#include <stdio.h>
#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_BFace.h"
#include "GR_Vertex.h"
#include "GR_Mesh2D.h"
#include "GR_BdryPatch2D.h"

void
writeFEA(Mesh2D& OutMesh, const char strBaseFileName[], const int iOrder,
	 const char strExtraFileSuffix[])
{
  assert(OutMesh.isSimplicial());
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s%s.fea", strBaseFileName, strExtraFileSuffix);

  /*newfile fea*/
  printf("Opening output file %s\n", strFileName);
  FILE *pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "2d mesh output");

  /*nverts ncells nfaces nbfaces*/
  fprintf(pFOutFile, "%u %u %u %u\n", OutMesh.iNumVerts(),
	  OutMesh.getNumCells(), OutMesh.getNumFaces(),
	  OutMesh.getNumBdryFaces());

  for (GR_index_t iV = 0; iV < OutMesh.iNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.16g %.16g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y());

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    Cell* pC = OutMesh.getCell(iC);
    for (int iV = 0; iV < pC->getNumVerts(); iV++) {
      fprintf(pFOutFile, "%u ", OutMesh.iVertIndex(pC->getVert(iV)));
    }
    /*cells: verts*/
    fprintf(pFOutFile, "\n");
  }

  fclose(pFOutFile);

  // Need to write the bdry data
  if (iOrder > 2) {
    sprintf(strFileName, "%s%s.bpts", strBaseFileName, strExtraFileSuffix);

    logMessage(1, "Opening boundary point output file %s\n", strFileName);
    pFOutFile = fopen(strFileName, "w");
    if (NULL == pFOutFile)
      vFatalError("Couldn't open boundary point output file for writing",
		  "2d mesh output");

    int iNumPoints = iOrder + 1;
    double adS[iNumPoints];
    for (int ii = 0; ii < iNumPoints; ii++) {
      adS[ii] = double(ii) / (iNumPoints - 1);
    }

    for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); iBF++) {
      BFace* pBF = OutMesh.getBFace(iBF);
      assert(
	  (pBF->getType() == Cell::eBdryEdge)
	      || (pBF->getType() == Cell::eIntBdryEdge));

      BdryPatch2D *pBP = dynamic_cast<BdryPatch2D*>(pBF->getPatchPointer());

      double adPointLoc[2];

      const double* const adVertLocA = pBF->getVert(0)->getCoords();
      const double* const adVertLocB = pBF->getVert(1)->getCoords();

      fprintf(pFOutFile, "%u\t%u\t%u\t%s\n", iBF,
	      OutMesh.iVertIndex(pBF->getVert(0)),
	      OutMesh.iVertIndex(pBF->getVert(1)),
	      pBF->getBdryCondition() == 1 ? "wall" : "far field");

      for (int ii = 0; ii < iNumPoints; ii++) {
	pBP->vGeodesicDistanceRatioBetweenTwoPoints(adS[ii], adVertLocA,
						    adVertLocB, adPointLoc);
	fprintf(pFOutFile, " %.16g %.16g\n", adPointLoc[0], adPointLoc[1]);
      }
      fprintf(pFOutFile, "\n");
    }
    fclose(pFOutFile);
  }

}
