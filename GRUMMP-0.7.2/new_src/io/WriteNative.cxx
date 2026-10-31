#include <stdio.h>
#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_BFace.h"
#include "GR_Vertex.h"
#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"

void
writeNative(Mesh2D& OutMesh, const char strBaseFileName[],
	    const bool qHighOrder, const char strExtraFileSuffix[])
{
  if (OutMesh.getNumIntBdryFaces() == 0) {
    writeNativeMESH(OutMesh, strBaseFileName, qHighOrder, strExtraFileSuffix);
  }
  else {
    writeNativeGRMESH(OutMesh, strBaseFileName, qHighOrder, strExtraFileSuffix);
  }
}

void
writeNative(VolMesh& OutMesh, const char strBaseFileName[],
	    const bool /*qHighOrder*/, const char strExtraFileSuffix[])
{
  if (OutMesh.getNumIntBdryFaces() == 0) {
    writeNativeMESH(OutMesh, strBaseFileName, strExtraFileSuffix);
  }
  else {
    writeNativeGRMESH(OutMesh, strBaseFileName, strExtraFileSuffix);
  }
}

void
writeNativeMESH(Mesh2D& OutMesh, const char strBaseFileName[],
		const bool qHighOrder, const char strExtraFileSuffix[])
{
  //  assert(OutMesh.qSimplicial());
  //assert(OutMesh.getNumIntBdryFaces() == 0);
  FILE *pFOutFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s%s.mesh", strBaseFileName, strExtraFileSuffix);

  /*newfile mesh*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "2d mesh output");

  /*ncells nfaces nbfaces nverts*/
  fprintf(pFOutFile, "%u %u %u %u\n", OutMesh.getNumCells(),
	  OutMesh.getNumFaces(), OutMesh.getNumBdryFaces(),
	  OutMesh.iNumVerts());

  for (GR_index_t iV = 0; iV < OutMesh.iNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%20.16lf %20.16lf\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y());

  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
    Face* pF = OutMesh.getFace(iF);
    GR_sindex_t iCA, iCB;
    Cell *pC = pF->getLeftCell();

    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
	iCA = -(dynamic_cast<BFace*>(pC))->getBdryCondition();

	break;

      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	Face * pF2;
	pF2 = pC->getFace(0);
	if (pF == pF2)
	  pF2 = pC->getFace(1);
	Cell * pC2;
	pC2 = pF2->getOppositeCell(pC);
	iCA = OutMesh.getCellIndex(pC2);
	break;

      default:
	iCA = OutMesh.getCellIndex(pC);
	break;
      }

    pC = pF->getRightCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:

	iCB = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:

	Face * pF2;
	pF2 = pC->getFace(0);
	if (pF == pF2)
	  pF2 = pC->getFace(1);
	Cell * pC2;
	pC2 = pF2->getOppositeCell(pC);
	iCB = OutMesh.getCellIndex(pC2);
	break;

      default:
	iCB = OutMesh.getCellIndex(pC);
	break;
      }
    /*faces: cells verts*/
    fprintf(pFOutFile, "%d %d %u %u\n", iCA, iCB,
	    OutMesh.iVertIndex(pF->getVert(0)),
	    OutMesh.iVertIndex(pF->getVert(1)));
  }

  for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); iBF++) {
    BFace* pBF = OutMesh.getBFace(iBF);
    /*bdryfaces: face bc verts*/
    fprintf(pFOutFile, "%u %d %u %u\n", OutMesh.getFaceIndex(pBF->getFace(0)),
	    pBF->getBdryCondition(), OutMesh.iVertIndex(pBF->getVert(0)),
	    OutMesh.iVertIndex(pBF->getVert(1)));
  }

//  for (GR_index_t iBF = 0; iBF < OutMesh.getNumIntBdryFaces(); iBF++) {
//    IntBdryEdge *pIBE = OutMesh.getIntBdryEdge(iBF);
//    /*intbdryfaces: faces verts*/
//    fprintf(pFOutFile, "%u %u %u %u\n", OutMesh.getFaceIndex(pIBE->getFace(0)),
//            OutMesh.getFaceIndex(pIBE->getFace(1)),
//            OutMesh.iVertIndex(pIBE->getVert(0)),
//            OutMesh.iVertIndex(pIBE->getVert(1)));
//  }

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    Cell* pC = OutMesh.getCell(iC);
    /*cells: region*/
    fprintf(pFOutFile, "%d\n", pC->getRegion());
  }

  fclose(pFOutFile);

  if (qHighOrder) {
    OutMesh.writeBdryPatchFile(strBaseFileName, strExtraFileSuffix);
  }

}
void
writeNativeGRMESH(Mesh2D& OutMesh, const char strBaseFileName[],
		  const bool qHighOrder, const char strExtraFileSuffix[])
{
  //  assert(OutMesh.qSimplicial());
  FILE *pFOutFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s%s.grmesh", strBaseFileName, strExtraFileSuffix);

  /*newfile mesh*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "2d mesh output");

  /*ncells nfaces nbfaces nibfaces nverts*/
  fprintf(pFOutFile, "%u %u %u %u %u\n", OutMesh.getNumCells(),
	  OutMesh.getNumFaces(), OutMesh.getNumBdryFaces(),
	  OutMesh.getNumIntBdryFaces(), OutMesh.iNumVerts());

  for (GR_index_t iV = 0; iV < OutMesh.iNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.16g %.16g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y());

  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
    Face* pF = OutMesh.getFace(iF);
    GR_sindex_t iCA, iCB;
    Cell *pC = pF->getLeftCell();

    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
	iCA = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:

	Face * pF2;
	pF2 = pC->getFace(0);
	if (pF == pF2)
	  pF2 = pC->getFace(1);
	Cell * pC2;
	pC2 = pF2->getOppositeCell(pC);
	iCA = OutMesh.getCellIndex(pC2);
	break;

      default:
	iCA = OutMesh.getCellIndex(pC);
	break;
      }

    pC = pF->getRightCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
	iCB = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:

	Face * pF2;
	pF2 = pC->getFace(0);
	if (pF == pF2)
	  pF2 = pC->getFace(1);
	Cell * pC2;
	pC2 = pF2->getOppositeCell(pC);
	iCB = OutMesh.getCellIndex(pC2);
	break;

      default:
	iCB = OutMesh.getCellIndex(pC);
	break;
      }
    /*faces: cells verts*/
    fprintf(pFOutFile, "%d %d %u %u\n", iCA, iCB,
	    OutMesh.iVertIndex(pF->getVert(0)),
	    OutMesh.iVertIndex(pF->getVert(1)));
  }

  for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); iBF++) {
    BFace* pBF = OutMesh.getBFace(iBF);
    /*bdryfaces: face bc verts*/
    fprintf(pFOutFile, "%u %d %u %u\n", OutMesh.getFaceIndex(pBF->getFace(0)),
	    pBF->getBdryCondition(), OutMesh.iVertIndex(pBF->getVert(0)),
	    OutMesh.iVertIndex(pBF->getVert(1)));
  }

  for (GR_index_t iBF = 0; iBF < OutMesh.getNumIntBdryFaces(); iBF++) {
    IntBdryEdge *pIBE = OutMesh.getIntBFace(iBF);
    /*intbdryfaces: faces verts*/
    fprintf(pFOutFile, "%u %u %u %u\n", OutMesh.getFaceIndex(pIBE->getFace(0)),
	    OutMesh.getFaceIndex(pIBE->getFace(1)),
	    OutMesh.iVertIndex(pIBE->getVert(0)),
	    OutMesh.iVertIndex(pIBE->getVert(1)));
  }

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    Cell* pC = OutMesh.getCell(iC);
    /*cells: region*/
    fprintf(pFOutFile, "%d\n", pC->getRegion());
  }

  fclose(pFOutFile);

  if (qHighOrder) {
    OutMesh.writeBdryPatchFile(strBaseFileName, strExtraFileSuffix);
  }

}

void
Mesh2D::writeBdryPatchFile(const char strBaseFileName[],
			   const char strExtraFileSuffix[]) const
{

  // I can't imagine why this assertion is needed.
  //  assert(qSimplicial());

  FILE *output_file = NULL;
  char filename[1024];
  sprintf(filename, "%s.bpatch%s", strBaseFileName, strExtraFileSuffix);

  /*newfile mesh*/
  logMessage(0, "Opening boundary patch output file %s\n", filename);

  if (output_file)
    fclose(output_file);

  output_file = fopen(filename, "w");
  if (NULL == output_file)
    vFatalError("Could not open boundary patch output file for writing",
		"Mesh2D::vWriteBdryPatchFile()");

  int num_gauss_points = 18;
  double beg_ratio[num_gauss_points], mid_ratio[num_gauss_points],
      end_ratio[num_gauss_points];

  // useful constant.
  double S = 0.5 * (1. + sqrt(1. / 3.));

  // one gauss point
  beg_ratio[0] = 0.;
  mid_ratio[0] = 0.5;
  end_ratio[0] = 1.;

  // two gauss points
  beg_ratio[1] = beg_ratio[0];
  mid_ratio[1] = 1. - S;
  end_ratio[1] = 0.5;

  beg_ratio[2] = end_ratio[1];
  mid_ratio[2] = S;
  end_ratio[2] = end_ratio[0];

  // three gauss points
  beg_ratio[3] = beg_ratio[0];
  mid_ratio[3] = 0.5 - sqrt(0.15);
  end_ratio[3] = 5. / 18.;

  beg_ratio[4] = end_ratio[3];
  mid_ratio[4] = 0.5;
  end_ratio[4] = 13. / 18.;

  beg_ratio[5] = end_ratio[4];
  mid_ratio[5] = 0.5 + sqrt(0.15);
  end_ratio[5] = end_ratio[0];

  // two times one gauss points
  beg_ratio[6] = beg_ratio[0];
  mid_ratio[6] = 0.25;
  end_ratio[6] = 0.5;

  beg_ratio[7] = end_ratio[6];
  mid_ratio[7] = 0.75;
  end_ratio[7] = end_ratio[0];

  // two times two gauss points
  beg_ratio[8] = beg_ratio[0];
  mid_ratio[8] = 0.5 * (1. - S);
  end_ratio[8] = 0.25;

  beg_ratio[9] = end_ratio[8];
  mid_ratio[9] = 0.5 * S;
  end_ratio[9] = 0.5;

  beg_ratio[10] = end_ratio[9];
  mid_ratio[10] = 1. - mid_ratio[9];
  end_ratio[10] = 0.75;

  beg_ratio[11] = end_ratio[10];
  mid_ratio[11] = 0.5 * (1. + S);
  end_ratio[11] = end_ratio[0];

  // two times three gauss points
  beg_ratio[12] = beg_ratio[0];
  mid_ratio[12] = 0.5 * (0.5 - sqrt(0.15));
  end_ratio[12] = 5. / 36.;

  beg_ratio[13] = end_ratio[12];
  mid_ratio[13] = 0.25;
  end_ratio[13] = 13. / 36.;

  beg_ratio[14] = end_ratio[13];
  mid_ratio[14] = 0.5 * (0.5 + sqrt(0.15));
  end_ratio[14] = 0.5;

  beg_ratio[15] = end_ratio[14];
  mid_ratio[15] = 0.5 * (1.5 - sqrt(0.15));
  end_ratio[15] = 0.5 + 5. / 36.;

  beg_ratio[16] = end_ratio[15];
  mid_ratio[16] = 0.75;
  end_ratio[16] = 0.5 + 13. / 36.;

  beg_ratio[17] = end_ratio[16];
  mid_ratio[17] = 0.5 * (1.5 + sqrt(0.15));
  end_ratio[17] = end_ratio[0];

  //   for(int i = 0; i < 18; i++) {
  //     printf("i = %d, beg = %e, mid = %e, end = %e\n",
  // 	   i, beg_ratio[i], mid_ratio[i], end_ratio[i]);
  //   }

  ///

  int num_faces = getNumFaces();
  double beg_param, end_param, gauss_param, face_length;

  GRCurve* curve = NULL;
  Face* face = NULL;
  BdryEdgeBase* bdry_edge = NULL;

  CubitVector gauss_coord, normal;
  bool reverse_normal;

  for (int i = 0; i < num_faces; i++) {

    face = getFace(i);
    if (face->isDeleted())
      continue;

    reverse_normal = false;

    if (face->isBdryFace()) {

      bdry_edge = dynamic_cast<BdryEdgeBase*>(face->getLeftCell());
      if (!bdry_edge) {
	bdry_edge = dynamic_cast<BdryEdgeBase*>(face->getRightCell());
	assert(bdry_edge);
	if (!bdry_edge->isForward())
	  reverse_normal = true;
      }
      else {
	if (!bdry_edge->isForward())
	  reverse_normal = true;
      }

      curve = bdry_edge->getCurve();
      if (bdry_edge->isForward()) {
	beg_param = bdry_edge->getVert0Param();
	end_param = bdry_edge->getVert1Param();
      }
      else {
	beg_param = bdry_edge->getVert1Param();
	end_param = bdry_edge->getVert0Param();
      }
      assert(iFuzzyComp(beg_param, end_param) == -1);

      face_length = curve->get_curve_geom()->arc_length(beg_param, end_param);

    }

    else
      continue;

    fprintf(output_file, "%d\n", i);

    for (int j = 0; j < num_gauss_points; j++) {

      gauss_param = curve->get_curve_geom()->param_at_arc_length(
	  beg_param, mid_ratio[j] * face_length);

      curve->get_curve_geom()->coord_at_param(gauss_param, gauss_coord);
      curve->get_curve_geom()->unit_normal(gauss_param, normal);
      if (reverse_normal)
	normal *= -1.;

      fprintf(output_file, "%.16g %.16g %.16g %.16g %.16g\n", gauss_coord.x(),
	      gauss_coord.y(), normal.x(), normal.y(),
	      (end_ratio[j] - beg_ratio[j]) * face_length);

    }

    fprintf(output_file, "\n");

  }

  fclose(output_file);

}
void
writeNativeMESH(VolMesh& OutMesh, const char strBaseFileName[],
		const char strExtraFileSuffix[])
{
  // assert(OutMesh.isSimplicial());

  FILE *pFOutFile = nullptr;
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s%s.vmesh", strBaseFileName, strExtraFileSuffix);

  /*newfile vmesh*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d mesh output");

  /*ncells nfaces nbfaces nibfaces nverts*/
  fprintf(pFOutFile, "%d %d %d %d\n", OutMesh.getNumCells(),
	  OutMesh.getNumFaces(), OutMesh.getNumBdryFaces(),
	  OutMesh.getNumVerts());

  for (GR_index_t iV = 0; iV < OutMesh.getNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.14g %.14g %.14g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y(), OutMesh.getVert(iV)->z());

  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
    Face* pF = OutMesh.getFace(iF);
    GR_index_t iCA, iCB;
    Cell *pC = pF->getLeftCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
//    case Cell::eIntBdryEdge:
//    case Cell::eIntTriBFace:
//    case Cell::eIntQuadBFace:
	iCA = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	{
	  Face *pFOther = pC->getFace(0);
	  if (pF == pFOther)
	    pFOther = pC->getFace(1);
	  Cell *pCOther = pFOther->getOppositeCell(pC);
	  iCA = OutMesh.getCellIndex(pCOther);
	}
	break;
      default:
	iCA = OutMesh.getCellIndex(pC);
	break;
      }

    pC = pF->getRightCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
//    case Cell::eIntBdryEdge:
//    case Cell::eIntTriBFace:
//    case Cell::eIntQuadBFace:
	iCB = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	{
	  Face *pFOther = pC->getFace(0);
	  if (pF == pFOther)
	    pFOther = pC->getFace(1);
	  Cell *pCOther = pFOther->getOppositeCell(pC);
	  iCB = OutMesh.getCellIndex(pCOther);
	}
	break;
      default:
	iCB = OutMesh.getCellIndex(pC);
	break;
      }
    /*faces: cells verts*/
    fprintf(pFOutFile, "%d %d ", iCA, iCB);
    for (int ii = 0; ii < pF->getNumVerts(); ii++) {
      fprintf(pFOutFile, "%d ", OutMesh.getVertIndex(pF->getVert(ii)));
    }
    fprintf(pFOutFile, "\n");
  }

  for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); iBF++) {
    BFace* pBF = OutMesh.getBFace(iBF);
    /*bdryfaces: face bc verts*/
    fprintf(pFOutFile, "%d %d ", OutMesh.getFaceIndex(pBF->getFace(0)),
	    pBF->getBdryCondition());
    for (int ii = 0; ii < pBF->getNumVerts(); ii++) {
      fprintf(pFOutFile, "%d ", OutMesh.getVertIndex(pBF->getVert(ii)));
    }
    fprintf(pFOutFile, "\n");
  }

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    fprintf(pFOutFile, "%d\n", OutMesh.getCell(iC)->getRegion());
  }
  fclose(pFOutFile);
}
void
writeNativeGRMESH(VolMesh& OutMesh, const char strBaseFileName[],
		  const char strExtraFileSuffix[])
{
  assert(OutMesh.isSimplicial());

  FILE *pFOutFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s%s.grvmesh", strBaseFileName, strExtraFileSuffix);

  /*newfile vmesh*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d mesh output");

  /*ncells nfaces nbfaces nibfaces nverts*/
  fprintf(pFOutFile, "%d %d %d %d %d\n", OutMesh.getNumCells(),
	  OutMesh.getNumFaces(), OutMesh.getNumBdryFaces(),
	  OutMesh.getNumIntBdryFaces(), OutMesh.getNumVerts());

  for (GR_index_t iV = 0; iV < OutMesh.getNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.14g %.14g %.14g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y(), OutMesh.getVert(iV)->z());

  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
    Face* pF = OutMesh.getFace(iF);
    GR_index_t iCA, iCB;
    Cell *pC = pF->getLeftCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
//    case Cell::eIntBdryEdge:
//    case Cell::eIntTriBFace:
//    case Cell::eIntQuadBFace:
	iCA = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	{
	  Face *pFOther = pC->getFace(0);
	  if (pF == pFOther)
	    pFOther = pC->getFace(1);
	  Cell *pCOther = pFOther->getOppositeCell(pC);
	  iCA = OutMesh.getCellIndex(pCOther);
	}
	break;
      default:
	iCA = OutMesh.getCellIndex(pC);
	break;
      }

    pC = pF->getRightCell();
    switch (pC->getType())
      {
      case Cell::eBdryEdge:
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
//    case Cell::eIntBdryEdge:
//    case Cell::eIntTriBFace:
//    case Cell::eIntQuadBFace:
	iCB = -(dynamic_cast<BFace*>(pC))->getBdryCondition();
	break;
      case Cell::eIntBdryEdge:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	{
	  Face *pFOther = pC->getFace(0);
	  if (pF == pFOther)
	    pFOther = pC->getFace(1);
	  Cell *pCOther = pFOther->getOppositeCell(pC);
	  iCB = OutMesh.getCellIndex(pCOther);
	}
	break;
      default:
	iCB = OutMesh.getCellIndex(pC);
	break;
      }
    /*faces: cells verts*/
    fprintf(pFOutFile, "%d %d %d %d %d\n", iCA, iCB,
	    OutMesh.getVertIndex(pF->getVert(0)),
	    OutMesh.getVertIndex(pF->getVert(1)),
	    OutMesh.getVertIndex(pF->getVert(2)));
  }

  for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); iBF++) {
    BFace* pBF = OutMesh.getBFace(iBF);
    /*bdryfaces: face bc verts*/
    fprintf(pFOutFile, "%d %d %d %d %d\n",
	    OutMesh.getFaceIndex(pBF->getFace(0)), pBF->getBdryCondition(),
	    OutMesh.getVertIndex(pBF->getVert(0)),
	    OutMesh.getVertIndex(pBF->getVert(1)),
	    OutMesh.getVertIndex(pBF->getVert(2)));
  }
  for (GR_index_t iBF = 0; iBF < OutMesh.getNumIntBdryFaces(); iBF++) {
    BFace *pIBE = OutMesh.getIntBFace(iBF);
    /*intbdryfaces: faces verts*/
    fprintf(pFOutFile, "%d %d %d %d %d\n",
	    OutMesh.getFaceIndex(pIBE->getFace(0)),
	    OutMesh.getFaceIndex(pIBE->getFace(1)),
	    OutMesh.getVertIndex(pIBE->getVert(0)),
	    OutMesh.getVertIndex(pIBE->getVert(1)),
	    OutMesh.getVertIndex(pIBE->getVert(2)));
  }
  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    fprintf(pFOutFile, "%d\n", OutMesh.getCell(iC)->getRegion());
  }
  fclose(pFOutFile);
}

