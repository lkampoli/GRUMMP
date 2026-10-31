#include "GR_config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GR_misc.h"
#include "GR_EntContainer.h"
#define TOPO_DIM 2
#define SPACE_DIM 3
#include "GR_SurfMesh.h"

void readFile_Surface(
                       const char * const strBaseFileName,
                       GR_index_t& iNumVerts,
                       GR_index_t& iNumFaces,
                       GR_index_t& iNumCells,
                       GR_index_t& iNumIntBdryFaces,
                       bool& qFaceVert,
                       bool& qCellVert,
                       bool& qFaceCell,
                       bool& qCellFace,
                       bool& qCellRegion,
                       bool& qIntBFaceFace,
                       bool& qIntBFaceVert,
                       bool& qIntBFaceBC,
                       EntContainer<Vert>& ECVerts,
                       GR_index_t (*&a2iFaceVert)[TOPO_DIM],
                       GR_sindex_t (*&a2iFaceCell)[2],
                       GR_index_t (*&a2iCellVert)[TOPO_DIM+1],
                       GR_index_t (*&a2iCellFace)[TOPO_DIM+1],
                       int *&aiCellRegion,
                       GR_index_t (*&a2iIntBFaceFace)[2],
                       int *&aiIntBFaceBC,
                       GR_index_t (*&a2iIntBFaceVert)[TOPO_DIM]
                       )
{
  iNumVerts = 0;
  iNumFaces = 0;
  iNumCells = 0;
  iNumIntBdryFaces = 0;
  qFaceVert = false;
  qCellVert = false;
  qFaceCell = false;
  qCellFace = false;
  qCellRegion = false;
  a2iFaceCell = reinterpret_cast<GR_sindex_t(*)[2]> (NULL);
  a2iFaceVert = reinterpret_cast<GR_index_t(*)[TOPO_DIM]> (NULL);
  a2iCellVert = reinterpret_cast<GR_index_t(*)[TOPO_DIM+1]> (NULL);
  a2iCellFace = reinterpret_cast<GR_index_t(*)[TOPO_DIM+1]> (NULL);
  aiCellRegion = reinterpret_cast<int *> (NULL);
  qIntBFaceFace = false;
  qIntBFaceVert = false;
  qIntBFaceBC   = false;
  a2iIntBFaceFace = reinterpret_cast<GR_index_t(*)[2]> (NULL);
  aiIntBFaceBC    = reinterpret_cast<int *> (NULL);
  a2iIntBFaceVert = reinterpret_cast<GR_index_t(*)[TOPO_DIM]> (NULL);
 
  int iVertOffset = 0, iFaceOffset = 0;
  int iCellOffset = 0, iBFaceOffset = 1;
  int res = -1;
  FILE *pFInFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  char strTemp[1024];
  snprintf(strTemp, 1024, "pwm");
  const char *pcWhere = strstr(strBaseFileName, strTemp);
  if (pcWhere) {
    snprintf(strFileName, 1024, "%s", strBaseFileName);
  } else {
    snprintf(strFileName, 1024, "%s.pwm", strBaseFileName);
  }

  /*newfile pwm*/
  printf("Opening input file %s\n", strFileName);
  if (pFInFile) fclose(pFInFile);
  pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile)
    vFatalError("Couldn't open input file for reading",
                "surface mesh input");

  /*1994*/
  res = fscanf(pFInFile, "1994\n");
  assert(res >= 0);  /*PWI Mesh exported from Loom*/
  res = fscanf(pFInFile, "PWI Mesh exported from Loom\n");
  assert(res >= 0);  /*drop_int*/
  res = fscanf(pFInFile, "%*d\n");
  assert(res >= 0);  /*drop_int*/
  res = fscanf(pFInFile, "%*d\n");
  assert(res >= 0);  /*nverts*/
  res = fscanf(pFInFile, "%u\n", &iNumVerts);
  assert(res >= 0);  /*drop_int nfaces drop_int drop_int drop_int drop_int drop_int */
  res = fscanf(pFInFile, "%*d %u %*d %*d %*d %*d %*d \n", &iNumFaces);
  assert(res >= 0);
  if (iNumVerts == 0)
    vFatalError("Tried to read vertex data without specifying number of vertices",
                "surface mesh input");
  ECVerts.reserve(iNumVerts);
  for (GR_index_t iV = 0; iV < iNumVerts; iV++) {
    GR_index_t iRealVert = iV;
    double dXX, dYY;
    double dZZ;
    /*verts: coords*/
    res = fscanf(pFInFile, "%lf%lf%lf\n", &dXX, &dYY, &dZZ);
    assert(res >= 0);    {
      double adCoord[] = {dXX, dYY, dZZ};
      (ECVerts.getEntry(iRealVert))->setCoords(SPACE_DIM, adCoord);
    }
  }


  if (iNumFaces == 0)
    vFatalError("Tried to read face data without specifying number of faces",
                "surface mesh input");
  /* Allocate both the face-vert and face-cell arrays if they haven't been.
     One or both may have to be deleted later. */
  if (! a2iFaceVert)
    a2iFaceVert = new GR_index_t[iNumFaces][TOPO_DIM];
  if (! a2iFaceCell)
    a2iFaceCell = new GR_sindex_t[iNumFaces][2];
  if (iNumVerts == 0)
    vWarning("Number of verts unknown while reading faces; sanity checking may be crippled.");
  if (iNumCells == 0)
    vWarning("Number of cells unknown while reading faces; sanity checking may be crippled.");
  for (GR_index_t iF = 0; iF < iNumFaces; iF++) {
    GR_index_t iRealFace = iF;
    GR_index_t iVertA = 2*iNumVerts, iVertB = 2*iNumVerts;
    /*faces: verts*/
    res = fscanf(pFInFile, "%u%u\n", &iVertA, &iVertB);
    assert(res >= 0);    /* Remove an offset from the vertex indices if needed. */
    iVertA -= iVertOffset;
    iVertB -= iVertOffset;
    if (((iNumVerts > 0) && 
          (iVertA >= iNumVerts || 
           iVertB >= iNumVerts) ) ) 
      vFatalError("Vertex index out of range",
                "surface mesh input");
    (a2iFaceVert)[iRealFace][0] = iVertA;
    (a2iFaceVert)[iRealFace][1] = iVertB;

  }

  qFaceVert = true;
  /* Get rid of the face-vert array if it isn't used. */
  if (!(qFaceVert) && a2iFaceVert) {
    delete [] a2iFaceVert;
    a2iFaceVert = reinterpret_cast<GR_index_t(*)[TOPO_DIM]> (NULL);
  }

  /* Get rid of the face-cell array if it isn't used. */
  if (!(qFaceCell) && (a2iFaceCell)){
    delete [] a2iFaceCell;
    a2iFaceCell = reinterpret_cast<GR_sindex_t(*)[2]> (NULL);
  }

  /* Get rid of the cell-vert array if it isn't used. */
  if (!(qCellVert) && (a2iCellVert)) {
    delete [] a2iCellVert;
    a2iCellVert = reinterpret_cast<GR_index_t(*)[TOPO_DIM+1]> (NULL);
  }

  /* Get rid of the cell-face array if it isn't used. */
  if (!(qCellFace) && (a2iCellFace)) {
    delete [] a2iCellFace;
    a2iCellFace = reinterpret_cast<GR_index_t(*)[TOPO_DIM+1]> (NULL);
  }

  /* Get rid of the cell-region array if it isn't used. */
  if (!(qCellRegion) && (aiCellRegion)) {
    delete [] aiCellRegion;
    aiCellRegion = reinterpret_cast<int *> (NULL);
  }

  /* The following lines prevent compiler complaints. */
  iBFaceOffset++;
  iFaceOffset++;
  iCellOffset++;
  iVertOffset++;
  fclose(pFInFile);
}

void
writeFile_Surface(SurfMesh& OutMesh, const char strBaseFileName[],
		  const char strExtraFileSuffix[])
{
  int iVertOffset = 0, iFaceOffset = 0;
  int iCellOffset = 0, iBFaceOffset = 1;
  FILE *pFOutFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s.pwm%s", strBaseFileName, strExtraFileSuffix);

  /*newfile pwm*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "surface mesh output");

  /*1994*/
  fprintf(pFOutFile, "1994\n");
  /*PWI Mesh exported from Loom*/
  fprintf(pFOutFile, "PWI Mesh exported from Loom\n");
  /*drop_int*/
  fprintf(pFOutFile, "drop_int\n");
  /*drop_int*/
  fprintf(pFOutFile, "drop_int\n");
  /*nverts*/
  fprintf(pFOutFile, "%d\n", OutMesh.getNumVerts());
  /*drop_int nfaces drop_int drop_int drop_int drop_int drop_int */
  fprintf(pFOutFile,
	  "drop_int %d drop_int drop_int drop_int drop_int drop_int \n",
	  OutMesh.getNumFaces());

  for (GR_index_t iV = 0; iV < OutMesh.getNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.14g %.14g %.14g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y(), OutMesh.getVert(iV)->z());

  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
    Face* pF = OutMesh.getFace(iF);
    /*faces: verts*/
    fprintf(pFOutFile, "%d %d\n",
	    OutMesh.getVertIndex(pF->getVert(0)) + iVertOffset,
	    OutMesh.getVertIndex(pF->getVert(1)) + iVertOffset);
  }

  /* Make sure compilers won't complain that these variables */
  /* aren't used. */
  iVertOffset++;
  iFaceOffset++;
  iBFaceOffset++;
  iCellOffset++;
  fclose(pFOutFile);
}
