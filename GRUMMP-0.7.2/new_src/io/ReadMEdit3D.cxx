#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "GR_misc.h"
#include "GR_VolMesh.h"

void
readMEdit3D(const char * const fileName, VolMesh * const pVM)
{
  char strFileName[1024];
  const char *pcWhere = strstr(fileName, "mesh");
  if (pcWhere) {
    snprintf(strFileName, 1024, "%s", fileName);
  }
  else {
    snprintf(strFileName, 1024, "%s.mesh", fileName);
  }

  /*newfile mesh*/
  printf("Opening input file %s\n", strFileName);
  FILE *pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile)
    vFatalError("Couldn't open input file for reading", "3d mesh input");

  /*MeshVersionFormatted 1*/
  int res = fscanf(pFInFile, "MeshVersionFormatted 1\n");
  assert(res >= 0);
  /*Dimension*/
  res = fscanf(pFInFile, "Dimension\n");
  assert(res >= 0);
  /*3*/
  res = fscanf(pFInFile, "3\n");
  assert(res >= 0);
  /*Vertices*/
  res = fscanf(pFInFile, "Vertices\n");
  assert(res >= 0);
  /*nverts*/
  GR_index_t iNumVerts = 0;
  res = fscanf(pFInFile, "%u\n", &iNumVerts);
  assert(res >= 0);

  if (iNumVerts == 0)
    vFatalError(
	"Tried to read vertex data without specifying number of vertices",
	"3d mesh input");
  for (GR_index_t iV = 0; iV < iNumVerts; iV++) {
    double dXX, dYY;
    double dZZ;
    /*verts: x y z discard_int*/
    res = fscanf(pFInFile, "%lf%lf%lf%*d\n", &dXX, &dYY, &dZZ);
    assert(res >= 0);
    pVM->createVert(dXX, dYY, dZZ);
  }

  /*Triangles*/
  res = fscanf(pFInFile, "Triangles\n");
  assert(res >= 0);
  /*nbfaces*/
  GR_index_t nBdryTris = 0;
  res = fscanf(pFInFile, "%u\n", &nBdryTris);
  assert(res >= 0);

  if (nBdryTris == 0)
    vFatalError(
	"Tried to read bdry face data without specifying number of bdry faces",
	"3d mesh input");
  GR_index_t (*bdryTriToVert)[3] = new GR_index_t[nBdryTris][3];
  GR_index_t *bdryCond = new GR_index_t[nBdryTris];
  for (GR_index_t iBF = 0; iBF < nBdryTris; iBF++) {
    GR_index_t iVertA = 2 * iNumVerts, iVertB = 2 * iNumVerts;
    GR_index_t iVertC = 2 * iNumVerts;
    GR_index_t thisBC = INT_MAX;
    /*bdryfaces: verts bc*/
    res = fscanf(pFInFile, "%u%u%u %u\n", &iVertA, &iVertB, &iVertC, &thisBC);
    assert(res >= 0);
    /* Remove offset from the vertex indices */
    iVertA -= 1;
    iVertB -= 1;
    iVertC -= 1;
    if (((iNumVerts > 0)
	&& (iVertA >= iNumVerts || iVertC >= iNumVerts || iVertB >= iNumVerts)))
      vFatalError("Vertex index out of range", "3d mesh input");
    bdryTriToVert[iBF][0] = iVertA;
    bdryTriToVert[iBF][1] = iVertB;
    bdryTriToVert[iBF][2] = iVertC;
    bdryCond[iBF] = thisBC;
  }

  /*Tetrahedra*/
  res = fscanf(pFInFile, "Tetrahedra\n");
  assert(res >= 0);
  /*ncells*/
  GR_index_t nTet = 0;
  res = fscanf(pFInFile, "%u\n", &nTet);
  assert(res >= 0);
  GR_index_t (*tetToVert)[4] = new GR_index_t[nTet][4];
  GR_index_t (*regions) = new GR_index_t[nTet];
  if (nTet == 0)
    vFatalError("Tried to read cell data without specifying number of cells",
		"3d mesh input");
  for (GR_index_t iC = 0; iC < nTet; iC++) {
    GR_index_t iVertA = 2 * iNumVerts, iVertB = 2 * iNumVerts, iVertC = 2
	* iNumVerts, iVertD = 2 * iNumVerts;
    int iReg = 0;
    /*cells: verts region*/
    res = fscanf(pFInFile, "%u%u%u%u %u\n", &iVertA, &iVertB, &iVertC, &iVertD,
		 &iReg);
    assert(res >= 0);
    /* Store the region that the cell belongs to. */
    iReg += iDefaultRegion;
    if (iReg <= 0 || iReg >= 32)
      vFatalError("Cell tagged with invalid region.", "3d mesh input");
    iVertA -= 1;
    iVertB -= 1;
    iVertC -= 1;
    iVertD -= 1;
    if (((iNumVerts > 0)
	&& (iVertA >= iNumVerts || iVertB >= iNumVerts || iVertD >= iNumVerts
	    || iVertC >= iNumVerts)))
      vFatalError("Vertex index out of range", "3d mesh input");
    tetToVert[iC][0] = iVertA;
    tetToVert[iC][1] = iVertB;
    tetToVert[iC][2] = iVertC;
    tetToVert[iC][3] = iVertD;
    regions[iC] = iReg;
  }
  // Unfortunately, the bdry face info as read from Intel medit files isn't
  // useful to us, not least because it's got two copies of each internal face.
  pVM->convertFromCellVert(nTet, nBdryTris, tetToVert, regions, bdryTriToVert);
  if (bdryTriToVert)
    delete[] bdryTriToVert;
  if (bdryCond)
    delete[] bdryCond;
  if (tetToVert)
    delete[] tetToVert;
  if (regions)
    delete[] regions;
}
