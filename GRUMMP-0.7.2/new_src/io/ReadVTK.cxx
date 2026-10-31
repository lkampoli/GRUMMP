#include "GR_config.h"
#include <stdio.h>
#include <stdlib.h>
#include "GR_misc.h"
#define TOPO_DIM 3
#define SPACE_DIM 3
#include "GR_VolMesh.h"

#ifndef VTK_TETRA
#define VTK_LINE 3
#define VTK_TRIANGLE 5
#define VTK_QUAD 9
#define VTK_TETRA 10
#define VTK_HEXAHEDRON 12
#define VTK_WEDGE 13
#define VTK_PYRAMID 14
#endif

#define DISCARD_LINE do {			\
    lineNum++;					\
    int i__i = fscanf(pFInFile, "%*[^\n]\n");	\
    i__i++;					\
  } while(0)

#define CHECK_READ(b, a) do {						\
    int _res_ = (a);							\
    lineNum++;								\
    if (_res_ != b) {							\
      logMessage(0, "File format error in function %s\n", __func__);	\
      logMessage(0, "File %s, line %d: \n", fileName, lineNum);		\
      return;					\
    }									\
  } while(0)

void
readVTK_ASCII(const char * const fileName, VolMesh * const pVM)
{
  FILE *pFInFile = fopen(fileName, "r");

  GR_index_t lineNum = 1;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;

  GR_index_t nVerts;
  CHECK_READ(1, (fscanf(pFInFile, "POINTS %u %*s\n", &nVerts)));

  double coordsTemp[3];
  for (GR_index_t i = 0; i < nVerts; i++) {
    CHECK_READ(
	3,
	(fscanf(pFInFile, "%lf%lf%lf\n", coordsTemp, coordsTemp + 1,
		coordsTemp + 2)));
    pVM->createVert(coordsTemp);
  }

  // Read in the actual connectivity info
  int iConnectSize, iNEnts;
  CHECK_READ(2, (fscanf(pFInFile, "CELLS %d %d\n", &iNEnts, &iConnectSize)));
  iConnectSize -= iNEnts;

  // aiConnect and aiOffset can be leaked if there's a read error.
  int *aiConnect = new int[iConnectSize];
  int *aiOffset = new int[iNEnts + 1];
  int iIndex = 0;
  aiOffset[0] = 0;
  for (int i = 0; i < iNEnts; i++) {
    int iVertsThisEnt;
    CHECK_READ(1, (fscanf(pFInFile, "%d", &iVertsThisEnt)));
    aiOffset[i + 1] = aiOffset[i] + iVertsThisEnt;
    for (int j = 0; j < iVertsThisEnt; j++) {
      CHECK_READ(1, (fscanf(pFInFile, "%d", &(aiConnect[iIndex]))));
      iIndex++;
    }
    assert(iIndex == aiOffset[i + 1]);
  }

  // aiTopo can be leaked if there's a read error.
  int *aiTopo = new int[iNEnts];
  // Read the cell type data as well
  {
    int iCheck;
    CHECK_READ(1, (fscanf(pFInFile, "%*s %d\n", &iCheck)));
    assert(iCheck == iNEnts);
  }

  int edgesRead = 0, trisRead = 0, quadsRead = 0, tetsRead = 0, pyrsRead = 0,
      prismsRead = 0, hexesRead = 0;
  for (int i = 0; i < iNEnts; i++) {
    CHECK_READ(1, (fscanf(pFInFile, "%d", &(aiTopo[i]))));
    switch (aiTopo[i])
      {
      case VTK_LINE:
	edgesRead++;
	break;
      case VTK_TRIANGLE:
	trisRead++;
	break;
      case VTK_QUAD:
	quadsRead++;
	break;
      case VTK_TETRA:
	tetsRead++;
	break;
      case VTK_PYRAMID:
	pyrsRead++;
	break;
      case VTK_WEDGE:
	prismsRead++;
	break;
      case VTK_HEXAHEDRON:
	hexesRead++;
	break;
      default:
	break;
      }
  }
  logMessage(
      2,
      "Read %d edges, %d tris, %d quads, %d tets, %d pyrs, %d prisms, %d hexes\n",
      edgesRead, trisRead, quadsRead, tetsRead, pyrsRead, prismsRead,
      hexesRead);

  // Create entities.  If any creation fails, forget about the whole thing.

  Vert* verts[8];
  bool qExist;
  for (int i = 0; i < iNEnts; i++) {
    if (aiTopo[i] == VTK_LINE)
      continue; // Don't try to create these.
    int ii = 0, iVertsThisEnt = aiOffset[i + 1] - aiOffset[i];
    for (; ii < iVertsThisEnt; ii++) {
      verts[ii] = pVM->getVert(aiConnect[aiOffset[i] + ii]);
    }
    switch (aiTopo[i])
      {
      case VTK_TRIANGLE:
      case VTK_QUAD:
	pVM->createFace(verts, iVertsThisEnt);
	break;
      case VTK_TETRA:
	pVM->createTetCell(qExist, verts[0], verts[1], verts[2], verts[3], 1);
	assert(!qExist);
	break;
      case VTK_PYRAMID:
	pVM->createPyrCell(qExist, verts[0], verts[1], verts[2], verts[3],
			   verts[4], 1);
	assert(!qExist);
	break;
      case VTK_WEDGE:
	pVM->createPrismCell(qExist, verts[0], verts[1], verts[2], verts[3],
			     verts[4], verts[5], 1);
	assert(!qExist);
	break;
      case VTK_HEXAHEDRON:
	pVM->createHexCell(qExist, verts[0], verts[1], verts[2], verts[3],
			   verts[4], verts[5], verts[6], verts[7], 1);
	assert(!qExist);
	break;
      default:
	break;
      }
  }

  delete[] aiTopo;
  delete[] aiOffset;
  delete[] aiConnect;

  // Create BFace's where needed for completeness.
  for (int i = pVM->getNumFaces() - 1; i >= 0; i--) {
    Face *pF = pVM->getFace(i);
    if (!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid()) {
      pVM->createBFace(pF);
    }
  }

  fclose(pFInFile);
}
