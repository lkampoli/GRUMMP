/*
 * WriteTetGen.cxx
 *
 *  Created on: Apr 11, 2017
 *      Author: cfog
 */

#include <cstdio>
#include <assert.h>

#include "GR_VolMesh.h"

//-----------------------------------------------------------------------------
void
readTetGenFile(const char strBaseFileName[], VolMesh* const pVM, int*& aiBCList,
	       int* p_numFields, double** p_fieldData)
//-----------------------------------------------------------------------------
{
  FILE *pFInFile = NULL;
  char strFileName[1024];
  pVM->purgeAllEntities();

  // This file format actually writes four different files, so any suffix
  // on the base file name is irrelevant.

  // Need to identify shock faces, including which of the two faces
  // incident on the bdry face  is actually on the upstream side.  Likewise,
  // need to tag verts that need to be duplicated and actually duplicate
  // them.; eventually, I'll need a map between the upstream and downstream
  // copies of a vert (those on the perimeter can map to themselves).

  // Also, in the process, the internal bdry faces can become regular bdry
  // faces with regular bdry conditions.

  // Verts that have faces with more than one BC (ends of shocks) can get
  // tagged as interior.

  // Then it should be possible to write as normal.

  ///////////////////////////////////////////////////////////////
  // First file: node locations (no tags for now).
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.node", strBaseFileName);
  pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile) {
    vFatalError("Couldn't open input node file for reading",
		"3D Tetgen mesh input");
  }

  GR_index_t nVerts = 0;
  int hasTags = 0;
  fscanf(pFInFile, "%u 3 %d %d\n", &nVerts, p_numFields, &hasTags);
  assert(!*p_fieldData);
  (*p_fieldData) = new double[(*p_numFields) * nVerts];
  for (GR_index_t iV = 0; iV < nVerts; iV++) {
    GR_index_t actualVert = INT_MAX;
    double x = LARGE_DBL, y = LARGE_DBL, z = LARGE_DBL;
    fscanf(pFInFile, "%u %lf %lf %lf", &actualVert, &x, &y, &z);
    assert(actualVert > 0 && actualVert <= nVerts);
    actualVert--;
    Vert *pV = pVM->createVert(x, y, z);

    for (int ii = 0; ii < *p_numFields; ii++) {
      fscanf(pFInFile, "%lf", &((*p_fieldData)[ii + iV * (*p_numFields)]));
    }
    if (hasTags) {
      int tag = INT_MIN;
      fscanf(pFInFile, "%d", &tag);
      switch (tag)
	{
	case -1:
	  // Phantom
	  pV->markAsPhantom(true);
	  pV->setType(Vert::ePhantom);
	  break;
	case 0:
	  // Interior
	  pV->setType(Vert::eInterior);
	  break;
	case 10:
	case 11:
	case 12:
	case 13:
	  // Shock point of some sort; map directly to vert types
	  pV->setType(tag);
	  break;
	default:
	  // The rest are bdry variants; these will be handled separately.
	  pV->setType(Vert::eBdry);
	  break;
	}
    } // Done reading a tag, if there is one.
  } // Done reading all node data.
  fclose(pFInFile);
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Next: element connectivity (.ele)
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.ele", strBaseFileName);
  pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile) {
    vFatalError("Couldn't open input element file for reading",
		"3D Tetgen mesh input");
  }

  GR_index_t nTets = 0;
  fscanf(pFInFile, "%u 4 0\n", &nTets);
  for (GR_index_t iC = 0; iC < nTets; iC++) {
    GR_index_t iV0, iV1, iV2, iV3;
    fscanf(pFInFile, "%*u %u %u %u %u\n", &iV0, &iV1, &iV2, &iV3);
    assert(iV0 > 0 && iV0 <= nVerts);
    assert(iV1 > 0 && iV1 <= nVerts);
    assert(iV2 > 0 && iV2 <= nVerts);
    assert(iV3 > 0 && iV3 <= nVerts);
    Vert *pV0 = pVM->getVert(iV0 - 1);
    Vert *pV1 = pVM->getVert(iV1 - 1);
    Vert *pV2 = pVM->getVert(iV2 - 1);
    Vert *pV3 = pVM->getVert(iV3 - 1);
    bool existed = false;
    pVM->createTetCell(existed, pV0, pV1, pV2, pV3);
    assert(!existed);
  }
  fclose(pFInFile);
  ///////////////////////////////////////////////////////////////

  // Don't need to read the .neigh file.  All it can do is confirm
  // what we already know anyway.

  ///////////////////////////////////////////////////////////////
  // Next: bdry face connectivity (.face)
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.face", strBaseFileName);
  pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile) {
    vFatalError("Couldn't open output boundary face file for writing",
		"3D Tetgen mesh output");
  }

  GR_index_t nBdryTris = 0;
  fscanf(pFInFile, "%u 1\n", &nBdryTris);
  if (aiBCList != nullptr) {
    delete[] aiBCList;
  }

  struct faceStruct {
    Face* actualFace;
    Vert *pV0, *pV1, *pV2;
    int BC;
    bool
    operator<(const struct faceStruct& fs) const
    {
      return actualFace < fs.actualFace;
    }
  };
  std::multiset<struct faceStruct> facesToSplit;
  aiBCList = new int[nBdryTris];
  GR_index_t iBFSoFar = 0;
  for (GR_index_t iBF = 0; iBF < nBdryTris; iBF++) {
    GR_index_t iV0, iV1, iV2;
    int BC;
    fscanf(pFInFile, "%*u %u %u %u %d\n", &iV0, &iV1, &iV2, &BC);
    assert(iV0 > 0 && iV0 <= nVerts);
    assert(iV1 > 0 && iV1 <= nVerts);
    assert(iV2 > 0 && iV2 <= nVerts);
    Vert *pV0 = pVM->getVert(iV0 - 1);
    Vert *pV1 = pVM->getVert(iV1 - 1);
    Vert *pV2 = pVM->getVert(iV2 - 1);
    Face* face = findCommonFace(pV0, pV1, pV2, nullptr, true);
    if (face->getRightCell()->isValid() && face->getLeftCell()->isValid()) {
      // This face will need to be split.  The tricky bit is
      // going to be keeping the BC's straight.
      struct faceStruct fs;
      fs.actualFace = face;
      fs.pV0 = pV0;
      fs.pV1 = pV1;
      fs.pV2 = pV2;
      fs.BC = BC;
      facesToSplit.insert(fs);
    }
    else {
      pVM->createBFace(face);
      aiBCList[iBFSoFar] = BC;
      iBFSoFar++; // This way, when we skip one, we keep aiBCList in sync.
    }
  }

  // Now split any faces that need it.
  while (!facesToSplit.empty()) {
    assert(facesToSplit.size() % 2 == 0);
    // Grab the top two, which had better match.
    std::set<struct faceStruct>::iterator iterFirst = facesToSplit.begin();
    std::set<struct faceStruct>::iterator iterSecond = iterFirst;
    iterSecond++;
    struct faceStruct fs1 = *iterFirst;
    struct faceStruct fs2 = *iterSecond;

    assert(!(fs1 < fs2) && !(fs2 < fs1));
    // First, create an internal bdry face.
    BFace* pIBF = pVM->createBFace(fs1.actualFace);

    // Now, identify which faceStruct corresponds to which mesh face
    // incident on the new int bdry face.

    Face *pF0 = pIBF->getFace(0);
    Face *pF1 = pIBF->getFace(1);
    double normal[3], thisNormal[3];
    pF0->calcUnitNormal(normal);
    if (pF0->getRightCell() == pIBF)
      vSCALE3D(normal, -1);

    calcNormal3D(fs1.pV0->getCoords(), fs1.pV1->getCoords(),
		 fs1.pV2->getCoords(), thisNormal);
    double dot = dDOT3D(normal, thisNormal);

    int iBC0 = fs1.BC, iBC1 = fs2.BC;
    if (dot) {
      // fs1 corresponds to pF0; BC's are right.
    }
    else {
      // fs2 corresponds to pF0
      iBC0 = fs2.BC;
      iBC1 = fs1.BC;
    }
    pVM->deleteBFace(pIBF);
    BFace* bdryFace0 = pVM->createBFace(pF0);
    bdryFace0->setBdryCondition(iBC0);
    aiBCList[iBFSoFar++] = iBC0;

    BFace* bdryFace1 = pVM->createBFace(pF1);
    bdryFace1->setBdryCondition(iBC1);
    aiBCList[iBFSoFar++] = iBC1;

    facesToSplit.erase(iterFirst);
    facesToSplit.erase(iterSecond);
  }
  assert(iBFSoFar == nBdryTris);
  fclose(pFInFile);
  ///////////////////////////////////////////////////////////////
  pVM->purgeAllEntities();
}

//-----------------------------------------------------------------------------
void
writeFileTetGen(VolMesh& OutMesh, const char strBaseFileName[],
		const char /*unused*/[], const int numVars,
		const double* const solnData)
//-----------------------------------------------------------------------------
{
  FILE *pFOutFile = NULL;
  char strFileName[1024];
  OutMesh.purgeAllEntities();

  // This file format actually writes four different files, so any suffix
  // on the base file name is irrelevant.

  // Need to identify shock faces, including which of the two faces
  // incident on the bdry face  is actually on the upstream side.  Likewise,
  // need to tag verts that need to be duplicated and actually duplicate
  // them.; eventually, I'll need a map between the upstream and downstream
  // copies of a vert (those on the perimeter can map to themselves).

  // Also, in the process, the internal bdry faces can become regular bdry
  // faces with regular bdry conditions.

  // Verts that have faces with more than one BC (ends of shocks) can get
  // tagged as interior.

  // Then it should be possible to write as normal.

  ///////////////////////////////////////////////////////////////
  // First file: node locations (no tags for now).
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.node", strBaseFileName);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile) {
    vFatalError("Couldn't open output node file for writing",
		"3D Tetgen mesh output");
  }

  // TODO: Currently the data written here is very much shock fitting specific.
  // This needs to be generalized.
  GR_index_t nVerts = OutMesh.getNumVerts();
  fprintf(pFOutFile, "%u 3 %d 1\n", nVerts, numVars);
  for (GR_index_t iV = 0; iV < nVerts; iV++) {
    Vert *pV = OutMesh.getVert(iV);
    int tag = 0;
    switch (pV->getVertType())
      {
      case Vert::eBdryApex:
      case Vert::eBdryCurve:
      case Vert::eBdry:
	// Tag should be the bdry condition.  But does this really matter
	// if BC's are applied to faces?
	tag = 1;
	break;
      case Vert::ePhantom:
	tag = -1;
	break;
      case Vert::eBdryPhantom:
	tag = -2;
	break;
      case Vert::eShockUpstream:
      case Vert::eShockDownstream:
      case Vert::eShockPerimeter:
      case Vert::eShockInteraction:
	tag = pV->getVertType();
	break;
      case Vert::eInterior:
      case Vert::eInteriorFixed:
	tag = 0;
	break;
      default:
	// eBdryTwoSide, ePseudoSurface, eBBox, eUnknown, eInvalid:
	tag = 0;
	assert(0); // Should do something smarter than label them as interior.
	break;

      }
    fprintf(pFOutFile, "%8u %20.16G %20.16G %20.16G ", iV + 1, pV->x(), pV->y(),
	    pV->z());
    for (int ii = 0; ii < numVars; ii++) {
      fprintf(pFOutFile, "%20G ", solnData[iV * numVars + ii]);
    }
    fprintf(pFOutFile, "%d\n", tag);
  }
  fclose(pFOutFile);
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Next: element connectivity (.ele)
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.ele", strBaseFileName);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile) {
    vFatalError("Couldn't open output element file for writing",
		"3D Tetgen mesh output");
  }

  GR_index_t nTets = OutMesh.getNumCells();
  fprintf(pFOutFile, "%u 4 0\n", nTets);
  for (GR_index_t iC = 0; iC < nTets; iC++) {
    Cell *pC = OutMesh.getCell(iC);
    fprintf(pFOutFile, "%8u %7u %7u %7u %7u\n", iC + 1,
	    OutMesh.getVertIndex(pC->getVert(0)) + 1,
	    OutMesh.getVertIndex(pC->getVert(1)) + 1,
	    OutMesh.getVertIndex(pC->getVert(2)) + 1,
	    OutMesh.getVertIndex(pC->getVert(3)) + 1);
  }
  fclose(pFOutFile);
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Next: neighbor file (.neigh)
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.neigh", strBaseFileName);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile) {
    vFatalError("Couldn't open output neighbor file for writing",
		"3D Tetgen mesh output");
  }
  fprintf(pFOutFile, "%u 4\n", nTets);
  for (GR_index_t iC = 0; iC < nTets; iC++) {
    Cell *pC = OutMesh.getCell(iC);
    TetCell *pTC = dynamic_cast<TetCell*>(pC);
    assert(pTC);
    GR_sindex_t index[4];
    for (int ii = 0; ii < 4; ii++) {
      Face *pF = pTC->getOppositeFace(pTC->getVert(ii));
      Cell *pCOpp = pF->getOppositeCell(pTC);
      if (pCOpp->getType() == CellSkel::eTet) {
	index[ii] = OutMesh.getCellIndex(pCOpp) + 1;
      }
      else {
	index[ii] = 0;
      }
    }
    fprintf(pFOutFile, "%8u %8d %8d %8d %8d\n", iC + 1, index[0], index[1],
	    index[2], index[3]);
  }
  fclose(pFOutFile);
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Next: bdry face connectivity (.face)
  ///////////////////////////////////////////////////////////////
  snprintf(strFileName, 1024, "%s.face", strBaseFileName);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile) {
    vFatalError("Couldn't open output boundary face file for writing",
		"3D Tetgen mesh output");
  }

  GR_index_t nBdryTris = OutMesh.getNumBdryFaces();
  fprintf(pFOutFile, "%u 1\n", nBdryTris);
  for (GR_index_t iBF = 0; iBF < nBdryTris; iBF++) {
    BFace *pBF = OutMesh.getBFace(iBF);
    fprintf(pFOutFile, "%6u %7u %7u %7u %2d\n", iBF + 1,
	    OutMesh.getVertIndex(pBF->getVert(0)) + 1,
	    OutMesh.getVertIndex(pBF->getVert(1)) + 1,
	    OutMesh.getVertIndex(pBF->getVert(2)) + 1, pBF->getBdryCondition());
  }
  fclose(pFOutFile);
  ///////////////////////////////////////////////////////////////
}
