#include <queue>

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Constrain.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Geometry.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_SurfMesh.h"
#include "GR_VolMesh.h"

bool
isBdryConstrained(VolMesh* pVM, SurfMesh* pSM,
		  std::set<ConstrainEdge>& sCESubSegs,
		  std::set<ConstrainEdge>& sCEInFacets,
		  std::set<ConstrainFace>& sCFSubFacets)
{
  assert(pVM->isSimplicial());
  assert(pSM->isSimplicial());

  //   //@@ Make sure that are the vertex hints are valid in both meshes
  //   pVM->vSetAllHintFaces();
  //   pSM->vSetAllHintFaces();

  //@@ Clear the lists of missing edges and faces
  sCESubSegs.clear();
  sCEInFacets.clear();
  sCFSubFacets.clear();

  // The old code for this went through a vertex at a time, if you can
  // imagine.  Looking at the problem now, the obvious choice seems to
  // be to make sets of subsements, in-facet edges, and subfacets from
  // the surface mesh.  Then make sets of edges and faces in the volume
  // mesh.  Finally, subtract these sets to get the real list of
  // unconstrained edges and faces.
  std::set<ConstrainEdge> sCESurfSegs, sCESurfInFacets, sCEVolEdges;
  std::set<ConstrainFace> sCFSurfTris, sCFVolTris;

  for (int iEdge = pSM->getNumEdgeFaces() - 1; iEdge >= 0; iEdge--) {
    Face *pF = pSM->getFace(iEdge);
    if (pF->isDeleted())
      continue;
    Cell *pCR = pF->getRightCell();
    Cell *pCL = pF->getLeftCell();

    BdryPatch3D *pBPLeft = NULL, *pBPRight = NULL;

    if (pCL->getType() != CellSkel::eBdryEdge) {
      pBPLeft = pSM->getBdryPatch(pSM->getCellIndex(pCL));
    }
    if (pCR->getType() != CellSkel::eBdryEdge) {
      pBPRight = pSM->getBdryPatch(pSM->getCellIndex(pCR));
    }

    int iV0 = pSM->iVertIndex(pF->getVert(0));
    int iV1 = pSM->iVertIndex(pF->getVert(1));

    if (pBPLeft == pBPRight) {
      // Interior to a facet
      sCESurfInFacets.insert(ConstrainEdge(iV0, iV1));
    }
    else {
      // On a segment
      sCESurfSegs.insert(ConstrainEdge(iV0, iV1));
    }
  }

  for (int iMulti = pSM->getNumMultiEdges() - 1; iMulti >= 0; iMulti--) {
    Face *pF = pSM->getMultiEdge(iMulti);

    int iV0 = pSM->iVertIndex(pF->getVert(0));
    int iV1 = pSM->iVertIndex(pF->getVert(1));

    // Multiedges are always On a segment
    sCESurfSegs.insert(ConstrainEdge(iV0, iV1));
  }

  // Now for the triangles in the surface list.
  for (int iCell = pSM->getNumCells() - 1; iCell >= 0; iCell--) {
    Cell *pC = pSM->getCell(iCell);
    if (pC->isDeleted())
      continue;

    int iV0 = pSM->iVertIndex(pC->getVert(0));
    int iV1 = pSM->iVertIndex(pC->getVert(1));
    int iV2 = pSM->iVertIndex(pC->getVert(2));

    //     vMessage(2, "Surf tri: %3d %3d %3d\n", iV0, iV1, iV2);
    sCFSurfTris.insert(ConstrainFace(iV0, iV1, iV2));
  }

  // Build comparable sets from the volume mesh.
  for (int iFace = pVM->getNumFaces() - 1; iFace >= 0; iFace--) {
    Face *pF = pVM->getFace(iFace);
    if (pF->isDeleted())
      continue;

    int iV0 = pVM->getVertIndex(pF->getVert(0));
    int iV1 = pVM->getVertIndex(pF->getVert(1));
    int iV2 = pVM->getVertIndex(pF->getVert(2));

    sCFVolTris.insert(ConstrainFace(iV0, iV1, iV2));
    sCEVolEdges.insert(ConstrainEdge(iV0, iV1));
    sCEVolEdges.insert(ConstrainEdge(iV1, iV2));
    sCEVolEdges.insert(ConstrainEdge(iV2, iV0));
  }

  // Now use set_difference to find out what's missing.  Yes, there will
  // be lots of entries from the volume mesh that don't match anything,
  // but this construction is still only order of the size of the volume
  // mesh....
  std::insert_iterator<std::set<ConstrainEdge> > iiSubSegs(sCESubSegs,
							   sCESubSegs.begin()),
      iiInFacets(sCEInFacets, sCEInFacets.begin());
  std::set_difference(sCESurfSegs.begin(), sCESurfSegs.end(),
		      sCEVolEdges.begin(), sCEVolEdges.end(), iiSubSegs);
  std::set_difference(sCESurfInFacets.begin(), sCESurfInFacets.end(),
		      sCEVolEdges.begin(), sCEVolEdges.end(), iiInFacets);

  std::insert_iterator<std::set<ConstrainFace> > iiFaces(sCFSubFacets,
							 sCFSubFacets.begin());
  std::set_difference(sCFSurfTris.begin(), sCFSurfTris.end(),
		      sCFVolTris.begin(), sCFVolTris.end(), iiFaces);

#ifndef NDEBUG
  //   vMessage(1, "Missing %zd input edges, %zd/%d total edges, %zd/%d faces.\n",
  // 	   sCESubSegs.size(), sCESubSegs.size()+sCEInFacets.size(),
  // 	   pSM->iNumFaces(), sCFSubFacets.size(), pSM->iNumCells());

  //   for (std::set<ConstrainFace>::iterator iter = sCFSubFacets.begin();
  //        iter != sCFSubFacets.end(); ++iter) {
  //     ConstrainFace CF = *iter;
  // //     vMessage(2, "Missing tri: %3d %3d %3d\n", CF.iV0, CF.iV1, CF.iV2);
  //   }

  // Check that all locked volume faces are in fact ones that are
  // present on the boundary in the surface mesh.
  for (int ii = pVM->getNumFaces() - 1; ii >= 0; ii--) {
    Face *pF = pVM->getFace(ii);
    if (pF->isDeleted())
      continue;
    if (pF->isLocked()) {
      int iV0 = pVM->getVertIndex(pF->getVert(0));
      int iV1 = pVM->getVertIndex(pF->getVert(1));
      int iV2 = pVM->getVertIndex(pF->getVert(2));
      ConstrainFace CF(iV0, iV1, iV2);
      assert(sCFSurfTris.count(CF) == 1);
    }
  }
#endif
  // Lock all the volume faces that have been found.
  std::set<ConstrainFace> sCFFound;
  std::insert_iterator<std::set<ConstrainFace> > iiFound(sCFFound,
							 sCFFound.begin());
  std::set_difference(sCFSurfTris.begin(), sCFSurfTris.end(),
		      sCFSubFacets.begin(), sCFSubFacets.end(), iiFound);

  std::set<ConstrainFace>::iterator iter = sCFFound.begin(), iterEnd =
      sCFFound.end();
  for (; iter != iterEnd; ++iter) {
    ConstrainFace CF = *iter;
    Vert *pVA = pVM->getVert(CF.iV0);
    Vert *pVB = pVM->getVert(CF.iV1);
    Vert *pVC = pVM->getVert(CF.iV2);
    Face *pF = findCommonFace(pVA, pVB, pVC, NULL, true);
    assert(pF);
    pF->lock(); // It might already have been, but now it is for sure.
  }

  // Also lock all the surface mesh edges that have been matched, so
  // there won't be a swap to take these away.
  for (int ii = pSM->getNumFaces() - 1; ii >= 0; ii--) {
    Face *pF = pSM->getFace(ii);
    if (pF->isDeleted())
      continue;
    int iV0 = pSM->iVertIndex(pF->getVert(0));
    int iV1 = pSM->iVertIndex(pF->getVert(1));
    ConstrainEdge CE(iV0, iV1);
    if (sCESubSegs.count(CE) == 0 && sCEInFacets.count(CE) == 0) {
      //       vMessage(2, "Locking edge %d %d\n", iV0, iV1);
      pF->lock(); // It might already have been, but now it is for sure.
    }
    else {
      //       vMessage(2, "Unlocking edge %d %d\n", iV0, iV1);
      pF->unlock(); // It might already have been, but now it is for sure.
    }
  }

  //@@ Issue a report on the results
  if (sCESubSegs.empty() && sCEInFacets.empty() && sCFSubFacets.empty()) {
    logMessage(1, "All edges and faces are properly constrained.\n");
    return true;
  }
  else {
    logMessage(2,
	       "Missing %zu input edges, %zu/%u total edges, %zu/%u faces.\n",
	       sCESubSegs.size(), sCESubSegs.size() + sCEInFacets.size(),
	       pSM->getNumFaces(), sCFSubFacets.size(), pSM->getNumCells());
    return false;
  }
}

void
VolMesh::identifyBadEdgesAndFaces(const Vert* const pV0, const Vert* const pV1,
				  const Vert* const pV2,
				  std::set<BadEdge>& sBEBadEdges,
				  std::set<Face*>& spFBadFaces)
{
  // Now proceed to recover the face.  First, find a cell containing an
  // edge of the required face.
  Cell *pC = beginPipe(pV0, pV1);
  assert(pC->hasVert(pV0));
  assert(pC->hasVert(pV1));

  sBEBadEdges.clear();
  spFBadFaces.clear();

  // If the face already exists within this cell, great.
  if (pC->hasVert(pV2))
    return;

  // Walk around edge to find cell that crosses the surface face.
  int qDone;
  do {
    qDone = 1;
    for (int iF = 0; iF < 4; iF++) {
      Face *pF = pC->getFace(iF);
      if (!(pF->hasVert(pV0) && pF->hasVert(pV1)))
	continue;
      assert(pF->getType() == Face::eTriFace);
      // Walk through this face if the search vertex is on the opposite
      // side of the plane of the face from the other vertex in pC.
      int iOrient0 = checkOrient3D(pF->getVert(0), pF->getVert(1),
				   pF->getVert(2), pV2);
      int iOrient1 = checkOrient3D(
	  pF->getVert(0), pF->getVert(1), pF->getVert(2),
	  dynamic_cast<TetCell*>(pC)->getOppositeVert(pF));
      if (iOrient0 * iOrient1 == -1) { // Don't walk in case of ties.
	qDone = 0;
	pC = pF->getOppositeCell(pC);
	if (pC->getType() == Cell::eIntTriBFace) {
	  pF = (pC->getFace(0) == pF ? pC->getFace(1) : pC->getFace(0));
	  pC = pF->getOppositeCell(pC);
	  assert(pC->getType() == Cell::eTet);
	}
	assert(pC->getType() == Cell::eTet);
	break;
      }
    } // end loop over faces of a tet
  }
  while (!qDone);
  // Now we have a cell containing both pV0 and pV1 and also spanning
  // the direction from pV0 towards pV2.  Again, we may already be
  // done.
  if (pC->hasVert(pV2))
    return;

  //@@@ Identify the bad edge and faces
  Vert *pVNorth = pVInvalidVert, *pVSouth = pVInvalidVert, *pVOther =
  pVInvalidVert;
  Face *pF;
  for (int iF = 0; iF < 4; iF++) {
    // Get a face which has only one of pV0 and pV1
    pF = pC->getFace(iF);
    if (pF->hasVert(pV0) && pF->hasVert(pV1))
      continue;
    // Identify the verts at either end of the bad edge
    if (pF->getVert(0) == pV0 || pF->getVert(0) == pV1) {
      pVOther = pF->getVert(0);
      if (checkOrient3D(pV0, pV1, pV2, pF->getVert(1)) == 1) {
	pVNorth = pF->getVert(1);
	pVSouth = pF->getVert(2);
      }
      else {
	pVNorth = pF->getVert(2);
	pVSouth = pF->getVert(1);
      }
    }
    else if (pF->getVert(1) == pV0 || pF->getVert(1) == pV1) {
      pVOther = pF->getVert(1);
      if (checkOrient3D(pV0, pV1, pV2, pF->getVert(0)) == 1) {
	pVNorth = pF->getVert(0);
	pVSouth = pF->getVert(2);
      }
      else {
	pVNorth = pF->getVert(2);
	pVSouth = pF->getVert(0);
      }
    }
    else {
      pVOther = pF->getVert(2);
      assert(pF->getVert(2) == pV0 || pF->getVert(2) == pV1);
      if (checkOrient3D(pV0, pV1, pV2, pF->getVert(0)) == 1) {
	pVNorth = pF->getVert(0);
	pVSouth = pF->getVert(1);
      }
      else {
	pVNorth = pF->getVert(1);
	pVSouth = pF->getVert(0);
      }
    }
    break;
  }

  sBEBadEdges.insert(BadEdge(pVNorth, pVSouth, pVOther, pF));
  spFBadFaces.insert(pF);

  std::queue<BadEdge> qBE;
  qBE.push(BadEdge(pVNorth, pVSouth, pVOther, pF));
  //@@@ Walk around the edge, identifying all bad faces and edges
  while (!qBE.empty()) {
    // Set up starting face, cell, and other vert of starting face
    BadEdge BE = qBE.front();
    qBE.pop();
    Face *pFStart;
    pFStart = pF = BE.pF;
    pVNorth = BE.pVN;
    pVSouth = BE.pVS;
    Vert *pV = BE.pVO;
    pC = pF->getLeftCell();
    if (pC->getType() != Cell::eTet)
      pC = pF->getRightCell();
    do {
      // Grab new face incident on edge and add to bad face list
      assert(pC->getType() == Cell::eTet);
      TetCell* pTC = dynamic_cast<TetCell*>(pC);
      Face *pFPrev = pF;
      pF = pTC->getOppositeFace(pV);
      spFBadFaces.insert(pF);
      // Grab other vertex of that face
      pV = pTC->getOppositeVert(pFPrev);
      pC = pF->getOppositeCell(pC);
      // If vertex is not one of pV0, pV1, pV2, there's another bad edge
      // here. Identify it and store it.
      if (pV != pV0 && pV != pV1 && pV != pV2) {
	int iOrient = checkOrient3D(pV0, pV1, pV2, pV);
	assert(iOrient != 0);
	if (iOrient == 1) {
	  // pV is on the "north" side of the surface face
	  BadEdge BETmp(pV, pVSouth, pVNorth, pF);
	  if (sBEBadEdges.count(BETmp) == 0) {
	    sBEBadEdges.insert(BETmp);
	    qBE.push(BETmp);
	  }
	}
	else {
	  // pV is on the "south" side of the surface face
	  BadEdge BETmp(pVNorth, pV, pVSouth, pF);
	  if (sBEBadEdges.count(BETmp) == 0) {
	    sBEBadEdges.insert(BETmp);
	    qBE.push(BETmp);
	  }
	}
      } // Done adding new bad edge
    }
    while (pF != pFStart);
  } // Done adding all bad edges
  logMessage(4, "Found %zd edges intersecting desired face.\n",
	     sBEBadEdges.size());
  logMessage(4, "Found %zd faces intersecting desired face.\n",
	     spFBadFaces.size());
}

//@ Recover a surface face that's missing from the volume mesh
// This is done by swapping if possible and by point insertion if necessary

int
VolMesh::recoverFace(const int aiVerts[3], int * const piPointsInserted,
		     const bool qAllowEdgeSwap)
{

  //   The following assertion chokes when called from
  //   VolMesh::VolMesh(Bdry3D&)
  //   assert(qValid());
  logMessage(3, "%s %d, %d, and %d.\n",
	     "Trying to recover surface face composed of vertices ", aiVerts[0],
	     aiVerts[1], aiVerts[2]);
  //@@ Find all faces and edges violating the given surface face.
  Vert *pV0 = getVert(aiVerts[0]);
  Vert *pV1 = getVert(aiVerts[1]);
  Vert *pV2 = getVert(aiVerts[2]);

  assert(pV0);
  assert(pV1);
  assert(pV2);

  return recoverFace(pV0, pV1, pV2, piPointsInserted, qAllowEdgeSwap);

}
int
VolMesh::recoverFace(const Vert * pV0, const Vert * pV1, const Vert * pV2,
		     int * const piPointsInserted, const bool qAllowEdgeSwap,
		     int iRecurDepth)
{

  if (iRecurDepth > 5)
    return 0;
  int iRetVal = 0;
  // First, make sure that all the edges of the surface face are
  // present.  If they aren't, this means that they are interior edges
  // to a surface facet (and so weren't recovered earlier when doing
  // bdry subsegments); and that there is a vertex encroaching on them
  // (which makes it so that the Delaunay tetrahedralization doesn't
  // match the constrained DT).  If the edge is missing, recover it
  // without point insertion (this should always be possible, but note
  // that there's still a check to be sure that it works).
  int iEdgeAOK = recoverEdge(pV0, pV1);
  int iEdgeBOK = recoverEdge(pV1, pV2);
  int iEdgeCOK = recoverEdge(pV2, pV0);
  if ((iEdgeAOK != 1) || (iEdgeBOK != 1) || (iEdgeCOK != 1)) {
    //		printf("can't recover edges\n");
    //		if (iEdgeAOK != 1){
    //			printf("edge 1\n");
    //			pV0->printVertInfo();
    //			pV1->printVertInfo();
    //		}
    //		if (iEdgeBOK != 1){
    //			printf("edge 2\n");
    //			pV1->printVertInfo();
    //			pV2->printVertInfo();
    //		}
    //		if (iEdgeCOK != 1){
    //			printf("edge 3\n");
    //			pV2->printVertInfo();
    //			pV0->printVertInfo();
    //		}
    return iRetVal;
  }

  std::set<BadEdge> sBEBadEdges;
  std::set<Face*> spFBadFaces;
  identifyBadEdgesAndFaces(pV0, pV1, pV2, sBEBadEdges, spFBadFaces);

  //@@ For each face in the list, swap it if possible.
  // Swap only if the face is of type T32 and all of its equatorial
  // points are on the same side of or in the plane of the surface face.
  bool qChange, qGlobalChange = false;
  do {
    qChange = false;
    std::queue<Face*> qpFBad;
    for (std::set<Face*>::iterator iterF = spFBadFaces.begin();
	iterF != spFBadFaces.end(); ++iterF) {
      qpFBad.push(*iterF);
    }

    while (!qpFBad.empty()) {
      Face* pF = qpFBad.front();
      qpFBad.pop();
      assert(pF->getType() == Face::eTriFace);
      TriFace *pTF = dynamic_cast<TriFace*>(pF);
      if (pTF->isDeleted())
	continue;
      Vert *pVVertA, *pVVertB, *pVVertC, *pVVertD, *pVVertE;
      TetCell* apTCTets[4];
      int iNTets;
      Vert *pVJunk0, *pVJunk1, *pVJunk2;
      eFaceCat eFC = pTF->categorizeFace(pVVertA, pVVertB, pVVertC, pVVertD,
					 pVVertE, apTCTets, iNTets, pVJunk0,
					 pVJunk1, pVJunk2);

      if (eFC == eT32) {
	assert(iNTets == 3);
	Face *apFDel[3];
	apFDel[0] = findCommonFace(apTCTets[0], apTCTets[1]);
	apFDel[1] = findCommonFace(apTCTets[1], apTCTets[2]);
	apFDel[2] = findCommonFace(apTCTets[2], apTCTets[0]);
	assert(apFDel[0]->isValid());
	assert(apFDel[1]->isValid());
	assert(apFDel[2]->isValid());
	int iOrientD = checkOrient3D(pV0, pV1, pV2, pVVertD);
	int iOrientE = checkOrient3D(pV0, pV1, pV2, pVVertE);
	if (iOrientD * iOrientE != -1)
	  continue;
	assert(pTF == apFDel[0] || pTF == apFDel[1] || pTF == apFDel[2]);

	qChange = true;
	qGlobalChange = true;
#ifndef NDEBUG
	int iSwaps =
#endif
	    reconfTet32_deprecated(apTCTets, pTF, pVVertA, pVVertB, pVVertC,
				   pVVertD, pVVertE);
	assert(iSwaps == 1);

	iRetVal++;
	logMessage(3, "Swapped 3-to-2\n");
	// Check each of the three faces to be sure that either
	// they've been deleted or all of their verts lie on the
	// same side of the target face (including ties).
	for (int ii = 0; ii < 3; ii++) {
	  Face *pFTest = apFDel[ii];
	  if (pFTest->isDeleted()) {
	    spFBadFaces.erase(pFTest);
	  }
	  else {
	    Vert *pVTest0 = pFTest->getVert(0);
	    Vert *pVTest1 = pFTest->getVert(1);
	    Vert *pVTest2 = pFTest->getVert(2);
	    int iOrient0 = checkOrient3D(pV0, pV1, pV2, pVTest0);
	    int iOrient1 = checkOrient3D(pV0, pV1, pV2, pVTest1);
	    int iOrient2 = checkOrient3D(pV0, pV1, pV2, pVTest2);
	    if ((iOrient0 * iOrient1 != -1) && (iOrient1 * iOrient2 != -1)
		&& (iOrient2 * iOrient0 != -1)) {
	      spFBadFaces.erase(pFTest);
	    }
	  }
	}
      } // T32 face
    } // End loop over all bad faces
    identifyBadEdgesAndFaces(pV0, pV1, pV2, sBEBadEdges, spFBadFaces);
    // Do this until no more progress can be made by simple swapping or
    // until all offending edges are removed
  }
  while (qChange && !spFBadFaces.empty());

  // If no swapping has been done at all, then we can do edge swapping
  // based on the original bad edge data.  Also, don't allow too much
  // explosion of bad edges before we give up entirely.

  // This implementation is surprisingly effective, considering how crude it is.
  //   bool qGotOne = false;
  if (!spFBadFaces.empty() && !qGlobalChange && sBEBadEdges.size() < 40
      && qAllowEdgeSwap) {
    logMessage(3, "Trying edge swapping...\n");
    bool qSuccess = false;
    for (std::set<BadEdge>::iterator iter = sBEBadEdges.begin();
	iter != sBEBadEdges.end() && !qSuccess; ++iter) {
      BadEdge BE = *iter;
      qSuccess = (0 != edgeSwap3D(BE.pF, BE.pVN, BE.pVS, BE.pVO, true));
      logMessage(3, "Face: %d  Verts: %d %d %d   %s\n", getFaceIndex(BE.pF),
		 getVertIndex(BE.pVN), getVertIndex(BE.pVS),
		 getVertIndex(BE.pVO), qSuccess ? "success" : "failure");
    }
    // If an edge has been removed, then recurse.
    if (qSuccess) {
      int iRetRecur = recoverFace(pV0, pV1, pV2, piPointsInserted,
				  qAllowEdgeSwap, iRecurDepth + 1);
      if (iRetRecur)
	return 1000 + iRetVal + iRetRecur;
      else
	return 0;
    }
  }

  if (spFBadFaces.empty()) {
    logMessage(3, "All bad faces removed by swapping\n");
  }
  else {
    logMessage(3, "After swapping, %d bad faces and %d bad edges remain\n",
	       static_cast<int>(spFBadFaces.size()),
	       static_cast<int>(sBEBadEdges.size()));

  }

  //@@ If swapping fails, insert a point to split an offending face.
  // The bad face should be split at its circumcenter, so long as this
  // doesn't violate encroachment rules.  Not yet implemented,
  // obviously.

  return iRetVal;
}

void
VolMesh::markCleanNeighbors(Cell* pC)
{
  int iReg = pC->getRegion();
  for (int iF = 0; iF < pC->getNumFaces(); iF++) {
    Face *pF = pC->getFace(iF);
    // If the face has been deleted, then we need not do anything.
    if (pF->isValid() && !pF->isDeleted()) {
      Cell *pCCand = pF->getOppositeCell(pC);
      assert(pCCand->isValid());

      // Be sure never to deal with this face again
      if (pCCand->getRegion() == iReg)
	continue;
      else if (pCCand->getRegion() == iInvalidRegion) {
	// If the opposite cell hasn't been marked, mark it the same as
	// this one.
	pCCand->setRegion(iReg);
	markCleanNeighbors(pCCand);
      } // Candidate cell needs to be tagged.
    } // Face is OK
  } // Loop over faces
}

// Find the cell incident on pV0 that lies astride the direction to pV1.
Cell*
VolMesh::beginPipe(const Vert* const pV0, const Vert* const pV1) const
{

  Face *pF = const_cast<Face*>(pV0->getFace(0));
  Cell *pC = pF->getLeftCell();
  if (pC->getType() == Cell::eTriBFace || pC->getType() == Cell::eIntTriBFace) {
    pC = pF->getRightCell();
    assert(pC->getType() == Cell::eTet);
  }
  bool qFoundIt = false;
#ifndef NDEBUG
  int iTrials = 0;
#endif
  while (!qFoundIt) {
    qFoundIt = true;
#ifndef NDEBUG
    iTrials++;
    assert(iTrials < 1000);
#endif
    Face *pFOpp = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV0);
    for (int iF = 0; iF < 4; iF++) {
      Face *pFTrial = pC->getFace(iF);
      if (pFTrial == pFOpp)
	continue;
      int iSign = pFTrial->getLeftCell() == pC ? -1 : 1;
      if (checkOrient3D(pFTrial->getVert(0), pFTrial->getVert(1),
			pFTrial->getVert(2), pV1) == -iSign) {
	// Incorrect orientation; walk through this face and try again
	pC = pFTrial->getOppositeCell(pC);
	if (pC->getType() == Cell::eIntTriBFace) {
	  pFTrial =
	      (pC->getFace(0) == pFTrial ? pC->getFace(1) : pC->getFace(0));
	  pC = pFTrial->getOppositeCell(pC);
	  assert(pC->getType() == Cell::eTet);
	}
	qFoundIt = false;
	break;
      }
    }
  }
  assert(pC->getType() == Cell::eTet);
  return pC;
}

// Find the face (other than pFLast) through which the line from pV0 to
// pV1 crosses the surface of Cell pC.
static void
vWhichFace(const Cell* const pC, const Face* const pFLast,
	   const Vert* const pV0, const Vert* const pV1, const Face*& pF,
	   bool& qHitEdge)
{
  qHitEdge = false;
  for (int iF = 0; iF < 4; iF++) {
    pF = pC->getFace(iF);
    // Don't give back the same face that was previously found
    if (pF == pFLast)
      continue;
    // Don't give back a face containing one of the endpoints
    if (pF->hasVert(pV0) || pF->hasVert(pV1))
      continue;
    const Vert *pVA = pF->getVert(0);
    const Vert *pVB = pF->getVert(1);
    const Vert *pVC = pF->getVert(2);
    // If pV0 and pV1 aren't on opposite sides of the face, then you'll
    // have to try another.
    if (checkOrient3D(pVA, pVB, pVC, pV0) != -checkOrient3D(pVA, pVB, pVC, pV1))
      continue;
    int iOrientAB = checkOrient3D(pVA, pVB, pV0, pV1);
    if (iOrientAB == 0) {
      qHitEdge = true;
      return;
    }  // Hit an edge
    int iOrientBC = checkOrient3D(pVB, pVC, pV0, pV1);
    if (iOrientBC == 0) {
      qHitEdge = true;
      return;
    }  // Hit an edge
    int iOrientCA = checkOrient3D(pVC, pVA, pV0, pV1);
    if (iOrientCA == 0) {
      qHitEdge = true;
      return;
    }  // Hit an edge
    if (iOrientAB == iOrientBC && iOrientAB == iOrientCA)
      return;
  }
}

//@ Work your way through a pipe trying to remove faces.
int
VolMesh::clearPipe(const Vert* const pVInit, const Vert* const pVEnd,
		   bool& qHitEdge, Face* apFSkipped[], int& iNSkip)
// Return 1 if the pipe is cleared, 0 if progress is made but the
// pipe still has internal faces, and -1 if no changes were made.
{
  // A new version of this routine, including edge swapping internally.

  // Find the first cell at the pVInit end of the pipe from pVInit to pVEnd
  Cell* pCStart0 = beginPipe(pVInit, pVEnd);

  logMessage(3, "Pipe begins: %5u  Pipe ends: %5u\n", getVertIndex(pVInit),
	     getVertIndex(pVEnd));

  // Walk from pVInit towards pVEnd until either (a) you hit an edge or (b)
  // you get to pVEnd.  Keep track of how many faces you had to skip
  // over. Any time you swap a face, you must begin the ream again,
  // because the identities of cells and faces change with swapping, and
  // you could end up holding an invalid cell or face pointer, and
  // cratering the code.

  // Find which face the edge (pVInit,pVEnd) passes through
  Face *pFCand = pFInvalidFace, *pFLast = pFInvalidFace;
  Cell *pC = pCStart0;
  iNSkip = 0;
  TetCell* apTCTets[4];
  int iNTets;
  bool qChange = false, qReset = false, qEndOfPipe = false;
  int iCount = 0;
  eFaceCat ePrevFC = eOther;
  Vert *pVPrevD = pVInvalidVert, *pVPrevE = pVInvalidVert;

  do {
    if (pC->hasVert(pVEnd)) {
      logMessage(3, "  Pipe already reamed.\n");
      return 1;
    }
    {
      const Face* pFTmp = pFInvalidFace;
      vWhichFace(pC, pFLast, pVInit, pVEnd, pFTmp, qHitEdge);
      pFCand = const_cast<Face*>(pFTmp);
    }

    if (qHitEdge) {
      // Try to swap away this edge in the volume mesh.  First
      // identify the edge.
      Vert *pVA = pFCand->getVert(0);
      Vert *pVB = pFCand->getVert(1);
      Vert *pVC = pFCand->getVert(2);

      Vert *pVNorth, *pVSouth, *pVOther;

      if (checkOrient3D(pVA, pVB, pVInit, pVEnd) == 0) {
	pVNorth = pVA;
	pVSouth = pVB;
	pVOther = pVC;
      }
      else if (checkOrient3D(pVB, pVC, pVInit, pVEnd) == 0) {
	pVNorth = pVB;
	pVSouth = pVC;
	pVOther = pVA;
      }
      else {
	assert(checkOrient3D(pVC, pVA, pVInit, pVEnd) == 0);
	pVNorth = pVC;
	pVSouth = pVA;
	pVOther = pVB;
      }
      int iSwaps = edgeSwap3D(pFCand, pVNorth, pVSouth, pVOther, true, pVInit,
			      pVEnd);
      if (iSwaps > 0) {
	// Success!  Well, at least we got rid of the offending edge.
	// So restart through the pipe.
	logMessage(3, "  Successfully flipped an edge.\n");
	qReset = qChange = true;
	qHitEdge = false;
      }
      else {
	qHitEdge = true;
	logMessage(3, "   Failed to flip edge.\n");
      }
    }
    else { // The normal case, where you haven't hit an edge
      assert(pFCand->getType() == Face::eTriFace);
      TriFace *pTF = dynamic_cast<TriFace*>(pFCand);
      Vert *pVPivot0, *pVPivot1, *pVOther;
      Vert *pVVertA, *pVVertB, *pVVertC, *pVVertD, *pVVertE;
      eFaceCat eFC = pTF->categorizeFace(pVVertA, pVVertB, pVVertC, pVVertD,
					 pVVertE, apTCTets, iNTets, pVPivot0,
					 pVPivot1, pVOther);
      if (eFC == eT23 || eFC == eT32 || eFC == eT44 || eFC == eN32) {
	logMessage(3, " FaceCat = %d %8u %8u %8u %8u %8u\n", eFC,
		   getVertIndex(pVVertA), getVertIndex(pVVertB),
		   getVertIndex(pVVertC), getVertIndex(pVVertD),
		   getVertIndex(pVVertE));
      }
      if ((eFC == eT44)
	  && (pVInit == pVVertD || pVEnd == pVVertE || pVInit == pVVertE
	      || pVEnd == pVVertD)) {
	// Do T44 swaps whenever the T44 face is first or last in the
	// pipe.
	logMessage(3, "    Swapped T44\n");
	reconfTet44_deprecated(apTCTets, pTF, pVVertA, pVVertB, pVVertC,
			       pVVertD, pVVertE);
	qChange = qReset = true;
      }
      else if (eFC == eT23) {
	// Do T23 swaps at the beginning or end of the pipe
	// OR
	// If the previous and next faces in the pipe share an edge.
	bool qSwap = (pVInit == pVVertD || pVInit == pVVertE || pVEnd == pVVertD
	    || pVEnd == pVVertE);
	if (!qSwap) {
	  // For interior faces of the pipe, check to see whether the
	  // preceding and following pipe faces share an edge.
	  assert(pFLast->isValid());
	  Cell *pCTemp = pFCand->getOppositeCell(pC);
	  assert(!pCTemp->hasVert(pVEnd));
	  const Face *pFNext;
	  bool qHitTemp = false;
	  vWhichFace(pCTemp, pFCand, pVInit, pVEnd, pFNext, qHitTemp);
	  assert(pFNext->isValid());
	  if (!qHitTemp) {
	    int iV, iMatch = 0;
	    for (iV = 0; iV < 3; iV++)
	      if (pFNext->hasVert(pFLast->getVert(iV)))
		iMatch++;
	    assert(iMatch == 1 || iMatch == 2);
	    if (iMatch == 2)
	      qSwap = true;
	  }
	}
	if (qSwap) {
	  logMessage(3, "    Swapped T23\n");
	  reconfTet23_deprecated(apTCTets, pTF, pVVertA, pVVertB, pVVertC,
				 pVVertD, pVVertE);
	  qChange = qReset = true;
	}
      } // Done with T23 case
      else if (eFC == eT32) {
	// Swap if this is the first face at one end or the other.
	// OR
	// If the previous face was also T32 about the same edge.
	// OR
	// If the previous and subsequent faces share an edge.
	bool qSwap = (pVInit == pVVertA || pVInit == pVVertB
	    || pVInit == pVVertC || pVEnd == pVVertA || pVEnd == pVVertB
	    || pVEnd == pVVertC)
	    || (ePrevFC == eT32
		&& ((pVPrevD == pVVertD && pVPrevE == pVVertE)
		    || (pVPrevD == pVVertE && pVPrevE == pVVertD)));
	if (!qSwap) {
	  // For interior faces of the pipe, check to see whether the
	  // preceding and following pipe faces share an edge.
	  assert(pFLast->isValid());
	  Cell *pCTemp = pFCand->getOppositeCell(pC);
	  assert(!pCTemp->hasVert(pVEnd));
	  const Face *pFNext;
	  bool qHitTemp = false;
	  vWhichFace(pCTemp, pFCand, pVInit, pVEnd, pFNext, qHitTemp);
	  assert(pFNext->isValid());
	  if (!qHitTemp) {
	    int iV, iMatch = 0;
	    for (iV = 0; iV < 3; iV++)
	      if (pFNext->hasVert(pFLast->getVert(iV)))
		iMatch++;
	    assert(iMatch == 1 || iMatch == 2);
	    if (iMatch == 2)
	      qSwap = true;
	  }
	} // Done checking for previous/next faces sharing an edge
	if (qSwap) {
	  logMessage(3, "    Swapped T32\n");
	  reconfTet32_deprecated(apTCTets, pTF, pVVertA, pVVertB, pVVertC,
				 pVVertD, pVVertE);
	  qChange = qReset = true;
	}
      } // Done with T32 case
      else if (eFC == eN32) {
	// Swap if this is the first face at one end or the other.
	// OR
	// If the previous and subsequent faces share an edge.
	bool qSwap = (pVInit == pVVertD || pVInit == pVVertE || pVEnd == pVVertD
	    || pVEnd == pVVertE);
	// 	if (!qSwap) {
	// 	  // For interior faces of the pipe, check to see whether the
	// 	  // preceding and following pipe faces share an edge.
	// 	  assert(pFLast->qValid());
	// 	  Cell *pCTemp = pFCand->pCCellOpposite(pC);
	// 	  assert(!pCTemp->qHasVert(pVEnd));
	// 	  const Face *pFNext;
	// 	  bool qHitTemp = false;
	// 	  vWhichFace(pCTemp, pFCand, pVInit, pVEnd, pFNext, qHitTemp);
	// 	  assert(pFNext->qValid());
	// 	  if (!qHitTemp) {
	// 	    int iV, iMatch = 0;
	// 	    for (iV = 0; iV < 3; iV++)
	// 	      if (pFNext->qHasVert(pFLast->pVVert(iV))) iMatch++;
	// 	    assert (iMatch == 1 || iMatch == 2);
	// 	    if (iMatch == 2)
	// 	      qSwap = true;
	// 	  }
	// 	} // Done checking for previous/next faces sharing an edge
	if (qSwap) {
	  // Target edge connects verts D and E; if that can be managed,
	  // then all is well.  Otherwise, forget it.
	  assert(pVVertD != pVPivot0 && pVVertE != pVPivot0);
	  assert(pVVertD != pVPivot1 && pVVertE != pVPivot1);
	  assert(pVVertD != pVOther && pVVertE != pVOther);
	  int iRetVal = edgeSwap3D(pTF, pVPivot0, pVPivot1, pVOther, false,
				   pVVertD, pVVertE);
	  if (iRetVal != 0) {
	    logMessage(3, "    Swapped N32\n");
	    iCount += iRetVal - 1; // It'll be incremented later...
	    qChange = qReset = true;
	  }
	}
      } // Done with N32 case
      else {
	qReset = false;
      }

      ePrevFC = eFC;
      pVPrevD = pVVertD;
      pVPrevE = pVVertE;
    }
    if (qReset) {
      pC = beginPipe(pVInit, pVEnd);
      pFLast = pFInvalidFace;
      qReset = false;
      iNSkip = 0;
      iCount++;
    }
    else {
      // No swapping, just go on
      apFSkipped[iNSkip] = pFCand;
      iNSkip = iNSkip + 1;
      // Grab the next cell in the pipe
      pFLast = pFCand;
      pC = pFCand->getOppositeCell(pC);
      if (pC->getType() == Cell::eIntTriBFace) {
	pFCand = (pC->getFace(0) == pFCand ? pC->getFace(1) : pC->getFace(0));
	pC = pFCand->getOppositeCell(pC);
      }
    }
    // If you made it through and no faces were skipped, you're done

    qEndOfPipe = (pC->hasVert(pVEnd) || pC->getType() == Cell::eTriBFace);
    if (qEndOfPipe && iNSkip == 0) {
      logMessage(3, "   Successfully reamed pipe; %d swaps.\n", iCount);
      return 1;
    }
  }
  while ((!qEndOfPipe) && (iCount < 10) && (iNSkip < 10));
  if (qChange)
    return 0;
  else {
    logMessage(4, "   Failed to ream pipe; %d swaps.\n", iCount);
    return -1;
  }
}

//@ Delete all entities lieing outside the surface mesh
void
VolMesh::deleteExternal(const SurfMesh* pSM, const int iNumOrigVerts,
			const int iNumBBoxVerts)
{
  assert(isSimplicial());
  assert(pSM->isSimplicial());
  GR_index_t iCell;
  // Clear the domain tag for each cell in the volume mesh.  At this
  // point, cells are being tagged only as internal or external;
  // multiple subdomains will come later.
  for (iCell = 0; iCell < getNumCells(); iCell++)
    getCell(iCell)->setRegion(iInvalidRegion);

  struct sDetach {
    Cell *pCOutside;
    Face *pFBdry;
    BdryPatch3D *pBP;
  };
  struct sDetach asBdryData[pSM->getNumCells()];

  //@@ Identify (at least most of) the internal floating faces.
  // This algorithm relies on working from bdry edges in the surface
  // mesh, which are hanging edges in the interior of the domain.  So
  // this code won't be sufficient if there's a surface bounded
  // completely by non-manifold edges.

  // First, set all surface cells as region 1; those that are not part
  // of the boundary of the domain will be marked as iOutsideRegion ---
  // a misnomer, but it makes it possible to propagate markings using
  // existing code.
  for (iCell = 0; iCell < pSM->getNumCells(); iCell++) {
    pSM->getCell(iCell)->setRegion(iInvalidRegion);
  }

  GR_index_t iNBFaces = pSM->getNumBdryFaces();
  if (iNBFaces != 0) {
    for (GR_index_t iBF = 0; iBF < iNBFaces; iBF++) {
      BFace *pBF = pSM->getBFace(iBF);
      Face *pF = pBF->getFace();
      Cell *pC = pF->getLeftCell();
      if (!pC->isValid())
	pC = pF->getRightCell();
      assert(pC->isValid());
      pC->setRegion(iOutsideRegion);
    }

    // For each marked cell in the volume mesh, recursively mark unmarked
    // neighbors
    for (iCell = 0; iCell < pSM->getNumCells(); iCell++) {
      Cell *pC = pSM->getCell(iCell);
      if (pC->getRegion() != iInvalidRegion)
	markCleanNeighbors(pC);
    }
  }

  // Now the ones you -don't- want to use for finding the exterior of
  // the volume mesh are marked as iOutside, while other are iInvalid.

  //@@ For each face in the surface mesh, find the matching volume face
  for (iCell = 0; iCell < pSM->getNumCells(); iCell++) {
    Cell *pCS = pSM->getCell(iCell);
    // If the cell's deleted or it's a hanging face in the interior,
    // don't do anything with this one.
    if (pCS->isDeleted() || pCS->getRegion() == iOutsideRegion) {
      asBdryData[iCell].pFBdry = pFInvalidFace;
      continue;
    }
    // Note: These are triangular surface mesh cells
    // Find the identities of the vertices in the volume mesh
    Vert* pVV0 = getVert(pSM->iVertIndex(pCS->getVert(0)));
    Vert* pVV1 = getVert(pSM->iVertIndex(pCS->getVert(1)));
    Vert* pVV2 = getVert(pSM->iVertIndex(pCS->getVert(2)));

    //@@@ Find the face that contains all of these verts
    Face *pF = findCommonFace(pVV0, pVV1, pVV2, pVInvalidVert, true);

    //     // First, find an edge containing two of them and the info needed to
    //     // walk around that edge
    //     Cell *pC = pCBeginPipe(pVV0, pVV1);
    //     assert(pC->qHasVert(pVV0));
    //     assert(pC->qHasVert(pVV1));
    //     Face *pF = pFInvalidFace;
    //     Vert *pV = pVInvalidVert;
    //     for (int iF = 0; iF < 4; iF++) {
    //       pF = pC->pFFace(iF);
    //       if (pF->qHasVert(pVV0) && pF->qHasVert(pVV1)) break;
    //     }
    //     for (int iV = 0; iV < 3; iV++) {
    //       pV = pF->pVVert(iV);
    //       if (pV != pVV0 && pV != pVV1) break;
    //     }

    //     // Walk around the edge until you hit the correct face
    //     Face *pFStart = pFInvalidFace;
    // #ifndef NDEBUG
    //     bool qTurned = false;
    // #endif
    //     while (!pF->qHasVert(pVV2) && pF != pFStart) {
    //       if (pFStart == pFInvalidFace) pFStart = pF;
    //       assert(pC->eType() == Cell::eTet);
    //       TetCell* pTC = dynamic_cast<TetCell*>(pC);
    //       Face *pFPrev = pF;
    //       pF = pTC->pFFaceOpposite(pV);
    //       pV = pTC->pVVertOpposite(pFPrev);
    //       pC = pF->pCCellOpposite(pC);
    //       if (pC->eType() == Cell::eTriBFace && !pF->qHasVert(pVV2)) {
    // 	// If you've hit the boundary without finding the vertex you
    // 	// were after, reverse direction and go back the other way.
    // #ifndef NDEBUG
    // 	assert(!qTurned);
    // 	qTurned = true;
    // #endif
    // 	pFPrev = pFStart = pF;
    // 	pF = pTC->pFFaceOpposite(pV);
    // 	pV = pTC->pVVertOpposite(pFPrev);
    // 	pC = pF->pCCellOpposite(pTC);
    //       }
    //     }

    assert(pF->hasVert(pVV0));
    assert(pF->hasVert(pVV1));
    assert(pF->hasVert(pVV2));

    // Is the cyclic ordering for the volume face and the surface cell
    // the same, or opposite?
    double adNormVol[3];
    pF->calcNormal(adNormVol);
    NORMALIZE3D(adNormVol);
    assert(pCS->getType() == Cell::eTriCell);
    double adVec01[] = adDIFF3D(pCS->getVert(1)->getCoords(),
	pCS->getVert(0)->getCoords());
      double adVec02[] = adDIFF3D(pCS->getVert(2)->getCoords(),
	  pCS->getVert(0)->getCoords());
      double adNormSurf[3];
      vCROSS3D(adVec01, adVec02, adNormSurf);
      NORMALIZE3D(adNormSurf);

      double adNormPatch[3];
      BdryPatch3D *pBP3D = pSM->getBdryPatch(iCell);
      pBP3D->vUnitNormal(pCS->getVert(0)->getCoords(), adNormPatch);

      double dPatchDot = -dDOT3D(adNormVol, adNormPatch);
      double dSurfDot = -dDOT3D(adNormVol, adNormSurf);

      // Mark the interior and exterior cells adjacent to this face
	int iRightHandPatch = iFuzzyComp(dPatchDot,0.);
	int iRightHandSurf = iFuzzyComp(dSurfDot,0.);
	assert(abs(iRightHandPatch) == 1);
	assert(abs(iRightHandSurf ) == 1);

      // These next few lines are relevant eventually for getting
      // multi-region support up and going again.
	int iRightRegion = pBP3D->iRightRegion();
	int iLeftRegion = pBP3D->iLeftRegion ();

	bool qIsInteriorBoundary = false;
	if (iRightRegion != iOutsideRegion &&
	    iLeftRegion != iOutsideRegion)
	qIsInteriorBoundary = true;

	asBdryData[iCell].pFBdry = pF;
	asBdryData[iCell].pBP = pBP3D;

      // Topology-based tagging isn't working as well as I'd like, so I'm
      // going to try geometry-based tagging.

	double adFaceCent[3],
    adCellCent[3];
    pF->calcCentroid(adFaceCent);
    pF->getLeftCell()->calcCentroid(adCellCent);
    double adLeftVec[] = adDIFF3D(adCellCent, adFaceCent);
    NORMALIZE3D(adLeftVec);
    double dDot = dDOT3D(adLeftVec, adNormSurf);
    if (iFuzzyComp(dDot, 0) == 0) {
      // Collapsed tet; try the other side.
      pF->getRightCell()->calcCentroid(adCellCent);
      double adRightVec[] = adDIFF3D(adCellCent, adFaceCent);
      NORMALIZE3D(adRightVec);
      dDot = -dDOT3D(adRightVec, adNormSurf);
    }
    if (qIsInteriorBoundary) {
      // For now, use the old topology-based scheme, which -seemed- to
      // work okay.
      assert(iRightRegion != iOutsideRegion);
      assert(iLeftRegion != iOutsideRegion);
      asBdryData[iCell].pCOutside = pCInvalidCell;
      if (iRightHandPatch == 1) {
	pF->getRightCell()->setRegion(iRightRegion);
	pF->getLeftCell()->setRegion(iLeftRegion);
      }
      else {
	pF->getRightCell()->setRegion(iLeftRegion);
	pF->getLeftCell()->setRegion(iRightRegion);
      }
    }
    else {
      // If dDot < 0, then for sure this cell lies inside, because the
      // surface mesh is oriented this way on purpose.
      assert(iFuzzyComp(dDot, 0) != 0);
      int iIntRegion = (
	  (iRightRegion == iOutsideRegion) ? iLeftRegion : iRightRegion);
      logMessage(4, "Face: %3u.  Verts %3u, %3u, %3u. Surf: %s.  Patch %s.\n",
		 getFaceIndex(pF), getVertIndex(pF->getVert(0)),
		 getVertIndex(pF->getVert(1)), getVertIndex(pF->getVert(2)),
		 iRightHandSurf == 1 ? "y" : "n",
		 iRightHandPatch == 1 ? "y" : "n");
      logMessage(4, "  Vert %3u @ (%10.6G, %10.6G, %10.6G)\n",
		 getVertIndex(pF->getVert(0)), pF->getVert(0)->x(),
		 pF->getVert(0)->y(), pF->getVert(0)->z());
      logMessage(4, "  Vert %3u @ (%10.6G, %10.6G, %10.6G)\n",
		 getVertIndex(pF->getVert(1)), pF->getVert(1)->x(),
		 pF->getVert(1)->y(), pF->getVert(1)->z());
      logMessage(4, "  Vert %3u @ (%10.6G, %10.6G, %10.6G)\n",
		 getVertIndex(pF->getVert(2)), pF->getVert(2)->x(),
		 pF->getVert(2)->y(), pF->getVert(2)->z());
      logMessage(4, "Surface mesh normal: (%10.6G, %10.6G, %10.6G)\n",
		 adNormSurf[0], adNormSurf[1], adNormSurf[2]);
      logMessage(4,
		 "Volume mesh normal:  (%10.6G, %10.6G, %10.6G) Dot: %8.3G\n",
		 adNormVol[0], adNormVol[1], adNormVol[2],
		 (dDOT3D(adNormSurf, adNormVol)));
      logMessage(4,
		 "Patch normal:        (%10.6G, %10.6G, %10.6G) Dot: %8.3G\n",
		 adNormPatch[0], adNormPatch[1], adNormPatch[2],
		 (dDOT3D(adNormSurf, adNormPatch)));
      if (dDot < 0) {
	Cell *pC = pF->getRightCell();
	logMessage(4, "Detaching cell %3u (<). Verts %3u, %3u, %3u, %3u.\n",
		   getCellIndex(pC), getVertIndex(pC->getVert(0)),
		   getVertIndex(pC->getVert(1)), getVertIndex(pC->getVert(2)),
		   getVertIndex(pC->getVert(3)));
	asBdryData[iCell].pCOutside = pF->getRightCell();
	pF->getRightCell()->setRegion(iOutsideRegion);
	pF->getLeftCell()->setRegion(iIntRegion);
      }
      else {
	Cell *pC = pF->getLeftCell();
	logMessage(4, "Detaching cell %3u (>). Verts %3u, %3u, %3u, %3u.\n",
		   getCellIndex(pC), getVertIndex(pC->getVert(0)),
		   getVertIndex(pC->getVert(1)), getVertIndex(pC->getVert(2)),
		   getVertIndex(pC->getVert(3)));
	asBdryData[iCell].pCOutside = pF->getLeftCell();
	pF->getLeftCell()->setRegion(iOutsideRegion);
	pF->getRightCell()->setRegion(iIntRegion);
      }
    }

    // Old topology-based tagging; fails mysteriously for at least
    // the hex-multi case; more if you try to fix it. ;-)
    //     // Tag cells as inside/outside based on the patch orientation
    //     if (iRightHandPatch == 1) {
    //       pF->pCCellRight()->vSetRegion(iRightRegion);
    //       pF->pCCellLeft() ->vSetRegion(iLeftRegion);
    //     }
    //     else {
    //       pF->pCCellRight()->vSetRegion(iLeftRegion);
    //       pF->pCCellLeft() ->vSetRegion(iRightRegion);
    //     }

    //     // Tag cell to be detached from the mesh based on surface mesh
    //     // orientation
    //     if (qIsInteriorBoundary) {
    //       assert(iRightRegion != iOutsideRegion);
    //       assert(iLeftRegion  != iOutsideRegion);
    //       apCDetach[iCell] = pCInvalidCell;
    //     }
    //     else {
    //       // There's something wrong with the following condition.  It works
    //       // nearly all the time, but for the hex-multi case, it causes
    //       // problems.
    //       if (iFuzzyComp(dPatchDot, dSurfDot) == 0) {
    // 	// When the surface mesh and the patch are oriented the same
    // 	// way, there are no problems.
    // 	if (iRightHandPatch == 1) {
    // 	  apCDetach[iCell] = pF->pCCellLeft();
    // 	}
    // 	else {
    // 	  apCDetach[iCell] = pF->pCCellRight();
    // 	}
    //       }
    //       else {
    // 	// In this case, do a geometric check.
    // 	double adFaceCent[3], adCellCent[3];
    // 	pF->vCentroid(adFaceCent);
    // 	pF->pCCellLeft()->vCentroid(adCellCent);
    // 	double adLeftVec[] = adDIFF3D(adCellCent, adFaceCent);
    // 	NORMALIZE3D(adLeftVec);
    // 	double dDot = dDOT3D(adLeftVec, adNormSurf);
    // 	// If dDot > 0, then for sure this cell lies inside, because the
    // 	// surface mesh is oriented this way on purpose.
    // 	assert(iFuzzyComp(dDot, 0) != 0);
    // 	if (dDot > 0) {
    // 	  apCDetach[iCell] = pF->pCCellRight();
    // 	}
    // 	else {
    // 	  apCDetach[iCell] = pF->pCCellLeft();
    // 	}
    //       }
    //     }
  } // End of loop associating surface faces with volume faces.

#ifndef NDEBUG
  logMessage(4, "Initial markings...\n");
  double dSum = 0;
  double adSum[iOutsideRegion + 1];
  for (int iiR = 0; iiR <= iOutsideRegion; adSum[iiR++] = 0) {
  }
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    double dSize = pC->calcSize();
    int iReg = pC->getRegion();
    if (iReg != iOutsideRegion)
      dSum += dSize;
    adSum[iReg] += dSize;
    logMessage(
	4,
	"Cell %3u. Verts %3u, %3u, %3u, %3u.  Size: %10.4g (%10.4g) Region %d (%10.4g).\n",
	iCell, getVertIndex(pC->getVert(0)), getVertIndex(pC->getVert(1)),
	getVertIndex(pC->getVert(2)), getVertIndex(pC->getVert(3)), dSize, dSum,
	iReg, adSum[iReg]);
  }
#endif

  // For each marked cell in the volume mesh, recursively mark unmarked
  // neighbors
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->getRegion() != iInvalidRegion)
      markCleanNeighbors(pC);
  }

  {
    int aiRegionCounts[iMaxRegionLabel];
    for (int i = 0; i < iMaxRegionLabel; i++)
      aiRegionCounts[i] = 0;
    for (iCell = 0; iCell < getNumCells(); iCell++) {
      int iReg = getCell(iCell)->getRegion();
      aiRegionCounts[iReg]++;
    }
    for (int i = 0; i < iMaxRegionLabel; i++) {
      if (aiRegionCounts[i] != 0)
	logMessage(2, "Region %d has %d tets\n", i, aiRegionCounts[i]);
    }
    logMessage(2, "(Region %d is outside, %d is invalid)\n",
    iOutsideRegion,
	       iInvalidRegion);
  }

#ifndef NDEBUG
  logMessage(4, "Intermediate markings...\n");
  for (int iiR = 0; iiR <= iOutsideRegion; adSum[iiR++] = 0) {
  }
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    double dSize = pC->calcSize();
    int iReg = pC->getRegion();
    if (iReg != iOutsideRegion)
      dSum += dSize;
    adSum[iReg] += dSize;
    logMessage(
	4,
	"Cell %3u. Verts %3u, %3u, %3u, %3u.  Size: %10.4g (%10.4g) Region %d (%10.4g).\n",
	iCell, getVertIndex(pC->getVert(0)), getVertIndex(pC->getVert(1)),
	getVertIndex(pC->getVert(2)), getVertIndex(pC->getVert(3)), dSize, dSum,
	iReg, adSum[iReg]);
  }
#endif

  // Delete the old bdry faces.
  for (GR_index_t iBFace = 0; iBFace < getNumBdryFaces(); iBFace++) {
    BFace *pBF = getBFace(iBFace);
    deleteBFace(pBF);
  }

  // Create the actual boundary that you'll need.
  for (iCell = 0; iCell < pSM->getNumCells(); iCell++) {
    Face *pF = asBdryData[iCell].pFBdry;
    if (!pF->isValid())
      continue;
    Cell *pC = asBdryData[iCell].pCOutside;

    if (pC->isValid()) {
      // This is the regular boundary case.
      if (pC->getType() != Cell::eTriBFace
	  && pC->getType() != Cell::eQuadBFace) {
	// Detach the outside cell from the face on the boundary
	deleteCell(pC);

	// Create a new boundary face and attach it to the face
	BFace *pBF = createBFace(pF);
	assert(pBF->getType() == Cell::eTriBFace);
	pBF->setPatch(asBdryData[iCell].pBP);
      }
      else {
	// Already have a BFace as the cell.
	assert(pF->getFaceLoc() == Face::eBdryFace);
	assert(!pC->isDeleted());
	// pC->vMarkNotDeleted(); // This should be redundant.
	pC->setRegion(iInsideRegion);
      }
    }
    else {
      // Internal boundary.  Need to duplicate the face and add an
      // internal bdry face between the two copies.
      assert(pF->getType() == Face::eTriFace);

      // FIX ME  CFOG, 19 Oct, 2010
      // This is a bit of a hack:  create a fake internal BFace to
      // attach the patch to, then clone this face and stick an internal
      // bdry face in between.  Really should be able to create a bface,
      // then classify it onto a piece of the geometry; or pass the
      // geometry info to the creation routine directly instead of
      // attached to a bdry face.

      // I may have unhacked this.. Dan 2013

      //      BFace *pBFInt = getNewBFace(3, true);
      //      pBFInt->setPatch(pSM->getBdryPatch(iCell));
      //      createBFace(pF, pBFInt);
      //      deleteBFace(pBFInt);
      BFace * pBFInt = createBFace(pF);
      pBFInt->setPatch(pSM->getBdryPatch(iCell));
    }
  }

#ifndef NDEBUG
  logMessage(4, "Final markings...\n");
  dSum = 0;
  for (int iiR = 0; iiR <= iOutsideRegion; adSum[iiR++] = 0) {
  }
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    double dSize = pC->calcSize();
    int iReg = pC->getRegion();
    if (iReg != iOutsideRegion)
      dSum += dSize;
    adSum[iReg] += dSize;
    logMessage(
	4,
	"Cell %3u. Verts %3u, %3u, %3u, %3u.  Size: %10.4g (%10.4g) Region %d (%10.4g).\n",
	iCell, getVertIndex(pC->getVert(0)), getVertIndex(pC->getVert(1)),
	getVertIndex(pC->getVert(2)), getVertIndex(pC->getVert(3)), dSize, dSum,
	iReg, adSum[iReg]);
  }
#endif

  // Mark exterior cells for deletion.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    if (pC->getRegion() == iInvalidRegion || pC->getRegion() == iOutsideRegion)
      deleteCell(pC);
  }

  // Doublecheck to be sure that no kept cells have a deleted vertex.
  // If so, this is an input error.  Further diagnosis requires throwing
  // an exception.
  // FIX ME: 0.2.2.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (!pC->isDeleted()) {
      for (int iVert = 0; iVert < pC->getNumVerts(); iVert++) {
	if (pC->getVert(iVert)->isDeleted()) {
	  vFatalError("Problem with orientation in input file.",
		      "VolMesh::vDeleteExternal()");
	}
      }
    }
  }

  // Mark verts at the corners of the original domain for deletion
  for (int iVert = iNumOrigVerts; iVert < iNumOrigVerts + iNumBBoxVerts;
      iVert++) {
    deleteVert(getVert(iVert));
  }

  // Verify that cells aren't screwed up; this can happen with bad
  // boundary orientation.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    if (!(pC->getFace(0)->isValid() && pC->getFace(1)->isValid()
	&& pC->getFace(2)->isValid() && pC->getFace(3)->isValid()))
      vFatalError("Orientation of an internal boundary is incorrect.",
		  "3D initial meshing from boundary geometry.");
  }

  // This purge messes up attempts to remove vertices that were inserted
  // for surface recovery.
  //   vPurge();
  assert(isValid());

  // #ifndef NDEBUG
  //   for (iCell = 0; iCell < iNumCells(); iCell++)
  //     assert(pCCell(iCell)->iRegion() != 0);
  // #endif
}

