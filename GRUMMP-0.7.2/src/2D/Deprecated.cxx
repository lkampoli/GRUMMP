#include "GR_Geometry.h"
#include "GR_InsertionManager.h"
#include "GR_List.h"
#include "GR_Mesh2D.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"

#ifndef FACET_MERGED
#include "GR_BdryPatch2D.h"
#include "GR_List.h"
#include "GR_Bdry2D.h"

// This type is used to keep track of edges created
// when doing the initial discretization of patches
typedef struct {
  int iPatch;
  int iV0;
  int iV1;
  Face *pF;
  bool qTriedRecover;
} structTrackEdges;

// Someplace to stash the info we need to correctly disconnect the
// exterior part of the mesh and replace it with boundary conditions.
struct FaceInfo {
  Face *pF;
  BFace *pBF;
  Cell *pC;
  FaceInfo() :
      pF(pFInvalidFace), pBF(pBFInvalidBFace), pC(pCInvalidCell)
  {
  }
};

void
Mesh2D::getNeighborVertsAndFaces(Vert * pV, List<Vert *>& LpVNear,
				 List<Face *>& LpFNear)
{
  // Returns a list of neighbor verts and originating faces..
  // Both can be done at the same time..

  bool qWentAround = false;
  bool qHitFace = false;
  bool qCanAdd = true;

  Face *pF, *pOriginalF;
  Cell *pCLeft, *pCRight;
  Face *OldpF;

  pOriginalF = pV->getHintFace();
  pF = pOriginalF;
  pCLeft = pF->getLeftCell();

  //  vMessage(3, "Going left\n");
  //  vMessage(3, "Original Face is %p\n", pF);

  // Go left...
  while ((!qWentAround) && (!qHitFace)) {
    //		vMessage(4, "Now looking at face %d\n", pF);
    if (qCanAdd) {
      // Add the face to the list..
      LpFNear.vAppendItem(pF);
      //      vMessage(4, "Adding %d to face list\n", pF);
      // Add the other vertex to list
      Vert *pVOtherV = pF->getVert(0) == pV ? pF->getVert(1) : pF->getVert(0);
      LpVNear.vAppendItem(pVOtherV);
    }
    else {
      // Not allowed to add because it was an interior boundary
    }				// if qCanAdd

    // Get the left cell
    // Check if it's a BFace we just went across..
    if (pCLeft->getType() == Cell::eBdryEdge) {
      qHitFace = true;
      //      vMessage(4, "Face was a boundary face\n");
      // No need to get another face..
    }
    else {
      // Must find the other face connected to the same vertex
      // If you hit an internal boundary, go through it..
      if (pCLeft->getType() == Cell::eIntBdryEdge) {
	// vMessage(4, "Face is an internal boundary face\n");
	// Only want to have one of the two faces around an internal
	// boundary in the list..
	OldpF = pF;
	pF = dynamic_cast<IntBdryEdge *>(pCLeft)->getOtherFace(OldpF);
	// vMessage(4, "Changing face %d to face %d\n", OldpF, pF);
	// Only want one of the two faces in the list,
	// and only want the vertex once..
	// so this one won't be added
	qCanAdd = false;
      }
      else {
	//		vMessage(4, "Regular face, getting new one...\n");
	for (int iF = pCLeft->getNumFaces() - 1; iF >= 0; iF--) {
	  Face *pFNew = pCLeft->getFace(iF);
	  if (pFNew != pF && pFNew->hasVert(pV)) {
	    pF = pFNew;
	    break;
	  }
	}
	//				vMessage(4, "New face is %d\n", pF);
	qCanAdd = true;
      }				// if it's an eIntBdryEdge

      Cell *OldpC = pCLeft;
      pCLeft = pF->getOppositeCell(OldpC);
      // check if it's the same..
      if (pF == pOriginalF) {
	qWentAround = true;
	//				vMessage(4, "Same as original face -- done\n");
	// If the first face was an interior boundary face, then
	// the two faces around the interior boundary would be first
	// and last, and we only want one of them in the list. So remove
	// the last one inserted

	if (OldpC->getType() == Cell::eIntBdryEdge) {
	  //					vMessage(4, "Last cell was an interior boundary, removing last face inserted\n");
	  LpFNear.vDelete(OldpF);
	  Vert *pVDel =
	      OldpF->getVert(0) == pV ? OldpF->getVert(1) : OldpF->getVert(0);
	  LpVNear.vDelete(pVDel);
	  //					vMessage(4, "Removed last vertex too\n");
	}

      }				// if it's the original face
    }				// if it's an eBdryEdge
  }				// while

  //  vMessage(4, "Now out of first loop\n");

  qHitFace = false;
  qCanAdd = true;

  pF = pOriginalF;
  pCRight = pF->getRightCell();

  //  vMessage(4, "Starting right with OriginalpF = %p\n", pF);

  if (pCRight->getType() == Cell::eBdryEdge) {
    qHitFace = true;
    //    vMessage(4, "Hit Bdry right away, done\n");
  }
  else {
    if (pCRight->getType() == Cell::eIntBdryEdge) {
      // vMessage(4, "Int Bdry right away....\n");
      OldpF = pF;
      pF = dynamic_cast<IntBdryEdge *>(pCRight)->getOtherFace(OldpF);
      // vMessage(4, "Changing %d with %d\n", OldpF, pF);
      qCanAdd = false;
    }
    else {
      //    	vMessage(4, "Getting adjacent face to the right\n");
      for (int iF = pCRight->getNumFaces() - 1; iF >= 0; iF--) {
	Face *pFNew = pCRight->getFace(iF);
	if (pFNew != pOriginalF && pFNew->hasVert(pV)) {
	  pF = pFNew;
	  break;
	}
      }
      qCanAdd = true;
    }
    pCRight = pF->getOppositeCell(pCRight);
  }

  //  vMessage(4, "Starting right with pF = %d\n", pF);

  // Go right...
  while ((!qWentAround) && (!qHitFace)) {
    if (qCanAdd) {
      // Add the face to the list..
      LpFNear.vAppendItem(pF);
      //      vMessage(4, "Adding face %d to list\n");
      // Add the other vertex to list
      Vert *pVOtherV = pF->getVert(0) == pV ? pF->getVert(1) : pF->getVert(0);
      LpVNear.vAppendItem(pVOtherV);
    }
    else {
      // Not allowed to add..
    }
    // Get the left cell
    // Check if it's a BFace we just went across..
    if (pCRight->getType() == Cell::eBdryEdge) {
      //      vMessage(4, "Face is a bdry\n");
      qHitFace = true;
      // No need to get another face..
    }
    else {
      // Must find the other face connected to the same vertex
      // Make sure to go through internal boundaries..
      if (pCRight->getType() == Cell::eIntBdryEdge) {
	OldpF = pF;
	pF = dynamic_cast<IntBdryEdge *>(pCRight)->getOtherFace(OldpF);

	//	vMessage(4, "Face is an internal bdry\n");
	//	vMessage(4, "Swapping %p for %p\n", OldpF, pF);
	qCanAdd = false;
      }
      else {
	for (int iF = pCRight->getNumFaces() - 1; iF >= 0; iF--) {
	  Face *pFNew = pCRight->getFace(iF);
	  if (pFNew != pF && pFNew->hasVert(pV)) {
	    pF = pFNew;
	    break;
	  }
	}
	//	vMessage(4, "Getting new face %d (%p)\n", iFaceIndex(pF), pF);
	qCanAdd = true;
      }

      Cell *OldpC = pCRight;
      pCRight = pF->getOppositeCell(OldpC);
      // check if it's the same..
      if (pF == pOriginalF) {
	//		vMessage(4, "Same as original face..\n");
	qWentAround = true;
      }
    }
  }
  //  vMessage(4, "Exiting from GetNeighborVertsAndFaces()\n");
}

//////////////////////////////////////////////////////////////
//
//@ Triangulation of a generalized boundary geometry
//
//////////////////////////////////////////////////////////////
Mesh2D::Mesh2D(const Bdry2D & B2DInfo, const int iQualMeas) :
    Mesh(), m_ECEdgeF(), m_ECEdgeBF(), m_ECTri(), m_ECQuad()
{
  m_qual = new Quality(this, iQualMeas);

  // List of edges to recover
  List<structTrackEdges*> LsTrackEdges;
  structTrackEdges *pEdge;
  int EdgeiV0 = -1, EdgeiV1 = -1, EdgeiPatch;

  // Concurrent list to keep track of what patch a vert is associated
  // with.
  // LsVertsPatch[0] refers to vert 0, etc...
  List<int> LsVertsPatch;

  m_LipschitzAlpha = 1.0;
  m_resolution = 1.0;

  // Set up an outer triangulation of a big square
  double dXMin, dXMax, dYMin, dYMax;
  dXMin = dYMin = 1.e100;
  dXMax = dYMax = -1.e100;

  // Find max/min of segments
  int iPatch, iNBC = B2DInfo.iNumPatches();
  for (iPatch = 0; iPatch < iNBC; iPatch++) {
    double dSegXMin = B2DInfo[iPatch]->dXMin();
    double dSegXMax = B2DInfo[iPatch]->dXMax();
    double dSegYMin = B2DInfo[iPatch]->dYMin();
    double dSegYMax = B2DInfo[iPatch]->dYMax();

    dXMin = std::min(dXMin, dSegXMin);
    dXMax = std::max(dXMax, dSegXMax);
    dYMin = std::min(dYMin, dSegYMin);
    dYMax = std::max(dYMax, dSegYMax);
  }

  double dXLo = dXMin - (dXMax - dXMin);
  double dXHi = dXMax + (dXMax - dXMin);
  double dYLo = dYMin - (dYMax - dYMin);
  double dYHi = dYMax + (dYMax - dYMin);

  logMessage(3, "dXLo: %f dXHi: %f dYLo: %f dYHi: %f\n", dXLo, dXHi, dYLo,
	     dYHi);

  GR_index_t iNumOrigVerts = iNumVerts();

  // add dummy entries..
  for (int iDum = 0; iDum < 4; iDum++) {
    LsVertsPatch.vAppendItem(-100);
  }

  createMeshInBox(*this, dXLo, dYLo, dXHi, dYHi);

  // Do Delaunay insertion of all the existing verts in the mesh
  setSwapType(eDelaunay);
  //	vMessage(2, "Number of verts now: %d\n", iNumVerts());
  //	for (int kk = 0; kk < iNumVerts(); kk++) {
  //		vMessage(2, "  %d -- (%f, %f)\n", kk, pVVert(kk)->dX(), pVVert(kk)->dY());
  //	}

  logMessage(3,
	     "Going through patches to determine original discretization... ");

  int iNBdryPatches = B2DInfo.iNumPatches();
  int iOriginalV = -1;
  BdryPatch2D *pPatch;
  List<int> LiInsertedVerts;
  List<int> LiMeshToBdryVerts;
  LiInsertedVerts.vClear();
  LiInsertedVerts.vDisallowSorting();
  LiMeshToBdryVerts.vClear();
  LiMeshToBdryVerts.vDisallowSorting();
  GRUMMP::DelaunaySwapDecider2D Del2D;
  GRUMMP::SwapManager2D SM2D(&Del2D, this);
  GRUMMP::SwappingInserter2D SI2D(this, &SM2D);
  SI2D.setForcedInsertion(true);
  for (iPatch = 0; iPatch < iNBdryPatches; iPatch++) {
    // Get the patch
    pPatch = B2DInfo[iPatch];
    // Call the original discretization info
    List<double> La2dOtherV;
    bool qInsertPointers;
    bool qLoop;
    double adCoords[2];

    La2dOtherV.vClear();

    La2dOtherV = pPatch->LdOriginalDiscretization(&qInsertPointers, &qLoop);

    // Insert the necessary verts, and build edge tracking list
    // at same time...

    bool qOriginalV = true;
    int iVert;
    int iMeshVert;
    Vert *pV;

    if (qInsertPointers) {
      iVert = pPatch->iBeginPointIndex();
      iMeshVert = iNumVerts();
      //      pV = pVVert(iVert);

      // We know iVert is associated to iPatch
      LsVertsPatch.vAppendItem(iPatch);
      pPatch->vGetBeginPoint(adCoords);

      //      adCoords[0] = pV->dX();
      //      adCoords[1] = pV->dY();

      if (qOriginalV) {
	iOriginalV = iMeshVert;
	qOriginalV = false;
      }

      // Insert pV0 for this patch
      // pV->vSetType(Vert::eBdryApex);

      //			vMessage(2, "iVert is: %d\n", iVert);
      if (!LiInsertedVerts.qInList(iVert)) {
	logMessage(4, "Inserting vert (%f, %f)\n", adCoords[0], adCoords[1]);
	// Must create the vertex
	pV = createVert(adCoords);
	pV->setType(Vert::eBdryApex);
	if (B2DInfo.qVertSmallAngle(iVert)) {
	  pV->markAsSmallAngleVert(true);
	}

	Cell *pCGuess = getCell(0);
	for (int ii = 0; pCGuess->isDeleted(); ii++) {
	  pCGuess = getCell(ii);
	}
#ifndef NDEBUG
	Vert* newVert = SI2D.insertPoint(adCoords, pCGuess, pV);
	assert(newVert->isValid());
#else
	(void) SI2D.insertPoint(adCoords, pCGuess, pV);
#endif

	LiInsertedVerts.vAppendItem(iVert);
	LiMeshToBdryVerts.vAppendItem(iMeshVert);

	//	vMessage(2, "Number of verts now: %d\n", iNumVerts());
	//	for (int kk = 0; kk < iNumVerts(); kk++) {
	//	  vMessage(2, "  %d -- (%f, %f)\n", kk, pVVert(kk)->dX(), pVVert(kk)->dY());
	//	}
	EdgeiV0 = iMeshVert;

      }
      else {
	// iV0 was already inserted. Find out which one it is in the mesh
	int iTarget = LiInsertedVerts.iSeek(iVert);
	EdgeiV0 = LiMeshToBdryVerts[iTarget];
      }

    } // if pointers are supposed to be inserted

    // now insert the middle "extra" verts specificied in La2dOtherV

    if (La2dOtherV.iLength() > 0) {
      for (int i = 0; i < La2dOtherV.iLength() / 2; i++) {

	// insert point
	adCoords[0] = La2dOtherV[2 * i];
	adCoords[1] = La2dOtherV[(2 * i) + 1];

	iMeshVert = iNumVerts();

	// Create new vert
	pV = createVert(adCoords);
	// we know the patch for this vert...
	LsVertsPatch.vAppendItem(iPatch);

	// Set vert info
	pV->setType(Vert::eBdryApex);

	if (!qOriginalV) {
	  EdgeiV1 = iNumVerts() - 1;
	  EdgeiPatch = iPatch;

	  pEdge = new structTrackEdges;
	  pEdge->iV0 = EdgeiV0;
	  pEdge->iV1 = EdgeiV1;
	  pEdge->iPatch = EdgeiPatch;
	  pEdge->pF = pFInvalidFace;
	  pEdge->qTriedRecover = false;
	  LsTrackEdges.vAppendItem(pEdge);

	  // prepare next one..
	  EdgeiV0 = EdgeiV1;

	}
	else {
	  qOriginalV = false;
	  EdgeiV0 = iNumVerts() - 1;
	  iOriginalV = EdgeiV0;
	}

	logMessage(4, "Inserting vert (%f, %f)\n", adCoords[0], adCoords[1]);

	Cell *pCGuess = getCell(0);
	for (int ii = 0; pCGuess->isDeleted(); ii++) {
	  pCGuess = getCell(ii);
	}
#ifndef NDEBUG
	Vert* newVert = SI2D.insertPoint(adCoords, pCGuess, pV);
	assert(newVert->isValid());
#else
	(void) SI2D.insertPoint (adCoords, pCGuess, pV);
#endif
	// This vertex is not in the Bdry2D list of points..
	int iDum = -1;
	LiInsertedVerts.vAppendItem(iDum);
	LiMeshToBdryVerts.vAppendItem(iMeshVert);

      } // loop over "middle" verts
    } // if there were any middle verts

    if (qInsertPointers) {
      // insert the "end" vert now
      iVert = pPatch->iEndPointIndex();
      iMeshVert = iNumVerts();
      //      pV = pVVert(iVert);

      LsVertsPatch.vAppendItem(iPatch);

      pPatch->vGetEndPoint(adCoords);
      //      adCoords[0] = pV->dX();
      //      adCoords[1] = pV->dY();

      //      vMessage(2, "Adding Patch: %d, iV0: %d, iV1: %d to list\n",
      //				pEdge->iPatch, pEdge->iV0, pEdge->iV1);
      //			vMessage(2, "(%f, %f) -- (%f, %f)\n",
      //				pVVert(pEdge->iV0)->dX(), pVVert(pEdge->iV0)->dY(),
      //				pVVert(pEdge->iV1)->dX(), pVVert(pEdge->iV1)->dY());

      // Insert pV1 for this patch
      //      pV->vSetType(Vert::eBdryApex);

      if (!LiInsertedVerts.qInList(iVert)) {
	logMessage(4, "Inserting vert (%f, %f)\n", adCoords[0], adCoords[1]);

	// Must create the vertex
	pV = createVert(adCoords);
	pV->setType(Vert::eBdryApex);
	if (B2DInfo.qVertSmallAngle(iVert)) {
	  pV->markAsSmallAngleVert(true);
	}

	Cell *pCGuess = getCell(0);
	for (int ii = 0; pCGuess->isDeleted(); ii++) {
	  pCGuess = getCell(ii);
	}
#ifndef NDEBUG
	Vert* newVert = SI2D.insertPoint(adCoords, pCGuess, pV);
	assert(newVert->isValid());
#else
	(void) SI2D.insertPoint(adCoords, pCGuess, pV);
#endif

	LiInsertedVerts.vAppendItem(iVert);
	LiMeshToBdryVerts.vAppendItem(iMeshVert);

	//
	EdgeiV1 = iMeshVert;

      }
      else {
	// This vertex was already inserted. Use the
	// mesh index stored in LiMeshToBdryVerts
	int iTarget = LiInsertedVerts.iSeek(iVert);
	EdgeiV1 = LiMeshToBdryVerts[iTarget];

      }

      EdgeiPatch = iPatch;

      pEdge = new structTrackEdges;
      pEdge->iV0 = EdgeiV0;
      pEdge->iV1 = EdgeiV1;
      pEdge->iPatch = EdgeiPatch;
      pEdge->pF = pFInvalidFace;
      pEdge->qTriedRecover = false;
      LsTrackEdges.vAppendItem(pEdge);

    } // if end pointer was to be inserted...

    // Now check if we must loop
    if (qLoop) {
      EdgeiV0 = EdgeiV1;
      EdgeiV1 = iOriginalV;
      EdgeiPatch = iPatch;

      pEdge = new structTrackEdges;
      pEdge->iV0 = EdgeiV0;
      pEdge->iV1 = EdgeiV1;
      pEdge->iPatch = EdgeiPatch;

      pEdge->pF = pFInvalidFace;
      pEdge->qTriedRecover = false;
      LsTrackEdges.vAppendItem(pEdge);
    }
  } // loop over all patches

  SM2D.swapAllFaces();
  logMessage(3, "done.\nNow trying to recover all of the edges... ");
  logMessage(3, "List of edges to be recovered contains %d items\n",
	     LsTrackEdges.iLength());
  for (int j = 0; j < LsTrackEdges.iLength(); j++) {
    pEdge = LsTrackEdges[j];
    logMessage(4, "  %d) iV0: %d, iV1: %d, iPatch: %d\n", j, pEdge->iV0,
	       pEdge->iV1, pEdge->iPatch);
  }
  setAllHintFaces();
  setVertFaceNeighbors();

  assert(isValid());
  // Another step must be performed here.
  // We will make sure that all the boundary edges can be recovered
  // by using something similar to stitching.
  //
  // If an edge that is in the list of edges does NOT exist in the
  // triangulation then insert a vertex at its midpoint (making sure the
  // list gets updated properly).
  //
  // This fixes problems that arise when two boundary edges (wrongly) cross
  // each other.
  //
  // Now, everytime an edge that IS present in the mesh is found,
  // check for the following:
  //  - check for encroachment (2 opposite verts)
  //  - If it is encroached on, check distance to that bdry patch
  //  - If distance is positive, split.
  //
  // When this is over, the encroachment check step can be skipped
  // altogether..!
  bool qAllEdgesPresent = false;
  int iEdge = 0;
  int iNumEdgesPresent = 0;
  while (!qAllEdgesPresent) {
    // Look at iEdge
    // if it is OK, iNumEdgesPresent++
    // iEdge++ <-- Make sure it does not go too far

    // Figure out if the edge is present in mesh
    bool qEdgePresent = false;
    // Get the edge info
    pEdge = LsTrackEdges[iEdge];
    int iV0, iV1;
    Vert *pV0, *pV1;
    // Indices
    iV0 = pEdge->iV0;
    iV1 = pEdge->iV1;
    iPatch = pEdge->iPatch;
    // Get pointers
    pV0 = getVert(iV0);
    pV1 = getVert(iV1);
    logMessage(4, "Looking at edge %d -- (%f, %f)->(%f, %f) -- Patch #%d\n",
	       iEdge, pV0->x(), pV0->y(), pV1->x(), pV1->y(), iPatch);

    // Can save some time if we already have a face pointer
    if (pEdge->pF != pFInvalidFace) {
      // Edge was at some point present. Make sure it
      // still is.
      if (!pEdge->pF->isDeleted()) {
	Vert *pFV0 = pEdge->pF->getVert(0);
	Vert *pFV1 = pEdge->pF->getVert(1);

	// Make sure vertices are the same
	qEdgePresent = (((pFV0 == pV0) && (pFV1 == pV1))
	    || ((pFV0 == pV1) && (pFV1 == pV0)));
      }
      else {
	// Edge was deleted..
	qEdgePresent = false;
      }

    }
    else {
      // No pointer to the face -- Have to determine if
      // the edge is present in the mesh

      // Neighbor lists
      List<Vert*> LpVNear;
      List<Face*> LpFNear;
      // Get the neighbors
      getNeighborVertsAndFaces(pV0, LpVNear, LpFNear);
      // Go through vertex list -- is pV1 in there?
      int iV = 0;
      while ((iV < LpVNear.iLength()) && (qEdgePresent == false)) {
	if (LpVNear[iV] == pV1)
	  qEdgePresent = true;
	iV++;
      }

      if (qEdgePresent) {
	// Set face pointer -- will save time later
	pEdge->pF = LpFNear[iV - 1];
      }
    } // whether or not a pointer to the face existed already

    // Check for encroachment if edge exists
    bool qEncroach = false;
    if (qEdgePresent) {
      // Make sure it is not encroached on if it's a curved boundary edge
      //				if (B2DInfo[iPatch]->eGeom2DType != BdryPatch2D::ePolyLine) {
      // Check for encroachment. Find 2 opposite verts.
      Cell *pCL, *pCR;
      Vert *pVL, *pVR;

      double adC0[2], adC1[2];

      adC0[0] = pV0->x();
      adC0[1] = pV0->y();
      adC1[0] = pV1->x();
      adC1[1] = pV1->y();

      // Get the midpoint of the edge...
      double adMid[2];
      B2DInfo[iPatch]->vRatioBetweenTwoPoints(0.5, adC0, adC1, adMid);

      // Is the left cell *really* to the left
      // of the patch..?

      pCL = pEdge->pF->getLeftCell();
      pCR = pEdge->pF->getRightCell();
      pVL = pCL->getOppositeVert(pEdge->pF);
      pVR = pCR->getOppositeVert(pEdge->pF);

      // Get orientations...
      int iOldL, iOldR, iNewL, iNewR;
      iOldL = checkOrient2D(pV0, pV1, pVL);
      iOldR = checkOrient2D(pV0, pV1, pVR);
      iNewL = checkOrient2D(adC0, adMid, pVL->getCoords());
      iNewR = checkOrient2D(adC0, adMid, pVR->getCoords());

      // Check difference in orientations
      if ((iOldL != iNewL) || (iOldR != iNewR)) {
	// Hmm.. Something wrong, I think....
	//				vMessage(3, "NEW METHOD: Wrong orientations!\n");
	qEncroach = true;
      }
      /*
       iPL = LsVertsPatch[iVertIndex(pVL)];
       iPR = LsVertsPatch[iVertIndex(pVR)];

       B2DInfo[iPatch]->vDistWhichSide(pVL->adCoords(), &dDL, &iSideL);
       B2DInfo[iPatch]->vDistWhichSide(pVR->adCoords(), &dDR, &iSideR);

       vMessage(3, "Left vert: (%f, %f) -- Dist: %f -- Side: %d\n",
       pVL->dX(), pVL->dY(), dDL, iSideL);
       vMessage(3, "Right vert: (%f, %f) -- Dist: %f -- Side: %d\n",
       pVR->dX(), pVR->dY(), dDR, iSideR);


       // Make sure other vert is not part of the same patch..
       if ( (iFuzzyComp(dDR, 0.0) > 0) && (iFuzzyComp(dDL, 0.0) > 0) ) {
       // If the side is the same... there is something wrong!
       if (iSideL == iSideR) {
       qEncroach = true;
       vMessage(3, "OLD METHOD: Wrong orientations!\n");
       }
       }
       */
    }

    // if these two are OK, then increment, otherwise, split!
    if ((qEdgePresent) && (!qEncroach)) {
      iNumEdgesPresent++;
    }
    else {
      // Reset counter
      iNumEdgesPresent = 0;
      // Make sure the face pointer is invalid
      pEdge->pF = pFInvalidFace;
      // Get the patch
      pPatch = B2DInfo[iPatch];

      if ((pPatch->eGeom2DType == BdryPatch2D::ePolyLine)
	  || ((!pEdge->qTriedRecover) && (!qEncroach))) {
	// Try to recover the old way...
	logMessage(4, "Trying to recover that edge...\n");
	recoverEdge(pEdge->iV0, pEdge->iV1);
	pEdge->qTriedRecover = true;
      }
      else {
	// Already tried to recover the "old" way, so do something else...
	// If it's not a line,then insert a point...
	double adLoc1[2], adLoc2[2], adNewPoint[2];
	adLoc1[0] = pV0->x();
	adLoc1[1] = pV0->y();
	adLoc2[0] = pV1->x();
	adLoc2[1] = pV1->y();
	// Get the point
	pPatch->vRatioBetweenTwoPoints(0.5, adLoc1, adLoc2, adNewPoint);
	// insert it!
	int iVNew = iNumVerts();
	Vert *pV = createVert(adNewPoint);
	// Associate patch with it
	LsVertsPatch.vAppendItem(iPatch);

	// Set type;
	pV->setType(Vert::eBdryApex);
	//
	Cell *pCGuess = getCell(0);
	for (int ii = 0; pCGuess->isDeleted(); ii++) {
	  pCGuess = getCell(ii);
	}
#ifndef NDEBUG
	Vert* newVert = SI2D.insertPoint(adNewPoint, pCGuess, pV);
	assert(newVert->isValid());
#else
	(void) SI2D.insertPoint(adNewPoint, pCGuess, pV);
#endif
	// Delete the old edge from the tracking list
	LsTrackEdges.vDelete(pEdge);

	// add the two new ones
	EdgeiV0 = iV0;
	EdgeiV1 = iVNew;
	EdgeiPatch = iPatch;

	pEdge = new structTrackEdges;
	pEdge->iV0 = EdgeiV0;
	pEdge->iV1 = EdgeiV1;
	pEdge->iPatch = EdgeiPatch;
	pEdge->pF = pFInvalidFace;
	pEdge->qTriedRecover = false;
	LsTrackEdges.vAppendItem(pEdge);

	EdgeiV0 = iVNew;
	EdgeiV1 = iV1;

	pEdge = new structTrackEdges;
	pEdge->iV0 = EdgeiV0;
	pEdge->iV1 = EdgeiV1;
	pEdge->iPatch = EdgeiPatch;
	pEdge->pF = pFInvalidFace;
	pEdge->qTriedRecover = false;
	LsTrackEdges.vAppendItem(pEdge);
      }
    }

    // Do the iEdge increment properly
    if (iEdge < LsTrackEdges.iLength() - 1) {
      // Can increment
      iEdge++;
    }
    else {
      // Loop back
      iEdge = 0;
    }

    // Are we there yet?
    if (iNumEdgesPresent == LsTrackEdges.iLength())
      qAllEdgesPresent = true;
  }
  //   vMessage(4, "Now marking all verts at the ends of curved patches as apexes\n");
  //   for (iPatch = 0; iPatch < iNBdryPatches; iPatch++) {
  //     // Get the patch
  //     pPatch = B2DInfo[iPatch];
  //     // Do nothing for polylines and patches of unknown type (the latter
  //     // are impossible anyway).
  //     if (pPatch->eGeom2DType == BdryPatch2D::eUnknown ||
  //	pPatch->eGeom2DType == BdryPatch2D::ePolyLine) continue;
  //     // For other patches, mark their endpoints as apexes to prevent
  //     // inappropriate surface smoothing.  The "if" is there because
  //     // circles have -1 as their begin and end indices.
  //     int iBeginVert = pPatch->iBeginPointIndex();
  //     if (iBeginVert >= 0) pVVert(iBeginVert)->vSetType(Vert::eBdryApex);
  //     int iEndVert   = pPatch->iEndPointIndex();
  //     if (iEndVert >= 0) pVVert(iEndVert)->vSetType(Vert::eBdryApex);
  //   }
  // Mark all cells as having an invalid region, since we don't yet know
  // where any of them are.
  GR_index_t iCell;
  for (iCell = 0; iCell < getNumCells(); iCell++)
    getCell(iCell)->setRegion(iInvalidRegion);

  struct FaceInfo *aFI = new struct FaceInfo[LsTrackEdges.iLength()];

  // For each boundary patch:
  //   1. Mark cells that are in the interior of the domain to be meshed
  //   with their region.
  //   2. Mark cells that are in the exterior for deletion.
  int iNumOneBdryPatches = 0;
  logMessage(4, "Assign proper boundary patches data to neighbor cells\n");
  for (iEdge = 0; iEdge < LsTrackEdges.iLength(); iEdge++) {

    logMessage(4, "Edge %d\n", iEdge);
    pEdge = LsTrackEdges[iEdge];

    int iV0 = pEdge->iV0;
    iPatch = pEdge->iPatch;
    Vert *pV0 = getVert(iV0);
    BFace *pBF;

    // Edge was already found before
    Face *pF = pEdge->pF;
    //     vMessage(4, "Face pointer is %p\n", pF);

    assert(pF->isValid());
#ifndef NDEBUG
    int iV1 = pEdge->iV1;
    Vert *pV1 = getVert(iV1);
    assert(pF->getVert(0) == pV0 || pF->getVert(0) == pV1);
    assert(pF->getVert(1) == pV0 || pF->getVert(1) == pV1);
#endif

    logMessage(4, "   pVVert(0): (%f, %f) -- pVVert(1) (%f, %f)\n",
	       pF->getVert(0)->x(), pF->getVert(0)->y(), pF->getVert(1)->x(),
	       pF->getVert(1)->y());

    // Determine whether the edge and the bdry patch are pointing the
    // same direction.
    bool qSameOrient = (pF->getVert(0) == pV0);

    logMessage(4, "Marking regions..\n");

    BdryPatch2D *pBPatch = B2DInfo[iPatch];
    // 1. Mark regions.
    if (pBPatch->iLeftRegion() != iInvalidRegion) {
      if (qSameOrient) {
	pF->getLeftCell()->setRegion(pBPatch->iLeftRegion());
	logMessage(4, "Marking left cell with left region number\n");
      }
      else {
	pF->getRightCell()->setRegion(pBPatch->iLeftRegion());
	logMessage(4, "Marking right cell with left region number\n");
      }
    }
    if (pBPatch->iRightRegion() != iInvalidRegion) {
      if (qSameOrient) {
	pF->getRightCell()->setRegion(pBPatch->iRightRegion());
	logMessage(4, "Marking right cell with right region number\n");
      }
      else {
	pF->getLeftCell()->setRegion(pBPatch->iRightRegion());
	logMessage(4, "Marking left cell with right region number\n");
      }
    }

    // 2. Mark boundaries.
    // Is one side of the edge a boundary?
    if (pBPatch->iBCWhichSide() != BC_SIDE_NONE) {
      // Make sure this is only done when it's an external boundary

      logMessage(4, "   This patch is an external boundary\n");
      Cell *pC;
      if (pBPatch->iBCWhichSide() == BC_SIDE_LEFT) {
	// BC is left side of the patch
	if (qSameOrient)
	  pC = pF->getLeftCell();
	else
	  pC = pF->getRightCell();
      }
      else {
	// BC is right side of the patch...
	if (qSameOrient)
	  pC = pF->getRightCell();
	else
	  pC = pF->getLeftCell();
      }

      pBF = getNewBFace();

      // Make sure the face is set properly.
      pBF->assignFaces(pF);
      pF->setFaceLoc(Face::eBdryFace);

      // This is used later on to delete the unused exterior cells
      aFI[iNumOneBdryPatches].pC = pC;
      aFI[iNumOneBdryPatches].pF = pF;
      aFI[iNumOneBdryPatches].pBF = pBF;
      iNumOneBdryPatches++;
    }
    else {
      // This is an internal boundary

      logMessage(4, "   This patch is an internal boundary\n");

      // Must create...
      //  1. An internal boundary face..
      //  2. Another face. The IBF is located "between" the two faces..

      logMessage(4, "   Creating the extra face\n");

      Face *pF2 = getNewFace();

      logMessage(4, "   Extra face is: %u (%p)\n", getFaceIndex(pF2), pF2);

      logMessage(4, "   Creating the new IntBFace..\n");

      pBF = getNewIntBdryEdge();

      // The new face will be the "right" face of the IBF.
      // This means that the IBF is located on the left side
      // of the new face..

      logMessage(4, "   Assigning the cells and verts to the faces\n");

      Cell *pCOther;
      if (qSameOrient) {
	pCOther = pF->getRightCell();
	pF2->assign(pBF, pCOther, pF->getVert(0), pF->getVert(1));
      }
      else {
	pCOther = pF->getLeftCell();
	pF2->assign(pCOther, pBF, pF->getVert(0), pF->getVert(1));
      }

      pF->replaceCell(pCOther, pBF);
      pCOther->replaceFace(pF, pF2);

      logMessage(4, "   Setting faces types\n");

      // Make sure they know they are internal boundaries
      pF2->setFaceLoc(Face::eBdryTwoSide);
      pF->setFaceLoc(Face::eBdryTwoSide);

      logMessage(4, "   Assigning faces to the IBF\n");

      // Make sure IBF knows about the two faces...
      //
      // Is the order really important here..?
      //
      dynamic_cast<IntBdryEdge*>(pBF)->assignFaces(pF, pF2);
    }

    // Set the patch index
    dynamic_cast<BdryEdgeBase *>(pBF)->setPatch(B2DInfo[iPatch]);

  }				// for loop over all patches

  logMessage(3, "done.\nNow deleting unnecessary cells and faces... ");

  // For every edge separating interior and exterior, replace the
  // exterior cell with a boundary face (which has a pointer to the
  // correct boundary patch so that geometric info is not lost).

  logMessage(4, "Now deleting the exterior cells of Bdry Faces...\n");
  for (iPatch = 0; iPatch < iNumOneBdryPatches; iPatch++) {
    // Only do this stuff if this is -really- a boundary face.
    logMessage(4, "   Checking exterior patch %d\n", iPatch);
    if (aFI[iPatch].pF->isValid()) {
      Cell *pC = aFI[iPatch].pC;
      Face *pF = aFI[iPatch].pF;
      BFace *pBF = aFI[iPatch].pBF;

      // Disconnect the old cell.
      pC->removeFace(pF);

      // Replace the old cell with the new bdry face in the face
      // connectivity.
      pF->replaceCell(pC, pBF);

    }
  }
  for (iPatch = 0; iPatch < iNumOneBdryPatches; iPatch++) {
    // Only do this stuff if this is -really- a boundary face.
    logMessage(4, "   Checking exterior patch %d\n", iPatch);
    if (aFI[iPatch].pF->isValid()) {
      Cell *pC = aFI[iPatch].pC;
      deleteCell(pC);
    }
  }
  delete[] aFI;

  logMessage(4, "Marking regions\n");

  // March out from marked cells, filling all unmarked neighbors with
  // correct region tag (interior) or marking for deletion (exterior) in
  // a second pass.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    // Can't march out from cells that don't know where they are, and
    // don't need to march out from cells in the exterior.
    if (pC->getRegion() == iInvalidRegion)
      continue;
    markCleanNeighbors(pC);
  }

  logMessage(4, "Marking exterior cells for deletion\n");

  // Mark exterior cells for deletion.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (!pC->isDeleted() && pC->getRegion() == iInvalidRegion) {
      //       pC->vMarkDeleted();
      deleteCell(pC);
    }
  }

  logMessage(4, "Marking exterior verts for deletion\n");

  // Mark farfield boundary faces for deletion.  These are always the
  // first four BFaces.
  GR_index_t iBFace;
  for (iBFace = 0; iBFace < 4; iBFace++) {
    deleteBFace(getBFace(iBFace));
    //     pBFBFace (iBFace)->vMarkDeleted ();
  }

  // Mark all faces bounded by two deleted cells for deletion. This must
  // come -after- bdry face deletion so that the farfield faces are
  // deleted properly.
  GR_index_t iFace;
  for (iFace = 0; iFace < getNumFaces(); iFace++) {
    Face *pF = getFace(iFace);
    if (pF->isDeleted())
      continue;
    if ((!pF->getLeftCell()->isValid() || pF->getLeftCell()->isDeleted())
	&& (!pF->getRightCell()->isValid() || pF->getRightCell()->isDeleted())) {
      //       pF->vMarkDeleted();
      deleteFace(pF);
    }
  }

  // Verify that markings are consistent in that all interior regions
  // are bounded away from each other by region boundaries.
  for (iFace = 0; iFace < getNumFaces(); iFace++) {
    Face *pF = getFace(iFace);
    if (pF->isDeleted())
      continue;			// Don't bother about deleted edges.

    if (pF->getFaceLoc() == Face::eBdryFace)
      continue;			// Ditto for boundaries.

    int iRegL = pF->getLeftCell()->getRegion();
    int iRegR = pF->getRightCell()->getRegion();
    if ((iRegL != iRegR) && (pF->getFaceLoc() != Face::eBdryTwoSide))
      vFatalError("Orientation of an internal boundary is incorrect.",
		  "2D initial meshing from boundary geometry.");
  }

  // Mark four farfield vertices for deletion.
  GR_index_t iVert;
  for (iVert = iNumOrigVerts; iVert < iNumOrigVerts + 4; iVert++) {
    getVert(iVert)->markAsDeleted();
    //     deleteVert(pVVert(iVert));
  }

  // Doublecheck to be sure that no kept cells have a deleted vertex.
  // If so, this is an input error.  Further diagnosis requires throwing
  // an exception.
  // FIX ME: 0.2.2.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (!pC->isDeleted()) {
      for (int iiV = 0; iiV < pC->getNumVerts(); iiV++) {
	if (pC->getVert(iiV)->isDeleted()) {
	  vFatalError("Problem with polygon orientation in input file.",
		      "Mesh2D::Mesh2D(Bdry2D)");
	}
      }
    }
  }

  // Verify that cells aren't screwed up; this can happen with bad
  // boundary orientation.
  for (iCell = 0; iCell < getNumCells(); iCell++) {
    Cell *pC = getCell(iCell);
    if (pC->isDeleted())
      continue;
    if (!(pC->getFace(0)->isValid() && pC->getFace(1)->isValid()
	&& pC->getFace(2)->isValid()))
      vFatalError("Orientation of an internal boundary is incorrect.",
		  "2D initial meshing from boundary geometry.");
  }

  logMessage(4, "About to purge mesh\n");

  // Clean up and exit

  purgeAllEntities();
  if (3 * getNumCells() + getNumBdryFaces() + 2 * getNumIntBdryFaces()
      != 2 * getNumFaces()) {
    logMessage(
	1,
	"NumCells: %u, NumIntBdryFaces: %u, NumBdryFaces: %u, NumFaces: %u\n",
	getNumCells(), getNumIntBdryFaces(), getNumBdryFaces(), getNumFaces());
    vFatalError(
	"Incompatible number of cells, faces, and boundary faces.\nProbably a boundary is oriented incorrectly.",
	"2D initial meshing from boundary geometry.");
  }

  // Mark vertices on internal boundaries.
  logMessage(4, "Now marking vertices on internal boundaries\n");
  for (iFace = 0; iFace < getNumFaces(); iFace++) {
    Face *pF = getFace(iFace);
    if (pF->getFaceLoc() == Face::eBdryTwoSide) {
      logMessage(4, "  ... for Face %u\n", iFace);
      pF->getVert(0)->setType(Vert::eBdryTwoSide);
      pF->getVert(1)->setType(Vert::eBdryTwoSide);
    }
  }

  // Mark vertices on external boundaries.  Must be done after internal
  // boundaries.
  logMessage(4, "Now marking vertices on external boundaries...\n");
  for (iBFace = 0; iBFace < getNumBdryFaces(); iBFace++) {
    logMessage(4, "  ... for BFace %u\n", iBFace);
    Face *pF = getBFace(iBFace)->getFace(0);
    //      pF->pVVert(0)->vSetType(Vert::eBdryApex);
    //      pF->pVVert(1)->vSetType(Vert::eBdryApex);
    switch (pF->getVert(0)->getVertType())
      {
      case Vert::eBdryTwoSide:
	pF->getVert(0)->setType(Vert::eBdryApex);
	break;
      case Vert::eBdryApex:
	break;
      default:
	pF->getVert(0)->setType(Vert::eBdry);
	break;
      }
    switch (pF->getVert(1)->getVertType())
      {
      case Vert::eBdryTwoSide:
	pF->getVert(1)->setType(Vert::eBdryApex);
	break;
      case Vert::eBdryApex:
	break;
      default:
	pF->getVert(1)->setType(Vert::eBdry);
	break;
      }
  }

  logMessage(4, "Getting rid of list of edges to be recovered\n");
  for (int kk = LsTrackEdges.iLength() - 1; kk >= 0; kk--) {
    pEdge = LsTrackEdges[kk];
    delete pEdge;
  }
  LsTrackEdges.vClear();

  // Don't know if the following is necessary...?
  setAllHintFaces();
  setVertFaceNeighbors();

  assert(isValid());

  //   vMessage(3, "done.\n Swapping to make sure everything is Delaunay...\n");

  // try to swap
  int iRes = SM2D.swapAllFaces();
  logMessage(4, "Finished swapping, and it swapped %d cells.\n", iRes);

}

#endif

// This is an old constructor that is used only in the scat2d executable.

//@ Delaunay triangulation of a point set, with optional recovery of perimeter
Mesh2D::Mesh2D(const EntContainer<Vert> & EC) :
    Mesh(EC), m_ECEdgeF(), m_ECEdgeBF(), m_ECTri(), m_ECQuad(), m_ECIntEdgeBF()
{

  m_qual = new Quality(this, 2);

  GR_index_t iVert;

  // Set up an outer triangulation of a big square
  double dXMin, dXMax, dYMin, dYMax;
  dXMin = dYMin = 1.e100;
  dXMax = dYMax = -1.e100;

  for (iVert = 0; iVert < iNumVerts(); iVert++) {
    double dX = getVert(iVert)->x();
    double dY = getVert(iVert)->y();
    dXMin = min(dXMin, dX);
    dXMax = max(dXMax, dX);
    dYMin = min(dYMin, dY);
    dYMax = max(dYMax, dY);
  }

  double dXLo = dXMin - (dXMax - dXMin);
  double dXHi = dXMax + (dXMax - dXMin);
  double dYLo = dYMin - (dYMax - dYMin);
  double dYHi = dYMax + (dYMax - dYMin);

  GR_index_t iNumOrigVerts = iNumVerts();

  Vert *pVLoLeft = createVert(dXLo, dYLo);
  Vert *pVLoRight = createVert(dXHi, dYLo);
  Vert *pVHiLeft = createVert(dXLo, dYHi);
  Vert *pVHiRight = createVert(dXHi, dYHi);

  double adLoc[2];
  adLoc[0] = dXLo;
  adLoc[1] = dYLo;
  pVLoLeft->setCoords(2, adLoc);
  adLoc[1] = dYHi;
  pVHiLeft->setCoords(2, adLoc);
  adLoc[0] = dXHi;
  pVHiRight->setCoords(2, adLoc);
  adLoc[1] = dYLo;
  pVLoRight->setCoords(2, adLoc);

  // Set up initial mesh
  //
  //    HL-------0--------HR
  //     |              _/|
  //     |     0     __/  |
  //     |         _/     |
  //     1      _4/       3
  //     |    _/     1    |
  //     | __/            |
  //     |/               |
  //    LL-------2--------LR
  //

  bool qExist;
  Face *pF0 = createFace(qExist, pVHiRight, pVHiLeft);
  assert(!qExist);
  Face *pF1 = createFace(qExist, pVHiLeft, pVLoLeft);
  assert(!qExist);
  Face *pF2 = createFace(qExist, pVLoLeft, pVLoRight);
  assert(!qExist);
  Face *pF3 = createFace(qExist, pVLoRight, pVHiRight);
  assert(!qExist);
  Face *pF4 = createFace(qExist, pVLoLeft, pVHiRight);
  assert(!qExist);

  Cell *pC0 = createTriCell(pF0, pF1, pF4, iDefaultRegion);
  Cell *pC1 = createTriCell(pF2, pF3, pF4, iDefaultRegion);

  BFace *pBF0 = createBFace(pF0);
  BFace *pBF1 = createBFace(pF1);
  BFace *pBF2 = createBFace(pF2);
  BFace *pBF3 = createBFace(pF3);

  assert(pF0->doFullCheck());
  assert(pF1->doFullCheck());
  assert(pF2->doFullCheck());
  assert(pF3->doFullCheck());
  assert(pF4->doFullCheck());

  assert(pBF0->doFullCheck());
  assert(pBF1->doFullCheck());
  assert(pBF2->doFullCheck());
  assert(pBF3->doFullCheck());

  assert(pC0->doFullCheck());
  assert(pC1->doFullCheck());

  assert(isValid());

  GRUMMP::DelaunaySwapDecider2D Del2D;
  GRUMMP::SwapManager2D SM2D(&Del2D, this);
  GRUMMP::SwappingInserter2D SI2D(this, &SM2D);
  SI2D.setForcedInsertion(true);
  // Do Delaunay insertion of all the existing verts in the mesh
  setSwapType(eDelaunay);
  for (iVert = 0; iVert < iNumOrigVerts; iVert++) {
    Vert *pV = getVert(iVert);
    double adCoords[2];
    adCoords[0] = pV->x();
    adCoords[1] = pV->y();
    // Insert, then swap.  Force insertion at the specified point and
    // use the existing vertex location.
    Cell *pCGuess = getCell(0);
    for (int ii = 0; pCGuess->isDeleted(); ii++) {
      pCGuess = getCell(ii);
    }
#ifndef NDEBUG
    Vert* newVert = SI2D.insertPoint(adCoords, pCGuess, pV);
    assert(newVert->isValid());
#else
    (void) SI2D.insertPoint (adCoords, pCGuess, pV);
#endif
  }
  setAllHintFaces();
  setVertFaceNeighbors();
  assert(isValid());

  // No perimeter specified
  // Delete those four outer points and anything incident on them
  Vert *apV[] =
    { pVLoLeft, pVLoRight, pVHiLeft, pVHiRight };
  for (int i = 0; i < 4; i++) {
    std::set<Cell*> spCInc;
    std::set<Vert*> spVTmp;
    std::set<BFace*> spBFaceInc;
    bool isBdry = false;
    findNeighborhoodInfo(apV[i], spCInc, spVTmp, &spBFaceInc, &isBdry);
    for (std::set<Cell*>::iterator iter = spCInc.begin(); iter != spCInc.end();
	++iter) {
      Cell *pC = *iter;
      deleteCell(pC);
    }
    for (std::set<BFace*>::iterator iter = spBFaceInc.begin();
	iter != spBFaceInc.end(); ++iter) {
      BFace *pBF = *iter;
      deleteBFace(pBF);
    }
    deleteVert(apV[i]);
  }
  //
  //  pVLoLeft->vMarkDeleted ();
  //  pVLoRight->vMarkDeleted ();
  //  pVHiLeft->vMarkDeleted ();
  //  pVHiRight->vMarkDeleted ();
  //
  //  pBFBFace (0)->vMarkDeleted ();
  //  pBFBFace (1)->vMarkDeleted ();
  //  pBFBFace (2)->vMarkDeleted ();
  //  pBFBFace (3)->vMarkDeleted ();
  //  // Do this the easy way, in which every cell is checked for incidence
  //  // on one of the deleted verts, then local marking is done.
  //  GR_index_t iCell, iFace;
  //  for (iCell = 0; iCell < iNumCells (); iCell++) {
  //    Cell *pC = pCCell (iCell);
  //    if (pC->pVVert (0)->qDeleted () ||
  //	pC->pVVert (1)->qDeleted () ||
  //	pC->pVVert (2)->qDeleted ()) {
  //      Face *pF;
  //      // Find the face that is -not- incident on a deleted vert.  Set it
  //      // up as a new boundary face.
  //      for (iFace = 0; iFace < 3; iFace++) {
  //	pF = pC->pFFace (iFace);
  //	if (!pF->pVVert (0)->qDeleted () &&
  //	    !pF->pVVert (1)->qDeleted ()) {
  //	  break;
  //	}
  //      }
  //      deleteCell(pC);
  //      if (iFace != 3) createBFace(pF);
  //    }
  //  }
  //
  //  for (iFace = 0; iFace < iNumFaces (); iFace++) {
  //    Face *pF = pFFace (iFace);
  //    if (pF->qDeleted()) continue;
  //    if (pF->pVVert (0)->qDeleted () ||
  //	pF->pVVert (1)->qDeleted ())
  //      deleteFace(pF);
  //  }

  // Clean up and exit
  purgeAllEntities();
  assert(isValid());
}
