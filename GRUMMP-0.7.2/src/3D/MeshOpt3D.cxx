#include <set>
#include <map>

#include <stdlib.h>
#include "GR_ADT.h"
#include "GR_Geometry.h"
#include "GR_SmoothingManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_VolMesh.h"

extern int aiEdgeReq[], aiEdgeDone[];

////@ Optimization driver for face swapping.  Swap measure is set elsewhere.
int
VolMesh::iSwap_deprecated(const int iMaxPasses, const bool /* qAlreadyMarked */)
{
  int iSwaps = 1;
  assert(m_swapMeasure != eNoSwap);

  int iTotal = 0;
  int iPasses = 0;

  while (iSwaps && iPasses < iMaxPasses) {
    ++iPasses;
    if (iMaxPasses > 1)
      logMessage(2, "   Pass %d:", iPasses);
    iSwaps = 0;
    GR_index_t iMessageInt = std::max(getNumTriFaces() / 10 + 1,
				      GR_index_t(100000));
    for (int i = getNumTriFaces() - 1; i >= 0; i--) {
      if (i % iMessageInt == 0) {
	logMessage(2, "Still have %d faces to go this pass; %d swaps so far\n",
		   i, iSwaps);
      }
      Face *pF = getFace(i);
      // Don't swap if the face no longer exists or if it doesn't need
      // to be swapped.
      if (pF->isSwapAllowed()) {
	// Always swap recursively
	iSwaps += iFaceSwap_deprecated(pF);
      }
    }
    iTotal += iSwaps;
    logMessage(2, "%6d swaps", iSwaps);
    if (iMaxPasses > 1)
      logMessage(2, ",%6d total", iTotal);
    logMessage(2, "\n");
    purgeAllEntities();
  }
  if (iMaxPasses > 1)
    logMessage(2, "Face swapping complete.\n");
  for (int ii = 2; ii < 10; ii++) {
    logMessage(
	3,
	"    Edge swapping from %d tets to %2d: %5d done out of %5d requests\n",
	ii, 2 * ii - 4, aiEdgeDone[ii], aiEdgeReq[ii]);
  }

  if (!(isValid()))
    vFoundBug("3D mesh reconnection");
  return iTotal;
}

static bool
qEdgeHitsBFace(const Vert * const pV0, const Vert * const pV1,
	       const ADT& BFaceTree, const VolMesh& VM)
{
  double adBBox[] =
    { min(pV0->x(), pV1->x()), max(pV0->x(), pV1->x()), min(pV0->y(), pV1->y()),
	max(pV0->y(), pV1->y()), min(pV0->z(), pV1->z()), max(pV0->z(),
							      pV1->z()) };
  std::vector<GR_index_t> veciQueryResult = BFaceTree.veciRangeQuery(adBBox);
  if (veciQueryResult.size() >= 20) {
    logMessage(3, "Warning: query returned %d possible bdry faces.\n",
	       static_cast<int>(veciQueryResult.size()));
  }

  for (int ii = veciQueryResult.size() - 1; ii >= 0; ii--) {
    Face *pF = VM.getBFace(veciQueryResult[ii])->getFace();
    // Now check whether the edge actually penetrates the face.  If so,
    // return true.  (Actually, return this even if the edge hits an
    // edge or corner of the face.)  For the face to be okay (edge
    // outside) then there must be one orientation of each sign (one
    // plus and one minus).  This is so regardless of which sign is the
    // inside orientation and which is the outside orientation.
    bool qHavePlus = false, qHaveMinus = false;
    switch (checkOrient3D(pF->getVert(0), pF->getVert(1), pV0, pV1))
      {
      case 1:
	qHavePlus = true;
	break;
      case -1:
	qHaveMinus = true;
	break;
      default:
	break;
      }
    switch (checkOrient3D(pF->getVert(1), pF->getVert(2), pV0, pV1))
      {
      case 1:
	qHavePlus = true;
	break;
      case -1:
	qHaveMinus = true;
	break;
      default:
	break;
      }
    switch (checkOrient3D(pF->getVert(2), pF->getVert(0), pV0, pV1))
      {
      case 1:
	qHavePlus = true;
	break;
      case -1:
	qHaveMinus = true;
	break;
      default:
	break;
      }
    if (!(qHavePlus && qHaveMinus))
      return true;
  }
  return false;
}

static bool
qIsVisible(const Cell * const pC, const Vert * const pV, const ADT& BFaceTree,
	   const VolMesh& VM)
{
  assert(pC->getType() == Cell::eTet);
  // Visibility is defined here as each vertex of the tet being visible
  // to the vertex inside the circumsphere.  Technically, one could have
  // situations where part of an edge was obscured, for example, but not
  // any vertex.  The general case would be, I think, very difficult to
  // test correctly.
  for (int iV = 0; iV < 4; iV++) {
    const Vert *pVSource = pC->getVert(iV);
    if (qEdgeHitsBFace(pV, pVSource, BFaceTree, VM))
      return false;
  }
  return true;
}

int
identifyNonDelaunayCells(VolMesh& VM, bool aqDelaunay[])
{
  int iRetVal = 0;

  // Create a vertex tree so that insphere searching can be done
  // efficiently.
  ADT VT(VM.pECGetVerts());

  // Also create a tree of bdry faces to make visibility checking
  // easier.
  double (*a2dBFace)[6] = new double[VM.getNumBdryFaces()][6];
  for (GR_index_t iBFace = 0; iBFace < VM.getNumBdryFaces(); iBFace++) {
    BFace *pBF = VM.getBFace(iBFace);
    if (pBF->isDeleted())
      continue;

    a2dBFace[iBFace][0] = MIN3(pBF->getVert(0)->x(), pBF->getVert(1)->x(),
			       pBF->getVert(2)->x());
    a2dBFace[iBFace][1] = MAX3(pBF->getVert(0)->x(), pBF->getVert(1)->x(),
			       pBF->getVert(2)->x());
    a2dBFace[iBFace][2] = MIN3(pBF->getVert(0)->y(), pBF->getVert(1)->y(),
			       pBF->getVert(2)->y());
    a2dBFace[iBFace][3] = MAX3(pBF->getVert(0)->y(), pBF->getVert(1)->y(),
			       pBF->getVert(2)->y());
    a2dBFace[iBFace][4] = MIN3(pBF->getVert(0)->z(), pBF->getVert(1)->z(),
			       pBF->getVert(2)->z());
    a2dBFace[iBFace][5] = MAX3(pBF->getVert(0)->z(), pBF->getVert(1)->z(),
			       pBF->getVert(2)->z());
  }
  ADT BFaceTree(VM.getNumBdryFaces(), 3, ADT::eBBoxes,
		reinterpret_cast<double*>(a2dBFace));
  delete[] a2dBFace;

  // Now check each cell to see how many verts fall into its
  // circumsphere.  Keep track of the max as well.
  int iMaxEncroach = 0;
  int *aiEncroach = new int[VM.getNumCells()];
  std::vector<GR_index_t> veciQueryResult;
  double dEps = 1.e-10;
  GR_index_t iCell;
  for (iCell = 0; iCell < VM.getNumCells(); iCell++) {
    Cell *pCCell = VM.getCell(iCell);
    if (pCCell->isDeleted())
      continue;
    assert(pCCell->getType() == Cell::eTet);
    TetCell *pTC = dynamic_cast<TetCell*>(pCCell);
    double adCircCent[3], dCircRad;
    dCircRad = pTC->calcCircumradius();
    pTC->calcCircumcenter(adCircCent);
    // Get all verts in the bounding box of the circumsphere
    double adRange[] =
      { adCircCent[0] - dCircRad - dEps, adCircCent[0] + dCircRad + dEps,
	  adCircCent[1] - dCircRad - dEps, adCircCent[1] + dCircRad + dEps,
	  adCircCent[2] - dCircRad - dEps, adCircCent[2] + dCircRad + dEps };
    veciQueryResult = VT.veciRangeQuery(adRange);
    assert(veciQueryResult.size() >= 4);
    Vert *pV0 = pTC->getVert(0);
    Vert *pV1 = pTC->getVert(1);
    Vert *pV2 = pTC->getVert(2);
    Vert *pV3 = pTC->getVert(3);
    int iNumEnc = 0;
    for (int ii = veciQueryResult.size() - 1; ii >= 0; ii--) {
      Vert *pV = VM.getVert(veciQueryResult[ii]);
      if ((pV == pV0) || (pV == pV1) || (pV == pV2) || (pV == pV3))
	continue;
      int iIsEncroached = isInsphere(pTC->getVert(0), pTC->getVert(1),
				     pTC->getVert(2), pTC->getVert(3), pV);
#ifndef NDEBUG
      if (iIsEncroached != 1) {
	// On or outside the sphere
	const double *adLoc = pV->getCoords();
	double dTmpDist = dDIST3D(adLoc, adCircCent);
	assert(iFuzzyComp(dTmpDist, dCircRad) >= 0);
      }
#endif
      if (iIsEncroached == 1 && qIsVisible(pTC, pV, BFaceTree, VM)) {
	iNumEnc++;
      }
    }
    aiEncroach[iCell] = iNumEnc;
    if (iNumEnc == 0) {
      aqDelaunay[iCell] = true;
    }
    else {
      iRetVal++;
    }
    if (iMaxEncroach <= iNumEnc) {
      iMaxEncroach = iNumEnc;

    }
//    iMaxEncroach = max(iNumEnc, iMaxEncroach);
  }
  logMessage(1, "Maximum number of extra verts in a circumsphere is %d\n",
	     iMaxEncroach);
  logMessage(1, "Num verts encroaching      Num tets\n");
  int *aiCount = new int[iMaxEncroach + 1];
  for (int ii = 0; ii < iMaxEncroach + 1; aiCount[ii++] = 0) {
  }
  for (iCell = 0; iCell < VM.getNumCells(); iCell++)
    aiCount[aiEncroach[iCell]]++;
  for (int ii = 0; ii <= iMaxEncroach; ii++)
    if (aiCount[ii] != 0)
      logMessage(1, "%12d%20d\n", ii, aiCount[ii]);
  delete[] aiEncroach;
  delete[] aiCount;
  return (iRetVal);
}

bool
VolMesh::makeDelaunay()
{
  assert(isSimplicial());
  // First, do insphere face swapping; this often gives a Delaunay mesh
  // right away.
  GRUMMP::SwapDecider3D *pSD3D = new GRUMMP::DelaunaySwapDecider3D(
      areSurfaceEdgeChangesAllowed(), getMaxAngleForSurfaceEdgeChanges());
  GRUMMP::SwapManager3D Swapper(pSD3D, this);
  (void) Swapper.swapAllFaces(false);
  purgeAllEntities();

  GRUMMP::OptMSSmoothingManager3D Smoother(this);

  int iOldNumNonDelaunay = 0, iNumNonDelaunay;
  // Identify non-Delaunay tets.
  bool *aqDelaunay = new bool[getNumCells()];
  iNumNonDelaunay = identifyNonDelaunayCells(*this, aqDelaunay);

  while (iNumNonDelaunay != 0 && iNumNonDelaunay != iOldNumNonDelaunay) {
    // Tag all verts as not in need of smoothing; for those in
    // non-Delaunay tets, this will change.

    // Tag (only) the verts of non-Delaunay tets for smoothing.
    for (GR_index_t iCell = 0; iCell < getNumCells(); iCell++) {
      if (!aqDelaunay[iCell]) {
	Cell *pC = getCell(iCell);
	moveVertEvent(pC->getVert(0));
	moveVertEvent(pC->getVert(1));
	moveVertEvent(pC->getVert(2));
	moveVertEvent(pC->getVert(3));
      }
    }
//    Swapper.parallelswapping = false;
    sendEvents();

    // Smooth, then swap only the affected faces.
    logMessage(2, "Smoothing...\n");
    Smoother.smoothAllQueuedVerts();
    logMessage(2, "Swapping...\n");
    (void) Swapper.swapAllQueuedFaces();

    // Identify non-Delaunay tets.
    iOldNumNonDelaunay = iNumNonDelaunay;
    iNumNonDelaunay = identifyNonDelaunayCells(*this, aqDelaunay);
    logMessage(2, "Had %d tets non-Delaunay last time; now %d\n",
	       iOldNumNonDelaunay, iNumNonDelaunay);

    // While still making progress, repeat.
  }

  GR_index_t iFace;
  Vert *pVVertA, *pVVertB, *pVVertC, *pVVertD, *pVVertE;
  Vert *pVPivot0, *pVPivot1, *pVOther;
  TetCell* apTCTets[4];
  int iNTets;
  GR_index_t iN44Cases, iAttempts, iSuccesses, iPasses = 0;
  do {
    iN44Cases = iAttempts = iSuccesses = 0;
    iPasses++;
    for (iFace = 0; iFace < getNumFaces(); iFace++) {
      Face *pF = getFace(iFace);
      if (pF->isDeleted() || pF->isBdryFace())
	continue;
      int iCL = getCellIndex(pF->getLeftCell());
      int iCR = getCellIndex(pF->getRightCell());
      if (aqDelaunay[iCL] || aqDelaunay[iCR])
	continue;
      assert(pF->getType() == Face::eTriFace);
      eFaceCat eFC = dynamic_cast<TriFace*>(pF)->categorizeFace(pVVertA,
								pVVertB,
								pVVertC,
								pVVertD,
								pVVertE,
								apTCTets,
								iNTets,
								pVPivot0,
								pVPivot1,
								pVOther);
      if (eFC == eOther || eFC == eBdry)
	continue;
      int iPointInCircsphere = isInsphere(pVVertA, pVVertB, pVVertC, pVVertD,
					  pVVertE);
      switch (eFC)
	{
	default:
	  break;
	case eN44:
	case eN32:
	  if (eFC == eN44)
	    iN44Cases++;
	  if (iPointInCircsphere == 1) {
	    assert(
		pVVertD != pVPivot0 && pVVertD != pVPivot1
		    && pVVertD != pVOther);
	    assert(
		pVVertE != pVPivot0 && pVVertE != pVPivot1
		    && pVVertE != pVOther);
	    logMessage(
		2,
		"Attempting to eliminate edge (%4u %4u) in favor of (%4u %4u)\n",
		getVertIndex(pVPivot0) + 1, getVertIndex(pVPivot1) + 1,
		getVertIndex(pVVertD) + 1, getVertIndex(pVVertE) + 1);
	    int iRes = edgeSwap3D(pF, pVPivot0, pVPivot1, pVOther, true,
				  pVVertD, pVVertE);
	    iAttempts++;
	    if (iRes > 0) {
	      iSuccesses++;
	    }
	  }
	  break;
	}
    }
    logMessage(2, "Number of N44 cases attempted: %d\n", iN44Cases);
    logMessage(
	2,
	"Succeeded on %d of %d attempts to remove pivot edges in N32 and N44 faces\n",
	iSuccesses, iAttempts);
  }
  while (iSuccesses > 0 && iPasses < 10);
  purgeAllEntities();

  delete[] aqDelaunay;
  return (iNumNonDelaunay == 0);
}
