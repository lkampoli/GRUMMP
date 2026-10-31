#include <set>

#include "GR_CellCV.h"
#include "GR_Geometry.h"
#include "GR_VolMesh.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"

class EdgeToBFace;

struct RemCell {
  Cell *pC;
  Face *pFOnlyV, *pFOnlyV1;
  RemCell() :
      pC(pCInvalidCell), pFOnlyV(pFInvalidFace), pFOnlyV1(pFInvalidFace)
  {
  }
  RemCell(Cell * const pCIn, Face * const pF0, Face * const pF1) :
      pC(pCIn), pFOnlyV(pF0), pFOnlyV1(pF1)
  {
  }
};

static int iDelReq = 0, iBdryDelReq = 0, iBdryCurveDelReq = 0;
static int iBdryApexDelReq = 0;
static int iDelDone = 0, iBdryDelDone = 0, iBdryCurveDelDone = 0;

// static int iStarRemeshFailed = 0, iStarRemeshSucceeded = 0,
//   iStarRestoreFailed = 0;

// class FakeTet : public TetCell {
// private:
//   /// Copy construction disallowed.
//   FakeTet(const FakeTet&) : TetCell() {assert(0);};
//   Vert *apV[4];
// public:
//   eCellType eType() const {return eFakeTet;};
//   FakeTet(Vert* const pV0, Vert* const pV1,
//           Vert* const pV2, Vert* const pV3) : TetCell(NULL)
//     {
//        apV[0] = pV0; apV[1] = pV1; apV[2] = pV2; apV[3] = pV3;
//     };
//   Vert * pVVert(const int i) const {assert(i>=0 && i<4); return apV[i];};
// };

// To remove a given vertex V from a mesh, there must be a vertex V1
// adjacent to V which is visible from all other vertices adjacent to
// V.  (Note that this specifically does NOT require that the
// neighborhood be convex.)  If this requirement is met, then all faces
// and all cells incident on both V and V1 can be removed from the mesh.
// The faces in the deleted cells which are incident on only one of V
// and V1 occur in pairs and can be matched up and reconnected to ensure
// a valid mesh on output.

bool
VolMesh::removeVert(Vert * const pV, int& iNewSwaps,
				 Vert * apVPreferred[], const int iNPref)
// Remove vertex *pV from the mesh.  Decisions about whether *pV
// SHOULD be removed are made elsewhere.  The vertices listed in
// apVPreferred are given a weak preference in determining to which
// vertex pV should be moved.  This preference takes the form of
// checking these vertices first in making removal choices.
{

  // Only works for tetrahedral meshes at this point.
  iNewSwaps = 0;
  iDelReq++;

  if (pV->getVertType() == Vert::eBdryApex) {
    logMessage(3, "Deletion of a boundary apex requested! Denied.\n");
    iBdryApexDelReq++;
    return false;
  }

  //  if (pV->getVertType() == Vert::eBdryTwoSide) {
  //    logMessage(4, "Deletion of verts on internal boundaries not yet supported.\n");
  //    return false;
  //  }

  // List all incident cells on *pV and all vertices in those cells.  As
  // a side effect, find out whether this is a boundary vertex.
  std::set<Cell*> spCIncident;
  std::set<Vert*> spVCand;
  std::set<BFace*> spBFIncident;
  std::set<Face*> spFUpdate;

  std::set<Vert*>::iterator iterV;
  std::set<Cell*>::iterator iterC;
  std::set<Face*>::iterator iterF;
  std::set<BFace*>::iterator iterBF;

  bool qBdryVert = false;
  bool qBdryCurveVert = false;

  findNeighborhoodInfo(pV, spCIncident, spVCand, &spBFIncident, &qBdryVert,
		       &spFUpdate);

  for (Cell *pC : spCIncident) {
    if (pC->getType() != CellSkel::eTet) {
      logMessage(3, "Can't delete a vertex with a non-tet incident\n");
      return false;
    }
  }

  std::set<Vert*> spVAllBdryCands;
  if (qBdryVert) {
    // This case shouldn't really happen, but it sometimes does anyway,
    // and it's fatal.
    if (pV->getVertType() == Vert::eInterior)
      return false;

    // BdryApex calls have already been bounced.
    assert(
	pV->getVertType() == Vert::eBdryCurve
	    || pV->getVertType() == Vert::eBdryTwoSide
	    || pV->getVertType() == Vert::eBdry);
    // Abort if boundaries are precious
    if (!areBdryChangesAllowed())
      return false;

    // Permanently disqualify all vertices which don't meet boundary rules
    // for appropriate contraction.  Those rules are:
    //   1.  BdryApex verts never get removed.
    //   2.  BdryCurve verts only go onto BdryCurve and BdryApex verts
    //       that fall on the same curve (i.e., the connecting edge must
    //       lie on a curve, so the two BFaces shared by the two verts
    //       must lie on different patches.
    //   3.  Bdry verts that are structured only go onto BdryCurve or
    //       BdryApex verts, or onto kept, structured Bdry verts.
    //   4.  Other Bdry verts can go onto any flavor of bdry vert.

    // First, find all of the faces incident on pV (faces in spFUpdate)
    // that are boundary faces.
    spVCand.clear();
    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
	iterBF++) {
      BFace *pBF = *iterBF;
      Face *pFBdry = pBF->getFace();
      assert(pFBdry->hasVert(pV));
      spVCand.insert(pFBdry->getVert(0));
      spVCand.insert(pFBdry->getVert(1));
      spVCand.insert(pFBdry->getVert(2));
    }
    spVCand.erase(pV);
    //    assert(spVCand.size() == spBFIncident.size());

    spVAllBdryCands.insert(spVCand.begin(), spVCand.end());

    if (iNPref > 0) {
      // If there are any preferred verts, they will be the only ones
      // considered as contraction targets.  If there are none, proceed as
      // usual.
      spVCand.clear();
      spVCand.insert(apVPreferred, apVPreferred + iNPref);
    }

    switch (pV->getVertType())
      {
      case Vert::eBdryCurve:
	{
	  qBdryCurveVert = true;
//			logMessage(1, "Bdry curve: %lu candidates", spVCand.size());
	  for (iterV = spVCand.begin(); iterV != spVCand.end();) {
	    Vert *pVCand = *iterV;
	    std::set<Vert*>::iterator iterTmp = iterV;
	    ++iterV;
	    // Eliminate all verts that aren't bdry curves or bdry apexes.
	    logMessage(1, " %d", pVCand->getVertType());
	    if (pVCand->getVertType() != Vert::eBdryApex
		&& pVCand->getVertType() != Vert::eBdryCurve) {
	      spVCand.erase(iterTmp);
	    }
//				} else if (hasSubsegs()
//						&& areSubsegVertsConnected(pVCand, pV) == NULL) {
//					spVCand.erase(iterTmp);
//				}

	    // Now for each remaining candidate, find the pair of bdry faces
	    // that connect it to the vertex to be removed, and eliminate
	    // any candidates that aren't on the same bdry curve.

	    else {
	      BFace* apBFTemp[20];
	      int ii = 0;
	      for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
		  iterBF++) {
		BFace *pBFCand = *iterBF;
		assert(pBFCand->hasVert(pV));
		if (pBFCand->hasVert(pVCand))
		  apBFTemp[ii++] = pBFCand;
	      }
	      if (ii != 2) {
		spVCand.erase(iterTmp);
	      }
	      else {
		assert(ii == 2);
		BFace *pBF0 = apBFTemp[0];
		BFace *pBF1 = apBFTemp[1];
		double adNorm0[3], adNorm1[3];
		pBF0->calcUnitNormal(adNorm0);
		pBF1->calcUnitNormal(adNorm1);
		double dDot = dDOT3D(adNorm0, adNorm1);
		if (dDot > cos(m_maxAngleForSurfSwap * M_PI / 180)) {
		  spVCand.erase(iterTmp);
		}
	      }
	    } // Done with the case where the candidate vertex is on a
	      // bdry curve somewhere
	  } // Done checking all candidate verts when the vert to be
	    // removed lies on a bdry

	} // end case
	logMessage(1, " %lu candidates left\n", spVCand.size());
	iBdryCurveDelReq++;
	break;
      case Vert::eBdryTwoSide:
      case Vert::eBdry:
	{
	  // Any candidate verts must be connected to the vert to be
	  // removed by an edge that actually lies on the boundary.
	  for (iterV = spVCand.begin(); iterV != spVCand.end();) {
	    Vert *pVCand = *iterV;
	    std::set<Vert*>::iterator iterTmp = iterV;
	    ++iterV;
	    // For each candidate, find the pair of bdry faces that
	    // connect it to the vertex to be removed, and eliminate any
	    // candidates that aren't on the same bdry curve.
	    int ii = 0;
	    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
		iterBF++) {
	      BFace *pBFCand = *iterBF;
	      assert(pBFCand->hasVert(pV));
	      if (pBFCand->hasVert(pVCand))
		ii++;
	    }
	    if (ii != 2) {
	      logMessage(
		  2,
		  "Blocking possible attempt to remove a bdry vert along a bad edge.\n");
	      logMessage(2, "# bdry faces shared by verts: %d\n", ii);
	      spVCand.erase(iterTmp);
	    }
	  } // Done checking all candidate verts when the vert to be
	    // removed lies on a bdry

	    // If the vertex is structured, then it should only be slid onto
	    // another structured vertex or a vert on a curve or apex.
	  if (pV->isStructured()) {
	    for (iterV = spVCand.begin(); iterV != spVCand.end();) {
	      Vert *pVCand = *iterV;
	      std::set<Vert*>::iterator iterTmp = iterV;
	      ++iterV;
	      if (!(pVCand->isStructured()
		  || pVCand->getVertType() == Vert::eBdryApex
		  || pVCand->getVertType() == Vert::eBdryCurve)) {
		spVCand.erase(iterTmp);
	      }
	    }
	  }
	}
	iBdryDelReq++;
	break;
      case Vert::eBdryApex:
      default:
	// Should never get here.
	assert(0);
	break;
      }
  }

  //  assert(spVCand.size() == size_t(iNPref) || iNPref == 0);

  // List all faces in incident cells that are opposite *pV
  std::vector<Face*> vecpFOpposite;
  for (iterC = spCIncident.begin(); iterC != spCIncident.end(); iterC++) {
    Cell *pC = *iterC;
    assert(pC->getType() == Cell::eTet);
    assert(pC->hasVert(pV));
    Face *pFOpp = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV);
    vecpFOpposite.push_back(pFOpp);
    spFUpdate.erase(pFOpp);
  }

  //   // Used only in assertions just before the actual deletion.
  // #ifndef NDEBUG
  //   int iTotalIncFaces = spFUpdate.size();
  // #endif

  // For each neighboring vertex, determine whether it is a valid target
  // for edge contraction.  Keep track of which faces affect which
  // vertices.  If no vertices are okay, start with the one with the
  // fewest interfering faces and try to remove that face from the hull
  // by using edge removal to break the connection from *pV to any of
  // the vertices of the face.  If that succeeds, call this routine
  // recursively to do the contraction; otherwise continue until an
  // attempt has been made to remove each problematic face.

  bool qVertOK;
  bool qAnyVertOK = false;
  double dMinDihed;
  double dMaxMinDihed = 0.; // Force a minimum dihedral angle before
  // attempting to remove a vertex.
  Vert *pV1 = pVInvalidVert, *pVBest = pVInvalidVert;
  for (iterV = spVCand.begin(); iterV != spVCand.end(); ++iterV) {
    pV1 = *iterV;
    qVertOK = true;
    dMinDihed = 180;
    for (unsigned int iF = 0; iF < vecpFOpposite.size() && qVertOK; iF++) {
      Face *pF = vecpFOpposite[iF];
      if (!pF->hasVert(pV1)) {
	// Determine whether a given face has the correct orientation
	// for the candidate contraction target.  If the orientation
	// test fails, pV1 can't see the correct side of the face pF.
	bool qResult;
	if ((pF->getRightCell())->hasVert(pV)) {
	  // For this case, the right-hand rule applies for orientation.
	  assert(checkOrient3D(pF, pV) == 1);
	  qResult = (checkOrient3D(pF, pV1) == 1);
	  if (qResult) {
	    Vert *apV[] =
	      { pF->getVert(0), pF->getVert(1), pF->getVert(2), pV1 };
	    TetCellCV FT(apV);
	    double dMin = FT.calcMinDihed();
	    double dMax = FT.calcMaxDihed();
	    dMinDihed = min(dMinDihed, min(dMin, 180 - dMax));
	    logMessage(4, "Min dihed (face): %10.6f  (vert): %10.6f\n",
		       FT.calcMinDihed(), dMinDihed);
	  }
	}
	else {
	  // Otherwise, the left-hand rule applies for orientation.
	  assert((pF->getLeftCell())->hasVert(pV));
	  assert(checkOrient3D(pF, pV) == -1);
	  qResult = (checkOrient3D(pF, pV1) == -1);
	  if (qResult) {
	    Vert *apV[] =
	      { pF->getVert(1), pF->getVert(0), pF->getVert(2), pV1 };
	    TetCellCV FT(apV);
	    double dMin = FT.calcMinDihed();
	    double dMax = FT.calcMaxDihed();
	    dMinDihed = min(dMinDihed, min(dMin, 180 - dMax));
	    logMessage(4, "Min dihed (face): %10.6f  (vert): %10.6f\n",
		       FT.calcMinDihed(), dMinDihed);
	  }
	}
	if (!qResult) {
	  qVertOK = false;
	}
	if (dMinDihed < dMaxMinDihed) {
	  qVertOK = false;
	}
      } // Done checking a face
    } // Done checking a vertex

    if (qVertOK && qBdryVert) {

      // For bdry verts, make certain that there aren't three (or
      // more) boundary verts shared by pV and pV1.  If there are,
      // then contracting this edge will collapse two boundary faces
      // onto each other.  While having BFaces collapse
      // "front-to-front" (interior sides together) is recoverable,
      // the "back-to-back" collapse is not recoverable.  And even the
      // front-to-front collapse may potentially delete a vertex that
      // wasn't tagged for deletion, so it isn't a great idea anyway.

      // This check is made easier because spVCand already holds a
      // list of all bdry verts incident on pV.

      // May have to handle special case of three, with one internal boundary face
      // and two external ones
      std::set<Vert*> spVNeigh1;
      std::set<Cell*> spCJunk;
      findNeighborhoodInfo(pV1, spCJunk, spVNeigh1);

      assert(spVNeigh1.count(pV) == 1);
      assert(spVCand.count(pV1) == 1);
      assert(spVCand.count(pV) == 0);

      int iCommon = 0;
      for (std::set<Vert*>::iterator iterV2 = spVNeigh1.begin();
	  iterV2 != spVNeigh1.end(); iterV2++) {
	Vert *pVTest = *iterV2;
	if (pVTest->isBdryVert() && spVAllBdryCands.count(pVTest) == 1)
	  iCommon++;
      }
      switch (iCommon)
	{
	case 0:
	case 1:
	  assert(0);
	  vFoundBug("bizarre result for bdry vertex removal");
	  break;
	case 2:
	  // Everything's okay.
	  break;
	case 3:
	  {
	    // check specific case of two external boundary faces
	    // one internal boundary faces

	    int iNumBdryFaces = 0, iNumIntBdryFaces = 0;

	    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
		iterBF++) {
	      BFace *pBF = *iterBF;

	      assert(pBF->hasVert(pV));
	      if (pBF->hasVert(pV1)) {
		if (pBF->getType() == Cell::eTriBFace)
		  iNumBdryFaces++;
		if (pBF->getType() == Cell::eIntTriBFace)
		  iNumIntBdryFaces++;
	      }
	    }
	    if (iNumBdryFaces != 2 || iNumIntBdryFaces != 1)
	      qVertOK = false;

	  }
	  break;
	default:
	  // Too many bdry verts incident!
	  logMessage(
	      4,
	      "%d common bdry vert neighbors in bdry edge contraction; 2 is max allowed\n",
	      iCommon);
	  qVertOK = false;
	  break;
	} // Done with switch
    } // Done with check for weird bdry circumstances.

    if (qVertOK) {
      if (dMinDihed > dMaxMinDihed) {
	// The following line has been moved inside the 'if' so that the
	// oddball case of adMinDihed[] == 0 can't cause a problem.
	// Hey, it happened once!
	qAnyVertOK = true;
	logMessage(4, "Best min dihedral yet.  Old: %10.6f  New: %10.6f\n",
		   dMaxMinDihed, dMinDihed);
	dMaxMinDihed = dMinDihed;
	pVBest = pV1;
      }
      if (dMinDihed > 0) {
	// Check to see whether this is a preferred vertex
	for (int iPref = 0; iPref < iNPref; iPref++) {
	  if (apVPreferred[iPref] == pV1) {
	    pVBest = pV1;
	    break;
	  }
	}
      } // Done with preferred vert check
    } // Done updating best vertex
  } // Done checking all candidate targets for edge contraction

  if (!qAnyVertOK) {
    return false;
  }

  // If we get to here, we know we CAN remove *pV by contracting edge
  // (*pV, *pVBest), and that this gives the largest minimum dihedral
  // angle of any possible contraction.
  assert(qAnyVertOK);
  assert(pVBest->isValid());
  // If the best you can do is really lousy, don't do it.  This doesn't
  // apply for structured verts, where traditional measures of mesh
  // quality are (at best) difficult to interpret.
  // Change angle threshhold here!
  if (dMaxMinDihed < 3. && !pV->isStructured() && iNPref == 0)
    return false;
  logMessage(3, "Contracting vert %5u onto %5u; min dihed after: %5.2f deg.\n",
	     getVertIndex(pV), getVertIndex(pVBest), dMaxMinDihed);
  logMessage(4, "Contracting vert %5u at (%10g, %10g, %10g)\n",
	     getVertIndex(pV), pV->x(), pV->y(), pV->z());
  logMessage(4, "       onto vert %5u at (%10g, %10g, %10g)\n",
	     getVertIndex(pVBest), pVBest->x(), pVBest->y(), pVBest->z());

  // At this point, we're sure that we want to do the removal, and what
  // vertex the target is.

  std::set<EdgeToBFace> sEBFSingleEdges;
  if (qBdryVert) {
    std::multiset<EdgeToBFace> sEBFAllEdges;
    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
	iterBF++) {
      BFace *pBF = *iterBF;
      assert(pBF->hasVert(pV));
      if (!pBF->hasVert(pVBest)) {
	// Need to keep track of edges that will combine with pVBest to
	// create a new BFace after contraction.
	Face *pF = pBF->getFace(0);

	std::pair<std::set<EdgeToBFace>::iterator, bool> pairRes;
	if (pF->getVert(0) == pV) {
	  pairRes = sEBFSingleEdges.insert(
	      EdgeToBFace(pF->getVert(1), pF->getVert(2), pBF));
	  assert(pairRes.second);
	}
	else if (pF->getVert(1) == pV) {
	  pairRes = sEBFSingleEdges.insert(
	      EdgeToBFace(pF->getVert(2), pF->getVert(0), pBF));
	  assert(pairRes.second);
	}
	else {
	  assert(pF->getVert(2) == pV);
	  pairRes = sEBFSingleEdges.insert(
	      EdgeToBFace(pF->getVert(0), pF->getVert(1), pBF));
	  assert(pairRes.second);
	}
      }
    }
  }

  // Before it's too late, we need to make a note of what the actual
  // regions are for faces in the multiregion case.
  std::map<Face*, int> mFaceRegion;

  std::vector<Face*>::iterator iterFVec;
  for (iterFVec = vecpFOpposite.begin(); iterFVec != vecpFOpposite.end();
      iterFVec++) {
    Face *pF = *iterFVec;
    int iThisReg = iInvalidRegion;

    Cell *pCL = pF->getLeftCell();
    Cell *pCR = pF->getRightCell();
    assert(pCL->isValid() && pCR->isValid());
    assert(pCL->hasVert(pV) || pCR->hasVert(pV));
    if (pCL->getType() == CellSkel::eTet) {
      iThisReg = pCL->getRegion();
      assert(pCR->getType() != CellSkel::eTet || iThisReg == pCR->getRegion());
    }
    else {
      assert(pCR->getType() == CellSkel::eTet);
      iThisReg = pCR->getRegion();
    }
    mFaceRegion.insert(std::pair<Face*, int>(pF, iThisReg));
  }

  // Delete all the cells.
  for (iterC = spCIncident.begin(); iterC != spCIncident.end(); iterC++) {
    deleteCell(*iterC);
  }
  // Add new cells.
  std::map<Face*, int>::iterator mIter, mEnd = mFaceRegion.end();
  std::set<Face*> spFToCheck;
  for (mIter = mFaceRegion.begin(); mIter != mEnd; mIter++) {
    Face *pF = mIter->first;
    // There are some faces that need to be skipped here.
    if (pF->hasVert(pVBest))
      continue;
    assert(!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid());
    // Get the appropriate region...
    int iThisReg = mIter->second;
    assert(iThisReg != iInvalidRegion);
    bool qExist;
//		for(int iV = 0; iV < 3; iV++){
//			pF->getVert(iV)->setHintFace(pF);
//		}
    Cell *pC = createTetCell(qExist, pVBest, pF, iThisReg);
//		for(int iF = 0; iF < 4; iF++){
//			if(pC->getFace(iF) != pF){
//				pVBest->setHintFace(pC->getFace(iF));
//				break;
//			}
//		}

    assert(!qExist);
    spFToCheck.insert(pC->getFace(0));
    spFToCheck.insert(pC->getFace(1));
    spFToCheck.insert(pC->getFace(2));
    spFToCheck.insert(pC->getFace(3));

#ifndef NDEBUG
    assert(pC->doFullCheck());
    assert(
	checkOrient3D(pC->getVert(0), pC->getVert(1), pC->getVert(2),
		      pC->getVert(3)) == 1);
#endif

  }

  std::set<EdgeToBFace>::iterator iterEBF, endEBF = sEBFSingleEdges.end();
  for (iterEBF = sEBFSingleEdges.begin(); iterEBF != endEBF; iterEBF++) {
    EdgeToBFace EBF = *iterEBF;
    Face *pF = findCommonFace(EBF.pVVert(0), EBF.pVVert(1), pVBest, NULL, true);
    assert(pF->isValid());
    BFace *pBF = createBFace(pF, EBF.pBFace());
    assert(pBF->doFullCheck());
  }

  // Now it's finally safe to delete the old BFace's.
  for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end(); iterBF++) {
    BFace* pBF = *iterBF;
    deleteBFace(pBF);
  }

  // Delete *pV.
  deleteVert(pV);

#ifndef NDEBUG
  // Here's a post-check to be sure that bdry faces haven't collapsed
  // onto each other.  For now, it just spews output.
  if (qBdryVert) {
    std::set<Cell*> spCJunk;
    std::set<Vert*> spVJunk;
    std::set<BFace*> spBFInc;
    bool qBdryChk = false;
    findNeighborhoodInfo(pVBest, spCJunk, spVJunk, &spBFInc, &qBdryChk);

    assert(qBdryChk);
    spVJunk.clear();
    logMessage(4, "pVBest has %zd incident BFaces.\n", spBFInc.size());
    for (iterBF = spBFInc.begin(); iterBF != spBFInc.end(); iterBF++) {
      BFace *pBF = *iterBF;
      logMessage(4, "  Bdryface:%6u.  Verts:%7u %7u %7u\n", getBFaceIndex(pBF),
		 getVertIndex(pBF->getVert(0)), getVertIndex(pBF->getVert(1)),
		 getVertIndex(pBF->getVert(2)));
      spVJunk.insert(pBF->getVert(0));
      spVJunk.insert(pBF->getVert(1));
      spVJunk.insert(pBF->getVert(2));
    }
    spVJunk.erase(pVBest);
    //    assert(spVJunk.size() == spBFInc.size());
  }
#endif

//	GRUMMP::QualMeasure<3> *pQM = new GRUMMP::ShortestToRadius3D();
//	GRUMMP::SwapDecider3D *SwapDec3D = new GRUMMP::MetricSwapDecider3D(true,true,pQM);
//	GRUMMP::SwapManager3D *SwapMan3D = new GRUMMP::SwapManager3D(SwapDec3D, this);

  // Try swapping all the faces that were changed.
  // Helpful for 3D isotropic mesh quality; makes no difference for 3D
  // anisotropic removal success, which is a little surprising.
  // Don't swap if you were aiming for a particular edge, because you
  // might swap it away, after going to a lot of trouble to get it.
//	if (iNPref == 0) {
//		for (iterF = spFToCheck.begin(); iterF != spFToCheck.end(); ++iterF) {
//			Face *pF = *iterF;
//
//			iNewSwaps += iFaceSwap_deprecated(pF);
//			iNewSwaps += SwapMan3D->swapFace(pF);
//		}
//	}

//	delete SwapMan3D;
//	delete SwapDec3D;
//	delete pQM;

  // Success!
  sendEvents();
  iDelDone++;
  if (qBdryCurveVert) {
    iBdryCurveDelDone++;
  }
  else if (qBdryVert) {
    iBdryDelDone++;
  }
  return true;
}

