//@ Include files, #define's, etc.
#include <math.h>
#include "GR_assert.h"
#include "GR_misc.h"
#include "GR_Geometry.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_CellCV.h"
#include "GR_BFace.h"
#include "GR_Quality.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"
#include "GR_Vec.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
static eFaceCat eFC;
static int iSwapDepth = 0;
static const int iMaxSwapDepth = 75;

int aiEdgeReq[10], aiEdgeDone[10];

//@ Decide which type of triangular face this is, anyway.
eFaceCat
TriFace::categorizeFace(Vert*& pVVertA, Vert*& pVVertB, Vert*& pVVertC,
			Vert*& pVVertD, Vert*& pVVertE, TetCell* apTCTets[],
			int& iNTets, Vert*& pVPivot0, Vert*& pVPivot1,
			Vert*& pVOther, const double dMaxAngle) const
{
  //@@ Case: one or both cells is not a tet, including boundaries
  Cell* pCLeft = getLeftCell();
  Cell* pCRight = getRightCell();
  if (!pCLeft->isValid() || !pCRight->isValid())
    return (eOther);
  if (isBdryFace())
    return (eBdry);

  // The first clause is aimed primarily at faces that are locked for
  // swapping.  A primary application of this is in surface recovery.
  if (!isSwapAllowed() || pCLeft->getType() != Cell::eTet
      || pCRight->getType() != Cell::eTet)
    return (eOther);

  //@@ List all tets formed from five verts defining two cells sharing iFace
  iNTets = 0;
  apTCTets[0] = dynamic_cast<TetCell*>(pCLeft);
  apTCTets[1] = dynamic_cast<TetCell*>(pCRight);

  // The (Face*)this is really const, because nothing gets done here.
  // But the Vert*'s really should be non-const, because the data is
  // being returned to a caller who might, potentially, decide to swap
  // the mesh.
  pVVertA = const_cast<Vert*>(getVert(0));
  pVVertB = const_cast<Vert*>(getVert(1));
  pVVertC = const_cast<Vert*>(getVert(2));
  pVVertD = apTCTets[0]->getOppositeVert(this);
  pVVertE = apTCTets[1]->getOppositeVert(this);
  iNTets = 2;

  // Any add'l cell must be a neighbor of both of the first two.
  int i;
  for (i = 0; i < 4; ++i) {
    Face* pFCandFace = apTCTets[0]->getFace(i);
    Cell* pCCandTet = pFCandFace->getOppositeCell(apTCTets[0]);
    // Don't bother with cells that aren't tets; these are unswappable
    // anyway
    if (pCCandTet->getType() != Cell::eTet)
      continue;
    if (pCCandTet != apTCTets[1] && pCCandTet->hasVert(pVVertE)) {
      // If this face internal to the configuration is locked, then forget it.
      if (pFCandFace->isLocked())
	return (eOther);
      apTCTets[iNTets++] = dynamic_cast<TetCell*>(pCCandTet);
    }
  }

  assert(iNTets <= 4);

  switch (iNTets)
    {
    case 2:
      {
	//@@ Two tets
	assert(checkOrient3D(pVVertA, pVVertC, pVVertB, pVVertD) == 1);
	assert(checkOrient3D(pVVertA, pVVertB, pVVertC, pVVertE) == 1);
	//       if (checkOrient3D(pVVertA, pVVertC, pVVertB, pVVertD) != 1 ||
	// 	  checkOrient3D(pVVertA, pVVertB, pVVertC, pVVertE) != 1)
	// 	vFatalError("Orientation error", "FaceSwap");

	// The following code for eN32 face identification uses only as
	// many orientation primitives as required to determine the actual
	// face type.  This is very useful, as this is a critical code
	// path.

	// Initialize these to ludicrous values.
	int iOrientA = 100, iOrientB = 100, iOrientC = 100;
	bool qApproxFlatA = false, qApproxFlatB = false, qApproxFlatC = false;
	bool qBadOrientation = false;

	// There is always the possibility that the configuration has one
	// or two correctly-oriented but nearly flat (within the tolerance
	// set by VolMesh::dMaxAngleForSurfSwap) faces.  So if the
	// orientation is positive, check whether it should perhaps be
	// tagged as zero instead.

	iOrientA = checkOrient3D(pVVertB, pVVertC, pVVertD, pVVertE);
	if (iOrientA == -1) {
	  pVPivot0 = pVVertB;
	  pVPivot1 = pVVertC;
	  pVOther = pVVertA;
	  qBadOrientation = true;
	}
	else {
	  if (iOrientA == 1) {
	    // Check whether it's nearly flat
	    double adNorm0[3], adNorm1[3];
	    calcNormal3D(pVVertB->getCoords(), pVVertC->getCoords(),
			 pVVertD->getCoords(), adNorm0);
	    calcNormal3D(pVVertC->getCoords(), pVVertB->getCoords(),
			 pVVertE->getCoords(), adNorm1);
	    NORMALIZE3D(adNorm0);
	    NORMALIZE3D(adNorm1);
	    double dDot = dDOT3D(adNorm0, adNorm1);
	    if (dDot > cos(dMaxAngle * M_PI / 180)) {
	      qApproxFlatA = true;
	    }
	  }
	  iOrientB = checkOrient3D(pVVertC, pVVertA, pVVertD, pVVertE);
	  if (iOrientB == -1) {
	    pVPivot0 = pVVertC;
	    pVPivot1 = pVVertA;
	    pVOther = pVVertB;
	    qBadOrientation = true;
	  }
	  else if (iOrientB == 1) {
	    // Check whether it's nearly flat
	    double adNorm0[3], adNorm1[3];
	    calcNormal3D(pVVertC->getCoords(), pVVertA->getCoords(),
			 pVVertD->getCoords(), adNorm0);
	    calcNormal3D(pVVertA->getCoords(), pVVertC->getCoords(),
			 pVVertE->getCoords(), adNorm1);
	    NORMALIZE3D(adNorm0);
	    NORMALIZE3D(adNorm1);
	    double dDot = dDOT3D(adNorm0, adNorm1);
	    if (dDot > cos(dMaxAngle * M_PI / 180)) {
	      qApproxFlatB = true;
	    }
	  }
	  if (iOrientA + iOrientB == 0)
	    return eN20;
	  if (!qBadOrientation) {
	    iOrientC = checkOrient3D(pVVertA, pVVertB, pVVertD, pVVertE);
	    if (iOrientC == -1) {
	      pVPivot0 = pVVertA;
	      pVPivot1 = pVVertB;
	      pVOther = pVVertC;
	      qBadOrientation = true;
	    }
	    else if (iOrientC == 1) {
	      // Check whether it's nearly flat
	      double adNorm0[3], adNorm1[3];
	      calcNormal3D(pVVertA->getCoords(), pVVertB->getCoords(),
			   pVVertD->getCoords(), adNorm0);
	      calcNormal3D(pVVertB->getCoords(), pVVertA->getCoords(),
			   pVVertE->getCoords(), adNorm1);
	      NORMALIZE3D(adNorm0);
	      NORMALIZE3D(adNorm1);
	      double dDot = dDOT3D(adNorm0, adNorm1);
	      if (dDot > cos(dMaxAngle * M_PI / 180)) {
		qApproxFlatC = true;
	      }
	    }
	  }
	}
	if (qBadOrientation) {
	  // Need to distinguish between eN32 cases and boundary cases.
	  // The current code does this incompletely, in that one of the
	  // two cells must have a bdry face for the boundary case to be
	  // recognized.  Nevertheless, the main point here is to weed out
	  // the cases where two bdry tets have a slightly reflex edge
	  // between them, and correctly -not- identify these as eN32.
	  Face *pF0 = apTCTets[0]->getOppositeFace(pVOther);
	  Face *pF1 = apTCTets[1]->getOppositeFace(pVOther);

	  if (pF0->isBdryFace() || pF1->isBdryFace())
	    return eBdryReflex;
	  else
	    return eN32;
	}

	switch (iOrientA + iOrientB + iOrientC)
	  {
	  case 0: // Impossible, but included for completeness.
	    return eOther;
	  case 1: // Two orientations are zero (three points co-linear).
		  // Hopelessly unswappable.  Bail out.
	    return eN20;
	  case 2:
	    {
	      //@@@ Four points are coplanar; T22, T44, N44
	      // One orientation must be 0; verts are re-labeled to make it
	      // iOrientC.  This implies that AB is the coplanar edge.
	      assert(!(iOrientA && iOrientB && iOrientC));
	      if (iOrientA == 0) {
		Vert* pVTemp;
		pVTemp = pVVertA;
		pVVertA = pVVertB;
		pVVertB = pVVertC;
		pVVertC = pVTemp;
	      }
	      else if (iOrientB == 0) {
		Vert* pVTemp;
		pVTemp = pVVertA;
		pVVertA = pVVertC;
		pVVertC = pVVertB;
		pVVertB = pVTemp;
	      }
	      else if (iOrientC == 0) {
	      }
	      // The following assertion has been removed deliberately, to
	      // accomodate the near-flat cases tested for with dot products
	      // above.  Prior to adding support for near-flat swapping, the
	      // orientation was in fact re-aligned correctly, so there's no
	      // cause for concern.
	      // Verify that the relabeling was done correctly
	      //	assert( checkOrient3D(pVVertA, pVVertB, pVVertD, pVVertE) ==
	      //	0 );

	      // The configuration of these two tets is classified based on the
	      // properties of the coplanar faces.
	      // 1 If both are BFaces with the same boundary condition,
	      //   swappable two tets to two.
	      // 2 If both are BFaces with different bdry cond, not swappable.
	      // 3 If one is a BFace and the other not, not swappable.
	      // 4 If neither is a BFace, both opposite cells are tets, and the
	      //   tets have the same fourth vert, swappable four to four.
	      // 5 If neither is a BFace, both opposite cells are tets, and the
	      //   tets do not have the same fourth vert, not swappable,
	      //   although some non-local transformations might make it so.
	      // 6 If neither is a BFace and one or both opposite cells is not a
	      //   tet, not swappable.
	      TetCell *pTCTetA = apTCTets[0], *pTCTetB = apTCTets[1];
	      // These are the faces that are coplanar
	      Face* pFFaceA = pTCTetA->getOppositeFace(pVVertC);
	      Face* pFFaceB = pTCTetB->getOppositeFace(pVVertC);
	      Cell* pCOppCellA = pFFaceA->getOppositeCell(pTCTetA);
	      Cell* pCOppCellB = pFFaceB->getOppositeCell(pTCTetB);

	      if (pCOppCellA->getType() == Cell::eTriBFace
		  && pCOppCellB->getType() == Cell::eTriBFace) {
		// Both faces are on the boundary
		TriBFace* pTBFA = dynamic_cast<TriBFace*>(pCOppCellA);
		TriBFace* pTBFB = dynamic_cast<TriBFace*>(pCOppCellB);
		// These are already -known- to be co-planar, so checking bdry
		// patches is sufficient.
		// FIX ME: 0.3  This is a spot where checking patches, not
		// BC's is the right thing to do.  But until patch
		// optimization gets written, the BC check is the better
		// choice, because then you can swap coplanar faces in meshes
		// that you read in.
		// 	  if (pTBFA->pPatchPointer() == pTBFB->pPatchPointer())
		if (pTBFA->getBdryCondition() == pTBFB->getBdryCondition())
		  return eT22; // case 1
		else
		  return eOther;  // case 2
	      }
	      else if (pCOppCellA->getType() == Cell::eTet
		  && pCOppCellB->getType() == Cell::eTet) {
		// Bail if any of the internal faces are locked.
		if (!pFFaceA->isSwapAllowed() || !pFFaceB->isSwapAllowed())
		  return (eOther);
		// Both faces are internal to the mesh and separate pairs of tets
		TetCell* pTCA = dynamic_cast<TetCell*>(pCOppCellA);
		TetCell* pTCB = dynamic_cast<TetCell*>(pCOppCellB);
		assert(pTCA && pTCB);
		if (pTCA->getOppositeVert(pFFaceA)
		    == pTCB->getOppositeVert(pFFaceB)) {
		  Face *pFCommon = findCommonFace(pTCA, pTCB);
		  assert(pFCommon); // Has to exist or we wouldn't be here.
		  if (pFCommon->isLocked()) {
		    return eOther;
		  }
		  else {
		    return eT44; // case 4
		  }
		}
		else {
		  pVPivot0 = pVVertA;
		  pVPivot1 = pVVertB;
		  pVOther = pVVertC;
		  return eN44; // case 5
		}
	      }
	      else
		// Exactly one face on the boundary or internal faces with
		// cells other than tets
		return eOther; // cases 3,6
	    } // End two-tet cases with four verts coplanar
	  case 3:
	    {
	      // This is probably a T23 case, but it might be a perturbed T22
	      // case, which should be tagged so that removal of surface verts
	      // will work better.  (Actually, it could even be a perturbed
	      // N20, which is a disaster during removal if you try to run
	      // with it as a T23! so nuke this case.)
	      if ((qApproxFlatA && qApproxFlatB)
		  || (qApproxFlatB && qApproxFlatC)
		  || (qApproxFlatC && qApproxFlatA))
		return eN20;
	      if (qApproxFlatA || qApproxFlatB || qApproxFlatC) {
		// It is a perturbed T22
		if (qApproxFlatA) {
		  Vert* pVTemp;
		  pVTemp = pVVertA;
		  pVVertA = pVVertB;
		  pVVertB = pVVertC;
		  pVVertC = pVTemp;
		}
		else if (qApproxFlatB) {
		  Vert* pVTemp;
		  pVTemp = pVVertA;
		  pVVertA = pVVertC;
		  pVVertC = pVVertB;
		  pVVertB = pVTemp;
		}
		else if (qApproxFlatC) {
		}

		TetCell *pTCTetA = apTCTets[0], *pTCTetB = apTCTets[1];
		// These are the faces that are coplanar
		Face* pFFaceA = pTCTetA->getOppositeFace(pVVertC);
		Face* pFFaceB = pTCTetB->getOppositeFace(pVVertC);
		Cell* pCOppCellA = pFFaceA->getOppositeCell(pTCTetA);
		Cell* pCOppCellB = pFFaceB->getOppositeCell(pTCTetB);

		if (pCOppCellA->getType() == Cell::eTriBFace
		    && pCOppCellB->getType() == Cell::eTriBFace) {
		  // Both faces are on the boundary
		  TriBFace* pTBFA = dynamic_cast<TriBFace*>(pCOppCellA);
		  TriBFace* pTBFB = dynamic_cast<TriBFace*>(pCOppCellB);
		  // These are already -known- to be co-planar, so checking bdry
		  // patches is sufficient.
		  // FIX ME: 0.3  This is a spot where checking patches, not
		  // BC's is the right thing to do.  But until patch
		  // optimization gets written, the BC check is the better
		  // choice, because then you can swap coplanar faces in meshes
		  // that you read in.
		  // 	  if (pTBFA->pPatchPointer() == pTBFB->pPatchPointer())
		  if (pTBFA->getBdryCondition() == pTBFB->getBdryCondition())
		    return eT22; // case 1
		  else
		    return eOther;  // case 2
		}
	      } // If one is nearly flat
		// Otherwise, the configuration is convex; swappable two for three.
	      return eT23;
	    }
	  default:
	    // No other cases should be possible
	    assert2(0, "Impossible convexity case");
	    return eOther;
	  }
	break;
      } // End of two tet cases
    case 3:
      {
	//@@ Three tet cases
	int iOrient1, iOrient2;
	pVVertA = pVVertE;
	pVVertB = pVVertD;
	pVVertC = const_cast<Vert*>(getVert(0));
	pVVertD = const_cast<Vert*>(getVert(1));
	pVVertE = const_cast<Vert*>(getVert(2));
	if (!apTCTets[2]->hasVert(pVVertC)) {
	}
	else if (!apTCTets[2]->hasVert(pVVertE)) {
	  Vert* pVTemp = pVVertC;
	  pVVertC = pVVertE;
	  pVVertE = pVVertD;
	  pVVertD = pVTemp;
	}
	else {
	  Vert* pVTemp = pVVertC;
	  pVVertC = pVVertD;
	  pVVertD = pVVertE;
	  pVVertE = pVTemp;
	}

	// The first of these returns will save a few orientation
	// primitive evaluations
	iOrient1 = checkOrient3D(pVVertA, pVVertC, pVVertB, pVVertD);
	if (iOrient1 == -1)
	  return eOther;
	else if (iOrient1 == 0)
	  return eN30;
	iOrient2 = checkOrient3D(pVVertA, pVVertB, pVVertC, pVVertE);
	if (iOrient2 == -1)
	  return eOther;
	else if (iOrient2 == 0)
	  return eN30;
	assert(iOrient1 == 1);
	assert(iOrient2 == 1);
	//@@@@ T32
	assert(!apTCTets[0]->hasVert(pVVertA));
	assert(!apTCTets[1]->hasVert(pVVertB));
	assert(!apTCTets[2]->hasVert(pVVertC));
	Face *pF01 = findCommonFace(apTCTets[0], apTCTets[1]);
	Face *pF12 = findCommonFace(apTCTets[1], apTCTets[2]);
	Face *pF20 = findCommonFace(apTCTets[2], apTCTets[0]);
	assert(
	    pF01 != pFInvalidFace && pF12 != pFInvalidFace && pF20 != pFInvalidFace);
	if (pF01->isLocked() || pF12->isLocked() || pF20->isLocked())
	  return eOther;
	else
	  return eT32;
      } // End of three tet cases
    case 4:
      //@@@ Four tets => N40
      return eN40;
      // End of four tet cases
    } // End of classification/weeding
      // Should never reach here, but just in case...
  return eOther;
}

//@ Quality function for purposes of edge swapping. 
// Would also work for face swapping, but there are more efficient
// approaches for that case.  dQualFunc must return a value which is
// larger for more desirable sets of four verts.
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

double
VolMesh::evalQualFunc(Vert* const pV0, Vert* const pV1, Vert* const pV2,
		      Vert* const pV3) const
{
  Vert *apV[] =
    { pV0, pV1, pV2, pV3 };
  TetCellCV FT(apV);
  switch (m_swapMeasure)
    {
    case eDelaunay:  // Seek to maximize shortest edge to circumradius.
      {
	double adRes[1];
	int i;
	calcShortestToRadius(&FT, &i, adRes);
	return adRes[0];
      }
    case eMaxDihed:  // Seek to minimize max dihedral angle.  Return
      // negative of max dihedral.
      return -FT.calcMaxDihed();
    case eMinSine:   // Seek to maximize min sine of dihedral angle.
      // Return min difference of dihedral from 0 or 180
      // degrees.
      {
	double adDihed[6];
	int iDum;
	FT.calcAllDihed(adDihed, &iDum);
	assert(iDum == 6);
	double dMinDiff = 90;
	for (iDum = 0; iDum < 6; iDum++) {
	  double dDiff = 90 - fabs(90 - adDihed[iDum]);
	  if (dDiff < dMinDiff)
	    dMinDiff = dDiff;
	}
	return dMinDiff;
      }
    case eFromQual:
      {
	double dValue = m_qual->dEvaluate(&FT);
	bool qMaximizeMin = m_qual->qMaximizeMinValue();
	return qMaximizeMin ? dValue : -dValue;
      }
    default:
      assert(0);
      break;
    }
  return (-1000);
}

// Pull in the >1200 lines of code that initializes data for the final
// configurations after edge swapping.  
#include "InitCanon.cxx"

//@ Locally reconfigure tets incident on an edge to improve mesh quality
// The edge given by pVNorth, pVSouth is a candidate for removal.  If
// the edge doesn't lie on a boundary and a reconfiguration which
// improves mesh quality exists, then a swap is done, replacing 4 tets
// with 4, 5 with 6, 6 with 8, or 7 with 10.  Generally N tets are
// replaced by 2(N-2).  At this point, the largest transformation
// considered is 7 tets to 10.  The number of possible configurations
// for 8 to 12 swapping is somewhere around 150, so it's not
// considered, at least for now.
//
// After the edge swap is complete, all new faces are face swapped.
// The return value is the total number of face swaps done, with edge
// swaps counting as N-2 face swaps.

int
VolMesh::edgeSwap3D(Face* const pF, Vert* pVNorth, Vert* pVSouth,
		    Vert* const pVOther, const bool qAllowAnyValid,
		    const Vert* const pVWant0, const Vert* const pVWant1)
{

  // The last three args are optional and used only for forced edge
  // swapping during surface recovery.  qAllowAnyValid (default:0) is
  // set to 1 to indicate that any submesh which removes the specified
  // edge is preferred.  This implies that the initial quality should be
  // set to some artifically low value and initial quality computation
  // skipped.  If there's some edge that we want in the mesh,
  // the verts of that edge are given as pVWant0 and pVWant1 (default:
  // pVInvalidVert).  A quality bonus is artificially given to
  // configurations which have the magic edge.

  // TODO: Go directly to a swap decider
  GRUMMP::MinMaxDihedSwapDecider3D SwapDec3D(true);
  GRUMMP::SwapManager3D SwapMan3D(&SwapDec3D, this);

  assert(pF->isValid());
  assert(pVNorth->isValid());
  assert(pVSouth->isValid());
  assert(pVOther->isValid());
  assert(pF->hasVert(pVNorth));
  assert(pF->hasVert(pVSouth));
  assert(pF->hasVert(pVOther));
  // This call only needs to appear once, but that's okay, because
  // subsequent calls short-circuit and exit immediately.
  vInitCanonConfigs();

  const double dInvalidPenalty = -10000;
  const double dBonus = 100;

  Vert *apVOuter[10];
  Cell *apCOrig[10];
  Face *apFNorth[10];
  Face *apFSouth[10];
  Face *apFOrig[10];
  // Determine the size of the orbit around the pivot edge and list
  // all the verts, faces, and cells involved.
  apFOrig[0] = pF;
  apVOuter[0] = pVOther;
  apCOrig[0] = pF->getLeftCell();

  if (apCOrig[0]->getType() == Cell::eTriBFace) {

    return bdryEdgeSwap3D(pF, pVNorth, pVSouth, pVOther);
  }
  if (apCOrig[0]->getType() == Cell::eIntTriBFace) {
    if (areSubsegVertsConnected(pVNorth, pVSouth))
      return false;

    Face * pFOther = (
	apCOrig[0]->getFace(0) == pF ?
	    apCOrig[0]->getFace(1) : apCOrig[0]->getFace(0));
    assert(pFOther != pF);
    GRUMMP::FaceSwapInfo3D FC1(pF, pVNorth, pVSouth, pVOther);
    GRUMMP::FaceSwapInfo3D FC2(pFOther, pVNorth, pVSouth, pVOther);

    GRUMMP::EdgeSwapInfo ES1(FC1);
    GRUMMP::EdgeSwapInfo ES2(FC2);
    // going to do a check to see they agree on boundary faces;
    BFace * bface11 = ES1.getBFace(0);
    BFace * bface12 = ES1.getBFace(1);
    BFace * bface21 = ES2.getBFace(0);
    BFace * bface22 = ES2.getBFace(1);
    if (!((bface11 == bface21 || bface11 == bface22)
	&& (bface12 == bface21 || bface12 == bface22)))
      return false;
    //  check if both sides can do the swap, if they both can, do the swap, with a forced removal
    bool doSwap1 = SwapDec3D.doEdgeSwap(ES1, const_cast<Vert*>(pVWant0),
					const_cast<Vert*>(pVWant1), false);
    bool doSwap2 = SwapDec3D.doEdgeSwap(ES2, const_cast<Vert*>(pVWant0),
					const_cast<Vert*>(pVWant1), false);

    if (doSwap1 && doSwap2) {
      SwapMan3D.reconfigureIntBdryEdge(ES1, ES2);
      return ES1.getNumOrigVerts() - 2 + ES2.getNumOrigVerts() - 2;
    }
    else
      return 0;
  }

  if (apCOrig[0]->getType() != Cell::eTet)
    return 0;

  int iReg = apCOrig[0]->getRegion();
  apFNorth[0] = apCOrig[0]->getOppositeFace(pVSouth);
  apFSouth[0] = apCOrig[0]->getOppositeFace(pVNorth);
  Face *pFNew = apCOrig[0]->getOppositeFace(apVOuter[0]);
  assert(pFNew != pF);
  int i = 1;
  while (i < 10 && pF != pFNew) {

    apFOrig[i] = pFNew;
    apVOuter[i] = apCOrig[i - 1]->getOppositeVert(apFOrig[i - 1]);
    apCOrig[i] = pFNew->getOppositeCell(apCOrig[i - 1]);
    // Hit a boundary or non-tet cell.
    if (apCOrig[i]->getType() == Cell::eTriBFace) {
      return bdryEdgeSwap3D(pFNew, pVNorth, pVSouth, apVOuter[i]);
    }
    if (apCOrig[i]->getType() == Cell::eIntTriBFace) {
      if (areSubsegVertsConnected(pVNorth, pVSouth))
	return false;

      Face * pFOther = (
	  apCOrig[i]->getFace(0) == pFNew ?
	      apCOrig[i]->getFace(1) : apCOrig[i]->getFace(0));
      assert(pFOther != pFNew);
      GRUMMP::FaceSwapInfo3D FC1(pFNew, pVNorth, pVSouth, apVOuter[i]);
      GRUMMP::FaceSwapInfo3D FC2(pFOther, pVNorth, pVSouth, apVOuter[i]);

      GRUMMP::EdgeSwapInfo ES1(FC1);
      GRUMMP::EdgeSwapInfo ES2(FC2);
      // going to do a check to see they agree on boundary faces;

      BFace * bface11 = ES1.getBFace(0);
      BFace * bface12 = ES1.getBFace(1);
      BFace * bface21 = ES2.getBFace(0);
      BFace * bface22 = ES2.getBFace(1);
      if (!((bface11 == bface21 || bface11 == bface22)
	  && (bface12 == bface21 || bface12 == bface22)))
	return false;

      bool doSwap1 = SwapDec3D.doEdgeSwap(ES1, const_cast<Vert*>(pVWant0),
					  const_cast<Vert*>(pVWant1), false);
      bool doSwap2 = SwapDec3D.doEdgeSwap(ES2, const_cast<Vert*>(pVWant0),
					  const_cast<Vert*>(pVWant1), false);
      if (doSwap1 && doSwap2) {
	SwapMan3D.reconfigureIntBdryEdge(ES1, ES2);

	return ES1.getNumOrigVerts() - 2 + ES2.getNumOrigVerts() - 2;
      }
      else
	return 0;
    }
    if (apCOrig[i]->getType() != Cell::eTet)
      return 0;
    // Give up if more than one region is incident on the edge.
    if (iReg != apCOrig[i]->getRegion())
      return 0;
    apFNorth[i] = apCOrig[i]->getOppositeFace(pVSouth);
    apFSouth[i] = apCOrig[i]->getOppositeFace(pVNorth);
    pFNew = apCOrig[i]->getOppositeFace(apVOuter[i]);
    i++;
  }
//#ifndef NDEBUG
//					std::set<Cell*> spCInc;
//					std::set<Vert*> spVTemp;
//	findNeighborhoodInfo(pVNorth,spCInc,spVTemp);
//	int iAdj = 0;
//	for(std::set<Cell*>::iterator it = spCInc.begin(); it != spCInc.end(); ++it)
//		iAdj += ((*it)->hasVert(pVSouth));
//	assert(iAdj == i);
//#endif
  if (pFNew != pF)
    return 0; // Ran out of space

  int iNOrig = std::min(i, 9);
  aiEdgeReq[iNOrig]++;
  if (iNOrig > 7) {
    logMessage(4, "Too many tets around edge: %d\n", iNOrig);
    return 0;
  }
  double dOrigQual;
  if (qAllowAnyValid)
    dOrigQual = dInvalidPenalty + 1;
  else {
    // Compute the quality of the current configuration.  This is the
    // smallest quality for any of the existing cells.
    dOrigQual = 1000;
    for (i = 0; i < iNOrig; i++) {
      Cell *pC = apCOrig[i];
      // Evaluate quality, depending on the swap criterion in use.
      double dQual = evalQualFunc(pC->getVert(0), pC->getVert(1),
				  pC->getVert(2), pC->getVert(3));
      dOrigQual = min(dQual, dOrigQual);
    }
  }

  // Set up info about the possible new configurations.

  // Compute quality for each possible tet pair (one pair per possible
  // triangle).
  double a3dNewQual[10][10][10];
  int iSign = checkOrient3D(apVOuter[0], apVOuter[1], pVSouth, pVNorth);
  int iV0, iV1, iV2;
#ifndef NDEBUG
  for (iV0 = 0; iV0 < 10; iV0++)
    for (iV1 = 0; iV1 < 10; iV1++)
      for (iV2 = 0; iV2 < 10; iV2++)
	a3dNewQual[iV0][iV1][iV2] = 0;
  for (int iDum = 0; iDum < iNOrig; iDum++)
    assert(
	iSign
	    == checkOrient3D(apVOuter[iDum], apVOuter[(iDum + 1) % iNOrig],
			     pVSouth, pVNorth));
#endif
  for (iV0 = 0; iV0 < iNOrig; iV0++)
    for (iV1 = iV0 + 1; iV1 < iNOrig; iV1++)
      for (iV2 = iV1 + 1; iV2 < iNOrig; iV2++) {
	if (checkOrient3D(apVOuter[iV0], apVOuter[iV1], apVOuter[iV2], pVNorth)
	    == iSign)
	  // Evaluate quality
	  a3dNewQual[iV0][iV1][iV2] = evalQualFunc(apVOuter[iV0], apVOuter[iV1],
						   apVOuter[iV2], pVNorth);
	else
	  a3dNewQual[iV0][iV1][iV2] = dInvalidPenalty;

	// Don't bother evaluating for the other tet if this one is
	// worse that the worst in the original config, because this tet
	// pair can never be in a final config anyway.
	if (a3dNewQual[iV0][iV1][iV2] >= dOrigQual) {
	  double dDummyQual;
	  if (checkOrient3D(apVOuter[iV0], apVOuter[iV2], apVOuter[iV1],
			    pVSouth) == iSign)
	    dDummyQual = evalQualFunc(apVOuter[iV0], apVOuter[iV2],
				      apVOuter[iV1], pVSouth);
	  else
	    dDummyQual = dInvalidPenalty;
	  a3dNewQual[iV0][iV1][iV2] = min(a3dNewQual[iV0][iV1][iV2],
					  dDummyQual);
	}

	// Copy things around symmetrically
	a3dNewQual[iV2][iV1][iV0] = a3dNewQual[iV2][iV0][iV1] =
	    a3dNewQual[iV1][iV2][iV0] = a3dNewQual[iV1][iV0][iV2] =
		a3dNewQual[iV0][iV2][iV1] = a3dNewQual[iV0][iV1][iV2];

      }

  // Find the new configuration with the highest quality.

  double dMaxQual = dOrigQual;
  int iCanon, iPerm;
  int iBestCanon = -1, iBestPerm = -1;
  bool qBonusPossible = (pVWant0 != pVInvalidVert);

  for (iCanon = 0; iCanon < aCS[iNOrig].iNCanon; iCanon++)
    for (iPerm = 0; iPerm < aCS[iNOrig].aConf[iCanon].iNPerm; iPerm++) {
      int iBonusGiven = 0;
      double dConfQual = 1000;
      // Once this config is known to be worse than the best one we know
      // about, we can give up on it.
      for (int iTri = 0;
	  iTri < iNOrig - 2 && (dConfQual >= dMaxQual || qBonusPossible);
	  iTri++) {
	iV0 = (aCS[iNOrig].aConf[iCanon].a2iFaceVert[iTri][0] + iPerm) % iNOrig;
	iV1 = (aCS[iNOrig].aConf[iCanon].a2iFaceVert[iTri][1] + iPerm) % iNOrig;
	iV2 = (aCS[iNOrig].aConf[iCanon].a2iFaceVert[iTri][2] + iPerm) % iNOrig;
	// Give a bonus if a particular edge is desired and that edge
	// is included in the configuration.
	if (qBonusPossible && !(iBonusGiven == 2)) {
	  int iThisBonus = 0;
	  Vert * pV0 = apVOuter[iV0];
	  Vert * pV1 = apVOuter[iV1];
	  Vert * pV2 = apVOuter[iV2];
	  if (pVWant0 == pV0 || pVWant0 == pV1 || pVWant0 == pV2)
	    iThisBonus++;
	  if (pVWant1 == pV0 || pVWant1 == pV1 || pVWant1 == pV2)
	    iThisBonus++;
	  if (iThisBonus > iBonusGiven)
	    iBonusGiven = iThisBonus;
	}
	dConfQual = min(dConfQual, a3dNewQual[iV0][iV1][iV2]);

	// Only give a bonus to valid configurations
	if (qBonusPossible && dConfQual > dOrigQual)
	  dConfQual += iBonusGiven * dBonus;

      }

      if (dConfQual > dMaxQual) {
	dMaxQual = dConfQual;
	iBestCanon = iCanon;
	iBestPerm = iPerm;
      }
    }

  // If the highest-quality alternative is worse than the current
  // configuration, do nothing.  The check for iBestCanon is there so
  // that some stupid roundoff in comparison of dMaxQual and dOrigQual
  // won't cause a reconfiguration when there shouldn't be one.
  if (dMaxQual <= dOrigQual || iBestCanon == -1) {
    return 0;
  }
  else {
    // If the highest-quality alternative is better than the current
    // configuration:
    aiEdgeDone[iNOrig]++;

    // Delete all of the interior tets and faces.  These must be
    // separate loops, or forced.
    for (i = 0; i < iNOrig; i++) {
      deleteCell(apCOrig[i]);
    }

    // Set up an array to hold all the new cells
    Cell** apCNew = new Cell*[2 * iNOrig - 4];

    // Renumber all the old face, vertex and cell data, offsetting it
    // by iBestPerm to make the desired configuration match the stored
    // canonical configuration.
    if (iBestPerm != 0) {
      Cell **apCTemp = new Cell*[iNOrig];

      Face **apFNorthTemp = new Face*[iNOrig];
      Face **apFSouthTemp = new Face*[iNOrig];
      Vert **apVOuterTemp = new Vert*[iNOrig];
      // Copy data to temporary arrays.
      for (i = 0; i < iNOrig; i++) {
	apCTemp[i] = apCOrig[i];
	apFNorthTemp[i] = apFNorth[i];
	apFSouthTemp[i] = apFSouth[i];
	apVOuterTemp[i] = apVOuter[i];
      }
      // Copy back with an offset.
      for (i = 0; i < iNOrig; i++) {
	apCOrig[i] = apCTemp[(i + iBestPerm) % iNOrig];
	apFNorth[i] = apFNorthTemp[(i + iBestPerm) % iNOrig];
	apFSouth[i] = apFSouthTemp[(i + iBestPerm) % iNOrig];
	apVOuter[i] = apVOuterTemp[(i + iBestPerm) % iNOrig];
      }
      delete[] apCTemp;
      delete[] apFNorthTemp;
      delete[] apFSouthTemp;
      delete[] apVOuterTemp;
    }

    // Re-connect the mesh locally
    Config *pConf = &(aCS[iNOrig].aConf[iBestCanon]);
    // Assign data to new faces
    int iOrient = 0;

    // There are iNOrig - 2 triangles on the equator, two tets apiece,
    // four verts per tet.
    for (i = 0; i < iNOrig - 2; i++) {
      // These are the equatorial tris, so all these indices should be
      // positive.
      assert(pConf->a2iFaceVert[i][0] >= 0);
      assert(pConf->a2iFaceVert[i][1] >= 0);
      assert(pConf->a2iFaceVert[i][2] >= 0);
      Vert* pVert0 = apVOuter[pConf->a2iFaceVert[i][0]];
      Vert* pVert1 = apVOuter[pConf->a2iFaceVert[i][1]];
      Vert* pVert2 = apVOuter[pConf->a2iFaceVert[i][2]];

      if (iOrient == 0)
	iOrient = checkOrient3D(pVert0, pVert1, pVert2, pVNorth);

      if (iOrient == 1) {
	assert(checkOrient3D(pVert0, pVert1, pVert2, pVNorth) == 1);
	assert(checkOrient3D(pVert0, pVert2, pVert1, pVSouth) == 1);
	bool qExist;
	apCNew[2 * i] = createTetCell(qExist, pVert0, pVert1, pVert2, pVNorth,
				      iReg);
	assert(!qExist);
	apCNew[2 * i + 1] = createTetCell(qExist, pVert0, pVert2, pVert1,
					  pVSouth, iReg);
	assert(!qExist);
      }
      else {
	assert(checkOrient3D(pVert0, pVert2, pVert1, pVNorth) == 1);
	assert(checkOrient3D(pVert0, pVert1, pVert2, pVSouth) == 1);
	bool qExist;
	apCNew[2 * i] = createTetCell(qExist, pVert0, pVert2, pVert1, pVNorth,
				      iReg);
	assert(!qExist);
	apCNew[2 * i + 1] = createTetCell(qExist, pVert0, pVert1, pVert2,
					  pVSouth, iReg);
	assert(!qExist);
      }
    }

#ifndef OMIT_VERTEX_HINTS
    // Reset vertex hints for polar verts
    pVNorth->setHintFace(apFNorth[0]);
    pVSouth->setHintFace(apFSouth[0]);
#endif

    // Validate the result
    for (i = 0; i < iNOrig * 2 - 4; i++) {
      // Check everything about the tet
      assert(apCNew[i]->doFullCheck());
      // Every tet has to hit one pole or the other
      assert(apCNew[i]->hasVert(pVNorth) || apCNew[i]->hasVert(pVSouth));
      // Assert tetrahedron orientation
      assert(
	  checkOrient3D(apCNew[i]->getVert(0), apCNew[i]->getVert(1),
			apCNew[i]->getVert(2), apCNew[i]->getVert(3)) == 1);
    }

    for (i = 0; i < iNOrig; i++) {
      assert(apFSouth[i]->doFullCheck());
      assert(apFNorth[i]->doFullCheck());
    }

    // If recursion is permitted, check all the new faces for possible
    // swaps, as well as all of the faces on the outside of the orbit.
    int iRetVal = iNOrig - 2;
    iSwapDepth++;
    if (iSwapDepth > iMaxSwapDepth)
      logMessage(3, "Swap recursion depth exceeded!\n");

    Face** apFRecur = new Face*[(2 * iNOrig - 4) * 4];
    // Grab all the faces for all the tets that were just created...
    Face **ppFCurr = apFRecur;
    for (i = 0; i < 2 * iNOrig - 4; i++) {
      Cell *pC = apCNew[i];
      for (int j = 0; j < 4; j++) {
	(*ppFCurr) = pC->getFace(j);
	ppFCurr++;
      }
    }
    assert(ppFCurr - apFRecur == (2 * iNOrig - 4) * 4);

    // Now sort them uniquely.
    std::sort(apFRecur, ppFCurr);
    Face **ppFNewEnd = std::unique(apFRecur, ppFCurr);
    assert(ppFNewEnd - apFRecur == 5 * iNOrig - 8);

    // And now swap them all.
    for (i = 0; i < 5 * iNOrig - 8; i++) {
      if (m_swapRecursively && iSwapDepth <= iMaxSwapDepth) {
	iRetVal += iFaceSwap_deprecated(apFRecur[i]);
      }
    }
    delete[] apFRecur;
    delete[] apCNew;

    // Return the total (equivalent) number of swaps performed
    iSwapDepth--;
    return iRetVal;
  }
}

//@ Locally reconfigure tets incident on a bdry edge to improve mesh quality
// The bdry edge given by pVNorth, pVSouth is a candidate for removal.
// pF is required to be a boundary face, and pVOther is the third vertex
// on that face.
// 
// If the bdry edge has two incident bdry faces that are coplanar and
// have the same boundary condition, and if a reconfiguration which
// improves mesh quality exists, then a swap is done, replacing 2 tets
// with 2, 3 with 4, ... or 6 with 10.  Generally N tets are replaced by
// 2(N-1).  At this point, the largest transformation considered is 6
// tets to 10.  The number of possible configurations for 7 to 12
// swapping is somewhere around 150, so it's not considered, at least
// for now. 
//
// All of these reconfigurations are close cousins to interior edge
// swaps.  A 2-to-2 boundary swap is just a 3-to-2 swap with one tet
// flattened into oblivion.  This means that all of the canonical
// configuration info from interior can be re-used here.
//
// After the edge swap is complete, all new interior faces are face
// swapped. The return value is the total number of face swaps done,
// with boundary edge swaps counting as N-3 face swaps. 
int
VolMesh::bdryEdgeSwap3D(Face* const pF, Vert* pVNorth, Vert* pVSouth,
			Vert* const pVOther)
{
  // Note the optional args used for forced edge swapping during surface
  // recovery aren't used here because they aren't needed.
  //
  if (areSubsegVertsConnected(pVNorth, pVSouth))
    return 0;

  assert(isSimplicial());
  assert(pF->getType() == Face::eTriFace);
  assert(pF->isValid());
  assert(pVNorth->isValid());
  assert(pVSouth->isValid());
  assert(pVOther->isValid());
  assert(pF->hasVert(pVNorth));
  assert(pF->hasVert(pVSouth));
  assert(pF->hasVert(pVOther));
  // The face must be on the boundary.
  assert(
      pF->getLeftCell()->getType() == Cell::eTriBFace
	  || pF->getRightCell()->getType() == Cell::eTriBFace);

  // This call only needs to appear once, but that's okay, because
  // subsequent calls short-circuit and exit immediately.
  vInitCanonConfigs();
  if (!m_OKToSwapSurfaceEdges // Surface edge swapping turned off
  || !m_allowBdryChanges)   // All surface changes turned off
    return 0;

  const double dInvalidPenalty = -10000;

  //@@ Identify all cells and faces incident on the bad edge
  Vert *apVOuter[11];
  Cell *apCOrig[11];
  Face *apFNorth[11];
  Face *apFSouth[11];
  Face *apFOrig[11];
  BFace *pBFBFace0, *pBFBFace1;

  // Determine the size of the orbit around the pivot edge and list
  // all the verts, faces, and cells involved.
  apFOrig[0] = pF;
  apVOuter[0] = pVOther;
  // ID the boundary associated with the original face *pF
  {
    Cell *pCTemp = pF->getLeftCell();
    if (pCTemp->getType() == Cell::eTriBFace) {
      pBFBFace0 = dynamic_cast<BFace*>(pCTemp);
      apCOrig[0] = pF->getRightCell();
    }
    else {
      pBFBFace0 = dynamic_cast<BFace*>(pF->getRightCell());
      apCOrig[0] = pCTemp;
    }
  }
  apFNorth[0] = apCOrig[0]->getOppositeFace(pVSouth);
  apFSouth[0] = apCOrig[0]->getOppositeFace(pVNorth);
  Face *pFNew = apCOrig[0]->getOppositeFace(apVOuter[0]);
  int i = 1;
  int iReg = apCOrig[0]->getRegion();
  while (i <= 10) {
    // Abort if there is more than one region involved.
    if (pFNew->getFaceLoc() == Face::eBdryTwoSide)
      return 0;
    apFOrig[i] = pFNew;
    apVOuter[i] = apCOrig[i - 1]->getOppositeVert(apFOrig[i - 1]);
    apCOrig[i] = pFNew->getOppositeCell(apCOrig[i - 1]);
    // Hit a boundary or non-tet cell.
    if (apCOrig[i]->getType() == Cell::eTriBFace)
      break;
    assert(apCOrig[i]->getType() == Cell::eTet);
    // Give up now if there is more than one region incident on the edge.
    if (apCOrig[i]->getRegion() != iReg)
      return 0;
    apFNorth[i] = apCOrig[i]->getOppositeFace(pVSouth);
    apFSouth[i] = apCOrig[i]->getOppositeFace(pVNorth);
    pFNew = apCOrig[i]->getOppositeFace(apVOuter[i]);
    i++;
  }
  int iNOrigCells = std::min(9, i);
  int iNOrigVerts = iNOrigCells + 1;
  aiEdgeReq[iNOrigCells]++; // Record the attempt even if we're bailing.

  // Do we have more cells incident on the edge than we're prepared to
  // handle?
  if (iNOrigCells > 6) {
    logMessage(4, "Too many tets around bdry edge: %d\n", iNOrigCells);
    return 0;
  }

  // Assign data for the second boundary face
  pBFBFace1 = dynamic_cast<BFace*>(apCOrig[i]);

  // A couple more easy exits for cases where surface edge swapping is not
  // legal.

  // Did we run out of space while identifying the entities around the
  // edge we'd like to swap?
  if (apCOrig[i]->getType() != Cell::eTriBFace)
    return 0;

  // Do the boundary faces have the same BC?
  if (!areBdryFacesOKToSwap(dynamic_cast<TriBFaceBase*>(pBFBFace0),
			    dynamic_cast<TriBFaceBase*>(pBFBFace1)))
    return 0;

  // Are the boundary faces near enough to being coplanar?

  // Strict coplanarity test for surface edge swapping.  Tested by
  // checking whether the four surface points in the configuration form
  // a degenerate tet.

  //   if (checkOrient3D(pVNorth, pVSouth, pVOther, apVOuter[iNOrigVerts-1]) != 0)
  //     return 0;

  // Allow some slightly non-planar surface edge swapping.
  {
    double adNorm1[3], adNorm2[3];
    calcUnitNormal(pVNorth->getCoords(), pVSouth->getCoords(),
		   pVOther->getCoords(), adNorm1);
    calcUnitNormal(pVNorth->getCoords(), apVOuter[iNOrigVerts - 1]->getCoords(),
		   pVSouth->getCoords(), adNorm2);
    double dCos = dDOT3D(adNorm1, adNorm2);
    double dMinCos = cos(m_maxAngleForSurfSwap * M_PI / 180.);
    if (iFuzzyComp(dCos, dMinCos) == -1)
      return 0;
  }

  // Do the two boundary triangles form a convex quad?  If so, the
  // orientation of the tets formed by these triangles and any
  // non-coplanar point (i.e., any other member of apVOuter) will be the
  // same.  This exit is redundant (no swapping would occur for this
  // case anyway), but saves a lot of quality computations.
  if (checkOrient3D(pVOther, apVOuter[iNOrigVerts - 1], pVNorth, apVOuter[1])
      != checkOrient3D(pVOther, pVSouth, apVOuter[iNOrigVerts - 1],
		       apVOuter[1]))
    return 0;

  double dOrigQual;
  //@@ Compute the quality of the current configuration.
  // This is the smallest quality for any of the existing cells.
  dOrigQual = 1000;
  for (i = 0; i < iNOrigCells; i++) {
    Cell *pC = apCOrig[i];
    // Evaluate quality, depending on the swap criterion in use.
    double dQual = evalQualFunc(pC->getVert(0), pC->getVert(1), pC->getVert(2),
				pC->getVert(3));
    dOrigQual = min(dQual, dOrigQual);
  }

  // Set up info about the possible new configurations.

  // Compute quality for each possible tet pair (one pair per possible
  // triangle).
  double a3dNewQual[10][10][10];
  int iSign = checkOrient3D(apVOuter[0], apVOuter[1], pVSouth, pVNorth);
  int iV0, iV1, iV2;
#ifndef NDEBUG
  for (iV0 = 0; iV0 < 10; iV0++)
    for (iV1 = 0; iV1 < 10; iV1++)
      for (iV2 = 0; iV2 < 10; iV2++)
	a3dNewQual[iV0][iV1][iV2] = 0;
  for (int iDum = 0; iDum < iNOrigCells; iDum++)
    assert(
	iSign
	    == checkOrient3D(apVOuter[iDum], apVOuter[iDum + 1], pVSouth,
			     pVNorth));
#endif
  for (iV0 = 0; iV0 < iNOrigVerts; iV0++)
    for (iV1 = iV0 + 1; iV1 < iNOrigVerts; iV1++)
      for (iV2 = iV1 + 1; iV2 < iNOrigVerts; iV2++) {
	if (checkOrient3D(apVOuter[iV0], apVOuter[iV1], apVOuter[iV2], pVNorth)
	    == iSign)
	  // Evaluate quality
	  a3dNewQual[iV0][iV1][iV2] = evalQualFunc(apVOuter[iV0], apVOuter[iV1],
						   apVOuter[iV2], pVNorth);
	else
	  a3dNewQual[iV0][iV1][iV2] = dInvalidPenalty;

	// Don't bother evaluating for the other tet if this one is
	// worse that the worst in the original config, because this tet
	// pair can never be in a final config anyway.
	if (a3dNewQual[iV0][iV1][iV2] >= dOrigQual) {
	  double dDummyQual;
	  if (checkOrient3D(apVOuter[iV0], apVOuter[iV2], apVOuter[iV1],
			    pVSouth) == iSign)
	    dDummyQual = evalQualFunc(apVOuter[iV0], apVOuter[iV2],
				      apVOuter[iV1], pVSouth);
	  else
	    dDummyQual = dInvalidPenalty;
	  a3dNewQual[iV0][iV1][iV2] = min(a3dNewQual[iV0][iV1][iV2],
					  dDummyQual);
	}

	// Copy things around symmetrically
	a3dNewQual[iV2][iV1][iV0] = a3dNewQual[iV2][iV0][iV1] =
	    a3dNewQual[iV1][iV2][iV0] = a3dNewQual[iV1][iV0][iV2] =
		a3dNewQual[iV0][iV2][iV1] = a3dNewQual[iV0][iV1][iV2];

      }

  // Find the new configuration with the highest quality.
  double dMaxQual = dOrigQual;
  int iCanon, iPerm;
  int iBestCanon = -1, iBestPerm = -1;
  for (iCanon = 0; iCanon < aCS[iNOrigVerts].iNCanon; iCanon++)
    for (iPerm = 0; iPerm < aCS[iNOrigVerts].aConf[iCanon].iNPerm; iPerm++) {
      double dConfQual = 1000;
      // Once this config is known to be worse than the best one we know
      // about, we can give up on it.
      for (int iTri = 0; iTri < iNOrigVerts - 2 && (dConfQual >= dMaxQual);
	  iTri++) {
	iV0 = (aCS[iNOrigVerts].aConf[iCanon].a2iFaceVert[iTri][0] + iPerm)
	    % iNOrigVerts;
	iV1 = (aCS[iNOrigVerts].aConf[iCanon].a2iFaceVert[iTri][1] + iPerm)
	    % iNOrigVerts;
	iV2 = (aCS[iNOrigVerts].aConf[iCanon].a2iFaceVert[iTri][2] + iPerm)
	    % iNOrigVerts;
	dConfQual = min(dConfQual, a3dNewQual[iV0][iV1][iV2]);
      }
      if (dConfQual > dMaxQual) {
	dMaxQual = dConfQual;
	iBestCanon = iCanon;
	iBestPerm = iPerm;
      }
    }

  // If the highest-quality alternative is worse than the current
  // configuration, do nothing.  The check for iBestCanon is there so
  // that some stupid roundoff in comparison of dMaxQual and dOrigQual
  // won't cause a reconfiguration when there shouldn't be one.
  if (dMaxQual <= dOrigQual || iBestCanon == -1) {
    logMessage(4, "Bdry edge failed: %d\n", iNOrigCells);
    return 0;
  }
  else {
    // If the highest-quality alternative is better than the current
    // configuration:
    aiEdgeDone[iNOrigCells]++;
    logMessage(2, "Bdry edge swap: %d\n", iNOrigCells);

    // Delete all of the interior tets and faces.  These must be
    // separate loops, or forced.  Outer faces are disconnected
    // automatically.
    for (i = 0; i < iNOrigCells; i++) {
      deleteCell(apCOrig[i]);
    }

    // Renumber all the old vertex data, offsetting it by iBestPerm to
    // make the desired configuration match the stored canonical
    // configuration.
    if (iBestPerm != 0) {
      Vert *apVOuterTemp[11];
      // Copy data to temporary arrays.
      for (i = 0; i < iNOrigVerts; i++) {
	apVOuterTemp[i] = apVOuter[i];
      }
      // Copy back with an offset.
      for (i = 0; i < iNOrigVerts; i++) {
	apVOuter[i] = apVOuterTemp[(i + iBestPerm) % iNOrigVerts];
      }
    }

    // Re-connect the mesh locally
    Config *pConf = &(aCS[iNOrigVerts].aConf[iBestCanon]);

    // Set up an array to hold the new cells
    Cell** apCNew = new Cell*[2 * (iNOrigVerts - 2)];

    int iOrient = 0;
    // There are iNOrigVerts - 2 triangles on the equator, two tets apiece.
    for (i = 0; i < iNOrigVerts - 2; i++) {
      // These are the equatorial tris, so all these indices should be
      // positive.
      assert(pConf->a2iFaceVert[i][0] >= 0);
      assert(pConf->a2iFaceVert[i][1] >= 0);
      assert(pConf->a2iFaceVert[i][2] >= 0);
      Vert* pVert0 = apVOuter[pConf->a2iFaceVert[i][0]];
      Vert* pVert1 = apVOuter[pConf->a2iFaceVert[i][1]];
      Vert* pVert2 = apVOuter[pConf->a2iFaceVert[i][2]];

      if (iOrient == 0)
	iOrient = checkOrient3D(pVert0, pVert1, pVert2, pVNorth);

      if (iOrient == 1) {
	assert(checkOrient3D(pVert0, pVert1, pVert2, pVNorth) == 1);
	assert(checkOrient3D(pVert0, pVert2, pVert1, pVSouth) == 1);
	bool qExist = false;
	apCNew[2 * i] = createTetCell(qExist, pVert0, pVert1, pVert2, pVNorth,
				      iReg);
	assert(!qExist);
	apCNew[2 * i + 1] = createTetCell(qExist, pVert0, pVert2, pVert1,
					  pVSouth, iReg);
	assert(!qExist);
      }
      else {
	assert(checkOrient3D(pVert0, pVert2, pVert1, pVNorth) == 1);
	assert(checkOrient3D(pVert0, pVert1, pVert2, pVSouth) == 1);
	bool qExist = false;
	apCNew[2 * i] = createTetCell(qExist, pVert0, pVert2, pVert1, pVNorth,
				      iReg);
	assert(!qExist);
	apCNew[2 * i + 1] = createTetCell(qExist, pVert0, pVert1, pVert2,
					  pVSouth, iReg);
	assert(!qExist);
      }
    }

    // Retrieve the faces that were just created on the boundary.
    {
      bool qExisted = true;
      Vert *pV0, *pV1;
      if (iBestPerm == 0) {
	pV0 = apVOuter[0];
	pV1 = apVOuter[iNOrigVerts - 1];
      }
      else {
	pV0 = apVOuter[iNOrigVerts - iBestPerm];
	pV1 = apVOuter[iNOrigVerts - iBestPerm - 1];
      }
      apFNorth[iNOrigVerts - 1] = createFace(qExisted, pV0, pV1, pVNorth);
      assert(qExisted); // Just created them; they'd -better- exist!

      apFSouth[iNOrigVerts - 1] = createFace(qExisted, pV1, pV0, pVSouth);
      assert(qExisted); // Just created them; they'd -better- exist!
    }

    // Attach the new faces to the bdry faces.
    BFace* pBFBFaceNew0 = createBFace(apFNorth[iNOrigVerts - 1], pBFBFace0);
    BFace* pBFBFaceNew1 = createBFace(apFSouth[iNOrigVerts - 1], pBFBFace1);
    assert(
	!m_strictPatchChecking
	    || (pBFBFaceNew0->getPatchPointer()
		== pBFBFaceNew1->getPatchPointer()));

    // Delete the old faces on the boundary.
    deleteBFace(pBFBFace0);
    deleteBFace(pBFBFace1);

#ifndef OMIT_VERTEX_HINTS
    // Reset vertex hints for polar verts
    pVNorth->setHintFace(apFNorth[0]);
    pVSouth->setHintFace(apFSouth[0]);
#endif

    // Validate the result
    for (i = 0; i < iNOrigVerts * 2 - 4; i++) {
      // Check everything about the tet
      assert(apCNew[i]->doFullCheck());
      // Every tet has to hit one pole or the other
      assert(apCNew[i]->hasVert(pVNorth) || apCNew[i]->hasVert(pVSouth));
      // Assert tetrahedron orientation
      assert(
	  checkOrient3D(apCNew[i]->getVert(0), apCNew[i]->getVert(1),
			apCNew[i]->getVert(2), apCNew[i]->getVert(3)) == 1);
    }

    for (i = 0; i < iNOrigVerts; i++) {
      assert(apFSouth[i]->doFullCheck());
      assert(apFNorth[i]->doFullCheck());
    }

    // If recursion is permitted, check all the new faces for possible
    // swaps, as well as all of the faces on the outside of the orbit.
    int iRetVal = iNOrigVerts - 2;
    iSwapDepth++;
    if (iSwapDepth > iMaxSwapDepth)
      logMessage(3, "Swap recursion depth exceeded!\n");

    Face** apFRecur = new Face*[(2 * iNOrigVerts - 4) * 4];
    // Grab all the faces for all the tets that were just created...
    Face **ppFCurr = apFRecur;
    for (i = 0; i < 2 * iNOrigVerts - 4; i++) {
      Cell *pC = apCNew[i];
      for (int j = 0; j < 4; j++) {
	(*ppFCurr) = pC->getFace(j);
	ppFCurr++;
      }
    }
    assert(ppFCurr - apFRecur == (2 * (iNOrigVerts - 2) * 4));

    // Now sort them uniquely.
    std::sort(apFRecur, ppFCurr);
    Face **ppFNewEnd = std::unique(apFRecur, ppFCurr);
    assert(ppFNewEnd - apFRecur == (5 * iNOrigVerts - 8));

    // And now swap them all.
    for (i = 0; i < 5 * iNOrigVerts - 8; i++) {
      if (m_swapRecursively && iSwapDepth <= iMaxSwapDepth) {
	iRetVal += iFaceSwap_deprecated(apFRecur[i]);
      }
    }
    delete[] apFRecur;
    delete[] apCNew;

    // Return the total (equivalent) number of swaps performed
    iSwapDepth--;
    return iRetVal;
  }
}

//@ Swap away face if this is legal and improves the mesh quality measure
// Swapping typically propagates through the mesh, and the return value from
// iFaceSwap is the total number of swaps done during this invocation.
int
VolMesh::iFaceSwap_deprecated(Face*& pF)
{
  if (!pF->isSwapAllowed())
    return (0);

  // If the face has been removed, don't bother
  assert(pF->doFullCheck()); // This should always be true, I'm pretty sure.

  if (iSwapDepth > iMaxSwapDepth)
    return (0);
  iSwapDepth++;

  Vert *pVVertA, *pVVertB, *pVVertC, *pVVertD, *pVVertE;
  Vert *pVPivot0, *pVPivot1, *pVOther;
  TetCell* apTCTets[4];
  int iNTets;
  TriFace* pTF = (dynamic_cast<TriFace*>(pF));
  if (!pTF->isValid())
    return 0;
  eFC = pTF->categorizeFace(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE,
			    apTCTets, iNTets, pVPivot0, pVPivot1, pVOther);

  //@@ Determine the preferable configuration and swap if necessary
  switch (eFC)
    {
    default: // Catch-all for bad cases
      vWarning("Shouldn't be here: face wasn't categorized");
      iSwapDepth--;
      return 0;
      // These cases are handled by edge swapping.
    case eN44:
    case eN32:
      {
	int iRetVal = 0;
	if (getSwapType() != eDelaunay && isEdgeSwappingAllowed())
	  iRetVal = edgeSwap3D(pF, pVPivot0, pVPivot1, pVOther);
	iSwapDepth--;
	return iRetVal;
      }
      // These cases are definitely unswappable
    case eBdry:
    case eN40:
    case eN30:
    case eN20:
    case eBdryReflex:
    case eOther:
      iSwapDepth--;
      return (0);
    case eT44:
    case eT22:
    case eT23:
    case eT32:
      if (qDoSwap(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE)) {
	int iRes = reconfigure_deprecated(apTCTets, dynamic_cast<TriFace*>(pF),
					  pVVertA, pVVertB, pVVertC, pVVertD,
					  pVVertE);
	iSwapDepth--;
	return iRes;
      }
      else {
	iSwapDepth--;
	return (0);
      }
    }
}

//@@ Determines whether swapping is needed for local Delaunay-ness.
// Returns true to swap, false to remain the same.
static bool
qDoSwapDelaunay(const Vert* pVVertA, const Vert* pVVertB, const Vert* pVVertC,
		const Vert* pVVertD, const Vert* pVVertE)
{
  int iInSphereRes = isInsphere(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE);

  switch (eFC)
    {
    case eT32:
      // Current tets are ABDE, BCDE, CADE.  Swap if pVVertE is on (0) or
      // outside (-1) the circumsphere of ABCD.  This is an arbitrary
      // tie-break; configurations which are co-spherical go to the two
      // tet configuration.
      return (iInSphereRes != 1);
    case eT23:
    case eT22:
    case eT44: // This one is the same as eT22 for all situations.
	       // Current tets are ABCD, CBAE.  Swap if pVVertE is inside (1) the
	       // circumsphere of ABCD.  Again, configurations which are
	       // co-spherical go to the two-tet configuration.
      return (iInSphereRes == 1);
    default:
      return (false);
    }
}

//@@ Determines whether swapping will improve the maximum face angle.
// Returns true to swap, false to remain the same.
static bool
qDoSwapMaxDihed(const Vert* pVVertA, const Vert* pVVertB, const Vert* pVVertC,
		const Vert* pVVertD, const Vert* pVVertE)
{
  // This algorithm finds face angles between adjacent faces by dotting
  // their unit normals.  The largest magnitude loses.
  //
  // To prevent pairs from flopping back and forth, a swap is only
  // authorized if the inequality is bigger than eps.

  double dEps = 1.e-10;
  switch (eFC)
    {
    case eT23:
    case eT32:
      {
	//@@@ The case where the five pVVerts form a genuinely convex set
	double adNormABD[3], adNormABE[3], adNormADE[3];
	double adNormBCD[3], adNormBCE[3], adNormBDE[3];
	double adNormCAD[3], adNormCAE[3], adNormCDE[3];
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertD->getCoords(), adNormABD);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);
	calcUnitNormal(pVVertA->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormADE);

	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormBCE);
	calcUnitNormal(pVVertB->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormBDE);

	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertD->getCoords(), adNormCAD);
	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertE->getCoords(), adNormCAE);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);

	double dLastDot, dWorst2 = -1., dWorst3 = -1.;

	// Care must be taken to be certain that I'm taking the dot
	// product of two vectors that both point into the enclosed angle
	// or both point out, or correct the sign.

	dLastDot = dDOT3D(adNormABD, adNormBCD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormCAD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormCAD, adNormABD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);

	dLastDot = dDOT3D(adNormABE, adNormBCE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormBCE, adNormCAE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormCAE, adNormABE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);

	dLastDot = -dDOT3D(adNormABD, adNormABE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormBCD, adNormBCE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormCAD, adNormCAE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);

	dLastDot = -dDOT3D(adNormADE, adNormBDE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormBDE, adNormCDE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormCDE, adNormADE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);

	if (dWorst3 > dWorst2 + dEps)
	  return (eFC == eT32);
	else if (dWorst2 > dWorst3 + dEps)
	  return (eFC == eT23);
	else
	  return (false);
      }
    case eT22:
      {
	double adNormACE[3], adNormECB[3], adNormBCD[3], adNormDCA[3];
	double adNormABC[3], adNormCDE[3], adNormABE[3];
	calcUnitNormal(pVVertA->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormACE);
	calcUnitNormal(pVVertE->getCoords(), pVVertC->getCoords(),
		       pVVertB->getCoords(), adNormECB);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertD->getCoords(), pVVertC->getCoords(),
		       pVVertA->getCoords(), adNormDCA);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);

	double dLastDot, dWorstNow = -1., dWorstSwapped = -1.;

	dLastDot = dDOT3D(adNormACE, adNormECB);
	dWorstNow = max(dWorstNow, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormDCA);
	dWorstNow = max(dWorstNow, dLastDot);
	dLastDot = fabs(dDOT3D(adNormABC, adNormABE));
	dWorstNow = max(dWorstNow, dLastDot);

	dLastDot = dDOT3D(adNormACE, adNormDCA);
	dWorstSwapped = max(dWorstSwapped, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormECB);
	dWorstSwapped = max(dWorstSwapped, dLastDot);
	dLastDot = fabs(dDOT3D(adNormCDE, adNormABE));
	dWorstSwapped = max(dWorstSwapped, dLastDot);

	if (dWorstNow > dWorstSwapped + dEps)
	  return (true);
	else
	  return (false);
      }
    case eT44:
      return (false);
    default:
      assert2(0, "Bad face type in finding max dihedral angle");
      return (false);
    }
}

//@@ Determines whether swapping will improve the maximum face angle.
// Returns true to swap, false to remain the same.
static bool
qDoSwapMinSine(const Vert* pVVertA, const Vert* pVVertB, const Vert* pVVertC,
	       const Vert* pVVertD, const Vert* pVVertE)
{
  // This algorithm finds sine of face angles between adjacent faces by
  // crossing their unit normals.  The idea is to pick the configuration
  // which maximizes this value.
  //
  // To prevent pairs from flopping back and forth, a swap is only
  // authorized if the inequality is bigger than eps.

  double dEps = 1.e-10;
  switch (eFC)
    {
    case eT23:
    case eT32:
      {
	double adNormABD[3], adNormABE[3], adNormADE[3], adNormABC[3];
	double adNormBCD[3], adNormBCE[3], adNormBDE[3];
	double adNormCAD[3], adNormCAE[3], adNormCDE[3];
	//@@@ The case where the five pVVerts form a genuinely convex set
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertD->getCoords(), adNormABD);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);
	calcUnitNormal(pVVertA->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormADE);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC); // Points to E

	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormBCE);
	calcUnitNormal(pVVertB->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormBDE);

	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertD->getCoords(), adNormCAD);
	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertE->getCoords(), adNormCAE);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);

	double dWorst2 = 1., dWorst3 = 1.;

	// Because I'm taking the MAGNITUDE of the cross product, the
	// actual sign of the result doesn't matter (vector can be
	// reversed).  This means that normal orientation doesn't matter.

	// Sines of dihedral angles of tet ABCD
	double adTemp[3], dMag;

	vCROSS3D(adNormABD, adNormBCD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCD, adNormCAD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAD, adNormABD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormABD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	// Sines of dihedral angles of tet ABCE
	vCROSS3D(adNormABE, adNormBCE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCE, adNormCAE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAE, adNormABE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormABE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	// Now have the smallest sine of dihedral angle in tet ABCD and ABCE

	//Sines of dihedral angles of tet ABDE
	vCROSS3D(adNormABD, adNormABE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABD, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABD, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormADE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	//Sines of dihedral angles of tet BCDE
	vCROSS3D(adNormBCD, adNormBCE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCD, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBDE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	//Sines of dihedral angles of tet CADE
	vCROSS3D(adNormCAD, adNormCAE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAD, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCDE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);
	// Now have the smallest sine of dihedral in tets ABDE, BCDE, CADE

	if (dWorst3 > dWorst2 + dEps)
	  return (eFC == eT23);
	else if (dWorst2 > dWorst3 + dEps)
	  return (eFC == eT32);
	else
	  return (false);
      }
    case eT22:
      {
	double adNormACE[3], adNormECB[3], adNormBCD[3], adNormDCA[3];
	double adNormABC[3], adNormCDE[3], adNormABE[3];
	calcUnitNormal(pVVertA->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormACE);
	calcUnitNormal(pVVertE->getCoords(), pVVertC->getCoords(),
		       pVVertB->getCoords(), adNormECB);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertD->getCoords(), pVVertC->getCoords(),
		       pVVertA->getCoords(), adNormDCA);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);

	double dWorstNow = 1., dWorstSwap = 1.;
	double adTemp[3], dMag;

	// Sines of dihedral angles of tet ABCD, except those that are the
	// same both ways
	vCROSS3D(adNormBCD, adNormDCA, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormBCD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormDCA, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Sines of dihedral angles of tet ABCE, except those that are the
	// same both ways
	vCROSS3D(adNormECB, adNormACE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormECB, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormACE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Only check the center face once; the sine is the same on both
	// sides
	vCROSS3D(adNormABE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Now have the smallest sine of dihedral angle in tet ABCD and ABCE

	// Sines of dihedral angles of tet DECA, except those that are the
	// same both ways
	vCROSS3D(adNormDCA, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormACE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormDCA, adNormACE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Sines of dihedral angles of tet DECB, except those that are the
	// same both ways
	vCROSS3D(adNormBCD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormECB, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormBCD, adNormECB, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Only check the center face once; the sine is the same on both
	// sides
	vCROSS3D(adNormABE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Now have the smallest sine of dihedral angle in tet DECA and DECB

	if (dWorstNow + dEps < dWorstSwap)
	  return (true);
	else
	  return (false);
      }
    case eT44:
      return (false);
    default:
      assert2(0, "Bad face type in finding max dihedral angle");
      return (false);
    }
}

bool
VolMesh::qDoSwap(const Vert* const pVVertA, const Vert* const pVVertB,
		 const Vert* const pVVertC, const Vert* const pVVertD,
		 const Vert* const pVVertE) const
{
  assert(pVVertA->isValid());
  assert(pVVertB->isValid());
  assert(pVVertC->isValid());
  assert(pVVertD->isValid());
  assert(pVVertE->isValid());

  switch (getSwapType())
    {
    case eDelaunay:
      return qDoSwapDelaunay(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE);
    case eMaxDihed:
      return qDoSwapMaxDihed(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE);
    case eMinSine:
      return qDoSwapMinSine(pVVertA, pVVertB, pVVertC, pVVertD, pVVertE);
    case eFromQual:
      {
	if (eFC != eT23 && eFC != eT32)
	  return false;
	assert(m_qual != NULL);
	TetCellCV TCCV11(const_cast<Vert*>(pVVertA), const_cast<Vert*>(pVVertC),
			 const_cast<Vert*>(pVVertB),
			 const_cast<Vert*>(pVVertD));
	TetCellCV TCCV12(const_cast<Vert*>(pVVertA), const_cast<Vert*>(pVVertB),
			 const_cast<Vert*>(pVVertC),
			 const_cast<Vert*>(pVVertE));
	double dQual11 = m_qual->dEvaluate(&TCCV11);
	double dQual12 = m_qual->dEvaluate(&TCCV12);

	TetCellCV TCCV21(const_cast<Vert*>(pVVertD), const_cast<Vert*>(pVVertE),
			 const_cast<Vert*>(pVVertA),
			 const_cast<Vert*>(pVVertB));
	TetCellCV TCCV22(const_cast<Vert*>(pVVertD), const_cast<Vert*>(pVVertE),
			 const_cast<Vert*>(pVVertB),
			 const_cast<Vert*>(pVVertC));
	TetCellCV TCCV23(const_cast<Vert*>(pVVertD), const_cast<Vert*>(pVVertE),
			 const_cast<Vert*>(pVVertC),
			 const_cast<Vert*>(pVVertA));
	double dQual21 = m_qual->dEvaluate(&TCCV21);
	double dQual22 = m_qual->dEvaluate(&TCCV22);
	double dQual23 = m_qual->dEvaluate(&TCCV23);

	bool qMaximize = m_qual->qMaximizeMinValue();
	static const double dEps = 1.e-10;
	if (qMaximize) {
	  double dQual1 = min(dQual11, dQual12);
	  double dQual2 = MIN3(dQual21, dQual22, dQual23);
	  if (dQual1 > dQual2 + dEps)
	    return (eFC == eT32);
	  else if (dQual2 > dQual1 + dEps)
	    return (eFC == eT23);
	  else
	    return false;
	}
	else {
	  double dQual1 = max(dQual11, dQual12);
	  double dQual2 = MAX3(dQual21, dQual22, dQual23);
	  if (dQual1 < dQual2 - dEps)
	    return (eFC == eT32);
	  else if (dQual2 < dQual1 - dEps)
	    return (eFC == eT23);
	  else
	    return false;
	}
      }
    default:
      assert2(0, "Bad swap measure specified");
      return (false);
    }
}

int
VolMesh::reconfigure_deprecated(TetCell* apTCTets[], TriFace* const pF,
				Vert* const pVVertA, Vert* const pVVertB,
				Vert* const pVVertC, Vert* const pVVertD,
				Vert* const pVVertE)
{
  assert(pF->doFullCheck());
  assert(pVVertA->isValid());
  assert(pVVertB->isValid());
  assert(pVVertC->isValid());
  assert(pVVertD->isValid());
  assert(pVVertE->isValid());

  switch (eFC)
    {
    case eT22:
      return (reconfTet22_deprecated(apTCTets, pF, pVVertA, pVVertB, pVVertC,
				     pVVertD, pVVertE));
    case eT23:
      return (reconfTet23_deprecated(apTCTets, pF, pVVertA, pVVertB, pVVertC,
				     pVVertD, pVVertE));
    case eT32:
      return (reconfTet32_deprecated(apTCTets, pF, pVVertA, pVVertB, pVVertC,
				     pVVertD, pVVertE));
    case eT44:
      return (reconfTet44_deprecated(apTCTets, pF, pVVertA, pVVertB, pVVertC,
				     pVVertD, pVVertE));
    default:
      assert2(0, "Bad face type");
      return (0);
    }
}

//@@ Swap two tets for three.
int
VolMesh::reconfTet23_deprecated(TetCell* apTCTets[], TriFace* const pF,
				Vert* const pVVertA, Vert* const pVVertB,
				Vert* const pVVertC, Vert* const pVVertD,
				Vert* const pVVertE)
{
  if (pF->isLocked())
    return 0;
  // Identify both cells and all seven faces in the two-tet config.
  TetCell* pTCTetD = apTCTets[0];
  TetCell* pTCTetE = apTCTets[1];
  if (pTCTetD->getRegion() != pTCTetE->getRegion())
    return 0;
  assert(pTCTetD->hasVert(pVVertD));
  assert(pTCTetE->hasVert(pVVertE));
  assert(pTCTetD == pF->getLeftCell());
  assert(pTCTetE == pF->getRightCell());
  //   Face* pFFaceD = pF;

  // Top half
  Face* pFFaceA = pTCTetE->getOppositeFace(pVVertA);
  Face* pFFaceB = pTCTetE->getOppositeFace(pVVertB);
  Face* pFFaceC = pTCTetE->getOppositeFace(pVVertC);

  // Bottom half
  Face* pFFaceA1 = pTCTetD->getOppositeFace(pVVertA);
  Face* pFFaceB1 = pTCTetD->getOppositeFace(pVVertB);
  Face* pFFaceC1 = pTCTetD->getOppositeFace(pVVertC);

  int iRegion = pTCTetD->getRegion();
  // Delete old tets
  deleteCell(pTCTetD);
  deleteCell(pTCTetE);

  // Create three new tets
  bool qExist;
  TetCell *pTCTetA = createTetCell(qExist, pVVertB, pVVertC, pVVertD, pVVertE,
				   iRegion);
  assert(!qExist);
  TetCell *pTCTetB = createTetCell(qExist, pVVertC, pVVertA, pVVertD, pVVertE,
				   iRegion);
  assert(!qExist);
  TetCell *pTCTetC = createTetCell(qExist, pVVertA, pVVertB, pVVertD, pVVertE,
				   iRegion);
  assert(!qExist);

  assert(!pTCTetA->hasVert(pVVertA));
  assert(!pTCTetB->hasVert(pVVertB));
  assert(!pTCTetC->hasVert(pVVertC));
  assert(checkOrient3D(pVVertA, pVVertB, pVVertD, pVVertE) == 1);
  assert(checkOrient3D(pVVertB, pVVertC, pVVertD, pVVertE) == 1);
  assert(checkOrient3D(pVVertC, pVVertA, pVVertD, pVVertE) == 1);

  assert(pTCTetA->doFullCheck());
  assert(pTCTetB->doFullCheck());
  assert(pTCTetC->doFullCheck());

#ifndef NDEBUG
  Face *pFFaceA2 = findCommonFace(pTCTetB, pTCTetC);
  Face *pFFaceB2 = findCommonFace(pTCTetC, pTCTetA);
  Face *pFFaceC2 = findCommonFace(pTCTetA, pTCTetB);
  assert(pFFaceA->doFullCheck());
  assert(pFFaceB->doFullCheck());
  assert(pFFaceC->doFullCheck());
  assert(pFFaceA1->doFullCheck());
  assert(pFFaceB1->doFullCheck());
  assert(pFFaceC1->doFullCheck());
  assert(pFFaceA2->doFullCheck());
  assert(pFFaceB2->doFullCheck());
  assert(pFFaceC2->doFullCheck());
#endif

  int iRetVal = 1;
  if (m_swapRecursively && iSwapDepth <= iMaxSwapDepth) {
    iRetVal += iFaceSwap_deprecated(pFFaceA);
    iRetVal += iFaceSwap_deprecated(pFFaceB);
    iRetVal += iFaceSwap_deprecated(pFFaceC);
    iRetVal += iFaceSwap_deprecated(pFFaceA1);
    iRetVal += iFaceSwap_deprecated(pFFaceB1);
    iRetVal += iFaceSwap_deprecated(pFFaceC1);
  }
  return (iRetVal);
}
//@@ Swap three tets for two.
int
VolMesh::reconfTet32_deprecated(TetCell* apTCTets[], TriFace* const pF,
				Vert* const pVVertA, Vert* const pVVertB,
				Vert* const pVVertC, Vert* const pVVertD,
				Vert* const pVVertE)
{
  // Identify the three cells of the three-tet config.  The cell ID's
  // require that the cell NOT containing pVVertA is pTCTetA, etc.
  TetCell* pTCTetA = apTCTets[0];
  TetCell* pTCTetB = apTCTets[1];
  TetCell* pTCTetC = apTCTets[2];
  assert(!pTCTetA->hasVert(pVVertA));
  assert(!pTCTetB->hasVert(pVVertB));
  assert(!pTCTetC->hasVert(pVVertC));
  int iRegion = pTCTetA->getRegion();
  if ((iRegion != pTCTetB->getRegion()) || (iRegion != pTCTetC->getRegion()))
    return 0;

  assert(pF->hasVert(pVVertA) || pF->hasVert(pVVertB) || pF->hasVert(pVVertC));
  assert(pF->hasVert(pVVertD));
  assert(pF->hasVert(pVVertE));

  // Identify the nine faces in the three-tet configuration. pFFaceA and
  // pFFaceA1 are faces in the shell not incident on pVVertA. pFFaceA2 IS
  // incident on pVVertA, because that's the easiest way to define it
  // uniquely.
  Face* pFFaceA = pTCTetA->getOppositeFace(pVVertE);
  Face* pFFaceB = pTCTetB->getOppositeFace(pVVertE);
  Face* pFFaceC = pTCTetC->getOppositeFace(pVVertE);

  Face* pFFaceA1 = pTCTetA->getOppositeFace(pVVertD);
  Face* pFFaceB1 = pTCTetB->getOppositeFace(pVVertD);
  Face* pFFaceC1 = pTCTetC->getOppositeFace(pVVertD);

  Face* pFFaceA2 = pTCTetC->getOppositeFace(pVVertB);
  Face* pFFaceB2 = pTCTetA->getOppositeFace(pVVertC);
  Face* pFFaceC2 = pTCTetB->getOppositeFace(pVVertA);

  if (pFFaceA2->isLocked() || pFFaceB2->isLocked() || pFFaceC2->isLocked())
    return 0;

  // Delete tets
  deleteCell(pTCTetA);
  deleteCell(pTCTetB);
  deleteCell(pTCTetC);

  // Create new tets
  bool qExist;
  TetCell *pTCTetD = createTetCell(qExist, pVVertA, pVVertC, pVVertB, pVVertD,
				   iRegion);
  assert(!qExist);
  TetCell *pTCTetE = createTetCell(qExist, pVVertA, pVVertB, pVVertC, pVVertE,
				   iRegion);
  assert(!qExist);

  assert(pTCTetD->hasVert(pVVertD));
  assert(pTCTetE->hasVert(pVVertE));
  assert(checkOrient3D(pVVertA, pVVertC, pVVertB, pVVertD) == 1);
  assert(checkOrient3D(pVVertA, pVVertB, pVVertC, pVVertE) == 1);

  assert(pTCTetD->doFullCheck());
  assert(pTCTetE->doFullCheck());

  assert(pFFaceA->doFullCheck());
  assert(pFFaceB->doFullCheck());
  assert(pFFaceC->doFullCheck());
#ifndef NDEBUG
  Face *pFFaceD = findCommonFace(pTCTetD, pTCTetE);
  assert(pFFaceD->doFullCheck());
#endif
  assert(pFFaceA1->doFullCheck());
  assert(pFFaceB1->doFullCheck());
  assert(pFFaceC1->doFullCheck());

  int iRetVal = 1;
  if (m_swapRecursively && iSwapDepth <= iMaxSwapDepth) {
    iRetVal += iFaceSwap_deprecated(pFFaceA);
    iRetVal += iFaceSwap_deprecated(pFFaceB);
    iRetVal += iFaceSwap_deprecated(pFFaceC);
    iRetVal += iFaceSwap_deprecated(pFFaceA1);
    iRetVal += iFaceSwap_deprecated(pFFaceB1);
    iRetVal += iFaceSwap_deprecated(pFFaceC1);
  }
  return (iRetVal);
}

//@@ Swap two tets for two, in the case where two faces are coplanar.
int
VolMesh::reconfTet22_deprecated(TetCell* apTCTets[], TriFace* const pF,
				Vert* const pVVertA, Vert* const pVVertB,
				Vert* const pVVertC, Vert* const pVVertD,
				Vert* const pVVertE)
{
  // This routine does only the 2-to-2 case.  4-to-4 is done elsewhere.
  // The cells and verts are assumed to have been set up properly already.

  /*      	        C
   /\
        	      /||\
        	     / |  \                 A,D are on the base;
   /  | | \                B,E are front, C,F are back.
   /   |    \               Cell A is on the left, B on
   /    |  |  \	            the right.
   /     |      \
  		/      |   |   \
  	       /     C |     F  \
  	      /	       |  E |    \
  	     /   B     | _ - B    \
  	    /  _ _ - - |    /  - _ \
  	   D -         |           -E
   -__	    A  |  /  D  __--
   -__     |     __-
   -___ |/__--
   -A-

   */
  if (!areSurfaceEdgeChangesAllowed() || !areBdryChangesAllowed())
    return 0;
  TetCell* pTCTetA = apTCTets[0];
  TetCell* pTCTetB = apTCTets[1];
  assert(pTCTetA->hasVert(pVVertD));
  assert(pTCTetB->hasVert(pVVertE));
  assert(pTCTetA == pF->getLeftCell());
  assert(pTCTetB == pF->getRightCell());
  assert(pTCTetA->getRegion() == pTCTetB->getRegion());

  Face* pFFaceA = pTCTetA->getOppositeFace(pVVertC);
  Face* pFFaceB = pTCTetA->getOppositeFace(pVVertB);
  Face* pFFaceC = pTCTetA->getOppositeFace(pVVertA);

  Face* pFFaceD = pTCTetB->getOppositeFace(pVVertC);
  Face* pFFaceE = pTCTetB->getOppositeFace(pVVertB);
  Face* pFFaceF = pTCTetB->getOppositeFace(pVVertA);

  if (pF->isLocked() || pFFaceA->isLocked() || pFFaceD->isLocked())
    return 0;

  // Check that faces A and D are bdry faces.  Then, if they belong to
  // different patches, bail out, so that we won't end up losing a
  // subsegment (from A to B, in this case) when we swap.
  Cell* pCOppA = pFFaceA->getOppositeCell(pTCTetA);
  Cell* pCOppD = pFFaceD->getOppositeCell(pTCTetB);
  assert(pCOppA->getType() == Cell::eTriBFace);
  assert(pCOppD->getType() == Cell::eTriBFace);
  TriBFace* pTBFA = dynamic_cast<TriBFace*>(pCOppA);
  TriBFace* pTBFD = dynamic_cast<TriBFace*>(pCOppD);
  if (!areBdryFacesOKToSwap(pTBFA, pTBFD))
    return 0;

  // Do the actual two-for-two swap here.  The faces all stay in the
  // same places except for pFFace, the dividing face, which gets turned.
  // Cell A ends up in front, cell B in back.  The only vert hints
  // affected are for verts 1 and 2, which can be set explicitly.

  int iRegion = pTCTetA->getRegion();
  assert(iRegion == pTCTetB->getRegion());

  // Delete old cells
  deleteCell(pTCTetA);
  deleteCell(pTCTetB);

  bool qExist;
  pTCTetA = createTetCell(qExist, pVVertD, pVVertA, pVVertE, pVVertC, iRegion);
  assert(!qExist);
  pTCTetB = createTetCell(qExist, pVVertE, pVVertB, pVVertD, pVVertC, iRegion);
  assert(!qExist);

  Face *pFBdryA = pTCTetA->getOppositeFace(pVVertC);
  Face *pFBdryB = pTCTetB->getOppositeFace(pVVertC);
  createBFace(pFBdryA, pTBFA);
  createBFace(pFBdryB, pTBFA);

  // Delete bdry faces now
  deleteBFace(pTBFA);
  deleteBFace(pTBFD);

  assert(pTCTetA->doFullCheck());
  assert(pTCTetB->doFullCheck());

#ifndef NDEBUG
  Face *pFNewMid = findCommonFace(pTCTetA, pTCTetB);
  pFFaceA = pTCTetA->getOppositeFace(pVVertC);
  pFFaceD = pTCTetB->getOppositeFace(pVVertC);
  assert(pFNewMid->doFullCheck());
  assert(pFFaceA->doFullCheck());
  assert(pFFaceB->doFullCheck());
  assert(pFFaceC->doFullCheck());
  assert(pFFaceD->doFullCheck());
  assert(pFFaceE->doFullCheck());
  assert(pFFaceF->doFullCheck());

  assert(pTCTetB->getOppositeVert(pFNewMid) == pVVertB);
  assert(checkOrient3D(pVVertC, pVVertD, pVVertE, pVVertB) == -1);
  assert(pTCTetA->getOppositeVert(pFNewMid) == pVVertA);
  assert(checkOrient3D(pVVertC, pVVertD, pVVertE, pVVertA) == 1);
#endif

  int iRetVal = 1;
  if (m_swapRecursively) {
    iRetVal += iFaceSwap_deprecated(pFFaceB);
    iRetVal += iFaceSwap_deprecated(pFFaceC);
    iRetVal += iFaceSwap_deprecated(pFFaceE);
    iRetVal += iFaceSwap_deprecated(pFFaceF);
  }
  return (iRetVal);
}

//@@ Swap four tets for four, in the case where two faces are coplanar.
int
VolMesh::reconfTet44_deprecated(TetCell* apTCTets[], TriFace* const pF,
				Vert* const pVVertA, Vert* const pVVertB,
				Vert* const pVVertC, Vert* const pVVertD,
				Vert* const pVVertE)
{
  // The cells and verts are assumed to have been set up properly already.

  // Here's the before picture:
  /*      	        C
   /\
        	      /||\
        	     / |  \                 A,D are on the base;
   /  | | \                B,E are front, C,F are back.
   /   |    \               Cell A is on the left, B on
   /    |  |  \	            the right.
   /     |      \
  		/      |   |   \
  	       /     C |     F  \
  	      /	       |  E |    \
  	     /   B     | _ - B    \
  	    /  _ _ - - |    /  - _ \
  	   D -         |           -E
   -__	    A  |  /  D  __--
   -__     |     __-
   -___ |/__--
   -A-

   ___---B_
   ____----	    /  --__
   D---            /       -E
   - _	    A     / |D  _ -
   \  - _       /   _ -   /
   \     -_   / _-      /
   \       -A-  |     /
   \     C1|     F1 /
   \      |  |    /
   \ B1  |  E1  /
   \    |     /	            A,D are on the base;
   \   | |  /               B1,E1 are front, C1,F1 are back.
   \  |   /                Cell A1 is on the left, B1 on
   \ || /                 the right.
   \| /                  Face G divides these two cells.
   \/
   C1

   */
  TetCell* pTCTetA = apTCTets[0];
  TetCell* pTCTetB = apTCTets[1];
  assert(pTCTetA->hasVert(pVVertD));
  assert(pTCTetB->hasVert(pVVertE));
  assert(pTCTetA == pF->getLeftCell());
  assert(pTCTetB == pF->getRightCell());

  Face* pFFaceA = pTCTetA->getOppositeFace(pVVertC);
  Face* pFFaceB = pTCTetA->getOppositeFace(pVVertB);
  Face* pFFaceC = pTCTetA->getOppositeFace(pVVertA);

  Face* pFFaceD = pTCTetB->getOppositeFace(pVVertC);
  Face* pFFaceE = pTCTetB->getOppositeFace(pVVertB);
  Face* pFFaceF = pTCTetB->getOppositeFace(pVVertA);

  // Now identify the cells, faces and vert on the opposite side of the plane
  TetCell* pTCTetA1 = dynamic_cast<TetCell*>(pFFaceA->getOppositeCell(pTCTetA));
  TetCell* pTCTetB1 = dynamic_cast<TetCell*>(pFFaceD->getOppositeCell(pTCTetB));

  // Check to make sure that there isn't something fishy going on with
  // region boundaries.
  if ((pTCTetA->getRegion() != pTCTetB->getRegion())
      || (pTCTetA1->getRegion() != pTCTetB1->getRegion()))
    return 0;

  assert(pTCTetA1->hasVert(pVVertD));
  assert(pTCTetB1->hasVert(pVVertE));

  Face* pFFaceB1 = pTCTetA1->getOppositeFace(pVVertB);
  Face* pFFaceC1 = pTCTetA1->getOppositeFace(pVVertA);

  Face* pFFaceE1 = pTCTetB1->getOppositeFace(pVVertB);
  Face* pFFaceF1 = pTCTetB1->getOppositeFace(pVVertA);
  Face* pFFaceG = pTCTetB1->getOppositeFace(pVVertE);

  Vert* pVVertC1 = pTCTetB1->getOppositeVert(pFFaceD);

  if (pF->isLocked() || pFFaceA->isLocked() || pFFaceD->isLocked()
      || pFFaceG->isLocked())
    return 0;

  // Do the actual four-for-four swap here.  The faces all stay in the
  // same places except for pF, the dividing face, which gets turned.
  // Cell A ends up in front, cell B in back.  The only vert hints
  // affected are for verts A and B, which are set explicitly.
  // Below the plane, pFFaceG gets moves, cell A1 ends up in front, cell
  // B1 in back.

  int iReg = pTCTetA->getRegion();
  int iReg1 = pTCTetA1->getRegion();
  // Delete all four tets
  deleteCell(pTCTetA);
  deleteCell(pTCTetB);
  deleteCell(pTCTetA1);
  deleteCell(pTCTetB1);

  // Now create four new tets directly from verts.
  bool qExist;
  pTCTetA = createTetCell(qExist, pVVertD, pVVertA, pVVertE, pVVertC, iReg);
  assert(!qExist);
  pTCTetB = createTetCell(qExist, pVVertE, pVVertB, pVVertD, pVVertC, iReg);
  assert(!qExist);

  pTCTetA1 = createTetCell(qExist, pVVertE, pVVertA, pVVertD, pVVertC1, iReg1);
  assert(!qExist);
  pTCTetB1 = createTetCell(qExist, pVVertD, pVVertB, pVVertE, pVVertC1, iReg1);
  assert(!qExist);

  assert(pTCTetA->doFullCheck());
  assert(pTCTetB->doFullCheck());
  assert(pTCTetA1->doFullCheck());
  assert(pTCTetB1->doFullCheck());

#ifndef NDEBUG
  Face *pFNewMid = findCommonFace(pTCTetA, pTCTetB);
  pFFaceA = findCommonFace(pTCTetA, pTCTetA1);
  pFFaceD = findCommonFace(pTCTetB, pTCTetB1);
  pFFaceG = findCommonFace(pTCTetA1, pTCTetB1);
  assert(pFNewMid->doFullCheck());
  assert(pFFaceA->doFullCheck());
  assert(pFFaceB->doFullCheck());
  assert(pFFaceC->doFullCheck());
  assert(pFFaceD->doFullCheck());
  assert(pFFaceE->doFullCheck());
  assert(pFFaceF->doFullCheck());

  assert(pFFaceB1->doFullCheck());
  assert(pFFaceC1->doFullCheck());
  assert(pFFaceE1->doFullCheck());
  assert(pFFaceF1->doFullCheck());
  assert(pFFaceG->doFullCheck());

  assert(pTCTetB->getOppositeVert(pFNewMid) == pVVertB);
  assert(checkOrient3D(pVVertC, pVVertD, pVVertE, pVVertB) == -1);
  assert(pTCTetA->getOppositeVert(pFNewMid) == pVVertA);
  assert(checkOrient3D(pVVertC, pVVertD, pVVertE, pVVertA) == 1);

  assert(pTCTetB1->getOppositeVert(pFFaceG) == pVVertB);
  assert(checkOrient3D(pVVertC1, pVVertD, pVVertE, pVVertB) == 1);
  assert(pTCTetA1->getOppositeVert(pFFaceG) == pVVertA);
  assert(checkOrient3D(pVVertC1, pVVertD, pVVertE, pVVertA) == -1);
#endif

  int iRetVal = 2;
  if (m_swapRecursively) {
    iRetVal += iFaceSwap_deprecated(pFFaceB);
    iRetVal += iFaceSwap_deprecated(pFFaceC);
    iRetVal += iFaceSwap_deprecated(pFFaceE);
    iRetVal += iFaceSwap_deprecated(pFFaceF);
    iRetVal += iFaceSwap_deprecated(pFFaceB1);
    iRetVal += iFaceSwap_deprecated(pFFaceC1);
    iRetVal += iFaceSwap_deprecated(pFFaceE1);
    iRetVal += iFaceSwap_deprecated(pFFaceF1);
  }
  return (iRetVal);
}

