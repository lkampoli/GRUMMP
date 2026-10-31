#include "GR_BFace3D.h"
#include "GR_Cell.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Face.h"
#include "GR_Geometry.h"
#include "GR_Util.h"

namespace GRUMMP
{
  FaceSwapInfo3D::FaceSwapInfo3D(Face* const face, Vert* const north,
				 Vert* const south, Vert* const other) :
      m_maxAngle(0), m_pF(face), m_fc(eN32), m_nTets(0)
  {
    m_pVerts[0] = north;
    m_pVerts[1] = south;
    m_pVerts[2] = other;
  }

// Notes
//
// There are several places where categorizeFace pre-checks tet
// orientation.  These tests should always pass for untangled meshes, of
// course.  The assumption built into the swapping code --- but never
// checked for up front --- is that the meshes it's handed will be
// untangled in advance.

  void
  FaceSwapInfo3D::getTets(int &nTets, Cell* pTets[4]) const
  {
    nTets = m_nTets;
    int i;
    for (i = 0; i < m_nTets; i++) {
      pTets[i] = m_pTets[i];
    }
    for (; i < 4; i++) {
      pTets[i] = NULL;
    }
  }

//@ Decide which type of triangular face this is, anyway.
  FaceSwapInfo3D::faceCat
  FaceSwapInfo3D::categorizeFace(Face* pF)
  {
    VALIDATE_INPUT(pF->isValid() && !pF->isDeleted());
    m_pF = pF;
    m_nTets = 0;
    m_fc = eOther;
    m_pVerts[0] = m_pVerts[1] = m_pVerts[2] = m_pVerts[3] = m_pVerts[4] =
	m_pVerts[5] = NULL;
    m_pBFaces[0] = m_pBFaces[1] = NULL;

    // Don't do anything silly with quads or locked faces ...
    if (pF->getType() != Face::eTriFace || pF->isLocked())
      return m_fc = eOther;
    // ... or boundaries
    if (pF->isBdryFace())
      return m_fc = (eBdry);

    //@@ Case: one or both cells is not a tet, including boundaries
    Cell* pCLeft = pF->getLeftCell();
    Cell* pCRight = pF->getRightCell();
    if (!pCLeft->isValid() || !pCRight->isValid())
      return (eOther);

    if (pCLeft->getType() != Cell::eTet || pCRight->getType() != Cell::eTet
	|| pCLeft->getRegion() != pCRight->getRegion())
      return (eOther);
    int iReg = pCLeft->getRegion();

    //@@ List all tets formed from five verts defining two cells sharing iFace
    m_pTets[0] = dynamic_cast<TetCell*>(pCLeft);
    m_pTets[1] = dynamic_cast<TetCell*>(pCRight);
    m_nTets = 2;
    assert(m_pTets[0] && m_pTets[1]);

    m_pVerts[0] = pF->getVert(0);
    m_pVerts[1] = pF->getVert(1);
    m_pVerts[2] = pF->getVert(2);
    m_pVerts[3] = m_pTets[0]->getOppositeVert(pF);
    m_pVerts[4] = m_pTets[1]->getOppositeVert(pF);

    // Any add'l cell must be a neighbor of both of the first two.
    int i;
    for (i = 0; i < 4 && m_nTets < 4; ++i) {
      Face* pFCandFace = m_pTets[0]->getFace(i);
      Cell* pCCandTet = pFCandFace->getOppositeCell(m_pTets[0]);
      // Don't bother with cells that aren't tets; these are unswappable
      // anyway
      if (pCCandTet->getType() != Cell::eTet)
	continue;
      if (pCCandTet != m_pTets[1] && pCCandTet->hasVert(m_pVerts[4])) {
	// If this face internal to the configuration is locked, then forget it.
	if (pFCandFace->isLocked() || pCCandTet->getRegion() != iReg)
	  return m_fc = (eOther);
	m_pTets[m_nTets++] = dynamic_cast<TetCell*>(pCCandTet);
      }
    }

    assert(m_nTets <= 4);

    switch (m_nTets)
      {
      case 2:
	{
	  //@@ Two tets
#ifndef NDEBUG
	  int check1 = checkOrient3D(m_pVerts[0], m_pVerts[2], m_pVerts[1],
				     m_pVerts[3]);
	  int check2 = checkOrient3D(m_pVerts[0], m_pVerts[1], m_pVerts[2],
				     m_pVerts[4]);
	  assert(check1 >= 0);
	  assert(check2 >= 0);
	  assert(check1 + check2 > 0); // At least one has to be positively oriented.
#endif

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

	  iOrientA = checkOrient3D(m_pVerts[1], m_pVerts[2], m_pVerts[3],
				   m_pVerts[4]);
	  if (iOrientA == -1) {
	    Vert *pVTmp = m_pVerts[0];
	    m_pVerts[0] = m_pVerts[1];
	    m_pVerts[1] = m_pVerts[2];
	    m_pVerts[2] = pVTmp;
	    qBadOrientation = true;
	  }
	  else {
	    if (iOrientA == 1) {
	      // Check whether it's nearly flat
	      double adNorm0[3], adNorm1[3];
	      calcNormal3D(m_pVerts[1]->getCoords(), m_pVerts[2]->getCoords(),
			   m_pVerts[3]->getCoords(), adNorm0);
	      calcNormal3D(m_pVerts[2]->getCoords(), m_pVerts[1]->getCoords(),
			   m_pVerts[4]->getCoords(), adNorm1);
	      NORMALIZE3D(adNorm0);
	      NORMALIZE3D(adNorm1);
	      double dDot = dDOT3D(adNorm0, adNorm1);
	      if (dDot > cos(m_maxAngle * M_PI / 180)) {
		qApproxFlatA = true;
	      }
	    }
	    iOrientB = checkOrient3D(m_pVerts[2], m_pVerts[0], m_pVerts[3],
				     m_pVerts[4]);
	    if (iOrientB == -1) {
	      Vert *pVTmp = m_pVerts[2];
	      m_pVerts[2] = m_pVerts[1];
	      m_pVerts[1] = m_pVerts[0];
	      m_pVerts[0] = pVTmp;
	      qBadOrientation = true;
	    }
	    else if (iOrientB == 1) {
	      // Check whether it's nearly flat
	      double adNorm0[3], adNorm1[3];
	      calcNormal3D(m_pVerts[2]->getCoords(), m_pVerts[0]->getCoords(),
			   m_pVerts[3]->getCoords(), adNorm0);
	      calcNormal3D(m_pVerts[0]->getCoords(), m_pVerts[2]->getCoords(),
			   m_pVerts[4]->getCoords(), adNorm1);
	      NORMALIZE3D(adNorm0);
	      NORMALIZE3D(adNorm1);
	      double dDot = dDOT3D(adNorm0, adNorm1);
	      if (dDot > cos(m_maxAngle * M_PI / 180)) {
		qApproxFlatB = true;
	      }
	    }
	    if (iOrientA + iOrientB == 0)
	      return m_fc = eN20;
	    if (!qBadOrientation) {
	      iOrientC = checkOrient3D(m_pVerts[0], m_pVerts[1], m_pVerts[3],
				       m_pVerts[4]);
	      if (iOrientC == -1) {
		// Already in this order
		// m_pVerts[0] = m_pVerts[0];
		// m_pVerts[1] = m_pVerts[1];
		// m_pVerts[2] = m_pVerts[2];
		qBadOrientation = true;
	      }
	      else if (iOrientC == 1) {
		// Check whether it's nearly flat
		double adNorm0[3], adNorm1[3];
		calcNormal3D(m_pVerts[0]->getCoords(), m_pVerts[1]->getCoords(),
			     m_pVerts[3]->getCoords(), adNorm0);
		calcNormal3D(m_pVerts[1]->getCoords(), m_pVerts[0]->getCoords(),
			     m_pVerts[4]->getCoords(), adNorm1);
		NORMALIZE3D(adNorm0);
		NORMALIZE3D(adNorm1);
		double dDot = dDOT3D(adNorm0, adNorm1);
		if (dDot > cos(m_maxAngle * M_PI / 180)) {
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
	    Face *pF0 = m_pTets[0]->getOppositeFace(m_pVerts[2]);
	    Face *pF1 = m_pTets[1]->getOppositeFace(m_pVerts[2]);

	    if (pF0->isBdryFace() || pF1->isBdryFace())
	      return m_fc = eBdryReflex;
	    else
	      return m_fc = eN32;
	  }

	  switch (iOrientA + iOrientB + iOrientC)
	    {
	    case 0: // Impossible, but included for completeness.
	      return m_fc = eOther;
	    case 1: // Two orientations are zero (three points co-linear).
		    // Hopelessly unswappable.  Bail out.
	      return m_fc = eN20;
	    case 2:
	      {
		//@@@ Four points are coplanar; T22, T44, N44
		// One orientation must be 0; verts are re-labeled to make it
		// iOrientC.  This implies that AB (01) is the coplanar edge.
		assert(!(iOrientA && iOrientB && iOrientC));
		if (iOrientA == 0) {
		  Vert* pVTemp;
		  pVTemp = m_pVerts[0];
		  m_pVerts[0] = m_pVerts[1];
		  m_pVerts[1] = m_pVerts[2];
		  m_pVerts[2] = pVTemp;
		}
		else if (iOrientB == 0) {
		  Vert* pVTemp;
		  pVTemp = m_pVerts[0];
		  m_pVerts[0] = m_pVerts[2];
		  m_pVerts[2] = m_pVerts[1];
		  m_pVerts[1] = pVTemp;
		}
		else if (iOrientC == 0) {
		}
		// The following assertion has been removed deliberately, to
		// accomodate the near-flat cases tested for with dot products
		// above.  Prior to adding support for near-flat swapping, the
		// orientation was in fact re-aligned correctly, so there's no
		// cause for concern.
		// Verify that the relabeling was done correctly
		//	assert( checkOrient3D(m_pVerts[0], m_pVerts[1], m_pVerts[3], m_pVerts[4]) ==
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
		TetCell *pTCTetA = m_pTets[0];
		TetCell *pTCTetB = m_pTets[1];
		// These are the faces that are coplanar
		Face* pFFaceA = pTCTetA->getOppositeFace(m_pVerts[2]);
		Face* pFFaceB = pTCTetB->getOppositeFace(m_pVerts[2]);
		// Bail if either of these faces are locked.
		if (pFFaceA->isLocked() || pFFaceB->isLocked())
		  return m_fc = eOther;

		Cell* pCOppCellA = pFFaceA->getOppositeCell(pTCTetA);
		Cell* pCOppCellB = pFFaceB->getOppositeCell(pTCTetB);

		if (pCOppCellA->getType() == Cell::eTriBFace
		    && pCOppCellB->getType() == Cell::eTriBFace) {
		  // Both faces are on the boundary
		  TriBFace* pTBFA = dynamic_cast<TriBFace*>(pCOppCellA);
		  TriBFace* pTBFB = dynamic_cast<TriBFace*>(pCOppCellB);
		  // These are already -known- to be co-planar, so checking bdry
		  // patches is sufficient.
		  // FIX ME:  This is a spot where checking patches, not
		  // BC's, is the right thing to do.  This needs a careful
		  // look in the context of the curved bdry meshing stuff,
		  // though.
		  // 	  if (pTBFA->pPatchPointer() == pTBFB->pPatchPointer())
		  if (pTBFA->getBdryCondition() == pTBFB->getBdryCondition()) {
		    m_pBFaces[0] = pTBFA;
		    m_pBFaces[1] = pTBFB;
		    return m_fc = eT22; // case 1
		  }
		  else
		    return m_fc = eOther;  // case 2
		}
		else if (pCOppCellA->getType() == Cell::eTet
		    && pCOppCellB->getType() == Cell::eTet) {
		  // Both faces are internal to the mesh and separate pairs of tets

		  // Find out if they have a common face
		  Face *pFCommon = findCommonFace(pCOppCellA, pCOppCellB);

		  if (pFCommon->isValid()) {
		    if (!pFCommon->isLocked()) {
		      m_pTets[2] = dynamic_cast<TetCell*>(pCOppCellA);
		      m_pTets[3] = dynamic_cast<TetCell*>(pCOppCellB);
		      m_nTets = 4;
		      m_pVerts[5] = m_pTets[2]->getOppositeVert(pFFaceA);
		      return m_fc = eT44; // case 4
		    }
		    else {
		      return m_fc = eOther;
		    }
		  }
		  else {
		    // Didn't hit a common face; this face is N44
		    return m_fc = eN44; // case 5
		  }
		}
		else
		  // Exactly one face on the boundary or internal faces with
		  // cells other than tets
		  return m_fc = eOther; // cases 3,6
	      } // End two-tet cases with four verts coplanar
	    case 3:
	      {
		// This is probably a T23 case, but it might be a perturbed T22
		// case, which should be tagged so that removal of surface verts
		// will work better.  (Actually, it could even be a perturbed
		// N20, which is a disaster during removal if you try to run
		// with it as a T23! so nuke this case.)
		//
		// FIX ME  Surely there must be some way to factor the
		// following clause an the near-identical one above to get
		// rid of the repetition.
		if ((qApproxFlatA && qApproxFlatB)
		    || (qApproxFlatB && qApproxFlatC)
		    || (qApproxFlatC && qApproxFlatA))
		  return m_fc = eN20;
		if (qApproxFlatA || qApproxFlatB || qApproxFlatC) {
		  // It is a perturbed T22
		  if (qApproxFlatA) {
		    Vert* pVTemp;
		    pVTemp = m_pVerts[0];
		    m_pVerts[0] = m_pVerts[1];
		    m_pVerts[1] = m_pVerts[2];
		    m_pVerts[2] = pVTemp;
		  }
		  else if (qApproxFlatB) {
		    Vert* pVTemp;
		    pVTemp = m_pVerts[0];
		    m_pVerts[0] = m_pVerts[2];
		    m_pVerts[2] = m_pVerts[1];
		    m_pVerts[1] = pVTemp;
		  }
		  else if (qApproxFlatC) {
		  }

		  TetCell *pTCTetA = m_pTets[0];
		  TetCell *pTCTetB = m_pTets[1];
		  // These are the faces that are coplanar
		  Face* pFFaceA = pTCTetA->getOppositeFace(m_pVerts[2]);
		  Face* pFFaceB = pTCTetB->getOppositeFace(m_pVerts[2]);
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
		    if (pTBFA->getBdryCondition()
			== pTBFB->getBdryCondition()) {
		      m_pBFaces[0] = pTBFA;
		      m_pBFaces[1] = pTBFB;
		      return m_fc = eT22; // case 1
		    }
		    else
		      return m_fc = eOther;  // case 2
		  }
		} // If one is nearly flat
		  // Otherwise, the configuration is convex; swappable two for three.
		return m_fc = eT23;
	      }
	    default:
	      // No other cases should be possible
	      assert2(0, "Impossible convexity case");
	      return m_fc = eOther;
	    }
	} // End of two tet cases
	assert(0); // Should never get here.
	break;
      case 3:
	{
	  //@@ Three tet cases
	  int iOrient1, iOrient2;
	  m_pVerts[0] = m_pVerts[4];
	  m_pVerts[1] = m_pVerts[3];
	  m_pVerts[2] = pF->getVert(0);
	  m_pVerts[3] = pF->getVert(1);
	  m_pVerts[4] = pF->getVert(2);
	  if (!m_pTets[2]->hasVert(m_pVerts[2])) {
	  }
	  else if (!m_pTets[2]->hasVert(m_pVerts[4])) {
	    Vert* pVTemp = m_pVerts[2];
	    m_pVerts[2] = m_pVerts[4];
	    m_pVerts[4] = m_pVerts[3];
	    m_pVerts[3] = pVTemp;
	  }
	  else {
	    Vert* pVTemp = m_pVerts[2];
	    m_pVerts[2] = m_pVerts[3];
	    m_pVerts[3] = m_pVerts[4];
	    m_pVerts[4] = pVTemp;
	  }

	  // The first of these returns will save a few orientation
	  // primitive evaluations
	  iOrient1 = checkOrient3D(m_pVerts[0], m_pVerts[2], m_pVerts[1],
				   m_pVerts[3]);
	  if (iOrient1 == -1)
	    return m_fc = eOther;
	  else if (iOrient1 == 0)
	    return m_fc = eN30;
	  iOrient2 = checkOrient3D(m_pVerts[0], m_pVerts[1], m_pVerts[2],
				   m_pVerts[4]);
	  if (iOrient2 == -1)
	    return m_fc = eOther;
	  else if (iOrient2 == 0)
	    return m_fc = eN30;
	  assert(iOrient1 == 1);
	  assert(iOrient2 == 1);
	  //@@@@ T32
	  assert(!m_pTets[0]->hasVert(m_pVerts[0]));
	  assert(!m_pTets[1]->hasVert(m_pVerts[1]));
	  assert(!m_pTets[2]->hasVert(m_pVerts[2]));
	  Face *pF01 = findCommonFace(m_pTets[0], m_pTets[1]);
	  Face *pF12 = findCommonFace(m_pTets[1], m_pTets[2]);
	  Face *pF20 = findCommonFace(m_pTets[2], m_pTets[0]);
	  // These have to exist, or we wouldn't be here.
	  if (!pF01 || !pF12 || !pF20) {
	    // We hit the jackpot: we're wrapped around a baffle interior to
	    // the mesh.
	    return m_fc = eOther;
	  }
	  assert(pF01 && pF12 && pF20);
	  if (pF01->isLocked() || pF12->isLocked() || pF20->isLocked())
	    return m_fc = eOther;
	  else
	    return m_fc = eT32;
	} // End of three tet cases
      case 4:
	//@@@ Four tets => N40
	return m_fc = eN40;
	// End of four tet cases
      } // End of classification/weeding
	// Should never reach here, but just in case...
    return m_fc = eOther;
  } // FaceSwapInfo3D::categorizeFace

}
