//@ Include files, #define's, etc.
#include <stdlib.h>
#include <queue>

#include <math.h>
#include "GR_assert.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_BFace.h"
#include "GR_Vertex.h"
#include "GR_Geometry.h"
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include "GR_Vec.h"
#include "GR_Aniso.h"

//@ Swap away iFace if this is legal and improves the mesh quality measure
// Swapping typically propagates through the mesh, and the return value from
// iFaceSwap is the total number of swaps done during this invocation.
int
Mesh2D::iFaceSwap_deprecated(Face*& pF)
{
  /*if(pAnisoRef!=NULL){
   return pAnisoRef->iAniSwap(pF); //Anistropic replacement for iFaceSwap
   }*/
  if (!pF->doFullCheck())
    return (0);
  if (!pF->isSwapAllowed())
    return 0;

  assert(pF->doFullCheck());
  assert(getSwapType() != eNoSwap);

  //@@ Case: one or both cells is not a tri, including boundaries
  Cell *pCLeft = pF->getLeftCell();
  Cell *pCRight = pF->getRightCell();

  if (!pCLeft->isValid() || !pCRight->isValid())
    return (0);
  if (pCLeft->getType() != Cell::eTriCell
      || pCRight->getType() != Cell::eTriCell)
    return (0);

  Vert *pVVertA = pF->getVert(0);
  Vert *pVVertB = pF->getVert(1);
  TriCell *pTCLeft = dynamic_cast<TriCell*>(pCLeft);
  TriCell *pTCRight = dynamic_cast<TriCell*>(pCRight);
  assert(pTCLeft != NULL);
  assert(pTCRight != NULL);

  Vert *pVVertC = pTCLeft->getOppositeVert(pF);
  Vert *pVVertD = pTCRight->getOppositeVert(pF);

  int iOrientA = checkOrient2D(pVVertA, pVVertD, pVVertC);
  int iOrientB = checkOrient2D(pVVertB, pVVertC, pVVertD);
  if (iOrientA != 1 || iOrientB != 1)
    return 0;
  if (doSwap(pVVertA, pVVertB, pVVertC, pVVertD)) {
    return reconfigure(pF);
  }
  else
    return (0);
}

//@ Swap away iFace if this is legal and improves the mesh quality measure
// Swapping typically propagates through the mesh, and the return value from
// iFaceSwap is the total number of swaps done during this invocation.
bool
Mesh2D::removeEdge(Face*& pF)
{
  //printf("Getting rid of an edge\n");
  //  assert(0);
  if (!pF->isValid() || pF->isDeleted() || !pF->isSwapAllowed())
    return (false);

  assert(pF->doFullCheck());

  //@@ Case: one or both cells is not a tri, including boundaries
  Cell *pCL = pF->getLeftCell();
  Cell *pCR = pF->getRightCell();
  if (!pCL->isValid() || !pCR->isValid())
    return (false);
  if (pCL->getType() != Cell::eTriCell || pCR->getType() != Cell::eTriCell)
    return (false);

  TriCell *pTCL = dynamic_cast<TriCell*>(pCL);
  TriCell *pTCR = dynamic_cast<TriCell*>(pCR);
  assert(pTCL != NULL);
  assert(pTCR != NULL);

  Vert *pVVertA = pF->getVert(0);
  Vert *pVVertB = pF->getVert(1);
  Vert *pVVertC = pTCL->getOppositeVert(pF);
  Vert *pVVertD = pTCR->getOppositeVert(pF);

  if (checkOrient2D(pVVertA, pVVertD, pVVertC) != 1
      || checkOrient2D(pVVertB, pVVertC, pVVertD) != 1)
    return false;

  disallowSwapRecursion();
  reconfigure(pF);
  allowSwapRecursion();
  return true;
}

// This tie breaker is used only to pick a unique configuration in the
// event that a tie occurs trying to decide  how to swap.  The only
// thing which this choice of tiebreaker has to recommend it is that it
// gives a unique answer for any two pairs of verts (AB and CD).
static bool
qTieBreak(const Vert * pVVertA, const Vert * pVVertB, const Vert * pVVertC,
	  const Vert * pVVertD)
{
  unsigned long ulA = reinterpret_cast<unsigned long>(pVVertA);
  unsigned long ulB = reinterpret_cast<unsigned long>(pVVertB);
  unsigned long ulC = reinterpret_cast<unsigned long>(pVVertC);
  unsigned long ulD = reinterpret_cast<unsigned long>(pVVertD);
  if (min(ulA, ulB) > min(ulC, ulD))
    return true;
  else
    return false;
}

//@@ Determines whether swapping is needed for local Delaunay-ness.
static bool
qDoSwapDelaunay(const Vert * pVVertA, const Vert * pVVertB,
		const Vert * pVVertC, const Vert * pVVertD)
{
  switch (isInCircle(pVVertA, pVVertB, pVVertC, pVVertD))
    {
    case -1:
      return false;		// pVVertD is outside

    case 1:
      return true;		// pVVertD is inside

    case 0:			// Need some sort of tie-breaker for uniqueness
      return false;		// Tie breaker design is a little tricky....
      //     return qTieBreak(pVVertA, pVVertB, pVVertC, pVVertD);

    default:
      assert2(0, "Illegal return value from iIncircle");
      return false;
    }
}

//@@ Determines whether swapping will improve the maximum face angle.
static bool
qDoSwapMaxDihed(const Vert * pVVertA, const Vert * pVVertB,
		const Vert * pVVertC, const Vert * pVVertD)
{
  // This algorithm finds face angles between adjacent faces by dotting
  // their unit normals.  The largest magnitude loses.
  //
  // To prevent pairs from flopping back and forth, the tie-breaker is
  // invoked if the inequality is small.

  double dEps = 1.e-12;
  // All of these normals are intended to point inwards.
  double adNormDA[2], adNormBD[2], adNormCB[2], adNormAC[2];
  calcNormal2D(pVVertA->getCoords(), pVVertD->getCoords(), adNormDA);
  calcNormal2D(pVVertD->getCoords(), pVVertB->getCoords(), adNormBD);
  calcNormal2D(pVVertB->getCoords(), pVVertC->getCoords(), adNormCB);
  calcNormal2D(pVVertC->getCoords(), pVVertA->getCoords(), adNormAC);

  NORMALIZE2D(adNormDA);
  NORMALIZE2D(adNormBD);
  NORMALIZE2D(adNormCB);
  NORMALIZE2D(adNormAC);

  double dDotBDA = dDOT2D(adNormBD, adNormDA);
  double dDotACB = dDOT2D(adNormAC, adNormCB);
  double dDotDAC = dDOT2D(adNormDA, adNormAC);
  double dDotCBD = dDOT2D(adNormCB, adNormBD);

  double dMaxThis = max(dDotBDA, dDotACB);
  double dMaxOther = max(dDotDAC, dDotCBD);

//   // Use the commented code only to create bad meshes to test swapping.
// #warning Deliberately broken swapping decision
//   if (dMaxThis > 0.97 || dMaxOther > 0.97) {
//     //     // This is the -right- decision
  if (dMaxThis > dMaxOther + dEps)
    return true;
  else if (dMaxThis + dEps < dMaxOther)
    return false;
  else
    return qTieBreak(pVVertA, pVVertB, pVVertC, pVVertD);
//     // Use the following clause only to create bad meshes to test swapping.
//   }
//   else {
//     if (dMaxThis < dMaxOther - dEps)
//       return true;
//     else if (dMaxThis - dEps > dMaxOther)
//       return false;
//     else
//       return qTieBreak (pVVertA, pVVertB, pVVertC, pVVertD);
//   }
}

bool
Mesh2D::doSwap(const Vert* const pVVertA, const Vert* const pVVertB,
	       const Vert* const pVVertL, const Vert* const pVVertR,
	       const Vert* const pVVertE) const
{
  assert(pVVertA->isValid());
  assert(pVVertB->isValid());
  assert(pVVertL->isValid());
  assert(pVVertR->isValid());
  assert(!pVVertE->isValid());

  if (pAnisoRef != NULL) {
    double dCircRad1, dCircRad2, dTemp;
    //   Determine which configuration has the smallest curvature
    dCircRad1 = pAnisoRef->dAnisoCircumRadius(pVVertA, pVVertB, pVVertL);
    //  dCircRad1*=dCircRad1;
    dTemp = pAnisoRef->dAnisoCircumRadius(pVVertA, pVVertB, pVVertR);
    dCircRad1 += dTemp;

    // printf("Circrad1 is %.3e\n",dCircRad1);
    dCircRad2 = pAnisoRef->dAnisoCircumRadius(pVVertL, pVVertR, pVVertA);
    //  dCircRad2*=dCircRad2;
    dTemp = pAnisoRef->dAnisoCircumRadius(pVVertL, pVVertR, pVVertB);
    dCircRad2 += dTemp;

    //   double dCurvOld,dCurvNew;
//     dCurvOld=pAnisoRef->dCurvature(pVVertA,pVVertB);
//     dCurvNew=pAnisoRef->dCurvature(pVVertL,pVVertR);
//     return(dCurvNew<dCurvOld);
    // printf("Circrad2 is %.3e\n",dCircRad2);
    return (dCircRad2 < dCircRad1);
  }

  //first check to see if swapping is possible
  if (checkOrient2D(pVVertB, pVVertL, pVVertA) != 1
      || checkOrient2D(pVVertL, pVVertA, pVVertR) != 1
      || checkOrient2D(pVVertA, pVVertR, pVVertB) != 1
      || checkOrient2D(pVVertR, pVVertB, pVVertL) != 1)
    return false;

  switch (getSwapType())
    {
    case eMinSine:
    case eDelaunay:
      return qDoSwapDelaunay(pVVertA, pVVertB, pVVertL, pVVertR);
    case eMaxDihed:
      return qDoSwapMaxDihed(pVVertA, pVVertB, pVVertL, pVVertR);
    case eRandom:
      return (drand48() > 0.5);
    case eFromQual:
      {
	assert(m_qual != NULL);
	TriCellCV TCCV11(const_cast<Vert*>(pVVertA), const_cast<Vert*>(pVVertB),
			 const_cast<Vert*>(pVVertL));
	TriCellCV TCCV12(const_cast<Vert*>(pVVertB), const_cast<Vert*>(pVVertA),
			 const_cast<Vert*>(pVVertR));
	double dQual11 = m_qual->dEvaluate(&TCCV11);
	double dQual12 = m_qual->dEvaluate(&TCCV12);

	TriCellCV TCCV21(const_cast<Vert*>(pVVertL), const_cast<Vert*>(pVVertR),
			 const_cast<Vert*>(pVVertB));
	TriCellCV TCCV22(const_cast<Vert*>(pVVertR), const_cast<Vert*>(pVVertL),
			 const_cast<Vert*>(pVVertA));
	double dQual21 = m_qual->dEvaluate(&TCCV21);
	double dQual22 = m_qual->dEvaluate(&TCCV22);
	bool qMaximize = m_qual->qMaximizeMinValue();
	if (qMaximize) {
	  double dQual1 = min(dQual11, dQual12);
	  double dQual2 = min(dQual21, dQual22);
	  if (dQual1 > dQual2)
	    return true;
	  else
	    return false;
	}
	else {
	  double dQual1 = max(dQual11, dQual12);
	  double dQual2 = max(dQual21, dQual22);
	  if (dQual1 < dQual2)
	    return true;
	  else
	    return false;
	}
      }
    default:
      return (false);
    }
}

int
Mesh2D::reconfigure(Face *& pFEdge)
{
  if (!pFEdge->isSwapAllowed())
    return 0;

  assert(pFEdge->doFullCheck());
  assert(!((m_swapMeasure == eNoSwap) && m_swapRecursively));
  if (pFEdge->getFaceLoc() == Face::eBdryTwoSide)
    return (0);

  // Identify all the parts
  TriCell *pTri0 = dynamic_cast<TriCell*>(pFEdge->getLeftCell());
  TriCell *pTri1 = dynamic_cast<TriCell*>(pFEdge->getRightCell());

  if (pCellQueue != NULL) {
    pCellQueue->vRemoveCell(pFEdge->getCell(0));
    pCellQueue->vRemoveCell(pFEdge->getCell(1));
  }
  Vert *pVA = pFEdge->getVert(0);
  Vert *pVB = pFEdge->getVert(1);
  Vert *pV0 = pTri0->getOppositeVert(pFEdge);
  Vert *pV1 = pTri1->getOppositeVert(pFEdge);

  Face *pFA = pTri0->getOppositeFace(pVB);
  Face *pFB = pTri0->getOppositeFace(pVA);
  Face *pFC = pTri1->getOppositeFace(pVA);
  Face *pFD = pTri1->getOppositeFace(pVB);
  int iReg = pTri0->getRegion();
  deleteCell(pTri0);
  deleteCell(pTri1);

  bool autoRotate = true;
  Cell *pC0 = createTriCell(pV1, pV0, pVA, iReg, autoRotate);
  Cell *pC1 = createTriCell(pV0, pV1, pVB, iReg, autoRotate);

  // Identify the edge because surface reconnection requires that
  // the new edge between the cells be accessible to it.
  pFEdge = findCommonFace(pC0, pC1);

#ifndef OMIT_VERTEX_HINTS
  pVA->setHintFace(pFA);
  pVB->setHintFace(pFB);
#endif
  assert(pC0->doFullCheck());
  assert(pC1->doFullCheck());

  assert(pFA->doFullCheck());
  assert(pFB->doFullCheck());
  assert(pFC->doFullCheck());
  assert(pFD->doFullCheck());

  // if(pAnisoRef!=NULL){
  if (pCellQueue != NULL) {
    pCellQueue->vAddCell(pC0);
    pCellQueue->vAddCell(pC1);
  }
  //  }
  int iRetVal = 1;
  if (m_swapRecursively) {
    if (pAnisoRef == NULL) {
      iRetVal += iFaceSwap_deprecated(pFA);
      iRetVal += iFaceSwap_deprecated(pFB);
      iRetVal += iFaceSwap_deprecated(pFC);
      iRetVal += iFaceSwap_deprecated(pFD);
    }
    else {
      //pAnisoRef->vSmoothMetric(pFEdge->getVert(0));
      //pAnisoRef->vSmoothMetric(pFEdge->getVert(1));
      iRetVal += pAnisoRef->iAniSwap(pFA);
      iRetVal += pAnisoRef->iAniSwap(pFB);
      iRetVal += pAnisoRef->iAniSwap(pFC);
      iRetVal += pAnisoRef->iAniSwap(pFD);
    }
  }

  return iRetVal;
}

//@ Make certain that a particular edge is included in a 2D mesh.
bool
Mesh2D::recoverEdge(const int iV0, const int iV1)
// Force swaps of all edges that are crossed by the edge between
// the verts indicated by iV0 and iV1.
{
  assert(iV0 != iV1);
  Vert *pV0 = getVert(iV0);
  Vert *pV1 = getVert(iV1);
  Face *pFDummy;
  return recoverEdge(pV0, pV1, pFDummy);
}

bool
Mesh2D::recoverEdge(Vert* const pV0, Vert* const pV1)
{
  Face *pFTmp;
  return recoverEdge(pV0, pV1, pFTmp);
}

bool
Mesh2D::recoverEdge(Vert* const pV0, Vert* const pV1, Face*& pFRet)
// Force swaps of all edges that are crossed by the edge between
// the verts indicated by pV0 and pV1.
{
  assert(pV0->isValid());
  assert(pV1->isValid());
  bool qOldRecursion = isSwapRecursionAllowed();
  disallowSwapRecursion();
  pFRet = pFInvalidFace;

//  vMessage(3, "Inside qRecoverEdge: %d, %d : (%f, %f) - (%f, %f)\n",
//	   iV0, iV1, pV0->dX (), pV0->dY (), pV1->dX (), pV1->dY ());

  // Find the cell across which to begin walking from pV0 to pV1.
  std::set<Cell*> spCInc;
  std::set<Vert*> spVNeigh;
  findNeighborhoodInfo(pV0, spCInc, spVNeigh);

  // For each cell, find out which, if any, of its non-incident faces
  // are crossed by pV0, pV1.  In the process, identify the face
  // connecting pV0 to pV1 if it exists.
  std::set<Cell*>::iterator iter;
  Face *pFToCross = pFInvalidFace;
  Cell *pC = pCInvalidCell;
  bool qGotIt = false;
  for (iter = spCInc.begin(); !qGotIt && (iter != spCInc.end()); ++iter) {
    Cell *pCCand = *iter;
    assert(pCCand->hasVert(pV0));
    for (int iF = pCCand->getNumFaces() - 1; iF >= 0; iF--) {
      Face *pFCand = pCCand->getFace(iF);
      if (pFCand->hasVert(pV0)) {
	if (pFCand->hasVert(pV1)) {
	  pFRet = pFCand;
	  return true;
	}
	// If a face incident on pV0 isn't also incident on pV1, then
	// don't bother.
      }

      else if ((checkOrient2D(pV0, pV1, pFCand->getVert(0))
	  == checkOrient2D(pV1, pV0, pFCand->getVert(1)))
	  && (checkOrient2D(pFCand->getVert(0), pFCand->getVert(1), pV0)
	      == checkOrient2D(pFCand->getVert(1), pFCand->getVert(0), pV1))) {
	// This edge crosses the segment pV0, pV1; it's the place to
	// start.
	pFToCross = pFCand;
	pC = pCCand;
	qGotIt = true;
	break;

      }
    }
  }
  assert(spVNeigh.find(pV1) == spVNeigh.end());
#ifndef NDEBUG
  if (pFToCross == pFInvalidFace) {
    pV0->printVertInfo();
    pV1->printVertInfo();
  }
#endif
  assert(pFToCross != pFInvalidFace);
  assert(pC != pCInvalidCell);
  if (pC->getType() != Cell::eTriCell) {
    // Needed to swap an edge that bounds a quad.
    return false;
  }

  // Walk from pV0 to pV1, listing all the edges that cross (pV0, pV1).
  // The walk should always succeed, hitting no vertices square on and
  // running into no boundaries.
  std::queue<Face*> qpFBadEdges;
  Face *pF = pFToCross;
  Vert *pV = pVInvalidVert;
  Vert *pVLeft, *pVRight;
  if (checkOrient2D(pV0, pV1, pF->getVert(0)) == 1) {
    pVLeft = pF->getVert(0);
    pVRight = pF->getVert(1);
  }
  else {
    pVLeft = pF->getVert(1);
    pVRight = pF->getVert(0);
  }
  assert(checkOrient2D(pV0, pV1, pVLeft) == 1);
  assert(checkOrient2D(pV1, pV0, pVRight) == 1);
  assert(checkOrient2D(pVLeft, pVRight, pV1) == 1);
  assert(checkOrient2D(pVRight, pVLeft, pV0) == 1);

  do {
    // Add the current to the list of those which need removing
    qpFBadEdges.push(pF);
    // Grab the next cell and vertex
    pC = pF->getOppositeCell(pC);
    if (pC->getType() != Cell::eTriCell) {
      return false;
    }

    // This cast will succeed; all other cases just exited.
    TriCell* pTC = dynamic_cast<TriCell*>(pC);
    assert(pTC != NULL);
    pV = pTC->getOppositeVert(pF);
    int iOrient = checkOrient2D(pV, pV0, pV1);
    switch (iOrient)
      {
      case 1:
	pF = pTC->getOppositeFace(pVLeft);
	pVLeft = pV;
	break;			// end left side case

      case -1:
	pF = pTC->getOppositeFace(pVRight);
	pVRight = pV;
	break;			// end right side case

      default:
	// This assertion fails if a vert lies exactly on the segment
	// connecting pV0 and pV1.  This is an error in the existing
	// triangulation, or in the edge recovery request.
	assert2(pV == pV1, "Ties should be impossible in forced edge swapping");
	break;
      }				// end switch

  }
  while (pV != pV1);		// end walk across mesh

  // Attempt to swap each edge in turn.  If the replacement would cross
  // the edge to be recovered, don't bother.  Also don't swap in convex
  // quad cases.  Eventually, the list of bad edges will be empty.
  while (!qpFBadEdges.empty()) {
    pF = qpFBadEdges.front();
    qpFBadEdges.pop();

    Cell *pCL = pF->getLeftCell();
    Cell *pCR = pF->getRightCell();
    // These assertions should be okay, even for mixed-element meshes,
    // because all pipes with quads in them should have returned before
    // this.
    assert(pCL->getType() == Cell::eTriCell);
    assert(pCR->getType() == Cell::eTriCell);
    TriCell* pTCL = dynamic_cast<TriCell*>(pCL);
    assert(pTCL != NULL);
    TriCell* pTCR = dynamic_cast<TriCell*>(pCR);
    assert(pTCR != NULL);

    Vert *pVA = pF->getVert(0);
    Vert *pVB = pF->getVert(1);
    Vert *pVC = pTCL->getOppositeVert(pF);
    Vert *pVD = pTCR->getOppositeVert(pF);
    assert(checkOrient2D(pVA, pVB, pVC) == 1);
    assert(checkOrient2D(pVA, pVB, pVD) == -1);
    assert(pF->getType() == Face::eEdgeFace);
    if (checkOrient2D(pVA, pVD, pVC) == 1
	&& checkOrient2D(pVB, pVC, pVD) == 1) {
      reconfigure(pF);	// No swap recursion
      // If the new diagonal crosses the line connecting pV0 and pV1,
      // pF still needs to be forced out, but it needs a different
      // hull around it.  In this case, leave pF in the list and  go
      // on to the next edge.
      // This test passes if pVC or pVD is identical to pV0 or pV1, or
      // if pVC and pVD are on the same side of pV0-pV1.

      if (checkOrient2D(pV0, pV1, pVC) * checkOrient2D(pV0, pV1, pVD) == -1) {
	qpFBadEdges.push(pF);
      }
    }
    else {
      // Re-queue this one so that you can try it again when it's turn
      // comes up.
      qpFBadEdges.push(pF);
    }
  }
  if (qOldRecursion)
    allowSwapRecursion();
  pFRet = pF;
  return true;
}

