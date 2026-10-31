/*
 * BFace3D.cxx
 *
 *  Created on: 2013-04-21
 *      Author: cfog
 */

#include "GR_BFace3D.h"
#include "GR_Geometry.h"

#include "CubitVector.hpp"
#include "FacetEvalTool.hpp"
#include "FacetSurface.hpp"
#include "RefFace.hpp"
#include "Surface.hpp"

bool
TriBFaceBase::isCoplanarWith(const double pointLoc[]) const
{
  return (checkOrient3D(getFace(), pointLoc) == 0);
}

enum eEncroachResult
TriBFaceBase::isEncroached(const enum eEncroachType eET) const
// In a constrained Delaunay tetrahedralization, the point opposite
// from a bdry face must encroach on it if any point is going to, so
// check this point only.
{

  // New version of the function that takes into account both
  // interior and regular boundary faces. Modified by SG 06/2008.

  TetCell* pC;
  const Face *pF;
  Vert *pV;

  switch (getNumFaces())
    { //Must consider two cases: interior and regular BFace.
    case 2:
      {
	pF = getFace(1);
	pC = dynamic_cast<TetCell*>(pF->getOppositeCell(this));
	assert(pC);
	pV = pC->getOppositeVert(pF);
	eEncroachResult enc_result = isPointEncroaching(pV->getCoords(), eET,
							false);
	if (enc_result != eClean)
	  return enc_result;
      }
      // falls through
      // no break
    case 1:
      {
	pF = getFace(0);
	pC = dynamic_cast<TetCell*>(pF->getOppositeCell(this));
	assert(pC);
	pV = pC->getOppositeVert(pF);
	return isPointEncroaching(pV->getCoords(), eET, false);
      }
      break;
    default:
      vFatalError("Triangular boundary faces must have one or two faces.",
		  "TriBFaceBase::eIsEncroached(const eEncroachType)");
      break;
    }
  // Should never get here.
  return eClean;

  // Old version commented out but kept here for reference.

//   const Face *pF = getFace();
//   Cell *pC = pF->pCCellOpposite(this);
//   assert(pC->eType() == Cell::eTet);
//   Vert *pVOpp = dynamic_cast<TetCell*>(pC)->pVVertOpposite(pF);
//   return eIsPointEncroaching(pVOpp->adCoords(), eET, false);

}

enum eEncroachResult
TriBFaceBase::isPointEncroachingBall(const double adPoint[3],
				     const bool qTieBreaker) const
				     // This is the easy version of encroachment: a bdry face is encroached
				     // iff a point lies within its circumsphere.
{
  const Face *pF = getFace();
  double adCircCent[3];
  dynamic_cast<const TriFace*>(pF)->calcCircumcenter(adCircCent);
  double dRadius = dDIST3D(adCircCent, pF->getVert(0)->getCoords());
  double dDistance = dDIST3D(adCircCent, adPoint);
  switch (iFuzzyComp(dDistance, dRadius))
    {
    case -1:
      // Distance smaller than radius -> inside
      return eSplit;
    case 1:
      // Distance larger than radius -> outside
      return eClean;
    case 0:
      return qTieBreaker ? eSplit : eClean;
    default:
      assert(0);
      return eClean;
    }
}

void
TriBFaceBase::calcSplitPoint(double split_pt[3]) const
{
  calcCircumcenter(split_pt);
  if (m_bdryPatch)
    return; //flat surface, nothing more to do.

  else if (m_surface) {
//#ifndef NDEBUG
//    Surface* this_surf= m_surface->get_surface_ptr();
//    FacetSurface* this_facet_surf = dynamic_cast<FacetSurface*>(this_surf);
//
//    if(this_facet_surf) {
//      if(this_facet_surf->get_eval_tool()->compare_tol() > 1.e-12) {
//        printf("compare tol = %e\n", this_facet_surf->get_eval_tool()->compare_tol());
//        assert(0);
//      }
//    }
//    else
//      vFatalError("Unsupported geometry type",
//                  "TriBFaceBase::compute_split_point");
//#endif
    CubitVector point(split_pt), closest_point;
    m_surface->find_closest_point_trimmed(point, closest_point);
    closest_point.get_xyz(split_pt);
  }

  else
    return;

}

void
TriBFaceBase::findClosestPointOnBFace(const double pt[3], double close[3])
{

  assert(getNumVerts() == 3);
  closestOnTriangle(pt, getVert(0)->getCoords(), getVert(1)->getCoords(),
		    getVert(2)->getCoords(), close);

  if (m_bdryPatch)
    return; //flat surface, nothing more to do.

  else if (m_surface) {

#ifndef NDEBUG
    Surface* this_surf = m_surface->get_surface_ptr();
    FacetSurface* this_facet_surf = dynamic_cast<FacetSurface*>(this_surf);
    if (this_facet_surf) {
      if (this_facet_surf->get_eval_tool()->compare_tol() > 1.e-12) {
	printf("compare tol = %e\n",
	       this_facet_surf->get_eval_tool()->compare_tol());
	assert(0);
      }
    }
    else
      vFatalError("Unsupported geometry type",
		  "TriBFaceBase::compute_split_point");
#endif

    CubitVector point(close), closest_point;
    m_surface->find_closest_point_trimmed(point, closest_point);
    closest_point.get_xyz(close);

  }
  else {
    // must be flat
    return;
  }

//    else vFatalError("No geometry attached to this patch",
//                     "TriBFaceBase::closest_on_bface");

}

static enum eEncroachResult
eIsPointEncroachingGeneralLens(const double adPoint[3], const bool qTieBreaker,
			       const TriBFaceBase * const pBTF,
			       const enum eEncroachType eET)
// This version of encroachment is more complicated, because we need
// to find the centers and radii of the spheres that make up the lens.
// Then the point encroaches iff it's inside both spheres.  Also, if
// the encroachment type eET == eNewLens, check to see whether the
// point happens to hit the jackpot and requires an offset point to be
// inserted instead of merely splitting the edge.
{
  const Face *pF = pBTF->getFace();
  double adCircCent[3];
  dynamic_cast<const TriFace*>(pF)->calcCircumcenter(adCircCent);
  double dFaceRadius = dDIST3D(adCircCent, pF->getVert(0)->getCoords());
  double dSphereRadius = dFaceRadius * 2 / sqrt(3.);
  double dEncDist = dSphereRadius / 2.;

  double adNorm[3], adCent1[3], adCent2[3];
  calcUnitNormal(pF->getVert(0)->getCoords(), pF->getVert(1)->getCoords(),
		 pF->getVert(2)->getCoords(), adNorm);
  adCent1[0] = adCircCent[0] + adNorm[0] * dEncDist;
  adCent1[1] = adCircCent[1] + adNorm[1] * dEncDist;
  adCent1[2] = adCircCent[2] + adNorm[2] * dEncDist;

  adCent2[0] = adCircCent[0] - adNorm[0] * dEncDist;
  adCent2[1] = adCircCent[1] - adNorm[1] * dEncDist;
  adCent2[2] = adCircCent[2] - adNorm[2] * dEncDist;

  double dDist1 = dDIST3D(adCent1, adPoint);
  double dDist2 = dDIST3D(adCent2, adPoint);

  int iFC1 = iFuzzyComp(dDist1, dSphereRadius);
  int iFC2 = iFuzzyComp(dDist2, dSphereRadius);

  enum eEncroachResult eRetVal = eClean;
  // If both are -1, then it definitely encroaches
  if (iFC1 == -1 && iFC2 == -1)
    eRetVal = eSplit;
  // If one is negative and one zero, then use the tiebreaker
  if ((iFC1 == -1 && iFC2 == 0) || (iFC1 == 0 && iFC2 == -1))
    eRetVal = (qTieBreaker ? eSplit : eClean);

  // Note that if an internal bdry face has had offset insertion done in
  // one region, a request for offset insertion in the other region will
  // instead cause the face to be split.  This is a theoretical problem,
  // but not expected to be an issue in practice.
  if (eRetVal == eSplit && eET == eNewLens
      && pBTF->isOffsetInsertionAllowed()) {
    // Now check whether to do the offset point thing.  The offset point
    // is inserted iff the proposed new point lies too far (farther than
    // dRadius) from each of the corners.  Also, if the offset point
    // lies outside the domain, it isn't inserted.  In principle, the
    // offset point also shouldn't be inserted if it's across an
    // internal bdry from the cell that triggered the insertion.  In
    // practice, once encroached edges are fixed, either the bdry edge
    // with the incorrectly added offset pt will have been split -or-
    // the original triggering cell will still be there.  On the second
    // try, the BFace will be split.

    double dDistA = dDIST3D(adPoint, pF->getVert(0)->getCoords());
    double dDistB = dDIST3D(adPoint, pF->getVert(1)->getCoords());
    double dDistC = dDIST3D(adPoint, pF->getVert(2)->getCoords());
    // The following gives a provable angle bound of 30 degrees
    // (shortest edge and circumradius equal).

    int iFCA = iFuzzyComp(dDistA, dFaceRadius);
    int iFCB = iFuzzyComp(dDistB, dFaceRadius);
    int iFCC = iFuzzyComp(dDistC, dFaceRadius);
    if (iFCA == 1 && iFCB == 1 && iFCC == 1) {
      eRetVal = eOffsetPoint;
      pBTF->requestOffsetInsertion();

      double adDiff[] = adDIFF3D(adPoint, adCircCent);
      double dDot = dDOT3D(adNorm, adDiff);

      // Identify which apex is needed; it's the one in the interior, if
      // only one is.
      if (dDot < 0) {
	// Inside domain
	pBTF->setInsertionLocation(adCent2);
      }
      else if (pBTF->getType() == Cell::eIntTriBFace) {
	// Pick the closer apex
	assert(iFuzzyComp(dDist1, dDist2) != 0);
	if (dDist1 < dDist2) {
	  pBTF->setInsertionLocation(adCent1);
	}
	else {
	  pBTF->setInsertionLocation(adCent2);
	}
      }
      else {
	// Don't try to insert this thing outside the domain!
	pBTF->clearOffsetInsertionRequest();
      }
    }
  }
  return eRetVal;
}

enum eEncroachResult
TriBFaceBase::isPointEncroachingNewLens(const double adPoint[3],
					const bool qTieBreaker) const
{
  return eIsPointEncroachingGeneralLens(adPoint, qTieBreaker, this, eNewLens);
}

enum eEncroachResult
TriBFaceBase::isPointEncroachingLens(const double adPoint[3],
				     const bool qTieBreaker) const
{
  return eIsPointEncroachingGeneralLens(adPoint, qTieBreaker, this, eLens);
}

bool
TriBFaceBase::doesPointProjectInside(const double adPointIn[3]) const
{
  // Find out where the given point projects onto the plane of the bdry
  // face.
  const Face *pF = getFace();
  double adPoint[] =
    { adPointIn[0], adPointIn[1], adPointIn[2] };
  pF->projectOntoFace(adPoint);

  // Now find the barycentric coordinates of the projection in the tet
  // opposite this bdry face
  assert(pF->getOppositeCell(this)->getType() == Cell::eTet);
  TetCell *pTC = dynamic_cast<TetCell*>(pF->getOppositeCell(this));
  double adBary[4];
  pTC->calcBarycentricCoords(adPoint, adBary);

  // If all the barycentrics are zero or positive, great.  If one (or
  // more) is (are) negative, then the projection falls outside the
  // face.
  int aiBarySign[] =
    { iFuzzyComp(adBary[0], 0), iFuzzyComp(adBary[1], 0), iFuzzyComp(adBary[2],
								     0),
	iFuzzyComp(adBary[3], 0) };
  int iNumNeg = ((aiBarySign[0] == -1 ? 1 : 0) + (aiBarySign[1] == -1 ? 1 : 0)
      + (aiBarySign[2] == -1 ? 1 : 0) + (aiBarySign[3] == -1 ? 1 : 0));
  if (iNumNeg > 0)
    return false;

  // Also, projection onto a corner is worthless...
  int iNumZero = ((aiBarySign[0] == 0 ? 1 : 0) + (aiBarySign[1] == 0 ? 1 : 0)
      + (aiBarySign[2] == 0 ? 1 : 0) + (aiBarySign[3] == 0 ? 1 : 0));
  if (iNumZero == 3)
    return false;
  return true;
}

