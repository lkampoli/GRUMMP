#include <math.h>
#include "GR_Geometry.h"
#include "GR_BdryPatch2D.h"
#include "GR_Bdry2D.h"
#include "GR_Bdry3D.h"

// *******************************************************************
//
// BdryPolygon
// -----------
//
// This is the simplest descendant of BdryPatch3D.
// It implements features for planar polygons, possibly with holes.
//
// *******************************************************************

BdryPolygon::BdryPolygon(const int iBCL, const int iBCR, const int iRL,
			 const int iRR, const int iNV, Vert *apVV[]) :
    BdryPatch3D(iBCL, iBCR, iRL, iRR, iNV, apVV)
{
  // Must be careful to actually determine the polygon orientation, as
  // just checking a single angle can cause serious problems.  Start by
  // finding a vector perpendicular to the plan of the polygon, though,
  // so you can have a reference for the sign of the cross-product.
  Vert *pVA = apVV[iNVerts - 2];
  Vert *pVB = apVV[iNVerts - 1];
  Vert *pVC = apVV[0];
  calcNormal3D(pVA->getCoords(), pVB->getCoords(), pVC->getCoords(), adNormal);

  int i = 0;
  while (i < iNVerts && iFuzzyComp(dMAG3D(adNormal), 0) == 0) {
    pVA = pVB;
    pVB = pVC;
    pVC = apVV[i++];

    calcNormal3D(pVA->getCoords(), pVB->getCoords(), pVC->getCoords(),
		 adNormal);
  }
  NORMALIZE3D(adNormal);

  double dTotalAngle = 0;
  for (i = 0; i < iNVerts; i++) {
    pVA = apVV[i];
    pVB = apVV[(i + 1) % iNVerts];
    pVC = apVV[(i + 2) % iNVerts];

    double adVecAtoB[] = adDIFF3D(pVB->getCoords(), pVA->getCoords());
    double adVecBtoC[] = adDIFF3D(pVC->getCoords(), pVB->getCoords());

    double adCross[3];
    vCROSS3D(adVecAtoB, adVecBtoC, adCross);
    double dDot = dDOT3D(adVecAtoB, adVecBtoC);

    double dCross = dDOT3D(adCross, adNormal);

    double dThisAngle = atan2(dCross, dDot);
    dTotalAngle += dThisAngle;
    logMessage(
	4, "Point %3d. Normalized angle: %6f; normalized total so far %6f.\n",
	i, dThisAngle / (2 * M_PI), dTotalAngle / (2 * M_PI));
  }
  // This had better be one full rotation, one way or the other.
  assert(iFuzzyComp(fabs(dTotalAngle/(2*M_PI)), 1) == 0);
  if (dTotalAngle < 0) {
    vSCALE3D(adNormal, -1);
  }
  // Ensure that the normal always points out for non-internal patches.
//   if (iRL == iOutsideRegion)
  vSCALE3D(adNormal, -1);
}

void
BdryPolygon::vComputeMapping(const BdryPatch* const /* pBP */, int& /* iDim */,
			     double& /* dScale */, double /* adTrans */[3],
			     double /* a2dRot */[3][3]) const
{
  assert(0);
}

void
BdryPolygon::vFindProjectionBasis(double adOffset[3], double adBasis1[3],
				  double adBasis2[3]) const
{
  // Find the normal to the BdryPolygon
  double adNorm[3];
  adOffset[0] = pVVert(0)->x();
  adOffset[1] = pVVert(0)->y();
  adOffset[2] = pVVert(0)->z();
  vUnitNormal(adOffset, adNorm); // Could use any location

  // Now find two more vectors that are perpendicular to each other and
  // to the normal.

  // Find a single basis vector by using cross products with the normal
  // and trial vecs.
  {
    double adTry[] =
      { 1, .3, .7 }; // Unit coordinate vectors are too
    // commonly present.
    vCROSS3D(adNorm, adTry, adBasis1);
    if (dMAG3D(adBasis1) < 1.e-12) {
      // Hit the jackpot!  adNorm and adTry are parallel!
      adTry[0] += 1;
      // Now they won't be.
      vCROSS3D(adNorm, adTry, adBasis1);
    }

    // We want Basis1 x Basis2 = Norm.  So Norm x Basis1 = Basis2.
    vCROSS3D(adNorm, adBasis1, adBasis2);

    // Now normalize these things.
    NORMALIZE3D(adBasis1);
    NORMALIZE3D(adBasis2);
    assert(0 == iFuzzyComp(dDOT3D(adNorm, adBasis1), 0.));
    assert(0 == iFuzzyComp(dDOT3D(adNorm, adBasis2), 0.));
    assert(0 == iFuzzyComp(dDOT3D(adBasis1, adBasis2), 0.));
  } // Got the basis vectors
}

void
BdryPolygon::B2DProjection(const Bdry3D& B3D, Bdry2D& B2D) const
{
  double adOffset[3], adBasis1[3], adBasis2[3];
  vFindProjectionBasis(adOffset, adBasis1, adBasis2);

  // Set up a Bdry2D with enough vertex storage for all the Bdry3D
  // storage, and enough patch storage for all the edges in the
  // BdryPolygon.
  B2D.vSetSize(iNumSegments(), B3D.iNumVerts());

  // Project all points in the B3D into the plane where the B2D is being
  // set up.
  for (GR_index_t iV = 0; iV < B3D.iNumVerts(); iV++) {
    Vert *pV = B3D.pVVert(iV);
    const double * const adLoc3D = pV->getCoords();

    // To project, subtract the offset distance (adOffset).  Then take
    // dot products between that vector (which lies in a plane through
    // the origin) and the two basis vectors (both of which lie in that
    // plane).
    double adInPlane[] = adDIFF3D(adLoc3D, adOffset);
    double adLoc2D[] =
      { dDOT3D(adInPlane, adBasis1), dDOT3D(adInPlane, adBasis2) };

    B2D.vAddPoint(adLoc2D);
    B2D.vSetBdryPoint(iV, false);
  }

  // Copy BdryPolygon connectivity to the B2D, with each segment
  // becoming a patch.  Don't bother with BC and region data, as it
  // isn't relevant.  At the same time,
  for (int iP = 0; iP < iNumSegments(); iP++) {
    Vert *pV0, *pV1;
    vGetSegmentVerts(iP, pV0, pV1);
    int iGlobalVert0 = B3D.iVertIndex(pV0);
    int iGlobalVert1 = B3D.iVertIndex(pV1);

    BdryPatch2D *pPatch = new BdrySeg(&B2D, iDefaultBC, iInvalidBC,
    iInvalidRegion,
				      iDefaultRegion, iGlobalVert0,
				      iGlobalVert1);
    B2D.vAddPatchToList(pPatch);

    double adPt[2];
    double adVertPt[2];
    double adNorm;
    double dTanX, dTanY;
    structVertConnect sVConnect;

    // Might as well build the connectivity list as we go along, right?
    sVConnect.iPatch = B2D.iNumPatches() - 1;
    sVConnect.iVert = iGlobalVert0;
    B2D.vSetBdryPoint(iGlobalVert0, true);

    // Get the tangent vector for iGlobalVert0
    pPatch->vPointAtParameterRatio(0.01, adPt);
    // Get iGlobalVert0...
    B2D.vGetPoint(iGlobalVert0, adVertPt);
    dTanX = adPt[0] - adVertPt[0];
    dTanY = adPt[1] - adVertPt[1];
    adNorm = sqrt(dTanX * dTanX + dTanY * dTanY);
    dTanX /= adNorm;
    dTanY /= adNorm;
    sVConnect.dAngle = atan2(dTanY, dTanX);
    // Add it to list
    B2D.vAddToConnectivityList(sVConnect);

    // Do the same for iGlobalVert1 too
    sVConnect.iPatch = B2D.iNumPatches() - 1;
    sVConnect.iVert = iGlobalVert1;
    B2D.vSetBdryPoint(iGlobalVert1, true);
    // Get the tangent vector for iGlobalVert1
    pPatch->vPointAtParameterRatio(0.99, adPt);
    B2D.vGetPoint(iGlobalVert1, adVertPt);
    dTanX = adPt[0] - adVertPt[0];
    dTanY = adPt[1] - adVertPt[1];
    adNorm = sqrt(dTanX * dTanX + dTanY * dTanY);
    dTanX /= adNorm;
    dTanY /= adNorm;
    sVConnect.dAngle = atan2(dTanY, dTanX);
    // Add it to list
    B2D.vAddToConnectivityList(sVConnect);

  }
}

int
BdryPolygon::iLocalVertIndex(const double adLoc[3]) const
// Find which vertex in the polygon matches the given location.
{
  int iMin = -1;
  double dMinDist = 1.e300, d2ndMin = 1.e300;
  for (int i = 0; i < iNVerts; i++) {
    Vert *pV = pVVert(i);
    double dThisDist = dDIST3D(adLoc, pV->getCoords());
    if (dThisDist < dMinDist) {
      d2ndMin = dMinDist;
      dMinDist = dThisDist;
      iMin = i;
    }
    else if (dThisDist < d2ndMin) {
      d2ndMin = dThisDist;
    }
    if (dMinDist < 1.e-12)
      return iMin;

  }
  // If we've gotten here, then we should check any children that exist
  BdryPolygon* pBPolySib = dynamic_cast<BdryPolygon*>(pB3DChild);
  if (pBPolySib) {
    int iCumTotal = 0;
    const BdryPolygon* pBPolyPrev = this;
    int iRetVal = -1;
    do {
      iCumTotal += pBPolyPrev->iNVerts;
      iRetVal = pBPolySib->iLocalVertIndex(adLoc);
      pBPolyPrev = pBPolySib;
      pBPolySib = dynamic_cast<BdryPolygon*>(pBPolySib->pB3DNextSib);
    }
    while (iRetVal == -1 && pBPolySib);
    if (iRetVal >= 0)
      return (iCumTotal + iRetVal);
  }
  if (dMinDist / d2ndMin < 1.e-2)
    return iMin;
  else
    return -1;
}

double
BdryPolygon::dDistFromBdry(const double adLoc[], int * const piSide) const
{
  double adNorm[3];
  vUnitNormal(adLoc, adNorm);
  // piSide is used only for checking whether a pt is in front or behind
  // a patch for length scale calculations.

  double adDiff[] = adDIFF3D(adLoc, apVVerts[0]->getCoords());
  double dSignedNormDist = dDOT3D(adDiff, adNorm);

  // Try to handle cases where the distance is negligible compared with
  // local dimensions in the BdryPatch.
  double dDot1 = dDOT3D(adLoc, adNorm);
  double dDot2 = dDOT3D(apVVerts[0]->getCoords(), adNorm);
  if (iFuzzyComp(dDot1, dDot2) == 0)
    dSignedNormDist = 0;
  // Normal is directed outwards, at least for BdryPolygon, so flip the
  // sign to get positive inside.
  if (piSide)
    *piSide = -iFuzzyComp(dSignedNormDist, 0);

  double adProj[3];
  adProj[XDIR] = adLoc[XDIR] - dSignedNormDist * adNorm[XDIR];
  adProj[YDIR] = adLoc[YDIR] - dSignedNormDist * adNorm[YDIR];
  adProj[ZDIR] = adLoc[ZDIR] - dSignedNormDist * adNorm[ZDIR];

  double dNetAngle = 0;
  // This loop will work with -all- segments of the polygon and its
  // children.
  int iSeg;
  for (iSeg = 0; iSeg < iNumSegments(); iSeg++) {
    Vert *pV0, *pV1;
    vGetSegmentVerts(iSeg, pV0, pV1);
    double adDiff0[] = adDIFF3D(adProj, pV0->getCoords());
    double adDiff1[] = adDIFF3D(adProj, pV1->getCoords());
    double dLen0 = dMAG3D(adDiff0);
    double dLen1 = dMAG3D(adDiff1);
    if (iFuzzyComp(dLen0, 0) == 0 || iFuzzyComp(dLen1, 0) == 0) {
      continue;
    }
    double dCos = dDOT3D(adDiff0, adDiff1);
    double adCross[3];
    vCROSS3D(adDiff0, adDiff1, adCross);
    double dSin = dMAG3D(adCross);
    // Take orientation into account
    if (dDOT3D(adCross, adNorm) < 0)
      dSin = -dSin;
    double dThisAngle = atan2(dSin, dCos);
    // The following conditional is true if adProj is on an edge of the
    // polygon.
    if (iFuzzyComp(fabs(dThisAngle), M_PI) == 0)
      return fabs(dSignedNormDist);
    dNetAngle += dThisAngle;
  }
  double dRotations = fabs(dNetAngle / (2 * M_PI));

  // If the projection falls inside the poly or on its boundary, return
  // the distance from the plane.  In principle, the winding number is
  // 0.5 for a point on the boundary, 1 if inside, and 0 if outside.
  // The following inequality leaves lots of wiggle room.
  if (dRotations > 0.25) {
    return (fabs(dSignedNormDist));
  }

  // Otherwise, do it the hard way:  find the distance to each
  // individual segment and take the min.
  else {
    double dMinDist = LARGE_DBL;
    for (iSeg = 0; iSeg < iNumSegments(); iSeg++) {
      Vert *pV0, *pV1;
      vGetSegmentVerts(iSeg, pV0, pV1);
#if 0
      // This is the old code for computing distance to an edge.
      double adDiffProj[] = adDIFF3D(adProj, pV0->getCoords());
      double adVec[] = adDIFF3D(pV0->getCoords(),
	  pV1->getCoords());

      // Find an in-plane normals to the edge from the cross product of
      // the edge with the poly normal.
      double adSecondNorm[3];
      vCROSS3D(adNorm, adVec, adSecondNorm);
      NORMALIZE3D(adSecondNorm);
      double dDotProj = dDOT3D(adDiffProj, adSecondNorm);

      // Project adLoc onto the line containing iPrev and iNext.
      double adEdgeProj[] =
      { adProj[XDIR] - dDotProj*adSecondNorm[XDIR],
	adProj[YDIR] - dDotProj*adSecondNorm[YDIR],
	adProj[ZDIR] - dDotProj*adSecondNorm[ZDIR]};

      // Find the closest point along the edge to adLoc --- the
      // projection if inside; otherwise an endpoint.
      double adDiff0[] = adDIFF3D(adEdgeProj,
	  pV0->getCoords());
      double adDiff1[] = adDIFF3D(adEdgeProj,
	  pV1->getCoords());
      const double *adNearestPoint;
      if (dDOT3D(adDiff0, adDiff1) < 0) {
	// Point projects inside edge!
	adNearestPoint = adEdgeProj;
      }
      else {
	// Point projects outside edge.
	if (dMAG3D(adDiff0) < dMAG3D(adDiff1)) {
	  adNearestPoint = pV0->getCoords();
	}
	else {
	  adNearestPoint = pV1->getCoords();
	}
      }

      // Now (finally) compute the distance from adLoc to the nearest
      // point on the edge iNext, iPrev.
      double dDistance = dDIST3D(adNearestPoint, adLoc);
#else
      // Here's the new code.
      double adNewDiff[] = adDIFF3D(adLoc, pV0->getCoords());
      double adVec[] = adDIFF3D(pV1->getCoords(),
	  pV0->getCoords());

      double dDot = dDOT3D(adNewDiff, adVec);
      double dEdgeLenSq = dDOT3D(adVec, adVec);

      // T \equiv diff.vec/(vec.vec) lies between 0 and 1 iff the
      // projection of adLoc into the line of the edge lies inside the
      // edge.  Don't need to compute this explicitly.

      double dDistance;
      if (dDot < 0) {
	dDistance = dDIST3D(adLoc, pV0->getCoords());
      }
      else if (dDot > dEdgeLenSq) {
	dDistance = dDIST3D(adLoc, pV1->getCoords());
      }
      else {
	double dOtherLenSq = dDOT3D(adNewDiff, adNewDiff);
	dDistance = sqrt(dOtherLenSq - (dDot * dDot) / dEdgeLenSq);
      }
#endif
      dMinDist = min(dDistance, dMinDist);
    } // Done checking all edges
      // Sometimes causes problems in conjunction with patch optimization.
//     assert(dMinDist > 0);
    return dMinDist;
  } // Projection of adLoc falls outside the polygon
}

void
BdryPolygon::vUnitNormal(const double[], double adNorm[]) const
{
  // Always return the surface normal.
  if (qReversed) {
    adNorm[0] = -adNormal[0];
    adNorm[1] = -adNormal[1];
    adNorm[2] = -adNormal[2];
  }
  else {
    adNorm[0] = adNormal[0];
    adNorm[1] = adNormal[1];
    adNorm[2] = adNormal[2];
  }
}

void
BdryPolygon::vFlipNormal()
{
  adNormal[0] *= -1;
  adNormal[1] *= -1;
  adNormal[2] *= -1;
}
