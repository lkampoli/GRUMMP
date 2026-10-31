#include <math.h>
#include "GR_BdryPatch2D.h"
#include "GR_Bdry2D.h"
#include "GR_List.h"

// *******************************************************************
//
// BdrySeg
// -------
//
// This is the simplest descendant of BdryPatch2D.
// It implements features for linear boundary segments.
//
// *******************************************************************

BdrySeg::BdrySeg() :
    BdryPatch2D(NULL, iInvalidBC, iInvalidBC, iInvalidRegion,
    iInvalidRegion,
		-1, -1)
{
  eGeom2DType = ePolyLine;
}

BdrySeg::BdrySeg(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
		 const int iRL, const int iRR, const int iBPt, const int iEPt) :
    BdryPatch2D(pBdryObject, iBCL, iBCR, iRL, iRR, iBPt, iEPt)
{

  // Set boundary type
  eGeom2DType = ePolyLine;
  // Calculate information about segment...
  vRefresh();
}

void
BdrySeg::vAssignData(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
		     const int iRL, const int iRR, const int iBPt,
		     const int iEPt)
{
  BdryPatch2D::vAssignData(pBdryObject, iBCL, iBCR, iRL, iRR, iBPt, iEPt);
  // Compute data..
  vRefresh();
}

void
BdrySeg::vComputeMapping(BdryPatch const * pBP, int& iDim, double& dScaling,
			 double adTrans[3], double a2dRot[3][3]) const
{
  assert(dynamic_cast<const BdrySeg*>(pBP));
  // Space dimensions
  iDim = 2;
  const BdrySeg *pBS = static_cast<const BdrySeg*>(pBP);

  double adThisStart[2], adThisEnd[2];
  this->vGetBeginPoint(adThisStart);
  this->vGetEndPoint(adThisEnd);

  double dThisLen = hypot((adThisStart[0] - adThisEnd[0]),
			  (adThisStart[1] - adThisEnd[1]));
  double dThisAngle = atan2((adThisEnd[1] - adThisStart[1]),
			    (adThisEnd[0] - adThisStart[0]));

  double adOtherStart[2], adOtherEnd[2];
  pBS->vGetBeginPoint(adOtherEnd);
  pBS->vGetEndPoint(adOtherStart);

  double dOtherLen = hypot((adOtherStart[0] - adOtherEnd[0]),
			   (adOtherStart[1] - adOtherEnd[1]));
  double dOtherAngle = atan2((adOtherEnd[1] - adOtherStart[1]),
			     (adOtherEnd[0] - adOtherStart[0]));

  dScaling = dOtherLen / dThisLen;

  // Rotation is assumed to be an identity on input.
  double dAngleDiff = dOtherAngle - dThisAngle;
  double dCos = cos(dAngleDiff);
  double dSin = sin(dAngleDiff);
  a2dRot[0][0] = dCos;
  a2dRot[0][1] = -dSin;
  a2dRot[1][0] = dSin;
  a2dRot[1][1] = dCos;

  // Finally the translation, which is assumed to be zero initially.
  double adTestPt[] =
    { (adThisStart[0] * a2dRot[0][0] + adThisStart[1] * a2dRot[0][1])
	* dScaling, (adThisStart[0] * a2dRot[1][0]
	+ adThisStart[1] * a2dRot[1][1]) * dScaling };
  adTrans[0] = adOtherStart[0] - adTestPt[0];
  adTrans[1] = adOtherStart[1] - adTestPt[1];
}

void
BdrySeg::vComputeMaxParameter()
{
  // Compute the length of the segment.
  // For a linear segment, not very hard..
  double adE[2], adB[2];

  vGetBeginPoint(adB);
  vGetEndPoint(adE);

  double dX = adE[0] - adB[0];
  double dY = adE[1] - adB[1];

  dMaxParam = sqrt(dX * dX + dY * dY);
}

void
BdrySeg::vComputeNormal()
{
  // Since the normal is the same all along the
  // segment.. Can calculate it only once.

  // Must make sure that normal points in the direction of
  // the REGION, not the boundary..

  if (iRegionSide == REGION_SIDE_LEFT) {
    // the region is to the left

    // Can use adSlopeVec since it has been calculated before...
    adNormal[0] = -adSlopeVec[1];
    adNormal[1] = adSlopeVec[0];
  }
  else {
    // The region is to the right, or it's both..

    adNormal[0] = adSlopeVec[1];
    adNormal[1] = -adSlopeVec[0];
  }

}

void
BdrySeg::vComputeSlopeVec()
{
  // return the value of the slope of the segment...
  double dMag;
  double adB[2], adE[2];

  vGetBeginPoint(adB);
  vGetEndPoint(adE);

  adSlopeVec[0] = adE[0] - adB[0];
  adSlopeVec[1] = adE[1] - adB[1];
  dMag = dMAG2D(adSlopeVec);
  adSlopeVec[0] /= dMag;
  adSlopeVec[1] /= dMag;

}

void
BdrySeg::vPointAtParameterRatio(const double dRatio, double adLoc[2]) const
{
  // Return the point located at the ratio between pVBegin and pVEnd..
  // Ratio must be between 0 and 1 ...

  double adB[2], adE[2];

  // Make sure ratio is correct..
  assert((iFuzzyComp(dRatio, 0.0) >= 0) && (iFuzzyComp(dRatio, 1.0) <= 0));
  // Get the coords
  vGetBeginPoint(adB);
  vGetEndPoint(adE);
  // For a linear segment, this is rather simple...
  adLoc[0] = adB[0] + dRatio * (adE[0] - adB[0]);
  adLoc[1] = adB[1] + dRatio * (adE[1] - adB[1]);

}

void
BdrySeg::vPointAtDistFromStart(const double dDist, double adLoc[2]) const
{
  // For a linear segment, the distance is equal to the parameter so...
  vPointAtParameterRatio(dDist / dMaxParam, adLoc);
}

void
BdrySeg::vPointAtDistFromEnd(const double dDist, double adLoc[2]) const
{
  // Easy for linear segments...
  vPointAtParameterRatio((dMaxParam - dDist) / dMaxParam, adLoc);
}

double
BdrySeg::dCurvature(const double[2]) const
{
  // Calculate the curvature of a linear segment.
  // Actually, no calculation is needed since the curvature
  // of a segment is 0 everywhere...

  return (0.0);
}

double
BdrySeg::dParameterClosestToPoint(const double adLoc[2]) const
{
  // Returns the parameter from pVBegin to the point on the line segment
  // closest to the point given by adLoc[2]. The parameter might be negative
  // or greater than dMaxParam, depending on the location of adLoc,
  // i.e. the location given by the parameter calculated in here is not
  // restricted to be between pVBegin and pVEnd.
  double dDelX, dDelY, dNormRatio, ds;
  double adB[2];

  vGetBeginPoint(adB);

  dDelX = adLoc[0] - adB[0];
  dDelY = adLoc[1] - adB[1];

  if (fabs(adNormal[0]) < 1e-15) {
    // If the Normal is almost pointing straight up or down,
    // then the parameter is simply the x-distance...
    dNormRatio = 0.0;
    ds = fabs(dDelX);

  }
  else {
    dNormRatio = adNormal[1] / adNormal[0];

    ds = (dDelY - (dDelX * dNormRatio))
	/ (adSlopeVec[1] - (adSlopeVec[0] * dNormRatio));
  }
  return ds;
}

void
BdrySeg::vBdryPointClosestToPoint(const double adLoc[2],
				  double adPoint[2]) const
{
  // Returns the point ON THE LINE SEGMENT closest to the point
  // given in adLoc[2].
  double ds;
  double adB[2], adE[2];

  vGetBeginPoint(adB);
  vGetEndPoint(adE);

  ds = dParameterClosestToPoint(adLoc);
  // ds can be greater than dMaxParam or negative...
  // Must take that into account
  if (ds >= dMaxParam) {
    adPoint[0] = adE[0];
    adPoint[1] = adE[1];
    return;
  }
  if (ds <= 0.0) {
    adPoint[0] = adB[0];
    adPoint[1] = adB[1];
    return;
  }
  // Point must lie on segment..
  vPointAtParameter(ds, adPoint);
}

double
BdrySeg::dDistFromBdry(const double adLoc[2], int * const piSide) const
{
  // Return the SHORTEST distance from the line segment.
  // Also return if the point is located to the right or
  // to the left of the segment.

  double adPoint[2], dX, dY, dDot, adN[2];

  // Find closest point location..
  vBdryPointClosestToPoint(adLoc, adPoint);
  // Calculate distance..
  dX = adPoint[0] - adLoc[0];
  dY = adPoint[1] - adLoc[1];
  double dDist = sqrt(dX * dX + dY * dY);
  // This normal always points to the right of the curve
  adN[0] = adSlopeVec[1];
  adN[1] = -adSlopeVec[0];

  // Calculate if point is inside (where normal points) or outside
  // Use Dot product
  dDot = dX * adN[0] + dY * adN[1];

  // -ve dot product means to the right
  if (piSide) {
    if (dDot <= 0.0) {
      *piSide = SIDE_RIGHT;
    }
    else {
      *piSide = SIDE_LEFT;
    }
  }
  return dDist;
}

double
BdrySeg::dParameterAtPoint(const double adLoc[2]) const
{
  // Returns the parameter for the point...
  // This will actually return the parameter for the
  // on the boundary closest to the point given.

  return dParameterClosestToPoint(adLoc);
}

void
BdrySeg::vUnitNormal(const double[], double adNorm[]) const
{
  // Returns the normal at location adLoc[].

  vNormalAtParameterRatio(0.0, adNorm);
}

void
BdrySeg::vNormalAtParameterRatio(const double /*dRatio*/,
				 double adNorm[2]) const
{
  // Trivial for a linear segment...
  // Normal is supposed to be up-to-date
  adNorm[0] = adNormal[0];
  adNorm[1] = adNormal[1];
}

bool
BdrySeg::qIsColinear(BdryPatch2D *, Vert * /*pVCommon*/, int */*iCase*/)
{
  // Determines if the patch in pPatch is
  // colinear with this one or not...
  assert(0);
  return (false);
  /*
   bool qComVertEnd1, qComVertEnd2;
   double adSlope[2];
   double dDot;
   // First of all, if the type is not the same, definitely not colinear!
   if (pPatch->eGeom2DType != eGeom2DType)
   return false;


   // Look at which vertices are the common vertices..
   if (pVEndVert () == pVCommon)
   qComVertEnd1 = true;
   else
   qComVertEnd1 = false;
   if (pPatch->pVEndVert () == pVCommon)
   qComVertEnd2 = true;
   else
   qComVertEnd2 = false;

   if ((qComVertEnd1) && (qComVertEnd2)) {
   // Common vertex is on both end vertices...
   *iCase = 2;
   // Dot product should be -1
   ((BdrySeg *) pPatch)->vGetSlopeVec (adSlope);
   dDot = adSlope[0] * adSlopeVec[0] + adSlope[1] * adSlopeVec[1];
   if (iFuzzyComp (dDot, -1.0) == 0) {
   // Last check: are the regions and BCs the same?
   // Must look at opposite ones since they're not pointing
   // in the same direction
   if (
   (iBCLeft == pPatch->iRightBdryCond ()) &&
   (iBCRight == pPatch->iLeftBdryCond ()) &&
   (iRegionRight == pPatch->iLeftRegion ()) &&
   (iRegionLeft == pPatch->iRightRegion ())
   )
   return true;
   else
   return false;
   }
   else
   return false;
   }
   else {
   if ((qComVertEnd1) || (qComVertEnd2)) {
   // Common vertex is at the end for one, and at the beginning for the other..
   *iCase = 1;
   // Dot product should be 1
   ((BdrySeg *) pPatch)->vGetSlopeVec (adSlope);
   dDot = adSlope[0] * adSlopeVec[0] + adSlope[1] * adSlopeVec[1];
   if (iFuzzyComp (dDot, 1.0) == 0) {
   // Last check: are the regions and BCs the same?
   if (
   (iBCLeft == pPatch->iLeftBdryCond ()) &&
   (iBCRight == pPatch->iRightBdryCond ()) &&
   (iRegionRight == pPatch->iRightRegion ()) &&
   (iRegionLeft == pPatch->iLeftRegion ())
   )
   return true;
   else
   return false;
   }
   else
   return false;
   }
   else {
   // Common vertex is the beginning vertex for both patches
   *iCase = 3;
   // dot product should be -1
   ((BdrySeg *) pPatch)->vGetSlopeVec (adSlope);
   dDot = adSlope[0] * adSlopeVec[0] + adSlope[1] * adSlopeVec[1];
   if (iFuzzyComp (dDot, -1.0) == 0) {
   // Last check: are the regions and BCs the same?
   // Must look at opposite ones since they're not pointing
   // in the same direction
   if (
   (iBCLeft == pPatch->iRightBdryCond ()) &&
   (iBCRight == pPatch->iLeftBdryCond ()) &&
   (iRegionRight == pPatch->iLeftRegion ()) &&
   (iRegionLeft == pPatch->iRightRegion ())
   )
   return true;
   else
   return false;
   }
   else
   return false;
   }
   }

   */
}

List<double>
BdrySeg::LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop)
{
  // This returns the original discretization info
  // For a BdrySeg, there are no extra verts to insert,
  // so leave LsNewVerts intact.
  // There is only one edge to insert also...
  List<double> LdList;
  LdList.vClear();

  *qInsertVertPointers = true;
  *qCloseLoop = false;
  return LdList;
}

void
BdrySeg::vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
				const double adLoc2[2], double adPoint[2]) const
{
  // Return the location of the point on the boundary that is located
  // at ratio dRatio between adLoc1[] and adLoc2[].

  // Both these points "should" be on the boundary...	
#ifndef NDEBUG
  double dDist1 = dDistFromBdry(adLoc1);
  double dDist2 = dDistFromBdry(adLoc2);

  assert(iFuzzyComp(dDist1 / 100.0, 0.0) == 0);
  assert(iFuzzyComp(dDist2 / 100.0, 0.0) == 0);
#endif

  double dAL1, dAL2, dALRatio;

  // Get first parameter
  dAL1 = dParameterAtPoint(adLoc1);
  // Get the second one
  dAL2 = dParameterAtPoint(adLoc2);

  // Get the parameter we want
  dALRatio = dAL1 + dRatio * (dAL2 - dAL1);

  // Return the point
  vPointAtParameter(dALRatio, adPoint);
}

void
BdrySeg::vGeodesicDistanceRatioBetweenTwoPoints(double dRatio,
						const double adLoc1[2],
						const double adLoc2[2],
						double adPoint[2]) const
{
  vRatioBetweenTwoPoints(dRatio, adLoc1, adLoc2, adPoint);
}

double
BdrySeg::dXMax()
{
  double adB[2], adE[2];
  vGetBeginPoint(adB);
  vGetEndPoint(adE);

  // Returns the maximum X coordinate of the patch..
  if (adB[0] > adE[0])
    return adB[0];
  else
    return adE[0];
}

double
BdrySeg::dYMax()
{
  double adB[2], adE[2];
  vGetBeginPoint(adB);
  vGetEndPoint(adE);

  // Returns the maximum Y coordinate of the patch..
  if (adB[1] > adE[1])
    return adB[1];
  else
    return adE[1];
}

double
BdrySeg::dXMin()
{
  double adB[2], adE[2];
  vGetBeginPoint(adB);
  vGetEndPoint(adE);
  // Returns the minimum X coordinate of the patch..
  if (adB[0] < adE[0])
    return adB[0];
  else
    return adE[0];
}

double
BdrySeg::dYMin()
{
  double adB[2], adE[2];
  vGetBeginPoint(adB);
  vGetEndPoint(adE);
  // Returns the minimum Y coordinate of the patch..
  if (adB[1] < adE[1])
    return adB[1];
  else
    return adE[1];
}

void
BdrySeg::vRefresh()
{
  // Make sure the pre-calculated items are up-to-date...

  // Make sure the length is computed properly...                 
  vComputeMaxParameter();
  // Compute the slope vector
  vComputeSlopeVec();
  // Calculate the normal too
  vComputeNormal();
}

BdrySeg::~BdrySeg()
{
  // Destructor
}
