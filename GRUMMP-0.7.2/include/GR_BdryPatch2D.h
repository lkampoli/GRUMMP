#ifndef GR_BDRYPATCH_2D_H
#define GR_BDRYPATCH_2D_H 1

#include "GR_config.h"
#include <math.h>
#include "GR_misc.h"
#include "GR_Bdry2D.h"
#include "GR_List.h"

#define BC_SIDE_LEFT 1
#define BC_SIDE_RIGHT 2
#define BC_SIDE_INT_BDRY 3
#define BC_SIDE_NONE 0
#define REGION_SIDE_LEFT 1
#define REGION_SIDE_RIGHT 2
#define REGION_SIDE_BOTH 0
#define SIDE_LEFT 1
#define SIDE_RIGHT 2
#define SIDE_UNDETERMINED 3

// *******************************************************************
//
// BdryPatch2D
// -----------
//
// This is the base class for all 2D boundary data (including internal
// boundaries between regions in a multi-region domain).  Eventually
// this will have a flock of derived classes hanging from it.
//
// *******************************************************************

class BdryPatch2D : public BdryPatch {
  /// Copy construction and operator= disallowed.
  BdryPatch2D(const BdryPatch2D&) :
      BdryPatch(), iBCSide(0), iRegionSide(0), pB2D(NULL), iBeginPt(), iEndPt(), dMaxParam(), eGeom2DType()
  {
    assert(0);
  }
  BdryPatch2D&
  operator=(const BdryPatch2D&)
  {
    assert(0);
    return (*this);
  }
protected:
  // Only one left and one right tag will be valid at any given time.
  //  int iBCLeft, iBCRight;
  //  int iRegionLeft, iRegionRight;
  // 1 = Left; 2 = Right; 0 = none.
  int iBCSide;
  // 1 = Left; 2 = Right; 0 = both
  int iRegionSide;

  // I'd really prefer not to have to store a Bdry2D object in the first
  // place, but vSetPatch(BdryPatch*) is problematical at present
  // without this call.  Hopefully this can all be eliminated in v0.2.1.
  // FIX ME
  // Pointer to the Bdry2D object
  Bdry2D* pB2D;
public:
  Bdry2D*
  pBdry2D() const
  {
    return pB2D;
  }
protected:
  //  Vert *pVBegin, *pVEnd;
  // use indices instead
  int iBeginPt, iEndPt;
  // Length of the patch
  double dMaxParam;
  // Virtual function to compute the length. Must be called in constructor.
  virtual void
  vComputeMaxParameter()
  {
    dMaxParam = 0.0;
  }
public:
  /// Constructor and destructor
  BdryPatch2D() :
      BdryPatch(), iBCSide(0), iRegionSide(0), pB2D(NULL), iBeginPt(-1), iEndPt(
	  -1), dMaxParam(0), eGeom2DType(eUnknown)
  {
    //    pVBegin = pVEnd = pVInvalidVert;
  }
  BdryPatch2D(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt) :
      BdryPatch(iBCL, iBCR, iRL, iRR), iBCSide(0), iRegionSide(0), pB2D(NULL), iBeginPt(
	  -1), iEndPt(-1), dMaxParam(0), eGeom2DType(eUnknown)
  {
    // Each side must be either a boundary condition or a region, but
    // not both.
    vAssignData(pBdryObject, iBCL, iBCR, iRL, iRR, iBPt, iEPt);
  }

  virtual
  ~BdryPatch2D()
  {
  }

  virtual void
  vComputeMapping(const BdryPatch* const pBP, int& iDim, double& dScale,
		  double adTrans[3], double a2dRot[3][3]) const = 0;

  enum eGeom2D {
    eUnknown,
    ePolyLine,
    eCircle,
    eCircularArc,
    eLongArc,
    eCubicParam,
    eBezier,
    eInterpolated
  };

  eGeom2D eGeom2DType;
  /// Inquiry functions

  int
  iBeginPointIndex() const
  {
    return iBeginPt;
  }
  int
  iEndPointIndex() const
  {
    return iEndPt;
  }
  int
  iLeftBdryCond() const
  {
    return iBCLeft;
  }
  int
  iRightBdryCond() const
  {
    return iBCRight;
  }
  int
  iLeftRegion() const
  {
    return iRegionLeft;
  }
  int
  iRightRegion() const
  {
    return iRegionRight;
  }
  int
  iBCWhichSide() const
  {
    return iBCSide;
  }
  int
  iRegionWhichSide() const
  {
    return iRegionSide;
  }

  void
  vGetBeginPoint(double adLoc[2]) const
  {
    pB2D->vGetPoint(iBeginPt, adLoc);
  }
  void
  vGetEndPoint(double adLoc[2]) const
  {
    pB2D->vGetPoint(iEndPt, adLoc);
  }
  int
  iBdryCond() const
  {
    if (iBCSide == 0)
      return iInvalidBC;
    if (iBCSide == 1)
      return iLeftBdryCond();
    if (iBCSide == 2)
      return iRightBdryCond();
    // unreachable..
    return iInvalidBC;
  }

  double
  dMaxParameter()
  {
    return dMaxParam;
  }
  void
  vPointAtParameter(const double dParamIn, double adLoc[2]) const
  {
    // Return the point located at parameter dParameter..
    // Parameter is measure from pVBegin to pVEnd.
    double dRatio;
    // Compute the ratio for this parameter...
    dRatio = dParamIn / dMaxParam;
    // Use ratio function...
    vPointAtParameterRatio(dRatio, adLoc);
  }

  virtual void
  vAssignData(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt)
  {
    // This is for patches initialized with the empty constructor

    // Each side must be either a boundary condition or a region, but
    // not both.

    assert(
	( (iBCL > 0) && (iBCL != iInvalidBC) ) || ( (iRL >= 0) && (iRL < iMaxRegionLabel) ));
    assert(
	( (iBCR > 0) && (iBCR != iInvalidBC) ) || ( (iRR >= 0) && (iRR < iMaxRegionLabel) ));
    // Both sides can't be boundaries.
    assert(iBCL == iInvalidBC || iBCR == iInvalidBC);
    iBCLeft = iBCL;
    iBCRight = iBCR;
    iRegionLeft = iRL;
    iRegionRight = iRR;
    // Make sure the Bdry2D object is not null
    assert(pBdryObject != NULL);
    pB2D = pBdryObject;
    // The verts had better be valid.
    // Not necessary here, circles don't need them...
    //    assert(pVB != pVInvalidVert);
    //    assert(pVE != pVInvalidVert);
    iBeginPt = iBPt;
    iEndPt = iEPt;

    // Check what side the region is on..
    if (iRL != iInvalidRegion) {
      if (iRR != iInvalidRegion) {
	iRegionSide = REGION_SIDE_BOTH;
      }
      else {
	iRegionSide = REGION_SIDE_LEFT;
      }
    }
    else {
      if (iRR != iInvalidRegion) {
	iRegionSide = REGION_SIDE_RIGHT;
      }
      else {
	assert2(0, "ERROR: 2 invalid regions in BdryPatch2D!\n");
      }
    }
    // Check what side the BC is...
    if (iBCL == iInvalidBC) {
      if (iBCR == iInvalidBC) {
	// seems to be 2 regions...
	iBCSide = BC_SIDE_NONE;
      }
      else {
	// Left invalid, right not..
	iBCSide = BC_SIDE_RIGHT;
      }
    }
    else {
      // Just to make sure...
      if (iBCR == iInvalidBC) {
	iBCSide = BC_SIDE_LEFT;
      }
      else {
	// Both are not invalid?!
	assert2(0, "ERROR: 2 non-invalid BCs in BdryPatch2D!\n");
      }
    }
    // Compute data that needs to be computed..
    //    vRefresh();
  }

  virtual double
  dGeodesicDistance(const double dParam1, const double dParam2) const;
  virtual void
  vPointAtParameterRatio(const double dRatio, double adLoc[2]) const = 0;
  virtual void
  vPointAtDistFromStart(const double dDist, double adLoc[2]) const = 0;
  virtual void
  vPointAtDistFromEnd(const double dDist, double adLoc[2]) const = 0;
  // The default version of dSpeed assumes arclength parametrization.
  virtual double
  dSpeed(const double) const
  {
    return 1;
  }
  virtual double
  dCurvature(const double[2]) const
  {
    return 0.0;
  }
  virtual double
  dParameterAtPoint(const double[2]) const
  {
    return 0.0;
  }
  virtual void
  vBdryPointClosestToPoint(const double adLoc[2], double adPoint[2]) const = 0;
  virtual double
  dDistFromBdry(const double[], int * const = NULL) const
  {
    return 0.0;
  }
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const = 0;
  virtual void
  vNormalAtParameterRatio(const double dRatio, double adNorm[2]) const = 0;
  virtual void
  vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
			 const double adLoc2[2], double adPoint[2]) const = 0;
  virtual void
  vRefresh() = 0;
  virtual double
  dXMax() = 0;
  virtual double
  dYMax() = 0;
  virtual double
  dXMin() = 0;
  virtual double
  dYMin() = 0;
  // Call to find out if the two verts "see" each other
  // (does not take into account OTHER patches that might
  //  block the view)
  virtual bool
  qVertsVisible() = 0;
  virtual List<double>
  LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop) = 0;
  // iCase (the int arg) is modified this way:
  //    - iCase = 1:  If both patches are pointing in the same direction
  //    - iCase = 2:  If both patches are pointing towards the common vertex
  //    - iCase = 3:  If both patches are pointing towars the "other" vertex
  virtual bool
  qIsColinear(BdryPatch2D*, Vert*, int*)
  {
    return false;
  }
  virtual bool
  qDisjointFrom(const BdryPatch* const pBP) const;
  virtual bool
  qDisjointFrom(const int iPt) const;
  virtual bool
  qDisjointFrom(const Vert* const pV) const;
  //  virtual bool   qIsPointOnBdry(const double[2]) {return false;}

  virtual void
  vGeodesicDistanceRatioBetweenTwoPoints(double dRatio, const double adLoc1[2],
					 const double adLoc2[2],
					 double adPoint[2]) const;

};

// *******************************************************************
//
// BdrySeg
// -------
//
// This is the simplest descendant of BdryPatch2D.
// It implements features for linear boundary segments.
//
// *******************************************************************

class BdrySeg : public BdryPatch2D {
  // Disable the default constructor and copy constructor
  BdrySeg(const BdrySeg&) :
      BdryPatch2D()
  {
    assert(0);
  }
protected:
  // Normal is the same everywhere..
  double adNormal[2];
  // Slope also...
  double adSlopeVec[2];
  // Compute the length...
  virtual void
  vComputeMaxParameter();
  // The next two are not virtual because
  // slope and normals will not always be the same

  // Compute the slope...
  void
  vComputeSlopeVec();
  // Compute the normal
  void
  vComputeNormal();

  // Compute the stuff needed to map this patch to some other patch.
  void
  vComputeMapping(const BdryPatch* const pBP, int& iDim, double& dScale,
		  double adTrans[3], double a2dRot[3][3]) const;

  // Used to compute the closest point on the boundary..
  virtual double
  dParameterClosestToPoint(const double adLoc[2]) const;
protected:
  void
  vGetSlopeVec(double adSlope[2])
  {
    // Returns the slope vector..
    adSlope[0] = adSlopeVec[0];
    adSlope[1] = adSlopeVec[1];
  }
public:
  BdrySeg();
  BdrySeg(Bdry2D *pBdryObject, const int iBCL, const int iBCR, const int iRL,
	  const int iRR, const int iBPt, const int iEPt);
  virtual
  ~BdrySeg();

  virtual void
  vAssignData(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt);

  virtual double
  dGeodesicDistance(const double dParam1, const double dParam2) const
  {
    return (dParam2 - dParam1);
  } // Arclength parametrization
  virtual void
  vPointAtParameterRatio(const double dRatio, double adLoc[2]) const;
  virtual void
  vPointAtDistFromStart(const double dDist, double adLoc[2]) const;
  virtual void
  vPointAtDistFromEnd(const double dDist, double adLoc[2]) const;
  virtual double
  dCurvature(const double adLoc[2]) const;
  virtual void
  vBdryPointClosestToPoint(const double adLoc[2], double adPoint[2]) const;
  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const;
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const;
  virtual void
  vNormalAtParameterRatio(const double dRatio, double adNorm[2]) const;
  virtual void
  vRefresh();
  virtual double
  dParameterAtPoint(const double adLoc[2]) const;
  virtual void
  vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
			 const double adLoc2[2], double adPoint[2]) const;
  virtual double
  dXMax();
  virtual double
  dYMax();
  virtual double
  dXMin();
  virtual double
  dYMin();
  // Linear segments -- verts see each other
  virtual bool
  qVertsVisible()
  {
    return true;
  }
  ;

  virtual List<double>
  LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop);

  // iCase is modified this way:
  //    - iCase = 1:  If both patches are pointing in the same direction
  //    - iCase = 2:  If both patches are pointing towards the common vertex
  //    - iCase = 3:  If both patches are pointing towars the "other" vertex
  virtual bool
  qIsColinear(BdryPatch2D* pPatch, Vert* pVCommon, int* iCase);
  virtual void
  vGeodesicDistanceRatioBetweenTwoPoints(double dRatio, const double adLoc1[2],
					 const double adLoc2[2],
					 double adPoint[2]) const;

};

// *******************************************************************
//
// BdryArc
// -------
//
// Yet another descendant of BdryPatch2D.
// It implements features for circular arcs AND circles.
//
// *******************************************************************

class BdryArc : public BdryPatch2D {
  // Disable the copy constructor
  BdryArc(const BdryArc&) :
      BdryPatch2D(), iCenterPt(), dRadius(), dAngleStart(), dAngleEnd()
  {
    assert(0);
  }
protected:
  // Need a center
  // Not an index, because arcs don't get a center passed to them
  double adCenter[2];
  // Store if for circles however...
  int iCenterPt;
  // Need a radius too
  double dRadius;
  // Need to know the start/end orientations
  // These angles are stored in radians even though
  // the user supplies them in degrees.
  // Start = 0, End = 2Pi for a full circle.
  double dAngleStart;
  double dAngleEnd;
  // Compute the length...
  virtual void
  vComputeMaxParameter();
  // The next two are not virtual because
  // slope and normals will not always be the same

  // Used to compute the closest point on the boundary..
  virtual double
  dParameterClosestToPoint(const double adLoc[2]) const;
  virtual double
  dAngleClosestToPoint(const double adLoc[2]) const;

public:
  BdryArc();
  // Constructor for circles
  BdryArc(Bdry2D *pBdryObject, const int iBCL, const int iBCR, const int iRL,
	  const int iRR, const double dRad, const int iCPt);
  // constructor for arcs
  BdryArc(Bdry2D *pBdryObject, const int iBCL, const int iBCR, const int iRL,
	  const int iRR, const int iBPt, const int iEPt, const double dRad,
	  const bool qIsLongArc = false);
  // destructor
  virtual
  ~BdryArc();

  virtual void
  vAssignData(Bdry2D *pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt);

  void
  vComputeMapping(const BdryPatch* const, int&, double&, double[3],
		  double[3][3]) const
  {
    assert(0);
  }

  // Should be private or protected; need copy constructor to do this.
  int
  iCenterPointIndex()
  {
    return iCenterPt;
  }
  double
  dGetRadius()
  {
    return dRadius;
  }
  ;
  void
  vGetCenter(double adLoc[2])
  {
    adLoc[0] = adCenter[0];
    adLoc[1] = adCenter[1];
  }
  ;
  double
  dGetStartAngle()
  {
    return dAngleStart;
  }
  ;
  double
  dGetEndAngle()
  {
    return dAngleEnd;
  }
  ;
  // End should be private or protected

  virtual double
  dGeodesicDistance(const double dParam1, const double dParam2) const
  // Arclength parametrized; allow for wraparounds with more than a
  // difference of 2 pi in parameter.
  {
    return min(fabs(dParam2 - dParam1),
	       fabs(fabs(dParam2 - dParam1) - 2 * M_PI * dRadius));
  }
  void
  vPointAtParameter(const double dParameter, double adLoc[2]) const;
  virtual void
  vPointAtParameterRatio(const double dRatio, double adLoc[2]) const;
  virtual void
  vPointAtDistFromStart(const double dDist, double adLoc[2]) const;
  virtual void
  vPointAtDistFromEnd(const double dDist, double adLoc[2]) const;
  virtual double
  dCurvature(const double adLoc[2]) const;
  virtual void
  vBdryPointClosestToPoint(const double adLoc[2], double adPoint[2]) const;
  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const;
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const;
  virtual void
  vNormalAtParameterRatio(const double dRatio, double adNorm[2]) const;
  virtual void
  vRefresh();
  virtual double
  dParameterAtPoint(const double adLoc[2]) const;
  virtual void
  vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
			 const double adLoc2[2], double adPoint[2]) const;
  virtual double
  dXMax();
  virtual double
  dYMax();
  virtual double
  dXMin();
  virtual double
  dYMin();
  // Only true if the region is NOT to the right...
  virtual bool
  qVertsVisible()
  {
    if (iRegionWhichSide() != REGION_SIDE_RIGHT)
      return true;
    else
      return false;
  }
  ;
  virtual List<double>
  LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop);

  // iCase is modified this way:
  //    - iCase = 1:  If both patches are pointing in the same direction
  //    - iCase = 2:  If both patches are pointing towards the common vertex
  //    - iCase = 3:  If both patches are pointing towars the "other" vertex
  virtual bool
  qIsColinear(BdryPatch2D* pPatch, Vert* pVCommon, int* iCase);

  virtual void
  vGeodesicDistanceRatioBetweenTwoPoints(double dRatio, const double adLoc1[2],
					 const double adLoc2[2],
					 double adPoint[2]) const;

};

// *******************************************************************
//
// BdryCubicParam2D
// ----------------
//
// Another BdryPatch2D descendent, this class is the base class used
// to implement different types of cubic parametric curves in 2D.
//
// Different constructor can be used for different type of representations.
//
// *******************************************************************

class BdryCubicParam2D : public BdryPatch2D {
  // Disable the copy constructor
  BdryCubicParam2D(const BdryCubicParam2D&) :
      BdryPatch2D(), iEqnTypeX(), iEqnTypeY(), iXNumFirstRoots(), iYNumFirstRoots(), iXNumSecondRoots(), iYNumSecondRoots(), dOrientationChange(), LdOrderedRoots(), LdOrientationChangeAtRoots()
  {
    assert(0);
  }
private:

protected:
  // Internal variables needed for the type of boundary patch

  // Cubic parametric curves are of the form:
  // x(t) = ax t^3 + bx t^2 + cx t + dx
  // y(t) = ay t^3 + by t^2 + cy t + dy
  // Note that 0 <= t <= 1 ONLY

  double adCoeffs[2][4];
  // Info about the equations
  int iEqnTypeX, iEqnTypeY; // 0 = const, 1 = linear, 2 = quadratic, 3 = cubic

  int iXNumFirstRoots;  // Number of roots to the 1st derivative of the X eqn
  int iYNumFirstRoots;  // Number of roots to the 1st derivative of the Y eqn
  int iXNumSecondRoots; // Number of roots to the 2nd derivative of the X eqn
  int iYNumSecondRoots; // Number of roots to the 2nd derivative of the Y eqn

  double adFirstRoots[2][2];
  double adSecondRoots[2];

  double dOrientationChange;

  List<double> LdOrderedRoots;
  List<double> LdOrientationChangeAtRoots;

  // Compute the length...
  //  virtual void vComputeMaxParameter();

  // Functions associated with the internal variables

public:
  BdryCubicParam2D();
protected:
  // Constructor needed for the child classes
  BdryCubicParam2D(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
		   const int iRL, const int iRR, const int iBPt,
		   const int iEPt);
public:
  // Constructor for cubic parametrics curves
  // There is redundancy in the parameters given here
  // (i.e. some coefficients can be calculated from the endpoints),
  // but include them all as a general case
  BdryCubicParam2D(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
		   const int iRL, const int iRR, const int iBPt, const int iEPt,
		   const double dax, const double dbx, const double dcx,
		   const double ddx, const double day, const double dby,
		   const double dcy, const double ddy);
  // destructor
  virtual
  ~BdryCubicParam2D();

  //
  // --------------------------------------------------------------------
  // Parametric curve data-access functions
  //

private:
  //
  // --------------------------------------------------------------------
  // Parametric curve specific functions
  //

  void
  vDetermineEqnTypes();
  void
  vSolveQuadratic(const double dA, const double dB, const double dC,
		  int *iNumRoots, double *dRoot1, double *dRoot2) const;
  void
  vSolveCubic(double dA, double dB, double dC, double dD, int *iNumRoots,
	      double adRoots[3]) const;
  void
  vCalculateRoots();
  void
  vCalculateOrientationChange();
  virtual void
  vGetBezierPoints(double dP1[2], double dP2[2], double dP3[2],
		   double dP4[2]) const;

public:
  void
  vComputeMapping(const BdryPatch* const, int&, double&, double[3],
		  double[3][3]) const
  {
    assert(0);
  }

  // Should be private or protected; need copy constructor to do this.
  void
  vGetAllCoefficients(double *dAx, double *dBx, double *dCx, double *dDx,
		      double *dAy, double *dBy, double *dCy, double *dDy);
  void
  vSetAllCoefficients(const double dAx, const double dBx, const double dCx,
		      const double dDx, const double dAy, const double dBy,
		      const double dCy, const double dDy);
  double
  dGetOrientationChange();
  int
  iNumOrderedRoots()
  {
    return LdOrderedRoots.iLength();
  }
  double
  dGetOrderedRoot(const int iIndex)
  {
    assert((iIndex >= 0) && (iIndex < LdOrderedRoots.iLength() + 2));
    if (iIndex == 0)
      return 0.0;
    if (iIndex == LdOrderedRoots.iLength() + 1)
      return 1.0;
    return LdOrderedRoots[iIndex - 1];
  }
  double
  dGetOrientationChangeAtRoot(const int iIndex)
  {
    assert((iIndex >= 0) && (iIndex < LdOrderedRoots.iLength() + 2));
    if (iIndex == 0)
      return 0.0;
    if (iIndex == LdOrderedRoots.iLength() + 1)
      return dOrientationChange;
    return LdOrientationChangeAtRoots[iIndex - 1];
  }
  // End should be private or protected.
  virtual double
  dSpeed(const double dParam) const;
  double
  dOrientationChangeAtParameter(const double dt) const;
  double
  dParameterAtOrientationChange(const double dChange,
				const double dStart) const;
  void
  vTangentAtParameterRatio(const double dRatio, double a2dTangent[2]) const;
  void
  vParameterClosestToPoint(const double adLoc[2], double *dLen) const;

  //
  // --------------------------------------------------------------------
  // Re-implementation of functions defined in base class
  //

  virtual double
  dCurvature(const double adLoc[2]) const;
  //
  virtual double
  dParameterAtPoint(const double adLoc[2]) const;
  //
  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const;
  //
  virtual void
  vPointAtParameterRatio(const double dRatio, double adLoc[2]) const;
  //
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const;
  //
  virtual void
  vNormalAtParameterRatio(const double dRatio, double adNorm[2]) const;
  //
  virtual void
  vAssignData(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt);
  //
  virtual void
  vRefresh();
  //
  virtual void
  vPointAtDistFromStart(const double dDist, double adLoc[2]) const;
  //
  virtual void
  vPointAtDistFromEnd(const double dDist, double adLoc[2]) const;
  //
  virtual void
  vBdryPointClosestToPoint(const double adLoc[2], double adPoint[2]) const;
  //
  virtual void
  vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
			 const double adLoc2[2], double adPoint[2]) const;
  //
  virtual double
  dXMax();
  //
  virtual double
  dYMax();
  //
  virtual double
  dXMin();
  //
  virtual double
  dYMin();

  virtual List<double>
  LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop);
  // Only true if the region is NOT to the right...
  virtual bool
  qVertsVisible()
  {
    if (iRegionWhichSide() != REGION_SIDE_RIGHT)
      return true;
    else
      return false;
  }
  ;

  // iCase is modified this way:
  //    - iCase = 1:  If both patches are pointing in the same direction
  //    - iCase = 2:  If both patches are pointing towards the common vertex
  //    - iCase = 3:  If both patches are pointing towars the "other" vertex
  virtual bool
  qIsColinear(BdryPatch2D* pPatch, Vert* pVCommon, int* iCase);

  /*
   virtual double dCurvature(const double adLoc[2]);
   virtual double dDistFromBdry(const double adLoc[2]);

   */
};

// *******************************************************************
//
// BdryBezier2D
// ----------------
//
// This class is a BdryCubicParam2D descendant. The only thing that has
// to be overloaded is the way the cubic coefficients are calculated
//
// *******************************************************************

class BdryBezier2D : public BdryCubicParam2D {
  // Disable the copy constructor
  BdryBezier2D(const BdryBezier2D&) :
      BdryCubicParam2D(), iCPt1(), iCPt2()
  {
    assert(0);
  }
private:
  // Control points
  int iCPt1, iCPt2;
  double adControl1[2], adControl2[2];

  //
  // Access functions
  //

  virtual void
  vGetBezierPoints(double dP1[2], double dP2[2], double dP3[2],
		   double dP4[2]) const;
public:
  BdryBezier2D();
  // Constructor for Bezier curves. Only need the end
  // points and two control points (adC1, adC2). adC1 is related to pVB,
  // and adC2 is related to pVE
  BdryBezier2D(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
	       const int iRL, const int iRR, const int iBPt, const int iEPt,
	       const int iC1, const int iC2);
  // destructor
  virtual
  ~BdryBezier2D();

  void
  vComputeMapping(const BdryPatch* const, int&, double&, double[3],
		  double[3][3]) const
  {
    assert(0);
  }

  // Should be private or protected; need copy constructor to do this.
  int
  iControlPointIndex(const int iWhich)
  {
    assert((iWhich > 0) && (iWhich < 3));
    if (iWhich == 1)
      return iCPt1;
    else
      return iCPt2;
  }
  // End should be private or protected.

  //
  // Base class functions that are overriden
  //

  virtual double
  dSpeed(const double dParam) const
  {
    return BdryCubicParam2D::dSpeed(dParam);
  }
  virtual void
  vAssignData(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt);

};

// *******************************************************************
//
// BdryInterpolated2D
// ----------------
//
// Another BdryPatch2D descendent, this class is the base class used
// to implement different types of cubic parametric curves in 2D.
//
// Different constructor can be used for different type of representations.
//
// *******************************************************************

class BdryInterpolated2D : public BdryPatch2D {
private:
  // Disable the copy constructor
  BdryInterpolated2D(const BdryInterpolated2D&) :
      BdryPatch2D(), LpBCP2D(), LiPoints(), LdOrderedRoots(), LdOrientationChangeAtRoots(), dOrientationChange(), qUseLastCoeffs(), iLastNewPoint()
  {
    assert(0);
  }

protected:
  // Internal variables needed for the type of boundary patch

  // We need a list of pointers to cubic patches
  List<BdryCubicParam2D*> LpBCP2D;
  // We also need a list of points
  List<int> LiPoints;
  // Orientation change stuff
  List<double> LdOrderedRoots;
  List<double> LdOrientationChangeAtRoots;
  // Total orientation change
  double dOrientationChange;

  // Used when computing a new spline
  bool qUseLastCoeffs;
  double dLastCoeffs[2][4];
  int iLastNewPoint;

public:
  BdryInterpolated2D();
  // Constructor for the interpolated curves.
  // LiThePoints contains ALL the points necessary to build the curves
  // LiThePoints[0] is iBPt, and LiThePoints[LiThePoints.iLength()-1] is iEPt
  BdryInterpolated2D(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
		     const int iRL, const int iRR, const List<int> LiThePoints);
  // destructor
  virtual
  ~BdryInterpolated2D();

  void
  vComputeMapping(const BdryPatch* const, int&, double&, double[3],
		  double[3][3]) const
  {
    assert(0);
  }

  //
  // --------------------------------------------------------------------
  // Interpolated curve data-access functions
  //
  //
  BdryCubicParam2D*
  pBCP2DGetPatch(const int iIndex)
  {
    assert((iIndex >= 0) && (iIndex < LpBCP2D.iLength()));
    return LpBCP2D[iIndex];
  }

private:
  //
  // --------------------------------------------------------------------
  // Interpolated curve specific functions
  //

  void
  vFindBCoeffs(const int iNumRows, double *adR);
  void
  vComputeAllSplines();
  void
  vGetSplineFromParameter(const double dAL, int *iSpline,
			  double *dRemainingAL) const;
  double
  dCalculateNewShortSplinesChange();
public:
  // Should be private or protected; need copy constructor to do this.
  int
  iNumSplines()
  {
    return LpBCP2D.iLength();
  }
  int
  iNumPoints()
  {
    return LiPoints.iLength();
  }
  int
  iGetPoint(const int iIndex)
  {
    assert((iIndex >= 0) && (iIndex < LiPoints.iLength()));
    return LiPoints[iIndex];
  }
  void
  vGetSplineCoefficients(const int iSpline, double *dAx, double *dBx,
			 double *dCx, double *dDx, double *dAy, double *dBy,
			 double *dCy, double *dDy);
  void
  vSetSplineCoefficients(const int iSpline, const double dAx, const double dBx,
			 const double dCx, const double dDx, const double dAy,
			 const double dBy, const double dCy, const double dDy);
  double
  dTotalOrientationChange()
  {
    return dOrientationChange;
  }
  ;
  void
  vSplitSplineAtParameter(const int iSpline, const double dSplit, double *dAx1,
			  double *dBx1, double *dCx1, double *dDx1,
			  double *dAy1, double *dBy1, double *dCy1,
			  double *dDy1, double *dAx2, double *dBx2,
			  double *dCx2, double *dDx2, double *dAy2,
			  double *dBy2, double *dCy2, double *dDy2);
  int
  iNumOrderedRoots()
  {
    return LdOrderedRoots.iLength();
  }
  ;
  double
  dGetOrderedRoot(const int iIndex)
  {
    assert((iIndex >= 0) && (iIndex < LdOrderedRoots.iLength()));
    return LdOrderedRoots[iIndex];
  }
  // End should be private or protected.

public:
  double
  dOrientationChangeAtParameter(const double dt) const;
  double
  dParameterAtOrientationChange(const double dChange,
				const double dStart) const;
  void
  vTangentAtParameterRatio(const double dRatio, double a2dTangent[2]);
  //
  // --------------------------------------------------------------------
  // Re-implementation of functions defined in base class
  //

  // This one must be overridden so that we can integrate each separate
  // cubic parametric separately.  Necessary to get the right answer.
  virtual double
  dGeodesicDistance(const double dParam1, const double dParam2) const;
  virtual double
  dSpeed(const double dParam) const;
  virtual double
  dCurvature(const double adLoc[2]) const;
  //
  virtual double
  dParameterAtPoint(const double adLoc[2]) const;
  //
  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const;
  //
  virtual void
  vPointAtParameterRatio(const double dRatio, double adLoc[2]) const;
  //
  virtual void
  vPointAtParameter(const double dAL, double adLoc[2]) const;
  //
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const;
  //
  virtual void
  vNormalAtParameterRatio(const double dRatio, double adNorm[2]) const;
  //
  virtual void
  vAssignData(Bdry2D* pBdryObject, const int iBCL, const int iBCR,
	      const int iRL, const int iRR, const int iBPt, const int iEPt);
  //
  virtual void
  vRefresh();

  //
  virtual void
  vPointAtDistFromStart(const double dDist, double adLoc[2]) const;
  //
  virtual void
  vPointAtDistFromEnd(const double dDist, double adLoc[2]) const;
  //
  virtual void
  vBdryPointClosestToPoint(const double adLoc[2], double adPoint[2]) const;
  //
  virtual void
  vRatioBetweenTwoPoints(const double dRatio, const double adLoc1[2],
			 const double adLoc2[2], double adPoint[2]) const;
  //
  virtual double
  dXMax();
  //
  virtual double
  dYMax();
  //
  virtual double
  dXMin();
  //
  virtual double
  dYMin();

  virtual List<double>
  LdOriginalDiscretization(bool *qInsertVertPointers, bool *qCloseLoop);
  // Only true if the region is NOT to the right...
  virtual bool
  qVertsVisible()
  {
    if (iRegionWhichSide() != REGION_SIDE_RIGHT)
      return true;
    else
      return false;
  }
  ;

  // iCase is modified this way:
  //    - iCase = 1:  If both patches are pointing in the same direction
  //    - iCase = 2:  If both patches are pointing towards the common vertex
  //    - iCase = 3:  If both patches are pointing towars the "other" vertex
  virtual bool
  qIsColinear(BdryPatch2D* pPatch, Vert* pVCommon, int* iCase);
};

#endif
