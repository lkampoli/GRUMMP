#ifndef GR_Bdry
#define GR_Bdry 

#include "GR_config.h"
#include "GR_Classes.h"

#include "GR_misc.h"

// This is the base class for all boundary data in both two and three
// dimensions.   
class BdryPatch {
  /// Copy construction disallowed
  BdryPatch(const BdryPatch&); // Never defined
  /// operator= disallowed; return value for picky compilers only
  BdryPatch&
  operator=(const BdryPatch&); // Never defined
protected:
  // Only one left and one right tag will be valid at any given time.
  int iBCLeft, iBCRight;
  int iRegionLeft, iRegionRight;
  bool qDel, qReversed;
public:
/// Constructor and destructor  
  BdryPatch() :
      iBCLeft(iInvalidBC), iBCRight(iInvalidBC), iRegionLeft(
      iInvalidRegion), iRegionRight(iInvalidRegion), qDel(false), qReversed(
	  false)
  {
  }
  BdryPatch(const int iBCL, const int iBCR, const int iRL, const int iRR) :
      iBCLeft(iBCL), iBCRight(iBCR), iRegionLeft(iRL), iRegionRight(iRR), qDel(
	  false), qReversed(false)
  {
    // Each side must be either a boundary condition or a region, but
    // not both.
    assert(
	( (iBCL > 0) && (iBCL != iInvalidBC) ) || ( (iRL >= 0) && (iRL < iMaxRegionLabel) ));
    assert(
	( (iBCR > 0) && (iBCR != iInvalidBC) ) || ( (iRR >= 0) && (iRR < iMaxRegionLabel) ));
    // Both sides can't be boundaries.
    assert(iBCL == iInvalidBC || iBCR == iInvalidBC);
  }
  virtual
  ~BdryPatch()
  {
  }
/// Inquiry functions
  int
  iBdryCond() const
  {
    assert(iBCLeft == iInvalidBC || iBCRight == iInvalidBC);
    if (iBCLeft == iInvalidBC)
      return iBCRight;
    else
      return iBCLeft;
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
  bool
  qIsBoundary() const
  {
    return iRegionLeft == iOutsideRegion || iRegionRight == iOutsideRegion;
  }
  bool
  qValid() const
  {
    return (this != NULL);
  }

  virtual void
  vComputeMapping(const BdryPatch* const pBP, int& iDim, double& dScale,
		  double adTrans[3], double a2dRot[3][3]) const = 0;

  virtual bool
  qDisjointFrom(const BdryPatch* const) const = 0;
  virtual bool
  qDisjointFrom(const Vert* const) const = 0;

  void
  vMarkDeleted()
  {
    qDel = true;
  }
  bool
  qDeleted() const
  {
    return qDel;
  }

  void
  vSwitchOrientation()
  {
    int iBCTmp = iBCLeft;
    iBCLeft = iBCRight;
    iBCRight = iBCTmp;
    int iRegTmp = iRegionLeft;
    iRegionLeft = iRegionRight;
    iRegionRight = iRegTmp;
    qReversed = true;
  }
  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const = 0;
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const = 0;
  virtual double
  dCurvature(const double[]) const
  {
    return 0;
  }
  virtual void
  vCurvature(const double[], double[2][2]) const
  {
    assert(0);
  }
};

typedef struct _location_ {
  double adLoc[3];
} Location;

class BdryRep {
protected:
  int iNPatches, iMaxNPatches;
public:
  BdryRep() :
      iNPatches(0), iMaxNPatches(0)
  {
  }
  virtual
  ~BdryRep()
  {
  }
  ;
  virtual int
  iNumPatches() const
  {
    return iNPatches;
  }
  virtual int
  iMaxNumPatches() const
  {
    return iMaxNPatches;
  }
};

#endif
