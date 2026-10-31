#ifndef GR_BDRY_2D_H

#define GR_BDRY_2D_H 1

#include <math.h>
#include "GR_config.h"

#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Vertex.h"
#include "GR_Bdry.h"
#include "GR_List.h"

//  typedef struct {
//      bool    qIsUsed;
//      int     iWasReplacedBy;
//  } structPatchUsed;

// Used for the list of verts -> patches connected to them
class structVertConnect {
public:
  int iVert;    // Vertex number in Bdry2D
  int iPatch;   // Patch number in Bdry2D;
  double dAngle;   // Angle of the tangent vector originating from iVert;
  bool
  operator<(const structVertConnect& sVC) const
  {
    return (iVert < sVC.iVert || (iVert == sVC.iVert && iPatch < sVC.iPatch)
	|| (iVert == sVC.iVert && iPatch == sVC.iPatch && dAngle < sVC.dAngle));
  }
};

// *******************************************************************
// Bdry2D
// ------
//
// Bdry2D is the collection of all 2D boundary patches.
// These segments are all descendants of BdryPatch2D.
//
// *******************************************************************

class Bdry2D : public BdryRep {
private:
  // Size of the array to hold patches...
  int iSz;
  // Make this a list since sometimes, a patch yields a number of other
  // patches...  
  List<BdryPatch2D*> LpBPatches;
  // Here too, need a list, because we might add some as we go
  // CFO-G, 4/01:  This should actually be a vector, or something.  This
  // will be part of The Great List Purge.
  List<double> LdPoints;
  // Size of the locations array
  int iSizePoints;
  Bdry2D(const Bdry2D&);
  Bdry2D&
  operator=(const Bdry2D&);

protected:
  // Actual number of patches assigned...
  int iCount;
  // Actual number of Locations assigned..
  int iNumPointsUsed;
  // To determine if point i is a bdry point or not
  List<bool> LqBdryPoint;
  // To determine if the patches are used or not (when optimized)
  int iNum2DPatchesUsed;
  List<structVertConnect> LsVertConnectivity;
  int iNumVertsConnectFound;
  int iMaxVertsConnectFound;
  List<int> LiNumPatchesToPoint;
  List<int> LiStartVert;
  List<bool> LqVertSmallAngle;
  List<bool> LqVertUsed;
//    List<structPatchUsed> LsPatchUsed;

public:
  Bdry2D() :
      BdryRep(), iSz(0), LpBPatches(), LdPoints(), iSizePoints(0), iCount(0), iNumPointsUsed(
	  0), LqBdryPoint(), iNum2DPatchesUsed(0), LsVertConnectivity(), iNumVertsConnectFound(
	  0), iMaxVertsConnectFound(0), LiNumPatchesToPoint(), LiStartVert(), LqVertSmallAngle(), LqVertUsed()
  {
    LqVertUsed.vClear();
//        LsPatchUsed.vClear();
    iNumVertsConnectFound = 0;
  }
  Bdry2D(const char strBaseFileName[], const bool qOptimizePatches = false);
//  Bdry2D(const ECVertex ECVertIn);
  ~Bdry2D();

  virtual int
  iNumPatches() const
  {
    return LpBPatches.iLength();
  }
  int
  iNumPatchesUsed() const
  {
    return iNum2DPatchesUsed;
  }
  int
  iPatchIndex(BdryPatch2D* const pBP) const
  {
    assert(pBP);
    int iRetVal = LpBPatches.iSeek(pBP);
    assert(iRetVal >= 0);
    return iRetVal;
  }

  int
  iNumPointsInUse() const
  {
    return iNumPointsUsed;
  }
  int
  iNumPoints() const
  {
    return LdPoints.iLength() / 2;
  }
  void
  vClearPoints()
  {
    LdPoints.vClear();
  }
  void
  vSetSize(const int iNumPatchesIn, const int iNVerts)
  {

    // Must set internal size too...
    iSz = iNumPatchesIn;
    // Set counter to 0...
    iCount = 0;

    LpBPatches.vClear();
    iNum2DPatchesUsed = 0;
    iSizePoints = iNVerts;
    iNumPointsUsed = 0;
    // Setup the used/unused vertices array
    LqBdryPoint.vClear();
//      LsPatchUsed.vClear();
    LsVertConnectivity.vClear();
    iNumVertsConnectFound = 0;
    LiStartVert.vClear();
    LiNumPatchesToPoint.vClear();
    LqVertSmallAngle.vClear();

  }

  void
  vAddPatchToList(BdryPatch2D* pPatch);

  // Access functions to the used/unused states of vertices and patches
  bool
  qIsBdryPoint(int iV) const;
  void
  vSetBdryPoint(int iV, bool qBdry);
//    bool qIsPatchUsed(int iChk) const;
//    int  iPatchWasReplacedBy(int iChk);

  BdryPatch2D*
  pBdry(const int i) const
  {
    assert(i >= 0 && i < LpBPatches.iLength());
    return (LpBPatches[i]);
  }

  BdryPatch2D*
  operator[](const int i) const
  {
    return pBdry(i);
  }

  void
  vGetPoint(const int i, double adLoc[2]) const
  {
    assert(i >= 0 && i < LdPoints.iLength() / 2);
    int iBase = 2 * i;
    adLoc[0] = LdPoints[iBase];
    adLoc[1] = LdPoints[iBase + 1];
  }

  void
  vSetPoint(const int i, const double adLoc[2])
  {
    assert(i >= 0 && i < LdPoints.iLength() / 2);
    int iBase = 2 * i;
    LdPoints[iBase] = adLoc[0];
    LdPoints[iBase + 1] = adLoc[1];
  }

  void
  vAddPoint(const double adLoc[2])
  {
    LdPoints.vAppendItem(adLoc[0]);
    LdPoints.vAppendItem(adLoc[1]);
    LiNumPatchesToPoint.vAppendItem(0);
    LqVertSmallAngle.vAppendItem(false);
    LqBdryPoint.vAppendItem(false);
    LiStartVert.vAppendItem(-1);
    iNumPointsUsed++;
  }

  void
  vOrderConnectivityInfo();
  int
  iNumPatchesConnectedToPoint(int iVertIndex) const
  {
    assert((iVertIndex >= 0) && (iVertIndex < LiNumPatchesToPoint.iLength()));
    return LiNumPatchesToPoint[iVertIndex];
  }
  void
  vSetNumPatchesConnectedToPoint(int iVertIndex, int iNumPatches);

  int
  iConnectivityListLength() const
  {
    return LsVertConnectivity.iLength();
  }
  ;
  structVertConnect
  sVCConnectivityListContent(int iIndex) const
  {
    assert(iIndex < LsVertConnectivity.iLength());
    return LsVertConnectivity[iIndex];
  }
  ;
  void
  vAddToConnectivityList(structVertConnect sVC)
  {
    LsVertConnectivity.vAppendItem(sVC);
    // Increment counter
    iNumVertsConnectFound++;
  }
  ;
  BdryPatch2D*
  pB2DPatchConnectedToPoint(int iVertIndex, int iPatchIndex);
  void
  vFindSmallAngles();
  bool
  qVertSmallAngle(int iV) const
  {
    assert((iV >= 0) && (iV < LqVertSmallAngle.iLength()));
    return LqVertSmallAngle[iV];
  }
  ;
  void
  vSetSmallAngle(int iV, bool qAngle)
  {
    assert((iV >= 0) && (iV < LqVertSmallAngle.iLength()));
    LqVertSmallAngle[iV] = qAngle;
  }
  ;
  void
  vSetStartVert(int iV, int iStart)
  {
    assert((iV >= 0) && (iV < LiStartVert.iLength()));
    LiStartVert[iV] = iStart;
  }
  int
  iStartVert(int iV) const
  {
    assert((iV >= 0) && (iV < LiStartVert.iLength()));
    return LiStartVert[iV];
  }
};

#endif
