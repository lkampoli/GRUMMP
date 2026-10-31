#ifndef GR_Sort
#define GR_Sort 1

// Define several sortable classes that are each useful in several
// places in GRUMMP.  Also, the STL set definition is included, so that
// all of the power and convenience of that implementation can be used
// to simplify coding later on.

#include <algorithm>
using std::sort;
#include <set>
using std::multiset;
using std::set;
using std::pair;

#include "GR_config.h"
#include "GR_Vertex.h"

// All sortable classes must define operator= and a comparison operator
// (operator<) which allows ties only for truly duplicate entries, at
// least, in addition to a constructor and access functions.

// This class is used when determining heuristically which vertices in a
// mesh are of type BdryApex and BdryCurve (both two and three dimensions)
class BdryNorm {
private:
  Vert *pV;
  double adNorm[3];
public:
  BdryNorm() :
      pV(pVInvalidVert)
  {
  }
  BdryNorm(const BdryNorm& BN) :
      pV(BN.pV)
  {
    adNorm[0] = BN.adNorm[0];
    adNorm[1] = BN.adNorm[1];
    adNorm[2] = BN.adNorm[2];
  }
  BdryNorm(Vert * const pVIn, double adN[]) :
      pV(pVIn)
  {
    assert(pVIn->isValid());
    adNorm[XDIR] = adN[XDIR];
    adNorm[YDIR] = adN[YDIR];
    if (pV->getSpaceDimen() == 3)
      adNorm[ZDIR] = adN[ZDIR];
    else
      adNorm[ZDIR] = 0;
  }
  BdryNorm&
  operator=(const BdryNorm& BN)
  {
    if (this != &BN) {
      pV = BN.pV;
      adNorm[XDIR] = BN.adNorm[XDIR];
      adNorm[YDIR] = BN.adNorm[YDIR];
      adNorm[ZDIR] = BN.adNorm[ZDIR];
    }
    return (*this);
  }
  Vert *
  pVVert() const
  {
    return pV;
  }
  const double *
  adNormal() const
  {
    return adNorm;
  }
  // These are intended for use in sorting by vertex only, so they don't
  // need to compare the normal vectors themselves.
  bool
  operator<(const BdryNorm& BN) const
  {
    return (pV < BN.pV);
  }
  bool
  operator==(const BdryNorm& BN) const
  {
    return (pV == BN.pV);
  }
};
// End of class BdryNorm

// The following class is used primarily for building face-cell
// connectivity from cell-vertex data.  Can be used for both 2D and 3D
// applications.
class FaceSort {
private:
  Vert *apV[4];
  Face *pF; // This data is used in 3D Watson insertion.
  int rank[4];
  int iCLeft, iFaceInCell, iMax;
public:
  FaceSort() :
      pF(pFInvalidFace), iCLeft(-1), iFaceInCell(-1), iMax(-1)
  {
    apV[0] = apV[1] = apV[2] = apV[3] = pVInvalidVert;
  }
  FaceSort(Vert * const pV0, Vert * const pV1, const int iCIn,
	   Face * const pFIn = pFInvalidFace) :
      pF(pFIn), iCLeft(iCIn), iFaceInCell(-1), iMax(2)
  {
    apV[0] = pV0;
    apV[1] = pV1;
    apV[2] = apV[3] = pVInvalidVert;
    vRankTwoWay();
  }
  FaceSort(Vert * const pV0, Vert * const pV1, Vert * const pV2, const int iCIn) :
      pF(pFInvalidFace), iCLeft(iCIn), iFaceInCell(-1), iMax(3)
  {
    apV[0] = pV0;
    apV[1] = pV1;
    apV[2] = pV2;
    apV[3] = pVInvalidVert;
    vRankThreeWay();
  }
  FaceSort(Vert * const pV0, Vert * const pV1, Vert * const pV2,
	   Vert * const pV3, const int iCIn, const int iFICIn) :
      pF(pFInvalidFace), iCLeft(iCIn), iFaceInCell(iFICIn), iMax(4)
  {
    apV[0] = pV0;
    apV[1] = pV1;
    apV[2] = pV2;
    apV[3] = pV3;
    for (iMax = 4; iMax >= 1; iMax--) {
      if (apV[iMax - 1]->isValid())
	break;
    }
    assert(iMax >= 2);
    switch (iMax)
      {
      case 2:
	vRankTwoWay();
	break;
      case 3:
	vRankThreeWay();
	break;
      case 4:
	vRankFourWay();
	break;
      }
  }
private:
  void
  vRankTwoWay()
  {
    if (apV[0] < apV[1]) {
      rank[0] = 0;
      rank[1] = 1;
    }
    else {
      rank[0] = 1;
      rank[1] = 0;
    }
  }
  void
  vRankThreeWay()
  {
    int iWinner = (apV[0] < apV[1]) ? 0 : 1;
    int iLoser = 1 - iWinner;
    if (apV[2] < apV[iWinner]) {
      rank[0] = 2; // smallest
      rank[1] = iWinner;
      rank[2] = iLoser;
    }
    else {
      rank[0] = iWinner;
      if (apV[2] < apV[iLoser]) {
	rank[1] = 2;
	rank[2] = iLoser;
      }
      else {
	rank[1] = iLoser;
	rank[2] = 2;
      }
    }
  }
  void
  vRankFourWay()
  {
    int iWinner01 = (apV[0] < apV[1]) ? 0 : 1;
    int iLoser01 = 1 - iWinner01;
    int iWinner23 = (apV[2] < apV[3]) ? 2 : 3;
    int iLoser23 = 5 - iWinner23;

    int iWL, iLW;
    if (apV[iWinner01] < apV[iWinner23]) {
      rank[0] = iWinner01;
      iWL = iWinner23;
    }
    else {
      rank[0] = iWinner23;
      iWL = iWinner01;
    }

    if (apV[iLoser01] < apV[iLoser23]) {
      rank[3] = iLoser23;
      iLW = iLoser01;
    }
    else {
      rank[3] = iLoser01;
      iLW = iLoser23;
    }

    if (apV[iLW] < apV[iWL]) {
      rank[1] = iLW;
      rank[2] = iWL;
    }
    else {
      rank[1] = iWL;
      rank[2] = iLW;
    }
  }
public:
  FaceSort&
  operator=(const FaceSort& FS)
  {
    if (this != &FS) {
      iMax = FS.iMax;
      for (int i = 0; i < iMax; i++) {
	apV[i] = FS.apV[i];
	rank[i] = FS.rank[i];
      }
      pF = FS.pF;
      iCLeft = FS.iCLeft;
      iFaceInCell = FS.iFaceInCell;
    }
    return (*this);
  }
  Vert *
  getVert(const int i) const
  {
    assert(i >= 0 && i < iMax);
    return apV[i];
  }
  Face *
  getFace() const
  {
    return pF;
  }
  int
  iLeftCell() const
  {
    return iCLeft;
  }
  int
  iWhichFaceInCell() const
  {
    return iFaceInCell;
  }
  int
  iRank(int i) const
  {
    assert(i >= 0 && i < iMax);
    return rank[i];
  }
  bool
  operator<(const FaceSort& FS) const
  {
    for (int ii = 0; ii < iMax; ii++) {
      if (apV[rank[ii]] != FS.apV[FS.rank[ii]]) {
	return (apV[rank[ii]] < FS.apV[FS.rank[ii]]);
      }
    }
    // The old slower (?) way
    /*       Vert *apV1[] = {apV[0], apV[1], apV[2], apV[3]}; */
    /*       Vert *apV2[] = {FS.apV[0], FS.apV[1], FS.apV[2], FS.apV[3]}; */
    /*       sort(apV1, apV1+iMax); */
    /*       sort(apV2, apV2+iMax); */

    /*       for (int ii = 0; ii < iMax; ii++) { */
    /* 	if (apV1[ii] != apV2[ii]) { */
    /* 	  return (apV1[ii] < apV2[ii]); */
    /* 	} */
    /*       } */
    // All verts are matched up.  Whether this next line has < or >
    // alters whether BFaces are on the left or the right of their faces.
    return (iCLeft > FS.iCLeft);
  }
  bool
  operator==(const FaceSort& FS) const
  {
    for (int ii = 0; ii < iMax; ii++) {
      if (apV[rank[ii]] != FS.apV[FS.rank[ii]]) {
	return false;
      }
    }
    // The old slower (?) way
    /*       Vert *apV1[] = {apV[0], apV[1], apV[2], apV[3]}; */
    /*       Vert *apV2[] = {FS.apV[0], FS.apV[1], FS.apV[2], FS.apV[3]}; */
    /*       sort(apV1, apV1+iMax); */
    /*       sort(apV2, apV2+iMax); */

    /*       for (int ii = 0; ii < iMax; ii++) { */
    /* 	if (apV1[ii] != apV2[ii]) { */
    /* 	  return false; */
    /* 	} */
    /*       } */
    // All verts are matched up.
    return (iCLeft == FS.iCLeft);
  }
  friend bool
  qIdenticalVerts(const FaceSort& FS1, const FaceSort& FS2)
  {
    for (int ii = 0; ii < FS1.iMax; ii++) {
      if (FS1.apV[FS1.rank[ii]] != FS2.apV[FS2.rank[ii]]) {
	return false;
      }
    }
    /*       Vert *apV1[] = {FS1.apV[0], FS1.apV[1], FS1.apV[2], FS1.apV[3]}; */
    /*       Vert *apV2[] = {FS2.apV[0], FS2.apV[1], FS2.apV[2], FS2.apV[3]}; */
    /*       sort(apV1, apV1+FS1.iMax); */
    /*       sort(apV2, apV2+FS2.iMax); */

    /*       for (int ii = 0; ii < FS1.iMax; ii++) { */
    /* 	if (apV1[ii] != apV2[ii]) { */
    /* 	  return false; */
    /* 	} */
    /*       } */
    // All verts are matched up.
    return (true);
  }
};
// End of class FaceSort

#endif
