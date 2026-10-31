#ifndef GR_Constrain
#define GR_Constrain 1
#include "GR_config.h"

class ConstrainEdge {
public:
  int iV0, iV1;
  ConstrainEdge(const int i0 = -1, const int i1 = -1) :
      iV0(min(i0, i1)), iV1(max(i0, i1))
  {
  }
  bool
  operator==(const ConstrainEdge& CE) const
  {
    return (iV0 == CE.iV0 && iV1 == CE.iV1);
  }
  bool
  operator<(const ConstrainEdge& CE) const
  {
    return ((iV0 < CE.iV0) || ((iV0 == CE.iV0) && (iV1 < CE.iV1)));
  }
  bool
  operator!=(const ConstrainEdge& CE) const
  {
    return (!(*this == CE));
  }
};

class ConstrainFace {
public:
  int iV0, iV1, iV2;
  ConstrainFace(const int i0 = -1, const int i1 = -1, const int i2 = -1) :
      iV0(i0), iV1(i1), iV2(i2)
  {
    // Must put the entries in order so that sorting will work.  The
    // following three swaps will always order these three (although not
    // all may be needed.
    if (iV0 > iV1) {
      int iTmp = iV1;
      iV1 = iV0;
      iV0 = iTmp;
    }
    if (iV1 > iV2) {
      int iTmp = iV1;
      iV1 = iV2;
      iV2 = iTmp;
    }
    if (iV0 > iV1) {
      int iTmp = iV1;
      iV1 = iV0;
      iV0 = iTmp;
    }
    assert(
	(iV0 < iV1 && iV1 < iV2)
	    || ((iV0 == -1) && (iV1 == -1) && (iV2 == -1)));
  }
  bool
  operator==(const ConstrainFace& CF) const
  {
    return (iV0 == CF.iV0 && iV1 == CF.iV1 && iV2 == CF.iV2);
  }
  bool
  operator<(const ConstrainFace& CF) const
  {
    return ((iV0 < CF.iV0) || ((iV0 == CF.iV0) && (iV1 < CF.iV1))
	|| ((iV0 == CF.iV0) && (iV1 == CF.iV1) && (iV2 < CF.iV2)));
  }
  bool
  operator!=(const ConstrainFace& CF) const
  {
    return (!(*this == CF));
  }
};

class BadEdge {
public:
  // Default operator= and copy constructor are okay here, because the
  // pointer isn't allocated in this class anyway; the data is owned
  // elsewhere. 
  Vert *pVN, *pVS, *pVO;
  Face *pF;
  BadEdge(Vert * const pVNorth = pVInvalidVert, Vert * const pVSouth =
  pVInvalidVert,
	  Vert * const pVOther = pVInvalidVert, Face * const pFace =
	  pFInvalidFace) :
      pVN(pVNorth), pVS(pVSouth), pVO(pVOther), pF(pFace)
  {
  }
  BadEdge(const BadEdge& BE) :
      pVN(BE.pVN), pVS(BE.pVS), pVO(BE.pVO), pF(BE.pF)
  {
  }
  BadEdge&
  operator=(const BadEdge& BE)
  {

    if (this != &BE) {
      pVN = BE.pVN;
      pVS = BE.pVS;
      pVO = BE.pVO;
      pF = BE.pF;
    }
    return (*this);
  }
  int
  operator==(const BadEdge& BE) const
  {
    return (pVN == BE.pVN && pVS == BE.pVS);
  }
  int
  operator<(const BadEdge& BE) const
  {
    return ((pVN < BE.pVN) || ((pVN == BE.pVN) && (pVS < BE.pVS)));
  }
  int
  operator!=(const BadEdge& BE) const
  {
    return (!(*this == BE));
  }
};

#endif

