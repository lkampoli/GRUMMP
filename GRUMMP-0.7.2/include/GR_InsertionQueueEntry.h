#ifndef GR_InsertionQueueEntry
#define GR_InsertionQueueEntry 1

#include "GR_config.h"
#include "GR_QualMeasure.h"

class InsertionQueueEntry {
public:
  enum eEntryType {
    eInvalid,
    eTetCell,
    eTriBFace,
    eSubSeg,
    eTriCell,
    eBdryEdge,
    eFace,
    e2DFace,
    e3DFace
  };

  // Note that multimaps with double keys sort from smallest to largest
  // by default, so big negative values go to the the head of the queue.
  enum dPriorityValues {
    dSubSegPriority = INT_MIN,
    dSubSegOffsetPriority = INT_MIN / 2,
    dBdryFacePriority = INT_MIN / 4,
    dBdryFaceOffsetPriority = INT_MIN / 8,
    dWorstCellPriority = INT_MIN / 16,
    dNeverQueuedPriority = INT_MAX,
    dInvalidPriority = INT_MAX
  };
private:

  mutable double dPri;
  double adNewPtLoc[3];
  union {
    Cell *pC;
    BFace *pBF; // Actually derived from Cell, so unnecessary, but clearer.
    Face *pF;
  };
  Vert *apV[4];
  enum eEntryType eET;
  int filterPT; // This variable stores the length scale ratio "sigma = circumradius/LeafSideSize"

  // Because bit-wise copy is correct for this class, the default
  // operator= and copy constructor don't need to be overridden.
public:
  InsertionQueueEntry() :
      dPri(dInvalidPriority), pC(pCInvalidCell), eET(eInvalid), filterPT(-1)
  {
    // An obviously bogus location value...
    adNewPtLoc[0] = adNewPtLoc[1] = adNewPtLoc[2] = -1e300;

  }
  InsertionQueueEntry(Cell * const);
  InsertionQueueEntry(BFace * const);
  InsertionQueueEntry(Face * const);
  InsertionQueueEntry(Vert * const, Vert * const);
//  double dEvaluatePriority() const;

  bool
  qInsertOffsetPoint() const
  {
    return ((dPri >= dSubSegOffsetPriority && dPri < dBdryFacePriority)
	|| (dPri >= dBdryFaceOffsetPriority && dPri < dWorstCellPriority));
  }

  void
  vsetNewCoords(double Coord, int i)
  {
    adNewPtLoc[i] = Coord;
  }

  double
  vgetNewCoords(int i)
  {
    return adNewPtLoc[i];
  }

  void
  vSetPointLocation(double adLoc[])
  {
    adNewPtLoc[0] = adLoc[0];
    adNewPtLoc[1] = adLoc[1];
    // The last component may be garbage, which is okay, but it may
    // also be illegal to read that location, which is not.
    if (eET == eTetCell || eET == eTriBFace || eET == eSubSeg)
      adNewPtLoc[2] = adLoc[2];
  }
  /// This function and the one after that are for the parallel insertion
  /// they deal with the length scale filter variable for quad tree.
  void
  updateFilter(int filter)
  {
    filterPT = filter;
  }
  int
  getFilter()
  {
    return filterPT;
  }
  const double*
  adPointLocation() const
  {
    return adNewPtLoc;
  }
  Cell *
  pCCell() const
  {
    assert(eET == eTetCell || eET == eTriCell);
    return pC;
  }
  Face *
  pFFace() const
  {
    assert(eET == e3DFace || eET == e2DFace);
    return pF;
  }
  BFace *
  pBFBFace() const
  {
    assert(eET == eTriBFace || eET == eBdryEdge);
    return pBF;
  }
  Vert *
  pVVert(const int i) const
  {
    assert(i >= 0 && i <= 3);
    assert(
	((eET == eSubSeg || eET == eBdryEdge || eET == e2DFace) && i <= 1)
	    || ((eET == eTriBFace || eET == eTriCell || eET == e3DFace)
		&& i <= 2) || (eET == eTetCell));
    return (apV[i]);
  }
  enum eEntryType
  eType() const
  {
    return eET;
  }
  double
  dPriority() const
  {
    return dPri;
  }
  void
  setPriority(const double dPrio) const
  {
    dPri = dPrio;
  }
  bool
  qStillInMesh() const;
  bool
  operator==(const InsertionQueueEntry& IQE) const;
  bool
  operator!=(const InsertionQueueEntry& IQE)
  {
    return (!(*this == IQE));
  }
  friend bool
  operator<(const InsertionQueueEntry& IQE1, const InsertionQueueEntry& IQE2);

};

#endif

