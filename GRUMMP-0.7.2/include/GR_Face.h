#ifndef GR_Face
#define GR_Face 1

#include "GR_config.h"
#include "GR_Entity.h"
/* #include "GR_Vertex.h" */
#include "GR_Cell.h"

class CubitVector;

// Categories into which tetrahedral faces are placed for purposes of
// swapping.  T indicates that the face is transformable, N indicates
// non-transformable.  For "T" faces, the digits are the number of tets
// before and after the transformation.  For "N" faces, the digits are
// the number of tets that would be present for a transformable case; a
// zero as the second digit indicates that the case is irrevocably
/// unswappable.
enum eFaceCat {
  eT23 = 0,
  eT32 = 1,
  eT22 = 2,
  eT44 = 3,
  eN32 = 4,
  eN44 = 5,
  eN40 = 6,
  eN30 = 7,
  eN20 = 8,
  eBdry = 9,
  eBdryReflex = 10,
  eOther = 11
};

///
class Face : public Entity {
protected:
  /// Copy construction disallowed.
  Face(const Face&) :
      Entity(), m_leftCell(NULL), m_rightCell(NULL), m_length(-1), m_badFace(
	  false), m_doNotDel(false), m_inUse(true), dQualInc(-1), bInsertion(
	  false), m_faceLoc(eInterior), m_restricted(false), m_sameShape(false), m_mergeable(
	  true), m_deleted(false), m_locked(false), bCombineToQuad(false)
  {
    assert(0);
  }
  ///
  Cell* m_leftCell;
  ///
  Cell* m_rightCell;
  ///
  Vert* m_verts[4];
  ///
  Face(Cell* const pCA, Cell* const pCB, const int iNV, Vert* const pV0 =
  pVInvalidVert,
       Vert* const pV1 = pVInvalidVert, Vert* const pV2 =
       pVInvalidVert,
       Vert* const pV3 = pVInvalidVert);

  double m_length;
  bool m_badFace;
  bool m_doNotDel;
  bool m_inUse;
  double dQualInc;
  bool bInsertion;

public:
  ///
  enum eFaceType {
    eEdgeFace, eMultiEdge, eTriFace, eQuadFace
  };
  ///
  // This enum must begin with eUnknown and end with eInvalidFaceLoc
  enum eFaceLoc {
    eUnknown = 0,
    eBdryFace,
    eBdryTwoSide,
    ePseudoSurface,
    eInterior,
    eExterior,
    eInvalidFaceLoc
  };
  // TODO These public variables must be fixed!
  int ordernumber = 0;
  bool Done = false;
  bool ToDo = false;
  int layerindex = 1e9;
protected:
  ///
  mutable unsigned int m_faceLoc :5;
  bool m_restricted :1, m_sameShape :1, m_mergeable :1, m_deleted :1,
      m_locked :1, bCombineToQuad :1;
public:
  ///
  virtual
  ~Face()
  {
  }
  ///
  Face&
  operator=(Face& F);
  /// Ensure a sane initial state for face data
  void
  resetAllData();
  /// Returns the type of face.
  /// TODO: Try to eliminate this function and replace switches on its return
  /// value with virtual functions.
  virtual eFaceType
  getType() const = 0;
  ///
  virtual int
  getNumVerts() const = 0;
  ///
  void
  setDefaultFlags();
  ///
  void
  copyAllFlags(const Face &F);
  ///
  int
  getFaceLoc() const;

  ///
  void
  setFaceLoc(const int newLoc)
  {
    assert(newLoc >= eUnknown && newLoc <= eInvalidFaceLoc);
    m_faceLoc = (newLoc & 0xF);
  }

  ///
  void
  setFaceLoc();
  bool
  isMergeable() const
  {
    return m_mergeable;
  }
  void
  markAsNotMergeable()
  {
    m_mergeable = false;
  }
  void
  markAsMergeable()
  {
    m_mergeable = true;
  }
  ///
  bool
  isDeleted() const
  {
    return m_deleted;
  }

  ///
  bool
  isNotDeleted() const
  {
    return !m_deleted;
  }

  ///
  void
  markAsDeleted()
  {
    m_deleted = true;
  }

  ///
  void
  markAsNotDeleted()
  {
    m_deleted = false;
  }

  ///
  bool
  isLocked() const
  {
    return m_locked;
  }

  ///
  void
  lock()
  {
    m_locked = true;
  }

  ///
  void
  unlock()
  {
    m_locked = false;
  }

  /// Functions to "play" with the restricted Delaunay tag.
  void
  setRestrictedDelaunay()
  {
    m_restricted = true;
  }

  ///
  void
  unsetRestrictedDelaunay()
  {
    m_restricted = false;
  }

  /// This returns the restricted Delaunay tag.
  bool
  isRestrictedDelaunay() const
  {
    return m_restricted;
  }

  bool
  isCombineToQuad()
  {
    return bCombineToQuad;
  }

  void
  CombineToQuad(bool Merge)
  {
    bCombineToQuad = Merge;
  }

  ///
  bool
  isSwapAllowed() const;
  ///
  bool
  isBdryFace() const;

  ///
  bool
  isValid() const
  {
    return (this != pFInvalidFace);
  }

  ///
  virtual int
  doFullCheck() const;

  ///
  virtual bool
  isLocallyDelaunay() const
  {
    return false;
  }

  inline virtual bool
  hasVert(const Vert* pV) const;

  ///
  virtual bool
  hasCell(const Cell* const pC) const
  {
    return (pC == m_leftCell || pC == m_rightCell);
  }

  /// Data member reporting
  Cell*
  getLeftCell() const
  {
    assert(isValid());
    return m_leftCell;
  }

  ///
  Cell*
  getRightCell() const
  {
    assert(isValid());
    return m_rightCell;
  }

  ///
  virtual Cell*
  getOppositeCell(const Cell* const pC) const
  {
    assert(isValid());
    assert(
	m_leftCell == pC || m_rightCell == pC || (m_leftCell == pCInvalidCell && m_rightCell == pCInvalidCell));
    return ((m_leftCell == pC) ? m_rightCell : m_leftCell);
  }

  ///
  virtual int
  getNumCells() const
  {
    int iN = 0;
    if (m_leftCell->isValid() && !m_leftCell->isBdryCell())
      iN++;

    if (m_rightCell->isValid() && !m_rightCell->isBdryCell())
      iN++;

    return iN;
  }
  virtual Vert*
  getOppositeVert(Vert* pV) const
  {
    assert(isValid());
    assert(getNumVerts() == 2);
    assert(hasVert(pV));
    if (pV == m_verts[0])
      return m_verts[1];
    else
      return m_verts[0];
  }
  int
  getNumFaces() const
  {
    return 0;
  }

  const Vert*
  getVert(const int i) const
  {
    assert(isValid());
    assert(i >= 0 && i < getNumVerts());
    return m_verts[i];
  }

  Vert*
  getVert(int i)
  {
    assert(isValid());
    assert(i >= 0 && i < getNumVerts());
    return m_verts[i];
  }

  const Face*
  getFace(const int) const
  {
    return this;
  }

  Face*
  getFace(const int)
  {
    return this;
  }

  Cell*
  getCell(const int i)
  {
    assert(i == 0 || i == 1);
    if (i == 0) {
      if (m_leftCell->isValid()) {
	int eCTL = m_leftCell->getType();
	if (eCTL == Cell::eTriBFace || eCTL == Cell::eIntTriBFace
	    || eCTL == Cell::eQuadBFace || eCTL == Cell::eIntQuadBFace) {
	  return (m_rightCell);
	}
	else
	  return (m_leftCell);
      }
      else { // pCL is invalid
	return (m_rightCell);
      }
    }
    else {
      if (m_leftCell->isValid() && m_rightCell->isValid()) {
	int eCTL = m_leftCell->getType();
	int eCTR = m_rightCell->getType();
	if (eCTL == Cell::eTriBFace || eCTR == Cell::eTriBFace
	    || eCTL == Cell::eIntTriBFace || eCTR == Cell::eIntTriBFace
	    || eCTL == Cell::eQuadBFace || eCTR == Cell::eQuadBFace
	    || eCTL == Cell::eIntQuadBFace || eCTR == Cell::eIntQuadBFace) {
	  return (pCInvalidCell);
	}
	else
	  return (m_rightCell);
      }
      else
	return pCInvalidCell; // If one cell's invalid, then the
      // valid one is returned for i=0, so
      // return garbage here.
    }
    // Can never get here.
    assert(0);
    return pCInvalidCell;
  }

  const Cell*
  getCell(const int i) const
  {
    assert(i == 0 || i == 1);
    if (i == 0) {
      if (m_leftCell->isValid()) {
	int eCTL = m_leftCell->getType();
	if (eCTL == Cell::eTriBFace || eCTL == Cell::eIntTriBFace
	    || eCTL == Cell::eQuadBFace || eCTL == Cell::eIntQuadBFace) {
	  return (m_rightCell);
	}
	else
	  return (m_leftCell);
      }
      else { // pCL is invalid
	return (m_rightCell);
      }
    }
    else {
      if (m_leftCell->isValid() && m_rightCell->isValid()) {
	int eCTL = m_leftCell->getType();
	int eCTR = m_rightCell->getType();
	if (eCTL == Cell::eTriBFace || eCTR == Cell::eTriBFace
	    || eCTL == Cell::eIntTriBFace || eCTR == Cell::eIntTriBFace
	    || eCTL == Cell::eQuadBFace || eCTR == Cell::eQuadBFace
	    || eCTL == Cell::eIntQuadBFace || eCTR == Cell::eIntQuadBFace) {
	  return (pCInvalidCell);
	}
	else
	  return (m_rightCell);
      }
      else
	return pCInvalidCell; // If one cell's invalid, then the
      // valid one is returned for i=0, so
      // return garbage here.
    }
    // Can never get here.
    assert(0);
    return pCInvalidCell;
  }

  void
  getAllFaceHandles(GRUMMP_Entity*[]) const
  /*aHandle*/
  {
  }
  void
  getAllCellHandles(GRUMMP_Entity* aHandle[]) const;

/// Data member manipulation
  virtual void
  addCell(Cell* const pCNew);
  ///
  virtual void
  removeCell(const Cell* const pCOld);
  ///
  virtual void
  replaceCell(const Cell* const pCOld, Cell* const pCNew);
  ///
  void
  replaceVert(const Vert* const pVOld, Vert* const pVNew);
  ///
  void
  interchangeCells(void);
  ///
  void
  interchangeCellsAndVerts(void);
  ///
  void
  assign(Cell* const pCA, Cell* const pCB, Vert* const pV0 =
  pVInvalidVert,
	 Vert* const pV1 = pVInvalidVert, Vert* const pV2 =
	 pVInvalidVert,
	 Vert* const pV3 = pVInvalidVert);
  /// Change right cell connectivity data.  No checks are done, so
  /// caveat programmer. This function should be called only from
  /// functions defined in MeshMod2D.cxx or MeshMod3D.cxx.
  void
  setRightCell(Cell* pC)
  {
    m_rightCell = pC;
  }
  /// Change left cell connectivity data.  No checks are done, so
  /// caveat programmer. This function should be called only from
  /// functions defined in MeshMod2D.cxx or MeshMod3D.cxx.
  void
  setLeftCell(Cell* pC)
  {
    m_leftCell = pC;
  }
  void
  setVerts(Vert* pV0, Vert* pV1, Vert* pV2 = pVInvalidVert, Vert* pV3 =
  pVInvalidVert);
/// Geometric computations
  virtual void
  calcNormal(double * const adNorm) const = 0;
  ///
  virtual void
  calcUnitNormal(double * const adNorm) const = 0;
  ///
  virtual void
  projectOntoFace(double adPoint[]) const;
  ///
  virtual double
  calcSize() const = 0;
  ///
  virtual void
  calcCentroid(double adCent[]) const;

  // The following functions are used only in 2D aniso refinement.  When that
  // is rewritten as a service, these functions will almost certainly be
  // removed, if only because it's a lot of memory use for something that's
  // needed only rarely.
  double
  getLength() const
  {
    return m_length;
  }
  void
  setLength(double dLengthIn)
  {
    assert(dLengthIn > 0);
    m_length = dLengthIn;
  }

  double
  getdQualInc() const
  {
    return dQualInc;
  }
  void
  setdQualInc(double d_QualInc)
  {
    assert(d_QualInc > 0);
    dQualInc = d_QualInc;
  }

  bool
  getbInsertion() const
  {
    return bInsertion;
  }
  void
  setbInsertion(double b_Insertion)
  {
    bInsertion = b_Insertion;
  }

  // void vSetAniLength(double dLength) {assert(dLength>0);dAniLength=dLength;}
  // double dGetAniLength() const {assert(dAniLength>=0);return dAniLength;}
  bool
  isBadFace() const
  {
#pragma omp critical(setBadFace)
    {
    }
    return m_badFace;
  }
  inline bool
  isBadFaceWOC() const
  {
    return m_badFace;
  }
  void
  setBadFace(bool qInput)
  {
#pragma omp critical(setBadFace)
    {
      m_badFace = qInput;
    }
  }
  inline void
  setBadFaceWOC(bool qInput)
  {
    m_badFace = qInput;
  }
  bool
  isInUse() const
  {
#pragma omp critical(SetInUseFace)
    {
    }
    return m_inUse;
  }
  void
  setInUse(bool qInput)
  {
#pragma omp critical(SetInUseFace)
    {
      m_inUse = qInput;
    }
  }
  inline void
  setInUseWOC(bool qInput)
  {
    m_inUse = qInput;
  }
  bool
  qDoNotDelete() const
  {
    return m_doNotDel;
  }
  void
  vSetDoNotDelete(bool qDoNotDelete_)
  {
    m_doNotDel = qDoNotDelete_;
  }
  /// uses printVertInfo()
  void
  printFaceInfo() const;
};

///
class QuadFace : public Face {
  /// Copy construction disallowed.
  QuadFace(const QuadFace& QF) :
      Face(QF)
  {
    assert(0);
  }
public:
  ///
  QuadFace(Cell* const pCA = pCInvalidCell, Cell* const pCB = pCInvalidCell,
	   Vert* const pV0 = pVInvalidVert, Vert* const pV1 = pVInvalidVert,
	   Vert* const pV2 = pVInvalidVert, Vert* const pV3 = pVInvalidVert) :
      Face(pCA, pCB, 4, pV0, pV1, pV2, pV3)
  {
  }
  ///
  int
  doFullCheck() const;
  ///
  void
  calcNormal(double * const adNorm) const;
  ///
  void
  calcUnitNormal(double * const adNorm) const;
  ///
  eFaceType
  getType() const
  {
    return eQuadFace;
  }
  ///
  int
  getNumVerts() const
  {
    return 4;
  }
  ///
  double
  calcSize() const;

  int
  getEntType() const
  {
    return iBase_FACE;
  }
  int
  getEntTopology() const
  {
    return iMesh_QUADRILATERAL;
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandle[]) const
  {
    aHandle[0] = static_cast<Entity*>(m_verts[0]);
    aHandle[1] = m_verts[1];
    aHandle[2] = m_verts[2];
    aHandle[3] = m_verts[3];
  }
};

///
class TriFace : public Face {

  /// Copy construction disallowed.
  TriFace(const TriFace& TF) :
      Face(TF)
  {
    assert(0);
  }
public:
  ///
  TriFace(Cell* const pCA = pCInvalidCell, Cell* const pCB = pCInvalidCell,
	  Vert* const pV0 = pVInvalidVert, Vert* const pV1 = pVInvalidVert,
	  Vert* const pV2 = pVInvalidVert) :
      Face(pCA, pCB, 3, pV0, pV1, pV2)
  {
  }
  ///
  int
  doFullCheck() const;
  ///
  void
  calcNormal(double * const adNorm) const;
  ///
  void
  calcUnitNormal(double * const adNorm) const;
  ///
  void
  calcCircumcenter(double adCircCent[]) const;
  ///SG added the following.
  CubitVector
  calcCircumcenter() const;
  ///
  double
  calcCircumradius() const;
  ///
  double
  calcSize() const;

  ///
  eFaceType
  getType() const
  {
    return eTriFace;
  }
  ///
  int
  getNumVerts() const
  {
    return 3;
  }
  ///
  inline virtual bool
  hasVert(const Vert* const pV) const;
  ///
  eFaceCat
  categorizeFace(Vert*& pVVertA, Vert*& pVVertB, Vert*& pVVertC, Vert*& pVVertD,
		 Vert*& pVVertE, TetCell* apTCTets[], int& iNTets,
		 Vert*& pVPivot0, Vert*& pVPivot1, Vert*& pVOther,
		 const double dMaxAngle = 0) const;
  ///
  bool
  isLocallyDelaunay() const;
  int
  getEntType() const
  {
    return iBase_FACE;
  }
  int
  getEntTopology() const
  {
    return iMesh_TRIANGLE;
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandle[]) const
  {
    aHandle[0] = m_verts[0];
    aHandle[1] = m_verts[1];
    aHandle[2] = m_verts[2];
  }
  bool
  isRightHanded(const Vert * const pV0, const Vert * const pV1,
		const Vert * const pV2)
  {
    if (pV0 == m_verts[0]) {
      assert(
	  (pV1 == m_verts[1] && pV2 == m_verts[2])
	      || (pV2 == m_verts[1] && pV1 == m_verts[2]));
      return (pV1 == m_verts[1] && pV2 == m_verts[2]);
    }
    else if (pV0 == m_verts[1]) {
      assert(
	  (pV1 == m_verts[2] && pV2 == m_verts[0])
	      || (pV2 == m_verts[2] && pV1 == m_verts[0]));
      return (pV1 == m_verts[2] && pV2 == m_verts[0]);
    }
    else {
      assert(pV0 == m_verts[2]);
      assert(
	  (pV1 == m_verts[0] && pV2 == m_verts[1])
	      || (pV2 == m_verts[0] && pV1 == m_verts[1]));
      return (pV1 == m_verts[0] && pV2 == m_verts[1]);
    }
  }
};

///
class EdgeFace : public Face {

  /// Copy construction disallowed.
protected:
  EdgeFace(const EdgeFace& EF) :
      Face(EF)
  {
    assert(0);
  }

public:
  ///
  EdgeFace(Cell* const pCA = pCInvalidCell, Cell* const pCB = pCInvalidCell,
	   Vert* const pV0 = pVInvalidVert, Vert* const pV1 = pVInvalidVert) :
      Face(pCA, pCB, 2, pV0, pV1)
  {
  }
  ///
  int
  doFullCheck() const;
  ///
  void
  calcNormal(double * const adNorm) const;
  ///
  void
  calcUnitNormal(double * const adNorm) const;
  ///
  double
  calcSize() const;
  ///
  eFaceType
  getType() const
  {
    return eEdgeFace;
  }
  ///
  int
  getNumVerts() const
  {
    return 2;
  }
  ///
  bool
  isLocallyDelaunay() const;
  ///
  void
  calcCentroid(double adLoc[]) const;

  int
  getEntType() const
  {
    return iBase_EDGE;
  }
  int
  getEntTopology() const
  {
    return iMesh_LINE_SEGMENT;
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandle[]) const
  {
    aHandle[0] = m_verts[0];
    aHandle[1] = m_verts[1];
  }
};

class MultiEdge : public EdgeFace {
  Cell *apCExtra[8];
  int iNCells;
  MultiEdge(MultiEdge& ME) :
      EdgeFace(ME), iNCells(0)
  {
    assert(0);
  }
public:
  MultiEdge&
  operator=(MultiEdge& ME)
  {
    if (this != &ME) {
      this->EdgeFace::operator=(ME);
      iNCells = ME.iNCells;
      for (int ii = 0; ii < iNCells; ii++) {
	apCExtra[ii] = ME.apCExtra[ii];
      }
    }
    return (*this);
  }
public:
  MultiEdge() :
      EdgeFace(), iNCells(0)
  {
  }
  // This constructor is actually never used.
//   MultiEdge(Vert* const pV0, Vert* const pV1, Cell **apCExtraIn,
//	    int iNCellsIn)
//     : EdgeFace(apCExtraIn[0], apCExtraIn[1], pV0, pV1), iNCells(iNCellsIn)
//     {
//       assert(iNCells > 2);
//       for (int i = 2; i < iNCells; i++)
//	apCExtra[i-2] = apCExtraIn[i];
//     }
  void
  vAssign(Vert* const pV0, Vert* const pV1, Cell **apCExtraIn, int iNCellsIn);
  void
  addCell(Cell* const pCNew);
  void
  removeCell(const Cell* const pCOld);
  void
  replaceCell(const Cell* const pCOld, Cell* const pCNew);
  const Cell *
  getCell(const int i) const
  {
    assert(i >= 0 && i < iNCells);
    switch (i)
      {
      case 0:
	return m_leftCell;
      case 1:
	return m_rightCell;
      default:
	return apCExtra[i - 2];
      }
  }
  Cell *
  getCell(const int i)
  {
    return const_cast<Cell*>(static_cast<const Face*>(this)->getCell(i));
  }
  bool
  hasCell(const Cell* const pC) const;
  int
  getNumCells() const
  {
    return iNCells;
  }
  eFaceType
  getType() const
  {
    return eMultiEdge;
  }
  Cell*
  getOppositeCell(const Cell* const pC) const;
  ///
  bool
  isLocallyDelaunay() const
  {
    assert(0);
    return false;
  } // Unsafe to use.
  int
  getEntType() const
  {
    return iBase_EDGE;
  }
  int
  getEntTopology() const
  {
    return iMesh_LINE_SEGMENT;
  }
};

inline int
Face::getFaceLoc() const
{
  return m_faceLoc;
}

inline bool
Face::hasVert(const Vert* pV) const
{
//  fprintf(stdout, "Inside Face::qHasVert: %d\n", this);
  assert(isValid());
  for (int i = getNumVerts() - 1; i >= 0; i--)
    if (pV == m_verts[i])
      return (true);
  return (false);
}

inline bool
TriFace::hasVert(const Vert* const pV) const
{

  assert(isValid());

  return (pV == m_verts[0] ? true : pV == m_verts[1] ? true :
	  pV == m_verts[2] ? true : false);

}

///
Cell*
findCommonCell(const Face* const pF0, const Face* const pF1);
Vert*
findCommonVert(Face* pF0, Face* pF1);
Vert*
findCommonVert(const Face* const pF0, const Face* const pF1,
	       const Face* const pF2);

#endif
