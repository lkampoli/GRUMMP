#ifndef GR_Vertex
#define GR_Vertex 1

#include <vector>

#include "GR_config.h"
#include "GR_Classes.h"
#include "GR_Entity.h"
#include "GR_Vec.h"

class BasicTopologyEntity;
class GeometryEntity;
class CubitVector;
class CubitBox;

class FaceList {
  static const int s_blockSize = 10;
  int m_numUsed;
  Face *m_faces[s_blockSize];
  FaceList *m_nextFaceList;
  FaceList(FaceList&) :
      m_numUsed(0), m_nextFaceList(NULL)
  {
    assert(0);
  }
public:
  FaceList() :
      m_numUsed(0), m_nextFaceList(NULL)
  {
  }
  ~FaceList()
  {
    if (m_nextFaceList) {
      delete m_nextFaceList;
      m_nextFaceList = NULL;
    }
  }
  bool
  addFace(Face * pF)
  {
    assert(m_numUsed <= s_blockSize);
    for (int ii = 0; ii < m_numUsed; ii++) {
      if (m_faces[ii] == pF)
	return false;
    }
    if (m_numUsed == s_blockSize) {
      if (!m_nextFaceList) {
	m_nextFaceList = new FaceList();
      }
      return (m_nextFaceList->addFace(pF));
    }
    else {
      m_faces[m_numUsed++] = pF;
    }
    return true;
  }
  bool
  removeFace(Face * pF)
  {
    int i;
    FaceList *pFLCand;
    FaceList *pFLLast;
    FaceList *pFLPrev;

    //#pragma omp threadprivate (pFLCand,pFLLast,pFLPrev)

    pFLCand = this;
    pFLLast = this;
    pFLPrev = NULL;
    assert(m_numUsed > 0);
    do {
      for (i = 0; i < m_numUsed; i++) {
	if (pFLCand->m_faces[i] == pF) {
	  // Get the last face...
	  while (pFLLast->m_nextFaceList) {
	    pFLPrev = pFLLast;
	    pFLLast = pFLLast->m_nextFaceList;
	  }
	  assert(pFLLast);
	  // Put it here instead
	  pFLCand->m_faces[i] = pFLLast->m_faces[pFLLast->m_numUsed - 1];
	  // And drop the last one
	  pFLLast->m_numUsed--;
	  assert(pFLLast->m_numUsed >= 0);
	  if (pFLLast->m_numUsed == 0 && !(pFLLast == this)) {
	    assert(pFLPrev);
	    delete pFLLast;
	    pFLPrev->m_nextFaceList = NULL;
	  }
	  return true;
	}
      }
      // If you got here, this block doesn't have it.
      pFLPrev = pFLCand;
      pFLCand = pFLCand->m_nextFaceList;
    }
    while (pFLCand);
    return false;
  }
  Face*
  getFace(int i) const
  {
    if (i >= s_blockSize)
      return m_nextFaceList->getFace(i - s_blockSize);
    else {
      if (i >= m_numUsed)
	return pFInvalidFace;
      else
	return m_faces[i];
    }
  }
  bool
  hasFace(const Face* const pF) const
  {
    int i;
    for (i = 0; i < m_numUsed; i++) {
      if (m_faces[i] == pF) {
	return true;
      }
    }
    // If you got here, this block doesn't have it.
    if (m_nextFaceList)
      return m_nextFaceList->hasFace(pF);
    else
      return false;
  }
  void
  clear()
  {
    if (m_nextFaceList)
      delete m_nextFaceList;
    m_nextFaceList = NULL;
    m_numUsed = 0;
  }
  int
  size()
  {
    if (m_nextFaceList) {
      assert(m_numUsed == s_blockSize);
      return s_blockSize + m_nextFaceList->size();
    }
    else
      return m_numUsed;
  }
public:
  FaceList&
  operator=(const FaceList& FL)
  {
    if (this != &FL) {
      m_numUsed = FL.m_numUsed;
      for (int ii = 0; ii < s_blockSize; ii++) {
	m_faces[ii] = FL.m_faces[ii];
      }
      if (m_nextFaceList != NULL) {
	delete m_nextFaceList;
	m_nextFaceList = NULL;
      }
      if (FL.m_nextFaceList != NULL) {
	m_nextFaceList = new FaceList;
	(*m_nextFaceList) = *(FL.m_nextFaceList);
      }
    }
    return *this;
  }
};

Face*
findCommonFace(const Vert * const pV0, const Vert * const pV1,
	       bool qCanBeConnectedBothSides = false);

Face*
findCommonFace(const Vert * const pV0, const Vert * const pV1,
	       const Vert * const pV2, const Vert * const pV3 = pVInvalidVert,
	       bool qCanBeConnectedBothSides = false);

class Vert : public Entity {
  ///
  double m_loc[3];
  ///
#ifndef OMIT_VERTEX_HINTS
  mutable Face *m_faceHint;
#endif
  FaceList *m_faceList;
  enum {
    eDefaultVecSize = 8
  };
  ///PUT THE FOLLOWING TWO IN AN UNION.
  GeometryEntity* parent_entity;        //Parent geometric entity (2D)
  BasicTopologyEntity* parent_topology; //Parent topological entity (3D)
  ///

  // This is a packed bitfield containing a set of flags, including 10
  // bits for user data. The mechanics of the packing are handled
  // automatically by the compiler.
  unsigned int m_dim :2, m_vertType :6;
  bool m_del :1, m_delRequested :1, m_active :1, m_struct :1, m_small :1,
      m_shell :1, m_phantom :1;

private:

  ///Length-scale for this vertex
  double m_lenScale;
  ///
  double m_metric[6];
  bool m_useMetric;
  bool m_inUse;
  int m_coarsentest;

  Vert(const Vert&); // Never defined anywhere, deliberately.

public:
  const Face *
  getFace(const int i) const
  {
    return (m_faceList->getFace(i));
  }
  Face *
  getFace(const int i)
  {
    return (m_faceList->getFace(i));
  }

private:
  Cell *
  getCell(const int)
  {
    assert(0);
    return pCInvalidCell;
  }
  const Cell *
  getCell(const int) const
  {
    assert(0);
    return pCInvalidCell;
  }

public:
  ///
  enum VertType {
    eUnknown = 0,
    eBdryApex,
    eBdryCurve,
    eBdryTwoSide,
    eBdry,
    ePseudoSurface,
    eInterior,
    eInteriorFixed,
    eBBox,
    eShockUpstream = 10,
    eShockDownstream = 11,
    eShockPerimeter = 12,
    eShockInteraction = 13,
    ePhantom = 14,
    eBdryPhantom = 15,
    eInvalid
  };
  ///
  Vert(const VertType VT = eInterior) :
      Entity(), m_faceHint(pFInvalidFace), m_faceList(NULL), parent_entity(), parent_topology(), m_dim(), m_vertType(
	  (VT & 0xF)), m_del(false), m_delRequested(false), m_active(false), m_struct(
	  false), m_small(false), m_shell(false), m_lenScale(0), m_useMetric(
	  false), m_inUse(false), m_coarsentest(0)
  // Others set by SetDefaultFlags
  {
    // By default, a vertex has eInterior type, 0 space dimensions, and
    // is smoothable.  All other flags are cleared.
    m_faceList = new FaceList();
    setDefaultFlags();
    m_loc[0] = m_loc[1] = m_loc[2] = 0.;
    clearLengthScale();
    for (int i = 0; i < 5; i++)
      m_metric[i] = 0;
  }

  ///
  ~Vert()
  {
    if (m_faceList)
      delete m_faceList;
  }
  ///
  Vert&
  operator=(const Vert& V);
  /// Ensure a sane initial state for face data
  void
  resetAllData();
public:
  // ************ Functions used in OpenMP swapping ********
  void
  isInUse(bool &IsInUse)
  {
#pragma omp critical(setInUseVert)
    {
      IsInUse = m_inUse;
    }
    //	  return m_inUse;
  }
  inline void
  isInUseWOC(bool &IsInUse)
  {
    IsInUse = m_inUse;
  }
  void
  setInUse(bool qInput)
  {
#pragma omp critical(setInUseVert)
    {
      m_inUse = qInput;
    }
  }
  inline void
  setInUseWOC(bool qInput)
  {
    m_inUse = qInput;
  }

  /* TSTT Begin TSTT Begin TSTT Begin TSTT Begin TSTT Begin TSTT Begin */
  // The following functions are very useful for GRUMMP support of TSTT.
  const Vert *
  getVert(const int /*iV*/) const
  {
    return this;
  }
  Vert *
  getVert(const int /*iV*/)
  {
    return this;
  }
  void
  getAllVertHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  void
  getAllCellHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  void
  getAllFaceHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }

  int
  getNumVerts() const
  {
    return 0;
  }
  int
  getNumFaces() const
  {
    return m_faceList->size();
  }
  int
  getNumCells() const
  {
    assert(0);
    return (0);
  }
  int
  getEntType() const
  {
    return iBase_VERTEX;
  }
  int
  getEntTopology() const
  {
    return iMesh_POINT;
  }

  ///
  Face *
  getHintFace() const;
  ///
  void
  setHintFace(Face *pF);
  ///
  Face *
  setHintFace(const Mesh* const pM);
  ///
  void
  clearHintFace()
  {
    m_faceHint = pFInvalidFace;
  }
  void
  updateHintFace();
  ///
  void
  setParentTopology(BasicTopologyEntity* const parent)
  {
    parent_topology = parent;
  }

  void
  addFace(Face * pF)
  {
    //#ifndef NDEBUG				// Reza has comment this
    //#pragma omp critical(voidaddFace)
    //{
    bool success = m_faceList->addFace(pF);
    assert(success);
    //}
    //#endif						// Reza has comment this
  }
  ///
  void
  removeFace(Face * pF)
  {
    //#ifndef NDEBUG				// Reza has comment this
    //#pragma omp critical(pgrF)
    //{
    //    bool success =
    m_faceList->removeFace(pF);
    //    assert(success);
    //}
    //#endif						// Reza has comment this
  }
  void
  clearFaceConnectivity()
  {
    //#pragma omp critical(clearFaceConnectivity)
    //	  {
    m_faceList->clear();
    //	  }
  }
  bool
  hasFace(const Face * const pF) const
  {
    return m_faceList->hasFace(pF);
  }
  ///
  BasicTopologyEntity*
  getParentTopology() const
  {
    return parent_topology;
  }
  ///
  void
  setDefaultFlags();
  ///
  void
  copyAllFlagsFrom(const Vert * const pV);
  void
  copyConnectivityFrom(const Vert& V);
  void
  copyLengthInfoFrom(const Vert& V);
  ///
  void
  setType(const unsigned int VT)
  {
    assert(VT <= eInvalid);
    m_vertType = (VT & 0xF);
  }
  ///
  unsigned int
  getVertType() const
  {
    return (m_vertType);
  }
  ///
  bool
  isBdryVert() const;
  /// Returns the number of space dimensions for the vertex (number of
  /// coordinates).  For surface meshes, this is different than the
  /// topological dimension of the mesh.
  int
  getSpaceDimen() const
  {
    return (m_dim);
  }
  ///
  void
  setSpaceDimen(const unsigned iSD)
  {
    assert(m_dim == 0);
    assert(iSD == 2 || iSD == 3);
    m_dim = (iSD & 0x3);
  }
  ///
  void
  markAsDeleted()
  {
    m_del = true;
  }
  ///
  bool
  isDeleted() const
  {
    return (m_del);
  }
  ///
  void
  markActive()
  {
    m_active = true;
  }
  ///
  void
  markInactive()
  {
    m_active = false;
  }
  ///
  bool
  isActive() const
  {
    return (m_active);
  }
  ///
  void
  markToDelete()
  {
    m_delRequested = true;
  }
  ///
  void
  markToKeep()
  {
    m_delRequested = false;
  }
  ///
  bool
  isDeletionRequested() const
  {
    return (m_delRequested);
  }
  ///
  void
  markStructured()
  {
    m_struct = true;
  }
  ///
  void
  clearStructured()
  {
    m_struct = false;
  }
  ///
  bool
  isStructured() const
  {
    return m_struct;
  }
  /// A phantom vertex isn't connected to any other entities,
  /// but is still in the database.
  bool
  isPhantom() const
  {
    return m_phantom;
  }

  void
  markAsPhantom(bool phantom)
  {
    m_phantom = phantom;
    if (m_phantom) {
      if (isBdryVert()) {
	setType(Vert::eBdryPhantom);
      }
      else {
	setType(Vert::ePhantom);
      }
    }
  }
  ///
  bool
  isSmoothable() const;
  ///
  double
  x() const
  {
    return m_loc[0];
  }
  ///
  double
  y() const
  {
    return m_loc[1];
  }
  ///
  double
  z() const
  {
    return (getSpaceDimen() == 3) ? m_loc[2] : 0;
  }
  ///
  void
  setCoords(const int iNumDim, const double adNewLoc[]);
  ///
  void
  setCoords(const char * const pcBuffer);
  ///
  void
  setCoords(const int iNumDim, const CubitVector& coord);
  ///

  //Added by SG.
  void
  setParentEntity(GeometryEntity* const parent)
  {
    parent_entity = parent;
  }
  ///
  GeometryEntity*
  getParentEntity() const
  {
    return parent_entity;
  }

  //The following four added on Sept. 11 2006 SG.//
  bool
  isSmallAngleVert() const
  {
    return m_small;
  }
  bool
  isShellVert() const
  {
    return m_shell;
  }
  void
  markAsSmallAngleVert(const bool small)
  {
    m_small = small;
  }
  void
  markAsShellVert(const bool shell)
  {
    m_shell = shell;
  }

  //The following added on Jan 24th 2007 SG //
  CubitBox
  bounding_box() const;
  const double*
  getCoords() const
  {
    assert(isValid());
    return m_loc;
  }

  bool
  isMetricUsed() const
  {
    return m_useMetric;
  }
  const double*
  getMetricArray() const
  {
    return m_metric;
  }

  bool
  isValid() const
  {
    return (this != pVInvalidVert);
  }

  //The following code has to do with length scale calculation.
  //Potential rewrite (Doug?) in a near future (Sept. 11 2006)

  double
  getLengthScale() const
  {
    return m_lenScale;
  }

  void
  setLengthScale(double dLength)
  {
    assert(dLength > 0);
    m_lenScale = dLength;
  }

  void
  clearLengthScale()
  {
    m_lenScale = LARGE_DBL;
//		m_isClosestVert = m_isSecondVert = false;
//		m_vert1 = m_vert2 = pVInvalidVert;
  }

  void
  setMetric(const double * const adMetricNew)
  {
    m_useMetric = true;
    if (m_dim == 2) {

      double dDet;
      dDet = adMetricNew[0] * adMetricNew[2] - adMetricNew[1] * adMetricNew[1];

      if (dDet <= 0 || adMetricNew[0] <= 0) {
	printf("Metric [%.3e,%.3e,%.3e] is not positive definite\n",
	       adMetricNew[0], adMetricNew[1], adMetricNew[2]);
	assert(0);
      }

      for (int i = 0; i < 3; i++) {
	m_metric[i] = adMetricNew[i];
      }
    }

    if (m_dim == 3) {

      double dDet3D;
      dDet3D = adMetricNew[0]
	  * (adMetricNew[2] * adMetricNew[5] - adMetricNew[4] * adMetricNew[4])
	  - adMetricNew[1]
	      * (adMetricNew[1] * adMetricNew[5]
		  - adMetricNew[4] * adMetricNew[3])
	  + adMetricNew[3]
	      * (adMetricNew[1] * adMetricNew[4]
		  - adMetricNew[2] * adMetricNew[3]);

      if (dDet3D <= 0) {
	printf(
	    "Metric [%.3e,%.3e,%.3e,%.3e,%.3e,%.3e] is not positive definite = %.3e \n",
	    adMetricNew[0], adMetricNew[1], adMetricNew[2], adMetricNew[3],
	    adMetricNew[4], adMetricNew[5], dDet3D);
	assert(0);
      }

      for (int i = 0; i < 6; i++) {
	m_metric[i] = adMetricNew[i];
      }
    }
    return;
  }
  double
  getMetric(int i) const
  {
    return m_metric[i];
  }

  ///
  friend Face*
  findCommonFace(const Vert * const pV0, const Vert * const pV1,
		 bool qCanBeConnectedBothSides);

  friend Face*
  findCommonFace(const Vert * const pV0, const Vert * const pV1,
		 const Vert * const pV2, const Vert * const pV3,
		 bool qCanBeConnectedBothSides);

  void
  printVertInfo() const;
  void
  copyCoordsFrom(const Vert& V);
};

double
calcDistanceBetween(const Vert * const pV1, const Vert * const pV2);

#endif
