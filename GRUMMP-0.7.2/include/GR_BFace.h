#ifndef GR_BFace
#define GR_BFace 1

class BFace;

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Cell.h"

#include <set>

class BasicTopologyEntity;
class GeometryEntity;
class CubitBox;

enum eEncroachType {
  eBall, eLens, eNewLens, eNoEncroach
};
enum eEncroachResult {
  eClean, eSplit, eOffsetPoint
};

// A boundary face is a Cell so that a Face can point to it just like
// any other Cell.  This gives a Face access to its BC if that's needed,
// for example, while making all Faces look exactly the same for stuff
// other than BC's.  At the same time, BFace's can walk right into the
// mesh by way of the actual Face to which they're attached.  The Cell
// dSize and vAllDihed member functions are also useful here; the former
// is overloaded to evaluated surface area of a BFace and the latter to
// evaluate angles in the surface mesh.

class BFace : public Cell {
protected:
  // Some data required for point insertion at/near boundaries.
  mutable bool m_offsetInsRequested :1, m_neverDoOffsetIns :1;
  mutable double m_insertionLoc[3];

  // In case of emergency (like reading a file with no geometry
  // present), a place to stash the associated bdry condition.
  int m_bdryCond;

  // Copy construction disallowed, but must be declared protected so
  // that derived classes can pretend to use it.
  BFace(const BFace&);
  // Called only by more-derived classes; all args provided.
  BFace(const bool qIsInternal, Face* const pFLeft, Face* const pFRight) :
      Cell(qIsInternal, pFLeft, pFRight), m_offsetInsRequested(false), m_neverDoOffsetIns(
	  false), m_bdryCond(iDefaultBC)
  {
    m_insertionLoc[0] = m_insertionLoc[1] = m_insertionLoc[2] = 1.e300;
  }
private:
  // These are functions needed for TSTT support, and TSTT doesn't
  // recognize the existence of such a thing as BFace's, so these can
  // all be private and garbage.
  int
  getNumCells() const
  {
    assert(0);
    return 1;
  }
  int
  getEntType() const
  {
    return iBase_ALL_TYPES;
  }
  int
  getEntTopology() const
  {
    return iMesh_ALL_TOPOLOGIES;
  }

  void
  getAllVertHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  void
  getAllFaceHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  void
  getAllCellHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
public:
  // The destructor has nothing to do.
  virtual
  ~BFace()
  {
  }
  // Necessary when (for example) removing deleted BFace's from EntContainer's
  BFace&
  operator=(const BFace& BF);
  // This is implemented differently in 2D and 3D.
  virtual int
  getBdryCondition() const = 0;

  // The new interface to access boundary conditions (Sept. 11 2006)
//	virtual int getBdryCondition(const bool left) const = 0;

  void
  setBdryCondition(const int bc)
  {
    assert(bc >= 0);
    m_bdryCond = bc;
  }
  // Is there any point in the mesh that encroaches on this BFace?
  virtual eEncroachResult
  isEncroached(const eEncroachType eET = eLens) const = 0;
  eEncroachResult
  isPointEncroaching(const double adPoint[], const eEncroachType eET = eLens,
		     const bool qTieBreaker = false) const
  {
    switch (eET)
      {
      case eBall:
	return isPointEncroachingBall(adPoint, qTieBreaker);
	break;
      case eLens:
	return isPointEncroachingLens(adPoint, qTieBreaker);
	break;
      case eNewLens:
	return isPointEncroachingNewLens(adPoint, qTieBreaker);
	break;
      default:
	assert(0);
	return (eClean);
      }
  }

  virtual bool
  isCoplanarWith(const double pointLoc[]) const = 0;
  virtual eEncroachResult
  isPointEncroachingBall(const double adPoint[],
			 const bool qTieBreak = false) const = 0;
  virtual eEncroachResult
  isPointEncroachingLens(const double adPoint[],
			 const bool qTieBreak = false) const = 0;
  virtual eEncroachResult
  isPointEncroachingNewLens(const double adPoint[],
			    const bool qTieBreak = false) const = 0;
  virtual double
  calcInsertionPriority() const
  {
    assert(0);
    return 0;
  }
  virtual void
  findOffsetInsertionPoint(const double[], double[]) const
  {
    assert(0);
  }
  const double*
  getInsertionLocation() const
  {
    return m_insertionLoc;
  }
  virtual void
  setInsertionLocation(const double[]) const
  {
    assert(0);
  }
  bool
  isOffsetInsertionRequested() const
  {
    return m_offsetInsRequested;
  }
  void
  requestOffsetInsertion() const
  {
    m_offsetInsRequested = true;
  }
  void
  clearOffsetInsertionRequest() const
  {
    m_offsetInsRequested = false;
  }
  bool
  isOffsetInsertionAllowed() const
  {
    return !m_neverDoOffsetIns;
  }
  void
  forbidOffsetInsertion() const
  {
    m_neverDoOffsetIns = true;
    m_offsetInsRequested = false;
  }
  void
  allowOffsetInsertion() const
  {
    m_neverDoOffsetIns = false;
  }
  virtual void
  calcCircumcenter(double[]) const = 0;
  // Closure makes no sense for BFace's
  bool
  isClosed() const
  {
    return true;
  }
  // It is a bdry face/cell...
  virtual bool
  isBdryCell() const
  {
    return true;
  }
  Face *
  getFace(const int i = 0)
  {
    assert(i < getNumFaces());
    return (m_faces[i]);
  }
  const Face *
  getFace(const int i = 0) const
  {
    assert(i < getNumFaces());
    return (m_faces[i]);
  }
  const Vert*
  getVert(const int i) const;
  Vert*
  getVert(const int i)
  {
    return const_cast<Vert*>(static_cast<const BFace*>(this)->getVert(i));
  }
  double
  calcSize() const;
  void
  calcVecSize(double adRes[]) const;
  void
  calcUnitNormal(double adRes[]) const;
  virtual void
  setPatch(BdryPatch* const pBP) = 0;
  virtual BdryPatch*
  getPatchPointer() const = 0;
  virtual void
  calcCentroid(double[]) const
  {
    assert(0);
  }

  //Split point for Delaunay refinement.
  virtual void
  calcSplitPoint(double*) const
  {
    assert(0);
  }

  //Computes and returns the face's bounding box.
  virtual CubitBox
  bounding_box() const = 0;

  //Computes the closest point on the BFace 
  //(projected to the attached surface if it is a curved surface)
  //Not implemented for QuadBFaceBase for now.
  virtual void
  findClosestPointOnBFace(const double pt[3], double close[3])
  {
    assert(0);
    assert(pt);
    close[0] = close[1] = close[2] = -LARGE_DBL;
  }

  //FIX ME: SG 09/2007.
  //The following was added to store geometry definition with the boundary
  //face. Initially, the 2D code was designed to store these as
  //GeometryEntity pointers. The 3D tetrahedral meshing code uses
  //BasicTopologyEntity. Given the difference, I am not sure the pure
  //virtual calls should remain, especially if they only apply to one
  //particular derived class...

  //Sets and returns the geometric entity associated with this boundary face.
  virtual void
  setGeometry(GeometryEntity* const geom_ent) = 0;
  virtual GeometryEntity*
  getGeometry() = 0;

  virtual void
  setTopologicalParent(BasicTopologyEntity* const topo_ent) = 0;
  virtual BasicTopologyEntity*
  getTopologicalParent() = 0;

};
// End of BFace base class declaration

#endif
