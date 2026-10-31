#ifndef GR_BFace3D
#define GR_BFace3D 1

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_BFace.h"
#include "GR_Bdry3D.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_InsertionQueueEntry.h"

#include "BasicTopologyEntity.hpp"
#include "CubitBox.hpp"
#include "RefFace.hpp"

class TriBFaceBase : public BFace {
private:
  BdryPatch3D* m_bdryPatch;
  RefFace* m_surface;
protected:
  /// Copy construction disallowed.  operator= inherited from BFace
  TriBFaceBase(const TriBFaceBase&);
public:
  TriBFaceBase(const bool qIsInternal, Face* const pFLeft, Face* const pFRight =
  pFInvalidFace) :
      BFace(qIsInternal, pFLeft, pFRight), m_bdryPatch(NULL), m_surface(
      NULL)
  {
  }
  TriBFaceBase&
  operator=(const TriBFaceBase& TBF)
  {
    if (this != &TBF) {
      this->BFace::operator=(TBF);
      m_bdryPatch = TBF.m_bdryPatch;
      m_surface = TBF.m_surface;
    }
    return (*this);
  }

  int
  getNumVerts() const
  {
    return 3;
  }
  // Now we know enough to give the boundary condition

  virtual int
  getBdryCondition(const bool) const
  {
    assert(0);
    return 0;
  }
  virtual int
  getBdryCondition() const
  {
    if (m_bdryPatch)
      return m_bdryPatch->iBdryCond();
    else
      return m_bdryCond;
  }

  virtual CubitBox
  bounding_box() const
  {
    assert(0);
    return CubitBox();
  }

  virtual bool
  isCoplanarWith(const double pointLoc[]) const;

  enum eEncroachResult
  isEncroached(const enum eEncroachType = eBall) const;
  // Tri specific versions of standard encroachment routines
  enum eEncroachResult
  isPointEncroachingBall(const double adPoint[],
			 const bool qTieBreak = false) const;
  // Lens checking isn't currently supported.
  enum eEncroachResult
  isPointEncroachingLens(const double adPoint[],
			 const bool qTieBreak = false) const;
  enum eEncroachResult
  isPointEncroachingNewLens(const double adPoint[],
			    const bool qTieBreak = false) const;
  double
  calcInsertionPriority() const
  // To prevent infinitely recursive insertion, don't split
  // TriBFace's that are too small relative to the local length
  // scale (at this point, 15 times smaller is counted as "too
  // small".  This is easily accomplished by giving such over-small 
  // triangles a fake priority, much higher than the cutoff for
  // splitting. 
  {
    double dActualSize = 2 / sqrt(2.)
	* (dynamic_cast<const TriFace*>(getFace()))->calcCircumradius();
    double dSizeLimit = 0;
    for (int ii = 0; ii < 3; ii++) {
      dSizeLimit += getVert(ii)->getLengthScale();
    }
    dSizeLimit /= (15 * 3); // 3 = num verts; 15 = size factor
    return (
	dActualSize > dSizeLimit ?
	    InsertionQueueEntry::dBdryFacePriority :
	    InsertionQueueEntry::dNeverQueuedPriority);
  }
  void
  findOffsetInsertionPoint(const double /*adTriggerLoc*/[],
			   double /*adOffsetLoc*/[]) const
  {
    assert(0);
  }
  void
  setInsertionLocation(const double adNewLoc[]) const
  {
    m_insertionLoc[0] = adNewLoc[0];
    m_insertionLoc[1] = adNewLoc[1];
    m_insertionLoc[2] = adNewLoc[2];
  }

  //Split point for Delaunay refinement.
  virtual void
  calcSplitPoint(double split_pt[3]) const;

  void
  calcCircumcenter(double adLoc[]) const
  {
    (dynamic_cast<const TriFace*>(getFace()))->calcCircumcenter(adLoc);
  }
  // Code needed for finding point locations for Shewchuk insertion;
  // these are all 3D-specific.
  bool
  doesPointProjectInside(const double adPoint[3]) const;
  int
  splitBdryFace(VolMesh * const pVM, const int iRecursionDepth);

  // Code for boundary connectivity in 3D; 2D is handled separately
  // because it needs to reference Bdry2D and BdryPatch2D objects.
  void
  setPatch(BdryPatch* const pBP)
  {
    setPatch(dynamic_cast<BdryPatch3D*>(pBP));
  }
  void
  setPatch(BdryPatch3D* const pBP3DIn)
  {
    assert(pBP3DIn != NULL);
    m_bdryPatch = pBP3DIn;
  }
  void
  clearGeom()
  {
    m_bdryPatch = nullptr;
    m_surface = nullptr;
  }
  BdryPatch*
  getPatchPointer() const
  {
    return m_bdryPatch;
  }
  RefFace*
  getRefFacePointer() const
  {
    return m_surface;
  }
  void
  calcCentroid(double adCent[]) const
  {
    adCent[0] = adCent[1] = adCent[2] = 0;
    for (int i = 0; i < 3; i++) {
      adCent[0] += getVert(i)->x();
      adCent[1] += getVert(i)->y();
      adCent[2] += getVert(i)->z();
    }
    adCent[0] /= 3;
    adCent[1] /= 3;
    adCent[2] /= 3;
  }

  void
  findClosestPointOnBFace(const double pt[3], double close[3]);

  virtual void
  setGeometry(GeometryEntity* const)
  {
    assert(0);
  }
  virtual GeometryEntity*
  getGeometry()
  {
    assert(0);
    return NULL;
  }

  virtual void
  setTopologicalParent(BasicTopologyEntity* const topo_ent)
  {
    RefFace* surface = (dynamic_cast<RefFace*>(topo_ent));
    m_surface = surface;
  }

  virtual BasicTopologyEntity*
  getTopologicalParent()
  {
    return m_surface;
  }

};

// Boundary edges in 3D that are true boundaries of the domain
class TriBFace : public TriBFaceBase {
private:
  // Disallow copy constructor
  TriBFace(TriBFace& TBF) :
      TriBFaceBase(TBF)
  {
    assert(0);
  }
  ;
public:
  TriBFace(Face* const pF = pFInvalidFace) :
      TriBFaceBase(false, pF)
  {
    setType(Cell::eTriBFace);
  }
  ~TriBFace()
  {
  }
  int
  getNumFaces() const
  {
    return 1;
  }
};
// End of TriBFace declaration

// Boundary edges in 3D that are true boundaries of the domain
class IntTriBFace : public TriBFaceBase {
private:
  // Disallow copy constructor
  IntTriBFace(IntTriBFace& ITBF) :
      TriBFaceBase(ITBF)
  {
    assert(0);
  }
public:
  IntTriBFace(Face* const pFLeft = pFInvalidFace, Face* const pFRight =
  pFInvalidFace) :
      TriBFaceBase(true, pFLeft, pFRight)
  {
    setType(Cell::eIntTriBFace);
  }
  ~IntTriBFace()
  {
  }
  int
  getNumFaces() const
  {
    return 2;
  }
  const Face *
  getOtherFace(const Face * const pF) const
  {
    if (pF == getFace(0))
      return getFace(1);
    else
      return getFace(0);
  }
};
// End of IntTriBFace declaration

class QuadBFaceBase : public BFace {
private:
  BdryPatch3D* m_bdryPatch;
  BasicTopologyEntity* m_surface;
protected:
  /// Copy construction disallowed.  operator= inherited from BFace
  QuadBFaceBase(const QuadBFaceBase& TBF) :
      BFace(TBF), m_bdryPatch(NULL), m_surface(NULL)
  {
    assert(0);
  }
public:
  QuadBFaceBase(const bool qIsInternal, Face* const pFLeft,
		Face* const pFRight = pFInvalidFace) :
      BFace(qIsInternal, pFLeft, pFRight), m_bdryPatch(NULL), m_surface(
      NULL)
  {
  }
  QuadBFaceBase&
  operator=(const QuadBFaceBase& QBF)
  {
    if (this != &QBF) {
      this->BFace::operator=(QBF);
      m_bdryPatch = QBF.m_bdryPatch;
      m_surface = QBF.m_surface;
    }
    return (*this);
  }
  int
  getNumVerts() const
  {
    return 4;
  }
  // Now we know enough to give the boundary condition

  virtual int
  getBdryCondition(const bool) const
  {
    assert(0);
    return 0;
  }
  virtual int
  getBdryCondition() const
  {
    return m_bdryPatch->iBdryCond();
  }

  virtual CubitBox
  bounding_box() const
  {
    assert(0);
    return CubitBox();
  }

  virtual bool
  isCoplanarWith(const double[]) const
  {
    assert(0);
    return false;
  }

  // Return value below only for the benefit of picky compilers.
  enum eEncroachResult
  isEncroached(const eEncroachType = eBall) const
  {
    assert(0);
    return eClean;
  }
  // Quad specific versions of standard encroachment routines
  enum eEncroachResult
  isPointEncroachingBall(const double[], const bool = false) const
  {
    assert(0);
    return eClean;
  }
  enum eEncroachResult
  isPointEncroachingLens(const double[], const bool = false) const
  {
    assert(0);
    return eClean;
  }
  enum eEncroachResult
  isPointEncroachingNewLens(const double[], const bool = false) const
  {
    assert(0);
    return eClean;
  }
  void
  calcCircumcenter(double[]) const
  {
    assert(0);
  }

  // Code for boundary connectivity in 3D; 2D is handled separately
  // because it needs to reference Bdry2D and BdryPatch2D objects.
  void
  setPatch(BdryPatch* const pBP)
  {
    vSetPatch(dynamic_cast<BdryPatch3D*>(pBP));
  }
  void
  vSetPatch(BdryPatch3D* const pBP3DIn)
  {
    assert(pBP3DIn != NULL);
    m_bdryPatch = pBP3DIn;
  }

  BdryPatch*
  getPatchPointer() const
  {
    return m_bdryPatch;
  }

  virtual void
  setGeometry(GeometryEntity* const)
  {
    assert(0);
  }
  virtual GeometryEntity*
  getGeometry()
  {
    assert(0);
    return NULL;
  }
  void
  clearGeom()
  {
    m_bdryPatch = nullptr;
    m_surface = nullptr;
  }

  virtual void
  setTopologicalParent(BasicTopologyEntity* const topo_ent)
  {
    RefFace* surface = (dynamic_cast<RefFace*>(topo_ent));
    assert(surface);
    m_surface = surface;
  }

  virtual BasicTopologyEntity*
  getTopologicalParent()
  {
    assert(m_surface);
    return m_surface;
  }

};
// End of QuadBFaceBase declaration

// Boundary edges in 3D that are true boundaries of the domain
class QuadBFace : public QuadBFaceBase {
private:
  // Disallow copy constructor
  QuadBFace(QuadBFace& QBF) :
      QuadBFaceBase(QBF)
  {
    assert(0);
  }
public:
  QuadBFace(Face* const pF = pFInvalidFace) :
      QuadBFaceBase(false, pF)
  {
    setType(Cell::eQuadBFace);
  }
  ~QuadBFace()
  {
  }
  int
  getNumFaces() const
  {
    return 1;
  }
};
// End of QuadBFace declaration

// Boundary edges in 3D that are true boundaries of the domain
class IntQuadBFace : public QuadBFaceBase {
private:
  // Disallow copy constructor
  IntQuadBFace(IntQuadBFace& IQBF) :
      QuadBFaceBase(IQBF)
  {
    assert(0);
  }
public:
  IntQuadBFace(Face* const pFLeft = pFInvalidFace, Face* const pFRight =
  pFInvalidFace) :
      QuadBFaceBase(true, pFLeft, pFRight)
  {
    setType(Cell::eIntQuadBFace);
  }
  ~IntQuadBFace()
  {
  }
  int
  getNumFaces() const
  {
    return 2;
  }
  const Face *
  getOtherFace(const Face * const pF) const
  {
    if (pF == getFace(0))
      return getFace(1);
    else
      return getFace(0);
  }
};
// End of IntQuadBFace declaration

#endif
