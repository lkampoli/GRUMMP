#ifndef GR_BFACE2D_H
#define GR_BFACE2D_H

#include "GR_config.h"
#include "GR_BFace.h"
#include "GR_Face.h"
#include "GR_GRCurve.h"

class Cell;
class CubitVector;
class CubitBox;
class GeometryEntity;
class Vert;

class BdryEdgeBase : public BFace {

public:

  enum SplitType {
    UNKNOWN, MID_TVT, MID_PARAM
  };

protected:

  //The curve associated with this boundary edge.
  GRCurve* m_curve;

  //The parameters of the beg and end verts on the curve.
  //Storing this here causes redundancy storage, however,
  //it is better than storing the param in the Vert class
  //since interior vertices have nothing to do with a parameter.
  double m_vert0_param, m_vert1_param;

  //Default constructor.
  BdryEdgeBase() :
      BFace(false, pFInvalidFace, pFInvalidFace), m_curve(
      NULL), m_vert0_param(-LARGE_DBL), m_vert1_param(LARGE_DBL), m_bdryPatch(
      NULL)
  {
  }
  ;

  // Typical constructor; called from BdryEdge and IntBdryEdge
  // constructors with first two args guaranteed to exist.
  BdryEdgeBase(const bool qIsInternal, Face* const pFLeft, Face* const pFRight =
  pFInvalidFace,
	       GRCurve* const parent_curve =
	       NULL,
	       double vert0_param = -LARGE_DBL, double vert1_param =
	       LARGE_DBL);

  // Copy construction disallowed, but declared protected so that
  // derived classes can pretend to use it.
  BdryEdgeBase(BdryEdgeBase& BEB); // Never defined
  BdryEdgeBase&
  operator=(const BdryEdgeBase& BEB); // Never defined

public:

  virtual
  ~BdryEdgeBase()
  {
  }

  //Sets the geometric curve associated with this boundary edge.
  void
  setCurve(GRCurve* const curve);

  //Returns a pointer to the geometric curve object associated to the BE.
  GRCurve*
  getCurve() const;
  bool
  doesCurveExist() const;

  //Sets the parameters on curve at begin and end points.
  void
  setVert0Param(double param);
  void
  setVert1Param(double param);

  //Returns parameter on curve at begin and end points.
  double
  getVert0Param();
  double
  getVert1Param();

  //Returns the boundary condition on left or right.
  int
  getBdryCondition(const bool left) const;
  int
  getLeftBdryCondition() const;
  int
  getRightBdryCondition() const;

  //Return the distance between pVVert(0) and pVVert(1)
  double
  calcLength() const;

  //Returns the edge's midpoint
  CubitVector
  calcMidpoint() const;

  //Is pVVert(0) before pVVert(1) along the curve?
  //The mesh is generated in such a way that the faces located at the boundary
  //are always in the same direction as the underlying curve. If the boundary
  //edge is in the same direction as its face, then it is also in the same direction
  //as the underlying curve. MIGHT NOT HOLD AFTER Mesh2D::vReorder() IS CALLED!
  bool
  isForward() const
  {
    return (getVert(0) == getFace(0)->getVert(0));
  }

  //Returns the number of verts. Always returns 2.
  int
  getNumVerts() const
  {
    return 2;
  }

  //This is needed by the refinement code. Returns the split parameter
  //and set split_coord. concentric_shell_split set to true if a shell
  //split is required.
  double
  getSplitData(const SplitType split_type, CubitVector& split_coord,
	       bool& concentric_shell_split);

  //Returns the axis-aligned bounding box for this boundary edge.
  CubitBox
  bounding_box() const;

  virtual bool
  isCoplanarWith(const double pointLoc[]) const;
  //Functions to call to verify if the boundary edge is encroached.
  virtual eEncroachResult
  isEncroached(const enum eEncroachType enc_type) const;

  virtual eEncroachResult
  isPointEncroachingBall(const double coord[],
			 const bool tie_break = false) const;

  virtual eEncroachResult
  isPointEncroachingLens(const double coord[],
			 const bool tie_break = false) const;

  virtual eEncroachResult
  isPointEncroachingNewLens(const double coord[],
			    const bool tie_break = false) const;

  //Offset insertion is disabled for now.
  /*   CubitVector get_offset_loc(const CubitVector& trigger) const; */
  /*   CubitVector get_offset_loc(bool left = true) const; */

  //These two return the center poit of this edge.
  void
  calcCentroid(double centroid[]) const;
  void
  calcCircumcenter(double circumcenter[]) const;

  //Sets and returns the geometric entity associated with this boundary face.
  //Pure virtual functions defined in BFace.
  void
  setGeometry(GeometryEntity* const geom_ent);
  GeometryEntity*
  getGeometry();

private:

  eEncroachResult
  isEncroachingLens(const double coord[], const bool tie_break,
		    const eEncroachType enc_type) const;

public:

  /// Set bdry condition to a given integer.  This is a temporary hack
  /// to deal with reading a mesh from native format.
  void
  setBdryCondition(const int i)
  {
    assert(i >= 0);
    m_bdryCond = i;
  }
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  //The following are no longer necessary but are defined//
  //because base class declare them as pure virtual.     //
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  int
  getBdryCondition() const
  {
    //This is the original function called by the BFace
    //interface. It only allows to store and return a
    //boundary on one side of a curve. Still maintained for
    //compatibility with IO functions generated by lex.c
    if (m_bdryCond >= 0)
      return m_bdryCond;

    int left = getLeftBdryCondition();
    int right = getRightBdryCondition();

#ifndef NDEBUG

    assert(
	getType() == CellSkel::eBdryEdge
	    || getType() == CellSkel::eIntBdryEdge /*|| getType() == CellSkel::ePseudoSurface*/ );

    switch (getType())
      {
      case CellSkel::eBdryEdge:
	assert(left == 0 || right == 0);
	assert(!(left == 0 && right == 0));
	break;
      case CellSkel::eIntBdryEdge:
	assert(left == iInvalidBC && right == iInvalidBC);
	break;
      }

#endif

    return left == 0 ? right : left;

  }

  //The following has been re-instated because it is necessary
  //for the "original" 3D mesh generation code.

private:

  BdryPatch* m_bdryPatch;

public:

  //This should not be necessary.
  int
  getPatchIndex()
  {
    //Patches do not exist in 2D anymore. However, do not want
    //to remove from BFace interface just yet.
    assert(0);
    return 0;
  }

  void
  setPatch(BdryPatch* const pBPBdryPatchIn)
  {
    //Patches do not exist in 2D anymore. However, do not want
    //to remove from BFace interface just yet.
    assert(pBPBdryPatchIn);
    m_bdryPatch = pBPBdryPatchIn;
  }

  BdryPatch*
  getPatchPointer() const
  {
    // assert(0);
    // assert(m_bdryPatch);
    return m_bdryPatch;
  }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  void
  setTopologicalParent(BasicTopologyEntity* const topo_ent)
  {
    if (topo_ent) {
    }
    assert(0);
    abort();
  }
  ;
  BasicTopologyEntity*
  getTopologicalParent()
  {
    return static_cast<BasicTopologyEntity*>(NULL);
  }
  virtual void
  findClosestPointOnBFace(const double pt[3], double close[3]);
};

class BdryEdge : public BdryEdgeBase {

private:

  // Disallow copy constructor
  BdryEdge(BdryEdge& BE) :
      BdryEdgeBase(BE)
  {
    assert(0);
  }

public:

  BdryEdge(Face* const pF = pFInvalidFace, GRCurve* const parent_curve = NULL,
	   double param_at_beg = -LARGE_DBL, double param_at_end = LARGE_DBL) :
      BdryEdgeBase(false, pF, pFInvalidFace, parent_curve, param_at_beg,
		   param_at_end)
  {
    setType(Cell::eBdryEdge);
  }

  BdryEdge&
  operator=(BdryEdge& BE)
  {

    if (this != &BE)
      this->BdryEdgeBase::operator=(BE);

    return *this;

  }
public:
  virtual
  ~BdryEdge()
  {
  }

  int
  getNumFaces() const
  {
    return 1;
  }

};

class IntBdryEdge : public BdryEdgeBase {

private:

  // Disallow copy constructor
  IntBdryEdge(IntBdryEdge& IBE) :
      BdryEdgeBase(IBE)
  {
    assert(0);
  }

public:

  IntBdryEdge(Face* const pFLeft = pFInvalidFace, Face* const pFRight =
  pFInvalidFace,
	      GRCurve* const parent_curve = NULL, double param_at_beg =
		  -LARGE_DBL,
	      double param_at_end = LARGE_DBL) :
      BdryEdgeBase(true, pFLeft, pFRight, parent_curve, param_at_beg,
		   param_at_end)
  {
    setType(Cell::eIntBdryEdge);
  }

  IntBdryEdge&
  operator=(IntBdryEdge& IBE)
  {

    if (this != &IBE)
      this->BdryEdgeBase::operator=(IBE);

    return *this;

  }

  virtual
  ~IntBdryEdge()
  {
  }

  int
  getNumFaces() const
  {
    return 2;
  }

  Face *
  getOtherFace(const Face * const pF) const
  {
    if (pF == getFace(0))
      return const_cast<Face*>(getFace(1));
    else
      return const_cast<Face*>(getFace(0));
  }

};

#endif
