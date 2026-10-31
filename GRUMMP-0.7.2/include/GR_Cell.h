#ifndef GR_Cell
#define GR_Cell 1

#include "GR_config.h"
#include "GR_CellSkel.h"
#include "GR_Geometry.h"
class CubitBox;
class CubitVector;

class Cell : public CellSkel {
private:
  int indexx = 1024;
protected:
  /// Copy construction disallowed.
  Cell(const Cell&);
  ///
  Face** m_faces;

protected:
  /// Constructor for all types of bdry faces; all args are provided.
  Cell(const bool isInternalBdry, Face* const leftFace, Face* const rightFace);
  ///
  Cell(const int numFaces, const int numVerts, Face** const facesIn = 0);
public:
  ///
  virtual
  ~Cell();
  ///
  Cell&
  operator=(const Cell& C);
  ///
  /// Ensure a sane initial state for cell data
  void
  resetAllData();
  /// Asks if it is a bdry cell.. return false by default
  virtual bool
  isBdryCell() const
  {
    return false;
  }
  ///
  virtual void
  setIndex();
  virtual int
  getIndex();
  ///
  virtual int
  doFullCheck() const;
  ///
  virtual int
  getNumFaces() const = 0;
  ///
  virtual int
  getNumVerts() const = 0;
  ///
  int
  getNumCells() const
  {
    return 0;
  }

  ///
  virtual const Vert*
  getVert(const int i) const = 0;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }

  ///
  const Face*
  getFace(const int i) const;

  Face*
  getFace(const int i);
  void
  getAllFaceHandles(GRUMMP_Entity* aHandles[]) const;

  ///
  const Cell*
  getCell(const int /*i*/) const
  {
    return this;
  }
  Cell*
  getCell(const int /*i*/)
  {
    return this;
  }
  void
  getAllCellHandles(GRUMMP_Entity*[]) const
  {/* Do nothing */
  }

  ///
  bool
  hasFace(const Face* pF) const;
  ///
  virtual bool
  hasVert(const Vert* pV) const;
  ///
  virtual const Face*
  getOppositeFace(const Vert* const pVVert) const;
  virtual Face*
  getOppositeFace(const Vert* const pVVert);
  ///
  virtual Vert*
  getOppositeVert(const Face* const pFFace) const;
  ///
  virtual void
  calcVecSize(double[]) const
  {
    assert(0);
  }
  ///
  virtual double
  calcSize() const = 0;
  /// The boolean here is true for results in degrees, false for radians.
  virtual void
  calcAllDihed(double[], int* const piN, bool = true) const
  {
    *piN = 0;
  } // Fail silently; pyrs, prisms, and hexes currently get skipped
  virtual void
  calcAllSolid(double[], int* const piN, bool = true) const
  {
    *piN = 0;
  } // Fail silently; pyrs, prisms, and hexes currently get skipped
  virtual bool
  isClosed() const = 0;
  ///
  virtual double
  calcShortestEdgeLength() const
  {
    assert(0);
    return -1.e300;
  }
  ///
  virtual double
  calcCircumradius() const
  {
    assert(0);
    return -1.e300;
  }

  ///
  virtual void
  addFace(Face* const pFNew, const int iFaceIndex = -1);

  ///
  virtual double
  calcInradius() const
  {
    assert(0);
    return -1.e300;
  }
  ///
  void
  replaceFace(const Face* const pFOld, Face* const pFNew);
  ///
  void
  removeFace(const Face* const pFOld);
  ///
  void
  assignFaces(Face* apF[]);
  // Must allow as few as one valid face, because bdry faces have only
  // one face attached to them, and bdry faces think they're cells.
  void
  assignFaces(Face* pF0, Face* pF1 = pFInvalidFace, Face* pF2 =
  pFInvalidFace,
	      Face* pF3 =
	      pFInvalidFace,
	      Face* pF4 = pFInvalidFace, Face* pF5 =
	      pFInvalidFace);
  /// This next one should probably be private, since it's only going to
  /// be called by the cell assignment function...
  virtual void
  canonicalizeFaceOrder()
  {
  }
  ///
  virtual void
  calcCentroid(double adCent[]) const;
  /// The following function name can't be changed because it's required for
  /// a CGM container class that tets are put into somewhere.
  /// TODO: Another CGM dependency to break.
  virtual CubitBox
  bounding_box() const;

  virtual bool
  isCircumcircleEmpty() const;

private:
  double m_anisoCircRad;
  bool bCombineToQuad;

public:
  /// These functions are used only in the aniso refinement queue.  They're
  /// actually dealing with aniso circumradius.
  /// TODO: Migrate these to the TriCell class.
  /// TODO: Generalize circumradius calculation to apply in metric or regular
  /// spaces as needed.
  double
  getAnisoCircumradius() const
  {
//    assert(m_anisoCircRad>0);
    return m_anisoCircRad;
  }
  void
  setAnisoCircumradius(double dCircRadIn)
  {
    m_anisoCircRad = dCircRadIn;
  }

  void
  CombineToQuad(bool Quad)
  {
    bCombineToQuad = Quad;
  }

  bool
  isCombineToQuad()
  {
    return bCombineToQuad;
  }

};

///
class SimplexCell : public Cell {
protected:
  /// Copy construction disallowed.
  SimplexCell(const SimplexCell& SC) :
      Cell(SC)
  {
    assert(0);
  }
  ///
  SimplexCell(const int iNFaces, const int iNVerts, Face** const ppFIn) :
      Cell(iNFaces, iNVerts, ppFIn)
  {
  }
public:
  ///
  const Face*
  getOppositeFace(const Vert* const pVVert) const;
  Face*
  getOppositeFace(const Vert* const pVVert);
  ///
  Vert*
  getOppositeVert(const Face* const pFFace) const;
  ///
  virtual void
  calcCircumcenter(double adRes[]) const = 0;
  ///
  virtual double
  calcCircumradius() const = 0;
  ///
  virtual bool
  circumscribesPoint(const double adPoint[]) const = 0;
  ///
  virtual void
  calcContainmentCenter(double adRes[]) const = 0;
  ///
  virtual double
  calcContainmentRadius() const = 0;
  ///
  virtual void
  calcBarycentricCoords(const double adPoint[], double adBary[]) const = 0;
  virtual bool
  isSmallAngleAffected() const = 0;

  virtual bool
  isPointOnEdge(const double adPoint[], GR_index_t & iV0, GR_index_t & iV1) = 0;

};

///
class TetCell : public SimplexCell {
  /// Copy construction disallowed.
  TetCell(const TetCell& TC) :
      SimplexCell(TC)
  {
    assert(0);
  }
public:
  ///
  TetCell(Face** const ppFIn = 0) :
      SimplexCell(4, 4, ppFIn)
  {
    setType(Cell::eTet);
  }
  /// Returns number of faces
  int
  getNumFaces() const
  {
    return 4;
  }
  /// Returns number of vertices
  int
  getNumVerts() const
  {
    return 4;
  }

  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandles[]) const;

  ///
  bool
  isClosed() const;
  ///
  virtual double
  calcSize() const;
  ///
  void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const;
  ///
  void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const;
/// Declared in SimplexCell
  void
  calcCircumcenter(double adRes[]) const;
  ///
  double
  calcCircumradius() const;
  ///
  double
  calcInradius() const;
  ///
  bool
  circumscribesPoint(const double adPoint[]) const;
  ///
  void
  calcContainmentCenter(double adRes[]) const;
  ///
  double
  calcContainmentRadius() const;
  ///
  void
  calcBarycentricCoords(const double adPoint[], double adBary[]) const;
  ///
  virtual bool
  isSmallAngleAffected() const
  {
    assert(0);
    return false;
  }
  ///
  bool
  isPointOnEdge(const double adPoint[], GR_index_t & iV0, GR_index_t & iV1);
  ///
  double
  calcShortestEdgeLength() const;

  //Checks if vertex is located inside tetra's circumsphere using
  //adaptive precision predicates.
  bool
  circumscribesPointAdaptive(const Vert* const vertex);
  int
  getEntType() const
  {
    return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    return iMesh_TETRAHEDRON;
  }
};

///
class TriCell : public SimplexCell {
  /// Copy construction disallowed.
  TriCell(const TriCell& TC) :
      SimplexCell(TC)
  {
    assert(0);
  }
//private:
//  int indexx = 1024;
public:
  ///
  TriCell(Face** const ppFIn = 0) :
      SimplexCell(3, 3, ppFIn)
  {
    setType(Cell::eTriCell);
  }
/// Declared in Cell

  ///
  int
  getNumFaces() const
  {
    return 3;
  }
  ///
  int
  getNumVerts() const
  {
    return 3;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandles[]) const;

  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  ///
  void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const;
  ///
  void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const;
/// Declared in SimplexCell
  void
  calcCircumcenter(double adRes[]) const;
  ///
  double
  calcCircumradius() const;
  ///
  bool
  circumscribesPoint(const double adPoint[]) const;
  ///
  void
  calcContainmentCenter(double adRes[]) const;
  ///
  double
  calcContainmentRadius() const;
  ///
  void
  calcBarycentricCoords(const double adPoint[], double adBary[]) const;
  ///
  virtual bool
  isSmallAngleAffected() const;

  double
  calcShortestEdgeLength() const;
  int
  getLargestLayernumber();
  int
  getSmallestLayernumber();
  bool
  isPointOnEdge(const double adPoint[], GR_index_t & iV0, GR_index_t & iV1);

  /// Things peculiar to surface cells
  void
  calcVecSize(double adRes[]) const;
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
};

///
class QuadCell : public Cell {
  /// Copy construction disallowed.
  QuadCell(const QuadCell& QC) :
      Cell(QC)
  {
    assert(0);
  }
public:
  ///
  QuadCell(Face** const ppFIn = 0) :
      Cell(4, 4, ppFIn)
  {
    setType(Cell::eQuadCell);
  }
/// Declared in Cell
  int
  getNumFaces() const
  {
    return 4;
  }
  ///
  int
  getNumVerts() const
  {
    return 4;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* /*aHandles*/[]) const;
  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  ///
  void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const;
  ///
  void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const;
  /// Things peculiar to surface cells
  void
  calcVecSize(double adRes[]) const;
  // iMesh support calls.
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
  // Get the faces into cyclic order!
  virtual void
  canonicalizeFaceOrder();
};

///
class FlakeCell : public Cell {
  /// Copy construction disallowed.
  FlakeCell(const FlakeCell&);
public:
  ///
  FlakeCell(Face** const ppFIn = 0) :
      Cell(3, 4, ppFIn)
  {
    setType(Cell::eFlake);
  }
/// Declared in Cell
  int
  getNumFaces() const
  {
    return 3;
  }
  ///
  int
  getNumVerts() const
  {
    return 4;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* /*aHandles*/[]) const;
  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  int
  getEntType() const
  {
    return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    return iMesh_FLAKE;
  }
};

///
class PyrCell : public Cell {
  /// Copy construction disallowed.
  PyrCell(const PyrCell& PCIn) :
      Cell(PCIn)
  {
    assert(0);
  }
public:
  ///
  PyrCell(Face** const ppFIn = 0) :
      Cell(5, 5, ppFIn)
  {
    setType(Cell::ePyr);
  }
/// Declared in Cell
  int
  getNumFaces() const
  {
    return 5;
  }
  ///
  int
  getNumVerts() const
  {
    return 5;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* /*aHandles*/[]) const;
  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  int
  getEntType() const
  {
    return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    return iMesh_PYRAMID;
  }
};

///
class PrismCell : public Cell {
  /// Copy construction disallowed.
  PrismCell(const PrismCell& PCIn) :
      Cell(PCIn)
  {
    assert(0);
  }
public:
  ///
  PrismCell(Face** const ppFIn = 0) :
      Cell(5, 6, ppFIn)
  {
    setType(Cell::ePrism);
  }
/// Declared in Cell
  int
  getNumFaces() const
  {
    return 5;
  }
  ///
  int
  getNumVerts() const
  {
    return 6;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* /*aHandles*/[]) const;
  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  int
  getEntType() const
  {
    return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    return iMesh_PRISM;
  }
};

///
class HexCell : public Cell {
  /// Copy construction disallowed.
  HexCell(const HexCell& HC) :
      Cell(HC)
  {
    assert(0);
  }
public:
  ///
  HexCell(Face** const ppFIn = 0) :
      Cell(6, 8, ppFIn)
  {
    setType(Cell::eHex);
  }
/// Declared in Cell
  int
  getNumFaces() const
  {
    return 6;
  }
  ///
  int
  getNumVerts() const
  {
    return 8;
  }
  ///
  virtual const Vert*
  getVert(const int i) const;
  virtual Vert*
  getVert(const int i)
  {
    return this->Entity::getVert(i);
  }
  void
  getAllVertHandles(GRUMMP_Entity* /*aHandles*/[]) const;
  ///
  double
  calcSize() const;
  ///
  bool
  isClosed() const;
  int
  getEntType() const
  {
    return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    return iMesh_HEXAHEDRON;
  }
};

///
Face*
findCommonFace(Cell* const pC0, Cell* const pC1);

#include "GR_AdaptPred.h"

inline bool
TetCell::circumscribesPointAdaptive(const Vert* const vertex)
{

  assert(vertex);

  //Must have the correct orientation for the insphere predicate to work.
//  assert(
//      orient3d_shew(getVert(0)->getCoords(), getVert(1)->getCoords(), getVert(2)->getCoords(), getVert(3)->getCoords()) < 0.);
//  return (insphere_shew(getVert(0)->getCoords(), getVert(1)->getCoords(),
//                        getVert(2)->getCoords(), getVert(3)->getCoords(),
//                        vertex->getCoords()) < 0.);
  assert(
      checkOrient3D(getVert(0)->getCoords(), getVert(1)->getCoords(),
		    getVert(2)->getCoords(), getVert(3)->getCoords()) > 0);

  return (isInsphere(getVert(0)->getCoords(), getVert(1)->getCoords(),
		     getVert(2)->getCoords(), getVert(3)->getCoords(),
		     vertex->getCoords()) > 0);

}

void
evalAllDihedsFromCoords(int& nDiheds, double *result, const double ptA[],
			const double ptB[], const double ptC[],
			const double ptD[], const double ptE[],
			const double ptF[]);

#endif
