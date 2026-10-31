#ifndef GR_BARYCENTRICLENGTHSCALE_H_
#define GR_BARYCENTRICLENGTHSCALE_H_

#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"

/// Really simple object, carries the mesh and a guess of the cell

/// Functions currently in SurfaceInsertion2D/3D.cxx

/// Basic idea is that for any coordinate in the mesh, it gives the length scale
/// using barycentric interpolation on the length scales of the vertices of the
/// cell. Uses the previous cell to traverse around.

/// Stores a guess for where it found its last point, cell wise
/// just an easier way when working through the mesh to avoid having to go aimlessly
/// assumes points are close by when searching

/// Might make use of search tree in the future, but right now
/// seems to move around the mesh well enough,

/// In general, one global search upon creation for the first point,
/// Then traversal

/// 2D has anisotropic ability, 3D does not yet.
/// 2D should really be barycentric and length scale 2D,
class BaryLengthScale2D {

protected:
  const Mesh2D * m_mesh;

  // the last cell we were in
  Cell * m_cellguess;

  // if we are in a curve, used for anisotropic only
  GRCurve * m_curve;

  // true if anisotropic
  bool m_aniso;

public:

  BaryLengthScale2D(const Mesh2D *mesh, bool aniso = false) :
      m_mesh(mesh), m_cellguess(m_mesh->getCell(0)), m_curve(NULL), m_aniso(
	  aniso)
  {
  }
  ;
  ~BaryLengthScale2D()
  {
  }
  ;

  Cell *
  getCell()
  {
    return m_cellguess;
  }

  double
  getLengthScale(const double & param, Cell * guess = NULL);
  double
  getAnisoLengthScale(const double & param, Cell * guess = NULL,
		      CubitVector * dir = NULL);

  void
  setCurve(GRCurve * curve)
  {
    assert(curve);
    m_curve = curve;
  }
  ;
  GRCurve *
  getCurve()
  {
    return m_curve;
  }

  void
  getMetric(const CubitVector & point, double metric[], Cell * guess =
  NULL);
  // taking CubitVectors in here takes some of the burden off of surface insertion
  double
  getAnisoLength(const CubitVector & pointA, const CubitVector & pointB,
		 Cell * guessA, Cell * guessB);

  // a modification on the mesh version
  Cell *
  findCell(const CubitVector & point, Cell * guess);
  Cell *
  findCell(const double * coords, Cell * guess);

};

class BaryLengthScale3D {

  const VolMesh * m_mesh;
  std::map<const Vert* const, double> m_lengthLookup;
  Cell * m_cellguess;
public:
  BaryLengthScale3D(const VolMesh *mesh) :
      m_mesh(mesh), m_cellguess(m_mesh->getCell(0))
  {
  }
  ;
  ~BaryLengthScale3D()
  {
  }
  ;
private:
  double
  getLengthScale(const double * coords);
public:
  double
  getLengthScale(const Vert * const vert);
  double
  getLengthScale(const CubitVector & cv)
  {
    double coords[3];
    cv.get_xyz(coords);
    return getLengthScale(coords);
  }
  Cell *
  getCell()
  {
    return m_cellguess;
  }
};

#endif /* GR_BARYCENTRICLENGTHSCALE_H_ */
