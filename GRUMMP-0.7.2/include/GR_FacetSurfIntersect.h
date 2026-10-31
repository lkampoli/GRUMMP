#ifndef GR_FACETSURFINTERSECT_H
#define GR_FACETSURFINTERSECT_H 1

#include <set>
#include <map>

#include "GR_SurfIntersect.h"

//Pre-declarations:

//CGM classes
class CubitFacet;
class CubitFacetEdge;
class CubitPoint;
class CubitVector;
class RefFace;
template<class X>
  class DLILIst;

//GRUMMP classes
class FacetSurfIntersect;
class Vert;

class FacetSurfTri : public SurfTri {

  CubitFacet* m_facet;

public:

  friend class FacetSurfIntersect;

  FacetSurfTri(RefFace* const surface, CubitFacet* const facet,
	       Vert* const vert1, Vert* const vert2, Vert* const vert3);

  FacetSurfTri(RefFace* const surface, CubitFacet* const facet,
	       Vert* const vert1, double* const curv1, Vert* const vert2,
	       double* const curv2, Vert* const vert3, double* const curv3);

  virtual
  ~FacetSurfTri()
  {
  }

  CubitFacet*
  get_facet() const
  {
    return m_facet;
  }

  //Projects coord located on m_facet to the curved patch of the facet.
  virtual void
  project_to_surface(const CubitVector& coord, CubitVector* const surf_coord,
		     CubitVector* const surf_normal = NULL) const;

private:

  //For the sake of defining them...
  FacetSurfTri() :
      SurfTri(), m_facet(NULL)
  {
    assert(0);
  }
  FacetSurfTri(const FacetSurfTri&) :
      SurfTri(), m_facet(NULL)
  {
    assert(0);
  }
  FacetSurfTri&
  operator=(const FacetSurfTri&)
  {
    assert(0);
    return *this;
  }

};

//Surface intersection tool designed to be used with
//CGM's Facet-based modeler.

class FacetSurfIntersect : public SurfIntersect {

  std::set<Vert*> m_facetsurf_verts;
public:

  //Creates an empty intersect object.
  FacetSurfIntersect();

  //Creates an intersect object for surface
  FacetSurfIntersect(RefFace* const surface);

  //Creates an intersect object for contents of surfaces.
  FacetSurfIntersect(const std::set<RefFace*>& surfaces);

  virtual
  ~FacetSurfIntersect();

private:

  void
  init_intersect(RefFace* const surface);

  void
  init_surf_pts(const std::set<CubitPoint*>& surf_pts,
		const std::set<CubitPoint*>& curv_pts,
		const std::set<CubitPoint*>& apex_pts,
		const std::set<CubitFacet*>& admissible_facets,
		std::map<CubitPoint*, Vert*>& vert_map,
		std::map<CubitPoint*, double*>& curv_map);

  void
  init_curv_pts(const std::set<CubitPoint*>& surf_pts,
		const std::set<CubitPoint*>& curv_pts,
		const std::set<CubitPoint*>& apex_pts,
		const std::set<CubitFacetEdge*>& admissible_edges,
		std::map<CubitPoint*, Vert*>& vert_map,
		std::map<CubitPoint*, double*>& curv_map);

  void
  init_apex_pts(const std::set<CubitPoint*>& surf_pts,
		const std::set<CubitPoint*>& curv_pts,
		const std::set<CubitPoint*>& apex_pts,
		const std::set<CubitPoint*>& admissible_pts,
		std::map<CubitPoint*, Vert*>& vert_map,
		std::map<CubitPoint*, double*>& curv_map);

  void
  init_point_handles(const std::map<CubitPoint*, Vert*>& vert_map);

  void
  init_curv_handles(const std::map<CubitPoint*, double*>& curv_map);

  void
  init_surf_tris(RefFace* const surface, DLIList<CubitFacet*>& facets,
		 const std::map<CubitPoint*, Vert*>& vert_map,
		 const std::map<CubitPoint*, double*>& curv_map);

  Vert*
  new_vert_at_point(CubitPoint* const my_point);

  //Returns points based 
  void
  get_points(DLIList<CubitPoint*>& the_points,
	     std::set<CubitPoint*>& apex_pts) const;

  void
  get_points(DLIList<CubitFacetEdge*>& the_edges,
	     std::set<CubitPoint*>& curve_pts) const;

  void
  get_points(DLIList<CubitFacet*>& the_facets,
	     std::set<CubitPoint*>& surface_pts) const;

  bool
  points_valid(const std::set<CubitPoint*>& apex_pts,
	       const std::set<CubitPoint*>& curve_pts,
	       const std::set<CubitPoint*>& surface_pts) const;

};

#endif
