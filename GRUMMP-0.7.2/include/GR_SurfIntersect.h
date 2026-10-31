#ifndef GR_SURFINTERSECT_H
#define GR_SURFINTERSECT_H

#include "GR_config.h"
#include "KDDTree.hpp"
#include "CubitVector.hpp"
#include <vector>
#include <list>
#include <set>
#include <map>

class CubitBox;
class RefFace;
class SurfIntersect;
class Vert;
class VoronoiEdge;

class SurfTri {

  struct PointData {
    PointData() :
	vertex(NULL), curvature(NULL)
    {
    }
    Vert* vertex;
    double* curvature;
  };

  //The RefFace on which this triangle lies
  RefFace* m_surface;

  //Pointer to the triangle's points and curvature
  PointData m_point_data[3];

  //The AreaOrder functor requires to know the area of the triangle 
  //Since it is used rarely, this variable defaults to NULL to save memory.
  //However, before calling anything using the functor, make sure it is
  //instantiated.
  //(could probably use a size measure slightly cheaper to compute, but hey
  //it is only done once and probably only when meshing smooth surface...)
  double* m_area;

  SurfTri(const SurfTri&);
  SurfTri&
  operator=(const SurfTri&);

  //Computes the mapping of coord in the triangle's parametric space
  void
  param_at_coord(const CubitVector& coord, double& s, double& t) const;

protected:

  SurfTri();

  //constructor without specified curvature
  SurfTri(RefFace* const surface, Vert* const vert1, Vert* const vert2,
	  Vert* const vert3);

  //constructor with specified curvature
  SurfTri(RefFace* const surface, Vert* const vert1, double* const curv1,
	  Vert* const vert2, double* const curv2, Vert* const vert3,
	  double* const curv3);

public:

  friend class SurfIntersect;

  virtual
  ~SurfTri();

  void
  set_curvature_at_vert(Vert* const vertex, double* const curvature);

  //Computes the area based on the vertices' coords.
  //Returns the area. If area already set, then returns this area without recomputing.
  double
  set_area();

  //Accessor functions.
  RefFace*
  get_surface() const
  {
    return m_surface;
  }
  Vert*
  get_vertex(int i) const;
  double*
  get_curvature(int i) const;
  double*
  get_area() const
  {
    return m_area;
  }

  //Axis-aligned bounding box
  CubitBox
  bounding_box() const;

  //
  double
  get_curvature_val(int i) const;
  CubitVector
  get_centroid() const;

  //Linear interpolation of the curvature over the triangle.
  double
  curvature_at_coord(const CubitVector& coord) const;

  //Projects coord to curved surface. In this general case,
  //only calls the corresponding function in the RefFace class.
  //Derived implementation makes better use of available data.
  virtual void
  project_to_surface(const CubitVector& coord, CubitVector* const surf_coord,
		     CubitVector* const surf_normal = NULL) const;

  //Used to order the triangles in terms of their area.
  struct AreaOrder {
    bool
    operator()(SurfTri* const tri1, SurfTri* const tri2) const
    {
      if (tri1 == tri2)
	return false;
      return tri1->set_area() >= tri2->set_area();
    }
  };

};

class SurfIntersect {

  typedef std::map<Vert*, std::list<SurfTri*> > VertFaceMap;

  //Inits the surface of m_surf_tree given facets.
  virtual void
  init_intersect(RefFace* const surface);

  //Computes and inits discrete approximation of curvature at vertices.
  virtual void
  init_curvature(const VertFaceMap& vert_to_faces);

  SurfIntersect(const SurfIntersect&);
  SurfIntersect&
  operator=(const SurfIntersect&);

protected:

  //Flag to know if curvature is stored.
  bool m_has_curvature_data;

  //To store SurfTri's and curvatures.

  std::set<Vert*> m_all_points;
  std::vector<double*> m_all_curvatures;

  //Spatial indexing tree for intersection search
  KDDTree<SurfTri*> m_surf_tree;
public:

  struct IntersectCompare;
  typedef std::set<SurfTri*, SurfTri::AreaOrder> SurfTriSet;

  struct Intersection {

    SurfTri* m_triangle;
    CubitVector m_intersect_coord;

    Intersection() :
	m_triangle(NULL), m_intersect_coord()
    {
    }
    Intersection(const Intersection& intersect) :
	m_triangle(intersect.m_triangle), m_intersect_coord(
	    intersect.m_intersect_coord)
    {
    }
    Intersection&
    operator=(const Intersection& intersect)
    {
      if (&intersect != this) {
	m_triangle = intersect.m_triangle;
	m_intersect_coord.set(intersect.m_intersect_coord);
      }
      return *this;
    }
    virtual
    ~Intersection()
    {
    }

    void
    reset_data()
    {
      m_triangle = static_cast<SurfTri*>(NULL);
      m_intersect_coord.set(0., 0., 0.);
    }

    void
    set_triangle(SurfTri* const triangle)
    {
      m_triangle = triangle;
    }
    const SurfTri*
    get_triangle() const
    {
      return m_triangle;
    }

    void
    set_intersect_coord(const CubitVector& coord)
    {
      m_intersect_coord.set(coord);
    }
    const CubitVector&
    get_intersect_coord() const
    {
      return m_intersect_coord;
    }

  };

  typedef std::set<Intersection, IntersectCompare> IntersectSet;

  //Constructor / destructor, setting the with_discrete_curvature switch to true
  //will cause the function to estimate the surface curvature at every vertex using
  //discrete differential geometry operators (only the maximum principal curvature is
  //stored - extending to include other curvature data should not be hard for the
  //motivated coder...).

  //Creates an empty SurfIntersect object.
  SurfIntersect(bool with_curvature = false);

  //Creates a SurfIntersect object for surface.
  SurfIntersect(RefFace* const surface, bool with_curvature = false);

  //Creates a SurfIntersect object for members of surfaces. 
  SurfIntersect(const std::set<RefFace*>& surfaces,
		bool with_curvature = false);

  virtual
  ~SurfIntersect();

  bool
  contains_curvature_data() const
  {
    return m_has_curvature_data;
  }

  void
  get_triangles(SurfTriSet& triangles);

  //Finds the points where a line going from beg_pt to end_pt
  //intersects the DISCRETE approximation of the surface.
  void
  intersect_with_line(const CubitVector& beg_pt, const CubitVector& end_pt,
		      IntersectSet& intersect);

  //Intersects DISCRETE surface with a Voronoi edge.
  void
  intersect_with_voronoi_edge(const VoronoiEdge& vor_edge,
			      IntersectSet& intersect);

  //A functor to compare intersection coordinates. operator() returns true when
  //iFuzzyComp returns 0 for the three components of the vectors, false otherwise.
  struct IntersectCompare {
    bool
    operator()(const SurfIntersect::Intersection& intersect1,
	       const SurfIntersect::Intersection& intersect2) const;
  };

  //This is deprecated: A functor that determines the farthest location from m_point.
  struct FarthestPoint {
    double m_max_dist;
    CubitVector *m_point, *m_farthest;
    double *m_curvature;

    FarthestPoint(CubitVector* const point, CubitVector* const farthest,
		  double* const curvature);
    void
    operator()(const Intersection& intersect);
  };

  //Replaced by this:
  struct FarthestIntersect {
    double m_max_dist;
    CubitVector m_from_point;
    Intersection* m_intersect;
    FarthestIntersect(const CubitVector& from_coord,
		      Intersection* const intersect);
    void
    operator()(const Intersection& intersect);
  };

  //Given two CubitVectors, finds the minimal and maximal components
  static CubitVector
  min_vec(const CubitVector& vec1, const CubitVector& vec2);
  static CubitVector
  max_vec(const CubitVector& vec1, const CubitVector& vec2);
  static CubitVector
  min_vec(const CubitVector* const vec1, const CubitVector* const vec2);
  static CubitVector
  max_vec(const CubitVector* const vec1, const CubitVector* const vec2);

  //Finds intersection between a line (passing by l0 and l1) and
  //the triangle (p0, p1, p2).
  static bool
  intersect_line_triangle(const CubitVector& l0, const CubitVector& l1,
			  const CubitVector& p0, const CubitVector& p1,
			  const CubitVector& p2, const bool restrict_to_ray,
			  const bool restrict_to_segment,
			  CubitVector& intersection_coord);

  //Returns true if a LINE defined by l0 and l1 intersects a triangle 
  //defined by p0, p1 and p2 (the intersection coord is also returned).
  static bool
  intersect_line_triangle(const CubitVector& l0, const CubitVector& l1,
			  const CubitVector& p0, const CubitVector& p1,
			  const CubitVector& p2,
			  CubitVector& intersection_coord);

  //Returns true if a LINE defined by l0 and l1 intersects a triangle 
  //defined by p0, p1 and p2 (the intersection coord is also returned).
  static bool
  intersect_ray_triangle(const CubitVector& l0, const CubitVector& l1,
			 const CubitVector& p0, const CubitVector& p1,
			 const CubitVector& p2,
			 CubitVector& intersection_coord);

  //Returns true if a LINE defined by l0 and l1 intersects a triangle 
  //defined by p0, p1 and p2 (the intersection coord is also returned).
  static bool
  intersect_segment_triangle(const CubitVector& l0, const CubitVector& l1,
			     const CubitVector& p0, const CubitVector& p1,
			     const CubitVector& p2,
			     CubitVector& intersection_coord);

};

inline CubitVector
SurfIntersect::min_vec(const CubitVector& vec1, const CubitVector& vec2)
{

  return CubitVector(std::min(vec1.x(), vec2.x()), std::min(vec1.y(), vec2.y()),
		     std::min(vec1.z(), vec2.z()));

}

inline CubitVector
SurfIntersect::min_vec(const CubitVector* const vec1,
		       const CubitVector* const vec2)
{

  return min_vec(*vec1, *vec2);

}

inline CubitVector
SurfIntersect::max_vec(const CubitVector& vec1, const CubitVector& vec2)
{

  return CubitVector(std::max(vec1.x(), vec2.x()), std::max(vec1.y(), vec2.y()),
		     std::max(vec1.z(), vec2.z()));

}

inline CubitVector
SurfIntersect::max_vec(const CubitVector* const vec1,
		       const CubitVector* const vec2)
{

  return CubitVector(*vec1, *vec2);

}

#endif
