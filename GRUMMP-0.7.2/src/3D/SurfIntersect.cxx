#include "GR_SurfIntersect.h"
#include "GR_misc.h"
#include "GR_DiscreteCurvTool.h"
#include "GR_Geometry.h"
#include "GR_SurfaceSampler.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "GeometryQueryTool.hpp"
#include "GMem.hpp"
#include "RefFace.hpp"

#include <algorithm>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <utility>

using std::deque;
using std::list;
using std::make_pair;
using std::map;
using std::vector;
using std::set;

SurfTri::SurfTri() :
    m_surface(NULL), m_area(NULL)
{
}

SurfTri::SurfTri(const SurfTri&) :
    m_surface(NULL), m_area(NULL)
{
  assert(0);
}

SurfTri&
SurfTri::operator=(const SurfTri&)
{
  assert(0);
  return *this;
}

SurfTri::SurfTri(RefFace* const surface, Vert* const vert1, Vert* const vert2,
		 Vert* const vert3) :
    m_surface(surface), m_area(NULL)
{

  m_point_data[0].vertex = vert1;
  m_point_data[1].vertex = vert2;
  m_point_data[2].vertex = vert3;

}

SurfTri::SurfTri(RefFace* const surface, Vert* const vert1, double* const curv1,
		 Vert* const vert2, double* const curv2, Vert* const vert3,
		 double* const curv3) :
    m_surface(surface), m_area(NULL)
{

  m_point_data[0].vertex = vert1;
  m_point_data[1].vertex = vert2;
  m_point_data[2].vertex = vert3;

  m_point_data[0].curvature = curv1;
  m_point_data[1].curvature = curv2;
  m_point_data[2].curvature = curv3;

}

void
SurfTri::project_to_surface(const CubitVector& coord,
			    CubitVector* const surf_coord,
			    CubitVector* const surf_normal) const
{

  assert(surf_coord);
  m_surface->find_closest_point_trimmed(coord, *surf_coord);

  if (surf_normal)
    surf_normal->set(m_surface->normal_at(*surf_coord));

}

SurfTri::~SurfTri()
{
  if (m_area)
    delete m_area;
}

void
SurfTri::set_curvature_at_vert(Vert* const vertex, double* const curvature)
{

  assert(curvature);

  if (m_point_data[0].vertex == vertex) {
    m_point_data[0].curvature = curvature;
    return;
  }
  if (m_point_data[1].vertex == vertex) {
    m_point_data[1].curvature = curvature;
    return;
  }
  if (m_point_data[2].vertex == vertex) {
    m_point_data[2].curvature = curvature;
    return;
  }

  vFatalError("Vertex is not found in this triangle\n",
	      "SurfTri::set_curvature_at_vert()");

}

double
SurfTri::set_area()
{

  assert(
      m_point_data[0].vertex && m_point_data[1].vertex
	  && m_point_data[2].vertex);

  if (m_area)
    return *m_area;

  double vec1[] = adDIFF3D(m_point_data[0].vertex->getCoords(),
      m_point_data[1].vertex->getCoords()), vec2[] =
  adDIFF3D(m_point_data[0].vertex->getCoords(),
      m_point_data[2].vertex->getCoords()), vec3[3];

  vCROSS3D(vec1, vec2, vec3);
  m_area = new double(0.5 * dMAG3D(vec3));

  return *m_area;

}

Vert*
SurfTri::get_vertex(int i) const
{

  assert(i >= 0 && i <= 2);
  return m_point_data[i].vertex;

}

double*
SurfTri::get_curvature(int i) const
{

  assert(i >= 0 && i <= 2);
  return m_point_data[i].curvature;

}

double
SurfTri::get_curvature_val(int i) const
{

  assert(i >= 0 && i <= 2);
  assert(m_point_data[i].curvature);
  return *m_point_data[i].curvature;

}

CubitVector
SurfTri::get_centroid() const
{

  assert(
      m_point_data[0].vertex && m_point_data[1].vertex
	  && m_point_data[2].vertex);

  return ((CubitVector(m_point_data[0].vertex->getCoords())
      + CubitVector(m_point_data[1].vertex->getCoords())
      + CubitVector(m_point_data[2].vertex->getCoords())) / 3.);

}

CubitBox
SurfTri::bounding_box() const
{

  Vert* vert[] =
    { m_point_data[0].vertex, m_point_data[1].vertex, m_point_data[2].vertex };

  CubitVector min_vec(MIN3(vert[0]->x(), vert[1]->x(), vert[2]->x()),
		      MIN3(vert[0]->y(), vert[1]->y(), vert[2]->y()),
		      MIN3(vert[0]->z(), vert[1]->z(), vert[2]->z()));
  CubitVector max_vec(MAX3(vert[0]->x(), vert[1]->x(), vert[2]->x()),
		      MAX3(vert[0]->y(), vert[1]->y(), vert[2]->y()),
		      MAX3(vert[0]->z(), vert[1]->z(), vert[2]->z()));

  return CubitBox(min_vec, max_vec);

}

//Computes the s-t coordinates in the triangle's parameteric space
//(triangle formed by points p0, p1 and p2) at coord. Note that coord
//must be coplanar with p0, p1 and p2 to get significant values.
//If coord is located inside the triangle, then 0 <= t <= 1,
//0 <= t  and t + s <= 1.

void
SurfTri::param_at_coord(const CubitVector& coord, double& s, double& t) const
{

  Vert* vert0 = get_vertex(0);
  Vert* vert1 = get_vertex(1);
  Vert* vert2 = get_vertex(2);
//vert0->printVertInfo();
//vert1->printVertInfo();
//vert2->printVertInfo();
//printf("coord %f %f %f\n",coord.x(),coord.y(),coord.z());
//#ifndef NDEBUG
//  double this_coord[] = { coord.x(), coord.y(), coord.z() };
//  assert( checkOrient3D(vert0->getCoords(), vert1->getCoords(),
//		    vert2->getCoords(), this_coord) == 0 );
//#endif

  double u[] =
	{ vert1->x() - vert0->x(), vert1->y() - vert0->y(), vert1->z()
	    - vert0->z() },

  v[] =
	{ vert2->x() - vert0->x(), vert2->y() - vert0->y(), vert2->z()
	    - vert0->z() },

  w[] =
    { coord.x() - vert0->x(), coord.y() - vert0->y(), coord.z() - vert0->z() };

  double uu = dDOT3D(u, u), uv = dDOT3D(u, v), vv = dDOT3D(v, v), wu = dDOT3D(
      w, u), wv = dDOT3D(w, v), den = uv * uv - uu * vv;

  s = (uv * wv - vv * wu) / den;
  t = (uv * wu - uu * wv) / den;
}

double
SurfTri::curvature_at_coord(const CubitVector& coord) const
{

  double s, t;
  param_at_coord(coord, s, t);

  //Coord must be located inside triangle
  assert(iFuzzyComp(s, 0.) != -1 && iFuzzyComp(s, 1.) != 1);
  assert(iFuzzyComp(t, 0.) != -1 && iFuzzyComp(s + t, 1.) != 1);

  return (get_curvature_val(0)
      + s * (get_curvature_val(1) - get_curvature_val(0))
      + t * (get_curvature_val(2) - get_curvature_val(0)));

}

SurfIntersect::SurfIntersect(bool with_curvature) :
    m_has_curvature_data(with_curvature), m_all_points(), m_all_curvatures(), m_surf_tree()
{
}

SurfIntersect::SurfIntersect(RefFace* const surface, bool with_curvature) :
    m_has_curvature_data(with_curvature), m_all_points(), m_all_curvatures(), m_surf_tree()
{

  init_intersect(surface);

}

SurfIntersect::SurfIntersect(const set<RefFace*>& surfaces, bool with_curvature) :
    m_has_curvature_data(with_curvature), m_all_points(), m_all_curvatures(), m_surf_tree()
{

  std::for_each(
      surfaces.begin(), surfaces.end(),
      std::bind1st(std::mem_fun(&SurfIntersect::init_intersect), this));

}

SurfIntersect::SurfIntersect(const SurfIntersect& SI) :
    m_has_curvature_data(SI.m_has_curvature_data), m_all_points(
	SI.m_all_points), m_all_curvatures(SI.m_all_curvatures), m_surf_tree(
	SI.m_surf_tree)
{
  assert(0);
}

SurfIntersect&
SurfIntersect::operator=(const SurfIntersect& SI)
{

  m_has_curvature_data = SI.m_has_curvature_data;
  m_all_points = SI.m_all_points;
  m_all_curvatures = SI.m_all_curvatures;
  m_surf_tree = SI.m_surf_tree;

  assert(0);
  return *this;

}

SurfIntersect::~SurfIntersect()
{

  DLIList<SurfTri*> tri_list;
  CubitVector min(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
  CubitVector max( LARGE_DBL, LARGE_DBL, LARGE_DBL);
  m_surf_tree.find(CubitBox(min, max), tri_list);

  for (int i = 0; i < tri_list.size(); i++)
    if (tri_list.next(i))
      delete tri_list.next(i);

  for (std::set<Vert*>::iterator it = m_all_points.begin();
      it != m_all_points.end(); ++it) {
    if ((*it))
      delete (*it);
  }
  for (std::vector<double*>::iterator it2 = m_all_curvatures.begin();
      it2 != m_all_curvatures.end(); ++it2) {
    if ((*it2))
      delete (*it2);
  }

  m_all_points.clear();
  m_all_curvatures.clear();

}

void
SurfIntersect::get_triangles(SurfTriSet& triangles)
{

  DLIList<SurfTri*> tri_list;
  CubitVector min(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
  CubitVector max( LARGE_DBL, LARGE_DBL, LARGE_DBL);
  m_surf_tree.find(CubitBox(min, max), tri_list);
  int tri_list_size = tri_list.size();

  for (int i = 0; i < tri_list_size; i++)
    triangles.insert(tri_list.get_and_step());

}

//In the base class, the intersect data is built using
//a surface mesh for graphics output. The most "modeler-independent"
//mesh available. Only problem is that it is given in the form
//of a polygon soup, so no connectivity or topological "affiliation".

//If better data is available, the derived classes will do a better
//job at this, especially to compute curvature info.

void
SurfIntersect::init_intersect(RefFace* const surface)
{

  GMem facet_container;
  surface->get_graphics(facet_container);

  int num_points = facet_container.point_list_size();
  GPoint* point_array = facet_container.point_list();

  Vert* point_coord[num_points];
  std::fill(point_coord, point_coord + num_points, static_cast<Vert*>(NULL));

  int num_facets = facet_container.fListCount;
  int* facet_array = facet_container.facet_list();

  //only used if m_has_curvature_data is true
  VertFaceMap vert_to_faces;

  for (int j = 0; j < num_facets;) {

    int num_verts = facet_array[j++];
    assert(num_verts == 3);
    Vert* verts[] =
      { NULL, NULL, NULL };

    for (int k = 0; k < num_verts; k++) {

      int index = facet_array[j++];

      if (!point_coord[index]) {

	Vert* this_vert = new Vert();
	double this_vert_coord[] =
	  { point_array[index].x, point_array[index].y, point_array[index].z };
	this_vert->setCoords(3, this_vert_coord);

	m_all_points.insert(this_vert);
	point_coord[index] = this_vert;

      }

      verts[k] = point_coord[index];

    }

    SurfTri* triangle = new SurfTri(surface, verts[0], verts[1], verts[2]);

    if (m_has_curvature_data) {

      list<SurfTri*> dummy_list;
      VertFaceMap::iterator it_vert;

      for (int i = 0; i < 3; i++) {
	it_vert = vert_to_faces.insert(make_pair(verts[i], dummy_list)).first;
	assert(it_vert != vert_to_faces.end());
	it_vert->second.push_back(triangle);
      }

    }

    m_surf_tree.add(triangle);

  }

  if (m_has_curvature_data)
    init_curvature(vert_to_faces);

}

//Note that the differential geometry operators only work on
//valid surface vertices i.e. vertices having a 1-star forming
//a topological circle.

//The boundary vertices are given a curvature of 0 for the moment
//(at the moment no need to implement, better support in derived classes).

void
SurfIntersect::init_curvature(const VertFaceMap& vert_to_faces)
{

  assert(!vert_to_faces.empty());

  vector<Vert*> bdry_verts, int_verts;
  VertFaceMap::const_iterator it = vert_to_faces.begin(), it_end =
      vert_to_faces.end();

  for (; it != it_end; ++it) {

    deque<Vert*> vertex_ring;

    Vert* vertex = it->first;
    list<SurfTri*> facets = it->second;
    assert(!facets.empty());
    list<SurfTri*>::iterator itf = facets.begin();

    while (!facets.empty()) {

      SurfTri* surf_triangle = *itf;
      assert(
	  surf_triangle->get_vertex(0) == vertex
	      || surf_triangle->get_vertex(1) == vertex
	      || surf_triangle->get_vertex(2) == vertex);

      Vert *vert1, *vert2;

      switch (surf_triangle->get_vertex(0) == vertex ? 1 :
	      surf_triangle->get_vertex(1) == vertex ? 2 : 3)
	{
	case 1:
	  vert1 = surf_triangle->get_vertex(1);
	  vert2 = surf_triangle->get_vertex(2);
	  break;
	case 2:
	  vert1 = surf_triangle->get_vertex(0);
	  vert2 = surf_triangle->get_vertex(2);
	  break;
	case 3:
	  vert1 = surf_triangle->get_vertex(0);
	  vert2 = surf_triangle->get_vertex(1);
	  break;
	}

      assert(vert1 != vertex && vert2 != vertex);

      if (vertex_ring.empty()) {
	vertex_ring.push_back(vert1);
	vertex_ring.push_back(vert2);
      }

      else {

	Vert *front_vert = vertex_ring.front(), *back_vert = vertex_ring.back();

	if (vert1 == front_vert) {
	  vertex_ring.push_front(vert2);
	}
	else if (vert1 == back_vert) {
	  vertex_ring.push_back(vert2);
	}
	else if (vert2 == front_vert) {
	  vertex_ring.push_front(vert1);
	}
	else if (vert2 == back_vert) {
	  vertex_ring.push_back(vert1);
	}
	else {
	  ++itf;
	  assert(itf != facets.end());
	  continue;
	}

      }

      facets.erase(itf);
      itf = facets.begin();

    }

    if (vertex_ring.front() == vertex_ring.back())
      int_verts.push_back(vertex);
    else
      bdry_verts.push_back(vertex);

  }

  //Computing the discrete curvature for interior vertices:
  //As this is designed for isotropic meshing,
  //only largest principal curvature considered.

  DiscreteCurvTool* DCT = DiscreteCurvTool::instance();

  vector<Vert*>::iterator itv = int_verts.begin(), itv_end = int_verts.end();
  list<SurfTri*>::const_iterator itf, itf_end;

  for (; itv != itv_end; ++itv) {

    Vert* vertex = *itv;

    it = vert_to_faces.find(vertex);
    assert(it != vert_to_faces.end());

    //As this is designed for isotropic meshing,
    //only largest principal curvature considered.
    DCT->init_tool(vertex, it->second);
    double* curvature = new double(DCT->principal_curvatures().first);

    itf = it->second.begin();
    itf_end = it->second.end();
    for (; itf != itf_end; ++itf)
      (*itf)->set_curvature_at_vert(vertex, curvature);
    delete[] curvature;
  }

  //Setting apex and curve vertices' curvature to 0 for now.

  itv = bdry_verts.begin();
  itv_end = bdry_verts.end();

  for (; itv != itv_end; ++itv) {

    Vert* vertex = *itv;

    it = vert_to_faces.find(vertex);
    assert(it != vert_to_faces.end());

    itf = it->second.begin();
    itf_end = it->second.end();
    for (; itf != itf_end; ++itf)
      (*itf)->set_curvature_at_vert(vertex, new double(0.));

  }

}

//The following function is used to intersect a line with a
//discrete surface. The initial version tried to split this
//line in multiple subsegment such that the number of faces
//returned by the tree was reasonable. It turns out that doing
//so ends up being more expensive because of the overhead
//resulting from the added complexity.

void
SurfIntersect::intersect_with_line(const CubitVector& beg_pt,
				   const CubitVector& end_pt,
				   IntersectSet& intersect)
{

  DLIList<SurfTri*> int_candidates;
  m_surf_tree.find(CubitBox(min_vec(beg_pt, end_pt), max_vec(beg_pt, end_pt)),
		   int_candidates);

  int size = int_candidates.size();

  if (size == 0)
    return;

  else {

    SurfTri* tri;
    CubitVector p0, p1, p2, int_pt;
    Vert *v0, *v1, *v2;

    for (int i = 0; i < size; i++) {

      tri = int_candidates.next(i);
      v0 = tri->get_vertex(0);
      v1 = tri->get_vertex(1);
      v2 = tri->get_vertex(2);

      p0 = CubitVector(v0->x(), v0->y(), v0->z());
      p1 = CubitVector(v1->x(), v1->y(), v1->z());
      p2 = CubitVector(v2->x(), v2->y(), v2->z());

      if (intersect_segment_triangle(beg_pt, end_pt, p0, p1, p2, int_pt)) {

	Intersection this_intersect;
	this_intersect.m_triangle = tri;
	this_intersect.m_intersect_coord = int_pt;

#ifndef NDEBUG
	if (m_has_curvature_data)
	  assert(
	      tri->get_curvature(0) && tri->get_curvature(1)
		  && tri->get_curvature(2));
#endif

	intersect.insert(this_intersect);

      }

    }

  }

}

void
SurfIntersect::intersect_with_voronoi_edge(const VoronoiEdge& vor_edge,
					   IntersectSet& intersect)
{

  intersect_with_line(vor_edge.beg_pt(), vor_edge.end_pt(), intersect);

}

bool
SurfIntersect::intersect_line_triangle(const CubitVector& l0,
				       const CubitVector& l1,
				       const CubitVector& p0,
				       const CubitVector& p1,
				       const CubitVector& p2,
				       const bool restrict_to_ray,
				       const bool restrict_to_segment,
				       CubitVector& intersection_coord)
{

  if (l0.about_equal(l1))
    return false;

  CubitVector u = p1 - p0;   //Parametric plane axis 1
  CubitVector v = p2 - p0;   //Parametric plane axis 2
  CubitVector n = u * v;     //Normal to plane (triangle)
  //n = n;
  if (iFuzzyComp(n.length_squared(), 0.) == 0)
    return false; //Degenerate

  CubitVector dir = l1 - l0; //Line direction

  double a = -n % (l0 - p0);
  double b = n % dir;

  if (iFuzzyComp(b, 0.) == 0)
    return false; //Line parallel to plane

  double r = a / b; //Intersection parameter along line;

  if (restrict_to_ray) {
    if (iFuzzyComp(r, 0.) == -1)
      return false;
  }
  else if (restrict_to_segment) {
    if (iFuzzyComp(r, 0.) == -1 || iFuzzyComp(r, 1.) == 1)
      return false;
  }

  CubitVector I = l0 + r * dir; //Intersection point between plane and line

  //Is this point inside triangle?
  //Use parametric plane equations to find out.

  CubitVector w = I - p0;

  double uu = u % u;
  double uv = u % v;
  double vv = v % v;
  double wu = w % u;
  double wv = w % v;
  double den = uv * uv - uu * vv;

  double s = (uv * wv - vv * wu) / den;
  if (iFuzzyComp(s, 0.) == -1 || iFuzzyComp(s, 1.) == 1)
    return false;

  double t = (uv * wu - uu * wv) / den;
  if (iFuzzyComp(t, 0.) == -1 || iFuzzyComp(s + t, 1.) == 1)
    return false;

  intersection_coord = I;

  return true;

}

bool
SurfIntersect::intersect_line_triangle(const CubitVector& l0,
				       const CubitVector& l1,
				       const CubitVector& p0,
				       const CubitVector& p1,
				       const CubitVector& p2,
				       CubitVector& intersection_coord)
{

  return intersect_line_triangle(l0, l1, p0, p1, p2, false, false,
				 intersection_coord);

}

bool
SurfIntersect::intersect_ray_triangle(const CubitVector& l0,
				      const CubitVector& l1,
				      const CubitVector& p0,
				      const CubitVector& p1,
				      const CubitVector& p2,
				      CubitVector& intersection_coord)
{

  return intersect_line_triangle(l0, l1, p0, p1, p2, true, false,
				 intersection_coord);

}

bool
SurfIntersect::intersect_segment_triangle(const CubitVector& l0,
					  const CubitVector& l1,
					  const CubitVector& p0,
					  const CubitVector& p1,
					  const CubitVector& p2,
					  CubitVector& intersection_coord)
{

  return intersect_line_triangle(l0, l1, p0, p1, p2, false, true,
				 intersection_coord);

}

bool
SurfIntersect::IntersectCompare::operator()(
    const Intersection& intersect1, const Intersection& intersect2) const
{

  int comp_x = iFuzzyComp(intersect1.m_intersect_coord.x(),
			  intersect2.m_intersect_coord.x());
  if (comp_x == -1)
    return true;
  else if (comp_x == 1)
    return false;
  else {
    int comp_y = iFuzzyComp(intersect1.m_intersect_coord.y(),
			    intersect2.m_intersect_coord.y());
    if (comp_y == -1)
      return true;
    else if (comp_y == 1)
      return false;
    else {
      int comp_z = iFuzzyComp(intersect1.m_intersect_coord.z(),
			      intersect2.m_intersect_coord.z());
      if (comp_z == -1)
	return true;
      else if (comp_z == 1)
	return false;
    }
  }

  return false;

}

SurfIntersect::FarthestPoint::FarthestPoint(CubitVector* const point,
					    CubitVector* const farthest,
					    double* curvature) :
    m_max_dist(-LARGE_DBL), m_point(point), m_farthest(farthest), m_curvature(
	curvature)
{
}

void
SurfIntersect::FarthestPoint::operator()(const Intersection& intersect)
{

  double dist = ((*m_point) - intersect.m_intersect_coord).length_squared();

  if (dist > m_max_dist) {
    m_farthest->set(intersect.m_intersect_coord);
    SurfTri* triangle = intersect.m_triangle;
    *m_curvature = triangle->curvature_at_coord(intersect.m_intersect_coord);
    m_max_dist = dist;
  }

}

SurfIntersect::FarthestIntersect::FarthestIntersect(
    const CubitVector& from_point, Intersection* const intersect) :
    m_max_dist(-LARGE_DBL), m_from_point(from_point), m_intersect(intersect)
{
}

void
SurfIntersect::FarthestIntersect::operator()(const Intersection& intersect)
{

  double dist = ((m_from_point) - intersect.m_intersect_coord).length_squared();
  if (dist > m_max_dist) {
    m_max_dist = dist;
    *m_intersect = intersect;
  }

}
