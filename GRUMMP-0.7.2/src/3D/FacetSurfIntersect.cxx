#include "GR_SurfIntersect.h"
#include "GR_DiscreteCurvTool.h"
#include "GR_FacetSurfIntersect.h"
#include "GR_misc.h"
#include "GR_Vertex.h"

#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"
#include "DLIList.hpp"
#include "FacetEvalTool.hpp"
#include "FacetCurve.hpp"
#include "FacetPoint.hpp"
#include "FacetSurface.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

#include <algorithm>
#include <deque>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <utility>

using std::deque;
using std::list;
using std::map;
using std::pair;
using std::set;

FacetSurfTri::FacetSurfTri(RefFace* const surface, CubitFacet* const facet,
			   Vert* const vert1, Vert* const vert2,
			   Vert* const vert3) :
    SurfTri(surface, vert1, vert2, vert3), m_facet(facet)
{
}

FacetSurfTri::FacetSurfTri(RefFace* const surface, CubitFacet* const facet,
			   Vert* const vert1, double* const curv1,
			   Vert* const vert2, double* const curv2,
			   Vert* const vert3, double* const curv3) :
    SurfTri(surface, vert1, curv1, vert2, curv2, vert3, curv3), m_facet(facet)
{
}

void
FacetSurfTri::project_to_surface(const CubitVector& coord,
				 CubitVector* const surf_coord,
				 CubitVector* const surf_normal) const
{

  assert(surf_coord);

  CubitVector area_coord, this_coord(coord);
  CubitBoolean outside, project_ok, normal_ok;
  double tolerance = 1.e-15;

  FacetEvalTool::facet_area_coordinate(m_facet, this_coord, area_coord);
  project_ok = FacetEvalTool::project_to_facet(m_facet, this_coord, area_coord,
					       *surf_coord, outside, tolerance);
  assert(project_ok == CUBIT_TRUE);
  if (surf_normal) {
    normal_ok = FacetEvalTool::eval_facet_normal(m_facet, area_coord,
						 *surf_normal);
    assert(normal_ok == CUBIT_TRUE);
  }

}

FacetSurfIntersect::FacetSurfIntersect() :
    SurfIntersect(false), m_facetsurf_verts()
{
}

FacetSurfIntersect::FacetSurfIntersect(RefFace* const surface) :
    SurfIntersect(true), m_facetsurf_verts()
{

  init_intersect(surface);

}

FacetSurfIntersect::FacetSurfIntersect(const set<RefFace*>& surfaces) :
    SurfIntersect(true), m_facetsurf_verts()
{

  std::for_each(
      surfaces.begin(), surfaces.end(),
      std::bind1st(std::mem_fun(&FacetSurfIntersect::init_intersect), this));

}

FacetSurfIntersect::~FacetSurfIntersect()
{

  for (std::set<Vert*>::iterator it = m_facetsurf_verts.begin();
      it != m_facetsurf_verts.end(); ++it) {
    if (m_all_points.count(*it) == 0)
      delete *it;
  }

}
void
FacetSurfIntersect::init_intersect(RefFace* const surface)
{

  set<CubitFacet*> all_facets;
  set<CubitFacetEdge*> all_edges;
  set<CubitPoint*> all_pts, bdry_pts, surf_pts, curv_pts, apex_pts;

  DLIList<CubitFacet*> facet_list;
  DLIList<CubitFacetEdge*> edge_list;
  DLIList<CubitPoint*> point_list;

  FacetSurface* facet_surface;
  FacetCurve* facet_curve;
  FacetPoint* facet_point;
  CubitPoint* cubit_point;

  facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());
  assert(facet_surface);

  //Get all the facets of this surface in a set.
  facet_surface->tris(facet_list);
  for (int i = 0; i < facet_list.size(); i++)
    all_facets.insert(facet_list.get_and_step());

  //Get all the points on this surface (including on curves and apexes).
  facet_surface->get_my_points(point_list);
  for (int i = 0; i < point_list.size(); i++)
    all_pts.insert(point_list.get_and_step());

  //Get all apex points in a set.
  DLIList<RefVertex*> my_points;
  surface->ref_vertices(my_points);

  for (int i = 0; i < my_points.size(); i++) {
    facet_point =
	dynamic_cast<FacetPoint*>(my_points.get_and_step()->get_point_ptr());
    assert(facet_point);
    cubit_point = facet_point->get_cubit_point();
    assert(cubit_point->is_feature());
    apex_pts.insert(cubit_point);
  }

  //Get all curve points and all edges in sets.
  DLIList<RefEdge*> my_curves;
  surface->ref_edges(my_curves);

  for (int i = 0; i < my_curves.size(); i++) {

    edge_list.clean_out();
    facet_curve =
	dynamic_cast<FacetCurve*>(my_curves.get_and_step()->get_curve_ptr());
    assert(facet_curve);
    facet_curve->get_facets(edge_list);
    for (int j = 0; j < edge_list.size(); j++)
      all_edges.insert(edge_list.get_and_step());

    point_list.clean_out();
    facet_curve->get_points(point_list);
    for (int j = 0; j < point_list.size(); j++) {
      cubit_point = point_list.get_and_step();
      if (!cubit_point->is_feature())
	curv_pts.insert(cubit_point);
    }

  }

  //From the three initialized sets, deduce interior surface points.
  std::set_union(apex_pts.begin(), apex_pts.end(), curv_pts.begin(),
		 curv_pts.end(), std::inserter(bdry_pts, bdry_pts.begin()));
  std::set_difference(all_pts.begin(), all_pts.end(), bdry_pts.begin(),
		      bdry_pts.end(),
		      std::inserter(surf_pts, surf_pts.begin()));

  assert(points_valid(apex_pts, curv_pts, surf_pts));

  map<CubitPoint*, Vert*> vert_map;
  map<CubitPoint*, double*> curv_map;

  init_surf_pts(surf_pts, curv_pts, apex_pts, all_facets, vert_map, curv_map);
  init_curv_pts(surf_pts, curv_pts, apex_pts, all_edges, vert_map, curv_map);
  init_apex_pts(surf_pts, curv_pts, apex_pts, all_pts, vert_map, curv_map);
  assert(vert_map.size() == curv_map.size());
  assert(
      vert_map.size() == apex_pts.size() + curv_pts.size() + surf_pts.size());
  init_point_handles(vert_map);
  init_curv_handles(curv_map);

  //Finally, initiate the surface triangles in the intersection tool.
  DLIList<CubitFacet*> facets;
  facet_surface->tris(facets);

  init_surf_tris(surface, facets, vert_map, curv_map);
//  for(set<CubitPoint*>::iterator it = surf_pts.begin(); it != surf_pts.end(); ++it)
//  	delete *it;
}

void
FacetSurfIntersect::init_surf_pts(const set<CubitPoint*>& surf_pts,
				  const set<CubitPoint*>& curv_pts,
				  const set<CubitPoint*>& apex_pts,
				  const set<CubitFacet*>& admissible_facets,
				  map<CubitPoint*, Vert*>& vert_map,
				  map<CubitPoint*, double*>& curv_map)
{

  if (surf_pts.empty())
    return;

  CubitPoint* my_point;
  DLIList<CubitFacet*> facet_list;
  DiscreteCurvTool* DCT = DiscreteCurvTool::instance();

  typedef map<CubitPoint*, Vert*> PointMap;
  PointMap point_map;
  PointMap::iterator it_map;

  set<CubitPoint*>::const_iterator it = surf_pts.begin(), it_end =
      surf_pts.end();

  for (; it != it_end; ++it) {

    my_point = *it;

    assert(vert_map.count(my_point) == 0);
    assert(curv_map.count(my_point) == 0);
    Vert * my_vert = new_vert_at_point(my_point);
    vert_map.insert(std::make_pair(my_point, my_vert));

    facet_list.clean_out();
    my_point->facets(facet_list);

#ifndef NDEBUG

    //checks if there is a valid ring around the point.

    deque<CubitPoint*> ring;
    DLIList<CubitFacet*> facet_list_copy(facet_list);

    while (facet_list_copy.size() != 0) {

      facet_list_copy.reset();
      int list_size = facet_list_copy.size();

      for (int i = 0; i < list_size; i++) {

	CubitFacet* facet = facet_list_copy.get();
	assert(admissible_facets.count(facet) == 1);

	CubitPoint* point0, *point1;
	facet->opposite_edge(my_point, point0, point1);

	assert(point0 && point1);
	assert(
	    surf_pts.count(point0) == 1 || curv_pts.count(point0) == 1
		|| apex_pts.count(point0) == 1);
	assert(
	    surf_pts.count(point1) == 1 || curv_pts.count(point1) == 1
		|| apex_pts.count(point1) == 1);

	if (ring.empty()) {
	  ring.push_back(point0);
	  ring.push_back(point1);
	}
	else {
	  CubitPoint *front = ring.front(), *back = ring.back();
	  if (point0 == front) {
	    ring.push_front(point1);
	  }
	  else if (point0 == back) {
	    ring.push_back(point1);
	  }
	  else if (point1 == front) {
	    ring.push_front(point0);
	  }
	  else if (point1 == back) {
	    ring.push_back(point0);
	  }
	  else {
	    assert(!facet_list_copy.is_at_end());
	    facet_list_copy.step();
	    continue;
	  }
	}

	facet_list_copy.remove();
	break;

      }

    }

    assert(facet_list_copy.size() == 0);
    assert(ring.front() == ring.back());

#endif

    point_map.clear();
    facet_list.reset();

    int list_size = facet_list.size();

    Vert* this_vert;
    CubitFacet* this_facet;
    CubitPoint* this_facet_pts[3];

    DiscreteCurvTool::VertTriplet triplet;
    list<DiscreteCurvTool::VertTriplet> triplet_list;

    for (int i = 0; i < list_size; i++) {

      this_facet = facet_list.get_and_step();
      this_facet_pts[0] = this_facet->point(0);
      this_facet_pts[1] = this_facet->point(1);
      this_facet_pts[2] = this_facet->point(2);

      for (int j = 0; j < 3; j++) {

	it_map = point_map.find(this_facet_pts[j]);

	if (it_map == point_map.end()) {
	  this_vert = new_vert_at_point(this_facet_pts[j]);
	  point_map.insert(std::make_pair(this_facet_pts[j], this_vert));
	}
	else {
	  this_vert = it_map->second;
	}

	triplet.set(j, this_vert);

      }

      triplet_list.push_back(triplet);

    }

    assert(point_map.count(my_point) == 1);
    DCT->init_tool(point_map.find(my_point)->second, triplet_list);
    curv_map.insert(
	std::make_pair(my_point,
		       new double(DCT->principal_curvatures().first)));

//     printf("surface curvature = %lf\n", DCT->principal_curvatures().first);

  }

}

void
FacetSurfIntersect::init_curv_pts(const set<CubitPoint*>& surf_pts,
				  const set<CubitPoint*>& curv_pts,
				  const set<CubitPoint*>& apex_pts,
				  const set<CubitFacetEdge*>& admissible_edges,
				  map<CubitPoint*, Vert*>& vert_map,
				  map<CubitPoint*, double*>& curv_map)
{

  if (curv_pts.empty())
    return;

  CubitPoint* point0, *point1, *point2;
  CubitFacetEdge *edge1, *edge2;
  DLIList<CubitFacetEdge*> my_edges;

  set<CubitPoint*>::const_iterator it = curv_pts.begin(), it_end =
      curv_pts.end();

  for (; it != it_end; ++it) {

    edge1 = edge2 = NULL;

    point0 = *it;
    assert(vert_map.count(point0) == 0);
    assert(curv_map.count(point0) == 0);

    my_edges.clean_out();
    point0->edges(my_edges);

    for (int i = 0; i < my_edges.size(); i++) {

      CubitFacetEdge* this_edge = my_edges.get_and_step();
      if (admissible_edges.count(this_edge) == 0)
	continue;

      if (this_edge->is_feature()) {
#ifndef NDEBUG
	if (!edge1) {
	  edge1 = this_edge;
	}
	else if (!edge2) {
	  edge2 = this_edge;
	}
	else {
	  assert(0);
	}
#else
	if(!edge1) {edge1 = this_edge;}
	else if(!edge2) {edge2 = this_edge; break;}
#endif	
      }

    }

    assert(edge1 && edge2);

    point1 = edge1->point(0) == point0 ? edge1->point(1) : edge1->point(0);
    point2 = edge2->point(0) == point0 ? edge2->point(1) : edge2->point(0);

    assert(
	surf_pts.count(point1) == 1 || curv_pts.count(point1) == 1
	    || apex_pts.count(point1) == 1);
    assert(
	surf_pts.count(point2) == 1 || curv_pts.count(point2) == 1
	    || apex_pts.count(point2) == 1);

    vert_map.insert(std::make_pair(point0, new_vert_at_point(point0)));

    //circumradius of triangle (p0,p1,p2) = good estimate of radius of curvature.
    //code lifted from TriFace class.
    double coord0[3], coord1[3], coord2[3];
    point0->coordinates(coord0);
    point1->coordinates(coord1);
    point2->coordinates(coord2);

    double vec1[] =
      { coord1[0] - coord0[0], coord1[1] - coord0[1], coord1[2] - coord0[2] };
    double vec2[] =
      { coord2[0] - coord0[0], coord2[1] - coord0[1], coord2[2] - coord0[2] };
    double vec3[3], circ_center[3];
    vCROSS3D(vec1, vec2, vec3);

    if (iFuzzyComp(dMAG3D_SQ(vec3), 0.) == 0) { //Three points are colinear
      curv_map.insert(std::make_pair(point0, new double(0.)));
      continue;
    }

    double dot1 = 0.5 * dDOT3D(vec1, vec1);
    double dot2 = 0.5 * dDOT3D(vec2, vec2);
    solve3By3(vec1, vec2, vec3, dot1, dot2, 0., circ_center);

    circ_center[0] += coord0[0];
    circ_center[1] += coord0[1];
    circ_center[2] += coord0[2];

    double circ_radius = dDIST3D(coord0, circ_center);
    assert(iFuzzyComp(circ_radius, 0.) == 1);
    double curvature = 1. / circ_radius;
    //double curvature = 0.;
    curv_map.insert(std::make_pair(point0, new double(curvature)));

//     double pt_coord[3];
//     point0->coordinates(pt_coord);
//     printf("curve curvature at %lf %lf %lf = %e\n", 
// 	   pt_coord[0], pt_coord[1], pt_coord[2], curvature);

  }

}

void
FacetSurfIntersect::init_apex_pts(const set<CubitPoint*>& surf_pts,
				  const set<CubitPoint*>& curv_pts,
				  const set<CubitPoint*>& apex_pts,
				  const set<CubitPoint*>& admissible_pts,
				  map<CubitPoint*, Vert*>& vert_map,
				  map<CubitPoint*, double*>& curv_map)
{

  if (apex_pts.empty())
    return;

  CubitPoint* my_point, *other_point;
  DLIList<CubitPoint*> neigh_points;
  map<CubitPoint*, double*>::iterator it_pt, itp;

  set<CubitPoint*>::const_iterator it = apex_pts.begin(), it_end =
      apex_pts.end();

  for (; it != it_end; ++it) {

    my_point = *it;
    assert(my_point->is_feature());
    assert(vert_map.count(my_point) == 0);
    assert(curv_map.count(my_point) == 0);
    assert(admissible_pts.count(my_point) == 1);
    vert_map.insert(std::make_pair(my_point, new_vert_at_point(my_point)));
    my_point->adjacent_points(neigh_points);

    double max_curv = 0.;

    for (int i = 0; i < neigh_points.size(); i++) {

      other_point = neigh_points.get_and_step();

      //If other_point is not in any of the sets, then it is on another RefFace.
      //If it is not already in curv_map, then it is an apex point having curvature = 0.
      if (surf_pts.count(other_point) == 1 || curv_pts.count(other_point) == 1
	  || apex_pts.count(other_point) == 1) {

	it_pt = curv_map.find(other_point);
	if (it_pt != curv_map.end())
	  max_curv = std::max(max_curv, *(it_pt->second));
      }

    }

    curv_map.insert(std::make_pair(my_point, new double(max_curv)));

//     double pt_coord[3];
//     my_point->coordinates(pt_coord);
//     printf("apex curvature at %lf %lf %lf = %e\n", 
// 	   pt_coord[0], pt_coord[1], pt_coord[2], min_curv);    

  }

}

void
FacetSurfIntersect::init_point_handles(const map<CubitPoint*, Vert*>& vert_map)
{

  map<CubitPoint*, Vert*>::const_iterator it = vert_map.begin(), it_end =
      vert_map.end();

  for (; it != it_end; ++it)
    m_all_points.insert(it->second);

}

void
FacetSurfIntersect::init_curv_handles(const map<CubitPoint*, double*>& curv_map)
{

  map<CubitPoint*, double*>::const_iterator it = curv_map.begin(), it_end =
      curv_map.end();

  for (; it != it_end; ++it)
    m_all_curvatures.push_back(it->second);

}

void
FacetSurfIntersect::init_surf_tris(RefFace* const surface,
				   DLIList<CubitFacet*>& facets,
				   const map<CubitPoint*, Vert*>& vert_map,
				   const map<CubitPoint*, double*>& curv_map)
{

  int num_tris = facets.size();
  CubitFacet* my_facet;
  CubitPoint* point0, *point1, *point2;

  facets.reset();

  for (int i = 0; i < num_tris; i++) {

    my_facet = facets.get_and_step();

    point0 = my_facet->point(0);
    point1 = my_facet->point(1);
    point2 = my_facet->point(2);
    assert(vert_map.count(point0) == 1 && curv_map.count(point0) == 1);
    assert(vert_map.count(point1) == 1 && curv_map.count(point1) == 1);
    assert(vert_map.count(point2) == 1 && curv_map.count(point2) == 1);

    FacetSurfTri* surf_tri = new FacetSurfTri(surface, my_facet,
					      vert_map.find(point0)->second,
					      curv_map.find(point0)->second,
					      vert_map.find(point1)->second,
					      curv_map.find(point1)->second,
					      vert_map.find(point2)->second,
					      curv_map.find(point2)->second);

    m_surf_tree.add(surf_tri);

  }

}

Vert*
FacetSurfIntersect::new_vert_at_point(CubitPoint* const my_point)
{

  Vert* my_vert = new Vert();
  m_facetsurf_verts.insert(my_vert);
  double my_coord[3];
  my_point->coordinates(my_coord);
  my_vert->setCoords(3, my_coord);

  return my_vert;

}

void
FacetSurfIntersect::get_points(DLIList<CubitPoint*>& the_points,
			       set<CubitPoint*>& apex_pts) const
{

  CubitPoint* my_point;
  the_points.reset();

  for (int i = 0; i < the_points.size(); i++) {
    my_point = the_points.get_and_step();
    if (my_point->is_feature())
      apex_pts.insert(my_point);
  }

}

void
FacetSurfIntersect::get_points(DLIList<CubitFacetEdge*>& the_edges,
			       set<CubitPoint*>& curve_pts) const
{

  CubitPoint* my_point0, *my_point1;
  CubitFacetEdge* my_edge;
  the_edges.reset();

  for (int i = 0; i < the_edges.size(); i++) {

    my_edge = the_edges.get_and_step();
    if (!my_edge->is_feature())
      continue;

    my_point0 = my_edge->point(0);
    my_point1 = my_edge->point(1);
    if (!my_point0->is_feature()) {
      curve_pts.insert(my_point0);
    }
    if (!my_point1->is_feature()) {
      curve_pts.insert(my_point1);
    }

  }

}

void
FacetSurfIntersect::get_points(DLIList<CubitFacet*>& the_facets,
			       set<CubitPoint*>& surface_pts) const
{

  the_facets.reset();

  CubitPoint* my_point0, *my_point1, *my_point2;
  CubitFacetEdge* my_edge0, *my_edge1, *my_edge2;
  CubitFacet* my_facet;

  set<CubitPoint*> all_pts;
  set<CubitPoint*> exclude_pts;

  //Have to proceed by inserting all the points and then 
  //removing those that are considered feature points.

  for (int i = 0; i < the_facets.size(); i++) {

    my_facet = the_facets.get_and_step();

    my_point0 = my_facet->point(0);
    my_point1 = my_facet->point(1);
    my_point2 = my_facet->point(2);
    all_pts.insert(my_point0);
    all_pts.insert(my_point1);
    all_pts.insert(my_point2);
    if (my_point0->is_feature())
      exclude_pts.insert(my_point0);
    if (my_point1->is_feature())
      exclude_pts.insert(my_point1);
    if (my_point2->is_feature())
      exclude_pts.insert(my_point2);

    my_edge0 = my_facet->edge(0);
    my_edge1 = my_facet->edge(1);
    my_edge2 = my_facet->edge(2);
    if (my_edge0->is_feature()) {
      exclude_pts.insert(my_edge0->point(0));
      exclude_pts.insert(my_edge0->point(1));
    }
    if (my_edge1->is_feature()) {
      exclude_pts.insert(my_edge1->point(0));
      exclude_pts.insert(my_edge1->point(1));
    }
    if (my_edge2->is_feature()) {
      exclude_pts.insert(my_edge2->point(0));
      exclude_pts.insert(my_edge2->point(1));
    }

  }

  std::set_difference(all_pts.begin(), all_pts.end(), exclude_pts.begin(),
		      exclude_pts.end(),
		      std::inserter(surface_pts, surface_pts.begin()));

}

bool
FacetSurfIntersect::points_valid(const set<CubitPoint*>& apex_pts,
				 const set<CubitPoint*>& curve_pts,
				 const set<CubitPoint*>& surface_pts) const
{

  set<CubitPoint*> dummy;

  std::set_intersection(apex_pts.begin(), apex_pts.end(), curve_pts.begin(),
			curve_pts.end(), std::inserter(dummy, dummy.begin()));

  if (!dummy.empty())
    return false;

  std::set_intersection(apex_pts.begin(), apex_pts.end(), surface_pts.begin(),
			surface_pts.end(), std::inserter(dummy, dummy.begin()));

  if (!dummy.empty())
    return false;

  std::set_intersection(curve_pts.begin(), curve_pts.end(), surface_pts.begin(),
			surface_pts.end(), std::inserter(dummy, dummy.begin()));

  if (!dummy.empty())
    return false;

  return true;

}
