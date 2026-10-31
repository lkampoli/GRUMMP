#include "GR_TriMeshBuilder.h"
#include "GR_BFace2D.h"
#include "GR_Geometry.h"
#include "GR_GRGeom2D.h"
#include "GR_GRCurve.h"
#include "GR_GRPoint.h"
#include "GR_Mesh2D.h"
#include "GR_Vertex.h"
#include "GR_VertexPair.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "Point.hpp"
#include "RTree.hpp"

#include <deque>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include <utility>

using std::deque;
using std::list;
using std::map;
using std::multimap;
using std::make_pair;
using std::pair;
using std::queue;
using std::set;
using std::vector;

TriMeshBuilder::TriMeshBuilder(const GRGeom2D* const geometry_to_mesh,
			       Mesh2D* const tri_mesh_to_build) :
    m_geometry(geometry_to_mesh), m_mesh(tri_mesh_to_build), m_WI(m_mesh), m_vert_tree(
	new RTree<Vert*>()), m_edge_tree(new RTree<VertexPair*>()), m_apex_vert_map(), m_boundary_edges()
{
}

TriMeshBuilder::~TriMeshBuilder()
{
  for (std::set<VertexPair*>::iterator it = m_boundary_edges.begin();
      it != m_boundary_edges.end(); ++it) {
    if ((*it))
      delete (*it);
  }

  m_boundary_edges.clear();
  m_apex_vert_map.clear();

  delete m_vert_tree;
  delete m_edge_tree;

}

bool
TriMeshBuilder::buildTerribleConstrainedDelaunayTriangulation(
    const CurveSamplingType sample_type)
{
  double xMin, yMin, xMax, yMax;
  findBoundingBox(xMin, yMin, xMax, yMax);
  createMeshInBox(*m_mesh, 2 * xMin - xMax, 2 * yMin - yMax, 2 * xMax - xMin,
		  2 * yMax - yMin);
//#ifndef NDEBUG
//			m_mesh->writeTempMesh();
//#endif
  //m_mesh->writeTempMesh();
  init_apex_verts();
  //m_mesh->writeTempMesh();
  //output_boundary_edges_to_file("temp1.mesh");

  init_curve_verts(sample_type);
  //m_mesh->writeTempMesh();
  //prepare_small_angles();

  //make_boundary_constrained_delaunay();
  //output_boundary_edges_to_file("temp.mesh");

  recover_boundary_flip_only();
  m_mesh->writeTempMesh();

  clean_up_domain();

  logMessage(0, "Finished building constrained Delaunay triangulation\n");

  return true;

}

bool
TriMeshBuilder::buildConstrainedDelaunayTriangulation(
    const CurveSamplingType sample_type)
{
  double xMin, yMin, xMax, yMax;
  findBoundingBox(xMin, yMin, xMax, yMax);
  createMeshInBox(*m_mesh, 2 * xMin - xMax, 2 * yMin - yMax, 2 * xMax - xMin,
		  2 * yMax - yMin);
//#ifndef NDEBUG
//			m_mesh->writeTempMesh();
//#endif
  init_apex_verts();
  //output_boundary_edges_to_file("temp1.mesh");

  init_curve_verts(sample_type);

  prepare_small_angles();

  make_boundary_constrained_delaunay();
  //output_boundary_edges_to_file("temp.mesh");
  recover_boundary_flip_only();

  clean_up_domain();

  logMessage(0, "Finished building constrained Delaunay triangulation\n");

  return true;

}

void
TriMeshBuilder::get_vertices(vector<Vert*>& vertices) const
{

  int num_verts = m_mesh->iNumVerts();

  for (int i = 0; i < num_verts; i++) {

    Vert* vertex = m_mesh->getVert(i);

    if (!vertex->isDeleted())
      vertices.push_back(vertex);

  }

}

void
TriMeshBuilder::get_vertices(deque<Vert*>& vertices) const
{

  int num_verts = m_mesh->iNumVerts();

  for (int i = 0; i < num_verts; i++) {

    Vert* vertex = m_mesh->getVert(i);

    if (!vertex->isDeleted())
      vertices.push_back(vertex);

  }

}

void
TriMeshBuilder::get_boundary_edges(
    std::vector<VertexPair*>& boundary_edges) const
{

  std::remove_copy_if(m_boundary_edges.begin(), m_boundary_edges.end(),
		      std::back_inserter(boundary_edges),
		      std::mem_fun(&VertexPair::is_deleted));

}

void
TriMeshBuilder::findBoundingBox(double& xMin, double& yMin, double& xMax,
				double& yMax) const
{

  xMin = LARGE_DBL;
  yMin = LARGE_DBL;
  xMax = -LARGE_DBL;
  yMax = -LARGE_DBL;

  list<GRPoint*> point_list;
  m_geometry->get_points(point_list);

  list<GRPoint*>::iterator itp = point_list.begin();
  list<GRPoint*>::iterator itp_end = point_list.end();

  for (; itp != itp_end; ++itp) {

    GRPoint* point = *itp;
    CubitVector coord = point->coordinates();

    xMin = std::min(xMin, coord.x());
    yMin = std::min(yMin, coord.y());
    xMax = std::max(xMax, coord.x());
    yMax = std::max(yMax, coord.y());

  }

  list<GRCurve*> curve_list;
  m_geometry->get_curves(curve_list);

  list<GRCurve*>::iterator itc = curve_list.begin();
  list<GRCurve*>::iterator itc_end = curve_list.end();

  for (; itc != itc_end; ++itc) {

    GRCurve* curve = *itc;

    DLIList<CubitVector*> int_ext;
    curve->get_curve_geom()->interior_extrema(int_ext);
    int size = int_ext.size();

    for (int i = 0; i < size; i++) {

      CubitVector* coord = int_ext.next(i);

      xMin = std::min(xMin, coord->x());
      yMin = std::min(yMin, coord->y());
      xMax = std::max(xMax, coord->x());
      yMax = std::max(yMax, coord->y());

      delete coord;

    }

  }
}

void
TriMeshBuilder::init_apex_verts()
{

  list<GRPoint*> point_list;
  m_geometry->get_points(point_list);

  list<GRPoint*>::iterator it = point_list.begin();
  list<GRPoint*>::iterator it_end = point_list.end();

  for (; it != it_end; ++it) {
    GRPoint* point = *it;
    insert_apex_vertex(point);
  }
}

void
TriMeshBuilder::init_curve_verts(const CurveSamplingType sample_type)
{

  list<GRCurve*> curve_list;
  m_geometry->get_curves(curve_list);

  list<GRCurve*>::iterator it = curve_list.begin();
  list<GRCurve*>::iterator it_end = curve_list.end();

  for (; it != it_end; ++it) {

    GRCurve* curve = *it;
    sample_curve(curve, sample_type);

  }

}

void
TriMeshBuilder::sample_curve(GRCurve* const curve,
			     CurveSamplingType sample_type)
{

  assert(curve);

  switch (sample_type)
    {

    case TVT_BASED:
      sample_curve_based_on_tvt(curve);
      break;

    default:
      vFatalError("Invalid sampling type",
		  "TriMeshBuilder::compute_curve_sampling()");
      break;

    }

}

void
TriMeshBuilder::sample_curve_based_on_tvt(GRCurve* const curve,
					  double max_tvt_allowed)
{

  double tvt = LARGE_DBL;

  GRPoint* beg_point = dynamic_cast<GRPoint*>(curve->get_start_point());
  GRPoint* end_point = dynamic_cast<GRPoint*>(curve->get_end_point());

  assert(beg_point && end_point);
  assert(m_apex_vert_map.find(beg_point) != m_apex_vert_map.end());
  assert(m_apex_vert_map.find(end_point) != m_apex_vert_map.end());

  Vert* beg_vert = m_apex_vert_map.find(beg_point)->second;
  Vert* end_vert = m_apex_vert_map.find(end_point)->second;

  VertexPair boundary_edge = VertexPair(curve, beg_vert, end_vert,
					curve->get_curve_geom()->min_param(),
					curve->get_curve_geom()->max_param());

  queue<VertexPair> edge_queue;
  edge_queue.push(boundary_edge);

  while (!edge_queue.empty()) {

    boundary_edge = edge_queue.front();

    tvt = boundary_edge.tangent_variation();

    if (iFuzzyComp(tvt, max_tvt_allowed) == 1) {
      VertexPair new_edge_1, new_edge_2;
      split_boundary_edge_at_mid_tvt(boundary_edge, new_edge_1, new_edge_2);
      edge_queue.push(new_edge_1);
      edge_queue.push(new_edge_2);
    }
    else {
      VertexPair* new_edge = new VertexPair(boundary_edge);
      m_boundary_edges.insert(new_edge);
      m_edge_tree->add(new_edge);
    }

    edge_queue.pop();

  }

}

void
TriMeshBuilder::prepare_small_angles()
{

  typedef pair<VertexPair*, VertexPair*> EdgePair;
  typedef pair<double, VertexPair*> AngleEdgePair;
  typedef map<GRPoint*, Vert*> ApexVertMap;

  //Every apex vertex potentially harbors one or many small angle complexes.
  //Find them and protect them by a shell.

  set<VertexPair*> small_angle_complex;

  ApexVertMap::const_iterator itv, itv_end = m_apex_vert_map.end();
  set<VertexPair*>::const_iterator ite, ite_end;

  for (itv = m_apex_vert_map.begin(); itv != itv_end; ++itv) {

    small_angle_complex.clear();

    //Find the boundary edges attached to apex_vert. Construct a map which
    //orders these edges as function of the plane angle formed by the tangent.

    set<AngleEdgePair> angle_edge_pairs;

    Vert* apex_vert = itv->second;
    CubitVector apex_coords(apex_vert->x(), apex_vert->y(), 0.);

    //Since edge to vert connectivity is not stored, we have to go through all the edges
    //to find those bounded by apex_vert. Fortunately, at this point, we have a small amount
    //of edges so this should not be too expensive.

    ite_end = m_boundary_edges.end();
    double min_length = LARGE_DBL, check_length;

    for (ite = m_boundary_edges.begin(); ite != ite_end; ++ite) {

      VertexPair* edge = *ite;
      if (edge->is_deleted())
	continue;

      if (edge->beg_vert() == apex_vert || edge->end_vert() == apex_vert) {
	check_length = edge->length_squared();
	min_length = check_length < min_length ? check_length : min_length;
	small_angle_complex.insert(edge);
      }

    }

    ite_end = small_angle_complex.end();
    double check_dist = 0.01 * sqrt(min_length);

    for (ite = small_angle_complex.begin(); ite != ite_end; ++ite) {

      VertexPair* edge = *ite;
      CubitVector check_coords;
      assert(edge->has_vert(apex_vert));
      assert(edge->beg_param() < edge->end_param());

      GRCurve* curve = edge->get_curve();

      if (edge->beg_vert() == apex_vert)
	curve->get_curve_geom()->coord_at_dist(check_dist, true, check_coords);
      else if (edge->end_vert() == apex_vert)
	curve->get_curve_geom()->coord_at_dist(check_dist, false, check_coords);
      else
	assert(0);

      double angle = atan2(check_coords.y() - apex_coords.y(),
			   check_coords.x() - apex_coords.x());

      angle_edge_pairs.insert(make_pair(angle, edge));

    }

    if (angle_edge_pairs.size() < 2)
      continue;

    //Add the first entry to the end of the map (with an offset of -2*M_PI)
    //to be able to complete a cycle around the apex_vertex when looking for
    //small angles.

    angle_edge_pairs.insert(
	make_pair(angle_edge_pairs.begin()->first + 2. * M_PI,
		  angle_edge_pairs.begin()->second));

    //Go through the map and identify pairs of boundary edges forming small angles.

    deque<EdgePair> edges_around_small_angle;

    set<AngleEdgePair>::iterator ita, ita_end = --angle_edge_pairs.end();

    for (ita = angle_edge_pairs.begin(); ita != ita_end;) {

      AngleEdgePair pair1 = *ita, pair2 = *(++ita);
      double angle1 = pair1.first, angle2 = pair2.first;

      //If the boundary edges are separated by an angle of more that 60 degrees,
      //then they cannot form a small angle. Otherwise, we need to protect the small angle.

      if ((angle2 - angle1) > M_PI / 3.)
	continue;

      //If the boundary edges do not share a common region where the small angle
      //is located, then they do not form a small angle.

      VertexPair *edge1 = pair1.second, *edge2 = pair2.second;
      int reg1 =
	  (edge1->beg_vert() == apex_vert) ?
	      edge1->get_curve()->region_left() :
	      edge1->get_curve()->region_right();
      int reg2 =
	  (edge2->beg_vert() == apex_vert) ?
	      edge2->get_curve()->region_right() :
	      edge2->get_curve()->region_left();

      if (reg1 != reg2)
	vFatalError("Two curves, sharing a small angle, "
		    "bound regions with different indices\n"
		    "There is probably an error in the boundary file.\n"
		    "I can't mesh the geometry. I need to abort\n",
		    "TriMeshBuilder::prepare_small_angles()");

      if (reg1 == 0)
	continue;

      //Found a pair of boundary edges forming a small angle. Insert them in
      //a small angle map (make sure no duplicate entries exist in the map).
      //Second entry is always on CW side of first entry.

      edges_around_small_angle.push_back(make_pair(edge1, edge2));

    }

    //edges_of_small_angle contains pairs of boundary edges sharing a small angle.
    //cross product of tangent(first) and tangent(second) is always positive.
    //Now time to form small angle complexes and split edges on concentric shells.
    while (!edges_around_small_angle.empty()) {

      small_angle_complex.clear();

      VertexPair* edge[] =
	{ static_cast<VertexPair*>(NULL), static_cast<VertexPair*>(NULL) };

      do {

	edge[0] = edges_around_small_angle.front().first;
	edge[1] = edges_around_small_angle.front().second;
	edges_around_small_angle.pop_front();
	if (!(edge[0]->is_not_deleted() && edge[1]->is_not_deleted()))
	  continue;
	assert(edge[0]->is_not_deleted() && edge[1]->is_not_deleted());

	small_angle_complex.insert(edge[0]);
	small_angle_complex.insert(edge[1]);

      }
      while (edge[1] == edges_around_small_angle.front().first);

      if (!small_angle_complex.empty()) {
	apex_vert->markAsSmallAngleVert(true);
	small_angle_split(apex_vert, small_angle_complex);
      }

    }

  }

}

void
TriMeshBuilder::make_boundary_constrained_delaunay()
{

  //Initialize most variables here to save time in loops.

  int isize = 0, jsize = 0;

  deque<Vert*> verts_to_check;
  deque<VertexPair*> edges_to_check;
  DLIList<Vert*> vert_candidates;
  DLIList<VertexPair*> edge_candidates;

  Vert* vertex;
  VertexPair* edge;
  EdgeSplitData ESD;

  get_vertices(verts_to_check);

  while (!verts_to_check.empty()) {

    vertex = verts_to_check.front();
    assert(!vertex->isDeleted());

    verts_to_check.pop_front();

    edge_candidates.clean_out();
    m_edge_tree->find(vertex->bounding_box(), edge_candidates);
    isize = edge_candidates.size();

    for (int i = 0; i < isize; i++) {

      edge = edge_candidates.next(i);
      assert(!edge->is_deleted());

      if (edge->has_vert(vertex))
	continue;

      if (vert_inside_or_on_edge_ball(vertex, edge)
	  && vert_is_visible_from_edge(vertex, edge)) {

	verts_to_check.push_front(vertex);

	ESD = split_boundary_edge(edge, VertexPair::MID_PARAM);

	verts_to_check.push_back(ESD.new_vertex);
	edges_to_check.push_back(ESD.new_edge1);
	edges_to_check.push_back(ESD.new_edge2);

	while (!edges_to_check.empty()) {

	  edge = edges_to_check.front();
	  edges_to_check.pop_front();

	  if (edge->is_deleted())
	    continue;

	  vert_candidates.clean_out();
	  m_vert_tree->find(edge->bounding_box(), vert_candidates);
	  jsize = vert_candidates.size();

	  for (int j = 0; j < jsize; j++) {

	    vertex = vert_candidates.next(j);
	    assert(!vertex->isDeleted());

	    if (edge->has_vert(vertex))
	      continue;

	    if (vert_inside_or_on_edge_ball(vertex, edge)
		&& vert_is_visible_from_edge(vertex, edge)) {

	      ESD = split_boundary_edge(edge, VertexPair::MID_TVT);

	      verts_to_check.push_back(ESD.new_vertex);
	      edges_to_check.push_back(ESD.new_edge1);
	      edges_to_check.push_back(ESD.new_edge2);

	      break;

	    }

	  }

	}

	break;

      }

    }

  }

}

bool
TriMeshBuilder::vert_inside_or_on_edge_ball(Vert* const vertex,
					    VertexPair* const edge) const
{

  CubitVector point(vertex->x(), vertex->y(), 0.);

  int comp_result = iFuzzyComp((point - edge->mid_point()).length(),
			       edge->radius());

  return (comp_result < 1);

}

bool
TriMeshBuilder::vert_is_visible_from_edge(Vert* const vertex,
					  VertexPair* const edge) const
{

  GRCurve* curve = edge->get_curve();

  //Is visibility occluded by an unmeshed region?

  if (curve->region_right() != 0 && curve->region_left() != 0)
    return true;

  else {

    int left_of_edge = checkOrient2D(edge->beg_vert(), edge->end_vert(),
				     vertex);
    if (left_of_edge == 0)
      return true; //right on the edge

    CubitVector vertex_coord(vertex->x(), vertex->y(), 0.), closest_coord;

    double closest_param = curve->get_curve_geom()->closest_param_on_curve(
	vertex_coord, edge->beg_param(), edge->end_param());
    curve->get_curve_geom()->coord_at_param(closest_param, closest_coord);

    CubitVector tangent_at_closest, from_closest_to_vertex(
	vertex_coord.x() - closest_coord.x(),
	vertex_coord.y() - closest_coord.y(), 0.);
    curve->get_curve_geom()->unit_tangent(closest_param, tangent_at_closest);

    double cross = tangent_at_closest.x() * from_closest_to_vertex.y()
	- tangent_at_closest.y() * from_closest_to_vertex.x();

    int left_of_curve = iFuzzyComp(cross, 0.);

    if (left_of_curve == 0)
      vFatalError("Your attempt at splitting an edge will insert a vertex"
		  "on a curve that is not the edge's parent."
		  "This suggests an error in the bondary file."
		  "I can't mesh that... I am aborting!",
		  "TriMeshBuilder::vert_is_visible_from_edge");

    if (left_of_edge != left_of_curve)
      return true;

    switch (left_of_curve)
      {

      case 1:
	if (curve->region_left() == 0)
	  return false;
	break;
      case -1:
	if (curve->region_right() == 0)
	  return false;
	break;
      default:
	assert(0);
	break;
      }

  }

  //If we got this far, it means that edge:
  // - has vertex in its diametral circle and
  // - has a meshed region on the side where vertex lies.

  //We need to test if visibility is occluded by other edges.

  CubitVector mini, maxi;

  Vert* vert[] =
    { edge->beg_vert(), edge->end_vert() };
  bool intersected[] =
    { false, false };

  mini.x(MIN3(vertex->x(), vert[0]->x(), vert[1]->x()));
  mini.y(MIN3(vertex->y(), vert[0]->y(), vert[1]->y()));
  mini.z(-1.e-12);
  maxi.x(MAX3(vertex->x(), vert[0]->x(), vert[1]->x()));
  maxi.y(MAX3(vertex->y(), vert[0]->y(), vert[1]->y()));
  maxi.z(1.e-12);

  CubitBox box(mini, maxi);

  DLIList<VertexPair*> occlude_candidates;
  m_edge_tree->find(box, occlude_candidates);
  int size = occlude_candidates.size();

  for (int i = 0; i < size; i++) {

    VertexPair* edge_test = occlude_candidates.next(i);

    if (edge_test == edge)
      continue;

    if (edge_test->beg_vert() == vertex || edge_test->end_vert() == vertex)
      continue;

    int intersect[] =
      { intersectSegments(vertex, vert[0], edge_test->beg_vert(),
			  edge_test->end_vert()), intersectSegments(
	  vertex, vert[1], edge_test->beg_vert(), edge_test->end_vert()) };

    if (intersect[0] == 1 && intersect[1] == 1)
      return false;

    for (int j = 0; j < 2; j++) {

      if (intersect[j] == 1) {

	assert(vert[j]->getVertType() != Vert::eInterior);

	if (!edge_test->has_vert(vert[j]))
	  intersected[j] = true;

	else if (vert[j]->getVertType() == Vert::eBdryApex) {

	  Vert* vert1 = vertex;

	  Vert* vert2 = vert[j];

	  Vert* vert3 = j == 0 ? vert[1] : vert[0];

	  Vert* vert4 =
	      edge_test->beg_vert() == vert[0] ?
		  edge_test->end_vert() : edge_test->beg_vert();

	  int orient1 = checkOrient2D(vert1, vert2, vert4);
	  int orient2 = checkOrient2D(vert2, vert3, vert4);
	  int orient3 = checkOrient2D(vert3, vert1, vert4);

	  if (orient1 == 1 && orient2 == 1 && orient3 == 1)
	    return false;
	  if (orient1 == -1 && orient2 == -1 && orient3 == -1)
	    return false;

	}

      }

    }

    if (intersected[0] && intersected[1])
      return false;

  }

  return true;

}

TriMeshBuilder::EdgeSplitData
TriMeshBuilder::split_boundary_edge(VertexPair* const edge_to_split,
				    const VertexPair::SplitType split_type,
				    const double shell_radius)
{

  assert(edge_to_split->is_not_deleted());

  EdgeSplitData ESD;
  ESD.new_edge1 = new VertexPair();
  ESD.new_edge2 = new VertexPair();
  ESD.new_vertex = m_mesh->getNewVert();

  if (split_type == VertexPair::INIT_SHELL)
    edge_to_split->initial_shell_split(shell_radius, ESD.new_edge1,
				       ESD.new_edge2, ESD.new_vertex);

  else
    edge_to_split->split(ESD.new_edge1, ESD.new_edge2, ESD.new_vertex,
			 split_type);

  set<Cell*> seed_guesses;
  find_seed_guesses(edge_to_split, seed_guesses);

  insert_vertex_in_mesh(ESD.new_vertex, &seed_guesses);
  ESD.new_vertex->setType(Vert::eBdry);

  edge_to_split->mark_deleted();
  m_boundary_edges.insert(ESD.new_edge1);
  m_boundary_edges.insert(ESD.new_edge2);

  m_edge_tree->remove(edge_to_split);
  m_edge_tree->add(ESD.new_edge1);
  m_edge_tree->add(ESD.new_edge2);

  return ESD;

}

void
TriMeshBuilder::split_boundary_edge_at_mid_tvt(VertexPair& split_edge,
					       VertexPair& new_edge1,
					       VertexPair& new_edge2)
{

  double new_param = split_edge.mid_tvt_param();

  CubitVector vert_coord;
  GRCurve* curve = split_edge.get_curve();
  curve->get_curve_geom()->coord_at_param(new_param, vert_coord);

  Vert* new_vertex = insert_curve_vertex(&split_edge, vert_coord);

  new_edge1.set_curve(curve);
  new_edge2.set_curve(curve);

  new_edge1.set_beg_vert(split_edge.beg_vert());
  new_edge1.set_end_vert(new_vertex);
  new_edge2.set_beg_vert(new_vertex);
  new_edge2.set_end_vert(split_edge.end_vert());

  new_edge1.set_beg_param(split_edge.beg_param());
  new_edge1.set_end_param(new_param);
  new_edge2.set_beg_param(new_param);
  new_edge2.set_end_param(split_edge.end_param());

}

void
TriMeshBuilder::small_angle_split(Vert* const apex_vertex,
				  const set<VertexPair*>& small_angle_complex)
{

  vector<VertexPair*> all_edges;
  std::copy(small_angle_complex.begin(), small_angle_complex.end(),
	    std::back_inserter(all_edges));

  double shell_radius = LARGE_DBL;

  VertexPair* edge = NULL;
  EdgeSplitData ESD;

  set<VertexPair*>::const_iterator it = small_angle_complex.begin(), it_end =
      small_angle_complex.end();

  for (; it != it_end; ++it) {

    edge = *it;
    assert(!edge->is_deleted());

    if (edge->has_two_apex_verts()) {

      ESD = split_boundary_edge(edge, VertexPair::MID_TVT);

      if (ESD.new_edge1->beg_vert() == apex_vertex
	  || ESD.new_edge1->end_vert() == apex_vertex)
	edge = ESD.new_edge1;
      else
	edge = ESD.new_edge2;

      all_edges.push_back(edge);

    }

    assert(!edge->is_deleted());

    shell_radius = std::min(edge->length_squared(), shell_radius);

  }

  vector<VertexPair*>::iterator itv = all_edges.begin(), itv_end =
      all_edges.end();

  for (; itv != itv_end; ++itv) {

    edge = *itv;
    if (edge->is_deleted())
      continue;

    assert(edge->beg_is_small() || edge->end_is_small());

    split_boundary_edge(edge, VertexPair::INIT_SHELL, 0.5 * sqrt(shell_radius));

  }

}

void
TriMeshBuilder::insert_vertex_in_mesh(Vert* const vertex,
				      const set<Cell*>* const seed_guesses)
{
  //TODO: Upgrade Mesh iterators so that the "if" and "else" code here can
  // collapse down to being identical and get factored.

  m_vert_tree->add(vertex);

  Cell* seed_cell = NULL;
  if (seed_guesses != NULL) {
    set<Cell*>::iterator iter = seed_guesses->begin(), iterEnd =
	seed_guesses->end();
    while (iter != iterEnd) {
      TriCell* cand = dynamic_cast<TriCell*>(*iter);
      assert(cand->isValid());
      if (cand->circumscribesPoint(vertex->getCoords())) {
	seed_cell = cand;
	break;
      }
      ++iter;
    }
  }
  else {
    EntContainer<TriCell>::iterator iter = m_mesh->triCell_begin(), iterEnd =
	m_mesh->triCell_end();
    for (; iter != iterEnd; ++iter) {
      TriCell *cand = &(*iter);
      assert(cand->isValid());
      if (cand->circumscribesPoint(vertex->getCoords())) {
	seed_cell = cand;
	break;
      }
    }
  }
  m_WI.computeHull(vertex->getCoords(), seed_cell);
  m_WI.insertPointInHull(vertex);
//#ifndef NDEBUG
//			m_mesh->writeTempMesh();
//#endif
}

Vert*
TriMeshBuilder::insert_vertex_in_mesh(const CubitVector& insert_location,
				      const set<Cell*>* const seed_guesses)
{

  Vert* vertex = m_mesh->createVert(insert_location.x(), insert_location.y());
  logMessage(2, "Inserting point at: %f %f\n", insert_location.x(),
	     insert_location.y());
  insert_vertex_in_mesh(vertex, seed_guesses);

  return vertex;

}

Vert*
TriMeshBuilder::insert_apex_vertex(GRPoint* const point)
{

  CubitVector insert_location = point->coordinates();

  Vert* vertex = insert_vertex_in_mesh(insert_location);

  vertex->setType(Vert::eBdryApex);
  vertex->setParentEntity(point);

  m_apex_vert_map.insert(make_pair(point, vertex));

  return vertex;

}

Vert*
TriMeshBuilder::insert_curve_vertex(VertexPair* const edge_getting_split,
				    const CubitVector& insert_location)
{

  std::set<Cell*> seed_guesses;
  find_seed_guesses(edge_getting_split, seed_guesses);

  Vert* vertex = insert_vertex_in_mesh(insert_location, &seed_guesses);

  vertex->setType(Vert::eBdry);
  vertex->setParentEntity(edge_getting_split->get_curve());

  return vertex;

}

void
TriMeshBuilder::find_seed_guesses(VertexPair* const edge_to_split,
				  set<Cell*>& guesses) const
{

  Vert* verts[] =
    { edge_to_split->beg_vert(), edge_to_split->end_vert() };

  std::set<Cell*> neigh_cells;
  std::set<Vert*> dummy_set;

  for (int i = 0; i < 2; i++) {

    findNeighborhoodInfo(verts[i], neigh_cells, dummy_set);
    std::copy(neigh_cells.begin(), neigh_cells.end(),
	      std::inserter(guesses, guesses.begin()));

  }

}

void
TriMeshBuilder::recover_boundary_flip_only()
{

  set<VertexPair*>::iterator it = m_boundary_edges.begin(), it_end =
      m_boundary_edges.end();

  for (; it != it_end; ++it) {

    VertexPair* boundary_edge = *it;

    if (boundary_edge->is_deleted())
      continue;

    Vert* vert1 = boundary_edge->beg_vert();
    Vert* vert2 = boundary_edge->end_vert();

#ifndef NDEBUG
    assert(m_mesh->recoverEdge(vert1, vert2));
    assert(findCommonFace(vert1, vert2, true) != pFInvalidFace);
#else
    m_mesh->recoverEdge(vert1, vert2);
#endif

  }

}

void
TriMeshBuilder::clean_up_domain()
{
  // This function combines the old init_boundary_entities and
  // delete_outside_of_domain functions, because creating a bdry edge for
  // a face that has two cells is ... tricky.

  //Marking all cells as being in an invalid region.
  int num_cells = m_mesh->getNumCells();

  for (int i = 0; i < num_cells; i++)
    m_mesh->getCell(i)->setRegion(iInvalidRegion);

  //Finding boundary faces and create BdryEdge or IntBdryEdge
  set<VertexPair*>::iterator it = m_boundary_edges.begin(), it_end =
      m_boundary_edges.end();

  for (; it != it_end; ++it) {

    VertexPair* edge = *it;
    if (edge->is_deleted())
      continue;

    GRCurve* curve = edge->get_curve();

    Vert* beg_vert = edge->beg_vert();
    Vert* end_vert = edge->end_vert();

    double beg_param = edge->beg_param();
    double end_param = edge->end_param();

    Face* face = findCommonFace(beg_vert, end_vert, true);

    if (face == pFInvalidFace) {
      cout << beg_vert->x() << ":" << beg_vert->y() << " " << end_vert->x()
	  << ":" << end_vert->y() << endl;
      vFatalError(
	  "Trying to build a boundary at a face not found in the mesh.\n"
	  "Boundary recovery probably failed.\n"
	  "I must stop here... sorry!\n",
	  "TriMeshBuilder::init_mesh_boundary()");
    }

    if (face->getVert(0) == beg_vert) {
      assert(face->getVert(1) == end_vert);
      face->getLeftCell()->setRegion(curve->region_left());
      face->getRightCell()->setRegion(curve->region_right());
    }
    else {
      assert(face->getVert(0) == end_vert);
      assert(face->getVert(1) == beg_vert);
      face->getLeftCell()->setRegion(curve->region_right()); // Yes, really
      face->getRightCell()->setRegion(curve->region_left()); // Yes, really
    }

    if (curve->region_left() == 0 || curve->region_right() == 0) {
      // Regular bdry edge
      m_mesh->createBdryEdge(beg_vert, end_vert, curve, beg_param, end_param);
    }
    else {
      // Internal bdry edge
      m_mesh->createIntBdryEdge(beg_vert, end_vert, curve, beg_param,
				end_param);
    }
  }

  //Tagging cells with its region.
  for (int i = 0; i < num_cells; i++) {

    Cell *cell = m_mesh->getCell(i);

    if (cell->isDeleted())
      continue;
    if (cell->getRegion() == iInvalidRegion)
      continue;

    m_mesh->markCleanNeighbors(cell);
  }

  //Deleting the cells lying outside of the domain as well as the four
  //bounding box vertices, the bounding box boundary faces and finally
  //the faces lying between two deleted cells.

  for (int i = 0; i < num_cells; i++) {
    Cell *cell = m_mesh->getCell(i);
    if (cell->isDeleted())
      continue;
    if (cell->getRegion() == iInvalidRegion
	|| cell->getRegion() == iOutsideRegion)
      m_mesh->deleteCell(cell);

  }

  //Mark farfield boundary faces for deletion. (the first four BFaces.)
  for (int i = 0; i < 4; i++)
    m_mesh->deleteBFace(m_mesh->getBFace(i));

  // Faces are deleted automagically at this point.

  // Mark four farfield vertices for deletion.
  for (int i = 0; i < 4; i++)
    m_mesh->deleteVert(m_mesh->getVert(i));

  purge_mesh();

}

void
TriMeshBuilder::purge_boundary_edges()
{

  set<VertexPair*> temp_set;

//  std::remove_copy_if(m_boundary_edges.begin(), m_boundary_edges.end(),
//                      std::inserter(temp_set, temp_set.end()),
//                      std::mem_fun(&VertexPair::is_deleted));
  for (std::set<VertexPair*>::iterator it = m_boundary_edges.begin();
      it != m_boundary_edges.end(); ++it) {
    if (!(*it)->is_deleted())
      temp_set.insert(*it);
    else
      delete (*it);

  }
  m_boundary_edges = temp_set;

}

void
TriMeshBuilder::purge_mesh()
{

  purge_boundary_edges();
  m_mesh->purgeAllEntities();

  assert(m_mesh->isValid());

}

void
TriMeshBuilder::output_boundary_edges_to_file(const char* const filename) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"TriMeshBuilder::output_verts_to_file");

  fprintf(out_file, "MeshVersionFormatted 1\n");
  fprintf(out_file, "Dimension 2\n");

  vector<Vert*> all_verts;
  get_vertices(all_verts);

  for (int i = 0; i < 4; i++)
    all_verts.erase(all_verts.begin());

  vector<Vert*>::iterator itv, itv_beg = all_verts.begin(), itv_end =
      all_verts.end();

  fprintf(out_file, "Vertices\n");
  fprintf(out_file, "%d\n", static_cast<int>(itv_end - itv_beg));

  for (itv = itv_beg; itv != itv_end; ++itv) {

    Vert* vertex = *itv;
    assert(!vertex->isDeleted());

    fprintf(out_file, "%lf %lf 1\n", vertex->x(), vertex->y());

  }

  vector<VertexPair*> all_edges;
  get_boundary_edges(all_edges);

  vector<VertexPair*>::iterator ite, ite_beg = all_edges.begin(), ite_end =
      all_edges.end();

  fprintf(out_file, "Edges\n");
  fprintf(out_file, "%d\n", static_cast<int>(ite_end - ite_beg));

  int color = 0;

  for (ite = ite_beg; ite != ite_end; ++ite) {

    VertexPair* edge = *ite;
    assert(!edge->is_deleted());

    Vert* vert1 = edge->beg_vert();
    Vert* vert2 = edge->end_vert();

    unsigned int index1 = std::find(itv_beg, itv_end, vert1) - itv_beg;
    unsigned int index2 = std::find(itv_beg, itv_end, vert2) - itv_beg;
    fprintf(out_file, "%d %d %d\n", index1 + 1, index2 + 1, color++);

  }

  fprintf(out_file, "End\n");
  fclose(out_file);
}

void
TriMeshBuilder::output_triangulation_to_file(const char* const filename) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"TriMeshBuilder::output_verts_to_file");

  fprintf(out_file, "MeshVersionFormatted 1\n");
  fprintf(out_file, "Dimension 2\n");

  vector<Vert*> all_verts;
  get_vertices(all_verts);

  vector<Vert*>::iterator itv, itv_beg = all_verts.begin(), itv_end =
      all_verts.end();

  fprintf(out_file, "Vertices\n");
  fprintf(out_file, "%d\n", static_cast<int>(itv_end - itv_beg));

  for (itv = itv_beg; itv != itv_end; ++itv) {

    Vert* vertex = *itv;
    assert(!vertex->isDeleted());

    fprintf(out_file, "%lf %lf 1\n", vertex->x(), vertex->y());

  }

  vector<Cell*> all_cells;
  int total_num_cells = m_mesh->getNumCells();

  for (int i = 0; i < total_num_cells; ++i) {

    Cell* cell = m_mesh->getCell(i);

    if (cell->isDeleted())
      continue;

    all_cells.push_back(m_mesh->getCell(i));

  }

  vector<Cell*>::iterator itc, itc_beg = all_cells.begin(), itc_end =
      all_cells.end();

  fprintf(out_file, "Triangles\n");
  fprintf(out_file, "%d\n", static_cast<int>(itc_end - itc_beg));

  for (itc = itc_beg; itc != itc_end; ++itc) {

    Cell* cell = *itc;
    assert(!cell->isDeleted());
    assert(cell->getNumVerts() == 3);

    int index1 = std::find(itv_beg, itv_end, cell->getVert(0)) - itv_beg;
    int index2 = std::find(itv_beg, itv_end, cell->getVert(1)) - itv_beg;
    int index3 = std::find(itv_beg, itv_end, cell->getVert(2)) - itv_beg;

    fprintf(out_file, "%d %d %d 1\n", index1 + 1, index2 + 1, index3 + 1);

  }

  fprintf(out_file, "End\n");

  fclose(out_file);

}

void
TriMeshBuilder::output_cells_to_file(const char* const filename,
				     const set<Cell*>& cells) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"TriMeshBuilder::output_verts_to_file");

  fprintf(out_file, "MeshVersionFormatted 1\n");
  fprintf(out_file, "Dimension 2\n");

  std::vector<Vert*> all_verts;

  set<Cell*>::iterator itc = cells.begin(), itc_end = cells.end();

  for (; itc != itc_end; ++itc) {

    Cell* cell = *itc;

    for (int i = 0; i < 3; i++)
      all_verts.push_back(cell->getVert(i));

  }

  std::vector<Vert*>::iterator itv, itv_end;
  std::sort(all_verts.begin(), all_verts.end());
  itv_end = std::unique(all_verts.begin(), all_verts.end());

  Vert* vertex = NULL;

  fprintf(out_file, "Vertices\n");
  fprintf(out_file, "%d\n", static_cast<int>(itv_end - all_verts.begin()));
  for (itv = all_verts.begin(); itv != itv_end; ++itv) {
    vertex = *itv;
    fprintf(out_file, "%lf %lf 1\n", vertex->x(), vertex->y());
  }

  fprintf(out_file, "Triangles\n");
  fprintf(out_file, "%d\n", static_cast<int>(cells.size()));

  for (itc = cells.begin(); itc != itc_end; ++itc) {

    Cell* cell = *itc;
    assert(!cell->isDeleted());
    assert(cell->getNumVerts() == 3);

    int index1 = std::find(all_verts.begin(), itv_end, cell->getVert(0))
	- all_verts.begin();
    int index2 = std::find(all_verts.begin(), itv_end, cell->getVert(1))
	- all_verts.begin();
    int index3 = std::find(all_verts.begin(), itv_end, cell->getVert(2))
	- all_verts.begin();

    fprintf(out_file, "%d %d %d 1\n", index1 + 1, index2 + 1, index3 + 1);

  }

  fprintf(out_file, "End\n");

  fclose(out_file);

}
