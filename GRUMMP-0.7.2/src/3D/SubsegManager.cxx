#include "GR_SubsegManager.h"
#include "GR_BaryCentricLengthScale.h"
#include "GR_misc.h"
#include "GR_CurveDiscretizer.h"
#include "GR_SmallAngleTool.h"
#include "GR_Subseg.h"
#include "GR_Vertex.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "KDDTree.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "GR_SurfMesh.h"
#include <deque>
#include <map>
#include <vector>
#include <iterator>
#include <functional>
#include <algorithm>
#include <utility>

using std::deque;
using std::map;
using std::vector;

SubsegManager::SubsegManager(BasicTopologyEntity * entity,
			     const double allowed_TVT,
			     const VolMesh * backgroundMesh,
			     CurveDiscretizer::CDType samplingType) :
    m_entity(entity), vertices(), subsegments(), vert_tree(
	new KDDTree<Vert*>()), subseg_tree(new KDDTree<Subseg*>()), apex_vert_map(), edge_subseg_map(), m_backgroundMesh(
	backgroundMesh), m_samplingType(samplingType)
{

  assert(entity);
  assert(
      strcmp(entity->class_name(), "Volume") == 0
	  || strcmp(entity->class_name(), "Surface") == 0);

  init_subsegs(allowed_TVT);

}
SubsegManager::SubsegManager(SurfMesh * surf, BasicTopologyEntity * entity) :
    m_entity(entity), vertices(), subsegments(), vert_tree(
	new KDDTree<Vert*>()), subseg_tree(new KDDTree<Subseg*>()), apex_vert_map(), edge_subseg_map(), m_backgroundMesh(
    NULL), m_samplingType(CurveDiscretizer::TANGENT)
{

  // lets create new verts in place of old ones, because its scary
  // to leave memory management to the original surf mesh
  // map from surf vert to newly created vert,
  // needed for connecting up subsegs
  std::map<Vert *, Vert *> vertMap;
  // handle vertices first.
  for (GR_index_t iV = 0; iV < surf->getNumVerts(); iV++) {
    Vert * oldVert = surf->getVert(iV);
    if (!oldVert->isDeleted() && oldVert->getVertType() != Vert::eInterior) {
      Vert * newVert = new Vert();
      newVert->setCoords(3, oldVert->getCoords());
      newVert->setType(oldVert->getVertType());
      newVert->setParentTopology(oldVert->getParentTopology());
      vertMap.insert(std::make_pair(oldVert, newVert));
    }
  }
  std::set<Subseg*> subsegSet;
  surf->getSubsegs(subsegSet);
  for (std::set<Subseg*>::iterator it = subsegSet.begin();
      it != subsegSet.end(); ++it) {
    Subseg * oldSubseg = *it;
    if (!oldSubseg->is_deleted()) {
      Vert * begVert = vertMap.find(oldSubseg->get_beg_vert())->second;
      Vert * endVert = vertMap.find(oldSubseg->get_end_vert())->second;
      if (oldSubseg->get_ref_edge()) {
	Subseg * newSubseg = new Subseg(begVert, endVert,
					oldSubseg->get_beg_param(),
					oldSubseg->get_end_param(),
					oldSubseg->get_ref_edge());
	subseg_tree->add(newSubseg);
	subsegments.push_back(newSubseg);
	EdgeSubsegMap::iterator esm_it = edge_subseg_map.find(
	    oldSubseg->get_ref_edge());
	if (esm_it != edge_subseg_map.end()) {
	  esm_it->second.push_back(newSubseg);
	}
	else {
	  vector<Subseg*> subsegs_of_ref_edge;
	  subsegs_of_ref_edge.push_back(newSubseg);
	  edge_subseg_map.insert(
	      std::make_pair(oldSubseg->get_ref_edge(), subsegs_of_ref_edge));
	}
      }
      else {
	if (begVert->getVertType() == Vert::eBdryCurve)
	  begVert->setParentTopology(oldSubseg->get_ref_face());
	if (endVert->getVertType() == Vert::eBdryCurve)
	  endVert->setParentTopology(oldSubseg->get_ref_face());
	Subseg * newSubseg = new Subseg(begVert, endVert,
					oldSubseg->get_ref_face());
	subseg_tree->add(newSubseg);
	subsegments.push_back(newSubseg);
      }
    }
  }

}

SubsegManager::~SubsegManager()
{

  if (vert_tree)
    delete vert_tree;
  if (subseg_tree)
    delete subseg_tree;

  vert_tree = static_cast<KDDTree<Vert*>*>(NULL);
  subseg_tree = static_cast<KDDTree<Subseg*>*>(NULL);

  for (std::vector<Vert*>::iterator it = vertices.begin(); it != vertices.end();
      ++it) {
    if ((*it))
      delete (*it);
  }
  for (std::vector<Subseg*>::iterator it = subsegments.begin();
      it != subsegments.end(); ++it) {
    if ((*it))
      delete (*it);
  }

  apex_vert_map.clear();
  edge_subseg_map.clear();

  vertices.clear();
  subsegments.clear();

}

void
SubsegManager::get_vertices(std::vector<Vert*>& verts) const
{

  std::remove_copy_if(vertices.begin(), vertices.end(),
		      std::back_inserter(verts),
		      std::mem_fun(&Vert::isDeleted));

}

void
SubsegManager::get_subsegments(vector<Subseg*>& subsegs) const
{

  std::remove_copy_if(subsegments.begin(), subsegments.end(),
		      std::back_inserter(subsegs),
		      std::mem_fun(&Subseg::is_deleted));

}

Vert*
SubsegManager::get_apex_vertex(RefVertex* const ref_vert) const
{

  ApexVertMap::const_iterator it = apex_vert_map.find(ref_vert);

  if (it == apex_vert_map.end())
    return static_cast<Vert*>(NULL);
  else {
    assert(
	dynamic_cast<RefVertex*>(it->second->getParentTopology()) == it->first);
    return it->second;
  }

}

void
SubsegManager::get_subsegments_of_edge(RefEdge* const ref_edge,
				       vector<Subseg*>& subsegs) const
{

  EdgeSubsegMap::const_iterator it = edge_subseg_map.find(ref_edge);

  if (it == edge_subseg_map.end())
    return;

  std::remove_copy_if(it->second.begin(), it->second.end(),
		      std::back_inserter(subsegs),
		      std::mem_fun(&Subseg::is_deleted));

}

Subseg*
SubsegManager::point_encroaching_on_subseg(const CubitVector& point,
					   const RefFace* const ref_face) const
{

  DLIList<Subseg*> return_candidates;
  subseg_tree->find(CubitBox(point, point), return_candidates);

  DLIList<RefFace*> ref_face_list;
  int size = return_candidates.size();

  for (int i = 0; i < size; i++) {

    Subseg* subseg = return_candidates.next(i);
    if (ref_face && subseg->get_ref_edge()) {
      ref_face_list.clean_out();
      subseg->get_ref_edge()->ref_faces(ref_face_list);
      if (ref_face_list.where_is_item(const_cast<RefFace*>(ref_face)) == -1)
	continue;
    }

    if (point_inside_or_on_subseg_ball(point, subseg))
      return subseg;

  }

  return static_cast<Subseg*>(NULL);

}

bool
SubsegManager::point_inside_subseg_ball(const CubitVector& point,
					const Subseg* const subseg)
{

  int comp_result = iFuzzyComp((point - subseg->mid_point()).length(),
			       subseg->ball_radius());

  return (comp_result == -1);

}

bool
SubsegManager::point_inside_or_on_subseg_ball(const CubitVector& point,
					      const Subseg* const subseg)
{

  int comp_result = iFuzzyComp((point - subseg->mid_point()).length(),
			       subseg->ball_radius());

  return (comp_result < 1);

}

bool
SubsegManager::vert_inside_subseg_ball(const Vert* const vertex,
				       const Subseg* const subseg)
{

  CubitVector point(vertex->x(), vertex->y(), vertex->z());

  return point_inside_subseg_ball(point, subseg);

}

bool
SubsegManager::vert_inside_or_on_subseg_ball(const Vert* const vertex,
					     const Subseg* const subseg)
{

  CubitVector point(vertex->x(), vertex->y(), vertex->z());

  return point_inside_or_on_subseg_ball(point, subseg);

}

void
SubsegManager::split_subseg(Subseg* const subseg_to_split,
			    vector<Vert*>* const new_verts,
			    vector<Subseg*>* const new_subsegs,
			    vector<Subseg*>* const old_subsegs)
{

  if (SmallAngleTool::instance()->subseg_has_small_angle(subseg_to_split)) {

    Vert* beg_vert = subseg_to_split->get_beg_vert();
    Vert* end_vert = subseg_to_split->get_end_vert();
    assert(beg_vert && end_vert);

    split_small_angle_subsegs(beg_vert, new_verts, new_subsegs, old_subsegs);
    split_small_angle_subsegs(end_vert, new_verts, new_subsegs, old_subsegs);

    return;

  }

  split_subseg_at_mid_param(subseg_to_split, false, new_verts, new_subsegs,
			    old_subsegs);

}

void
SubsegManager::init_subsegs(const double allowed_TVT)
{

  logMessage(1, "\nInitializing Subsegment Manager\n");

  create_apex_verts();
  create_edge_verts(allowed_TVT);

  if (subsegments.empty())
    return;

  if (m_samplingType == CurveDiscretizer::LENGTHSCALE)
    init_small_angle_tool(60.);
  else
    init_small_angle_tool();

  //	strongly_delaunayize();

  //Garbage collection: removing deleted entities from the object.
  //	purge();

  //Output some stats to screen.
  printf("  number of verts = %lu\n", vertices.size());
  printf("  number of subsegs = %lu\n", subsegments.size());

//	output_subsegs_to_file("subsegs.vtk");

}

void
SubsegManager::create_apex_verts()
{

  assert(apex_vert_map.empty());

  DLIList<RefVertex*> ref_vert_list;

  m_entity->ref_vertices(ref_vert_list);

  int size = ref_vert_list.size();

  for (int i = 0; i < size; i++) {

    RefVertex* ref_vert = ref_vert_list.next(i);
    CubitVector vert_coord = ref_vert->coordinates();

    Vert* vertex = new Vert();
    vertex->setType(Vert::eBdryApex);
    vertex->setCoords(3, vert_coord);
    vertex->setParentTopology(ref_vert);

    vert_tree->add(vertex);
    vertices.push_back(vertex);
    apex_vert_map.insert(std::make_pair(ref_vert, vertex));

  }

  assert(
      apex_vert_map.size() == static_cast<unsigned int>(ref_vert_list.size()));

}

void
SubsegManager::create_edge_verts(const double allowed_TVT)
{

  assert(allowed_TVT > 0.);
  assert(allowed_TVT < M_PI);

  assert(
      apex_vert_map.size()
	  == static_cast<unsigned int>(m_entity->num_ref_vertices()));

  DLIList<RefEdge*> ref_edge_list;
  m_entity->ref_edges(ref_edge_list);

  int num_ref_edges = ref_edge_list.size();
  for (int i = 0; i < num_ref_edges; i++) {
    //Compute discretization data for the ref_edge
    RefEdge* ref_edge = ref_edge_list.next(i);
    CurveDiscretizer::DiscreteCurve discrete_curve;
    CurveDiscretizer * CD;
    if (m_samplingType == CurveDiscretizer::LENGTHSCALE)
      CD = new LengthScaleDiscretizer(ref_edge, m_backgroundMesh);
    else if (m_samplingType == CurveDiscretizer::QUASISTRUCT)
      CD = new QuasiStructDiscretizer(ref_edge, m_backgroundMesh);
    else
      CD = new TangentDiscretizer(ref_edge, allowed_TVT);
    CD->discretize(discrete_curve);
    RefVertex* beg_ref_vert = ref_edge->start_vertex();
    RefVertex* end_ref_vert = ref_edge->end_vertex();

    assert(apex_vert_map.find(beg_ref_vert) != apex_vert_map.end());
    assert(apex_vert_map.find(end_ref_vert) != apex_vert_map.end());
    Vert* beg_vert = (apex_vert_map.find(beg_ref_vert))->second;
    Vert* end_vert = (apex_vert_map.find(end_ref_vert))->second;

    //Creating the vertex in order from start to end of ref_edge
    vector<Subseg::VertData> edge_verts;
    Subseg::VertData vert_data;

    vert_data.param = ref_edge->start_param();
    vert_data.vertex = beg_vert;
    edge_verts.push_back(vert_data);

    //Skips first and last vertices... initialized separately.
    CurveDiscretizer::DiscreteCurve::const_iterator it = discrete_curve.begin();
    CurveDiscretizer::DiscreteCurve::const_iterator it_end =
	--discrete_curve.end();

    while (++it != it_end) {

      double vert_param(it->get_param());
      CubitVector vert_coord(it->get_coord());

      Vert* vertex = new Vert();
      vertex->setCoords(3, vert_coord);
      vertex->setParentTopology(ref_edge);
      if (CD->get_type() == CurveDiscretizer::LENGTHSCALE) {
	const CurveDiscretizer::CurvePoint * cp = &(*it);
	if (dynamic_cast<LengthScaleDiscretizer*>(CD)->isIntBdryIntersection(
	    cp)) {
	  vertex->setType(Vert::eBdryTwoSide);
	}
	else
	  vertex->setType(Vert::eBdryCurve);
      }
      else
	vertex->setType(Vert::eBdryCurve);

      vertices.push_back(vertex);
      vert_tree->add(vertex);

      vert_data.param = vert_param;
      vert_data.vertex = vertex;
      edge_verts.push_back(vert_data);

    }

    vert_data.param = ref_edge->end_param();
    vert_data.vertex = end_vert;
    edge_verts.push_back(vert_data);
    assert(edge_verts.size() == discrete_curve.size());
    //Joining the curve verts and forming subsegments
    vector<Subseg::VertData>::iterator itv = edge_verts.begin();
    vector<Subseg::VertData>::iterator itv_end = --(edge_verts.end());

    vector<Subseg*> subsegs_of_ref_edge;

    for (; itv != itv_end;) {
      Subseg* subseg = new Subseg(*itv, *(++itv), ref_edge);

      subsegments.push_back(subseg);
      subsegs_of_ref_edge.push_back(subseg);
      subseg_tree->add(subseg);

    }

    //Initializing a entry for ref_edge in edge_subseg_map.
#ifndef NDEBUG
    assert(
	edge_subseg_map.insert(std::make_pair(ref_edge, subsegs_of_ref_edge)).second);
#else
    edge_subseg_map.insert( std::make_pair(ref_edge, subsegs_of_ref_edge) );
#endif
    delete CD;
  }
}

void
SubsegManager::init_small_angle_tool(const double minimum_angle)
{

  SmallAngleTool* small_angle_tool = SmallAngleTool::instance();
  //Initializing small angle data.
  small_angle_tool->clear_small_angle_data();
  small_angle_tool->set_minimum_angle(minimum_angle);
  small_angle_tool->initialize(subsegments);

  //First, splitting small angle subsegs bounded by two apex verts.
  std::set<Subseg*> multiple_entries;
  small_angle_tool->get_multiple_entries(multiple_entries);
  std::set<Subseg*>::iterator its = multiple_entries.begin();
  std::set<Subseg*>::iterator its_end = multiple_entries.end();

  for (; its != its_end; ++its) {
    split_subseg_at_mid_param(*its, true);

  }

#ifndef NDEBUG
  multiple_entries.clear();
  small_angle_tool->get_multiple_entries(multiple_entries);
  assert(multiple_entries.empty());
#endif
  //Splitting all subsegs part of a small angle complex
  vector<Vert*> all_verts;
  get_vertices(all_verts);
  vector<Vert*>::iterator itv = all_verts.begin();
  vector<Vert*>::iterator itv_end = all_verts.end();

  for (; itv != itv_end; ++itv) {
    split_small_angle_subsegs(*itv);
  }

}

void
SubsegManager::strongly_delaunayize(std::vector<Vert*>* verts)
{

#ifndef NDEBUG
  int small_angle_size = SmallAngleTool::instance()->size();
#endif

  //Initialize most variables here to save time
  //during the indented loops.

  int isize = 0, jsize = 0;

  deque<Vert*> verts_to_check;
  deque<Subseg*> subsegs_to_check;

  vector<Vert*> new_verts;
  vector<Subseg*> new_subsegs;

  DLIList<Vert*> vert_candidates;
  DLIList<Subseg*> subseg_candidates;

  Vert* vertex;
  Subseg* subseg;

  //The following depends on the status of verts.
  //If verts is NULL, then all vertices are checked.
  //If verts is empty, then no vertex is checked.
  //Otherwise, all the vertices of verts are checked.

  if (!verts) {
    std::copy(vertices.begin(), vertices.end(),
	      std::back_inserter(verts_to_check));
  }
  else {
    if (!verts->empty())
      std::copy(verts->begin(), verts->end(),
		std::back_inserter(verts_to_check));
    else
      return;
  }

  //The Delaunayization itself. None of the subsegs will have
  //any of verts_to_check inside their diamteral sphere at exit.

  while (!verts_to_check.empty()) {

    vertex = verts_to_check.front();
    verts_to_check.pop_front();
    assert(!vertex->isDeleted());

    subseg_candidates.clean_out();
    subseg_tree->find(vertex->bounding_box(), subseg_candidates);
    isize = subseg_candidates.size();

    for (int i = 0; i < isize; i++) {

      subseg = subseg_candidates.next(i);
      assert(!subseg->is_deleted());

      if (subseg->has_vert(vertex))
	continue;

      if (vert_inside_or_on_subseg_ball(vertex, subseg)) {

	verts_to_check.push_back(vertex);

	new_verts.clear();
	new_subsegs.clear();
	split_subseg(subseg, &new_verts, &new_subsegs);

	std::copy(new_verts.begin(), new_verts.end(),
		  std::back_inserter(verts_to_check));
	std::copy(new_subsegs.begin(), new_subsegs.end(),
		  std::back_inserter(subsegs_to_check));

	while (!subsegs_to_check.empty()) {

	  subseg = subsegs_to_check.front();
	  subsegs_to_check.pop_front();
	  if (subseg->is_deleted())
	    continue;

	  vert_candidates.clean_out();
	  vert_tree->find(subseg->bounding_box(), vert_candidates);
	  jsize = vert_candidates.size();

	  for (int j = 0; j < jsize; j++) {

	    vertex = vert_candidates.next(j);
	    assert(!vertex->isDeleted());

	    if (subseg->has_vert(vertex))
	      continue;

	    if (vert_inside_or_on_subseg_ball(vertex, subseg)) {

	      new_verts.clear();
	      new_subsegs.clear();
	      split_subseg(subseg, &new_verts, &new_subsegs);

	      std::copy(new_verts.begin(), new_verts.end(),
			std::back_inserter(verts_to_check));
	      std::copy(new_subsegs.begin(), new_subsegs.end(),
			std::back_inserter(subsegs_to_check));

	      break;

	    }

	  }

	}

	break;

      }

    }

  }

#ifndef NDEBUG
  assert(small_angle_size == SmallAngleTool::instance()->size());
#endif

#ifndef NDEBUG

  //Some debugging code that verifies if all subsegments are
  //strongly Delaunay as they should.

  vector<Subseg*>::iterator it = subsegments.begin();
  vector<Subseg*>::iterator it_end = subsegments.end();

  int num_missed = 0;

  for (; it != it_end; ++it) {

    subseg = *it;
    if (subseg->is_deleted())
      continue;

    vert_candidates.clean_out();
    vert_tree->find(subseg->bounding_box(), vert_candidates);
    isize = vert_candidates.size();

    for (int i = 0; i < isize; i++) {
      vertex = vert_candidates.next(i);
      if (subseg->has_vert(vertex))
	continue;
      if (vert_inside_or_on_subseg_ball(vertex, subseg))
	num_missed++;
    }

  }

  if (num_missed != 0) {
    printf("%d subsegments are not strongly Delaunay... aborting\n",
	   num_missed);
    abort();
  }

#endif

}

void
SubsegManager::split_small_angle_subsegs(Vert* const vertex,
					 vector<Vert*>* const new_verts,
					 vector<Subseg*>* const new_subsegs,
					 vector<Subseg*>* const old_subsegs)
{

  SmallAngleTool* small_angle_tool = SmallAngleTool::instance();

  vector<Subseg*> subsegs_to_split;
  small_angle_tool->get_subsegs_of_vert(vertex, subsegs_to_split);

  if (subsegs_to_split.empty())
    return;
  // Let's also actually mark these as small angle verts.
  vertex->markAsSmallAngleVert(true);

  //Get the split radius, then split the subsegs.
  double split_radius = small_angle_tool->compute_split_radius(
      subsegs_to_split);
  double split_param = -LARGE_DBL;

  vector<Subseg*>::iterator it = subsegs_to_split.begin();
  vector<Subseg*>::iterator it_end = subsegs_to_split.end();

  Subseg* subseg = NULL;

  for (; it != it_end; ++it) {

    subseg = *it;
    assert(!subseg->is_deleted());
    assert(subseg->get_beg_vert() != subseg->get_end_vert());

    if (subseg->get_beg_vert() == vertex)
      split_param = subseg->param_at_distance_from_beg(split_radius);
    else if (subseg->get_end_vert() == vertex)
      split_param = subseg->param_at_distance_from_end(split_radius);
    else
      assert(0);

    assert(
	iFuzzyComp(split_param, subseg->get_beg_param()) >= 0
	    && iFuzzyComp(split_param, subseg->get_end_param()) <= 0);

    split_subseg_at_param(subseg, split_param, true, new_verts, new_subsegs,
			  old_subsegs);

  }

}

void
SubsegManager::split_subseg_at_param(Subseg* const subseg_to_split,
				     const double split_param,
				     const bool small_angle_split,
				     vector<Vert*>* const new_verts,
				     vector<Subseg*>* const new_subsegs,
				     vector<Subseg*>* const old_subsegs)
{

  assert(!subseg_to_split->is_deleted());

  Vert* new_vert = new Vert(Vert::eUnknown);

  Subseg* new_subseg1 = new Subseg();
  Subseg* new_subseg2 = new Subseg();
  RefEdge* ref_edge = subseg_to_split->get_ref_edge();

  subseg_to_split->split_at_param(split_param, new_subseg1, new_subseg2,
				  new_vert);

  vertices.push_back(new_vert);
  vert_tree->add(new_vert);

  subsegments.push_back(new_subseg1);
  subsegments.push_back(new_subseg2);

  subseg_tree->add(new_subseg1);
  subseg_tree->add(new_subseg2);
  subseg_tree->remove(subseg_to_split);
  if (ref_edge) {
    EdgeSubsegMap::iterator it = edge_subseg_map.find(ref_edge);

    assert(it != edge_subseg_map.end());
    it->second.push_back(new_subseg1);
    it->second.push_back(new_subseg2);
  }
  if (small_angle_split) {
    SmallAngleTool* SAT = SmallAngleTool::instance();
    SAT->remove_subseg(subseg_to_split);
    SAT->add_subseg(new_subseg1);
    SAT->add_subseg(new_subseg2);
  }

  if (new_verts) {
    new_verts->push_back(new_vert);
  }

  if (new_subsegs) {
    new_subsegs->push_back(new_subseg1);
    new_subsegs->push_back(new_subseg2);
  }

  if (old_subsegs) {
    old_subsegs->push_back(subseg_to_split);
  }

  assert(subseg_to_split->is_deleted());
  assert(!SmallAngleTool::instance()->subseg_has_small_angle(subseg_to_split));

}

void
SubsegManager::split_subseg_at_mid_param(Subseg* const subseg_to_split,
					 const bool small_angle_split,
					 vector<Vert*>* const new_verts,
					 vector<Subseg*>* const new_subsegs,
					 vector<Subseg*>* const old_subsegs)
{

  double param1 = subseg_to_split->get_beg_param();
  double param2 = subseg_to_split->get_end_param();
  double split_param = 0.5 * (param1 + param2);

  split_subseg_at_param(subseg_to_split, split_param, small_angle_split,
			new_verts, new_subsegs, old_subsegs);

}

void
SubsegManager::purge_verts()
{

  //At the moment, this function is not useful since
  //no vertex ever gets deleted from the Subseg Manager.
  //However, function is here, just in case...

  Vert* vertex;
  vector<Vert*> undeleted_verts;

  vector<Vert*>::iterator itv = vertices.begin(), itv_end = vertices.end();

  for (; itv != itv_end; ++itv) {

    vertex = *itv;

    if (vertex->isDeleted()) {
      //An apex vert should never get deleted in the subseg manager.
      assert(vertex->getVertType() != Vert::eBdryApex);
      vert_tree->remove(vertex);
      delete vertex;
      continue;
    }

    undeleted_verts.push_back(vertex);

  }

  vertices.swap(undeleted_verts);

}

void
SubsegManager::purge_subsegs()
{

  vector<Subseg*> subseg_vector;

  //First, clearing the deleted subsegs from edge_subseg_map
  EdgeSubsegMap::iterator it = edge_subseg_map.begin(), it_end =
      edge_subseg_map.end();

  for (; it != it_end; ++it) {

    std::remove_copy_if(it->second.begin(), it->second.end(),
			std::back_inserter(subseg_vector),
			std::mem_fun(&Subseg::is_deleted));

    it->second.swap(subseg_vector);
    subseg_vector.clear();

  }

  //Finding deleted subsegs in the subsegments member vector.
  std::remove_copy_if(subsegments.begin(), subsegments.end(),
		      std::back_inserter(subseg_vector),
		      std::mem_fun(&Subseg::is_not_deleted));

  //Erasing deleted subsegs from subsegments.
  subsegments.erase(
      std::remove_if(subsegments.begin(), subsegments.end(),
		     std::mem_fun(&Subseg::is_deleted)),
      subsegments.end());

  //Deallocating the deleted subsegs.

  for (std::vector<Subseg*>::iterator its = subseg_vector.begin();
      its != subseg_vector.end(); ++its) {
    if ((*its))
      delete (*its);
  }

#ifndef NDEBUG

  //Deleted subsegs are removed from the tree as soon as
  //they are deleted. This is a sanity check to make sure
  //this truly happens.

  CubitVector min_vec(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
  CubitVector max_vec( LARGE_DBL, LARGE_DBL, LARGE_DBL);
  DLIList<Subseg*> tree_content;

  subseg_tree->find(CubitBox(min_vec, max_vec), tree_content);

  int size = tree_content.size();
  assert(size == static_cast<int>(subsegments.size()));

  for (int i = 0; i < size; i++)
    assert(tree_content.next(i)->is_not_deleted());

#endif

}

void
SubsegManager::output_verts_to_file(const char* const filename) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"SubsegManager::output_verts_to_file");

  vector<Vert*> all_verts;
  get_vertices(all_verts);

  vector<Vert*>::iterator it;
  vector<Vert*>::iterator it_beg = all_verts.begin();
  vector<Vert*>::iterator it_end = all_verts.end();
  unsigned int num_verts = it_end - it_beg;

  fprintf(out_file, "MeshVersionFormatted 1\n");
  fprintf(out_file, "Dimension 3\n");
  fprintf(out_file, "Vertices\n");
  fprintf(out_file, "%d\n", num_verts);
  for (it = it_beg; it != it_end; ++it) {
    Vert* vertex = *it;
    fprintf(out_file, "%lf %lf %lf 1\n", vertex->x(), vertex->y(), vertex->z());
  }
  fprintf(out_file, "End\n");

  fclose(out_file);

}

void
SubsegManager::output_subsegs_to_file(const char* const filename) const
{

  FILE *pFOutFile = NULL;

  pFOutFile = fopen(filename, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  vector<Vert*> all_verts;
  get_vertices(all_verts);

  vector<Vert*>::iterator itv;
  vector<Vert*>::iterator itv_beg = all_verts.begin();
  vector<Vert*>::iterator itv_end = all_verts.end();
  unsigned int num_verts = itv_end - itv_beg;

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", num_verts);
  for (itv = itv_beg; itv != itv_end; ++itv) {
    Vert* vertex = *itv;
    fprintf(pFOutFile, "%lf %lf %lf\n", vertex->x(), vertex->y(), vertex->z());
  }

  vector<Subseg*> all_subsegs;
  get_subsegments(all_subsegs);

  vector<Subseg*>::iterator its;
  vector<Subseg*>::iterator its_beg = all_subsegs.begin();
  vector<Subseg*>::iterator its_end = all_subsegs.end();
  unsigned int num_edges = its_end - its_beg;

  fprintf(pFOutFile, "CELLS %d %d\n", num_edges, num_edges * 3);
  for (its = its_beg; its != its_end; its++) {
    Subseg* subseg = *its;
    Vert* v1 = subseg->get_beg_vert();
    Vert* v2 = subseg->get_end_vert();
    unsigned int index1 = std::find(itv_beg, itv_end, v1) - itv_beg;
    unsigned int index2 = std::find(itv_beg, itv_end, v2) - itv_beg;
    fprintf(pFOutFile, "2 %d %d\n", index1, index2);
  }
  fprintf(pFOutFile, "CELL_TYPES %d\n", num_edges);

  for (its = its_beg; its != its_end; its++) {

    fprintf(pFOutFile, "3\n");
  }
  fclose(pFOutFile);
  //	Subseg* subseg;
  //		SubsegMapType::const_iterator it = m_map.begin(), it_end = m_map.end();
  //		SubsegSet::const_iterator its, its_end;
  //		GR_index_t nSubsegs = 0;
  //		Vert * pV;
  //		for( ; it != it_end; ++it) {
  //			its = it->second.begin();
  //			its_end = it->second.end();
  //			for( ; its != its_end; ++its) {
  //				nSubsegs++;
  //			}
  //		}
  //		FILE *pFOutFile = NULL;
  //
  //		pFOutFile = fopen("SubsegMap.vtk", "w");
  //		if (NULL == pFOutFile)
  //			vFatalError("Couldn't open output file for writing",
  //					"3d Vtk mesh output");
  //
  //		GR_index_t nno = nSubsegs*2;   // number of nodes
  //		//-------------------------------------------------------
  //		// Write the VTK header details
  //		//-------------------------------------------------------
  //
  //		fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  //		fprintf(pFOutFile, "GRUMMP Tetra example\n");
  //		fprintf(pFOutFile, "ASCII\n");
  //		fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  //		fprintf(pFOutFile, "POINTS %d float\n", nno);
  //		it = m_map.begin(), it_end = m_map.end();
  //		for( ; it != it_end; ++it) {
  //			its = it->second.begin();
  //			its_end = it->second.end();
  //			for( ; its != its_end; ++its) {
  //				subseg = *its;
  //				pV = subseg->get_beg_vert();
  //				fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
  //				pV = subseg->get_end_vert();
  //				fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
  //
  //			}
  //		}
  //		int counter = -1;
  //		fprintf(pFOutFile, "CELLS %d %d\n", nSubsegs,nSubsegs*3);
  //		it = m_map.begin(), it_end = m_map.end();
  //
  //		for( ; it != it_end; ++it) {
  //			its = it->second.begin();
  //			its_end = it->second.end();
  //			for( ; its != its_end; ++its) {
  //				fprintf(pFOutFile, "2 %d %d\n",counter+1,counter+2);
  //				counter = counter+2;
  //
  //			}
  //		}
  //
  //
  //
  //
  //		//-------------------------------------
  //		// write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  //		// VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  //		// VTK_PYRAMID = 14)
  //		//-------------------------------------
  //		fprintf(pFOutFile, "CELL_TYPES %d\n", nSubsegs);
  //		it = m_map.begin(), it_end = m_map.end();
  //		for( ; it != it_end; ++it) {
  //			its = it->second.begin();
  //			its_end = it->second.end();
  //			for( ; its != its_end; ++its) {
  //				fprintf(pFOutFile, "3\n");
  //			}
  //		}
  //		fclose(pFOutFile);
}
