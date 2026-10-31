#include "GR_SmallAngleTool.h"
#include "GR_Subseg.h"
#include "GR_Vertex.h"
#include "DLIList.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

using std::inserter;
using std::make_pair;
using std::pair;
using std::set;
using std::vector;

SmallAngleTool* SmallAngleTool::_instance_ = NULL;

SmallAngleTool*
SmallAngleTool::instance()
{

  if (_instance_ == NULL)
    _instance_ = new SmallAngleTool();

  return _instance_;

}

SmallAngleTool::SmallAngleTool() :
    vert_map(), minimum_angle(90.)
{
}

SmallAngleTool::~SmallAngleTool()
{

  vert_map.clear();

  delete _instance_;
  _instance_ = NULL;

}

double
SmallAngleTool::compute_angle(const Vert* const vertex,
			      const Subseg* const subseg1,
			      const Subseg* const subseg2) const
{

#ifndef NDEBUG
  Vert* v1 = subseg1->get_beg_vert();
  Vert* v2 = subseg1->get_end_vert();
  assert(v1 == vertex || v2 == vertex);
  v1 = subseg2->get_beg_vert();
  v2 = subseg2->get_end_vert();
  assert(v1 == vertex || v2 == vertex);
#endif

  CubitVector vert_coord(vertex->x(), vertex->y(), vertex->z());
  CubitVector tangent1(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
  CubitVector tangent2(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);

  RefVertex* ref_vert = dynamic_cast<RefVertex*>(vertex->getParentTopology());
  RefEdge* edge1 = subseg1->get_ref_edge();
  RefEdge* edge2 = subseg2->get_ref_edge();
  edge1->tangent(vert_coord, tangent1);
  edge2->tangent(vert_coord, tangent2);

  if (edge1 == edge2) {
    assert(edge1->is_periodic() && edge2->is_periodic());
    assert(edge1->start_vertex() == edge2->start_vertex());
    assert(edge1->end_vertex() == edge2->end_vertex());
    tangent2 = -tangent2;
  }
  else {
    if (edge1->start_vertex() != ref_vert)
      tangent1 = -tangent1;
    if (edge2->start_vertex() != ref_vert)
      tangent2 = -tangent2;
  }

  return tangent1.interior_angle(tangent2);

}

void
SmallAngleTool::remove_acceptable_angle_entries()
{

  if (vert_map.empty())
    return;

  typedef pair<ItMap, bool> MyPair;

  for (ItMap it = vert_map.begin(); it != vert_map.end();) {

    Vert* vertex = it->first;
    vector<MyPair> to_check(0);

    do {
      to_check.push_back(make_pair(it, false));
      if (++it == vert_map.end())
	break;
    }
    while (it->first == vertex);

    for (unsigned int i = 0; i < to_check.size(); i++) {
      Subseg* subseg1 = to_check[i].first->second;
      for (unsigned int j = i + 1; j < to_check.size(); j++) {
	if (to_check[i].second && to_check[j].second)
	  continue;
	Subseg* subseg2 = to_check[j].first->second;
	if (compute_angle(vertex, subseg1, subseg2) < minimum_angle) {
	  to_check[i].second = true;
	  to_check[j].second = true;
	}
      }
    }

    for (unsigned int i = 0; i < to_check.size(); i++)
      if (!to_check[i].second)
	vert_map.erase(to_check[i].first);

  }

}

void
SmallAngleTool::set_minimum_angle(const double angle)
{

  assert(angle > 0.);
//  assert(angle < 180.);
  minimum_angle = angle;

}

void
SmallAngleTool::initialize(const vector<Subseg*>& all_subsegs)
{

  //Clear the contents of the map
  vert_map.clear();

  //Add all the subsegments attached to a RefVertex
  for (unsigned int i = 0; i < all_subsegs.size(); i++) {

    Subseg* subsegment = all_subsegs[i];
    RefVertex* refvert1 = subsegment->get_beg_ref_vert();
    RefVertex* refvert2 = subsegment->get_end_ref_vert();

    if (refvert1) {
      vert_map.insert(make_pair(subsegment->get_beg_vert(), subsegment));
    }
    if (refvert2) {
      vert_map.insert(make_pair(subsegment->get_end_vert(), subsegment));
    }

  }

  //Remove all the subsegments that are not part of a small angle complex.
  remove_acceptable_angle_entries();

}

void
SmallAngleTool::get_multiple_entries(set<Subseg*>& multiple_entries) const
{

  for (ItcMap it = vert_map.begin(); it != vert_map.end(); it++) {
    Subseg* subseg = it->second;
    if (subseg->get_beg_ref_vert() && subseg->get_end_ref_vert())
      multiple_entries.insert(subseg);
  }

//  std::sort(multiple_entries.begin(), multiple_entries.end());
//  vector<Subseg*>::iterator it_end
//    = std::unique(multiple_entries.begin(), multiple_entries.end());
//  multiple_entries.erase(it_end, multiple_entries.end());

}

void
SmallAngleTool::add_subseg(Subseg* const subseg, const bool check_angle)
{

  RefVertex* refvert1 = subseg->get_beg_ref_vert();
  RefVertex* refvert2 = subseg->get_end_ref_vert();

  assert(!(refvert1 && refvert2));

  Vert* vertex = NULL;
  if (refvert1)
    vertex = subseg->get_beg_vert();
  else if (refvert2)
    vertex = subseg->get_end_vert();
  else
    return;

  if (vert_map.count(vertex) == 0)
    return;

  if (check_angle) {
    vector<Subseg*> subsegs_of_vert;
    get_subsegs_of_vert(vertex, subsegs_of_vert);

    for (unsigned int i = 0; i < subsegs_of_vert.size(); i++) {
      if (compute_angle(vertex, subseg, subsegs_of_vert[i]) < minimum_angle) {
	vert_map.insert(make_pair(vertex, subseg));
	return;
      }
    }
    return;
  }
  else
    vert_map.insert(make_pair(vertex, subseg));

}

void
SmallAngleTool::remove_subseg(Subseg* const subseg)
{

  RefVertex* refvert1 = subseg->get_beg_ref_vert();

  if (refvert1) {
    Vert* vertex = subseg->get_beg_vert();
    pair<ItMap, ItMap> it_pair = vert_map.equal_range(vertex);
    for (ItMap it = it_pair.first; it != it_pair.second; it++) {
      if (it->second == subseg) {
	vert_map.erase(it);
	break;
      }
    }
  }

  RefVertex* refvert2 = subseg->get_end_ref_vert();

  if (refvert2) {
    Vert* vertex = subseg->get_end_vert();
    pair<ItMap, ItMap> it_pair = vert_map.equal_range(vertex);
    for (ItMap it = it_pair.first; it != it_pair.second; it++) {
      if (it->second == subseg) {
	vert_map.erase(it);
	break;
      }
    }
  }

}

void
SmallAngleTool::get_subsegs_of_vert(Vert* const vertex,
				    vector<Subseg*>& subsegments) const
{

  pair<ItcMap, ItcMap> it_pair = vert_map.equal_range(vertex);

  for (ItcMap it = it_pair.first; it != it_pair.second; it++)
    subsegments.push_back(it->second);

}

void
SmallAngleTool::get_all_subsegs(std::vector<Subseg*>& subsegments) const
{

  for (ItcMap it = vert_map.begin(); it != vert_map.end(); it++)
    subsegments.push_back(it->second);

}

double
SmallAngleTool::compute_split_radius(
    const vector<Subseg*>& subsegs_to_split) const
{

  double split_radius = LARGE_DBL;

  for (unsigned int i = 0; i < subsegs_to_split.size(); i++) {
    Subseg* subsegment = subsegs_to_split[i];
    double radius = 0.5 * subsegment->get_length();
    if (radius < split_radius)
      split_radius = radius;
  }

  return split_radius;

}

bool
SmallAngleTool::vertex_has_small_angle(Vert* const vertex) const
{

  if (vert_map.count(vertex) == 0)
    return false;
  else
    return true;

}

bool
SmallAngleTool::subseg_has_small_angle(const Subseg* const subseg) const
{

  Vert* vert = subseg->get_beg_vert();
  pair<ItcMap, ItcMap> it_pair = vert_map.equal_range(vert);

  for (ItcMap it = it_pair.first; it != it_pair.second; it++)
    if (it->second == subseg)
      return true;

  vert = subseg->get_end_vert();
  it_pair = vert_map.equal_range(vert);

  for (ItcMap it = it_pair.first; it != it_pair.second; it++)
    if (it->second == subseg)
      return true;

  return false;

}

void
SmallAngleTool::print_map_to_screen() const
{

  printf("--- Vertex --- --- Subseg ---\n");

  for (ItcMap it = vert_map.begin(); it != vert_map.end(); it++) {
    printf("   %p      %p %d\n", it->first, it->second,
	   it->second->is_deleted());
  }

}
