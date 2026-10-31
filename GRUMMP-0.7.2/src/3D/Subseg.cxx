#include "GR_Subseg.h"
#include "GR_misc.h"
#include "GR_Vertex.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include <algorithm>

SubsegBridge::~SubsegBridge()
{
  m_parent = static_cast<Subseg*>(NULL);
}

Subseg::~Subseg()
{

  if (m_vd1)
    delete m_vd1;
  if (m_vd2)
    delete m_vd2;

  for (std::vector<SubsegBridge*>::iterator it = m_children.begin();
      it != m_children.end(); ++it) {
    if ((*it))
      delete (*it);
  }
}

double
Subseg::get_length() const
{

  double vector[] =
    { m_vd2->vertex->x() - m_vd1->vertex->x(), m_vd2->vertex->y()
	- m_vd1->vertex->y(), m_vd2->vertex->z() - m_vd1->vertex->z() };
  return dMAG3D(vector);

}

double
Subseg::param_at_distance_from_vert(const double distance,
				    const Vert* const vert) const
{

  assert(distance >= 0.);
  assert(distance <= get_length());
  assert(vert == m_vd1->vertex || vert == m_vd2->vertex);
  assert(m_vd1->param < m_vd2->param);
  assert(m_parent_edge);

  bool from_lo = vert == m_vd1->vertex;
  CubitVector from_vert_coord(vert->x(), vert->y(), vert->z());

  double param_lo = m_vd1->param, param_up = m_vd2->param;
  CubitVector coord_lo, coord_up, current_coord;
  m_parent_edge->position_from_u(param_lo, coord_lo);
  m_parent_edge->position_from_u(param_up, coord_up);
  double dist_lo = (coord_lo - from_vert_coord).length();
  double dist_up = (coord_up - from_vert_coord).length();

  assert(
      (dist_lo <= distance && dist_up >= distance)
	  || (dist_lo >= distance && dist_up <= distance));

  double current_param = -LARGE_DBL, current_dist = LARGE_DBL;

  // FIX ME!
  //Having a fixed tolerance like this is likely to cause problems.
  while (fabs(current_dist - distance) > 1.e-12) {

    current_param = 0.5 * (param_lo + param_up);
    m_parent_edge->position_from_u(current_param, current_coord);
    current_dist = (current_coord - from_vert_coord).length();

    if (from_lo) {
      if (current_dist > distance)
	param_up = current_param;
      else
	param_lo = current_param;
    }
    else {
      if (current_dist > distance)
	param_lo = current_param;
      else
	param_up = current_param;
    }

    //     printf("current_param = %lf, error = %.5e\n",
    // 	   current_param, fabs(current_dist - distance));

  }

  //   printf("current param = %lf - m_vd1->param = %lf - m_vd2->param = %lf - distance = %e - current_dist = %e\n",
  // 	 current_param, m_vd1->param, m_vd2->param, distance, current_dist);

  return current_param;

}

RefVertex*
Subseg::get_beg_ref_vert() const
{

  return dynamic_cast<RefVertex*>(get_beg_vert()->getParentTopology());

}

RefVertex*
Subseg::get_end_ref_vert() const
{

  return dynamic_cast<RefVertex*>(get_end_vert()->getParentTopology());

}

void
Subseg::compute_refinement_split_data(double split_coords[3],
				      double& split_param, bool qForceCoord)
{

  assert(m_vd1 && m_vd2);

  Vert* vert1 = m_vd1->vertex, *vert2 = m_vd2->vertex;
  assert(vert1 && !vert1->isDeleted() && vert1->getSpaceDimen() == 3);
  assert(vert2 && !vert2->isDeleted() && vert2->getSpaceDimen() == 3);

  assert(
      vert1->getVertType() == Vert::eBdryApex
	  || vert1->getVertType() == Vert::eBdryCurve);
  assert(
      vert2->getVertType() == Vert::eBdryApex
	  || vert2->getVertType() == Vert::eBdryCurve);
  if (qForceCoord) {
    if (m_parent_edge) {
      CubitVector temp(split_coords);
      split_param = m_parent_edge->u_from_position(temp);
    }
    else {
      double dLength = get_length();
      double dist = dDIST3D(split_coords, vert1->getCoords());
      double tempParam = dist / dLength;
      split_param = m_vd1->param + tempParam * (m_vd2->param - m_vd1->param);
    }
    return;
  }
  //This is the regular split case. The split point will
  //be at the mid parameter (or at the mid point if there is no parent curve)
  if (vert1->getVertType() == vert2->getVertType()) {

    //No parent curve defined, split at mid point.
    if (!m_parent_edge) {

      //Compute the split coordinates
      split_coords[0] = 0.5 * (vert1->x() + vert2->x());
      split_coords[1] = 0.5 * (vert1->y() + vert2->y());
      split_coords[2] = 0.5 * (vert1->z() + vert2->z());

      //Set the parameter manually to some dummy value.
      split_param = -LARGE_DBL;

    }

    //Split at middle parameter between
    else {

      //Compute the split parameter.
      assert(
	  m_vd1->param >= m_parent_edge->start_param()
	      && m_vd1->param <= m_parent_edge->end_param());
      assert(
	  m_vd2->param >= m_parent_edge->start_param()
	      && m_vd2->param <= m_parent_edge->end_param());

      split_param = 0.5 * (m_vd1->param + m_vd2->param);

      assert(split_param > m_vd1->param && split_param < m_vd2->param);
      assert(
	  split_param >= m_parent_edge->start_param()
	      && split_param <= m_parent_edge->end_param());

      //Compute the split coordinates
      CubitVector tmp_coord;
      m_parent_edge->position_from_u(split_param, tmp_coord);
      tmp_coord.get_xyz(split_coords);

    }

  }

  //In this case, exactly one of the verts is a BdryApex.
  //To make sure all edges with small incident surface angles get
  //split safely, split all such segments in power of 2 lengths,
  //using the biggest power of 2 less than 2/3 of current edge
  //length.  This will ensure a split between 1/3 and 2/3.
  else {

    assert(
	(vert1->getVertType() == Vert::eBdryApex
	    && vert2->getVertType() == Vert::eBdryCurve)
	    || (vert2->getVertType() == Vert::eBdryApex
		&& vert1->getVertType() == Vert::eBdryCurve));

    //Compute the log_e of 2/3 the length, then convert to log_2 by
    //multiplying by log_2(e) (equivalent to dividing by
    //log_e(2)).  Then take the largest integer below this.

    double length = get_length();
    double split_power = floor(log(length / 1.5) * M_LOG2E);
    double split_dist = pow(2., split_power);
    assert(length <= 3. * split_dist && 1.5 * split_dist <= length);

    //No parent curve defined, split at split_dist from apex vertex.
    if (!m_parent_edge) {

      double subseg_vec[] =
      adDIFF3D(vert2->getCoords(), vert1->getCoords());
      NORMALIZE3D(subseg_vec);

      if (vert1->getVertType() == Vert::eBdryApex) {
	split_coords[0] = vert1->x() + split_dist * subseg_vec[0];
	split_coords[1] = vert1->y() + split_dist * subseg_vec[1];
	split_coords[2] = vert1->z() + split_dist * subseg_vec[2];
      }
      else {
	split_coords[0] = vert2->x() - split_dist * subseg_vec[0];
	split_coords[1] = vert2->y() - split_dist * subseg_vec[1];
	split_coords[2] = vert2->z() - split_dist * subseg_vec[2];
      }

      //Set the parameter manually to some dummy value.
      split_param = -LARGE_DBL;

    }

    //We have a parent curve attached, so find the parameter on the curve
    //that is at split_dist from the apex vertex.
    else {

      //Compute the split parameter.
      assert(
	  m_vd1->param >= m_parent_edge->start_param()
	      && m_vd1->param <= m_parent_edge->end_param());
      assert(
	  m_vd2->param >= m_parent_edge->start_param()
	      && m_vd2->param <= m_parent_edge->end_param());

      split_param =
	  vert1->getVertType() == Vert::eBdryApex ?
	      param_at_distance_from_vert(split_dist, vert1) :
	      param_at_distance_from_vert(split_dist, vert2);

      assert(split_param > m_vd1->param && split_param < m_vd2->param);
      assert(
	  split_param >= m_parent_edge->start_param()
	      && split_param <= m_parent_edge->end_param());

      //Compute the split coordinates
      CubitVector tmp_coord;
      m_parent_edge->position_from_u(split_param, tmp_coord);
      tmp_coord.get_xyz(split_coords);

    }

  }

}

void
Subseg::refinement_split(const double split_coords[3], const double split_param,
			 Vert* const new_vert, Subseg* const new_subseg1,
			 Subseg* const new_subseg2)
{

  //Set the new subsegments param.
  new_subseg1->m_vd1->param = m_vd1->param;
  new_subseg1->m_vd2->param = new_subseg2->m_vd1->param = split_param;
  new_subseg2->m_vd2->param = m_vd2->param;

  //Perform the split
  split(split_coords, new_vert, new_subseg1, new_subseg2);

#ifndef NDEBUG

  if (new_subseg1->m_parent_edge) {
    assert(
	new_subseg1->m_vd1->param >= m_parent_edge->start_param()
	    && new_subseg1->m_vd1->param <= m_parent_edge->end_param());
    assert(
	new_subseg1->m_vd2->param >= m_parent_edge->start_param()
	    && new_subseg1->m_vd2->param <= m_parent_edge->end_param());
  }
  if (new_subseg2->m_parent_edge) {
    assert(
	new_subseg2->m_vd1->param >= m_parent_edge->start_param()
	    && new_subseg2->m_vd1->param <= m_parent_edge->end_param());
    assert(
	new_subseg2->m_vd2->param >= m_parent_edge->start_param()
	    && new_subseg2->m_vd2->param <= m_parent_edge->end_param());
  }

#endif

}

void
Subseg::refinement_split(Vert* const new_vert, Subseg* const new_subseg1,
			 Subseg* const new_subseg2)
{

  double split_coords[3], split_param;
  compute_refinement_split_data(split_coords, split_param);
  refinement_split(split_coords, split_param, new_vert, new_subseg1,
		   new_subseg2);

}

void
Subseg::split_at_param(const double split_param, Subseg* const new_subseg1,
		       Subseg* const new_subseg2, Vert* const new_vert)
{

  //	assert( m_parent_edge );
  assert(is_not_deleted());

  //Make sure this param is within a valid range.
  assert(split_param > m_vd1->param && split_param < m_vd2->param);
  double split_coords[3];
  // allows splitting of subsegs without background geometry - assumes straight lines
  if (m_parent_edge) {
    CubitVector tmp_coords;

    m_parent_edge->position_from_u(split_param, tmp_coords);
    //		printf("Dist %.16f %.16f %.16f\n",(tmp_coords-CubitVector(get_beg_vert()->getCoords())).length(),
    //				(tmp_coords-CubitVector(get_end_vert()->getCoords())).length(),
    //				(CubitVector(get_beg_vert()->getCoords())-CubitVector(get_end_vert()->getCoords())).length());
    tmp_coords.get_xyz(split_coords);

  }
  else {
    double split_fraction = 0.5;    //split_param/(m_vd2->param - m_vd1->param);
    split_coords[0] = m_vd1->vertex->x()
	+ split_fraction * (m_vd2->vertex->x() - m_vd1->vertex->x());
    split_coords[1] = m_vd1->vertex->y()
	+ split_fraction * (m_vd2->vertex->y() - m_vd1->vertex->y());
    split_coords[2] = m_vd1->vertex->z()
	+ split_fraction * (m_vd2->vertex->z() - m_vd1->vertex->z());
  }
  //Obtain the split coordinates from the param.

  //Perform the split with the info we gathered.
  split(split_coords, new_vert, new_subseg1, new_subseg2);

  //We know the split param, so set it in the new subsegments.
  new_subseg1->m_vd1->param = m_vd1->param;
  new_subseg1->m_vd2->param = split_param;
  new_subseg2->m_vd1->param = split_param;
  new_subseg2->m_vd2->param = m_vd2->param;

}

void
Subseg::split_at_coords(const CubitVector& split_coords, Vert* const new_vert,
			Subseg* const new_subseg1, Subseg* const new_subseg2)
{

  double split_loc[3];
  split_coords.get_xyz(split_loc);

  split(split_loc, new_vert, new_subseg1, new_subseg2);

}

void
Subseg::split(const double split_coords[3], Vert* const new_vert,
	      Subseg* const new_subseg1, Subseg* const new_subseg2)
{

  assert(new_vert && !new_vert->isDeleted());
  assert(new_subseg1 && new_subseg1->is_not_deleted());
  assert(new_subseg2 && new_subseg2->is_not_deleted());

  mark_deleted();

  new_vert->setCoords(3, split_coords);
  new_vert->setType(Vert::eBdryCurve);

  new_vert->setParentTopology(m_parent_edge);

  assert(new_subseg1->m_vd2);
  assert(new_subseg2->m_vd1);
  assert(new_subseg2->m_vd2);

  new_subseg1->m_vd1->vertex = m_vd1->vertex;
  new_subseg1->m_vd2->vertex = new_vert;
  new_subseg2->m_vd1->vertex = new_vert;
  new_subseg2->m_vd2->vertex = m_vd2->vertex;

  new_subseg1->m_parent_edge = m_parent_edge;
  new_subseg2->m_parent_edge = m_parent_edge;

  new_subseg1->m_deleted = false;
  new_subseg2->m_deleted = false;

}

CubitVector
Subseg::mid_point() const
{

  return CubitVector(0.5 * (m_vd1->vertex->x() + m_vd2->vertex->x()),
		     0.5 * (m_vd1->vertex->y() + m_vd2->vertex->y()),
		     0.5 * (m_vd1->vertex->z() + m_vd2->vertex->z()));

}

double
Subseg::ball_radius() const
{

  return 0.5 * calcDistanceBetween(m_vd1->vertex, m_vd2->vertex);

}

CubitBox
Subseg::bounding_box() const
{

  CubitVector midpoint = mid_point();
  double radius = ball_radius();

  double x1 = midpoint.x() - radius, x2 = midpoint.x() + radius;
  double y1 = midpoint.y() - radius, y2 = midpoint.y() + radius;
  double z1 = midpoint.z() - radius, z2 = midpoint.z() + radius;

  CubitVector max_coord(std::max(x1, x2), std::max(y1, y2), std::max(z1, z2));
  CubitVector min_coord(std::min(x1, x2), std::min(y1, y2), std::min(z1, z2));

  return CubitBox(min_coord, max_coord);

}

void
Subseg::closest_on_subseg(const double pt[3], double close[3]) const
{

  Vert *vert0 = get_beg_vert(), *vert1 = get_end_vert();

  double diff0[] = adDIFF3D(vert0->getCoords(), pt);
  double diff1[] = adDIFF3D(vert1->getCoords(), vert0->getCoords());
  double param = -dDOT3D(diff0, diff1) / dMAG3D_SQ(diff1);

  if (param <= 0.) {
    close[0] = vert0->x();
    close[1] = vert0->y();
    close[2] = vert0->z();
    return; //no need to project.
  }

  else if (param >= 1.) {
    close[0] = vert1->x();
    close[1] = vert1->y();
    close[2] = vert1->z();
    return; //no need to project.
  }

  else {

    close[0] = vert0->x() + param * diff1[0];
    close[1] = vert0->y() + param * diff1[1];
    close[2] = vert0->z() + param * diff1[2];

    if (m_parent_edge) { //project onto curve
      CubitVector close_pt(close), closest;
      m_parent_edge->closest_point(close_pt, closest);
      closest.get_xyz(close);
    }

  }

}

SubsegMap::~SubsegMap()
{

  //If the m_to_release vector is not empty, start by releasing its content.
  release_deleted_subsegs();

  //The subsegments left in the map must be deleted.
  //They were instantiated elsewhere. Subseg entries appear
  //in more than one key, we need to filter these first.

  SubsegSet to_delete;
  SubsegMapType::iterator it = m_map.begin(), it_end = m_map.end();

  for (; it != it_end; ++it)
    std::copy(it->second.begin(), it->second.end(),
	      std::inserter(to_delete, to_delete.begin()));

  //All remaining subsegs are in the to_delete set. Let's delete.
  for (SubsegSet::iterator it2 = to_delete.begin(); it2 != to_delete.end();
      ++it2) {
    if ((*it2))
      delete (*it2);
  }

  m_map.clear();

}

void
SubsegMap::release_deleted_subsegs()
{
  std::vector<Subseg*>::iterator it = m_to_release.begin(), it_end =
      m_to_release.end();

#ifndef NDEBUG
  for (; it != it_end; ++it)
    assert((*it)->is_deleted());
#endif

  for (it = m_to_release.begin(); it != it_end; ++it)
    if ((*it))
      delete (*it);

  m_to_release.clear();

}

void
SubsegMap::add_subseg(Subseg* const subseg)
{

  assert(subseg);
  Vert *vert0 = subseg->get_beg_vert(), *vert1 = subseg->get_end_vert();

#ifndef NDEBUG

  assert(vert0 && vert1);
  assert(vert0 != vert1);
  assert(!vert0->isDeleted());
  assert(!vert1->isDeleted());
  //	assert( vert0->getVertType() == Vert::eBdryApex ||
  //			vert0->getVertType() == Vert::eBdryCurve );
  //	assert( vert1->getVertType() == Vert::eBdryApex ||
  //			vert1->getVertType() == Vert::eBdryCurve );

  //	std::set<Cell*> neigh_cells;
  //	findNeighborhoodInfo(vert0, neigh_cells, neigh_verts);
  std::set<Vert*> neigh_verts;
  for (int iF = 0; iF < vert0->getNumFaces(); iF++) {
    Face * pF = vert0->getFace(iF);
    for (int iV = 0; iV < pF->getNumVerts(); iV++) {
      neigh_verts.insert(pF->getVert(iV));
    }
  }
  assert(neigh_verts.count(vert1) == 1);

#endif

  SubsegMapType::iterator it0 =
      m_map.insert(std::make_pair(vert0, SubsegSet())).first, it1 =
      m_map.insert(std::make_pair(vert1, SubsegSet())).first;

  it0->second.insert(subseg);
  it1->second.insert(subseg);

  assert(it0->first == vert0);
  assert(it1->first == vert1);
  assert(it0->second.count(subseg) == 1);
  assert(it1->second.count(subseg) == 1);

}

void
SubsegMap::remove_subseg(Subseg* const subseg)
{

  assert(subseg);
  Vert *vert[] =
    { subseg->get_beg_vert(), subseg->get_end_vert() };

  assert(vert[0] && vert[1]);
  assert(vert[0] != vert[1]);
  assert(vert[0]->isBdryVert());
  assert(vert[1]->isBdryVert());

  //Erase the two entries:
  SubsegMapType::iterator it;

  for (int i = 0; i < 2; ++i) {

    it = m_map.find(vert[i]);
    assert(it != m_map.end());             //vertex must be in map.
    assert(it->second.count(subseg) == 1); //subseg must be connected to vert.

    it->second.erase(subseg);
    if (it->second.empty())
      m_map.erase(it);

  }

  //Finally, move the subseg into the m_to_release vector. We do
  //not want to release the memory right now as the Subseg might still
  //be stored in the TetMeshRefiner::m_subseg_queue.
  assert(subseg->is_deleted());
  m_to_release.push_back(subseg);

  //   for(vector<Subseg*>::iterator its = m_to_release.begin(); its != m_to_release.end(); ++its) {
  //     Subseg* my_sub = *its;
  //     printf("Subseg %p, deleted = %d\n", my_sub, my_sub->is_deleted());
  //   }

}

bool
SubsegMap::vert_has_subseg(Vert* const vertex, Subseg* const subseg) const
{

  assert(vertex);
  assert(!vertex->isDeleted());
  assert(
      vertex->getVertType() == Vert::eBdryApex
	  || vertex->getVertType() == Vert::eBdryCurve);

  SubsegMapType::const_iterator it = m_map.find(vertex);
  assert(it != m_map.end());

  return (it->second.count(subseg) == 1);

}

Subseg*
SubsegMap::verts_are_connected(Vert* const vert0, Vert* const vert1) const
{

  assert(vert0 && vert1);
  assert(vert0 != vert1);
  assert(!vert0->isDeleted());
  assert(!vert1->isDeleted());

  //If either of the vertices is not at an apex or on
  //a curve, then it definitely cannot be part of a subsegment.
  if (vert0->getVertType() != Vert::eBdryApex
      && vert0->getVertType() != Vert::eBdryCurve) {
    assert(m_map.find(vert0) == m_map.end());
    return static_cast<Subseg*>(NULL);
  }
  else if (vert1->getVertType() != Vert::eBdryApex
      && vert1->getVertType() != Vert::eBdryCurve) {
    assert(m_map.find(vert1) == m_map.end());
    return static_cast<Subseg*>(NULL);
  }

  //If we get to this point, then both vertices must be at
  //an apex or on a curve.
  assert(
      (vert0->getVertType() == Vert::eBdryApex
	  || vert0->getVertType() == Vert::eBdryCurve)
	  && (vert1->getVertType() == Vert::eBdryApex
	      || vert1->getVertType() == Vert::eBdryCurve));

  SubsegMapType::const_iterator it = m_map.find(vert0);

  //If vert0 is not in the map, then for sure there are
  //no subsegment connecting vert0 to vert1;
  if (it == m_map.end()) {
#ifndef NDEBUG
    it = m_map.find(vert1);
    if (it != m_map.end())
      assert(
	  std::find_if(it->second.begin(), it->second.end(),
		       SubsegHasVert(vert0)) == it->second.end());
#endif
    return static_cast<Subseg*>(NULL);
  }

  SubsegSet::const_iterator its = std::find_if(it->second.begin(),
					       it->second.end(),
					       SubsegHasVert(vert1));

  Subseg* ret_subseg =
      its != it->second.end() ? *its : static_cast<Subseg*>(NULL);

#ifndef NDEBUG

  if (ret_subseg) {

    assert(ret_subseg->has_vert(vert0));
    assert(ret_subseg->has_vert(vert1));
    //Make sure these two verts are connected in the mesh.
    int iMax0 = vert0->getNumFaces();

    bool foundIt = false;
    for (int i0 = 0; i0 < iMax0; i0++) {
      Face *pFCand = vert0->getFace(i0);
      assert(pFCand->isValid());
      assert(!pFCand->isDeleted());
      if (pFCand->hasVert(vert1)) {
	foundIt = true;
	break;
      }
    }
    assert(foundIt);
//		std::set<Cell*> neigh_cells;
//		std::set<Vert*> neigh_verts;
//		findNeighborhoodInfo(vert0, neigh_cells, neigh_verts);
//		assert(neigh_verts.count(vert1) == 1);

  }

#endif

  return ret_subseg;

}

const std::set<Subseg*>&
SubsegMap::get_subsegs(Vert* const vert0) const
{

  assert(vert0);
  assert(!vert0->isDeleted());
  assert(
      vert0->getVertType() == Vert::eBdryApex
	  || vert0->getVertType() == Vert::eBdryCurve);

  assert(m_map.count(vert0) == 1);
  return m_map.find(vert0)->second;

}
bool
SubsegMap::isSubsegVert(Vert * vert) const
{
  if (vert->isDeleted()
      || !(vert->getVertType() == Vert::eBdryApex
	  || vert->getVertType() == Vert::eBdryCurve))
    return false;
  return m_map.count(vert) > 0;
}
std::pair<SubsegMap::SubsegSet::const_iterator,
    SubsegMap::SubsegSet::const_iterator>
SubsegMap::get_subseg_range(Vert* const vert0) const
{

  assert(vert0);
  assert(isSubsegVert(vert0));

  SubsegMapType::const_iterator it = m_map.find(vert0);
  assert(it != m_map.end());
  return std::make_pair(it->second.begin(), it->second.end());

}

void
SubsegMap::get_all_subsegs(std::set<Subseg*>& all_subsegs) const
{

  all_subsegs.clear();

  SubsegMapType::const_iterator it = m_map.begin(), it_end = m_map.end();

  for (; it != it_end; ++it)
    all_subsegs.insert(it->second.begin(), it->second.end());

}
GR_index_t
SubsegMap::num_subsegs()
{

  SubsegMapType::const_iterator it = m_map.begin(), it_end = m_map.end();
  GR_index_t num = 0;
  for (; it != it_end; ++it)
    num += it->second.size();
  return num;
}

void
SubsegMap::get_connected_verts(Vert* const vert0,
			       std::vector<Vert*>& verts) const
{

  assert(vert0);
  assert(!vert0->isDeleted());
  assert(
      vert0->getVertType() == Vert::eBdryApex
	  || vert0->getVertType() == Vert::eBdryCurve);

  assert(m_map.count(vert0) == 1);

  Subseg* subseg;
  const SubsegSet subseg_set = get_subsegs(vert0);
  SubsegSet::const_iterator it = subseg_set.begin(), it_end = subseg_set.end();

  for (; it != it_end; ++it) {

    subseg = *it;
    assert(subseg->has_vert(vert0));

    verts.push_back(
	vert0 == subseg->get_beg_vert() ?
	    subseg->get_end_vert() : subseg->get_beg_vert());

  }
}
void
SubsegMap::replaceVert(Vert * oldVert, Vert * newVert)
{

  // first, for every subseg on this vertex
  SubsegMapType::iterator itm = m_map.find(oldVert);

  for (std::set<Subseg*>::const_iterator its = itm->second.begin();
      its != itm->second.end(); ++its) {
    Subseg * subseg = *its;
    subseg->replace_vert(oldVert, newVert);
  }
  m_map.insert(std::make_pair(newVert, itm->second));
  m_map.erase(itm);
}
//void SubsegMap::update_vert_data_after_purge(const std::map<Vert*, Vert*>& vert_map) {
//
//	//Start by releasing the deleted subsegments, if any.
//	release_deleted_subsegs();
//
//	//Since all (or most of) the keys in the map changed,
//	//we might as well rebuild the entire map.
//	Subseg* subseg;
//	Vert* old_vert, *new_vert;
//	SubsegSet::iterator its, its_end;
//	SubsegMapType updated_map;
//	SubsegMapType::iterator itm = m_map.begin(), itm_end = m_map.end();
//
//	for( ; itm != itm_end; ++itm) {
//
//		old_vert = itm->first;
//		assert(vert_map.count(old_vert) == 1);
//		new_vert = vert_map.find(old_vert)->second;
//
//		its = itm->second.begin();
//		its_end = itm->second.end();
//
//		for( ; its != its_end; ++its) {
//			subseg = *its;
//			assert(subseg->is_not_deleted());
//			subseg->replace_vert(old_vert, new_vert);
//		}
//
//		updated_map.insert( std::make_pair(new_vert, itm->second) );
//
//	}
//
//	//Swap the updated map with the old map and check validity.
//	updated_map.swap(m_map);
//	assert(valid());
//
//}

//Checks the validity of the map entries.
bool
SubsegMap::valid() const
{

  Vert *vert0, *vert1;
  Subseg* subseg;
  SubsegMapType::const_iterator it = m_map.begin(), it_end = m_map.end(),
      it_tmp;
  SubsegSet::const_iterator its, its_end;

  std::set<Cell*> neigh_cells;
  std::set<Vert*> neigh_verts;

  //Go through all vertices stored in the map.

  for (; it != it_end; ++it) {

    vert0 = it->first;

    switch (vert0->getVertType())
      {
      case Vert::eBdryApex:
	break;
      case Vert::eBdryCurve:
	if (it->second.size() == 1) {
	  // Termination of a string of subsegs; relabel this as an apex
	  vert0->setType(Vert::eBdryApex);
	}
//			else if (it->second.size() != 2)
//				return false;
	break;
      default:
	return false;
      }
    findNeighborhoodInfo(vert0, neigh_cells, neigh_verts);

    ///Check all subsegs connected to vert0 for validity.

    its = it->second.begin();
    its_end = it->second.end();

    for (; its != its_end; ++its) {

      subseg = *its;
      assert(subseg->get_beg_vert());
      assert(subseg->get_end_vert());

      if (!subseg->has_vert(vert0))
	return false;

      //Obviously, the two vertices stored in the subseg
      //must be different.
      vert1 =
	  vert0 == subseg->get_beg_vert() ?
	      subseg->get_end_vert() : subseg->get_beg_vert();
      if (vert0 == vert1)
	return false;

      //The same subseg must have an entry in the map both
      //for vert0 and vert1.
      it_tmp = m_map.find(vert1);
      if (it_tmp == m_map.end())
	return false;
      if (it_tmp->second.count(subseg) != 1)
	return false;

      //If there is a subseg, then the two verts must
      //also be connected in the mesh.
      if (neigh_verts.count(vert1) != 1)
	return false;

    }

  }

  return true;

}

void
SubsegMap::print() const
{

  Subseg* subseg;
  SubsegMapType::const_iterator it = m_map.begin(), it_end = m_map.end();
  SubsegSet::const_iterator its, its_end;
  GR_index_t nSubsegs = 0;
  Vert * pV;
  for (; it != it_end; ++it) {
    its = it->second.begin();
    its_end = it->second.end();
    for (; its != its_end; ++its) {
      nSubsegs++;
    }
  }
  FILE *pFOutFile = NULL;

  pFOutFile = fopen("SubsegMap.vtk", "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  GR_index_t nno = nSubsegs * 2;   // number of nodes
  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nno);
  it = m_map.begin(), it_end = m_map.end();
  for (; it != it_end; ++it) {
    its = it->second.begin();
    its_end = it->second.end();
    for (; its != its_end; ++its) {
      subseg = *its;
      pV = subseg->get_beg_vert();
      fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
      pV = subseg->get_end_vert();
      fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());

    }
  }
  int counter = -1;
  fprintf(pFOutFile, "CELLS %d %d\n", nSubsegs, nSubsegs * 3);
  it = m_map.begin(), it_end = m_map.end();

  for (; it != it_end; ++it) {
    its = it->second.begin();
    its_end = it->second.end();
    for (; its != its_end; ++its) {
      fprintf(pFOutFile, "2 %d %d\n", counter + 1, counter + 2);
      counter = counter + 2;

    }
  }

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "CELL_TYPES %d\n", nSubsegs);
  it = m_map.begin(), it_end = m_map.end();
  for (; it != it_end; ++it) {
    its = it->second.begin();
    its_end = it->second.end();
    for (; its != its_end; ++its) {
      fprintf(pFOutFile, "3\n");
    }
  }
  fclose(pFOutFile);
}
