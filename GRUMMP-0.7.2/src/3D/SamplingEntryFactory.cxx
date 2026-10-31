#include "GR_SamplingEntryFactory.h"
#include "GR_Face.h"
#include "GR_SamplingQueue.h"
#include "GR_SubsegManager.h"
#include "GR_SurfaceSampler.h"
#include "GR_SurfIntersect.h"
#include "GR_Vertex.h"

#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

using std::vector;
using std::set;

SamplingEntryFactory*
SamplingEntryFactory::instance()
{
  static SamplingEntryFactory factory;
  return &factory;
}
SamplingEntryFactory::~SamplingEntryFactory()
{

//	for(FaceQueueEntrySet::iterator it = m_used_faces.begin(); it != m_used_faces.end(); ++it){
//		if((*it)) delete (*it);
//	}
}
void
SamplingEntryFactory::build_entries(
    const SurfaceSampler* const sampler,
    const vector<VertWithTag>& sample_verts,
    const vector<SurfaceSampler::RestrictFace>& restrict_faces,
    vector<SamplingEntry*>* new_entries)
{

  if (sample_verts.empty() && restrict_faces.empty())
    return;

  assert(new_entries);
  reset_factory();

  m_sampler = sampler;
  m_new_entries = new_entries;
  assert(valid_data_passed(sample_verts));
  assert(valid_data_passed(restrict_faces));

  for (vector<SurfaceSampler::RestrictFace>::const_iterator it =
      restrict_faces.begin(); it != restrict_faces.end(); ++it)
    initialize_face(*it);

  for (vector<VertWithTag>::const_iterator it = sample_verts.begin();
      it != sample_verts.end(); ++it)
    inspect_vertex(*it);

  for (FaceQueueEntrySet::iterator it = m_test_faces.begin();
      it != m_test_faces.end(); ++it)
    inspect_face(*it);

}

void
SamplingEntryFactory::build_entries(
    const SurfaceSampler* const sampler, const VertWithTag& vertex_and_tag,
    const vector<SurfaceSampler::RestrictFace>& restrict_faces,
    vector<SamplingEntry*>* new_entries)
{

  if (restrict_faces.empty())
    return;

  assert(new_entries);
  reset_factory();

  m_sampler = sampler;
  m_new_entries = new_entries;
  assert(valid_data_passed(vertex_and_tag));
//  assert( valid_data_passed(restrict_faces) );

  for (vector<SurfaceSampler::RestrictFace>::const_iterator it =
      restrict_faces.begin(); it != restrict_faces.end(); ++it)
    initialize_face(*it);

  //assert( test_faces_have_vertex(vertex_and_tag.first) );

  if (inspect_vertex(vertex_and_tag)) {
    //Not sure if the following is really necessary.
//     if( !check_peripheral_edges(vertex_and_tag) )
//       build_topo_violation(vertex_and_tag.first, vertex_and_tag.second);
  }

  for (FaceQueueEntrySet::iterator it = m_test_faces.begin();
      it != m_test_faces.end(); ++it)
    inspect_face(*it);

}

void
SamplingEntryFactory::build_entries(const SurfaceSampler* const sampler,
				    const VertWithTag& vertex_and_tag,
				    vector<SamplingEntry*>* new_entries)
{

  assert(new_entries);
  reset_factory();

  m_sampler = sampler;
  m_new_entries = new_entries;
  assert(valid_data_passed(vertex_and_tag));

  inspect_vertex(vertex_and_tag);

}

void
SamplingEntryFactory::reset_factory()
{

  m_sampler = static_cast<SurfaceSampler*>(NULL);

  m_test_faces.clear();
  m_used_faces.clear();
  m_checked_edges.clear();
  m_checked_bridges.clear();
  m_bad_bridges.clear();

}

void
SamplingEntryFactory::initialize_face(
    const SurfaceSampler::RestrictFace& restrict_face)
{

  FaceQueueEntry* fqe_to_test = restrict_face.first;

  //This case might happen when 2 subsegs are split at the same time
  //due to the presence of a small angle. The splits are done first,
  //and the sampling entries are built. It is possible to two splits
  //on a same face cause either a face to get deleted, or a face to
  //loose it restricted Delaunay status. Therefore, we need to perform
  //the following check.

  if (fqe_to_test->isDeleted())
    return;

  if (!fqe_to_test->getFace()->isRestrictedDelaunay())
    return;

  m_test_faces.insert(fqe_to_test);

  switch (restrict_face.second.size())
    {
    case 1: //usual case, does nothing.
      break;
    default: //multiple surf intersection, build entry.
      build_multiple_intersect(fqe_to_test);
      break;
    }

}

bool
SamplingEntryFactory::inspect_vertex(const VertWithTag& vertex_and_tag)
{

  Vert* vertex = vertex_and_tag.first;
  SampleVertTag* vertex_tag = vertex_and_tag.second;
  switch (vertex->getVertType())
    {

    case Vert::eBdryApex:
      return inspect_apex_vert(vertex, vertex_tag);
      break;
    case Vert::eBdryCurve:
    case Vert::eBdryTwoSide:
      return inspect_curv_vert(vertex, vertex_tag);
      break;
    case Vert::ePseudoSurface:
      return inspect_surf_vert(vertex, vertex_tag);
      break;
    default:
      vertex->printVertInfo();
      vFatalError("", "SamplingEntryFactory::inspect_vertex");
      break;
    }
  return false;
}

void
SamplingEntryFactory::inspect_face(FaceQueueEntry* const fqe)
{

  if (m_used_faces.find(fqe) != m_used_faces.end())
    return;

  if (!inspect_face_quality(fqe))
    return;

}

bool
SamplingEntryFactory::inspect_face_quality(FaceQueueEntry* const fqe)
{
  Face * face = fqe->getFace();
  if (!m_qual_samp)
    return true; //no quality sampling wanted.

//  assert( m_used_faces.find(face) == m_used_faces.end() );

  //If the triangle is at a small angle,
  //we do not want to touch it.
  if (face->getVert(0)->isSmallAngleVert()
      || face->getVert(1)->isSmallAngleVert()
      || face->getVert(2)->isSmallAngleVert()) {
    return true;
  }

  double shortest =
      MIN3(
	  dDIST_SQ_3D(face->getVert(1)->getCoords(), face->getVert(0)->getCoords()),
	  dDIST_SQ_3D(face->getVert(2)->getCoords(), face->getVert(0)->getCoords()),
	  dDIST_SQ_3D(face->getVert(2)->getCoords(), face->getVert(1)->getCoords()));
  TriFace* tface = dynamic_cast<TriFace*>(face);
  assert(tface);
  double radius_edge_ratio = tface->calcCircumradius() / sqrt(shortest);

  if (radius_edge_ratio > m_radius_edge_bound) {
    build_quality(fqe, radius_edge_ratio);
    return false; //need to improve face quality.
  }

  return true; //face passes quality requirement.

}

bool
SamplingEntryFactory::inspect_apex_vert(Vert* const vertex,
					SampleVertTag* const vertex_tag)
{

//   printf("parent entity = %p %p %p\n",
// 	 dynamic_cast<RefVertex*>(vertex->get_parent_topology()),
// 	 dynamic_cast<RefEdge*>(vertex->get_parent_topology()),
// 	 dynamic_cast<RefFace*>(vertex->get_parent_topology()));

//  assert( dynamic_cast<RefVertex*>(vertex->getParentTopology()) );
  assert(vertex_tag->get_vertex() == vertex);

  //Check if the subsegs are present.
  if (!check_subsegs(vertex, vertex_tag))
    return false;

  //Also need to check the vertices on the ring are correctly connected to subsegs.
  set<Vert*> bdry_ring_verts;
  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		BdryRingMapper(vertex, &bdry_ring_verts));

  if (static_cast<int>(bdry_ring_verts.size())
      != std::count_if(
	  bdry_ring_verts.begin(), bdry_ring_verts.end(),
	  GRUMMP::extern_mem_fun(&SamplingEntryFactory::check_subsegs, this)))
    return false;

  //For apex vertices, only proper connections to the subsegements is checked.
  //This is deemed sufficient. Disk correctness will eventually be
  //checked when testing the curve or surface vertices attached to the apex.

  //Because of possible non-manifold conditions at apices, it is difficult
  //to correctly establish how to test the disk condition anyways.

  //Future testing will tell whether this is sufficient or not.

  //We can still check the disk in "normal" conditions:

//   int num_bridges = vertex_tag->num_bridges();
//   assert( supported_apex_vertex(num_bridges, vertex, vertex_tag) );

//   switch (num_bridges) {

//   case 0:
//   case 1:
//     if( !check_disk(vertex, vertex_tag) )
//       build_topo_violation(vertex, vertex_tag);
//     break;

//   case 2:
//     set<SubsegBridge*>::const_iterator it = vertex_tag->get_bridges().begin();
//     SubsegBridge* bridge1 = *it, *bridge2 = *(++it);

//     if( !check_disk(vertex, vertex_tag, bridge1, bridge2) )
//       build_topo_violation(vertex, vertex_tag);
//     break;

//   }

  return true;

}

bool
SamplingEntryFactory::inspect_curv_vert(Vert* const vertex,
					SampleVertTag* const vertex_tag)
{

  assert(vertex_tag->get_vertex() == vertex);
//  assert( dynamic_cast<RefEdge*>(vertex->getParentTopology()) );
//  assert( vertex_tag->get_bridges().size() == 2 );

  //Check if subsegs are present on the ring faces' perimeter.
  //If so, make sure this vertex does not encroach the subsegs.
  //If so, mark the subsegment for split (it is encroached after all).
  vector<SubsegBridge*> ring_bridges;
  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		RingBridgeMapper(vertex, m_sampler, &ring_bridges));

  vector<SubsegBridge*> encroached_bridges;
  std::remove_copy_if(ring_bridges.begin(), ring_bridges.end(),
		      std::back_inserter(encroached_bridges),
		      BridgeNotEncroached(vertex));

  //treat the encroached bridge as a missing subsegment.
  std::for_each(
      encroached_bridges.begin(),
      encroached_bridges.end(),
      GRUMMP::extern_mem_fun(&SamplingEntryFactory::build_missing_subseg,
			     this));

  //Check if the subsegs attached to this vertex are present.
  if (!check_subsegs(vertex, vertex_tag))
    return false;

  //Also need to check the vertices on the ring are correctly connected to subsegs.
  set<Vert*> bdry_ring_verts;
  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		BdryRingMapper(vertex, &bdry_ring_verts));

  if (static_cast<int>(bdry_ring_verts.size())
      != std::count_if(
	  bdry_ring_verts.begin(), bdry_ring_verts.end(),
	  GRUMMP::extern_mem_fun(&SamplingEntryFactory::check_subsegs, this)))
    return false;

  //Now check the disk around this curve vertex. Disk will start a
  //one bridge and stop at the other.
  set<SubsegBridge*>::const_iterator it = vertex_tag->get_bridges().begin();
  SubsegBridge* bridge1 = *it, *bridge2 = *(++it);
  assert(++it == vertex_tag->get_bridges().end());

  if (!check_disk(vertex, vertex_tag, bridge1, bridge2)) {
    build_topo_violation(vertex, vertex_tag);
    return false;
  }

  return true;

}

bool
SamplingEntryFactory::inspect_surf_vert(Vert* const vertex,
					SampleVertTag* const vertex_tag)
{

  assert(vertex_tag->get_vertex() == vertex);
  assert(vertex_tag->get_bridges().empty());
  assert(
      dynamic_cast<RefFace*>(vertex->getParentTopology())
	  == m_sampler->m_surface);

  //Need to check the vertices on the ring are correctly connected to subsegs.
  set<Vert*> bdry_ring_verts;
  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		BdryRingMapper(vertex, &bdry_ring_verts));

  if (static_cast<int>(bdry_ring_verts.size())
      != std::count_if(
	  bdry_ring_verts.begin(), bdry_ring_verts.end(),
	  GRUMMP::extern_mem_fun(&SamplingEntryFactory::check_subsegs, this)))
    return false;

  if (!check_disk(vertex, vertex_tag)) {
    build_topo_violation(vertex, vertex_tag);
    return false;
  }

  return true;

}

bool
SamplingEntryFactory::check_subsegs(Vert* const vertex)
{

  SampleVertTag* vertex_tag = m_sampler->get_vert_tag(vertex);
  return check_subsegs(vertex, vertex_tag);

}

bool
SamplingEntryFactory::check_subsegs(Vert* const vertex,
				    SampleVertTag* const vertex_tag)
{

  assert(vertex_tag->get_vertex() == vertex);
//  assert( dynamic_cast<RefVertex*>(vertex->getParentTopology()) ||
//	  dynamic_cast<RefEdge*>  (vertex->getParentTopology()) ||
//	  dynamic_cast<RefFace*>  (vertex->getParentTopology()));

  //The vertex lies on the boundary of the surface, therefore, it has
  //subsegment bridges attached to it. The set of bridges attached to
  //vertex is stored in the vertex tag. Here we inspect the restricted 1-star
  //around vertex to see whether these bridges are all present in the sampled
  //surface. If they are not, then we need to split the subsegment.

  Vert* other_vertex;
  SubsegBridge* bridge;

  set<SubsegBridge*>::const_iterator itb = vertex_tag->get_bridges().begin(),
      itb_end = vertex_tag->get_bridges().end();

  bool ret_value = true;

  for (; itb != itb_end; ++itb) {

    bridge = *itb;
    assert(bridge->get_parent()->is_not_deleted());
    assert(bridge->get_vert(0) == vertex || bridge->get_vert(1) == vertex);

    if (m_bad_bridges.find(bridge) != m_bad_bridges.end())
      return false;
    if (!m_checked_bridges.insert(bridge).second)
      continue;

    other_vertex =
	vertex == bridge->get_vert(0) ?
	    bridge->get_vert(1) : bridge->get_vert(0);
    assert(other_vertex != vertex);

    //Check if the edge represented by the bridge is present in the ring.
    if (!edge_in_ring(vertex, other_vertex, vertex_tag->get_ring_faces())) {
      m_bad_bridges.insert(bridge);
      build_missing_subseg(bridge);
      ret_value = false;
    }

  } //END CHECKING BRIDGES.

  return ret_value;

}

bool
SamplingEntryFactory::check_disk(Vert* const vertex,
				 SampleVertTag* const vertex_tag,
				 SubsegBridge* beg_bridge,
				 SubsegBridge* end_bridge)
{

  if (vertex_tag->get_ring_faces().empty())
    return true;

  RingMap ring_map;
  RingMap::iterator it_ring;
  RingRange ring_range;

  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		RingMapper(vertex, &ring_map));

  Vert *beg_vert, *end_vert, *this_vert;
  FaceQueueEntry *this_face, *prev_face;

  if (beg_bridge) { //walk initialization for bdry vert.

    assert(end_bridge);
//    assert(vertex->getVertType() == Vert::eBdryApex ||
//	   vertex->getVertType() == Vert::eBdryCurve);

    beg_vert =
	vertex == beg_bridge->get_vert(0) ?
	    beg_bridge->get_vert(1) : beg_bridge->get_vert(0);
    end_vert =
	vertex == end_bridge->get_vert(0) ?
	    end_bridge->get_vert(1) : end_bridge->get_vert(0);
    this_vert = beg_vert;

    ring_range = ring_map.equal_range(beg_vert);
    if (std::distance(ring_range.first, ring_range.second) != 1)
      return false;

    prev_face = this_face = ring_range.first->second;

  }

  else { //walk initialization for interior vert.

    assert(!end_bridge);
    assert(vertex->getVertType() == Vert::ePseudoSurface);

    beg_vert = end_vert = this_vert = ring_map.begin()->first;

    ring_range = ring_map.equal_range(beg_vert);
    if (std::distance(ring_range.first, ring_range.second) != 2)
      return false;

    it_ring = ring_range.first;
    this_face = it_ring->second;
    prev_face = (++it_ring)->second;

  }

  ring_map.erase(ring_range.first);

  while (true) { //the walk around the vertex.

    if (ring_map.empty())
      return false; //did not reach end_vert, but no more faces.

    if (!edge_checked(vertex, this_vert)) { //do not do this again if done previously.
      if (!valid_dihedral_angle(prev_face, this_face, vertex, this_vert))
	return false;
      if (!valid_normal_deviation(vertex, this_vert, vertex_tag))
	return false;
    }

    if (!valid_edge_length(vertex, this_vert, vertex_tag))
      return false;

    assert(vertex_tag->get_vertex() == vertex);
    this_vert = face_other_vert(this_face->getFace(), vertex, this_vert);

    if (this_vert == end_vert) {
      ring_range = ring_map.equal_range(this_vert);
      if (std::distance(ring_range.first, ring_range.second) != 1)
	return false;
      ring_map.erase(ring_range.first, ring_range.second);
      break;
    }

    ring_range = ring_map.equal_range(this_vert);
    if (std::distance(ring_range.first, ring_range.second) != 2)
      return false;

#ifndef NDEBUG
    it_ring = ring_range.first;
    assert(it_ring->second == this_face || (++it_ring)->second == this_face);
#endif

    it_ring = ring_range.first;
    prev_face = this_face;
    this_face =
	it_ring->second != prev_face ? it_ring->second : (++it_ring)->second;

    ring_map.erase(ring_range.first, ring_range.second);

  }

  // Long-winded, but possible to set breakpoints more easily.
  if (ring_map.empty())
    return true;
  else
    return false;
//  return (ring_map.empty());

}

bool
SamplingEntryFactory::check_peripheral_edges(const VertWithTag& vert_and_tag)
{

  return check_peripheral_edges(vert_and_tag.first, vert_and_tag.second);

}

bool
SamplingEntryFactory::check_peripheral_edges(Vert* const vertex,
					     SampleVertTag* vertex_tag)
{

  set<FactoryEdge> periph_edges;
  std::for_each(vertex_tag->get_ring_faces().begin(),
		vertex_tag->get_ring_faces().end(),
		PeriphEdgeID(vertex, &periph_edges));

  Vert* vert0, *vert1;
  SampleVertTag* tag0, *tag1;
  vector<FaceQueueEntry*> edge_faces;

  set<FactoryEdge>::iterator it = periph_edges.begin(), it_end =
      periph_edges.end();

  for (; it != it_end; ++it) {

    assert(m_checked_edges.count(*it) == 0);

    edge_faces.clear();

    vert0 = it->get_vert(0);
    vert1 = it->get_vert(1);

    //Edge coincides with a subseg. Nothing more to do.
    if ((vert0->getVertType() == Vert::eBdryApex
	|| vert0->getVertType() == Vert::eBdryCurve)
	&& (vert1->getVertType() == Vert::eBdryApex
	    || vert1->getVertType() == Vert::eBdryCurve))
      continue;

    tag0 = m_sampler->get_vert_tag(vert0);
    tag1 = m_sampler->get_vert_tag(vert1);

    std::set_intersection(tag0->get_ring_faces().begin(),
			  tag0->get_ring_faces().end(),
			  tag1->get_ring_faces().begin(),
			  tag1->get_ring_faces().end(),
			  std::inserter(edge_faces, edge_faces.begin()));

    if (std::distance(
	edge_faces.begin(),
	std::remove_if(edge_faces.begin(), edge_faces.end(),
		       std::mem_fun(&FaceQueueEntry::isDeleted))) != 2)
      return false;

    if (!valid_dihedral_angle(*(edge_faces.begin()), *(edge_faces.begin() + 1),
			      vert0, vert1))
      return false;

  }

  return true;

}

bool
SamplingEntryFactory::valid_normal_deviation(Vert* const vert0,
					     Vert* const vert1,
					     SampleVertTag* tag0,
					     SampleVertTag* tag1) const
{

  if (!m_norm_samp)
    return true;

#ifndef NDEBUG
  if (tag0)
    assert(tag0->get_vertex() == vert0);
  if (tag1)
    assert(tag1->get_vertex() == vert1);
#endif

  if (!tag0)
    tag0 = m_sampler->get_vert_tag(vert0);
  if (!tag1)
    tag1 = m_sampler->get_vert_tag(vert1);
  assert(tag0->normal_is_set());
  assert(tag1->normal_is_set());

//   double norm_var = *(tag0->get_normal()) % *(tag1->get_normal());
//   if(norm_var <= 0.)
//     printf("invalid norm deviation: %lf\n", norm_var);
#ifndef NDEBUG
  // This clause is a much less compact way of saying what's in the #else, but
  // with the advantage that a debugger can set breakpoints in branches of it.
  CubitVector norm0 = *(tag0->get_normal());
  CubitVector norm1 = *(tag1->get_normal());
  double dot = norm0 % norm1;
  if (dot > 0) {
    return true;
  }
  else {
    return false;
  }
#else
  return ( *(tag0->get_normal()) % *(tag1->get_normal()) > 0. );
#endif

}

bool
SamplingEntryFactory::valid_dihedral_angle(FaceQueueEntry* const fqe0,
					   FaceQueueEntry* const fqe1,
					   Vert* const vert0,
					   Vert* const vert1) const
{

  //A dihedral angle between faces is considered "topologically"
  //valid if it is under 90 degrees (says I!).

  if (!m_dihed_samp)
    return true;
  if (fqe0->isEqual(fqe1))
    return true;

  assert(fqe0->hasVert(vert0) && fqe1->hasVert(vert1));
  assert(fqe1->hasVert(vert0) && fqe1->hasVert(vert1));

  //The face normal is pre-computed and stored in the face tags.
  RestrictedFaceTag *face_tag0 = m_sampler->get_face_tag(fqe0), *face_tag1 =
      m_sampler->get_face_tag(fqe1);

  //We have to get correct normal orientation based on edge from vert0 to vert1.
  double orient0, orient1;

  switch (fqe0->getVert(0) == vert0 ? 1 : fqe0->getVert(1) == vert0 ? 2 : 3)
    {
    case 1:
      orient0 = vert1 == fqe0->getVert(1) ? 1. : -1.;
      break;
    case 2:
      orient0 = vert1 == fqe0->getVert(2) ? 1. : -1.;
      break;
    case 3:
      orient0 = vert1 == fqe0->getVert(0) ? 1. : -1.;
      break;
    default:
      assert(0);
      break;
    }

  switch (fqe1->getVert(0) == vert0 ? 1 : fqe1->getVert(1) == vert0 ? 2 : 3)
    {
    case 1:
      orient1 = vert1 == fqe1->getVert(1) ? -1. : 1.;
      break;
    case 2:
      orient1 = vert1 == fqe1->getVert(2) ? -1. : 1.;
      break;
    case 3:
      orient1 = vert1 == fqe1->getVert(0) ? -1. : 1.;
      break;
    default:
      assert(0);
      break;
    }

//   double dihed =
//     (orient0 * face_tag0->get_face_normal()) %
//     (orient1 * face_tag1->get_face_normal());
//   if(dihed <= 0.)
//     printf("invalid dihed: %lf\n", dihed);

#ifndef NDEBUG
  CubitVector norm0 = face_tag0->get_face_normal();
  CubitVector norm1 = face_tag1->get_face_normal();
  double result = orient0 * orient1 * (norm0 % norm1);
  if (result > 0.) {
    return true;
  }
  else {
    return false;
  }
#else
  return ( (orient0 * face_tag0->get_face_normal()) %
      (orient1 * face_tag1->get_face_normal()) > 0. );
#endif
}

bool
SamplingEntryFactory::valid_edge_length(Vert* const vert0, Vert* const vert1,
					SampleVertTag* const tag0) const
{

  if (!m_size_samp)
    return true;

  assert(vert0 && vert1);
  assert(tag0);
  assert(tag0->get_vertex() == vert0);

  double curv = tag0->get_curvature();
  double edge_length = dDIST3D(vert0->getCoords(), vert1->getCoords());
  double radius = iFuzzyComp(curv, 0.) == 0 ? LARGE_DBL : 1 / curv;
  assert(iFuzzyComp(curv, 0.) != -1);

//   double edge_bound = m_radius_curv_ratio * edge_length;
//   if(radius <= edge_bound)
//     printf("invalid edge length: radius = %lf, edge bound = %lf\n", radius, edge_bound);

  return (radius > m_radius_curv_ratio * edge_length);

}

bool
SamplingEntryFactory::face_has_edge(const Face* const face,
				    const Vert* const vert0,
				    const Vert* const vert1) const
{

  return (face->hasVert(vert0) && face->hasVert(vert1));

}

bool
SamplingEntryFactory::edge_checked(Vert* const vert0, Vert* const vert1)
{

  return !m_checked_edges.insert(FactoryEdge(vert0, vert1)).second;

}

bool
SamplingEntryFactory::edge_in_ring(const Vert* const vert0,
				   const Vert* const vert1,
				   const FaceQueueEntrySet& ring_faces) const
{

  const FaceQueueEntry* fqe;
  FaceQueueEntrySet::const_iterator itf = ring_faces.begin(), itf_end =
      ring_faces.end();

  for (; itf != itf_end; ++itf) {
    fqe = *itf;
    if (fqe->isDeleted())
      continue;
    if (fqe->getFace()->hasVert(vert0) && fqe->getFace()->hasVert(vert1))
      return true;
  }

  return false;

}

Face*
SamplingEntryFactory::largest_face_in_ring(
    const FaceQueueEntrySet& ring_faces) const
{

  FaceQueueEntry* large_face = NULL;
  std::for_each(ring_faces.begin(), ring_faces.end(), LargestFace(large_face));

  assert(large_face);
  return large_face->getFace();

}

FaceQueueEntry*
SamplingEntryFactory::largest_unused_face_in_ring(
    const FaceQueueEntrySet& ring_faces) const
{

  assert(!ring_faces.empty());
  FaceQueueEntry* large_face = NULL;
  std::for_each(ring_faces.begin(), ring_faces.end(),
		LargestFace(large_face, &m_used_faces));
  return large_face;
}

Vert*
SamplingEntryFactory::face_other_vert(Face* const face, Vert* const vert0,
				      Vert* const vert1) const
{

  assert(face && vert0 && vert1);
  assert(vert0 != vert1);
  assert(face->getNumVerts() == 3);
  assert(face->hasVert(vert0) && face->hasVert(vert1));

  return (
      vert0 == face->getVert(0) ?
	  vert1 == face->getVert(1) ? face->getVert(2) : face->getVert(1)
      :
      vert0 == face->getVert(1) ?
	  vert1 == face->getVert(0) ? face->getVert(2) : face->getVert(0)
      :
      vert1 == face->getVert(1) ? face->getVert(0) : face->getVert(1));

}

void
SamplingEntryFactory::RingBridgeMapper::operator()(FaceQueueEntry* const fqe)
{

  assert(fqe);

  if (fqe->isDeleted())
    return;
  assert(fqe->hasVert(m_vertex));
  Vert *v0, *v1;

  switch (fqe->getVert(0) == m_vertex ? 1 : fqe->getVert(1) == m_vertex ? 2 : 3)
    {
    case 1:
      v0 = fqe->getVert(1);
      v1 = fqe->getVert(2);
      break;
    case 2:
      v0 = fqe->getVert(0);
      v1 = fqe->getVert(2);
      break;
    case 3:
      v0 = fqe->getVert(0);
      v1 = fqe->getVert(1);
      break;

    }

// #ifndef NDEBUG
//   if(v0->iVertType() == Vert::ePseudoSurface)
//     assert(v0->get_bridges.empty());
//   else if(v1->iVertType() == Vert::ePseudoSurface)
//     assert(v1->get_bridges.empty());
// #endif

  if (v0->getVertType() != Vert::ePseudoSurface
      && v1->getVertType() != Vert::ePseudoSurface) {

//    assert(v0->getVertType() == Vert::eBdryApex ||
//	   v0->getVertType() == Vert::eBdryCurve);
//    assert(v1->getVertType() == Vert::eBdryApex ||
//	   v1->getVertType() == Vert::eBdryCurve);

    SampleVertTag* v0_tag = m_sampler->get_vert_tag(v0), *v1_tag =
	m_sampler->get_vert_tag(v1);

    std::set_intersection(v0_tag->get_bridges().begin(),
			  v0_tag->get_bridges().end(),
			  v1_tag->get_bridges().begin(),
			  v1_tag->get_bridges().end(),
			  std::back_inserter(*m_bridges));

  }

}

bool
SamplingEntryFactory::BridgeNotEncroached::operator()(
    SubsegBridge* const bridge)
{

  assert(m_vertex);

  return !SubsegManager::vert_inside_or_on_subseg_ball(m_vertex,
						       bridge->get_parent());

}

void
SamplingEntryFactory::BdryRingMapper::operator()(FaceQueueEntry* const fqe)
{

  Vert* my_verts[2];

  switch (fqe->getVert(0) == m_vertex ? 1 : fqe->getVert(1) == m_vertex ? 2 : 3)
    {
    case 1:
      my_verts[0] = fqe->getVert(1);
      my_verts[1] = fqe->getVert(2);
      break;
    case 2:
      my_verts[0] = fqe->getVert(0);
      my_verts[1] = fqe->getVert(2);
      break;
    case 3:
      my_verts[0] = fqe->getVert(0);
      my_verts[1] = fqe->getVert(1);
      break;

    }

  for (int i = 0; i < 2; i++) {

    switch (my_verts[i]->getVertType())
      {
      case Vert::eBdryApex:
      case Vert::eBdryCurve:
	m_bdry_ring_verts->insert(my_verts[i]);
	break;
      default:
	break;
      }

  }

}

void
SamplingEntryFactory::PeriphEdgeID::operator()(FaceQueueEntry* const fqe)
{

  if (fqe->isDeleted())
    return;

  assert(fqe);
  assert(fqe->getNumVerts() == 3);
  assert(fqe->hasVert(m_vertex));

  switch (fqe->getVert(0) == m_vertex ? 1 : fqe->getVert(1) == m_vertex ? 2 : 3)
    {

    case 1:
      m_edges->insert(FactoryEdge(fqe->getVert(1), fqe->getVert(2)));
      break;
    case 2:
      m_edges->insert(FactoryEdge(fqe->getVert(0), fqe->getVert(2)));
      break;
    case 3:
      m_edges->insert(FactoryEdge(fqe->getVert(0), fqe->getVert(1)));
      break;

    }

}

void
SamplingEntryFactory::RingMapper::operator()(FaceQueueEntry* const fqe)
{

  if (fqe->isDeleted())
    return;

  assert(fqe->hasVert(m_vertex));
  assert(fqe->getNumVerts() == 3);

  switch (fqe->getVert(0) == m_vertex ? 1 : fqe->getVert(1) == m_vertex ? 2 : 3)
    {
    case 1:
      m_ring_map->insert(std::make_pair(fqe->getVert(1), fqe));
      m_ring_map->insert(std::make_pair(fqe->getVert(2), fqe));
      break;
    case 2:
      m_ring_map->insert(std::make_pair(fqe->getVert(0), fqe));
      m_ring_map->insert(std::make_pair(fqe->getVert(2), fqe));
      break;
    case 3:
      m_ring_map->insert(std::make_pair(fqe->getVert(0), fqe));
      m_ring_map->insert(std::make_pair(fqe->getVert(1), fqe));
      break;

    }

}

void
SamplingEntryFactory::LargestFace::operator()(FaceQueueEntry* const fqe)
{

  if (fqe->isDeleted())
    return;

  if (m_exclude_faces->find(fqe) != m_exclude_faces->end())
    return;

  double face_size = fqe->getFace()->calcSize();
  if (face_size > m_large_size) {
    m_large_size = face_size;
    m_fqe = fqe;
  }

}

//Following four are just debugging code.

bool
SamplingEntryFactory::valid_data_passed(const VertWithTag& vert_and_tag) const
{

  if (vert_and_tag.second->get_vertex() != vert_and_tag.first)
    return false;
  if (m_sampler->get_vert_tag(vert_and_tag.first) != vert_and_tag.second)
    return false;
  return true;

}

bool
SamplingEntryFactory::valid_data_passed(
    const vector<VertWithTag>& sample_verts) const
{

  assert(m_sampler);
  assert(m_sampler->restrict_delaunay_initialized());

  vector<VertWithTag>::const_iterator it = sample_verts.begin(), it_end =
      sample_verts.end();

  for (; it != it_end; ++it)
    if (!valid_data_passed(*it))
      return false;

  return true;

}

bool
SamplingEntryFactory::valid_data_passed(
    const vector<SurfaceSampler::RestrictFace>& restrict_faces) const
{

  assert(m_sampler);
  assert(m_sampler->restrict_delaunay_initialized());

  vector<SurfaceSampler::RestrictFace>::const_iterator it =
      restrict_faces.begin(), it_end = restrict_faces.end();

  for (; it != it_end; ++it) {

    FaceQueueEntry* fqe = it->first;

    if (fqe->isDeleted())
      continue;
    if (!fqe->getFace()->isRestrictedDelaunay()) {
      assert(!m_sampler->face_has_tag(fqe));
      continue;
    }
    Vert* v0 = fqe->getVert(0);
    Vert* v1 = fqe->getVert(1);
    Vert* v2 = fqe->getVert(2);

    if (!m_sampler->face_has_tag(fqe))
      return false;

    assert(fqe->getNumVerts() == 3);

    if (!m_sampler->vert_has_tag(v0))
      return false;
    if (!m_sampler->vert_has_tag(v1))
      return false;
    if (!m_sampler->vert_has_tag(v2))
      return false;

  }

  return true;

}

bool
SamplingEntryFactory::test_faces_have_vertex(Vert* const vertex) const
{

  FaceQueueEntrySet::const_iterator it = m_test_faces.begin(), it_end =
      m_test_faces.end();

  for (; it != it_end; ++it)
    if (!(*it)->getFace()->hasVert(vertex))
      return false;

  return true;

}

// bool SamplingEntryFactory::
// all_subsegs_present(const SurfaceSampler* sampler,
// 		    const vector<Vert*>& verts_to_check) {

//   if(verts_to_check.empty()) return false;
//   assert(sampler);

//   reset_factory();

//   return std::find_if(verts_to_check.begin(), verts_to_check.end(),
// 		      GRUMMP::extern_mem_fun(subsegs_are_present, this) );

// }

// bool SamplingEntryFactory::
// subsegs_are_present(const Vert* const vertex) {

//   SampleVertTag* vertex_tag = m_sampler->get_vert_tag(vertex);

//   return ( vertex_tag->get_bridges().end() ==
// 	   std::find_if( vertex_tag->get_bridges.begin(), vertex_tag->get_bridges.end(),
// 			 RingMissBridge(vertex_tag->get_ring_faces()) ) );

// }

// bool SamplingEntryFactory::RingMissBridge::
// operator (const SubsegBridge* const bridge) {

//   m_it     = m_ring_faces.begin();
//   m_it_end = m_ring_faces.end();

//   for( ; m_it != m_it_end; ++m_it) {
//     m_face = *m_it;
//     if( m_face->qDeleted() ) continue;
//     if( !m_face->qHasTwoVerts(bridge->get_beg_vert(),
// 			      bridge->get_end_vert()) ) return true;
//   }

//   return false;

// }
