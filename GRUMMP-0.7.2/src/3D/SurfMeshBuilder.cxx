#include "GR_SurfMeshBuilder.h"
#include "GR_SurfaceSampler.h"
#include "GR_SamplingEntryFactory.h"
#include "GR_SamplingQueue.h"
#include "GR_Subseg.h"
#include "GR_SubsegManager.h"
#include "GR_VolMesh.h"
#include "BasicTopologyEntity.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

#include <stdio.h>

//not so portable...
#include <algorithm>
#include <functional>
#include <map>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

using std::vector;
using std::set;
using std::map;
using std::multimap;
using std::pair;

SurfaceSampler*
SamplerContainer::get_sampler(int index) const
{

  assert(index >= 0 && index < num_samplers());
  return m_samplers[index];

}

SurfaceSampler*
SamplerContainer::get_sampler_of_surface(RefFace* const surface) const
{

  assert(m_sampler_indices.count(surface) == 1);
  int index = m_sampler_indices.find(surface)->second;
  assert(index >= 0 && index < static_cast<int>(m_samplers.size()));

  return m_samplers[index];

}

void
SamplerContainer::get_all_samplers(vector<SurfaceSampler*>& face_samplers) const
{

  assert(face_samplers.empty());
  std::copy(m_samplers.begin(), m_samplers.end(),
	    std::back_inserter(face_samplers));

}

void
SamplerContainer::get_all_samplers(
    std::vector<SamplerFacePair>& face_samplers) const
{

  assert(face_samplers.empty());
  assert(m_samplers.size() == m_sampler_indices.size());

  face_samplers.resize(m_samplers.size());

  SamplerIndexMap::const_iterator it = m_sampler_indices.begin(), it_end =
      m_sampler_indices.end();
  unsigned int index = 0;

  for (; it != it_end; ++it) {
    index = it->second;
    face_samplers[index] = std::make_pair(m_samplers[index], it->first);
  }

}

void
SamplerContainer::destroy_face_samplers()
{

  if (m_samplers.empty())
    return;

  for (GR_index_t iSampler = 0; iSampler < m_samplers.size(); ++iSampler)
    if (m_samplers[iSampler])
      delete m_samplers[iSampler];

  m_samplers.clear();
  m_sampler_indices.clear();

}

SurfMeshBuilder::SurfMeshBuilder(const SurfMeshBuilder&) :
    m_entity(NULL), m_modeler_type(UNKNOWN), m_subseg_manager(NULL), m_sampling_queue(
    NULL), m_face_samplers(), m_bls(NULL)
{
  assert(0);
}

SurfMeshBuilder&
SurfMeshBuilder::operator=(const SurfMeshBuilder&)
{
  assert(0);
  return (*this);
}

SurfMeshBuilder::SurfMeshBuilder(
    BasicTopologyEntity * entity, const ModelerType mod_type,
    const double allowed_TVT, const VolMesh * backgroundMesh,
    const CurveDiscretizer::CDType subsegSamplingType) :
    m_entity(entity), m_modeler_type(mod_type), m_sampling_queue(
	new SamplingQueue()), m_face_samplers(), m_bls(NULL)
{

  assert(mod_type != UNKNOWN);

  //	Initializing the mesh builder and the samplers.
  if (backgroundMesh)
    m_bls = new BaryLengthScale3D(backgroundMesh);

  m_subseg_manager = new SubsegManager(entity, allowed_TVT, backgroundMesh,
				       subsegSamplingType);

  init_surface_samplers(m_modeler_type);
  init_singular_points();
  init_subseg_bridges();
  init_sampling_queue();
}

SurfMeshBuilder::SurfMeshBuilder(SurfMesh * surf, BasicTopologyEntity * entity,
				 const ModelerType mod_type) :
    m_entity(entity), m_modeler_type(mod_type), m_subseg_manager(
	new SubsegManager(surf, entity)), m_sampling_queue(new SamplingQueue()), m_face_samplers(), m_bls(
    NULL)
{
  assert(mod_type != UNKNOWN);

  init_surface_samplers(m_modeler_type);
  init_singular_points();
  init_subseg_bridges();
  init_sampling_queue();

}

SurfMeshBuilder::~SurfMeshBuilder()
{
  m_sampling_queue->clean_queue();
  m_face_samplers.destroy_face_samplers();

  //Delete other stuff.
  if (m_subseg_manager)
    delete m_subseg_manager;
  if (m_sampling_queue)
    delete m_sampling_queue;
  if (m_bls)
    delete m_bls;
}

void
SurfMeshBuilder::get_face_samplers(vector<SurfaceSampler*>& face_samplers) const
{

  m_face_samplers.get_all_samplers(face_samplers);

}

void
SurfMeshBuilder::get_subsegs(vector<Subseg*>& subsegs) const
{

  assert(m_subseg_manager);
  m_subseg_manager->get_subsegments(subsegs);

}
void
SurfMeshBuilder::get_subsegs(set<Subseg*>& subsegs) const
{
  vector<Subseg*> vec_subsegs;
  assert(m_subseg_manager);
  m_subseg_manager->get_subsegments(vec_subsegs);
  for (std::vector<Subseg*>::iterator it = vec_subsegs.begin();
      it != vec_subsegs.end(); ++it) {
    subsegs.insert(*it);
  }

}

void
SurfMeshBuilder::get_subsegs_from_bridges(set<Subseg*>& subsegs) const
{
  set<Subseg*> tempSubsegs;
  get_subsegs(tempSubsegs);
  std::set<Subseg*>::iterator it_sub = tempSubsegs.begin(), it_sub_end =
      tempSubsegs.end();

  for (; it_sub != it_sub_end; ++it_sub) {

    Subseg * oldSubseg = *it_sub;
    std::vector<SubsegBridge*> subseg_bridges;
    oldSubseg->get_children(subseg_bridges);
    std::vector<SubsegBridge*>::iterator it_bridge = subseg_bridges.begin(),
	it_bridge_end = subseg_bridges.end();

    for (; it_bridge != it_bridge_end; ++it_bridge) {
      SubsegBridge* bridge = *it_bridge;
      Subseg * newSubseg;
      Vert* v0 = bridge->get_vert(0);
      Vert* v1 = bridge->get_vert(1);

      newSubseg = new Subseg(v0, v1, oldSubseg->get_beg_param(),
			     oldSubseg->get_end_param(),
			     oldSubseg->get_ref_edge());
      subsegs.insert(newSubseg);
    }
  }

}
void
SurfMeshBuilder::get_subseg_verts(vector<Vert*>& subseg_verts) const
{

  assert(m_subseg_manager);
  m_subseg_manager->get_vertices(subseg_verts);

}

void
SurfMeshBuilder::get_sample_verts(vector<Vert*>& sample_verts,
				  RefFace* const ref_face) const
{

  SurfaceSampler* sampler;

  if (ref_face) {
    sampler = m_face_samplers.get_sampler_of_surface(ref_face);
    sampler->get_sample_verts(sample_verts);
  }

  else {

    int num_samplers = m_face_samplers.num_samplers();

    for (int i = 0; i < num_samplers; i++) {
      sampler = m_face_samplers.get_sampler(i);
      sampler->get_sample_verts(sample_verts);
    }

  }

}

void
SurfMeshBuilder::get_restricted_faces(
    vector<FaceQueueEntry*>& restrict_faces) const
{

  int num_samplers = m_face_samplers.num_samplers();

  SurfaceSampler* sampler;
  vector<FaceQueueEntry*> fqes;

  for (int i = 0; i < num_samplers; i++) {
    sampler = m_face_samplers.get_sampler(i);
    sampler->get_restricted_faces(fqes);
  }

  std::copy(fqes.begin(), fqes.end(), std::back_inserter(restrict_faces));

#ifndef NDEBUG

  //Just making sure elements in restrict_faces are unique.
  FaceQueueEntrySet tmp_faces;
  std::copy(restrict_faces.begin(), restrict_faces.end(),
	    std::inserter(tmp_faces, tmp_faces.begin()));
  assert(tmp_faces.size() == restrict_faces.size());

#endif

}

void
SurfMeshBuilder::get_restricted_faces(
    std::vector<std::pair<FaceQueueEntry*, RefFace*>>& restrict_faces) const
{

  typedef vector<pair<SurfaceSampler*, RefFace*> > FaceSamplerVec;

  FaceSamplerVec all_samplers;
  m_face_samplers.get_all_samplers(all_samplers);

  FaceSamplerVec::iterator its = all_samplers.begin(), its_end =
      all_samplers.end();

  SurfaceSampler* sampler;
  RefFace* surface;
  vector<FaceQueueEntry*> fqes;

  for (; its != its_end; ++its) {

    fqes.clear();

    sampler = its->first;
    surface = its->second;
    //		sampler->print_surface_sampler();
    sampler->get_restricted_faces(fqes);

    for (vector<FaceQueueEntry*>::iterator itf = fqes.begin();
	itf != fqes.end(); ++itf)
      restrict_faces.push_back(std::make_pair(*itf, surface));
  }
  fqes.clear();

#ifndef NDEBUG

  FaceQueueEntrySet tmp_faces;
  vector<pair<FaceQueueEntry*, RefFace*> >::iterator it_tmp;

  for (it_tmp = restrict_faces.begin(); it_tmp != restrict_faces.end();
      ++it_tmp)
    tmp_faces.insert(it_tmp->first);

  assert(tmp_faces.size() == restrict_faces.size());

#endif

}
void
SurfMeshBuilder::print_restricted_faces() const
{

  typedef vector<pair<SurfaceSampler*, RefFace*> > FaceSamplerVec;

  FaceSamplerVec all_samplers;
  m_face_samplers.get_all_samplers(all_samplers);

  FaceSamplerVec::iterator its = all_samplers.begin(), its_end =
      all_samplers.end();

  SurfaceSampler* sampler;

  vector<FaceQueueEntry*> fqes;

  for (; its != its_end; ++its) {

    fqes.clear();

    sampler = its->first;

    sampler->print_surface_sampler();

  }

}

// TODO: A bunch of diagnostic output in here that doesn't use logMessage.
void
SurfMeshBuilder::sample_all_surfaces()
{
  logMessage(1, "Sampling all surfaces with %d samplers\n",
	     m_face_samplers.num_samplers());
  logMessage(1, "Sampling queue length: %d\n", m_sampling_queue->queue_size());
  GR_index_t num_insertions = 0;
  SamplingEntry* sampling_entry = NULL;
  vector<SamplingEntry*> new_entries;

  FaceTopoSamplingEntry* topo_entry = NULL;
  while (true) {

    //Grab the top entry from the sampling queue.
    sampling_entry = m_sampling_queue->top_entry();

    //Queue is empty. Terminates here!
    if (!sampling_entry) {
      purge();
      break;
    }

    if (sampling_entry->remove_entry()) {
      m_sampling_queue->pop_top_entry();
      continue;
    }

    if (sampling_entry->entry_deleted()) {
      assert(sampling_entry->get_entry_type() == SamplingEntry::BAD_TOPO);
      logMessage(3, "Bad topo\n");
      topo_entry = dynamic_cast<FaceTopoSamplingEntry*>(sampling_entry);
      assert(topo_entry);
      assert(topo_entry->get_vertex());
      assert(topo_entry->get_vertex_tag());
      FaceQueueEntrySet ring_faces =
	  topo_entry->get_vertex_tag()->get_ring_faces();
      FaceQueueEntrySet::iterator it_f = ring_faces.begin(), it_e =
	  ring_faces.end();
      for (; it_f != it_e; ++it_f) {
	FaceQueueEntry* fqe = (*it_f);
	if (!fqe->getFace()->hasVert(topo_entry->get_vertex()))
	  topo_entry->get_vertex_tag()->remove_ring_face(fqe);

      }

      new_entries.clear();
      SamplingEntryFactory::instance()->build_entries(
	  m_face_samplers.get_sampler_of_surface(topo_entry->get_surface()),
	  std::make_pair(topo_entry->get_vertex(),
			 topo_entry->get_vertex_tag()),
	  &new_entries);
      std::for_each(
	  new_entries.begin(),
	  new_entries.end(),
	  std::bind1st(std::mem_fun(&SurfMeshBuilder::add_to_sampling_queue),
		       this));

      topo_entry->mark_for_removal();
      {
	char fileName[1024];
	snprintf(fileName, 1023, "sampled-surf-%d.vtk", num_insertions);
	output_stitched_surface(fileName);
      }
      continue;

    }

    // hacky way to get it to stop, but seems to work if the queue gets stuck
    // basically shifts the burden of the surface mesh onto refinement.
    num_insertions += split_sampling_entry(sampling_entry);
    //		if(sampling_entry->entry_deleted()){
    //			//			m_sampling_queue->pop_top_entry();
    //		} else {
    sampling_entry->incrementAttemptNum();
    //		}
    //		if(sampling_entry->getAttemptNum() == 4){
    //			face_entry =
    //					dynamic_cast<FaceSamplingEntry*>(sampling_entry);
    //			if(face_entry)
    //				face_entry->getFace()->lock();
    //		}
    //		if(sampling_entry->getAttemptNum() == 5){
    //			m_sampling_queue->pop_top_entry();
    //		}

    if (num_insertions % 100 == 0
	&& num_insertions > 0 /*&& num_insertions < 1000*/) {
      logMessage(1,
		 "number of sampling vertices inserted = %d, queue size %d\n",
		 num_insertions, m_sampling_queue->queue_size());
      char fileName[1024];
      snprintf(fileName, 1023, "sampled-surf-%d.vtk", num_insertions);
      output_stitched_surface(fileName);
    }
    //output_stitched_surface("final-sampled.vtk");
    //			writeVTKLegacyWithoutPurge(*m_face_samplers.get_sampler(0)->m_mesh,"tempmesh");

//		}

  }

  logMessage(1, "number of sampling vertices inserted = %d\n", num_insertions);
}

void
SurfMeshBuilder::refine_by_background_mesh(const VolMesh * backgroundMesh)
{
  BaryLengthScale3D bls(backgroundMesh);
  m_sampling_queue->empty_queue();

  int nInsertions;
  do {
    nInsertions = 0;
    FaceVec restricted_faces;
    get_restricted_faces(restricted_faces);
    for (FaceVec::iterator it = restricted_faces.begin();
	it != restricted_faces.end(); ++it) {

      FaceQueueEntry * fqe = it->first;
      if (!fqe->isValid() || fqe->isDeleted())
	continue;
      Face * face = fqe->getFace();
      if (!face->isRestrictedDelaunay())
	continue;
      double radius, current_size, target_size;
      double adCircCent[3];
      calcCircumcenter3D(face->getVert(0), face->getVert(1), face->getVert(2),
			 adCircCent);
      radius = dDIST3D(adCircCent, face->getVert(0)->getCoords());
      current_size = sqrt(2.) * radius;

      target_size = (1. / 3.)
	  * (bls.getLengthScale(face->getVert(0))
	      + bls.getLengthScale(face->getVert(1))
	      + bls.getLengthScale(face->getVert(2)));

      if (current_size > target_size) {
	// we need to split the face
	nInsertions += split_face_and_sample(fqe, it->second);
      }
    }
    logMessage(1, "Inserted %d verts this pass\n", nInsertions);
  }
  while (nInsertions > 0);
}
void
SurfMeshBuilder::sample_all_surfaces_from_surf()
{

  unsigned int num_insertions = 0;

  SamplingEntry* sampling_entry = NULL;
  vector<SamplingEntry*> new_entries;

  FaceTopoSamplingEntry* topo_entry = NULL;

  while (true) {
    //Grab the top entry from the sampling queue.
    sampling_entry = m_sampling_queue->top_entry();
//		printf("queue size %d\n",m_sampling_queue->queue_size());
    //Queue is empty. Terminates here!
    if (!sampling_entry) {
      purge();
      break;
    }

    //If the top entry is deleted, then we need to pop it from the queue.
    if (sampling_entry->remove_entry()) {
      m_sampling_queue->pop_top_entry();
      continue;
    }

    if (sampling_entry->entry_deleted()) {

      assert(sampling_entry->get_entry_type() == SamplingEntry::BAD_TOPO);

      topo_entry = dynamic_cast<FaceTopoSamplingEntry*>(sampling_entry);
      assert(topo_entry);
      assert(topo_entry->get_vertex());
      assert(topo_entry->get_vertex_tag());
      FaceQueueEntrySet ring_faces =
	  topo_entry->get_vertex_tag()->get_ring_faces();
      FaceQueueEntrySet::iterator it_f = ring_faces.begin(), it_e =
	  ring_faces.end();
      for (; it_f != it_e; ++it_f) {
	FaceQueueEntry* fqe = (*it_f);
	if (!fqe->getFace()->hasVert(topo_entry->get_vertex()))
	  topo_entry->get_vertex_tag()->remove_ring_face(fqe);

      }

      new_entries.clear();
      SamplingEntryFactory::instance()->build_entries(
	  m_face_samplers.get_sampler_of_surface(topo_entry->get_surface()),
	  std::make_pair(topo_entry->get_vertex(),
			 topo_entry->get_vertex_tag()),
	  &new_entries);

      std::for_each(
	  new_entries.begin(),
	  new_entries.end(),
	  std::bind1st(std::mem_fun(&SurfMeshBuilder::add_to_sampling_queue),
		       this));

      topo_entry->mark_for_removal();

      continue;

    }
    // intentially avoid bad_topology entries.
    if (sampling_entry->get_entry_type() == SamplingEntry::BAD_TOPO)
      m_sampling_queue->pop_top_entry();
    else
      num_insertions += split_sampling_entry(sampling_entry);

  }

  logMessage(1, "number of sampling vertices inserted = %d\n", num_insertions);
}

void
SurfMeshBuilder::insert_interior_points_into_sampled_surf(
    std::vector<CubitVector> & points)
{
  // given a point we need to figure out where to insert it
  // which potentially means going between surface samplers
  // the subseg thing is different, thats going to be harder.
  // need to check subsegs first, and if not, then split
  // TODO OCTREE THIS
  SurfaceSampler * sampler = m_face_samplers.get_sampler(0);

  for (GR_index_t iPoint = 0; iPoint < points.size(); iPoint++) {
    Vert vertex;
    vertex.setCoords(3, points[iPoint]);
    vertex.setParentTopology(sampler->m_surface);
    vertex.setType(Vert::ePseudoSurface);
    sampler->insert_arbitrary_point(&vertex, pCInvalidCell, NULL);

  }

}
void
SurfMeshBuilder::insert_boundary_points_into_sampled_surf(
    std::vector<CubitVector> & points)
{
  // given a point we need to figure out where to insert it
  // which potentially means going between surface samplers
  // the subseg thing is different, thats going to be harder.
  // need to check subsegs first, and if not, then split
  SurfaceSampler * ss = m_face_samplers.get_sampler(0);
  RefFace *surface = ss->m_surface;
  DLIList<RefEdge*> refedges;
  surface->ref_edges(refedges);
  DLIList<RefVertex*> refvertices;
  surface->ref_vertices(refvertices);
  for (GR_index_t iPoint = 0; iPoint < points.size(); iPoint++) {
    CubitVector & point = points[iPoint];
    bool qRefVert = false;
    for (int iRefVert = 0; iRefVert < refvertices.size(); iRefVert++) {
      RefVertex * refVert = refvertices.get_and_step();

      if ((refVert->coordinates() - point).length() < 1.e-6) {
	printf("ref point %f %f %f\n", point[0], point[1], point[2]);
	qRefVert = true;
	break;
      }

    }
    if (!qRefVert) {
      for (int iRefEdge = 0; iRefEdge < refedges.size(); iRefEdge++) {
	RefEdge * refEdge = refedges.get_and_step();

	// need to find the face and force this point into it.
	//find face;

	CubitVector closest;
	refEdge->closest_point_trimmed(point, closest);
	if ((closest - point).length() < 1.e-6) {
	  printf("curve point %f %f %f\n", point[0], point[1], point[2]);

	  break;
	}
      }
    }
  }

}
int
SurfMeshBuilder::split_subseg_and_sample(Subseg* const subseg_to_split,
					 vector<SubsegSamplePair>* new_verts)
{

  assert(m_sampling_queue->is_empty());

  SamplingEntry* sampling_entry = NULL;
  int num_entries = 0;

  num_entries += split_subseg(subseg_to_split, new_verts);

  while (1) {

    sampling_entry = m_sampling_queue->top_entry();

    if (!sampling_entry) {
      break;
    }
    if (sampling_entry->entry_deleted()) {
      m_sampling_queue->pop_top_entry();
      continue;
    }

    num_entries += split_sampling_entry(sampling_entry, new_verts);

  }

  return num_entries;

}

int
SurfMeshBuilder::split_face_and_sample(FaceQueueEntry* const fqe,
				       RefFace* const insert_surface,
				       vector<SubsegSamplePair>* new_verts)
{

  assert(m_sampling_queue->is_empty());

  SamplingEntry* sampling_entry = NULL;

  int num_entries = 0;

  num_entries += split_restricted_face(fqe, insert_surface, 0, new_verts);

  while (1) {

    sampling_entry = m_sampling_queue->top_entry();

    if (!sampling_entry) {
      break;
    }
    if (sampling_entry->entry_deleted()) {
      m_sampling_queue->pop_top_entry();
      continue;
    }

    num_entries += split_sampling_entry(sampling_entry, new_verts);

  }

  return num_entries;

}

int
SurfMeshBuilder::split_sampling_entry(SamplingEntry* const sampling_entry,
				      vector<SubsegSamplePair>* new_verts)
{

  assert(!sampling_entry->entry_deleted());

  switch (sampling_entry->get_entry_type())
    {

    case SamplingEntry::SUBSEG:
      {
	SubsegSamplingEntry* subseg_entry =
	    dynamic_cast<SubsegSamplingEntry*>(sampling_entry);
	assert(subseg_entry);
	return split_subseg(subseg_entry, new_verts);

      }
    default:
      {

	FaceSamplingEntry* face_entry =
	    dynamic_cast<FaceSamplingEntry*>(sampling_entry);
	assert(face_entry);

	if (face_entry->getFace()->isLocked()) {
	  return 0;
	}
	return split_restricted_face(face_entry, new_verts);

      }

    }

}

int
SurfMeshBuilder::split_subseg(SubsegSamplingEntry* const subseg_entry,
			      vector<SubsegSamplePair>* new_verts)
{

#ifndef NDEBUG
  int num_verts = split_subseg(subseg_entry->get_bridge()->get_parent(),
			       new_verts);
  assert(subseg_entry->entry_deleted());
  return num_verts;
#else
  return split_subseg(subseg_entry->get_bridge()->get_parent(), new_verts);
#endif

}

int
SurfMeshBuilder::split_subseg(Subseg* const subseg,
			      vector<SubsegSamplePair>* new_verts)
{
  typedef vector<Vert*> VertexVec;
  typedef vector<Subseg*> SubsegVec;
  typedef vector<SubsegBridge*> BridgeVec;
  //	typedef vector<SamplingEntry*> SamEntVec;
  typedef vector<SurfaceSampler::RestrictFace> RestFaceVec;
  //First thing first, split the subsegment in m_subseg_manager. This split
  //might cause other splits (due to forced Delaunayization for instance).
  //Therefore, we need to keep track of what was created and what was deleted.

  assert(!subseg->is_deleted());

  VertexVec new_vertices;
  SubsegVec new_subsegs, old_subsegs;
  m_subseg_manager->split_subseg(subseg, &new_vertices, &new_subsegs,
				 &old_subsegs);
  //The following lacks elegance. When I designed the SubsegManager class, I
  //did not expect to have to use the new verts and subsegs outside the class.
  //The interface does not allow to know directly which new vertex is associated
  //to which new subsegment (should use a vector of struct or something of the
  //sort to do that). However, if one looks into SubsegManager::split_subseg_at_param(),
  //it is possible to see that the entries are inserted at the back of the vectors.
  //Therefore, the first entry in new_verts splits the first entry in old_subsegs and
  //is attached to the TWO first entries in new_subsegs. The vectors are all "sync-ed"
  //together and data is retrieved here as such. It is f-ugly, but it works well.
  //If someone wants to change the SubsegManager interface, be my guest! ;-)

  assert(!new_subsegs.empty());
  assert(2 * new_vertices.size() == new_subsegs.size());
  assert(2 * old_subsegs.size() == new_subsegs.size());

  //That being said, now that we have the new vertices, new and old subsegs, we
  //need to update the SubsegBridges. We will iterate through the containers, so
  //we have to assume that the assumption on the container's structure (above) holds.

  //Variables for split
  Subseg *new_sub1, *new_sub2, *old_sub;
  SubsegBridge *new_bridge1, *new_bridge2, *old_bridge;
  BridgeVec subseg_bridges;

  VertexVec::iterator it_vert = new_vertices.begin(), it_vert_end =
      new_vertices.end();
  SubsegVec::iterator it_new_sub = new_subsegs.begin(), it_old_sub =
      old_subsegs.begin();
  BridgeVec::iterator it_bridge, it_bridge_end;

  //Variables for rebuild
  int ii, index = -1;
  double curvature;
  Vert *new_vert;
  SampleVertTag *vert_tag;

  vector<RestFaceVec> new_rest_faces;
  new_rest_faces.reserve(30);
  vector<SamplingEntry*> new_entries;
  new_entries.reserve(200);

  vector<SurfaceSampler*> samplers;
  samplers.reserve(30);
  vector<Vert*> inserted_verts;
  inserted_verts.reserve(30);
  vector<SampleVertTag*> inserted_tags;
  inserted_tags.reserve(30);

  //That was a bunch of variables. Now the idea is to use the old subseg i.e. the one that
  //was split and the old subseg bridges to update the Subseg/SubsegBridge data structure.
  //The orientation of the new subsegs and new bridges is kept consistent with the
  //orientation of the old subseg and bridges. The code itself should (I say should) be
  //clear enough to understand what is going on.

  do {

    new_vert = *(it_vert++);
    new_sub1 = *(it_new_sub++);
    new_sub2 = *(it_new_sub++);
    old_sub = *(it_old_sub++);

    assert(old_sub->is_deleted());

    assert(new_sub1->has_vert(new_vert));
    assert(new_sub2->has_vert(new_vert));
    assert(new_sub1->has_vert(old_sub->get_beg_vert()));
    assert(new_sub2->has_vert(old_sub->get_end_vert()));

    subseg_bridges.clear();
    old_sub->get_children(subseg_bridges);

    assert(new_sub1->get_end_vert() == new_vert);
    assert(new_sub2->get_beg_vert() == new_vert);
    assert(new_sub1->get_ref_edge() == new_sub2->get_ref_edge());
    curvature = compute_curvature(m_modeler_type, new_sub1->get_ref_edge(),
				  new_sub1->get_end_param());

    it_bridge = subseg_bridges.begin();
    it_bridge_end = subseg_bridges.end();

    for (; it_bridge != it_bridge_end; ++it_bridge) {

      ++index;
      new_rest_faces.resize(index + 1);
      assert(new_rest_faces[index].empty());

      old_bridge = *it_bridge;
      new_bridge1 = new SubsegBridge();
      new_bridge2 = new SubsegBridge();

      new_bridge1->set_parent(new_sub1);
      new_bridge2->set_parent(new_sub2);

      new_bridge1->set_ref_face(old_bridge->get_ref_face());
      new_bridge2->set_ref_face(old_bridge->get_ref_face());

      samplers.push_back(
	  m_face_samplers.get_sampler_of_surface(old_bridge->get_ref_face()));

      //			assert(old_bridge->get_vert(0)->getHintFace());
      //			assert(old_bridge->get_vert(1)->getHintFace());

      inserted_verts.push_back(
	  samplers[index]->insert_sample_point(new_vert,
					       old_bridge->get_vert(0),
					       old_bridge->get_vert(1),
					       NULL,
					       &new_rest_faces[index]));

      if (new_verts)
	new_verts->push_back(std::make_pair(new_vert, inserted_verts[index]));

      new_bridge1->set_vert(0, old_bridge->get_vert(0));
      new_bridge1->set_vert(1, inserted_verts[index]);
      new_bridge2->set_vert(0, inserted_verts[index]);
      new_bridge2->set_vert(1, old_bridge->get_vert(1));

      new_sub1->add_bridge(new_bridge1);
      new_sub2->add_bridge(new_bridge2);

      vert_tag = samplers[index]->get_vert_tag(old_bridge->get_vert(0));
      assert(vert_tag);
      vert_tag->remove_bridge(old_bridge);
      vert_tag->add_bridge(new_bridge1);

      vert_tag = samplers[index]->get_vert_tag(old_bridge->get_vert(1));
      assert(vert_tag);
      vert_tag->remove_bridge(old_bridge);
      vert_tag->add_bridge(new_bridge2);

      vert_tag = samplers[index]->get_vert_tag(inserted_verts[index]);
      assert(vert_tag);
      vert_tag->add_bridge(new_bridge1);
      vert_tag->add_bridge(new_bridge2);
      vert_tag->set_curvature(curvature);
      inserted_tags.push_back(vert_tag);

      //Compute the normal to the surface at the
      //insertion location and store it in the vertex tag.
      vert_tag->set_normal(
	  samplers[index]->normal_at(inserted_verts[index], 1.e-5));

    }

  }
  while (it_vert != it_vert_end);

  for (std::vector<RestFaceVec>::iterator itt = new_rest_faces.begin();
      itt != new_rest_faces.end(); ++itt) {

  }
  assert(it_new_sub == new_subsegs.end() && it_old_sub == old_subsegs.end());

  //That's it. The necessary vertices were inserted in the corresponding samplers.
  //Now build the sampling entries.

  for (ii = 0; ii <= index; ii++) {

    new_entries.clear();

    SamplingEntryFactory::instance()->build_entries(
	samplers[ii], std::make_pair(inserted_verts[ii], inserted_tags[ii]),
	new_rest_faces[ii], &new_entries);

    std::for_each(
	new_entries.begin(),
	new_entries.end(),
	std::bind1st(std::mem_fun(&SurfMeshBuilder::add_to_sampling_queue),
		     this));

  }

  return static_cast<int>(new_vertices.size());

}

int
SurfMeshBuilder::split_restricted_face(FaceSamplingEntry* const face_entry,
				       vector<SubsegSamplePair>* new_verts)
{

  //Splitting a restricted face is relatively straightforward. The face_entry
  //contains the information necessary to split the face. Notably, it allows to
  //access the RestrictedFaceTag for this entry, which contains the insertion
  //location.
  return split_restricted_face(face_entry->getFaceQueueEntry(),
			       face_entry->get_surface(),
			       face_entry->getAttemptNum(), new_verts);

}

int
SurfMeshBuilder::split_restricted_face(FaceQueueEntry* const fqe,
				       RefFace* const insert_surface,
				       int attemptNum,
				       vector<SubsegSamplePair>* new_verts)
{

  Face * face_to_split = fqe->getFace();
  assert(face_to_split->isRestrictedDelaunay());
  //Get the intersection data for this face.
  SurfaceSampler* face_sampler = m_face_samplers.get_sampler_of_surface(
      insert_surface);

  RestrictedFaceTag* face_tag = face_sampler->get_face_tag(fqe);
  assert(face_tag);
  //Project the intersection point with the facetted
  //model to the real curved surface.
  CubitVector insert_coord, insert_normal;
  face_tag->get_intersection().get_triangle()->project_to_surface(
      face_tag->get_intersection().get_intersect_coord(), &insert_coord,
      &insert_normal);

  //Before inserting the point, we need to make sure it is not too close to the
  //boundary. Since we do not have a mesh that can be used as a search structure,
  //we have to rely on the trees contained in the subseg manager.
  Subseg* encroach_subseg = m_subseg_manager->point_encroaching_on_subseg(
      insert_coord, insert_surface);

  if (encroach_subseg)
    return split_subseg(encroach_subseg, new_verts);

  //We are not encroaching, therefore we split.

  //Because SurfaceSampler::insert_sample_point inserts a copy of a known vertex,
  //we have to create a dummy vertex here and pass it to the sampler.
  Vert new_vertex;
  new_vertex.setCoords(3, insert_coord);
  new_vertex.setParentTopology(insert_surface);
  new_vertex.setType(Vert::ePseudoSurface);
  //Find the curvature at the insertion location. This curvature
  //will be stored in the vertex tag once the vertex is inserted.
  double curvature =
      face_tag->get_intersection().get_triangle()->curvature_at_coord(
	  face_tag->get_intersection().get_intersect_coord());
  vector<SurfaceSampler::RestrictFace> new_restrict_faces;
  Vert* insert_vert = face_sampler->insert_sample_point(&new_vertex,
							face_to_split,
							attemptNum,
							&new_restrict_faces);
  if (!insert_vert->isValid())
    return 0;
  //This is useful in the edge / face recovery stage.
  //	insert_vert->printVertInfo();
  if (new_verts)
    new_verts->push_back(std::make_pair(static_cast<Vert*>(NULL), insert_vert));

  //Set the proper data in the vertex tag. We store the normal to the surface
  //and the curvature of the surface at the vertex coords.
  SampleVertTag* vert_tag = face_sampler->get_vert_tag(insert_vert);
  assert(vert_tag);
  vert_tag->set_normal(insert_normal);
  vert_tag->set_curvature(curvature);

  //Build the new sampling entries and store them in the priority queue.
  vector<SamplingEntry*> new_sampling_entries;
  SamplingEntryFactory::instance()->build_entries(
      face_sampler, std::make_pair(insert_vert, vert_tag), new_restrict_faces,
      &new_sampling_entries);

  for (vector<SamplingEntry*>::iterator it = new_sampling_entries.begin();
      it != new_sampling_entries.end(); ++it)
    add_to_sampling_queue(*it);

  //Only one vertex can be inserted.
  return 1;

}
void
SurfMeshBuilder::add_to_sampling_queue(SamplingEntry* const sampling_entry)
{

  assert(sampling_entry);

  if (!sampling_entry->entry_deleted())
    m_sampling_queue->add_to_queue(sampling_entry);
  else
    delete sampling_entry;

}

void
SurfMeshBuilder::init_surface_samplers(const ModelerType modeler_type)
{

  if (!m_face_samplers.empty())
    m_face_samplers.destroy_face_samplers();
  assert(m_face_samplers.empty());

  //Get all the surfaces.
  assert(m_entity);
  DLIList<RefFace*> all_surfaces;
  if (strcmp(m_entity->get_class_name(), "Surface") == 0)
    all_surfaces.insert(dynamic_cast<RefFace*>(m_entity));
  else
    m_entity->ref_faces(all_surfaces);
  unsigned int i, num_surfaces = all_surfaces.size();
  RefFace* surface;

  //Build and store the sampler for each surface.
  for (i = 0; i < num_surfaces; i++) {

    surface = all_surfaces.next(i);
    m_face_samplers.m_sampler_indices.insert(std::make_pair(surface, i));
    m_face_samplers.m_samplers.push_back(
	new SurfaceSampler(surface, modeler_type));

  }

}

void
SurfMeshBuilder::init_singular_points()
{

  vector<Vert*> all_vertices;
  m_subseg_manager->get_vertices(all_vertices);

  vector<Vert*>::iterator it = all_vertices.begin(), it_end =
      all_vertices.end();

  for (; it != it_end; ++it) {

    Vert* vertex = *it;
    RefVertex* ref_vertex =
	dynamic_cast<RefVertex*>(vertex->getParentTopology());

    if (ref_vertex) {

      DLIList<RefEdge*> curves_attached_to_vert;
      ref_vertex->ref_edges(curves_attached_to_vert);
      if (curves_attached_to_vert.size() == 0)
	vFatalError(
	    "GRUMMP does not currently support geometries with singular points\n"
	    "i.e. points that are attached to less than two curves\n"
	    "e.g. the tip of a cone\n"
	    "I am aborting here, sorry!",
	    "SurfMeshBuilder:init_singular_points()");

    }

  }

}

void
SurfMeshBuilder::init_subseg_bridges()
{

  //Need to store what has been inserted to avoid inserting twice.
  typedef map<RefFace*, Vert*> SurfaceToVert;
  SurfaceToVert::iterator it_beg_vert, it_end_vert;

  typedef map<Vert*, SurfaceToVert> InsertMap;
  InsertMap inserted_verts;

  //Define other variables here.
  int num_surf;
  double curv_beg, curv_end;

  RefEdge* curve;
  RefFace* surface;
  DLIList<RefFace*> neigh_surfaces;

  Subseg* subseg;
  SubsegBridge* subseg_bridge;
  Vert* subseg_beg_vert, *subseg_end_vert, *bridge_beg_vert, *bridge_end_vert;
  SurfaceSampler* sampler;
  SampleVertTag* vert_tag;

  //Insert all subseg vertices, only once.
  //Also set the normal to the surface at the inserted location (store in tag).
  //Also set the curvature of the curve at the insertion location (store in tag).

  vector<Subseg*> all_subsegs;
  m_subseg_manager->get_subsegments(all_subsegs);
  vector<Subseg*>::iterator it_sub = all_subsegs.begin(), it_sub_end =
      all_subsegs.end();

  for (; it_sub != it_sub_end; ++it_sub) {

    subseg = *it_sub;
    curve = subseg->get_ref_edge();
    neigh_surfaces.clean_out();
    if (curve)
      curve->ref_faces(neigh_surfaces);
    else
      neigh_surfaces.append(subseg->get_ref_face());
    num_surf = neigh_surfaces.size();
    curv_beg = curv_end = -LARGE_DBL;

    subseg_beg_vert = subseg->get_beg_vert();
    subseg_end_vert = subseg->get_end_vert();

    for (int i = 0; i < num_surf; i++) {

      surface = neigh_surfaces.next(i);
      sampler = m_face_samplers.get_sampler_of_surface(surface);

      bridge_beg_vert = bridge_end_vert = NULL;

      SurfaceToVert &beg_vertices = inserted_verts[subseg_beg_vert],
	  &end_vertices = inserted_verts[subseg_end_vert];

      it_beg_vert = beg_vertices.find(surface);
      it_end_vert = end_vertices.find(surface);

      if (it_beg_vert == beg_vertices.end()) { //beg vert not inserted.

	if (it_end_vert == end_vertices.end()) { //end vert not inserted.
	    //					subseg_end_vert->printVertInfo();
	  bridge_end_vert = sampler->insert_sample_point(subseg_end_vert);
	  assert(sampler->m_mesh->isValid());
	  end_vertices.insert(std::make_pair(surface, bridge_end_vert));
	  if (curv_end < -1.e100)
	    curv_end = compute_curvature(m_modeler_type, curve,
					 subseg->get_end_param());
	  vert_tag = sampler->get_vert_tag(bridge_end_vert);
	  vert_tag->set_normal(sampler->normal_at(bridge_end_vert, 1.e-5));
	  vert_tag->set_curvature(curv_end);
	}
	else {
	  bridge_end_vert = it_end_vert->second;
	} //end vert inserted.
	assert(bridge_end_vert);
	bridge_beg_vert = sampler->insert_sample_point(subseg_beg_vert,
						       bridge_end_vert);
	assert(sampler->m_mesh->isValid());
	beg_vertices.insert(std::make_pair(surface, bridge_beg_vert));
	if (curv_beg < -1.e100)
	  curv_beg = compute_curvature(m_modeler_type, curve,
				       subseg->get_beg_param());
	vert_tag = sampler->get_vert_tag(bridge_beg_vert);
	vert_tag->set_normal(sampler->normal_at(bridge_beg_vert, 1.e-5));
	vert_tag->set_curvature(curv_beg);

      }

      else { //beg vert inserted.

	bridge_beg_vert = it_beg_vert->second;

	if (it_end_vert == end_vertices.end()) { //end vert not inserted.

	  bridge_end_vert = sampler->insert_sample_point(subseg_end_vert,
							 bridge_beg_vert);
	  assert(sampler->m_mesh->isValid());
	  end_vertices.insert(std::make_pair(surface, bridge_end_vert));
	  if (curv_end < -1.e100)
	    curv_end = compute_curvature(m_modeler_type, curve,
					 subseg->get_end_param());
	  vert_tag = sampler->get_vert_tag(bridge_end_vert);
	  vert_tag->set_normal(sampler->normal_at(bridge_end_vert, 1.e-5));
	  vert_tag->set_curvature(curv_end);
	}
	else {
	  bridge_end_vert = it_end_vert->second;
	} //end vert inserted.
      }

      assert(bridge_beg_vert && bridge_end_vert);

      //Now that both vertices have been inserted in the sampler,
      //update the subseg and subseg bridge for this surface.

      subseg_bridge = new SubsegBridge();
      subseg_bridge->set_ref_face(surface);
      subseg_bridge->set_parent(subseg);
      subseg_bridge->set_vert(0, bridge_beg_vert);
      subseg_bridge->set_vert(1, bridge_end_vert);

      subseg->add_bridge(subseg_bridge);

      sampler->attach_bridge(bridge_beg_vert, subseg_bridge);
      sampler->attach_bridge(bridge_end_vert, subseg_bridge);
    }
  }
}

void
SurfMeshBuilder::init_sampling_queue()
{
  if (!m_sampling_queue->is_empty())
    m_sampling_queue->empty_queue();

  int num_samplers = m_face_samplers.num_samplers();
  SurfaceSampler* sampler;

  vector<SamplingEntry*> queue_entries;
  vector<SurfaceSampler::RestrictFace> restrict_faces;
  vector<SurfaceSampler::VertWithTag> sample_verts;

  for (int i = 0; i < num_samplers; i++) {
    sampler = m_face_samplers.get_sampler(i);
//                printf("%d %p %d\n", i, sampler, sampler->m_mesh->getNumVerts());
    restrict_faces.clear();
    sample_verts.clear();
//                sampler->output_restricted_delaunay("rdel_before.vtk");

    sampler->init_restricted_delaunay(sample_verts, restrict_faces);
//                sampler->output_restricted_delaunay("rdel_after.vtk");
    queue_entries.clear();
    SamplingEntryFactory::instance()->build_entries(sampler, sample_verts,
						    restrict_faces,
						    &queue_entries);
    std::for_each(
	queue_entries.begin(),
	queue_entries.end(),
	std::bind1st(std::mem_fun(&SurfMeshBuilder::add_to_sampling_queue),
		     this));
  }
}

void
SurfMeshBuilder::purge()
{

  //Before doing anything else, we need
  //Must first remove deleted entries from the sampling queue.

  //
  //	//Purge contents of subseg manager.
  //	m_subseg_manager->purge();

  //	//Purge face samplers and update sampling queue.
  int num_samplers = m_face_samplers.num_samplers();
  SurfaceSampler* sampler;
  //
  //	map<Vert*, Vert*> vert_map, global_vert_map;
  //	map<Face*, Face*> face_map, global_face_map;
  //	map<Cell*, Cell*> cell_map, global_cell_map;
  //
  for (int i = 0; i < num_samplers; i++) {

    sampler = m_face_samplers.get_sampler(i);
    sampler->purge_sampler();
  }

  //
  //		std::copy(vert_map.begin(), vert_map.end(),
  //				std::inserter(global_vert_map, global_vert_map.begin()));
  //		std::copy(face_map.begin(), face_map.end(),
  //				std::inserter(global_face_map, global_face_map.begin()));
  //		std::copy(cell_map.begin(), cell_map.end(),
  //				std::inserter(global_cell_map, global_cell_map.begin()));
  //
  //	}
  //
  //	//Finally, modify the queue entries using the face and cell maps.
  //	m_sampling_queue->rebuild_queue(global_vert_map, global_face_map, global_cell_map);

}

double
SurfMeshBuilder::compute_curvature(ModelerType mod_type, RefEdge* const curve,
				   double param)
{

  switch (mod_type)
    {

    case FACET:
      {
	if (curve)
	  return SurfMeshBuilder::discrete_curvature_estimate(curve, param);
	else
	  return 0.;
	break;
      }
    case UNKNOWN:
      vFatalError("Modeler type not set properly.",
		  "SurfMeshBuilder::curvature()");
      break;

    case ACIS:
    case OCC:
    default:
      vFatalError("Only support for facet modeler is currently implemented.",
		  "SurfMeshBuilder::curvature()");
      break;

    }
  return DBL_MAX;
}

double
SurfMeshBuilder::discrete_curvature_estimate(RefEdge* const curve, double param)
{

  //We select two points on the curve:
  //one slightly before param and one slightly after.
  //Then compute the circumradius at param to approximate
  //the radius of curvature.

  double min_param = curve->start_param();
  double max_param = curve->end_param();
  assert(iFuzzyComp(min_param, max_param) == -1);
  assert(
      iFuzzyComp(param, min_param) != -1 && iFuzzyComp(param, max_param) != 1);

  double delta = 0.05 * (max_param - min_param);
  double param_lo, param_hi;

  if (iFuzzyComp(param - delta, min_param) == -1) {
    param_lo = min_param;
    param = param_lo + delta;
    param_hi = param + delta;
  }
  else if (iFuzzyComp(param + delta, max_param) == 1) {
    param_hi = max_param;
    param = param_hi - delta;
    param_lo = param - delta;
  }
  else {
    param_lo = param - delta;
    param_hi = param + delta;
  }

  CubitVector coord0, coord1, coord2;
  curve->position_from_u(param_lo, coord0);
  curve->position_from_u(param, coord1);
  curve->position_from_u(param_hi, coord2);

  double vec0[] =
    { coord0.x(), coord0.y(), coord0.z() };
  double vec1[] =
	{ coord1.x() - coord0.x(), coord1.y() - coord0.y(), coord1.z()
	    - coord0.z() };
  double vec2[] =
	{ coord2.x() - coord0.x(), coord2.y() - coord0.y(), coord2.z()
	    - coord0.z() };
  double vec3[3], circ_center[3];
  vCROSS3D(vec1, vec2, vec3);

  //Three points are colinear
  if (iFuzzyComp(dMAG3D_SQ(vec3), 0.) == 0)
    return 0.;

  double dot1 = 0.5 * dDOT3D(vec1, vec1);
  double dot2 = 0.5 * dDOT3D(vec2, vec2);
  solve3By3(vec1, vec2, vec3, dot1, dot2, 0., circ_center);

  circ_center[0] += coord0.x();
  circ_center[1] += coord0.y();
  circ_center[2] += coord0.z();

  double circ_radius = dDIST3D(vec0, circ_center);
  assert(iFuzzyComp(circ_radius, 0.) == 1);

  return (1. / circ_radius);

}

void
SurfMeshBuilder::output_stitched_surface(
    const char* const strBaseFileName) const
{

  FILE *pFOutFile = NULL;
  char strFileName[1024];
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s.vtk", strBaseFileName);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s", strBaseFileName);
  }
  pFOutFile = fopen(strFileName, "w");
  if (!pFOutFile)
    vFatalError("Cannot open output file for writing",
		"SurfMeshBuilder::output_restricted_delaunay");

  typedef vector<FaceQueueEntry*> FacePrint;
  typedef vector<Vert*> VertPrint;
  typedef map<Vert*, int> VertIndex;

  int counter = 1;

  FacePrint faces_to_print;
  FacePrint::iterator itf, itf_end;
  VertPrint verts_to_print;
  VertPrint::iterator itv, itv_end;
  VertIndex vert_indices;

  int num_samplers = m_face_samplers.num_samplers();
  SurfaceSampler* sampler;

  for (int i = 0; i < num_samplers; i++) {
    sampler = m_face_samplers.get_sampler(i);
    sampler->get_restricted_faces(faces_to_print);
  }

  itf = faces_to_print.begin();
  itf_end = faces_to_print.end();

  for (; itf != itf_end; ++itf) {

    Face* face = (*itf)->getFace();

    for (int i = 0; i < face->getNumVerts(); i++) {
      Vert* vert = face->getVert(i);
      if (vert_indices.insert(std::make_pair(vert, counter)).second) {
	verts_to_print.push_back(vert);
	++counter;
      }
    }
  }

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "Geom tri data\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n",
	  static_cast<int>(verts_to_print.size()));

  itv = verts_to_print.begin();
  itv_end = verts_to_print.end();
  logMessage(3, "Verts %lu to print \n", verts_to_print.size());
  for (; itv != itv_end; ++itv) {
    Vert* vert = *itv;
    fprintf(pFOutFile, "%lf %lf %lf\n", vert->x(), vert->y(), vert->z());
  }

  itf = faces_to_print.begin();

  unsigned int nTri = static_cast<int>(itf_end - itf);
  fprintf(pFOutFile, "Cells %d %d\n", nTri, nTri * 4);

  for (; itf != itf_end; ++itf) {
    Face* face = (*itf)->getFace();
    fprintf(pFOutFile, "3 %d %d %d\n",
	    vert_indices.find(face->getVert(0))->second - 1,
	    vert_indices.find(face->getVert(1))->second - 1,
	    vert_indices.find(face->getVert(2))->second - 1);
  }

  for (; itf != itf_end; ++itf) {
    FaceQueueEntry * fqe = *itf;
    if (fqe)
      delete fqe;
  }
  faces_to_print.clear();
  fprintf(pFOutFile, "CELL_TYPES %d\n", nTri);
  for (unsigned int i = 0; i < nTri; i++) {
    fprintf(pFOutFile, "5\n");
  }

  fclose(pFOutFile);

}
