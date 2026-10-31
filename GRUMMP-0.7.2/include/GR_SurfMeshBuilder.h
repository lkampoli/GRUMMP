#ifndef GR_SURFMESH_BUILDER_H
#define GR_SURFMESH_BUILDER_H 1

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include "GR_BaryCentricLengthScale.h"
#include "GR_CurveDiscretizer.h"
#include "GR_VolMesh.h"
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <utility>

class BasicTopologyEntity;
class CubitBox;
class CubitVector;
class FaceSamplingEntry;
class RefEdge;
class RefFace;
class SurfaceSampler;
class SamplingQueue;
class SamplingEntry;
class Subseg;
class SubsegManager;
class SubsegSamplingEntry;
class TetCell;
class Vert;
class VolMesh;

class SurfMeshBuilder;
//A class that defines a container to store surface samplers.

class SamplerContainer {

public:

  friend class SurfMeshBuilder;

  typedef std::pair<SurfaceSampler*, RefFace*> SamplerFacePair;

  typedef std::vector<SurfaceSampler*> SamplerVec;
  typedef GR_map<RefFace*, int> SamplerIndexMap;

private:

  SamplerVec m_samplers;
  SamplerIndexMap m_sampler_indices;

  //If build_all_samplers is true, then automatically build a sampler
  //for each surface in the geometric model. Otherwise, samplers need
  //to be built manually by calling build_face_sampler.
  SamplerContainer() :
  m_samplers(), m_sampler_indices() {
  }
  ~SamplerContainer() {
    destroy_face_samplers();
  }

  bool empty() const {
    return m_samplers.empty();
  }

  int num_samplers() const {
    return static_cast<int>(m_samplers.size());
  }

  SurfaceSampler* get_sampler(int index) const;
  SurfaceSampler* get_sampler_of_surface(RefFace* const surface) const;

  void get_all_samplers(std::vector<SurfaceSampler*>& face_samplers) const;
  void get_all_samplers(std::vector<SamplerFacePair>& face_samplers) const;

  void destroy_face_samplers();

};

class SurfMeshBuilder {

public:

  //Potentially available modelers.
  enum ModelerType {
    ACIS, FACET, OCC, GRUMMP, UNKNOWN
  };
  //Not used now, but a potential idea I do not want to forget...
  enum SamplingType {
    CURVATURE, NONE
  };

  //This typedef is used to return new vertices when a split is called from the outside.
  //First entry is a subseg vertex, second entry is the associated sampling vertex.
  //Obviously, when splitting a restricted face, no subseg vertex is inserted, the
  //first entry will therefore be NULL.
  typedef std::pair<Vert*, Vert*> SubsegSamplePair;

private:
  // the volume we are interested in sampling,
  // using this general form allows for arbitrary
  // surfaces as well (RefVolume or RefFace)
  BasicTopologyEntity * m_entity;

  //The type of modeler used.
  ModelerType m_modeler_type;

  //The object to manage the subsegments discretization of RefEdges.
  SubsegManager* m_subseg_manager;

  //Sampling queue. Made of many priority queues (different sampling criteria).
  SamplingQueue* m_sampling_queue;

  //Container for the face samplers.
  SamplerContainer m_face_samplers;

  // should I sample based on an underlying mesh
  BaryLengthScale3D * m_bls;

public:
  //Constructor - Initializes all the necessary stuff.
  //curv_sample true for curvature-based sampling,
  //sample_radius is a coefficient that controls curvature-based sampling.
  SurfMeshBuilder(BasicTopologyEntity * entity, const ModelerType mod_type =
		      UNKNOWN,
		  const double allowed_TVT = 0.25 * M_PI,
		  const VolMesh * backgroundMesh = NULL,
		  const CurveDiscretizer::CDType subsegSamplingType =
		      CurveDiscretizer::TANGENT);

  // An experiment, to build one off of a surface mesh
  SurfMeshBuilder(SurfMesh * surf, BasicTopologyEntity * entity,
		  const ModelerType mod_type);

  //Destructor
  ~SurfMeshBuilder();

  //Once the sampling queue is initialized, add sample vertices until
  //the sampling queue becomes empty.
  void
  sample_all_surfaces();

  // given a background mesh and a sampled surface,
  // refine the sampled surface based on a background mesh
  void
  refine_by_background_mesh(const VolMesh * backgroundMesh);
  // this is used if you initialize from a surface mesh
  void
  sample_all_surfaces_from_surf();

  // insert points into sampled surface
  void
  insert_interior_points_into_sampled_surf(std::vector<CubitVector> & points);
  void
  insert_boundary_points_into_sampled_surf(std::vector<CubitVector> & points);

  //Splits subseg_to_split and samples until m_sampling_queue becomes empty.
  int
  split_subseg_and_sample(Subseg* const subseg_to_split,
			  std::vector<SubsegSamplePair>* new_verts = NULL);

  //Split face_to_split (must be a restricted face) and samples
  //until m_sampling_queue becomes empty.
  int
  split_face_and_sample(FaceQueueEntry* const fqe_to_split,
			RefFace* const insert_surface,
			std::vector<SubsegSamplePair>* new_verts = NULL);

  //Returns the face sampler pointers from m_face_samplers in the face_samplers vector.
  void
  get_face_samplers(std::vector<SurfaceSampler*>& face_samplers) const;

  //Populates subseg_verts with the subseg vertices as stored in m_subseg_manager.
  void
  get_subseg_verts(std::vector<Vert*>& subseg_verts) const;

  //Populates the set with all the undeleted subsegments in m_subseg_manager.
  //Two versions
  void
  get_subsegs(std::vector<Subseg*>& subsegs) const;
  void
  get_subsegs(std::set<Subseg*>& subsegs) const;

  // this is the actual mapping, to the verts in the surface mesh
  void
  get_subsegs_from_bridges(std::set<Subseg*>& subsegs) const;
  //Populates sample_verts with the sample vertices of ref_face.
  //If the ref_face pointer is non-NULL, then the vertices for all sampled surfaces are returned.
  void
  get_sample_verts(std::vector<Vert*>& sample_verts, RefFace* const ref_face =
  NULL) const;

  //Populates the vector with all the restricted faces defined in the face samplers.
  void
  get_restricted_faces(std::vector<FaceQueueEntry*>& restrict_faces) const;
  //Populates the vector with Restricted Face / RefFace pairs.
  void
  get_restricted_faces(
      std::vector<std::pair<FaceQueueEntry*, RefFace*> >& restrict_faces) const;

  void
  print_restricted_faces() const;
  //Purges the meshes in all the samplers, and updates data structures.
  //Does not seem to work (entirely) correctly at the moment.
  void
  purge();

  //Outputs the complete surface mesh (coming from all samplers).
  void
  output_stitched_surface(const char* const filename) const;

  //The curvature at a given param on a curve.
  static double
  compute_curvature(ModelerType mod_type, RefEdge* const curve, double param);

  //Estimates the curvature of a curve point by computing the circumradius
  //of three points the neighborhood of param.
  static double
  discrete_curvature_estimate(RefEdge* const curve, double param);

  // gets the entity
  BasicTopologyEntity *
  getEntity() const
  {
    return m_entity;
  }

private:

  //Copy constructor and operator= private.
  SurfMeshBuilder(const SurfMeshBuilder&);
  SurfMeshBuilder&
  operator=(const SurfMeshBuilder&);

  //Adds sampling_entry to m_sampling_queue
  void
  add_to_sampling_queue(SamplingEntry* const sampling_entry);

  //Returns the top valid (undeleted) entry in the sampling queue.
  //If the queue is empty, returns NULL.
  SamplingEntry*
  get_top_sampling_entry() const;

  //Returns the number of new vertices created.
  int
  split_sampling_entry(SamplingEntry* const sampling_entry,
		       std::vector<SubsegSamplePair>* new_verts = NULL);

  int
  split_subseg(Subseg* const subseg_to_split,
	       std::vector<SubsegSamplePair>* new_verts = NULL);

  int
  split_restricted_face(FaceQueueEntry* const fqe, RefFace* const insert_surf,
			int attemptNum,
			std::vector<SubsegSamplePair>* new_verts = NULL);

  int
  split_subseg(SubsegSamplingEntry* const subseg_entry,
	       std::vector<SubsegSamplePair>* new_verts = NULL);

  int
  split_restricted_face(FaceSamplingEntry* const face_entry,
			std::vector<SubsegSamplePair>* new_verts = NULL);

  //Initializes the surface samplers. One sampler per surface.
  void
  init_surface_samplers(const ModelerType modeler_type);

  //Builds all the subseg bridges and insert their vertices in
  //the corresponding face samplers. Initializes the data structure. (SG, c 2007)
  //
  // Subseg bridges relate versions of a subsegment that are found on
  // different surfaces.  This is essential so that the samples on the two
  // surfaces remain consistent with each other.
  //
  // In the process of setting these up, steps are taken to protect small angles.
  void
  init_subseg_bridges();

  //Initializes the sampling queue with data coming from the SurfaceSamplers.
  void
  init_sampling_queue();

  //Right now does nothing (except calling vFatalError...).
  void
  init_singular_points();

  //STUFF TO IMPLEMENT.. maybe here...
  void
  protect_small_dihedral()
  {
    assert(0);
  }

  //Used to set the properties of the surface samplers.
  void
  sampler_curv_based(SurfaceSampler* const sampler, const bool curv_based,
		     const double ratio);

};

#endif
