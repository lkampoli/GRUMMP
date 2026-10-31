#ifndef GR_SAMPLINGENTRYFACTORY_H
#define GR_SAMPLINGENTRYFACTORY_H

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include "GR_SamplingQueue.h"
#include "GR_SurfaceSampler.h"
#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

class Face;
class SampleVertTag;
class SamplingEntry;
class SubsegBridge;
class SurfaceSampler;
class Vert;

class FactoryEdge {

  Vert *m_v0, *m_v1;

  FactoryEdge() :
      m_v0(NULL), m_v1(NULL)
  {
  }

public:

  typedef std::pair<Vert*, SampleVertTag*> VertAndTag;

  FactoryEdge(Vert* v0, Vert* v1) :
      m_v0(std::min(v0, v1)), m_v1(std::max(v0, v1))
  {
    assert(m_v0 != m_v1);
  }
  FactoryEdge(const FactoryEdge& edge) :
      m_v0(edge.m_v0), m_v1(edge.m_v1)
  {
  }
  FactoryEdge&
  operator=(const FactoryEdge& edge)
  {
    if (&edge != this) {
      m_v0 = edge.m_v0;
      m_v1 = edge.m_v1;
    }
    return *this;
  }
  virtual
  ~FactoryEdge()
  {
  }

  bool
  operator<(const FactoryEdge& edge) const
  {
    assert(m_v0 < m_v1 && edge.m_v0 < edge.m_v1);
    return (m_v0 == edge.m_v0 ? m_v1 < edge.m_v1 : m_v0 < edge.m_v0);

  }
  bool
  operator==(const FactoryEdge& edge) const
  {
    assert(m_v0 < m_v1 && edge.m_v0 < edge.m_v1);
    return (m_v0 == edge.m_v0 && m_v1 == edge.m_v1);
  }

  Vert*
  get_vert(int i) const
  {
    assert(i == 0 || i == 1);
    return (i == 0 ? m_v0 : m_v1);
  }

  void
  set_verts(Vert* const v0, Vert* const v1)
  {
    assert(v0 != v1);
    m_v0 = std::min(v0, v1);
    m_v1 = std::max(v0, v1);
  }

  bool
  face_has_edge(const Face* const face) const
  {
    return (face->hasVert(m_v0) && face->hasVert(m_v1));
  }

  bool
  valid() const
  {
    if (!this)
      return false;
    if (!m_v0 || !m_v1)
      return false;
    if (m_v0 >= m_v1)
      return false;
    return true;
  }

};

class SamplingEntryFactory {

private:

  const SurfaceSampler* m_sampler;
  std::vector<SamplingEntry*>* m_new_entries;

  std::set<FactoryEdge> m_checked_edges;   //edges already checked.
  std::set<SubsegBridge*> m_bad_bridges; //bridges already checked and found to be bad.
  std::set<SubsegBridge*> m_checked_bridges; //ALL bridges already checked.
  FaceQueueEntrySet m_test_faces; //faces to test (all new restrict faces).
  FaceQueueEntrySet m_used_faces; //faces with entry already build.

  //parameters controlling sampling:

  bool m_dihed_samp, m_norm_samp;
  bool m_size_samp, m_qual_samp;
  double m_radius_edge_bound; //coefficient for quality sampling.
  double m_radius_curv_ratio; //coefficient for size sampling based on curvature.
  double m_max_dihed; //maximum dihedral angle between faces.
  double m_max_angle; //maximum angle deviation along an edge.
  //The (only) constructor, and it is private.
  SamplingEntryFactory() :
      m_sampler(NULL), m_new_entries(NULL), m_checked_edges(), m_bad_bridges(), m_checked_bridges(), m_test_faces(), m_used_faces(), m_dihed_samp(
	  true), m_norm_samp(true), m_size_samp(true), m_qual_samp(true), m_radius_edge_bound(
	  1.), m_radius_curv_ratio(1.), m_max_dihed(90.), m_max_angle(90.)
  {
  }

  //Some useful typedefs.
  typedef std::multimap<Vert*, FaceQueueEntry*> RingMap;
  typedef std::pair<RingMap::iterator, RingMap::iterator> RingRange;

public:

  typedef SurfaceSampler::VertWithTag VertWithTag;

  //Access to the single instance of this class.
  static SamplingEntryFactory*
  instance();

  //Functions to call to build sampling entries from a set of restricted faces.
  //When many vertices and many faces need to be checked. Typically at initialization
  //of the restricted triangulation.
  void
  build_entries(const SurfaceSampler* const sampler,
		const std::vector<VertWithTag>& sample_verts,
		const std::vector<SurfaceSampler::RestrictFace>& restrict_faces,
		std::vector<SamplingEntry*>* new_entries);

  //Typically, the one to call after an insertion. vertex MUST be appearing
  //in ALL the restricted faces. Failure to do so leads to a fatal error.
  void
  build_entries(const SurfaceSampler* const sampler,
		const VertWithTag& vertex_and_tag,
		const std::vector<SurfaceSampler::RestrictFace>& restrict_faces,
		std::vector<SamplingEntry*>* new_entries);

  //Checks the topological correctness of the neighborhood of vertex.
  void
  build_entries(const SurfaceSampler* const sampler,
		const VertWithTag& vertex_and_tag,
		std::vector<SamplingEntry*>* new_entries);

  virtual
  ~SamplingEntryFactory();

private:
  //Functions used to actually build the sampling entries.
  void
  build_missing_subseg(SubsegBridge* const bridge);

  void
  build_multiple_intersect(FaceQueueEntry* const fqe);

  void
  build_topo_violation(Vert* const vertex, SampleVertTag* vertex_tag);

  void
  build_size(Face* const face);

  void
  build_quality(FaceQueueEntry* const fqe, double radius_edge_ratio);

  //Initialization functions
  void
  reset_factory();

  void
  initialize_face(const SurfaceSampler::RestrictFace& restrict_face);

  //Insepction functions
  bool
  inspect_vertex(const VertWithTag& vertex_and_tag);
  bool
  inspect_apex_vert(Vert* const vertex, SampleVertTag* const vertex_tag);
  bool
  inspect_curv_vert(Vert* const vertex, SampleVertTag* const vertex_tag);
  bool
  inspect_surf_vert(Vert* const vertex, SampleVertTag* const vertex_tag);

  void
  inspect_face(FaceQueueEntry* const fqe);
  bool
  inspect_face_quality(FaceQueueEntry* const fqe);

  //Functions to verify topology.
  bool
  check_subsegs(Vert* const vertex);
  bool
  check_subsegs(Vert* const vertex, SampleVertTag* const vertex_tag);

  bool
  check_disk(Vert* const vertex, SampleVertTag* const vertex_tag,
	     SubsegBridge* beg_bridge = NULL, SubsegBridge* end_bridge = NULL);

  bool
  check_peripheral_edges(const VertWithTag& vert_and_tag);
  bool
  check_peripheral_edges(Vert* const vertex, SampleVertTag* const vertex_tag);

  bool
  valid_normal_deviation(Vert* const vert0, Vert* const vert1,
			 SampleVertTag* tag0 = NULL,
			 SampleVertTag* tag1 = NULL) const;

  bool
  valid_dihedral_angle(FaceQueueEntry* const face0, FaceQueueEntry* const face1,
		       Vert* const vert0, Vert* const vert1) const;

  bool
  valid_edge_length(Vert* const vert0, Vert* const vert1,
		    SampleVertTag* const tag0) const;

  //
  bool
  edge_checked(Vert* const vert0, Vert* const vert1);

  bool
  edge_in_ring(const Vert* const vert0, const Vert* const vert1,
	       const FaceQueueEntrySet& ring_faces) const;

  bool
  face_has_edge(const Face* const face, const Vert* const vert0,
		const Vert* const vert1) const;

  Face*
  largest_face_in_ring(const FaceQueueEntrySet& ring_faces) const;
  FaceQueueEntry*
  largest_unused_face_in_ring(const FaceQueueEntrySet& ring_faces) const;

  Vert*
  face_other_vert(Face* const face, Vert* const vert0, Vert* const vert1) const;

  //Functor to transform the ring faces of a vertex into a more usable form.
  struct RingMapper {

    const Vert* m_vertex;
    SamplingEntryFactory::RingMap* m_ring_map;

    RingMapper(const Vert* const vertex, RingMap* ring_map) :
	m_vertex(vertex), m_ring_map(ring_map)
    {
    }
    void
    operator()(FaceQueueEntry* const fqe);

  };

  struct RingBridgeMapper {

    const Vert* m_vertex;
    const SurfaceSampler* m_sampler;
    std::vector<SubsegBridge*>* m_bridges;

    RingBridgeMapper(const Vert* const vertex,
		     const SurfaceSampler* const sampler,
		     std::vector<SubsegBridge*>* bridges) :
	m_vertex(vertex), m_sampler(sampler), m_bridges(bridges)
    {
    }

    void
    operator()(FaceQueueEntry* const fqe);

  };

  struct BridgeNotEncroached {

    const Vert* const m_vertex;

    BridgeNotEncroached(const Vert* const vertex) :
	m_vertex(vertex)
    {
    }
    bool
    operator()(SubsegBridge* const bridge);

  };

  //Given the ring faces, finds the ring vertices that are
  //located on a boundary entity and returns them in the set (if any).
  struct BdryRingMapper {

    const Vert* m_vertex;
    std::set<Vert*>* m_bdry_ring_verts;

    BdryRingMapper(const Vert* const vertex, std::set<Vert*>* bdry_ring_verts) :
	m_vertex(vertex), m_bdry_ring_verts(bdry_ring_verts)
    {
    }
    void
    operator()(FaceQueueEntry* const fqe);

  };

  //Functor to find the triangle edge not connected to m_vertex
  struct PeriphEdgeID {

    const Vert* m_vertex;
    std::set<FactoryEdge>* m_edges;

    PeriphEdgeID(const Vert* vertex, std::set<FactoryEdge>* edges) :
	m_vertex(vertex), m_edges(edges)
    {
    }
    void
    operator()(FaceQueueEntry* const fqe);

  };

  //Functor to find the largest face among a group of faces.
  struct LargestFace {

    FaceQueueEntry*& m_fqe;
    const FaceQueueEntrySet* m_exclude_faces;
    double m_large_size;

    LargestFace(FaceQueueEntry*& fqe,
		const FaceQueueEntrySet* const exclude_faces = NULL) :
	m_fqe(fqe), m_exclude_faces(exclude_faces), m_large_size(-LARGE_DBL)
    {
    }
    ;
    void
    operator()(FaceQueueEntry* const fqe);

  };

  //Only used in debugging mode.

  bool
  valid_data_passed(const VertWithTag& vert_and_tag) const;
  bool
  valid_data_passed(const std::vector<VertWithTag>& sample_verts) const;
  bool
  valid_data_passed(
      const std::vector<SurfaceSampler::RestrictFace>& restrict_faces) const;

  bool
  test_faces_have_vertex(Vert* const vertex) const;

  //Just for the sake of defining them.
  SamplingEntryFactory(const SamplingEntryFactory&);
  SamplingEntryFactory&
  operator=(const SamplingEntryFactory&);

};

inline void
SamplingEntryFactory::build_missing_subseg(SubsegBridge* const bridge)
{

  m_new_entries->push_back(
      new SubsegSamplingEntry(m_sampler->m_surface, bridge,
			      bridge->get_parent()->get_length()));

}

inline void
SamplingEntryFactory::build_multiple_intersect(FaceQueueEntry* const fqe)
{

  if (m_used_faces.insert(fqe).second)
    m_new_entries->push_back(
	new FaceSamplingEntry(SamplingEntry::MULTIPLE, m_sampler->m_surface,
			      fqe, fqe->getFace()->calcSize()));

}

inline void
SamplingEntryFactory::build_topo_violation(Vert* const vertex,
					   SampleVertTag* vertex_tag)
{
  FaceQueueEntry* fqe = largest_unused_face_in_ring(
      vertex_tag->get_ring_faces());
  if (!fqe)
    return; //all possible faces are already used.

#ifndef NDEBUG
  bool inserted = m_used_faces.insert(fqe).second;
  assert(inserted);
#else
  m_used_faces.insert(fqe);
#endif

  m_new_entries->push_back(
      new FaceTopoSamplingEntry(m_sampler->m_surface, fqe, vertex, vertex_tag,
				fqe->getFace()->calcSize()));

}

inline void
SamplingEntryFactory::build_size(Face* const face)
{

  assert(face);
  FaceQueueEntry * fqe = new FaceQueueEntry(face);
#ifndef NDEBUG
  bool inserted = m_used_faces.insert(fqe).second;
  assert(inserted);
#else
  m_used_faces.insert(fqe);
#endif

  m_new_entries->push_back(
      new FaceSamplingEntry(SamplingEntry::SIZE, m_sampler->m_surface, fqe,
			    face->calcSize()));

}

inline void
SamplingEntryFactory::build_quality(FaceQueueEntry* const fqe,
				    double radius_edge_ratio)
{

  assert(fqe);
  assert(radius_edge_ratio > m_radius_edge_bound);
#ifndef NDEBUG
  bool inserted = m_used_faces.insert(fqe).second;
  assert(inserted);
#else
  m_used_faces.insert(fqe);
#endif

  m_new_entries->push_back(
      new FaceSamplingEntry(SamplingEntry::QUALITY, m_sampler->m_surface, fqe,
			    radius_edge_ratio));

}
namespace GRUMMP
{

  template<typename ReturnType, typename ClassType, typename ArgType>
    class extern_mem_fun_t : public std::unary_function<ArgType, ReturnType> {
    private:
      ReturnType
      (ClassType::*m_func_ptr)(ArgType);
      ClassType* m_ptr_obj;
    public:
      extern_mem_fun_t(ReturnType
      (ClassType::*ptr)(ArgType),
		       ClassType* obj) :
	  m_func_ptr(ptr), m_ptr_obj(obj)
      {
      }
      ReturnType
      operator()(ArgType arg)
      {
	return ((m_ptr_obj->*m_func_ptr)(arg));
      }
    };

  template<typename ReturnType, typename ClassType, typename ArgType>
    extern_mem_fun_t<ReturnType, ClassType, ArgType>
    extern_mem_fun(ReturnType
    (ClassType::*ptr_mfunc)(ArgType),
		   ClassType* type)
    {
      return (extern_mem_fun_t<ReturnType, ClassType, ArgType>(ptr_mfunc, type));
    }

}
#endif
