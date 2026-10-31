#ifndef GR_SURFACESAMPLER_H
#define GR_SURFACESAMPLER_H

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include "GR_Face.h"
#include "GR_Subseg.h"
#include "GR_SurfIntersect.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_InsertionManager.h"
#include <map>
#include <utility>

//Defined elsewhere
class Cell;
class CubitVector;
class RefFace;
class SamplingEntry;
class TriFace;
class Vert;
class VoronoiEdge;

//Defined in this header file
class SampleVertTag;
class RestrictedFaceTag;
class SamplingEntryFactory;

class SurfaceSampler : public GRUMMP::Observer {

public:

  friend class SurfMeshBuilder;
  friend class SamplingEntryFactory;

  typedef std::map<Vert*, Vert*> VertMap;
  typedef std::map<Face*, Face*> FaceMap;
  typedef std::map<Cell*, Cell*> CellMap;

  //Might want to relace these with TR1::unordered_map
  typedef std::map<Vert*, SampleVertTag*> SampleVertDict;
  typedef std::map<FaceQueueEntry*, RestrictedFaceTag*, FaceQueueEntryComp> RestrictedFaceDict;

  typedef std::pair<FaceQueueEntry*, SurfIntersect::IntersectSet> RestrictFace;
  typedef std::pair<Vert*, SampleVertTag*> VertWithTag;

private:

  //The modeler type we are using.
  SurfMeshBuilder::ModelerType m_modeler_type;

  VolMesh* m_mesh;
  GRUMMP::WatsonInserter* m_WI;
  RefFace* m_surface;
  SurfIntersect* m_surf_intersect;

  SampleVertDict m_sample_verts;

  RestrictedFaceDict m_restricted_faces;

  //-- m_mesh: Contains the mesh used to sample the surface.
  //-- m_surface: Contains the geometric surface we want to sample.
  //-- m_surf_intersect: Used to compute intersections with m_surface.
  //-- m_entry_builder: Used to construct sampling entries from restricted triangles.

  //-- m_sample_verts: Dictionary associating every sample vert in m_mesh with its tag.
  //-- m_restricted_faces: Dictionnary associating a restricted face with tag.

  //Sanity check to make sure everything is initialized.
  bool m_restrict_initialized;

  //Observable stuff
  std::vector<Cell*> m_new_cells;
  // I really dislike this, it exists to ensure that I clean up facequeueentries correctly.
  std::set<FaceQueueEntry*> m_temporary_fqes;
public:

  //(The only public) constructor.
  //RefFace* surface is the surface to sample
  SurfaceSampler(RefFace* const surface,
		 const SurfMeshBuilder::ModelerType modeler_type);

  //Destructor.
  ~SurfaceSampler();

public:

  void
  print_surface_sampler();

  //A bunch of interfaces are defined to insert a new sample points, each
  //allow to pass some sort of guess to expedite Watson Hull computation :
  //** All return a pointer to the newly inserted vertex.
  //** Note that the subseg bridge pointers are not updated in the vertex tags.

  //- Inserts a new vertex at vertex's location.
  //- Optionally, specify verts to find a seed guess.
  //- Return a pointer to the new vertex. Optionally BUILDS (instantiates)
  //  the new sampling entries computed from the restricted faces.
  Vert*
  insert_sample_point(const Vert* const vertex, const Vert* const guess1 = NULL,
		      const Vert* const guess2 = NULL,
		      const Vert* const guess3 = NULL,
		      std::vector<RestrictFace>* new_restrict_faces = NULL);

  //- Inserts a new_vertex at vertex's location.
  //- Uses the cells adacent to face as likely seed guesses.
  //- Return a pointer to the new vertex. Optionally BUILDS (instantiates)
  //  the new sampling entries computed from the restricted faces.
  Vert*
  insert_sample_point(const Vert* const vertex, const Face* const face,
		      int attemptNum,
		      std::vector<RestrictFace>* new_restrict_faces = NULL);

  //- Inserts a new_vertex at vertex's location.
  //- Uses the contents of vert_guesses to build a list of potential seed cells.
  //- Return a pointer to the new vertex. Optionally BUILDS (instantiates)
  //  the new sampling entries computed from the restricted faces.
  Vert*
  insert_sample_point(const Vert* const vertex,
		      const std::set<const Vert*>& vert_guesses,
		      std::vector<RestrictFace>* new_restrict_faces = NULL);

  //- Inserts a new_vertex at vertex's location.
  //- Uses the contents of seed_guesses as potential seed cells.
  //- Return a pointer to the new vertex. Optionally BUILDS (instantiates)
  //  the new sampling entries computed from the restricted faces.
  //- If bail_on_bad_guesses == true, will return NULL if no valid
  //  Watson seed cells are found in seed_guesses.
  Vert*
  insert_sample_point(const Vert* const vertex,
		      const std::set<Cell*>& seed_guesses,
		      const bool bail_on_bad_guesses,
		      std::vector<RestrictFace>* new_restrict_faces = NULL);
  // use this if we don't have any idea where the seed is
  Vert*
  insert_arbitrary_point(const Vert* const vertex, Cell * guess,
			 std::vector<RestrictFace>* new_restrict_faces);
  //- Populates the vector with the sample vertices contained in m_sample_verts.
  void
  get_sample_verts(std::vector<Vert*>& sample_verts) const;

  //- Populates the vector with the faces part of the restricted
  //  Delaunay triangulation (identified with qRestricted == 1 in the Face class).
  //** Note that this is not necessarily all the faces in m_restricted_faces.
  //   Some of them might have been deleted without being removed from the dict.
  void
  get_restricted_faces(std::vector<FaceQueueEntry*>& restricted_faces) const;

  // this takes an array of FaceQueueEntries with valid faces but not valid Entries
  // (different vertices from the face) and resets the default verts
  void
  reset_restricted_faces() const;
  //- Returns the SampleVertTag of vertex.
  SampleVertTag*
  get_vert_tag(Vert* const vertex) const;

  //- Returns the RestrictedFaceTag of face.
  RestrictedFaceTag*
  get_face_tag(FaceQueueEntry* fqe) const;

  bool
  vert_has_tag(Vert* const vertex) const
  {
    return m_sample_verts.count(vertex) == 1;
  }
  bool
  face_has_tag(FaceQueueEntry* const fqe) const
  {
    return m_restricted_faces.count(fqe) == 1;
  }

private:

  //- The default constructor made private. No use for it now.
  SurfaceSampler() :
      m_modeler_type(SurfMeshBuilder::UNKNOWN), m_mesh(NULL), m_WI(NULL), m_surface(
      NULL), m_surf_intersect(NULL), m_sample_verts(), m_restricted_faces(), m_restrict_initialized(
	  false)
  {
    assert(0);
  }

  //- The copy constructor made private. No use for it now.
  SurfaceSampler(const SurfaceSampler&) :
      m_modeler_type(SurfMeshBuilder::UNKNOWN), m_mesh(NULL), m_WI(NULL), m_surface(
      NULL), m_surf_intersect(NULL), m_sample_verts(), m_restricted_faces(), m_restrict_initialized(
	  false)
  {
    assert(0);
  }

  //- The operator= made private. No use for it now.
  SurfaceSampler&
  operator=(const SurfaceSampler&)
  {
    assert(0);
    return *this;
  }

  //The following methods are for sampler initialization:

  //- Instantiates the VolMesh contained in m_mesh.
  //- Builds a tetrahedral mesh inside a rectangular prism
  //  big enough to contain containing m_ref_face.
  void
  initialize_mesh();

  //- Computes the restricted triangulation data.
  //- Builds m_restricted_faces and tags the faces appropriately.
  //- Returns the restricted faces and the sample verts (with tags)
  //  in the vectors.
  //  Called from SurfMeshBuilder.
  void
  init_restricted_delaunay(std::vector<VertWithTag>& sample_verts,
			   std::vector<RestrictFace>& restricted_faces);

  //- Determines if the face associated to voronoi_edge is restricted.
  //  A face is said to be restricted if its Voronoi edge intersects the surface.
  //- IF intersection(s) are found, returns the restricted faces and their
  //  interesction(s) coords in the vector.
  void
  compute_restricted_data(const VoronoiEdge& voronoi_edge,
			  std::vector<RestrictFace>& restricted_faces);

  // It is possible, althought very rare, that a face be found to
  // be restricted Delaunay while being attached to the tetrahedralization's
  // bounding box. In this case, the culprits must be split until
  // the situation is corrected. The restrict_faces vector gets updates accordingly.
  // This is never called TODO: Remove This function
//  void split_box_attached(std::vector<RestrictFace>& restrict_faces);

  //- IF no VoronoiEdge intersects the surface, it needs to be seeded
  //  arbitrarily with vertices by calling this functions.
  void
  seed_surface(std::vector<RestrictFace>& restricted_faces);

public:

  //-Returns the normal at a point on m_surface
  // (projects the point on the surface first).
  //**Note: The tolerance typically does not need
  //  to be set as low as for a project-to-surface operation.
  CubitVector
  normal_at(const Vert* const vertex, const double tolerance = 1.e-3) const;

  //-Finds the coord on the surface closest to point (point_on_surf).
  // Also computes the normal at the closest point (normal_at_point).
  void
  project_to_surface(const CubitVector& point, CubitVector& point_on_surf,
		     CubitVector& normal_at_point, const double tolerance =
			 1.e-14) const;

  //-Obtains the 1-star of vertex in the restricted triangulation.
  //**Note: this 1-star does not necessarily form a cycle around the vertex,
  //  either because of a topology infraction or because vertex is on the boundary.
//  void restricted_neighborhood(const Vert* const vertex,
//			       std::vector<Face*>* const neigh_faces,
//			       std::vector<Vert*>* const neigh_verts = NULL) const;

  bool
  restrict_delaunay_initialized() const
  {
    return m_restrict_initialized;
  }

  //The following three methods mostly used for debugging purposes.

  //-Outputs the restricted Delaunay triangulation to file in medit format.
  void
  output_restricted_delaunay(const char* const filename) const;

  //-Outputs the faces_to_print in medit format
  void
  output_faces(const char* const filename,
	       const std::set<Face*>& faces_to_print) const;
  void
  output_faces(const char* const filename,
	       const std::vector<Face*>& faces_to_print) const;

  // Observer details
  void
  receiveCreatedCells(std::vector<Cell*>& createdCells);
private:

  //-After an insertion in m_mesh, the new cells are passed to this
  // function. It builds the Voronoi diagram accordingly, populating
  // new_edges with the newly created Voronoi edges.
  void
  build_voronoi_edges(const std::vector<Cell*>& new_cells,
		      std::vector<VoronoiEdge>& new_edges);

  //- Attaches a subseg bridge to a vertex (in its tag).
  void
  attach_bridge(Vert* const vertex, SubsegBridge* const bridge);

  //- Removes a subseg bridge from the vertex (in its tag).
  void
  remove_bridge(Vert* const vertex, SubsegBridge* const bridge);

  //- Removes the tags associated to deleted
  //  sample verts from m_sample_verts.
  void
  clean_vert_tags();

  //- Removes the tags associated to deleted
  //  restricted faces from m_restricted_faces.
  void
  clean_face_tags();

  //Garbage collection.
  void
  purge_sampler();

};

//A Voronoi edge is a line segment connecting the circumcenters
//of two adjacent Delaunay tetrahedra. As such, we store the coords
//of these circumcenters and the face shared by the cells.

class VoronoiEdge {

  Face* m_face;
  CubitVector m_beg, m_end;

public:

  //Slighlty faster if we know the face.
  VoronoiEdge(Face* const face, Cell* const cell1, Cell* const cell2);

  //If we don't know the face.
  VoronoiEdge(Cell* const cell1, Cell* const cell2);

  //If we know the face and the coords.
  VoronoiEdge(Face* const face, const CubitVector& beg_coord,
	      const CubitVector& end_coord) :
      m_face(face), m_beg(beg_coord), m_end(end_coord)
  {
  }

  //Copy constructor.
  VoronoiEdge(const VoronoiEdge& VE) :
      m_face(VE.m_face), m_beg(VE.m_beg), m_end(VE.m_end)
  {
  }

  VoronoiEdge&
  operator=(const VoronoiEdge& VE)
  {
    if (&VE != this) {
      m_face = VE.m_face;
      m_beg = VE.m_beg;
      m_end = VE.m_end;
    }
    return *this;
  }

  ~VoronoiEdge()
  {
  }

  //Returns face associated with this edge.
  Face*
  get_face() const
  {
    return m_face;
  }

  //Returns Voronoi vertex coords
  const CubitVector&
  beg_pt() const
  {
    return m_beg;
  }
  const CubitVector&
  end_pt() const
  {
    return m_end;
  }

  //Computes the Voronoi vertex coord associated to cell. Cell type
  //must either be eTet (coord = cell's circumcenter) or eTriBFace
  //(coord = triangle's circumcenter).
  static void
  voronoi_vert_coord(Cell* const cell, CubitVector& coord);

  void
  print_voronoi() const
  {
    logMessage(1, "%f %f %f\n", m_beg.x(), m_beg.y(), m_beg.z());
    logMessage(1, "%f %f %f\n", m_end.x(), m_end.y(), m_end.z());
  }

private:

  VoronoiEdge();
};

class SampleVertTag {

  //Info stored in tag:
  const Vert* m_vertex;  //The vertex.
  CubitVector* m_normal;  //Normal to   the surface at m_vertex.
  double m_curv;    //Surface max principal curvature at m_vertex.

  std::set<SubsegBridge*> m_bridges;    //Subseg bridges attached to m_vertex.
  FaceQueueEntrySet m_ring_faces; //The 1-ring of restricted faces.

  //Private copy constructor and operator=.
  SampleVertTag(const SampleVertTag&) :
      m_vertex(NULL), m_normal(NULL), m_curv(-LARGE_DBL), m_bridges(), m_ring_faces()
  {
    assert(0);
  }
  SampleVertTag&
  operator=(const SampleVertTag&)
  {
    assert(0);
    return *this;
  }

public:

  SampleVertTag(const Vert* const vertex) :
      m_vertex(vertex), m_normal(NULL), m_curv(-LARGE_DBL), m_bridges(), m_ring_faces()
  {
  }

  ~SampleVertTag()
  {
    if (m_normal)
      delete m_normal;
  }

  //Stuff to deal with the vertex.
  void
  set_vertex(Vert* const new_vertex)
  {
    m_vertex = new_vertex;
  }
  const Vert*
  get_vertex() const
  {
    return m_vertex;
  }

  //Stuff to deal with normal info.
  void
  set_normal(const CubitVector& normal)
  {
    if (m_normal)
      delete m_normal;
    m_normal = new CubitVector(normal);
  }
  const CubitVector*
  get_normal() const
  {
    return m_normal;
  }
  bool
  normal_is_set() const
  {
    return (m_normal != NULL);
  }

  //Stuff to deal with the curvature info.
  void
  set_curvature(const double curvature)
  {
    m_curv = curvature;
  }
  double
  get_curvature()
  {
    return m_curv;
  }

  //Stuff to deal with bridges.
  void
  add_bridge(SubsegBridge* const bridge)
  {
    assert(bridge->has_vert(const_cast<Vert*>(m_vertex)));
    m_bridges.insert(bridge);
  }

  void
  remove_bridge(SubsegBridge* const bridge)
  {
    m_bridges.erase(bridge);
  }

  const std::set<SubsegBridge*>&
  get_bridges()
  {
    assert(bridges_valid());
    return m_bridges;
  }

  void
  get_bridges(std::vector<SubsegBridge*>& bridges) const
  {
    assert(bridges_valid());
    std::copy(m_bridges.begin(), m_bridges.end(), std::back_inserter(bridges));
  }

  int
  num_bridges() const
  {
    assert(bridges_valid());
    return static_cast<int>(m_bridges.size());
  }

  bool
  has_bridge(SubsegBridge* const bridge) const
  {
    assert(bridges_valid());
    return (m_bridges.count(bridge) == 1);
  }

  bool
  bridges_valid() const
  {
    //Just checks if the set contains a bridge attached to a deleted subseg. Should not!
    return (std::find_if(m_bridges.begin(), m_bridges.end(),
			 std::mem_fun(&SubsegBridge::parent_deleted))
	== m_bridges.end());
  }

  //Stuff to deal with ring faces.
  FaceQueueEntrySet&
  get_ring_faces()
  {
    return m_ring_faces;
  }

  void
  add_ring_face(FaceQueueEntry* const fqe)
  {
    m_ring_faces.insert(fqe);
  }

  void
  remove_ring_face(FaceQueueEntry* const fqe)
  {
    m_ring_faces.erase(fqe);
  }

  bool
  has_ring_face(FaceQueueEntry* const fqe) const
  {
    return m_ring_faces.count(fqe) == 1;
  }

  int
  num_ring_faces() const
  {
    return std::count_if(m_ring_faces.begin(), m_ring_faces.end(),
			 std::mem_fun(&FaceQueueEntry::isNotDeleted));
  }

  void
  replace_ring_face(FaceQueueEntry* const old_fqe,
		    FaceQueueEntry* const new_fqe)
  {
    FaceQueueEntrySet::iterator itf = m_ring_faces.find(old_fqe);
    if (itf == m_ring_faces.end())
      return;
    m_ring_faces.erase(itf);
    m_ring_faces.insert(new_fqe);
  }

  void
  replace_all_ring_faces(const FaceQueueEntrySet& new_ring_faces)
  {
    m_ring_faces = new_ring_faces;
  }

  //Clears the set of ring faces.
  void
  clear_ring()
  {
    m_ring_faces.clear();
  }

  void
  clean_ring()
  { //Removed deleted faces from ring.
    FaceQueueEntrySet dummy;
    std::remove_copy_if(m_ring_faces.begin(), m_ring_faces.end(),
			std::inserter(dummy, dummy.begin()),
			std::mem_fun(&FaceQueueEntry::isDeleted));
    m_ring_faces.swap(dummy);
  }

  void
  clean_tag()
  {
    clean_ring();
    assert(m_vertex && m_vertex->isValid());
    assert(!m_vertex->isDeleted());
    assert(bridges_valid());
    assert(
	std::find_if(m_ring_faces.begin(), m_ring_faces.end(),
		     std::mem_fun(&FaceQueueEntry::isDeleted))
	    == m_ring_faces.end());
  }

};

class RestrictedFaceTag {

  //The tagged face.
  FaceQueueEntry* m_fqe;
  //The normal to the face.
  CubitVector m_face_normal;
  //Intersection data. Stores only one intersection per face. In cases of
  //multi intersection, must decide which one to keep (done elsewhere).
  SurfIntersect::Intersection m_intersect;

  //Private constructor, can only be called by friend class.
  RestrictedFaceTag(FaceQueueEntry* const fqe) :
      m_fqe(fqe), m_face_normal(), m_intersect()
  {
    assert(m_fqe && m_fqe->isValid());
    assert(!m_fqe->isDeleted());
    assert(m_fqe->getNumVerts() == 3);
    assert(m_fqe->getFace()->isRestrictedDelaunay());
    compute_face_normal();
  }

public:

  friend class SurfaceSampler;

  ~RestrictedFaceTag()
  {
    reset_tag();
  }

  void
  reset_tag()
  {
    m_intersect.reset_data();
  }

  void
  set_face(Face* const face)
  {
    m_fqe->setFace(face);
  }

  //Sets intersection data.
  void
  set_intersect(const SurfIntersect::Intersection& intersect)
  {
    m_intersect = intersect;
  }

  //Accessor functions
  FaceQueueEntry*
  get_fqe() const
  {
    return m_fqe;
  }
  const Face*
  get_face() const
  {
    return m_fqe->getFace();
  }
  const CubitVector&
  get_face_normal() const
  {
    return m_face_normal;
  }
  const SurfIntersect::Intersection&
  get_intersection() const
  {
    return m_intersect;
  }

  void
  clean_tag()
  {
    assert(m_fqe && m_fqe->getFace()->isValid());
    assert(!m_fqe->isDeleted());
    assert(m_fqe->getNumVerts() == 3);
    assert(m_fqe->getFace()->isRestrictedDelaunay());
  }

private:

  void
  compute_face_normal()
  {
    assert(m_fqe);
    double normal[3];
    ::calcNormal(m_fqe->getVert(0)->getCoords(), m_fqe->getVert(1)->getCoords(),
		 m_fqe->getVert(2)->getCoords(), normal);
    m_face_normal.set(normal);
  }

  RestrictedFaceTag(const RestrictedFaceTag&);
  RestrictedFaceTag&
  operator=(const RestrictedFaceTag&);

};

#endif
