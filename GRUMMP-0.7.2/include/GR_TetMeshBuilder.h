#ifndef GR_TET_MESH_BUILDER_H
#define GR_TET_MESH_BUILDER_H

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_InsertionManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_RefinementManager3D.h"
#include <map>
#include <set>

class CubitVector;
class RefFace;
class Subseg;
class TriFace;
class Vert;
class VolMesh;

class RecoverEdge {

  Vert *m_v0, *m_v1;

  //Store extra info for insertion on curved boundaries.
  //The RecoverEdge stores either a Subseg or a TriFace and a RefFace:
  //Don't bother with a union on the first draft of this code.

  //If the edge is a subsegment (subcurve), store the Subseg.
  Subseg* m_subseg;

  //If the edge is on a surface (RefFace), store the TriFace and the RefFace
  TriFace* m_face;
  RefFace* m_surface;

public:

  RecoverEdge() :
      m_v0(NULL), m_v1(NULL), m_subseg(NULL), m_face(NULL), m_surface(
      NULL)
  {
  }
  RecoverEdge(const RecoverEdge& edge) :
      m_v0(edge.m_v0), m_v1(edge.m_v1), m_subseg(edge.m_subseg), m_face(
	  edge.m_face), m_surface(edge.m_surface)
  {
  }
  RecoverEdge&
  operator=(const RecoverEdge& edge)
  {
    if (&edge != this) {
      m_v0 = edge.m_v0;
      m_v1 = edge.m_v1;
      m_subseg = edge.m_subseg;
      m_face = edge.m_face;
      m_surface = edge.m_surface;
    }
    return *this;
  }

  RecoverEdge(Vert* const v0, Vert* const v1) :
      m_v0(std::min(v0, v1)), m_v1(std::max(v0, v1)), m_subseg(NULL), m_face(
      NULL), m_surface(NULL)
  {
    assert(m_v0 != m_v1);
    assert(m_v0 < m_v1);
  }

  //For RecoverEdge objects coming from a surface mesh, need to keep track
  //of the subsegsegment to know where the edge must be split
  void
  set_subseg(Subseg* const subseg)
  {
    m_subseg = subseg;
  }

  void
  set_face(TriFace* const face)
  {
    m_face = face;
  }
  void
  set_surface(RefFace* const surface)
  {
    m_surface = surface;
  }

  Vert*
  get_v0() const
  {
    return m_v0;
  }
  Vert*
  get_v1() const
  {
    return m_v1;
  }

  Subseg*
  get_subseg() const
  {
    return m_subseg;
  }

  TriFace*
  get_face() const
  {
    return m_face;
  }
  RefFace*
  get_surface() const
  {
    return m_surface;
  }

  bool
  operator<(const RecoverEdge& edge) const
  {
    assert(edge.valid());
    return ((m_v0 < edge.m_v0) || ((m_v0 == edge.m_v0) && (m_v1 < edge.m_v1)));
  }

  bool
  operator==(const RecoverEdge& edge) const
  {
    assert(edge.valid());
    return (edge.get_v0() == m_v0 && edge.get_v1() == m_v1);
  }

  bool
  operator!=(const RecoverEdge& edge) const
  {
    return !(*this == edge);
  }

  bool
  valid() const
  {
    if (!m_v0 || !m_v1)
      return false;
    return (m_v0 < m_v1);
  }

};

class RecoverFace {

  //The face to recover as defined by its three verts.
  Vert *m_v0, *m_v1, *m_v2;

  //Can store the face to recover and the surface associated
  TriFace* m_face;
  RefFace* m_surface;

public:

  RecoverFace() :
      m_v0(NULL), m_v1(NULL), m_v2(NULL), m_face(NULL), m_surface(NULL)
  {
  }
  RecoverFace(const RecoverFace& face) :
      m_v0(face.m_v0), m_v1(face.m_v1), m_v2(face.m_v2), m_face(face.m_face), m_surface(
	  face.m_surface)
  {
  }
  RecoverFace&
  operator=(const RecoverFace& face)
  {
    if (&face != this) {
      m_v0 = face.m_v0;
      m_v1 = face.m_v1;
      m_v2 = face.m_v2;
      m_face = face.m_face;
      m_surface = face.m_surface;
    }
    return *this;
  }

  RecoverFace(Vert* const v0, Vert* const v1, Vert* const v2) :
      m_v0(MIN3(v0, v1, v2)), m_v1(MID3(v0, v1, v2)), m_v2(MAX3(v0, v1, v2)), m_face(
      NULL), m_surface(NULL)
  {
    assert(m_v0 != m_v1 && m_v0 != m_v2 && m_v1 != m_v2);
    assert(m_v0 < m_v1 && m_v1 < m_v2);
  }

  void
  set_face(TriFace* const face)
  {
    assert(face);
    m_face = face;
  }
  void
  set_surface(RefFace* const surface)
  {
    assert(surface);
    m_surface = surface;
  }

  Vert*
  get_v0() const
  {
    return m_v0;
  }
  Vert*
  get_v1() const
  {
    return m_v1;
  }
  Vert*
  get_v2() const
  {
    return m_v2;
  }

  TriFace*
  get_face() const
  {
    return m_face;
  }
  RefFace*
  get_surface() const
  {
    return m_surface;
  }

  bool
  operator<(const RecoverFace& face) const
  {
    assert(face.valid());
    return (((m_v0 < face.m_v0)) || ((m_v0 == face.m_v0) && (m_v1 < face.m_v1))
	|| ((m_v0 == face.m_v0) && (m_v1 == face.m_v1) && (m_v2 < face.m_v2)));
  }

  bool
  operator==(const RecoverFace& face) const
  {
    assert(face.valid());
    return ((m_v0 == face.m_v0) && (m_v1 == face.m_v1) && (m_v2 == face.m_v2));
  }

  bool
  operator!=(const RecoverFace& face) const
  {
    return !(*this == face);
  }

  bool
  valid() const
  {
    if (!m_v0 || !m_v1)
      return false;
    return (m_v0 < m_v1 && m_v1 < m_v2);
  }

};

class TetMeshBuilder {

  VolMesh* m_mesh_init;
  //The surface mesh (in sampler, as defined by restricted faces).
  SurfMeshBuilder* m_surf_mesh;

  // During surface recovery, insert-and-swap must be used; this
  // insertion manager handles that.  The swapping stuff is required to
  // get Delaunay swapping to occur (as opposed to nothing).
  GRUMMP::DelaunaySwapDecider3D* m_pDel3D;
  GRUMMP::SwapManager3D* m_pSM3D;
  GRUMMP::SwappingInserter3D* m_pSI3D;
  GRUMMP::WatsonInserter * m_WI;
  //A vertex dictionary to maintain correspondance from surface to tet mesh.
  //Changed from a hash_map to a std::map for a more deterministic behavior.
  //Switching the typedef is enough to move back to the hash_map.

  typedef std::map<Vert*, Vert*> VertexDict;

  VertexDict m_sample_vert_dict;
  VertexDict m_subseg_vert_dict;

  //Typedef used in the class to store "special" vertex pairs.
  typedef std::vector<SurfMeshBuilder::SubsegSamplePair> VertVec;

public:

  typedef std::vector<std::pair<Vert*, Vert*> > SubsegVec;

  //Constructor and destructor.
  //Optionally, the constructor can return a set containing
  //a copy of the subsegments, as contained in the SubsegManager of
  //m_surf_mesh, BUT with the vertex pointers of the tetrahedral mesh.
  //NOTE: These Subseg pointers are instantiated here, the memory must
  //be freed when the go out of scope.
  TetMeshBuilder(SurfMeshBuilder* const surf_mesh);

  ~TetMeshBuilder();

  // Deletes the initial mesh created internally, when needed.  Sometimes that
  // mesh pointer is also used externally, so it can't just be destroyed at
  // TMB destruction time.
  void
  destroyInitMesh();

  // builds the mesh, can return false if unable to recover watertight exterior
  bool
  buildMesh(std::set<Subseg*>* subsegs = NULL);
  //Gets the tetrahedral mesh as contained in m_mesh_init.
  VolMesh*
  get_mesh() const;

  /*   //Outputs the content of the set in medit format. */
  /*   void output_faces(const char* const filename, */
  /* 		    const std::set<Face*>& faces) const; */

private:

  //Default constructor, copy constructor and operator= made private.
  //No use for them at the moment.
  TetMeshBuilder();
  TetMeshBuilder(const TetMeshBuilder&);
  TetMeshBuilder&
  operator=(const TetMeshBuilder&);

  //This instantiates m_mesh_init. It creates a tet mesh of
  //a large bounding box containing all the sampling vertices.
  void
  init_mesh(RefVolume * const refVolume);

  //Inserts the sample vertices in the tet mesh.
  //This will also initialize the two vert dictionaries.
  void
  insert_sample_verts();

  //After building the tet mesh from the sampling data, attempts
  //to recover the boundary edges and faces.
  bool
  recover_boundary();

  //Recovers all edges. Called from recover_boundary();
  //Returns true if the entire boundary is constrained (edges and faces)
  bool
  recover_all_edges(std::set<RecoverEdge>& bdry_edges_to_recover,
		    std::set<RecoverEdge>& face_edges_to_recover);

  //Attempts to recover an edge by swapping.
  int
  recover_edge_flip(const RecoverEdge& edge_to_recover);

  //Attempts to recover an edge by splitting.
  int
  recover_edge_split(const RecoverEdge& edge_to_recover);

  //Attempts to recover a face by swapping.
  int
  recover_face_flip(const RecoverFace& face_to_recover);

  //Attemps to recover a face by splitting.
  int
  recover_face_split(const RecoverFace& face_to_recover);

  //Inserts a copy of the verts contained in new_verts
  //(obtained from sampling algorithm) into the tetrahedral mesh
  void
  insert_new_verts(const VertVec& new_verts);

  //Inserts new_vertex in the tet mesh,
  //then swaps to retain Delaunay property.
  bool
  insert_vert_and_swap(Vert* const vertex);

  //Verifies if all edges and faces defined in the surface mesh
  //are present in the Delaunay tetrahedralization.
  bool
  is_bdry_constrained(std::set<RecoverEdge>& bdry_edges_to_recover,
		      std::set<RecoverEdge>& face_edges_to_recover,
		      std::set<RecoverFace>& faces_to_recover) const;

  //Deletes the cells that are outside the meshed domain.
  bool
  clean_exterior();

  bool
  mesh_is_const_delaunay() const;

  void
  update_vert_dict(const std::map<Vert*, Vert*>& verts_old_to_new);

  //Returns the normal on the surface at a location closest to point.
  CubitVector
  surf_normal_at_point(RefFace* const surface, const CubitVector& point,
		       const double tolerance = 1.e-3) const;

  void
  set_surf_compare_tol(RefFace* const surface, const double tolerance) const;

  Cell*
  getSeedCellFromVerts(const Vert* const vertex, const Vert* const vert_guess1 =
  NULL,
		       const Vert* const vert_guess2 = NULL,
		       const Vert* const vert_guess3 = NULL);

  void
  print_markings() const;

  //Creates a copy of the subsegments contained in the SubsegManager
  //object of m_surf_mesh. Replaces the vertices of the SubsegManager
  //by those of the tetrahedral mesh.
  //NOTE: The Subseg pointers are allocated here, the memory must be
  //freed when they are no longer used.
  void
  get_subseg_copy(std::set<Subseg*>* subseg_copy) const;

  //A functor to help construct the mesh bounding box
  //from a set of vertices.
  struct ComputeBoundingBox {
    CubitVector *m_mini, *m_maxi;
    ComputeBoundingBox(CubitVector* const mini, CubitVector* const maxi) :
	m_mini(mini), m_maxi(maxi)
    {
    }
    void
    operator()(const Vert* const vertex);
  };

};

#endif
