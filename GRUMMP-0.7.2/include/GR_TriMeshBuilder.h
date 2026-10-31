#ifndef GR_TRIMESHBUILDER_H
#define GR_TRIMESHBUILDER_H

#include "GR_config.h"
#include "GR_VertexPair.h"
#include "GR_InsertionManager.h"

#include <deque>
#include <map>
#include <set>
#include <vector>

class CubitBox;
class CubitVector;
class GRCurve;
class GRPoint;
class GRGeom2D;
template<class X>
  class RTree;

// TODO This class only ever does one thing: it builds a Delaunay mesh for
// an input geometry, including discretizing that geometry.  There's no reason
// not to convert this to a function like
//
// buildConstrainedDelaunayTriangulation(Mesh2D& mesh, GRGeom2D& geom, const
//                                       CurveSamplingType sampleType);
//
// Member data and member functions are all private and can be file static (or
// part of a struct that's file static).
//
// This is a stylistic thing, not functional, so it's low priority.

class TriMeshBuilder {

  //The geometric domain's boundary definition.
  const GRGeom2D* m_geometry;

  //The mesh to be generated.
  Mesh2D* m_mesh;

  // A Watson insertion object to handle the actual insertion of vertices.
  GRUMMP::WatsonInserter m_WI;

  //The spatial indexing tree for verts.
  RTree<Vert*>* m_vert_tree;

  //The spatial indexing tree for edges.
  RTree<VertexPair*>* m_edge_tree;

  //A map between apex points and their associated vertex.
  std::map<GRPoint*, Vert*> m_apex_vert_map;

  //The boundary edges that needs to be recovered.
  std::set<VertexPair*> m_boundary_edges;

public:

  enum CurveSamplingType {
    TVT_BASED, UNKNOWN
  };

  //Public constructor
  TriMeshBuilder(const GRGeom2D* const geometry_to_mesh,
		 Mesh2D* const tri_mesh_to_build);
  //Destructor
  ~TriMeshBuilder();

  //The interface through which a constrained Delaunay triangulation is build.
  bool
  buildConstrainedDelaunayTriangulation(const CurveSamplingType sample_type =
      TVT_BASED);
  bool
  buildTerribleConstrainedDelaunayTriangulation(
      const CurveSamplingType sample_type = TVT_BASED);

private:

  //Populates vertices with the UNDELETED vertices in mesh. 
  void
  get_vertices(std::vector<Vert*>& vertices) const;
  void
  get_vertices(std::deque<Vert*>& vertices) const;

  //Populates boundary_edges with the contents of m_boundary_edges.
  void
  get_boundary_edges(std::vector<VertexPair*>& boundary_edges) const;

  //Returns a bounding box of the geometric model
  void
  findBoundingBox(double& xMin, double& yMin, double& xMax, double& yMax) const;

  //These are not used, so make private.
  TriMeshBuilder();
  TriMeshBuilder(const TriMeshBuilder&);
  TriMeshBuilder&
  operator=(const TriMeshBuilder&);

  struct EdgeSplitData {

    Vert* new_vertex;
    VertexPair* new_edge1;
    VertexPair* new_edge2;

  };

  //Initializes the vertices whose parent entity is a point.
  void
  init_apex_verts();

  //Initializes the vertices whose parent entity is a curve.
  void
  init_curve_verts(const CurveSamplingType sample_type);

  //Samples a curve according to sample_type.
  void
  sample_curve(GRCurve* const curve, CurveSamplingType sample_type);

  //Samples a curve according to the total variation of the tangent subtended to the curve.
  void
  sample_curve_based_on_tvt(GRCurve* const curve,
			    double max_tvt_allowed = M_PI / 2.);

  //Prepares the neighborhood of small angles so boundary edges sharing the 
  //apex of a small angle do not encroach on each other.
  //This is done by inserting vertices on concentric shells.
  void
  prepare_small_angles();

  //Delaunayizes the boundary edges.
  void
  make_boundary_constrained_delaunay();

  //Called from make_boundary_constrained_delaunay.
  //Checks if a vertex is visible from an edge.
  bool
  vert_inside_or_on_edge_ball(Vert* const vertex, VertexPair* const edge) const;

  //Called from make_boundary_constrained_delaunay.
  //Checks if a vertex is visible from an edge.
  bool
  vert_is_visible_from_edge(Vert* const vertex, VertexPair* const edge) const;

  //Splits edge_to_split at split_param along the edge's parent curve.
  //shell_radius used iif split_type == VertexPair::INIT_SHEL,
  EdgeSplitData
  split_boundary_edge(VertexPair* const edge_to_split,
		      const VertexPair::SplitType split_type =
			  VertexPair::MID_TVT,
		      const double shell_radius = -LARGE_DBL);

  //The same as above, but does not modify m_boundary_edges. However, the vertex
  //created by the split is inserted on the curve. Useful for sampling curves.
  void
  split_boundary_edge_at_mid_tvt(VertexPair& edge_to_split,
				 VertexPair& new_edge_1,
				 VertexPair& new_edge_2);

  //Split all the edges contained in small_angle_complex on a 
  //concentric shell centered at apex_vert.
  void
  small_angle_split(Vert* const apex_vertex,
		    const std::set<VertexPair*>& small_angle_complex);

  //Functions to operate on m_mesh.

  // TODO:  The two variants of insert_vertex_in_mesh have weirdly different
  // signatures.  Both should return the vertex.  Getting down to only one of
  // these would be nice.  At least they're private, so they don't provide a
  // bad API to users.

  //Inserts vertex in m_mesh using Watson insertion.
  //The coordinates of vertex must be properly initialized.
  void
  insert_vertex_in_mesh(Vert* const vertex,
			const std::set<Cell*>* const seed_guesses = NULL);

  //Inserts a vertex at insert_location in m_mesh.
  Vert*
  insert_vertex_in_mesh(const CubitVector& insert_location,
			const std::set<Cell*>* const seed_guesses = NULL);

  //Inserts a vertex having point as parent entity.
  Vert*
  insert_apex_vertex(GRPoint* const point);

  //Inserts a vertex at insert_location (located on a curve).
  //The first version is to be used when splitting an edge. It is able
  //to find a guess seed cell for Watson insertion, therefore insertion is faster.
  Vert*
  insert_curve_vertex(VertexPair* const edge_getting_split,
		      const CubitVector& insert_location);

  //Finds guesses for a seed cell to be used in Watson insertion
  //when splitting edge_to_split.
  void
  find_seed_guesses(VertexPair* const edge_to_split,
		    std::set<Cell*>& guesses) const;

  //Purge deleted boundary edges and adjust data structures.
  void
  purge_boundary_edges();

  //Purge m_mesh of deleted entities and adjust data structures.
  void
  purge_mesh();

  void
  recover_boundary_flip_only();

  void
  clean_up_domain();

public:

  //Outputs in medit format

  void
  output_boundary_edges_to_file(const char* const filename) const;

  /// TODO This function should become a separte Mesh2D output function, not a member func here.
  void
  output_triangulation_to_file(const char* const filename) const;

  void
  output_cells_to_file(const char* const filename,
		       const std::set<Cell*>& cells) const;
};

#endif
