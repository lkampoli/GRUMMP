#ifndef GR_INITVOLMESH_H
#define GR_INITVOLMESH_H

#include "GR_config.h"
#include "KDDTree.hpp"
#include <map>
#include <utility>
#include <vector>

class Cell;
class CubitBox;
class RefEdge;
class RefFace;
class RefVertex;
class SamplingQueue;
class Subseg;
class Vert;
class VolMesh;
class Voronoi;
class VoronoiEdge;

//This class is used to initialize volume meshes.
//It actually creates the vertices located on curves and at
//the points where these curves meet. At the moment
//it is set up to initialize constrained Delaunay tetrahedral
//meshes, however the class could probably be extended to other
//mesh types if need be.

class InitVolMesh {

  //The mesh
  VolMesh* mesh;

  //Voronoi data (only edges and vertices) from that mesh
  Voronoi* voronoi;

  //The vertex pairs that must be connected.
  std::vector<Subseg*> subsegments;

  //Spatial indexing tree containing the subsegments' bounding boxes
  KDDTree<Subseg*>* subseg_tree;

  //Spatial indexing tree
  KDDTree<Vert*>* vert_tree;

  //The surface sampling queue
  SamplingQueue* sampling_queue;

public:

  //This constructor is used to initialize the RefEdges based on
  //total variation of the tangent. Theoretical guarantees should not
  //hold if this value is set higher than pi / 4.
  InitVolMesh(const double allowed_TVT = 0.245 * M_PI);

  InitVolMesh(RefFace* const ref_face, const double allowed_TVT = 0.245 * M_PI);

  ~InitVolMesh();

  //Populates the vertices vector with the UNDELETED mesh verts.
  void
  get_vertices(std::vector<Vert*>& vertices) const;

  //Populates the subsegs vector with all the contents of the subsegments set.
  //The subsegs marked as deleted will NOT be added.
  void
  get_subsegments(std::vector<Subseg*>& subsegs) const;

  //A series of static functions to know if a point lies inside the
  //diamteral ball of a subsegment.
  static bool
  point_inside_subseg_ball(const CubitVector& point,
			   const Subseg* const subseg);
  static bool
  point_inside_or_on_subseg_ball(const CubitVector& point,
				 const Subseg* const subseg);
  static bool
  vert_inside_subseg_ball(const Vert* const vertex, const Subseg* const subseg);
  static bool
  vert_inside_or_on_subseg_ball(const Vert* const vertex,
				const Subseg* const subseg);

  //Public interface to split a subsegment. Should this be public or private?
  //Will split subseg_to_split:
  //Returns true if subseg_to_split was part of a small angle complex.
  //Can optionally return the newly created vertices and subsegs
  void
  split_subseg(Subseg* const subseg_to_split,
	       std::vector<Vert*>* const new_verts = NULL,
	       std::vector<Subseg*>* const new_subsegs = NULL);

  //The following are functions to output relevant initialization
  //data for visualization. The output is in medit format.

  //Outputs all vertices in mesh (including bounding box) to file.
  void
  output_verts_to_file(const char* const filename) const;

  //Outputs the restricted Delaunay triangulation (faces of the 
  //Delaunay tetrahedralization with their Voronoi edge intersecting
  //the geometric surface) to file.
  void
  output_restricted_delaunay(const char* const filename) const;

  //Outputs subsegments to file
  void
  output_subsegments(const char* const filename) const;

  //Outputs tetrahedral mesh to file
  void
  output_mesh_to_file(const char* const filename) const;

  //Outputs Voronoi edges to file
  void
  output_voronoi_to_file(const char* const filename) const;

private:

  //Should not be used.
  InitVolMesh(const InitVolMesh&) :
      mesh(NULL), voronoi(NULL), subsegments(), subseg_tree(NULL), vert_tree(
      NULL), sampling_queue(NULL)
  {
    assert(0);
  }

  //Should not be used.
  InitVolMesh&
  operator=(const InitVolMesh&)
  {
    assert(0);
    return *this;
  }

  //Creates the vertices corresponding to topology vertices.
  void
  init_vert_verts(std::map<RefVertex*, Vert*>& refvert_to_vert);

  //Creates the vertices lying on topological edges.
  void
  init_edge_verts(const double allowed_TVT,
		  const std::map<RefVertex*, Vert*>& refvert_to_vert);

  //Grooms the neighborhood of small angles and initializes SmallAngleTool::instance.
  void
  init_small_angle_tool(const double minimum_angle = 90.);

  //Will split all the subsgments having, in their diametral sphere, 
  //one (or many) vertex contained in new_vertices. After executing 
  //this function, all subsegments will be strongly Delaunay with respect 
  //to the vertices of new_vertices. To make ALL subsegs strongly Delaunay, 
  //just put all known vertices in the vector.
  void
  strongly_delaunayize(std::vector<Vert*> new_vertices);

  //Splits all the small angle subsegs attached to vertex (vertex must
  //be a key in SmallAngleTool::instance otherwise nothing will happen).
  void
  split_small_angle_subsegs(Vert* const vertex,
			    std::vector<Vert*>* const new_verts = NULL,
			    std::vector<Subseg*>* const new_subsegs = NULL);

  //Carries the operations required for a subseg split.
  //Returns the newly created subsegments in a pair.
  void
  split_subseg_at_param(Subseg* const subseg_to_split, const double split_param,
			const bool small_angle_split,
			std::vector<Vert*>* const new_verts = NULL,
			std::vector<Subseg*>* const new_subsegs = NULL);

  //Same as above but splits automatically at mid param point.
  void
  split_subseg_at_mid_param(Subseg* const subseg_to_split,
			    const bool small_angle_split,
			    std::vector<Vert*>* const new_verts =
			    NULL,
			    std::vector<Subseg*>* const new_subsegs = NULL);

  //Inserts insert_vertex in mesh if insert_vertex splits a subseg.
  void
  insert_on_subseg(Vert* const insert_vertex, Subseg* const insert_subseg);

  //Inserts insert_vertex in mesh if insert_vertex is located on a RefFace
  void
  insert_on_surface(const CubitVector& insert_coord,
		    RefFace* const insert_ref_face,
		    VoronoiEdge* const voronoi_edge);

  //Inserts vertex in mesh. This should only be called if a seed cannot
  //be easily found. Finding the seed will take longer.
  void
  insert_in_tet_mesh(Vert* const vertex);
  //Inserts vertex in mesh using Watson insertion (interior insertion).
  //The Voronoi verts and edges as well as sampling_queue get updated. 
  void
  insert_in_tet_mesh(Vert* const vertex, Cell* const seed_cell);

  //Updates data contained in voronoi after an insertion.
  //Populates new_voronoi_edges with the newly created edges.
  void
  update_voronoi(const std::vector<Cell*>& new_cells,
		 std::vector<VoronoiEdge*>& new_voronoi_edges);

  //Find a seed cell (tet cell with vertex in its circumcenter) for insertion.
  Cell*
  find_seed(Vert* const vertex) const;

  //Deallocates subsegs marked deleted. Alors removes them from 
  //the subsegments vector (member of this class).
  void
  purge_subsegs();

  //Builds a Delaunay tet mesh (without any boundaries) from the vertices
  //contained in all_vertices.
  void
  build_mesh();

  //Curved surface sampling.
  void
  sample_surface();

  //Building stuff... obviously.
  void
  build_voronoi();
  void
  build_sampling_queue();

  //Prints info about the current mesh to screen.
  void
  print_mesh() const;

  //This is a predicate determining if vertex is inside cell's
  //cirumsphere (cell must be a TetCell).
  class InCircumsphere {
    Vert* vertex;
  public:
    InCircumsphere(Vert* _vertex_) :
	vertex(_vertex_)
    {
    }
    bool
    operator()(Cell* const cell) const;
  };

  class VoronoiEdgeBuilder {
    Voronoi* m_voronoi;
    std::vector<VoronoiEdge*>* new_edges;
  public:
    VoronoiEdgeBuilder(Voronoi* const _voronoi_,
		       std::vector<VoronoiEdge*>* new_voronoi_edges) :
	m_voronoi(_voronoi_), new_edges(new_voronoi_edges)
    {
    }
    void
    operator()(Face* const face);
  };

};

#endif
