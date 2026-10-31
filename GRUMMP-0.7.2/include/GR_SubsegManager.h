#ifndef GR_SUBSEGMANAGER_H
#define GR_SUBSEGMANAGER_H

#include "GR_config.h"
#include "GR_BaryCentricLengthScale.h"
#include "GR_CurveDiscretizer.h"
#include <vector>
#include <map>
#include "KDDTree.hpp"

class BasicTopologyEntity;
class CubitVector;
class RefEdge;
class RefFace;
class RefVertex;
class Subseg;
class Vert;

class SubsegManager {

  //POSSIBLE IMPROVEMENT: Currently storing Vert* and Subseg* in a vector
  //and in the trees. It should be possible to only store the pointers in 
  //the trees. The subsegs are also stored a third time in edge_subseg_map.
  //Triple storage == not good!
  BasicTopologyEntity * m_entity;

  //Vector containing all instantiated vertices.
  std::vector<Vert*> vertices;

  //Vector containing all instantiated subsegments.
  std::vector<Subseg*> subsegments;

  //Contains the bounding boxes of all verts.
  KDDTree<Vert*>* vert_tree;

  //Contains the bounding boxes of all subsegs.
  KDDTree<Subseg*>* subseg_tree;

  //A "dictionary" allowing to retrieve the vertex of a given RefVertex.
  typedef std::map<RefVertex*, Vert*> ApexVertMap;
  ApexVertMap apex_vert_map;

  //A "dictionary" allowing to retrieve the subsegments of a given RefEdge. 
  typedef std::map<RefEdge*, std::vector<Subseg*> > EdgeSubsegMap;
  EdgeSubsegMap edge_subseg_map;

  // allows for sampling over length scale
  const VolMesh * m_backgroundMesh;

  // allows for more subseg flexibility
  CurveDiscretizer::CDType m_samplingType;

public:

  //Constructor, builds a set of subsegments for the model.
  //Curve discretization is based on tangent variation.
  SubsegManager(BasicTopologyEntity * entity,
		const double allowed_TVT = (60. / 180.) * M_PI,
		const VolMesh * backgroundMesh = NULL,
		CurveDiscretizer::CDType samplingType =
		    CurveDiscretizer::TANGENT);

  //Constructor, creates it off of subsegments from a surface mesh
  SubsegManager(SurfMesh * surf, BasicTopologyEntity * entity);

  //Destructor
  ~SubsegManager();

  //Populates the vector with the UNDELETED entries in vertices.
  void
  get_vertices(std::vector<Vert*>& verts) const;

  //Populates the vector with the UNDELETED entries in subsegments.
  void
  get_subsegments(std::vector<Subseg*>& subsegs) const;

  //Returns the vertex associated to ref_vert. If ref_vert is not
  //present in apex_vert_map, returns a NULL pointer.
  Vert*
  get_apex_vertex(RefVertex* const ref_vert) const;

  //Populates the subsegs vector with the subsegments attached to ref_edge.
  //The subsegments marked as deleted ARE NOT added to subsegs. If ref_edge
  //is not contained in the map, adds nothing to the vector.
  void
  get_subsegments_of_edge(RefEdge* const ref_edge,
			  std::vector<Subseg*>& subsegs) const;

  //Returns a subseg having coord in its diametral sphere, if one exists. 
  //Otherwise returns NULL. 
  //Passing a non-NULL ref_face pointer will restrict the search 
  //to the subsegments of ref_face
  Subseg*
  point_encroaching_on_subseg(const CubitVector& point,
			      const RefFace* const ref_face = NULL) const;

  //A series of static functions to know if a point lies inside the
  //diametral ball of a subsegment.
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

  //Public interface to split a subsegment. 
  //Will split subseg_to_split:
  //Can optionally return the newly created vertices and subsegs
  void
  split_subseg(Subseg* const subseg_to_split,
	       std::vector<Vert*>* const new_verts = NULL,
	       std::vector<Subseg*>* const new_subsegs = NULL,
	       std::vector<Subseg*>* const old_subsegs = NULL);

  //Garbage collection.
  void
  purge()
  {
    purge_verts();
    purge_subsegs();
  }

  //The following are functions to output some relevant
  //data for visualization. The output is in medit format.

  //Outputs vertices to file in medit format.
  void
  output_verts_to_file(const char* const filename) const;

  //Outputs subsegments to file in medit format.
  void
  output_subsegs_to_file(const char* const filename) const;

private:

  //Copy constructor and operator= are private. There should
  //only be one instance of this class anyway (make it singleton?)
  SubsegManager(const SubsegManager&) :
      m_entity(NULL), vertices(), subsegments(), vert_tree(NULL), subseg_tree(
      NULL), apex_vert_map(), edge_subseg_map(), m_backgroundMesh(
      NULL), m_samplingType(CurveDiscretizer::TANGENT)
  {
    assert(0);
  }

  SubsegManager&
  operator=(const SubsegManager&)
  {
    assert(0);
    return *this;
  }

  //This function will create apex and edge vertices and then
  //connect them appropriately to form the initial set of subsegments.
  void
  init_subsegs(const double allowed_TVT);

  //Creates the vertices corresponding to topological vertices.
  void
  create_apex_verts();

  //Creates the vertices lying on topological edges.
  void
  create_edge_verts(const double allowed_TVT);
  //Initializes SmallAngleTool::instance.
  void
  init_small_angle_tool(const double minimum_angle = 90.);

  //Will split all the subsgments having, in their diametral sphere, 
  //one (or many) vertex contained in new_vertices. After executing 
  //this function, all subsegments will be strongly Delaunay with respect 
  //to the vertices of new_vertices. To make ALL subsegs strongly Delaunay, 
  //just put all known vertices in the vector.
  void
  strongly_delaunayize(std::vector<Vert*>* verts = NULL);

  //Splits all the small angle subsegs attached to vertex (vertex must
  //be a key in SmallAngleTool::instance otherwise nothing will happen).
  void
  split_small_angle_subsegs(Vert* const vertex,
			    std::vector<Vert*>* const new_verts = NULL,
			    std::vector<Subseg*>* const new_subsegs = NULL,
			    std::vector<Subseg*>* const old_subsegs = NULL);

  //Carries the operations required for a subseg split.
  //Returns the newly created subsegments in a pair.
  void
  split_subseg_at_param(Subseg* const subseg_to_split, const double split_param,
			const bool small_angle_split,
			std::vector<Vert*>* const new_verts = NULL,
			std::vector<Subseg*>* const new_subsegs = NULL,
			std::vector<Subseg*>* const old_subsegs = NULL);

  //Split subsegs at mid parameter ( 0.5 * (param_beg + param_end) ).
  void
  split_subseg_at_mid_param(Subseg* const subseg_to_split,
			    const bool small_angle_split,
			    std::vector<Vert*>* const new_verts =
			    NULL,
			    std::vector<Subseg*>* const new_subsegs = NULL,
			    std::vector<Subseg*>* const old_subsegs = NULL);

  //Deallocates verts marked deleted. Also removes them from 
  //the vertices vector (member of this class).
  void
  purge_verts();

  //Deallocates subsegs marked deleted. Also removes them from 
  //the subsegments vector (member of this class).
  void
  purge_subsegs();

};

#endif
