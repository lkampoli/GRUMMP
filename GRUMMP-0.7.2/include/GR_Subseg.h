#ifndef GR_SUBSEG_H
#define GR_SUBSEG_H

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include "GR_Mesh.h"
#include <vector>

class CubitBox;
class CubitVector;
class RefEdge;
class RefFace;
class RefVertex;
class Subseg;
class Vert;

class SubsegBridge {

  RefFace* m_ref_face; //The surface to which we are bridging with this entity
  Subseg* m_parent;    //The parent subseg on a neighboring curve.
  Vert* m_verts[2]; //The two SURFACE vertices (not the same as the Subseg verts).

  SubsegBridge(const SubsegBridge&);
  SubsegBridge&
  operator=(const SubsegBridge&);

public:

  SubsegBridge() :
      m_ref_face(NULL), m_parent(NULL)
  {
    m_verts[0] = NULL;
    m_verts[1] = NULL;
  }

  ~SubsegBridge();

  void
  set_ref_face(RefFace* const ref_face)
  {
    m_ref_face = ref_face;
  }

  RefFace*
  get_ref_face() const
  {
    return m_ref_face;
  }

  void
  set_parent(Subseg* const parent)
  {
    m_parent = parent;
  }

  Subseg*
  get_parent() const
  {
    return m_parent;
  }

  void
  set_vert(int i, Vert* const vertex)
  {
    assert(i == 0 || i == 1);
    m_verts[i] = vertex;
  }

  Vert*
  get_vert(int i) const
  {
    assert(i == 0 || i == 1);
    return m_verts[i];
  }

  bool
  has_vert(Vert* const vertex) const
  {
    return (m_verts[0] == vertex || m_verts[1] == vertex);
  }

  void
  replace_vert(Vert* const old_vert, Vert* const replace_by)
  {
    assert(has_vert(old_vert));
    if (m_verts[0] == old_vert)
      m_verts[0] = replace_by;
    else
      m_verts[1] = replace_by;
  }

  bool
  parent_deleted() const;

};

class Subseg {

public:

  //We are storing the vertex and its parameter on the RefEdge
  struct VertData {

    Vert* vertex;
    double param;

    VertData() :
	vertex(NULL), param(-LARGE_DBL)
    {
    }

    ~VertData()
    {
    }

    VertData(Vert* const _vertex_, double _param_) :
	vertex(_vertex_), param(_param_)
    {
    }

    VertData(const VertData& VD) :
	vertex(VD.vertex), param(VD.param)
    {
    }

    VertData&
    operator=(const VertData& VD)
    {
      if (&VD != this) {
	vertex = VD.vertex;
	param = VD.param;
      }
      return *this;
    }

  };

private:

  VertData *m_vd1, *m_vd2;

  RefEdge* m_parent_edge;

  // potential for creating new subsegs when splitting
  // geometry, this is used to indicate the subseg lies on this point.
  RefFace* m_ref_face;

  std::vector<SubsegBridge*> m_children;

  bool m_deleted, m_locked;

public:

  Subseg(const VertData& vd1, const VertData& vd2, RefEdge* const edge) :
      m_vd1(new VertData(vd1.vertex, vd1.param)), m_vd2(
	  new VertData(vd2.vertex, vd2.param)), m_parent_edge(edge), m_ref_face(
      NULL), m_children(), m_deleted(false), m_locked(false)
  {
  }

  Subseg(Vert* const v1, Vert* const v2, double param1, double param2,
	 RefEdge* const edge) :
      m_vd1(new VertData(v1, param1)), m_vd2(new VertData(v2, param2)), m_parent_edge(
	  edge), m_ref_face(NULL), m_children(), m_deleted(false), m_locked(
	  false)
  {
  }

  Subseg(Vert* const v1, Vert* const v2) :
      m_vd1(new VertData(v1, 0.)), m_vd2(new VertData(v2, 1.)), m_parent_edge(
      NULL), m_ref_face(NULL), m_children(), m_deleted(false), m_locked(false)
  {
  }

  Subseg(Vert* const v1, Vert* const v2, RefFace * ref_face) :
      m_vd1(new VertData(v1, 0.)), m_vd2(new VertData(v2, 1.)), m_parent_edge(
      NULL), m_ref_face(ref_face), m_children(), m_deleted(false), m_locked(
	  false)
  {
  }

  Subseg() :
      m_vd1(new VertData()), m_vd2(new VertData()), m_parent_edge(NULL), m_ref_face(
      NULL), m_children(), m_deleted(false), m_locked(false)
  {
  }

  Subseg(const Subseg& SS) :
      m_vd1(new VertData(SS.m_vd1->vertex, SS.m_vd1->param)), m_vd2(
	  new VertData(SS.m_vd2->vertex, SS.m_vd2->param)), m_parent_edge(
	  SS.m_parent_edge), m_ref_face(SS.m_ref_face), m_children(
	  SS.m_children), m_deleted(SS.m_deleted), m_locked(SS.m_locked)
  {
  }

  Subseg&
  operator=(const Subseg& SS)
  {
    if (&SS != this) {
      m_vd1->vertex = SS.m_vd1->vertex;
      m_vd2->vertex = SS.m_vd2->vertex;

      m_vd1->param = SS.m_vd1->param;
      m_vd2->param = SS.m_vd2->param;

      m_parent_edge = SS.m_parent_edge;
      m_ref_face = SS.m_ref_face;
      m_children = SS.m_children;

      m_deleted = SS.m_deleted;
      m_locked = SS.m_locked;
    }
    return *this;

  }

  ~Subseg();

  //Marks this subseg as deleted
  void
  mark_deleted()
  {
    m_deleted = true;
  }

  //Checks if the subseg is marked as deleted (or not).
  bool
  is_deleted() const
  {
    return m_deleted;
  }
  bool
  is_not_deleted() const
  {
    return !m_deleted;
  }

  //
  void
  lock()
  {
    m_locked = true;
  }
  void
  unlock()
  {
    m_locked = false;
  }
  bool
  locked()
  {
    return m_locked;
  }

  //Returns the vertices defining this subseg.
  Vert*
  get_beg_vert() const
  {
    return m_vd1->vertex;
  }
  Vert*
  get_end_vert() const
  {
    return m_vd2->vertex;
  }

  //Returns the vertices' parametric location
  double
  get_beg_param() const
  {
    return m_vd1->param;
  }
  double
  get_end_param() const
  {
    return m_vd2->param;
  }

  double
  get_vert_param(const Vert* const vert0) const
  {
    assert(has_vert(vert0));
    return vert0 == m_vd1->vertex ? m_vd1->param : m_vd2->param;
  }

  void
  get_other_vert(const Vert* const vert0, Vert*& vert1) const
  {
    assert(has_vert(vert0));
    vert1 = vert0 == m_vd1->vertex ? m_vd2->vertex : m_vd1->vertex;
  }

  void
  get_other_vert_data(const Vert* const vert0, Vert*& vert1,
		      double& param1) const
  {
    assert(has_vert(vert0));
    if (vert0 == m_vd1->vertex) {
      vert1 = m_vd2->vertex;
      param1 = m_vd2->param;
    }
    else {
      vert1 = m_vd1->vertex;
      param1 = m_vd1->param;
    }
  }

  //Returns the attached RefEdge.
  RefEdge*
  get_ref_edge() const
  {
    return m_parent_edge;
  }
  //Returns the surface this is on
  RefFace*
  get_ref_face() const
  {
    return m_ref_face;
  }
  //Adds a child to this subseg.
  void
  add_bridge(SubsegBridge* const child)
  {
    if (child->get_parent() == NULL)
      child->set_parent(this);
    assert(child->get_parent() == this);
    m_children.push_back(child);
  }

  //Populates the children vector with the children of this subseg.
  void
  get_children(std::vector<SubsegBridge*>& children) const
  {
    children = m_children;
  }

  //Computes the length of this subsegment
  double
  get_length() const;

  //Checks if Subseg has Vert vert;
  bool
  has_vert(const Vert* const vertex) const
  {
    return (m_vd1->vertex == vertex || m_vd2->vertex == vertex);
  }

  //If m_vd1 or m_vd2 have have a RefVertex as parent
  //topological entity, returns the pointer to this RefVertex.
  //Otherwise, a NULL pointer is returned.
  RefVertex*
  get_beg_ref_vert() const;
  RefVertex*
  get_end_ref_vert() const;

  //
  void
  replace_vert(Vert* const old_vert, Vert* const new_vert)
  {
    assert(has_vert(old_vert));
    if (m_vd1->vertex == old_vert)
      m_vd1->vertex = new_vert;
    else
      m_vd2->vertex = new_vert;
  }

  double
  param_at_distance_from_vert(const double distance,
			      const Vert* const vert) const;

  //Returns the parameter ON m_edge at a given distance for an end point.
  double
  param_at_distance_from_beg(const double distance) const
  {
    return param_at_distance_from_vert(distance, get_beg_vert());
  }
  double
  param_at_distance_from_end(const double distance) const
  {
    return param_at_distance_from_vert(distance, get_end_vert());
  }

  //This call is used when performing refinement. It decides between the various
  //types of splits which one is appropriate for this subsegment, computes the
  //insertion location and then performs the split:
  // - REGULAR CASE: Splits at mid param.
  // - NO PARENT CURVE DEFINED: Splits at mid point.
  // - SMALL ANGLE SPLIT: Computes the concentric shell point and splits at this point.

  //This computes the split data only, does not do or assign anything.
  void
  compute_refinement_split_data(double split_coords[3], double& split_param,
				bool qForceCoord = false);
  //Performs the split if split_coords and split_param are known (for function above).
  void
  refinement_split(const double split_coords[3], const double split_param,
		   Vert* const new_vert, Subseg* const new_subseg1,
		   Subseg* const new_subseg2);
  //Performs the split if split_coords and split_param are not yet known.
  void
  refinement_split(Vert* const new_vert, Subseg* const new_subseg1,
		   Subseg* const new_subseg2);

  //Splits the subsegment at param.
  //new_subseg1, new_subseg2 and new_vert must be declared
  //and instantiated elsewhere before calling the function
  void
  split_at_param(const double split_param, Subseg* const new_subseg1,
		 Subseg* const new_subseg2, Vert* const new_vert);

  //Same as above expect that split is done at middle parameter
  //between m_vd1->param and m_vd2->param.
  void
  split_mid_param(Subseg* const new_subseg1, Subseg* const new_subseg2,
		  Vert* const new_vert)
  {
    split_at_param(0.5 * (m_vd1->param + m_vd2->param), new_subseg1,
		   new_subseg2, new_vert);
  }

  //Assuming that we know where to split the subsegment,
  //then go ahead and split it.
  //NOTE: this call does not set the parameter at the subseg's
  //vertices. It needs to be done manually somewhere else is this road is taken.
  void
  split_at_coords(const double split_coords[3], Vert* const new_vert,
		  Subseg* const new_subseg1, Subseg* const new_subseg2)
  {
    split(split_coords, new_vert, new_subseg1, new_subseg2);
  }

  void
  split_at_coords(const CubitVector& split_coords, Vert* const new_vert,
		  Subseg* const new_subseg1, Subseg* const new_subseg2);

  //Returns the subseg midpoint
  CubitVector
  mid_point() const;

  //Returns the ball radius = 0.5 * subseg length
  double
  ball_radius() const;

  //Returns a coordinate-aligned box containing a ball of radius
  //0.5 * subsegment's length and centered at the subsegment's mid point 
  CubitBox
  bounding_box() const;

  //Computes closest point from pt on the subseg 
  //(projects to the underlying curve if it is curved)
  void
  closest_on_subseg(const double pt[3], double close[3]) const;

private:

  void
  split(const double split_coords[3], Vert* const new_vert,
	Subseg* const new_subseg1, Subseg* const new_subseg2);

};

inline bool
SubsegBridge::parent_deleted() const
{
  return m_parent->is_deleted();
}

//This class is basically a wrapper around a TR1 unordered_map
//used to store vertex to subsegment connectivity. This is needed by curved
//boundary as a surface face can actually be connected to 3 vertices
//on the same curve. In such a situation, it is impossible to
//determine which of the face's three edges are actually subsegments.

class SubsegMap {

  typedef std::set<Subseg*> SubsegSet;
  typedef GR_map<Vert*, SubsegSet> SubsegMapType;

  SubsegMapType m_map;

  std::vector<Subseg*> m_to_release;

  //Forbid copy construction and operator= since not needed right now.
  SubsegMap(const SubsegMap&);
  SubsegMap& operator=(const SubsegMap&);

public:

  //Creates and empty map;
  SubsegMap() :
  m_map(), m_to_release() {
  }

  ~SubsegMap();

  // these functions are either internal, or never currently called
private:
  // Release the deleted subsegs in m_to_release.
  // used in destructor
  void release_deleted_subsegs();

  // Clears the map.
  // never called right now
  void clear() {
    m_map.clear();
  }

  //Checks if vertex is connected to subseg.
  bool vert_has_subseg(Vert* const vertex, Subseg* const subseg) const;

public:
  //Checks if the map is empty.
  bool empty() const {
    return m_map.empty();
  }

  //How many entries in the map?
  int size() const {
    return static_cast<int>(m_map.size());
  }

  //Returns the set of subsegs attached to a given vert. Moved from private to public for use in Length3D.
  const std::set<Subseg*>& get_subsegs(Vert* const vert0) const;

  //Add connectivity info for subseg vert0-vert1 to the map.
  void add_subseg(Subseg* const subseg);

  //Removes subseg vert0-vert1 from the map.
  void remove_subseg(Subseg* const subseg);

  //Returns the subseg connecting the two verts (if there is one).
  //Otherwise, returns NULL.
  Subseg* verts_are_connected(Vert* const vert0, Vert* const vert1) const;

  bool isSubsegVert(Vert * vert) const;
  //Returns the set iterator range containing the subsegs of vert0.
  std::pair<SubsegSet::const_iterator, SubsegSet::const_iterator> get_subseg_range(
      Vert* const vert0) const;

  //Populates all_subsegs with all subsegs in the map.
  void get_all_subsegs(std::set<Subseg*>& all_subsegs) const;

  // never called right now
  GR_index_t num_subsegs();

  //Gets a vector of vertices connected to vert0.
  void get_connected_verts(Vert* const vert0,
      std::vector<Vert*>& verts) const;

  // used during purging, replaces a vertex with a new one
  void replaceVert(Vert * oldVert, Vert * newVert);

//  void update_vert_data_after_purge(const std::map<Vert*, Vert*>& vert_map);

  //Checks the validity of the map entries.
  bool valid() const;

  //Prints the content of the map;
  void print() const;

private:

  struct SubsegHasVert {
    const Vert* m_vert;
    SubsegHasVert(const Vert* const vert) :
    m_vert(vert) {
    }
    bool operator()(const Subseg* const subseg) const {
      return subseg->has_vert(m_vert);
    }
  };

};

//This is used to identify edges that already have been checked
//for encroachment while computing the Watson cavity around
//a new vertex.

class TetRefinerEdge {

protected:

  Vert *m_v0, *m_v1;
  BFace* m_bface;

public:

  TetRefinerEdge() :
      m_v0(NULL), m_v1(NULL), m_bface(NULL)
  {
  }
  TetRefinerEdge(Vert* const v0, Vert* const v1) :
      m_v0(v0 < v1 ? v0 : v1), m_v1(v0 < v1 ? v1 : v0), m_bface(NULL)
  {
    assert(valid());
  }
  TetRefinerEdge(Vert* const v0, Vert* const v1, BFace* const bface) :
      m_v0(v0 < v1 ? v0 : v1), m_v1(v0 < v1 ? v1 : v0), m_bface(bface)
  {
    assert(valid());
  }
  TetRefinerEdge(const TetRefinerEdge& edge) :
      m_v0(edge.m_v0), m_v1(edge.m_v1), m_bface(edge.m_bface)
  {
    assert(valid());
  }
  TetRefinerEdge&
  operator=(const TetRefinerEdge& edge)
  {
    if (&edge != this) {
      m_v0 = edge.m_v0;
      m_v1 = edge.m_v1;
      m_bface = edge.m_bface;
      assert(valid());
    }
    return *this;
  }
  virtual
  ~TetRefinerEdge()
  {
  }

  Vert*
  vert(int i) const
  {
    assert(i == 0 || i == 1);
    return i == 0 ? m_v0 : m_v1;
  }
  Vert*
  beg_vert() const
  {
    return m_v0;
  }
  Vert*
  end_vert() const
  {
    return m_v1;
  }

  BFace*
  get_bface() const
  {
    return m_bface;
  }
  int
  is_reflex() const;

  bool
  operator<(const TetRefinerEdge& edge) const
  {
    assert(valid() && edge.valid());
    return ((m_v0 < edge.m_v0) || ((m_v0 == edge.m_v0) && (m_v1 < edge.m_v1)));
  }
  bool
  operator==(const TetRefinerEdge& edge) const
  {
    assert(valid() && edge.valid());
    return (m_v0 == edge.m_v0 && m_v1 == edge.m_v1);
  }
  bool
  operator!=(const TetRefinerEdge& edge) const
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

//This derived class is used to infer subsegments given
//a tetrahedral mesh using GRUMMP's native geometry. Note
//that this native geometry does not allow orphaned edges
//or surfaces not bounding at least on meshed region.

class TetRefinerEdgeBFace : public TetRefinerEdge {

  std::set<BFace*> m_bfaces;

public:

  TetRefinerEdgeBFace() :
      TetRefinerEdge(), m_bfaces()
  {
  }

  TetRefinerEdgeBFace(Vert* const v0, Vert* const v1) :
      TetRefinerEdge(v0, v1), m_bfaces()
  {
    assert(valid());
  }
  TetRefinerEdgeBFace(const TetRefinerEdgeBFace& edge) :
      TetRefinerEdge(edge), m_bfaces(edge.m_bfaces)
  {
    assert(valid());
  }
  TetRefinerEdgeBFace&
  operator=(const TetRefinerEdgeBFace& edge)
  {
    m_v0 = edge.m_v0;
    m_v1 = edge.m_v1;
    m_bfaces = edge.m_bfaces;
    assert(valid());
    return *this;
  }
  virtual
  ~TetRefinerEdgeBFace()
  {
  }

  const std::set<BFace*>&
  get_bfaces() const
  {
    return m_bfaces;
  }

  void
  add_bface(BFace* const bface)
  {
    m_bfaces.insert(bface);
  }

  void
  remove_bface(BFace* const bface)
  {
    m_bfaces.erase(bface);
  }

  void
  clear_bfaces()
  {
    m_bfaces.clear();
  }

};

#endif
