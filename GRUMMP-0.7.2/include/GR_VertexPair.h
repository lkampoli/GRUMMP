#ifndef GR_VERTEXPAIR_H
#define GR_VERTEXPAIR_H

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Vertex.h"
#include <utility>

#include "CubitVector.hpp"

class GRCurve;
class GRPoint;
class CubitBox;

//This class is used to build the 2D mesh.
//WILL ONLY WORK IN 2D because GRCurves are only implemented in 2D...

class VertexPair {

public:

  struct VertData {

    Vert* vertex;
    double vert_param;

  };

private:

  bool m_deleted;

  GRCurve* m_curve;

  std::pair<VertData, VertData> m_vert_pair;

public:

  enum SplitType {
    MID_TVT, MID_PARAM, EQUIDISTANT, INIT_SHELL, UNKNOWN
  };

  //Constructors and destructor

  VertexPair();

  VertexPair(GRCurve* const curve, std::pair<VertData, VertData>& vert_pair);

  VertexPair(GRCurve* const curve, const VertData& vert_data_1,
	     const VertData& vert_data_2);

  VertexPair(GRCurve* const curve, Vert* const vertex_1, Vert* const vertex_2,
	     double v1_param = LARGE_DBL, double v2_param = LARGE_DBL);

  VertexPair(const VertexPair& vertex_pair);

  VertexPair&
  operator=(const VertexPair& vertex_pair);

  ~VertexPair();

  //Marks this vertex pair as deleted.
  void
  mark_deleted()
  {
    m_deleted = true;
  }

  //Checks if this vertex pair is deleted or not.
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

  //Sets the parent curve
  void
  set_curve(GRCurve* const curve)
  {
    m_curve = curve;
  }

  //Returns the parent curve
  GRCurve*
  get_curve() const
  {
    return m_curve;
  }

  //Sets bounding vertices
  void
  set_beg_vert(Vert* const vertex)
  {
    m_vert_pair.first.vertex = vertex;
  }
  void
  set_end_vert(Vert* const vertex)
  {
    m_vert_pair.second.vertex = vertex;
  }

  //Returns bounding vertices
  Vert*
  beg_vert() const
  {
    return m_vert_pair.first.vertex;
  }
  Vert*
  end_vert() const
  {
    return m_vert_pair.second.vertex;
  }

  //Returns the point associated with
  //beg_vert and end_vert if these are apex vertices.
  GRPoint*
  beg_point() const;
  GRPoint*
  end_point() const;

  //Sets bounding vertices parameter
  void
  set_beg_param(double param)
  {
    assert(iFuzzyComp(param, end_param()) == -1);
    m_vert_pair.first.vert_param = param;
  }
  void
  set_end_param(double param)
  {
    assert(iFuzzyComp(param, beg_param()) == 1);
    m_vert_pair.second.vert_param = param;
  }

  //Returns bounding vertices parameter
  double
  beg_param() const
  {
    return m_vert_pair.first.vert_param;
  }
  double
  end_param() const
  {
    return m_vert_pair.second.vert_param;
  }

  //Sets/unsets bounding vertices as small angle vert or shell vert
  void
  set_beg_small(bool is_small)
  {
    m_vert_pair.first.vertex->markAsSmallAngleVert(is_small);
  }
  void
  set_beg_shell(bool is_shell)
  {
    m_vert_pair.first.vertex->markAsShellVert(is_shell);
  }
  void
  set_end_small(bool is_small)
  {
    m_vert_pair.second.vertex->markAsSmallAngleVert(is_small);
  }
  void
  set_end_shell(bool is_shell)
  {
    m_vert_pair.second.vertex->markAsShellVert(is_shell);
  }

  //Returns bounding vertices small / shell status
  bool
  beg_is_small() const
  {
    return m_vert_pair.first.vertex->isSmallAngleVert();
  }
  bool
  beg_is_shell() const
  {
    return m_vert_pair.first.vertex->isShellVert();
  }
  bool
  end_is_small() const
  {
    return m_vert_pair.second.vertex->isSmallAngleVert();
  }
  bool
  end_is_shell() const
  {
    return m_vert_pair.second.vertex->isShellVert();
  }

  //Determines if vertex is one of the verts in the pair.
  bool
  has_vert(const Vert* const vertex) const
  {
    return (vertex == beg_vert() || vertex == end_vert());
  }

  //Determines if beg_vert and end_vert are apex vertices.
  bool
  has_two_apex_verts() const
  {
    return (beg_point() && end_point());
  }

  //Returns half the length of this
  double
  radius() const
  {
    return 0.5 * length();
  }

  //Return the mid point of this
  CubitVector
  mid_point() const
  {
    return CubitVector(0.5 * (beg_vert()->x() + end_vert()->x()),
		       0.5 * (beg_vert()->y() + end_vert()->y()), 0.);
  }

  //Returns the Euclidian distance between bounding vertices.
  double
  length() const;
  double
  length_squared() const;

  //Computes and returns the total tangent variation between vertices.
  double
  tangent_variation() const;

  //Computes and returns the param of the curve such that it splits
  //the vertex pair into two segments spanning equal TVT.
  double
  mid_tvt_param() const;

  //Computes a coordinate-aligned bounding box
  CubitBox
  bounding_box() const;

  //Public interface call to split this VertexPair segment
  void
  split(VertexPair* const new_pair1, VertexPair* const new_pair2,
	Vert* const new_vert, const SplitType split_type = MID_TVT,
	const bool allow_concentric_shell_split = true);

  //Function to call to perform the INITIAL concentric shell split.
  //ONE of the vertex must be a small angle vertex, the other
  //must be a regular vertex.
  //distance is the variable to specify the distance of the new vertex
  //from the small angle vert. It must obviously be smaller than length().
  void
  initial_shell_split(double distance, VertexPair* const new_pair1,
		      VertexPair* const new_pair2, Vert* const new_vert);

private:

  //Functions performing the necessary operations to perform a split
  //of the vertex-bounded segment.

  void
  split(const CubitVector& split_coord, const double split_param,
	const bool concentric_shell_split, VertexPair* const new_pair1,
	VertexPair* const new_pair2, Vert* const new_vert);

  void
  split_at_coord(const CubitVector& coord, VertexPair* const new_pair1,
		 VertexPair* const new_pair2, Vert* const new_vert,
		 const bool allow_concetric_shell_split);

  void
  compute_split_data(const SplitType split_type,
		     const bool allow_concentric_shell_split,
		     CubitVector& split_coord, double& split_param,
		     bool& concentric_shell_split) const;

};

#endif
