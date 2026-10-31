#include "GR_VertexPair.h"
#include "GR_misc.h"
#include "GR_GRCurve.h"
#include "GR_GRPoint.h"
#include "GR_Vertex.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"

#include <utility>
#include <algorithm>

using std::pair;
using std::make_pair;

VertexPair::VertexPair() :
    m_deleted(false), m_curve(NULL), m_vert_pair()
{

  VertData vd1, vd2;
  vd1.vertex = vd2.vertex = NULL;
  vd1.vert_param = vd2.vert_param = LARGE_DBL;

  m_vert_pair = make_pair(vd1, vd2);

}

VertexPair::VertexPair(GRCurve* const curve,
		       std::pair<VertData, VertData>& vert_pair) :
    m_deleted(false), m_curve(curve), m_vert_pair(vert_pair)
{
}

VertexPair::VertexPair(GRCurve* const curve, const VertData& vert_data_1,
		       const VertData& vert_data_2) :
    m_deleted(false), m_curve(curve), m_vert_pair(
	make_pair(vert_data_1, vert_data_2))
{
  assert(iFuzzyComp(vert_data_1.vert_param, vert_data_2.vert_param) == -1);
}

VertexPair::VertexPair(GRCurve* const curve, Vert* const vert1,
		       Vert* const vert2, double vert1_param,
		       double vert2_param) :
    m_deleted(false), m_curve(curve), m_vert_pair()
{
  assert(iFuzzyComp(vert1_param, vert2_param) == -1);

  VertData vd1, vd2;

  vd1.vertex = vert1;
  vd2.vertex = vert2;

  vd1.vert_param = vert1_param;
  vd2.vert_param = vert2_param;

  m_vert_pair = make_pair(vd1, vd2);
}

VertexPair::VertexPair(const VertexPair& vertex_pair) :
    m_deleted(vertex_pair.m_deleted), m_curve(vertex_pair.m_curve), m_vert_pair(
	vertex_pair.m_vert_pair)
{
}

VertexPair&
VertexPair::operator=(const VertexPair& vertex_pair)
{
  if (&vertex_pair != this) {
    m_deleted = vertex_pair.m_deleted;
    m_curve = vertex_pair.m_curve;
    m_vert_pair = vertex_pair.m_vert_pair;
  }
  return *this;

}

VertexPair::~VertexPair()
{
}

GRPoint*
VertexPair::beg_point() const
{

  return dynamic_cast<GRPoint*>(beg_vert()->getParentEntity());

}

GRPoint*
VertexPair::end_point() const
{

  return dynamic_cast<GRPoint*>(end_vert()->getParentEntity());

}

double
VertexPair::length() const
{

  Vert* v1 = m_vert_pair.first.vertex;
  Vert* v2 = m_vert_pair.second.vertex;

  double vector[] =
    { v2->x() - v1->x(), v2->y() - v1->y() };

  return dMAG2D(vector);

}

double
VertexPair::length_squared() const
{

  Vert* v1 = m_vert_pair.first.vertex;
  Vert* v2 = m_vert_pair.second.vertex;

  double diffx = v2->x() - v1->x();
  double diffy = v2->y() - v1->y();

  return (diffx * diffx + diffy * diffy);

}

double
VertexPair::tangent_variation() const
{

  double param1 = m_vert_pair.first.vert_param;
  double param2 = m_vert_pair.second.vert_param;

  return m_curve->get_curve_geom()->TVT(param1, param2);

}

double
VertexPair::mid_tvt_param() const
{

  return m_curve->get_curve_geom()->mid_TVT(beg_param(), end_param());

}

// CubitBox VertexPair::
// bounding_box() const {

//   Vert* v1 = m_vert_pair.first.vertex;
//   Vert* v2 = m_vert_pair.second.vertex;

//   CubitVector mini( std::min(v1->dX(), v2->dX()),
// 		    std::min(v1->dY(), v2->dY()),
// 		    0. );

//   CubitVector maxi( std::max(v1->dX(), v2->dX()),
// 		    std::max(v1->dY(), v2->dY()),
// 		    0. );

//   return CubitBox(mini, maxi);

// }

CubitBox
VertexPair::bounding_box() const
{

  CubitVector midpoint = mid_point();
  double rad = radius();

  double x1 = midpoint.x() - rad, x2 = midpoint.x() + rad;
  double y1 = midpoint.y() - rad, y2 = midpoint.y() + rad;

  CubitVector max_coord(std::max(x1, x2), std::max(y1, y2), 0.);
  CubitVector min_coord(std::min(x1, x2), std::min(y1, y2), 0.);

  return CubitBox(min_coord, max_coord);

}

void
VertexPair::split(VertexPair* const new_pair1, VertexPair* const new_pair2,
		  Vert* const new_vert, const SplitType split_type,
		  const bool allow_concentric_shell_split)
{

  assert(new_pair1);
  assert(new_pair2);
  assert(new_vert);

  CubitVector split_coord;
  double split_param;
  bool concentric_shell_split;

  compute_split_data(split_type, allow_concentric_shell_split, split_coord,
		     split_param, concentric_shell_split);

  split(split_coord, split_param, concentric_shell_split, new_pair1, new_pair2,
	new_vert);

}

void
VertexPair::initial_shell_split(double distance, VertexPair* const new_pair1,
				VertexPair* const new_pair2,
				Vert* const new_vert)
{

  assert(distance < length());
  assert(new_pair1);
  assert(new_pair2);
  assert(new_vert);

  assert(beg_is_small() || end_is_small());
  assert(!(beg_is_small() && end_is_small()));
  assert(!beg_is_shell());
  assert(!end_is_shell());

  CubitVector split_coord;
  double split_param;

  if (beg_is_small()) {
    assert(
	iFuzzyComp(beg_param(), m_curve->get_curve_geom()->min_param()) == 0);
    split_param = m_curve->get_curve_geom()->coord_at_dist(distance, true,
							   split_coord);
  }
  else {
    assert(
	iFuzzyComp(end_param(), m_curve->get_curve_geom()->max_param()) == 0);
    split_param = m_curve->get_curve_geom()->coord_at_dist(distance, false,
							   split_coord);
  }

  split(split_coord, split_param, true, new_pair1, new_pair2, new_vert);

}

void
VertexPair::split(const CubitVector& split_coord, const double split_param,
		  const bool concentric_shell_split,
		  VertexPair* const new_pair1, VertexPair* const new_pair2,
		  Vert* const new_vert)
{

  new_vert->setCoords(2, split_coord);
  new_vert->setParentEntity(m_curve);
  new_vert->setType(Vert::eBdry);

  if (concentric_shell_split)
    new_vert->markAsShellVert(true);

  new_pair1->set_curve(m_curve);
  new_pair2->set_curve(m_curve);

  new_pair1->set_beg_vert(beg_vert());
  new_pair1->set_end_vert(new_vert);
  new_pair2->set_beg_vert(new_vert);
  new_pair2->set_end_vert(end_vert());

  assert(iFuzzyComp(beg_param(), split_param) == -1);
  assert(iFuzzyComp(split_param, end_param()) == -1);

  new_pair1->set_beg_param(beg_param());
  new_pair1->set_end_param(split_param);
  new_pair2->set_beg_param(split_param);
  new_pair2->set_end_param(end_param());

}

void
VertexPair::compute_split_data(const SplitType split_type,
			       const bool allow_concentric_shell_split,
			       CubitVector& split_coord, double& split_param,
			       bool& concentric_shell_split) const
{

  split_coord.set(-LARGE_DBL, -LARGE_DBL, 0.);
  split_param = LARGE_DBL;
  concentric_shell_split = false;

  if (allow_concentric_shell_split) {

    if (beg_is_small() && end_is_shell()) {
      split_param = m_curve->get_curve_geom()->coord_at_dist(0.5 * length(),
							     true, split_coord);
      concentric_shell_split = true;
      return;
    }
    if (end_is_small() && beg_is_shell()) {
      split_param = m_curve->get_curve_geom()->coord_at_dist(0.5 * length(),
							     false,
							     split_coord);
      concentric_shell_split = true;
      return;
    }
    if (beg_is_shell() && end_is_shell()) {
      split_param = m_curve->get_curve_geom()->coord_at_mid_dist(beg_param(),
								 end_param(),
								 split_coord);
      concentric_shell_split = true;
      return;
    }

  }

  switch (split_type)
    {

    case (MID_TVT):
      split_param = mid_tvt_param();
      break;

    case (MID_PARAM):
      split_param = 0.5 * (beg_param() + end_param());
      break;

    case (EQUIDISTANT):
      vFatalError(
	  "Equidistant splitting has not been implemented yet... sorry!",
	  "VertexPair::compute_split_data()");
      break;

    case (INIT_SHELL):
      vFatalError(
	  "To initialize shell splitting, use VertexPair::initial_shell_split()",
	  "VertexPair::compute_split_data()");
      break;
    default:
      assert(0);
      break;

    }

  m_curve->get_curve_geom()->coord_at_param(split_param, split_coord);

}

