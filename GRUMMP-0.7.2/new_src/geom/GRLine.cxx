#include "GR_GRCurveGeom.h"
#include <iostream>
#include <algorithm>
#include "CubitBox.hpp"
using namespace std;

GRLine::GRLine() :
    GRCurveGeom()
{

  n_pts = 2;
  n_dim = 0;
  points = new CubitVector[2];
  points[0].x(-LARGE_DBL);
  points[0].y(-LARGE_DBL);
  points[0].z(-LARGE_DBL);
  points[1].x(-LARGE_DBL);
  points[1].y(-LARGE_DBL);
  points[1].z(-LARGE_DBL);

}

GRLine::GRLine(const CubitVector& coord_start, const CubitVector& coord_end) :
    GRCurveGeom()
{

  n_pts = 2;
  points = new CubitVector[n_pts];
  points[0].set(coord_start.x(), coord_start.y(), coord_start.z());
  points[1].set(coord_end.x(), coord_end.y(), coord_end.z());

  if (iFuzzyComp(points[0].z(), 0.) == 0
      && iFuzzyComp(points[1].z(), 0.) == 0) {
    n_dim = 2;
  }
  else {
    n_dim = 3;
  }

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    vFatalError("Both points of the line are at the same coordinate.",
		"GRLine::GRLine(const CubitVector&, const CubitVector&)");
  }

}

GRLine::GRLine(const GRLine& L) :
    GRCurveGeom()
{

  n_pts = L.n_pts;
  n_dim = L.n_dim;
  points = new CubitVector[n_pts];
  points[0].set(L.points[0].x(), L.points[0].y(), L.points[0].z());
  points[1].set(L.points[1].x(), L.points[1].y(), L.points[1].z());

}

GRLine&
GRLine::operator=(const GRLine& L)
{

  n_pts = L.n_pts;
  n_dim = L.n_dim;
  points[0].set(L.points[0].x(), L.points[0].y(), L.points[0].z());
  points[1].set(L.points[1].x(), L.points[1].y(), L.points[1].z());

  return *this;

}

GRLine::~GRLine()
{
}

CubitBox
GRLine::bounding_box() const
{

  assert(isValid());

  double xmin = std::min(points[0].x(), points[1].x());
  double ymin = std::min(points[0].y(), points[1].y());
  double zmin = std::min(points[0].z(), points[1].z());
  double xmax = std::max(points[0].x(), points[1].x());
  double ymax = std::max(points[0].y(), points[1].y());
  double zmax = std::max(points[0].z(), points[1].z());

  CubitVector v_min(xmin, ymin, zmin);
  CubitVector v_max(xmax, ymax, zmax);
  CubitBox bbox(v_min, v_max);

  return bbox;

}

double
GRLine::param_at_coord(const CubitVector& coord) const
{

  assert(isValid());
  double param_x, param_y, param_z;
  bool on_curve = coord_on_curve(coord, &param_x, &param_y, &param_z);
  assert(on_curve);
#ifndef NDEBUG
  switch (n_dim)
    {
    case 2:
      assert(iFuzzyComp(param_x, param_y) == 0);
      break;
    case 3:
      assert(iFuzzyComp(param_x, param_y) == 0);
      assert(iFuzzyComp(param_x, param_z) == 0);
      assert(iFuzzyComp(param_y, param_z) == 0);
      break;
    default:
      assert(0);
      break;
    }
#endif

  return param_x;

}

double
GRLine::param_at_arc_length(const double param_start, const double length) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param_start, min_param()) >= 0
	  && iFuzzyComp(param_start, max_param()) <= 0);
  assert(iFuzzyComp(length, 0.) >= 0);

  if (iFuzzyComp(length, arc_length(param_start, max_param())) > 0) {
    return max_param();
  }

  CubitVector startCoord;
  coord_at_param(param_start, startCoord);
  CubitVector dir = points[1] - points[0];
  dir.normalize();
  dir *= length;

  CubitVector coord = startCoord + dir;
  return param_at_coord(coord);

}

void
GRLine::coord_at_param(const double param, CubitVector& coord) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  coord.x(points[0].x() + param * (points[1].x() - points[0].x()));
  coord.y(points[0].y() + param * (points[1].y() - points[0].y()));
  if (n_dim == 3)
    coord.z(points[0].z() + param * (points[1].z() - points[0].z()));
  else
    coord.z(0.);

}

double
GRLine::coord_at_dist(const double dist, const bool start,
		      CubitVector& coord) const
{

  assert(isValid());
  assert(dist >= 0.);

  CubitVector v(points[1].x() - points[0].x(), points[1].y() - points[0].y(),
		0.);

  double length = v.length();
  assert(iFuzzyComp(length, 0.) == 1);

  if (dist > length) {
    if (start) {
      coord.set(points[1].x(), points[1].y(), 0.);
      return max_param();
    }
    else {
      coord.set(points[0].x(), points[0].y(), 0.);
      return min_param();
    }
  }

  v.length(dist);

  if (start) {
    coord = points[0] + v;
    return (dist / length);
  }
  else {
    coord = points[1] - v;
    return ((length - dist) / length);
  }

}

double
GRLine::coord_at_mid_dist(const double param1, const double param2,
			  CubitVector& middle) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  double mid_param = 0.5 * (param1 + param2);
  coord_at_param(mid_param, middle);

  return mid_param;

}

//double GRLine::
//coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
//		  CubitVector& middle) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  middle.set( 0.5 * (coord1.x() + coord2.x()),
//	      0.5 * (coord1.y() + coord2.y()),
//	      0. );
//
//  return param_at_coord(middle);
//
//}

void
GRLine::closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			       const double param1, const double param2) const
{

  assert(isValid());
  if (n_dim == 2)
    assert(iFuzzyComp(coord.z(), 0.) == 0);

  double p1 = min_param(), p2 = max_param();
  if (iFuzzyComp(param1, min_param()) >= 0
      && iFuzzyComp(param2, max_param()) <= 0) {
    assert(iFuzzyComp(param1, param2) <= 0);
    assert(
	iFuzzyComp(param1, min_param()) >= 0
	    && iFuzzyComp(param1, max_param()) <= 0);
    assert(
	iFuzzyComp(param2, min_param()) >= 0
	    && iFuzzyComp(param2, max_param()) <= 0);
    p1 = param1;
    p2 = param2;
  }

  CubitVector c1, c2;
  coord_at_param(p1, c1);
  coord_at_param(p2, c2);

  if (coord.about_equal(c1)) {
    closest.set(c1);
    return;
  }

  CubitVector V1, V2;
  V1.set(coord.x() - c1.x(), coord.y() - c1.y(), coord.z() - c1.z());
  V2.set(c2.x() - c1.x(), c2.y() - c1.y(), c2.z() - c1.z());

  double mag = V2.length();
  assert(mag > 0);
  V2.normalize();

  double dot = V1 % V2;

  if (dot < 0.) { //Closest point is c1.
    closest.set(c1);
  }
  else if (dot > mag) { //Closest point is c2.
    closest.set(c2);
  }
  else {
    closest.set(c1.x() + dot * V2.x(), c1.y() + dot * V2.y(),
		c1.z() + dot * V2.z());
  }

}

double
GRLine::closest_param_on_curve(const CubitVector& coord, const double param1,
			       const double param2) const
{

  assert(isValid());

  CubitVector closest;
  closest_coord_on_curve(coord, closest, param1, param2);
  return (param_at_coord(closest));

}

bool
GRLine::coord_on_curve(const CubitVector& coord, double* const param_x,
		       double* const param_y, double* const param_z) const
{
	//cout<<"@@@@"<<endl;
  assert(isValid());
  double xparam, yparam, zparam;
  bool xequal = false, yequal = false, zequal = false;

  if (iFuzzyComp(points[0].x(), points[1].x()) != 0)
    xparam = (coord.x() - points[0].x()) / (points[1].x() - points[0].x());
  else {
    if (iFuzzyComp(coord.x(), points[0].x()) == 0)
      xequal = true;
    else
      return false;
  }

  if (iFuzzyComp(points[0].y(), points[1].y()) != 0)
    yparam = (coord.y() - points[0].y()) / (points[1].y() - points[0].y());
  else {
    if (iFuzzyComp(coord.y(), points[0].y()) == 0)
      yequal = true;
    else
      return false;
  }

  if (iFuzzyComp(points[0].z(), points[1].z()) != 0)
    zparam = (coord.z() - points[0].z()) / (points[1].z() - points[0].z());
  else {
    if (iFuzzyComp(coord.z(), points[0].z()) == 0)
      zequal = true;
    else
      return false;
  }

  switch (n_dim)
    {

    case 2:
      {

	if (xequal && yequal) {
	  vFatalError("The two points defining the line are the same",
		      "GRLine::coord_on_curve(...)");
	  return false;
	}

	if (xequal && !yequal)
	  xparam = yparam;
	if (yequal && !xequal)
	  yparam = xparam;

	if (iFuzzyComp(xparam, yparam) == 0 && iFuzzyComp(xparam, 0.) >= 0
	    && iFuzzyComp(xparam, 1.) <= 0) {
	  if (param_x != NULL) {
	    *param_x = xparam;
	  }
	  if (param_y != NULL) {
	    *param_y = yparam;
	  }
	  return true;
	}
	else
	  return false;
	break;
      }

    case 3:
      {
	if (xequal && yequal && zequal) {
	  vFatalError("The two points defining the line are the same",
		      "GRLine::coord_on_curve(...)");
	  return false;
	}

	if (xequal && !yequal)
	  xparam = yparam;
	else if (xequal && !zequal)
	  xparam = zparam;
	if (yequal && !xequal)
	  yparam = xparam;
	else if (yequal && !zequal)
	  yparam = zparam;
	if (zequal && !xequal)
	  zparam = xparam;
	else if (zequal && !yequal)
	  zparam = yparam;

	if (iFuzzyComp(xparam, yparam) == 0 && iFuzzyComp(xparam, zparam) == 0
	    && iFuzzyComp(yparam, zparam) == 0 && iFuzzyComp(xparam, 0.) >= 0
	    && iFuzzyComp(xparam, 1.) <= 0) {
	  if (param_x != NULL) {
	    *param_x = xparam;
	  }
	  if (param_y != NULL) {
	    *param_y = yparam;
	  }
	  if (param_z != NULL) {
	    *param_z = zparam;
	  }
	  return true;
	}
	else
	  return false;
	break;
      }

    }
  assert(0);
  return false;

}

void
GRLine::first_deriv(const double param, CubitVector& FD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  FD.set(points[1].x() - points[0].x(), points[1].y() - points[0].y(),
	 points[1].z() - points[0].z());

}

//void GRLine::
//first_deriv(const CubitVector& coord,
//	    CubitVector& FD) const {
//
//  first_deriv(param_at_coord(coord), FD);
//
//}

void
GRLine::unit_normal(const double param, CubitVector& normal) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  if (n_dim == 3) {
    vFatalError("The normal to a curve is only defined in 2D",
		"GRLine::vUnitNormal(const double, CubitVector&)");
    assert(0);
  }

  CubitVector V;
  V.set(points[1] - points[0]);
  V.length(1.);

  normal.set(-V.y(), V.x(), 0.);

}

//void GRLine::
//unit_normal(const CubitVector& coord,
//	    CubitVector& normal) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  //Normal is the same everywhere, parameter does not matter.
//  unit_normal(0.5, normal);
//
//}

void
GRLine::unit_tangent(const double param, CubitVector& tangent) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  tangent.set(points[1].x() - points[0].x(), points[1].y() - points[0].y(),
	      points[1].z() - points[0].z());

  tangent.normalize();

}

//void GRLine::
//unit_tangent(const CubitVector& coord,
//	     CubitVector& tangent) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  //Tangent is the same everywhere, parameter does not matter.
//  unit_tangent(0.5, tangent);
//
//}

double
GRLine::arc_length() const
{

  assert(isValid());

  CubitVector V = points[1] - points[0];
  return V.length();

}

double
GRLine::arc_length(const double param1, const double param2) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  CubitVector V1;
  CubitVector V2;

  coord_at_param(param1, V1);
  coord_at_param(param2, V2);

  V2 -= V1;

  return V2.length();

}

//double GRLine::
//arc_length(const CubitVector& coord1,
//	   const CubitVector& coord2) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  CubitVector V = coord2 - coord1;
//
//  return V.length();
//
//}

//void GRLine::
//center_point(CubitVector& mid_pt) const {
//
//  assert(isValid());
//  mid_pt = 0.5*(points[0] + points[1]);
//
//}

void
GRLine::center_point(const double param1, const double param2,
		     CubitVector& mid_pt) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  CubitVector V1, V2;
  coord_at_param(param1, V1);
  coord_at_param(param2, V2);

  V2 -= V1;
  V2.length(0.5 * V2.length());
  mid_pt = points[0] + V2;
}

//void GRLine::
//center_point(const CubitVector& coord1,
//	     const CubitVector& coord2,
//	     CubitVector& mid_pt) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  CubitVector V = coord2 - coord1;
//  V.length(0.5 * V.length());
//  mid_pt = points[0] + V;
//
//}

bool
GRLine::isValid() const
{
  assert(this);
  if (n_pts != 2) {
    logMessage(0, "ERROR: The GRLine object is not defined with two points");
    return false;
  }
  if (points == NULL) {
    logMessage(0, "ERROR: The GRLine points array is uninitialized");
    return false;
  }
  if (n_dim < 2 || n_dim > 3) {
    logMessage(0, "ERROR: The number of dimension n_dim must be 2 or 3");
    return false;
  }
  if (points[0].x() < -1.e200) {
    logMessage(0, "ERROR: The GRLine first point x-coordinate is invalid");
    return false;
  }
  if (points[0].y() < -1.e200) {
    logMessage(0, "ERROR: The GRLine first point y-coordinate is invalid");
    return false;
  }
  if (points[0].z() < -1.e200) {
    logMessage(0, "ERROR: The GRLine first point z-coordinate is invalid");
    return false;
  }
  if (points[1].x() < -1.e200) {
    logMessage(0, "ERROR: The GRLine second point x-coordinate is invalid");
    return false;
  }
  if (points[1].y() < -1.e200) {
    logMessage(0, "ERROR: The GRLine second point y-coordinate is invalid");
    return false;
  }
  if (points[1].z() < -1.e200) {
    logMessage(0, "ERROR: The GRLine second point z-coordinate is invalid");
    return false;
  }

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    logMessage(0, "ERROR: The GRLine two points are at the same coordinate");
    return false;
  }

  return true;

}
