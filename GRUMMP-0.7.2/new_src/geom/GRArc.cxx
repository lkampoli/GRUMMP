#include "GR_GRCurveGeom.h"
#include "CubitBox.hpp"
#include <iostream>

GRArc::GRArc() :
    GRCurveGeom(), radius(-1.), angle_beg(-1.), angle_span(-1.)
{
  n_pts = 3;
  points = new CubitVector[n_pts];
  n_dim = 2;

  points[0].set(-LARGE_DBL, -LARGE_DBL, 0.);
  points[1].set(-LARGE_DBL, -LARGE_DBL, 0.);
  points[2].set(-LARGE_DBL, -LARGE_DBL, 0.);

}

GRArc::GRArc(const CubitVector& center, const CubitVector& beg_pt,
	     const CubitVector& end_pt) :
    GRCurveGeom(), radius(-1.), angle_beg(-1.), angle_span(-1.)
{
  n_pts = 3;
  points = new CubitVector[n_pts];
  n_dim = 2;

  points[0] = center;
  points[1] = beg_pt;
  points[2] = end_pt;

  int i;
  for (i = 0; i < n_pts; i++)
    points[i].z(0.);

  radius = (points[1] - points[0]).length();
  if (iFuzzyComp((points[2] - points[0]).length(), radius) != 0)
    vFatalError(
	"Arc start and end points are not at the same distance from the center.",
	"GRArc::GRArc(const CubitVector&, const CubitVector&, const CubitVector&)");

  if (iFuzzyComp(points[1].x(), points[2].x()) == 0
      && iFuzzyComp(points[1].y(), points[2].y()) == 0)
    vFatalError(
	"Arc start and end point are the same. Use circle instead.",
	"GRArc::GRArc(const CubitVector&, const CubitVector&, const CubitVector&)");

  double angle_end;
  angle_beg = atan2(points[1].y() - points[0].y(),
		    points[1].x() - points[0].x());
  angle_end = atan2(points[2].y() - points[0].y(),
		    points[2].x() - points[0].x());

  if (iFuzzyComp(angle_beg, 0.) == 0)
    angle_beg = 0.;
  if (iFuzzyComp(angle_beg, 0.) < 0)
    angle_beg += 2. * M_PI;
  if (iFuzzyComp(angle_end, 0.) < 0)
    angle_end += 2. * M_PI;

  if (iFuzzyComp(angle_end, angle_beg) == 1)
    angle_span = angle_end - angle_beg;
  else
    angle_span = 2. * M_PI - angle_beg + angle_end;

  assert(isValid());

}

GRArc::GRArc(const CubitVector& beg_pt, const CubitVector& end_pt,
	     const double rad, const bool short_arc) :
    GRCurveGeom(), radius(rad), angle_beg(-1.), angle_span(-1.)
{

  n_pts = 3;
  n_dim = 2;
  points = new CubitVector[n_pts];

  if (iFuzzyComp(radius, 0.) <= 0) {
    vFatalError("Arc has an invalid radius.",
		"GRArc::GRArc(CubitVector&, CubitVector&, double,  bool)");
  }

  points[1] = beg_pt;
  points[1].z(0.);
  points[2] = end_pt;
  points[2].z(0.);

  if (iFuzzyComp(points[1].x(), points[2].x()) == 0
      && iFuzzyComp(points[1].y(), points[2].y()) == 0) {
    vFatalError(
	"Both input points are the same. Cannot determine on which side to put the center.",
	"GRArc::GRArc(CubitVector&, CubitVector&, double,  bool)");
  }

  CubitVector V;
  if (short_arc)
    V = points[2] - points[1];
  else
    V = points[1] - points[2];

  double L = V.length();
  double D = (radius * radius) - (0.25 * L * L);
  if (D < 0.) {
    if (fabs(D) < 1.e-8)
      D = 0.;
    else
      vFatalError("Unable to determine arc center. Error in boundary file?",
		  "GRArc::GRArc(CubitVector&, CubitVector&, double,  bool)");
  }
  else
    D = sqrt(D);

  V /= L;

  if (short_arc)
    points[0].set(points[1].x() + 0.5 * L * V.x() - D * V.y(),
		  points[1].y() + 0.5 * L * V.y() + D * V.x(), 0.);
  else
    points[0].set(points[2].x() + 0.5 * L * V.x() - D * V.y(),
		  points[2].y() + 0.5 * L * V.y() + D * V.x(), 0.);

  double angle_end;
  angle_beg = atan2(points[1].y() - points[0].y(),
		    points[1].x() - points[0].x());
  angle_end = atan2(points[2].y() - points[0].y(),
		    points[2].x() - points[0].x());

  //This seems redundant, but I am having problem with floating point precision here.
  if (iFuzzyComp(angle_beg, 0.) == 0
      || iFuzzyComp(angle_beg + 2. * M_PI, 2. * M_PI) == 0)
    angle_beg = 0.;

  if (iFuzzyComp(angle_beg, 0.) < 0)
    angle_beg += 2. * M_PI;
  if (iFuzzyComp(angle_end, 0.) < 0)
    angle_end += 2. * M_PI;

  if (iFuzzyComp(angle_end, angle_beg) == 1)
    angle_span = angle_end - angle_beg;
  else
    angle_span = 2. * M_PI - angle_beg + angle_end;

#ifndef NDEBUG

  //printf("Angle span = %lf\n", angle_span);

  if (short_arc)
    assert(iFuzzyComp(angle_span, 0.) > 0 && iFuzzyComp(angle_span, M_PI) <= 0);
  else
    assert(
	iFuzzyComp(angle_span, M_PI) > 0 && iFuzzyComp(angle_span, 2. * M_PI) < 0);
#endif

  assert(isValid());

}

GRArc::GRArc(const GRArc& arc) :
    GRCurveGeom(), radius(arc.radius), angle_beg(arc.angle_beg), angle_span(
	arc.angle_span)
{
  assert(arc.isValid());

  n_pts = arc.n_pts;
  points = new CubitVector[n_pts];
  assert(n_pts == 3);
  n_dim = arc.n_dim;

  points[0] = arc.points[0];
  points[1] = arc.points[1];
  points[2] = arc.points[2];

  assert(isValid());

}

GRArc&
GRArc::operator=(const GRArc& arc)
{

  assert(arc.isValid());

  n_pts = arc.n_pts;
  assert(n_pts == 3);
  n_dim = arc.n_dim;

  points[0] = arc.points[0];
  points[1] = arc.points[1];
  points[2] = arc.points[2];

  radius = arc.radius;
  angle_beg = arc.angle_beg;
  angle_span = arc.angle_span;

  assert(isValid());

  return *this;

}

GRArc::~GRArc()
{
}

CubitBox
GRArc::bounding_box() const
{

  assert(isValid());

  DLIList<CubitVector*> coord_list;
  interior_extrema(coord_list);
  coord_list.reset();

  double x_max = std::max(points[1].x(), points[2].x());
  double x_min = std::min(points[1].x(), points[2].x());
  double y_max = std::max(points[1].y(), points[2].y());
  double y_min = std::min(points[1].y(), points[2].y());
  CubitVector v_min, v_max;

  int i;
  for (i = 0; i < coord_list.size(); i++) {
    CubitVector* coord = coord_list.next(i);
    x_max = std::max(x_max, coord->x());
    x_min = std::min(x_min, coord->x());
    y_max = std::max(y_max, coord->y());
    y_min = std::min(y_min, coord->y());
    delete coord;
  }
  coord_list.clean_out();

  v_min.set(x_min, y_min, 0.);
  v_max.set(x_max, y_max, 0.);

  return CubitBox(v_min, v_max);

}

double
GRArc::param_at_coord(const CubitVector& coord) const
{

  assert(isValid());

  double param;
  bool on_curve = coord_on_curve(coord, &param);
  assert(on_curve);

  return param;

}

double
GRArc::param_at_arc_length(const double param_start, const double length) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param_start, min_param()) >= 0
	  && iFuzzyComp(param_start, max_param()) <= 0);
  assert(iFuzzyComp(length, 0.) >= 0);
  assert(iFuzzyComp(length, arc_length(param_start, max_param())) <= 0);

  double this_angle_span = length / radius;
  return (param_start + this_angle_span / angle_span);

}

void
GRArc::coord_at_param(const double param, CubitVector& coord) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  double arg = angle_beg + param * angle_span;
  coord.set(points[0].x() + radius * cos(arg),
	    points[0].y() + radius * sin(arg), 0.);

}

double
GRArc::coord_at_dist(const double dist, const bool start,
		     CubitVector& coord) const
{

  assert(isValid());
  assert(iFuzzyComp(dist, 0.) >= 0);
  assert(iFuzzyComp(dist, arc_length()) <= 0);

  //this_angle_span is the angle subtended by a chord
  //of length dist on the arc. Simply obtained using
  //law of cosines.

  double arg = 1. - 0.5 * ((dist * dist) / (radius * radius));
  assert(iFuzzyComp(arg, 0.) > 0 && iFuzzyComp(arg, M_PI) < 0);
  double this_angle_span = acos(arg);

  double span_from_beg = LARGE_DBL;

  if (start)
    span_from_beg = this_angle_span;
  else
    span_from_beg = angle_span - this_angle_span;

  arg = angle_beg + span_from_beg;

  coord.set(points[0].x() + radius * cos(arg),
	    points[0].y() + radius * sin(arg), 0.);

  return (span_from_beg / angle_span);

}

double
GRArc::coord_at_mid_dist(const double param1, const double param2,
			 CubitVector& middle) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);
  assert(iFuzzyComp(param1, param2) < 0);

  double param = 0.5 * (param1 + param2);
  coord_at_param(param, middle);

  return param;

}

//double GRArc::
//coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
//		  CubitVector& middle) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//
//  assert(iFuzzyComp(param1, param2) < 0);
//
//  double param  = 0.5 * (param1 + param2);
//  coord_at_param(param, middle);
//
//  return param;
//
//}

inline void
GRArc::closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			      const double param1, const double param2) const
{

  assert(isValid());

  coord_at_param(closest_param_on_curve(coord, param1, param2), closest);

}

double
GRArc::closest_param_on_curve(const CubitVector& coord, const double param1,
			      const double param2) const
{

  assert(isValid());

  double p1, p2;
  if (iFuzzyComp(param1, min_param()) < 0
      || iFuzzyComp(param2, min_param()) < 0) {
    p1 = min_param();
    p2 = max_param();
  }
  else {
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

  //If coord is the center point of the arc:
  if (iFuzzyComp(points[0].x(), coord.x()) == 0
      && iFuzzyComp(points[0].y(), coord.y()) == 0)
    return (0.5 * (p1 + p2));

  double angle = atan2(coord.y() - points[0].y(), coord.x() - points[0].x());
  if (angle < 0.)
    angle += 2. * M_PI;

  assert(iFuzzyComp(angle, 0.) >= 0 && iFuzzyComp(angle, 2. * M_PI) <= 0);

  double angle_end = angle_beg + angle_span;

  if (iFuzzyComp(angle_end, 2 * M_PI) <= 0) {
    if (iFuzzyComp(angle, angle_end) <= 0
	&& iFuzzyComp(angle, angle_beg) >= 0) {
      double p = (angle - angle_beg) / angle_span;
      if (p < p1)
	return p1;
      if (p > p2)
	return p2;
      return p;
    }
  }
  else {
    if (iFuzzyComp(angle, angle_beg) < 0)
      angle += 2 * M_PI;
    if (iFuzzyComp(angle, angle_end) <= 0
	&& iFuzzyComp(angle, angle_beg) >= 0) {
      double p = (angle - angle_beg) / angle_span;
      if (p < p1)
	return p1;
      if (p > p2)
	return p2;
      return p;
    }
  }

  double angle_avg = 0.5 * (angle_beg + angle_end) + M_PI;

  if (iFuzzyComp(angle, angle_end) == 1 && iFuzzyComp(angle, angle_avg) <= 0)
    return p2;
  else
    return p1;

  assert(0);
  return -LARGE_DBL;

}

bool
GRArc::coord_on_curve(const CubitVector& coord, double* const param_X,
		      double* const param_Y, double* const param_Z) const
{

  assert(isValid());

  double rad = (points[0] - coord).length();

  if (iFuzzyComp(rad, radius) != 0) {
    if (param_X != NULL)
      *param_X = - LARGE_DBL;
    if (param_Y != NULL)
      *param_Y = - LARGE_DBL;
    if (param_Z != NULL)
      *param_Z = - LARGE_DBL;
    return false;
  }

  double angle = atan2(coord.y() - points[0].y(), coord.x() - points[0].x());
  if (iFuzzyComp(angle, 0.) < 0)
    angle += 2. * M_PI;

  double angle_end = angle_beg + angle_span;

  //   printf("before angle = %lf, angle_beg = %lf, angle_end = %lf\n",
  // 	 angle, angle_beg, angle_end);

  if (iFuzzyComp(angle, 0.) >= 0 && iFuzzyComp(angle, angle_beg) < 0)
    angle += 2. * M_PI;

  //    printf("after angle = %lf, angle_beg = %lf, angle_end = %lf\n\n",
  //  	 angle, angle_beg, angle_end);

  assert(
      iFuzzyComp(angle, angle_beg) >= 0 && iFuzzyComp(angle, angle_beg + 2 * M_PI) <= 0);

  if (iFuzzyComp(angle, angle_beg) >= 0 && iFuzzyComp(angle, angle_end) <= 0) {
    if (param_X != NULL)
      *param_X = (angle - angle_beg) / angle_span;
    if (param_Y != NULL)
      *param_Y = (angle - angle_beg) / angle_span;
    if (param_Z != NULL)
      *param_Z = -LARGE_DBL;
    return true;
  }
  assert(iFuzzyComp(angle, angle_end) > 0);

  if (param_X != NULL)
    *param_X = - LARGE_DBL;
  if (param_Y != NULL)
    *param_Y = - LARGE_DBL;
  if (param_Z != NULL)
    *param_Z = - LARGE_DBL;
  return false;

}

void
GRArc::first_deriv(const double param, CubitVector& FD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  double arg = angle_beg + param * angle_span;
  FD.set(-radius * angle_span * sin(arg), radius * angle_span * cos(arg), 0.);

}

//inline void GRArc::
//first_deriv(const CubitVector& coord,
//	    CubitVector& FD) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  first_deriv(param_at_coord(coord), FD);
//
//}

void
GRArc::second_deriv(const double param, CubitVector& SD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  double span_sq = angle_span * angle_span;
  double arg = angle_beg + param * angle_span;
  SD.set(-radius * span_sq * cos(arg), -radius * span_sq * sin(arg), 0.);

}

//void GRArc::
//second_deriv(const CubitVector& coord,
//	     CubitVector& SD) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  second_deriv(param_at_coord(coord), SD);
//
//}

double
GRArc::curvature(const double param) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  return 1. / radius;

}

//double GRArc::
//curvature(const CubitVector& coord) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  return 1. / radius;
//
//}

void
GRArc::unit_normal(const double param, CubitVector& normal) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  CubitVector tangent;
  unit_tangent(param, tangent);

  normal.set(-tangent.y(), tangent.x(), 0.);

}

//void GRArc::
//unit_normal(const CubitVector& coord,
//	    CubitVector& normal) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  unit_normal(param_at_coord(coord), normal);
//
//}

void
GRArc::unit_tangent(const double param, CubitVector& tangent) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  first_deriv(param, tangent);
  tangent.normalize();

}

//void GRArc::
//unit_tangent(const CubitVector& coord,
//	     CubitVector& tangent) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  unit_tangent(param_at_coord(coord), tangent);
//
//}

double
GRArc::arc_length() const
{

  assert(isValid());
  return arc_length(min_param(), max_param());

}

double
GRArc::arc_length(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, 0.) >= 0 && iFuzzyComp(param1, 1.) <= 0);
  assert(iFuzzyComp(param2, 0.) >= 0 && iFuzzyComp(param2, 1.) <= 0);

  assert(iFuzzyComp(param1, param2) <= 0);
  if (iFuzzyComp(param1, param2) == 0)
    return 0.;
  return radius * (param2 - param1) * angle_span;

}

//double GRArc::
//arc_length(const CubitVector& coord1,
//	   const CubitVector& coord2) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return arc_length(param_at_coord(coord1),
//		    param_at_coord(coord2));
//
//}

//void GRArc::
//center_point(CubitVector& mid_pt) const {
//
//  assert(isValid());
//  double param_avg = 0.5 * (min_param() + max_param());
//  coord_at_param(param_avg, mid_pt);
//
//}

void
GRArc::center_point(const double param1, const double param2,
		    CubitVector& mid_pt) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, 0.) >= 0 && iFuzzyComp(param1, 1.) <= 0);
  assert(iFuzzyComp(param2, 0.) >= 0 && iFuzzyComp(param2, 1.) <= 0);
  //assert(iFuzzyComp(param1, param2) < 0);

  double param_avg = 0.5 * (param1 + param2);
  coord_at_param(param_avg, mid_pt);

}

//void GRArc::
//center_point(const CubitVector& coord1, const CubitVector& coord2,
//	     CubitVector& mid_pt) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param_avg = 0.5 * (param_at_coord(coord1) + param_at_coord(coord2));
//  coord_at_param(param_avg, mid_pt);
//
//}

//double GRArc::
//TVT() const {
//
//  assert(isValid());
//  return TVT(min_param(), max_param());
//
//}

double
GRArc::TVT(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, 0.) >= 0 && iFuzzyComp(param1, 1.) <= 0);
  assert(iFuzzyComp(param2, 0.) >= 0 && iFuzzyComp(param2, 1.) <= 0);
  assert(iFuzzyComp(param1, param2) < 0);

  return (param2 - param1) * angle_span;

}

//double GRArc::
//TVT(const CubitVector& coord1,
//    const CubitVector& coord2) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return TVT(param_at_coord(coord1), param_at_coord(coord2));
//
//}

double
GRArc::mid_TVT(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, 0.) >= 0 && iFuzzyComp(param1, 1.) <= 0);
  assert(iFuzzyComp(param2, 0.) >= 0 && iFuzzyComp(param2, 1.) <= 0);
  assert(iFuzzyComp(param1, param2) < 0);

  return 0.5 * (param1 + param2);

}

//double GRArc::
//mid_TVT(const CubitVector& coord1,
//	const CubitVector& coord2) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return 0.5 * (param_at_coord(coord1) + param_at_coord(coord2));
//
//}

void
GRArc::interior_extrema(DLIList<CubitVector*>& point_list) const
{

  assert(isValid());
  assert(point_list.size() == 0);
  assert(
      iFuzzyComp(angle_beg, 0.) >= 0 && iFuzzyComp(angle_beg, 2. * M_PI) < 0);

  double angle_end = angle_beg + angle_span;

  assert(
      iFuzzyComp(angle_end, 0.) >= 0 && iFuzzyComp(angle_end, 4. * M_PI) < 0);

  double angle = 0.5 * M_PI;
  while (iFuzzyComp(angle, 4. * M_PI) < 0) {
    if (iFuzzyComp(angle_beg, angle) < 0 && iFuzzyComp(angle_end, angle) > 0) {
      CubitVector* cv = new CubitVector(points[0].x() + radius * cos(angle),
					points[0].y() + radius * sin(angle),
					0.);
      point_list.insert(cv);
    }
    angle += 0.5 * M_PI;
  }

}

bool
GRArc::line_intersect(const double */*coords1[3]*/,
		      const double */*coords2[3]*/, const double /*begParam*/,
		      const double /*endParam*/,
		      std::vector<double> &/* params*/) const
{

  vFatalError("Line Intersection not set up for Arcs.\n",
	      "GRArc::line_intersect");
  return false;
}

double
GRArc::x(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].x();

}

double
GRArc::y(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].y();

}

double
GRArc::z(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].z();

}

void
GRArc::get_coord(const int pt, CubitVector& coord) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  coord = points[pt];

}

const CubitVector&
GRArc::get_coord(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  return points[pt];

}

void
GRArc::set_x(const int pt, const double XX)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  double change = XX - points[pt].x();
  points[pt].x(XX);

  switch (pt)
    {
    case 0:
      {
	points[1].x(points[1].x() + change);
	points[2].x(points[2].x() + change);
	break;
      }
    case 1:
      {
	points[0].x(points[0].x() + change);
	points[2].x(points[2].x() + change);
	break;
      }
    case 2:
      {
	points[0].x(points[0].x() + change);
	points[1].x(points[1].x() + change);
	break;
      }
    default:
      {
	assert(0);
	break;
      }
    }
}

void
GRArc::set_y(const int pt, const double YY)
{
  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  double change = YY - points[pt].y();
  points[pt].y(YY);

  switch (pt)
    {
    case 0:
      {
	points[1].y(points[1].y() + change);
	points[2].y(points[2].y() + change);
	break;
      }
    case 1:
      {
	points[0].y(points[0].y() + change);
	points[2].y(points[2].y() + change);
	break;
      }
    case 2:
      {
	points[0].y(points[0].y() + change);
	points[1].y(points[1].y() + change);
	break;
      }
    default:
      {
	assert(0);
	break;
      }
    }
}

void
GRArc::set_z(const int pt, const double ZZ)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  assert(0);
  if (iFuzzyComp(ZZ, 0.) != 0) {
  }

}

void
GRArc::set_coord(const int pt, const CubitVector& new_coord)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  CubitVector change = new_coord - points[pt];
  assert(iFuzzyComp(change.z(), 0.) == 0);

  points[0] += change;
  points[1] += change;
  points[2] += change;

}

bool
GRArc::isValid() const
{

  assert(this);

  if (n_pts != 3)
    return false;
  if (n_dim != 2)
    return false;
  if (points == NULL)
    return false;
  if (points[0].x() < -1.e200)
    return false;
  if (points[0].y() < -1.e200)
    return false;
  if (iFuzzyComp(points[0].z(), 0.) != 0)
    return false;
  if (points[1].x() < -1.e200)
    return false;
  if (points[1].y() < -1.e200)
    return false;
  if (iFuzzyComp(points[1].z(), 0.) != 0)
    return false;
  if (points[2].x() < -1.e200)
    return false;
  if (points[2].y() < -1.e200)
    return false;
  if (iFuzzyComp(points[2].z(), 0.) != 0)
    return false;
  if (radius <= 0.)
    return false;
  if (iFuzzyComp((points[1] - points[0]).length(), radius) != 0
      && iFuzzyComp((points[2] - points[0]).length(), radius) != 0)
    return false;
  if (iFuzzyComp(angle_beg, 0.) == -1 || iFuzzyComp(angle_beg, 2 * M_PI) >= 0)
    return false;
  if (iFuzzyComp(angle_span, 0.) <= 0)
    return false;
  if (iFuzzyComp(angle_span, 2. * M_PI) == 1)
    return false;

  return true;

}
