#include "GR_GRCurveGeom.h"
#include "CubitBox.hpp"
#include <map>

GRSpline::GRSpline() :
    GRCurveGeom(), num_cubics(0), cubics(NULL), closed(false)
{

  n_pts = 0;
  points = NULL;
  n_dim = 2;

}

GRSpline::~GRSpline()
{

  int i;
  for (i = 0; i < num_cubics; i++)
    if (cubics[i])
      delete cubics[i];
  delete[] cubics;

}

CubitBox
GRSpline::bounding_box() const
{

  assert(isValid());

  DLIList<CubitVector*> int_extrema;
  interior_extrema(int_extrema);
  int_extrema.reset();

  CubitVector v1, v2;
  coord_at_param(min_param(), v1);
  coord_at_param(max_param(), v2);

  CubitVector v_min(std::min(v1.x(), v2.x()), std::min(v1.y(), v2.y()), 0.);
  CubitVector v_max(std::max(v1.x(), v2.x()), std::max(v1.y(), v2.y()), 0.);

  int i;
  for (i = 0; i < int_extrema.size(); i++) {
    CubitVector* v = int_extrema.next(i);

    v_min.set(std::min(v_min.x(), v->x()), std::min(v_min.y(), v->y()), 0.);
    v_max.set(std::max(v_max.x(), v->x()), std::max(v_max.y(), v->y()), 0.);
    delete v;
  }
  int_extrema.clean_out();

  return CubitBox(v_min, v_max);

}

double
GRSpline::param_at_coord(const CubitVector& coord) const
{

  assert(isValid());

  double param;
  bool on_curve = coord_on_curve(coord, &param);

  if (!on_curve) {
    assert(on_curve);
  }

  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  return param;

}

double
GRSpline::param_at_arc_length(const double param_start,
			      const double length) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param_start, min_param()) >= 0
	  && iFuzzyComp(param_start, max_param()) <= 0);
  assert(iFuzzyComp(length, 0.) >= 0);
  assert(iFuzzyComp(length, arc_length(param_start, max_param())) <= 0);

  int i = static_cast<int>(param_start);
  if (i == num_cubics)
    i--;
  double this_length = 0., prev_length;
  double this_param = param_start - static_cast<double>(i);

  while (1) {
    prev_length = this_length;
    this_length += cubics[i]->arc_length(this_param, cubics[i]->max_param());

    if (this_length >= length)
      break;
    if (++i == num_cubics)
      return max_param();

    this_param = cubics[i]->min_param();
  }

  double param_lo = this_param;
  double param_up = cubics[i]->max_param();
  double param_new = 0.5 * (param_lo + param_up);
  double target_length = length - prev_length;
  double current_length = cubics[i]->arc_length(param_lo, param_new);

  int num_iters = 0;

  while (fabs(target_length - current_length) > cg_tolerance && ++num_iters < 75) {

    if (current_length <= target_length) {
      param_lo = param_new;
    }
    else {
      param_up = param_new;
    }

    param_new = 0.5 * (param_lo + param_up);
    current_length = cubics[i]->arc_length(this_param, param_new);

  }

  double return_param = param_new + static_cast<double>(i);
  assert(param_start < return_param);

  return return_param;

}

void
GRSpline::coord_at_param(const double param, CubitVector& coord) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double param_on_cubic;
  if (i == num_cubics)
    param_on_cubic = cubics[--i]->max_param();
  else
    param_on_cubic = param - static_cast<double>(i);

  cubics[i]->coord_at_param(param_on_cubic, coord);

}

double
GRSpline::coord_at_dist(const double target_dist, const bool start,
			CubitVector& coord) const
{

  assert(isValid());
  assert(iFuzzyComp(target_dist, 0.) >= 0);
  assert(target_dist <= arc_length()); //necessary but not sufficient for validity.

  double current_dist = LARGE_DBL;
  CubitVector end_coord, dum_coord;

  char string1[] =
      "The distance specified is greater than any distance on the curve";
  char string2[] =
      "GRSpline::coord_at_dist(const double, const bool, CubitVector&)";

  int i;
  if (start) {
    i = -1;
    coord_at_param(min_param(), end_coord);
    while (++i != num_cubics) {
      cubics[i]->coord_at_param(cubics[i]->max_param(), dum_coord);
      current_dist = (dum_coord - end_coord).length();
      if (current_dist > target_dist)
	break;
    }

    if (i == num_cubics)
      vFatalError(string1, string2);
  }
  else {
    i = num_cubics;
    coord_at_param(max_param(), end_coord);
    while (--i != -1) {
      cubics[i]->coord_at_param(cubics[i]->min_param(), dum_coord);
      current_dist = (dum_coord - end_coord).length();
      if (current_dist > target_dist)
	break;
    }

    if (i == -1)
      vFatalError(string1, string2);
  }

  //We found in which cubic curve the target_coord is located.
  //Now find exactly where.

  double param_lo = cubics[i]->min_param();
  double param_up = cubics[i]->max_param();
  double param_new = 0.5 * (param_lo + param_up);

  cubics[i]->coord_at_param(param_new, coord);
  current_dist = (coord - end_coord).length();

  while (fabs(current_dist - target_dist) > cg_tolerance) {

    if (current_dist <= target_dist) {
      if (start)
	param_lo = param_new;
      else
	param_up = param_new;
    }
    else {
      if (start)
	param_up = param_new;
      else
	param_lo = param_new;
    }

    param_new = 0.5 * (param_lo + param_up);
    cubics[i]->coord_at_param(param_new, coord);
    current_dist = (coord - end_coord).length();

  }

  return (static_cast<double>(i) + param_new);

}

double
GRSpline::coord_at_mid_dist(const double param1, const double param2,
			    CubitVector& middle) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);
  assert(iFuzzyComp(param1, param2) <= 0);

  double param_lo = param1;
  double param_up = param2;

  CubitVector coord1, coord2;
  coord_at_param(param1, coord1);
  coord_at_param(param2, coord2);

  double param = 0.5 * (param_lo + param_up);
  coord_at_param(param, middle);
  double dist1 = (coord1 - middle).length();
  double dist2 = (coord2 - middle).length();

  while (fabs(dist1 - dist2) > cg_tolerance) {

    if (dist1 > dist2)
      param_up = param;
    else
      param_lo = param;

    param = 0.5 * (param_lo + param_up);
    coord_at_param(param, middle);
    dist1 = (coord1 - middle).length();
    dist2 = (coord2 - middle).length();

  }

  return param;

}

//double GRSpline::
//coord_at_mid_dist(const CubitVector& coord1,
//		  const CubitVector& coord2,
//		  CubitVector& middle) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return coord_at_mid_dist(param_at_coord(coord1),
//			   param_at_coord(coord2),
//			   middle);
//
//}

void
GRSpline::closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
				 const double param1, const double param2) const
{

  assert(isValid());

  double param = closest_param_on_curve(coord, param1, param2);
  coord_at_param(param, closest);

  //assert(coord_on_curve(closest));

}

double
GRSpline::closest_param_on_curve(const CubitVector& coord, const double param1,
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
    p1 = param1;
    p2 = param2;
  }

  assert(iFuzzyComp(p1, min_param()) >= 0 && iFuzzyComp(p1, max_param()) <= 0);
  assert(iFuzzyComp(p2, min_param()) >= 0 && iFuzzyComp(p2, max_param()) <= 0);
  assert(iFuzzyComp(p1, p2) <= 0);

  //Determining which cubic curves contain
  //param1 and param2.
  int beg_cubic = static_cast<int>(p1);
  int end_cubic = static_cast<int>(p2);

  //If p2 == max_param, then we have a special case.
  if (end_cubic == num_cubics) {
    if (beg_cubic == end_cubic)
      return max_param();
    else
      --end_cubic;
  }

  //Initializing min_dist and return_param
  double return_param, min_dist;

  CubitVector coord1, coord2;
  coord_at_param(p1, coord1);
  coord_at_param(p2, coord2);
  double length1 = (coord - coord1).length();
  double length2 = (coord - coord2).length();
  if (length1 < length2) {
    return_param = p1;
    min_dist = length1;
  }
  else {
    return_param = p2;
    min_dist = length2;
  }

  //Testing first cubic in the range using lower param bound
  GRCubic* cub = cubics[beg_cubic];
  double param = static_cast<double>(beg_cubic);
  double beg_param = p1 - param;
  double end_param = p2 - param;

  if (beg_cubic == end_cubic)
    param += cub->closest_param_on_curve(coord, beg_param, end_param);
  else
    param += cub->closest_param_on_curve(coord, beg_param, cub->max_param());

  coord_at_param(param, coord1);
  length1 = (coord - coord1).length();
  if (length1 < min_dist) {
    return_param = param;
    min_dist = length1;
  }

  //If there is only one cubic to check, return here
  if (beg_cubic == end_cubic)
    return return_param;

  //Testing last cubic in the range using upper param bound
  cub = cubics[end_cubic];
  param = static_cast<double>(end_cubic);
  end_param = p2 - param;
  param += cub->closest_param_on_curve(coord, cub->min_param(), end_param);
  coord_at_param(param, coord1);
  length1 = (coord - coord1).length();
  if (length1 < min_dist) {
    return_param = param;
    min_dist = length1;
  }

  //Finally, checking the rest of the cubic curves.
  for (int i = beg_cubic + 1; i < end_cubic; i++) {
    cub = cubics[i];
    param = static_cast<double>(i);
    param += cub->closest_param_on_curve(coord, cub->min_param(),
					 cub->max_param());
    coord_at_param(param, coord1);
    length1 = (coord - coord1).length();
    if (length1 < min_dist) {
      return_param = param;
      min_dist = length1;
    }
  }

  assert(
      iFuzzyComp(return_param, p1) >= 0 && iFuzzyComp(return_param, p2) <= 0);

  return return_param;

}

bool
GRSpline::coord_on_curve(const CubitVector& coord, double* const param_X,
			 double* const param_Y, double* const param_Z) const
{

  assert(isValid());

  int i;
  for (i = 0; i < num_cubics; i++) {
    if (cubics[i]->coord_on_curve(coord, param_X, param_Y, param_Z)) {
      if (param_X != NULL) {
	*param_X += static_cast<double>(i);
      }
      if (param_Y != NULL) {
	*param_Y += static_cast<double>(i);
      }
      if (param_Z != NULL) {
	*param_Z = - LARGE_DBL;
      }
      return true;
    }
  }

  return false;

}

void
GRSpline::first_deriv(const double param, CubitVector& FD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double dum_param;
  if (i == num_cubics) {
    i--;
    dum_param = cubics[i]->max_param();
  }
  else {
    dum_param = param - static_cast<double>(i);
  }

  cubics[i]->first_deriv(dum_param, FD);

}

//void GRSpline::
//first_deriv(const CubitVector& coord,
//	    CubitVector& FD) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  first_deriv(param, FD);
//
//}

void
GRSpline::second_deriv(const double param, CubitVector& SD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double dum_param;
  if (i == num_cubics) {
    i--;
    dum_param = cubics[i]->max_param();
  }
  else {
    dum_param = param - static_cast<double>(i);
  }

  cubics[i]->second_deriv(dum_param, SD);

}
//
//void GRSpline::
//second_deriv(const CubitVector& coord,
//	     CubitVector& SD) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  second_deriv(param, SD);
//
//}
//
double
GRSpline::curvature(const double param) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double dum_param;
  if (i == num_cubics) {
    i--;
    dum_param = cubics[i]->max_param();
  }
  else {
    dum_param = param - static_cast<double>(i);
  }

  return cubics[i]->curvature(dum_param);

}

//double GRSpline::
//curvature(const CubitVector& coord) const {
//
//  assert(qValid());
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  return curvature(param);
//
//}

void
GRSpline::unit_normal(const double param, CubitVector& normal) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double dum_param;
  if (i == num_cubics) {
    i--;
    dum_param = cubics[i]->max_param();
  }
  else {
    dum_param = param - static_cast<double>(i);
  }

  cubics[i]->unit_normal(dum_param, normal);

}

//void GRSpline::
//unit_normal(const CubitVector& coord,
//	    CubitVector& normal) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  unit_normal(param, normal);
//
//}
//
void
GRSpline::unit_tangent(const double param, CubitVector& tangent) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  int i = static_cast<int>(param);
  double dum_param;
  if (i == num_cubics) {
    i--;
    dum_param = cubics[i]->max_param();
  }
  else {
    dum_param = param - static_cast<double>(i);
  }

  cubics[i]->unit_tangent(dum_param, tangent);

}

//void GRSpline::
//unit_tangent(const CubitVector& coord,
//	     CubitVector& tangent) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  unit_tangent(param, tangent);
//
//}

double
GRSpline::arc_length() const
{

  return arc_length(min_param(), max_param());

}

double
GRSpline::arc_length(const double param1, const double param2) const
{

  assert(isValid());
  //assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  double length = 0;

  int i = static_cast<int>(param1);
  int j = static_cast<int>(param2);
  double dum_param1 = param1 - static_cast<double>(i);
  double dum_param2 = param2 - static_cast<double>(j);

  if (i == num_cubics) {
    i--;
    dum_param1 = cubics[i]->max_param();
  }
  if (j == num_cubics) {
    j--;
    dum_param2 = cubics[j]->max_param();
  }

  assert(i >= 0 && i < num_cubics);
  assert(j >= 0 && j < num_cubics);
  assert(i <= j);
  assert(
      iFuzzyComp(dum_param1, cubics[i]->min_param()) >= 0
	  && iFuzzyComp(dum_param2, cubics[i]->max_param()) <= 0);
  assert(
      iFuzzyComp(dum_param1, cubics[j]->min_param()) >= 0
	  && iFuzzyComp(dum_param2, cubics[j]->max_param()) <= 0);

  if (i == j) {
    length += cubics[i]->arc_length(dum_param1, dum_param2);
  }
  else {
    length += cubics[i]->arc_length(dum_param1, cubics[i]->max_param());
    length += cubics[j]->arc_length(cubics[j]->min_param(), dum_param2);

    int k;
    for (k = i + 1; k < j; k++)
      length += cubics[k]->arc_length();

  }

  return length;

}

//double GRSpline::
//arc_length(const CubitVector& coord1,
//	   const CubitVector& coord2) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//  return arc_length(param1, param2);
//
//}

//void GRSpline::
//center_point(CubitVector& mid_pt) const {
//
//  assert(isValid());
//
//  center_point(min_param(), max_param(), mid_pt);
//
//}

void
GRSpline::center_point(const double param1, const double param2,
		       CubitVector& mid_pt) const
{

  assert(isValid());
  // assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  double param_lo = min(param1, param2);
  double param_up = max(param1, param2);
  // double total_length = arc_length(param1, param2);

  double total_length = arc_length(param_lo, param_up);
  double target_length = 0.5 * total_length;

//  double param_lo = param1;
//  double param_up = param2;
  double param_new = 0.5 * (param_lo + param_up);
  //total_length = arc_length(param1, param_new);
  total_length = arc_length(param_lo, param_new);
  double paramTemp = param_lo;
  while (fabs(total_length - target_length) > cg_tolerance) {

    if (total_length < target_length)
      param_lo = param_new;
    else
      param_up = param_new;

    param_new = 0.5 * (param_lo + param_up);
    // total_length = arc_length(param1, param_new);
    total_length = arc_length(paramTemp, param_new);

  }

  coord_at_param(param_new, mid_pt);

}

//void GRSpline::
//center_point(const CubitVector& coord1,
//	     const CubitVector& coord2,
//	     CubitVector& mid_pt) const {
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//  center_point(param1, param2, mid_pt);
//
//}
//
//double GRSpline::
//TVT() const {
//
//  assert(isValid());
//  return TVT(min_param(), max_param());
//
//}

double
GRSpline::TVT(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  double tot_var = 0.;

  int i = static_cast<int>(param1);
  int j = static_cast<int>(param2);
  double dum_param1 = param1 - static_cast<double>(i);
  double dum_param2 = param2 - static_cast<double>(j);

  if (i == num_cubics) {
    i--;
    dum_param1 = cubics[i]->max_param();
  }
  if (j == num_cubics) {
    j--;
    dum_param2 = cubics[j]->max_param();
  }

  assert(i >= 0 && i < num_cubics);
  assert(j >= 0 && j < num_cubics);
  assert(i <= j);
  assert(
      iFuzzyComp(dum_param1, cubics[i]->min_param()) >= 0
	  && iFuzzyComp(dum_param2, cubics[i]->max_param()) <= 0);
  assert(
      iFuzzyComp(dum_param1, cubics[j]->min_param()) >= 0
	  && iFuzzyComp(dum_param2, cubics[j]->max_param()) <= 0);

  if (i == j)
    tot_var += cubics[i]->TVT(dum_param1, dum_param2);
  else {
    tot_var += cubics[i]->TVT(dum_param1, cubics[i]->max_param());
    tot_var += cubics[j]->TVT(cubics[j]->min_param(), dum_param2);

    int k;
    for (k = i + 1; k < j; k++)
      tot_var += cubics[k]->TVT(cubics[k]->min_param(), cubics[k]->max_param());

  }

  return tot_var;

}

//double GRSpline::
//TVT(const CubitVector& coord1,
//    const CubitVector& coord2) const {
//
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//  return TVT(param1, param2);
//
//}

double
GRSpline::mid_TVT(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  //All the critical points located between param1 and param2
  //are first put into a set. This set will be used to compute
  //the TVT between param1 and param2.

  int beg_cubic = static_cast<int>(param1);
  int end_cubic = static_cast<int>(param2);
  if (beg_cubic == num_cubics)
    --beg_cubic;
  if (end_cubic == num_cubics)
    --end_cubic;
  assert(beg_cubic <= end_cubic);

  std::set<double> cubic_crit_pts;
  for (int i = beg_cubic; i < end_cubic + 1; i++)
    cubics[i]->get_crit_pts(cubic_crit_pts, static_cast<double>(i));

  std::set<double>::iterator itL = cubic_crit_pts.upper_bound(param1);
  std::set<double>::iterator itU = cubic_crit_pts.upper_bound(param2);

  std::set<double, param_compare> critical_pts;
  std::copy(itL, itU, std::inserter(critical_pts, critical_pts.begin()));
  critical_pts.insert(param1);
  critical_pts.insert(param2);

  assert(iFuzzyComp(*critical_pts.begin(), param1) == 0);
  assert(iFuzzyComp(*(--critical_pts.end()), param2) == 0);

  //The TVT between param1 and param2 is computed. The accumulated
  //TVT at each critical point is stored into a map (along with the
  //critical point parameter).

  double total_TVT = 0.;

  std::map<double, double> TVT_map;
  TVT_map.insert(std::make_pair(total_TVT, *critical_pts.begin()));

  // Klocwork complains, but loop is never entered if itL == end().
  for (itL = critical_pts.begin(); itL != (--critical_pts.end());) {
    double p1 = *itL, p2 = *(++itL);
    total_TVT += TVT(p1, p2);
    TVT_map.insert(std::make_pair(total_TVT, *itL));
  }

  if (fabs(total_TVT) < 1.e-3)
    return 0.5 * (param1 + param2);

  //Once the TVT is know on the interval, compute the target TVT
  //and find between which critical point pair this target TVT
  //is found.

  double target_TVT = 0.5 * total_TVT;
  std::map<double, double>::iterator itm = TVT_map.upper_bound(target_TVT);

  assert(itm != TVT_map.end() && itm != TVT_map.begin());
  double TVT_lo = (--itm)->first;
  double param_lo = itm->second;
  double TVT_up = (++itm)->first;
  double param_up = itm->second;

  assert(iFuzzyComp(TVT_lo, target_TVT) <= 0);
  assert(iFuzzyComp(TVT_up, target_TVT) >= 0);

  //Compute the tangent at these lower and upper critical points.

  CubitVector tan_lo, tan_up;
  unit_tangent(param_lo, tan_lo);
  unit_tangent(param_up, tan_up);

  //Finding the angle needed to achieve the target TVT. Use the cross
  //product to find if this angle must be computed CW or CCW.

  double cross = tan_lo.x() * tan_up.y() - tan_lo.y() * tan_up.x();
  double angle_rotate = target_TVT - TVT_lo;
  assert(iFuzzyComp(angle_rotate, 0.) >= 0);
  if (iFuzzyComp(cross, 0.) < 0)
    angle_rotate *= -1.;

  //Computing the target tangent using a 2D rotation matrix

  double cos_rotate = cos(angle_rotate);
  double sin_rotate = sin(angle_rotate);
  CubitVector tan_target(tan_lo.x() * cos_rotate - tan_lo.y() * sin_rotate,
			 tan_lo.y() * cos_rotate + tan_lo.x() * sin_rotate, 0.);
  assert(iFuzzyComp(tan_target.length(), 1.) == 0);

  int i = static_cast<int>(param_lo);

#ifndef NDEBUG

  int j = static_cast<int>(param_up);

  if (i != j) {
    assert(j == i + 1);
    if (j < num_cubics)
      assert(
	  iFuzzyComp(param_up - static_cast<double>(j), cubics[j]->min_param())
	      == 0);
    else if (j == num_cubics)
      assert(iFuzzyComp(param_up, max_param()) == 0);
    else
      assert(0);
  }

#endif

  //Compute the coefficient at which the tangent vector is 
  //parallel to the unit tangent previously computed.

  int num_roots;
  double roots[2];
  const double* const x_coeffs = cubics[i]->get_x_coeff();
  const double* const y_coeffs = cubics[i]->get_y_coeff();

  double a = 3. * (tan_target.x() * y_coeffs[3] - tan_target.y() * x_coeffs[3]);
  double b = 2. * (tan_target.x() * y_coeffs[2] - tan_target.y() * x_coeffs[2]);
  double c = tan_target.x() * y_coeffs[1] - tan_target.y() * x_coeffs[1];

  solveQuadratic(a, b, c, &num_roots, &roots[0], &roots[1]);

  double ret_param = -LARGE_DBL;
  for (int k = 0; k < num_roots; k++) {
    roots[k] += static_cast<double>(i);
    if (iFuzzyComp(roots[k], param_lo) >= 0
	&& iFuzzyComp(roots[k], param_up) <= 0) {
      ret_param = roots[k];
      break;
    }
  }

  if (iFuzzyComp(ret_param, min_param()) < 0
      || iFuzzyComp(ret_param, max_param()) > 0)
    ret_param = cubics[i]->mid_TVT(param_lo, param_up);

  assert(
      iFuzzyComp(ret_param, min_param()) >= 0
	  && iFuzzyComp(ret_param, max_param()) <= 0);
  //assert( fabs( (TVT(param1, ret_param) / TVT(param1, param2)) - 0.5 ) < 1.e-5 );

  return ret_param;

}

//double GRSpline::
//mid_TVT(const CubitVector& coord1,
//	const CubitVector& coord2) const {
//
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//  return mid_TVT(param1, param2);
//
//}

void
GRSpline::interior_extrema(DLIList<CubitVector*>& point_list) const
{

  int i;
  for (i = 0; i < num_cubics; i++) {
    DLIList<CubitVector*> this_cubic;
    cubics[i]->interior_extrema(point_list);
    point_list += this_cubic;
  }
}
bool
GRSpline::line_intersect(const double coords1[3], const double coords2[3],
			 const double begParam, const double endParam,
			 std::vector<double> & params) const
{
  std::vector<double> parameters;
  bool qRetVal = false;
  for (int i = 0; i < num_cubics; i++) {

    qRetVal = cubics[i]->line_intersect(coords1, coords2,
					begParam - static_cast<double>(i),
					endParam - static_cast<double>(i),
					parameters) || qRetVal;
    for (int j = 0; j < static_cast<int>(parameters.size()); j++) {
      parameters[j] = parameters[j] + static_cast<double>(i);
      params.push_back(parameters[j]);
    }
    parameters.clear();
  }
  return qRetVal;

}
