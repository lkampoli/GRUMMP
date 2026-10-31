#include "GR_GRCurveGeom.h"
#include "CubitBox.hpp"
#include <iostream>

GRCubic::GRCubic() :
    GRCurveGeom(), x_coeff(NULL), y_coeff(NULL), crit_pts()
{

  n_pts = 0;
  points = NULL;
  n_dim = 2;

  x_coeff = new double[4];
  y_coeff = new double[4];

  unsigned int i;
  for (i = 0; i < 4; i++)
    x_coeff[i] = y_coeff[i] = -LARGE_DBL;

}

GRCubic::GRCubic(const double x_coefficient[], const double y_coefficient[]) :
    GRCurveGeom(), x_coeff(NULL), y_coeff(NULL), crit_pts()
{

  n_pts = 0;
  points = NULL;
  n_dim = 2;

  x_coeff = new double[4];
  y_coeff = new double[4];

  unsigned int i;
  for (i = 0; i < 4; i++) {
    x_coeff[i] = x_coefficient[i];
    y_coeff[i] = y_coefficient[i];
  }

  compute_crit_pts();

}

GRCubic::GRCubic(const GRCubic& cubic) :
    GRCurveGeom(), x_coeff(NULL), y_coeff(NULL), crit_pts()
{

  assert(cubic.isValid());

  n_pts = cubic.n_pts;
  points = NULL;
  n_dim = cubic.n_dim;

  x_coeff = new double[4];
  y_coeff = new double[4];

  for (int i = 0; i < 4; i++) {
    x_coeff[i] = cubic.x_coeff[i];
    y_coeff[i] = cubic.y_coeff[i];
  }

  compute_crit_pts();

}

GRCubic&
GRCubic::operator=(const GRCubic& cubic)
{
  if (&cubic != this) {
    assert(cubic.isValid());

    n_pts = cubic.n_pts;
    points = NULL;
    n_dim = cubic.n_dim;

    for (int i = 0; i < 4; i++) {
      x_coeff[i] = cubic.x_coeff[i];
      y_coeff[i] = cubic.y_coeff[i];
    }

    compute_crit_pts();
  }
  return *this;

}

GRCubic::~GRCubic()
{

  if (x_coeff)
    delete[] x_coeff;
  if (y_coeff)
    delete[] y_coeff;
  crit_pts.clear();

}

void
GRCubic::set_coeffs(const double xcoeff[4], const double ycoeff[4])
{

  int i;
  for (i = 0; i < 4; i++) {
    x_coeff[i] = xcoeff[i];
    y_coeff[i] = ycoeff[i];
  }

}
void
GRCubic::set_xcoeffs(const double xcoeff[4])
{

  int i;
  for (i = 0; i < 4; i++)
    x_coeff[i] = xcoeff[i];

  compute_crit_pts();
}

void
GRCubic::set_ycoeffs(const double ycoeff[4])
{

  int i;
  for (i = 0; i < 4; i++)
    y_coeff[i] = ycoeff[i];

  compute_crit_pts();
}

CubitBox
GRCubic::bounding_box() const
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
GRCubic::param_at_coord(const CubitVector& coord) const
{

  assert(isValid());

  double param;
  bool on_curve = coord_on_curve(coord, &param);

  assert(on_curve);
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  return param;

}

double
GRCubic::param_at_arc_length(const double param_start,
			     const double length) const
{

  //This function is not necessary.

  assert(isValid());
  assert(
      iFuzzyComp(param_start, min_param()) >= 0
	  && iFuzzyComp(param_start, max_param()) <= 0);
  assert(iFuzzyComp(length, 0.) >= 0);
  assert(iFuzzyComp(length, arc_length(param_start, max_param())) <= 0);

  vFatalError("Function not defined for this class",
	      "GRCubic::param_at_arc_length()");

  return -LARGE_DBL;

}

void
GRCubic::coord_at_param(const double param, CubitVector& coord) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  double xx = (x_coeff[0] + x_coeff[1] * param + x_coeff[2] * param * param
      + x_coeff[3] * param * param * param);
  double yy = (y_coeff[0] + y_coeff[1] * param + y_coeff[2] * param * param
      + y_coeff[3] * param * param * param);

  coord.set(xx, yy, 0.);

}

double
GRCubic::coord_at_dist(const double dist, const bool start,
		       CubitVector& coord) const
{

  //Function not necessary for this class
  assert(isValid());
  assert(iFuzzyComp(dist, 0.) >= 0);
  if (start) {
  }

  coord.set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);

  vFatalError("Function not defined for this class",
	      "GRCubic::coord_at_dist()");

  return -LARGE_DBL;

}

double
GRCubic::coord_at_mid_dist(const double param1, const double param2,
			   CubitVector& middle) const
{

  //Function not necessary for this class
  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  middle.set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);

  vFatalError("Function not defined for this class",
	      "GRCubic::coord_at_dist()");

  return -LARGE_DBL;

}

//double GRCubic::coord_at_mid_dist(const CubitVector& coord1,
//                                  const CubitVector& coord2,
//                                  CubitVector& middle) const
//{
//
//  //Function not necessary for this class
//  assert(qValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  middle.set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
//
//  vFatalError("Function not defined for this class",
//      "GRCubic::coord_at_dist()");
//
//  return -LARGE_DBL;
//
//}

void
GRCubic::closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
				const double param1, const double param2) const
{

  assert(isValid());

  double param = closest_param_on_curve(coord, param1, param2);

  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  coord_at_param(param, closest);
  assert(coord_on_curve(closest));

}

double
GRCubic::closest_param_on_curve(const CubitVector& coord, const double param1,
				const double param2) const
{

  assert(isValid());
  double p1, p2;
  if (iFuzzyComp(param1, min_param()) < 0
      && iFuzzyComp(param2, min_param()) < 0) {
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

  //Between critical points, curvature change is monotone. Furthermore, the first and
  //second derivatives of the distance are known. Hence a Newton's method appears
  //to be a good choice to find closest parameter. In the case of a negative second
  //derivative, the Newton's method will not converge, then switch to bisection.

  std::set<double, param_compare> param_set(crit_pts);
  std::set<double>::iterator it1, it2;

  //Erasing params outside of the search range.

  it1 = param_set.begin();
  it2 = param_set.upper_bound(p1);
  if (it2 != param_set.begin())
    param_set.erase(it1, it2);
  param_set.insert(p1);

  it1 = param_set.lower_bound(p2);
  it2 = param_set.end();
  if (it1 != param_set.end())
    param_set.erase(it1, it2);
  param_set.insert(p2);

  //param_set now only contains parameters in the range we are
  //interested in. Here, we initialize min_dist and return_param
  //to the begining of the search range.

  CubitVector normal, lower_coord, upper_coord, lower_vec, upper_vec;
  double lower_param, upper_param, lower_cross, upper_cross;

  lower_param = *(param_set.begin());
  unit_normal(lower_param, normal);
  coord_at_param(lower_param, lower_coord);
  lower_vec = coord - lower_coord;

  double min_dist = lower_vec.length();
  double return_param = lower_param;

  if (iFuzzyComp(min_dist, 0.) == 0)
    return return_param;
  else {
    lower_vec /= min_dist;
    lower_cross = normal.x() * lower_vec.y() - lower_vec.x() * normal.y();
  }

  //Testing every search interval separately to find the minimal distance
  //to the curve.
  for (it1 = ++(param_set.begin()); it1 != param_set.end(); it1++) {

    upper_param = *it1;
    unit_normal(upper_param, normal);
    coord_at_param(upper_param, upper_coord);
    upper_vec = coord - upper_coord;
    double this_dist = upper_vec.length();

    if (iFuzzyComp(this_dist, 0.) == 0)
      return upper_param;
    else {
      upper_vec /= this_dist;
      upper_cross = normal.x() * upper_vec.y() - upper_vec.x() * normal.y();
    }

    //If the closest param is not in the current range, the two
    //cross products will have the same sign. If there is no
    //minimum in the interval, there is not need to invoke the optimizer.
    if (iFuzzyComp((upper_cross * lower_cross), 0.) > 0) {

      CubitVector v = upper_coord - lower_coord;
      v.normalize();
      double angle = GR_acos(lower_vec % v);

      if (angle < 0.5 * M_PI) {
	v *= -1.;
	angle = GR_acos(upper_vec % v);
	if (angle < 0.5 * M_PI) {
	  std::pair<std::set<double>::iterator, bool> P;
	  P = param_set.insert(0.5 * (lower_param + upper_param));
	  if (P.second)
	    it1 = --(P.first);
	  continue;
	}
      }

      if (this_dist < min_dist) {
	min_dist = this_dist;
	return_param = upper_param;
      }

    }

    else { //We might have a minimum, call the optimizer.

      double param;
      bool converged = newton(coord, lower_param, upper_param, param);
      if (!converged)
	converged = bisection(coord, lower_param, upper_param, param);
      assert(converged);

      CubitVector closest;
      coord_at_param(param, closest);
      double dist = (coord - closest).length();

      if (dist < min_dist) {
	min_dist = dist;
	return_param = param;
      }

    }

    lower_cross = upper_cross;
    lower_param = upper_param;
    lower_coord = upper_coord;
    lower_vec = upper_vec;

  }

  return return_param;

}

static double
calcDerivDist_wrt_Param(const GRCubic* const cubic, const CubitVector& target,
			const double param)
{
  CubitVector grad, current;
  cubic->coord_at_param(param, current);
  cubic->first_deriv(param, grad);
  double deriv = 2
      * ((current.x() - target.x()) * grad.x()
	  + (current.y() - target.y()) * grad.y());
  return deriv;
}

bool
GRCubic::coord_on_curve(const CubitVector& coord, double* const param_X,
			double* const param_Y, double* const param_Z) const
{
  // Initialize params so that if we fail, they're already set for that.
  if (param_X != NULL)
    *param_X = - LARGE_DBL;
  if (param_Y != NULL)
    *param_Y = - LARGE_DBL;
  if (param_Z != NULL)
    *param_Z = - LARGE_DBL;

  // The old code (based on root finding for the x coordinate) fails when the
  // spline is vertical, because the roots for x basically dependent only on
  // the small noise in the non-constant coefficients.

  // Instead, try doing both solves.  If they agree, we're done.  Also, if
  // there are no roots for either that are in range, we're also done.
  int num_roots_x, num_roots_y;
  double roots_x[3], roots_y[3];

  solveCubic(x_coeff[3], x_coeff[2], x_coeff[1], x_coeff[0] - coord.x(),
	     &num_roots_x, roots_x);
  solveCubic(y_coeff[3], y_coeff[2], y_coeff[1], y_coeff[0] - coord.y(),
	     &num_roots_y, roots_y);

  bool anyValidX = false, anyValidY = false;

  for (int iy = 0; iy < num_roots_y; iy++) {
    if (iFuzzyComp(roots_y[iy], min_param()) >= 0
	&& iFuzzyComp(roots_y[iy], max_param()) <= 0) {
      anyValidY = true;
    }
  }
  for (int ix = 0; ix < num_roots_x; ix++) {
    if (iFuzzyComp(roots_x[ix], min_param()) >= 0
	&& iFuzzyComp(roots_x[ix], max_param()) <= 0) {
      anyValidX = true;
      for (int iy = 0; iy < num_roots_y; iy++) {
	if (iFuzzyComp(roots_x[ix], roots_y[iy]) == 0) {
	  // I don't like the distance check, because it fails if you scale
	  // the geometry.
	  //    if (iFuzzyComp(roots[i], min_param()) >= 0
	  //        && iFuzzyComp(roots[i], max_param()) <= 0) {
	  //      double yy = (y_coeff[0] + y_coeff[1] * roots[i]
	  //          + y_coeff[2] * pow(roots[i], 2) + y_coeff[3] * pow(roots[i], 3));
	  //      if (fabs(yy - coord.y()) < 1e-6) {
	  if (param_X != NULL)
	    *param_X = roots_x[ix];
	  if (param_Y != NULL)
	    *param_Y = roots_x[ix];
	  if (param_Z != NULL)
	    *param_Z = -LARGE_DBL;
	  return true;
	} // if roots match
      } // loop over y roots
    } // if x root is in range
  } // loop over x roots

    // If neither parametric cubic has a valid root, then forget it.
  if (!anyValidX and !anyValidY)
    return false;

  // If the two solves give inconsistent results, try it the hard way:
  // minimize the distance.

  // The square of the distance error.
  static const double funcTol = 1.e-16;

  // Since the objective function (dist^2) is quadratic, this is surprisingly
  // easy to meet.
  static const double derivTol = 1.e-16;

  // The new (Feb 2015) implementation of this function works by minimizing the
  // distance from the coordinate to this parametric curve, with limits to keep
  // the resulting parameter within the parameter limits.  In practice, this is
  // done by finding a zero of the derivative, then confirming that the distance
  // between that point and the target is very small.
  double paramA = min_param();
  double paramB = max_param();

  // Use Brent's method for root finding.  This implementation follows
  // pseudocode in the Brent's method wikipedia page.

  double derivA = calcDerivDist_wrt_Param(this, coord, paramA);
  double derivB = calcDerivDist_wrt_Param(this, coord, paramB);

  if (iFuzzyComp(derivA * derivB, 0) != -1) {
    if (iFuzzyComp(derivA, 0) == 0 && iFuzzyComp(derivB, 0) == 0) {
      // Bad luck: surface normals for the end points happen to meet at the
      // point in question.

      // Check the slope just off one endpoint:
      double param1 = paramA * 0.9 + paramB * 0.1;
      double deriv1 = calcDerivDist_wrt_Param(this, coord, param1);
      if (deriv1 > 0) {
	// Headed back to A
	paramB = param1;
	derivB = deriv1;
      }
      else {
	// Check the other end
	double param2 = paramA * 0.1 + paramB * 0.9;
	double deriv2 = calcDerivDist_wrt_Param(this, coord, param2);
	if (deriv2 < 0) {
	  // Headed back to B
	  paramA = param2;
	  derivA = deriv2;
	}
	else {
	  // Try the middle and hope for the best.
	  paramA = param1;
	  derivA = deriv1;
	  paramB = param2;
	  derivB = deriv2;
	}
      }
    }

    // See if one of these points matches the input point.  If so, rig things
    // so that we return easily.
    CubitVector testLoc;
    coord_at_param(paramB, testLoc);
    double funcB = (testLoc.x() - coord.x()) * (testLoc.x() - coord.x())
	+ (testLoc.y() - coord.y()) * (testLoc.y() - coord.y());

    coord_at_param(paramA, testLoc);
    double funcA = (testLoc.x() - coord.x()) * (testLoc.x() - coord.x())
	+ (testLoc.y() - coord.y()) * (testLoc.y() - coord.y());
    if (funcA < funcTol) {
      funcB = funcA;
      std::swap(paramA, paramB);
      derivB = 0;
      derivA = 1;
    }
    else if (funcB < funcTol) {
      // Do nothing; just don't do the else clause.
    }
    else {
      // Need to see whether there's some interior point that has opposite sign.
      bool foundOne = false;
      for (int ii = 1; ii < 7; ii++) {
	double param = (paramA * ii + paramB * (7 - ii)) / 7;
	double deriv = calcDerivDist_wrt_Param(this, coord, param);
	if (iFuzzyComp(derivA * deriv, 0) != 1) {
	  // These values are at least signed and zero; more likely opposite signs
	  paramB = param;
	  derivB = deriv;
	  foundOne = true;
	  break;
	}
      }
      if (!foundOne)
	return false; // No zero of the deriv, so no minimum
    }
  }

  if (fabs(derivA) < fabs(derivB)) {
    // Switch so that B has the smaller derivative (better soln guess)
    std::swap(paramA, paramB);
    std::swap(derivA, derivB);
  }
  double paramC = paramA, derivC = derivA;
  double paramS = -10, paramD = -10, delta = 1.e-10;
  bool mflag = true;
  int step = 0;
  while (fabs(derivB) >= derivTol && fabs(paramB - paramA) > delta && step < 20) {
    step++;
    if (iFuzzyComp(derivA, derivC) != 0 && iFuzzyComp(derivB, derivC) != 0) {
      // Inverse quadratic interpolation
//      logMessage(1, "IQI\n");
      paramS = paramA * derivB * derivC
	  / ((derivA - derivB) * (derivA - derivC))
	  + paramB * derivA * derivC / ((derivB - derivA) * (derivB - derivC))
	  + paramC * derivA * derivB / ((derivC - derivA) * (derivC - derivB));
    }
    else {
      // Linear interpolation
//      logMessage(1, "LI\n");
      paramS = paramB - derivB * (paramB - paramA) / (derivB - derivA);
    }
    if ((iFuzzyComp(
	(paramS - (0.75 * paramA + 0.25 * paramB)) * (paramB - paramS), 0) == -1)
	or (mflag and fabs(paramS - paramB) >= 0.5 * fabs(paramB - paramC))
	or (mflag and fabs(paramB - paramC) < delta)
	or (!mflag and fabs(paramS - paramB) >= 0.5 * fabs(paramC - paramD))
	or (!mflag and fabs(paramC - paramD) < delta)) {
      // Bisection
//      logMessage(1, "Bisect\n");
      paramS = 0.5 * (paramA + paramB);
      mflag = true;
    }
    else {
      mflag = false;
    }
    double derivS = calcDerivDist_wrt_Param(this, coord, paramS);
//    logMessage(1, "Step %d. A = (%10.6F %14.8G)\n"
//                  "         B = (%10.6F %14.8G)\n"
//                  "         S = (%10.6F %14.8G)\n",
//                  step, paramA, derivA, paramB, derivB, paramS, derivS);
    paramD = paramC;
    paramC = paramB;
    derivC = derivB;
    if (derivA * derivS < 0) {
      paramB = paramS;
      derivB = derivS;
    }
    else {
      paramA = paramS;
      derivA = derivS;
    }
    if (fabs(derivA) < fabs(derivB)) {
      // Switch so that B has the smaller derivative (better soln guess)
      std::swap(paramA, paramB);
      std::swap(derivA, derivB);
    }
  }

  CubitVector thisLoc;
  coord_at_param(paramB, thisLoc);
  double func = (thisLoc.x() - coord.x()) * (thisLoc.x() - coord.x())
      + (thisLoc.y() - coord.y()) * (thisLoc.y() - coord.y());
  if (func > funcTol) {
//    logMessage(1, "Failing. Param: %10.6G  Dist^2: %10.6G (%10.5f, %10.5f)\n",
//               paramB, func, thisLoc.x(), thisLoc.y());
    return false;
  }
  if (param_X != NULL)
    *param_X = paramB;
  return true;

}

void
GRCubic::first_deriv(const double param, CubitVector& FD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  FD.set(x_coeff[1] + 2. * param * x_coeff[2] + 3. * param * param * x_coeff[3],
	 y_coeff[1] + 2. * param * y_coeff[2] + 3. * param * param * y_coeff[3],
	 0.);

}

//void GRCubic::first_deriv(const CubitVector& coord, CubitVector& FD) const
//{
//
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  first_deriv(param, FD);
//
//}

void
GRCubic::second_deriv(const double param, CubitVector& SD) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  SD.set(2. * x_coeff[2] + 6. * param * x_coeff[3],
	 2. * y_coeff[2] + 6. * param * y_coeff[3], 0.);

}

//void GRCubic::second_deriv(const CubitVector& coord, CubitVector& SD) const
//{
//
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  second_deriv(param, SD);
//
//}

double
GRCubic::curvature(const double param) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  CubitVector V(
      x_coeff[1] + 2. * x_coeff[2] * param + 3. * x_coeff[3] * param * param,
      y_coeff[1] + 2. * y_coeff[2] * param + 3. * y_coeff[3] * param * param,
      0.);
  CubitVector A(2. * x_coeff[2] + 6. * x_coeff[3] * param,
		2. * y_coeff[2] + 6. * y_coeff[3] * param, 0.);

  double cross = V.x() * A.y() - A.x() * V.y();
  double mag = V.length();

  return (fabs(cross) / pow(mag, 3.));

}

//double GRCubic::curvature(const CubitVector& coord) const
//{
//
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  return curvature(param);
//
//}

void
GRCubic::unit_normal(const double param, CubitVector& normal) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  CubitVector tangent;
  unit_tangent(param, tangent);

  normal.set(-tangent.y(), tangent.x(), 0.);

}

//void GRCubic::unit_normal(const CubitVector& coord, CubitVector& normal) const
//{
//
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  unit_normal(param, normal);
//
//}

void
GRCubic::unit_tangent(const double param, CubitVector& tangent) const
{

  assert(isValid());
  assert(
      iFuzzyComp(param, min_param()) >= 0
	  && iFuzzyComp(param, max_param()) <= 0);

  first_deriv(param, tangent);
  tangent.normalize();

}

//void GRCubic::unit_tangent(const CubitVector& coord, CubitVector& tangent) const
//{
//
//  assert(coord_on_curve(coord));
//
//  double param = param_at_coord(coord);
//  unit_tangent(param, tangent);
//
//}

double
GRCubic::arc_length() const
{

  assert(isValid());
  return arc_length(min_param(), max_param());

}

double
GRCubic::arc_length(const double param1, const double param2) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  // Integrates (sqrt((dx/dt)^2 + (dy/dt)^2)) dt from param1 to param2
  // using Simpson's rule with three levels of refinement, plus
  // Richardson extrapolaton to estimate the result for an infinitely
  // large number of integration intervals.

  // Based on numerical experimentation, a precision equals to square
  // root of machine precision is obtained using interval values of
  // about param2 - param1 = 0.1 (although this is not necessarly
  // always true since parameterization is not linear).

  double interval = param2 - param1;
  const int maxPoints = 1024;
  double speed[maxPoints + 1]; // That had better be enough points!

  //   while(interval > 0.1) {
  //     interval *= static_cast<double>(num_int);
  //     ++num_int;
  //     interval /= static_cast<double>(num_int);
  //   }

  // Start with values @ 0 (param1) and 256 (param2).
  // First loop: j = 1: compute @ 128; [0] + 4*[128] + [256]; *interval/6
  // Second loop: j = 1: compute @ 64; [0] + 4*[64] + [128];
  //              j = 3: compute @ 192; [128] + 4*[192] + [256]; *interval/12
  CubitVector FD;
  first_deriv(param1, FD);
  speed[0] = FD.length();

  first_deriv(param2, FD);
  speed[maxPoints] = FD.length();

  double length = 0, oldLength = -1;
  int jMax = 1;
  do {
    jMax *= 2;
    int stride = maxPoints / jMax;
    oldLength = length;
    length = 0;
    for (int j = 1; j < jMax; j += 2) {
      first_deriv((param1 + (interval * j) / jMax), FD);
      speed[j * stride] = FD.length();
      length += speed[(j - 1) * stride] + 4 * speed[(j) * stride]
	  + speed[(j + 1) * stride];
    }
    length *= (interval / (3 * jMax));
    logMessage(4, "jMax: %3d  Length: %18.16F\n", jMax, length);
  }
  while (jMax < maxPoints && iFuzzyComp(oldLength, length) != 0);
  return length;

}

//double GRCubic::arc_length(const CubitVector& coord1, const CubitVector& coord2) const
//{
//
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  double param1 = param_at_coord(coord1);
//  double param2 = param_at_coord(coord2);
//  return arc_length(param1, param2);
//
//}

//void GRCubic::center_point(CubitVector& mid_pt) const
//{
//
//  assert(isValid());
//  center_point(min_param(), max_param(), mid_pt);
//
//}

void
GRCubic::center_point(const double param1, const double param2,
		      CubitVector& mid_pt) const
{

  assert(isValid());
  assert(iFuzzyComp(param1, param2) <= 0);
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);

  double total_length = arc_length(param1, param2);
  double target_length = 0.5 * total_length;

  double param_lo = param1;
  double param_up = param2;
  double param_new = 0.5 * (param_lo + param_up);

  total_length = arc_length(param1, param_new);

  while (iFuzzyComp(total_length, target_length) != 0) {

    if (total_length <= target_length) {
      param_lo = param_new;
    }
    else {
      param_up = param_new;
    }

    double split_ratio = 0.5;

    //This could speed up convergence by taking advantage of
    //first order info... or not since there are two more
    //arc_length calls involved.

    //     double length1 = arc_length(param1, param_lo);
    //     double length2 = length1 + arc_length(param_lo, param_up);
    //     double split_ratio = (length2 - target_length) / (length2 - length1);

    param_new = param_lo + split_ratio * (param_up - param_lo);
    total_length = arc_length(param1, param_new);

  }

  coord_at_param(param_new, mid_pt);

}

//void GRCubic::center_point(const CubitVector& coord1,
//                           const CubitVector& coord2, CubitVector& mid_pt) const
//{
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  center_point(param_at_coord(coord1), param_at_coord(coord2), mid_pt);
//
//}

//double GRCubic::TVT() const
//{
//
//  assert(isValid());
//  return TVT(min_param(), max_param());
//
//}

double
GRCubic::TVT(const double param1, const double param2) const
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
  CubitVector tangent_beg, tangent_end;

  std::set<double>::const_iterator itL = crit_pts.lower_bound(param1);
  std::set<double>::const_iterator itU = crit_pts.lower_bound(param2);

  if (itL == itU) {
    unit_tangent(param1, tangent_beg);
    unit_tangent(param2, tangent_end);
    tot_var += fabs(GR_acos(tangent_beg % tangent_end));
  }
  else {
    --itU;
    unit_tangent(param1, tangent_beg);
    unit_tangent(*itL, tangent_end);
    tot_var += fabs(GR_acos(tangent_beg % tangent_end));

    while (itL != itU) {
      unit_tangent(*itL, tangent_beg);
      unit_tangent(*(++itL), tangent_end);
      tot_var += fabs(GR_acos(tangent_beg % tangent_end));
    }

    unit_tangent(*itU, tangent_beg);
    unit_tangent(param2, tangent_beg);
    tot_var += fabs(GR_acos(tangent_beg % tangent_end));
  }

  return tot_var;

}

//double GRCubic::TVT(const CubitVector& coord1, const CubitVector& coord2) const
//{
//
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return TVT(param_at_coord(coord1), param_at_coord(coord2));
//
//}

double
GRCubic::mid_TVT(const double param1, const double param2) const
{

  double total_TVT = TVT(param1, param2);
  double target_TVT = 0.5 * total_TVT;

  double param_lo = param1;
  double param_up = param2;
  double param_new = 0.5 * (param_lo + param_up);

  total_TVT = TVT(param1, param_new);

  if (fabs(total_TVT) < 1.e-3)
    return 0.5 * (param1 + param2);

  while (fabs(total_TVT - target_TVT) > cg_tolerance) {

    if (total_TVT <= target_TVT) {
      param_lo = param_new;
    }
    else {
      param_up = param_new;
    }

    double split_ratio = 0.5;

    //     double TVT1 = TVT(param1, param_lo);
    //     double TVT2 = TVT1 + TVT(param_lo, param_up);
    //     double split_ratio = (TVT2 - target_TVT) / (TVT2 - TVT1);

    param_new = param_lo + split_ratio * (param_up - param_lo);
    total_TVT = TVT(param1, param_new);

  }

  return param_new;

}

//double GRCubic::mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const
//{
//
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//
//  return mid_TVT(param_at_coord(coord1), param_at_coord(coord2));
//
//}

void
GRCubic::interior_extrema(DLIList<CubitVector*>& point_list) const
{

  assert(isValid());

  std::set<double>::const_iterator it;
  for (it = crit_pts.begin(); it != crit_pts.end(); it++) {
    double param = *it;

    double x_ext = (3. * x_coeff[3] * param * param + 2. * x_coeff[2] * param
	+ x_coeff[1]);
    double y_ext = (3. * y_coeff[3] * param * param + 2. * y_coeff[2] * param
	+ y_coeff[1]);

    if (iFuzzyComp(x_ext, 0.) == 0 || iFuzzyComp(y_ext, 0.) == 0) {
      CubitVector coord;
      coord_at_param(param, coord);
      point_list.insert(new CubitVector(coord));
    }
  }
}

bool
GRCubic::line_intersect(const double coords1[3], const double coords2[3],
			const double begParam, const double endParam,
			std::vector<double> & params) const
{

  double v0x = coords1[0];
  double v0y = coords1[1];
  double v1x = coords2[0];
  double v1y = coords2[1];
  int iNRoots;
  double adRoots[3];

  CubitBox box = bounding_box();
  if (box.max_x() < min(v0x, v1x) || box.min_x() > max(v0x, v1x)
      || box.max_y() < min(v0y, v1y) || box.min_y() > max(v0y, v1y)) {
    return false;
  }
  if (fabs(v1x - v0x) > 1.e-12) {
    double m = (v1y - v0y) / (v1x - v0x);
    // going to solve y(t) - m*(x(t)-x0)-y0 = 0

    solveCubic(y_coeff[3] - m * x_coeff[3], y_coeff[2] - m * x_coeff[2],
	       y_coeff[1] - m * x_coeff[1],
	       y_coeff[0] - m * x_coeff[0] + m * v0x - v0y, &iNRoots, adRoots);
  }
  else {
    // if the face is vertical, need to be more clever

    solveCubic(x_coeff[3], x_coeff[2], x_coeff[1],
	       x_coeff[0] - 0.5 * (v1x + v0x), &iNRoots, adRoots);
  }
  bool qRetVal = false;
  for (int i = 0; i < iNRoots; i++) {
    if (iFuzzyComp(adRoots[i], 1.) < 1 && iFuzzyComp(adRoots[i], 0.) > -1) {
      CubitVector coord;
      coord_at_param(adRoots[i], coord);
      if (iFuzzyComp(coord.y(), max(v0y, v1y)) < 1
	  && iFuzzyComp(coord.y(), min(v0y, v1y)) > -1
	  && iFuzzyComp(coord.x(), max(v0x, v1x)) < 1
	  && iFuzzyComp(coord.x(), min(v0x, v1x)) > -1
	  && iFuzzyComp(adRoots[i], endParam) < 1
	  && iFuzzyComp(adRoots[i], begParam) > -1) {
	params.push_back(adRoots[i]);
	qRetVal = true;
      }
    }
  }
  return qRetVal;
}

double
GRCubic::x(const int) const
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::X(const int)");
  return -DBL_MAX;
}

double
GRCubic::y(const int) const
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::Y(const int)");
  return -DBL_MAX;
}

double
GRCubic::z(const int) const
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::Z(const int)");
  return -DBL_MAX;
}

void
GRCubic::get_coord(const int, CubitVector&) const
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::get_coord(const int, CubitVector&)");
}

const CubitVector&
GRCubic::get_coord(const int) const
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::get_coord(const int)");
  // Execution can't get here; the following lines silence compiler warnings.
  static CubitVector silencer(0, 0, 0);
  return silencer;
}

void
GRCubic::set_x(const int, const double)
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::set_X(const int, const double)");
}

void
GRCubic::set_y(const int, const double)
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::set_Y(const int, const double)");
}

void
GRCubic::set_z(const int, const double)
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::set_Z(const int, const double)");
}

void
GRCubic::set_coord(const int, const CubitVector&)
{
  vFatalError("A cubic curve is not defined with points.",
	      "GRCubic::set_coord(const int, const CubitVector&)");
}

bool
GRCubic::isValid() const
{

  if (n_pts != 0) {
    return false;
  }
  if (points != NULL) {
    return false;
  }
  if (crit_pts.size() < 2) {
    return false;
  }
  if (iFuzzyComp(*(crit_pts.begin()), min_param()) != 0) {
    return false;
  }
  if (iFuzzyComp(*(--(crit_pts.end())), max_param()) != 0) {
    return false;
  }

  if (x_coeff == NULL)
    return false;
  if (y_coeff == NULL)
    return false;
  if (x_coeff[0] < -1.e200)
    return false;
  if (x_coeff[1] < -1.e200)
    return false;
  if (x_coeff[2] < -1.e200)
    return false;
  if (x_coeff[3] < -1.e200)
    return false;
  if (y_coeff[0] < -1.e200)
    return false;
  if (y_coeff[1] < -1.e200)
    return false;
  if (y_coeff[2] < -1.e200)
    return false;
  if (y_coeff[3] < -1.e200)
    return false;

  return true;

}

void
GRCubic::compute_crit_pts()
{

  crit_pts.clear();

  //Firstly, begin and end points of the cubic curve.

  crit_pts.insert(min_param());
  crit_pts.insert(max_param());

  //Secondly, roots of x''(t), y''(t) and z''(t).

  double one_third = 1. / 3.;
  if (iFuzzyComp(x_coeff[3], 0.) != 0) {
    double root = -one_third * (x_coeff[2] / x_coeff[3]);
    if (root > (min_param() + 5.e-3) && root < (max_param() - 5.e-3))
      crit_pts.insert(root);
  }
  if (iFuzzyComp(y_coeff[3], 0.) != 0) {
    double root = -one_third * (y_coeff[2] / y_coeff[3]);
    if (root > (min_param() + 5.e-3) && root < (max_param() - 5.e-3))
      crit_pts.insert(root);
  }

  //Thirdly, roots of x'(t) and y'(t).

  double A, B, C, root1, root2;
  int num_roots;

  A = 3. * x_coeff[3];
  B = 2. * x_coeff[2];
  C = x_coeff[1];
  solveQuadratic(A, B, C, &num_roots, &root1, &root2);
  insert_crit_pts(num_roots, root1, root2);

  A = 3 * y_coeff[3];
  B = 2 * y_coeff[2];
  C = y_coeff[1];
  solveQuadratic(A, B, C, &num_roots, &root1, &root2);
  insert_crit_pts(num_roots, root1, root2);

  //Finally, parameters where curvature = 0.

  A = x_coeff[2] * y_coeff[3] - y_coeff[2] * x_coeff[3];
  B = x_coeff[1] * y_coeff[3] - y_coeff[1] * x_coeff[3];
  C = one_third * (x_coeff[1] * y_coeff[2] - y_coeff[1] * x_coeff[2]);
  solveQuadratic(A, B, C, &num_roots, &root1, &root2);
  insert_crit_pts(num_roots, root1, root2);

}

void
GRCubic::insert_crit_pts(const int num_roots, const double root1,
			 const double root2)
{

  switch (num_roots)
    {
    case 1:
      {
	if (root1 > (min_param() + 5.e-3) && root1 < (max_param() - 5.e-3))
	  crit_pts.insert(root1);
	break;
      }
    case 2:
      {
	if (root1 > (min_param() + 5.e-3) && root1 < (max_param() - 5.e-3))
	  crit_pts.insert(root1);
	if (root2 > (min_param() + 5.e-3) && root2 < (max_param() - 5.e-3))
	  crit_pts.insert(root2);
	break;
      }
    default:
      break;
    }

}

bool
GRCubic::newton(const CubitVector& coord, const double lower_param,
		const double upper_param, double& return_param) const
{

  assert(iFuzzyComp(lower_param, upper_param) < 0);

  return_param = 0.5 * (lower_param + upper_param);

  int num_iters = 0;
  CubitVector loc, FD, SD;
  double dist, grad, scale_grad, hess, step, direction;

  coord_at_param(return_param, loc);
  first_deriv(return_param, FD);
  second_deriv(return_param, SD);

  dist = (loc - coord).length_squared();
  grad = 2. * ((loc.x() - coord.x()) * FD.x() + (loc.y() - coord.y()) * FD.y());
  hess = 2.
      * (pow(FD.x(), 2.) + pow(FD.y(), 2.) + SD.x() * (loc.x() - coord.x())
	  + SD.y() * (loc.y() - coord.y()));

  { //Compute a scale gradient to define stopping criterion.
    coord_at_param(lower_param, loc);
    first_deriv(lower_param, FD);
    double grad_up = 2.
	* ((loc.x() - coord.x()) * FD.x() + (loc.y() - coord.y()) * FD.y());
    coord_at_param(upper_param, loc);
    first_deriv(upper_param, FD);
    double grad_lo = 2.
	* ((loc.x() - coord.x()) * FD.x() + (loc.y() - coord.y()) * FD.y());

    scale_grad = MAX3(fabs(grad_lo), fabs(grad_up), fabs(grad));
  }

  do {

    if ((fabs(grad) / scale_grad) < cg_tolerance)
      return true;

    if (hess <= 0.) { //hessian is negative if extremum is a max.
      CubitVector coord1, coord2;
      coord_at_param(lower_param, coord1);
      coord_at_param(upper_param, coord2);
      double dist1 = (coord1 - coord).length_squared();
      double dist2 = (coord2 - coord).length_squared();
      return_param = dist1 < dist2 ? lower_param : upper_param;
      return true;
    }

    direction = -grad / hess; //Taking a Newton step.
    step = 1.;

    do {

      if (step < cg_tolerance)
	return false;

      double new_param = return_param + (step * direction);

      if (new_param > upper_param) { //return_param cannot be greater than upper_param
	if (iFuzzyComp(return_param, upper_param) == 0)
	  return true;
	step = (return_param - upper_param) * hess / grad;
	assert(step > 0.);
	new_param = upper_param;
      }
      else if (new_param < lower_param) { //return_param cannot be smaller than lower_param
	if (iFuzzyComp(return_param, lower_param) == 0)
	  return true;
	step = (return_param - lower_param) * hess / grad;
	assert(step > 0.);
	new_param = lower_param;
      }

      assert(iFuzzyComp(lower_param, upper_param) < 0);
      assert(
	  iFuzzyComp(new_param, upper_param) <= 0
	      && iFuzzyComp(new_param, lower_param) >= 0);

      coord_at_param(new_param, loc);
      double new_dist = (loc - coord).length_squared();

      if (iFuzzyComp(new_dist, dist + (1.e-4 * step * grad * direction)) <= 0) {

	first_deriv(new_param, FD);
	double new_grad = 2.
	    * ((loc.x() - coord.x()) * FD.x() + (loc.y() - coord.y()) * FD.y());

	if (iFuzzyComp(fabs(new_grad * direction), 0.9 * fabs(grad * direction))
	    <= 0) {
	  second_deriv(return_param, SD);
	  return_param = new_param;
	  grad = new_grad;
	  dist = new_dist;
	  hess = 2.
	      * (pow(FD.x(), 2.) + pow(FD.y(), 2.)
		  + SD.x() * (loc.x() - coord.x())
		  + SD.y() * (loc.y() - coord.y()));
	  break;
	}

	else
	  step *= 0.5;

      }

      else
	step *= 0.5;

    }
    while (1.);

  }
  while (++num_iters < 15);

  return false;

}

bool
GRCubic::bisection(const CubitVector& coord, const double lower_param,
		   const double upper_param, double& return_param) const
{

  assert(iFuzzyComp(lower_param, upper_param) < 0);

  bool converged = true;
  double param_up = upper_param, param_lo = lower_param;

  CubitVector coord_lo, coord_up;
  coord_at_param(param_up, coord_up);
  coord_at_param(param_lo, coord_lo);

  CubitVector normal_lo, normal_up;
  unit_normal(param_up, normal_up);
  unit_normal(param_lo, normal_lo);

  CubitVector diff = coord - coord_up;
  double cross_up = normal_up.x() * diff.y() - diff.x() * normal_up.y();
  diff = coord - coord_lo;
  double cross_lo = normal_lo.x() * diff.y() - diff.x() * normal_lo.y();

  assert(iFuzzyComp(cross_up * cross_lo, 0.) <= 0);

  int iters = 0;
  double cross, scale_cross = std::max(fabs(cross_lo), fabs(cross_up));
  CubitVector this_normal, this_coord;

  do {

    //double ratio = fabs(cross_lo) / (fabs(cross_lo) + fabs(cross_up));
    double ratio = 0.5;
    return_param = param_lo + ratio * (param_up - param_lo);

    unit_normal(return_param, this_normal);
    coord_at_param(return_param, this_coord);

    diff = coord - this_coord;
    cross = (this_normal.x() * diff.y()) - (diff.x() * this_normal.y());

    if (cross_up * cross <= 0.) {
      param_lo = return_param;
      cross_lo = cross;
    }
    else {
      assert(cross_lo * cross <= 0.);
      param_up = return_param;
      cross_up = cross;
    }

  }
  while ((fabs(cross) > 1.e-16) && (fabs(cross) / scale_cross) > cg_tolerance
      && ++iters < 50);

  if (iters == 50) {
    converged = false;
  }

  assert(iFuzzyComp(return_param, lower_param) >= 0);
  assert(iFuzzyComp(return_param, upper_param) <= 0);

  return converged;

}
