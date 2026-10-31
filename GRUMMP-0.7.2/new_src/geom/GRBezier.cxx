#include "GR_GRCurveGeom.h"
#include <iostream>

GRBezier::GRBezier() :
    GRSpline()
{

  n_pts = 4;
  n_dim = 2;
  points = new CubitVector[4];
  cg_tolerance = 1e-12;

  int i;
  for (i = 0; i < n_pts; i++)
    points[i].set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);

  num_cubics = 1;
  cubics = new GRCubic*[num_cubics];

}

GRBezier::GRBezier(const CubitVector input_points[4]) :
    GRSpline()
{

  n_pts = 4;
  n_dim = 2;
  points = new CubitVector[4];
  cg_tolerance = 1e-12;

  int i;
  for (i = 0; i < n_pts; i++)
    points[i] = input_points[i];

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0)
    points[1] = points[3];

  if (iFuzzyComp(points[2].x(), points[3].x()) == 0
      && iFuzzyComp(points[2].y(), points[3].y()) == 0)
    points[2] = points[0];

  double x_coeff[4], y_coeff[4];
  compute_coeffs(x_coeff, y_coeff);

  num_cubics = 1;
  cubics = new GRCubic*[num_cubics];
  cubics[0] = new GRCubic(x_coeff, y_coeff);

  assert(isValid());

}

GRBezier::GRBezier(const GRBezier& B) :
    GRSpline()
{

  assert(B.isValid());

  n_pts = B.n_pts;
  n_dim = B.n_dim;
  cg_tolerance = B.cg_tolerance;
  points = new CubitVector[4];
  for (int i = 0; i < 4; i++)
    points[i] = B.points[i];

  num_cubics = B.num_cubics;
  cubics = new GRCubic*[num_cubics];
  cubics[0] = new GRCubic(*(B.cubics[0]));

  assert(isValid());

}

GRBezier&
GRBezier::operator=(const GRBezier& B)
{

  assert(B.isValid());

  n_pts = B.n_pts;
  n_dim = B.n_dim;

  for (int i = 0; i < 4; i++)
    points[i] = B.points[i];

  cubics[0] = new GRCubic(*(B.cubics[0]));

  assert(isValid());
  return *this;

}

GRBezier::~GRBezier()
{
}

double
GRBezier::x(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].x();

}

double
GRBezier::y(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].y();

}

double
GRBezier::z(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt].z();

}

void
GRBezier::get_coord(const int pt, CubitVector& coord) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  coord.set(points[pt].x(), points[pt].y(), 0.);

}

const CubitVector&
GRBezier::get_coord(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  return points[pt];

}

void
GRBezier::set_x(const int pt, const double X_in)
{

  assert(pt >= 0 && pt < n_pts);
  points[pt].x(X_in);
  set_coeffs();

}

void
GRBezier::set_y(const int pt, const double Y_in)
{

  assert(pt >= 0 && pt < n_pts);
  points[pt].y(Y_in);
  set_coeffs();

}

void
GRBezier::set_z(const int, const double)
{
  vFatalError("Modifying z-coord is not allowed on this 2D curve.",
	      "GRBezier::set_Z(const int, const double)");
}

void
GRBezier::set_coord(const int pt, const CubitVector& new_coord)
{

  assert(pt >= 0 && pt < n_pts);
  assert(iFuzzyComp(new_coord.z(), 0.) == 0);

  points[pt] = new_coord;
  set_coeffs();

}

bool
GRBezier::isValid() const
{

  if (points == NULL) {
    return false;
  }
  if (points[0].x() < -1.e200) {
    return false;
  }
  if (points[0].y() < -1.e200) {
    return false;
  }
  if (points[1].x() < -1.e200) {
    return false;
  }
  if (points[1].y() < -1.e200) {
    return false;
  }
  if (points[2].x() < -1.e200) {
    return false;
  }
  if (points[2].y() < -1.e200) {
    return false;
  }
  if (points[3].x() < -1.e200) {
    return false;
  }
  if (points[3].y() < -1.e200) {
    return false;
  }

  if (cubics == NULL) {
    return false;
  }
  if (num_cubics != 1) {
    return false;
  }
  if (!cubics[0]->isValid()) {
    return false;
  }
  if (n_pts != 4) {
    return false;
  }
  if (n_dim != 2) {
    return false;
  }

  return true;

}

void
GRBezier::compute_coeffs(double x_coeff[4], double y_coeff[4])
{

  x_coeff[0] = points[0].x();
  x_coeff[1] = -3 * points[0].x() + 3 * points[1].x();
  x_coeff[2] = 3 * points[0].x() - 6 * points[1].x() + 3 * points[2].x();
  x_coeff[3] = -points[0].x() + 3 * points[1].x() - 3 * points[2].x()
      + points[3].x();

  y_coeff[0] = points[0].y();
  y_coeff[1] = -3 * points[0].y() + 3 * points[1].y();
  y_coeff[2] = 3 * points[0].y() - 6 * points[1].y() + 3 * points[2].y();
  y_coeff[3] = -points[0].y() + 3 * points[1].y() - 3 * points[2].y()
      + points[3].y();

}

void
GRBezier::set_coeffs()
{

  double x_coeff[4], y_coeff[4];
  compute_coeffs(x_coeff, y_coeff);
  cubics[0]->set_coeffs(x_coeff, y_coeff);

}
