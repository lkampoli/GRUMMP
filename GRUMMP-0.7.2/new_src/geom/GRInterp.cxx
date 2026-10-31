#include "GR_GRCurveGeom.h"

GRInterp::GRInterp() :
    GRSpline()
{

  num_cubics = 0;
  cubics = NULL;
  closed = false;

}

GRInterp::GRInterp(const CubitVector input_pts[], const int num_pts,
		   const bool is_closed) :
    GRSpline()
{

  if (num_pts < 4)
    vFatalError(
	"At least four points are needed to build a notAKnot interpolated spline.",
	"GRInterp::GRInterp(const CubitVector, const int, const bool)");

  n_pts = num_pts;
  points = new CubitVector[n_pts];
  n_dim = 2;

  num_cubics = n_pts - 1;
  closed = is_closed;

  int i;
  for (i = 0; i < n_pts; i++)
    points[i] = input_pts[i];

  double (*x_coeff)[4] = new double[num_cubics][4];
  double (*y_coeff)[4] = new double[num_cubics][4];

//  if(!closed) compute_coeffs_natural(x_coeff, y_coeff);
  if (!closed)
    compute_coeffs_notAKnot(x_coeff, y_coeff);
  else
    compute_coeffs_closed(x_coeff, y_coeff);

  cubics = new GRCubic*[num_cubics];

  for (i = 0; i < num_cubics; i++)
    cubics[i] = new GRCubic(x_coeff[i], y_coeff[i]);
//  for(i = 0; i < num_cubics; i++){
//  	for(GR_index_t j = 0; j < 4; j++){
//  		printf("%e ",x_coeff[i][j]);
//  	}
//  	printf("cubic %d x\n",i);
//  	for(GR_index_t j = 0; j < 4; j++){
//  		printf("%e ",y_coeff[i][j]);
//  	}
//  	printf("cubic %d y\n",i);
//  }
//

  delete[] x_coeff;
  delete[] y_coeff;

}

GRInterp::GRInterp(const GRInterp& interp) :
    GRSpline()
{

  assert(interp.isValid());

  n_pts = interp.n_pts;
  n_dim = 2;
  num_cubics = interp.num_cubics;

  points = new CubitVector[n_pts];
  int i;
  for (i = 0; i < n_pts; i++)
    points[i] = interp.points[i];

  cubics = new GRCubic*[num_cubics];
  for (i = 0; i < n_pts - 1; i++)
    cubics[i] = new GRCubic(*(interp.cubics[i]));

  closed = interp.closed;

  assert(isValid());

}

GRInterp&
GRInterp::operator=(const GRInterp& interp)
{

  assert(interp.isValid());

  n_pts = interp.n_pts;
  n_dim = 2;
  num_cubics = interp.num_cubics;

  points = new CubitVector[n_pts];
  int i;
  for (i = 0; i < n_pts; i++)
    points[i] = interp.points[i];

  cubics = new GRCubic*[num_cubics];
  for (i = 0; i < n_pts - 1; i++)
    cubics[i] = new GRCubic(*(interp.cubics[i]));

  closed = interp.closed;

  assert(isValid());

  return *this;

}

GRInterp::~GRInterp()
{
}

double
GRInterp::x(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  return points[pt].x();

}

double
GRInterp::y(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  return points[pt].y();

}

double
GRInterp::z(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  return points[pt].z();

}

void
GRInterp::get_coord(const int pt, CubitVector& coord) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  coord = points[pt];

}

const CubitVector&
GRInterp::get_coord(const int pt) const
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  return points[pt];

}

void
GRInterp::set_x(const int pt, const double X_in)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  if (pt == 0 || pt == n_pts - 1) {
    if (iFuzzyComp(points[0].x(), points[n_pts - 1].x()) == 0
	&& iFuzzyComp(points[0].y(), points[n_pts - 1].y()) == 0
	&& iFuzzyComp(points[0].z(), points[n_pts - 1].z()) == 0) {
      points[0].x(X_in);
      points[n_pts - 1].x(X_in);
    }
    else
      points[pt].x(X_in);
  }

  double (*x_coeff)[4] = new double[num_cubics][4];
  double (*y_coeff)[4] = new double[num_cubics][4];

  if (!closed)
    compute_coeffs_natural(x_coeff, y_coeff);
  else
    compute_coeffs_closed(x_coeff, y_coeff);

  int i;
  for (i = 0; i < num_cubics; i++)
    cubics[i]->set_coeffs(x_coeff[i], y_coeff[i]);

  delete[] x_coeff;
  delete[] y_coeff;

}

void
GRInterp::set_y(const int pt, const double Y_in)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);

  if (pt == 0 || pt == n_pts - 1) {
    if (iFuzzyComp(points[0].x(), points[n_pts - 1].x()) == 0
	&& iFuzzyComp(points[0].y(), points[n_pts - 1].y()) == 0
	&& iFuzzyComp(points[0].z(), points[n_pts - 1].z()) == 0) {
      points[0].y(Y_in);
      points[n_pts - 1].y(Y_in);
    }
    else
      points[pt].y(Y_in);
  }

  double (*x_coeff)[4] = new double[num_cubics][4];
  double (*y_coeff)[4] = new double[num_cubics][4];

  if (!closed)
    compute_coeffs_natural(x_coeff, y_coeff);
  else
    compute_coeffs_closed(x_coeff, y_coeff);

  int i;
  for (i = 0; i < num_cubics; i++)
    cubics[i]->set_coeffs(x_coeff[i], y_coeff[i]);

  delete[] x_coeff;
  delete[] y_coeff;
}

void
GRInterp::set_z(const int, const double)
{
  vFatalError("Cannot modify z-coord for a 2D curve.",
	      "GRInterp::set_z(const int, const double)");
}

void
GRInterp::set_coord(const int pt, const CubitVector& new_coord)
{

  assert(isValid());
  assert(pt >= 0 && pt < n_pts);
  assert(iFuzzyComp(new_coord.z(), 0.) == 0);

  if (pt == 0 || pt == n_pts - 1) {
    if (iFuzzyComp(points[0].x(), points[n_pts - 1].x()) == 0
	&& iFuzzyComp(points[0].y(), points[n_pts - 1].y()) == 0
	&& iFuzzyComp(points[0].z(), points[n_pts - 1].z()) == 0)
      points[0] = points[n_pts - 1] = new_coord;
    else
      points[pt] = new_coord;
  }

  double (*x_coeff)[4] = new double[num_cubics][4];
  double (*y_coeff)[4] = new double[num_cubics][4];

  if (!closed)
    compute_coeffs_natural(x_coeff, y_coeff);
  else
    compute_coeffs_closed(x_coeff, y_coeff);

  int i;
  for (i = 0; i < num_cubics; i++)
    cubics[i]->set_coeffs(x_coeff[i], y_coeff[i]);

  delete[] x_coeff;
  delete[] y_coeff;

}

bool
GRInterp::isValid() const
{

  if (n_pts < 3) {
    return false;
  }
  if (points == NULL) {
    return false;
  }
  if (cubics == NULL) {
    return false;
  }

  return true;

}

void
GRInterp::compute_coeffs_natural(double (*x_coeff)[4], double (*y_coeff)[4])
{

  //Computes the coefficients of the natural spline.
  //That is the spline with zero second derivative at both ends.

  double** LHS;
  double* RHS = new double[n_pts];

  LHS = new double*[n_pts];

  int i;
  for (i = 0; i < n_pts; i++)
    LHS[i] = new double[3];

  //Computing x coefficients
  LHS[0][0] = 0;
  LHS[0][1] = 2;
  LHS[0][2] = 1;
  RHS[0] = 3. * (points[1].x() - points[0].x());

  // This is a re-derivation of Serge's BC.  CFOG 6/2014.
  LHS[n_pts - 2][0] = 2.;
  LHS[n_pts - 2][1] = 7.;
  LHS[n_pts - 2][2] = 0.;
  RHS[n_pts - 2] = 3. * (points[n_pts - 1].x() - points[n_pts - 2].x())
      + 6 * (points[n_pts - 2].x() - points[n_pts - 3].x());

  for (i = 1; i < n_pts - 2; i++) {
    LHS[i][0] = 1.;
    LHS[i][1] = 4.;
    LHS[i][2] = 1.;
    RHS[i] = 3. * (points[i + 1].x() - points[i - 1].x());
  }

//  for (i = 0; i < n_pts - 1; i++) {
//    logMessage(2, "%d LHS: %5f %5f %5f RHS: %5f\n", i,
//               LHS[i][0], LHS[i][1], LHS[i][2], RHS[i]);
//  }

  solveTriDiag(LHS, RHS, n_pts - 1);

  for (i = 0; i < n_pts - 2; i++) {
    x_coeff[i][0] = points[i].x();
    x_coeff[i][1] = RHS[i];
    x_coeff[i][2] = 3. * (points[i + 1].x() - points[i].x()) - 2. * RHS[i]
	- RHS[i + 1];
    x_coeff[i][3] = 2. * (points[i].x() - points[i + 1].x()) + RHS[i]
	+ RHS[i + 1];
//    logMessage(2, "x-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//             i, x_coeff[i][0], x_coeff[i][1], x_coeff[i][2], x_coeff[i][3]);
  }
  x_coeff[i][0] = points[i].x();
  x_coeff[i][1] = RHS[i];
  x_coeff[i][2] = x_coeff[i - 1][2] + 3 * x_coeff[i - 1][3];
  x_coeff[i][3] = points[i + 1].x() - x_coeff[i][0] - x_coeff[i][1]
      - x_coeff[i][2];
//  logMessage(2, "x-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//           i, x_coeff[i][0], x_coeff[i][1], x_coeff[i][2], x_coeff[i][3]);

  //Getting y coefficients
  LHS[0][0] = 0;
  LHS[0][1] = 2;
  LHS[0][2] = 1;
  RHS[0] = 3. * (points[1].y() - points[0].y());

  // This is a re-derivation of Serge's BC.  CFOG 6/2014.
  LHS[n_pts - 2][0] = 2.;
  LHS[n_pts - 2][1] = 7.;
  LHS[n_pts - 2][2] = 0.;
  RHS[n_pts - 2] = 3. * (points[n_pts - 1].y() - points[n_pts - 2].y())
      + 6 * (points[n_pts - 2].y() - points[n_pts - 3].y());

  for (i = 1; i < n_pts - 2; i++) {
    LHS[i][0] = 1.;
    LHS[i][1] = 4.;
    LHS[i][2] = 1.;
    RHS[i] = 3. * (points[i + 1].y() - points[i - 1].y());
  }

  for (i = 0; i < n_pts - 1; i++) {
//    logMessage(2, "%d LHS: %5f %5f %5f RHS: %5f\n", i,
//               LHS[i][0], LHS[i][1], LHS[i][2], RHS[i]);
  }

  solveTriDiag(LHS, RHS, n_pts - 1);

  for (i = 0; i < n_pts - 2; i++) {
    y_coeff[i][0] = points[i].y();
    y_coeff[i][1] = RHS[i];
    y_coeff[i][2] = 3. * (points[i + 1].y() - points[i].y()) - 2. * RHS[i]
	- RHS[i + 1];
    y_coeff[i][3] = 2. * (points[i].y() - points[i + 1].y()) + RHS[i]
	+ RHS[i + 1];
//    logMessage(2, "y-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//             i, y_coeff[i][0], y_coeff[i][1], y_coeff[i][2], y_coeff[i][3]);
  }
  y_coeff[i][0] = points[i].y();
  y_coeff[i][1] = RHS[i];
  y_coeff[i][2] = y_coeff[i - 1][2] + 3 * y_coeff[i - 1][3];
  y_coeff[i][3] = points[i + 1].y() - y_coeff[i][0] - y_coeff[i][1]
      - y_coeff[i][2];
//  logMessage(2, "y-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//           i, y_coeff[i][0], y_coeff[i][1], y_coeff[i][2], y_coeff[i][3]);

  for (i = 0; i < n_pts; i++)
    delete[] LHS[i];
  delete[] LHS;
  delete[] RHS;
}

void
GRInterp::compute_coeffs_notAKnot(double (*x_coeff)[4], double (*y_coeff)[4])
{

  // Computes the coefficients of the "not-a-knot" spline.
  // That is the spline with continuous third derivatives at the first interior
  // knot.

  double** LHS;
  double* RHS = new double[n_pts];

  LHS = new double*[n_pts];

  int i;
  for (i = 0; i < n_pts; i++)
    LHS[i] = new double[3];

  // This system is actually set up to solve for the linear coefficients.
  // Each spline is of the form a + b*t + c*t^2 + d*t^3.  There are two
  // conditions for each segment (f is the function being splined):
  //
  // A: f(i) = a_i
  // B: f(i+1) - f(i) = a_i + b_i + c_i + d_i - a_i
  //
  // and two shared conditions at each end:
  //
  // C: f'(i+1) = b_i + 2*c_i + 3*d_i = b_i+1
  // D: f''(i+1)/2 = c_i + 3*d_i = c_i+1
  //
  // Using eqn D to eliminate d_i:
  //
  // E: 3*(f(i+1)-f(i)) = 3*b_i + 2*c_i + c_i+1
  // F:         b_i+1 =   b_i +   c_i + c_i+1
  //
  // Writing eqn E for the next segment and adding to eqn E yields:
  //
  // G: 3*(f(i+2)-f(i)) = 3*b_i + 3*b_i+1 + 2*c_i + 3*c_i+1 + c_i+2
  //
  // Using eqn F for both i and i+1 to eliminate the c_i's, we finally get:
  //
  // 3*(f(i+2)-f(i)) = b_i + 4*b_i+1 + b_i+2
  //
  // The boundary condition at the left can be derived by combining eqn B@1 and
  // @2; eqn C@1; eqn D@2; and the condition that d_1 = d_2.  After eliminating
  // the c_i and d_i, we get:
  //
  // 2*b_1 + 4*b_2 = f(2) + 4*f(1) - 5*f(0)
  //
  // and similarly at the right side.
  //
  // Once we know all the b_i, we can easily use eqns E and F to get c_i and
  // then eqn C to get d_i.

  //Computing x coefficients
  LHS[0][0] = 0;
  LHS[0][1] = 2;
  LHS[0][2] = 4;
  RHS[0] = points[2].x() + 4 * points[1].x() - 5 * points[0].x();

  LHS[n_pts - 2][0] = 2.;
  LHS[n_pts - 2][1] = 4.;
  LHS[n_pts - 2][2] = 0.;
  RHS[n_pts - 2] = points[n_pts - 1].x() + 4 * points[n_pts - 2].x()
      - 5 * points[n_pts - 3].x();

  for (i = 1; i < n_pts - 2; i++) {
    LHS[i][0] = 1.;
    LHS[i][1] = 4.;
    LHS[i][2] = 1.;
    RHS[i] = 3. * (points[i + 1].x() - points[i - 1].x());
  }

//  for (i = 0; i < n_pts - 1; i++) {
//    logMessage(2, "%d LHS: %5f %5f %5f RHS: %5f\n", i,
//               LHS[i][0], LHS[i][1], LHS[i][2], RHS[i]);
//  }

  solveTriDiag(LHS, RHS, n_pts - 1);

  for (i = 0; i < n_pts - 2; i++) {
    x_coeff[i][0] = points[i].x();
    x_coeff[i][1] = RHS[i];
    x_coeff[i][2] = 3. * (points[i + 1].x() - points[i].x()) - 2. * RHS[i]
	- RHS[i + 1];
    x_coeff[i][3] = 2. * (points[i].x() - points[i + 1].x()) + RHS[i]
	+ RHS[i + 1];
//    logMessage(2, "x-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//             i, x_coeff[i][0], x_coeff[i][1], x_coeff[i][2], x_coeff[i][3]);
  }
  // Here we take advantage of the boundary condition and previous data.
  x_coeff[i][0] = points[i].x();
  x_coeff[i][1] = RHS[i];
  x_coeff[i][2] = x_coeff[i - 1][2] + 3 * x_coeff[i - 1][3];
  x_coeff[i][3] = x_coeff[i - 1][3];
//  logMessage(2, "x-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//           i, x_coeff[i][0], x_coeff[i][1], x_coeff[i][2], x_coeff[i][3]);

  //Getting y coefficients
  LHS[0][0] = 0;
  LHS[0][1] = 2;
  LHS[0][2] = 4;
  RHS[0] = points[2].y() + 4 * points[1].y() - 5 * points[0].y();

  LHS[n_pts - 2][0] = 2.;
  LHS[n_pts - 2][1] = 4.;
  LHS[n_pts - 2][2] = 0.;
  RHS[n_pts - 2] = points[n_pts - 1].y() + 4 * points[n_pts - 2].y()
      - 5 * points[n_pts - 3].y();

  for (i = 1; i < n_pts - 2; i++) {
    LHS[i][0] = 1.;
    LHS[i][1] = 4.;
    LHS[i][2] = 1.;
    RHS[i] = 3. * (points[i + 1].y() - points[i - 1].y());
  }

//  for (i = 0; i < n_pts - 1; i++) {
//    logMessage(2, "%d LHS: %5f %5f %5f RHS: %5f\n", i,
//               LHS[i][0], LHS[i][1], LHS[i][2], RHS[i]);
//  }

  solveTriDiag(LHS, RHS, n_pts - 1);

  for (i = 0; i < n_pts - 2; i++) {
    y_coeff[i][0] = points[i].y();
    y_coeff[i][1] = RHS[i];
    y_coeff[i][2] = 3. * (points[i + 1].y() - points[i].y()) - 2. * RHS[i]
	- RHS[i + 1];
    y_coeff[i][3] = 2. * (points[i].y() - points[i + 1].y()) + RHS[i]
	+ RHS[i + 1];
//    logMessage(2, "y-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//             i, y_coeff[i][0], y_coeff[i][1], y_coeff[i][2], y_coeff[i][3]);
  }
  // Here we take advantage of the boundary condition and previous data.
  y_coeff[i][0] = points[i].y();
  y_coeff[i][1] = RHS[i];
  y_coeff[i][2] = y_coeff[i - 1][2] + 3 * y_coeff[i - 1][3];
  y_coeff[i][3] = y_coeff[i - 1][3];
//  logMessage(2, "y-coeffs: %2d %20.12g %20.12g %20.12g %20.12g\n",
//           i, y_coeff[i][0], y_coeff[i][1], y_coeff[i][2], y_coeff[i][3]);

  for (i = 0; i < n_pts; i++)
    delete[] LHS[i];
  delete[] LHS;
  delete[] RHS;

}

void
GRInterp::compute_coeffs_closed(double (*x_coeff)[4], double (*y_coeff)[4])
{

  if (x_coeff && y_coeff) {
  }
  assert(0);

  //Not supported yet, but have to write a solver for
  //a matrix of the form
  //[4 1  ...   1]
  //[1 4 1       ]
  //[0 1 4 1     ]
  //[     ...    ]
  //[1    ... 1 4]
  //The RHS is also slightly different

}
