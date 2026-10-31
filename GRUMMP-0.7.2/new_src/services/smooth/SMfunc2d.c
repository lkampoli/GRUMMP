#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

void
SMcomputeQuadCosines(double *vtx1, double *vtx2, double *vtx3, double *function,
		     int *num_values, double *vtx4)
{
  ////////////QUAD/////////////////
  /*
   *                                    5
   * 		vtx1 -------------------------------- vtx4
   *                  /   \                           /     /
   *                 /       \                  /          /
   *                /           \          /6             /
   * 				 /               \   /                 /
   *           3  /              /    \                 /4
   *             /           /           \2            /
   *            /        /                  \         /
   *           /     /                         \     /
   *          /  /                                \ /
   *        vtx2 -------------------------------- vtx3
   *                             1
   *
   *
   *
   */
  double xlen, ylen;
  double L1, L2, L3, L4, L5, L6;

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;

  xlen = vtx3[XDIR] - vtx4[XDIR];
  ylen = vtx3[YDIR] - vtx4[YDIR];
  L4 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx4[XDIR];
  ylen = vtx1[YDIR] - vtx4[YDIR];
  L5 = xlen * xlen + ylen * ylen;

  xlen = vtx2[XDIR] - vtx4[XDIR];
  ylen = vtx2[YDIR] - vtx4[YDIR];
  L6 = xlen * xlen + ylen * ylen;

  function[0] = .5 * (L3 + L5 - L6) / sqrt(L3 * L5); //v1
  function[1] = .5 * (L1 + L3 - L2) / sqrt(L1 * L3); //v2
  function[2] = .5 * (L1 + L4 - L6) / sqrt(L1 * L4); //v3
  function[3] = .5 * (L5 + L4 - L2) / sqrt(L5 * L4); //v4
  *num_values = 4;
}

void
SMcomputeTriCosines(double *vtx1, double *vtx2, double *vtx3, double *function,
		    int *num_values)
{
  double xlen, ylen;
  double L1, L2, L3;

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;

  function[0] = .5 * (L2 + L3 - L1) / (sqrt(L2 * L3));
  function[1] = .5 * (L1 + L3 - L2) / (sqrt(L1 * L3));
  function[2] = .5 * (L1 + L2 - L3) / (sqrt(L1 * L2));
  *num_values = 3;
}

void
SMcomputeInteriorTriCosines(double *vtx1, double *vtx2, double *vtx3,
			    double *function, int *num_values)
{
  double xlen, ylen;
  double L1, L2, L3;

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;

  function[0] = .5 * (L2 + L3 - L1) / (sqrt(L2 * L3));
  *num_values = 1;
}

void
SMcomputeNegTriCosines(double *vtx1, double *vtx2, double *vtx3,
		       double *function, int *num_values)
{
  int i;
  SMcomputeTriCosines(vtx1, vtx2, vtx3, function, num_values);
  for (i = 0; i < *num_values; i++)
    function[i] = -function[i];
}

void
SMcomputeTriAngles(double *vtx1, double *vtx2, double *vtx3, double *function,
		   int *num_values)
{
  int i;

  SMcomputeTriCosines(vtx1, vtx2, vtx3, function, num_values);
  for (i = 0; i < *num_values; i++)
    function[i] = SAFE_ACOS(function[i]);
}

void
SMcomputeInteriorTriAngles(double *vtx1, double *vtx2, double *vtx3,
			   double *function, int *num_values)
{
  int i;

  SMcomputeInteriorTriCosines(vtx1, vtx2, vtx3, function, num_values);
  for (i = 0; i < *num_values; i++)
    function[i] = SAFE_ACOS(function[i]);
}

void
SMcomputeNegTriAngles(double *vtx1, double *vtx2, double *vtx3,
		      double *function, int *num_values)
{
  int i;
  SMcomputeTriAngles(vtx1, vtx2, vtx3, function, num_values);
  for (i = 0; i < *num_values; i++)
    function[i] = -function[i];
}

void
SMcomputeQuadSines(double *vtx1, double *vtx2, double *vtx3, double *function,
		   int *num_values, double *vtx4)
{
  double xlen, ylen;
  double L1, L2, L3, L4, L5, L6;
  double sL1, sL3, sL4, sL5;
  double cos0;
////////////QUAD/////////////////
  /*
   *                                    5
   * 					vtx1 -------------------------------- vtx4
   *                  /   \                           /     /
   *                 /       \                  /          /
   *                /           \          /6             /
   * 				 /               \   /                 /
   *           3  /              /    \                 /4
   *             /           /           \2            /
   *            /        /                  \         /
   *           /     /                         \     /
   *          /  /                                \ /
   *        vtx2 -------------------------------- vtx3
   *                             1
   *
   *
   *
   */
  /* Computes the sine of the angles of the triangle, but only returns
   positive values - as it stands this function cannot be used for mesh
   untangling */

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;
  sL1 = sqrt(L1);

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);

  xlen = vtx3[XDIR] - vtx4[XDIR];
  ylen = vtx3[YDIR] - vtx4[YDIR];
  L4 = xlen * xlen + ylen * ylen;
  sL4 = sqrt(L4);

  xlen = vtx1[XDIR] - vtx4[XDIR];
  ylen = vtx1[YDIR] - vtx4[YDIR];
  L5 = xlen * xlen + ylen * ylen;
  sL5 = sqrt(L5);

  xlen = vtx2[XDIR] - vtx4[XDIR];
  ylen = vtx2[YDIR] - vtx4[YDIR];
  L6 = xlen * xlen + ylen * ylen;

  cos0 = .5 * (L3 + L5 - L6) / (sL3 * sL5); //v1
  function[0] = sqrt(1.0 - cos0 * cos0);
  cos0 = .5 * (L1 + L3 - L2) / (sL1 * sL3); //v2
  function[1] = sqrt(1.0 - cos0 * cos0);
  cos0 = .5 * (L1 + L4 - L6) / (sL1 * sL4); //v3
  function[2] = sqrt(1.0 - cos0 * cos0);
  cos0 = .5 * (L5 + L4 - L2) / (sL5 * sL4); //v4
  function[3] = sqrt(1.0 - cos0 * cos0);
  *num_values = 4;
}

void
SMcomputeTriSines(double *vtx1, double *vtx2, double *vtx3, double *function,
		  int *num_values)
{
  double xlen, ylen;
  double L1, L2, L3;
  double sL1, sL2, sL3;
  double cos0;

  /* Computes the sine of the angles of the triangle, but only returns 
   positive values - as it stands this function cannot be used for mesh
   untangling */

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;
  sL1 = sqrt(L1);

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);

  cos0 = .5 * (L2 + L3 - L1) / (sL2 * sL3);
  function[0] = sqrt(1.0 - cos0 * cos0);
  function[1] = (sL2 * function[0]) / sL1;
  function[2] = (sL3 * function[0]) / sL1;
  *num_values = 3;
}

void
SMcomputeInteriorEdgeratio(double *vtx1, double *vtx2, double *vtx3,
			   double *function, int *num_values)
{
  double xlen, ylen;
  double L2, L3;
  double sL2, sL3;
//	double cos0;
  double skewness;
//	double sin0;
  *num_values = 1;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);
//    cos0 = .5*(L2+L3-L1)/(sL2*sL3);
//    skewness = L2*L3*2./(L2*L2+L3*L3);
  skewness = sL2 * sL3 * 2. / (L2 + L3);
  function[0] = skewness;
//    function[0] = skewness;
}
void
SMcomputeInteriorSquareness(double *vtx1, double *vtx2, double *vtx3,
			    double *function, int *num_values)
{
  double xlen, ylen;
  double L1, L2, L3;
  double sL2, sL3;
  double cos0;
  double skewness;
  double sin0;
  *num_values = 1;
  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);
  cos0 = .5 * (L2 + L3 - L1) / (sL2 * sL3);
//    skewness = L2*L3*2./(L2*L2+L3*L3);
  skewness = sL2 * sL3 * 2. / (L2 + L3);
  sin0 = sqrt(1.0 - cos0 * cos0);
  function[0] = sin0 * skewness;
//    function[0] = skewness;
}
void
SMcomputeInteriorTriSines(double *vtx1, double *vtx2, double *vtx3,
			  double *function, int *num_values)
{
  double xlen, ylen;
  double L1, L2, L3;
  double sL2, sL3;
  double cos0;

  /* Computes the sine of the angles of the triangle, but only returns 
   positive values - as it stands this function cannot be used for mesh
   untangling */

  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
  L1 = xlen * xlen + ylen * ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);

  cos0 = .5 * (L2 + L3 - L1) / (sL2 * sL3);
  function[0] = sqrt(1.0 - cos0 * cos0);
//    function[0] = function[0]*function[0];
  *num_values = 1;
}

void
SMcomputeTriJacobians(double *vtx1, double *vtx2, double *vtx3,
		      double *function, int *num_values)
{
  double a, b, c;
  double equal_jac;
  double L1_2;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* Computes the square of the difference between the actual jacbian
   and the jacobian of an equilateral triangle based on the edge opposite
   the free vertex.  We need to minimize the maximum difference and 
   the function is therefore negated.  The corresponing gradient code 
   is in SMcomputeJacobianGradients.  There is one value per triangle */

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L1_2 = (xx2 - xx3) * (xx2 - xx3) + (yy2 - yy3) * (yy2 - yy3);
  equal_jac = .5 * sqrt(3.) * L1_2;

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;
  function[0] = a + b * xx1 + c * yy1;
  function[0] = -((function[0] - equal_jac) / equal_jac)
      * ((function[0] - equal_jac) / equal_jac);

  *num_values = 1;
}

void
SMcomputeTriArea(double *vtx1, double *vtx2, double *vtx3, double *function,
		 int *num_values)
{
  double a, b, c;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;

  function[0] = .5 * (a + b * xx1 + c * yy1);
  *num_values = 1;
}

void
SMcomputeScaledTriJacobians(double *vtx1, double *vtx2, double *vtx3,
			    double *function, int *num_values)
{
  double a, b, c;
  double L1, L2, L3;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* This function computes the scaled jacobian for each angle of the triangle
   -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
   Note that this is equivalent to the sine of the angle but is an efficient way to compute
   both positive and negative sines for mesh untangling.  There are 3 function values 
   per triangle and we need to maximize the minimum scaled jacobian.  The 
   corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L1 = sqrt((xx2 - xx3) * (xx2 - xx3) + (yy2 - yy3) * (yy2 - yy3));
  L2 = sqrt((xx1 - xx3) * (xx1 - xx3) + (yy1 - yy3) * (yy1 - yy3));
  L3 = sqrt((xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1));

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;
  function[0] = a + b * xx1 + c * yy1;
  function[0] = function[0] / (L2 * L3);

  a = xx3 * yy1 - xx1 * yy3;
  b = yy3 - yy1;
  c = xx1 - xx3;
  function[1] = a + b * xx2 + c * yy2;
  function[1] = function[1] / (L1 * L3);

  a = xx1 * yy2 - xx2 * yy1;
  b = yy1 - yy2;
  c = xx2 - xx1;
  function[2] = a + b * xx3 + c * yy3;
  function[2] = function[2] / (L1 * L2);

  *num_values = 3;
}

void
SMcomputeInteriorScaledTriJacobians(double *vtx1, double *vtx2, double *vtx3,
				    double *function, int *num_values)
{
  double a, b, c;
  double L2, L3;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* This function computes the scaled jacobian for each angle of the triangle
   -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
   Note that this is equivalent to the sine of the angle but is an efficient way to compute
   both positive and negative sines for mesh untangling.  There are 3 function values 
   per triangle and we need to maximize the minimum scaled jacobian.  The 
   corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L2 = sqrt((xx1 - xx3) * (xx1 - xx3) + (yy1 - yy3) * (yy1 - yy3));
  L3 = sqrt((xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1));

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;
  function[0] = a + b * xx1 + c * yy1;
  function[0] = function[0] / (L2 * L3);

  *num_values = 1;
}

void
SMcomputeAreaLengthRatio(double *vtx1, double *vtx2, double *vtx3,
			 double *function, int *num_values)
{
  double a, b, c;
  double L1, L2, L3;
  double L1_xdiff, L2_xdiff, L3_xdiff;
  double L1_ydiff, L2_ydiff, L3_ydiff;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* comptes the SMcomputeAreaLengthRatio function which
   gives the ratio of the the area of the triangle and the sum of the squares of 
   the length of the edges */

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L1_xdiff = xx2 - xx3;
  L1_ydiff = yy2 - yy3;
  L2_xdiff = xx1 - xx3;
  L2_ydiff = yy1 - yy3;
  L3_xdiff = xx1 - xx2;
  L3_ydiff = yy1 - yy2;

  L1 = L1_xdiff * L1_xdiff + L1_ydiff * L1_ydiff;
  L2 = L2_xdiff * L2_xdiff + L2_ydiff * L2_ydiff;
  L3 = L3_xdiff * L3_xdiff + L3_ydiff * L3_ydiff;

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;
  function[0] = .5 * (a + b * xx1 + c * yy1);
  function[0] = 12. / sqrt(3.) * function[0] / (L1 + L2 + L3);

  *num_values = 1;
}

void
SMcomputeLengthAreaRatio(double *vtx1, double *vtx2, double *vtx3,
			 double *function, int *num_values)
{
  int i;
  SMcomputeAreaLengthRatio(vtx1, vtx2, vtx3, function, num_values);
  for (i = 0; i < (*num_values); i++)
    function[i] = -1 / function[i];
}

void
SMcomputeAngGradients(double *vtx1, double *vtx2, double *vtx3,
		      double (*gradient)[2], int *num_values)
{
  int i;
  int num_vert = 3;
  double *cosines, tmp;

  MY_MALLOC(cosines, (double *), sizeof(double) * num_vert, 1);

  SMcomputeTriCosines(vtx1, vtx2, vtx3, cosines, num_values);
  SMcomputeCosGradients(vtx1, vtx2, vtx3, gradient, num_values);

  /* compute the gradient of the angle */
  for (i = 0; i < *num_values; i++) {
    tmp = -1. / sqrt(1 - cosines[i] * cosines[i]);
    gradient[i][XDIR] *= tmp;
    gradient[i][YDIR] *= tmp;
  }
  *num_values = num_vert;
  MY_FREE(cosines);
}

void
SMcomputeInteriorAngGradients(double *vtx1, double *vtx2, double *vtx3,
			      double (*gradient)[2], int *num_values)
{
  int i;
  int num_vert = 3;
  double *cosines, tmp;

  MY_MALLOC(cosines, (double *), sizeof(double) * num_vert, 1);

  SMcomputeInteriorTriCosines(vtx1, vtx2, vtx3, cosines, num_values);
  SMcomputeInteriorCosGradients(vtx1, vtx2, vtx3, gradient, num_values);

  /* compute the gradient of the angle */
  for (i = 0; i < *num_values; i++) {
    tmp = -1. / sqrt(1 - cosines[i] * cosines[i]);
    gradient[i][XDIR] *= tmp;
    gradient[i][YDIR] *= tmp;
  }
  *num_values = num_vert;
  MY_FREE(cosines);
}

void
SMcomputeNegAngGradients(double *vtx1, double *vtx2, double *vtx3,
			 double (*gradient)[2], int *num_values)
{
  int i, j;

  SMcomputeAngGradients(vtx1, vtx2, vtx3, gradient, num_values);
  for (i = 0; i < *num_values; i++) {
    for (j = 0; j < 2; j++)
      gradient[i][j] = -gradient[i][j];
  }
}

void
SMcomputeCosGradientsQuad(double *vtx1, double *vtx2, double *vtx3,
			  double (*gradient)[2], int *num_values, double *vtx4)
{
//    int i;
  int num_vert = 4;
  double x_1, x_2, x_3, x_4, y_1, y_2, y_3, y_4;
  double a, c, d, e;
  double a2, b2, c2, d2, e2, f2;
  double gx, gy;
  double term1, term2;

  x_1 = vtx1[XDIR];
  x_2 = vtx2[XDIR];
  x_3 = vtx3[XDIR];
  x_4 = vtx4[XDIR];
  y_1 = vtx1[YDIR];
  y_2 = vtx2[YDIR];
  y_3 = vtx3[YDIR];
  y_4 = vtx4[XDIR];
  ////////////QUAD/////////////////
  /*
   *                                    e
   * 					vtx1 -------------------------------- vtx4
   *                  /   \                           /     /
   *                 /       \                  /          /
   *                /           \          /f             /
   * 				 /               \   /                 /
   *           a  /              /    \                 /d
   *             /           /           \b            /
   *            /        /                  \         /
   *           /     /                         \     /
   *          /  /                                \ /
   *        vtx2 -------------------------------- vtx3
   *                             c
   *
   *
   */
  a2 = (x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1);
  b2 = (x_1 - x_3) * (x_1 - x_3) + (y_1 - y_3) * (y_1 - y_3);
  c2 = (x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2) * (y_3 - y_2);
  d2 = (x_3 - x_4) * (x_3 - x_4) + (y_3 - y_4) * (y_3 - y_4);
  e2 = (x_1 - x_4) * (x_1 - x_4) + (y_1 - y_4) * (y_1 - y_4);
  f2 = (x_4 - x_2) * (x_4 - x_2) + (y_4 - y_2) * (y_4 - y_2);

  a = sqrt(a2);
  c = sqrt(c2);
  d = sqrt(d2);
  e = sqrt(e2);

  /* compute the gradient of the cosine */

  //gradient of cosine 1
  term1 = (a2 + f2 - e2) / (2 * a2 * e) * (1 / a);
  term2 = (e2 + f2 - a2) / (2 * a * e2) * (1 / e);

  gx = term1 * (x_1 - x_2) + term2 * (x_1 - x_4);
  gy = term1 * (y_1 - y_2) + term2 * (y_1 - y_4);

  gradient[0][XDIR] = gx;
  gradient[0][YDIR] = gy;

  //gradient of cosine 2
  term1 = (a2 + b2 - c2) / (2 * a2 * c) * (1 / a);
  term2 = (c2 + b2 - a2) / (2 * a * c2) * (1 / c);

  gx = term1 * (x_2 - x_1) + term2 * (x_2 - x_3);
  gy = term1 * (y_2 - y_1) + term2 * (y_2 - y_3);

  gradient[1][XDIR] = gx;
  gradient[1][YDIR] = gy;

  //gradient of cosine 3
  term1 = (c2 + f2 - d2) / (2 * c2 * d) * (1 / c);
  term2 = (d2 + f2 - c2) / (2 * c * d2) * (1 / d);

  gx = term1 * (x_3 - x_2) + term2 * (x_3 - x_4);
  gy = term1 * (y_3 - y_2) + term2 * (y_3 - y_4);

  gradient[2][XDIR] = gx;
  gradient[2][YDIR] = gy;
  //gradient of cosine 4
  term1 = (d2 + b2 - e2) / (2 * d2 * e) * (1 / d);
  term2 = (e2 + b2 - d2) / (2 * d * e2) * (1 / e);

  gx = term1 * (x_4 - x_3) + term2 * (x_4 - x_1);
  gy = term1 * (y_4 - y_3) + term2 * (y_4 - y_1);

  gradient[3][XDIR] = gx;
  gradient[3][YDIR] = gy;

  *num_values = num_vert;
}

void
SMcomputeCosGradients(double *vtx1, double *vtx2, double *vtx3,
		      double (*gradient)[2], int *num_values)
{
  int i;
  int num_vert = 3;
  double x_1, x_2, x_3, y_1, y_2, y_3;
  double a, b, c;
  double a2, b2, c2;
  double gx, gy;
  double term1, term2;

  x_1 = vtx1[XDIR];
  x_2 = vtx2[XDIR];
  x_3 = vtx3[XDIR];
  y_1 = vtx1[YDIR];
  y_2 = vtx2[YDIR];
  y_3 = vtx3[YDIR];

  a2 = (x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1);
  b2 = (x_1 - x_3) * (x_1 - x_3) + (y_1 - y_3) * (y_1 - y_3);
  c2 = (x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2) * (y_3 - y_2);

  a = sqrt(a2);
  b = sqrt(b2);
  c = sqrt(c2);

  /* compute the gradient of the cosine */
  for (i = 0; i < num_vert; i++) {
    if (i == 0) {
      term1 = (a2 + c2 - b2) / (2 * a2 * b) * (1 / a);
      term2 = (b2 + c2 - a2) / (2 * a * b2) * (1 / b);
      /* SM_LOG_FLOPS(__SM_GRADIENT__,14); */
    }
    else if (i == 1) {
      term1 = (a2 + b2 - c2) / (2 * a2 * c) * (1 / a);
      term2 = -1 / (a * c);
      /* SM_LOG_FLOPS(__SM_GRADIENT__,9); */
    }
    else {
      term1 = -1 / (b * c);
      term2 = (b2 + a2 - c2) / (2 * b2 * c) * (1 / b);
      /* SM_LOG_FLOPS(__SM_GRADIENT__,9); */
    }

    gx = term1 * (x_1 - x_2) + term2 * (x_1 - x_3);
    gy = term1 * (y_1 - y_2) + term2 * (y_1 - y_3);
    /* SM_LOG_FLOPS(__SM_GRADIENT__,10); */

    gradient[i][XDIR] = gx;
    gradient[i][YDIR] = gy;
  }
  *num_values = num_vert;
}

void
SMcomputeInteriorCosGradients(double *vtx1, double *vtx2, double *vtx3,
			      double (*gradient)[2], int *num_values)
{
  double x_1, x_2, x_3, y_1, y_2, y_3;
  double a, b;
  double a2, b2, c2;
  double gx, gy;
  double term1, term2;

  x_1 = vtx1[XDIR];
  x_2 = vtx2[XDIR];
  x_3 = vtx3[XDIR];
  y_1 = vtx1[YDIR];
  y_2 = vtx2[YDIR];
  y_3 = vtx3[YDIR];

  a2 = (x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1);
  b2 = (x_1 - x_3) * (x_1 - x_3) + (y_1 - y_3) * (y_1 - y_3);
  c2 = (x_3 - x_2) * (x_3 - x_2) + (y_3 - y_2) * (y_3 - y_2);

  a = sqrt(a2);
  b = sqrt(b2);

  /* compute the gradient of the cosine */
  term1 = (a2 + c2 - b2) / (2 * a2 * b) * (1 / a);
  term2 = (b2 + c2 - a2) / (2 * a * b2) * (1 / b);

  gx = term1 * (x_1 - x_2) + term2 * (x_1 - x_3);
  gy = term1 * (y_1 - y_2) + term2 * (y_1 - y_3);

  gradient[0][XDIR] = gx;
  gradient[0][YDIR] = gy;
  *num_values = 1;
}

void
SMcomputeNegCosGradients(double *vtx1, double *vtx2, double *vtx3,
			 double (*gradient)[2], int *num_values)
{
  int i, j;

  SMcomputeCosGradients(vtx1, vtx2, vtx3, gradient, num_values);
  for (i = 0; i < *num_values; i++) {
    for (j = 0; j < 2; j++) {
      gradient[i][j] = -gradient[i][j];
    }
  }
}

void
SMcomputeSineGradientsQuad(double *vtx1, double *vtx2, double *vtx3,
			   double (*gradient)[2], int *num_values, double *vtx4)
{
  double sqrtcos0, sqrtcos1, sqrtcos2, sqrtcos3;
  double *cosine;
  double (*cos_gradient)[2];
  int num_vert;

  num_vert = 4;

  MY_MALLOC(cosine, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(cos_gradient, (double(*)[2]), sizeof(double) * num_vert * 2, 1);

  SMcomputeQuadCosines(vtx1, vtx2, vtx3, cosine, num_values, vtx4);
  SMcomputeCosGradientsQuad(vtx1, vtx2, vtx3, cos_gradient, num_values, vtx4);

  sqrtcos0 = sqrt(1. - cosine[0] * cosine[0]);
  sqrtcos1 = sqrt(1. - cosine[1] * cosine[1]);
  sqrtcos2 = sqrt(1. - cosine[2] * cosine[2]);
  sqrtcos3 = sqrt(1. - cosine[3] * cosine[3]);
  gradient[0][XDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][XDIR];
  gradient[0][YDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][YDIR];
  gradient[1][XDIR] = (-1. / sqrtcos1) * cosine[1] * cos_gradient[1][XDIR];
  gradient[1][YDIR] = (-1. / sqrtcos1) * cosine[1] * cos_gradient[1][YDIR];
  gradient[2][XDIR] = (-1. / sqrtcos2) * cosine[2] * cos_gradient[2][XDIR];
  gradient[2][YDIR] = (-1. / sqrtcos2) * cosine[2] * cos_gradient[2][YDIR];
  gradient[3][XDIR] = (-1. / sqrtcos3) * cosine[3] * cos_gradient[3][XDIR];
  gradient[3][YDIR] = (-1. / sqrtcos3) * cosine[3] * cos_gradient[3][YDIR];
  MY_FREE(cosine);
  MY_FREE(cos_gradient);

  *num_values = 4;
}

void
SMcomputeSineGradients(double *vtx1, double *vtx2, double *vtx3,
		       double (*gradient)[2], int *num_values)
{
  double sqrtcos0, sqrtcos1, sqrtcos2;
  double *cosine;
  double (*cos_gradient)[2];
  int num_vert;

  num_vert = 3;

  MY_MALLOC(cosine, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(cos_gradient, (double(*)[2]), sizeof(double) * num_vert * 2, 1);

  SMcomputeTriCosines(vtx1, vtx2, vtx3, cosine, num_values);
  SMcomputeCosGradients(vtx1, vtx2, vtx3, cos_gradient, num_values);

  sqrtcos0 = sqrt(1. - cosine[0] * cosine[0]);
  sqrtcos1 = sqrt(1. - cosine[1] * cosine[1]);
  sqrtcos2 = sqrt(1. - cosine[2] * cosine[2]);
  gradient[0][XDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][XDIR];
  gradient[0][YDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][YDIR];
  gradient[1][XDIR] = (-1. / sqrtcos1) * cosine[1] * cos_gradient[1][XDIR];
  gradient[1][YDIR] = (-1. / sqrtcos1) * cosine[1] * cos_gradient[1][YDIR];
  gradient[2][XDIR] = (-1. / sqrtcos2) * cosine[2] * cos_gradient[2][XDIR];
  gradient[2][YDIR] = (-1. / sqrtcos2) * cosine[2] * cos_gradient[2][YDIR];

  MY_FREE(cosine);
  MY_FREE(cos_gradient);

  *num_values = 3;
}

void
SMcomputeInteriorEdgeratioGradients(double *vtx1, double *vtx2, double *vtx3,
				    double (*gradient)[2], int *num_values)
{
//    double sqrtcos0;
  double *cosine;
  double *squareness;
  double (*cos_gradient)[2];
  int num_vert;

  double xlen, ylen;
  double L2, L3;
  double sL2, sL3;
  xlen = vtx2[XDIR] - vtx3[XDIR];
  ylen = vtx2[YDIR] - vtx3[YDIR];
//    L1 = xlen*xlen + ylen*ylen;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);
  num_vert = 1;

  MY_MALLOC(cosine, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(squareness, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(cos_gradient, (double(*)[2]), sizeof(double) * num_vert * 2, 1);

//    double term1 = (4*sL2*pow(sL3,6)-4*pow(sL2,5)*pow(sL3,2))/(pow((L2*L2+L3*L3),2));
//    double term2 = (4*sL3*pow(sL2,6)-4*pow(sL3,5)*pow(sL2,2))/(pow((L2*L2+L3*L3),2));
  double term1 = (2 * pow(sL3, 3) - 2 * pow(sL2, 2) * sL3)
      / (pow((L2 + L3), 2));
  double term2 = (2 * pow(sL2, 3) - 2 * pow(sL3, 2) * sL2)
      / (pow((L2 + L3), 2));
  gradient[0][XDIR] = (term1 / sL2 * (vtx1[XDIR] - vtx3[XDIR])
      + term2 / sL3 * (vtx1[XDIR] - vtx2[XDIR]));
  gradient[0][YDIR] = (term1 / sL2 * (vtx1[YDIR] - vtx3[YDIR])
      + term2 / sL3 * (vtx1[YDIR] - vtx2[YDIR]));
//    gradient[0][XDIR] = (term1/sL2*(vtx1[XDIR] - vtx3[XDIR])+term2/sL3*(vtx1[XDIR] - vtx2[XDIR]));
//    gradient[0][YDIR] = (term1/sL2*(vtx1[YDIR] - vtx3[YDIR])+term2/sL3*(vtx1[YDIR] - vtx2[YDIR]));
//    printf("squareness: %f\n", squareness[0]);
  MY_FREE(cosine);
  MY_FREE(cos_gradient);
  MY_FREE(squareness);

  *num_values = 1;
}

void
SMcomputeInteriorSquarenessGradients(double *vtx1, double *vtx2, double *vtx3,
				     double (*gradient)[2], int *num_values)
{
  double sqrtcos0;
  double *cosine;
  double *squareness;
  double (*cos_gradient)[2];
  int num_vert;

  double xlen, ylen;
  double L2, L3;
  double sL2, sL3;

  xlen = vtx1[XDIR] - vtx3[XDIR];
  ylen = vtx1[YDIR] - vtx3[YDIR];
  L2 = xlen * xlen + ylen * ylen;
  sL2 = sqrt(L2);

  xlen = vtx1[XDIR] - vtx2[XDIR];
  ylen = vtx1[YDIR] - vtx2[YDIR];
  L3 = xlen * xlen + ylen * ylen;
  sL3 = sqrt(L3);
  num_vert = 1;

  MY_MALLOC(cosine, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(squareness, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(cos_gradient, (double(*)[2]), sizeof(double) * num_vert * 2, 1);

  SMcomputeInteriorTriCosines(vtx1, vtx2, vtx3, cosine, num_values);
  SMcomputeInteriorSquareness(vtx1, vtx2, vtx3, squareness, num_values);
  SMcomputeInteriorCosGradients(vtx1, vtx2, vtx3, cos_gradient, num_values);

  sqrtcos0 = sqrt(1. - cosine[0] * cosine[0]);
//    double term1 = (4*sL2*pow(sL3,6)-4*pow(sL2,5)*pow(sL3,2))/(pow((L2*L2+L3*L3),2));
//    double term2 = (4*sL3*pow(sL2,6)-4*pow(sL3,5)*pow(sL2,2))/(pow((L2*L2+L3*L3),2));
  double term1 = (2 * pow(sL3, 3) - 2 * pow(sL2, 2) * sL3)
      / (pow((L2 + L3), 2));
  double term2 = (2 * pow(sL2, 3) - 2 * pow(sL3, 2) * sL2)
      / (pow((L2 + L3), 2));
  gradient[0][XDIR] = squareness[0] / sqrtcos0 * (-1. / sqrtcos0) * cosine[0]
      * cos_gradient[0][XDIR]
      + sqrtcos0
	  * (term1 / sL2 * (vtx1[XDIR] - vtx3[XDIR])
	      + term2 / sL3 * (vtx1[XDIR] - vtx2[XDIR]));
  gradient[0][YDIR] = squareness[0] / sqrtcos0 * (-1. / sqrtcos0) * cosine[0]
      * cos_gradient[0][YDIR]
      + sqrtcos0
	  * (term1 / sL2 * (vtx1[YDIR] - vtx3[YDIR])
	      + term2 / sL3 * (vtx1[YDIR] - vtx2[YDIR]));
//    gradient[0][XDIR] = (term1/sL2*(vtx1[XDIR] - vtx3[XDIR])+term2/sL3*(vtx1[XDIR] - vtx2[XDIR]));
//    gradient[0][YDIR] = (term1/sL2*(vtx1[YDIR] - vtx3[YDIR])+term2/sL3*(vtx1[YDIR] - vtx2[YDIR]));
//    printf("squareness: %f\n", squareness[0]);
  MY_FREE(cosine);
  MY_FREE(cos_gradient);
  MY_FREE(squareness);

  *num_values = 1;
}

void
SMcomputeInteriorSineGradients(double *vtx1, double *vtx2, double *vtx3,
			       double (*gradient)[2], int *num_values)
{
  double sqrtcos0;
  double *cosine;
  double (*cos_gradient)[2];
  int num_vert;

  num_vert = 1;

  MY_MALLOC(cosine, (double *), sizeof(double) * num_vert, 1);
  MY_MALLOC(cos_gradient, (double(*)[2]), sizeof(double) * num_vert * 2, 1);

  SMcomputeInteriorTriCosines(vtx1, vtx2, vtx3, cosine, num_values);
  SMcomputeInteriorCosGradients(vtx1, vtx2, vtx3, cos_gradient, num_values);

  sqrtcos0 = sqrt(1. - cosine[0] * cosine[0]);
  gradient[0][XDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][XDIR];
  gradient[0][YDIR] = (-1. / sqrtcos0) * cosine[0] * cos_gradient[0][YDIR];
  MY_FREE(cosine);
  MY_FREE(cos_gradient);

  *num_values = 1;
}

void
SMcomputeJacobianGradients(double *vtx1, double *vtx2, double *vtx3,
			   double (*gradient)[2], int *num_values)
{
  double a, b, c;
  double equal_jac;
  double L1;

  L1 = (vtx2[XDIR] - vtx3[XDIR]) * (vtx2[XDIR] - vtx3[XDIR]);
  L1 = L1 + (vtx2[YDIR] - vtx3[YDIR]) * (vtx2[YDIR] - vtx3[YDIR]);
  equal_jac = .5 * sqrt(3.) * L1;

  a = vtx2[XDIR] * vtx3[YDIR] - vtx3[XDIR] * vtx2[YDIR];
  b = vtx2[YDIR] - vtx3[YDIR];
  c = vtx3[XDIR] - vtx2[XDIR];

  gradient[0][XDIR] = -2 * (a + b * vtx1[XDIR] + c * vtx1[YDIR] - equal_jac) * b
      / equal_jac;
  gradient[0][YDIR] = -2 * (a + b * vtx1[XDIR] + c * vtx1[YDIR] - equal_jac) * c
      / equal_jac;

  *num_values = 1;
}

void
SMcomputeScaledJacobianGradients(double *vtx1, double *vtx2, double *vtx3,
				 double (*gradient)[2], int *num_values)
{
  double f;
  double L1, L2, L3;
  double L1_2, L2_2, L3_2;
  double dL2, dL3;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* This function computes the scaled jacobian for each angle of the triangle
   -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
   Note that this is equivalent to the sine of the angle but is an efficient way to compute
   both positive and negative sines for mesh untangling.  There are 3 function values 
   per triangle and we need to maximize the minimum scaled jacobian.  The 
   corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L1_2 = (xx2 - xx3) * (xx2 - xx3) + (yy2 - yy3) * (yy2 - yy3);
  L1 = sqrt(L1_2);

  L2_2 = (xx1 - xx3) * (xx1 - xx3) + (yy1 - yy3) * (yy1 - yy3);
  L2 = sqrt(L2_2);
  dL2 = L2_2 * L2;

  L3_2 = (xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1);
  L3 = sqrt(L3_2);
  dL3 = L3_2 * L3;

  f = xx2 * yy3 - xx3 * yy2 + (yy2 - yy3) * xx1 + (xx3 - xx2) * yy1;
  gradient[0][XDIR] = (yy2 - yy3) / (L3 * L2) - f * (xx1 - xx2) / (dL3 * L2)
      - f * (xx1 - xx3) / (L3 * dL2);
  gradient[0][YDIR] = (xx3 - xx2) / (L3 * L2) - f * (yy1 - yy2) / (dL3 * L2)
      - f * (yy1 - yy3) / (L3 * dL2);

  f = xx3 * yy1 - xx1 * yy3 + (yy3 - yy1) * xx2 + (xx1 - xx3) * yy2;
  gradient[1][XDIR] = (yy2 - yy3) / (L1 * L3) - f * (xx1 - xx2) / (L1 * dL3);
  gradient[1][YDIR] = (xx3 - xx2) / (L1 * L3) - f * (yy1 - yy2) / (L1 * dL3);

  f = xx1 * yy2 - xx2 * yy1 + (yy1 - yy2) * xx3 + (xx2 - xx1) * yy3;
  gradient[2][XDIR] = (yy2 - yy3) / (L2 * L1) - f * (xx1 - xx3) / (dL2 * L1);
  gradient[2][YDIR] = (xx3 - xx2) / (L2 * L1) - f * (yy1 - yy3) / (dL2 * L1);

  *num_values = 3;
}

void
SMcomputeInteriorScaledJacobianGradients(double *vtx1, double *vtx2,
					 double *vtx3, double (*gradient)[2],
					 int *num_values)
{
  double f;
  /* Don't ever need L1 or L1_2, so don't declare or compute them. */
  double L2, L3;
  double L2_2, L3_2;
  double dL2, dL3;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* This function computes the scaled jacobian for each angle of the triangle
   -- Jac/(L1*L2) where L1 and L2 are the length of the edges containing the angle.
   Note that this is equivalent to the sine of the angle but is an efficient way to compute
   both positive and negative sines for mesh untangling.  There are 3 function values 
   per triangle and we need to maximize the minimum scaled jacobian.  The 
   corresponding gradient code is in the function SMcomputeScaledJacobianGradients*/

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L2_2 = (xx1 - xx3) * (xx1 - xx3) + (yy1 - yy3) * (yy1 - yy3);
  L2 = sqrt(L2_2);
  dL2 = L2_2 * L2;

  L3_2 = (xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1);
  L3 = sqrt(L3_2);
  dL3 = L3_2 * L3;

  f = xx2 * yy3 - xx3 * yy2 + (yy2 - yy3) * xx1 + (xx3 - xx2) * yy1;
  gradient[0][XDIR] = (yy2 - yy3) / (L3 * L2) - f * (xx1 - xx2) / (dL3 * L2)
      - f * (xx1 - xx3) / (L3 * dL2);
  gradient[0][YDIR] = (xx3 - xx2) / (L3 * L2) - f * (yy1 - yy2) / (dL3 * L2)
      - f * (yy1 - yy3) / (L3 * dL2);

  *num_values = 1;
}

void
SMcomputeAreaLengthRatioGradients(double *vtx1, double *vtx2, double *vtx3,
				  double (*gradient)[2], int *num_values)
{
  double a, b, c;
  double L1, L2, L3;
  double L1_xdiff, L2_xdiff, L3_xdiff;
  double L1_ydiff, L2_ydiff, L3_ydiff;
  double xx1, xx2, xx3, yy1, yy2, yy3;

  /* comptes the gradient of the SMcomputeAreaLengthRatio function which
   gives the ratio of the the area of the triangle and the sum of the squares of 
   the length of the edges */

  xx1 = vtx1[XDIR];
  yy1 = vtx1[YDIR];
  xx2 = vtx2[XDIR];
  yy2 = vtx2[YDIR];
  xx3 = vtx3[XDIR];
  yy3 = vtx3[YDIR];

  L1_xdiff = xx2 - xx3;
  L1_ydiff = yy2 - yy3;
  L2_xdiff = xx1 - xx3;
  L2_ydiff = yy1 - yy3;
  L3_xdiff = xx1 - xx2;
  L3_ydiff = yy1 - yy2;

  L1 = L1_xdiff * L1_xdiff + L1_ydiff * L1_ydiff;
  L2 = L2_xdiff * L2_xdiff + L2_ydiff * L2_ydiff;
  L3 = L3_xdiff * L3_xdiff + L3_ydiff * L3_ydiff;

  a = xx2 * yy3 - xx3 * yy2;
  b = yy2 - yy3;
  c = xx3 - xx2;
  gradient[0][XDIR] = .5 * b / (L1 + L2 + L3)
      - (.5 * (a + b * xx1 + c * yy1) * (4 * xx1 - 2 * xx3 - 2 * xx2))
	  / ((L1 + L2 + L3) * (L1 + L2 + L3));
  gradient[0][YDIR] = .5 * c / (L1 + L2 + L3)
      - (.5 * (a + b * xx1 + c * yy1) * (4 * yy1 - 2 * yy3 - 2 * yy2))
	  / ((L1 + L2 + L3) * (L1 + L2 + L3));
  gradient[0][XDIR] = 12. / sqrt(3.) * gradient[0][XDIR];
  gradient[0][YDIR] = 12. / sqrt(3.) * gradient[0][YDIR];

  *num_values = 1;
}

void
SMNormJacSquared2D(double *vtx1, double *vtx2, double *vtx3, double *function,
		   int *num_values)
{
  double a1[2], a2[2];

  /* the a matrix contains the Jacobian */
  a1[0] = vtx2[0] - vtx1[0];
  a1[1] = vtx2[1] - vtx1[1];

  a2[0] = vtx3[0] - vtx1[0];
  a2[1] = vtx3[1] - vtx1[1];

  function[0] = SMfrobenius_norm_squared2x2(a1, a2);
  function[0] *= -1.;

  *num_values = 1;
}

void
SMcomputeNormJacSquaredGradients2D(double *vtx1, double *vtx2, double *vtx3,
				   double (*gradient)[2], int *num_values)
{
  double a1[2], a2[2];

  /* the a matrix contains the Jacobian */
  a1[0] = vtx2[0] - vtx1[0];
  a1[1] = vtx2[1] - vtx1[1];

  a2[0] = vtx3[0] - vtx1[0];
  a2[1] = vtx3[1] - vtx1[1];

  gradient[0][0] = -2 * (a1[0] + a2[0]);
  gradient[0][1] = -2 * (a1[1] + a2[1]);
  *num_values = 1;
}

