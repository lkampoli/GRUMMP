#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "SMsmooth.h"
#include "SMdihed_func.h"

static void
vCross(const double adVecA[3], const double adVecB[3], double adResult[3])
{
  adResult[0] = adVecA[1] * adVecB[2] - adVecA[2] * adVecB[1];
  adResult[1] = adVecA[2] * adVecB[0] - adVecA[0] * adVecB[2];
  adResult[2] = adVecA[0] * adVecB[1] - adVecA[1] * adVecB[0];
}

static double
dMagnitude(const double adVec[3])
{
  return (sqrt(adVec[0] * adVec[0] + adVec[1] * adVec[1] + adVec[2] * adVec[2]));
}

static double
dDot(const double adVecA[3], const double adVecB[3])
{
  return adVecA[0] * adVecB[0] + adVecA[1] * adVecB[1] + adVecA[2] * adVecB[2];
}

static void
vUnitNormal(const double adCoord0[3], const double adCoord1[3],
	    const double adCoord2[3], double adResult[3])
{
  double adVecA[3], adVecB[3], dMag;
  int i;
  for (i = 0; i < 3; i++) {
    adVecA[i] = adCoord1[i] - adCoord0[i];
    adVecB[i] = adCoord2[i] - adCoord0[i];
  }

  vCross(adVecA, adVecB, adResult);

  dMag = dMagnitude(adResult);

  adResult[0] /= dMag;
  adResult[1] /= dMag;
  adResult[2] /= dMag;
}

void
vSineDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double adResult[6], int* const piNResult)
{
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  double adTemp[3];

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);

  vCross(adNormA, adNormB, adTemp);
  adResult[5] = dMagnitude(adTemp); /* Edge 23 */

  vCross(adNormA, adNormC, adTemp);
  adResult[4] = dMagnitude(adTemp); /* Edge 13 */

  vCross(adNormA, adNormD, adTemp);
  adResult[3] = dMagnitude(adTemp); /* Edge 12 */

  vCross(adNormB, adNormC, adTemp);
  adResult[2] = dMagnitude(adTemp); /* Edge 03 */

  vCross(adNormB, adNormD, adTemp);
  adResult[1] = dMagnitude(adTemp); /* Edge 02 */

  vCross(adNormC, adNormD, adTemp);
  adResult[0] = dMagnitude(adTemp); /* Edge 01 */
}

void
vCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
		 const double adCoord2[3], const double adCoord3[3],
		 double adResult[6], int* const piNResult)
{
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);

  adResult[5] = -dDot(adNormA, adNormB); /* Edge 23 */
  adResult[4] = -dDot(adNormA, adNormC); /* Edge 13 */
  adResult[3] = -dDot(adNormA, adNormD); /* Edge 12 */
  adResult[2] = -dDot(adNormB, adNormC); /* Edge 03 */
  adResult[1] = -dDot(adNormB, adNormD); /* Edge 02 */
  adResult[0] = -dDot(adNormC, adNormD); /* Edge 01 */
}

#define iFuzzyOne(a) (fabs((a)-1) < 1.e-12)

void
vSolids(const double adCoord0[3], const double adCoord1[3],
	const double adCoord2[3], const double adCoord3[3], double adResult[4],
	int* const piNResult)
{
  double adDiheds[6];
  int iNDiheds;
  vDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adDiheds, &iNDiheds);
  assert(iNDiheds == 6);

  /* Solid angle at vertex 0 */
  adResult[0] = adDiheds[0] + adDiheds[1] + adDiheds[2] - M_PI;
  /* Solid angle at vertex 1 */
  adResult[1] = adDiheds[0] + adDiheds[3] + adDiheds[4] - M_PI;
  /* Solid angle at vertex 2 */
  adResult[2] = adDiheds[1] + adDiheds[3] + adDiheds[5] - M_PI;
  /* Solid angle at vertex 3 */
  adResult[3] = adDiheds[2] + adDiheds[4] + adDiheds[5] - M_PI;

  *piNResult = 4;
}

void
vDihedrals(const double adCoord0[3], const double adCoord1[3],
	   const double adCoord2[3], const double adCoord3[3],
	   double adResult[6], int* const piNResult)
{
  double adNormA[3], adNormB[3], adNormC[3], adNormD[3];
  double dTemp;

  *piNResult = 6;
  /* Unit normals (inward pointing) for all four faces */
  vUnitNormal(adCoord2, adCoord1, adCoord3, adNormA);
  vUnitNormal(adCoord0, adCoord2, adCoord3, adNormB);
  vUnitNormal(adCoord1, adCoord0, adCoord3, adNormC);
  vUnitNormal(adCoord0, adCoord1, adCoord2, adNormD);

  assert(iFuzzyOne(dMagnitude(adNormA)));
  assert(iFuzzyOne(dMagnitude(adNormB)));
  assert(iFuzzyOne(dMagnitude(adNormC)));
  assert(iFuzzyOne(dMagnitude(adNormD)));

  dTemp = -dDot(adNormA, adNormB);
  adResult[5] = SAFE_ACOS(dTemp); /* Edge 23 */

  dTemp = -dDot(adNormA, adNormC);
  adResult[4] = SAFE_ACOS(dTemp); /* Edge 13 */

  dTemp = -dDot(adNormA, adNormD);
  adResult[3] = SAFE_ACOS(dTemp); /* Edge 12 */

  dTemp = -dDot(adNormB, adNormC);
  adResult[2] = SAFE_ACOS(dTemp); /* Edge 03 */

  dTemp = -dDot(adNormB, adNormD);
  adResult[1] = SAFE_ACOS(dTemp); /* Edge 02 */

  dTemp = -dDot(adNormC, adNormD);
  adResult[0] = SAFE_ACOS(dTemp); /* Edge 01 */
}

void
vNegateDihedrals(const double adCoord0[3], const double adCoord1[3],
		 const double adCoord2[3], const double adCoord3[3],
		 double adResult[6], int* const piNResult)
{
  int i;
  vDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adResult, piNResult);
  for (i = 0; i < *piNResult; i++)
    adResult[i] = -adResult[i];
}

void
vNegateCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
		       const double adCoord2[3], const double adCoord3[3],
		       double adResult[6], int* const piNResult)
{
  int i;
  vCosineDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adResult, piNResult);
  for (i = 0; i < *piNResult; i++)
    adResult[i] = -adResult[i];
}

void
vScaledJacobian(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
  double xx0, xx1, xx2, xx3;
  double yy0, yy1, yy2, yy3;
  double zz0, zz1, zz2, zz3;
  double L1, L2, L3, L4, L5, L6;
  double L1_squared, L2_squared, L3_squared;
  double L4_squared, L5_squared, L6_squared;
  double jacobian;

  xx0 = adCoord0[XDIR];
  yy0 = adCoord0[YDIR];
  zz0 = adCoord0[ZDIR];
  xx1 = adCoord1[XDIR];
  yy1 = adCoord1[YDIR];
  zz1 = adCoord1[ZDIR];
  xx2 = adCoord2[XDIR];
  yy2 = adCoord2[YDIR];
  zz2 = adCoord2[ZDIR];
  xx3 = adCoord3[XDIR];
  yy3 = adCoord3[YDIR];
  zz3 = adCoord3[ZDIR];

  L1_squared = (xx1 - xx0) * (xx1 - xx0) + (yy1 - yy0) * (yy1 - yy0)
      + (zz1 - zz0) * (zz1 - zz0);
  L2_squared = (xx2 - xx0) * (xx2 - xx0) + (yy2 - yy0) * (yy2 - yy0)
      + (zz2 - zz0) * (zz2 - zz0);
  L3_squared = (xx3 - xx0) * (xx3 - xx0) + (yy3 - yy0) * (yy3 - yy0)
      + (zz3 - zz0) * (zz3 - zz0);
  L4_squared = (xx1 - xx2) * (xx1 - xx2) + (yy1 - yy2) * (yy1 - yy2)
      + (zz1 - zz2) * (zz1 - zz2);
  L5_squared = (xx3 - xx1) * (xx3 - xx1) + (yy3 - yy1) * (yy3 - yy1)
      + (zz3 - zz1) * (zz3 - zz1);
  L6_squared = (xx3 - xx2) * (xx3 - xx2) + (yy3 - yy2) * (yy3 - yy2)
      + (zz3 - zz2) * (zz3 - zz2);

  L1 = sqrt(L1_squared);
  L4 = sqrt(L4_squared);
  L2 = sqrt(L2_squared);
  L5 = sqrt(L5_squared);
  L3 = sqrt(L3_squared);
  L6 = sqrt(L6_squared);

  jacobian = dComputeTetJacobian(adCoord0, adCoord1, adCoord2, adCoord3);

  adResult[0] = jacobian / (L1 * L2 * L3);
  adResult[1] = jacobian / (L1 * L4 * L5);
  adResult[2] = jacobian / (L2 * L4 * L6);
  adResult[3] = jacobian / (L3 * L5 * L6);

  *piNResult = 4;
}

void
vSMRSVolumeRatio(const double adCoord0[3], const double adCoord1[3],
		 const double adCoord2[3], const double adCoord3[3],
		 double adResult[6], int* const piNResult)
{
  double jacobian, volume, srms;
  double xx0, xx1, xx2, xx3;
  double yy0, yy1, yy2, yy3;
  double zz0, zz1, zz2, zz3;
  double L1_squared, L2_squared, L3_squared;
  double L4_squared, L5_squared, L6_squared;

  /*   computes the ratio of (1/6 sum_i L_i^2) ^ (3/2) / (8.47*Volume) */
  /*   this needs to be negated as it ranges from 1 to infinity (equilateral */
  /*   triangle is 1, zero volume is infinity) and we therefore need to */
  /*   minimize the maximum ratio */

  xx0 = adCoord0[XDIR];
  yy0 = adCoord0[YDIR];
  zz0 = adCoord0[ZDIR];
  xx1 = adCoord1[XDIR];
  yy1 = adCoord1[YDIR];
  zz1 = adCoord1[ZDIR];
  xx2 = adCoord2[XDIR];
  yy2 = adCoord2[YDIR];
  zz2 = adCoord2[ZDIR];
  xx3 = adCoord3[XDIR];
  yy3 = adCoord3[YDIR];
  zz3 = adCoord3[ZDIR];

  jacobian = dComputeTetJacobian(adCoord0, adCoord1, adCoord2, adCoord3);
  volume = fabs(jacobian) / 6.;

  L1_squared = (xx1 - xx0) * (xx1 - xx0) + (yy1 - yy0) * (yy1 - yy0)
      + (zz1 - zz0) * (zz1 - zz0);
  L2_squared = (xx2 - xx0) * (xx2 - xx0) + (yy2 - yy0) * (yy2 - yy0)
      + (zz2 - zz0) * (zz2 - zz0);
  L3_squared = (xx3 - xx0) * (xx3 - xx0) + (yy3 - yy0) * (yy3 - yy0)
      + (zz3 - zz0) * (zz3 - zz0);
  L4_squared = (xx1 - xx2) * (xx1 - xx2) + (yy1 - yy2) * (yy1 - yy2)
      + (zz1 - zz2) * (zz1 - zz2);
  L5_squared = (xx3 - xx1) * (xx3 - xx1) + (yy3 - yy1) * (yy3 - yy1)
      + (zz3 - zz1) * (zz3 - zz1);
  L6_squared = (xx3 - xx2) * (xx3 - xx2) + (yy3 - yy2) * (yy3 - yy2)
      + (zz3 - zz2) * (zz3 - zz2);

  srms = sqrt(
      1. / 6.
	  * (L1_squared + L2_squared + L3_squared + L4_squared + L5_squared
	      + L6_squared));

  adResult[0] = -pow(srms, 3.) / (8.47967 * volume);
  *piNResult = 1;
}

void
vComputeTetVolume(const double adCoord0[3], const double adCoord1[3],
		  const double adCoord2[3], const double adCoord3[3],
		  double adResult[6], int* const piNResult)
{
  double jacobian, volume;

  jacobian = dComputeTetJacobian(adCoord0, adCoord1, adCoord2, adCoord3);
  volume = jacobian / 6.;
  adResult[0] = volume;
  *piNResult = 1;
}

double
dComputeTetJacobian(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3])
{
  double a, b, c, d;
  double xx0, xx1, xx2, xx3;
  double yy0, yy1, yy2, yy3;
  double zz0, zz1, zz2, zz3;

  xx0 = adCoord0[XDIR];
  yy0 = adCoord0[YDIR];
  zz0 = adCoord0[ZDIR];
  xx1 = adCoord1[XDIR];
  yy1 = adCoord1[YDIR];
  zz1 = adCoord1[ZDIR];
  xx2 = adCoord2[XDIR];
  yy2 = adCoord2[YDIR];
  zz2 = adCoord2[ZDIR];
  xx3 = adCoord3[XDIR];
  yy3 = adCoord3[YDIR];
  zz3 = adCoord3[ZDIR];

  a = xx1 * (yy2 * zz3 - zz2 * yy3) - xx2 * (yy1 * zz3 - yy3 * zz1)
      + xx3 * (yy1 * zz2 - zz1 * yy2);
  b = (yy2 * zz3 - zz2 * yy3) - (yy1 * zz3 - yy3 * zz1)
      + (yy1 * zz2 - yy2 * zz1);
  c = xx1 * (zz3 - zz2) - xx2 * (zz3 - zz1) + xx3 * (zz2 - zz1);
  d = xx1 * (yy2 - yy3) - xx2 * (yy1 - yy3) + xx3 * (yy1 - yy2);

  return (a - b * xx0 - c * yy0 - d * zz0);
}

/*******************************************************************/
/*                  GRADIENT FUNCTIONS                             */
/*******************************************************************/

void
vGradSineDihedrals(const double adCoord0[3], const double adCoord1[3],
		   const double adCoord2[3], const double adCoord3[3],
		   double (*adGradient)[3], int* const piNGradient)
{
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);

  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);

  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */

  /* call the function */
  g_ad_vSineDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);

  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    register double *res = DERIV_GRAD(DTResult[i11]);
    adGradient[i11][0] = res[0];
    adGradient[i11][1] = res[1];
    adGradient[i11][2] = res[2];
  }
}

void
vGradCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
		     const double adCoord2[3], const double adCoord3[3],
		     double (*adGradient)[3], int* const piNGradient)
{
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);

  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);

  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */

  /* call the function */
  g_ad_vCosineDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);

  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    adGradient[i11][0] = DERIV_GRAD(DTResult[i11])[0];
    adGradient[i11][1] = DERIV_GRAD(DTResult[i11])[1];
    adGradient[i11][2] = DERIV_GRAD(DTResult[i11])[2];
  }
}

void
vGradSolids(const double adCoord0[3], const double adCoord1[3],
	    const double adCoord2[3], const double adCoord3[3],
	    double (*adGradient)[3], int* const piNGradient)
{
  /* First, compute the gradients of the dihedral angles. */
  double adGradDihed[6][3];
  int iNGrad;
  vGradDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adGradDihed, &iNGrad);
  assert(iNGrad == 6);

  /* Gradient of solid angle at 0 */
  adGradient[0][0] = adGradDihed[0][0] + adGradDihed[1][0] + adGradDihed[2][0];
  adGradient[0][1] = adGradDihed[0][1] + adGradDihed[1][1] + adGradDihed[2][1];
  adGradient[0][2] = adGradDihed[0][2] + adGradDihed[1][2] + adGradDihed[2][2];

  /* Gradient of solid angle at 1 */
  adGradient[1][0] = adGradDihed[0][0] + adGradDihed[3][0] + adGradDihed[4][0];
  adGradient[1][1] = adGradDihed[0][1] + adGradDihed[3][1] + adGradDihed[4][1];
  adGradient[1][2] = adGradDihed[0][2] + adGradDihed[3][2] + adGradDihed[4][2];

  /* Gradient of solid angle at 2 */
  adGradient[2][0] = adGradDihed[1][0] + adGradDihed[3][0] + adGradDihed[5][0];
  adGradient[2][1] = adGradDihed[1][1] + adGradDihed[3][1] + adGradDihed[5][1];
  adGradient[2][2] = adGradDihed[1][2] + adGradDihed[3][2] + adGradDihed[5][2];

  /* Gradient of solid angle at 3 */
  adGradient[3][0] = adGradDihed[2][0] + adGradDihed[4][0] + adGradDihed[5][0];
  adGradient[3][1] = adGradDihed[2][1] + adGradDihed[4][1] + adGradDihed[5][1];
  adGradient[3][2] = adGradDihed[2][2] + adGradDihed[4][2] + adGradDihed[5][2];

  *piNGradient = 4;
}

void
vGradDihedrals(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double (*adGradient)[3], int* const piNGradient)
{
  int i11;
  DERIV_TYPE DTVecA[3], DTVecB[3], DTVecC[3], DTVecD[3];
  DERIV_TYPE DTResult[6];

  ADALLOCATE(3);

  /* Put indep. values into the DERIV_TYPE structures. */
  COPY_VAL_TO_DERIV_ARRAY(DTVecA, adCoord0, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecB, adCoord1, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecC, adCoord2, 3);
  COPY_VAL_TO_DERIV_ARRAY(DTVecD, adCoord3, 3);

  /* Begin with all zeros for partials */
  ADZERO_ARRAY(DTVecA, 3);
  ADZERO_ARRAY(DTVecB, 3);
  ADZERO_ARRAY(DTVecC, 3);
  ADZERO_ARRAY(DTVecD, 3);

  /* Set the appropriate partials. */
  ADINIT_ARRAY(1.0, 1, DTVecA, 3); /* Indep #1-3 */

  /* call the function */
  g_ad_vDihedrals(DTVecA, DTVecB, DTVecC, DTVecD, DTResult, piNGradient);

  /* Unpack the gradients and return them */
  for (i11 = 0; i11 < 6; i11++) {
    adGradient[i11][0] = DERIV_GRAD(DTResult[i11])[0];
    adGradient[i11][1] = DERIV_GRAD(DTResult[i11])[1];
    adGradient[i11][2] = DERIV_GRAD(DTResult[i11])[2];
  }
}

void
vNegateGradDihedrals(const double adCoord0[3], const double adCoord1[3],
		     const double adCoord2[3], const double adCoord3[3],
		     double (*adGradient)[3], int* const piNGradient)
{
  int i, j;
  vGradDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adGradient,
		 piNGradient);
  for (i = 0; i < *piNGradient; i++) {
    for (j = 0; j < 3; j++)
      adGradient[i][j] = -adGradient[i][j];
  }
}

void
vNegateGradCosineDihedrals(const double adCoord0[3], const double adCoord1[3],
			   const double adCoord2[3], const double adCoord3[3],
			   double (*adGradient)[3], int* const piNGradient)
{
  int i, j;
  vGradCosineDihedrals(adCoord0, adCoord1, adCoord2, adCoord3, adGradient,
		       piNGradient);
  for (i = 0; i < *piNGradient; i++) {
    for (j = 0; j < 3; j++)
      adGradient[i][j] = -adGradient[i][j];
  }
}

void
vGradScaledJacobian(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
		    double (*adGradient)[3], int* const piNGradient)
{
  double xx0, xx1, xx2, xx3;
  double yy0, yy1, yy2, yy3;
  double zz0, zz1, zz2, zz3;
  double L1, L2, L3, L4, L5, L6;
  double L1_squared, L2_squared, L3_squared;
  double L4_squared, L5_squared, L6_squared;
  double L1_3_2, L2_3_2, L3_3_2;
  double a, b, c, d, J;
  double factor1, factor2, factor3, factor4;

  xx0 = adCoord0[XDIR];
  yy0 = adCoord0[YDIR];
  zz0 = adCoord0[ZDIR];
  xx1 = adCoord1[XDIR];
  yy1 = adCoord1[YDIR];
  zz1 = adCoord1[ZDIR];
  xx2 = adCoord2[XDIR];
  yy2 = adCoord2[YDIR];
  zz2 = adCoord2[ZDIR];
  xx3 = adCoord3[XDIR];
  yy3 = adCoord3[YDIR];
  zz3 = adCoord3[ZDIR];

  L1_squared = (xx1 - xx0) * (xx1 - xx0) + (yy1 - yy0) * (yy1 - yy0)
      + (zz1 - zz0) * (zz1 - zz0);
  L2_squared = (xx2 - xx0) * (xx2 - xx0) + (yy2 - yy0) * (yy2 - yy0)
      + (zz2 - zz0) * (zz2 - zz0);
  L3_squared = (xx3 - xx0) * (xx3 - xx0) + (yy3 - yy0) * (yy3 - yy0)
      + (zz3 - zz0) * (zz3 - zz0);
  L4_squared = (xx1 - xx2) * (xx1 - xx2) + (yy1 - yy2) * (yy1 - yy2)
      + (zz1 - zz2) * (zz1 - zz2);
  L5_squared = (xx3 - xx1) * (xx3 - xx1) + (yy3 - yy1) * (yy3 - yy1)
      + (zz3 - zz1) * (zz3 - zz1);
  L6_squared = (xx3 - xx2) * (xx3 - xx2) + (yy3 - yy2) * (yy3 - yy2)
      + (zz3 - zz2) * (zz3 - zz2);

  L1 = sqrt(L1_squared);
  L4 = sqrt(L4_squared);
  L2 = sqrt(L2_squared);
  L5 = sqrt(L5_squared);
  L3 = sqrt(L3_squared);
  L6 = sqrt(L6_squared);

  L1_3_2 = L1 * L1_squared;
  L2_3_2 = L2 * L2_squared;
  L3_3_2 = L3 * L3_squared;

  a = xx1 * (yy2 * zz3 - zz2 * yy3) - xx2 * (yy1 * zz3 - yy3 * zz1)
      + xx3 * (yy1 * zz2 - zz1 * yy2);
  b = (yy2 * zz3 - zz2 * yy3) - (yy1 * zz3 - yy3 * zz1)
      + (yy1 * zz2 - yy2 * zz1);
  c = xx1 * (zz3 - zz2) - xx2 * (zz3 - zz1) + xx3 * (zz2 - zz1);
  d = xx1 * (yy2 - yy3) - xx2 * (yy1 - yy3) + xx3 * (yy1 - yy2);

  J = a - b * xx0 - c * yy0 - d * zz0;

  factor1 = 1 / (L1 * L2 * L3);
  factor2 = J / (L1_3_2 * L2 * L3);
  factor3 = J / (L1 * L2_3_2 * L3);
  factor4 = J / (L1 * L2 * L3_3_2);

  adGradient[0][0] = -b * factor1 - (xx0 - xx1) * factor2
      - (xx0 - xx2) * factor3 - (xx0 - xx3) * factor4;
  adGradient[0][1] = -c * factor1 - (yy0 - yy1) * factor2
      - (yy0 - yy2) * factor3 - (yy0 - yy3) * factor4;
  adGradient[0][2] = -d * factor1 - (zz0 - zz1) * factor2
      - (zz0 - zz2) * factor3 - (zz0 - zz3) * factor4;

  factor1 = 1 / (L1 * L4 * L5);
  adGradient[1][0] = -b * factor1 - J * (xx0 - xx1) / (L1_3_2 * L4 * L5);
  adGradient[1][1] = -c * factor1 - J * (yy0 - yy1) / (L1_3_2 * L4 * L5);
  adGradient[1][2] = -d * factor1 - J * (zz0 - zz1) / (L1_3_2 * L4 * L5);

  factor1 = 1 / (L2 * L4 * L6);
  adGradient[2][0] = -b * factor1 - J * (xx0 - xx1) / (L2_3_2 * L4 * L6);
  adGradient[2][1] = -c * factor1 - J * (yy0 - yy1) / (L2_3_2 * L4 * L6);
  adGradient[2][2] = -d * factor1 - J * (zz0 - zz1) / (L2_3_2 * L4 * L6);

  factor1 = 1 / (L3 * L5 * L6);
  adGradient[3][0] = -b * factor1 - J * (xx0 - xx1) / (L3_3_2 * L5 * L6);
  adGradient[3][1] = -c * factor1 - J * (yy0 - yy1) / (L3_3_2 * L5 * L6);
  adGradient[3][2] = -d * factor1 - J * (zz0 - zz1) / (L3_3_2 * L5 * L6);

  *piNGradient = 4;
}

void
vGradSMRSVolumeRatio(const double adCoord0[3], const double adCoord1[3],
		     const double adCoord2[3], const double adCoord3[3],
		     double (*adGradient)[3], int* const piNGradient)
{
  /* double srms; */
  double xx0, yy0, zz0;
  /* double xx1, xx2, xx3; */
  /* double yy1, yy2, yy3; */
  /* double zz1, zz2, zz3; */
  /* double L1_squared, L2_squared, L3_squared; */
  /* double L4_squared, L5_squared, L6_squared; */
  /* double a, b, c, d; */

  xx0 = adCoord0[XDIR];
  yy0 = adCoord0[YDIR];
  zz0 = adCoord0[ZDIR];
  /* xx1 = adCoord1[XDIR];  yy1 = adCoord1[YDIR];  zz1 = adCoord1[ZDIR]; */
  /* xx2 = adCoord2[XDIR];  yy2 = adCoord2[YDIR];  zz2 = adCoord2[ZDIR]; */
  /* xx3 = adCoord3[XDIR];  yy3 = adCoord3[YDIR];  zz3 = adCoord3[ZDIR]; */

  /* a = xx1*(yy2*zz3-zz2*yy3) - xx2*(yy1*zz3-yy3*zz1) + xx3*(yy1*zz2 - zz1*yy2); */
  /* b = (yy2*zz3-zz2*yy3) - (yy1*zz3-yy3*zz1) + (yy1*zz2-yy2*zz1); */
  /* c = xx1*(zz3-zz2) - xx2*(zz3-zz1)+xx3*(zz2-zz1); */
  /* d = xx1*(yy2-yy3) - xx2*(yy1-yy3)+xx3*(yy1-yy2); */

  //  jacobian = - a + b*xx0 + c*yy0 + d*zz0;
  /* L1_squared = (xx1-xx0)*(xx1-xx0) + (yy1-yy0)*(yy1-yy0) + (zz1-zz0)*(zz1-zz0); */
  /* L2_squared = (xx2-xx0)*(xx2-xx0) + (yy2-yy0)*(yy2-yy0) + (zz2-zz0)*(zz2-zz0); */
  /* L3_squared = (xx3-xx0)*(xx3-xx0) + (yy3-yy0)*(yy3-yy0) + (zz3-zz0)*(zz3-zz0); */
  /* L4_squared = (xx1-xx2)*(xx1-xx2) + (yy1-yy2)*(yy1-yy2) + (zz1-zz2)*(zz1-zz2); */
  /* L5_squared = (xx3-xx1)*(xx3-xx1) + (yy3-yy1)*(yy3-yy1) + (zz3-zz1)*(zz3-zz1); */
  /* L6_squared = (xx3-xx2)*(xx3-xx2) + (yy3-yy2)*(yy3-yy2) + (zz3-zz2)*(zz3-zz2); */

  /* srms = sqrt(1./6.*(L1_squared+L2_squared +L3_squared+ */
  /* 		     L4_squared+L5_squared+L6_squared)); */
  //  srms3 = pow(srms,3.);
  {

    /*     double A, B, C, D; */
    /*     A = yy1*zz2 + zz1*yy3 + yy2*zz3 - yy3*zz2 - zz3*yy1 - yy2*zz1; */
    /*     B = xx1*zz2 + zz1*xx3 + xx2*zz3 - xx3*zz2 - zz3*xx1 - xx2*zz1; */
    /*     C = xx1*yy2 + yy1*xx3 + xx2*yy3 - xx3*yy2 - yy3*xx1 - xx2*yy1; */
    /*     D = xx1*yy2*zz3 + yy1*zz2*xx3 + zz1*xx2*yy3 - xx3*yy2*zz1 - yy3*zz2*xx1 - zz3*xx2*yy1; */

    /*     double E = A * xx0 - B * yy0 + C * zz0 - D; //this is 6 * volume. */
    /*     if(E < 0.) { A *= -1.; B *= -1; C *= -1.; E *= -1; } */

    /*     double sqrt_sum_sq_len  = sqrt(L1_squared+L2_squared +L3_squared+ */
    /* 				  L4_squared+L5_squared+L6_squared); */
    /*     double sqrt_sum_sq_len3 = sqrt_sum_sq_len * sqrt_sum_sq_len * sqrt_sum_sq_len; */

    /*     double scale = - sqrt(3) / 36; */

    /*     double grad = scale * (((A * sqrt_sum_sq_len3) / (E * E)) + ((1.5 * sqrt_sum_sq_len * (6. * xx0 - 2. * (xx1+xx2+xx3))) / E)); */

    double coordshi[] =
      { xx0 + 1.e-12, yy0, zz0 };
    double coordslo[] =
      { xx0 - 1.e-12, yy0, zz0 };

    double result[6];
    int num_result;
    double res1, res2;

    vSMRSVolumeRatio(coordshi, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res1 = result[0];
    vSMRSVolumeRatio(coordslo, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res2 = result[0];

    double gradx = (res1 - res2) / 2.e-12;

    coordshi[0] = xx0;
    coordshi[1] += 1.e-12;
    coordslo[0] = xx0;
    coordslo[1] -= 1.e-12;
    vSMRSVolumeRatio(coordshi, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res1 = result[0];
    vSMRSVolumeRatio(coordslo, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res2 = result[0];

    double grady = (res1 - res2) / 2.e-12;

    coordshi[1] = yy0;
    coordshi[2] += 1.e-12;
    coordslo[1] = yy0;
    coordslo[2] -= 1.e-12;
    vSMRSVolumeRatio(coordshi, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res1 = result[0];
    vSMRSVolumeRatio(coordslo, adCoord1, adCoord2, adCoord3, result,
		     &num_result);
    res2 = result[0];

    double gradz = (res1 - res2) / 2.e-12;

    /*     printf("finite diff grad = %e, analytic grad = %e\n", gradx, grad); */

    /*     abort(); */

    adGradient[0][0] = gradx;
    adGradient[0][1] = grady;
    adGradient[0][2] = gradz;

  }

  /*   adGradient[0][0] = -( 3./2.*srms*(-2*xx1+6*xx0-2*xx2 -2*xx3)/ (8.47967 * jacobian) */
  /* 			+ 6*srms3*b/(8.47967*jacobian*jacobian) ); */
  /*   adGradient[0][1] = -(3./2.*srms*(-2*yy1+6*yy0-2*yy2 -2*y3)/ (8.47967 * jacobian) */
  /* 		       + 6*srms3*c/(8.47967*jacobian*jacobian) ); */
  /*   adGradient[0][2] = -(3./2.*srms*(-2*zz1+6*zz0-2*z2 -2*zz3)/ (8.47967 * jacobian) */
  /* 		       + 6*srms3*d/(8.47967*jacobian*jacobian) ); */
  *piNGradient = 1;
}

void
vNormJacSquared(const double adCoord0[3], const double adCoord1[3],
		const double adCoord2[3], const double adCoord3[3],
		double adResult[6], int* const piNResult)
{
  double a1[3], a2[3], a3[3];

  /* the a matrix contains the Jacobian */
  a1[0] = adCoord1[0] - adCoord0[0];
  a1[1] = adCoord1[1] - adCoord0[1];
  a1[2] = adCoord1[2] - adCoord0[2];

  a2[0] = adCoord2[0] - adCoord0[0];
  a2[1] = adCoord2[1] - adCoord0[1];
  a2[2] = adCoord2[2] - adCoord0[2];

  a3[0] = adCoord3[0] - adCoord0[0];
  a3[1] = adCoord3[1] - adCoord0[1];
  a3[2] = adCoord3[2] - adCoord0[2];

  adResult[0] = SMfrobenius_norm_squared3x3(a1, a2, a3);
  adResult[0] *= -1.;

  *piNResult = 1;
}

void
vGradNormJacSquared(const double adCoord0[3], const double adCoord1[3],
		    const double adCoord2[3], const double adCoord3[3],
		    double (*adGradient)[3], int* const piNGradient)
{
  double a1[3], a2[3], a3[3];
  /* the a matrix contains the Jacobian */
  a1[0] = adCoord1[0] - adCoord0[0];
  a1[1] = adCoord1[1] - adCoord0[1];
  a1[2] = adCoord1[2] - adCoord0[2];

  a2[0] = adCoord2[0] - adCoord0[0];
  a2[1] = adCoord2[1] - adCoord0[1];
  a2[2] = adCoord2[2] - adCoord0[2];

  a3[0] = adCoord3[0] - adCoord0[0];
  a3[1] = adCoord3[1] - adCoord0[1];
  a3[2] = adCoord3[2] - adCoord0[2];

  adGradient[0][0] = -2 * (a1[0] + a2[0] + a3[0]);
  adGradient[0][1] = -2 * (a1[1] + a2[1] + a3[1]);
  adGradient[0][2] = -2 * (a1[2] + a2[2] + a3[2]);
  *piNGradient = 1;
}

void
vCondition(const double adCoord0[3], const double adCoord1[3],
	   const double adCoord2[3], const double adCoord3[3],
	   double adResult[6], int* const piNResult)
{
  double norm_a2, norm_adj2;
  double a1[3], a2[3], a3[3];
  double t1[3], t2[3], t3[3];
  double alpha;
  double f1, f1_hat;

  /* the a matrix contains the Jacobian */
  a1[0] = adCoord1[0] - adCoord0[0];
  a1[1] = adCoord1[1] - adCoord0[1];
  a1[2] = adCoord1[2] - adCoord0[2];

  a2[0] = adCoord2[0] - adCoord0[0];
  a2[1] = adCoord2[1] - adCoord0[1];
  a2[2] = adCoord2[2] - adCoord0[2];

  a3[0] = adCoord3[0] - adCoord0[0];
  a3[1] = adCoord3[1] - adCoord0[1];
  a3[2] = adCoord3[2] - adCoord0[2];

  alpha = SMdeterminant3x3(a1, a2, a3);
  norm_a2 = SMfrobenius_norm_squared3x3(a1, a2, a3);
  SMadjoint3x3(a1, a2, a3, t1, t2, t3);
  norm_adj2 = SMfrobenius_norm_squared3x3(t1, t2, t3);
  /*    norm_adj2 = SMfrobenius_norm_squared_adj3x3(a1,a2,a3); */

  f1 = norm_a2 / pow(fabs(alpha), 0.6666667);
  f1_hat = norm_adj2 / pow(fabs(alpha), 1.333333);

  adResult[0] = -sqrt(f1 * f1_hat);

  *piNResult = 1;
}

void
vGradCondition(const double adCoord0[3], const double adCoord1[3],
	       const double adCoord2[3], const double adCoord3[3],
	       double (*adGradient)[3], int* const piNGradient)
{
  double norma, norm_adja;
  double kappa, tmp, alpha;
  double column_1[3], column_2[3], column_3[3];
  double a1[3], a2[3], a3[3];
  double r1[3], r2[3], r3[3];
  double q1[3], q2[3], q3[3];
  double s1[3], s2[3], s3[3];
  double t1[3], t2[3], t3[3];

  /* the a matrix contains the Jacobian */
  a1[0] = adCoord1[0] - adCoord0[0];
  a1[1] = adCoord1[1] - adCoord0[1];
  a1[2] = adCoord1[2] - adCoord0[2];

  a2[0] = adCoord2[0] - adCoord0[0];
  a2[1] = adCoord2[1] - adCoord0[1];
  a2[2] = adCoord2[2] - adCoord0[2];

  a3[0] = adCoord3[0] - adCoord0[0];
  a3[1] = adCoord3[1] - adCoord0[1];
  a3[2] = adCoord3[2] - adCoord0[2];

  alpha = SMdeterminant3x3(a1, a2, a3);
  norma = SMfrobenius_norm_squared3x3(a1, a2, a3);
  SMadjoint3x3(a1, a2, a3, t1, t2, t3);
  norm_adja = SMfrobenius_norm_squared3x3(t1, t2, t3);
  /*   norm_adja = SMfrobenius_norm_squared_adj3x3( a1, a2, a3 ); */

  SMtranspose3x3(a1, a2, a3, q1, q2, q3);
  SMmultiply3x3(q1, q2, q3, a1, a2, a3, r1, r2, r3);
  SMmultiply3x3(a1, a2, a3, r1, r2, r3, s1, s2, s3);

  s1[0] *= norma;
  s1[1] *= norma;
  s1[2] *= norma;
  s2[0] *= norma;
  s2[1] *= norma;
  s2[2] *= norma;
  s3[0] *= norma;
  s3[1] *= norma;
  s3[2] *= norma;

  SMadjoint3x3(a1, a2, a3, q1, q2, q3);

  SMtranspose3x3(q1, q2, q3, r1, r2, r3);
  kappa = sqrt(norma * norm_adja) / alpha;

  tmp = alpha * kappa * kappa;
  s1[0] += tmp * r1[0];
  s1[1] += tmp * r1[1];
  s1[2] += tmp * r1[2];
  s2[0] += tmp * r2[0];
  s2[1] += tmp * r2[1];
  s2[2] += tmp * r2[2];
  s3[0] += tmp * r3[0];
  s3[1] += tmp * r3[1];
  s3[2] += tmp * r3[2];

  tmp = norma * norma + norm_adja;
  column_1[0] = tmp * a1[0] - s1[0];
  column_1[1] = tmp * a1[1] - s1[1];
  column_1[2] = tmp * a1[2] - s1[2];
  column_2[0] = tmp * a2[0] - s2[0];
  column_2[1] = tmp * a2[1] - s2[1];
  column_2[2] = tmp * a2[2] - s2[2];
  column_3[0] = tmp * a3[0] - s3[0];
  column_3[1] = tmp * a3[1] - s3[1];
  column_3[2] = tmp * a3[2] - s3[2];

  tmp = kappa * alpha * alpha;
  column_1[0] /= tmp;
  column_1[1] /= tmp;
  column_1[2] /= tmp;
  column_2[0] /= tmp;
  column_2[1] /= tmp;
  column_2[2] /= tmp;
  column_3[0] /= tmp;
  column_3[1] /= tmp;
  column_3[2] /= tmp;

  adGradient[0][0] = (column_1[0] + column_2[0] + column_3[0]);
  adGradient[0][1] = (column_1[1] + column_2[1] + column_3[1]);
  adGradient[0][2] = (column_1[2] + column_2[2] + column_3[2]);
  *piNGradient = 1;

  /* compare the analytic solution to finite differences */
  /*      { */
  /*        double function[6], function_pert[6]; */
  /*        double diff_grad[3]; */
  /*        double epsilon; */
  /*        int piNResult; */
  /*        double free_vtx[3]; */

  /*        free_vtx[0] = adCoord0[0]; */
  /*        free_vtx[1] = adCoord0[1]; */
  /*        free_vtx[2] = adCoord0[2]; */

  /*        epsilon = 1E-15; */
  /*        vCondition(adCoord0,adCoord1,adCoord2,adCoord3,function,&piNResult); */

  /*        free_vtx[0]+=epsilon; */
  /*        vCondition(free_vtx,adCoord1,adCoord2,adCoord3,function_pert,&piNResult); */
  /*        diff_grad[0] = (function_pert[0]-function[0])/epsilon; */

  /*        free_vtx[0]-=epsilon; */
  /*        free_vtx[1]+=epsilon; */
  /*        vCondition(free_vtx,adCoord1,adCoord2,adCoord3,function_pert,&piNResult); */
  /*        diff_grad[1] = (function_pert[0]-function[0])/epsilon; */

  /*        free_vtx[1]-=epsilon; */
  /*        free_vtx[2]+=epsilon; */
  /*        vCondition(free_vtx,adCoord1,adCoord2,adCoord3,function_pert,&piNResult); */
  /*        diff_grad[2] = (function_pert[0]-function[0])/epsilon; */

  /*        printf("Condition Gradient X: anal %f diff %f difference %f\n", */
  /*  	     adGradient[0][0], diff_grad[0], adGradient[0][0]-diff_grad[0]);  */
  /*        printf("Condition Gradient Y: anal %f diff %f difference %f\n", */
  /*  	     adGradient[0][1], diff_grad[1], adGradient[0][1]-diff_grad[1]); */
  /*        printf("Condition Gradient Z: anal %f diff %f difference %f\n", */
  /*  	     adGradient[0][2], diff_grad[2], adGradient[0][2]-diff_grad[2]); */
  /*      } */

}

