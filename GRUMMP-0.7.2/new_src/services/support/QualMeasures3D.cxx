#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_QualMeasure.h"
#include "GR_Quality.h"
#include "GR_Geometry.h"
#include "GR_TMOPQual.h"

#ifdef HAVE_MESQUITE
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "TShapeB1.hpp"
#include "IdealShapeTarget.hpp"
#include "TQualityMetric.hpp"
#endif

using namespace GRUMMP;

#include <math.h>

#ifdef HAVE_MESQUITE
/*
 * The next set of functions are used for getting an eigendecomposition of a metric for the Aniso adapt.
 */

// Symmetric Householder reduction to tridiagonal form.
void
tred(double V[3][3], double d[3], double e[3])
{

  int n = 3;
  for (int j = 0; j < n; j++) {
    d[j] = V[n - 1][j];
  }

  // Householder reduction to tridiagonal form.
  for (int i = n - 1; i > 0; i--) {

    // Scale to avoid under/overflow.
    double scale = 0.0;
    double h = 0.0;
    for (int k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i - 1];
      for (int j = 0; j < i; j++) {
	d[j] = V[i - 1][j];
	V[i][j] = 0.0;
	V[j][i] = 0.0;
      }
    }
    else {

      // Generate Householder vector.
      for (int k = 0; k < i; k++) {
	d[k] /= scale;
	h += d[k] * d[k];
      }
      double f = d[i - 1];
      double g = sqrt(h);
      if (f > 0) {
	g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i - 1] = f - g;
      for (int j = 0; j < i; j++) {
	e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.
      for (int j = 0; j < i; j++) {
	f = d[j];
	V[j][i] = f;
	g = e[j] + V[j][j] * f;
	for (int k = j + 1; k <= i - 1; k++) {
	  g += V[k][j] * d[k];
	  e[k] += V[k][j] * f;
	}
	e[j] = g;
      }
      f = 0.0;
      for (int j = 0; j < i; j++) {
	e[j] /= h;
	f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (int j = 0; j < i; j++) {
	e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
	f = d[j];
	g = e[j];
	for (int k = j; k <= i - 1; k++) {
	  V[k][j] -= (f * e[k] + g * d[k]);
	}
	d[j] = V[i - 1][j];
	V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.
  for (int i = 0; i < n - 1; i++) {
    V[n - 1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i + 1];
    if (h != 0.0) {
      for (int k = 0; k <= i; k++) {
	d[k] = V[k][i + 1] / h;
      }
      for (int j = 0; j <= i; j++) {
	double g = 0.0;
	for (int k = 0; k <= i; k++) {
	  g += V[k][i + 1] * V[k][j];
	}
	for (int k = 0; k <= i; k++) {
	  V[k][j] -= g * d[k];
	}
      }
    }
    for (int k = 0; k <= i; k++) {
      V[k][i + 1] = 0.0;
    }
  }
  for (int j = 0; j < n; j++) {
    d[j] = V[n - 1][j];
    V[n - 1][j] = 0.0;
  }
  V[n - 1][n - 1] = 1.0;
  e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.
void
tql(double V[3][3], double d[3], double e[3])
{

  int n = 3;
  for (int i = 1; i < n; i++) {
    e[i - 1] = e[i];
  }
  e[n - 1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0, -52.0);
  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element
    tst1 = std::max(tst1, fabs(d[l]) + fabs(e[l]));
    int m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps * tst1) {
	break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.
    if (m > l) {
      int iter = 0;
      do {
	iter = iter + 1;  // (Could check iteration count here.)

	// Compute implicit shift
	double g = d[l];
	double p = (d[l + 1] - g) / (2.0 * e[l]);

	double r = hypot(p, 1.0);

	if (p < 0) {
	  r = -r;
	}
	d[l] = e[l] / (p + r);
	d[l + 1] = e[l] * (p + r);
	double dl1 = d[l + 1];
	double h = g - d[l];
	for (int i = l + 2; i < n; i++) {
	  d[i] -= h;
	}
	f = f + h;

	// Implicit QL transformation.
	p = d[m];
	double c = 1.0;
	double c2 = c;
	double c3 = c;
	double el1 = e[l + 1];
	double s = 0.0;
	double s2 = 0.0;
	for (int i = m - 1; i >= l; i--) {
	  c3 = c2;
	  c2 = c;
	  s2 = s;
	  g = c * e[i];
	  h = c * p;
	  r = hypot(p, e[i]);
	  e[i + 1] = s * r;
	  s = e[i] / r;
	  c = p / r;
	  p = c * d[i] - s * g;
	  d[i + 1] = h + s * (c * g + s * d[i]);

	  // Accumulate transformation.
	  for (int k = 0; k < n; k++) {
	    h = V[k][i + 1];
	    V[k][i + 1] = s * V[k][i] + c * h;
	    V[k][i] = c * V[k][i] - s * h;
	  }
	}
	p = -s * s2 * c3 * el1 * e[l] / dl1;
	e[l] = s * p;
	d[l] = c * p;

	// Check for convergence.
      }
      while (fabs(e[l]) > eps * tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  // Sort eigenvalues and corresponding vectors.
  for (int i = 0; i < n - 1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i + 1; j < n; j++) {
      if (d[j] < p) {
	k = j;
	p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
	p = V[j][i];
	V[j][i] = V[j][k];
	V[j][k] = p;
      }
    }
  }
}

void
EigenDecomposition(double A[3][3], double V[3][3], double d[3])
{
  int n = 3;
  double e[3] =
    { 0 };
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred(V, d, e);
  tql(V, d, e);
}

/*
 * Computes the coordianates of a cell in metric space based on edge length.
 * The first point is set at the origin, the next one along an axis. If the length of all the other
 * edges is know, there is a unique triangle and tetrahedron. NOTE: orientation is not considered in this function.
 */
static void
convertTetToMetricSpace(const Vert* pV0, double *adLoc0, const Vert* pV1,
			double *adLoc1, const Vert* pV2, double *adLoc2,
			const Vert* pV3, double *adLoc3)
{

  double adMetric[6], a, b, c, d, e, f, r, s, t, u, v;

  for (unsigned i = 0; i < 6; i++) {
    adMetric[i] = (pV0->getMetric(i) + pV1->getMetric(i) + pV2->getMetric(i)
	+ pV3->getMetric(i)) / 4.;
  }

  GRUMMP::TMOPQual TQ;

  a = TQ.MLength3D(pV1->getCoords(), pV2->getCoords(), adMetric);
  b = TQ.MLength3D(pV0->getCoords(), pV2->getCoords(), adMetric);
  c = TQ.MLength3D(pV0->getCoords(), pV1->getCoords(), adMetric);
  d = TQ.MLength3D(pV1->getCoords(), pV3->getCoords(), adMetric);
  e = TQ.MLength3D(pV2->getCoords(), pV3->getCoords(), adMetric);
  f = TQ.MLength3D(pV0->getCoords(), pV3->getCoords(), adMetric);

  if (c == 0)
    c = 1e-15;
  r = (pow(c, 2) + pow(b, 2) - pow(a, 2)) / (2 * c);
  s = sqrt(fabs(pow(b, 2) - pow(r, 2)));
  t = (pow(c, 2) + pow(f, 2) - pow(d, 2)) / (2 * c);

  if (s == 0)
    s = 1e-15;
  u = (pow(e, 2) - pow(d, 2) - pow(s, 2) - pow(r, 2) + 2 * t * r + pow(c, 2)
      - 2 * t * c) / (-2 * s);
  v = sqrt(fabs(pow(f, 2) - pow(u, 2) - pow(t, 2)));

  adLoc0[0] = 0;
  adLoc0[1] = 0;
  adLoc0[2] = 0;

  adLoc1[0] = c;
  adLoc1[1] = 0;
  adLoc1[2] = 0;

  adLoc2[0] = r;
  adLoc2[1] = s;
  adLoc2[2] = 0;

  adLoc3[0] = t;
  adLoc3[1] = u;
  adLoc3[2] = v;
}
#endif

void
SineDihedralAngles::evalAllFromVerts(int& nValues, double* values,
				     const Vert* vert0, const Vert* vert1,
				     const Vert* vert2, const Vert* vert3,
				     const Vert* vert4, const Vert* vert5,
				     const Vert* vert6, const Vert* vert7) const
{
  // For now, only tets.
  assert(vert3);
  assert(!vert4 && !vert5 && !vert6 && !vert7);
#ifdef HAVE_MESQUITE
  if (bEvalMetric == true) {
    double adMA[3];
    double adMB[3];
    double adMC[3];
    double adMD[3];
    convertTetToMetricSpace(vert0, adMA, vert1, adMB, vert2, adMC, vert3, adMD);
    evalAllFromCoords(nValues, values, adMA, adMB, adMC, adMD);
  }
  else {
#endif

    const double * const adA = vert0->getCoords();
    const double * const adB = vert1->getCoords();
    const double * const adC = vert2->getCoords();
    const double * const adD = vert3->getCoords();

    evalAllFromCoords(nValues, values, adA, adB, adC, adD);
#ifdef HAVE_MESQUITE
  }
#endif
}

void
SineDihedralAngles::evalAllFromCoords(int& nValues, double* values,
				      const double *adA, const double *adB,
				      const double *adC, const double *adD,
				      const double *adE, const double *adF,
				      const double *adG,
				      const double *adH) const
{
  if (!adE && !adF && !adG && !adH) {
    nValues = 6;
    double adNorm0[3], adNorm1[3], adNorm2[3], adNorm3[3];
    calcUnitNormal(adC, adB, adD, adNorm0);
    calcUnitNormal(adA, adC, adD, adNorm1);
    calcUnitNormal(adB, adA, adD, adNorm2);
    calcUnitNormal(adA, adB, adC, adNorm3);

    // The order in which the dihedrals are assigned matters for
    // computation of solid angles.  The way they're currently set up,
    // combining them as (0,1,2), (0,3,4), (1,3,5), (2,4,5) gives (in
    // order) solid angles at vertices A, B, C, D

    double cross[3];

    // Always positive, even if the sine is negative, so no good for untangling.
    vCROSS3D(adNorm0, adNorm1, cross);
    values[0] = dMAG3D(cross); // Edge CD

    vCROSS3D(adNorm0, adNorm2, cross);
    values[1] = dMAG3D(cross); // Edge BD

    vCROSS3D(adNorm0, adNorm3, cross);
    values[2] = dMAG3D(cross); // Edge BC

    vCROSS3D(adNorm1, adNorm2, cross);
    values[3] = dMAG3D(cross); // Edge AD

    vCROSS3D(adNorm1, adNorm3, cross);
    values[4] = dMAG3D(cross); // Edge AC

    vCROSS3D(adNorm2, adNorm3, cross);
    values[5] = dMAG3D(cross); // Edge AB
  }
  else if (!adF && !adG && !adH) {
    // Pyramid
    nValues = 0; // actually 8
  }
  else if (!adG && !adH) {
    // Prism
    nValues = 9;
    evalAllDihedsFromCoords(nValues, values, adA, adB, adC, adD, adE, adF);
    // It's probably possible to do this more efficiently, but
    // until it's a performance issue, who cares?
    for (int ii = 0; ii < 9; ii++) {
      values[ii] = sin(values[ii]);
    }
  }
  else {
    // Hex
    nValues = 0; // actually 12
  }
}

void
DihedralAngles::evalAllFromVerts(int& nValues, double* values,
				 const Vert* vert0, const Vert* vert1,
				 const Vert* vert2, const Vert* vert3,
				 const Vert* vert4, const Vert* vert5,
				 const Vert* vert6, const Vert* vert7) const
{
  // For now, only tets.
  assert(vert3);

#ifdef HAVE_MESQUITE
  if (bEvalMetric == true) {
    assert(!vert4 && !vert5 && !vert6 && !vert7);
    double adMA[3];
    double adMB[3];
    double adMC[3];
    double adMD[3];
    convertTetToMetricSpace(vert0, adMA, vert1, adMB, vert2, adMC, vert3, adMD);
    evalAllFromCoords(nValues, values, adMA, adMB, adMC, adMD);

  }
  else {
#endif

    const double * const ad0 = vert0->getCoords();
    const double * const ad1 = vert1->getCoords();
    const double * const ad2 = vert2->getCoords();
    const double * const ad3 = vert3->getCoords();
    const double *ad4 = nullptr, *ad5 = nullptr, *ad6 = nullptr, *ad7 = nullptr;
    if (vert4)
      ad4 = vert4->getCoords();
    if (vert5)
      ad5 = vert5->getCoords();
    if (vert6)
      ad6 = vert6->getCoords();
    if (vert7)
      ad7 = vert7->getCoords();
    evalAllFromCoords(nValues, values, ad0, ad1, ad2, ad3, ad4, ad5, ad6, ad7);
#ifdef HAVE_MESQUITE
  }
#endif
}

void
DihedralAngles::evalAllFromCoords(int& nValues, double* values,
				  const double *adA, const double *adB,
				  const double *adC, const double *adD,
				  const double *adE, const double *adF,
				  const double *adG, const double *adH) const
{
  assert(adA && adB && adC && adD);
  if (!adE && !adF && !adG && !adH) {
    nValues = 6;
    double adNorm0[3], adNorm1[3], adNorm2[3], adNorm3[3];
    calcUnitNormal(adC, adB, adD, adNorm0);
    calcUnitNormal(adA, adC, adD, adNorm1);
    calcUnitNormal(adB, adA, adD, adNorm2);
    calcUnitNormal(adA, adB, adC, adNorm3);

    // The order in which the dihedrals are assigned matters for
    // computation of solid angles.  The way they're currently set up,
    // combining them as (0,1,2), (0,3,4), (1,3,5), (2,4,5) gives (in
    // order) solid angles at vertices A, B, C, D

    double dot = -dDOT3D(adNorm0, adNorm1);
    values[5] = GR_acos(dot); // Edge CD

    dot = -dDOT3D(adNorm0, adNorm2);
    values[4] = GR_acos(dot); // Edge BD

    dot = -dDOT3D(adNorm0, adNorm3);
    values[3] = GR_acos(dot); // Edge BC

    dot = -dDOT3D(adNorm1, adNorm2);
    values[2] = GR_acos(dot); // Edge AD

    dot = -dDOT3D(adNorm1, adNorm3);
    values[1] = GR_acos(dot); // Edge AC

    dot = -dDOT3D(adNorm2, adNorm3);
    values[0] = GR_acos(dot); // Edge AB
  }
  else if (!adF && !adG && !adH) {
    // Pyramid; eight dihedrals
    nValues = 0; // For now
  }
  else if (!adG && !adH) {
    // Prism; nine dihedrals
    nValues = 9;
    evalAllDihedsFromCoords(nValues, values, adA, adB, adC, adD, adE, adF);
  }
  else {
    // Hex; twelve dihedrals
    nValues = 0; // For now
  }
}

void
SolidAngles::evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
			      const Vert* vert1, const Vert* vert2,
			      const Vert* vert3, const Vert* vert4,
			      const Vert* vert5, const Vert* vert6,
			      const Vert* vert7) const
{
  // For now, only tets.
  assert(vert3);
  assert(!vert4 && !vert5 && !vert6 && !vert7);

#ifdef HAVE_MESQUITE
  if (bEvalMetric == true) {
    double adMA[3];
    double adMB[3];
    double adMC[3];
    double adMD[3];
    convertTetToMetricSpace(vert0, adMA, vert1, adMB, vert2, adMC, vert3, adMD);
    evalAllFromCoords(nValues, values, adMA, adMB, adMC, adMD);

  }
  else {
#endif

    const double * const adA = vert0->getCoords();
    const double * const adB = vert1->getCoords();
    const double * const adC = vert2->getCoords();
    const double * const adD = vert3->getCoords();
    evalAllFromCoords(nValues, values, adA, adB, adC, adD);

#ifdef HAVE_MESQUITE
  }
#endif
}

void
SolidAngles::evalAllFromCoords(int& nValues, double* values, const double *adA,
			       const double *adB, const double *adC,
			       const double *adD, const double *adLoc4,
			       const double *adLoc5, const double *adLoc6,
			       const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  int nDiheds = -1;
  double diheds[12];
  m_DA.evalAllFromCoords(nDiheds, diheds, adA, adB, adC, adD, adLoc4, adLoc5,
			 adLoc6, adLoc7);
  assert(nDiheds == 6); // Only tets for now.
  nValues = 4;

  // The order in which the dihedrals are assigned matters for
  // computation of solid angles.  The way they're currently set up,
  // combining them as (0,1,2), (0,3,4), (1,3,5), (2,4,5) gives (in
  // order) solid angles at vertices A, B, C, D

  values[0] = diheds[0] + diheds[1] + diheds[2] - M_PI;
  values[1] = diheds[0] + diheds[3] + diheds[4] - M_PI;
  values[2] = diheds[1] + diheds[3] + diheds[5] - M_PI;
  values[3] = diheds[2] + diheds[4] + diheds[5] - M_PI;
}

///////////////////////////////////////////////////////////////////////////////
void
ShortestToRadius3D::evalAllFromVerts(int& nValues, double* values,
				     const Vert* vert0, const Vert* vert1,
				     const Vert* vert2, const Vert* vert3,
				     const Vert* vert4, const Vert* vert5,
				     const Vert* vert6, const Vert* vert7) const
{
  assert(!vert4 && !vert5 && !vert6 && !vert7);

#ifdef HAVE_MESQUITE
  if (bEvalMetric == true) {
    double adMA[3];
    double adMB[3];
    double adMC[3];
    double adMD[3];
    convertTetToMetricSpace(vert0, adMA, vert1, adMB, vert2, adMC, vert3, adMD);

    nValues = 1;
    int iN = 1;
    int* piN = &iN;
    calcShortestToRadius(adMA, adMB, adMC, adMD, piN, values);
  }
  else {
#endif

    nValues = 1;
    int iN = 1;
    int* piN = &iN;
    calcShortestToRadius(vert0->getCoords(), vert1->getCoords(),
			 vert2->getCoords(), vert3->getCoords(), piN, values);

#ifdef HAVE_MESQUITE
  }
#endif
}

void
ShortestToRadius3D::evalAllFromCoords(int& nValues, double* values,
				      const double *adA, const double *adB,
				      const double *adC, const double *adD,
				      const double *adLoc4,
				      const double *adLoc5,
				      const double *adLoc6,
				      const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  nValues = 1;
  int iN = 1;
  int* piN = &iN;
  calcShortestToRadius(adA, adB, adC, adD, piN, values);
}

void
RadiusRatio3D::evalAllFromCoords(int& nValues, double* values,
				 const double *adA, const double *adB,
				 const double *adC, const double *adD,
				 const double *adLoc4, const double *adLoc5,
				 const double *adLoc6,
				 const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  nValues = 1;
  // It would be really nice not to have to go through all this;
  // importing the calcs for this quality measure might be easier...
  Vert VA, VB, VC, VD;
  VA.setCoords(3, adA);
  VB.setCoords(3, adB);
  VC.setCoords(3, adC);
  VD.setCoords(3, adD);
  double circRad = dCircumradius(&VA, &VB, &VC, &VD);
  double inRad = dInradius(&VA, &VB, &VC, &VD);
//	logMessage(MSG_MANAGER, "Inradius: %f  Circradius: %f\n", inRad, circRad);

  values[0] = 3 * inRad / circRad;
}

void
MarcumSkewness3D::evalAllFromCoords(int& nValues, double* values,
				    const double *adA, const double *adB,
				    const double *adC, const double *adD,
				    const double *adLoc4, const double *adLoc5,
				    const double *adLoc6,
				    const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  nValues = 1;
  static constexpr double scaling = 9. / 8 * sqrt(3);
  // It would be really nice not to have to go through all this;
  // importing the calcs for this quality measure might be easier...
  Vert VA, VB, VC, VD;
  VA.setCoords(3, adA);
  VB.setCoords(3, adB);
  VC.setCoords(3, adC);
  VD.setCoords(3, adD);
  double vol = dVolTet(&VA, &VB, &VC, &VD);
  double radius = dCircumradius(&VA, &VB, &VC, &VD);
  logMessage(MSG_MANAGER, "Volume: %f  Radius: %f\n", vol, radius);
  values[0] = scaling * vol / (radius * radius * radius);
}

void
VolumeArea3D::evalAllFromCoords(int& nValues, double* values, const double *adA,
				const double *adB, const double *adC,
				const double *adD, const double *adLoc4,
				const double *adLoc5, const double *adLoc6,
				const double *adLoc7) const
{
  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  nValues = 1;
  static constexpr double scaling = 2 * M_SQRT2 * pow(3., 1.75); // about 20

  // It would be really nice not to have to go through all this;
  // importing the calcs for this quality measure might be easier...
  Vert VA, VB, VC, VD;
  VA.setCoords(3, adA);
  VB.setCoords(3, adB);
  VC.setCoords(3, adC);
  VD.setCoords(3, adD);
  double vol = dVolTet(&VA, &VB, &VC, &VD);

  double areaABC = dAreaTri(&VA, &VB, &VC);
  double areaABD = dAreaTri(&VA, &VB, &VD);
  double areaBCD = dAreaTri(&VB, &VC, &VD);
  double areaCAD = dAreaTri(&VC, &VD, &VA);

  double totalArea = areaABC + areaABD + areaBCD + areaCAD;

  values[0] = scaling * vol / pow(totalArea, 1.5);
}

void
VolumeLength3D::evalAllFromVerts(int& nValues, double* values,
				 const Vert* vert0, const Vert* vert1,
				 const Vert* vert2, const Vert* vert3,
				 const Vert* vert4, const Vert* vert5,
				 const Vert* vert6, const Vert* vert7) const
{
  assert(!vert4 && !vert5 && !vert6 && !vert7);
#ifdef HAVE_MESQUITE
  if (bEvalMetric == true) {

    double adMA[3], adMB[3], adMC[3], adMD[3], dRMS3D, dVol3D, adMetric[6];

    convertTetToMetricSpace(vert0, adMA, vert1, adMB, vert2, adMC, vert3, adMD);
    for (unsigned i = 0; i < 6; i++) {
      adMetric[i] = (vert0->getMetric(i) + vert1->getMetric(i)
	  + vert2->getMetric(i) + vert3->getMetric(i)) / 4.;
    }

    double a, b, c, d, e, f;

    GRUMMP::TMOPQual TQ;

    a = TQ.MLength3D(vert1->getCoords(), vert2->getCoords(), adMetric);
    b = TQ.MLength3D(vert0->getCoords(), vert2->getCoords(), adMetric);
    c = TQ.MLength3D(vert0->getCoords(), vert1->getCoords(), adMetric);
    d = TQ.MLength3D(vert1->getCoords(), vert3->getCoords(), adMetric);
    e = TQ.MLength3D(vert2->getCoords(), vert3->getCoords(), adMetric);
    f = TQ.MLength3D(vert0->getCoords(), vert3->getCoords(), adMetric);

    dRMS3D = pow(
	sqrt(
	    (pow(a, 2) + pow(b, 2) + pow(c, 2) + pow(d, 2) + pow(e, 2)
		+ pow(f, 2)) / 6.),
	3);

    nValues = 1;

    dVol3D = fabs(
	(1 / 6.)
	    * ((adMA[0] - adMD[0])
		* ((adMB[1] - adMD[1]) * (adMC[2] - adMD[2])
		    - (adMC[1] - adMD[1]) * (adMB[2] - adMD[2]))
		- (adMB[0] - adMD[0])
		    * ((adMA[1] - adMD[1]) * (adMC[2] - adMD[2])
			- (adMC[1] - adMD[1]) * (adMA[2] - adMD[2]))
		+ (adMC[0] - adMD[0])
		    * ((adMA[1] - adMD[1]) * (adMB[2] - adMD[2])
			- (adMB[1] - adMD[1]) * (adMA[2] - adMD[2]))));

    values[0] = dVol3D / dRMS3D;
    values[0] = values[0] * 6 * sqrt(2);
    values[0] = fabs(1 - values[0]);
  }
  else {
#endif
    const double * const adA = vert0->getCoords();
    const double * const adB = vert1->getCoords();
    const double * const adC = vert2->getCoords();
    const double * const adD = vert3->getCoords();
    evalAllFromCoords(nValues, values, adA, adB, adC, adD);
#ifdef HAVE_MESQUITE
  }
#endif
}

void
VolumeLength3D::evalAllFromCoords(int& nValues, double* values,
				  const double *adA, const double *adB,
				  const double *adC, const double *adD,
				  const double *adLoc4, const double *adLoc5,
				  const double *adLoc6,
				  const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  nValues = 1;

  double a = dDIST_SQ_3D(adB, adC);
  double b = dDIST_SQ_3D(adA, adC);
  double c = dDIST_SQ_3D(adA, adB);
  double d = dDIST_SQ_3D(adB, adD);
  double e = dDIST_SQ_3D(adC, adD);
  double f = dDIST_SQ_3D(adA, adD);

  double dRMS3D = sqrt((a + b + c + d + e + f) / 6.);

  double dVol3D = fabs(
      (1 / 6.)
	  * ((adA[0] - adD[0])
	      * ((adB[1] - adD[1]) * (adC[2] - adD[2])
		  - (adC[1] - adD[1]) * (adB[2] - adD[2]))
	      - (adB[0] - adD[0])
		  * ((adA[1] - adD[1]) * (adC[2] - adD[2])
		      - (adC[1] - adD[1]) * (adA[2] - adD[2]))
	      + (adC[0] - adD[0])
		  * ((adA[1] - adD[1]) * (adB[2] - adD[2])
		      - (adB[1] - adD[1]) * (adA[2] - adD[2]))));

  values[0] = dVol3D / (dRMS3D * dRMS3D * dRMS3D);
  values[0] = values[0] * 6 * sqrt(2);
  values[0] = fabs(1 - values[0]);

}

