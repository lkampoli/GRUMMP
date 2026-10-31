/*
 * TMOPQual.cxx
 *
 *  Created on: Jun 30, 2017
 *      Author: cfog
 */

#include "GR_config.h"
#include "GR_Cell.h"
#include "GR_TMOPQual.h"
#include "GR_Vec.h"

#ifdef HAVE_MESQUITE
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "IdealShapeTarget.hpp"
#include "TQualityMetric.hpp"

#include "TShapeB1.hpp"
#include "TShapeOrientB1.hpp"
#include "TShapeSizeB1.hpp"
#include "TSizeB1.hpp"
#include "TShapeSizeOrientB1.hpp"

#include "MeshImpl.hpp"
#include "ReferenceMesh.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "TQualityMetric.hpp"
#include "GR_Aniso.h"
#endif

////////////************/////////////
#ifdef HAVE_MESQUITE

using namespace Mesquite;

using GRUMMP::TMOPQual;

void
TMOPQual::CheckTypeCell(CellSkel::eCellType Type)
{

  switch (Type)
    {
    case CellSkel::eTriCell:
    case CellSkel::eTriCV:
      iVs = 3;
      iNM = 3;
      break;
    case CellSkel::eQuadCell:
    case CellSkel::eQuadCV:
      iVs = 4;
      iNM = 3;

      break;
    case CellSkel::eTet:
    case CellSkel::eTetCV:
      iVs = 4;
      iNM = 6;
      break;
    case CellSkel::ePyr:
    case CellSkel::ePyrCV:
      iVs = 5;
      iNM = 6;
      break;
    case CellSkel::ePrism:
    case CellSkel::ePrismCV:
      iVs = 6;
      iNM = 6;
      break;
    case CellSkel::eHex:
    case CellSkel::eHexCV:
      iVs = 8;
      iNM = 6;
      break;
    default:
      assert(0);
      break;
    }
  iT = iVs * 3;
}

bool
TMOPQual::CheckMetric(const Vert *apV[])
{
  bool bM = true;
  for (int i = 0; i < iVs; i++) {
    assert(apV[i]->isValid());
    for (int ii = 0; ii < iNM; ii++)
      if (!apV[i]->getMetric(ii))
	bM = false;
  }
  return bM;
}

bool
TMOPQual::isMetricAvailable(CellSkel::eCellType Type, const Vert* pV0,
			    const Vert* pV1, const Vert*pV2, const Vert* pV3,
			    const Vert* pV4, const Vert* pV5, const Vert* pV6,
			    const Vert* pV7)
{
  bool bM;
  const Vert *apV[8] =
    { pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7 };
  CheckTypeCell(Type);
  bM = CheckMetric(apV);
  return bM;
}

Mesquite::TMetric*
TMOPQual::SingleQualMetric(enum eQMetric eQM)
{
  switch (eQM)
    {
    case (eShape):
      return new Mesquite::TShapeB1;
      break;
    case (eSize):
      return new Mesquite::TSizeB1;
      break;
    case (eShapeSize):
      return new Mesquite::TShapeSizeB1;
      break;
    case (eShapeOrient):
      return new Mesquite::TShapeOrientB1;
      break;
    case (eShapeSizeOrient):
      return new Mesquite::TShapeSizeOrientB1;
      break;
    default:
      return new Mesquite::TShapeB1;
      break;
    }
  return NULL;
}

void
TMOPQual::EigenDecomposition2D(double *dA, double *dVec, double *dEV)
{

  double dRoot = sqrt((dA[0] - dA[2]) * (dA[0] - dA[2]) + 4 * dA[1] * dA[1]);
  dEV[0] = 0.5 * ((dA[0] + dA[2]) + dRoot);                  //Larger eigenvalue
  dEV[1] = 0.5 * ((dA[0] + dA[2]) - dRoot);                 //smaller eigenvalue

  if (dA[1] < 1e-6 && dA[1] > -1e-6) {

    dVec[0] = 1;
    dVec[1] = 0;
    dVec[2] = 0;
    dVec[3] = 1;
    double adRes[2], adDiff[2];
    adRes[0] = dA[0] * dVec[0] + dA[1] * dVec[1];

    adRes[1] = dA[1] * dVec[0] + dA[2] * dVec[1];

    adDiff[0] = fabs(adRes[0] - dEV[0] * dVec[0]);
    adDiff[1] = fabs(adRes[1] - dEV[0] * dVec[1]);
    if (adDiff[0] > 1e-6 || adDiff[1] > 1e-6) {
      dVec[0] = 0;
      dVec[1] = 1;
      dVec[2] = 1;
      dVec[3] = 0;
    }
  }
  else {

    dVec[2] = dA[0] - dEV[0];
    dVec[0] = -dA[1];

    double dMag = sqrt(dVec[0] * dVec[0] + dVec[2] * dVec[2]);
    dVec[2] /= dMag;
    dVec[0] /= dMag;

    //EigenVector corresponding to smaller eigenvalue
    dVec[3] = dA[0] - dEV[1];
    dVec[1] = -dA[1];
    dMag = sqrt(dVec[1] * dVec[1] + dVec[3] * dVec[3]);
    dVec[3] /= dMag;
    dVec[1] /= dMag;
  }
}

void
TMOPQual::tred(double V[3][3], double d[3], double e[3])
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
TMOPQual::tql(double V[3][3], double d[3], double e[3])
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
TMOPQual::EigenDecomposition3D(double *dA, double *dVec, double *d)
{

  double A[3][3];
  A[0][0] = dA[0];
  A[0][1] = dA[1];
  A[0][2] = dA[3];
  A[1][0] = dA[1];
  A[1][1] = dA[2];
  A[1][2] = dA[4];
  A[2][0] = dA[3];
  A[2][1] = dA[4];
  A[2][2] = dA[5];
  double V[3][3];
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

  double dNorm0 = sqrt(
      V[0][0] * V[0][0] + V[1][0] * V[1][0] + V[2][0] * V[2][0]);
  double dNorm1 = sqrt(
      V[0][1] * V[0][1] + V[1][1] * V[1][1] + V[2][1] * V[2][1]);
  double dNorm2 = sqrt(
      V[0][2] * V[0][2] + V[1][2] * V[1][2] + V[2][2] * V[2][2]);

  dVec[0] = V[0][0] / dNorm0;
  dVec[1] = V[1][0] / dNorm0;
  dVec[2] = V[2][0] / dNorm0;

  dVec[3] = V[0][1] / dNorm1;
  dVec[4] = V[1][1] / dNorm1;
  dVec[5] = V[2][1] / dNorm1;

  dVec[6] = V[0][2] / dNorm2;
  dVec[7] = V[1][2] / dNorm2;
  dVec[8] = V[2][2] / dNorm2;
}

void
TMOPQual::getEigenDecomposition(double *dA, double *dVec, double *dEV)
{

  if (iNM == 3)
    EigenDecomposition2D(dA, dVec, dEV);
  else if (iNM == 6) {
    EigenDecomposition3D(dA, dVec, dEV);
  }
  else
    assert(0);
}

void
TMOPQual::SetTargetfromEigen(CellSkel::eCellType Type, double *dT, double *dVec,
			     double *dEV)
{

  double dTest[iNM == 3 ? 4 : 9];
  if (iNM == 3) {

    dTest[0] = dVec[0] * (1 / sqrt(dEV[0]));
    dTest[1] = dVec[1] * (1 / sqrt(dEV[0]));

    dTest[2] = dVec[2] * (1 / sqrt(dEV[1]));
    dTest[3] = dVec[3] * (1 / sqrt(dEV[1]));

    double loc0[2] =
      { 0, 0 };
    double loc1[2] =
      { dTest[0], dTest[1] };
    double loc2[2] =
      { dTest[2], dTest[3] };

    switch (Type)
      {
      case CellSkel::eTriCell:
      case CellSkel::eTriCV:
	if (checkOrient2D(loc0, loc1, loc2) == 1) {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = 0;
	  dT[3] = loc1[0];
	  dT[4] = loc1[1];
	  dT[5] = 0;
	  dT[6] = loc2[0];
	  dT[7] = loc2[1];
	  dT[8] = 0;
	}
	else {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = 0;
	  dT[3] = loc2[0];
	  dT[4] = loc2[1];
	  dT[5] = 0;
	  dT[6] = loc1[0];
	  dT[7] = loc1[1];
	  dT[8] = 0;
	}
	break;
      case CellSkel::eQuadCell:
      case CellSkel::eQuadCV:
	if (checkOrient2D(loc0, loc1, loc2) == 1) {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = 0;
	  dT[3] = loc1[0];
	  dT[4] = loc1[1];
	  dT[5] = 0;
	  dT[6] = loc1[0] + loc2[0];
	  dT[7] = loc1[1] + loc2[1];
	  dT[8] = 0;
	  dT[9] = loc2[0];
	  dT[10] = loc2[1];
	  dT[11] = 0;
	}
	else {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = 0;
	  dT[3] = loc2[0];
	  dT[4] = loc2[1];
	  dT[5] = 0;
	  dT[6] = loc1[0] + loc2[0];
	  dT[7] = loc1[1] + loc2[1];
	  dT[8] = 0;
	  dT[9] = loc1[0];
	  dT[10] = loc1[1];
	  dT[11] = 0;
	}
	break;
      default:
	assert(0);
	break;
      }
  }
  if (iNM == 6) {

    dTest[0] = dVec[0] * (1 / sqrt(dEV[0]));
    dTest[1] = dVec[1] * (1 / sqrt(dEV[0]));
    dTest[2] = dVec[2] * (1 / sqrt(dEV[0]));

    dTest[3] = dVec[3] * (1 / sqrt(dEV[1]));
    dTest[4] = dVec[4] * (1 / sqrt(dEV[1]));
    dTest[5] = dVec[5] * (1 / sqrt(dEV[1]));

    dTest[6] = dVec[6] * (1 / sqrt(dEV[2]));
    dTest[7] = dVec[7] * (1 / sqrt(dEV[2]));
    dTest[8] = dVec[8] * (1 / sqrt(dEV[2]));

    double loc0[3] =
      { 0, 0, 0 };
    double loc1[3] =
      { dTest[0], dTest[1], dTest[2] };
    double loc2[3] =
      { dTest[3], dTest[4], dTest[5] };
    double loc3[3] =
      { dTest[6], dTest[7], dTest[8] };

    switch (Type)
      {
      case CellSkel::eTet:
      case CellSkel::eTetCV:
	if (checkOrient3D(loc0, loc1, loc2, loc3) == 1) {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = loc0[2];
	  dT[3] = loc1[0];
	  dT[4] = loc1[1];
	  dT[5] = loc1[2];
	  dT[6] = loc2[0];
	  dT[7] = loc2[1];
	  dT[8] = loc2[2];
	  dT[9] = loc3[0];
	  dT[10] = loc3[1];
	  dT[11] = loc3[2];
	}
	else {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = loc0[2];
	  dT[3] = loc2[0];
	  dT[4] = loc2[1];
	  dT[5] = loc2[2];
	  dT[6] = loc1[0];
	  dT[7] = loc1[1];
	  dT[8] = loc1[2];
	  dT[9] = loc3[0];
	  dT[10] = loc3[1];
	  dT[11] = loc3[2];
	}
	break;
      case CellSkel::ePyr:
      case CellSkel::ePyrCV:
	assert(0);
	break;
      case CellSkel::ePrism:
      case CellSkel::ePrismCV:
	assert(0);
	break;
      case CellSkel::eHex:
      case CellSkel::eHexCV:
	if (checkOrient3D(loc0, loc1, loc2, loc3) == 1) {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = loc0[2];
	  dT[3] = loc1[0];
	  dT[4] = loc1[1];
	  dT[5] = loc1[2];
	  dT[6] = loc1[0] + loc2[0];
	  dT[7] = loc1[1] + loc2[1];
	  dT[8] = loc1[2] + loc2[2];
	  dT[9] = loc2[0];
	  dT[10] = loc2[1];
	  dT[11] = loc2[2];
	  dT[12] = loc3[0];
	  dT[13] = loc3[1];
	  dT[14] = loc3[2];
	  dT[15] = loc3[0] + loc1[0];
	  dT[16] = loc3[1] + loc1[1];
	  dT[17] = loc3[2] + loc1[2];
	  dT[18] = loc3[0] + loc1[0] + loc2[0];
	  dT[19] = loc3[1] + loc1[1] + loc2[1];
	  dT[20] = loc3[2] + loc1[2] + loc2[2];
	  dT[21] = loc3[0] + loc2[0];
	  dT[22] = loc3[1] + loc2[1];
	  dT[23] = loc3[2] + loc2[2];
	}
	else {
	  dT[0] = loc0[0];
	  dT[1] = loc0[1];
	  dT[2] = loc0[2];
	  dT[3] = loc2[0];
	  dT[4] = loc2[1];
	  dT[5] = loc2[2];
	  dT[6] = loc1[0] + loc2[0];
	  dT[7] = loc1[1] + loc2[1];
	  dT[8] = loc1[2] + loc2[2];
	  dT[9] = loc1[0];
	  dT[10] = loc1[1];
	  dT[11] = loc1[2];
	  dT[12] = loc3[0];
	  dT[13] = loc3[1];
	  dT[14] = loc3[2];
	  dT[15] = loc3[0] + loc2[0];
	  dT[16] = loc3[1] + loc2[1];
	  dT[17] = loc3[2] + loc2[2];
	  dT[18] = loc3[0] + loc1[0] + loc2[0];
	  dT[19] = loc3[1] + loc1[1] + loc2[1];
	  dT[20] = loc3[2] + loc1[2] + loc2[2];
	  dT[21] = loc3[0] + loc1[0];
	  dT[22] = loc3[1] + loc1[1];
	  dT[23] = loc3[2] + loc1[2];
	}
	break;
      default:
	assert(0);
	break;
      }
  }
}

void
TMOPQual::SetTargetfromMetricNewVertice(double *adTarget,
					CellSkel::eCellType Type,
					double *dMetric, const Vert* pV0,
					const Vert* pV1, const Vert* pV2,
					const Vert* pV3, const Vert* pV4,
					const Vert* pV5, const Vert* pV6,
					const Vert* pV7)
{

  assert(!pV3 && !pV4 && !pV5 && !pV6 && !pV7);

  double adM[iNM];
  double adM0[iNM];
  double adM1[iNM];
  double aT[iT];
  double adEV[iNM == 3 ? 2 : 3];
  double adEVec[iNM == 3 ? 4 : 9];

  for (int i = 0; i < iNM; i++)
    adM[i] = 0;

  if (iNM == 3) {
    assert(Type == CellSkel::eTriCell);
    for (int i = 0; i < iNM; i++) {
      adM0[i] = pV0->getMetric(i);
      adM1[i] = pV1->getMetric(i);
    }
    dLogMetricCell(adM0, adM1, dMetric, adM);
  }

  if (iNM == 6) {
    assert(Type == CellSkel::eTet);
    for (int i = 0; i < iNM; i++) {
      adM[i] = (pV0->getMetric(i) + pV1->getMetric(i) + pV2->getMetric(i)
	  + dMetric[i]) / 4.;
    }
  }

  getEigenDecomposition(adM, adEVec, adEV);
  SetTargetfromEigen(Type, aT, adEVec, adEV);

  for (int i = 0; i < iT; i++)
    adTarget[i] = aT[i];
}

void
TMOPQual::SetTargetfromMetric(double *adTarget, CellSkel::eCellType Type,
			      const Vert* pV0, const Vert* pV1, const Vert* pV2,
			      const Vert* pV3, const Vert* pV4, const Vert* pV5,
			      const Vert* pV6, const Vert* pV7)
{

  isMetricAvailable(Type, pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7);

  double adM[iNM];
  double aT[iT];
  double adEV[iNM == 3 ? 2 : 3];
  double adEVec[iNM == 3 ? 4 : 9];

  for (int i = 0; i < iNM; i++)
    adM[i] = 0;

  if (iNM == 3)
    dLogMetricCell(pV0, pV1, pV2, adM);

  if (iNM == 6) {
    for (int i = 0; i < iNM; i++) {
      adM[i] = (pV0->getMetric(i) + pV1->getMetric(i) + pV2->getMetric(i)
	  + pV3->getMetric(i)) / 4.;
    }
  }

  getEigenDecomposition(adM, adEVec, adEV);
  SetTargetfromEigen(Type, aT, adEVec, adEV);

  for (int i = 0; i < iT; i++)
    adTarget[i] = aT[i];
}

double
TMOPQual::dEval(enum eQMetric eQM, CellSkel::eCellType Type, const Vert* pV0,
		const Vert* pV1, const Vert*pV2, const Vert* pV3,
		const Vert* pV4, const Vert* pV5, const Vert* pV6,
		const Vert* pV7)
{

  double dQ;
  switch (Type)
    {
    case CellSkel::eTriCell:
    case CellSkel::eTriCV:
      double dTriTarget[9];
      SetTargetfromMetric(dTriTarget, Type, pV0, pV1, pV2);
      dQ = dEval(eQM, Type, dTriTarget, pV0->getCoords(), pV1->getCoords(),
		 pV2->getCoords());
      break;
    case CellSkel::eQuadCell:
    case CellSkel::eQuadCV:
      double dQuadTarget[12];
      SetTargetfromMetric(dQuadTarget, Type, pV0, pV1, pV2, pV3);
      dQ = dEval(eQM, Type, dQuadTarget, pV0->getCoords(), pV1->getCoords(),
		 pV2->getCoords(), pV3->getCoords());
      break;
    case CellSkel::eTet:
    case CellSkel::eTetCV:
      double dTetTarget[12];
      SetTargetfromMetric(dTetTarget, Type, pV0, pV1, pV2, pV3);
      dQ = dEval(eQM, Type, dTetTarget, pV0->getCoords(), pV1->getCoords(),
		 pV2->getCoords(), pV3->getCoords());
      break;
    case CellSkel::ePyr:
    case CellSkel::ePyrCV:
      assert(0);
      break;
    case CellSkel::ePrism:
    case CellSkel::ePrismCV:
      assert(0);
      break;
    case CellSkel::eHex:
    case CellSkel::eHexCV:
      double dHexTarget[24];
      SetTargetfromMetric(dHexTarget, Type, pV0, pV1, pV2, pV3, pV4, pV5, pV6,
			  pV7);
      dQ = dEval(eQM, Type, dHexTarget, pV0->getCoords(), pV1->getCoords(),
		 pV2->getCoords(), pV3->getCoords(), pV4->getCoords(),
		 pV5->getCoords(), pV6->getCoords(), pV7->getCoords());
      break;
    default:
      assert(0);
      break;
    }
  return dQ;

}

double
TMOPQual::dEvalgetBestTarget(enum eQMetric eQM, CellSkel::eCellType Type,
			     double* dT, double *dTbest, const double *adLoc0,
			     const double *adLoc1, const double *adLoc2)
{

  assert(Type == CellSkel::eTriCell || Type == CellSkel::eTet);

  const int connect[3] =
    { 0, 1, 2 };
  const unsigned long int Connect[3] =
    { 0, 1, 2 };
  const bool Fixed[3] =
    { 0, 0, 0 };
  double dCell0[9];
  double dTarget1[9];
  double dTarget2[9];
  double dTarget3[9];
  double dTarget4[9];
  double dQual;
  double dMinQual = 1e30;

  dCell0[0] = adLoc0[0];
  dCell0[1] = adLoc0[1];
  dCell0[2] = 0;
  dCell0[3] = adLoc1[0];
  dCell0[4] = adLoc1[1];
  dCell0[5] = 0;
  dCell0[6] = adLoc2[0];
  dCell0[7] = adLoc2[1];
  dCell0[8] = 0;

  for (unsigned i = 0; i < 9; i++)
    dTarget1[i] = dT[i];

  Mesquite::TMetric* pM_QMetric;
  pM_QMetric = SingleQualMetric(eQM);

  Mesquite::MsqError Error;
  Mesquite::PatchData CellPatch;
  CellPatch.fill(3, dCell0, 1, Mesquite::TRIANGLE, Connect, Fixed, Error);

  Mesquite::MeshImpl Mesh(3, 1, Mesquite::TRIANGLE, Fixed, dTarget1, connect);
  Mesquite::ReferenceMesh Ref_Mesh(&Mesh);
  Mesquite::RefMeshTargetCalculator Tc(&Ref_Mesh, true);

  dTarget2[0] = -dTarget1[0];
  dTarget2[3] = -dTarget1[3];
  dTarget2[6] = -dTarget1[6];

  dTarget2[1] = -dTarget1[1];
  dTarget2[4] = -dTarget1[4];
  dTarget2[7] = -dTarget1[7];

  Mesquite::MeshImpl Mesh2(3, 1, Mesquite::TRIANGLE, Fixed, dTarget2, connect);
  Mesquite::ReferenceMesh Ref_Mesh2(&Mesh2);
  Mesquite::RefMeshTargetCalculator Tc2(&Ref_Mesh2, true);

//	dTarget3[0] = dT[0];
//	dTarget3[1] = dT[1];
//	dTarget3[2] = dT[2];
//	dTarget3[3] = dT[3];
//	dTarget3[4] = dT[4] + dT[7];
//	dTarget3[5] = dT[5];
//	dTarget3[6] = dT[6];
//	dTarget3[7] = dT[7];
//	dTarget3[8] = dT[8];

  dTarget3[0] = dT[0];
  dTarget3[1] = dT[1];
  dTarget3[2] = dT[2];
  dTarget3[3] = dT[6];
  dTarget3[4] = dT[7];
  dTarget3[5] = dT[8];
  dTarget3[6] = -dT[3];
  dTarget3[7] = -dT[4];
  dTarget3[8] = -dT[5];

  Mesquite::MeshImpl Mesh3(3, 1, Mesquite::TRIANGLE, Fixed, dTarget3, connect);
  Mesquite::ReferenceMesh Ref_Mesh3(&Mesh3);
  Mesquite::RefMeshTargetCalculator Tc3(&Ref_Mesh3, true);

//	dTarget4[0] = dT[0];
//	dTarget4[1] = dT[1];
//	dTarget4[2] = dT[2];
//	dTarget4[3] = -dT[3];
//	dTarget4[4] = -(dT[4] + dT[7]);
//	dTarget4[5] = -dT[5];
//	dTarget4[6] = -dT[6];
//	dTarget4[7] = -dT[7];
//	dTarget4[8] = -dT[8];

  dTarget4[0] = dT[0];
  dTarget4[1] = dT[1];
  dTarget4[2] = dT[2];
  dTarget4[3] = -dT[6];
  dTarget4[4] = -dT[7];
  dTarget4[5] = -dT[8];
  dTarget4[6] = dT[3];
  dTarget4[7] = dT[4];
  dTarget4[8] = dT[5];

  Mesquite::MeshImpl Mesh4(3, 1, Mesquite::TRIANGLE, Fixed, dTarget4, connect);
  Mesquite::ReferenceMesh Ref_Mesh4(&Mesh4);
  Mesquite::RefMeshTargetCalculator Tc4(&Ref_Mesh4, true);

  Mesquite::TQualityMetric TQMetric(&Tc, pM_QMetric);
  Mesquite::TQualityMetric TQMetric2(&Tc2, pM_QMetric);
  Mesquite::TQualityMetric TQMetric3(&Tc3, pM_QMetric);
  Mesquite::TQualityMetric TQMetric4(&Tc4, pM_QMetric);

  Mesquite::TQualityMetric *pTQ[4] =
    { &TQMetric, &TQMetric2, &TQMetric3, &TQMetric4 };

  pTQ[0]->evaluate(CellPatch, 0, dQual, Error);

  dMinQual = dQual;

  dTbest[3] = dTarget1[3];
  dTbest[6] = dTarget1[6];
  dTbest[4] = dTarget1[4];
  dTbest[7] = dTarget1[7];

  pTQ[1]->evaluate(CellPatch, 0, dQual, Error);

  if (dQual < dMinQual) {
    dTbest[3] = dTarget2[3];
    dTbest[6] = dTarget2[6];
    dTbest[4] = dTarget2[4];
    dTbest[7] = dTarget2[7];
    dMinQual = dQual;
  }

  pTQ[2]->evaluate(CellPatch, 0, dQual, Error);

  if (dQual < dMinQual) {
    dTbest[3] = dTarget3[3];
    dTbest[6] = dTarget3[6];
    dTbest[4] = dTarget3[4];
    dTbest[7] = dTarget3[7];
    dMinQual = dQual;
  }

  pTQ[3]->evaluate(CellPatch, 0, dQual, Error);

  if (dQual < dMinQual) {
    dTbest[3] = dTarget4[3];
    dTbest[6] = dTarget4[6];
    dTbest[4] = dTarget4[4];
    dTbest[7] = dTarget4[7];
    dMinQual = dQual;
  }

  return dMinQual;
}

double
TMOPQual::dEval(enum eQMetric eQM, CellSkel::eCellType Type, double* dT,
		const double *adLoc0, const double *adLoc1,
		const double *adLoc2, const double *adLoc3,
		const double *adLoc4, const double *adLoc5,
		const double *adLoc6, const double *adLoc7)
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  double dQ;
  switch (Type)
    {
    case CellSkel::eTriCell:
    case CellSkel::eTriCV:
      {
	const int connect[3] =
	  { 0, 1, 2 };
	const unsigned long int Connect[3] =
	  { 0, 1, 2 };
	const bool Fixed[3] =
	  { 0, 0, 0 };
	double dCell0[9];
	double dTarget[9];
	double dQual;
	double dMinQual = 1e30;

	dCell0[0] = adLoc0[0];
	dCell0[1] = adLoc0[1];
	dCell0[2] = 0;
	dCell0[3] = adLoc1[0];
	dCell0[4] = adLoc1[1];
	dCell0[5] = 0;
	dCell0[6] = adLoc2[0];
	dCell0[7] = adLoc2[1];
	dCell0[8] = 0;

	for (unsigned i = 0; i < 9; i++)
	  dTarget[i] = dT[i];

	Mesquite::TMetric* pM_QMetric;
	pM_QMetric = SingleQualMetric(eQM);

	Mesquite::MsqError Error;
	Mesquite::PatchData CellPatch;
	CellPatch.fill(3, dCell0, 1, Mesquite::TRIANGLE, Connect, Fixed, Error);

	if (false) {
	  Mesquite::IdealShapeTarget TcIso;
	  Mesquite::TQualityMetric TQMetricIso(&TcIso, pM_QMetric);
	  TQMetricIso.evaluate(CellPatch, 0, dQ, Error);
	  return dQ;
	}

	Mesquite::MeshImpl Mesh(3, 1, Mesquite::TRIANGLE, Fixed, dTarget,
				connect);
	Mesquite::ReferenceMesh Ref_Mesh(&Mesh);
	Mesquite::RefMeshTargetCalculator Tc(&Ref_Mesh, true);

	dTarget[0] = -dTarget[0];
	dTarget[3] = -dTarget[3];
	dTarget[6] = -dTarget[6];

	dTarget[1] = -dTarget[1];
	dTarget[4] = -dTarget[4];
	dTarget[7] = -dTarget[7];

	Mesquite::MeshImpl Mesh1(3, 1, Mesquite::TRIANGLE, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh1(&Mesh1);
	Mesquite::RefMeshTargetCalculator Tc1(&Ref_Mesh1, true);

	dTarget[0] = dT[3];
	dTarget[1] = dT[4];
	dTarget[2] = dT[5];
	dTarget[3] = dT[0];
	dTarget[4] = dT[1];
	dTarget[5] = dT[2];
	dTarget[6] = -dT[6];
	dTarget[7] = -dT[7];
	dTarget[8] = -dT[8];

	Mesquite::MeshImpl Mesh2(3, 1, Mesquite::TRIANGLE, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh2(&Mesh2);
	Mesquite::RefMeshTargetCalculator Tc2(&Ref_Mesh2, true);

	dTarget[0] = -dTarget[0];
	dTarget[3] = -dTarget[3];
	dTarget[6] = -dTarget[6];

	dTarget[1] = -dTarget[1];
	dTarget[4] = -dTarget[4];
	dTarget[7] = -dTarget[7];

	Mesquite::MeshImpl Mesh3(3, 1, Mesquite::TRIANGLE, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh3(&Mesh3);
	Mesquite::RefMeshTargetCalculator Tc3(&Ref_Mesh3, true);

	Mesquite::TQualityMetric TQMetric(&Tc, pM_QMetric);
	Mesquite::TQualityMetric TQMetric1(&Tc1, pM_QMetric);
	Mesquite::TQualityMetric TQMetric2(&Tc2, pM_QMetric);
	Mesquite::TQualityMetric TQMetric3(&Tc3, pM_QMetric);

	Mesquite::PatchData CellPatch1;
	Mesquite::PatchData CellPatch2;

	double dCell1[9] =
	  { dCell0[3], dCell0[4], dCell0[5], dCell0[6], dCell0[7], dCell0[8],
	      dCell0[0], dCell0[1], dCell0[2] };

	CellPatch1.fill(3, dCell1, 1, Mesquite::TRIANGLE, Connect, Fixed,
			Error);

	double dCell2[9] =
	  { dCell0[6], dCell0[7], dCell0[8], dCell0[0], dCell0[1], dCell0[2],
	      dCell0[3], dCell0[4], dCell0[5] };

	CellPatch2.fill(3, dCell2, 1, Mesquite::TRIANGLE, Connect, Fixed,
			Error);

	Mesquite::TQualityMetric *pTQ[4] =
	  { &TQMetric, &TQMetric1, &TQMetric2, &TQMetric3 };
	Mesquite::PatchData *pCP[3] =
	  { &CellPatch, &CellPatch1, &CellPatch2 };

	for (int i = 0; i < 3; i++) {
	  for (int ii = 0; ii < 4; ii++) {
	    pTQ[ii]->evaluate(*pCP[i], 0, dQual, Error);
	    dMinQual = dQual < dMinQual ? dQual : dMinQual;
	  }
	}

	dQ = dMinQual;
	delete pM_QMetric;
	return dQ;
      }
      break;
    case CellSkel::eQuadCell:
    case CellSkel::eQuadCV:
      {
	const int connect[4] =
	  { 0, 1, 2, 3 };
	const unsigned long int Connect[4] =
	  { 0, 1, 2, 3 };
	const bool Fixed[4] =
	  { 0, 0, 0, 0 };
	double dCell0[12];
	double dTarget[12];
	double dQual;
	double dMinQual = 1e30;

	dCell0[0] = adLoc0[0];
	dCell0[1] = adLoc0[1];
	dCell0[2] = 0;
	dCell0[3] = adLoc1[0];
	dCell0[4] = adLoc1[1];
	dCell0[5] = 0;
	dCell0[6] = adLoc2[0];
	dCell0[7] = adLoc2[1];
	dCell0[8] = 0;
	dCell0[9] = adLoc3[0];
	dCell0[10] = adLoc3[1];
	dCell0[11] = 0;

	for (unsigned i = 0; i < 12; i++)
	  dTarget[i] = dT[i];

	Mesquite::TMetric* pM_QMetric;
	pM_QMetric = SingleQualMetric(eQM);

	Mesquite::MsqError Error;
	Mesquite::PatchData CellPatch;
	CellPatch.fill(4, dCell0, 1, Mesquite::QUADRILATERAL, Connect, Fixed,
		       Error);

	if (false) {
	  Mesquite::IdealShapeTarget TcIso;
	  Mesquite::TQualityMetric TQMetricIso(&TcIso, pM_QMetric);
	  TQMetricIso.evaluate(CellPatch, 0, dQ, Error);
	  return dQ;
	}

	Mesquite::MeshImpl Mesh(4, 1, Mesquite::QUADRILATERAL, Fixed, dTarget,
				connect);
	Mesquite::ReferenceMesh Ref_Mesh(&Mesh);
	Mesquite::RefMeshTargetCalculator Tc(&Ref_Mesh, true);

	dTarget[0] = -dTarget[0];
	dTarget[3] = -dTarget[3];
	dTarget[6] = -dTarget[6];
	dTarget[9] = -dTarget[9];

	dTarget[1] = -dTarget[1];
	dTarget[4] = -dTarget[4];
	dTarget[7] = -dTarget[7];
	dTarget[10] = -dTarget[10];

	Mesquite::MeshImpl Mesh2(4, 1, Mesquite::QUADRILATERAL, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh2(&Mesh2);
	Mesquite::RefMeshTargetCalculator Tc2(&Ref_Mesh2, true);

	dTarget[0] = dT[3];
	dTarget[1] = dT[4];
	dTarget[2] = dT[5];

	dTarget[3] = dT[0];
	dTarget[4] = dT[1];
	dTarget[5] = dT[2];

	dTarget[6] = -dT[9];
	dTarget[7] = -dT[10];
	dTarget[8] = -dT[11];

	dTarget[9] = -dT[6];
	dTarget[10] = -dT[7];
	dTarget[11] = -dT[8];

	Mesquite::MeshImpl Mesh1(4, 1, Mesquite::QUADRILATERAL, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh1(&Mesh1);
	Mesquite::RefMeshTargetCalculator Tc1(&Ref_Mesh1, true);

	dTarget[0] = -dTarget[0];
	dTarget[3] = -dTarget[3];
	dTarget[6] = -dTarget[6];
	dTarget[9] = -dTarget[9];

	dTarget[1] = -dTarget[1];
	dTarget[4] = -dTarget[4];
	dTarget[7] = -dTarget[7];
	dTarget[10] = -dTarget[10];

	Mesquite::MeshImpl Mesh3(4, 1, Mesquite::QUADRILATERAL, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh3(&Mesh3);
	Mesquite::RefMeshTargetCalculator Tc3(&Ref_Mesh3, true);

	Mesquite::TShapeB1 Tshape;

	Mesquite::TQualityMetric TQMetric(&Tc, pM_QMetric);
	Mesquite::TQualityMetric TQMetric1(&Tc1, pM_QMetric);
	Mesquite::TQualityMetric TQMetric2(&Tc2, pM_QMetric);
	Mesquite::TQualityMetric TQMetric3(&Tc3, pM_QMetric);

	Mesquite::PatchData CellPatch1;
	Mesquite::PatchData CellPatch2;
	Mesquite::PatchData CellPatch3;

	double dCell1[12] =
	      { dCell0[3], dCell0[4], dCell0[5], dCell0[6], dCell0[7],
		  dCell0[8], dCell0[9], dCell0[10], dCell0[11], dCell0[0],
		  dCell0[1], dCell0[2] };

	CellPatch1.fill(4, dCell1, 1, Mesquite::QUADRILATERAL, Connect, Fixed,
			Error);

	double dCell2[12] =
	  { dCell0[6], dCell0[7], dCell0[8], dCell0[9], dCell0[10], dCell0[11],
	      dCell0[0], dCell0[1], dCell0[2], dCell0[3], dCell0[4], dCell0[5] };

	CellPatch2.fill(4, dCell2, 1, Mesquite::QUADRILATERAL, Connect, Fixed,
			Error);

	double dCell3[12] =
	  { dCell0[9], dCell0[10], dCell0[11], dCell0[0], dCell0[1], dCell0[2],
	      dCell0[3], dCell0[4], dCell0[5], dCell0[6], dCell0[7], dCell0[8] };

	CellPatch3.fill(4, dCell3, 1, Mesquite::QUADRILATERAL, Connect, Fixed,
			Error);

	Mesquite::TQualityMetric *pTQ[4] =
	  { &TQMetric, &TQMetric1, &TQMetric2, &TQMetric3 };
	Mesquite::PatchData *pCP[4] =
	  { &CellPatch, &CellPatch1, &CellPatch2, &CellPatch3 };

	for (int i = 0; i < 4; i++) {
	  for (int ii = 0; ii < 4; ii++) {
	    pTQ[ii]->evaluate(*pCP[i], 0, dQual, Error);
	    dMinQual = dQual < dMinQual ? dQual : dMinQual;
	  }
	}

	dQ = dMinQual;
	delete pM_QMetric;

	return dQ;
      }
      break;
    case CellSkel::eTet:
    case CellSkel::eTetCV:
      {
	const int connect[4] =
	  { 0, 1, 2, 3 };
	const unsigned long int Connect[4] =
	  { 0, 1, 2, 3 };
	const bool Fixed[4] =
	  { 0, 0, 0, 0 };
	double dCell0[12];
	double dTarget[12];
	double dQual;
	double dMinQual = 1e30;

	dCell0[0] = adLoc0[0];
	dCell0[1] = adLoc0[1];
	dCell0[2] = adLoc0[2];
	dCell0[3] = adLoc1[0];
	dCell0[4] = adLoc1[1];
	dCell0[5] = adLoc1[2];
	dCell0[6] = adLoc2[0];
	dCell0[7] = adLoc2[1];
	dCell0[8] = adLoc2[2];
	dCell0[9] = adLoc3[0];
	dCell0[10] = adLoc3[1];
	dCell0[11] = adLoc3[2];

	for (unsigned i = 0; i < 12; i++)
	  dTarget[i] = dT[i];

	Mesquite::TMetric* pM_QMetric;
	pM_QMetric = SingleQualMetric(eQM);

	Mesquite::MsqError Error;
	Mesquite::PatchData CellPatch;
	CellPatch.fill(4, dCell0, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
		       Error);

	if (false) {
	  Mesquite::IdealShapeTarget TcIso;
	  Mesquite::TQualityMetric TQMetricIso(&TcIso, pM_QMetric);
	  TQMetricIso.evaluate(CellPatch, 0, dQ, Error);
	  return dQ;
	}

	Mesquite::MeshImpl Mesh(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				connect);
	Mesquite::ReferenceMesh Ref_Mesh(&Mesh);
	Mesquite::RefMeshTargetCalculator Tc(&Ref_Mesh, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = -dT[6];
	dTarget[4] = -dT[7];
	dTarget[5] = -dT[8];
	dTarget[6] = dT[3];
	dTarget[7] = dT[4];
	dTarget[8] = dT[5];
	dTarget[9] = dT[9];
	dTarget[10] = dT[10];
	dTarget[11] = dT[11];

	Mesquite::MeshImpl Mesh1(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh1(&Mesh1);
	Mesquite::RefMeshTargetCalculator Tc1(&Ref_Mesh1, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = -dT[3];
	dTarget[4] = -dT[4];
	dTarget[5] = -dT[5];
	dTarget[6] = -dT[6];
	dTarget[7] = -dT[7];
	dTarget[8] = -dT[8];
	dTarget[9] = dT[9];
	dTarget[10] = dT[10];
	dTarget[11] = dT[11];

	Mesquite::MeshImpl Mesh2(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh2(&Mesh2);
	Mesquite::RefMeshTargetCalculator Tc2(&Ref_Mesh2, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = dT[6];
	dTarget[4] = dT[7];
	dTarget[5] = dT[8];
	dTarget[6] = -dT[3];
	dTarget[7] = -dT[4];
	dTarget[8] = -dT[5];
	dTarget[9] = dT[9];
	dTarget[10] = dT[10];
	dTarget[11] = dT[11];

	Mesquite::MeshImpl Mesh3(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh3(&Mesh3);
	Mesquite::RefMeshTargetCalculator Tc3(&Ref_Mesh3, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = dT[6];
	dTarget[4] = dT[7];
	dTarget[5] = dT[8];
	dTarget[6] = dT[3];
	dTarget[7] = dT[4];
	dTarget[8] = dT[5];
	dTarget[9] = -dT[9];
	dTarget[10] = -dT[10];
	dTarget[11] = -dT[11];

	Mesquite::MeshImpl Mesh4(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh4(&Mesh4);
	Mesquite::RefMeshTargetCalculator Tc4(&Ref_Mesh4, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = dT[3];
	dTarget[4] = dT[4];
	dTarget[5] = dT[5];
	dTarget[6] = -dT[6];
	dTarget[7] = -dT[7];
	dTarget[8] = -dT[8];
	dTarget[9] = -dT[9];
	dTarget[10] = -dT[10];
	dTarget[11] = -dT[11];

	Mesquite::MeshImpl Mesh5(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh5(&Mesh5);
	Mesquite::RefMeshTargetCalculator Tc5(&Ref_Mesh5, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = -dT[6];
	dTarget[4] = -dT[7];
	dTarget[5] = -dT[8];
	dTarget[6] = -dT[3];
	dTarget[7] = -dT[4];
	dTarget[8] = -dT[5];
	dTarget[9] = -dT[9];
	dTarget[10] = -dT[10];
	dTarget[11] = -dT[11];

	Mesquite::MeshImpl Mesh6(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh6(&Mesh6);
	Mesquite::RefMeshTargetCalculator Tc6(&Ref_Mesh6, true);

	dTarget[0] = dT[0];
	dTarget[1] = dT[1];
	dTarget[2] = dT[2];
	dTarget[3] = -dT[3];
	dTarget[4] = -dT[4];
	dTarget[5] = -dT[5];
	dTarget[6] = dT[6];
	dTarget[7] = dT[7];
	dTarget[8] = dT[8];
	dTarget[9] = -dT[9];
	dTarget[10] = -dT[10];
	dTarget[11] = -dT[11];

	Mesquite::MeshImpl Mesh7(4, 1, Mesquite::TETRAHEDRON, Fixed, dTarget,
				 connect);
	Mesquite::ReferenceMesh Ref_Mesh7(&Mesh7);
	Mesquite::RefMeshTargetCalculator Tc7(&Ref_Mesh7, true);

	Mesquite::TShapeB1 Tshape;

	Mesquite::TQualityMetric TQMetric(&Tc, pM_QMetric), TQMetric1(
	    &Tc1, pM_QMetric), TQMetric2(&Tc2, pM_QMetric), TQMetric3(
	    &Tc3, pM_QMetric), TQMetric4(&Tc4, pM_QMetric), TQMetric5(
	    &Tc5, pM_QMetric), TQMetric6(&Tc6, pM_QMetric), TQMetric7(
	    &Tc7, pM_QMetric);

	Mesquite::PatchData CellPatch1, CellPatch2, CellPatch3, CellPatch4,
	    CellPatch5, CellPatch6, CellPatch7, CellPatch8, CellPatch9,
	    CellPatch10, CellPatch11;

	double dCell1[12] =
	  { dCell0[0], dCell0[1], dCell0[2], dCell0[9], dCell0[10], dCell0[11],
	      dCell0[3], dCell0[4], dCell0[5], dCell0[6], dCell0[7], dCell0[8] };

	CellPatch1.fill(4, dCell1, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell2[12] =
	      { dCell0[3], dCell0[4], dCell0[5], dCell0[0], dCell0[1],
		  dCell0[2], dCell0[9], dCell0[10], dCell0[11], dCell0[6],
		  dCell0[7], dCell0[8] };

	CellPatch2.fill(4, dCell2, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell3[12] =
	  { dCell0[9], dCell0[10], dCell0[11], dCell0[3], dCell0[4], dCell0[5],
	      dCell0[0], dCell0[1], dCell0[2], dCell0[6], dCell0[7], dCell0[8] };

	CellPatch3.fill(4, dCell3, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell4[12] =
	      { dCell0[0], dCell0[1], dCell0[2], dCell0[6], dCell0[7],
		  dCell0[8], dCell0[9], dCell0[10], dCell0[11], dCell0[3],
		  dCell0[4], dCell0[5] };

	CellPatch4.fill(4, dCell4, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell5[12] =
	  { dCell0[9], dCell0[10], dCell0[11], dCell0[0], dCell0[1], dCell0[2],
	      dCell0[6], dCell0[7], dCell0[8], dCell0[3], dCell0[4], dCell0[5] };

	CellPatch5.fill(4, dCell5, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell6[12] =
	  { dCell0[6], dCell0[7], dCell0[8], dCell0[9], dCell0[10], dCell0[11],
	      dCell0[0], dCell0[1], dCell0[2], dCell0[3], dCell0[4], dCell0[5] };

	CellPatch6.fill(4, dCell6, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell7[12] =
	  { dCell0[3], dCell0[4], dCell0[5], dCell0[9], dCell0[10], dCell0[11],
	      dCell0[6], dCell0[7], dCell0[8], dCell0[0], dCell0[1], dCell0[2] };

	CellPatch7.fill(4, dCell7, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell8[12] =
	      { dCell0[6], dCell0[7], dCell0[8], dCell0[3], dCell0[4],
		  dCell0[5], dCell0[9], dCell0[10], dCell0[11], dCell0[0],
		  dCell0[1], dCell0[2] };

	CellPatch8.fill(4, dCell8, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell9[12] =
	  { dCell0[9], dCell0[10], dCell0[11], dCell0[6], dCell0[7], dCell0[8],
	      dCell0[3], dCell0[4], dCell0[5], dCell0[0], dCell0[1], dCell0[2] };

	CellPatch9.fill(4, dCell9, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			Error);

	double dCell10[12] =
	      { dCell0[6], dCell0[7], dCell0[8], dCell0[0], dCell0[1],
		  dCell0[2], dCell0[3], dCell0[4], dCell0[5], dCell0[9],
		  dCell0[10], dCell0[11] };

	CellPatch10.fill(4, dCell10, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			 Error);

	double dCell11[12] =
	      { dCell0[3], dCell0[4], dCell0[5], dCell0[6], dCell0[7],
		  dCell0[8], dCell0[0], dCell0[1], dCell0[2], dCell0[9],
		  dCell0[10], dCell0[11] };

	CellPatch11.fill(4, dCell11, 1, Mesquite::TETRAHEDRON, Connect, Fixed,
			 Error);

	Mesquite::TQualityMetric *pTQ[8] =
	  { &TQMetric, &TQMetric1, &TQMetric2, &TQMetric3, &TQMetric4,
	      &TQMetric5, &TQMetric6, &TQMetric7 };
	Mesquite::PatchData *pCP[12] =
	  { &CellPatch, &CellPatch1, &CellPatch2, &CellPatch3, &CellPatch4,
	      &CellPatch5, &CellPatch6, &CellPatch7, &CellPatch8, &CellPatch9,
	      &CellPatch10, &CellPatch11 };

	for (int i = 0; i < 12; i++) {
	  for (int ii = 0; ii < 8; ii++) {
	    pTQ[ii]->evaluate(*pCP[i], 0, dQual, Error);
	    dMinQual = dQual < dMinQual ? dQual : dMinQual;
	  }
	}
	dQ = dMinQual;
	delete pM_QMetric;

	return dQ;
      }
      break;
    case CellSkel::ePyr:
    case CellSkel::ePyrCV:
      assert(0);  ///pending
      break;
    case CellSkel::ePrism:
    case CellSkel::ePrismCV:
      assert(0);  ///pending
      break;
    case CellSkel::eHex:
    case CellSkel::eHexCV:
      assert(0);  ///pending
      break;
    default:
      assert(0);
      break;
    }
  return dQ;
}

void
TMOPQual::vRLamdaRt2D(double *adVec, double *adEig, double *adM)
{

  double dA[4];

  dA[0] = adVec[0] * adEig[0] + adVec[1] * adEig[2];
  dA[1] = adVec[0] * adEig[1] + adVec[1] * adEig[3];
  dA[2] = adVec[2] * adEig[0] + adVec[3] * adEig[2];
  dA[3] = adVec[2] * adEig[1] + adVec[3] * adEig[3];

  adM[0] = dA[0] * adVec[0] + dA[1] * adVec[1];
  adM[1] = dA[0] * adVec[2] + dA[1] * adVec[3];
  adM[2] = dA[2] * adVec[0] + dA[3] * adVec[1];
  adM[3] = dA[2] * adVec[2] + dA[3] * adVec[3];
}

void
TMOPQual::vRLamdaRt3D(double *adVec, double *adEig, double *adM)
{

  double dA[9];

  dA[0] = adVec[0] * adEig[0] + adVec[1] * adEig[3] + adVec[2] * adEig[6];
  dA[1] = adVec[0] * adEig[1] + adVec[1] * adEig[4] + adVec[2] * adEig[7];
  dA[2] = adVec[0] * adEig[2] + adVec[1] * adEig[5] + adVec[2] * adEig[8];
  dA[3] = adVec[3] * adEig[0] + adVec[4] * adEig[3] + adVec[5] * adEig[6];
  dA[4] = adVec[3] * adEig[1] + adVec[4] * adEig[4] + adVec[5] * adEig[7];
  dA[5] = adVec[3] * adEig[2] + adVec[4] * adEig[5] + adVec[5] * adEig[8];
  dA[6] = adVec[6] * adEig[0] + adVec[7] * adEig[3] + adVec[8] * adEig[6];
  dA[7] = adVec[6] * adEig[1] + adVec[7] * adEig[4] + adVec[8] * adEig[7];
  dA[8] = adVec[6] * adEig[2] + adVec[7] * adEig[5] + adVec[8] * adEig[8];

  adM[0] = dA[0] * adVec[0] + dA[1] * adVec[1] + dA[2] * adVec[2];
  adM[1] = dA[0] * adVec[3] + dA[1] * adVec[4] + dA[2] * adVec[5];
  adM[2] = dA[0] * adVec[6] + dA[1] * adVec[7] + dA[2] * adVec[8];
  adM[3] = dA[3] * adVec[0] + dA[4] * adVec[1] + dA[5] * adVec[2];
  adM[4] = dA[3] * adVec[3] + dA[4] * adVec[4] + dA[5] * adVec[5];
  adM[5] = dA[3] * adVec[6] + dA[4] * adVec[7] + dA[5] * adVec[8];
  adM[6] = dA[6] * adVec[0] + dA[7] * adVec[1] + dA[8] * adVec[2];
  adM[7] = dA[6] * adVec[3] + dA[7] * adVec[4] + dA[8] * adVec[5];
  adM[8] = dA[6] * adVec[6] + dA[7] * adVec[7] + dA[8] * adVec[8];
}

void
TMOPQual::vLog2D(double *adM, double *adMLog)
{

  double adVec[4];
  double adEig_[2];
  double adEig[4];
  double adMLog_[4];

  EigenDecomposition2D(adM, adVec, adEig_);

  adEig[0] = log(adEig_[0]);
  adEig[1] = 0;
  adEig[2] = 0;
  adEig[3] = log(adEig_[1]);

  vRLamdaRt2D(adVec, adEig, adMLog_);

  adMLog[0] = adMLog_[0];
  adMLog[1] = adMLog_[1];
  adMLog[2] = adMLog_[3];
}

void
TMOPQual::vExp2D(double *adM, double *adMExp)
{

  double adVec[4];
  double adEig_[2];
  double adEig[4];
  double adMExp_[4];

  EigenDecomposition2D(adM, adVec, adEig_);

  adEig[0] = exp(adEig_[0]);
  adEig[1] = 0;
  adEig[2] = 0;
  adEig[3] = exp(adEig_[1]);

  vRLamdaRt2D(adVec, adEig, adMExp_);

  adMExp[0] = adMExp_[0];
  adMExp[1] = adMExp_[1];
  adMExp[2] = adMExp_[3];

}

void
TMOPQual::vLog3D(double *adM, double *adMLog)
{

  double adVec[9];
  double adEig_[3];
  double adEig[9];
  double adMLog_[9];

  getEigenDecomposition(adM, adVec, adEig_);

  adEig[0] = log(adEig_[0]);
  adEig[1] = 0;
  adEig[2] = 0;
  adEig[3] = 0;
  adEig[4] = log(adEig_[1]);
  adEig[5] = 0;
  adEig[6] = 0;
  adEig[7] = 0;
  adEig[8] = log(adEig_[2]);

  vRLamdaRt3D(adVec, adEig, adMLog_);

  adMLog[0] = adMLog_[0];
  adMLog[1] = adMLog_[3];
  adMLog[2] = adMLog_[4];
  adMLog[3] = adMLog_[6];
  adMLog[4] = adMLog_[7];
  adMLog[5] = adMLog_[8];

}

void
TMOPQual::vExp3D(double *adM, double *adMExp)
{

  double adVec[9];
  double adEig_[3];
  double adEig[9];
  double adMExp_[9];

  getEigenDecomposition(adM, adVec, adEig_);

  adEig[0] = exp(adEig_[0]);
  adEig[1] = 0;
  adEig[2] = 0;
  adEig[3] = 0;
  adEig[4] = exp(adEig_[1]);
  adEig[5] = 0;
  adEig[6] = 0;
  adEig[7] = 0;
  adEig[8] = exp(adEig_[2]);

  vRLamdaRt3D(adVec, adEig, adMExp_);

  adMExp[0] = adMExp_[0];
  adMExp[1] = adMExp_[3];
  adMExp[2] = adMExp_[4];
  adMExp[3] = adMExp_[6];
  adMExp[4] = adMExp_[7];
  adMExp[5] = adMExp_[8];
}

void
TMOPQual::LogEuclideanInterpMetric(Vert* pV, std::set<Vert*> setpV)
{

  std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
  double dDist[setpV.size()];
  double dTDist = 0;
  double adM[iNMetricValues];
  double adNewM[iNMetricValues];
  double dVM[iNMetricValues];
  GR_index_t i = 0;

  for (; iter != iterEnd; ++iter) {
    Vert* psetV = (*iter);
    dDist[i] = calcDistanceBetween(pV, psetV);
    dTDist += dDist[i];
    i++;
  }

  iter = setpV.begin();
  iterEnd = setpV.end();

  i = 0;
  for (int ii = 0; ii < iNMetricValues; ii++) {
    dVM[ii] = 0;
    adM[ii] = 0;
  }

  for (; iter != iterEnd; ++iter) {
    Vert* psetV = (*iter);

    double dLogVM[iNMetricValues];

    for (int ii = 0; ii < iNMetricValues; ii++) {
      dVM[ii] = psetV->getMetric(ii);
    }

    if (iNMetricValues == 3)
      vLog2D(dVM, dLogVM);
    if (iNMetricValues == 6)
      vLog3D(dVM, dLogVM);

    for (int ii = 0; ii < iNMetricValues; ii++) {
      adM[ii] += (dDist[i] * 1. / dTDist) * dLogVM[ii];
    }
    i++;
  }

  if (iNMetricValues == 3)
    vExp2D(adM, adNewM);
  if (iNMetricValues == 6)
    vExp3D(adM, adNewM);

  pV->setMetric(adNewM);

  return;
}

void
TMOPQual::LogEuclideanInterpMetric(double *dCoords, double *adNewM,
				   std::set<Vert*> setpV)
{

  std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
  double dDist[setpV.size()];
  double dTDist = 0;
  double adM[iNMetricValues];
  GR_index_t i = 0;

  for (; iter != iterEnd; ++iter) {
    Vert* psetV = (*iter);

    if (iNMetricValues == 3)
      dDist[i] = dDIST2D(dCoords, psetV->getCoords());
    if (iNMetricValues == 6)
      dDist[i] = dDIST3D(dCoords, psetV->getCoords());
    dTDist += dDist[i];
    i++;
  }

  double dVM[iNMetricValues];
  for (int ii = 0; ii < iNMetricValues; ii++) {
    dVM[ii] = 0;
    adM[ii] = 0;
  }

  iter = setpV.begin();
  iterEnd = setpV.end();
  i = 0;
  for (; iter != iterEnd; ++iter) {
    Vert* psetV = (*iter);

    double dLogVM[iNMetricValues];
    for (int ii = 0; ii < iNMetricValues; ii++) {
      dVM[ii] = psetV->getMetric(ii);
    }

    if (iNMetricValues == 3)
      vLog2D(dVM, dLogVM);
    if (iNMetricValues == 6)
      vLog3D(dVM, dLogVM);

    for (int ii = 0; ii < iNMetricValues; ii++) {
      adM[ii] += (dDist[i] * 1. / dTDist) * dLogVM[ii];
    }
    i++;
  }

  if (iNMetricValues == 3)
    vExp2D(adM, adNewM);
  if (iNMetricValues == 6)
    vExp3D(adM, adNewM);
  return;
}

double
TMOPQual::dMLength2D(const double adLoc0[2], const double adLoc1[2],
		     const double adMet[3])
{
  int i;
  int j;
  double dRetLength = 0;
  double adVec[2];
  double a2dM[2][2];

  a2dM[0][0] = adMet[0];
  a2dM[1][1] = adMet[2];
  a2dM[0][1] = adMet[1];
  a2dM[1][0] = adMet[1];

  for (i = 0; i < 2; i++)
    adVec[i] = adLoc1[i] - adLoc0[i];

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dRetLength += a2dM[i][j] * adVec[i] * adVec[j];
    }
  }

  assert(dRetLength >= 0);
  dRetLength = sqrt(dRetLength);

  return dRetLength;
}

double
TMOPQual::dMLength2D(const Face* const pF)
{

  double adMet[iNMetricValues];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet[i] = pF->getVert(0)->getMetric(i) + pF->getVert(1)->getMetric(i);
    adMet[i] /= 2.0;
  }
  return dMLength2D(pF->getVert(0)->getCoords(), pF->getVert(1)->getCoords(),
		    adMet);
}

double
TMOPQual::dLogMetricLength(double *Coords0, double *adMet0, double *Coords1,
			   double *adMet1)
{

  double dML;
  double dML0;
  double dML1;
  double dA;

  dML0 = dMLength2D(Coords0, Coords1, adMet0);
  dML1 = dMLength2D(Coords0, Coords1, adMet1);
  if (dML0 == dML1)
    dML = dML0;
  else {
    dA = std::max(dML0, dML1) / std::min(dML0, dML1);
    dML = std::max(dML0, dML1) * ((dA - 1) / (dA * log(dA)));
  }
  return dML;
}

double
TMOPQual::dLogMetricLength(const Vert *pV0, const Vert *pV1)
{

  double adMet0[3];
  double adMet1[3];
  double dCoords0[2];
  double dCoords1[2];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet0[i] = pV0->getMetric(i);
  }

  for (int i = 0; i < iNMetricValues; i++) {
    adMet1[i] = pV1->getMetric(i);
  }

  dCoords0[0] = pV0->x();
  dCoords0[1] = pV0->y();

  dCoords1[0] = pV1->x();
  dCoords1[1] = pV1->y();

  return dLogMetricLength(dCoords0, adMet0, dCoords1, adMet1);
}

double
TMOPQual::dLogMetricLength(const Face* const pF)
{
  return dLogMetricLength(pF->getVert(0), pF->getVert(1));
}

void
TMOPQual::dLogMetricCell(double *adMet0, double *adMet1, double *adMet2,
			 double *adM)
{

  double adLogM0[3];
  double adLogM1[3];
  double adLogM2[3];
  double adM_[3];

  vLog2D(adMet0, adLogM0);
  vLog2D(adMet1, adLogM1);
  vLog2D(adMet2, adLogM2);

  for (int i = 0; i < iNMetricValues; i++) {
    adM_[i] = (adLogM0[i] + adLogM1[i] + adLogM2[i]) / 3.;
  }

  vExp2D(adM_, adM);
}

void
TMOPQual::dLogMetricCell(const Vert *pV0, const Vert *pV1, const Vert *pV2,
			 double *adM)
{

  double adMet0[3];
  double adMet1[3];
  double adMet2[3];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet0[i] = pV0->getMetric(i);
    adMet1[i] = pV1->getMetric(i);
    adMet2[i] = pV2->getMetric(i);
  }

  dLogMetricCell(adMet0, adMet1, adMet2, adM);
}

void
TMOPQual::dLogMetricCell(const Cell *pC, double *adM)
{

  dLogMetricCell(pC->getVert(0), pC->getVert(1), pC->getVert(2), adM);
}

double
TMOPQual::MLength3D(const double adLoc0[3], const double adLoc1[3],
		    const double adMet[6])
{
  int i;
  int j;
  double dRetLength = 0;
  double adVec[3];
  double a3dM[3][3];

  a3dM[0][0] = adMet[0];
  a3dM[0][1] = adMet[1];
  a3dM[0][2] = adMet[3];
  a3dM[1][0] = adMet[1];
  a3dM[1][1] = adMet[2];
  a3dM[1][2] = adMet[4];
  a3dM[2][0] = adMet[3];
  a3dM[2][1] = adMet[4];
  a3dM[2][2] = adMet[5];

  for (i = 0; i < 3; i++)
    adVec[i] = adLoc1[i] - adLoc0[i];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dRetLength += a3dM[i][j] * adVec[i] * adVec[j];
    }
  }
  assert(dRetLength >= 0);
  dRetLength = sqrt(dRetLength);

  return dRetLength;
}

////////////////////////////////////////////////////////////////////////////////
// For Merging Triangular cells into Quad cells

void
TMOPQual::MergeTriMeshToQuadMesh(bool MaxQual, Mesh2D *pM2D, double *dQRange)
{

  bMaxQualperCell = MaxQual;
  EvalMeshTritoQuad(pM2D, dQRange);

  Face *pF;
  for (unsigned int i = 0; i < pM2D->getNumFaces(); i++) {

    pF = pM2D->getFace(i);
    if (pF->isCombineToQuad()) {
      pF->CombineToQuad(false);
      QuadCell *pQC;
      pQC = CombineTrisToQuad(pM2D, pF);
      assert(pQC->doFullCheck());
    }
  }
}

void
TMOPQual::EvalMeshTritoQuad(Mesh2D *pM2D, double *QualRange)
{

  Cell *pC;
  int NumQuads;

  do {
    NumQuads = 0;
    for (unsigned int i = 0; i < pM2D->getNumCells(); i++) {
      pC = pM2D->getCell(i);
      if (pC->isCombineToQuad())
	continue;
      if (EvalCellToQuad(pC, QualRange))
	NumQuads++;
    }
  }
  while (NumQuads > 0);
}

bool
TMOPQual::EvalCellToQuad(Cell *pC, double *dQRange)
{

  double dQualRange;
  bool quad = 0;
  double EvalFace[3];
  Face *pF[3];

  if (!dQRange)
    dQualRange = 1e9;
  else
    dQualRange = *dQRange;

  for (int i = 0; i < 3; i++) {
    pF[i] = pC->getFace(i);

    EvalFace[i] = EvalTriToQuad(pF[i], NULL, false);

    if (bMaxQualperCell == false) {

      if (EvalFace[i] <= dQualRange) {
	pF[i]->getLeftCell()->CombineToQuad(1);
	pF[i]->getRightCell()->CombineToQuad(1);
	pF[i]->CombineToQuad(1);
	quad = true;
	continue;
      }
    }
  }

  if ((EvalFace[0] > dQualRange) && (EvalFace[1] > dQualRange)
      && (EvalFace[2] > dQualRange)) {
    return false;
  }

  if (bMaxQualperCell) {
    if (EvalFace[0] <= EvalFace[1] && EvalFace[0] <= EvalFace[2]) {
      if (CheckOppositeCellToQuad(pF[0], pC)) {
	pF[0]->getLeftCell()->CombineToQuad(1);
	pF[0]->getRightCell()->CombineToQuad(1);
	quad = true;
      }
    }

    if (EvalFace[1] < EvalFace[0] && EvalFace[1] <= EvalFace[2]) {
      if (CheckOppositeCellToQuad(pF[1], pC)) {
	pF[1]->getLeftCell()->CombineToQuad(1);
	pF[1]->getRightCell()->CombineToQuad(1);
	quad = true;
      }
    }

    if (EvalFace[2] < EvalFace[0] && EvalFace[2] < EvalFace[1]) {
      if (CheckOppositeCellToQuad(pF[2], pC)) {
	pF[2]->getLeftCell()->CombineToQuad(1);
	pF[2]->getRightCell()->CombineToQuad(1);
	quad = true;
      }
    }
  }

  return quad;
}

bool
TMOPQual::EvalFaceToMerge(Face *pFtoMerge, Face *pF1, Face *pF2)
{

  bool quad = false;

  double EvalFtoMerge = EvalTriToQuad(pFtoMerge, NULL, false);
  double EvalFace1 = EvalTriToQuad(pF1, NULL, false);
  double EvalFace2 = EvalTriToQuad(pF2, NULL, false);

  if (EvalFtoMerge <= EvalFace1 && EvalFtoMerge <= EvalFace2) {
    pFtoMerge->CombineToQuad(1);
    quad = true;
  }

  return quad;
}

bool
TMOPQual::CheckOppositeCellToQuad(Face *pF, Cell *pC)
{

  Face *pF1 = pF->getOppositeCell(pC)->getOppositeFace(pF->getVert(0));
  Face *pF2 = pF->getOppositeCell(pC)->getOppositeFace(pF->getVert(1));
  bool Quad = EvalFaceToMerge(pF, pF1, pF2);

  return Quad;
}

double
TMOPQual::EvalQuadCell(QuadCell* pCell, double *NewCellCoords, bool ForNewCell)
{

  double Qual;
  Vert *pV0 = pCell->getVert(0);
  Vert *pV1 = pCell->getVert(1);
  Vert *pV2 = pCell->getVert(2);
  Vert *pV3 = pCell->getVert(2);

  Qual = EvalTriToQuad(pV0, pV1, pV2, pV3, NewCellCoords, ForNewCell);

  return Qual;
}

double
TMOPQual::EvalTriToQuad(Face *pFace, double *NewCellCoords, bool ForNewCell)
{

  assert(pFace->doFullCheck());
  if (pFace->isBdryFace())
    return 1e10;

  double Qual;
  Cell* pCell0 = pFace->getLeftCell();
  Cell* pCell1 = pFace->getRightCell();

  assert(pCell0->doFullCheck());
  assert(pCell1->doFullCheck());
  assert(pCell0->getRegion() == pCell1->getRegion());

  if (pCell0->isCombineToQuad() || pCell1->isCombineToQuad()) {
    Qual = 1e10;
    return Qual;
  }

  Vert *pV0 = pFace->getVert(0);
  Vert *pV2 = pFace->getVert(1);
  Vert *pV3 = pCell0->getOppositeVert(pFace);
  Vert *pV1 = pCell1->getOppositeVert(pFace);

  //checking old triangles
  assert(checkOrient2D(pV0, pV2, pV3) == 1);
  assert(checkOrient2D(pV0, pV1, pV2) == 1);

  if (!CheckConvexQuad(pV0, pV1, pV2, pV3)) {
    Qual = 1e10;
    return Qual;
  }

  Qual = EvalTriToQuad(pV0, pV1, pV2, pV3, NewCellCoords, ForNewCell);

  return Qual;

}

double
TMOPQual::EvalTriToQuad(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3,
			double *NewCellCoords, bool ForNewCell)
{

  double CellCoord[12];
  double Coord[12];
  double Target[12];
  double TargetLenght;

  SetQuadTargetfromMetric(pV0, pV1, pV2, pV3, Target);
  TargetLenght = sqrt(
      pow((Target[0] - Target[3]), 2) + pow((Target[1] - Target[4]), 2));

  if (ForNewCell) {
    for (int i = 0; i < 12; i++) {
      Coord[i] = NewCellCoords[i];
    }
  }
  else {
    Coord[0] = pV0->x();
    Coord[1] = pV0->y();
    Coord[2] = pV0->z();
    Coord[3] = pV1->x();
    Coord[4] = pV1->y();
    Coord[5] = pV1->z();
    Coord[6] = pV2->x();
    Coord[7] = pV2->y();
    Coord[8] = pV2->z();
    Coord[9] = pV3->x();
    Coord[10] = pV3->y();
    Coord[11] = pV3->z();
  }

  double lenghtf0 = fabs(
      TargetLenght
	  - sqrt(
	      pow((Coord[0] - Coord[3]), 2) + pow((Coord[1] - Coord[4]), 2)));
  double lenghtf1 = fabs(
      TargetLenght
	  - sqrt(
	      pow((Coord[3] - Coord[6]), 2) + pow((Coord[4] - Coord[7]), 2)));
  double lenghtf2 = fabs(
      TargetLenght
	  - sqrt(
	      pow((Coord[6] - Coord[9]), 2) + pow((Coord[7] - Coord[10]), 2)));
  double lenghtf3 = fabs(
      TargetLenght
	  - sqrt(
	      pow((Coord[9] - Coord[0]), 2) + pow((Coord[10] - Coord[1]), 2)));

  int index = 0;
  int count = 0;

  if (lenghtf0 <= lenghtf1 || lenghtf0 <= lenghtf2 || lenghtf0 <= lenghtf3) {
    index = 0;
  }
  else if (lenghtf1 < lenghtf0 || lenghtf1 <= lenghtf2
      || lenghtf1 <= lenghtf3) {
    index = 3;
  }
  else if (lenghtf2 < lenghtf0 || lenghtf2 < lenghtf1 || lenghtf2 <= lenghtf3) {
    index = 6;
  }
  else if (lenghtf3 < lenghtf0 || lenghtf3 < lenghtf1 || lenghtf3 < lenghtf2) {
    index = 9;
  }

  do {
    CellCoord[count] = Coord[index];
    index++;
    if (index >= 12) {
      index = 0;
    }
    count++;
  }
  while (count < 12);

  const int connect[4] =
    { 3, 0, 1, 2 };
  const unsigned long int Connect0[4] =
    { 0, 1, 2, 3 };
  const unsigned long int Connect1[4] =
    { 1, 2, 3, 0 };
  const unsigned long int Connect2[4] =
    { 2, 3, 0, 1 };
  const unsigned long int Connect3[4] =
    { 3, 0, 1, 2 };
  const bool fixed[4] =
    { 0, 0, 0, 0 };

  Mesquite::PatchData CellPatch0;
  Mesquite::PatchData CellPatch1;
  Mesquite::PatchData CellPatch2;
  Mesquite::PatchData CellPatch3;
  Mesquite::MsqError Error;

  Mesquite::MeshImpl Mesh(4, 1, Mesquite::QUADRILATERAL, fixed, Target,
			  connect);
  Mesquite::ReferenceMesh Ref_Mesh(&Mesh);
  Mesquite::RefMeshTargetCalculator Tc(&Ref_Mesh, false);

  //For Isotropic only
  Mesquite::IdealShapeTarget TcIsoIdeal;
  Mesquite::TShapeB1 TShapeOrient;
  Mesquite::TQualityMetric PQMetric(&Tc, &TShapeOrient);

  CellPatch0.fill(4, CellCoord, 1, Mesquite::QUADRILATERAL, Connect0, fixed,
		  Error);
  if (Error)
    std::cout << Error << std::endl;

  CellPatch1.fill(4, CellCoord, 1, Mesquite::QUADRILATERAL, Connect1, fixed,
		  Error);
  if (Error)
    std::cout << Error << std::endl;

  CellPatch2.fill(4, CellCoord, 1, Mesquite::QUADRILATERAL, Connect2, fixed,
		  Error);
  if (Error)
    std::cout << Error << std::endl;

  CellPatch3.fill(4, CellCoord, 1, Mesquite::QUADRILATERAL, Connect3, fixed,
		  Error);
  if (Error)
    std::cout << Error << std::endl;

  double dQual0;
  double dQual1;
  double dQual2;
  double dQual3;
  double dQual;

  PQMetric.evaluate(CellPatch0, 0, dQual0, Error);
  PQMetric.evaluate(CellPatch1, 0, dQual1, Error);
  PQMetric.evaluate(CellPatch2, 0, dQual2, Error);
  PQMetric.evaluate(CellPatch3, 0, dQual3, Error);

  dQual = std::min(std::min(std::min(dQual0, dQual1), dQual2), dQual3);

  return dQual;
}

void
TMOPQual::SetQuadTargetfromMetric(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3,
				  double *Target) const
{

  std::vector<double> lid(4);
  double adMetric[3];
  double adEigVals[2];
  double adTest[4], adTestOri[4];
  double adEigVec1[2], adEigVec2[2];

  adMetric[0] = (pV0->getMetric(0) + pV1->getMetric(0) + pV2->getMetric(0)
      + pV3->getMetric(0)) / 4.0;
  adMetric[1] = (pV0->getMetric(1) + pV1->getMetric(1) + pV2->getMetric(1)
      + pV3->getMetric(1)) / 4.0;
  adMetric[2] = (pV0->getMetric(2) + pV1->getMetric(2) + pV2->getMetric(2)
      + pV3->getMetric(2)) / 4.0;

  //Computing the eigenvalues of the Metric
  double dRoot = sqrt(
      (adMetric[0] - adMetric[2]) * (adMetric[0] - adMetric[2])
	  + 4 * adMetric[1] * adMetric[1]);
  adEigVals[0] = 0.5 * ((adMetric[0] + adMetric[2]) + dRoot); //Larger eigenvalue
  adEigVals[1] = 0.5 * ((adMetric[0] + adMetric[2]) - dRoot); //smaller eigenvalue

  //Computing the normalized eigenvectors of Metric
  if (adMetric[1] < 1e-6 && adMetric[1] > -1e-6) {

    adEigVec1[0] = 1;
    adEigVec1[1] = 0;
    adEigVec2[0] = 0;
    adEigVec2[1] = 1;

    double adRes[2], adDiff[2];
    adRes[0] = adMetric[0] * adEigVec1[0] + adMetric[1] * adEigVec1[1];
    adRes[1] = adMetric[1] * adEigVec1[0] + adMetric[2] * adEigVec1[1];

    adDiff[0] = fabs(adRes[0] - adEigVals[0] * adEigVec1[0]);
    adDiff[1] = fabs(adRes[1] - adEigVals[0] * adEigVec1[1]);

    if (adDiff[0] > 1e-6 || adDiff[1] > 1e-6) {
      adEigVec1[0] = 0;
      adEigVec1[1] = 1;
      adEigVec2[0] = 1;
      adEigVec2[1] = 0;
    }
  }

  else {

    // EigenVector corresponding to larger eigenvalue
    adEigVec1[0] = adMetric[0] - adEigVals[0];
    adEigVec1[1] = -adMetric[1];

    double dNorm = sqrt(
	adEigVec1[0] * adEigVec1[0] + adEigVec1[1] * adEigVec1[1]);
    adEigVec1[0] /= dNorm;
    adEigVec1[1] /= dNorm;

    //EigenVector corresponding to smaller eigenvalue
    adEigVec2[0] = adMetric[0] - adEigVals[1];
    adEigVec2[1] = -adMetric[1];
    dNorm = sqrt(adEigVec2[0] * adEigVec2[0] + adEigVec2[1] * adEigVec2[1]);
    adEigVec2[0] /= dNorm;
    adEigVec2[1] /= dNorm;

  }

  //setting the Target based on the eigenvectors and eigenvalues
  adTest[0] = adEigVec1[0] * (1. / sqrt(adEigVals[0]));
  adTest[2] = adEigVec1[1] * (1. / sqrt(adEigVals[0]));

  adTest[1] = adEigVec2[0] * (1. / sqrt(adEigVals[1]));
  adTest[3] = adEigVec2[1] * (1. / sqrt(adEigVals[1]));

  //testing the orientation of cell
  if (adTest[0] * adTest[3] > adTest[1] * adTest[2]) {

    adTestOri[0] = adTest[0];
    adTestOri[1] = adTest[1];
    adTestOri[2] = adTest[2];
    adTestOri[3] = adTest[3];
  }
  else {

    adTestOri[0] = adTest[1];
    adTestOri[1] = adTest[0];
    adTestOri[2] = adTest[3];
    adTestOri[3] = adTest[2];
  }

  //testing the norm of the first column
  if (adTestOri[1] * adTestOri[1] + adTestOri[3] * adTestOri[3]
      > adTestOri[0] * adTestOri[0] + adTestOri[2] * adTestOri[2]) {

    lid[0] = -adTestOri[1];
    lid[1] = adTestOri[0];
    lid[2] = -adTestOri[3];
    lid[3] = adTestOri[2];
  }
  else {

    lid[0] = adTestOri[0];
    lid[1] = adTestOri[1];
    lid[2] = adTestOri[2];
    lid[3] = adTestOri[3];

  }

  Target[0] = 0;
  Target[1] = 0;
  Target[2] = 0;
  Target[3] = lid[0];
  Target[4] = lid[2];
  Target[5] = 0;
  Target[6] = lid[0] + lid[1];
  Target[7] = lid[2] + lid[3];
  Target[8] = 0;
  Target[9] = lid[1];
  Target[10] = lid[3];
  Target[11] = 0;

}

bool
TMOPQual::CheckConvexQuad(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3)
{

  Vert *pV[4] =
    { pV0, pV1, pV2, pV3 };
  bool sign = 0;
  int n = 4; //number of vertices
  for (int i = 0; i < n; i++) {
    double dx1 = pV[(i + 2) % n]->x() - pV[(i + 1) % n]->x();
    double dy1 = pV[(i + 2) % n]->y() - pV[(i + 1) % n]->y();
    double dx2 = pV[i]->x() - pV[(i + 1) % n]->x();
    double dy2 = pV[i]->y() - pV[(i + 1) % n]->y();
    double zcrossproduct = dx1 * dy2 - dy1 * dx2;
    if (i == 0)
      sign = zcrossproduct > 0;
    else {
      if (sign != (zcrossproduct > 0))
	return false;
    }
  }
  return sign;

}

QuadCell *
TMOPQual::CombineTrisToQuad(Mesh2D *pM2D, Face *pF)
{

  Cell* cell[2];
  cell[0] = pF->getLeftCell();
  cell[1] = pF->getRightCell();

  Vert *vertA = pF->getVert(0);
  Vert *vertB = pF->getVert(1);
  Vert *vertC = cell[0]->getOppositeVert(pF);
  Vert *vertD = cell[1]->getOppositeVert(pF);

  //checking old triangles
  assert(checkOrient2D(vertA, vertB, vertC) == 1);
  assert(checkOrient2D(vertA, vertD, vertB) == 1);
  assert(cell[0]->getRegion() == cell[1]->getRegion());

  QuadCell *pQCell;
  int iReg = cell[0]->getRegion();

  pM2D->deleteCell(cell[0]);
  pM2D->deleteCell(cell[1]);

  pQCell = pM2D->createQuadCell(vertA, vertD, vertB, vertC, iReg);

  if (pQCell->isValid())
    pQCell->canonicalizeFaceOrder();
  else
    assert(0);

  return pQCell;
}

///////////////////////////////////////////////////////////////////////////////////////////

#endif

