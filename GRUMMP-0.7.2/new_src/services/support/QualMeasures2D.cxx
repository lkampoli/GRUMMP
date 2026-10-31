#include "GR_QualMeasure.h"
#include "GR_Geometry.h"
#include "GR_Quality.h"
#include "GR_Mesh2D.h"
#include "GR_TMOPQual.h"
#include "GR_Vec.h"
#include <math.h>

using namespace GRUMMP;

double
QualMeasure::eval(const Vert* vert0, const Vert* vert1, const Vert* vert2,
		  const Vert* vert3, const Vert* vert4, const Vert* vert5,
		  const Vert* vert6, const Vert* vert7) const
{
  VALIDATE_INPUT(vert0->getSpaceDimen() == m_dim);
  int nValues = -1;
  double values[12];
  evalAllFromVerts(nValues, values, vert0, vert1, vert2, vert3, vert4, vert5,
		   vert6, vert7);
  if (shouldMinimize())
    return -(*(std::max_element(values, values + nValues)));
  else
    return (*(std::min_element(values, values + nValues)));
}

double
QualMeasure::eval(const Cell* pC) const
{
  int nValues = -1;
  double values[12];
  evalAllFromCell(nValues, values, pC);
  if (shouldMinimize())
    return -(*(std::max_element(values, values + nValues)));
  else
    return (*(std::min_element(values, values + nValues)));
}

void
QualMeasure::evalAllFromCell(int& nValues, double* values, const Cell* pC) const
{
  const Vert* apVert[] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
  for (int ii = 0; ii < pC->getNumVerts(); ii++) {
    apVert[ii] = pC->getVert(ii);
  }
  evalAllFromVerts(nValues, values, apVert[0], apVert[1], apVert[2], apVert[3],
		   apVert[4], apVert[5], apVert[6], apVert[7]);
}

void
QualMeasure::evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
			      const Vert* vert1, const Vert* vert2,
			      const Vert* vert3, const Vert* vert4,
			      const Vert* vert5, const Vert* vert6,
			      const Vert* vert7) const
{
  const double* adLoc[] =
    { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
  if (vert0)
    adLoc[0] = vert0->getCoords();
  if (vert1)
    adLoc[1] = vert1->getCoords();
  if (vert2)
    adLoc[2] = vert2->getCoords();
  if (vert3)
    adLoc[3] = vert3->getCoords();
  if (vert4)
    adLoc[4] = vert4->getCoords();
  if (vert5)
    adLoc[5] = vert5->getCoords();
  if (vert6)
    adLoc[6] = vert6->getCoords();
  if (vert7)
    adLoc[7] = vert7->getCoords();

  evalAllFromCoords(nValues, values, adLoc[0], adLoc[1], adLoc[2], adLoc[3],
		    adLoc[4], adLoc[5], adLoc[6], adLoc[7]);
}

void
Angle2D::evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
			  const Vert* vert1, const Vert* vert2,
			  const Vert* vert3, const Vert* vert4,
			  const Vert* vert5, const Vert* vert6,
			  const Vert* vert7) const
{
  assert(!vert4 && !vert5 && !vert6 && !vert7);
  if (vert3->isValid()) {
    // The quad case
    const double * const adA = vert0->getCoords();
    const double * const adB = vert1->getCoords();
    const double * const adC = vert2->getCoords();
    const double * const adD = vert3->getCoords();
    evalAllFromCoords(nValues, values, adA, adB, adC, adD);

  }

#ifdef HAVE_MESQUITE
  else {
    if (bEvalMetric == true) {

      double dML0;
      double dML1;
      double dML2;
      double adMetric[3];

      GRUMMP::TMOPQual TQ(2);
      TQ.dLogMetricCell(vert0, vert1, vert2, adMetric);

//		  for(unsigned int i=0; i<3; i++){
//			  adMetric[i] = (vert0->getMetric(i) + vert1->getMetric(i) +  vert2->getMetric(i) )/3.;
//		  }

      double dCoordsV0[2];
      double dCoordsV1[2];
      double dCoordsV2[2];

      dCoordsV0[0] = vert0->x();
      dCoordsV0[1] = vert0->y();

      dCoordsV1[0] = vert1->x();
      dCoordsV1[1] = vert1->y();

      dCoordsV2[0] = vert2->x();
      dCoordsV2[1] = vert2->y();

      dML0 = TQ.dLogMetricLength(dCoordsV0, adMetric, dCoordsV1, adMetric);
      dML1 = TQ.dLogMetricLength(dCoordsV1, adMetric, dCoordsV2, adMetric);
      dML2 = TQ.dLogMetricLength(dCoordsV2, adMetric, dCoordsV0, adMetric);

//		  AnisoRefinement *m_pAR;
//		  dML0 = m_pAR->dMetricLength(vert0->getCoords(),vert1->getCoords(),adMetric);
//		  dML1 = m_pAR->dMetricLength(vert1->getCoords(),vert2->getCoords(),adMetric);
//		  dML2 = m_pAR->dMetricLength(vert2->getCoords(),vert0->getCoords(),adMetric);

      nValues = 3;
      values[0] = GR_acos(
	  (pow(dML1, 2) + pow(dML2, 2) - pow(dML0, 2)) / (2.0 * dML1 * dML2)); // angle @ 0
      values[1] = GR_acos(
	  (pow(dML0, 2) + pow(dML2, 2) - pow(dML1, 2)) / (2.0 * dML0 * dML2)); // angle @ 1
      values[2] = GR_acos(
	  (pow(dML0, 2) + pow(dML1, 2) - pow(dML2, 2)) / (2.0 * dML0 * dML1)); // angle @ 2

    }
#endif
    else {
      const double * const adA = vert0->getCoords();
      const double * const adB = vert1->getCoords();
      const double * const adC = vert2->getCoords();
      evalAllFromCoords(nValues, values, adA, adB, adC);
    }
#ifdef HAVE_MESQUITE
  }
#endif

}

void
Angle2D::evalAllFromCoords(int& nValues, double* values, const double *adA,
			   const double *adB, const double *adC,
			   const double *adD, const double *adLoc4,
			   const double *adLoc5, const double *adLoc6,
			   const double *adLoc7) const
{

  assert(!adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  if (!adD) {
    // The tri case
    nValues = 3;
    double adNorm0[] =
      { +adB[1] - adA[1], -adB[0] + adA[0] };
    double adNorm1[] =
      { +adC[1] - adB[1], -adC[0] + adB[0] };
    double adNorm2[] =
      { +adA[1] - adC[1], -adA[0] + adC[0] };

    NORMALIZE2D(adNorm0);
    NORMALIZE2D(adNorm1);
    NORMALIZE2D(adNorm2);

    double dArg0 = -dDOT2D(adNorm0, adNorm2);
    double dArg1 = -dDOT2D(adNorm0, adNorm1);
    double dArg2 = -dDOT2D(adNorm1, adNorm2);

    values[0] = GR_acos(dArg0); // angle @ 0
    values[1] = GR_acos(dArg1); // angle @ 1
    values[2] = GR_acos(dArg2); // angle @ 2
  }
  else {
    // The quad case
    nValues = 4;
    double adNorm0[] =
      { +adB[1] - adA[1], -adB[0] + adA[0] };
    double adNorm1[] =
      { +adC[1] - adB[1], -adC[0] + adB[0] };
    double adNorm2[] =
      { +adD[1] - adC[1], -adD[0] + adC[0] };
    double adNorm3[] =
      { +adA[1] - adD[1], -adA[0] + adD[0] };

    NORMALIZE2D(adNorm0);
    NORMALIZE2D(adNorm1);
    NORMALIZE2D(adNorm2);
    NORMALIZE2D(adNorm3);

    double dArg0 = -dDOT2D(adNorm0, adNorm3);
    double dArg1 = -dDOT2D(adNorm0, adNorm1);
    double dArg2 = -dDOT2D(adNorm1, adNorm2);
    double dArg3 = -dDOT2D(adNorm2, adNorm3);

    values[0] = GR_acos(dArg0); // angle @ 0
    values[1] = GR_acos(dArg1); // angle @ 1
    values[2] = GR_acos(dArg2); // angle @ 2
    values[3] = GR_acos(dArg3); // angle @ 3
  }
}

void
EdgeLengthToRadius2D::evalAllFromCoords(int& nValues, double* values,
					const double *adA, const double *adB,
					const double *adC, const double *adLoc3,
					const double *adLoc4,
					const double *adLoc5,
					const double *adLoc6,
					const double *adLoc7) const
{

  assert(!adLoc3 && !adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);

  nValues = 3;

  double adNorm0[] =
    { +adB[1] - adA[1], -adB[0] + adA[0] };
  double adNorm1[] =
    { +adC[1] - adB[1], -adC[0] + adB[0] };
  double adNorm2[] =
    { +adA[1] - adC[1], -adA[0] + adC[0] };

  NORMALIZE2D(adNorm0);
  NORMALIZE2D(adNorm1);
  NORMALIZE2D(adNorm2);

  double dArg0 = -dDOT2D(adNorm0, adNorm2);
  double dArg1 = -dDOT2D(adNorm0, adNorm1);
  double dArg2 = -dDOT2D(adNorm1, adNorm2);
  values[0] = 2 * sin(GR_acos(dArg0));
  values[1] = 2 * sin(GR_acos(dArg1));
  values[2] = 2 * sin(GR_acos(dArg2));
}

void
RadiusRatio2D::evalAllFromCoords(int& nValues, double* values,
				 const double *adA, const double *adB,
				 const double *adC, const double *adLoc3,
				 const double *adLoc4, const double *adLoc5,
				 const double *adLoc6,
				 const double *adLoc7) const
{

  assert(!adLoc3 && !adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  // The tri case only
  nValues = 1;

  double dArea;
  calcNormal2D(adA, adB, adC, &dArea);
  dArea /= 2;
  double dLenA = dDIST2D(adA, adB);
  double dLenB = dDIST2D(adB, adC);
  double dLenC = dDIST2D(adC, adA);

  values[0] = (16 * dArea * dArea
      / ((dLenA + dLenB + dLenC) * (dLenA * dLenB * dLenC)));
}

void
AreaPerim2D::evalAllFromCoords(int& nValues, double* values, const double *adA,
			       const double *adB, const double *adC,
			       const double *adLoc3, const double *adLoc4,
			       const double *adLoc5, const double *adLoc6,
			       const double *adLoc7) const
{

  assert(!adLoc3 && !adLoc4 && !adLoc5 && !adLoc6 && !adLoc7);
  // The tri case only
  nValues = 1;
  // Range:  Normalized to be 0 to 1

  const double RConst = 12. * sqrt(3.); // about 20
  double dArea, dLenA, dLenB, dLenC;
  calcNormal2D(adA, adB, adC, &dArea);
  dArea /= 2;
  dLenA = dDIST2D(adA, adB);
  dLenB = dDIST2D(adB, adC);
  dLenC = dDIST2D(adC, adA);
  double dPerim = dLenA + dLenB + dLenC;
  values[0] = RConst * dArea / (dPerim * dPerim);
}
