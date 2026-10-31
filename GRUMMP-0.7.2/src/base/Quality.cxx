#include <math.h>
#include "GR_assert.h"
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Geometry.h"
#include "GR_Face.h"
#include "GR_Mesh.h"
#include "GR_Quality.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"
#include "GR_Aniso.h"

void
calcCircumcenter(const double *adLoc0, const double *adLoc1,
		 const double *adLoc2, double adRes[])
{

  double adRow0[] =
    { adLoc0[0] - adLoc1[0], adLoc0[1] - adLoc1[1] };
  double adRow1[] =
    { adLoc0[0] - adLoc2[0], adLoc0[1] - adLoc2[1] };

  double dRHS0 = 0.5 * (dDOT2D (adLoc0, adLoc0) - dDOT2D(adLoc1, adLoc1));
  double dRHS1 = 0.5 * (dDOT2D (adLoc0, adLoc0) - dDOT2D(adLoc2, adLoc2));

  double dDet = adRow0[0] * adRow1[1] - adRow1[0] * adRow0[1];
  double dDet0 = dRHS0 * adRow1[1] - dRHS1 * adRow0[1];
  double dDet1 = dRHS1 * adRow0[0] - dRHS0 * adRow1[0];

  adRes[0] = dDet0 / dDet;
  adRes[1] = dDet1 / dDet;
}

//@ Cell quality functions

//@@ Normalized ratio of shortest edge length to circumradius
void
calcShortestToRadius(const CellSkel* const pC, int * const piN, double adRes[])
// Range:  Normalized to be 0 to 1
{
  assert(
      pC->getType() == Cell::eTet || pC->getType() == Cell::eTetCV
	  || pC->getType() == Cell::eTriCell || pC->getType() == Cell::eTriCV);
  assert(piN != NULL);
  assert(adRes != NULL);
  *piN = 1;
  if (pC->getType() == Cell::eTriCell || pC->getType() == Cell::eTriCV) {
    calcShortestToRadius(pC->getVert(0)->getCoords(),
			 pC->getVert(1)->getCoords(),
			 pC->getVert(2)->getCoords(), piN, adRes);
  }
  else {
//      double dLen01_Sq = dDIST_SQ_3D(pC->pVVert(0)->adCoords(),
//  				   pC->pVVert(1)->adCoords());
//      double dLen02_Sq = dDIST_SQ_3D(pC->pVVert(0)->adCoords(),
//  				   pC->pVVert(2)->adCoords());
//      double dLen03_Sq = dDIST_SQ_3D(pC->pVVert(0)->adCoords(),
//  				   pC->pVVert(3)->adCoords());
//      double dLen12_Sq = dDIST_SQ_3D(pC->pVVert(1)->adCoords(),
//  				   pC->pVVert(2)->adCoords());
//      double dLen13_Sq = dDIST_SQ_3D(pC->pVVert(1)->adCoords(),
//  				   pC->pVVert(3)->adCoords());
//      double dLen23_Sq = dDIST_SQ_3D(pC->pVVert(2)->adCoords(),
//  				   pC->pVVert(3)->adCoords());
    // For an equilateral tetrahedron with edge length 1, the circumradius
    // is sqrt(3/8), so we need this factor to normalize.
    assert(pC->isValid());
    const double * const adA = pC->getVert(0)->getCoords();
    const double * const adB = pC->getVert(1)->getCoords();
    const double * const adC = pC->getVert(2)->getCoords();
    const double * const adD = pC->getVert(3)->getCoords();

    calcShortestToRadius(adA, adB, adC, adD, piN, adRes);
  }
}

void
calcShortestToRadius(const double *adLoc0, const double *adLoc1,
		     const double *adLoc2, int * const piN, double adRes[])
// Range:  Normalized to be 0 to 1
{
  assert(piN != NULL);
  assert(adRes != NULL);
  *piN = 1;

  double dLen01 = dDIST2D(adLoc0, adLoc1);
  double dLen12 = dDIST2D(adLoc1, adLoc2);
  double dLen20 = dDIST2D(adLoc2, adLoc0);
  double dShortest = min(min(dLen01, dLen12), dLen20);
  // For an equilateral triangle with edge length 1, the circumradius
  // is sqrt(1/3), so we need this factor to normalize.
  double dRadius;
  if (checkOrient2D(adLoc0, adLoc1, adLoc2) == 0) {
    // Tri is totally flat.
    dRadius = LARGE_DBL;
  }
  else {
    double adTemp[2];
    calcCircumcenter(adLoc0, adLoc1, adLoc2, adTemp);
    dRadius = dDIST2D(adTemp, adLoc0);
  }
  adRes[0] = sqrt(1. / 3.) * dShortest / dRadius;
}

void
calcShortestToRadius(const double *adA, const double *adB, const double *adC,
		     const double *adD, int * const piN, double adRes[])
{

  assert(piN != NULL);
  assert(adRes != NULL);
  *piN = 1;

  double adCircCent[3];

  double adRow1[] =
    { adA[0] - adB[0], adA[1] - adB[1], adA[2] - adB[2] };
  double adRow2[] =
    { adA[0] - adC[0], adA[1] - adC[1], adA[2] - adC[2] };
  double adRow3[] =
    { adA[0] - adD[0], adA[1] - adD[1], adA[2] - adD[2] };

  double dRHS1 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adB, adB));
  double dRHS2 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adC, adC));
  double dRHS3 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adD, adD));

  solve3By3(adRow1, adRow2, adRow3, dRHS1, dRHS2, dRHS3, adCircCent);

  double dLenSq_AB = dDIST_SQ_3D(adA, adB);
  double dLenSq_AC = dDIST_SQ_3D(adA, adC);
  double dLenSq_AD = dDIST_SQ_3D(adA, adD);
  double dLenSq_BC = dDIST_SQ_3D(adB, adC);
  double dLenSq_BD = dDIST_SQ_3D(adB, adD);
  double dLenSq_CD = dDIST_SQ_3D(adC, adD);

  double dShortestSq = min(
      min(min(dLenSq_AB, dLenSq_AC), min(dLenSq_AD, dLenSq_BC)),
      min(dLenSq_BD, dLenSq_CD));
  double dShortestLen = sqrt(dShortestSq);
  double dRadius = dDIST3D(adCircCent, adA);

  adRes[0] = sqrt(0.375) * dShortestLen / dRadius;
}

static double
dQKolSmirn(const double dLam)
// Computes the Kolmogorov-Smirnov significance, Q_KS
// 1 - Q_KS is the confidence with which we can reject the null
// hypothesis that the distributions are the same.
{
  // If the argument lambda is small, then the two distributions are
  // statistically the same.
  if (dLam < 0.1)
    return 1;
  // If the argument is large, the two are unquestionably different.
  else if (dLam > 4)
    return 0;

  double dSum = 0, d2LamSq = dLam * dLam * 2.;
  // Compute enough terms that the last will be smaller than 1.e-10
  int iMaxTerms = 1 + int(ceil(sqrt(25. / d2LamSq)));
  for (int i = 1; i <= iMaxTerms; i += 2)
    dSum += exp(-i * i * d2LamSq) - exp(-(i + 1) * (i + 1) * d2LamSq);
  return (2 * dSum);
}

Quality::Quality(Mesh* const pMIn, const int iMeas) :
    pM(pMIn), pQM(nullptr), pvfQualityFunc(NULL), iQualEvals(0), iMaxEvals(20), iNBins(
	0), dMaxQual(0), dMinQual(0), dScale(1), aiNumQualVals(
    NULL), adSignifImp(NULL), a2dBinFrac(NULL), adMinQual(NULL), adMaxQual(
    NULL), eET(eAll)
{
  static const int aiNBins[] =
    { 30, 36, 36, 30, 16, 30, 25, 25, 25, 25, 25 };
  static const double adMax[] =
    { 180, 72, 180, 360, 32, 360, 1, 1, 1, 1, 1, 1 };
  static const double adMin[] =
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//  static const bool aqMaxMin[] =
//    { false, true, true, false, true, true, true, true, true, true, true, true };
  static const double adScale[] =
    { 180 / M_PI, 180 / M_PI, 180 / M_PI, 180 / M_PI, 180 / M_PI, 180 / M_PI, 1,
	1, 1, 1, 1, 1 };
  assert(iMeas >= 0 && iMeas < 12);
  int dim = (pM->getType() == Mesh::eMesh2D ? 2 : 3);
  if (dim == 2) {
    switch (iMeas)
      {
      case 0: // Max dihedral
      case 3: // Max solid
	pQM = new GRUMMP::Angle2D();
	setEvalType(eMax);
	break;
      case 1: // Min dihedral
      case 4: // Min solid
	pQM = new GRUMMP::Angle2D();
	setEvalType(eMin);
	break;
      case 2: // All dihedrals
      case 5: // All solid
	pQM = new GRUMMP::Angle2D();
	break;
      case 6: // Inradius to circumradius
	pQM = new GRUMMP::RadiusRatio2D();
	break;
      case 7: // Marcum's sliver measure
      case 9: // volume to area ratio
      case 11:
	pQM = new GRUMMP::AreaPerim2D();
	break;
      case 8: // Marcum's skewness measure, which is a ratio of volume to
	      // circumradius cubed
	pQM = new GRUMMP::MarcumSkewness3D();
	break;
      case 10: // shortest edge to circumradius
	pQM = new GRUMMP::EdgeLengthToRadius2D();
	setEvalType(eMin);
	break;
      default:
	vFatalError("Invalid quality measure selected", "Quality::Quality()");
	break;
      }
  }
  else {
    switch (iMeas)
      {
      case 0: // Max dihedral
	pQM = new GRUMMP::DihedralAngles();
	setEvalType(eMax);
	break;
      case 1: // Min dihedral
	pQM = new GRUMMP::DihedralAngles();
	setEvalType(eMin);
	break;
      case 2: // All dihedrals
	pQM = new GRUMMP::DihedralAngles();
	break;
      case 3: // Max solid
	pQM = new GRUMMP::SolidAngles();
	setEvalType(eMax);
	break;
      case 4: // Min solid
	pQM = new GRUMMP::SolidAngles();
	setEvalType(eMin);
	break;
      case 5: // All solids
	pQM = new GRUMMP::SolidAngles();
	break;
      case 6: // Inradius to circumradius
	pQM = new GRUMMP::RadiusRatio3D();
	break;
      case 7: // Marcum's sliver measure, which is a ratio of volume to
	      // the cube of RMS edge length
      case 11:
	pQM = new GRUMMP::VolumeLength3D();
	break;
      case 8: // Marcum's skewness measure, which is a ratio of volume to
	      // circumradius cubed
	pQM = new GRUMMP::MarcumSkewness3D();
	break;
      case 9: // volume to area ratio
	pQM = new GRUMMP::VolumeArea3D();
	break;
      case 10: // shortest edge to circumradius
	pQM = new GRUMMP::ShortestToRadius3D();
	break;
      default:
	vFatalError("Invalid quality measure selected", "Quality::Quality()");
	break;
      }
  }
  iNBins = aiNBins[iMeas];
  dMaxQual = adMax[iMeas];
  dMinQual = adMin[iMeas];
  dScale = adScale[iMeas];

  a2dBinFrac = new double*[iMaxEvals];
  adMinQual = new double[iMaxEvals];
  adMaxQual = new double[iMaxEvals];
  adSignifImp = new double[iMaxEvals];
  adSignifImp[0] = 1;
  aiNumQualVals = new int[iMaxEvals];
//	qMaxMin = aqMaxMin[iMeas];
  for (int i = 0; i < iMaxEvals; i++) {
    adMinQual[i] = adMax[iMeas]; // This pre-sets the min/max measured
    adMaxQual[i] = adMin[iMeas]; // quality appropriately for use.
    a2dBinFrac[i] = new double[iNBins];
    for (int ii = 0; ii < iNBins; ii++)
      a2dBinFrac[i][ii] = 0.;
  }
}

Quality::~Quality()
{
  for (int i = 0; i < iMaxEvals; i++)
    delete[] a2dBinFrac[i];
  delete[] a2dBinFrac;
  delete[] adSignifImp;
  delete[] aiNumQualVals;
  delete[] adMinQual;
  delete[] adMaxQual;
  if (pQM)
    delete pQM;
}

int
Quality::iComputeBin(const double dQual) const
{
  int iBin = int(iNBins * (dQual - dMinQual) / (dMaxQual - dMinQual));
  assert(iBin <= iNBins && iBin >= 0);
  if (iBin == iNBins)
    iBin = iNBins - 1;
  return (iBin);
}

void
Quality::vEvaluate()
{
  if (iQualEvals >= iMaxEvals) {
    logMessage(0, "%s\n %s %d\n",
	       "Too many evaluations of mesh mesh quality; out of data space",
	       "Try initializing the Quality object with iMaxEvals > ",
	       iMaxEvals);
    return;
  }
  logMessage(1, "Evaluating quality of mesh - Qual Eval %d ... ", iQualEvals);
  int iSum = 0;
  double dMin = adMinQual[iQualEvals];
  double dMax = adMaxQual[iQualEvals];
  for (int i = pM->getNumCells() - 1; i >= 0; i--) {
    Cell *pC = pM->getCell(i);
    assert(pC->isValid());
    if (pC->isDeleted())
      continue;
    assert(pC->doFullCheck() == 1);
    int iN;
    double adResult[24];
    if (pQM) {
      pQM->evalAllFromCell(iN, adResult, pC);
    }
    else {
      pvfQualityFunc(pC, &iN, adResult);
    }
    if (eET != eAll && iN > 1) {
      if (eET == eMin) {
	adResult[0] = *(min_element(adResult, adResult + iN));
      }
      else {
	adResult[0] = *(max_element(adResult, adResult + iN));
      }
      iN = 1;
    }
    for (int ii = 0; ii < iN; ii++) {
      double dQual = adResult[ii] * dScale;
      dMin = min(dMin, dQual);
      dMax = max(dMax, dQual);
      int iBin = iComputeBin(dQual);
      a2dBinFrac[iQualEvals][iBin]++;
      iSum++;
    }
  }
  aiNumQualVals[iQualEvals] = iSum;

  double dFrac = 1 / double(iSum);
  int ii;
  for (ii = 0; ii < iNBins; ii++)
    a2dBinFrac[iQualEvals][ii] *= dFrac;
  adMinQual[iQualEvals] = dMin;
  adMaxQual[iQualEvals] = dMax;

  // Assess the statistical significance of the change in the quality
  // distribution between this evaluation and the previous one using the
  // Kolmogorov-Smirnov test for binned data.
  if (iQualEvals >= 1) {
    double dDiff = 0, dMaxDiff = 0;
    // Compute the difference in the cumulative distribution bin-by-bin,
    // and keep track of the maximum value.
    for (ii = 0; ii < iNBins; ii++) {
      dDiff = dDiff + a2dBinFrac[iQualEvals - 1][ii] // Old value
      - a2dBinFrac[iQualEvals][ii]; // New value
      if (fabs(dDiff) > dMaxDiff)
	dMaxDiff = fabs(dDiff);
    }
    int iN1 = aiNumQualVals[iQualEvals - 1];
    int iN2 = aiNumQualVals[iQualEvals];
    adSignifImp[iQualEvals] = 1
	- dQKolSmirn(dMaxDiff * sqrt(double(iN1 * iN2) / double(iN1 + iN2)));
  }

  iQualEvals++;

  logMessage(1, "done.\n");
}

double
Quality::dEvaluate(const CellSkel * const pC) const
{
  assert(pC->isValid());
  assert(!pC->isDeleted());
  assert(pC->doFullCheck() == 1);
  int iN;
  double adResult[24];
  pvfQualityFunc(pC, &iN, adResult);

  assert(iN == 1);
  return (adResult[0] * dScale);
}

void
Quality::vWriteToFile(const char strQualFileName[]) const
{
  FILE *pF = fopen(strQualFileName, "w");
  if (pF == NULL)
    vFatalError("Couldn't open output file for quality data.",
		"writing quality data to disk");

  // Begin with a header telling which quality measure and what the
  // min/max were for each pass.
  fprintf(pF, "# Quality measure:  %s\n", pQM->getName().c_str());

  fprintf(pF, "#\n");
  fprintf(pF, "# Qual eval        Min qual      Max qual    Conf imp\n");
//   fprintf(pF, "# Qual eval        Min qual      Max qual\n");
  int ii;
  for (ii = 0; ii < iQualEvals; ii++) {
    fprintf(pF, "#    %3d       %11.6f   %11.6f   %11.6f\n", ii, adMinQual[ii],
	    adMaxQual[ii], adSignifImp[ii]);
//     fprintf(pF, "#    %3d       %11.6f   %11.6f\n",
// 	    ii, adMinQual[ii], adMaxQual[ii]);
  }

  fprintf(pF, "#\n");
  for (int i = 0; i < iNBins; i++) {
    fprintf(pF, "%8.3f",
	    dMinQual + (double(i) + 0.5) / iNBins * (dMaxQual - dMinQual));
    for (ii = 0; ii < iQualEvals; ii++)
      fprintf(pF, "%13.6f", a2dBinFrac[ii][i]);
    fprintf(pF, "\n");
  }
  fclose(pF);
}

void
Quality::vImportQualityData(const double adData[])
{
  if (iQualEvals >= iMaxEvals) {
    logMessage(0, "%s\n %s %d\n",
	       "Too many evaluations of mesh mesh quality; out of data space",
	       "Try initializing the Quality object with iMaxEvals > ",
	       iMaxEvals);
    return;
  }
  logMessage(1, "Importing quality data...");
  int iSum = 0;
  double dMin = adMinQual[iQualEvals];
  double dMax = adMaxQual[iQualEvals];
  for (int i = pM->getNumCells() - 1; i >= 0; i--) {
    double dQual = adData[i];
    dMin = min(dMin, dQual);
    dMax = max(dMax, dQual);
    int iBin = iComputeBin(dQual);
    a2dBinFrac[iQualEvals][iBin]++;
    iSum++;
  }
  aiNumQualVals[iQualEvals] = iSum;

  double dFrac = 1 / double(iSum);
  int ii;
  for (ii = 0; ii < iNBins; ii++)
    a2dBinFrac[iQualEvals][ii] *= dFrac;
  adMinQual[iQualEvals] = dMin;
  adMaxQual[iQualEvals] = dMax;

  // Assess the statistical significance of the change in the quality
  // distribution between this evaluation and the previous one using the
  // Kolmogorov-Smirnov test for binned data.
  if (iQualEvals >= 1) {
    double dDiff = 0, dMaxDiff = 0;
    // Compute the difference in the cumulative distribution bin-by-bin,
    // and keep track of the maximum value.
    for (ii = 0; ii < iNBins; ii++) {
      dDiff = dDiff + a2dBinFrac[iQualEvals - 1][ii] // Old value
      - a2dBinFrac[iQualEvals][ii]; // New value
      if (fabs(dDiff) > dMaxDiff)
	dMaxDiff = fabs(dDiff);
    }
    int iN1 = aiNumQualVals[iQualEvals - 1];
    int iN2 = aiNumQualVals[iQualEvals];
    adSignifImp[iQualEvals] = 1
	- dQKolSmirn(dMaxDiff * sqrt(double(iN1 * iN2) / double(iN1 + iN2)));
  }

  iQualEvals++;

  logMessage(1, "done.\n");
}

