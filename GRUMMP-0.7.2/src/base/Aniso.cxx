#include <set>
#include <vector>
#include <string.h>
#include <math.h>

#include "GR_ADT.h"
#include "GR_Aniso.h"
#include "GR_BdryPatch2D.h"
#include "GR_Face.h"
#include "GR_Geometry.h"
#include "GR_GRCurve.h"
#include "GR_LengthAniso2D.h"
#include "GR_TMOPQual.h"

#define MAXCELLS 3

typedef struct MetricInfo__ {
  double dX, dY;
  double dMet[5];
} MetricInfo;

AnisoRefinement::AnisoRefinement(Mesh2D* pM2D_,
				 char strMetricFileName_[FILE_NAME_LEN],
				 double dAlpha_) :
    pM2D(pM2D_), pSI2D(new SwappingInserter2D(pM2D, NULL, NULL)), pvSmoothData(
    NULL), iOrder(2), dAlpha(dAlpha_), iNMetricValues(3)
{
  sprintf(strMetricFileName, "%s", strMetricFileName_);
  FILE *pFMetricFile = fopen(strMetricFileName_, "r");
  vAssignMetricsToVerts(pFMetricFile);
//	for(GR_index_t iV = 0; iV < pM2D_->getNumVerts(); iV++){
//
//		Vert * vert = pM2D_->getVert(iV);
//
//		double metric[5] = {1.0,0.0,1.0,0.0,0.0};
//		//		double xc = 1.0032;
//		//		double yc = 0.5;
//		//		double x = vert->x()-xc;
//		//		double y = vert->y()-yc;
//		//		double theta = atan2(y,x);
//		//		double r = sqrt(x*x + y*y);
//		//		double a = 0.1/(r*r+0.01);
//		//		double b = 1.0/(r*r+0.01);
//		//		metric[0] = a+(b-a)*sin(theta)*sin(theta);
//		//		metric[1] = 0.5*sin(2.*theta)*(a-b);
//		//		metric[2] = b+(a-b)*sin(theta)*sin(theta);
//		//		double x = vert->x();
//		//		double y = vert->y();
//		double a = 100.0;
//		double b = 1.0;
//		double theta = M_PI_4;
//		//		if(x < 0.5){
//		//			metric[0] = a;
//		//			metric[2] = b;
//		//		}
//		//		if ( y < 0.7 ){
//		metric[0] = a+(b-a)*sin(theta)*sin(theta);
//		metric[1] = 0.5*sin(2.*theta)*(a-b);
//		metric[2] = b+(a-b)*sin(theta)*sin(theta);
//		//		}
//		printf("metric %f %f %f\n",metric[0],metric[1],metric[2]);
//		vert->setMetric(metric);
//	}
  vSmoothMetric();
//	 vAssignLengthsToFaces();

  printf("Metrics Assigned\n");
}
// use this one for when the metric has already been assigned
AnisoRefinement::AnisoRefinement(Mesh2D* pM2D_, double dAlpha_) :
    pM2D(pM2D_), m_length(nullptr), pSI2D(
	new SwappingInserter2D(pM2D, NULL, NULL)), pvSmoothData(
    NULL), iOrder(2), dAlpha(dAlpha_), iNMetricValues(3)
{
//  m_length = new LengthAniso2D(LengthAniso::eMesh, 1, 1);
//  m_length->initLength(pM2D);
}
AnisoRefinement::AnisoRefinement(Mesh2D* pM2D_, LengthAniso2D* length_,
				 int iOrder_, double dAlpha_) :
    pM2D(pM2D_), m_length(length_), pSI2D(
	new SwappingInserter2D(pM2D, NULL, length_)), pvSmoothData(
    NULL), iOrder(iOrder_), dAlpha(dAlpha_), iNMetricValues(3)
{
  sprintf(strMetricFileName, "%s", "./METRIC.metric");
  SMinitSmoothing(2, 'f', OPTMS_DEFAULT, OPTMS_DEFAULT, &pvSmoothData);

  m_length->initLength(pM2D);
  //vAssignMetricsToVerts(vFunc);
  // vSmoothMetric();
  // vAssignLengthsToFaces();
  iNMetricValues = 3;
  printf("Metrics Assigned\n");
}
;

void
AnisoRefinement::vAssignMetricsToVerts(FILE *pFMetricFile)
{
//	FILE *pFMetricFile = fopen(strMetricFileName, "r");
  if (pFMetricFile == NULL)
    assert(0);

  MetricInfo MetInfo;
  std::vector<MetricInfo> vecMetricInfo;
  int iDataRead = 0;

  iDataRead = fscanf(pFMetricFile, "%d", &(iOrder));
  assert(iOrder > 1 && iOrder < 4);
  iOrder = 2;

  switch (iNMetricValues)
    {
    case 3:
      while (1) {
	iDataRead = fscanf(pFMetricFile, "%lf %lf %lf %lf %lf ", &(MetInfo.dX),
			   &(MetInfo.dY), &(MetInfo.dMet[0]),
			   &(MetInfo.dMet[1]), &(MetInfo.dMet[2]));
	MetInfo.dMet[3] = 0;
	MetInfo.dMet[4] = 0;
	if (iDataRead != 5)
	  break;
	if (MetInfo.dMet[0] * MetInfo.dMet[2]
	    - MetInfo.dMet[1] * MetInfo.dMet[1] < 0) {
	  printf(
	      "[%le %le %le] det=%le is invalid\n",
	      MetInfo.dMet[0],
	      MetInfo.dMet[1],
	      MetInfo.dMet[2],
	      MetInfo.dMet[0] * MetInfo.dMet[2]
		  - MetInfo.dMet[1] * MetInfo.dMet[1]);
	  assert(0);
	  MetInfo.dMet[0] = 1e-3; //this should never happen
	  MetInfo.dMet[1] = 0;
	  MetInfo.dMet[2] = 1e-3;
	}
	vecMetricInfo.push_back(MetInfo);
      }
      break;
    case 5:
      while (1) {
	iDataRead = fscanf(pFMetricFile, "%le %le %le %le %le %le %le ",
			   &(MetInfo.dX), &(MetInfo.dY), &(MetInfo.dMet[0]),
			   &(MetInfo.dMet[1]), &(MetInfo.dMet[2]),
			   &(MetInfo.dMet[3]), &(MetInfo.dMet[4]));
	if (iDataRead != 7)
	  break;
	vecMetricInfo.push_back(MetInfo);
      }
      break;
    default:
      assert(0);
      break;
    }

  fclose(pFMetricFile);

  //Metric has now been read and is stored in the vector;

  std::set<int> setAlreadyRefined;
  ADT VTree(pM2D->pECGetVerts());
  std::vector<GR_index_t> veciQueryResult;
  for (std::vector<MetricInfo>::iterator Iter = vecMetricInfo.begin();
      Iter != vecMetricInfo.end(); ++Iter) {
    double dEps = 1e-10;
    do {
      double adBBox[] =
	{ Iter->dX - dEps, Iter->dX + dEps, Iter->dY - dEps, Iter->dY + dEps };
      veciQueryResult = VTree.veciRangeQuery(adBBox);
      dEps *= 10;
    }
    while (veciQueryResult.empty());
    size_t resSize = veciQueryResult.size();
    for (size_t i = 0; i < resSize; i++) {
      int iV = veciQueryResult[i];
      if (setAlreadyRefined.find(iV) != setAlreadyRefined.end()) {
	continue;
      }
      setAlreadyRefined.insert(iV);
      Vert *pV = pM2D->getVert(iV);

      logMessage(
	  MSG_DEBUG,
	  "Setting metric from vAssignMetricsToVerts. Terrible mistake.\n");
      pV->setMetric(Iter->dMet);
    }
  }
  size_t iSize = vecMetricInfo.size();
  assert(iSize == size_t(pM2D->iNumVerts()));
  //this can be removed later

  logMessage(1, "Read Metric for %zd locations\n", vecMetricInfo.size());
  logMessage(1, "Mesh has %d vertices\n", pM2D->iNumVerts());

  return;
}

void
AnisoRefinement::vAssignMetricsToVertsTest(enum eMetricType Type)
{
  for (GR_index_t i = 0; i < pM2D->getNumVerts(); i++) {

    Vert *pV = pM2D->getVert(i);
    assert(pV->isValid());

    double adM[3];
    double d_hmin;
    double d_hmax;
    double d_alphaX;
    double d_alphaY;
    double d_alpha;
    double d_h1;
    double d_h2;
    double d_Theta;
    double d_X = pV->x();
    double d_Y = pV->y();

    switch (Type)
      {

      case e2DTestCross:

	d_hmin = 0.005;
	d_hmax = 0.1;
	d_alphaX = 20 * (fabs(d_X - 0.5));
	d_alphaY = 20 * (fabs(d_Y - 0.5));
	d_h1 = std::min((pow(2, d_alphaX)) * d_hmin, d_hmax);
	d_h2 = std::min((pow(2, d_alphaY)) * d_hmin, d_hmax);

	adM[0] = pow(d_h1, -2);
	adM[1] = 0;
	adM[2] = pow(d_h2, -2);

	break;

      case e2DTestQuarterCircle:

	d_Theta = atan2(d_Y, d_X);
	d_hmax = 0.1;
	d_alpha = 10 * (fabs(0.75 - sqrt(pow(d_X, 2) + pow(d_Y, 2))));
	d_h1 = std::min(0.002 * pow(5, d_alpha), d_hmax);
	d_h2 = std::min(0.05 * pow(2, d_alpha), d_hmax);

	adM[0] = pow(d_h1, -2) * pow(cos(d_Theta), 2)
	    + pow(d_h2, -2) * pow(sin(d_Theta), 2);
	adM[1] = (pow(d_h1, -2) - pow(d_h2, -2)) * cos(d_Theta) * sin(d_Theta);
	adM[2] = pow(d_h1, -2) * pow(sin(d_Theta), 2)
	    + pow(d_h2, -2) * pow(cos(d_Theta), 2);

	break;

      default:
	assert(0);
	break;
      }

    pV->setMetric(adM);
  }
}

void
AnisoRefinement::vAssignMetricsToVerts(void
(vFunc)(const double adLoc[], double adRetVal[]))
{
  //Assign the metric based on a function

  double adNewMetric[3];
  Vert* pV;
  GR_index_t i;

  for (i = 0; i < pM2D->iNumVerts(); i++) {
    pV = pM2D->getVert(i);
    vFunc(pV->getCoords(), adNewMetric);
    // printf("vFunc returns [%.3e,%.3e,%.3e]\n",adNewMetric[0],adNewMetric[1],adNewMetric[2]);
    logMessage(
	MSG_DEBUG,
	"Setting metric from vAssignMetricsToVerts. Terrible mistake.\n");
    pV->setMetric(adNewMetric);
  }

  //   for(i=0;i<pM2D->iNumVerts();i++){
  //     pV=pM2D->pVVert(i);
  //     printf("(%.3e,%.3e) ->[%.3e,%.3e,%.3e]\n",pV->dX(),pV->dY(),pV->dGetMetric(0),pV->dGetMetric(1),pV->dGetMetric(2));
  //   }
}

void
AnisoRefinement::vSmoothMetric(Vert* pVVert)
{
  std::set<Cell*> spCInc; //just to make vneighborhood work
  std::set<Vert*> spVNeigh; //
  std::set<Vert*>::iterator iterV;
  std::set<Cell*>::iterator iterC;

  findNeighborhoodInfo(pVVert, spCInc, spVNeigh);

  double adCurrMet[5] =
    { 0, 0, 0, 0, 0 };
  double adAvgMet[5] =
    { 0, 0, 0, 0, 0 };
  int i;
  double dSmoothRat = dAlpha;

  assert(pVVert->isValid());

  for (iterV = spVNeigh.begin(); iterV != spVNeigh.end(); ++iterV) {
    assert((*iterV)->isValid());
    for (i = 0; i < iNMetricValues; i++) {
      adAvgMet[i] += (*iterV)->getMetric(i);
    }
  }
  for (i = 0; i < iNMetricValues; i++) {
    adCurrMet[i] = pVVert->getMetric(i);
    adAvgMet[i] /= double(spVNeigh.size());
  }

  for (i = 0; i < iNMetricValues; i++) {
    adCurrMet[i] = adCurrMet[i] + dSmoothRat * (adAvgMet[i] - adCurrMet[i]);
  }

  if (pM2D->pCellQueue != NULL) {
    for (iterC = spCInc.begin(); iterC != spCInc.end(); iterC++) {
      pM2D->pCellQueue->vRemoveCell(*iterC);
    }
  }

  logMessage(MSG_DEBUG,
	     "Setting metric from qSmoothMetric. Terrible mistake.\n");
  pVVert->setMetric(adCurrMet);

  if (pM2D->pCellQueue != NULL) {
    for (iterC = spCInc.begin(); iterC != spCInc.end(); iterC++) {
      pM2D->pCellQueue->vAddCell(*iterC);
    }
  }
  return;
}

bool
AnisoRefinement::qSmoothMetric(Vert* pV0, const Vert* const pV1)
{
  //Smooth pV0 to match pV1+alpha*len  return whether or not things changed

  double dFudge = 0.99;
  //Given two vertices, increase pV0 to satisfy alpha-lipchitz
  double adMet0[2][2];
  double adMet1[2][2];
  double adMetNew[5] =
    { 0, 0, 0, 0, 0 };
  double adEig0[2];
  double adEig1[2];
  double dSize;
  double dDiag0, dDiag1;
  double dBound0 = 0;
  double dBound1 = 0;
  int i, j;
  double dLen = calcDistanceBetween(pV0, pV1);
  //  double dLen=dMetricLength(pV0->adCoords(),pV1->adCoords(),pV0->adMetricArray());
  // dLen+=dMetricLength(pV0->adCoords(),pV1->adCoords(),pV1->adMetricArray());
  //  dLen*=0.5;
  bool qChange = false;
  double dMinEigValue = 1e-3;
  adMet0[0][0] = pV0->getMetric(0);
  adMet0[0][1] = pV0->getMetric(1);
  adMet0[1][0] = pV0->getMetric(1);
  adMet0[1][1] = pV0->getMetric(2);

  adMet1[0][0] = pV1->getMetric(0);
  adMet1[0][1] = pV1->getMetric(1);
  adMet1[1][0] = pV1->getMetric(1);
  adMet1[1][1] = pV1->getMetric(2);

  printf("Starting Metric is [%.3e,%.3e,%.3e]\n", adMet0[0][0], adMet0[0][1],
	 adMet0[1][1]);
  printf("Compared to [%.3e,%.3e,%.3e] at %.3e away\n", adMet1[0][0],
	 adMet1[0][1], adMet1[1][1], dLen);

  if (iFuzzyComp(adMet0[0][1], 0) == 0) {
    //  printf("Awkward...\n");
    if (iFuzzyComp(adMet1[0][1], 0) == 0) {
      //    printf("But okay in the end\n");
      adEig0[0] = 1;
      adEig0[1] = 0;
      adEig1[0] = 0;
      adEig1[1] = 1;
    }
    else {
      dSize = adMet1[1][1] * adMet1[1][1] - 2 * adMet1[0][0] * adMet1[1][1]
	  + adMet1[0][0] * adMet1[0][0] + 4 * adMet1[0][1] * adMet1[0][1];
      //    printf("dSize^2=%.3e\n",dSize);
      dSize = sqrt(dSize);
      adEig0[0] = 0.5 * (adMet1[0][0] - adMet1[1][1] + dSize) / adMet1[0][1];
      adEig0[1] = 1;
      adEig1[0] = 0.5 * (adMet1[0][0] - adMet1[1][1] - dSize) / adMet1[0][1];
      adEig1[1] = 1;
    }
    // adMet0[0][1]=1e-8;
    //  adMet0[1][0]=1e-8; //to avoid division by zero. It can probably be replaced with something more satisfactory in the future
  }
  else {
    dSize = adMet0[1][1] * adMet0[1][1] - 2 * adMet0[0][0] * adMet0[1][1]
	+ adMet0[0][0] * adMet0[0][0] + 4 * adMet0[0][1] * adMet0[0][1];
    dSize = sqrt(dSize);

    adEig0[0] = 0.5 * (adMet0[0][0] - adMet0[1][1] + dSize) / adMet0[0][1];
    adEig0[1] = 1;

    adEig1[0] = 0.5 * (adMet0[0][0] - adMet0[1][1] - dSize) / adMet0[0][1];
    adEig1[1] = 1;
  }
  dDiag0 = sqrt(adEig0[0] * adEig0[0] + adEig0[1] * adEig0[1]);
  dDiag1 = sqrt(adEig1[0] * adEig1[0] + adEig1[1] * adEig1[1]);

  for (i = 0; i < 2; i++) {
    adEig0[i] /= dDiag0;
    adEig1[i] /= dDiag1;
  } // now the eigenvectors are normalized
    //  printf("Eigenvectors [%.3e,%.3e] and [%.3e,%.3e]\n",adEig0[0],adEig0[1],adEig1[0],adEig1[1]);

  dDiag0 = 0;
  dDiag1 = 0;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dDiag0 += adMet0[i][j] * adEig0[i] * adEig0[j];
    }
  }
  dDiag0 = sqrt(dDiag0);

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dDiag1 += adMet0[i][j] * adEig1[i] * adEig1[j];
    }
  }
  dDiag1 = sqrt(dDiag1);

  if (dDiag0 < dMinEigValue)
    dDiag0 = dMinEigValue;
  if (dDiag1 < dMinEigValue)
    dDiag1 = dMinEigValue;

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dBound0 += adMet1[i][j] * adEig0[i] * adEig0[j];
    }
  }
  dBound0 = sqrt(dBound0);
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dBound1 += adMet1[i][j] * adEig1[i] * adEig1[j];
    }
  }
  dBound1 = sqrt(dBound1);

  //  printf("Comparing (%.3e,%.3e) to (%.3e,%.3e)\n",dDiag0,dDiag1,dBound0,dBound1);
  if ((dBound0 - dDiag0) > (dAlpha * dLen)) {
    dDiag0 = dBound0 - dAlpha * dLen * dFudge;

    qChange = true;
  }

  if ((dBound1 - dDiag1) > (dAlpha * dLen)) {
    dDiag1 = dBound1 - dAlpha * dLen * dFudge;

    qChange = true;
  }
  //  printf("New Diag0 = %.3e\n",dDiag0);
  //  printf("New Diag1 = %.3e\n",dDiag1);

  if (!qChange) {
    //   printf("\n");
    return false;
  }

  assert((dBound0 - dDiag0) <= (dAlpha * dLen));
  assert((dBound1 - dDiag1) <= (dAlpha * dLen));

  double t1, t2, t6, t9, t10, t14, t16, t18;
  t1 = dDiag0 * dDiag0;
  t2 = adEig0[0] * t1;
  t6 = 1 / (adEig0[0] * adEig1[1] - adEig1[0] * adEig0[1]);
  t9 = dDiag1 * dDiag1;
  t10 = adEig1[0] * t9;
  t14 = adEig1[0] * t6;
  t16 = adEig0[0] * t6;
  t18 = -t2 * t14 + t10 * t16;
  adMetNew[0] = t2 * adEig1[1] * t6 - t10 * adEig0[1] * t6;
  adMetNew[1] = t18;
  adMetNew[2] = -adEig0[1] * t1 * t14 + adEig1[1] * t9 * t16;
  // printf("New Metric is [%.3e,%.3e,%.3e]\n\n",adMetNew[0],adMetNew[1],adMetNew[2]);

  logMessage(MSG_DEBUG,
	     "Setting metric from qSmoothMetric. Terrible mistake.\n");
  pV0->setMetric(adMetNew);

  //  assert(!qSmoothMetric(pV0,pV1));
  return true;
}

void
AnisoRefinement::vSmoothMetric()
{
  GR_index_t iV;
  printf("Smoothing Metrics with alpha=%.3e.....", dAlpha);
  for (iV = 0; iV < pM2D->iNumVerts(); iV++) {
    vSmoothMetric(pM2D->getVert(iV));
  }
  printf("Done\n");

  return;
}

double
dSpaceDensity(const Vert* const pV, const double adNorm[2])
{
  double dRetVal;
  const double* adM = pV->getMetricArray();

  dRetVal = adM[0] * adNorm[0] * adNorm[0] + 2 * adNorm[0] * adNorm[1] * adM[1]
      + adM[2] * adNorm[1] * adNorm[1];
  dRetVal = sqrt(dRetVal);
  return dRetVal;
}

//This is a vertex that is already in the mesh so it is interpolated using it's neighbo-u-rs
void
AnisoRefinement::vInterpolateMetric(Vert* pVVert)
{
  assert(pVVert->isValid());

  std::set<Cell*> spCInc;
  std::set<Vert*> spVNeigh;
  std::set<Vert*>::iterator iterV;

  findNeighborhoodInfo(pVVert, spCInc, spVNeigh);
  if (pVVert->isBdryVert()) {
    for (iterV = spVNeigh.begin(); iterV != spVNeigh.end(); ++iterV) {
    	Vert* pV = *iterV;
    	if (!pV->isBdryVert())
    		spVNeigh.erase(iterV);
    }
  }
  // printf("using neighbours to interpolate (%.3e,%.3e)\n",pVVert->dX(),pVVert->dY());
  vInterpolateMetric(pVVert, &(spVNeigh));

  return;
}

void
vLinearInterpolate(const double adLoc1[2], const double adVal1[],
		   const double adLoc2[2], const double adVal2[],
		   const double adLoc3[2], const double adVal3[],
		   const double adNewPoint[2], double adNewValue[],
		   int iNumValues)
{
  // vMessage(1,"Interpolating with 3 vertices\n");
  //  printf("(%.3e,%.3e)[%.3e,%.3e,%.3e]-(%.3e,%.3e)[%.3e,%.3e,%.3e]-(%.3e,%.3e)[%.3e,%.3e,%.3e]\n",adLoc1[0],adLoc1[1],adVal1[0],adVal1[1],adVal1[2],adLoc2[0],adLoc2[1],adVal2[0],adVal2[1],adVal2[2],adLoc3[0],adLoc3[1],adVal3[0],adVal3[1],adVal3[2]);

  double dDetM;
  double dX[3] =
    { adLoc1[0], adLoc2[0], adLoc3[0] };
  double dY[3] =
    { adLoc1[1], adLoc2[1], adLoc3[1] };
  double adNewLoc[3] =
    { 1, adNewPoint[0], adNewPoint[1] };
  double adSol[3];
  int i, j;
  double a2dMi[3][3];

  //   printf("adNewLoc is [%.3e,%.3e,%.3e]\n",adNewLoc[0],adNewLoc[1],adNewLoc[2]);
  dDetM = dX[1] * dY[2] - dY[1] * dX[2] - dX[0] * dY[2] + dY[0] * dX[2]
      + dX[0] * dY[1] - dY[0] * dX[1];
  //   printf("Det is %.3e\n",dDetM);
  if (fabs(dDetM) < 1e-10) {
    //Just use the average
    for (i = 0; i < iNumValues; i++) {
      adNewValue[i] = 1.0 / 3.0 * (adVal1[i] + adVal2[i] + adVal3[i]);
    }
  }
  else {
    a2dMi[0][0] = dX[1] * dY[2] - dY[1] * dX[2];
    a2dMi[0][1] = dX[2] * dY[0] - dX[0] * dY[2];
    a2dMi[0][2] = dX[0] * dY[1] - dY[0] * dX[1];
    a2dMi[1][0] = dY[1] - dY[2];
    a2dMi[1][1] = dY[2] - dY[0];
    a2dMi[1][2] = dY[0] - dY[1];
    a2dMi[2][0] = dX[2] - dX[1];
    a2dMi[2][1] = dX[0] - dX[2];
    a2dMi[2][2] = dX[1] - dX[0];
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) {
	a2dMi[i][j] /= dDetM;
      }
    }
    //       for(i=0;i<3;i++){
    //        printf("a2dMi(%d) is [%.3e,%.3e,%.3e]\n",i,a2dMi[i][0],a2dMi[i][1],a2dMi[i][2]);
    //             }

    for (j = 0; j < 3; j++) {
      adSol[j] = 0;
      for (i = 0; i < 3; i++) {
	adSol[j] += adNewLoc[i] * a2dMi[i][j];
      }
    }
    for (j = 0; j < iNumValues; j++) {
      // adNewValue[j]=0;
      adNewValue[j] = adSol[0] * adVal1[j];
      adNewValue[j] += adSol[1] * adVal2[j];
      adNewValue[j] += adSol[2] * adVal3[j];
    }
  }
  //  printf("[%.3e,%.3e,%.3e] at (%.3e,%.3e)\n",adNewValue[0],adNewValue[1],adNewValue[2],adNewPoint[0],adNewPoint[1]);

  return;
}

//This vertex may not have connectivity yet and so is interpolated as an average of possible containing cells
void
AnisoRefinement::vInterpolateMetric(Vert* pVVert, std::set<Vert*>* psetpVVerts)
{
  std::set<Vert*>::iterator iter;
  Vert* pV0 = NULL;
  Vert* pV1 = NULL;
  Vert* pV2 = NULL;
  Vert* pV3 = NULL;
  Vert* pV4 = NULL;
  Vert* pVa = NULL;
  Vert* pVb = NULL;
  Vert* pVc = NULL;
  double dDist0, dDist1, dDist2, dDist3;
  const double* adM0;
  const double* adM1;
  const double* adM2;
  const double* adM3;

  double dNumToAvg = 0;
  double adNewMetric[iNMetricValues];
  double adTempMetric[iNMetricValues];
  int i;
  int iInTri;
  //Just in case
  for (i = 0; i < iNMetricValues; i++) {
    adNewMetric[i] = 0;
  }

  switch (psetpVVerts->size())
    {
    case 1:
      assert(0);
      //I shouldn't be calling this
      break;
    case 2:
      //Weighted average
      //   assert(0); //Nor this
      //  vMessage(1,"Interpolating along a line\n");
      iter = psetpVVerts->begin();
      pV0 = (*iter);
      ++iter;
      pV1 = (*iter);
      dDist0 = calcDistanceBetween(pV0, pVVert);
      dDist1 = calcDistanceBetween(pV1, pVVert);

      adM0 = pV0->getMetricArray();
      adM1 = pV1->getMetricArray();

      for (i = 0; i < iNMetricValues; i++) {
	adNewMetric[i] = (dDist1 * adM0[i] + dDist0 * adM1[i])
	    / (dDist1 + dDist0);
      }

      logMessage(MSG_DEBUG,
		 "Setting metric from vInterpolateMetric. Terrible mistake.\n");
      pVVert->setMetric(adNewMetric);
      //  printf("New metric is [%.3e,%.3e,%.3e]\n",pVVert->dGetMetric(0),pVVert->dGetMetric(1),pVVert->dGetMetric(2));
      return;

      break;
    case 3:
      //    printf("Interpolating with 3\n");
      iter = psetpVVerts->begin();
      pV1 = *iter;
      ++iter;
      pV2 = *iter;
      ++iter;
      pV3 = *iter;
      if (isInTri(pVVert, pV1, pV2, pV3) < 1) {
	assert(pVVert->isBdryVert());
	dDist1 = calcDistanceBetween(pV1, pVVert);
	dDist2 = calcDistanceBetween(pV2, pVVert);
	dDist3 = calcDistanceBetween(pV3, pVVert);

	adM1 = pV1->getMetricArray();
	adM2 = pV2->getMetricArray();
	adM3 = pV3->getMetricArray();

	for (i = 0; i < iNMetricValues; i++) {
	  adNewMetric[i] = (1.0 / dDist1 * adM1[i] + 1.0 / dDist2 * adM2[i]
	      + 1.0 / dDist3 * adM3[i])
	      / (1.0 / dDist1 + 1.0 / dDist2 + 1.0 / dDist3);
	}
      }
      else {
	vLinearInterpolate(pV1->getCoords(), pV1->getMetricArray(),
			   pV2->getCoords(), pV2->getMetricArray(),
			   pV3->getCoords(), pV3->getMetricArray(),
			   pVVert->getCoords(), adNewMetric, iNMetricValues);
      }
      break;
    case 4:
      // printf("Interpolating with 4\n");
      iter = psetpVVerts->begin();
      pV1 = *iter;
      ++iter;
      pV2 = *iter;
      ++iter;
      pV3 = *iter;
      ++iter;
      pV4 = *iter;
      //   printf("(%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n",pV1->dX(),pV1->dY(),pV2->dX(),pV2->dY(),pV3->dX(),pV3->dY(),pV4->dX(),pV4->dY());
      //   printf("newpoint=(%.3e,%.3e)\n",pVVert->dX(),pVVert->dY());

      dNumToAvg = 0;
      //This could be cleaned up
      for (int iNumPossTri = 0; iNumPossTri < 4; iNumPossTri++) {
	if (iNumPossTri == 0) {
	  pVa = pV1;
	  pVb = pV2;
	  pVc = pV3;
	}
	if (iNumPossTri == 1) {
	  pVa = pV1;
	  pVb = pV2;
	  pVc = pV4;
	}
	if (iNumPossTri == 2) {
	  pVa = pV1;
	  pVb = pV3;
	  pVc = pV4;
	}
	if (iNumPossTri == 3) {
	  pVa = pV2;
	  pVb = pV3;
	  pVc = pV4;
	}
	//	printf("(%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n",pVa->dX(),pVa->dY(),pVb->dX(),pVb->dY(),pVc->dX(),pVc->dY());
	iInTri = isInTri(pVVert, pVa, pVb, pVc);
	//	printf("IsInTri returns %d\n",iInTri);
	if (iInTri == 3 || iInTri == 2) {
	  vLinearInterpolate(pVa->getCoords(), pVa->getMetricArray(),
			     pVb->getCoords(), pVb->getMetricArray(),
			     pVc->getCoords(), pVc->getMetricArray(),
			     pVVert->getCoords(), adTempMetric, iNMetricValues);
	  for (i = 0; i < iNMetricValues; i++) {
	    if (iInTri == 2)
	      adTempMetric[i] *= 0.5; //this is a quick hack so that there isn't a bias towards a configuration with the vertex on a diagonal
	    adNewMetric[i] += adTempMetric[i];
	  }
	  if (iInTri == 2)
	    dNumToAvg += 0.5;
	  else
	    dNumToAvg += 1.0;
	}
      }

      assert(dNumToAvg > 0);
      for (i = 0; i < iNMetricValues; i++)
	adNewMetric[i] /= (dNumToAvg * 1.0);
      //  printf("%.1f possible cells contain the vertex\n",dNumToAvg);
      //    for(i=0;i<iNMetricValues;i++) adNewMetric[i]/=(dNumToAvg*1.0);
      //       printf("New Metric is [%.3e,%.3e,%.3e]\n",adNewMetric[0],adNewMetric[1],adNewMetric[2]);

      break;
    default:
      assert(0);
      break;
    }

  //   printf("New Metric is [%.3e,%.3e,%.3e] at (%.3e,%.3e)\n",adNewMetric[0],adNewMetric[1],adNewMetric[2],pVVert->dX(),pVVert->dY());
  // printf("\n");
  if (adNewMetric[0] < 0
      || (adNewMetric[0] * adNewMetric[2] - adNewMetric[1] * adNewMetric[1])
	  < 0) {

    printf("Interpolating with %zd\n", psetpVVerts->size());
    for (iter = psetpVVerts->begin(); iter != psetpVVerts->end(); ++iter) {
      printf("(%.6e,%.6e)  with [%.6e,%.6e,%.6e]\n", (*iter)->x(), (*iter)->y(),
	     (*iter)->getMetric(0), (*iter)->getMetric(1),
	     (*iter)->getMetric(2));
    }
    printf("New Metric is [%.6e,%.6e,%.6e] at (%.6e,%.6e)\n", adNewMetric[0],
	   adNewMetric[1], adNewMetric[2], pVVert->x(), pVVert->y());
    assert(0);
  }

  logMessage(MSG_DEBUG,
	     "Setting metric from vInterpolateMetric. Terrible mistake.\n");
  pVVert->setMetric(adNewMetric);
  //  vSmoothMetric(pVVert);
  return;
}

void
AnisoRefinement::vAniMidPoint(Face* pF, double adMidPoint[2])
{
  //This function is used to "fix" the area around pF
  //It finds the anisotropic midpoint of the face

  const double* adP1 = pF->getVert(0)->getCoords();
  const double* adP2 = pF->getVert(1)->getCoords();

  double dLen = calcDistanceBetween(pF->getVert(0), pF->getVert(1));
  double dRho0, dRho1, dM0, dM1, dtemp, dMid;
  int i;
  double adNorm[2];

  for (i = 0; i < 2; i++) {
    adNorm[i] = 1.0 * (adP2[i] - adP1[i]) / dLen;
  }

  dRho0 = dSpaceDensity(pF->getVert(0), adNorm);
  dRho1 = dSpaceDensity(pF->getVert(1), adNorm);

  // printf("dRho0 is %.3e, dRho1 is %.3e, dLen is %.3e\n",dRho0,dRho1,dLen);

  if (fabs(dRho0 - dRho1) < 1e-9)
    dMid = 0.5 * dLen;
  else {
    dM0 = dRho0 * dRho0;
    dM1 = dRho1 * dRho1;

    dtemp = pow(4.0, 0.6666666666666667)
	* pow(dRho0 * dM0 + dRho1 * dM1, 0.6666666666666667);
    dMid = dLen * (dM0 - dtemp / 4.0) / (-dM1 + dM0);
  }
  // printf("dMid is %.3e\n",dMid);
  for (i = 0; i < 2; i++) {
    adMidPoint[i] = adP1[i] + dMid * adNorm[i];
  }
  // printf("The aniMidPoint is (%.3e,%.3e)\n",adMidPoint[0],adMidPoint[1]);

  return;
}

double
AnisoRefinement::dMetricLength(const Face* const pF)
{
  // WARNING: this function returns the maximum length of a face, the
  // actual length depends on which cell you are looking from
  double adMet[iNMetricValues];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet[i] = pF->getVert(0)->getMetric(i) + pF->getVert(1)->getMetric(i);
    adMet[i] /= 2.0;
  }

  return dMetricLength(pF->getVert(0)->getCoords(), pF->getVert(1)->getCoords(),
		       adMet);
}

double
AnisoRefinement::dMetricLength(const double adLoc0[2], const double adLoc1[2],
			       const double adMet[5])
{
  int i;
  int j;
  double dRetLength = 0;
  double adVec[2];
  double a2dM[2][2];

//	before, changed for testing
//	a2dM[0][0] = adMet[0];
//	a2dM[1][1] = adMet[2];

// swapped for tests... [1][1] with [0][0].... im not sure why it changes during the full aniso steps... but this was easier for testing
  a2dM[0][0] = adMet[0];
  a2dM[1][1] = adMet[2];
  a2dM[0][1] = adMet[1];
  a2dM[1][0] = adMet[1];

//	 printf("(%.3e,%.3e) and (%.3e,%.3e)\n",adLoc0[0],adLoc0[1],adLoc1[0],adLoc1[1]);
//	 printf("with metric [%.3e,%.3e,%.3e]\n",adMet[0],adMet[1],adMet[2]);

  for (i = 0; i < 2; i++)
    adVec[i] = adLoc1[i] - adLoc0[i];

  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      dRetLength += a2dM[i][j] * adVec[i] * adVec[j];
    }
  }

  assert(dRetLength >= 0);
  dRetLength = sqrt(dRetLength);

//	    printf("      are %.4f apart\n",dRetLength);

  return dRetLength;
}

bool
AnisoRefinement::qAniEncroachFace(BFace* pBF, double adLoc[2])
{
  double dLenAB, dLenAC, dLenBC;
  double adCv[2];
  double adMet[iNMetricValues];
  int i;
  double adC1[2], adC2[2];
  double adCurVec[2];
  double adNewLoc[2];

  //Determine the Average metric of the two bdry vertices
  for (i = 0; i < iNMetricValues; i++) {
    adMet[i] = pBF->getVert(0)->getMetric(i) + pBF->getVert(1)->getMetric(i);
    adMet[i] *= 0.5;
  }
  //Calculate the offset that should account for curvature

  //First the projection ratio of Vec1 on Vec0
  //   if(false){
  //   double adBisec[2];
  //   double adRealMid[2];
  //   double dRatio;
  //   double dMaxTheta;
  //   double adVec0[2];
  //   double adVec1[2];
  //   double adPerpVec[2];
  //     for(i=0;i<2;i++){
  //       adVec0[i]=pBF->pVVert(1)->adCoords()[i] - pBF->pVVert(0)->adCoords()[i];
  //       adVec1[i]=adLoc[i]-pBF->pVVert(0)->adCoords()[i];
  //     }
  //     dRatio=dDOT2D(adVec0,adVec1);
  //     dRatio/=adVec0[0]*adVec0[0]+adVec0[1]*adVec0[1];
  //     // printf("dRatio=%f\n",dRatio);
  //     if(dRatio<0 || dRatio>1) return false;

  //     adRealMid[0]=pBF->pVVert(0)->dX()+dRatio*(-(pBF->pVVert(0)->dX())+pBF->pVVert(1)->dX());
  //     adRealMid[0]=pBF->pVVert(0)->dY()+dRatio*(-(pBF->pVVert(0)->dY())+pBF->pVVert(1)->dY());

  //     for(i=0;i<2;i++){
  //       adPerpVec[i]=adLoc[i]-adRealMid[i];
  //     }

  //     double dRange=0.1;
  //     double iNumInt=20;
  //     double dTheta=0;
  //     double dBestRatio;
  //     double dTestRatio;
  //     dMaxTheta=0;
  //     dBestRatio=dRatio;
  //     //printf("\ndRatio= %.4f\n",dRatio);
  //     int ii;
  //     for(ii=0;ii<=iNumInt;ii++){
  //       dTestRatio=(dRatio-dRange)+ii*(2.0*dRange)/(1.0*iNumInt);
  //       if(dTestRatio>=0 && dTestRatio<=1){
  // 	dynamic_cast<BdryEdge *>(pBF)->vGetRatioPoint(dTestRatio,adBisec);
  // 	adCurVec[0]=adBisec[0]-adRealMid[0];
  // 	adCurVec[1]=adBisec[1]-adRealMid[1];

  // 	dTheta=dDOT2D(adPerpVec,adCurVec);
  // 	dTheta/=sqrt(adPerpVec[0]*adPerpVec[0]+adPerpVec[1]*adPerpVec[1]);
  // 	dTheta/=sqrt(adCurVec[0]*adCurVec[0]+adCurVec[1]*adCurVec[1]);
  // 	// printf("dTheta = %.4f\n",dTheta);
  // 	if (fabs(dTheta)>fabs(dMaxTheta)){
  // 	  dMaxTheta=dTheta;
  // 	  dBestRatio=dTestRatio;
  // 	}
  //       }
  //     }
  //     //   if(!(fabs(dMaxTheta)>0.9)){
  //     //     printf("SHIT = %.4f\n",dMaxTheta);
  //     //     printf("dBestRatio= %.4f\n",dBestRatio);
  //     //     //  assert(0);
  //     //   }
  //     //  printf("dMaxTheta = %.4f\n",dMaxTheta);

  //     dynamic_cast<BdryEdge *>(pBF)->vGetRatioPoint(dBestRatio,adBisec);
  //     adCurVec[0]=adBisec[0]-adRealMid[0];
  //     adCurVec[1]=adBisec[1]-adRealMid[1];

  //     printf("adCurVec = (%.3e,%.3e)\n",adCurVec[0],adCurVec[1]);
  //   }

  adCurVec[0] = 0;
  adCurVec[1] = 0;

  //Reverse that offset for the location
  adNewLoc[0] = adLoc[0] - adCurVec[0];
  adNewLoc[1] = adLoc[1] - adCurVec[1];

  dLenAB = dMetricLength(pBF->getVert(0)->getCoords(),
			 pBF->getVert(1)->getCoords(), adMet);
  dLenAC = dMetricLength(pBF->getVert(0)->getCoords(), adNewLoc, adMet);
  dLenBC = dMetricLength(adNewLoc, pBF->getVert(1)->getCoords(), adMet);

  adCv[0] = dLenAB * dLenAB - dLenBC * dLenBC + dLenAC * dLenAC;

  adCv[0] /= (2.0 * dLenAB);
  adCv[1] = dLenAC * dLenAC - adCv[0] * adCv[0];
  adCv[1] = sqrt(adCv[1]);

  //Now we have adLoc's virtual location, we can check encroachment easily

  adC1[0] = 0.5 * dLenAB;
  adC1[1] = 0.5 / sqrt(3.0) * dLenAB; //+dShift;

  adC2[0] = 0.5 * dLenAB;
  adC2[1] = -0.5 / sqrt(3.0) * dLenAB; //+dShift;

  double dDist1 = dDIST2D(adCv, adC1);
  double dDist2 = dDIST2D(adCv, adC2);

  //  int iFC1 = iFuzzyComp(dDist1, dLenAB/sqrt(3));
  // int iFC2 = iFuzzyComp(dDist2, dLenAB/sqrt(3));
  double dNewRad = sqrt(adC2[1] * adC2[1] + adC2[0] * adC2[0]);
  int iFC1 = iFuzzyComp(dDist1, dNewRad);
  int iFC2 = iFuzzyComp(dDist2, dNewRad);

  if (iFC1 < 0 && iFC2 < 0) {
    //    assert(0);
    return true; //The boundary is encroached
  }
  else {
    //boundary not encroached, check gauss points
    return false;
    return (!qIsCurvedBdryOkay(pBF->getFace(0), adLoc)); //test for GaussPoints
  }
}

//Figure out if any of the input vertices encroaches on a boundary in the anisotropic sense. This is probably inefficient, but it's a small part of the overall process.
void
AnisoRefinement::vHealBdry()
{
  bool qEncroachmentFound;
  GR_index_t iBF, i, ii;
  double adTestLoc[2];
  BFace* pBdryFace = NULL;
  std::set<Cell*> spCInc; //just to make vNeighborhood work
  std::set<Vert*> spVNeigh; //
  std::set<Vert*> spVPossEncVerts;
  std::set<Vert*> spVTemp;
  std::set<Vert*>::iterator iterV, iterV2;
  Vert* pVertStart;
  printf("Healing %d bdry faces\n", pM2D->getNumBdryFaces());

  for (iBF = 0; iBF < pM2D->getNumBdryFaces(); iBF++) {
    spVPossEncVerts.clear();
    pBdryFace = pM2D->getBFace(iBF);
    if (pBdryFace->isDeleted())
      continue;
    //first get the opposite vertex of the boundary face
    pVertStart = pBdryFace->getFace(0)->getCell(0)->getOppositeVert(
	pBdryFace->getFace(0));
    spVPossEncVerts.insert(pVertStart);

    //then for a given number of iterations,
    for (i = 0; i < 2; i++) {
      spVTemp.clear();
      //find the neighbors of each poss encroaching vertex
      for (ii = 0; ii < 2; ii++) {
	for (iterV = spVPossEncVerts.begin(); iterV != spVPossEncVerts.end();
	    ++iterV) {
	  spVNeigh.clear();
	  findNeighborhoodInfo(*iterV, spCInc, spVNeigh);
	  for (iterV2 = spVNeigh.begin(); iterV2 != spVNeigh.end(); iterV2++) {
	    //if they are not boundary verts, insert them.
	    if (!(*iterV2)->isBdryVert())
	      spVTemp.insert(*iterV2);
	  }
	}

      }
      //copy the temp set into the possenc set
      for (iterV = spVTemp.begin(); iterV != spVTemp.end(); ++iterV) {
	spVPossEncVerts.insert(*iterV);
      }
    }
    //Now you have some possibly encroaching verts, check if they should be removed

    for (iterV = spVPossEncVerts.begin(); iterV != spVPossEncVerts.end();
	++iterV) {
      qEncroachmentFound = false;
      //iterate through vertices
      adTestLoc[0] = ((*iterV)->x());
      adTestLoc[1] = ((*iterV)->y());
      // printf("Testing (%.3e,%.3e)\n",adTestLoc[0],adTestLoc[1]);
      if (qAniEncroachFace(pBdryFace, adTestLoc)) {
	qEncroachmentFound = true;
	//  printf("Boundary Encroached\n");
	//pFEncroached=pM2D->pBFBFace(iBF)->pFFace(0);
	//assert(pFEncroached->qIsBdryFace());
	//pCEncroached=pFEncroached->pCCell(0);
	//	Vert*  pVNew=pM2D->pVNewVert();
	//	  pM2D->iInsertOnBdryFace (pVNew,pFEncroached,pCEncroached,false,false);
	int iTest;
	bool qRemoved;
	// printf("removing Vert\n");
	qRemoved = pM2D->removeVertCheaply(*iterV, iTest);
	//qRemoved=pM2D->qRemoveVert(pM2D->pVVert(iV),iTest);
	if (qRemoved) {
	  // printf("Gone\n");
	  // qEncroachmendFound=true;
	  // setpVVerts.erase(iter);
	}
	else {
	  printf("Can't remove Vert (%.3e,%.3e) inserting instead\n",
		 (*iterV)->x(), (*iterV)->y());
	  Face * pFEncroached = pM2D->getBFace(iBF)->getFace(0);
	  double newPointLoc[2];
	  pM2D->getBFace(iBF)->calcCentroid(newPointLoc);
	  pSI2D->insertPointOnBdryFace(newPointLoc, pFEncroached);
	  break;
	}
      }
    }
    if (qEncroachmentFound) {
      iBF--;
      break;
    }
  }

  pM2D->purgeAllEntities();
  return;
}

Cell*
AnisoRefinement::pCAnisoCircumCenter(Cell *pC, double adCircCen[2])
{
  //create the virtual mesh around the cell
  VirtualTransMesh2D VTM2D(this, pC);

  //if the virtual mesh was made without running in to a boundary
  if (VTM2D.pEndOfMeshBdryFace == NULL) {
    //    	  if(!pC->isCircumcircleEmpty())
    //        if((!pC->isCircumcircleEmpty() && (VTM2D.vecpCCells.size()>1)) || VTM2D.vecpCCells.size()>=MAXCELLS ){
    //          //the new location may be too close to another vertex so special procedures are required
    //          VTM2D.vAdaptToNonEmptyCC();
    //
    //        }
    VTM2D.vUnTransformDest(adCircCen);
    return VTM2D.vecpCCells.back();
  }

  // TODO: Is the following code correct?  Can it even ever be executed for
  // possible input meshes?

  if (VTM2D.pEndOfMeshBdryFace->isBdryFace()) {
    //the virtual mesh hit a boundary so that face needs to be split
    adCircCen[0] = LARGE_DBL - 1;
    adCircCen[1] = LARGE_DBL - 1;
  }

  //  printf("need to output the boundary face associate with the last face used in contruction\n");
  //  assert(VTM2D.pEndOfMeshBdryFace->qIsBdryFace());
  Cell * cell =
      VTM2D.pEndOfMeshBdryFace->getCell(0)->getType() == Cell::eTriCell ?
	  VTM2D.pEndOfMeshBdryFace->getCell(1) :
	  VTM2D.pEndOfMeshBdryFace->getCell(0);
  assert(cell->getType() != Cell::eTriCell);

  return cell;
}

int
AnisoRefinement::iRefine(CellQueue* pCellQueue, int iMaxVertsAdded = 10)
{
  Cell* pCelltoFix;
  Cell* pCellforInsert;
  double adCircCen[2];
  int iNumVertsAdded = 0;
  int iBdryEncroached;
  int iCounter = 0;
  int iWhentoPrint = 100;
  int i;
  // assert2(pM2D->iSwap()==0,"FAILED:Mesh wasn't properly prepared\n");

  Vert* apVOfCell[3];

  while (1) {
    pCelltoFix = pCellQueue->pCTopCell();
    if (pCelltoFix == NULL)
      break;
    for (i = 0; i < 3; i++) {
      apVOfCell[i] = pCelltoFix->getVert(i);
    }

    //   printf("\nCell to refine is (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",pCelltoFix->pVVert(0)->dX(),pCelltoFix->pVVert(0)->dY(),pCelltoFix->pVVert(1)->dX(),pCelltoFix->pVVert(1)->dY(),pCelltoFix->pVVert(2)->dX(),pCelltoFix->pVVert(2)->dY());
    //  printf("CellQueue qHasCell returns %d (%d=true)\n",pCellQueue->qHasCell(pCelltoFix),true);
    //		pCelltoFix->printCellInfo();
    pCellforInsert = pCAnisoCircumCenter(pCelltoFix, adCircCen);
    //		printf("Circumcenter is (%.3e,%.3e)\n",adCircCen[0],adCircCen[1]);

    if (pCellforInsert->isBdryCell()) {
      //adCircCen[0]=0.5*(pCellforInsert->pVVert(0)->dX()+pCellforInsert->pVVert(1)->dX());
      //adCircCen[1]=0.5*(pCellforInsert->pVVert(0)->dY()+pCellforInsert->pVVert(1)->dY());
      //  printf("(%.3e,%.3e) on face (%.3e,%.3e)-(%.3e,%.3e)\n",pVNewAniVert->dX(),pVNewAniVert->dY(),pCellforInsert->pVVert(0)->dX(),pCellforInsert->pVVert(0)->dY(),pCellforInsert->pVVert(1)->dY(),pCellforInsert->pVVert(1)->dY());
      pSI2D->insertPointOnBdryFace(adCircCen, pCellforInsert->getFace(0));
      iNumVertsAdded++;
      iCounter++;
      //     printf("Inserted on the boundary\n")
    }
    else {
      //       printf("pCellforInsert is (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",pCellforInsert->pVVert(0)->dX(),pCellforInsert->pVVert(0)->dY(),pCellforInsert->pVVert(1)->dX(),pCellforInsert->pVVert(1)->dY(),pCellforInsert->pVVert(2)->dX(),pCellforInsert->pVVert(2)->dY());
      //           printf("Need to test for encroachment\n");
      // iBdryEncroached=iCheckBoundary(adCircCen);
      iBdryEncroached = iCheckBoundary(adCircCen, pCellforInsert);
      iNumVertsAdded += iBdryEncroached;
      iCounter += iBdryEncroached;
      //          printf("iCheckBoundary returns %d\n",iBdryEncroached);
      //    if(iBdryEncroached!=0) printf("inserted %d vertices on bdry\n",iBdryEncroached);

      if (iBdryEncroached == 0) {
	//             printf("No Encroachment, inserting at (%.3e,%.3e)\n",adCircCen[0],adCircCen[1]);
	//       printf("inserting in interior\n");
	// Find the actual cell containing the new point and insert it.
	pSI2D->insertPoint(adCircCen, pCellforInsert);
	// printf("inserted in interior\n");
	//      int itest;

	iNumVertsAdded++;
	iCounter++;
      }
    }
    bool qTest = true;
    for (i = 0; i < 3; i++) {
      qTest &= pCelltoFix->hasVert(apVOfCell[i]);
    }
    if (qTest) {
      printf("iFixFace needed\n");
      //assert2(0,"Didn't change cell\n");
      iFixFace(pCelltoFix);
    }

    //   else{
    //     printf("Interior, No problem \n");
    //      printf("Inserting at (%.3e,%.3e)\n",adCircCen[0],adCircCen[1]);
    //     //set the new vertex to the circumcenter coordinates
    //     pVNewAniVert->vSetCoords(2,adCircCen);
    //     pM2D->iInsertInInterior(pVNewAniVert,pCellforInsert,true);
    //     }

    //  printf("%d vertices added\n",iNumVertsAdded);

    // pCellQueue->vPopBiggest();
    if (iCounter >= iWhentoPrint) {
      printf("Inserted %d Vertices so far: Big Cells: %zd\n", iNumVertsAdded,
	     pCellQueue->iBig());

      iCounter = 0;
    }

    //    assert(iRepairMesh()==0);
//		writeVTK(*pM2D, "After_DuringInsertion");
    if (iNumVertsAdded >= iMaxVertsAdded) {
      printf("\nNot allowed to insert anymore\n\n");
      break;
    }
  }

  printf("Added %d Vertices in total\n\n", iNumVertsAdded);
  return iNumVertsAdded;
}

//This function takes a location for insertion. It then checks all boundary faces for encroachment. Any encroached boundary faces are split and the return value is the number of splits done (Inefficient but definitely works)
int
AnisoRefinement::iCheckBoundary(double adLoc[2])
{
  GR_index_t i, iRetVal;
  Face* pFEncroached;
  double adMidPoint[2];

  iRetVal = 0;
  for (i = 0; i < pM2D->getNumBdryFaces(); i++) {
    if (qAniEncroachFace(pM2D->getBFace(i), adLoc)) {
      pFEncroached = pM2D->getBFace(i)->getFace(0);
      vAniMidPoint(pFEncroached, adMidPoint);
      //  printf("Boundary split at (%.3e,%.3e)\n",pVNew->dX(),pVNew->dY());
      assert(pFEncroached->isBdryFace());
      //  qAniLengthsFixed=!qRepairingMesh; //this will determine whether or not inserton bdry face checks the new cells for validity
      pSI2D->insertPointOnBdryFace(adMidPoint, pFEncroached);

      //  printf("Split Edge with (%.3e,%.3e)\n",pVNew->dX(),pVNew->dY());
      iRetVal++;
      //  break;
    }
  }

  return iRetVal;
}

//This function takes a cell, and finds all the boundary faces around it. It then checks encroachment of the location ->This has no guarantees of working, but is quick.
int
AnisoRefinement::iCheckBoundary(double adLoc[2], Cell* pCellwLoc)
{
  // printf("Checking boundary\n");
  //  int iRecur=4;
  int iRetVal;
  Face* pFEncroached;
  double adMidPoint[2];

  std::set<BFace*> setpBF;
  std::set<BFace*>::iterator iter;
  Cell* pCellNeigh_i;
  Cell* pCellNeigh_ii;
  Cell* pCellNeigh_iii;

  assert(!pCellwLoc->isBdryCell());
  // printf("Building Set of Boundary Faces\n");
  int i, ii, iii;
  for (i = 0; i < pCellwLoc->getNumFaces(); i++) {
    //  printf("i=%d ",i);
    pCellNeigh_i = pCellwLoc->getFace(i)->getOppositeCell(pCellwLoc);
    if (pCellNeigh_i->isBdryCell())
      setpBF.insert(dynamic_cast<BFace*>(pCellNeigh_i));
    else {
      for (ii = 0; ii < 3; ii++) {
	// printf("ii=%d ",ii);
	pCellNeigh_ii = pCellNeigh_i->getFace(i)->getOppositeCell(pCellNeigh_i);
	if (pCellNeigh_ii->isBdryCell())
	  setpBF.insert(dynamic_cast<BFace*>(pCellNeigh_ii));
	else {
	  for (iii = 0; iii < 3; iii++) {
	    // printf("iii=%d ",iii);
	    pCellNeigh_iii = pCellNeigh_i->getFace(i)->getOppositeCell(
		pCellNeigh_i);
	    if (pCellNeigh_iii->isBdryCell())
	      setpBF.insert(dynamic_cast<BFace*>(pCellNeigh_iii));
	  }
	}
      }
    }
  }
  // printf("\n");
  //  printf("Searching %d boundaries\n",setpBF.size());

  iRetVal = 0;
  for (iter = setpBF.begin(); iter != setpBF.end(); ++iter) {
    if (qAniEncroachFace(*iter, adLoc)) {
      pFEncroached = (*iter)->getFace(0);
      adMidPoint[0] = 0.5
	  * (pFEncroached->getVert(0)->x() + pFEncroached->getVert(1)->x());
      adMidPoint[1] = 0.5
	  * (pFEncroached->getVert(0)->y() + pFEncroached->getVert(1)->y());
      assert(pFEncroached->isBdryFace());
      //    printf("Split Edge at (%.3e,%.3e)\n",adMidPoint[0],adMidPoint[1]);
      pSI2D->insertPointOnBdryFace(adMidPoint, pFEncroached);

      iRetVal++;
      //  if(iRetVal>1) printf("Multiple Split\n");
      // break;
    }
  }
  return iRetVal;
}

VirtualTransMesh2D::VirtualTransMesh2D(AnisoRefinement* pAnisoRef_,
				       Cell* pCStart) :
    vecpCCells(), vecvvVerts(), pAnisoRef(pAnisoRef_), pEndOfMeshBdryFace(
    NULL)
{
  virtualvert vvNewVert[3];
  double dFopp[3];
  int i;
  assert2(!pCStart->isBdryCell(), "Bdry Cell in TransMesh??\n");

  vecpCCells.push_back(pCStart);
  ///the lengths of edges
  for (i = 0; i < 3; i++) {
    vvNewVert[i].pVRealVert = pCStart->getVert(i);
  }

  double adMetABC[pAnisoRef->iNMetricValues];
  for (i = 0; i < pAnisoRef->iNMetricValues; i++) {
    //   printf("Metrics %d->%.3e, %3e, %3e\n",i,pCStart->getVert(0)->dGetMetric(i),pCStart->getVert(1)->dGetMetric(i),pCStart->getVert(2)->dGetMetric(i));
    adMetABC[i] = pCStart->getVert(0)->getMetric(i)
	+ pCStart->getVert(1)->getMetric(i) + pCStart->getVert(2)->getMetric(i);
    adMetABC[i] /= 3.0;
  }
  //  printf("Metric is [%.3e,%.3e,%.3e]\n",adMetABC[0],adMetABC[1],adMetABC[2]);

  dFopp[0] = pAnisoRef->dMetricLength(pCStart->getVert(1)->getCoords(),
				      pCStart->getVert(2)->getCoords(),
				      adMetABC);
  dFopp[1] = pAnisoRef->dMetricLength(pCStart->getVert(2)->getCoords(),
				      pCStart->getVert(0)->getCoords(),
				      adMetABC);
  dFopp[2] = pAnisoRef->dMetricLength(pCStart->getVert(0)->getCoords(),
				      pCStart->getVert(1)->getCoords(),
				      adMetABC);

  //These points are chosen so as to satisfy the length requirements
  //reference is 0,0
  // printf("Lengths are [BC%.3e,CA%.3e,AB%.3e]\n",dFopp[0],dFopp[1],dFopp[2]);
  vvNewVert[0].dvX = 0;
  vvNewVert[0].dvY = 0;
  //Next point is at (dFopp2,0)
  vvNewVert[1].dvX = dFopp[2];
  vvNewVert[1].dvY = 0;
  //Third point is at (more complicated)
  for (i = 0; i < 3; i++) {
    dFopp[i] *= dFopp[i];
  }
  vvNewVert[2].dvX = 0.5 * (dFopp[2] + dFopp[1] - dFopp[0]) / sqrt(dFopp[2]);
  vvNewVert[2].dvY = 2 * sqrt(dFopp[2]) * vvNewVert[2].dvX;
  vvNewVert[2].dvY += dFopp[0] - dFopp[2];
  vvNewVert[2].dvY -= vvNewVert[2].dvX * vvNewVert[2].dvX;
  vvNewVert[2].dvY = sqrt(vvNewVert[2].dvY);
  //  printf("Locations are (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",vvNewVert[0].dvX,vvNewVert[0].dvY,vvNewVert[1].dvX,vvNewVert[1].dvY,vvNewVert[2].dvX,vvNewVert[2].dvY);
  //load the vertices in order
  for (i = 0; i < 3; i++) {
    vecvvVerts.push_back(vvNewVert[i]);
  }

  //  printf("Locations used to be (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",vvNewVert[0].pVRealVert->dX(),vvNewVert[0].pVRealVert->dY(),vvNewVert[1].pVRealVert->dX(),vvNewVert[1].pVRealVert->dY(),vvNewVert[2].pVRealVert->dX(),vvNewVert[2].pVRealVert->dY());

  //The start of this mesh is the centroid of the cell because it's guaranteed to be interior
  adStart[0] = (vvNewVert[0].dvX + vvNewVert[1].dvX + vvNewVert[2].dvX) / 3;
  adStart[1] = (vvNewVert[0].dvY + vvNewVert[1].dvY + vvNewVert[2].dvY) / 3;

  //  printf("Virtual Centroid is (%.3e,%.3e)\n",adStart[0],adStart[1]);
  //Then we compute the adDest of the mesh
  double adA[2], adB[2], adC[2];

  adA[0] = vvNewVert[0].dvX;
  adA[1] = vvNewVert[0].dvY;

  adB[0] = vvNewVert[1].dvX;
  adB[1] = vvNewVert[1].dvY;

  adC[0] = vvNewVert[2].dvX;
  adC[1] = vvNewVert[2].dvY;
  vTriCircumcenter(adA, adB, adC, adDest);
  //  printf("Virtual Circumcenter is (%.3e,%.3e)\n",adDest[0],adDest[1]);

  vBuildMeshToDest();
  // if(vecpCCells.size()>2)  printf(" #=%d ",vecpCCells.size());
}

void
VirtualTransMesh2D::vAdaptToNonEmptyCC()
{

  pEndOfMeshBdryFace = NULL;
  int iF;
  std::vector<Cell *>::iterator iterCells;
  std::vector<virtualvert>::iterator iterVerts;
  iterCells = vecpCCells.begin();
  Cell* pCFirst = *iterCells;
  iterCells++;
  Cell* pCSecond = *iterCells;
  Face* pFShared = NULL;
  virtualvert vvVert0;
  virtualvert vvVert1;

  // printf("CC was (%.3e,%.3e)\n",adDest[0],adDest[1]);

  //first, put the adCircCen in the middle of the first edge
  for (iF = 0; iF < pCFirst->getNumFaces(); iF++) {
    if (pCSecond->hasFace(pCFirst->getFace(iF))) {
      pFShared = pCFirst->getFace(iF);
      break;
    }
  }
  assert(pFShared!=NULL);
  vvVert0 = vvFindVirtualVert(pFShared->getVert(0));
  vvVert1 = vvFindVirtualVert(pFShared->getVert(1));
  adDest[0] = (vvVert0.dvX + vvVert1.dvX) * 0.5;
  adDest[1] = (vvVert0.dvY + vvVert1.dvY) * 0.5;
  //then , call the smooth procedure with all the vertex locations from the virtual mesh
  //  printf("Starting with point %.3e,%.3e\n",adDest[0],adDest[1]);
  int iMaxVertNeigh = 4;
  double **a2dNeighPts = new double *[iMaxVertNeigh];
  int **a2iFaceConn = new int *[iMaxVertNeigh];
  int i;
  int iNeigh = int(vecvvVerts.size());
  //printf("iNeigh = %d\n",iNeigh);

  for (i = 0; i < iMaxVertNeigh; i++) {
    a2dNeighPts[i] = new double[2];
    a2iFaceConn[i] = new int[2];
  }

  i = 0;
  int iNV = int(vecvvVerts.size());
  for (iterVerts = vecvvVerts.begin(); iterVerts != vecvvVerts.end();
      iterVerts++) {
    if (i == iMaxVertNeigh) {
      iNV = iMaxVertNeigh;
      iNeigh = iMaxVertNeigh;
      //   assert(0);
      break;
    }
    a2dNeighPts[i][0] = iterVerts->dvX;
    a2dNeighPts[i][1] = iterVerts->dvY;
    i++;
    //    printf("(%.3e,%.3e) theta=%.3e\n",a2dNeighPts[i][0],a2dNeighPts[i][1],atan2(a2dNeighPts[i][1]-adDest[1],a2dNeighPts[i][0]-adDest[0]));
  }

  //Simple bubblesort to order a2dNeighPts by angle
  double dX, dY;
  double dAng, dAng2;

  //  assert(iNV>1);
  bool qSorted = false;
  virtualvert vvVert2;
  while (qSorted == false) {
    qSorted = true;
    for (i = 0; i < iNV - 1; i++) {
      dX = a2dNeighPts[i][0] - adDest[0];
      dY = a2dNeighPts[i][1] - adDest[1];
      dAng = atan2(dY, dX);
      dAng2 = atan2(a2dNeighPts[i + 1][1] - adDest[1],
		    a2dNeighPts[i + 1][0] - adDest[0]);
      //   printf("dAng(i)=%.3e\n",dAng);
      //    printf("dAng(i+1)=%.3e\n",dAng2);
      if (dAng < dAng2) {
	//      printf("switch\n");
	qSorted = false;
	a2dNeighPts[i][0] = a2dNeighPts[i + 1][0];
	a2dNeighPts[i][1] = a2dNeighPts[i + 1][1];
	a2dNeighPts[i + 1][0] = dX + adDest[0];
	a2dNeighPts[i + 1][1] = dY + adDest[1];
      }
    }
    //   printf("\n");
    // 	    for(int ii=0;ii<iNV;ii++){
    // 	      printf("(%.3e,%.3e) theta=%.3e\n",a2dNeighPts[ii][0],a2dNeighPts[ii][1],atan2(a2dNeighPts[ii][1]-adDest[1],a2dNeighPts[ii][0]-adDest[0]));
    // 	    }
  }
  //    printf("\nfinal:\n");
  //       for(i=0;i<iNV;i++){
  // 	printf("(%.3e,%.3e) theta=%.3e\n",a2dNeighPts[i][0],a2dNeighPts[i][1],atan2(a2dNeighPts[i][1]-adDest[1],a2dNeighPts[i][0]-adDest[0]));
  //       }

  a2iFaceConn[0][0] = 0;
  a2iFaceConn[0][1] = iNeigh - 1;
  int ii;
  for (ii = 1; ii < iNeigh; ii++) {
    a2iFaceConn[ii][0] = ii;
    a2iFaceConn[ii][1] = ii - 1;
  }

  SMinitSmoothStats(pAnisoRef->pvSmoothData);
  // SMuntangle (iNeigh, iNeigh, adDest, a2dNeighPts, a2iFaceConn,pAnisoRef->pM2D->pvSmoothData);  //Currently not available

  SMsmooth(iNeigh, iNeigh, adDest, a2dNeighPts, a2iFaceConn,
	   pAnisoRef->pvSmoothData, 0);
  //  printf("Smoothed Location is (%.3e,%.3e)\n",adDest[0],adDest[1]);

  for (i = 0; i < iMaxVertNeigh; i++) {
    delete[] a2dNeighPts[i];
    delete[] a2iFaceConn[i];
  }
  delete[] a2dNeighPts;
  delete[] a2iFaceConn;

  //erase the mesh until the new location is within the last cell
  double adV0[2];
  double adV1[2];
  double adV2[2];
  bool qLastCell = false;
  while (!qLastCell) {
    iterCells = vecpCCells.end();
    iterCells--;
    Cell* pCLast = *iterCells;
    vvVert0 = vvFindVirtualVert(pCLast->getVert(0));
    vvVert1 = vvFindVirtualVert(pCLast->getVert(1));
    vvVert2 = vvFindVirtualVert(pCLast->getVert(2));

    adV0[0] = vvVert0.dvX;
    adV0[1] = vvVert0.dvY;
    adV1[0] = vvVert1.dvX;
    adV1[1] = vvVert1.dvY;
    adV2[0] = vvVert2.dvX;
    adV2[1] = vvVert2.dvY;

    if (isInTri(adDest, adV0, adV1, adV2) >= 2) {
      qLastCell = true;
    }
    else
      vecpCCells.erase(iterCells);
  }
  //  printf("Found the last cell\n");

  //the resulting smooth location is where the point will be inserted
  //    assert(0);
}

virtualvert
VirtualTransMesh2D::vvFindVirtualVert(Vert* pVReal)
{
  std::vector<virtualvert>::iterator iter;
  virtualvert vvRetVirtualVert;
  bool qFound = false;
  // printf("Searching through %d verts\n",vecvvVerts.size());
  for (iter = vecvvVerts.begin(); iter != vecvvVerts.end(); ++iter) {
    if ((*iter).pVRealVert == pVReal) {
      vvRetVirtualVert = *iter;
      qFound = true;
    }
  }
  assert(qFound);
  return vvRetVirtualVert;
}

void
VirtualTransMesh2D::vAddCell(Face* pFforAdd)
{
  //This fuction adds the appropriate cell to the vecpCCells
  //and computes the imaginary location of the new vertex
  virtualvert vvVert0, vvVert1; //the two locations already known
  virtualvert vvVert2; //the new vertex
  double dLen[3];
  double dTheta;
  bool qLeftOfFace = false;

  //the new cell
  Cell* pCNewCell;
  pCNewCell = pFforAdd->getOppositeCell(vecpCCells.back());
  //   printf("Adding Cell(%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n",pCNewCell->getVert(0)->dX(),pCNewCell->getVert(0)->dY(),pCNewCell->getVert(1)->dX(),pCNewCell->getVert(1)->dY(),pCNewCell->getVert(2)->dX(),pCNewCell->getVert(2)->dY());

  vecpCCells.push_back(pCNewCell);
  qLeftOfFace = (pCNewCell == pFforAdd->getLeftCell());

  double adMetNew[pAnisoRef->iNMetricValues];
  int i;
  for (i = 0; i < pAnisoRef->iNMetricValues; i++) {
    adMetNew[i] = pCNewCell->getVert(0)->getMetric(i)
	+ pCNewCell->getVert(1)->getMetric(i)
	+ pCNewCell->getVert(2)->getMetric(i);
    adMetNew[i] /= 3.0;
  }
  // printf("New cell metric =[%.3e,%.3e,%.3e]\n",adMetNew[0],adMetNew[1],adMetNew[2]);

  Face* pF;
  pF = pCNewCell->getOppositeFace(pFforAdd->getVert(0));
  dLen[0] = pAnisoRef->dMetricLength(pF->getVert(0)->getCoords(),
				     pF->getVert(1)->getCoords(), adMetNew);

  pF = pCNewCell->getOppositeFace(pFforAdd->getVert(1));
  dLen[1] = pAnisoRef->dMetricLength(pF->getVert(0)->getCoords(),
				     pF->getVert(1)->getCoords(), adMetNew);

  pF = pFforAdd;
  dLen[2] = pAnisoRef->dMetricLength(pF->getVert(0)->getCoords(),
				     pF->getVert(1)->getCoords(), adMetNew);

  vvVert0 = vvFindVirtualVert(pFforAdd->getVert(0));
  vvVert1 = vvFindVirtualVert(pFforAdd->getVert(1));
  vvVert2.pVRealVert = pCNewCell->getOppositeVert(pFforAdd);

  //  printf("New cell lengths =%.3e, %.3e, %.3e\n",dLen[0],dLen[1],dLen[2]);

  double dOldLen = (vvVert0.dvX - vvVert1.dvX) * (vvVert0.dvX - vvVert1.dvX);
  dOldLen += (vvVert0.dvY - vvVert1.dvY) * (vvVert0.dvY - vvVert1.dvY);
  dOldLen = sqrt(dOldLen);

  double dScale = 1.0 * dOldLen / dLen[0];
  //  printf("Scale = %.3e\n",dScale);
  dLen[0] *= dScale;
  dLen[1] *= dScale;
  dLen[2] *= dScale;
  //the pointer for the new vertex
  //     printf("Converted to =%.3e, %.3e, %.3e\n",dLen[0],dLen[1],dLen[2]);

  //  printf("The old vertices are at (%.3e,%.3e) and (%.3e,%.3e)\n",vvVert0.dvX,vvVert0.dvY,vvVert1.dvX,vvVert1.dvY);
  //  printf("The edge lengths are %.3e,%.3e,%.3e\n",dLen[0],dLen[1],dLen[2]);
  //the computation

  dTheta = 0.5 * (dLen[1] * dLen[1] + dLen[2] * dLen[2] - dLen[0] * dLen[0])
      / dLen[1] / dLen[2];
  dTheta = acos(dTheta);
  if (qLeftOfFace)
    dTheta = -dTheta; //depends on which side

  vvVert2.dvX = cos(dTheta) * (vvVert1.dvX - vvVert0.dvX)
      + sin(dTheta) * (vvVert1.dvY - vvVert0.dvY);
  vvVert2.dvX = vvVert0.dvX + dLen[1] / dLen[2] * vvVert2.dvX;

  vvVert2.dvY = -sin(dTheta) * (vvVert1.dvX - vvVert0.dvX)
      + cos(dTheta) * (vvVert1.dvY - vvVert0.dvY);
  vvVert2.dvY = vvVert0.dvY + dLen[1] / dLen[2] * vvVert2.dvY;

  //  printf("The new Vertex is at (%.3e,%.3e)\n",vvVert2.dvX,vvVert2.dvY);

  vecvvVerts.push_back(vvVert2);

  return;
}

bool
VirtualTransMesh2D::qDoesLineIntersectFace(double adP1[2], double adP2[2],
					   Face* pF)
//return 1 if yes
//return -1 if no
//return 0 if it is on the line

{
  double adP3[2], adP4[2];
  double dTest;

  //  printf("Faces real coordinate (%.3e,%.3e)-(%.3e,%.3e)\n",pF->pVVert(0)->dX(),pF->pVVert(0)->dY(),pF->pVVert(1)->dX(),pF->pVVert(1)->dY());
  adP3[0] = vvFindVirtualVert(pF->getVert(0)).dvX;
  adP3[1] = vvFindVirtualVert(pF->getVert(0)).dvY;

  adP4[0] = vvFindVirtualVert(pF->getVert(1)).dvX;
  adP4[1] = vvFindVirtualVert(pF->getVert(1)).dvY;

  //  printf("line (%.3e,%.3e)-(%.3e,%.3e) and (%.3e,%.3e)-(%.3e,%.3e)\n",adP1[0],adP1[1],adP2[0],adP2[1],adP3[0],adP3[1],adP4[0],adP4[1]);

  dTest = (adP4[0] - adP3[0]) * (adP1[1] - adP3[1])
      - (adP4[1] - adP3[1]) * (adP1[0] - adP3[0]);
  dTest /= (adP4[1] - adP3[1]) * (adP2[0] - adP1[0])
      - (adP4[0] - adP3[0]) * (adP2[1] - adP1[1]);
  //  printf("First test give %.6e\n",dTest);
  if (dTest < 0 || dTest > 1)
    return false;

  dTest = (adP2[0] - adP1[0]) * (adP3[1] - adP1[1])
      - (adP2[1] - adP1[1]) * (adP3[0] - adP1[0]);
  dTest /= (adP2[1] - adP1[1]) * (adP4[0] - adP3[0])
      - (adP2[0] - adP1[0]) * (adP4[1] - adP3[1]);
  //  printf("Second test give %.6e\n",dTest);
  if (dTest < 0 || dTest > 1)
    return false;
  else
    return true;
}

void
VirtualTransMesh2D::vBuildMeshToDest()
{

  int i;
  assert(vecpCCells.size() == 1);
  //make sure the mesh hasn't been made yet (this can be changed later
  bool qDone = false;
  int iisintri;

  Face* pF = NULL;
  Face* pFLastFace = NULL;
  Cell* pCLast = NULL;

  double adP1[2];
  double adP2[2];
  double adP3[2];

  while (qDone == false) {
    pCLast = vecpCCells.back();

    adP1[0] = vvFindVirtualVert(pCLast->getVert(0)).dvX;
    adP1[1] = vvFindVirtualVert(pCLast->getVert(0)).dvY;

    adP2[0] = vvFindVirtualVert(pCLast->getVert(1)).dvX;
    adP2[1] = vvFindVirtualVert(pCLast->getVert(1)).dvY;

    adP3[0] = vvFindVirtualVert(pCLast->getVert(2)).dvX;
    adP3[1] = vvFindVirtualVert(pCLast->getVert(2)).dvY;

    iisintri = isInTri(adDest, adP1, adP2, adP3);

    if (iisintri == 2 || iisintri == 3) {
      qDone = true;
      break;
    }

    for (i = 0; i < 3; i++) {
      qDone = true;
      if (pFLastFace != vecpCCells.back()->getFace(i)) { //only check the faces not equal to the last face
	pF = vecpCCells.back()->getFace(i);
	if (qDoesLineIntersectFace(adStart, adDest, pF)) {

	  pFLastFace = pF;
	  //     printf("Found intersection with face (%.3e,%.3e)-(%.3e,%.3e)\n",pF->pVVert(0)->dX(),pF->pVVert(0)->dY(),pF->pVVert(1)->dX(),pF->pVVert(1)->dY());
	  qDone = false;
	  if (pF->isBdryFace()) {
	    //	  printf("Hit bdry\n");
	    pEndOfMeshBdryFace = pF;
	    return;
	  }
	  if (vecpCCells.size() == MAXCELLS) {
	    return;
	  }
	  // 		printf("%d\n",vecpCCells.size());
	  //  		pEndOfMeshBdryFace=pF;
	  // 		vAdaptToNonEmptyCC();

	  //  		return;
	  //  	      }
	  vAddCell(pF);
	  break;
	}
      }
    }
  }
  return;
}

void
VirtualTransMesh2D::vUnTransformLoc(double adVirtLoc[2], double adRealLoc[2])
{
  //Should assert the destination is within the final cell
  int i;
  Cell* pC = vecpCCells.back();
  virtualvert vvVert[3]; //the three virtual vertices
  Vert* pVRealVert[3]; //the three real vertices

  for (i = 0; i < 3; i++) {
    pVRealVert[i] = pC->getVert(i);
    vvVert[i] = vvFindVirtualVert(pVRealVert[i]);
  }
  double R0[2], R1[2], R2[2];
  double V0[2], V1[2], V2[2];
  double t1, t2, t11, t13, t14, t18, t20, t23, t27, t30, t33;

  R0[0] = pVRealVert[0]->x(); //this should all be cleaned up
  R0[1] = pVRealVert[0]->y();

  R1[0] = pVRealVert[1]->x();
  R1[1] = pVRealVert[1]->y();

  R2[0] = pVRealVert[2]->x();
  R2[1] = pVRealVert[2]->y();

  V0[0] = vvVert[0].dvX;
  V0[1] = vvVert[0].dvY;

  V1[0] = vvVert[1].dvX;
  V1[1] = vvVert[1].dvY;

  V2[0] = vvVert[2].dvX;
  V2[1] = vvVert[2].dvY;
  //	printf("The real triangle has vertices (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",R0[0],R0[1],R1[0],R1[1],R2[0],R2[1]);
  //	printf("The virtual triangle has vertices (%.3e,%.3e) (%.3e,%.3e) (%.3e,%.3e)\n",V0[0],V0[1],V1[0],V1[1],V2[0],V2[1]);

//	assert(isInTri(adDest,V0,V1,V2)>1);

  t1 = R2[0] - R0[0];
  t2 = -V1[1] + V0[1];
  t11 = 1
      / (-V2[0] * V1[1] + V2[0] * V0[1] + V0[0] * V1[1] + V1[0] * V2[1]
	  - V1[0] * V0[1] - V0[0] * V2[1]);
  t13 = R1[0] - R0[0];
  t14 = V2[1] - V0[1];
  t18 = adVirtLoc[0] - V0[0];
  t20 = -V1[0] + V0[0];
  t23 = V2[0] - V0[0];
  t27 = adVirtLoc[1] - V0[1];
  t30 = R2[1] - R0[1];
  t33 = R1[1] - R0[1];
  adRealLoc[0] = (t1 * t2 * t11 + t13 * t14 * t11) * t18
      + (-t1 * t20 * t11 - t13 * t23 * t11) * t27 + R0[0];
  adRealLoc[1] = (t30 * t2 * t11 + t33 * t14 * t11) * t18
      + (-t30 * t20 * t11 - t33 * t23 * t11) * t27 + R0[1];
  //  printf("The virtual point (%.3e,%.3e) is really (%.3e,%.3e)\n",adVirtLoc[0],adVirtLoc[1],adRealLoc[0],adRealLoc[1]);

  return;
}

// void VirtualTransMesh::vPrintMesh()
// {
//   std::vector<virtualvert>::iterator iterVV;
//   std::vector<Cell *>::iterator iterC;

//   FILE* fpt;

//   fpt = fopen("virtmesh.tec","w");

//   fprintf(fpt,"TITLE = \"TECplot of VirtualMesh\"\n");
//   fprintf(fpt,"VARIABLES = \"X\", \"Y\"\n");
//   fprintf(fpt,"ZONE N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",vecvvVerts.size(),vecpCCells.size());
//   for (iterVV=vecvvVerts.begin(),iterVV!=vecvvVerts.end();iterVV++){
//     fprintf(fpt,"%lf  %lf\n",(*iterVV).dvX,(*iterVV).dvY);
//   }
//   for (iter
// 	 fclose(fpt);
// }

void
vTriCircumcenter(double adA[2], double adB[2], double adC[2], double adRes[2])
{
  double adRow0[] =
    { adA[0] - adB[0], adA[1] - adB[1] };
  double adRow1[] =
    { adA[0] - adC[0], adA[1] - adC[1] };

  double dRHS0 = 0.5 * (dDOT2D (adA, adA) - dDOT2D(adB, adB));
  double dRHS1 = 0.5 * (dDOT2D (adA, adA) - dDOT2D(adC, adC));

  double dDet = adRow0[0] * adRow1[1] - adRow1[0] * adRow0[1];
  double dDet0 = dRHS0 * adRow1[1] - dRHS1 * adRow0[1];
  double dDet1 = dRHS1 * adRow0[0] - dRHS0 * adRow1[0];

  adRes[0] = dDet0 / dDet;
  adRes[1] = dDet1 / dDet;

  double dComp = dDIST2D(adRes, adA);
  double dTest = dDIST2D(adRes, adB);
  assert(0 == iFuzzyComp(dComp, dTest));
  dTest = dDIST2D(adRes, adC);
  assert(0 == iFuzzyComp(dComp, dTest));
}

void
vGetVirtualLengths(double adAr[2], double adBr[2], double adCr[2],
		   double adNr[2], double dLAB, double dLBC, double dLCA,
		   double dRetLengths[3])
{
  //this function takes 3 vertex locations, a point within the triangle, and the associated lengths of the triangle. It returns a set of lengths A-N,B-N,C-N used for point insertion

  double t1, t2, t5, t8, t10, t11, t16, t18, t19, t20, t21, t23, t24, t30, t31,
      t35, t41, t42, t44, t45, t47, t51, t54, t57, t60, t61, lAN, t64, t65, lBN,
      t69, t70, t73, t74, lCN;
  //this is generated by maple:

  t1 = adCr[1];
  t2 = adAr[1];
  t5 = adBr[0];
  t8 = adAr[0];
  t10 = adCr[0];
  t11 = adBr[1];
  t16 = 1 / (-t5 * t1 + t5 * t2 + t8 * t1 + t10 * t11 - t10 * t2 - t8 * t11);
  t18 = dLAB * dLAB;
  t19 = dLBC * dLBC;
  t20 = dLCA * dLCA;
  t21 = t18 - t19 + t20;
  t23 = t21 / dLAB;
  t24 = t11 - t2;
  t30 = adNr[0] - t8;
  t31 = (dLAB * (-t1 + t2) * t16 + t23 * t24 * t16 / 2.0) * t30;
  t35 = t5 - t8;
  t41 = adNr[1] - t2;
  t42 = (-dLAB * (-t10 + t8) * t16 - t23 * t35 * t16 / 2.0) * t41;
  t44 = fabs(t31 + t42);
  t45 = t44 * t44;
  t47 = t21 * t21;
  t51 = sqrt(4.0 * t20 - t47 / t18);
  t54 = t51 * t24 * t16 * t30;
  t57 = t51 * t35 * t16 * t41;
  t60 = fabs(t54 / 2.0 - t57 / 2.0);
  t61 = t60 * t60;
  lAN = sqrt(t45 + t61);
  t64 = fabs(t31 + t42 - dLAB);
  t65 = t64 * t64;
  lBN = sqrt(t65 + t61);
  t69 = fabs(t31 + t42 - t23 / 2.0);
  t70 = t69 * t69;
  t73 = fabs(t54 / 2.0 - t57 / 2.0 - t51 / 2.0);
  t74 = t73 * t73;
  lCN = sqrt(t70 + t74);

  dRetLengths[0] = lAN;
  dRetLengths[1] = lBN;
  dRetLengths[2] = lCN;

  for (int i = 0; i < 3; i++) {
    assert(dRetLengths[i] > 0);
  }
  // printf("Returning lengths [%.3e,%.3e,%.3e]\n",dRetLengths[0],dRetLengths[1],dRetLengths[2]);
  return;
}
double
AnisoRefinement::dAvgArea()
{
  GR_index_t i;
  double dArea = 0;

  pM2D->purgeAllEntities(); //just in case...

  for (i = 0; i < pM2D->getNumCells(); i++) {
    dArea += dAnisoArea(pM2D->getCell(i));
  }
  logMessage(1, "Total Ani-Area is %.3e\n", dArea);

  dArea /= pM2D->getNumCells();
  logMessage(1, "Avg Ani-Area per cell is %.3e\n", dArea);

  assert(dArea > 0);
  return dArea;
}

double
AnisoRefinement::dAvgCircRadius()
{
  double dRetVal = 0;
  double dCircRad;
  GR_index_t i;

  for (i = 0; i < pM2D->getNumCells(); i++) {
    dCircRad = dAnisoCircumRadius(pM2D->getCell(i));
    dRetVal += dCircRad * dCircRad;

    //    if(dCircRad>1e6){
    //        Cell* pC=pM2D->pCCell(i);
    //        //    printf("Vertices (%.3e,%.3e),(%.3e,%.3e),(%.3e,%.3e)\n",pC->pVVert(0)->dX(),pC->pVVert(0)->dY(),pC->pVVert(1)->dX(),pC->pVVert(1)->dY(),pC->pVVert(2)->dX(),pC->pVVert(2)->dY());
    //        for(int j=0;j<3;j++){
    // 	 printf("Metric [%.3e,%.3e,%.3e]\n",pC->pVVert(j)->dGetMetric(0),pC->pVVert(j)->dGetMetric(1),pC->pVVert(j)->dGetMetric(2));
    //        }
    //        assert(0);
    //     }

  }
  dRetVal /= 1.0 * pM2D->getNumCells();
  assert(dRetVal > 0);
  assert(dRetVal < 1e12);
  return sqrt(dRetVal);
}

double
AnisoRefinement::dAnisoArea(const Cell* const pC)
{
  return dAnisoArea(pC->getVert(0), pC->getVert(1), pC->getVert(2));
}

double
AnisoRefinement::dAnisoCircumRadius(const Cell* const pC)
{
  return dAnisoCircumRadius(pC->getVert(0), pC->getVert(1), pC->getVert(2));
}

double
AnisoRefinement::dAnisoArea(const Vert* const pVVert1,
			    const Vert* const pVVert2,
			    const Vert* const pVVert3)
{
  double adMet[3];
  int i;

  for (i = 0; i < iNMetricValues; i++) {
    adMet[i] = pVVert1->getMetric(i) + pVVert2->getMetric(i)
	+ pVVert3->getMetric(i);
    adMet[i] /= 3.0;
  }
  double dLen12 = dMetricLength(pVVert1->getCoords(), pVVert2->getCoords(),
				adMet);
  double dLen23 = dMetricLength(pVVert2->getCoords(), pVVert3->getCoords(),
				adMet);
  double dLen31 = dMetricLength(pVVert3->getCoords(), pVVert1->getCoords(),
				adMet);
  //  printf("Lengths %.3e,%3e,%.3e\n",dLen12,dLen23,dLen31);
  return dArea(dLen12, dLen23, dLen31);
}

//double AnisoRefinement::dLogAnisoArea(const Vert* const pVVert1,
//		const Vert* const pVVert2,
//		const Vert* const pVVert3)
//{
//	double adMet1[3];
//	double adMet2[3];
//	double adMet3[3];
//	double adMet[3];
//
//	int i;
//
//	GRUMMP::TMOPQual TQ;
//
//	for (i = 0; i < iNMetricValues; i++) {
//		adMet1[i] = pVVert1->getMetric(i);
//		adMet2[i] = pVVert2->getMetric(i);
//		adMet3[i] = pVVert3->getMetric(i);
//	}
//
//	TQ.dLogMetricCell(adMet1,adMet2,adMet3,adMet);
//
//	double dLen12 = dMetricLength(pVVert1->getCoords(), pVVert2->getCoords(), adMet);
//	double dLen23 = dMetricLength(pVVert2->getCoords(), pVVert3->getCoords(), adMet);
//	double dLen31 = dMetricLength(pVVert3->getCoords(), pVVert1->getCoords(), adMet);
//	//  printf("Lengths %.3e,%3e,%.3e\n",dLen12,dLen23,dLen31);
//	return dArea(dLen12, dLen23, dLen31);
//}

double
AnisoRefinement::dAnisoCircumRadius(const Vert* const pVVert1,
				    const Vert* const pVVert2,
				    const Vert* const pVVert3)
{
  double adMet[iNMetricValues];
  int i;

  for (i = 0; i < iNMetricValues; i++) {
    adMet[i] = pVVert1->getMetric(i) + pVVert2->getMetric(i)
	+ pVVert3->getMetric(i);
    adMet[i] /= 3.0;
  }
  double dLen12 = dMetricLength(pVVert1->getCoords(), pVVert2->getCoords(),
				adMet);
  double dLen23 = dMetricLength(pVVert2->getCoords(), pVVert3->getCoords(),
				adMet);
  double dLen31 = dMetricLength(pVVert3->getCoords(), pVVert1->getCoords(),
				adMet);
  //  printf("Lengths %.3e,%3e,%.3e\n",dLen12,dLen23,dLen31);
  return dCircumRadius(dLen12, dLen23, dLen31);

}

//returns  the circumradius based only on length
double
dCircumRadius(double dL1, double dL2, double dL3)
{
  double dRetVal;

  dRetVal = (dL1 + dL2 + dL3) * (-dL1 + dL2 + dL3) * (dL1 - dL2 + dL3)
      * (dL1 + dL2 - dL3);

  if (dRetVal <= 0) {
    return LARGE_DBL;
    assert(0);
  }

  dRetVal = dL1 * dL2 * dL3 / sqrt(dRetVal);

  return dRetVal;
}

double
dArea(double dL1, double dL2, double dL3)
{
  double dRetVal;

  dRetVal = dL1 + dL2 + dL3;
  dRetVal *= (-dL1 + dL2 + dL3);
  dRetVal *= (dL1 - dL2 + dL3);
  dRetVal *= (dL1 + dL2 - dL3);
  assert(dRetVal > 0);
  dRetVal = 0.25 * sqrt(dRetVal);

  return dRetVal;
}

double
dAniBiggness(const Cell* const pC)
{
  assert(0);
  double dL1, dL2, dL3;
  double dRetVal;
  assert(pC->isValid());
  //  dL1=pC->pFFace(0)->dGetAniLength();
  //  dL2=pC->pFFace(1)->dGetAniLength();
  // dL3=pC->pFFace(2)->dGetAniLength();
  //  return dCircumRadius(dL1,dL2,dL3);

  // printf("Lengths [%.3e,%.3e,%.3e]\n",dL1,dL2,dL3);
  //returns the area of the cell
  dRetVal = dL1 + dL2 + dL3;
  dRetVal *= (-dL1 + dL2 + dL3);
  dRetVal *= (dL1 - dL2 + dL3);
  dRetVal *= (dL1 + dL2 - dL3);
  if (dRetVal <= 0) {
    //  printf("Bad Lengths [%.3e,%.3e,%.3e]\n",dL1,dL2,dL3);
    return (LARGE_DBL - 1);

    // assert(0);
  }
  //assert(dRetVal>0);
  // printf("Size is1/4*sqrt(%.3e)\n",dRetVal);
  dRetVal = 0.25 * sqrt(dRetVal);
  return dRetVal;

}

double
dAniPoorness(const Cell* const pC)
{
  assert(0);
  assert(pC->getNumFaces() == 3);
  //  double dLA = pC->pFFace(0)->dGetAniLength();
  // double dLB = pC->pFFace(1)->dGetAniLength();
  // double dLC = pC->pFFace(2)->dGetAniLength();
  double dLA, dLB, dLC;
  //  printf("\nLengths %.3e,%.3e,%.3e\n",dLA,dLB,dLC);
  double da, db, dc;

  da = (dLB * dLB + dLC * dLC - dLA * dLA);
  da /= 2.0 * dLB * dLC;
  if (fabs(da) > 1)
    return LARGE_DBL;
  da = acos(da);

  db = (dLA * dLA + dLC * dLC - dLB * dLB);
  db /= 2.0 * dLA * dLC;
  if (fabs(db) > 1)
    return LARGE_DBL;
  db = acos(db);

  dc = (dLB * dLB + dLA * dLA - dLC * dLC);
  dc /= 2.0 * dLB * dLA;
  if (fabs(dc) > 1)
    return LARGE_DBL;
  dc = acos(dc);

  //  printf("Angles %.3e,%.3e,%.3e\n",da,db,dc);
  //printf("return %.3e\n",max(da,max(db,dc)));
  return (1.0 / sin(min(da, min(db, dc))));
}

// void AnisoRefinement::vWriteLengths()
// {
//   int iC,iF;
//   Cell* pC;
//   Face* pF;

//   for(iC=0;iC<pM2D->iNumCells();iC++){
//     pC=pM2D->pCCell(iC);
//       printf("Cell %d (%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n has lengths ",iC,pC->pVVert(0)->dX(),pC->pVVert(0)->dY(),pC->pVVert(1)->dX(),pC->pVVert(1)->dY(),pC->pVVert(2)->dX(),pC->pVVert(2)->dY());
//     for(iF=0;iF<pC->iNumFaces();iF++){
//       pF=pC->pFFaceOpposite(pC->pVVert(iF)));
//      printf("%.3e ",dMetricLength(pF->pVVert(0)->adCoords);
//     }
//     printf("\n\n");
//   }
//   return;
// }

int
AnisoRefinement::iAniSwap(Face* pF)
{

  //The purpose of this function is to empty the circumcircles of the two cells adjacent to pF. If this is not possible, here come the heuristics
  int iRetVal = 0;

  //  if(pF->qIsBadFace()) printf(".");

  if (!pF->doFullCheck())
    return (0);
  //   return (0);
  // Never swap a face that divides two regions in a multi-region mesh.
  if (pF->getFaceLoc() == Face::eBdryTwoSide
      || pF->getType() == Face::eMultiEdge)
    return (0);
  assert(pM2D->getSwapType() != eNoSwap);

  //@@ Case: one or both cells is not a tri, including boundaries
  Cell *pCLeft = pF->getLeftCell();
  Cell *pCRight = pF->getRightCell();
  if (pF->isBdryFace())
    return 0;

  if (!pCLeft->isValid() || !pCRight->isValid())
    return (0);
  if (pCLeft->getType() != Cell::eTriCell
      || pCRight->getType() != Cell::eTriCell)
    return (0);

  Vert *pVVertA = pF->getVert(0);
  Vert *pVVertB = pF->getVert(1);
  TriCell *pTCLeft = dynamic_cast<TriCell*>(pCLeft);
  TriCell *pTCRight = dynamic_cast<TriCell*>(pCRight);
  assert(pTCLeft != NULL);
  assert(pTCRight != NULL);

  Vert *pVVertC = pTCLeft->getOppositeVert(pF);
  Vert *pVVertD = pTCRight->getOppositeVert(pF);
  // printf("Wants to swap (%.3e,%.3e)-(%.3e,%.3e)\n",pVVertA->dX(),pVVertA->dY(),pVVertB->dX(),pVVertB->dY());
  double dLenAB, dLenBC, dLenAC, dLenBD, dLenAD;
  double adMetABC[iNMetricValues];
  double adMetABD[iNMetricValues];
  double dRat;
  int i;
  //  bool qIsEmptyCC;
  bool qDoSwap = qDoesSwapImprove(pF);

  if (!qDoSwap) {
    for (i = 0; i < iNMetricValues; i++) {
      adMetABC[i] = pVVertA->getMetric(i) + pVVertB->getMetric(i)
	  + pVVertC->getMetric(i);
      adMetABC[i] /= 3.0;
    }
    for (i = 0; i < iNMetricValues; i++) {
      adMetABD[i] = pVVertA->getMetric(i) + pVVertB->getMetric(i)
	  + pVVertD->getMetric(i);
      adMetABD[i] /= 3.0;
    }
    dLenAB = dMetricLength(pVVertA->getCoords(), pVVertB->getCoords(),
			   adMetABC);
    dLenAC = dMetricLength(pVVertA->getCoords(), pVVertC->getCoords(),
			   adMetABC);
    dLenBC = dMetricLength(pVVertB->getCoords(), pVVertC->getCoords(),
			   adMetABC);
    dRat = 1.0 * dLenAB
	/ dMetricLength(pVVertA->getCoords(), pVVertB->getCoords(), adMetABD);
    dLenAD = dRat
	* dMetricLength(pVVertA->getCoords(), pVVertD->getCoords(), adMetABD);
    dLenBD = dRat
	* dMetricLength(pVVertB->getCoords(), pVVertD->getCoords(), adMetABD);
    pF->setBadFace(
	!qCircumCirclesEmpty(dLenAB, dLenAC, dLenBC, dLenAD, dLenBD));
  }
  if (qDoSwap) {
    //  if(pF->qIsBadFace()) return 0;

    //Repeat process except A->C B->D C->B D->A

    for (i = 0; i < iNMetricValues; i++) {
      adMetABC[i] = pVVertC->getMetric(i) + pVVertD->getMetric(i)
	  + pVVertB->getMetric(i);
      adMetABC[i] /= 3.0;
    }
    for (i = 0; i < iNMetricValues; i++) {
      adMetABD[i] = pVVertC->getMetric(i) + pVVertD->getMetric(i)
	  + pVVertA->getMetric(i);
      adMetABD[i] /= 3.0;
    }

    dLenAB = dMetricLength(pVVertC->getCoords(), pVVertD->getCoords(),
			   adMetABC);
    dLenAC = dMetricLength(pVVertC->getCoords(), pVVertB->getCoords(),
			   adMetABC);
    dLenBC = dMetricLength(pVVertD->getCoords(), pVVertB->getCoords(),
			   adMetABC);

    dRat = 1.0 * dLenAB
	/ dMetricLength(pVVertC->getCoords(), pVVertD->getCoords(), adMetABD);
    dLenAD = dRat
	* dMetricLength(pVVertC->getCoords(), pVVertA->getCoords(), adMetABD);
    dLenBD = dRat
	* dMetricLength(pVVertD->getCoords(), pVVertA->getCoords(), adMetABD);
    pF->setBadFace(
	!qCircumCirclesEmpty(dLenAB, dLenAC, dLenBC, dLenAD, dLenBD));

    int iOrientA = checkOrient2D(pVVertA, pVVertD, pVVertC);
    int iOrientB = checkOrient2D(pVVertB, pVVertC, pVVertD);
    if (iOrientA == -1 || iOrientB == -1 || iOrientA == 0 || iOrientB == 0) {
      //The circumcircle(s) are non-empty but the face can't be swapped
      //   printf("Wants to swap (%.3e,%.3e)-(%.3e,%.3e) but can't\n",pVVertA->dX(),pVVertA->dY(),pVVertB->dX(),pVVertB->dY());
      //    assert2(0,"Can't Swap... but want's to\n");

      pF->setBadFace(true);
      //  iFixFace(pF);
      return 0;
    }

    iRetVal += pM2D->reconfigure(pF);
    pM2D->sendEvents();
  }

  return iRetVal;
}

bool
AnisoRefinement::qDoesSwapImprove(const Face* const pF)
{
  //Choose the configuration with the smallest maximum circumradius
  assert(!pF->isBdryFace());
  const Vert* pVA = pF->getVert(0);
  const Vert* pVB = pF->getVert(1);
  const Vert* pVC = pF->getLeftCell()->getOppositeVert(pF);
  const Vert* pVD = pF->getRightCell()->getOppositeVert(pF);

  int i;

  double dOldCircRad1 = dAnisoCircumRadius(pVA, pVB, pVC);
  double dOldCircRad2 = dAnisoCircumRadius(pVA, pVB, pVD);

  double dNewCircRad1 = dAnisoCircumRadius(pVA, pVC, pVD);
  double dNewCircRad2 = dAnisoCircumRadius(pVB, pVC, pVD);

  double dLenAB = dMetricLength(pF);
  double dLenAC = dMetricLength(pF->getLeftCell()->getOppositeFace(pVB));
  double dLenBC = dMetricLength(pF->getLeftCell()->getOppositeFace(pVA));
  double dLenBD = dMetricLength(pF->getRightCell()->getOppositeFace(pVA));
  double dLenAD = dMetricLength(pF->getRightCell()->getOppositeFace(pVB));

  double dMet1[3];
  double dMet2[3];

  for (i = 0; i < 3; i++) {
    dMet1[i] = (pVB->getMetric(i) + pVC->getMetric(i) + pVD->getMetric(i))
	/ 3.0;
    dMet2[i] = (pVA->getMetric(i) + pVC->getMetric(i) + pVD->getMetric(i))
	/ 3.0;
  }

  double dLenCD = 0.5
      * (dMetricLength(pVC->getCoords(), pVD->getCoords(), dMet1)
	  + dMetricLength(pVC->getCoords(), pVD->getCoords(), dMet2));

  double dMinOldL1 = min(dLenAB, min(dLenAC, dLenBC)) / dOldCircRad1;
  double dMinOldL2 = min(dLenAB, min(dLenAD, dLenBD)) / dOldCircRad2;

  double dMinNewL1 = min(dLenAD, min(dLenAC, dLenCD)) / dNewCircRad1;
  double dMinNewL2 = min(dLenCD, min(dLenBD, dLenBC)) / dNewCircRad2;

  return (min(dMinOldL1, dMinOldL2) < min(dMinNewL1, dMinNewL2));

}

int
AnisoRefinement::iFixFace(Cell* pC)
{
  //iFixFace the longest edge of the cell
  //assert(0);
  Face* pFLongest;
  Face* pF;
  double dLen;
  double dLongest = 0;
  int i;
  double adMet[iNMetricValues];

  for (i = 0; i < iNMetricValues; i++) {
    adMet[i] = pC->getVert(0)->getMetric(i) + pC->getVert(1)->getMetric(i)
	+ pC->getVert(2)->getMetric(i);
    adMet[i] /= 3.0;
  }

  for (i = 0; i < pC->getNumFaces(); i++) {
    pF = pC->getFace(i);
    dLen = dMetricLength(pF->getVert(0)->getCoords(),
			 pF->getVert(1)->getCoords(), adMet);
    if (dLen > dLongest) {
      dLongest = dLen;
      pFLongest = pF;
    }
  }
  // printf("Longest face is (%.3e)\n",dLen);
  return iFixFace(pFLongest);
}

int
AnisoRefinement::iFixFace(Face* pF)
{

  assert(pF->isValid());
  if (pF->isBdryFace())
    assert(0);
  //  int i;
  static int iNumAdded = 0;

  Face* pFLongest = pF;
  //   dLongLength=dMetricLength(pF);

  //   //Split that face at the Anisotropic Midpoint, or on the boundary if it encroaches

  //  printf("Face (%.3e,%.3e)-(%.3e,%.3e) with shell (%.3e,%.3e) and (%.3e,%.3e)\n",pF->pVVert(0)->dX(),pF->pVVert(0)->dY(),pF->pVVert(1)->dX(),pF->pVVert(1)->dY(),pF->pCCell(0)->pVVertOpposite(pF)->dX(),pF->pCCell(0)->pVVertOpposite(pF)->dY(),pF->pCCell(1)->pVVertOpposite(pF)->dX(),pF->pCCell(1)->pVVertOpposite(pF)->dY());

  double adMidPoint[2];
  vAniMidPoint(pFLongest, adMidPoint);
  //  printf("Want to insert at (%.3e,%.3e)\n",adMidPoint[0],adMidPoint[1]);

  if (!(pFLongest->isBdryFace())) {
    int iBdryEncr = iCheckBoundary(adMidPoint, pFLongest->getCell(0));
    if (iBdryEncr > 0) {
      iNumAdded += iBdryEncr;
      //           printf("inserted %d for fixing +bdry\n",iNumAdded);
      return iNumAdded;
      //   //   printf("Interior\n");
    }
  }
  // std::set<Vert*> setpVVerts;
  // setpVVerts.insert(pFLongest->pVVert(0));
  // setpVVerts.insert(pFLongest->pVVert(1));
  // vInterpolateMetric(pVNew,&setpVVerts);
  iNumAdded++;
  //  printf("Point (%.3e,%.3e) being inserted to fix face\n",pVNew->dX(),pVNew->dY());

  pSI2D->insertPointOnFace(adMidPoint, pFLongest);
  //  printf("inserted %d for fixing\n",iNumAdded);

  return iNumAdded;
}

bool
qCircumCirclesEmpty(double dLenDiag, double dLenL1, double dLenL2,
		    double dLenR1, double dLenR2)
{
  double imA[2] =
    { 0, 0 };
  double imB[2] =
    { dLenDiag, 0 };
  double imC[2];
  double imD[2];

  imC[0] = 0.5 * (dLenDiag * dLenDiag - dLenL2 * dLenL2 + dLenL1 * dLenL1)
      / dLenDiag;
  imC[1] = (dLenL1 * dLenL1 - imC[0] * imC[0]);
  if (imC[1] < 0)
    return false;
  imC[1] = sqrt(imC[1]);

  imD[0] = 0.5 * (dLenDiag * dLenDiag - dLenR2 * dLenR2 + dLenR1 * dLenR1)
      / dLenDiag;
  imD[1] = (dLenR1 * dLenR1 - imD[0] * imD[0]);
  if (imD[1] < 0)
    return false;
  imD[1] = -sqrt(imD[1]);

  int iOrientA = checkOrient2D(imA, imD, imC);
  int iOrientB = checkOrient2D(imB, imC, imD);
  if (iOrientA == -1 || iOrientB == -1 || (iOrientA == 0 || iOrientB == 0)) //changed to or .
    return true;

  else
    return (isInCircle(imA, imB, imC, imD) != 1);
  //iIncircle returns one if the imD is inside the circumcircle of imA,imB,imC//If the vertex is on the circumcircle (returns 0), that's close enough to being outside.
}

int
AnisoRefinement::iAniCoarsen(CellQueue* pCellQueue, int iMaxNewBig = 10)
{

  if (iMaxNewBig == 0)
    iMaxNewBig = 10;
  //  bool qKeepBdry=false;  //no longer works
  bool qVerbose = false;
  //  int iStartVert=pM2D->iNumVerts();
  size_t iNumStartBig = pCellQueue->iBig();
  printf("Coarsening verts unless %d big cells are created\n", iMaxNewBig);

  int iRemoved = 0;
  bool qRemoved = false;
  //  double adMet[iNMetricValues];
  //  double adMidPoint[2];
  Face* pFtoDel;
  int test;
  while (pCellQueue->iSmall() != 0
      && pCellQueue->iBig() < (iNumStartBig + iMaxNewBig)) {
    //while (iRemoved<1100){
    qRemoved = false;
    //   pCellQueue->vTest();
    if (qVerbose)
      printf("CellQueue has %zd:big %zd:small\n", pCellQueue->iBig(),
	     pCellQueue->iSmall());

    pFtoDel = pCellQueue->pFSmallestFace();
    if (!pFtoDel)
      continue;
    assert(pFtoDel->getVert(0)->isValid());
    assert(pFtoDel->getVert(1)->isValid());
    if (qVerbose) {
      printf("Collapsing (%.6e,%.6e)-(%.6e,%.6e)\n", pFtoDel->getVert(0)->x(),
	     pFtoDel->getVert(0)->y(), pFtoDel->getVert(1)->x(),
	     pFtoDel->getVert(1)->y());
    }
    assert(!pFtoDel->getVert(0)->isDeleted());
    assert(!pFtoDel->getVert(1)->isDeleted());

    //   for(int i=0;i<iNMetricValues;i++){
    //       adMet[i]=pFtoDel->pVVert(0)->dGetMetric(i)+pFtoDel->pVVert(1)->dGetMetric(i);
    //       adMet[i]*=0.5;
    //     }

    //     adMidPoint[0]=0.5*(pFtoDel->pVVert(0)->dX()+pFtoDel->pVVert(1)->dX());
    //     adMidPoint[1]=0.5*(pFtoDel->pVVert(0)->dY()+pFtoDel->pVVert(1)->dY());
    assert(pFtoDel->doFullCheck());
    Vert* pVVerttoDel = NULL;
    Vert* pVVertOther = NULL;
    Vert* pVTemp = NULL;
    double dDet0 = 0.0;
    double dDet1 = 0.0;
    //Rules:
    //Delete the vertex with the smallest metric determinant OR
    //the non-boundary vertex if only one vertex is on the boundary

    dDet0 =
	(pFtoDel->getVert(0)->getMetric(0) * pFtoDel->getVert(0)->getMetric(2))
	    - (pFtoDel->getVert(0)->getMetric(1)
		* pFtoDel->getVert(0)->getMetric(1));
    dDet1 =
	(pFtoDel->getVert(1)->getMetric(0) * pFtoDel->getVert(1)->getMetric(2))
	    - (pFtoDel->getVert(1)->getMetric(1)
		* pFtoDel->getVert(1)->getMetric(1));

    //   printf(" fabs(log(dDet0/dDet1)) is %.3e\n",fabs(log(dDet0/dDet1)))

    if (fabs(log10(dDet0 / dDet1)) < 2) { //The metric doesn't change too much

      if (dDet0 > dDet1) {
	pVVerttoDel = pFtoDel->getVert(1);
	pVVertOther = pFtoDel->getVert(0);
      }
      else {
	pVVerttoDel = pFtoDel->getVert(0);
	pVVertOther = pFtoDel->getVert(1);
      }

      if (pVVerttoDel->isBdryVert()) {
	if (!pVVertOther->isBdryVert()) {
	  pVTemp = pVVerttoDel;
	  pVVerttoDel = pVVertOther;
	  pVVertOther = pVTemp;
	}
      }
      // This clause looks like it's a geometry-specific hack for a
      // particular case.  This is VERY BAD!  Don't do it!
      //       if((pVVerttoDel->y()==0 && fabs(pVVerttoDel->x())==500 )||(pVVertOther->y()==0 && fabs(pVVertOther->x())==500 ))
      //       {
      //           pFtoDel->vSetDoNotDelete(true);
      //           pCellQueue->vRemoveFace(pFtoDel);
      //           goto leap;
      //       }

      //   vAniMidPoint(pFtoDel,adMidPoint);
      if (qVerbose)
	printf("Calling qRemoveVertCheaply for (%.3e,%.3e)\n", pVVerttoDel->x(),
	       pVVerttoDel->y());

      qRemoved = (pM2D->removeVertCheaply(pVVerttoDel, test));

      if (!qRemoved) {
	qRemoved = (pM2D->removeVertCheaply(pVVertOther, test));
      }
    }
    if (qRemoved) {
      // if(qVerbose) printf("Successfully removed a Vertex\n");
      //	printf("To (%.6e,%.6e)\n",adMidPoint[0],adMidPoint[1]);
      //	pFtoDel->pVVert(0)->vSetCoords(2,adMidPoint);
      //	pFtoDel->pVVert(0)->vSetMetric(adMet);
      //	printf("(second) Vertex removed\n");
      iRemoved++;
      if (iRemoved % 100 == 0) {
	printf("Removed %d\n", iRemoved);

	//	pCellQueue->vTest();
      }
      //	if(iRemoved>2200) qVerbose=true;
    }
    else {
      if (qVerbose)
	printf("Couldn't delete Face (%.3e,%.3e)-(%.3e,%.3e)\n",
	       pFtoDel->getVert(0)->x(), pFtoDel->getVert(0)->y(),
	       pFtoDel->getVert(1)->x(), pFtoDel->getVert(1)->y());
      pFtoDel->vSetDoNotDelete(true);
      //	printf("Couldn't remove this one\n");
      pCellQueue->vRemoveFace(pFtoDel);
    }

//		writeVTK(*pM2D, "After_DuringCoarsen");
    if (qVerbose)
      printf("Removed %d\n", iRemoved);

  }

  pM2D->purgeAllEntities();
//	 pM2D->iSwap_deprecated(10); /////////////////////////////////////////////////////////////////////////////////////////////////
  // pCellQueue->vReset();

  printf("There are now %zd big, and %zd small\n", pCellQueue->iBig(),
	 pCellQueue->iSmall());
  printf("Removed a net total of %d vertices\n", iRemoved);
  //  qKeepBdry=true;
  //  pCellQueue->vTest();
  return iRemoved;
}

void
AnisoRefinement::vImproveMesh(CellQueue* pCellQueue, int iMaxVerts,
			      double dRate = 0.9, bool qVerbose = false)
{

  assert(dRate < 1 && dRate >= 0);
  pM2D->purgeAllEntities();

  int iNInnerLoops = 1;
  int iNVCurrent = pM2D->iNumVerts();
  int iNVAdded = 0;
  int iNVCoarsened = 0;

  char strOutName[FILE_NAME_LEN];
  if (qVerbose) {
    sprintf(strOutName, "initial.dat");
    vWriteDistribution(strOutName);
  }
  //First, add as many verts as you need to/can
  if (iNVCurrent < iMaxVerts)
    iNVAdded = iRefine(pCellQueue, iMaxVerts - iNVCurrent);
  //return;
  if (qVerbose) {
    sprintf(strOutName, "ref0.dat");
    vWriteDistribution(strOutName);
  }
  writeVTK(*pM2D, "After_DuringInsertion");
  //Then, coarsen until you create iRate*iNVAdded new big cells
  if (iNVAdded < 100)
    iNVAdded = 100;
  iNVCoarsened = iAniCoarsen(pCellQueue, int(dRate * iNVAdded));
  if (qVerbose) {
    sprintf(strOutName, "cut0.dat");
    vWriteDistribution(strOutName);
  }

  int iLoop = 0;
  while (iNVCoarsened != 0 && iLoop < iNInnerLoops) {
    iLoop++;

    iNVAdded = iRefine(pCellQueue, iNVCoarsened);
    if (qVerbose) {
      sprintf(strOutName, "ref%d.dat", iLoop);
      vWriteDistribution(strOutName);
    }
    iNVCoarsened = iAniCoarsen(pCellQueue, int(dRate * iNVAdded));
    if (qVerbose) {
      sprintf(strOutName, "cut%d.dat", iLoop);
      vWriteDistribution(strOutName);
    }
  }
  double dTest = 0.3;
  //int iFinalAdd =int (max(1.0*iNVCoarsened,dTest*(iMaxVerts-iNVCurrent)));
  int iFinalAdd = iNVCoarsened;
  iNVAdded = iRefine(pCellQueue, iFinalAdd * dRate * dTest);
  if (qVerbose) {
    sprintf(strOutName, "final.dat");
    vWriteDistribution(strOutName);
  }

  return;

}

void
AnisoRefinement::vWriteDistribution(char strOutName[])
{
  printf("Writing Distribution to:  %s\n", strOutName);
  FILE *fpt;
  fpt = fopen(strOutName, "w");
  pM2D->purgeAllEntities();

  Cell* pC = NULL;
  double dX, dY;

  for (GR_index_t i = 0; i < pM2D->getNumCells(); i++) {
    pC = pM2D->getCell(i);
    dY = sqrt(3) * dAnisoCircumRadius(pC);
    dX = LARGE_DBL;
    for (int iF = 0; iF < pC->getNumFaces(); iF++) {
      dX = min(dX, dMetricLength(pC->getFace(iF)));
    }
    fprintf(fpt, "%e  %e\n", dX, dY);
  }
  fclose(fpt);
}

bool
AnisoRefinement::qIsCurvedBdryOkay(Face* pF, const double *adVertLocC)
{
  bool qOkay;
  assert(pF->isBdryFace());
  int iNumGaussPoints = 18;
  double adS[iNumGaussPoints];
  // double adBeginRatio[iNumGaussPoints];
  // double adEndRatio[iNumGaussPoints];

  double dS = 0.5 * (1 + sqrt(1. / 3.));

  adS[0] = 0.5; // one gauss point
  // adBeginRatio[0] = 0;
  // adEndRatio[0] = 1.0;

  adS[1] = 1 - dS; // two gauss points
  // adBeginRatio[1] = 0;
  // adEndRatio[1] = 0.5;
  adS[2] = dS;
  // adBeginRatio[2] = 0.5;
  // adEndRatio[2] = 1.0;

  adS[3] = 0.5 - sqrt(0.15); // three gauss points
  // adBeginRatio[3] = 0;
  // adEndRatio[3] = 5./18.;
  adS[4] = 0.5;
  // adBeginRatio[4] = 5./18.;
  // adEndRatio[4] = 13./18.;
  adS[5] = 0.5 + sqrt(0.15);
  // adBeginRatio[5] = 13./18.;
  // adEndRatio[5] = 1.;

  adS[6] = 0.25; // two times one gauss points
  // adBeginRatio[6] = 0;
  // adEndRatio[6] = 0.5;
  adS[7] = 0.75;
  // adBeginRatio[7] = 0.5;
  // adEndRatio[7] = 1.0;

  adS[8] = 0.5 - 0.5 * dS; // two times two gauss points
  // adBeginRatio[8] = 0;
  // adEndRatio[8] = 0.25;
  adS[9] = 0.5 * dS;
  // adBeginRatio[9] = 0.25;
  // adEndRatio[9] = 0.5;
  adS[10] = 1 - 0.5 * dS;
  // adBeginRatio[10] = 0.5;
  // adEndRatio[10] = 0.75;
  adS[11] = 0.5 + 0.5 * dS;
  // adBeginRatio[11] = 0.75;
  // adEndRatio[11] = 1.;

  adS[12] = 0.5 * (0.5 - sqrt(0.15));
  // adBeginRatio[12] =  0;
  // adEndRatio[12] =  5./36.;
  adS[13] = 0.25;
  // adBeginRatio[13] =  5./36.;
  // adEndRatio[13] =  13./36.;
  adS[14] = 0.5 * (0.5 + sqrt(0.15));
  // adBeginRatio[14] =  13./36.;
  // adEndRatio[14] =  0.5;
  adS[15] = 0.5 + 0.5 * (0.5 - sqrt(0.15)); // two times three gauss points
  // adBeginRatio[15] = 0.5 + 0;
  // adEndRatio[15] = 0.5 + 5./36.;
  adS[16] = 0.5 + 0.25;
  // adBeginRatio[16] = 0.5 +5./36.;
  // adEndRatio[16] = 0.5 +13./36.;
  adS[17] = 0.5 + 0.5 * (0.5 + sqrt(0.15));
  // adBeginRatio[17] = 0.5 + 13./36.;
  // adEndRatio[17] = 1.;

  Cell *pC = pF->getLeftCell();
  if ((pC->getType() != Cell::eBdryEdge) && (pC->getType() != Cell::eTriBFace)
      && (pC->getType() != Cell::eQuadBFace)
      && (pC->getType() != Cell::eIntBdryEdge)
      && (pC->getType() != Cell::eIntTriBFace)
      && (pC->getType() != Cell::eIntQuadBFace)) {
    pC = pF->getRightCell();
  }
  assert(
      (pC->getType() == Cell::eBdryEdge) || (pC->getType() == Cell::eTriBFace)
	  || (pC->getType() == Cell::eQuadBFace)
	  || (pC->getType() == Cell::eIntBdryEdge)
	  || (pC->getType() == Cell::eIntTriBFace)
	  || (pC->getType() == Cell::eIntQuadBFace));
  double adGaussLoc[2];
  //     double adBeginLoc[2];
  //     double adEndLoc[2];
  //     double dLength;
  //double adNormal[2];
  const double *adVertLocA;
  const double *adVertLocB;

  double beg_param, end_param, gauss_param, face_length;
  adVertLocA = pC->getVert(0)->getCoords();
  adVertLocB = pC->getVert(1)->getCoords();
  //double adLocA=*adVertLocA;
  // double adLocB=*adVertLocB;
  // Identify bdry geometry info
  BdryEdgeBase* bdry_edge = NULL;
  bool reverse_normal;
  GRCurve* curve = NULL;
  CubitVector gauss_coord, normal;
  reverse_normal = false;
  if (pF->isBdryFace()) {
    bdry_edge = dynamic_cast<BdryEdgeBase*>(pF->getLeftCell());
    if (!bdry_edge) {
      bdry_edge = dynamic_cast<BdryEdgeBase*>(pF->getRightCell());
      assert(bdry_edge);
      if (!bdry_edge->isForward())
	reverse_normal = true;
    }
    else {
      if (!bdry_edge->isForward())
	reverse_normal = true;
    }
    curve = bdry_edge->getCurve();
    //  BdryPatch2D *pBP = dynamic_cast<BdryPatch2D*>((dynamic_cast<BFace*>(pC))->pPatchPointer());
    // Done

    //		if (bdry_edge->isForward()) {
    //			beg_param = bdry_edge->getVert0Param();
    //			end_param = bdry_edge->getVert1Param();
    //		}
    //		else {
    //			beg_param = bdry_edge->getVert1Param();
    //			end_param = bdry_edge->getVert0Param();
    //		}
    beg_param = min(bdry_edge->getVert0Param(), bdry_edge->getVert1Param());
    end_param = max(bdry_edge->getVert0Param(), bdry_edge->getVert1Param());

    assert(iFuzzyComp(beg_param, end_param) == -1);

    face_length = curve->get_curve_geom()->arc_length(beg_param, end_param);
  }
  qOkay = true;
  for (int iGP = 0; iGP < iNumGaussPoints; iGP++) {

    //      pBP->vGeodesicDistanceRatioBetweenTwoPoints(adBeginRatio[iGP], adVertLocA, adVertLocB, adBeginLoc);
    //      pBP->vGeodesicDistanceRatioBetweenTwoPoints(adEndRatio[iGP], adVertLocA, adVertLocB, adEndLoc);

    //      dLength = fabs(pBP->dGeodesicDistance(pBP->dParameterAtPoint(adBeginLoc),pBP->dParameterAtPoint(adEndLoc)));

    // Find correct parametric location of Gauss point
    // Find physical location of Gauss point
    //  pBP->vGeodesicDistanceRatioBetweenTwoPoints(adS[iGP], adVertLocA, adVertLocB, adGaussLoc);
    gauss_param = curve->get_curve_geom()->param_at_arc_length(
	beg_param, adS[iGP] * face_length);
    curve->get_curve_geom()->coord_at_param(gauss_param, gauss_coord);
    //gauss_coord.x()=-gauss_coord[0];
    //gauss_coord.y()=-gauss_coord[1];
    adGaussLoc[0] = gauss_coord[0];
    adGaussLoc[1] = gauss_coord[1];
    // Find normal
    //pBP->vUnitNormal(adGaussLoc, adNormal);
    curve->get_curve_geom()->unit_normal(gauss_param, normal);
    if (reverse_normal)
      normal *= -1.;
    int iOrACB, iOrACg, iOrCBg;

    // Use Gauss point info
    iOrACB = checkOrient2D(adVertLocA, adVertLocC, adVertLocB);
    iOrACg = checkOrient2D(adVertLocA, adVertLocC, adGaussLoc);
    iOrCBg = checkOrient2D(adVertLocC, adVertLocB, adGaussLoc);

    if (!(iOrACg == iOrACB && iOrCBg == iOrACB)) {
      qOkay = false;
      // 	printf("VertA = %e   %e\n",adVertLocA[0],adVertLocA[1]);
      // 	printf("VertB = %e   %e\n",adVertLocB[0],adVertLocB[1]);
      // 	printf("VertC = %e   %e\n",adVertLocC[0],adVertLocC[1]);
      // 	printf("Gauss = %e   %e\n",adGaussLoc[0],adGaussLoc[1]);

      // 	printf("iOrACB,iOrACg,iOrCBg = (%d,%d,%d)\n",iOrACB,iOrACg,iOrCBg);

      //	assert2(0,"Boundary Protection didn't work\n");
    }
  }
  return qOkay;
}

void
AnisoRefinement::vFixForCurve()
{
  pM2D->purgeAllEntities();
  GR_index_t iBF;
  Face* pF = NULL;
  int iTest;
  for (iBF = 0; iBF < pM2D->getNumBdryFaces(); iBF++) {
    pF = pM2D->getBFace(iBF)->getFace(0);
    if (!(qIsCurvedBdryOkay(pF,
			    pF->getCell(0)->getOppositeVert(pF)->getCoords()))) {
      assert(
	  pM2D->removeVertCheaply(pF->getCell(0)->getOppositeVert(pF), iTest));
      iBF = 0;
    }
  }
  return;
}
#ifdef HAVE_MESQUITE
bool
AnisoRefinement::bRemoveVertexwithQual(Face* const InteriorFaceVertRemove)
{

  std::set<Vert*>::iterator iter, iterE;
  std::set<Cell*>::iterator iterC, iterCE;

  std::set<Cell*> spCJunk;
  std::set<Vert*> spVJunk;
  std::set<Cell*> spC0;
  std::set<Cell*> spC1;
  std::set<Vert*> spV0;
  std::set<Vert*> spV1;
  Vert *pV0 = InteriorFaceVertRemove->getVert(0);
  Vert *pV1 = InteriorFaceVertRemove->getVert(1);
  Vert *pVDel;
  Vert *pVOther;

//	if there's no point with a request for deletion, do nothing
  if (pV0->isDeletionRequested() == false
      && pV1->isDeletionRequested() == false)
    return false;

//	Decide which vert to delete with TMOP qual, depending on an average quality of its neighbors
  else {
    double Qual0 = 0;
    double Qual1 = 0;
    double dCount = 0;
    findNeighborhoodInfo(pV0, spC0, spVJunk);
    iterC = spC0.begin();
    iterCE = spC0.end();

    for (; iterC != iterCE; iterC++) {
      Cell *spC = (*iterC);

      double dT[9];
      GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
      GRUMMP::TMOPQual TQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, spC->getVert(0),
			     spC->getVert(1), spC->getVert(2));
      Qual0 += TQ.dEval(QM, CellSkel::eTriCell, dT,
			spC->getVert(0)->getCoords(),
			spC->getVert(1)->getCoords(),
			spC->getVert(2)->getCoords());
      dCount++;
    }

    Qual0 = Qual0 / dCount;

    dCount = 0;

    findNeighborhoodInfo(pV1, spC1, spVJunk);
    iterC = spC1.begin();
    iterCE = spC1.end();

    for (; iterC != iterCE; iterC++) {
      Cell *spC = (*iterC);

      double dT[9];
      GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
      GRUMMP::TMOPQual TQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, spC->getVert(0),
			     spC->getVert(1), spC->getVert(2));
      Qual1 += TQ.dEval(QM, CellSkel::eTriCell, dT,
			spC->getVert(0)->getCoords(),
			spC->getVert(1)->getCoords(),
			spC->getVert(2)->getCoords());
      dCount++;
    }

    Qual1 = Qual1 / dCount;

    if (Qual0 >= Qual1) {
      pVDel = pV0;
      pVOther = pV1;
    }
    else {
      pVDel = pV1;
      pVOther = pV0;
    }
  }

// Delete the point with edge collapse, if fail try removing the point cheaply by swapping
  int iNewSwaps;
  bool bDeleted = false;

  std::set<Cell*> spC;
  std::set<Vert*> spV;
  findNeighborhoodInfo(pVDel, spC, spV);

  spV.erase(pVDel);
  if (pVDel->getVertType() != Vert::eBdryApex && pVDel->getVertType() != Vert::ePseudoSurface) { // && pVDel->getVertType() != Vert::eBdry){
    if (pVDel->isDeletionRequested()) {
      bDeleted = pM2D->removeVert(pVDel, iNewSwaps);
    }
  }

//	if removing the point is successful, set all it's neigbouring verties to keep
  if (bDeleted == true) {
    iter = spV.begin();
    iterE = spV.end();
    for (; iter != iterE; ++iter) {
      Vert* psV = (*iter);
      psV->markToKeep();
    }
  }

//	if removing the point failed, try removing the other point
  if (bDeleted == false) {
    if (pVOther->getVertType() != Vert::eBdryApex && pVDel->getVertType() != Vert::ePseudoSurface) { // && pVOther->getVertType() != Vert::eBdry){
      findNeighborhoodInfo(pVOther, spC, spV);
      spV.erase(pVOther);

      if (pVOther->isDeletionRequested() == true) {
	bDeleted = pM2D->removeVert(pVOther, iNewSwaps);
      }
    }
  }

//	if removing the other point is successful, set all it's neighbouring points to keep
  if (bDeleted == true) {
    iter = spV.begin();
    iterE = spV.end();
    for (; iter != iterE; ++iter) {
      Vert* psV = (*iter);
      psV->markToKeep();
    }
  }

//	if removing both points fails, mark the points to keep
  if (bDeleted == false) {
    pV0->markToKeep();
    pV1->markToKeep();
  }

  return bDeleted;
}

bool
AnisoRefinement::bRemoveVertex(Vert *pVert)
{

  int iNewSwaps;
  bool bDeleted = false;

//	Try deleting vert by edge collapse, if fails try deleting cheaply by swapping
  bDeleted = pM2D->removeVert(pVert, iNewSwaps);

// if failed to delete vert, mark point to keep
  if (bDeleted == false) {
    pVert->markToKeep();
  }

  return bDeleted;
}

AnisoRefinement3D::AnisoRefinement3D(VolMesh* pVMesh_) :
    pVMesh(pVMesh_), iOrder(2), iNMetricValues(5)
{
}

bool
AnisoRefinement3D::bRemoveVertexwithQual(Vert *pV0, Vert *pV1)
{

  std::set<Vert*>::iterator iter, iterE;
  std::set<Cell*>::iterator iterC, iterCE;

  std::set<Cell*> spCJunk;
  std::set<Vert*> spVJunk;
  std::set<Cell*> spC0;
  std::set<Cell*> spC1;
  std::set<Vert*> spV0;
  std::set<Vert*> spV1;
  Vert *pVDel;
  Vert *pVOther;

//	if none of the points can be deleted...
  if (pV0->isDeletionRequested() == false
      && pV1->isDeletionRequested() == false)
    return false;

//	computes average neighbouring quality for both points, chooses to delete the one with the worst quality
  else {
    double Qual0 = 0;
    double Qual1 = 0;
    double dCount = 0;
    findNeighborhoodInfo(pV0, spC0, spVJunk);
    iterC = spC0.begin();
    iterCE = spC0.end();

    for (; iterC != iterCE; iterC++) {
      Cell *spC = (*iterC);

      double dT[9];
      GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
      GRUMMP::TMOPQual TQ;

      if (checkOrient3D(spC->getVert(0), spC->getVert(1), spC->getVert(2),
			spC->getVert(3)) == -1) {
	TQ.SetTargetfromMetric(dT, CellSkel::eTet, spC->getVert(0),
			       spC->getVert(2), spC->getVert(1),
			       spC->getVert(3));
	Qual0 += TQ.dEval(QM, CellSkel::eTet, dT, spC->getVert(0)->getCoords(),
			  spC->getVert(2)->getCoords(),
			  spC->getVert(1)->getCoords(),
			  spC->getVert(3)->getCoords());
      }
      else {
	TQ.SetTargetfromMetric(dT, CellSkel::eTet, spC->getVert(0),
			       spC->getVert(1), spC->getVert(2),
			       spC->getVert(3));
	Qual0 += TQ.dEval(QM, CellSkel::eTet, dT, spC->getVert(0)->getCoords(),
			  spC->getVert(1)->getCoords(),
			  spC->getVert(2)->getCoords(),
			  spC->getVert(3)->getCoords());
      }
      dCount++;
    }

    Qual0 = Qual0 / dCount;
    dCount = 0;

    findNeighborhoodInfo(pV1, spC1, spVJunk);
    iterC = spC1.begin();
    iterCE = spC1.end();

    for (; iterC != iterCE; iterC++) {
      Cell *spC = (*iterC);

      double dT[9];
      GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
      GRUMMP::TMOPQual TQ;

      if (checkOrient3D(spC->getVert(0), spC->getVert(2), spC->getVert(1),
			spC->getVert(3)) == -1) {
	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, spC->getVert(0),
			       spC->getVert(2), spC->getVert(1),
			       spC->getVert(3));
	Qual1 += TQ.dEval(QM, CellSkel::eTriCell, dT,
			  spC->getVert(0)->getCoords(),
			  spC->getVert(2)->getCoords(),
			  spC->getVert(1)->getCoords(),
			  spC->getVert(3)->getCoords());
      }
      else {
	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, spC->getVert(0),
			       spC->getVert(1), spC->getVert(2),
			       spC->getVert(3));
	Qual1 += TQ.dEval(QM, CellSkel::eTriCell, dT,
			  spC->getVert(0)->getCoords(),
			  spC->getVert(1)->getCoords(),
			  spC->getVert(2)->getCoords(),
			  spC->getVert(3)->getCoords());
      }
      dCount++;
    }

    Qual1 = Qual1 / dCount;

    if (Qual0 >= Qual1) {
      pVDel = pV0;
      pVOther = pV1;
    }
    else {
      pVDel = pV1;
      pVOther = pV0;
    }
  }

  int iNewSwaps;
  bool bDeleted = false;

  std::set<Cell*> spC;
  std::set<Vert*> spV;
  findNeighborhoodInfo(pVDel, spC, spV);
  spV.erase(pVDel);

//	Delete the point byu edge collapse
  if ((pVDel->getVertType() != Vert::eBdryApex
      || pVDel->getVertType() != Vert::eBdry) && (pVDel->isDeletionRequested()))
    bDeleted = pVMesh->removeVert(pVDel, iNewSwaps);

//	if succes in deleting the point, set all other neighbouring point to keep... this to avoid over coarsening an area
  if (bDeleted == true) {
    iter = spV.begin();
    iterE = spV.end();
    for (; iter != iterE; ++iter) {
      Vert* psV = (*iter);
      psV->markToKeep();
    }
  }

//	if failed to delete the point, try to delete the other one
  if (bDeleted == false) {
    findNeighborhoodInfo(pVOther, spC, spV);
    spV.erase(pVOther);
    if ((pVOther->getVertType() == Vert::eBdryApex
	|| pVDel->getVertType() != Vert::eBdry)
	&& (pVOther->isDeletionRequested() == true))
      return false;
    bDeleted = pVMesh->removeVert(pVOther, iNewSwaps);
  }
// if succes in deleting the other point, keep all the other neigbouring point to avoid over coarsening
  if (bDeleted == true) {
    iter = spV.begin();
    iterE = spV.end();
    for (; iter != iterE; ++iter) {
      Vert* psV = (*iter);
      psV->markToKeep();
    }
  }

  //if false up to here... mark not to remove...
  if (bDeleted == false) {
    pV0->markToKeep();
    pV1->markToKeep();
  }

  return bDeleted;
}

bool
AnisoRefinement3D::bRemoveVertex(Vert *pVert)
{

  int iNewSwaps;
  bool bDeleted = false;

  bDeleted = pVMesh->removeVert(pVert, iNewSwaps);

  if (bDeleted == false)
    pVert->markToKeep();

  return bDeleted;
}

// At the moment only sets two different analytic metric functions: quarter circle along z axis and cross section along z axis
void
AnisoRefinement3D::vAssignMetricsToVertsTest(enum eMetricType Type)
{
  for (GR_index_t i = 0; i < pVMesh->getNumVerts(); i++) {

    Vert *pV = pVMesh->getVert(i);
    assert(pV->isValid());

    double adM[6];
    double d_hmin;
    double d_hmax;
    double d_alphaX;
    double d_alphaY;
    double d_alpha;
    double d_h1;
    double d_h2;
    double d_Theta;
    double d_X = pV->x();
    double d_Y = pV->y();

    switch (Type)
      {

      case e3DTestCross:

	d_hmin = 0.005;
	d_hmax = 0.1;
	d_alphaX = 20 * (fabs(d_X - 0.5));
	d_alphaY = 20 * (fabs(d_Y - 0.5));
	d_h1 = std::min((pow(2, d_alphaX)) * d_hmin, d_hmax);
	d_h2 = std::min((pow(2, d_alphaY)) * d_hmin, d_hmax);

	adM[0] = pow(d_h1, -2);
	adM[1] = 0;
	adM[2] = pow(d_h2, -2);
	adM[3] = 0;
	adM[4] = 0;
	adM[5] = 100;

	break;

      case e3DTestQuarterCircle:

	d_Theta = atan2(d_Y, d_X);
	d_hmax = 0.1;
	d_alpha = 10 * (fabs(0.75 - sqrt(pow(d_X, 2) + pow(d_Y, 2))));
	d_h1 = std::min(0.002 * pow(5, d_alpha), d_hmax);
	d_h2 = std::min(0.05 * pow(2, d_alpha), d_hmax);

	adM[0] = pow(d_h1, -2) * pow(cos(d_Theta), 2)
	    + pow(d_h2, -2) * pow(sin(d_Theta), 2);
	adM[1] = (pow(d_h1, -2) - pow(d_h2, -2)) * cos(d_Theta) * sin(d_Theta);
	adM[2] = pow(d_h1, -2) * pow(sin(d_Theta), 2)
	    + pow(d_h2, -2) * pow(cos(d_Theta), 2);
	adM[3] = 0;
	adM[4] = 0;
	adM[5] = 100;

	break;

      default:
	assert(0);
	break;
      }

    pV->setMetric(adM);
  }
}
#endif
