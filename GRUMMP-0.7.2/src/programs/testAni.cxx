#define TESTMESH
//#define TESTVERT
#define MAXVERTS 20000
#define REF_RATIO 1.5 //1.3
//the approximate size of the new mesh relative to the old one
//#define ALPHA 0
#define ANGLE_TOLERANCE 28*M_PI/180.0 //The smallest angle considered okay

//   0   Maximum angle
//   1   Minimum angle
//   2   All angles
//   6   Aspect ratio I (inscribed to circum radius)
//   9   Aspect ratio II (based on area versus perimeter)
//   10  Aspect ratio III (shortest edge to circumradius)

//#define BdryFILENAME "/home/pagnutti/GRUMMP/doug/src/programs/Meshes/test"
//#define FILENAME "/home/pagnutti/GRUMMP/doug/src/programs/Meshes/test"
#define OutFILENAME

#include <iostream>
#include "GR_Vertex.h"
#include "GR_Aniso.h"
#include "math.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "GR_config.h"
#include "GR_AdaptPred.h"
#include "GR_GRGeom2D.h"
#include "GR_InsertionQueue.h"
#include "GR_Mesh2D.h"
#include "GR_LengthAniso2D.h"
#include "GR_LengthAniso3D.h"
#include "GR_misc.h"
#include "GR_ADT.h"
#include "GR_Cell.h"
#include "GR_CellCV.h"
#include "GR_CellSkel.h"
#include <vector>
#include <set>
#include "GR_CellQueue.h"
#include "GR_Aniso.h"
#include "GR_Geometry.h"

using namespace std;

// static void vIdentity(const double /*adLoc*/[], double adRetVal[])

// {
//   adRetVal[0]=1.0;
//   adRetVal[1]=0;
//   adRetVal[2]=1.0;
//   return;
// }

static void
vTest(const double adLoc[], double adRetVal[])
{

//  double angle = atan2(adLoc[1]-0.5,adLoc[0]-0.5);
//  double length = sqrt((adLoc[0]-0.5)*(adLoc[0]-0.5) + (adLoc[1]-0.5)*(adLoc[1]-0.5));
//  /*
//  double garb = 0.009;
//
//  //adRetVal[0]=garb+(1-garb)*sin(angle)*sin(angle);
//  //adRetVal[1]=(1-garb)*cos(angle)*sin(angle);
//  //adRetVal[2]=garb+(1-garb)*cos(angle)*cos(angle);
//
//  adRetVal[0]=garb+(1-garb)*sqrt(2)*length*sin(angle)*sin(angle);
//  adRetVal[1]=(1-garb)*sqrt(2)*length*cos(angle)*sin(angle);
//  adRetVal[2]=garb+(1-garb)*sqrt(2)*length*cos(angle)*cos(angle);
//
//  // Multiplying all elements in tensor by a scalar doesn't do anything, currently.
//  double factor = 1;
//  adRetVal[0] *= factor;
//  adRetVal[1] *= factor;
//  adRetVal[2] *= factor;
//  */
//
//	double t = angle+M_PI_2;
//	double D = 0.1;
//	double d = 0.1*(1-0.84*sin(M_PI*length*sqrt(2)));
//	double s = sin(t);
//	double c = cos(t);
//	double q = (c/D)*(c/D)+(s/d)*(s/d);
//	adRetVal[0] = (c*c)/(D*D)+(s*s)/(d*d);
//	adRetVal[1] = (c*s)/(D*D)-(c*s)/(d*d);
//	adRetVal[2] = (s*s)/(D*D)+(c*c)/(d*d);
//
////  if(adLoc[1]<0.05){
////    adRetVal[2]=1e4;
////  }
////
////  if(adLoc[0]>0.475 && adLoc[0]<0.525){
////    adRetVal[0]=1e4;
////  }
//
//  return;

  adRetVal[0] = 1;
  adRetVal[1] = 0;
  adRetVal[2] = 1;

  if (adLoc[1] < 0.05) {
    adRetVal[2] = 1e4;
  }

  if (adLoc[0] > 0.475 && adLoc[0] < 0.525) {
    adRetVal[0] = 1e4;
  }

  return;
}

int
main(int argc, char *argv[])
{
  assert(argc == 4);

  void
  (*pvfMetric)(const double adLoc[], double adRetVal[])=vTest;

  char strBaseFileName[FILE_NAME_LEN];
  char strMeshInFileName[FILE_NAME_LEN];
  char strMeshOutFileName[FILE_NAME_LEN];
  char strQualityFileName[FILE_NAME_LEN];
  char strMetricFileName[FILE_NAME_LEN];
  char strMeshInFileNameWithSuffix[FILE_NAME_LEN];

  snprintf(strBaseFileName, FILE_NAME_LEN, "%s", argv[1]);
  snprintf(strMeshInFileName, FILE_NAME_LEN, "%s", argv[2]);
  snprintf(strMeshOutFileName, FILE_NAME_LEN, "%s", argv[3]);

  char strPreSmooth[FILE_NAME_LEN];
  char strSwap[FILE_NAME_LEN];
  char strCoarsen[FILE_NAME_LEN];

  //  char strSwap1[FILE_NAME_LEN];
  //  char strSwap2[FILE_NAME_LEN];
  char strHeal[FILE_NAME_LEN];
  char strRefine[FILE_NAME_LEN];
  char strInput[FILE_NAME_LEN];

  snprintf(strPreSmooth, FILE_NAME_LEN, "%s_PreSmooth", strMeshOutFileName);
  snprintf(strSwap, FILE_NAME_LEN, "%s_swap", strMeshOutFileName);
  snprintf(strCoarsen, FILE_NAME_LEN, "%s_coarsen", strMeshOutFileName);

  // sprintf(strSwap1,"%s_swap1",strMeshOutFileName);
  // sprintf(strSwap2,"%s_swap2",strMeshOutFileName);
  snprintf(strHeal, FILE_NAME_LEN, "%s_heal", strMeshOutFileName);
  snprintf(strRefine, FILE_NAME_LEN, "%s_refine", strMeshOutFileName);
  snprintf(strInput, FILE_NAME_LEN, "%s_input", strMeshOutFileName);

  printf("Base file name: %s\n", strBaseFileName);
  printf("Mesh Input file name: %s.mesh\n", strMeshInFileName);
  printf("Mesh Output file name: %s.mesh\n", strMeshOutFileName);

  sprintf(strQualityFileName, "%s.qual", strMeshOutFileName);
  sprintf(strMetricFileName, "%s.metric", strMeshInFileName);
  sprintf(strMeshInFileNameWithSuffix, "%s.mesh", strMeshInFileName);

  GRUMMPInit("testAni");

  Mesh2D *pM2D;
  GRGeom2D *pG2D = NULL;

  pG2D = new GRGeom2D(strBaseFileName);
  assert(pG2D!=NULL);
  // Reading an existing mesh file
  pM2D = new Mesh2D(2);
  pM2D->readFromFile(strMeshInFileNameWithSuffix);
  pM2D->relateMeshToGeometry(pG2D);
  assert(pM2D->isValid());
  logMessage(1, "Number of cells read: %d\n", pM2D->getNumCells());
  // vWriteFile_Mesh2D(*pM2D,strInput); 

  LengthAniso2D* pLength = new LengthAniso2D(LengthAniso2D::eFunc, 1, 1,
  LARGE_DBL);
  pLength->provideFunc(pvfMetric);
  AnisoRefinement ARMesh2D(pM2D, pLength, 2, 0);

  pM2D->attachAnisoRef(&ARMesh2D);
  ARMesh2D.vHealBdry();
  //  vWriteFile_Mesh2D(*pM2D,strHeal);
  printf("Finishing Swapping...\n");
  pM2D->iSwap_deprecated(3);
  //  vWriteFile_Mesh2D(*pM2D,strSwap);

  int iNumNonEmpty = 0;
  GR_index_t i;
  for (i = 0; i < pM2D->getNumCells(); i++) {
    if (!pM2D->getCell(i)->isCircumcircleEmpty())
      iNumNonEmpty++;
  }
  printf(" %d/%d cells have non-empty CC\n", iNumNonEmpty, pM2D->getNumCells());

  double dRho = REF_RATIO;
  double dThetaMin = ANGLE_TOLERANCE;
  double dAdes = ARMesh2D.dAvgArea() / dRho;
  double dIdeal = 2 * sqrt(dAdes / sqrt(3));
  double dTau = dIdeal
      * (7.0 - 4.0 * sqrt(3.0) * sin(dThetaMin)
	  - 4.0 * cos(dThetaMin) * cos(dThetaMin))
      / (4.0 * cos(dThetaMin) * cos(dThetaMin) - 1.0);

  double dMaxCircRad = 1.0 / sqrt(3) * (dIdeal + (dTau));
  double dMinLength = dIdeal - dTau;

  printf(
      "RefRat=%.3e,dThetaMin=%.3e dAdes=%.3e, dIdeal=%.3e, dTau=%.3e, dMaxCR=%.3e, dMinL=%.3e\n",
      dRho, dThetaMin, dAdes, dIdeal, dTau, dMaxCircRad, dMinLength);

  assert(dMaxCircRad > 0);
  assert(dMinLength > 0);

  CellQueue* pCellQueue = new CellQueue(dynamic_cast<Mesh*>(pM2D), dMaxCircRad,
					dMinLength, true);
  int iMaxVerts = MAXVERTS;

  printf("Improving Mesh\n");
  ARMesh2D.vImproveMesh(pCellQueue, iMaxVerts, 0.1, false); //0.2 works well
  //  ARMesh2D.vHealBdry();
  // vWriteFile_Mesh2D(*pM2D,strPreSmooth); 

  pM2D->iSwap_deprecated(10);

  ARMesh2D.vFixForCurve();
//  writeFile_Mesh2D(*pM2D,strMeshOutFileName);
  pM2D->writeBdryPatchFile(strMeshOutFileName, "");
  writeVTK(*pM2D, strMeshOutFileName);

  delete pM2D;
  delete pG2D;

  return 1;
}
