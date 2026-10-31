#include <assert.h>
#include <math.h>

#include "GR_Aniso.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Geometry.h"
#include "GR_SwapDecider.h"
#include "GR_Util.h"

using std::min;

using namespace GRUMMP;

#ifdef SIM_ANNEAL_TEST
extern double simAnnealTemp;
static bool simAnnealComp(const double a, const double b)
{
  assert(0 <= a && a <= 1);
  assert(0 <= b && b <= 1);
  double diff = 1/b - 1/a + 1.e-10;
  diff /= (simAnnealTemp);
  diff = (diff + 1) / 2;
  double rand = drand48();
  return (rand < diff);
}
#endif

double
CircRadfromLength(double a, double b, double c)
{
  assert(a >= 0 && b >= 0 && c >= 0);

  //if(!qValidTri(a,b,c)) return LARGE_DBL;

  double dCircRad = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c);
  if (dCircRad <= 0) {
    // vMessage(1,"\nZERO, ahhhhh!\n\n");
    return LARGE_DBL;
  }
  // assert(dCircRad > 0.00000);
  // vMessage(1,"CircRad is first %f\n",dCircRad);
  dCircRad = 1.0 * a * b * c / sqrt(dCircRad);
  //vMessage(1,"CircRad is %f\n\n",dCircRad);

  //   double  dCircRad =(a*a+b*b+c*c)/2;
  //   vMessage(1,"CircRad is first %f\n",dCircRad);
  //   dCircRad = dCircRad - a/b/b/c/c - b/a/a/c/c - c/a/a/b/b;
  //   vMessage(1,"CircRad is %f\n\n",dCircRad);
  //   if(dCircRad <= 0){
  //     vMessage(1,"\nZERO, ahhhhh!\n\n");
  //   }
  //   dCircRad = sqrt(1.0/dCircRad);

  return dCircRad;
}

double
AniLength(const double adVert1[2], const double adMetric1[3],
	  const double adVert2[2], const double adMetric2[3])
{

  double dLength, a, b;
  double t3, t5, t10, t12, t18, t22;

  //  vMessage(1,"The distance between (%.3f,%.3f) and (%.3f,%.3f) is to be computed\n",adVert1[0],adVert1[1],adVert2[0],adVert2[1]);

  //   vMessage(1,"The metric for Vert1 is [%.4f,%.4f,%.4f,%.4f]\n",adMetric1[0],adMetric1[2],adMetric1[1],adMetric1[3]);
  //   vMessage(1,"The metric for Vert1 is [%.4f,%.4f,%.4f,%.4f]\n",adMetric2[0],adMetric2[2],adMetric2[1],adMetric2[3]);
  //See maple sheet
  t3 = adVert2[0] - adVert1[0];
  t5 = adMetric1[0];
  t10 = adVert2[1] - adVert1[1];
  t12 = adMetric1[1];
  t18 = adMetric1[1];
  t22 = adMetric1[2];
  a = (t3 * (adMetric2[0] - t5) + t10 * (adMetric2[1] - t12)) * t3
      + (t3 * (adMetric2[1] - t18) + t10 * (adMetric2[2] - t22)) * t10;
  b = (t3 * t5 + t10 * t12) * t3 + (t3 * t18 + t10 * t22) * t10;

  //  vMessage(1,"t3 = %.3f, t4 = %.3f, t6 = %.3f, t10 = %.3f, t15 = %.3f, t17 = %.3f,a is %.3f,b is %.3f\n",t3,t4,t6,t10,t15,t17,a,b);

  if (fabs(a) < 1e-10 && b > 1e-10)
    dLength = sqrt(b);
  else if (a > 1e-10 && b < 1e-10)
    dLength = 2.0 / 3 * sqrt(a);
  else if (a < 1e-10 && b < 1e-10)
    dLength = 0;

  else {
    dLength = 2.0 / 3.0 * (sqrt(pow(a + b, 3.0)) - sqrt(b * b * b)) / a;
  }

  //  vMessage(1,"The coeff are [%.3f,%.3f,%.3f]\n",adMetric1[0],adMetric1[2],adMetric1[3]);
  //  vMessage(1,"And [%.3f,%.3f,%.3f]\n",adMetric2[0],adMetric2[2],adMetric2[3]);
  //  vMessage(1,"The poly constants are a=%.3f,b=%.3f\n",a,b);

  if (dLength <= 0) {
    logMessage(1, "\nThe poly constants are a=%.3f,b=%.3f\n", a, b);
    logMessage(1, "The distance between (%.3f,%.3f) and (%.3f,%.3f) is %.4f\n",
	       adVert1[0], adVert1[1], adVert2[0], adVert2[1], dLength);
    logMessage(1, "With metrics [%.3f,%.3f,%.3f] and [%.3f,%.3f,%.3f]\n\n",
	       adMetric1[0], adMetric1[1], adMetric1[2], adMetric2[0],
	       adMetric2[1], adMetric2[2]);
  }
  //      vMessage(1,"The distance between (%.3f,%.3f) and (%.3f,%.3f) is %.4f\n",adVert1[0],adVert1[1],adVert2[0],adVert2[1],dLength);

  assert(dLength >= 0);

  return dLength;
}

double
AniLength(const Vert * const pV1, const Vert * const pV2)
{

  double P1[2], P2[2], T1[4], T2[4];

  // vMessage(1,"The distance between (%.3f,%.3f) and (%.3f,%.3f) is to be computed\n",pV1->x(),pV1->dY(),pV2->x(),pV2->dY());

  if (pV2->getSpaceDimen() != 2) {
    logMessage(
	1,
	"The distance between (%.3f,%.3f)%d and (%.3f,%.3f)%d is to be computed\n",
	pV1->x(), pV1->y(), pV1->getSpaceDimen(), pV2->x(), pV2->y(),
	pV2->getSpaceDimen());
  }
  assert(pV1->getSpaceDimen() == 2);
  assert(pV2->getSpaceDimen() == 2);
  P1[0] = pV1->x();
  P1[1] = pV1->y();
  P2[0] = pV2->x();
  P2[1] = pV2->y();
  // for(i=0;i<4;i++){
  T1[0] = pV1->getMetric(0);
  T1[1] = pV1->getMetric(1);
  T1[2] = pV1->getMetric(2);
  T2[0] = pV2->getMetric(0);
  T2[1] = pV2->getMetric(1);
  T2[2] = pV2->getMetric(2);

  //  }

  return AniLength(P1, T1, P2, T2);
}

/// Unless overridden in a derived class, this function decides whether
/// the current configuration is preferable to the alternate
/// configuration by seeking to maximize the minimum value of the
/// quality function for the current or speculative tets.
bool
SwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FC) const
{
  if (!m_pQM) {
    // How could we have gotten here?  The subclasses that don't set a
    // specific shape quality measure are supposed to override this
    // function.  This is definitely a library bug, if it occurs.
    assert(0);
    return false;
  }
  // If we've gotten this far, we know for sure that it's topologically
  // valid to do the reconfiguration.  But is it advisable in terms of
  // mesh quality?

  // First, make sure that the verts that are definitely supposed to be
  // here actually are.
  Vert *pVA = FC.getVertA();
  Vert *pVB = FC.getVertB();
  Vert *pVC = FC.getVertC();
  Vert *pVD = FC.getVertD();
  VALIDATE_INPUT(
      pVA->isValid() && pVB->isValid() && pVC->isValid() && pVD->isValid());
  // Verts are ordered like this: A and B are the ends of the existing
  // face; C and D are the opposite verts in the two tris incident on
  // the face.

  // If these values never get changed, there'll never be any face swaps.
  double qualBefore = 1, qualAfter = 0;

  // Old configuration
  double qualADB = m_pQM->eval(pVA, pVD, pVB);
  double qualABC = m_pQM->eval(pVA, pVB, pVC);
  qualBefore = min(qualADB, qualABC);

  // New configuration
  double qualADC = m_pQM->eval(pVA, pVD, pVC);
  double qualCDB = m_pQM->eval(pVC, pVD, pVB);
  qualAfter = min(qualADC, qualCDB);

  // In cases of (near) tie, don't do anything; this will help termination.
#ifdef SIM_ANNEAL_TEST
  return simAnnealComp(qualAfter, qualBefore);
#else
  return (fuzzyGreater(qualAfter, qualBefore));
#endif
}

bool
DelaunaySwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FC) const
{
  if (!FC.isSwappable())
    return false;
  int iIncircleRes = isInCircle(FC.getVertA(), FC.getVertB(), FC.getVertC(),
				FC.getVertD());
  return (iIncircleRes == 1);
}

bool
UniformDegreeSwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FC) const
{
  if (!FC.isSwappable())
    return false;
  Vert *verts[4];
  verts[0] = FC.getVertA();
  verts[1] = FC.getVertB();
  verts[2] = FC.getVertC();
  verts[3] = FC.getVertD();

  //   if (pVA->qIsBdryVert() || pVB->qIsBdryVert() ||
  //       pVC->qIsBdryVert() || pVD->qIsBdryVert()) {
  //     return false;
  //   }

  int degDiff[4];
  for (int ii = 0; ii < 4; ii++) {
    if (verts[ii]->getVertType() == Vert::eBdryApex)
      return false;
    int deg = verts[ii]->getNumFaces();
    bool qBdry = verts[ii]->isBdryVert();
    if (qBdry) {
      degDiff[ii] = 2 * (deg - 1) - 6;
    }
    else {
      degDiff[ii] = deg - 6;
    }
  }

  int before = (degDiff[0] * degDiff[0] + degDiff[1] * degDiff[1]
      + degDiff[2] * degDiff[2] + degDiff[3] * degDiff[3]);

  degDiff[0] -= (verts[0]->isBdryVert() ? 2 : 1);
  degDiff[1] -= (verts[1]->isBdryVert() ? 2 : 1);
  degDiff[2] += (verts[2]->isBdryVert() ? 2 : 1);
  degDiff[3] += (verts[3]->isBdryVert() ? 2 : 1);

  int after = (degDiff[0] * degDiff[0] + degDiff[1] * degDiff[1]
      + degDiff[2] * degDiff[2] + degDiff[3] * degDiff[3]);

  return (after < before);

}
bool
AnisoSwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FC) const
{
  if (!FC.isSwappable())
    return false;
  Face *pF;
  pF = FC.getFace();
  if (!pF->doFullCheck())
    return (0);
  assert(pF->doFullCheck());
  // assert (eSwapType () != eNone);
  Cell *pCLeft = pF->getLeftCell();
  Cell *pCRight = pF->getRightCell();
  Vert *pVVertA = pF->getVert(0);
  Vert *pVVertB = pF->getVert(1);
  TriCell *pTCLeft = dynamic_cast<TriCell*>(pCLeft);
  TriCell *pTCRight = dynamic_cast<TriCell*>(pCRight);
  assert(pTCLeft != NULL);
  assert(pTCRight != NULL);

  Vert *pVVertL = pTCLeft->getOppositeVert(pF);
  Vert *pVVertR = pTCRight->getOppositeVert(pF);

  int iOrientA = checkOrient2D(pVVertA, pVVertR, pVVertL);
  int iOrientB = checkOrient2D(pVVertB, pVVertL, pVVertR);
  if (iOrientA != 1 || iOrientB != 1)
    return 0;

  assert(pVVertA->isValid());
  assert(pVVertB->isValid());
  assert(pVVertL->isValid());
  assert(pVVertR->isValid());

  double Area1, Area2;
  /*  if (pVVertA->getMetric(0)==1 || pVVertB->getMetric(0)==1 || pVVertL->getMetric(0)==1 || pVVertR->getMetric(0)==1)
   {
   double vectorAB[] = {pVVertA->x() - pVVertB->x(),
   pVVertA->y() - pVVertB->y() };
   double AB= dMAG2D(vectorAB);
   double vectorLR[] = {pVVertL->x() - pVVertR->x(),
   pVVertL->y() - pVVertR->y() };
   double LR= dMAG2D(vectorLR);
   double  vectorBR[] = {pVVertB->x() - pVVertR->x(),
   pVVertB->y() - pVVertR->y() };
   double BR= dMAG2D(vectorBR);
   double vectorBL[] = {pVVertB->x() - pVVertL->x(),
   pVVertB->y() - pVVertL->y() };
   double BL= dMAG2D(vectorBL);
   double vectorAR[] = {pVVertA->x() - pVVertR->x(),
   pVVertA->y() - pVVertR->y() };
   double AR= dMAG2D(vectorAR);
   double vectorAL[] = {pVVertA->x() - pVVertL->x(),
   pVVertA->y() - pVVertL->y() };
   double AL= dMAG2D(vectorAL);
   Area1 = CircRadfromLength(AL,AR,LR)+CircRadfromLength(BL,BR,LR);
   Area2 = CircRadfromLength(AL,AB,BL)+CircRadfromLength(AB,AR,BR);
   }
   else
   {*/
  double AL, AR, BL, BR, LR, AB;
  AL = AniLength(pVVertA, pVVertL);
  AR = AniLength(pVVertA, pVVertR);
  BL = AniLength(pVVertB, pVVertL);
  BR = AniLength(pVVertB, pVVertR);
  LR = AniLength(pVVertL, pVVertR);
  AB = AniLength(pVVertA, pVVertB);

  Area1 = CircRadfromLength(AL, AR, LR) + CircRadfromLength(BL, BR, LR);
  Area2 = CircRadfromLength(AL, AB, BL) + CircRadfromLength(AB, AR, BR);
  // }
  return (Area1 < Area2);
  /* double FrobeniusNormOld=0 , FrobeniusNormNew=0;

   double vectorAB[] = {pVVertA->x() - pVVertB->x(),
   pVVertA->y() - pVVertB->y() };
   double vectorLR[] = {pVVertR->x() - pVVertL->x(),
   pVVertR->y() - pVVertL->y() };
   double vectorBR[] = {pVVertR->x() - pVVertB->x(),
   pVVertR->y() - pVVertB->y() };
   double vectorBL[] = {pVVertB->x() - pVVertL->x(),
   pVVertB->y() - pVVertL->y() };
   double vectorAR[] = {pVVertA->x() - pVVertR->x(),
   pVVertA->y() - pVVertR->y() };
   double vectorAL[] = {pVVertL->x() - pVVertA->x(),
   pVVertL->y() - pVVertA->y() };
   double vectorAB[]={pTCLeft->pVVert(1)->x()-pTCLeft->pVVert(0)->x(), pTCLeft->pVVert(1)->y()-pTCLeft->pVVert(0)->y()};
   double vectorAL[]={pTCLeft->pVVert(2)->x()-pTCLeft->pVVert(1)->x(), pTCLeft->pVVert(2)->y()-pTCLeft->pVVert(1)->y()};
   double vectorBL[]={pTCLeft->pVVert(0)->x()-pTCLeft->pVVert(2)->x(), pTCLeft->pVVert(0)->y()-pTCLeft->pVVert(2)->y()};
   double vectorBR[]={pTCRight->pVVert(2)->x()-pTCRight->pVVert(0)->x(), pTCRight->pVVert(2)->y()-pTCRight->pVVert(0)->y()};
   double vectorAR[]={pTCRight->pVVert(0)->x()-pTCRight->pVVert(2)->x(), pTCRight->pVVert(0)->y()-pTCRight->pVVert(2)->y()};
   double vectorLR[]={pTCRight->pVVert(2)->x()-pTCLeft->pVVert(1)->x(), pTCRight->pVVert(2)->y()-pTCLeft->pVVert(1)->dY()};
   double w[2][2];
   w[0][0]=(vectorAL[1]-0.5*vectorAB[1])/(vectorAB[0]*vectorAL[1]-vectorAB[1]*vectorAL[0]);
   w[0][1]=(vectorAL[0]-0.5*vectorAB[0])/((vectorAB[1]*vectorAL[0]-vectorAL[1]*vectorAB[0]));
   w[1][0]=(sqrt(3)*vectorAB[1])/(2*(vectorAB[1]*vectorAL[0]-vectorAL[1]*vectorAB[0]));
   w[1][1]=(sqrt(3)*vectorAB[0])/(2*(vectorAL[1]*vectorAB[0]-vectorAB[1]*vectorAL[0]));

   double T[2][2];
   double Target[2][2];

   Target[0][0]=0;
   Target[0][1]=-sqrt((pVVertL->getMetric(2))-(pVVertL->getMetric(1)));
   Target[1][0]=sqrt(pVVertL->getMetric(0));
   Target[1][1]=pVVertL->getMetric(1);

   T[0][0]=-Target[0][1];//1;
   T[0][1]=2*(Target[0][0]-0.5*Target[0][1])/sqrt(3);//sqrt(3)/3;
   T[1][0]=-Target[1][1];//0;
   T[1][1]=2*(Target[1][0]-0.5*Target[1][1])/sqrt(3);//20*sqrt(3)/3;

   double WInvbyT[2][2];
   WInvbyT[0][0]=T[0][0]*w[0][0]+T[1][0]*w[0][1];
   WInvbyT[0][1]=T[0][1]*w[0][0]+T[1][1]*w[0][1];
   WInvbyT[1][0]=T[0][0]*w[1][0]+T[1][0]*w[1][1];
   WInvbyT[1][1]=T[0][1]*w[1][0]+T[1][1]*w[1][1];

   FrobeniusNormOld+=sqrt(WInvbyT[0][0]*WInvbyT[0][0]+WInvbyT[0][1]*WInvbyT[0][1]+WInvbyT[1][0]*WInvbyT[1][0]+WInvbyT[1][1]*WInvbyT[1][1]);

   w[0][0]=(vectorBR[1]+0.5*vectorAB[1])/(-vectorAB[0]*vectorBR[1]+vectorAB[1]*vectorBR[0]);
   w[0][1]=(vectorBR[0]+0.5*vectorAB[0])/((-vectorAB[1]*vectorBR[0]+vectorBR[1]*vectorAB[0]));
   w[1][0]=(-sqrt(3)*vectorAB[1])/(2*(-vectorAB[1]*vectorBR[0]+vectorBR[1]*vectorAB[0]));
   w[1][1]=(-sqrt(3)*vectorAB[0])/(2*(-vectorBR[1]*vectorAB[0]+vectorAB[1]*vectorBR[0]));

   Target[0][0]=0;
   Target[0][1]=-sqrt((pVVertR->getMetric(2))-(pVVertR->getMetric(1)));
   Target[1][0]=sqrt(pVVertR->getMetric(0));
   Target[1][1]=pVVertR->getMetric(1);

   T[0][0]=-Target[0][1];//1;
   T[0][1]=2*(Target[0][0]-0.5*Target[0][1])/sqrt(3);//sqrt(3)/3;
   T[1][0]=-Target[1][1];//0;
   T[1][1]=2*(Target[1][0]-0.5*Target[1][1])/sqrt(3);//20*sqrt(3)/3;

   WInvbyT[0][0]=T[0][0]*w[0][0]+T[1][0]*w[0][1];
   WInvbyT[0][1]=T[0][1]*w[0][0]+T[1][1]*w[0][1];
   WInvbyT[1][0]=T[0][0]*w[1][0]+T[1][0]*w[1][1];
   WInvbyT[1][1]=T[0][1]*w[1][0]+T[1][1]*w[1][1];

   FrobeniusNormOld+=sqrt(WInvbyT[0][0]*WInvbyT[0][0]+WInvbyT[0][1]*WInvbyT[0][1]+WInvbyT[1][0]*WInvbyT[1][0]+WInvbyT[1][1]*WInvbyT[1][1]);
   /////////////////////
   w[0][0]=(vectorAR[1]-0.5*vectorLR[1])/(vectorLR[0]*vectorAR[1]-vectorLR[1]*vectorAR[0]);
   w[0][1]=(vectorAR[0]-0.5*vectorLR[0])/((vectorLR[1]*vectorAR[0]-vectorAR[1]*vectorLR[0]));
   w[1][0]=(sqrt(3)*vectorLR[1])/(2*(vectorLR[1]*vectorAR[0]-vectorAR[1]*vectorLR[0]));
   w[1][1]=(sqrt(3)*vectorLR[0])/(2*(vectorAR[1]*vectorLR[0]-vectorLR[1]*vectorAR[0]));

   Target[0][0]=0;
   Target[0][1]=-sqrt((pVVertA->getMetric(2))-(pVVertA->getMetric(1)));
   Target[1][0]=sqrt(pVVertA->getMetric(0));
   Target[1][1]=pVVertA->getMetric(1);

   T[0][0]=-Target[0][1];//1;
   T[0][1]=2*(Target[0][0]-0.5*Target[0][1])/sqrt(3);//sqrt(3)/3;
   T[1][0]=-Target[1][1];//0;
   T[1][1]=2*(Target[1][0]-0.5*Target[1][1])/sqrt(3);//20*sqrt(3)/3;

   FrobeniusNormNew+=sqrt(WInvbyT[0][0]*WInvbyT[0][0]+WInvbyT[0][1]*WInvbyT[0][1]+WInvbyT[1][0]*WInvbyT[1][0]+WInvbyT[1][1]*WInvbyT[1][1]);

   w[0][0]=(vectorBL[1]+0.5*vectorLR[1])/(-vectorLR[0]*vectorBL[1]+vectorLR[1]*vectorBL[0]);
   w[0][1]=(vectorBL[0]+0.5*vectorLR[0])/((-vectorLR[1]*vectorBL[0]+vectorBL[1]*vectorLR[0]));
   w[1][0]=(-sqrt(3)*vectorLR[1])/(2*(-vectorLR[1]*vectorBL[0]+vectorBL[1]*vectorLR[0]));
   w[1][1]=(-sqrt(3)*vectorLR[0])/(2*(-vectorBL[1]*vectorLR[0]+vectorLR[1]*vectorBL[0]));

   Target[0][0]=0;
   Target[0][1]=-sqrt((pVVertB->getMetric(2))-(pVVertB->getMetric(1)));
   Target[1][0]=sqrt(pVVertB->getMetric(0));
   Target[1][1]=pVVertB->getMetric(1);

   T[0][0]=-Target[0][1];//1;
   T[0][1]=2*(Target[0][0]-0.5*Target[0][1])/sqrt(3);//sqrt(3)/3;
   T[1][0]=-Target[1][1];//0;
   T[1][1]=2*(Target[1][0]-0.5*Target[1][1])/sqrt(3);//20*sqrt(3)/3;

   FrobeniusNormNew+=sqrt(WInvbyT[0][0]*WInvbyT[0][0]+WInvbyT[0][1]*WInvbyT[0][1]+WInvbyT[1][0]*WInvbyT[1][0]+WInvbyT[1][1]*WInvbyT[1][1]);

   return (FrobeniusNormOld < FrobeniusNormNew);*/
}
//////////////////****************************/////////////////////////////
// This tie breaker is used only to pick a unique configuration in the
// event that a tie occurs trying to decide  how to swap.  The only
// thing which this choice of tiebreaker has to recommend it is that it
// gives a unique answer for any two pairs of verts (AB and CD).

static bool
qTieBreak(const Vert * pVVertA, const Vert * pVVertB, const Vert * pVVertC,
	  const Vert * pVVertD)
{
  unsigned long ulA = reinterpret_cast<unsigned long>(pVVertA);
  unsigned long ulB = reinterpret_cast<unsigned long>(pVVertB);
  unsigned long ulC = reinterpret_cast<unsigned long>(pVVertC);
  unsigned long ulD = reinterpret_cast<unsigned long>(pVVertD);
  if (min(ulA, ulB) > min(ulC, ulD))
    return true;
  else
    return false;
}

//@@ Determines whether swapping will improve the maximum face angle.
// Returns true to swap, false to remain the same.
bool
MinMaxAngleSwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FSI2D) const
{
  if (!FSI2D.isSwappable())
    return false;
  // This algorithm finds face angles between adjacent faces by dotting
  // their unit normals.  The largest magnitude loses.
  //
  // To prevent pairs from flopping back and forth, the tie-breaker is
  // invoked if the inequality is small.

  Vert *pVVertA = FSI2D.getVertA();
  Vert *pVVertB = FSI2D.getVertB();
  Vert *pVVertC = FSI2D.getVertC();
  Vert *pVVertD = FSI2D.getVertD();

  double dEps = 1.e-12;
  // All of these normals are intended to point inwards.
  double adNormDA[2], adNormBD[2], adNormCB[2], adNormAC[2];
  calcNormal2D(pVVertA->getCoords(), pVVertD->getCoords(), adNormDA);
  calcNormal2D(pVVertD->getCoords(), pVVertB->getCoords(), adNormBD);
  calcNormal2D(pVVertB->getCoords(), pVVertC->getCoords(), adNormCB);
  calcNormal2D(pVVertC->getCoords(), pVVertA->getCoords(), adNormAC);

  NORMALIZE2D(adNormDA);
  NORMALIZE2D(adNormBD);
  NORMALIZE2D(adNormCB);
  NORMALIZE2D(adNormAC);

  double dDotBDA = dDOT2D(adNormBD, adNormDA);
  double dDotACB = dDOT2D(adNormAC, adNormCB);
  double dDotDAC = dDOT2D(adNormDA, adNormAC);
  double dDotCBD = dDOT2D(adNormCB, adNormBD);

  double dMaxThis = max(dDotBDA, dDotACB);
  double dMaxOther = max(dDotDAC, dDotCBD);

  if (dMaxThis > dMaxOther + dEps)
    return true;
  else if (dMaxThis + dEps < dMaxOther)
    return false;
  else
    return qTieBreak(pVVertA, pVVertB, pVVertC, pVVertD);
}
bool
DelaunaySwapDeciderSurf::doFaceSwap(const FaceSwapInfo2D& FC) const
{
  if (!FC.isSwappable())
    return false;

  double adCircCent[3];
  calcCircumcenter3D(FC.getVertA(), FC.getVertB(), FC.getVertC(), adCircCent);

  double dRadius = dDIST3D(adCircCent, FC.getVertA()->getCoords());
  double dDistance = dDIST3D(adCircCent, FC.getVertD()->getCoords());

  switch (iFuzzyComp(dDistance, dRadius))
    {
    case -1:
      // Distance smaller than radius -> inside
      return true;
    case 1:
      // Distance larger than radius -> outside
      return false;
    case 0:
      return qTieBreak(FC.getVertA(), FC.getVertB(), FC.getVertC(),
		       FC.getVertD());
    default:
      assert(0);
      return eClean;
    }
}

bool
MetricSwapDecider2D::isMetricOK(Vert *pV) const
{

  double dDet;
  double adMetric[3];

  for (int i = 0; i < 3; i++) {
    adMetric[i] = pV->getMetric(i);
  }

  dDet = adMetric[0] * adMetric[2] - adMetric[1] * adMetric[1];

  if (dDet <= 0 || adMetric[0] <= 0) {
    printf("Metric [%.3e,%.3e,%.3e] is not positive definite\n", adMetric[0],
	   adMetric[1], adMetric[2]);
    return false;
  }
  return true;
}

bool
MetricSwapDecider2D::doFaceSwap(const FaceSwapInfo2D& FSI2D) const
{

  if (!m_pQM || m_pQM->getName() != "angles (2D)") {
    logMessage(0, "Need to set a MinMax Angle 2D as Quality Measure...\n");
    assert(0);
    return false;
  }

  if (!FSI2D.isSwappable())
    return false;

  Vert* pVA = FSI2D.getVertA();
  Vert* pVB = FSI2D.getVertB();
  Vert* pVC = FSI2D.getVertC();
  Vert* pVD = FSI2D.getVertD();

  assert(isMetricOK(pVA));
  assert(isMetricOK(pVB));
  assert(isMetricOK(pVC));
  assert(isMetricOK(pVD));

  assert(pVA->isValid());
  assert(pVB->isValid());
  assert(pVC->isValid());
  assert(pVD->isValid());

  int iOrientABC = checkOrient2D(pVA, pVB, pVC);
  int iOrientADB = checkOrient2D(pVA, pVD, pVB);
  int iOrientADC = checkOrient2D(pVA, pVD, pVC);
  int iOrientBCD = checkOrient2D(pVB, pVC, pVD);

  if (iOrientABC != 1 && iOrientADB != 1 && iOrientADC != 1 && iOrientBCD != 1)
    return 0;

  m_pQM->vEvalInMetric(bMetric);

  double qualBefore = 1, qualAfter = 0;

  // Old configuration
  double qualADB = m_pQM->eval(pVA, pVD, pVB);
  double qualABC = m_pQM->eval(pVA, pVB, pVC);
  qualBefore = min(qualADB, qualABC);

  // New configuration
  double qualADC = m_pQM->eval(pVA, pVD, pVC);
  double qualCDB = m_pQM->eval(pVC, pVD, pVB);
  qualAfter = min(qualADC, qualCDB);

  return (qualAfter > qualBefore);
}

//bool MinMaxAngleSwapDeciderSurf::doFaceSwap(const FaceSwapInfo2D& FSI2D) const
//{
//  if (!FSI2D.isSwappable()) return false;
//  // This algorithm finds face angles between adjacent faces by dotting
//  // their unit normals.  The largest magnitude loses.
//  //
//  // To prevent pairs from flopping back and forth, the tie-breaker is
//  // invoked if the inequality is small.
//
//  Vert *pVVertA = FSI2D.getVertA();
//  Vert *pVVertB = FSI2D.getVertB();
//  Vert *pVVertC = FSI2D.getVertC();
//  Vert *pVVertD = FSI2D.getVertD();
//  pVVertA->printVertInfo();
//  pVVertB->printVertInfo();
//  pVVertC->printVertInfo();
//  pVVertD->printVertInfo();
//  double dEps = 1.e-12;
//  // All of these normals are intended to point inwards.
//  // project to 2D. Orientation irrelevant, I think
//  double dA[3],dB[3],dC[3],dD[3];
//  double adNormDA[2], adNormBD[2], adNormCB[2], adNormAC[2];
//
//  projectTo2D(pVVertA,pVVertD,pVVertB,dA,dD,dB);
//  calcNormal2D(dA, dD, adNormDA);
//  calcNormal2D(dD, dB, adNormBD);
//
//  projectTo2D(pVVertB,pVVertC,pVVertA,dB,dC,dA);
//  calcNormal2D(dB, dC, adNormCB);
//  calcNormal2D(dC, dA, adNormAC);
//
//  NORMALIZE2D(adNormDA);
//  NORMALIZE2D(adNormBD);
//  NORMALIZE2D(adNormCB);
//  NORMALIZE2D(adNormAC);
//
//  double dDotBDA = (dDOT2D (adNormBD, adNormDA));
//  double dDotACB = (dDOT2D (adNormAC, adNormCB));
//  double dDotDAC = (dDOT2D (adNormDA, adNormAC));
//  double dDotCBD = (dDOT2D (adNormCB, adNormBD));
//
//  double dMaxThis = max(dDotBDA, dDotACB);
//  double dMaxOther = max(dDotDAC, dDotCBD);
//
//  if (dMaxThis > dMaxOther + dEps)
//    return true;
//  else if (dMaxThis + dEps < dMaxOther)
//    return false;
//  else
//    return qTieBreak(pVVertA, pVVertB, pVVertC, pVVertD);
//}

