#include <math.h>

#include "GR_EdgeCanon.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Geometry.h"
#include "GR_SwapDecider.h"
#include "GR_Util.h"

using std::min;

using namespace GRUMMP;

#undef SIM_ANNEAL_TEST
#ifdef SIM_ANNEAL_TEST
double simAnnealTemp = 1;
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

/// Unless overridden in a derived class, this function decides whether
/// the current configuration is preferable to the alternate
/// configuration by seeking to maximize the minimum value of the
/// quality function for the current or speculative tets.
bool
SwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
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
  Vert *pVE = FC.getVertE();
  VALIDATE_INPUT(
      pVA->isValid() && pVB->isValid() && pVC->isValid() && pVD->isValid()
	  && pVE->isValid());
  // Can't check whether they belong to a mesh; we don't have one, and
  // they might be verts that aren't in a mesh yet anyway.

  // Verts are ordered like this:
  //
  // 2->3, 3->2 swaps: verts A, B, C are always the ones on the internal
  // triangle in the two tet state.
  //
  // 2->2, 4->4 swaps: verts A, B, C are always the ones on the (a, for
  // 4->4) current internal triangle.  A and B are two of the four
  // coplanar verts (D and E are the others).
  //
  // Vert D is to the left of tri ABC (orientation(ABCD) is inverted).
  // Vert E is to the right (orientation(ABCE) is correct).
  //
  // Vert F (only for 4->4 swaps) is the analog to vert C, but on the
  // other side of plane ABDE.

  // If these values never get changed, there'll never be any face swaps.
  double qualBefore = 1, qualAfter = 0;
  switch (FC.getFaceCat())
    {
    case FaceSwapInfo3D::eT23:
      {
	// Old configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	qualBefore = min(qualACBD, qualABCE);

	// New configuration
	double qualDEAB = m_pQM->eval(pVD, pVE, pVA, pVB);
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	qualAfter = min(qualDEAB, min(qualDEBC, qualDECA));
	break;
      }
    case FaceSwapInfo3D::eT32:
      {
	// Old configuration
	double qualDEAB = m_pQM->eval(pVD, pVE, pVA, pVB);
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	qualBefore = min(qualDEAB, min(qualDEBC, qualDECA));

	// New configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	qualAfter = min(qualACBD, qualABCE);
	break;
      }
    case FaceSwapInfo3D::eT22:
      {
	// Can't swap T22 if the bdry is precious.
	if (!m_allowBdryChanges)
	  return false;
	// Old configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	qualBefore = min(qualACBD, qualABCE);

	// New configuration
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	qualAfter = min(qualDEBC, qualDECA);
	break;
      }
    case FaceSwapInfo3D::eT44:
      {
	Vert *pVF = FC.getVertF();
	VALIDATE_INPUT(pVF->isValid());
	// Old configuration
	double qualABDC = m_pQM->eval(pVA, pVB, pVD, pVC);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	double qualABEF = m_pQM->eval(pVA, pVB, pVE, pVF);
	double qualABFD = m_pQM->eval(pVA, pVB, pVF, pVD);
	qualBefore = min(min(qualABDC, qualABCE), min(qualABEF, qualABFD));

	// New configuration
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	double qualDEAF = m_pQM->eval(pVD, pVE, pVA, pVF);
	double qualDEFB = m_pQM->eval(pVD, pVE, pVF, pVB);
	qualAfter = min(min(qualDEBC, qualDECA), min(qualDEAF, qualDEFB));
	break;
      }
    default:
      return false;
    }
  // In cases of (near) tie, don't do anything; this will help termination. 
#ifdef SIM_ANNEAL_TEST
  return simAnnealComp(qualAfter, qualBefore);
#else
  return (fuzzyGreater(qualAfter, qualBefore));
#endif
}

bool
DelaunaySwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
{
  int iInSphereRes = isInsphere(FC.getVertA(), FC.getVertB(), FC.getVertC(),
				FC.getVertD(), FC.getVertE());

  switch (FC.getFaceCat())
    {
    case FaceSwapInfo3D::eT32:
      // Current tets are ABDE, BCDE, CADE.  Swap if pVVertE is on (0) or
      // outside (-1) the circumsphere of ABCD.  This is an arbitrary
      // tie-break; configurations which are co-spherical go to the two
      // tet configuration.
      return (iInSphereRes != 1);
    case FaceSwapInfo3D::eT23:
    case FaceSwapInfo3D::eT22:
    case FaceSwapInfo3D::eT44: // This one is the same as eT22 for all situations.
      // Current tets are ABCD, CBAE.  Swap if pVVertE is inside (1) the
      // circumsphere of ABCD.  Again, configurations which are
      // co-spherical go to the two-tet configuration.  For T22 and T44,
      // ties remain unchanged.
      return (iInSphereRes == 1);
    default:
      return (false);
    }
}

bool
UniformDegreeSwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
{
  // This first cut at the code requires that all verts be interior.
  if (FC.getVertA()->isBdryVert() || FC.getVertB()->isBdryVert()
      || FC.getVertC()->isBdryVert() || FC.getVertD()->isBdryVert()
      || FC.getVertE()->isBdryVert())
    return false;

  int iEdgeDegA = (FC.getVertA()->getNumFaces() + 8) / 3;
  int iEdgeDegB = (FC.getVertB()->getNumFaces() + 8) / 3;
  int iEdgeDegC = (FC.getVertC()->getNumFaces() + 8) / 3;
  int iEdgeDegD = (FC.getVertD()->getNumFaces() + 8) / 3;
  int iEdgeDegE = (FC.getVertE()->getNumFaces() + 8) / 3;

  int iEdgeDiffA = iEdgeDegA - 13;
  int iEdgeDiffB = iEdgeDegB - 13;
  int iEdgeDiffC = iEdgeDegC - 13;
  int iEdgeDiffD = iEdgeDegD - 13;
  int iEdgeDiffE = iEdgeDegE - 13;

  int iScoreBefore = (iEdgeDiffA * iEdgeDiffA + iEdgeDiffB * iEdgeDiffB
      + iEdgeDiffC * iEdgeDiffC + iEdgeDiffD * iEdgeDiffD
      + iEdgeDiffE * iEdgeDiffE);
  int iScoreAfter = iScoreBefore;

  switch (FC.getFaceCat())
    {
    case FaceSwapInfo3D::eT32:
      // Removes an edge DE:
      iScoreAfter += -2 * iEdgeDiffD - 2 * iEdgeDiffE + 2;
      break;
    case FaceSwapInfo3D::eT23:
      // Adds an edge DE:
      iScoreAfter += 2 * iEdgeDiffD + 2 * iEdgeDiffE + 2;
      break;
    case FaceSwapInfo3D::eT22:
    case FaceSwapInfo3D::eT44:
    default:
      // Do nothing.
      break;
    }
  return ((iScoreAfter < iScoreBefore)
      || ((iScoreAfter == iScoreBefore)
	  && (FC.getFaceCat() == FaceSwapInfo3D::eT32)));
}

//@@ Determines whether swapping will improve the maximum face angle.
// Returns true to swap, false to remain the same.
bool
MinMaxDihedSwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
{
#ifdef SIM_ANNEAL_TEST
  return SwapDecider3D::doFaceSwap(FC);
#endif
  Vert* pVVertA = FC.getVertA();
  Vert* pVVertB = FC.getVertB();
  Vert* pVVertC = FC.getVertC();
  Vert* pVVertD = FC.getVertD();
  Vert* pVVertE = FC.getVertE();

  FaceSwapInfo3D::faceCat eFC = FC.getFaceCat();

  // This algorithm finds face angles between adjacent faces by dotting
  // their unit normals.  The largest magnitude loses.
  //
  // To prevent pairs from flopping back and forth, a swap is only
  // authorized if the inequality is bigger than eps.

  double dEps = 1.e-10;
  switch (eFC)
    {
    case FaceSwapInfo3D::eT23:
    case FaceSwapInfo3D::eT32:
      {
	//@@@ The case where the five pVVerts form a genuinely convex set
	double adNormABD[3], adNormABE[3], adNormADE[3];
	double adNormBCD[3], adNormBCE[3], adNormBDE[3];
	double adNormCAD[3], adNormCAE[3], adNormCDE[3];
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertD->getCoords(), adNormABD);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);
	calcUnitNormal(pVVertA->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormADE);

	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormBCE);
	calcUnitNormal(pVVertB->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormBDE);

	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertD->getCoords(), adNormCAD);
	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertE->getCoords(), adNormCAE);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);

	double dLastDot, dWorst2 = -1., dWorst3 = -1.;

	// Care must be taken to be certain that I'm taking the dot
	// product of two vectors that both point into the enclosed angle
	// or both point out, or correct the sign.

	dLastDot = dDOT3D(adNormABD, adNormBCD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormCAD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormCAD, adNormABD); /* Both in */
	dWorst2 = max(dWorst2, dLastDot);

	dLastDot = dDOT3D(adNormABE, adNormBCE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormBCE, adNormCAE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);
	dLastDot = dDOT3D(adNormCAE, adNormABE); /* Both out */
	dWorst2 = max(dWorst2, dLastDot);

	dLastDot = -dDOT3D(adNormABD, adNormABE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormBCD, adNormBCE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormCAD, adNormCAE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);

	dLastDot = -dDOT3D(adNormADE, adNormBDE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormBDE, adNormCDE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);
	dLastDot = -dDOT3D(adNormCDE, adNormADE); /* One in, one out */
	dWorst3 = max(dWorst3, dLastDot);

	if (dWorst3 > dWorst2 + dEps)
	  return (eFC == FaceSwapInfo3D::eT32);
	else if (dWorst2 > dWorst3 + dEps)
	  return (eFC == FaceSwapInfo3D::eT23);
	else
	  return (false);
      }
    case FaceSwapInfo3D::eT22:
      {
	double adNormACE[3], adNormECB[3], adNormBCD[3], adNormDCA[3];
	double adNormABC[3], adNormCDE[3], adNormABE[3];
	calcUnitNormal(pVVertA->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormACE);
	calcUnitNormal(pVVertE->getCoords(), pVVertC->getCoords(),
		       pVVertB->getCoords(), adNormECB);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertD->getCoords(), pVVertC->getCoords(),
		       pVVertA->getCoords(), adNormDCA);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);

	double dLastDot, dWorstNow = -1., dWorstSwapped = -1.;

	dLastDot = dDOT3D(adNormACE, adNormECB);
	dWorstNow = max(dWorstNow, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormDCA);
	dWorstNow = max(dWorstNow, dLastDot);
	dLastDot = fabs(dDOT3D(adNormABC, adNormABE));
	dWorstNow = max(dWorstNow, dLastDot);

	dLastDot = dDOT3D(adNormACE, adNormDCA);
	dWorstSwapped = max(dWorstSwapped, dLastDot);
	dLastDot = dDOT3D(adNormBCD, adNormECB);
	dWorstSwapped = max(dWorstSwapped, dLastDot);
	dLastDot = fabs(dDOT3D(adNormCDE, adNormABE));
	dWorstSwapped = max(dWorstSwapped, dLastDot);

	if (dWorstNow > dWorstSwapped + dEps)
	  return (true);
	else
	  return (false);
      }
    case FaceSwapInfo3D::eT44:
      return (false);
    default:
      assert2(0, "Bad face type in finding max dihedral angle");
      return (false);
    }
}

//@@ Determines whether swapping will improve the maximum face angle.
// Returns true to swap, false to remain the same.
bool
MaxMinSineSwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
{
#ifdef SIM_ANNEAL_TEST
  return SwapDecider3D::doFaceSwap(FC);
#endif
  Vert* pVVertA = FC.getVertA();
  Vert* pVVertB = FC.getVertB();
  Vert* pVVertC = FC.getVertC();
  Vert* pVVertD = FC.getVertD();
  Vert* pVVertE = FC.getVertE();

  FaceSwapInfo3D::faceCat eFC = FC.getFaceCat();

  // This algorithm finds sine of face angles between adjacent faces by
  // crossing their unit normals.  The idea is to pick the configuration
  // which maximizes this value.
  //
  // To prevent pairs from flopping back and forth, a swap is only
  // authorized if the inequality is bigger than eps.

  double dEps = 1.e-10;
  switch (eFC)
    {
    case FaceSwapInfo3D::eT23:
    case FaceSwapInfo3D::eT32:
      {

//    	double TestCoordA[3];
//    	double TestCoordB[3];
//    	double TestCoordC[3];
//    	double TestCoordD[3];
//    	double TestCoordE[3];
//
//    	dMetricCoords_3D(pVVertA,TestCoordA,pVVertB,TestCoordB,pVVertC,TestCoordC,pVVertD,TestCoordD,pVVertE,TestCoordE);

	double adNormABD[3], adNormABE[3], adNormADE[3], adNormABC[3];
	double adNormBCD[3], adNormBCE[3], adNormBDE[3];
	double adNormCAD[3], adNormCAE[3], adNormCDE[3];
	//@@@ The case where the five pVVerts form a genuinely convex set

//      calcUnitNormal(TestCoordA, TestCoordB,
//    		  TestCoordD, adNormABD);
//      calcUnitNormal(TestCoordA, TestCoordB,
//    		  TestCoordE, adNormABE);
//      calcUnitNormal(TestCoordA, TestCoordD,
//    		  TestCoordE, adNormADE);
//
//      calcUnitNormal(TestCoordA, TestCoordB,
//    		  TestCoordC, adNormABC); // Points to E
//
//      calcUnitNormal(TestCoordB, TestCoordC,
//    		  TestCoordD, adNormBCD);
//      calcUnitNormal(TestCoordB, TestCoordC,
//    		  TestCoordE, adNormBCE);
//      calcUnitNormal(TestCoordB, TestCoordD,
//    		  TestCoordE, adNormBDE);
//
//      calcUnitNormal(TestCoordC, TestCoordA,
//    		  TestCoordD, adNormCAD);
//      calcUnitNormal(TestCoordC, TestCoordA,
//    		  TestCoordE, adNormCAE);
//      calcUnitNormal(TestCoordC, TestCoordD,
//    		  TestCoordE, adNormCDE);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertD->getCoords(), adNormABD);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);
	calcUnitNormal(pVVertA->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormADE);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC); // Points to E

	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormBCE);
	calcUnitNormal(pVVertB->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormBDE);

	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertD->getCoords(), adNormCAD);
	calcUnitNormal(pVVertC->getCoords(), pVVertA->getCoords(),
		       pVVertE->getCoords(), adNormCAE);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);

	double dWorst2 = 1., dWorst3 = 1.;

	// Because I'm taking the MAGNITUDE of the cross product, the
	// actual sign of the result doesn't matter (vector can be
	// reversed).  This means that normal orientation doesn't matter.

	// Sines of dihedral angles of tet ABCD
	double adTemp[3], dMag;

	vCROSS3D(adNormABD, adNormBCD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCD, adNormCAD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAD, adNormABD, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormABD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	// Sines of dihedral angles of tet ABCE
	vCROSS3D(adNormABE, adNormBCE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCE, adNormCAE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAE, adNormABE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormABE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormBCE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	vCROSS3D(adNormCAE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst2 = min(dWorst2, dMag);

	// Now have the smallest sine of dihedral angle in tet ABCD and ABCE

	//Sines of dihedral angles of tet ABDE
	vCROSS3D(adNormABD, adNormABE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABD, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABD, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormABE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormADE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	//Sines of dihedral angles of tet BCDE
	vCROSS3D(adNormBCD, adNormBCE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCD, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCE, adNormBDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBCE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormBDE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	//Sines of dihedral angles of tet CADE
	vCROSS3D(adNormCAD, adNormCAE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAD, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCAE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);

	vCROSS3D(adNormCDE, adNormADE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorst3 = min(dWorst3, dMag);
	// Now have the smallest sine of dihedral in tets ABDE, BCDE, CADE

	if (dWorst3 > dWorst2 + dEps)
	  return (eFC == FaceSwapInfo3D::eT23);
	else if (dWorst2 > dWorst3 + dEps)
	  return (eFC == FaceSwapInfo3D::eT32);
	else
	  return (false);
      }
    case FaceSwapInfo3D::eT22:
      {

//    	double TestCoordA[3];
//    	double TestCoordB[3];
//    	double TestCoordC[3];
//    	double TestCoordD[3];
//    	double TestCoordE[3];
//
//    	dMetricCoords_3D(pVVertA,TestCoordA,pVVertB,TestCoordB,pVVertC,TestCoordC,pVVertD,TestCoordD,pVVertE,TestCoordE);

	double adNormACE[3], adNormECB[3], adNormBCD[3], adNormDCA[3];
	double adNormABC[3], adNormCDE[3], adNormABE[3];
//      calcUnitNormal(TestCoordA, TestCoordC,
//    		  TestCoordE, adNormACE);
//      calcUnitNormal(pVVertE->getCoords(), TestCoordC,
//    		  TestCoordB, adNormECB);
//      calcUnitNormal(TestCoordB, TestCoordC,
//    		  TestCoordD, adNormBCD);
//      calcUnitNormal(TestCoordD, TestCoordC,
//    		  TestCoordA, adNormDCA);
//
//      calcUnitNormal(TestCoordA, TestCoordB,
//    		  TestCoordC, adNormABC);
//      calcUnitNormal(TestCoordC, TestCoordD,
//    		  TestCoordE, adNormCDE);
//      calcUnitNormal(TestCoordA, TestCoordB,
//    		  TestCoordE, adNormABE);

	calcUnitNormal(pVVertA->getCoords(), pVVertC->getCoords(),
		       pVVertE->getCoords(), adNormACE);
	calcUnitNormal(pVVertE->getCoords(), pVVertC->getCoords(),
		       pVVertB->getCoords(), adNormECB);
	calcUnitNormal(pVVertB->getCoords(), pVVertC->getCoords(),
		       pVVertD->getCoords(), adNormBCD);
	calcUnitNormal(pVVertD->getCoords(), pVVertC->getCoords(),
		       pVVertA->getCoords(), adNormDCA);

	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertC->getCoords(), adNormABC);
	calcUnitNormal(pVVertC->getCoords(), pVVertD->getCoords(),
		       pVVertE->getCoords(), adNormCDE);
	calcUnitNormal(pVVertA->getCoords(), pVVertB->getCoords(),
		       pVVertE->getCoords(), adNormABE);

	double dWorstNow = 1., dWorstSwap = 1.;
	double adTemp[3], dMag;

	// Sines of dihedral angles of tet ABCD, except those that are the
	// same both ways
	vCROSS3D(adNormBCD, adNormDCA, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormBCD, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormDCA, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Sines of dihedral angles of tet ABCE, except those that are the
	// same both ways
	vCROSS3D(adNormECB, adNormACE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormECB, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	vCROSS3D(adNormACE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Only check the center face once; the sine is the same on both
	// sides
	vCROSS3D(adNormABE, adNormABC, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstNow = min(dWorstNow, dMag);

	// Now have the smallest sine of dihedral angle in tet ABCD and ABCE

	// Sines of dihedral angles of tet DECA, except those that are the
	// same both ways
	vCROSS3D(adNormDCA, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormACE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormDCA, adNormACE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Sines of dihedral angles of tet DECB, except those that are the
	// same both ways
	vCROSS3D(adNormBCD, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormECB, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	vCROSS3D(adNormBCD, adNormECB, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Only check the center face once; the sine is the same on both
	// sides
	vCROSS3D(adNormABE, adNormCDE, adTemp);
	dMag = dMAG3D(adTemp);
	dWorstSwap = min(dWorstSwap, dMag);

	// Now have the smallest sine of dihedral angle in tet DECA and DECB

	if (dWorstNow + dEps < dWorstSwap)
	  return (true);
	else
	  return (false);
      }
    case FaceSwapInfo3D::eT44:
      return (false);
    default:
      assert2(0, "Bad face type in finding max dihedral angle");
      return (false);
    }
}

bool
MetricSwapDecider3D::doFaceSwap(const FaceSwapInfo3D& FC) const
{

  Vert* pVA = FC.getVertA();
  Vert* pVB = FC.getVertB();
  Vert* pVC = FC.getVertC();
  Vert* pVD = FC.getVertD();
  Vert* pVE = FC.getVertE();

  assert(pVA->isValid());
  assert(pVB->isValid());
  assert(pVC->isValid());
  assert(pVD->isValid());
  assert(pVE->isValid());

  FaceSwapInfo3D::faceCat eFC = FC.getFaceCat();
  m_pQM->vEvalInMetric(true);
  // If these values never get changed, there'll never be any face swaps.
  double MinqualBefore = 1, MinqualAfter = 0;
  double MaxqualBefore = 1, MaxqualAfter = 0;
  switch (eFC)
    {
    case FaceSwapInfo3D::eT23:
      {
	// Old configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	MinqualBefore = min(qualACBD, qualABCE);
	MaxqualBefore = max(qualACBD, qualABCE);

	// New configuration
	double qualDEAB = m_pQM->eval(pVD, pVE, pVA, pVB);
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	MinqualAfter = min(qualDEAB, min(qualDEBC, qualDECA));
	MaxqualAfter = max(qualDEAB, max(qualDEBC, qualDECA));

	break;
      }

    case FaceSwapInfo3D::eT32:
      {
	// Old configuration
	double qualDEAB = m_pQM->eval(pVD, pVE, pVA, pVB);
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	MinqualBefore = min(qualDEAB, min(qualDEBC, qualDECA));
	MaxqualBefore = max(qualDEAB, max(qualDEBC, qualDECA));

	// New configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	MinqualAfter = min(qualACBD, qualABCE);
	MaxqualAfter = max(qualACBD, qualABCE);

	break;
      }

    case FaceSwapInfo3D::eT22:
      {
	// Can't swap T22 if the bdry is precious.
	if (!m_allowBdryChanges)
	  return false;
	// Old configuration
	double qualACBD = m_pQM->eval(pVA, pVC, pVB, pVD);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	MinqualBefore = min(qualACBD, qualABCE);
	MaxqualBefore = max(qualACBD, qualABCE);

	// New configuration
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	MinqualAfter = min(qualDEBC, qualDECA);
	MaxqualAfter = max(qualDEBC, qualDECA);

	break;
      }
    case FaceSwapInfo3D::eT44:
      {
	Vert *pVF = FC.getVertF();
	VALIDATE_INPUT(pVF->isValid());
	// Old configuration
	double qualABDC = m_pQM->eval(pVA, pVB, pVD, pVC);
	double qualABCE = m_pQM->eval(pVA, pVB, pVC, pVE);
	double qualABEF = m_pQM->eval(pVA, pVB, pVE, pVF);
	double qualABFD = m_pQM->eval(pVA, pVB, pVF, pVD);
	MinqualBefore = min(min(qualABDC, qualABCE), min(qualABEF, qualABFD));
	MaxqualBefore = max(max(qualABDC, qualABCE), max(qualABEF, qualABFD));

	// New configuration
	double qualDEBC = m_pQM->eval(pVD, pVE, pVB, pVC);
	double qualDECA = m_pQM->eval(pVD, pVE, pVC, pVA);
	double qualDEAF = m_pQM->eval(pVD, pVE, pVA, pVF);
	double qualDEFB = m_pQM->eval(pVD, pVE, pVF, pVB);
	MinqualAfter = min(min(qualDEBC, qualDECA), min(qualDEAF, qualDEFB));
	MaxqualAfter = max(max(qualDEBC, qualDECA), max(qualDEAF, qualDEFB));

	break;
      }

    default:
      {
	assert2(0, "Bad Face Cat for swapping");
	return false;
      }
    }

  bool bSwap = false;

  if (m_pQM->getName() == "VolumeLength3D") {
    bSwap = (MaxqualAfter < MaxqualBefore);
  }

  else if (m_pQM->getName() == "sine of dihedral angles") {
    bSwap = (MinqualAfter > MinqualBefore);
  }

  else {
    bSwap = (MinqualAfter > MinqualBefore);
  }

  return bSwap;
}

double
SwapDecider3D::calcOrigQuality(const EdgeSwapInfo& ES) const
{
  double origQual = DBL_MAX;

  // Compute the quality of the current configuration.  This is the
  // smallest quality for any of the existing cells.
  for (int i = ES.getNumOrigCells() - 1; i >= 0; --i) {
    Cell *pC = ES.getOrigCell(i);
    // Evaluate quality, depending on the swap criterion in use.
    assert(
	1
	    == checkOrient3D(pC->getVert(0), pC->getVert(1), pC->getVert(2),
			     pC->getVert(3)));

    double qual = m_pQM->eval(pC->getVert(0), pC->getVert(1), pC->getVert(2),
			      pC->getVert(3));
    origQual = std::min(qual, origQual);
  }
  return origQual;
}
// experiment in recovering edges. qWant = true if I want the edge,
// qWant = false if i don't want the edge
bool
SwapDecider3D::doEdgeSwap(EdgeSwapInfo& ES, Vert *pVWant0, Vert *pVWant1,
			  bool qWant) const
{
#ifdef SIM_ANNEAL_TEST
  return false;
#endif

  if ((!m_pQM && !pVWant0 && !pVWant1) || !ES.isSwappableEdge()) {
    // If there's no cell-by-cell quality measure, call the whole thing off.
    // unless we have a target edge, in which case, just give bonuses to those configs
    // with the edge
    return false;
  }
  // If there's a particular edge that's desired, to force it into the
  // mesh, then all configurations with this edge get a bonus that
  // will make it so they're better than the original config, while
  // configs without the new edge get a penalty that ensures they'll
  // be worse.  Since we're looking for the best config that has the
  // given edge, the original quality is irrelevant and that
  // calculation gets skipped.
  bool qEdgePresent = false, qBonusPossible = (pVWant0->isValid()
      && pVWant1->isValid());
  double origQual = qBonusPossible ? 0 : calcOrigQuality(ES);
  if (!qWant)
    origQual = -5000;
  int numOrigVerts = ES.getNumOrigVerts();
  if (numOrigVerts < 3)
    return false;
  int numOrigCells = ES.getNumOrigCells();
  const double invalidPenalty = -10000;
  const double edgeBonus = 1000;
  bool qReturnVal = false;

  Vert *northVert = ES.getNorthVert();
  Vert *southVert = ES.getSouthVert();

  // Compute quality for each possible tet pair (one pair per possible
  // triangle). 
  double a3dNewQual[10][10][10];
  int iSign = checkOrient3D(ES.getVert(0), ES.getVert(1), southVert, northVert);
  int iV0, iV1, iV2;
#ifndef NDEBUG
  for (iV0 = 0; iV0 < 10; iV0++)
    for (iV1 = 0; iV1 < 10; iV1++)
      for (iV2 = 0; iV2 < 10; iV2++)
	a3dNewQual[iV0][iV1][iV2] = 0;
  for (int iDum = 0; iDum < numOrigCells; iDum++)
    assert(
	iSign
	    == checkOrient3D(ES.getVert(iDum),
			     ES.getVert((iDum + 1) % numOrigVerts), southVert,
			     northVert));
#endif
  for (iV0 = 0; iV0 < numOrigVerts; iV0++) {
    Vert *pV0 = ES.getVert(iV0);
    for (iV1 = iV0 + 1; iV1 < numOrigVerts; iV1++) {
      Vert *pV1 = ES.getVert(iV1);
      for (iV2 = iV1 + 1; iV2 < numOrigVerts; iV2++) {
	Vert *pV2 = ES.getVert(iV2);
	if (checkOrient3D(pV0, pV1, pV2, northVert) == iSign)
	  // Evaluate quality
	  if (m_pQM)
	    a3dNewQual[iV0][iV1][iV2] = m_pQM->eval(pV0, pV1, pV2, northVert);
	  else
	    a3dNewQual[iV0][iV1][iV2] = 1.0; // set to 1 so that it checks the other side
	else
	  a3dNewQual[iV0][iV1][iV2] = invalidPenalty;

	// Don't bother evaluating for the other tet if this one is
	// worse that the worst in the original config, because this tet
	// pair can never be in a final config anyway.
	if (a3dNewQual[iV0][iV1][iV2] >= origQual) {
	  double dDummyQual = invalidPenalty;
	  if (checkOrient3D(pV0, pV2, pV1, southVert) == iSign) {
	    if (m_pQM)
	      dDummyQual = m_pQM->eval(pV0, pV2, pV1, southVert);
	    else
	      dDummyQual = 1.0;
	  }
	  a3dNewQual[iV0][iV1][iV2] = std::min(a3dNewQual[iV0][iV1][iV2],
					       dDummyQual);

	}

	// Copy things around symmetrically
	a3dNewQual[iV2][iV1][iV0] = a3dNewQual[iV2][iV0][iV1] =
	    a3dNewQual[iV1][iV2][iV0] = a3dNewQual[iV1][iV0][iV2] =
		a3dNewQual[iV0][iV2][iV1] = a3dNewQual[iV0][iV1][iV2];

      } // Loop over third vert
    } // Loop over second vert
  } // Loop over first vert

    // Find the new configuration with the highest quality.
  double dMaxQual = origQual;
  int iCanon, iPerm;
  const ConfSet& CS = EdgeConfigs::getInstance().getConfSet(numOrigVerts);
  for (iCanon = 0; iCanon < CS.getNumCanon(); iCanon++) {
    const EdgeConfig& Conf = CS.getConfig(iCanon);
    int iNPerm = Conf.getNumPerm();
    for (iPerm = 0; iPerm < iNPerm; iPerm++) {
      double dConfQual = 1000;
      // Once this config is known to be worse than the best one we know
      // about, we can give up on it.
      for (int iTri = 0;
	  iTri < numOrigVerts - 2 && (dConfQual >= dMaxQual || qBonusPossible);
	  iTri++) {
	iV0 = (Conf.getFaceVert(iTri, 0) + iPerm) % numOrigVerts;
	iV1 = (Conf.getFaceVert(iTri, 1) + iPerm) % numOrigVerts;
	iV2 = (Conf.getFaceVert(iTri, 2) + iPerm) % numOrigVerts;
	// Give a bonus if a particular edge is desired and that edge
	// is included in the configuration.
	if (qBonusPossible && !qEdgePresent) {
	  Vert * pV0 = ES.getVert(iV0);
	  Vert * pV1 = ES.getVert(iV1);
	  Vert * pV2 = ES.getVert(iV2);
	  if ((pVWant0 == pV0 || pVWant0 == pV1 || pVWant0 == pV2)
	      && (pVWant1 == pV0 || pVWant1 == pV1 || pVWant1 == pV2))
	    qEdgePresent = true;
	}
	dConfQual = min(dConfQual, a3dNewQual[iV0][iV1][iV2]);
      }
      // Only give a bonus to valid configurations
      if (qBonusPossible)
	dConfQual += qEdgePresent ? edgeBonus : -edgeBonus;
      if (dConfQual > dMaxQual) {
	ES.setNewConfig(iCanon, iPerm);
	qReturnVal = true;
      }
    } // Loop over permutations
  } // Loop over canonical configurations
  return qReturnVal;
}
// For now this ONLY checks for delaunay criteria before swapping
// where Delaunay is based on being in the circumsphere of the face
// some assumption of planarity
bool
SwapDecider3D::doBdryEdgeSwap(EdgeSwapInfo& ES) const
{
  // If strict patch checking is on, the following call checks that
  // the two bfaces are on the same patch.  Otherwise, it uses a
  // weaker check:  are the bdry conditions the same and the faces
  // nearly coplanar. 

  BFace *bdryFace0 = ES.getBFace(0);
  BFace *bdryFace1 = ES.getBFace(1);
  if (!areBdryFacesOKToSwap(bdryFace0, bdryFace1)) {
    return false;
  }
  // Do the two boundary triangles form a convex quad?  If so, the
  // orientation of the tets formed by these triangles and any
  // non-coplanar point (i.e., any other member of m_outerVerts) will be the
  // same.  This exit is redundant (no swapping would occur for this
  // case anyway), but saves a lot of quality computations.
  unsigned numOrigVerts = ES.getNumOrigVerts();

  // Grab the first and last verts on the "equator".
  Vert *firstVert = ES.getVert(0);
  Vert *lastVert = ES.getVert(numOrigVerts - 1);
  Vert *northVert = ES.getNorthVert();
  Vert *southVert = ES.getSouthVert();
  Vert *interVert = ES.getVert(1);
//  bdryFace0->printCellInfo();
//  bdryFace1->printCellInfo();
//  printf("original verts\n");
//  for(int i = 0; i < numOrigVerts; i++)
//  	ES.getVert(i)->printVertInfo();
  if (checkOrient3D(firstVert, lastVert, northVert, interVert)
      != checkOrient3D(firstVert, southVert, lastVert, interVert)) {
    return false;
  }
  if (strcmp(m_name.c_str(), "Delaunay") == 0) {
    // check delaunay property of existing orientation
    double adCircCent[3];
    bdryFace0->calcCircumcenter(adCircCent);

    double dRadius = dDIST3D(adCircCent, bdryFace0->getVert(0)->getCoords());
    double dDistance = dDIST3D(adCircCent, lastVert->getCoords());

    switch (iFuzzyComp(dDistance, dRadius))
      {
      case -1:
	// Distance smaller than radius -> inside
	return doEdgeSwap(ES, ES.getNorthVert(), ES.getSouthVert(), false);
	break;
      case 1:
	// Distance larger than radius -> outside
	return false;
	break;
      case 0:
	return doEdgeSwap(ES, ES.getNorthVert(), ES.getSouthVert(), false);
	break;
      default:
	assert(0);
	return true;
      }
  }
  return doEdgeSwap(ES);
}

bool
SwapDecider3D::areBdryFacesOKToSwap(const BFace * const pBF0,
				    const BFace * const pBF1) const
{
  //     if (qStrictPatchChecking) {
  //       return (pBF0->pPatchPointer() == pBF1->pPatchPointer());
  //     }
  //     else {
  // Check that the BC's are the same and that the normals are more or
  // less aligned (within the tolerances of dMaxAngleForSurfSwap).

  if (pBF0->getBdryCondition() == pBF1->getBdryCondition()) {
    double adNorm0[3], adNorm1[3];
    pBF0->calcUnitNormal(adNorm0);
    pBF1->calcUnitNormal(adNorm1);
    NORMALIZE3D(adNorm0);
    NORMALIZE3D(adNorm1);
    double dDot = dDOT3D(adNorm0, adNorm1);
    return (dDot >= cos(m_maxBdryAngle * M_PI / 180));
  }
  else {
    return false;
  }
  //     }
}

double
SwapDecider3D::evalTetQual(const Cell* const pC) const
{
  VALIDATE_INPUT(pC->isValid() && !pC->isDeleted());
  VALIDATE_INPUT(pC->getType() == CellSkel::eTet);
  if (!m_pQM) {
    return 0;
  }
  return m_pQM->eval(pC->getVert(0), pC->getVert(1), pC->getVert(2),
		     pC->getVert(3));
}
