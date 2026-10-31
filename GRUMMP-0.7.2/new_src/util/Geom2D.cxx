#include <math.h>
#include "GR_misc.h"
#include "GR_Geometry.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_AdaptPred.h"
#include "GR_Vec.h"

//@ Compute the vector normal for a set of three verts
void
calcNormal2D(const double adA[2], const double adB[2], const double adC[2],
	     double * const pdRes)
{
  double adTempB[] =
    { adB[0] - adA[0], adB[1] - adA[1] };
  double adTempC[] =
    { adC[0] - adA[0], adC[1] - adA[1] };
  *pdRes = adTempB[0] * adTempC[1] - adTempB[1] * adTempC[0];
}

//@ Compute the vector normal for a set of two verts
// This is always rotated 90 degrees counterclockwise from the vector
// running from vertex A to vertex B.
void
calcNormal2D(const double adA[2], const double adB[2], double adRes[2])
{
  adRes[0] = adA[1] - adB[1];
  adRes[1] = adB[0] - adA[0];
}

//@ Determine whether a vertex lies inside, on or outside a triangle
int
isInTri(const Vert * const pV, const Vert * const pVertA,
	const Vert * const pVertB, const Vert * const pVertC)
// For pV to be inside the tri formed by pVertA-C, the orientation
// has to be right wrt every face.  The return value is -1 if the
// vert lies outside the tri, 0 if it's colinear with all three
// faces (impossible for valid tris), 1 if it coincides with a vert
// in the tri, 2 if it lies on an edge, and 3 if it lies strictly
// inside.

{
  assert(pV->isValid() && pV->getSpaceDimen() == 2);
  assert(pVertA->isValid() && pVertA->getSpaceDimen() == 2);
  assert(pVertB->isValid() && pVertB->getSpaceDimen() == 2);
  assert(pVertC->isValid() && pVertC->getSpaceDimen() == 2);

  if (checkOrient2D(pVertA, pVertB, pVertC) == (-1)) {
    //The vertices are clockwise and need to be changed
    return isInTri(pV, pVertC, pVertB, pVertA);
  }
  //    assert(iOrient2D(pVertA,pVertB,pVertC)==1);

  int iRetVal = 0;
  int iCheck = checkOrient2D(pVertA, pVertB, pV);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient2D(pVertB, pVertC, pV);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient2D(pVertC, pVertA, pV);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  if (iRetVal == 0)
    printf(
	"(%.3e,%.3e) is coincident to all three faces of (%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n",
	pV->x(), pV->y(), pVertA->x(), pVertA->y(), pVertB->x(), pVertB->y(),
	pVertC->x(), pVertC->y());

  return (iRetVal);
}

int
isInTri(const double adVert[2], const double adVertA[2],
	const double adVertB[2], const double adVertC[2])
// For Vert to be inside the tri formed by pVertA-C, the orientation
// has to be right wrt every face.  The return value is -1 if the
// vert lies outside the tri, 0 if it's colinear with all three
// faces (impossible for valid tris), 1 if it coincides with a vert
// in the tri, 2 if it lies on an edge, and 3 if it lies strictly
// inside.

{
  if (checkOrient2D(adVertA, adVertB, adVertC) == (-1)) {
    //The vertices are clockwise and need to be changed
    return isInTri(adVert, adVertC, adVertB, adVertA);
  }
  assert(checkOrient2D(adVertA, adVertB, adVertC) == 1);

  int iRetVal = 0;
  int iCheck = checkOrient2D(adVertA, adVertB, adVert);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient2D(adVertB, adVertC, adVert);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient2D(adVertC, adVertA, adVert);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  return (iRetVal);
}

//@ 2D orientation test
int
checkOrient2D(const Vert * const pVertA, const Vert * const pVertB,
	      const Vert * const pVertC)
// Computes the orientation of three verts in 2D
{
  assert(pVertA->isValid());
  assert(pVertB->isValid());
  assert(pVertC->isValid());
  assert(pVertA->getSpaceDimen() == 2);
  assert(pVertB->getSpaceDimen() == 2);
  assert(pVertC->getSpaceDimen() == 2);
  return checkOrient2D(pVertA->getCoords(), pVertB->getCoords(),
		       pVertC->getCoords());
}

int
checkOrient2D(const double adA[2], const double adB[2], const double adC[2])
{
  if (qAdaptPred) {
    double dResult = orient2d_shew(adA, adB, adC);
    if (dResult > 0)
      return 1;
    else if (dResult < 0)
      return -1;
    else
      return 0;
  }
  else {
    double xa = adA[0];
    double ya = adA[1];

    double xb = adB[0];
    double yb = adB[1];

    double xc = adC[0];
    double yc = adC[1];

    // Currently no exact arithmetic
    double dDet = ((xb - xa) * (yc - ya) - (xc - xa) * (yb - ya));
    // Scale the determinant by the mean of the magnitudes of vectors:
    double dScale = (dDIST2D(adA, adB) + dDIST2D(adB, adC) + dDIST2D(adC, adA))
	/ 3;
    dDet /= (dScale * dScale);
    if (dDet > 1.e-10)
      return 1;
    else if (dDet < -1.e-10)
      return -1;
    else
      return 0;
  }
}

//@ Incircle test
int
isInCircle(const Vert * const pVertA, const Vert * const pVertB,
	   const Vert * const pVertC, const Vert * const pVertD)
// Determines whether pVertD is inside or outside the circumcircle
// of the tri formed by pVertA-C, using as nearly exact arithmetic
// as required.  Returns 1 if pVertD is inside the
// circumcircle of the tri, -1 if outside, 0 if on.
{
  assert(pVertA->isValid());
  assert(pVertB->isValid());
  assert(pVertC->isValid());
  assert(pVertD->isValid());
  return isInCircle(pVertA->getCoords(), pVertB->getCoords(),
		    pVertC->getCoords(), pVertD->getCoords());
}

int
isInCircle(const double adA[2], const double adB[2], const double adC[2],
	   const double adD[2])
// Determines whether pVertD is inside or outside the circumcircle
// of the tri formed by pVertA-C, using as nearly exact arithmetic
// as required.  Returns 1 if pVertD is inside the
// circumcircle of the tri, -1 if outside, 0 if on.
{
  if (qAdaptPred) {
    double dResult = incircle_shew(adA, adB, adC, adD);
    if (dResult > 0)
      return 1;
    else if (dResult < 0)
      return -1;
    else
      return 0;
  }
  else {
    double xa = adA[0];
    double ya = adA[1];
    double wa2 = xa * xa + ya * ya;

    double xb = adB[0];
    double yb = adB[1];
    double wb2 = xb * xb + yb * yb;

    double xc = adC[0];
    double yc = adC[1];
    double wc2 = xc * xc + yc * yc;

    double xd = adD[0];
    double yd = adD[1];
    double wd2 = xd * xd + yd * yd;

    // Currently no exact arithmetic
    double dDet1 = ((xb - xa)
	* ((yc - ya) * (wd2 - wa2) - (yd - ya) * (wc2 - wa2))
	+ (xc - xa) * ((yd - ya) * (wb2 - wa2) - (yb - ya) * (wd2 - wa2))
	+ (xd - xa) * ((yb - ya) * (wc2 - wa2) - (yc - ya) * (wb2 - wa2)));
    int iOrient = checkOrient2D(adA, adB, adC);
    double dProd = dDet1 * iOrient;
    // Scale by the mean of the magnitudes of vectors:
    double dScale = (dDIST2D(adB, adA) + dDIST2D(adC, adA) + dDIST2D(adB, adC))
	/ 3;
    // Unclear exactly how to scale this...
    dProd /= (dScale * dScale * dScale * dScale);
    if (dProd > 1.e-10)
      return -1;
    else if (dProd < -1.e-10)
      return 1;
    else
      return 0;
  }
}

bool
calcLineIntersection(const Vert* const pVA, const Vert* const pVB,
		     const Vert* const pVC, const Vert* const pVD,
		     double intersection[2])
{

  assert(pVA->isValid() && pVB->isValid() && pVC->isValid() && pVD->isValid());
  assert(
      pVA->getSpaceDimen() == 2 && pVB->getSpaceDimen() == 2
	  && pVC->getSpaceDimen() == 2 && pVD->getSpaceDimen() == 2);

  double det_x, det_y, det_denom, det1, det2;
  double x1 = pVA->x(), x2 = pVB->x(), x3 = pVC->x(), x4 = pVD->x();
  double y1 = pVA->y(), y2 = pVB->y(), y3 = pVC->y(), y4 = pVD->y();

  det1 = x1 * y2 - x2 * y1;
  det2 = x3 * y4 - x4 * y3;
  det_x = det1 * (x3 - x4) - det2 * (x1 - x2);
  det_y = det1 * (y3 - y4) - det2 * (y1 - y2);
  det_denom = (x1 - x2) * (y3 - y4) - (x3 - x4) * (y1 - y2);

  if (iFuzzyComp(det_denom, 0.) == 0)
    return false; //parallel
  else {
    intersection[0] = det_x / det_denom;
    intersection[1] = det_y / det_denom;
    return true;
  }

}

//int intersectRays(const Vert* const pVA, const Vert* pVB,
//		const Vert* const pVC, const Vert* pVD,
//		double* intersection) {
//
//	double inter[2];
//	if( !calcLineIntersection(pVA, pVB, pVC, pVD, inter) ) return -1;
//
//	double x = inter[0], y = inter[1];
//
//	if(intersection != NULL)
//	{ intersection[0] = x; intersection[1] = y; }
//
//	double x1 = pVC->x(), x2 = pVD->x();
//	double y1 = pVC->y(), y2 = pVD->y();
//
//	//   printf("intersection: %lf %lf\n", x, y);
//	//   printf("x1 = %lf, y1 = %lf\n", x1, y1);
//	//   printf("x2 = %lf, y2 = %lf\n", x2, y2);
//
//	double vec_vAvB[]    = { pVB->x() - pVA->x(), pVB->y() - pVA->y() };
//	double vec_vAinter[] = { inter[0]  - pVA->x(), inter[1]  - pVA->y() };
//	if( iFuzzyComp(dDOT2D(vec_vAvB, vec_vAinter), 0.) <= 0) return -1;
//
//	if( iFuzzyComp(std::min(x1, x2), x) <= 0 &&
//			iFuzzyComp(std::max(x1, x2), x) >= 0 &&
//			iFuzzyComp(std::min(y1, y2), y) <= 0 &&
//			iFuzzyComp(std::max(y1, y2), y) >= 0 ) return 1;
//
//	//   if( iFuzzyComp(std::min(x1, x2), x) == 0 &&
//	//       iFuzzyComp(std::max(x1, x2), x) == 0 &&
//	//       iFuzzyComp(std::min(y1, y2), y) == 0 &&
//	//       iFuzzyComp(std::max(y1, y2), y) == 0 ) return 0;
//
//	return -1;
//
//}

int
intersectSegments(const Vert* const pVA, const Vert* pVB, const Vert* const pVC,
		  const Vert* pVD)
{

  //The ends of both lines are included in intersection search.

  double inter[2];
  if (!calcLineIntersection(pVA, pVB, pVC, pVD, inter))
    return -1;

  double x = inter[0], y = inter[1];
  double x1 = pVA->x(), x2 = pVB->x(), x3 = pVC->x(), x4 = pVD->x();
  double y1 = pVA->y(), y2 = pVB->y(), y3 = pVC->y(), y4 = pVD->y();

  if (iFuzzyComp(std::min(x1, x2), x) <= 0
      && iFuzzyComp(std::max(x1, x2), x) >= 0
      && iFuzzyComp(std::min(y1, y2), y) <= 0
      && iFuzzyComp(std::max(y1, y2), y) >= 0
      && iFuzzyComp(std::min(x3, x4), x) <= 0
      && iFuzzyComp(std::max(x3, x4), x) >= 0
      && iFuzzyComp(std::min(y3, y4), y) <= 0
      && iFuzzyComp(std::max(y3, y4), y) >= 0)
    return 1;

  return -1;

}

double
dAreaTri(const Vert* const pV0, const Vert* pV1, const Vert* const pV2)
{
  double e0 = dDIST3D(pV0->getCoords(), pV1->getCoords());
  double e1 = dDIST3D(pV1->getCoords(), pV2->getCoords());
  double e2 = dDIST3D(pV2->getCoords(), pV0->getCoords());
  return dAreaTriFromEdgeLengths(e0, e1, e2);
}

double
dAreaTriFromEdgeLengths(double e0, double e1, double e2)
{
  // This specific ordering is required for numerical stability.
  if (e0 < e1) {
    std::swap(e0, e1);
  }
  // Now e0 > e1
  if (e1 < e2) {
    std::swap(e1, e2);
    // Now e1 > e2, but we don't know about e0 vs the new e1
    if (e0 < e1) {
      std::swap(e0, e1);
    }
  }
  assert(e0 >= e1);
  assert(e1 >= e2);
  // The parentheses enforce a numerically stable order of operations.
  double Area = 0.25
      * sqrt(
	  (e0 + (e1 + e2)) * (e2 - (e0 - e1)) * (e2 + (e0 - e1))
	      * (e0 + (e1 - e2)));
  return Area;
}

//void projectTo2D(const Vert* const pVertA, const Vert* const pVertB,	const Vert* const pVertC,
//		double * pA, double * pB, double * pC, double * pD){
//
//	projectTo2D(pVertA->getCoords(),pVertB->getCoords(), pVertC->getCoords(),pA,pB,pC,pD);
//
//}
//void projectTo2D(const double adA[3],const double adB[3], const double adC[3],
//		double * pA, double * pB, double * pC, double * pD){
//
//	double adNorm[3];
//	calcNormal3D(adA,adB,adC,adNorm);
//	NORMALIZE3D(adNorm);
//	//	printf("normal %f %f %f\n",adNorm[0],adNorm[1],adNorm[2]);
//	double adBasisX[3];
//	for(int i = 0; i<3; i++)
//		adBasisX[i] = adA[i]-adB[i];
//	NORMALIZE3D(adBasisX);
//	//	printf("X Basis %f %f %f\n",adBasisX[0],adBasisX[1], adBasisX[2]);
//	double adBasisY[3];
//	vCROSS3D(adNorm, adBasisX, adBasisY);
//
//	//	printf("Y Basis %f %f %f\n",adBasisY[0],adBasisY[1], adBasisY[2]);
//
//	pA[0] = dDOT3D(adA, adBasisX);
//	pA[1] = dDOT3D(adA, adBasisY);
//	pA[2] = 0.;
//	pB[0] = dDOT3D(adB, adBasisX);
//	pB[1] = dDOT3D(adB, adBasisY);
//	pB[2] = 0.;
//
//	pC[0] = dDOT3D(adC, adBasisX);
//	pC[1] = dDOT3D(adC, adBasisY);
//	pC[2] = 0.;
//
//	if(pD){
//		double pDTemp[2];
//		pD[2] = 0.;
//		pDTemp[0] = dDOT3D(pD, adBasisX);
//		pDTemp[1] = dDOT3D(pD, adBasisY);
//		pD[0] = pDTemp[0];
//		pD[1] = pDTemp[1];
//	}
//}
void
calcCovarianceMatrix(const Vert* const pV1, const Vert* pV2,
		     const Vert* const pV3, double * adMoment)
{

  // lets do this explicitly for now
  double x1 = pV1->x();
  double y1 = pV1->y();

  double x2 = pV2->x();
  double y2 = pV2->y();

  double x3 = pV3->x();
  double y3 = pV3->y();

  double areaConst = dAreaTri(pV1, pV2, pV3) / 18.;

  // using the algorithm from INRIA's principal component analysis for this
  adMoment[0] = areaConst
      * (x1 * x1 - x1 * x2 - x1 * x3 + x2 * x2 - x2 * x3 + x3 * x3);
  adMoment[1] = areaConst
      * (x1 * y1 - (x1 * y2) / 2 - (x2 * y1) / 2 - (x1 * y3) / 2 + x2 * y2
	  - (x3 * y1) / 2 - (x2 * y3) / 2 - (x3 * y2) / 2 + x3 * y3);
  adMoment[2] = areaConst
      * (y1 * y1 - y1 * y2 - y1 * y3 + y2 * y2 - y2 * y3 + y3 * y3);

  //

}

bool
calcSegmentCircleIntersections(const double adA[2], const double adB[2],
			       const double adCent[2], const double dRadius2,
			       double adIntersect[])
{
  // we define a point on the line by p = a+(b-a)*t and the circle by (x-xc)^2 + (y-yc)^2 - r^2 = 0
  double a = dDIST_SQ_2D(adA, adB);
  double b = 2.0
      * ((adB[0] - adA[0]) * (adA[0] - adCent[0])
	  + (adB[1] - adA[1]) * (adA[1] - adCent[1]));
  double c = (adA[0] - adCent[0]) * (adA[0] - adCent[0])
      + (adA[1] - adCent[1]) * (adA[1] - adCent[1]) - dRadius2;
  int iNRoots;
  double dt1, dt2;
  solveQuadratic(a, b, c, &iNRoots, &dt1, &dt2);

  // should always be this case
  // should be impossible for both intersection points to be on the line, as one endpoint
  // should be inside the circle
//	printf("%d roots\n",iNRoots);
//	printf("%f %f | %f %f\n",adA[0],adA[1],adB[0],adB[1]);
//	printf("dR2 %f | %f %f\n",dRadius2,adCent[0],adCent[1]);
//	printf("iNRoots %d %f %f\n",iNRoots,dt1,dt2);
  if (iNRoots == 2) {
//		printf("t1 t2 %f %f\n",dt1,dt2);
    // extra condition at end handles bounding case, which is technically false if one point lies
    // on boundary, and other outside
    if (iFuzzyComp(dt1, 0.0) > -1 && iFuzzyComp(dt1, 1.0) < 1) {
      adIntersect[0] = adA[0] + (adB[0] - adA[0]) * dt1;
      adIntersect[1] = adA[1] + (adB[1] - adA[1]) * dt1;
//			printf("intersect %f %f\n",adIntersect[0],adIntersect[1]);

      return true;
    }
    else if (iFuzzyComp(dt2, 0.0) > -1 && iFuzzyComp(dt2, 1.0) < 1) {
      adIntersect[0] = adA[0] + (adB[0] - adA[0]) * dt2;
      adIntersect[1] = adA[1] + (adB[1] - adA[1]) * dt2;
//			printf("intersect %f %f\n",adIntersect[0],adIntersect[1]);
      return true;
    }
    return false;
  }
  else if (iNRoots == 1) {

    if (iFuzzyComp(dt1, 0.0) > -1 && iFuzzyComp(dt1, 1.0) < 0) {
      adIntersect[0] = adA[0] + (adB[0] - adA[0]) * dt1;
      adIntersect[0] = adA[1] + (adB[1] - adA[1]) * dt1;
      return true;
    }

  }

//	assert(0);
  return false;

}
//@ End of file
//@ End of file
