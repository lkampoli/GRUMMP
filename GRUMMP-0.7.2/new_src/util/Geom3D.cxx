#include <math.h>
#include "GR_misc.h"
#include "GR_Geometry.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_AdaptPred.h"
#include "GR_Vec.h"

#define MINMAX4(da, db, dc, dd, dMin, dMax) { \
		double dMax1, dMax2, dMin1, dMin2; \
		if (da > db) {dMax1 = da; dMin1 = db;} \
		else         {dMax1 = db; dMin1 = da;} \
		if (dc > dd) {dMax2 = dc; dMin2 = dd;} \
		else         {dMax2 = dd; dMin2 = dc;} \
		if (dMax1 > dMax2) dMax = dMax1; \
		else               dMax = dMax2; \
		if (dMin1 < dMin2) dMin = dMin1; \
		else               dMin = dMin2; \
}

void
calcNormal3D(const double adA[3], const double adB[3], const double adC[3],
	     double * const pdRes)
{
  assert(pdRes != NULL);
  double adDiffBA[3] =
  adDIFF3D(adB, adA);
  double adDiffCA[3] =
  adDIFF3D(adC, adA);
  vCROSS3D(adDiffBA, adDiffCA, pdRes);

}

void
calcQuadNormal(const double pt0[3], const double pt1[3], const double pt2[3],
	       const double pt3[3], double norm[3])
{
  // Need midsides
  double mid01[] =
	{ 0.5 * (pt0[0] + pt1[0]), 0.5 * (pt0[1] + pt1[1]), 0.5
	    * (pt0[2] + pt1[2]) };
  double mid12[] =
	{ 0.5 * (pt1[0] + pt2[0]), 0.5 * (pt1[1] + pt2[1]), 0.5
	    * (pt1[2] + pt2[2]) };
  double mid23[] =
	{ 0.5 * (pt2[0] + pt3[0]), 0.5 * (pt2[1] + pt3[1]), 0.5
	    * (pt2[2] + pt3[2]) };
  double mid30[] =
	{ 0.5 * (pt3[0] + pt0[0]), 0.5 * (pt3[1] + pt0[1]), 0.5
	    * (pt3[2] + pt0[2]) };

  // Now the differences across the center
  double diffRL[] = adDIFF3D(mid01, mid23);
  double diffTB[] = adDIFF3D(mid12, mid30);
  vCROSS3D(diffRL, diffTB, norm);
}

//@ 3D orientation test
int
checkOrient3D(const Vert* const pVertA, const Vert* const pVertB,
	      const Vert* const pVertC, const Vert* const pVertD)
{
  assert(pVertA->isValid());
  assert(pVertB->isValid());
  assert(pVertC->isValid());
  assert(pVertD->isValid());

  return (checkOrient3D(pVertA->getCoords(), pVertB->getCoords(),
			pVertC->getCoords(), pVertD->getCoords()));
}

int
checkOrient3D(const Face* const pF, const double adA[3])
{
  assert(pF->isValid());
  return (checkOrient3D(pF->getVert(0)->getCoords(),
			pF->getVert(1)->getCoords(),
			pF->getVert(2)->getCoords(), adA));
}

int
checkOrient3D(const Face* const pF, const Vert * const pV)
{
  assert(pF->isValid());
  return (checkOrient3D(pF->getVert(0)->getCoords(),
			pF->getVert(1)->getCoords(),
			pF->getVert(2)->getCoords(), pV->getCoords()));
}

//@@ 3D orientation predicate, with adaptive precision arithmetic
int
checkOrient3D(const double adA[3], const double adB[3], const double adC[3],
	      const double adD[3])
// Computes the orientation of four verts in 3D, using as nearly
// exact arithmetic as required.
{
  if (qAdaptPred) {
    double dResult = orient3d_shew(adA, adB, adC, adD);
    if (dResult > 0)
      return (-1);
    else if (dResult < 0)
      return (1);
    else
      return 0;
  }
  else {
    double dXa = adA[0];
    double dYa = adA[1];
    double dZa = adA[2];

    double dXb = adB[0];
    double dYb = adB[1];
    double dZb = adB[2];

    double dXc = adC[0];
    double dYc = adC[1];
    double dZc = adC[2];

    double dXd = adD[0];
    double dYd = adD[1];
    double dZd = adD[2];

    double dMaxX, dMinX, dMaxY, dMinY, dMaxZ, dMinZ;

    double dDX2 = dXb - dXa;
    double dDX3 = dXc - dXa;
    double dDX4 = dXd - dXa;
    MINMAX4(dXa, dXb, dXc, dXd, dMinX, dMaxX);

    double dDY2 = dYb - dYa;
    double dDY3 = dYc - dYa;
    double dDY4 = dYd - dYa;
    MINMAX4(dYa, dYb, dYc, dYd, dMinY, dMaxY);

    double dDZ2 = dZb - dZa;
    double dDZ3 = dZc - dZa;
    double dDZ4 = dZd - dZa;
    MINMAX4(dZa, dZb, dZc, dZd, dMinZ, dMaxZ);
    double dMax = max(max(dMaxX - dMinX, dMaxY - dMinY), dMaxZ - dMinZ);

    // dDet is proportional to the cell volume
    double dDet = (dDX2 * (dDY3 * dDZ4 - dDY4 * dDZ3)
	+ dDX3 * (dDY4 * dDZ2 - dDY2 * dDZ4)
	+ dDX4 * (dDY2 * dDZ3 - dDY3 * dDZ2));

    //   // Compute a length scale based on edge lengths.
    //   double dScale = ( dDIST3D(adA, adB) + dDIST3D(adA, adC) +
    // 		    dDIST3D(adA, adD) + dDIST3D(adB, adC) +
    // 		    dDIST3D(adB, adD) + dDIST3D(adC, adD) ) / 6.;

    //   dDet /= (dScale*dScale*dScale);

    //   double dError = 1.e-13;

    // Compute an upper bound on the error bound.

    const double dMachEps = 2.22044605e-13; // about 2^-52 * 1000;

    double dError = dMachEps * dMax * dMax * dMax;

    if (dDet > dError)
      return (1);
    else if (dDet < -dError)
      return (-1);

    //     If neither of those two worked, compute a more accurate error bound.
    //     The stuff in parentheses is the result of perturbing each term in
    //     the determinant by tweaking each multiplicand by the same amount.
    //     That amount is a multiple of machine zero.
    //       dError = dMachEps *
    //         (fabs(dDX2*dDY3) + fabs(dDX2*dDZ4) + fabs(dDY3*dDZ4) +
    //          fabs(dDX3*dDY4) + fabs(dDX3*dDZ2) + fabs(dDY4*dDZ2) +
    //          fabs(dDX4*dDY2) + fabs(dDX4*dDZ3) + fabs(dDY2*dDZ3) -
    //          fabs(dDZ2*dDY3) + fabs(dDZ2*dDX4) + fabs(dDY3*dDX4) +
    //          fabs(dDZ3*dDY4) + fabs(dDZ3*dDX2) + fabs(dDY4*dDX2) +
    //          fabs(dDZ4*dDY2) + fabs(dDZ4*dDX3) + fabs(dDY2*dDX3));
    //       if (dDet > dError)
    //         return (1);
    //       else if (dDet < -dError)
    //         return (-1);
    //       else
    return (0);
  }
}

//@@ Insphere primitive with adaptive arithmetic
int
isInsphere(const Vert* const pVertA, const Vert* const pVertB,
	   const Vert* const pVertC, const Vert* const pVertD,
	   const Vert* const pVertE)
// Determines whether pVertE is inside or outside the circumsphere
// of the tet formed by pVertA-D, using as nearly exact arithmetic
// as required.  Returns 1 if pVertE is inside the circumsphere of
// the tet, -1 if outside, 0 if on.
{
  assert(pVertA->isValid());
  assert(pVertB->isValid());
  assert(pVertC->isValid());
  assert(pVertD->isValid());
  assert(pVertE->isValid());

  return isInsphere(pVertA->getCoords(), pVertB->getCoords(),
		    pVertC->getCoords(), pVertD->getCoords(),
		    pVertE->getCoords());
}

int
isInsphere(const double adA[3], const double adB[3], const double adC[3],
	   const double adD[3], const double adE[3])
{
  if (qAdaptPred) {
    double dResult = insphere_shew(adA, adB, adC, adD, adE)
	* orient3d_shew(adA, adB, adC, adD);
    if (dResult > 0)
      return 1;
    else if (dResult < 0)
      return -1;
    else
      return 0;
  }
  else {
    double dXa = adA[0];
    double dYa = adA[1];
    double dZa = adA[2];

    double dXb = adB[0];
    double dYb = adB[1];
    double dZb = adB[2];

    double dXc = adC[0];
    double dYc = adC[1];
    double dZc = adC[2];

    double dXd = adD[0];
    double dYd = adD[1];
    double dZd = adD[2];

    double dXe = adE[0];
    double dYe = adE[1];
    double dZe = adE[2];

    double a2dInSphMat[4][4], dWa;
    dWa = dXa * dXa + dYa * dYa + dZa * dZa;

    a2dInSphMat[0][0] = dXb - dXa;
    a2dInSphMat[0][1] = dYb - dYa;
    a2dInSphMat[0][2] = dZb - dZa;
    a2dInSphMat[0][3] = dXb * dXb + dYb * dYb + dZb * dZb - dWa;

    a2dInSphMat[1][0] = dXc - dXa;
    a2dInSphMat[1][1] = dYc - dYa;
    a2dInSphMat[1][2] = dZc - dZa;
    a2dInSphMat[1][3] = dXc * dXc + dYc * dYc + dZc * dZc - dWa;

    a2dInSphMat[2][0] = dXd - dXa;
    a2dInSphMat[2][1] = dYd - dYa;
    a2dInSphMat[2][2] = dZd - dZa;
    a2dInSphMat[2][3] = dXd * dXd + dYd * dYd + dZd * dZd - dWa;

    a2dInSphMat[3][0] = dXe - dXa;
    a2dInSphMat[3][1] = dYe - dYa;
    a2dInSphMat[3][2] = dZe - dZa;
    a2dInSphMat[3][3] = dXe * dXe + dYe * dYe + dZe * dZe - dWa;

    // Set up a scale that is the average of the distance from A to each
    // of the other verts.  Use some trickery to take advantage of
    // already knowing what the differences in coordinates are.
    double dAveLen = 0.25
	* (dMAG3D(a2dInSphMat[0]) + dMAG3D(a2dInSphMat[1])
	    + dMAG3D(a2dInSphMat[2]) + dMAG3D(a2dInSphMat[3]));
    double dScale = pow(dAveLen, 5);
    double dDet = calcDet4By4(a2dInSphMat) / dScale
	* checkOrient3D(adA, adB, adC, adD);
    double dEps = 1.e-12;
    if (dDet > dEps)
      return (-1);
    else if (dDet < -dEps)
      return (1);
    else
      return (0);
  }
}

//@ Determine whether a vertex lies inside, on or outside a tet
int
isInTet(const double testPt[3], const double ptA[3], const double ptB[3],
	const double ptC[3], const double ptD[3])
// For pV to be inside the tet formed by pVertA-D, the orientation
// has to be right wrt every face.  The return value is -1 if the
// vert lies outside the tet, 0 if it's coplanar with all four
// faces (impossible for valid tets), 1 if it coincides with a vert
// in the tet, 2 if it lies on an edge, 3 if it lies on a face, and
// a 4 if it lies strictly inside.
{
  int iRetVal = 0;
  int iCheck = checkOrient3D(ptB, ptD, ptC, testPt);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient3D(ptC, ptD, ptA, testPt);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient3D(ptA, ptD, ptB, testPt);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  iCheck = checkOrient3D(ptA, ptB, ptC, testPt);
  if (iCheck == -1)
    return -1;
  iRetVal += iCheck;

  return (iRetVal);
}

int
isInTet(const double testPt[3], const Cell* const cell)
// For pV to be inside the tet formed by pVertA-D, the orientation
// has to be right wrt every face.  The return value is -1 if the
// vert lies outside the tet, 0 if it's coplanar with all four
// faces (impossible for valid tets), 1 if it coincides with a vert
// in the tet, 2 if it lies on an edge, 3 if it lies on a face, and
// a 4 if it lies strictly inside.
{
  assert(cell->doFullCheck() && cell->getType() == CellSkel::eTet);
  return isInTet(testPt, cell->getVert(0)->getCoords(),
		 cell->getVert(1)->getCoords(), cell->getVert(2)->getCoords(),
		 cell->getVert(3)->getCoords());
}

////@ Determine whether a vertex lies inside, on or outside a tet
//int isCoplanarInTet(const Vert* const pV, const Vert* const pVertA,
//		const Vert* const pVertB, const Vert* const pVertC,
//		const Vert* const pVertD)
//// For pV to be inside the tet formed by pVertA-D, the orientation
//// has to be right wrt every face.  The return value is -1 if the
//// vert lies outside the tet, 0 if it's coplanar with all four
//// faces (impossible for valid tets), 1 if it coincides with a vert
//// in the tet, 2 if it lies on an edge, 3 if it lies on a face, and
//// a 4 if it lies strictly inside.
//{
//	assert(pV->isValid() && pV->getSpaceDimen() == 3);
//	assert(pVertA->isValid() && pVertA->getSpaceDimen() == 3);
//	assert(pVertB->isValid() && pVertB->getSpaceDimen() == 3);
//	assert(pVertC->isValid() && pVertC->getSpaceDimen() == 3);
//	assert(pVertD->isValid() && pVertD->getSpaceDimen() == 3);
//
//	// Find a normal and two vectors in the plane.
//	double adNorm[3];
//	calcNormal(pVertA->getCoords(), pVertB->getCoords(), pVertC->getCoords(), adNorm);
//	NORMALIZE3D(adNorm);
//
//	double adBasisX[] =
//	{pVertA->x() - pVertB->x(), pVertA->y() - pVertB->y(), pVertA->z()
//			- pVertB->z()};
//	NORMALIZE3D(adBasisX);
//
//	double adBasisY[3];
//	vCROSS3D(adNorm, adBasisX, adBasisY);
//
//	assert(iFuzzyComp(dMAG3D(adBasisY), 1.) == 0);
//	assert(qFuzzyPerp3D(adBasisX, adBasisY));
//	assert(qFuzzyPerp3D(adNorm, adBasisY));
//	assert(qFuzzyPerp3D(adBasisX, adNorm));
//
//	// Project onto a plane
//	Vert VProjA, VProjB, VProjC, VProjD, VProj;
//	double adTemp[2];
//
//	adTemp[0] = dDOT3D(pVertA->getCoords(), adBasisX);
//	adTemp[1] = dDOT3D(pVertA->getCoords(), adBasisY);
//	VProjA.setCoords(2, adTemp);
//
//	adTemp[0] = dDOT3D(pVertB->getCoords(), adBasisX);
//	adTemp[1] = dDOT3D(pVertB->getCoords(), adBasisY);
//	VProjB.setCoords(2, adTemp);
//
//	adTemp[0] = dDOT3D(pVertC->getCoords(), adBasisX);
//	adTemp[1] = dDOT3D(pVertC->getCoords(), adBasisY);
//	VProjC.setCoords(2, adTemp);
//
//	adTemp[0] = dDOT3D(pVertD->getCoords(), adBasisX);
//	adTemp[1] = dDOT3D(pVertD->getCoords(), adBasisY);
//	VProjD.setCoords(2, adTemp);
//
//	adTemp[0] = dDOT3D(pV ->getCoords(), adBasisX);
//	adTemp[1] = dDOT3D(pV ->getCoords(), adBasisY);
//	VProj .setCoords(2, adTemp);
//
//	// Check to see if pV lies inside any of the triangles implied by
//	// verts A-D
//	int iCheck = isInTri(&VProjA, &VProjB, &VProjC, &VProj);
//	if (iCheck != -1) return iCheck;
//
//	iCheck = isInTri(&VProjA, &VProjD, &VProjB, &VProj);
//	if (iCheck != -1) return iCheck;
//
//	iCheck = isInTri(&VProjB, &VProjD, &VProjC, &VProj);
//	if (iCheck != -1) return iCheck;
//
//	iCheck = isInTri(&VProjC, &VProjD, &VProjA, &VProj);
//	if (iCheck != -1) return iCheck;
//
//	return -1;
//}

//int isReflex(const EdgeFace* const pF)
//{
//	assert(pF->getType() == Face::eEdgeFace);
//	assert(pF->getLeftCell()->getType()  == Cell::eTriCell);
//	assert(pF->getRightCell()->getType() == Cell::eTriCell);
//	TriCell* pTC0 = dynamic_cast<TriCell*>(pF->getLeftCell ());
//	TriCell* pTC1 = dynamic_cast<TriCell*>(pF->getRightCell());
//	const Vert* pVA = pF->getVert(0);
//	const Vert* pVB = pF->getVert(1);
//	const Vert* pVC = pTC0->getOppositeVert(pF);
//	const Vert* pVD = pTC1->getOppositeVert(pF);
//	return (-checkOrient3D(pVA, pVB, pVC, pVD));
//}

int
calcSegmentTriangleIntersection(const Vert * const pV0, const Vert * const pV1,
				const Vert* const pFV0, const Vert* const pFV1,
				const Vert* const pFV2, double tol,
				double * intersection)
{
  return calcSegmentTriangleIntersection(pV0->getCoords(), pV1->getCoords(),
					 pFV0, pFV1, pFV2, tol, intersection);
}
int
calcSegmentTriangleIntersection(const double l0[3], const double l1[3],
				const Vert* const pFV0, const Vert* const pFV1,
				const Vert* const pFV2, double tol,
				double * intersection)
{
  if (iFuzzyComp(dDIST3D(l0, l1), 0.) == 0)
    return -1;

  const double * p0 = pFV0->getCoords();
  const double * p1 = pFV1->getCoords();
  const double * p2 = pFV2->getCoords();

  const double u[] =
    { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] }; //Parametric plane axis 1
  const double v[] =
    { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] }; //Parametric plane axis 2
  double normal[3];  //Normal to plane (triangle)
  calcNormal3D(p0, p1, p2, normal);

  if (iFuzzyComp(dMAG3D_SQ(normal), 0.) == 0)
    return -1; //Degenerate

  double dir1[] =
    { l1[0] - l0[0], l1[1] - l0[1], l1[2] - l0[2] };
  double dir2[] =
    { l0[0] - p0[0], l0[1] - p0[1], l0[2] - p0[2] };

  double a = -dDOT3D(normal, dir2);
  double b = dDOT3D(normal, dir1);

  // need to check if line is in the plane
  if (iFuzzyComp(b, 0.) == 0) {
    if (iFuzzyComp(a, 0.) == 0) // on the same plane
      return 0;
    else
      return -1; // Line parallel to plane
  }

  double r = a / b; //Intersection parameter along line;

  if (iFuzzyComp(r, 0.) == -1 || iFuzzyComp(r, 1.) == 1)
    return -1;
  double I[] =
    { l0[0] + r * dir1[0], l0[1] + r * dir1[1], l0[2] + r * dir1[2] };
  //Intersection point between plane and line
  if (intersection != NULL) {
    intersection[0] = I[0];
    intersection[1] = I[1];
    intersection[2] = I[2];
  }
  //Is this point inside triangle?
  //Use parametric plane equations to find out.

  double w[] =
    { I[0] - p0[0], I[1] - p0[1], I[2] - p0[2] };

  double uu = dDOT3D(u, u);
  double uv = dDOT3D(u, v);
  double vv = dDOT3D(v, v);
  double wu = dDOT3D(w, u);
  double wv = dDOT3D(w, v);
  double den = uv * uv - uu * vv;

  double s = (uv * wv - vv * wu) / den;
  if (iFuzzyComp(s, -tol) == -1 || iFuzzyComp(s, 1. + tol) == 1)
    return -1;

  double t = (uv * wu - uu * wv) / den;
  if (iFuzzyComp(t, -tol) == -1 || iFuzzyComp(s + t, 1. + tol) == 1)
    return -1;
  return 1;
}
void
calcCircumcenter3D(const Vert* const pVA, const Vert* pVB,
		   const Vert* const pVC, double * adCircCent)
{
#ifndef NDEBUG
  double adA[3] =
    { 0, 0, 0 };
#endif
  double adB[3] =
    { pVB->x() - pVA->x(), pVB->y() - pVA->y(), pVB->z() - pVA->z() };
  double adC[3] =
    { pVC->x() - pVA->x(), pVC->y() - pVA->y(), pVC->z() - pVA->z() };

  double adRow1[3], adRow2[3], adRow3[3];

  adRow1[0] = adB[0];
  adRow1[1] = adB[1];
  adRow1[2] = adB[2];

  adRow2[0] = adC[0];
  adRow2[1] = adC[1];
  adRow2[2] = adC[2];

  vCROSS3D(adRow1, adRow2, adRow3);

  double dRHS1 = 0.5 * dDOT3D(adB, adB);
  double dRHS2 = 0.5 * dDOT3D(adC, adC);
  double dRHS3 = 0;

  solve3By3(adRow1, adRow2, adRow3, dRHS1, dRHS2, dRHS3, adCircCent);

#ifndef NDEBUG
  double dMagA = dDIST3D(adCircCent, adA);
  double dMagB = dDIST3D(adCircCent, adB);
  double dMagC = dDIST3D(adCircCent, adC);
  assert(0 == iFuzzyComp(dMagA, dMagB));
  assert(0 == iFuzzyComp(dMagA, dMagC));
#endif

  adCircCent[0] += pVA->x();
  adCircCent[1] += pVA->y();
  adCircCent[2] += pVA->z();
}

int
closestOnTriangle(const double pt[3], const double t0[3], const double t1[3],
		  const double t2[3], double close[3])
{

  double u[] = adDIFF3D(t1, t0), v[] = adDIFF3D(t2, t0), w[] =
  adDIFF3D(t0, pt),

  a = dDOT3D(u, u), b = dDOT3D(u, v), c = dDOT3D(v, v), d = dDOT3D(u, w), e =
      dDOT3D(v, w),

  det = fabs(a * c - b * b), s = b * e - c * d, t = b * d - a * e;

  int zone;

  if (s + t <= det) {
    if (s < 0.) {
      if (t < 0.)
	zone = 4;
      else
	zone = 3;
    }
    else if (t < 0.)
      zone = 5;
    else
      zone = 0;
  }
  else {
    if (s < 0.)
      zone = 2;
    else if (t < 0.)
      zone = 6;
    else
      zone = 1;
  }

  switch (zone)
    {

    case 0:
      { //closest point inside the triangle.
	double inv_det = 1. / det;
	s *= inv_det;
	t *= inv_det;
	break;
      }

    case 1:
      { //closest point along s + t = 1;
	//min dist at (c + e - b - d) / (a - 2b + c)
	double num = c + e - b - d;
	if (num <= 0.)
	  s = 0;
	else {
	  double den = a - 2. * b + c; //always positive.
	  s = num >= den ? 1. : num / den;
	}
	t = 1. - s;
	break;
      }

    case 2:
      { //closest point either along s + t = 1 or s = 0.
	double C0 = b + d, C1 = c + e;
	if (C1 > C0) { //closest point alogn s + t = 1. (like case 1).
	  double num = C1 - C0;
	  double den = a - 2. * b + c;  //always positive.
	  s = num >= den ? 1. : num / den;
	  t = 1. - s;
	}
	else { //closest point along s = 0. (like case 3)
	  s = 0.;
	  t = e >= 0. ? 0. : -e >= c ? 1. : -e / c;
	}
	break;
      }

    case 3:
      { //closest point along s = 0. c must be positive (=v*v)
	//min dist at t' such that ct' + e = 0
	s = 0.;
	t = e >= 0. ? 0. : -e >= c ? 1. : -e / c;
	break;
      }

    case 4:
      { //closest point along s = 0 or t = 0
	if (d < 0.) { //closest point along t = 0.
	  s = -d >= a ? 1. : -d / a;
	  t = 0.;
	}
	else {
	  s = 0.;
	  t = e >= 0. ? 0. : -e >= c ? 1. : -e / c;
	}
	break;
      }

    case 5:
      { //closest point along t = 0.
	s = d >= 0. ? 0. : -d >= a ? 1. : -d / a;
	t = 0.;
	break;
      }

    case 6:
      { //closest point along s + t = 1 or along t = 0.
	double C0 = a + d;
	double C1 = b + e;
	if (C1 < C0) { //closest point along s + t = 1
	  double num = C0 - C1;
	  double den = a - 2. * b + c;
	  t = num >= den ? 1. : num / den;
	  s = 1. - t;
	}
	else { //closest point along t = 0.
	  s = d >= 0. ? 0. : -d >= a ? 1. : -d / a;
	  t = 0.;
	}
	break;
      }

    default:
      assert(0);
      break;

    }

  close[0] = t0[0] + s * u[0] + t * v[0];
  close[1] = t0[1] + s * u[1] + t * v[1];
  close[2] = t0[2] + s * u[2] + t * v[2];
  return zone;
}

double
dVolTet(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	const Vert* pV3)
{
  double u = dDIST3D(pV0->getCoords(), pV1->getCoords());
  double U = dDIST3D(pV2->getCoords(), pV3->getCoords());
  double v = dDIST3D(pV1->getCoords(), pV2->getCoords());
  double V = dDIST3D(pV0->getCoords(), pV3->getCoords());
  double w = dDIST3D(pV2->getCoords(), pV0->getCoords());
  double W = dDIST3D(pV1->getCoords(), pV3->getCoords());

  double X = (w - U + v) * (U + v + w);
  double x = (U - v + w) * (v - w + U);
  double Y = (u - V + w) * (V + w + u);
  double y = (V - w + u) * (w - u + V);
  double Z = (v - W + u) * (W + u + v);
  double z = (W - u + v) * (u - v + W);

  double a = sqrt(x * Y * Z);
  double b = sqrt(y * Z * X);
  double c = sqrt(z * X * Y);
  double d = sqrt(x * y * z);

  double volume = sqrt(
      (-a + b + c + d) * (a - b + c + d) * (a + b - c + d) * (a + b + c - d))
      / (192 * u * v * w);

  return volume;
}

double
dCircumradius(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	      const Vert* pV3)
{
  double vol = dVolTet(pV0, pV1, pV2, pV3);
  double u = dDIST3D(pV0->getCoords(), pV1->getCoords());
  double U = dDIST3D(pV2->getCoords(), pV3->getCoords());
  double v = dDIST3D(pV1->getCoords(), pV2->getCoords());
  double V = dDIST3D(pV0->getCoords(), pV3->getCoords());
  double w = dDIST3D(pV2->getCoords(), pV0->getCoords());
  double W = dDIST3D(pV1->getCoords(), pV3->getCoords());

  double K = dAreaTriFromEdgeLengths(u * U, v * V, w * W);
  return K / (6 * vol);
}

double
dInradius(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	  const Vert* pV3)
{
  double vol = dVolTet(pV0, pV1, pV2, pV3);
  double u = dDIST3D(pV0->getCoords(), pV1->getCoords());
  double U = dDIST3D(pV2->getCoords(), pV3->getCoords());
  double v = dDIST3D(pV1->getCoords(), pV2->getCoords());
  double V = dDIST3D(pV0->getCoords(), pV3->getCoords());
  double w = dDIST3D(pV2->getCoords(), pV0->getCoords());
  double W = dDIST3D(pV1->getCoords(), pV3->getCoords());

  double area_uVW = dAreaTriFromEdgeLengths(u, V, W);
  double area_UvW = dAreaTriFromEdgeLengths(U, v, W);
  double area_UVw = dAreaTriFromEdgeLengths(U, V, w);
  double area_uvw = dAreaTriFromEdgeLengths(u, v, w);
  logMessage(MSG_DETAIL, "Areas: %12f %12f %12f %12f\n", area_uVW, area_UvW,
	     area_UVw, area_uvw);

  return 3 * vol / (area_uVW + area_UvW + area_UVw + area_uvw);
}
//@ End of file
