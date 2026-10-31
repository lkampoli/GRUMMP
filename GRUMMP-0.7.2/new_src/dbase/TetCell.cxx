#include <math.h>
#include "GR_Geometry.h"
#include "GR_Cell.h"
#include "GR_CellCV.h"
#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"
#include "CubitVector.hpp"

const Vert*
TetCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert2(i >= 0 && i <= 3, "Vertex index out of range");
  const Face* const pFBase = m_faces[0];
  if (i == 3)
    return (getOppositeVert(pFBase));
  if (i == 0)
    return (pFBase->getVert(0));
  // For verts 1 and 2, return value depends on face orientation.
  if (pFBase->getRightCell() == static_cast<const Cell*>(this))
    return pFBase->getVert(i);
  else
    return pFBase->getVert(3 - i);
}

void
TetCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  const Face* const pFBase = m_faces[0];
  pFBase->getAllVertHandles(aHandles);
  // For verts 1 and 2, return value depends on face orientation, so
  // switch if needed.
  if (pFBase->getRightCell() != static_cast<const Cell*>(this)) {
    GRUMMP_Entity* pvTmp = aHandles[1];
    aHandles[1] = aHandles[2];
    aHandles[2] = pvTmp;
  }
  aHandles[3] = getOppositeVert(pFBase);
}

double
TetCell::calcSize() const
{
  assert(doFullCheck());
  const Vert *pV0, *pV1, *pV2, *pV3;
  const Face* const pFBase = m_faces[0];
  pV0 = pFBase->getVert(0);
  if (pFBase->getRightCell() == static_cast<const Cell*>(this)) {
    pV1 = pFBase->getVert(1);
    pV2 = pFBase->getVert(2);
  }
  else {
    pV1 = pFBase->getVert(2);
    pV2 = pFBase->getVert(1);
  }
  pV3 = getOppositeVert(pFBase);

  double adVec1[] =
    { pV1->x() - pV0->x(), pV1->y() - pV0->y(), pV1->z() - pV0->z() };
  double adVec2[] =
    { pV2->x() - pV0->x(), pV2->y() - pV0->y(), pV2->z() - pV0->z() };
  double adVec3[] =
    { pV3->x() - pV0->x(), pV3->y() - pV0->y(), pV3->z() - pV0->z() };
  double adCross[3];
  vCROSS3D(adVec2, adVec3, adCross);
  return dDOT3D(adCross, adVec1) / 6.;
}

bool
TetCell::isClosed() const
{
  // Check to be sure that the tet has only four vertices and that each
  // appears in exactly three faces.
  const Vert *apV[] =
    { pVInvalidVert, pVInvalidVert,
    pVInvalidVert, pVInvalidVert };
  // If the cell is deleted, it doesn't matter whether it's closed or
  // not.
  if (isDeleted())
    return true;
  int iFaceCount[4];
  for (int iF = 0; iF < 4; iF++) {
    const Face *pF = getFace(iF);
    assert(pF->isValid() && !pF->isDeleted());
    for (int iV = 0; iV < 3; iV++) {
      const Vert *pV = pF->getVert(iV);
      assert(pV->isValid() && !pV->isDeleted());
      int ii;
      for (ii = 0; ii < 4; ii++) {
	if (apV[ii] == pVInvalidVert) {
	  apV[ii] = pV;
	  iFaceCount[ii] = 1;
	  break;
	}
	else if (apV[ii] == pV) {
	  iFaceCount[ii]++;
	  break;
	}
      }
      if (ii == 4)
	return (false);
    } // Done with this vertex
  } // Done with this face
  if (iFaceCount[0] != 3 || iFaceCount[1] != 3 || iFaceCount[2] != 3
      || iFaceCount[3] != 3)
    return (false);

  // Is each vert shared by exactly three faces?  (Should be, based on
  // the previous result...
  static const int tetTuples[][3] =
    {
      { 0, 1, 2 },
      { 0, 1, 3 },
      { 0, 2, 3 },
      { 1, 2, 3 } };
  for (int ii = 0; ii < 4; ii++) {
    const Face* pFA = getFace(tetTuples[ii][0]);
    const Face* pFB = getFace(tetTuples[ii][1]);
    const Face* pFC = getFace(tetTuples[ii][2]);
    Vert* pV = findCommonVert(pFA, pFB, pFC);
    if (!pV->isValid())
      return false;
  }

  return (true);
}

void
TetCell::calcAllDihed(double adDihed[], int* const piNDihed,
		      const bool in_degrees) const
{
  *piNDihed = 6;
  const double * const adA = getVert(0)->getCoords();
  const double * const adB = getVert(1)->getCoords();
  const double * const adC = getVert(2)->getCoords();
  const double * const adD = getVert(3)->getCoords();

  double adNorm0[3], adNorm1[3], adNorm2[3], adNorm3[3];
  calcUnitNormal(adC, adB, adD, adNorm0);
  calcUnitNormal(adA, adC, adD, adNorm1);
  calcUnitNormal(adB, adA, adD, adNorm2);
  calcUnitNormal(adA, adB, adC, adNorm3);

  // The order in which the dihedrals are assigned matters for
  // computation of solid angles.  The way they're currently set up,
  // combining them as (0,1,2), (0,3,4), (1,3,5), (2,4,5) gives (in
  // order) solid angles at vertices A, B, C, D

  double dScale = in_degrees ? 180. / M_PI : 1;
  double dDot = -dDOT3D(adNorm0, adNorm1);
  adDihed[5] = GR_acos(dDot) * dScale; // Edge CD

  dDot = -dDOT3D(adNorm0, adNorm2);
  adDihed[4] = GR_acos(dDot) * dScale; // Edge BD

  dDot = -dDOT3D(adNorm0, adNorm3);
  adDihed[3] = GR_acos(dDot) * dScale; // Edge BC

  dDot = -dDOT3D(adNorm1, adNorm2);
  adDihed[2] = GR_acos(dDot) * dScale; // Edge AD

  dDot = -dDOT3D(adNorm1, adNorm3);
  adDihed[1] = GR_acos(dDot) * dScale; // Edge AC

  dDot = -dDOT3D(adNorm2, adNorm3);
  adDihed[0] = GR_acos(dDot) * dScale; // Edge AB
}

void
TetCell::calcAllSolid(double adSolid[], int* const piNSolid,
		      const bool in_degrees) const
{
  *piNSolid = 4;
  double adDihed[6];
  int iNDihed;
  calcAllDihed(adDihed, &iNDihed, in_degrees);

  double to_subtract = in_degrees ? 180 : M_PI;
  adSolid[0] = adDihed[0] + adDihed[1] + adDihed[2] - to_subtract;
  adSolid[1] = adDihed[0] + adDihed[3] + adDihed[4] - to_subtract;
  adSolid[2] = adDihed[1] + adDihed[3] + adDihed[5] - to_subtract;
  adSolid[3] = adDihed[2] + adDihed[4] + adDihed[5] - to_subtract;
}

void
TetCell::calcCircumcenter(double adCircCent[]) const
{

  //   assert(qValid());
  //   const double* const A = pVVert(0)->adCoords();
  //   const double* const B = pVVert(1)->adCoords();
  //   const double* const C = pVVert(2)->adCoords();
  //   const double* const D = pVVert(3)->adCoords();

  //   double BA[] = adDIFF3D(B, A), CA[] = adDIFF3D(C, A), DA[] = adDIFF3D(D, A);
  //   double BA_mag = dMAG3D(BA), CA_mag = dMAG3D(CA), DA_mag = dMAG3D(DA);

  //   double CD_cross[] = { CA[1] * DA[2] - DA[1] * CA[2],
  // 			CA[2] * DA[0] - DA[2] * CA[0],
  // 			CA[0] * DA[1] - DA[0] * CA[1] };
  //   double DB_cross[] = { DA[1] * BA[2] - BA[1] * DA[2],
  // 			DA[2] * BA[0] - BA[2] * DA[0],
  // 			DA[0] * BA[1] - BA[0] * DA[1] };
  //   double BC_cross[] = { BA[1] * CA[2] - CA[1] * BA[2],
  // 			BA[2] * CA[0] - CA[2] * BA[0],
  // 			BA[0] * CA[1] - CA[0] * BA[1] };

  //   double denom = 0.5 / orient3d_shew(B, C, D, A);

  //   adCircCent[0] = (BA_mag * CD_cross[0] + CA_mag * DB_cross[0] + DA_mag * BC_cross[0]) * denom;
  //   adCircCent[1] = (BA_mag * CD_cross[1] + CA_mag * DB_cross[1] + DA_mag * BC_cross[1]) * denom;
  //   adCircCent[2] = (BA_mag * CD_cross[2] + CA_mag * DB_cross[2] + DA_mag * BC_cross[2]) * denom;

  assert(isValid());
  double adB[3] =
    { getVert(1)->x() - getVert(0)->x(), getVert(1)->y() - getVert(0)->y(),
	getVert(1)->z() - getVert(0)->z() };
  double adC[3] =
    { getVert(2)->x() - getVert(0)->x(), getVert(2)->y() - getVert(0)->y(),
	getVert(2)->z() - getVert(0)->z() };
  double adD[3] =
    { getVert(3)->x() - getVert(0)->x(), getVert(3)->y() - getVert(0)->y(),
	getVert(3)->z() - getVert(0)->z() };

  double *adRow1 = adB;
  double *adRow2 = adC;
  double *adRow3 = adD;

  double dRHS1 = 0.5 * dDOT3D(adB, adB);
  double dRHS2 = 0.5 * dDOT3D(adC, adC);
  double dRHS3 = 0.5 * dDOT3D(adD, adD);

  solve3By3(adRow1, adRow2, adRow3, dRHS1, dRHS2, dRHS3, adCircCent);
#ifndef NDEBUG
  double adA[] =
    { 0, 0, 0 };
  double dComp = dDIST3D(adCircCent, adA);
  double dTest = dDIST3D(adCircCent, adB);
  assert(0 == iFuzzyComp(dComp, dTest));
  dTest = dDIST3D(adCircCent, adC);
  assert(0 == iFuzzyComp(dComp, dTest));
  dTest = dDIST3D(adCircCent, adD);
  assert(0 == iFuzzyComp(dComp, dTest));
#endif
  adCircCent[0] += getVert(0)->x();
  adCircCent[1] += getVert(0)->y();
  adCircCent[2] += getVert(0)->z();
  // The following assertion fails for flat tets.
  //assert(qCircumscribes(adCircCent));
}

double
TetCell::calcInradius() const
{

  //Inradius = 3 * Volume / sum(Triangle face Areas)

  double tot_area = 0.;
  for (int i = 0; i < 4; ++i)
    tot_area += getFace(i)->calcSize();

  return (3. * calcSize() / tot_area);

}

double
TetCell::calcCircumradius() const
{
  assert(isValid());
  if (checkOrient3D(getVert(0), getVert(1), getVert(2), getVert(3)) == 0) {
    // Tet is totally flat. Four vertices must be cocircular.
    return LARGE_DBL;
  }
  else {
    double adTemp[3];
    calcCircumcenter(adTemp);
    adTemp[0] -= getVert(0)->x();
    adTemp[1] -= getVert(0)->y();
    adTemp[2] -= getVert(0)->z();
    return (dMAG3D(adTemp));
  }
}

bool
TetCell::circumscribesPoint(const double adPoint[]) const
{
  // Returns true if adPoint lies in or on the circumsphere of *this
  return (-1
      != isInsphere(getVert(0)->getCoords(), getVert(1)->getCoords(),
		    getVert(2)->getCoords(), getVert(3)->getCoords(), adPoint));
}

void
TetCell::calcContainmentCenter(double adResult[]) const
{
  assert(isValid());
  calcCircumcenter(adResult);
  double adBary[4], REps = 1.e-10;
  int iTrips = 0;
  bool qReTry;
  do {
    qReTry = false;
    iTrips++;
    assert(iTrips <= 2);
    calcBarycentricCoords(adResult, adBary);
    int iLowBary = ((adBary[0] < REps ? 1 : 0) + (adBary[1] < REps ? 2 : 0)
	+ (adBary[2] < REps ? 4 : 0) + (adBary[3] < REps ? 8 : 0));
    if (iLowBary == 0)
      return;

    switch (iLowBary)
      {
      case 1: // Need containment center of face opposite vert 0
	{
	  if (fabs(adBary[0]) < REps)
	    break;
	  // Project onto the plane of that face.
	  const Face* pF = getOppositeFace(getVert(0));
	  pF->projectOntoFace(adResult);
	  // Re-cycle to make sure that the point is in or on the tet
	  qReTry = true;
	}
	break;
      case 2: // Need containment center of face opposite vert 1
	{
	  if (fabs(adBary[1]) < REps)
	    break;
	  // Project onto the plane of that face.
	  const Face* pF = getOppositeFace(getVert(1));
	  pF->projectOntoFace(adResult);
	  // Re-cycle to make sure that the point is in or on the tet
	  qReTry = true;
	}
	break;
      case 4: // Need containment center of face opposite vert 2
	{
	  if (fabs(adBary[2]) < REps)
	    break;
	  // Project onto the plane of that face.
	  const Face* pF = getOppositeFace(getVert(2));
	  pF->projectOntoFace(adResult);
	  // Re-cycle to make sure that the point is in or on the tet
	  qReTry = true;
	}
	break;
      case 8: // Need containment center of face opposite vert 3
	{
	  if (fabs(adBary[3]) < REps)
	    break;
	  // Project onto the plane of that face.
	  const Face* pF = getOppositeFace(getVert(3));
	  pF->projectOntoFace(adResult);
	  // Re-cycle to make sure that the point is in or on the tet
	  qReTry = true;
	}
	break;
      case 3: // Need midside of edge between verts 2 and 3
	adResult[0] = (getVert(2)->x() + getVert(3)->x()) * 0.5;
	adResult[1] = (getVert(2)->y() + getVert(3)->y()) * 0.5;
	adResult[2] = (getVert(2)->z() + getVert(3)->z()) * 0.5;
	break;
      case 5: // Need midside of edge between verts 1 and 3
	adResult[0] = (getVert(1)->x() + getVert(3)->x()) * 0.5;
	adResult[1] = (getVert(1)->y() + getVert(3)->y()) * 0.5;
	adResult[2] = (getVert(1)->z() + getVert(3)->z()) * 0.5;
	break;
      case 9: // Need midside of edge between verts 1 and 2
	adResult[0] = (getVert(1)->x() + getVert(2)->x()) * 0.5;
	adResult[1] = (getVert(1)->y() + getVert(2)->y()) * 0.5;
	adResult[2] = (getVert(1)->z() + getVert(2)->z()) * 0.5;
	break;
      case 6: // Need midside of edge between verts 0 and 3
	adResult[0] = (getVert(0)->x() + getVert(3)->x()) * 0.5;
	adResult[1] = (getVert(0)->y() + getVert(3)->y()) * 0.5;
	adResult[2] = (getVert(0)->z() + getVert(3)->z()) * 0.5;
	break;
      case 10: // Need midside of edge between verts 0 and 2
	adResult[0] = (getVert(0)->x() + getVert(2)->x()) * 0.5;
	adResult[1] = (getVert(0)->y() + getVert(2)->y()) * 0.5;
	adResult[2] = (getVert(0)->z() + getVert(2)->z()) * 0.5;
	break;
      case 12: // Need midside of edge between verts 0 and 1
	adResult[0] = (getVert(0)->x() + getVert(1)->x()) * 0.5;
	adResult[1] = (getVert(0)->y() + getVert(1)->y()) * 0.5;
	adResult[2] = (getVert(0)->z() + getVert(1)->z()) * 0.5;
	break;
      default:
	assert2(0, "Unexpected difficulty in finding containment center");
	break;
      }
  }
  while (qReTry);
}

double
TetCell::calcContainmentRadius() const
{
  double adTemp[3];
  calcContainmentCenter(adTemp);
  adTemp[0] -= getVert(0)->x();
  adTemp[1] -= getVert(0)->y();
  adTemp[2] -= getVert(0)->z();
  return (dMAG3D(adTemp));
}

void
TetCell::calcBarycentricCoords(const double adPoint[], double adBary[]) const
{
  assert(isValid());

  const Vert* pVA = getVert(0);
  const Vert* pVB = getVert(1);
  const Vert* pVC = getVert(2);
  const Vert* pVD = getVert(3);

  // The old version of this code (commented out below) has a robustness
  // problem with small tets far from the origin.

  // Solve:
  // / 1  1  1  1  \ / b0 \   / 1 \ .
  // | xa xb xc xd | | b1 |   | x |
  // | ya yb yc yd | | b2 | = | y |
  // \ za zb zc zd / \ b3 /   \ z /

  //   int i;
  //   double a2dLHS[4][4], adRHS[4], a2dMat[4][4];
  //   double dDet, dDet0, dDet1, dDet2, dDet3;

  //   a2dLHS[0][0] = 1.;
  //   a2dLHS[0][1] = pVA->dX();
  //   a2dLHS[0][2] = pVA->dY();
  //   a2dLHS[0][3] = pVA->dZ();

  //   a2dLHS[1][0] = 1.;
  //   a2dLHS[1][1] = pVB->dX();
  //   a2dLHS[1][2] = pVB->dY();
  //   a2dLHS[1][3] = pVB->dZ();

  //   a2dLHS[2][0] = 1.;
  //   a2dLHS[2][1] = pVC->dX();
  //   a2dLHS[2][2] = pVC->dY();
  //   a2dLHS[2][3] = pVC->dZ();

  //   a2dLHS[3][0] = 1.;
  //   a2dLHS[3][1] = pVD->dX();
  //   a2dLHS[3][2] = pVD->dY();
  //   a2dLHS[3][3] = pVD->dZ();

  //   adRHS[0] = 1.;
  //   adRHS[1] = adPoint[0];
  //   adRHS[2] = adPoint[1];
  //   adRHS[3] = adPoint[2];

  //   dDet = dDet4By4(a2dLHS);

  //   for (i = 0; i <= 3; ++i) {
  //     a2dMat[0][i] = adRHS    [i];
  //     a2dMat[1][i] = a2dLHS[1][i];
  //     a2dMat[2][i] = a2dLHS[2][i];
  //     a2dMat[3][i] = a2dLHS[3][i];
  //   }
  //   dDet0 = dDet4By4(a2dMat);

  //   for (i = 0; i <= 3; ++i) {
  //     a2dMat[0][i] = a2dLHS[0][i];
  //     a2dMat[1][i] = adRHS    [i];
  //     a2dMat[2][i] = a2dLHS[2][i];
  //     a2dMat[3][i] = a2dLHS[3][i];
  //   }
  //   dDet1 = dDet4By4(a2dMat);

  //   for (i = 0; i <= 3; ++i) {
  //     a2dMat[0][i] = a2dLHS[0][i];
  //     a2dMat[1][i] = a2dLHS[1][i];
  //     a2dMat[2][i] = adRHS    [i];
  //     a2dMat[3][i] = a2dLHS[3][i];
  //   }
  //   dDet2 = dDet4By4(a2dMat);

  //   for (i = 0; i <= 3; ++i) {
  //     a2dMat[0][i] = a2dLHS[0][i];
  //     a2dMat[1][i] = a2dLHS[1][i];
  //     a2dMat[2][i] = a2dLHS[2][i];
  //     a2dMat[3][i] = adRHS    [i];
  //   }
  //   dDet3 = dDet4By4(a2dMat);

  //   assert(fabs(dDet) > 1.e-12);
  //   adBary[0] = dDet0 / dDet;
  //   adBary[1] = dDet1 / dDet;
  //   adBary[2] = dDet2 / dDet;
  //   adBary[3] = dDet3 / dDet;

  // New version:

  // Solve:
  // / 1    1     1     1   \ / b0 \   /   1  \ .
  // | 0  xb-xa xc-xa xd-xa | | b1 |   | x-xa |
  // | 0  yb-ya yc-ya yd-ya | | b2 | = | y-ya |
  // \ 0  zb-za zc-za zd-za / \ b3 /   \ z-za /
  //
  // by solving the 3x3 for b1, b2, b3 and subtracting to get b0.
  // The solution for b1, b2, b3 has been immunized against some forms
  // of insignificance by this linear combination.  Each row can be
  // normalized to make it so that small determinants for the 3x3 are
  // genuinely bad news.

  double a2dLHS[3][3];

  a2dLHS[0][0] = pVB->x() - pVA->x();
  a2dLHS[0][1] = pVB->y() - pVA->y();
  a2dLHS[0][2] = pVB->z() - pVA->z();

  a2dLHS[1][0] = pVC->x() - pVA->x();
  a2dLHS[1][1] = pVC->y() - pVA->y();
  a2dLHS[1][2] = pVC->z() - pVA->z();

  a2dLHS[2][0] = pVD->x() - pVA->x();
  a2dLHS[2][1] = pVD->y() - pVA->y();
  a2dLHS[2][2] = pVD->z() - pVA->z();

  double adRHS[] =
    { adPoint[0] - pVA->x(), adPoint[1] - pVA->y(), adPoint[2] - pVA->z() };

  double adNorm[] =
    { a2dLHS[0][0] * a2dLHS[0][0] + a2dLHS[1][0] * a2dLHS[1][0]
	+ a2dLHS[2][0] * a2dLHS[2][0], a2dLHS[0][1] * a2dLHS[0][1]
	+ a2dLHS[1][1] * a2dLHS[1][1] + a2dLHS[2][1] * a2dLHS[2][1],
	a2dLHS[0][2] * a2dLHS[0][2] + a2dLHS[1][2] * a2dLHS[1][2]
	    + a2dLHS[2][2] * a2dLHS[2][2] };
  adNorm[0] = 1. / sqrt(adNorm[0]);
  adNorm[1] = 1. / sqrt(adNorm[1]);
  adNorm[2] = 1. / sqrt(adNorm[2]);

  a2dLHS[0][0] *= adNorm[0];
  a2dLHS[1][0] *= adNorm[0];
  a2dLHS[2][0] *= adNorm[0];
  adRHS[0] *= adNorm[0];

  a2dLHS[0][1] *= adNorm[1];
  a2dLHS[1][1] *= adNorm[1];
  a2dLHS[2][1] *= adNorm[1];
  adRHS[1] *= adNorm[1];

  a2dLHS[0][2] *= adNorm[2];
  a2dLHS[1][2] *= adNorm[2];
  a2dLHS[2][2] *= adNorm[2];
  adRHS[2] *= adNorm[2];

  double dDet = a2dLHS[0][0]
      * (a2dLHS[1][1] * a2dLHS[2][2] - a2dLHS[1][2] * a2dLHS[2][1])
      + a2dLHS[1][0]
	  * (a2dLHS[2][1] * a2dLHS[0][2] - a2dLHS[0][1] * a2dLHS[2][2])
      + a2dLHS[2][0]
	  * (a2dLHS[0][1] * a2dLHS[1][2] - a2dLHS[0][2] * a2dLHS[1][1]);

  assert(fabs(dDet) > 1.e-12);

  double dDet1 = adRHS[0]
      * (a2dLHS[1][1] * a2dLHS[2][2] - a2dLHS[1][2] * a2dLHS[2][1])
      + a2dLHS[1][0] * (a2dLHS[2][1] * adRHS[2] - adRHS[1] * a2dLHS[2][2])
      + a2dLHS[2][0] * (adRHS[1] * a2dLHS[1][2] - adRHS[2] * a2dLHS[1][1]);

  double dDet2 = a2dLHS[0][0]
      * (adRHS[1] * a2dLHS[2][2] - adRHS[2] * a2dLHS[2][1])
      + adRHS[0] * (a2dLHS[2][1] * a2dLHS[0][2] - a2dLHS[0][1] * a2dLHS[2][2])
      + a2dLHS[2][0] * (a2dLHS[0][1] * adRHS[2] - a2dLHS[0][2] * adRHS[1]);

  double dDet3 = a2dLHS[0][0]
      * (a2dLHS[1][1] * adRHS[2] - a2dLHS[1][2] * adRHS[1])
      + a2dLHS[1][0] * (adRHS[1] * a2dLHS[0][2] - a2dLHS[0][1] * adRHS[2])
      + adRHS[0] * (a2dLHS[0][1] * a2dLHS[1][2] - a2dLHS[0][2] * a2dLHS[1][1]);

  adBary[3] = dDet3 / dDet;
  adBary[2] = dDet2 / dDet;
  adBary[1] = dDet1 / dDet;
  adBary[0] = 1. - adBary[1] - adBary[2] - adBary[3];
}

double
TetCell::calcShortestEdgeLength() const
{
  register double dDX, dDY, dDZ;
  register const Vert *pV0 = getVert(0);
  register const Vert *pV1 = getVert(1);
  register const Vert *pV2 = getVert(2);
  register const Vert *pV3 = getVert(3);

  dDX = pV0->x() - pV1->x();
  dDY = pV0->y() - pV1->y();
  dDZ = pV0->z() - pV1->z();
  double dLen01_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  dDX = pV0->x() - pV2->x();
  dDY = pV0->y() - pV2->y();
  dDZ = pV0->z() - pV2->z();
  double dLen02_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  dDX = pV0->x() - pV3->x();
  dDY = pV0->y() - pV3->y();
  dDZ = pV0->z() - pV3->z();
  double dLen03_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  dDX = pV1->x() - pV2->x();
  dDY = pV1->y() - pV2->y();
  dDZ = pV1->z() - pV2->z();
  double dLen12_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  dDX = pV1->x() - pV3->x();
  dDY = pV1->y() - pV3->y();
  dDZ = pV1->z() - pV3->z();
  double dLen13_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  dDX = pV2->x() - pV3->x();
  dDY = pV2->y() - pV3->y();
  dDZ = pV2->z() - pV3->z();
  double dLen23_Sq = dDX * dDX + dDY * dDY + dDZ * dDZ;

  double dShortest_Sq = min(
      min(min(dLen01_Sq, dLen02_Sq), min(dLen03_Sq, dLen12_Sq)),
      min(dLen13_Sq, dLen23_Sq));
  return (sqrt(dShortest_Sq));
}
bool
TetCell::isPointOnEdge(const double adPoint[], GR_index_t & iV0,
		       GR_index_t & iV1)
{
  double adBary[4];
  calcBarycentricCoords(adPoint, adBary);
  int i, iLowBary = 0;
  for (i = 0; i <= 3; i++)
    if (adBary[i] < 1e-10 && adBary[i] > -1e-10) {
      adBary[i] = 0;
      iLowBary++;
    }

  if (iLowBary == 2) {

    if (adBary[0] < 1e-10) {
      if (adBary[1] < 1e-10) {
	iV0 = 2;
	iV1 = 3;
      }
      else if (adBary[2] < 1e-10) {
	iV0 = 1;
	iV1 = 3;
      }
      else {
	iV0 = 2;
	iV1 = 1;
      }
    }
    else if (adBary[1] < 1e-10) {
      if (adBary[2] < 1e-10) {
	iV0 = 0;
	iV1 = 3;
      }
      else {
	iV0 = 2;
	iV1 = 0;
      }
    }
    else {
      iV0 = 0;
      iV1 = 1;
    }
    assert(getVert(iV0)->isValid());
    assert(getVert(iV1)->isValid());
    return true;
  }
  return false;
}
void
TetCellCV::vCircumcenter(double adRes[]) const
{
  assert(qValid());
  const double * const adA = getVert(0)->getCoords();
  const double * const adB = getVert(1)->getCoords();
  const double * const adC = getVert(2)->getCoords();
  const double * const adD = getVert(3)->getCoords();

  double adRow1[] =
    { adA[0] - adB[0], adA[1] - adB[1], adA[2] - adB[2] };
  double adRow2[] =
    { adA[0] - adC[0], adA[1] - adC[1], adA[2] - adC[2] };
  double adRow3[] =
    { adA[0] - adD[0], adA[1] - adD[1], adA[2] - adD[2] };

  double dRHS1 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adB, adB));
  double dRHS2 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adC, adC));
  double dRHS3 = 0.5 * (dDOT3D(adA, adA) - dDOT3D(adD, adD));

  solve3By3(adRow1, adRow2, adRow3, dRHS1, dRHS2, dRHS3, adRes);

#ifndef NDEBUG
  double dComp = dDIST3D(adRes, adA);
  double dTest = dDIST3D(adRes, adB);
  assert(0 == iFuzzyComp(dComp, dTest));
  dTest = dDIST3D(adRes, adC);
  assert(0 == iFuzzyComp(dComp, dTest));
  dTest = dDIST3D(adRes, adD);
  assert(0 == iFuzzyComp(dComp, dTest));
#endif
// The following assertion fails for flat tets.
// assert(qCircumscribes(adRes));

}

void
TetCellCV::calcAllDihed(double adDihed[], int* const piNDihed,
			const bool in_degrees) const
{
  *piNDihed = 6;
  const double * const adA = getVert(0)->getCoords();
  const double * const adB = getVert(1)->getCoords();
  const double * const adC = getVert(2)->getCoords();
  const double * const adD = getVert(3)->getCoords();

  double adNorm0[3], adNorm1[3], adNorm2[3], adNorm3[3];
  calcUnitNormal(adC, adB, adD, adNorm0);
  calcUnitNormal(adA, adC, adD, adNorm1);
  calcUnitNormal(adB, adA, adD, adNorm2);
  calcUnitNormal(adA, adB, adC, adNorm3);

  // The order in which the dihedrals are assigned matters for
  // computation of solid angles.  The way they're currently set up,
  // combining them as (0,1,2), (0,3,4), (1,3,5), (2,4,5) gives (in
  // order) solid angles at vertices A, B, C, D

  double dScale = in_degrees ? 180. / M_PI : 1;
  double dDot = -dDOT3D(adNorm0, adNorm1);
  adDihed[5] = GR_acos(dDot) * dScale; // Edge CD

  dDot = -dDOT3D(adNorm0, adNorm2);
  adDihed[4] = GR_acos(dDot) * dScale; // Edge BD

  dDot = -dDOT3D(adNorm0, adNorm3);
  adDihed[3] = GR_acos(dDot) * dScale; // Edge BC

  dDot = -dDOT3D(adNorm1, adNorm2);
  adDihed[2] = GR_acos(dDot) * dScale; // Edge AD

  dDot = -dDOT3D(adNorm1, adNorm3);
  adDihed[1] = GR_acos(dDot) * dScale; // Edge AC

  dDot = -dDOT3D(adNorm2, adNorm3);
  adDihed[0] = GR_acos(dDot) * dScale; // Edge AB
}

void
TetCellCV::calcAllSolid(double adSolid[], int* const piNSolid,
			const bool in_degrees) const
{
  *piNSolid = 4;
  double adDihed[6];
  int iNDihed;
  calcAllDihed(adDihed, &iNDihed, in_degrees);

  double to_subtract = in_degrees ? 180 : M_PI;
  adSolid[0] = adDihed[0] + adDihed[1] + adDihed[2] - to_subtract;
  adSolid[1] = adDihed[0] + adDihed[3] + adDihed[4] - to_subtract;
  adSolid[2] = adDihed[1] + adDihed[3] + adDihed[5] - to_subtract;
  adSolid[3] = adDihed[2] + adDihed[4] + adDihed[5] - to_subtract;
}

double
TetCellCV::calcSize() const
{
  assert(iFullCheck());
  const Vert *pV0 = getVert(0);
  const Vert *pV1 = getVert(1);
  const Vert *pV2 = getVert(2);
  const Vert *pV3 = getVert(3);

  double adVec1[] =
    { pV1->x() - pV0->x(), pV1->y() - pV0->y(), pV1->z() - pV0->z() };
  double adVec2[] =
    { pV2->x() - pV0->x(), pV2->y() - pV0->y(), pV2->z() - pV0->z() };
  double adVec3[] =
    { pV3->x() - pV0->x(), pV3->y() - pV0->y(), pV3->z() - pV0->z() };
  double adCross[3];
  vCROSS3D(adVec2, adVec3, adCross);
  return fabs(dDOT3D(adCross, adVec1)) / 6.;
}
