#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Geometry.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"
#include "CubitVector.hpp"

const Vert *
TriCell::getVert(const int i) const
{
  //assert (doFullCheck ());
  assert2(i >= 0 && i <= 2, "Vertex index out of range");
  // New implementation aimed at speeding up this routine, since it's a
  // surprisingly large fraction of CPU time for TSTT-based swapping ---
  // like around 10% of the total before the switch, shockingly
  // (probably even bigger without Babel overhead, so GRUMMP-native
  // swapping could be spend an even bigger fraction of time here).
  if (i == 2) {
    Vert *pV0 = m_faces[0]->getVert(0);
    Vert *pV1 = m_faces[0]->getVert(1);
    Vert *pV2 = m_faces[1]->getVert(0);
    if (pV2 == pV0 || pV2 == pV1)
      return m_faces[1]->getVert(1);
    else
      return pV2;
  }
  else {
    return (
	(m_faces[0]->getLeftCell() == static_cast<const Cell*>(this)) ?
	    m_faces[0]->getVert(i) : m_faces[0]->getVert(1 - i));
  }
  // The following is the original code, retained because it's much
  // clearer in intent than the modified version.
//   const Face *const pFBase = ppFFaces[0];
//   assert (pFBase != NULL);
//   const Vert *pVRetVal;
//   switch (i)
//     {
//     case 0:
//       if (pFBase->pCCellLeft () == this)
// 	pVRetVal = (pFBase->pVVert (0));
//       else
// 	pVRetVal = (pFBase->pVVert (1));
//       break;
//     case 1:
//       if (pFBase->pCCellLeft () == this)
// 	pVRetVal = (pFBase->pVVert (1));
//       else
// 	pVRetVal = (pFBase->pVVert (0));
//       break;
//     case 2:
//       pVRetVal = (pVVertOpposite (pFBase));
//       break;
//     default:
//       assert2 (0, "Invalid vertex index");
//       pVRetVal = (pVInvalidVert);
//       break;
//     }
//   return (pVRetVal);
}

void
TriCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  // New implementation aimed at speeding up this routine, since it's a
  // surprisingly large fraction of CPU time for ITAPS-based swapping ---
  // like around 10% of the total before the switch, shockingly
  // (probably even bigger without Babel overhead, so GRUMMP-native
  // swapping could be spending an even bigger fraction of time here).
  Vert *pV0 = m_faces[0]->getVert(0);
  Vert *pV1 = m_faces[0]->getVert(1);
  Vert *pV2 = m_faces[1]->getVert(0);
  if (m_faces[0]->getLeftCell() == static_cast<const Cell*>(this)) {
    aHandles[0] = pV0;
    aHandles[1] = pV1;
  }
  else {
    aHandles[0] = pV1;
    aHandles[1] = pV0;
  }
  if (pV2 == pV0 || pV2 == pV1)
    aHandles[2] = m_faces[1]->getVert(1);
  else
    aHandles[2] = pV2;
}

double
TriCell::calcSize() const
{
  assert(doFullCheck());
  double adTemp[3];
  calcVecSize(adTemp);
  if (m_faces[0]->getVert(0)->getSpaceDimen() == 2)
    return (fabs(*adTemp));
  else
    return (dMAG3D(adTemp));
}

bool
TriCell::isClosed() const
{
  // Check to be sure that the tet has only four vertices and that each
  // appears in exactly three faces.
  const Vert *apV[] =
    { pVInvalidVert, pVInvalidVert, pVInvalidVert };
  // If the cell is deleted, it doesn't matter whether it's closed or
  // not. 
  if (isDeleted())
    return true;
  int iFaceCount[3];
  for (int iF = 0; iF < 3; iF++) {
    const Face *pF = getFace(iF);
    assert(pF->isValid() && !pF->isDeleted());
    for (int iV = 0; iV < 2; iV++) {
      const Vert *pV = pF->getVert(iV);
      assert(pV->isValid() && !pV->isDeleted());
      int ii;
      for (ii = 0; ii < 3; ii++) {
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
      if (ii == 3)
	return (false);
    }				// Done with this vertex

  }				// Done with this face

  if (iFaceCount[0] != 2 || iFaceCount[1] != 2 || iFaceCount[2] != 2)
    return (false);
  return (true);
}

void
TriCell::calcAllDihed(double adDihed[], int * const piNDihed,
		      const bool in_degrees) const
{
  assert(adDihed != NULL);
  assert(piNDihed != NULL);
  *piNDihed = 3;
  const double * const adCoord0 = getVert(0)->getCoords();
  const double * const adCoord1 = getVert(1)->getCoords();
  const double * const adCoord2 = getVert(2)->getCoords();

  double adNorm0[] =
    { +adCoord1[1] - adCoord0[1], -adCoord1[0] + adCoord0[0] };
  double adNorm1[] =
    { +adCoord2[1] - adCoord1[1], -adCoord2[0] + adCoord1[0] };
  double adNorm2[] =
    { +adCoord0[1] - adCoord2[1], -adCoord0[0] + adCoord2[0] };

  NORMALIZE2D(adNorm0);
  NORMALIZE2D(adNorm1);
  NORMALIZE2D(adNorm2);

  double dArg0 = -dDOT2D(adNorm0, adNorm2);
  double dArg1 = -dDOT2D(adNorm0, adNorm1);
  double dArg2 = -dDOT2D(adNorm1, adNorm2);

  double dScale = in_degrees ? 180. / M_PI : 1;
  adDihed[0] = GR_acos(dArg0) * dScale; // angle @ 0
  adDihed[1] = GR_acos(dArg1) * dScale; // angle @ 1
  adDihed[2] = GR_acos(dArg2) * dScale; // angle @ 2
}

void
TriCell::calcAllSolid(double adSolid[], int * const piNSolid,
		      const bool in_degrees) const
{
  calcAllDihed(adSolid, piNSolid, in_degrees);
}

void
TriCell::calcCircumcenter(double adRes[]) const
{
  const double * const adA = getVert(0)->getCoords();
  const double * const adB = getVert(1)->getCoords();
  const double * const adC = getVert(2)->getCoords();

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

  // The following is commented out because the code is known to work
  // properly, but on some machine/OS combinations, the assertions can
  // get triggered.
//  #ifndef NDEBUG
//    double dComp = dDIST2D (adRes, adA);
//    double dTest = dDIST2D (adRes, adB);
//    assert (0 == iFuzzyComp (dComp, dTest));
//    dTest = dDIST2D (adRes, adC);
//    assert (0 == iFuzzyComp (dComp, dTest));
//  #endif
}

double
TriCell::calcCircumradius() const
{
  // 14-06-99 - Charles Boivin
  // Modified this to use only 2 dimensions
  // since it's in TriCell.
  assert(isValid());
  if (checkOrient2D(getVert(0), getVert(1), getVert(2)) == 0) {
    // Tri is totally flat.
    return LARGE_DBL;
  }
  else {
    double adTemp[2];
    calcCircumcenter(adTemp);
    const double * const adCoord0 = getVert(0)->getCoords();
    return dDIST2D(adTemp, adCoord0);
  }
}

bool
TriCell::circumscribesPoint(const double adPoint[]) const
{
  // Somewhat to my surprise, using the exact predicate disambiguates cases
  // where the input includes collections of points that are co-circular but
  // whose locations have been rounded off.  (cyl and branch are two examples.)
  //   CFO-G, May 2013
  return (incircle_shew(getVert(0)->getCoords(), getVert(1)->getCoords(),
			getVert(2)->getCoords(), adPoint) > 0);
//  return (1 == isInCircle(getVert(0)->getCoords(), getVert(1)->getCoords(),
//			 getVert(2)->getCoords(), adPoint));
}

void
TriCell::calcContainmentCenter(double adRes[]) const
{
  calcCircumcenter(adRes);
}

double
TriCell::calcContainmentRadius() const
{
  assert(isValid());
  double adTemp[3];
  calcContainmentCenter(adTemp);
  adTemp[0] -= getVert(0)->x();
  adTemp[1] -= getVert(0)->y();
  adTemp[2] -= getVert(0)->z();
  return (dMAG3D(adTemp));
}

void
TriCell::calcBarycentricCoords(const double adPoint[], double adBary[]) const
// This is written for triangles in 2D.  Triangles in 3D should
// probably be projected into 2D, with assertions about coplanarity.
{
  assert(isValid());

  const Vert *pVA = getVert(0);
  const Vert *pVB = getVert(1);
  const Vert *pVC = getVert(2);

  assert(pVA->getSpaceDimen() == 2);

  double adRow1[] =
    { 1., 1., 1. };
  double adRow2[] =
    { pVA->x(), pVB->x(), pVC->x() };
  double adRow3[] =
    { pVA->y(), pVB->y(), pVC->y() };

  solve3By3(adRow1, adRow2, adRow3, 1., adPoint[0], adPoint[1], adBary);
}

bool
TriCell::isSmallAngleAffected() const
{

  int num_shell = 0;
  int num_small = 0;

  //If the triangle has a small angle corner, don't split
  for (int i = 0; i < 3; i++) {
    if (getVert(i)->isSmallAngleVert())
      num_small++;
    else if (getVert(i)->isShellVert())
      num_shell++;
  }

  if (num_small == 1 && num_shell == 2) {
    // check to see if we are actually at a small angle,
    // could be on other side (think of a meshing an interior Z shape or an X)
    // should only check if an internal boundary is affected
    for (int i = 0; i < 3; i++) {
      if (getVert(i)->getVertType() == Vert::eBdryTwoSide) {
	double angles[3];
	int numDihed;
	calcAllDihed(angles, &numDihed, false);
	for (int j = 0; j < 3; j++)
	  if (angles[j] > M_PI / 3. && getVert(j)->isSmallAngleVert())
	    return false;
      }
    }
    return true;
  }

  //If the three vertices are shell verts, don't split
  if (num_shell == 3)
    return true;

  //If two vertices are shell verts, then split if
  //the face formed by these shell verts subtend a
  //small angle.
  else if (num_shell == 2 && num_small == 0) {

    assert(
	!getVert(0)->isSmallAngleVert() && !getVert(1)->isSmallAngleVert()
	    && !getVert(2)->isSmallAngleVert());

    for (int i = 0; i < 3; i++) {
      const Face* face = getFace(i);
      if (face->getVert(0)->isShellVert() && face->getVert(1)->isShellVert()) {
	Cell* cell = face->getOppositeCell(this);
	if (dynamic_cast<TriCell*>(cell)) {
	  Vert* v0 = cell->getVert(0);
	  Vert* v1 = cell->getVert(1);
	  Vert* v2 = cell->getVert(2);
	  if (v0->isSmallAngleVert() || v1->isSmallAngleVert()
	      || v2->isSmallAngleVert())
	    return true;
	  else if (v0->isShellVert() && v1->isShellVert() && v2->isShellVert())
	    return true;
	  else
	    i = 3;
	}
      }
    }

  }

  return false;

}

int
TriCell::getLargestLayernumber()
{
  assert(this->isValid());
  int layermax = 0;

  layermax = getFace(0)->layerindex;
  layermax =
      layermax > getFace(1)->layerindex ? layermax : getFace(1)->layerindex;
  layermax =
      layermax > getFace(2)->layerindex ? layermax : getFace(2)->layerindex;
//			layermax = getVert(0)->layerindex;
//			layermax = layermax > getVert(1)->layerindex ? layermax : getVert(1)->layerindex;
//			layermax = layermax > getVert(2)->layerindex ? layermax : getVert(2)->layerindex;

  return layermax;
}

int
TriCell::getSmallestLayernumber()
{
  assert(this->isValid());
  int layermin = 1e5;
  layermin = getFace(0)->layerindex;
  layermin =
      layermin < getFace(1)->layerindex ? layermin : getFace(1)->layerindex;
  layermin =
      layermin < getFace(2)->layerindex ? layermin : getFace(2)->layerindex;

//		layermin = getVert(0)->layerindex;
//		layermin = layermin < getVert(1)->layerindex ? layermin : getVert(1)->layerindex;
//		layermin = layermin < getVert(2)->layerindex ? layermin : getVert(2)->layerindex;
  return layermin;
}

double
TriCell::calcShortestEdgeLength() const
{
  register double dDX, dDY;
  register const Vert *pV0 = getVert(0);
  register const Vert *pV1 = getVert(1);
  register const Vert *pV2 = getVert(2);

  dDX = pV0->x() - pV1->x();
  dDY = pV0->y() - pV1->y();
  double dLen01_Sq = dDX * dDX + dDY * dDY;

  dDX = pV0->x() - pV2->x();
  dDY = pV0->y() - pV2->y();
  double dLen02_Sq = dDX * dDX + dDY * dDY;

  dDX = pV1->x() - pV2->x();
  dDY = pV1->y() - pV2->y();
  double dLen12_Sq = dDX * dDX + dDY * dDY;

  double dShortest_Sq = min(min(dLen01_Sq, dLen02_Sq), dLen12_Sq);
  return (sqrt(dShortest_Sq));
}
bool
TriCell::isPointOnEdge(const double adPoint[], GR_index_t & iV0,
		       GR_index_t & iV1)
{
  double adBary[3];
  calcBarycentricCoords(adPoint, adBary);
  int i, iLowBary = 0;
  for (i = 0; i <= 2; i++)
    if (adBary[i] < 1e-10 && adBary[i] > -1e-10) {
      adBary[i] = 0;
      iLowBary++;
    }
  if (iLowBary == 1) {

    if (adBary[0] < 1e-10) {
      iV0 = 2;
      iV1 = 1;
    }
    else if (adBary[1] < 1e-10) {
      iV0 = 2;
      iV1 = 1;
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
TriCell::calcVecSize(double adRes[]) const
{
  if (getVert(0)->getSpaceDimen() == 2) {
    calcNormal2D(getVert(0)->getCoords(), getVert(1)->getCoords(),
		 getVert(2)->getCoords(), adRes);
    adRes[0] *= 0.5;
  }
  else {
    calcNormal(getVert(0)->getCoords(), getVert(1)->getCoords(),
	       getVert(2)->getCoords(), adRes);
    vSCALE3D(adRes, 0.5);
  }
}

//////////////////
