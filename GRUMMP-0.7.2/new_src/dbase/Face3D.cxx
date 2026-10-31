#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Geometry.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Vec.h"
#include "GR_VolMesh.h"

void
QuadFace::calcNormal(double * const adNorm) const
// Checks to see whether verts are in 2D or 3D.  For 2D, returns
// with the normal (directed area) in adNorm[0] (this allows
// vNormal to be called with the address of a double as its arg).
// For 3D, returns a true vector normal.
{
  assert(isValid());
  assert(m_verts[0]->getSpaceDimen() == 3);

  double adV1[3], adV2[3];
  ::calcNormal(m_verts[0]->getCoords(), m_verts[1]->getCoords(),
	       m_verts[2]->getCoords(), adV1);
  ::calcNormal(m_verts[0]->getCoords(), m_verts[2]->getCoords(),
	       m_verts[3]->getCoords(), adV2);
  adNorm[0] = (adV1[0] + adV2[0]) * 0.5;
  adNorm[1] = (adV1[1] + adV2[1]) * 0.5;
  adNorm[2] = (adV1[2] + adV2[2]) * 0.5;
}

void
QuadFace::calcUnitNormal(double * const adUnitNorm) const
// Checks to see whether verts are in 2D or 3D.  For 2D, returns
// with the unitNormal (directed area) in adUnitNorm[0] (this allows
// vUnitNormal to be called with the address of a double as its arg).
// For 3D, returns a true vector unitNormal.
{
  assert(isValid());
  assert(m_verts[0]->getSpaceDimen() == 3);

  double adV1[3], adV2[3];
  ::calcNormal(m_verts[0]->getCoords(), m_verts[1]->getCoords(),
	       m_verts[2]->getCoords(), adV1);
  ::calcNormal(m_verts[0]->getCoords(), m_verts[2]->getCoords(),
	       m_verts[3]->getCoords(), adV2);
  adUnitNorm[0] = (adV1[0] + adV2[0]);
  adUnitNorm[1] = (adV1[1] + adV2[1]);
  adUnitNorm[2] = (adV1[2] + adV2[2]);
  NORMALIZE3D(adUnitNorm);
}

double
QuadFace::calcSize() const
{
  double adTmp1[] = adDIFF3D(m_verts[0]->getCoords(),
      m_verts[2]->getCoords());
  double adTmp2[] = adDIFF3D(m_verts[1]->getCoords(),
      m_verts[3]->getCoords());
  double adTmp[3];
  vCROSS3D(adTmp1, adTmp2, adTmp);
  return (0.5 * dMAG3D(adTmp));
}

int
QuadFace::doFullCheck() const
{
  if (!Face::doFullCheck())
    return 0;
  // Check whether face location is correct.
  assert(m_rightCell->isValid() && m_leftCell->isValid());
  switch (m_faceLoc)
    {
    case eBdryFace:
      if (m_rightCell->getType() != Cell::eQuadBFace
	  && m_leftCell->getType() != Cell::eQuadBFace) {
	return 0;
      }
      break;
    case eBdryTwoSide:
      if (m_rightCell->getType() != Cell::eIntQuadBFace
	  && m_leftCell->getType() != Cell::eIntQuadBFace) {
	return 0;
      }
      break;
    case eExterior:
    case eInterior:
      if ((m_rightCell->getType() != Cell::eHex
	  && m_rightCell->getType() != Cell::ePyr
	  && m_rightCell->getType() != Cell::ePrism)
	  || (m_leftCell->getType() != Cell::eHex
	      && m_leftCell->getType() != Cell::ePyr
	      && m_leftCell->getType() != Cell::ePrism)
	  || m_leftCell->getRegion() != m_rightCell->getRegion()) {
	return 0;
      }
      break;
    default:
      // For all other face types, something is wrong.
      // Catches eUnknown, ePseudoSurface, and eInvalidFaceLoc.
      return 0;
    }
  // Exit is here rather than case-by-case so that other checks can be
  // added easily at a later date.
  return 1;
}

void
TriFace::calcNormal(double * const adNorm) const
// Checks to see whether verts are in 2D or 3D.  For 2D, returns
// with the normal (directed area) in adNorm[0] (this allows
// vNormal to be called with the address of a double as its arg).
// For 3D, returns a true vector normal.
{
  assert(isValid());
  assert(m_verts[0]->getSpaceDimen() == 3);

  ::calcNormal(m_verts[0]->getCoords(), m_verts[1]->getCoords(),
	       m_verts[2]->getCoords(), adNorm);
  vSCALE3D(adNorm, 0.5);
}

void
TriFace::calcUnitNormal(double * const adUnitNorm) const
{
  assert(isValid());
  assert(m_verts[0]->getSpaceDimen() == 3);

  ::calcNormal(m_verts[0]->getCoords(), m_verts[1]->getCoords(),
	       m_verts[2]->getCoords(), adUnitNorm);
  NORMALIZE3D(adUnitNorm);
}

void
TriFace::calcCircumcenter(double adCircCent[]) const
{
  assert(isValid());
  calcCircumcenter3D(getVert(0), getVert(1), getVert(2), adCircCent);
//#ifndef NDEBUG
//  double adA[3] = {0,0,0};
//#endif
//  double adB[3] = {getVert(1)->x() - getVert(0)->x(),
//		   getVert(1)->y() - getVert(0)->y(),
//		   getVert(1)->z() - getVert(0)->z()};
//  double adC[3] = {getVert(2)->x() - getVert(0)->x(),
//		   getVert(2)->y() - getVert(0)->y(),
//		   getVert(2)->z() - getVert(0)->z()};
//
//
//  double adRow1[3], adRow2[3], adRow3[3];
//
//  adRow1[0] = adB[0];
//  adRow1[1] = adB[1];
//  adRow1[2] = adB[2];
//
//  adRow2[0] = adC[0];
//  adRow2[1] = adC[1];
//  adRow2[2] = adC[2];
//
//  vCROSS3D(adRow1, adRow2, adRow3);
//
//  double dRHS1 = 0.5 * dDOT3D(adB, adB);
//  double dRHS2 = 0.5 * dDOT3D(adC, adC);
//  double dRHS3 = 0;
//
//  solve3By3(adRow1, adRow2, adRow3, dRHS1, dRHS2, dRHS3, adCircCent);
//
//#ifndef NDEBUG
//  double dMagA = dDIST3D(adCircCent, adA);
//  double dMagB = dDIST3D(adCircCent, adB);
//  double dMagC = dDIST3D(adCircCent, adC);
//  assert( 0 == iFuzzyComp(dMagA, dMagB));
//  assert( 0 == iFuzzyComp(dMagA, dMagC));
//#endif
//
//  adCircCent[0] += getVert(0)->x();
//  adCircCent[1] += getVert(0)->y();
//  adCircCent[2] += getVert(0)->z();

}

CubitVector
TriFace::calcCircumcenter() const
{

  double circum[] =
    { 0., 0., 0. };
  calcCircumcenter(circum);

  return CubitVector(circum[0], circum[1], circum[2]);

}

double
TriFace::calcCircumradius() const
{
  assert(isValid());
  double adCent[3];
  calcCircumcenter(adCent);

  double dDistance = dDIST3D(getVert(0)->getCoords(), adCent);
#ifndef NDEBUG
  double dDist1 = dDIST3D(getVert(1)->getCoords(), adCent);
  double dDist2 = dDIST3D(getVert(2)->getCoords(), adCent);
  assert(0 == iFuzzyComp(dDistance, dDist1));
  assert(0 == iFuzzyComp(dDistance, dDist2));
#endif
  return dDistance;
}

double
TriFace::calcSize() const
{
  double adTmp1[] = adDIFF3D(m_verts[0]->getCoords(),
      m_verts[1]->getCoords());
  double adTmp2[] = adDIFF3D(m_verts[0]->getCoords(),
      m_verts[2]->getCoords());
  double adTmp[3];
  vCROSS3D(adTmp1, adTmp2, adTmp);
  return (0.5 * dMAG3D(adTmp));
}

bool
TriFace::isLocallyDelaunay() const
// Return true if you would -not- swap this face when doing Delaunay
// face swapping.
{
  Vert *pVA, *pVB, *pVC, *pVD, *pVE, *pVPivot0, *pVPivot1, *pVOther;
  TetCell* apTCTets[4];
  int iNTets;
  eFaceCat eFC = categorizeFace(pVA, pVB, pVC, pVD, pVE, apTCTets, iNTets,
				pVPivot0, pVPivot1, pVOther);
  switch (eFC)
    {
    case eT23:
    case eT32:
    case eT22:
    case eT44:
      {
	int iRes = isInsphere(pVA, pVB, pVC, pVD, pVE);
	if (eFC == eT32)
	  iRes *= -1;
	return (iRes != 1);
      }
    default:
      return true;
    }
}

int
TriFace::doFullCheck() const
{
  if (!Face::doFullCheck())
    return 0;

  // Check whether face location is correct.
  assert(m_rightCell->isValid() && m_leftCell->isValid());
  switch (m_faceLoc)
    {
    case eBdryFace:
      if (m_rightCell->getType() != Cell::eTriBFace
	  && m_leftCell->getType() != Cell::eTriBFace) {
	return 0;
      }
      break;
    case eBdryTwoSide:
      if (m_rightCell->getType() != Cell::eIntTriBFace
	  && m_leftCell->getType() != Cell::eIntTriBFace) {
	return 0;
      }
      break;
    case eExterior:
    case eInterior:
      if ((m_rightCell->getType() != Cell::eTet
	  && m_rightCell->getType() != Cell::ePyr
	  && m_rightCell->getType() != Cell::ePrism)
	  || (m_leftCell->getType() != Cell::eTet
	      && m_leftCell->getType() != Cell::ePyr
	      && m_leftCell->getType() != Cell::ePrism)
	  || m_leftCell->getRegion() != m_rightCell->getRegion()) {
	return 0;
      }
      break;
    default:
      // For all other face types, something is wrong.
      // Catches eUnknown, ePseudoSurface, and eInvalidFaceLoc.
      return 0;
    }
  // Exit is here rather than case-by-case so that other checks can be
  // added easily at a later date.
  return 1;
}

Vert *
pVCommonVert(const Face* const pF0, const Face* const pF1,
	     const Face* const pF2)
{
  Vert *apVCand0[4], *apVCand1[4], *apVCand2[4];
  pF0->getAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apVCand0));
  pF1->getAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apVCand1));
  pF2->getAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apVCand2));

  int iNCands = pF0->getNumVerts();
  for (int j = iNCands - 1; j >= 0; j--) {
    Vert *pVCand = apVCand0[j];
    bool qHasCand = false;
    for (int i = pF1->getNumVerts() - 1; i >= 0; i--) {
      Vert *pVTest = apVCand1[i];
      if (pVTest == pVCand) {
	qHasCand = true;
	break;
      }
    }
    if (!qHasCand) {
      apVCand0[j] = apVCand0[iNCands - 1];
      iNCands--;
    }
  }

  assert(iNCands == 2 || iNCands == 0);

  for (int j = iNCands - 1; j >= 0; j--) {
    Vert *pVCand = apVCand0[j];
    bool qHasCand = false;
    for (int i = pF2->getNumVerts() - 1; i >= 0; i--) {
      Vert *pVTest = apVCand2[i];
      if (pVTest == pVCand) {
	qHasCand = true;
	break;
      }
    }
    if (!qHasCand) {
      apVCand0[j] = apVCand0[iNCands - 1];
      iNCands--;
    }
  }

  assert(iNCands == 1 || iNCands == 0);

  if (iNCands == 1)
    return apVCand0[0];
  else
    return pVInvalidVert;
}
