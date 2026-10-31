#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Geometry.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

const Vert*
QuadCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert2(i >= 0 && i <= 3, "Vertex index out of range");
  // Assume that faces are returned in cyclic order by pFFace().
  const Face * const pFBase = getFace(0);
  const Face * const pFOpp = getFace(2);
  switch (i)
    {
    case 0:
      if (pFBase->getLeftCell() == static_cast<const Cell*>(this)) {
	return pFBase->getVert(0);
      }
      else {
	return pFBase->getVert(1);
      }
      break;
    case 1:
      if (pFBase->getLeftCell() == static_cast<const Cell*>(this)) {
	return pFBase->getVert(1);
      }
      else {
	return pFBase->getVert(0);
      }
      break;
    case 2:
      if (pFOpp->getLeftCell() == static_cast<const Cell*>(this)) {
	return pFOpp->getVert(0);
      }
      else {
	return pFOpp->getVert(1);
      }
      break;
    case 3:
      if (pFOpp->getLeftCell() == static_cast<const Cell*>(this)) {
	return pFOpp->getVert(1);
      }
      else {
	return pFOpp->getVert(0);
      }
      break;
    default:
      return pVInvalidVert;
    }
}

void
QuadCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  for (int i = 0; i < 4; i++) {
    aHandles[i] = const_cast<Vert*>(getVert(i));
  }
}

double
QuadCell::calcSize() const
{
  assert(doFullCheck());
  double adVecSize[3];
  calcVecSize(adVecSize);
  if (getVert(1)->getSpaceDimen() == 2)
    return fabs(*adVecSize);
  else
    return dMAG3D(adVecSize);
}

void
QuadCell::calcAllDihed(double adDihed[], int * const piNDihed,
		       const bool in_degrees) const
{
  *piNDihed = 4;
  const double * const adCoord0 = getVert(0)->getCoords();
  const double * const adCoord1 = getVert(1)->getCoords();
  const double * const adCoord2 = getVert(2)->getCoords();
  const double * const adCoord3 = getVert(3)->getCoords();

  double adNorm0[] =
    { +adCoord1[1] - adCoord0[1], -adCoord1[0] + adCoord0[0] };
  double adNorm1[] =
    { +adCoord2[1] - adCoord1[1], -adCoord2[0] + adCoord1[0] };
  double adNorm2[] =
    { +adCoord3[1] - adCoord2[1], -adCoord3[0] + adCoord2[0] };
  double adNorm3[] =
    { +adCoord0[1] - adCoord3[1], -adCoord0[0] + adCoord3[0] };

  NORMALIZE2D(adNorm0);
  NORMALIZE2D(adNorm1);
  NORMALIZE2D(adNorm2);
  NORMALIZE2D(adNorm3);

  double dScale = in_degrees ? 180. / M_PI : 1;
  adDihed[0] = GR_acos(-dDOT2D(adNorm3, adNorm0)) * dScale;
  adDihed[1] = GR_acos(-dDOT2D(adNorm0, adNorm1)) * dScale;
  adDihed[2] = GR_acos(-dDOT2D(adNorm1, adNorm2)) * dScale;
  adDihed[3] = GR_acos(-dDOT2D(adNorm2, adNorm3)) * dScale;

}

void
QuadCell::calcAllSolid(double adSolid[], int * const piNSolid,
		       const bool in_degrees) const
{
  calcAllDihed(adSolid, piNSolid, in_degrees);

}

void
QuadCell::calcVecSize(double adRes[]) const
{
  if (getVert(0)->getSpaceDimen() == 2) {
    double adDiff1[] =
      { getVert(1)->x() - getVert(0)->x(), getVert(1)->y() - getVert(0)->y() };
    double adDiff2[] =
      { getVert(2)->x() - getVert(0)->x(), getVert(2)->y() - getVert(0)->y() };
    double adDiff3[] =
      { getVert(3)->x() - getVert(0)->x(), getVert(3)->y() - getVert(0)->y() };
    double dCross1 = dCROSS2D(adDiff1, adDiff2);
    double dCross2 = dCROSS2D(adDiff2, adDiff3);
    *adRes = (dCross1 + dCross2) * 0.5;
  }
  else {
    double adVec0[3], adVec1[3];
    calcNormal(getVert(0)->getCoords(), getVert(1)->getCoords(),
	       getVert(2)->getCoords(), adVec0);
    calcNormal(getVert(0)->getCoords(), getVert(2)->getCoords(),
	       getVert(3)->getCoords(), adVec1);
    adRes[0] = (adVec0[0] + adVec1[0]) * 0.5;
    adRes[1] = (adVec0[1] + adVec1[1]) * 0.5;
    adRes[2] = (adVec0[2] + adVec1[2]) * 0.5;
  }
}

bool
QuadCell::isClosed() const
{
  // Check to be sure that the tet has only four vertices and that each
  // appears in exactly three faces.
  const Vert *apV[] =
    { pVInvalidVert, pVInvalidVert, pVInvalidVert,
    pVInvalidVert };
  // If the cell is deleted, it doesn't matter whether it's closed or
  // not. 
  if (isDeleted())
    return true;
  int iFaceCount[4];
  for (int iF = 0; iF < 4; iF++) {
    const Face *pF = getFace(iF);
    assert(pF->isValid() && !pF->isDeleted());
    for (int iV = 0; iV < 2; iV++) {
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
    }				// Done with this vertex

  }				// Done with this face

  if (iFaceCount[0] != 2 || iFaceCount[1] != 2 || iFaceCount[2] != 2
      || iFaceCount[3] != 2)
    return (false);
  return (true);
}

#define vSwitchFaces(pFA, pFB) do { \
    Face *pFTmp = pFA; \
    pFA = pFB; \
    pFB = pFTmp; \
  } while (0);

void
QuadCell::canonicalizeFaceOrder()
{
  // This isn't particularly pretty, but it has to be done...
  Vert *apV[4];
  if (m_faces[0]->getLeftCell() == this
      || (m_faces[0]->getLeftCell() == pCInvalidCell
	  && m_faces[0]->getRightCell() == pCInvalidCell)) {
    apV[0] = m_faces[0]->getVert(0);
    apV[1] = m_faces[0]->getVert(1);
  }
  else {
    apV[0] = m_faces[0]->getVert(1);
    apV[1] = m_faces[0]->getVert(0);
  }

  if (!m_faces[1]->hasVert(apV[1])) {
    // Need to find out which face is actually the next in order
    if (m_faces[2]->hasVert(apV[1])) {
      vSwitchFaces(m_faces[1], m_faces[2]);
    }
    else {
      assert(m_faces[3]->hasVert(apV[1]));
      vSwitchFaces(m_faces[1], m_faces[3]);
    }
  }

  if (!m_faces[3]->hasVert(apV[0])) {
    vSwitchFaces(m_faces[2], m_faces[3]);
  }

  // Check correctness
  assert(m_faces[0]->hasVert(apV[0]) && m_faces[0]->hasVert(apV[1]));
  assert(!m_faces[1]->hasVert(apV[0]) && m_faces[1]->hasVert(apV[1]));
  assert(!m_faces[2]->hasVert(apV[0]) && !m_faces[2]->hasVert(apV[1]));
  assert(m_faces[3]->hasVert(apV[0]) && !m_faces[3]->hasVert(apV[1]));
}

