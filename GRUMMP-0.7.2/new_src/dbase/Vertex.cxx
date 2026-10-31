#include <math.h>
#include <stdlib.h>
#include "GR_Vertex.h"
#include "GR_BFace.h"
#include "GR_Mesh.h"
#include "GR_Bdry3D.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"

#if (defined(SOLARIS2) && !defined(NDEBUG))
// Solaris puts the definition of finite() in an odd place.  Need this
// only in some assertions, so don't include the file when assertions
// are null.
#include <ieeefp.h>
#endif

Face*
Vert::getHintFace() const
{
  assert(!isDeleted());
  if (m_faceHint->isValid() && m_faceHint->hasVert(this)) {
    return m_faceHint;
  }
  else {
#ifdef USE_VECTOR
    if (vecpFInc.size() > 0) {
      m_faceHint = vecpFInc[0];
#else
    if (m_faceList->size() > 0) {
      m_faceHint = m_faceList->getFace(0);
#endif
      assert(m_faceHint->isValid());
      assert(m_faceHint->hasVert(this));
      return m_faceHint;
    }
    else {
      return m_faceList->getFace(0);
//      return pFInvalidFace;
    }
  }
}

void
Vert::setHintFace(Face *pF)
{
  assert(isValid() && !isDeleted());
  assert(pF->isValid() && !pF->isDeleted());
#ifndef NDEBUG
  bool qOK = false;
  for (int ii = 0; ii < pF->getNumVerts(); ii++) {
    if (pF->getVert(ii) == this)
      qOK = true;
  }
  assert(qOK);
#endif
  m_faceHint = pF;
}

Face*
Vert::setHintFace(const Mesh* const pM)
{
  assert(isValid() && !isDeleted());

  if (m_faceHint->isValid() && !m_faceHint->isDeleted()
      && m_faceHint->hasVert(this))
    return m_faceHint;

  int i;
  // Check first among boundary faces.
  BFace *pBF;
  Face *pF;
//#pragma omp private (pBF,pF)
  for (i = pM->getNumBdryFaces() - 1; i >= 0; i--) {
    pBF = pM->getBFace(i);
    if (pBF->isValid() && !pBF->isDeleted() && pBF->hasVert(this)) {
      m_faceHint = pBF->getFace(0);
      return m_faceHint;
    }
  }
  // Failing that, try all the faces until you get a hit.
  for (i = pM->getNumFaces() - 1; i >= 0; i--) {
    pF = pM->getFace(i);
    if (pF->isValid() && !pF->isDeleted() && pF->hasVert(this)) {
      m_faceHint = pF;
      return m_faceHint;
    }
  }
  return pFInvalidFace;
}

void
Vert::updateHintFace()
{
  assert(isValid() && !isDeleted());

  Face *pF;
//#pragma omp threadprivate (pF)
  if (m_faceList->size() == 0) {
    pF = pFInvalidFace;
  }
  else {
    pF = m_faceList->getFace(0);
    assert(pF->isValid() && !pF->isDeleted() && pF->hasVert(this));
    if (pF == m_faceHint) {
      // This is bad, because the hint face is presumably about to be
      // deleted, or we wouldn't be here.
      pF = m_faceList->getFace(1);
      assert(pF->isValid() && !pF->isDeleted() && pF->hasVert(this));
    }
  }
  m_faceHint = pF;
}

Face*
findCommonFace(const Vert * const pV0, const Vert * const pV1,
	       bool qCanBeConnectedBothSides)
{
#ifdef USE_VECTOR
  std::vector<Face *>::const_iterator iter0 = pV0->vecpFInc.begin(),
  iter0_end = pV0->vecpFInc.end(),
  iter1, iter1_end = m_vert1->vecpFInc.end();

  for (; iter0 != iter0_end; iter0++) {
    for (iter1 = m_vert1->vecpFInc.begin(); iter1 != iter1_end; iter1++) {
      if (*iter0 == *iter1) return (*iter0);
    }
  }
  return pFInvalidFace;
#else
  int iMax0 = pV0->m_faceList->size();

  Face *result;
  Face *possible;
  Face *pFCand;
//#pragma omp critical(FacefindCommonFace)
//	{
  result = pFInvalidFace;
  possible = pFInvalidFace;
  for (int i0 = 0; i0 < iMax0; i0++) {
    pFCand = pV0->m_faceList->getFace(i0);
    assert(pFCand->isValid());
    assert(!pFCand->isDeleted());
    if (pFCand->hasVert(pV1)) {
      if (!pFCand->getLeftCell()->isValid()
	  || !pFCand->getRightCell()->isValid()
	  || pFCand->getType() == Face::eMultiEdge) {
	result = pFCand;
      }
      else {
	possible = pFCand;
      }
    }
  }

  if (!(result || !possible || qCanBeConnectedBothSides)) {
    std::cout << "findCommonFace Failed: " << pV0->x() << ":" << pV0->y() << " "
	<< pV1->x() << ":" << pV1->y() << std::endl;
    assert(result || !possible || qCanBeConnectedBothSides);
  }
  if (!result && qCanBeConnectedBothSides)
    result = possible;
  return result;
#endif
}

Face*
findCommonFace(const Vert * const pV0, const Vert * const pV1,
	       const Vert * const pV2, const Vert * const pV3,
	       bool qCanBeConnectedBothSides)
{
  assert(pV0->isValid());
  assert(pV1->isValid());
  assert(pV2->isValid());
#ifdef USE_VECTOR
  std::vector<Face *>::const_iterator iterA = pV0->vecpFInc.begin(),
  iterA_end = pV0->vecpFInc.end(),
  iterB, iterB_end = m_vert1->vecpFInc.end();
  Face *apFCom01[pV0->vecpFInc.size()];
  int iComCount = 0;

  Face *result = pFInvalidFace;

  for (; iterA != iterA_end; iterA++) {
    assert(*iterA != pFInvalidFace);
    if((*iterA)->isDeleted()) continue;
    for (iterB = m_vert1->vecpFInc.begin(); iterB != iterB_end; iterB++) {
      assert(*iterB != pFInvalidFace);
      if((*iterB)->isDeleted()) continue;
      if (*iterA == *iterB) {
	apFCom01[iComCount++] = *iterA;
      }
    }
  }

  // Now check which of those faces are also common to pV2.
  Face **iterCom_end = apFCom01 + iComCount;
  iterB_end = m_vert2->vecpFInc.end();
  for (Face **iterCom = apFCom01;
      iterCom != iterCom_end && !result;
      iterCom++) {
    assert(*iterCom != pFInvalidFace);
    if((*iterCom)->isDeleted()) continue;
    for (iterB = m_vert2->vecpFInc.begin();
	iterB != iterB_end && !result;
	iterB++) {
      if (*iterCom == *iterB) {
	result = (*iterCom);
      }
    }
  }
#else
  int iMax0 = pV0->m_faceList->size();
  int iMax1 = pV1->m_faceList->size();
  int iMax2 = pV2->m_faceList->size();

  Face *result = pFInvalidFace;
  Face *possible = pFInvalidFace;

  if (iMax0 < iMax1 && iMax0 < iMax2) {
    for (int i0 = 0; i0 < iMax0; i0++) {
      Face *pFCand = pV0->m_faceList->getFace(i0);
      if ((pFCand->getNumVerts() == 4 && !pV3->isValid())
	  || (pFCand->getNumVerts() == 3 && pV3->isValid()))
	continue;
      assert(pFCand != pFInvalidFace);
      if (!pFCand->isDeleted() && pFCand->hasVert(pV1)
	  && pFCand->hasVert(pV2)) {
	if (pFCand->getLeftCell()->isValid()
	    && pFCand->getRightCell()->isValid()) {
	  possible = pFCand;
	}
	else {
	  result = pFCand;
	}
      }
    }
  }
  else if (iMax1 < iMax2) {
    for (int i1 = 0; i1 < iMax1; i1++) {
      Face *pFCand = pV1->m_faceList->getFace(i1);
      if ((pFCand->getNumVerts() == 4 && !pV3->isValid())
	  || (pFCand->getNumVerts() == 3 && pV3->isValid()))
	continue;
      assert(pFCand != pFInvalidFace);
      if (!pFCand->isDeleted() && pFCand->hasVert(pV0)
	  && pFCand->hasVert(pV2)) {
	if (pFCand->getLeftCell()->isValid()
	    && pFCand->getRightCell()->isValid()) {
	  possible = pFCand;
	}
	else {
	  result = pFCand;
	}
      }
    }
  }
  else {
    // Vert 2 has the fewest neighbors.
    for (int i2 = 0; i2 < iMax2; i2++) {
      Face *pFCand = pV2->m_faceList->getFace(i2);
      if ((pFCand->getNumVerts() == 4 && !pV3->isValid())
	  || (pFCand->getNumVerts() == 3 && pV3->isValid()))
	continue;
      assert(pFCand != pFInvalidFace);
      if (!pFCand->isDeleted() && pFCand->hasVert(pV1)
	  && pFCand->hasVert(pV0)) {
	if (pFCand->getLeftCell()->isValid()
	    && pFCand->getRightCell()->isValid()) {
	  possible = pFCand;
	}
	else {
	  result = pFCand;
	}
      }
    }
  }

  // Prefer a result that has only one cell attached to it.
  if (!(result || !possible || qCanBeConnectedBothSides)) {
    std::cout << pV0->x() << ":" << pV0->y() << ":" << pV0->z() << std::endl;
    std::cout << pV1->x() << ":" << pV1->y() << ":" << pV1->z() << std::endl;
    std::cout << pV2->x() << ":" << pV2->y() << ":" << pV2->z() << std::endl;
    assert(result || !possible || qCanBeConnectedBothSides);
  }
  if (!result && qCanBeConnectedBothSides)
    result = possible;
#endif

  assert(
      !result->isValid()
	  || (result->hasVert(pV0) && result->hasVert(pV1)
	      && result->hasVert(pV2)));
  if (pV3->isValid() && result->isValid()) {
    bool qOK = result->hasVert(pV3);
    if (!qOK) {
      // Trying to create a quad, but hit a snag...
      switch (result->getType())
	{
	case Face::eTriFace:
	  vFatalError("Tried to identify a quad, but three verts form a tri.",
		      "findCommonFace");
	  break;
	case Face::eQuadFace:
	  vWarning(
	      "Tried to identify a quad, but another exists with three of these verts.");
	  break;
	case Face::eEdgeFace:
	case Face::eMultiEdge:
	  assert(0);
	  // Should never be able to ID an edge w/ three verts in common.
	  break;
	}
    }
  } // Done confirming face identity
  return result;
}

void
Vert::setDefaultFlags()
{
  m_dim = 0;
  m_vertType = eInterior;
  m_del = m_delRequested = m_active = m_struct = m_small = m_shell = false;
  m_lenScale = LARGE_DBL;
}

void
Vert::resetAllData()
{
//#pragma omp critical(VertresetAllDatacheck)
//{
  setDefaultFlags();
  m_faceHint = pFInvalidFace;
#ifdef USE_VECTOR
  vecpFInc.clear();
#else
  if (m_faceList)
    delete m_faceList;
  m_faceList = new FaceList();
#endif
  clearLengthScale();
  m_loc[0] = m_loc[1] = m_loc[2] = 0;
//}
}

void
Vert::copyAllFlagsFrom(const Vert * const V)
{
  m_dim = V->m_dim;
  m_vertType = V->m_vertType;
  m_del = V->m_del;
  m_delRequested = V->m_delRequested;
  m_active = V->m_active;
  m_struct = V->m_struct;
  m_small = V->m_small;
  m_shell = V->m_shell;
}

bool
Vert::isBdryVert() const
{
  bool qRetVal = false;
  qRetVal = ((m_vertType == eBdryApex) || (m_vertType == eBdryCurve)|| (m_vertType == eBdryTwoSide) || (m_vertType == eBdry)|| (m_vertType == eBdryPhantom));
  return qRetVal;
}

bool
Vert::isSmoothable() const
{
  return (!isStructured() && !isDeletionRequested() && !isDeleted() && !(getVertType() == eBdryApex) && !(getVertType() == ePseudoSurface) && !(getVertType() == eBdryCurve));
  // Technically, bdry curve points -are- smoothable, but that isn't
  // implemented yet.
}

void
Vert::setCoords(const int iNumDim, const double adNewLoc[])
{
  assert(iNumDim == 2 || iNumDim == 3);
  if (getSpaceDimen() == 0)
    setSpaceDimen(iNumDim);
  // The following assertion is false for ITAPS setting 2D coords...
#ifndef ITAPS
  assert(iNumDim == getSpaceDimen());
#endif

  assert(finite(adNewLoc[0]));
  m_loc[0] = adNewLoc[0];
  assert(finite(adNewLoc[1]));
  m_loc[1] = adNewLoc[1];
  if (iNumDim == 3) {
    assert(finite(adNewLoc[2]));
    m_loc[2] = adNewLoc[2];
  }
}

void
Vert::setCoords(const char * const pcBuffer)
{
  switch (getSpaceDimen())
    {
    case 2:
      sscanf(pcBuffer, "%lf %lf", &(m_loc[0]), &(m_loc[1]));
      assert(finite(m_loc[0]));
      assert(finite(m_loc[1]));
      break;
    case 3:
      sscanf(pcBuffer, "%lf %lf %lf", &(m_loc[0]), &(m_loc[1]), &(m_loc[2]));
      assert(finite(m_loc[0]));
      assert(finite(m_loc[1]));
      assert(finite(m_loc[2]));
      break;
    default:
      assert(0);
      break;
    }
}

void
Vert::setCoords(const int iNumDim, const CubitVector& coord)
{
  assert(iNumDim == 2 || iNumDim == 3);
  if (getSpaceDimen() == 0)
    setSpaceDimen(iNumDim);
  assert(iNumDim == getSpaceDimen());
  assert(finite(coord.x()));
  m_loc[0] = coord.x();
  assert(finite(coord.y()));
  m_loc[1] = coord.y();
  if (iNumDim == 3) {
    assert(finite(coord.z()));
    m_loc[2] = coord.z();
  }
}

CubitBox
Vert::bounding_box() const
{

  CubitVector coord1;
  CubitVector coord2;

  if (getSpaceDimen() == 2) {
    coord1.set(x(), y(), 0.);
    coord2.set(x(), y(), 0.);
  }
  else {
    coord1.set(x(), y(), z());
    coord2.set(x(), y(), z());
  }

  return CubitBox(coord1, coord2);

}

double
calcDistanceBetween(const Vert * const pV0, const Vert * const pV1)
{
  assert(pV0->getSpaceDimen() == pV1->getSpaceDimen());
  if (pV0->getSpaceDimen() == 2)
    return dDIST2D(pV0->getCoords(), pV1->getCoords());
  else {
    assert(pV0->getSpaceDimen() == 3);
    return dDIST3D(pV0->getCoords(), pV1->getCoords());
  }
}

void
Vert::copyConnectivityFrom(const Vert& V)
{
  m_faceHint = V.m_faceHint; // This hint will not be valid unless the face
  // is updated, too.
  (*m_faceList) = *(V.m_faceList);
  for (int i = getNumFaces() - 1; i >= 0; i--) {
    Face* pF = getFace(i);
    pF->replaceVert(&V, this);
  }
}

void
Vert::copyLengthInfoFrom(const Vert& V)
{
  // TODO Technically, should copy metric here.  But the metric may or may not
  // stay part of the Vert class.
  m_lenScale = V.m_lenScale;
  parent_entity = V.parent_entity;
  parent_topology = V.parent_topology;
  if (V.m_useMetric) {
    m_useMetric = true;
    for (int i = 0; i < m_dim * m_dim - 1; i++)
      m_metric[i] = V.m_metric[i];
  }
}

void
Vert::copyCoordsFrom(const Vert& V)
{
  m_loc[0] = V.m_loc[0];
  m_loc[1] = V.m_loc[1];
  m_loc[2] = V.m_loc[2];
}

Vert&
Vert::operator=(const Vert& V)
{
  if (this != &V) {
    copyCoordsFrom(V);
    copyConnectivityFrom(V);
    copyAllFlagsFrom(&V);
    // TODO Technically, should copy metric here.  But the metric may or may not
    // stay part of the Vert class.
    copyLengthInfoFrom(V);
  }
  return (*this);
}
void
Vert::printVertInfo() const
{
  char vtype[10][20] =
    { "Unknown", "BdryApex", "BdryCurve", "BdryTwoSide", "Bdry",
	"PseudoSurface", "Interior", "InteriorFixed", "BBox", "Invalid" };

  if (m_dim == 2)
    logMessage(1, "%p %.15f %.15f %s\n", this, m_loc[0], m_loc[1],
	       vtype[m_vertType]);
  else
    logMessage(1, "%p %.15f %.15f %.15f %d %s\n", this, m_loc[0], m_loc[1],
	       m_loc[2], isDeleted(), vtype[m_vertType]);

}
