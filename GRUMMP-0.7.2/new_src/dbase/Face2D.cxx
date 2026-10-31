#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Geometry.h"
#include "GR_Vec.h"

void
EdgeFace::calcNormal(double * const adNorm) const
{
  // For three space dimensions, the result is ambiguous, being
  // restricted only to a plane.  In this case, the vector returned will
  // have a zero Z-component if possible.

  // Two dimensions, the result is ambiguous only with respect to sign.

  assert(isValid());
  int iSpaceDim = m_verts[0]->getSpaceDimen();

  if (iSpaceDim == 2) {
    adNorm[0] = m_verts[0]->y() - m_verts[1]->y();
    adNorm[1] = m_verts[1]->x() - m_verts[0]->x();
  }
  else {
    double adTemp[] = adDIFF3D (m_verts[0]->getCoords (),
	m_verts[1]->getCoords ());
    if ((iFuzzyComp(adTemp[0], 0.) == 0) && (iFuzzyComp(adTemp[1], 0.) == 0)) {
      adNorm[0] = fabs(adTemp[2]);
      adNorm[1] = 0;
      adNorm[2] = 0;
    }
    else {
      adNorm[0] = adTemp[1];
      adNorm[1] = -adTemp[0];
      adNorm[2] = 0;
      double dMult = dMAG3D (adTemp) / dMAG3D(adNorm);
      vSCALE3D(adNorm, dMult);
    }
  }
}

void
EdgeFace::calcUnitNormal(double * const adNorm) const
{
  // For three space dimensions, the result is ambiguous, being
  // restricted only to a plane.  In this case, the vector returned will
  // have a zero Z-component if possible.

  // Two dimensions, the result is ambiguous only with respect to sign.

  assert(isValid());
  int iSpaceDim = m_verts[0]->getSpaceDimen();

  if (iSpaceDim == 2) {
    adNorm[0] = m_verts[0]->y() - m_verts[1]->y();
    adNorm[1] = m_verts[1]->x() - m_verts[0]->x();
    NORMALIZE2D(adNorm);
  }
  else {
    double adTemp[] = adDIFF3D (m_verts[0]->getCoords (),
	m_verts[1]->getCoords ());
    if (iFuzzyComp(adTemp[2], 0.) != 0) {
      adNorm[0] = adTemp[1];
      adNorm[1] = -adTemp[0];
      adNorm[2] = 0;
      NORMALIZE3D(adNorm);
    }
    else {
      adNorm[0] = 1;
      adNorm[1] = 0;
      adNorm[2] = 0;
    }
  }
}

int
EdgeFace::doFullCheck() const
{

  if (!Face::doFullCheck())
    return 0;

  // Check whether face location is correct.
  if (m_rightCell->isValid() && m_leftCell->isValid()) {
    switch (m_faceLoc)
      {
      case eBdryFace:
	if (m_rightCell->getType() != Cell::eBdryEdge
	    && m_leftCell->getType() != Cell::eBdryEdge) {
	  return 0;
	}
	break;
      case eBdryTwoSide:

	if (m_rightCell->getType() != Cell::eIntBdryEdge
	    && m_leftCell->getType() != Cell::eIntBdryEdge) {
	  return 0;
	}
	break;
      case eExterior:
      case eInterior:
      case ePseudoSurface:
	if ((m_rightCell->getType() != Cell::eTriCell
	    && m_rightCell->getType() != Cell::eQuadCell)
	    || (m_leftCell->getType() != Cell::eTriCell
		&& m_leftCell->getType() != Cell::eQuadCell)) {
	  printf("Type:   %d %d\n", m_leftCell->getType(),
		 m_rightCell->getType());
	  printf("Region: %d %d\n", m_leftCell->getRegion(),
		 m_rightCell->getRegion());
	  return 0;
	}

	if (m_leftCell->getRegion() != m_rightCell->getRegion())
	  return -1;
	break;
      default:
	// For all other face types, something is wrong.
	// Catches eUnknown, ePseudoSurface, and eInvalidFaceLoc.
	return 0;
      }
  }
  // Exit is here rather than case-by-case so that other checks can be
  // added easily at a later date.
  return 1;
}

void
EdgeFace::calcCentroid(double adLoc[]) const
{
  if (isBdryFace()) {
    if (m_leftCell->getType() == Cell::eBdryEdge
	|| m_leftCell->getType() == Cell::eIntBdryEdge) {
      m_leftCell->calcCentroid(adLoc);
    }
    else {
      assert(
	  m_rightCell->getType() == Cell::eBdryEdge
	      || m_rightCell->getType() == Cell::eIntBdryEdge);
      m_rightCell->calcCentroid(adLoc);
    }
  }
  else
    Face::calcCentroid(adLoc);
}

// The following routine is not yet needed, so it hasn't been written.
bool
EdgeFace::isLocallyDelaunay() const
{
  assert(0);
  return false;
}

double
EdgeFace::calcSize() const
{
  double dRetVal;
  if (m_verts[0]->getSpaceDimen() == 2) {
    double adTmp[] = adDIFF2D(m_verts[0]->getCoords(),
	m_verts[1]->getCoords());
    dRetVal = dMAG2D(adTmp);
  }
  else {
    double adTmp[] = adDIFF3D(m_verts[0]->getCoords(),
	m_verts[1]->getCoords());
    dRetVal = dMAG3D(adTmp);
  }
  return dRetVal;
}

