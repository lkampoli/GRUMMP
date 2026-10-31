#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Geometry.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

Face::Face(Cell* const pCA, Cell* const pCB, const int iNV, Vert* const pV0,
	   Vert* const pV1, Vert* const pV2, Vert* const pV3) :
    m_leftCell(pCA), m_rightCell(pCB), m_length(-1), m_badFace(false), m_doNotDel(
	false), m_faceLoc(eInterior), m_restricted(false), m_sameShape(false), m_deleted(
	false), m_locked(false)
// Other stuff initialized in SetDefaultFlags and SetFaceLoc
{
  assert(iNV > 0 && iNV <= 4);
  m_verts[0] = pV0;
  m_verts[1] = pV1;
  m_verts[2] = pV2;
  m_verts[3] = pV3;
  setFaceLoc();
}

/// This function makes it so that a face has no connectivity and
/// default flags and location.  This is needed when re-using a Face
/// stored in an array, which would otherwise inherit garbage data.
void
Face::resetAllData()
{
//#pragma omp critical(resetAllDataFace)
//{
  m_leftCell = m_rightCell = pCInvalidCell;
  m_verts[0] = m_verts[1] = m_verts[2] = m_verts[3] = pVInvalidVert;
  setDefaultFlags();
  setFaceLoc();
//}
}

void
Face::copyAllFlags(const Face &F)
{
  m_faceLoc = F.m_faceLoc;
  m_restricted = F.m_restricted;
  m_sameShape = F.m_sameShape;
  m_deleted = F.m_deleted;
  m_locked = F.m_locked;
}

void
Face::setDefaultFlags()
{
  m_faceLoc = eInterior;
  m_sameShape = m_deleted = m_restricted = m_locked = false;
}

Face&
Face::operator=(Face& F)
{
  assert(getType() == F.getType());
  if (this != &F) {
    if (!isDeleted()) {
      // The face *this will no longer be connected to its previous
      // cells. Only do this if the face hasn't been marked as deleted
      // from the mesh already.
      if (m_leftCell->isValid() && m_leftCell->hasFace(this))
	m_leftCell->removeFace(this);
      if (m_rightCell->isValid() && m_rightCell->hasFace(this))
	m_rightCell->removeFace(this);
    }

    // Copy all the flag information for the face.
    copyAllFlags(F);

    // Not only copy the cells, but also make sure that the face
    // pointers for the cells are updated properly.
    m_leftCell = F.m_leftCell;
    if (m_leftCell->isValid())
      m_leftCell->replaceFace(&F, this);
    m_rightCell = F.m_rightCell;
    if (m_rightCell->isValid())
      m_rightCell->replaceFace(&F, this);

    // Copy the vertex pointers and make sure that the vertex hints are
    // updated properly.
    for (int i = 0; i < getNumVerts(); i++) {
      m_verts[i] = F.m_verts[i];
      if (m_verts[i]->isValid() && !m_verts[i]->isDeleted()) {
	m_verts[i]->setHintFace(this);
	m_verts[i]->removeFace(&F);
	m_verts[i]->addFace(this);
      }
    }
  }
  assert(!isDeleted());
  m_length = F.m_length;
  m_badFace = F.m_badFace;
  m_doNotDel = F.m_doNotDel;
  return (*this);
}

int
Face::doFullCheck() const
{

  if (!isValid() || isDeleted())
    return 0;
  int checking = 1;
//#pragma omp critical (dofullcheck)
//  {
  for (int iV = getNumVerts() - 1; iV >= 0; iV--) {
    const Vert *pV = getVert(iV);
    if (!pV->hasFace(this))
      checking *= 0;
  }
//  }
  if (checking == 0)
    return 0;
  if ((m_leftCell->isValid()
      && (m_leftCell->isDeleted() || !m_leftCell->hasFace(this)))
      || (m_rightCell->isValid()
	  && (m_rightCell->isDeleted() || !m_rightCell->hasFace(this)))) {
    return 0;
  }
  else
    return 1;
}

bool
Face::isBdryFace() const
{
  assert(isValid());

  // Modified 07-Sep-99 by Charles Boivin
  // Face location flags should now be set properly, so no need to
  // compute it anymore..
  //
  // 28-Jan-00:  Does this work properly for 3D?  CFO-G

  if (m_faceLoc == Face::eUnknown) {
    if (m_leftCell->doFullCheck()) {
      switch (m_leftCell->getType())
	{
	case Cell::eBdryEdge:
	case Cell::eTriBFace:
	case Cell::eQuadBFace:
	  m_faceLoc = Face::eBdryFace;
	  break;
	case Cell::eIntBdryEdge:
	case Cell::eIntTriBFace:
	case Cell::eIntQuadBFace:
	  m_faceLoc = Face::eBdryTwoSide;
	  break;
	default:
	  // Check the other cell
	  if (m_rightCell->doFullCheck()) {
	    switch (m_rightCell->getType())
	      {
	      case Cell::eBdryEdge:
	      case Cell::eTriBFace:
	      case Cell::eQuadBFace:
		m_faceLoc = Face::eBdryFace;
		break;
	      case Cell::eIntBdryEdge:
	      case Cell::eIntTriBFace:
	      case Cell::eIntQuadBFace:
		m_faceLoc = Face::eBdryTwoSide;
		break;
	      default:
		m_faceLoc = Face::eInterior;
		break;
	      } // end of inner switch
	  } // end of pCR handling
	  break;
	} // end of outer switch
    } // end of pCL handling
  }

  // FIX ME: This shouldn't need to be computed every time.
  //	int iLoc = iFaceLoc();

  if ((m_faceLoc == Face::eBdryFace) || (m_faceLoc == Face::eBdryTwoSide))
    return true;
  else
    return false;
}

bool
Face::isSwapAllowed() const
{
  assert(isValid());
  // These clauses are ordered so that the  most likely hits are first
  // (I think?).  Of course, they're mostly inline functions that access
  // boolean member functions, so they're all quick.

  if (isDeleted() || isBdryFace() || isLocked() ||(getFaceLoc() == Face::ePseudoSurface) || (getType() != Face::eTriFace && getType() != Face::eEdgeFace))
    return false;
  else
    return true;
}

/// Deprecated for external use, though it may survive in Cell::operator=()
void
Face::addCell(Cell* const pCNew)
{
  assert(isValid());
  assert2(
      (m_leftCell == pCInvalidCell || m_rightCell == pCInvalidCell) && (m_leftCell != pCInvalidCell || m_rightCell != pCInvalidCell),
      "Need exactly one invalid cell");
  if (m_leftCell == pCInvalidCell)
    m_leftCell = pCNew;
  else
    m_rightCell = pCNew;
  setFaceLoc();
}

/// Deprecated for external use, though it may survive in Cell::operator=()
void
Face::removeCell(const Cell* const pCOld)
{

  assert2((m_leftCell == pCOld || m_rightCell == pCOld),
	  "Cell must bound the face");
  if (m_leftCell == pCOld)
    m_leftCell = pCInvalidCell;
  else {
    assert(m_rightCell == pCOld);
    m_rightCell = pCInvalidCell;
  }
}

/// Deprecated for external use, though it may survive in Cell::operator=()
void
Face::replaceCell(const Cell* const pCOld, Cell* const pCNew)
{
  assert(isValid());
  assert(pCNew->isValid());
  assert2((m_leftCell == pCOld || m_rightCell == pCOld),
	  "Cell must bound the face");
  if (m_leftCell == pCOld)
    m_leftCell = pCNew;
  else
    m_rightCell = pCNew;
  setFaceLoc();
}

void
Face::replaceVert(const Vert* const pVOld, Vert* const pVNew)
{
  assert(isValid());
  assert(pVNew->isValid());
  int i;
  for (i = 0; i < getNumVerts(); i++) {
    if (pVOld == m_verts[i])
      break;
  }
  assert(i < getNumVerts());
  assert(m_verts[i] == pVOld);
  m_verts[i] = pVNew;
}

void
Face::interchangeCells(void)
{
  assert(isValid());
  Cell * const pCTemp = m_leftCell;
  m_leftCell = m_rightCell;
  m_rightCell = pCTemp;
}

void
Face::interchangeCellsAndVerts(void)
{
  assert(isValid());
  interchangeCells();
  Vert * apVTmp[4];
  int iV, iNV = getNumVerts();
  for (iV = 0; iV < iNV; iV++) {
    apVTmp[iV] = m_verts[iV];
  }
  for (iV = 0; iV < iNV; iV++) {
    m_verts[iV] = apVTmp[iNV - iV - 1];
  }
}

void
Face::setFaceLoc()
{
  int iLType = Cell::eInvalid;
  int iRType = Cell::eInvalid;
  if (m_leftCell->isValid())
    iLType = m_leftCell->getType();
  if (m_rightCell->isValid())
    iRType = m_rightCell->getType();
  if (iLType == Cell::eIntBdryEdge || iLType == Cell::eIntTriBFace
      || iLType == Cell::eIntQuadBFace || iRType == Cell::eIntBdryEdge
      || iRType == Cell::eIntTriBFace || iRType == Cell::eIntQuadBFace) {
    m_faceLoc = eBdryTwoSide;
  }
  else if (iLType == Cell::eBdryEdge || iLType == Cell::eTriBFace
      || iLType == Cell::eQuadBFace || iRType == Cell::eBdryEdge
      || iRType == Cell::eTriBFace || iRType == Cell::eQuadBFace) {
    m_faceLoc = eBdryFace;
  }
  else {
    m_faceLoc = eInterior;
  }
}

void
Face::projectOntoFace(double adPoint[]) const
{
  double adNorm[3];
  calcUnitNormal(adNorm);

  int iDim = getVert(0)->getSpaceDimen();
  if (iDim == 2) {
    double adTemp[] =
      { adPoint[0] - getVert(0)->x(), adPoint[1] - getVert(0)->y() };
    double dDistToFace = dDOT2D(adTemp, adNorm);
    adPoint[0] -= adNorm[0] * dDistToFace;
    adPoint[1] -= adNorm[1] * dDistToFace;
  }
  else {
    double adTemp[] =
      { adPoint[0] - getVert(0)->x(), adPoint[1] - getVert(0)->y(), adPoint[2]
	  - getVert(0)->z() };
    double dDistToFace = dDOT3D(adTemp, adNorm);
    adPoint[0] -= adNorm[0] * dDistToFace;
    adPoint[1] -= adNorm[1] * dDistToFace;
    adPoint[2] -= adNorm[2] * dDistToFace;
  }
}

void
Face::calcCentroid(double adCent[]) const
{
  assert(isValid());
  if (isBdryFace()) {
    // Use the BFace centroid routine so that the result lies on the
    // boundary.
    if (m_leftCell->isBdryCell()) {
      m_leftCell->calcCentroid(adCent);
    }
    else {
      assert(m_rightCell->isBdryCell());
      m_rightCell->calcCentroid(adCent);
    }
  }
  else {
    int iDim = getVert(0)->getSpaceDimen();
    switch (iDim)
      {
      case 2:
	// Always okay in 2D.
	{
	  adCent[0] = adCent[1] = 0;
	  for (int i = 0; i < getNumVerts(); i++) {
	    adCent[0] += getVert(i)->x();
	    adCent[1] += getVert(i)->y();
	  }
	  adCent[0] /= getNumVerts();
	  adCent[1] /= getNumVerts();
	}
	break;
      case 3:
	{
	  assert(getType() == eTriFace);
	  adCent[0] = adCent[1] = adCent[2] = 0;
	  for (int i = 0; i < getNumVerts(); i++) {
	    adCent[0] += getVert(i)->x();
	    adCent[1] += getVert(i)->y();
	    adCent[2] += getVert(i)->z();
	  }
	  adCent[0] /= getNumVerts();
	  adCent[1] /= getNumVerts();
	  adCent[2] /= getNumVerts();
	}
	break;
      default:
	assert2(0, "Unknown number of dimensions");
	break;
      }
  }
}

Cell*
findCommonCell(const Face* const pF0, const Face* const pF1)
{
  // Modified 1/27/11 so that it will work properly with MultiEdges.
  if (pF0 == pFInvalidFace || pF1 == pFInvalidFace)
    return (pCInvalidCell);
  int iNC0 = pF0->getNumCells();
  int iNC1 = pF1->getNumCells();
  for (int i0 = 0; i0 < iNC0; i0++) {
    const Cell *pC0 = pF0->getCell(i0);
    if (!pC0->isValid())
      continue;
    for (int i1 = 0; i1 < iNC1; i1++) {
      const Cell *pC1 = pF1->getCell(i1);
      if (pC0 == pC1)
	return const_cast<Cell*>(pC0);
    }
  }
  return (pCInvalidCell);
}

Vert*
findCommonVert(Face* pF0, Face* pF1)
{
  if (pF0 == pFInvalidFace || pF1 == pFInvalidFace)
    return (pVInvalidVert);
  GRUMMP_Entity *apV0[4], *apV1[4];
  int iNV0 = pF0->getNumVerts(), iNV1 = pF1->getNumVerts();
  pF0->getAllVertHandles(apV0);
  pF1->getAllVertHandles(apV1);
  for (int i = 0; i < iNV0; i++) {
    for (int j = 0; j < iNV1; j++) {
      if (apV0[i] == apV1[j])
	return static_cast<Vert*>(apV0[i]);
    }
  }
  return (pVInvalidVert);
}

Vert *
findCommonVert(const Face* const pF0, const Face* const pF1,
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

void
Face::getAllCellHandles(GRUMMP_Entity* apEnt[]) const
{
  int i = 0;
  if (m_leftCell->isValid() && !m_leftCell->isBdryCell()) {
    apEnt[i++] = m_leftCell;
  }
  if (m_rightCell->isValid() && !m_rightCell->isBdryCell()) {
    apEnt[i] = m_rightCell;
  }
}

void
Face::printFaceInfo() const
{
  char faceloc[7][20] =
    { "Unknown", "BdryFace", "BdryTwoSide", "PseudoSurface", "Interior",
	"Exterior", "InvalidFaceLoc" };

  logMessage(1, "Face with %d verts of type %s\n", getNumVerts(),
	     faceloc[m_faceLoc]);

  for (int i = 0; i < getNumVerts(); i++) {
    m_verts[i]->printVertInfo();
  }

}

