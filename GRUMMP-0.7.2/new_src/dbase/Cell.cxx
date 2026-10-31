#include "GR_assert.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"
#include "CubitBox.hpp"

Cell::Cell(const bool qIsInternalBdry, Face* const pFLeft, Face* const pFRight) :
    CellSkel(), m_faces(NULL), m_anisoCircRad(-1), bCombineToQuad(false)
// A constructor for bdry faces of all stripes.
{
  if (qIsInternalBdry) {
    m_faces = new Face*[2];
    m_faces[0] = pFLeft;
    m_faces[1] = pFRight;
  }
  else {
    m_faces = new Face*[1];
    m_faces[0] = pFLeft;
    assert(pFRight == pFInvalidFace);
  }
  // Set flags to default values.
  // Cell type will be explicitly set by derived class constructor.
}

Cell::Cell(const int iNFaces, const int iNVerts, Face** const ppFIn) :
    CellSkel(), m_faces(NULL), m_anisoCircRad(-1)
{
  assert(iNFaces >= 2 && iNFaces <= 6);
#ifndef NDEBUG
  switch (iNFaces)
    {
    case 2:
      assert(iNVerts == 2);
      break;
    case 3:
      assert(iNVerts == 3 || iNVerts == 4);  // Tri or flake
      break;
    case 4:
      assert(iNVerts == 4);
      break;
    case 5:
      assert(iNVerts == 5 || iNVerts == 6);
      break;
    case 6:
      assert(iNVerts == 8);
      break;
    default:
      assert2(0, "Illegal number of faces for cell");
      break;
    }
#endif
  // Cell type will be explicitly set by derived class constructor.
  m_faces = new Face*[iNFaces];
  int i;
  if (ppFIn) {
    for (i = 0; i < iNFaces; i++) {
      m_faces[i] = ppFIn[i];
      assert2(m_faces[i] != pFInvalidFace, "Uninitialized face");
    }
    canonicalizeFaceOrder();
  }
  else
    for (i = 0; i < iNFaces; i++)
      m_faces[i] = pFInvalidFace;
  setRegion(1);
}

Cell::~Cell()
{
  if (m_faces != NULL)
    delete[] m_faces;
}

Cell&
Cell::operator=(const Cell& C)
{
  assert(getType() == C.getType());
  if (this != &C) {
    int i;
    // The cell *this will no longer be connected to its previous faces,
    // if the cell still exists.  If it doesn't, then either the face is
    // deleted or the face has been modified so that it isn't attached
    // to the cell anymore.
    if (!isDeleted()) {
      for (i = getNumFaces() - 1; i >= 0; i--)
	if (m_faces[i]->isValid() && m_faces[i]->hasCell(this))
	  m_faces[i]->removeCell(this);
    }
    // Copy the face pointers and connect them to this cell
    for (i = getNumFaces() - 1; i >= 0; i--) {
      m_faces[i] = C.m_faces[i];
      if (m_faces[i]->isValid()) {
	m_faces[i]->replaceCell(&C, this);
	assert(m_faces[i]->hasCell(this));
      }
    }
    copyAllFlags(C);
  }
  return (*this);
}

int
Cell::doFullCheck() const
{
  if (!isValid() || isDeleted())
    return 0;
  bool qOkay = true;

  int iInvalid = 0;
//#pragma omp critical (doFullCheckCell)
//  {
  int nFaces = getNumFaces(); // Speeds up debug-mode code a surprising amount
  for (int i = 0; qOkay && i < nFaces; i++) {
//		const Face* pFThisFace = getFace(i);
    // All of the bounds and validity checks are already done.
    // Again, this speeds up debug-mode code a surprising amount in some cases.
    const Face* pFThisFace = m_faces[i];

    if (pFThisFace->isValid()) {
      qOkay = qOkay && pFThisFace->hasCell(this);
    }

    else
      iInvalid++;
  }
//  }

  if (iInvalid == 0 && qOkay)
    return 1;
  else if (iInvalid == getNumFaces() && qOkay)
    return 2;
  else
    return 0;
}

const Face*
Cell::getFace(const int i) const
{
  assert(isValid());
  assert2(i >= 0 && i < getNumFaces(), "Index out of range");
  return m_faces[i];
}

Face*
Cell::getFace(const int i)
{
  assert(isValid());
  assert2(i >= 0 && i < getNumFaces(), "Index out of range");
  return m_faces[i];
}

bool
Cell::hasFace(const Face* pF) const
{
  assert(isValid());
  assert(pF->isValid());
  for (int i = 0; i < getNumFaces(); i++)
    if (m_faces[i] == pF)
      return (true);
  return (false);
}

const Face*
Cell::getOppositeFace(const Vert* const pV) const
{
  assert(0 && pV->isValid());
  return (pFInvalidFace);
}
Face*
Cell::getOppositeFace(const Vert* const pV)
{
  assert(0 && pV->isValid());
  return (pFInvalidFace);
}

Vert*
Cell::getOppositeVert(const Face* const pF) const
{
  assert(0 && pF->isValid());
  return (pVInvalidVert);
}

void
Cell::addFace(Face* const pFNew, const int iInd)
{
  assert(isValid());
  assert2(pFNew->isValid(), "Tried to add non-existent face");
  if (iInd == -1) { // No face index specified
    int iNull = -1;
    for (int i = 0; i < getNumFaces() && iNull == -1; i++) {
      assert2(pFNew != m_faces[i], "Face already exists for this cell");
      if (!m_faces[i]->isValid()) {
	iNull = i;
      }
    }
    assert2(iNull != -1, "No blank faces for this cell.");
    m_faces[iNull] = pFNew;
  }
  else {
    assert(iInd >= 0 && iInd < getNumFaces());
    assert2(m_faces[iInd] == pFInvalidFace,
	    "Face with that index exists already");
    m_faces[iInd] = pFNew;
  }
}

void
Cell::replaceFace(const Face* const pFOld, Face* const pFNew)
{
  assert(isValid());
  assert2(pFOld->isValid(), "Tried to replace non-existent face");
  assert2(pFNew->isValid(), "Tried to add non-existent face");
  int iRepl = -1;
  for (int i = 0; i < getNumFaces(); i++) {
    if (m_faces[i] == pFOld)
      iRepl = i;
  }
  assert2(iRepl >= 0, "Tried to replace face which doesn't bound this cell");
  m_faces[iRepl] = pFNew;
}

void
Cell::removeFace(const Face* const pFOld)
{
  if (!isValid())
    return; // A useful short-circuit when destroying part
  // of a mesh.  This way there is no need to be
  // careful about the order in which objects are
  // killed.
  assert2(pFOld->isValid(), "Tried to remove non-existent face");
  int iRepl = -1;
  for (int i = 0; i < getNumFaces(); i++) {
    if (m_faces[i] == pFOld)
      iRepl = i;
  }
  assert2(iRepl >= 0, "Tried to remove face which doesn't bound this cell");
  m_faces[iRepl] = pFInvalidFace;
}

void
Cell::assignFaces(Face* apF[])
{
  assert(isValid());
  assert(apF != NULL);
  for (int i = 0; i < getNumFaces(); i++) {
    m_faces[i] = apF[i];
    assert2(m_faces[i] != pFInvalidFace, "Uninitialized face");
  }
  canonicalizeFaceOrder();
}

void
Cell::assignFaces(Face* pF0, Face* pF1, Face* pF2, Face* pF3, Face* pF4,
		  Face* pF5)
// All alternate interface to vAssign, so that stuff doesn't have
// to be loaded into an array by the user.
{
//#pragma omp critical(CellassignFaces)
//{
  assert(isValid());
  Face *apF[6];
  apF[0] = pF0;
  apF[1] = pF1;
  apF[2] = pF2;
  apF[3] = pF3;
  apF[4] = pF4;
  apF[5] = pF5;
  assignFaces(apF);
//}
}

void
Cell::resetAllData()
{
//#pragma omp critical(resetAllDataCell)
//{
  setDefaultFlags();
  for (int i = 0; i < getNumFaces(); i++)
    m_faces[i] = pFInvalidFace;
//}
}

const Face*
SimplexCell::getOppositeFace(const Vert* const pV) const
{
  assert(doFullCheck());
  assert(pV->isValid());

  for (int i = 0; i < getNumFaces(); i++) {
    assert(getFace(i)->isValid());
    if (!getFace(i)->hasVert(pV))
      return (getFace(i));
  }
  assert2(0, "Strange. All faces of a simplex are incident on a single vert.");
  return (pFInvalidFace);
}

Face*
SimplexCell::getOppositeFace(const Vert* const pV)
{
  return const_cast<Face*>(static_cast<const SimplexCell*>(this)->getOppositeFace(
      pV));
}

Vert*
SimplexCell::getOppositeVert(const Face* const pFFaceIn) const
{
  assert(doFullCheck());
  assert(pFFaceIn->isValid());

  int iNVerts = getNumVerts();
  Face* pF = (m_faces[0] == pFFaceIn) ? m_faces[1] : m_faces[0];
  GRUMMP_Entity* apV[3];
  pFFaceIn->getAllVertHandles(apV);
  for (int i = 0; i < iNVerts - 1; i++) {
    Vert *pVCand = pF->getVert(i);
    bool qFoundMatch = false;
    for (int j = 0; j < iNVerts - 1 && !qFoundMatch; j++) {
      qFoundMatch = (apV[j] == pVCand);
    }
    if (!qFoundMatch)
      return pVCand;
  }
  assert2(0, "Strange. All faces of a simplex are incident on a single vert.");
  return (pVInvalidVert);
}

bool
Cell::hasVert(const Vert* pV) const
{
  assert(isValid());
  for (int i = 0; i < getNumFaces(); i++)
    if (m_faces[i]->hasVert(pV))
      return true;
  return false;
}

void
Cell::calcCentroid(double adCent[]) const
{
  // This is a really crappy centroid calculation.  More precisely, it's
  // correct for tris and tets, and bogus for other cell types.  So this
  // should probably be moved to SimplexCell::vCentroid.  Even then,
  // it'll be obsolete for general mesh geometry...
  assert(isValid());
  // assert(iFullCheck() == 1);
  int iDim = getVert(0)->getSpaceDimen();
  switch (iDim)
    {
    case 2:
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

bool
Cell::isCircumcircleEmpty() const
{
  int iF;
  for (iF = 0; iF < this->getNumFaces(); iF++) {
    if (this->getFace(iF)->isBadFace())
      return false;
  }
  return true;

}

CubitBox
Cell::bounding_box() const
{

  CubitVector mini( LARGE_DBL, LARGE_DBL, LARGE_DBL), maxi(-LARGE_DBL,
							   -LARGE_DBL,
							   -LARGE_DBL);

  assert(getNumVerts() > 0);

  int space_dimen = getVert(0)->getSpaceDimen();
  assert(space_dimen == 2 || space_dimen == 3);

  double x, y, z;

  for (int i = 0; i < getNumVerts(); i++) {

    const Vert* vertex = getVert(i);
    assert(vertex->getSpaceDimen() == space_dimen);

    x = vertex->x();
    y = vertex->y();
    z = space_dimen == 3 ? vertex->z() : 0.;

    mini.set(std::min(mini.x(), x), std::min(mini.y(), y),
	     std::min(mini.z(), z));
    maxi.set(std::max(maxi.x(), x), std::max(maxi.y(), y),
	     std::max(maxi.z(), z));

  }

  return CubitBox(mini, maxi);

}

Face *
findCommonFace(Cell* const pC0, Cell* const pC1)
{
  if (pC0->doFullCheck() != 1 || pC1->doFullCheck() != 1)
    return (pFInvalidFace);
  for (int iF0 = 0; iF0 < pC0->getNumFaces(); iF0++) {
    Face *pF0 = pC0->getFace(iF0);
    for (int iF1 = 0; iF1 < pC1->getNumFaces(); iF1++)
      if (pF0 == pC1->getFace(iF1))
	return pF0;
  }
  return pFInvalidFace;
}

void
Cell::getAllFaceHandles(GRUMMP_Entity* aHandles[]) const
{
  for (int i = getNumFaces() - 1; i >= 0; i--) {
    aHandles[i] = m_faces[i];
  }
}
int
Cell::getIndex()
{
  return Cell::indexx;
}

void
Cell::setIndex()
{
  if (this->getFace(0)->isBdryFace() or this->getFace(1)->isBdryFace()
      or this->getFace(2)->isBdryFace()) {
    Cell::indexx = 0;
  }
  else {
    Cell* Cell1 = (this->getFace(0)->getOppositeCell(this));
    assert(Cell1);
    int index1 = Cell1->getIndex();
    Cell* Cell2 = (this->getFace(1)->getOppositeCell(this));
    assert(Cell2);
    int index2 = Cell2->getIndex();
    Cell* Cell3 = (this->getFace(2)->getOppositeCell(this));
    assert(Cell3);
    int index3 = Cell3->getIndex();
    int mindex = index1;
    mindex = mindex < index2 ? mindex : index2;
    mindex = mindex < index3 ? mindex : index3;
    Cell::indexx = mindex + 1;
  }
}

