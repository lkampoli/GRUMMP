#include "GR_Face.h"
#include "GR_BFace.h"
#include "GR_Cell.h"
#include "GR_Vec.h"

BFace&
BFace::operator=(const BFace& BF)
{
  assert(getType() == BF.getType());
  assert(!BF.isDeleted());
  if (this != &BF) {
    assert(getNumFaces() == 1 || getNumFaces() == 2);
    assert(BF.getNumFaces() == getNumFaces());
    for (int i = 0; i < getNumFaces(); i++) {
      if (!isDeleted() && m_faces[i]->isValid() && m_faces[i]->hasCell(this))
	m_faces[i]->removeCell(this);
      m_faces[i] = BF.m_faces[i];
      if (m_faces[i]->isValid())
	m_faces[i]->replaceCell(&BF, this);
    }
    copyAllFlags(BF);
    m_offsetInsRequested = BF.m_offsetInsRequested;
    m_neverDoOffsetIns = BF.m_neverDoOffsetIns;
    m_insertionLoc[0] = BF.m_insertionLoc[0];
    m_insertionLoc[1] = BF.m_insertionLoc[1];
    m_insertionLoc[2] = BF.m_insertionLoc[2];
    m_bdryCond = BF.m_bdryCond;
  }
  return (*this);
}

void
BFace::calcVecSize(double adRes[]) const
{
  assert(getFace(0)->isValid());
  int iDim = getFace(0)->getVert(0)->getSpaceDimen();
  if (isDeleted()) {
    adRes[0] = adRes[1] = 0;
    if (iDim == 3)
      adRes[2] = 0;
    return;
  }
  Face *pF = m_faces[0];
  int iSign = (pF->getLeftCell() == dynamic_cast<const Cell*>(this)) ? 1 : -1;
  pF->calcNormal(adRes);
  if (iDim == 2) {
    vSCALE2D(adRes, iSign);
  }
  else {
    vSCALE3D(adRes, iSign);
  }
}

void
BFace::calcUnitNormal(double adRes[]) const
{
  assert(getFace(0)->isValid());
  adRes[0] = adRes[1] = adRes[2] = 0;
  int iDim = getFace(0)->getVert(0)->getSpaceDimen();
  if (isDeleted()) {
    return;
  }
  Face *pF = m_faces[0];
  int iSign = (dynamic_cast<BFace*>(pF->getLeftCell()) == this) ? 1 : -1;
  pF->calcNormal(adRes);

  if (iDim == 2) {
    vSCALE2D(adRes, iSign);
    NORMALIZE2D(adRes);
  }
  else {
    vSCALE3D(adRes, iSign);
    NORMALIZE3D(adRes);
  }
}

const Vert*
BFace::getVert(const int i) const
{
  assert(i >= 0 && i < getNumVerts());
  assert(
      dynamic_cast<BFace*>(getFace(0)->getLeftCell()) == this
	  || dynamic_cast<BFace*>(getFace(0)->getRightCell()) == this);
  bool qReverse = true;
  // Traverse the bdry face in the direction so that the face is
  // traversed in the opposite direction than it is for the interior
  // cell.
  if (getType() == eBdryEdge)
    qReverse = false;
  if (getType() == eIntBdryEdge)
    qReverse = false;
  if (dynamic_cast<BFace*>(getFace(0)->getLeftCell()) == this)
    qReverse = !qReverse;
  if (qReverse)
    return (getFace(0)->getVert(getNumVerts() - i - 1));
  else
    return (getFace(0)->getVert(i));
}

double
BFace::calcSize() const
{
  if (isDeleted())
    return (double(0));
  else
    return m_faces[0]->calcSize();
}

