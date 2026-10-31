#include "GR_Face.h"

void
MultiEdge::vAssign(Vert* const pV0, Vert* const pV1, Cell **apCExtraIn,
		   int iNCellsIn)
{
  m_verts[0] = pV0;
  m_verts[1] = pV1;
  m_leftCell = apCExtraIn[0];
  m_rightCell = apCExtraIn[1];
  iNCells = iNCellsIn;
  assert(iNCells > 2);
  for (int i = 2; i < iNCellsIn; i++)
    apCExtra[i - 2] = apCExtraIn[i];
}

bool
MultiEdge::hasCell(const Cell* const pC) const
{
  bool qRetVal = (pC == m_leftCell || pC == m_rightCell);
  for (int i = 2; i < iNCells && !qRetVal; i++)
    qRetVal = qRetVal || (apCExtra[i - 2] == pC);
  return qRetVal;
}

Cell*
MultiEdge::getOppositeCell(const Cell* const pC) const
// This operator traverses all the cells incident on the face if
// called often enough.
{
  if (pC == m_leftCell)
    return m_rightCell;
  if (pC == m_rightCell)
    return apCExtra[0];
  for (int i = 2; i < iNCells - 1; i++)
    if (pC == apCExtra[i - 2])
      return apCExtra[i - 1];
  assert(pC == apCExtra[iNCells - 3]);
  return m_leftCell;
}

void
MultiEdge::removeCell(const Cell* const pCOld)
{
  // The assumption here is that a cell will eventually replace this
  // one.  That is, this call assumes that the actual non-manifold mesh
  // topology isn't being changed, just the identities of triangles on
  // the surface.
  assert(iNCells > 2 || apCExtra[0] == pCInvalidCell);
  if (pCOld == m_leftCell)
    m_leftCell = apCExtra[iNCells - 3];
  else if (pCOld == m_rightCell)
    m_rightCell = apCExtra[iNCells - 3];
  else {
    int i;
    for (i = 2; i < iNCells; i++) {
      if (pCOld == apCExtra[i - 2]) {
	apCExtra[i - 2] = apCExtra[iNCells - 3];
	break;
      }
    }
    assert(i < iNCells);
  }
  apCExtra[iNCells - 3] = pCInvalidCell;
  iNCells--;
}

void
MultiEdge::replaceCell(const Cell* const pCOld, Cell* const pCNew)
{
  assert(pCNew->isValid());
  if (pCOld == m_leftCell)
    m_leftCell = pCNew;
  else if (pCOld == m_rightCell)
    m_rightCell = pCNew;
  else {
    int i;
    for (i = 2; i < iNCells; i++) {
      if (pCOld == apCExtra[i - 2]) {
	apCExtra[i - 2] = pCNew;
	break;
      }
    }
    assert(i < iNCells);
  }
}

void
MultiEdge::addCell(Cell* const pCNew)
{
  assert(pCNew != NULL);
  switch (iNCells)
    {
    case 0:
      m_leftCell = pCNew;
      break;
    case 1:
      m_rightCell = pCNew;
      break;
    default:
      {
	int ii;
	for (ii = 2; ii < iNCells; ii++) {
	  if (apCExtra[ii - 2] == NULL)
	    continue;
	}
	apCExtra[ii - 2] = pCNew;
      }
      break;
    }
  iNCells++;
  assert(iNCells < 10);
}

