#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

const Vert*
HexCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert(i >= 0 && i <= 7);
  switch (i)
    {
    case 0:
      return findCommonVert(m_faces[0], m_faces[1], m_faces[4]);
    case 1:
      return findCommonVert(m_faces[0], m_faces[2], m_faces[1]);
    case 2:
      return findCommonVert(m_faces[0], m_faces[3], m_faces[2]);
    case 3:
      return findCommonVert(m_faces[0], m_faces[4], m_faces[3]);
    case 4:
      return findCommonVert(m_faces[5], m_faces[1], m_faces[4]);
    case 5:
      return findCommonVert(m_faces[5], m_faces[2], m_faces[1]);
    case 6:
      return findCommonVert(m_faces[5], m_faces[3], m_faces[2]);
    case 7:
      return findCommonVert(m_faces[5], m_faces[4], m_faces[3]);
    default:
      // Should never get here in debug mode.
      assert(0);
      return pVInvalidVert;
    }
}

void
HexCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  for (int i = 0; i < 8; i++) {
    aHandles[i] = const_cast<Vert*>(getVert(i));
  }
}

double
HexCell::calcSize() const
{
  assert2(0, "Volume calculation not supported for hexes");
  assert(doFullCheck());
//   const Vert *pV0, *pV1, *pV2;
//   const Face* const pFBase = ppFFaces[0];
//   pV0 = pFBase->pVVert(0);
//   if (pFBase->pCCellRight() == (const Cell*)this) {
//     pV1 = pFBase->pVVert(1);
//     pV2 = pFBase->pVVert(2);
//   }
//   else {
//     pV1 = pFBase->pVVert(2);
//     pV2 = pFBase->pVVert(1);
//   }

  return (0);
}

bool
HexCell::isClosed() const
{
  // Are all verts part of the hex?
  for (int iF = getNumFaces() - 1; iF >= 0; iF--) {
    const Face *pF = getFace(iF);
    for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
      if (!hasVert(pF->getVert(iV)))
	return false;
    }
  }

  // Is each vert shared by exactly three faces?  (Should be, based on
  // the previous result...
  static const int hexTuples[][4] =
    {
      { 0, 1, 2, 1 },
      { 0, 2, 3, 2 },
      { 0, 3, 4, 3 },
      { 0, 4, 1, 0 },
      { 5, 1, 2, 5 },
      { 5, 2, 3, 6 },
      { 5, 3, 4, 7 },
      { 5, 4, 1, 4 } };
  for (int ii = 0; ii < 8; ii++) {
    const Face* pFA = getFace(hexTuples[ii][0]);
    const Face* pFB = getFace(hexTuples[ii][1]);
    const Face* pFC = getFace(hexTuples[ii][2]);
    Vert* pV = findCommonVert(pFA, pFB, pFC);
    const Vert* pVExpected = getVert(hexTuples[ii][3]);
    if (!pV->isValid() || pV != pVExpected)
      return false;
  }

  return true;
}
