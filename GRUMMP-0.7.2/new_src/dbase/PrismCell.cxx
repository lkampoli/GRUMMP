#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Geometry.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

void
evalAllDihedsFromCoords(int& nDiheds, double *result, const double ptA[],
			const double ptB[],
			const double ptC[], const double ptD[],
			const double ptE[], const double ptF[])
{
  // ABC are cyclic at the bottom (RH inward)
  // DEF are cyclic at the top (RH outward)
  // Quads on the sides are (RH outward)
  // ABDE, BCEF CAFD

  double normABC[3], normDEF[3], normABED[3], normBCFE[3], normCADF[3];
  calcQuadNormal(ptA, ptB, ptE, ptD, normABED);
  calcQuadNormal(ptB, ptC, ptF, ptE, normBCFE);
  calcQuadNormal(ptC, ptA, ptD, ptF, normCADF);

  calcNormal3D(ptA, ptC, ptB, normABC);
  calcNormal3D(ptD, ptE, ptF, normDEF);

  NORMALIZE3D(normABC);
  NORMALIZE3D(normDEF);
  NORMALIZE3D(normABED);
  NORMALIZE3D(normBCFE);
  NORMALIZE3D(normCADF);

  double dots[9];
  // First the bottom edges (ABC)
  dots[0] = -dDOT3D(normABC, normABED);  // AB
  dots[1] = -dDOT3D(normABC, normBCFE);  // BC
  dots[2] = -dDOT3D(normABC, normCADF);  // CA
  // Now the top edges (DEF)
  dots[3] = -dDOT3D(normDEF, normABED);  // DE
  dots[4] = -dDOT3D(normDEF, normBCFE);  // EF
  dots[5] = -dDOT3D(normDEF, normCADF);  // FD
  // Now the lateral edges
  dots[6] = -dDOT3D(normABED, normCADF); // AD
  dots[7] = -dDOT3D(normBCFE, normABED); // BE
  dots[8] = -dDOT3D(normCADF, normBCFE); // CF

  // Finally, get the dihedral angles from these.
  nDiheds = 9;
  for (int ii = 0; ii < 9; ii++) {
    result[ii] = GR_acos(dots[ii]);
  }
}

const Vert*
PrismCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert(i >= 0 && i <= 5); // Vertex index out of range
  assert(doFullCheck());
  assert(i >= 0 && i <= 5); // Vertex index out of range
  switch (i)
    {
    case 0:
      return findCommonVert(m_faces[0], m_faces[2], m_faces[3]);
    case 1:
      return findCommonVert(m_faces[0], m_faces[1], m_faces[3]);
    case 2:
      return findCommonVert(m_faces[1], m_faces[2], m_faces[3]);
    case 3:
      return findCommonVert(m_faces[0], m_faces[2], m_faces[4]);
    case 4:
      return findCommonVert(m_faces[0], m_faces[1], m_faces[4]);
    case 5:
      return findCommonVert(m_faces[1], m_faces[2], m_faces[4]);
    default:
      assert(0);
      return pVInvalidVert;
    }
}

void
PrismCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  for (int i = 0; i < 6; i++) {
    aHandles[i] = const_cast<Vert*>(getVert(i));
  }
}

double
PrismCell::calcSize() const
{
  assert2(0, "No volume calculation yet for prisms");
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
PrismCell::isClosed() const
{
  for (int iF = getNumFaces() - 1; iF >= 0; iF--) {
    const Face *pF = getFace(iF);
    for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
      if (!hasVert(pF->getVert(iV)))
	return false;
    }
  }

  // Is each vert shared by exactly three faces?  (Should be, based on
  // the previous result...
  static const int prismTuples[][4] =
    {
      { 3, 2, 0, 0 },
      { 3, 0, 1, 1 },
      { 3, 1, 2, 2 },
      { 4, 2, 0, 3 },
      { 4, 0, 1, 4 },
      { 4, 1, 2, 5 } };
  for (int ii = 0; ii < 6; ii++) {
    const Face* pFA = getFace(prismTuples[ii][0]);
    const Face* pFB = getFace(prismTuples[ii][1]);
    const Face* pFC = getFace(prismTuples[ii][2]);
    Vert* pV = findCommonVert(pFA, pFB, pFC);
    const Vert* pVExpected = getVert(prismTuples[ii][3]);
    if (!pV->isValid() || pV != pVExpected)
      return false;
  }

  return true;
}

