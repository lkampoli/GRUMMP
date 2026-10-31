#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

const Vert*
PyrCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert(i >= 0 && i <= 4); // Vertex index out of range
  const Face* const pFBase = m_faces[4];
  switch (i)
    {
    case 0:
    case 2:
      return pFBase->getVert(i);
    case 1:
    case 3:
      if (pFBase->getRightCell() == static_cast<const Cell*>(this))
	return pFBase->getVert(i);
      else
	return pFBase->getVert(4 - i);
    case 4:
      {
	const Face *pF0 = getFace(0);
	for (int ii = 0; ii < 3; ii++) {
	  const Vert *pV = pF0->getVert(ii);
	  if (!pFBase->hasVert(pV))
	    return pV;
	}
      }
      assert(0); // Should never get here.
      break;
    }
  return pVInvalidVert;
}

double
PyrCell::calcSize() const
{
  // This calculation is exact for a pyramid with a bi-linear base.
  // Straightforward in principle: a coordinate transformation from a canonical
  // pyramid to the real one, then to integrate 1 over the physical pyramid,
  // you integrate det(jacobian) over the canonical pyramid.  In the end,
  // it works out to be exactly the triple product of three vectors.

  const double* coords0 = getVert(0)->getCoords();
  const double* coords1 = getVert(1)->getCoords();
  const double* coords2 = getVert(2)->getCoords();
  const double* coords3 = getVert(3)->getCoords();
  const double* coords4 = getVert(4)->getCoords();

  double vecA[] =
    { (coords0[0] + coords1[0] + coords2[0] + coords3[0]) * 0.25, (coords0[1]
	+ coords1[1] + coords2[1] + coords3[1]) * 0.25, (coords0[2] + coords1[2]
	+ coords2[2] + coords3[2]) * 0.25 };
  double vecB[] =
    { (coords0[0] - coords1[0] - coords2[0] + coords3[0]) * 0.25, (coords0[1]
	- coords1[1] - coords2[1] + coords3[1]) * 0.25, (coords0[2] - coords1[2]
	- coords2[2] + coords3[2]) * 0.25 };
  double vecC[] =
    { (coords0[0] + coords1[0] - coords2[0] - coords3[0]) * 0.25, (coords0[1]
	+ coords1[1] - coords2[1] - coords3[1]) * 0.25, (coords0[2] + coords1[2]
	- coords2[2] - coords3[2]) * 0.25 };
  double vecE[] =
    { coords4[0] - vecA[0], coords4[1] - vecA[1], coords4[2] - vecA[2] };

  double result = (+vecB[0] * vecC[1] * vecE[2] + vecB[1] * vecC[2] * vecE[0]
      + vecB[2] * vecC[0] * vecE[1] - vecB[0] * vecC[2] * vecE[1]
      - vecB[1] * vecC[0] * vecE[2] - vecB[2] * vecC[1] * vecE[0]) * 4. / 3.;

  return result;
}

void
PyrCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  for (int i = 0; i < 5; i++) {
    aHandles[i] = const_cast<Vert*>(getVert(i));
  }
}

bool
PyrCell::isClosed() const
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
  static const int pyrTuples[][4] =
    {
      { 4, 3, 0, 0 },
      { 4, 0, 1, 1 },
      { 4, 1, 2, 2 },
      { 4, 2, 3, 3 },
      { 0, 1, 2, 4 } };
  for (int ii = 0; ii < 5; ii++) {
    const Face* pFA = getFace(pyrTuples[ii][0]);
    const Face* pFB = getFace(pyrTuples[ii][1]);
    const Face* pFC = getFace(pyrTuples[ii][2]);
    Vert* pV = findCommonVert(pFA, pFB, pFC);
    const Vert* pVExpected = getVert(pyrTuples[ii][3]);
    if (!pV->isValid() || pV != pVExpected)
      return false;
  }

  return true;
}

