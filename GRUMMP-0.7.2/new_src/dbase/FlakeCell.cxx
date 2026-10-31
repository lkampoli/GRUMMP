#include <math.h>
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_Vec.h"

const Vert*
FlakeCell::getVert(const int i) const
{
  assert(doFullCheck());
  assert(i >= 0 && i <= 3); // Vertex index out of range
  const Face* const pFQuad = m_faces[0];
  switch (i)
    {
    case 0:
    case 2:
      return pFQuad->getVert(i);
    case 1:
    case 3:
      if (pFQuad->getRightCell() == static_cast<const Cell*>(this))
	return pFQuad->getVert(i);
      else
	return pFQuad->getVert(4 - i);
    }
  return pVInvalidVert;
}

double
FlakeCell::calcSize() const
{
  // Could do this calculation assuming a curved quad and two flat tris, but
  // at this point at least that isn't needed.
  return 0;
}

void
FlakeCell::getAllVertHandles(GRUMMP_Entity* aHandles[]) const
{
  assert(doFullCheck());
  for (int i = 0; i < 4; i++) {
    aHandles[i] = const_cast<Vert*>(getVert(i));
  }
}

bool
FlakeCell::isClosed() const
{
  for (int iF = getNumFaces() - 1; iF >= 0; iF--) {
    const Face *pF = getFace(iF);
    for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
      if (!hasVert(pF->getVert(iV)))
	return false;
    }
  }
  return true;
}

