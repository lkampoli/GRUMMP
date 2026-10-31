#include "GR_assert.h"
#include "GR_CellSkel.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"

int
CellSkel::doFullCheck() const
{
  if (!isValid() || isDeleted() ) {
    return 0;
  }
  bool qOkay = true;
  for (int i = 0; i < getNumVerts() && qOkay; i++) {
    qOkay = getVert(i)->isValid();
  }
  return (qOkay);
}

double
CellSkel::calcMaxDihed() const
{
  assert(isValid());
  int iNDihed;
  double adDihed[12];
  calcAllDihed(adDihed, &iNDihed);
  double dRetVal = 0;
  for (int i = 0; i < iNDihed; i++)
    dRetVal = (dRetVal > adDihed[i]) ? dRetVal : adDihed[i];
  return dRetVal;
}

double
CellSkel::calcMinDihed() const
{
  assert(isValid());
  int iNDihed;
  double adDihed[12];
  calcAllDihed(adDihed, &iNDihed);
  double dRetVal = 180;
  for (int i = 0; i < iNDihed; i++)
    dRetVal = (dRetVal < adDihed[i]) ? dRetVal : adDihed[i];
  return dRetVal;
}

double
CellSkel::calcMaxSolid() const
{
  assert(isValid());
  int iNSolid;
  double adSolid[12];
  calcAllSolid(adSolid, &iNSolid);
  double dRetVal = 0;
  for (int i = 0; i < iNSolid; i++)
    dRetVal = (dRetVal > adSolid[i]) ? dRetVal : adSolid[i];
  return dRetVal;
}

double
CellSkel::calcMinSolid() const
{
  assert(isValid());
  int iNSolid;
  double adSolid[12];
  calcAllSolid(adSolid, &iNSolid);
  double dRetVal = 360;
  for (int i = 0; i < iNSolid; i++)
    dRetVal = (dRetVal < adSolid[i]) ? dRetVal : adSolid[i];
  return dRetVal;
}

void
CellSkel::printCellInfo() const
{

  char cellType[20][20] =
    { "TriCell", "QuadCell", "Tet", "Pyr", "Prism", "Hex", "BdryEdge",
	"IntBdryEdge", "TriBFace", "IntTriBFace", "QuadBFace", "IntQuadBFace",
	"TriCV", "QuadCV", "TetCV", "PyrCV", "PrismCV", "HexCV", "Invalid" };

  logMessage(1, "Cell %p with %d verts of type %s\n", this, getNumVerts(),
	     cellType[m_cellType]);

  for (int i = 0; i < getNumVerts(); i++)
    getVert(i)->printVertInfo();
}
