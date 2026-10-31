#include "GR_config.h"
#include "GR_assert.h"
#include "GR_CellCV.h"

double
CellCV::calcSize() const
{
  // Different stuff for different types of CellCV!
  int iNF = getNumFaces();
  switch (iNV)
    {
    case 3:
      // Triangle.
      {
	assert(iNF == 3);
	const Vert *pV0 = getVert(0);
	const Vert *pV1 = getVert(1);
	const Vert *pV2 = getVert(2);
	int iDim = pV0->getSpaceDimen();
	double adEdge01[] =
	  { 0, 0, 0 };
	double adEdge12[] =
	  { 0, 0, 0 };
	for (int i = 0; i < iDim; i++) {
	  adEdge01[i] = pV1->getCoords()[i] - pV0->getCoords()[i];
	  adEdge12[i] = pV2->getCoords()[i] - pV1->getCoords()[i];
	}
	if (iDim == 2) {
	  return (0.5 * fabs(dCROSS2D(adEdge01, adEdge12)));
	}
	else {
	  double adRes[3];
	  vCROSS3D(adEdge01, adEdge12, adRes);
	  return 0.5 * dMAG3D(adRes);
	} // Triangle in 2D or 3D
      }
    case 4:
      assert(iNF == 4);
      // Tet or quad.  No way to tell, for quads in a surface mesh!
      // Assume it's a tet if it's 3D.
      if (getVert(0)->getSpaceDimen() == 2) {
	assert(0);
	return -1;
      }
      else {
	const Vert *pV0 = getVert(0);
	const Vert *pV1 = getVert(1);
	const Vert *pV2 = getVert(2);
	const Vert *pV3 = getVert(3);

	double adVec1[] =
	  { pV1->x() - pV0->x(), pV1->y() - pV0->y(), pV1->z() - pV0->z() };
	double adVec2[] =
	  { pV2->x() - pV0->x(), pV2->y() - pV0->y(), pV2->z() - pV0->z() };
	double adVec3[] =
	  { pV3->x() - pV0->x(), pV3->y() - pV0->y(), pV3->z() - pV0->z() };
	double adCross[3];
	vCROSS3D(adVec2, adVec3, adCross);
	return dDOT3D(adCross, adVec1) / 6.;
      }
    case 5:
      // Pyramid
      assert(iNF == 5);
      assert(0);
      return -1;
    case 6:
      // Prism
      assert(iNF == 5);
      assert(0);
      return -1;
    case 8:
      // Hex
      assert(iNF == 6);
      assert(0);
      return -1;
    default:
      assert(0);
      return -1;
    }
  // Should never reach here, but some compilers whine anyway.
//   assert(0);
//   return -1;
}

void
CellCV::calcAllDihed(double adDihed[], int* const piNDihed,
		     const bool in_degrees) const
{
  // Different stuff for different types of CellCV!
  int iNF = getNumFaces();
  if (iNV == 3) {
    assert(iNF == 3);
    // Triangle.
    const Vert *pV0 = getVert(0);
    const Vert *pV1 = getVert(1);
    const Vert *pV2 = getVert(2);
    int iDim = pV0->getSpaceDimen();
    double adEdge01[] =
      { 0, 0, 0 };
    double adEdge12[] =
      { 0, 0, 0 };
    double adEdge20[] =
      { 0, 0, 0 };
    for (int i = 0; i < iDim; i++) {
      adEdge01[i] = pV1->getCoords()[i] - pV0->getCoords()[i];
      adEdge12[i] = pV2->getCoords()[i] - pV1->getCoords()[i];
      adEdge20[i] = pV0->getCoords()[i] - pV2->getCoords()[i];
    }
    NORMALIZE3D(adEdge01);
    NORMALIZE3D(adEdge12);
    NORMALIZE3D(adEdge20);
    double dArg0 = -dDOT3D(adEdge20, adEdge01);
    double dArg1 = -dDOT3D(adEdge01, adEdge12);
    double dArg2 = -dDOT3D(adEdge12, adEdge20);

    double dScale = in_degrees ? 180. / M_PI : 1;
    adDihed[0] = GR_acos(dArg0) * dScale; // angle @ 0
    adDihed[1] = GR_acos(dArg1) * dScale; // angle @ 1
    adDihed[2] = GR_acos(dArg2) * dScale; // angle @ 2
    *piNDihed = 3;
  } // Triangle in 2D or 3D
  else if (iNV == 4) {
    // Tet or quad.  No way to tell, for quads in a surface mesh!  So
    // assume tets in 3D.
    if (getVert(0)->getSpaceDimen() == 2) {
      // Quad in 2D
    }
    else {
      // Tet in 3D
      const TetCellCV *pTCCV = reinterpret_cast<const TetCellCV*>(this);
      //      (pVVert(0), pVVert(1), pVVert(2), pVVert(3));
      pTCCV->TetCellCV::calcAllDihed(adDihed, piNDihed);
      assert(*piNDihed == 6);
    }
  }
  else if (iNV == 5) {
    // Pyramid
    assert(0);
  }
  else if (iNV == 6) {
    // Prism
    assert(0);
  }
  else if (iNV == 8) {
    // Hex
    assert(0);
  }
}

void
CellCV::calcAllSolid(double adDihed[], int* const piNDihed,
		     const bool in_degrees) const
{
  // Different stuff for different types of CellCV, and we have to
  // deduce which type we have!
  int iNF = getNumFaces();
  if (iNV == 3) {
    assert(iNF == 3);
    // Triangle.
    calcAllDihed(adDihed, piNDihed, in_degrees);
  } // Triangle in 2D or 3D
  else if (iNV == 4) {
    // Tet or quad.  No way to tell, for quads in a surface mesh!
    assert(0);
  }
  else if (iNV == 5) {
    // Pyramid
    assert(0);
  }
  else if (iNV == 6) {
    // Prism
    assert(0);
  }
  else if (iNV == 8) {
    // Hex
    assert(0);
  }
}
