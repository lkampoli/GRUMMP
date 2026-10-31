#include "GR_config.h"
#include "GR_assert.h"
#include "GR_VertConnect.h"

int
VertConnect::iCountKeptNeigh() const
{
  int ii, mySize = int(iSize());
  int iRetVal = 0;
  for (ii = 0; ii < mySize; ii++) {
    if (!pVNeigh(ii)->isDeletionRequested())
      iRetVal++;
  }
  return iRetVal;
}

int
iCountKept2nd(const int iV, const VertConnect aVV[], const Mesh* const pMesh)
{
  std::set<Vert*> spVSec;
  size_t iN;
  for (iN = 0; iN < aVV[iV].iSize(); iN++) {
    int iNeigh = pMesh->getVertIndex(aVV[iV].pVNeigh(iN));
    for (size_t iN2 = 0; iN2 < aVV[iNeigh].iSize(); iN2++)
      spVSec.insert(aVV[iNeigh].pVNeigh(iN2));
  }
  for (iN = 0; iN < aVV[iV].iSize(); iN++)
    spVSec.erase(aVV[iV].pVNeigh(iN));
  spVSec.erase(pMesh->getVert(iV));

  int iRetVal = 0;
  for (std::set<Vert*>::iterator iterV = spVSec.begin(); iterV != spVSec.end();
      ++iterV) {
    Vert *pV = *iterV;
    if (!pV->isDeletionRequested())
      iRetVal++;
  }
  return iRetVal;
}

double
VertConnect::dClosestLegitSurfPoint(const Vert *pV,
				    const Vert::VertType VTLegitVertType) const
				    // Find the nearest point along a common line/surface where we want to
				    // start marching into a higher-dimensional (surface/volume)
				    // pseudo-structured anisotropic mesh.
{
  size_t iV;
  ssize_t iNearest = -1;
  double dMinLen = 1.e20;
  for (iV = 0; iV < iSize(); iV++) {
    Vert *pVN = pVNeigh(iV);
    // This conditional makes it so that only BdryCurve and BdryApex
    // points are checked to find the nearest point for BdryCurve
    // points, and so on.
    if (pVN->getVertType() <= VTLegitVertType) {
      double dDist = calcDistanceBetween(pV, pVN);
      dMinLen = min(dMinLen, dDist);
      iNearest = iV;
    }
  }
  if (iNearest == -1)
    return 0;
  else
    return dMinLen;
}

bool
VertConnect::qNormalNeighborExists(const Vert * const pVCurr,
				   const Vert * const pVPrev,
				   Vert * &pVNext) const
{
  // Allow about a 20 deg (acos(dDotTol)) variation in angle from one step to
  // the next.
  static const double dDotTol = 0.94;
  double dMinDistInMarchingDir = LARGE_DBL;
  Vert *pVBest = pVInvalidVert;

  // This works even for 2D, because the z-coordinate is always zero.
  double adVecPrevCurr[] = adDIFF3D(pVCurr->getCoords(),
      pVPrev->getCoords());

  NORMALIZE3D(adVecPrevCurr);

  for (GR_index_t iV = 0; iV < iSize(); iV++) {
    Vert *pVCand = pVNeigh(iV);
    if (pVCand == pVPrev)
      continue;
    double adVecCurrCand[] = adDIFF3D(pVCand->getCoords(),
	pVCurr->getCoords());
    double dDist = dMAG3D(adVecCurrCand);
    NORMALIZE3D(adVecCurrCand);

    double dDot = dDOT3D(adVecCurrCand, adVecPrevCurr);
    if (dDot > dDotTol && dDist < dMinDistInMarchingDir) {
      dMinDistInMarchingDir = dDist;
      pVBest = pVCand;
    }
  } // Done looping over all neighbors
  pVNext = pVBest;
  return pVNext->isValid();
}
