#include <assert.h>
#include <stdlib.h>
#include "GR_Geometry.h"
#include "GR_Mesh2D.h"
#include "GR_BFace.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Vec.h"

#define iMaxVertNeigh MAX_NUM_PTS


//@ Optimization driver for face swapping.  Swap measure is set elsewhere.
int
Mesh2D::iSwap_deprecated(const int iMaxPasses, const bool /* qAlreadyMarked */)
{
  int iSwaps = 1;
  assert(m_swapMeasure != eNoSwap);

  int iTotal = 0;
#ifndef NDEBUG
  int iPasses = 0;

  while (iSwaps && iPasses < iMaxPasses) {
    ++iPasses;
    if (iMaxPasses > 1)
      logMessage(2, "   Pass %d:", iPasses);
#endif
    iSwaps = 0;
    for (GR_index_t i = 0; i < getNumFaces(); i++) {
      // Always swap recursively
      Face *pF = getFace(i);
      iSwaps += iFaceSwap_deprecated(pF);
    }
    iTotal += iSwaps;
    logMessage(2, "%6d swaps", iSwaps);
#ifndef NDEBUG
    if (iMaxPasses > 1)
      logMessage(2, ",%6d total", iTotal);
#endif
    logMessage(2, "\n");
#ifndef NDEBUG
  }
  if (iMaxPasses > 1)
    logMessage(2, "Face swapping complete.\n");
#endif

  // No need to purge faces and cells after 2D swapping
  setAllHintFaces();
  setVertFaceNeighbors();
  return iTotal;
}

bool
Mesh2D::makeDelaunay()
{
  // This assertion shouldn't matter; face swapping should do nothing for faces
  // adjacent to a quad.
  // assert(isSimplicial());
  // Do incircle face swapping; this gives a Delaunay mesh without
  // further ado.
  GRUMMP::SwapDecider2D *pSD2D = new GRUMMP::DelaunaySwapDecider2D();
  GRUMMP::SwapManager2D Swapper(pSD2D, this);
  (void) Swapper.swapAllFaces();

  return true;
}
