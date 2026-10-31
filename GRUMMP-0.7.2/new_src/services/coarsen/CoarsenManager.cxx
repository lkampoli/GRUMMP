/*
 * CoarsenManager.cxx
 *
 *  Created on: 2016-08-08
 */

#include <stdlib.h>

#include "GR_assert.h"
#include "GR_CoarsenManager.h"
#include "GR_Geometry.h"
#include "GR_Length.h"
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include "GR_misc.h"
#include "GR_SmoothingManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Vec.h"
#include "GR_VertConnect.h"

#define STRATEGY 2

namespace GRUMMP
{
  int
  CoarsenManager::markChain(Vert * const pVSurf, Vert * const pVFirst,
			    const VertConnect aVC[], const double dMaxDist,
			    const int iSkipNum)
  {
    pVSurf->markToKeep(); // Be sure -not- to delete this vertex.
    bool qChainOK = true;
    int iRetVal = 0;
    int iCurr = m_mesh->getVertIndex(pVFirst);
    Vert *pVLastKept = pVSurf;
    Vert *pVCurr = pVFirst;
    Vert *pVPrev = pVSurf;
    Vert *pVNext = pVInvalidVert;
    Vert *apVSkipped[8]; // Store the intermediate skipped verts, so that
    // those after the last kept vertex can be restored.

    if (pVSurf->getSpaceDimen() == 2)
      logMessage(4, "Chain starting from location (%10g, %10g), type %d\n",
		 pVSurf->x(), pVSurf->y(), pVSurf->getVertType());
    else
      logMessage(4,
		 "Chain starting from location (%10g, %10g, %10g), type %d\n",
		 pVSurf->x(), pVSurf->y(), pVSurf->z(), pVSurf->getVertType());

    while (1) { // Exit from this loop is via the return, below.
      int iSkipped = 0;
      while (iSkipped < iSkipNum && qChainOK) {
	assert(!pVPrev->isActive());
	assert(pVCurr->isActive());
	qChainOK = aVC[iCurr].qNormalNeighborExists(pVCurr, pVPrev, pVNext);
	pVCurr->markStructured();
	// Make sure nothing happens to this vert later.
	pVCurr->markInactive();
	apVSkipped[iSkipped] = pVCurr;
	if (pVCurr->getSpaceDimen() == 2)
	  logMessage(4, "Chain skipping location (%10g, %10g), type %d\n",
		     pVCurr->x(), pVCurr->y(), pVCurr->getVertType());
	else
	  logMessage(4, "Chain skipping location (%10g, %10g, %10g), type %d\n",
		     pVCurr->x(), pVCurr->y(), pVCurr->z(),
		     pVCurr->getVertType());

	if (!pVNext->isValid() || !pVNext->isActive())
	  qChainOK = false;
	if (qChainOK) {
	  pVPrev = pVCurr;
	  pVCurr = pVNext;
	  iCurr = m_mesh->getVertIndex(pVCurr);
	  iSkipped++;
	}
      }
      if (!qChainOK || (calcDistanceBetween(pVLastKept, pVCurr) > dMaxDist)) {
	logMessage(4, "Chain done\n");
	// Set length scale for last vertex and tag all its neighbors to
	// have their LS checked.  This prevents undesirable clustering
	// near the ends of chains.

	//pVCurr->setLengthScale(dMaxDist);
	m_length->setLengthScale(pVCurr);

	std::set<Vert*> spVNeigh;
	std::set<Cell*> spCTmp;
	findNeighborhoodInfo(pVCurr, spCTmp, spVNeigh);

	//m_vertsToUpdateLenScale.insert(spVNeigh.begin(), spVNeigh.end());
	for (Vert* pVneigh : spVNeigh) {
	  m_length->setLengthScale(pVneigh);
	}

	// Restore the last batch of skipped verts to a state where they
	// can be included in the MIS for this dimension.
	for (int ii = std::max(0, iSkipped - 1); ii > (iSkipNum + 1) / 2;
	    ii++) {
	  apVSkipped[ii]->markActive();
	}
	return iRetVal;
      }
      // For chains lying -on- the boundary, vertices are kept.  Chains
      // in the interior are deleted and re-inserted afterwards.
      pVCurr->markToKeep(); // Keep this vertex

      pVLastKept = pVCurr;
      iRetVal++;
      if (pVCurr->getSpaceDimen() == 2)
	logMessage(4, "Chain keeping location (%10g, %10g), type %d\n",
		   pVCurr->x(), pVCurr->y(), pVCurr->getVertType());
      else
	logMessage(4, "Chain keeping location (%10g, %10g, %10g), type %d\n",
		   pVCurr->x(), pVCurr->y(), pVCurr->z(),
		   pVCurr->getVertType());

      pVCurr->markStructured();
      // Make sure nothing happens to this vert later.
      pVCurr->markInactive();

      // Need to be pointing at the first vert to skip before returning to
      // the top of the loop.
      qChainOK = aVC[iCurr].qNormalNeighborExists(pVCurr, pVPrev, pVNext);
      if (!pVNext->isValid() || !pVNext->isActive())
	qChainOK = false;
      if (qChainOK) {
	pVPrev = pVCurr;
	pVCurr = pVNext;
	iCurr = m_mesh->getVertIndex(pVCurr);
      }
      else {
	logMessage(4, "Chain done\n");
	// Set length scale for last vertex and tag all its neighbors to
	// have their LS checked.

	//pVCurr->setLengthScale(dMaxDist);
	m_length->setLengthScale(pVCurr);

	std::set<Vert*> spVNeigh;
	std::set<Cell*> spCTmp;
	findNeighborhoodInfo(pVCurr, spCTmp, spVNeigh);

	//m_vertsToUpdateLenScale.insert(spVNeigh.begin(), spVNeigh.end());
	for (Vert* pVneigh : spVNeigh) {
	  m_length->setLengthScale(pVneigh);
	}

	// Restore the last batch of skipped verts to a state where they
	// can be included in the MIS for this dimension.
	for (int ii = std::max(0, iSkipped - 1); ii > (iSkipNum + 1) / 2;
	    ii++) {
	  apVSkipped[ii]->markActive();
	}
	return iRetVal;
      }
    }
    return iRetVal;
  }

  void
  CoarsenManager::deactivateChain(Vert * const pVSurf, Vert * const pVFirst,
				  const VertConnect aVC[],
				  const double dMaxDist) const
  {
    bool qChainOK = true;
    int iCurr = m_mesh->getVertIndex(pVFirst);
    Vert *pVCurr = pVFirst;
    Vert *pVPrev = pVSurf;
    Vert *pVNext;

    do {
      assert(!pVPrev->isActive());
      pVCurr->markStructured();
      pVCurr->markInactive();
      qChainOK = aVC[iCurr].qNormalNeighborExists(pVCurr, pVPrev, pVNext);
      if (qChainOK) {
	pVPrev = pVCurr;
	pVCurr = pVNext;
	iCurr = m_mesh->getVertIndex(pVCurr);
      }
    }
    while (pVCurr->isActive() && qChainOK
	&& calcDistanceBetween(pVCurr, pVPrev) < dMaxDist);
  }

  void
  CoarsenManager::addNormalVerts(const VertConnect aVC[],
				 const Vert::VertType VTVertTypeToStartFrom)
  {
    GR_index_t iV;
    int iAdded = 0;
    static const double dDistRatioLimit = 0.1;

    // Coarsen by choosing every iNSkip+1 vertex along valid lines.
    for (iV = 0; iV < m_mesh->getNumVerts(); ++iV) {
      Vert *pV = m_mesh->getVert(iV);
      // This vertex can be marched from only if (a) it is a sufficiently
      // important vertex (part of a sufficiently low dimension feature of
      // the mesh) and (b) the part of the mesh being marched into really
      // is anisotropic; this rule is actually a bit more complicated for
      // verts on BdryCurves.

      // Is this vertex part of a low-dimensional enough feature?
      if (pV->getVertType() <= VTVertTypeToStartFrom) {
	// Yes, it is.  Now find the distance to the nearest vertex that
	// is also low-D enough.  This is the long axis of the anisotropic
	// mesh, locally.
	double dLongDist = aVC[iV].dClosestLegitSurfPoint(
	    pV, VTVertTypeToStartFrom);

	int iNSkip = getInteriorSkip();
	bool qInFill = fillBoundary();

	if (VTVertTypeToStartFrom <= Vert::eBdryCurve) {
	  // Find the distance to the two closest surface points and the
	  // closest interior point.
	  double dDistS1 = LARGE_DBL, dDistS2 = LARGE_DBL, dDistInt =
	  LARGE_DBL;

	  if (dDistInt > dDistS1 && dDistS1 < dDistRatioLimit * dLongDist) {
	    // This is a case where an anisotropic mesh is present along
	    // one or both surfaces that meet at a fold.  March along
	    // lines on the surface keeping every fourth vertex to match
	    // what will happen in the interior of the mesh.  Only march
	    // from curve verts that are marked to be kept.
	    iNSkip = 3;
	    qInFill = false;
	  }
	  // For all other cases where the nearest surface vertex is
	  // close, we're going to keep every second vertex on the
	  // surface.
	  iNSkip = 1;
	  // If the two surfaces meeting at this curve are both
	  // anisotropically meshed, we'll fill in along the curve.
	  // Otherwise, we won't.
	  qInFill = (dDistS2 < dDistRatioLimit * dLongDist);
	}

	// When qInFill is set, keep all verts that have a chain coming
	// out from them.
	bool qInitialDeleteRequest = pV->isDeletionRequested();
	GR_index_t iNewlyAdded = 0;
	if (qInFill && qInitialDeleteRequest) {
	  pV->markToKeep();
	  iNewlyAdded++;
	}

	// Only check neighbors as starting points of chains.
	std::set<Vert*> spVNeigh;
	{
	  std::set<Cell*> spCTmp;
	  findNeighborhoodInfo(pV, spCTmp, spVNeigh);
	}

	std::set<Vert*>::iterator iter;
	for (iter = spVNeigh.begin(); iter != spVNeigh.end(); ++iter) {
	  Vert *pVCand = *iter;
	  if (!pVCand->isActive())
	    continue;
	  int iAddedThisChain = 0;
	  double dDistRatio =
	  dDIST3D(pVCand->getCoords(), pV->getCoords()) / dLongDist;
	  if (dDistRatio < dDistRatioLimit) {
	    pV->markStructured();
	    if (!pV->isDeletionRequested()) {
	      iAddedThisChain = markChain(pV, pVCand, aVC, dLongDist, iNSkip);
	      iNewlyAdded += iAddedThisChain;
	    }
	  } // Done running a chain out in one direction
	    // The following clause obliterates chains that didn't manage
	    // to build anything (including those where spacing was too
	    // large initially), as well as chains that weren't meant to
	    // build anything in the first place.
	  if (iAddedThisChain == 0 && dDistRatio < 0.5) {
	    deactivateChain(pV, pVCand, aVC, dLongDist);
	  }
	} // Done checking all active neighbors of this source vert
	  // If the vertex was tagged for inclusion in the fine mesh only
	  // because of a (non-existant) chain that could be built from it,
	  // then re-mark the vertex for deletion.

	if (iNewlyAdded == 1 && qInitialDeleteRequest) {
	  pV->markToDelete();
	  iNewlyAdded--;
	}
	logMessage(4, "From vertex %u, added %u verts in chain.\n", iV,
		   iNewlyAdded);
	iAdded += iNewlyAdded;
      } // Done checking verts that can have a PSAM line built from them
    } // Done checking all verts
    logMessage(2, "Added %d verts in pseudo-structured anisotropic mesh\n",
	       iAdded);
  }

  void
  CoarsenManager::addToCoarseVerts(const VertConnect aVC[]) const
  {
    GR_index_t iV;
    int iChanges = 0;

    // The initial MIS is set up using a frontal advance through the mesh,
    // beginning from all kept verts, and working outward a layer at a
    // time.  For BdryCurves in 3D and Bdrys in 2D, this will produce
    // near-maximal size on the first pass, which lexicographic ordering
    // of verts does not accomplish and which subsequent improvement can
    // not achieve because there isn't enough wiggle room to move vertices
    // around.

    {
      set<Vert*> spVChecked;
      for (iV = 0; iV < m_mesh->getNumVerts(); ++iV) {
	Vert *pV = m_mesh->getVert(iV);
	if (!pV->isDeletionRequested()) {
	  spVChecked.insert(pV);
	}
      }

      // Now do the frontal advance...
      //    int iOldLength = 0;
      set<Vert*> spVCurrentLayer;
      // Make an array copy of the latest layer of verts checked; this is
      // convenient to have, especially as the identity of new layers is
      // easiest to produce as an array.
      Vert** apVPrevLayer = new Vert*[spVChecked.size()];
      Vert** ppVEnd = std::copy(spVChecked.begin(), spVChecked.end(),
				apVPrevLayer);
      assert(
	  static_cast<unsigned int>(ppVEnd - apVPrevLayer)
	      == spVChecked.size());
      size_t iPrevSize = spVChecked.size();
      do {
//        int iStop = iOldLength;
//        iOldLength = LpVChecked.iLength();
//        for (iV = LpVChecked.iLength() - 1; iV >= iStop; iV--) {
//	Vert *pV = LpVChecked[iV];
	spVCurrentLayer.clear();
	for (size_t ii = 0; ii < iPrevSize; ii++) {
	  Vert *pV = apVPrevLayer[ii];
	  int iVert = m_mesh->getVertIndex(pV);
	  // Add all active neighbors of the current vert to the list, and
	  // mark each of these verts for inclusion in the coarse mesh if
	  // there isn't a conflicting vert already included.
	  for (int iii = 0; iii < int(aVC[iVert].iSize()); iii++) {
	    Vert *pVNeigh = aVC[iVert].pVNeigh(iii);
	    if (pVNeigh->isActive() && pVNeigh->isDeletionRequested()) {
//	    LpVChecked.vAddItem(pVNeigh);
	      // Add to the list of all checked vertices and to the list
	      // of vertices checked on this pass.
	      spVCurrentLayer.insert(pVNeigh);
	      int iNeigh = m_mesh->getVertIndex(pVNeigh);
	      if (aVC[iNeigh].iCountKeptNeigh() == 0) {
		pVNeigh->markToKeep();
		iChanges++;
	      }
	    }
	  } // Done looping over one listed verts neighbors
	} // Done looping over all verts inserted in the last layer of the
	  // advance

	  // Put a list of the newest layer of verts into apVPrevLayer.
	delete[] apVPrevLayer;
	apVPrevLayer = new Vert*[spVCurrentLayer.size()];
	ppVEnd = set_difference(spVCurrentLayer.begin(), spVCurrentLayer.end(),
				spVChecked.begin(), spVChecked.end(),
				apVPrevLayer);
	iPrevSize = ppVEnd - apVPrevLayer;

	// Now update the list of all verts that have been checked.
	spVChecked.insert(apVPrevLayer, ppVEnd);
      }
      while (iPrevSize > 0);
    } // Done with first insertion pass over vertices.
    logMessage(2,
	       "Initial site insertion pass for this phase added %7d verts\n",
	       iChanges);

    int iIter = 0;
    do {
      iIter++;
      iChanges = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); ++iV) {
	Vert *pV = m_mesh->getVert(iV);
	if (pV->isActive() && pV->isDeletionRequested()) {
	  int iNumKeptNeigh = aVC[iV].iCountKeptNeigh();

	  if (iNumKeptNeigh == 0) {
	    pV->markToKeep();
	    ++iChanges;
	  }
	  else if (iNumKeptNeigh == 1) {
	    int iKept2nd = iCountKept2nd(iV, aVC, m_mesh);
	    for (int iNeigh = 0; iNeigh < int(aVC[iV].iSize()); iNeigh++) {
	      Vert *pVNeigh = aVC[iV].pVNeigh(iNeigh);
	      if (pVNeigh->isActive() && !pVNeigh->isDeletionRequested()) {
		int iDiff = iKept2nd
		    - iCountKept2nd(m_mesh->getVertIndex(pVNeigh), aVC, m_mesh);
		double dThresh = atan((double(iDiff) - 0.5) * 5.) / M_PI + 0.5;
#ifndef WIN_VC32
		if (drand48() < dThresh) {
#else
		  if ( rand()/double(RAND_MAX) < dThresh ) {
#endif
		  assert(pV->isActive());
		  assert(pVNeigh->isActive());
		  pVNeigh->markToDelete();
		  pV->markToKeep();
		}
		break;
	      } /* neighbor marked for keeping in the coarse mesh */
	    } /* loop over first neighbors */
	  } /* exactly one kept neighbor */
	} /* current vertex to be deleted */
      } /* loop over all vertices */
      logMessage(2, "Site insertion pass%3d added%7d verts\n", iIter, iChanges);
    }
    while (iIter < 20 && iChanges); /* iteration loop */

    // End with a pass over all vertices to ensure that the marked
    // vertices are a maximal independent set.  There is a small chance
    // that movement will result in a non-maximal set otherwise.
    for (iV = 0; iV < m_mesh->getNumVerts(); ++iV) {
      Vert *pV = m_mesh->getVert(iV);
      if (pV->isActive() && pV->isDeletionRequested()) {
	int iNumKeptNeigh = aVC[iV].iCountKeptNeigh();
	if (iNumKeptNeigh == 0) {
	  pV->markToKeep();
	  ++iChanges;
	}
      }
    } // Done with final insertion pass over vertices.
    logMessage(2, "Final site insertion pass for this phase added %7d verts\n",
	       iChanges);
  }

  void
  CoarsenManager::deactivateAll() const
  {
    for (GR_index_t iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      m_mesh->getVert(iV)->markInactive();
    }
  }

  int
  CoarsenManager::activate(const Vert::VertType VT) const
  {
    int iRetVal = 0;
    for (GR_index_t iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      Vert *pV = m_mesh->getVert(iV);
      if (pV->getVertType() == VT) {
	m_mesh->getVert(iV)->markActive();
	iRetVal++;
      }
    }
    return iRetVal;
  }

  void
  CoarsenManager::selectCoarseVerts(VertConnect aVC[])
// This routine is designed to choose which vertices in a mesh should
// remain when the mesh is coarsened.  These vertices are tagged using
// Vert::vMarkToKeep; those that will be deleted are tagged using
// Vert::vMarkToDelete.
//
// Selection is heirarchical, working first on bdry curves, then on
// the entire boundary, then in the interior of the domain.  At each
// level, first a check is made for stretched (pseudo-structured) mesh
// regions, by building chains from verts already tagged for retention
// into the next higher-dimensional entity.
  {
    logMessage(1, "Deciding which verts to keep for the coarse mesh...\n");
    // Tag every vertex as not ACTIVE and DELETE_REQUESTED (except for
    // bdry apexes, which must always be keep for geometric integrity).
    GR_index_t iV;
    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      Vert *pV = m_mesh->getVert(iV);
      pV->markInactive();
      if (pV->getVertType() == Vert::eBdryApex || pV->getVertType() == Vert::ePseudoSurface) {
	pV->markToKeep();
	logMessage(4, "Vertex %u is an apex; kept.\n", iV);
      }
    }

    // First, find verts to add along bdry curves.  This will apply only
    // for 3D (and possible surface; that isn't clear yet), and checking
    // the number of Verts found on BdryCurves will short-circuit the real
    // work here when appropriate.
    int iActive = activate(Vert::eBdryCurve);
    if (iActive > 0) {
      logMessage(2, " Deciding which verts to keep along bdry curves...\n");
      // Conflict graph must be updated, because otherwise we don't take
      // length scale into account.
      updateConflictGraph(aVC);
      // Conflict graph must be updated, because otherwise we don't take
      // length scale into account.
      updateConflictGraph(aVC);
      // Add a MIS of maximal size within the active verts.
      addToCoarseVerts(aVC);

      // Lock all verts.
      deactivateAll();
    }
    {
      int iCurrentNumberToDelete = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	if (m_mesh->getVert(iV)->isDeletionRequested()
	    && !m_mesh->getVert(iV)->isDeleted()) {
	  iCurrentNumberToDelete++;
	}
      }
      logMessage(3, "After bdry curves: %d to remove.\n",
		 iCurrentNumberToDelete);
    }

    // Now work with the rest of the boundary, including internal
    // boundaries and pseudo-surfaces.
    iActive = (activate(Vert::eBdryTwoSide) + activate(Vert::eBdry)
	+ activate(Vert::ePseudoSurface));
    if (iActive > 0) {
      logMessage(2, " Deciding which verts to keep on boundaries...\n");
      addNormalVerts(aVC, Vert::eBdryCurve);
//updateLengthScale();
      m_length->setLengthScale(m_mesh);
      updateConflictGraph(aVC);
      addToCoarseVerts(aVC);
      deactivateAll();
    }
    {
      int iCurrentNumberToDelete = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	if (m_mesh->getVert(iV)->isDeletionRequested()
	    && !m_mesh->getVert(iV)->isDeleted()) {
	  iCurrentNumberToDelete++;
	}
      }
      logMessage(3, "After all bdrys: %d to remove.\n", iCurrentNumberToDelete);
    }

    // In regions where the volume mesh is structured (STRUCTURED),
    // choose every fourth vertex along structured lines.
    iActive = activate(Vert::eInterior);
    if (iActive > 0) {
      logMessage(2, " Deciding which verts to keep in the interior...\n");
      // March out from surface-like stuff, skipping three vertices between
      // kept verts (keeping 1 out of 4) to equalize aspect ratios sooner.
      addNormalVerts(aVC, Vert::ePseudoSurface);
//updateLengthScale();
      m_length->setLengthScale(m_mesh);
      updateConflictGraph(aVC);
      addToCoarseVerts(aVC);
      deactivateAll();
    }
    {
      int iCurrentNumberToDelete = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	if (m_mesh->getVert(iV)->isDeletionRequested()
	    && !m_mesh->getVert(iV)->isDeleted()) {
	  iCurrentNumberToDelete++;
	}
      }
      logMessage(3, "After interior: %d to remove.\n", iCurrentNumberToDelete);
    }

//    vMessage(2, " Cleaning up along boundaries...\n");
//    // To improve coarse mesh quality near boundaries, retain any deleted
//    // non-STRUCTURED boundary vertex that is in conflict with a kept,
//    // non-STRUCTURED interior vertex.
//    for (iV = 0; iV < iNumVerts(); iV++) {
//      Vert *pV = pVVert(iV);
//      if (pV->qIsBdryVert() && pV->qDeletionRequested() &&
//	!pV->qIsStructured()) {
//        for (int ii = 0; ii < aVC[iV].iSize(); ii++) {
//	Vert *pVNeigh = aVC[iV].pVNeigh(ii);
//	if (!pVNeigh->qIsBdryVert() && !pVNeigh->qDeletionRequested() &&
//	    !pVNeigh->qIsStructured()) {
//	  pV->vMarkToKeep();
//	}
//        }
//      }
//    }
    logMessage(1, "Done deciding which verts to keep for the coarse mesh...\n");
  }

  void
  CoarsenManager::removeTaggedVerts(bool smooth, std::set<Vert*> *cantRemove)
  {
    // Remove vertices from the mesh, with intermittent swapping to ease
    // the process and get good mesh quality.
    GR_index_t iDeletedThisPass, iDeletedSincePurge, iTotalDeleted;
    GR_index_t iNumSwaps, iNMarkedForDeletion, iV;
    int iPasses;

    iDeletedSincePurge = iTotalDeleted = iPasses = iNumSwaps =
	iNMarkedForDeletion = 0;
    // TODO: This is only going to work for 2D.  Unless the swap managers can be made
    // to share a common base class, this function may have to be split.  That
    // basically means that there'll have to be separate coarsening managers for
    // 2D and 3D.

    GRUMMP::DelaunaySwapDecider2D * SD2D = NULL;
    GRUMMP::DelaunaySwapDecider3D * SD3D = NULL;
    GRUMMP::SwapManager2D * SM2D = NULL;
    GRUMMP::SwapManager3D * SM3D = NULL;

    if (m_mesh->getType() == Mesh::eMesh2D) {
      SD2D = new GRUMMP::DelaunaySwapDecider2D();
      SM2D = new GRUMMP::SwapManager2D(SD2D, dynamic_cast<Mesh2D*>(m_mesh));
    }
    else {
      // Better be a VolMesh
      assert(m_mesh->getType() == Mesh::eVolMesh);
      SD3D = new GRUMMP::DelaunaySwapDecider3D(true); // Bdry is changeable
      SM3D = new GRUMMP::SwapManager3D(SD3D, dynamic_cast<VolMesh*>(m_mesh));
    }

//   If one were doing Delaunay vertex removal (which is non-robust),
//   then Delaunay swapping would be the right thing to do.  Otherwise,
//   it just wrecks the connectivity so that you can't remove verts as
//   well as you otherwise would.
//   vSetSwapType(eDelaunay);
    if (m_mesh->getSwapType() != eMinSine) {
      if (SM2D != NULL) {
	SM2D->swapAllFaces();
      }
      else {
	assert(SM3D != NULL);
	SM3D->swapAllFaces();
      }
    }

    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      logMessage(3, "Vertex %u at (%f, %f) is%s tagged for deletion.\n", iV,
		 m_mesh->getVert(iV)->x(), m_mesh->getVert(iV)->y(),
		 m_mesh->getVert(iV)->isDeletionRequested() ? "" : " not");
      if (m_mesh->getVert(iV)->isDeletionRequested())
	iNMarkedForDeletion++;

    }

    logMessage(1, "Removing vertices from the mesh using strategy %d.\n",
    STRATEGY);
    bool qBdryOnly = false;
    GRUMMP::OptMSSmoothingManager* OMS = GRUMMP::OptMSSmoothingManager::Create(
	m_mesh);

    do {
      iPasses++;
      if (iPasses >= 3)
	qBdryOnly = false;

      logMessage(1, "  Pass %d:  ", iPasses);
      iDeletedThisPass = 0;
      GR_index_t iPrintInt = max(GR_index_t(25),
				 m_mesh->getNumVerts() / 10 + 1);

      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	logMessage(4, "Checking vertex %u for removal\n", iV);
	Vert *pV = m_mesh->getVert(iV);

#if (STRATEGY == 3)
	if (pV->isDeleted() ||
	    (!pV->isDeletionRequested()) ||
	    (!pV->isBdryVert() && qBdryOnly)) continue;
	// Contract only onto the nearest kept vertex.
	std::set<Cell*> spCDummy;
	std::set<Vert*> spVNeigh;
	findNeighborhoodInfo(pV, spCDummy, spVNeigh);

	// Find the nearest kept vertex.
	double dMinDist = LARGE_DBL;
	Vert *pVClosest = pVInvalidVert;
	for (std::set<Vert*>::iterator iterV = spVNeigh.begin();
	    iterV != spVNeigh.end(); ++iterV) {
	  Vert *pVCand = *iterV;
	  if (!pVCand->isDeletionRequested()) {
	    double dDist = dDIST3D(pVCand->getCoords(), pV->getCoords());
	    if (dDist < dMinDist) {
	      dMinDist = dDist;
	      pVClosest = pVCand;
	    }
	  }
	}
	// If there -is- a kept neighbor, try the deletion.
	int iNewSwaps;
	if (pVClosest->isValid() &&
	    removeVertByContraction(pV, iNewSwaps, &pVClosest, 1)) {
	  iNumSwaps += iNewSwaps;
	  iTotalDeleted++;
	  iDeletedThisPass++;
	  iDeletedSincePurge++;
	} // Deleted a vert.
#elif (STRATEGY == 2)
	if (pV->isDeleted() || (!pV->isBdryVert() && qBdryOnly))
	  continue;
	int iNewSwaps;
	if (pV->isDeletionRequested()
	    && m_mesh->removeVert(pV, iNewSwaps)) {
	  m_mesh->sendEvents();
	  iNumSwaps += iNewSwaps;
	  iTotalDeleted++;
	  iDeletedThisPass++;
	  iDeletedSincePurge++;
	} // Deleted a vert.
#elif (STRATEGY == 1)
	// A new approach: for each kept vertex, try to delete all deletable
	// neighbors.
	// Under odd circumstances, a vertex whose deletion was not
	// requested may actually be removed from the mesh, so we must
	// check for this...
	if (!pV->isDeletionRequested() && !pV->isDeleted()) {
	  // Find all neighboring verts, then throw out the ones that
	  // don't need to be deleted.
	  std::set<Cell*> spCDummy;
	  std::set<Vert*> spVNeigh;
	  findNeighborhoodInfo(pV, spCDummy, spVNeigh);

	  std::set<Vert*>::iterator iterV;
	  for (iterV = spVNeigh.begin(); iterV != spVNeigh.end(); ) {
	    Vert *pVCand = *iterV;
	    assert(!pVCand->isDeleted());
	    if (!pVCand->isDeletionRequested()) {
	      std::set<Vert*>::iterator iterTmp = iterV;
	      ++iterV;
	      spVNeigh.erase(iterTmp);
	    }
	    else ++iterV;
	  }

	  // Go through the list repeatedly until no more deletions are
	  // made onto this vertex.
	  bool qChange;
	  do {
	    qChange = false;
	    for (iterV = spVNeigh.begin(); iterV != spVNeigh.end(); ) {
	      Vert *pVCand = *iterV;
	      std::set<Vert*>::iterator iterTmp = iterV;
	      ++iterV;
	      int iNewSwaps = 0;
	      if (pVCand->isDeleted() ||
		  removeVertByContraction(pVCand, iNewSwaps, &pV, 1)) {
		qChange = true;
		spVNeigh.erase(iterTmp);
		iNumSwaps += iNewSwaps;
		iTotalDeleted++;
		iDeletedThisPass++;
		iDeletedSincePurge++;
	      } // Deleted a vert.
	    } // Loop over neighboring verts.
	  }while (qChange && !spVNeigh.empty());
	} // Is this vert being kept?
#endif
	if (iV % iPrintInt == 0 && iDeletedThisPass != 0) {
	  logMessage(2, "At vertex %u out of %u; %u verts removed, %u swaps\n",
		     iV, m_mesh->getNumVerts(), iDeletedThisPass, iNumSwaps);
	}
      } // Loop over verts.
      assert(m_mesh->isValid());

      logMessage(2, "%6u vertices removed, %6u total.  %6u still to remove\n",
		 iDeletedThisPass, iTotalDeleted,
		 iNMarkedForDeletion - iTotalDeleted);
      m_mesh->sendEvents();
      // If nothing was deleted, this would be a waste of time.
      if (iDeletedThisPass != 0) {
	logMessage(2, "  Smoothing to improve cell shapes...\n");
	if (smooth)
	  OMS->smoothAllVerts(1);
	logMessage(2, "  Swapping to improve connectivity...\n");
	if (SM2D != NULL) {
	  iNumSwaps += SM2D->swapAllFaces();
	}
	else {
	  assert(SM3D != NULL);
	  iNumSwaps += SM3D->swapAllFaces();
	}
	// Using 3D bad cell repair techniques does no good here.
      }
    }
    while (iTotalDeleted < iNMarkedForDeletion && iDeletedThisPass > 0);
    logMessage(2, "Deleted %u vertices out of %u marked for deletion.\n",
	       iTotalDeleted, iNMarkedForDeletion);
    logMessage(2, "A total of %u swaps were performed %s\n", iNumSwaps,
	       "to allow deletion and improve mesh quality.");

    // Because the data structures are now re-using old space, the need to
    // purge on the fly is much less, and purging really messes up the whole
    // Observer thing.
//  if (iDeletedSincePurge != 0) purgeAllEntities();
    delete OMS;

    // Info on undeleted verts.
    logMessage(2, "The following verts could not be deleted:\n");
    if (iTotalDeleted < iNMarkedForDeletion) {
      GR_index_t iBdryCurve = 0, iBdryStruct = 0, iBdry = 0, iBdryTotal = 0;
      GR_index_t iIntStruct = 0, iInt = 0, iIntTotal = 0, iGrandTotal = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	Vert *pV = m_mesh->getVert(iV);
	if (pV->isDeleted())
	  continue;
//      assert(!pV->isDeleted());
	if (pV->isDeletionRequested()) {
	  switch (pV->getVertType())
	    {
	    case Vert::eBdryCurve:
	      iBdryCurve++;
	      iBdryTotal++;
	      iGrandTotal++;
	      break;
	    case Vert::eBdry:
	    case Vert::eBdryTwoSide:
	      if (pV->isStructured())
		iBdryStruct++;
	      else
		iBdry++;
	      iBdryTotal++;
	      iGrandTotal++;
	      break;
	    case Vert::eInterior:
	      if (pV->isStructured())
		iIntStruct++;
	      else
		iInt++;
	      iIntTotal++;
	      iGrandTotal++;
	      break;
	    default:
	      assert(0);
	      break;
	    }
	  if (cantRemove)
	    cantRemove->insert(pV);
	  if (pV->getSpaceDimen() == 2) {
	    logMessage(2, "Vert index %6u at (%13.6g, %13.6g). Type: %d\n", iV,
		       pV->x(), pV->y(), pV->getVertType());
	  }
	  else {
	    logMessage(2,
		       "Vert index %6u at (%13.6g, %13.6g, %13.6g). Type: %d\n",
		       iV, pV->x(), pV->y(), pV->z(), pV->getVertType());
	  }
	}
      }
      assert(iGrandTotal == iNMarkedForDeletion - iTotalDeleted);
      assert(iGrandTotal == iIntTotal + iBdryTotal);
      assert(iBdryTotal == iBdry + iBdryStruct + iBdryCurve);
      assert(iIntTotal == iInt + iIntStruct);
      logMessage(1, "Total verts not removed: %u\n", iGrandTotal);
      logMessage(
	  1,
	  " On the bdry: %3u on curves, %3u structured, %3u others. Total: %u\n",
	  iBdryCurve, iBdryStruct, iBdry, iBdryTotal);
      logMessage(
	  1,
	  " In interior:                %3u structured, %3u others. Total: %u\n",
	  iIntStruct, iInt, iIntTotal);
    }

    // Finish up.
    if (!m_mesh->isValid()) {
      vFoundBug("removing verts to create a coarse mesh\n");
    }
    if (SM2D) {
      delete SM2D;
      assert(SD2D);
      delete SD2D;
    }
    else {
      assert(SM3D);
      delete SM3D;
      assert(SD3D);
      delete SD3D;
    }

  }

  void
  CoarsenManager::updateConflictGraph(VertConnect aVC[]) const
// The goal is to create a list, for each vertex i, of vertices j for
// which the distance between V_i and V_j is smaller than the sum of
// their length scales. An empirical scaling factor is used here to
// make sure approximately the right number of verts remain after
// coarsening.
//
// To accomplish this, a search marches outward from each vertex,
// seeking all points that are topologically near a given vertex and
// also geometrically nearby.  A strictly geometric search pierces
// thin boundaries, which is bad, so that approach has been
// abandoned.
//
// Only vertices that are current marked active are considered; others
// are either already taken care of or not yet relevant.
  {
    logMessage(
	1, "Updating vertex-to-vertex conflict graph using length scale...\n");
    GR_index_t iV, iInfoInterval = max(GR_index_t(25),
				       m_mesh->getNumVerts() / 10);
    double dBubbleSizeConst = m_mesh->getCoarseningConstant();
    GR_index_t iUpdated = 0;
    GR_index_t uiMaxCand = 100;

    Vert **apVCand = new Vert*[uiMaxCand];

    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      if (iV % iInfoInterval == 0) {
	logMessage(2, "  Updating conflicts for vertex %u...\n", iV);
      }
      Vert *pV = m_mesh->getVert(iV);
      if (!pV->isActive())
	continue;

      iUpdated++;

      // Initially, create a conflict set containing only pV, and a
      // candidate array containing the current conflict list for pV.
      // Also, initialize the set of tested verts to pV.

      set<Vert*> spVConflict, spVTested;
      spVConflict.insert(pV);
      spVTested.insert(pV);

      size_t iNCand = aVC[iV].iSize();
      if (uiMaxCand < iNCand) {
	uiMaxCand = iNCand;
	delete[] apVCand;
	apVCand = new Vert*[uiMaxCand];
      }
      Vert **ppVFirst = aVC[iV].ppVVerts();
      Vert **ppVLast = std::copy(ppVFirst, ppVFirst + iNCand, apVCand);
      assert(size_t(ppVLast - apVCand) == iNCand);

      int iPass = 0;
      set<Vert*> spVNewCand;
      do {
	spVNewCand.clear();
	spVTested.insert(apVCand, apVCand + iNCand);
	iPass++;

	for (size_t i = 0; i < iNCand; i++) {
	  Vert *pVCand = apVCand[i];
	  if (!pVCand->isActive())
	    continue;

	  double dDist = calcDistanceBetween(pV, pVCand);

	  // The definition of dMinLegalDist has been tweaked so that, in
	  // 2D, the mesh is about one-fourth as dense as before
	  // coarsening and the point distribution for coarse meshes is
	  // esthetically reasonable.
	  double dMinLegalDist = dBubbleSizeConst
	      * (pV->getLengthScale() + pVCand->getLengthScale());

	  if (dDist <= dMinLegalDist) {
	    // Conflict!  Add to conflict list and be sure to check
	    // anything in conflict with it.
	    logMessage(
		3,
		"Verts %u and %u. Dist = %12.6f, min legal = %12.6f, ratio = %12.6f\n",
		m_mesh->getVertIndex(pV), m_mesh->getVertIndex(pVCand), dDist,
		dMinLegalDist, dDist / dMinLegalDist);
	    spVConflict.insert(pVCand);
	    int iCand = m_mesh->getVertIndex(pVCand);
	    Vert **ppVBegin = aVC[iCand].ppVVerts();
	    Vert **ppVEnd = ppVBegin + aVC[iCand].iSize();
	    spVNewCand.insert(ppVBegin, ppVEnd);
	  } // Done with conflicting candidate
	} // Done checking all candidates
	if (spVNewCand.size() > uiMaxCand) {
	  delete[] apVCand;
	  uiMaxCand = spVNewCand.size() * 2;
	  apVCand = new Vert*[uiMaxCand];
	}
	ppVLast = std::set_difference(spVNewCand.begin(), spVNewCand.end(),
				      spVTested.begin(), spVTested.end(),
				      apVCand);
	iNCand = ppVLast - apVCand;
      }
      while (iNCand > 0); // No more candidates to check!
      spVConflict.erase(pV);

      // Construct a VertConnect entry for vertex iV
      logMessage(3, "Vertex %d: number of conflicts = %zd\n", iV,
		 spVConflict.size());
      aVC[iV].vSetSize(spVConflict.size());
      std::set<Vert*>::iterator it;
      size_t ii = 0;
      for (it = spVConflict.begin(); it != spVConflict.end(); it++, ii++) {
	aVC[iV].vSetNeigh(ii, *it);
      }
      assert(ii == aVC[iV].iSize());
    } // End of loop over verts
    logMessage(1, "Updated conflict graph for %u out of %u verts.\n", iUpdated,
	       m_mesh->getNumVerts());
    logMessage(2, "Maximum number of candidates considered at once was %d\n",
	       uiMaxCand);
    delete[] apVCand;
  }

  void
  CoarsenManager::initConflictGraph(VertConnect aVC[]) const
// The full conflict graph is overkill for coarsening boundary curves
// and for anisotropic coarsening.  In those cases, neighbor info is
// all that's needed.  So compute neighbor info only, and later on
// compute the true conflict graph incrementally using
// vUpdateConflictGraph.
  {
    logMessage(
	1,
	"Creating vertex-to-vertex conflict graph using neighbor data only...\n");
    GR_index_t iV, iInfoInterval = max(GR_index_t(25),
				       m_mesh->getNumVerts() / 10);

    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      if (iV % iInfoInterval == 0) {
	logMessage(2, "  Checking conflicts for vertex %u...\n", iV);
      }
      Vert *pV = m_mesh->getVert(iV);
      std::set<Vert*> spVConflict;
      {
	std::set<Cell*> spCTmp;
	findNeighborhoodInfo(pV, spCTmp, spVConflict);
      }
      spVConflict.erase(pV);

      // Construct a VertConnect entry for vertex iV
      size_t newSize = spVConflict.size();
      aVC[iV].vSetSize(newSize);
      std::set<Vert*>::iterator iter;
      int ii = 0;
      for (iter = spVConflict.begin(); iter != spVConflict.end();
	  ++iter, ii++) {
	aVC[iV].vSetNeigh(ii, *iter);
      }
    } // End of loop over verts
    logMessage(1, "Done initializing conflict graph.\n");
  }

// Create a coarse unstructured mesh.  dLengthRatio is the ratio by which
// to increase the length scale.
  void
  CoarsenManager::coarsen(const double dLengthRatio)
  {
    assert(m_mesh->isSimplicial());

    m_length->scale(dLengthRatio);
    m_length->setLengthScale(m_mesh);

    logMessage(2, "Done increasing length scale.\n");

    coarsenToLengthScale();
  }

  void
  CoarsenManager::coarsenToLengthScale(const bool qHeal)
  {
    assert(m_mesh->isSimplicial());
    bool qRebuildPatchesAfter = (m_mesh->getType() == Mesh::eVolMesh)
	&& (m_mesh->getNumBdryFaces() == m_mesh->getNumBdryPatches());
    logMessage(1, "Coarsening mesh...\n");

    GR_index_t iV;
    /*
     * Grading should already have been checked by Length object when setting the coarsened length scale.
     */
//  // Now make sure that grading, etc, is okay...
//  logMessage(2, "Updating length scale to take grading into account...\n");
//  for (iV = 0; iV < getNumVerts(); iV++) {
//    m_vertsToUpdateLenScale.insert(getVert(iV));
//  }
//  updateLengthScale();
//	logMessage(2, "Done with grading check.\n");
    // Classify verts on the basis of geometry only; this should prevent
    // degeneration of holes with multiple coarsening.
    // The following code is a little bit dubious, in that it could cause
    // problems with points that lie on patch boundaries.  So for now
    // just make sure that anything definitively ID'd as an apex from the
    // bdry patches isn't reclassified.
    // FIX ME 0.5.0
    m_mesh->identifyVertexTypesGeometrically();
    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      logMessage(4, "Vert: %5u  Type: %d\n", iV,
		 m_mesh->getVert(iV)->getVertType());
    }

    // Always keep bdry apexes.  Also, for the mesh healing case, keep all
    // boundary verts.  This should be okay unless the boundary itself has
    // changed shape rapidly, but it is a bit of a hack. CFO-G 10/17/02
    for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
      Vert *pV = m_mesh->getVert(iV);
      if (pV->getVertType() == Vert::eBdryApex) {
	pV->markToKeep();
      }
      else if (pV->isBdryVert() && qHeal) {
	pV->markToKeep();
      }
      else {
	pV->markToDelete();
      }
    }

    VertConnect *aVC = new VertConnect[m_mesh->getNumVerts()];
    initConflictGraph(aVC);
    {
      int iCurrentNumberToDelete = 0;
      for (iV = 0; iV < m_mesh->getNumVerts(); iV++) {
	if (m_mesh->getVert(iV)->isDeletionRequested()
	    && !m_mesh->getVert(iV)->isDeleted()) {
	  iCurrentNumberToDelete++;
	}
      }
      logMessage(3, "Before selection: %d to remove.\n",
		 iCurrentNumberToDelete);
    }

    selectCoarseVerts(aVC);

    // Remove connectivity information
    delete[] aVC;

    removeTaggedVerts();

    if (qRebuildPatchesAfter)
      m_mesh->buildBdryPatches();

    logMessage(1, "Done coarsening mesh\n");
  }
}
