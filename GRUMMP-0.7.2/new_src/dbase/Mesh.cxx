#include <queue>

#include "GR_config.h"
#include "GR_assert.h"
#include "GR_Mesh.h"
#include "GR_Vertex.h"
#include "GR_Face.h"

void
Mesh::reorderVertices(GR_index_t*& aRReorderInt)
{
  int iNVerts = getNumVerts(), iVert;
  //@@ Move vertex data
  // If we copy this to a new array and copy it back, all the pointers from
  // faces get updated automagically (well, actually, there's explicit code in
  // Vert::operator=).
  EntContainer<Vert> ECTemp;
  ECTemp.reserve(iNVerts);

  for (iVert = 0; iVert < iNVerts; iVert++) {
    ECTemp[iVert] = (*getVert(iVert));
  }

  for (iVert = 0; iVert < iNVerts; iVert++) {
    // Copy back into re-ordered locations
    Vert* pVSource = &(ECTemp[aRReorderInt[iVert]]);
    Vert* pVTarget = getVert(iVert);
    (*pVTarget) = (*pVSource);
  }
}

struct _reorder_Verts_Lex_ {
  // Pointers move through index of vertices and returns true if first argument is less than the second in x, then y, and then z

  bool
  operator()(const Vert * v1, const Vert * v2) const
  {
    return ((v1->x() < v2->x()) || (v1->x() == v2->x() && v1->y() < v2->y())
	|| (v1->x() == v2->x() && v1->y() == v2->y() && v1->z() < v2->z()));
  }
} reorderVertsLex;

void
Mesh::reorderVerts_Lex()
{
  std::set<const Vert *, _reorder_Verts_Lex_> Verts_Lex;
  std::set<const Vert *, _reorder_Verts_Lex_>::iterator iterLex;

  for (GR_index_t i = 0; i < getNumVerts(); i++) {
    Vert* vert = getVert(i);
    if (vert->isDeleted())
      continue;
    Verts_Lex.insert(vert);
  }

  GR_index_t * aRReorder = new GR_index_t[Verts_Lex.size()];

  GR_index_t iterVert = 0;

  for (iterLex = Verts_Lex.begin(); iterLex != Verts_Lex.end();
      ++iterLex, ++iterVert) {
    const Vert* thisVert = *iterLex;
    GR_index_t index = getVertIndex(thisVert);
    aRReorder[iterVert] = index;
  }
  reorderVertices(aRReorder);

  delete[] aRReorder;

  //@@ Reset vertex hints
  setAllHintFaces();
  setVertFaceNeighbors();

  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Lexicographic vertex re-ordering");
  }
  logMessage(2, "done.\n");
}

//@@ Compute hint faces for all vertices
void
Mesh::clearAllHintFaces()
{
  GR_index_t i;
  for (i = 0; i < getNumVerts(); i++)
    getVert(i)->clearHintFace();
}

//@@ Compute hint faces for all vertices
void
Mesh::setAllHintFaces()
{
  // Have to clear vertex hints -before- checking validity, else bad
  // hints will cause an irrelevant assertion.
  // logMessage(2, "Setting all hint faces...");
#ifndef NDEBUG
  clearAllHintFaces();
  //  assert(qValid());
#endif

  GR_index_t i;
  for (i = 0; i < getNumFaces(); i++) {
    Face* pF = getFace(i);
    if (pF->isDeleted() || !(pF->doFullCheck()))
      continue;
    for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
      Vert* pV = pF->getVert(ii);
      if (pV->isValid()) {
	pV->setHintFace(pF);
      }
    }
  }

  // Bdry faces are prefered as vertex hints, so overwrite as needed.
  for (i = 0; i < getNumBdryFaces(); i++) {
    BFace* pBF = getBFace(i);
    if (pBF->isDeleted())
      continue;
    Face* pF = pBF->getFace(0);
    if (pF->isDeleted())
      continue;
    for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
      Vert* pV = pF->getVert(ii);
      if (pV->isValid()) {
	pV->setHintFace(pF);
      }
    }
  }
  // logMessage(2, "done\n");
}

void
Mesh::setVertFaceNeighbors()
{
  // logMessage(2, "Setting all vert-face connectivity...");
  clearVertFaceNeighbors();
  // Have to clear vertex hints -before- checking validity, else bad
  // hints will cause an irrelevant assertion.
  //  assert(qValid());

  GR_index_t i;
  for (i = 0; i < getNumFaces(); i++) {
    Face* pF = getFace(i);
    if (pF->isDeleted())
      continue;
    for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
      Vert* pV = pF->getVert(ii);
      if (pV != pVInvalidVert) {
	pV->addFace(pF);
      }
    }
  }
  // logMessage(2, "done\n");
}

void
Mesh::clearVertFaceNeighbors()
{
  GR_index_t i;
  for (i = 0; i < getNumVerts(); i++)
    getVert(i)->clearFaceConnectivity();
}

void
findNeighborhoodInfo(const Vert* const pVert, std::set<Cell*>& spCInc,
		     std::set<Vert*>& spVNeigh, std::set<BFace*>* pspBFInc,
		     bool *pqBdryVert, std::set<Face*>* pspFNearby)
{
  spCInc.clear();
  spVNeigh.clear();
  if (pspBFInc) {
    pspBFInc->clear();
    assert(pqBdryVert != NULL);
    (*pqBdryVert) = false;
  }
  if (pspFNearby)
    pspFNearby->clear();

  Face * pF = pFInvalidFace;
  Cell * pC = pCInvalidCell;

  for (int iF = pVert->getNumFaces() - 1; iF >= 0; iF--) {
    pF = const_cast<Face*>(pVert->getFace(iF));
    assert(!pF->isDeleted());

    pC = pF->getLeftCell();
    if (pC->isValid() && !pC->isBdryCell()) {
      // In addition to doing the insertion, check to see
      // whether this cell had already been added.
      bool wasInserted = spCInc.insert(pC).second;
      if (wasInserted) {
	if (pspFNearby) {
	  for (int iCF = pC->getNumFaces() - 1; iCF >= 0; iCF--) {
	    pspFNearby->insert(pC->getFace(iCF));
	  }
	}
	for (int iV = pC->getNumVerts() - 1; iV >= 0; iV--) {
	  Vert *pVAdd = pC->getVert(iV);
	  if (!pVAdd->isDeleted())
	    spVNeigh.insert(pVAdd);
	}
      }
    }
    else if (pspBFInc) {
      pspBFInc->insert(dynamic_cast<BFace*>(pC));
    }
    pC = pF->getRightCell();
    if (pC->isValid() && !pC->isBdryCell()) {
      // In addition to doing the insertion, check to see
      // whether this cell had already been added.
      bool wasInserted = spCInc.insert(pC).second;
      if (wasInserted) {
	if (pspFNearby) {
	  for (int iCF = pC->getNumFaces() - 1; iCF >= 0; iCF--) {
	    pspFNearby->insert(pC->getFace(iCF));
	  }
	}
	for (int iV = pC->getNumVerts() - 1; iV >= 0; iV--) {
	  Vert *pVAdd = pC->getVert(iV);
	  if (!pVAdd->isDeleted())
	    spVNeigh.insert(pVAdd);
	}
      }
    }
    else if (pspBFInc) {
      pspBFInc->insert(dynamic_cast<BFace*>(pC));
    }
    for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
      Vert *pVAdd = pF->getVert(iV);
      if (!pVAdd->isDeleted())
	spVNeigh.insert(pVAdd);
    }
  }
  spVNeigh.erase(const_cast<Vert*>(pVert));

  if (pspBFInc) {
    assert(pqBdryVert);
    *pqBdryVert = (!pspBFInc->empty());
  }
  //
  //
  //  Face *pF = pVert->getHintFace();
  //  assert(pF->isValid() && pF->hasVert(pVert));
  //
  //  spCInc.clear();
  //  spVNeigh.clear();
  //  if (pspBFInc) {
  //    pspBFInc->clear();
  //    assert(pqBdryVert != NULL);
  //    (*pqBdryVert) = false;
  //  }
  //  if (pspFNearby) pspFNearby->clear();
  //
  //  Cell *pC = pF->getLeftCell();
  //  if (!pC->isValid() || pC->isBdryCell()) {
  //    pC = pF->getOppositeCell(pC);
  //  }
  //  spCInc.insert(pC);
  //  std::queue<Cell*> qpC;
  //  qpC.push(pC);
  //
  //  // Keep adding cells incident on pVert until they're all there.
  //  // spCInc.iLength() increases as cells are added, but the process
  //  // is guaranteed to exit.  In the process, create a list of pointers
  //  // to all the neighboring verts of pVert and all the faces incident on
  //  // pVert. The latter list will have to be trimmed later by removing faces
  //  // opposite to pVert.
  //  while (!qpC.empty()) {
  //    pC = qpC.front();
  //    qpC.pop();
  //    if(pC->isDeleted()) continue;
  //    // Add the parts of this cell to the appropriate lists.
  //    for (int iV = pC->getNumVerts() - 1; iV >= 0; iV--) {
  //	Vert *pVAdd = pC->getVert(iV);
  //	if (!pVAdd->isDeleted())
  //	  spVNeigh.insert(pVAdd);
  //    }
  //
  //    for (int iF = pC->getNumFaces() - 1; iF >= 0; iF--) {
  //      pF = pC->getFace(iF);
  //      if (pspFNearby) {
  //	pspFNearby->insert(pF);
  //      }
  //      if (pF->hasVert(pVert)) {
  //	// This Face has the right vertex, so the cell behind it does, too.
  //	Cell *pCCand = pF->getOppositeCell(pC);
  //	if (!pCCand->isValid()) continue;
  //	assert(pCCand->hasVert(pVert));
  //	if (pCCand->isBdryCell() && pspBFInc) {
  //	  pspBFInc->insert(dynamic_cast<BFace*>(pCCand));
  //	}
  //	// If BFaces  get added to this list, then we can walk right
  //	// across internal bdry faces.
  //	if (spCInc.count(pCCand) == 0) {
  //	  // The insertion is safe anyway, but the queue push is not!
  //	  spCInc.insert(pCCand);
  //	  qpC.push(pCCand);
  //	}
  //      }
  //    } // Done checking all faces for this cell
  //  } // Done checking all cells
  //
  //  // Clean up vertex list
  //  spVNeigh.erase(const_cast<Vert*>(pVert));
  //
  //  // Clean up cell list by removing any BFaces from it
  //  std::set<Cell*>::iterator iter;
  //  for (iter = spCInc.begin(); iter != spCInc.end(); ) {
  //    pC = *iter;
  //    if (pC->getType() == Cell::eBdryEdge ||
  //	pC->getType() == Cell::eIntBdryEdge ||
  //	pC->getType() == Cell::eTriBFace ||
  //	pC->getType() == Cell::eQuadBFace ||
  //	pC->getType() == Cell::eIntTriBFace ||
  //	pC->getType() == Cell::eIntQuadBFace) {
  //      std::set<Cell*>::iterator iterTmp = ++iter;
  //      spCInc.erase(iterTmp);
  //    }
  //    else ++iter;
  //  }
  //
  //  // Identify whether this is a bdry vert
  //  if (pspBFInc) {
  //    assert(pqBdryVert);
  //    *pqBdryVert = (!pspBFInc->empty());
  //  }
}

bool
Mesh::hasVert(const Vert* const pV) const
{

  int i, num_verts = getNumVerts();

  for (i = 0; i < num_verts; i++) {
    if (pV == getVert(i))
      return true;
  }

  return false;

}

//int Mesh::queueEncroachedBdryEntities(InsertionQueue& IQ) const
//{
//	// Make a set of all BFaces.
//	std::set<BFace*> spBF;
//	GR_index_t i;
//	for (i = 0; i < getNumTotalBdryFaces(); i++) {
//		spBF.insert(getBFace(i));
//	}
//	return queueEncroachedBdryEntities(IQ, spBF);
//}

int
Mesh::calcNumFacesNonDelaunay() const
{
  int retVal = 0;
  for (GR_index_t iF = 0; iF < getNumFaces(); iF++) {
    Face *pF = getFace(iF);
    if (!pF->isDeleted() && !pF->isLocallyDelaunay())
      retVal++;
  }
  return retVal;
}

//void Mesh::setLengthScale(const double dLen) {
//	m_LipschitzAlpha = 1;
//	m_resolution = 1;
//	for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
//		Vert *pV = getVert(iV);
//		pV->setLengthScale(dLen);
//	}
//}

void
Mesh::checkEdgeLengths() const
{
  int iNEdges = 0;
  double dMin = LARGE_DBL;
  double dMax = 0;
  double dSum = 0;

  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
    Vert *pV = getVert(iV);
    if (pV->isDeleted())
      continue;
    std::set<Vert*> spVNeigh;
    {
      std::set<Cell*> spCTmp;
      findNeighborhoodInfo(pV, spCTmp, spVNeigh);
    }
    double dLenScale;

    // Use the Vertex-stored length scale
    dLenScale = pV->getLengthScale();

    std::set<Vert*>::iterator iterV;
    for (iterV = spVNeigh.begin(); iterV != spVNeigh.end(); ++iterV) {
      Vert *pVNeigh = *iterV;
      double dDist = calcDistanceBetween(pVNeigh, pV);
      double dRelLen = dDist / dLenScale;
      dMax = max(dMax, dRelLen);
      dMin = min(dMin, dRelLen);
      dSum += dRelLen;
      iNEdges++;
    }
  }
  double dAve = dSum / iNEdges;
  logMessage(
      2, "Checked relative length scale for %d edges, counting duplicates.\n",
      iNEdges);
  logMessage(2, "Average: %f.  Min: %f.  Max: %f.\n", dAve, dMin, dMax);
}

void
Mesh::checkForITAPSBdryErrors()
{
  // This routine is a hack to make sure that all faces on the bdry have
  // a BFace attached to them.  Any faces that have a missing cell are
  // automatically assigned a BFace with a default BC instead.  If
  // called at the wrong time, this routine could easily break a mesh in
  // a way that ITAPS can't fix, so it should only be called from vPurge,
  // which in turn is called from the ITAPS implementation only when
  // allowing entity handles to change.
  for (GR_index_t iF = 0; iF < getNumFaces(); iF++) {
    Face *pF = getFace(iF);
    if (pF->isDeleted())
      continue;
    Cell *pCL = pF->getLeftCell();
    Cell *pCR = pF->getRightCell();
    if (!pCL->isValid() && !pCR->isValid()) {
      deleteFace(pF);
      continue;
    }
    if (!pCL->isValid() || !pCR->isValid()) {
      BFace *pBF = createBFace(pF);
      assert(pBF->doFullCheck());
      assert(pF->doFullCheck());
    }
  }
}

#ifdef HAVE_MPI
const GRUMMP::Ghost&
Mesh::getGhostEntry(Entity* const ent, bool& status)
{
  std::map<Entity*, GRUMMP::Ghost>::const_iterator iter = m_ghostMap.find(ent);
  status = (iter != m_ghostMap.end());
  if (!status) {
    // Create a new entry for this entity.
    addGhostEntry(ent, getPartID());
    iter = m_ghostMap.find(ent);
    status = (iter != m_ghostMap.end());
    assert(status);
  }
  status = status && iter->second.getLocalID() == getPartID();
  return iter->second;
}

bool
Mesh::addGhostEntry(Entity* const ent, const GRUMMP::PartGID& ownerID,
		    const GR_index_t nCopies)
{
  GRUMMP::Ghost ghost(getPartID(), nCopies);
  ghost.setOwnerID(ownerID);
  std::pair<std::map<Entity*, GRUMMP::Ghost>::iterator, bool> result =
      m_ghostMap.insert(std::make_pair(ent, ghost));
  return result.second;
}

void
Mesh::removeGhostEntry(Entity* const ent)
{
  std::map<Entity*, GRUMMP::Ghost>::iterator result = m_ghostMap.find(ent);
  if (result != m_ghostMap.end()) {
    m_ghostMap.erase(result);
  }
  if (ent->getEntStatus() == Entity::eGhosted) {
    ent->setEntStatus(Entity::eInterior);
  }
}

// For part bdry entities, remove all copies but the owner, or a stub if
// we know we're not the owner.
void
Mesh::resetGhostEntry(Entity* const ent)
{
  std::map<Entity*, GRUMMP::Ghost>::iterator result = m_ghostMap.find(ent);
  assert(result != m_ghostMap.end());
  // Reset the number of copies the entity knows about to the original number
  // it has when the only copies are entities shared on the part bdry.
  (result->second).reset();
}

bool
Mesh::addCopyToGhostEntry(Entity* const ent, const GRUMMP::EntityGID& remote,
			  const bool isBdry)
{
  bool status = true;
  GRUMMP::Ghost ghost(getGhostEntry(ent, status));
  if (status) {
    status = ghost.addCopy(remote, isBdry);
    m_ghostMap[ent] = ghost;
    assert(m_ghostMap[ent].getNumCopies() <= m_ghostMap[ent].getMaxNumCopies());
    assert(ghost.getNumCopies() <= ghost.getMaxNumCopies());
  }
  return status;
}
#endif

#ifndef HAVE_OPENMP
int
omp_in_parallel()
{
  return 0;
}
#endif
