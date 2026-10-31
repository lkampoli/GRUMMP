#include <map>

using std::map;
using std::pair;

#include "GR_VolMesh.h"

//@ Weed out empty entries in the connectivity tables
void
VolMesh::purgeFaces(map<Face*, Face*>* face_map)
// Get rid of all unused faces in the connectivity table.
{
  //   ECTriF.purge();
  //   ECQuadF.purge();

  if (getNumTriFaces() > 0) {
    TriFace *pTFFor = m_ECTriF.pTFirst();
    TriFace *pTFBack = m_ECTriF.pTLast();

    GR_index_t iNewSize = 0;
    assert(pTFFor != pTFBack);

    while (pTFFor != pTFBack) {

      while (pTFFor != pTFBack && !pTFFor->isDeleted()) {
//				if(face_map) face_map->insert( std::make_pair(pTFFor, pTFFor) );
	pTFFor = m_ECTriF.pTNext(pTFFor);
	iNewSize++;
      }

      while (pTFBack != pTFFor && pTFBack->isDeleted())
	pTFBack = m_ECTriF.pTPrev(pTFBack);

      if (pTFFor != pTFBack) {
	assert(pTFFor->isDeleted());
	if (face_map)
	  face_map->insert(std::make_pair(pTFBack, pTFFor));
	iNewSize++;
	(*pTFFor) = (*pTFBack);
	pTFBack->markAsDeleted();
	pTFFor = m_ECTriF.pTNext(pTFFor);
      }
      else if (!pTFFor->isDeleted()) {
	// This case only arises if there are no deleted faces to be purged.
//				if(face_map) face_map->insert( std::make_pair(pTFFor, pTFFor) );
	iNewSize++;
      }
    } // Loop over all faces
      // This resets the number of entries in use to iNewSize.
    m_ECTriF.reserve(iNewSize);
  } // Done with tris

  if (getNumQuadFaces() > 0) {
    QuadFace *pQFFor = m_ECQuadF.pTFirst();
    QuadFace *pQFBack = m_ECQuadF.pTLast();

    GR_index_t iNewSize = 0;
    assert(pQFFor != pQFBack);

    while (pQFFor != pQFBack) {
      while (pQFFor != pQFBack && !pQFFor->isDeleted()) {
//				if(face_map) face_map->insert( std::make_pair(pQFFor, pQFFor) );
	pQFFor = m_ECQuadF.pTNext(pQFFor);
	iNewSize++;
      }

      while (pQFBack != pQFFor && pQFBack->isDeleted())
	pQFBack = m_ECQuadF.pTPrev(pQFBack);

      if (pQFFor != pQFBack) {
	assert(pQFFor->isDeleted());
	if (face_map)
	  face_map->insert(std::make_pair(pQFBack, pQFFor));
	iNewSize++;
	(*pQFFor) = (*pQFBack);
	pQFBack->markAsDeleted();
	pQFFor = m_ECQuadF.pTNext(pQFFor);
      }
      else if (!pQFFor->isDeleted()) {
	// This case only arises if there are no deleted faces to be purged.
//				if(face_map) face_map->insert( std::make_pair(pQFFor, pQFFor) );
	iNewSize++;
      }
    } // Loop over all faces
      // This resets the number of entries in use to iNewSize.
    m_ECQuadF.reserve(iNewSize);
  } // Done with quads
}

GR_index_t
VolMesh::purgeCellRange(GR_index_t iFor, GR_index_t iBack,
			map<Cell*, Cell*>* cell_map)
{
  GR_index_t iNewSize = 0;
  // The case of exactly one valid cell
  if ((iFor == iBack) && getCell(iFor)->isValid()
      && !getCell(iFor)->isDeleted()) {
    //			if(cell_map) cell_map->insert( std::make_pair(getCell(iFor), getCell(iFor)) );
    iNewSize = 1;
    iFor = 1; // Necessary because of an assertion later.
  }
  while (iFor < iBack) {
    while (!getCell(iFor)->isDeleted() && iFor < iBack) {
      //				if(cell_map) cell_map->insert( std::make_pair(getCell(iFor), getCell(iFor)) );
      iFor++;
      iNewSize++;
    }
    while (getCell(iBack)->isDeleted() && iFor < iBack)
      iBack--;
    if (iFor != iBack) {
      Cell* pCFor = getCell(iFor);
      Cell* pCBack = getCell(iBack);
      assert(pCFor->isDeleted());
      assert(!pCBack->isDeleted());
      if (cell_map)
	cell_map->insert(std::make_pair(pCBack, pCFor));

      (*pCFor) = (*pCBack);
      pCBack->markAsDeleted();
      // #ifndef NDEBUG
      //	for (int i = 0; i < pCFor->iNumFaces(); i++)
      //	  assert(pCFor->pFFace(i)->iFullCheck());
      // #endif
      iNewSize++;
      iFor++;
    }
    else if (!getCell(iFor)->isDeleted()) {
      //				if(cell_map) cell_map->insert( std::make_pair(getCell(iFor), getCell(iFor)) );
      iNewSize++;
      iFor++;
    }
  }
  return iNewSize;
}

void
VolMesh::purgeCells(map<Cell*, Cell*>* cell_map)
// Get rid of all unused cells in the connectivity table.
{
  // Purge tetrahedra
  if (getNumTetCells() > 0) {
    GR_index_t iFor = 0, iBack = getNumTetCells() - 1;

    GR_index_t iNewSize = purgeCellRange(iFor, iBack, cell_map);
    // Reset the number of tet cells in use
    m_ECTet.reserve(iNewSize);
  } // Tetrahedron iteration

  // Purge pyramids
  if (getNumPyrCells() > 0) {
    GR_index_t iFor = getNumTetCells(), iBack = iFor + getNumPyrCells() - 1;

    GR_index_t iNewSize = purgeCellRange(iFor, iBack, cell_map);
    // Reset the number of pyramid cells in use
    m_ECPyr.reserve(iNewSize);
  } // Pyramid iteration

  // Purge prisms
  if (getNumPrismCells() > 0) {
    GR_index_t iFor = getNumTetCells() + getNumPyrCells(), iBack = iFor
	+ getNumPrismCells() - 1;

    GR_index_t iNewSize = purgeCellRange(iFor, iBack, cell_map);
    // Reset the number of prism cells in use
    m_ECPrism.reserve(iNewSize);
  } // Prism iteration

  // Purge hexahedra
  if (getNumHexCells() > 0) {
    GR_index_t iFor = getNumTetCells() + getNumPyrCells() + getNumPrismCells(),
	iBack = iFor + getNumHexCells() - 1;

    GR_index_t iNewSize = purgeCellRange(iFor, iBack, cell_map);
    m_ECHex.reserve(iNewSize);
  } // Hexahedron iteration
}

void
VolMesh::purgeBdryFaces(map<BFace*, BFace*>* bface_map)
// Get rid of all unused cells in the connectivity table.
{
  {
    TriBFace *pTBFFor = m_ECTriBF.pTFirst();
    TriBFace *pTBFBack = m_ECTriBF.pTLast();

    GR_index_t iNewSize = 0;

    if (m_ECTriBF.size() != 0) {
      while (pTBFFor != pTBFBack) {
	while (!pTBFFor->isDeleted() && pTBFFor != pTBFBack) {
//				if(bface_map) bface_map->insert( std::make_pair(pTBFFor, pTBFFor) );
	  pTBFFor = m_ECTriBF.pTNext(pTBFFor);
	  iNewSize++;
	}

	while (pTBFBack->isDeleted() && pTBFFor != pTBFBack)
	  pTBFBack = m_ECTriBF.pTPrev(pTBFBack);

	if (pTBFFor != pTBFBack) {
	  assert(pTBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pTBFBack, pTBFFor));
	  iNewSize++;
	  (*pTBFFor) = (*pTBFBack);
	  pTBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pTBFFor->getNumFaces(); i++)
	    assert(pTBFFor->getFace(i)->doFullCheck());
#endif
	  pTBFFor = m_ECTriBF.pTNext(pTBFFor);
	}
	else if (!pTBFFor->isDeleted()) {
//				if(bface_map) bface_map->insert( std::make_pair(pTBFFor, pTBFFor) );
	  iNewSize++;
	}

      }
    }
    // Reset the number of tri bdry faces in use
    m_ECTriBF.reserve(iNewSize);
  } // Triangular boundary face iteration
  {
    IntTriBFace *pITBFFor = m_ECIntTriBF.pTFirst();
    IntTriBFace *pITBFBack = m_ECIntTriBF.pTLast();

    GR_index_t iNewSize = 0;

    if (m_ECIntTriBF.size() != 0) {
      while (pITBFFor != pITBFBack) {
	while (!pITBFFor->isDeleted() && pITBFFor != pITBFBack) {
//				if(bface_map) bface_map->insert( std::make_pair(pITBFFor, pITBFFor) );
	  pITBFFor = m_ECIntTriBF.pTNext(pITBFFor);
	  iNewSize++;
	}

	while (pITBFBack->isDeleted() && pITBFFor != pITBFBack)
	  pITBFBack = m_ECIntTriBF.pTPrev(pITBFBack);

	if (pITBFFor != pITBFBack) {
	  assert(pITBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pITBFBack, pITBFFor));
	  iNewSize++;
	  (*pITBFFor) = (*pITBFBack);
	  pITBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pITBFFor->getNumFaces(); i++)
	    assert(pITBFFor->getFace(i)->doFullCheck());
#endif
	  pITBFFor = m_ECIntTriBF.pTNext(pITBFFor);
	}
	else if (!pITBFFor->isDeleted()) {
//				if(bface_map) bface_map->insert( std::make_pair(pITBFFor, pITBFFor) );
	  iNewSize++;
	}
      }
    }
    // Reset the number of internal tri bdry faces in use
    m_ECIntTriBF.reserve(iNewSize);
  } // Internal triangular boundary face iteration
  {
    QuadBFace *pQBFFor = m_ECQuadBF.pTFirst();
    QuadBFace *pQBFBack = m_ECQuadBF.pTLast();

    int iNewSize = 0;

    if (m_ECQuadBF.size() != 0) {
      while (pQBFFor != pQBFBack) {
	while (!pQBFFor->isDeleted() && pQBFFor != pQBFBack) {
//				if(bface_map) bface_map->insert( std::make_pair(pQBFFor, pQBFFor) );
	  pQBFFor = m_ECQuadBF.pTNext(pQBFFor);
	  iNewSize++;
	}

	while (pQBFBack->isDeleted() && pQBFFor != pQBFBack)
	  pQBFBack = m_ECQuadBF.pTPrev(pQBFBack);

	if (pQBFFor != pQBFBack) {
	  assert(pQBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pQBFBack, pQBFFor));
	  iNewSize++;
	  (*pQBFFor) = (*pQBFBack);
	  pQBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pQBFFor->getNumFaces(); i++)
	    assert(pQBFFor->getFace(i)->doFullCheck());
#endif
	  pQBFFor = m_ECQuadBF.pTNext(pQBFFor);
	}
	else if (!pQBFFor->isDeleted()) {
//				if(bface_map) bface_map->insert( std::make_pair(pQBFFor, pQBFFor) );
	  iNewSize++;
	}
      }
    }
    // Reset the number of quad bdry faces in use
    m_ECQuadBF.reserve(iNewSize);
  } // Quadangular boundary face iteration
  {
    IntQuadBFace *pITBFFor = m_ECIntQuadBF.pTFirst();
    IntQuadBFace *pITBFBack = m_ECIntQuadBF.pTLast();

    GR_index_t iNewSize = 0;

    if (m_ECIntQuadBF.size() != 0) {
      while (pITBFFor != pITBFBack) {
	while (!pITBFFor->isDeleted() && pITBFFor != pITBFBack) {
//				if(bface_map) bface_map->insert( std::make_pair(pITBFFor, pITBFFor) );
	  pITBFFor = m_ECIntQuadBF.pTNext(pITBFFor);
	  iNewSize++;
	}

	while (pITBFBack->isDeleted() && pITBFFor != pITBFBack)
	  pITBFBack = m_ECIntQuadBF.pTPrev(pITBFBack);

	if (pITBFFor != pITBFBack) {
	  assert(pITBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pITBFBack, pITBFFor));
	  iNewSize++;
	  (*pITBFFor) = (*pITBFBack);
	  pITBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pITBFFor->getNumFaces(); i++)
	    assert(pITBFFor->getFace(i)->doFullCheck());
#endif
	  pITBFFor = m_ECIntQuadBF.pTNext(pITBFFor);
	}
	else if (!pITBFFor->isDeleted()) {
//				if(bface_map) bface_map->insert( std::make_pair(pITBFFor, pITBFFor) );
	  iNewSize++;
	}
      }
    }
    // Reset the number of internal quad bdry faces in use
    m_ECIntQuadBF.reserve(iNewSize);
  } // Internal quadrilateral boundary face iteration
}

void
VolMesh::purgeVerts(map<Vert*, Vert*>* vert_map)
// Renumber vertices, removing those that have been marked as deleted
{
  logMessage(2, "Purging verts\n");
  GR_index_t iDelCount = 0;
  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
    Vert *pV = getVert(iV);
    if (pV->isDeleted())
      iDelCount++;
  }
  if (iDelCount != 0) {
    setAllHintFaces();
    setVertFaceNeighbors();

    GR_index_t *aiVertInd = new GR_index_t[getNumVerts()];
    GR_index_t iVActive = 0;

    // Copy the vertex data down within the list and establish a lookup so
    // that the faces can figure out where their verts are now stored.
    for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
      Vert *pV = getVert(iV);
      if (pV->isDeleted() || (!pV->getHintFace()->isValid())) {
	// If vert is no longer part of the mesh, copy over it.
	pV->markAsDeleted(); // Yes, this may have zero effect, but it's
	// an inline setting of a single bit, so who
	// cares?
	aiVertInd[iV] = INT_MAX;
      }
      else {
	aiVertInd[iV] = iVActive;
	Vert *pVA = getVert(iVActive);
	// Copy the data from its old to its new location
	if (pVA != pV) {
	  if (vert_map)
	    vert_map->insert(std::make_pair(pV, pVA));

	  (*pVA) = (*pV);
	  // update the subsegments if we have them
	  if (hasSubsegs() && isSubsegVert(pV)) {
	    replaceSubsegVert(pV, pVA);
	  }
	  pV->markAsDeleted();
	}

	iVActive++;
      }
    }

    //     vClearVertFaceNeighbors();
    //     vClearAllHintFaces();
    //     for (GR_index_t iF = 0; iF < iNumFaces(); iF++) {
    //       Face *pF = pFFace(iF);
    //       if (pF->qDeleted()) continue;
    //       Vert *apVNew[4];
    //       for (int ii = pF->iNumVerts() - 1; ii >= 0; ii--)
    // 	apVNew[ii] = pVVert(aiVertInd[iVertIndex(pF->pVVert(ii))]);
    //       // Passing the extra arguments does no harm, as they are ignored.
    //       pF->vAssign(pF->pCCellLeft(), pF->pCCellRight(), apVNew[0],
    // 		  apVNew[1], apVNew[2], apVNew[3]);
    //     }
    delete[] aiVertInd;
    m_ECVerts.reserve(iVActive);
    setAllHintFaces();
    setVertFaceNeighbors();
  }
//	else if(vert_map) {
//		for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
//			Vert *pV = getVert(iV);
//			vert_map->insert( std::make_pair(pV, pV) );
//		}
//	}
  // If there are no verts deleted, then skip all the hard work.
  return;
}

void
VolMesh::purgeAllEntities(map<Vert*, Vert*>* vert_map,
			  map<Face*, Face*>* face_map,
			  map<Cell*, Cell*>* cell_map,
			  map<BFace*, BFace*>* bface_map)
{
  clearEventQueues();
#ifdef ITAPS
  checkForITAPSBdryErrors();
#endif
  purgeCells(cell_map);
  purgeFaces(face_map);
  purgeBdryFaces(bface_map);
  purgeVerts(vert_map);

}

typedef struct _vert_vert_ {
  int iV0, iV1;
  _vert_vert_(const int i0 = -1, const int i1 = -1) :
      iV0(i0), iV1(i1)
  {
  }
  bool
  operator<(const struct _vert_vert_& v2) const
  {
    return (iV0 < v2.iV0 || (iV0 == v2.iV0 && iV1 < v2.iV1));
  }
} VertVert;

typedef struct _renum_ {
  int iV, iDeg;
  _renum_(const int iVIn = -1, const int iDegIn = 0) :
      iV(iVIn), iDeg(iDegIn)
  {
  }
} Renum;

extern "C"
{
  static int
  qRenumComp3D(const void *pv0, const void *pv1)
  {
    const Renum * const pR0 = (const Renum *) pv0;
    const Renum * const pR1 = (const Renum *) pv1;
    return ((pR0->iDeg < pR1->iDeg) ? -1 : ((pR0->iDeg == pR1->iDeg) ? 0 : 1));
  }
}

//@ Reorder vertices in a 2D mesh using the reverse Cuthill-McKee algorithm
void
VolMesh::reorderVerts_RCM(void)
{
  logMessage(2, "Renumbering vertices...");
  assert(isSimplicial());
  int iVert, iFace, iNVerts = getNumVerts(), iNFaces = getNumFaces();
  //@@ Construct the vert-vert connectivity
  // Make a list containing each edge (once in each direction) and sort
  int iNEdges = 0;
  VertVert *aVV = NULL;
  {
    std::set<VertVert> sVV;
    for (iFace = iNFaces - 1; iFace >= 0; iFace--) {
      Face *pF = getFace(iFace);
      int iV0 = getVertIndex(pF->getVert(0));
      int iV1 = getVertIndex(pF->getVert(1));
      int iV2 = getVertIndex(pF->getVert(2));
      sVV.insert(VertVert(iV1, iV0));
      sVV.insert(VertVert(iV0, iV1));
      sVV.insert(VertVert(iV1, iV2));
      sVV.insert(VertVert(iV2, iV1));
      sVV.insert(VertVert(iV2, iV0));
      sVV.insert(VertVert(iV0, iV2));
    }
    // This is actually twice the number of unique edges, because both
    // directions appear here.  iNEdges shouldn't be roughly the same as
    // iNFaces.
    iNEdges = sVV.size();
    if (iNEdges > 0) {
      aVV = new VertVert[iNEdges];
      VertVert *pVV = copy(sVV.begin(), sVV.end(), aVV);
      assert(pVV - aVV == iNEdges);
    }
  }

  // Tabulate starting locations and degree for each vertex.
  std::vector<int> viStart(iNVerts + 1), viDegree(iNVerts);
  std::vector<bool> vqUsed(iNVerts);
  int iStart = 0, iNAttachedVerts = 0;
  for (iVert = 0; iVert < iNVerts; iVert++) {
    viStart[iVert] = iStart;
    for (; iStart < iNEdges && aVV[iStart].iV0 == iVert; iStart++) {
    }
    if (iStart != viStart[iVert])
      iNAttachedVerts++;
    viDegree[iVert] = iStart - viStart[iVert];
    vqUsed[iVert] = false;
  }
  viStart[iNVerts] = iNEdges;

  // Set up space to compute re-ordering
  Renum *aRReorder = new Renum[iNVerts];

  int iNumRenumbered = 0;
  while (iNumRenumbered < iNAttachedVerts) {
    //@@ Initialize the front to be a single unrenumbered vertex with lowest degree
    // Find lowest degree vertex
    int iMinDeg = 100, iMinVert = -1;
    for (iVert = 0; iVert < iNAttachedVerts; iVert++) {
      int iDeg = viStart[iVert + 1] - viStart[iVert];
      if (iDeg < iMinDeg && iDeg > 0 && !vqUsed[iVert]) {
	iMinDeg = iDeg;
	iMinVert = iVert;
	if (iMinDeg == 2)
	  break;	// Can never be lower than this, so why go on?
      }
    }
    assert(iMinVert != -1);

    // Initialize
    vqUsed[iMinVert] = true;
    aRReorder[iNumRenumbered] = Renum(iMinVert, iMinDeg);
    iNumRenumbered++;
    int iEnd = iNumRenumbered, iBegin = iNumRenumbered - 1;

    //@@ Propagate across the mesh
    do {
      //@@@ Add all the next layer of vertices
      int iEnt;
      for (iEnt = iBegin; iEnt < iEnd; iEnt++) {
	int iV = aRReorder[iEnt].iV;
	for (int iNeigh = viStart[iV]; iNeigh < viStart[iV + 1]; iNeigh++) {
	  // Grab a candidate
	  assert(aVV[iNeigh].iV0 == iV);
	  int iCand = aVV[iNeigh].iV1;
	  if (!vqUsed[iCand]) {
	    // Add it if it isn't already in the list
	    aRReorder[iNumRenumbered].iV = iCand;
	    aRReorder[iNumRenumbered].iDeg = viDegree[iCand];
	    vqUsed[iCand] = true;
	    iNumRenumbered++;
	  }
	}			// Loop over neighbors of a given vertex
      }				// Loop over all verts in the most recent layer

				//@@@ Sort the new layer in order of increasing degree
      iBegin = iEnd;
      iEnd = iNumRenumbered;
      qsort(&(aRReorder[iBegin]), static_cast<GR_index_t>(iEnd - iBegin),
	    sizeof(Renum), qRenumComp3D);
    }
    while (iEnd != iBegin);
  }

  for (iVert = 0; iNumRenumbered != iNVerts; iVert++) {
    if (!vqUsed[iVert]) {
      aRReorder[iNumRenumbered].iV = iVert;
      vqUsed[iVert] = true;
      iNumRenumbered++;
    }
  }

  //@@ Reverse the ordering
  for (iVert = 0; iVert < iNVerts / 2; iVert++) {
    int iTemp = aRReorder[iVert].iV;
    aRReorder[iVert].iV = aRReorder[iNVerts - iVert - 1].iV;
    aRReorder[iNVerts - iVert - 1].iV = iTemp;
  }

  //@@ Move vertex data
  // Create a new copy
  EntContainer<Vert> ECTemp;
  ECTemp.reserve(iNVerts);
  for (iVert = 0; iVert < iNVerts; iVert++) {
    // Copy both coordinates and vertex tags
    ECTemp[iVert] = (*getVert(iVert));
  }
  for (iVert = 0; iVert < iNVerts; iVert++) {
    // Copy back into re-ordered locations
    Vert *pVSource = &(ECTemp[aRReorder[iVert].iV]);
    Vert *pVTarget = getVert(iVert);
    (*pVTarget) = (*pVSource);
  }
  delete[] aVV;
  delete[] aRReorder;

  //@@ Reset vertex hints
  setAllHintFaces();
  setVertFaceNeighbors();

  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Vertex re-ordering for 3D meshes");
  }
  logMessage(2, "done.\n");
}

void
VolMesh::reorderCells()
// Reorder cells according to vertex indices, so that cells with the
// lowest sum of vertex indices are first in order.  The hardest part
// is actually to get the faces and cells reconnected afterwards.
{
  logMessage(2, "Renumbering cells...");
  assert(isSimplicial());
  std::multimap<GR_index_t, GR_index_t> mOrderedCells;
  GR_index_t iC;
  for (iC = 0; iC < getNumCells(); iC++) {
    Cell *pC = getCell(iC);
    GR_index_t iSum = 0;
    int iNV = pC->getNumVerts();
    for (int iV = 0; iV < iNV; iV++) {
      iSum += getVertIndex(pC->getVert(iV));
    }
    mOrderedCells.insert(pair<GR_index_t, GR_index_t>(iSum, iC));
  }
  // Now we have a list of all the cells, sorted in order by smallest
  // sum of vertex indices!  Use this to create a map that tells
  // which old cell goes into which slot in the new numbering.
  std::vector<GR_index_t> aiCellMapping(getNumCells());
  std::multimap<GR_index_t, GR_index_t>::iterator iter;
  iC = 0;
  for (iter = mOrderedCells.begin(); iter != mOrderedCells.end(); ++iter) {
    logMessage(4, "Old cell: %u  New cell: %u  Sum: %u\n", iter->second, iC,
	       iter->first);
    aiCellMapping[iC] = iter->second;
    iC++;
  }
  assert(iC == getNumCells());

  // Finally, update the cell information.  Fundamentally, re-ordering
  // the cells results in one or more circular permutations of cells
  // (A->B->Q->L->A, etc).  So follow each of these.
  std::vector<bool> aqUpdated(getNumCells());
  for (iC = 0; iC < getNumCells(); iC++) {
    aqUpdated[iC] = false;
  }

  GR_index_t iNumUpdated = 0;
  GR_index_t iFirst = 0;
  while (iNumUpdated != getNumCells()) {
    while (aqUpdated[iFirst])
      iFirst++;
    assert(iFirst < getNumCells());
    if (aiCellMapping[iFirst] == iFirst) {
      aqUpdated[iFirst] = true;
      iNumUpdated++;
    }
    else {
      assert(getCell(iFirst)->getType() == Cell::eTet);
      TetCell TC;
      TC = *(dynamic_cast<TetCell*>(getCell(iFirst)));
      GR_index_t iNext = iFirst;
      while (aiCellMapping[iNext] != iFirst) {
	(*getCell(iNext)) = (*getCell(aiCellMapping[iNext]));
	aqUpdated[iNext] = true;
	iNumUpdated++;
	iNext = aiCellMapping[iNext];
      }
      (*getCell(iNext)) = TC;
      aqUpdated[iNext] = true;
      iNumUpdated++;
    }
  }
  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Cell re-ordering for 3D meshes");
  }
  logMessage(2, "done.\n");
}

void
VolMesh::reorderFaces()
// Reorder faces according to vertex indices, so that faces with the
// lowest sum of vertex indices are first in order.  The hardest part
// is actually to get the faces and cells reconnected afterwards.
{
  logMessage(2, "Renumbering faces...");
  assert(isSimplicial());
  std::multimap<GR_index_t, GR_index_t> mOrderedFaces;
  GR_index_t iF, iNF = getNumFaces();
  for (iF = 0; iF < iNF; iF++) {
    Face *pF = getFace(iF);

    GR_index_t iSum = getVertIndex(pF->getVert(0))
	+ getVertIndex(pF->getVert(1));
    mOrderedFaces.insert(pair<GR_index_t, GR_index_t>(iSum, iF));
  }
  // Now we have a list of all the faces, sorted in order by smallest
  // average of vertex indices!  Use this to create a map of which -old-
  // face becomes which -new- face.
  std::vector<GR_index_t> aiFaceMapping(iNF);
  std::multimap<GR_index_t, GR_index_t>::iterator iter;
  iF = 0;
  for (iter = mOrderedFaces.begin(); iter != mOrderedFaces.end(); ++iter) {
    logMessage(4, "Old face: %u  New face: %u  Sum: %u\n", iter->second, iF,
	       iter->first);
    aiFaceMapping[iF] = iter->second;
    iF++;
  }
  assert(iF == iNF);

  // Finally, update the cell information.  Fundamentally, re-ordering
  // the cells results in one or more circular permutations of cells
  // (A->B->Q->L->A, etc).  So follow each of these.
  std::vector<bool> aqUpdated(iNF);
  for (iF = 0; iF < iNF; iF++) {
    aqUpdated[iF] = false;
  }

  GR_index_t iNumUpdated = 0;
  GR_index_t iFirst = 0;
  while (iNumUpdated != iNF) {
    while (aqUpdated[iFirst])
      iFirst++;
    assert(iFirst < iNF);
    if (aiFaceMapping[iFirst] == iFirst) {
      aqUpdated[iFirst] = true;
      iNumUpdated++;
    }
    else {
      assert(getFace(iFirst)->getType() == Face::eTriFace);
      TriFace TF;
      TF = *(dynamic_cast<TriFace*>(getFace(iFirst)));
      GR_index_t iNext = iFirst;
      while (aiFaceMapping[iNext] != iFirst) {
	(*getFace(iNext)) = (*getFace(aiFaceMapping[iNext]));
	aqUpdated[iNext] = true;
	iNumUpdated++;
	iNext = aiFaceMapping[iNext];
      }
      (*getFace(iNext)) = TF;
      aqUpdated[iNext] = true;
      iNumUpdated++;
    }
  }

  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Face re-ordering for 3D meshes");
  }
  logMessage(2, "done.\n");
}

