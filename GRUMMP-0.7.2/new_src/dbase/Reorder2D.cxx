#include <map>

#include "GR_Mesh2D.h"

using std::map;
using std::multimap;
using std::pair;

//@ Purge unused entries
void
Mesh2D::purgeAllEntities(map<Vert*, Vert*>* vert_map,
			 map<Face*, Face*>* face_map,
			 map<Cell*, Cell*>* cell_map,
			 map<BFace*, BFace*>* bface_map)
{
  clearEventQueues();
#ifdef ITAPS
  checkForITAPSBdryErrors();
#endif
  purgeCells(cell_map);
  purgeBdryFaces(bface_map);
  purgeFaces(face_map);
  purgeVerts(vert_map);

  if (pCellQueue != NULL)
    pCellQueue->vReset();
}

void
Mesh2D::purgeVerts(map<Vert*, Vert*>* vert_map)
// Renumber vertices, removing those that have been marked as deleted
{
  GR_index_t *aiVertInd = new GR_index_t[iNumVerts()];
  GR_index_t iV, iVActive;

  setAllHintFaces();
  setVertFaceNeighbors();
  iVActive = 0;
  // Copy the vertex data down within the list and establish a lookup so
  // that the faces can figure out where their verts are now stored.
  for (iV = 0; iV < iNumVerts(); iV++) {
    Vert *pV = getVert(iV);
    if (pV->isDeleted()) {  // If vert is deleted, copy over it.
      aiVertInd[iV] = GR_index_t(INT_MAX) * 2;
    }
    else {
      aiVertInd[iV] = iVActive;
      Vert *pVA = getVert(iVActive);
      // Copy the data from its old to its new location
      if (pVA != pV) {
	if (vert_map)
	  vert_map->insert(std::make_pair(pV, pVA));

	(*pVA) = (*pV);
	pV->markAsDeleted();
      }
      iVActive++;
    }
  }

  // This should happen automatically now.
//   // Change the vertices for all the faces.
//   vClearVertFaceNeighbors();
//   vClearAllHintFaces();
//   for (int iF = iNumFaces() - 1; iF >= 0; iF--) {
//     Face *pF = pFFace(iF);
//     if (pF->qDeleted()) continue;
//     Vert *apVNew[4];
//     for (int ii = pF->iNumVerts() - 1; ii >= 0; ii--) {
//       GR_index_t tmp = aiVertInd[iVertIndex(pF->pVVert(ii))];

//       apVNew[ii] = pVVert(
//     // Passing the extra arguments does no harm, as they are ignored.
//     pF->vSetVerts(apVNew[0], apVNew[1], apVNew[2], apVNew[3]);
//   }
  setAllHintFaces();
  setVertFaceNeighbors();
//  assert(isValid());

  // Now change the influencing verts, where they occur, to match the
  // new vertex numbering.

  for (iV = 0; iV < iVActive; iV++) {
    Vert *pV = getVert(iV);
    assert(pV->isValid() && !pV->isDeleted());
    if (pV->getLengthScale() < 0) {
      // Don't bother copying data if no length scale exists yet.
      pV->clearLengthScale();
      continue;
    }
  }
  delete[] aiVertInd;

  setupVerts(iVActive);
//  for (iV = 0; iV < iVActive; iV++)					// these are the codes Reza added solely for debugging
//  {
//	  Vert *pV = getVert(iV);
//	  std::cout<<pV<<"\t"<<pV->getNumFaces()<<std::endl;
//  }
}

void
Mesh2D::purgeFaces(map<Face*, Face*>* face_map)
// Get rid of all unused faces in the connectivity table.
{
  EdgeFace *pFFor = m_ECEdgeF.pTFirst();
  EdgeFace *pFBack = m_ECEdgeF.pTLast();

  int iNewSize = 0;
  // If one's valid, they both are, but this keeps Klocwork happy...
  if (pFFor != pFInvalidFace && pFBack != pFInvalidFace) {
    // The case of exactly one valid face
    if (pFFor == pFBack && !pFFor->isDeleted()) {
//      if(face_map) face_map->insert( std::make_pair(pFFor, pFFor) );
      iNewSize++;
    }

    while (pFFor != pFBack) {
      while (!pFFor->isDeleted() && pFFor != pFBack) {
//	if(face_map) face_map->insert( std::make_pair(pFFor, pFFor) );
	pFFor = m_ECEdgeF.pTNext(pFFor);
	iNewSize++;
      }

      while (pFBack->isDeleted() && pFFor != pFBack)
	pFBack = m_ECEdgeF.pTPrev(pFBack);

      if (pFFor != pFBack) {
	assert(pFFor->isDeleted());
	if (face_map)
	  face_map->insert(std::make_pair(pFBack, pFFor));
	iNewSize++;
	(*pFFor) = (*pFBack);
	assert(!pFFor->isDeleted());
	pFBack->markAsDeleted();
	pFFor = m_ECEdgeF.pTNext(pFFor);
      }
      // Allow for the case where nothing has been deleted
      else if (!pFFor->isDeleted()) {
//	if(face_map) face_map->insert( std::make_pair(pFFor, pFFor) );
	iNewSize++;
      }
    }				// Loop over all faces
  }

  // This resets the number of entries in use to iNewSize.
  m_ECEdgeF.reserve(iNewSize);
  resetVertexConnectivity();

}

void
Mesh2D::purgeCells(map<Cell*, Cell*>* cell_map)
// Get rid of all unused cells in the connectivity table.
{
  {
    TriCell *pCFor = m_ECTri.pTFirst();
    TriCell *pCBack = m_ECTri.pTLast();
    int iNewSize = 0;

    // If one's valid, they both are, but this keeps Klocwork happy...
    if (pCFor != pCInvalidCell && pCBack != pCInvalidCell) {
      // The case of exactly one valid cell
      if (pCFor == pCBack && pCFor != pCInvalidCell && !pCFor->isDeleted()) {
	iNewSize = 1;
//	if(cell_map) cell_map->insert( std::make_pair(pCFor, pCFor) );
      }

      while (pCFor != pCBack) {
	while (!pCFor->isDeleted() && pCFor != pCBack) {
//	  if(cell_map) cell_map->insert( std::make_pair(pCFor, pCFor) );
	  pCFor = m_ECTri.pTNext(pCFor);
	  iNewSize++;
	}

	while (pCBack->isDeleted() && pCFor != pCBack)
	  pCBack = m_ECTri.pTPrev(pCBack);

	if (pCFor != pCBack) {
	  assert(pCFor->isDeleted());
	  if (cell_map)
	    cell_map->insert(std::make_pair(pCBack, pCFor));
	  iNewSize++;
	  (*pCFor) = (*pCBack);
	  pCBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pCFor->getNumFaces(); i++) {
	    assert(pCFor->getFace(i)->doFullCheck());
	  }
#endif
	  pCFor = m_ECTri.pTNext(pCFor);
	}
	// Allow for the case where nothing has been deleted
	else if (!pCFor->isDeleted()) {
//	  if(cell_map) cell_map->insert( std::make_pair(pCFor, pCFor) );
	  iNewSize++;
	}
      }				// Loop over triangles
      // Reset the number of tri cells in use
    }
    m_ECTri.reserve(iNewSize);
  }				// Contex for loop over triangles

  {
    QuadCell *pCFor = m_ECQuad.pTFirst();
    QuadCell *pCBack = m_ECQuad.pTLast();
    int iNewSize = 0;

    // If one's valid, they both are, but this keeps Klocwork happy...
    if (pCFor != pCInvalidCell && pCBack != pCInvalidCell) {
      // The case of exactly one valid cell
      if (pCFor == pCBack && pCFor != pCInvalidCell && !pCFor->isDeleted())
	iNewSize++;
      while (pCFor != pCBack) {
	while (!pCFor->isDeleted() && pCFor != pCBack) {
	  pCFor = m_ECQuad.pTNext(pCFor);
	  iNewSize++;
	}

	while (pCBack->isDeleted() && pCFor != pCBack)
	  pCBack = m_ECQuad.pTPrev(pCBack);

	if (pCFor != pCBack) {
	  assert(pCFor->isDeleted());
	  iNewSize++;

	  (*pCFor) = (*pCBack);
	  pCBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pCFor->getNumFaces(); i++)
	    assert(pCFor->getFace(i)->doFullCheck());
#endif
	  pCFor = m_ECQuad.pTNext(pCFor);
	}
	// Allow for the case where nothing has been deleted
	else if (!pCFor->isDeleted())
	  iNewSize++;
      }				// Loop over quads
      // Reset the number of quad cells in use
    }
    m_ECQuad.reserve(iNewSize);
  }				// Context for loop over quads
}

void
Mesh2D::purgeBdryFaces(map<BFace*, BFace*>* bface_map)
{

  // Get rid of all unused bdry faces in the connectivity table.
  // First the "real" bdry edges
  {
    BdryEdge *pEBFFor = m_ECEdgeBF.pTFirst();
    BdryEdge *pEBFBack = m_ECEdgeBF.pTLast();

    int iNewSize = 0;

    // If one's valid, they both are, but this keeps Klocwork happy...
    if (pEBFFor != pBFInvalidBFace && pEBFBack != pBFInvalidBFace) {
      while (pEBFFor != pEBFBack) {

	while (!pEBFFor->isDeleted() && pEBFFor != pEBFBack) {
//	  if(bface_map) bface_map->insert( std::make_pair(pEBFFor, pEBFFor) );
	  pEBFFor = m_ECEdgeBF.pTNext(pEBFFor);
	  iNewSize++;
	}

	while (pEBFBack->isDeleted() && pEBFFor != pEBFBack)
	  pEBFBack = m_ECEdgeBF.pTPrev(pEBFBack);

	if (pEBFFor != pEBFBack) {
	  assert(pEBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pEBFBack, pEBFFor));
	  iNewSize++;
	  (*pEBFFor) = (*pEBFBack);
	  pEBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pEBFFor->getNumFaces(); i++) {
	    assert(pEBFFor->getFace(i)->doFullCheck());
	  }
#endif
	  pEBFFor = m_ECEdgeBF.pTNext(pEBFFor);
	}
	else if (!pEBFFor->isDeleted()) {
//	  if(bface_map) bface_map->insert( std::make_pair(pEBFFor, pEBFFor) );
	  iNewSize++;
	}
      }
    }
    //Set the number of bdry faces in use
    m_ECEdgeBF.reserve(iNewSize);
  }

  // Now the internal bdry faces
  {
    int iNewSize = 0;
    IntBdryEdge *pIEBFFor = m_ECIntEdgeBF.pTFirst();
    IntBdryEdge *pIEBFBack = m_ECIntEdgeBF.pTLast();

    // If one's valid, they both are, but this keeps Klocwork happy...
    if (pIEBFFor != pBFInvalidBFace && pIEBFBack != pBFInvalidBFace) {
      // The case of exactly one valid internal bdry edge
      if (pIEBFFor == pIEBFBack && !pIEBFFor->isDeleted()) {
//	if(bface_map) bface_map->insert( std::make_pair(pIEBFFor, pIEBFFor) );
	iNewSize = 1;
      }

      while (pIEBFFor != pIEBFBack) {
	while (!pIEBFFor->isDeleted() && pIEBFFor != pIEBFBack) {
//	  if(bface_map) bface_map->insert( std::make_pair(pIEBFFor, pIEBFFor) );
	  pIEBFFor = m_ECIntEdgeBF.pTNext(pIEBFFor);
	  iNewSize++;
	}

	while (pIEBFBack->isDeleted() && pIEBFFor != pIEBFBack)
	  pIEBFBack = m_ECIntEdgeBF.pTPrev(pIEBFBack);

	if (pIEBFFor != pIEBFBack) {
	  assert(pIEBFFor->isDeleted());
	  if (bface_map)
	    bface_map->insert(std::make_pair(pIEBFBack, pIEBFFor));
	  iNewSize++;
	  (*pIEBFFor) = (*pIEBFBack);
	  pIEBFBack->markAsDeleted();
#ifndef NDEBUG
	  for (int i = 0; i < pIEBFFor->getNumFaces(); i++)
	    assert(pIEBFFor->getFace(i)->doFullCheck());
#endif
	  pIEBFFor = m_ECIntEdgeBF.pTNext(pIEBFFor);
	}
	else if (!pIEBFFor->isDeleted()) {
//	  if(bface_map) bface_map->insert( std::make_pair(pIEBFFor, pIEBFFor) );
	  iNewSize++;
	}
      }
    }
    // Reset the number of bdry faces in use
    m_ECIntEdgeBF.reserve(iNewSize);
  }
}

typedef struct _vert_vert_ {
  int iV0, iV1;
  _vert_vert_(const int i0 = -1, const int i1 = -1) :
      iV0(i0), iV1(i1)
  {
  }
} VertVert;

extern "C"
{
  static int
  qVVComp(const void *pv0, const void *pv1)
  {
    const VertVert * const pVV0 = (const VertVert *) pv0;
    const VertVert * const pVV1 = (const VertVert *) pv1;
    return
	(pVV0->iV0 < pVV1->iV0
	    || (pVV0->iV0 == pVV1->iV0 && pVV0->iV1 < pVV1->iV1)) ? -1 : 1;
  }
}

extern "C"
{
  static int
  qRenumComp2D(const void *pv0, const void *pv1)
  {
    const Renum * const pR0 = (const Renum *) pv0;
    const Renum * const pR1 = (const Renum *) pv1;
    return ((pR0->iDeg < pR1->iDeg) ? -1 : ((pR0->iDeg == pR1->iDeg) ? 0 : 1));
  }
}

//@ Reorder vertices in a 2D mesh using the reverse Cuthill-McKee algorithm
void
Mesh2D::reorderVerts_RCM(void)
{
  logMessage(2, "Renumbering vertices...");
  int iVert, iFace, iNVerts = iNumVerts(), iNFaces = getNumFaces();
  //@@ Construct the vert-vert connectivity
  // Make a list containing each edge (once in each direction) and sort
  VertVert *aVV = new VertVert[2 * iNFaces];
  for (iFace = iNFaces - 1; iFace >= 0; iFace--) {
    Face *pF = getFace(iFace);
    int iV0 = iVertIndex(pF->getVert(0));
    int iV1 = iVertIndex(pF->getVert(1));
    aVV[2 * iFace + 1] = VertVert(iV0, iV1);
    aVV[2 * iFace] = VertVert(iV1, iV0);
  }

  qsort(aVV, 2 * iNFaces, sizeof(VertVert), qVVComp);
  // Tabulate starting locations and degree for each vertex.
  int *aiStart = new int[iNVerts + 1];
  int *aiDegree = new int[iNVerts];
  bool *aqUsed = new bool[iNVerts];
  int iStart = 0, iNAttachedVerts = 0;
  for (iVert = 0; iVert < iNVerts; iVert++) {
    aiStart[iVert] = iStart;
    for (; (iStart < 2 * iNFaces) && (aVV[iStart].iV0 == iVert); iStart++) {
    }
    if (iStart != aiStart[iVert])
      iNAttachedVerts++;
    aiDegree[iVert] = iStart - aiStart[iVert];
    aqUsed[iVert] = false;
  }

  aiStart[iNVerts] = 2 * iNFaces;

  // Set up space to compute re-ordering

  // Making array of structs.. don't want this
  Renum *aRReorder = new Renum[iNVerts];

  int iNumRenumbered = 0;
  while (iNumRenumbered < iNAttachedVerts) {
    //@@ Initialize the front to be a single unrenumbered vertex with lowest degree
    // Find lowest degree vertex
    int iMinDeg = 100, iMinVert = -1;
    for (iVert = 0; iVert < iNAttachedVerts; iVert++) {
      int iDeg = aiStart[iVert + 1] - aiStart[iVert];
      if (iDeg < iMinDeg && iDeg > 0 && !aqUsed[iVert]) {
	iMinDeg = iDeg;
	iMinVert = iVert;
	if (iMinDeg == 2)
	  break;	// Can never be lower than this, so why go on?
      }
    }
    assert(iMinVert != -1);

    // Initialize
    aqUsed[iMinVert] = true;
    aRReorder[iNumRenumbered] = Renum(iMinVert, iMinDeg);
    iNumRenumbered++;
    int iEnd = iNumRenumbered, iBegin = iNumRenumbered - 1;

    //@@ Propagate across the mesh
    do {
      //@@@ Add all the next layer of vertices
      int iEnt;
      for (iEnt = iBegin; iEnt < iEnd; iEnt++) {
	int iV = aRReorder[iEnt].iV;
	for (int iNeigh = aiStart[iV]; iNeigh < aiStart[iV + 1]; iNeigh++) {
	  // Grab a candidate
	  assert(aVV[iNeigh].iV0 == iV);
	  int iCand = aVV[iNeigh].iV1;
	  if (!aqUsed[iCand]) {
	    // Add it if it isn't already in the list
	    aRReorder[iNumRenumbered].iV = iCand;
	    aRReorder[iNumRenumbered].iDeg = aiDegree[iCand];
	    aqUsed[iCand] = true;
	    iNumRenumbered++;
	  }
	}			// Loop over neighbors of a given vertex
      }				// Loop over all verts in the most recent layer

				//@@@ Sort the new layer in order of increasing degree
      iBegin = iEnd;
      iEnd = iNumRenumbered;
      qsort(&(aRReorder[iBegin]), static_cast<GR_index_t>(iEnd - iBegin),
	    sizeof(Renum), qRenumComp2D);
    }
    while (iEnd != iBegin);
  }
  delete[] aVV;
  delete[] aiStart;
  delete[] aiDegree;
  delete[] aqUsed;

  for (iVert = 0; iNumRenumbered != iNVerts; iVert++) {
    if (!aqUsed[iVert]) {
      aRReorder[iNumRenumbered].iV = iVert;
      aqUsed[iVert] = true;
      iNumRenumbered++;
    }
  }

  GR_index_t * aRReorderInt = new GR_index_t[iNVerts];

  for (iVert = 0; iVert < iNVerts; iVert++) {

    aRReorderInt[iVert] = aRReorder[iNVerts - iVert - 1].iV;

  }

  //@@ Move vertex data
  // If we copy this to a new array and copy it back, all the pointers from
  // faces get updated automagically (well, actually, there's explicit code in
  // Vert::operator=).
  reorderVertices(aRReorderInt);

  delete[] aRReorder;
  delete[] aRReorderInt;

  //@@ Reset vertex hints
  setAllHintFaces();
  setVertFaceNeighbors();

  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Vertex re-ordering for 2D meshes");
  }
  logMessage(2, "done.\n");
}

void
Mesh2D::reorderCells()
// Reorder cells according to vertex indices, so that cells with the
// lowest sum of vertex indices are first in order.  The hardest part
// is actually to get the faces and cells reconnected afterwards.
{
  logMessage(2, "Renumbering cells...");
  assert(isSimplicial());
  multimap<GR_index_t, GR_index_t> mOrderedCells;
  GR_index_t iC;
  for (iC = 0; iC < getNumCells(); iC++) {
    Cell *pC = getCell(iC);
    GR_index_t iSum = 0;
    GR_index_t iNV = pC->getNumVerts();
    for (GR_index_t iV = 0; iV < iNV; iV++) {
      iSum += iVertIndex(pC->getVert(iV));
    }
    mOrderedCells.insert(pair<GR_index_t, GR_index_t>(iSum, iC));
  }
  // Now we have a list of all the cells, sorted in order by smallest
  // sum of vertex indices!  Use this to create a map that tells
  // which old cell goes into which slot in the new numbering.  Note that this
  // ordering is not absolutely guaranteed to be unique.
  GR_index_t *aiCellMapping = new GR_index_t[getNumCells()];
  multimap<GR_index_t, GR_index_t>::iterator iter;
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
  bool *aqUpdated = new bool[getNumCells()];
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
      assert(getCell(iFirst)->getType() == Cell::eTriCell);
      TriCell TC;
      TC = *(dynamic_cast<TriCell*>(getCell(iFirst)));
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
    vFoundBug("Cell re-ordering for 2D meshes");
  }
  logMessage(2, "done.\n");
  delete[] aiCellMapping;
  delete[] aqUpdated;
}

void
Mesh2D::reorderFaces()
// Reorder faces according to vertex indices, so that faces with the
// lowest sum of vertex indices are first in order.  The hardest part
// is actually to get the faces and cells reconnected afterwards.
{
  logMessage(2, "Renumbering faces...");
  assert(isSimplicial());
  multimap<GR_index_t, GR_index_t> mOrderedFaces;
  GR_index_t iF, iNF = getNumFaces();
  for (iF = 0; iF < iNF; iF++) {
    Face *pF = getFace(iF);

    GR_index_t iSum = iVertIndex(pF->getVert(0)) + iVertIndex(pF->getVert(1));
    mOrderedFaces.insert(pair<GR_index_t, GR_index_t>(iSum, iF));
  }
  // Now we have a list of all the faces, sorted in order by smallest
  // average of vertex indices!  Use this to create a map of which -old-
  // face becomes which -new- face.
  GR_index_t *aiFaceMapping = new GR_index_t[iNF];
  multimap<GR_index_t, GR_index_t>::iterator iter;
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
  bool *aqUpdated = new bool[iNF];
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
      assert(getFace(iFirst)->getType() == Face::eEdgeFace);
      EdgeFace EF;
      EF = *(dynamic_cast<EdgeFace*>(getFace(iFirst)));
      GR_index_t iNext = iFirst;
      while (aiFaceMapping[iNext] != iFirst) {
	(*getFace(iNext)) = (*getFace(aiFaceMapping[iNext]));
	aqUpdated[iNext] = true;
	iNumUpdated++;
	iNext = aiFaceMapping[iNext];
      }
      (*getFace(iNext)) = EF;
      aqUpdated[iNext] = true;
      iNumUpdated++;
    }
  }

  resetVertexConnectivity();

  if (!(isValid())) {
    logMessage(2, "aborted!\n");
    vFoundBug("Face re-ordering for 2D meshes");
  }
  logMessage(2, "done.\n");
  delete[] aqUpdated;
  delete[] aiFaceMapping;
}

