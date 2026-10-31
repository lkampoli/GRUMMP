#include <algorithm>
#include "GR_EntContainer.h"
#include "GR_Vertex.h"
#include "GR_Mesh.h"
#include "GR_ADT.h"

#undef VERIFY_QUERY

ADT::ADT(const GR_index_t iNumObj, const int iNumDim, const eDataT eDType,
	 const double* const pdDataIn) :
    pADT_Left(NULL), pADT_Right(NULL), pADT_Parent(NULL), iLevel(0), iDataInd(
	0), iNDim(iNumDim), eDT(eDType), pdData(NULL)
{
  sADTSort *asADTS = new sADTSort[iNumObj];
  int iNV = iNumVars();
  for (GR_index_t i = 0; i < iNumObj; i++) {
    asADTS[i].vSetData(iNV, pdDataIn[i * iNV], &pdDataIn[i * iNV], i);
  }
  vBuildADT(asADTS, iNumObj);
  for (GR_index_t i = 0; i < iNumObj; i++) {
    asADTS[i].vFreeData();
  }
  delete[] asADTS;
}

ADT::ADT(const Mesh* const pMesh) :
    pADT_Left(NULL), pADT_Right(NULL), pADT_Parent(NULL), iLevel(0), iDataInd(
	0), iNDim(pMesh->getVert(0)->getSpaceDimen()), eDT(ePoints), pdData(
    NULL)
{
  int iNumObj = pMesh->getNumVerts();
  sADTSort *asADTS = new sADTSort[iNumObj];
  double *adDataIn = new double[iNumObj * iNDim];
  for (GR_index_t ii = 0; ii < pMesh->getNumVerts(); ii++) {
    Vert *pV = pMesh->getVert(ii);
    adDataIn[ii * iNDim] = pV->x();
    adDataIn[ii * iNDim + 1] = pV->y();
    if (iNDim == 3) {
      adDataIn[ii * iNDim + 2] = pV->z();
    }
  }
  for (int ii = 0; ii < iNumObj; ii++) {
    asADTS[ii].vSetData(iNDim, adDataIn[ii * iNDim], &adDataIn[ii * iNDim], ii);
  }
  vBuildADT(asADTS, iNumObj);
  for (int i = 0; i < iNumObj; i++) {
    asADTS[i].vFreeData();
  }
  delete[] adDataIn;
  delete[] asADTS;
}

ADT::ADT(const EntContainer<Vert>* const pVAVerts) :
    pADT_Left(NULL), pADT_Right(NULL), pADT_Parent(NULL), iLevel(0), iDataInd(
	0), iNDim(pVAVerts->getEntry(0)->getSpaceDimen()), eDT(ePoints), pdData(
    NULL)
{
  int iNumObj = pVAVerts->lastEntry();
  sADTSort *asADTS = new sADTSort[iNumObj];
  double *adDataIn = new double[iNumObj * iNDim];
  int i;
  for (i = 0; i < iNumObj; i++) {
    adDataIn[i * iNDim] = pVAVerts->getEntry(i)->x();
    adDataIn[i * iNDim + 1] = pVAVerts->getEntry(i)->y();
    if (iNDim == 3) {
      adDataIn[i * iNDim + 2] = pVAVerts->getEntry(i)->z();
    }
  }
  for (i = 0; i < iNumObj; i++) {
    asADTS[i].vSetData(iNDim, adDataIn[i * iNDim], &adDataIn[i * iNDim], i);
  }
  vBuildADT(asADTS, iNumObj);
  delete[] adDataIn;
  delete[] asADTS;
}

void
ADT::vBuildADT(sADTSort asADTS[], const GR_index_t iNumObj)
{
  assert(iNumObj > 0);
  std::sort(asADTS, asADTS + iNumObj);
  // A balanced split, with an extra element into the left subtree if needed.
  int iMedian = iNumObj / 2;

  if (pdData != NULL)
    delete pdData;
  pdData = new double[iNumVars()];

  int i;
  iDataInd = asADTS[iMedian].iInd;
  for (i = iNumVars() - 1; i >= 0; i--) {
    pdData[i] = asADTS[iMedian].pdData[i];
  }

  // Update the keys so that we can sort again at the next level.
  int iNextKey = iDataKey(iLevel + 1);
  for (GR_index_t is = 0; is < iNumObj; is++) {
    asADTS[is].dKey = asADTS[is].pdData[iNextKey];
  }

  GR_index_t iLeftSize = iMedian;
  GR_index_t iRightSize = iNumObj - iMedian - 1;
  if (iLeftSize > 0) {
    pADT_Left = new ADT();
    pADT_Left->pADT_Parent = this;
    pADT_Left->iNDim = iNDim;
    pADT_Left->iLevel = iLevel + 1;
    pADT_Left->eDT = eDT;
    pADT_Left->vBuildADT(asADTS, iLeftSize);
  }
  if (iRightSize > 0) {
    pADT_Right = new ADT();
    pADT_Right->pADT_Parent = this;
    pADT_Right->iNDim = iNDim;
    pADT_Right->iLevel = iLevel + 1;
    pADT_Right->eDT = eDT;
    pADT_Right->vBuildADT(asADTS + iMedian + 1, iRightSize);
  }
}

const ADT*
ADT::pADTFindNodeWithMaxValue(const int iKeyIn) const
{
  const ADT* pADTResult = this;
  int iKey = iDataKey();

  if (iKeyIn != iKey) {
    // In this case, start by searching the left subtree, if it exists.
    if (pADTLeftChild() != NULL) {
      const ADT* pADTRes2 = pADTLeftChild()->pADTFindNodeWithMaxValue(iKeyIn);

      // Now compare to this node.
      if (pADTResult->pdData[iKeyIn] < pADTRes2->pdData[iKeyIn]) {
	pADTResult = pADTRes2;
      }
    }
  }
  // Always check the right subtree, if it exists.
  if (pADTRightChild() != NULL) {
    const ADT* pADTRes2 = pADTRightChild()->pADTFindNodeWithMaxValue(iKeyIn);

    // Now compare to this node.
    if (pADTResult->pdData[iKeyIn] < pADTRes2->pdData[iKeyIn])
      pADTResult = pADTRes2;
  }
  return pADTResult;
}

const ADT*
ADT::pADTFindNodeWithMinValue(const int iKeyIn) const
{
  const ADT* pADTResult = this;
  int iKey = iDataKey();

  if (iKeyIn != iKey) {
    // In this case, start by searching the right subtree, if it exists.
    if (pADTRightChild() != NULL) {
      const ADT* pADTRes2 = pADTRightChild()->pADTFindNodeWithMinValue(iKeyIn);

      // Now compare to this node.
      if (pADTResult->pdData[iKeyIn] > pADTRes2->pdData[iKeyIn]) {
	pADTResult = pADTRes2;
      }
    }
  }
  // Always check the left subtree, if it exists.
  if (pADTLeftChild() != NULL) {
    const ADT* pADTRes2 = pADTLeftChild()->pADTFindNodeWithMinValue(iKeyIn);

    // Now compare to this node.
    if (pADTResult->pdData[iKeyIn] > pADTRes2->pdData[iKeyIn])
      pADTResult = pADTRes2;
  }
  return pADTResult;
}

void
ADT::vBalancingPromotion()
{
  if (pADT_Left == NULL && pADT_Right == NULL) {
    // This node has no descendants.  It needs to be unlinked from the
    // tree.
    if (pADT_Parent == NULL) {
      // The whole tree is gone!
      iLevel = -1;
      if (pdData)
	delete[] pdData;
      pdData = NULL;
    }
    else {
      // Unlink from parent.
      if (pADT_Parent->pADT_Left == this) {
	pADT_Parent->pADT_Left = NULL;
      }
      else {
	assert(pADT_Parent->pADT_Right == this);
	pADT_Parent->pADT_Right = NULL;
      }
      delete this;
    }
    return;
  }

  // Promote from the side that's larger.
  GR_index_t iLeftSize = pADT_Left == NULL ? 0 : pADT_Left->iTreeSize();
  GR_index_t iRightSize = pADT_Right == NULL ? 0 : pADT_Right->iTreeSize();

  int iKey = iDataKey();
  ADT *pADTPromote;
  if (iLeftSize >= iRightSize) {
    // Find the node to promote from the left side.
    assert(iLeftSize > 0);
    pADTPromote = const_cast<ADT*>(pADT_Left->pADTFindNodeWithMaxValue(iKey));
  }
  else {
    // Find the node to promote from the right side.
    assert(iRightSize > 0);
    pADTPromote = const_cast<ADT*>(pADT_Right->pADTFindNodeWithMinValue(iKey));
  }

  iDataInd = pADTPromote->iDataInd;
  assert(iNDim == pADTPromote->iNDim);
  assert(eDT == pADTPromote->eDT);
  for (int i = iNumVars() - 1; i >= 0; i--) {
    pdData[i] = pADTPromote->pdData[i];
  }
  pADTPromote->vBalancingPromotion();
}

GR_index_t
ADT::iTreeSize() const
{
  int iLeftSize = 0, iRightSize = 0;
  if (pADTLeftChild() != NULL) {
    iLeftSize = pADTLeftChild()->iTreeSize();
  }
  if (pADTRightChild() != NULL) {
    iRightSize = pADTRightChild()->iTreeSize();
  }
  return 1 + iLeftSize + iRightSize;
}

// The following #defines can be used as statements.
#define RECUR_LEFT(data, ind) \
      if (pADT_Left) \
	pADT_Left->vAddNode(data, ind); \
      else \
	pADT_Left = new ADT(this, iLevel+1, ind, iNDim, eDT, data)
#define RECUR_RIGHT(data, ind) \
      if (pADT_Right) \
	pADT_Right->vAddNode(data, ind); \
      else \
	pADT_Right = new ADT(this, iLevel+1, ind, iNDim, eDT, data)

void
ADT::vAddNode(const double pdNewData[], const GR_index_t iNewInd)
{
  if (iLevel == -1) {
    // This should be a completely empty tree.
    assert(pADT_Left == NULL);
    assert(pADT_Right == NULL);
    assert(pADT_Parent == NULL);
    iLevel = 0;
    iDataInd = iNewInd;
    assert(iNDim != 0);
    assert(eDT == ePoints || eDT == eBBoxes);
    assert(pdData == NULL);
    pdData = new double[iNumVars()];
    for (int i = 0; i < iNumVars(); i++) {
      pdData[i] = pdNewData[i];
    }
    return;
  }

  int iLeftSize = 0, iRightSize = 0;
  if (pADT_Left) iLeftSize = pADT_Left->iTreeSize();
  if (pADT_Right) iRightSize = pADT_Right->iTreeSize();
  int iKey = iDataKey();
  double dDiv = pdData[iKey];
  if (iLeftSize == iRightSize) {
    // The tree is currently balanced, so who cares where we insert?
    if (pdNewData[iKey] < dDiv) {
      RECUR_LEFT(pdNewData, iNewInd);
    }
    else {
      RECUR_RIGHT(pdNewData, iNewInd);
    }
  }
  else if (iLeftSize > iRightSize) {
    assert(iLeftSize >= 1 && pADT_Left);
    // Need to add a node to the right subtree.
    if (pdNewData[iKey] > dDiv) {
      // New data belongs on right (the easy case)
      RECUR_RIGHT(pdNewData, iNewInd);
    }
    else {
      // First, take the data from the current node and put it into the
      // right subtree.
      RECUR_RIGHT(pdData, iDataInd);

      // Grab the biggest entry in the left subtree.
      ADT *pADT_Cand = const_cast<ADT*>(pADT_Left->pADTFindNodeWithMaxValue(
	  iKey));
      if (pADT_Cand->pdData[iKey] < pdNewData[iKey]) {
	// The new data will now divide the two subtrees.
	iDataInd = iNewInd;
	for (int ii = 0; ii < iNumVars(); ii++) {
	  pdData[ii] = pdNewData[ii];
	}
      }
      else {
	// Promote the candidate you just found, then insert the new
	// data in the left subtree.
	iDataInd = pADT_Cand->iDataInd;
	for (int ii = 0; ii < iNumVars(); ii++) {
	  pdData[ii] = pADT_Cand->pdData[ii];
	}
	pADT_Cand->vBalancingPromotion();
	RECUR_LEFT(pdNewData, iNewInd);
      }
    }
  } // Done with forcing a new node into the right subtree
  else {
    assert(iRightSize >= 1 && pADT_Right);
    // Need to add a node to the left subtree.
    if (pdNewData[iKey] < dDiv) {
      // New data belongs on left (the easy case)
      RECUR_LEFT(pdNewData, iNewInd);
    }
    else {
      // First, take the data from the current node and put it into the
      // left subtree.
      RECUR_LEFT(pdData, iDataInd);

      // Grab the smallest entry in the right subtree.
      ADT *pADT_Cand = const_cast<ADT*>(pADT_Right->pADTFindNodeWithMinValue(
	  iKey));
      if (pADT_Cand->pdData[iKey] > pdNewData[iKey]) {
	// The new data will now divide the two subtrees.
	iDataInd = iNewInd;
	for (int ii = 0; ii < iNumVars(); ii++) {
	  pdData[ii] = pdNewData[ii];
	}
      }
      else {
	// Promote the candidate you just found, then insert the new
	// data in the right subtree.
	iDataInd = pADT_Cand->iDataInd;
	for (int ii = 0; ii < iNumVars(); ii++) {
	  pdData[ii] = pADT_Cand->pdData[ii];
	}
	pADT_Cand->vBalancingPromotion();
	RECUR_RIGHT(pdNewData, iNewInd);
      }
    }
  } // Done with forcing a new node into the left subtree
}

void
ADT::vPointQuery(const double adRange[],
		 std::vector<GR_index_t>& veciResult) const
{
  int iKey = iDataKey();
  double dThisDiv = pdData[iKey];

  // Check the point for membership in the test range.
  bool qAccept = true;
  for (int i = 0; i < iNDim && qAccept; i++)
    qAccept = ((pdData[i] <= adRange[2 * i + 1])
	&& (pdData[i] >= adRange[2 * i]));

  if (qAccept)
    // Add this point.
    veciResult.push_back(iDataInd);

  // Point falls above the bottom of the range.  Check left subtree.
  if ((dThisDiv >= adRange[iKey * 2]) && pADTLeftChild())
    pADTLeftChild()->vPointQuery(adRange, veciResult);

  // Point falls below the top of the range.  Check right subtree.
  if ((dThisDiv <= adRange[iKey * 2 + 1]) && pADTRightChild())
    pADTRightChild()->vPointQuery(adRange, veciResult);
}

void
ADT::vBBoxQuery(const double adRange[],
		std::vector<GR_index_t>& veciResult) const
{
  int iKey = iDataKey();
  double dThisDiv = pdData[iKey];

  // Check the current node for overlap with the test bbox.
  bool qAccept = true;
  for (int i = 0; i < iNDim && qAccept; i++)
    qAccept = ((pdData[2 * i] <= adRange[2 * i + 1])
	&& (pdData[2 * i + 1] >= adRange[2 * i]));

  // The node overlaps; add it to the result list.
  if (qAccept)
    // Add this BBox.
    veciResult.push_back(iDataInd);

  ADT *pADTLeft = pADTLeftChild();
  ADT *pADTRight = pADTRightChild();

  // Left subtree must be checked if it exists and either
  //   (a) the last divider was a minimum or
  //   (b) the low end of the range is below the divider
  if (pADTLeft && (iKey % 2 == 0 || dThisDiv >= adRange[(iKey / 2) * 2]))
    pADTLeft->vBBoxQuery(adRange, veciResult);

  // Right subtree must be checked if it exists and either
  //   (a) the last divider was a maximum or
  //   (b) the high end of the range is above the divider
  if ((pADTRight) && (iKey % 2 == 1 || dThisDiv <= adRange[(iKey / 2) * 2 + 1]))
    pADTRight->vBBoxQuery(adRange, veciResult);
}

ADT::iterator
ADT::begin()
{
  ADT *pADT = this;
  while (pADT->pADT_Left) {
    pADT = pADT->pADT_Left;
  }
  return iterator(pADT);
}

ADT::iterator
ADT::end()
{
  return iterator(NULL);
}

const ADT::iterator
ADT::begin() const
{
  const ADT *pADT = this;
  while (pADT->pADT_Left) {
    pADT = pADT->pADT_Left;
  }
  return iterator(const_cast<ADT*>(pADT));
}

const ADT::iterator
ADT::end() const
{
  return iterator(NULL);
}

#ifdef VERIFY_QUERY
void
ADT::vPointExhaust(const double adRange[],
    std::vector<GR_index_t>& veciResult) const
{
  for (ADT::iterator iter = begin(); iter != end(); ++iter) {
    int iIndex = iter->iDataIndex();
    bool qAccept = true;
    for (int i = 0; i < iNDim && qAccept; i++)
    qAccept = ((iter->pdData[i] <= adRange[2 * i + 1])
	&& (iter->pdData[i] >= adRange[2 * i]));

    // The node overlaps; add it to the result list.
    if (qAccept)
    veciResult.push_back(iIndex);
  }
}

void
ADT::vBBoxExhaust(const double adRange[],
    std::vector<GR_index_t>& veciResult) const
{
  for (ADT::iterator iter = begin(); iter != end(); ++iter) {
    int iIndex = iter->iDataIndex();
    bool qAccept = true;
    for (int i = 0; i < iNDim && qAccept; i++)
    qAccept = ((iter->pdData[2 * i] <= adRange[2 * i + 1])
	&& (iter->pdData[2 * i + 1] >= adRange[2 * i]));
    if (qAccept)
    veciResult.push_back(iIndex);
  }
}
#endif

std::vector<GR_index_t>
ADT::veciRangeQuery(double adRange[]) const
{
  int iDimRange = 2 * iNDim;
  int iCheck;
  for (iCheck = 0; iCheck < iDimRange; iCheck += 2) {
    assert(adRange[iCheck] <= adRange[iCheck + 1]);
    if (iFuzzyComp(adRange[iCheck], adRange[iCheck + 1]) == 0) {
      adRange[iCheck] -= 1.e-8;
      adRange[iCheck + 1] += 1.e-8;
    }
  }
  std::vector<GR_index_t> veciResult;
#ifdef VERIFY_QUERY
  std::vector<GR_index_t> veciCheck;
#endif

  switch (eDT)
    {
    case ePoints:
#ifdef VERIFY_QUERY
      vPointExhaust(adRange, veciCheck);
#endif
      vPointQuery(adRange, veciResult);
      break;
    case eBBoxes:
#ifdef VERIFY_QUERY
      vBBoxExhaust(adRange, veciCheck);
#endif
      vBBoxQuery(adRange, veciResult);
      break;
#ifndef NDEBUG
    default:
      assert(0);
      break;
#endif
    }
#ifdef VERIFY_QUERY
  logMessage(2, "Exhaustive search gave %d items\n", veciCheck.size());
  int ii;
  std::sort(veciCheck.begin(), veciCheck.end());
  for (ii = 0; ii < veciCheck.size(); ii++) {
    logMessage(2, "  %d\n", veciCheck[ii]);
  }
  logMessage(2, "Tree search gave %d items\n", veciResult.size());
  std::sort(veciResult.begin(), veciResult.end());
  for (ii = 0; ii < veciResult.size(); ii++) {
    logMessage(2, "  %d\n", veciResult[ii]);
  }
  return(veciCheck);
#else
  return (veciResult);
#endif
}

void
ADT::vSpew() const
{
  ADT::iterator iter = begin();
  for (; iter != end(); ++iter) {
    logMessage(3, "  Level: %2d  Index:  %3u\n", iter->iTreeLevel(),
	       iter->iDataIndex());
  }
}

bool
ADT::qValid() const
{
  if (iLevel == -1) {
    // This is the special case of an empty tree.
    assert(pADT_Parent == NULL);
    return true;
  }

  if (pADT_Parent) {
    if ((pADT_Parent->pADT_Left != this && pADT_Parent->pADT_Right != this)
	|| (pADT_Parent->iLevel != iLevel - 1) || (pADT_Parent->eDT != eDT)
	|| (pADT_Parent->iNDim != iNDim))
      return false;
  }

  // For a tree to be valid, it must correctly divide its children, and
  // its children must also be valid.
  int iKey = iDataKey();
  double dDiv = pdData[iKey];

  ADT::iterator iter;
  // Check the left subtree
  if (pADT_Left) {
    iter = pADT_Left->begin();
    for (; &(*iter) != this; ++iter) {
      if (iter->pdData[iKey] > dDiv)
	return false;
    }
    if (!pADT_Left->qValid())
      return false;
  }

  // Check the right subtree
  if (pADT_Right) {
    iter = pADT_Right->begin();
    for (; iter->iLevel < iLevel; ++iter) {
      if (iter->pdData[iKey] < dDiv)
	return false;
    }
    if (!pADT_Right->qValid())
      return false;
  }

  return true;
}

Vert*
findVertex(Mesh* const pMesh, const ADT& tree, const double coords[])
{
  double offset = 1.e-10;
  double range[] =
    { coords[0] - offset, coords[0] + offset, coords[1] - offset, coords[1]
	+ offset, coords[2] - offset, coords[2] + offset };
  std::vector<GR_index_t> queryResult = tree.veciRangeQuery(range);
  assert(queryResult.size() <= 1);
  if (queryResult.size() == 1) {
    GR_index_t vertIndex = queryResult[0];
    Vert *result = pMesh->getVert(vertIndex);

    assert(result->x() == coords[0]);
    assert(result->y() == coords[1]);
    assert(result->getSpaceDimen() == 2 || result->z() == coords[2]);
    return result;
  }
  else {
    return pVInvalidVert;
  }
}

void
createADTFromAllCells(Mesh* const pM, ADT*& m_BBoxTree)
{
  // Create a tree of the cells from that mesh to make interpolation
  // efficient
  const int dim = pM->getVert(0)->getSpaceDimen();

  double *a2dBBox = new double[pM->getNumCells() * 2 * dim];
  GR_index_t iCell;
  for (iCell = 0; iCell < pM->getNumCells(); iCell++) {
    for (int ii = 0; ii < dim; ii++) {
      a2dBBox[iCell * dim * 2 + 2 * ii] = LARGE_DBL;
      a2dBBox[iCell * dim * 2 + 2 * ii + 1] = -LARGE_DBL;
    }
    Cell* pC = pM->getCell(iCell);
    for (int iV = 0; iV < pC->getNumVerts(); iV++) {
      const double* adCoord = pC->getVert(iV)->getCoords();
      a2dBBox[iCell * 2 * dim + 0] = min(a2dBBox[iCell * 2 * dim + 0],
					 adCoord[0]);
      a2dBBox[iCell * 2 * dim + 1] = max(a2dBBox[iCell * 2 * dim + 1],
					 adCoord[0]);
      a2dBBox[iCell * 2 * dim + 2] = min(a2dBBox[iCell * 2 * dim + 2],
					 adCoord[1]);
      a2dBBox[iCell * 2 * dim + 3] = max(a2dBBox[iCell * 2 * dim + 3],
					 adCoord[1]);
      if (dim == 3) {
	a2dBBox[iCell * 6 + 4] = min(a2dBBox[iCell * 6 + 4], adCoord[2]);
	a2dBBox[iCell * 6 + 5] = max(a2dBBox[iCell * 6 + 5], adCoord[2]);
      }
    }
  }
  m_BBoxTree = new ADT(pM->getNumCells(), dim, ADT::eBBoxes,
		       reinterpret_cast<double*>(a2dBBox));
  delete[] a2dBBox;
}

Cell*
findCellContainingPoint(const Mesh& M, ADT* BBoxTree, const double adCoord[3],
			double adBestBary[4])
{
  std::vector<GR_index_t> veciQueryResult;
  Cell *pC;
  // Find cell (nearly?) containing the target point
  double dEps = 1.e-8;
  do {
    double adBB[] =
      { adCoord[0] - dEps, adCoord[0] + dEps, adCoord[1] - dEps, adCoord[1]
	  + dEps, adCoord[2] - dEps, adCoord[2] + dEps };
    veciQueryResult = BBoxTree->veciRangeQuery(adBB);
    dEps *= 10;
  }
  while (veciQueryResult.empty());
  // Find out which cell contains the vert, or which one comes closest
  // if none contain it.  "Closest" is defined in the sense of having
  // barycentric coordinates that violate the range (0,1) by as little
  // as possible.
  double dBestValue = 1.e50;
  int iBestCell = -1;
  for (int iQ = veciQueryResult.size() - 1; iQ >= 0 && dBestValue > 0.5; iQ--) {
    // Get barycentric coords.  If best_value <(=) 0.5 on completion
    // then the vert lies within (on an edge of) the cell in list
    int iCell = veciQueryResult[iQ];
    pC = M.getCell(iCell);
    double adBary[] =
      { 0.3, 0.3, 0.4, 0 }; // Always valid
    dynamic_cast<SimplexCell*>(pC)->calcBarycentricCoords(adCoord, adBary);
    double dThisValue = max(max(fabs(adBary[0] - 0.5), fabs(adBary[1] - 0.5)),
			    max(fabs(adBary[2] - 0.5), fabs(adBary[3] - 0.5)));
    if (dThisValue < dBestValue) {
      iBestCell = iCell;
      dBestValue = dThisValue;
      adBestBary[0] = adBary[0];
      adBestBary[1] = adBary[1];
      adBestBary[2] = adBary[2];
      adBestBary[3] = adBary[3];
    }
  } // Done checking all cells which have bboxes overlapping with the
    // vertex
  assert(iBestCell >= 0);
  pC = M.getCell(iBestCell);

  return pC;
}

