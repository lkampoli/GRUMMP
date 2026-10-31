#ifndef GR_ADT

#define GR_ADT 1

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_EntContainer.h"

struct sADTSort {
public:
  double dKey, *pdData;
  GR_index_t iInd;
  sADTSort() :
      dKey(0), pdData(NULL), iInd(0)
  {
  }
  sADTSort(const sADTSort& S) :
      dKey(S.dKey), pdData(S.pdData), iInd(S.iInd)
  {
  }
  sADTSort(const int iDataSize, const double dKeyIn, const double adDataIn[],
	   const GR_index_t iIndIn) :
      dKey(dKeyIn), pdData(0), iInd(iIndIn)
  {
    assert(iDataSize >= 0);
    pdData = new double[iDataSize];
    for (int i = 0; i < iDataSize; i++) {
      pdData[i] = adDataIn[i];
    }
  }
  ~sADTSort()
  {
  }
  sADTSort&
  operator=(const sADTSort& S)
  {
    if (&S != this) {
      dKey = S.dKey;
      pdData = S.pdData;
      iInd = S.iInd;
    }
    return (*this);
  }
  void
  vSetData(const int iDataSize, const double dKeyIn, const double adDataIn[],
	   const GR_index_t iIndIn)
  {
    assert(iDataSize >= 0);
    dKey = dKeyIn;
    iInd = iIndIn;
    if (pdData != NULL)
      delete[] pdData;
    pdData = new double[iDataSize];
    for (int i = 0; i < iDataSize; i++) {
      pdData[i] = adDataIn[i];
    }
  }
  void
  vFreeData()
  {
    delete[] pdData;
  }
  void
  vSetKey(const double dKeyIn)
  {
    dKey = dKeyIn;
  }
  bool
  operator<(const sADTSort& sADTS) const
  {
    int iFC = iFuzzyComp(dKey, sADTS.dKey);
    return (iFC == -1 || (iFC == 0 && iInd < sADTS.iInd));
  }
};

class ADTiterator;

class ADT {
public:
  enum eDataT {
    ePoints = 1, eBBoxes = 2
  };
private:
  // Pointers to children and parent.
  ADT *pADT_Left, *pADT_Right, *pADT_Parent;
  // Level within the tree, original index of the data in the input array.
  int iLevel;
  GR_index_t iDataInd;
  // Number of spatial dimensions of the data
  int iNDim;
  // Data is points or bounding boxes?
  eDataT eDT;
  // A dynamic array storing the data associated with this node.
  double *pdData;
  ADT(const ADT&); // Not defined
  const ADT&
  operator=(const ADT&);  // Not defined
  void
  vCloneTree(const ADT* const pADTIn);
  void
  vBuildADT(sADTSort asADTS[], const GR_index_t iNumObj);
public:
  ADT(ADT * const pADT_ParentIn = NULL, const int iLevIn = -1,
      const GR_index_t iDataIndIn = 0, const int iNDimIn = 0,
      const eDataT eDTIn = ePoints, const double * const pdDataIn = NULL) :
      pADT_Left(NULL), pADT_Right(NULL), pADT_Parent(pADT_ParentIn), iLevel(
	  iLevIn), iDataInd(iDataIndIn), iNDim(iNDimIn), eDT(eDTIn), pdData(
      NULL)
  {
    if (pdDataIn && iNumVars() != 0) {
      pdData = new double[iNumVars()];
      for (int i = iNumVars() - 1; i >= 0; i--) {
	pdData[i] = pdDataIn[i];
      }
    }
  }
  ///
  ADT(const GR_index_t iNumObj, const int iNumDim, const eDataT eDType,
      const double* const pdData);
  ///
  ADT(const Mesh* const pMesh);
  ADT(const EntContainer<Vert>* const ECVerts);
  ///
  ~ADT()
  {
    if (pdData) {
      delete[] pdData;
      pdData = NULL;
    }
    if (pADT_Left)
      delete pADT_Left;
    if (pADT_Right)
      delete pADT_Right;
  }
  ///
  GR_index_t
  iDataIndex() const
  {
    return (iDataInd);
  }
  ///
  int
  iNumVars() const
  {
    return (eDT == ePoints) ? iNDim : iNDim * 2;
  }
  ///
  int
  iDataKey() const
  {
    return iLevel % iNumVars();
  }
  ///
  int
  iDataKey(const int iLev) const
  {
    return iLev % iNumVars();
  }
  ///
  void
  vSetDataType(const eDataT eDTIn)
  {
    eDT = eDTIn;
  }
  ///
  void
  vSetDimensions(const int iNumDim)
  {
    assert(iNumDim >= 1 && iNumDim <= 3);
    iNDim = iNumDim;
  }
  ///
  void
  vCopyData(const ADT& ADT_In)
  {
    assert(eDT == ADT_In.eDT);
    assert(iNDim == ADT_In.iNDim);
    assert(pdData != NULL);
    assert(ADT_In.pdData != NULL);
    iLevel = ADT_In.iLevel;
    iDataInd = ADT_In.iDataInd;
    for (int i = iNumVars() - 1; i >= 0; i--) {
      pdData[i] = ADT_In.pdData[i];
    }
  }
  int
  iTreeLevel() const
  {
    return iLevel;
  }
  ///
  ADT*
  pADTLeftChild() const
  {
    return pADT_Left;
  }
  ///
  ADT*
  pADTRightChild() const
  {
    return pADT_Right;
  }
  ///
  ADT*
  pADTParent() const
  {
    return pADT_Parent;
  }
  ///
  void
  vSetNodeData(const GR_index_t iIndex, const double adData[]);
  ///
  int
  iNumDimensions() const
  {
    return iNDim;
  }
  ///
  void
  vAddNode(const double adNewData[], const GR_index_t iInd);
  ///
  std::vector<GR_index_t>
  veciRangeQuery(double adRange[]) const;
private:
  void
  vPointQuery(const double adRange[], std::vector<GR_index_t>& LiResult) const;
  void
  vBBoxQuery(const double adRange[], std::vector<GR_index_t>& LiResult) const;
#ifdef VERIFY_QUERY
  void
  vPointExhaust(const double adRange[],
		std::vector<GR_index_t>& LiResult) const;
  void
  vBBoxExhaust(const double adRange[], std::vector<GR_index_t>& LiResult) const;
#endif
public:
  ///
  const ADT*
  pADTFindNodeWithMaxValue(const int iKey) const;
  ///
  const ADT*
  pADTFindNodeWithMinValue(const int iKey) const;
  ///
  void
  vRemoveNode(const double adData[], const GR_index_t iInd);
  ///
  void
  vBalancingPromotion();
  ///
  GR_index_t
  iTreeSize() const;
  /// Needed for iterators to work properly.
  typedef ADTiterator iterator;
  iterator
  begin();
  iterator
  end();
  const iterator
  begin() const;
  const iterator
  end() const;
  void
  vSpew() const;
  bool
  qValid() const;
};
// End of ADT declaration

class ADTiterator {
private:
  ADT *pADT;
public:
  ADTiterator(ADT * const pADTIn = NULL) :
      pADT(pADTIn)
  {
  }
  ADTiterator(const ADTiterator& Ai) :
      pADT(Ai.pADT)
  {
  }
  ~ADTiterator()
  {
  }

  ADTiterator&
  operator=(const ADTiterator& Ai)
  {
    if (&Ai != this) {
      pADT = Ai.pADT;
    }
    return (*this);
  }
  ADT&
  operator*() const
  {
    return *pADT;
  }
  ADT*
  operator->() const
  {
    return pADT;
  }
  ADTiterator
  operator++(int)
  {
    if (pADT->pADTRightChild()) {
      *this = pADT->pADTRightChild()->begin();
    }
    else {
      while (pADT->pADTParent() != NULL
	  && pADT->pADTParent()->pADTRightChild() == pADT) {
	pADT = pADT->pADTParent();
      }
      pADT = pADT->pADTParent();
    }
    return *this;
  }
  ADTiterator&
  operator++()
  {
    if (pADT->pADTRightChild()) {
      *this = pADT->pADTRightChild()->begin();
    }
    else {
      while (pADT->pADTParent() != NULL
	  && pADT->pADTParent()->pADTRightChild() == pADT) {
	pADT = pADT->pADTParent();
      }
      pADT = pADT->pADTParent();
    }
    return *this;
  }
  friend bool
  operator==(const ADTiterator& ADTi0, const ADTiterator& ADTi1)
  {
    return ADTi0.pADT == ADTi1.pADT;
  }
  friend bool
  operator!=(const ADTiterator& ADTi0, const ADTiterator& ADTi1)
  {
    return ADTi0.pADT != ADTi1.pADT;
  }
};

Vert*
findVertex(Mesh* const pMesh, const ADT& tree, const double coords[]);

void
createADTFromAllCells(Mesh* const pM, ADT*& m_BBoxTree);

Cell*
findCellContainingPoint(const Mesh& M, ADT* BBoxTree, const double adCoord[3],
		       double adBestBary[4]);

#endif
