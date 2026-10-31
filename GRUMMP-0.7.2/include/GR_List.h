#ifndef GR_List
#define GR_List 1

#include "GR_config.h"
#include "GR_misc.h"
#include <algorithm>
#include <stdlib.h>

template<class T>
  class List {
    int iMaxElem, iNumElem;
    bool qSorted, qAllowSorting;
    T* pCData;
    inline void
    vResize();
  public:
    List() :
	iMaxElem(10), iNumElem(0), qSorted(false), qAllowSorting(true), pCData(
	NULL)
    {
      pCData = new T[iMaxElem];
    }
    List(const int iSize) :
	iMaxElem(iSize), iNumElem(0), qSorted(false), qAllowSorting(true), pCData(
	NULL)
    {
      assert(iSize > 0);
      pCData = new T[iMaxElem];
    }
    inline
    List(const List<T>& CListIn);
    ~List()
    {
      if (pCData) {
	delete[] pCData;
	pCData = NULL;
      }
    }
    inline List<T>&
    operator=(const List<T>& CRHS);
    void
    vDisallowSorting()
    {
      qAllowSorting = false;
    }
    inline T&
    operator[](int iIndex) const;
    int
    iLength() const
    {
      return iNumElem;
    }
    inline bool
    qInList(const T& CTarget);
    inline bool
    qInList(const T& CTarget) const;
    inline int
    iSeek(const T& CTarget) const;
    inline void
    vAddItem(const T& CNewElem);
    inline void
    vAppendItem(const T& CUniqElem);
    inline void
    vAddList(const List<T>& CList);
    inline void
    vAppendList(const List<T>& CList);
    inline void
    vOrderList();
    inline void
    vUniquifyList();
    inline void
    vDelete(const T& CElem);
    inline void
    vDeleteByIndex(const int i);
    inline void
    vClear();
    // The following function is a crude hack, necessary because sometimes
    // the default constructor doesn't seem to get called.
    inline void
    vReset()
    {
      iMaxElem = 10;
      iNumElem = 0;
      qSorted = false, qAllowSorting = true;
      pCData = new T[iMaxElem];
    }
  };

template<class T>
  inline
  List<T>::List(const List<T>& CListIn) :
      iMaxElem(CListIn.iMaxElem), iNumElem(CListIn.iNumElem), qSorted(
	  CListIn.qSorted), qAllowSorting(CListIn.qAllowSorting), pCData()
  {
    assert(CListIn.iNumElem >= 0);
    assert(CListIn.iMaxElem >= CListIn.iNumElem);
    assert(CListIn.pCData != NULL);
    pCData = new T[iMaxElem];
    for (int i_ = 0; i_ < iNumElem; i_++)
      pCData[i_] = CListIn.pCData[i_];
  }

template<class T>
  inline List<T>&
  List<T>::operator=(const List<T>& CRHS)
  {
    if (this != &CRHS) {
      // Clean up current contents, copy contents from CRHS.  Max size 
      // is shrunk to the smallest possible multiple of iSizeInc.
      delete[] pCData;
      iNumElem = CRHS.iNumElem;
      iMaxElem = iNumElem + 10;
      pCData = new T[iMaxElem];
      for (int i_ = 0; i_ < iNumElem; i_++)
	pCData[i_] = CRHS.pCData[i_];
      qSorted = CRHS.qSorted;
      qAllowSorting = CRHS.qAllowSorting;
    }
    return (*this);
  }

template<class T>
  inline T&
  List<T>::operator[](int iIndex) const
  {
    assert(iIndex >= 0);
    assert(iIndex < iNumElem);
    return (pCData[iIndex]);
  }

template<class T>
  inline bool
  List<T>::qInList(const T& CTarget)
  {
    // If the list is big, order it to allow binary searching
    if (!qSorted && iNumElem > 100 && qAllowSorting)
      vOrderList();
    int iBotOfRange = 0;
    int iTopOfRange = iNumElem - 1;
    // Do binary search to narrow the list to be checked, then check
    // sequentially at the end.  If the list is unsorted, just check
    // sequentially. 
    if (qSorted) {
      while (iBotOfRange < iTopOfRange - 8) {
	int iTestElem = (iTopOfRange + iBotOfRange) / 2;
	if (pCData[iTestElem] == CTarget)
	  return (true);
	else if (pCData[iTestElem] < CTarget)
	  iBotOfRange = iTestElem;
	else
	  iTopOfRange = iTestElem;
      }
    }
    for (int i_ = iBotOfRange; i_ <= iTopOfRange; i_++)
      if (pCData[i_] == CTarget)
	return (true);
    return (false);
  }

template<class T>
  inline bool
  List<T>::qInList(const T& CTarget) const
  {
    // If the list is big, order it to allow binary searching
    int iBotOfRange = 0;
    int iTopOfRange = iNumElem - 1;
    // Do binary search to narrow the list to be checked, then check
    // sequentially at the end.  If the list is unsorted, just check
    // sequentially. 
    if (qSorted) {
      while (iBotOfRange < iTopOfRange - 8) {
	int iTestElem = (iTopOfRange + iBotOfRange) / 2;
	if (pCData[iTestElem] == CTarget)
	  return (true);
	else if (pCData[iTestElem] < CTarget)
	  iBotOfRange = iTestElem;
	else
	  iTopOfRange = iTestElem;
      }
    }
    for (int i_ = iBotOfRange; i_ <= iTopOfRange; i_++)
      if (pCData[i_] == CTarget)
	return (true);
    return (false);
  }

template<class T>
  inline int
  List<T>::iSeek(const T& CTarget) const
  {
    for (int i_ = 0; i_ < iNumElem; i_++)
      if (pCData[i_] == CTarget)
	return (i_);
    return (-1);
  }

template<class T>
  inline void
  List<T>::vAddItem(const T& CNewElem)
  {
    for (int i_ = 0; i_ < iNumElem; i_++)
      if (CNewElem == pCData[i_])
	return;
    if (iNumElem == iMaxElem)
      vResize();
    pCData[iNumElem] = CNewElem;
    iNumElem++;
    qSorted = false;
    return;
  }

template<class T>
  inline void
  List<T>::vAppendItem(const T& CUniqElem)
  {
    if (iNumElem == iMaxElem)
      vResize();
    pCData[iNumElem] = CUniqElem;
    iNumElem++;
    qSorted = false;
    return;
  }

template<class T>
  inline void
  List<T>::vAppendList(const List<T>& CList)
  {
    for (int i_ = 0; i_ < CList.iNumElem; i_++)
      vAppendItem(CList[i_]);
    qSorted = false;
  }

template<class T>
  inline void
  List<T>::vAddList(const List<T>& CList)
  {
    for (int i_ = 0; i_ < CList.iNumElem; i_++)
      vAddItem(CList[i_]);
    qSorted = false;
  }

template<class T>
  inline void
  List<T>::vOrderList()
  {
    std::sort(pCData, pCData + iNumElem);
    //  qsort((void*) pCData, (size_t) iNumElem, sizeof(T), iCompPtr);
    qSorted = true;
  }

template<class T>
  inline void
  List<T>::vUniquifyList()
  {
    assert(qAllowSorting);
    vOrderList();
    for (int i_ = iNumElem - 1; i_ > 0; i_--)
      if (pCData[i_] == pCData[i_ - 1]) {
	pCData[i_] = pCData[iNumElem - 1];
	iNumElem--;
      }
  }

template<class T>
  inline void
  List<T>::vDelete(const T& CElem)
  {
    int i_;
    for (i_ = 0; pCData[i_] != CElem && i_ < iNumElem; ++i_) {
    }
    if (i_ < iNumElem) {
      pCData[i_] = pCData[iNumElem - 1];
      --iNumElem;
      qSorted = false;
    }
  }

template<class T>
  inline void
  List<T>::vDeleteByIndex(const int iToDel)
  {
    assert(iToDel >= 0 && iToDel < iNumElem);
    pCData[iToDel] = pCData[iNumElem - 1];
    --iNumElem;
    qSorted = false;
  }

template<class T>
  inline void
  List<T>::vClear()
  {
    iNumElem = 0;
    iMaxElem = 10;
    if (pCData == NULL)
      // Kludge to fix a case where constructor wasn't called for a static
      // List variable.
      qAllowSorting = true;
    else
      delete[] pCData;
    pCData = new T[iMaxElem];
    qSorted = false;
  }

template<class T>
  inline void
  List<T>::vResize()
  {
    iMaxElem *= 2;
    T* pCNewData = new T[iMaxElem];
    for (int i_RS = 0; i_RS < iNumElem; i_RS++)
      pCNewData[i_RS] = pCData[i_RS];
    delete[] pCData;
    pCData = pCNewData;
  }

template<class T>
  inline int
  operator==(const List<T>& CA, const List<T>&CB)
  {
    if (CA.iLength() != CB.iLength())
      return (0);
    for (int i_ = CA.iLength() - 1; i_ >= 0; i_--)
      if (!CB.qInList(CA[i_]))
	return (0);
    return (1);
  }

#undef inline_loop

#endif 

