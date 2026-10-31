#ifndef GR_EntContainer
#define GR_EntContainer

#include "GR_assert.h"
#include "GR_misc.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <list>
#include <omp.h>

// The following definitions allow "only" 2^31 entities of a particular
// type.  That translates to a mesh with 2^31 (2G) faces as the limiting
// case.  That mesh is going to require O(100GB) of memory for faces,
// O(200GB) or so for the whole mesh.  Yes, it's possible to exceed this
// in parallel, but not likely on a single processor.
#define MAX_LISTS 20
#define FIRST_BITS 11
#define FIRST_SIZE (2 << (FIRST_BITS -1 ))
#define FIRST_SIZE_P 8192//(2 << (FIRST_BITS-5))
#define FIRST_SIZE_P_SMALL 8
#define MAX_ENTITIES size_t((size_t(2)) << (FIRST_BITS + MAX_LISTS - 2))
#define MAX_NUM_THREADS 64

#ifndef NUM_PROCS
#define NUM_PROCS 4//omp_get_num_procs()
#endif
///
template<typename T>
  class EntContainer {
    ///
    enum {
      m_maxLists = MAX_LISTS,
      m_firstListSize = FIRST_SIZE,
      p_firstListSize = FIRST_SIZE_P,
      p_firstListSize_small = FIRST_SIZE_P_SMALL,
      p_maxNumThreads = MAX_NUM_THREADS
    };
    /// An array of allocatable arrays of objects of type T.
    T *m_data[m_maxLists];

    std::list<T*> empties;
    typename std::list<T*>::iterator itEmpties;
    GR_index_t Empties_Size;

    /// The allocated size of each list.
    GR_index_t m_listSize[m_maxLists];

    /// The cumulative size of all lists so far, including the current one.
    GR_index_t m_cumSize[m_maxLists];

    /// The number of lists actually in use.
    GR_index_t m_nLists;

    /// The number of entries in use.  This should always be m_lastEntry -
    /// the size of m_empties.
    GR_index_t m_liveEntries;

    /// The index of the last entry in use.
    GR_index_t m_lastEntry;

    /// The overall number of T's allocated; identical to m_cumSize[m_nLists].
    GR_index_t m_totalSize;

    // ****** OMP Private Variables ******
    /// to see whether we are in parallel section or not
    bool isOMPEnabled;

    GR_index_t accumulative_empties;

    // Private member functions

    /// Copy constructor is -really- dubious, because of pointers to data
    /// that's being copied.
    EntContainer(const EntContainer& EC) :
	m_nLists(1), m_liveEntries(0), m_lastEntry(0), m_totalSize(
	    m_firstListSize), empties()
    {
      assert(0);
    }
  public:

    //********** OMP Public Variables and member functions*************
    // this array is to store all the entities that have been modified in
    // the parallel section to be used later in the syncronization

    std::list<T*> *empties_pa;

    double
    sync_entry(bool SizeFree = false)
    {
      double paralleltime = 0;
      if (NUM_PROCS <= 1)
	return paralleltime;
      if (!SizeFree) {
	// This function should be called within a serial section before entering the parallel section
	GR_index_t chunk = empties.size();
	while (chunk < p_firstListSize * NUM_PROCS) {
	  resize(m_totalSize + 1);
	  chunk = empties.size();
	}
	GR_index_t chunkSize = p_firstListSize; //static_cast<GR_index_t> (chunk/double(omp_get_num_procs()));
	GR_index_t chunkLeft = chunk - chunkSize * NUM_PROCS;
	if (m_lastEntry == m_liveEntries)
	  m_lastEntry += chunkSize * NUM_PROCS;
	else if (m_lastEntry - m_liveEntries <= chunkSize * NUM_PROCS)
	  m_lastEntry = m_lastEntry + chunkSize * NUM_PROCS
	      - (m_lastEntry - m_liveEntries);
	//m_lastEntry = m_totalSize-chunkLeft;
	paralleltime = omp_get_wtime();
	omp_set_dynamic(0);
	omp_set_num_threads(NUM_PROCS);
	typename std::list<T*>::iterator iterStart = empties.begin();
	typename std::list<T*>::iterator iterEnd = empties.begin();
#pragma omp parallel default(shared) private (iterStart,iterEnd)
	{
	  /// spreading some initial empty memory slots to all threadprivate arrays

	  int id = omp_get_thread_num();
	  iterStart = empties.begin();
	  iterEnd = empties.begin();
	  std::advance(iterStart, chunkLeft + chunkSize * id);
	  std::advance(iterEnd, chunkLeft + chunkSize * (id + 1));
	  empties_pa[id].insert(empties_pa[id].begin(), iterStart, iterEnd);
	}
	paralleltime = omp_get_wtime() - paralleltime;
	iterStart = empties.begin();
	iterEnd = empties.begin();
	std::advance(iterStart, chunkLeft);
	std::advance(iterEnd, chunk);
	empties.erase(iterStart, iterEnd);
	Empties_Size = empties.size();
	itEmpties = empties.end();
      }
      else {
// This function should be called within a serial section before entering the parallel section
	GR_index_t chunk = empties.size();
	while (chunk < p_firstListSize_small * NUM_PROCS) {
	  resize(m_totalSize + 1);
	  chunk = empties.size();
	}
	GR_index_t chunkSize = p_firstListSize_small; //static_cast<GR_index_t> (chunk/double(omp_get_num_procs()));
	GR_index_t chunkLeft = chunk - chunkSize * NUM_PROCS;
	if (m_lastEntry == m_liveEntries)
	  m_lastEntry += chunkSize * NUM_PROCS;
	else if (m_lastEntry - m_liveEntries <= chunkSize * NUM_PROCS)
	  m_lastEntry = m_lastEntry + chunkSize * NUM_PROCS
	      - (m_lastEntry - m_liveEntries);
	//m_lastEntry = m_totalSize-chunkLeft;
	paralleltime = omp_get_wtime();
	omp_set_dynamic(0);
	omp_set_num_threads(NUM_PROCS);
	typename std::list<T*>::iterator iterStart = empties.begin();
	typename std::list<T*>::iterator iterEnd = empties.begin();
#pragma omp parallel default(shared) private (iterStart,iterEnd)
	{
	  int id = omp_get_thread_num();
	  iterStart = empties.begin();
	  iterEnd = empties.begin();
	  std::advance(iterStart, chunkLeft + chunkSize * id);
	  std::advance(iterEnd, chunkLeft + chunkSize * (id + 1));
	  empties_pa[id].insert(empties_pa[id].begin(), iterStart, iterEnd);
	  /// spreading some initial empty memory slots to all threadprivate arrays
	}
	paralleltime = omp_get_wtime() - paralleltime;
	iterStart = empties.begin();
	iterEnd = empties.begin();
	std::advance(iterStart, chunkLeft);
	std::advance(iterEnd, chunk);
	empties.erase(iterStart, iterEnd);
	Empties_Size = empties.size();
	itEmpties = empties.end();
      }
      return paralleltime;
    }

    void
    sync_exit()
    {
      if (NUM_PROCS > 1) {
	assert(omp_in_parallel() == 0);
	GR_index_t sum = 0;
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  sum += empties_pa[ID].size();
	  empties.splice(empties.end(), empties_pa[ID]);
	}
	//purge();
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Iterator functionality
    //////////////////////////////////////////////////////////////////////

    class iterator {
    private:
      const EntContainer<T>* m_pEC;
      GR_index_t m_whichList;
      T *m_data, *m_listEnd, *m_listBegin;
//    iterator& operator=(const typename EntContainer<T>::iterator&)
//    { assert(0); return this; }
    public:
      iterator&
      operator=(typename EntContainer<T>::iterator& iter)
      {
	m_pEC = iter.m_pEC;
	m_whichList = iter.m_whichList;
	m_data = iter.m_data;
	m_listEnd = iter.m_listEnd;
	m_listBegin = iter.m_listBegin;
	return *this;
      }
      iterator(const typename EntContainer<T>::iterator& iter) :
	  m_pEC(iter.m_pEC), m_whichList(iter.m_whichList), m_data(iter.m_data), m_listEnd(
	      iter.m_listEnd), m_listBegin(iter.m_listBegin)
      {
      }
      iterator(const EntContainer<T>& EC, const bool begin = true) :
	  m_pEC(&EC), m_whichList(0), m_data(NULL), m_listEnd(NULL), m_listBegin(
	  NULL)
      {
	GR_index_t index;
	if (begin) {
	  m_whichList = 0;
	  index = 0;
	}
	else {
	  // Initializing at the end.
	  m_whichList = m_pEC->m_nLists - 1;
	  index = m_pEC->m_lastEntry
	      - (m_whichList == 0 ? 0 : m_pEC->m_cumSize[m_whichList - 1]);
	}
	m_data = m_pEC->m_data[m_whichList] + index;
	m_listBegin = m_pEC->m_data[m_whichList];

	if (begin && index != m_pEC->m_listSize[m_whichList]
	    && m_data->isDeleted()) {
	  increment();
	}

	GR_index_t lastIndex = m_pEC->m_listSize[m_whichList];
	if (m_whichList == m_pEC->m_nLists - 1) {
	  // This can happen at the beginning, too.
	  lastIndex = m_pEC->m_lastEntry
	      - (m_whichList == 0 ? 0 : m_pEC->m_cumSize[m_whichList - 1]);
	}

	m_listEnd = m_pEC->m_data[m_whichList] + lastIndex;
      }
      ~iterator()
      {
      }
      // The compiler-default operator= and copy constructor are fine.
      T&
      operator*() const
      {
	return (*m_data);
      }
      T*
      operator->() const
      {
	return m_data;
      }
      iterator
      operator++(int)
      {
	// Post-increment
	// Return the value that hasn't yet been incremented
	iterator __tmp = *this;
	increment();
	return __tmp;
      }
      iterator&
      operator++()
      {
	// Pre-increment
	// Increment, then return.
	increment();
	return *this;
      }
      bool
      operator==(const iterator& iter) const
      {
	return m_data == iter.m_data;
      }
      bool
      operator!=(const iterator& iter) const
      {
	return m_data != iter.m_data;
      }
      friend std::ostream&
      operator<<(std::ostream& os, const iterator& iter)
      {
	os << "EntContainer: " << iter.m_pEC << " Data: " << iter.m_data;
	return (os);
      }
      GR_index_t
      getindexit()
      {
	return (m_pEC->getIndex(m_data));
      }

      void
      jump(GR_index_t num)
      {
	assert(m_pEC->isValid(num));
	GR_index_t index = this->getindexit() + num;
	assert(m_pEC->isValid(index));
	m_data = m_pEC->getEntry(index);
	m_whichList = m_pEC->getList(index);
	m_listBegin = m_pEC->m_data[m_whichList];
	GR_index_t lastIndex = m_pEC->m_listSize[m_whichList];
	if (m_whichList == m_pEC->m_nLists - 1) {
	  // This can happen at the beginning, too.
	  lastIndex = m_pEC->m_lastEntry
	      - (m_whichList == 0 ? 0 : m_pEC->m_cumSize[m_whichList - 1]);
	}
	m_listEnd = m_pEC->m_data[m_whichList] + lastIndex;
	if (m_data->isDeleted())
	  increment();

      }
    private:
      void
      increment()
      {
	// Warning!  This function cheerfully increments past the end, so
	// you'd better be careful to check the end condition every time...
	// Carefully, so that wraparound will work properly...
//       if (pEC->size() == 0) return;
	do {
	  m_data++;
	  if (m_data == m_listEnd) {
	    if (m_whichList == m_pEC->m_nLists - 1) {
	      return;
	    }
	    else {
	      // Switch over to the next list
	      m_whichList++;
	      m_listBegin = m_data = m_pEC->m_data[m_whichList];
	      GR_index_t lastIndex = m_pEC->m_listSize[m_whichList];
	      if (m_whichList == m_pEC->m_nLists - 1) {
		// This can happen at the beginning, too.
		lastIndex =
		    m_pEC->m_lastEntry
			- (m_whichList == 0 ?
			    0 : m_pEC->m_cumSize[m_whichList - 1]);
	      }

	      m_listEnd = m_pEC->m_data[m_whichList] + lastIndex;
	    }
	  }
	}
	while (m_data->isDeleted());
      }
    };

    //////////////////////////////////////////////////////////////////////
    // Done with iterator stuff.
    //////////////////////////////////////////////////////////////////////

    ///
    ~EntContainer()
    {
      for (GR_index_t i = 0; i < m_nLists; i++) {
	assert(m_data[i] != NULL);
	delete[] m_data[i];
      }
      m_nLists = 0;
      for (int i = 0; i < NUM_PROCS; i++)
	empties_pa[i].clear();
      empties.clear();

      delete[] empties_pa;
      //e_nLists = 0;
    }

    ///
    EntContainer() :
	m_nLists(1), m_liveEntries(0), m_lastEntry(0), m_totalSize(
	    m_firstListSize), isOMPEnabled(false)
    {
      for (GR_index_t ii = 0; ii < m_maxLists; ii++) {
	m_listSize[ii] = 0;
	m_cumSize[ii] = 0;
	m_data[ii] = NULL;
      }
      m_listSize[0] = m_totalSize;
      m_cumSize[0] = m_totalSize;
      m_data[0] = new T[m_totalSize];
      empties.clear();
      accumulative_empties = p_firstListSize;
      for (GR_index_t iter = 0; iter < m_totalSize; ++iter)
	empties.push_back(&(m_data[0][m_totalSize - 1 - iter]));
      itEmpties = empties.end();
      Empties_Size = empties.size();

      empties_pa = new std::list<T*>[NUM_PROCS];

    }

    ///
    iterator
    begin() const
    {
      return iterator(*this, true);
    }
    iterator
    end() const
    {
      return iterator(*this, false);
    }

    iterator
    advance(typename EntContainer<T>::iterator& iter, GR_index_t num)
    {
      assert(isValid(num));
      //iterator iter(*this,true);
      for (GR_index_t i = 0; i < num; i++)
	++iter;
      return iter;
    }

    ///
    void
    clear()
    {
      // Free up all data, then set up an empty container with
      // m_firstListSize entries allocated.
      m_liveEntries = 0;
      m_lastEntry = 0;
      m_nLists = 1;    // e_nLists = 1;
      //e_lastList=0;
      m_totalSize = m_firstListSize;
      for (GR_index_t ii = 0; ii < m_maxLists; ii++) {
	m_listSize[ii] = 0;
	m_cumSize[ii] = 0;
	if (m_data[ii])
	  delete[] m_data[ii];
	m_data[ii] = NULL;
      }
      m_listSize[0] = m_totalSize;
      m_cumSize[0] = m_totalSize;
      m_data[0] = new T[m_totalSize];
      empties.clear();
      for (GR_index_t iter = 0; iter < m_totalSize; ++iter)
	empties.push_back(&(m_data[0][m_totalSize - 1 - iter]));
    }
    ///

    /// \brief The number of entries actually in use.
    ///
    /// Note that there may be holes in the data structure, so the last
    /// entry in use may be at a higher index than this.  So use the
    /// iterators to loop over containers.
    GR_index_t
    size() const
    {
      return m_liveEntries;
    }

    /// \brief The index of the last entry ever used.
    ///
    /// This is a backward compatibility function so that it's possible to
    /// write the old-style integer loops over entities.  This function
    /// will go away, so its use in new code is a very bad idea.
    GR_index_t
    lastEntry() const
    {
      return m_lastEntry;
    }

    /// This function could be used later to get the number of available empty memory slots
    GR_index_t
    emptyavailable() const
    {
      return empties.size();
    }

    /// \brief The number of entries currently allocated.
    GR_index_t
    capacity() const
    {
      return m_cumSize[m_nLists - 1];
    }

    /// \brief The most entities one could store in the data structures as
    /// currently compiled.
    ///
    /// Note that actually using all of this hypothetical total is likely
    /// to blow the virtual memory on your machine.  Caveat programmer.
    GR_index_t
    max_size() const
    {
      return GR_index_t(m_firstListSize) << (m_maxLists);
    }

    bool
    isEmptiesUnique(GR_index_t &checksize)
    {
      std::list<T*> empties1;
      //std::copy(empties.begin(), empties.end(), empties1.begin());
      typename std::list<T*>::iterator iterS;
      typename std::list<T*>::iterator iterE = empties.end();
      for (iterS = empties.begin(); iterS != iterE; iterS++)
	empties1.push_back(*iterS);
      empties1.sort();
      empties1.unique();
      checksize = empties1.size();
      return (empties1.size() == empties.size());
    }

    /// This function is deprecated; new code shouldn't use it.
    void
    reserve(GR_index_t newSize)
    {
      if (newSize < m_totalSize) {
	m_lastEntry = newSize;
	m_liveEntries = 0;
	empties.clear();
	for (GR_index_t counter = newSize; counter < m_totalSize; counter++) {
	  T* bufferCheck = getEntryChecking(counter);
	  empties.push_front(bufferCheck);
	}

	for (GR_index_t i = 0; i < newSize; i++) {
	  GR_index_t iBuff = i;
	  i = newSize - 1 - i;
	  T* pT = getEntry(i);
	  if (pT->isDeleted())
	    empties.push_back(pT);
	  else
	    m_liveEntries++;
	  i = iBuff;
	}
	resize(newSize);
	assert(empties.size() + m_liveEntries == m_totalSize);

      }
      else {
	resize(newSize);
	m_lastEntry = newSize;
	m_liveEntries = 0;
	empties.clear();
	for (GR_index_t counter = newSize; counter < m_totalSize; counter++) {
	  T* bufferCheck = getEntryChecking(counter);
	  empties.push_front(bufferCheck);
	}

	for (GR_index_t i = 0; i < newSize; i++) {
	  GR_index_t iBuff = i;
	  i = newSize - 1 - i;
	  T* pT = getEntry(i);
	  if (pT->isDeleted())
	    empties.push_back(pT);
	  else
	    m_liveEntries++;
	  i = iBuff;
	}
	assert(empties.size() + m_liveEntries == m_totalSize);
      }
    }
    /// \brief Force a change in size for the container.
    ///
    /// Unlike the std::vector reserve function, this one also reduces
    /// size when appropriate.

    void
    purge()
    {
      if (m_lastEntry > 0) {
	T* pTFor = pTFirst();
	T* pTBack = pTLast();

	GR_index_t newSize = 0;
	while (pTFor != pTBack) {
	  while (pTFor != pTBack && !pTFor->isDeleted()) {
	    pTFor = pTNext(pTFor);
	    newSize++;
	  }

	  while (pTBack != pTFor && pTBack->isDeleted()) {
	    pTBack = pTPrev(pTBack);
	  }

	  if (pTFor != pTBack) {
	    assert(pTFor->isDeleted());
	    newSize++;
	    (*pTFor) = (*pTBack);
	    pTBack->markAsDeleted();
	    pTFor = pTNext(pTFor);
	  }
	  else if (!pTFor->isDeleted()) {
	    newSize++;
	  }
	}

	//GR_index_t indexList=0 , indexSum=0;
	empties.clear();
	for (GR_index_t counter = newSize; counter < m_totalSize; counter++) {
	  T*buffPointer = getEntryChecking(counter);
	  empties.push_front(buffPointer);
	}
	//m_liveEntries = m_lastEntry = newSize;
	resize(newSize);
	m_liveEntries = m_lastEntry = newSize;
	assert(m_liveEntries + empties.size() == m_totalSize);
	assert(m_liveEntries < capacity());
      }
    }

    GR_index_t
    purgeLiveEmpties(GR_index_t Numbers)
    {
      for (GR_index_t iterE = 0; iterE < Numbers; iterE++)
	empties.pop_front();
      return (empties.size());
    }

    /// \brief Returns a pointer to an unused entity in the container.
    ///
    /// This entity may be at the end, or we may be re-using one that's
    /// been deleted.

    // ******* getNewEntry() for OMP ********

    T*
    getNewEntry(int ID = 0)
    {
      if (omp_in_parallel() == 0) {
	T*buffer;
	assert(omp_in_parallel() == 0);
	if (empties.empty()) {
	  resize(m_totalSize + 1);
	  m_lastEntry++;
	  m_liveEntries++;
	  buffer = empties.back();
	  empties.pop_back();
	  return (buffer);
	}
	else if (!empties.empty() && m_lastEntry == m_totalSize) {
	  m_liveEntries++;
	  buffer = empties.back();
	  empties.pop_back();
	  return (buffer);
	}
	else if (!empties.empty() && m_lastEntry < m_totalSize) {
	  if (m_lastEntry - m_liveEntries == 0)
	    m_lastEntry++;
	  m_liveEntries++;
	  buffer = empties.back();
	  empties.pop_back();
	  return (buffer);
	}
	else
	  assert(0);
      }
      else {
	T*buff;
	if (empties_pa[ID].empty()) {
	  resize(m_totalSize + 1, ID);
	  buff = empties_pa[ID].back();
	  empties_pa[ID].pop_back();
#pragma omp atomic
	  m_liveEntries++;
	  return (buff);
	}
	else {
	  buff = empties_pa[ID].back();
	  empties_pa[ID].pop_back();
#pragma omp atomic
	  m_liveEntries++;
	  return (buff);

	}
      }
      assert(0);
      return NULL;
    }

    //******* deleteEntry for OMP **********

    void
    deleteEntry(T* dead, int ID = 0)
    {
      if (omp_in_parallel() == 0) {
	assert(contains(dead) && dead->isDeleted());
	empties.push_back(dead);
	m_liveEntries--;
//    		assert(empties.size() + m_liveEntries == m_totalSize);
      }
      else {
	assert(contains(dead) && dead->isDeleted());
	empties_pa[ID].push_back(dead);
	if (m_liveEntries == 0)
	  assert(0);
#pragma omp atomic
	m_liveEntries--;
      }
    }

    //******* resize for OMP*********
    void
    resize(GR_index_t newSize, int ID = 0)
    {
      if (omp_in_parallel() == 0) {
	assert(omp_in_parallel() == 0);
	if (newSize == m_totalSize)
	  return;
	else if (newSize > m_totalSize) {
	  GR_index_t m_nListsOld = m_nLists;
	  GR_index_t iMax = m_cumSize[m_nLists - 1];
	  while (newSize > m_cumSize[m_nLists - 1] && m_nLists < m_maxLists) {
	    m_data[m_nLists] = new T[iMax];
	    m_listSize[m_nLists] = iMax;
	    // At this point we need to allocate the pointers of empty new slots to
	    // the empty container (empties) which stores all the pointers to empty slots.
	    //for (GR_index_t counter=0 ; counter<iMax ; counter++) {fillEmpties(&(m_data[m_nLists][iMax-1-counter]));}
	    iMax = m_totalSize = m_cumSize[m_nLists] = 2 * iMax;
	    m_nLists++;
	  }
	  for (GR_index_t countE2 = m_nListsOld; countE2 < m_nLists;
	      countE2++) {
	    GR_index_t EndSize = m_listSize[countE2];
	    for (GR_index_t countE3 = 0; countE3 < EndSize; countE3++)
	      empties.push_front(&(m_data[countE2][countE3]));
	  }
	  assert(iMax >= newSize && iMax < 2 * newSize);
	  if (iMax < newSize)
	    vFatalError("Out of space", "resizing variable-sized array");
	  assert(empties.size() + m_liveEntries == m_totalSize);
	}
	else {
	  GR_index_t Numbers = m_totalSize;
	  assert(newSize < m_totalSize);
	  assert(empties.size() + m_liveEntries == m_totalSize);
	  while (m_totalSize - m_listSize[m_nLists - 1] >= newSize
	      && m_nLists > 1) {
	    // m_nLists >= 2 or we wouldn't be here.
	    assert(
		m_totalSize - m_listSize[m_nLists - 1]
		    == m_cumSize[m_nLists - 2]);
	    m_nLists--;
	    delete[] m_data[m_nLists];
	    m_totalSize -= m_listSize[m_nLists];
	    m_listSize[m_nLists] = m_cumSize[m_nLists] = 0;
	    m_data[m_nLists] = NULL;
	  }
	  Numbers -= m_totalSize;
	  if (Numbers > 0)
	    purgeLiveEmpties(Numbers);

	  assert(empties.size() + m_liveEntries == m_totalSize);
	}
      }
      else // **************** Parallel resize *****************
      {
#pragma omp critical(resize_pa)
	{
	  if (Empties_Size == 0) {
	    empties.clear();
	    GR_index_t m_nListsOld = m_nLists;
	    GR_index_t iMax = m_cumSize[m_nLists - 1];
	    GR_index_t SumSize = 0;
	    while (newSize > m_cumSize[m_nLists - 1] && m_nLists < m_maxLists) {
	      m_data[m_nLists] = new T[iMax];
	      m_listSize[m_nLists] = iMax;
	      iMax = m_totalSize = m_cumSize[m_nLists] = 2 * iMax;
	      m_nLists++;
	    }
	    for (GR_index_t countE2 = m_nListsOld; countE2 < m_nLists;
		countE2++) {
	      GR_index_t EndSize = m_listSize[countE2];
	      SumSize += EndSize;
	      for (GR_index_t countE3 = 0; countE3 < EndSize; countE3++)
		empties.push_front(&(m_data[countE2][countE3]));
	    }
	    assert(iMax >= newSize && iMax < 2 * newSize);
	    if (iMax < newSize)
	      vFatalError("Out of space", "resizing variable-sized array");
	    Empties_Size = SumSize;
	    itEmpties = empties.end();
	    accumulative_empties *= 2;
	  }
	}
	GR_index_t SizeofEmpty;
	GR_index_t counterStart = 0;
#pragma omp critical(resize_pa)
	{
	  SizeofEmpty = Empties_Size;
	  counterStart = 0;
	  if (p_firstListSize < SizeofEmpty)
	    counterStart = SizeofEmpty - p_firstListSize;
	  Empties_Size = counterStart;
	  m_lastEntry += SizeofEmpty - counterStart;
	}
	typename std::list<T*>::iterator itC = empties.begin();
	typename std::list<T*>::iterator itE;
	std::advance(itC, counterStart);
	itE = itC;
	std::advance(itE, SizeofEmpty - counterStart);
	empties_pa[ID].insert(empties_pa[ID].begin(), itC, itE);
      }
    }

    /// \brief Finds the pointer for a particular index
    T*
    getEntry(GR_index_t i) const
    {
      //assert(omp_in_parallel() == 0);
      assert(isValid(i));
      GR_index_t iList = 0, iSum = 0;
      while (i >= (iSum += m_listSize[iList]) && iList < m_nLists)
	iList++;
      assert(iList < m_nLists);
      iSum -= m_listSize[iList];
      assert(iSum <= i && iSum + m_listSize[iList] > i);
      return (m_data[iList] + (i - iSum));
    }

    GR_index_t
    getList(GR_index_t i) const
    {
      assert(isValid(i));
      GR_index_t iList = 0, iSum = 0;
      while (i >= (iSum += m_listSize[iList]) && iList < m_nLists)
	iList++;
      assert(iList < m_nLists);
      iSum -= m_listSize[iList];
      assert(iSum <= i && iSum + m_listSize[iList] > i);
      return (iList);
    }

    T*
    getEntryChecking(const GR_index_t i) const
    {
      assert(isValid(i, m_totalSize));
      GR_index_t iList = 0, iSum = 0;
      while (i >= (iSum += m_listSize[iList]) && iList < m_nLists)
	iList++;
      assert(iList < m_nLists);
      iSum -= m_listSize[iList];
      assert(iSum <= i && iSum + m_listSize[iList] > i);
      return (m_data[iList] + (i - iSum));
    }

    T*
    getOnlyEntryEmpty(const GR_index_t i)
    {
      assert(isValid(i, empties.size()));
      typename std::list<T*>::iterator iterStart = empties.begin();
      std::advance(iterStart, i);
      T*buff = *iterStart;
      return (buff);
    }

    /// \brief Finds the entity for a particular index
    T&
    operator[](const GR_index_t i) const
    {
      return (*getEntry(i));
    }

    /// \brief Turn a pointer in the container into an index.
    ///
    /// Really useful in two places:  output and debugging.
    GR_index_t
    getIndex(const T* const pT) const
    {
      // This function's runtime depends on the number of lists, so it's
      // technically log(n), with a very small constant.
      GR_index_t iList;
      for (iList = 0;
	  (!((const_cast<T*>(pT) >= m_data[iList])
	      && (const_cast<T*>(pT) < m_data[iList] + m_listSize[iList]))
	      && iList < m_nLists); iList++) {
      }
      if (iList == m_nLists)
	return (~0); // Not a valid pointer.
      assert(iList < m_nLists);
      assert(const_cast<T*>(pT) >= m_data[iList]);
      assert(const_cast<T*>(pT) < m_data[iList] + m_listSize[iList]);
      GR_index_t prevSize = (iList == 0) ? 0 : m_cumSize[iList - 1];
      assert(
	  (reinterpret_cast<size_t>(pT)
	      - reinterpret_cast<size_t>(m_data[iList])) % sizeof(T) == 0);
      return prevSize + (const_cast<T*>(pT) - m_data[iList]);
    }

    /// \brief op= is not safe for containers of mesh entities because of
    /// the pointer interconnections between entities.
  private:
    EntContainer&
    operator=(const EntContainer& EC)
    {
      assert(EC.m_liveEntries >= 0);
      assert(EC.m_totalSize > 0);
      // Do nothing if the source says EC=EC.
      if (this != &EC) {
	// Clean up the existing data.
	int i;
	for (i = 0; i < m_nLists; i++)
	  delete[] m_data[i];

	m_liveEntries = EC.m_liveEntries;
	m_lastEntry = EC.m_lastEntry;
	m_nLists = EC.m_nLists;
	m_totalSize = EC.m_totalSize;
	int iSum = 0;
	for (int ii = 0; ii < m_nLists; ii++) {
	  assert(EC.m_listSize[ii] >= 0);
	  assert(EC.m_data[ii] != NULL);
	  m_listSize[ii] = EC.m_listSize[ii];
	  m_cumSize[ii] = EC.m_cumSize[ii];
	  m_data[ii] = new T[m_listSize[ii]];
	  iSum += m_listSize[ii];
	  for (i = 0; i < m_listSize[ii]; i++)
	    m_data[ii][i] = EC.m_data[ii][i];
	}
	//    assert(iSum == iMax);
      }
      return (*this);
    }
  public:
    ///
    T*
    pTFirst() const
    {
      if (m_liveEntries == 0)
	return (0);
      else
	return (&(m_data[0][0]));
    }
    ///
    T*
    pTLast() const
    {
      if (m_lastEntry == 0)
	return (0);
      else
	return (getEntry(m_lastEntry - 1));
    }
    ///
    T*
    pTNext(T* const pTCurr) const
    {
      assert(contains(const_cast<const T*>(pTCurr)));
      GR_index_t iList = 0, iSum = 0;
      while (!((pTCurr >= m_data[iList])
	  && (pTCurr < m_data[iList] + m_listSize[iList])) && iList < m_nLists)
	iSum += m_listSize[iList++];
      assert(iList < m_nLists);
      assert(pTCurr >= m_data[iList]);
      assert(pTCurr < m_data[iList] + m_listSize[iList]);
      assert(
	  (reinterpret_cast<unsigned long>(pTCurr)
	      - reinterpret_cast<unsigned long>(m_data[iList])) % sizeof(T)
	      == 0);
      // This seems to be wrong.  iSum -= m_listSize[iList];
      if (pTCurr == pTLast()
	  || (pTCurr - m_data[iList] == std::ptrdiff_t(m_listSize[iList] - 1)
	      && iList == m_nLists - 1))
	return 0;
      else if (pTCurr - m_data[iList] == std::ptrdiff_t(m_listSize[iList] - 1))
	return m_data[iList + 1];
      else
	return pTCurr + 1;
    }

    ///
    T*
    pTPrev(T* const pTCurr) const
    {
      assert(contains(const_cast<const T*>(pTCurr)));
      GR_index_t iList = 0, iSum = 0;
      while (!((pTCurr >= m_data[iList])
	  && (pTCurr < m_data[iList] + m_listSize[iList])) && iList < m_nLists)
	iSum += m_listSize[iList++];
      assert(iList < m_nLists);
      assert(pTCurr >= m_data[iList]);
      assert(pTCurr < m_data[iList] + m_listSize[iList]);
      assert(
	  (reinterpret_cast<unsigned long>(pTCurr)
	      - reinterpret_cast<unsigned long>(m_data[iList])) % sizeof(T)
	      == 0);
      // This seems to be wrong.  iSum -= m_listSize[iList];
      if (iList == 0 && pTCurr - m_data[iList] == 0)
	return 0;
      else if (pTCurr - m_data[iList] == 0)
	return m_data[iList - 1] + m_listSize[iList - 1] - 1;
      else
	return pTCurr - 1;
    }
    ///
    bool
    contains(const T* const pT) const
    {
      GR_index_t ind = getIndex(pT);
      return isValid(ind);
    }
  private:
    ///
    bool
    isValid(GR_index_t i) const
    {
      return (i < m_lastEntry);
    }

    bool
    isValid(const GR_index_t i, const GR_index_t endList) const
    {
      return (i < endList);
    }
  };

#endif
