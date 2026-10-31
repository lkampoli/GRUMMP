#define BOOST_TEST_MODULE TEST_REZA
#include <boost/test/unit_test.hpp>
#include <GR_EntContainer.h>
#include <GR_Mesh2D.h>
#include <omp.h>

struct TestContainer {
  GR_index_t LiveEntries;
  GR_index_t LastEntry;
  GR_index_t TotalSize;
  GR_index_t NEmpties;
  GR_index_t first_size;
  GR_index_t ENum_Pointers;

  //************ declaration of empty containers *************;
  EntContainer<TriCell> m_TECTri;

  void
  CheckEntries(GR_index_t i, GR_index_t j, GR_index_t k)
  {
    BOOST_CHECK_EQUAL(i, LiveEntries);
    BOOST_CHECK_EQUAL(j, LastEntry);
    BOOST_CHECK_EQUAL(k, NEmpties);
  }

  void
  UpdateEntries()
  {
    LiveEntries = m_TECTri.size();
    LastEntry = m_TECTri.lastEntry();
    TotalSize = m_TECTri.capacity();
    NEmpties = LastEntry - LiveEntries;
    first_size = FIRST_SIZE;
    ENum_Pointers = m_TECTri.emptyavailable();
  }

  //************** putting something into containers*********
  void
  InsertingData(int ii, const int ID = 0)
  {
    if (omp_in_parallel() == 0) {
      for (int i = 0; i < ii; i++) {
	TriCell *pTC = m_TECTri.getNewEntry();
	pTC->setDefaultFlags();
	//pTC=dynamic_cast<TriCell*> (i+1);
      }
    }
    else {
      assert(omp_in_parallel() != 0);
      //			#pragma omp critical (UnitTestInsertingData)
      //			{
      for (int i = 0; i < ii; i++) {
	TriCell *pTC = m_TECTri.getNewEntry(ID);
	pTC->setDefaultFlags();
      }
      //			}

    }
  }

  void
  DeleteFirstData()
  {
    TriCell *pTC = m_TECTri.pTFirst();
    if (pTC == 0) {
      BOOST_WARN_MESSAGE(pTC != 0, "the containers are already empty");
    }
    else {
      //logMessage(3,"the pointer to the first slot of memory is: "<<pTC<<std::endl;
      pTC->markAsDeleted();
      m_TECTri.deleteEntry(static_cast<TriCell*>(pTC));
      //pTC=NULL;//dynamic_cast<TriCell*>pCInvalidCell;
      //Mesh2D::deleteCellEvent(pTC);
    }
  }

  void
  DeleteData(const GR_index_t i, const int ID = 0)
  {
    if (omp_in_parallel() == 0) {
      TriCell *pTC1 = m_TECTri.getEntry(i);
      assert(m_TECTri.contains(pTC1));
      pTC1->markAsDeleted();
      m_TECTri.deleteEntry(static_cast<TriCell*>(pTC1));
    }
    else {
      assert(omp_in_parallel() != 0);
      //			#pragma omp critical (UnitTestDeleteData)
      //			{
      TriCell *pTC1 = m_TECTri.getEntry(i);
      assert(m_TECTri.contains(pTC1));
      pTC1->markAsDeleted();
      m_TECTri.deleteEntry(static_cast<TriCell*>(pTC1), ID);
      //			}
    }
  }

  void
  SyncEntry()
  {
    m_TECTri.sync_entry();
  }

  void
  SyncExit()
  {
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      std::list<TriCell*>::iterator iterS = m_TECTri.empties_pa[ID].begin();
      std::list<TriCell*>::iterator iterE = m_TECTri.empties_pa[ID].end();
      for (; iterS != iterE; ++iterS)
	(*iterS)->markAsDeleted();
    }
    m_TECTri.sync_exit();
  }

  void
  likeReconfigure(TriCell * const pTC, const int ID)
  {
    GR_index_t index = m_TECTri.getIndex(pTC);
    DeleteData(index, ID);
    InsertingData(1, ID);
  }

  bool
  IsListUnique(std::list<TriCell*> CheckList)
  {
    std::list<TriCell*> CheckBuffer;
    std::list<TriCell*>::iterator iterS = CheckList.begin();
    std::list<TriCell*>::iterator iterE = CheckList.end();
    for (; iterS != iterE; iterS++)
      CheckBuffer.push_back(*iterS);
    CheckBuffer.sort();
    CheckBuffer.unique();
    return (CheckList.size() == CheckBuffer.size());
  }

  bool
  isDeleted(const GR_index_t i)
  {
    TriCell *pTC1 = m_TECTri.getEntry(i);
    assert(m_TECTri.contains(pTC1));
    return (pTC1->isDeleted());
  }

  void
  Resize(GR_index_t i)
  {
    m_TECTri.resize(i);
  }

  void
  clear()
  {
    m_TECTri.clear();
  }

  TriCell *
  First()
  {
    return m_TECTri.pTFirst();
  }

  TriCell *
  Last()
  {
    return m_TECTri.pTLast();
  }

  TriCell *
  Next(TriCell *pTCPre)
  {
    return m_TECTri.pTNext(pTCPre);
  }

  TriCell *
  Previous(TriCell *pTCNow)
  {
    return m_TECTri.pTPrev(pTCNow);
  }

  void
  purge()
  {
    m_TECTri.purge();
  }

  void
  Reserve(GR_index_t newSize)
  {
    m_TECTri.reserve(newSize);
  }

  bool
  OpEqu(EntContainer<TriCell>::iterator iter1,
	EntContainer<TriCell>::iterator iter2)
  {
    return (iter1 == iter2);
  }

  bool
  OpNEqu(EntContainer<TriCell>::iterator iter1,
	 EntContainer<TriCell>::iterator iter2)
  {
    return (iter1 != iter2);
  }

  EntContainer<TriCell>::iterator
  OpPreInc(EntContainer<TriCell>::iterator& iter)
  {
    EntContainer<TriCell>::iterator iterbuffer = ++iter;
    return (iterbuffer);
  }

  EntContainer<TriCell>::iterator
  OpPostInc(EntContainer<TriCell>::iterator& iter)
  {
    EntContainer<TriCell>::iterator iterbuffer = iter++;
    return (iterbuffer);
  }

  EntContainer<TriCell>::iterator
  Entity_begin() const
  {
    return m_TECTri.begin();
  }
  EntContainer<TriCell>::iterator
  Entity_end() const
  {
    return m_TECTri.end();
  }

  GR_index_t
  getIndex(const TriCell *pTC) const
  {
    return m_TECTri.getIndex(pTC);
  }

  TriCell*
  getEntry(const GR_index_t i) const
  {
    return m_TECTri.getEntry(i);
  }

  TriCell*
  getEntryChecking(const GR_index_t i) const
  {
    return m_TECTri.getEntryChecking(i);
  }

  TriCell*
  getEmEntry(const GR_index_t i)
  {
    return m_TECTri.getOnlyEntryEmpty(i);
  }
};

BOOST_AUTO_TEST_CASE(INITIAL_VALUES)
{
  logMessage(3, "Testing the default EntContainer Constructor... ");

  TestContainer TCell1;
  TCell1.UpdateEntries();
  TCell1.CheckEntries(0, 0, 0);
  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(INSERT_DELETE_DATA)
{
  logMessage(
      3,
      "Testing simple insertion and deletion via getNewEntry and deleteEntry... ");

  TestContainer TCell;
  TCell.InsertingData(4);
  TCell.UpdateEntries();
  BOOST_CHECK_MESSAGE(TCell.LiveEntries == 4,
		      "there's a conflict with the LiveEntries in containers");
  BOOST_CHECK_MESSAGE(TCell.LastEntry == 4,
		      "there's a conflict with the LastEntry in containers");
  BOOST_CHECK_MESSAGE(TCell.NEmpties == 0,
		      "Empties' vector should be empty but it's not");

  TCell.DeleteFirstData();
  TCell.UpdateEntries();
  BOOST_CHECK_MESSAGE(
      TCell.LiveEntries == 3,
      "there's a conflict with the LiveEntries because it should be 4-1=3");
  BOOST_CHECK_MESSAGE(
      TCell.LastEntry == 4,
      "after deletion of first entity the last entity should not change");
  BOOST_CHECK_MESSAGE(
      TCell.NEmpties == 1,
      "Empties' vector should possess one slot of memory right now but it's not");

  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(ITERATORS)
{
  logMessage(3, "Testing iterator class... ");

  TestContainer TCellI;
  TCellI.InsertingData(10);
  TCellI.DeleteFirstData();
  EntContainer<TriCell>::iterator TCiter = TCellI.Entity_begin(), TCiterEnd =
      TCellI.Entity_end();

  TriCell *pTCIT = &(*TCiter);
  GR_index_t ind = TCellI.getIndex(pTCIT);
  BOOST_CHECK_MESSAGE(ind == 1,
		      "This error means that getIndex is working wrong");
  TriCell *pTCIT1 = TCellI.getEntry(ind);
  BOOST_CHECK_MESSAGE(pTCIT == pTCIT1,
		      "This error means that getEntry is working wrong");

  for (GR_index_t i = 1; TCiter != TCiterEnd; ++TCiter, ++i) {
    TriCell *pTCITl = &(*TCiter);
    TriCell *pTCITl1 = TCellI.getEntry(i);
    BOOST_CHECK_EQUAL(pTCITl, pTCITl1);
  }

  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(TWO_LISTS)
{
  logMessage(
      3,
      "Testing insertion and deletion when there are more than two lists... ");

  TestContainer TCell2;
  TCell2.UpdateEntries();
  GR_index_t numindex = TCell2.first_size;
  TCell2.InsertingData(numindex + 1);
  TCell2.UpdateEntries();
  BOOST_CHECK_EQUAL(2 * TCell2.first_size, TCell2.TotalSize); // we have created the new list and the size of the new list should be m_cumSize[0]=firstlist..
  // which in total the size of all the allocated slots should be (1+1)*firstlist
  BOOST_CHECK_EQUAL(TCell2.LastEntry - TCell2.first_size, 1); // the last entry in use should be the first slot of memory in the second list.
  TCell2.InsertingData(numindex);
  TCell2.UpdateEntries();
  BOOST_CHECK_EQUAL(4 * TCell2.first_size, TCell2.TotalSize); // we should have jumped to the third list by now.
  // and the size of the third list should be m_cumSize[1]=2*firstlist...
  // the total size should now be (2+1+1)*firstlist
  BOOST_CHECK_EQUAL(TCell2.LiveEntries, 2 * numindex + 1);
  GR_index_t OldTotalSize = TCell2.TotalSize;

  TCell2.DeleteData(10);
  TCell2.DeleteData(11);
  TCell2.DeleteData(12);
  TCell2.DeleteData(numindex + 10);
  TCell2.DeleteData(numindex + 11);
  TCell2.DeleteData(numindex + 12);
  TCell2.UpdateEntries();
  BOOST_CHECK_EQUAL(OldTotalSize, TCell2.TotalSize); // we should still be in the same list so the number of totalsize should not be changed.
  TCell2.InsertingData(1);
  TCell2.UpdateEntries();
  BOOST_CHECK_EQUAL(TCell2.NEmpties, 5); // it checks if the getNewEntry has reused one of the already-vacated slots of memory...
  // or instead using a new one.

  TriCell* OpTC = TCell2.getEntry(numindex + 12);
  BOOST_CHECK(!OpTC->isDeleted());// Right now the way getNewEntry is working is Last In First Out (LIFO) but it should change..
  // to FIFO
  TCell2.clear();
  TCell2.UpdateEntries();
  //BOOST_CHECK(TCell2);
  BOOST_CHECK_MESSAGE(TCell2.LiveEntries == 0, "Clear has not done it's job");
  BOOST_CHECK_MESSAGE(TCell2.NEmpties == 0,
		      "there should be nothing in the empty vectors so "
		      "Clear has not done its job");
  BOOST_CHECK_MESSAGE(
      TCell2.TotalSize == numindex,
      "there should only be one list so the Clear has not done its job");

  TCell2.InsertingData(numindex + 1);
  TCell2.UpdateEntries();
  OldTotalSize = TCell2.TotalSize;
  GR_index_t OldLastEntry = TCell2.LastEntry;
  TCell2.DeleteData(numindex);
  TCell2.UpdateEntries();
  BOOST_CHECK_EQUAL(OldLastEntry, TCell2.LastEntry); // The last entry has been deleted but the value of the last entry should not change
  BOOST_CHECK_EQUAL(OldTotalSize, TCell2.TotalSize); // The total size should not change even there's no entity in the second list
  TriCell* OpTC1 = TCell2.getEntry(TCell2.LastEntry - 1);
  BOOST_CHECK(OpTC1->isDeleted());

  logMessage(3, "done\n");
}

// the following test case checks whether "pTNext", "pTFirst", "pTPrev", and "pTLast" are working or not.
// All these four functions just iterate amongst the entities of the container and point to the slots of ...
// contiguous memory in the m_data no matter whether they are empty, deleted or not.
BOOST_AUTO_TEST_CASE(NEIGHBORS)
{
  logMessage(3, "Testing pTFirst, pTLast, pTNext, ptPrev... ");

  // *********** Checking "pTFirst" and "pTLast"**************
  TestContainer TCellN;
  TCellN.UpdateEntries();
  TriCell *pTCStart = TCellN.First();
  TriCell *pTCEnd = TCellN.Last();
  BOOST_CHECK_MESSAGE(
      pTCStart==NULL,
      "When the container is empty pointer to the first slot of memory should refer to nothing");
  BOOST_CHECK_MESSAGE(
      pTCEnd==NULL,
      "When the container is empty pointer to the first slot of memory should refer to nothing");

  // ************ Checking "pTNext" *************************
  GR_index_t numindex = TCellN.first_size;
  TCellN.InsertingData(2 * numindex);
  TCellN.DeleteData(10);
  TCellN.DeleteData(numindex);
  TCellN.DeleteData(numindex + 10);
  TCellN.UpdateEntries();

  pTCEnd = TCellN.Last();
  pTCStart = TCellN.First();
  BOOST_CHECK(pTCEnd!=NULL);
  BOOST_CHECK(pTCStart!=NULL);
  TriCell *pTCCheck;
  GR_index_t ic = TCellN.getIndex(pTCStart);
  for (; pTCStart != pTCEnd && ic < TCellN.LastEntry;) {
    pTCCheck = TCellN.getEntry(ic);
    BOOST_CHECK_EQUAL(pTCCheck, pTCStart);

    ic++;
    pTCStart = TCellN.Next(pTCStart);
    BOOST_REQUIRE(ic <= TCellN.LastEntry);
  }

  // ************ Checking "pTPrev" **************
  pTCEnd = TCellN.Last();
  pTCStart = TCellN.First();
  ic = TCellN.getIndex(pTCEnd);
  for (; pTCEnd != pTCStart && ic > 0;) {
    pTCCheck = TCellN.getEntry(ic);
    BOOST_CHECK_EQUAL(pTCCheck, pTCEnd);
    ic--;
    pTCEnd = TCellN.Previous(pTCEnd);
  }

  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(Purge)
{
  logMessage(3, "Testing Purge and Resize... ");

  /* first test is for the case when there's only one list
   * so after purging and sorting the undeleted entites at
   * the end */
  TestContainer TCellP;
  TCellP.UpdateEntries();
  GR_index_t numindex = TCellP.first_size;
  TCellP.InsertingData(numindex);
  GR_index_t NumDelete = static_cast<GR_index_t>(numindex / 2.0); // number of entities that are set to be deleted.
  GR_index_t icP = 10; // some arbitrary number that shows where to start deleting
  for (GR_index_t id = 0; id < NumDelete; icP++, id++) {
    TCellP.DeleteData(icP);
  }
  TCellP.UpdateEntries();
  GR_index_t NumLives = TCellP.LiveEntries;
  GR_index_t OldTotal = TCellP.TotalSize;
  GR_index_t OldLastEntry = TCellP.LastEntry;

  TCellP.purge();
  TCellP.UpdateEntries();

  BOOST_CHECK_EQUAL(NumLives, TCellP.LiveEntries);
  BOOST_CHECK_EQUAL(TCellP.NEmpties, 0);
  BOOST_CHECK_EQUAL(OldTotal, TCellP.TotalSize); // if purging is working fine the number of TotalSize...
  //should not change because there's only one list.
  BOOST_CHECK_EQUAL(OldLastEntry, TCellP.LastEntry + NumDelete);// LastEntry=OldLastEntry-number-of-deleted-entities

  /* Second test is for the case that there are more than one list (here three)
   * and the number of deleted entries are more than m_totalSize - m_listSize[m_nLists - 1].
   * so, we should at least get rid of one list */

  TCellP.clear();
  TCellP.UpdateEntries();
  TCellP.InsertingData(3 * numindex);		// 3 lists should be created
  NumDelete = static_cast<GR_index_t>(1.5 * numindex);
  for (GR_index_t id = 0, icPl = 10; id < NumDelete; icPl++, id++)  // ic=100
      {
    TCellP.DeleteData(icPl); // we have deleted enough to be in the second list
  }
  TCellP.UpdateEntries();
  NumLives = TCellP.LiveEntries;
  OldTotal = TCellP.TotalSize;
  OldLastEntry = TCellP.LastEntry;

  BOOST_CHECK_MESSAGE(NumLives < 2 * numindex,
		      "you are not doing a good testing \n"
		      "NumLives<2*numindex otherwise resize would do nothing");

  TCellP.purge();// so by purging we should now get rid of the third list.
  TCellP.UpdateEntries();

  BOOST_CHECK_EQUAL(NumLives, TCellP.LiveEntries);
  BOOST_CHECK_EQUAL(TCellP.NEmpties, 0);
  BOOST_CHECK_EQUAL(TCellP.TotalSize, 2 * numindex);

  //******* check whether after purging the live entries are contiguous or not********

  TriCell *pTCStart = TCellP.First();
  TriCell *pTCEnd = TCellP.getEntry(TCellP.LastEntry - 1);

  for (; pTCStart != pTCEnd;) {
    BOOST_CHECK(!pTCStart->isDeleted());
    pTCStart = TCellP.Next(pTCStart);
  }

  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(OPERATORS)
{
  logMessage(3, "Testing Operators Overloading... ");

  TestContainer TCellO;
  TCellO.InsertingData(10);
  EntContainer<TriCell>::iterator TCiter1 = TCellO.Entity_begin(), TCiter2 =
      TCellO.Entity_begin(), TCiterEnd = TCellO.Entity_end();

  BOOST_CHECK(TCellO.OpEqu(TCiter1, TCiter2));
  BOOST_CHECK(TCellO.OpNEqu(TCiter1, TCiterEnd));

  EntContainer<TriCell>::iterator TCiterCheck1 = TCellO.OpPreInc(TCiter1);
  EntContainer<TriCell>::iterator TCiterCheck2 = TCellO.OpPostInc(TCiter2);

  BOOST_CHECK(TCellO.OpEqu(TCiterCheck1, TCiter1));
  BOOST_CHECK(TCellO.OpNEqu(TCiterCheck2, TCiter2));

  TriCell *pTCheck1 = TCellO.getEntry(0);
  TriCell *pTCheck2 = TCellO.getEntry(1);
  TriCell *pTCheck3 = &(*TCiter2);
  TriCell *pTCheck4 = &(*TCiterCheck2);

  BOOST_CHECK_EQUAL(pTCheck1, pTCheck4);
  BOOST_CHECK_EQUAL(pTCheck2, pTCheck3);

  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(RESERVE)
{
  logMessage(3, "Testing reserve function which is really trouble maker... ");
  TestContainer TCellR;
  TCellR.UpdateEntries();
  TCellR.InsertingData(3 * TCellR.first_size);
  TCellR.UpdateEntries();
  BOOST_CHECK_EQUAL(TCellR.TotalSize, 4 * TCellR.first_size);

  for (GR_index_t id = 0, ic = 5; id < 10; ic++, id++)  // ic=100
      {
    TCellR.DeleteData(ic); // we have deleted enough to be in the second list
  }
  TCellR.Reserve(3 * TCellR.first_size - 10);

  // Reserve for the case when the size is bigger than m_totalSize

  TCellR.clear();
  TCellR.UpdateEntries();
  TCellR.InsertingData(TCellR.first_size - 10);
  TCellR.Reserve(3 * TCellR.first_size);
  TCellR.UpdateEntries();
  BOOST_CHECK_EQUAL(TCellR.LiveEntries, 3 * TCellR.first_size);
  BOOST_CHECK_EQUAL(TCellR.ENum_Pointers, TCellR.first_size);
  logMessage(3, "done\n");
}

BOOST_AUTO_TEST_CASE(EMPTIES)
{
  logMessage(3, "Testing the compatibility of the entities of the empties... ");
  TestContainer TCellE;
  TCellE.InsertingData(20);
  TCellE.UpdateEntries();
  for (GR_index_t i = TCellE.LastEntry; i < TCellE.TotalSize; i++) {
    TriCell* pTE1 = TCellE.getEntryChecking(i);
    TriCell* pTE2 = TCellE.getEmEntry(TCellE.TotalSize - 1 - i);
    BOOST_CHECK_EQUAL(pTE1, pTE2);
  }
  TCellE.DeleteData(5);
  TCellE.DeleteData(19);
  TCellE.UpdateEntries();

  TriCell* pTEM1 = TCellE.getEmEntry(TCellE.TotalSize - 20);
  TriCell* pTEM2 = TCellE.getEmEntry(TCellE.TotalSize - 19);
  TriCell* pTEM3 = TCellE.getEntry(5);
  TriCell* pTEM4 = TCellE.getEntry(19);

  BOOST_CHECK_EQUAL(pTEM1, pTEM3);
  BOOST_CHECK_EQUAL(pTEM2, pTEM4);
  BOOST_CHECK_EQUAL(TCellE.ENum_Pointers, TCellE.TotalSize - 18);

  TCellE.purge();
  TCellE.UpdateEntries();
  BOOST_CHECK_EQUAL(TCellE.ENum_Pointers, TCellE.TotalSize - 18);
  for (GR_index_t i = TCellE.LastEntry; i < TCellE.TotalSize; i++) {
    TriCell* pTE1 = TCellE.getEntryChecking(i);
    TriCell* pTE2 = TCellE.getEmEntry(TCellE.TotalSize - 1 - i);
    BOOST_CHECK_EQUAL(pTE1, pTE2);
  }

  //	for (GR_index_t i=TCellE.ENum_Pointers ; i < TCellE.TotalSize ; i++)
  //	{
  //		TriCell* pTE2 = TCellE.getEmEntry(i);
  //		BOOST_CHECK(pTE2==NULL);
  //	}

  TCellE.DeleteData(10);
  pTEM1 = TCellE.getEmEntry(TCellE.TotalSize - 18);
  pTEM2 = TCellE.getEntry(10);
  BOOST_CHECK_EQUAL(pTEM1, pTEM2);
  BOOST_CHECK(TCellE.isDeleted(10));

  TCellE.InsertingData(1);
  pTEM2 = TCellE.getEntry(10);
  BOOST_CHECK(!TCellE.isDeleted(10));

  TCellE.InsertingData(10);
  TCellE.UpdateEntries();

  BOOST_CHECK_EQUAL(TCellE.ENum_Pointers, TCellE.TotalSize - 28);
  for (GR_index_t i = TCellE.LastEntry; i < TCellE.TotalSize; i++) {
    TriCell* pTEMl1 = TCellE.getEntryChecking(i);
    TriCell* pTEMl2 = TCellE.getEmEntry(TCellE.TotalSize - 1 - i);
    BOOST_CHECK_EQUAL(pTEMl1, pTEMl2);
  }

  logMessage(3, "Done.");
}
BOOST_AUTO_TEST_CASE(PARALLEL_Container)
{
//	logMessage(3,"Testing parallel insertion and deletion...");
//	if (NUM_PROCS > 1)
//	{
//		int num = FIRST_SIZE_PT+2;
//		TestContainer TCellP;
//		TCellP.InsertingData(20);
//		TCellP.UpdateEntries();
//		TCellP.SyncEntry();
//		GR_index_t Sizes[NUM_PROCS];
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//			Sizes[Id] = TCellP.m_TECTri.empties_pa[Id].size();
//		omp_set_dynamic(0);
//		omp_set_num_threads(NUM_PROCS);
//#pragma omp parallel shared(num,TCellP) 			//insertion without resizing
//		{
//			assert(omp_get_num_threads() == NUM_PROCS);
//			int ID = omp_get_thread_num();
//			TCellP.InsertingData(num,ID);
//#pragma omp barrier
//#pragma omp single
//			TCellP.UpdateEntries();
//		}
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//		{
//			BOOST_CHECK_EQUAL(Sizes[Id],TCellP.m_TECTri.empties_pa[Id].size()+(num-FIRST_SIZE_PT));
//			Sizes[Id] = TCellP.m_TECTri.empties_pa[Id].size();
//		}
//		omp_set_dynamic(0);
//		omp_set_num_threads(NUM_PROCS);
//#pragma omp parallel default(none) shared(num,TCellP)			//deletion...
//		{
//			//int ID = omp_get_thread_num();
//#pragma omp for schedule(static,2)
//			for (int ii=0 ; ii < NUM_PROCS*2 ; ii++ )
//			{
//				int ID = omp_get_thread_num();
//				TCellP.DeleteData(ii,ID);
//			}
//		}
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//		{
//			BOOST_CHECK_EQUAL(Sizes[Id],TCellP.m_TECTri.empties_pa[Id].size()-2);
//			Sizes[Id] = TCellP.m_TECTri.empties_pa[Id].size();
//		}
//		omp_set_dynamic(0);
//		omp_set_num_threads(NUM_PROCS);
//#pragma omp parallel default(none) shared(Sizes,TCellP)		//insertion with resize
//		{
//			int ID = omp_get_thread_num();
//			TCellP.InsertingData(Sizes[ID]+1,ID);
//		}
//
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//			BOOST_CHECK_EQUAL(FIRST_SIZE_PT-1,TCellP.m_TECTri.empties_pa[Id].size());
//
//		TCellP.SyncExit();
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//			BOOST_CHECK_EQUAL(0,TCellP.m_TECTri.empties_pa[Id].size());
//		GR_index_t checksize = 0;
//		BOOST_CHECK(TCellP.m_TECTri.isEmptiesUnique(checksize));
//		BOOST_CHECK_EQUAL(checksize,TCellP.m_TECTri.emptyavailable());
//
//		TCellP.InsertingData(400);
//		TCellP.UpdateEntries();
//		TCellP.SyncEntry();
//		for (int Id=0 ; Id < NUM_PROCS ; Id++)
//		{
//			BOOST_CHECK_EQUAL(FIRST_SIZE_PT,TCellP.m_TECTri.empties_pa[Id].size());
//			TCellP.IsListUnique(TCellP.m_TECTri.empties_pa[Id]);
//		}
//		GR_index_t num1 = 20;
//		omp_set_dynamic(0);
//		omp_set_num_threads(NUM_PROCS);
//#pragma omp parallel default(none) shared(TCellP,num1)
//		{
//#pragma omp for schedule(dynamic)
//			for (GR_index_t i=0 ; i < num1 ; i++)
//			{
//				int ID = omp_get_thread_num();
//				TriCell* pTC = TCellP.getEntry(i);
//				TCellP.likeReconfigure(pTC,ID);
//			}
//#pragma omp barrier
//#pragma omp single
//			{
//				for (int Id=0 ; Id < NUM_PROCS ; Id++)
//					TCellP.IsListUnique(TCellP.m_TECTri.empties_pa[Id]);
//			}
//		}
//		TCellP.SyncExit();
//		checksize = 0;
//		BOOST_CHECK(TCellP.m_TECTri.isEmptiesUnique(checksize));
//		BOOST_CHECK_EQUAL(checksize,TCellP.m_TECTri.emptyavailable());
//	}
//
//	logMessage(3," done\n");
}

