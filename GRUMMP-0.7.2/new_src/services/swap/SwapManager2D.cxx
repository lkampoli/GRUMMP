#include "GR_FaceSwapInfo.h"
#include "GR_Geometry.h"
#include "GR_SwapManager.h"
#include "GR_SwapDecider.h"
#include "GR_Mesh2D.h"
#include <omp.h>

#ifdef SIM_ANNEAL_TEST
extern double simAnnealTemp;
#endif

namespace GRUMMP
{
  void
  FaceSwapInfo2D::calculateFaceInfo(Face* pF)
  {
    m_nonSwappable = false;
    if (!pF->isSwapAllowed() /*|| pF->isLocked()*/) {
      m_nonSwappable = true;
      return;
    }

    //@@ Case: one or both cells is not a tri, including boundaries
    Cell *pCLeft = pF->getLeftCell();
    Cell *pCRight = pF->getRightCell();
    if (!pCLeft->isValid() || !pCRight->isValid()
	|| pCLeft->getType() != Cell::eTriCell
	|| pCRight->getType() != Cell::eTriCell) {
      m_nonSwappable = true;
      return;
    }

    m_pTris[0] = dynamic_cast<TriCell*>(pCLeft);
    m_pTris[1] = dynamic_cast<TriCell*>(pCRight);
    // Already checked that these are tri's, so the dynamic casts should
    // always work unless the database is corrupt somehow.
    assert(m_pTris[0] != NULL);
    assert(m_pTris[1] != NULL);

    m_pVerts[0] = pF->getVert(0);
    m_pVerts[1] = pF->getVert(1);
    m_pVerts[2] = m_pTris[0]->getOppositeVert(pF);
    m_pVerts[3] = m_pTris[1]->getOppositeVert(pF);
    if (m_pVerts[0]->getSpaceDimen() == 2) {
      int iOrientA = checkOrient2D(m_pVerts[0], m_pVerts[3], m_pVerts[2]);
      int iOrientB = checkOrient2D(m_pVerts[1], m_pVerts[2], m_pVerts[3]);

      m_orientOK = (iOrientA == 1 && iOrientB == 1);
    }
    else {    // TODO:ASK CARL HOW DO TO THIS
      double unitNormalA[3], unitNormalB[3];
      calcUnitNormal(m_pVerts[0]->getCoords(), m_pVerts[3]->getCoords(),
		     m_pVerts[2]->getCoords(), unitNormalA);
      calcUnitNormal(m_pVerts[1]->getCoords(), m_pVerts[2]->getCoords(),
		     m_pVerts[3]->getCoords(), unitNormalB);
      // take the dot of the unit normals, assert < 90 degrees.
      double dot = dDOT3D(unitNormalA, unitNormalB);
      if (iFuzzyComp(dot, 0.) > -1)
	m_orientOK = true;
      else
	m_orientOK = false;

    }
  }

  inline void
  FaceSwapInfo2D::setTris(Face* pF)
  {
    m_pTris[0] = dynamic_cast<TriCell*>(pF->getLeftCell());
    m_pTris[1] = dynamic_cast<TriCell*>(pF->getRightCell());

    m_pVerts[0] = pF->getVert(0);
    m_pVerts[1] = pF->getVert(1);
    m_pVerts[2] = m_pTris[0]->getOppositeVert(pF);
    m_pVerts[3] = m_pTris[1]->getOppositeVert(pF);

  }

  SwapManager2D::~SwapManager2D()
  {
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      m_successes += m_successesThreads[ID][0];
      m_attempts += m_attemptsThreads[ID][0];
    }

    m_pM2D->removeObserver(this);
    logMessage(2, "SwapManager exiting.\n");
    if (m_pSD)
      logMessage(2, "  Used a %s SwapDecider.\n", m_pSD->getName().c_str());
    logMessage(2, "  Face swapping stats: \n");
    logMessage(2, "    %6d swaps out of %6d attempts\n", m_successes,
	       m_attempts);
  }
// This function is solely used for the serial version.
  int
  SwapManager2D::swapFace(FaceSwapInfo2D& FSI2D)
  {
    assert(omp_in_parallel() == 0);
    m_attempts++;
    if (m_pSD->doFaceSwap(FSI2D)) {
      reconfigure(FSI2D);
      m_pM2D->sendEvents();
      m_successes++;
#ifndef NDEBUG
      m_pM2D->writeTempMesh();
#endif
      return 1;
    }
    else
      return (0);
  }

  template<class myType>
    void
    SwapManager2D::LoadBalancing(std::vector<myType> VecToBalanced[],
				 int num_procs)
    {
#if NUM_PROCS>1
      int VecSize[num_procs];
      int TotSize = 0;
      for (int ID = 0; ID < num_procs; ID++) {
	VecSize[ID] = VecToBalanced[ID].size();
	TotSize += VecSize[ID];
      }
      int AveSize = int(floor(double(TotSize) / num_procs));
      for (int i = 0; i < num_procs - 1; i++) {
	int *MaxPos = std::max_element(VecSize, VecSize + num_procs);
	int *MinPos = std::min_element(VecSize, VecSize + num_procs);
	typename std::vector<myType>::iterator it = VecToBalanced[MaxPos
	    - VecSize].begin();
	std::advance(it, AveSize);
	VecToBalanced[MinPos - VecSize].insert(
	    VecToBalanced[MinPos - VecSize].end(), it,
	    VecToBalanced[MaxPos - VecSize].end());
	VecToBalanced[MaxPos - VecSize].erase(
	    it, VecToBalanced[MaxPos - VecSize].end());
	*MinPos += *MaxPos - AveSize;
	*MaxPos = AveSize;
      }
#endif
    }
//****Defining new function for extracting conflicting faces
//****and creating new face queue without conflicts ready to parallel
  double
  SwapManager2D::faceQueueManager(int passes, bool SwapAll)
  {
    double paralleltime;
    double FaceQueue_Time = omp_get_wtime();
    if (!SwapAll) {
      SwapAllowedFaceQueue.clear();
      NConflictFaceQueue.clear();
      std::vector<Face*> SwapAllowedFaceQueueLists[NUM_PROCS];
      struct IteratorPairs {
	std::set<Face*>::iterator first;
	std::set<Face*>::iterator second;
      };
      GR_index_t faceQueueSize = faceQueue.size();
      IteratorPairs chunks[NUM_PROCS];
      GR_index_t chunk_size = GR_index_t(faceQueueSize / double(NUM_PROCS));
      std::set<Face*>::iterator iter_begin = faceQueue.begin();
      for (int i = 0; i < NUM_PROCS - 1; i++) {
	std::set<Face*>::iterator iter_end = iter_begin;
	for (GR_index_t k = 0; k < chunk_size; k++)
	  ++iter_begin;
	chunks[i].first = iter_end;
	chunks[i].second = iter_begin;
	SwapAllowedFaceQueueLists[i].reserve(
	    GR_index_t(faceQueueSize / double(NUM_PROCS)) + 1);
      }
      chunks[NUM_PROCS - 1].first = iter_begin;
      chunks[NUM_PROCS - 1].second = faceQueue.end();
      SwapAllowedFaceQueueLists[NUM_PROCS - 1].reserve(
	  GR_index_t(faceQueueSize / double(NUM_PROCS)) + 1);

      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel default(none) shared(SwapAllowedFaceQueueLists,faceQueueSize,chunks)// private(FSI2D)
      {
	int ID = omp_get_thread_num();
	FaceSwapInfo2D FSI2D;
	for (std::set<Face*>::iterator it = chunks[ID].first;
	    it != chunks[ID].second; ++it) {
	  Face *pF = *it;
	  FSI2D.calculateFaceInfo(pF);
	  if (FSI2D.isSwappable() && m_pSD->doFaceSwap(FSI2D))
	    SwapAllowedFaceQueueLists[ID].push_back(pF);
	}
      }

      double ConflictReTime = omp_get_wtime();
      std::vector<Face*>::iterator iterNew = NConflictFaceQueue.begin(),
	  iterSoFar = NConflictFaceQueue.end();
      std::set<Cell*> CellLists;
      std::set<Face*> FaceLists;
      std::set<Vert*> VertLists;

      paralleltime = omp_get_wtime();
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(iterSoFar,iterNew,SwapAllowedFaceQueueLists) private(CellLists,FaceLists,VertLists)
      {
	int IDT = omp_get_thread_num();
	for (GR_index_t iterCount = 0;
	    iterCount < SwapAllowedFaceQueueLists[IDT].size(); iterCount++) {
	  FaceLists.clear();
	  CellLists.clear();
	  VertLists.clear();
	  Face *pF = (SwapAllowedFaceQueueLists[IDT][iterCount]);

	  //************ finding 4 vertexes of Quad surrounding the corresponded Face
	  Cell *pCLeft = pF->getLeftCell();
	  Cell *pCRight = pF->getRightCell();
	  if (pCRight->getType() != CellSkel::eTriCell)
	    continue;
	  for (int i = pCRight->getNumVerts() - 1; i >= 0; i--) {
	    Vert* pV1 = pCRight->getVert(i);
	    assert(pV1->isValid());
	    VertLists.insert(pV1);
	  }
	  for (int i = pCLeft->getNumVerts() - 1; i >= 0; i--) {
	    Vert* pV2 = pCLeft->getVert(i);
	    assert(pV2->isValid());
	    VertLists.insert(pV2);

	  }
	  for (std::set<Vert*>::iterator NV = VertLists.begin();
	      NV != VertLists.end(); ++NV) {
	    for (int NF = 0; NF < (*NV)->getNumFaces(); ++NF) {
	      Face* VPF = (*NV)->getFace(NF);
	      for (int NC = 0; NC < VPF->getNumCells(); NC++) {
		Cell* VPC = VPF->getCell(NC);
		CellLists.insert(VPC);
	      }
	    }
	  }

	  for (std::set<Cell*>::iterator itC = CellLists.begin();
	      itC != CellLists.end(); ++itC) {
	    Cell* NewC = *itC;
	    for (int NumF = 0; NumF < NewC->getNumFaces(); ++NumF) {
	      Face* VCF = NewC->getFace(NumF);
	      FaceLists.insert(VCF);
	    }
	  }
	  bool no_conflict = false;
	  double Critical_Time = omp_get_wtime();
#pragma omp critical(NConflictFaceQueueinsert1)
	  {
	    for (std::set<Face*>::iterator it2 = FaceLists.begin();
		it2 != FaceLists.end(); ++it2) {
	      for (std::vector<Face*>::iterator iterNew1 = iterNew;
		  iterNew1 != iterSoFar; ++iterNew1) {
		if ((*iterNew1) == (*it2)) {
		  no_conflict = true;
		  break;
		}
	      }
	      if (no_conflict)
		break;
	      else
		continue;
	    }
	  } //****End of Critical Section
	  Critical_Time = omp_get_wtime() - Critical_Time;
	  Time_Spent_in_Critical += Critical_Time;
	  if (!no_conflict) {
	    NConflictFaceQueue.push_back(pF);
	    iterNew = NConflictFaceQueue.begin();
	    iterSoFar = NConflictFaceQueue.end();
	  }
	  else
	    RemainingFaces[IDT][0].push_back(pF);
	  //*************** Done with conflicts removal *******
	} //******* End of for *****
      } //**** End of parallel Section
      paralleltime = omp_get_wtime() - paralleltime;
      Parallel_Time += paralleltime;
      ConflictReTime = omp_get_wtime() - ConflictReTime;
      RemovingConflicts += ConflictReTime;
    }

    //************ for the cases that "swapAllFaces" has been called *******
    logMessage(3, "passes: %d\n", passes);
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      FaceStructs[ID].clear();
      RemainingFaces[ID][0].clear();
    }
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_PROCS);
#pragma omp parallel
    {
      GR_index_t attemptsThreads = 0;
      GR_index_t successesThreads = 0;
      double SwappingCheck = omp_get_wtime();
      // *** these are timing for critical and non-critical sections ***
//    	double CTSum = 0;
//    	double NCTSum = 0;
//    	double CT = 0;
//    	double NCT = 0;
//    	NCT = omp_get_wtime();
      bool IsInUse = false;
      int ID = omp_get_thread_num();
      FaceSwapInfo2D FS2D;

      for (GR_index_t fN = 0; fN < faceQueueVecs[ID].size(); ++fN) {
	Face *pF = faceQueueVecs[ID][fN];
	FS2D.calculateFaceInfo(pF);
	attemptsThreads++;
	if (!FS2D.isSwappable() || !m_pSD->doFaceSwap(FS2D)) {
	  continue;
	}
	else
	  RemainingFaces[ID][0].push_back(pF);
      }
//      NCT = omp_get_wtime() - NCT;
//      NCTSum += NCT;
      faceQueueVecs[ID].clear();
      logMessage(3, "Number of faces in thread %d is %lu \n", ID,
		 RemainingFaces[ID][0].size());
//      CT = omp_get_wtime();
      for (GR_index_t fN = 0; fN < RemainingFaces[ID][0].size(); ++fN) {
	IsInUse = false;
	Face *pF = RemainingFaces[ID][0][fN];
	Vert* pV[4];
	FS2D.setTris(pF);
	for (int i = 0; i < 4 && !IsInUse; i++)
	  pV[i] = FS2D.getVert(i);
#pragma omp critical(CheckingIfVertIsInUseFQM)
	{
	  if (!IsInUse)
	    pV[0]->isInUseWOC(IsInUse);
	  if (!IsInUse)
	    pV[1]->isInUseWOC(IsInUse);
	  if (!IsInUse)
	    pV[2]->isInUseWOC(IsInUse);
	  if (!IsInUse)
	    pV[3]->isInUseWOC(IsInUse);
	  if (!IsInUse) {
	    for (int i = 0; i < 4; i++)
	      pV[i]->setInUseWOC(true);
	  }
	} // end of critical

	if (IsInUse) {
	  Face* pFif2 = FS2D.getFace();
	  faceQueueVecs[ID].push_back(pFif2);
	}
	else {
	  FaceStructs[ID].push_back(FS2D);
	  m_successesThreads[ID][0]++;
	  successesThreads++;
	}
      } // end of for
//      CT = omp_get_wtime() - CT;
//      CTSum += CT;
      RemainingFaces[ID][0].clear();
      RemainingFaces[ID][0].swap(faceQueueVecs[ID]);
//      CriticalTime[ID][0] += CTSum;
//      SwappingCheckTime[ID][0] += NCTSum;
      SwappingCheck = omp_get_wtime() - SwappingCheck;
      FQMThread[ID][0] += SwappingCheck;
      m_attemptsThreads[ID][0] += attemptsThreads;
      m_successesThreads[ID][0] += successesThreads;
    } // end of parallel

    FaceQueue_Time = omp_get_wtime() - FaceQueue_Time;
    return FaceQueue_Time;
  } //****** End of Function*****

  int
  SwapManager2D::swapAllFaces(bool OMP_Flag)
  {
    double SwappingStartTime = omp_get_wtime();
    if (!OMP_Flag) {
      // Initialize the queue with all the faces in the mesh.
      EntContainer<EdgeFace>::iterator ECiter = m_pM2D->edgeFace_begin(),
	  ECiterEnd = m_pM2D->edgeFace_end();
      for (; ECiter != ECiterEnd; ++ECiter) {
	Face *pF = &(*ECiter);
	faceQueue.insert(pF);
      }
      int NumTotSwaps = swapAllQueuedFacesRecursive();
      SwappingStartTime = omp_get_wtime() - SwappingStartTime;
      TotalSwappingTime += SwappingStartTime;
      return NumTotSwaps;
    }
    else {
      //     Initialize the queue with all the faces in the mesh.
      double swapAllFacess = omp_get_wtime();
      double paralleltime = omp_get_wtime();
      GR_index_t chunk_size = GR_index_t(
	  double(m_pM2D->getNumFaces()) / NUM_PROCS);
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(chunk_size)
      {
	EntContainer<EdgeFace>::iterator ECiter = m_pM2D->edgeFace_begin();
	EntContainer<EdgeFace>::iterator ECiterEnd = m_pM2D->edgeFace_begin();

	int ID = omp_get_thread_num();
	if (ID < NUM_PROCS - 1) {
	  ECiter.jump(chunk_size * ID);
	  ECiterEnd.jump(chunk_size * (ID + 1));

	}
	else {
	  EntContainer<EdgeFace>::iterator ECiterEndBuff =
	      m_pM2D->edgeFace_end();

	  ECiter.jump(chunk_size * ID);
	  ECiterEnd = ECiterEndBuff;

	}
	for (; ECiter != ECiterEnd; ++ECiter) {
	  Face *pF = &(*ECiter);
	  faceQueueVecs[ID].push_back(pF);
	}
      }
      paralleltime = omp_get_wtime() - paralleltime;
      Parallel_Time += paralleltime;
      swapAllFacess = omp_get_wtime() - swapAllFacess;
      swapAllFacess_Time += swapAllFacess;

      paralleltime = 0;
      Sync_Entry_Time += m_pM2D->Mesh2D_sync_entry(paralleltime, true);
//      Sync_Entry_Time += paralleltime;
      Parallel_Time += paralleltime;

      int NumSwaps = swapAllQueuedFacesRecursive(OMP_Flag, true);

      Sync_Exit_Time += m_pM2D->Mesh2D_sync_exit(true);

      for (int i = 0; i < NUM_PROCS; i++)
	faceQueueVecs[i].clear();

      SwappingStartTime = omp_get_wtime() - SwappingStartTime;
      TotalSwappingTime += SwappingStartTime;

      return NumSwaps;
    }
  }

  int
  SwapManager2D::swapAllQueuedFacesRecursive(bool OMP_Flag, bool SwapAll)
  {
    if (!OMP_Flag) {
      int numSwaps = 0, swapsThisPass = 0, passes = 0;
      do {
	int attempts = faceQueue.size();
	swapsThisPass = swapAllQueuedFaces();
	numSwaps += swapsThisPass;
	passes++;

	logMessage(2, "Pass %d, %d swaps out of %d attempts\n", passes,
		   swapsThisPass, attempts);
	m_pM2D->sendEvents();
	logMessage(2, "  Master swap queue now has %zu faces.\n",
		   faceQueue.size());
#ifdef SIM_ANNEAL_TEST
	logMessage(2, "  Sim anneal temp is %f\n", simAnnealTemp);
	simAnnealTemp /= 1.1;
      }while (swapsThisPass != 0 && simAnnealTemp > 1.e-6);
#else
      }
      while (swapsThisPass != 0);

#endif
      return numSwaps;
    }
    else {
      int RFZ = 0;
      assert(OMP_Flag);
      parallelswapping = OMP_Flag;
      m_pM2D->IsInParallel(true);
      double SendTime = omp_get_wtime();
      if (!SwapAll)
	m_pM2D->sendEvents();
      SendTime = omp_get_wtime() - SendTime;
      SendEventsTime += SendTime;
      int numSwaps = 0, swapsThisPass = 0, passes = 0;
      do {
	RFZ = 0;
	LoadBalancing(faceQueueVecs, NUM_PROCS);
	Quary_Time += faceQueueManager(passes, SwapAll);
	LoadBalancing(FaceStructs, NUM_PROCS);
	swapsThisPass = swapAllQueuedFaces(OMP_Flag);
	numSwaps += swapsThisPass;

	logMessage(2, "   %d swaps at pass %d.\n", swapsThisPass, passes);

	if (!SwapAll) {
	  faceQueue.clear();
	  for (int ID = 0; ID < NUM_PROCS; ID++) {
	    for (std::vector<Face*>::iterator it =
		RemainingFaces[ID][0].begin();
		it != RemainingFaces[ID][0].end(); ++it)
	      faceQueue.insert(*it);
	    RemainingFaces[ID][0].clear();
	  }
	}
	else {
	  omp_set_dynamic(0);
	  omp_set_num_threads(NUM_PROCS);
#pragma omp parallel
	  {
	    int ID = omp_get_thread_num();
	    faceQueueVecs[ID].clear();
	    for (GR_index_t fB = 0; fB < RemainingFaces[ID][0].size(); ++fB) {
	      Face* pF = RemainingFaces[ID][0][fB];
	      if (pF->isValid() && !pF->isDeleted()) {
		faceQueueVecs[ID].push_back(pF);
	      }
	    }
	    RemainingFaces[ID][0].clear();
	  }
	}

	SendTime = omp_get_wtime();
	m_pM2D->sendEvents();
	SendTime = omp_get_wtime() - SendTime;
	SendEventsTime += SendTime;

	for (GR_index_t ID = 0; ID < NUM_PROCS; ID++)
	  RFZ += faceQueueVecs[ID].size();
	logMessage(3, "   Number of remaining faces is: %d.\n", RFZ);

	for (GR_index_t ID = 0; ID < NUM_PROCS; ID++)
	  logMessage(3, "   faceQueueVecs[%u] now has %zu faces.\n", ID,
		     faceQueueVecs[ID].size());

	passes++;

#ifdef SIM_ANNEAL_TEST
	logMessage(2, "  Sim anneal temp is %f\n", simAnnealTemp);
	simAnnealTemp /= 1.1;
      }while (swapsThisPass != 0 && simAnnealTemp > 1.e-6);
#else
      }
      while (swapsThisPass != 0);
#endif

      m_pM2D->IsInParallel(false);
      return numSwaps;
    }
  }

  int
  SwapManager2D::swapAllQueuedFaces(bool OMP_Flag)
  {
    parallelswapping = OMP_Flag;
    if (!OMP_Flag) {
      m_pM2D->sendEvents();
      std::set<Face*> liveFaceQueue;
      // Initialize the queue with all the faces in the queue.
      swap(faceQueue, liveFaceQueue);
      logMessage(2, "  Master swap queue now has %zu faces.\n",
		 faceQueue.size());
      std::set<Face*>::iterator iter = liveFaceQueue.begin(), iterEnd =
	  liveFaceQueue.end();
      int numLocalSwaps = 0;
      FaceSwapInfo2D FSI2D;
      for (; iter != iterEnd; ++iter) {
	Face *pF = *iter;
	if (!pF)
	  continue;
//        if(!pF->doFullCheck()) continue;
//        assert(pF->doFullCheck());
//        assert(m_pM2D->isValid());
	FSI2D.calculateFaceInfo(pF);
	numLocalSwaps += swapFace(FSI2D);
      } // for
      return numLocalSwaps;
    }
    else {
      int numTotalSwaps = 0;
      double paralleltime = 0;
      int numSwapslocal[NUM_PROCS][20] =
	{
	  { 0 } };
      paralleltime = omp_get_wtime();
      double Time_SwapAll = omp_get_wtime();

      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel default(none) shared(numSwapslocal)
      {
	double ThreadReconTime = omp_get_wtime();
	int ID = omp_get_thread_num();
	GR_index_t FSSize = FaceStructs[ID].size();
	numSwapsG[ID][0] += FSSize;
	numSwapslocal[ID][0] += FSSize;
	for (GR_index_t fN = 0; fN < FSSize; ++fN)
	  reconfigure(FaceStructs[ID][fN], ID, true);
	ThreadReconTime = omp_get_wtime() - ThreadReconTime;
	ReconTime[ID][0] += ThreadReconTime;
      } //end of parallel region

      Time_SwapAll = omp_get_wtime() - Time_SwapAll;
      swapAllQueuedFacesTime += Time_SwapAll;
      paralleltime = omp_get_wtime() - paralleltime;
      Parallel_Time += paralleltime;
      for (int ID = 0; ID < NUM_PROCS; ID++)
	numTotalSwaps += numSwapslocal[ID][0];
      return numTotalSwaps;
    }
  }

  bool
  SwapManager2D::reconfigure(FaceSwapInfo2D& FSI2D, Cell* pTri[], Cell *pC0,
			     Cell *pC1, Vert *pV0, Vert *pV1, Vert *pV2,
			     Vert *pV3) const
  {
    FSI2D.getTris(pTri);
    int iReg = pTri[0]->getRegion();
    int ID = omp_get_thread_num();
    assert(omp_get_num_threads() == NUM_PROCS);

    m_pM2D->deleteCellRe(pTri[0], ID);
    m_pM2D->deleteCellRe(pTri[1], ID);

    pV0 = FSI2D.getVertA();
    pV1 = FSI2D.getVertB();
    pV2 = FSI2D.getVertC();
    pV3 = FSI2D.getVertD();

    pC0 = m_pM2D->createTriCell(pV3, pV2, pV0, iReg, true, ID);
    pC1 = m_pM2D->createTriCell(pV2, pV3, pV1, iReg, true, ID);

    assert(pC0->doFullCheck());
    assert(pC1->doFullCheck());

    return true;
  }

  bool
  SwapManager2D::reconfigure(FaceSwapInfo2D& FSI2D, int ID, bool OMP_Flag) const
  {
    if (!OMP_Flag) {
      Cell* pTri[2];
      FSI2D.getTris(pTri);
      int iReg = pTri[0]->getRegion();
      m_pM2D->deleteCell(pTri[0]);
      m_pM2D->deleteCell(pTri[1]);

      Cell *pC0 = m_pM2D->createTriCell(FSI2D.getVertD(), FSI2D.getVertC(),
					FSI2D.getVertA(), iReg, true);
      Cell *pC1 = m_pM2D->createTriCell(FSI2D.getVertC(), FSI2D.getVertD(),
					FSI2D.getVertB(), iReg, true);

      assert(pC0->doFullCheck());
      assert(pC1->doFullCheck());

      return true;
    }
    else {
      Cell* pTri[2];
      FSI2D.getTris(pTri);
      int iReg = pTri[0]->getRegion();

      m_pM2D->deleteCellRe(pTri[0], ID);
      m_pM2D->deleteCellRe(pTri[1], ID);

      Cell *pC0 = m_pM2D->createTriCell(FSI2D.getVertD(), FSI2D.getVertC(),
					FSI2D.getVertA(), iReg, false, ID);
      Cell *pC1 = m_pM2D->createTriCell(FSI2D.getVertC(), FSI2D.getVertD(),
					FSI2D.getVertB(), iReg, false, ID);

      assert(pC0->doFullCheck());
      assert(pC1->doFullCheck());

      for (GR_index_t i = 0; i < 4; i++) {
	Vert* pV = FSI2D.getVert(i);
	pV->setInUseWOC(false);
      }

      return true;

    }
  }

  void
  SwapManager2D::receiveDeletedFaces(std::vector<Face*>& deletedFaces)
  {
    if (!parallelswapping) {
      logMessage(4, "SwapManager2D received word about %zu deleted faces.\n",
		 deletedFaces.size());
      std::vector<Face*>::iterator iter = deletedFaces.begin(), iterEnd =
	  deletedFaces.end();
      for (; iter != iterEnd; ++iter) {
	Face *pF = *iter;
	faceQueue.erase(pF);
      }
    }
  }

  void
  SwapManager2D::receiveCreatedCells(std::vector<Cell*>& createdCells)
  {
    double CCTime = omp_get_wtime();
    logMessage(4, "SwapManager2D received word about %zu new cells.\n",
	       createdCells.size());
    if (!parallelswapping) {
      std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
	  createdCells.end();
      for (; iter != iterEnd; ++iter) {
	Cell *pC = *iter;
	if (pC->getType() != CellSkel::eTriCell)
	  continue;
	for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pC->getFace(i);
	  faceQueue.insert(pF);
	}
      }
    }
    else {
#pragma omp parallel
      {
	int ID = omp_get_thread_num();
#pragma omp for schedule(guided)
	for (GR_index_t iter = 0; iter < createdCells.size(); ++iter) {
	  Cell *pC = createdCells[iter];
	  if (pC->getType() != CellSkel::eTriCell)
	    continue;
	  for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
	    Face *pF = pC->getFace(i);
	    pF->setBadFaceWOC(false);
	    faceQueueVecs[ID].push_back(pF);
	  }
	}
      } // end of parallel
    }
    CCTime = omp_get_wtime() - CCTime;
    CreatedCells += CCTime;
  }

  void
  SwapManager2D::receiveMovedVerts(std::vector<Vert*>& movedVerts)
  {
    logMessage(4, "SwapManager2D received word about %zu moved verts.\n",
	       movedVerts.size());
    std::vector<Vert*>::iterator iter = movedVerts.begin(), iterEnd =
	movedVerts.end();
    int counter = 0;
    if (!parallelswapping) {
      for (; iter != iterEnd; ++iter) {
	Vert *pV = *iter;
	// Need to be sure to get not just incident faces but also
	// opposite faces.
	for (int i = pV->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pV->getFace(i);
	  Cell *pC = pF->getLeftCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--)
	    faceQueue.insert(pC->getFace(ii));

	  pC = pF->getRightCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--)
	    faceQueue.insert(pC->getFace(ii));
	}
      }
    }
    else {
      for (; iter != iterEnd; ++iter) {
	Vert *pV = *iter;
	// Need to be sure to get not just incident faces but also
	// opposite faces.
	for (int i = pV->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pV->getFace(i);
	  Cell *pC = pF->getLeftCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    int ID = counter % NUM_PROCS;
	    counter++;
	    pC->getFace(ii)->setBadFace(false);
	    faceQueueVecs[ID].push_back(pC->getFace(ii));
	  }

	  pC = pF->getRightCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    int ID = counter % NUM_PROCS;
	    counter++;
	    pC->getFace(ii)->setBadFace(false);
	    faceQueueVecs[ID].push_back(pC->getFace(ii));
	  }
	}
      }
    } // End of if
  }

  void
  FaceSwapInfo2D::getTris(Cell* pTris[2]) const
  {
    pTris[0] = m_pTris[0];
    pTris[1] = m_pTris[1];
  }

  bool
  SwapManagerSurf::reconfigure(FaceSwapInfo2D& FSI2D, Cell* pTri[], Cell *pC0,
			       Cell *pC1, Vert *pV0, Vert *pV1, Vert *pV2,
			       Vert *pV3) const
  {
    FSI2D.getTris(pTri);
    int iReg = pTri[0]->getRegion();
    int ID = omp_get_thread_num();
    assert(omp_get_num_threads() == NUM_PROCS);

    RefFace * surface0 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[0]);
    RefFace * surface1 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[1]);
    // going to do the check here, since it is the easiest
    if (surface0 != surface1)
      return false; // shouldn't happen
    m_pM2D->deleteCellRe(pTri[0], ID);
    m_pM2D->deleteCellRe(pTri[1], ID);

    pV0 = FSI2D.getVertA();
    pV1 = FSI2D.getVertB();
    pV2 = FSI2D.getVertC();
    pV3 = FSI2D.getVertD();

    pC0 = m_pM2D->createTriCell(pV3, pV2, pV0, iReg, false, ID);
    pC1 = m_pM2D->createTriCell(pV2, pV3, pV1, iReg, false, ID);

    assert(pC0->doFullCheck());
    assert(pC1->doFullCheck());

    dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC0, surface0);
    dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC1, surface1);
    return true;
  }

  bool
  SwapManagerSurf::reconfigure(FaceSwapInfo2D& FSI2D, int ID,
			       bool OMP_Flag) const
  {
    if (!OMP_Flag) {
      Cell* pTri[2];
      FSI2D.getTris(pTri);
      int iReg = pTri[0]->getRegion();
      RefFace * surface0 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[0]);
      RefFace * surface1 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[1]);

      m_pM2D->deleteCell(pTri[0]);
      m_pM2D->deleteCell(pTri[1]);

      Cell *pC0 = m_pM2D->createTriCell(FSI2D.getVertD(), FSI2D.getVertC(),
					FSI2D.getVertA(), iReg);
      Cell *pC1 = m_pM2D->createTriCell(FSI2D.getVertC(), FSI2D.getVertD(),
					FSI2D.getVertB(), iReg);

      dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC0, surface0);
      dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC1, surface1);
      assert(pC0->doFullCheck());
      assert(pC1->doFullCheck());

      return true;
    }
    else {
      Cell* pTri[2];
      FSI2D.getTris(pTri);
      int iReg = pTri[0]->getRegion();
      assert(ID == omp_get_thread_num());
      assert(omp_get_num_threads() == NUM_PROCS);
      RefFace * surface0 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[0]);
      RefFace * surface1 = dynamic_cast<SurfMesh*>(m_pM2D)->getSurface(pTri[1]);
      m_pM2D->deleteCellRe(pTri[0], ID);
      m_pM2D->deleteCellRe(pTri[1], ID);

      Cell *pC0 = m_pM2D->createTriCell(FSI2D.getVertD(), FSI2D.getVertC(),
					FSI2D.getVertA(), iReg, false, ID);
      Cell *pC1 = m_pM2D->createTriCell(FSI2D.getVertC(), FSI2D.getVertD(),
					FSI2D.getVertB(), iReg, false, ID);

      assert(pC0->doFullCheck());
      assert(pC1->doFullCheck());
      dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC0, surface0);
      dynamic_cast<SurfMesh*>(m_pM2D)->setSurface(pC1, surface1);
      for (GR_index_t i = 0; i < 4; i++) {
	Vert* pV = FSI2D.getVert(i);
	pV->setInUse(false);
      }

      return true;

    }
  }

} // namespace GRUMMP
