#include "GR_EdgeCanon.h"
#include "GR_EdgeSwapInfo.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Geometry.h"
#include "GR_SwapManager.h"
#include "GR_SwapDecider.h"
#include "GR_VolMesh.h"
#include <omp.h>
#include <algorithm>

#ifdef SIM_ANNEAL_TEST
extern double simAnnealTemp;
#endif

namespace GRUMMP
{

  SwapManager3D::~SwapManager3D()
  {
    for (int i = 0; i < 11; i++) {
      for (int ID = 0; ID < NUM_PROCS; ++ID) {
	m_edgeAttempts[i] += m_edgeAttemptsThreads[i][ID];
	m_edgeSuccesses[i] += m_edgeSuccessesThreads[i][ID];
	m_faceAttempts[i] += m_faceAttemptsThreads[i][ID];
	m_faceSuccesses[i] += m_faceSuccessesThreads[i][ID];
      }
    }
    if (m_pVM)
      m_pVM->removeObserver(this);
    logMessage(2, "SwapManager exiting.\n");
    logMessage(2, "  Used a %s SwapDecider.\n", m_pSD->getName().c_str());
    logMessage(2, "  Face swapping stats: \n");
    logMessage(2, "    T22: %6d swaps out of %6d attempts\n",
	       m_faceSuccesses[FaceSwapInfo3D::eT22],
	       m_faceAttempts[FaceSwapInfo3D::eT22]);
    logMessage(2, "    T23: %6d swaps out of %6d attempts\n",
	       m_faceSuccesses[FaceSwapInfo3D::eT23],
	       m_faceAttempts[FaceSwapInfo3D::eT23]);
    logMessage(2, "    T32: %6d swaps out of %6d attempts\n",
	       m_faceSuccesses[FaceSwapInfo3D::eT32],
	       m_faceAttempts[FaceSwapInfo3D::eT32]);
    logMessage(2, "    T44: %6d swaps out of %6d attempts\n",
	       m_faceSuccesses[FaceSwapInfo3D::eT44],
	       m_faceAttempts[FaceSwapInfo3D::eT44]);
    logMessage(2, "  Edge swapping stats: \n");
    for (int i = 3; i < 11; i++) {
      logMessage(2, "    %2d tets for %2d:  %6d out of %6d attempts\n", i,
		 2 * i - 4, m_edgeSuccesses[i], m_edgeAttempts[i]);
    }
  }

  template<class myType>
    void
    SwapManager3D::LoadBalancing(std::vector<myType> VecToBalanced[],
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

  void
  SwapManager3D::faceQueueManager3D()
  {
    double FQM_Time = omp_get_wtime();
    std::vector<FaceEdgeInfo> FSBuffer[NUM_PROCS];
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      FaceStructs[ID].clear();
      RemainingFaces[ID].clear();
      FSBuffer[ID].clear();
    }

    omp_set_dynamic(0);
    omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(FSBuffer)
    {
      double SwappingCheck = omp_get_wtime();
      bool IsInUse = false;
      int ID = omp_get_thread_num();
      for (GR_index_t fN = 0; fN < faceQueueVecs[ID].size(); ++fN) {
	Face *pF = faceQueueVecs[ID][fN];
	VALIDATE_INPUT(pF->isValid());
	FaceSwapInfo3D FC;
	FC.categorizeFace(pF);
	FaceSwapInfo3D::faceCat fc = FC.getFaceCat();
	switch (fc)
	  {
	  case FaceSwapInfo3D::eT22:
	  case FaceSwapInfo3D::eT23:
	  case FaceSwapInfo3D::eT32:
	  case FaceSwapInfo3D::eT44:
	    if (!m_pSD->doFaceSwap(FC))
	      continue;
	    else {
	      FaceEdgeInfo FEI1;
	      FEI1.FCS = FC;
	      FEI1.FaceOrEdge = true;
	      assert(FEI1.FCS.getFaceCat() != FaceSwapInfo3D::eOther);
	      FSBuffer[ID].push_back(FEI1);
	    }
	    break;
	  case FaceSwapInfo3D::eN32:
	  case FaceSwapInfo3D::eN44:
	    {
	      EdgeSwapInfo ES(FC);
	      if (!m_pSD->doEdgeSwap(ES))
		continue;
	      else {
		FaceEdgeInfo FEI2;
		FEI2.FCS = FC;
		FEI2.ESS = ES;
		FEI2.FaceOrEdge = false;
		assert(FEI2.FCS.getFaceCat() != FaceSwapInfo3D::eOther);
		FSBuffer[ID].push_back(FEI2);
	      }
	    }
	    break;
	  default:
	    break;
	  }
      }
      faceQueueVecs[ID].clear();
      for (GR_index_t fN = 0; fN < FSBuffer[ID].size(); ++fN) {
	IsInUse = false;
	FaceEdgeInfo FEI3 = FSBuffer[ID][fN];
	FaceSwapInfo3D FC = FEI3.FCS;
	FaceSwapInfo3D::faceCat fc = FC.getFaceCat();

	switch (fc)
	  {
	  case FaceSwapInfo3D::eT22:
	  case FaceSwapInfo3D::eT23:
	  case FaceSwapInfo3D::eT32:
	  case FaceSwapInfo3D::eT44:
	    {
	      m_faceAttemptsThreads[fc][ID]++;
	      int ENDV = 0;
	      if (fc == FaceSwapInfo3D::eT44)
		ENDV = 6;
	      else
		ENDV = 5;
	      double CT = omp_get_wtime();
#pragma omp critical(setInUseVert111)
	      {
		Vert* pV = FC.getVertA();
		pV->isInUseWOC(IsInUse);
		for (int i = 1; i < ENDV && !IsInUse; i++) {
		  pV = FC.getVert(i);
		  pV->isInUseWOC(IsInUse);
		}
		if (!IsInUse) {
		  pV = FC.getVertA();
		  pV->setInUseWOC(true);
		  for (int i = 1; i < ENDV; i++) {
		    pV = FC.getVert(i);
		    pV->setInUseWOC(true);
		  }
		}
	      }
	      CT = omp_get_wtime() - CT;
	      CriticalTime[ID][0] += CT;
	      if (IsInUse) {
		Face* pFinStudy = FC.getFace();
		RemainingFaces[ID].push_back(pFinStudy);
		break;
	      }
	      else {
		FaceStructs[ID].push_back(FEI3);
	      }
	    }
	    break;
	  case FaceSwapInfo3D::eN32:
	  case FaceSwapInfo3D::eN44:
	    {
	      EdgeSwapInfo ES(FC);
	      m_edgeAttemptsThreads[ES.getNumOrigCells()][ID]++;
	      double CT = omp_get_wtime();
#pragma omp critical(setInUseVert111)
	      {
		for (std::set<Vert*>::iterator it = ES.m_TotalVerts.begin();
		    it != ES.m_TotalVerts.end() && !IsInUse; ++it) {
		  Vert* pV = *it;
		  pV->isInUseWOC(IsInUse);
		}
		if (!IsInUse) {
		  for (std::set<Vert*>::iterator it = ES.m_TotalVerts.begin();
		      it != ES.m_TotalVerts.end(); ++it) {
		    Vert* pVl = *it;
		    pVl->setInUseWOC(true);
		  }
		}
	      }
	      CT = omp_get_wtime() - CT;
	      CriticalTime[ID][0] += CT;
	      if (IsInUse) {
		Face* pFif1 = FC.getFace();
		RemainingFaces[ID].push_back(pFif1);
		break;
	      }
	      FaceStructs[ID].push_back(FEI3);
	    }
	    break;
	  default:
	    break;
	  } // End of Switch

      } // End of for
      SwappingCheck = omp_get_wtime() - SwappingCheck;
      FQMThread[ID][0] += SwappingCheck;
    } // End of parallel

    FQM_Time = omp_get_wtime() - FQM_Time;
    FaceQueueManagerTime += FQM_Time;
  }

  int
  SwapManager3D::swapFace(const FaceEdgeInfo& FEI, int ID, bool OMP_Flag)
  {
    FaceSwapInfo3D::faceCat fc = FEI.FCS.getFaceCat();
    if (FEI.FaceOrEdge) // ******* swapping faces *****
    {
      int ENDV = 0;
      if (fc == FaceSwapInfo3D::eT44)
	ENDV = 6;
      else
	ENDV = 5;
      bool reconfigured = reconfigure(FEI.FCS, ID, OMP_Flag);
      if (reconfigured)
	m_faceSuccessesThreads[fc][ID]++;
      Vert* pV = FEI.FCS.getVertA();
      pV->setInUse(false);
      for (int i = 1; i < ENDV; i++) {
	pV = FEI.FCS.getVert(i);
	pV->setInUse(false);
      }
      return 1 + (fc == FaceSwapInfo3D::eT44 ? 1 : 0);
    }
    else //****** swapping edges ********
    {
      reconfigureEdge(FEI.ESS, ID, OMP_Flag);
      m_edgeSuccessesThreads[FEI.ESS.getNumOrigCells()][ID]++;
      for (std::set<Vert*>::iterator it = FEI.ESS.m_TotalVerts.begin();
	  it != FEI.ESS.m_TotalVerts.end(); ++it) {
	Vert* pV = *it;
	pV->setInUse(false);
      }
      return (FEI.ESS.getNumOrigVerts() - 2);
    }
  }

  int
  SwapManager3D::swapFace(Face* pF, int ID, bool OMP_Flag)
  {
    if (!OMP_Flag) {
      VALIDATE_INPUT(pF->isValid());
      if (pF->isDeleted())
	return 0;
      FaceSwapInfo3D FC;
      FC.categorizeFace(pF);
      return swapFace(FC);
    }
    else { // This Section is only for SwapInsertor
      assert(SwapInsertor);
      VALIDATE_INPUT(pF->isValid());
      if (pF->isDeleted())
	return 0;
      FaceSwapInfo3D FC;
      FC.categorizeFace(pF);
      return swapFace(FC, ID, OMP_Flag);
    }
  }

  int
  SwapManager3D::swapFace(const FaceSwapInfo3D& FC, int ID, bool OMP_Flag)
  {
    if (!OMP_Flag) {
      FaceSwapInfo3D::faceCat fc = FC.getFaceCat();
      switch (fc)
	{
	case FaceSwapInfo3D::eT22:
	case FaceSwapInfo3D::eT23:
	case FaceSwapInfo3D::eT32:
	case FaceSwapInfo3D::eT44:
	  m_faceAttempts[fc]++;
	  if (m_pSD->doFaceSwap(FC) && reconfigure(FC)) {
	    m_faceSuccesses[fc]++;
	    return 1 + (fc == FaceSwapInfo3D::eT44 ? 1 : 0);
	  }
	  else
	    return 0;
	  break;
	case FaceSwapInfo3D::eN32:
	case FaceSwapInfo3D::eN44:
	  {
	    EdgeSwapInfo ES(FC);

	    m_edgeAttempts[ES.getNumOrigCells()]++;
	    if (m_pSD->doEdgeSwap(ES)) {
	      reconfigureEdge(ES);
	      m_edgeSuccesses[ES.getNumOrigCells()]++;
	      return (ES.getNumOrigVerts() - 2);
	    }
	    else
	      return 0;
	  }
	  break;
//        // Attempt to put bdry edge swapping for quality
//        // Can't swap bdry faces, but can swap bdry edges
//        // provided they don't lie on a subseg and that
//        // the angle formed by the normals of adjacent bdry faces is
//        // sufficiently small
//        case FaceSwapInfo3D::eBdry: {
//        	// do I need to check both faces?
//        	EdgeSwapInfo ES(FC);
//        	ES.
//        }
//          break;
	default:
	  return 0;
	  break;
	}
    }
    else	// This Section is only for SwapInsertor
    {
      assert(SwapInsertor);
      FaceSwapInfo3D::faceCat fc = FC.getFaceCat();
      switch (fc)
	{
	case FaceSwapInfo3D::eT22:
	case FaceSwapInfo3D::eT23:
	case FaceSwapInfo3D::eT32:
	case FaceSwapInfo3D::eT44:
	  m_faceAttempts[fc]++;
	  if (m_pSD->doFaceSwap(FC) && reconfigure(FC, ID, OMP_Flag)) {
	    m_faceSuccesses[fc]++;
	    return 1 + (fc == FaceSwapInfo3D::eT44 ? 1 : 0);
	  }
	  else
	    return 0;
	  break;
	case FaceSwapInfo3D::eN32:
	case FaceSwapInfo3D::eN44:
	  {
	    EdgeSwapInfo ES(FC);

	    m_edgeAttempts[ES.getNumOrigCells()]++;
	    if (m_pSD->doEdgeSwap(ES)) {
	      reconfigureEdge(ES, ID, OMP_Flag);
	      m_edgeSuccesses[ES.getNumOrigCells()]++;
	      return (ES.getNumOrigVerts() - 2);
	    }
	    else
	      return 0;
	  }
	  break;
//        // Attempt to put bdry edge swapping for quality
//        // Can't swap bdry faces, but can swap bdry edges
//        // provided they don't lie on a subseg and that
//        // the angle formed by the normals of adjacent bdry faces is
//        // sufficiently small
//        case FaceSwapInfo3D::eBdry: {
//        	// do I need to check both faces?
//        	EdgeSwapInfo ES(FC);
//        	ES.
//        }
//          break;
	default:
	  return 0;
	  break;
	}
    }
    // we should not end up here anymore. This function is obsolete

    // Should never reach here.
    return 0;
  }
// takes in a boundary face with edge to swap and two verts connected to edge
  bool
  SwapManager3D::swapBdryEdge(BFace * bface, Vert * v0, Vert * v1,
			      Vert * vOther)
  {

    // first make sure if one is a curves, that they do NOT lie on a subseg
    if ((v0->getVertType() == Vert::eBdryCurve
	|| v1->getVertType() == Vert::eBdryCurve)
	&& m_pVM->areSubsegVertsConnected(v0, v1) != NULL)
      return false;

    Face * pF = bface->getFace(0);

    FaceSwapInfo3D FC(pF, v0, v1, vOther);
    EdgeSwapInfo ES(FC);
    if (!ES.isSwappableEdge())
      return false;
    m_edgeAttempts[ES.getNumOrigCells()]++;
    assert(ES.isBdryEdge());
    if (bface->getType() == Cell::eTriBFace) {
      if (m_pSD->doBdryEdgeSwap(ES)) {
	reconfigureEdge(ES);
	m_edgeSuccesses[ES.getNumOrigCells()]++;
	return true;
      }
      return false;
    }
    else {

      Face * pF2 = bface->getFace(1);
      FaceSwapInfo3D FC2(pF2, v0, v1, vOther);
      EdgeSwapInfo ES2(FC2);
      if (!ES2.isSwappableEdge())
	return false;
      // going to do a check to see they agree on boundary faces;
//    	BFace * bface11 = ES.getBFace(0);
//    	BFace * bface12 = ES.getBFace(1);
//    	BFace * bface21 = ES2.getBFace(0);
//    	BFace * bface22 = ES2.getBFace(1);
//    	logMessage(1,"bface ptrs %p %p | %p %p\n",bface11,bface12,bface21,bface22);
//    	 return false;
      //  check if both sides can do the swap, if they both can, do the swap, with a forced removal
      bool doSwap1 = m_pSD->doBdryEdgeSwap(ES);
      bool doSwap2 = m_pSD->doBdryEdgeSwap(ES2);

      if (doSwap1 && doSwap2) {
	reconfigureIntBdryEdge(ES, ES2);
	return true;
      }
    }
    return false;
  }
  int
  SwapManager3D::swapAllBdryEdges(bool OMP_Flag)
  {
    int iSwaps = 0;
    // collect all edges.
    // need a way to pick edges out but edges not directly stored
    if (!OMP_Flag) {
      for (GR_index_t iBF = 0; iBF < m_pVM->getNumTotalBdryFaces(); iBF++) {
	BFace * bface = m_pVM->getBFace(iBF);
	if (bface->isDeleted())
	  continue;
	// lets swap each edge of these faces
	Face * pF = bface->getFace(0);
	bool qSwapped = false;
	// find the first nonApex vert to try. Preference of edges is
	// Bdry-Bdry, Bdry-BdryCurve, Bdry-BdryApex, BdryCurve-BdryCurve (if not subseg)
	Vert *verts[] =
	  { pVInvalidVert, pVInvalidVert, pVInvalidVert };
	for (int iV = 0; iV < 3; iV++) {
	  if (pF->getVert(iV)->getVertType() != Vert::eBdryApex) {
	    verts[0] = pF->getVert(iV);
	    verts[1] = pF->getVert((iV + 1) % 3);
	    verts[2] = pF->getVert((iV + 2) % 3);
	    break;
	  }
	}
	// Haven't found a vert that I can attempt to swap its edge from so skip this face and move on
	if (!verts[0])
	  continue;
	// pick off pairs of verts. Try 0 and 1
	{
	  qSwapped = swapBdryEdge(bface, verts[0], verts[1], verts[2]);
	}
	if (!qSwapped) {
	  qSwapped = swapBdryEdge(bface, verts[0], verts[2], verts[1]);
	}
	// Finally, try verts[1] and verts[2] if they are swappable
	if (!qSwapped
	    && (verts[1]->getVertType() == Vert::eBdry
		|| verts[1]->getVertType() == Vert::eBdryTwoSide
		|| verts[2]->getVertType() == Vert::eBdry
		|| verts[2]->getVertType() == Vert::eBdryTwoSide)) {
	  qSwapped = swapBdryEdge(bface, verts[1], verts[2], verts[0]);

	}
	iSwaps += qSwapped;
      }
    }
    return iSwaps;
  }
  int
  SwapManager3D::swapAllFaces(bool OMP_Flag)
  {
    double SwappingTime = omp_get_wtime();
    if (!OMP_Flag) {
      faceQueue.clear();
      // Initialize the queue with all the faces in the mesh.
      EntContainer<TriFace>::iterator ECiter = m_pVM->triFace_begin(),
	  ECiterEnd = m_pVM->triFace_end();
      for (; ECiter != ECiterEnd; ++ECiter) {
	Face *pF = &(*ECiter);
	faceQueue.insert(pF);
      }
      int NumTotSwaps = swapAllQueuedFacesRecursive();
      SwappingTime = omp_get_wtime() - SwappingTime;
      TotalSwappingTime += SwappingTime;
      return NumTotSwaps;
    }
    else {
      double initial_query_time = omp_get_wtime();
      double paralleltime = omp_get_wtime();
      GR_index_t chunk_size = GR_index_t(
	  double(m_pVM->getNumFaces()) / NUM_PROCS);
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(chunk_size)
      {
	EntContainer<TriFace>::iterator ECiter = m_pVM->triFace_begin();
	EntContainer<TriFace>::iterator ECiterEnd = m_pVM->triFace_begin();

	int ID = omp_get_thread_num();
	if (ID < NUM_PROCS - 1) {
	  ECiter.jump(chunk_size * ID);
	  ECiterEnd.jump(chunk_size * (ID + 1));

	}
	else {
	  EntContainer<TriFace>::iterator ECiterEndBuff = m_pVM->triFace_end();

	  ECiter.jump(chunk_size * ID);
	  ECiterEnd = ECiterEndBuff;

	}
	for (; ECiter != ECiterEnd; ++ECiter) {
	  Face *pF = &(*ECiter);
	  faceQueueVecs[ID].push_back(pF);
	}
      }

      initial_query_time = omp_get_wtime() - initial_query_time;
      Quary_Time += initial_query_time;

      paralleltime = omp_get_wtime() - paralleltime;
      Parallel_Time += paralleltime;
      paralleltime = 0;

      Sync_Entry_Time += m_pVM->Mesh3D_sync_entry(paralleltime, true);

      Parallel_Time += paralleltime;

      int NumSwaps = swapAllQueuedFacesRecursive(OMP_Flag, false);

      Sync_Exit_Time += m_pVM->Mesh3D_sync_Exit(true);

      for (int i = 0; i < NUM_PROCS; i++)
	faceQueueVecs[i].clear();
      SwappingTime = omp_get_wtime() - SwappingTime;
      TotalSwappingTime += SwappingTime;
      return NumSwaps;
    }
    //return 0;
  }

  int
  SwapManager3D::swapAllQueuedFacesRecursive(bool OMP_Flag, bool SwapIn)
  {
    if (!OMP_Flag) {
      parallelswapping = OMP_Flag;
      // m_pVM->sendEvents();
      int numSwaps = 0, swapsThisPass = 0, passes = 0;
      do {
	int attempts = faceQueue.size();
	swapsThisPass = swapAllQueuedFaces(OMP_Flag);
	numSwaps += swapsThisPass;
	passes++;

	logMessage(2, "Pass %d, %d swaps out of %d attempts\n", passes,
		   swapsThisPass, attempts);
	m_pVM->sendEvents();
	logMessage(2, "  Master swap queue now has %zu faces.\n",
		   faceQueue.size());
#ifdef SIM_ANNEAL_TEST
	logMessage(2, "  Sim anneal temp is %f\n", simAnnealTemp);
	simAnnealTemp /= 1.1;
      }while (swapsThisPass != 0 && simAnnealTemp > 1.e-6);
#else
      }
      while (swapsThisPass != 0 && passes < 20);
#endif

      return numSwaps;
    }
    else if (!SwapIn) {
      assert(OMP_Flag);
      m_pVM->IsInParallel(true);
      int RFZ = 0;
      double SendTime;
      parallelswapping = OMP_Flag;
      assert(OMP_Flag);
      int numSwaps = 0, swapsThisPass = 0, passes = 0;
      do {
	RFZ = 0;
	LoadBalancing(faceQueueVecs, NUM_PROCS);
	faceQueueManager3D();
	LoadBalancing(FaceStructs, NUM_PROCS);
	swapsThisPass = swapAllQueuedFaces(OMP_Flag);
	numSwaps += swapsThisPass;

	if (SwapIn) {
	  faceQueue.clear();
	  for (int ID = 0; ID < NUM_PROCS; ID++) {
	    for (std::vector<Face*>::iterator it = RemainingFaces[ID].begin();
		it != RemainingFaces[ID].end(); ++it)
	      faceQueue.insert(*it);
	    RemainingFaces[ID].clear();
	  }
	}
	else {
	  for (int ID = 0; ID < NUM_PROCS; ID++) {
	    faceQueueVecs[ID].clear();
	    //faceQueueVecs[ID].reserve(RemainingFaces[ID].size());
	  }
	  omp_set_dynamic(0);
	  omp_set_num_threads(NUM_PROCS);
#pragma omp parallel
	  {
	    int ID = omp_get_thread_num();
	    for (std::vector<Face*>::iterator it = RemainingFaces[ID].begin();
		it != RemainingFaces[ID].end(); ++it) {
	      if (!(*it)->isDeleted()) {
		faceQueueVecs[ID].push_back(*it);
	      }
	    }
	    RemainingFaces[ID].clear();
	  }
	}

	SendTime = omp_get_wtime();
	m_pVM->sendEvents();
	SendTime = omp_get_wtime() - SendTime;
	SendEventsTime += SendTime;

	for (GR_index_t ID = 0; ID < NUM_PROCS; ID++)
	  RFZ += faceQueueVecs[ID].size();
	logMessage(2, "   %d swaps at pass %d.\n", swapsThisPass, passes);
	logMessage(2, "   Number of remaining faces is: %d.\n", RFZ);

	for (GR_index_t ID = 0; ID < NUM_PROCS; ID++)
	  logMessage(2, "   faceQueueVecs[%u] now has %lu faces.\n", ID,
		     faceQueueVecs[ID].size());

	passes++;

#ifdef SIM_ANNEAL_TEST
	logMessage(2, "  Sim anneal temp is %f\n", simAnnealTemp);
	simAnnealTemp /= 1.1;
      }while ((swapsThisPass != 0 || RFZ >0) && simAnnealTemp > 1.e-6);
#else
      }
      while (swapsThisPass != 0 || RFZ > 0);
#endif

      m_pVM->IsInParallel(false);
      return numSwaps;
    }
    else {   // ***** This section has been added for SwapInsertor
      assert(SwapIn);
      parallelswapping = OMP_Flag;
      SwapInsertor = SwapIn;
      // m_pVM->sendEvents();
      int numSwaps = 0, swapsThisPass = 0, passes = 0;
      do {
	int attempts = faceQueue.size();
	swapsThisPass = swapAllQueuedFaces(OMP_Flag);
	numSwaps += swapsThisPass;
	passes++;

	logMessage(2, "Pass %d, %d swaps out of %d attempts\n", passes,
		   swapsThisPass, attempts);
	m_pVM->sendEvents();
	logMessage(2, "  Master swap queue now has %zu faces.\n",
		   faceQueue.size());
#ifdef SIM_ANNEAL_TEST
	logMessage(2, "  Sim anneal temp is %f\n", simAnnealTemp);
	simAnnealTemp /= 1.1;
      }while (swapsThisPass != 0 && simAnnealTemp > 1.e-6);
#else
      }
      while (swapsThisPass != 0 && passes < 20);
#endif

      return numSwaps;
    }
  }
  int
  SwapManager3D::swapAllQueuedFaces(bool OMP_Flag)
  {
    parallelswapping = OMP_Flag;
    if (!OMP_Flag) {
      m_pVM->sendEvents();
      std::set<Face*> liveFaceQueue;

      // Initialize the queue with all the faces in the queue.
      swap(faceQueue, liveFaceQueue);

      std::set<Face*>::iterator iter = liveFaceQueue.begin(), iterEnd =
	  liveFaceQueue.end();
      int numSwaps = 0;
      for (; iter != iterEnd; ++iter) {
	Face *pF = *iter;
	numSwaps += swapFace(pF);
      } // for

      return numSwaps;
    }
    else if (!SwapInsertor) {
      assert(OMP_Flag);
      //LoadBalancing <FaceEdgeInfo> (FaceStructs,NUM_PROCS);
      int numSwapps[NUM_PROCS][100] =
	{
	  { 0 } };
      int numTotalSwapps = 0;
      double swappingtime1 = omp_get_wtime();
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(numSwapps)
      {
	double ThreadReconTime = omp_get_wtime();
	int ID = omp_get_thread_num();
	for (GR_index_t fN = 0; fN < FaceStructs[ID].size(); ++fN)
	  numSwapps[ID][0] += swapFace(FaceStructs[ID][fN], ID, true);
	ThreadReconTime = omp_get_wtime() - ThreadReconTime;
	ReconTime[ID][0] += ThreadReconTime;
      }
      swappingtime1 = omp_get_wtime() - swappingtime1;
      Parallel_Time += swappingtime1;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	numTotalSwapps += numSwapps[ID][0];
	NumSwapps[ID] += numSwapps[ID][0];
      }
      return numTotalSwapps;
    }
    else {
      assert(SwapInsertor);
      int ID = omp_get_thread_num();
      m_pVM->sendEvents();
      std::set<Face*> liveFaceQueue;

      // Initialize the queue with all the faces in the queue.
      swap(faceQueue, liveFaceQueue);

      std::set<Face*>::iterator iter = liveFaceQueue.begin(), iterEnd =
	  liveFaceQueue.end();
      int numSwaps = 0;
      for (; iter != iterEnd; ++iter) {
	Face *pF = *iter;
	numSwaps += swapFace(pF, ID, OMP_Flag);
      } // for

      return numSwaps;
    }
    assert(0); // you should not get here!
    return 0;
  }

  bool
  SwapManager3D::reconfigure(const FaceSwapInfo3D& FC, int ID, bool OMP_Flag)
  {
    assert(!FC.getFace()->isLocked());

    Cell* pTets[4];
    int nTets;
    FC.getTets(nTets, pTets);

    Vert *pVA = FC.getVertA();
    Vert *pVB = FC.getVertB();
    Vert *pVC = FC.getVertC();
    Vert *pVD = FC.getVertD();
    Vert *pVE = FC.getVertE();
    Vert *pVF = FC.getVertF();

    assert(pVA->isValid());
    assert(pVB->isValid());
    assert(pVC->isValid());

    switch (FC.getFaceCat())
      {
      case FaceSwapInfo3D::eT22:
	return (reconfTet22(pTets, pVA, pVB, pVC, pVD, pVE, pVF, FC.getBFaceA(),
			    FC.getBFaceB(), ID, OMP_Flag));
      case FaceSwapInfo3D::eT23:
	return (reconfTet23(pTets, pVA, pVB, pVC, pVD, pVE, pVF, ID, OMP_Flag));
      case FaceSwapInfo3D::eT32:
	return (reconfTet32(pTets, pVA, pVB, pVC, pVD, pVE, pVF, ID, OMP_Flag));
      case FaceSwapInfo3D::eT44:
	return (reconfTet44(pTets, pVA, pVB, pVC, pVD, pVE, pVF, ID, OMP_Flag));
      default:
	assert2(0, "Bad face type in reconfigure");
	return false;
      }
  }

//@@ Swap two tets for three.
  bool
  SwapManager3D::reconfTet23(Cell* pTets[4], Vert * const pVVertA,
			     Vert* const pVVertB, Vert * const pVVertC,
			     Vert * const pVVertD, Vert * const pVVertE,
			     Vert * const pVVertF, int ID, bool OMP_Flag)
  {
    assert(
	pTets[0] != NULL && pTets[1] != NULL && pVVertA != NULL && pVVertB != NULL && pVVertC != NULL && pVVertD != NULL && pVVertE != NULL);
    assert(pTets[2] == NULL && pTets[3] == NULL && pVVertF == NULL);
    // Identify both cells and all seven faces in the two-tet config.
    Cell* pTetD = pTets[0];
    Cell* pTetE = pTets[1];
    int iRegion = pTetD->getRegion();

    // Delete old tets
    if (!OMP_Flag) {
      m_pVM->deleteCell(pTetD);
      m_pVM->deleteCell(pTetE);
    }
    else {
      m_pVM->deleteCellParallel(pTetD, ID);
      m_pVM->deleteCellParallel(pTetE, ID);
    }

    // Create three new tets
    bool qExist;

    TetCell *pTCTetA = m_pVM->createTetCell(qExist, pVVertB, pVVertC, pVVertD,
					    pVVertE, iRegion, OMP_Flag, ID);
    assert(!qExist);

    TetCell *pTCTetB = m_pVM->createTetCell(qExist, pVVertC, pVVertA, pVVertD,
					    pVVertE, iRegion, OMP_Flag, ID);
    assert(!qExist);

    TetCell *pTCTetC = m_pVM->createTetCell(qExist, pVVertA, pVVertB, pVVertD,
					    pVVertE, iRegion, OMP_Flag, ID);
    assert(!qExist);

    assert(!pTCTetA->hasVert(pVVertA));
    assert(!pTCTetB->hasVert(pVVertB));
    assert(!pTCTetC->hasVert(pVVertC));
    // Orientation and correctness tests for the tets are done at
    // creation time; it's a waste of time to repeat the assertions.

    return true;
  }

//@@ Swap three tets for two.
  bool
  SwapManager3D::reconfTet32(Cell* pTets[4], Vert * const pVVertA,
			     Vert* const pVVertB, Vert * const pVVertC,
			     Vert * const pVVertD, Vert * const pVVertE,
			     Vert * const pVVertF, int ID, bool OMP_Flag)
  {
    assert(
	pTets[0] != NULL && pTets[1] != NULL && pTets[2] != NULL && pVVertA != NULL && pVVertB != NULL && pVVertC != NULL && pVVertD != NULL && pVVertE != NULL);
    assert(pTets[3] == NULL && pVVertF == NULL);
    Cell* pTetA = pTets[0];
    Cell* pTetB = pTets[1];
    Cell* pTetC = pTets[2];
    int iRegion = pTetA->getRegion();
    if ((iRegion != pTetB->getRegion()) || (iRegion != pTetC->getRegion()))
      return false;

    // Delete tets
    if (!OMP_Flag) {
      m_pVM->deleteCell(pTetA);
      m_pVM->deleteCell(pTetB);
      m_pVM->deleteCell(pTetC);
    }
    else {
      m_pVM->deleteCellParallel(pTetA, ID);
      m_pVM->deleteCellParallel(pTetB, ID);
      m_pVM->deleteCellParallel(pTetC, ID);
    }

    // Create new tets
    bool qExist;
    TetCell *pTetD = m_pVM->createTetCell(qExist, pVVertA, pVVertC, pVVertB,
					  pVVertD, iRegion, OMP_Flag, ID);
    assert(!qExist);

    TetCell *pTetE = m_pVM->createTetCell(qExist, pVVertA, pVVertB, pVVertC,
					  pVVertE, iRegion, OMP_Flag, ID);
    assert(!qExist);

    assert(pTetD->hasVert(pVVertD));
    assert(pTetE->hasVert(pVVertE));

    return true;
  }

//@@ Swap two tets for two, in the case where two faces are coplanar.
  bool
  SwapManager3D::reconfTet22(Cell* pTets[4], Vert * const pVVertA,
			     Vert* const pVVertB, Vert * const pVVertC,
			     Vert * const pVVertD, Vert * const pVVertE,
			     Vert * const pVVertF, TriBFace * const pTBFA,
			     TriBFace * const pTBFB, int ID, bool OMP_Flag)
  {
    // This routine does only the 2-to-2 case.  4-to-4 is done elsewhere.
    // The cells and verts are assumed to have been set up properly already.

    /*      	        C
     /\
        	      /||\
        	     / |  \                 A,D are on the base;
     /  | | \                B,E are front, C,F are back.
     /   |    \               Cell A is on the left, B on
     /    |  |  \	            the right.
     /     |      \
  		/      |   |   \
  	       /     C |     F  \
  	      /	       |  E |    \
  	     /   B     | _ - B    \
  	    /  _ _ - - |    /  - _ \
  	   D -         |           -E
     -__	    A  |  /  D  __--
     -__     |     __-
     -___ |/__--
     -A-

     */
    assert(
	pTets[0] != NULL && pTets[1] != NULL && pVVertA != NULL && pVVertB != NULL && pVVertC != NULL && pVVertD != NULL && pVVertE != NULL);
    assert(pTets[2] == NULL && pTets[3] == NULL && pVVertF == NULL);
    Cell* pTCTetA = pTets[0];
    Cell* pTCTetB = pTets[1];

    // Do the actual two-for-two swap here.  The faces all stay in the
    // same places except for pFFace, the dividing face, which gets turned.
    // Cell A ends up in front, cell B in back.  The only vert hints
    // affected are for verts 1 and 2, which can be set explicitly.

    int iRegion = pTCTetA->getRegion();
    assert(iRegion == pTCTetB->getRegion());

    // Delete old cells
    if (!OMP_Flag) {
      m_pVM->deleteCell(pTCTetA);
      m_pVM->deleteCell(pTCTetB);
    }
    else {
      m_pVM->deleteCellParallel(pTCTetA, ID);
      m_pVM->deleteCellParallel(pTCTetB, ID);
    }

    bool qExist;
    pTCTetA = m_pVM->createTetCell(qExist, pVVertD, pVVertA, pVVertE, pVVertC,
				   iRegion, OMP_Flag, ID);
    assert(!qExist);

    pTCTetB = m_pVM->createTetCell(qExist, pVVertE, pVVertB, pVVertD, pVVertC,
				   iRegion, OMP_Flag, ID);
    assert(!qExist);

    Face *pFBdryA = pTCTetA->getOppositeFace(pVVertC);
    Face *pFBdryB = pTCTetB->getOppositeFace(pVVertC);
    if (!OMP_Flag) {
      m_pVM->createBFace(pFBdryA, pTBFA);
      m_pVM->createBFace(pFBdryB, pTBFA);

      m_pVM->deleteBFace(pTBFA);
      m_pVM->deleteBFace(pTBFB);

    }
    else {
      m_pVM->createBFaceParallel(pFBdryA, pTBFA, ID);
      m_pVM->createBFaceParallel(pFBdryB, pTBFA, ID);

      m_pVM->deleteBFaceParallel(pTBFA, ID);
      m_pVM->deleteBFaceParallel(pTBFB, ID);

    }
    // Delete bdry faces now

    return true;
  }

//@@ Swap four tets for four, in the case where two faces are coplanar.
  bool
  SwapManager3D::reconfTet44(Cell* pTets[4], Vert * const pVVertA,
			     Vert* const pVVertB, Vert * const pVVertC,
			     Vert * const pVVertD, Vert * const pVVertE,
			     Vert * const pVVertF, int ID, bool OMP_Flag)
  {
    // The cells and verts are assumed to have been set up properly already.

    // Here's the before picture:
    /*      	        C
     /\
        	      /||\
        	     / |  \                 A,D are on the base;
     /  | | \                B,E are front, C,F are back.
     /   |    \               Cell A is on the left, B on
     /    |  |  \	            the right.
     /     |      \
  		/      |   |   \
  	       /     C |     F  \
  	      /	       |  E |    \
  	     /   B     | _ - B    \
  	    /  _ _ - - |    /  - _ \
  	   D -         |           -E
     -__	    A  |  /  D  __--
     -__     |     __-
     -___ |/__--
     -A-

     ___---B_
     ____----	    /  --__
     D---            /       -E
     - _	    A     / |D  _ -
     \  - _       /   _ -   /
     \     -_   / _-      /
     \       -A-  |     /
     \     C1|     F1 /
     \      |  |    /
     \ B1  |  E1  /
     \    |     /	            A,D are on the base;
     \   | |  /               B1,E1 are front, C1,F1 are back.
     \  |   /                Cell A1 is on the left, B1 on
     \ || /                 the right.
     \| /                  Face G divides these two cells.
     \/
     C1

     */
    assert(
	pTets[0] != NULL && pTets[1] != NULL && pTets[2] != NULL && pTets[3] != NULL && pVVertA != NULL && pVVertB != NULL && pVVertC != NULL && pVVertD != NULL && pVVertE != NULL);
    Cell* pTCTetA = pTets[0];
    Cell* pTCTetB = pTets[1];
    Cell* pTCTetC = pTets[2];
    Cell* pTCTetD = pTets[3];

    int iReg = pTCTetA->getRegion();
    int iReg1 = pTCTetC->getRegion();
    // Delete all four tets
    if (!OMP_Flag) {
      m_pVM->deleteCell(pTCTetA);
      m_pVM->deleteCell(pTCTetB);
      m_pVM->deleteCell(pTCTetC);
      m_pVM->deleteCell(pTCTetD);
    }
    else {
      m_pVM->deleteCellParallel(pTCTetA, ID);
      m_pVM->deleteCellParallel(pTCTetB, ID);
      m_pVM->deleteCellParallel(pTCTetC, ID);
      m_pVM->deleteCellParallel(pTCTetD, ID);
    }

    // Now create four new tets directly from verts.
    bool qExist;
    pTCTetA = m_pVM->createTetCell(qExist, pVVertD, pVVertA, pVVertE, pVVertC,
				   iReg, OMP_Flag, ID);
    assert(!qExist);
    pTCTetB = m_pVM->createTetCell(qExist, pVVertE, pVVertB, pVVertD, pVVertC,
				   iReg, OMP_Flag, ID);
    assert(!qExist);

    pTCTetC = m_pVM->createTetCell(qExist, pVVertE, pVVertA, pVVertD, pVVertF,
				   iReg1, OMP_Flag, ID);
    assert(!qExist);
    pTCTetD = m_pVM->createTetCell(qExist, pVVertD, pVVertB, pVVertE, pVVertF,
				   iReg1, OMP_Flag, ID);
    assert(!qExist);
    return true;
  }

  void
  SwapManager3D::receiveDeletedFaces(std::vector<Face*>& deletedFaces)
  {
    logMessage(4, "SwapManager3D received word about %zu deleted faces.\n",
	       deletedFaces.size());
//    if (!parallelswapping)
//    {
    std::vector<Face*>::iterator iter = deletedFaces.begin(), iterEnd =
	deletedFaces.end();
    for (; iter != iterEnd; ++iter) {
      Face *pF = *iter;
      faceQueue.erase(pF);
    }
//    }
  }

  void
  SwapManager3D::receiveCreatedCells(std::vector<Cell*>& createdCells)
  {
    logMessage(4, "SwapManager3D received word about %zu new cells.\n",
	       createdCells.size());
    std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
	createdCells.end();
    if (!parallelswapping || SwapInsertor) {
      for (; iter != iterEnd; ++iter) {
	Cell *pC = *iter;
	if (pC->getType() != CellSkel::eTet)
	  continue;
	for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pC->getFace(i);
	  faceQueue.insert(pF);
	}
      }
    }
    else {
      int counter = 0;
      for (; iter != iterEnd; ++iter) {
	Cell *pC = *iter;
	if (pC->getType() != CellSkel::eTet)
	  continue;
	for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pC->getFace(i);
	  int ID = counter % NUM_PROCS;
	  counter++;
	  pF->setBadFaceWOC(false);
	  faceQueueVecs[ID].push_back(pF);
	}
      }
    }
  }

  void
  SwapManager3D::receiveMovedVerts(std::vector<Vert*>& movedVerts)
  {
    logMessage(4, "SwapManager3D received word about %zu moved verts.\n",
	       movedVerts.size());
    std::vector<Vert*>::iterator iter = movedVerts.begin(), iterEnd =
	movedVerts.end();
    if (!parallelswapping || SwapInsertor) {
      for (; iter != iterEnd; ++iter) {
	Vert *pV = *iter;
	// Need to be sure to get not just incident faces but also
	// opposite faces.
	for (int i = pV->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pV->getFace(i);
	  Cell *pC = pF->getLeftCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    faceQueue.insert(pC->getFace(ii));
	  }

	  pC = pF->getRightCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    faceQueue.insert(pC->getFace(ii));
	  }
	}
      }
    }
    else {
      int counter = 0;
      for (; iter != iterEnd; ++iter) {
	Vert *pV = *iter;
	for (int i = pV->getNumFaces() - 1; i >= 0; i--) {
	  Face *pF = pV->getFace(i);
	  Cell *pC = pF->getLeftCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    int ID = counter % NUM_PROCS;
	    counter++;
	    pC->getFace(ii)->setBadFaceWOC(false);
	    faceQueueVecs[ID].push_back(pC->getFace(ii));
	  }

	  pC = pF->getRightCell();
	  for (int ii = pC->getNumFaces() - 1; ii >= 0; ii--) {
	    int ID = counter % NUM_PROCS;
	    counter++;
	    pC->getFace(ii)->setBadFaceWOC(false);
	    faceQueueVecs[ID].push_back(pC->getFace(ii));
	  }
	}
      }
    }
  }

  void
  SwapManager3D::reconfigureEdge(const EdgeSwapInfo& ES, int ID, bool OMP_Flag)
  {
    // The interior version; the bdry version calls this, then does a
    // little post-processing to get bfaces set up properly.
    // If the highest-quality alternative is better than the current
    // configuration:
    unsigned numOrigCells = ES.getNumOrigCells();
    unsigned numOrigVerts = ES.getNumOrigVerts();

    // Get the region identifier for this collection of tets.
    int iReg = ES.getOrigCell(0)->getRegion();

    // Delete all of the interior tets and faces.  These must be
    // separate loops, or forced.
    if (!OMP_Flag) {
      assert(omp_in_parallel() == 0);
      for (unsigned i = 0; i < numOrigCells; i++) {
	m_pVM->deleteCell(ES.getOrigCell(i));
	assert(iReg == ES.getOrigCell(i)->getRegion());
      }
    }
    else {
      if (NUM_PROCS > 1)
	assert(omp_in_parallel() != 0);
      for (unsigned i = 0; i < numOrigCells; i++) {
	m_pVM->deleteCellParallel(ES.getOrigCell(i), ID);
	assert(iReg == ES.getOrigCell(i)->getRegion());
      }
    }

    // The old version stored all the newly-created cells for two
    // reasons:
    //   1. So that their faces could be swapped recursively.  With
    //      Observer, that's now handled, including automatically
    //      handling cases where recursion isn't desired.
    //   2. To post-check cells for correctness.  This is done in the
    //      tet cell creation code, so it doesn't need to be done again.

    // Also, old face data isn't needed anymore (it was only used for
    // post-validation, which is handled internal to the new creation
    // calls.

    // Get the new configuration data.
    unsigned bestCanon, bestPerm;
    ES.getNewConfig(bestCanon, bestPerm);
    EdgeConfig & Conf =
	EdgeConfigs::getInstance().getConfSet(numOrigVerts).getConfig(
	    bestCanon);

    int iOrient = 0;

    // There are numOrigVerts - 2 triangles on the equator, two tets apiece,
    // four verts per tet.
    Vert *northVert = ES.getNorthVert();
    Vert *southVert = ES.getSouthVert();

    for (unsigned i = 0; i < numOrigVerts - 2; i++) {
      // These are the equatorial tris, so all these indices should be
      // positive.
      assert(Conf.getFaceVert(i, 0) >= 0);
      assert(Conf.getFaceVert(i, 1) >= 0);
      assert(Conf.getFaceVert(i, 2) >= 0);
      Vert* pVert0 = ES.getVert(
	  (Conf.getFaceVert(i, 0) + bestPerm) % numOrigVerts);
      Vert* pVert1 = ES.getVert(
	  (Conf.getFaceVert(i, 1) + bestPerm) % numOrigVerts);
      Vert* pVert2 = ES.getVert(
	  (Conf.getFaceVert(i, 2) + bestPerm) % numOrigVerts);

      // Calculate the orientation exactly once.
      if (iOrient == 0)
	iOrient = checkOrient3D(pVert0, pVert1, pVert2, northVert);

      if (iOrient == 1) {
	assert(checkOrient3D(pVert0, pVert1, pVert2, northVert) == 1);
	assert(checkOrient3D(pVert0, pVert2, pVert1, southVert) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, northVert,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, southVert,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
      else {
	assert(checkOrient3D(pVert0, pVert2, pVert1, northVert) == 1);
	assert(checkOrient3D(pVert0, pVert1, pVert2, southVert) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, northVert,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, southVert,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
    }

    // Old code marked vertices for smoothing; Observer makes this
    // unneccessary.

    // Tet and face validation is done when they're created.

    // Swap recursion has been replaced by batch queueing of faces that
    // would otherwise be swapped.
    if (ES.isBdryEdge()) {

      // Replace the two bdry faces.
      BFace* oldBFace0 = ES.getBFace(0);
      BFace* oldBFace1 = ES.getBFace(1);

      // Retrieve the faces that were just created on the boundary.
      Face *faceOnBdry0, *faceOnBdry1;
      {
	const Vert *pV0 = ES.getVert(0);
	const Vert *pV1 = ES.getVert(numOrigVerts - 1);

	faceOnBdry0 = findCommonFace(pV0, pV1, northVert);
	assert(faceOnBdry0->isValid());
	// Just created them; they'd -better- exist!
	// Exactly one cell should exist
	assert(
	    faceOnBdry0->getLeftCell()->isValid()
		^ faceOnBdry0->getRightCell()->isValid());

	faceOnBdry1 = findCommonFace(pV1, pV0, southVert);
	assert(faceOnBdry1->isValid());
	// Just created them; they'd -better- exist!
	// Exactly one cell should exist
	assert(
	    faceOnBdry1->getLeftCell()->isValid()
		^ faceOnBdry1->getRightCell()->isValid());
      }

      // Attach the new faces to the bdry faces.
      BFace* newBFace0;
      BFace* newBFace1;
      if (!OMP_Flag) {
	newBFace0 = m_pVM->createBFace(faceOnBdry0, oldBFace0);
	newBFace1 = m_pVM->createBFace(faceOnBdry1, oldBFace1);

	m_pVM->deleteBFace(oldBFace0);
	m_pVM->deleteBFace(oldBFace1);
      }
      else {
	newBFace0 = m_pVM->createBFaceParallel(faceOnBdry0, oldBFace0, ID);
	newBFace1 = m_pVM->createBFaceParallel(faceOnBdry1, oldBFace1, ID);

	m_pVM->deleteBFaceParallel(oldBFace0, ID);
	m_pVM->deleteBFaceParallel(oldBFace1, ID);
      }
      assert(newBFace0->isValid());
      assert(newBFace1->isValid());

      // Delete the old faces on the boundary.
    }
  }
  void
  SwapManager3D::reconfigureIntBdryEdge(const EdgeSwapInfo& ES1,
					const EdgeSwapInfo& ES2, int ID,
					bool OMP_Flag)
  {
    // The interior version; the bdry version calls this, then does a
    // little post-processing to get bfacs set up properly.
    // If the highest-quality alternative is better than the current
    // configuration:
    unsigned numOrigCells1 = ES1.getNumOrigCells();
    unsigned numOrigVerts1 = ES1.getNumOrigVerts();
    assert(ES1.isBdryEdge() && ES2.isBdryEdge()); // should be same edge
	// Get the region identifier for this collection of tets.
    int iReg = ES1.getOrigCell(0)->getRegion();

    // Delete all of the interior tets and facES1.  ThES1e must be
    // separate loops, or forced.
    if (!OMP_Flag) {
      assert(omp_in_parallel() == 0);
      for (unsigned i = 0; i < numOrigCells1; i++) {
	m_pVM->deleteCell(ES1.getOrigCell(i));
	assert(iReg == ES1.getOrigCell(i)->getRegion());
      }
    }
    else {
      if (NUM_PROCS > 1)
	assert(omp_in_parallel() != 0);
      for (unsigned i = 0; i < numOrigCells1; i++) {
	m_pVM->deleteCellParallel(ES1.getOrigCell(i), ID);
	assert(iReg == ES1.getOrigCell(i)->getRegion());
      }
    }

    // Get the new configuration data.
    unsigned bestCanon1, bestPerm1;
    ES1.getNewConfig(bestCanon1, bestPerm1);
    EdgeConfig & Conf1 =
	EdgeConfigs::getInstance().getConfSet(numOrigVerts1).getConfig(
	    bestCanon1);

    int iOrient = 0;

    // There are numOrigVerts - 2 triangles on the equator, two tets apiece,
    // four verts per tet.
    Vert *northVert1 = ES1.getNorthVert();
    Vert *southVert1 = ES1.getSouthVert();

    for (unsigned i = 0; i < numOrigVerts1 - 2; i++) {
      // These are the equatorial tris, so all these indices should be
      // positive.
      assert(Conf1.getFaceVert(i, 0) >= 0);
      assert(Conf1.getFaceVert(i, 1) >= 0);
      assert(Conf1.getFaceVert(i, 2) >= 0);
      Vert* pVert0 = ES1.getVert(
	  (Conf1.getFaceVert(i, 0) + bestPerm1) % numOrigVerts1);
      Vert* pVert1 = ES1.getVert(
	  (Conf1.getFaceVert(i, 1) + bestPerm1) % numOrigVerts1);
      Vert* pVert2 = ES1.getVert(
	  (Conf1.getFaceVert(i, 2) + bestPerm1) % numOrigVerts1);

      // Calculate the orientation exactly once.
      if (iOrient == 0)
	iOrient = checkOrient3D(pVert0, pVert1, pVert2, northVert1);

      if (iOrient == 1) {
	assert(checkOrient3D(pVert0, pVert1, pVert2, northVert1) == 1);
	assert(checkOrient3D(pVert0, pVert2, pVert1, southVert1) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, northVert1,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, southVert1,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
      else {
	assert(checkOrient3D(pVert0, pVert2, pVert1, northVert1) == 1);
	assert(checkOrient3D(pVert0, pVert1, pVert2, southVert1) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, northVert1,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, southVert1,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
    }

    unsigned numOrigCells2 = ES2.getNumOrigCells();
    unsigned numOrigVerts2 = ES2.getNumOrigVerts();

    // Get the region identifier for this collection of tets.
    iReg = ES2.getOrigCell(0)->getRegion();

    // Delete all of the interior tets and faces.  These must be
    // separate loops, or forced.
    if (!OMP_Flag) {
      assert(omp_in_parallel() == 0);
      for (unsigned i = 0; i < numOrigCells2; i++) {
	m_pVM->deleteCell(ES2.getOrigCell(i));
	assert(iReg == ES2.getOrigCell(i)->getRegion());
      }
    }
    else {
      if (NUM_PROCS > 1)
	assert(omp_in_parallel() != 0);
      for (unsigned i = 0; i < numOrigCells2; i++) {
	m_pVM->deleteCellParallel(ES2.getOrigCell(i), ID);
	assert(iReg == ES2.getOrigCell(i)->getRegion());
      }
    }

    // Get the new configuration data.
    unsigned bestCanon2, bestPerm2;
    ES2.getNewConfig(bestCanon2, bestPerm2);
    EdgeConfig & Conf2 =
	EdgeConfigs::getInstance().getConfSet(numOrigVerts2).getConfig(
	    bestCanon2);

    // There are numOrigVerts - 2 triangles on the equator, two tets apiece,
    // four verts per tet.
    Vert *northVert2 = ES2.getNorthVert();
    Vert *southVert2 = ES2.getSouthVert();

    for (unsigned i = 0; i < numOrigVerts2 - 2; i++) {
      // These are the equatorial tris, so all these indices should be
      // positive.
      assert(Conf2.getFaceVert(i, 0) >= 0);
      assert(Conf2.getFaceVert(i, 1) >= 0);
      assert(Conf2.getFaceVert(i, 2) >= 0);
      Vert* pVert0 = ES2.getVert(
	  (Conf2.getFaceVert(i, 0) + bestPerm2) % numOrigVerts2);
      Vert* pVert1 = ES2.getVert(
	  (Conf2.getFaceVert(i, 1) + bestPerm2) % numOrigVerts2);
      Vert* pVert2 = ES2.getVert(
	  (Conf2.getFaceVert(i, 2) + bestPerm2) % numOrigVerts2);

      // Calculate the orientation exactly once.
      if (iOrient == 0)
	iOrient = checkOrient3D(pVert0, pVert1, pVert2, northVert2);

      if (iOrient == 1) {
	assert(checkOrient3D(pVert0, pVert1, pVert2, northVert2) == 1);
	assert(checkOrient3D(pVert0, pVert2, pVert1, southVert2) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, northVert2,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, southVert2,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
      else {
	assert(checkOrient3D(pVert0, pVert2, pVert1, northVert2) == 1);
	assert(checkOrient3D(pVert0, pVert1, pVert2, southVert2) == 1);
	bool qExist;
	Cell *pC;
	pC = m_pVM->createTetCell(qExist, pVert0, pVert2, pVert1, northVert2,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
	pC = m_pVM->createTetCell(qExist, pVert0, pVert1, pVert2, southVert2,
				  iReg, OMP_Flag, ID);
	assert(pC->isValid() && !qExist);
      }
    }

    // Replace the two bdry faces. Shouldn't matter which edge I take them from
    BFace* oldBFace0 = ES1.getBFace(0);
    BFace* oldBFace1 = ES1.getBFace(1);

    // Retrieve the faces that were just created on the boundary.
    Face *faceOnBdry0, *faceOnBdry1;
    {
      const Vert *pV0 = ES1.getVert(0);
      const Vert *pV1 = ES1.getVert(numOrigVerts1 - 1);

      faceOnBdry0 = findCommonFace(pV0, pV1, northVert1, NULL, true);
      assert(faceOnBdry0->isValid());
      // Just created them; they'd -better- exist!
      // Both Cells should exist
      assert(
	  faceOnBdry0->getLeftCell()->isValid()
	      && faceOnBdry0->getRightCell()->isValid());

      faceOnBdry1 = findCommonFace(pV1, pV0, southVert1, NULL, true);
      assert(faceOnBdry1->isValid());

      assert(
	  faceOnBdry1->getLeftCell()->isValid()
	      && faceOnBdry1->getRightCell()->isValid());
    }

    // Attach the new faces to the bdry faces.
    BFace* newBFace0;
    BFace* newBFace1;
    if (!OMP_Flag) {
      newBFace0 = m_pVM->createBFace(faceOnBdry0, oldBFace0);
      newBFace1 = m_pVM->createBFace(faceOnBdry1, oldBFace1);

      m_pVM->deleteBFace(oldBFace0);
      m_pVM->deleteBFace(oldBFace1);
    }
    else {
      newBFace0 = m_pVM->createBFaceParallel(faceOnBdry0, oldBFace0, ID);
      newBFace1 = m_pVM->createBFaceParallel(faceOnBdry1, oldBFace1, ID);

      m_pVM->deleteBFaceParallel(oldBFace0, ID);
      m_pVM->deleteBFaceParallel(oldBFace1, ID);
    }
    assert(newBFace0->isValid());
    assert(newBFace1->isValid());

    // Delete the old faces on the boundary.
  }
} // namespace GRUMMP
