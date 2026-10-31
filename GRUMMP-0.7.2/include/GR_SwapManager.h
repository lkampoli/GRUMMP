#ifndef GR_SwapManager_h
#define GR_SwapManager_h

#include <set>
#include <fstream>

#include "GR_Cell.h"
#include "GR_EdgeSwapInfo.h"
#include "GR_Face.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Mesh2D.h"
#include "GR_Observer.h"
//#include "GR_SwapDecider.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"
#include "GR_SurfMesh.h"
namespace GRUMMP
{
  class SwapDecider2D;
  class SwapDecider3D;
  class SwapManager2D : public Observer {
    SwapDecider2D* m_pSD;
  protected:
    Mesh2D* m_pM2D;
  private:
    GR_index_t m_attempts, m_successes;
    GR_index_t m_attemptsThreads[NUM_PROCS][20];
    GR_index_t m_successesThreads[NUM_PROCS][20];

    bool parallelswapping;
    /// Timer Variables
    double TotalSwappingTime;
    double Sync_Entry_Time;
    double Sync_Exit_Time;
    double Quary_Time;
    double Time_Spent_in_Critical;
    double SendEventsTime;
    double swapAllQueuedFacesTime;
    double swapAllFacess_Time;
    double Parallel_Time;
    double RemovingConflicts;
    double SwappingCheckTime[NUM_PROCS][100];
    double CriticalTime[NUM_PROCS][100];
    double FQMThread[NUM_PROCS][100];
    double ReconTime[NUM_PROCS][100];
    double CreatedCells;
    int numSwapsG[NUM_PROCS][200];
    ///

    /// Default construction is a really, really bad idea.
    SwapManager2D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SwapManager2D(SwapManager2D&);
    SwapManager2D&
    operator=(SwapManager2D&);
    // Need to have a queue and ways to access it, but that can wait for
    // specializations --- OpenMP, MPI, etc --- and/or Observer giving
    // feedback about faces that should be looked at.
    //
    // The following set is the simple serial version.  For parallel
    // (either version) it isn't going to be adequate.
    std::set<Face*> faceQueue;
    std::vector<Face*> SwapAllowedFaceQueue; //**********NEW for OMP************
    std::vector<Face*> NConflictFaceQueue;   //**********NEW for OMP************
    std::vector<Face*> RemainingFaces[NUM_PROCS][200];
    std::vector<Face*> faceQueueVecs[NUM_PROCS];
    template<class myType>
      void
      LoadBalancing(std::vector<myType> VecToBalanced[], int rows);
  public:
    /// Set the swap criterion and mesh at construction time.
    SwapManager2D(SwapDecider2D* pSD, Mesh2D* pMeshToSwap) :
	m_pSD(pSD), m_pM2D(pMeshToSwap), m_attempts(0), m_successes(0), faceQueue()
    {
      m_pM2D->addObserver(this,
			  Observable::cellCreated | Observable::faceDeleted);
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	m_successesThreads[ID][0] = 0;
	m_attemptsThreads[ID][0] = 0;
	numSwapsG[ID][0] = 0;
	CriticalTime[ID][0] = 0;
	FQMThread[ID][0] = 0;
	SwappingCheckTime[ID][0] = 0;
	ReconTime[ID][0] = 0;
      }
      parallelswapping = false;
      TotalSwappingTime = 0;
      Sync_Entry_Time = 0;
      Sync_Exit_Time = 0;
      Quary_Time = 0;
      Time_Spent_in_Critical = 0;
      SendEventsTime = 0;
      swapAllQueuedFacesTime = 0;
      swapAllFacess_Time = 0;
      Parallel_Time = 0;
      RemovingConflicts = 0;
      CreatedCells = 0;
    }
    virtual
    ~SwapManager2D();

    SwapDecider2D *
    changeSwapDecider(SwapDecider2D * swapDec)
    {
      SwapDecider2D * retVal = m_pSD;
      m_pSD = swapDec;
      return retVal;
    }
    // Observer function overrides
    void
    receiveDeletedFaces(std::vector<Face*>& deletedFaces);
    void
    receiveCreatedCells(std::vector<Cell*>& createdCells);
    void
    receiveMovedVerts(std::vector<Vert*>& movedVerts);
    void
    Set_Swap2D_Timer()
    {
      TotalSwappingTime = 0;
      Sync_Entry_Time = 0;
      Sync_Exit_Time = 0;
      Quary_Time = 0;
      Time_Spent_in_Critical = 0;
      SendEventsTime = 0;
      swapAllQueuedFacesTime = 0;
      swapAllFacess_Time = 0;
      Parallel_Time = 0;
      RemovingConflicts = 0;
      CreatedCells = 0;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	numSwapsG[ID][0] = 0;
	CriticalTime[ID][0] = 0;
	FQMThread[ID][0] = 0;
	SwappingCheckTime[ID][0] = 0;
	ReconTime[ID][0] = 0;
      }
    }
    void
    Print_Swap2D_Timer(int n = 2)
    {
      logMessage(n, "The time elapsed in Swapping is %3.10f\n",
		 TotalSwappingTime - SendEventsTime);
      logMessage(n, "The time elapsed in initial faceQueue is %3.10f\n",
		 swapAllFacess_Time);
      logMessage(n, "The time elapsed in Sync_Entry is %3.10f\n",
		 Sync_Entry_Time);
      logMessage(n, "The time elapsed in Sync_Exit is %3.10f\n",
		 Sync_Exit_Time);
      logMessage(n, "The time elapsed in FaceQueueManager function is %3.10f\n",
		 Quary_Time);
//		logMessage(n,"The time elapsed in Critical Section of FaceQueueManager function is %3.10f\n",Time_Spent_in_Critical);
      logMessage(n, "The time elapsed in sendEvents function is %3.10f\n",
		 SendEventsTime);
      logMessage(n, "The time elapsed in Reconfiguration is %3.10f\n",
		 swapAllQueuedFacesTime);
      logMessage(n, "The time elapsed in Parallel Regions is %3.10f\n",
		 Parallel_Time);
      logMessage(n, "The time elapsed in Sending created cells is %3.10f\n",
		 CreatedCells);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(n, "The number of faces swapped by thread %d : %d faces.\n",
		   ID, numSwapsG[ID][0]);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n,
	    "The time elapsed in CT and NCT for RC by thread %d are %3.10f && %3.10f && %3.10f\n",
	    ID, CriticalTime[ID][0], SwappingCheckTime[ID][0],
	    FQMThread[ID][0]);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n, "The time spent by thread %d for reconfiguration is %3.10f\n",
	    ID, ReconTime[ID][0]);
//		logMessage(n,"The time elapsed in Removing Conflicts from SwapAllowedFaceQueue is %3.10f\n",RemovingConflicts);
//		logMessage(n,"The time elapsed in checking which faces are eligible for swapping in faceQueue is %3.10f\n",SwappingCheckTime);
    }

    // If a function is added to change swap criteria, that's fine, but
    // that function must be sure to add all (interior) faces to the
    // swap queue -before- it returns.
    int
    swapAllFaces(bool OMP_Flag = false);
    int
    swapAllQueuedFaces(bool OMP_Flag = false);
    int
    swapAllQueuedFacesRecursive(bool OMP_Flag = false, bool SwapAll = false);
    int
    swapFace(FaceSwapInfo2D& FSI2D);
    void
    clearQueue()
    {
      faceQueue.clear();
      SwapAllowedFaceQueue.clear();	//**********NEW for OMP************
      NConflictFaceQueue.clear();
      //RemainingFaces.clear();
      for (int i = 0; i < NUM_PROCS; i++) {
	faceQueueVecs[i].clear();
	RemainingFaces[i][0].clear();
      }
    }					//**********NEW for OMP************
    std::vector<FaceSwapInfo2D> FaceStructs[NUM_PROCS];

    virtual bool
    reconfigure(FaceSwapInfo2D& FSI2D, int ID = 0, bool OMP_Flag = false) const;
    virtual bool
    reconfigure(FaceSwapInfo2D& FSI2D, Cell* pTri[], Cell *pC0, Cell *pC1,
		Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3) const;
    double
    faceQueueManager(int passes, bool SwapAll);	//**********NEW for OMP************
  };

  class SwapManagerSurf : public SwapManager2D {
    /// Default construction is a really, really bad idea.
    SwapManagerSurf();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SwapManagerSurf(SwapManagerSurf&);
    SwapManagerSurf&
    operator=(SwapManagerSurf&);
    // Need to have a queue and ways to access it, but that can wait for
    // specializations --- OpenMP, MPI, etc --- and/or Observer giving
    // feedback about faces that should be looked at.
    //

  public:

    /// Set the swap criterion and mesh at construction time.
    SwapManagerSurf(SwapDecider2D* pSD, Mesh2D* pMeshToSwap) :
	SwapManager2D(pSD, pMeshToSwap)
    {
    }
    ;

    virtual
    ~SwapManagerSurf()
    {
    }
    ;
    // only two functions I need;
    virtual bool
    reconfigure(FaceSwapInfo2D& FSI2D, int ID = 0, bool OMP_Flag = false) const;
    virtual bool
    reconfigure(FaceSwapInfo2D& FSI2D, Cell* pTri[], Cell *pC0, Cell *pC1,
		Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3) const;
  };

  class SwapManager3D : public Observer {
    SwapDecider3D* m_pSD;
    VolMesh* m_pVM;
    // The following set is going to eventually need to be replaced, at
    // least for some situations, with a more general queue.  If nothing
    // else, the set won't work well for parallel.
    std::set<Face*> faceQueue;
    std::vector<Face*> faceQueueVecs[NUM_PROCS];
    std::vector<Face*> RemainingFaces[NUM_PROCS];
    std::vector<Face*> NConflictFaceQueue[NUM_PROCS];
    /// Keep track of the number of swaps attempted and succeeded.
    GR_index_t m_edgeAttempts[11], m_edgeSuccesses[11];
    GR_index_t m_faceAttempts[11], m_faceSuccesses[11];

    GR_index_t m_edgeAttemptsThreads[11][NUM_PROCS],
	m_edgeSuccessesThreads[11][NUM_PROCS];
    GR_index_t m_faceAttemptsThreads[11][NUM_PROCS],
	m_faceSuccessesThreads[11][NUM_PROCS];
    int NumSwapps[NUM_PROCS];

    struct FaceEdgeInfo {
      FaceSwapInfo3D FCS;
      EdgeSwapInfo ESS;
      //Face* face = FCS.getFace();
      bool FaceOrEdge;
    };
    std::vector<FaceEdgeInfo> FaceStructs[NUM_PROCS];

    /// Time Variables
    double Sync_Entry_Time;
    double Sync_Exit_Time;
    double swapAllFacess_Time;
    double Parallel_Time;
    double SendEventsTime;
    double Quary_Time;
    double swapAllQueuedFacesTime;
    double FaceQueueManagerTime;
    double CriticalTime[NUM_PROCS][100];
    double FQMThread[NUM_PROCS][100];
    double ReconTime[NUM_PROCS][100];
    double NotCT[NUM_PROCS][100];
    double TotalSwappingTime;
    ///
    bool parallelswapping;
    bool SwapInsertor;

    /// Default construction is a really, really bad idea.
    SwapManager3D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    template<class myType>
      void
      LoadBalancing(std::vector<myType> VecToBalanced[], int rows);
    SwapManager3D(SwapManager3D&);
    SwapManager3D&
    operator=(SwapManager3D&);
    // Need to have a queue and ways to access it, but that can wait for
    // specializations --- OpenMP, MPI, etc --- and/or Observer giving
    // feedback about faces that should be looked at.
  public:
    /// Set the swap criterion and mesh at construction time.
    SwapManager3D(SwapDecider3D* pSD, VolMesh* pMeshToSwap) :
	m_pSD(pSD), m_pVM(pMeshToSwap), faceQueue(), parallelswapping(false), SwapInsertor(
	    false)
    {
      m_pVM->addObserver(
	  this,
	  Observable::cellCreated | Observable::faceDeleted
	      | Observable::vertMoved);
      for (int i = 0; i < 11; i++) {
	m_edgeAttempts[i] = 0;
	m_edgeSuccesses[i] = 0;
	m_faceAttempts[i] = 0;
	m_faceSuccesses[i] = 0;
      }
      for (int i = 0; i < 11; i++) {
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  m_edgeAttemptsThreads[i][ID] = 0;
	  m_edgeSuccessesThreads[i][ID] = 0;
	  m_faceAttemptsThreads[i][ID] = 0;
	  m_faceSuccessesThreads[i][ID] = 0;
	}
      }
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	NumSwapps[ID] = 0;
	CriticalTime[ID][0] = 0;
	FQMThread[ID][0] = 0;
	ReconTime[ID][0] = 0;
	NotCT[ID][0] = 0;
      }

      Sync_Entry_Time = 0;
      Sync_Exit_Time = 0;
      swapAllFacess_Time = 0;
      Parallel_Time = 0;
      SendEventsTime = 0;
      Quary_Time = 0;
      swapAllQueuedFacesTime = 0;
      FaceQueueManagerTime = 0;
      TotalSwappingTime = 0;
      parallelswapping = false;
    }
    void
    Set_Swap3D_Timer()
    {
      Sync_Entry_Time = 0;
      Sync_Exit_Time = 0;
      swapAllFacess_Time = 0;
      Parallel_Time = 0;
      SendEventsTime = 0;
      Quary_Time = 0;
      swapAllQueuedFacesTime = 0;
      FaceQueueManagerTime = 0;
      TotalSwappingTime = 0;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	NumSwapps[ID] = 0;
	CriticalTime[ID][0] = 0;
	FQMThread[ID][0] = 0;
	ReconTime[ID][0] = 0;
	NotCT[ID][0] = 0;
      }
    }
    void
    Print_Swap3D_Timer(int n = 2)
    {
      logMessage(n, "The time elapsed in Swapping is %3.10f\n",
		 TotalSwappingTime - SendEventsTime);
      logMessage(n, "The time elapsed in Reconfiguration is %3.10f\n",
		 Parallel_Time);
      logMessage(n, "The time elapsed in faceQueueManager function is %3.10f\n",
		 FaceQueueManagerTime);
      logMessage(n, "The time elapsed in SendEvents function is %3.10f\n",
		 SendEventsTime);
      logMessage(n, "The time elapsed in Sync_Entry function is %3.10f\n",
		 Sync_Entry_Time);
      logMessage(n, "The time elapsed in Sync_Exit function is %3.10f\n",
		 Sync_Exit_Time);
      logMessage(n, "The time elapsed in initial query is %3.10f\n",
		 Quary_Time);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(n, "The number of faces swapped by thread %d : %d faces.\n",
		   ID, NumSwapps[ID]);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n,
	    "The time elapsed in CT for RC by thread %d are %3.10f && %3.10f\n",
	    ID, CriticalTime[ID][0], FQMThread[ID][0]);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n, "The time spent by thread %d for reconfiguration is %3.10f\n",
	    ID, ReconTime[ID][0]);
    }
    ~SwapManager3D();

    SwapDecider3D *
    changeSwapDecider(SwapDecider3D * swapDec)
    {
      SwapDecider3D * retVal = m_pSD;
      m_pSD = swapDec;
      return retVal;
    }

    // TODO  The following is a hack and should be removed once it's no longer
    // useful to prevent a double free on exit after Watson insertion.
    void
    detachFromObservable()
    {
      m_pVM = NULL;
    } // No longer usable; destroyed later

      // Observer function overrides
    void
    receiveDeletedFaces(std::vector<Face*>& deletedFaces);
    void
    receiveCreatedCells(std::vector<Cell*>& createdCells);
    void
    receiveMovedVerts(std::vector<Vert*>& movedVerts);

    // If a function is added to change swap criteria, that's fine, but
    // that function must be sure to add all (interior) faces to the
    // swap queue -before- it returns.

    int
    swapAllFaces(bool OMP_Flag = false);
    int
    swapAllQueuedFaces(bool OMP_Flag = false);
    int
    swapAllQueuedFacesRecursive(bool OMP_Flag = false, bool SwapAll = false);
    int
    swapAllBdryEdges(bool OMP_Flag = false);
    bool
    swapBdryEdge(BFace * bface, Vert * v0, Vert * v1, Vert * vOther);
    int
    swapFace(Face* face, int ID = 0, bool OMP_Flag = false);
    int
    swapFace(const FaceSwapInfo3D& FC, int ID = 0, bool OMP_Flag = false);
    int
    swapFace(const FaceEdgeInfo& FEI, int ID, bool OMP_Flag = false);
    void
    clearQueue()
    {
      faceQueue.clear();
    }
    void
    AddFacetoQueue(Face *pF)
    {
      faceQueue.insert(pF);
    }
    bool
    reconfigure(const FaceSwapInfo3D& FC, int ID = 0, bool OMP_Flag = false);
    void
    reconfigureEdge(const EdgeSwapInfo& ES, int ID = 0, bool OMP_Flag = false);
    void
    reconfigureIntBdryEdge(const EdgeSwapInfo& ES1, const EdgeSwapInfo& ES2,
			   int ID = 0, bool OMP_Flag = false);
    void
    faceQueueManager3D();
  private:
    // These are the functions that do the actual reconfiguration.
    bool
    reconfTet22(Cell* pTets[4], Vert * const pVVertA, Vert* const pVVertB,
		Vert * const pVVertC, Vert * const pVVertD,
		Vert * const pVVertE, Vert * const pVVertF,
		TriBFace * const pTBFA, TriBFace * const pTBFB, int ID = 0,
		bool OMP_Flag = false);
    bool
    reconfTet23(Cell* pTets[4], Vert * const pVVertA, Vert* const pVVertB,
		Vert * const pVVertC, Vert * const pVVertD,
		Vert * const pVVertE, Vert * const pVVertF, int ID = 0,
		bool OMP_Flag = false);
    bool
    reconfTet32(Cell* pTets[4], Vert * const pVVertA, Vert* const pVVertB,
		Vert * const pVVertC, Vert * const pVVertD,
		Vert * const pVVertE, Vert * const pVVertF, int ID = 0,
		bool OMP_Flag = false);
    bool
    reconfTet44(Cell* pTets[4], Vert * const pVVertA, Vert* const pVVertB,
		Vert * const pVVertC, Vert * const pVVertD,
		Vert * const pVVertE, Vert * const pVVertF, int ID = 0,
		bool OMP_Flag = false);
  };

}
;

#endif
