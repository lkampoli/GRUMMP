#ifndef GR_REFINEMENTMANAGER3D_H
#define GR_REFINEMENTMANAGER3D_H 1
#include "GR_BaseQueueEntry.h"
#include "GR_BFace.h"
#include "GR_Face.h"
#include "GR_SwapManager.h"
#include "GR_SmoothingManager.h"
#include "GR_SwapDecider.h"
#include "GR_Subseg.h"
#include "GR_TMOPQual.h"
#include "GR_TreeConstructor.h"

#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <functional>
#include <tr1/unordered_map>

//Class in which all the functions to compute
//a cell's priority are defined. Should be easily expandable.

class OldPriorityCalculator {

public:

  enum PrioType {
    EDGE_RADIUS,
    SIN_DIHED,
    BIAS_SIN_DIHED,
    VOLUME_LENGTH,
    MIN_SOLID,
    NORMALIZED_SHAPE_RATIO,
    MEAN_RATIO,
    VOLUME_AREA_RATIO,
    INRAD_ALTITUDE,
    INRAD_LENGTH,
    LONGEDGE_MINALT,
	SIZE_ONLY,
    UNKNOWN
  };

private:

  PrioType m_type;

  //The object to compute the length scale.
  //Instantiated elsewhere, pointer is passed to the constructor.
  std::shared_ptr<GRUMMP::Length> m_length;

  OldPriorityCalculator(const OldPriorityCalculator&) :
      m_type(UNKNOWN), m_length(NULL)
  {
    assert(0);
  }

  OldPriorityCalculator&
  operator=(const OldPriorityCalculator&)
  {
    assert(0);
    return *this;
  }

public:

  OldPriorityCalculator(PrioType type = UNKNOWN) :
      m_type(type), m_length(NULL)
  {
  }

  virtual
  ~OldPriorityCalculator();

  void
  set_length_calculator(std::shared_ptr<GRUMMP::Length> const length);

  void
  set_priority_type(PrioType type)
  {
    m_type = type;
  }

  //Return the shape priority and the length priority separately.

  void
  compute_priority(const Cell* const cell, double& length_prio,
		   double& shape_prio) const;

  //This function is set up to return:
  // . 0 if the cell is small enough with respect to the sizing field.
  // . -LARGE_DBL if the cell MUST NOT BE SPLIT (useful for small angles)
  // . A positive size priority if the cell must be split.

  double
  compute_length_priority(const Cell* const cell) const;

  // computes the shape part..
  double
  compute_shape_priority(const Cell* const cell) const;

  //Static functions computing different quality measures for a tet.

  static double
  edge_radius_ratio(const Cell* const cell);

  static double
  min_sin_dihed(const Cell* const cell);

  static double
  bias_min_sin_dihed(const Cell* const cell);

  static double
  volume_length(const Cell* const cell);

  static double
  min_solid(const Cell* const cell);

  static double
  normalized_shape_ratio(const Cell* const cell);

  static double
  mean_ratio(const Cell* const cell);

  static double
  volume_area_ratio(const Cell* const cell);

  static double
  inrad_altitude_ratio(const Cell* const cell);

  static double
  inrad_longest_edge_ratio(const Cell* const cell);

  static double
  long_edge_min_altitude_ratio(const Cell* const cell);
};

typedef std::priority_queue<CellQueueEntry*, std::vector<CellQueueEntry*>,
    CellQueueEntry::PrioOrder> RefinementPrioQueue;

//Priority queue to control refinement.
class RefinementQueue : public RefinementPrioQueue {

  typedef std::vector<CellQueueEntry*> EntryVec;

  OldPriorityCalculator* m_calculator;

  double m_min_shape_priority, m_max_shape_priority;

  // Never defined
  RefinementQueue(const RefinementQueue&);
  RefinementQueue&
  operator=(const RefinementQueue&);

public:

  RefinementQueue() :
      RefinementPrioQueue(), m_calculator(
	  new OldPriorityCalculator(OldPriorityCalculator::UNKNOWN)), m_min_shape_priority(
	  -LARGE_DBL), m_max_shape_priority(LARGE_DBL)
  {
  }

  RefinementQueue(OldPriorityCalculator::PrioType prio_type, bool max_prio_top =
		      true,
		  double min_shape_priority = 0., double max_shape_priority =
		  LARGE_DBL) :
      RefinementPrioQueue(
	  max_prio_top ?
	      CellQueueEntry::PrioOrder() : CellQueueEntry::PrioOrder(true)), m_calculator(
	  new OldPriorityCalculator(prio_type)), m_min_shape_priority(
	  min_shape_priority), m_max_shape_priority(max_shape_priority)
  {
  }

  virtual
  ~RefinementQueue();

  //Access to the vector upon which the queue is built.
  //EntryVec& impl() const { return c; }

  //Returns the first UNDELETED cell at the top of the queue.
  //All the deleted entries at the top will get removed.
  //Returns NULL if the queue no longer contains valid cells.
  Cell*
  get_top_valid_cell();

  //Pops the top entry from the queue. Returns false is the queue is empty.
  bool
  pop_top_entry();

  //Adds a cell to the queue.
  //(builds a CellQueueEntry and computes priority in the background)
  //Returns true if the cell is actually inserted (it is within the
  //refinement bounds of the refinement criterion).
  void
  add_entry(Cell* const cell);
  // used for quality refinement, to ensure we only refine for quality,
  // ignoring length scale
  void
  add_quality_entry(Cell* const cell);

  //Removes all deleted entries from the queue.
  void
  clean_queue();

  //Empties the queue completely
  void
  empty_queue()
  {
    while (pop_top_entry()) {
    }
    assert(empty());
  }

  //
  void
  set_priority_type(OldPriorityCalculator::PrioType type)
  {
    assert(m_calculator);
    m_calculator->set_priority_type(type);
  }

  void
  set_min_shape_priority(double min_shape_priority)
  {
    assert(min_shape_priority >= 0.);
    m_min_shape_priority = min_shape_priority;
  }

  void
  set_max_shape_priority(double max_shape_priority)
  {
    assert(max_shape_priority >= 0.);
    m_max_shape_priority = max_shape_priority;
  }

  void
  set_length_calculator(std::shared_ptr<GRUMMP::Length> const length)
  {
    assert(length && m_calculator);
    m_calculator->set_length_calculator(length);
  }

  //Prints the content of the queue and aborts.
  void
  print_queue();

  int
  queue_size()
  {
    return size();
  }
  ;

};

namespace GRUMMP
{
/// This is the class formerly known as TetMeshRefiner, renamed to fit in with
/// everything else.
  class RefinementManager3D : public Observer {

    enum SplitType {
      NONE, CELL, LOCKED
    };
    //The tetrahedral mesh to refine.
    VolMesh* m_mesh;

    // InsertionPointCalculator
    InsertionPointCalculator *m_IPC;

    // Watson Inserter
    WatsonInserter3D* m_WI3D;

    // Insertion Swapper and other stuff
    // Used when Watson insertion gives a bad hull
    // Always inserts and swaps for delaunay
    DelaunaySwapDecider3D *m_Del3D;
    SwapManager3D *m_DelSwap;
    SwappingInserter3D * m_SI3D;

    //The length scale calculator.
    std::shared_ptr<GRUMMP::Length> m_isoLength;

    // Thrown in here due to queuing issues.
    // Normal Queue does size + quality queuing
    bool m_refineForQuality;

    // The queues "controlling" the refinement.
    std::deque<Subseg*> m_subseg_queue;
    std::deque<InsertionQueueEntry> m_bdryFacesToSplit;
    RefinementQueue* m_cell_queue;

    /// When refining for size, queue boundary faces for size only
    bool m_refineBdryForSize;

    // Sets of unsplittable cells (point encroaches)
    // boundary faces/subsegs (swapping inserter fails due to curves)
    // subsegs
    std::set<InsertionQueueEntry> m_lockedCells;
    std::set<InsertionQueueEntry> m_lockedBdryFaces;
    std::set<Subseg*> m_lockedSubsegs;

    //Octree structure:
    GRUMMP::OctTree *OctT;

    //***** Parallel Varaibles *******
    bool IsInParallel;
    Rect3D* ThreadLeaves[NUM_PROCS];
    std::set<Rect3D*> IndependentSet;

    WatsonInserter3D* Para_m_WI3D[NUM_PROCS];
    SwappingInserter3D * Para_m_SI3D[NUM_PROCS];
    std::deque<Subseg*> Para_m_subseg_queue[NUM_PROCS];
    std::deque<InsertionQueueEntry> Para_m_bdryFacesToSplit[NUM_PROCS];
    RefinementQueue* Para_m_cell_queue[NUM_PROCS];
//	Length3D* Para_m_length[NUM_PROCS];
    // These containers are introduced to store all the entities ...
    // that have been inherited from other places and at the time of ...
    // sending events, each thread could not put it into its working leaves.
    std::vector<TetCell*> Remaining_Cells;
    std::deque<BFace*> Remaining_m_bdryFacesToSplit;
    std::deque<Subseg*> Remaining_m_subseg;
    std::vector<BFace*> Para_m_newBdryFaces[NUM_PROCS];

    std::set<InsertionQueueEntry> Para_m_lockedCells[NUM_PROCS];
    std::set<InsertionQueueEntry> Para_m_lockedBdryFaces[NUM_PROCS];
    std::set<Subseg*> Para_m_lockedSubsegs[NUM_PROCS];

    GR_index_t NumCreatedCells;
    omp_lock_t *lck;
    bool *CheckLock;
    int FilterLevel;
    int Granularity;
    int Filter;
    int ThreadID;
    bool OKTOGO;

    double TotalRefineTime;
    double JustParallelTime;
    double SendEventsTime;
    double SendCCTime;
    double SendCBTime;
    double CumWRTime;
    double SecondPassesTime;
    double FIndSetsTime;
    double ThreadTime[NUM_PROCS][100];
    int ThreadInNum[NUM_PROCS][100];

    //Make these private (for now).
    // These next three functions are deliberately never defined
    RefinementManager3D();
    RefinementManager3D(const RefinementManager3D&);
    RefinementManager3D&
    operator=(const RefinementManager3D&);

    //To initialize the queues with encroached boundary and
    //cells that do not meet refinement criteria.
    void
    queueEncroachedSubsegs();
    void
    queueEncroachedBdryFaces();
    void
    queueAllCells(bool OMP_FLAG = false, int NQLayers = 3);
    GR_index_t
    queueAllCellsThreadWise(int ID, int Filter);

    //Clears the deleted entries at the top of m_subseg_queue and m_subface_queue
    //Returns false if the queue is empty after clearing, true otherwise.
    bool
    clearTopSubsegQueue();
    bool
    clearTopBdryFaceQueue();

  public:

    //The constructor...
    RefinementManager3D(VolMesh* const mesh,
			InsertionPointCalculator * const IPC,
			const std::set<Subseg*>* const subsegs = NULL,
			const eEncroachType encroach_type = eBall,
			std::shared_ptr<GRUMMP::Length> length3D = NULL,
			const OldPriorityCalculator::PrioType =
			    OldPriorityCalculator::EDGE_RADIUS);

    //... and destructor.
    ~RefinementManager3D();

    // left from SG's work. could be useful so keeping it, but not
    // actively called
    void
    printQualityData() const;

    //Refine the mesh for quality (and size).
    void
    refineForQuality(bool Para_Flag = false, int NOL = 3);

    void
    refineMesh(bool Para_Flag = false, int NOL = 3);

    // Similar to 2D version, used for surface insertion
    void
    refineAllQueuedCellsForQuality(double angle);

    // This refines all the boundaries, subsegs + bfaces
    // If another mesh is given, use that for the length scale
    int
    refineBoundaryForLength();

    // queue won't queue based on length scale, ONLY based on quality
    // dangerous, but used for adaptation.
    void
    setQueueForQuality()
    {
      m_refineForQuality = true;
    }

    inline void
    Set_Refine3D_Timer()
    {
      TotalRefineTime = 0;
      JustParallelTime = 0;
      SendEventsTime = 0;
      SendCCTime = 0;
      SendCBTime = 0;
      CumWRTime = 0;
      SecondPassesTime = 0;
      FIndSetsTime = 0;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	ThreadTime[ID][100] = 0;
	ThreadInNum[ID][100] = 0;
      }
    }
    inline void
    Print_Refine3D_Timer(int n = 2)
    {
      logMessage(n, "The total time elapsed in Refining is %3.10f\n",
		 TotalRefineTime - SendEventsTime);
      logMessage(
	  n,
	  "The time elapsed in all parallel portions of the code is %3.10f\n",
	  JustParallelTime);
      logMessage(n, "The time elapsed in sending events is %3.10f\n",
		 SendEventsTime);
      logMessage(n, "The time elapsed in sending Created Cells is %3.10f\n",
		 SendCCTime);
      logMessage(
	  n, "The time elapsed in sending Created Boundary Faces is %3.10f\n",
	  SendCBTime);
      logMessage(n, "The time elapsed in calling WhichRect Faces is %3.10f\n",
		 CumWRTime);
      logMessage(
	  n,
	  "The time elapsed in creating second layers of OctTree is %3.10f\n",
	  SecondPassesTime);
      logMessage(n, "The time elapsed in finding independent sets is %3.10f\n",
		 FIndSetsTime);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n,
	    "The time elapsed in thread %d : %3.10f seconds, inserted %d verticies.\n",
	    ID, ThreadTime[ID][100], ThreadInNum[ID][100]);
    }

  public:

    //Functions to split a mesh entity.
    // subseg and bdry face return true or false
    bool
    splitSubseg(Subseg* const subeg_to_split,
		std::pair<Subseg*, Subseg*>* new_subsegs = NULL);
    bool
    splitBdryFace(BFace* bface_to_split, bool splitAtEdge = false);

    // this can return LOCKED (stop trying to split due to encroachment)
    // CELL (success), or NONE (encroached)
    SplitType
    splitTet(Cell* const cellToSplit, bool splitAtCentroid = false);

    bool
    splitSubsegInParallel(Subseg* const subeg_to_split,
			  std::pair<Subseg*, Subseg*>* new_subsegs = NULL);
    bool
    splitBdryFaceInParallel(BFace* bface_to_split,
			    bool splitAtCentroid = false);
    SplitType
    splitTetInParallel(Cell* const cellToSplit, bool splitAtCentroid = false);

    void
    onlyRefineForQuality()
    {
      m_refineForQuality = true;
    }
    ;

  private:
    //Refines the boundary entities until they satisfy the sizing field.
    // This is called by the function below, which handles all the boundaries
    int
    refineSubsegsForLength();

    //Runs a refinement loop, inserts all vertices it can
    //until all queues are empty.
    //Returns total number of vertices inserted.
    int
    runRefinementLoop();
    // Running refinement loop by one thread for one leaf.
    int
    runRefinementLoopForEachLeaf(int ID, int FLevel, GR_index_t& RVerts);
    // Refine encroached subsegs and boundary edges before getting to the main loop of interior refinement.
    // to ease the problem of parallel refinement
    int
    runRefinementForEncroachedSubsegs();
    int
    runRefinementForEncroachedBdryFaces();

    //Check for encroachment of boundary entities.
    bool
    isBdryFaceEncroached(const BFace* const bface) const;
    bool
    isSubsegEncroached(const Subseg* const subseg) const;
    void
    checkNewBdryFacesForEncroachment();

    // Observable information
    void
    receiveCreatedCells(std::vector<Cell*>& createdCells);
    void
    receiveCreatedBFaces(std::vector<BFace*>& createdBFaces);

    //This function initializes m_subseg_map with either a set of
    //subsegs that was obtain somewhere else, or by infering the
    //subsegs from m_mesh (this last case should only be used with
    //GRUMMP's native geometry).
    void
    initializeSubsegs(const std::set<Subseg*>* const subsegs);

    //Finds the seed cells for a subsegment.
    bool
    findSeedCellsForSubseg(double * const newLoc, Subseg* const subseg,
			   std::vector<Cell*>& seed_cells,
			   std::vector<BFace*>& seed_bfaces);

  };

#ifdef HAVE_MESQUITE
  /*
   * Class for Aniso Edge Adaptation in 3D
   *
   * Extension of the 2D Aniso edge adaptation with TMOP. This is still an unfinished code.
   * The main problem of the algorithm was in swapping cells. Two different quality criterias are
   * used for swapping: SineDihedralAngles and VolumeLength3D. VolumeLentgh is used to force swapping cells
   * by improving the best cell into a better quality, which also creates a lot of bad cells. SineDihedral angles
   * is used to fix as much of the bad cells created by the previous swapper. This gave good alignment but only at surfaces,
   * and at interior cells it created a lot of bad cells.
   * For this reason both swapping criteria are set to improve the worst quality (as it should be... assuring bad cells are not created), but
   * as a results, good alignment isn't obtained even at the surface. A solution to fix this is still needed.
   */

  class AnisoEdgeRefinementManager3D : public Observer {

    enum eInterpType {
      eDefault, eAverage, eLinear, eLogEuclidean
    };
    enum eInsertType {
      eQualityEdge, eMetricLengthEdge, eCoarsening
    };

    TMOPQual *m_TQ;
    VolMesh* m_VMesh;
    double dMaxMLength;
    double dMinMLength;
    InsertionPointCalculator *m_IPC;
    GRUMMP::SmoothingManagerMesquite3D *m_pSMM;
    AnisoRefinement3D *m_pAR;
    GRUMMP::SwapManager3D *m_pSM3D;
    GRUMMP::SwapManager3D *m_pSM3D_2;
    SwappingInserter3D* m_SI3D;
    InsertionQueue m_IQBE;
    InsertionQueue m_IQIE;

    int iNMetricValues;
    bool bSmooth;
    bool bReceiveCreated = false;
    eInterpType m_Type = eDefault;

    AnisoEdgeRefinementManager3D();
    AnisoEdgeRefinementManager3D(const RefinementManager3D&);
    RefinementManager3D&
    operator=(const RefinementManager3D&);

  public:

    AnisoEdgeRefinementManager3D(VolMesh* meshToRefine, double dMax_MLength,
				 double dMin_MLength,
				 GRUMMP::InsertionPointCalculator * pIPC,
				 GRUMMP::SmoothingManagerMesquite3D *pSMM,
				 AnisoRefinement3D *pAR, SwapManager3D *pSM3D,
				 SwapManager3D *pSM3D_2,
				 SwappingInserter3D *pSI3D) :
	m_TQ(new TMOPQual), m_VMesh(meshToRefine), dMaxMLength(dMax_MLength), dMinMLength(
	    dMin_MLength), m_IPC(pIPC), m_pSMM(pSMM), m_pAR(pAR), m_pSM3D(
	    pSM3D), m_pSM3D_2(pSM3D_2), m_SI3D(pSI3D),

	m_IQBE(m_VMesh, InsertionQueueEntry::dInvalidPriority, true), m_IQIE(
	    m_VMesh, InsertionQueueEntry::dInvalidPriority, true), iNMetricValues(
	    6), bSmooth(false)
    {

      m_VMesh->addObserver(this,
			   Observable::bfaceCreated | Observable::faceCreated);
      m_TQ->SetMetricValue(6);
      m_TQ->CheckTypeCell(CellSkel::eTet);
      m_IQBE.vSetMinMLength(dMin_MLength);
      m_IQBE.vSetMaxMLength(dMax_MLength);
      m_IQBE.vSetMetricTol(1e-6);
      m_IQIE.vSetMinMLength(dMin_MLength);
      m_IQIE.vSetMaxMLength(dMax_MLength);
      m_IQIE.vSetMetricTol(1e-6);
      m_VMesh->allowBdryChanges();
      SetInterpType(eDefault);
    }

    ~AnisoEdgeRefinementManager3D()
    {
      delete m_SI3D;
      delete m_TQ;
      if (m_VMesh)
	m_VMesh->removeObserver(this);
    }

    void
    receiveCreatedFaces(std::vector<Face*>& createdFaces);
    void
    receiveCreatedBFaces(std::vector<BFace*>& createdBFaces);

//	Runs the algorithm to aniso adapt a mesh, a metric must be previously set in the mesh
    GR_index_t
    AnisoRefineMesh();

//	Modifies the mesh according to the queues
    GR_index_t
    AnisoRefineAllQueued();

//	Splits edges if the quality of the configuration after inserting a point is improved
    GR_index_t
    AnisoQualityEdgeInsertion();

//	Splits edges which metric length is bigger than a max metric length allowed
    GR_index_t
    AnisoMetricLengthEdgeInsertion();

//	Coarses the mesh based on metric length. Edges are queued for deleting if any of its vertices has at least two edges smaller
//	than a min metric length. Wich point to delete with edge collapse is decided by quality
    GR_index_t
    AnisoMetricEdgeCoarsening();

//	Uses Blobal Smoothing with mesquite based on targets defined by an eigen decomposition of a metric
    void
    GlobalSmoothing();

    void
    emptyQueue(InsertionQueue IQ)
    {
      IQ.emptyQueue();
      assert(IQ.iQueueLength() == 0);
    }

  private:

//	Decides where to split an edge, exact location is computed with InsertionPointCalculator
//	at the moment only insert point at midsections of an edge
    bool
    splitEdge(Vert * const pV0, Vert* const pV1);

//	Decides where to split a boundary edge, exact location is computed with InsertionPointCalculator
//	at the moment only insert point at midsections of an edge
    bool
    splitBdryEdge(Vert * const pV0, Vert* const pV1);

//	Functions add to queue all the edges to be modified by priority
    GR_index_t
    queueAll();
    GR_index_t
    queueEdges();
    GR_index_t
    queueBEdges();

//	Adds elements to the queue
    bool
    addToPriorityQueueEdge(Vert* const pV0, Vert* const pV1/*Edge*/);
    bool
    addToPriorityQueueBdryEdge(Vert* const pV0, Vert* const pV1/*Edge*/);

//	Functions used to interpolated metrics
    void
    InterpMetric(Vert* pV, std::set<Vert*> setpV, enum eInterpType Type =
		     eDefault);
    void
    AverageMetric(Vert* pV, std::set<Vert*> setpV);
    void
    LinearInterpMetric(Vert* pV, std::set<Vert*> setpV);
    void
    SetInterpType(eInterpType eIT)
    {
      m_Type = eIT;
    }
    ;

//	Sets priorities for each step of the process
    void
    SetQueuesPriority(enum eInsertType Type);

//	If set, evaluate all the new cells created
    void
    SetReceiveCreated(bool bRC)
    {
      bReceiveCreated = bRC;
    }
    ;

  };
#endif

}
#endif
