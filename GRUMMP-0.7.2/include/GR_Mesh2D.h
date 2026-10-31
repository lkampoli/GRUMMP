#ifndef GR_Mesh2D
#define GR_Mesh2D 1

#include "GR_config.h"
#include "GR_Classes.h"
#include "GR_BFace2D.h"
#include "GR_Cell.h"
#include "GR_CellCV.h"
#include "GR_Face.h"
#include "GR_GRCurve.h"
#include "GR_Mesh.h"
#include "GR_EntContainer.h"
#include <map>
#include <omp.h>

class GRGeom2D;

typedef struct _renum_ {
  int iV, iDeg;
  _renum_(const int iVIn = -1, const int iDegIn = 0) :
      iV(iVIn), iDeg(iDegIn)
  {
  }
} Renum;

/**\brief The 2D specialization of the Mesh base class.
 *
 * Functions inherited from Mesh will not be documented again here.
 */
class Mesh2D : public Mesh {
protected:
  /// Variable size array containing faces in the 2D mesh
  EntContainer<EdgeFace> m_ECEdgeF;

  /// Variable size array containing bdry faces in the 2D mesh
  EntContainer<BdryEdge> m_ECEdgeBF;

  /// Variable sized array containing triangular cells in the 2D mesh
  EntContainer<TriCell> m_ECTri;

  /// Variable sized array containing quadrilateral cells in the 2D mesh
  EntContainer<QuadCell> m_ECQuad;

  /// Variable size array containing internal bdry faces in the 2D mesh
  EntContainer<IntBdryEdge> m_ECIntEdgeBF;

private:
  /// The copy constructor currently appears to be unnecessary.
  Mesh2D(const Mesh2D&);

  /// operator= has probably always been unnecessary.
  Mesh2D&
  operator=(const Mesh2D&);

public:
  double
  Mesh2D_sync_entry(double& paralleltime, const bool IsSwapping)
  {
    double Sync_Entry_StartTime = omp_get_wtime();
    if (NUM_PROCS > 1) {
      paralleltime = 0;
      if (IsSwapping) {
	paralleltime += m_ECEdgeF.sync_entry(IsSwapping);
	paralleltime += m_ECTri.sync_entry(IsSwapping);
      }
      else {
	paralleltime += m_ECEdgeF.sync_entry();
	paralleltime += m_ECTri.sync_entry();
	//paralleltime += m_ECQuad.sync_entry();
	paralleltime += m_ECIntEdgeBF.sync_entry();
	paralleltime += m_ECVerts.sync_entry();
	paralleltime += m_ECEdgeBF.sync_entry();
      }
    }
    Sync_Entry_StartTime = omp_get_wtime() - Sync_Entry_StartTime;
    return Sync_Entry_StartTime;
  }

  double
  Mesh2D_sync_exit(bool IsSwapping)
  {
    double Sync_Exit_StartTime = 0;
    if (NUM_PROCS > 1) {
      if (IsSwapping) {
	Sync_Exit_StartTime = omp_get_wtime();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<EdgeFace*>::iterator iterS =
	      m_ECEdgeF.empties_pa[ID].begin();
	  std::list<EdgeFace*>::iterator iterE = m_ECEdgeF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECEdgeF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<TriCell*>::iterator iterS = m_ECTri.empties_pa[ID].begin();
	  std::list<TriCell*>::iterator iterE = m_ECTri.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTri.sync_exit();
      }
      else {
	Sync_Exit_StartTime = omp_get_wtime();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<EdgeFace*>::iterator iterS =
	      m_ECEdgeF.empties_pa[ID].begin();
	  std::list<EdgeFace*>::iterator iterE = m_ECEdgeF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECEdgeF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<TriCell*>::iterator iterS = m_ECTri.empties_pa[ID].begin();
	  std::list<TriCell*>::iterator iterE = m_ECTri.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTri.sync_exit();
//			  for (int ID=0 ; ID<omp_get_num_procs() ; ID++)
//		  	  {
//				  std::list<QuadCell*>::iterator iterS1=m_ECQuad.empties_pa[ID].begin();
//				  std::list<QuadCell*>::iterator iterE=m_ECQuad.empties_pa[ID].end();
//				  if (m_ECQuad.empties_pa[ID].empty()) assert(0);
//				  for ( ; iterS1 != iterE ; ++iterS1)
//				  {
////					  if (&(*iterS1) == NULL) assert(0);
////					  if (!(*iterS1)->isValid()) assert(0);
//					  (*iterS1)->markAsDeleted();
//				  }
//			  }
//			   m_ECQuad.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<IntBdryEdge*>::iterator iterS =
	      m_ECIntEdgeBF.empties_pa[ID].begin();
	  std::list<IntBdryEdge*>::iterator iterE =
	      m_ECIntEdgeBF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECIntEdgeBF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<Vert*>::iterator iterS = m_ECVerts.empties_pa[ID].begin();
	  std::list<Vert*>::iterator iterE = m_ECVerts.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECVerts.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<BdryEdge*>::iterator iterS =
	      m_ECEdgeBF.empties_pa[ID].begin();
	  std::list<BdryEdge*>::iterator iterE =
	      m_ECEdgeBF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECEdgeBF.sync_exit();
      }
      Sync_Exit_StartTime = omp_get_wtime() - Sync_Exit_StartTime;
    }
    return Sync_Exit_StartTime;
  }
  /// Read mesh data from a file.
  void
  readFromFile(const char strFileName[], int* numVars = nullptr, double** data =
		   nullptr);
  void
  readFromFileTerrible(const char strFileName[]);
  /// Associate this mesh with a geometry.
  void
  relateMeshToGeometry(GRGeom2D* const geometry);
  void PseudoFaces();
  ///Default constructor called by TriMeshBuilder.
  Mesh2D(const int iQualMeas = 2);

  ///This constructor is deprecated, but still required for 3D initial
  ///tetrahedralization until Serge's facet branch merges in -and- can
  ///read .bdry and .smesh files.
  Mesh2D(const Bdry2D& B2D, const int iQualMeas = 2);

  //Creating a Delaunay mesh from a set of vertices.
  //THIS SHOULD BE REIMPLEMENTED IN THE TriMeshBuilder CLASS.
  // This is only used in the scat2d executable, apparently, so it can wait.
  Mesh2D(const EntContainer<Vert>& EC);

  ///Destructor.
  virtual
  ~Mesh2D()
  {
  }

  /// Inquiry functions

  //////////////////////////////////////////////////////////////////////
  /// Functions that directly manipulate or query the mesh database. ///
  //////////////////////////////////////////////////////////////////////

  GR_index_t
  iVertIndex(const Vert* const pV) const
  {
    return (m_ECVerts.getIndex(pV));
  }
  GR_index_t
  iNumVerts() const
  {
    return m_ECVerts.lastEntry();
  }

  /// Given a face pointer, find the corresponding integer index.
  virtual GR_index_t
  getFaceIndex(const Face* const pF) const;
  virtual GR_index_t
  getNumFaces() const
  {
    return m_ECEdgeF.lastEntry();
  }
  virtual Face*
  getFace(const GR_index_t i) const;
  ///
  virtual GR_index_t
  getBFaceIndex(const BFace* const pBF) const;
  virtual GR_index_t
  getNumBdryFaces() const
  {
    return m_ECEdgeBF.lastEntry();
  }
  virtual GR_index_t
  getNumIntBdryFaces() const
  {
    return m_ECIntEdgeBF.lastEntry();
  }
  virtual BFace*
  getBFace(const GR_index_t i) const;
  virtual GR_index_t
  getNumBdryPatches()
  {
    assert(0);
    return (~0);/*B2D.iNumPatches();*/
  }
  virtual void
  buildBdryPatches();
  ///
  virtual GR_index_t
  getCellIndex(const Cell* const pC) const;
  GR_index_t
  getNumTriCells() const
  {
    return m_ECTri.lastEntry();
  }
  GR_index_t
  getNumQuadCells() const
  {
    return m_ECQuad.lastEntry();
  }
  virtual GR_index_t
  getNumCells() const
  {
    return getNumTriCells() + getNumQuadCells();
  }

  GR_index_t
  getNumEdgesInUse() const
  {
    return m_ECEdgeF.size();
  }
  GR_index_t
  getNumBdryEdgesInUse() const
  {
    return m_ECEdgeBF.size();
  }
  GR_index_t
  getNumIntBdryEdgesInUse() const
  {
    return m_ECIntEdgeBF.size();
  }
  GR_index_t
  getNumTrisInUse() const
  {
    return m_ECTri.size();
  }
  GR_index_t
  getNumQuadsInUse() const
  {
    return m_ECQuad.size();
  }
  GR_index_t
  getNumTetsInUse() const
  {
    return 0;
  }
  GR_index_t
  getNumPyrsInUse() const
  {
    return 0;
  }
  GR_index_t
  getNumPrismsInUse() const
  {
    return 0;
  }
  GR_index_t
  getNumHexesInUse() const
  {
    return 0;
  }

  virtual Cell*
  getCell(const GR_index_t i) const;
  QuadCell*
  getQuadCell(const GR_index_t i) const;
  ///
  virtual bool
  isValidEntHandle(void* pvEnt) const;

  ///
  virtual eMeshType
  getType() const
  {
    return Mesh::eMesh2D;
  }
  ///
  /// perimeter length or  area:
  virtual double
  calcBoundarySize() const;
  /// interior   area  or volume:
  virtual double
  calcInteriorSize() const;
  virtual bool
  isWatertight() const;
  ///
  /// Validity checking
  virtual bool
  isValid() const;

public:
  double gradeparameter = 1.;
  double determineparameter = 1.;
  double minangle = 0.447832396928932;
  ///
  //Determine whether the edge should be added to the edge-based queue.
  bool
  DeterGoodTri(TriCell* cell, Face* face);
  bool
  isLocallyDelaunay() const;
  /////////////////////////////////////////////////////////////////
  /// Primitive functions that modify mesh geometry or topology ///
  /////////////////////////////////////////////////////////////////
  // The generic, ITAPS-style calls
  /*   Face* createFace(Vert * const apV[], const int iNVerts); */
  /*   Cell* createCell(Face * const apF[], const int iNFaces, const int iReg); */
  /*   Cell* createCell(Vert * const apV[], const int iNFaces, const int iReg); */

  // The size-specific (and more efficient) calls; the ITAPS-style calls
  // will either call these or switch on number of verts and execute the
  // same code.
  virtual Vert*
  createVert(const double dX, const double dY, const double dZ = 0);
  virtual Vert*
  createVert(const double adCoords[]);
  BFace*
  createBFace(Face * const pF, BFace * const pBF = pBFInvalidBFace);
  /// This form is deprecated, but currently needed to read native format.
  BFace*
  createBFace(Face * const pF, const int i);
  Face*
  createFace(int minlayerindex, bool& qAlreadyExisted, Vert * const pV0,
	     Vert * const pV1, const bool qForceDuplicate = false, int ID = 0);
  Face*
  createFace(bool& qAlreadyExisted, Vert * const pV0, Vert * const pV1,
	     const bool qForceDuplicate = false, int ID = 0);
  /// This form will almost certainly be deprecated in future
  IntBdryEdge*
  createIntBdryEdge(Face * const pFA, Face * const pFB, BFace * const pBF =
  pBFInvalidBFace);
  /// Create a new internal bdry edge, including setting geometry info.
  IntBdryEdge*
  createIntBdryEdge(Vert* const pVBeg, Vert* const pVEnd, GRCurve* const curve,
		    double dBegParam, double dEndParam);
  /// Create a new bdry edge, including setting geometry info.
  BdryEdge*
  createBdryEdge(Vert* const pVBeg, Vert* const pVEnd, GRCurve* const curve,
		 double dBegParam, double dEndParam);
  virtual TriCell*
  createTriCell(Face * const pF0, Face * const pF1, Face * const pF2,
		const int iReg = iDefaultRegion, int ID = 0);
  QuadCell*
  createQuadCell(Face* const pF0, Face* const pF1, Face* const pF2,
		 Face* const pF3, const int iReg = iDefaultRegion);

  virtual TriCell*
  createTriCell(int minlayerindex, Vert * const pV0, Vert * const pV1,
		Vert * const pV2, const int iReg = iDefaultRegion,
		const bool autoRotate = false, int ID = 0);

  virtual TriCell*
  createTriCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
		const int iReg = iDefaultRegion, const bool autoRotate = true,
		int ID = 0);

  QuadCell*
  createQuadCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
		 Vert * const pV3, const int iReg = iDefaultRegion);

  Cell*
  createCell(int minlayerindex, Vert * const vert, Face * const face,
	     const int iReg)
  {
    if (omp_in_parallel() == 0)
      return createTriCell(minlayerindex, vert, face->getVert(0),
			   face->getVert(1), iReg, false, 0);
    else {
      int ID = omp_get_thread_num();
      return createTriCell(minlayerindex, vert, face->getVert(0),
			   face->getVert(1), iReg, false, ID);
    }
  }
  Cell*
  createCell(Vert * const vert, Face * const face, const int iReg)
  {
    if (omp_in_parallel() == 0)
      return createTriCell(vert, face->getVert(0), face->getVert(1), iReg,
			   false);
    else {
      int ID = omp_get_thread_num();
      return createTriCell(vert, face->getVert(0), face->getVert(1), iReg,
			   false, ID);
    }
  }
  bool
  deleteVert(Vert * const pV);
  bool
  deleteFace(Face * const pF);
  bool
  deleteBFace(BFace * const pBF);
  bool
  deleteCell(Cell * const pC);
  bool
  deleteCellWithoutDeletingFaces(Cell * const pC);

  bool
  deleteCellRe(Cell* pC, int ID = 0);
  bool
  deleteFaceRe(Face* pF, int ID = 0);

  //bool deleteCell(Cell *pC);
  ///////////////////////////////////////
  /// Database maintainence functions ///
  ///////////////////////////////////////
  ///
  virtual void
  purgeBdryFaces(std::map<BFace*, BFace*>* bface_map = NULL);
  ///
  virtual void
  purgeFaces(std::map<Face*, Face*>* face_map = NULL);
  ///
  virtual void
  purgeCells(std::map<Cell*, Cell*>* cell_map = NULL);
  ///
  virtual void
  purgeVerts(std::map<Vert*, Vert*>* vert_map = NULL);
  /// Return the number of triangular cells.
  virtual void
  purgeAllEntities(std::map<Vert*, Vert*>* vert_map = NULL,
		   std::map<Face*, Face*>* face_map = NULL,
		   std::map<Cell*, Cell*>* cell_map = NULL,
		   std::map<BFace*, BFace*>* bface_map = NULL);
  void
  resetVertexConnectivity();
  ///
protected:
  virtual void
  reorderVerts_RCM(void);
  virtual void
  reorderCells(void);
  virtual void
  reorderFaces(void);

public:

  ////////////////////////////////////////////////////////
  /// Functions that affect the behavior of algorithms ///
  ////////////////////////////////////////////////////////

  virtual void
  makeSimplicial();

  // Insertion:
  // The worst shape quality measure (shortest edge to circumradius) allowed.
  // The worst shape quality measure allowed.
  virtual double
  getWorstAllowedCellShape() const
  {
    // These should correspond to minimum angles of 20.7, 25.7, and 30 deg,
    // respectively.
    switch (m_encroachType)
      {
      case eBall:
	return sqrt(1. / 6.); // = sin(20.7 deg) / sin(60 deg)
      case eLens:
	return 0.5; // = sin(25.7 deg) / sin(60 deg)
      case eNewLens:
	return sqrt(1. / 3.); // = sin(30 deg) / sin(60 deg)
      default:
	assert(0);
	return 0;
      }
  }
  // Coarsen:
  virtual double
  getCoarseningConstant() const
  {
    return 0.45;
  }
  // End coarsen

  ///////////////////////////////////////////////////////
  /// Low-level algorithms that change mesh geom/topo ///
  ///////////////////////////////////////////////////////
  // Swapping
private:
  virtual bool
  doSwap(const Vert* const pVVertA, const Vert* const pVVertB,
	 const Vert* const pVVertC, const Vert* const pVVertD,
	 const Vert* const pVVertE = pVInvalidVert) const;

public:
  virtual int
  iFaceSwap_deprecated(Face*& pF);

  // Vertex removal
  virtual bool
  removeVert(Vert* const pV, int& iNewSwaps, Vert* apVPreferred[] =
  NULL,
			  const int iNPref = 0);

  bool
  removeVertByRemeshing(Vert* const pV, int& iNewSwaps);
  void
  getNeighborVertsAndFaces(Vert* pV, List<Vert*>& LpVNear,
			   List<Face*>& LpFNear);
  friend class TriMeshBuilder;
  friend class TriMeshRefiner;

  //
  typedef std::pair<BdryEdgeBase*, BdryEdgeBase*> BEdgePair;

  //////////////////////////////////////////////////
  /// High-level algorithms that change the mesh ///
  //////////////////////////////////////////////////
  // This call is deprecated; don't use it in new code.  It's still used
  // in a couple of places, but those are slated for change: coarsening,
  // aniso refinement in 2D, and the STL mesher.
  virtual int
  iSwap_deprecated(const int iMaxPasses = 2, const bool qAlreadyMarked = false);
  virtual bool
  makeDelaunay();

  ////////////////////////////////////////////////////////
  /// Internal routines to assist with mesh coarsening ///
  ////////////////////////////////////////////////////////
public:
  virtual void
  identifyVertexTypesGeometrically() const;
  void
  markCleanNeighbors(Cell* pC);

public:
  //This function is no longer defined in 2D, it is replaced by
  //delete_interior_verts_in_ball (accepting different arguments).
  //Since it is a pure virtual function in the base class, a definition
  //is still needed.
  virtual bool
  deleteInteriorPointsInBall(const Cell*, const double[], const double, int&)
  {
    vFatalError("This function is not defined in 2D",
		"Mesh2D::qDeleteInteriorPointsInBall");
    return false;
  }

public:
  //The following two are called by InsertionQueue to add encroached entities
  //at the top of the priority queue. The new mesh refinement implementation in
  //TriMeshRefiner no longer queues encroached entities. These functions are
  //no longer necessary. The implementation is commented out in Mesh2D.cxx
//  virtual int queueEncroachedBdryEntities(InsertionQueue&,
//                                          const std::set<BFace*>&) const
//  {
//    vFatalError("Function no longer defined in 2D",
//                "Mesh2D::iQueueEncroachedBdryEntities()");
//    return 0;
//  }

//  virtual int queueEncroachedBdryEntities(InsertionQueue&,
//					   const std::set<BFace*>&,
//					   const double[],
//					   const double) const {
//    vFatalError("Function no longer defined in 2D",
//		"Mesh2D::iQueueEncroachedBdryEntities()");
//    return 0;
//  }

  IntBdryEdge*
  getIntBFace(const GR_index_t i) const;

  void
  setupVerts(const GR_index_t i)
  {
    m_ECVerts.reserve(i);
  }
  void
  setupFaces(const GR_index_t i)
  {
    m_ECEdgeF.reserve(i);
  }
  void
  setupBFaces(const GR_index_t i)
  {
    m_ECEdgeBF.reserve(i);
  }
  void
  setupIntBFaces(const GR_index_t i)
  {
    m_ECIntEdgeBF.reserve(i);
  }
  void
  setupTriCells(const GR_index_t i)
  {
    m_ECTri.reserve(i);
  }
  void
  setupQuadCells(const GR_index_t i)
  {
    m_ECQuad.reserve(i);
  }

  /// Add entries to connectivity tables
  Face*
  getNewFace(const int iNV = 2, int ID = 0);
  BFace*
  getNewBFace(const int iNBV = 2, int ID = 0);
  IntBdryEdge*
  getNewIntBdryEdge(const int iNBV = 2, int ID = 0);
  Cell*
  getNewCell(const int iNF = 3, int ID = 0);
  /// Prints a summary of the patches...
  void
  checkPatches();
  // Mesh optimization, point insertion, etc
  int
  reconfigure(Face*& pF);
  ///
  bool
  recoverEdge(const int iV0, const int iV1);
  bool
  recoverEdge(Vert* const pV0, Vert* const pV1);
  bool
  recoverEdge(Vert* const pV0, Vert* const pV1, Face*& pF);
  ///
  bool
  removeEdge(Face*& pF);
  bool
  removeVertCheaply(Vert* const pV, int& iNewSwaps);

  Cell*
  beginPipe(const Vert* const, const Vert* const) const
  {
    assert(0);
    return 0;
  }
public:

  /// PreProcessing for Ruppert - 1
  void
  checkForLargeAngles(const double dMaxAngle = M_PI_2); // 90 deg
  /// PreProcessing for Ruppert - 2
  void
  checkForSmallAngles();
  //
  // The following three functions are used in the experimental
  // anisotropic meshing code.
  void
  generateAnisotropic(const int iNStructBC, const int aiStructBC[],
		      const double dScale, const double dGrading);
  void
  refineBoundaryToLengthScale(const int iNStructBC, const int aiStructBC[]);

  // both of these assume the cells exist
  bool
  createQuadFromTris(Cell* const pCA, Cell* const pCB);
  bool
  createQuadFromVerts(Vert* const pVA, Vert* const pVB, Vert* const pVC,
		      Vert* const pVD);
  // the opposite, so I can unit test more easily
  bool
  createTrisFromQuad(Cell * const pC);
  void
  setLSChecking(bool qCheck);
  void
  TriMergetoQuad();
  void
  TriMergetoQuad2();
  double
  calcsquareness(Vert* const v0, Vert* const v1, Vert* const v2,
		 Vert* const v3);
  bool
  updateQueue(InsertionQueue MergeQ, Cell* const pCA, Cell* const pCB);
  bool
  Orientation(TriCell* Tri1, TriCell* Tri2);

protected:
  void
  getBdry2DInfo(const Bdry2D& B2DInfo);

  ///
public:
  void
  writeBdryPatchFile(const char strBaseFileName[],
		     const char strExtraFileSuffix[]) const;

public:
  /// Iterator functions
  EntContainer<EdgeFace>::iterator
  edgeFace_begin() const
  {
    return m_ECEdgeF.begin();
  }
  EntContainer<EdgeFace>::iterator
  edgeFace_end() const
  {
    return m_ECEdgeF.end();
  }
  EntContainer<EdgeFace>::iterator
  edgeFace_advance(EntContainer<EdgeFace>::iterator& iter, GR_index_t num)
  {
    return m_ECEdgeF.advance(iter, num);
  }

  EntContainer<Vert>::iterator
  Vert_begin() const
  {
    return m_ECVerts.begin();
  }
  EntContainer<Vert>::iterator
  Vert_end() const
  {
    return m_ECVerts.end();
  }
  EntContainer<Vert>::iterator
  Vert_advance(EntContainer<Vert>::iterator& iter, GR_index_t num)
  {
    return m_ECVerts.advance(iter, num);
  }

  EntContainer<TriCell>::iterator
  triCell_begin() const
  {
    return m_ECTri.begin();
  }
  EntContainer<TriCell>::iterator
  triCell_end() const
  {
    return m_ECTri.end();
  }

  friend class AnisoRefinement;

  Cell*
  findCell(const double newPtLoc[], Cell* const cellGuess, int &numViolated,
	   Face *facesViolated[3], int &numTies, Face *facesTied[3]) const;

  virtual void
  writeTempMesh();

////To perturbe a mesh by a percentaje of its edge length. Used for testing only
//  void PerturbeMesh(double dPertMeasure);
//
////Adds a line of quad cells outside of the boundary of the mesh. Used for testing only
//  void QuadOutsideBdry(double dLength);

////Turns a quad cell into two triangular cells
//  void QuadintoTriMesh(Cell *pC);
};

/// Write a file in GRUMMP native format. Picks whether to write with or
/// without internal bdry faces.
void
writeNative(Mesh2D& OutMesh, const char strBaseFileName[],
	    const bool qHighOrder = false,
	    const char strExtraFileSuffix[] = "");
/// Write a file in GRUMMP native format without interior bdry faces.
void
writeNativeMESH(Mesh2D& OutMesh, const char strBaseFileName[],
		const bool qHighOrder = false, const char strExtraFileSuffix[] =
		    "");
/// Write a file in GRUMMP native format with interior bdry faces.
void
writeNativeGRMESH(Mesh2D& OutMesh, const char strBaseFileName[],
		  const bool qHighOrder = false,
		  const char strExtraFileSuffix[] = "");
/// Write a file in generic finite element (cell-vertex) format
void
writeFEA(Mesh2D& OutMesh, const char strBaseFileName[], const int iOrder = 2,
	 const char strExtraFileSuffix[] = "");
/// Write a file in the VTK legacy format (.vtk, uncompressed)
void
writeVTK(Mesh2D& OutMesh, const char strBaseFileName[],
	       const char strExtraFileSuffix[] = "",
	       const enum Mesh::ExtraData ED = Mesh::eStandard);

void
writeVTKMetric(Mesh2D& OutMesh, const char strBaseFileName[],
	       const char strExtraFileSuffix[] = "",
	       const int *cellData = NULL);

void
writeTriangle(Mesh2D& OutMesh, const char strBaseFileName[]);

/// Read a file in GRUMMP native format (has Internal Boundaries)
void
readNativeGRMESH(Mesh2D& OutMesh, const char strFileName[]);

/// Read a file in Triangle format
void
readTriangle(Mesh2D& OutMesh, const char strFileName[]);

/// Read in a file without Internal Boundaries
void
readNativeMESH(Mesh2D& OutMesh, const char strFileName[]);

/// Read a file in generic finite element (cell-vertex) format
void
readFEA(Mesh2D& OutMesh, const char strFileName[]);

/// Read a file in the VTK legacy format (.vtk, uncompressed)
void
readVTK(Mesh2D& OutMesh, const char strFileName[]);

// TODO: This function should eventually live outside the core database, with either
// mesh builders or inserters.
void
createMeshInBox(Mesh2D& Mesh, double dXLo, double dYLo, double dXHi,
		double dYHi);

#endif
