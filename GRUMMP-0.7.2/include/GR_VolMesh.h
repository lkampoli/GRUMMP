#ifndef GR_VolMesh
#define GR_VolMesh 1

#include "GR_config.h"
#include "GR_BFace3D.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_List.h"
#include "GR_Mesh.h"
#include "GR_EntContainer.h"
#include "GR_Subseg.h"

class CubitBox;
class TetRefinerEdge;

///
class VolMesh : public Mesh {
private:
  ///
  bool m_OKToSwapSurfaceEdges :1, m_doEdgeSwaps :1, m_strictPatchChecking :1;
  ///
  // This angle is used to determine when swapping of surface edges is
  // allowed, and also to decide when a surface edge is actually on a
  // bdry curve or not (in the absence of other info and of surface
  // patch recovery code)
  double m_maxAngleForSurfSwap; // In degrees
  /// Three-dimensional topology
  EntContainer<TriFace> m_ECTriF;      /// Faces
  EntContainer<QuadFace> m_ECQuadF;     /// Faces
  EntContainer<TriBFace> m_ECTriBF;     /// BFaces
  EntContainer<QuadBFace> m_ECQuadBF;    /// BFaces
  EntContainer<TetCell> m_ECTet;       /// Cells
  EntContainer<PyrCell> m_ECPyr;       /// Cells
  EntContainer<PrismCell> m_ECPrism;     /// Cells
  EntContainer<HexCell> m_ECHex;       /// Cells
  EntContainer<FlakeCell> m_ECFlake;     /// Cells
  EntContainer<IntTriBFace> m_ECIntTriBF;  /// Internal bdry faces
  EntContainer<IntQuadBFace> m_ECIntQuadBF; /// Internal bdry faces
  Bdry3D *m_Bdry3D;       /// Bdry info
  bool m_bdryFromMesh;
  SubsegMap *m_SubsegMap; // contains subseg info
private:
  /// Operator= is disabled, because I can't see a use case for it.
  VolMesh&
  operator=(const VolMesh&)
  {
    assert(0);
    return *this;
  }
  /// No use case for this, either.
  VolMesh(const VolMesh&);
public:

  double
  Mesh3D_sync_entry(double& paralleltime, bool IsSwapping)
  {
    double Sync_Entry_StartTime = omp_get_wtime();
    if (NUM_PROCS > 1) {
      if (IsSwapping) {
	paralleltime += m_ECTriF.sync_entry(IsSwapping);
	paralleltime += m_ECTriBF.sync_entry(IsSwapping);
	paralleltime += m_ECTet.sync_entry(IsSwapping);
      }
      else {
	paralleltime += m_ECVerts.sync_entry();
	paralleltime += m_ECTriF.sync_entry();
	paralleltime += m_ECTriBF.sync_entry();
	paralleltime += m_ECTet.sync_entry();
//			  paralleltime +=m_ECQuadF.sync_entry();
//			  paralleltime +=m_ECQuadBF.sync_entry();
//			  paralleltime +=m_ECPyr.sync_entry();
//			  paralleltime +=m_ECPrism.sync_entry();
//			  paralleltime +=m_ECHex.sync_entry();
//	paralleltime += m_ECFlake.sync_entry();
//			  paralleltime +=m_ECIntTriBF.sync_entry();
//			  paralleltime +=m_ECIntQuadBF.sync_entry();
      }
    }
    Sync_Entry_StartTime = omp_get_wtime() - Sync_Entry_StartTime;
    return Sync_Entry_StartTime;
  }

  double
  Mesh3D_sync_Exit(bool IsSwapping)
  {
    double Sync_Exit_StartTime = omp_get_wtime();
    if (NUM_PROCS > 1) {
      if (IsSwapping) {
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TriFace*>::iterator iterS = m_ECTriF.empties_pa[ID].begin();
	  std::list<TriFace*>::iterator iterE = m_ECTriF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTriF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TriBFace*>::iterator iterS =
	      m_ECTriBF.empties_pa[ID].begin();
	  std::list<TriBFace*>::iterator iterE = m_ECTriBF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTriBF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TetCell*>::iterator iterS = m_ECTet.empties_pa[ID].begin();
	  std::list<TetCell*>::iterator iterE = m_ECTet.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTet.sync_exit();
      }
      else {
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  std::list<Vert*>::iterator iterS = m_ECVerts.empties_pa[ID].begin();
	  std::list<Vert*>::iterator iterE = m_ECVerts.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECVerts.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TriFace*>::iterator iterS = m_ECTriF.empties_pa[ID].begin();
	  std::list<TriFace*>::iterator iterE = m_ECTriF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTriF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TriBFace*>::iterator iterS =
	      m_ECTriBF.empties_pa[ID].begin();
	  std::list<TriBFace*>::iterator iterE = m_ECTriBF.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTriBF.sync_exit();
	for (int ID = 0; ID < NUM_PROCS; ++ID) {
	  std::list<TetCell*>::iterator iterS = m_ECTet.empties_pa[ID].begin();
	  std::list<TetCell*>::iterator iterE = m_ECTet.empties_pa[ID].end();
	  for (; iterS != iterE; ++iterS)
	    (*iterS)->markAsDeleted();
	}
	m_ECTet.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<QuadFace*>::iterator iterS=m_ECQuadF.empties_pa[ID].begin();
//				  std::list<QuadFace*>::iterator iterE=m_ECQuadF.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECQuadF.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<QuadBFace*>::iterator iterS=m_ECQuadBF.empties_pa[ID].begin();
//				  std::list<QuadBFace*>::iterator iterE=m_ECQuadBF.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECQuadBF.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<PyrCell*>::iterator iterS=m_ECPyr.empties_pa[ID].begin();
//				  std::list<PyrCell*>::iterator iterE=m_ECPyr.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECPyr.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<PrismCell*>::iterator iterS=m_ECPrism.empties_pa[ID].begin();
//				  std::list<PrismCell*>::iterator iterE=m_ECPrism.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECPrism.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<HexCell*>::iterator iterS=m_ECHex.empties_pa[ID].begin();
//				  std::list<HexCell*>::iterator iterE=m_ECHex.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECHex.sync_exit();
	//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
	//			  {
	//				  std::list<FlakeCell*>::iterator iterS=m_ECFlake.empties_pa[ID].begin();
	//				  std::list<FlakeCell*>::iterator iterE=m_ECFlake.empties_pa[ID].end();
	//				  for ( ; iterS != iterE ; ++iterS)
	//					  (*iterS)->markAsDeleted();
	//			  }
	//			  m_ECFlake.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<IntTriBFace*>::iterator iterS=m_ECIntTriBF.empties_pa[ID].begin();
//				  std::list<IntTriBFace*>::iterator iterE=m_ECIntTriBF.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECIntTriBF.sync_exit();
//			  for (int ID=0 ; ID<NUM_PROCS ; ++ID)
//			  {
//				  std::list<IntQuadBFace*>::iterator iterS=m_ECIntQuadBF.empties_pa[ID].begin();
//				  std::list<IntQuadBFace*>::iterator iterE=m_ECIntQuadBF.empties_pa[ID].end();
//				  for ( ; iterS != iterE ; ++iterS)
//					  (*iterS)->markAsDeleted();
//			  }
//			  m_ECIntQuadBF.sync_exit();
      }
    }
    Sync_Exit_StartTime = omp_get_wtime() - Sync_Exit_StartTime;
    return Sync_Exit_StartTime;
  }

  ///
//  VolMesh(const char strBaseFileName[], const int iQualMeas = 2,
//	  const double dMaxAngle = 0, Bdry3D * const pB3DIn = NULL);
  ///Added March 2007 by SG.
  VolMesh(const CubitBox& box, const int iQualMeas = 2, const double dMaxAngle =
	      0);
  ///Modified March 2007 by SG.
  VolMesh(const EntContainer<Vert>& EC, const int iQualMeas = 2,
	  const double dMaxAngle = 0, const CubitBox* const = NULL);
  ///
  VolMesh(SurfMesh &SM, const int iQualMeas = 2, const double dRes = 1,
	  const double dGrade = 1, const double dMaxAngle = 0);
  ///
  VolMesh(const int iQualMeas = 2, const double dMaxAngle = 0);
  ///
  virtual
  ~VolMesh();
  ///
  void
  readFromFile(const char strBaseFileName[], int* numVars = nullptr,
	       double** data = nullptr);
  ///

  void
  setBdryData(Bdry3D* const bdryData);
public:
  void
  createPatchesFromBdryFaces(const int aiBCList[]);
public:
  void
  initMeshInBrick(const double dXMin, const double dXMax, const double dYMin,
		  const double dYMax, const double dZMin, const double dZMax);
  ///
  void
  convertFromCellVert(const GR_index_t iNCells, const GR_index_t iNBdryFaces,
		      const GR_index_t a2iCellVert[][4],
		      const GR_index_t aiRegions[],
		      const GR_index_t a2iBFaceVert[][3]);
  ///
  void
  convertFromCellVert(const GR_index_t nTet, const GR_index_t nPyr,
		      const GR_index_t nPrism, const GR_index_t nHex,
		      const GR_index_t nBdryTri, const GR_index_t nBdryQuad,
		      const GR_index_t tetToVert[][4],
		      const GR_index_t pyrToVert[][5],
		      const GR_index_t prismToVert[][6],
		      const GR_index_t hexToVert[][8],
		      const GR_index_t bdryTriToVert[][3],
		      const GR_index_t bdryQuadToVert[][4]);
  ///
  eMeshType
  getType() const
  {
    return Mesh::eVolMesh;
  }
  ///
  GR_index_t
  getNumFaces() const
  {
    return m_ECTriF.lastEntry() + m_ECQuadF.lastEntry();
  }

  GR_index_t
  getNumEdgesInUse() const
  {
    return 0;
  }
  ///
  GR_index_t
  getNumTriFaces() const
  {
    return m_ECTriF.lastEntry();
  }
  GR_index_t
  getNumTrisInUse() const
  {
    return m_ECTriF.size();
  }
  ///
  GR_index_t
  getNumQuadFaces() const
  {
    return m_ECQuadF.lastEntry();
  }
  GR_index_t
  getNumQuadsInUse() const
  {
    return m_ECQuadF.size();
  }
  ///
  GR_index_t
  getNumBdryFaces() const
  {
    return m_ECTriBF.lastEntry() + m_ECQuadBF.lastEntry();
  }
  GR_index_t
  getNumBdryFacesInUse() const
  {
    return m_ECTriBF.size() + m_ECQuadBF.size();
  }
  ///
  GR_index_t
  getNumTriBdryFaces() const
  {
    return m_ECTriBF.lastEntry();
  }
  ///
  GR_index_t
  getNumQuadBdryFaces() const
  {
    return m_ECQuadBF.lastEntry();
  }
  ///
  GR_index_t
  getNumIntBdryFaces() const
  {
    return m_ECIntTriBF.lastEntry() + m_ECIntQuadBF.lastEntry();
  }
  ///
  GR_index_t
  getNumIntTriBdryFaces() const
  {
    return m_ECIntTriBF.lastEntry();
  }
  ///
  GR_index_t
  getNumIntQuadBdryFaces() const
  {
    return m_ECIntQuadBF.lastEntry();
  }
  ///
  GR_index_t
  getNumCells() const
  {
    return m_ECTet.lastEntry() + m_ECPyr.lastEntry() + m_ECPrism.lastEntry()
	+ m_ECHex.lastEntry() + m_ECFlake.lastEntry();
  }
  ///
  GR_index_t
  getNumTetCells() const
  {
    return m_ECTet.lastEntry();
  }
  GR_index_t
  getNumTetsInUse() const
  {
    return m_ECTet.size();
  }
  ///
  GR_index_t
  getNumPyrCells() const
  {
    return m_ECPyr.lastEntry();
  }
  GR_index_t
  getNumPyrsInUse() const
  {
    return m_ECPyr.size();
  }
  ///
  GR_index_t
  getNumPrismCells() const
  {
    return m_ECPrism.lastEntry();
  }
  GR_index_t
  getNumPrismsInUse() const
  {
    return m_ECPrism.size();
  }
  ///
  GR_index_t
  getNumHexCells() const
  {
    return m_ECHex.lastEntry();
  }
  GR_index_t
  getNumHexesInUse() const
  {
    return m_ECHex.size();
  }
  ///
  GR_index_t
  getNumFlakeCells() const
  {
    return m_ECFlake.lastEntry();
  }
  GR_index_t
  getNumFlakesInUse() const
  {
    return m_ECFlake.size();
  }
  ///
  Face*
  getFace(const GR_index_t i) const;
  ///
  BFace*
  getBFace(const GR_index_t i) const;
  ///
  BFace*
  getIntBFace(const GR_index_t i) const;
  ///
  Cell*
  getCell(const GR_index_t i) const;
  /// Index retrieval
  virtual GR_index_t
  getFaceIndex(const Face* const pF) const;
  ///
  virtual GR_index_t
  getBFaceIndex(const BFace* const pBF) const;
  ///
  virtual GR_index_t
  getCellIndex(const Cell* const pC) const;
  /// TODO: These are here only for a clunky way of reading .vmesh files
  void
  reserveTriFaces(GR_index_t size)
  {
    m_ECTriF.reserve(size);
  }
  void
  reserveTets(GR_index_t size)
  {
    m_ECTet.reserve(size);
  }
  void
  reserveTriBFaces(GR_index_t size)
  {
    m_ECTriBF.reserve(size);
  }
#ifdef ITAPS
  ///
  virtual bool
  isValidEntHandle(void* pvEnt) const;
#endif
  /// Add entries to connectivity tables
  Face*
  getNewFace(const int iNV = 3, int ID = 0);
  ///
  BFace*
  getNewBFace(const int iNBV = 3, const bool qIsInternal = false, int ID = 0);
  ///
  Cell*
  getNewCell(const int iNF = 4, const int iNV = 4, int ID = 0);
  ///
  Cell*
  getNewCell(const Cell* const pC); // Returned cell has
  // same type and
  // data as old one.
  /// New topology modifiers, a la ITAPS
  virtual Face*
  createFace(Vert * const /*apV*/[], const int /*iNVerts*/)
  {
    return pFInvalidFace;
  }
  /// Note the changes made by Dan, to have it identify IntBdryFaces vs
  /// BdryFaces by number of attached cells.
  virtual BFace*
  createBFace(Face * const pF, BFace * const pBF =
  pBFInvalidBFace);
  BFace*
  createBFaceParallel(Face * const pF, BFace * const pBF =
  pBFInvalidBFace,
		      int ID = 0);
  /*   virtual Cell* createCell(Face * const /\* apF *\/[], const int /\* iNFaces *\/, */
  /* 			   const int /\* iReg = iDefaultRegion *\/) */
  /*   {return pCInvalidCell;} */
  /*   virtual Cell* createCell(Vert * const /\* apV *\/[], const int /\* iNVerts *\/, */
  /* 			   const int /\* iReg = iDefaultRegion *\/) */
  /*   {return pCInvalidCell;} */

  Vert*
  createVert(const double dX, const double dY, const double dZ);
  Vert*
  createVert(const double adCoords[]);

  /// Create a new vertex and connect faces to the new one
  Vert*
  createVert(const Vert* pVIn);
  /// Create a new vertex without face connectivity
  Vert*
  duplicateVert(const Vert* pVIn);

  Face*
  createFace(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
	     Vert * const pVC, Vert * const pVD = pVInvalidVert,
	     bool qForceDuplicate = false, int ID = 0);
  TetCell*
  createTetCell(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
		Vert * const pVC, Vert * const pVD, const int iRegion =
		iDefaultRegion,
		bool OMP_Flag = false, int ID = 0);
  /// The following version should only be used if the correct vertex
  /// ordering isn't known in advance.
  TetCell*
  createTetCell(bool &qExistedAlready, Vert * const pV, Face * const pF,
		const int iRegion = iDefaultRegion);
  TetCell*
  createTetCell(bool &qExistedAlready, Face * const pFA, Face * const pFB,
		Face * const pFC, Face * const pFD, const int iRegion =
		iDefaultRegion,
		int ID = 0);
  Cell*
  createCell(Vert * const vert, Face * const face, const int region)
  {
    assert(face->getType() == Face::eTriFace);
    bool temp;
    return createTetCell(temp, vert, face, region);
  }
  Cell*
  createCell(int minlayerindex, Vert * const vert, Face * const face,
	     const int region)
  {
    (void) minlayerindex;
    assert(face->getType() == Face::eTriFace);
    bool temp;
    return createTetCell(temp, vert, face, region);
  }

  // Vertex order: cyclic around the quad.
  FlakeCell*
  createFlakeCell(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
		  Vert * const pVC, Vert * const pVD, const int iRegion =
		  iDefaultRegion);
  // Vertex order: cyclic around bottom, then top.  RH normal points in
  // for bottom.
  PyrCell*
  createPyrCell(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
		Vert * const pVC, Vert * const pVD, Vert * const pVE,
		const int iRegion = iDefaultRegion);
  // Faces are: tris in cyclic order, then quad.  Verts, if given, are
  // the verts in the same order as above.  See implementation for a
  // sketch.
  PyrCell*
  createPyrCell(bool &qExistedAlready, Face* const pF014, Face* const pF124,
		Face* const pF234, Face* const pF304, Face* const pF0123,
		const int iRegion = iDefaultRegion, Vert * pVA =
		pVInvalidVert,
		Vert * pVB = pVInvalidVert, Vert * pVC =
		pVInvalidVert,
		Vert * pVD = pVInvalidVert, Vert * pVE =
		pVInvalidVert);
  // Vertex order: cyclic around bottom, then top.  RH normal points in
  // for bottom, out for top.
  PrismCell*
  createPrismCell(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
		  Vert * const pVC, Vert * const pVD, Vert * const pVE,
		  Vert * const pVF, const int iRegion =
		  iDefaultRegion);
  // Faces are: quads in cyclic order, bottom tri, top tri. Verts, if
  // given, are the verts in the same order as above.
  PrismCell*
  createPrismCell(bool &qExistedAlready, Face * const pFA, Face * const pFB,
		  Face * const pFC, Face * const pFD, Face * const pFE,
		  const int iRegion = iDefaultRegion, Vert * pVA =
		  pVInvalidVert,
		  Vert * pVB = pVInvalidVert, Vert * pVC =
		  pVInvalidVert,
		  Vert * pVD = pVInvalidVert, Vert * pVE =
		  pVInvalidVert,
		  Vert * pVF = pVInvalidVert);
protected:
  // Vertex order: cyclic at one end, then cyclic at the other (same
  // orientation, same starting corner.  ABC is the bottom, DEF is the top.
  PrismCell*
  replaceCellsWithPrism(Vert* pVA, Vert* pVB, Vert* pVC, Vert* pVD, Vert* pVE,
			Vert* pVF);
public:
  // Vertex order: cyclic around bottom, then top.  RH normal points in
  // for bottom, out for top.
  HexCell*
  createHexCell(bool &qExistedAlready, Vert * const pVA, Vert * const pVB,
		Vert * const pVC, Vert * const pVD, Vert * const pVE,
		Vert * const pVF, Vert * const pVG, Vert * const pVH,
		const int iRegion = iDefaultRegion);
  // Faces are: sides in cyclic order, bottom, top. Verts, if
  // given, are the verts in the same order as above.
  HexCell*
  createHexCell(bool &qExistedAlready, Face* const pFA, Face* const pFB,
		Face* const pFC, Face* const pFD, Face* const pFE,
		Face* const pFF, const int iRegion = iDefaultRegion,
		Vert * pVA =
		pVInvalidVert,
		Vert * pVB = pVInvalidVert, Vert * pVC =
		pVInvalidVert,
		Vert * pVD = pVInvalidVert, Vert * pVE =
		pVInvalidVert,
		Vert * pVF = pVInvalidVert, Vert * pVG =
		pVInvalidVert,
		Vert * pVH = pVInvalidVert);

  /*   Cell* createCell(bool &qExistedAlready, */
  /* 		   Face * const pFA, Face * const pFB, Face * const pFC, */
  /* 		   Face * const pFD, Face * const pFE, */
  /* 		   const int iRegion = iDefaultRegion); */

  bool
  deleteVert(Vert * const pV);
  bool
  deleteFace(Face * const pF);
  bool
  deleteBFace(BFace * const pBF);
  bool
  deleteCell(Cell * const pC);

  bool
  deleteFaceParallel(Face * const pF, int ID = 0);
  bool
  deleteBFaceParallel(BFace * const pBF, int ID = 0);
  bool
  deleteCellParallel(Cell * const pC, int ID = 0);

  /// Clean up connectivity tables
  void
  purgeCells(std::map<Cell*, Cell*>* cell_map = NULL);
  ///
  void
  purgeBdryFaces(std::map<BFace*, BFace*>* bface_map = NULL);
  ///
  void
  purgeFaces(std::map<Face*, Face*>* face_map = NULL);
  ///
  void
  purgeVerts(std::map<Vert*, Vert*>* vert_map = NULL);
  ///
  void
  purgeAllEntities(std::map<Vert*, Vert*>* vert_map = NULL,
		   std::map<Face*, Face*>* face_map = NULL,
		   std::map<Cell*, Cell*>* cell_map = NULL,
		   std::map<BFace*, BFace*>* bface_map = NULL);

  /*   /// Ensure that the correct number of verts exist */
  /*   void vSetupVerts(const int i, */
  /* 		   const bool qExact = false) {ECVerts.vSetup(i, qExact);} */
  ///
  void
  assignInitialLengthScales(double dRes, double dAlpha);
//protected:
  virtual void
  reorderVerts_RCM(void);
  virtual void
  reorderCells(void);
  virtual void
  reorderFaces(void);
public:
  double
  calcBoundarySize() const;
  double
  calcInteriorSize() const;
  bool
  isWatertight() const;
  bool
  isMeshOrientationOK(const bool qFix = true);
/// Conversion to/from simplicial form
  ///
  void
  makeSimplicial()
  {
    assert(getNumCells() == getNumTetCells());
    m_simplexMesh = true;
  }

public:

  /// Validity checking
  bool
  isValid() const;

//  // Bdry encroachment check.
//  virtual int queueEncroachedBdryEntities(InsertionQueue& IQ,
//					   const std::set<BFace*>& spBF) const;
//  virtual int queueEncroachedBdryEntities(InsertionQueue& IQ,
//					   const std::set<BFace*>& spBF,
//					   const double adLoc[],
//					   const double dOffset = 0) const;
public:
  // The worst shape quality measure allowed.
  virtual double
  getWorstAllowedCellShape() const
  {
    switch (m_encroachType)
      {
      case eBall:
	return sqrt(0.09375); // sqrt(3) sqrt(2) / 8
      case eLens:
	return sqrt(0.1875);  // sqrt(3) / 4
      case eNewLens:
	return sqrt(0.375);   // sqrt(3) sqrt(2) / 4
	// The mythical best-case scenario would be to use sqrt(3) / 2,
	// because then slivers would be impossible.
      default:
	assert(0);
	return 1;
      }
  }
// Mesh optimization routines
/// Sliver fix
public:
  int
  repairBadCells(bool qAllowInsertion = false);
  ///
  int
  breakBadBdryCells(bool qAllowInsertion = false);
  ///
  int
  splitLargeAngle(Cell *pC, const double dSplitAngle,
		  const double * const adDihed = 0, const bool qOnlyInterior =
		      false);

  ///
private:
  virtual bool
  qDoSwap(const Vert* const pVVertA, const Vert* const pVVertB,
	  const Vert* const pVVertC, const Vert* const pVVertD,
	  const Vert* const pVVertE) const;

public:
  // This call is deprecated; don't use it in new code.  It's still used
  // in a couple of places, but those are slated for change: coarsening,
  // aniso refinement in 2D, and the STL mesher.
  virtual int
  iSwap_deprecated(const int iMaxPasses = 2, const bool qAlreadyMarked = false);
  ///
  int
  iFaceSwap_deprecated(Face*& pF);
  ///
  bool
  areSurfaceEdgeChangesAllowed() const
  {
    return m_OKToSwapSurfaceEdges;
  }
  ///
  void
  allowSurfaceEdgeChanges(const double dMaxAngle = 0)
  {
    assert(dMaxAngle >= 0 && dMaxAngle < 180.);
    m_OKToSwapSurfaceEdges = true;
    m_maxAngleForSurfSwap = dMaxAngle;
  }
  ///
  double
  getMaxAngleForSurfaceEdgeChanges() const
  {
    return m_maxAngleForSurfSwap;
  }
  ///
  void
  disallowSurfaceEdgeChanges()
  {
    m_OKToSwapSurfaceEdges = false;
  }
  ///
  bool
  isStrictPatchCheckingUsed() const
  {
    return m_strictPatchChecking;
  }
  ///
  void
  useStrictPatchChecking()
  {
    m_strictPatchChecking = true;
  }

  bool
  isBdryFromMesh() const
  {
    return m_bdryFromMesh;
  }
  ///
  void
  useLenientPatchChecking()
  {
    m_strictPatchChecking = false;
  }
  virtual GR_index_t
  getNumBdryPatches()
  {
    return m_Bdry3D->iNumPatches();
  }
  ///
  bool
  areBdryFacesOKToSwap(const TriBFaceBase * const pBF0,
		       const TriBFaceBase * const pBF1) const;
  ///
  bool
  isEdgeSwappingAllowed() const
  {
    return m_doEdgeSwaps;
  }
  ///
  void
  allowEdgeSwapping()
  {
    m_doEdgeSwaps = true;
  }
  ///
  void
  disallowEdgeSwapping()
  {
    m_doEdgeSwaps = false;
  }

  ///
  void
  refineToLengthScale(const bool qSameSize, const int iMaxPasses = 40);
  ///
  virtual bool
  removeVert(Vert * const pV, int& iNewSwaps,
			  Vert * apVPreferred[] = NULL, const int iNPref = 0);
  virtual void
  identifyVertexTypesGeometrically() const;
  virtual double
  getCoarseningConstant() const
  {
    return 0.36;
  }
  ///
  bool
  makeDelaunay();
  ///
  int
  reconfigure_deprecated(TetCell* apTCTets[], TriFace* const pFFace,
			 Vert* const pVVertA, Vert* const pVVertB,
			 Vert* const pVVertC, Vert* const pVVertD,
			 Vert* const pVVertE);
  ///
  int
  reconfTet23_deprecated(TetCell* apTCTets[], TriFace* const pFFace,
			 Vert* const pVVertA, Vert* const pVVertB,
			 Vert* const pVVertC, Vert* const pVVertD,
			 Vert* const pVVertE);
  ///
  int
  reconfTet32_deprecated(TetCell* apTCTets[], TriFace* const pFFace,
			 Vert* const pVVertA, Vert* const pVVertB,
			 Vert* const pVVertC, Vert* const pVVertD,
			 Vert* const pVVertE);
  ///
  int
  reconfTet22_deprecated(TetCell* apTCTets[], TriFace* const pFFace,
			 Vert* const pVVertA, Vert* const pVVertB,
			 Vert* const pVVertC, Vert* const pVVertD,
			 Vert* const pVVertE);
  ///
  int
  reconfTet44_deprecated(TetCell* apTCTets[], TriFace* const pFFace,
			 Vert* const pVVertA, Vert* const pVVertB,
			 Vert* const pVVertC, Vert* const pVVertD,
			 Vert* const pVVertE);
  ///
  int
  edgeSwap3D(Face* const pF, Vert* pVNorth, Vert* pVSouth, Vert* const pVOther,
	     const bool qAllowAnyValid = false, const Vert* const pVWant0 =
	     pVInvalidVert,
	     const Vert* const pVWant1 = pVInvalidVert);
  ///
  int
  bdryEdgeSwap3D(Face* const pF, Vert* pVNorth, Vert* pVSouth,
		 Vert* const pVOther);
private:
  double
  evalQualFunc(Vert* const pV0, Vert* const pV1, Vert* const pV2,
	       Vert* const pV3) const;

  GR_index_t
  purgeCellRange(GR_index_t iFor, GR_index_t iBack,
		 std::map<Cell*, Cell*>* cell_map);
  void
  handleFlakeCell(Cell* pFCABED, Vert* pVA, Vert* pVB, Vert* pVE, Vert* pVD);

public:
  // Boundary recovery routines
  void
  recoverConstrainedSurface(SurfMesh& SM, const int iNumOrigVerts,
			    const int iNumBBoxVerts = 8);
  void
  markCleanNeighbors(Cell* pC);
  int
  recoverFace(const int iVerts[3], int * const piPointsInserted = NULL,
	      const bool qAllowEdgeSwap = true);
  int
  recoverFace(const Vert * pV0, const Vert * pV1, const Vert * pV2,
	      int * const piPointsInserted = NULL, const bool qAllowEdgeSwap =
		  true,
	      const int iRecurDepth = 0);
  int
  recoverFaceByRemoval(const Vert * pV0, const Vert * pV1, const Vert * pV2,
		       const int iRecurDepth = 0);
  int
  recoverEdge(const Vert* const pV0, const Vert* const pV1,
	      const bool qInsertionOK = false);
  int
  recoverEdge(const Vert* const pV0, const Vert* const pV1, SurfMesh& SM,
	      const bool qInsertionOK = false);
  void
  splitEdge(const Vert* const pV0, const Vert* const pV1, SurfMesh& SM);
  // Find the cell incident on pV0 that covers the ray towards pV1 (may be
  // a multiple choice question if pV0 and pV1 are adjacent, or if there
  // is a face containing one but not the other and coplanar with both.
  Cell*
  beginPipe(const Vert* const pV0, const Vert* const pV1) const;
  int
  clearPipe(const Vert* const pVInit, const Vert* const pVEnd, bool& qHitEdge,
	    Face* apFSkipped[], int& iNSkip);
  void
  deleteExternal(const SurfMesh* pSM, const int iNumOrigVerts,
		 const int iNumBBoxVerts = 8);
  void
  identifyBadEdgesAndFaces(const Vert* const pV0, const Vert* const pV1,
			   const Vert* const pV2,
			   std::set<BadEdge>& sBEBadEdges,
			   std::set<Face*>& spFBadFaces);
public:
  /// Mesh quality assessment
  void
  classifyTets() const;
  ///
  friend class VolQual;

public:
  /// Iterator functions
  EntContainer<TriFace>::iterator
  triFace_begin() const
  {
    return m_ECTriF.begin();
  }
  EntContainer<TriFace>::iterator
  triFace_end() const
  {
    return m_ECTriF.end();
  }
  EntContainer<TriFace>::iterator
  triFace_advance(EntContainer<TriFace>::iterator& iter, GR_index_t num)
  {
    return m_ECTriF.advance(iter, num);
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

  Cell*
  findCell(const double newPtLoc[], Cell* const cellGuess, bool& status) const;

  virtual void
  writeTempMesh();

  // Subseg handling functions go here
  // For now just call subsegMap function.
  // everything here just calls the map directly

  // checks whether the map has been initialized, used in one place for now
  bool
  hasSubsegs()
  {
    return (m_SubsegMap != NULL);
  }

  // called in RefinementManager3D for now
  bool
  isSubsegMapEmpty() const
  {
    return m_SubsegMap->empty();
  }

  bool
  isSubsegVert(Vert * vert) const
  {
    return m_SubsegMap->isSubsegVert(vert);
  }
  // used for debugging purposes, mostly
  int
  subsegMapSize() const
  {
    return m_SubsegMap->size();
  }

  void
  addSubseg(Subseg* const subseg)
  {
    m_SubsegMap->add_subseg(subseg);
  }
  void
  removeSubseg(Subseg* const subseg)
  {
    m_SubsegMap->remove_subseg(subseg);
  }

  Subseg *
  areSubsegVertsConnected(Vert* const vert0, Vert* const vert1) const
  {
    if (vert0->getVertType() == Vert::eInterior
	|| vert1->getVertType() == Vert::eInterior)
      return NULL;
    return m_SubsegMap->verts_are_connected(vert0, vert1);
  }

  // useful when reconnecting subsegs, will be refactored elsewhere
  std::pair<std::set<Subseg*>::const_iterator, std::set<Subseg*>::const_iterator>
  getSubsegRange(Vert* const vert0) const
  {
    return m_SubsegMap->get_subseg_range(vert0);
  }

  void
  getAllSubsegs(std::set<Subseg*>& all_subsegs) const
  {
    all_subsegs.clear();
    m_SubsegMap->get_all_subsegs(all_subsegs);
  }
  //
  void
  replaceSubsegVert(Vert * oldVert, Vert * newVert)
  {
    m_SubsegMap->replaceVert(oldVert, newVert);
  }
  //  void updateSubsegVertDataAfterPurge(const std::map<Vert*, Vert*>& vert_map) { if(m_SubsegMap) m_SubsegMap->update_vert_data_after_purge(vert_map);}

  // standard check for validity
  // fairly extensive.
  bool
  subsegMapIsValid() const
  {
    return m_SubsegMap->valid();
  }

  // used to take in an old subseg and create new subsegs based on the given information
  // performs the split, etc.
  std::pair<Subseg*, Subseg*>
  createNewSubsegsBySplit(Subseg * oldSubseg, const double splitCoords[3],
			  const double splitParam, Vert* const newVert);

  // formerly a RefinementManager3D member function, moved here because this is required for Length3D.
  void
  initializeSubsegs(const set<Subseg*>* const given_subsegs = NULL);

  // strictly a debugging thing, names it SubsegMap.vtk
  void
  printSubsegMap() const
  {
#ifndef NDEBUG
    m_SubsegMap->print();
#endif
  }
  // only new function, potentially dangerous,
  // Only used to pass the subsegMap Pointer to the length3D object
  // TODO figure out a less intrusive, safer way to call this
  SubsegMap *
  getSubsegMap()
  {
    assert(m_SubsegMap);
    return m_SubsegMap;
  }

	void initVolMeshFromSurfMesh(SurfMesh& SM);
};

class EdgeToBFace {
private:
  Vert *pV0, *pV1;
  BFace* pBF;
public:
  EdgeToBFace() :
      pV0(pVInvalidVert), pV1(pVInvalidVert), pBF(NULL)
  {
  }
  EdgeToBFace(const EdgeToBFace& EBF) :
      pV0(EBF.pV0), pV1(EBF.pV1), pBF(EBF.pBF)
  {
  }
  EdgeToBFace(Vert * const pV0In, Vert * const pV1In, BFace * const pBFIn =
  NULL) :
      pV0(pV0In), pV1(pV1In), pBF(pBFIn)
  {
    assert(pV0In->isValid() && pV1In->isValid());
    if (pV0 > pV1) {
      Vert *pVTmp = pV1;
      pV1 = pV0;
      pV0 = pVTmp;
    }
  }
  EdgeToBFace&
  operator=(const EdgeToBFace& EBF)
  {
    if (this != &EBF) {
      pV0 = EBF.pV0;
      pV1 = EBF.pV1;
      pBF = EBF.pBF;
    }
    return (*this);
  }
  Vert *
  pVVert(const int i) const
  {
    assert(i == 0 || i == 1);
    return (i == 0) ? pV0 : pV1;
  }
  //  BdryPatch* pPatchPointer() const {return pBF;}
  BFace *
  pBFace() const
  {
    return pBF;
  }
  bool
  operator!=(const EdgeToBFace& EBF) const
  {
    return !(*this == EBF);
  }
  bool
  operator==(const EdgeToBFace& EBF) const
  {
    return ((pV0 == EBF.pV0 && pV1 == EBF.pV1)
	|| (pV0 == EBF.pV1 && pV1 == EBF.pV0));
  }
  bool
  operator<(const EdgeToBFace& EBF) const
  {
    return ((pV0 < EBF.pV0) || (pV0 == EBF.pV0 && pV1 < EBF.pV1));
  }
};

// Protect small dihedral angles in a volume mesh before beginning insertion,
// to prevent runaway insertion in the small angle.
//void protectSmallAngles(VolMesh* const pVM);

// Check to see whether a VolMesh correctly conforms to a surface mesh.
bool
isBdryConstrained(VolMesh* pVM, SurfMesh* pSM,
		  std::set<ConstrainEdge>& sCESubSegs,
		  std::set<ConstrainEdge>& sCEInFacets,
		  std::set<ConstrainFace>& sCFSubFacets);

void
readUGridBinary(const char * const fileName, VolMesh * const pVM,
		int *& bdryCond);
void
readVGridBinary(const char * const fileName, VolMesh * const pVM);
void
readVTK_ASCII(const char * const fileName, VolMesh * const pVM);
void
readMEdit3D(const char * const fileName, VolMesh * const pVM);
void
readVMeshFile(const char strBaseFileName[], VolMesh * const pVM,
	      int *& bdryCond);
void
readTetGenFile(const char strBaseFileName[], VolMesh * const pVM,
	       int*& aiBCList, int* numFields, double** fieldData);

void
writeVTKLegacy(VolMesh& OutMesh, const char strBaseFileName[],
	       const char strExtraFileSuffix[] = "",
	       const int *cellData = NULL);
void
readTecPlotASCII(const char* const fileName, VolMesh * const pVM);

#ifdef HAVE_MESQUITE
void
writeVTKLegacy3DTMOPQual(VolMesh& OutMesh, const char strBaseFileName[],
			 const char strExtraFileSuffix[] = "",
			 const int *cellData = NULL);
#endif

/// These need to be moved at some point to somewhere, since they are fairly dimensionless
void
writeVTK_CellSet(std::vector<Cell*> & cell_list, const char strBaseFileName[]);
void
writeVTK_SingleFace(const Face * face, const char strBaseFileName[]);

// writes a point set stored in a vector of CubitVectors
void
writeVTKPointSet(std::vector<CubitVector> &cv_list, double * data,
		 const char strBaseFileName[], const char strExtraFileSuffix[]);

/// Write a file in GRUMMP native format. Picks whether to write with or
/// without internal bdry faces.
void
writeNative(VolMesh& OutMesh, const char strBaseFileName[],
	    const bool qHighOrder = false,
	    const char strExtraFileSuffix[] = "");
void
writeNativeMESH(VolMesh & Mesh, const char strBaseFileName[],
		const char strExtraFileSuffix[] = "");
void
writeNativeGRMESH(VolMesh & Mesh, const char strBaseFileName[],
		  const char strExtraFileSuffix[] = "");
void
writeFile_PWM(VolMesh & Mesh, const char strBaseFileName[],
	      const char strExtraFileSuffix[] = "");
void
writeFileTetGen(VolMesh& OutMesh, const char strBaseFileName[],
		const char strInfix[], const int numVars,
		const double* const solnData);

VolMesh* createVolMeshFromBdryFile(char strBaseFileName[FILE_NAME_LEN],
		int iBdryOptLevel, int iQualMeasure, double dResolution,
		double dGrading, Bdry3D*& pB3D);

void
readFile_Native(const char * const strBaseFileName, GR_index_t& iNumVerts,
		GR_index_t& iNumFaces, GR_index_t& iNumCells,
		GR_index_t& iNumBdryFaces, GR_index_t& iNumIntBdryFaces,
		bool& qFaceVert, bool& qFaceCell, bool& qCellRegion,
		bool& qBFaceFace, bool& qBFaceVert, bool& qBFaceBC,
		bool& qIntBFaceFace, bool& qIntBFaceVert, bool& qIntBFaceBC,
		EntContainer<Vert>& EC, GR_index_t (*&a2iFaceVert)[3],
		GR_sindex_t (*&a2iFaceCell)[2], int *&aiCellRegion,
		GR_index_t *&aiBFaceFace, int *&aiBFaceBC,
		GR_index_t (*&a2iBFaceVert)[3],
		GR_index_t (*&a2iIntBFaceFace)[2], int *&aiIntBFaceBC,
		GR_index_t (*&a2iIntBFaceVert)[3]);

void
repairBadCells(VolMesh* pVM);
void
recoverSurface(VolMesh* const pVM, SurfMesh* const pSM);
int
identifyNonDelaunayCells(VolMesh& VM, bool aqDelaunay[]);
#endif
