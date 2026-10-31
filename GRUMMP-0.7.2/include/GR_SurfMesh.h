#ifndef GR_SurfMesh
#define GR_SurfMesh 1

#include "GR_config.h"
#include "GR_Mesh2D.h"
#include "GR_EntContainer.h"
#include "GR_Subseg.h"
///
typedef std::vector<std::pair<FaceQueueEntry*, RefFace*> > FaceVec;

class SurfMesh : public Mesh2D {
  // MultiEdge's are needed to allow non-manifold surface meshes.
  EntContainer<MultiEdge> m_ECME;
  /// surface geometry data
  Bdry3D *m_Bdry3D;
  /// Tie-in from cells in surface mesh to boundary patches.  This array
  // should be variable-sized, but I couldn't be bothered on the first
  // cut at this code.  (CFO-G, 11/2000)
  BdryPatch3D **m_bdryPatches;
  GR_index_t m_numBdryPatches;

  /// this will handle curved surfaces
  SubsegMap * m_SubsegMap;
  // all boundary info will be passed in somehow
  // need a way to convert
  // going to carry one to one mapping between surfaces and
  // TriCells. Each TriCell has its own entry in this
  RefFace **m_surfaces;
  GR_index_t m_numSurfaces;
private:
  /// This constructor might be helpful for volume meshing if we cared
  /// about the surface mesh afterwards; the volume initialization
  /// doesn't currently clean up the surface mesh if points are inserted
  /// in it.
  SurfMesh(const SurfMesh&);
  /// Op= isn't needed either
  SurfMesh&
  operator=(const SurfMesh&);

public:
  /// This is the default constructor
  SurfMesh(const int iQualMeas = 2) :
      Mesh2D(iQualMeas), m_ECME(), m_Bdry3D(NULL), m_bdryPatches(NULL), m_numBdryPatches(
	  0), m_SubsegMap(NULL), m_surfaces(NULL), m_numSurfaces(0)
  {
  }
  ///
  SurfMesh(Bdry3D * const pB3D, const int iQualMeas = 2);
  ///
  SurfMesh(const char strBaseFileName[], const int iQualMeas = 2);
  ///
  SurfMesh(FaceVec faces, std::set<Subseg*> subsegs, const int iQualMeas = 2);

  SurfMesh(EntContainer<Vert>& EC, const int iQualMeas = 2) :
      Mesh2D(iQualMeas), m_ECME(), m_Bdry3D(NULL), m_bdryPatches(NULL), m_numBdryPatches(
	  0), m_SubsegMap(NULL), m_surfaces(), m_numSurfaces(0)
  {
    // Don't pass the vertex list to the Mesh2D constructor, because
    // it will try to triangulate the point set as if it were 2D.
    for (GR_index_t i = 0; i < EC.lastEntry(); i++) {
      Vert *pVOld = EC.getEntry(i);
      Vert *pVNew = m_ECVerts.getNewEntry();
      pVNew->setCoords(3, pVOld->getCoords());
      pVNew->copyAllFlagsFrom(pVOld);
      // Don't copy vertex data any deeper than this.
    }

    m_bdryPatches = new BdryPatch3D*[50];
    m_numBdryPatches = 50;
    // Can't set up faces and cells here, because there's not enough
    // info to do this correctly.
  }
  ///
  virtual
  ~SurfMesh()
  {
    if (m_bdryPatches) {
      delete[] m_bdryPatches;
      m_bdryPatches = NULL;
    }
    if (m_SubsegMap)
      delete m_SubsegMap;
    if (m_surfaces) {
      delete[] m_surfaces;
      m_surfaces = NULL;
    }
  }
  ///
  void
  setVertFaceNeighbors();
  virtual Vert*
  createVert(const double adCoords[]);
  virtual Vert*
  createVert(const double dX, const double dY, const double dZ);
  MultiEdge*
  createMultiEdge(Face* const pF);
  MultiEdge*
  createMultiEdge(Vert* const pV0, Vert* const pV1);
  using Mesh2D::createTriCell;
//  virtual TriCell* createTriCell(int minlayerindex, Vert * const pV0, Vert * const pV1,
//				 Vert * const pV2, const int iReg = iDefaultRegion,
//				 const bool autoRotate = false, int ID=0);
  virtual TriCell*
  createTriCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
		const int iReg = iDefaultRegion, const bool autoRotate = false,
		int ID = 0);
  virtual TriCell*
  createTriCell(Face * const pF0, Face * const pF1, Face * const pF2,
		const int iReg = iDefaultRegion, int ID = 0);
  virtual QuadCell*
  createQuadCell(Face * const pF0, Face * const pF1, Face * const pF2,
		 Face * const pF3, const int iReg = iDefaultRegion);
  virtual QuadCell*
  createQuadCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
		 Vert * const pV3, const int iReg = iDefaultRegion);

  virtual GR_index_t
  getNumFaces() const
  {
    return m_ECEdgeF.lastEntry() + m_ECME.lastEntry();
  }
  virtual GR_index_t
  getFaceIndex(const Face* const pF) const;
  virtual Face*
  getFace(const GR_index_t i) const;
  ///
  GR_index_t
  getNumEdgeFaces() const
  {
    return m_ECEdgeF.lastEntry();
  }
  GR_index_t
  getNumMultiEdges() const
  {
    return m_ECME.lastEntry();
  }
  ///
  MultiEdge*
  getNewMultiEdge();
  ///
  MultiEdge*
  getMultiEdge(const GR_index_t i) const
  {
    assert(i < m_ECME.lastEntry());
    return m_ECME.getEntry(i);
  }
  ///
  void
  convertFromCellVert(EntContainer<TriCellCV>& ECTriCV);
  ///
  Bdry3D *
  getBoundaryData() const
  {
    return m_Bdry3D;
  }
  ///
  BdryPatch3D *
  getBdryPatch(const GR_index_t iCell) const
  {
    assert(iCell < getNumCells());
    return m_bdryPatches[iCell];
  }
  BdryPatch3D *
  getBdryPatch(const Cell* const pC) const
  {
    GR_index_t iCell = getCellIndex(pC);
    return m_bdryPatches[iCell];
  }
  void
  setBdryPatch(const GR_index_t iCell, BdryPatch3D* const pBP)
  {
    assert(pBP != NULL);
    if (m_numBdryPatches == 0) {
      assert(m_bdryPatches == NULL);
      m_numBdryPatches = getNumCells() < 50 ? 50 : getNumCells();
      m_bdryPatches = new BdryPatch3D*[m_numBdryPatches];
    }
    if (iCell >= m_numBdryPatches) {
      assert(m_numBdryPatches > 0);
      BdryPatch3D **apBPTemp = m_bdryPatches;
      GR_index_t iNewSize = m_numBdryPatches * 2;
      while (iNewSize < iCell)
	iNewSize *= 2;
      m_bdryPatches = new BdryPatch3D*[iNewSize];
      // Copy the data into the new array
      for (GR_index_t iTmp = 0; iTmp < m_numBdryPatches; iTmp++) {
	m_bdryPatches[iTmp] = apBPTemp[iTmp];
      }
      m_numBdryPatches = iNewSize;
      delete[] apBPTemp;
    }
    m_bdryPatches[iCell] = pBP;
  }
  void
  setBdryPatch(const Cell* const pC, BdryPatch3D* const pBP)
  {
    GR_index_t iCell = getCellIndex(pC);
    setBdryPatch(iCell, pBP);
  }
  ///
//  void
//  createPatches();
  void
  setSurface(const GR_index_t iCell, RefFace * const refFace)
  {
    assert(refFace != NULL);
    if (m_numSurfaces == 0) {
      assert(m_surfaces == NULL);
      m_numSurfaces = getNumCells() < 50 ? 50 : getNumCells();
      m_surfaces = new RefFace*[m_numSurfaces];
    }
    if (iCell >= m_numSurfaces) {
      assert(m_numSurfaces > 0);
      RefFace **surfaceTemp = m_surfaces;
      GR_index_t iNewSize = m_numSurfaces * 2;
      while (iNewSize < iCell)
	iNewSize *= 2;
      m_surfaces = new RefFace*[iNewSize];
      // Copy the data into the new array
      for (GR_index_t iTmp = 0; iTmp < m_numSurfaces; iTmp++) {
	m_surfaces[iTmp] = surfaceTemp[iTmp];
      }
      m_numSurfaces = iNewSize;
      delete[] surfaceTemp;
    }
    m_surfaces[iCell] = refFace;
  }

  void
  setSurface(const Cell * const pC, RefFace * const refFace)
  {
    GR_index_t iCell = getCellIndex(pC);
    setSurface(iCell, refFace);
  }
  RefFace *
  getSurface(const GR_index_t iCell) const
  {
    assert(iCell < getNumCells());
    return m_surfaces[iCell];
  }
  RefFace *
  getSurface(const Cell* const pC) const
  {
    GR_index_t iCell = getCellIndex(pC);
    return m_surfaces[iCell];
  }
  void
  getSubsegs(std::set<Subseg*>& all_subsegs)
  {
    m_SubsegMap->get_all_subsegs(all_subsegs);
  }
  void
  addSubseg(Subseg * subseg)
  {
    m_SubsegMap->add_subseg(subseg);
  }
  Subseg*
  areSubsegVertsConnected(Vert* const vert0, Vert* const vert1) const
  {
    if (vert0->getVertType() == Vert::eInterior
	|| vert1->getVertType() == Vert::eInterior)
      return NULL;
    return m_SubsegMap->verts_are_connected(vert0, vert1);
  }
  void
  writeSubsegMap()
  {
#ifndef NDEBUG
    m_SubsegMap->print();
#endif
  }
public:
  ///
  eMeshType
  getType() const
  {
    return eSurfMesh;
  }
  ///
  int
  reconfigure(Face*& pF);
  /*   /// */
  /*   bool qInsertPoint(const double adPoint[3], Cell* const pC, */
  /* 		   int * const piSwaps, const bool qSwap = true, */
  /* 		   const bool qForce = false, Vert* pVNew = pVInvalidVert); */
  ///
  int
  insertOnFace(Vert* const pVNew, Face* const pF, Cell* const pC,
	       const bool qSwap, bool* qInsertOnEdge = NULL);
  ///
  int
  insertOnMultiEdge(Vert* const pVNew, MultiEdge* const pME, const bool qSwap);
  /*   /// */
  /*   int iInsertInInterior(Vert* const pVNew, Cell* const pC, */
  /* 			const bool qSwap); */
  ///
  void
  markSmallAngles();
  ///
//	void protectSmallAngles();
  ///
  int
  iFaceSwap_deprecated(Face*& pF);
  ///
  virtual bool
  doSwap(const Vert* const pVVertA, const Vert* const pVVertB,
	 const Vert* const pVVertC, const Vert* const pVVertD,
	 const Vert* const pVVertE = pVInvalidVert) const;

  virtual bool
  doSwap(double dLenDiag, double dLenCL1, double dLenCL2, double dLenCR1,
	 double dLenCR2);

public:
  ///
  double
  calcBoundarySize() const;
  ///
  double
  calcInteriorSize() const;
  ///
  bool
  isWatertight() const;
  ///
  void
  purgeVerts(std::map<Vert*, Vert*>* vert_map = NULL);
  /// need to rearrange cells.
  virtual void
  purgeAllEntities(std::map<Vert*, Vert*>* vert_map = NULL,
		   std::map<Face*, Face*>* face_map = NULL,
		   std::map<Cell*, Cell*>* cell_map = NULL,
		   std::map<BFace*, BFace*>* bface_map = NULL);

};

void
writeFile_Surface(SurfMesh & Mesh, const char strBaseFileName[],
		  const char strExtraFileSuffix[] = "");

void
readFile_Surface(const char * const strBaseFileName, GR_index_t& iNumVerts,
		 GR_index_t& iNumFaces, GR_index_t& iNumCells,
		 GR_index_t& iNumIntBdryFaces, bool& qFaceVert, bool& qCellVert,
		 bool& qFaceCell, bool& qCellFace, bool& qCellRegion,
		 bool& qIntBFaceFace, bool& qIntBFaceVert, bool& qIntBFaceBC,
		 EntContainer<Vert>& ECVerts, GR_index_t (*&a2iFaceVert)[2],
		 GR_sindex_t (*&a2iFaceCell)[2], GR_index_t (*&a2iCellVert)[3],
		 GR_index_t (*&a2iCellFace)[3], int *&aiCellRegion,
		 GR_index_t (*&a2iIntBFaceFace)[2], int *&aiIntBFaceBC,
		 GR_index_t (*&a2iIntBFaceVert)[2]);

#endif
