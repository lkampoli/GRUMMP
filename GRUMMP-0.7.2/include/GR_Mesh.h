#ifndef GR_Mesh
#define GR_Mesh 1

#include <map>
#include <set>
//#include <omp.h>

#include "GR_CellQueue.h"
#include "GR_config.h"
#include "GR_BFace.h"
#include "GR_Cell.h"
#include "GR_EntContainer.h"
#include "GR_Face.h"
#include "GR_Ghost.h"
#include "GR_Length.h"
#include "GR_MPI_Defs.h"
#include "GR_Observable.h"
#include "GR_Quality.h"
#include "GR_Vertex.h"
#include <cstring>

enum eSwap {
	eNoSwap, eDelaunay, eMaxDihed, eMinSine, eFromQual, eRandom
};

class CellQueue;
class AnisoRefinement;

/**\brief This is the base class for all mesh data.
 *
 * Vertex data is the primary thing stored internally; connectivity data is
 * declared in the derived classes.
 */
// The publishing side of the Observer paradigm is implemented in the
// Observable class.
class Mesh: public GRUMMP::Observable {
public:
	enum ExtraData {
		eStandard = 0, eMetric = 1, eTMOPQual = 2
	};
private:
	/// Operator= disallowed.
	Mesh&
	operator=(const Mesh&) {
		assert(0);
		return (*this);
	}
protected:
	std::set<Vert*> m_vertsToUpdateLenScale;
	/// Variable-sized array for storing vertex data.
	EntContainer<Vert> m_ECVerts;
#ifdef HAVE_MPI
	/// A container for all the entities (all types) that have a copy here and
	/// somewhere else.  All of these entities will be stored in the normal
	/// containers.
	std::map<Entity*, GRUMMP::Ghost> m_ghostMap;
#endif // HAVE_MPI
	//  void m_ECVerts_sync_entry () {m_ECVerts.sync_entry();}
	//  void m_ECVerts_sync_exit () {m_ECVerts.sync_exit();}
	/// Size factor: this is the approx number of edges per
	/// feature size.
	double m_resolution;
	/// Grading constant
	double m_LipschitzAlpha;
	/// Quality measure.  Will eventually not be a Mesh member.
	mutable Quality *m_qual;
protected:
	///  If true, deleted verts are kept in the database.
	bool m_keepPhantomVerts :1;
	///
	bool m_swapRecursively :1;
	///
	bool m_allowBdryChanges :1;
	///
	bool m_isLengthScaleInitialized :1;
	bool m_isLengthScaleFromCellSizes :1;
	bool m_interpolateLengthScale :1;
	///
	bool m_simplexMesh :1;
	///
	eSwap m_swapMeasure :3;
	///
	eEncroachType m_encroachType :3;
	///
	/// Pointer to Cell Queue
public:
	// These next two variables are strictly for 2D aniso insertion.
	// TODO: Factor these out of Mesh, and even out of Mesh2D.
	CellQueue* pCellQueue;          //A pointer to the meshes InsertionQueue
	AnisoRefinement* pAnisoRef; //A pointer to the meshes Anisotropic Refinement
	//Pointer to a set of Faces
	//  std::set<Face*>* psetPFaces;  //used in anisotropic meshing
/// Constructor and destructor

protected:
	// an experiment in visualization,
	bool m_writeTempMeshes :1;
	char m_tempMeshFileName[FILE_NAME_LEN];
	GR_index_t m_tempMeshIndex;
public:
	enum eMeshType {
		eMesh2D, eSurfMesh, eVolMesh
	};
	enum eInsertType {
		eBadCell, eBdryFaceEncroach, eBdrySegEncroach
	};
	///
	Mesh() :
			m_vertsToUpdateLenScale(), m_ECVerts(), m_resolution(1), m_LipschitzAlpha(1), m_qual(NULL), m_keepPhantomVerts(
					false), m_swapRecursively(true), m_allowBdryChanges(true), m_isLengthScaleInitialized(
					false), m_isLengthScaleFromCellSizes(false), m_interpolateLengthScale(
					false), m_simplexMesh(true), m_swapMeasure(eDelaunay), m_encroachType(
					eLens), pCellQueue(NULL), pAnisoRef(NULL), m_writeTempMeshes(
					false), m_tempMeshIndex(1)
#ifdef HAVE_MPI
					, m_PartID(), m_partition()
#endif
	{
		strcpy(m_tempMeshFileName, "tempMesh");
		// The value of m_interiorSkip is the correct default for all meshes.
		// The value of m_fillBdry is the proper default for 2D and surface;
		// the 3D constructors override this.
		// Can't initialize Quality at the Mesh level.
	}
	///
	Mesh(const EntContainer<Vert>& ECPts) :
			m_vertsToUpdateLenScale(), m_ECVerts(), m_resolution(1), m_LipschitzAlpha(1),
				m_qual(NULL), m_keepPhantomVerts(
					false), m_swapRecursively(true), m_allowBdryChanges(true), m_isLengthScaleInitialized(
					false), m_isLengthScaleFromCellSizes(false), m_interpolateLengthScale(
					false), m_simplexMesh(true), m_swapMeasure(eDelaunay), m_encroachType(
					eLens), pCellQueue(NULL), pAnisoRef(NULL), m_writeTempMeshes(
					false), m_tempMeshIndex(1)
#ifdef HAVE_MPI
					, m_PartID(), m_partition()
#endif
	{
		// The value of m_interiorSkip is the correct default for all meshes.
		// The value of m_fillBdry is the proper default for 2D and surface;
		// the 3D constructors override this.
		copyVertData(ECPts);
		strcpy(m_tempMeshFileName, "tempMesh");
	}
	///
	Mesh(const Mesh& M);
	///
	virtual ~Mesh() {
		if (m_qual)
			delete m_qual;
	}
	///
	virtual void
	readFromFile(const char strBaseFileName[], int* numVars = nullptr,
			double** data = nullptr) = 0;
/// Inquiry functions
public:

	////////////////////////////////////////////////////////
	/// Functions that directly query the mesh database. ///
	////////////////////////////////////////////////////////

	// This function returns a pointer to internal data.  This isn't a great idea.
	// TODO: Come up with a better alternative method to replace the three calls
	// to this function.
	EntContainer<Vert>*
	pECGetVerts() {
		return (&m_ECVerts);
	}
	///checks all verts until (if) it finds pV. slow...
	bool
	hasVert(const Vert* const pV) const;
	///
	GR_index_t getVertIndex(const Vert* const pV) const {
		return (m_ECVerts.getIndex(pV));
	}
	///
	virtual GR_index_t
	getFaceIndex(const Face* const pF) const = 0;
	///
	virtual GR_index_t
	getBFaceIndex(const BFace* const pBF) const = 0;
	///
	virtual GR_index_t
	getCellIndex(const Cell* const pC) const = 0;
#ifdef ITAPS
	///
	virtual bool
	isValidEntHandle(void* pvEnt) const = 0;
#endif
	///
	virtual eMeshType
	getType() const = 0;
	///
	/// perimeter length or  area:
	virtual double calcBoundarySize() const {
		assert(0);
		return 0;
	}
	/// interior   area  or volume:
	virtual double calcInteriorSize() const {
		assert(0);
		return 0;
	}
	virtual bool isWatertight() const {
		assert(0);
		return 0;
	}
	///
	/// Validity checking
	virtual bool
	isValid() const = 0;
//
	  int
	  calcNumFacesNonDelaunay() const;
	/////////////////////////////////////////////////////////////////
	/// Primitive functions that modify mesh geometry or topology ///
	/////////////////////////////////////////////////////////////////
	/// Create a vertex at the given location.
	virtual Vert*
	createVert(const double dX, const double dY, const double dZ) = 0;
	virtual Vert*
	createVert(const double adCoords[]) = 0;
	/*   virtual Face* createFace(Vert * const apV[], const int iNVerts) = 0; */
	virtual BFace*
	createBFace(Face * const pF, BFace * const pBF =
	pBFInvalidBFace) = 0;
	/*   virtual Cell* createCell(Face * const apF[], const int iNFaces, */
	/* 			   const int iReg = iDefaultRegion) = 0; */
	/*   virtual Cell* createCell(Vert * const apV[], const int iNVerts, */
	/* 			   const int iReg = iDefaultRegion) = 0; */
	virtual Cell*
	createCell(Vert* const vert, Face* const face, const int iReg =
			iDefaultRegion) = 0;
	virtual Cell*
	createCell(int minlayerindex, Vert* const vert, Face* const face,
			const int iReg = iDefaultRegion) = 0;
	virtual bool
	deleteVert(Vert * const pV) = 0;
	virtual bool
	deleteFace(Face * const pF) = 0;
	virtual bool
	deleteBFace(BFace * const pBF) = 0;
	virtual bool
	deleteCell(Cell * const pC) = 0;

	///////////////////////////////////////
	/// Database maintainence functions ///
	///////////////////////////////////////

	/// The following routine adds bdry faces where needed after ITAPS
	/// modification.  Doesn't put BFaces on the right patches, though...
	void
	checkForITAPSBdryErrors();
	///
	virtual void
	purgeBdryFaces(std::map<BFace*, BFace*>* bface_map = NULL) = 0;
	///
	virtual void
	purgeFaces(std::map<Face*, Face*>* face_map = NULL) = 0;
	///
	virtual void
	purgeCells(std::map<Cell*, Cell*>* cell_map = NULL) = 0;
	///
	virtual void
	purgeVerts(std::map<Vert*, Vert*>* vert_map = NULL) = 0;
	///
	virtual void
	purgeAllEntities(std::map<Vert*, Vert*>* vert_map = NULL,
			std::map<Face*, Face*>* face_map = NULL,
			std::map<Cell*, Cell*>* cell_map = NULL,
			std::map<BFace*, BFace*>* bface_map = NULL) = 0;

	enum eReorderingType {
		RCM, Lexicographic
	};
	void reorderAllEntities(const eReorderingType type = RCM) {
		purgeAllEntities();
		switch (type) {
		case RCM:
		default:
			reorderVerts_RCM();
			break;
		case Lexicographic:
			reorderVerts_Lex();
			break;
		}
		//purgeAllEntities();
		reorderCells();
		reorderFaces();

		//Modified by SG July 2007:
		//The following seems to be used to keep consistent face normal orientation.
		//However, the most recent iteration of the 2D tri meshing code with curved boundaries
		//needs to go through BdryEdge and IntBdryEdge objects to obtain normals. To achieve this,
		//the faces must be oriented in the SAME direction as their underlying curve.
		//Interchanging the vertices will in obviously change the face orientation and
		//therefore, the normal direction will become inconsistent. Hence, the skipping of this
		//piece of code for 2D triangular meshes.

		if (getType() == Mesh::eMesh2D)
			return;
		for (GR_index_t iF = 0; iF < getNumFaces(); iF++) {
			Face *pF = getFace(iF);
			if (pF->isBdryFace() && pF->getLeftCell()->isBdryCell()) {
				pF->interchangeCellsAndVerts();
			}
		}

	}

private:
	void
	reorderVerts_Lex(void);
public:
	void
	reorderVertices(GR_index_t*& aRReorderInt);
protected:
	virtual void
	reorderVerts_RCM(void) = 0;
	virtual void
	reorderCells(void) = 0;
	virtual void
	reorderFaces(void) = 0;

	void copyVertData(const EntContainer<Vert>& ECPts) {
		// The value of m_interiorSkip is the correct default for all meshes.
		// The value of m_fillBdry is the proper default for 2D and surface;
		// the 3D constructors override this.
		int iDim = ECPts.getEntry(0)->getSpaceDimen();
		for (GR_index_t i = 0; i < ECPts.lastEntry(); i++) {
			Vert* pVOld = ECPts.getEntry(i);
			Vert* pVNew = m_ECVerts.getNewEntry();
			pVNew->setCoords(iDim, pVOld->getCoords());
			pVNew->copyAllFlagsFrom(pVOld);
			// Don't copy vertex data any deeper than this.
		}
	}

public:
	///
	void
	setAllHintFaces();
	void
	clearAllHintFaces();

	virtual void
	setVertFaceNeighbors();
	void
	clearVertFaceNeighbors();

	////////////////////////////////////////////////////////
	/// Functions that affect the behavior of algorithms ///
	////////////////////////////////////////////////////////

	///
	eSwap getSwapType() const {
		return m_swapMeasure;
	}
	///
	void setSwapType(const eSwap eSwapIn) {
		m_swapMeasure = eSwapIn;
	}
	///
	bool areBdryChangesAllowed() const {
		return m_allowBdryChanges;
	}
	///
	void allowBdryChanges() {
		m_allowBdryChanges = true;
	}
	///
	void disallowBdryChanges() {
		m_allowBdryChanges = false;
	}
	///
	bool isSwapRecursionAllowed() const {
		return m_swapRecursively;
	}
	///
	void allowSwapRecursion() {
		m_swapRecursively = true;
	}
	///
	void disallowSwapRecursion() {
		m_swapRecursively = false;
	}
	///
	void setLengthScaleFromCellSizes() {
		m_isLengthScaleFromCellSizes = true;
	}
	///
	void unsetLengthScaleFromCellSizes() {
		m_isLengthScaleFromCellSizes = false;
	}
	///

//	void initLengthScale(const double dRes = 1, const double dGrade = 1);
	void setInterpolatedLengthScale() {
		m_interpolateLengthScale = true;
	}
	bool isInterpolatedLengthScaleUsed() {
		return m_interpolateLengthScale;
	}
	void
	checkEdgeLengths() const;
	/// Return the number of vertices in the database; some may be deleted.
	GR_index_t getNumVerts() const {
		return m_ECVerts.lastEntry();
	}

	GR_index_t getNumVertsInUse() const {
		return m_ECVerts.size();
	}
	virtual GR_index_t
	getNumEdgesInUse() const = 0;
	virtual GR_index_t
	getNumTrisInUse() const = 0;
	virtual GR_index_t
	getNumQuadsInUse() const = 0;
	virtual GR_index_t
	getNumTetsInUse() const = 0;
	virtual GR_index_t
	getNumPyrsInUse() const = 0;
	virtual GR_index_t
	getNumPrismsInUse() const = 0;
	virtual GR_index_t
	getNumHexesInUse() const = 0;

	/// Return the number of faces in the database; some may be deleted.
	virtual GR_index_t
	getNumFaces() const = 0;
	///
	virtual GR_index_t
	getNumBdryFaces() const = 0;
	///
	virtual GR_index_t
	getNumIntBdryFaces() const = 0;

	///
	GR_index_t getNumTotalBdryFaces() const {
		return getNumBdryFaces() + getNumIntBdryFaces();
	}

	/// Return the number of cells in the database; some may be deleted.
	virtual GR_index_t
	getNumCells() const = 0;

	/// Return a cell by index.
	virtual Cell*
	getCell(const GR_index_t i) const = 0;

	/// Return a vertex by index.
	Vert*
	getVert(const GR_index_t i) const {
		return m_ECVerts.getEntry(i);
	}

	/// Return a face by index.
	virtual Face*
	getFace(const GR_index_t i) const = 0;

	/// Return the handle of an unused vertex in the database.
	Vert*
	getNewVert(int ID = 0) {
		Vert *pV;
		if (omp_in_parallel() == 0)
			pV = m_ECVerts.getNewEntry();
		else {
			assert(ID == omp_get_thread_num());
			pV = m_ECVerts.getNewEntry(ID);
		}
		pV->setDefaultFlags();
		pV->clearFaceConnectivity();
		return (pV);
	}

	///
	virtual BFace*
	getBFace(const GR_index_t i) const = 0;
	virtual BFace*
	getIntBFace(const GR_index_t i) const = 0;

	/// Conversion to/from simplicial form
	bool isSimplicial() const {
		return m_simplexMesh;
	}
	///
	virtual void makeSimplicial() {
		assert(0);
	}
	///
	void allowNonSimplicial() {
		m_simplexMesh = false;
	}

	bool arePhantomVertsKept() const {
		return m_keepPhantomVerts;
	}

	void keepPhantomVerts(bool keepPhantoms) {
		m_keepPhantomVerts = keepPhantoms;
	}

	///
	void setQualityMeasure(const int iQual) {
		delete m_qual;
		m_qual = new Quality(this, iQual);
	}
	void evaluateQuality() const {
		m_qual->vEvaluate();
	}
	double evaluateQuality(const CellSkel* const pCS) const {
		return m_qual->dEvaluate(pCS);
	}
	void writeQualityToFile(const char strQualFileName[]) const {
		m_qual->vWriteToFile(strQualFileName);
	}

	// The next three functions are used only for Doug's aniso insertion.
	void attachCellQueue(CellQueue* pCellQueueIn) {
		pCellQueue = pCellQueueIn;
	}
	void removeCellQueue() {
		pCellQueue = NULL;
	}

	void attachAnisoRef(AnisoRefinement* pAnisoRefIn) {
		pAnisoRef = pAnisoRefIn;
	}

	/*   void vAttachFaceSet(std::set<Face*>* psetPFacesIn){assert(psetPFaces==NULL);psetPFaces=psetPFacesIn;} */
	/*   void vRemoveFaceSet(){psetPFaces=NULL;} */

public:
	///
	void setEncroachmentType(const eEncroachType eET) {
		m_encroachType = eET;
	}
	eEncroachType getEncroachmentType() const {
		return m_encroachType;
	}

	void meltIntBdryFaces(std::vector<Face*>& to_lock);

	void write_vertInfo_for_pFLeft_faces(std::vector<Face*>& faces);
	void lock_faces(std::vector<Face*>& faces);

	// The worst shape quality measure (shortest edge to circumradius) allowed.
	virtual double
	getWorstAllowedCellShape() const = 0;

//#ifdef HAVE_MESQUITE
//  void
//  addToCoarseVertsAnisoLength(const VertConnect aVC[]) const;
//  void
//  selectCoarseVertsAnisoLength(VertConnect aVC[], double dMinMLength);
//#endif
//
//  void
//  deactivateAll() const;
//  int
//  activate(const Vert::VertType VT) const;
	virtual Cell*
	beginPipe(const Vert* const pVLast, const Vert* const pVNew) const = 0;

public:
	virtual double
	getCoarseningConstant() const = 0;
	virtual void
	identifyVertexTypesGeometrically() const = 0;
//  void
//  initConflictGraph(VertConnect aVC[]) const;
//#ifdef HAVE_MESQUITE
//  void
//  vConflictGraphAnisoLength(double dMinMLength);
//  void
//  updateConflictGraphAnisoLength(VertConnect aVC[], double dMinMLength) const;
//#endif

	void
	removeTaggedVerts(bool smooth = true, std::set<Vert*> *cantRemove =
	NULL);
	virtual GR_index_t
	getNumBdryPatches() = 0;
	virtual void buildBdryPatches() {
	}
	virtual bool
	removeVert(Vert * const pV, int& iNewSwaps, Vert * apVPreferred[] = NULL,
			const int iNPref = 0) = 0;

// Mesh adaptation and related routines
	virtual bool
	makeDelaunay() = 0;
	virtual void
	adaptToLengthScale(const eEncroachType eET,
			std::shared_ptr<GRUMMP::Length> const sizing_field = NULL);
	// Iterator functions
	EntContainer<Vert>::iterator vert_begin() const {
		return m_ECVerts.begin();
	}
	EntContainer<Vert>::iterator vert_end() const {
		return m_ECVerts.end();
	}

	// This will output an unPurged version of the mesh
	virtual void
	writeTempMesh() = 0;
	void setWriteTempMesh(bool writeTempMeshes) {
		m_writeTempMeshes = writeTempMeshes;
	}

#ifdef HAVE_MPI

	// Functions for MPI meshing support
protected:
	GRUMMP::PartGID m_PartID;
	GRUMMP::Partition* m_partition;
public:
	void setPartID(int rank, GR_index_t localPartID) {
		m_PartID = GRUMMP::PartGID(rank, localPartID);
	}
	GRUMMP::PartGID getPartID() const {
		return m_PartID;
	}
	void setPartition(GRUMMP::Partition* partition) {
		m_partition = partition;
	}
	GRUMMP::Partition*
	getPartition() const {
		return m_partition;
	}

	bool
	addGhostEntry(Entity* const ent, const GRUMMP::PartGID& ownerID,
			const GR_index_t nCopies = 1);
	// Removes an entity from the Ghost map entirely
	void
	removeGhostEntry(Entity* const ent);
	// Resets info about copies for a ghosted entity
	void
	resetGhostEntry(Entity* const ent);
	bool
	addCopyToGhostEntry(Entity* const ent, const GRUMMP::EntityGID& remote,
			const bool isBdry = false);
	const GRUMMP::Ghost&
	getGhostEntry(Entity* const ent, bool& status);
	friend class GRUMMP::Partition;
	void
	matchPartBdrys();
	void
	setupGhostData();
	void
	clearGhostData();

private:
	void
	matchPartBdryVertices();
	void
	matchPartBdryFaces();

	void
	sendVertInfo(const enum GRUMMP::GhostMsgTags_t tag,
			const GRUMMP::PartGID& targetID, const std::set<Entity*>& verts);
	void
	replyWithVertIDs(const enum GRUMMP::GhostMsgTags_t tag,
			const std::set<Vert*>& vertsReceived);
#endif
};

/// Return info about entities in the one-cell ball around the vertex.
void
findNeighborhoodInfo(const Vert* const pVert, std::set<Cell*>& spCInc,
		std::set<Vert*>& spVNeigh, std::set<BFace*>* pspBFInc =
		NULL, bool *pqBdryVert = NULL, std::set<Face*>* pspFNearby = NULL);
// this function only should be used in debug mode, its meant to
// output the mesh without changing its ordering
// works by repeating vertices, and only outputting valid cells
// so each vertex is potentially written many times
// resulting in much larger files, but again, debug only
void
writeVTKLegacyWithoutPurge(const Mesh& OutMesh, const char strBaseFileName[]);

#endif // Done with class Mesh
