/*
 * GR_InsertionManager.h
 *
 *  Created on: May 16, 2012
 *      Author: cfog
 */

#ifndef GR_INSERTION_MANAGER_H_
#define GR_INSERTION_MANAGER_H_

// There are three levels of abstraction for insertion, because there are (at
// least) three important activities:
//
// 1.  At a high level, there's a driver that makes decisions about why and
//     where to add points to the mesh. (Generic name: Refiners; these are
//     the ones that are Observers of Mesh objects.)
// 2.  At a low level, there's code that does the actual insertion (valid mesh
//     to valid mesh).  (Generic name: Inserters)
// 3.  Finally, at the middle level, there's a need to be able to insert a point
//     with some additional guarantees on the properties of the mesh other than
//     just validity (that the mesh is always Delaunay, for instance).
//     (Generic name: InsertionManagers)
// An inserter needs to be able to do several things:
//
// insertPoint   This is the obvious, general purpose call.
//
// After some internal dispatch, as well as for special purposes, specific
// calls are helpful:
//
// insertPointInCell
// insertPointOnFace
// insertPointOnEdge
// insertPointOnBdryFace
// insertPointOnBdryEdge
//
// In 2D, the face and edge calls will be the same; if possible, make the edge
// calls private.  For Watson, it may be possible to make all the interior calls
// the same and the two bdry calls the same.
//
// Each also needs some internal state data:  the mesh it's operating on and a
// (typically Delaunay) SwapManager, at a minimum.

#include <map>
#include <queue>
#include "GR_InsertionQueue.h"
#include "GR_Length.h"
#include "GR_LengthAniso2D.h"
#include "GR_QualMeasure.h"
#include "GR_SwapManager.h"
#include "GR_Subseg.h"
#include "GR_SmoothingManager.h"
#include "GR_TMOPQual.h"
#include "GR_TreeConstructor.h"
#include <omp.h>

namespace GRUMMP
{
  class Inserter {
  protected:
    Mesh* m_pMesh;
    GR_index_t m_nInCell, m_nOnTri, m_nOnEdge, m_nOnBdryTri, m_nOnBdryEdge,
	m_nOnIntBdryTri, m_nOnIntBdryEdge;
    BdryEdgeBase::SplitType m_splitType;
    bool m_force;
  private:
    /// Default construction is a really, really bad idea.
    Inserter();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    Inserter(Inserter&);
    /// Operator= isn't obviously helpful, either.
    Inserter&
    operator=(Inserter&);
    //
  public:
    /// Constructor; the mesh must be supplied at construction time.
    Inserter(Mesh* pMesh) :
	m_pMesh(pMesh), m_nInCell(0), m_nOnTri(0), m_nOnEdge(0), m_nOnBdryTri(
	    0), m_nOnBdryEdge(0), m_nOnIntBdryTri(0), m_nOnIntBdryEdge(0), m_splitType(
	    BdryEdgeBase::MID_TVT), m_force(false)
    {
    }
    virtual
    ~Inserter()
    {
    }

    /// True iff points are required to be inserted at precisely the specified
    /// location.
    bool
    forceInsertion() const
    {
      return m_force;
    }

    /// Set whether to force insertion at the specified point; returns the old
    /// value of this parameter, in case you need to restore the value later.
    bool
    setForcedInsertion(const bool newForce)
    {
      bool retVal = m_force;
      m_force = newForce;
      return retVal;
    }

    virtual Vert*
    insertPoint(const double newPtLoc[], Cell* const cellGuess,
		Vert* const existingVert = NULL) = 0;

    virtual Vert*
    insertPointInCell(const double newPtLoc[], Cell* const cell,
		      Vert* const existingVert = NULL) = 0;
    virtual Vert*
    insertPointOnFace(const double newPtLoc[], Face* const face,
		      Vert* const existingVert = NULL) = 0;
    /*     virtual bool insertPointOnEdge(const double newPtLoc[], Face* const face, */
    /*                                    Vert* const vert0, Vert* const vert1) = 0; */

    /*     virtual bool insertPointOnBdry(const double newPtLoc[], */
    /*                                    Face*  const faceGuess) = 0; */
    virtual Vert*
    insertPointOnBdryFace(const double newPtLoc[], Face* const face,
			  Vert* const existingVert = NULL)
    {
      return insertPointOnFace(newPtLoc, face, existingVert);
    }
    virtual Vert*
    insertPointOnInternalBdryFace(const double newPtLoc[], Face* const face,
				  Vert* const existingVert = NULL)
    {
      return insertPointOnFace(newPtLoc, face, existingVert);
    }
    /*     virtual bool insertPointOnBdryEdge(const double newPtLoc[], */
    /* 				       Face* const face, Vert* const vert0, */
    /* 				       Vert* const vert1) = 0; */
  };

  class SwappingInserter2D : public Inserter {
    SwapManager2D* m_swapper;
    LengthAniso2D* m_length;
    GR_index_t m_swapsDone;
  private:
    /// Default construction is a really, really bad idea.
    SwappingInserter2D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SwappingInserter2D(SwappingInserter2D&);
    /// Operator= isn't obviously helpful, either.
    SwappingInserter2D&
    operator=(SwappingInserter2D&);
    //
  public:
    /// Constructor; the mesh must be supplied at construction time.
    SwappingInserter2D(Mesh2D* mesh, SwapManager2D* swapper,
		       LengthAniso2D* length = NULL) :
	Inserter(mesh), m_swapper(swapper), m_length(length), m_swapsDone(0)
    {
    }
    virtual
    ~SwappingInserter2D()
    {
    }

    virtual Vert*
    insertPoint(const double newPtLoc[], Cell* const cellGuess,
		Vert* const existingVert = NULL);

    /*     virtual bool insertPointInInterior(const double newPtLoc[], */
    /*                                        Cell* const cellGuess); */
    virtual Vert*
    insertPointInCell(const double newPtLoc[], Cell* const cell,
		      Vert* const existingVert = NULL);
    virtual Vert*
    insertPointOnFace(const double newPtLoc[], Face* const face,
		      Vert* const existingVert = NULL);
    /*     virtual bool insertPointOnEdge(const double newPtLoc[], Face* const face, */
    /*                                    Vert* const vert0, Vert* const vert1); */

    /*     virtual bool insertPointOnBdry(const double newPtLoc[], */
    /*                                    Face* const faceGuess); */
    virtual Vert*
    insertPointOnBdryFace(const double newPtLoc[], Face* const face,
			  Vert* const existingVert = NULL);

    virtual Vert*
    insertPointOnInternalBdryFace(const double newPtLoc[], Face* const face,
				  Vert* const existingVert = NULL);

    GR_index_t
    getNumSwapsDone() const
    {
      return m_swapsDone;
    }
    GR_index_t
    getNumInCell() const
    {
      return m_nInCell;
    }
    GR_index_t
    getNumOnEdge() const
    {
      return m_nOnEdge;
    }
    GR_index_t
    getNumOnBdryEdge() const
    {
      return m_nOnBdryEdge;
    }
    GR_index_t
    getNumOnIntBdryEdge() const
    {
      return m_nOnIntBdryEdge;
    }

  };

  class SwappingInserter3D : public Inserter {
    SwapManager3D* m_swapper;
    GR_index_t m_swapsDone;
  private:
    /// Default construction is a really, really bad idea.
    SwappingInserter3D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SwappingInserter3D(SwappingInserter3D&);
    /// Operator= isn't obviously helpful, either.
    SwappingInserter3D&
    operator=(SwappingInserter3D&);
    //
  public:
    /// Constructor; the mesh must be supplied at construction time.
    SwappingInserter3D(VolMesh* mesh, SwapManager3D* swapper) :
	Inserter(mesh), m_swapper(swapper), m_swapsDone(0)
    {
    }
    virtual
    ~SwappingInserter3D()
    {
    }

    SwapManager3D*
    changeSwapManager(SwapManager3D* newSwapper)
    {
      SwapManager3D* retVal = m_swapper;
      m_swapper = newSwapper;
      return retVal;
    }
    virtual Vert*
    insertPoint(const double newPtLoc[], Cell* const cellGuess,
		Vert* const existingVert = NULL);

    /*     virtual bool insertPointInInterior(const double newPtLoc[], */
    /*                                        Cell* const cellGuess); */
    virtual Vert*
    insertPointInCell(const double newPtLoc[], Cell* const cell,
		      Vert* const existingVert = NULL);
    virtual Vert*
    insertPointOnFace(const double newPtLoc[], Face* const face,
		      Vert* const existingVert = NULL);
    virtual Vert*
    insertPointOnEdge(const double newPtLoc[], Vert* const vert0,
		      Vert* const vert1, Vert* const existingVert = NULL);

    /// Lets test this...
    Vert*
    insertPointOnSubseg(const double newPtLoc[], Vert* const vert0,
			Vert* const vert1, Vert* const existingVert = NULL);
    // an experimental function, made to handle the case where watson fails
    // but the point to insert is "on a face" but on a curved face, and
    // does not lie in either of the cells adjacent to the face
    // really
//    Vert* insertPointOnCurvedIntBFace(const double newPtLoc[], Face* const face,
//                                    Vert* const existingVert = NULL);
    /*     virtual bool insertPointOnBdry(const double newPtLoc[], */
    /*                                    Face* const faceGuess); */
//    virtual bool insertPointOnBdryFace(const double newPtLoc[],
//                                       Face* const face);
//    virtual bool insertPointOnInternalBdryFace(const double newPtLoc[],
//                                               Face* const face);
//    virtual bool insertPointOnBdryEdge(const double newPtLoc[],
//                                       Face* const face, Vert* const vert0,
//                                       Vert* const vert1);
    GR_index_t
    getNumSwapsDone() const
    {
      return m_swapsDone;
    }
    GR_index_t
    getNumInCell() const
    {
      return m_nInCell;
    }
    GR_index_t
    getNumOnEdge() const
    {
      return m_nOnEdge;
    }
    GR_index_t
    getNumOnBdryEdge() const
    {
      return m_nOnBdryEdge;
    }
    GR_index_t
    getNumOnIntBdryEdge() const
    {
      return m_nOnIntBdryEdge;
    }
    GR_index_t
    getNumOnFace() const
    {
      return m_nOnTri;
    }
    GR_index_t
    getNumOnBdryFace() const
    {
      return m_nOnBdryTri;
    }
    GR_index_t
    getNumOnIntBdryFace() const
    {
      return m_nOnIntBdryTri;
    }

  };

/// This class stores all the data needed for inserting a point into a known
/// cavity in the mesh.
///
/// Cavity insertion removes all cells within the cavity, then "spraypaints"
/// the hull faces of that cavity with new cells, all connected to a given
/// vertex.
///
/// In the interior, this is *very* easy.
///
/// On boundaries, it's harder, because of the bdry face split.  Some old bdry
/// faces are removed and replaced with others.  This is essentially the same
/// process at one dimension lower, but that's not as well / generically
/// supported by GRUMMP data structures.
///
/// This functionality has obvious usage in Watson insertion, but it's also
/// useful when removing a vertex from the mesh by edge collapse, which
/// is why this functionality and the Watson cavity formation are segregated
/// into separate classes.
///
/// The CavityInserter also keeps track of faces that should be removed,
/// but these are auto-deleted when cells are deleted.
///
/// Responsibility for the semantics of Ruppert / Shewchuk insertion are
/// divided into several phases.
///
/// First, WatsonInfo computes the insertion hull for a candidate point.
///
/// Second, the Watson insertion driver handles bdry encroachment, and
/// may decide to abort an interior or bdry face insertion in favor of
/// bdry face  or bdry segment insertion.  For some bdry
/// protection schemes, bdry insertion may require removal of some
/// interior vertices, which happens at this stage as well.  Also, especially
/// for multicore implementation, some filtering of points that are too
/// close together may happen at this stage.
///
/// Finally, the WatsonInfo objects for points that are actually being
/// inserted are passed to a CavityInserter for action.
  class CavityInserter {
  protected:
    typedef std::set<Cell*> cellContainer_t;
    typedef std::set<Face*> faceContainer_t;
    typedef std::set<BFace*> bfaceContainer_t;

    typedef cellContainer_t::iterator cellIter_t;
    typedef faceContainer_t::iterator faceIter_t;
    typedef bfaceContainer_t::iterator bfaceIter_t;

    /// Cells in the hull that will be deleted on insertion.
    cellContainer_t m_cellsToRemove;

    /// Bdry faces on the part of the bdry being remeshed
    bfaceContainer_t m_bfacesToRemove;

    /// Faces on the insertion hull will have new cells attached to them.
    faceContainer_t m_facesOnHull;

    /// Faces inside the hull are deleted on insertion.
    faceContainer_t m_hullInteriorFaces;

    /// The mesh into which points will be inserted.
    Mesh* m_mesh;

    /// Need to know which region the new cell associated with a face
    /// should be in.
    std::vector<std::pair<Face*, int> > m_hullFaceData;
  public:
    CavityInserter(Mesh* const mesh);
    virtual
    ~CavityInserter()
    {
    }

    void
    clear()
    {
      m_cellsToRemove.clear();
      m_bfacesToRemove.clear();
      m_facesOnHull.clear();
      m_hullInteriorFaces.clear();
    }

    Mesh*
    getMesh()
    {
      return m_mesh;
    }
    void
    writeCavityVTK(const char strBaseFileName[]);
    void
    writeHullVTK(const char strBaseFileName[]);
    void
    makeCavityConvex(const double adPoint[]);
    bool
    isHullValid(const double adPoint[]);
    double
    dMinLengthOfNewFaces(double vert[]);
  protected:
    /// Record info about the regions of cells incident onto hull faces, which
    /// is going to be needed when it's time to reconnect the mesh later.
    ///
    /// Data is stored in a vector of Face*-int pairs: the int tells what
    /// region the new cell incident on the Face will be in.
    ///
    /// Pre-condition: The collection of hull faces must be up to date.
    /// Post-condition: The regions for all new cells will be recorded.
    virtual void
    recordRegions();
  private:
    /// Delete all entities inside the hull.
    ///
    /// Cells are deleted directly.  Faces are auto-deleted if/when they no
    /// longer have any incident cells.  (This is handled automatically by the
    /// mesh database.)
    ///
    /// Pre-condition:  The collection of cells in the cavity is complete.
    /// Post-condition:  Cavity is empty, with no entities inside the given
    /// hull.

    // Sets the length-scale from interpolated, if possible
    void
    deleteHullInterior() const;

    /// Reconnect all the faces in the hull to a given vertex.
    ///
    /// For every face in the hull that doesn't contain vert, a new cell is
    /// created connecting them.  Other entities inside the cavity are
    /// auto-created.
    ///
    /// Pre-condition: The cavity must be empty, and all its faces must be
    /// visible from vert and non-degenerate.
    /// Post-condition: The cavity is filled with fully-connected cells with
    /// positive volumes.  All faces of these cells are also properly connected.
    virtual void
    connectHullToVert(Vert* const vert);

  protected:
    /// Remesh parts of the bdry that must be changed during bdry insertion.
    /// This is virtual because it's done differently for topologically 2D and
    /// 3D meshes.
    virtual void
    updateMeshBdry(Vert* const newVert, const double splitParam) const = 0;
  public:
    /// Retriangulate/retetrahedralize a cavity with all faces connected to a
    /// new vertex created at the given location.  This form is especially
    /// useful as the execution phase of Watson vertex insertion.  The hull
    /// must be initialized in advance.
    ///
    /// This routine creates a vertex, then calls the vertex version of this
    /// overloaded function.
    ///
    /// Arguments:
    ///   candPtLoc[in]   The location of the new vert to which all faces on
    ///                     the hull of the cavity will be connected.
    ///   vertType[in]    Tells what flavor of vert this is (interior, bdry,
    ///                     bdry curve, etc).
    ///   splitParam[in]  For bdry faces, gives the parametric location of the
    ///                     new vertex so that the new bdry faces can be set
    ///                     properly.  For interior faces, this arg is never
    ///                     used, so any value will do.
    ///
    /// Pre-condition:  candPtLoc must be in or on the cavity and all faces of the
    ///   hull must be visible from it, except faces on the bdry that will be
    ///   deleted and replaced.
    /// Post-condition: The cavity is filled with fully-connected cells with
    /// positive volumes.  All faces of these cells are also properly connected.
    Vert*
    insertPointInHull(const double candPtLoc[], const Vert::VertType vertType,
		      const double splitParam = -1);

    /// Retriangulate/retetrahedralize a cavity with all (non-incident) faces
    /// connected to a particular vert.  This form is especially useful for
    /// vertex removal.
    ///
    /// Arguments:
    ///   targetVert[in]  The vertex to which faces of the hull will be
    ///      connected.
    ///
    /// Pre-condition / post-condition: same as point version of this overloaded
    ///   function.
    Vert*
    insertPointInHull(Vert* targetVert, const Vert::VertType vertType =
			  Vert::eInterior,
		      const double splitParam = -1);

    /// Add a cell to the cavity.  In the process, maintain the collections of
    /// interior and hull skin faces.
    ///
    /// If a face of the cell is currently listed as on the skin (because the
    /// opposite cell is in the cavity, this call will remove it from the skin
    /// and add it to the interior.  Otherwise, the face is added to the skin.
    ///
    /// Arguments:
    ///   cellToAdd[in]  The cell to be added to the cavity.
    ///
    /// Pre-condition: None.
    /// Post-condition:  cellToAdd has been added to the list of cells inside
    ///   the cavity, and the hull and interior face lists are up to date.
    /// Returns true if and only if the cell wasn't already part of the cavity.
    bool
    addCellToCavity(Cell* cellToAdd);
    bool
    isCellInCavity(Cell* const cell)
    {
      return m_cellsToRemove.count(cell);
    }

    /// Add a bdry face to the hull of the cavity; these bfaces are going to be
    /// deleted and replaced.
    ///
    /// Arguments:
    ///   bfaceToAdd[in]  The BFace to be added to the cavity.
    ///
    /// Pre-condition: None.
    /// Post-condition:  bfaceToAdd has been added to the list of BFaces on
    /// the hull of the cavity that must be replaced.  No changes to hull and
    /// interior face lists.
    void
    addBFaceToCavity(BFace* bfaceToAdd);
    bool
    isBFaceInCavity(BFace* const bface)
    {
      return m_bfacesToRemove.count(bface);
    }
    void
    fixIntBFacesInCavity();
    /// gets number of cells to remove -> used in 3D check
    int
    getNumCellsToRemove()
    {
      return m_cellsToRemove.size();
    }
  };

  class CavityInserter2D : public CavityInserter {
  protected:
    /// Remesh parts of the bdry that must be changed during bdry insertion. (2D version)
    ///
    /// In 2D, only one bdry face can be the one on which the vertex is
    /// being inserted.  Two new bdry faces are created: one connected
    /// to each end point of that edge.  This works for multi-region interior
    /// bdry cases as well.
    ///
    /// Pre-condition: Exactly one bdry face to replace.  The interior of the
    ///   cavity has already been properly re-connected.
    /// Post-condition: Old bface has been deleted, two new bfaces take its
    ///   place.
    virtual void
    updateMeshBdry(Vert* const newVert, const double splitParam) const;
  public:
    CavityInserter2D(Mesh* const mesh) :
	CavityInserter(mesh)
    {
    }
  };

  class CavityInserter3D : public CavityInserter {
  protected:
    /// Remesh parts of the bdry that must be changed during bdry insertion.
    /// This is virtual because it's done differently for topologically 2D and
    /// 3D meshes.
    ///
    /// In 3D, we can replace an arbitrary number of bdry faces N.  Typically
    /// these will be coplanar and on a single bdry patch.  However, insertion
    /// on a bdry curve (that is, adding a point on an edge that separates two
    /// parts of the surface) is supported, too, as is the case of curved bdrys.
    /// Curved bdrys are a geometric flourish on what is, here, a strictly
    /// topological process.  The N old bfaces are replaced by N+2 new ones.
    /// New bdry tris will fall on the same piece of the geometric bdry as the
    /// existing edge they are connected to.  Old bfaces are deleted.
    /// This works for multi-region interior bdry cases as well.
    ///
    /// Pre-condition: N bdry face to replace.  The interior of the
    ///   cavity has already been properly re-connected.
    /// Post-condition: Old bfaces have been deleted, new bfaces take their
    ///   place.
    virtual void
    updateMeshBdry(Vert* const newVert, const double /*splitParam*/) const;
  public:
    CavityInserter3D(Mesh* const mesh) :
	CavityInserter(mesh)
    {
    }
  };

/// This class stores all the data needed for Watson insertion of a point.
///
/// Watson insertion uses the CavityInserter for actual mesh modification.
///
/// When a point is inserted near or on the mesh boundary, the point may fall
/// inside the protection zone for a face on the boundary ("encroach" on that
/// bdry face).  In this case, the boundary faces must be split before the
/// interior point can be inserted; for bdry points that encroach bdry faces
/// that aren't being replaced, those encroached bdry faces should be split as
/// well, but this should happen after inserting the current point, to avoid
/// infinite loops.
///
/// WatsonInfo also keeps track of faces that should be removed, but these are
/// auto-deleted when cells are deleted.
///
/// Responsibility for the semantics of Ruppert / Shewchuk insertion are
/// divided into several phases.
///
/// First, WatsonInfo computes the insertion hull for a candidate point.
///
/// Second, the Watson insertion driver handles bdry encroachment, and
/// may decide to abort an interior or bdry face insertion in favor of
/// bdry face (interior only) or bdry segment insertion.  For some bdry
/// protection schemes, bdry insertion may require removal of some
/// interior vertices, which happens at this stage as well.  Also, especially
/// for multicore implementation, some filtering of points that are too
/// close together may happen at this stage.
///
/// Finally, the CavityInserter is called to do the actual insertion.

  class WatsonInserter {
  protected:
    /// The bulk of information about mesh entities to be modified is stored in
    /// in the CavityInserter.
    CavityInserter *m_CI;

    /// All bdry faces for which the new point lies in their protection zone.
    std::set<BFace*> m_bfacesEncroached;

    /// The location of the current candidate point for insertion.
    double m_candPtLoc[3];

    /// The parametric location of the candidate point.
    double m_splitParam;

    /// The type of vertex (interior, bdry, bdry curve, etc)
    Vert::VertType m_vertType;

    /// How (and whether) to check for bdry encroachment.
    eEncroachType m_encroachType;

    /// Retrieve the mesh pointer, which is stored by the cavity inserter.
    Mesh*
    getMesh() const
    {
      return m_CI->getMesh();
    }

  public:
    /// This constructor deliberately contains no info about the point to
    /// insert.  Instead, a new point can be specified via a separate member
    /// function and the object re-used for multiple insertions.
    WatsonInserter(Mesh* mesh) :
	m_CI(NULL), m_bfacesEncroached(), m_splitParam(-1.), m_vertType(
	    Vert::eUnknown), m_encroachType(mesh->getEncroachmentType())
    {
      if (mesh->getType() == Mesh::eVolMesh) {
	m_CI = new CavityInserter3D(mesh);
      }
      else {
	m_CI = new CavityInserter2D(mesh);
      }
      m_candPtLoc[0] = m_candPtLoc[1] = m_candPtLoc[2] = 0;

    }

    virtual
    ~WatsonInserter()
    {
      delete m_CI;
    }

    /// Inserts the point into its hull.  This call assumes that any
    /// encroachment fixes or hull repair is done elsewhere.  The point
    /// is specified through a previous call to computeHull.
    Vert*
    insertPointInHull(Vert* newVert = NULL);

    /// Gathers all the necessary information about a potential insertion hull.
    ///
    /// Starting from a seed cell, which must contain the vertex in its circum-
    /// ball, find all cells that contain the vertex in their circumballs.
    void
    computeHull(const double adPoint[], Cell* const seedCell,
		BFace* const onBdryFace = NULL, const double splitParam = -1);

    double
    minLengthOfNewFaces(double newVert[]);
    //////////////////////
    void
    FrontalComputeHull(const double adPoint[], Cell* const seedCell,
		       BFace* const onBdryFace = NULL, const double splitParam =
			   -1);

    /// Are there any encroached bdry faces for this Watson insertion?
    bool
    areBFacesEncroached() const
    {
      return !(m_bfacesEncroached.empty());
    }

    /// Get the encroached bdry faces
    ///
    /// Data is returned in the form of a std::pair of iterators.
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator>
    getEncroachedBFaces() const
    {
      return std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator>(
	  m_bfacesEncroached.begin(), m_bfacesEncroached.end());
    }
    ///
    int
    getNumCellsToRemove()
    {
      return m_CI->getNumCellsToRemove();
    }
    void
    writeCavityVTK(const char strBaseFileName[])
    {
      m_CI->writeCavityVTK(strBaseFileName);
    }
    void
    writeHullVTK(const char strBaseFileName[])
    {
      m_CI->writeHullVTK(strBaseFileName);
    }
    bool
    isHullValid();

  protected:

    /// Find a cell that circumscribes a given location by checking all
    /// cells in the mesh until it finds one.
    SimplexCell*
    findSeedNaively(const double newLoc[]) const;

    /// Find a cell that circumscribes a given location.
    ///
    /// First, checks the cells in the given set.  If none of them
    /// circumscribes the location, then call findSeedNaively to check
    /// every cell in the mesh.
    SimplexCell*
    findSeedFromGuesses(const double newLoc[],
			const std::set<SimplexCell*>* const guesses) const;
  };
// Used for 3D surfaces, has some things different
  class WatsonInserter3D : public WatsonInserter {

    /// All bdry faces (3D only) which are touched during the hull computation.
    std::set<BFace*> m_bfacesHit;

    /// All Subsegs (3D only) for which the new point lies in their protection zone.
    std::set<Subseg*> m_subsegsEncroached;

  public:
    WatsonInserter3D(Mesh* mesh) :
	WatsonInserter(mesh)
    {
    }
    ;

    virtual
    ~WatsonInserter3D()
    {
    }
    ;

    /// for now these are two distinct functions, distinct from the 2D calls
    /// Major difference is the presence of subsegs, which mean more checks
    /// need to be made. Could be merged with 2D, but would heavily weigh
    /// 2D code down...
    void
    compute3DHull(const double adPoint[], Cell* const seedCell,
		  BFace* const onBdryFace = NULL);
    void
    computeSubsegHull(const double adPoint[], std::vector<Cell*>& seed_cells,
		      std::vector<BFace*>& seed_bfaces);

    /// Are there any encroached subsegs for this Watson insertion?
    bool
    areSubsegsEncroached() const
    {
      return !(m_subsegsEncroached.empty());
    }

    /// Get the encroached subsegs
    std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator>
    getEncroachedSubsegs() const
    {

      return std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator>(
	  m_subsegsEncroached.begin(), m_subsegsEncroached.end());
    }
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator>
    getHitBFaces() const
    {
      return std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator>(
	  m_bfacesHit.begin(), m_bfacesHit.end());
    }
    bool
    areBFacesHit() const
    {
      return !(m_bfacesHit.empty());
    }
    /// Supplementary function
    /// Finds encroached subsegments and returns whether it found one or not
    /// stores in m_subsegsEncroached
    /// Copied Entirely from TetMeshRefiner
    bool
    hasEncroachedSubseg(Cell* const cell, const double newLoc[],
			std::set<TetRefinerEdge>* edges_visited);
  private:

  };

/// This is meant to choose point insertion locations,
/// based on a cell or boundary face. Should allow for
/// other methods to be used. At this point, finding the point
/// can be done either by using the cell, or intelligently traversing
/// based on its information
  class InsertionPointCalculator {
  public:
    // For Optimal (Erten and Ungor) we have Type I is OFFCENTER
    // Type IV is circumcenter, Type II is OTHERCIRCUMCENTER, Type III is Sink
    // such that we have [IV, I, II, III] as far as types go

    enum ePointType {
      eCircumcenter,
      eOffcenter,
      eSink,
      eOthercircumcenter,
      eShapeOptimal,
      eSizeOptimal,
      eAF,
      eCentroid,
      eAnisoEdge,
      eInvalid
    };

  protected:
    // these two are used for storing statistics
    // for example, in offcenter, stores whether the circumcenter or offcenter
    // is calculated, small overhead for just one type, but worth it for me
    GR_index_t m_numTypes;
    GR_index_t * m_numCalcs;

#ifdef HAVE_MESQUITE
    GRUMMP::TMOPQual::eQMetric Iso_QM;
    GRUMMP::TMOPQual::eQMetric Aniso_QM;
    double dMetricTol;
#endif
    // name
    std::string m_name;
  private:
    InsertionPointCalculator(const InsertionPointCalculator& IPC);
    InsertionPointCalculator&
    operator=(const InsertionPointCalculator&);
  public:
    InsertionPointCalculator(const std::string name, GR_index_t numTypes) :
	m_numTypes(numTypes), m_numCalcs(new GR_index_t[numTypes]), m_name(name)
    {
      for (GR_index_t i = 0; i < m_numTypes; i++)
	m_numCalcs[i] = 0;
#ifdef HAVE_MESQUITE
      Iso_QM = GRUMMP::TMOPQual::eShapeOrient;
      Aniso_QM = GRUMMP::TMOPQual::eShapeOrient;
      dMetricTol = 1e-6;
#endif
    }
    virtual
    ~InsertionPointCalculator()
    {
      delete[] m_numCalcs;
    }
    // Given a cell, compute the insertion point.
    // return the type of point
    virtual ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const = 0;

    // Only for aniso edge adaptment so far...
    virtual ePointType
    calcInsertionPoint(const Face * const pF, double adPointLoc[]) const
    {
      pF->doFullCheck();
      adPointLoc[0] = 0;
      assert(0);
      return eInvalid;
    }
    ;
    virtual ePointType
    calcInsertionPoint(const BFace * const pBF, double adPointLoc[]) const
    {
      pBF->doFullCheck();
      adPointLoc[0] = 0;
      assert(0);
      return eInvalid;
    }
    ;
    virtual ePointType
    calcInsertionPoint(const Vert* const pV0, const Vert* const pV1,
		       double adPointLoc[]) const
    {
      pV0->isValid();
      pV1->isValid();
      adPointLoc[0] = 0;
      assert(0);
      return eInvalid;
    }
    ;

    // Tell me what stats you want
    virtual void
    writeStats() = 0;

    void
    updateStats(ePointType pointType)
    {
      m_numCalcs[pointType]++;
    }
    const std::string&
    getName() const
    {
      return m_name;
    }
#ifdef HAVE_MESQUITE
    void
    SetIsoQM(GRUMMP::TMOPQual::eQMetric QM)
    {
      Iso_QM = QM;
    }
    void
    SetAnisoQM(GRUMMP::TMOPQual::eQMetric QM)
    {
      Aniso_QM = QM;
    }
    GRUMMP::TMOPQual::eQMetric
    getIsoQM()
    {
      return Iso_QM;
    }
    GRUMMP::TMOPQual::eQMetric
    getAnisoQM()
    {
      return Aniso_QM;
    }
    void
    SetMetricTol(double dMT)
    {
      dMetricTol = dMT;
    }
    double
    SetMetricTol()
    {
      return dMetricTol;
    }
#endif

  };
  class CircumcenterIPC : public InsertionPointCalculator {
  private:
  public:
    CircumcenterIPC() :
	InsertionPointCalculator("circumcenter", 7)
    {
    }
    ~CircumcenterIPC()
    {
      writeStats();
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      if (m_numCalcs[eCircumcenter] > 0)
	logMessage(1, "Calculated %u Circumcenters using %s\n",
		   m_numCalcs[eCircumcenter], m_name.c_str());
    }
  };

  class OffcenterIPC2D : public InsertionPointCalculator {
  private:
    // alpha is minimum angle
    double m_alpha;
    // store the coefficient used for computing offcenters/petals
    // sqrt((0.5/(1-cos(m_alpha))-0.25)
    double m_coeff;

  public:
    OffcenterIPC2D(double alpha = 0.447832396928932) :
	InsertionPointCalculator("offcenter", 7), m_alpha(alpha)
    {
      // note, triangle uses 0.475, which is smaller than mine, so their offcenter is smaller
      m_coeff = 0.5 * sqrt(1. / sin(alpha / 2.) / sin(alpha / 2.) - 1.0);
    }
    ~OffcenterIPC2D()
    {
      writeStats();
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(1, "Calculated %u Circumcenters using %s\n",
		 m_numCalcs[eCircumcenter], m_name.c_str());
      logMessage(1, "Calculated %u Offcenters using %s\n",
		 m_numCalcs[eOffcenter], m_name.c_str());
    }
  };

#ifdef HAVE_MESQUITE
  /*
   *  Class for Aniso Edge Adaptation
   *
   *  Computes the location along an edge where a point should be inserted
   *  based on TMOP quality, a metric must be previously set.
   */

  class AnisoEdgeIPC2D : public InsertionPointCalculator {

  protected:

    double dMinLength;
    double dMaxLength;

  public:
    AnisoEdgeIPC2D(double dMinML, double dMaxML) :
	InsertionPointCalculator("AnisoEdge", 5), dMinLength(dMinML), dMaxLength(
	    dMaxML)
    {
    }

    ~AnisoEdgeIPC2D()
    {
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;
    ePointType
    calcInsertionPoint(const Face * const pF, double adPointLoc[]) const;
    ePointType
    calcInsertionPoint(const BFace * const pBF, double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(4, "Calculated %u splits using %s\n", m_numCalcs[eAnisoEdge],
		 m_name.c_str());
    }
  };

  /*
   *  Class for Aniso Edge Adaptation in 3D
   *
   *  This class is intended to compute the location for inserting a point used by the anis edge adaptation in 3D.
   *  Due to the complexity of 3D, at the moment it only computes midsection of edges.
   */
  class AnisoEdgeIPC3D : public InsertionPointCalculator {

  public:
    AnisoEdgeIPC3D() :
	InsertionPointCalculator("AnisoTarget", 5)
    {
    }

    ~AnisoEdgeIPC3D()
    {
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;
    ePointType
    calcInsertionPoint(const Vert* const pV0, const Vert* const pV1,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(4, "Calculated %u splits using %s\n", m_numCalcs[eAnisoEdge],
		 m_name.c_str());
    }
  };
#endif

  class OptimalSteinerIPC2D : public InsertionPointCalculator {
  private:
    struct VoronoiNode {
      VoronoiNode(double node[2])
      {
	m_node[0] = node[0];
	m_node[1] = node[1];
      }
      double m_node[2];
      double&
      operator[](GR_index_t i)
      {
	return m_node[i];
      }
    };
    // alpha is minimum acceptable angle
    double m_alpha;
    // these are all coefficients used frequently
    double m_petalRad2;
    double m_coeffOff;
    double m_coeffRad;

  public:
    OptimalSteinerIPC2D(double alpha = 0.447832396928932) :
	InsertionPointCalculator("optimal", 7), m_alpha(alpha)
    { //, m_offcenter(new OffcenterIPC2D(alpha)) {
      m_petalRad2 = 0.25 / sin(alpha) / sin(alpha);
      m_coeffRad = sqrt(m_petalRad2 - 0.25);
      m_coeffOff = 0.5 * sqrt(1 / sin(alpha / 2) / sin(alpha / 2) - 1.0);
      assert(m_coeffRad < m_coeffOff);
    }
    ~OptimalSteinerIPC2D()
    {
      writeStats();
    }  	//delete m_offcenter;}

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(1, "Calculated %u Offcenters using %s\n",
		 m_numCalcs[eOffcenter], m_name.c_str());
      logMessage(1, "Calculated %u Sinks using %s\n", m_numCalcs[eSink],
		 m_name.c_str());
      logMessage(1, "Calculated %u OtherCircumcenters using %s\n",
		 m_numCalcs[eOthercircumcenter], m_name.c_str());
      logMessage(1, "Calculated %u Circumcenters using %s\n",
		 m_numCalcs[eCircumcenter], m_name.c_str());
    }
  };

////
////
  class OffcenterEngwirda2D : public InsertionPointCalculator {
  private:
    // alpha is minimum angle
    double m_alpha;
    // store the coefficient used for computing offcenters/petals
    // sqrt((0.5/(1-cos(m_alpha))-0.25)
    double m_coeff;
    double m_g;

    double
    FindAltitude(const Vert* pVA, const Vert* pVB, double adPointLoc[],
		 const Vert* pVC) const;

  public:
    OffcenterEngwirda2D(double grade, double alpha = 0.447832396928932) :
	InsertionPointCalculator("EngwirdaOffcenter", 7), m_alpha(alpha)
    {
      /// note, triangle uses 0.475, which is smaller than mine, so their offcenter is smaller
      m_coeff = 0.5 * sqrt(1. / sin(alpha / 2.) / sin(alpha / 2.) - 1.0);
      /// Why not just use this?
      //m_coeff = 0.5*cos(alpha/2.)/sin(alpha/2.);
      m_g = grade;

    }
    ~OffcenterEngwirda2D()
    {
      writeStats();
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(1, "Calculated %u Circumcenters using %s\n",
		 m_numCalcs[eCircumcenter], m_name.c_str());
      logMessage(1, "Calculated %u Shape-Optimal offcenters using %s\n",
		 m_numCalcs[eShapeOptimal], m_name.c_str());
      logMessage(1, "Calculated %u Size-Optimal offcenters using %s\n",
		 m_numCalcs[eSizeOptimal], m_name.c_str());
    }
  };
///////
///////
  class AF2D : public InsertionPointCalculator {
  private:
    // alpha is minimum angle
    double m_alpha;
    // store the coefficient used for computing offcenters/petals
    // sqrt((0.5/(1-cos(m_alpha))-0.25)
    double m_coeff;
    double m_g;

  public:
    AF2D(double grade, double alpha = 0.447832396928932) :
	InsertionPointCalculator("AdvancingFront", 7), m_alpha(alpha)
    {
      /// note, triangle uses 0.475, which is smaller than mine, so their offcenter is smaller
      m_coeff = 0.5 * sqrt(1. / sin(alpha / 2.) / sin(alpha / 2.) - 1.0);
      /// Why not just use this?
      //m_coeff = 0.5*cos(alpha/2.)/sin(alpha/2.);
      m_g = grade;

    }
    ~AF2D()
    {
      writeStats();
    }

    using InsertionPointCalculator::calcInsertionPoint;
    ePointType
    calcInsertionPoint(const SimplexCell * const cell,
		       double adPointLoc[]) const;

    void
    writeStats()
    {
      logMessage(1, "Calculated %u Circumcenters using %s\n",
		 m_numCalcs[eCircumcenter], m_name.c_str());
      logMessage(1, "Calculated %u AdvancingFront points using %s\n",
		 m_numCalcs[eAF], m_name.c_str());
    }
  };

///
///
///
///

  class RefinementManager2D : public Observer {
    /// Quality measure to use in prioritizing insertions
    GRUMMP::QualMeasure *m_QM;
    /// Mesh to refine
    Mesh2D* m_mesh;
    // TODO Instead of a Watson inserter, this should just be an inserter
    // passed in as an arg, so that, for instance, insert with aniso swapping
    // Just Works.
    /// The WatsonInserter does the heavy lifting when the time comes.
    WatsonInserter* m_WI;
    // a length scale sort of object, so that cell size calculations
    // can be done, too.
    std::shared_ptr<GRUMMP::Length> m_sizingField;

    // This will be useful both for aniso refinement and for things like
    // offcenter or enhanced offcenter insertion.
    InsertionPointCalculator *m_IPC;
    //
    //OffcenterCalculator *m_OC;

    /// A priority queue of cells that should be split by insertion
    InsertionQueue m_IQ;
    /// A (non-priority) queue of bdry faces that should be split by insertion
    std::deque<InsertionQueueEntry> m_bdryFacesToSplit;

    /// When refining for size, queue boundary faces for size only
    bool m_refineBdryForSize;
    int keybdryconditionNum;

    ///*********** Parallel variables *******
    InsertionQueue Remaining_m_IQ;
    std::deque<BdryEdgeBase*> Remaining_m_bdryFacesToSplit;
    std::set<Rect*> IndependentSet;
    int FilterLevel;
    int Granularity;
    int Filter;
    int ThreadID;
    GRUMMP::Quadtree *QuadT;
    bool IsInParallel;
    InsertionQueue Para_m_IQ[NUM_PROCS];
    Rect* ThreadLeaves[NUM_PROCS];
    WatsonInserter* Para_m_WI[NUM_PROCS];
//	Length2D* Para_m_sizing_field[NUM_PROCS];
//    std::vector<Vert*>VertsToSetLengthScale[NUM_PROCS];
//	std::ofstream SendEventsFile;

    double ThreadTime[NUM_PROCS][101];
    int ThreadInNum[NUM_PROCS][101];
    double TotalRefineTime;
    double JustParallelTime;
    double QueueingTime;
    double SendEventsTime;
    double SendEventsDCTime;
    double SendEventsCCTime;
    double SendEventsCBTime;
    double SecondPassesTime;
    double FIndSetsTime;
    double SyncEntryTime;
    double SyncExitTime;
    bool OkToGo;
    int Queueonface;
    std::deque<BdryEdgeBase*> Para_m_bdryFacesToSplit[NUM_PROCS];
    // These containers are introduced to store all the entities ...
    // that have been inherited from other places and at the time of ...
    // sending events, each thread could not put it into its working leaves.
//	std::vector<TriCell*> Remaining_Cells;

    /// Default construction is a really, really bad idea.
    RefinementManager2D();
    /// Copy construction and op= aren't much better.
    RefinementManager2D(RefinementManager2D&);
    RefinementManager2D&
    operator=(RefinementManager2D&);
    //Deletes the INTERIOR vertices located in the boundary edge's
    //diametral ball from the mesh.
    void
    deleteVertsInEdgeBall(BdryEdgeBase* const edge_to_split);

  public:
    /// Set the quality measure (normalized to 0..1) and mesh at construction time.
    /// Default insertion set to Circumcenters
    RefinementManager2D(Mesh2D* meshToRefine,
			GRUMMP::InsertionPointCalculator * pIPC,
			std::shared_ptr<GRUMMP::Length> sizing_field = NULL,
			QualMeasure* pQM = NULL) :
	m_QM(pQM), m_mesh(meshToRefine), m_WI(new WatsonInserter(m_mesh)), m_sizingField(
	    sizing_field), m_IPC(pIPC), m_IQ(
	    m_mesh, InsertionQueueEntry::dInvalidPriority, true), m_bdryFacesToSplit(), m_refineBdryForSize(
	    false), Remaining_m_IQ(m_mesh, true,
				   InsertionQueueEntry::dInvalidPriority, true), Remaining_m_bdryFacesToSplit(), IndependentSet(), FilterLevel(
	    1000), Granularity(10000), Filter(0), ThreadID(0)
    {
      m_mesh->addObserver(
	  this,
	  Observable::cellCreated | Observable::cellDeleted
	      | Observable::bfaceCreated);
      if (strcmp(pIPC->getName().c_str(), "circumcenter") == 0) {
	m_IQ.setPriority(InsertionQueue::eWorstQuality);
      }
      else if (strcmp(pIPC->getName().c_str(), "offcenter") == 0) {
	m_IQ.setPriority(InsertionQueue::eLongestShortest);
      }
      else if (strcmp(pIPC->getName().c_str(), "optimal") == 0) {
	m_IQ.setPriority(InsertionQueue::eShortestEdge);
      }
      QuadT = new Quadtree(m_mesh);
      ///*********** Parallel Variable Initialization ****
      IsInParallel = false;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	Para_m_IQ[ID] = InsertionQueue(m_mesh, true,
				       InsertionQueueEntry::dInvalidPriority,
				       true);
	ThreadLeaves[ID] = NULL;
	Para_m_WI[ID] = new WatsonInserter(m_mesh);
//    	  Para_m_sizing_field[ID] = new Length2D(*m_sizingField);
	ThreadTime[ID][100] = 0;
	ThreadTime[ID][0] = 0;
	ThreadTime[ID][26] = 0;
	ThreadTime[ID][25] = 0;
	ThreadTime[ID][51] = 0;
	ThreadTime[ID][50] = 0;
	ThreadInNum[ID][100] = 0;
      }
//      SendEventsFile.open("SendEventsTraking.txt");
      TotalRefineTime = 0;
      JustParallelTime = 0;
      QueueingTime = 0;
      SendEventsTime = 0;
      SendEventsDCTime = 0;
      SendEventsCCTime = 0;
      SendEventsCBTime = 0;
      SecondPassesTime = 0;
      FIndSetsTime = 0;
      SyncEntryTime = 0;
      SyncExitTime = 0;
      OkToGo = false;
      Queueonface = 0;
      keybdryconditionNum = 1;
    }
    ~RefinementManager2D()
    {
      delete m_WI;
      delete QuadT;
      for (int i = 0; i < NUM_PROCS; i++) {
	delete Para_m_WI[i];
//    	  delete Para_m_sizing_field[i];
      }
//      delete[] lck;
      if (m_mesh)
	m_mesh->removeObserver(this);
//      SendEventsFile.close();
    }
    /// Pass the parameter from user to RM2D
    void
    queuetype(int Queuetype)
    {
      if (Queuetype == 0) {
	Queueonface = 0;
      }
      else {
	Queueonface = 1;
      }
    }
    void
    setDifferentQueue(int DelaunayorEngwirda)
    {
      if (DelaunayorEngwirda == 2)
	m_IQ.iDelaunayorEngwirda = 1;
      else
	m_IQ.iDelaunayorEngwirda = 0;
    }
    void
    setKeyBdryCondition(int KBC)
    {
      keybdryconditionNum = KBC;

    }
    /// Receive and handle cells added to the database
    void
    receiveDeletedCells(std::vector<Cell*>& deletedCells);
    /// Receive and handle cells deleted from the database
    void
    receiveCreatedCells(std::vector<Cell*>& createdCells);
    /// Receive and check new bfaces for encroachment
    void
    receiveCreatedBFaces(std::vector<BFace*>& createdBFaces);

    /// Refine to match required mesh quality
    ///
    /// With a Watson inserter, circumcenter point placement,
    /// shortest-edge-to-circumradius as a quality measure, and no length scale,
    /// this is Ruppert's/Shewchuk's 2D insertion scheme.
    ///
    /// Other choices lead to other insertion schemes.  Possibilities include
    /// but definitely aren't limited to:
    ///  - adding a length scale gets you the scheme of Ollivier-Gooch and Boivin
    ///  - offcenter point placement gets you Ungor's scheme
    ///  - changes to quality, point placement, and to an aniso swapping
    ///    inserter will give aniso refinement.  Pagnutti's scheme could be
    ///    framed this way; a TMOP-based scheme is also possible.
    ///
    /// Return value: number of vertices added.
    /// where all the magic is done
    GR_index_t
    refineForQuality(const GR_index_t maxMeshVerts = INT_MAX, bool ParaFlag =
			 false,
		     int NQLayer = 5);

    GR_index_t
    refineForQualityonfaces(const GR_index_t maxMeshVerts = INT_MAX,
			    bool ParaFlag = false, int NQLayer = 5);
    void
    InsertionCalculatorForAF(Face* face, TriCell* cellToSplit,
			     double insertionPoint[]) const;
    double
    FindAltitude(const Vert* pVA, const Vert* pVB, double circumcenter[],
		 const Vert* pVC) const;

    /// Both functions below use the same "refineForQuality" function, however,
    /// refineAllQueuedCells lets RefinementManager watch
    GR_index_t
    refineMesh(const GR_index_t maxMeshVerts = INT_MAX, bool ParaFlag = false,
	       int NQLayer = 5);
    GR_index_t
    refineMeshonfaces(const GR_index_t maxMeshVerts = INT_MAX, bool ParaFlag =
			  false,
		      int NQLayer = 5);
    GR_index_t
    refineAllQueuedCells(const GR_index_t maxMeshVerts = INT_MAX);

    /// only looks at boundaries
    GR_index_t
    refineBoundaryForLength(const GR_index_t maxMeshVerts = INT_MAX);
    /// This function simply does what refinForQuality does with the addition of ...
    /// some parallel features. Also, it only refines to the defined level by the Quadtree filter.
    GR_index_t
    refinForQualityFilteredbyQuadtree(int ID, int Filter, GR_index_t &RVerts,
				      const GR_index_t maxMeshVerts = INT_MAX);
    GR_index_t
    queueAllCellsThreadWise(int ID, int FilterLevel);
    inline void
    Set_Refine2D_Timer()
    {
      TotalRefineTime = 0;
      JustParallelTime = 0;
      QueueingTime = 0;
      SendEventsTime = 0;
      SendEventsDCTime = 0;
      SendEventsCCTime = 0;
      SendEventsCBTime = 0;
      SecondPassesTime = 0;
      FIndSetsTime = 0;
      SyncEntryTime = 0;
      SyncExitTime = 0;
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	ThreadTime[ID][100] = 0;
	ThreadTime[ID][0] = 0;
	ThreadTime[ID][25] = 0;
	ThreadTime[ID][26] = 0;
	ThreadTime[ID][51] = 0;
	ThreadTime[ID][50] = 0;
	ThreadInNum[ID][100] = 0;
      }
    }

    inline void
    Print_Refine2D_Timer(int n = 2)
    {
      logMessage(n, "The total time elapsed in Refining is %3.10f\n",
		 TotalRefineTime - SendEventsTime);
      logMessage(
	  n,
	  "The time elapsed in all parallel portions of the code is %3.10f\n",
	  JustParallelTime);
      logMessage(n, "The time elapsed in queueing cells into  is %3.10f\n",
		 QueueingTime);
      logMessage(n, "The time elapsed in SendingEvents  is %3.10f\n",
		 SendEventsTime);
      logMessage(n, "The time elapsed in Sending Deleted Cells  is %3.10f\n",
		 SendEventsDCTime);
      logMessage(n, "The time elapsed in Sending Created Cells  is %3.10f\n",
		 SendEventsCCTime);
      logMessage(n, "The time elapsed in Sending Created BdryCells is %3.10f\n",
		 SendEventsCBTime);
      logMessage(n,
		 "The time elapsed in Creating second layer of QT is %3.10f\n",
		 SecondPassesTime);
      logMessage(n, "The time elapsed in Finding Independent Sets is %3.10f\n",
		 FIndSetsTime);
      logMessage(n, "The time elapsed in Synchronisation at Entry is %3.10f\n",
		 SyncEntryTime);
      logMessage(n, "The time elapsed in Synchronisation at Exit is %3.10f\n",
		 SyncExitTime);
      for (int ID = 0; ID < NUM_PROCS; ID++)
	logMessage(
	    n,
	    "The time elapsed in thread %d : %3.10f seconds, inserted %d verticies.\n",
	    ID, ThreadTime[ID][100], ThreadInNum[ID][100]);
    }

    GR_index_t
    refineBoundaryForLengthGeom(const GR_index_t maxMeshVerts =
    INT_MAX);

    /// Split any bdry edges that are encroached until they aren't any more
    /// Public so that I can fix meshes without refining.
    void
    fixBoundaryEncroachment();

    /// resets the queue, used in surface insertion
    void
    emptyQueue()
    {
      m_IQ.emptyQueue();
      assert(m_IQ.iQueueLength() == 0);
    }

    void
    setQualityTarget(double edgeToRadius)
    {
      m_IQ.vSetQualityTarget(edgeToRadius);
    }
    // controls cell refinement
    void
    setRefineForShape(bool shape)
    {
      m_IQ.qSetShape(shape);
    }
    void
    setRefineForSize(bool size)
    {
      m_IQ.qSetSize(size);
    }

    void
    updateEncroachedBFaceQueue();
    bool
    Insert(double insertionPoint[3], bool success,
	   InsertionPointCalculator::ePointType pointType,
	   TriCell* const cellToSplit);
    Vert*
    InsertEntity(InsertionPointCalculator::ePointType pointType,
		 GR_index_t& vertsAdded);
    int
    BdryConditionFromFace(Face* face);

  private:
    /// Find the location of and insert a new point to split this triangle.
    ///
    /// Depending on the choice of Inserter, the resulting mesh may be Delaunay
    /// or it may have been swapped to match some other quality criterion.
    ///
    /// Returns true iff the triangle is split successfully, so it'll return
    /// false if the insertion point encroaches on a bdry edge.
    bool
    splitTriangle(TriCell* const cellToSplit, bool Para_Flag = false,
		  int ID = 0);

    /// Find the location of and insert a new point to split this bdry edge.
    ///
    /// The location choice is currently based on total variation of the
    /// tangent.
    ///
    /// Returns true iff the edge is split successfully, which had better be
    /// always.
    bool
    splitBdryEdge(BdryEdgeBase* const bdryEdgeToSplit, bool Para_Flag = false,
		  int ID = 0);

    /// Add all encroached bdry edges to the queue so that they can be split
    /// first.
    void
    queueEncroachedBdryEdges();

    /// Add all cells that need splitting to the queue.
    GR_index_t
    queueAllCells(bool Para_Flag = false, int NQLayers = 5);
    /// Add all faces at the frontal
    GR_index_t
    queueAllFaces(bool Para_Flag = false, int NQLayers = 5);

    /// Add this cell to the priority queue
    ///
    /// For iso meshes, priority will be based on size and shape.  The same will
    /// likely be true for aniso meshes, but the measures will be different
    bool
    addToPriorityQueue(Cell* const cell);
    bool
    addToFrontal(Face* const face);

  };

#ifdef HAVE_MESQUITE
  /*
   * Class for Aniso Edge Adaptation
   *
   * This class is for aniso adapting a mesh by modifying edges, based on
   * TMOP quality and metric length. A metric must be set in the mesh.
   */
  class AnisoEdgeRefinementManager2D : public Observer {

    enum eInterpType {
      eDefault, eAverage, eLinear, eLogEuclidean
    };
    enum eInsertType {
      eQualityEdge, eMetricLengthEdge, eCoarsening
    };
    TMOPQual *m_TQ;
    Mesh2D* m_mesh;
    GRUMMP::SwapDecider2D* m_pSD2D;
    GRUMMP::SwapManager2D *m_pSM2D;
    SwappingInserter2D* m_SI2D;
    GRUMMP::SmoothingManagerMesquite2D *m_pSMM;
    InsertionPointCalculator *m_IPC;
    AnisoRefinement *m_pAR;
    eInterpType m_Type = eDefault;
    InsertionQueue m_IQBF;
    InsertionQueue m_IQIF;
    int iNMetricValues;
    double m_maxMetricLength;
    double m_minMetricLength;
    bool bReceiveCreated, m_useGlobalSmoothing;

    AnisoEdgeRefinementManager2D();
    AnisoEdgeRefinementManager2D(RefinementManager2D&);
    RefinementManager2D&
    operator=(RefinementManager2D&);

  public:

    AnisoEdgeRefinementManager2D(Mesh2D* meshToRefine,
				 const double dMax_MLength,
				 const double dMin_MLength, const double global,
				 GRUMMP::InsertionPointCalculator * pIPC,
				 GRUMMP::SmoothingManagerMesquite2D *pSMM,
				 AnisoRefinement *pAR, SwapDecider2D *pSD2D =
				 NULL) :
	m_TQ(new GRUMMP::TMOPQual(2)), m_mesh(meshToRefine), m_pSD2D(pSD2D), m_pSM2D(
	    new GRUMMP::SwapManager2D(m_pSD2D, m_mesh)), m_SI2D(
	    new SwappingInserter2D(m_mesh, NULL)), m_pSMM(pSMM), m_IPC(pIPC), m_pAR(
	    pAR), m_IQBF(m_mesh, InsertionQueueEntry::dInvalidPriority, true), m_IQIF(
	    m_mesh, InsertionQueueEntry::dInvalidPriority, true), iNMetricValues(
	    3), m_maxMetricLength(dMax_MLength), m_minMetricLength(
	    dMin_MLength), bReceiveCreated(false), m_useGlobalSmoothing(global)
    {
      m_mesh->addObserver(this,
			  Observable::bfaceCreated | Observable::faceCreated);
      m_IQBF.vSetMinMLength(dMin_MLength);
      m_IQBF.vSetMaxMLength(dMax_MLength);
      m_IQBF.vSetMetricTol(1e-6);
      m_IQIF.vSetMinMLength(dMin_MLength);
      m_IQIF.vSetMaxMLength(dMax_MLength);
      m_IQIF.vSetMetricTol(1e-6);
      m_mesh->allowBdryChanges();
      m_SI2D->setForcedInsertion(true);
      m_mesh->attachAnisoRef(m_pAR);
      SetInterpType(eDefault);
    }
    ~AnisoEdgeRefinementManager2D()
    {
      delete m_SI2D;
      delete m_pSM2D;
      delete m_TQ;
      if (m_mesh)
	m_mesh->removeObserver(this);
    }

    void
    receiveCreatedFaces(std::vector<Face*>& createdFaces);
    void
    receiveCreatedBFaces(std::vector<BFace*>& createdBFaces);

//	  Runs the algorithm to aniso adapt a mesh, a metric must be previously set in the mesh
    GR_index_t
    AnisoRefineMesh();

//	  Modifies the mesh according to the queues
    GR_index_t
    AnisoRefineAllQueued();

//	  Splits edges if the quality of the configuration after inserting a point is improved
    GR_index_t
    AnisoQualityEdgeInsertion();

//	  Splits edges which metric length is bigger than a max metric length allowed
    GR_index_t
    AnisoMetricLengthEdgeInsertion();

//	  Coarses the mesh based on metric length. Edges are queued for deleting if any of its vertices has at least two edges smaller
//	  than a min metric length. Wich point to delete with edge collapse is decided by quality
    GR_index_t
    AnisoMetricEdgeCoarsening();

//	  Uses Local smoothing with mesquite based on targets defined by an eigen decomposition of a metric
//	  iIter = number of passes of local smoothing
    void
    LocalSmoothing(int iIter);

//	  Uses Global Smoothing with mesquite based on targets defined by an eigen decomposition of a metric
    void
    GlobalSmoothing();

    void
    emptyQueue(InsertionQueue IQ)
    {
      IQ.emptyQueue();
      assert(IQ.iQueueLength() == 0);
    }

  private:

//	  Decides where to split an edge, exact location is computed with InsertionPointCalculator
    bool
    splitFace(Face* const FaceToSplit);

//	  Decides where to split a boundary edge, exact location is computed with InsertionPointCalculator
    bool
    splitBdryEdge(BdryEdgeBase* bdryEdgeToSplit);

//	  Functions add to queue all the edges to be modified by priority
    GR_index_t
    queueAll();
    GR_index_t
    queueFaces();
    GR_index_t
    queueBFaces();

//	  Adds elements to the queue
    bool
    addToPriorityQueue(Face* const pF);
    bool
    addToPriorityQueue(BdryEdgeBase* const pBE);

//	  Functions used to interpolated metrics
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

//	  Sets priorities for each step of the process
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

///  FOR NOW, IN ITS OWN HEADER

//  class RefinementManager3D : public Observer {
//    /// Quality measure to use in prioritizing insertions
//    GRUMMP::QualMeasure<3> *m_QM;
//    /// Mesh to refine
//    Mesh2D* m_mesh;
//    //  Instead of a Watson inserter, this should just be an inserter
//    // passed in as an arg, so that, for instance, insert with aniso swapping
//    // Just Works.
//    /// The WatsonInserter does the heavy lifting when the time comes.
//    WatsonInserter* m_WI;
//    // Add a length scale sort of object, so that cell size calculations
//    // can be done, too.
//
//    // Eventually will want an InsertionPointCalculator (or some such name)
//    // This will be useful both for aniso refinement and for things like
//    // offcenter or enhanced offcenter insertion.
//
//    /// A priority queue of cells that should be split by insertion
//    InsertionQueue m_IQ;
//    /// A (non-priority) queue of bdry faces that should be split by insertion
//    std::deque<BdryEdgeBase*> m_bdryFacesToSplit;
//
//    /// Default construction is a really, really bad idea.
//    RefinementManager3D();
//    /// Copy construction and op= aren't much better.
//    RefinementManager3D(RefinementManager3D&);
//    RefinementManager3D& operator=(RefinementManager3D&);
//
//  public:
//    /// Set the quality measure (normalized to 0..1) and mesh at construction time.
//    RefinementManager3D(VolMesh* meshToRefine, QualMeasure<3>* pQM = NULL) :
//      m_QM(pQM), m_mesh(meshToRefine), m_WI(new WatsonInserter(m_mesh)),
//      m_IQ(m_mesh, InsertionQueueEntry::dInvalidPriority, true)
//    {
//      m_mesh->addObserver(this, Observable::cellCreated |
//                          Observable::cellDeleted);
//    }
//    ~RefinementManager3D() {
//      delete m_WI;
//      if (m_mesh) m_mesh->removeObserver(this);
//    }
//    /// Receive and handle cells added to the database
//    void receiveDeletedCells(std::vector<Cell*>& deletedCells);
//    /// Receive and handle cells deleted from the database
//    void receiveCreatedCells(std::vector<Cell*>& deletedCells);
//
//    /// Refine to match required mesh quality
//    ///
//    /// With a Watson inserter, circumcenter point placement,
//    /// shortest-edge-to-circumradius as a quality measure, and no length scale,
//    /// this is Shewchuk's 3D insertion scheme.
//    ///
//    /// Other choices lead to other insertion schemes.  Possibilities include
//    /// but definitely aren't limited to:
//    ///  - adding a length scale gets you the scheme of Ollivier-Gooch and Boivin
//    ///  - offcenter point placement gets you Ungor's scheme
//    ///  - changes to quality, point placement, and to an aniso swapping
//    ///    inserter will give aniso refinement.  Pagnutti's scheme could be
//    ///    framed this way; a TMOP-based scheme is also possible.
//    ///
//    /// Return value: number of vertices added.
//    GR_index_t refineForQuality(const GR_index_t maxMeshVerts = INT_MAX);
//
//  private:
//    /// Find the location of and insert a new point to split this triangle.
//    ///
//    /// Depending on the choice of Inserter, the resulting mesh may be Delaunay
//    /// or it may have been swapped to match some other quality criterion.
//    ///
//    /// Returns true iff the triangle is split successfully, so it'll return
//    /// false if the insertion point encroaches on a bdry edge.
//    bool splitTet(TetCell* const cellToSplit);
//
//    /// Find the location of and insert a new point to split this bdry tri.
//    ///
//    /// The location choice is currently the circumcenter of the tri.
//    ///
//    /// Returns true iff the tri is split successfully, which had better be
//    /// always, unless the split point would encroach on a bdry (sub)segment.
//    bool splitBdryFace(BdryTriBase* const bdryTriToSplit);
//
//    /// Add all encroached bdry edges to the queue so that they can be split
//    /// first.
//    void queueEncroachedBdryFaces();
//
//    /// Split any bdry edges that are encroached until they aren't any more
//    void fixBoundaryEncroachment();
//
//    /// Add all cells that need splitting to the queue.
//    GR_index_t queueAllCells();
//
//    /// Add this cell to the priority queue
//    ///
//    /// For iso meshes, priority will be based on size and shape.  The same will
//    /// likely be true for aniso meshes, but the measures will be different
//    bool addToPriorityQueue(Cell* const cell);
//  };

}
// End namespace

#endif /* GR_INSERTIONMANAGER_H_ */
