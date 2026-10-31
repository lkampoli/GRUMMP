#ifndef GR_SURFACEINSERTION_H
#define GR_SURFACEINSERTION_H 1

#include "GR_BaryCentricLengthScale.h"
#include "GR_FacetSurfIntersect.h"
#include "GR_GRGeom2D.h"
#include "GR_GRPoint.h"
#include "GR_InsertionManager.h"
#include "GR_Length.h"
#include "GR_LengthAniso2D.h"
#include "GR_RefinementManager3D.h"
#include "GR_SmoothingManager.h"
#include "GR_SubsegManager.h"
#include "GR_SurfaceSampler.h"
#include "GR_SurfIntersect.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_SurfMesh.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_TetMeshBuilder.h"

#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "DLIList.hpp"
//#include "FacetEvalTool.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "RefVertex.hpp"

/// This is the header for all Surface/Curve Insertion,
/// 2D is really a subset of 3D, so all function naming
/// based on 3D class, even though 2D written first
/// Major difference is in the sampling, as 3D uses
/// existing sampling code written by S. Gosselin

namespace GRUMMP
{

/// SurfaceInsertion Base Class
/// Base Class should never explicitly created (abstract), but handles a lot of
/// information/function calls both 2D/3D classes need
/// SurfaceInserter is an Observer so it can determine
/// what vertices are left alone in the whole process

  class SurfaceInserter : public Observer {
  private:
    /// Default construction is a really, really bad idea.
    SurfaceInserter();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SurfaceInserter(SurfaceInserter&);
    /// Operator= isn't obviously helpful, either.
    SurfaceInserter&
    operator=(SurfaceInserter&);
  protected:
    // Place to store the lengthscale
    std::shared_ptr<GRUMMP::Length> m_isoLength;

    // This is the insertionPointCalculat
    InsertionPointCalculator *m_IPC;

    /// Initial LengthScale Inputs.
    /// Free parameters (mesherments)
    /// Used for LengthScale calculations
    double m_scale;
    double m_grading;

    // Length scale and Metric management objects
    LengthAniso* m_anisoLength;

    /// If true, lengthscales are computed as the minimum
    /// of the mesh lengthscale and the radius of curvature
    /// Is the mesh anisotropic, is it structured?
    /// All the flags for mesh type

    // sampling by curvature removed, wasn't being used at all
//	bool m_sampleByCurvature;
    bool m_isAnisotropic;
    bool m_isStructured;
    bool m_useOffsets;
    bool m_offsetByNormals;

    // these are used for handling offsets
    // used for linear translational offsets
    std::pair<int, CubitVector> m_offsetCoordinate;
    // collects pairs, based on "curve number"
    std::vector<std::pair<int, int>> m_offsetPairs;
    GR_index_t m_numberOfOffsets;

    char m_params[256];
    /// Statistics to keep, all positive
    GR_index_t numInitVerts, numInsertedVerts, numRemovedVerts,
	numSmoothedVerts, numRefinedVerts, numInitCells;

    /// set of Verts we could not remove from removal
    std::set<Vert*> m_leftoverVertsToRemove;
    /// set of Verts we could not remove during insertion
    /// Or set of Verts that don't encroach. Only used in 3D for now
    std::set<Vert*> m_vertsToRemoveAfterInsertion;

    /// Queue for keeping track of unchanged verts
    /// basically follows verts, and if they are smoothed
    /// marks them as moved (using active/inactive)
    /// Used for statistic collecting
    std::set<Vert*> vertQueue;
    void
    receiveMovedVerts(std::vector<Vert*>& movedVerts);

    /// Given a cell and a region number, set all cells
    /// within that boundary to that value
    /// cell is initial cell
    /// region is region number to set to
    /// mesh included so it can be outputted everytime new region number set,
    /// really only used to make pretty pictures
    void
    setRegionFromCell(Cell * cell, const int region, Mesh * mesh = NULL) const;

    // updates the verts after purging done
    void
    updateVertsToRemoveAfterInsertion(const std::map<Vert*, Vert*>& vert_map);

  protected:

    /// straightforward constructor/destructor
    SurfaceInserter(std::shared_ptr<Length> len, double grading, double scale,
		    const char strSurfParams[] =
		    NULL);

    /// these could be public, but keep them here
    /// both collect information on the mesh
    void
    initializeStats(Mesh* mesh);
    void
    finalizeStats(Mesh & mesh);

    // outputs a lot of the meshes in debug mode, avoids cluttering

  public:
    /// Destructor is public and virtual
    virtual
    ~SurfaceInserter()
    {
      if (m_IPC)
	delete m_IPC;
      if (m_anisoLength)
	delete m_anisoLength;
    }
  };

/// stores the point parameter, the cell its in, and the vert, if its been
/// inserted into the mesh.
/// insertionPoint makes life easier in 2D, its basically a discrete curve,
/// as in a surface mesh in 2D.
/// by storing the cell the point is already in, this makes lookup way easier,
/// and traversal based
  class insPoint2D {
  private:
    // since only the parameter is used for comparison, the others can be
    // happily modified if need be, without adjusting the order of the set
    double m_param;
    mutable Cell * m_cell;
    mutable Vert * m_vert;
    mutable double m_metric[5];
  public:
    insPoint2D() :
	m_param(-1.), m_cell(NULL), m_vert(NULL)
    {
    }
    ;
    insPoint2D(double param, Cell * cell = NULL, Vert * vert = NULL) :
	m_param(param), m_cell(cell), m_vert(vert)
    {
    }
    ;

    // very basic, needs ability to get itself and such
    double
    getParam() const
    {
      return m_param;
    }
    ;
    void
    setParam(double param)
    {
      m_param = param;
    }
    Vert *
    getVert() const
    {
      return m_vert;
    }
    ;
    void
    setVert(Vert * vert) const
    {
      m_vert = vert;
    }
    ;
    Cell *
    getCell() const
    {
      return m_cell;
    }
    ;
    void
    setCell(Cell * cell) const
    {
      m_cell = cell;
    }
    ;
    double
    getMetric(int i) const
    {
      return m_metric[i];
    }
    void
    setMetric(double metric[5]) const
    {
      for (int i = 0; i < 5; i++) {
	m_metric[i] = metric[i];
      }
    }
    ;
    // this allows for sorting in sets, and such.
    // Rounds to the 13th decimal place, which is fine
    // for comparison
    bool
    operator<(const insPoint2D & cPBF) const
    {
      return round(m_param * 1.e13) / 1.e13
	  < round(cPBF.m_param * 1.e13) / 1.e13;
    }
  };
// definitions of data structure used in 2D
// Every curve has its own sorted set of points, sorted by parameter
// and the whole object has a map to the set up,
// really only set up for potentially offsetting/layering,
// otherwise, the order of curves is irrelevant
  typedef std::set<insPoint2D> insPointSet;
  typedef std::set<insPoint2D>::iterator insPointSet_it;
  typedef std::map<GRCurve*, insPointSet*> insertionCurveMap;
  typedef std::map<GRCurve*, insPointSet*>::iterator insertionCurveMap_it;

/// Here is how the input to surfaceInserter2D works

/// this is a bit.. complicated, to allow for the versatility of surface insertion
/// it works similar to the settings in the original GRUMMP
/// 'd' corresponds to the swapdecider, in this case theres only one choice for now
/// 'dd' -> Use Delaunay SwapDecider
/// 'a' corresponds to use anisotropic
/// anisotropic refinement is theoretically possible
/// but currently not working well
/// 'o' is one way to do offsets
/// it takes a given curve, and offsets it by a prescribed vector
/// a certain number of times.
/// sample: o1x0.0y0.1n100, offsets curve 1 by (0,0.1) 100 times
/// if curvenum < 0, offsets all curves by normal amount
/// o-1x0.0y0.1n10 offsets all curves in a normal distance by the magnitude
/// of the vector, will include examples
/// 'p' offsets by pairs, in the corresponding boundary file
/// projects from 'p1' to 'p2'
/// sample: o2c1c2c2c3, offsets curve 1 to curve 2, and curve 2 to curve 3
/// in a normal direction
/// lastly 'q' is used for quasi-structured meshes
/// should always finish a string with offset info

/// theres an idea that the anisotropic version and structured
/// and offset could be higher level classes on a base 2D class
/// but realistically, its just easier to use flags to direct it

/// The obvious next change is to use a quad tree to find cells for
/// point insertion
  class SurfaceInserter2D : public SurfaceInserter {
  private:
    /// Default construction is a really, really bad idea.
    SurfaceInserter2D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SurfaceInserter2D(SurfaceInserter2D&);
    /// Operator= isn't obviously helpful, either.
    SurfaceInserter2D &
    operator=(SurfaceInserter2D&);

  public:
    SurfaceInserter2D(Mesh2D* mesh, std::shared_ptr<GRUMMP::Length> len,
		      const double scale, const double grading,
		      const char strSurfParams[]);

    ~SurfaceInserter2D();

    /// does all the work
    void
    insertCurveInMesh(const char strCurveFilename[], const char strParams[]);
    void
    refineBoundaryForLength(GRGeom2D* geom);
	void equidistributeArcLength(
			std::vector<double>& paramToInsert,
			const std::vector<double>& newPointsArclength,
			GRCurveGeom* curveGeom, insPointSet& pointSet);

  private:

    /// Basically carries one of everything
    /// Mesh cleanup handled in tri.cxx
    Mesh2D * m_mesh;

    /// Geometry only used internally, handled internally
    GRGeom2D * m_geom;

    /// keeps track of points I want to insert
    insertionCurveMap m_pointsToInsert;

    // if we want to offset
    std::multimap<GRCurve *, GRCurve*> m_offsetCurveMap;

    /// these are declared here so they can observe
    /// everything that goes on, and know where to locally
    /// make changes, SwapDecider is a pointer, because of the way
    /// that class was defined, and used by other classes

    SwapDecider2D *SwapDec2D;
    SwapManager2D *SwapMan2D;
    OptMSSmoothingManager2D SmoothMan2D;

    union {
      RefinementManager2D * RefMan2D;
      AnisoRefinement * ARMesh2D;
    };
    /// all the insertion, recovery, and creation of new boundary faces is done in here
    int
    insertPointsIntoMesh();

    int
    insertBoundaryPointsIntoMesh();
    /// sampleSurface samples boundary curves and optimizes using Jacobi
    /// and an equidistribution principle
    int
    sampleSurface();
    int
    sampleSurfaceForBoundaryPoints();
    // this does the offset part of things
    void
    sampleOffsetCurves();

    /// an experimental idea
    int
    sampleQuasiStructured();

    /// Removes vertices that encroach on sampled edges and are too close to sampled points
    int
    removeVertsNearSurface();

    /// updates the containing cell for each sample point
    /// makes finding intersections with internal boundaries easier
    /// as well as point placement
    void
    updateSamplePoints();

    // post processing, based on input string
    void
    refineAndSmooth(const char strParams[]);

    // rarely ever called, but there in case you decide to go crazy and for some reason
    // get yourself into a situation where you can't remove everything you want to.
    int
    removeVertsAfterInsertion();

    /// used to relabel regions
    void
    setRegions();

    /// Determines intersections between internal boundaries and a particular
    /// curve, and stores in a set with the internal edge and the parameter number
    /// its a lot of information

    /// requires that internal boundaries have distinct regions on each side
    void
    sampleForIntBdryIntersections();

    /// checks all consecutive new internal boundaries for a small angle,
    /// then inserts additional points on the concentric circle
    /// based off of small angle protection done for TriMeshBuilder
    /// protects small angles for future refinement, often is overkill
    /// but guarantees that the refinement algorithm will terminate
    void
    identifyAndFixSmallAngle(Vert * apex_vert);

    /// used for anisotropic refinement
    /// two metrics are used, one for length scale
    /// and one to ensure we can swap back to the original anisotropy
    void
    setVertLengthScaleMetrics();
    void
    setVertOriginalMetrics();
    void
    setVertMetric(Vert * vert);

    // creates quads between offset layers
    void
    createOffsetQuadCells();

    // useful for debugging
    void
    writeSamplePointsToFile(const char * filename);
  };

/// 3D version of the class

/// Currently, the major algorithmic difference is that "refVolumes" are
/// inserted one at a time, and the whole process is done first
/// If a closed surface is inserted, it is assumed that it does not contact
/// the boundary all. If open surface is inserted,
/// it is assumed that the boundaries of the open surface
/// match up with the boundaries of the background mesh
  class SurfaceInserter3D : public SurfaceInserter {
  private:
    /// Default construction is a really, really bad idea.
    SurfaceInserter3D();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SurfaceInserter3D(SurfaceInserter3D&);
    /// Operator= isn't obviously helpful, either.
    SurfaceInserter3D &
    operator=(SurfaceInserter3D&);

  public:

    ///  only valid constructor
    SurfaceInserter3D(VolMesh * pVMesh, std::shared_ptr<Length> len,
		      double grading, double scale, const char strSurfParams[]);

    /// destructor
    ~SurfaceInserter3D();

    /// calls all the base functions in order
    void
    insertSurfaceInMesh(const char strSurface[], const char strParams[]);

  private:
    /// initial mesh
    VolMesh * m_mesh;

    /// basic idea is that each connected surface has its own surface mesh
    /// if this is an array of more than size 0, then it means we are inserting
    /// multiple disjoint surfaces in at a time
    /// surface mesh(es) used after sampling
    /// IMPORTANT DISTINCTION: surface meshes use "Cells" to make the mesh,
    /// These cells are 3D TriCells, however since they are inserted
    /// as 'faces' into the mesh, there is a naming issue
    std::vector<SurfMesh*> m_surfMeshes;

    /// holds on to the "current" refVolume we are working on at a time
    /// could probably pass into each function
    RefVolume * m_currentRefVolume;

    /// algorithm is slightly different for a closed surface than a singly connected one
    /// this is due to the surface sampling algorithm already in existence
    bool m_isSurfaceClosed;

    // same as in 2D, all operate on the initial mesh
    DelaunaySwapDecider3D* SwapDec3D;
    SwapManager3D* SwapMan3D;
    OptMSSmoothingManager3D * SmoothMan3D;
    RefinementManager3D * RefMan3D;

    std::multimap<RefFace*, RefFace*> m_offsetSurfaceMap;

    // base functions used in 3D
    // all return a number of points
    int
    sampleSurface();

    // removal is different than in 2D
    // For recovery, only vertices that encroach are removed initially
    // after insertion, vertices that would be too close are removed
    // these are also used to help recovery, giving an extra set of vertices
    // to use to help with removal
    int
    removeVertsNearSurface();
    int
    removeVertsAfterInsertion();
    int
    removeLeftoverVerts();
    // this handles the point insertion and recovery
    // recovery is the big challenge
    int
    insertPointsIntoMesh();

    // post processing, more or less identical to 2D in structure
    void
    refineAndSmooth(const char strParams[]);

    /// needed for closed surfaces, Could be eliminated with proper
    /// search structure
    BFace*
    findNearestBoundaryCell(const double * coords, double * projectedCoords);

    /// slightly different than 2D, in that the geometry data structure is different
    void
    setRegions();

    // copying surface meshes isn't allowed, nor is it trivial to implement
    // so this ensures distinct surfaces are created, off the base information
    void
    sampleOffsetSurfaces(FaceVec & fqes, std::set<Subseg*> & subsegSet);

    /// all of these are functions are additional to 3D,
    /// required to handle the extra topology not present
    /// in two dimensions

    // this overrides the normal one, and fixes subsegments
    // eventually subsegments should be better handled by the base mesh
    bool
    removeVert(Vert * vert);

    // Resample for intersections with internal boundaries
    // Works by taking the first sampling of the mesh
    // Determining the discrete intersection between
    // the surface and the internal boundary
    // using the first sampling and moving points
    // to form new subsegments on the mesh
    // the new subsegments are added to the old
    // and the sampling process is redone
    // potentially expensive, but stronger guarantees
    void
    sampleForIntBdryIntersections();

    // experimental function
    int
    sampleQuasiStructured();

    // attempts to recover face by swapping for quality
    // in the neighborhood to change the topology enough
    Face *
    recoverFaceBySwapping(Vert * pV0, Vert * pV1, Vert * pV2);

    bool
    recoverFacePairByEdgeSwap(Face * sharedFace, Vert * sharedVerts[],
			      Vert * vertsA[], Vert * vertsB[],
			      std::pair<Face*, Face*> & facePair);
    Face*
    attemptFaceRecovery(std::vector<Face*>& newFaces, const GR_index_t iC,
			Vert* const faceVerts[3],
			std::vector<Cell*>& facesUnableToRecover,
			Cell* const faceToInsert, GR_index_t& iFacesRecovered,
			Cell*& cellGuess);
  };

}

void
smoothBdryRefine(Mesh2D* pM2D, GRGeom2D* pG2D, double dBGrading, double dBScale,
		 char strSurfInsParams[]);

#endif
