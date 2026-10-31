#ifndef GR_SmoothingManager_h
#define GR_SmoothingManager_h

#include <set>

#include "GR_Cell.h"
#include "GR_Mesh.h"
#include "GR_Observer.h"
#include "GR_Vertex.h"
#include "SMsmooth.h"

#ifdef HAVE_MESQUITE
#include "Mesquite.hpp"
#include "MsqIMesh.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "AnisoWrapper.hpp"
#include "InstructionQueue.hpp"
#include "PlanarDomain.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "MeshWriter.hpp"
#include "IdealWeightInverseMeanRatio.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "iMesh_GRUMMP.hh"
#endif

namespace GRUMMP
{
  class SmoothingManager : public Observer {
  protected:
    Mesh* m_pMesh;
  private:
    /// Default construction is a really, really bad idea.
    SmoothingManager();
    /// Copy construction isn't harmful, but it isn't clear that it's
    /// helpful, either.
    SmoothingManager(SmoothingManager&);
    /// Operator= isn't obviously helpful, either.
    SmoothingManager&
    operator=(SmoothingManager&);
    // Need to have a queue and ways to access it, but that can wait for
    // specializations --- OpenMP, MPI, etc --- and/or Observer giving
    // feedback about faces that should be looked at.
    //
    // The following set is the simple serial version.  For parallel
    // (either version) it isn't going to be adequate.
    std::set<Vert*> vertQueue;
  public:
    /// Set the mesh at construction time.
    SmoothingManager(Mesh* pMeshToSmooth) :
	m_pMesh(pMeshToSmooth), vertQueue()
    {
      m_pMesh->addObserver(
	  this,
	  Observable::cellCreated | Observable::vertDeleted
	      | Observable::vertMoved);
    }
    virtual
    ~SmoothingManager();

    // Observer function overrides
    void
    receiveDeletedVerts(std::vector<Vert*>& deletedVerts);
    void
    receiveMovedVerts(std::vector<Vert*>& movedVerts);
    void
    receiveCreatedCells(std::vector<Cell*>& createdCells);
    void
    emptyQueue()
    {
      vertQueue.clear();
      assert(vertQueue.size() == 0);
    }

    virtual int
    smoothAllVerts(const int passes);
    int
    smoothAllQueuedVerts();
    virtual int
    smoothVert(Vert* vert) = 0;
  };

#ifdef HAVE_MESQUITE
  class SmoothingManagerMesquite2D : public SmoothingManager {
  protected:
    Mesquite::AnisoWrapper Queue;
    Mesquite::MsqIMesh* msqMesh;

  public:
    SmoothingManagerMesquite2D(Mesh2D * const pM2D);
    ~SmoothingManagerMesquite2D();
    void
    vSetTargetfromMetric();

//	  Smooths all the points in the mesh
    int
    smoothAllVerts(const int passes);

//	  Creates a patch the neighbouring cells of the point and passes it to Mesquite for smoothing
//	  For the local smoothing function of SmoothingManagerMesquite2D, the maxTime at the AnisWrapper.cxx
//	  for 2D smoothing must be set really low (1 or 2)

    int
    smoothVert(Vert* pV);
  };

  /*
   * Similar approach for smoothing with TMOP based on target matrices as in 2D.
   * At the moment only rigth angle tetrahedra are set as targets,
   * due to the complexity of 3D, a higher range of targets is needed for this smoothing
   * approach to work accordingly.
   */
  class SmoothingManagerMesquite3D : public SmoothingManager {
  protected:

    Mesquite::AnisoWrapper3D Queue;
    Mesquite::Mesh* msqMesh;
    bool bITarget;
  public:
    SmoothingManagerMesquite3D(VolMesh * const pM3D, bool bIsoTarget);
    void
    vSetTargetfromMetric();
    int
    smoothAllVerts(const int passes);
    int
    smoothVert(Vert* /*vert*/)
    {
      assert(0);
      return 0;
    }
  };
#endif

  class OptMSSmoothingManager : public SmoothingManager {
  protected:
    /// An opaque struct that contains context data for mesh smoothing.
    SMsmooth_data *m_smoothData;
    double **m_neighCoords;
    int **m_faceConn;
  public:
    OptMSSmoothingManager(Mesh * const pMeshToSmooth, const int iDim);
    ~OptMSSmoothingManager()
    {
      SMfinalizeSmoothing(m_smoothData);
    }
    static OptMSSmoothingManager*
    Create(Mesh * const pMeshToSmooth);
    virtual int
    smoothVert(Vert* vert) = 0;
    virtual int
    smoothAllVerts(const int passes = 5)
    {
      SMinitSmoothStats(m_smoothData);
      int retVal = SmoothingManager::smoothAllVerts(passes);
      SMprintSmoothStats(m_smoothData);
      resetThreshold();
      return retVal;
    }
    // Set a new smoothing measure
    virtual void
    setSmoothingGoal(const int goal) = 0;
    // Set a new smoothing technique (default: floating threshold, which is
    // quite good).
    void
    setSmoothingTechnique(const int technique)
    {
      SMsetSmoothTechnique(m_smoothData, technique);
    }
    void
    setFunctionID(const int FunctionID)
    {
      SMsetFunctionID(m_smoothData, FunctionID);
    }
    // Set a cutoff below which optimization will be applied.
    void
    setSmoothingThreshold(const double threshold)
    {
      SMsetSmoothThreshold(m_smoothData, threshold);
    }
    // Reset the cutoff based on the results from the last pass of smoothing.
    void
    resetThreshold()
    {
      // For floating threshold after at least one smoothing pass, this
      // will reset the threshold using the global minimum value from
      // the last pass.  In other cases, the threshold value remains
      // unchanged.
      double threshold = m_smoothData->smooth_param->lap_accept_value;
      SMconvertToDegrees(m_smoothData->smooth_param->function_id, &threshold);
      SMsetSmoothThreshold(m_smoothData, threshold);
    }
    // Access info about quality
    // Spew smoothing stats
  };

  class OptMSSmoothingManager2D : public OptMSSmoothingManager {
  private:
    int
    smoothBdryVert(Vert* const vert, Face * const face,
		   BdryEdge* const bdryEdgeInit);
    int
    smoothIntBdryVert(Vert* const vert, Face * const face,
		      IntBdryEdge* const bdryEdgeInit);
  public:
    OptMSSmoothingManager2D(Mesh2D * const pMeshToSmooth, int qualMeasure =
    OPTMS_DEFAULT);
    virtual
    ~OptMSSmoothingManager2D();
    void
    setSmoothingGoal(const int goal);
    int
    smoothVert(Vert* vert);
  };

  class OptMSSmoothingManager3D : public OptMSSmoothingManager {
  private:
    std::set<Cell*> m_incCells;
    std::set<Vert*> m_neighVerts;
    std::set<BFace*> m_incBdryFaces;
    std::set<Face*> m_nearFaces; // All faces in tets incident on pV.
  public:
    OptMSSmoothingManager3D(VolMesh * const pMeshToSmooth, int qualMeasure =
    OPTMS_DEFAULT);
    virtual
    ~OptMSSmoothingManager3D();
    void
    setSmoothingGoal(const int goal);
    int
    smoothVert(Vert* vert);
  };
}
;

#endif
