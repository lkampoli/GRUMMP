#include "GR_Cell.h"
#include "GR_SmoothingManager.h"
#include "GR_Vertex.h"
#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"
#include "GR_Face.h"
namespace GRUMMP
{
  SmoothingManager::~SmoothingManager()
  {
    m_pMesh->removeObserver(this);
    // Smoothing stats will eventually be reported here.
  }

  void
  SmoothingManager::receiveCreatedCells(std::vector<Cell*>& createdCells)
  {
    logMessage(3, "SmoothingManager received word about %zu new cells.\n",
	       createdCells.size());
    std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
	createdCells.end();
    for (; iter != iterEnd; ++iter) {
      Cell *pC = *iter;
      if (pC->isDeleted()
	  || (pC->getType() != CellSkel::eTet
	      && pC->getType() != CellSkel::eTriCell))
	continue;
      for (int i = pC->getNumVerts() - 1; i >= 0; i--) {
	Vert *pV = pC->getVert(i);
	logMessage(4, "  Queueing vert %p\n", pV);
	vertQueue.insert(pV);
      }
    }
  }

  void
  SmoothingManager::receiveDeletedVerts(std::vector<Vert*>& deletedVerts)
  {
    logMessage(3, "SmoothingManager received word about %zu deleted verts.\n",
	       deletedVerts.size());
    std::vector<Vert*>::iterator iter = deletedVerts.begin(), iterEnd =
	deletedVerts.end();
    for (; iter != iterEnd; ++iter) {
      Vert *pV = *iter;
      logMessage(4, "  De-queueing vert %p\n", pV);
      vertQueue.erase(pV);
    }
  }

  void
  SmoothingManager::receiveMovedVerts(std::vector<Vert*>& movedVerts)
  {
    logMessage(3, "SmoothingManager received word about %zu moved verts.\n",
	       movedVerts.size());
    std::vector<Vert*>::iterator iter = movedVerts.begin(), iterEnd =
	movedVerts.end();
    for (; iter != iterEnd; ++iter) {
      Vert *pV = *iter;
      for (int i = pV->getNumFaces() - 1; i >= 0; i--) {
	Face *pF = pV->getFace(i);
	for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
	  logMessage(4, "  Queueing vert %p\n", pF->getVert(ii));
	  vertQueue.insert(pF->getVert(ii));
	}
      }
    }
  }

  int
  SmoothingManager::smoothAllQueuedVerts()
  {

    // Initialize the active queue with all the verts in the master queue.
    std::set<Vert*> liveVertQueue;
    swap(vertQueue, liveVertQueue);
    vertQueue.clear();
    logMessage(3, "Live vert queue contains %zu entries.\n",
	       liveVertQueue.size());
    std::set<Vert*>::iterator iter = liveVertQueue.begin(), iterEnd =
	liveVertQueue.end();
    int numSmoothed = 0;
    for (; iter != iterEnd; ++iter) {
      Vert *pV = *iter;
      numSmoothed += smoothVert(pV);
    }
    return numSmoothed;
  }

  int
  SmoothingManager::smoothAllVerts(const int maxPasses)
  {
    // Initialize the active queue with all the verts in the mesh.
    vertQueue.clear();
    EntContainer<Vert>::iterator ECiter = m_pMesh->vert_begin(), ECiterEnd =
	m_pMesh->vert_end();
    for (; ECiter != ECiterEnd; ++ECiter) {
      Vert *pV = &(*ECiter);
      if (pV->getNumFaces() != 0)
	vertQueue.insert(pV);

    }

    int numSmoothed = 0, smoothedThisPass = 0, passes = 0;
    do {
      int attempts = vertQueue.size();
      smoothedThisPass = smoothAllQueuedVerts();
      numSmoothed += smoothedThisPass;
      passes++;

      logMessage(2, "Pass %d, %d smoothed out of %d attempts\n", passes,
		 smoothedThisPass, attempts);
      m_pMesh->sendEvents();
      logMessage(2, "  Master smoothing queue now has %zu verts.\n",
		 vertQueue.size());
    }
    while (smoothedThisPass != 0 && passes < maxPasses);
    return numSmoothed;
  }

  OptMSSmoothingManager::OptMSSmoothingManager(Mesh * const pMeshToSmooth,
					       const int iDim) :
      SmoothingManager(pMeshToSmooth), m_smoothData(NULL), m_neighCoords(
      NULL), m_faceConn(NULL)
  {
    SMinitSmoothing(iDim, 'f', OPTMS_DEFAULT, OPTMS_DEFAULT,
		    reinterpret_cast<void**>(&m_smoothData));
  }

  OptMSSmoothingManager*
  OptMSSmoothingManager::Create(Mesh * const pMesh)
  {
    switch (pMesh->getType())
      {
      case Mesh::eMesh2D:
	return new OptMSSmoothingManager2D(dynamic_cast<Mesh2D*>(pMesh));
	break;
      case Mesh::eVolMesh:
	return new OptMSSmoothingManager3D(dynamic_cast<VolMesh*>(pMesh));
	break;
      default:
	return NULL;
      }
  }
}
