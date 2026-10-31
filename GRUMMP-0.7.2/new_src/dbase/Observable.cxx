#include <algorithm>

#include "GR_Observable.h"
#include "GR_Observer.h"
#include "GR_misc.h"

namespace GRUMMP
{

  Observable::Observable() :
      m_deletedCells(), m_createdCells(), m_deletedBFaces(), m_createdBFaces(), m_deletedFaces(), m_createdFaces(), m_deletedVerts(), m_createdVerts(), m_movedVerts(), m_observers(), m_allObservedEvents(
	  0), m_eventCount(0), m_eventThreshold(100000), IsParallel(false), m_OKToSend(
	  true), IsSwapping2D(false), IsSwapping3D(false)
  {
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      Para_m_eventCount[ID] = 0;
//		  Para_m_deletedCells[ID] = 0;	Para_m_createdCells[ID] = 0;
//		  Para_m_deletedBFaces[ID] = 0;	Para_m_createdBFaces[ID] = 0;
//		  Para_m_deletedFaces[ID] = 0;	Para_m_createdFaces[ID] = 0;
//		  Para_m_deletedVerts[ID] = 0;	Para_m_createdVerts[ID] = 0;
//		  Para_m_movedVerts[ID] = 0;
    }
  }

  Observable::~Observable()
  {
    ObsContIter iter = m_observers.begin(), iterEnd = m_observers.end();
    for (; iter != iterEnd; ++iter) {
      Observer* obs = iter->first;
      obs->detachFromObservable();
    }
  }

  size_t
  Observable::addObserver(Observer* obs, const unsigned events)
  {
    // An Observer have more than one entry in the map.
    m_observers.insert(ObsPair(obs, events));
    m_allObservedEvents |= events;
    logMessage(3, "Added observer of events %X; all observed events: %X\n",
	       events, m_allObservedEvents);
    return m_allObservedEvents;
  }

  size_t
  Observable::removeObserver(Observer* obs)
  {
    m_observers.erase(obs);
    m_allObservedEvents = 0;
    ObsContIter iter = m_observers.begin(), iterEnd = m_observers.end();
    for (; iter != iterEnd; ++iter) {
      unsigned int oe = iter->second;
      m_allObservedEvents |= oe;
    }
    return m_allObservedEvents;
  }

  template<class T>
    void
    Observable::deleteEvent(T* t, std::vector<T*>& deletedEnts,
			    std::vector<T*>& createdEnts)
    {
      if (!IsParallel) {
	typename std::vector<T*>::iterator iter = std::find(createdEnts.begin(),
							    createdEnts.end(),
							    t);

	if (iter != createdEnts.end()) {
	  createdEnts.erase(iter);
	  m_eventCount--;
	}
	else {
	  deletedEnts.push_back(t);
	  m_eventCount++;
	}
	//if (m_eventCount >= m_eventThreshold) sendEvents();
      }
      else {
	int ID = omp_get_thread_num();
//		    typename std::vector<T*>::iterator iter =
//		      std::find(createdEnts.begin(), createdEnts.end(), t);
//
//		    if (iter != createdEnts.end()) {
//		      createdEnts.erase(iter);
//		      Para_m_eventCount[ID]--;
//		    }
//		    else {
	deletedEnts.push_back(t);
	Para_m_eventCount[ID]++;
//		    }
	//if (m_eventCount >= m_eventThreshold) sendEvents();
      }
    }

  void
  Observable::deleteVertEvent(Vert* v)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (vertDeleted | vertCreated | vertMoved)) {
	deleteEvent(v, m_deletedVerts, m_createdVerts);
	if (m_allObservedEvents & vertMoved) {
	  std::vector<Vert*>::iterator iter = std::find(m_movedVerts.begin(),
							m_movedVerts.end(), v);

	  if (iter != m_movedVerts.end()) {
	    m_movedVerts.erase(iter);
	    m_eventCount--;
	  }
	}
      }
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (vertDeleted | vertCreated | vertMoved)) {
	deleteEvent(v, Para_m_deletedVerts[ID], Para_m_createdVerts[ID]);
	if (m_allObservedEvents & vertMoved) {
	  std::vector<Vert*>::iterator iter = std::find(
	      Para_m_movedVerts[ID].begin(), Para_m_movedVerts[ID].end(), v);

	  if (iter != Para_m_movedVerts[ID].end()) {
	    Para_m_movedVerts[ID].erase(iter);
	    Para_m_eventCount[ID]--;
	  }
	}
      }
    }
  }

  void
  Observable::deleteFaceEvent(Face* f)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (faceDeleted | faceCreated))
	deleteEvent(f, m_deletedFaces, m_createdFaces);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (faceDeleted | faceCreated))
	deleteEvent(f, Para_m_deletedFaces[ID], Para_m_createdFaces[ID]);
    }
  }

  void
  Observable::deleteBFaceEvent(BFace* bf)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (bfaceDeleted | bfaceCreated))
	deleteEvent(bf, m_deletedBFaces, m_createdBFaces);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (bfaceDeleted | bfaceCreated))
	deleteEvent(bf, Para_m_deletedBFaces[ID], Para_m_createdBFaces[ID]);
    }
  }

  void
  Observable::deleteCellEvent(Cell* c)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (cellDeleted | cellCreated))
	deleteEvent(c, m_deletedCells, m_createdCells);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (cellDeleted | cellCreated))
	deleteEvent(c, Para_m_deletedCells[ID], Para_m_createdCells[ID]);
    }
  }

  template<class T>
    void
    Observable::createEvent(T* t, std::vector<T*>& createdEnts)
    {
      if (!IsParallel) {
	createdEnts.push_back(t);
	m_eventCount++;
      }
      else {
	int ID = omp_get_thread_num();
	createdEnts.push_back(t);
	Para_m_eventCount[ID]++;
      }

      //if (m_eventCount >= m_eventThreshold) sendEvents();
    }

  void
  Observable::createVertEvent(Vert* v)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (vertCreated | vertDeleted))
	createEvent(v, m_createdVerts);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (vertCreated | vertDeleted))
	createEvent(v, Para_m_createdVerts[ID]);
    }
  }

  void
  Observable::createFaceEvent(Face* f)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (faceCreated | faceDeleted))
	createEvent(f, m_createdFaces);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (faceCreated | faceDeleted))
	createEvent(f, Para_m_createdFaces[ID]);
    }
  }

  void
  Observable::createBFaceEvent(BFace* bf)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (bfaceCreated | bfaceDeleted))
	createEvent(bf, m_createdBFaces);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (bfaceCreated | bfaceDeleted))
	createEvent(bf, Para_m_createdBFaces[ID]);
    }
  }

  void
  Observable::createCellEvent(Cell* c)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & (cellCreated | cellDeleted))
	createEvent(c, m_createdCells);
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & (cellCreated | cellDeleted))
	createEvent(c, Para_m_createdCells[ID]);
    }
  }

  void
  Observable::moveVertEvent(Vert* v)
  {
    if (!IsParallel) {
      if (m_allObservedEvents & vertMoved) {
	m_movedVerts.push_back(v);
	m_eventCount++;
	if (m_eventCount >= m_eventThreshold)
	  sendEvents();
      }
    }
    else {
      int ID = omp_get_thread_num();
      if (m_allObservedEvents & vertMoved) {
	Para_m_movedVerts[ID].push_back(v);
	Para_m_eventCount[ID]++;
//		      if (m_eventCount >= m_eventThreshold) sendEvents();
      }
    }
  }

  void
  Observable::clearEventQueues()
  {
    m_deletedVerts.clear();
    m_deletedFaces.clear();
    m_deletedCells.clear();
    m_deletedBFaces.clear();
    m_createdVerts.clear();
    m_createdFaces.clear();
    m_createdCells.clear();
    m_createdBFaces.clear();
    m_movedVerts.clear();

    if (IsParallel) {
      for (int ID = 0; ID < NUM_PROCS; ID++) {
	Para_m_deletedVerts[ID].clear();
	Para_m_deletedFaces[ID].clear();
	Para_m_deletedCells[ID].clear();
	Para_m_deletedBFaces[ID].clear();

	Para_m_createdVerts[ID].clear();
	Para_m_createdFaces[ID].clear();
	Para_m_createdCells[ID].clear();
	Para_m_createdBFaces[ID].clear();
      }
    }

    m_eventCount = 0;
  }

  void
  Observable::sendEvents()
  {
    if (!m_OKToSend)
      return;
    ObsContIter iter = m_observers.begin(), iterEnd = m_observers.end();
    for (; iter != iterEnd; ++iter) {
      Observer* obs = iter->first;
      unsigned int oe = iter->second;

      logMessage(4,
		 "Sending info about mesh changes.  Recorded the following:\n");
      logMessage(4, "  Verts: %8zu deleted, %8zu created, %8zu moved.\n",
		 m_deletedVerts.size(), m_createdVerts.size(),
		 m_movedVerts.size());
      logMessage(4, "  Faces: %8zu deleted, %8zu created.\n",
		 m_deletedFaces.size(), m_createdFaces.size());
      logMessage(4, "  Cells: %8zu deleted, %8zu created.\n",
		 m_deletedCells.size(), m_createdCells.size());
      logMessage(4, "  BFaces: %8zu deleted, %8zu created.\n",
		 m_deletedBFaces.size(), m_createdBFaces.size());
      logMessage(4, "  Observer asking for %u\n", oe);

      if (!IsParallel) {
	// Send deleted stuff first.
	if (oe & vertDeleted)
	  obs->receiveDeletedVerts(m_deletedVerts);
	if (oe & faceDeleted)
	  obs->receiveDeletedFaces(m_deletedFaces);
	if (oe & cellDeleted)
	  obs->receiveDeletedCells(m_deletedCells);
	if (oe & bfaceDeleted)
	  obs->receiveDeletedBFaces(m_deletedBFaces);

	if (oe & vertCreated)
	  obs->receiveCreatedVerts(m_createdVerts);
	if (oe & faceCreated)
	  obs->receiveCreatedFaces(m_createdFaces);
	if (oe & cellCreated)
	  obs->receiveCreatedCells(m_createdCells);
	if (oe & bfaceCreated)
	  obs->receiveCreatedBFaces(m_createdBFaces);

	if (oe & vertMoved)
	  obs->receiveMovedVerts(m_movedVerts);

      }
      else {
	for (int ID = 0; ID < NUM_PROCS; ID++) {
	  if (oe & vertDeleted) {
	    std::vector<Vert*>::iterator it1 = Para_m_deletedVerts[ID].begin();
	    m_deletedVerts.insert(m_deletedVerts.end(), it1,
				  Para_m_deletedVerts[ID].end());
	  }
	  if (oe & faceDeleted) {
	    std::vector<Face*>::iterator it2 = Para_m_deletedFaces[ID].begin();
	    m_deletedFaces.insert(m_deletedFaces.end(), it2,
				  Para_m_deletedFaces[ID].end());
	  }
	  if (oe & cellDeleted) {
	    std::vector<Cell*>::iterator it3 = Para_m_deletedCells[ID].begin();
	    m_deletedCells.insert(m_deletedCells.end(), it3,
				  Para_m_deletedCells[ID].end());
	  }
	  if (oe & bfaceDeleted) {
	    std::vector<BFace*>::iterator it4 =
		Para_m_deletedBFaces[ID].begin();
	    m_deletedBFaces.insert(m_deletedBFaces.end(), it4,
				   Para_m_deletedBFaces[ID].end());
	  }
	  if (oe & vertCreated) {
	    std::vector<Vert*>::iterator it5 = Para_m_createdVerts[ID].begin();
	    m_createdVerts.insert(m_createdVerts.end(), it5,
				  Para_m_createdVerts[ID].end());
	  }
	  if (oe & faceCreated) {
	    std::vector<Face*>::iterator it6 = Para_m_createdFaces[ID].begin();
	    m_createdFaces.insert(m_createdFaces.end(), it6,
				  Para_m_createdFaces[ID].end());
	  }
	  if (oe & cellCreated) {
	    std::vector<Cell*>::iterator it = Para_m_createdCells[ID].begin();
	    m_createdCells.insert(m_createdCells.end(), it,
				  Para_m_createdCells[ID].end());
	  }
	  if (oe & bfaceCreated) {
	    std::vector<BFace*>::iterator it7 =
		Para_m_createdBFaces[ID].begin();
	    m_createdBFaces.insert(m_createdBFaces.end(), it7,
				   Para_m_createdBFaces[ID].end());
	  }
	}
	// Send deleted stuff first.
	if (oe & vertDeleted)
	  obs->receiveDeletedVerts(m_deletedVerts);
	if (oe & faceDeleted)
	  obs->receiveDeletedFaces(m_deletedFaces);
	if (oe & cellDeleted)
	  obs->receiveDeletedCells(m_deletedCells);
	if (oe & bfaceDeleted)
	  obs->receiveDeletedBFaces(m_deletedBFaces);

	if (oe & vertCreated)
	  obs->receiveCreatedVerts(m_createdVerts);
	if (oe & faceCreated)
	  obs->receiveCreatedFaces(m_createdFaces);
	if (oe & cellCreated)
	  obs->receiveCreatedCells(m_createdCells);
	if (oe & bfaceCreated)
	  obs->receiveCreatedBFaces(m_createdBFaces);

	if (oe & vertMoved)
	  obs->receiveMovedVerts(m_movedVerts);
      }

    } // Loop over observers
    clearEventQueues();
  } // sendEvents
} // namespace GRUMMP
