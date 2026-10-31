#ifndef GR_Observable_h
#define GR_Observable_h

#include <vector>
#include <map>

#include "GR_config.h"
#include <omp.h>
#include <GR_EntContainer.h>

namespace GRUMMP
{
  class Observer;

  class Observable {
  public:
    // The possible events that a subscriber might want to know about.
    // They're numbered the way they are so that an observer can
    // subscribe to all that it wants at one time.
    enum ObsEvents {
      vertCreated = 0x1, edgeCreated = 0x2, // For future use in 3D
      faceCreated = 0x4,
      cellCreated = 0x8,
      bfaceCreated = 0x10,
      vertDeleted = 0x100,
      edgeDeleted = 0x200, // For future use in 3D
      faceDeleted = 0x400,
      cellDeleted = 0x800,
      bfaceDeleted = 0x1000,
      vertMoved = 0x10000
    };
    Observable();

    // Events that haven't been sent before the Observable is destroyed
    // don't get forwarded on.  But since the Observable -is- the mesh,
    // that doesn't matter, because any remaining Observers wouldn't be
    // able to do anything with the mesh after that anyway.
    virtual
    ~Observable();
    // Return value: all observed events.
    size_t
    addObserver(Observer* obs, const unsigned events);
    // Return value: all observed events.
    size_t
    removeObserver(Observer* obs);

    void
    deleteVertEvent(Vert* v);
    void
    deleteFaceEvent(Face* f);
    void
    deleteCellEvent(Cell* c);
    void
    deleteBFaceEvent(BFace* bf);

    void
    createVertEvent(Vert* v);
    void
    createFaceEvent(Face* f);
    void
    createCellEvent(Cell* c);
    void
    createBFaceEvent(BFace* bf);

    void
    moveVertEvent(Vert* v);

    void
    sendEvents();
    void
    clearEventQueues();
    void
    isOKToSend(const bool isOK)
    {
      m_OKToSend = isOK;
    }
    void
    IsInParallel(const bool ParFlag)
    {
      IsParallel = ParFlag;
    }
  private:
    template<class T>
      void
      deleteEvent(T* t, std::vector<T*>& deletedEnts,
		  std::vector<T*>& createdEnts);
    template<class T>
      void
      createEvent(T* t, std::vector<T*>& createdEnts);

    typedef std::pair<Observer*, unsigned int> ObsPair;
    typedef std::multimap<Observer*, unsigned int> ObsContainer;
    typedef std::multimap<Observer*, unsigned int>::iterator ObsContIter;

    std::vector<Cell*> m_deletedCells, m_createdCells;
    std::vector<BFace*> m_deletedBFaces, m_createdBFaces;
    std::vector<Face*> m_deletedFaces, m_createdFaces;
    std::vector<Vert*> m_deletedVerts, m_createdVerts;
    std::vector<Vert*> m_movedVerts;

    //*********** parallel Variables ************
    std::vector<Cell*> Para_m_deletedCells[NUM_PROCS],
	Para_m_createdCells[NUM_PROCS];
    std::vector<BFace*> Para_m_deletedBFaces[NUM_PROCS],
	Para_m_createdBFaces[NUM_PROCS];
    std::vector<Face*> Para_m_deletedFaces[NUM_PROCS],
	Para_m_createdFaces[NUM_PROCS];
    std::vector<Vert*> Para_m_deletedVerts[NUM_PROCS],
	Para_m_createdVerts[NUM_PROCS];
    std::vector<Vert*> Para_m_movedVerts[NUM_PROCS];
    unsigned int Para_m_eventCount[NUM_PROCS];

    //*******************************************
    ObsContainer m_observers;
    unsigned int m_allObservedEvents, m_eventCount, m_eventThreshold;

    bool IsParallel;
    bool m_OKToSend;
    bool IsSwapping2D;
    bool IsSwapping3D;
  };
}

#endif
