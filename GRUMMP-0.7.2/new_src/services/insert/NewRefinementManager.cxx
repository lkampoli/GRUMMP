/*
 * NewRefinementManager.cxx
 *
 *  Created on: Dec 9, 2015
 *      Author: cfog
 */
#include <map>
#include <queue>
#include <memory>

#include "GR_InsertionManager.h"
#include "GR_PriorityCalculator.h"

namespace GRUMMP
{
  class QueueEntry {
  };
  class Queue {
  public:
    Queue();
    virtual
    ~Queue()
    {
    }
    virtual QueueEntry*
    top() = 0;
    virtual void
    pop() = 0;
    virtual void
    addToQueue(QueueEntry*) = 0;
    virtual void
    initialize(Mesh*) = 0;
    virtual bool
    getInsertionPoint(QueueEntry*, double candLoc[], Mesh*) = 0;
    virtual void
    insertPoint(QueueEntry*, double loc[], Mesh*) = 0;
    virtual bool
    empty() = 0;
  };
  class CompositeQueue : public Queue {
  private:
    std::map<int, std::shared_ptr<Queue> > m_queueList;
  public:
    QueueEntry*
    top()
    {
      std::map<int, std::shared_ptr<Queue> >::iterator iter, iterEnd =
	  m_queueList.end();
      for (iter = m_queueList.begin(); iter != iterEnd; ++iter) {
	std::shared_ptr<Queue> temp = iter->second;
	if (!temp->empty()) {
	  return temp->top();
	}
      }
      return nullptr;
    }
  };
  template<class T>
    class RefinementPriorityQueue : Queue {
    private:
      std::priority_queue<T> m_prioQueue;
      std::shared_ptr<PriorityCalculator> m_PC;
      double m_threshold;
    public:
      RefinementPriorityQueue(std::shared_ptr<PriorityCalculator>& ap_PC,
			      std::shared_ptr<InsertionPointCalculator>& ap_IPC,
			      std::shared_ptr<Inserter>& ap_I,
			      const double threshold);

      double
      getThreshold() const
      {
	return m_threshold;
      }

      void
      setThreshold(double threshold)
      {
	m_threshold = threshold;
      }

      QueueEntry*
      top()
      {
	return m_prioQueue.top();
      }
      void
      pop()
      {
	m_prioQueue.pop_front();
      }
      void
      addToQueue(T t)
      {/* calc priority and insert */
      }
    };

  template<class T>
    class RefinementQueue {
    private:
      std::queue<T> m_queue;
    public:
      RefinementQueue(std::shared_ptr<InsertionPointCalculator>& ap_IPC,
		      std::shared_ptr<Inserter>& ap_I);
      T&
      top()
      {
	return m_queue.top();
      }
      void
      pop()
      {
	m_queue.pop_front();
      }
      void
      addToQueue(T t)
      {
	m_queue.push(t);
      }
    };

// The container for these needs to have the same basic set of calls
// with internal dispatch to actual queues.  T needs to be some object
// with an operator<() and enough info to keep PC and IPC happy.
  class RefManager {
    std::shared_ptr<Queue> m_queue;
    Mesh *m_pMesh;
  public:
    RefManager(std::shared_ptr<Queue> queue, Mesh* mesh) :
	m_queue(queue), m_pMesh(mesh)
    {
    }
    void
    refineMesh();
  };

  void
  RefManager::refineMesh()
  {
    // Queue all entities: each subqueue iterates over the parts of the
    // mesh it cares about.  Not sure about subsegs in 3D....
    m_queue->initialize(m_pMesh);
    do {
      // Grab an entry from the top of the queue; don't need to know what kind
      QueueEntry* pQE = m_queue->top();
      double candLoc[3];
      // Calculate the insertion point and decide whether to accept or reject
      bool OKtoInsert = m_queue->getInsertionPoint(pQE, candLoc, m_pMesh);
      // Do the actual insertion
      if (OKtoInsert)
	m_queue->insertPoint(pQE, candLoc, m_pMesh);
    }
    while (!m_queue->empty());
  }
}

