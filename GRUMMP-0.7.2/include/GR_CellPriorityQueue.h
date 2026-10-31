#ifndef GR_CellPriorityQueue
#define GR_CellPriorityQueue 1

#include <set>
#include <map>
#include <memory>
#include <bits/shared_ptr.h> // To make eclipse happy
#include <limits.h>

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_Mesh.h"
#include "GR_PriorityCalculator.h"

namespace GRUMMP
{
  class CellPriorityQueueEntry {
  private:
    double m_newPtLoc[3];
    Cell *m_cell;
    Vert *m_verts[4];
    // Because bit-wise copy is correct for this class, the default
    // operator= and copy constructor don't need to be overridden.
  public:
    CellPriorityQueueEntry(Cell * const cell) :
	m_cell(cell)
    {
      // An obviously bogus location value...
      m_newPtLoc[0] = m_newPtLoc[1] = m_newPtLoc[2] = -1e300;
      m_verts[3] = pVInvalidVert; // May get over-written.
      for (int i = 0; i < cell->getNumVerts(); i++) {
	m_verts[i] = cell->getVert(i);
      }
    }
    void
    setPointLocation(double newLoc[])
    {
      m_newPtLoc[0] = newLoc[0];
      m_newPtLoc[1] = newLoc[1];
      // The last component may be garbage, which is okay, but it may
      // also be illegal to read that location, which is not.
      if (m_cell->getType() == Cell::eTet)
	m_newPtLoc[2] = newLoc[2];
    }

    const double*
    getPointLocation() const
    {
      return m_newPtLoc;
    }
    Cell *
    getCell() const
    {
      return m_cell;
    }
    Vert *
    getVert(const int i) const
    {
      assert(i >= 0 && i < m_cell->getNumVerts());
      return (m_verts[i]);
    }
    bool
    isStillInMesh() const;
    bool
    operator==(const CellPriorityQueueEntry& IQE) const;
    bool
    operator!=(const CellPriorityQueueEntry& IQE)
    {
      return (!(*this == IQE));
    }
  };

  class CellPriorityQueue {
  private:
    std::multimap<double, CellPriorityQueueEntry> m_MMap;
    double m_worstAllowable;
    std::shared_ptr<PriorityCalculator> m_prioCalc;
    // Don't see any good reason to allow default or copy construction, or
    // operator=.
    CellPriorityQueue(); // Never defined
    CellPriorityQueue(const CellPriorityQueue&); // Never defined
    CellPriorityQueue
    operator=(const CellPriorityQueue&); // Never defined
  public:
    // Note that multimaps with double keys sort from smallest to largest
    // by default, so big negative values go to the the head of the queue.
    enum dPriorityValues {
      dWorstCellPriority = INT_MIN / 16,
      dNeverQueuedPriority = INT_MAX,
      dInvalidPriority = INT_MAX
    };
    CellPriorityQueue(std::shared_ptr<PriorityCalculator> prioCalc,
		      const double target = CellPriorityQueue::dInvalidPriority) :
	m_MMap(), m_worstAllowable(target), m_prioCalc(prioCalc)
    {
    }
    ~CellPriorityQueue()
    {
    }

    bool
    getNextCell(CellPriorityQueueEntry &entry)
    {
      // Grab the first entry in the multimap, and return the data.
      // m_MMap.begin() returns an iterator which (through overloading of
      // operator ->) can access the contents of the pair<double,
      // CellPriorityQueueEntry> that constitutes the first mmap entry.
      cleanQueue();
      if (isQueueEmpty()) {
	return false;
      }
      else {
	entry = m_MMap.begin()->second;
	logMessage(3, "Grabbed queue entry with quality measure %f\n",
		   m_MMap.begin()->first);
	return true;
      }
    }
  private:
    void
    cleanQueue()
    // Remove all null entries from the top of the queue.
    {
      if (isQueueEmpty())
	return;
      logMessage(4, "Cleaning queue of length %u...", getQueueLength());
      do {
	CellPriorityQueueEntry entry = m_MMap.begin()->second;
	if (entry.isStillInMesh()) {
	  logMessage(4, "now of length %u.\n", getQueueLength());
	}
	m_MMap.erase(m_MMap.begin());
      }
      while (!isQueueEmpty());
    }

    double
    getTopPriority()
    {
      assert(!isQueueEmpty());
      // Grab the first entry in the multimap, and return its priority.
      // m_MMap.begin() returns an iterator which (through overloading of
      // operator ->) can access the contents of the pair<double,

      // CellPriorityQueueEntry> that constitutes the first mmap entry.
      return m_MMap.begin()->first;
    }
    void
    printQueue()
    {
      std::multimap<double, CellPriorityQueueEntry>::iterator it;
      for (it = m_MMap.begin(); it != m_MMap.end(); it++) {
	printf("%p: Still in mesh = %d, deleted = %d, priority = %lf\n",
	       it->second.getCell(), it->second.isStillInMesh(),
	       it->second.getCell()->isDeleted(), it->first);
      }
    }

    GR_index_t
    getQueueLength() const
    {
      return m_MMap.size();
    }
    bool
    isQueueEmpty()
    {
      return m_MMap.empty();
    }
  public:
    bool
    addEntry(SimplexCell* const cell)
    {
      // Return value of insert is an iterator; not currently checked.
      double key = m_prioCalc->calcPriority(cell);

      if (key < m_worstAllowable) {
	/* This hack is designed to fix interior tets.  It works.  But
	 most of the bad ones tend to be at the boundaries, so it
	 doesn't have much net effect. */
	/* || */
	/* 	  (IQE.eType() == CellPriorityQueueEntry::eTetCell && */
	/* 	   (!IQE.pCCell()->pVVert(0)->qIsBdryVert() && */
	/* 	    !entry.getCell()->pVVert(1)->qIsBdryVert() && */
	/* 	    !entry.getCell()->pVVert(2)->qIsBdryVert() && */
	/* 	    !entry.getCell()->pVVert(3)->qIsBdryVert()) && */
	/* 	   dKey < m_worstAllowable*2) */
	/* 	if (entry.eType() == CellPriorityQueueEntry::eTetCell || */
	/* 	    entry.eType() == CellPriorityQueueEntry::eTriCell) { */
	/* 	  double adCent[3]; adCent[2] = 0; */
	/* 	  IQE.getCell()->vCentroid(adCent); */
	/* 	  vMessage(3, "Key: %10.6f.  Centroid: (%10.6f, %10.6f, %10.6f)\n", */
	/* 		   dKey, adCent[0], adCent[1], adCent[2]); */
	/* 	} */
	CellPriorityQueueEntry entry(cell);
	m_MMap.insert(std::make_pair(key, entry));
	return true;
      }
      else {
	return false;
      }
    }
    double
    setQualityTarget(const double newTarget)
    {
      double temp = m_worstAllowable;
      m_worstAllowable = newTarget;
      return temp;
    }
    double
    getQualityTarget() const
    {
      return m_worstAllowable;
    }
  };
} // End namespace GRUMMP
#endif

