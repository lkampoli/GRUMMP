#ifndef GR_TETREFINERQUEUE_H
#define GR_TETREFINERQUEUE_H

#include <float.h>

#include <deque>
#include <queue>
#include <vector>

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_BFace.h"
#include "GR_Cell.h"
#include "GR_Subseg.h"

//TODO: REMOVE THIS

//class TetRefinerEntry {
//
//  TetRefinerEntry(const TetRefinerEntry&) :
//      m_priority(LARGE_DBL)
//  {
//    assert(0);
//  }
//  TetRefinerEntry& operator=(const TetRefinerEntry&)
//  {
//    assert(0);
//    return *this;
//  }
//
//protected:
//
//  double m_priority;
//
//  TetRefinerEntry(double priority = LARGE_DBL) :
//      m_priority(priority)
//  {
//  }
//
//public:
//
//  virtual ~TetRefinerEntry()
//  {
//  }
//
//  //Types of entries (derived classes).
//  enum EntryType {
//    UNKNOWN, SUBSEG, SUBFACE, CELL
//  };
//
//  virtual EntryType type() const = 0;
//
//  virtual bool deleted() const = 0;
//
//  bool not_deleted() const
//  {
//    return !deleted();
//  }
//
//  void set_priority(const double priority)
//  {
//    m_priority = priority;
//  }
//
//  double get_priority() const
//  {
//    return m_priority;
//  }
//
//  //These three functions are redefined in the base classes.
//  virtual Subseg* get_subseg() const
//  {
//    assert(0);
//    vFatalError("", "TetRefinerEntry::get_subseg()");
//    return NULL;
//  }
//  virtual BFace* get_subface() const
//  {
//    assert(0);
//    vFatalError("", "TetRefinerEntry::get_subface()");
//    return NULL;
//  }
//  virtual Cell* get_cell() const
//  {
//    assert(0);
//    vFatalError("", "TetRefinerEntry::get_cell()");
//    return NULL;
//  }
//
//  struct PriorityOrder {
//    bool operator()(const TetRefinerEntry* const entry1,
//                    const TetRefinerEntry* const entry2)
//    {
//      return (entry1->m_priority < entry2->m_priority);
//    }
//  };
//
//};
//
//class SubsegRefinerEntry: public TetRefinerEntry {
//
//  Subseg* m_subseg;
//
//  SubsegRefinerEntry() :
//      TetRefinerEntry(), m_subseg(NULL)
//  {
//  }
//
//  SubsegRefinerEntry& operator=(const SubsegRefinerEntry&)
//  {
//    assert(0);
//    return *this;
//  }
//
//public:
//
//  SubsegRefinerEntry(Subseg* const subseg, double priority = LARGE_DBL) :
//      TetRefinerEntry(priority), m_subseg(subseg)
//  {
//    assert(m_subseg->is_not_deleted());
//  }
//
//  SubsegRefinerEntry(const SubsegRefinerEntry& SRE) :
//      TetRefinerEntry(SRE.m_priority), m_subseg(SRE.m_subseg)
//  {
//    assert(m_subseg->is_not_deleted());
//  }
//
//  virtual ~SubsegRefinerEntry()
//  {
//  }
//
//  EntryType type() const
//  {
//    return TetRefinerEntry::SUBSEG;
//  }
//
//  bool deleted() const
//  {
//    return m_subseg->is_deleted();
//  }
//
//  Subseg* get_subseg() const
//  {
//    return m_subseg;
//  }
//
//};
//
//class SubfaceRefinerEntry: public TetRefinerEntry {
//
//  BFace* m_subface;
//  Vert * m_verts[3];
//  SubfaceRefinerEntry& operator=(const SubfaceRefinerEntry&)
//  {
//    assert(0);
//    return *this;
//  }
//
//public:
//
//  SubfaceRefinerEntry(BFace* const subface, double priority = LARGE_DBL) :
//      TetRefinerEntry(priority), m_subface(subface)
//  {
//    assert(!m_subface->isDeleted());
//    assert(
//        m_subface->getType() == Cell::eTriBFace || m_subface->getType() == Cell::eIntTriBFace);
//    for(int iV = 0; iV < 3; iV++)
//    	m_verts[iV] =  subface->getFace(0)->getVert(iV);
//  }
//
//  SubfaceRefinerEntry(const SubfaceRefinerEntry& SRE) :
//      TetRefinerEntry(SRE.m_priority), m_subface(SRE.m_subface)
//  {
//    assert(!m_subface->isDeleted());
//    assert(
//        m_subface->getType() == Cell::eTriBFace || m_subface->getType() == Cell::eIntTriBFace);
//    for(int iV = 0; iV < 3; iV++)
//    	m_verts[iV] =  m_subface->getFace(0)->getVert(iV);
//  }
//
//  virtual ~SubfaceRefinerEntry()
//  {
//  }
//
//  EntryType type() const
//  {
//    return TetRefinerEntry::SUBFACE;
//  }
//
//  bool deleted() const
//  {
//    if(m_subface->isDeleted()) return true;
//    for(int iV = 0; iV < 3; iV++)
//    	if(m_verts[iV] != m_subface->getFace(0)->getVert(iV)) return true;
//    return false;
//
//  }
//
//  BFace* get_subface() const
//  {
//    return m_subface;
//  }
//
//};
//
//class CellRefinerEntry: public TetRefinerEntry {
//
//  Cell* m_cell;
//
//  CellRefinerEntry& operator=(const CellRefinerEntry&)
//  {
//    assert(0);
//    return *this;
//  }
//
//public:
//
//  CellRefinerEntry(const CellRefinerEntry& CRE) :
//      TetRefinerEntry(CRE.m_priority), m_cell(CRE.m_cell)
//  {
//    assert(!m_cell->isDeleted());
//    assert(m_cell->getType() == Cell::eTet);
//  }
//
//  CellRefinerEntry(Cell* const cell, double priority) :
//      TetRefinerEntry(priority), m_cell(cell)
//  {
//    assert(!m_cell->isDeleted());
//    assert(m_cell->getType() == Cell::eTet);
//  }
//
//  virtual ~CellRefinerEntry()
//  {
//  }
//
//  EntryType type() const
//  {
//    return TetRefinerEntry::CELL;
//  }
//
//  bool deleted() const
//  {
//    return m_cell->isDeleted();
//  }
//
//  Cell* get_cell() const
//  {
//    return m_cell;
//  }
//
//};
//
//class PriorityCalculator {
//
//public:
//
//  enum PriorityType {
//    NONE, RADIUS_EDGE
//  };
//
//private:
//
//  PriorityType m_type;
//
//  PriorityCalculator() :
//      m_type(NONE)
//  {
//    assert(0);
//  }
//  PriorityCalculator(const PriorityCalculator&) :
//      m_type(NONE)
//  {
//    assert(0);
//  }
//  PriorityCalculator& operator=(const PriorityCalculator&)
//  {
//    assert(0);
//    return *this;
//  }
//
//  double radius_edge_priority(const Cell* const cell) const;
//
//  bool known_type()
//  {
//    return (m_type == RADIUS_EDGE);
//  }
//
//public:
//
//  PriorityCalculator(PriorityType type) :
//      m_type(type)
//  {
//    assert(known_type());
//  }
//  virtual ~PriorityCalculator()
//  {
//  }
//
//  double compute_priority(const Cell* const entry) const;
//
//};
//
//class TetRefinerPrioQueue: public std::priority_queue<TetRefinerEntry*,
//    std::vector<TetRefinerEntry*>, TetRefinerEntry::PriorityOrder> {
//
//  typedef std::vector<TetRefinerEntry*> EntryVec;
//
//public:
//
//  //Returns a reference to the vector upon which the queue is built.
//  EntryVec& impl()
//  {
//    return c;
//  }
//
//  void clean_queue();
//
//};
//
//class TetRefinerQueue {
//
//  //Do not need to use priority queues for subsegs and subfaces.
//  typedef std::deque<SubsegRefinerEntry*> SubsegQueue;
//  typedef std::deque<SubfaceRefinerEntry*> SubFaceQueue;
//
//  //The three queues that will control refinement;
//  SubsegQueue m_subseg_queue;
//  SubFaceQueue m_subface_queue;
//  TetRefinerPrioQueue m_cell_queue;
//
//  //Used to compute priority of cells.
//  PriorityCalculator* m_calculator;
//
//  TetRefinerQueue() :
//      m_subseg_queue(), m_subface_queue(), m_cell_queue(), m_calculator(NULL)
//  {
//    assert(0);
//  }
//  TetRefinerQueue(const TetRefinerQueue&) :
//      m_subseg_queue(), m_subface_queue(), m_cell_queue(), m_calculator(NULL)
//  {
//    assert(0);
//  }
//  TetRefinerQueue& operator=(const TetRefinerQueue&)
//  {
//    assert(0);
//    return *this;
//  }
//
//public:
//
//  TetRefinerQueue(PriorityCalculator::PriorityType prio_type) :
//      m_subseg_queue(), m_subface_queue(), m_cell_queue(), m_calculator(
//          new PriorityCalculator(prio_type))
//  {
//  }
//
//  virtual ~TetRefinerQueue()
//  {
//    if (m_calculator) delete m_calculator;
//    empty_queue();
//  }
//
//  void insert_subseg_front(Subseg* const subseg);
//  void insert_subseg_back(Subseg* const subseg);
//  void insert_subface_front(BFace* const subface);
//  void insert_subface_back(BFace* const subface);
//  void insert_cell(Cell* const cell);
//
//  TetRefinerEntry* get_top_entry() const;
//
//  bool pop_top_entry();
//
//  void empty_queue()
//  {
//    while (pop_top_entry()) {
//    }
//  }
//
//  void clean_queue();
//
//  bool empty() const
//  {
//    return (m_subseg_queue.empty() && m_subface_queue.empty()
//        && m_cell_queue.empty());
//  }
//
//};
//
//inline double PriorityCalculator::compute_priority(const Cell* const cell) const
//{
//
//  switch (m_type) {
//
//    case RADIUS_EDGE:
//      return radius_edge_priority(cell);
//    default:
//      assert(0);
//      break;
//  }
//  return DBL_MAX;
//}
//
//inline void TetRefinerQueue::insert_subseg_front(Subseg* const subseg)
//{
//  m_subseg_queue.push_front(new SubsegRefinerEntry(subseg));
//}
//
//inline void TetRefinerQueue::insert_subseg_back(Subseg* const subseg)
//{
//  m_subseg_queue.push_back(new SubsegRefinerEntry(subseg));
//}
//
//inline void TetRefinerQueue::insert_subface_front(BFace* const bface)
//{
//  m_subface_queue.push_front(new SubfaceRefinerEntry(bface));
//}
//
//inline void TetRefinerQueue::insert_subface_back(BFace* const bface)
//{
//  m_subface_queue.push_back(new SubfaceRefinerEntry(bface));
//}
//
//inline void TetRefinerQueue::insert_cell(Cell* const cell)
//{
//  double priority = m_calculator->compute_priority(cell);
//  //if( priority )
//  m_cell_queue.push(new CellRefinerEntry(cell, priority));
//}

#endif
