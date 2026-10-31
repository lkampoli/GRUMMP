#ifndef GR_SAMPLINGQUEUE_H
#define GR_SAMPLINGQUEUE_H

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include <queue>
#include <vector>
#include <map>

class Cell;
class Face;
class RefFace;
class SampleVertTag;
class SamplingEntryFactory;
class SubsegBridge;
class SurfaceSampler;

class SamplingEntry {

  int m_insertionAttempt;
  RefFace* m_ref_face;
  double m_priority;

  // Never used
  SamplingEntry&
  operator=(const SamplingEntry&);

protected:
  SamplingEntry() :
      m_insertionAttempt(0), m_ref_face(NULL), m_priority(-LARGE_DBL)
  {
  }
  SamplingEntry(const SamplingEntry& SE) :
      m_insertionAttempt(0), m_ref_face(SE.m_ref_face), m_priority(
	  SE.m_priority)
  {
  }

  SamplingEntry(RefFace* const ref_face, const double priority = -LARGE_DBL) :
      m_insertionAttempt(0), m_ref_face(ref_face), m_priority(priority)
  {
  }

public:

  struct QueueOrder {
    bool
    operator()(const SamplingEntry* const entry1,
	       const SamplingEntry* const entry2) const
    {
      return (entry1->m_priority < entry2->m_priority);
    }
  };

  enum EntryType {
    SUBSEG, MULTIPLE, BAD_TOPO, SIZE, QUALITY, UNKNOWN
  };

  virtual
  ~SamplingEntry()
  {
  }

  void
  set_priority(const double priority)
  {
    m_priority = priority;
  }

  double
  get_priority() const
  {
    return m_priority;
  }

  virtual EntryType
  get_entry_type() const = 0;

  virtual RefFace*
  get_surface() const
  {
    return m_ref_face;
  }

  virtual bool
  entry_deleted() const = 0;

  virtual bool
  remove_entry() const
  {
    return entry_deleted();
  }

  virtual void
  print_entry() const = 0;

  int
  getAttemptNum() const
  {
    return m_insertionAttempt;
  }
  void
  incrementAttemptNum()
  {
    m_insertionAttempt++;
  }

};

class SubsegSamplingEntry : public SamplingEntry {

  SubsegBridge* m_bridge;

  SubsegSamplingEntry(RefFace* const ref_face, SubsegBridge* const bridge,
		      const double priority = -LARGE_DBL) :
      SamplingEntry(ref_face, priority), m_bridge(bridge)
  {
  }

  SubsegSamplingEntry();
  SubsegSamplingEntry(const SubsegSamplingEntry&);
  SubsegSamplingEntry&
  operator=(const SubsegSamplingEntry&);

public:

  friend class SamplingEntryFactory;

  virtual
  ~SubsegSamplingEntry()
  {
  }

  EntryType
  get_entry_type() const
  {
    return SUBSEG;
  }

  SubsegBridge*
  get_bridge() const
  {
    return m_bridge;
  }

  bool
  entry_deleted() const;

  void
  print_entry() const;

};

class FaceSamplingEntry : public SamplingEntry {

  EntryType m_entry_type;

  FaceSamplingEntry(const FaceSamplingEntry&);

  FaceSamplingEntry&
  operator=(const FaceSamplingEntry&);

protected:

  FaceQueueEntry *m_fqe;
  Cell *m_cell_left, *m_cell_right;
  FaceSamplingEntry(EntryType entry_type, RefFace* const ref_face,
		    FaceQueueEntry* const fqe, const double priority =
			-LARGE_DBL);

  FaceSamplingEntry() :
      SamplingEntry(), m_entry_type(UNKNOWN), m_fqe(NULL), m_cell_left(
      NULL), m_cell_right(NULL)
  {
    assert(0);
  }

public:

  friend class SamplingEntryFactory;

  virtual
  ~FaceSamplingEntry()
  {
  }

  EntryType
  get_entry_type() const
  {
    return m_entry_type;
  }
  FaceQueueEntry*
  getFaceQueueEntry() const
  {
    return m_fqe;
  }
  Face*
  getFace() const
  {
    return m_fqe->getFace();
  }
  Cell*
  get_cell_left() const
  {
    return m_cell_left;
  }
  Cell*
  get_cell_right() const
  {
    return m_cell_right;
  }

  void
  set_face(Face* const face)
  {
    m_fqe->setFace(face);
  }
  void
  set_cell_left(Cell* const cell)
  {
    m_cell_left = cell;
  }
  void
  set_cell_right(Cell* const cell)
  {
    m_cell_right = cell;
  }
  virtual bool
  entry_deleted() const;

  void
  print_entry() const;

};

class FaceTopoSamplingEntry : public FaceSamplingEntry {

  Vert* m_offending_vert;
  SampleVertTag* m_offending_vert_tag;

  bool m_remove;

  FaceTopoSamplingEntry(RefFace* const ref_face, FaceQueueEntry* const fqe,
			Vert* const offending_vert,
			SampleVertTag* const offending_vert_tag,
			const double priority = -LARGE_DBL) :
      FaceSamplingEntry(SamplingEntry::BAD_TOPO, ref_face, fqe, priority), m_offending_vert(
	  offending_vert), m_offending_vert_tag(offending_vert_tag), m_remove(
	  false)
  {
  }

//  FaceTopoSamplingEntry()
//    : FaceSamplingEntry(), m_offending_vert(NULL),
//      m_offending_vert_tag(NULL), m_remove(true) { assert(0); }
  FaceTopoSamplingEntry(const FaceTopoSamplingEntry&) :
      FaceSamplingEntry(), m_offending_vert(NULL), m_offending_vert_tag(
      NULL), m_remove(true)
  {
    assert(0);
  }
  FaceTopoSamplingEntry&
  operator=(const FaceTopoSamplingEntry&)
  {
    assert(0);
    return *this;
  }

public:

  friend class SamplingEntryFactory;

  virtual
  ~FaceTopoSamplingEntry()
  {
  }

  void
  set_vertex(Vert* const vertex);

  Vert*
  get_vertex() const
  {
    return m_offending_vert;
  }
  SampleVertTag*
  get_vertex_tag() const
  {
    return m_offending_vert_tag;
  }

  void
  mark_for_removal()
  {
    m_remove = true;
  }

  virtual bool
  remove_entry() const
  {
    return m_remove;
  }

  virtual bool
  entry_deleted() const;

};

class SamplingPriorityQueue : public std::priority_queue<SamplingEntry*,
    std::vector<SamplingEntry*>, SamplingEntry::QueueOrder> {

  typedef std::vector<SamplingEntry*> EntryVec;

public:

  //The variable c is the vector upon which the queue is built,
  //and it is protected in the std::priority_queue.
  //impl() acts as an accessor.
  EntryVec
  impl()
  {
    return c;
  }

  //Removes deleted entries from the queue.
  void
  clean_queue();

  void
  rebuild_queue(const std::map<Vert*, Vert*>& vert_map,
		const std::map<Face*, Face*>& face_map,
		const std::map<Cell*, Cell*>& cell_map);

};

class SamplingQueue {

  SamplingPriorityQueue m_multiple;
  SamplingPriorityQueue m_subseg;
  SamplingPriorityQueue m_bad_topo;
  SamplingPriorityQueue m_size;
  SamplingPriorityQueue m_quality;

  //These are made private and should not be used.
  SamplingQueue(const SamplingQueue&);
  SamplingQueue&
  operator=(const SamplingQueue&);

public:

  //Creates an empty queue.
  SamplingQueue() :
      m_multiple(), m_subseg(), m_bad_topo(), m_size(), m_quality()
  {
  }

  //Destructor.
  ~SamplingQueue()
  {
    empty_queue();
  }

  //Returns the top entry in the queue.
  //CAUTION: The top entry might not be in the mesh anymore, so check before using.
  SamplingEntry*
  top_entry() const;

  //Deletes and removes the top entry from the queue.
  bool
  pop_top_entry();

  //Returns true if the queue is empty. False otherwise.
  bool
  is_empty() const
  {
    return (m_multiple.empty() && m_subseg.empty() && m_bad_topo.empty()
	&& m_size.empty() && m_quality.empty());
  }

  //Add a SamplingQueueEntry to the queue...
  void
  add_to_queue(SamplingEntry* const sampling_entry);

  int
  queue_size() const
  {
//    logMessage(1,"Sampling Entry Queue Size %lu %lu %lu %lu %lu \n",
//	       m_multiple.size(),m_subseg.size(), m_bad_topo.size(),
//	       m_size.size(),m_quality.size());

    return m_multiple.size() + m_subseg.size() + m_bad_topo.size()
	+ m_size.size() + m_quality.size();
  }

  void
  clean_queue()
  {

    m_multiple.clean_queue();
    m_subseg.clean_queue();
    m_bad_topo.clean_queue();
    m_size.clean_queue();
    m_quality.clean_queue();

  }

  void
  rebuild_queue(const std::map<Vert*, Vert*>& vert_map,
		const std::map<Face*, Face*>& face_map,
		const std::map<Cell*, Cell*>& cell_map)
  {

    m_multiple.rebuild_queue(vert_map, face_map, cell_map);
    m_subseg.rebuild_queue(vert_map, face_map, cell_map);
    m_bad_topo.rebuild_queue(vert_map, face_map, cell_map);
    m_size.rebuild_queue(vert_map, face_map, cell_map);
    m_quality.rebuild_queue(vert_map, face_map, cell_map);

  }

  //Deletes all the entries in the queue and empties it.
  void
  empty_queue()
  {
    while (pop_top_entry()) {
    }
    assert(is_empty());
  }

  void
  print_queue();

};

#endif
