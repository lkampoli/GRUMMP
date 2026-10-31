#include "GR_SamplingQueue.h"
#include "GR_Subseg.h"
#include "GR_SurfaceSampler.h"

#include <algorithm>
#include <map>
#include <queue>
#include <vector>

using std::map;
using std::vector;

bool
SubsegSamplingEntry::entry_deleted() const
{
  return m_bridge->get_parent()->is_deleted();
}
void
SubsegSamplingEntry::print_entry() const
{

  Subseg * ss = get_bridge()->get_parent();
  printf("ss \n");
  ss->get_beg_vert()->printVertInfo();
  ss->get_end_vert()->printVertInfo();
}
FaceSamplingEntry::FaceSamplingEntry(EntryType entry_type,
				     RefFace* const ref_face,
				     FaceQueueEntry* const fqe,
				     const double priority) :
    SamplingEntry(ref_face, priority), m_entry_type(entry_type), m_fqe(fqe), m_cell_left(
	fqe->getFace()->getLeftCell()), m_cell_right(
	fqe->getFace()->getRightCell())
{
  for (GR_index_t iV = 0; iV < fqe->getNumVerts(); iV++)
    assert(m_entry_type != SUBSEG);
  assert(!m_fqe->isDeleted());
  assert(!m_cell_left->isDeleted());
  assert(!m_cell_right->isDeleted());

}

bool
FaceSamplingEntry::entry_deleted() const
{

  bool qRetVal = ((m_fqe->isDeleted()
      || (m_cell_left != m_fqe->getFace()->getLeftCell())
      || (m_cell_right != m_fqe->getFace()->getRightCell())));
  if (qRetVal)
    return qRetVal;
  return !m_fqe->getFace()->isRestrictedDelaunay();

}
void
FaceSamplingEntry::print_entry() const
{
  printf("face %p %d\n", getFace(), getFace()->isRestrictedDelaunay());
  getFace()->printFaceInfo();

}
void
FaceTopoSamplingEntry::set_vertex(Vert* const vertex)
{

  m_offending_vert = vertex;
  assert(m_offending_vert == m_offending_vert_tag->get_vertex());

}

bool
FaceTopoSamplingEntry::entry_deleted() const
{

  if (m_fqe == NULL || m_cell_left == pCInvalidCell
      || m_cell_right == pCInvalidCell)
    return true;

  bool qRetVal = ((m_fqe->isDeleted()
      || (m_cell_left != m_fqe->getFace()->getLeftCell())
      || (m_cell_right != m_fqe->getFace()->getRightCell())));
  if (qRetVal)
    return qRetVal;
  return !m_fqe->getFace()->isRestrictedDelaunay();

}

void
SamplingPriorityQueue::clean_queue()
{

  if (empty())
    return;
  assert(std::__is_heap(c.begin(), c.end(), SamplingEntry::QueueOrder()));

#ifndef NDEBUG
  GR_index_t num_total = size();
#endif

  SamplingEntry* entry;
  EntryVec valid_entries, deleted_entries;
  EntryVec::iterator it = c.begin(), it_end = c.end();

  for (; it != it_end; ++it) {
    entry = *it;
    if (entry->remove_entry())
      deleted_entries.push_back(entry);
    else
      valid_entries.push_back(entry);
  }

#ifndef NDEBUG
  GR_index_t num_valid = valid_entries.size(), num_deleted =
      deleted_entries.size();
  assert(num_total == num_valid + num_deleted);
#endif

  c.swap(valid_entries);
  std::make_heap(c.begin(), c.end(), SamplingEntry::QueueOrder());
  assert(std::__is_heap(c.begin(), c.end(), SamplingEntry::QueueOrder()));

#ifndef NDEBUG
  assert(size() == num_valid);
#endif

  for (EntryVec::iterator it2 = deleted_entries.begin();
      it2 != deleted_entries.end(); ++it2) {
    if ((*it2))
      delete (*it2);
  }

}

void
SamplingPriorityQueue::rebuild_queue(const map<Vert*, Vert*>& vert_map,
				     const map<Face*, Face*>& face_map,
				     const map<Cell*, Cell*>& cell_map)
{

  SamplingEntry::EntryType entry_type;

  SamplingEntry *entry;
  FaceSamplingEntry *face_entry;
  FaceTopoSamplingEntry *face_topo_entry;

  map<Vert*, Vert*>::const_iterator itv;
  map<Face*, Face*>::const_iterator itf;
  map<Cell*, Cell*>::const_iterator itc1, itc2;

  vector<SamplingEntry*>::iterator it = c.begin(), it_end = c.end();

  for (; it != it_end; ++it) {

    entry = *it;
    entry_type = entry->get_entry_type();

    //nothing to change for a subseg entry
    if (entry_type == SamplingEntry::SUBSEG)
      continue;

    face_entry = dynamic_cast<FaceSamplingEntry*>(entry);
    assert(face_entry);
    itf = face_map.find(face_entry->getFace());
    itc1 = cell_map.find(face_entry->get_cell_left());
    itc2 = cell_map.find(face_entry->get_cell_right());

    //update the data in bad topo entries.
    if (entry_type == SamplingEntry::BAD_TOPO) {

      face_topo_entry = dynamic_cast<FaceTopoSamplingEntry*>(face_entry);
      assert(face_topo_entry);

      itv = vert_map.find(face_topo_entry->get_vertex());
      assert(itv != vert_map.end());
      face_topo_entry->set_vertex(itv->second);

      if (itf == face_map.end())
	face_topo_entry->set_face(pFInvalidFace);
      else
	face_topo_entry->set_face(itf->second);

      if (itc1 == cell_map.end())
	face_topo_entry->set_cell_left(pCInvalidCell);
      else
	face_topo_entry->set_cell_left(itc1->second);
      if (itc2 == cell_map.end())
	face_topo_entry->set_cell_right(pCInvalidCell);
      else
	face_topo_entry->set_cell_right(itc2->second);

    }

    //update all other face sampling entries
    else {

      assert(itf != face_map.end());
      face_entry->set_face(itf->second);

      assert(itc1 != cell_map.end());
      face_entry->set_cell_left(itc1->second);
      assert(itc2 != cell_map.end());
      face_entry->set_cell_right(itc2->second);

      assert(!face_entry->entry_deleted());

    }

#ifndef NDEBUG
    if (!face_entry->entry_deleted()) {

      assert(face_entry->getFace()->isValid());
      assert(!face_entry->getFace()->isDeleted());

      assert(face_entry->get_cell_left()->isValid());
      assert(face_entry->get_cell_right()->isValid());
      assert(!face_entry->get_cell_left()->isDeleted());
      assert(!face_entry->get_cell_right()->isDeleted());

      assert(
	  face_entry->get_cell_left() == face_entry->getFace()->getLeftCell());
      assert(
	  face_entry->get_cell_right()
	      == face_entry->getFace()->getRightCell());

    }
#endif

  }

}

SamplingEntry*
SamplingQueue::top_entry() const
{

  if (!m_multiple.empty())
    return (m_multiple.top());
  if (!m_subseg.empty())
    return (m_subseg.top());
  if (!m_bad_topo.empty())
    return (m_bad_topo.top());
  if (!m_size.empty())
    return (m_size.top());
  if (!m_quality.empty())
    return (m_quality.top());

  return static_cast<SamplingEntry*>(NULL);

}

bool
SamplingQueue::pop_top_entry()
{

  if (!m_multiple.empty()) {
    SamplingEntry* dummy = m_multiple.top();
    m_multiple.pop();
    if (dummy)
      delete dummy;
    return true;
  }
  if (!m_subseg.empty()) {
    SamplingEntry* dummy = m_subseg.top();
    m_subseg.pop();
    if (dummy)
      delete dummy;
    return true;
  }
  if (!m_bad_topo.empty()) {
    SamplingEntry* dummy = m_bad_topo.top();
    m_bad_topo.pop();
    if (dummy)
      delete dummy;
    return true;
  }
  if (!m_size.empty()) {
    SamplingEntry* dummy = m_size.top();
    m_size.pop();
    if (dummy)
      delete dummy;
    return true;
  }
  if (!m_quality.empty()) {
    SamplingEntry* dummy = m_quality.top();
    m_quality.pop();
    if (dummy)
      delete dummy;
    return true;
  }

  return false;

}

void
SamplingQueue::add_to_queue(SamplingEntry* const sampling_entry)
{

  switch (sampling_entry->get_entry_type())
    {

    case SamplingEntry::MULTIPLE:
      m_multiple.push(sampling_entry);
      break;
    case SamplingEntry::SUBSEG:
      m_subseg.push(sampling_entry);
      break;
    case SamplingEntry::BAD_TOPO:
      m_bad_topo.push(sampling_entry);
      break;
    case SamplingEntry::SIZE:
      m_size.push(sampling_entry);
      break;
    case SamplingEntry::QUALITY:
      m_quality.push(sampling_entry);
      break;
    default:
      vFatalError("Unknown sampling entry type",
		  "SamplingQueue::add_to_queue()");
      break;

    }

}

void
SamplingQueue::print_queue()
{

  SamplingEntry* entry = top_entry();
  printf("Printing Sampling Queue\n");
  while (entry) {
    printf("type = %d, priority = %lf, deleted = %d\n", entry->get_entry_type(),
	   entry->get_priority(), entry->entry_deleted());
    if (entry->get_entry_type() != SamplingEntry::SUBSEG) {
      FaceSamplingEntry* face_entry = dynamic_cast<FaceSamplingEntry*>(entry);
      printf("face %d %d %d\n", face_entry->getFaceQueueEntry()->isDeleted(),
	     face_entry->getFace()->isDeleted(),
	     face_entry->getFace()->isRestrictedDelaunay());
      face_entry->getFace()->printFaceInfo();

    }
    else {
      SubsegSamplingEntry* ss_entry = dynamic_cast<SubsegSamplingEntry*>(entry);
      Subseg * ss = ss_entry->get_bridge()->get_parent();
      printf("ss \n");
      ss->get_beg_vert()->printVertInfo();
      ss->get_end_vert()->printVertInfo();
    }
    pop_top_entry();
    entry = top_entry();
  }

  assert(0);

}
