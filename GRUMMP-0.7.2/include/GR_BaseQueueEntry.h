/*
 * GR_BaseQueueEntry.h
 *
 *  Created on: 2013-09-13
 *      Author: zaide
 *      Serves to hold information about objects (cells or faces) and check
 *      if they are still in the mesh. Meant to go in lots of places.
 */

#ifndef GR_BASEQUEUEENTRY_H_
#define GR_BASEQUEUEENTRY_H_
#include "GR_config.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_BFace.h"
#include <map>

/// Currently mostly contains FaceQueueEntry and its Comparator.
/// Used in TetMeshBuilder for curved surfaces, stripped down version
/// of insertion queue. Challenge is that a lot of searching needs to be done,
/// hence the more complicated comparator

/// Also contains CellQueueEntry used in RefinementManager3D. It acts like a
/// lean version of InsertionQueueEntry
class FaceQueueEntry {

private:
  Face * m_face;
  Vert ** m_verts;

  Vert **
  getVerts() const
  {
    return m_verts;
  }
public:
  friend struct FaceQueueEntryComp;
  FaceQueueEntry() :
      m_face(NULL), m_verts(NULL)
  {
    ;
  }
  FaceQueueEntry(Face * face) :
      m_face(face)
  {
    assert(face->getNumVerts() > 0);
    m_verts = new Vert *[face->getNumVerts()];
    for (int iV = 0; iV < face->getNumVerts(); iV++)
      m_verts[iV] = m_face->getVert(iV);
  }
  ~FaceQueueEntry()
  {
    delete[] m_verts;
  }
  void
  setFace(Face * face)
  {
    m_face = face;
    reset();
  }
  bool
  isDeleted() const
  {
    if (m_face->isDeleted())
      return true;
    for (int iV = 0; iV < m_face->getNumVerts(); iV++)
      if (!m_face->hasVert(m_verts[iV]))
	return true;

    return false;
  }
  bool
  isNotDeleted() const
  {
    return !isDeleted();
  }
  bool
  isEqual(FaceQueueEntry* entry) const
  {
    if (entry == NULL)
      return false;
    if (m_face != entry->m_face)
      return false;
    for (int iV = 0; iV < m_face->getNumVerts(); iV++)
      if (!entry->hasVert(m_verts[iV]))
	return false;
    return true;
  }

  Face *
  getFace() const
  {
    return m_face;
  }
  void
  reset() const
  {
    for (int iV = 0; iV < m_face->getNumVerts(); iV++)
      m_verts[iV] = m_face->getVert(iV);
  }
  GR_index_t
  getNumVerts() const
  {
    assert(m_face->getNumVerts() > 0);
    return m_face->getNumVerts();
  }
  bool
  isValid()
  {
    return m_face->isValid();
  }
  Vert *
  getVert(const int i) const
  {
    assert(i >= 0 && i < m_face->getNumVerts());
    assert(m_face->getVert(i));
    return (m_face->getVert(i));
  }
//	Vert *getStoredVert(const int i) const
//	{
//		assert(i >= 0 && i < m_face->getNumVerts());
//		return (m_verts[i]);
//	}
  bool
  hasVert(const Vert * vert) const
  {
    return m_face->hasVert(vert);
  }

};

struct FaceQueueEntryComp {
  bool
  operator()(const FaceQueueEntry * fqe1, const FaceQueueEntry * fqe2) const
  {

    if (fqe1->getFace() != fqe2->getFace()) {
      return (fqe1->getFace() < fqe2->getFace());
    }
    else {

      // Puts deleted greater, so that if the faces are equal, deleted facequeueentries show up after
      if (fqe1->isDeleted() || fqe2->isDeleted())
	return fqe1->isDeleted() > fqe2->isDeleted();

      Vert * temp_vert;
      Vert ** verts1 = fqe1->getVerts();
      Vert ** verts2 = fqe2->getVerts();
      assert(fqe1->getNumVerts() == 3);
      // sort verts1
      if (verts1[0] > verts1[1]) {
	temp_vert = verts1[1];
	verts1[1] = verts1[0];
	verts1[0] = temp_vert;
      }
      if (verts1[1] > verts1[2]) {
	temp_vert = verts1[1];
	verts1[1] = verts1[2];
	verts1[2] = temp_vert;
      }
      if (verts1[0] > verts1[1]) {
	temp_vert = verts1[1];
	verts1[1] = verts1[0];
	verts1[0] = temp_vert;
      }
      // sort verts2
      if (verts2[0] > verts2[1]) {
	temp_vert = verts2[1];
	verts2[1] = verts2[0];
	verts2[0] = temp_vert;
      }
      if (verts2[1] > verts2[2]) {
	temp_vert = verts2[1];
	verts2[1] = verts2[2];
	verts2[2] = temp_vert;
      }
      if (verts2[0] > verts2[1]) {
	temp_vert = verts2[1];
	verts2[1] = verts2[0];
	verts2[0] = temp_vert;
      }

      return verts1[0] < verts2[0] || verts1[1] < verts2[1]
	  || verts1[2] < verts2[2];
    }
  }
};

typedef std::set<FaceQueueEntry*, FaceQueueEntryComp> FaceQueueEntrySet;

class CellQueueEntry {

  Cell* m_cell;
  double m_priority;
  Vert* m_verts[4];
  CellQueueEntry() :
      m_cell(NULL), m_priority(LARGE_DBL)
  {
    assert(0);
  }
  CellQueueEntry(const CellQueueEntry&) :
      m_cell(NULL), m_priority(LARGE_DBL)
  {
    assert(0);
  }
  CellQueueEntry&
  operator=(const CellQueueEntry&)
  {
    assert(0);
    return *this;
  }

public:

  CellQueueEntry(Cell* const cell, double priority = LARGE_DBL) :
      m_cell(cell), m_priority(priority)
  {
    for (int iV = 0; iV < 4; iV++)
      m_verts[iV] = m_cell->getVert(iV);
  }

  virtual
  ~CellQueueEntry()
  {
  }

  Cell*
  get_cell() const
  {
    assert(m_cell);
    return m_cell;
  }
  Vert*
  get_vert(int i)
  {
    assert(m_verts);
    assert(i >= 0 && i < 4);
    return m_verts[i];
  }
  double
  get_priority() const
  {
    assert(m_cell);
    return m_priority;
  }

  void
  set_cell(Cell* cell)
  {
    m_cell = cell;
    for (int iV = 0; iV < 4; iV++)
      m_verts[iV] = m_cell->getVert(iV);
  }

  bool
  isDeleted() const
  {
    assert(m_cell);
    if (m_cell->isDeleted())
      return true;
    for (int iV = 0; iV < m_cell->getNumVerts(); iV++) {
      if (m_verts[iV] != m_cell->getVert(iV))
	return true;
    }
    return false;
  }
  bool
  is_equal(CellQueueEntry * const entry)
  {
    if (entry == this)
      return true;
    if (entry == NULL)
      return false;
    if (m_cell != entry->get_cell())
      return false;
    for (int iV = 0; iV < m_cell->getNumVerts(); iV++)
      if (!entry->get_cell()->hasVert(m_verts[iV]))
	return false;
    return true;
  }
  class PrioOrder {
    bool m_reverse;
  public:
    PrioOrder(const bool& reverse = false) :
	m_reverse(reverse)
    {
    }
    bool
    max_on_top() const
    {
      return !m_reverse;
    }
    bool
    operator()(const CellQueueEntry* const entry0,
	       const CellQueueEntry* const entry1) const
    {
      int iCond = iFuzzyComp(entry0->m_priority, entry1->m_priority);

      if (m_reverse) {
	if (iCond == 1)
	  return true;
	else if (iCond == -1)
	  return false;
	else
	  return entry0 > entry1;
      }
      else {
	if (iCond == 1)
	  return false;
	else if (iCond == -1)
	  return true;
	else
	  return entry0 < entry1;
      }
    }
  };

};
#endif /* GR_BASEQUEUEENTRY_H_ */
