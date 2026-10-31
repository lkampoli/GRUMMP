/*
 * Ghost.cxx
 *
 *  Created on: Aug 8, 2016
 *      Author: cfog
 */

#include "GR_config.h"

#ifdef HAVE_MPI

#include "GR_misc.h"
#include "GR_Ghost.h"
#include "GR_MPI_Defs.h"

using GRUMMP::Ghost;

Ghost::Ghost(const PartGID& localID, const GR_index_t nCopies) :
    m_ownerID(-1), m_localID(localID), m_copyID(nullptr), m_maxCopies(nCopies), m_copies(
	0), m_partBdryCopies(0)
{
  m_copyID = new EntityGID[m_maxCopies];
  assert(m_copies <= m_maxCopies);
  assert(m_partBdryCopies <= m_copies);
}

Ghost::Ghost(const Ghost& g) :
    m_ownerID(g.m_ownerID), m_localID(g.m_localID), m_copyID(
	new EntityGID[g.m_maxCopies]), m_maxCopies(g.m_maxCopies), m_copies(
	g.m_copies), m_partBdryCopies(g.m_partBdryCopies)
{
  for (GR_index_t ii = 0; ii < m_maxCopies; ii++) {
    m_copyID[ii] = g.m_copyID[ii];
  }
  assert(m_copies <= m_maxCopies);
  assert(m_partBdryCopies <= m_copies);
}

Ghost&
Ghost::operator=(const Ghost& g)
{
  if (&g != this) {
    m_localID = g.m_localID;
    m_ownerID = g.m_ownerID;
    if (m_maxCopies != g.m_maxCopies) {
      this->resetMaxNumCopies(g.m_maxCopies);
    }
    m_copies = g.m_copies;
    for (GR_index_t ii = 0; ii < m_maxCopies; ii++) {
      m_copyID[ii] = g.m_copyID[ii];
    }
    m_partBdryCopies = g.m_partBdryCopies;
  }
  assert(m_copies <= m_maxCopies);
  assert(m_partBdryCopies <= m_copies);
  return *this;
}

Ghost::Ghost(const Ghost& g, PartGID& target, PartGID& source, Entity* localEnt) :
    m_ownerID(target), m_localID(g.m_localID), m_copyID(
	new EntityGID[g.m_maxCopies]), m_maxCopies(g.m_maxCopies), m_copies(
	g.m_copies)
{
  for (GR_index_t ii = 0; ii < m_copies; ii++) {
    EntityGID thisID = g.m_copyID[ii];
#ifndef NDEBUG
    int numChanged = 0;
#endif
    if (thisID.getPart() == target) {
#ifndef NDEBUG
      numChanged++;
#endif
      m_copyID[ii] = EntityGID(source, localEnt);
    }
    else {
      m_copyID[ii] = thisID;
    }
#ifndef NDEBUG
    assert(numChanged <= 1);
#endif
  }
  m_partBdryCopies = g.m_partBdryCopies;
  assert(m_copies <= m_maxCopies);
  assert(m_partBdryCopies <= m_copies);
}

Ghost::~Ghost()
{
  if (m_copyID)
    delete[] m_copyID;
}

void
Ghost::resetMaxNumCopies(const GR_index_t newMaxNumCopies)
{
  EntityGID *tmpData = new EntityGID[newMaxNumCopies];
  for (GR_index_t ii = 0; ii < m_maxCopies; ii++) {
    tmpData[ii] = m_copyID[ii];
  }
  m_maxCopies = newMaxNumCopies;
  delete[] m_copyID;
  m_copyID = tmpData;
  assert(m_copies <= m_maxCopies);
  assert(m_partBdryCopies <= m_copies);
}

void
Ghost::setNumPartBdryCopies(const GR_index_t input)
{
  // Could be bigger than the number of copies we currently have.
  assert(input <= m_maxCopies);
  m_partBdryCopies = input;
}

bool
Ghost::addCopy(const EntityGID& egid, const bool isBdry)
{
  assert(m_copies <= m_maxCopies);
//    logMessage(1, "In addCopy; part %d\n", m_ownerID.getUniqueID());
  if (egid.getPart() == m_localID)
    return false;
  // There may be a copy with the correct part ID (null entity).  If
  // so, overwrite that entry.
  for (GR_index_t ii = 0; ii < m_copies; ++ii) {
    if (m_copyID[ii].getPart() == egid.getPart()) {
      m_copyID[ii] = egid;
      return true;
    }
  }
  // Otherwise, try to add a copy.
  if (m_copies == m_maxCopies) {
    resetMaxNumCopies((m_maxCopies + 1) * 3 / 2);
    assert(m_maxCopies > m_copies);
  }
  m_copyID[m_copies] = egid;
  m_copies++;
  assert(m_copies <= m_maxCopies);
  if (isBdry) {
    ++m_partBdryCopies;
    assert(m_partBdryCopies <= m_copies);
  }
  return true;
}

void
Ghost::reset()
{
  m_copies = m_partBdryCopies;
}

int
Ghost::merge(const Ghost& g, const bool isBdry)
{
  int origCopies = m_copies;
  for (GR_index_t ii = 0; ii < g.m_copies; ++ii) {
    EntityGID egid = g.getCopyID(ii);
    addCopy(egid, isBdry);
  }
  return (m_copies - origCopies);
}

void
Ghost::serialize(char buffer[], int& offset) const
{
  assert(m_copies <= m_maxCopies);
  transcribeToBuffer(buffer, offset, &m_copies, sizeof(GR_index_t));
  transcribeToBuffer(buffer, offset, m_copyID, sizeof(EntityGID) * m_copies);
  transcribeToBuffer(buffer, offset, &m_partBdryCopies, sizeof(GR_index_t));
}

void
Ghost::deserialize(const char* const buffer, int& offset)
{
  transcribeFromBuffer(buffer, offset, &m_copies, sizeof(GR_index_t));

  // Make sure there's enough space
  if (m_copies > m_maxCopies) {
    EntityGID *tmpData = new EntityGID[m_copies];
    m_maxCopies = m_copies;
    delete[] m_copyID;
    m_copyID = tmpData;
  }

  transcribeFromBuffer(buffer, offset, m_copyID, sizeof(EntityGID) * m_copies);
  assert(m_copies <= m_maxCopies);

  GR_index_t junk;
  transcribeFromBuffer(buffer, offset, &junk, sizeof(GR_index_t));
  assert(m_partBdryCopies <= m_copies);
}
#endif
