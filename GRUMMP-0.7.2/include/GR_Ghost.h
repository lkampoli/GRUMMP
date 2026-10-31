/*
 * GR_Ghost.h
 *
 *  Created on: Jul 27, 2016
 *      Author: cfog
 */

#ifndef GR_GHOST_H_
#define GR_GHOST_H_

#include "GR_config.h"

#ifdef HAVE_MPI

#include "GR_MPI_Defs.h"

namespace GRUMMP
{
  class Ghost {
    PartGID m_ownerID, m_localID;
    EntityGID *m_copyID;
    // m_maxCopies is the size of the array
    // m_copies is the number stored so far
    GR_index_t m_maxCopies, m_copies, m_partBdryCopies;
  public:
    // CFOG: I'm not entirely sure that the default value for the first arg
    // is safe.  Will add assertions to Mesh ghost retrieval functions...
    Ghost(const PartGID& localID = PartGID(), const GR_index_t nCopies = 1);
    Ghost(const Ghost&);

    // This constructor changes the owner part ID, and replaces the Entity
    // GID of the copy of this entity on target with the EGID of the entity
    // on source, in preparation for sending data back to other parts.
    Ghost(const Ghost& g, PartGID& target, PartGID& source, Entity* localEnt);

    Ghost&
    operator=(const Ghost&);

    virtual
    ~Ghost();
    const PartGID&
    getLocalID() const
    {
      return m_localID;
    }
    const PartGID&
    getOwnerID() const
    {
      return m_ownerID;
    }
    void
    setOwnerID(const PartGID& owner)
    {
      m_ownerID = owner;
    }
    const EntityGID&
    getCopyID(const GR_index_t copy = 0) const
    {
      assert(copy < m_maxCopies);
      return m_copyID[copy];
    }
    GR_index_t
    getMaxNumCopies() const
    {
      return m_maxCopies;
    }
    GR_index_t
    getNumCopies() const
    {
      return m_copies;
    }
    GR_index_t
    getNumPartBdryCopies() const
    {
      return m_partBdryCopies;
    }
    void
    setNumPartBdryCopies(const GR_index_t numPBCopies);
    void
    resetMaxNumCopies(const GR_index_t newMaxNumCopies);
    bool
    addCopy(const EntityGID& egid, const bool isBdry = false);

    // Remove all copy info except the owner, and a stub for non-owning
    // host.
    void
    reset();

    // Combine the contents of g into the contents of *this.
    int
    merge(const Ghost& g, const bool isBdry = false);

    // How much space will this Ghost take up in an MPI transport buffer?
    int
    MPIdataSize() const
    {
      return 2 * sizeof(int) + m_copies * sizeof(EntityGID);
    }
    // Write to a buffer for MPI transport; the buffer is allocated by the caller
    void
    serialize(char buffer[], int& offset) const;
    // Read from an MPI buffer; the buffer i
    void
    deserialize(const char* const buffer, int& offset);
  };
} // End namespace GRUMMP

#endif // HAVE_MPI

#endif /* GR_GHOST_H_ */
