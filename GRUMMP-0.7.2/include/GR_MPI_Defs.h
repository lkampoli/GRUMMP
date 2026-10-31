/*
 * GR_MPI_Defs.h
 *
 *  Created on: Jul 28, 2016
 *      Author: cfog
 */

#ifndef GR_MPI_DEFS_H_
#define GR_MPI_DEFS_H_

#include "GR_config.h"

#ifdef HAVE_MPI

#include <limits>

namespace GRUMMP
{
  static const int PART_BITS = 12;
  static const int PROC_BITS = sizeof(int) * 8 - PART_BITS;

// Default values for Process ID, Part local ID and ENTITY local ID
  static const int PROC_UNDEFINED = -2;
  static const int PROC_ANY = -3;
  static const int PROC_ALL = -4;

  static const int PART_UNDEFINED = 0;
  static const int PART_ALL = PART_UNDEFINED - 1;
  static const int PART_ANY = PART_UNDEFINED - 2;

//  static const int ENTITY_UNDEFINED =
//      std::numeric_limits<int>::max();
//  static const int ENTITY_ALL = ENTITY_UNDEFINED-1;
//  static const int ENTITY_ANY = ENTITY_UNDEFINED-2;

  enum EntityStatus_t {
    eInt, eExt, eBdryInt, eBdryExt
  };
  enum EntType_t {
    eVert, eFace, eCell, eAll, eAny
  };

  enum GhostMsgTags_t {
    PartBdryVertCoords = 32768,
    PartBdryVertCopies,
    PartBdryFaces,
    GhostVertCoords,
    GhostVertCopies,
    GhostEntities
  };

  typedef GR_index_t PartLID_t;
  typedef GR_index_t EntityLID_t;

// Global handle for a mesh part. Valid on all processes and parts
  class PartGID {
    union {
      struct {
	GR_sindex_t iProcID :PROC_BITS;
	PartLID_t iPartLID :PART_BITS;
      };
      unsigned int iUnique;
    };
  public:
    PartGID(GR_sindex_t _iRank = PROC_UNDEFINED, PartLID_t _iPartLID =
		PART_UNDEFINED) :
	iProcID(_iRank), iPartLID(_iPartLID)
    {
    }
    ;

    PartGID(const PartGID& other) :
	iProcID(other.iProcID), iPartLID(other.iPartLID)
    {
    }

    bool
    operator ==(const PartGID &other) const
    {
      return (iProcID == other.iProcID && iPartLID == other.iPartLID);
    }

    bool
    operator !=(const PartGID& other) const
    {
      return !(*this == other);
    }

    bool
    operator <(const PartGID &other) const
    {
      if (iProcID < other.iProcID)
	return true;
      if (iProcID > other.iProcID)
	return false;
      return (iPartLID < other.iPartLID);
    }
    int
    getPartLocalID() const
    {
      return iPartLID;
    }
    int
    getRank() const
    {
      return iProcID;
    }
    unsigned int
    getUniqueID() const
    {
      return iUnique;
    }
  };

  class EntityGID {
    Entity *m_ent;
    PartGID m_part;
  public:
    EntityGID(PartGID part = PartGID(), Entity *ent = nullptr) :
	m_ent(ent), m_part(part)
    {
    }

    const Entity*
    getEnt() const
    {
      return m_ent;
    }

    void
    setEnt(Entity* const ent)
    {
      m_ent = ent;
    }

    const PartGID&
    getPart() const
    {
      return m_part;
    }

    void
    setPart(const PartGID& part)
    {
      m_part = part;
    }
    bool
    operator==(const EntityGID& other) const
    {
      return (m_ent == other.m_ent && m_part == other.m_part);
    }
    bool
    operator<(const EntityGID& other) const
    {
      return (m_part < other.m_part
	  || (m_part == other.m_part && m_ent < other.m_ent));
    }
  };

//! Write data to an MPI buffer
//
// buffer [inout]   Data is written to here; memory handled by caller
// offset [inout]   Location within buffer to write; updated to point after
//                     data that is written
// data [in]        Stuff to be transcribed; memory handled by caller
// blockSize [in]   Amount of stuff to be transcribed
  void
  transcribeToBuffer(char* const buffer, int& offset, const void* const data,
		     const size_t blockSize);

//! Read data from an MPI buffer
//
// buffer [in]      Data is read from here; memory handled by caller
// offset [inout]   Location within buffer to read; updated to point after
//                     data that is read
// data [inout]     Location to copy data to; memory handled by caller
// blockSize [in]   Amount of stuff to be transcribed
  void
  transcribeFromBuffer(const char* const buffer, int& offset, void* const data,
		       const size_t blockSize);
}

#endif // HAVE_MPI

#endif /* GR_MPI_DEFS_H_ */
