#ifndef GR_Partition
#define GR_Partition 1

#include "GR_config.h"
#ifdef HAVE_MPI

#include <mpi.h>
#include <set>
#include <vector>

#include "GR_Mesh.h"
#include "GR_MPI_Defs.h"

namespace GRUMMP
{
  class Partition {
  private:
    int m_numLocalParts;
    GR_index_t m_numGlobalParts;
    int m_rank;
    MPI::Intercomm m_comm;
    std::vector<Mesh*> m_meshes;

    Partition(const Partition&);
    Partition&
    operator=(const Partition&);
  public:
    Partition(const int MeshType, const int numParts = 1,
	      MPI::Intercomm communicator = MPI::COMM_WORLD);
    ~Partition();
    int
    getNumLocalParts() const
    {
      return m_numLocalParts;
    }
    int
    getNumGlobalParts() const
    {
      return m_numGlobalParts;
    }
    Mesh*
    getPart(const int i)
    {
      if (i < 0 || i >= m_numLocalParts)
	return NULL;
      else
	return m_meshes[i];
    }
    PartGID
    getPartGID(PartLID_t localID)
    {
      return PartGID(m_rank, localID);
    }
    // Do any behind the scenes work needed so that global partition queries and
    // remote handles are up to date (and therefore accessible without further
    // communication).
    void
    sync();

    //! How many parts have entities shared with this one?
    //
    // Communication: None (pre-computed by sync).
    GR_index_t
    getNumPartNeighbors(const int partIndex);

    //! What parts have entities shared with this one?
    //
    // Communication: None (pre-computed by sync).
    void
    getPartNeighbors(std::vector<Mesh*>& partList);

    // Part bdry entities: number and identity
    // Entities in a part: via Mesh member functions
    // Adjacencies: via entity member functions; but need to have ghosts!
    // Get entity owner
    // Do I own this entity?
    // Entity status (internal, part bdry, ghost)
    // Number of copies of an entity (owned: accurate; ghost: 2)
    //

    void
    updateGhosts(Mesh* pM);

//    void sendGhostInfo(const PartGID targetID, const GhostInfo GI);
//    void receiveGhostInfo(const PartGID& myID, const int mesg);
//    void relayGhostInfo(const int mesg);
//    void receiveRelayedGhostInfo(const int mesg);

//    void getAllGhostData(const GhostInfo& GI, int& numVerts, double*& coords,
//                         int& numEntities, int*& offsets, int& numIndices,
//                         int*& indices, Entity**& handles);
  public:
    void
    sendVertexData(const GhostMsgTags_t messageTag, const PartGID& target,
		   const int numVerts, const double coords[], Vert* handles[]);
    void
    receiveVertexData(const GhostMsgTags_t messageTag, PartGID& senderID,
		      int& numVerts, double*& coords, Vert**& handles);
    void
    sendConnectivityData(const GhostMsgTags_t messageTag,
			 const PartGID& sourceID, const PartGID& targetID,
			 const int numEntities, const Entity* handles[],
			 const int offsets[], const int numVerts,
			 const Vert* verts[]);
    void
    receiveConnectivityData(const GhostMsgTags_t messageTag, PartGID& sourceID,
			    int& numEntities, Entity**& handles, int*& offsets,
			    int& numVerts, Vert**& verts);
    void
    receiveCopyData(const GhostMsgTags_t messageTag, PartGID& senderID,
		    Mesh* pMesh);
    const MPI::Intercomm&
    getCommunicator() const
    {
      return m_comm;
    }
  };

#endif

} // End namespace GRUMMP
#endif

