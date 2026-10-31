/*
 * Partition.cxx
 *
 *  Created on: Jul 28, 2016
 *      Author: cfog
 */

#include "GR_config.h"

#ifdef HAVE_MPI

#include <map>
#include <set>
#include "GR_ADT.h"
#include "GR_Partition.h"
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include "GR_SurfMesh.h"
#include "GR_VolMesh.h"

using std::map;
using std::set;
using namespace GRUMMP;

Partition::Partition(const int MeshType, const int numParts,
		     MPI::Intercomm communicator) :
    m_numLocalParts(numParts), m_numGlobalParts(-1), m_rank(-1), m_comm(
	communicator), m_meshes(numParts)
{
  m_rank = m_comm.Get_rank();
  for (int ii = 0; ii < m_numLocalParts; ii++) {
    switch (MeshType)
      {
      default: // nobreak
      case Mesh::eMesh2D:
	m_meshes[ii] = new Mesh2D;
	break;
      case Mesh::eSurfMesh:
	m_meshes[ii] = new SurfMesh;
	break;
      case Mesh::eVolMesh:
	m_meshes[ii] = new VolMesh;
	break;
      }
    m_meshes[ii]->setPartition(this);
    m_meshes[ii]->setPartID(m_rank, ii);
  }
  m_comm.Allreduce(&m_numLocalParts, &m_numGlobalParts, 1, MPI::INT, MPI::SUM);
}

Partition::~Partition()
{
  for (int ii = 0; ii < m_numLocalParts; ii++) {
    delete m_meshes[ii];
  }
}

class GhostEnts {
  typedef std::map<PartGID, std::set<Entity*> > map_t;
  map_t m_data;
public:
  typedef std::map<PartGID, std::set<Entity*> >::iterator GEiter;
  typedef std::map<PartGID, std::set<Entity*> >::const_iterator GE_const_iter;
  GhostEnts() :
      m_data()
  {
  }
  ~GhostEnts()
  {
  }
  size_t
  size() const
  {
    return m_data.size();
  }
  bool
  addEntity(Entity* newEnt, const PartGID& target)
  {
    assert(newEnt);
    map_t::iterator iter = m_data.find(target);
    if (iter == m_data.end()) {
      std::set<Entity*> tmp;
      std::pair<map_t::iterator, bool> result = m_data.insert(
	  std::make_pair(target, tmp));
      if (!result.second)
	return false;
      iter = result.first;
    }
    return (iter->second.insert(newEnt)).second;
  }
  GE_const_iter
  begin() const
  {
    return m_data.begin();
  }
  GE_const_iter
  end() const
  {
    return m_data.end();
  }
};

void
Mesh::sendVertInfo(const enum GhostMsgTags_t tag, const PartGID& targetID,
		   const std::set<Entity*>& verts)
{
  PartGID sourceID = getPartID();
  logMessage(3, "Should send a message from part %x to part %x; %lu verts\n",
	     sourceID.getUniqueID(), targetID.getUniqueID(), verts.size());

  std::set<Entity*>::const_iterator iter = verts.begin();
  const std::set<Entity*>::const_iterator iterEnd = verts.end();

  // Sending only verts, so the number of other entities is zero.
  int numVerts = verts.size();
  double *coords = new double[numVerts * 3];
  Vert **handles = new Vert*[numVerts];
  int iV = 0;
  for (; iter != iterEnd; ++iter, iV++) {
    Vert* vert = dynamic_cast<Vert*>(*iter);
    assert(vert);
    coords[iV * 3 + 0] = vert->x();
    coords[iV * 3 + 1] = vert->y();
    coords[iV * 3 + 2] = vert->z();
    handles[iV] = vert;
    logMessage(4, "Sending vert %p, coords (%f, %f, %f)\n", vert, vert->x(),
	       vert->y(), vert->z());
  }

  Partition* ptn = getPartition();
  ptn->sendVertexData(tag, targetID, numVerts, coords, handles);

  // Must delete memory that was allocated elsewhere.
  if (coords)
    delete[] coords;
  if (handles)
    delete[] handles;
}

// This function packages up connectivity info for transport; the actual
// send is done by a partition member function, including handling the data
// buffer.
static void
sendConnectInfo(const enum GhostMsgTags_t tag, const PartGID& targetID,
		Mesh* const pMesh, const std::set<Entity*>& ents)
{
  PartGID sourceID = pMesh->getPartID();
  logMessage(3, "Should send a message from part %x to part %x\n",
	     sourceID.getUniqueID(), targetID.getUniqueID());

  std::set<Entity*>::const_iterator iter = ents.begin();
  const std::set<Entity*>::const_iterator iterEnd = ents.end();

  // Sending only faces.
  int numEnts = ents.size();
  const Entity** handles = new const Entity*[numEnts];
  int *offsets = new int[numEnts + 1];

  // Two passes through the data: one to count the verts, another to put them
  // into the array (after it's been created).
  int numVerts = 0, ii = 0;
  for (; iter != iterEnd; ++iter, ++ii) {
    Entity *pEnt = *iter;
    handles[ii] = pEnt;
    offsets[ii] = numVerts;
    numVerts += pEnt->getNumVerts();
  }
  assert(ii == numEnts);
  offsets[ii] = numVerts; // The end of the last list of vert connectivity

  Vert const **verts = new const Vert*[numVerts];
  // Now again, to transcribe the verts.
  ii = 0;
  for (iter = ents.begin(); iter != iterEnd; ++iter) {
    Entity *pEnt = *iter;
    for (int iV = 0; iV < pEnt->getNumVerts(); ++iV, ++ii) {
      Vert *pVLocal = pEnt->getVert(iV);
      bool success = false;
      verts[ii] = nullptr;
      Ghost g = pMesh->getGhostEntry(pVLocal, success);
      for (GR_index_t iCopy = 0; iCopy < g.getNumCopies(); ++iCopy) {
	EntityGID egid = g.getCopyID(iCopy);
	logMessage(4, "Part: %d.  Copy #%d.  EGID: (%d, %p)\n",
		   sourceID.getUniqueID(), iCopy, egid.getPart().getUniqueID(),
		   egid.getEnt());
	if (egid.getPart() == targetID) {
	  // These entity handles aren't necessarily valid pointers,
	  // nor necessarily in the process address space if they
	  // were, so dynamic_cast can fail if the system is checking
	  // whether it's a valid pointer.  reinterpret_cast = "These
	  // bits means what I say they mean because I say so" is
	  // a pretty risky thing to do, unless there's no other
	  // choice.
	  verts[ii] = reinterpret_cast<const Vert*>(egid.getEnt());
	  break; // No point finishing the loop.
	}
      }
      assert(verts[ii] != nullptr);
    }
  }
  assert(ii == numVerts);

  Partition* ptn = pMesh->getPartition();
  ptn->sendConnectivityData(tag, sourceID, targetID, numEnts, handles, offsets,
			    numVerts, verts);

  if (handles)
    delete[] handles;
  if (offsets)
    delete[] offsets;
  if (verts)
    delete[] verts;
}

static void
receivePBFaceInfo(Mesh * const pMesh)
{
  Partition *ptn = pMesh->getPartition();
  int numEntities = 0, numVerts = 0, *offsets = nullptr;
  Vert ** verts = nullptr;
  Entity ** handles = nullptr;
  PartGID source;

  ptn->receiveConnectivityData(PartBdryFaces, source, numEntities, handles,
			       offsets, numVerts, verts);

  // The faces are foreign handles; the verts are local handles.
  for (int ii = 0; ii < numEntities; ++ii) {
    int numVertsThisFace = offsets[ii + 1] - offsets[ii];
    assert(numVertsThisFace >= 2 && numVertsThisFace <= 4);
    Face *pF = nullptr;
    switch (numVertsThisFace)
      {
      case 2:
	pF = findCommonFace(verts[offsets[ii]], verts[offsets[ii] + 1]);
	break;
      case 3:
	pF = findCommonFace(verts[offsets[ii]], verts[offsets[ii] + 1],
			    verts[offsets[ii] + 2]);
	break;
      case 4:
	pF = findCommonFace(verts[offsets[ii]], verts[offsets[ii] + 1],
			    verts[offsets[ii] + 2], verts[offsets[ii] + 3]);
	break;
      default:
	assert(0);
	break;
      }
    assert(pF);
    EntityGID egid(source, handles[ii]);
    bool success = pMesh->addCopyToGhostEntry(pF, egid);
    assert(success);
  }

  logMessage(3, "Part %d: Successfully processed %d entities from part %d\n",
	     pMesh->getPartID().getUniqueID(), numEntities,
	     source.getUniqueID());
  delete[] offsets;
  delete[] verts;
  delete[] handles;
}

static void
receiveCellInfo(Mesh * const pMesh)
{
  Partition *ptn = pMesh->getPartition();
  // Is this a topologically 2D or 3D mesh?
  bool is3D = (pMesh->getType() == Mesh::eVolMesh);

  int numEntities = 0, numVerts = 0, *offsets = nullptr;
  Vert ** verts = nullptr;
  Entity ** handles = nullptr;
  PartGID source;

  ptn->receiveConnectivityData(GhostEntities, source, numEntities, handles,
			       offsets, numVerts, verts);

  // The faces are foreign handles; the verts are local handles.
  for (int ii = 0; ii < numEntities; ++ii) {
    int numVertsThisCell = offsets[ii + 1] - offsets[ii];
    assert(numVertsThisCell >= 2 && numVertsThisCell <= 4);
    Cell *pC = nullptr;
    switch (numVertsThisCell)
      {
      case 3:
	{
	  // This is a triangle
	  Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pMesh);
	  pC = pM2D->createTriCell(verts[offsets[ii]], verts[offsets[ii] + 1],
				   verts[offsets[ii] + 2]);

	  break;
	}
      case 4:
	// Could be a quad or a tet.
	if (is3D) {
	  // It's a tet.
	  VolMesh *pVolMesh = dynamic_cast<VolMesh*>(pMesh);
	  bool qAlreadyExisted = false;
	  pC = pVolMesh->createTetCell(qAlreadyExisted, verts[offsets[ii]],
				       verts[offsets[ii] + 1],
				       verts[offsets[ii] + 2],
				       verts[offsets[ii] + 3]);
	  assert(!qAlreadyExisted);
	}
	else {
	  Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pMesh);
	  pC = pM2D->createQuadCell(verts[offsets[ii]], verts[offsets[ii] + 1],
				    verts[offsets[ii] + 2],
				    verts[offsets[ii] + 3]);
	}
	break;
      case 5:
	{
	  // It's a pyramid
	  VolMesh *pVolMesh = dynamic_cast<VolMesh*>(pMesh);
	  bool qAlreadyExisted = false;
	  pC = pVolMesh->createPyrCell(qAlreadyExisted, verts[offsets[ii]],
				       verts[offsets[ii] + 1],
				       verts[offsets[ii] + 2],
				       verts[offsets[ii] + 3],
				       verts[offsets[ii] + 4]);
	  assert(!qAlreadyExisted);
	  break;
	}
      case 6:
	{
	  // It's a prism
	  VolMesh *pVolMesh = dynamic_cast<VolMesh*>(pMesh);
	  bool qAlreadyExisted = false;
	  pC = pVolMesh->createPrismCell(qAlreadyExisted, verts[offsets[ii]],
					 verts[offsets[ii] + 1],
					 verts[offsets[ii] + 2],
					 verts[offsets[ii] + 3],
					 verts[offsets[ii] + 4],
					 verts[offsets[ii] + 5]);
	  assert(!qAlreadyExisted);
	  break;
	}
      case 8:
	{
	  // It's a hex
	  VolMesh *pVolMesh = dynamic_cast<VolMesh*>(pMesh);
	  bool qAlreadyExisted = false;
	  pC = pVolMesh->createHexCell(qAlreadyExisted, verts[offsets[ii]],
				       verts[offsets[ii] + 1],
				       verts[offsets[ii] + 2],
				       verts[offsets[ii] + 3],
				       verts[offsets[ii] + 4],
				       verts[offsets[ii] + 5],
				       verts[offsets[ii] + 6],
				       verts[offsets[ii] + 7]);
	  assert(!qAlreadyExisted);
	  break;
	}
      default:
	assert(0);
	break;
      }
    assert(pC);
    pC->setEntStatus(Entity::eForeign);
    EntityGID egid(source, handles[ii]);
    bool success = pMesh->addCopyToGhostEntry(pC, egid);
    assert(success);
  }

  // Now tag faces that have no non-foreign cells adjacent to them
  // as foreign.
  for (GR_index_t ii = 0; ii < pMesh->getNumFaces(); ++ii) {
    Face *pF = pMesh->getFace(ii);
    if (pF->isDeleted())
      continue;
    bool isLeftForeign(true), isRightForeign(true);
    Cell *pC = pF->getLeftCell();
    if (pC != pCInvalidCell) {
      isLeftForeign = (pC->getEntStatus() == Entity::eForeign);
    }
    pC = pF->getRightCell();
    if (pC != pCInvalidCell) {
      isRightForeign = (pC->getEntStatus() == Entity::eForeign);
    }
    if (isLeftForeign && isRightForeign) {
      pF->setEntStatus(Entity::eForeign);
    }
  }

  logMessage(3, "Part %d: Successfully processed %d entities from part %d\n",
	     pMesh->getPartID().getUniqueID(), numEntities,
	     source.getUniqueID());
  delete[] offsets;
  delete[] verts;
  delete[] handles;
}

static int
receiveGhostVertMessage(const enum GhostMsgTags_t tag, Mesh* const pMesh,
			ADT& tree, std::set<Vert*>& vertsReceived)
{
  bool isBdry = (tag == PartBdryVertCoords);
  Partition *ptn = pMesh->getPartition();
  int numVerts = 0;
  double *coords = nullptr;
  Vert **handles = nullptr;
  PartGID senderID, myID = pMesh->getPartID();

  // Grab the actual message.
  ptn->receiveVertexData(tag, senderID, numVerts, coords, handles);
  logMessage(3, "On part %d, got %d verts. %p\n", myID.getUniqueID(), numVerts,
	     coords);

  int matchedVerts = 0;
  for (int ii = 0; ii < numVerts; ii++) {
    logMessage(4, "Part %d.  Coords: (%f, %f, %f)\n", myID.getUniqueID(),
	       coords[ii * 3], coords[ii * 3 + 1], coords[ii * 3 + 2]);
    Vert* localVert = findVertex(pMesh, tree, coords + ii * 3);
    if (localVert) {
      logMessage(4, "Part %d.  Found vert %p.   Coords: (%f, %f, %f)\n",
		 myID.getUniqueID(), localVert, localVert->x(), localVert->y(),
		 localVert->z());
      assert(localVert);
      assert(localVert->x() == *(coords + ii * 3));
      assert(localVert->y() == *(coords + ii * 3 + 1));
      matchedVerts++;
    }
    else {
      // Create a new vertex
      localVert = pMesh->createVert(coords + ii * 3);
      localVert->setEntStatus(Entity::eForeign);
      // Setting the owner as sender is potentially incorrect:  this
      // could actually be a ghost copy, especially on a part bdry.  It's
      // possible that this will matter someday.
      pMesh->addGhostEntry(localVert, senderID);
      tree.vAddNode(coords + ii * 3, pMesh->getNumVerts() - 1);
    }
    EntityGID remoteEGID(senderID, handles[ii]);
    vertsReceived.insert(localVert);
    bool success = pMesh->addCopyToGhostEntry(localVert, remoteEGID, isBdry);
    assert(success);
  }
  logMessage(3, "On part %d, matched %d verts out of %d.\n", myID.getUniqueID(),
	     matchedVerts, numVerts);

  if (coords)
    delete[] coords;
  if (handles)
    delete[] handles;
  return numVerts;
}

void
Mesh::replyWithVertIDs(const enum GhostMsgTags_t tag,
		       const std::set<Vert*>& vertsReceived)
{
  PartGID myID = getPartID();
  // Note that the size of vertsReceived won't necessarily match numTotalGhosts,
  // because for some verts, we're getting incoming info about multiple copies.
  // Stage 3:  Send copy info back
  //
  // Send messages back to all parts you just received messages from, with
  // a list of all copy info for every vert they sent info about.
  std::map<EntityGID, Ghost> allGhosts;
  {
    const std::set<Vert*>::iterator iterEnd = vertsReceived.end();
    for (std::set<Vert*>::iterator iter = vertsReceived.begin();
	iter != iterEnd; ++iter) {
      // No need to bother for verts with no copies, or ones we don't own.
      Vert *pV = *iter;
      assert(pV->getEntStatus() != Entity::eInterior);
      // This vert will be owned elsewhere if it's a ghost.

      bool OK = false;
      Ghost g = getGhostEntry(pV, OK);
      assert(OK);
      if (!OK)
	continue; // Something horrible has happened.
      for (GR_index_t ii = 0; ii < g.getNumCopies(); ii++) {
	EntityGID egid = g.getCopyID(ii);
	PartGID target = egid.getPart();
	Ghost foreign(g, target, myID, pV);
	std::pair<std::map<EntityGID, Ghost>::iterator, bool> result =
	    allGhosts.insert(std::make_pair(egid, foreign));
	assert(result.second);
      }
    }
  }

  MPI::Intercomm comm = getPartition()->getCommunicator();
  // Now we have a map of EGID's to all of their copies.  The map is sorted
  // by PartGID's first, so all data that's bound for a given remote part is
  // in order.
  {
    const std::map<EntityGID, Ghost>::iterator iterEnd = allGhosts.end();
    std::map<EntityGID, Ghost>::iterator iter0, iter1;
    iter0 = iter1 = allGhosts.begin();
    while (iter0 != iterEnd) {
      PartGID thisPart = iter0->first.getPart();
      iter1 = iter0;
      ++iter1;
      while (iter1 != iterEnd && thisPart == iter1->first.getPart()) {
	++iter1;
      }
      // Now we have a range from iter0 to iter1.
      // In a first pass over this range, we accumulate the size of the
      // message we're about to send.
      int messageSize = sizeof(int);
      int numEntities = 0;
      for (std::map<EntityGID, Ghost>::iterator iter = iter0; iter != iter1;
	  ++iter, ++numEntities) {
	messageSize += sizeof(Entity*) + iter->second.MPIdataSize();
      }

      // Now create a buffer to put the data into
      char *buffer = new char[messageSize];
      int offset = 0, targetRank = -1;
      transcribeToBuffer(buffer, offset, &numEntities, sizeof(int));
      for (std::map<EntityGID, Ghost>::iterator iter = iter0; iter != iter1;
	  ++iter, ++numEntities) {
	const Entity* ent = iter->first.getEnt();
	transcribeToBuffer(buffer, offset, &ent, sizeof(Entity*));
	iter->second.serialize(buffer, offset);
	targetRank = iter->first.getPart().getRank();
      }
      assert(offset == messageSize);
      logMessage(
	  3, "Sending message of length %d from part %d to rank %d (tag %d)\n",
	  messageSize, myID.getUniqueID(), targetRank, PartBdryVertCopies);
      comm.Isend(buffer, messageSize, MPI::CHAR, targetRank, tag);
      iter0 = iter1;
      delete [] buffer;
    }
  }
}

void
Mesh::matchPartBdryVertices()
{
  Partition* ptn = getPartition();
  MPI::Intercomm comm = ptn->getCommunicator();
  if (!ptn || ptn->getNumGlobalParts() == 1)
    return;

  PartGID myID = getPartID();
  GhostEnts GV;

  // Stage 1:  Send local info to owners
  //
  // Assemble and send messages with the global handles and coords of
  // all foreign verts, one to each neighboring part.

  // 1a:  Create a data structure to hold all this info, sorted by the
  // part to which to send data.
  for (GR_index_t ii = 0; ii < getNumVerts(); ii++) {
    Vert* pV = getVert(ii);
    if (!pV->isOwned()) {
      // Look up the ghost info for this vertex
      bool success = false;
      Ghost ghost = getGhostEntry(pV, success);
      success = GV.addEntity(pV, ghost.getOwnerID());
      assert(success);
    }
  }
  // 1b:  Create and send a message for each.
  // How many will there be?  We care because this tells us how many to
  // expect to receive.
  int numPartsSentTo = GV.size();

  GhostEnts::GE_const_iter GVend = GV.end();
  for (GhostEnts::GE_const_iter iter = GV.begin(); iter != GVend; ++iter) {
    logMessage(4, "Part %d: target PID %d, %lu verts\n", myID.getUniqueID(),
	       iter->first.getUniqueID(), iter->second.size());
    sendVertInfo(PartBdryVertCoords, iter->first, iter->second);
  }

  // Stage 2:  Assemble info about copies
  //
  // Receive messages with foreign EntityGID's for verts we own; store the
  // EGID's.  Keep receiving until there have been enough messages per vert.
  int numTotalGhosts = 0;
  for (GR_index_t ii = 0; ii < getNumVerts(); ii++) {
    Vert* pV = getVert(ii);
    if (!pV->isOwned() || pV->getEntStatus() == Entity::eInterior)
      continue;
    bool success = false;
    Ghost ghost = getGhostEntry(pV, success);
    assert(success);
    logMessage(4, "Part %d: adding %d copies\n", myID.getUniqueID(),
	       ghost.getMaxNumCopies());
    numTotalGhosts += ghost.getMaxNumCopies();
  }

  std::set<Vert*> vertsReceived;
  {
    ADT tree(this);
    int numPendingGhosts = numTotalGhosts;
    int numMessagesReceived = 0;
    logMessage(4, "Part %d needs to receive %d verts\n", myID.getUniqueID(),
	       numPendingGhosts);
    while (numPendingGhosts > 0) {
      int numGhostsReceived = receiveGhostVertMessage(PartBdryVertCoords, this,
						      tree, vertsReceived);
      numMessagesReceived++;
      numPendingGhosts -= numGhostsReceived;
      logMessage(4, "Part %d received %d verts; %d to go.\n",
		 myID.getUniqueID(), numGhostsReceived, numPendingGhosts);
    }
  }
  // Note that the size of vertsReceived won't necessarily match numTotalGhosts,
  // because for some verts, we're getting incoming info about multiple copies.

  // Stage 3:  Send copy info back
  //
  // Send messages back to all parts you just received messages from, with
  // a list of all copy info for every vert they sent info about.
  replyWithVertIDs(PartBdryVertCopies, vertsReceived);

  // Stage 4:  Receive copy info
  //
  // Receive and store info about all copies of foreign verts this part has
  // copies of.  # of messages is the same as Stage 1.
  for (int ii = 0; ii < numPartsSentTo; ++ii) {
    PartGID senderID;
    ptn->receiveCopyData(PartBdryVertCopies, senderID, this);
  }
}

void
Mesh::matchPartBdryFaces()
{
  // Stage 1:  Send local PBF info to others
  //
  // Assemble and send messages with the global handles and coords of
  // all foreign verts, one to each neighboring part.

  GhostEnts GE;

  // 1a:  Create a data structure to hold all this info, sorted by the
  // part to which to send data.
  for (GR_index_t ii = 0; ii < getNumFaces(); ii++) {
    Face* pF = getFace(ii);
    if (pF->getEntStatus() == Entity::ePartBdry) {
      // Look up the ghost info for this face
      bool success = false;
      Ghost ghost = getGhostEntry(pF, success);
      PartGID target = ghost.getOwnerID();
      if (target == getPartID()) {
	// Need to get the remote ID instead.  There can only be one
	// copy.
	assert(ghost.getNumCopies() == 1);
	EntityGID egid = ghost.getCopyID(0);
	assert(!(egid.getPart() == getPartID()));
	target = egid.getPart();
      }
      success = GE.addEntity(pF, target);
      assert(success);
    }
  }
  // 1b:  Create and send a message for each.
  // Need to know how many we're sending, because that's the number we'll
  // need to receive shortly.
  GR_index_t numPartsSentTo = GE.size();

  GhostEnts::GE_const_iter GEend = GE.end();
  for (GhostEnts::GE_const_iter iter = GE.begin(); iter != GEend; ++iter) {
    sendConnectInfo(PartBdryFaces, iter->first, this, iter->second);
  }

  int numConnectedVerts = 0, numPartBdryFaces = 0;
  for (GR_index_t ii = 0; ii < getNumFaces(); ++ii) {
    Face *pF = getFace(ii);
    if (pF->getEntStatus() == Entity::ePartBdry) {
      numConnectedVerts += pF->getNumVerts();
      ++numPartBdryFaces;
    }
  }

  // Now catch the data that others sent to you.  The number of messages
  // received will match the number sent.
  for (GR_index_t ii = 0; ii < numPartsSentTo; ++ii) {
    receivePBFaceInfo(this);
  }
}

void
Mesh::matchPartBdrys()
{
  matchPartBdryVertices();
  matchPartBdryFaces();
}

template<typename C>
  static void
  addGhostedNeighbors(const std::set<C*>& ents, const Ghost& g,
		      GhostEnts& ghosts)
  {
    GR_index_t nCopies = g.getNumCopies();
    for (typename std::set<C*>::iterator iter = ents.begin();
	iter != ents.end(); ++iter) {
      C* tmp = *iter;
      for (GR_index_t ii = 0; ii < nCopies; ++ii) {
	PartGID pgid = g.getCopyID(ii).getPart();
	ghosts.addEntity(tmp, pgid);
      }
      if (tmp->getEntStatus() != Entity::ePartBdry) {
	tmp->setEntStatus(Entity::eGhosted);
      }
    }
  }

void
Mesh::clearGhostData()
{
  // For verts, clear ghosting info and info about foreign entities.
  for (GR_index_t ii = 0; ii < getNumVerts(); ++ii) {
    Vert *pV = getVert(ii);
    if (pV->isDeleted())
      continue;
    Entity::Entity_Status_t stat = pV->getEntStatus();
    switch (stat)
      {
      case Entity::eForeign:
      case Entity::eGhosted:
	removeGhostEntry(pV);
	break;
      case Entity::ePartBdry:
	resetGhostEntry(pV);
	break;
      default:
	break;
      }
  }

  // For faces, clear ghosting info and info about foreign entities.
  for (GR_index_t ii = 0; ii < getNumFaces(); ++ii) {
    Face *pF = getFace(ii);
    if (pF->isDeleted())
      continue;
    Entity::Entity_Status_t stat = pF->getEntStatus();
    if ((stat == Entity::eForeign) || (stat == Entity::eGhosted)) {
      removeGhostEntry(pF);
    }
  }

  // For cells, clear ghosting info and info about foreign entities.  Also,
  // delete foreign cells.  This will lead to auto-deletion of foreign faces.
  for (GR_index_t ii = 0; ii < getNumCells(); ++ii) {
    Cell *pC = getCell(ii);
    if (pC->isDeleted())
      continue;
    Entity::Entity_Status_t stat = pC->getEntStatus();
    if (stat == Entity::eForeign) {
      removeGhostEntry(pC);
      deleteCell(pC);
    }
    else if (stat == Entity::eGhosted) {
      removeGhostEntry(pC);
    }
  }

  // Now remove the foreign verts, which don't auto-remove
  for (GR_index_t ii = 0; ii < getNumVerts(); ++ii) {
    Vert *pV = getVert(ii);
    if (pV->isDeleted())
      continue;
    Entity::Entity_Status_t stat = pV->getEntStatus();
    if (stat == Entity::eForeign) {
      bool result = deleteVert(pV);
      assert(result);
    }
  }

  purgeCells();
  purgeFaces();
  purgeVerts();
}

void
Mesh::setupGhostData()
{
  // Step 1: Collect neighbor data, by parts.  Specifically, we need vertices and
  // cells.
  GhostEnts verts, cells;
  for (GR_index_t ii = 0; ii < getNumVerts(); ++ii) {
    Vert *pV = getVert(ii);
    if (pV->getEntStatus() == Entity::ePartBdry) {
      std::set<Cell*> spCInc;
      std::set<Vert*> spVNeigh;
      findNeighborhoodInfo(pV, spCInc, spVNeigh);

      bool status = false;
      Ghost g(getGhostEntry(pV, status));
      assert(status);
      logMessage(
	  4,
	  "Part %u checking vert at (%f %f %f) %lu verts, %lu cells; %d remote copies\n",
	  getPartID().getUniqueID(), pV->x(), pV->y(), pV->z(), spVNeigh.size(),
	  spCInc.size(), g.getNumCopies());

      // Every part that has a copy of this vert needs to get a copy of
      // these neighbor verts and cells.  Also, the verts and cells being
      // sent are marked with status "ghosted" unless they're already
      // part bdry entities.
      addGhostedNeighbors(spVNeigh, g, verts);
      addGhostedNeighbors(spCInc, g, cells);
    }
  }
  assert(verts.size() == cells.size());
  GR_index_t numPartsSentTo = verts.size();

  // Stage 2:  Sending info about verts that will be ghosted to all neighboring parts.
  GhostEnts::GE_const_iter GVend = verts.end();
  for (GhostEnts::GE_const_iter iter = verts.begin(); iter != GVend; ++iter) {
    sendVertInfo(GhostVertCoords, iter->first, iter->second);
  }

  ADT tree(this);
  assert(tree.qValid());
  assert(tree.iTreeSize() > 0);
  std::set<Vert*> vertsReceived;

  for (GR_index_t ii = 0; ii < numPartsSentTo; ++ii) {
    receiveGhostVertMessage(GhostVertCoords, this, tree, vertsReceived);
  }

  // Stage 3:  Send vertex copy info back
  //
  // Send messages back to all parts you just received messages from, with
  // a list of all copy info for every vert they sent info about.
  replyWithVertIDs(GhostVertCopies, vertsReceived);

  // Stage 4:  Receive vertex copy info
  //
  // Receive and store info about all copies of foreign verts this part has
  // copies of.  # of messages is the same as Stage 1.
  for (GR_index_t ii = 0; ii < numPartsSentTo; ++ii) {
    PartGID senderID;
    getPartition()->receiveCopyData(GhostVertCopies, senderID, this);
  }

  // Stage 5: Send cell ghost info
  GhostEnts::GE_const_iter GEend = cells.end();
  for (GhostEnts::GE_const_iter iter = cells.begin(); iter != GEend; ++iter) {
    sendConnectInfo(GhostEntities, iter->first, this, iter->second);
  }

  // Stage 6: Receive cell ghost info
  // Now catch the data that others sent to you.  The number of messages
  // received will match the number sent.
  for (GR_index_t ii = 0; ii < numPartsSentTo; ++ii) {
    receiveCellInfo(this);
  }
}
// Update all ghost information based on the current ghosting rules.
// TODO: Initially, the only ghosting rule (hard-coded) is to share everything
// incident on verts on the part bdry, plus the closures of those entities.

void
GRUMMP::transcribeToBuffer(char* const buffer, int& offset,
			   const void* const data, const size_t blockSize)
{
  memcpy(buffer + offset, data, blockSize);
  offset += blockSize;
}

void
GRUMMP::transcribeFromBuffer(const char* const buffer, int& offset,
			     void* const data, const size_t blockSize)
{
  memcpy(data, buffer + offset, blockSize);
  offset += blockSize;
}

void
Partition::sendVertexData(const GhostMsgTags_t messageTag,
			  const PartGID& target, const int numVerts,
			  const double coords[], Vert* verts[])
{
  // First calculate the total message size
  size_t bufferSize = sizeof(int) // numVerts
  + sizeof(double) * numVerts * 3 // coords
  + sizeof(Entity**) * numVerts; // handles
  char* buffer = new char[bufferSize];

  logMessage(3, "Sending to %d: %d verts\n", target.getRank(), numVerts);
  // Now transcribe data, in the order the arguments are given.  memcpy is
  // low-level, but faster than looping over anything by hand.
  int currLoc = 0;

  // numVerts
  transcribeToBuffer(buffer, currLoc, &numVerts, sizeof(int));

  // coords
  transcribeToBuffer(buffer, currLoc, coords, sizeof(double) * numVerts * 3);

  // handles
  transcribeToBuffer(buffer, currLoc, verts, sizeof(Vert*) * (numVerts));

  assert(size_t(currLoc) == bufferSize);
  m_comm.Isend(buffer, bufferSize, MPI::CHAR, target.getRank(), messageTag);

  delete[] buffer;
}

void
Partition::sendConnectivityData(const GhostMsgTags_t messageTag,
				const PartGID& source, const PartGID& target,
				const int numEntities, const Entity* handles[],
				const int offsets[], const int numVerts,
				const Vert* verts[])
{
  // First calculate the total message size
  size_t bufferSize = sizeof(PartGID) + sizeof(int) // numEntities
      + sizeof(Entity*) * numEntities // handles
      + sizeof(int) * (numEntities + 1) // offsets
      + sizeof(int) // numVerts
      + sizeof(Vert*) * numVerts; // indices
  char* buffer = new char[bufferSize];

  logMessage(3, "Sending to %d: %d entities, %d connectivity entries\n",
	     target.getRank(), numEntities, numVerts);
  // Now transcribe data, in the order the arguments are given.  memcpy is
  // low-level, but faster than looping over anything by hand.
  int currLoc = 0;

  // Source part GID
  transcribeToBuffer(buffer, currLoc, &source, sizeof(PartGID));

  // numEntities
  transcribeToBuffer(buffer, currLoc, &numEntities, sizeof(int));

  // handles
  transcribeToBuffer(buffer, currLoc, handles, sizeof(Entity*) * numEntities);

  // offsets
  transcribeToBuffer(buffer, currLoc, offsets, sizeof(int) * (numEntities + 1));

  // numVerts
  transcribeToBuffer(buffer, currLoc, &numVerts, sizeof(int));

  // handles
  transcribeToBuffer(buffer, currLoc, verts, sizeof(Vert*) * numVerts);

  assert(size_t(currLoc) == bufferSize);
  m_comm.Isend(buffer, bufferSize, MPI::CHAR, target.getRank(), messageTag);

  delete[] buffer;
}

void
Partition::receiveCopyData(const GhostMsgTags_t messageTag, PartGID& senderID,
			   Mesh* pMesh)
{
  bool isBdry = (messageTag == PartBdryVertCopies);
  MPI::Status stat;
  m_comm.Probe(MPI::ANY_SOURCE, messageTag, stat);
  size_t messageSize = stat.Get_count(MPI::CHAR);
  // Get a message and unpack it; make sure it's the one you just got info
  // about with Probe.
  char *buffer = new char[messageSize];
  m_comm.Recv(buffer, messageSize, MPI::CHAR, stat.Get_source(), stat.Get_tag(),
	      stat);

  senderID = stat.Get_source();

  int numGhostedEntities = -1, offset = 0, numTotalCopies = 0;
  transcribeFromBuffer(buffer, offset, &numGhostedEntities, sizeof(int));

  // Now loop over the entities and unpack their data.
  for (int ii = 0; ii < numGhostedEntities; ++ii) {
    Entity* ent;
    transcribeFromBuffer(buffer, offset, &ent, sizeof(Entity*));

    Ghost g(pMesh->m_ghostMap[ent]);
    Ghost g2;
    g2.deserialize(buffer, offset);
    g.merge(g2, isBdry);
    pMesh->m_ghostMap[ent] = g;
    numTotalCopies += g.getNumCopies();
  }
  logMessage(3, "Part %d received %d copies of %d entities\n",
	     pMesh->getPartID().getUniqueID(), numTotalCopies,
	     numGhostedEntities);

  delete [] buffer;
}

void
Partition::receiveVertexData(const GhostMsgTags_t messageTag, PartGID& senderID,
			     int& numVerts, double*& coords, Vert**& handles)
{
  MPI::Status stat;
  m_comm.Probe(MPI::ANY_SOURCE, messageTag, stat);
  size_t messageSize = stat.Get_count(MPI::CHAR);
  // Get a message and unpack it; make sure it's the one you just got info
  // about with Probe.
  char *buffer = new char[messageSize];
  m_comm.Recv(buffer, messageSize, MPI::CHAR, stat.Get_source(), stat.Get_tag(),
	      stat);

  senderID = stat.Get_source();

  // Assemble local ghost data from that.
  // Send a reply with info about local handles for those entities.

  // Now read out data, in the order the arguments are given.  memcpy is
  // low-level, but faster than looping over anything by hand.
  int currLoc = 0;

  // numVerts
  transcribeFromBuffer(buffer, currLoc, &numVerts, sizeof(int));

  // coords
  coords = new double[numVerts * 3];
  transcribeFromBuffer(buffer, currLoc, coords, sizeof(double) * numVerts * 3);

  // handles
  handles = new Vert*[numVerts];
  transcribeFromBuffer(buffer, currLoc, handles, sizeof(Vert*) * (numVerts));

  assert(size_t(currLoc) == messageSize);
  size_t sizeCheck = sizeof(int) // numVerts
  + sizeof(double) * numVerts * 3 // coords
  + sizeof(Vert*) * numVerts; // handles
  assert(size_t(sizeCheck) == messageSize);
  logMessage(3, "Received from %d: %d verts\n", senderID.getRank(), numVerts);
  delete[] buffer;
}

void
Partition::receiveConnectivityData(const GhostMsgTags_t messageTag,
				   PartGID& source, int& numEntities,
				   Entity**& handles, int*& offsets,
				   int& numVerts, Vert**& verts)
{
  MPI::Status stat;
  m_comm.Probe(MPI::ANY_SOURCE, messageTag, stat);
  size_t messageSize = stat.Get_count(MPI::CHAR);
  // Get a message and unpack it; make sure it's the one you just got info
  // about with Probe.
  char *buffer = new char[messageSize];
  m_comm.Recv(buffer, messageSize, MPI::CHAR, stat.Get_source(), stat.Get_tag(),
	      stat);

  // Now read out data, in the order the arguments are given.  memcpy is
  // low-level, but faster than looping over anything by hand.
  int currLoc = 0;

  // Source part GID
  transcribeFromBuffer(buffer, currLoc, &source, sizeof(PartGID));

  // numEntities
  transcribeFromBuffer(buffer, currLoc, &numEntities, sizeof(int));

  // handles
  handles = new Entity*[numEntities];
  transcribeFromBuffer(buffer, currLoc, handles,
		       sizeof(Entity*) * (numEntities));

  // offsets
  offsets = new int[numEntities + 1];
  transcribeFromBuffer(buffer, currLoc, offsets,
		       sizeof(int) * (numEntities + 1));

  // numVerts
  transcribeFromBuffer(buffer, currLoc, &numVerts, sizeof(int));

  // indices
  verts = new Vert*[numVerts];
  transcribeFromBuffer(buffer, currLoc, verts, sizeof(Vert*) * numVerts);

  assert(size_t(currLoc) == messageSize);
  logMessage(3, "Received from %d: %d entities, %d connectivity entries\n",
	     source.getRank(), numEntities, numVerts);
  delete[] buffer;
}

#endif
