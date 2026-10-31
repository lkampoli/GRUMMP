#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"
#include "iMesh_GRUMMP_misc.hh"

#include <vector>
#include <algorithm>
#include <iostream>

#include "iBase.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"

namespace ITAPS_GRUMMP
{

//////////////////////////////////////////////////////////////////////
////  EntityList (non-unique, ordered)                            ////
//////////////////////////////////////////////////////////////////////

  EntityList::EntityList(const int iDim) :
      EntitySetBase(iDim), lpEntities(), sVRAux(), spIters(), spWorksets(), iVerts(
	  0), iEdges(0), iTriFaces(0), iQuadFaces(0), iTriCells(0), iQuadCells(
	  0), iTetCells(0), iPyrCells(0), iPrismCells(0), iHexCells(0)
  {
  }

  EntityList::EntityList(const EntityList& EL) :
      EntitySetBase(EL), lpEntities(EL.lpEntities), sVRAux(EL.sVRAux), spIters(
	  EL.spIters), spWorksets(EL.spWorksets), iVerts(EL.iVerts), iEdges(
	  EL.iEdges), iTriFaces(EL.iTriFaces), iQuadFaces(EL.iQuadFaces), iTriCells(
	  EL.iTriCells), iQuadCells(EL.iQuadCells), iTetCells(EL.iTetCells), iPyrCells(
	  EL.iPyrCells), iPrismCells(EL.iPrismCells), iHexCells(EL.iHexCells)
  {
  }

  void
  EntityList::vGetEntities(const int entity_type, const int entity_topology,
			   iBase_EntityHandle** entity_handles,
			   int *entity_handles_allocated,
			   int *entity_handles_size) const
  {
    // Call the overloaded method to get all the entities of that type/topology
    // into a vector
    std::vector<Entity*> vpE;
    vGetEntities(entity_type, entity_topology, vpE);

    // Then stick them into the array
    *entity_handles_size = vpE.size();
    vEnsureArraySize<iBase_EntityHandle>(entity_handles,
					 entity_handles_allocated,
					 *entity_handles_size);
    std::vector<Entity*>::iterator iter = vpE.begin(), iterEnd = vpE.end();
    for (int iEH = 0; iter != iterEnd; ++iter, iEH++) {
      (*entity_handles)[iEH] = *iter;
    }
  }

  void
  EntityList::vGetEntities(const int entity_type, const int entity_topology,
			   std::vector<Entity*>& vpE) const
  {
    // Don't clear vpE here, because then chaining these together for
    // multiple entity sets won't work.
    std::list<Entity*>::const_iterator iter = lpEntities.begin(), iterEnd =
	lpEntities.end();

    if (entity_topology == iMesh_ALL_TOPOLOGIES
	&& entity_type == iBase_ALL_TYPES) {
      // Insert all the entities into the vector
      vpE.insert(vpE.end(), iter, iterEnd);
    }
    else if (entity_type == iBase_ALL_TYPES) {
      // Insert those entities with the desired topology
      for (; iter != iterEnd; ++iter) {
	if ((*iter)->getEntTopology() == entity_topology)
	  vpE.push_back(*iter);
      }
    }
    else if (entity_topology == iMesh_ALL_TOPOLOGIES) {
      // Insert those entities with the desired type
      for (; iter != iterEnd; ++iter) {
	if ((*iter)->getEntType() == entity_type)
	  vpE.push_back(*iter);
      }
    }
    else {
      // Insert those entities with the desired type and topology
      for (; iter != iterEnd; ++iter) {
	if ((*iter)->getEntType() == entity_type
	    && (*iter)->getEntTopology() == entity_topology)
	  vpE.push_back(*iter);
      }
    }
  }

  void
  EntityList::vAddVertex(Vert* const pV)
  {
    // Unlike the entity set class, this method is only called for auxilliary
    // vertices (which are not contained explicitly in the list)
    VertRec VR(pV);
    std::set<VertRec>::iterator iter = sVRAux.find(VR);
    // To keep track of how many copies of each vertex there are, must see if
    // the vertex was already added
    if (iter == sVRAux.end()) {
      sVRAux.insert(VR);
    }
    else {
      // Cannot modify the existing VertRec, so one must delete the old one
      // and insert a new one with an incremented count
      VR.vSetCount(iter->iCount());
      VR.iIncrementCount();
      sVRAux.erase(iter);
      sVRAux.insert(VR);
    }
  }

  void
  EntityList::vRemoveVertex(Vert* const pV)
  {
    // Only called for auxilliary vertices
    qVertNumsOK = false;
    VertRec VR(pV);
    std::set<VertRec>::iterator iter = sVRAux.find(VR);
    if (iter == sVRAux.end()) {
      // This should only happen if a user tries to remove an entity
      // that was never added to the ES.
      // GRTS_throw(INVALID_ARGUMENT);
      return;
    }
    else {
      // Cannot modify the existing VertRec, so one must delete the old one
      // and insert a new one with a decremented count
      int iCount = iter->iCount();
      iCount--;
      assert(iCount >= 0);
      sVRAux.erase(iter);
      if (iCount != 0) {
	VR.vSetCount(iCount);
	sVRAux.insert(VR);
      }
    }
  }

  void
  EntityList::vAddEntities(const iBase_EntityHandle entity_handles[],
			   const int entity_handles_size)

  {
    for (int iEH = 0; iEH < entity_handles_size; iEH++) {
      vAddEntity(entity_handles[iEH]);
    }
  }

  void
  EntityList::vAddEntity(void* const entity_handle)

  {
    Entity *pE = static_cast<Entity*>(entity_handle);

    // Ensure that the verts for each entity will be included, at
    // least in shadow form.
    if (pE->getEntType() != iBase_VERTEX) {
      int iStop = pE->getNumVerts();
      for (int i = 0; i < iStop; i++)
	vAddVertex(pE->getVert(i));
    }

    switch (pE->getEntTopology())
      {
      case iMesh_POINT:
	assert(dynamic_cast<Vert*>(pE));
	iVerts++;
	if (!sVRAux.empty())
	  qVertNumsOK = false;
	break;
      case iMesh_LINE_SEGMENT:
	assert(dynamic_cast<EdgeFace*>(pE));
	iEdges++;
	break;
      case iMesh_TRIANGLE:
	if (iMeshTopoDim == 3) {
	  assert(dynamic_cast<TriFace*>(pE));
	  iTriFaces++;
	}
	else {
	  assert(dynamic_cast<TriCell*>(pE));
	  iTriCells++;
	}
	break;
      case iMesh_QUADRILATERAL:
	if (iMeshTopoDim == 3) {
	  assert(dynamic_cast<QuadFace*>(pE));
	  iQuadFaces++;
	}
	else {
	  assert(dynamic_cast<QuadCell*>(pE));
	  iQuadCells++;
	}
	break;
      case iMesh_TETRAHEDRON:
	assert(dynamic_cast<TetCell*>(pE));
	iTetCells++;
	break;
      case iMesh_PYRAMID:
	assert(dynamic_cast<PyrCell*>(pE));
	iPyrCells++;
	break;
      case iMesh_PRISM:
	assert(dynamic_cast<PrismCell*>(pE));
	iPrismCells++;
	break;
      case iMesh_HEXAHEDRON:
	assert(dynamic_cast<HexCell*>(pE));
	iHexCells++;
	break;
      case iMesh_POLYGON:
      case iMesh_POLYHEDRON:
      case iMesh_SEPTAHEDRON:
	throw(iBase_NOT_SUPPORTED);
	break;
      default:
	// How did somebody create an entity of this topology?
	throw(iBase_INVALID_ARGUMENT);
      }

    lpEntities.push_back(pE);

    // Inform all iterators for this list of the addition.
    std::set<EntityListIter*>::iterator itersIter = spIters.begin(),
	itersIterEnd = spIters.end();
    for (; itersIter != itersIterEnd; itersIter++)
      (*itersIter)->vAdded(pE);

    // Inform all worksets for this list of the addition.
    std::set<EntityListWorkset*>::iterator setsIter = spWorksets.begin(),
	setsIterEnd = spWorksets.end();
    for (; setsIter != setsIterEnd; setsIter++)
      (*setsIter)->vAdded(pE);
  }

  void
  EntityList::vRemoveEntity(Entity* const pE)
  {
    std::list<Entity*>::iterator iter, iterFound, iterEnd = lpEntities.end();
    iterFound = iterEnd;
    // Need to find the last copy of the entity in the list (the spec says
    // "the most recently added copy is removed") so cannot use std::find
    for (iter = lpEntities.begin(); iter != iterEnd; ++iter)
      if ((*iter) == pE)
	iterFound = iter;
    // If no copies of the entity exist, better quit
    if (iterFound == iterEnd)
      return;

    // Inform all iterators for this list of the removal.
    std::set<EntityListIter*>::iterator itersIter = spIters.begin(),
	itersIterEnd = spIters.end();
    for (; itersIter != itersIterEnd; itersIter++)
      (*itersIter)->vDelete(iterFound);

    // Inform all worksets for this list of the removal.
    std::set<EntityListWorkset*>::iterator setsIter = spWorksets.begin(),
	setsIterEnd = spWorksets.end();
    for (; setsIter != setsIterEnd; setsIter++)
      (*setsIter)->vDelete(iterFound);

    // Ensure that the verts for each entity will be removed, even
    // though they may not be members of the entity set.
    if (pE->getEntTopology() != iMesh_POINT) {
      int iStop = pE->getNumVerts();
      for (int i = 0; i < iStop; i++)
	vRemoveVertex(pE->getVert(i));
    }

    switch (pE->getEntTopology())
      {
      case iMesh_POINT:
	assert(dynamic_cast<Vert*>(pE));
	iVerts--;
	if (!sVRAux.empty())
	  qVertNumsOK = false;
	break;
      case iMesh_LINE_SEGMENT:
	assert(dynamic_cast<EdgeFace*>(pE));
	iEdges--;
	break;
      case iMesh_TRIANGLE:
	if (dynamic_cast<TriFace*>(pE) != NULL)
	  iTriFaces--;
	else
	  iTriCells--;
	break;
      case iMesh_QUADRILATERAL:
	if (dynamic_cast<QuadFace*>(pE) != NULL)
	  iQuadFaces--;
	else
	  iQuadCells--;
	break;
      case iMesh_TETRAHEDRON:
	assert(dynamic_cast<TetCell*>(pE));
	iTetCells--;
	break;
      case iMesh_PYRAMID:
	assert(dynamic_cast<PyrCell*>(pE));
	iPyrCells--;
	break;
      case iMesh_PRISM:
	assert(dynamic_cast<PrismCell*>(pE));
	iPrismCells--;
	break;
      case iMesh_HEXAHEDRON:
	assert(dynamic_cast<HexCell*>(pE));
	iHexCells--;
	break;
      case iMesh_POLYGON:
      case iMesh_POLYHEDRON:
      case iMesh_SEPTAHEDRON:
	throw(iBase_NOT_SUPPORTED);
	break;
      default:
	// How did somebody create an entity of this topology?
	throw(iBase_INVALID_ARGUMENT);
      }

    lpEntities.erase(iterFound);
  }

  void
  EntityList::vRemoveEntities(const iBase_EntityHandle entity_handles[],
			      const int iNumEnt)

  {
    for (int iEH = 0; iEH <= iNumEnt; iEH++) {
      Entity *pE = static_cast<Entity*>(entity_handles[iEH]);
      vRemoveEntity(pE);
    }
  }

  bool
  EntityList::qIsEntityInSet(void* entity_handle) const

  {
    Entity *pE = static_cast<Entity*>(entity_handle);
    std::list<Entity*>::const_iterator iter = std::find(lpEntities.begin(),
							lpEntities.end(), pE);
    return (iter != lpEntities.end());
  }

  void
  EntityList::vRebuildAuxiliaryVertList()
  {
    std::list<Entity*>::iterator iter = lpEntities.begin(), iterEnd =
	lpEntities.end();
    int iVEnd;
    Entity* pETemp;
    while (iter != iterEnd) {
      pETemp = *iter;
      if (pETemp->getEntType() != iBase_VERTEX) {
	iVEnd = pETemp->getNumVerts();
	for (int iV = 0; iV < iVEnd; iV++)
	  vAddVertex(pETemp->getVert(iV));
      }
      ++iter;
    }
    qVertNumsOK = false;
  }
}
