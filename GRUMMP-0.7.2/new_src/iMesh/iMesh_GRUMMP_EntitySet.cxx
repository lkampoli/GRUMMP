#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"
#include "iMesh_GRUMMP_misc.hh"

#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"

template<class T>
  static bool inline
  qIsInSet(const std::set<T>& s, const T& TT)
  {
    typename std::set<T>::const_iterator iter = s.find(TT);
    return (!(iter == s.end()));
  }

namespace ITAPS_GRUMMP
{

//////////////////////////////////////////////////////////////////////
////  EntitySet  (unique, non-ordered)                            ////
//////////////////////////////////////////////////////////////////////

  EntitySet::EntitySet(const int iDim) :
      EntitySetBase(iDim), sVR(), sVRAux(), spEdgeF(), spTriF(), spQuadF(), spTriC(), spQuadC(), spTetC(), spPyrC(), spPrismC(), spHexC(), spIters(), spWorksets()
  {
  }

  EntitySet::EntitySet(const EntityList& EL) :
      EntitySetBase(EL.iMeshTopologicalDimension()), sVR(), sVRAux(), spEdgeF(), spTriF(), spQuadF(), spTriC(), spQuadC(), spTetC(), spPyrC(), spPrismC(), spHexC(), spIters(), spWorksets()
  {
    Entity* pE;
    std::list<Entity*>::const_iterator iter = EL.lpEntities.begin(), iterEnd =
	EL.lpEntities.end();
    for (; iter != iterEnd; ++iter) {
      pE = const_cast<Entity*>(*iter);
      switch (pE->getEntTopology())
	{
	case iMesh_POINT:
	  {
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV);
	    sVR.insert(VertRec(pV));
	    break;
	  }
	case iMesh_LINE_SEGMENT:
	  assert(dynamic_cast<EdgeFace*>(pE));
	  spEdgeF.insert(dynamic_cast<EdgeFace*>(pE));
	  break;
	case iMesh_TRIANGLE:
	  if (dynamic_cast<TriFace*>(pE) != NULL)
	    spTriF.insert(dynamic_cast<TriFace*>(pE));
	  else {
	    assert(dynamic_cast<TriCell*>(pE));
	    spTriC.insert(dynamic_cast<TriCell*>(pE));
	  }
	  break;
	case iMesh_QUADRILATERAL:
	  if (dynamic_cast<QuadFace*>(pE) != NULL)
	    spQuadF.insert(dynamic_cast<QuadFace*>(pE));
	  else {
	    assert(dynamic_cast<QuadCell*>(pE));
	    spQuadC.insert(dynamic_cast<QuadCell*>(pE));
	  }
	  break;
	case iMesh_TETRAHEDRON:
	  assert(dynamic_cast<TetCell*>(pE));
	  spTetC.insert(dynamic_cast<TetCell*>(pE));
	  break;
	case iMesh_PYRAMID:
	  assert(dynamic_cast<PyrCell*>(pE));
	  spPyrC.insert(dynamic_cast<PyrCell*>(pE));
	  break;
	case iMesh_PRISM:
	  assert(dynamic_cast<PrismCell*>(pE));
	  spPrismC.insert(dynamic_cast<PrismCell*>(pE));
	  break;
	case iMesh_HEXAHEDRON:
	  assert(dynamic_cast<HexCell*>(pE));
	  spHexC.insert(dynamic_cast<HexCell*>(pE));
	  break;
	case iMesh_POLYGON:
	case iMesh_POLYHEDRON:
	case iMesh_SEPTAHEDRON:
	  throw(iBase_NOT_SUPPORTED);
	default:
	  // How did somebody create an entity of this topology?
	  throw(iBase_INVALID_ARGUMENT);
	}
    }

    vRebuildAuxiliaryVertList();

    std::set<EntitySetBase*>::iterator setIter = EL.spESMembers.begin(),
	setIterEnd = EL.spESMembers.end();
    for (; setIter != setIterEnd; setIter++)
      vAddEntitySet(*setIter);
  }

  EntitySet::EntitySet(const EntitySet& ES) :
      EntitySetBase(ES), sVR(ES.sVR), sVRAux(ES.sVRAux), spEdgeF(ES.spEdgeF), spTriF(
	  ES.spTriF), spQuadF(ES.spQuadF), spTriC(ES.spTriC), spQuadC(
	  ES.spQuadC), spTetC(ES.spTetC), spPyrC(ES.spPyrC), spPrismC(
	  ES.spPrismC), spHexC(ES.spHexC), spIters(ES.spIters), spWorksets(
	  ES.spWorksets)
  {
    vRebuildAuxiliaryVertList();
  }

  void
  EntitySet::vGetEntities(const int entity_type, const int entity_topology,
			  iBase_EntityHandle** entity_handles,
			  int *entity_handles_allocated,
			  int *entity_handles_size) const

  {
    std::set<Entity*> spE;
    vGetEntities(entity_type, entity_topology, spE);
    *entity_handles_size = spE.size();

    vEnsureArraySize<iBase_EntityHandle>(entity_handles,
					 entity_handles_allocated,
					 *entity_handles_size);
    std::set<Entity*>::iterator iter = spE.begin();
    int iEH;
    for (iEH = 0; iter != spE.end(); ++iter, iEH++) {
      (*entity_handles)[iEH] = *iter;
    }
  }

  void
  EntitySet::vGetEntities(const int entity_type, const int entity_topology,
			  std::set<Entity*>& spE) const

  {
    // Don't clear spE here, because then chaining these together for
    // multiple entity sets won't work.
    if (entity_topology == iMesh_POINT
	|| (entity_topology == iMesh_ALL_TOPOLOGIES
	    && entity_type == iBase_VERTEX)) {
      spE.insert(sVR.begin(), sVR.end());
    }
    else if (entity_topology == iMesh_LINE_SEGMENT
	|| (entity_topology == iMesh_ALL_TOPOLOGIES && entity_type == iBase_EDGE)) {
      spE.insert(spEdgeF.begin(), spEdgeF.end());
    }
    else if (entity_topology == iMesh_TRIANGLE) {
      // IMO (CFO-G), a well-behaved entity set will never have both
      // tri cells and tri faces, but you never know...
      if (!spTriC.empty()) {
	spE.insert(spTriC.begin(), spTriC.end());
      }
      if (!spTriF.empty()) {
	spE.insert(spTriF.begin(), spTriF.end());
      }
    }
    else if (entity_topology == iMesh_QUADRILATERAL) {
      // IMO (CFO-G), a well-behaved entity set will never have both
      // quad cells and quad faces, but you never know...
      if (!spQuadC.empty()) {
	spE.insert(spQuadC.begin(), spQuadC.end());
      }
      if (!spQuadF.empty()) {
	spE.insert(spQuadF.begin(), spQuadF.end());
      }
    }
    else if (entity_topology == iMesh_ALL_TOPOLOGIES
	&& entity_type == iBase_FACE) {
      ////////////////////////////////////////////////
      // Returning all 2D entities.
      ////////////////////////////////////////////////
      // IMO (CFO-G), a well-behaved entity set will never have both
      // tri cells and tri faces, but you never know...
      if (!spTriC.empty()) {
	spE.insert(spTriC.begin(), spTriC.end());
      }
      if (!spTriF.empty()) {
	spE.insert(spTriF.begin(), spTriF.end());
      }
      // IMO (CFO-G), a well-behaved entity set will never have both
      // quad cells and quad faces, but you never know...
      if (!spQuadC.empty()) {
	spE.insert(spQuadC.begin(), spQuadC.end());
      }
      if (!spQuadF.empty()) {
	spE.insert(spQuadF.begin(), spQuadF.end());
      }
    }
    else if (entity_topology == iMesh_TETRAHEDRON) {
      spE.insert(spTetC.begin(), spTetC.end());
    }
    else if (entity_topology == iMesh_PYRAMID) {
      spE.insert(spPyrC.begin(), spPyrC.end());
    }
    else if (entity_topology == iMesh_PRISM) {
      spE.insert(spPrismC.begin(), spPrismC.end());
    }
    else if (entity_topology == iMesh_HEXAHEDRON) {
      spE.insert(spHexC.begin(), spHexC.end());
    }
    else if (entity_topology == iMesh_ALL_TOPOLOGIES
	&& entity_type == iBase_REGION) {
      spE.insert(spTetC.begin(), spTetC.end());
      spE.insert(spPyrC.begin(), spPyrC.end());
      spE.insert(spPrismC.begin(), spPrismC.end());
      spE.insert(spHexC.begin(), spHexC.end());
    }
    else if (entity_topology == iMesh_ALL_TOPOLOGIES
	&& entity_type == iBase_ALL_TYPES) {
      spE.insert(sVR.begin(), sVR.end());
      spE.insert(spEdgeF.begin(), spEdgeF.end());
      spE.insert(spTriF.begin(), spTriF.end());
      spE.insert(spQuadF.begin(), spQuadF.end());
      spE.insert(spTriC.begin(), spTriC.end());
      spE.insert(spQuadC.begin(), spQuadC.end());
      spE.insert(spTetC.begin(), spTetC.end());
      spE.insert(spPyrC.begin(), spPyrC.end());
      spE.insert(spPrismC.begin(), spPrismC.end());
      spE.insert(spHexC.begin(), spHexC.end());
    }
  }

  void
  EntitySet::vAddVertex(Vert* const pV, const bool qFromEntity)

  {
    if (qFromEntity) {
      VertRec VR(pV);
      std::set<VertRec>::iterator iter = sVRAux.find(VR);
      if (iter == sVRAux.end()) {
	sVRAux.insert(VR);
      }
      else {
	VR.vSetCount(iter->iCount());
	VR.iIncrementCount();
	sVRAux.erase(iter);
	sVRAux.insert(VR);
      }
    }
    else {
      VertRec VR(pV);
      // Always add it; the set can decide whether to include it.
      std::pair<std::set<VertRec>::iterator, bool> p = sVR.insert(VR);
      if (!sVRAux.empty() && p.second)
	qVertNumsOK = false;
    }
  }

  void
  EntitySet::vRemoveVertex(Vert* const pV, const bool qFromEntity)

  {
    qVertNumsOK = false;
    if (qFromEntity) {
      VertRec VR(pV);
      std::set<VertRec>::iterator iter = sVRAux.find(VR);
      if (iter == sVRAux.end()) {
	// This should only happen if a user tries to remove an entity
	// that was never added to the ES.

	throw(iBase_INVALID_ARGUMENT);
      }
      else {
	// Remove one copy of the vertex from the auxiliary set.
	int iCount = iter->iCount();
	iCount--;
	assert(iCount >= 0);
	sVRAux.erase(iter);
	if (iCount != 0) {
	  VR.vSetCount(iCount);
	  sVRAux.insert(VR);
	}
	qVertNumsOK = false;
      }
    }
    else {
      VertRec VR(pV);
      sVR.erase(VR);
      qVertNumsOK = false;
    }
  }

  void
  EntitySet::vAddEntities(const iBase_EntityHandle entity_handles[],
			  const int iNumEnt)

  {
    for (int iEH = 0; iEH < iNumEnt; iEH++) {
      vAddEntity(entity_handles[iEH]);
    }
  }

  void
  EntitySet::vAddEntity(void* const entity_handle)

  {
    Entity *pE = static_cast<Entity*>(entity_handle);
    switch (pE->getEntTopology())
      {
      case iMesh_POINT:
	assert(dynamic_cast<Vert*>(pE));
	vAddVertex(dynamic_cast<Vert*>(pE), false);
	return;
      case iMesh_LINE_SEGMENT:
	{
	  EdgeFace *pEF = dynamic_cast<EdgeFace*>(pE);
	  assert(pEF);
	  if (spEdgeF.find(pEF) == spEdgeF.end())
	    spEdgeF.insert(pEF);
	  else
	    return;
	  break;
	}
      case iMesh_TRIANGLE:
	{
	  TriFace *pTF = dynamic_cast<TriFace*>(pE);
	  if (pTF) {
	    if (spTriF.find(pTF) == spTriF.end())
	      spTriF.insert(pTF);
	    else
	      return;
	  }
	  else {
	    TriCell* pTC = dynamic_cast<TriCell*>(pE);
	    assert(pTC);
	    if (spTriC.find(pTC) == spTriC.end())
	      spTriC.insert(pTC);
	    else
	      return;
	  }
	  break;
	}
      case iMesh_QUADRILATERAL:
	{
	  QuadFace *pQF = dynamic_cast<QuadFace*>(pE);
	  if (pQF) {
	    if (spQuadF.find(pQF) == spQuadF.end())
	      spQuadF.insert(pQF);
	    else
	      return;
	  }
	  else {
	    QuadCell* pQC = dynamic_cast<QuadCell*>(pE);
	    assert(pQC);
	    if (spQuadC.find(pQC) == spQuadC.end())
	      spQuadC.insert(dynamic_cast<QuadCell*>(pE));
	    else
	      return;
	  }
	  break;
	}
      case iMesh_TETRAHEDRON:
	{
	  TetCell* pTC = dynamic_cast<TetCell*>(pE);
	  assert(pTC);
	  if (spTetC.find(pTC) == spTetC.end())
	    spTetC.insert(pTC);
	  else
	    return;
	  break;
	}
      case iMesh_PYRAMID:
	{
	  PyrCell* pPyC = dynamic_cast<PyrCell*>(pE);
	  assert(pPyC);
	  if (spPyrC.find(pPyC) == spPyrC.end())
	    spPyrC.insert(pPyC);
	  else
	    return;
	  break;
	}
      case iMesh_PRISM:
	{
	  PrismCell* pPrC = dynamic_cast<PrismCell*>(pE);
	  assert(pPrC);
	  if (spPrismC.find(pPrC) == spPrismC.end())
	    spPrismC.insert(pPrC);
	  else
	    return;
	  break;
	}
      case iMesh_HEXAHEDRON:
	{
	  HexCell* pHC = dynamic_cast<HexCell*>(pE);
	  assert(pHC);
	  if (spHexC.find(pHC) == spHexC.end())
	    spHexC.insert(pHC);
	  else
	    return;
	  break;
	}
      case iMesh_POLYGON:
      case iMesh_POLYHEDRON:
      case iMesh_SEPTAHEDRON:
      default:
	// How did somebody create an entity of this topology?
	throw(iBase_INVALID_ARGUMENT);
      }

    // Ensure that the verts for each entity will be included, at
    // least in shadow form.
    for (int i = 0; i < pE->getNumVerts(); i++)
      vAddVertex(pE->getVert(i), true);
  }

  void
  EntitySet::vRemoveEntity(Entity* const pE)

  {
    // It is possible that the entity was never added to the set, in
    // which case informing the iterators does nothing...

    // Inform all of the iterators over this set of the removal.
    std::set<EntitySetIter*>::iterator itersIter = spIters.begin(),
	itersIterEnd = spIters.end();
    for (; itersIter != itersIterEnd; itersIter++)
      (*itersIter)->vDelete(pE);

    // Inform all of the worksets over this set of the removal.
    std::set<EntitySetWorkset*>::iterator setsIter = spWorksets.begin(),
	setsIterEnd = spWorksets.end();
    for (; setsIter != setsIterEnd; setsIter++)
      (*setsIter)->vDelete(pE);

    switch (pE->getEntTopology())
      {
      case iMesh_POINT:
	{
	  assert(dynamic_cast<Vert*>(pE));
	  vRemoveVertex(dynamic_cast<Vert*>(pE), false);
	  return;
	}
      case iMesh_LINE_SEGMENT:
	{
	  EdgeFace* pEF = dynamic_cast<EdgeFace*>(pE);
	  assert(pEF);
	  std::set<EdgeFace*>::iterator EFIter = spEdgeF.find(pEF);
	  if (EFIter == spEdgeF.end())
	    return;
	  else
	    spEdgeF.erase(EFIter);
	  break;
	}
      case iMesh_TRIANGLE:
	{
	  TriFace *pTF = dynamic_cast<TriFace*>(pE);
	  if (pTF) {
	    std::set<TriFace*>::iterator TFIter = spTriF.find(pTF);
	    if (TFIter == spTriF.end())
	      return;
	    else
	      spTriF.erase(TFIter);
	  }
	  else {
	    TriCell* pTC = dynamic_cast<TriCell*>(pE);
	    assert(pTC);
	    std::set<TriCell*>::iterator TCIter = spTriC.find(pTC);
	    if (TCIter == spTriC.end())
	      return;
	    else
	      spTriC.erase(TCIter);
	  }
	  break;
	}
      case iMesh_QUADRILATERAL:
	{
	  QuadFace *pQF = dynamic_cast<QuadFace*>(pE);
	  if (pQF) {
	    std::set<QuadFace*>::iterator QFIter = spQuadF.find(pQF);
	    if (QFIter == spQuadF.end())
	      return;
	    else
	      spQuadF.erase(QFIter);
	  }
	  else {
	    QuadCell* pQC = dynamic_cast<QuadCell*>(pE);
	    assert(pQC);
	    std::set<QuadCell*>::iterator QCIter = spQuadC.find(pQC);
	    if (QCIter == spQuadC.end())
	      return;
	    else
	      spQuadC.erase(QCIter);
	  }
	  break;
	}
      case iMesh_TETRAHEDRON:
	{
	  TetCell* pTC = dynamic_cast<TetCell*>(pE);
	  assert(pTC);
	  std::set<TetCell*>::iterator TCIter = spTetC.find(pTC);
	  if (TCIter == spTetC.end())
	    return;
	  else
	    spTetC.erase(TCIter);
	  break;
	}
      case iMesh_PYRAMID:
	{
	  PyrCell* pPyC = dynamic_cast<PyrCell*>(pE);
	  assert(pPyC);
	  std::set<PyrCell*>::iterator PyCIter = spPyrC.find(pPyC);
	  if (PyCIter == spPyrC.end())
	    return;
	  else
	    spPyrC.erase(PyCIter);
	  break;
	}
      case iMesh_PRISM:
	{
	  PrismCell* pPrC = dynamic_cast<PrismCell*>(pE);
	  assert(pPrC);
	  std::set<PrismCell*>::iterator PrCIter = spPrismC.find(pPrC);
	  if (PrCIter == spPrismC.end())
	    return;
	  else
	    spPrismC.erase(PrCIter);
	  break;
	}
      case iMesh_HEXAHEDRON:
	{
	  HexCell* pHC = dynamic_cast<HexCell*>(pE);
	  assert(pHC);
	  std::set<HexCell*>::iterator HCIter = spHexC.find(pHC);
	  if (HCIter == spHexC.end())
	    return;
	  else
	    spHexC.erase(HCIter);
	  break;
	}
      case iMesh_POLYGON:
      case iMesh_POLYHEDRON:
      case iMesh_SEPTAHEDRON:
	// How did somebody create an entity of this topology?
	throw(iBase_NOT_SUPPORTED);
	break;
      default:
	// How did somebody create an entity of this topology?
	throw(iBase_INVALID_ARGUMENT);
      } // Done with switch

    // Ensure that the verts for each entity will be removed, even
    // though they may not be members of the entity set.
    for (int i = 0; i < pE->getNumVerts(); i++)
      vRemoveVertex(pE->getVert(i), true);
  }

  void
  EntitySet::vRemoveEntities(const iBase_EntityHandle entity_handles[],
			     const int iNumEnt)

  {
    for (int iEH = 0; iEH < iNumEnt; iEH++) {
      Entity *pE = static_cast<Entity*>(entity_handles[iEH]);
      vRemoveEntity(pE);
    }
  }

  bool
  EntitySet::qIsEntityInSet(void* entity_handle) const

  {
    Entity *pE = static_cast<Entity*>(entity_handle);
    switch (pE->getEntTopology())
      {
      case iMesh_POINT:
	{
	  assert(dynamic_cast<Vert*>(pE));
	  return ::qIsInSet(sVR, VertRec(dynamic_cast<Vert*>(pE)));
	  // no break needed
	}
      case iMesh_LINE_SEGMENT:
	{
	  assert(dynamic_cast<EdgeFace*>(pE));
	  return ::qIsInSet(spEdgeF, dynamic_cast<EdgeFace*>(pE));
	  // no break needed
	}
      case iMesh_TRIANGLE:
	{
	  TriFace *pTF = dynamic_cast<TriFace*>(pE);
	  if (pTF) {
	    return ::qIsInSet(spTriF, pTF);
	  }
	  else {
	    assert(dynamic_cast<TriCell*>(pE));
	    return ::qIsInSet(spTriC, dynamic_cast<TriCell*>(pE));
	  }
	  break;
	}
      case iMesh_QUADRILATERAL:
	{
	  QuadFace *pQF = dynamic_cast<QuadFace*>(pE);
	  if (pQF) {
	    return ::qIsInSet(spQuadF, pQF);
	  }
	  else {
	    assert(dynamic_cast<QuadCell*>(pE));
	    return ::qIsInSet(spQuadC, dynamic_cast<QuadCell*>(pE));
	  }
	  break;
	}
      case iMesh_TETRAHEDRON:
	{
	  assert(dynamic_cast<TetCell*>(pE));
	  return ::qIsInSet(spTetC, dynamic_cast<TetCell*>(pE));
	  // no break needed
	}
      case iMesh_PYRAMID:
	{
	  assert(dynamic_cast<PyrCell*>(pE));
	  return ::qIsInSet(spPyrC, dynamic_cast<PyrCell*>(pE));
	}
      case iMesh_PRISM:
	{
	  assert(dynamic_cast<PrismCell*>(pE));
	  return ::qIsInSet(spPrismC, dynamic_cast<PrismCell*>(pE));
	  break;
	}
      case iMesh_HEXAHEDRON:
	{
	  assert(dynamic_cast<HexCell*>(pE));
	  return ::qIsInSet(spHexC, dynamic_cast<HexCell*>(pE));
	  break;
	}
      case iMesh_POLYGON:
      case iMesh_POLYHEDRON:
      case iMesh_SEPTAHEDRON:
	// How did somebody create an entity of this topology?
	throw(iBase_INVALID_ENTITY_TYPE);
	return false;
	break;
      default:
	throw(iBase_INVALID_ARGUMENT);
	return false;
	break;
      }
  }

  void
  EntitySet::vRebuildAuxiliaryVertList()
  {
    {
      std::set<EdgeFace*>::iterator iter = spEdgeF.begin();
      std::set<EdgeFace*>::iterator iterEnd = spEdgeF.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 2; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<TriFace*>::iterator iter = spTriF.begin();
      std::set<TriFace*>::iterator iterEnd = spTriF.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 3; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<QuadFace*>::iterator iter = spQuadF.begin();
      std::set<QuadFace*>::iterator iterEnd = spQuadF.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 4; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<TriCell*>::iterator iter = spTriC.begin();
      std::set<TriCell*>::iterator iterEnd = spTriC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 3; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<QuadCell*>::iterator iter = spQuadC.begin();
      std::set<QuadCell*>::iterator iterEnd = spQuadC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 4; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<TetCell*>::iterator iter = spTetC.begin();
      std::set<TetCell*>::iterator iterEnd = spTetC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 4; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<PyrCell*>::iterator iter = spPyrC.begin();
      std::set<PyrCell*>::iterator iterEnd = spPyrC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 5; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<PrismCell*>::iterator iter = spPrismC.begin();
      std::set<PrismCell*>::iterator iterEnd = spPrismC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 6; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    {
      std::set<HexCell*>::iterator iter = spHexC.begin();
      std::set<HexCell*>::iterator iterEnd = spHexC.end();
      while (iter != iterEnd) {
	for (int iV = 0; iV < 8; iV++)
	  vAddVertex((*iter)->getVert(iV), true);
	++iter;
      }
    }
    qVertNumsOK = false;
  }
}
