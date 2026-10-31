#include <map>
#include <string>
#include <algorithm>
#include <unistd.h>

#include "iMesh_GRUMMP.hh"
#include "iMesh_GRUMMP_misc.hh"
#include "GR_Mesh.h"
#include "iMesh_GRUMMP_Tag.hh"
#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"

static bool
isRootSet(iMesh_Base*, iBase_EntitySetHandle entity_set_handle)
{
  return (entity_set_handle == GR_rootSet);
}

using namespace ITAPS_GRUMMP;

template<class T>
  static int
  iGetGenericData(iBase_Taggable* entity_handle, TagBase* ptbTheTag,
		  char*& genericData)

  {
    int numValues = ptbTheTag->iGetSize();
    genericData = new char[numValues * sizeof(T)];
    Tag<T>* properTag = dynamic_cast<Tag<T>*>(ptbTheTag);
    assert(properTag);
    if (!properTag->getTagValue(entity_handle,
				reinterpret_cast<T*>(genericData))) {
      delete[] genericData;
      genericData = NULL;
      throw(iBase_TAG_NOT_FOUND);
      return 0;
    }
    return numValues * sizeof(T);
  }

template<class T>
  static T
  tGetSpecifData(iBase_Taggable* entity_handle, iBase_TagHandle tag_handle,
		 iBase_TagValueType theType)

  {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    T result;
    if (ptbTheTag->iGetTagType() != theType || ptbTheTag->iGetSize() != 1) {
      throw(iBase_INVALID_TAG_HANDLE);
    }
    else {
      Tag<T>* properTag = dynamic_cast<Tag<T>*>(ptbTheTag);
      assert(properTag);
      if (!properTag->getTagValue(entity_handle, &result))
	throw(iBase_TAG_NOT_FOUND);
    }
    return result;
  }

template<class T>
  static void
  vSetSpecifData(iBase_Taggable* entity_handle, iBase_TagHandle tag_handle,
		 iBase_TagValueType theType, T& tag_value)

  {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    if (ptbTheTag->iGetTagType() != theType || ptbTheTag->iGetSize() != 1) {
      throw(iBase_INVALID_TAG_HANDLE);
    }
    else {
      Tag<T>* properTag = dynamic_cast<Tag<T>*>(ptbTheTag);
      assert(properTag);
      properTag->setTagValue(entity_handle, &tag_value);
    }
  }

template<class T>
  static int
  iGetSpecifArrData(const iBase_EntityHandle entity_handles[], int arrSize,
		    iBase_TagHandle tag_handle, iBase_TagValueType theType,
		    T*& values)

  {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    if (ptbTheTag->iGetTagType() != theType) {
      throw(iBase_INVALID_TAG_HANDLE);
    }
    else {
      int numValues = ptbTheTag->iGetSize();
      values = new T[numValues * arrSize];
      Tag<T>* properTag = dynamic_cast<Tag<T>*>(ptbTheTag);
      assert(properTag);
      int lastIndex = properTag->getTagValueArr(entity_handles, arrSize,
						values);
      if (lastIndex != arrSize) {
	delete[] values;
	throw(iBase_INVALID_ENTITY_HANDLE);
      }
      return numValues * arrSize;
    }
  }

template<class T>
  static void
  vSetSpecifArrData(const iBase_EntityHandle entity_handles[], int numEnts,
		    iBase_TagHandle tag_handle, iBase_TagValueType theType,
		    const T values[], int numValuesIn)

  {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    if (ptbTheTag->iGetTagType() != theType) {
      throw(iBase_INVALID_TAG_HANDLE);
    }
    else {
      int numValues = ptbTheTag->iGetSize();
      if (numEnts * numValues != numValuesIn) {
	throw(iBase_BAD_ARRAY_SIZE);
      }
      else {
	Tag<T>* properTag = dynamic_cast<Tag<T>*>(ptbTheTag);
	assert(properTag);
	properTag->setTagValueArr(entity_handles, numEnts, values);
      }
    }
  }

iMesh_Base::iMesh_Base(const enum eMesh_Type eMT_in) :
    pM(NULL), externalData(false), m_err(iBase_SUCCESS), eMT(eMT_in), spEntSets(), spSetWorksets(), spSetIters(), spTags()
{
  iMessageStdoutLevel = 1;
}

iMesh_Base*
iMesh_Base::Create(const enum eMesh_Type eGT_In)
{
  switch (eGT_In)
    {
    case ePlane:
      return new iMesh_Plane;
    case eSurface:
      return new iMesh_Surface;
    case eVolume:
      return new iMesh_Volume;
    default:
      // Impossible case
      assert(0);
      return NULL;
    }
}

void
iMesh_Base::vDestroySetsTagsIters()

{
  { // Delete all entity sets and clear the set
    std::set<EntitySetBase*>::iterator iter, iterEnd = spEntSets.end();
    for (iter = spEntSets.begin(); iter != iterEnd; ++iter)
      delete (*iter);
    spEntSets.clear();
  }
  { // Delete all entity iterators and clear the set
    std::set<EntIter*>::iterator iter, iterEnd = spSetIters.end();
    for (iter = spSetIters.begin(); iter != iterEnd; ++iter)
      delete (*iter);
    spSetIters.clear();
  }
  { // Delete all worksets and clear the set
    std::set<Workset*>::iterator iter, iterEnd = spSetWorksets.end();
    for (iter = spSetWorksets.begin(); iter != iterEnd; ++iter)
      delete (*iter);
    spSetWorksets.clear();
  }
  { // Delete all remaining tags
    std::set<TagBase*>::iterator iter, iterEnd = spTags.end();
    for (iter = spTags.begin(); iter != iterEnd; ++iter)
      delete (*iter);
  }
}

void
iMesh_Base::destroyEntSet(/*in*/iBase_EntitySetHandle entity_set)

{
  // First ensure the arg is not the 'root set'
  if (isRootSet(this, entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  // 'Else' is not necessary if the error-handling method is to throw, but
  // is necessary if the error-handling method is silent, etc.
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    spEntSets.erase(pESB);
    delete pESB;
  }
}

bool
iMesh_Base::isList(/*in*/iBase_EntitySetHandle entity_set)

{
  if (isRootSet(this, entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase* pESB = static_cast<EntitySetBase*>(entity_set);
    // The entity set itself can tell us if it is ordered
    return pESB->qIsOrdered();
  }
}

int
iMesh_Base::getNumEntSets(/*in*/iBase_EntitySetHandle entity_set_handle,
/*in*/int num_hops)

{
  if (isRootSet(this, entity_set_handle)) {
    return spEntSets.size();
  }
  else {
    // If we get here, it is not the root set
    EntitySetBase* pESB = static_cast<EntitySetBase*>(entity_set_handle);
    std::set<EntitySetBase*> spESRes;
    pESB->vGetEntitySets(spESRes, num_hops);
    return spESRes.size();
  }
}

void
iMesh_Base::getEntSets(
    /*in*/iBase_EntitySetHandle entity_set_handle, /*in*/
    int num_hops,
    /*inout*/iBase_EntitySetHandle** contained_entity_set_handles,
    /*inout*/int* contained_entity_set_handles_allocated,
    /*out*/int *contained_entity_set_handles_size)

{
  if (isRootSet(this, entity_set_handle)) {
    // The mesh case.  All entity sets are contained in the root set,
    // so the number of hops is irrelevant
    *contained_entity_set_handles_size = spEntSets.size();
    vEnsureArraySize<iBase_EntitySetHandle>(
	contained_entity_set_handles, contained_entity_set_handles_allocated,
	*contained_entity_set_handles_size);
    std::set<EntitySetBase*>::iterator iter = spEntSets.begin(), iterEnd =
	spEntSets.end();
    for (int iES = 0; iter != iterEnd; iES++, ++iter) {
      (*contained_entity_set_handles)[iES] = *iter;
    }
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set_handle);
    // Entity sets have their own method for this
    pESB->vGetEntitySets(contained_entity_set_handles,
			 contained_entity_set_handles_allocated,
			 contained_entity_set_handles_size, num_hops);
  }
}

void
iMesh_Base::addEntToSet(/*in*/iBase_EntityHandle entity_handle,
/*inout*/iBase_EntitySetHandle entity_set)

{
  if (isRootSet(this, entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    // Works for both sets and lists
    pESB->vAddEntity(entity_handle);
  }
}

void
iMesh_Base::addEntArrToSet(/*in*/const iBase_EntityHandle entity_handles[],
/*in*/int entity_handles_size,
			   /*inout*/iBase_EntitySetHandle entity_set)

{
#ifndef NDEBUG
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
#endif
  if (isRootSet(this, entity_set)) {
    // Cannot add an entity to the root set, is implicitly contained
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    pESB->vAddEntities(entity_handles, entity_handles_size);
  }
}

void
iMesh_Base::rmvEntFromSet(/*in*/iBase_EntityHandle entity_handle,
/*inout*/iBase_EntitySetHandle entity_set)

{
  if (isRootSet(this, entity_set)) {
    // Cannot remove an entity from the root set, must delete the entity instead
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    pESB->vRemoveEntity(entity_handle);
  }
}

void
iMesh_Base::rmvEntArrFromSet(
/*in*/const iBase_EntityHandle entity_handles[],
			     /*in*/int entity_handles_size,
			     /*inout*/iBase_EntitySetHandle entity_set)

{
#ifndef NDEBUG
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
#endif
  if (isRootSet(this, entity_set)) {
    // Same as for single entity removal
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    pESB->vRemoveEntities(entity_handles, entity_handles_size);
  }
}

void
iMesh_Base::addEntSet(/*in*/iBase_EntitySetHandle entity_set_to_add,
/*inout*/iBase_EntitySetHandle entity_set)

{
  if (isRootSet(this, entity_set) || isRootSet(this, entity_set_to_add)) {
    // Every entity set is implicitly contained in the root set
    throw(iBase_INVALID_ARGUMENT);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    pESB->vAddEntitySet(entity_set_to_add);
  }
}

void
iMesh_Base::rmvEntSet(/*in*/iBase_EntitySetHandle entity_set_to_remove,
/*inout*/iBase_EntitySetHandle entity_set)

{
  if (isRootSet(this, entity_set) || isRootSet(this, entity_set_to_remove)) {
    // Cannot remove an entity set, must delete it
    throw(iBase_INVALID_ARGUMENT);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    pESB->vRemoveEntitySet(entity_set_to_remove);
  }
}

bool
iMesh_Base::isEntContained(
/*in*/iBase_EntitySetHandle containing_entity_set,
			   /*in*/iBase_EntityHandle entity_handle)

{
  if (isRootSet(this, containing_entity_set)) {
    // The mesh interface contains all entities, so this had better
    // return true 100% of the time.
    return true;
  }
  else {
    EntitySetBase* pESBContainer =
	static_cast<EntitySetBase*>(containing_entity_set);
    return pESBContainer->qIsEntityInSet(entity_handle);
  }
}

bool
iMesh_Base::isEntSetContained(
/*in*/iBase_EntitySetHandle containing_entity_set,
			      /*in*/iBase_EntitySetHandle contained_entity_set)

{
  EntitySetBase* pESBContained =
      static_cast<EntitySetBase*>(contained_entity_set);
  if (isRootSet(this, containing_entity_set)) {
    // The mesh interface contains all entity sets, this had better
    // return true 100% of the time.
    return (spEntSets.find(pESBContained) != spEntSets.end());
  }
  else {
    EntitySetBase* pESBContainer =
	static_cast<EntitySetBase*>(containing_entity_set);
    return pESBContainer->qIsMyMember(pESBContained);
  }
}

void
iMesh_Base::addPrntChld(/*inout*/iBase_EntitySetHandle parent_entity_set,
/*inout*/iBase_EntitySetHandle child_entity_set)

{
  // Can't establish parent-child relationships with the root set
  if (isRootSet(this, parent_entity_set) || isRootSet(this, child_entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }

  EntitySetBase *pESB_P = static_cast<EntitySetBase*>(parent_entity_set);
  EntitySetBase *pESB_C = static_cast<EntitySetBase*>(child_entity_set);

  pESB_C->vAddParent(pESB_P);
  pESB_P->vAddChild(pESB_C);
}

void
iMesh_Base::rmvPrntChld(/*inout*/iBase_EntitySetHandle parent_entity_set,
/*inout*/iBase_EntitySetHandle child_entity_set)

{
  // Can't establish parent-child relationships with the root set,
  // so removing the relationship would make no sense
  if (isRootSet(this, parent_entity_set) || isRootSet(this, child_entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }

  EntitySetBase *pESB_P = static_cast<EntitySetBase*>(parent_entity_set);
  EntitySetBase *pESB_C = static_cast<EntitySetBase*>(child_entity_set);

  // If there was no relationship before, these calls do nothing
  pESB_P->vRemoveChild(pESB_C);
  pESB_C->vRemoveParent(pESB_P);
}

bool
iMesh_Base::isChildOf(/*in*/iBase_EntitySetHandle parent_entity_set,
/*in*/iBase_EntitySetHandle child_entity_set)

{

  // Can't establish parent-child relationships with the root set,
  // so checking the relationship would make no sense
  if (isRootSet(this, parent_entity_set) || isRootSet(this, child_entity_set)) {
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }

  EntitySetBase *pESB_P = static_cast<EntitySetBase*>(parent_entity_set);
  EntitySetBase *pESB_C = static_cast<EntitySetBase*>(child_entity_set);
  return (pESB_P->qIsMyChild(pESB_C));
}

int
iMesh_Base::getNumChld(/*in*/iBase_EntitySetHandle entity_set, /*in*/
		       int num_hops)

{
  if (isRootSet(this, entity_set)) {
    // Can't have parent/child relationships involving the root set.
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    return pESB->iGetNumChildren(num_hops);
  }
}

int
iMesh_Base::getNumPrnt(/*in*/iBase_EntitySetHandle entity_set, /*in*/
		       int num_hops)

{
  if (isRootSet(this, entity_set)) {
    // Can't have parent/child relationships involving the root set.
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set);
    return pESB->iGetNumParents(num_hops);
  }
}

void
iMesh_Base::getChldn(/*in*/iBase_EntitySetHandle from_entity_set,
/*in*/int num_hops,
		     /*inout*/iBase_EntitySetHandle** entity_set_handles,
		     /*inout*/int *entity_set_handles_allocated,
		     /*out*/int *entity_set_handles_size)

{
  if (isRootSet(this, from_entity_set)) {
    // Can't have parent/child relationships involving the root set.
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase* pESB = static_cast<EntitySetBase*>(from_entity_set);
    pESB->vGetChildren(num_hops, entity_set_handles,
		       entity_set_handles_allocated, entity_set_handles_size);
  }
}

void
iMesh_Base::getPrnts(/*in*/iBase_EntitySetHandle from_entity_set,
/*in*/int num_hops,
		     /*inout*/iBase_EntitySetHandle** entity_set_handles,
		     /*inout*/int* entity_set_handles_allocated,
		     /*out*/int* entity_set_handles_size)

{
  if (isRootSet(this, from_entity_set)) {
    // Can't have parent/child relationships involving the root set.
    throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  else {
    EntitySetBase* pESB = static_cast<EntitySetBase*>(from_entity_set);
    pESB->vGetParents(num_hops, entity_set_handles,
		      entity_set_handles_allocated, entity_set_handles_size);
  }
}

void
iMesh_Base::vCopyRootSet(/*out*/EntitySetBase*& pESBResult,
/*in*/const EntitySetBase* const pESBCompare)

{
  // Copies all of the entities in the mesh into pESBResult, except for the
  // entities that are already contained in pESBCompare
  pESBResult = new EntitySet();
  int iNVert_db = pM->getNumVerts();
  int iNFace_db = pM->getNumFaces();
  int iNCell_db = pM->getNumCells();
  iBase_EntityHandle entity_handle;
  // All verts
  int iE = 0;
  // Skip to the first non-deleted vert
  while (iE < iNVert_db && pM->getVert(iE)->isDeleted())
    iE++;

  while (iE < iNVert_db) {
    entity_handle = pM->getVert(iE);
    if ((pESBCompare == NULL) || (!pESBCompare->qIsEntityInSet(entity_handle)))
      pESBResult->vAddEntity(entity_handle);
    iE++;
    while (iE < iNVert_db && pM->getVert(iE)->isDeleted())
      iE++;
  }
  // All edge faces
  iE = 0;
  // Skip to the first non-deleted face
  while (iE < iNFace_db && pM->getFace(iE)->isDeleted())
    iE++;

  while (iE < iNFace_db) {
    entity_handle = pM->getFace(iE);
    if ((pESBCompare == NULL) || (!pESBCompare->qIsEntityInSet(entity_handle)))
      pESBResult->vAddEntity(entity_handle);
    iE++;
    while (iE < iNFace_db && pM->getFace(iE)->isDeleted())
      iE++;
  }
  // All cells
  iE = 0;
  // Skip to the first non-deleted cell
  while (iE < iNCell_db && pM->getCell(iE)->isDeleted())
    iE++;

  while (iE < iNCell_db) {
    entity_handle = pM->getCell(iE);
    if ((pESBCompare == NULL) || (!pESBCompare->qIsEntityInSet(entity_handle)))
      pESBResult->vAddEntity(entity_handle);
    iE++;
    while (iE < iNCell_db && pM->getCell(iE)->isDeleted())
      iE++;
  }
  // All entity sets
  std::set<EntitySetBase*>::iterator iter = spEntSets.begin(), iterEnd =
      spEntSets.end();
  for (; iter != iterEnd; ++iter)
    if ((pESBCompare == NULL) || (!pESBCompare->qIsMyMember(*iter)))
      pESBResult->vAddEntitySet(*iter);
}

void
iMesh_Base::vCopyRootSet( /*inout*/iBase_EntityHandle** entity_handles,
/*inout*/int* entity_handles_allocated,
			 /*out*/int* entity_handles_size)

{
  // Copies all of the entities in the mesh into the array
  int iNVert_db = pM->getNumVerts();
  int iNFace_db = pM->getNumFaces();
  int iNCell_db = pM->getNumCells();
  int iNEnt_db = iNVert_db + iNFace_db + iNCell_db;

  vEnsureArraySize<iBase_EntityHandle>(entity_handles, entity_handles_allocated,
				       iNEnt_db);

  int iEH = 0;
  { // All verts
    int iV = 0;
    // Skip to the first non-deleted vert
    while (iV < iNVert_db && pM->getVert(iV)->isDeleted())
      iV++;

    while (iV < iNVert_db) {
      (*entity_handles)[iEH] = pM->getVert(iV);
      iEH++;
      iV++;
      while (iV < iNVert_db && pM->getVert(iV)->isDeleted())
	iV++;
    }
  }

  { // All edge faces
    int iF = 0;
    // Skip to the first non-deleted face
    while (iF < iNFace_db && pM->getFace(iF)->isDeleted())
      iF++;

    while (iF < iNFace_db) {
      (*entity_handles)[iEH] = pM->getFace(iF);
      iEH++;
      iF++;
      while (iF < iNFace_db && pM->getFace(iF)->isDeleted())
	iF++;
    }
  }

  { // All cells
    int iCell = 0;
    // Skip to the first non-deleted cell
    while (iCell < iNCell_db && pM->getCell(iCell)->isDeleted())
      iCell++;

    while (iCell < iNCell_db) {
      (*entity_handles)[iEH] = pM->getCell(iCell);
      iEH++;
      iCell++;
      while (iCell < iNCell_db && pM->getCell(iCell)->isDeleted())
	iCell++;
    }
  }
  *entity_handles_size = iEH;
}

void
iMesh_Base::subtract(/*in*/iBase_EntitySetHandle entity_set_1,
/*in*/iBase_EntitySetHandle entity_set_2,
		     /*inout*/iBase_EntitySetHandle* entity_set_result)

{
  // Subtracts entity_set_2 from entity_set_1 and creates a new entity
  // set to store the result
  EntitySetBase *pESBResult = NULL;
  if (isRootSet(this, entity_set_1)) {
    // Copy every entity in the mesh except what is in entity_set_2 into the
    // resulting entity set
    EntitySetBase *pESB2 = static_cast<EntitySetBase*>(entity_set_2);
    vCopyRootSet(pESBResult, pESB2);
  }
  else if (isRootSet(this, entity_set_2)) {
    // The result has to be an empty entity set because the root set contains
    // all the entities, so at least the entities in entity_set_1
    pESBResult = new EntitySet();
  }
  else {
    // Call the method from SetBoolOps
    EntitySetBase *pESB1 = static_cast<EntitySetBase*>(entity_set_1);
    EntitySetBase *pESB2 = static_cast<EntitySetBase*>(entity_set_2);
    vSubtract(pESB1, pESB2, pESBResult);
  }
  // Keep track of the new set, makes validity checks easier
  spEntSets.insert(pESBResult);
  *entity_set_result = pESBResult;
}

void
iMesh_Base::intersect(/*in*/iBase_EntitySetHandle entity_set_1,
/*in*/iBase_EntitySetHandle entity_set_2,
		      /*inout*/iBase_EntitySetHandle* entity_set_result)

{
  // Intersects entity_set_2 with entity_set_1 and creates a new entity
  // set to store the result

  EntitySetBase *pESBResult = NULL;
  if (isRootSet(this, entity_set_1)) {
    // Create an unordered copy of entity_set_2 since the root set is unordered
    // any intersection involving a set will result in a set
    EntitySetBase *pESB2 = static_cast<EntitySetBase*>(entity_set_2);
    if (pESB2->qIsOrdered())
      pESBResult = new EntitySet(*dynamic_cast<EntityList*>(pESB2));
    else
      pESBResult = new EntitySet(*dynamic_cast<EntitySet*>(pESB2));
  }
  else if (isRootSet(this, entity_set_2)) {
    // Create an unordered copy of entity_set_1 since the root set is unordered
    // any intersection involving a set will result in a set
    EntitySetBase *pESB1 = static_cast<EntitySetBase*>(entity_set_1);
    if (pESB1->qIsOrdered())
      pESBResult = new EntitySet(*dynamic_cast<EntityList*>(pESB1));
    else
      pESBResult = new EntitySet(*dynamic_cast<EntitySet*>(pESB1));
  }
  else {
    // Call the method from SetBoolOps
    EntitySetBase *pESB1 = static_cast<EntitySetBase*>(entity_set_1);
    EntitySetBase *pESB2 = static_cast<EntitySetBase*>(entity_set_2);
    vIntersect(pESB1, pESB2, pESBResult);
  }
  // Keep track of the new set, makes validity checks easier
  spEntSets.insert(pESBResult);
  *entity_set_result = pESBResult;
}

void
iMesh_Base::unite(/*in*/iBase_EntitySetHandle entity_set_1,
/*in*/iBase_EntitySetHandle entity_set_2,
		  /*inout*/iBase_EntitySetHandle* entity_set_result)

{

  EntitySetBase* pESBResult = NULL;
  if (isRootSet(this, entity_set_1) || isRootSet(this, entity_set_2)) {
    // Copy all of the entities in the root set
    vCopyRootSet(pESBResult);
  }
  else {
    // Call the method from SetBoolOps
    EntitySetBase *pESB1 = static_cast<EntitySetBase*>(entity_set_1);
    EntitySetBase *pESB2 = static_cast<EntitySetBase*>(entity_set_2);
    vUnite(pESB1, pESB2, pESBResult);
  }
  // Keep track of the new set, makes validity checks easier
  spEntSets.insert(pESBResult);
  *entity_set_result = pESBResult;
}

void
iMesh_Base::vRemoveAllReferencesTo(Entity* const pE)

{
  // Intended to remove an entity from all entity sets.
  {
    ::std::set<EntitySetBase*>::iterator iter;
    for (iter = spEntSets.begin(); iter != spEntSets.end(); ++iter) {
      (*iter)->vRemoveEntity(pE);
    }
  }

  // Access and remove all tags on this entity
  iBase_TagHandle* tag_handles = NULL;
  int tag_handles_allocated = 0, tag_handles_size;
  getAllTags(pE, &tag_handles, &tag_handles_allocated, &tag_handles_size);
  for (int i = 0; i < tag_handles_size; i++) {
    rmvTag(pE, tag_handles[i]);
  }
  if (tag_handles)
    free(tag_handles);

  // Find all cases where this entity is a tag value and remove those, too.
  {
    std::set<TagBase*>::iterator iter = spTags.begin(), iterEnd = spTags.end();
    for (; iter != iterEnd; ++iter) {
      if (getTagType(*iter) == iBase_ENTITY_HANDLE) {
	Tag<iBase_EntityHandle>* pt =
	    dynamic_cast<Tag<iBase_EntityHandle>*>(*iter);
	assert(pt);
	iBase_EntityHandle EH = pE;
	pt->vRemoveTagValue(&EH);
      }
    }
  }
}

/// Order faces based on lexicographic ordering of their vertices.
class FaceComp : public std::binary_function<Face*, Face*, bool> {
public:
  bool
  operator()(const Face* const pFA, const Face* const pFB) const
  {
    const Vert *apVA[4], *apVB[4];
    int iNVA = pFA->getNumVerts();
    int iNVB = pFB->getNumVerts();
    int i;

    for (i = 0; i < iNVA; i++) {
      apVA[i] = pFA->getVert(i);
    }
    for (i = 0; i < iNVB; i++) {
      apVB[i] = pFB->getVert(i);
    }

    std::sort(apVA, apVA + iNVA);
    std::sort(apVB, apVB + iNVB);

    for (i = 0; i < std::min(iNVA, iNVB); i++) {
      if (apVA[i] < apVB[i])
	return true;
      else if (apVA[i] > apVB[i])
	return false;
    }
    // If the lengths are the same, A < B if A has -fewer- vertices
    // than B.  So for identical faces, A !< B and B !< A, as
    // required.
    return (iNVA < iNVB);
  }
};

#include "iMesh_GRUMMP_Tag.hh"

bool
iMesh_Base::qFindTag(const char name[], const int name_len, TagBase*& ptbResult)

// Finds a tag handle using its name, returns true if the tag is found
// and sets ptbResult to its memory location, otherwise returns false
{
  // All valid tag handles are store in spTags, so check there
  std::set<TagBase*>::iterator iter = spTags.begin(), iterEnd = spTags.end();
  for (; iter != iterEnd; ++iter) {
    char thisName[120];
    int thisNameLen = 120;
    (*iter)->vGetTagName(thisName, thisNameLen);
    if (strncmp(name, thisName, name_len) == 0) {
      ptbResult = (*iter);
      return true;
    }
  }
  return false;
}

void
iMesh_Base::createTag(const char* strName, const int tag_name_len,
		      const int tagSize, const int type,
		      iBase_TagHandle* tagHandle)

{
  try {
    TagBase* ptbTheTag;
    // Check to make sure a tag with that name does not already exist
    if (qFindTag(strName, tag_name_len, ptbTheTag)) {
      throw(iBase_TAG_ALREADY_EXISTS);
    }
    // Must have at least one value stored
    if (tagSize < 1) {
      throw(iBase_INVALID_ARGUMENT);
    }
    switch (type)
      {
      case iBase_INTEGER:
	{
	  ptbTheTag = new Tag<int>(strName, tagSize, type);
	}
	break;

      case iBase_DOUBLE:
	{
	  ptbTheTag = new Tag<double>(strName, tagSize, type);
	}
	break;

      case iBase_ENTITY_HANDLE:
	{
	  ptbTheTag = new Tag<iBase_EntityHandle>(strName, tagSize, type);
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  ptbTheTag = new Tag<iBase_EntitySetHandle>(strName, tagSize, type);
	}
	break;

      case iBase_BYTES:
	{
	  ptbTheTag = new Tag<char>(strName, tagSize, type);
	}
	break;

      default:
	{
	  throw(iBase_INVALID_ARGUMENT);
	}
	break;
      }
    // Keep track of the tag
    spTags.insert(ptbTheTag);
    *tagHandle = ptbTheTag;
  }
  // Just in case there wasn't enough memory to create the tag
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::destroyTag(iBase_TagHandle tagHandle, bool forced)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tagHandle);
  // Ensure the validity of the tag handle
  // Cannot delete a tag that is in use without setting forced to true
  if (!forced && !ptbTheTag->qIsEmpty()) {
    throw(iBase_TAG_IN_USE);
  }
  // First remove it from the set of valid tag handles
  spTags.erase(ptbTheTag);
  delete ptbTheTag;
}

void
iMesh_Base::getTagName(iBase_TagHandle tagHandle, char* name, int name_len)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tagHandle);
  ptbTheTag->vGetTagName(name, name_len);
}

int
iMesh_Base::getTagSizeValues(iBase_TagHandle tagHandle)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tagHandle);
  return ptbTheTag->iGetSize();
}

int
iMesh_Base::getTagSizeBytes(iBase_TagHandle tagHandle)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tagHandle);
  switch (ptbTheTag->iGetTagType())
    {
    case iBase_INTEGER:
      return sizeof(int) * (ptbTheTag->iGetSize());
    case iBase_DOUBLE:
      return sizeof(double) * (ptbTheTag->iGetSize());
    case iBase_ENTITY_HANDLE:
      return sizeof(iBase_EntityHandle) * (ptbTheTag->iGetSize());
    case iBase_ENTITY_SET_HANDLE:
      return sizeof(iBase_EntitySetHandle) * (ptbTheTag->iGetSize());
    case iBase_BYTES:
      return ptbTheTag->iGetSize();
    default:
      // Should never be able to get here, because you can't create a tag
      // with any other type.
      assert(0);
      throw(iBase_INVALID_TAG_HANDLE);
    }
}

iBase_TagHandle
iMesh_Base::getTagHandle(const char *name, const int name_len)

{
  TagBase* ptbTheTag;
  if (qFindTag(name, name_len, ptbTheTag))
    return ptbTheTag;
  else {
    throw(iBase_TAG_NOT_FOUND);
  }
}

int
iMesh_Base::getTagType(iBase_TagHandle tagHandle)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tagHandle);
  return ptbTheTag->iGetTagType();
}

void
iMesh_Base::getData(iBase_EntityHandle entity_handle,
		    iBase_TagHandle tag_handle, char** tag_value,
		    int* tag_value_allocated, int* tag_value_size)

{
  // Supposed to be able to handle memory allocation failures
  try {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    char* genericData;
    int type = ptbTheTag->iGetTagType();

    // The templatized iGetGenericData method allocates the array of chars
    // necessary and returns the resulting size
    switch (type)
      {
      case iBase_INTEGER:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<int>(entity_handle, ptbTheTag,
						 genericData);
	}
	break;

      case iBase_DOUBLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<double>(entity_handle, ptbTheTag,
						    genericData);
	}
	break;

      case iBase_ENTITY_HANDLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<iBase_EntityHandle>(entity_handle,
								ptbTheTag,
								genericData);
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<iBase_EntitySetHandle>(
	      entity_handle, ptbTheTag, genericData);
	}
	break;

      case iBase_BYTES:
	{
	  // The simple case, the number of chars is the tag size
	  *tag_value_size = ptbTheTag->iGetSize();
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericData = new char[*tag_value_size];
	  if (!genericTag->getTagValue(entity_handle, genericData)) {
	    delete[] genericData;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;
      }
    vEnsureArraySize<char>(tag_value, tag_value_allocated, *tag_value_size);
    for (int i = 0; i < *tag_value_size; i++) {
      (*tag_value)[i] = genericData[i];
    }
    delete[] genericData;
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

int
iMesh_Base::getIntData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<int>(entity_handle, tag_handle, iBase_INTEGER);
}

double
iMesh_Base::getDblData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<double>(entity_handle, tag_handle, iBase_DOUBLE);
}

iBase_EntityHandle
iMesh_Base::getEHData(iBase_EntityHandle entity_handle,
		      iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<iBase_EntityHandle>(entity_handle, tag_handle,
					    iBase_ENTITY_HANDLE);
}

iBase_EntitySetHandle
iMesh_Base::getESHData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<iBase_EntitySetHandle>(entity_handle, tag_handle,
					       iBase_ENTITY_SET_HANDLE);
}

void
iMesh_Base::setData(iBase_EntityHandle entity_handle,
		    iBase_TagHandle tag_handle, const char* tag_value,
		    int tag_size)

{
  try {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    int type = ptbTheTag->iGetTagType();
    int iNumValues = ptbTheTag->iGetSize();
    switch (type)
      {
      case iBase_INTEGER:
	{
	  if (static_cast<size_t>(tag_size) != iNumValues * sizeof(int)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  Tag<int>* properTag = dynamic_cast<Tag<int>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_handle,
				 reinterpret_cast<const int*>(tag_value));
	  break;
	}

      case iBase_DOUBLE:
	{
	  if (static_cast<size_t>(tag_size) != iNumValues * sizeof(double)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  Tag<double>* properTag = dynamic_cast<Tag<double>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_handle,
				 reinterpret_cast<const double*>(tag_value));
	  break;
	}

      case iBase_ENTITY_HANDLE:
	{
	  if (static_cast<size_t>(tag_size)
	      != iNumValues * sizeof(iBase_EntityHandle)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
#ifndef NDEBUG
	  const iBase_EntityHandle* entity_handles =
	      reinterpret_cast<const iBase_EntityHandle*>(tag_value);
	  for (int i = 0; i < iNumValues; i++) {
	    iBase_EntityHandle EH_tmp = entity_handles[i];
	    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
	      throw(iBase_INVALID_ENTITY_HANDLE);
	  }
#endif
	  Tag<iBase_EntityHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntityHandle>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(
	      entity_handle,
	      reinterpret_cast<const iBase_EntityHandle*>(tag_value));
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  if (static_cast<size_t>(tag_size)
	      != iNumValues * sizeof(iBase_EntitySetHandle)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
#ifndef NDEBUG
	  const iBase_EntitySetHandle* entity_set_handles =
	      reinterpret_cast<const iBase_EntitySetHandle*>(tag_value);
	  for (int i = 0; i < iNumValues; i++) {
	    iBase_EntitySetHandle ESH_tmp = entity_set_handles[i];
	    if (ESH_tmp == NULL
		|| !qHasSet(dynamic_cast<EntitySetBase*>(ESH_tmp)))
	      throw(iBase_INVALID_ENTITYSET_HANDLE);
	  }
#endif
	  Tag<iBase_EntitySetHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntitySetHandle>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(
	      entity_handle,
	      reinterpret_cast<const iBase_EntitySetHandle*>(tag_value));
	}
	break;

      case iBase_BYTES:
	{
	  if (tag_size != iNumValues) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericTag->setTagValue(entity_handle, tag_value);
	}
	break;
      }
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::setIntData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle, int tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<int>(entity_handle, tag_handle, iBase_INTEGER, tag_value);
}

void
iMesh_Base::setDblData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle, double tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<double>(entity_handle, tag_handle, iBase_DOUBLE, tag_value);
}

void
iMesh_Base::setEHData(iBase_EntityHandle entity_handle,
		      iBase_TagHandle tag_handle, iBase_EntityHandle tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<iBase_EntityHandle>(entity_handle, tag_handle,
				     iBase_ENTITY_HANDLE, tag_value);
}

void
iMesh_Base::setESHData(iBase_EntityHandle entity_handle,
		       iBase_TagHandle tag_handle,
		       iBase_EntitySetHandle tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<iBase_EntitySetHandle>(entity_handle, tag_handle,
					iBase_ENTITY_SET_HANDLE, tag_value);
}

void
iMesh_Base::getAllTags(iBase_EntityHandle entity_handle,
		       iBase_TagHandle** tag_handles,
		       int *tag_handles_allocated, int *tag_handles_size)

{
  try {
    std::set<TagBase*>::iterator iter = spTags.begin(), iterEnd = spTags.end();
    // Create a vector to store the tags temporarily, since it is
    // impossible to know how many there will be
    std::vector<iBase_TagHandle> tagVect;
    for (; iter != iterEnd; ++iter)
      if ((*iter)->qIsTagged(entity_handle))
	tagVect.push_back(*iter);
    *tag_handles_size = tagVect.size();
    vEnsureArraySize<iBase_TagHandle>(tag_handles, tag_handles_allocated,
				      *tag_handles_size);
    for (int i = 0; i < *tag_handles_size; i++)
      (*tag_handles)[i] = tagVect[i];
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::rmvTag(iBase_EntityHandle entity_handle, iBase_TagHandle tag_handle)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
  ptbTheTag->qUntagEntity(entity_handle);
}

void
iMesh_Base::getArrData(const iBase_EntityHandle* entity_handles,
		       int entity_handles_size, iBase_TagHandle tag_handle,
		       char** tag_value, int* tag_value_allocated,
		       int* tag_value_size)

{
  try {
    for (int i = 0; i < entity_handles_size; i++) {
      iBase_EntityHandle EH_tmp = entity_handles[i];
      if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
	throw(iBase_INVALID_ENTITY_HANDLE);
    }
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    int type = ptbTheTag->iGetTagType();
    int iNumValues = ptbTheTag->iGetSize();
    int dataSizeInBytes;
    char* genericValues = NULL;
    int lastIndex;
    switch (type)
      {

      case iBase_INTEGER:
	{
	  dataSizeInBytes = iNumValues * sizeof(int);
	  *tag_value_size = entity_handles_size * dataSizeInBytes;
	  vEnsureArraySize<char>(tag_value, tag_value_allocated,
				 *tag_value_size);
	  Tag<int>* properTag = dynamic_cast<Tag<int>*>(ptbTheTag);
	  assert(properTag);
	  genericValues = new char[*tag_value_size];
	  lastIndex = properTag->getTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<int*>(genericValues));
	  if (lastIndex != entity_handles_size) {
	    delete[] genericValues;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;

      case iBase_DOUBLE:
	{
	  dataSizeInBytes = iNumValues * sizeof(double);
	  *tag_value_size = entity_handles_size * dataSizeInBytes;
	  vEnsureArraySize<char>(tag_value, tag_value_allocated,
				 *tag_value_size);
	  Tag<double>* properTag = dynamic_cast<Tag<double>*>(ptbTheTag);
	  assert(properTag);
	  genericValues = new char[*tag_value_size];
	  lastIndex = properTag->getTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<double*>(genericValues));
	  if (lastIndex != entity_handles_size) {
	    delete[] genericValues;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;

      case iBase_ENTITY_HANDLE:
	{
	  dataSizeInBytes = iNumValues * sizeof(iBase_EntityHandle*);
	  *tag_value_size = entity_handles_size * dataSizeInBytes;
	  vEnsureArraySize<char>(tag_value, tag_value_allocated,
				 *tag_value_size);
	  Tag<iBase_EntityHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntityHandle>*>(ptbTheTag);
	  assert(properTag);
	  genericValues = new char[*tag_value_size];
	  lastIndex = properTag->getTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<iBase_EntityHandle*>(genericValues));
	  if (lastIndex != entity_handles_size) {
	    delete[] genericValues;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  dataSizeInBytes = iNumValues * sizeof(iBase_EntitySetHandle*);
	  *tag_value_size = entity_handles_size * dataSizeInBytes;
	  vEnsureArraySize<char>(tag_value, tag_value_allocated,
				 *tag_value_size);
	  Tag<iBase_EntitySetHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntitySetHandle>*>(ptbTheTag);
	  assert(properTag);
	  genericValues = new char[*tag_value_size];
	  lastIndex = properTag->getTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<iBase_EntitySetHandle*>(genericValues));
	  if (lastIndex != entity_handles_size) {
	    delete[] genericValues;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;

      case iBase_BYTES:
	{
	  dataSizeInBytes = iNumValues;
	  *tag_value_size = entity_handles_size * dataSizeInBytes;
	  vEnsureArraySize<char>(tag_value, tag_value_allocated,
				 *tag_value_size);
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericValues = new char[*tag_value_size];
	  lastIndex = genericTag->getTagValueArr(entity_handles,
						 entity_handles_size,
						 genericValues);
	  if (lastIndex != entity_handles_size) {
	    delete[] genericValues;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;
      }
    for (int iByte = 0; iByte < *tag_value_size; iByte++) {
      (*tag_value)[iByte] = genericValues[iByte];
    }
    // This conditional can only fail for a bad data type.
    if (genericValues)
      delete[] genericValues;
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::getIntArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  int** tag_value, int* tag_value_allocated,
			  int* tag_value_size)

{
  int* piValues;
  // Templatized method bodies are found in the header file
  *tag_value_size = iGetSpecifArrData(entity_handles, entity_handles_size,
				      tag_handle, iBase_INTEGER, piValues);
  vEnsureArraySize<int>(tag_value, tag_value_allocated, *tag_value_size);
  for (int i = 0; i < *tag_value_size; i++)
    (*tag_value)[i] = piValues[i];
  delete[] piValues;
}

void
iMesh_Base::getDblArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  double** tag_value, int *tag_value_allocated,
			  int *tag_value_size)

{
  double* pdValues;
  // Templatized method bodies are found in the header file
  *tag_value_size = iGetSpecifArrData(entity_handles, entity_handles_size,
				      tag_handle, iBase_DOUBLE, pdValues);
  vEnsureArraySize<double>(tag_value, tag_value_allocated, *tag_value_size);
  for (int i = 0; i < *tag_value_size; i++)
    (*tag_value)[i] = pdValues[i];
  delete[] pdValues;
}

void
iMesh_Base::getEHArrData(const iBase_EntityHandle entity_handles[],
			 int entity_handles_size, iBase_TagHandle tag_handle,
			 iBase_EntityHandle** tag_value,
			 int* tag_value_allocated, int* tag_value_size)

{
  iBase_EntityHandle* tmp_values;
  // Templatized method bodies are found in the header file
  *tag_value_size = iGetSpecifArrData(entity_handles, entity_handles_size,
				      tag_handle, iBase_ENTITY_HANDLE,
				      tmp_values);
  vEnsureArraySize<iBase_EntityHandle>(tag_value, tag_value_allocated,
				       *tag_value_size);
  for (int i = 0; i < *tag_value_size; i++)
    (*tag_value)[i] = tmp_values[i];
  delete[] tmp_values;
}

void
iMesh_Base::getESHArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  iBase_EntitySetHandle** tag_value,
			  int* tag_value_allocated, int* tag_value_size)

{
  iBase_EntitySetHandle* tmp_values;
  // Templatized method bodies are found in the header file
  *tag_value_size = iGetSpecifArrData(entity_handles, entity_handles_size,
				      tag_handle, iBase_ENTITY_SET_HANDLE,
				      tmp_values);
  vEnsureArraySize<iBase_EntitySetHandle>(tag_value, tag_value_allocated,
					  *tag_value_size);
  for (int i = 0; i < *tag_value_size; i++)
    (*tag_value)[i] = tmp_values[i];
  delete[] tmp_values;
}

void
iMesh_Base::setArrData(const iBase_EntityHandle entity_handles[],
		       int entity_handles_size, iBase_TagHandle tag_handle,
		       const char* tag_values, int tag_values_size)

{
  try {
    for (int i = 0; i < entity_handles_size; i++) {
      iBase_EntityHandle EH_tmp = entity_handles[i];
      if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
	throw(iBase_INVALID_ENTITY_HANDLE);
    }
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    int iNumValues = ptbTheTag->iGetSize();
    int type = ptbTheTag->iGetTagType();
    switch (type)
      {

      case iBase_INTEGER:
	{
	  if (entity_handles_size * iNumValues * sizeof(int)
	      != static_cast<size_t>(tag_values_size)) {
	    throw(iBase_BAD_ARRAY_SIZE);
	  }
	  Tag<int>* intTag = dynamic_cast<Tag<int>*>(ptbTheTag);
	  assert(intTag);
	  intTag->setTagValueArr(entity_handles, entity_handles_size,
				 reinterpret_cast<const int*>(tag_values));
	}
	break;

      case iBase_DOUBLE:
	{
	  if (entity_handles_size * iNumValues * sizeof(double)
	      != static_cast<size_t>(tag_values_size)) {
	    throw(iBase_BAD_ARRAY_SIZE);
	  }
	  Tag<double>* doubleTag = dynamic_cast<Tag<double>*>(ptbTheTag);
	  assert(doubleTag);
	  doubleTag->setTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<const double*>(tag_values));
	}
	break;

      case iBase_ENTITY_HANDLE:
	{
	  if (entity_handles_size * iNumValues * sizeof(iBase_EntityHandle)
	      != static_cast<size_t>(tag_values_size)) {
	    throw(iBase_BAD_ARRAY_SIZE);
	  }
#ifndef NDEBUG
	  const iBase_EntityHandle* aehTemp =
	      reinterpret_cast<const iBase_EntityHandle*>(tag_values);
	  for (int i = 0; i < entity_handles_size * iNumValues; i++) {
	    iBase_EntityHandle EH_tmp = aehTemp[i];
	    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
	      throw(iBase_INVALID_ENTITY_HANDLE);
	  }
#endif
	  Tag<iBase_EntityHandle>* EHTag =
	      dynamic_cast<Tag<iBase_EntityHandle>*>(ptbTheTag);
	  assert(EHTag);
	  EHTag->setTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<const iBase_EntityHandle*>(tag_values));
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  if (entity_handles_size * iNumValues * sizeof(iBase_EntitySetHandle)
	      != static_cast<size_t>(tag_values_size)) {
	    throw(iBase_BAD_ARRAY_SIZE);
	  }
#ifndef NDEBUG
	  const iBase_EntitySetHandle* aeshTemp =
	      reinterpret_cast<const iBase_EntitySetHandle*>(tag_values);
	  for (int i = 0; i < entity_handles_size * iNumValues; i++) {
	    iBase_EntitySetHandle ESH_tmp = aeshTemp[i];
	    if (ESH_tmp == NULL
		|| !qHasSet(dynamic_cast<EntitySetBase*>(ESH_tmp)))
	      throw(iBase_INVALID_ENTITYSET_HANDLE);
	  }
#endif
	  Tag<iBase_EntitySetHandle>* ESHTag = dynamic_cast<Tag<
	      iBase_EntitySetHandle>*>(ptbTheTag);
	  assert(ESHTag);
	  ESHTag->setTagValueArr(
	      entity_handles, entity_handles_size,
	      reinterpret_cast<const iBase_EntitySetHandle*>(tag_values));
	}
	break;

      case iBase_BYTES:
	{
	  if (entity_handles_size * iNumValues != tag_values_size) {
	    throw(iBase_BAD_ARRAY_SIZE);
	  }
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericTag->setTagValueArr(entity_handles, entity_handles_size,
				     tag_values);
	}
	break;
      }
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::setIntArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  const int tag_values[], int tag_values_size)

{
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  // Templatized method bodies are found in the header file
  vSetSpecifArrData(entity_handles, entity_handles_size, tag_handle,
		    iBase_INTEGER, tag_values, tag_values_size);
}

void
iMesh_Base::setDblArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  const double tag_values[], int tag_values_size)

{
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  // Templatized method bodies are found in the header file
  vSetSpecifArrData(entity_handles, entity_handles_size, tag_handle,
		    iBase_DOUBLE, tag_values, tag_values_size);
}

void
iMesh_Base::setEHArrData(const iBase_EntityHandle entity_handles[],
			 int entity_handles_size, iBase_TagHandle tag_handle,
			 const iBase_EntityHandle tag_values[],
			 int tag_values_size)

{
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  for (int i = 0; i < tag_values_size; i++) {
    iBase_EntityHandle EH_tmp = tag_values[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  // Templatized method bodies are found in the header file
  vSetSpecifArrData(entity_handles, entity_handles_size, tag_handle,
		    iBase_ENTITY_HANDLE, tag_values, tag_values_size);
}

void
iMesh_Base::setESHArrData(const iBase_EntityHandle entity_handles[],
			  int entity_handles_size, iBase_TagHandle tag_handle,
			  const iBase_EntitySetHandle tag_values[],
			  int tag_values_size)

{
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  for (int i = 0; i < tag_values_size; i++) {
    iBase_EntitySetHandle ESH_tmp = tag_values[i];
    if (ESH_tmp == NULL || !qHasSet(dynamic_cast<EntitySetBase*>(ESH_tmp)))
      throw(iBase_INVALID_ENTITYSET_HANDLE);
  }
  // Templatized method bodies are found in the header file
  vSetSpecifArrData(entity_handles, entity_handles_size, tag_handle,
		    iBase_ENTITY_HANDLE, tag_values, tag_values_size);
}

void
iMesh_Base::rmvArrTag(const iBase_EntityHandle entity_handles[],
		      int entity_handles_size, iBase_TagHandle tag_handle)

{
  for (int i = 0; i < entity_handles_size; i++) {
    iBase_EntityHandle EH_tmp = entity_handles[i];
    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
      throw(iBase_INVALID_ENTITY_HANDLE);
  }
  TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
  // Tag member function works well enough
  ptbTheTag->iUntagEntityArr(entity_handles, entity_handles_size);
}

void
iMesh_Base::setEntSetData(iBase_EntitySetHandle entity_set,
			  iBase_TagHandle tag_handle, const char* tag_value,
			  int tag_size)

{
  try {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    int type = ptbTheTag->iGetTagType();
    int iNumValues = ptbTheTag->iGetSize();
    // Could have created a templatized method, but there are only four lines
    switch (type)
      {
      case iBase_INTEGER:
	{
	  if (static_cast<size_t>(tag_size) != iNumValues * sizeof(int)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  Tag<int>* properTag = dynamic_cast<Tag<int>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_set,
				 reinterpret_cast<const int*>(tag_value));
	  break;
	}

      case iBase_DOUBLE:
	{
	  if (static_cast<size_t>(tag_size) != iNumValues * sizeof(double)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  Tag<double>* properTag = dynamic_cast<Tag<double>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_set,
				 reinterpret_cast<const double*>(tag_value));
	  break;
	}

      case iBase_ENTITY_HANDLE:
	{
	  if (static_cast<size_t>(tag_size)
	      != iNumValues * sizeof(iBase_EntityHandle)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  const iBase_EntityHandle* aEHInEnts =
	      reinterpret_cast<iBase_EntityHandle const *>(tag_value);
	  for (int i = 0; i < iNumValues; i++) {
	    iBase_EntityHandle EH_tmp = aEHInEnts[i];
	    if (EH_tmp == NULL || !pM->isValidEntHandle(EH_tmp))
	      throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	  Tag<iBase_EntityHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntityHandle>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_set, aEHInEnts);
	  break;
	}

      case iBase_ENTITY_SET_HANDLE:
	{
	  if (static_cast<size_t>(tag_size)
	      != iNumValues * sizeof(iBase_EntitySetHandle)) {
	    throw(iBase_INVALID_TAG_HANDLE);
	  }
	  const iBase_EntitySetHandle* aESHInSets =
	      reinterpret_cast<iBase_EntitySetHandle const *>(tag_value);
	  for (int i = 0; i < iNumValues; i++) {
	    iBase_EntitySetHandle ESH_tmp = aESHInSets[i];
	    if (ESH_tmp == NULL
		|| !qHasSet(dynamic_cast<EntitySetBase*>(ESH_tmp)))
	      throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	  Tag<iBase_EntitySetHandle>* properTag = dynamic_cast<Tag<
	      iBase_EntitySetHandle>*>(ptbTheTag);
	  assert(properTag);
	  properTag->setTagValue(entity_set, aESHInSets);
	  break;
	}

      case iBase_BYTES:
	{
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericTag->setTagValue(entity_set, tag_value);
	}
	break;
      }
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::setEntSetIntData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle, int tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<int>(entity_set, tag_handle, iBase_INTEGER, tag_value);
}

void
iMesh_Base::setEntSetDblData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle, double tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<double>(entity_set, tag_handle, iBase_DOUBLE, tag_value);
}

void
iMesh_Base::setEntSetEHData(iBase_EntitySetHandle entity_set,
			    iBase_TagHandle tag_handle,
			    iBase_EntityHandle tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<iBase_EntityHandle>(entity_set, tag_handle,
				     iBase_ENTITY_HANDLE, tag_value);
}

void
iMesh_Base::setEntSetESHData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle,
			     iBase_EntitySetHandle tag_value)

{
  // Templatized method bodies are found in the header file
  vSetSpecifData<iBase_EntitySetHandle>(entity_set, tag_handle,
					iBase_ENTITY_SET_HANDLE, tag_value);
}

void
iMesh_Base::getEntSetData(iBase_EntitySetHandle entity_set,
			  iBase_TagHandle tag_handle, char** tag_value,
			  int* tag_value_allocated, int* tag_value_size)

{
  char* genericData = NULL;
  try {
    TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
    int type = ptbTheTag->iGetTagType();
    // Can use the same templatized method as for getData, since its
    // operation does not depend on the key type
    switch (type)
      {
      case iBase_INTEGER:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<int>(entity_set, ptbTheTag,
						 genericData);
	}
	break;

      case iBase_DOUBLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<double>(entity_set, ptbTheTag,
						    genericData);
	}
	break;

      case iBase_ENTITY_HANDLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<iBase_EntityHandle>(entity_set,
								ptbTheTag,
								genericData);
	}
	break;

      case iBase_ENTITY_SET_HANDLE:
	{
	  // Templatized method bodies are found in the header file
	  *tag_value_size = iGetGenericData<iBase_EntitySetHandle>(entity_set,
								   ptbTheTag,
								   genericData);
	}
	break;

      case iBase_BYTES:
	{
	  *tag_value_size = ptbTheTag->iGetSize();
	  Tag<char>* genericTag = dynamic_cast<Tag<char>*>(ptbTheTag);
	  assert(genericTag);
	  genericData = new char[*tag_value_size];
	  if (!genericTag->getTagValue(entity_set, genericData)) {
	    delete[] genericData;
	    throw(iBase_INVALID_ENTITY_HANDLE);
	  }
	}
	break;
      }
    vEnsureArraySize<char>(tag_value, tag_value_allocated, *tag_value_size);
    for (int i = 0; i < *tag_value_size; i++)
      (*tag_value)[i] = genericData[i];
    if (genericData)
      delete[] genericData;
  }
  catch (std::bad_alloc&) {
    if (genericData)
      delete[] genericData;
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
  catch (...) {
    if (genericData)
      delete[] genericData;
    throw;
  }
}

int
iMesh_Base::getEntSetIntData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<int>(entity_set, tag_handle, iBase_INTEGER);
}

double
iMesh_Base::getEntSetDblData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<double>(entity_set, tag_handle, iBase_DOUBLE);
}

iBase_EntityHandle
iMesh_Base::getEntSetEHData(iBase_EntitySetHandle entity_set,
			    iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<iBase_EntityHandle>(entity_set, tag_handle,
					    iBase_ENTITY_HANDLE);
}

iBase_EntitySetHandle
iMesh_Base::getEntSetESHData(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle tag_handle)

{
  // Templatized method bodies are found in the header file
  return tGetSpecifData<iBase_EntitySetHandle>(entity_set, tag_handle,
					       iBase_ENTITY_SET_HANDLE);
}

void
iMesh_Base::getAllEntSetTags(iBase_EntitySetHandle entity_set,
			     iBase_TagHandle** tag_handles,
			     int *tag_handles_allocated, int *tag_handles_size)

{
  try {
    std::set<TagBase*>::iterator iter = spTags.begin(), iterEnd = spTags.end();
    *tag_handles_size = 0;
    // Create a vector to store the tags temporarily, since it is
    // impossible to know how many there will be
    std::vector<iBase_TagHandle> tagVect;
    for (; iter != iterEnd; ++iter)
      if ((*iter)->qIsTagged(entity_set))
	tagVect.push_back(*iter);
    *tag_handles_size = tagVect.size();
    vEnsureArraySize<iBase_TagHandle>(tag_handles, tag_handles_allocated,
				      *tag_handles_size);
    for (int i = 0; i < *tag_handles_size; i++)
      (*tag_handles)[i] = tagVect[i];
  }
  catch (std::bad_alloc&) {
    throw(iBase_MEMORY_ALLOCATION_FAILED);
  }
}

void
iMesh_Base::rmvEntSetTag(iBase_EntitySetHandle entity_set,
			 iBase_TagHandle tag_handle)

{
  TagBase* ptbTheTag = static_cast<TagBase*>(tag_handle);
  // Tag member function works well enough
  ptbTheTag->qUntagEntity(entity_set);
}

void
iMesh_Base::getEntArr2ndAdj(iBase_EntityHandle const* entity_handles,
			    const int entity_handles_size,
			    const int bridge_entity_type,
			    const int requested_entity_type,
			    iBase_EntityHandle** adj_entity_handles,
			    int* adj_entity_handles_allocated,
			    int* adj_entity_handles_size, int** offset,
			    int* offset_allocated, int* offset_size, int* err)
{
  m_err = *err = iBase_SUCCESS;
  TRY_ARRAY_SIZE(offset, offset_allocated, entity_handles_size + 1, int);
  *offset_size = entity_handles_size + 1;

  // It's a little crazy, but the only visibly feasible way to do this
  // is to compute the second adjacencies once to find out how big they
  // are, then size the output array, then compute them again for money.
  int ii;
  (*offset)[0] = 0;
  for (ii = 0; ii < entity_handles_size; ii++) {
    iBase_EntityHandle* temp_handles = NULL;
    int temp_handles_allocated = 0, temp_handles_size;
    getEnt2ndAdj(entity_handles[ii], bridge_entity_type, requested_entity_type,
		 &temp_handles, &temp_handles_allocated, &temp_handles_size,
		 err);
    if (*err != iBase_SUCCESS)
      return;
    (*offset)[ii + 1] = (*offset)[ii] + temp_handles_size;
    if (temp_handles)
      free(temp_handles);
  }

  *adj_entity_handles_size = (*offset)[entity_handles_size];
  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		 *adj_entity_handles_size, iBase_EntityHandle);

  for (ii = 0; ii < entity_handles_size; ii++) {
    int temp_handles_allocated = (*offset)[ii + 1] - (*offset)[ii],
	temp_handles_size;
    iBase_EntityHandle *temp_handles = (*adj_entity_handles) + (*offset)[ii];
    getEnt2ndAdj(entity_handles[ii], bridge_entity_type, requested_entity_type,
		 &temp_handles, &temp_handles_allocated, &temp_handles_size,
		 err);
    if (*err != iBase_SUCCESS)
      return;
    assert((*offset)[ii + 1] - (*offset)[ii] == temp_handles_size);
  }
}

void
iMesh_Base::writeTags(const char tagFileName[],
		      const iBase_EntityHandle allHandles[],
		      const int allHandles_size,
		      const iBase_EntitySetHandle allSetHandles[],
		      const int allSetHandles_size)
{
  int numTags = spTags.size();

  FILE* pFTagFile = fopen(tagFileName, "w");
  if (numTags == 0) {
    fclose(pFTagFile);
    unlink(tagFileName);
    return;
  }

  TagBase* *allTags = new TagBase*[numTags];
  std::copy(spTags.begin(), spTags.end(), allTags);

  // Write basic data for each tag.
  int *allTypes = new int[numTags];
  int *allSizes = new int[numTags];

  fprintf(pFTagFile, "Total number of tags: %d\n", numTags);
  for (int tag = 0; tag < numTags; tag++) {
    allSizes[tag] = getTagSizeValues(allTags[tag]);
    allTypes[tag] = getTagType(allTags[tag]);
    char namePtr[1024];
    getTagName(allTags[tag], namePtr, 1024);
    fprintf(pFTagFile, "Tag #%d: Type: %d Size: %d Name: %s \n", tag,
	    allTypes[tag], allSizes[tag], namePtr);
  }
  fprintf(pFTagFile, "\n");

  iBase_TagHandle *myTags = new iBase_TagHandle[numTags];
  int myTags_allocated = numTags, myTags_size = 0;

  // Now write data for each entity.
  for (int ent = 0; ent < allHandles_size; ent++) {
    iBase_EntityHandle entHandle = allHandles[ent];
    getAllTags(entHandle, &myTags, &myTags_allocated, &myTags_size);
    if (myTags_size != 0) {
      fprintf(pFTagFile, "Entity %d has %d tags.\n", ent, myTags_size);
      for (int ii = 0; ii < myTags_size; ii++) {
	iBase_TagHandle tagHandle = myTags[ii];
	int tagIndex = findArrayIndex(static_cast<TagBase*>(tagHandle), allTags,
				      allTags + numTags);

	fprintf(pFTagFile, "  Tag %d. Value(s): ", tagIndex);
	switch (allTypes[tagIndex])
	  {
	  case iBase_INTEGER:
	    if (allSizes[tagIndex] == 1) {
	      fprintf(pFTagFile, "%d\n", getIntData(entHandle, tagHandle));
	    }
	    else {
	      int *myData = NULL, myData_allocated = 0, myData_size;
	      getIntArrData(&entHandle, 1, tagHandle, &myData,
			    &myData_allocated, &myData_size);
	      assert(myData_size == allSizes[tagIndex]);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		fprintf(pFTagFile, "%d ", myData[iDat]);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  case iBase_DOUBLE:
	    if (allSizes[tagIndex] == 1) {
	      fprintf(pFTagFile, "%.12g\n", getDblData(entHandle, tagHandle));
	    }
	    else {
	      double *myData = NULL;
	      int myData_allocated = 0, myData_size;
	      getDblArrData(&entHandle, 1, tagHandle, &myData,
			    &myData_allocated, &myData_size);
	      assert(myData_size == allSizes[tagIndex]);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		fprintf(pFTagFile, "%.12g ", myData[iDat]);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  case iBase_ENTITY_HANDLE:
	    if (allSizes[tagIndex] == 1) {
	      iBase_EntityHandle dataHandle = getEHData(entHandle, tagHandle);
	      int entIndex = findArrayIndex(dataHandle, allHandles,
					    allHandles + allHandles_size);
	      fprintf(pFTagFile, "%d\n", entIndex);
	    }
	    else {
	      iBase_EntityHandle *myData = NULL;
	      int myData_allocated = 0, myData_size;
	      getEHArrData(&entHandle, 1, tagHandle, &myData, &myData_allocated,
			   &myData_size);
	      assert(myData_size == allSizes[tagIndex]);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		iBase_EntityHandle dataHandle = myData[iDat];
		int entIndex = findArrayIndex(dataHandle, allHandles,
					      allHandles + allHandles_size);
		fprintf(pFTagFile, "%d ", entIndex);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  case iBase_ENTITY_SET_HANDLE:
	    if (allSizes[tagIndex] == 1) {
	      iBase_EntitySetHandle dataHandle = getESHData(entHandle,
							    tagHandle);
	      int setIndex = findArrayIndex(dataHandle, allSetHandles,
					    allSetHandles + allSetHandles_size);
	      fprintf(pFTagFile, "%d\n", setIndex);
	    }
	    else {
	      iBase_EntitySetHandle *myData = NULL;
	      int myData_allocated = 0, myData_size;
	      getESHArrData(&entHandle, 1, tagHandle, &myData,
			    &myData_allocated, &myData_size);
	      assert(myData_size == allSizes[tagIndex]);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		iBase_EntitySetHandle dataHandle = myData[iDat];
		int setIndex = findArrayIndex(
		    dataHandle, allSetHandles,
		    allSetHandles + allSetHandles_size);
		fprintf(pFTagFile, "%d ", setIndex);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  case iBase_BYTES:
	    {
	      assert(allSizes[ii] >= 1);
	      char *myData = NULL;
	      int myData_allocated = 0, myData_size;
	      getData(entHandle, tagHandle, &myData, &myData_allocated,
		      &myData_size);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		fprintf(pFTagFile, "%c", myData[iDat]);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  }
      } // Done with all tags for this entity
      fprintf(pFTagFile, "\n");
    } // Done with if (tags exist for this entity)
  } // Done checking all entities

    // Now repeat for each set.
  for (int set = 0; set < allSetHandles_size; set++) {
    iBase_EntitySetHandle setHandle = allSetHandles[set];
    getAllEntSetTags(setHandle, &myTags, &myTags_allocated, &myTags_size);
    if (myTags_size != 0) {
      fprintf(pFTagFile, "Set %d has %d tags.\n", set, myTags_size);
      for (int ii = 0; ii < myTags_size; ii++) {
	iBase_TagHandle tagHandle = myTags[ii];
	int tagIndex = findArrayIndex(static_cast<TagBase*>(tagHandle), allTags,
				      allTags + numTags);

	fprintf(pFTagFile, "  Tag %d. Value(s): ", tagIndex);
	assert(allSizes[tagIndex] == 1 || allTypes[tagIndex] == iBase_BYTES);
	switch (allTypes[tagIndex])
	  {
	  case iBase_INTEGER:
	    fprintf(pFTagFile, "%d\n", getEntSetIntData(setHandle, tagHandle));
	    break;
	  case iBase_DOUBLE:
	    fprintf(pFTagFile, "%.12g\n",
		    getEntSetDblData(setHandle, tagHandle));
	    break;
	  case iBase_ENTITY_HANDLE:
	    {
	      iBase_EntityHandle dataHandle = getEntSetEHData(setHandle,
							      tagHandle);
	      int entIndex = findArrayIndex(dataHandle, allHandles,
					    allHandles + allHandles_size);
	      fprintf(pFTagFile, "%d\n", entIndex);
	    }
	    break;
	  case iBase_ENTITY_SET_HANDLE:
	    {
	      iBase_EntitySetHandle dataHandle = getEntSetESHData(setHandle,
								  tagHandle);
	      int setIndex = findArrayIndex(dataHandle, allSetHandles,
					    allSetHandles + allSetHandles_size);
	      fprintf(pFTagFile, "%d\n", setIndex);
	    }
	    break;
	  case iBase_BYTES:
	    {
	      assert(allSizes[ii] >= 1);
	      char *myData = NULL;
	      int myData_allocated = 0, myData_size;
	      getEntSetData(setHandle, tagHandle, &myData, &myData_allocated,
			    &myData_size);
	      for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
		fprintf(pFTagFile, "%c", myData[iDat]);
	      }
	      fprintf(pFTagFile, "\n");
	      if (myData)
		free(myData);
	    }
	    break;
	  }
      } // Done with all tags for this set
      fprintf(pFTagFile, "\n");
    } // Done with if (tags exist for this set)
  } // Done checking all sets

    // Finally the root set.
  getAllEntSetTags(GR_rootSet, &myTags, &myTags_allocated, &myTags_size);
  if (myTags_size != 0) {
    fprintf(pFTagFile, "Root set has %d tags.\n", myTags_size);
    for (int ii = 0; ii < myTags_size; ii++) {
      iBase_TagHandle tagHandle = myTags[ii];
      int tagIndex = findArrayIndex(static_cast<TagBase*>(tagHandle), allTags,
				    allTags + numTags);

      fprintf(pFTagFile, "  Tag %d. Value(s): ", tagIndex);
      assert(allSizes[tagIndex] == 1 || allTypes[tagIndex] == iBase_BYTES);
      switch (allTypes[tagIndex])
	{
	case iBase_INTEGER:
	  fprintf(pFTagFile, "%d\n", getEntSetIntData(GR_rootSet, tagHandle));
	  break;
	case iBase_DOUBLE:
	  fprintf(pFTagFile, "%.12g\n",
		  getEntSetDblData(GR_rootSet, tagHandle));
	  break;
	case iBase_ENTITY_HANDLE:
	  {
	    iBase_EntityHandle dataHandle = getEntSetEHData(GR_rootSet,
							    tagHandle);
	    int entIndex = findArrayIndex(dataHandle, allHandles,
					  allHandles + allHandles_size);
	    fprintf(pFTagFile, "%d\n", entIndex);
	  }
	  break;
	case iBase_ENTITY_SET_HANDLE:
	  {
	    iBase_EntitySetHandle dataHandle = getEntSetESHData(GR_rootSet,
								tagHandle);
	    int setIndex = findArrayIndex(dataHandle, allSetHandles,
					  allSetHandles + allSetHandles_size);
	    fprintf(pFTagFile, "%d\n", setIndex);
	  }
	  break;
	case iBase_BYTES:
	  {
	    assert(allSizes[ii] >= 1);
	    char *myData = NULL;
	    int myData_allocated = 0, myData_size;
	    getEntSetData(GR_rootSet, tagHandle, &myData, &myData_allocated,
			  &myData_size);
	    for (int iDat = 0; iDat < allSizes[tagIndex]; iDat++) {
	      fprintf(pFTagFile, "%c", myData[iDat]);
	    }
	    fprintf(pFTagFile, "\n");
	    if (myData)
	      free(myData);
	  }
	  break;
	}
    } // Done with all tags for this set
    fprintf(pFTagFile, "\n");
  } // Done with if (tags exist for this set)

  fclose(pFTagFile);

  delete[] allSizes;
  delete[] allTypes;
  delete[] allTags;
  delete[] myTags;
}
