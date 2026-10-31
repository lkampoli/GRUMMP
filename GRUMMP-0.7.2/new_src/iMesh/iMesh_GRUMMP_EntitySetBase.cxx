#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_misc.hh"

using std::set;

#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"

namespace ITAPS_GRUMMP
{

  EntitySetBase::~EntitySetBase()
  {
    set<EntitySetBase*>::iterator iter, iterEnd;

    // Sever ties with children
    iterEnd = spESChildren.end();
    for (iter = spESChildren.begin(); iter != iterEnd; ++iter)
      (*iter)->spESParents.erase(this);

    // Sever ties with parents
    iterEnd = spESParents.end();
    for (iter = spESParents.begin(); iter != iterEnd; ++iter)
      (*iter)->spESChildren.erase(this);

    // Sever ties with owners
    iterEnd = spESOwners.end();
    for (iter = spESOwners.begin(); iter != iterEnd; ++iter)
      (*iter)->spESMembers.erase(this);

    // Sever ties with members
    iterEnd = spESMembers.end();
    for (iter = spESMembers.begin(); iter != iterEnd; ++iter)
      (*iter)->spESOwners.erase(this);
  }

  int
  EntitySetBase::iGetNumEntSets(const int iNumHops) const

  {
    set<EntitySetBase*> spESRes;
    vGetEntitySets(spESRes, iNumHops);
    return spESRes.size();
  }

  int
  EntitySetBase::iGetNumChildren(const int iNumHops) const

  {
    set<EntitySetBase*> spESRes;
    vGetChildren(spESRes, iNumHops);
    return spESRes.size();
  }

  int
  EntitySetBase::iGetNumParents(const int iNumHops) const

  {
    set<EntitySetBase*> spESRes;
    vGetParents(spESRes, iNumHops);
    return spESRes.size();
  }

  void
  EntitySetBase::vGetEntitySets(set<EntitySetBase*>& spESRes,
				const int iNumHops) const

  {
    set<EntitySetBase*> spESToCheck(spESMembers);
    spESRes.insert(spESMembers.begin(), spESMembers.end());
    int iSize = spESRes.size(), iOldSize = 0;
    int iHops = 0;
    while (iHops != iNumHops && iSize > iOldSize) {
      iHops++;
      iOldSize = iSize;
      set<EntitySetBase*> spESNew;
      set<EntitySetBase*>::iterator iter = spESToCheck.begin(), iterEnd =
	  spESToCheck.end();
      for (; iter != iterEnd; ++iter) {
	spESNew.insert((*iter)->spESMembers.begin(),
		       (*iter)->spESMembers.end());
      }
      // Only need to check those that aren't already in the result.
      // 'insert' can be used with a range from a different set, but
      // 'erase' cannot, so the set must be iterated through
      iterEnd = spESRes.end();
      for (iter = spESRes.begin(); iter != iterEnd; ++iter)
	spESNew.erase(*iter);
      // Now put all of these into the result list:
      spESRes.insert(spESNew.begin(), spESNew.end());
      spESToCheck = spESNew;
      iSize = spESRes.size();
    }
  }

  void
  EntitySetBase::vGetEntitySets(iBase_EntitySetHandle** ent_set_handles,
				int *ent_set_handles_allocated,
				int *ent_set_handles_size,
				const int iNumHops) const

  {
    set<EntitySetBase*> spESResult;
    vGetEntitySets(spESResult, iNumHops);

    *ent_set_handles_size = spESResult.size();
    vEnsureArraySize<iBase_EntitySetHandle>(ent_set_handles,
					    ent_set_handles_allocated,
					    *ent_set_handles_size);

    set<EntitySetBase*>::iterator iter = spESResult.begin(), iterEnd =
	spESResult.end();
    int iES;
    for (iES = 0; iter != iterEnd; ++iter, iES++) {
      (*ent_set_handles)[iES] = *iter;
    }
    assert(iES == *ent_set_handles_size);
  }

  void
  EntitySetBase::vAddEntitySet(void* ent_set_handle)

  {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(ent_set_handle);
    spESMembers.insert(pESB);
    pESB->spESOwners.insert(this);
  }

  void
  EntitySetBase::vAddEntitySet(EntitySetBase* pESB)

  {
    spESMembers.insert(pESB);
    pESB->spESOwners.insert(this);
  }

  void
  EntitySetBase::vRemoveEntitySet(void* ent_set_handle)

  {
    EntitySetBase *pESB = static_cast<EntitySetBase*>(ent_set_handle);
    spESMembers.erase(pESB);
    pESB->spESOwners.erase(this);
  }

  void
  EntitySetBase::vGetChildren(set<EntitySetBase*>& spESResult,
			      const int iHops) const

  {
    spESResult.insert(spESChildren.begin(), spESChildren.end());
    set<EntitySetBase*> spESToCheck(spESChildren);
    int iH = 0, iSize = spESResult.size(), iOldSize = 0;
    while (iH != iHops && iSize > iOldSize) {
      iOldSize = iSize;
      set<EntitySetBase*> spESNew;
      set<EntitySetBase*>::iterator iter = spESToCheck.begin(), iterEnd =
	  spESToCheck.end();
      for (; iter != iterEnd; ++iter) {
	spESNew.insert((*iter)->spESChildren.begin(),
		       (*iter)->spESChildren.end());
      }
      // Only need to check those that aren't already in the result.
      iterEnd = spESResult.end();
      for (iter = spESResult.begin(); iter != iterEnd; ++iter)
	spESNew.erase(*iter);
      // Now put all of these into the result list:
      spESResult.insert(spESNew.begin(), spESNew.end());
      spESToCheck = spESNew;
      iSize = spESResult.size();
      iH++;
    }
  }

  void
  EntitySetBase::vGetChildren(const int iHops,
			      iBase_EntitySetHandle** ent_set_handles,
			      int *ent_set_handles_allocated,
			      int *ent_set_handles_size) const

  {
    set<EntitySetBase*> spESResult;
    vGetChildren(spESResult, iHops);

    *ent_set_handles_size = spESResult.size();
    vEnsureArraySize<iBase_EntitySetHandle>(ent_set_handles,
					    ent_set_handles_allocated,
					    *ent_set_handles_size);

    set<EntitySetBase*>::iterator iter = spESResult.begin();
    int iES;
    for (iES = 0; iter != spESResult.end(); ++iter, iES++) {
      (*ent_set_handles)[iES] = *iter;
    }
    assert(iES == *ent_set_handles_size);
  }

  void
  EntitySetBase::vGetParents(set<EntitySetBase*>& spESResult,
			     const int iHops) const

  {
    spESResult.insert(spESParents.begin(), spESParents.end());
    set<EntitySetBase*> spESToCheck(spESParents);
    int iH = 0, iSize = spESResult.size(), iOldSize = 0;
    while (iH != iHops && iSize > iOldSize) {
      iOldSize = iSize;
      set<EntitySetBase*> spESNew;
      set<EntitySetBase*>::iterator iter = spESToCheck.begin(), iterEnd =
	  spESToCheck.end();
      for (; iter != iterEnd; ++iter) {
	spESNew.insert((*iter)->spESParents.begin(),
		       (*iter)->spESParents.end());
      }
      // Only need to check those that aren't already in the result.
      iterEnd = spESResult.end();
      for (iter = spESResult.begin(); iter != iterEnd; ++iter)
	spESNew.erase(*iter);
      // Now put all of these into the result list:
      spESResult.insert(spESNew.begin(), spESNew.end());
      spESToCheck = spESNew;
      iSize = spESResult.size();
      iH++;
    }
  }

  void
  EntitySetBase::vGetParents(const int iHops,
			     iBase_EntitySetHandle** ent_set_handles,
			     int *ent_set_handles_allocated,
			     int *ent_set_handles_size) const

  {
    set<EntitySetBase*> spESResult;
    vGetParents(spESResult, iHops);

    *ent_set_handles_size = spESResult.size();
    vEnsureArraySize<iBase_EntitySetHandle>(ent_set_handles,
					    ent_set_handles_allocated,
					    *ent_set_handles_size);

    set<EntitySetBase*>::iterator iter = spESResult.begin();
    int iES;
    for (iES = 0; iter != spESResult.end(); ++iter, iES++) {
      (*ent_set_handles)[iES] = *iter;
    }
    assert(iES == *ent_set_handles_size);
  }

// Other methods named vUnite, vSubtract, and vIntersect are found in 
// 'iMesh_GRUMMP_SetBoolOps.cc'

  void
  vUnite(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	 EntitySetBase*& pESBRes)
  {
    bool qOrdered1 = pESB1->qIsOrdered();
    bool qOrdered2 = pESB2->qIsOrdered();
    if (qOrdered1) {
      const EntityList *pEL1 = dynamic_cast<const EntityList*>(pESB1);
      assert(pEL1);
      if (qOrdered2) {
	// Both are lists.
	const EntityList *pEL2 = dynamic_cast<const EntityList*>(pESB2);
	// Verify that the casts worked properly.
	assert(pEL2);
	EntityList *pELRes = new EntityList;
	vUnite(*pEL1, *pEL2, *pELRes);
	pESBRes = pELRes;
      }
      else {
	// Second is not a list, so the result will be a set.
	EntitySet ES1(*pEL1);
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vUnite(ES1, *pES2, *pESRes);
	pESBRes = pESRes;
      }
    }
    else {
      const EntitySet *pES1 = dynamic_cast<const EntitySet*>(pESB1);
      assert(pES1);
      if (qOrdered2) {
	// Second is a list; convert it to a set before operating.
	const EntityList* pEL2 = dynamic_cast<const EntityList*>(pESB2);
	assert(pEL2);
	EntitySet ES2(*pEL2);
	EntitySet *pESRes = new EntitySet;
	vUnite(*pES1, ES2, *pESRes);
	pESBRes = pESRes;
      }
      else {
	// Both are sets.
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vUnite(*pES1, *pES2, *pESRes);
	pESBRes = pESRes;
      }
    }
  }

  void
  vSubtract(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	    EntitySetBase*& pESBRes)
  {
    bool qOrdered1 = pESB1->qIsOrdered();
    bool qOrdered2 = pESB2->qIsOrdered();
    if (qOrdered1) {
      const EntityList *pEL1 = dynamic_cast<const EntityList*>(pESB1);
      assert(pEL1);
      if (qOrdered2) {
	// Both are lists.
	const EntityList *pEL2 = dynamic_cast<const EntityList*>(pESB2);
	// Verify that the casts worked properly.
	assert(pEL2);
	EntityList *pELRes = new EntityList;
	vSubtract(*pEL1, *pEL2, *pELRes);
	pESBRes = pELRes;
      }
      else {
	// Second is not a list, so the result will be a set.
	EntitySet ES1(*pEL1);
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vSubtract(ES1, *pES2, *pESRes);
	pESBRes = pESRes;
      }
    }
    else {
      const EntitySet *pES1 = dynamic_cast<const EntitySet*>(pESB1);
      assert(pES1);
      if (qOrdered2) {
	// Second is a list; convert it to a set before operating.
	const EntityList* pEL2 = dynamic_cast<const EntityList*>(pESB2);
	assert(pEL2);
	EntitySet ES2(*pEL2);
	EntitySet *pESRes = new EntitySet;
	vSubtract(*pES1, ES2, *pESRes);
	pESBRes = pESRes;
      }
      else {
	// Both are sets.
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vSubtract(*pES1, *pES2, *pESRes);
	pESBRes = pESRes;
      }
    }
  }

  void
  vIntersect(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	     EntitySetBase*& pESBRes)
  {
    bool qOrdered1 = pESB1->qIsOrdered();
    bool qOrdered2 = pESB2->qIsOrdered();
    if (qOrdered1) {
      const EntityList *pEL1 = dynamic_cast<const EntityList*>(pESB1);
      assert(pEL1);
      if (qOrdered2) {
	// Both are lists.
	const EntityList *pEL2 = dynamic_cast<const EntityList*>(pESB2);
	// Verify that the casts worked properly.
	assert(pEL2);
	EntityList *pELRes = new EntityList;
	vIntersect(*pEL1, *pEL2, *pELRes);
	pESBRes = pELRes;
      }
      else {
	// Second is not a list, so the result will be a set.
	EntitySet ES1(*pEL1);
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vIntersect(ES1, *pES2, *pESRes);
	pESBRes = pESRes;
	std::set<EntitySetBase*>::iterator iter = ES1.spESMembers.begin(),
	    iterEnd = ES1.spESMembers.end();
	for (; iter != iterEnd; ++iter)
	  (*iter)->vRemoveOwner(&ES1);
      }
    }
    else {
      const EntitySet *pES1 = dynamic_cast<const EntitySet*>(pESB1);
      assert(pES1);
      if (qOrdered2) {
	// Second is a list; convert it to a set before operating.
	const EntityList* pEL2 = dynamic_cast<const EntityList*>(pESB2);
	assert(pEL2);
	EntitySet ES2(*pEL2);
	EntitySet *pESRes = new EntitySet;
	vIntersect(*pES1, ES2, *pESRes);
	pESBRes = pESRes;
	std::set<EntitySetBase*>::iterator iter = ES2.spESMembers.begin(),
	    iterEnd = ES2.spESMembers.end();
	for (; iter != iterEnd; ++iter)
	  (*iter)->vRemoveOwner(&ES2);
      }
      else {
	// Both are sets.
	const EntitySet *pES2 = dynamic_cast<const EntitySet*>(pESB2);
	assert(pES2);
	EntitySet *pESRes = new EntitySet;
	vIntersect(*pES1, *pES2, *pESRes);
	pESBRes = pESRes;
      }
    }
  }
} // end namespace ITAPS_GRUMMP
