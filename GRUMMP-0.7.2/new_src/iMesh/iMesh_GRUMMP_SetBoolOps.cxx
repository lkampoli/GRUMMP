#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_misc.hh"
#include "GR_Entity.h"
#include "GR_Face.h"

#include <algorithm>

template<class _Key>
  inline static void
  vSubtract(const std::set<_Key>& sA, const std::set<_Key>& sB,
	    std::set<_Key>& sResult)
  {
    sResult = sA;
    typename std::set<_Key>::iterator iter = sB.begin(), iterEnd = sB.end();
    while (iter != iterEnd) {
      sResult.erase(*iter);
      ++iter;
    }
  }

template<class _Key>
  inline static void
  vIntersect(const std::set<_Key>& sA, const std::set<_Key>& sB,
	     std::set<_Key>& sResult)
  {
    int iMaxSize = sA.size();
    _Key aKTmp[iMaxSize];
    _Key *pKTmp = ::std::set_intersection(sA.begin(), sA.end(), sB.begin(),
					  sB.end(), aKTmp);
    sResult.clear();
    sResult.insert(aKTmp, pKTmp);
  }

template<class _Key>
  inline static void
  vUnite(const std::set<_Key>& sA, const std::set<_Key>& sB,
	 std::set<_Key>& sResult)
  {
    sResult = sA;
    sResult.insert(sB.begin(), sB.end());
  }

namespace ITAPS_GRUMMP
{
  void
  vSubtract(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESResult)
  {
    // The only way to do this, alas, is to subtract one set at a time.
    ::vSubtract(ES1.sVR, ES2.sVR, ESResult.sVR);
    ::vSubtract(ES1.spEdgeF, ES2.spEdgeF, ESResult.spEdgeF);
    ::vSubtract(ES1.spTriF, ES2.spTriF, ESResult.spTriF);
    ::vSubtract(ES1.spQuadF, ES2.spQuadF, ESResult.spQuadF);
    ::vSubtract(ES1.spTriC, ES2.spTriC, ESResult.spTriC);
    ::vSubtract(ES1.spQuadC, ES2.spQuadC, ESResult.spQuadC);
    ::vSubtract(ES1.spTetC, ES2.spTetC, ESResult.spTetC);
    ::vSubtract(ES1.spPyrC, ES2.spPyrC, ESResult.spPyrC);
    ::vSubtract(ES1.spPrismC, ES2.spPrismC, ESResult.spPrismC);
    ::vSubtract(ES1.spHexC, ES2.spHexC, ESResult.spHexC);

    // Don't forget to update the entity set data, too.  Compute the set
    // difference of the sets that are members of this one, and add
    // member-owner links for them.
    std::set<EntitySetBase*> spESMemberDiff;
    ::vSubtract(ES1.spESMembers, ES2.spESMembers, spESMemberDiff);
    std::set<EntitySetBase*>::iterator iter;
    for (iter = spESMemberDiff.begin(); iter != spESMemberDiff.end(); ++iter) {
      ESResult.vAddEntitySet(*iter);
    }

    // Build the auxiliary vertex list from scratch.
    ESResult.vRebuildAuxiliaryVertList();
  }

  void
  vIntersect(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESResult)
  {
    // The only way to do this, alas, is to intersect one set at a time.
    ::vIntersect(ES1.sVR, ES2.sVR, ESResult.sVR);
    ::vIntersect(ES1.spEdgeF, ES2.spEdgeF, ESResult.spEdgeF);
    ::vIntersect(ES1.spTriF, ES2.spTriF, ESResult.spTriF);
    ::vIntersect(ES1.spQuadF, ES2.spQuadF, ESResult.spQuadF);
    ::vIntersect(ES1.spTriC, ES2.spTriC, ESResult.spTriC);
    ::vIntersect(ES1.spQuadC, ES2.spQuadC, ESResult.spQuadC);
    ::vIntersect(ES1.spTetC, ES2.spTetC, ESResult.spTetC);
    ::vIntersect(ES1.spPyrC, ES2.spPyrC, ESResult.spPyrC);
    ::vIntersect(ES1.spPrismC, ES2.spPrismC, ESResult.spPrismC);
    ::vIntersect(ES1.spHexC, ES2.spHexC, ESResult.spHexC);

    // Don't forget to update the entity set data, too.  Compute the set
    // intersection of the sets that are members of this one, and add
    // member-owner links for them.
    std::set<EntitySetBase*> spESMemberInt;
    ::vIntersect(ES1.spESMembers, ES2.spESMembers, spESMemberInt);
    std::set<EntitySetBase*>::iterator iter;
    for (iter = spESMemberInt.begin(); iter != spESMemberInt.end(); ++iter) {
      ESResult.vAddEntitySet(*iter);
    }

    // Build the auxiliary vertex list from scratch.
    ESResult.vRebuildAuxiliaryVertList();
  }

  void
  vUnite(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESResult)
  {
    // The only way to do this, alas, is to unite one set at a time.
    ::vUnite(ES1.sVR, ES2.sVR, ESResult.sVR);
    ::vUnite(ES1.spEdgeF, ES2.spEdgeF, ESResult.spEdgeF);
    ::vUnite(ES1.spTriF, ES2.spTriF, ESResult.spTriF);
    ::vUnite(ES1.spQuadF, ES2.spQuadF, ESResult.spQuadF);
    ::vUnite(ES1.spTriC, ES2.spTriC, ESResult.spTriC);
    ::vUnite(ES1.spQuadC, ES2.spQuadC, ESResult.spQuadC);
    ::vUnite(ES1.spTetC, ES2.spTetC, ESResult.spTetC);
    ::vUnite(ES1.spPyrC, ES2.spPyrC, ESResult.spPyrC);
    ::vUnite(ES1.spPrismC, ES2.spPrismC, ESResult.spPrismC);
    ::vUnite(ES1.spHexC, ES2.spHexC, ESResult.spHexC);

    // Don't forget to update the entity set data, too.  Compute the set
    // union of the sets that are members of this one, and add
    // member-owner links for them.
    std::set<EntitySetBase*> spESMemberUnion;
    ::vUnite(ES1.spESMembers, ES2.spESMembers, spESMemberUnion);
    std::set<EntitySetBase*>::iterator iter;
    for (iter = spESMemberUnion.begin(); iter != spESMemberUnion.end();
	++iter) {
      ESResult.vAddEntitySet(*iter);
    }

    // Build the auxiliary vertex list from scratch.
    ESResult.vRebuildAuxiliaryVertList();
  }

  void
  vSubtract(const EntityList& EL1, const EntityList& EL2, EntityList& ELResult)
  {
    ELResult.lpEntities = EL1.lpEntities;

    ELResult.iVerts = EL1.iVerts;
    ELResult.iEdges = EL1.iEdges;
    ELResult.iTriFaces = EL1.iTriFaces;
    ELResult.iQuadFaces = EL1.iQuadFaces;
    ELResult.iTriCells = EL1.iTriCells;
    ELResult.iQuadCells = EL1.iQuadCells;
    ELResult.iTetCells = EL1.iTetCells;
    ELResult.iPyrCells = EL1.iPyrCells;
    ELResult.iPrismCells = EL1.iPrismCells;
    ELResult.iHexCells = EL1.iHexCells;

    std::list<Entity*>::iterator iter, iterFound, iterEnd;
    std::list<Entity*>::const_iterator iter2 = EL2.lpEntities.begin(),
	iter2End = EL2.lpEntities.end();
    for (; iter2 != iter2End; iter2++) {
      iter = ELResult.lpEntities.begin();
      iterFound = ELResult.lpEntities.end();
      iterEnd = ELResult.lpEntities.end();
      while (iter != iterEnd) {
	if (*iter == *iter2)
	  iterFound = iter;
	++iter;
      }
      if (iterFound != iterEnd) {
	Entity* pE = (*iterFound);
	switch (pE->getEntTopology())
	  {
	  case iMesh_POINT:
	    ELResult.iVerts--;
	    break;
	  case iMesh_LINE_SEGMENT:
	    ELResult.iEdges--;
	    break;
	  case iMesh_TRIANGLE:
	    if (EL1.iMeshTopoDim == 3)
	      ELResult.iTriFaces--;
	    else
	      ELResult.iTriCells--;
	    break;
	  case iMesh_QUADRILATERAL:
	    if (EL1.iMeshTopoDim == 3)
	      ELResult.iQuadFaces--;
	    else
	      ELResult.iQuadCells--;
	    break;
	  case iMesh_TETRAHEDRON:
	    ELResult.iTetCells--;
	    break;
	  case iMesh_PYRAMID:
	    ELResult.iPyrCells--;
	    break;
	  case iMesh_PRISM:
	    ELResult.iPrismCells--;
	    break;
	  case iMesh_HEXAHEDRON:
	    ELResult.iHexCells--;
	    break;
	  case iMesh_POLYGON:
	  case iMesh_POLYHEDRON:
	  case iMesh_SEPTAHEDRON:
	    throw(iBase_NOT_SUPPORTED);
	    break;
	  default:
	    // How did somebody create an entity of this topology?
	    throw(iBase_INVALID_ARGUMENT);
	  } // end switch

	ELResult.lpEntities.erase(iterFound);

      } // end if

    } // end for

      // Don't forget to update the entity set data, too.  Compute the set
      // difference of the sets that are members of this one, and add
      // member-owner links for them.
    std::set<EntitySetBase*> spESMemberDiff;
    ::vSubtract(EL1.spESMembers, EL2.spESMembers, spESMemberDiff);
    std::set<EntitySetBase*>::iterator setIter, setIterEnd =
	spESMemberDiff.end();
    for (setIter = spESMemberDiff.begin(); setIter != setIterEnd; setIter++) {
      ELResult.vAddEntitySet(*setIter);
    }

    // Build the auxiliary vertex list from scratch.
    ELResult.vRebuildAuxiliaryVertList();
  }

  void
  vIntersect(const EntityList& EL1, const EntityList& EL2, EntityList& ELResult)
  {
    std::list<Entity*> vpTemp2(EL2.lpEntities);
    std::list<Entity*>::const_iterator iter1 = EL1.lpEntities.begin(),
	iter1End = EL1.lpEntities.end();
    std::list<Entity*>::iterator iter2Found, iter2End;

    for (; iter1 != iter1End; iter1++) {
      iter2End = vpTemp2.end();
      iter2Found = std::find(vpTemp2.begin(), iter2End, *iter1);
      if (iter2Found != iter2End) {
	Entity* pE = *iter2Found;
	switch (pE->getEntTopology())
	  {
	  case iMesh_POINT:
	    ELResult.iVerts++;
	    break;
	  case iMesh_LINE_SEGMENT:
	    ELResult.iEdges++;
	    break;
	  case iMesh_TRIANGLE:
	    if (dynamic_cast<TriFace*>(pE) != NULL)
	      ELResult.iTriFaces++;
	    else
	      ELResult.iTriCells++;
	    break;
	  case iMesh_QUADRILATERAL:
	    if (dynamic_cast<QuadFace*>(pE) != NULL)
	      ELResult.iQuadFaces++;
	    else
	      ELResult.iQuadCells++;
	    break;
	  case iMesh_TETRAHEDRON:
	    ELResult.iTetCells++;
	    break;
	  case iMesh_PYRAMID:
	    ELResult.iPyrCells++;
	    break;
	  case iMesh_PRISM:
	    ELResult.iPrismCells++;
	    break;
	  case iMesh_HEXAHEDRON:
	    ELResult.iHexCells++;
	    break;
	  case iMesh_POLYGON:
	  case iMesh_POLYHEDRON:
	  case iMesh_SEPTAHEDRON:
	    throw(iBase_NOT_SUPPORTED);
	    break;
	  default:
	    // How did somebody create an entity of this topology?
	    throw(iBase_INVALID_ARGUMENT);
	  } // end switch
	ELResult.lpEntities.push_back(pE);
	vpTemp2.erase(iter2Found);
      } // end if

    } // end for

      // Don't forget to update the entity set data, too.  Compute the set
      // difference of the sets that are members of this one, and add
      // member-owner links for them.
    std::set<EntitySetBase*> spESMemberInt;
    ::vIntersect(EL1.spESMembers, EL2.spESMembers, spESMemberInt);
    std::set<EntitySetBase*>::iterator iter;
    for (iter = spESMemberInt.begin(); iter != spESMemberInt.end(); ++iter) {
      ELResult.vAddEntitySet(*iter);
    }

    // Build the auxiliary vertex list from scratch.
    ELResult.vRebuildAuxiliaryVertList();
  }

  void
  vUnite(const EntityList& EL1, const EntityList& EL2, EntityList& ELResult)
  {
    ELResult.lpEntities = EL1.lpEntities;
    ELResult.lpEntities.insert(ELResult.lpEntities.end(),
			       EL2.lpEntities.begin(), EL2.lpEntities.end());
    ELResult.iVerts = EL1.iVerts + EL2.iVerts;
    ELResult.iEdges = EL1.iEdges + EL2.iEdges;
    ELResult.iTriFaces = EL1.iTriFaces + EL2.iTriFaces;
    ELResult.iQuadFaces = EL1.iQuadFaces + EL2.iQuadFaces;
    ELResult.iTriCells = EL1.iTriCells + EL2.iTriCells;
    ELResult.iQuadCells = EL1.iQuadCells + EL2.iQuadCells;
    ELResult.iTetCells = EL1.iTetCells + EL2.iTetCells;
    ELResult.iPyrCells = EL1.iPyrCells + EL2.iPyrCells;
    ELResult.iPrismCells = EL1.iPrismCells + EL2.iPrismCells;
    ELResult.iHexCells = EL1.iHexCells + EL2.iHexCells;

    // Don't forget to update the entity set data, too.  Compute the set
    // union of the sets that are members of this one, and add
    // member-owner links for them.
    std::set<EntitySetBase*> spESMemberUnion;
    ::vUnite(EL1.spESMembers, EL2.spESMembers, spESMemberUnion);
    std::set<EntitySetBase*>::iterator iter;
    for (iter = spESMemberUnion.begin(); iter != spESMemberUnion.end();
	++iter) {
      ELResult.vAddEntitySet(*iter);
    }

    // Build the auxiliary vertex list from scratch.
    ELResult.vRebuildAuxiliaryVertList();
  }

}
