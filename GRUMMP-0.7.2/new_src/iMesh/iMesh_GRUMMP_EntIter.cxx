#include "iMesh_GRUMMP_misc.hh"
#include "iBase.h"
#include "iMesh_GRUMMP_EntIter.hh"
#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"
#include "iMesh_GRUMMP_EntitySet.hh"

namespace ITAPS_GRUMMP
{

  EntIter::EntIter(const int eEType, const int eETopo, const int iD) :
      iDim(iD), eDT(eUndefined), eExtraDT(eUndefined), pE(NULL)
  {
    {
      int err = iBase_SUCCESS;
      if (err != iBase_SUCCESS)
	throw(err);
    }

    switch (eETopo)
      {
      case iMesh_POINT:
	eDT = eVerts;
	break;
      case iMesh_LINE_SEGMENT:
	eDT = eEdgeFaces; // ITAPS / GRUMMP terminology clash.
	break;
      case iMesh_TRIANGLE:
	if (iDim == 3) {
	  eDT = eTriFaces;
	}
	else {
	  eDT = eTriCells;
	}
	break;
      case iMesh_QUADRILATERAL:
	// Quads are numbered sequentially beginning just after the
	// tris, so set counters appropriately.
	if (iDim == 3) {
	  eDT = eQuadFaces;
	}
	else {
	  eDT = eQuadCells;
	}
	break;
      case iMesh_POLYGON:
	// Quads are numbered sequentially beginning just after the
	// tris, so set counters appropriately.
	if (iDim == 3) {
	  eDT = ePolygonFaces;
	}
	else {
	  eDT = ePolygonCells;
	}
	break;
      case iMesh_TETRAHEDRON:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = eTetCells;
	break;
      case iMesh_PYRAMID:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = ePyrCells;
	break;
      case iMesh_PRISM:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = ePrismCells;
	break;
      case iMesh_HEXAHEDRON:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = eHexCells;
	break;
      case iMesh_SEPTAHEDRON:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = eSeptaCells;
	break;
      case iMesh_POLYHEDRON:
	if (iDim == 2)
	  eDT = eImpossible;
	else
	  eDT = ePolyhedra;
	break;
      case iMesh_ALL_TOPOLOGIES:
	{
	  // Choose the type of data from entity type, with no topology
	  // input.
	  switch (eEType)
	    {
	    case iBase_VERTEX:
	      eDT = eVerts;
	      break;
	    case iBase_EDGE:
	      {
		eDT = eEdgeFaces;
		break;
	      }
	    case iBase_FACE:
	      {
		if (iDim == 3) {
		  eDT = eAllFaces;
		}
		else {
		  eDT = eAllCells2D;
		}
		break;
	      }
	    case iBase_REGION:
	      {
		if (iDim == 3) {
		  eDT = eAllCells3D;
		}
		else {
		  eDT = eImpossible;
		}
		break;
	      }
	    case iBase_ALL_TYPES:
	      {
		eDT = eEverything;
		eExtraDT = eVerts;
	      }
	    } // Done with the switch over entity type
	  break;
	} // Done with the ALL_TOPOLOGIES case
      } // Done with the switch over entity topology
  }

  MeshIter::MeshIter(Mesh * const pMIn, const int eETypeIn, const int eETopoIn) :
      EntIter(eETypeIn, eETopoIn, (pMIn->getType() == Mesh::eVolMesh ? 3 : 2)), pM(
	  pMIn), iNext(0), iLast(0)
  {
    qReset();
  }

  bool
  MeshIter::qReset()
  {
    switch (eDT)
      {
      case eVerts:
	iNext = 0;
	iLast = pM->getNumVerts();
	while (iNext < iLast && pM->getVert(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case eEdgeFaces:
	if (pM->getType() == Mesh::eVolMesh) {
	  // In this case, no edges will exist, so the default bounds are
	  // okay.
	  iNext = iLast = 0;
	}
	else {
	  iNext = 0;
	  iLast = pM->getNumFaces();
	  while (iNext < iLast && pM->getFace(iNext)->isDeleted()) {
	    iNext++;
	  }
	}
	break;
      case eAllFaces:  // This is -always- 3D; see EntIter constructor
	iNext = 0;
	iLast = pM->getNumFaces();
	while (iNext < iLast && pM->getFace(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case eTriFaces:
	iNext = 0;
	iLast = dynamic_cast<const VolMesh*>(pM)->getNumTriFaces();
	while (iNext < iLast && pM->getFace(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case eQuadFaces:
	iNext = dynamic_cast<const VolMesh*>(pM)->getNumTriFaces();
	iLast = pM->getNumFaces();
	while (iNext < iLast && pM->getFace(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case eTriCells:
	{
	  const Mesh2D *pM2D = dynamic_cast<const Mesh2D*>(pM);
	  iNext = 0;
	  iLast = pM2D->getNumTriCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case eQuadCells:
	{
	  const Mesh2D *pM2D = dynamic_cast<const Mesh2D*>(pM);
	  iLast = pM2D->getNumCells();
	  iNext = pM2D->getNumTriCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case eTetCells:
	{
	  const VolMesh *pVM = dynamic_cast<const VolMesh*>(pM);
	  iNext = 0;
	  iLast = pVM->getNumTetCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case ePyrCells:
	{
	  const VolMesh *pVM = dynamic_cast<const VolMesh*>(pM);
	  iLast = pVM->getNumTetCells() + pVM->getNumPyrCells();
	  iNext = pVM->getNumTetCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case ePrismCells:
	{
	  const VolMesh *pVM = dynamic_cast<const VolMesh*>(pM);
	  iLast = pVM->getNumTetCells() + pVM->getNumPyrCells()
	      + pVM->getNumPrismCells();
	  iNext = pVM->getNumTetCells() + pVM->getNumPyrCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case eHexCells:
	{
	  const VolMesh *pVM = dynamic_cast<const VolMesh*>(pM);
	  iLast = pVM->getNumCells();
	  iNext = pVM->getNumTetCells() + pVM->getNumPyrCells()
	      + pVM->getNumPrismCells();
	  while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	    iNext++;
	  }
	  break;
	}
      case eAllCells2D:
      case eAllCells3D:
	iNext = 0;
	iLast = pM->getNumCells();
	while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case eEverything:
	// If there are no verts, there's nothing at all, so work with
	// those initially.
	eExtraDT = eVerts;
	iNext = 0;
	iLast = pM->getNumVerts();
	while (iNext < iLast && pM->getVert(iNext)->isDeleted()) {
	  iNext++;
	}
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	iNext = iLast = 0;
	break;
      case eUndefined:
	assert(0);
	break;
      } // Done with the switch over data type for workset
    pE = NULL;
    return (iNext < iLast);
  }

  bool
  MeshIter::qIncrement()
  {
    switch (eDT)
      {
      case eVerts:
	while (iNext < iLast && pM->getVert(iNext)->isDeleted()) {
	  iNext++;
	}
	if (iNext < iLast)
	  pE = pM->getVert(iNext);
	else
	  pE = NULL;
	iNext++;
	break;
      case eEdgeFaces:
      case eAllFaces:
      case eTriFaces:
      case eQuadFaces:
	while (iNext < iLast && pM->getFace(iNext)->isDeleted()) {
	  iNext++;
	}
	if (iNext < iLast)
	  pE = pM->getFace(iNext);
	else
	  pE = NULL;
	iNext++;
	break;
      case eAllCells2D:
      case eAllCells3D:
      case eTriCells:
      case eQuadCells:
      case eTetCells:
      case ePyrCells:
      case ePrismCells:
      case eHexCells:
	while (iNext < iLast && pM->getCell(iNext)->isDeleted()) {
	  iNext++;
	}
	if (iNext < iLast)
	  pE = pM->getCell(iNext);
	else
	  pE = NULL;
	iNext++;
	break;
      case eEverything:
	{
	  eDT = eExtraDT;
	  assert(
	      eDT == eVerts || eDT == eEdgeFaces || eDT == eAllFaces
		  || eDT == eAllCells2D || eDT == eAllCells3D);
	  bool qRetVal = qIncrement();
	  if (!qRetVal) {
	    if (eDT == eVerts && iDim == 3) {
	      eDT = eAllFaces; // No edges to begin with.
	      qRetVal = (qReset() && qIncrement());
	      if (!qRetVal)
		return false; // No faces means no regions.
	    }
	    else if (eDT == eVerts && iDim == 2) {
	      eDT = eEdgeFaces;
	      qRetVal = (qReset() && qIncrement());
	      if (!qRetVal)
		return false; // No edges means no faces.
	    }
	    else if (eDT == eEdgeFaces) {
	      eDT = eAllCells2D;
	      qRetVal = (qReset() && qIncrement());
	      // Nothing to go on to from here if there are no cells.
	    }
	    else if (eDT == eAllFaces) {
	      eDT = eAllCells3D;
	      qRetVal = (qReset() && qIncrement());
	      // Nothing to go on to from here if there are no cells.
	    }
	  }
	  eExtraDT = eDT;
	  eDT = eEverything;
	  break;
	}
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	// Just to be sure you don't get caught in a loop here....
	pE = NULL;
	break;
      case eUndefined:
	assert(0);
	break;
      } // Done with switch
    return (pE != NULL);
  } // qIncrement

  EntitySetIter::EntitySetIter(EntitySet * const pESIn, const int eETypeIn,
			       const int eETopoIn) :
      EntIter(eETypeIn, eETopoIn, pESIn->iMeshTopologicalDimension()), pES(
	  pESIn), itVR(), itEdgeF(), itTriF(), itQuadF(), itTriC(), itQuadC(), itTetC(), itPyrC(), itPrismC(), itHexC()
  {
    // The entity set that is being iterated over needs to be able to inform
    // the iterator class of additions or deletions which may invalidate the 
    // set iterators, so it has a set of pointers to iterators over itself
    pESIn->spIters.insert(this);
    qReset();
  }

  EntitySetIter::~EntitySetIter()
  {
    // After this iterator is ended (deleted), the entity set better not be
    // trying to inform it of additions to or deletions from the set
    pES->spIters.erase(this);
  }

  bool
  EntitySetIter::qReset()
  {
    // Setting the entity pointer to NULL tells the qIncrement method not
    // to increment the STL iterator before dereferencing to get a new pointer
    // The iterator class works this way so that it can recover if the item
    // it is pointing to is deleted and happens to be the last item
    pE = NULL;

    // Rather than nesting a bunch of if-else statements, the function
    // relies on returning to ensure that invalid steps won't be reached
    switch (eDT)
      {
      case eVerts:
	itVR = pES->sVR.begin();
	return (itVR != pES->sVR.end());
      case eEdgeFaces:
	itEdgeF = pES->spEdgeF.begin();
	return (itEdgeF != pES->spEdgeF.end());
      case eTriFaces:
	itTriF = pES->spTriF.begin();
	return (itEdgeF != pES->spEdgeF.end());
      case eQuadFaces:
	itQuadF = pES->spQuadF.begin();
	return (itQuadF != pES->spQuadF.end());
      case eAllFaces:
	eExtraDT = eTriFaces;
	itTriF = pES->spTriF.begin();
	if (itTriF != pES->spTriF.end())
	  return true;
	eExtraDT = eQuadFaces;
	itQuadF = pES->spQuadF.begin();
	return (itQuadF != pES->spQuadF.end());
      case eTriCells:
	itTriC = pES->spTriC.begin();
	return (itTriC != pES->spTriC.end());
      case eQuadCells:
	itQuadC = pES->spQuadC.begin();
	return (itQuadC != pES->spQuadC.end());
      case eTetCells:
	itTetC = pES->spTetC.begin();
	return (itTetC != pES->spTetC.end());
      case ePyrCells:
	itPyrC = pES->spPyrC.begin();
	return (itPyrC != pES->spPyrC.end());
      case ePrismCells:
	itPrismC = pES->spPrismC.begin();
	return (itPrismC != pES->spPrismC.end());
      case eHexCells:
	itHexC = pES->spHexC.begin();
	return (itHexC != pES->spHexC.end());
      case eAllCells2D:
	eExtraDT = eTriCells;
	itTriC = pES->spTriC.begin();
	if (itTriC != pES->spTriC.end())
	  return true;
	eExtraDT = eQuadCells;
	itQuadC = pES->spQuadC.begin();
	return (itQuadC != pES->spQuadC.end());
      case eAllCells3D:
	eExtraDT = eTetCells;
	itTetC = pES->spTetC.begin();
	if (itTetC != pES->spTetC.end())
	  return true;
	eExtraDT = ePyrCells;
	itPyrC = pES->spPyrC.begin();
	if (itPyrC != pES->spPyrC.end())
	  return true;
	eExtraDT = ePrismCells;
	itPrismC = pES->spPrismC.begin();
	if (itPrismC != pES->spPrismC.end())
	  return true;
	eExtraDT = eHexCells;
	itHexC = pES->spHexC.begin();
	return (itHexC != pES->spHexC.end());
      case eEverything:
	if (iDim == 2) {
	  itVR = pES->sVR.begin();
	  if (itVR != pES->sVR.end())
	    return true;

	  eExtraDT = eEdgeFaces;
	  itEdgeF = pES->spEdgeF.begin();
	  if (itEdgeF != pES->spEdgeF.end())
	    return true;

	  eExtraDT = eTriCells;
	  itTriC = pES->spTriC.begin();
	  if (itTriC != pES->spTriC.end())
	    return true;

	  eExtraDT = eQuadCells;
	  itQuadC = pES->spQuadC.begin();
	  return (itQuadC != pES->spQuadC.end());
	}
	else {
	  itVR = pES->sVR.begin();
	  if (itVR != pES->sVR.end())
	    return true;

	  eExtraDT = eTriFaces;
	  itTriF = pES->spTriF.begin();
	  if (itTriF != pES->spTriF.end())
	    return true;

	  eExtraDT = eQuadFaces;
	  itQuadF = pES->spQuadF.begin();
	  if (itQuadF != pES->spQuadF.end())
	    return true;

	  eExtraDT = eTetCells;
	  itTetC = pES->spTetC.begin();
	  if (itTetC != pES->spTetC.end())
	    return true;

	  eExtraDT = ePyrCells;
	  itPyrC = pES->spPyrC.begin();
	  if (itPyrC != pES->spPyrC.end())
	    return true;

	  eExtraDT = ePrismCells;
	  itPrismC = pES->spPrismC.begin();
	  if (itPrismC != pES->spPrismC.end())
	    return true;

	  eExtraDT = eHexCells;
	  itHexC = pES->spHexC.begin();
	  return (itHexC != pES->spHexC.end());
	}
	break;

      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	break;
      case eUndefined:
	assert(0);
	break;
      } // Done with the switch over data type for workset
    return false;
  }

  bool
  EntitySetIter::qIncrement()
  {
    // Rather than nesting a bunch of if-else statements, the function
    // relies on returning to ensure that invalid steps won't be reached
    switch (eDT)
      {
      // Only increment the STL iterator if the entity pointer is NULL
      case eVerts:
	if (itVR == pES->sVR.end())
	  return false;
	if (pE != NULL) {
	  itVR++;
	  if (itVR == pES->sVR.end())
	    return false;
	}
	pE = itVR->pVVert();
	return true;
      case eEdgeFaces:
	if (itEdgeF == pES->spEdgeF.end())
	  return false;
	if (pE != NULL) {
	  itEdgeF++;
	  if (itEdgeF == pES->spEdgeF.end())
	    return false;
	}
	pE = *itEdgeF;
	return true;
      case eTriFaces:
	if (itTriF == pES->spTriF.end())
	  return false;
	if (pE != NULL) {
	  itTriF++;
	  if (itTriF == pES->spTriF.end())
	    return false;
	}
	pE = *itTriF;
	return true;
      case eQuadFaces:
	if (itQuadF == pES->spQuadF.end())
	  return false;
	if (pE != NULL) {
	  itQuadF++;
	  if (itQuadF == pES->spQuadF.end())
	    return false;
	}
	pE = *itQuadF;
	return true;
      case eAllFaces:
	if (eExtraDT == eTriFaces) {
	  if (pE != NULL)
	    itTriF++;
	  if (itTriF != pES->spTriF.end()) {
	    pE = *itTriF;
	    return true;
	  }
	  // Move on to quad faces...
	  pE = NULL;
	  eExtraDT = eQuadFaces;
	  itQuadF = pES->spQuadF.begin();
	}
	if (itQuadF == pES->spQuadF.end())
	  return false;
	if (pE != NULL) {
	  itQuadF++;
	  if (itQuadF == pES->spQuadF.end())
	    return false;
	}
	pE = *itQuadF;
	return true;
      case eTriCells:
	if (itTriC == pES->spTriC.end())
	  return false;
	if (pE != NULL) {
	  itTriC++;
	  if (itTriC == pES->spTriC.end())
	    return false;
	}
	pE = *itTriC;
	return true;
      case eQuadCells:
	if (itQuadC == pES->spQuadC.end())
	  return false;
	if (pE != NULL) {
	  itQuadC++;
	  if (itQuadC == pES->spQuadC.end())
	    return false;
	}
	pE = *itQuadC;
	return true;
      case eTetCells:
	if (itTetC == pES->spTetC.end())
	  return false;
	if (pE != NULL) {
	  itTetC++;
	  if (itTetC == pES->spTetC.end())
	    return false;
	}
	pE = *itTetC;
	return true;
      case ePyrCells:
	if (itPyrC == pES->spPyrC.end())
	  return false;
	if (pE != NULL) {
	  itPyrC++;
	  if (itPyrC == pES->spPyrC.end())
	    return false;
	}
	pE = *itPyrC;
	return true;
      case ePrismCells:
	if (itPrismC == pES->spPrismC.end())
	  return false;
	if (pE != NULL) {
	  itPrismC++;
	  if (itPrismC == pES->spPrismC.end())
	    return false;
	}
	pE = *itPrismC;
	return true;
      case eHexCells:
	if (itHexC == pES->spHexC.end())
	  return false;
	if (pE != NULL) {
	  itHexC++;
	  if (itHexC == pES->spHexC.end())
	    return false;
	}
	pE = *itHexC;
	return true;
      case eAllCells2D:
	if (eExtraDT == eTriCells) {
	  if (pE == NULL) {
	    pE = *itTriC;
	    return true;
	  }
	  itTriC++;
	  if (itTriC != pES->spTriC.end()) {
	    pE = *itTriC;
	    return true;
	  }
	  // Move on to quad cells...
	  pE = NULL;
	  eExtraDT = eQuadCells;
	  itQuadC = pES->spQuadC.begin();
	}
	if (itQuadC == pES->spQuadC.end())
	  return false;
	if (pE != NULL) {
	  itQuadC++;
	  if (itQuadC == pES->spQuadC.end())
	    return false;
	}
	pE = *itQuadC;
	return true;
      case eAllCells3D:
	if (eExtraDT == eTetCells) {
	  if (pE == NULL) {
	    pE = *itTetC;
	    return true;
	  }
	  itTetC++;
	  if (itTetC != pES->spTetC.end()) {
	    pE = *itTetC;
	    return true;
	  }
	  // Move on to pyramids...
	  pE = NULL;
	  eExtraDT = ePyrCells;
	  itPyrC = pES->spPyrC.begin();
	}
	if (eExtraDT == ePyrCells) {
	  if (itPyrC != pES->spPyrC.end()) {
	    if (pE == NULL) {
	      pE = *itPyrC;
	      return true;
	    }
	    itPyrC++;
	    if (itPyrC != pES->spPyrC.end()) {
	      pE = *itPyrC;
	      return true;
	    }
	  }
	  // Move on to prisms...
	  pE = NULL;
	  eExtraDT = ePrismCells;
	  itPrismC = pES->spPrismC.begin();
	}
	if (eExtraDT == ePrismCells) {
	  if (itPrismC != pES->spPrismC.end()) {
	    if (pE == NULL) {
	      pE = *itPrismC;
	      return true;
	    }
	    itPrismC++;
	    if (itPrismC != pES->spPrismC.end()) {
	      pE = *itPrismC;
	      return true;
	    }
	  }
	  // Move on to hexes...
	  pE = NULL;
	  eExtraDT = eHexCells;
	  itHexC = pES->spHexC.begin();
	}
	// Must be eHexCells
	if (itHexC == pES->spHexC.end())
	  return false;
	if (pE != NULL) {
	  itHexC++;
	  if (itHexC == pES->spHexC.end())
	    return false;
	}
	pE = *itHexC;
	return true;
      case eEverything:
	if (iDim == 2) {
	  switch (eExtraDT)
	    {
	    case eVerts:
	      if (pE != NULL) {
		itVR++;
	      }
	      if (itVR != pES->sVR.end()) {
		pE = itVR->pVVert();
		return true;
	      }
	      // Move on to edges
	      pE = NULL;
	      eExtraDT = eEdgeFaces;
	      itEdgeF = pES->spEdgeF.begin();
	      // fall through
	    case eEdgeFaces:
	      if (pE != NULL)
		itEdgeF++;
	      if (itEdgeF != pES->spEdgeF.end()) {
		pE = *itEdgeF;
		return true;
	      }
	      // Move on to tri cells...
	      pE = NULL;
	      eExtraDT = eTriCells;
	      itTriC = pES->spTriC.begin();
	      // fall through
	    case eTriCells:
	      if (pE != NULL)
		itTriC++;
	      if (itTriC != pES->spTriC.end()) {
		pE = *itTriC;
		return true;
	      }
	      // Move on to quad cells...
	      pE = NULL;
	      eExtraDT = eQuadCells;
	      itQuadC = pES->spQuadC.begin();
	      // fall through
	    case eQuadCells:
	      if (pE != NULL)
		itQuadC++;
	      if (itQuadC != pES->spQuadC.end()) {
		pE = *itQuadC;
		return true;
	      }
	      return false;
	    case eTetCells:
	    case ePyrCells:
	    case ePrismCells:
	    case eHexCells:
	    default:
	      assert(0);
	    } // Done with switch to retrieve Everything in 2D
	} // Done with 2D
	else {
	  assert(iDim == 3);
	  switch (eExtraDT)
	    {
	    case eVerts:
	      if (pE != NULL) {
		itVR++;
	      }
	      if (itVR != pES->sVR.end()) {
		pE = itVR->pVVert();
		return true;
	      }
	      // Move on to tri faces...
	      pE = NULL;
	      eExtraDT = eTriFaces;
	      itTriF = pES->spTriF.begin();
	      // fall through
	    case eTriFaces:
	      if (pE != NULL)
		itTriF++;
	      if (itTriF != pES->spTriF.end()) {
		pE = *itTriF;
		return true;
	      }
	      // Move on to quad faces...
	      pE = NULL;
	      eExtraDT = eQuadFaces;
	      itQuadF = pES->spQuadF.begin();
	      // fall through
	    case eQuadFaces:
	      if (pE != NULL)
		itQuadF++;
	      if (itQuadF != pES->spQuadF.end()) {
		pE = *itQuadF;
		return true;
	      }
	      // Move on to tet cells...
	      pE = NULL;
	      eExtraDT = eTetCells;
	      itTetC = pES->spTetC.begin();
	      // fall through
	    case eTetCells:
	      if (pE != NULL)
		itTetC++;
	      if (itTetC != pES->spTetC.end()) {
		pE = *itTetC;
		return true;
	      }
	      // Move on to pyramids...
	      pE = NULL;
	      eExtraDT = ePyrCells;
	      itPyrC = pES->spPyrC.begin();
	      // fall through
	    case ePyrCells:
	      if (itPyrC != pES->spPyrC.end()) {
		if (pE == NULL) {
		  pE = *itPyrC;
		  return true;
		}
		itPyrC++;
		if (itPyrC != pES->spPyrC.end()) {
		  pE = *itPyrC;
		  return true;
		}
	      }
	      // Move on to prisms...
	      pE = NULL;
	      eExtraDT = ePrismCells;
	      itPrismC = pES->spPrismC.begin();
	      // fall through
	    case ePrismCells:
	      if (itPrismC != pES->spPrismC.end()) {
		if (pE == NULL) {
		  pE = *itPrismC;
		  return true;
		}
		itPrismC++;
		if (itPrismC != pES->spPrismC.end()) {
		  pE = *itPrismC;
		  return true;
		}
	      }
	      // Move on to hexes...
	      pE = NULL;
	      eExtraDT = eHexCells;
	      itHexC = pES->spHexC.begin();
	      // fall through
	    case eHexCells:
	      if (itHexC == pES->spHexC.end())
		return false;
	      if (pE != NULL) {
		itHexC++;
		if (itHexC == pES->spHexC.end())
		  return false;
	      }
	      pE = *itHexC;
	      return true;
	    default:
	      return false;
	    } // Done with switch to retrieve Everything in 3D
	} // Done with 3D
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	// Do absolutely nothing, as there won't ever be any of these.
	return false;
      case eUndefined:
	assert(0); // Some kind of data error must have occurred...
	return false;
      } // Done with switch
    return false;
  } // qIncrement

  void
  EntitySetIter::vDelete(Entity* pEToDelete)
  {
    // The STL iterator should only be invalidated if it happens to
    // be pointing to the entity to be deleted
    if (pE == pEToDelete) {
      // Flag the entity handle stored by the iterator and increment
      pE = NULL;
      switch (eDT)
	{
	case eVerts:
	  itVR++;
	  break;
	case eEdgeFaces:
	  itEdgeF++;
	  break;
	case eTriFaces:
	  itTriF++;
	  break;
	case eQuadFaces:
	  itQuadF++;
	  break;
	case eAllFaces:
	  if (eExtraDT == eTriFaces) {
	    itTriF++;
	    if (itTriF == pES->spTriF.end()) {
	      eExtraDT = eQuadFaces;
	      itQuadF = pES->spQuadF.begin();
	    }
	  }
	  else {
	    itQuadF++;
	  }
	  break;
	case eTriCells:
	  itTriC++;
	  break;
	case eQuadCells:
	  itQuadC++;
	  break;
	case eTetCells:
	  itTetC++;
	  break;
	case ePyrCells:
	  itPyrC++;
	  break;
	case ePrismCells:
	  itPrismC++;
	  break;
	case eHexCells:
	  itHexC++;
	  break;
	case eAllCells2D:
	  if (eExtraDT == eTriCells) {
	    itTriC++;
	    if (itTriC == pES->spTriC.end()) {
	      eExtraDT = eQuadCells;
	      itQuadC = pES->spQuadC.begin();
	    }
	  }
	  else {
	    itQuadC++;
	  }
	  break;
	case eAllCells3D:
	  if (eExtraDT == eTetCells) {
	    itTetC++;
	    if (itTetC == pES->spTetC.end()) {
	      eExtraDT = ePyrCells;
	      itPyrC = pES->spPyrC.begin();
	      if (itPyrC == pES->spPyrC.end()) {
		eExtraDT = ePrismCells;
		itPrismC = pES->spPrismC.begin();
		if (itPrismC == pES->spPrismC.end()) {
		  eExtraDT = eHexCells;
		  itHexC = pES->spHexC.begin();
		} // if ( itPrismC...
	      } // if ( itPyrC...
	    }
	  }
	  else if (eExtraDT == ePyrCells) {
	    itPyrC++;
	    if (itPyrC == pES->spPyrC.end()) {
	      eExtraDT = ePrismCells;
	      itPrismC = pES->spPrismC.begin();
	      if (itPrismC == pES->spPrismC.end()) {
		eExtraDT = eHexCells;
		itHexC = pES->spHexC.begin();
	      } // if ( itPrismC...
	    } // if ( itPyrC...
	  }
	  else if (eExtraDT == ePrismCells) {
	    itPrismC++;
	    if (itPrismC == pES->spPrismC.end()) {
	      eExtraDT = eHexCells;
	      itHexC = pES->spHexC.begin();
	    } // if ( itPrismC...
	  }
	  else if (eExtraDT == eHexCells) {
	    itHexC++;
	  }
	  break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
	  // Do absolutely nothing, as there won't ever be any of these.
	  break;
	case eEverything: // The sub-iterator will be in use, not the top one.
	case eUndefined:
	  assert(0);
	  break;
	} // Done with switch
    }  // Done with if
  } // vDelete

  EntityListIter::EntityListIter(EntityList * const pELIn, const int eETypeIn,
				 const int eETopoIn) :
      EntIter(eETypeIn, eETopoIn, pELIn->iMeshTopologicalDimension()), pEL(
	  pELIn), iFound(0), iTotal(0), iter()
  {
    // The entity list that is being iterated over needs to be able to inform
    // the iterator class of additions or deletions which may invalidate the 
    // list iterator, so it has a list of pointers to iterators over itself
    pELIn->spIters.insert(this);
    qReset();
  }

  EntityListIter::~EntityListIter()
  {
    // After this iterator is ended (deleted), the entity set better not be
    // trying to inform it of additions to or deletions from the set
    pEL->spIters.erase(this);
  }

  bool
  EntityListIter::qReset()
  {
    // As with entity set iterators, setting the entity pointer to NULL
    // tells the qIncrement method not to increment the STL list iterator
    // before dereferencing it
    // The iterator class works this way so that it can recover if the item
    // it is pointing to is deleted and happens to be the last item
    // List iterators also store the number of entities of that type/topology
    // to be found and the number already found so that it can stop when it
    // has found all that it will
    pE = NULL;
    iFound = 0;
    iter = pEL->lpEntities.begin();

    // The entity list also stores how many of each type/topology it contains
    switch (eDT)
      {
      case eVerts:
	iTotal = pEL->iVerts;
	break;
      case eEdgeFaces:
	iTotal = pEL->iEdges;
	break;
      case eTriFaces:
	iTotal = pEL->iTriFaces;
	break;
      case eQuadFaces:
	iTotal = pEL->iQuadFaces;
	break;
      case eAllFaces:
	iTotal = pEL->iTriFaces + pEL->iQuadFaces;
	break;
      case eAllCells2D:
	iTotal = pEL->iTriCells + pEL->iQuadCells;
	break;
      case eAllCells3D:
	iTotal = pEL->iTetCells + pEL->iPyrCells + pEL->iPrismCells
	    + pEL->iHexCells;
	break;
      case eTriCells:
	iTotal = pEL->iTriCells;
	break;
      case eQuadCells:
	iTotal = pEL->iQuadCells;
	break;
      case eTetCells:
	iTotal = pEL->iTetCells;
	break;
      case ePyrCells:
	iTotal = pEL->iPyrCells;
	break;
      case ePrismCells:
	iTotal = pEL->iPrismCells;
	break;
      case eHexCells:
	iTotal = pEL->iHexCells;
	break;
      case eEverything:
	iTotal = pEL->lpEntities.size();
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	iTotal = 0;
	return false;
	break;
      case eUndefined:
	assert(0);
	return false;
	break;
      } // Done with switch

    assert((iter != pEL->lpEntities.end()) || (iTotal == 0));
    return (iFound < iTotal);
  }

  bool
  EntityListIter::qIncrement()
  {
    // If all have been found, don't bother looking anymore
    if (iFound >= iTotal)
      return false;
    assert(iter != pEL->lpEntities.end());

    // Increment the STL iterator if the entity pointer is non-NULL
    // (i.e., if we haven't just reset).  
    if (pE != NULL)
      ++iter;

    // Find the next entity of the correct type/topology
    switch (eDT)
      {
      case eVerts:
	while ((*iter)->getEntType() != iBase_VERTEX)
	  ++iter;
	break;
      case eAllFaces:
      case eAllCells2D:
	while ((*iter)->getEntType() != iBase_FACE)
	  ++iter;
	break;
      case eEdgeFaces:
	while ((*iter)->getEntType() != iBase_EDGE)
	  ++iter;
	break;
      case eTriFaces:
      case eTriCells:
	while ((*iter)->getEntTopology() != iMesh_TRIANGLE)
	  ++iter;
	break;
      case eQuadFaces:
      case eQuadCells:
	while ((*iter)->getEntTopology() != iMesh_QUADRILATERAL)
	  ++iter;
	break;
      case eTetCells:
	while ((*iter)->getEntTopology() != iMesh_TETRAHEDRON)
	  ++iter;
	break;
      case ePyrCells:
	while ((*iter)->getEntTopology() != iMesh_PYRAMID)
	  ++iter;
	break;
      case ePrismCells:
	while ((*iter)->getEntTopology() != iMesh_PRISM)
	  ++iter;
	break;
      case eHexCells:
	while ((*iter)->getEntTopology() != iMesh_HEXAHEDRON)
	  ++iter;
	break;
      case eAllCells3D:
	while ((*iter)->getEntType() != iBase_REGION)
	  ++iter;
	break;
      case eEverything:
	// Actually do nothing here...
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	return false;
	break;
      case eUndefined:
	assert(0);
	return false;
	break;
      } // Done with switch

    iFound++;
    assert(iter != pEL->lpEntities.end());
    // Get the new entity pointer and return true
    pE = (*iter);
    return true;
  } // qIncrement

  void
  EntityListIter::vAdjustIterForDeletionByType(
      std::list<Entity*>::iterator itToDelete, const int iType)
  {
    std::list<Entity*>::iterator iterTemp, iterPrev;
    if (iter == itToDelete) {
      if (iter == pEL->lpEntities.begin()) {
	if ((*iter)->getEntType() == iType)
	  iTotal--;
	++iter;
      }
      else {
	iTotal--;
	iFound--;
	if (iFound == 0) {
	  iter = pEL->lpEntities.begin();
	  pE = NULL;
	}
	else {
	  for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter; iterTemp++)
	    if ((*iterTemp)->getEntType() == iType)
	      iterPrev = iterTemp;
	  iter = iterPrev;
	}
      }
    }
    else if ((*itToDelete)->getEntType() == iType) {
      iTotal--;
      iterTemp = pEL->lpEntities.begin();
      while (iterTemp != iter && iterTemp != itToDelete)
	iterTemp++;
      if (iterTemp != iter)
	iFound--;
    }
  }

  void
  EntityListIter::vAdjustIterForDeletionByTopo(
      std::list<Entity*>::iterator itToDelete, const int iTopo)
  {
    std::list<Entity*>::iterator iterTemp, iterPrev;
    if (iter == itToDelete) {
      if (iter == pEL->lpEntities.begin()) {
	if ((*iter)->getEntTopology() == iTopo)
	  iTotal--;
	++iter;
      }
      else {
	iTotal--;
	iFound--;
	if (iFound == 0) {
	  iter = pEL->lpEntities.begin();
	  pE = NULL;
	}
	else {
	  for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter; iterTemp++)
	    if ((*iterTemp)->getEntTopology() == iTopo)
	      iterPrev = iterTemp;
	  iter = iterPrev;
	}
      }
    }
    else if ((*itToDelete)->getEntTopology() == iTopo) {
      iTotal--;
      iterTemp = pEL->lpEntities.begin();
      while (iterTemp != iter && iterTemp != itToDelete)
	iterTemp++;
      if (iterTemp != iter)
	iFound--;
    }
  }

  void
  EntityListIter::vDelete(std::list<Entity*>::iterator itToDelete)
  {
    // May seem overly complicated, but if the entity pointed to by the STL iterator
    // is going to be deleted, must recover the previous entity of that type/topology
    // Two possibilities: EITHER 
    // (1) the entity is the type/topology or the iterator hasn't yet been incremented
    //     and the first item is the one to be removed ( and it might just happen to be 
    //     the right type/topology ) OR
    // (2) the entity is the correct type/topology in which case the previous entity
    //     of that type/topology needs to be recovered
    // Even if that entity is not the one to be deleted, need to check that the found
    // and total counts are up to date
    switch (eDT)
      {
      case eVerts:
	vAdjustIterForDeletionByType(itToDelete, iBase_VERTEX);
	break;
      case eAllFaces:
      case eAllCells2D:
	vAdjustIterForDeletionByType(itToDelete, iBase_FACE);
	break;
      case eEdgeFaces:
	vAdjustIterForDeletionByType(itToDelete, iBase_EDGE);
	break;
      case eTriFaces:
      case eTriCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_TRIANGLE);
	break;
      case eQuadFaces:
      case eQuadCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_QUADRILATERAL);
	break;
      case eTetCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_TETRAHEDRON);
	break;
      case ePyrCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_PYRAMID);
	break;
      case ePrismCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_PRISM);
	break;
      case eHexCells:
	vAdjustIterForDeletionByTopo(itToDelete, iMesh_HEXAHEDRON);
	break;
      case eAllCells3D:
	vAdjustIterForDeletionByType(itToDelete, iBase_REGION);
	break;
      case eEverything:
	{
	  std::list<Entity*>::iterator iterTemp, iterPrev;
	  if (iter == itToDelete) {
	    if (iter == pEL->lpEntities.begin()) {
	      iTotal--;
	      ++iter;
	    }
	    else {
	      iTotal--;
	      iFound--;
	      if (iFound == 0) {
		iter = pEL->lpEntities.begin();
		pE = NULL;
	      }
	      else {
		for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter;
		    iterTemp++)
		  iterPrev = iterTemp;
		iter = iterPrev;
	      }
	    }
	  }
	  else {
	    iTotal--;
	    iterTemp = pEL->lpEntities.begin();
	    while (iterTemp != iter && iterTemp != itToDelete)
	      iterTemp++;
	    if (iterTemp != iter)
	      iFound--;
	  }
	}
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	break;
      case eUndefined:
	assert(0);
	break;
      } // Done with switch

  } // vDelete

  void
  EntityListIter::vAdded(Entity* pEAdded)
  {
    // One special case where the STL iterator is pointing to the beginning of
    // an empty list (which will also be the end).  In that case the STL iterator
    // needs to be reset so it can at least pickup the entity that was added
    if (iFound == 0 && iter == pEL->lpEntities.end()) {
      iter = pEL->lpEntities.begin();
      pE = NULL;
    }

    // The entity will always be added at the end and as long as the STL iterator
    // isn't pointing to the end (see above) it won't miss the new entity, so just
    // check if the new entity is of the desired type/topology
    switch (eDT)
      {
      case eVerts:
	if (pEAdded->getEntType() == iBase_VERTEX)
	  iTotal++;
	break;
      case eEdgeFaces:
	if (pEAdded->getEntType() == iBase_EDGE)
	  iTotal++;
	break;
      case eTriFaces:
      case eTriCells:
	if (pEAdded->getEntTopology() == iMesh_TRIANGLE)
	  iTotal++;
	break;
      case eQuadFaces:
      case eQuadCells:
	if (pEAdded->getEntTopology() == iMesh_QUADRILATERAL)
	  iTotal++;
	break;
      case eAllFaces:
      case eAllCells2D:
	if (pEAdded->getEntType() == iBase_FACE)
	  iTotal++;
	break;
      case eAllCells3D:
	if (pEAdded->getEntType() == iBase_REGION)
	  iTotal++;
	break;
      case eTetCells:
	if (pEAdded->getEntTopology() == iMesh_TETRAHEDRON)
	  iTotal++;
	break;
      case ePyrCells:
	if (pEAdded->getEntTopology() == iMesh_PYRAMID)
	  iTotal++;
	break;
      case ePrismCells:
	if (pEAdded->getEntTopology() == iMesh_PRISM)
	  iTotal++;
	break;
      case eHexCells:
	if (pEAdded->getEntTopology() == iMesh_HEXAHEDRON)
	  iTotal++;
	break;
      case ePolygonCells:
      case ePolygonFaces:
      case ePolyhedra:
      case eSeptaCells:
      case eImpossible:
	break;
      case eEverything: // The sub-iterator will be in use, not the top one.
      case eUndefined:
	assert(0);
	break;
      } // Done with switch      
  } // vAdded

} // Close namespace

