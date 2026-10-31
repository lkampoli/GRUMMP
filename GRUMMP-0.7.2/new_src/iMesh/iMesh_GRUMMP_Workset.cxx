#include "iMesh_GRUMMP_misc.hh"
#include "iMesh_GRUMMP_EntitySet.hh"
#include "iMesh_GRUMMP_Workset.hh"
#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"

namespace ITAPS_GRUMMP {

Workset::Workset(const int eEType, const int eETopo, const int iSzIn,
		const int iD) :
		iSize(iSzIn), iDim(iD), eDT(eUndefined), eCurrentDT(eUndefined) {
	{
		int err = iBase_SUCCESS;
		if (err != iBase_SUCCESS)
			throw(err);
	}
	if (iSize < 0)
		throw(iBase_INVALID_ARGUMENT);

	switch (eETopo) {
	case iMesh_POINT:
		eDT = eCurrentDT = eVerts;
		break;
	case iMesh_LINE_SEGMENT:
		eDT = eCurrentDT = eEdgeFaces; // ITAPS / RefImpl terminology clash.
		break;
	case iMesh_TRIANGLE:
		if (iDim == 3) {
			eDT = eCurrentDT = eTriFaces;
		} else {
			eDT = eCurrentDT = eTriCells;
		}
		break;
	case iMesh_QUADRILATERAL:
		// Quads are numbered sequentially beginning just after the
		// tris, so set counters appropriately.
		if (iDim == 3) {
			eDT = eCurrentDT = eQuadFaces;
		} else {
			eDT = eCurrentDT = eQuadCells;
		}
		break;
	case iMesh_POLYGON:
		// Polygons are numbered sequentially beginning just after the
		// quads, so set counters appropriately.
		if (iDim == 3) {
			eDT = eCurrentDT = ePolygonFaces;
		} else {
			eDT = eCurrentDT = ePolygonCells;
		}
		break;
	case iMesh_TETRAHEDRON:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = eTetCells;
		break;
	case iMesh_PYRAMID:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = ePyrCells;
		break;
	case iMesh_PRISM:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = ePrismCells;
		break;
	case iMesh_HEXAHEDRON:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = eHexCells;
		break;
	case iMesh_SEPTAHEDRON:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = eSeptaCells;
		break;
	case iMesh_POLYHEDRON:
		if (iDim == 2)
			eDT = eCurrentDT = eImpossible;
		else
			eDT = eCurrentDT = ePolyhedra;
		break;
	case iMesh_ALL_TOPOLOGIES: {
		// Choose the type of data from entity type, with no topology
		// input.
		switch (eEType) {
		case iBase_VERTEX:
			eDT = eCurrentDT = eVerts;
			break;
		case iBase_EDGE: {
			eDT = eCurrentDT = eEdgeFaces;
			break;
		}
		case iBase_FACE: {
			if (iDim == 3) {
				eDT = eAllFaces;
				eCurrentDT = eTriFaces;
			} else {
				eDT = eAllCells2D;
				eCurrentDT = eTriCells;
			}
			break;
		}
		case iBase_REGION: {
			if (iDim == 3) {
				eDT = eAllCells3D;
				eCurrentDT = eTetCells;
			} else {
				eDT = eCurrentDT = eImpossible;
			}
			break;
		}
		case iBase_ALL_TYPES: {
			eDT = eEverything;
			eCurrentDT = eVerts;
		}
		} // Done with the switch over entity type
		break;
	} // Done with the ALL_TOPOLOGIES case
	} // Done with the switch over entity topology
}

MeshWorkset::MeshWorkset(Mesh * const pMIn, const int eETypeIn,
		const int eETopoIn, const int iSzIn) :
		Workset(eETypeIn, eETopoIn, iSzIn,
				(pMIn->getType() == Mesh::eVolMesh ? 3 : 2)), pM(pMIn), iNext(
				0), iTotal(0) {
	qReset();
}

bool MeshWorkset::qReset() {
	iNext = 0;
	switch (eDT) {
	case eVerts:
		iTotal = pM->getNumVerts();
		break;
	case eEdgeFaces:
		if (pM->getType() == Mesh::eVolMesh) {
			// In this case, no edges will exist, so the default bounds are
			// okay.
			break;
		} else {
			// Both surface and 2D meshes will be fine with this.
			iTotal = pM->getNumFaces();
			break;
		}
	case eAllFaces:
		iTotal = pM->getNumFaces();
		break;
	case eTriFaces:
		iTotal = dynamic_cast<VolMesh*>(pM)->getNumTriFaces();
		break;
	case eQuadFaces:
		iNext = dynamic_cast<VolMesh*>(pM)->getNumTriFaces();
		iTotal = pM->getNumFaces();
		break;
	case eTriCells: {
		Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pM);
		iTotal = pM2D->getNumTriCells();
		break;
	}
	case eQuadCells: {
		Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pM);
		iTotal = pM2D->getNumCells();
		iNext = pM2D->getNumTriCells();
		break;
	}
	case eTetCells: {
		VolMesh *pVM = dynamic_cast<VolMesh*>(pM);
		iTotal = pVM->getNumTetCells();
		break;
	}
	case ePyrCells: {
		VolMesh *pVM = dynamic_cast<VolMesh*>(pM);
		iTotal = pVM->getNumTetCells() + pVM->getNumPyrCells();
		iNext = pVM->getNumTetCells();
		break;
	}
	case ePrismCells: {
		VolMesh *pVM = dynamic_cast<VolMesh*>(pM);
		iTotal = pVM->getNumTetCells() + pVM->getNumPyrCells()
				+ pVM->getNumPrismCells();
		iNext = pVM->getNumTetCells() + pVM->getNumPyrCells();
		break;
	}
	case eHexCells: {
		VolMesh *pVM = dynamic_cast<VolMesh*>(pM);
		iTotal = pVM->getNumCells();
		iNext = pVM->getNumTetCells() + pVM->getNumPyrCells()
				+ pVM->getNumPrismCells();
		break;
	}
	case eAllCells2D:
	case eAllCells3D:
		iTotal = pM->getNumCells();
		break;
	case eEverything:
		// If there are no verts, there's nothing at all, so work with
		// those initially.
		eCurrentDT = eVerts;
		iTotal = pM->getNumVerts();
		break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
		iTotal = 0;
		break;
	case eUndefined:
		assert(0);
		break;
	} // Done with the switch over data type for workset
	return (iNext + iSize < iTotal);
}

bool MeshWorkset::qGetNextBlock(iBase_EntityHandle** entity_handles,
		int *entity_handles_allocated, int *entity_handles_size)

		{
	bool qRetVal = false;
	vEnsureArraySize<iBase_EntityHandle>(entity_handles,
			entity_handles_allocated, iSize);
	int iEH = 0;

	switch (eDT) {
	case eVerts:
		for (; iEH < iSize && iNext < iTotal; iNext++) {
			Vert *pV = pM->getVert(iNext);
			if (!pV->isDeleted())
				(*entity_handles)[iEH++] = pV;
		}
		qRetVal = (iNext < iTotal);
		break;
	case eEdgeFaces:
	case eAllFaces:
	case eTriFaces:
	case eQuadFaces:
		for (; iEH < iSize && iNext < iTotal; iNext++) {
			Face *pF = pM->getFace(iNext);
			if (!pF->isDeleted())
				(*entity_handles)[iEH++] = pF;
		}
		qRetVal = (iNext < iTotal);
		break;
	case eAllCells2D:
	case eAllCells3D:
	case eTriCells:
	case eQuadCells:
	case eTetCells:
	case ePyrCells:
	case ePrismCells:
	case eHexCells:
		for (; iEH < iSize && iNext < iTotal; iNext++) {
			Cell* pC = pM->getCell(iNext);
			if (!pC->isDeleted())
				(*entity_handles)[iEH++] = pC;
		}
		qRetVal = (iNext < iTotal);
		break;
	case eEverything:
		qRetVal = false;
		if (eCurrentDT == eVerts) {
			for (; iEH < iSize && iNext < iTotal; iNext++) {
				Vert *pV = pM->getVert(iNext);
				if (!pV->isDeleted())
					(*entity_handles)[iEH++] = pV;
			}
			qRetVal = (iNext < iTotal);
			if (!qRetVal) {
				// Try to find what's next.
				if (iDim == 2) {
					eCurrentDT = eEdgeFaces;
					iNext = 0;
					iTotal = pM->getNumFaces();
					if (iTotal == 0)
						return false; // No edges means no faces.
				} else {
					assert(iDim == 3);
					eCurrentDT = eAllFaces;
					iNext = 0;
					iTotal = pM->getNumFaces();
					if (iTotal == 0)
						return false; // No faces means no regions.
				}
			} // Done shifting from verts to edges/faces
		} // Done adding verts

		if (!qRetVal && (eCurrentDT == eEdgeFaces || eCurrentDT == eAllFaces)) {
			for (; iEH < iSize && iNext < iTotal; iNext++) {
				Face *pF = pM->getFace(iNext);
				if (!pF->isDeleted())
					(*entity_handles)[iEH++] = pF;
			}
			qRetVal = (iNext < iTotal);
			if (!qRetVal) {
				// Shift to whatever's next and continue.
				if (eCurrentDT == eEdgeFaces) {
					eCurrentDT = eAllCells2D;
					iNext = 0;
					iTotal = pM->getNumCells();
					if (iTotal == 0)
						return false; // There aren't any 2D faces.
				} else {
					assert(eCurrentDT == eAllFaces);
					eCurrentDT = eAllCells3D;
					iNext = 0;
					iTotal = pM->getNumCells();
					if (iTotal == 0)
						return false; // There aren't any 2D faces.
				}
			} // Done shifting from GRUMMP faces to GRUMMP cells
		} // Done adding GRUMMP faces

		if (!qRetVal
				&& (eCurrentDT == eAllCells2D || eCurrentDT == eAllCells3D)) {
			for (; iEH < iSize && iNext < iTotal; iNext++) {
				Cell* pC = pM->getCell(iNext);
				if (!pC->isDeleted())
					(*entity_handles)[iEH++] = pC;
			}
			qRetVal = (iNext < iTotal);
			// No possibility to go on to anything else.
		}
		break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
		// Do absolutely nothing, as there won't ever be any of these.
		*entity_handles_size = 0;
		qRetVal = false;
		break;
	case eUndefined:
		assert(0);
		*entity_handles_size = 0;
		qRetVal = false;
		break;
	} // Done with switch
	*entity_handles_size = iEH;
	return ((*entity_handles_size) > 0);
} // qGetNextBlock

EntitySetWorkset::EntitySetWorkset(EntitySet * const pESIn, const int eETypeIn,
		const int eETopoIn, const int iSzIn) :
		Workset(eETypeIn, eETopoIn, iSzIn, pESIn->iMeshTopologicalDimension()), pES(
				pESIn), itVR(), itEdgeF(), itTriF(), itQuadF(), itTriC(), itQuadC(), itTetC(), itPyrC(), itPrismC(), itHexC(), itVREnd(), itEdgeFEnd(), itTriFEnd(), itQuadFEnd(), itTriCEnd(), itQuadCEnd(), itTetCEnd(), itPyrCEnd(), itPrismCEnd(), itHexCEnd() {
	assert(pESIn);
	// The entity set that is being iterated over needs to be able to inform
	// the iterator class of additions or deletions which may invalidate the
	// set iterators, so it has a set of pointers to iterators over itself
	pESIn->spWorksets.insert(this);
	qReset();
}

EntitySetWorkset::~EntitySetWorkset() {
	// After this iterator is ended (deleted), the entity set better not be
	// trying to inform it of additions to or deletions from the set
	pES->spWorksets.erase(this);
}

bool EntitySetWorkset::qReset() {
	// Rather than nesting a bunch of if-else statements, the function
	// relies on returning to ensure that invalid steps won't be reached
	switch (eDT) {
	case eVerts:
		itVR = pES->sVR.begin();
		itVREnd = pES->sVR.end();
		return (itVR != itVREnd);
	case eEdgeFaces:
		itEdgeF = pES->spEdgeF.begin();
		itEdgeFEnd = pES->spEdgeF.end();
		return (itEdgeF != itEdgeFEnd);
	case eTriFaces:
		itTriF = pES->spTriF.begin();
		itTriFEnd = pES->spTriF.end();
		return (itTriF != itTriFEnd);
	case eQuadFaces:
		itQuadF = pES->spQuadF.begin();
		itQuadFEnd = pES->spQuadF.end();
		return (itQuadF != itQuadFEnd);
	case eAllFaces:
		eCurrentDT = eTriFaces;
		itTriF = pES->spTriF.begin();
		itTriFEnd = pES->spTriF.end();
		if (itTriF != itTriFEnd)
			return true;
		eCurrentDT = eQuadFaces;
		itQuadF = pES->spQuadF.begin();
		itQuadFEnd = pES->spQuadF.end();
		return (itQuadF != itQuadFEnd);
	case eAllCells2D:
		eCurrentDT = eTriCells;
		itTriC = pES->spTriC.begin();
		itTriCEnd = pES->spTriC.end();
		if (itTriC != itTriCEnd)
			return true;
		eCurrentDT = eQuadCells;
		itQuadC = pES->spQuadC.begin();
		itQuadCEnd = pES->spQuadC.end();
		return (itQuadC != itQuadCEnd);
	case eAllCells3D:
		eCurrentDT = eTetCells;
		itTetC = pES->spTetC.begin();
		itTetCEnd = pES->spTetC.end();
		if (itTetC != itTetCEnd)
			return true;
		eCurrentDT = ePyrCells;
		itPyrC = pES->spPyrC.begin();
		itPyrCEnd = pES->spPyrC.end();
		if (itPyrC != itPyrCEnd)
			return true;
		eCurrentDT = ePrismCells;
		itPrismC = pES->spPrismC.begin();
		itPrismCEnd = pES->spPrismC.end();
		if (itPrismC != itPrismCEnd)
			return true;
		eCurrentDT = eHexCells;
		itHexC = pES->spHexC.begin();
		itHexCEnd = pES->spHexC.end();
		return (itHexC != itHexCEnd);
	case eTriCells:
		itTriC = pES->spTriC.begin();
		itTriCEnd = pES->spTriC.end();
		return (itTriC != itTriCEnd);
	case eQuadCells:
		itQuadC = pES->spQuadC.begin();
		itQuadCEnd = pES->spQuadC.end();
		return (itQuadC != itQuadCEnd);
	case eTetCells:
		itTetC = pES->spTetC.begin();
		itTetCEnd = pES->spTetC.end();
		return (itTetC != itTetCEnd);
	case ePyrCells:
		itPyrC = pES->spPyrC.begin();
		itPyrCEnd = pES->spPyrC.end();
		return (itPyrC != itPyrCEnd);
	case ePrismCells:
		itPrismC = pES->spPrismC.begin();
		itPrismCEnd = pES->spPrismC.end();
		return (itPrismC != itPrismCEnd);
	case eHexCells:
		itHexC = pES->spHexC.begin();
		itHexCEnd = pES->spHexC.end();
		return (itHexC != itHexCEnd);
	case eEverything:
		eCurrentDT = eVerts;
		itVR = pES->sVR.begin();
		itVREnd = pES->sVR.end();
		return (itVR != itVREnd);
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
	return false;
}

bool EntitySetWorkset::qGetNextBlock(iBase_EntityHandle** entity_handles,
		int *entity_handles_allocated, int *entity_handles_size)

		{
	vEnsureArraySize<iBase_EntityHandle>(entity_handles,
			entity_handles_allocated, iSize);
	int iEH, i;

	switch (eDT) {
	case eVerts:
		for (i = 0, iEH = 0; itVR != itVREnd && i < iSize; i++, iEH++, itVR++) {
			(*entity_handles)[iEH] = itVR->pVVert();
		}
		break;
	case eAllFaces: {
		i = 0;
		iEH = 0;
		if (eCurrentDT == eTriFaces) {
			for (; itTriF != itTriFEnd && i < iSize; i++, iEH++, itTriF++) {
				(*entity_handles)[iEH] = *itTriF;
			}
			if (i < iSize) {
				eCurrentDT = eQuadFaces; // Switch to quads
				itQuadF = pES->spQuadF.begin();
				itQuadFEnd = pES->spQuadF.end();
			}
		}
		if (eCurrentDT == eQuadFaces) {
			for (; itQuadF != itQuadFEnd && i < iSize; i++, iEH++, itQuadF++) {
				(*entity_handles)[iEH] = *itQuadF;
			}
		}
	}
		break;
	case eEdgeFaces:
		for (i = 0, iEH = 0; itEdgeF != itEdgeFEnd && i < iSize;
				i++, iEH++, itEdgeF++) {
			(*entity_handles)[iEH] = *itEdgeF;
		}
		break;
	case eTriFaces:
		for (i = 0, iEH = 0; itTriF != itTriFEnd && i < iSize;
				i++, iEH++, itTriF++) {
			(*entity_handles)[iEH] = *itTriF;
		}
		break;
	case eQuadFaces:
		for (i = 0, iEH = 0; itQuadF != itQuadFEnd && i < iSize;
				i++, iEH++, itQuadF++) {
			(*entity_handles)[iEH] = *itQuadF;
		}
		break;
	case eAllCells2D:
		i = 0;
		iEH = 0;
		if (eCurrentDT == eTriCells) {
			for (; itTriC != itTriCEnd && i < iSize; i++, iEH++, itTriC++) {
				(*entity_handles)[iEH] = *itTriC;
			}
			if (i < iSize) {
				eCurrentDT = eQuadCells; // Switch to quads
				itQuadC = pES->spQuadC.begin();
				itQuadCEnd = pES->spQuadC.end();
			}
		}
		if (eCurrentDT == eQuadCells) {
			for (; itQuadC != itQuadCEnd && i < iSize; i++, iEH++, itQuadC++) {
				(*entity_handles)[iEH] = *itQuadC;
			}
		}
		break;
	case eAllCells3D:
		i = 0;
		iEH = 0;
		if (eCurrentDT == eTetCells) {
			for (; itTetC != itTetCEnd && i < iSize; i++, iEH++, itTetC++) {
				(*entity_handles)[iEH] = *itTetC;
			}
			if (i < iSize) {
				eCurrentDT = ePyrCells; // Switch to pyramids
				itPyrC = pES->spPyrC.begin();
				itPyrCEnd = pES->spPyrC.end();
			}
		}
		if (eCurrentDT == ePyrCells) {
			for (; itPyrC != itPyrCEnd && i < iSize; i++, iEH++, itPyrC++) {
				(*entity_handles)[iEH] = *itPyrC;
			}
			if (i < iSize) {
				eCurrentDT = ePrismCells; // Switch to prisms
				itPrismC = pES->spPrismC.begin();
				itPrismCEnd = pES->spPrismC.end();
			}
		}
		if (eCurrentDT == ePrismCells) {
			for (; itPrismC != itPrismCEnd && i < iSize;
					i++, iEH++, itPrismC++) {
				(*entity_handles)[iEH] = *itPrismC;
			}
			if (i < iSize) {
				eCurrentDT = eHexCells; // Switch to hexes
				itHexC = pES->spHexC.begin();
				itHexCEnd = pES->spHexC.end();
			}
		}
		if (eCurrentDT == eHexCells) {
			for (; itHexC != itHexCEnd && i < iSize; i++, iEH++, itHexC++) {
				(*entity_handles)[iEH] = *itHexC;
			}
		}
		break;
	case eEverything: {
		i = 0, iEH = 0;
		switch (eCurrentDT) {
		case eVerts:
			for (; itVR != itVREnd && i < iSize; i++, iEH++, itVR++)
				(*entity_handles)[iEH] = *itVR;
			if (i == iSize)
				break;
			// fallthrough
		case eEdgeFaces:
			for (; itEdgeF != itEdgeFEnd && i < iSize; i++, iEH++, itEdgeF++)
				(*entity_handles)[iEH] = *itEdgeF;
			if (i == iSize)
				break;
			// fallthrough
		case eTriFaces:
			for (; itTriF != itTriFEnd && i < iSize; i++, iEH++, itTriF++)
				(*entity_handles)[iEH] = *itTriF;
			if (i == iSize)
				break;
			// fallthrough
		case eQuadFaces:
			for (; itQuadF != itQuadFEnd && i < iSize; i++, iEH++, itQuadF++)
				(*entity_handles)[iEH] = *itQuadF;
			if (i == iSize)
				break;
			// fallthrough
		case eTriCells:
			for (; itTriC != itTriCEnd && i < iSize; i++, iEH++, itTriC++)
				(*entity_handles)[iEH] = *itTriC;
			if (i == iSize)
				break;
			// fallthrough
		case eQuadCells:
			for (; itQuadC != itQuadCEnd && i < iSize; i++, iEH++, itQuadC++)
				(*entity_handles)[iEH] = *itQuadC;
			if (i == iSize)
				break;
			// fallthrough
		case eTetCells:
			for (; itTetC != itTetCEnd && i < iSize; i++, iEH++, itTetC++)
				(*entity_handles)[iEH] = *itTetC;
			if (i == iSize)
				break;
			// fallthrough
		case ePyrCells:
			for (; itPyrC != itPyrCEnd && i < iSize; i++, iEH++, itPyrC++)
				(*entity_handles)[iEH] = *itPyrC;
			if (i == iSize)
				break;
			// fallthrough
		case ePrismCells:
			for (; itPrismC != itPrismCEnd && i < iSize; i++, iEH++, itPrismC++)
				(*entity_handles)[iEH] = *itPrismC;
			if (i == iSize)
				break;
			// fallthrough
		case eHexCells:
			for (; itHexC != itHexCEnd && i < iSize; i++, iEH++, itHexC++)
				(*entity_handles)[iEH] = *itHexC;
			break; // Nowhere to go from here.
		default:
			// Should never get here.
			break;
		}
		break;
	}
	case eTriCells:
		for (i = 0, iEH = 0; itTriC != itTriCEnd && i < iSize;
				i++, iEH++, itTriC++) {
			(*entity_handles)[iEH] = *itTriC;
		}
		break;
	case eQuadCells:
		for (i = 0, iEH = 0; itQuadC != itQuadCEnd && i < iSize;
				i++, iEH++, itQuadC++) {
			(*entity_handles)[iEH] = *itQuadC;
		}
		break;
	case eTetCells:
		for (i = 0, iEH = 0; itTetC != itTetCEnd && i < iSize;
				i++, iEH++, itTetC++) {
			(*entity_handles)[iEH] = *itTetC;
		}
		break;
	case ePyrCells:
		for (i = 0, iEH = 0; itPyrC != itPyrCEnd && i < iSize;
				i++, iEH++, itPyrC++) {
			(*entity_handles)[iEH] = *itPyrC;
		}
		break;
	case ePrismCells:
		for (i = 0, iEH = 0; itPrismC != itPrismCEnd && i < iSize;
				i++, iEH++, itPrismC++) {
			(*entity_handles)[iEH] = *itPrismC;
		}
		break;
	case eHexCells:
		for (i = 0, iEH = 0; itHexC != itHexCEnd && i < iSize;
				i++, iEH++, itHexC++) {
			(*entity_handles)[iEH] = *itHexC;
		}
		break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
		*entity_handles_size = 0;
		return false;
		break;
	case eUndefined:
		assert(0);
		*entity_handles_size = 0;
		return false;
		break;
	} // Done with switch
	*entity_handles_size = i;
	return (i > 0);
} // qGetNextBlock

void EntitySetWorkset::vDelete(Entity* pEToDelete)

{
	switch (eCurrentDT) {
	case eVerts:
		if (itVR != itVREnd && itVR->pVVert() == pEToDelete)
			itVR++;
		break;
	case eEdgeFaces:
		if (itEdgeF != itEdgeFEnd && (*itEdgeF) == pEToDelete)
			itEdgeF++;
		break;
	case eTriFaces:
		if (itTriF != itTriFEnd && (*itTriF) == pEToDelete)
			itTriF++;
		break;
	case eQuadFaces:
		if (itQuadF != itQuadFEnd && (*itQuadF) == pEToDelete)
			itQuadF++;
		break;
	case eTriCells:
		if (itTriC != itTriCEnd && (*itTriC) == pEToDelete)
			itTriC++;
		break;
	case eQuadCells:
		if (itQuadC != itQuadCEnd && (*itQuadC) == pEToDelete)
			itQuadC++;
		break;
	case eTetCells:
		if (itTetC != itTetCEnd && (*itTetC) == pEToDelete)
			itTetC++;
		break;
	case ePyrCells:
		if (itPyrC != itPyrCEnd && (*itPyrC) == pEToDelete)
			itPyrC++;
		break;
	case ePrismCells:
		if (itPrismC != itPrismCEnd && (*itPrismC) == pEToDelete)
			itPrismC++;
		break;
	case eHexCells:
		if (itHexC != itHexCEnd && (*itHexC) == pEToDelete)
			itHexC++;
		break;
	case eAllFaces:
		if (eCurrentDT == eTriFaces) {
			if (itTriF != itTriFEnd && (*itTriF) == pEToDelete)
				itTriF++;
			// This clause, and its analogs below, shouldn't be necessary,
			// but leave them for now anyway....
			if (itTriF == itTriFEnd) {
				eCurrentDT = eQuadFaces;
				itQuadF = pES->spQuadF.begin();
			}
		} else {
			if (itQuadF != itQuadFEnd && (*itQuadF) == pEToDelete)
				itQuadF++;
		}
		break;
	case eAllCells2D:
		if (eCurrentDT == eTriCells) {
			if (itTriC != itTriCEnd && (*itTriC) == pEToDelete)
				itTriC++;
			if (itTriC == itTriCEnd) {
				eCurrentDT = eQuadCells;
				itQuadC = pES->spQuadC.begin();
			}
		} else {
			if (itQuadC != itQuadCEnd && (*itQuadC) == pEToDelete)
				itQuadC++;
		}
		break;
	case eAllCells3D:
		if (eCurrentDT == eTetCells) {
			if (itTetC != itTetCEnd && (*itTetC) == pEToDelete)
				itTetC++;
			if (itTetC == itTetCEnd) {
				eCurrentDT = ePyrCells;
				itPyrC = pES->spPyrC.begin();
				if (itPyrC == itPyrCEnd) {
					eCurrentDT = ePrismCells;
					itPrismC = pES->spPrismC.begin();
					if (itPrismC == itPrismCEnd) {
						eCurrentDT = eHexCells;
						itHexC = pES->spHexC.begin();
					} // if ( itPrismC...
				} // if ( itPyrC...
			}
		} else if (eCurrentDT == ePyrCells) {
			if (itPyrC != itPyrCEnd && (*itPyrC) == pEToDelete)
				itPyrC++;
			if (itPyrC == itPyrCEnd) {
				eCurrentDT = ePrismCells;
				itPrismC = pES->spPrismC.begin();
				if (itPrismC == itPrismCEnd) {
					eCurrentDT = eHexCells;
					itHexC = pES->spHexC.begin();
				} // if ( itPrismC...
			} // if ( itPyrC...
		} else if (eCurrentDT == ePrismCells) {
			if (itPrismC != itPrismCEnd && (*itPrismC) == pEToDelete)
				itPrismC++;
			if (itPrismC == itPrismCEnd) {
				eCurrentDT = eHexCells;
				itHexC = pES->spHexC.begin();
			} // if ( itPrismC...
		} else if (eCurrentDT == eHexCells) {
			if (itHexC != itHexCEnd && (*itHexC) == pEToDelete)
				itHexC++;
		}
		break;
	case eEverything: {
		// This one, though much simpler than those above, should work
		// just fine.
		eDT = eCurrentDT;
		vDelete(pEToDelete); // This may update eCurrentDT along the way.
		eDT = eEverything;
	}
		break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
		// Do absolutely nothing, as there won't ever be any of these.
		break;
	case eUndefined:
		assert(0);
		break;
	} // Done with switch
} // vDelete

EntityListWorkset::EntityListWorkset(EntityList * const pELIn,
		const int eETypeIn, const int eETopoIn, const int iSzIn) :
		Workset(eETypeIn, eETopoIn, iSzIn, pELIn->iMeshTopologicalDimension()), pEL(
				pELIn), iFound(0), iTotal(0), qContinue(false), iter(), iterEnd() {
	assert(pELIn);
	pELIn->spWorksets.insert(this);
	qReset();
}

EntityListWorkset::~EntityListWorkset() {
	pEL->spWorksets.erase(this);
}

bool EntityListWorkset::qReset()

{
	iFound = 0;
	iter = pEL->lpEntities.begin();
	iterEnd = pEL->lpEntities.end();
	// 'qContinue' is a flag to tell the 'qGetNextBlock' method to
	// increment the STL iterator before getting the next entity
	qContinue = false;

	switch (eDT) {
	case eVerts:
		iTotal = pEL->iNumVerts();
		break;
	case eEdgeFaces:
		iTotal = pEL->iNumEdges();
		break;
	case eTriFaces:
		iTotal = pEL->iNumTriFaces();
		break;
	case eQuadFaces:
		iTotal = pEL->iNumQuadFaces();
		break;
	case eAllFaces:
		iTotal = pEL->iNumTriFaces() + pEL->iNumQuadFaces();
		break;
	case eAllCells2D:
		iTotal = pEL->iNumTriCells() + pEL->iNumQuadCells();
		break;
	case eAllCells3D:
		iTotal = pEL->iNumTetCells() + pEL->iNumPyrCells()
				+ pEL->iNumPrismCells() + pEL->iNumHexCells();
		break;
	case eTriCells:
		iTotal = pEL->iNumTriCells();
		break;
	case eQuadCells:
		iTotal = pEL->iNumQuadCells();
		break;
	case eTetCells:
		iTotal = pEL->iNumTetCells();
		break;
	case ePyrCells:
		iTotal = pEL->iNumPyrCells();
		break;
	case ePrismCells:
		iTotal = pEL->iNumPrismCells();
		break;
	case eHexCells:
		iTotal = pEL->iNumHexCells();
		break;
	case eEverything:
		iTotal = pEL->iNumVerts() + pEL->iNumEdges() + pEL->iNumTriFaces()
				+ pEL->iNumQuadFaces() + pEL->iNumTriCells()
				+ pEL->iNumQuadCells() + pEL->iNumTetCells()
				+ pEL->iNumPyrCells() + pEL->iNumPrismCells()
				+ pEL->iNumHexCells();
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
	assert(iter != iterEnd || iTotal == 0);
	return (iFound < iTotal);
}

void EntityListWorkset::vAddEntsByType(iBase_EntityHandle** entity_handles,
		int iType, int& iCurr) {
	Entity* pEnt;
	while (iCurr < iSize && iFound < iTotal) {
		if (qContinue)
			++iter;
		else
			qContinue = true;
		assert(iter != iterEnd);
		pEnt = *iter;
		if (pEnt->getEntType() == iType) {
			(*entity_handles)[iCurr++] = pEnt;
			iFound++;
		}
	}
}

void EntityListWorkset::vAddEntsByTopo(iBase_EntityHandle** entity_handles,
		int iTopo, int& iCurr) {
	Entity* pEnt;
	while (iCurr < iSize && iFound < iTotal) {
		if (qContinue)
			++iter;
		else
			qContinue = true;
		assert(iter != iterEnd);
		pEnt = *iter;
		if (pEnt->getEntTopology() == iTopo) {
			(*entity_handles)[iCurr++] = pEnt;
			iFound++;
		}
	}
}

bool EntityListWorkset::qGetNextBlock(iBase_EntityHandle** entity_handles,
		int *entity_handles_allocated, int *entity_handles_size)

		{
	int i;
	vEnsureArraySize<iBase_EntityHandle>(entity_handles,
			entity_handles_allocated, iSize);

	switch (eDT) {
	case eVerts:
		i = 0;
		vAddEntsByType(entity_handles, iBase_VERTEX, i);
		break;
	case eEdgeFaces:
		i = 0;
		vAddEntsByType(entity_handles, iBase_EDGE, i);
		break;
	case eAllFaces:
	case eAllCells2D:
		i = 0;
		vAddEntsByType(entity_handles, iBase_FACE, i);
		break;
	case eAllCells3D:
		i = 0;
		vAddEntsByType(entity_handles, iBase_REGION, i);
		break;
	case eTriFaces:
	case eTriCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_TRIANGLE, i);
		break;
	case eQuadFaces:
	case eQuadCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_QUADRILATERAL, i);
		break;
	case eTetCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_TETRAHEDRON, i);
		break;
	case ePyrCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_PYRAMID, i);
		break;
	case ePrismCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_PRISM, i);
		break;
	case eHexCells:
		i = 0;
		vAddEntsByTopo(entity_handles, iMesh_HEXAHEDRON, i);
		break;
	case eEverything: {
		Entity* pEnt;
		i = 0;
		while (i < iSize && iFound < iTotal) {
			if (qContinue)
				++iter;
			else
				qContinue = true;
			assert(iter != iterEnd);
			pEnt = *iter;
			(*entity_handles)[i++] = pEnt;
			iFound++;
		}
	}
		break;
	case ePolygonCells:
	case ePolygonFaces:
	case ePolyhedra:
	case eSeptaCells:
	case eImpossible:
		*entity_handles_size = 0;
		return false;
		break;
	case eUndefined:
		assert(0);
		*entity_handles_size = 0;
		return false;
		break;
	} // Done with switch
	*entity_handles_size = i;
	return (i > 0);
} // qGetNextBlock

void EntityListWorkset::vAdjustIterForDeletionByType(
		std::list<Entity*>::iterator itToDelete, const int iType) {
	std::list<Entity*>::iterator iterTemp, iterPrev;
	if (iter == itToDelete) {
		if (iter == pEL->lpEntities.begin()) {
			if ((*iter)->getEntType() == iType)
				iTotal--;
			++iter;
		} else {
			iTotal--;
			iFound--;
			if (iFound == 0) {
				iter = pEL->lpEntities.begin();
			} else {
				for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter;
						iterTemp++)
					if ((*iterTemp)->getEntType() == iType)
						iterPrev = iterTemp;
				iter = iterPrev;
			}
		}
	} else if ((*itToDelete)->getEntType() == iType) {
		iTotal--;
		iterTemp = pEL->lpEntities.begin();
		while (iterTemp != iter && iterTemp != itToDelete)
			iterTemp++;
		if (iterTemp != iter)
			iFound--;
	}
}

void EntityListWorkset::vAdjustIterForDeletionByTopo(
		std::list<Entity*>::iterator itToDelete, const int iTopo) {
	std::list<Entity*>::iterator iterTemp, iterPrev;
	if (iter == itToDelete) {
		if (iter == pEL->lpEntities.begin()) {
			if ((*iter)->getEntTopology() == iTopo)
				iTotal--;
			++iter;
		} else {
			iTotal--;
			iFound--;
			if (iFound == 0) {
				iter = pEL->lpEntities.begin();
			} else {
				for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter;
						iterTemp++)
					if ((*iterTemp)->getEntTopology() == iTopo)
						iterPrev = iterTemp;
				iter = iterPrev;
			}
		}
	} else if ((*itToDelete)->getEntTopology() == iTopo) {
		iTotal--;
		iterTemp = pEL->lpEntities.begin();
		while (iterTemp != iter && iterTemp != itToDelete)
			iterTemp++;
		if (iterTemp != iter)
			iFound--;
	}
}

void EntityListWorkset::vDelete(std::list<Entity*>::iterator itToDelete)

{
	switch (eDT) {
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
	case eEverything: {
		std::list<Entity*>::iterator iterTemp, iterPrev;
		if (iter == itToDelete) {
			if (iter == pEL->lpEntities.begin()) {
				iTotal--;
				++iter;
			} else {
				iTotal--;
				iFound--;
				if (iFound == 0) {
					iter = pEL->lpEntities.begin();
					qContinue = false;
				} else {
					for (iterTemp = pEL->lpEntities.begin(); iterTemp != iter;
							iterTemp++)
						iterPrev = iterTemp;
					iter = iterPrev;
				}
			}
		} else {
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

void EntityListWorkset::vAdded(Entity* pEAdded)

{
	if (iFound == 0 && iter == iterEnd) {
		iter = pEL->lpEntities.begin();
		qContinue = false;
	}

	switch (eDT) {
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
	case eEverything:
		iTotal++;
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
} // vAdded
} // Close namespace
