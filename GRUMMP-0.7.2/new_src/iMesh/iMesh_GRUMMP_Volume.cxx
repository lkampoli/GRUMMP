#include "iMesh.h"
#include "iMesh_GRUMMP.hh"
using namespace ITAPS_GRUMMP;

#include <iostream>
#include <vector>
#include "GR_Entity.h"
#include "GR_Geometry.h"
#include "GR_VolMesh.h"
#include "iMesh_GRUMMP_misc.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"

// Unless otherwise indicated, these funcs should be 3D-specific, in
// iMesh_Volume

iMesh_Volume::iMesh_Volume(VolMesh *pVM) :
    iMesh_Base(eVolume)
{
  pM = pVM;
  externalData = true;
}

void
iMesh_Volume::vCreateMesh()
{
  VolMesh *pVM = new VolMesh();
  pM = pVM;
}

void
iMesh_Volume::load(const char *name, const char *, int* err)
{
  m_err = *err = iBase_SUCCESS;
  VolMesh *pVM = new VolMesh;
  pVM->readFromFile(name);
  pM = pVM;
}

void
iMesh_Volume::save(const char *name, const char *, int* err)
{
  m_err = *err = iBase_SUCCESS;
  VolMesh *pVM = dynamic_cast<VolMesh*>(pM);
  if (pVM->isSimplicial()) {
    writeNative(*pVM, name);
  }
  else {
    writeVTKLegacy(*pVM, name);
  }
}

void
iMesh_Volume::setAdjTable(int* adjacency_table, int* err)
{
  *err = iBase_SUCCESS;
  if (adjacency_table[1] != iBase_UNAVAILABLE
      || adjacency_table[4] != iBase_UNAVAILABLE
      || adjacency_table[5] != iBase_UNAVAILABLE
      || adjacency_table[6] != iBase_UNAVAILABLE
      || adjacency_table[7] != iBase_UNAVAILABLE
      || adjacency_table[9] != iBase_UNAVAILABLE
      || adjacency_table[13] != iBase_UNAVAILABLE) {
    m_err = *err = iBase_NOT_SUPPORTED;
  }
  getAdjTable(adjacency_table);
}

void
iMesh_Volume::getAdjTable(int* adjacency_table)
{
  // Array is known to be okay before this call.

  // The contents of this adjacency table should be:
  //    V E F R
  //  V I U L L
  //  E U U U U
  //  F I U I I
  //  R I U I I
  // where I = immediate
  //       L = local traversal
  //       U = unavailable

  // The adjacency info available in GRUMMP is always the same.
  adjacency_table[0] = iBase_AVAILABLE;
  adjacency_table[1] = iBase_UNAVAILABLE;
  adjacency_table[2] = iBase_ALL_ORDER_1;
  adjacency_table[3] = iBase_ALL_ORDER_1;

  adjacency_table[4] = iBase_UNAVAILABLE;
  adjacency_table[5] = iBase_UNAVAILABLE;
  adjacency_table[6] = iBase_UNAVAILABLE;
  adjacency_table[7] = iBase_UNAVAILABLE;

  adjacency_table[8] = iBase_ALL_ORDER_1;
  adjacency_table[9] = iBase_UNAVAILABLE;
  adjacency_table[10] = iBase_AVAILABLE;
  adjacency_table[11] = iBase_ALL_ORDER_1;

  adjacency_table[12] = iBase_ALL_ORDER_1;
  adjacency_table[13] = iBase_UNAVAILABLE;
  adjacency_table[14] = iBase_ALL_ORDER_1;
  adjacency_table[15] = iBase_AVAILABLE;
}

void
iMesh_Volume::getEntities(const int entity_type, const int entity_topology,
			  iBase_EntityHandle** entity_handles,
			  int* entity_handles_allocated,
			  int* entity_handles_size, int *err)

{
  *err = iBase_SUCCESS;

  VolMesh *pVM = dynamic_cast<VolMesh*>(pM);

  if (entity_topology == iMesh_POINT
      || (entity_topology == iMesh_ALL_TOPOLOGIES && entity_type == iBase_VERTEX)) {
    ////////////////////////////////////////////////
    // Returning all vertices.
    ////////////////////////////////////////////////
    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, getNumVertices(),
		   iBase_EntityHandle);

    unsigned int iNV = pM->getNumVerts(); // This is >= iNVert
    unsigned int iNOrigV = iNV;
    unsigned int iV, iEH;
    for (iV = 0, iEH = 0; iV < iNOrigV; iV++) {
      Vert *pV = pM->getVert(iV);
      if (pV->isDeleted()) {
	iNV--;
      }
      else {
	(*entity_handles)[iEH] = pV;
	iEH++;
      }
    }
    assert(iEH == iNV);
    assert(iNV <= iNOrigV);
    *entity_handles_size = iEH;
  }
  else if (entity_topology == iMesh_POLYGON) {
    ////////////////////////////////////////////////
    // Returning all polygonal faces, of which there are zero.
    ////////////////////////////////////////////////
    int iNPoly = 0;
    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNPoly,
		   iBase_EntityHandle);
    *entity_handles_size = 0;
  }
  else if (entity_topology == iMesh_TRIANGLE) {
    ////////////////////////////////////////////////
    // Returning all triangular faces.
    ////////////////////////////////////////////////
    int iNTriTemp = 0;
    iNTriTemp = pVM->getNumTriFaces();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNTriTemp,
		   iBase_EntityHandle);

    int iTri, iEH;
    for (iTri = 0, iEH = 0; iTri < iNTriTemp; iTri++) {
      Face *pF = pM->getFace(iTri);
      if (pF->isDeleted()) {
	iNTriTemp--;
      }
      else {
	(*entity_handles)[iEH] = pF;
	iEH++;
      }
    }
    *entity_handles_size = iNTriTemp;
    assert(iEH == *entity_handles_size);
  } // Done with triangular faces.
  else if (entity_topology == iMesh_QUADRILATERAL) {
    ////////////////////////////////////////////////
    // Returning all quadrilateral faces.
    ////////////////////////////////////////////////
    int iNQuadTemp = 0, iNTriTemp = 0;
    iNTriTemp = pVM->getNumTriFaces();
    iNQuadTemp = pVM->getNumQuadFaces();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNQuadTemp,
		   iBase_EntityHandle);

    int iQuad, iEH;
    for (iQuad = iNTriTemp, iEH = 0; iQuad < iNQuadTemp + iNTriTemp; iQuad++) {
      Face *pF = pVM->getFace(iQuad);
      if (pF->isDeleted()) {
	iNQuadTemp--;
      }
      else {
	(*entity_handles)[iEH] = pF;
	iEH++;
      }
    }
    *entity_handles_size = iNQuadTemp;
    assert(iEH == *entity_handles_size);
  } // Done with quadrilateral cells.
  else if (entity_topology == iMesh_ALL_TOPOLOGIES
      && entity_type == iBase_FACE) {
    ////////////////////////////////////////////////
    // Returning all faces.
    ////////////////////////////////////////////////
    int iNFace = 0;
    iNFace = pVM->getNumFaces();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNFace,
		   iBase_EntityHandle);

    int iFace, iEH;
    for (iFace = 0, iEH = 0; iFace < iNFace; iFace++) {
      Face *pF = pVM->getFace(iFace);
      if (pF->isDeleted()) {
	iNFace--;
      }
      else {
	(*entity_handles)[iEH] = pF;
	iEH++;
      }
    }
    *entity_handles_size = iNFace;
    assert(iEH == *entity_handles_size);
  } // Done with all faces.
  else if (entity_topology == iMesh_TETRAHEDRON) {
    ////////////////////////////////////////////////
    // Returning all tets.
    ////////////////////////////////////////////////
    int iNTetTemp = pVM->getNumTetCells();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNTetTemp,
		   iBase_EntityHandle);

    int iTet, iEH;
    for (iTet = 0, iEH = 0; iTet < iNTetTemp; iTet++) {
      Cell *pC = pVM->getCell(iTet);
      if (pC->isDeleted()) {
	iNTetTemp--;
      }
      else {
	(*entity_handles)[iEH] = pC;
	iEH++;
      }
    }
    *entity_handles_size = iNTetTemp;
    assert(iEH == *entity_handles_size);
  }
  else if (entity_topology == iMesh_PYRAMID) {
    ////////////////////////////////////////////////
    // Returning all pyramids.
    ////////////////////////////////////////////////
    int iNPyrTemp = pVM->getNumPyrCells();
    int iBegin = pVM->getNumTetCells();
    int iEnd = iBegin + iNPyrTemp;

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNPyrTemp,
		   iBase_EntityHandle);

    int iPyr, iEH;
    for (iPyr = iBegin, iEH = 0; iPyr < iEnd; iPyr++) {
      Cell *pC = pVM->getCell(iPyr);
      if (pC->isDeleted()) {
	iNPyrTemp--;
      }
      else {
	(*entity_handles)[iEH] = pC;
	iEH++;
      }
    }
    *entity_handles_size = iNPyrTemp;
    assert(iEH == *entity_handles_size);
  }
  else if (entity_topology == iMesh_PRISM) {
    ////////////////////////////////////////////////
    // Returning all prisms.
    ////////////////////////////////////////////////
    int iNPrismTemp = pVM->getNumPrismCells();
    int iBegin = pVM->getNumTetCells() + pVM->getNumPyrCells();
    int iEnd = iBegin + iNPrismTemp;

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNPrismTemp,
		   iBase_EntityHandle);

    int iPrism, iEH;
    for (iPrism = iBegin, iEH = 0; iPrism < iEnd; iPrism++) {
      Cell *pC = pVM->getCell(iPrism);
      if (pC->isDeleted()) {
	iNPrismTemp--;
      }
      else {
	(*entity_handles)[iEH] = pC;
	iEH++;
      }
    }
    *entity_handles_size = iNPrismTemp;
    assert(iEH == *entity_handles_size);
  }
  else if (entity_topology == iMesh_HEXAHEDRON) {
    ////////////////////////////////////////////////
    // Returning all hexs.
    ////////////////////////////////////////////////
    int iNHexTemp = pVM->getNumHexCells();
    int iBegin = pVM->getNumTetCells() + pVM->getNumPyrCells()
	+ pVM->getNumPrismCells();
    int iEnd = iBegin + iNHexTemp;

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNHexTemp,
		   iBase_EntityHandle);

    int iHex, iEH;
    for (iHex = iBegin, iEH = 0; iHex < iEnd; iHex++) {
      Cell *pC = pVM->getCell(iHex);
      if (pC->isDeleted()) {
	iNHexTemp--;
      }
      else {
	(*entity_handles)[iEH] = pC;
	iEH++;
      }
    }
    *entity_handles_size = iNHexTemp;
    assert(iEH == *entity_handles_size);
  }
  else if (entity_topology == iMesh_ALL_TOPOLOGIES
      && entity_type == iBase_REGION) {
    ////////////////////////////////////////////////
    // Returning all cells!
    ////////////////////////////////////////////////
    int iNCell = pVM->getNumCells();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNCell,
		   iBase_EntityHandle);

    int iCell, iEH;
    for (iCell = 0, iEH = 0; iCell < iNCell; iCell++) {
      Cell *pC = pVM->getCell(iCell);
      if (pC->isDeleted()) {
	iNCell--;
      }
      else {
	(*entity_handles)[iEH] = pC;
	iEH++;
      }
    }
    *entity_handles_size = iNCell;
    assert(iEH == *entity_handles_size);
  } // Done with all cells.
  else if (entity_topology == iMesh_ALL_TOPOLOGIES
      && entity_type == iBase_ALL_TYPES) {
    ////////////////////////////////////////////////
    // Returning all entities!
    ////////////////////////////////////////////////

    vCopyRootSet(entity_handles, entity_handles_allocated, entity_handles_size);
    if (*err != iBase_SUCCESS)
      return;
  }
  else if (entity_topology == iMesh_POLYHEDRON
      || entity_topology == iMesh_SEPTAHEDRON
      || entity_topology == iMesh_LINE_SEGMENT
      || (entity_topology == iMesh_ALL_TOPOLOGIES && entity_type == iBase_EDGE)) {
    ////////////////////////////////////////////////
    // There aren't any of these.
    ////////////////////////////////////////////////
    int iN = 0;
    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iN,
		   iBase_EntityHandle);
    // Nothing else to do, except return iN when I can.
    *entity_handles_size = 0;
  }
}

void
iMesh_Volume::createEntArr(/*in*/const int new_entity_topology,
/*in*/const iBase_EntityHandle* lower_order_entity_handles,
			   /*in*/const int lower_order_entity_handles_size,
			   /*out*/iBase_EntityHandle** new_entity_handles,
			   /*out*/int* new_entity_handles_allocated,
			   /*out*/int* new_entity_handles_size,
			   /*inout*/int** status,
			   /*inout*/int* status_allocated,
			   /*out*/int* status_size, int *err)
{
  *err = iBase_SUCCESS;

  *new_entity_handles_size = 0;
  *status_size = 0;
  int vertsPerEnt;
  switch (new_entity_topology)
    {
    case iMesh_POINT:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      return;
      break;
    case iMesh_TRIANGLE:
      vertsPerEnt = 3;
      break;
    case iMesh_QUADRILATERAL:
      vertsPerEnt = 4;
      break;
    case iMesh_TETRAHEDRON:
      vertsPerEnt = 4;
      break;
    case iMesh_PRISM:
      vertsPerEnt = 6;
      break;
    case iMesh_PYRAMID:
      vertsPerEnt = 5;
      break;
    case iMesh_HEXAHEDRON:
      vertsPerEnt = 8;
      break;
    case iMesh_SEPTAHEDRON:
    case iMesh_POLYHEDRON:
    case iMesh_LINE_SEGMENT:
    case iMesh_POLYGON:
      m_err = *err = iBase_NOT_SUPPORTED;
      break;
    case iMesh_ALL_TOPOLOGIES:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      break;
    } // End switch over topologies

  // Now set up arrays and call createEnt to do the actual work
  if (lower_order_entity_handles_size % vertsPerEnt != 0) {
    m_err = *err = iBase_INVALID_ENTITY_COUNT;
    return;
  }

  (*status_size) = (*new_entity_handles_size) = lower_order_entity_handles_size
      / vertsPerEnt;
  TRY_ARRAY_SIZE(new_entity_handles, new_entity_handles_allocated,
		 *new_entity_handles_size, iBase_EntityHandle);
  TRY_ARRAY_SIZE(status, status_allocated, *status_size, int);

  if (*err != iBase_SUCCESS)
    return;

  for (int iLO = 0, iEH = 0, iStat = 0;
      iLO <= lower_order_entity_handles_size - 1; iEH++, iStat++, iLO +=
	  vertsPerEnt) {
    createEnt(new_entity_topology, lower_order_entity_handles + iLO,
	      vertsPerEnt, (*new_entity_handles) + iEH, (*status) + iStat, err);
  }
}

void
iMesh_Volume::createEnt(/*in*/const int new_entity_topology,
/*in*/const iBase_EntityHandle* lower_order_entity_handles,
			/*in*/const int lower_order_entity_handles_size,
			/*out*/iBase_EntityHandle* new_entity_handle,
			/*out*/int* status, int *err)
{
  VolMesh *pVM = dynamic_cast<VolMesh*>(pM);

  switch (new_entity_topology)
    {
    case iMesh_POINT:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      *new_entity_handle = NULL;
      *status = iBase_CREATION_FAILED;
      break;
    case iMesh_TRIANGLE:
      {
	// We're going to create a new TriFace!  That means three
	// vertices.
	if (lower_order_entity_handles_size != 3) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  return;
	}

	Vert *pV0 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	Vert *pV1 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	Vert *pV2 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));

	if (!pV0->isValid() || !pV1->isValid() || !pV2->isValid()) {
	  m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	  *new_entity_handle = NULL;
	  *status = iBase_CREATION_FAILED;
	  break;
	}

	bool qExistedAlready = false;
	Face* pF = pVM->createFace(qExistedAlready, pV0, pV1, pV2);
	if (!pF->isValid()) {
	  m_err = *err = iBase_ENTITY_CREATION_ERROR;
	  *new_entity_handle = NULL;
	  *status = iBase_CREATION_FAILED;
	  break;
	}
	assert(pF->isValid());
	if (qExistedAlready) {
	  // Return this face as a handle, with status ALREADY_EXISTED
	  *new_entity_handle = pF;
	  *status = iBase_ALREADY_EXISTED;
	}
	else {
	  // Created it.
	  *new_entity_handle = pF;
	  *status = iBase_NEW;
	}
      }
      break;
    case iMesh_QUADRILATERAL:
      {
	// We're going to create a new QuadFace!  That means four
	// vertices
	if (lower_order_entity_handles_size != 4) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  break;
	}

	Vert *pV0 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	Vert *pV1 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	Vert *pV2 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	Vert *pV3 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[3]));

	if (!pV0->isValid() || !pV1->isValid() || !pV2->isValid()
	    || !pV3->isValid()) {
	  m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	  *new_entity_handle = NULL;
	  *status = iBase_CREATION_FAILED;
	  break;
	}

	bool qExistedAlready = false;
	Face* pF = pVM->createFace(qExistedAlready, pV0, pV1, pV2, pV3);
	if (!pF->isValid()) {
	  m_err = *err = iBase_ENTITY_CREATION_ERROR;
	  *new_entity_handle = NULL;
	  *status = iBase_CREATION_FAILED;
	  break;
	}
	assert(pF->isValid());
	if (qExistedAlready) {
	  // Return this face as a handle, with status ALREADY_EXISTED
	  *new_entity_handle = pF;
	  *status = iBase_ALREADY_EXISTED;
	}
	else {
	  // Created it.
	  *new_entity_handle = pF;
	  *status = iBase_NEW;
	}
      }
      break;
    case iMesh_TETRAHEDRON:
      {
	// New tets can be built from faces or vertices.
	if (lower_order_entity_handles_size != 4) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  break;
	}

	Face *apF[4];
	apF[0] =
	    dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	apF[1] =
	    dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	apF[2] =
	    dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	apF[3] =
	    dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[3]));

	if (apF[0]->isValid() && apF[1]->isValid() && apF[2]->isValid()
	    && apF[3]->isValid()) {
	  bool qExist;
	  Cell *pC = pVM->createTetCell(qExist, apF[0], apF[1], apF[2], apF[3]);
	  if (pC->isValid()) {
	    m_err = *err = iBase_SUCCESS;
	    *new_entity_handle = pC;
	    *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	  }
	  else {
	    m_err = *err = iBase_ENTITY_CREATION_ERROR;
	    *new_entity_handle = NULL;
	    *status = iBase_CREATION_FAILED;
	    break;
	  }
	} // Created tet from faces
	else {
	  // Had better all be verts!
	  Vert *apV[4];
	  apV[0] =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	  apV[1] =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  apV[2] =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	  apV[3] =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[3]));

	  if (!apV[0]->isValid() || !apV[1]->isValid() || !apV[2]->isValid()
	      || !apV[3]->isValid()) {
	    m_err = *err = iBase_INVALID_ENTITY_TYPE;
	    *new_entity_handle = NULL;
	    *status = iBase_CREATION_FAILED;
	  }
	  else {
	    bool qExist;
	    Cell *pC = pVM->createTetCell(qExist, apV[0], apV[1], apV[2],
					  apV[3]);
	    if (pC->isValid()) {
	      m_err = *err = iBase_SUCCESS;
	      *new_entity_handle = pC;
	      *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	    }
	    else {
	      m_err = *err = iBase_ENTITY_CREATION_ERROR;
	      *new_entity_handle = NULL;
	      *status = iBase_CREATION_FAILED;
	      break;
	    }
	  }
	} // Created tet from verts.
      } // Done creating tets
      break;
    case iMesh_PRISM:
      {
	switch (lower_order_entity_handles_size)
	  {
	  case 5:
	    {
	      Face *pF0 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	      Face *pF1 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	      Face *pF2 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	      Face *pF3 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	      Face *pF4 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	      assert(pF0 && pF1 && pF2 && pF3 && pF4);
	      bool qExist;
	      *new_entity_handle = pVM->createPrismCell(qExist, pF0, pF1, pF2,
							pF3, pF4);
	      *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	    }
	    break;
	  case 6:
	    {
	      Vert *pV0 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	      Vert *pV1 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	      Vert *pV2 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	      Vert *pV3 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	      Vert *pV4 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	      Vert *pV5 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[5]));
	      assert(pV0 && pV1 && pV2 && pV3 && pV4 && pV5);
	      bool qExist;
	      *new_entity_handle = pVM->createPrismCell(qExist, pV0, pV1, pV2,
							pV3, pV4, pV5);
	      *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	    }
	    break;
	  default:
	    m_err = *err = iBase_INVALID_ENTITY_COUNT;
	    break;
	  } // end switch on entity count
      }
      break;
    case iMesh_PYRAMID:
      // This stuff will work once the GRUMMP mesh mod overhaul is merged
      // into the trunk.
      {
	if (lower_order_entity_handles_size != 5) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  break;
	}
	Vert *pV0 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	if (pV0 != NULL) {
	  Vert *pV1 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  Vert *pV2 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	  Vert *pV3 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	  Vert *pV4 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	  assert(pV0 && pV1 && pV2 && pV3 && pV4);
	  bool qExist;
	  *new_entity_handle = pVM->createPyrCell(qExist, pV0, pV1, pV2, pV3,
						  pV4);
	  *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	}
	else {
	  Face *pF0 =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	  Face *pF1 =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  Face *pF2 =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	  Face *pF3 =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	  Face *pF4 =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	  assert(pF0 && pF1 && pF2 && pF3 && pF4);
	  bool qExist;
	  *new_entity_handle = pVM->createPyrCell(qExist, pF0, pF1, pF2, pF3,
						  pF4);
	  *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	}
      }
      break;
    case iMesh_HEXAHEDRON:
      // This stuff will work once the GRUMMP mesh mod overhaul is merged
      // into the trunk.
      {
	switch (lower_order_entity_handles_size)
	  {
	  default:
	    m_err = *err = iBase_INVALID_ENTITY_COUNT;
	    break;
	  case 6:
	    {
	      // From faces
	      Face *pF0 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	      Face *pF1 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	      Face *pF2 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	      Face *pF3 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	      Face *pF4 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	      Face *pF5 =
		  dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[5]));
	      assert(pF0 && pF1 && pF2 && pF3 && pF4 && pF5);
	      bool qExist;
	      *new_entity_handle = pVM->createHexCell(qExist, pF0, pF1, pF2,
						      pF3, pF4, pF5);
	    }
	    break;
	  case 8:
	    {
	      // From verts
	      Vert *pV0 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	      Vert *pV1 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	      Vert *pV2 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	      Vert *pV3 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[3]));
	      Vert *pV4 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[4]));
	      Vert *pV5 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[5]));
	      Vert *pV6 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[6]));
	      Vert *pV7 =
		  dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[7]));
	      assert(pV0 && pV1 && pV2 && pV3 && pV4 && pV5 && pV6 && pV7);
	      bool qExist;
	      *new_entity_handle = pVM->createHexCell(qExist, pV0, pV1, pV2,
						      pV3, pV4, pV5, pV6, pV7);
	      *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	    }
	  } // Done switching on number of LO entities
      }
      break;
    case iMesh_SEPTAHEDRON:
    case iMesh_POLYHEDRON:
    case iMesh_LINE_SEGMENT:
    case iMesh_POLYGON:
      m_err = *err = iBase_NOT_SUPPORTED;
      break;
    case iMesh_ALL_TOPOLOGIES:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      break;
    }
}

void
iMesh_Volume::deleteEnt(iBase_EntityHandle entity_handle, int *err)
{
  *err = iBase_SUCCESS;

  Entity *pE = static_cast<Entity*>(entity_handle);
  switch (pE->getEntTopology())
    {
    case iMesh_TETRAHEDRON:
    case iMesh_PYRAMID:
    case iMesh_PRISM:
    case iMesh_HEXAHEDRON:
      {
	// In this case, there are no upward adjacencies, so it's always
	// safe to get rid of these.
	vRemoveAllReferencesTo(pE);
	// Remove connections from its edges.
	Cell *pC = static_cast<Cell*>(pE);
	pM->deleteCell(pC);
	break;
      }
    case iMesh_QUADRILATERAL:
    case iMesh_TRIANGLE:
      // These could have cells that depend on them.
      {
	Face *pF = dynamic_cast<Face*>(pE);
	assert(pF);
	bool qLeftOkay = !pF->getLeftCell()->isValid()
	    || pF->getLeftCell()->isDeleted();
	bool qRightOkay = !pF->getRightCell()->isValid()
	    || pF->getRightCell()->isDeleted();
	if (pF->isBdryFace()) {
	  BFace *pBF = pBFInvalidBFace;
	  // Identify the bdry face, if it's still attached.
	  if (qLeftOkay) {
	    // Only chance is on the right; if it's a BFace, the dynamic
	    // cast will put it in pBF, otherwise NULL.
	    pBF = dynamic_cast<BFace*>(pF->getRightCell());
	  }
	  else if (qRightOkay) {
	    // Same thing, other side.
	    pBF = dynamic_cast<BFace*>(pF->getLeftCell());
	  }
	  if (pBF) {
	    pM->deleteBFace(pBF);
	  }
	  qLeftOkay = qRightOkay = true;
	}

	if (qLeftOkay && qRightOkay) {
	  vRemoveAllReferencesTo(pE);
	  pM->deleteFace(pF);
	}
	else {
	  m_err = *err = iBase_FAILURE;
	}
	break;
      }
    case iMesh_POINT:
      // Every vertex that's connected to anything must have a valid
      // hint face.
      {
	Vert *pV = dynamic_cast<Vert*>(pE);
	assert(pV);
	if (!pV->isDeleted()) {
	  Face *pF = pV->getHintFace();
	  if (!pF->isValid() || !pF->hasVert(pV) || pF->isDeleted()) {
	    vRemoveAllReferencesTo(pE);
	    pM->deleteVert(pV);
	  }
	  else {
	    m_err = *err = iBase_FAILURE;
	  }
	}
	break;
      }
    case iMesh_LINE_SEGMENT:
      // This one will be supported eventually, but not yet.
    case iMesh_SEPTAHEDRON:
    case iMesh_POLYHEDRON:
    case iMesh_POLYGON:
      // These are never expected to be supported by GRUMMP.
      m_err = *err = iBase_NOT_SUPPORTED;
      break;
    case iMesh_ALL_TOPOLOGIES:
      // This one is illegal.
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      break;
    }
}

struct qNotHasVert : public std::binary_function<Face*, Vert*, bool> {
  bool
  operator()(const Face* pF, const Vert* pV) const
  {
    return !pF->hasVert(pV);
  }
};

struct qHasVert : public std::binary_function<Face*, Vert*, bool> {
  bool
  operator()(const Face* pF, const Vert* pV) const
  {
    return pF->hasVert(pV);
  }
};

void
iMesh_Volume::getEntAdj(const iBase_EntityHandle entity_handle,
			const int entity_type_requested,
			iBase_EntityHandle** adj_entity_handles,
			int* adj_entity_handles_allocated,
			int* adj_entity_handles_size, int *err)
{
  Entity *pE = static_cast<Entity*>(entity_handle);
  int eET = pE->getEntType();
  if (eET == entity_type_requested || iBase_EDGE == entity_type_requested) {
    // Not supposed to return any data in this case.
    *adj_entity_handles_size = 0;
    TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, 0,
		   iBase_EntityHandle);
  }
  else if (eET == iBase_VERTEX && entity_type_requested != iBase_VERTEX) {
    Vert *pV = dynamic_cast<Vert*>(pE);
    assert(pV->isValid());

    if (entity_type_requested == iBase_FACE) {
      int iCount = pV->getNumFaces();
      *adj_entity_handles_size = iCount;
      TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		     *adj_entity_handles_size, iBase_EntityHandle);
      for (int i = 0; i < iCount; i++) {
	(*adj_entity_handles)[i] = pV->getFace(i);
      }
    }
    else if (entity_type_requested == iBase_REGION) {
      std::set<Cell*> spC;
      int iNFaces = pV->getNumFaces();
      for (int i = 0; i < iNFaces; i++) {
	Face *pF = pV->getFace(i);
	Cell *pC = pF->getLeftCell();
	if (pC->isValid() && !pC->isBdryCell())
	  spC.insert(pC);
	pC = pF->getRightCell();
	if (pC->isValid() && !pC->isBdryCell())
	  spC.insert(pC);
      }
      *adj_entity_handles_size = spC.size();
      TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		     *adj_entity_handles_size, iBase_EntityHandle);
      std::copy(spC.begin(), spC.end(), *adj_entity_handles);
    }
    else {
      assert(entity_type_requested == iBase_ALL_TYPES);
      std::set<Cell*> spC;
      int iNFaces = pV->getNumFaces();
      Face **apF = new Face*[iNFaces];
      for (int i = 0; i < iNFaces; i++) {
	Face *pF = apF[i] = pV->getFace(i);
	Cell *pC = pF->getLeftCell();
	if (pC->isValid() && !pC->isBdryCell())
	  spC.insert(pC);
	pC = pF->getRightCell();
	if (pC->isValid() && !pC->isBdryCell())
	  spC.insert(pC);
      }

      *adj_entity_handles_size = spC.size() + iNFaces;
      TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		     *adj_entity_handles_size, iBase_EntityHandle);
      iBase_EntityHandle* dummy = std::copy(apF, apF + iNFaces,
					    *adj_entity_handles);
      std::copy(spC.begin(), spC.end(), dummy);
      delete[] apF;
    }
  }
  else {
    switch (entity_type_requested)
      {
      case iBase_VERTEX:
	{
	  *adj_entity_handles_size = pE->getNumVerts();
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  pE->getAllVertHandles(*adj_entity_handles);
	  break;
	}
      case iBase_FACE:
	{
	  *adj_entity_handles_size = pE->getNumFaces();
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  pE->getAllFaceHandles(*adj_entity_handles);
	  break;
	}
      case iBase_REGION:
	{
	  *adj_entity_handles_size = pE->getNumCells();
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  pE->getAllCellHandles(*adj_entity_handles);
	  if (*adj_entity_handles[0] == NULL) {
	    *adj_entity_handles[0] = *adj_entity_handles[1];
	    *adj_entity_handles_size = *adj_entity_handles_size - 1;
	  }
	  break;
	}
      case iBase_ALL_TYPES:
	{
	  *adj_entity_handles_size = pE->getNumCells() + pE->getNumVerts()
	      + pE->getNumFaces();
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  iBase_EntityHandle *tmp = *adj_entity_handles;
	  pE->getAllVertHandles(tmp);
	  tmp += pE->getNumVerts();
	  pE->getAllFaceHandles(tmp);
	  tmp += pE->getNumFaces();
	  pE->getAllCellHandles(tmp);
	  tmp += pE->getNumCells();
	  assert(tmp - *adj_entity_handles == *adj_entity_handles_size);
	  break;
	}
      case iBase_EDGE:
	{
	  // Do nothing, because there are no edges in the 3D mesh.
	  // In fact, applications with any sense have multiple ways to
	  // decide not to get here, but don't throw an error anyway,
	  // because the spec implies we shouldn't.
	  *adj_entity_handles_size = 0;
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  break;
	}
      default:
	// Can't get here, because you'll have already thrown an exception.
	assert(0);
	break;
      }
  } // The general case (not vertex nor requesting type of current entity
}

void
iMesh_Volume::getEntArrAdj(const iBase_EntityHandle* entity_handles,
			   const int entity_handles_size,
			   const int entity_type_requested,
			   iBase_EntityHandle** adj_entity_handles,
			   int* adj_entity_handles_allocated,
			   int* adj_entity_handles_size, int** offset,
			   int* offset_allocated, int* offset_size, int *err)
{
  // Even though it may be a bit awkward, put the switch on requested type
  // on the -outside-, with separate loops inside each case statement.  This
  // will make the code longer, but presumably faster...
  *err = iBase_SUCCESS;

  *offset_size = entity_handles_size + 1;
  TRY_ARRAY_SIZE(offset, offset_allocated, *offset_size, int);

  switch (entity_type_requested)
    {
    case iBase_VERTEX:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET != entity_type_requested) {
	    iEOut += pE->getNumVerts();
	  }
	}

	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET != entity_type_requested) {
	    pE->getAllVertHandles(*adj_entity_handles + (*offset)[iE]);
	  }
	} // The general case (not vertex nor requesting type of current entity
	break;
      }
    case iBase_EDGE:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  // Since GRUMMP currently stores no edges, there's nothing to do here.
	}
	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  // Do nothing, because there are no edges represented in the 3D mesh.
	  // In fact, applications with any sense have multiple ways to
	  // decide not to get here, but don't throw an error anyway,
	  // because the spec implies we shouldn't.
	}
	break;
      }
    case iBase_FACE:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else {
	    iEOut += pE->getNumFaces();
	  }
	}

	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX
	      && entity_type_requested != iBase_VERTEX) {
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    int iStart = (*offset)[iE];
	    int iCount = (*offset)[iE + 1] - iStart;
	    for (int i = 0; i < iCount; i++) {
	      (*adj_entity_handles)[i + iStart] = pV->getFace(i);
	    }
	  }
	  else {
	    pE->getAllFaceHandles(*adj_entity_handles + (*offset)[iE]);
	  }
	}
	break;
      }

    case iBase_REGION:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX) {
	    std::set<Cell*> spC;
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid()); // Should never fail, otherwise why are we here?
	    int iNFaces = pV->getNumFaces();
	    for (int i = 0; i < iNFaces; i++) {
	      Face *pF = pV->getFace(i);
	      Cell *pC = pF->getLeftCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	      pC = pF->getRightCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	    }
	    iEOut += spC.size();
	  }
	  else {
	    iEOut += pE->getNumCells();
	  }
	}
	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	iBase_EntityHandle* pEH = *adj_entity_handles;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX) {
	    std::set<Cell*> spC;
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid()); // Should never fail, otherwise why are we here?
	    int iNFaces = pV->getNumFaces();
	    for (int i = 0; i < iNFaces; i++) {
	      Face *pF = pV->getFace(i);
	      Cell *pC = pF->getLeftCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	      pC = pF->getRightCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	    }
	    pEH = std::copy(spC.begin(), spC.end(), pEH);
	    assert(pEH == (*adj_entity_handles) + (*offset)[iE + 1]);
	  }
	  else {
	    pE->getAllCellHandles(*adj_entity_handles + (*offset)[iE]);
	  }
	}
	break;
      }
    case iBase_ALL_TYPES:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == iBase_VERTEX) {
	    std::set<Cell*> spC;
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid()); // Should never fail, otherwise why are we here?
	    int iNFaces = pV->getNumFaces();
	    for (int i = 0; i < iNFaces; i++) {
	      Face *pF = pV->getFace(i);
	      Cell *pC = pF->getLeftCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	      pC = pF->getRightCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	    }
	    iEOut += spC.size() + iNFaces;
	  }
	  else {
	    iEOut += pE->getNumCells() + pE->getNumFaces() + pE->getNumVerts();
	  }
	}
	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	iBase_EntityHandle* pEH = *adj_entity_handles;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX) {
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    std::set<Cell*> spC;
	    int iNFaces = pV->getNumFaces();
	    Face **apF = new Face*[iNFaces];
	    for (int i = 0; i < iNFaces; i++) {
	      Face *pF = apF[i] = pV->getFace(i);
	      Cell *pC = pF->getLeftCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	      pC = pF->getRightCell();
	      if (pC->isValid() && !pC->isBdryCell())
		spC.insert(pC);
	    }
	    iEOut += spC.size() + iNFaces;
	    pEH = std::copy(apF, apF + iNFaces, pEH);
	    pEH = std::copy(spC.begin(), spC.end(), pEH);
	    delete[] apF;
	  }
	  else {
	    pE->getAllCellHandles(pEH);
	    pEH += pE->getNumCells();
	    pE->getAllFaceHandles(pEH);
	    pEH += pE->getNumFaces();
	    pE->getAllVertHandles(pEH);
	    pEH += pE->getNumVerts();
	  }
	}
	assert(pEH == *adj_entity_handles + *adj_entity_handles_size);
	break;
      }
    } // End of the switch
} // End of iMesh_getEntArrAdj

void
iMesh_Volume::getEnt2ndAdj(const iBase_EntityHandle entity_handle,
			   const int bridge_entity_type,
			   const int requested_entity_type,
			   iBase_EntityHandle** adjacent_entities,
			   int* adjacent_entities_allocated,
			   int* adjacent_entities_size, int* err)
{
  *err = iBase_SUCCESS;
  Entity *pE = static_cast<Entity*>(entity_handle);
  int input_type = pE->getEntType();

  // If the requested type is EDGE, there's no data to be returned.
  if (requested_entity_type == iBase_EDGE) {
    *adjacent_entities_size = 0;
    return;
  }

  // Look, the key here is pretty simple:  identify short-circuit cases
  // that can be handled by a single ITAPS adjacency request; this
  // includes things starting from a vertex, for which GRUMMP has a
  // built-in call that makes a second adjacency call always
  // unnecessary.  The short-circuit cases (in addition to vertex) are
  // those in which the bridge type matches one of the other two and
  // isn't iBase_ALL_TYPES; these are handled by a first-adjacency call
  // to requested_type.  Otherwise, call getEntAdj, then getEntArrAdj
  // and call it a day.

  if (input_type == iBase_VERTEX) {
    // Huge possibilities for optimization here compared with two
    // adjacency calls, because the second one is never needed.  The
    // vertex adjacency call would require a call to vNeighborhood
    // anyway, and that contains all the info we need.

    //  ->all->all   = everything from neighborhood calc
    //  ->all->regi  = 1st adj regions
    //  ->all->face  = all faces from neighborhood calc
    //  ->all->edge  = nothing
    //  ->all->vert  = all verts from neighborhood calc

    //  ->regi->all  = everything from neighborhood calc except regions
    //  ->regi->regi = nothing
    //  ->regi->face = all faces from neighborhood calc
    //  ->regi->edge = nothing
    //  ->regi->vert = all verts from neighborhood calc

    //  ->face->all  = all of neighborhood calc except opposite faces 
    //  ->face->regi = 1st adj regions
    //  ->face->face = nothing
    //  ->face->edge = nothing
    //  ->face->vert = most verts from neighborhood calc

    //  ->edge->all  = all of neighborhood calc except opposite faces
    //                    and verts not connected by an edge
    //  ->edge->regi = 1st adj regions
    //  ->edge->face = 1st adj faces
    //  ->edge->edge = nothing
    //  ->edge->vert = most verts from neighborhood calc

    //  ->vert->all  = nothing
    //  ->vert->regi = nothing
    //  ->vert->face = nothing
    //  ->vert->edge = nothing
    //  ->vert->vert = nothing

    if (bridge_entity_type == iBase_VERTEX
	|| (bridge_entity_type == requested_entity_type
	    && bridge_entity_type != iBase_ALL_TYPES)
	|| requested_entity_type == iBase_EDGE) {
      // This handles eleven of twenty-five cases by returning nothing!
      *adjacent_entities_size = 0;
      TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		     *adjacent_entities_size, iBase_EntityHandle);
    }
    else if ((requested_entity_type == iBase_REGION)
	|| (requested_entity_type == iBase_FACE
	    && bridge_entity_type == iBase_EDGE)) {
      // This handles another four cases (->e->r, ->f->r, ->a->r, ->e->f).
      getEntAdj(entity_handle, requested_entity_type, adjacent_entities,
		adjacent_entities_allocated, adjacent_entities_size, err);
    }
    else {
      // In all of these cases, we need to compute the neighborhood of
      // the vertex.
      Vert *pV = dynamic_cast<Vert*>(entity_handle);
      std::set<Cell*> spCInc;
      std::set<Vert*> spVNeigh;
      std::set<Face*> spFNearby;
      findNeighborhoodInfo(pV, spCInc, spVNeigh, NULL, NULL, &spFNearby);
      int iNCells = spCInc.size();
      int iNFaces = spFNearby.size();
      int iNVerts = spVNeigh.size();

      if (bridge_entity_type >= iBase_REGION) {
	// A total of six cases: bridge = region or all, with requested
	// = vertex, face, or all
	if (requested_entity_type == iBase_ALL_TYPES) {
	  // Two more cases (regi->all and ->all->all).  These are easy:
	  // they take everything

	  if (bridge_entity_type == iBase_REGION) {
	    *adjacent_entities_size = iNFaces + iNVerts;
	  }
	  else {
	    *adjacent_entities_size = iNCells + iNFaces + iNVerts;
	  }

	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  int ii = 0;
	  if (bridge_entity_type != iBase_REGION) {
	    std::copy(spCInc.begin(), spCInc.end(), (*adjacent_entities) + ii);
	    ii += iNCells;
	  }

	  std::copy(spVNeigh.begin(), spVNeigh.end(),
		    (*adjacent_entities) + ii);
	  ii += iNVerts;

	  std::copy(spFNearby.begin(), spFNearby.end(),
		    (*adjacent_entities) + ii);
	  ii += iNFaces;

	  assert(ii == *adjacent_entities_size);
	} // Done with requested_type = ALL_TYPES for region and all
	  // bridge
	else if (requested_entity_type == iBase_FACE) {
	  // Two cases (->all->face, ->regi->face); either way, we need
	  // all the faces.
	  *adjacent_entities_size = iNFaces;
	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  iBase_EntityHandle* ehResult = std::copy(spFNearby.begin(),
						   spFNearby.end(),
						   *adjacent_entities);
	  assert(ehResult - *adjacent_entities == *adjacent_entities_size);
	}
	else {
	  assert(requested_entity_type == iBase_VERTEX);
	  // Two cases (->all->vert, ->face->vert); either way, we need
	  // all the verts.
	  *adjacent_entities_size = iNVerts;
	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  iBase_EntityHandle* ehResult = std::copy(spVNeigh.begin(),
						   spVNeigh.end(),
						   *adjacent_entities);
	  assert(ehResult - *adjacent_entities == *adjacent_entities_size);
	}
      } // Done with bridging via region or all
      else {
	// The other combinations (->face->all, ->face->vert, edge->all,
	// edge->vert) are tricky, or at least tedious.  For now (read:
	// until it turns out to be an issue), these will get done as
	// two single adjacency calls.  For edges, this implies an empty
	// result, but that's fine.
	if (bridge_entity_type == iBase_EDGE) {
	  *adjacent_entities_size = 0;
	}
	else {
	  assert(
	      bridge_entity_type == iBase_FACE
		  && (requested_entity_type == iBase_ALL_TYPES
		      || requested_entity_type == iBase_VERTEX));
	  iBase_EntityHandle *intermediate_ents = NULL;
	  int intermediate_ents_allocated = 0, intermediate_ents_size = 0;
	  getEntAdj(entity_handle, bridge_entity_type, &intermediate_ents,
		    &intermediate_ents_allocated, &intermediate_ents_size, err);

	  iBase_EntityHandle *temp_ents = NULL;
	  int temp_ents_allocated = 0, temp_ents_size = 0;
	  int *offset = NULL, offset_allocated = 0, offset_size;
	  getEntArrAdj(intermediate_ents, intermediate_ents_size,
		       requested_entity_type, &temp_ents, &temp_ents_allocated,
		       &temp_ents_size, &offset, &offset_allocated,
		       &offset_size, err);
	  std::set<iBase_EntityHandle> setOutput(temp_ents,
						 temp_ents + temp_ents_size);
	  setOutput.erase(entity_handle);
	  *adjacent_entities_size = setOutput.size();
	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  std::copy(setOutput.begin(), setOutput.end(), *adjacent_entities);
	  if (intermediate_ents)
	    free(intermediate_ents);
	  if (temp_ents)
	    free(temp_ents);
	  if (offset)
	    free(offset);
	}
      } // Done with the face/edge->all/vert cases
    } // Done with vert->A->B cases that aren't just first adjacencies
  } // Done with vert->A->B
  else if ((bridge_entity_type == input_type
      || bridge_entity_type == requested_entity_type)
      && bridge_entity_type != iBase_ALL_TYPES) {
    *adjacent_entities_size = 0;
    TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		   *adjacent_entities_size, iBase_EntityHandle);
  }
  else {
    // No better alternative algorithmically than to make two
    // consecutive adjacency calls (one entity, one array).
    iBase_EntityHandle *intermediate_ents = NULL;
    int intermediate_ents_allocated = 0, intermediate_ents_size;
    getEntAdj(entity_handle, bridge_entity_type, &intermediate_ents,
	      &intermediate_ents_allocated, &intermediate_ents_size, err);

    iBase_EntityHandle *temp_ents = NULL;
    int temp_ents_allocated = 0, temp_ents_size = 0;
    int *offset = NULL, offset_allocated = 0, offset_size;
    getEntArrAdj(intermediate_ents, intermediate_ents_size,
		 requested_entity_type, &temp_ents, &temp_ents_allocated,
		 &temp_ents_size, &offset, &offset_allocated, &offset_size,
		 err);

    std::set<iBase_EntityHandle> setOutput(temp_ents,
					   temp_ents + temp_ents_size);
    setOutput.erase(entity_handle);
    *adjacent_entities_size = setOutput.size();
    TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		   *adjacent_entities_size, iBase_EntityHandle);
    std::copy(setOutput.begin(), setOutput.end(), *adjacent_entities);

    if (intermediate_ents)
      free(intermediate_ents);
    if (temp_ents)
      free(temp_ents);
    if (offset)
      free(offset);
  }
}

