#include "iMesh.h"

#include <vector>
#include "GR_Entity.h"
#include "GR_Geometry.h"
#include "GR_Mesh2D.h"
#include "iMesh_GRUMMP.hh"
#include "iMesh_GRUMMP_misc.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"

#undef SLOW_BUT_SURE

using namespace ITAPS_GRUMMP;

iMesh_Plane::iMesh_Plane(Mesh2D *pM2D) :
    iMesh_Base(ePlane)
{
  pM = pM2D;
  externalData = true;
}

void
iMesh_Plane::vCreateMesh()
{
  Mesh2D *pM2D = new Mesh2D();
  pM = pM2D;
}

void
iMesh_Plane::load(const char *name, const char* /*options*/, int* err)
{
  m_err = *err = iBase_SUCCESS;
  Mesh2D *pM2D = new Mesh2D();
  pM2D->readFromFile(name);
  pM = pM2D;
}

void
iMesh_Plane::save(const char *name, const char* /*options*/, int *err)
{
  m_err = *err = iBase_SUCCESS;
  Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pM);
  writeVTK(*pM2D, name);
}

void
iMesh_Plane::setAdjTable(int* adjacency_table, int* err)
{
  m_err = *err = iBase_SUCCESS;
  if (adjacency_table[3] != iBase_UNAVAILABLE
      || adjacency_table[7] != iBase_UNAVAILABLE
      || adjacency_table[11] != iBase_UNAVAILABLE
      || adjacency_table[12] != iBase_UNAVAILABLE
      || adjacency_table[13] != iBase_UNAVAILABLE
      || adjacency_table[14] != iBase_UNAVAILABLE
      || adjacency_table[15] != iBase_UNAVAILABLE) {
    m_err = *err = iBase_NOT_SUPPORTED;
  }
  getAdjTable(adjacency_table);
}

void
iMesh_Plane::getAdjTable(int* adjacency_table)
{
  // Array is known to be okay before this call.

  // The contents of this adjacency table should be:
  //    V E F R
  //  V I L L U
  //  E I I I U
  //  F I I I U
  //  R U U U U
  // where I = immediate
  //       L = local traversal
  //       U = unavailable

  // The adjacency info available in GRUMMP is always the same.
  adjacency_table[0] = iBase_AVAILABLE;
  adjacency_table[1] = iBase_ALL_ORDER_1;
  adjacency_table[2] = iBase_ALL_ORDER_1;
  adjacency_table[3] = iBase_UNAVAILABLE;

  adjacency_table[4] = iBase_ALL_ORDER_1;
  adjacency_table[5] = iBase_AVAILABLE;
  adjacency_table[6] = iBase_ALL_ORDER_1;
  adjacency_table[7] = iBase_UNAVAILABLE;

  adjacency_table[8] = iBase_ALL_ORDER_1;
  adjacency_table[9] = iBase_ALL_ORDER_1;
  adjacency_table[10] = iBase_AVAILABLE;
  adjacency_table[11] = iBase_UNAVAILABLE;

  adjacency_table[12] = iBase_UNAVAILABLE;
  adjacency_table[13] = iBase_UNAVAILABLE;
  adjacency_table[14] = iBase_UNAVAILABLE;
  adjacency_table[15] = iBase_UNAVAILABLE;
}

void
iMesh_Plane::getEntities(const int entity_type, const int entity_topology,
			 iBase_EntityHandle** entity_handles,
			 int* entity_handles_allocated,
			 int* entity_handles_size, int *err)
{
  m_err = *err = iBase_SUCCESS;
  Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pM);

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
    *entity_handles_size = getNumVertices();
  }
  else if (entity_topology == iMesh_LINE_SEGMENT
      || (entity_topology == iMesh_ALL_TOPOLOGIES && entity_type == iBase_EDGE)) {
    ////////////////////////////////////////////////
    // Returning all edges
    ////////////////////////////////////////////////
    int iNFace_db = pM2D->getNumFaces();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNFace_db,
		   iBase_EntityHandle);

    int iF = 0, iEH = 0;
    // Skip to the first non-deleted face
    while (iF < iNFace_db && pM2D->getFace(iF)->isDeleted())
      iF++;

    while (iF < iNFace_db) {
      (*entity_handles)[iEH] = pM2D->getFace(iF);
      iEH++;
      iF++;
      while (iF < iNFace_db && pM2D->getFace(iF)->isDeleted())
	iF++;
    }
    *entity_handles_size = iEH;
    assert(iEH <= (*entity_handles_allocated));
  } // Done with edges.
  else if (entity_topology == iMesh_TRIANGLE) {
    ////////////////////////////////////////////////
    // Returning all triangular cells.
    ////////////////////////////////////////////////
    int iNTri_db = pM2D->getNumTriCells();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNTri_db,
		   iBase_EntityHandle);

    int iTri = 0, iEH = 0;
    // Skip to the first non-deleted tri cell
    while (iTri < iNTri_db && pM2D->getCell(iTri)->isDeleted())
      iTri++;

    while (iTri < iNTri_db) {
      (*entity_handles)[iEH] = pM2D->getCell(iTri);
      iEH++;
      iTri++;
      while (iTri < iNTri_db && pM2D->getCell(iTri)->isDeleted())
	iTri++;
    }
    *entity_handles_size = iEH;
    assert(iEH <= (*entity_handles_allocated));
  } // Done with triangular cells.
  else if (entity_topology == iMesh_QUADRILATERAL) {
    ////////////////////////////////////////////////
    // Returning all quadrilateral cells.
    ////////////////////////////////////////////////
    int iNTri_db = pM2D->getNumTriCells();
    int iNQuad_db = pM2D->getNumQuadCells();
    int iNCell_db = iNTri_db + iNQuad_db;

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNQuad_db,
		   iBase_EntityHandle);

    int iQuad = iNTri_db, iEH = 0;
    // Skip to the first non-deleted quad cell
    while (iQuad < iNCell_db && pM2D->getCell(iQuad)->isDeleted())
      iQuad++;

    while (iQuad < iNCell_db) {
      (*entity_handles)[iEH] = pM2D->getCell(iQuad);
      iEH++;
      iQuad++;
      while (iQuad < iNCell_db && pM2D->getCell(iQuad)->isDeleted())
	iQuad++;
    }
    *entity_handles_size = iEH;
    assert(iEH <= (*entity_handles_allocated));
  } // Done with quadrilateral cells.
  else if (entity_topology == iMesh_ALL_TOPOLOGIES
      && entity_type == iBase_FACE) {
    ////////////////////////////////////////////////
    // Returning all cells.
    ////////////////////////////////////////////////
    int iNCell_db = pM2D->getNumCells();

    TRY_ARRAY_SIZE(entity_handles, entity_handles_allocated, iNCell_db,
		   iBase_EntityHandle);

    int iCell = 0, iEH = 0;
    // Skip to the first non-deleted cell
    while (iCell < iNCell_db && pM2D->getCell(iCell)->isDeleted())
      iCell++;

    while (iCell < iNCell_db) {
      (*entity_handles)[iEH] = pM2D->getCell(iCell);
      iEH++;
      iCell++;
      while (iCell < iNCell_db && pM2D->getCell(iCell)->isDeleted())
	iCell++;
    }
    *entity_handles_size = iEH;
    assert(iEH <= (*entity_handles_allocated));
  } // Done with all cells.
  else if (entity_topology == iMesh_ALL_TOPOLOGIES
      && entity_type == iBase_ALL_TYPES) {
    ////////////////////////////////////////////////
    // Returning all entities!
    ////////////////////////////////////////////////

    vCopyRootSet(entity_handles, entity_handles_allocated, entity_handles_size);
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
  else {
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
iMesh_Plane::createEnt(/*in*/const int new_entity_topology,
/*in*/const iBase_EntityHandle* lower_order_entity_handles,
		       /*in*/const int lower_order_entity_handles_size,
		       /*out*/iBase_EntityHandle* new_entity_handle,
		       /*out*/int* status, int *err)
{
  Mesh2D *pM2D = dynamic_cast<Mesh2D*>(pM);

  switch (new_entity_topology)
    {
    case iMesh_POINT:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      *new_entity_handle = NULL;
      *status = iBase_CREATION_FAILED;
      break;
    case iMesh_LINE_SEGMENT:
      {
	// We're going to create a new EdgeFace!  That means two vertices.
	if (lower_order_entity_handles_size != 2) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	  break;
	}

#ifdef NDEBUG
	Vert *pV0 = static_cast<Vert*>(lower_order_entity_handles[0]);
	Vert *pV1 = static_cast<Vert*>(lower_order_entity_handles[1]);
#else
	Vert *pV0 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	Vert *pV1 =
	    dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));

	if (!pV0->isValid() || !pV1->isValid()) {
	  m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	  break;
	}
#endif

	bool qExist;
	Face *pF = pM2D->createFace(qExist, pV0, pV1);

	if (pF->isValid()) {
	  // Return this face as a handle, with status ALREADY_EXISTED
	  m_err = *err = iBase_SUCCESS;
	  *new_entity_handle = pF;
	  *status = qExist ? iBase_ALREADY_EXISTED : iBase_NEW;
	}
	else {
	  m_err = *err = iBase_ENTITY_CREATION_ERROR;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	}
	break;
      }
    case iMesh_TRIANGLE:
      {
	// New triangles to be built from faces or vertices.
	if (lower_order_entity_handles_size != 3) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	  break;
	}

	Entity *pE =
	    static_cast<Entity*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	int ET = pE->getEntType();
	Face *apF[3];

	if (ET == iBase_VERTEX) {
	  // Create
#ifdef NDEBUG
	  Vert *pV0 = static_cast<Vert*>(lower_order_entity_handles[0]);
	  Vert *pV1 = static_cast<Vert*>(lower_order_entity_handles[1]);
	  Vert *pV2 = static_cast<Vert*>(lower_order_entity_handles[2]);
#else
	  Vert *pV0 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	  Vert *pV1 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  Vert *pV2 =
	      dynamic_cast<Vert*>(static_cast<Entity*>(lower_order_entity_handles[2]));

	  if (!pV0->isValid() || !pV1->isValid() || !pV2->isValid()) {
	    m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	    *status = iBase_CREATION_FAILED;
	    *new_entity_handle = NULL;
	    break;
	  }
#endif
	  bool qExist;
	  apF[0] = pM2D->createFace(qExist, pV0, pV1);
	  apF[1] = pM2D->createFace(qExist, pV1, pV2);
	  apF[2] = pM2D->createFace(qExist, pV2, pV0);

	} // Done creating/identifying edges from verts.
	else {
#ifdef NDEBUG
	  apF[0] = static_cast<Face*>(lower_order_entity_handles[0]);
	  apF[1] = static_cast<Face*>(lower_order_entity_handles[1]);
	  apF[2] = static_cast<Face*>(lower_order_entity_handles[2]);
#else
	  apF[0] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	  apF[1] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  apF[2] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
#endif

	  // The following also handles the case of getting handed weird
	  // stuff instead of legitimate edge faces.
	  if (!apF[0]->isValid() || !apF[1]->isValid() || !apF[2]->isValid()) {
	    m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	    *status = iBase_CREATION_FAILED;
	    *new_entity_handle = NULL;
	    break;
	  }
	}

	// Cell might exist already.
	Cell *pC01 = findCommonCell(apF[0], apF[1]);
	if (pC01->isValid()) {
	  Cell *pC12 = findCommonCell(apF[1], apF[2]);
	  if (pC12->isValid() && pC01 == pC12) {
	    // Hit a cell!
	    m_err = *err = iBase_SUCCESS;
	    *new_entity_handle = pC01;
	    *status = iBase_ALREADY_EXISTED;
	    return; // Nothing more to do.
	  }
	}

	TriCell* pTC = pM2D->createTriCell(apF[0], apF[1], apF[2]);
	if (pTC->isValid()) {
	  m_err = *err = iBase_SUCCESS;
	  *new_entity_handle = pTC;
	  *status = iBase_NEW;
	  pTC->canonicalizeFaceOrder();
	  assert(pTC->doFullCheck());
	} // Created the cell!
	else {
	  m_err = *err = iBase_ENTITY_CREATION_ERROR;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	} // Creation error.
	break;
      } // Done with tris
    case iMesh_QUADRILATERAL:
      {
	pM2D->allowNonSimplicial();
	// New quadrilaterals to be built from faces or vertices.
	if (lower_order_entity_handles_size != 4) {
	  m_err = *err = iBase_INVALID_ENTITY_COUNT;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	  break;
	}

	Entity *pE =
	    static_cast<Entity*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	int ET = pE->getEntType();
	Face *apF[4];

	if (ET == iBase_VERTEX) {
#ifdef NDEBUG
	  Vert *pV0 = static_cast<Vert*>(lower_order_entity_handles[0]);
	  Vert *pV1 = static_cast<Vert*>(lower_order_entity_handles[1]);
	  Vert *pV2 = static_cast<Vert*>(lower_order_entity_handles[2]);
	  Vert *pV3 = static_cast<Vert*>(lower_order_entity_handles[3]);
#else
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
	    *status = iBase_CREATION_FAILED;
	    *new_entity_handle = NULL;
	    break;
	  }
#endif

	  bool qExist;
	  apF[0] = pM2D->createFace(qExist, pV0, pV1);
	  apF[1] = pM2D->createFace(qExist, pV1, pV2);
	  apF[2] = pM2D->createFace(qExist, pV2, pV3);
	  apF[3] = pM2D->createFace(qExist, pV3, pV0);
	} // Done identifying and/or creating edges.
	else {
#ifdef NDEBUG
	  apF[0] = static_cast<Face*>(lower_order_entity_handles[0]);
	  apF[1] = static_cast<Face*>(lower_order_entity_handles[1]);
	  apF[2] = static_cast<Face*>(lower_order_entity_handles[2]);
	  apF[3] = static_cast<Face*>(lower_order_entity_handles[3]);
#else
	  apF[0] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[0]));
	  apF[1] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[1]));
	  apF[2] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[2]));
	  apF[3] =
	      dynamic_cast<Face*>(static_cast<Entity*>(lower_order_entity_handles[3]));
#endif

	  // The following also handles the case of getting handed weird
	  // stuff instead of legitimate edge faces.
	  if (!apF[0]->isValid() || !apF[1]->isValid() || !apF[2]->isValid()
	      || !apF[3]->isValid()) {
	    m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	    *status = iBase_CREATION_FAILED;
	    *new_entity_handle = NULL;
	    break;
	  }
	}

	Cell *pC01 = findCommonCell(apF[0], apF[1]);
	if (pC01->isValid()) {
	  Cell *pC23 = findCommonCell(apF[2], apF[3]);
	  if (pC23->isValid() && pC01 == pC23) {
	    m_err = *err = iBase_SUCCESS;
	    *new_entity_handle = pC01;
	    *status = iBase_ALREADY_EXISTED;
	    break;
	  }
	}

	QuadCell* pQC = pM2D->createQuadCell(apF[0], apF[1], apF[2], apF[3]);

	if (pQC->isValid()) {
	  m_err = *err = iBase_SUCCESS;
	  *new_entity_handle = pQC;
	  *status = iBase_NEW;
	  pQC->canonicalizeFaceOrder();
	  assert(pQC->doFullCheck());
	} // Got the cell
	else {
	  m_err = *err = iBase_ENTITY_CREATION_ERROR;
	  *status = iBase_CREATION_FAILED;
	  *new_entity_handle = NULL;
	} // Failed, one way or another.
	break;
      } // Done creating quadrilaterals

    case iMesh_TETRAHEDRON:
    case iMesh_PYRAMID:
    case iMesh_PRISM:
    case iMesh_HEXAHEDRON:
    case iMesh_SEPTAHEDRON:
    case iMesh_POLYHEDRON:
    case iMesh_POLYGON:
      // Not in 2D, thanks.
      m_err = *err = iBase_NOT_SUPPORTED;
      break;
    case iMesh_ALL_TOPOLOGIES:
      // Nonsense request.
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      break;
    }
}

void
iMesh_Plane::createEntArr(/*in*/const int new_entity_topology,
/*in*/const iBase_EntityHandle* lower_order_entity_handles,
			  /*in*/const int lower_order_entity_handles_size,
			  /*out*/iBase_EntityHandle** new_entity_handles,
			  /*out*/int* new_entity_handles_allocated,
			  /*out*/int* new_entity_handles_size,
			  /*inout*/int** status,
			  /*inout*/int* status_allocated,
			  /*out*/int* status_size, int *err)
{
  m_err = *err = iBase_SUCCESS;

  *new_entity_handles_size = 0;
  *status_size = 0;
  int vertsPerEnt;
  switch (new_entity_topology)
    {
    case iMesh_POINT:
      m_err = *err = iBase_INVALID_ENTITY_TOPOLOGY;
      return;
      break;
    case iMesh_LINE_SEGMENT:
      vertsPerEnt = 2;
      break;
    case iMesh_TRIANGLE:
      vertsPerEnt = 3;
      break;
    case iMesh_QUADRILATERAL:
      vertsPerEnt = 4;
      break;
    case iMesh_TETRAHEDRON:
    case iMesh_PRISM:
    case iMesh_PYRAMID:
    case iMesh_HEXAHEDRON:
    case iMesh_SEPTAHEDRON:
    case iMesh_POLYHEDRON:
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
iMesh_Plane::deleteEnt(iBase_EntityHandle entity_handle, int *err)
{
  m_err = *err = iBase_SUCCESS;

  Entity *pE = static_cast<Entity*>(entity_handle);
  switch (pE->getEntTopology())
    {
    case iMesh_QUADRILATERAL:
    case iMesh_TRIANGLE:
      {
	// In this case, there are no upward adjacencies, so it's always
	// safe to get rid of these.
	vRemoveAllReferencesTo(pE);
	// Remove connections from its edges.
#ifdef NDEBUG
	Cell *pC = static_cast<Cell*>(pE);
#else
	Cell *pC = dynamic_cast<Cell*>(pE);
#endif
	pM->deleteCell(pC);
      }
      break;
    case iMesh_LINE_SEGMENT:
      // These could have cells that depend on them.
      {
#ifdef NDEBUG
	EdgeFace *pEF = static_cast<EdgeFace*>(pE);
#else
	EdgeFace *pEF = dynamic_cast<EdgeFace*>(pE);
	assert(pEF);
#endif

	bool qLeftOkay = !pEF->getLeftCell()->isValid()
	    || pEF->getLeftCell()->isDeleted();
	bool qRightOkay = !pEF->getRightCell()->isValid()
	    || pEF->getRightCell()->isDeleted();
	if (pEF->isBdryFace()) {
	  BFace *pBF = pBFInvalidBFace;
	  // Identify the bdry face, if it's still attached.
	  if (qLeftOkay) {
	    // Only chance is on the right; if it's a BFace, the dynamic
	    // cast will put it in pBF, otherwise NULL.
#ifdef NDEBUG
	    pBF = static_cast<BFace*>(pEF->getRightCell());
#else
	    pBF = dynamic_cast<BFace*>(pEF->getRightCell());
#endif
	  }
	  else if (qRightOkay) {
	    // Same thing, other side.
#ifdef NDEBUG
	    pBF = static_cast<BFace*>(pEF->getLeftCell());
#else
	    pBF = dynamic_cast<BFace*>(pEF->getLeftCell());
#endif
	  }
	  pM->deleteBFace(pBF);
	  qLeftOkay = qRightOkay = true;
	}

	if (qLeftOkay && qRightOkay) {
	  vRemoveAllReferencesTo(pE);
	  pM->deleteFace(pEF);
	}
	else {
	  m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	}
	break;
      }
    case iMesh_POINT:
      // Every vertex that's connected to anything must have a valid
      // hint face.
      {
#ifdef NDEBUG
	Vert *pV = static_cast<Vert*>(pE);
#else
	Vert *pV = dynamic_cast<Vert*>(pE);
	assert(pV);
#endif
	if (!pV->isDeleted()) {
	  Face *pF = pV->getHintFace();
	  if (!pF->isValid() || !pF->hasVert(pV) || pF->isDeleted()) {
	    vRemoveAllReferencesTo(pE);
	    pM->deleteVert(pV);
	  }
	  else {
	    m_err = *err = iBase_INVALID_ENTITY_HANDLE;
	  }
	}
	break;
      }
    default:
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
iMesh_Plane::getEntAdj(const iBase_EntityHandle entity_handle,
		       const int entity_type_requested,
		       iBase_EntityHandle** adj_entity_handles,
		       int* adj_entity_handles_allocated,
		       int* adj_entity_handles_size, int *err)
{
  Entity *pE = static_cast<Entity*>(entity_handle);
  int eET = pE->getEntType();
  if (eET == entity_type_requested) {
    // Not supposed to return any data in this case.
    *adj_entity_handles_size = 0;
    TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, 0,
		   iBase_EntityHandle);
  }
  else if (eET == iBase_VERTEX && entity_type_requested != iBase_VERTEX) {
    // Now get fancy with a call that computes the neighborhood.  This
    // is the only expensive adjacency request, requiring a local mesh
    // traversal.
    Vert *pV = dynamic_cast<Vert*>(pE);
    assert(pV->isValid());
    std::set<Cell*> spCInc;
    std::set<Vert*> spVNeigh;
    std::set<BFace*> spBFInc;
    bool qBdryVert;
    std::set<Face*> spFNearby;
    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFInc, &qBdryVert,
			 &spFNearby);
    switch (entity_type_requested)
      {
      case iBase_EDGE:
	// Add faces incident on pV.
	*adj_entity_handles_size = std::count_if(spFNearby.begin(),
						 spFNearby.end(),
						 bind2nd(qHasVert(), pV));
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		       *adj_entity_handles_size, iBase_EntityHandle);
	std::remove_copy_if(spFNearby.begin(), spFNearby.end(),
			    *adj_entity_handles, bind2nd(qNotHasVert(), pV));
	break;
      case iBase_FACE:
	*adj_entity_handles_size = spCInc.size();
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		       *adj_entity_handles_size, iBase_EntityHandle);
	std::copy(spCInc.begin(), spCInc.end(), *adj_entity_handles);
	break;
      case iBase_REGION:
	*adj_entity_handles_size = 0;
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		       *adj_entity_handles_size, iBase_EntityHandle);
	break;
      case iBase_ALL_TYPES:
	*adj_entity_handles_size = spCInc.size()
	    + std::count_if(spFNearby.begin(), spFNearby.end(),
			    bind2nd(qHasVert(), pV));
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
		       *adj_entity_handles_size, iBase_EntityHandle);
	{
	  iBase_EntityHandle* tmpPtr = std::remove_copy_if(
	      spFNearby.begin(), spFNearby.end(), *adj_entity_handles,
	      bind2nd(qNotHasVert(), pV));
	  std::copy(spCInc.begin(), spCInc.end(), tmpPtr);
	}
	break;
      } // End switch
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
      case iBase_EDGE:
	{
	  *adj_entity_handles_size = pE->getNumFaces();
	  TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated,
			 *adj_entity_handles_size, iBase_EntityHandle);
	  pE->getAllFaceHandles(*adj_entity_handles);
	  break;
	}
      case iBase_FACE:
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
      case iBase_REGION:
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
	}
	break;
      } // End of switch over requested type
  } // The general case (not vertex nor requesting type of current entity
}

void
iMesh_Plane::getEntArrAdj(const iBase_EntityHandle* entity_handles,
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

  *offset_size = entity_handles_size + 1;
  TRY_ARRAY_SIZE(offset, offset_allocated, *offset_size, int);

  switch (entity_type_requested)
    {
    case iBase_VERTEX:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE = 0, iEOut = 0;
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
    case iBase_REGION:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE = 0, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  // Since there are no regions in 2D, do nothing here.
	}
	(*offset)[iE] = iEOut;

	// Now set up the output array and transcribe the data.
	TRY_ARRAY_SIZE(adj_entity_handles, adj_entity_handles_allocated, iEOut,
		       iBase_EntityHandle);
	*adj_entity_handles_size = iEOut;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  // Do nothing, because there are no regions in the 2D mesh.
	  // In fact, applications with any sense have multiple ways to
	  // decide not to get here, but don't throw an error anyway,
	  // because the spec implies we shouldn't.
	}
	break;
      }
    case iBase_EDGE:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE = 0, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX
	      && entity_type_requested != iBase_VERTEX) {
	    // Now get fancy with a call that computes the neighborhood.  This
	    // is the only expensive adjacency request, requiring a local mesh
	    // traversal.
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    std::set<Cell*> spCInc;
	    std::set<Vert*> spVNeigh;
	    std::set<BFace*> spBFInc;
	    bool qBdryVert;
	    std::set<Face*> spFNearby;
	    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFInc, &qBdryVert,
				 &spFNearby);
	    int iCount = std::count_if(spFNearby.begin(), spFNearby.end(),
				       bind2nd(qHasVert(), pV));
	    iEOut += iCount;
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

	iBase_EntityHandle* pEH = *adj_entity_handles;

	for (iE = 0; iE < entity_handles_size; iE++) {
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX
	      && entity_type_requested != iBase_VERTEX) {
	    // Now get fancy with a call that computes the neighborhood.  This
	    // is the only expensive adjacency request, requiring a local mesh
	    // traversal.
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    std::set<Cell*> spCInc;
	    std::set<Vert*> spVNeigh;
	    std::set<BFace*> spBFInc;
	    bool qBdryVert;
	    std::set<Face*> spFNearby;
	    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFInc, &qBdryVert,
				 &spFNearby);

	    // Add faces incident on pV.
	    pEH = std::remove_copy_if(spFNearby.begin(), spFNearby.end(), pEH,
				      bind2nd(qNotHasVert(), pV));
	    assert(pEH <= (*adj_entity_handles + iEOut));
	  }
	  else {
	    pE->getAllFaceHandles(*adj_entity_handles + (*offset)[iE]);
	  }
	}
	break;
      }

    case iBase_FACE:
      {
	// Pre-compute the output size, and put the correct values in the
	// offset array.
	int iE = 0, iEOut = 0;
	for (iE = 0; iE < entity_handles_size; iE++) {
	  (*offset)[iE] = iEOut;
	  Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	  int ET = pE->getEntType();
	  if (ET == entity_type_requested) {
	    // Not supposed to return any data in this case.
	  }
	  else if (ET == iBase_VERTEX) {
	    // Now get fancy with a call that computes the neighborhood.  This
	    // is the only expensive adjacency request, requiring a local mesh
	    // traversal.
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    std::set<Cell*> spCInc;
	    std::set<Vert*> spVNeigh;
	    std::set<BFace*> spBFInc;
	    bool qBdryVert;
	    std::set<Face*> spFNearby;
	    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFInc, &qBdryVert,
				 &spFNearby);
	    iEOut += spCInc.size();
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
	    // Now get fancy with a call that computes the neighborhood.  This
	    // is the only expensive adjacency request, requiring a local mesh
	    // traversal.
	    Vert *pV = dynamic_cast<Vert*>(pE);
	    assert(pV->isValid());
	    std::set<Cell*> spCInc;
	    std::set<Vert*> spVNeigh;
	    std::set<BFace*> spBFInc;
	    bool qBdryVert;
	    std::set<Face*> spFNearby;
	    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFInc, &qBdryVert,
				 &spFNearby);
	    pEH = std::copy(spCInc.begin(), spCInc.end(), pEH);
	    assert(pEH == (*adj_entity_handles) + (*offset)[iE + 1]);
	  }
	  else {
	    pE->getAllCellHandles(*adj_entity_handles + (*offset)[iE]);
	  }
	}
	break;
      }
    case iBase_ALL_TYPES:
      // Pre-compute the output size, and put the correct values in the
      // offset array.
      int iE = 0, iEOut = 0;
      for (iE = 0; iE < entity_handles_size; iE++) {
	(*offset)[iE] = iEOut;
	Entity *pE = static_cast<Entity*>(entity_handles[iE]);
	int ET = pE->getEntType();
	if (ET == iBase_VERTEX) {
	  // Now get fancy with a call that computes the neighborhood.  This
	  // is the only expensive adjacency request, requiring a local mesh
	  // traversal.
	  Vert *pV = dynamic_cast<Vert*>(pE);
	  assert(pV->isValid());
	  std::set<Cell*> spCInc;
	  std::set<Vert*> spVNeigh;
	  std::set<Face*> spFNearby;
	  findNeighborhoodInfo(pV, spCInc, spVNeigh, NULL, NULL, &spFNearby);
	  int iCount = std::count_if(spFNearby.begin(), spFNearby.end(),
				     bind2nd(qHasVert(), pV));
	  iEOut += iCount + spCInc.size();
	}
	else {
	  iEOut += pE->getNumFaces() + pE->getNumCells() + pE->getNumVerts();
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
	if (ET == iBase_VERTEX) {
	  // Now get fancy with a call that computes the neighborhood.  This
	  // is the only expensive adjacency request, requiring a local mesh
	  // traversal.
	  Vert *pV = dynamic_cast<Vert*>(pE);
	  assert(pV->isValid());
	  std::set<Cell*> spCInc;
	  std::set<Vert*> spVNeigh;
	  std::set<Face*> spFNearby;
	  findNeighborhoodInfo(pV, spCInc, spVNeigh, NULL, NULL, &spFNearby);

	  // Add faces incident on pV.
	  pEH = std::remove_copy_if(spFNearby.begin(), spFNearby.end(), pEH,
				    bind2nd(qNotHasVert(), pV));
	  pEH = std::copy(spCInc.begin(), spCInc.end(), pEH);
	  assert(pEH <= (*adj_entity_handles + iEOut));
	}
	else {
	  pE->getAllFaceHandles(pEH);
	  pEH += pE->getNumFaces();

	  pE->getAllCellHandles(pEH);
	  pEH += pE->getNumCells();

	  pE->getAllVertHandles(pEH);
	  pEH += pE->getNumVerts();
	}
      }
      break;
    } // End of the switch
} // End of iMesh_getEntArrAdj

void
iMesh_Plane::getEnt2ndAdj(const iBase_EntityHandle entity_handle,
			  const int bridge_entity_type,
			  const int requested_entity_type,
			  iBase_EntityHandle** adjacent_entities,
			  int* adjacent_entities_allocated,
			  int* adjacent_entities_size, int* err)
{
  m_err = *err = iBase_SUCCESS;
  Entity *pE = static_cast<Entity*>(entity_handle);
  int input_type = pE->getEntType();

  // If regions are involved at all, return no data.
  if (bridge_entity_type == iBase_REGION
      || requested_entity_type == iBase_REGION) {
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
    // adjacency calls, because the second one is never needed.

    // Possible cases:
    // FACE and ALL_TYPES as bridge type turn out to be identical.

    //  ->all->all   = everything from neighborhood calc
    //  ->all->face  = 1st adj faces
    //  ->all->edge  = all edges from neighborhood calc
    //  ->all->vert  = all verts from neighborhood calc

    //  ->face->all  = all of neighborhood calc except faces
    //  ->face->face = nothing
    //  ->face->edge = all edges from neighborhood calc
    //  ->face->vert = all verts from neighborhood calc

    //  ->edge->all  = all of neighborhood calc except edges
    //  ->edge->face = 1st adj faces
    //  ->edge->edge = nothing
    //  ->edge->vert = most verts from neighborhood calc

    //  ->vert->all  = nothing
    //  ->vert->face = nothing
    //  ->vert->edge = nothing
    //  ->vert->vert = nothing

    if (bridge_entity_type == iBase_VERTEX
	|| (bridge_entity_type == requested_entity_type
	    && bridge_entity_type != iBase_ALL_TYPES)) {
      // These six cases do nothing at all.
      *adjacent_entities_size = 0;
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
      int iNEdges = spFNearby.size();
      int iNVerts = spVNeigh.size();
      // The following is true both for interior and bdry verts.
      assert(iNCells + iNVerts == iNEdges);

      if (requested_entity_type == iBase_FACE) {
	// Two cases (->e->f, ->a->f): return all neighborhood cells
	*adjacent_entities_size = iNCells;
	TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		       *adjacent_entities_size, iBase_EntityHandle);
	std::copy(spCInc.begin(), spCInc.end(), (*adjacent_entities));
      }
      else if (requested_entity_type == iBase_EDGE) {
	// Two cases (->f->e, ->a->e): return all neighborhood edges
	*adjacent_entities_size = iNEdges;
	TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		       *adjacent_entities_size, iBase_EntityHandle);
	std::copy(spFNearby.begin(), spFNearby.end(), (*adjacent_entities));
      }
      else if (requested_entity_type == iBase_VERTEX) {
	if (bridge_entity_type == iBase_EDGE) {
	  // One case (->e->v): return some of the neighborhood verts
	  *adjacent_entities_size = std::count_if(spFNearby.begin(),
						  spFNearby.end(),
						  bind2nd(qHasVert(), pV));
	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  std::set<Face*>::iterator iter = spFNearby.begin(), iterEnd =
	      spFNearby.end();
	  int entIndex = 0;
	  for (; iter != iterEnd; ++iter) {
	    Face *pF = *iter;
	    if (pF->hasVert(pV)) {
	      if (pF->getVert(0) == pV) {
		(*adjacent_entities)[entIndex++] = pF->getVert(1);
	      }
	      else {
		assert(pF->getVert(1) == pV);
		(*adjacent_entities)[entIndex++] = pF->getVert(0);
	      }
	    }
	  }
	  assert(entIndex == *adjacent_entities_size);
	}
	else {
	  // Two cases (->f->v, ->a->v): return all neighborhood verts
	  *adjacent_entities_size = iNVerts;
	  TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
			 *adjacent_entities_size, iBase_EntityHandle);
	  std::copy(spVNeigh.begin(), spVNeigh.end(), (*adjacent_entities));
	}
      } // Done with VERTEX output
      else if (requested_entity_type == iBase_ALL_TYPES) {
	// Three more cases (->edge->all, ->face->all, and ->all->all)
	bool qReturnEdges = (bridge_entity_type != iBase_EDGE);
	bool qReturnCells = (bridge_entity_type != iBase_FACE);
	*adjacent_entities_size = iNVerts + (qReturnEdges ? iNEdges : 0)
	    + (qReturnCells ? iNCells : 0);

	TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		       *adjacent_entities_size, iBase_EntityHandle);
	int ii = 0;
	if (qReturnCells) {
	  std::copy(spCInc.begin(), spCInc.end(), (*adjacent_entities) + ii);
	  ii += iNCells;
	}

	if (qReturnEdges) {
	  std::copy(spFNearby.begin(), spFNearby.end(),
		    (*adjacent_entities) + ii);
	  ii += iNEdges;
	}

	// Include the verts regardless
	std::copy(spVNeigh.begin(), spVNeigh.end(), (*adjacent_entities) + ii);
	ii += iNVerts;

	assert(ii == *adjacent_entities_size);
      } // Done with requested_type = ALL_TYPES
    } // Done with all cases that require neighborhood calculation
  } // Done with all cases that have input_type == VERTEX
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

    iBase_EntityHandle *tmp_ents = NULL;
    int tmp_ents_allocated = 0, tmp_ents_size = 0;
    int *offset = NULL, offset_allocated = 0, offset_size;
    getEntArrAdj(intermediate_ents, intermediate_ents_size,
		 requested_entity_type, &tmp_ents, &tmp_ents_allocated,
		 &tmp_ents_size, &offset, &offset_allocated, &offset_size, err);
    // Need to uniquify this data.
    std::set<iBase_EntityHandle> tmpSet;
    tmpSet.insert(tmp_ents, tmp_ents + tmp_ents_size);
    tmpSet.erase(entity_handle);
    *adjacent_entities_size = tmpSet.size();
    TRY_ARRAY_SIZE(adjacent_entities, adjacent_entities_allocated,
		   *adjacent_entities_size, iBase_EntityHandle);
    std::copy(tmpSet.begin(), tmpSet.end(), (*adjacent_entities));

    if (intermediate_ents)
      free(intermediate_ents);
    if (offset)
      free(offset);
    if (tmp_ents)
      free(tmp_ents);
  }
}

