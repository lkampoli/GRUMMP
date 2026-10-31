#include "iMesh_unitTest.hh"

//@@ Test the core Mesh interface.

// Test stuff specific to interface iMesh_  Those functions are:
//    load, save, getRootSet, getGeometricDim, getDfltStorage,
//    getAdjTable, getNumOfType, getNumOfTopo, getAllVtxCoords,
//    getVtxCoordIndex, getEntities, getVtxArrCoords, and
//    getAdjEntities
//
//    Mark C. Miller, Mon Dec  6 10:39:08 PST 2010
//    Added logic to fail asap if input mesh doesn't appear to
//    have appropriate structure that the test assumes.

static void
setMinAdjTable3D(int aiAdjTable[16])
{
  // Minimal adjacency table for 3D: all major implementations can
  // support the entities requested here.
  aiAdjTable[0] = iBase_AVAILABLE;
  aiAdjTable[1] = iBase_UNAVAILABLE; // Don't care about any actual
				     // adjacencies for this purpose.
  aiAdjTable[2] = iBase_UNAVAILABLE;
  aiAdjTable[3] = iBase_UNAVAILABLE;
  aiAdjTable[4] = iBase_UNAVAILABLE;
  aiAdjTable[5] = iBase_UNAVAILABLE;
  aiAdjTable[6] = iBase_UNAVAILABLE;
  aiAdjTable[7] = iBase_UNAVAILABLE;
  aiAdjTable[8] = iBase_UNAVAILABLE;
  aiAdjTable[9] = iBase_UNAVAILABLE;
  aiAdjTable[10] = iBase_AVAILABLE;
  aiAdjTable[11] = iBase_UNAVAILABLE;
  aiAdjTable[12] = iBase_UNAVAILABLE;
  aiAdjTable[13] = iBase_UNAVAILABLE;
  aiAdjTable[14] = iBase_UNAVAILABLE;
  aiAdjTable[15] = iBase_AVAILABLE;
}

void
vTestSetAdjTable(string strFilename, int dim)
{
  int err;
  int aiAdjTable[16];

  cout << " Testing setAdjTable" << endl;
  // Test 1: Call setGeomDim first, then setAdjTable
  {
    cout << "  setGeomDim, then setAdjTable... " << endl;
    iMesh_Instance instance;
    iMesh_newMesh("", &instance, &err, 0);
    setMinAdjTable3D(aiAdjTable);
    iMesh_setGeometricDimension(instance, 3, &err);
    CHECK_ERR(err);
    iMesh_setAdjTable(instance, aiAdjTable, 16, &err);
    CHECK_ERR(err);
    iMesh_dtor(instance, &err);
    CHECK_ERR(err);
  }

  // Test 2: Call setAdjTable, then setGeomDim.
  {
    cout << "  setAdjTable, then setGeomDim... " << endl;
    iMesh_Instance instance;
    iMesh_newMesh("", &instance, &err, 0);
    setMinAdjTable3D(aiAdjTable);
    iMesh_setAdjTable(instance, aiAdjTable, 16, &err);
    CHECK_ERR(err);
    iMesh_setGeometricDimension(instance, 3, &err);
    CHECK_ERR(err);

    // Check whether the output is the same with no change to dimension.
    int *newTable, newTable_alloc = 0, newTable_size;
    iMesh_getAdjTable(instance, &newTable, &newTable_alloc, &newTable_size,
		      &err);

    bool qAdjTableUnchanged = true;
    for (int i = 0; i < 16; i++) {
      if (newTable[i] != aiAdjTable[i])
	qAdjTableUnchanged = false;
    }
    TEST(qAdjTableUnchanged);
    iMesh_dtor(instance, &err);
    CHECK_ERR(err);
    UT_FREE(newTable);
  }

  // TODO Figure out a way to have this test be compatible with surf meshes.
  // Test 3: Load a mesh, then call setAdjTable.
//  if (dim == 3) {
//    cout << "  load, then setAdjTable... " << endl;
//    iMesh_Instance instance;
//    iMesh_newMesh("", &instance, &err, 0);
//    iBase_EntitySetHandle tmpRootSet;
//    iMesh_getRootSet(instance, &tmpRootSet, &err);
//    iMesh_load(instance, tmpRootSet, strFilename.c_str(), "", &err,
//	       strFilename.length(), 0);
//    qIssueLoadMessage(strFilename.c_str(), err);
//
//    CHECK_ERR_CRITICAL(err);
//    setMinAdjTable3D(aiAdjTable);
//    iMesh_setAdjTable(instance, aiAdjTable, 16, &err);
//    CHECK_ERR(err);
//    bool qSetEntitiesProperly = true;
//    if (aiAdjTable[0] != iBase_AVAILABLE ||
//	aiAdjTable[10] != iBase_AVAILABLE ||
//	aiAdjTable[15] != iBase_AVAILABLE)
//      qSetEntitiesProperly = false;
//    TEST(qSetEntitiesProperly);
//    iMesh_dtor(instance, &err);
//    CHECK_ERR(err);
//  }
}

void
vTestMeshInterface(iMesh_Instance& instance, string strFilename)
{
  int err;
  // Modular: getRootSet must succeed
  iMesh_getRootSet(instance, &pvRootSet, &err);
  CHECK_ERR_CRITICAL(err);

  // Modular: load must succeed, or we're dead in the water.
  iMesh_load(instance, pvRootSet, strFilename.c_str(), "", &err,
	     strFilename.length(), 0);
  qIssueLoadMessage(strFilename.c_str(), err);
  CHECK_ERR_CRITICAL(err);

  // Modular: save is tested only to the extent that it shouldn't throw
  // an error.
//   iMesh_save(instance, pvRootSet, SAVE_NAME, "", &err, strlen(SAVE_NAME), 0);
//   CHECK_ERR(err);

  // Verify that the implementation returns a dimension.

  // Modular: getGeometricDim can be tested in isolation; with known
  // data, the answer can be tested specifically.
  cout << " Getting geometric dimension... ";
  int iGeomDim;
  iMesh_getGeometricDimension(instance, &iGeomDim, &err);
  CHECK_ERR(err);
  cout << iGeomDim << endl;
  TEST(iGeomDim > 0);
  if (iGeomDim < 2) {
    qIssueLoadMessage(strFilename.c_str(), iBase_FAILURE);
    CHECK_ERR_CRITICAL(iBase_FAILURE);
  }

  // Verify that the implementation returns a storage order preference.
  // Modular: getDfltStorage can be tested in isolation.
  cout << "Getting default storage order..." << endl;
  int SODef;
  iMesh_getDfltStorage(instance, &SODef, &err);
  CHECK_ERR(err);
  TEST(SODef == iBase_BLOCKED || SODef == iBase_INTERLEAVED);

  // Verify that the implementation returns a valid adjacency table.
  // Modular: getAdjTable can be tested in isolation.
  cout << " Getting adjacency table." << endl;
  int* aAdjTable;
  int iSize;
  int adj_table_allocated = 0;
  iMesh_getAdjTable(instance, &aAdjTable, &adj_table_allocated, &iSize, &err);
  CHECK_ERR(err);
  TESTEQ(iSize, 16);
  // Keep a copy for later reference.
  a2iAdjTable[0][0] = aAdjTable[0];
  a2iAdjTable[0][1] = aAdjTable[1];
  a2iAdjTable[0][2] = aAdjTable[2];
  a2iAdjTable[0][3] = aAdjTable[3];

  a2iAdjTable[1][0] = aAdjTable[4];
  a2iAdjTable[1][1] = aAdjTable[5];
  a2iAdjTable[1][2] = aAdjTable[6];
  a2iAdjTable[1][3] = aAdjTable[7];

  a2iAdjTable[2][0] = aAdjTable[8];
  a2iAdjTable[2][1] = aAdjTable[9];
  a2iAdjTable[2][2] = aAdjTable[10];
  a2iAdjTable[2][3] = aAdjTable[11];

  a2iAdjTable[3][0] = aAdjTable[12];
  a2iAdjTable[3][1] = aAdjTable[13];
  a2iAdjTable[3][2] = aAdjTable[14];
  a2iAdjTable[3][3] = aAdjTable[15];

  // Simple test for setAdjTable. Call it with whatever getAdjTable gave us.
  cout << " Setting adjacency table." << endl;
  iMesh_setAdjTable(instance, aAdjTable, iSize, &err);
  CHECK_ERR(err);

  // Check that impl. returned same values as in getAdjTable call
  // and that diagonal terms for verts and 2D/3D mesh ents are
  // properly set.
  bool qAdjTableOk = true;
  for (int i = 0; i < iSize && qAdjTableOk; i++) {
    // Should see only available/unavailable on diagonal
    // Should only see available for verts and highest dim entities.
    if (i == 0 || i == 5 || i == 10 || i == 15) {
      if ((i == 0 || (i == iGeomDim * 5)) && aAdjTable[i] != iBase_AVAILABLE)
	qAdjTableOk = false;
      else if (aAdjTable[i] != iBase_AVAILABLE
	  && aAdjTable[i] != iBase_UNAVAILABLE)
	qAdjTableOk = false;
    }
    if (aAdjTable[i] != a2iAdjTable[i / 4][i % 4])
      qAdjTableOk = false;
  }
  TEST(qAdjTableOk);

  UT_FREE(aAdjTable);

  vTestSetAdjTable(strFilename, iGeomDim);

  //@@ getEntities; entity type, topo
  bool qOK;
  int eEType;
  int eETopo;

  //@@@ Setup of legal argument combinations and expected return sizes.
  for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
    for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
      a2iEntCount[eEType][eETopo] = 0;
    }
  }

  //Modular: With known data, getNumOfTopo can be tested independently.
  int iVal, num2DElems = 0, num3DElems = 0;
  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_POINT, &iVal, &err);
  CHECK_ERR(err);
  if (iVal < 55) {
    qIssueLoadMessage(strFilename.c_str(), iBase_FAILURE);
    CHECK_ERR_CRITICAL(iBase_FAILURE);
  }
  a2iEntCount[iBase_VERTEX][iMesh_POINT] =
      a2iEntCount[iBase_VERTEX][iMesh_ALL_TOPOLOGIES] =
	  a2iEntCount[iBase_ALL_TYPES][iMesh_POINT] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_LINE_SEGMENT, &iVal, &err);
  CHECK_ERR(err);
  a2iEntCount[iBase_EDGE][iMesh_LINE_SEGMENT] =
      a2iEntCount[iBase_ALL_TYPES][iMesh_LINE_SEGMENT] =
	  a2iEntCount[iBase_EDGE][iMesh_ALL_TOPOLOGIES] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_TRIANGLE, &iVal, &err);
  CHECK_ERR(err);
  num2DElems += iVal;
  a2iEntCount[iBase_FACE][iMesh_TRIANGLE] =
      a2iEntCount[iBase_ALL_TYPES][iMesh_TRIANGLE] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_QUADRILATERAL, &iVal, &err);
  CHECK_ERR(err);
  num2DElems += iVal;
  a2iEntCount[iBase_FACE][iMesh_QUADRILATERAL] =
      a2iEntCount[iBase_ALL_TYPES][iMesh_QUADRILATERAL] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_POLYGON, &iVal, &err);
  CHECK_ERR(err);
  a2iEntCount[iBase_ALL_TYPES][iMesh_POLYGON] =
      a2iEntCount[iBase_FACE][iMesh_POLYGON] = iVal;

  // Total for faces
  a2iEntCount[iBase_FACE][iMesh_ALL_TOPOLOGIES] =
      a2iEntCount[iBase_FACE][iMesh_TRIANGLE]
	  + a2iEntCount[iBase_FACE][iMesh_QUADRILATERAL]
	  + a2iEntCount[iBase_FACE][iMesh_POLYGON];

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_TETRAHEDRON, &iVal, &err);
  CHECK_ERR(err);
  num3DElems += iVal;
  a2iEntCount[iBase_ALL_TYPES][iMesh_TETRAHEDRON] =
      a2iEntCount[iBase_REGION][iMesh_TETRAHEDRON] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_PYRAMID, &iVal, &err);
  CHECK_ERR(err);
  num3DElems += iVal;
  a2iEntCount[iBase_ALL_TYPES][iMesh_PYRAMID] =
      a2iEntCount[iBase_REGION][iMesh_PYRAMID] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_PRISM, &iVal, &err);
  CHECK_ERR(err);
  num3DElems += iVal;
  a2iEntCount[iBase_ALL_TYPES][iMesh_PRISM] =
      a2iEntCount[iBase_REGION][iMesh_PRISM] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_HEXAHEDRON, &iVal, &err);
  CHECK_ERR(err);
  num3DElems += iVal;
  a2iEntCount[iBase_ALL_TYPES][iMesh_HEXAHEDRON] =
      a2iEntCount[iBase_REGION][iMesh_HEXAHEDRON] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_SEPTAHEDRON, &iVal, &err);
  CHECK_ERR(err);
  a2iEntCount[iBase_ALL_TYPES][iMesh_SEPTAHEDRON] =
      a2iEntCount[iBase_REGION][iMesh_SEPTAHEDRON] = iVal;

  iMesh_getNumOfTopo(instance, pvRootSet, iMesh_POLYHEDRON, &iVal, &err);
  CHECK_ERR(err);
  a2iEntCount[iBase_ALL_TYPES][iMesh_POLYHEDRON] =
      a2iEntCount[iBase_REGION][iMesh_POLYHEDRON] = iVal;

  // Total for regions
  a2iEntCount[iBase_REGION][iMesh_ALL_TOPOLOGIES] =
      a2iEntCount[iBase_REGION][iMesh_TETRAHEDRON]
	  + a2iEntCount[iBase_REGION][iMesh_PYRAMID]
	  + a2iEntCount[iBase_REGION][iMesh_PRISM]
	  + a2iEntCount[iBase_REGION][iMesh_HEXAHEDRON]
	  + a2iEntCount[iBase_REGION][iMesh_SEPTAHEDRON]
	  + a2iEntCount[iBase_REGION][iMesh_POLYHEDRON];

  // Overall total
  a2iEntCount[iBase_ALL_TYPES][iMesh_ALL_TOPOLOGIES] =
      a2iEntCount[iBase_VERTEX][iMesh_ALL_TOPOLOGIES]
	  + a2iEntCount[iBase_EDGE][iMesh_ALL_TOPOLOGIES]
	  + a2iEntCount[iBase_FACE][iMesh_ALL_TOPOLOGIES]
	  + a2iEntCount[iBase_REGION][iMesh_ALL_TOPOLOGIES];

  int iTopoDim = iGeomDim;
  if (iGeomDim == 3 && a2iEntCount[iBase_REGION][iMesh_ALL_TOPOLOGIES] == 0)
    iTopoDim = 2;
  if (iTopoDim == 2 && num2DElems < 10) {
    qIssueLoadMessage(strFilename.c_str(), iBase_FAILURE);
    CHECK_ERR_CRITICAL(iBase_FAILURE);
  }
  if (iTopoDim == 3 && num3DElems < 10) {
    qIssueLoadMessage(strFilename.c_str(), iBase_FAILURE);
    CHECK_ERR_CRITICAL(iBase_FAILURE);
  }

  // Modular: Can also test those totals directly when there's data
  // available.

  // Modular: With data, can test getNumOfType directly.

  //@@@ Test loop
  // Modular: This loop can be done in isolation, because we'll know how
  // many to expect in advance.
  cout << " Global entity query tests" << endl;
  for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
      cout << "  Type: ";
      cout.width(9);
      cout << astrTypes[eEType] << " Topology: ";
      cout.width(14);
      cout << astrTopologies[eETopo];
      //@@@@ Get entities, type, topo, dim
      iBase_EntityHandle* aEntHandles = NULL;
      int iNumEnt, entity_handles_allocated = 0;
      iMesh_getEntities(instance, pvRootSet,
			static_cast<iBase_EntityType>(eEType),
			static_cast<iMesh_EntityTopology>(eETopo), &aEntHandles,
			&entity_handles_allocated, &iNumEnt, &err);
      switch (err)
	{
	case iBase_SUCCESS:
	  if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	    qOK = (iNumEnt == a2iEntCount[eEType][eETopo]);
	    cout << " OK" << endl;
	    iNPassed++;
	  }
	  else {
	    cout << " succeeded when it should have thrown an error." << endl;
	    iNFailed++;
	  }
	  break;
	case iBase_BAD_TYPE_AND_TOPO:
	  if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	    // This should only occur if a precondition test failed.
	    cout << endl;
	    cout << "Threw an error instead of returning! " << endl;
	    writeErrorDescr(err);
	    iNFailed++;
	  }
	  else {
	    cout << " threw an error, as expected." << endl;
	    iNPassed++;
	  }
	  break;
	case iBase_NOT_SUPPORTED:
	  // This should never occur, according to the spec.
	  cout << "; feasible though unsupported combination" << endl;
	  writeErrorDescr(err);
	  iNFailed++;
	  break;
	default:
	  cout << endl;
	  cout << "Unexpected error! " << endl;
	  writeErrorDescr(err);
	  iNFailed++;
	  break;
	}
      UT_FREE(aEntHandles);
    }
  } // Done with double loop to verify consistency of number of entities
    // fetched and number of entities that the mesh claims to have.

    //@@ getVtxArrCoords should insist on vertex handles
  cout << " Verifying that getVtxArrCoords requires vertex handles as input..."
      << endl;

  iBase_EntityHandle *aEntHandles = NULL;
  int iNumEnt, iEntHandles_allocated = 0;
  iMesh_getEntities(instance, pvRootSet, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
		    &aEntHandles, &iEntHandles_allocated, &iNumEnt, &err);
  CHECK_ERR(err);

  double *adCoords;
  int SO = iBase_BLOCKED;
  int coords_size, coords_allocated = 0;
  iMesh_getVtxArrCoords(instance, aEntHandles, iNumEnt, SO, &adCoords,
			&coords_allocated, &coords_size, &err);
  CHECK_ERR2(err, iBase_INVALID_ENTITY_TYPE);
  // The following should cause a break, unless the mesh happens to
  // contain only vertices.
  if (!err) {
    // Should only get here if all the entities in the entire mesh
    // happened to be vertices.
    int iNumVerts;
    iMesh_getNumOfType(instance, pvRootSet, iBase_VERTEX, &iNumVerts, &err);
    CHECK_ERR(err);
    TESTEQ(iNumVerts, iNumEnt);
    TEST(coords_size == iNumVerts * 3);
  }
  UT_FREE(aEntHandles);
  UT_FREE(adCoords);
}
