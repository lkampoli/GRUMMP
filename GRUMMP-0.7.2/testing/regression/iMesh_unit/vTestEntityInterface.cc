#include "iMesh_unitTest.hh"

// Test stuff specific to interface Entity.  Those functions are:
//    initEntIter, getNextEntIter, resetEntIter, endEntIter,
//    getEntTopo, getEntType, getVtxCoord, getEntAdj

void
vTestEntityInterface(iMesh_Instance& instance)
{
  // getEntTopo, getEntType, getVtxCoord, getEntAdj can all be tested
  // via consistency with the global calls.  For instance, all entities
  // retrieved by the global getEntities(FACE, TRIANGLE) call should be
  // of type face and topology triangle.

  int err;
  bool qOK;
  //@@@ Test loop for getEntTopo and getEntType.
  cout << " Verifying retrieved entity type and topology" << endl;
  int eEType, eETopo;
  for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
      qOK = true;
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	cout << "  Type: ";
	cout.width(9);
	cout << astrTypes[eEType] << " Topology: ";
	cout.width(14);
	cout << astrTopologies[eETopo] << endl;
	iBase_EntityHandle *aEntHandles = NULL;
	int iNumEnt, iEntHandlesAlloc = 0;

	// This call must succeed, or the Mesh tests wouldn't all pass.
	iMesh_getEntities(instance, pvRootSet,
			  static_cast<iBase_EntityType>(eEType),
			  static_cast<iMesh_EntityTopology>(eETopo),
			  &aEntHandles, &iEntHandlesAlloc, &iNumEnt, &err);
	CHECK_ERR(err);

	if (err == iBase_SUCCESS) {
	  for (int iEnt = 0; iEnt < iNumEnt; iEnt++) {
	    iBase_EntityHandle pEnt = aEntHandles[(iEnt)];
	    int ETy, ETo;

	    iMesh_getEntType(instance, pEnt, &ETy, &err);
	    CHECK_ERR(err);
	    iMesh_getEntTopo(instance, pEnt, &ETo, &err);
	    CHECK_ERR(err);

	    qOK = qOK
		&& ((ETy == eEType || eEType == iBase_ALL_TYPES)
		    && (ETo == eETopo || eETopo == iMesh_ALL_TOPOLOGIES));
	  }
	}
	UT_FREE(aEntHandles);
	TEST(qOK);
      }
    }
  } // Done with double loop to verify consistency of entity type and
    // topology.

    // Now test getVtxCoord for consistency with the getEntities(VERTEX)
    // -> getVtxArrCoord chain.
  cout << " Checking single vertex coordinate retrieval" << endl;
  {
    iBase_EntityHandle *aVtxHandles = NULL;
    double *adAllCoords = NULL;
    double adCoords[3];

    int iGeomDim;
    iMesh_getGeometricDimension(instance, &iGeomDim, &err);
    CHECK_ERR(err);

    // Get all the vertex coordinates.
    int iVtxHandlesAlloc = 0, iCoordsAlloc = 0;
    int iNVerts, iNCoords;
    iMesh_getEntities(instance, pvRootSet, iBase_VERTEX, iMesh_POINT,
		      &aVtxHandles, &iVtxHandlesAlloc, &iNVerts, &err);
    CHECK_ERR(err);
    int SO = iBase_INTERLEAVED;
    iMesh_getVtxArrCoords(instance, aVtxHandles, iNVerts, SO, &adAllCoords,
			  &iCoordsAlloc, &iNCoords, &err);
    CHECK_ERR(err);
    TEST(iNCoords == iGeomDim * iNVerts);

    // Now get them vertex by vertex and compare.
    qOK = true;
    for (int iV = 0, iC = 0; iV < iNVerts; iV++) {
      iBase_EntityHandle vert_handle = aVtxHandles[(iV)];
      iMesh_getVtxCoord(instance, vert_handle, &adCoords[0], &adCoords[1],
			&adCoords[2], &err);
      CHECK_ERR(err);
      for (int ii = 0; ii < iGeomDim; iC++, ii++) {
	qOK = qOK && (adAllCoords[(iC)] == adCoords[ii]);
      }
      if (iGeomDim == 2) {
	qOK = qOK && (adCoords[2] == 0);
      }
    }
    UT_FREE(adAllCoords);
    UT_FREE(aVtxHandles);
    TEST(qOK);
  }

  //@@ Test single entity adjacency call for consistency with the global
  // call.

  //@@@ Test loop
  cout << " Checking single entity adjacency calls" << endl;
  for (eEType = iBase_VERTEX; eEType < iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo < iMesh_ALL_TOPOLOGIES; eETopo++) {
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	qOK = true;
	cout << "  Type: ";
	cout.width(9);
	cout << astrTypes[eEType] << " Topology: ";
	cout.width(14);
	cout << astrTopologies[eETopo] << endl;

	// Now test for equivalence in adjacency operations.
	int eAdjType, eLastType = iBase_REGION;
	for (eAdjType = iBase_VERTEX; eAdjType <= eLastType; eAdjType++) {
	  // Skip any that are identities or unsupported.
	  if (eAdjType == eEType
	      || (a2iAdjTable[eEType][eAdjType] == iBase_UNAVAILABLE))
	    continue;
	  cout << "   Adjacent entity type:  ";
	  cout.width(9);
	  cout << astrTypes[eAdjType] << endl;

	  iBase_EntityHandle *aEnts = NULL, *aAdjEnt = NULL;
	  int *aiOffset = NULL, *aiIndex = NULL;
	  int iNumAdj, iOffsetSize, iIndexSize, iEntsSize;
	  int iAdjAlloc = 0, iOffsetAlloc = 0, iEntsAlloc = 0, iIndexAlloc = 0;
	  // Get adjacencies using the entityset call.

	  iMesh_getAdjEntIndices(instance, pvRootSet,
				 static_cast<iBase_EntityType>(eEType),
				 static_cast<iMesh_EntityTopology>(eETopo),
				 static_cast<iBase_EntityType>(eAdjType),
				 &aEnts, &iEntsAlloc, &iEntsSize, &aAdjEnt,
				 &iAdjAlloc, &iNumAdj, &aiIndex, &iIndexAlloc,
				 &iIndexSize, &aiOffset, &iOffsetAlloc,
				 &iOffsetSize, &err);
	  CHECK_ERR(err);

	  // Verify that the size of the output arrays are correct.
	  qOK = (iEntsSize == iOffsetSize - 1 && err == iBase_SUCCESS);

	  int iEnt, iOff;
	  for (iEnt = 0, iOff = 0; iEnt < iEntsSize && qOK; iEnt++, iOff++) {
	    iBase_EntityHandle handle = aEnts[iEnt];
	    int iNLocalAdj, iAdjEntHandlesAlloc = 0;
	    iBase_EntityHandle *aAdjEntHandles = NULL;
	    iMesh_getEntAdj(instance, handle,
			    static_cast<iBase_EntityType>(eAdjType),
			    &aAdjEntHandles, &iAdjEntHandlesAlloc, &iNLocalAdj,
			    &err);
	    CHECK_ERR(err);
	    // Number of adjacencies must match.
	    qOK = (iNLocalAdj == (aiOffset[(iOff + 1)] - aiOffset[(iOff)]));
	    // Now check that the adjacencies match.
	    int iOffset = aiOffset[iOff];
	    for (int iAdj = 0; iAdj < iNLocalAdj; iAdj++) {
	      iBase_EntityHandle global_handle =
		  aAdjEnt[aiIndex[iOffset + iAdj]];
	      iBase_EntityHandle local_handle = aAdjEntHandles[iAdj];
	      qOK = qOK && (global_handle == local_handle);
	    } // Loop over local adjacencies
	    UT_FREE(aAdjEntHandles);
	  } // Loop over entities
	  UT_FREE(aEnts);
	  UT_FREE(aAdjEnt);
	  UT_FREE(aiOffset);
	  UT_FREE(aiIndex);
	} // Loop over adjacency types
	TEST(qOK);
      }
    }
  } // Done with double loop to verify consistency of adjacencies

  // Test the scalar second adjacency call.
  vSecondAdjacencyTest(instance);

  // Iterators should be tested by iterating over every legal type and
  // topology combination to be sure that all the entities retrieved by
  // the global call are traversed by the iterator.  To make sure reset
  // works properly, run through each collection twice.
  cout << " Testing single entity iterators" << endl;
  for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
      // First, grab everything and put it in a set.
      int iMaxCount = a2iEntCount[eEType][eETopo];
      iBase_EntityHandle *aEnt = (iBase_EntityHandle*) calloc(
	  iMaxCount, sizeof(iBase_EntityHandle));
      int iEntAlloc = iMaxCount;
      int iActualCount = 0;
//       if (a2qLegalTypeAndTopo[eEType][eETopo]) {
      cout << "  Type: ";
      cout.width(9);
      cout << astrTypes[eEType] << " Topology: ";
      cout.width(14);
      cout << astrTopologies[eETopo] << endl;
      iMesh_getEntities(instance, pvRootSet,
			static_cast<iBase_EntityType>(eEType),
			static_cast<iMesh_EntityTopology>(eETopo), &aEnt,
			&iEntAlloc, &iActualCount, &err);
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	CHECK_ERR(err);
      }

      // Should see what happens with zero-sized iterators.
      // if (iActualCount == 0) continue;

      // Put all these entities into a big set.  This will need to be a
      // multiset when testing entity set queries, because the entity
      // set could have the same entity in it twice.
      std::set<iBase_EntityHandle> sEnt;
      sEnt.insert(aEnt, aEnt + iActualCount);
      std::set<iBase_EntityHandle> sEntCopy(sEnt);

      qOK = ((int) (sEnt.size()) == iActualCount);

      // Now iterate over these entities, and remove them one by one
      // from the set, hopefully succeeding every time.
      iBase_EntityIterator entIter;
      int qHasData = true;

      do {
	iMesh_initEntIter(instance, pvRootSet,
			  static_cast<iBase_EntityType>(eEType),
			  static_cast<iMesh_EntityTopology>(eETopo), 0,
			  &entIter, &err);
	if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	  BREAK_ON_ERR(err);
	}
	else {
	  break;
	}
	int iMissedCount = 0;
	int iExtraCount = 0;
	iBase_EntityHandle entHandle;
	iMesh_getNextEntIter(instance, entIter, &entHandle, &qHasData, &err);
	BREAK_ON_ERR(err);
	while (qHasData) {
	  int iRes = sEnt.erase(entHandle);
	  if (iRes > 1)
	    iExtraCount++;
	  else if (iRes == 0)
	    iMissedCount++;
	  iMesh_getNextEntIter(instance, entIter, &entHandle, &qHasData, &err);
	  BREAK_ON_ERR(err);
	}
	qOK = qOK && (iMissedCount == 0) && (iExtraCount == 0)
	    && (sEnt.size() == 0);

	iMesh_resetEntIter(instance, entIter, &err);
	BREAK_ON_ERR(err);

	// Now repeat to be sure the reset worked.
	iMissedCount = 0;
	iExtraCount = 0;
	iMesh_getNextEntIter(instance, entIter, &entHandle, &qHasData, &err);
	BREAK_ON_ERR(err);
	while (qHasData) {
	  int iRes = sEntCopy.erase(entHandle);
	  if (iRes > 1)
	    iExtraCount++;
	  else if (iRes == 0)
	    iMissedCount++;
	  iMesh_getNextEntIter(instance, entIter, &entHandle, &qHasData, &err);
	  BREAK_ON_ERR(err);
	}
	qOK = qOK && (iMissedCount == 0) && (iExtraCount == 0)
	    && (sEntCopy.size() == 0);

	iMesh_endEntIter(instance, entIter, &err);
	BREAK_ON_ERR(err);
      }while (0);

      switch (err)
	{
	case iBase_FAILURE:
	  // Could mean there aren't any of these; check for this.
	  if (iMaxCount == 0)
	    iNPassed++;
	  else {
	    cout << "Max count should be zero; is " << iMaxCount << endl;
	    qAllPassed = false;
	    iNFailed++;
	  }
	  break;
	case iBase_BAD_TYPE_AND_TOPO:
	  {
	    if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	      cout << "Didn't produce iterator for valid combinations." << endl;
	      qAllPassed = false;
	      iNFailed++;
	    }
	    else {
	      iNPassed++;
	    }
	  }
	  break;
	case iBase_INVALID_ARGUMENT:
	  // Can happen for 2D meshes requested to iterate over
	  // regions.
	  {
	    int iGeomDim;
	    iMesh_getGeometricDimension(instance, &iGeomDim, &err);
	    CHECK_ERR(err);
	    if (iGeomDim == 2
		&& ((eEType == iBase_REGION)
		    || ((eEType == iBase_ALL_TYPES)
			&& ((eETopo == iMesh_POLYHEDRON)
			    || (eETopo == iMesh_TETRAHEDRON)
			    || (eETopo == iMesh_PYRAMID)
			    || (eETopo == iMesh_PRISM)
			    || (eETopo == iMesh_HEXAHEDRON)
			    || (eETopo == iMesh_SEPTAHEDRON)))))
	      iNPassed++;
	    else {
	      cout << "Valid arguments incorrectly diagnosed." << endl;
	      qAllPassed = false;
	      iNFailed++;
	    }
	  }
	  break;
	case iBase_SUCCESS:
	  break;
	default:
	  cout << "Unexpected error type!" << endl;
	  TEST(false);
	  break;
	} // Done with catching exceptions
//       } // This was a legal type/topo combination
      UT_FREE(aEnt);
    }
  } // Done with double loop to test entity iterators.

  // Test reciprocity of entity adjacency calls.
  cout << " Testing single entity adjacency reciprocity" << endl;
  // Edge->vertex
  if (a2iAdjTable[1][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][1] != iBase_UNAVAILABLE) {
    cout << "  Edges and verts:";
    vEntityReciprocityTest(instance, iBase_EDGE, iBase_VERTEX);
    vEntityReciprocityTest(instance, iBase_VERTEX, iBase_EDGE);
    cout << endl;
  }

  // Face->vertex
  if (a2iAdjTable[2][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][2] != iBase_UNAVAILABLE) {
    cout << "  Faces and verts:";
    vEntityReciprocityTest(instance, iBase_FACE, iBase_VERTEX);
    vEntityReciprocityTest(instance, iBase_VERTEX, iBase_FACE);
    cout << endl;
  }

  // Region->vertex
  if (a2iAdjTable[3][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and verts:";
    vEntityReciprocityTest(instance, iBase_REGION, iBase_VERTEX);
    vEntityReciprocityTest(instance, iBase_VERTEX, iBase_REGION);
    cout << endl;
  }

  // Face->edge
  if (a2iAdjTable[2][1] != iBase_UNAVAILABLE
      && a2iAdjTable[1][2] != iBase_UNAVAILABLE) {
    cout << "  Faces and edges:";
    vEntityReciprocityTest(instance, iBase_FACE, iBase_EDGE);
    vEntityReciprocityTest(instance, iBase_EDGE, iBase_FACE);
    cout << endl;
  }

  // Region->edge
  if (a2iAdjTable[3][1] != iBase_UNAVAILABLE
      && a2iAdjTable[1][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and edges:";
    vEntityReciprocityTest(instance, iBase_REGION, iBase_EDGE);
    vEntityReciprocityTest(instance, iBase_EDGE, iBase_REGION);
    cout << endl;
  }

  // Region->face
  if (a2iAdjTable[3][2] != iBase_UNAVAILABLE
      && a2iAdjTable[2][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and faces:";
    vEntityReciprocityTest(instance, iBase_REGION, iBase_FACE);
    vEntityReciprocityTest(instance, iBase_FACE, iBase_REGION);
    cout << endl;
  }

  cout << "Entity adjacency identity test..." << endl;
  //@@ Ensure that identity adjacency calls return nothing.
  // Test to be sure that no VERTEX has an adjacent VERTEX, etc.  The
  // iMesh data model specifies that entities of the same type are never
  // adjacent, but are connected instead through entities of different
  // type incident on both.
  for (int iType = 0; iType <= 3; iType++) {
    if (a2iAdjTable[iType][iType] != iBase_UNAVAILABLE) {
      cout << "  Type: ";
      cout.width(9);
      cout << astrTypes[iType] << endl;
      iBase_EntityType eType = static_cast<iBase_EntityType>(iType);

      iBase_EntityHandle * aEntHand = NULL;
      int iNumEnt, iEntHandAlloc = 0;
      iMesh_getEntities(instance, pvRootSet, eType, iMesh_ALL_TOPOLOGIES,
			&aEntHand, &iEntHandAlloc, &iNumEnt, &err);
      CHECK_ERR(err);

      iBase_EntityHandle * aAdjHand = NULL;
      int iAdjHandAlloc = 0;
      qOK = true;
      for (int i = 0; i < iNumEnt; i++) {
	iBase_EntityHandle entity = aEntHand[i];
	int iNAdj;
	iMesh_getEntAdj(instance, entity, eType, &aAdjHand, &iAdjHandAlloc,
			&iNAdj, &err);
	CHECK_ERR(err);
	qOK = qOK && (iNAdj == 0);
      }
      UT_FREE(aAdjHand);
      UT_FREE(aEntHand);
      TEST(qOK);
    }
  }
}
