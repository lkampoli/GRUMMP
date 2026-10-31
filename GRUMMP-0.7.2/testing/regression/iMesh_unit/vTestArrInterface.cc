#include "iMesh_unitTest.hh"

void
vTestArrInterface(iMesh_Instance& instance)
{
  bool qOK;
  int eEType, eETopo, err;

  //@@@ Test loop
  cout << " Verifying retrieved entity type and topology (array version)"
      << endl << "    and global entity adjacency" << endl;
  for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	cout << "  Type: ";
	cout.width(9);
	cout << astrTypes[eEType] << " Topology: ";
	cout.width(14);
	cout << astrTopologies[eETopo] << endl;
	//@@@@ Get entities, type, topo
	iBase_EntityHandle * aEntHandles = NULL;
	int iNumEnt, iEntHandlesAlloc = 0;
	iMesh_getEntities(instance, pvRootSet,
			  static_cast<iBase_EntityType>(eEType),
			  static_cast<iMesh_EntityTopology>(eETopo),
			  &aEntHandles, &iEntHandlesAlloc, &iNumEnt, &err);
	CHECK_ERR(err);

	if (err == iBase_SUCCESS) {
	  int * aiEType = NULL;
	  int * aiETopo = NULL;

	  int topo_size, type_size, topo_alloc = 0, type_alloc = 0;
	  iMesh_getEntArrTopo(instance, aEntHandles, iNumEnt, &aiETopo,
			      &topo_alloc, &topo_size, &err);
	  CHECK_ERR(err);
	  iMesh_getEntArrType(instance, aEntHandles, iNumEnt, &aiEType,
			      &type_alloc, &type_size, &err);
	  CHECK_ERR(err);

	  qOK = (topo_size == iNumEnt) && (type_size == iNumEnt);

	  int iTopo, iType;
	  for (iTopo = 0, iType = 0; iTopo < topo_size; iTopo++, iType++) {
	    int eTo = aiETopo[(iType)];
	    int eTy = aiEType[(iType)];
	    qOK = qOK
		&& ((eTo == eETopo || eETopo == iMesh_ALL_TOPOLOGIES)
		    && (eTy == eEType || eEType == iBase_ALL_TYPES));
	  }
	  UT_FREE(aiEType);
	  UT_FREE(aiETopo);
	  assert(iTopo <= topo_size);
	  assert(iType <= type_size);
	}
	UT_FREE(aEntHandles);
	if (!qOK) {
	  cout << endl << "Inconsistency in entity type or topology.\n\n";
	  qAllPassed = false;
	  iNFailed++;
	}

	TEST(qOK);
	if (eEType == iBase_ALL_TYPES)
	  continue;

	qOK = true;
	//@@@@ Equivalence tests for entityGetAdj and entitysetGetAdj
	// Now test for equivalence in adjacency operations.
	int eAdjType, eLastType = iBase_REGION;
	for (eAdjType = iBase_VERTEX; eAdjType <= eLastType; eAdjType++) {
	  // Skip any that are identities or unsupported.
	  if (eAdjType == eEType || a2iAdjTable[eEType][eAdjType] == 0)
	    continue;
	  cout << "   Adjacent type: ";
	  cout.width(9);
	  cout << astrTypes[eAdjType] << endl;

	  // Get adjacencies using the entityset call.
	  iBase_EntityHandle *aEnts2 = NULL, *aAdjEnt2 = NULL;
	  int *aiIndex = NULL, *aiOffset2 = NULL;
	  int iNumAdj2, iOffset2Size, iIndexSize, iEnts2Size;
	  int iAdjEnt2Alloc = 0, iOffset2Alloc = 0, iIndexAlloc = 0,
	      iEnts2Alloc = 0;
	  iMesh_getAdjEntIndices(instance, pvRootSet,
				 static_cast<iBase_EntityType>(eEType),
				 static_cast<iMesh_EntityTopology>(eETopo),
				 static_cast<iBase_EntityType>(eAdjType),
				 &aEnts2, &iEnts2Alloc, &iEnts2Size, &aAdjEnt2,
				 &iAdjEnt2Alloc, &iNumAdj2, &aiIndex,
				 &iIndexAlloc, &iIndexSize, &aiOffset2,
				 &iOffset2Alloc, &iOffset2Size, &err);
	  CHECK_ERR(err);

	  if (err == iBase_SUCCESS) {
	    int *aiOffset1;
	    iBase_EntityHandle *aAdjEnt1 = NULL;
	    // Get adjacencies using entityGetAdjacencies.
	    int iNumAdj1, iOffset1Size;
	    int iAdjEnt1Alloc = 0, iOffset1Alloc = 0;
	    iMesh_getEntArrAdj(instance, aEnts2, iEnts2Size,
			       static_cast<iBase_EntityType>(eAdjType),
			       &aAdjEnt1, &iAdjEnt1Alloc, &iNumAdj1, &aiOffset1,
			       &iOffset1Alloc, &iOffset1Size, &err);
	    CHECK_ERR(err);

	    // Verify that the size of the output arrays are correct.
	    qOK = (iNumAdj1 >= iNumAdj2) && (iNumAdj1 == iIndexSize)
		&& (iNumEnt == iEnts2Size) && (iNumEnt == iOffset1Size - 1)
		&& (iNumEnt == iOffset2Size - 1);

	    // Verify that the contents are the same.
	    // Loop first over entities
	    int iAdj = 0;
	    for (int iEnt = 0; iEnt < iNumEnt; iEnt++) {
	      qOK = (aiOffset1[iEnt + 1] == aiOffset2[iEnt + 1]);
	      for (; iAdj < aiOffset1[iEnt + 1]; iAdj++) {
		iBase_EntityHandle entA = aAdjEnt1[iAdj];
		iBase_EntityHandle entB = aAdjEnt2[aiIndex[iAdj]];
		qOK = (entA == entB);
	      }
	    }
	    assert(!qOK || iAdj == iNumAdj1);
	    assert(!qOK || iAdj == iIndexSize);
	    UT_FREE(aAdjEnt1);
	    UT_FREE(aiOffset1);
	  }
	  UT_FREE(aEnts2);
	  UT_FREE(aAdjEnt2);
	  UT_FREE(aiIndex);
	  UT_FREE(aiOffset2);
	}
	TEST(qOK);
      }
    }
  } // Done with double loop to verify consistency of bulk adjacency calls.

  // Test the scalar second adjacency call.
  vSecondAdjacencyArrayTest(instance);

  // Iterators should be tested by iterating over every legal type and
  // topology combination to be sure that all the entities retrieved by
  // the global call are traversed by the iterator.  To make sure reset
  // works properly, run through each collection twice.
  cout << " Testing array iterators" << endl;
  for (eEType = iBase_VERTEX; eEType <= iBase_ALL_TYPES; eEType++) {
    for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
      cout << "  Type: ";
      cout.width(9);
      cout << astrTypes[eEType] << " Topology: ";
      cout.width(14);
      cout << astrTopologies[eETopo] << endl;
      // First, grab everything and put it in a set.
      int iMaxCount = a2iEntCount[eEType][eETopo];
      iBase_EntityHandle *aEnt = (iBase_EntityHandle*) calloc(
	  iMaxCount, sizeof(iBase_EntityHandle));
      int iEntAlloc = iMaxCount;
      int iActualCount;

      std::set<iBase_EntityHandle> sEnt;
      iMesh_getEntities(instance, pvRootSet,
			static_cast<iBase_EntityType>(eEType),
			static_cast<iMesh_EntityTopology>(eETopo), &aEnt,
			&iEntAlloc, &iActualCount, &err);
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	CHECK_ERR(err);

	// Put all these entities into a big set.  This will need to be a
	// multiset when testing entity set queries, because the entity
	// set could have the same entity in it twice.
	sEnt.insert(aEnt, aEnt + iActualCount);
	qOK = ((int) (sEnt.size()) == iActualCount);
      }
      else
	qOK = true;
      std::set<iBase_EntityHandle> sEntCopy(sEnt);

      // Now iterate over these entities, and remove them one by one
      // from the set, hopefully succeeding every time.
      iBase_EntityArrIterator arrIter;
      iBase_EntityHandle *aEntHandles = (iBase_EntityHandle*) calloc(
	  100, sizeof(iBase_EntityHandle*));
      int iEntHandlesAlloc = 100;
      int qHasData = true;

      do {
	iMesh_initEntArrIter(instance, pvRootSet,
			     static_cast<iBase_EntityType>(eEType),
			     static_cast<iMesh_EntityTopology>(eETopo), 100, 0,
			     &arrIter, &err);
	if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	  BREAK_ON_ERR(err);
	}
	else {
	  break;
	}
	int iMissedCount = 0;
	int iExtraCount = 0;
	int iNumHandles;
	iMesh_getNextEntArrIter(instance, arrIter, &aEntHandles,
				&iEntHandlesAlloc, &iNumHandles, &qHasData,
				&err);
	BREAK_ON_ERR(err);
	while (qHasData) {
	  for (int ii = 0; ii < iNumHandles; ii++) {
	    int iRes = sEnt.erase(aEntHandles[ii]);
	    if (iRes > 1)
	      iExtraCount++;
	    else if (iRes == 0)
	      iMissedCount++;
	  }
	  iMesh_getNextEntArrIter(instance, arrIter, &aEntHandles,
				  &iEntHandlesAlloc, &iNumHandles, &qHasData,
				  &err);
	  BREAK_ON_ERR(err);
	}
	qOK = qOK && (iMissedCount == 0) && (iExtraCount == 0)
	    && (sEnt.size() == 0);

	iMesh_resetEntArrIter(instance, arrIter, &err);
	BREAK_ON_ERR(err);

	// Now repeat to be sure the reset worked.
	iMissedCount = 0;
	iExtraCount = 0;
	iMesh_getNextEntArrIter(instance, arrIter, &aEntHandles,
				&iEntHandlesAlloc, &iNumHandles, &qHasData,
				&err);
	BREAK_ON_ERR(err);
	while (qHasData) {
	  for (int ii = 0; ii < iNumHandles; ii++) {
	    int iRes = sEntCopy.erase(aEntHandles[ii]);
	    if (iRes > 1)
	      iExtraCount++;
	    else if (iRes == 0)
	      iMissedCount++;
	  }
	  iMesh_getNextEntArrIter(instance, arrIter, &aEntHandles,
				  &iEntHandlesAlloc, &iNumHandles, &qHasData,
				  &err);
	  BREAK_ON_ERR(err);
	}
	qOK = qOK && (iMissedCount == 0) && (iExtraCount == 0)
	    && (sEntCopy.size() == 0);

	iMesh_endEntArrIter(instance, arrIter, &err);
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
      UT_FREE(aEntHandles);
      UT_FREE(aEnt);
    }
  } // Done with double loop to test entity iterators.

  // Test reciprocity of array adjacency calls.
  // Edge->vertex
  cout << " Testing entity array adjacency reciprocity" << endl;
  if (a2iAdjTable[1][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][1] != iBase_UNAVAILABLE) {
    cout << "  Edges and verts:";
    vArrayReciprocityTest(instance, iBase_EDGE, iBase_VERTEX);
    vArrayReciprocityTest(instance, iBase_VERTEX, iBase_EDGE);
    cout << endl;
  }

  // Face->vertex
  if (a2iAdjTable[2][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][2] != iBase_UNAVAILABLE) {
    cout << "  Faces and verts:";
    vArrayReciprocityTest(instance, iBase_FACE, iBase_VERTEX);
    vArrayReciprocityTest(instance, iBase_VERTEX, iBase_FACE);
    cout << endl;
  }

  // Region->vertex
  if (a2iAdjTable[3][0] != iBase_UNAVAILABLE
      && a2iAdjTable[0][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and verts:";
    vArrayReciprocityTest(instance, iBase_REGION, iBase_VERTEX);
    vArrayReciprocityTest(instance, iBase_VERTEX, iBase_REGION);
    cout << endl;
  }

  // Face->edge
  if (a2iAdjTable[2][1] != iBase_UNAVAILABLE
      && a2iAdjTable[1][2] != iBase_UNAVAILABLE) {
    cout << "  Faces and edges:";
    vArrayReciprocityTest(instance, iBase_FACE, iBase_EDGE);
    vArrayReciprocityTest(instance, iBase_EDGE, iBase_FACE);
    cout << endl;
  }

  // Region->edge
  if (a2iAdjTable[3][1] != iBase_UNAVAILABLE
      && a2iAdjTable[1][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and edges:";
    vArrayReciprocityTest(instance, iBase_REGION, iBase_EDGE);
    vArrayReciprocityTest(instance, iBase_EDGE, iBase_REGION);
    cout << endl;
  }

  // Region->face
  if (a2iAdjTable[3][2] != iBase_UNAVAILABLE
      && a2iAdjTable[2][3] != iBase_UNAVAILABLE) {
    cout << "  Regions and faces:";
    vArrayReciprocityTest(instance, iBase_REGION, iBase_FACE);
    vArrayReciprocityTest(instance, iBase_FACE, iBase_REGION);
    cout << endl;
  }

  cout << " Array adjacency identity test..." << endl;
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
      int * aiOffset = NULL;
      int iNAdj, iNOff;
      int iAdjHandAlloc = 0;
      int iOffsetAlloc = 0;
      iMesh_getEntArrAdj(instance, aEntHand, iNumEnt, eType, &aAdjHand,
			 &iAdjHandAlloc, &iNAdj, &aiOffset, &iOffsetAlloc,
			 &iNOff, &err);
      CHECK_ERR(err);
      UT_FREE(aAdjHand);
      UT_FREE(aiOffset);
      UT_FREE(aEntHand);
      TEST(iNAdj == 0 && iNOff == iNumEnt + 1);
    }
  }
}
