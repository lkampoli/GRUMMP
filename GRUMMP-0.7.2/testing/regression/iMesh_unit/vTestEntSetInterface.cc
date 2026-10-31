#include "iMesh_unitTest.hh"

void
vTestEntSetInterface(iMesh_Instance& instance)
{
  int err;
  //     createEntSet, destroyEntSet, isList, getNumEntSets, getEntSets,
  //     addEntToSet, rmvEntFromSet, addEntArrToSet, rmvEntArrFromSet,
  //     addEntSet, rmvEntSet, isContainedIn
  iBase_EntitySetHandle apvEntSets[11];
  int iNumSets, iOrigNumSets;
  bool qGotAllSets;
  int i;

  // Find out how many sets there were originally 
  iMesh_getNumEntSets(instance, pvRootSet, -1, &iOrigNumSets, &err);
  // Make a place to store retrieved sets
  int iSetsAlloc = 10 + iOrigNumSets;
  iBase_EntitySetHandle *aSets = (iBase_EntitySetHandle*) calloc(
      iSetsAlloc, sizeof(iBase_EntitySetHandle));

  cout << "Testing creation, basic methods, deletion...              ";
  // Try creating a number of sets and lists and see if they show up as members
  // of the root set...
  bool qOK = true;

  for (i = 0; i < 5; i++) {
    iMesh_createEntSet(instance, true, &apvEntSets[2 * i], &err);
    CHECK_ERR(err);
    if (err)
      qOK = false;
    int list;
    iMesh_isList(instance, apvEntSets[2 * i], &list, &err);
    CHECK_ERR(err);
    if (err || !list)
      qOK = false;
  }
  for (i = 0; i < 5; i++) {
    iMesh_createEntSet(instance, false, &apvEntSets[2 * i + 1], &err);
    CHECK_ERR(err);
    if (err)
      qOK = false;
    int list;
    iMesh_isList(instance, apvEntSets[2 * i + 1], &list, &err);
    CHECK_ERR(err);
    if (err || list)
      qOK = false;
  }
  TEST(qOK);

  // Whether the number of hops is -1, 0, 1 is irrelevant for the root set,
  // especially since none of the sets contain each other
  qGotAllSets = true;

  iMesh_getNumEntSets(instance, pvRootSet, -1, &iNumSets, &err);
  CHECK_ERR(err);
  qGotAllSets = qGotAllSets && ((iNumSets - iOrigNumSets) == 10);

  iMesh_getNumEntSets(instance, pvRootSet, 0, &iNumSets, &err);
  CHECK_ERR(err);
  qGotAllSets = qGotAllSets && ((iNumSets - iOrigNumSets) == 10);

  iMesh_getNumEntSets(instance, pvRootSet, 1, &iNumSets, &err);
  CHECK_ERR(err);
  qGotAllSets = qGotAllSets && ((iNumSets - iOrigNumSets) == 10);

  TEST(qGotAllSets);

  // Try to get all the entity sets from the root set
  iMesh_getEntSets(instance, pvRootSet, 0, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
  CHECK_ERR(err);
  qGotAllSets = false;
  if ((iNumSets - iOrigNumSets) == 10)
    GOT_ALL_ITEMS(apvEntSets, aSets, iNumSets, qGotAllSets);
  TEST(qGotAllSets);

  // Try to destroy all the sets
  for (i = 0; i < 10; i++) {
    iMesh_destroyEntSet(instance, apvEntSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing addition and removal for sets in sets...          ";

  // Create 10 sets, the even ones of which are lists ( irrelevant as far
  // as sets containing one another, but good for thoroughness )
  for (i = 0; i < 10; i++) {
    iMesh_createEntSet(instance, (i % 2 == 0), &apvEntSets[i], &err);
    CHECK_ERR(err);
  }

  // Create a hierarchy of sets where            [0]
  //     set 0 contains sets 1-3,           ______|_____
  //     set 1 contains sets 4-5,           [1], [2], [3]
  //     set 3 contains sets 6-7,         ___|___   ___|___
  // and set 4 contains sets 8-9          [4],[5]   [6],[7]
  //                                              ___|___
  // (and test the contain query)                 [8],[9]

  qGotAllSets = true;
  int iIsContained;
  for (i = 1; i <= 3; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[0], &err);
    CHECK_ERR(err);
    iMesh_isEntSetContained(instance, apvEntSets[0], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && iIsContained;
  }
  for (i = 4; i <= 5; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[1], &err);
    CHECK_ERR(err);
    iMesh_isEntSetContained(instance, apvEntSets[1], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && iIsContained;
  }
  for (i = 6; i <= 7; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[3], &err);
    CHECK_ERR(err);
    iMesh_isEntSetContained(instance, apvEntSets[3], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && iIsContained;
  }
  for (i = 8; i <= 9; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[6], &err);
    CHECK_ERR(err);
    iMesh_isEntSetContained(instance, apvEntSets[6], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && iIsContained;
  }
  TEST(qGotAllSets);
  vDoSaveLoad(&instance, &pvRootSet);

  if (qSaveLoad) {
    vRecoverSetsContain(instance, pvRootSet, apvEntSets);
  }

  // The numbers expected should be obvious from the diagram above, keeping in
  // mind that -1 hops really means infinite
  int aiTotal[] =
    { 9, 2, 0, 4, 0, 0, 2, 0, 0, 0 };
  int aiDirect[] =
    { 3, 2, 0, 2, 0, 0, 2, 0, 0, 0 };
  int aiTwoLayer[] =
    { 7, 2, 0, 4, 0, 0, 2, 0, 0, 0 };
  int aiOffset[] =
    { 1, 4, 10, 6, 10, 10, 8, 10, 10, 10 };
  qGotAllSets = true;
  int iSet;
  for (iSet = 0; iSet < 10; iSet++) {
    int iTotal, iDirect, iTwoLayer;
    int iTotalRet, iDirectRet, iTwoLayerRet;

    iMesh_getNumEntSets(instance, apvEntSets[iSet], -1, &iTotal, &err);
    CHECK_ERR(err);
    iMesh_getEntSets(instance, apvEntSets[iSet], -1, &aSets, &iSetsAlloc,
		     &iTotalRet, &err);
    CHECK_ERR(err);
    GOT_ALL_ITEMS((apvEntSets + aiOffset[iSet]), aSets, iTotalRet, qGotAllSets);
    TEST(
	(aiTotal[iSet] == iTotal) && (aiTotal[iSet] == iTotalRet)
	    && qGotAllSets);

    iMesh_getNumEntSets(instance, apvEntSets[iSet], 0, &iDirect, &err);
    CHECK_ERR(err);
    iMesh_getEntSets(instance, apvEntSets[iSet], 0, &aSets, &iSetsAlloc,
		     &iDirectRet, &err);
    CHECK_ERR(err);
    GOT_ALL_ITEMS((apvEntSets + aiOffset[iSet]), aSets, iDirectRet,
		  qGotAllSets);
    TEST(
	(aiDirect[iSet] == iDirect) && (aiDirect[iSet] == iDirectRet)
	    && qGotAllSets);

    iMesh_getNumEntSets(instance, apvEntSets[iSet], 1, &iTwoLayer, &err);
    CHECK_ERR(err);
    iMesh_getEntSets(instance, apvEntSets[iSet], 1, &aSets, &iSetsAlloc,
		     &iTwoLayerRet, &err);
    CHECK_ERR(err);
    GOT_ALL_ITEMS((apvEntSets + aiOffset[iSet]), aSets, iTwoLayerRet,
		  qGotAllSets);
    TEST(
	(aiTwoLayer[iSet] == iTwoLayer) && (aiTwoLayer[iSet] == iTwoLayerRet)
	    && qGotAllSets);
  }

  // Test the ability to remove sets from sets
  qGotAllSets = true;
  for (i = 1; i <= 3; i++) {
    iMesh_rmvEntSet(instance, apvEntSets[i], apvEntSets[0], &err);
    CHECK_ERR(err);
    int iIsContained;
    iMesh_isEntSetContained(instance, apvEntSets[0], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && !iIsContained;
  }
  for (i = 4; i <= 5; i++) {
    iMesh_rmvEntSet(instance, apvEntSets[i], apvEntSets[1], &err);
    CHECK_ERR(err);
    int iIsContained;
    iMesh_isEntSetContained(instance, apvEntSets[1], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && !iIsContained;
  }
  for (i = 6; i <= 7; i++) {
    iMesh_rmvEntSet(instance, apvEntSets[i], apvEntSets[3], &err);
    CHECK_ERR(err);
    int iIsContained;
    iMesh_isEntSetContained(instance, apvEntSets[3], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && !iIsContained;
  }
  for (i = 8; i <= 9; i++) {
    iMesh_rmvEntSet(instance, apvEntSets[i], apvEntSets[4], &err);
    CHECK_ERR(err);
    int iIsContained;
    iMesh_isEntSetContained(instance, apvEntSets[4], apvEntSets[i],
			    &iIsContained, &err);
    CHECK_ERR(err);
    qGotAllSets = qGotAllSets && !iIsContained;
  }
  TEST(qGotAllSets);
  cout << "DONE\n";

  cout << "Testing addition and removal for entities in sets...      ";
  iBase_EntityHandle aEnts[10];

  // Get some entities to work with (any generic entities will do)
  GET_ENTS(aEnts, 10);

  // Create a hierarchy of sets where            [0]
  //     set 0 contains sets 1-3,           ______|_____
  //     set 1 contains sets 4-5,           [1], [2], [3]
  //     set 3 contains sets 6-7,         ___|___   ___|___
  // and set 4 contains sets 8-9          [4],[5]   [6],[7]
  //                                              ___|___
  // (and test the contain query)                 [8],[9]
  //
  // This has already been done above, so presumably it'll just work.
  //
  for (i = 1; i <= 3; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[0], &err);
    CHECK_ERR(err);
  }
  for (i = 4; i <= 5; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[1], &err);
    CHECK_ERR(err);
  }
  for (i = 6; i <= 7; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[3], &err);
    CHECK_ERR(err);
  }
  for (i = 8; i <= 9; i++) {
    iMesh_addEntSet(instance, apvEntSets[i], apvEntSets[6], &err);
    CHECK_ERR(err);
  }

  // Ensure that none of the entities are already in the newly created sets
  bool qEmpty = true;
  for (i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      int iCont;
      iMesh_isEntContained(instance, apvEntSets[i], aEnts[j], &iCont, &err);
      CHECK_ERR(err);
      qEmpty = qEmpty && !iCont;
    }
  }
  TEST(qEmpty);

  // Add all the entities to each of the sets in two different ways
  // Since the first five sets are lists, they should end up with
  // two copies of each entity
  for (i = 0; i < 10; i++) {
    iMesh_addEntArrToSet(instance, aEnts, 10, apvEntSets[i], &err);
    CHECK_ERR(err);
    for (int j = 0; j < 10; j++) {
      iMesh_addEntToSet(instance, aEnts[j], apvEntSets[i], &err);
      CHECK_ERR(err);
    }
  }

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsContain(instance, pvRootSet, apvEntSets);
    iBase_EntityHandle *aEntsTmp = NULL;
    int aEntsTmpAlloc = 0, aEntsTmpSize;
    iMesh_getEntities(instance, apvEntSets[1], iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, &aEntsTmp, &aEntsTmpAlloc,
		      &aEntsTmpSize, &err);
    CHECK_ERR(err);
    TESTEQ(aEntsTmpSize, 10);
    for (int j = 0; j < 10; j++)
      aEnts[j] = aEntsTmp[j];
    UT_FREE(aEntsTmp);
  }

  // Test that all of the entities are in each of the sets
  bool qAllIn = true;
  for (i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      int iCont;
      iMesh_isEntContained(instance, apvEntSets[i], aEnts[j], &iCont, &err);
      CHECK_ERR(err);
      qAllIn = qAllIn && iCont;
    }
  }
  TEST(qAllIn);

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsContain(instance, pvRootSet, apvEntSets);
    iBase_EntityHandle *aEntsTmp = NULL;
    int aEntsTmpAlloc = 0, aEntsTmpSize;
    iMesh_getEntities(instance, apvEntSets[1], iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, &aEntsTmp, &aEntsTmpAlloc,
		      &aEntsTmpSize, &err);
    CHECK_ERR(err);
    TESTEQ(aEntsTmpSize, 10);
    for (int j = 0; j < 10; j++)
      aEnts[j] = aEntsTmp[j];
    UT_FREE(aEntsTmp);
  }

  // Remove one copy of the entities using two different function calls
  for (i = 0; i < 5; i++) {
    iMesh_rmvEntArrFromSet(instance, aEnts, 10, apvEntSets[i], &err);
    CHECK_ERR(err);
    for (int j = 0; j < 10; j++) {
      iMesh_rmvEntFromSet(instance, aEnts[j], apvEntSets[i + 5], &err);
      CHECK_ERR(err);
    }
  }

  // Check that the lists still contain a copy of each entity
  qAllIn = true;
  for (i = 0; i < 5; i++) {
    for (int j = 0; j < 10; j++) {
      int iCont;
      iMesh_isEntContained(instance, apvEntSets[2 * i], aEnts[j], &iCont, &err);
      CHECK_ERR(err);
      qAllIn = qAllIn && iCont;
    }
  }
  TEST(qAllIn);

  // Check that the sets no longer have any of the entities
  qEmpty = true;
  for (i = 0; i < 5; i++) {
    for (int j = 0; j < 10; j++) {
      int iCont;
      iMesh_isEntContained(instance, apvEntSets[2 * i + 1], aEnts[j], &iCont,
			   &err);
      CHECK_ERR(err);
      qEmpty = qEmpty && !iCont;
    }
  }
  TEST(qEmpty);

  // Remove the entities from the lists...
  for (i = 0; i < 5; i++) {
    iMesh_rmvEntArrFromSet(instance, aEnts, 10, apvEntSets[2 * i], &err);
    CHECK_ERR(err);
  }

  // ...and ensure that the lists and sets are indeed empty
  qEmpty = true;
  for (i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      int iCont;
      iMesh_isEntContained(instance, apvEntSets[i], aEnts[j], &iCont, &err);
      CHECK_ERR(err);
      qEmpty = qEmpty && !iCont;
    }
  }
  TEST(qEmpty);
  cout << "DONE\n";

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsContain(instance, pvRootSet, apvEntSets);
  }

  if (qErrTests) {
    cout << "Testing error throwing...                                 ";
    // Error throwing tests simply try to create a situation where the interface
    // should, according to the spec, throw an error, and verify that:
    // (a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reversal of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice versa, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.

    // Adding set to root set
    iMesh_addEntSet(instance, apvEntSets[0], pvRootSet, &err);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    // Adding root set to another set
    iMesh_addEntSet(instance, pvRootSet, apvEntSets[0], &err);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    // Removing set from root set
    iMesh_rmvEntSet(instance, apvEntSets[0], pvRootSet, &err);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    // Removing root set from set
    iMesh_rmvEntSet(instance, pvRootSet, apvEntSets[0], &err);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    cout << "DONE\n";
  }
  for (i = 0; i < 10; i++) {
    iMesh_destroyEntSet(instance, apvEntSets[i], &err);
    CHECK_ERR(err);
  }
  vDoSaveLoad(&instance, &pvRootSet);
  UT_FREE(aSets);
}
