#include "iMesh_unitTest.hh"

void
vTestSetRelationInterface(iMesh_Instance& instance)
{
  int err;
  //     addPrntChld, rmvPrntChld, isChildOf, getNumChld, getNumPrnt,
  //     getChldn, getPrnts
  cout << "Testing relations and parent/child queries...             ";
  iBase_EntitySetHandle apvSets[10];
  iBase_EntitySetHandle *aSets = (iBase_EntitySetHandle*) calloc(
      10, sizeof(iBase_EntitySetHandle));
  int i, iNumSets, iSetsAlloc = 10;
  int iIsChild;

  // Create ten sets; the even ones are lists
  for (i = 0; i < 10; i++) {
    iMesh_createEntSet(instance, (i % 2 == 0), &apvSets[i], &err);
    CHECK_ERR(err);
  }

  // Create a hierarchy of sets where                     [0]
  //     set 0 is the parent of sets 1-3,           _______|_______
  //     set 1 is the parent of sets 4-5,           [1],  [2],  [3]
  //     set 3 is the parent of sets 6-7,         ___|___     ___|___
  // and set 6 is the parent of sets 8-9          [4],[5]     [6],[7]
  //                                                        ___|___
  //                                                        [8],[9]
  // Also test the isChildOf query
  iIsChild = true;
  bool qOK = true;
  for (i = 1; i <= 3; i++) {
    iMesh_addPrntChld(instance, apvSets[0], apvSets[i], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[0], apvSets[i], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && iIsChild;
  }

  for (int j = 4; j <= 5; j++) {
    iMesh_addPrntChld(instance, apvSets[1], apvSets[j], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[1], apvSets[j], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && iIsChild;
  }

  for (int j = 6; j <= 7; j++) {
    iMesh_addPrntChld(instance, apvSets[3], apvSets[j], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[3], apvSets[j], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && iIsChild;
  }

  for (int j = 8; j <= 9; j++) {
    iMesh_addPrntChld(instance, apvSets[6], apvSets[j], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[6], apvSets[j], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && iIsChild;
  }
  TEST(qOK);
  cout << "DONE\n";

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsParentChild(instance, pvRootSet, apvSets);
  }

  cout << "Testing ability to get children/parents...                ";
  int aiTotalChld[] =
    { 9, 2, 0, 4, 0, 0, 2, 0, 0, 0 };
  int aiDirectChld[] =
    { 3, 2, 0, 2, 0, 0, 2, 0, 0, 0 };
  int aiTwoLayerChld[] =
    { 7, 2, 0, 4, 0, 0, 2, 0, 0, 0 };
  int aiTotalPrnt[] =
    { 0, 1, 1, 1, 2, 2, 2, 2, 3, 3 };
  int aiDirectPrnt[] =
    { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  int aiTwoLayerPrnt[] =
    { 0, 1, 1, 1, 2, 2, 2, 2, 2, 2 };

  iBase_EntitySetHandle a2hTotalChld[][9] =
    {
      { apvSets[1], apvSets[2], apvSets[3], apvSets[4], apvSets[5], apvSets[6],
	  apvSets[7], apvSets[8], apvSets[9] },
      { apvSets[4], apvSets[5], NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { apvSets[6], apvSets[7], apvSets[8], apvSets[9], NULL, NULL, NULL,
      NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { apvSets[8], apvSets[9], NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL } };
  iBase_EntitySetHandle a2hDirectChld[][3] =
    {
      { apvSets[1], apvSets[2], apvSets[3] },
      { apvSets[4], apvSets[5], NULL },
      { NULL, NULL, NULL },
      { apvSets[6], apvSets[7], NULL },
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
      { apvSets[8], apvSets[9], NULL },
      { NULL, NULL, NULL },
      { NULL, NULL, NULL },
      { NULL, NULL, NULL } };
  iBase_EntitySetHandle a2hTwoLayerChld[][7] =
    {
      { apvSets[1], apvSets[2], apvSets[3], apvSets[4], apvSets[5], apvSets[6],
	  apvSets[7] },
      { apvSets[4], apvSets[5], NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { apvSets[6], apvSets[7], apvSets[8], apvSets[9], NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { apvSets[8], apvSets[9], NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL },
      { NULL, NULL, NULL, NULL, NULL, NULL, NULL } };

  iBase_EntitySetHandle a2hTotalPrnt[][3] =
    {
      { NULL, NULL, NULL },
      { apvSets[0], NULL, NULL },
      { apvSets[0], NULL, NULL },
      { apvSets[0], NULL, NULL },
      { apvSets[0], apvSets[1], NULL },
      { apvSets[0], apvSets[1], NULL },
      { apvSets[0], apvSets[3], NULL },
      { apvSets[0], apvSets[3], NULL },
      { apvSets[0], apvSets[3], apvSets[6] },
      { apvSets[0], apvSets[3], apvSets[6] } };
  iBase_EntitySetHandle a2hDirectPrnt[][1] =
    {
      { NULL },
      { apvSets[0] },
      { apvSets[0] },
      { apvSets[0] },
      { apvSets[1] },
      { apvSets[1] },
      { apvSets[3] },
      { apvSets[3] },
      { apvSets[6] },
      { apvSets[6] } };
  iBase_EntitySetHandle a2hTwoLayerPrnt[][2] =
    {
      { NULL, NULL },
      { apvSets[0], NULL },
      { apvSets[0], NULL },
      { apvSets[0], NULL },
      { apvSets[0], apvSets[1] },
      { apvSets[0], apvSets[1] },
      { apvSets[0], apvSets[3] },
      { apvSets[0], apvSets[3] },
      { apvSets[3], apvSets[6] },
      { apvSets[3], apvSets[6] } };

  bool qNumsOK = true;
  for (int i = 0; i < 10; i++) {
    int iTotalChld, iDirectChld, iTwoLayerChld;
    int iTotalPrnt, iDirectPrnt, iTwoLayerPrnt;

    iMesh_getNumChld(instance, apvSets[i], -1, &iTotalChld, &err);
    CHECK_ERR(err);
    iMesh_getNumChld(instance, apvSets[i], 0, &iDirectChld, &err);
    CHECK_ERR(err);
    iMesh_getNumChld(instance, apvSets[i], 1, &iTwoLayerChld, &err);
    CHECK_ERR(err);

    iMesh_getNumPrnt(instance, apvSets[i], -1, &iTotalPrnt, &err);
    CHECK_ERR(err);
    iMesh_getNumPrnt(instance, apvSets[i], 0, &iDirectPrnt, &err);
    CHECK_ERR(err);
    iMesh_getNumPrnt(instance, apvSets[i], 1, &iTwoLayerPrnt, &err);
    CHECK_ERR(err);

    qNumsOK = ((iTotalChld == aiTotalChld[i])
	&& (iDirectChld == aiDirectChld[i])
	&& (iTwoLayerChld == aiTwoLayerChld[i])
	&& (iTotalPrnt == aiTotalPrnt[i]) && (iDirectPrnt == aiDirectPrnt[i])
	&& (iTwoLayerPrnt == aiTwoLayerPrnt[i]));
    TEST(qNumsOK);

    // Check that the sets can return the parents and children
    // counted above
    ///
    // Total children
    iMesh_getChldn(instance, apvSets[i], -1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    bool qGotAllSets = false;
    GOT_ALL_ITEMS(a2hTotalChld[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiTotalChld[i]);

    // Direct children
    iMesh_getChldn(instance, apvSets[i], 0, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    qGotAllSets = false;
    GOT_ALL_ITEMS(a2hDirectChld[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiDirectChld[i]);

    // TwoLayer children
    iMesh_getChldn(instance, apvSets[i], 1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    qGotAllSets = false;
    GOT_ALL_ITEMS(a2hTwoLayerChld[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiTwoLayerChld[i]);

    // Total parents
    iMesh_getPrnts(instance, apvSets[i], -1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    qGotAllSets = false;
    GOT_ALL_ITEMS(a2hTotalPrnt[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiTotalPrnt[i]);

    // Direct parents
    iMesh_getPrnts(instance, apvSets[i], 0, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    qGotAllSets = false;
    GOT_ALL_ITEMS(a2hDirectPrnt[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiDirectPrnt[i]);

    // TwoLayer parents
    iMesh_getPrnts(instance, apvSets[i], 1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR(err);
    qGotAllSets = false;
    GOT_ALL_ITEMS(a2hTwoLayerPrnt[i], aSets, iNumSets, qGotAllSets);
    TEST(qGotAllSets && iNumSets == aiTwoLayerPrnt[i]);
  }
  cout << "DONE\n";

  double dData;
  void* pvData = &dData;
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

    // Root set as parent
    iMesh_addPrntChld(instance, pvRootSet, apvSets[0], &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Root set as child
    iMesh_addPrntChld(instance, apvSets[0], pvRootSet, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Removing root as parent
    iMesh_rmvPrntChld(instance, pvRootSet, apvSets[0], &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Removing root as child
    iMesh_rmvPrntChld(instance, apvSets[0], pvRootSet, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Checking for root as parent
    iMesh_isChildOf(instance, pvRootSet, apvSets[0], &iIsChild, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Checking for root as child
    iMesh_isChildOf(instance, apvSets[0], pvRootSet, &iIsChild, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Getting child count for root set
    iMesh_getNumChld(instance, pvRootSet, -1, &iNumSets, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Getting parent count for root set
    iMesh_getNumPrnt(instance, pvRootSet, -1, &iNumSets, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Getting children for root set
    iMesh_getChldn(instance, pvRootSet, -1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    // Getting parents for root set
    iMesh_getPrnts(instance, pvRootSet, -1, &aSets, &iSetsAlloc, &iNumSets,
		   &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITYSET_HANDLE);

    cout << "DONE\n";
  }

  cout << "Testing relation removal...                               ";

  // Ensure that the parent-child relations can be removed without throwing
  // an error an that the query reflects this removal
  qOK = true;
  for (i = 1; i <= 3; i++) {
    iMesh_rmvPrntChld(instance, apvSets[0], apvSets[i], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[0], apvSets[i], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && (!iIsChild);
  }
  for (i = 4; i <= 5; i++) {
    iMesh_rmvPrntChld(instance, apvSets[1], apvSets[i], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[1], apvSets[i], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && (!iIsChild);
  }
  for (i = 6; i <= 7; i++) {
    iMesh_rmvPrntChld(instance, apvSets[3], apvSets[i], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[3], apvSets[i], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && (!iIsChild);
  }
  for (i = 8; i <= 9; i++) {
    iMesh_rmvPrntChld(instance, apvSets[6], apvSets[i], &err);
    CHECK_ERR(err);
    iMesh_isChildOf(instance, apvSets[6], apvSets[i], &iIsChild, &err);
    CHECK_ERR(err);
    qOK = qOK && (!iIsChild);
  }

  // Clean up the sets
  for (i = 0; i < 10; i++) {
    iMesh_destroyEntSet(instance, apvSets[i], &err);
    CHECK_ERR(err);
  }
  UT_FREE(aSets);
  vDoSaveLoad(&instance, &pvRootSet);
  cout << "DONE\n";
}
