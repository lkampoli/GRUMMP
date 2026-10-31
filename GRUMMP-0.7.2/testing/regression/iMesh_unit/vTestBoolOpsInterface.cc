#include "iMesh_unitTest.hh"

void
vTestBoolOpsInterface(iMesh_Instance& instance)
{
  int err;
  cout << "Testing methods with regular sets and lists...\n";
  iBase_EntitySetHandle aTestSets[12];
  iBase_EntityHandle aTestEnts[10];
  int i;

  // Get some entities to work with
  GET_ENTS(aTestEnts, 10);

  // Create a number of entity sets for testing
  // Call the entities in use a,b,c,d,e,f,g,h,i,j (indices 0-9)
  // Entity set 0 is a list and contains      : a,b,c,d,e,a,b,c,d,e,E1,E2,E3,E4,E5
  // Entity set 1 is a list and contains      : g,f,e,d,E2,E3,E4,E5
  // Entity set 2 is a list and contains      : f,f,g,g,h,h,i,i,j,j
  // Entity set 3 is not a list and contains  : a,b,c,d,e,E1
  // Entity set 4 is not a list and contains  : d,e,f,g,E1,E2
  // Entity set 5 is not a list and contains  : f,g,h,i,j
  for (i = 0; i < 6; i++) {
    iMesh_createEntSet(instance, i < 3, &aTestSets[i], &err);
    CHECK_ERR(err);
  }

  for (i = 0; i < 5; i++) {
    /* Add ents a,b,c,d,e to sets 0 and 3. */
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[0], &err);
    CHECK_ERR(err);
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[3], &err);
    CHECK_ERR(err);
    /* Add sets 1, 2, 3, 4, 5 to set 0 */
    iMesh_addEntSet(instance, aTestSets[i + 1], aTestSets[0], &err);
    CHECK_ERR(err);
  }

  /* Add ents a,b,c,d,e to sets 0 again. */
  for (i = 0; i < 5; i++) {
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[0], &err);
    CHECK_ERR(err);
  }

  for (i = 3; i < 7; i++) {
    // Add g,f,e,d to set 1
    iMesh_addEntToSet(instance, aTestEnts[9 - i], aTestSets[1], &err);
    CHECK_ERR(err);

    // Add sets 2,3,4,5 to set 1
    iMesh_addEntSet(instance, aTestSets[i - 1], aTestSets[1], &err);
    CHECK_ERR(err);

    // Add ents d, e, f, g to set 4
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[4], &err);
    CHECK_ERR(err);

    // Add sets 1 and 2 twice each to set 4.
    iMesh_addEntSet(instance, aTestSets[(i + 1) / 3], aTestSets[4], &err);
    CHECK_ERR(err);
  }

  for (i = 5; i < 10; i++) {
    // Add f,g,h,i,j twice in succession to both set 2 and set 5
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[2], &err);
    CHECK_ERR(err);
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[2], &err);
    CHECK_ERR(err);
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[5], &err);
    CHECK_ERR(err);
    iMesh_addEntToSet(instance, aTestEnts[i], aTestSets[5], &err);
    CHECK_ERR(err);
  }

  // Add set 1 to set 3.
  iMesh_addEntSet(instance, aTestSets[1], aTestSets[3], &err);
  CHECK_ERR(err);

  cout << "Testing intersect with regular lists...                   ";

  // E0 n E2; expect empty list
  iMesh_intersect(instance, aTestSets[0], aTestSets[2], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], true, aTestSets, 0, aTestEnts, 0);

  // E0 n E1; expect d,e,E2,E3,E4,E5
  iMesh_intersect(instance, aTestSets[0], aTestSets[1], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], true, aTestSets + 2, 4,
		     aTestEnts + 3, 2);

  // E1 n E0; expect e,d,E2,E3,E4,E5
  iMesh_intersect(instance, aTestSets[1], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[4], aTestEnts[3] };
    vCheckBoolOpResult(instance, aTestSets[8], true, aTestSets + 2, 4, aEnts,
		       2);
  }

  // E1 n E2; expect g,f
  iMesh_intersect(instance, aTestSets[1], aTestSets[2], &aTestSets[9], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[6], aTestEnts[5] };
    vCheckBoolOpResult(instance, aTestSets[9], true, aTestSets, 0, aEnts, 2);
  }

  // E2 n E0; expect empty list
  iMesh_intersect(instance, aTestSets[2], aTestSets[0], &aTestSets[10], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[10], true, aTestSets, 0, aTestEnts, 0);

  // E2 n E1; expect f,g
  iMesh_intersect(instance, aTestSets[2], aTestSets[1], &aTestSets[11], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[11], true, aTestSets, 0, aTestEnts + 5,
		     2);

  for (i = 6; i < 12; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing intersect with regular sets...                    ";

  // E3 n E4; expect d,e,E1 (unordered)
  iMesh_intersect(instance, aTestSets[3], aTestSets[4], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 1,
		     aTestEnts + 3, 2);

  // E4 n E3; expect d,e,E1 (unordered)
  iMesh_intersect(instance, aTestSets[4], aTestSets[3], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 1, 1,
		     aTestEnts + 3, 2);

  // E3 n E5; expect empty set
  iMesh_intersect(instance, aTestSets[3], aTestSets[5], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets, 0, aTestEnts, 0);

  // E5 n E3; expect empty set
  iMesh_intersect(instance, aTestSets[3], aTestSets[5], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets, 0, aTestEnts, 0);

  // E4 n E5; expect f,g (unordered)
  iMesh_intersect(instance, aTestSets[4], aTestSets[5], &aTestSets[10], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[10], false, aTestSets, 0,
		     aTestEnts + 5, 2);

  // E5 n E4; expect f,g (unordered)
  iMesh_intersect(instance, aTestSets[4], aTestSets[5], &aTestSets[11], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[11], false, aTestSets, 0,
		     aTestEnts + 5, 2);

  for (i = 6; i < 12; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing intersect with sets and lists...                  ";

  // E0 n E5; expect empty set
  iMesh_intersect(instance, aTestSets[0], aTestSets[5], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets, 0, aTestEnts, 0);

  // E0 n E4; expect d,e,E1,E2 (unordered)
  iMesh_intersect(instance, aTestSets[0], aTestSets[4], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 1, 2,
		     aTestEnts + 3, 2);

  // E4 n E0; expect d,e,E1,E2 (unordered)
  iMesh_intersect(instance, aTestSets[4], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets + 1, 2,
		     aTestEnts + 3, 2);

  // E4 n E2; expect f,g (unordered)
  iMesh_intersect(instance, aTestSets[4], aTestSets[2], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets, 0, aTestEnts + 5,
		     2);

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing subtract with regular lists...                    ";

  // E0 - E2; expect E0 = a,b,c,d,e,a,b,c,d,e,E1,E2,E3,E4,E5
  iMesh_subtract(instance, aTestSets[0], aTestSets[2], &aTestSets[6], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[10];
    for (int i_ = 0; i_ < 5; i_++) {
      aEnts[i_] = aEnts[i_ + 5] = aTestEnts[i_];
    }
    vCheckBoolOpResult(instance, aTestSets[6], true, aTestSets + 1, 5, aEnts,
		       10);
  }

  // E0 - E1; expect a,b,c,d,e,a,b,c,E1
  iMesh_subtract(instance, aTestSets[0], aTestSets[1], &aTestSets[7], &err);
  CHECK_ERR(err);
  {
    // This time, only eight of these ents should match...
    iBase_EntityHandle aEnts[10];
    for (int i_ = 0; i_ < 5; i_++) {
      aEnts[i_] = aEnts[i_ + 5] = aTestEnts[i_];
    }
    vCheckBoolOpResult(instance, aTestSets[7], true, aTestSets + 1, 1, aEnts,
		       8);
  }

  // E1 - E0; expect g,f
  iMesh_subtract(instance, aTestSets[1], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[6], aTestEnts[5] };
    vCheckBoolOpResult(instance, aTestSets[8], true, aTestSets, 0, aEnts, 2);
  }

  // E2 - E1; expect f,g,h,h,i,i,j,j
  iMesh_subtract(instance, aTestSets[2], aTestSets[1], &aTestSets[9], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[5], aTestEnts[6], aTestEnts[7], aTestEnts[7], aTestEnts[8],
	  aTestEnts[8], aTestEnts[9], aTestEnts[9] };
    vCheckBoolOpResult(instance, aTestSets[9], true, aTestSets, 0, aEnts, 8);
  }
  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing subtract with regular sets...                     ";

  // E3 - E5; expect E3 = a,b,c,d,e,E1 (unordered)
  iMesh_subtract(instance, aTestSets[3], aTestSets[5], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 1, aTestEnts,
		     5);

  // E3 - E4; expect a,b,c (unordered)
  iMesh_subtract(instance, aTestSets[3], aTestSets[4], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets, 0, aTestEnts, 3);

  // E4 - E3; expect f,g,E2 (unordered)
  iMesh_subtract(instance, aTestSets[4], aTestSets[3], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets + 2, 1,
		     aTestEnts + 5, 2);

  // E4 - E5; expect d,e,E1,E2 (unordered)
  iMesh_subtract(instance, aTestSets[4], aTestSets[5], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets + 1, 2,
		     aTestEnts + 3, 2);

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing subtract with lists and sets...                   ";

  // E0 - E5; expect a,b,c,d,e,E1,E2,E3,E4,E5 (unordered)
  iMesh_subtract(instance, aTestSets[0], aTestSets[5], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 5, aTestEnts,
		     5);

  // E0 - E4; expect a,b,c,E3,E4,E5 (unordered)
  iMesh_subtract(instance, aTestSets[0], aTestSets[4], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 3, 3, aTestEnts,
		     3);

  // E4 - E0; expect f,g (unordered)
  iMesh_subtract(instance, aTestSets[4], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets, 0, aTestEnts + 5,
		     2);

  // E4 - E2; expect d,e,E1,E2 (unordered)
  iMesh_subtract(instance, aTestSets[4], aTestSets[2], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets + 1, 2,
		     aTestEnts + 3, 2);

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing unite with regular lists...                       ";

  // E0 U E2; expect a,b,c,d,e,a,b,c,d,e,f,f,g,g,h,h,i,i,j,j,E1,E2,E3,E4,E5
  iMesh_unite(instance, aTestSets[0], aTestSets[2], &aTestSets[6], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[0], aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4],
	  aTestEnts[0], aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4],
	  aTestEnts[5], aTestEnts[5], aTestEnts[6], aTestEnts[6], aTestEnts[7],
	  aTestEnts[7], aTestEnts[8], aTestEnts[8], aTestEnts[9], aTestEnts[9] };
    vCheckBoolOpResult(instance, aTestSets[6], true, aTestSets + 1, 5, aEnts,
		       20);
  }

  // E0 U E1; expect a,b,c,d,e,a,b,c,d,e,g,f,e,d,E1,E2,E3,E4,E5
  iMesh_unite(instance, aTestSets[0], aTestSets[1], &aTestSets[7], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[0], aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4],
	  aTestEnts[0], aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4],
	  aTestEnts[6], aTestEnts[5], aTestEnts[4], aTestEnts[3] };
    vCheckBoolOpResult(instance, aTestSets[7], true, NULL, -1, aEnts, 14);
  }

  // E1 U E0; expect g,f,e,d,a,b,c,d,e,a,b,c,d,e,E1,E2,E3,E4,E5
  iMesh_unite(instance, aTestSets[1], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[6], aTestEnts[5], aTestEnts[4], aTestEnts[3], aTestEnts[0],
	  aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4], aTestEnts[0],
	  aTestEnts[1], aTestEnts[2], aTestEnts[3], aTestEnts[4] };
    vCheckBoolOpResult(instance, aTestSets[8], true, NULL, -1, aEnts, 14);
  }

  // E1 U E2; expect g,f,e,d,f,f,g,g,h,h,i,i,j,j,E2,E3,E4,E5
  iMesh_unite(instance, aTestSets[1], aTestSets[2], &aTestSets[9], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle aEnts[] =
      { aTestEnts[6], aTestEnts[5], aTestEnts[4], aTestEnts[3], aTestEnts[5],
	  aTestEnts[5], aTestEnts[6], aTestEnts[6], aTestEnts[7], aTestEnts[7],
	  aTestEnts[8], aTestEnts[8], aTestEnts[9], aTestEnts[9] };
    vCheckBoolOpResult(instance, aTestSets[9], true, aTestSets + 2, 4, aEnts,
		       14);
  }

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing unite with regular sets...                        ";

  // E3 U E5; expect a,b,c,d,e,f,g,h,i,j,E1 (unordered)
  iMesh_unite(instance, aTestSets[3], aTestSets[5], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 1, aTestEnts,
		     10);

  // E3 U E4; expect a,b,c,d,e,f,g,E1,E2 (unordered)
  iMesh_unite(instance, aTestSets[3], aTestSets[4], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 1, 2, aTestEnts,
		     7);

  // E4 U E3; expect a,b,c,d,e,f,g,E1,E2 (unordered)
  iMesh_unite(instance, aTestSets[4], aTestSets[3], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets + 1, 2, aTestEnts,
		     7);

  // E4 U E5; expect d,e,f,g,h,i,j,E1,E2 (unordered)
  iMesh_unite(instance, aTestSets[4], aTestSets[5], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets + 1, 2,
		     aTestEnts + 3, 7);

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing unite with lists and sets...                      ";

  // E0 U E5; expect a,b,c,d,e,f,g,h,i,j,E1,E2,E3,E4,E5 (unordered)
  iMesh_unite(instance, aTestSets[0], aTestSets[5], &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 5, aTestEnts,
		     10);

  // E0 U E4; expect a,b,c,d,e,f,g,E1,E2,E3,E4,E5 (unordered)
  iMesh_unite(instance, aTestSets[0], aTestSets[4], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 1, 5, aTestEnts,
		     7);

  // E4 U E0; expect a,b,c,d,e,f,g,E1,E2,E3,E4,E5 (unordered)
  iMesh_unite(instance, aTestSets[4], aTestSets[0], &aTestSets[8], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[8], false, aTestSets + 1, 5, aTestEnts,
		     7);

  // E4 U E2; expect d,e,f,g,h,i,j,E1,E2 (unordered)
  iMesh_unite(instance, aTestSets[4], aTestSets[2], &aTestSets[9], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[9], false, aTestSets + 1, 2,
		     aTestEnts + 3, 7);

  for (i = 6; i < 10; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing methods with the root set... may take a while...\n";

  cout << "Testing intersect with the root set...                    ";

  // E0 n R; expect a,b,c,d,e,E1,E2,E3,E4,E5 (unordered)
  iMesh_intersect(instance, aTestSets[0], pvRootSet, &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets + 1, 5, aTestEnts,
		     5);

  // R n E0; expect a,b,c,d,e,E1,E2,E3,E4,E5 (unordered)
  iMesh_intersect(instance, pvRootSet, aTestSets[0], &aTestSets[7], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[7], false, aTestSets + 1, 5, aTestEnts,
		     5);

  for (i = 6; i < 8; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";

  cout << "Testing subtract with the root set...                     ";

  // E0 - R; expect empty set
  iMesh_subtract(instance, aTestSets[0], pvRootSet, &aTestSets[6], &err);
  CHECK_ERR(err);
  vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets, 0, aTestEnts, 0);

  // R - E0; expect the complement of E0 (unordered; includes E0 and E6)
  iMesh_subtract(instance, pvRootSet, aTestSets[0], &aTestSets[7], &err);
  CHECK_ERR(err);

  {
    int iIsList;
    iMesh_isList(instance, aTestSets[7], &iIsList, &err);
    CHECK_ERR(err);

    iBase_EntitySetHandle *aGotSets = NULL;
    int iNumSets;
    int iGotSetsAlloc = 0;
    iMesh_getEntSets(instance, aTestSets[7], 0, &aGotSets, &iGotSetsAlloc,
		     &iNumSets, &err);
    CHECK_ERR(err);
    bool qGotAllItems = false;
    if (iNumSets == 2) {
      qGotAllItems = ((aGotSets[0] == aTestSets[0]
	  || aGotSets[1] == aTestSets[0])
	  && (aGotSets[0] == aTestSets[6] || aGotSets[1] == aTestSets[6]));
    }
    UT_FREE(aGotSets);

    int iNumEnts;
    iBase_EntityHandle * aGotEnts = NULL;
    int iGotEntsAlloc = 0;
    iMesh_getEntities(instance, aTestSets[7], iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, &aGotEnts, &iGotEntsAlloc,
		      &iNumEnts, &err);
    CHECK_ERR(err);
    int iNumAllEnts;
    iMesh_getNumOfType(instance, pvRootSet, iBase_ALL_TYPES, &iNumAllEnts,
		       &err);
    CHECK_ERR(err);

    bool qIsEmpty = false;
    if (iNumEnts == iNumAllEnts - 5) {
      // Size is right; check each entry in the result to make sure it's
      // not in the original set.
      qIsEmpty = true;
      for (i = 0; i < iNumEnts; i++) {
	int iIsCont;
	iMesh_isEntContained(instance, aTestSets[0], aGotEnts[i], &iIsCont,
			     &err);
	CHECK_ERR(err);
	qIsEmpty = (qIsEmpty && !iIsCont);
      }
    }
    TEST(qIsEmpty && qGotAllItems && !iIsList);

    for (i = 6; i < 8; i++) {
      iMesh_destroyEntSet(instance, aTestSets[i], &err);
      CHECK_ERR(err);
    }
    UT_FREE(aGotEnts);
  }
  cout << "DONE\n";

  cout << "Testing unite with the root set...                        ";
  // E0 U R; expect a copy of R (unordered)
  iMesh_unite(instance, aTestSets[0], pvRootSet, &aTestSets[6], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle *aAllEnts = NULL;
    int iAllEntsAlloc = 0, iNumAllEnts;
    iMesh_getEntities(instance, pvRootSet, iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, &aAllEnts, &iAllEntsAlloc,
		      &iNumAllEnts, &err);
    CHECK_ERR(err);
    vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets, 6, aAllEnts,
		       iNumAllEnts);
    UT_FREE(aAllEnts);
  }

  // R U E0; expect a copy of R (unordered)
  iMesh_unite(instance, pvRootSet, aTestSets[0], &aTestSets[7], &err);
  CHECK_ERR(err);
  {
    iBase_EntityHandle *aAllEnts = NULL;
    int iAllEntsAlloc = 0, iNumAllEnts;
    iMesh_getEntities(instance, pvRootSet, iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, &aAllEnts, &iAllEntsAlloc,
		      &iNumAllEnts, &err);
    CHECK_ERR(err);
    vCheckBoolOpResult(instance, aTestSets[6], false, aTestSets, 6, aAllEnts,
		       iNumAllEnts);
    UT_FREE(aAllEnts);
  }
  cout << "DONE\n";

  // Clean up the sets
  for (i = 0; i < 6; i++) {
    iMesh_destroyEntSet(instance, aTestSets[i], &err);
    CHECK_ERR(err);
  }
}
