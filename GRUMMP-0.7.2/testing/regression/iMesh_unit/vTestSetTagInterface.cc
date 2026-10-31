#include "iMesh_unitTest.hh"

void
vTestSetTagInterface(iMesh_Instance& instance)
{
  iBase_TagHandle *aGotTags = NULL;
  iBase_TagHandle aTags[7];
  iBase_EntitySetHandle apvSets[6] =
    { NULL };
  iBase_EntityHandle aEnts[3];
  int iNumChars, err;
  // Testing normal operation of SetTag
  //     getEntSetData, getEntSetIntData, getEntSetDblData,
  //     getEntSetEHData, setEntSetData, setEntSetIntData,
  //     setEntSetDblData, setEntSetEHData, getAllEntSetTags,
  //     rmvEntSetTag
  cout << "Testing entity set tagging functionality...\n";

  {
    iBase_EntityIterator pvIter;
    // Get at least 3 entities to work with
    iMesh_initEntIter(instance, pvRootSet, iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, 0, &pvIter, &err);
    CHECK_ERR(err);
    if (err != iBase_SUCCESS) {
      cout << "\nCould not get mesh data to work with... exiting...\n";
      TEST(data_preparation);
      return;
    }

    int qGotData, i = 0;
    iMesh_getNextEntIter(instance, pvIter, &aEnts[i], &qGotData, &err);
    CHECK_ERR(err);
    i++;
    while (qGotData && i < 3) {
      iMesh_getNextEntIter(instance, pvIter, &aEnts[i], &qGotData, &err);
      CHECK_ERR(err);
      i++;
    }
    if (!qGotData) {
      cout << "\nCould not get enough mesh data to work with... exiting...\n";
      TEST(data_preparation);
      return;
    }
    iMesh_endEntIter(instance, pvIter, &err);
    CHECK_ERR(err);
  }

  // Create a tag of each type and 3 entity sets to be tagged
  for (int i = 0; i < 3; i++) {
    iMesh_createTag(instance, astrTagNames[i], 1, aTagTypes[i], &aTags[i], &err,
		    aiTagNameLen[i]);
    CHECK_ERR(err);
    iMesh_createEntSet(instance, false, &apvSets[i], &err);
    CHECK_ERR(err);
    if (err != iBase_SUCCESS) {
      cout << "Could not create entity sets to work with... " << "exiting...\n";
      TEST(data_preparation);
      return;
    }
  }

  // Be sure to test tagging of the root set (mesh) as well.  Yes, this
  // leaks a set; so what?
  apvSets[0] = pvRootSet;
  // ... and a generic tag for completeness
  iMesh_createTag(instance, astrTagNames[3], 6, aTagTypes[3], &aTags[3], &err,
		  aiTagNameLen[3]);
  CHECK_ERR(err);

  cout << "Testing get and set methods with integers...              ";
  {
    // Try setting the data using the generic function...
    iNumChars = sizeof(int);
    int* piData = new int(25);
    char *aChars = (char*) piData;

    iMesh_setEntSetData(instance, apvSets[0], aTags[0], aChars, iNumChars,
			&err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getEntSetData(instance, apvSets[0], aTags[0], &aOutput, &iOutputAlloc,
			&iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(iNumChars, (int )sizeof(int));
    TESTEQ(*(int* )aOutput, 25);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setEntSetIntData(instance, apvSets[0], aTags[0], 30, &err);
    CHECK_ERR(err);
    int iVal;
    iMesh_getEntSetIntData(instance, apvSets[0], aTags[0], &iVal, &err);
    CHECK_ERR(err);
    TESTEQ(iVal, 30);
    UT_FREE(aOutput);
    delete piData;
  }
  cout << "DONE\n";

  cout << "Testing get and set methods with doubles...               ";
  {
    // Try setting the data using the generic function...
    iNumChars = sizeof(double);
    double* pdData = new double(1.8);
    char *aChars = (char*) pdData;

    iMesh_setEntSetData(instance, apvSets[0], aTags[1], aChars, iNumChars,
			&err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getEntSetData(instance, apvSets[0], aTags[1], &aOutput, &iOutputAlloc,
			&iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(iNumChars, (int )sizeof(double));
    TESTEQ(*(double* )aOutput, 1.8);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setEntSetDblData(instance, apvSets[0], aTags[1], 1.9, &err);
    CHECK_ERR(err);
    double dVal;
    iMesh_getEntSetDblData(instance, apvSets[0], aTags[1], &dVal, &err);
    CHECK_ERR(err);
    TESTEQ(dVal, 1.9);
    UT_FREE(aOutput);
    delete pdData;
  }
  cout << "DONE\n";

  cout << "Testing get and set methods with entity handles...        ";
  {
    // Try setting the data using the generic function...
    iNumChars = sizeof(void*);
    char *aChars = (char*) (&aEnts[0]);

    iMesh_setEntSetData(instance, apvSets[0], aTags[2], aChars, iNumChars,
			&err);
    CHECK_ERR(err);

    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getEntSetData(instance, apvSets[0], aTags[2], &aOutput, &iOutputAlloc,
			&iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(iNumChars, (int )sizeof(void*));
    TESTEQ(*(iBase_EntityHandle* )aOutput, aEnts[0]);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setEntSetEHData(instance, apvSets[0], aTags[2], aEnts[1], &err);
    CHECK_ERR(err);

    iBase_EntityHandle EH;
    iMesh_getEntSetEHData(instance, apvSets[0], aTags[2], &EH, &err);
    CHECK_ERR(err);
    UT_FREE(aOutput);
    TESTEQ(EH, aEnts[1]);
  }
  cout << "DONE\n";

  cout << "Testing get and set methods with generic data (bytes)...  ";
  {
    // Try setting the data using the generic function...
    iNumChars = 6;
    char acData[] = "abcdef";
    iMesh_setEntSetData(instance, apvSets[0], aTags[3], acData, iNumChars,
			&err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getEntSetData(instance, apvSets[0], aTags[3], &aOutput, &iOutputAlloc,
			&iNumChars, &err);
    CHECK_ERR(err);
    bool qGotData = (iNumChars == 6)
	&& !strncmp(reinterpret_cast<char*>(aOutput), acData, iNumChars);
    UT_FREE(aOutput);
    TEST(qGotData);
  }
  cout << "DONE\n";

  cout << "Testing tag finding...                                    ";
  {
    int iNumTags, iGotTagsAlloc = 0;
    iMesh_getAllEntSetTags(instance, apvSets[0], &aGotTags, &iGotTagsAlloc,
			   &iNumTags, &err);
    CHECK_ERR(err);
    bool qTagsPresent = false;
    if (iNumTags == 4)
      GOT_ALL_ITEMS(aTags, aGotTags, iNumTags, qTagsPresent);
    TEST(qTagsPresent);
  }
  cout << "DONE\n";
  vDoSaveLoad(&instance, &pvRootSet);

  if (qSaveLoad) {
    iBase_EntityHandle ent;
    vRecoverTags(instance, pvRootSet, aTags);
    apvSets[0] = pvRootSet;

    cout << "Retesting get methods after load/save...\n";
    cout << "Testing get method with integers...              ";
    {
      int iVal;
      iMesh_getEntSetIntData(instance, apvSets[0], aTags[0], &iVal, &err);
      CHECK_ERR(err);
      TESTEQ(iVal, 30);
    }
    cout << "DONE\n";

    cout << "Testing get method with doubles...               ";
    {
      double dVal;
      iMesh_getEntSetDblData(instance, apvSets[0], aTags[1], &dVal, &err);
      CHECK_ERR(err);
      TESTEQ(dVal, 1.9);
    }
    cout << "DONE\n";

    // No practical way to test entity handle tags:  what's the right
    // answer?

    cout << "Testing get method with generic data (bytes)...  ";
    {
      int iNumChars;
      char acData[] = "abcdef";
      void *aOutput = NULL;
      int iOutputAlloc = 0;
      iMesh_getEntSetData(instance, apvSets[0], aTags[3], &aOutput,
			  &iOutputAlloc, &iNumChars, &err);
      CHECK_ERR(err);
      bool qGotData = (iNumChars == 6)
	  && !strncmp(reinterpret_cast<char*>(aOutput), acData, iNumChars);
      UT_FREE(aOutput);
      TEST(qGotData);
    }
    cout << "DONE\n";
  }

  // Testing error behavior
  if (qErrTests) {
    // Error throwing tests simply try to create a situation where the interface
    // should, according to the spec, throw an error, and verify that:
    // (a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reversal of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice versa, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.
    {
      iMesh_Instance new_instance;
      iMesh_newMesh("", &new_instance, &err, 0);
      cout << "Testing error throwing for SetTag methods...              ";
      char data[] =
	{ 'a', 'b', 'c', 'd' };
      void *aOutput = data;
      int iOutputAlloc = 4;
      {
	for (int i = 0; i < 3; i++) {
	  iMesh_createTag(new_instance, astrTagNames[i], 1, aTagTypes[i],
			  &aTags[i + 4], &err, aiTagNameLen[i]);
	  CHECK_ERR(err);
	  iMesh_createEntSet(new_instance, false, &apvSets[i + 3], &err);
	  CHECK_ERR(err);
	  if (err != iBase_SUCCESS) {
	    cout << "Could not create entity sets to work with... "
		<< "exiting...\n";
	    TEST(data_preparation);
	    return;
	  }
	}
      }
      //setting tag data for a mesh with no data
      iMesh_setEntSetData(new_instance, apvSets[3], aTags[4], aOutput,
			  iOutputAlloc, &err);
      CHECK_ERR(err);

      // setting tag data for a mesh with no data
      iMesh_setEntSetIntData(new_instance, apvSets[3], aTags[4], 15, &err);
      CHECK_ERR(err);

      //  setting tag data for a mesh with no data
      iMesh_setEntSetDblData(new_instance, apvSets[3], aTags[5], 1.5, &err);
      CHECK_ERR(err);

      // setting tag data for a mesh with no data
      iMesh_setEntSetEHData(new_instance, apvSets[3], aTags[6], aEnts[0], &err);
      CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

      //  getting tag data for a mesh with no data
      iMesh_getEntSetData(new_instance, apvSets[3], aTags[4], &aOutput,
			  &iOutputAlloc, &iNumChars, &err);
      CHECK_ERR(err);

      //  getting tag data for a mesh with no data
      int iData;
      iMesh_getEntSetIntData(new_instance, apvSets[3], aTags[4], &iData, &err);
      CHECK_ERR(err);

      //  getting tag data for a mesh with no data
      double dData;
      iMesh_getEntSetDblData(new_instance, apvSets[3], aTags[5], &dData, &err);
      CHECK_ERR(err);

      //  getting tag data for a mesh with no data
      iBase_EntityHandle EH;
      iMesh_getEntSetEHData(new_instance, apvSets[3], aTags[6], &EH, &err);
      CHECK_ERR2(err, iBase_TAG_NOT_FOUND);

      //  removing tag data for a mesh with no data
      iMesh_rmvEntSetTag(new_instance, apvSets[3], aTags[6], &err);
      CHECK_ERR(err);

      iMesh_dtor(new_instance, &err);
    }

    //  getting tag data for an invalid output array size
    int iOutputAlloc = 1;
    char aOutput[] = "";
    iMesh_getEntSetData(instance, apvSets[0], aTags[0], &aOutput, &iOutputAlloc,
			&iNumChars, &err);
    CHECK_ERR2(err, iBase_BAD_ARRAY_SIZE);

    // setting integer tag data for a double tag handle
    iMesh_setEntSetIntData(instance, apvSets[0], aTags[1], 15, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    //  setting double tag data for an integer tag handle
    iMesh_setEntSetDblData(instance, apvSets[0], aTags[0], 1.5, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    if (!qSaveLoad) {
      // If you've done save/load, then you'd need to re-establish the
      //  identify of the entity in question.  Yuck.

      //  setting entity handle tag data for an integer tag handle
      iMesh_setEntSetEHData(instance, apvSets[0], aTags[0], aEnts[0], &err);
      CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);
    }

    //  getting integer tag data for a double tag handle
    int iData;
    iMesh_getEntSetIntData(instance, apvSets[0], aTags[1], &iData, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    //  getting double tag data for an integer tag handle
    double dData;
    iMesh_getEntSetDblData(instance, apvSets[0], aTags[0], &dData, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    //  getting entity handle tag data for an integer tag handle
    iBase_EntityHandle EH;
    iMesh_getEntSetEHData(instance, apvSets[0], aTags[0], &EH, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    //  getting all tags into an array of the wrong size
    int iGotTagsAlloc = 1;
    iMesh_getAllEntSetTags(instance, apvSets[0], &aGotTags, &iGotTagsAlloc,
			   &iNumChars, &err);
    CHECK_ERR2(err, iBase_BAD_ARRAY_SIZE);

    cout << "DONE\n";
  }
  cout << "Testing tag removal and forced/unforced destruction...    ";
  {
    bool qOK = true;
    iMesh_rmvEntSetTag(instance, apvSets[0], aTags[0], &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_rmvEntSetTag(instance, apvSets[0], aTags[1], &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[0], false, &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[1], false, &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    TEST(qOK);
  }
  // trying to destroy tag with data without using forced
  iMesh_destroyTag(instance, aTags[2], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  //  trying to destroy tag with data without using forced
  iMesh_destroyTag(instance, aTags[3], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  {
    bool qOK = true;
    iMesh_destroyTag(instance, aTags[2], true, &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[3], true, &err);
    qOK = (err == iBase_SUCCESS);
    CHECK_ERR(err);

    TEST(qOK);
  }
  cout << "DONE\n";
  UT_FREE(aGotTags);
}
