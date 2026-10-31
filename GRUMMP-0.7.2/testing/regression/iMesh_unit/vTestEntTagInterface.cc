#include "iMesh_unitTest.hh"

void
vTestEntTagInterface(iMesh_Instance& instance)
{
  iBase_TagHandle aTags[7];
  iBase_EntityHandle aEnts[3];
  iBase_TagHandle * aGotTags = NULL;
  int iNumChars, err;

  //     getData, getIntData, getDblData, getEHData, setData,
  //     setIntData, setDblData, setEHData, getAllTags, rmvTag
  cout << "Testing basic tagging functionality...\n";

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

  int iNumInitTags;
  iBase_TagHandle *aInitTags = NULL;
  {
    int iGotTagsAlloc = 0;
    iMesh_getAllTags(instance, aEnts[0], &aInitTags, &iGotTagsAlloc,
		     &iNumInitTags, &err);
    CHECK_ERR(err);
  }

  // Create tags of each type that store a single value
  for (int i = 0; i < 3; i++) {
    iMesh_createTag(instance, astrTagNames[i], 1, aTagTypes[i], &aTags[i], &err,
		    aiTagNameLen[i]);
    CHECK_ERR(err);
  }
  // ... and a generic tag for completeness
  iMesh_createTag(instance, astrTagNames[3], 6, aTagTypes[3], &aTags[3], &err,
		  aiTagNameLen[3]);
  CHECK_ERR(err);

  cout << "Testing get and set methods with integers...              ";
  // Try setting the data using the integer function...
  {
    iNumChars = sizeof(int);
    char aChars[iNumChars];

    {
      int* piData = new int(25);
      char *pcData = (char*) piData;

      for (int i = 0; i < iNumChars; i++)
	aChars[i] = pcData[i];
      delete piData;
    }
    iMesh_setData(instance, aEnts[0], aTags[0], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact

    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getData(instance, aEnts[0], aTags[0], &aOutput, &iOutputAlloc,
		  &iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(iNumChars, (int )sizeof(int));
    TESTEQ(*(int* )aOutput, 25);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setIntData(instance, aEnts[0], aTags[0], 30, &err);
    CHECK_ERR(err);
    int iData;
    iMesh_getIntData(instance, aEnts[0], aTags[0], &iData, &err);
    CHECK_ERR(err);
    TESTEQ(iData, 30);
    UT_FREE(aOutput);
  } // Integer tests
  cout << "DONE\n";

  cout << "Testing get and set methods with doubles...               ";
  {
    // Try setting the data using the generic function...
    iNumChars = sizeof(double);
    char aChars[iNumChars];
    {
      double* pdData = new double(1.8);
      char* pcData = (char*) pdData;
      for (int i = 0; i < iNumChars; i++)
	aChars[i] = pcData[i];
      delete pdData;
    }
    iMesh_setData(instance, aEnts[0], aTags[1], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getData(instance, aEnts[0], aTags[1], &aOutput, &iOutputAlloc,
		  &iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(iNumChars, (int )sizeof(double));
    TESTEQ(*(double* )aOutput, 1.8);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setDblData(instance, aEnts[0], aTags[1], 1.9, &err);
    CHECK_ERR(err);
    double dData;
    iMesh_getDblData(instance, aEnts[0], aTags[1], &dData, &err);
    CHECK_ERR(err);
    UT_FREE(aOutput);
    TESTEQ(dData, 1.9);
  } // Testing doubles
  cout << "DONE\n";

  cout << "Testing get and set methods with entity handles...        ";
  {
    // Try setting the data using the generic function...
    iNumChars = sizeof(iBase_EntityHandle);
    char aChars[iNumChars];
    {
      iBase_EntityHandle* ppvData = new iBase_EntityHandle(aEnts[1]);
      char *pcData = (char*) ppvData;
      for (int i = 0; i < iNumChars; i++)
	aChars[i] = pcData[i];
      delete ppvData;
    }
    iMesh_setData(instance, aEnts[0], aTags[2], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getData(instance, aEnts[0], aTags[2], &aOutput, &iOutputAlloc,
		  &iNumChars, &err);
    CHECK_ERR(err);
    TESTEQ(static_cast<unsigned>(iNumChars), sizeof(iBase_EntityHandle));
    TESTEQ(*(iBase_EntityHandle* )aOutput, aEnts[1]);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    iMesh_setEHData(instance, aEnts[0], aTags[2], aEnts[2], &err);
    CHECK_ERR(err);

    iBase_EntityHandle entity;
    iMesh_getEHData(instance, aEnts[0], aTags[2], &entity, &err);
    CHECK_ERR(err);
    TESTEQ(entity, aEnts[2]);
    UT_FREE(aOutput);
  } // Entity handle tag
  cout << "DONE\n";

  cout << "Testing get and set methods with generic data (bytes)...  ";
  {
    // Try setting the data using the generic function...
    iNumChars = 6;
    char aChars[] = "abcdef";
    iMesh_setData(instance, aEnts[0], aTags[3], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    void *aOutput = NULL;
    int iOutputAlloc = 0;
    iMesh_getData(instance, aEnts[0], aTags[3], &aOutput, &iOutputAlloc,
		  &iNumChars, &err);
    CHECK_ERR(err);
    bool qGotData = (iNumChars == 6);
    qGotData = qGotData
	&& !strncmp(aChars, reinterpret_cast<char*>(aOutput), iNumChars);
    TEST(qGotData);
    UT_FREE(aOutput);
  } // Done with generic tests
  cout << "DONE\n";

  cout << "Testing tag finding...                                    ";
  int iGotTagsAlloc = 0, iNumTags;
  iMesh_getAllTags(instance, aEnts[0], &aGotTags, &iGotTagsAlloc, &iNumTags,
		   &err);
  CHECK_ERR(err);
  TEST(iNumTags == 4 + iNumInitTags);
  bool qTagsPresent = true;
  if (iNumTags == 4 + iNumInitTags) {
    iBase_TagHandle aAllTags[iNumTags];
    int ii;
    for (ii = 0; ii < 4; ii++)
      aAllTags[ii] = aTags[ii];
    for (; ii < 4 + iNumInitTags; ii++)
      aAllTags[ii] = aInitTags[ii - 4];
    GOT_ALL_ITEMS(aGotTags, aAllTags, iNumTags, qTagsPresent);
  }
  UT_FREE(aGotTags);
  UT_FREE(aInitTags);
  TEST(qTagsPresent);
  cout << "DONE\n";

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverTags(instance, pvRootSet, aTags, &(aEnts[0]));

    cout << "Testing get again after save/load..." << endl;
    cout << "Testing get with integers...              ";
    {
      int iData;
      iMesh_getIntData(instance, aEnts[0], aTags[0], &iData, &err);
      CHECK_ERR(err);
      TESTEQ(iData, 30);
    } // Integer tests
    cout << "DONE\n";

    cout << "Testing get and set methods with doubles...               ";
    {
      double dData;
      iMesh_getDblData(instance, aEnts[0], aTags[1], &dData, &err);
      CHECK_ERR(err);
      TESTEQ(dData, 1.9);
    } // Testing doubles
    cout << "DONE\n";

    // No obvious way to test entity handle tags after save/load...
    cout << "Testing get and set methods with generic data (bytes)...  ";
    {
      // Try getting the data using the generic function...
      iNumChars = 6;
      char aChars[] = "abcdef";

      void *aOutput = NULL;
      int iOutputAlloc = 0;
      iMesh_getData(instance, aEnts[0], aTags[3], &aOutput, &iOutputAlloc,
		    &iNumChars, &err);
      CHECK_ERR(err);
      bool qGotData = (iNumChars == 6);
      qGotData = qGotData
	  && !strncmp(aChars, reinterpret_cast<char*>(aOutput), iNumChars);
      TEST(qGotData);
      UT_FREE(aOutput);
    } // Done with generic tests
    cout << "DONE\n";
  }

  if (qErrTests) {
    // Create a new instance and some tags within it.
    iMesh_Instance new_instance;
    iMesh_newMesh("", &new_instance, &err, 0);

    for (int i = 0; i < 3; i++) {
      iMesh_createTag(new_instance, astrTagNames[i], 1, aTagTypes[i],
		      &aTags[i + 4], &err, aiTagNameLen[i]);
      CHECK_ERR(err);
    }
    cout << "Testing error throwing for methods...              ";
    // Error throwing tests simply try to create a situation where the interface
    // should, according to the spec, throw an error, and verify that:
    // (instance,a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reversal of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice versa, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.

    // Setting tag data for a mesh with no data
    char aData[10];
    iMesh_setData(new_instance, aEnts[0], aTags[4], aData, iNumChars, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // Setting tag data for a mesh with no data
    iMesh_setIntData(new_instance, aEnts[0], aTags[4], 15, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    //  Setting tag data for a mesh with no data
    iMesh_setDblData(new_instance, aEnts[0], aTags[5], 1.5, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // Setting  tag data for a mesh with no data
    iMesh_setEHData(new_instance, aEnts[0], aTags[6], aEnts[1], &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // Getting tag data for a mesh with no data
    void *aOutput;
    int iOutputAlloc = 0;
    iMesh_getData(new_instance, aEnts[0], aTags[4], &aOutput, &iOutputAlloc,
		  &iNumChars, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    int iData;
    // Getting tag data for a mesh with no data
    iMesh_getIntData(new_instance, aEnts[0], aTags[4], &iData, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    double dData;
    // Getting tag data for a mesh with no data
    iMesh_getDblData(new_instance, aEnts[0], aTags[5], &dData, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    iBase_EntityHandle EH;
    // Getting tag data for a mesh with no data
    iMesh_getEHData(new_instance, aEnts[0], aTags[6], &EH, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // Removing tag data for a mesh with no data
    iMesh_rmvTag(new_instance, aEnts[0], aTags[4], &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    for (int i = 0; i < 3; i++) {
      iMesh_destroyTag(new_instance, aTags[i + 4], true, &err);
      CHECK_ERR(err);
    }
    iMesh_dtor(new_instance, &err);
    CHECK_ERR(err);

    // Back to using the original tags on the original mesh...
    // setting  integer tag data for a double tag handle
    iMesh_setIntData(instance, aEnts[0], aTags[1], 15, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // Setting double tag data for an integer tag handle
    iMesh_setDblData(instance, aEnts[0], aTags[0], 1.5, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    if (!qSaveLoad) {
      // If you've done save/load, then you'd need to re-establish the
      //  identify of the entity in question.  Yuck.  If not, this
      //  function could throw either of two errors.

      // Setting entity handle tag data for an integer tag handle
      iMesh_setEHData(instance, aEnts[0], aTags[0], aEnts[1], &err);
      CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);
    }

    // Getting integer tag data for a double tag handle
    iMesh_getIntData(instance, aEnts[0], aTags[1], &iData, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // Getting double tag data for an integer tag handle
    iMesh_getDblData(instance, aEnts[0], aTags[0], &dData, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // Getting entity handle tag data for an integer tag handle
    iMesh_getEHData(instance, aEnts[0], aTags[0], &EH, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    cout << "DONE\n";
  }

  cout << "Testing tag removal and forced/unforced destruction...    ";
  {
    bool qOK = true;
    // Untag the entity for the first two tags...
    iMesh_rmvTag(instance, aEnts[0], aTags[0], &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_rmvTag(instance, aEnts[0], aTags[1], &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    // ... and see if they can be destroyed without an error
    iMesh_destroyTag(instance, aTags[0], false, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[1], false, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    if (qOK)
      iNPassed++;
  }
  // Try to destroy the other tags without first untagging the entity
  // Trying to destroy tag with data without using 'forced'
  iMesh_destroyTag(instance, aTags[2], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  // Trying to destroy tag with data without using 'forced'
  iMesh_destroyTag(instance, aTags[3], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  //
  {
    bool qOK = true;
    // Force the destruction of the tags
    iMesh_destroyTag(instance, aTags[2], true, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[3], true, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    if (qOK)
      iNPassed++;
  }

  cout << "DONE\n";
}
