#include "iMesh_unitTest.hh"

void
vTestArrTagInterface(iMesh_Instance& instance)
{
  // Not testing save/load in this function; if entity tags are handled
  // properly with save/load (as tested in the entity tag test), they'll
  // be fine here, just a lot harder to test.

  //     getArrData, getIntArrData, getDblArrData, getEHArrData,
  //     setArrData, setIntArrData, setDblArrData, setEHArrData,
  //     rmvArrTag
  cout << "Testing basic array tagging functionality...\n";
  iBase_TagHandle aTags[7];
  int iOutputAlloc = 9 * sizeof(double);
  void *aOutput = calloc(iOutputAlloc, 1);
  iBase_EntityHandle aEnts[3];
  int aints[9];
  double aDoubles[9];
  iBase_EntityHandle aHandles[9];

  int iNumChars, iNumVals, err;

  // Begin the "normal" operation testing.

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

  // Create tags of each type that store a single value
  for (int i = 0; i < 3; i++) {
    iMesh_createTag(instance, astrTagNames[i], 3, aTagTypes[i], &aTags[i], &err,
		    aiTagNameLen[i]);
    CHECK_ERR(err);
  }
  // ... and a generic tag for completeness
  iMesh_createTag(instance, astrTagNames[3], 6, aTagTypes[3], &aTags[3], &err,
		  aiTagNameLen[3]);
  CHECK_ERR(err);

  cout << "Testing get and set methods with integers...              ";
  bool qGotData;
  {
    // Try setting the data using the generic function...
    iNumChars = 3 * 3 * sizeof(int);
    int aiData[] =
      { 4, 8, 12, 16, 20, 24, 28, 32, 36 };
    char *aChars = (char*) (aiData);

    iMesh_setArrData(instance, aEnts, 3, aTags[0], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    iMesh_getArrData(instance, aEnts, 3, aTags[0], &aOutput, &iOutputAlloc,
		     &iNumChars, &err);
    CHECK_ERR(err);
    qGotData = (iNumChars == 3 * 3 * (int) sizeof(int));
    if (qGotData) {
      int *piOutput = (int*) aOutput;
      for (int i = 0; i < 9; i++)
	qGotData = qGotData && (piOutput[i] == 4 * (i + 1));
    }
    TEST(qGotData);

    // Now try to change the data using the specific function
    // and ensure it comes back intact
    for (int i = 0; i < 9; i++)
      aints[i] = 5 * (i + 1);
    iMesh_setIntArrData(instance, aEnts, 3, aTags[0], aints, 9, &err);
    CHECK_ERR(err);
    int *aintsBack = NULL, iintsAlloc = 0;
    iMesh_getIntArrData(instance, aEnts, 3, aTags[0], &aintsBack, &iintsAlloc,
			&iNumVals, &err);
    CHECK_ERR(err);
    qGotData = (iNumVals == 9);
    if (qGotData) {
      for (int i = 0; i < 9; i++)
	qGotData = qGotData && (aintsBack[i] == 5 * (i + 1));
    }
    TEST(qGotData);
    UT_FREE(aintsBack);
  }
  cout << "DONE\n"; // Testing array int tags.

  cout << "Testing get and set methods with doubles...               ";
  {
    // Try setting the data using the generic function...
    iNumChars = 9 * sizeof(double);
    double adData[] =
      { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    char *aChars = (char*) (adData);
    iMesh_setArrData(instance, aEnts, 3, aTags[1], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact
    iNumChars = 0;
    iMesh_getArrData(instance, aEnts, 3, aTags[1], &aOutput, &iOutputAlloc,
		     &iNumChars, &err);
    CHECK_ERR(err);
    qGotData = (iNumChars == 9 * (int) sizeof(double));
    if (qGotData) {
      double *pdOutput = (double*) aOutput;
      for (int i = 0; i < 9; i++)
	qGotData = qGotData
	    && APPROX_EQUAL(pdOutput[i], ((double )(i + 1)) / 10);
    }
    TEST(qGotData);
    // Now try to change the data using the specific function
    // and ensure it comes back intact
    for (int i = 0; i < 9; i++)
      aDoubles[i] = ((double) (i + 1)) / 5;
    iMesh_setDblArrData(instance, aEnts, 3, aTags[1], aDoubles, 9, &err);
    CHECK_ERR(err);
    double *aDoublesBack = NULL;
    int iDoublesBackAlloc = 0;
    iMesh_getDblArrData(instance, aEnts, 3, aTags[1], &aDoublesBack,
			&iDoublesBackAlloc, &iNumVals, &err);
    CHECK_ERR(err);
    qGotData = (iNumVals == 9);
    if (qGotData) {
      for (int i = 0; i < 9; i++)
	qGotData = qGotData
	    && APPROX_EQUAL(aDoublesBack[i], ((double )(i + 1)) / 5);
    }
    TEST(qGotData);
    UT_FREE(aDoublesBack);
  }
  cout << "DONE\n";

  cout << "Testing get and set methods with entity handles...        ";
  {
    // Try setting the data using the generic function...
    iNumChars = 9 * sizeof(iBase_EntityHandle);

    iBase_EntityHandle aEntTemp[9];
    for (int i = 0; i < 9; i++)
      aEntTemp[i] = aEnts[i % 3];
    char *aChars = (char*) aEntTemp;

    iMesh_setArrData(instance, aEnts, 3, aTags[2], aChars, iNumChars, &err);
    CHECK_ERR(err);

    // ... then seeing if it comes back intact
    iMesh_getArrData(instance, aEnts, 3, aTags[2], &aOutput, &iOutputAlloc,
		     &iNumChars, &err);
    CHECK_ERR(err);
    qGotData = (iNumChars == 9 * (int) sizeof(iBase_EntityHandle));
    if (qGotData) {
      iBase_EntityHandle *aEntOut = (iBase_EntityHandle*) aOutput;
      for (int i = 0; i < 9; i++)
	qGotData = qGotData && (aEntOut[i] == aEntTemp[i]);
    }
    TEST(qGotData);

    // Now try to change the data using the specific function
    // and ensure it comes back intact
    for (int i = 0; i < 9; i++)
      aHandles[i] = aEnts[(i + 1) % 3];

    iMesh_setEHArrData(instance, aEnts, 3, aTags[2], aHandles, 9, &err);
    CHECK_ERR(err);
    iBase_EntityHandle *aEntOut = NULL;
    int iEntOutAlloc = 0;
    iMesh_getEHArrData(instance, aEnts, 3, aTags[2], &aEntOut, &iEntOutAlloc,
		       &iNumVals, &err);
    CHECK_ERR(err);

    qGotData = (iNumVals == 9);
    if (qGotData) {
      for (int i = 0; i < 9; i++)
	qGotData = qGotData && (aEntOut[i] == aEnts[(i + 1) % 3]);
    }
    TEST(qGotData);
    UT_FREE(aEntOut);
  }
  cout << "DONE\n"; // with handle tag testing

  cout << "Testing get and set methods with generic data (bytes)...  ";
  {
    // Try setting the data using the generic function...
    iNumChars = 18;
    char aChars[] = "abcdefghijklmnopqr";

    iMesh_setArrData(instance, aEnts, 3, aTags[3], aChars, iNumChars, &err);
    CHECK_ERR(err);
    // ... then seeing if it comes back intact

    iMesh_getArrData(instance, aEnts, 3, aTags[3], &aOutput, &iOutputAlloc,
		     &iNumChars, &err);
    CHECK_ERR(err);
    qGotData = (iNumChars == 18)
	&& !strncmp(aChars, reinterpret_cast<char*>(aOutput), iNumChars);
    TEST(qGotData);
  }
  cout << "DONE\n"; // with generic tag testing

  // Testing error throwing
  if (qErrTests) {
    iMesh_Instance new_instance;
    iMesh_newMesh("", &new_instance, &err, 0);

    cout << "Testing error throwing for array tag methods...              ";
    // Error throwing tests simply try to create a situation where the intefrace
    // should, according to the spec, throw an error, and verify that:
    // (a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reveral of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice vera, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.

    for (int i = 0; i < 3; i++) {
      iMesh_createTag(new_instance, astrTagNames[i], 3, aTagTypes[i],
		      &aTags[i + 4], &err, aiTagNameLen[i]);
      CHECK_ERR(err);
    }

    // Setting tag data for a mesh with no data
    iMesh_setArrData(new_instance, aEnts, 3, aTags[4], aOutput, iNumChars,
		     &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // setting tag data for a mesh with no data
    iMesh_setIntArrData(new_instance, aEnts, 3, aTags[4], aints, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // setting tag data for a mesh with no data
    iMesh_setDblArrData(new_instance, aEnts, 3, aTags[5], aDoubles, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // setting tag data for a mesh with no data
    iMesh_setEHArrData(new_instance, aEnts, 3, aTags[6], aHandles, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    //  getting tag data for a mesh with no data
    iMesh_getArrData(new_instance, aEnts, 3, aTags[4], &aOutput, &iOutputAlloc,
		     &iNumChars, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    int *aIntsOut;
    int iIntsAlloc = 0;
    // getting tag data for a mesh with no data
    iMesh_getIntArrData(new_instance, aEnts, 3, aTags[4], &aIntsOut,
			&iIntsAlloc, &iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // getting tag data for a  mesh with no data
    double *aDoublesOut;
    int iDoublesAlloc = 0;
    iMesh_getDblArrData(new_instance, aEnts, 3, aTags[5], &aDoublesOut,
			&iDoublesAlloc, &iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // getting tag data for a mesh with no data
    iBase_EntityHandle *aHandlesOut = NULL;
    int iHandlesAlloc = 0;
    iMesh_getEHArrData(new_instance, aEnts, 3, aTags[6], &aHandlesOut,
		       &iHandlesAlloc, &iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_HANDLE);

    // setting integer tag data for a double tag handle
    iMesh_setIntArrData(instance, aEnts, 3, aTags[1], aints, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // setting double tag data for an integer tag handle
    iMesh_setDblArrData(instance, aEnts, 3, aTags[0], aDoubles, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // setting entity handle tag data for integer tag handle
    iMesh_setEHArrData(instance, aEnts, 3, aTags[0], aHandles, 9, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // setting tag data for an invalid input array size
    iMesh_setIntArrData(instance, aEnts, 3, aTags[0], aints, 12, &err);
    CHECK_ERR2(err, iBase_BAD_ARRAY_SIZE);

    // Excercise error getting functions...
    int err2, err3;
#ifdef IMESH_VERSION_GE
#if IMESH_VERSION_GE(1,3,0)
    iMesh_getErrorType(instance, &err2);
#else
    iMesh_getErrorType(instance, &err2, &err3);
#endif
#else
    iMesh_getErrorType(instance, &err2, &err3);
#endif
    CHECK_ERR2(err2, iBase_BAD_ARRAY_SIZE);

    char descr[256];
    descr[0] = '\0';
#ifdef IMESH_VERSION_GE
#if IMESH_VERSION_GE(1,3,0)
    iMesh_getDescription(instance, descr, (int) sizeof(descr));
#else
    iMesh_getDescription(instance, descr, &err3, (int) sizeof(descr));
#endif
#else
    iMesh_getDescription(instance, descr, &err3, (int) sizeof(descr));
#endif
    TEST(descr != "");

    //  setting double tag data for an invalid input array size
    iMesh_setDblArrData(instance, aEnts, 3, aTags[1], aDoubles, 12, &err);
    CHECK_ERR2(err, iBase_BAD_ARRAY_SIZE);

    //  setting entity handle tag data for an invalid input array size
    iMesh_setEHArrData(instance, aEnts, 2, aTags[2], aHandles, 9, &err);
    CHECK_ERR2(err, iBase_BAD_ARRAY_SIZE);

    // getting integer tag data for a double tag handle
    iMesh_getIntArrData(instance, aEnts, 3, aTags[1], &aIntsOut, &iIntsAlloc,
			&iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // getting double tag data for an integer tag handle
    iMesh_getDblArrData(instance, aEnts, 3, aTags[0], &aDoublesOut,
			&iDoublesAlloc, &iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    // getting entity handle tag data for an integer tag handle
    iMesh_getEHArrData(instance, aEnts, 3, aTags[0], &aHandlesOut,
		       &iHandlesAlloc, &iNumVals, &err);
    CHECK_ERR2(err, iBase_INVALID_TAG_HANDLE);

    iMesh_dtor(new_instance, &err);
    cout << "DONE\n";
  } // Done with error checking
  cout << "Testing tag removal and forced/unforced destruction...    ";
  {
    bool qOK = true;
    iMesh_rmvArrTag(instance, aEnts, 3, aTags[0], &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_rmvArrTag(instance, aEnts, 3, aTags[1], &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[0], false, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[1], false, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    if (qOK)
      iNPassed++;
  }
  // trying to destroy tag with data without using forced
  iMesh_destroyTag(instance, aTags[2], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  iMesh_destroyTag(instance, aTags[3], false, &err);
  CHECK_ERR2(err, iBase_TAG_IN_USE);

  {
    bool qOK = true;
    iMesh_destroyTag(instance, aTags[2], true, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iMesh_destroyTag(instance, aTags[3], true, &err);
    qOK = qOK && (err == iBase_SUCCESS);
    CHECK_ERR(err);

    iNPassed++;
  }
  cout << "DONE\n";
  free(aOutput);
}
