//@ File header

#ifndef CREATE_ARG
#define CREATE_ARG ""
#endif

// This file is intended to test functionality of a iMesh implementation.
// Mesh query, modification, iterators, entity sets, and tags are all
// tested.  In cases where it's practical to define a priori a correct
// result for a query, that is done.  In other cases, consistency is the
// only thing checked.  In still others, the only possible correctness
// check is that the function took correct input and didn't choke
// (things like getAdjTable or getGeometricDim, for instance).

// It is worth noting that, although efforts have been made to make each
// test as independent as possible, it is unavoidable that a failure
// during one test may cause failures in subsequent tests, sometimes
// even across different interface tests.  In light of this, it is best
// to try to correct errors in the order they occur, or at least to
// exclude interface tests that cause errors to see if later errors are
// avoided.

// Modifications:
//    Mark C. Miller, Thu Aug 12 09:53:51 PDT 2010
//    Added Save/Load option and randomly sprinkled DoSaveLoad calls
//    throughout the tests to provide a means to ensure relevant
//    state is preserved across _save and _load operations.
//
//    Mark C. Miller, Mon Aug 16 09:37:27 PDT 2010
//    Deal with fact that iMesh_unitTest maintains a handle on the
//    root set and that may need to change during Save/Load. So, it
//    was added to vDoSaveLoad arglist.
//
//    Mark C. Miller, Wed Dec  1 09:51:42 PST 2010
//    Replaced use of '/tmp' with a directory in . named with pid
//    of process. Note that if iMesh_unitTest does not exit 'normally'
//    the resulting temporary directories can be left around.
//
//    Mark C. Miller, Wed Dec  8 00:01:49 PST 2010
//    Turn off save-load-clean if we are not even doing save-load.

#define MAIN_PROGRAM
#include "iMesh_unitTest.hh"
#include <omp.h>
#include <dirent.h>
#include <errno.h>
#include <limits>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <vector>

static char TMPDIR[64];

void
qFailed()
{
  iNFailed++;
}

void
qIssueLoadMessage(const char *filename, int err)
{
  if (err == iBase_SUCCESS)
    return;

  cout << endl;
  cout << "Unable to load file \"" << filename << "\"" << endl;
  cout << endl;

  int len = strlen(filename);
  bool isVTK = false;

  if (filename[len - 4] == '.' && filename[len - 3] == 'v'
      && filename[len - 2] == 't' && filename[len - 1] == 'k')
    isVTK = true;

  if (isVTK) {
    cout << "The file name extension suggests this is a VTK file.\n"
	"For testing purposes, all 'compliant' implementations\n"
	"are required to read VTK files. However, VTK files come\n"
	"in many flavors and may contain a lot of optional things.\n"
	"Therefore, for testing purposes, the flavor of VTK files\n"
	"an iMesh implementation is required to be able to read is\n"
	"constrained to the following:\n"
	"   - Must be ASCII format\n"
	"   - Must be DATASET UNSTRUCTURED_GRID\n"
	"   - May NOT contain optional VTK data\n"
	"   - Must be float data type\n"
	"   - Must have >55 verticies and >10 elements\n"
	"   - Must be 2D or 3D\n"
	"   - For 3D: Must consist only of Tets, Pyrmaids, Prisms and/or Hexs\n"
	"   - For 2D: Must consist only of Triangles and/or Quads\n";
    cout << endl;
  }
  else {
    cout << "The file name extension suggests this is NOT a VTK file.\n"
	"Not all implementations read the same file formats.\n"
	"You may wish to independently confirm that the iMesh\n"
	"implementation you have linked to can indeed read this file.\n"
	"In addition, also be aware that for purposes of testing\n"
	"the input mesh must have >55 vertices and >10 elements\n";
    cout << endl;
  }
}

bool
qCheckNumOfType(iMesh_Instance instance, iBase_EntitySetHandle SH, int iType,
		int iNumEntsExpected)
{
  int err, iNumEnts;
  iMesh_getNumOfType(instance, SH, iType, &iNumEnts, &err);
  CHECK_ERR(err);
  return (iNumEnts == iNumEntsExpected);
}

bool
qCheckNumOfTopo(iMesh_Instance instance, iBase_EntitySetHandle SH, int iTopo,
		int iNumEntsExpected)
{
  int err, iNumEnts;
  iMesh_getNumOfTopo(instance, SH, iTopo, &iNumEnts, &err);
  CHECK_ERR(err);
  return (iNumEnts == iNumEntsExpected);
}

void
vCheckEntities(iMesh_Instance instance, iBase_EntitySetHandle SH, int iType,
	       int iTopo, iBase_EntityHandle aEHExpected[],
	       int iNumEntsExpected)
{
  int err;
  iBase_EntityHandle *aGotEnts = NULL;
  int iNumEnts;
  int iGotEntsAlloc = 0;
  iMesh_getEntities(instance, SH, iType, iTopo, &aGotEnts, &iGotEntsAlloc,
		    &iNumEnts, &err);
  CHECK_ERR(err);

  int iIsList;
  iMesh_isList(instance, SH, &iIsList, &err);
  CHECK_ERR(err);

  // This used to be the test:  did you get precisely the entities
  // you expected?  Can't do this in the presence of save-load
  // pairs.
  if (!qSaveLoad) {
    bool qGotAllEnts = true;
    if (iNumEnts == iNumEntsExpected) {
      if (iIsList) {
	for (int ii = 0; ii < iNumEntsExpected; ii++) {
	  qGotAllEnts = qGotAllEnts && (aGotEnts[ii] == aEHExpected[ii]);
	}
      }
      else {
	GOT_ALL_ITEMS(aEHExpected, aGotEnts, iNumEntsExpected, qGotAllEnts);
      }
    }
    TEST((iNumEnts == iNumEntsExpected) && qGotAllEnts);
  }

  // The following tests are weaker, but work even if you've done a
  // save-load pair.
  bool qTypesOK = true, qToposOK = true;
  if (iType != iBase_ALL_TYPES) {
    // No point in testing if anything will pass.
    int *aiEType = NULL, type_alloc = 0, type_size;
    iMesh_getEntArrType(instance, aGotEnts, iNumEnts, &aiEType, &type_alloc,
			&type_size, &err);
    qTypesOK = (iNumEnts == type_size);
    for (int ii = 0; ii < type_size; ii++) {
      qTypesOK = qTypesOK && (aiEType[ii] == iType);
    }
    UT_FREE(aiEType);
  }
  if (iTopo != iMesh_ALL_TOPOLOGIES) {
    // No point in testing if anything will pass.
    int *aiETopo = NULL, topo_alloc = 0, topo_size;
    iMesh_getEntArrTopo(instance, aGotEnts, iNumEnts, &aiETopo, &topo_alloc,
			&topo_size, &err);
    qToposOK = (iNumEnts == topo_size);
    for (int ii = 0; ii < topo_size; ii++) {
      qToposOK = qToposOK && (aiETopo[ii] == iTopo);
    }
    UT_FREE(aiETopo);
  }
  UT_FREE(aGotEnts);

  TEST(qTypesOK && qToposOK);
}

void
vCheckBoolOpResult(iMesh_Instance instance, iBase_EntitySetHandle SH,
		   int iShouldBeList, iBase_EntitySetHandle aSHExpected[],
		   int iNumSetsExpected, iBase_EntityHandle aEntsExpected[],
		   int iNumEntsExpected)
{
  int err;
  int iIsList;
  iMesh_isList(instance, SH, &iIsList, &err);
  CHECK_ERR(err);
  TESTEQ(iIsList, iShouldBeList);

  if (iNumSetsExpected != -1) {
    iBase_EntitySetHandle *aGotSets = NULL;
    int iGotSetsAlloc = 0, iNumSets;
    iMesh_getEntSets(instance, SH, 0, &aGotSets, &iGotSetsAlloc, &iNumSets,
		     &err);
    CHECK_ERR(err);
    TESTEQ(iNumSets, iNumSetsExpected);

    bool qGotAllSets = true;
    if (iNumSets == iNumSetsExpected) {
      GOT_ALL_ITEMS(aSHExpected, aGotSets, iNumSetsExpected, qGotAllSets);
      TEST(qGotAllSets);
    }
    UT_FREE(aGotSets);
  }

  int iNumEnts;
  iBase_EntityHandle * aGotEnts = NULL;
  int iGotEntsAlloc = 0;
  iMesh_getEntities(instance, SH, iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
		    &aGotEnts, &iGotEntsAlloc, &iNumEnts, &err);
  CHECK_ERR(err);
  TESTEQ(iNumEnts, iNumEntsExpected);
  bool qGotAllEnts = (iNumEnts == iNumEntsExpected);
  if (qGotAllEnts) {
    if (iShouldBeList) {
      for (int ii = 0; ii < iNumEntsExpected; ii++) {
	qGotAllEnts = qGotAllEnts && (aGotEnts[ii] == aEntsExpected[ii]);
	TESTEQ(aGotEnts[ii], aEntsExpected[ii]);
      }
    }
    else {
      GOT_ALL_ITEMS(aEntsExpected, aGotEnts, iNumEntsExpected, qGotAllEnts);
    }
    TEST(qGotAllEnts);
  }
  UT_FREE(aGotEnts);
}

void
CleanSaveLoad()
{
  DIR *dir;
  struct dirent *dent;
  dir = opendir(TMPDIR);
  if (chdir(TMPDIR) == 0) {
    while (dir && (dent = readdir(dir))) {
      string dname = string(dent->d_name);
      if (dname == "." || dname == "..")
	continue;
      unlink(dname.c_str());
    }
    if (chdir("..") != 0)
      std::cerr << "Unable to cleanup temp data dir" << endl;
  }
  closedir(dir);
  rmdir(TMPDIR);
}

void
CleanSaveLoadOnSignal(int sig)
{
  CleanSaveLoad();
  exit(sig);
}

//
// Test save/load to ensure that the two operations result in
// a valid and expected iMesh_Instance object where relevant
// state is preserved. If a run of iMesh_unitTest can succeed
// without save/loads but fails with them, that means that
// the implementation's save/load operations are not working
// in that they fail to preserve the state iMesh_unitTest is
// expecting and testing for.
//
void
vDoSaveLoad(iMesh_Instance *instancep, iBase_EntitySetHandle *rootSet)
{
  static int numcalls = 0;
  int imErr;
  iMesh_Instance instance = *instancep;
  iMesh_Instance newinst;
  iBase_EntitySetHandle newroot;
  DIR *dir;
  struct dirent *dent;
  bool foundIt = false;

  if (!qSaveLoad)
    return;

//     cout << " ***Save/Load...";
  char imeshBasename[128];
  char imeshFilename[128];
  snprintf(imeshBasename, sizeof(imeshFilename), "iMesh_unitTest_data_%d",
	   numcalls);
  snprintf(imeshFilename, sizeof(imeshFilename), "%s/%s", TMPDIR,
	   imeshBasename);
  int len = strlen(imeshFilename);
  char options[1024];
  snprintf(options, sizeof(options), "grummp:silent grummp:vtk refimpl:silent");
  iMesh_save(*instancep, *rootSet, imeshFilename, options, &imErr, len,
	     strlen(options));
  CHECK_ERR(imErr);

  // Lets do a tiny bit of work to ensure we've actually got a file on disk now
  dir = opendir(TMPDIR);
  while (dir && (dent = readdir(dir))) {
    struct stat statbuf;
    errno = 0;
    memset(&statbuf, 0, sizeof(statbuf));
    if (strncmp(imeshBasename, dent->d_name, strlen(imeshBasename)) == 0) {
      char tmpname[256];
      snprintf(tmpname, sizeof(tmpname), "%s/%s", TMPDIR, dent->d_name);
      if (stat(tmpname, &statbuf) == 0 && errno == 0 && statbuf.st_size > 0)
	foundIt = true;
    }
  }
  closedir(dir);
  if (!foundIt) {
    CHECK_ERR(iBase_FILE_WRITE_ERROR);
    return;
  }

  // Ok, now try to read this file into a new instance and root set
  iMesh_newMesh("", &newinst, &imErr, 0);
  iMesh_getRootSet(newinst, &newroot, &imErr);
  iMesh_load(newinst, newroot, imeshFilename, options, &imErr, len,
	     strlen(options));
  if (imErr == iBase_SUCCESS) {
    // Destroy the old mesh now
    iMesh_dtor(*instancep, &imErr);
    CHECK_ERR(imErr);

    *instancep = newinst;
    *rootSet = newroot;
  }
  else {
    qIssueLoadMessage(imeshFilename, imErr);
    CHECK_ERR(imErr);
  }
//     cout << " Save/Load-DONE*** ";

  numcalls++;
}

void
vRecoverSetsContain(iMesh_Instance instance, iBase_EntitySetHandle pvRootSet,
		    iBase_EntitySetHandle apvEntSets[])
{
  // The save-load cycle probably invalidated our set handles.  Bummer.
  // So now we need to re-identify the sets.

  iBase_EntitySetHandle *aSets = NULL;
  int iTotalRet, iSetsAlloc = 0, err;
  iMesh_getEntSets(instance, pvRootSet, -1, &aSets, &iSetsAlloc, &iTotalRet,
		   &err);
  // Once sets are identified, these should be unique integers.
  int aiSetNum[] =
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

  // The complicated bit is that there may be other sets defined than
  // the ten that were just created.  Also, there's no guarantee at
  // all about the order in which sets are listed in the return data.
  // Distinguishing characteristics:
  // Set   Sets cont  Cont by   Set/List
  // 0     3/9         ---         L
  // 1     2/2          0          S
  // 2     0/0          0          L
  // 3     2/4          0          S
  // 4     0/0          1          L
  // 5     0/0          1          S
  // 6     2/2          3          L
  // 7     0/0          3          S
  // 8     0/0          6          L
  // 9     0/0          6          S

  // This code assumes that there won't be any remarkable coincidences
  // in the sets that the implementation auto-created.

  for (int ii = 0; ii < 10; apvEntSets[ii++] = NULL)
    ;

  std::vector<int> vecNoDescend;
  for (int iSet = 0; iSet < iTotalRet; iSet++) {
    int iNumSets = -1;
    iMesh_getNumEntSets(instance, aSets[iSet], -1, &iNumSets, &err);
    CHECK_ERR(err);

    switch (iNumSets)
      {
      case 0:
	vecNoDescend.push_back(iSet);
	break;
      case 2:
	int iIsList;
	iMesh_isList(instance, aSets[iSet], &iIsList, &err);
	CHECK_ERR(err);
	if (iIsList) {
	  apvEntSets[6] = aSets[iSet];
	}
	else {
	  apvEntSets[1] = aSets[iSet];
	}
	break;
      case 4:
	apvEntSets[3] = aSets[iSet];
	break;
      case 9:
	apvEntSets[0] = aSets[iSet];
	break;
      default:
	// Do nothing with these; they aren't ours.
	break;
      }
  }
  TEST(vecNoDescend.size() >= 6);
  for (int i = vecNoDescend.size() - 1; i >= 0; i--) {
    int iIsCont, iIsList;
    iBase_EntitySetHandle es = aSets[vecNoDescend[i]];
    iMesh_isList(instance, es, &iIsList, &err);
    CHECK_ERR(err);

    iMesh_isEntSetContained(instance, apvEntSets[0], es, &iIsCont, &err);
    CHECK_ERR(err);
    if (iIsCont) {
      TEST(iIsList);
      apvEntSets[2] = es;
    }
    else {
      iMesh_isEntSetContained(instance, apvEntSets[3], es, &iIsCont, &err);
      CHECK_ERR(err);
      if (iIsCont) {
	TEST(!iIsList);
	apvEntSets[7] = es;
      }
      else {
	iMesh_isEntSetContained(instance, apvEntSets[1], es, &iIsCont, &err);
	CHECK_ERR(err);
	if (iIsCont) {
	  if (iIsList)
	    apvEntSets[4] = es;
	  else
	    apvEntSets[5] = es;
	}
	else {
	  iMesh_isEntSetContained(instance, apvEntSets[6], es, &iIsCont, &err);
	  CHECK_ERR(err);
	  if (iIsCont) {
	    if (iIsList)
	      apvEntSets[8] = es;
	    else
	      apvEntSets[9] = es;
	  }
	} // Fell through to checking for membership in set 6
      } // Fell through to checking for membership in sets 1 or 6
    } // Fell through to checking for membership in sets 1, 3, or 6
  } // Checking what to do with each set.
  for (int ii = 0; ii < 10; ii++) {
    TEST(apvEntSets[ii] != NULL);
  }
  UT_FREE(aSets);
}

void
vRecoverSetsParentChild(iMesh_Instance instance,
			iBase_EntitySetHandle pvRootSet,
			iBase_EntitySetHandle apvEntSets[])
{
  // The save-load cycle probably invalidated our set handles.  Bummer.
  // So now we need to re-identify the sets.

  iBase_EntitySetHandle *aSets = NULL;
  int iTotalRet, iSetsAlloc = 0, err;
  iMesh_getEntSets(instance, pvRootSet, -1, &aSets, &iSetsAlloc, &iTotalRet,
		   &err);
  // Once sets are identified, these should be unique integers.
  int aiSetNum[] =
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

  // The complicated bit is that there may be other sets defined than
  // the ten that were just created.  Also, there's no guarantee at
  // all about the order in which sets are listed in the return data.
  // Distinguishing characteristics:
  // Set   TotChild  TotPar ImmPar  Set/List
  // 0     9          0       ---   L
  // 1     2          1        0    S
  // 2     0          1        0    L
  // 3     4          1        0    S
  // 4     0          2        1    L
  // 5     0          2        1    S
  // 6     2          2        3    L
  // 7     0          2        3    S
  // 8     0          3        6    L
  // 9     0          3        6    S

  // This code assumes that there won't be any remarkable coincidences
  // in the sets that the implementation auto-created.

  for (int ii = 0; ii < 10; apvEntSets[ii++] = NULL)
    ;

  std::vector<int> vecNoDescend;
  for (int iSet = 0; iSet < iTotalRet; iSet++) {
    int iNumChildren = -1;
    iMesh_getNumChld(instance, aSets[iSet], -1, &iNumChildren, &err);
    CHECK_ERR(err);

    switch (iNumChildren)
      {
      case 0:
	vecNoDescend.push_back(iSet);
	break;
      case 2:
	{
	  int iNumParents = -1;
	  iMesh_getNumPrnt(instance, aSets[iSet], -1, &iNumParents, &err);
	  CHECK_ERR(err);

	  if (iNumParents == 2) {
	    apvEntSets[6] = aSets[iSet];
	  }
	  else {
	    apvEntSets[1] = aSets[iSet];
	  }
	}
	break;
      case 4:
	apvEntSets[3] = aSets[iSet];
	break;
      case 9:
	apvEntSets[0] = aSets[iSet];
	break;
      default:
	// Do nothing with these; they aren't ours.
	break;
      }
  }
  TEST(vecNoDescend.size() >= 6);
  for (int i = vecNoDescend.size() - 1; i >= 0; i--) {
    int iIsList;
    iBase_EntitySetHandle es = aSets[vecNoDescend[i]];
    iMesh_isList(instance, es, &iIsList, &err);
    CHECK_ERR(err);

    iBase_EntitySetHandle *aesParent = NULL;
    int aesParentAlloc = 0, aesParentSize;

    iMesh_getPrnts(instance, es, 0, &aesParent, &aesParentAlloc, &aesParentSize,
		   &err);
    CHECK_ERR(err);
    if (aesParentSize != 1)
      continue; // This isn't one of ours.

    int iNumParents = -1;
    iMesh_getNumPrnt(instance, es, -1, &iNumParents, &err);
    CHECK_ERR(err);

    switch (iNumParents)
      {
      case 1:
	// There's only one of these.
	TEST(iIsList);
	apvEntSets[2] = es;
	break;
      case 2:
	if (iIsList) {
	  // Only one of these.
	  apvEntSets[4] = es;
	}
	else {
	  if (aesParent[0] == apvEntSets[1]) {
	    apvEntSets[5] = es;
	  }
	  else {
	    TESTEQ(aesParent[0], apvEntSets[3]);
	    apvEntSets[7] = es;
	  }
	}
	break;
      case 3:
	if (iIsList) {
	  apvEntSets[8] = es;
	}
	else {
	  apvEntSets[9] = es;
	}
	break;
      default:
	// Bad things are about to happen, somewhere...
	break;
      }
    UT_FREE(aesParent);
  } // Checking what to do with each set.
  for (int ii = 0; ii < 10; ii++) {
    TEST(apvEntSets[ii] != NULL);
  }
  UT_FREE(aSets);
}

void
vRecoverSetsFromEntCounts(iMesh_Instance instance,
			  iBase_EntitySetHandle pvRootSet,
			  iBase_EntitySetHandle *apvSets,
			  const int aiNumTypeEnts[], const int aiNumTopoEnts[])
{
  // In the original setup, we had:
  // Set #     List?   Matching type/topo
  //   0         Y         (Empty)
  //   1         Y         Type
  //   2         Y         Topo
  //   3         N         Type
  //   4         N         Topo
  //
  // Set 0 could actually be masqueraded for fairly easily, but that's
  // probably okay....
  int i, err;
  int iTotalEntsExpectedType = 0, iTotalEntsExpectedTopo = 0;
  for (i = 0; i < 5; i++) {
    apvSets[i] = NULL;
  }
  for (i = 0; i < iBase_ALL_TYPES; i++) {
    iTotalEntsExpectedType += aiNumTypeEnts[i];
  }
  for (i = 0; i < iMesh_ALL_TOPOLOGIES; i++) {
    iTotalEntsExpectedTopo += aiNumTopoEnts[i];
  }

  iBase_EntitySetHandle *allSets = NULL;
  int allSets_alloc = 0, allSets_size;
  iMesh_getEntSets(instance, pvRootSet, 0, &allSets, &allSets_alloc,
		   &allSets_size, &err);
  CHECK_ERR(err);

  for (i = 0; i < allSets_size; i++) {
    int iIsList;
    iBase_EntitySetHandle es = allSets[i];
    iMesh_isList(instance, es, &iIsList, &err);
    CHECK_ERR(err);

    int iNumAllEnts;
    iMesh_getNumOfType(instance, es, iBase_ALL_TYPES, &iNumAllEnts, &err);
    if (iNumAllEnts == 0) {
      apvSets[0] = es;
    }

    // It's possible to make a mistake with this simple test, but
    // another set would have to happen to hit the jackpot...
    if (iTotalEntsExpectedType == iNumAllEnts && !iIsList) {
      apvSets[3] = es;
    }
    else if (iTotalEntsExpectedTopo == iNumAllEnts && !iIsList) {
      apvSets[4] = es;
    }
    else if (iTotalEntsExpectedTopo < iNumAllEnts && iIsList) {
      apvSets[2] = es;
    }
    else if (iTotalEntsExpectedType < iNumAllEnts && iIsList) {
      apvSets[1] = es;
    }
  }
  bool qSetsOK = true;
  for (i = 0; i < 5; i++) {
    qSetsOK = qSetsOK && (apvSets[i] != NULL);
  }
  TEST(qSetsOK);
  UT_FREE(allSets);
}

void
vRecoverTags(iMesh_Instance instance, iBase_EntitySetHandle pvRootSet,
	     iBase_TagHandle aTags[], iBase_EntityHandle *aEnt)
{
  // This is easier than sets, because we can get the tags by name.
  int err;
  for (int i = 0; i < 4; i++) {
    iMesh_getTagHandle(instance, astrTagNames[i], &(aTags[i]), &err,
		       aiTagNameLen[i]);
    CHECK_ERR_CRITICAL(err);
  }

  if (aEnt) {
    // Identify the entity that's been tagged.
    iBase_EntityIterator iter;
    iMesh_initEntIter(instance, pvRootSet, iBase_ALL_TYPES,
		      iMesh_ALL_TOPOLOGIES, 0, &iter, &err);
    CHECK_ERR(err);
    int qGotData;
    char strTagName[100];
    iBase_EntityHandle ent;
    qGotData = true;
    while (qGotData) {
      iMesh_getNextEntIter(instance, iter, &ent, &qGotData, &err);
      CHECK_ERR(err);

      iBase_TagHandle *allTags = NULL;
      int allTags_alloc = 0, allTags_size;
      iMesh_getAllTags(instance, ent, &allTags, &allTags_alloc, &allTags_size,
		       &err);
      // Now loop over all of those tags.
      if (allTags_size < 4) {
	UT_FREE(allTags);
	continue;
      }
      int iNumFound = 0;
      for (int k = 0; k < 4; k++) {
	for (int j = 0; j < allTags_size; j++) {
	  if (aTags[k] == allTags[j]) {
	    iNumFound++;
	    break;
	  }
	}
      }
      UT_FREE(allTags);
      if (iNumFound == 4) {
	*aEnt = ent;
	break;
      }
    } // Done looping over entities.
    iMesh_endEntIter(instance, iter, &err);
  } // Done identifying the entity that has these tags on it.
}

void
vSetupTypeAndTopo()
{
  // Which combinations are okay doesn't depend on whether the mesh
  // database is 2D, 3D, etc.  Meshes are required to returned nothing
  // if they don't have something, that's all.

  int eEType;
  int eETopo;

  //@@ Setup data
  for (eETopo = iMesh_POINT; eETopo <= iMesh_ALL_TOPOLOGIES; eETopo++) {
    for (eEType = iBase_VERTEX; eEType < iBase_ALL_TYPES; eEType++) {
      a2qLegalTypeAndTopo[eEType][eETopo] = false;
    }
    a2qLegalTypeAndTopo[iBase_ALL_TYPES][eETopo] = true;
  }

  a2qLegalTypeAndTopo[iBase_VERTEX][iMesh_POINT] = true;

  a2qLegalTypeAndTopo[iBase_EDGE][iMesh_LINE_SEGMENT] = true;

  a2qLegalTypeAndTopo[iBase_FACE][iMesh_TRIANGLE] = true;
  a2qLegalTypeAndTopo[iBase_FACE][iMesh_QUADRILATERAL] = true;
  a2qLegalTypeAndTopo[iBase_FACE][iMesh_POLYGON] = true;

  a2qLegalTypeAndTopo[iBase_REGION][iMesh_TETRAHEDRON] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_PYRAMID] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_PRISM] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_HEXAHEDRON] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_SEPTAHEDRON] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_POLYHEDRON] = true;

  a2qLegalTypeAndTopo[iBase_VERTEX][iMesh_ALL_TOPOLOGIES] = true;
  a2qLegalTypeAndTopo[iBase_EDGE][iMesh_ALL_TOPOLOGIES] = true;
  a2qLegalTypeAndTopo[iBase_FACE][iMesh_ALL_TOPOLOGIES] = true;
  a2qLegalTypeAndTopo[iBase_REGION][iMesh_ALL_TOPOLOGIES] = true;
  a2qLegalTypeAndTopo[iBase_ALL_TYPES][iMesh_ALL_TOPOLOGIES] = true;
}

// Second adjacencies are tested by comparison to a chain of first
// adjacencies.
void
vSecondAdjacencyTest(iMesh_Instance instance)
{
  int err = iBase_SUCCESS;
  // Start by checking that first adjacencies don't return the entity itself.

  // Check all possible combinations of second adjacencies.
  for (int sourceType = iBase_VERTEX; sourceType < iBase_ALL_TYPES;
      sourceType++) {
    if (a2iAdjTable[sourceType][sourceType] != iBase_AVAILABLE)
      return;
    for (int bridgeType = iBase_VERTEX; bridgeType <= iBase_ALL_TYPES;
	bridgeType++) {
      for (int requestedType = iBase_VERTEX; requestedType <= iBase_ALL_TYPES;
	  requestedType++) {
	// Don't do this test if -either- of the required first
	// adjacencies is unavailable.
	if ((bridgeType != iBase_ALL_TYPES && requestedType != iBase_ALL_TYPES)
	    && ((a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_1
		&& a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_LOGN
		&& a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_N)
		|| (a2iAdjTable[bridgeType][requestedType] != iBase_ALL_ORDER_1
		    && a2iAdjTable[bridgeType][requestedType]
			!= iBase_ALL_ORDER_LOGN
		    && a2iAdjTable[bridgeType][requestedType]
			!= iBase_ALL_ORDER_N)))
	  continue;

	cout << "  Source type: ";
	cout.width(9);
	cout << astrTypes[sourceType] << "  Bridge type: ";
	cout.width(9);
	cout << astrTypes[bridgeType] << " Target type: ";
	cout.width(9);
	cout << astrTypes[requestedType] << endl;

	// Start an iterator and grab an entity.
	iBase_EntityIterator entIter;
	iMesh_initEntIter(instance, pvRootSet, sourceType, iMesh_ALL_TOPOLOGIES,
			  0, &entIter, &err);
	iBase_EntityHandle ent;
	int hasData, count = 0;
	bool qOK = true;
	iMesh_getNextEntIter(instance, entIter, &ent, &hasData, &err);
	while (hasData && qOK && count < 100) {
	  count++;
	  std::set<iBase_EntityHandle> setDirect;
	  {
	    // Grab second adjacency directly; the result goes into an
	    // STL set (setDirect) for later comparison.
	    iBase_EntityHandle *adjEnts = NULL;
	    int adjEnts_allocated = 0, adjEnts_size = 0;
	    iMesh_getEnt2ndAdj(instance, ent, bridgeType, requestedType,
			       &adjEnts, &adjEnts_allocated, &adjEnts_size,
			       &err);
	    for (int ii = 0; ii < adjEnts_size; ii++) {
	      setDirect.insert(adjEnts[ii]);
	    }
	    // Was the returned list unique?
	    qOK = (setDirect.size() == adjEnts_size);
	    UT_FREE(adjEnts);
	  }

	  // Spec says the original ent shouldn't be in the result.
	  qOK = (setDirect.find(ent) == setDirect.end());

	  // Some comparison cases are easy...
	  std::set<iBase_EntityHandle> setIndirect;
	  if ((bridgeType == requestedType || bridgeType == sourceType)
	      && bridgeType != iBase_ALL_TYPES) {
	    // The output should be empty, so do nothing.
	  }
	  else if ((sourceType < bridgeType && bridgeType < requestedType
	      && requestedType != iBase_ALL_TYPES)
	      || (sourceType > bridgeType && bridgeType > requestedType
		  && sourceType != iBase_ALL_TYPES)) {
	    // This case is simply a first adjacency call.   It might
	    // succeed even if the implementation doesn't support one of
	    // the requested adjacencies, but never gets tested in that
	    // case.
	    iBase_EntityHandle *adjEnts = NULL;
	    int adjEnts_allocated = 0, adjEnts_size = 0;
	    iMesh_getEntAdj(instance, ent, requestedType, &adjEnts,
			    &adjEnts_allocated, &adjEnts_size, &err);
	    for (int ii = 0; ii < adjEnts_size; ii++) {
	      setIndirect.insert(adjEnts[ii]);
	    }
	    UT_FREE(adjEnts);
	  }
	  else {
	    // Do two first adjacency calls, which should give the same
	    // result as the second adjacency call.

	    // Grab first adjacency...
	    iBase_EntityHandle *adjEnts = NULL;
	    int adjEnts_allocated = 0, adjEnts_size = 0;
	    iMesh_getEntAdj(instance, ent, bridgeType, &adjEnts,
			    &adjEnts_allocated, &adjEnts_size, &err);
	    CHECK_ERR(err);
	    // ... and first adjacencies of those
	    iBase_EntityHandle *adjEnts2 = NULL;
	    int adjEnts2_allocated = 0, adjEnts2_size;
	    int *offset, offset_allocated = 0, offset_size;
	    iMesh_getEntArrAdj(instance, adjEnts, adjEnts_size, requestedType,
			       &adjEnts2, &adjEnts2_allocated, &adjEnts2_size,
			       &offset, &offset_allocated, &offset_size, &err);
	    CHECK_ERR(err);
	    for (int ii = 0; ii < adjEnts2_size; ii++) {
	      setIndirect.insert(adjEnts2[ii]);
	    }
	    UT_FREE(adjEnts);
	    UT_FREE(offset);
	    UT_FREE(adjEnts2);
	  }
	  // Spec says the original ent shouldn't be in the result.
	  setIndirect.erase(ent);

	  // Compare results: STL set comparison should handle this
	  // flawlessly.
	  qOK = (setDirect == setIndirect);
// 	  if (!qOK) {
// 	    iBase_EntityHandle *aHandles =
// 	      new iBase_EntityHandle[setDirect.size()];
// 	    iBase_EntityHandle *result =
// 	      std::set_difference(setIndirect.begin(), setIndirect.end(),
// 				  setDirect.begin(), setDirect.end(), aHandles);
// 	    assert(result - aHandles > 0);
// 	  }
	  iMesh_getNextEntIter(instance, entIter, &ent, &hasData, &err);
	} // Loop over entities.
	iMesh_endEntIter(instance, entIter, &err);
	TEST(qOK);
      } // Loop over requested type
    } // Loop over bridge type
  } // Loop over source topologies
}

void
vSecondAdjacencyArrayTest(iMesh_Instance instance)
{
  int err = iBase_SUCCESS;
  // Start by checking that first adjacencies don't return the entity itself.

  // Check all possible combinations of second adjacencies.
  for (int sourceType = iBase_VERTEX; sourceType < iBase_ALL_TYPES;
      sourceType++) {
    for (int bridgeType = iBase_VERTEX; bridgeType <= iBase_ALL_TYPES;
	bridgeType++) {
      for (int requestedType = iBase_VERTEX; requestedType <= iBase_ALL_TYPES;
	  requestedType++) {
	if ((bridgeType != iBase_ALL_TYPES
	    && a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_1
	    && a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_LOGN
	    && a2iAdjTable[sourceType][bridgeType] != iBase_ALL_ORDER_N)
	    || (requestedType != iBase_ALL_TYPES
		&& a2iAdjTable[bridgeType][requestedType] != iBase_ALL_ORDER_1
		&& a2iAdjTable[bridgeType][requestedType]
		    != iBase_ALL_ORDER_LOGN
		&& a2iAdjTable[bridgeType][requestedType] != iBase_ALL_ORDER_N))
	  continue;

	cout << "  Source type: ";
	cout.width(9);
	cout << astrTypes[sourceType] << "  Bridge type: ";
	cout.width(9);
	cout << astrTypes[bridgeType] << " Target type: ";
	cout.width(9);
	cout << astrTypes[requestedType] << endl;

	// Start an iterator and grab an entity.
	iBase_EntityArrIterator wsIter;
	iMesh_initEntArrIter(instance, pvRootSet, sourceType,
			     iMesh_ALL_TOPOLOGIES, 100, 0, &wsIter, &err);
	iBase_EntityHandle *ents = NULL;
	int ents_allocated = 0, ents_size = 0;
	int hasData, count = 0;
	bool qOK = true;
	iMesh_getNextEntArrIter(instance, wsIter, &ents, &ents_allocated,
				&ents_size, &hasData, &err);
	while (hasData && qOK && count < 5) {
	  count++;

	  // Grab second adjacency for all entities directly
	  iBase_EntityHandle *adjEnts = NULL;
	  int adjEnts_allocated = 0, adjEnts_size = 0;
	  int *offsets = NULL, offsets_allocated = 0, offsets_size = 0;
	  iMesh_getEntArr2ndAdj(instance, ents, ents_size, bridgeType,
				requestedType, &adjEnts, &adjEnts_allocated,
				&adjEnts_size, &offsets, &offsets_allocated,
				&offsets_size, &err);

	  for (int i = 0; i < ents_size; i++) {
	    // Shove the array data for this ent into a set.
	    std::set<iBase_EntityHandle> setArray;
	    for (int ii = offsets[i]; ii < offsets[i + 1]; ii++) {
	      setArray.insert(adjEnts[ii]);
	    }

	    // Was the returned list unique?
	    qOK = (setArray.size() == offsets[i + 1] - offsets[i]);

	    // Now grab 2nd adj entity by entity
	    iBase_EntityHandle *adjEnts2 = NULL;
	    int adjEnts2_allocated = 0, adjEnts2_size = 0;
	    iMesh_getEnt2ndAdj(instance, ents[i], bridgeType, requestedType,
			       &adjEnts2, &adjEnts2_allocated, &adjEnts2_size,
			       &err);

	    // Shove that into a set as well.
	    std::set<iBase_EntityHandle> setSingle;
	    setSingle.insert(adjEnts2, adjEnts2 + adjEnts2_size);
	    UT_FREE(adjEnts2);
	    qOK = (setArray == setSingle);
	  }
	  UT_FREE(adjEnts);
	  UT_FREE(offsets);
	  iMesh_getNextEntArrIter(instance, wsIter, &ents, &ents_allocated,
				  &ents_size, &hasData, &err);
	} // Loop over entities.
	iMesh_endEntArrIter(instance, wsIter, &err);
	UT_FREE(ents);
	TEST(qOK);
      } // Loop over requested type
    } // Loop over bridge type
  } // Loop over source topologies
}

//@ Adjacency reciprocity test kernel

void
vEntityReciprocityTest(iMesh_Instance instance, iBase_EntityType eET1,
		       iBase_EntityType eET2)
{
  // Iterate over entities of type eET1
  int err;
  int qHasData, qOK = true;
  iBase_EntityIterator entIter;
  iMesh_initEntIter(instance, pvRootSet, eET1, iMesh_ALL_TOPOLOGIES, 0,
		    &entIter, &err);
  CHECK_ERR_CRITICAL(err);
  iBase_EntityHandle entity;
  iMesh_getNextEntIter(instance, entIter, &entity, &qHasData, &err);
  if (err != iBase_SUCCESS) {
    qOK = false;
  }
  while (qOK && qHasData) {
    // Grab the entities of type eET2 adjacent to entity.
    iBase_EntityHandle* aAdjEnt1 = NULL;
    int iNAdjEnt1, iAdjEnt_allocated = 0;
    iMesh_getEntAdj(instance, entity, eET2, &aAdjEnt1, &iAdjEnt_allocated,
		    &iNAdjEnt1, &err);
    if (err != iBase_SUCCESS) {
      qOK = false;
      break;
    }
    // Loop over those entities
    for (int iAdj1 = 0; iAdj1 < iNAdjEnt1 && qOK; iAdj1++) {
      iBase_EntityHandle adj_entity = aAdjEnt1[iAdj1];
      // Now get their adjacencies.
      iBase_EntityHandle* aAdjEnt2;
      int iNAdjEnt2, iAdjEnt2Alloc = 0;
      iMesh_getEntAdj(instance, adj_entity, eET1, &aAdjEnt2, &iAdjEnt2Alloc,
		      &iNAdjEnt2, &err);
      if (err != iBase_SUCCESS) {
	qOK = false;
	break;
      }
      // Loop over -those- to make sure the original entity is one of
      // them.
      qOK = false;
      for (int iAdj2 = 0; iAdj2 < iNAdjEnt2; iAdj2++) {
	iBase_EntityHandle adj_ent2 = aAdjEnt2[iAdj2];
	if (adj_ent2 == entity) {
	  qOK = true;
	  break;
	}
      } // Loop over entities of type eET1 adj to eET2's adj to the
	// original ent.
      UT_FREE(aAdjEnt2);
    } // Loop over entities of type eET2 adj to original ent
    UT_FREE(aAdjEnt1);
    iMesh_getNextEntIter(instance, entIter, &entity, &qHasData, &err);
    if (err != iBase_SUCCESS) {
      qOK = false;
    }
  } // Loop over entities of type eET1.
  iMesh_endEntIter(instance, entIter, &err);
  CHECK_ERR_CRITICAL(err);
  if (qOK) {
    cout << " passed";
    TEST(true);
  }
  else {
    cout << " failed";
    TEST(false);
  }
} // Done with entity reciprocity test

void
vArrayReciprocityTest(iMesh_Instance instance, iBase_EntityType eET1,
		      iBase_EntityType eET2)
{
  int err = 0;
  // First, get all eET1 -> eET2 adjacency info.
  iBase_EntityHandle *aHandles1 = NULL, *aHandles2 = NULL, *aAdjHandles1 =
  NULL;
  int *aiFirstOffset = NULL, *aiSecondOffset = NULL;
  int iHandles1_allocated = 0, iHandles2_allocated = 0;
  int iAdjHandles_allocated = 0, iOffset1_allocated = 0;
  int iOffset2_allocated = 0;
  int iNumEntities;
  iMesh_getEntities(instance, pvRootSet, eET1, iMesh_ALL_TOPOLOGIES, &aHandles1,
		    &iHandles1_allocated, &iNumEntities, &err);
  CHECK_ERR_CRITICAL(err);
  int iNum1st, iOffset1st, iNum2nd, iOffset2nd;
  iMesh_getEntArrAdj(instance, aHandles1, iNumEntities, eET2, &aHandles2,
		     &iHandles2_allocated, &iNum1st, &aiFirstOffset,
		     &iOffset1_allocated, &iOffset1st, &err);
  CHECK_ERR_CRITICAL(err);
  iMesh_getEntArrAdj(instance, aHandles2, iNum1st, eET1, &aAdjHandles1,
		     &iAdjHandles_allocated, &iNum2nd, &aiSecondOffset,
		     &iOffset2_allocated, &iOffset2nd, &err);
  CHECK_ERR_CRITICAL(err);

  bool qOK = (iNumEntities == iOffset1st - 1) && (iNum1st == iOffset2nd - 1)
      && (iNum1st == aiFirstOffset[iNumEntities])
      && (iNum2nd == aiSecondOffset[iNum1st]) && (aiFirstOffset[0] == 0)
      && (aiSecondOffset[0] == 0);

  // For every entity, check that every adjacent entity also thinks
  // the first entity is adjacent to it.
  int iEH1, iEH2;
  for (iEH1 = 0; qOK && iEH1 < iNumEntities; iEH1++) {
    void * Handle1 = aHandles1[iEH1];
    int iBegin = aiFirstOffset[iEH1];
    int iEnd = aiFirstOffset[iEH1 + 1] - 1;
    for (iEH2 = iBegin; iEH2 <= iEnd; iEH2++) {
      bool qFoundIt = false;
      int iBegin2 = aiSecondOffset[iEH2];
      int iEnd2 = aiSecondOffset[iEH2 + 1] - 1;
      for (int iAEH = iBegin2; !qFoundIt && iAEH <= iEnd2; iAEH++) {
	iBase_EntityHandle AdjHandle = aAdjHandles1[iAEH];
	if (Handle1 == AdjHandle)
	  qFoundIt = true;
      }
      qOK = qFoundIt;
    } // Done with inner loop
  } // Done with outer loop
  UT_FREE(aHandles1);
  UT_FREE(aHandles2);
  UT_FREE(aAdjHandles1);
  UT_FREE(aiFirstOffset);
  UT_FREE(aiSecondOffset);
  if (qOK) {
    cout << " passed";
    TEST(true);
  }
  else {
    cout << " failed";
    TEST(false);
  }
} // Done with array reciprocity test

void
vTestArrModInterface(iMesh_Instance& instance)
{
  int err;
  // setVtxArrCoords is trivial to test:  get vertices, change their coords,
  // and verify that the change worked.
  iBase_EntityHandle *aVerts = NULL;
  int iNVerts;
  int iVertsAlloc = 0;
  iMesh_getEntities(instance, pvRootSet, iBase_VERTEX, iMesh_POINT, &aVerts,
		    &iVertsAlloc, &iNVerts, &err);
  CHECK_ERR(err);

  double * adCoords1 = NULL;
  int iCoords1Size;
  int SO = iBase_INTERLEAVED;
  int iCoords1Alloc = 0;
  iMesh_getVtxArrCoords(instance, aVerts, iNVerts, SO, &adCoords1,
			&iCoords1Alloc, &iCoords1Size, &err);
  CHECK_ERR(err);
  int iGeomDim;
  iMesh_getGeometricDimension(instance, &iGeomDim, &err);
  CHECK_ERR(err);

  bool qOK;
  {
    cout << " Setting vertex coordinates..." << endl;
    double * adCoordsTmp = (double*) calloc(iNVerts * iGeomDim, sizeof(double));
    for (int iV = 0; iV < iNVerts; iV++) {
      int iOff = iGeomDim * iV;
      for (int i = 0; i < iGeomDim; i++) {
	adCoordsTmp[iOff + i] = adCoords1[iOff + i] + iV + i + 1;
      }
    }

    iMesh_setVtxArrCoords(instance, aVerts, iNVerts, SO, adCoordsTmp,
			  iCoords1Size, &err);
    CHECK_ERR(err);

    double * adCoords2 = NULL;
    int iCoords2Size;
    int iCoords2Alloc = 0;
    iMesh_getVtxArrCoords(instance, aVerts, iNVerts, SO, &adCoords2,
			  &iCoords2Alloc, &iCoords2Size, &err);
    CHECK_ERR(err);

    qOK = (iCoords2Size == iCoords1Size);
    for (int iV = 0; iV < iNVerts; iV++) {
      int iOff = iGeomDim * iV;
      for (int i = 0; i < iGeomDim; i++) {
	qOK = qOK
	    && (APPROX_EQUAL(adCoords1[(iOff + i)],
			     adCoords2[(iOff + i)] - iV - i - 1));
      }
    }
    TEST(qOK);
    // Now set them back.
    iMesh_setVtxArrCoords(instance, aVerts, iNVerts, SO, adCoords1,
			  iCoords1Size, &err);
    CHECK_ERR(err);
    UT_FREE(adCoordsTmp);
    UT_FREE(aVerts);
    UT_FREE(adCoords1);
    UT_FREE(adCoords2);
  }

  // Now create a flock of new vertices.  Then verify that their handles
  // and coordinates are returned properly .
  {
    qOK = true;
    cout << " Creating new vertices..." << endl;
    double * adCoords2 = (double*) calloc(10 * iGeomDim, sizeof(double));
    for (int i = 0; i < 10; i++) {
      switch (iGeomDim)
	{
	case 3:
	  adCoords2[i * iGeomDim + 0] = 3.14159 * i; // i*pi
	  adCoords2[i * iGeomDim + 1] = 2.81828 * i; // i*e
	  adCoords2[i * iGeomDim + 2] = 1.414 * i; // i*sqrt(2)
	  break;
	case 2:
	  adCoords2[i * iGeomDim + 0] = 3.14159 * i; // i*pi
	  adCoords2[i * iGeomDim + 1] = 2.81828 * i; // i*e
	  break;
	default:
	  assert(0);
	  break;
	}
    }
    vDoSaveLoad(&instance, &pvRootSet);

    iBase_EntityHandle *aNewVerts = NULL;
    int iNumNewVerts;
    int SO = iBase_INTERLEAVED;
    int iNewVertsAlloc = 0;
    iMesh_createVtxArr(instance, 10, iBase_INTERLEAVED, adCoords2,
		       10 * iGeomDim, &aNewVerts, &iNewVertsAlloc,
		       &iNumNewVerts, &err);
    CHECK_ERR(err);
    qOK = (iNumNewVerts == 10);

    // Now verify the coords.
    double * adCoords3 = NULL;
    int iCoords3Size;
    int iCoords3Alloc = 0;
    iMesh_getVtxArrCoords(instance, aNewVerts, iNumNewVerts, SO, &adCoords3,
			  &iCoords3Alloc, &iCoords3Size, &err);
    CHECK_ERR(err);
    qOK = qOK && (iCoords3Size == iGeomDim * 10);

    for (int iV = 0; iV < 10; iV++) {
      for (int i = 0; i < iGeomDim; i++) {
	qOK =
	    qOK
		&& (adCoords3[(iV * iGeomDim + i)]
		    == adCoords2[(iV * iGeomDim + i)]);
      }
    }
    UT_FREE(adCoords2);
    UT_FREE(adCoords3);
    TEST(qOK);
    iMesh_deleteEntArr(instance, aNewVerts, iNumNewVerts, &err);
    CHECK_ERR(err);
    vDoSaveLoad(&instance, &pvRootSet);
    UT_FREE(aNewVerts);

    // Check that they were really deleted.
    int iNewNVerts;
    iMesh_getNumOfTopo(instance, pvRootSet, iMesh_POINT, &iNewNVerts, &err);
    CHECK_ERR(err);
    TESTEQ(iNewNVerts, iNVerts);
  } // Done checking array vertex creation
}

void
vTestTagInterface(iMesh_Instance& instance)
{
  //     createTag, destroyTag, getTagName, getTagSizeValues,
  //     getTagSizeBytes, getTagHandle, getTagType
  iBase_TagHandle aTags[4], pvBadTag, theTag;
  char strTagName[100];
  int err;
  int iType, iSizeBytes, iSizeValues;

  cout << "Creating tags of each type...                             ";

  bool qNamesOK = true, qTypesOK = true, qSizeBytesOK = true, qSizeValuesOK =
      true, qRetrieved = true;
  // Create a tag of each type and check that the queries for each of the various
  // parameters are working properly
  for (int i = 0; i < 4; i++) {
    iMesh_createTag(instance, astrTagNames[i], aiTagSizes[i], aTagTypes[i],
		    &aTags[i], &err, aiTagNameLen[i]);
    CHECK_ERR(err);

    iMesh_getTagName(instance, aTags[i], strTagName, &err, 100);
    CHECK_ERR(err);
    qNamesOK = qNamesOK && (!strncmp(astrTagNames[i], strTagName, 100));

    iMesh_getTagType(instance, aTags[i], &iType, &err);
    CHECK_ERR(err);
    qTypesOK = qTypesOK && (iType == aTagTypes[i]);

    iMesh_getTagSizeBytes(instance, aTags[i], &iSizeBytes, &err);
    CHECK_ERR(err);
    qSizeBytesOK = qSizeBytesOK
	&& (iSizeBytes == aiTagBytes[i] * aiTagSizes[i]);

    iMesh_getTagSizeValues(instance, aTags[i], &iSizeValues, &err);
    CHECK_ERR(err);
    qSizeValuesOK = qSizeValuesOK && (iSizeValues == aiTagSizes[i]);

    iMesh_getTagHandle(instance, astrTagNames[i], &theTag, &err,
		       aiTagNameLen[i]);
    CHECK_ERR(err);
    qRetrieved = qRetrieved && (theTag == aTags[i]);
  }

  TEST(qNamesOK);
  TEST(qTypesOK);
  TEST(qSizeBytesOK);
  TEST(qSizeValuesOK);
  TEST(qRetrieved);
  cout << "DONE\n";

  if (qErrTests) {
    cout << "Testing error throwing for Tag methods...                 ";
    // Error throwing tests simply try to create a situation where the interface
    // should, according to the spec, throw an error, and verify that:
    // (a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reversal of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice versa, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.
    // Duplicate tag name not caught
    iMesh_createTag(instance, astrTagNames[0], 2, iBase_DOUBLE, &pvBadTag, &err,
		    aiTagNameLen[0]);
    CHECK_ERR2(err, iBase_TAG_ALREADY_EXISTS);

    // Invalid tag value size not caught
    iMesh_createTag(instance, "Bad", -1, iBase_ENTITY_HANDLE, &pvBadTag, &err,
		    3);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    // Invalid tag value type not caught
    iMesh_createTag(instance, "Bad", 1, iBase_TagValueType(5), &pvBadTag, &err,
		    3);
    CHECK_ERR2(err, iBase_INVALID_ARGUMENT);

    //
    cout << "DONE\n";
  }

  cout << "Testing tag destruction...                                ";
  // Try to destroy tags... can't test forced/unforced destruction just yet
  for (int i = 0; i < 4; i++) {
    iMesh_destroyTag(instance, aTags[i], false, &err);
    CHECK_ERR(err);
  }
  cout << "DONE\n";
}

// Test some of the version macros
static void
vTestVersionMacros()
{
  const bool have_version_ge = true;
  const bool version_ge_0_0_8 = true;
  const bool version_ge_99_3_5 = false;
  bool ge;

// It'd be nicer to use #if defined(IMESH_VERSION_GE) && IMESH_VERSION_GE(1,3,0)
// but I've seen CPP that don't handle the && operator correctly.
#ifdef IMESH_VERSION_GE

#if IMESH_VERSION_GE(0,0,8)
  ge = true;
#else
  ge = false;
#endif // #if IMESH_VERSION_GE(0,0,8)
  TESTEQ(ge, version_ge_0_0_8);

#if IMESH_VERSION_GE(99,3,5)
  ge = true;
#else
  ge = false;
#endif // #if IMESH_VERSION_GE(99,3,5)
  TESTEQ(ge, version_ge_99_3_5);

#else

  TESTEQ(false, have_version_ge);

#endif // #ifdef IMESH_VERSION_GE

  char tmp[256];
  sprintf(tmp, "iMesh_Version_%1d.%1d.%1d",
  IMESH_VERSION_MAJOR,
	  IMESH_VERSION_MINOR, IMESH_VERSION_PATCH);
  TESTEQ(IMESH_VERSION_STRING, string(tmp));
}

//@ Test iterators in the presence of deletion and addition of entities.
// Basically, getting iterators to always behave properly when entities
// are added or deleted during iteration is a non-trivial exercise.  The
// current test checks only iterators for edges and triangles (for 2D
// meshes), or triangles and tetrahedra (for 3D meshes), as those are
// what's needed for the swapping service.  Other iterators should be
// bomb-proofed in the same way, and should be checked eve

//@ Main program

int
main(int argc, char *argv[])
{
  // Check command line arg
  int iFile = 1, err;
  string sArg;
  string filename = "";
  const char* sOutFile;
  struct sigaction sa_struct;
  while (iFile < argc) {
    sArg = argv[iFile];
    if (sArg == "-h" || sArg == "--help") {
      TYPE_HELP
      ;
      TYPE_USAGE
      ;
      TYPE_USAGE_HELP
      ;
      return 0;
    }
    else if (sArg == "--no-error") {
      qErrTests = false;
    }
    else if (sArg == "-s" || sArg == "--silent") {
      qSilent = true;
    }
    else if (argv[iFile][0] == '-' && argv[iFile][1] == 'o') {
      if (argv[iFile][2] == '\0') {
	sOutFile = "out_test.txt";
      }
      else {
	sOutFile = (argv[iFile] + 2);
      }
      qOutFile = true;
    }
    else if (argv[iFile][0] == '-' && argv[iFile][1] == '-'
	&& argv[iFile][2] == 'o' && argv[iFile][3] == 'u'
	&& argv[iFile][4] == 't' && argv[iFile][5] == 'p'
	&& argv[iFile][6] == 'u' && argv[iFile][7] == 't') {
      if (argv[iFile][8] != '\0' && argv[iFile][8] != '=') {
	sArg = (argv[iFile] + 2);
	TYPE_USAGE
	;
	std::cerr << "Invalid option: '" << sArg << "'.  Try '" << argv[0]
	    << " --help'." << endl << endl;
	return 1;
      }
      else if (argv[iFile][8] == '=') {
	if (argv[iFile][9] == '\0') {
	  TYPE_USAGE
	  ;
	  std::cerr << "No output file given!  Try '" << argv[0] << " --help'."
	      << endl << endl;
	  return 1;
	}
	sOutFile = (argv[iFile] + 9);
      }
      else {
	sOutFile = "out_test.txt";
      }
      qOutFile = true;
    }
    else if (sArg == "--no-save-load") {
      qSaveLoad = false;
      qSaveLoadClean = false;
    }
    else if (sArg == "--no-save-load-clean") {
      qSaveLoadClean = false;
    }
    else if (argv[iFile][0] == '-') {
      if (argv[iFile][1] == '-')
	sArg = (argv[iFile] + 2);
      else
	sArg = (argv[iFile] + 1);
      TYPE_USAGE
      ;
      std::cerr << "Invalid option: '" << sArg << "'.  Try '" << argv[0]
	  << " --help'." << endl << endl;
      return 1;
    }
    else {
      filename = argv[iFile];
    }
    iFile++;
  }
  if (filename == "") {
    TYPE_USAGE
    ;
    std::cerr << "Need an input file!  Try '" << argv[0] << " --help'." << endl
	<< endl;
    return 1;
  }

  // This signal stuff is here solely to affect cleanup
  // of tmp dir and files on irregular exit for only the
  // most likely signals to encounter in normal development.
  memset(&sa_struct, 0, sizeof(sa_struct));
  sa_struct.sa_handler = CleanSaveLoadOnSignal;
  sigaction(SIGFPE, &sa_struct, 0);
  sigaction(SIGABRT, &sa_struct, 0);
  sigaction(SIGSEGV, &sa_struct, 0);
  sigaction(SIGTRAP, &sa_struct, 0);
  sigaction(SIGBUS, &sa_struct, 0);

#ifdef SETVBUF_REVERSED
  setvbuf(stdout, _IONBF, NULL, 1);
#else
  setvbuf(stdout, NULL, _IONBF, 1);
#endif

  std::streambuf* outbuf = cout.rdbuf();
  std::streambuf* errbuf = std::cerr.rdbuf();
  std::ofstream file, nullFile;
  if (qOutFile) {
    file.open(sOutFile);
    if (!file.is_open()) {
      std::cerr << endl << "Error: Could not open output file '" << sOutFile
	  << "'!" << endl << endl;
      return 1;
    }
    cout.rdbuf(file.rdbuf());
    std::cerr.rdbuf(file.rdbuf());
  }
  if (qSilent) {
    nullFile.open("/dev/null/");
    cout.rdbuf(nullFile.rdbuf());
    std::cerr.rdbuf(nullFile.rdbuf());
  }
  if (qSaveLoad) {
    snprintf(TMPDIR, sizeof(TMPDIR), "./tmpdata_%llu",
	     (unsigned long long) getpid());
    mkdir(TMPDIR, 0770);
  }

  vTestVersionMacros();

  iMesh_Instance instance;
  iMesh_newMesh("", &instance, &err, 0);
  if (err != iBase_SUCCESS) {
    cout << "Can't create an iMesh instance; giving up!" << endl;
    exit(1);
  }

  vSetupTypeAndTopo();

  {
    // Test stuff specific to interface iMesh_  Those functions are:
    //    load, save, getRootSet, getGeometricDimension, getDfltStorage,
    //    getAdjTable, getNumOfType, getNumOfTopo, getAllVtxCoords,
    //    getVtxCoordIndex, getEntities, getVtxArrCoords, and
    //    getAdjEntities

    cout << "Testing mesh interface..." << endl;
    qAllPassed = true;
    vTestMeshInterface(instance, filename);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface Entity.  Those functions are:
  //    initEntIter, getNextEntIter, resetEntIter, endEntIter,
  //    getEntTopo, getEntType, getVtxCoord, getEntAdj
  bool qEntityPassed = false;
  {
    qAllPassed = true;
    vTestEntityInterface(instance);
    if (qAllPassed) {
      qEntityPassed = true;
      iNIntPassed++;
    }
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface Arr.  Those functions are:
  //     initEntArrIter, getEntArrNextIter, resetEntArrIter,
  //     endEntArrIter, getEntArrTopo, getEntArrType, getEntArrAdj
  {
    qAllPassed = true;
    vTestArrInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface Modify.  Those functions are:
  //     setVtxCoords, createVtx, createEnt, deleteEnt
  bool qModifyPassed = false;
  {
    qAllPassed = true;
    vTestModifyInterface(instance);
    if (qAllPassed) {
      qModifyPassed = true;
      iNIntPassed++;
    }
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface ModArr.  Those functions are:
  //     setVtxArrCoords, createVtxArr, createEntArr, deleteEntArr
  {
    qAllPassed = true;
    vTestArrModInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Tag and entity set testing, appended June 2005 by Ben Kerby

  // Test stuff specific to interface EntSet.  Those functions are:
  //     createEntSet, destroyEntSet, isList, getNumEntSets, getEntSets,
  //     addEntToSet, rmvEntFromSet, addEntArrToSet, rmvEntArrFromSet,
  //     addEntSet, rmvEntSet, isContainedIn
  {
    cout << "Testing entity set interface..." << endl;
    qAllPassed = true;
    vTestEntSetInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface SetRelation.  Those functions are:
  //     addPrntChld, rmvPrntChld, isChildOf, getNumChld, getNumPrnt,
  //     getChldn, getPrnts
  {
    qAllPassed = true;
    vTestSetRelationInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface SetBoolOps.  Those functions are:
  //     subtract, intersect, unite
  {
    qAllPassed = true;
    vTestBoolOpsInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Now that entity sets have been tested, go back and re-test query
  // using entity sets.

  // Retest stuff specific to interface Mesh that can take entity_set
  // args.  Those functions are:
  //    getNumOfType, getNumOfTopo, getAllVtxCoords,
  //    getVtxCoordIndex, getEntities, and getAdjEntities

  {
    qAllPassed = true;
    vTestMeshWithSets(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test single-entity iterators over EntitySets. Those functions are:
  //    initEntIter, getNextEntIter, resetEntIter, endEntIter
  {
    cout << "Testing single entity iterators with entity sets..." << endl;
    qAllPassed = true;
    vTestEntityWithSets(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test array iterators over EntitySets. Those functions are:
  //     initEntArrIter, getEntArrNextIter, resetEntArrIter,
  //     endEntArrIter
  {
    cout << "Testing entity array iterators over entity sets..." << endl;
    qAllPassed = true;
    vTestArrayWithSets(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface Tag.  Those functions are:
  //     createTag, destroyTag, getTagName, getTagSizeValues,
  //     getTagSizeBytes, getTagHandle, getTagType
  {
    cout << "Testing tag interface..." << endl;
    qAllPassed = true;
    vTestTagInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface EntTag.  Those functions are:
  //     getData, getIntData, getDblData, getEHData, setData,
  //     setIntData, setDblData, setEHData, getAllTags, rmvTag
  {
    cout << "Testing entity tag interface..." << endl;
    qAllPassed = true;
    vTestEntTagInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface ArrTag.  Those functions are:
  //     getArrData, getIntArrData, getDblArrData, getEHArrData,
  //     setArrData, setIntArrData, setDblArrData, setEHArrData,
  //     rmvArrTag
  {
    cout << "Testing array tag interface..." << endl;
    qAllPassed = true;
    vTestArrTagInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test stuff specific to interface SetTag.  Those functions are:
  //     getEntSetData, getEntSetIntData, getEntSetDblData,
  //     getEntSetEHData, setEntSetData, setEntSetIntData,
  //     setEntSetDblData, setEntSetEHData, getAllEntSetTags,
  //     rmvEntSetTag
  {
    cout << "Testing entity set tag interface..." << endl;
    qAllPassed = true;
    vTestSetTagInterface(instance);
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // Test coordinate handling for 1 and 2D meshes
  {
    cout << "Testing lower geometric dimensioned meshes..." << endl;
    qAllPassed = true;
    vTestLowerDimensionedCoordinates();
    if (qAllPassed)
      iNIntPassed++;
    else
      iNIntFailed++;
    TEST_STATS
    ;
  }

  // These tests aren't yet ready for primetime, so they're commented out.

  // Now test the behavior of iterators in the presence of mesh
  // modification, provided that both Entity and Modify interfaces
  // passed their tests.
//   if (qEntityPassed && qModifyPassed) {
//     qAllPassed = true;
//     vTestIterWithModify(iMesh_getGeometricDimension(instance, ));
//     if ( qAllPassed )
//       iNIntPassed++;
//     else
//       iNIntFailed++;
//     TEST_STATS;
//   }
//   else {
//     iNNotImpl++;
//     cout << "Can't test iterators with Modify; either Entity or Modify failed a test.\n" << endl;
//   }

  if (qSilent || qOutFile) {
    cout.rdbuf(outbuf);
    std::cerr.rdbuf(errbuf);
  }

  std::stringstream stst;

  stst << endl << "FINAL RESULTS: " << std::setw(4) << iNPassed
      << " tests passed, " << endl << "               " << std::setw(4)
      << iNFailed << " tests failed in a total of " << endl << "               "
      << std::setw(4) << iNPassed + iNFailed << " tests, a pass rate of "
      << (int) 100 * iNPassed / (iNPassed + iNFailed) << "%." << endl << endl
      << "               " << std::setw(4)
      << iNIntPassed + iNIntFailed + iNNotImpl
      << " interfaces/interface combinations were tested: " << endl
      << "               " << std::setw(4) << iNIntPassed
      << " passed all interface-specific tests"
      << (iNNotImpl == 0 ? " and" : ", ") << endl << "               "
      << std::setw(4) << iNIntFailed << " failed at least one test";
  if (iNNotImpl != 0) {
    stst << ", and " << endl << "               " << std::setw(4) << iNNotImpl
	<< (iNNotImpl == 1 ? " was " : " were ") << "not implemented. " << endl
	<< endl;
  }
  else
    stst << "." << endl << endl;

  cout << stst.str();

  if (qOutFile) {
    file << stst.str();
    file.close();
  }
  if (qSaveLoadClean) {
    CleanSaveLoad();
  }

  iMesh_dtor(instance, &err);
  exit(iNFailed);
}
