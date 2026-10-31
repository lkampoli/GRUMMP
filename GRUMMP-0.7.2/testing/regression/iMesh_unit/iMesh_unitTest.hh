//#include "iMesh_unitTest_config.h"

#include "iBase.h"
#include "iMesh.h"

//-------------------------------------------------------

//@@ Standard include files

#include <algorithm>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

//@@ Macros and file static variables.

#define TYPE_HELP do {				\
    cout << endl;							\
    cout << "This file is intended to test functionality of a TSTT implementation." << endl; \
    cout << "Mesh query, modification, iterators, entity sets, and tags are all" << endl; \
    cout << "tested.  In cases where it's practical to define a priori a correct" << endl; \
    cout << "result for a query, that is done.  In other cases, consistency is the" << endl; \
    cout << "only thing checked.  In still others, the only possible correctness" << endl; \
    cout << "check is that the function took correct input and didn't choke" << endl; \
    cout << "(things like getAdjTable or getGeometricDim, for instance)." << endl; \
    cout << endl;							\
    cout << "It is worth noting that, although efforts have been made to make each" << endl; \
    cout << "test as independent as possible, it is unavoidable that a failure" << endl; \
    cout << "during one test may cause failures in subsequent tests, sometimes" << endl; \
    cout << "even across different interface tests.  In light of this, it is best" << endl; \
    cout << "to try to correct errors in the order they occur, or at least to" << endl;	\
    cout << "exclude interface tests that cause errors to see if later errors are" << endl; \
    cout << "avoided." << endl;						\
    cout << endl;							\
    cout << "Correct argument handling, including both pre- and post-conditions," << endl; \
    cout << "is checked elsewhere (iMesh_precondTest.cc and" << endl;	\
    cout << "iMesh_postcondTest.cc)." << endl;				\
  } while (0)

#define TYPE_USAGE do {\
    cerr << endl;							\
    cerr << "Usage: " << argv[0] << " [ -h | --help ]" << endl;	\
    cerr << "   or  " << argv[0] << " [arguments] <mesh_filename>" \
	      << endl << endl;						\
    cerr << "where help and arguments are:" << endl;		\
    cerr << " -h       | --help            Display usage information and exit." << endl; \
    cerr << " -s       | --silent          Display only final result." << endl; \
    cerr << " -o<file> | --output=<file>   Output the results to the filename given." << endl; \
    cerr << " --no-error                   Don't test all error throwing." << endl; \
    cerr << " --no-save-load               Don't attempt periodic saves and re-loads." << endl; \
    cerr << " --no-save-load-clean         Don't clean up temporary files from periodic saves " \
         << "and re-loads." << endl << endl; \
  } while (0)

#define TYPE_USAGE_HELP do {						\
    cerr << "All arguments must be given before the mesh filename, where:" << endl; \
    cerr << " 'help'    displays this message and ignores all other options," << endl; \
    cerr << " 'error'   adds more tests to more thoroughly test TSTT methods" << endl; \
    cerr << "           with regard to error-throwing circumstances," << endl; \
    cerr << " 'silent'  limits the output to only the final error count, and" << endl; \
    cerr << " 'output'  sends the program output to the filename given in the" << endl; \
    cerr << "           rest of the argument (the default is 'out_test.txt')." << endl << endl; \
  } while (0)

#define TEST(a) do {							\
    if (a) iNPassed++;							\
    else {								\
      qAllPassed = false;						\
      qFailed();							\
      cout << endl << __FILE__ << ":" << __LINE__ << ": " "Failure testing " << #a << endl; \
    }									\
  } while (0)

#define TEST2(a, b) do {							\
    if (a) iNPassed++;							\
    else {								\
      qAllPassed = false;						\
      iNFailed++;							\
      cout << endl << __FILE__ << ":" << __LINE__ << ": " "Failure: " << b << endl; \
    }									\
  } while (0)

#define TESTEQ(a,b) do {                                                        \
    if (a == b) iNPassed++;                                                     \
    else {                                                              \
      qAllPassed = false;                                               \
      qFailed();                                                       \
      cout << endl << __FILE__ << ":" << __LINE__ << ": Failure testing " \
           << #a << "(" << a << ") == " << #b << "(" <<b << ")" << endl; \
    }                                                                   \
  } while (0)

#define TEST_STATS do {							\
    cout << endl << "Test results so far:  " << iNPassed << " passed, " \
	 << iNFailed << " failed." << endl;				\
  } while(0)

// The following is currently not counted toward the pass total, as a
// successful return doesn't mean anything --- an empty function body
// can do that!
inline void
writeErrorDescr(int err)
{
  static const char* errorStrings[] =
    { "iBase_SUCCESS", "iBase_MESH_ALREADY_LOADED", "iBase_FILE_NOT_FOUND",
	"iBase_FILE_WRITE_ERROR", "iBase_NIL_ARRAY", "iBase_BAD_ARRAY_SIZE",
	"iBase_BAD_ARRAY_DIMENSION", "iBase_INVALID_ENTITY_HANDLE",
	"iBase_INVALID_ENTITY_COUNT", "iBase_INVALID_ENTITY_TYPE",
	"iBase_INVALID_ENTITY_TOPOLOGY", "iBase_BAD_TYPE_AND_TOPO",
	"iBase_ENTITY_CREATION_ERROR", "iBase_INVALID_TAG_HANDLE",
	"iBase_TAG_NOT_FOUND", "iBase_TAG_ALREADY_EXISTS", "iBase_TAG_IN_USE",
	"iBase_INVALID_ENTITYSET_HANDLE", "iBase_INVALID_ITERATOR_HANDLE",
	"iBase_INVALID_ARGUMENT", "iBase_MEMORY_ALLOCATION_FAILED",
	"iBase_NOT_SUPPORTED", "iBase_FAILURE" };
  cout << errorStrings[err] << endl;
}

#define CHECK_ERR_CRITICAL(err_) do {					\
    if ( err_ != iBase_SUCCESS) {					\
      writeErrorDescr(err_);						\
      qFailed();							\
      cerr << "Cannot continue interface test, exiting..." << endl;     \
      if ( qSaveLoadClean )                                             \
        CleanSaveLoad();                                                \
      exit(iNFailed);                                                   \
    }									\
    else								\
      iNPassed++;							\
  } while (0)

#define CHECK_ERR(err_) do {					\
    if (err_) {							\
      cout << "At line " << __LINE__ << " in " << __FILE__	\
	   << ", hit error:  ";					\
      writeErrorDescr(err_);					\
      qFailed();						\
    }								\
  } while (0)

#define CHECK_ERR2(err_, expected_err) do {			\
    if ( err_ != expected_err) {				\
      cout << "At line " << __LINE__ << " in " << __FILE__	\
	   << ", should have hit error:  ";			\
      writeErrorDescr(expected_err);				\
      cout << "Instead, hit error:  ";				\
      writeErrorDescr(err_);					\
      qFailed();						\
    }								\
    else							\
      iNPassed++;						\
  } while (0)

// The following macro is a little awkward: it can't be used in a
// single-line block (like in an if-else, for instance), because it's
// not a single line!
#define BREAK_ON_ERR(err_) 			\
    CHECK_ERR(err_);				\
    if (err_) break

//@ Testing for bad array arguments
#define TEST_BAD_ARRAY_SIZE CHECK_ERR2(err, iMesh_BAD_ARRAY_SIZE)
#define TEST_NIL_ARRAY CHECK_ERR2(err, iMesh_NIL_ARRAY)

#define GOT_ALL_ITEMS(_input, _output, _numItems, _qGotAll) do {	\
    _qGotAll = true;							\
    bool _gotItem;							\
    for ( int _i=0; _i<_numItems && _qGotAll; _i++ ) {			\
      _gotItem = false;							\
      for ( int _j=0; _j<_numItems && !_gotItem; _j++ )			\
	_gotItem = _gotItem || (_input[_i] == _output[_j]);		\
      _qGotAll = _qGotAll && _gotItem;					\
    }									\
  } while (0)

#define TOL 1.e-10
#define APPROX_EQUAL(a, b) (((a) < (b) + TOL) && ((b) < (a) + TOL))
#define UT_FREE(a) do {if (a) free(a); a = NULL;} while(0)

#define GET_ENTS(aEnts_, iTargetCount_) do {                            \
    iBase_EntityHandle *aEntsTmp = NULL;                                \
    int iEntsTmpAlloc = 0, iNEnts;                                      \
    iMesh_getEntities(instance, pvRootSet, iBase_ALL_TYPES,             \
                      iMesh_ALL_TOPOLOGIES,                             \
                      &aEntsTmp, &iEntsTmpAlloc, &iNEnts, &err);        \
    CHECK_ERR(err);                                                     \
    if (err || iNEnts < iTargetCount_) {                                \
      cout << "\nCould not get mesh data to work with... exiting...\n"; \
      TEST ( data_preparation );                                        \
      if ( qSaveLoadClean )                                             \
        CleanSaveLoad();                                                \
      exit(iNFailed);                                                   \
    }                                                                   \
    for (int iE = 0; iE < iTargetCount_; iE++) {                        \
      aEnts_[iE] = aEntsTmp[iE];                                        \
    }                                                                   \
    UT_FREE(aEntsTmp);                                                  \
  } while(0)

// Variables defined to make TEST statments more descriptive
static const bool error_throwing = false;
static const bool data_preparation = false;

#ifdef MAIN_PROGRAM
int iNPassed = 0, iNFailed = 0,
iNNotImpl = 0, iNIntPassed = 0, iNIntFailed = 0;
bool qAllPassed = true, qErrTests = true, qSilent = false, qOutFile = false;
bool qSaveLoad = true;
bool qSaveLoadClean = true;

int a2iAdjTable[4][4];

iBase_EntitySetHandle pvRootSet = NULL;
bool a2qLegalTypeAndTopo[5][12];
int a2iEntCount[5][12];
const char *astrTypes[] = {"vertex",
  "edge",
  "face",
  "region",
  "all types"};
const char *astrTopologies[] = {"point",
  "line segment",
  "polygon",
  "triangle",
  "quadrilateral",
  "polyhedron",
  "tetrahedron",
  "pyramid",
  "prism",
  "hexahedron",
  "septahedron",
  "all topologies"};

const char *astrTagNames[] = {"unitTest Int Tag", "unitTest Double Tag",
  "unitTest Ent Tag", "unitTest Generic Tag"};
int aiTagNameLen[] = {16, 19, 16, 20};
int aTagTypes[] = {iBase_INTEGER, iBase_DOUBLE, iBase_ENTITY_HANDLE, iBase_BYTES};
int aiTagSizes[] = {1, 2, 3, 4};
int aiTagBytes[] = {sizeof(int), sizeof(double),
  sizeof(void*), sizeof(char)};
#else
extern int iNPassed, iNFailed, iNNotImpl, iNIntPassed, iNIntFailed;
extern bool qAllPassed, qErrTests, qSilent, qOutFile;
extern bool qSaveLoad, qSaveLoadClean;

extern int a2iAdjTable[4][4];

extern iBase_EntitySetHandle pvRootSet;
extern bool a2qLegalTypeAndTopo[5][12];
extern int a2iEntCount[5][12];
extern char* astrTypes[];
extern char* astrTopologies[];
extern char* astrTagNames[];
extern int aiTagNameLen[];
extern int aiTagSizes[];
extern int aiTagBytes[];
extern int aTagTypes[];

#endif

void
qFailed();
bool
qCheckNumOfType(iMesh_Instance instance, iBase_EntitySetHandle SH, int iType,
		int iNumEntsExpected);
bool
qCheckNumOfTopo(iMesh_Instance instance, iBase_EntitySetHandle SH, int iTopo,
		int iNumEntsExpected);
void
vCheckBoolOpResult(iMesh_Instance instance, iBase_EntitySetHandle SH,
		   int iShouldBeList, iBase_EntitySetHandle aSHExpected[],
		   int iNumSetsExpected, iBase_EntityHandle aEHExpected[],
		   int iNumEntsExpected);
void
CleanSaveLoad();
void
CleanSaveLoadOnSignal(int);
void
vDoSaveLoad(iMesh_Instance *instancep, iBase_EntitySetHandle *rootSet);
void
vCheckEntities(iMesh_Instance instance, iBase_EntitySetHandle SH, int iType,
	       int iTopo, iBase_EntityHandle aEHExpected[],
	       int iNumEntsExpected);
void
vTestIterWithModify(const int iDim);
void
qIssueLoadMessage(const char *filename, int err);
void
vRecoverSetsParentChild(iMesh_Instance instance,
			iBase_EntitySetHandle pvRootSet,
			iBase_EntitySetHandle apvEntSets[]);
void
vRecoverSetsContain(iMesh_Instance instance, iBase_EntitySetHandle pvRootSet,
		    iBase_EntitySetHandle apvEntSets[]);
void
vRecoverSetsFromEntCounts(iMesh_Instance instance,
			  iBase_EntitySetHandle pvRootSet,
			  iBase_EntitySetHandle *apvSets,
			  const int aiNumTypeEnts[], const int aiNumTopoEnts[]);
void
vRecoverTags(iMesh_Instance instance, iBase_EntitySetHandle pvRootSet,
	     iBase_TagHandle aTags[], iBase_EntityHandle *aEnt = 0);
void
vSecondAdjacencyTest(iMesh_Instance instance);
void
vSecondAdjacencyArrayTest(iMesh_Instance instance);
void
vEntityReciprocityTest(iMesh_Instance instance, iBase_EntityType eET1,
		       iBase_EntityType eET2);
void
vArrayReciprocityTest(iMesh_Instance instance, iBase_EntityType eET1,
		      iBase_EntityType eET2);
void
vTestMeshInterface(iMesh_Instance& instance, string strFilename);
void
vTestEntityInterface(iMesh_Instance& instance);
void
vTestArrInterface(iMesh_Instance& instance);
void
vTestModifyInterface(iMesh_Instance& instance);
void
vTestEntSetInterface(iMesh_Instance& instance);
void
vTestMeshWithSets(iMesh_Instance& instance);
void
vTestEntityWithSets(iMesh_Instance& instance);
void
vTestArrayWithSets(iMesh_Instance& instance);
void
vTestEntTagInterface(iMesh_Instance& instance);
void
vTestArrTagInterface(iMesh_Instance& instance);
void
vTestSetRelationInterface(iMesh_Instance& instance);
void
vTestBoolOpsInterface(iMesh_Instance& instance);
void
vTestSetTagInterface(iMesh_Instance& instance);
void
vTestLowerDimensionedCoordinates();
