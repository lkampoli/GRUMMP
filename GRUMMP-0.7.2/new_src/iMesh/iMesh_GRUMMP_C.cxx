// This file contains the C API implementation of the iMesh interface
// with the GRUMMP backend.  Many calls are dispatched to the
// iMesh_GRUMMP class (which will probably eventually morph to
// GRUMMP::Mesh) and its derived classes.
#define IMESH_MAIN
#include "iMesh.h"
#include "iMesh_GRUMMP.hh"
#include "iMesh_GRUMMP_misc.hh"

using namespace ITAPS_GRUMMP;

#include <vector>
#include <string.h>
#include "GR_Entity.h"
#include "GR_Vertex.h"
#include "GR_Geometry.h"
#include "GR_Mesh2D.h"
#include "GR_SurfMesh.h"
#include "GR_VolMesh.h"
#include "GR_iMesh_Classes.h"
#include "iMesh_GRUMMP_misc.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"
#include "iMesh_GRUMMP_EntitySet.hh"

#ifndef VTK_TETRA
#define VTK_VERTEX 1
#define VTK_LINE 3
#define VTK_TRIANGLE 5
#define VTK_QUAD 9
#define VTK_TETRA 10
#define VTK_HEXAHEDRON 12
#define VTK_WEDGE 13
#define VTK_PYRAMID 14
#endif

#define DISCARD_LINE do {			\
    lineNum++;					\
    int i__i = fscanf(pFInFile, "%*[^\n]\n");	\
    i__i++;					\
  } while(0)

#define CHECK_READ(b, a) do {						\
    int _res_ = (a);							\
    lineNum++;								\
    if (_res_ != b) {							\
      logMessage(0, "File format error in function %s\n", __func__);	\
      logMessage(0, "File %s, line %zd: \n", name, lineNum);		\
      *err = iBase_FAILURE;						\
      instance->pGB->setError(*err); return;					\
    }									\
  } while(0)

#define CHECK_RANGE(val_, min_, max_) do {				\
    if (val_ < min_ || val_ > ssize_t(max_)) {				\
      vMessage(0, "Index range error in function %s\n", __func__);	\
      vMessage(0, "File %s, line %zd: \n", strFileName, lineNum);	\
      *err = iBase_FAILURE;						\
      instance->pGB->setError(*err); return;					\
    }									\
  } while(0)

#define CHECK_POS_RANGE(val_, max_) do {					\
    if (val_ > max_) {							\
      vMessage(0, "Index range error in function %s\n", __func__);	\
      vMessage(0, "File %s, line %zd: \n", strFileName, lineNum);	\
      *err = iBase_FAILURE;						\
      instance->pGB->setError(*err); return;					\
    }									\
  } while(0)

#define CHECK_VERTEX_RANGE(val_) CHECK_POS_RANGE(val_, nVerts)

static bool
isRootSet(iMesh_Instance, iBase_EntitySetHandle entity_set_handle)
{
  return (entity_set_handle == GR_rootSet);
}

static void
changeFileSuffix(const char strFileName[], const char strOldSuffix[],
		 const char strNewSuffix[], char strNewFileName[],
		 const int maxSize)
{
  char *tmpName = new char[maxSize];
  snprintf(tmpName, maxSize, "%s", strFileName);
  char* suffixLoc = strcasestr(tmpName, strOldSuffix);
  if (suffixLoc != NULL)
    *suffixLoc = '\0';
  snprintf(strNewFileName, maxSize, "%s%s", tmpName, strNewSuffix);
}

extern "C"
{

// Error functions
  void iMesh_getErrorType(iMesh_Instance instance, int *err)
  {
    *err = instance->getError();
  }

  void iMesh_getDescription(iMesh_Instance instance, char *descr, int descr_len)
  {
    const char *strDescriptions[] =
      { "Success!", "Mesh already loaded", "No mesh data", "File not found",
	  "File write error", "Nil array", "Bad array size",
	  "Bad array dimension", "Invalid entity handle",
	  "Invalid entity count", "Invalid entity type",
	  "Invalid entity topology", "Bad type and topo",
	  "Entity creation error", "Invalid tag handle", "Tag not found",
	  "Tag already exists", "Tag in use", "Invalid entityset handle",
	  "Invalid iterator handle", "Invalid argument",
	  "Memory allocation failed", "Not supported", "Failure" };
    snprintf(descr, descr_len - 1, "%s", strDescriptions[instance->getError()]);
  }

  void iMesh_newMesh(const char *options_in, iMesh_Instance *instance, int *err,
		     int options_len)
  {
    // From old Mesh factory
    *err = iBase_SUCCESS;
    bool qHas2D = false, qHas3D = true, qHasSurface = false;
    void *pNoWhere = NULL;
    char options[1024];
    NULL_TERMINATE(options_in, options_len, options);
    for (int i = strlen(options) - 1; i >= 0; i--) {
      options[i] = tolower(options[i]);
    }

    if (options != NULL && options_len != 0) {
      qHas2D = (pNoWhere != strstr(options, "grummp:2d"));
      qHas3D = (pNoWhere != strstr(options, "grummp:3d"));
      qHasSurface = (pNoWhere != strstr(options, "grummp:surface"));
    }
    if (qHas2D) {
      *instance = new iMesh_Instance_Private(iMesh_Base::ePlane);
    }
    else if (qHasSurface) {
      *instance = new iMesh_Instance_Private(iMesh_Base::eSurface);
    }
    else { // Default to 3D
      assert(qHas3D);
      *instance = new iMesh_Instance_Private(iMesh_Base::eVolume);
    }
  }

  void iMesh_dtor(iMesh_Instance instance, int *err)
  {
    *err = iBase_SUCCESS;
    delete instance;
  }

  void iMesh_load(iMesh_Instance instance,
		  const iBase_EntitySetHandle entity_set_handle,
		  const char *name_in, const char* options_in, int *err,
		  int name_len, int options_len)
  {
    iBase_EntitySetHandle root_set;
    *err = iBase_SUCCESS;
    iMesh_getRootSet(instance, &root_set, err);
    if (*err != iBase_SUCCESS)
      return;
    char name[1024], baseName[1024], options[1024];
    NULL_TERMINATE(name_in, name_len, baseName);
    NULL_TERMINATE(options_in, options_len, options);
    for (int i = strlen(options) - 1; i >= 0; i--) {
      options[i] = tolower(options[i]);
    }

    bool qSilent = (strstr(options, "grummp:silent") != NULL);

    // Check options to see whether to treat this as a VTK file.
    if ((strstr(options, "grummp:vtk") != NULL) // Supposed to be VTK
    && (strlen(baseName) >= 4)
	&& (strcasecmp(baseName + strlen(baseName) - 4, ".vtk")))
	// Doesn't end in .vtk already
	{
      sprintf(name, "%s.vtk", baseName);
    }
    else {
      sprintf(name, "%s", baseName);
    }

    if (entity_set_handle == root_set) {
      if (instance->pGB && instance->pGB->pMesh()
	  && instance->pGB->pMesh()->getNumVerts() != 0) {
	*err = iBase_MESH_ALREADY_LOADED;
	instance->pGB->setError(*err);
	return;
      }

      if (strcasecmp(name + strlen(name) - 4, ".vtk")) {
	// Call the GRUMMP file loader
	if (strstr(name, ".mesh") != NULL) {
	  // In this case, it's a 2D mesh.
	  iMesh_setGeometricDimension(instance, 2, err);
	  // Should never hit an error
	}
	else {
	  iMesh_setGeometricDimension(instance, 3, err);
	  // Should never hit an error
	}

	instance->pGB->pMesh()->readFromFile(name);
	return;
      }

      FILE *pFInFile = fopen(name, "r");
      *err = iBase_SUCCESS;
      if (NULL == pFInFile) {
	if (instance->pGB == NULL) {
	  // This makes sure that there's an object to set the
	  // error on.
	  iMesh_setGeometricDimension(instance, 3, err);
	}
	*err = iBase_FILE_NOT_FOUND;
	instance->pGB->setError(*err);
	return;
      }
      if (!qSilent) {
	logMessage(1, "Reading VTK legacy input\n");
      }

      // Discard the four line header
      size_t lineNum = 1;
      DISCARD_LINE
      ;
      DISCARD_LINE
      ;
      DISCARD_LINE
      ;
      DISCARD_LINE
      ;

      GR_index_t nVerts;
      CHECK_READ(1, (fscanf(pFInFile, "POINTS %u %*s\n", &nVerts)));
      if (nVerts <= 0) {
	*err = iBase_BAD_ARRAY_SIZE;
	instance->pGB->setError(*err);
	return;
      }

      double (*coordsTemp)[3], minZ = 1.e300, maxZ = -1.e300;
      // coordsTemp is leaked if there's a read error.
      coordsTemp = new double[nVerts][3];
      for (GR_index_t i = 0; i < nVerts; i++) {
	CHECK_READ(
	    3,
	    (fscanf(pFInFile, "%lf%lf%lf\n", coordsTemp[i], coordsTemp[i] + 1,
		    coordsTemp[i] + 2)));
	minZ = min(minZ, coordsTemp[i][2]);
	maxZ = max(maxZ, coordsTemp[i][2]);
      }

      int geomDim = -1;
      if (fabs(minZ - maxZ) < 1.e-10) {
	// All have the same z coord, so it's 2D
	iMesh_setGeometricDimension(instance, 2, err);
	geomDim = 2;
      }
      else {
	iMesh_setGeometricDimension(instance, 3, err);
	geomDim = 3;
      }

      // Read in the actual connectivity info
      int iConnectSize, iNEnts;
      CHECK_READ(2,
		 (fscanf(pFInFile, "CELLS %d %d\n", &iNEnts, &iConnectSize)));
      iConnectSize -= iNEnts;

      // aiConnect and aiOffset are leaked if there's a read error.
      int *aiConnect = new int[iConnectSize];
      int *aiOffset = new int[iNEnts + 1];
      int iIndex = 0;
      aiOffset[0] = 0;
      for (int i = 0; i < iNEnts; i++) {
	int iVertsThisEnt;
	CHECK_READ(1, (fscanf(pFInFile, "%d", &iVertsThisEnt)));
	aiOffset[i + 1] = aiOffset[i] + iVertsThisEnt;
	for (int j = 0; j < iVertsThisEnt; j++) {
	  CHECK_READ(1, (fscanf(pFInFile, "%d", &(aiConnect[iIndex]))));
	  iIndex++;
	}
	assert(iIndex == aiOffset[i + 1]);
      }

      // aiTopo is leaked if there's a read error
      iMesh_EntityTopology *aiTopo = new iMesh_EntityTopology[iNEnts];
      // Read the cell type data as well
      {
	int iCheck;
	CHECK_READ(1, (fscanf(pFInFile, "%*s %d\n", &iCheck)));
	assert(iCheck == iNEnts);
      }
      int edgesRead = 0, trisRead = 0, quadsRead = 0, tetsRead = 0,
	  pyrsRead = 0, prismsRead = 0, hexesRead = 0;
      for (int i = 0; i < iNEnts; i++) {
	int iVTKtype;
	CHECK_READ(1, (fscanf(pFInFile, "%d", &iVTKtype)));
	switch (iVTKtype)
	  {
	  case VTK_LINE:
	    aiTopo[i] = iMesh_LINE_SEGMENT;
	    edgesRead++;
	    break;
	  case VTK_TRIANGLE:
	    aiTopo[i] = iMesh_TRIANGLE;
	    trisRead++;
	    break;
	  case VTK_QUAD:
	    aiTopo[i] = iMesh_QUADRILATERAL;
	    quadsRead++;
	    break;
	  case VTK_TETRA:
	    aiTopo[i] = iMesh_TETRAHEDRON;
	    tetsRead++;
	    break;
	  case VTK_PYRAMID:
	    aiTopo[i] = iMesh_PYRAMID;
	    pyrsRead++;
	    break;
	  case VTK_WEDGE:
	    aiTopo[i] = iMesh_PRISM;
	    prismsRead++;
	    break;
	  case VTK_HEXAHEDRON:
	    aiTopo[i] = iMesh_HEXAHEDRON;
	    hexesRead++;
	    break;
	  default:
	    *err = iBase_NOT_SUPPORTED;
	    instance->pGB->setError(*err);
	    return;
	    break;
	  }
      }
      fclose(pFInFile);
      if (!qSilent) {
	logMessage(
	    1,
	    "Read %d edges, %d tris, %d quads, %d tets, %d pyrs, %d prisms, %d hexes\n",
	    edgesRead, trisRead, quadsRead, tetsRead, pyrsRead, prismsRead,
	    hexesRead);
      }
      if ((tetsRead + pyrsRead + prismsRead + hexesRead == 0)
	  && (geomDim == 3)) {
	// This looks like a surface mesh.
	geomDim = -3;
	iMesh_setGeometricDimension(instance, -3, err);
      }

      // Create entities.  Since this is going through the iMesh calling
      // sequence, 2D and 3D can be handled together.  Insert in order, even
      // though this constrains the way that set and tag files are laid out.

      // If any creation fails, forget about the whole thing.

      // The 2D code ignores the extra arg, and the base class declaration
      // requires three args
      iBase_EntityHandle *allVerts = new iBase_EntityHandle[nVerts];
      // allVerts is leaked if there's a read error
      for (GR_index_t i = 0; i < nVerts; i++) {
	iMesh_createVtx(instance, coordsTemp[i][0], coordsTemp[i][1],
			coordsTemp[i][2], &allVerts[i], err);
      }
      delete[] coordsTemp;

      iBase_EntityHandle verts[8];
      for (int i = 0; i < iNEnts && *err == iBase_SUCCESS; i++) {
	if (aiTopo[i] == iMesh_LINE_SEGMENT)
	  continue; // Don't try to
	// build edges
	int ii = 0, iVertsThisEnt = aiOffset[i + 1] - aiOffset[i];
	for (; ii < iVertsThisEnt; ii++) {
	  verts[ii] = allVerts[aiConnect[aiOffset[i] + ii]];
	}
	if (aiTopo[i] == iMesh_PRISM) {
	  // VTK flips orientation of triangular faces in its prisms
	  iBase_EntityHandle tmpVert = verts[2];
	  verts[2] = verts[1];
	  verts[1] = tmpVert;

	  tmpVert = verts[5];
	  verts[5] = verts[4];
	  verts[4] = tmpVert;
	}
	iBase_EntityHandle newEnt;
	int status;
	iMesh_createEnt(instance, aiTopo[i], verts, iVertsThisEnt, &newEnt,
			&status, err);
      }

      delete[] aiTopo;
      delete[] aiOffset;
      delete[] aiConnect;
      delete[] allVerts;

      // Need to create bdry faces, because the VTK file doesn't have any of
      // these.
      instance->pGB->pMesh()->checkForITAPSBdryErrors();

      if (!qSilent) {
	int iNLineSeg, iNTri, iNQuad, iNTet, iNPyr, iNPrism, iNHex;
	iMesh_getNumOfTopo(instance, root_set, iMesh_LINE_SEGMENT, &iNLineSeg,
			   err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_TRIANGLE, &iNTri, err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_QUADRILATERAL, &iNQuad,
			   err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_TETRAHEDRON, &iNTet, err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_PYRAMID, &iNPyr, err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_PRISM, &iNPrism, err);
	iMesh_getNumOfTopo(instance, root_set, iMesh_HEXAHEDRON, &iNHex, err);
	printf(
	    "Created %d edges, %d tris, %d quads, %d tets, %d pyrs, %d prisms, %d hexes\n",
	    iNLineSeg, iNTri, iNQuad, iNTet, iNPyr, iNPrism, iNHex);

      }
    }
    else {
      *err = iBase_NOT_SUPPORTED;
      instance->pGB->setError(*err);
      return;
    }

    // Now read the sets.

    // Get a file name for the set file and open it.
    char setName[1024];
    snprintf(setName, 1024, "%s.sets", baseName);
    FILE* pFSetFile = fopen(setName, "r");
    if (!pFSetFile) {
      // Strip off the suffix from the file name and replace it with
      // .sets; then see if this file exists.
      changeFileSuffix(name, "vtk", "sets", setName, 1024);
      pFSetFile = fopen(setName, "r");
    }

    // Get a file name for the tag file and open it.
    char tagFileName[1024];
    snprintf(tagFileName, 1024, "%s.tags", baseName);
    FILE* pFTagFile = fopen(tagFileName, "r");
    if (!pFTagFile) {
      // Strip off the suffix from the file name and replace it with
      // .tags; then see if this file exists.
      changeFileSuffix(name, "vtk", "tags", tagFileName, 1024);
      pFTagFile = fopen(tagFileName, "r");
    }

    if (pFSetFile || pFTagFile) {
      // Need a full list of all entity handles
      iBase_EntityHandle *allHandles;
      int allHandles_allocated = 0, allHandles_size;
      iMesh_getEntities(instance, root_set, iBase_ALL_TYPES,
			iMesh_ALL_TOPOLOGIES, &allHandles,
			&allHandles_allocated, &allHandles_size, err);

      int numSets = 0;
      iBase_EntitySetHandle *allSets = NULL;
      // Don't try to read the file if it doesn't exist.
      if (pFSetFile) {
	int res = fscanf(pFSetFile, "Total number of sets: %d\n", &numSets);
	assert(res == 1);
	allSets = new iBase_EntitySetHandle[numSets];
	// Header info
	for (int setInd = 0; setInd < numSets; setInd++) {
	  int junk;
	  char type[20];
	  res = fscanf(pFSetFile, "Set #%d; %20s\n", &junk, type);
	  assert(res == 2);
	  iMesh_createEntSet(instance, (strstr(type, "un") == NULL),
			     &(allSets[setInd]), err);
	  if (*err != iBase_SUCCESS) {
	    delete[] allSets;
	    return;
	  }
	}

	// Now read and add members.
	for (int setInd = 0; setInd < numSets; setInd++) {
	  iBase_EntitySetHandle set = allSets[setInd];
	  int junk;

	  // Start by reading entities
	  res = fscanf(pFSetFile, " Set #%d\nContains entities:", &junk);
	  assert(res == 1 && junk == setInd);
	  int entity, junk2;
	  while (2 == fscanf(pFSetFile, " (%d %d)", &entity, &junk2)) {
	    assert(entity >= 0 && entity < allHandles_size);
	    iMesh_addEntToSet(instance, allHandles[entity], set, err);
	  }

	  // Now read set members
	  int numMembers;
	  res = fscanf(pFSetFile, " Set members (%d):", &numMembers);
	  assert(res == 1);
	  for (int ii = 0; ii < numMembers; ii++) {
	    res = fscanf(pFSetFile, "%d", &junk);
	    assert(res == 1);
	    assert(junk >= 0 && junk < numSets);
	    iMesh_addEntSet(instance, allSets[junk], set, err);
	  }

	  // Now read children
	  int numChildren;
	  res = fscanf(pFSetFile, " Child sets (%d):", &numChildren);
	  assert(res == 1);
	  for (int ii = 0; ii < numChildren; ii++) {
	    res = fscanf(pFSetFile, "%d", &junk);
	    assert(res == 1);
	    assert(junk >= 0 && junk < numSets);
	    iMesh_addPrntChld(instance, set, allSets[junk], err);
	  }

	  // Now read parents, though nothing needs to be done.
	  int numParents;
	  res = fscanf(pFSetFile, " Parent sets (%d):", &numParents);
	  assert(res == 1);
	  for (int ii = 0; ii < numParents; ii++) {
	    res = fscanf(pFSetFile, "%d", &junk);
	    assert(res == 1);
	    assert(junk >= 0 && junk < numSets);
	  }
	}
	fclose(pFSetFile);
      }

      // Don't read the tag file if it doesn't exist
      if (pFTagFile) {
	// Read tag header info, and create tags.
	int numTags;
	int res = fscanf(pFTagFile, " Total number of tags: %d", &numTags);
	assert(res == 1);
	iBase_TagHandle allTags[numTags];
	int allTypes[numTags], allSizes[numTags];

	for (int tag = 0; tag < numTags; tag++) {
	  int type, size;
	  res = fscanf(pFTagFile, " Tag #%*d: Type: %d Size: %d Name: ", &type,
		       &size);
	  assert(res == 2);
	  char tagName[100];
	  char *chRes = fgets(tagName, 100, pFTagFile);
	  assert(chRes == tagName);
	  int len = strlen(tagName) - 1;
	  tagName[len - 1] = 0;
	  iMesh_createTag(instance, tagName, size, type, &allTags[tag], err,
			  len);
	  allTypes[tag] = type;
	  allSizes[tag] = size;
	  if (*err != iBase_SUCCESS)
	    return;
	}

	int entInd, entTags, scanRes, entsTagged = 0;
	scanRes = fscanf(pFTagFile, " Entity %d has %d tags.\n", &entInd,
			 &entTags);
	while (scanRes == 2) {
	  entsTagged++;
	  iBase_EntityHandle ent = allHandles[entInd];
	  for (int iTag = 0; iTag < entTags; iTag++) {
	    int tagInd;
	    res = fscanf(pFTagFile, " Tag %d. Value(s): ", &tagInd);
	    assert(res == 1);
	    iBase_TagHandle tag = allTags[tagInd];
	    int type = allTypes[tagInd];
	    int size = allSizes[tagInd];
	    switch (type)
	      {
	      case iBase_INTEGER:
		{
		  int data[size];
		  for (int ii = 0; ii < size; ii++) {
		    res = fscanf(pFTagFile, "%d", &(data[ii]));
		    assert(res == 1);
		  }
		  iMesh_setIntArrData(instance, &ent, 1, tag, data, size, err);
		}
		break;
	      case iBase_DOUBLE:
		{
		  double data[size];
		  for (int ii = 0; ii < size; ii++) {
		    res = fscanf(pFTagFile, "%lf", &(data[ii]));
		    assert(res == 1);
		  }
		  iMesh_setDblArrData(instance, &ent, 1, tag, data, size, err);
		}
		break;
	      case iBase_ENTITY_HANDLE:
		{
		  iBase_EntityHandle data[size];
		  for (int ii = 0; ii < size; ii++) {
		    int tagValInd;
		    res = fscanf(pFTagFile, "%d", &tagValInd);
		    assert(res == 1);
		    data[ii] = allHandles[tagValInd];
		  }
		  iMesh_setEHArrData(instance, &ent, 1, tag, data, size, err);
		}
		break;
	      case iBase_ENTITY_SET_HANDLE:
		{
		  if (allSets != NULL) {
		    iBase_EntitySetHandle data[size];
		    for (int ii = 0; ii < size; ii++) {
		      int tagValInd;
		      res = fscanf(pFTagFile, "%d", &tagValInd);
		      assert(res == 1);
		      data[ii] = allSets[tagValInd];
		    }
		    iMesh_setESHArrData(instance, &ent, 1, tag, data, size,
					err);
		  }
		  else {
		    // Not sure how we can be reading set tag handles
		    // without any sets.  Definitely implies an error in
		    // whatever wrote the files.
		    assert(0);
		  }
		}
		break;
	      case iBase_BYTES:
		{
		  char data[size];
		  for (int ii = 0; ii < size; ii++) {
		    res = fscanf(pFTagFile, "%c", &(data[ii]));
		    assert(res == 1);
		  }
		  iMesh_setData(instance, ent, tag, data, size, err);
		}
		break;
	      } // End of switch over tag type
	  } // End of loop over tags for this entity
	  scanRes = fscanf(pFTagFile, " Entity %d has %d tags.\n", &entInd,
			   &entTags);
	} // End reading tag data for all entities.
	if (!qSilent) {
	  logMessage(2, "Found %d entities with tags.\n", entsTagged);
	}

	// Don't bother with set tags if there are no sets
	if (numSets > 0) {
	  // Now repeat for sets.
	  int setInd, setTags, setsTagged = 0;
	  scanRes = fscanf(pFTagFile, "Set %d has %d tags.\n", &setInd,
			   &setTags);
	  while (scanRes == 2) {
	    setsTagged++;
	    iBase_EntitySetHandle set = allSets[setInd];
	    for (int iTag = 0; iTag < setTags; iTag++) {
	      int tagInd;
	      res = fscanf(pFTagFile, " Tag %d. Value(s): ", &tagInd);
	      assert(res == 1);
	      iBase_TagHandle tag = allTags[tagInd];
	      int type = allTypes[tagInd];
	      int size = allSizes[tagInd];
	      assert(size == 1 || type == iBase_BYTES);
	      switch (type)
		{
		case iBase_INTEGER:
		  {
		    int data;
		    res = fscanf(pFTagFile, "%d", &data);
		    assert(res == 1);
		    iMesh_setEntSetIntData(instance, set, tag, data, err);
		  }
		  break;
		case iBase_DOUBLE:
		  {
		    double data;
		    res = fscanf(pFTagFile, "%lf", &data);
		    assert(res == 1);
		    iMesh_setEntSetDblData(instance, set, tag, data, err);
		  }
		  break;
		case iBase_ENTITY_HANDLE:
		  {
		    int tagValInd;
		    res = fscanf(pFTagFile, "%d", &tagValInd);
		    assert(res == 1);
		    iBase_EntityHandle data = allHandles[tagValInd];
		    iMesh_setEntSetEHData(instance, set, tag, data, err);
		  }
		  break;
		case iBase_ENTITY_SET_HANDLE:
		  {
		    int tagValInd;
		    res = fscanf(pFTagFile, "%d", &tagValInd);
		    assert(res == 1);
		    iBase_EntitySetHandle data = allSets[tagValInd];
		    iMesh_setEntSetESHData(instance, set, tag, data, err);
		  }
		  break;
		case iBase_BYTES:
		  {
		    char data[size];
		    for (int ii = 0; ii < size; ii++) {
		      res = fscanf(pFTagFile, "%c", &(data[ii]));
		      assert(res == 1);
		    }
		    iMesh_setEntSetData(instance, set, tag, data, size, err);
		  }
		  break;
		} // End of switch over tag type
	    } // End of loop over tags for this setity
	    scanRes = fscanf(pFTagFile, " Set %d has %d tags.\n", &setInd,
			     &setTags);
	  } // End reading tag data for all sets
	    // Now repeat for root set
	  scanRes = fscanf(pFTagFile, "Root set has %d tags.\n", &setTags);
	  if (scanRes == 1) {
	    setsTagged++;
	    for (int iTag = 0; iTag < setTags; iTag++) {
	      int tagInd;
	      res = fscanf(pFTagFile, " Tag %d. Value(s): ", &tagInd);
	      iBase_TagHandle tag = allTags[tagInd];
	      int type = allTypes[tagInd];
	      int size = allSizes[tagInd];
	      assert(size == 1 || type == iBase_BYTES);
	      switch (type)
		{
		case iBase_INTEGER:
		  {
		    int data;
		    res = fscanf(pFTagFile, "%d", &data);
		    iMesh_setEntSetIntData(instance, root_set, tag, data, err);
		  }
		  break;
		case iBase_DOUBLE:
		  {
		    double data;
		    res = fscanf(pFTagFile, "%lf", &data);
		    iMesh_setEntSetDblData(instance, root_set, tag, data, err);
		  }
		  break;
		case iBase_ENTITY_HANDLE:
		  {
		    int tagValInd;
		    res = fscanf(pFTagFile, "%d", &tagValInd);
		    iBase_EntityHandle data = allHandles[tagValInd];
		    iMesh_setEntSetEHData(instance, root_set, tag, data, err);
		  }
		  break;
		case iBase_ENTITY_SET_HANDLE:
		  {
		    int tagValInd;
		    res = fscanf(pFTagFile, "%d", &tagValInd);
		    iBase_EntitySetHandle data = allSets[tagValInd];
		    iMesh_setEntSetESHData(instance, root_set, tag, data, err);
		  }
		  break;
		case iBase_BYTES:
		  {
		    char data[size];
		    for (int ii = 0; ii < size; ii++) {
		      res = fscanf(pFTagFile, "%c", &(data[ii]));
		    }
		    iMesh_setEntSetData(instance, root_set, tag, data, size,
					err);
		  }
		  break;
		} // End of switch over tag type
	    } // End of loop over tags for the root set
	  } // End reading tag data for root set
	  if (!qSilent) {
	    logMessage(2, "Found %d sets with tags.\n", setsTagged);
	  }

	} // Read set tags only if there are sets.
	fclose(pFTagFile);
      } // Done reading tag file
      if (allSets)
	delete[] allSets;
      if (allHandles)
	free(allHandles);
    } // Done reading sets and tags

    iBase_TagHandle TH;
    iMesh_getTagHandle(instance, "fixed", &TH, err, 5);
    if (*err == iBase_SUCCESS) {
      // In this case, we've already got the bloody thing; leave well
      // enough alone, if only because its semantics could be
      // different.
      return;
    }

    // Set up a tag for vertex type data.  Mesquite wants this
    // exposed, and it does no real harm.
    iBase_TagHandle TH2;
    iMesh_createTag(instance, "fixed", 1, iBase_INTEGER, &TH2, err, 5);
    if (*err != iBase_SUCCESS)
      return;

    // Now set data for it, by looping over verts and transcribing.
    for (unsigned int iV = 0; iV < instance->pGB->getNumVertices(); iV++) {
      Vert *pV = instance->pGB->pMesh()->getVert(iV);
      int iType = pV->getVertType();
      if (iType == Vert::eBdryApex || iType == Vert::eBdryCurve
	  || iType == Vert::eBdryTwoSide || iType == Vert::eBdry)
	iMesh_setIntData(instance, pV, TH2, 1, err);
      else
	iMesh_setIntData(instance, pV, TH2, 0, err);
      if (*err != iBase_SUCCESS)
	return;
    }
  }

  void iMesh_save(iMesh_Instance instance,
		  const iBase_EntitySetHandle entity_set_handle,
		  const char *name_in, const char* options_in, int *err,
		  const int name_len, int options_len)
  {
    // Base public, with switch
    iBase_EntitySetHandle root_set;
    *err = iBase_SUCCESS;
    iMesh_getRootSet(instance, &root_set, err);
    if (*err != iBase_SUCCESS)
      return;
    char name[1024], options[1024];
    NULL_TERMINATE(name_in, name_len, name);
    NULL_TERMINATE(options_in, options_len, options);
    for (int i = strlen(options) - 1; i >= 0; i--) {
      options[i] = tolower(options[i]);
    }

    // At this point, there's nothing to silence.
//     bool qSilent = (NULL != strstr(options, "grummp:silent"));

    if (entity_set_handle == root_set) {
      instance->pGB->save(name, options, err);

      // Writing sets.  This requires finding indices for all the
      // entities.  Basically, this uses a binary search in the
      // allHandles array to find where a particular handle is, then
      // pointer arithmetic to convert that to an integer.

      // Get a file name for the set file and open it.
      char setName[1024];
      snprintf(setName, 1024, "%s.sets", name);
      FILE *pFOutFile = fopen(setName, "w");

      iBase_EntitySetHandle *allSetHandles = NULL;
      int allSetHandles_allocated = 0, allSetHandles_size;
      iMesh_getEntSets(instance, root_set, 0, &allSetHandles,
		       &allSetHandles_allocated, &allSetHandles_size, err);

      // Start with some summary info about sets, so that they can be
      // created up front on read.
      fprintf(pFOutFile, "Total number of sets: %d\n", allSetHandles_size);

      for (int setInd = 0; setInd < allSetHandles_size; setInd++) {
	iBase_EntitySetHandle setHandle = allSetHandles[setInd];
	int isList;
	iMesh_isList(instance, setHandle, &isList, err);
	fprintf(pFOutFile, "Set #%d; %s\n", setInd,
		isList ? "ordered" : "unordered");
      }

      iBase_EntityHandle *allHandles = NULL;
      int allHandles_allocated = 0, allHandles_size;
      iMesh_getEntities(instance, root_set, iBase_ALL_TYPES,
			iMesh_ALL_TOPOLOGIES, &allHandles,
			&allHandles_allocated, &allHandles_size, err);
      assert(*err == iBase_SUCCESS);

      // Now write membership and parent-child data
      for (int setInd = 0; setInd < allSetHandles_size; setInd++) {
	iBase_EntitySetHandle setHandle = allSetHandles[setInd];
	fprintf(pFOutFile, "Set #%d\n", setInd);
	// Now iterate over the set to find its entity members
	iBase_EntityIterator entIter;
	iMesh_initEntIter(instance, setHandle, iBase_ALL_TYPES,
			  iMesh_ALL_TOPOLOGIES, 0, &entIter, err);

	iBase_EntityHandle entHandle;
	int hasData;
	iMesh_getNextEntIter(instance, entIter, &entHandle, &hasData, err);
	fprintf(pFOutFile, "Contains entities: ");
	while (hasData) {
	  int index = findArrayIndex(entHandle, allHandles,
				     allHandles + allHandles_size);
	  int topo;
	  iMesh_getEntTopo(instance, entHandle, &topo, err);
	  fprintf(pFOutFile, "(%d %d)", index, topo);
	  iMesh_getNextEntIter(instance, entIter, &entHandle, &hasData, err);
	}
	fprintf(pFOutFile, "\n");
	iMesh_endEntIter(instance, entIter, err);

	// Now its set members:
	{
	  iBase_EntitySetHandle *mySetHandles = NULL;
	  int mySetHandles_allocated = 0, mySetHandles_size;
	  iMesh_getEntSets(instance, setHandle, 0, &mySetHandles,
			   &mySetHandles_allocated, &mySetHandles_size, err);
	  fprintf(pFOutFile, "Set members (%d): ", mySetHandles_size);
	  if (mySetHandles_size > 0) {
	    for (int ii = 0; ii < mySetHandles_size; ii++) {
	      iBase_EntitySetHandle mySH = mySetHandles[ii];
	      int index = findArrayIndex(mySH, allSetHandles,
					 allSetHandles + allSetHandles_size);
	      fprintf(pFOutFile, "%d ", index);
	    }
	    free(mySetHandles);
	  }
	  fprintf(pFOutFile, "\n");
	}

	// Now child sets:
	{
	  iBase_EntitySetHandle *mySetHandles = NULL;
	  int mySetHandles_allocated = 0, mySetHandles_size;
	  iMesh_getChldn(instance, setHandle, 0, &mySetHandles,
			 &mySetHandles_allocated, &mySetHandles_size, err);
	  fprintf(pFOutFile, "Child sets (%d): ", mySetHandles_size);
	  if (mySetHandles_size > 0) {
	    for (int ii = 0; ii < mySetHandles_size; ii++) {
	      iBase_EntitySetHandle mySH = mySetHandles[ii];
	      int index = findArrayIndex(mySH, allSetHandles,
					 allSetHandles + allSetHandles_size);
	      fprintf(pFOutFile, "%d ", index);
	    }
	    free(mySetHandles);
	  }
	  fprintf(pFOutFile, "\n");
	}

	// Now parent sets:
	{
	  iBase_EntitySetHandle *mySetHandles = NULL;
	  int mySetHandles_allocated = 0, mySetHandles_size;
	  iMesh_getPrnts(instance, setHandle, 0, &mySetHandles,
			 &mySetHandles_allocated, &mySetHandles_size, err);
	  fprintf(pFOutFile, "Parent sets (%d): ", mySetHandles_size);
	  if (mySetHandles_size > 0) {
	    for (int ii = 0; ii < mySetHandles_size; ii++) {
	      iBase_EntitySetHandle mySH = mySetHandles[ii];
	      int index = findArrayIndex(mySH, allSetHandles,
					 allSetHandles + allSetHandles_size);
	      fprintf(pFOutFile, "%d ", index);
	    }
	    free(mySetHandles);
	  }
	  fprintf(pFOutFile, "\n");
	} // Done writing this set
	fprintf(pFOutFile, "\n");
      } // Done writing all sets.
      fclose(pFOutFile);

      // Tags have to be done from somewhere with more direct access to
      // the tag data

      // Get a file name for the tag file and open it.
      char tagName[1024];
      snprintf(tagName, 1024, "%s.tags", name);

      instance->pGB->writeTags(tagName, allHandles, allHandles_size,
			       allSetHandles, allSetHandles_size);

      // Clean up
      if (allHandles)
	free(allHandles);
      if (allSetHandles)
	free(allSetHandles);
    }
    else {
      *err = iBase_NOT_SUPPORTED;
      instance->pGB->setError(*err);
    }
  }

  void iMesh_getRootSet(iMesh_Instance, iBase_EntitySetHandle *root_set,
			int *err)
  {
    // Base public
    *err = iBase_SUCCESS;
    (*root_set) = GR_rootSet;
  }

  void iMesh_getGeometricDimension(iMesh_Instance instance, int *geom_dim,
				   int *err)
  {
    // Inline specialized func
    *err = iBase_SUCCESS;
    *geom_dim = instance->pGB->getGeometricDimension();
  }

  void iMesh_setGeometricDimension(iMesh_Instance instance, int geom_dim,
				   int *err)
  {
    *err = iBase_SUCCESS;
    if (instance->pGB == NULL || instance->pGB->pMesh() == NULL
	|| instance->pGB->pMesh()->getNumVerts() == 0) {
      // In this case, kill the old mesh database and create a new one.
      switch (geom_dim)
	{
	case -3:
	  instance->vInitiMesh(iMesh_Base::eSurface);
	  break;
	case 1:
	  *err = iBase_FAILURE;
	  // Make it so that there is at least -some-
	  // mesh available...
	  instance->vInitiMesh(iMesh_Base::ePlane);
	  instance->pGB->setError(*err);
	  break;
	case 2:
	  instance->vInitiMesh(iMesh_Base::ePlane);
	  break;
	case 3:
	  instance->vInitiMesh(iMesh_Base::eVolume);
	  break;
	}
    }
    else {
      // In this case, there's already mesh data, and it's too late.
      *err = iBase_FAILURE;
      instance->pGB->setError(*err);
    }
  }

  void iMesh_getDfltStorage(iMesh_Instance /*instance*/, int *order, int *err)
  {
    // Inline base func
    *err = iBase_SUCCESS;
    *order = iBase_INTERLEAVED;
  }

  void iMesh_setAdjTable(iMesh_Instance instance,
  /*intout*/int* adjacency_table,
			 /*in*/int adjacency_table_size, int *err)
  {
    // Derived types
    *err = iBase_SUCCESS;
    if (adjacency_table_size != 16) {
      *err = iBase_INVALID_ARGUMENT;
      instance->pGB->setError(*err);
    }

    if (!instance->pGB) {
      if (adjacency_table[15] == iBase_AVAILABLE
	  || adjacency_table[3] != iBase_UNAVAILABLE
	  || adjacency_table[7] != iBase_UNAVAILABLE
	  || adjacency_table[11] != iBase_UNAVAILABLE
	  || adjacency_table[12] != iBase_UNAVAILABLE
	  || adjacency_table[13] != iBase_UNAVAILABLE
	  || adjacency_table[14] != iBase_UNAVAILABLE) {
	iMesh_setGeometricDimension(instance, 3, err);
      }
      else {
	iMesh_setGeometricDimension(instance, 2, err);
      }
    }
    instance->pGB->setAdjTable(adjacency_table, err);
  }

  void iMesh_getAdjTable(iMesh_Instance instance, int** adjacency_table,
  /*inout*/int* adjacency_table_allocated,
			 /*out*/int* adjacency_table_size, int *err)
  {
    // Derived types
    *err = iBase_SUCCESS;

    *adjacency_table_size = 16;
    TRY_ARRAY_SIZE2(adjacency_table, adjacency_table_allocated,
		    *adjacency_table_size, int);

    instance->pGB->getAdjTable(*adjacency_table);
  }

  void iMesh_getNumOfType(iMesh_Instance instance,
  /*in*/const iBase_EntitySetHandle entity_set_handle,
			  /*in*/const int entity_type, int *num_type, int *err)
  {
    // Base public
    *err = iBase_SUCCESS;
    CHECK_ENTITY_SET(instance, entity_set_handle, err);

    if (isRootSet(instance, entity_set_handle)) {
      // The case for the entire mesh
      switch (entity_type)
	{
	case iBase_VERTEX:
	  *num_type = (instance->pGB->getNumVertices());
	  break;
	case iBase_FACE:
	  *num_type = (instance->pGB->getNumTriangles()
	      + instance->pGB->getNumQuadrilaterals());
	  break;
	case iBase_REGION:
	  *num_type = (instance->pGB->getNumTetrahedra()
	      + instance->pGB->getNumPyramids() + instance->pGB->getNumPrisms()
	      + instance->pGB->getNumHexahedra());
	  break;
	case iBase_EDGE:
	  *num_type = (instance->pGB->getNumLineSegments());
	  break;
	case iBase_ALL_TYPES:
	  *num_type = (instance->pGB->getNumVertices()
	      + instance->pGB->getNumLineSegments()
	      + instance->pGB->getNumTriangles()
	      + instance->pGB->getNumQuadrilaterals()
	      + instance->pGB->getNumTetrahedra()
	      + instance->pGB->getNumPyramids() + instance->pGB->getNumPrisms()
	      + instance->pGB->getNumHexahedra());
	  break;
	default:
	  *err = iBase_INVALID_ENTITY_TYPE;
	  instance->pGB->setError(*err);
	  *num_type = 0;
	  break;
	}
    }
    else {
      // Have an actual entity set argument
      EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set_handle);
      switch (entity_type)
	{
	case iBase_VERTEX:
	  *num_type = (pESB->iNumVerts());
	  break;
	case iBase_EDGE:
	  *num_type = (pESB->iNumEdges());
	  break;
	case iBase_FACE:
	  // Pretty cheesey, but works in both 2D and 3D
	  *num_type = (pESB->iNumTriFaces() + pESB->iNumQuadFaces()
	      + pESB->iNumTriCells() + pESB->iNumQuadCells());
	  break;
	case iBase_REGION:
	  *num_type = (pESB->iNumTetCells() + pESB->iNumPyrCells()
	      + pESB->iNumPrismCells() + pESB->iNumHexCells());
	  break;
	case iBase_ALL_TYPES:
	  *num_type = pESB->iNumVerts() + pESB->iNumEdges()
	      + pESB->iNumTriFaces() + pESB->iNumQuadFaces()
	      + pESB->iNumTriCells() + pESB->iNumQuadCells()
	      + pESB->iNumTetCells() + pESB->iNumPyrCells()
	      + pESB->iNumPrismCells() + pESB->iNumHexCells();
	  break;
	default:
	  *err = iBase_INVALID_ENTITY_TYPE;
	  instance->pGB->setError(*err);
	  *num_type = 0;
	  break;
	}
    }
  }

  void iMesh_getNumOfTopo(iMesh_Instance instance,
  /*in*/const iBase_EntitySetHandle entity_set_handle,
			  /*in*/const int entity_topology, int *num_topo,
			  int *err)
  {
    // Base public
    *err = iBase_SUCCESS;
    CHECK_ENTITY_SET(instance, entity_set_handle, err);

    if (isRootSet(instance, entity_set_handle)) {
      // The case for the entire mesh
      switch (entity_topology)
	{
	case iMesh_POINT:
	  *num_topo = instance->pGB->getNumVertices();
	  break;
	case iMesh_LINE_SEGMENT:
	  *num_topo = instance->pGB->getNumLineSegments();
	  break;
	case iMesh_TRIANGLE:
	  *num_topo = instance->pGB->getNumTriangles();
	  break;
	case iMesh_QUADRILATERAL:
	  *num_topo = instance->pGB->getNumQuadrilaterals();
	  break;
	case iMesh_TETRAHEDRON:
	  *num_topo = instance->pGB->getNumTetrahedra();
	  break;
	case iMesh_PYRAMID:
	  *num_topo = instance->pGB->getNumPyramids();
	  break;
	case iMesh_PRISM:
	  *num_topo = instance->pGB->getNumPrisms();
	  break;
	case iMesh_HEXAHEDRON:
	  *num_topo = instance->pGB->getNumHexahedra();
	  break;
	case iMesh_POLYGON:
	case iMesh_POLYHEDRON:
	case iMesh_SEPTAHEDRON:
	  *num_topo = 0;
	  break;
	case iMesh_ALL_TOPOLOGIES:
	  *num_topo = (instance->pGB->getNumVertices()
	      + instance->pGB->getNumLineSegments()
	      + instance->pGB->getNumTriangles()
	      + instance->pGB->getNumQuadrilaterals()
	      + instance->pGB->getNumTetrahedra()
	      + instance->pGB->getNumPyramids() + instance->pGB->getNumPrisms()
	      + instance->pGB->getNumHexahedra());
	  break;
	default:
	  *err = iBase_INVALID_ENTITY_TOPOLOGY;
	  instance->pGB->setError(*err);
	  *num_topo = 0;
	  break;
	}
    }
    else {
      // Have an actual entity set argument
      EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set_handle);
      switch (entity_topology)
	{
	case iMesh_POINT:
	  *num_topo = pESB->iNumVerts();
	  break;
	case iMesh_LINE_SEGMENT:
	  *num_topo = pESB->iNumEdges();
	  break;
	case iMesh_TRIANGLE:
	  // This is fairly cheesey, but works for both 2D and 3D, and the
	  // cost should be pretty small.
	  *num_topo = pESB->iNumTriFaces() + pESB->iNumTriCells();
	  break;
	case iMesh_QUADRILATERAL:
	  // This is fairly cheesey, but works for both 2D and 3D, and the
	  // cost should be pretty small.
	  *num_topo = pESB->iNumQuadFaces() + pESB->iNumQuadCells();
	  break;
	case iMesh_TETRAHEDRON:
	  *num_topo = pESB->iNumTetCells();
	  break;
	case iMesh_PYRAMID:
	  *num_topo = pESB->iNumPyrCells();
	  break;
	case iMesh_PRISM:
	  *num_topo = pESB->iNumPrismCells();
	  break;
	case iMesh_HEXAHEDRON:
	  *num_topo = pESB->iNumHexCells();
	  break;
	case iMesh_ALL_TOPOLOGIES:
	  *num_topo = (pESB->iNumVerts() + pESB->iNumEdges()
	      + pESB->iNumTriFaces() + pESB->iNumTriCells()
	      + pESB->iNumQuadFaces() + pESB->iNumQuadCells()
	      + pESB->iNumTetCells() + pESB->iNumPyrCells()
	      + pESB->iNumPrismCells() + pESB->iNumHexCells());
	  break;
	case iMesh_POLYGON:
	case iMesh_POLYHEDRON:
	case iMesh_SEPTAHEDRON:
	  *num_topo = 0;
	  break;
	default:
	  *err = iBase_INVALID_ENTITY_TOPOLOGY;
	  instance->pGB->setError(*err);
	  *num_topo = 0;
	  break;
	}
    }
  }

  void iMesh_optimize(iMesh_Instance instance, int *areHandlesInvalidated,
		      int *err)
  {
    // Base public
    *err = iBase_SUCCESS;
    instance->pGB->pMesh()->purgeAllEntities();
    *areHandlesInvalidated = true;
  }

  void iMesh_getEntities(iMesh_Instance instance,
			 const iBase_EntitySetHandle entity_set_handle,
			 const int entity_type, const int entity_topology,
			 iBase_EntityHandle** entity_handles,
			 int* entity_handles_allocated,
			 int* entity_handles_size, int *err)

  {
    *err = iBase_SUCCESS;
    // Zero this out in case there's an error.  This way the caller will
    // be sure to know there's no data, even if they allocated the array.
    *entity_handles_size = 0;
    EntitySetBase* pESB = static_cast<EntitySetBase*>(entity_set_handle);
    bool isRoot = isRootSet(instance, entity_set_handle);
    bool isOtherSet = instance->pGB->qHasSet(pESB);
    if (!isRoot && !isOtherSet) {
      *err = iBase_INVALID_ENTITYSET_HANDLE;
      instance->pGB->setError(*err);
      return;
    }
    CHECK_TYPE_AND_TOPO(entity_type, entity_topology, err);

    if (!isRoot) {
      pESB->vGetEntities(entity_type, entity_topology, entity_handles,
			 entity_handles_allocated, entity_handles_size);
      if (*err != iBase_SUCCESS)
	return;
    }
    else {
      instance->pGB->getEntities(entity_type, entity_topology, entity_handles,
				 entity_handles_allocated, entity_handles_size,
				 err);
    } // Done dealing with base mesh object
  }

  void iMesh_getVtxArrCoords(iMesh_Instance instance,
  /*in*/const iBase_EntityHandle* vertex_handles,
			     /*in*/const int vertex_handles_size,
			     /*in*/int storage_order,
			     /*inout*/double** coords,
			     /*inout*/int* coords_allocated,
			     /*out*/int* coords_size, int *err)
  {
    *err = iBase_SUCCESS;

    for (int i = 0; i < vertex_handles_size; i++) {
      CHECK_ENTITY(instance, vertex_handles[i], err);
    }

    // Start by figuring out how much data we'll need to store.
    int iGeomDim = instance->pGB->getGeometricDimension();
    *coords_size = vertex_handles_size * iGeomDim;
    TRY_ARRAY_SIZE2(coords, coords_allocated, *coords_size, double);

    // Now actually copy the data.
    if (storage_order == iBase_BLOCKED) {
      int iOut0 = 0;
      int iOut1 = iOut0 + vertex_handles_size;
      int iOut2 = iOut1 + vertex_handles_size;
      for (int iVert = 0; iVert < vertex_handles_size; iVert++) {
	Entity *pE = dynamic_cast<Entity*>(vertex_handles[iVert]);
	if (pE->getEntType() != iBase_VERTEX) {
	  *err = iBase_INVALID_ENTITY_TYPE;
	  instance->pGB->setError(*err);
	  return;
	}
	const double * const adCoords = pE->getVert(0)->getCoords();
	(*coords)[iOut0++] = adCoords[0];
	(*coords)[iOut1++] = adCoords[1];
	if (iGeomDim == 3)
	  (*coords)[iOut2++] = adCoords[2];
      }
    }
    else {
      int iOut = 0;
      for (int iVert = 0; iVert < vertex_handles_size; iVert++) {
	Entity *pE = static_cast<Entity*>(vertex_handles[iVert]);
	if (pE->getEntType() != iBase_VERTEX) {
	  *err = iBase_INVALID_ENTITY_TYPE;
	  instance->pGB->setError(*err);
	  return;
	}
	const double * const adCoords = pE->getVert(0)->getCoords();
	(*coords)[iOut++] = adCoords[0];
	(*coords)[iOut++] = adCoords[1];
	if (iGeomDim == 3)
	  (*coords)[iOut++] = adCoords[2];
      }
    }
  }

  void iMesh_initEntArrIter(iMesh_Instance instance,
  /*in*/const iBase_EntitySetHandle entity_set_handle,
			    /*in*/const int requested_entity_type,
			    /*in*/const int requested_entity_topology,
			    /*in*/const int requested_array_size,
			    /*in*/const int /*resilient*/,
			    /*out*/iBase_EntityArrIterator* entArr_iterator,
			    int *err)
  {
    // GRUMMP's iterators are (nearly) resilient now, and there's no support
    // for treating non-resilient operators differently, so ignore the arg.
    *err = iBase_SUCCESS;
    CHECK_ENTITY_SET(instance, entity_set_handle, err);
    CHECK_TYPE_AND_TOPO(requested_entity_type, requested_entity_topology, err);
    Workset* pWSNew = NULL;

    try {
      if (isRootSet(instance, entity_set_handle)) {
	pWSNew = new MeshWorkset(instance->pGB->pMesh(), requested_entity_type,
				 requested_entity_topology,
				 requested_array_size);
      }
      else {
	EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set_handle);
	if (pESB->qIsOrdered()) {
	  pWSNew = new EntityListWorkset(dynamic_cast<EntityList*>(pESB),
					 requested_entity_type,
					 requested_entity_topology,
					 requested_array_size);
	}
	else {
	  pWSNew = new EntitySetWorkset(dynamic_cast<EntitySet*>(pESB),
					requested_entity_type,
					requested_entity_topology,
					requested_array_size);
	}
      }
    }
    catch (iBase_ErrorType& eErr) {
      *err = eErr;
      instance->pGB->setError(*err);
      if (pWSNew)
	delete pWSNew;
      return;
    }

    assert(pWSNew);
    instance->pGB->addWorkset(pWSNew);
    *entArr_iterator = pWSNew;
  }

  void iMesh_getNextEntArrIter(iMesh_Instance instance,
  /*in*/iBase_EntityArrIterator entArr_iterator,
			       /*inout*/iBase_EntityHandle** entity_handles,
			       /*inout*/int* entity_handles_allocated,
			       /*out*/int* entity_handles_size, int *data_left,
			       int *err)
  {
    *err = iBase_SUCCESS;
    Workset * pWS = static_cast<Workset*>(entArr_iterator);
    CHECK_ENTARR_ITER(instance, pWS, err);
    TRY(*data_left = pWS->qGetNextBlock(entity_handles,
					entity_handles_allocated,
					entity_handles_size)
    ;
  );
}

void iMesh_resetEntArrIter(iMesh_Instance instance,
/*in*/iBase_EntityArrIterator entArr_iterator,
			   int *err)
{
  *err = iBase_SUCCESS;
  Workset * pWS = static_cast<Workset*>(entArr_iterator);
  CHECK_ENTARR_ITER(instance, pWS, err);
  TRY(pWS->qReset()
  ;
);
}

void iMesh_endEntArrIter(iMesh_Instance instance,
/*in*/iBase_EntityArrIterator entArr_iterator,
			 int *err)
{
*err = iBase_SUCCESS;
Workset * pWS = static_cast<Workset*>(entArr_iterator);
CHECK_ENTARR_ITER(instance, pWS, err);
instance->pGB->removeWorkset(pWS, err);
delete pWS;
}

void iMesh_getEntArrTopo(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*inout*/int** topology,
			 /*inout*/int* topology_allocated,
			 /*out*/int* topology_size,
			 int *err)
{
*err = iBase_SUCCESS;
#ifndef NDEBUG
for (int i = 0; i < entity_handles_size; i++) {
  CHECK_ENTITY(instance, entity_handles[i], err);
}
#endif

TRY_ARRAY_SIZE2(topology, topology_allocated, entity_handles_size, int);
for (int i = 0; i < entity_handles_size; i++) {
  Entity *pE = static_cast<Entity*>(entity_handles[i]);
  (*topology)[i] = pE->getEntTopology();
}
*topology_size = entity_handles_size;
}

void iMesh_getEntArrType(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*inout*/int** type,
			 /*inout*/int* type_allocated,
			 /*out*/int* type_size,
			 int *err)
{
*err = iBase_SUCCESS;
#ifndef NDEBUG
for (int i = 0; i < entity_handles_size; i++) {
  CHECK_ENTITY(instance, entity_handles[i], err);
}
#endif

TRY_ARRAY_SIZE2(type, type_allocated, entity_handles_size, int);
for (int i = 0; i < entity_handles_size; i++) {
  Entity *pE = static_cast<Entity*>(entity_handles[i]);
  (*type)[i] = pE->getEntType();
}
*type_size = entity_handles_size;
}

void iMesh_getEntArrAdj(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			/*in*/const int entity_handles_size,
			/*in*/const int entity_type_requested,
			/*inout*/iBase_EntityHandle** adj_entity_handles,
			/*inout*/int* adj_entity_handles_allocated,
			/*out*/int* adj_entity_handles_size,
			/*inout*/int** offset,
			/*inout*/int* offset_allocated,
			/*out*/int* offset_size,
			int *err)
{
*err = iBase_SUCCESS;

#ifndef NDEBUG
for (int i = 0; i < entity_handles_size; i++) {
  CHECK_ENTITY(instance, entity_handles[i], err);
}
#endif

instance->pGB->getEntArrAdj(entity_handles, entity_handles_size,
			    entity_type_requested, adj_entity_handles,
			    adj_entity_handles_allocated,
			    adj_entity_handles_size, offset, offset_allocated,
			    offset_size, err);
} // End of iMesh_getEntArrAdj

void iMesh_createEntSet(iMesh_Instance instance,
/*in*/const int isList,
			/*out*/iBase_EntitySetHandle* entity_set_created,
			int *err)
{
*err = iBase_SUCCESS;

EntitySetBase* pESB;
int iDim;
iMesh_getGeometricDimension(instance, &iDim, err);
if (iDim == 3) {
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(instance, &root_set, err);
  int nRegions = 0;
  iMesh_getNumOfType(instance, root_set, iBase_REGION, &nRegions, err);
  if (nRegions == 0)
    iDim = 2;
}

try {
  if (isList) {
    pESB = new EntityList(iDim);
  }
  else {
    pESB = new EntitySet(iDim);
  }
}
catch (...) {
// Hard to imagine anything but a memory allocation error...
  *err = iBase_MEMORY_ALLOCATION_FAILED;
  instance->pGB->setError(*err);
}
assert(pESB);
instance->pGB->addEntSet(pESB);
*entity_set_created = pESB;
}

void iMesh_destroyEntSet(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set,
			 int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(instance->pGB->destroyEntSet(entity_set)
;
);
}

void iMesh_isList(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
		  int *is_list, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
iBase_EntitySetHandle root_set;
iMesh_getRootSet(instance, &root_set, err);
if (*err != iBase_SUCCESS)
return;
if (root_set == entity_set)
*is_list = 0;
else
*is_list = instance->pGB->isList(entity_set);
}

void iMesh_getNumEntSets(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_handle,
			 /*in*/const int num_hops, int *num_sets, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
TRY(*num_sets = instance->pGB->getNumEntSets(entity_set_handle, num_hops)
;
);
}

void iMesh_getEntSets(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_handle,
		      /*in*/const int num_hops,
		      /*out*/iBase_EntitySetHandle** contained_set_handles,
		      /*out*/int* contained_set_handles_allocated,
		      /*out*/int* contained_set_handles_size,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
TRY(
instance->pGB->getEntSets(entity_set_handle, num_hops, contained_set_handles,
			  contained_set_handles_allocated,
			  contained_set_handles_size)
;
);
}

void iMesh_addEntToSet(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		       /*inout*/iBase_EntitySetHandle entity_set, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(instance->pGB->addEntToSet(entity_handle, entity_set)
;
);
}

void iMesh_rmvEntFromSet(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
			 /*inout*/iBase_EntitySetHandle entity_set, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(instance->pGB->rmvEntFromSet(entity_handle, entity_set)
;
);
}

void iMesh_addEntArrToSet(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			  /*in*/const int entity_handles_size,
			  /*inout*/iBase_EntitySetHandle entity_set,
			  int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(
instance->pGB->addEntArrToSet(entity_handles, entity_handles_size, entity_set)
;
);
}

void iMesh_rmvEntArrFromSet(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			    /*in*/const int entity_handles_size,
			    /*inout*/iBase_EntitySetHandle entity_set,
			    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(
instance->pGB->rmvEntArrFromSet(entity_handles, entity_handles_size, entity_set)
;
);
}

void iMesh_addEntSet(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set_to_add,
		     /*inout*/iBase_EntitySetHandle entity_set_handle, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_to_add, err);
CHECK_ENTITY_SET(instance, entity_set_handle, err);
TRY(instance->pGB->addEntSet(entity_set_to_add, entity_set_handle)
;
);
}

void iMesh_rmvEntSet(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_to_remove,
		     /*inout*/iBase_EntitySetHandle entity_set_handle, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_to_remove, err);
CHECK_ENTITY_SET(instance, entity_set_handle, err);
TRY(instance->pGB->rmvEntSet(entity_set_to_remove, entity_set_handle)
;
);
}

void iMesh_isEntContained(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle containing_entity_set,
			  /*in*/const iBase_EntityHandle contained_entity,
			  int *is_contained, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, containing_entity_set, err);
CHECK_ENTITY(instance, contained_entity, err);
TRY(
*is_contained = instance->pGB->isEntContained(containing_entity_set,
					      contained_entity)
;
);
}

void iMesh_isEntArrContained(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle containing_set,
			     /*in*/const iBase_EntityHandle* entity_handles,
			     /*in*/int num_entity_handles,
			     /*inout*/int** is_contained,
			     /*inout*/int* is_contained_allocated,
			     /*out*/int* is_contained_size,
			     /*out*/int* err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, containing_set, err);
TRY_ARRAY_SIZE2(is_contained, is_contained_allocated, num_entity_handles, int);
*is_contained_size = num_entity_handles;
for (int i = 0; i < num_entity_handles; i++) {
CHECK_ENTITY(instance, entity_handles[i], err);
TRY(
(*is_contained)[i] = instance->pGB->isEntContained(containing_set,
						   entity_handles[i])
;
);
}
}

void iMesh_isEntSetContained(
iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle containing_entity_set,
/*in*/const iBase_EntitySetHandle contained_entity_set, int *is_contained,
int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, containing_entity_set, err);
CHECK_ENTITY_SET(instance, contained_entity_set, err);
TRY(
*is_contained = instance->pGB->isEntSetContained(containing_entity_set,
						 contained_entity_set)
;
);
}

void iMesh_addPrntChld(iMesh_Instance instance,
/*inout*/iBase_EntitySetHandle parent_entity_set,
		       /*inout*/iBase_EntitySetHandle child_entity_set,
		       int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, parent_entity_set, err);
CHECK_ENTITY_SET(instance, child_entity_set, err);
TRY(instance->pGB->addPrntChld(parent_entity_set, child_entity_set)
;
);
}

void iMesh_rmvPrntChld(iMesh_Instance instance,
/*inout*/iBase_EntitySetHandle parent_entity_set,
		       /*inout*/iBase_EntitySetHandle child_entity_set,
		       int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, parent_entity_set, err);
CHECK_ENTITY_SET(instance, child_entity_set, err);
TRY(instance->pGB->rmvPrntChld(parent_entity_set, child_entity_set)
;
);
}

void iMesh_isChildOf(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle parent_entity_set,
		     /*in*/const iBase_EntitySetHandle child_entity_set,
		     int *is_child, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, parent_entity_set, err);
CHECK_ENTITY_SET(instance, child_entity_set, err);
TRY(*is_child = instance->pGB->isChildOf(parent_entity_set, child_entity_set)
;
);
}

void iMesh_getNumChld(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
		      /*in*/const int num_hops, int *num_child, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(*num_child = instance->pGB->getNumChld(entity_set, num_hops)
;
);
}

void iMesh_getNumPrnt(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
		      /*in*/const int num_hops, int *num_parent, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
TRY(*num_parent = instance->pGB->getNumPrnt(entity_set, num_hops)
;
);
}

void iMesh_getChldn(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle from_entity_set,
		    /*in*/const int num_hops,
		    /*out*/iBase_EntitySetHandle** entity_set_handles,
		    /*out*/int* entity_set_handles_allocated,
		    /*out*/int* entity_set_handles_size,
		    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, from_entity_set, err);
TRY(
instance->pGB->getChldn(from_entity_set, num_hops, entity_set_handles,
			entity_set_handles_allocated, entity_set_handles_size)
;
);
}

void iMesh_getPrnts(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle from_entity_set,
		    /*in*/const int num_hops,
		    /*out*/iBase_EntitySetHandle** entity_set_handles,
		    /*out*/int* entity_set_handles_allocated,
		    /*out*/int* entity_set_handles_size,
		    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, from_entity_set, err);
TRY(
instance->pGB->getPrnts(from_entity_set, num_hops, entity_set_handles,
			entity_set_handles_allocated, entity_set_handles_size)
;
);
}

void iMesh_setVtxArrCoords(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* vertex_handles,
			   /*in*/const int vertex_handles_size,
			   /*in*/const int storage_order,
			   /*in*/const double* new_coords,
			   /*in*/const int new_coords_size,
			   int *err)
{
*err = iBase_SUCCESS;

#ifndef NDEBUG
for (int i = 0; i < vertex_handles_size; i++) {
CHECK_ENTITY(instance, vertex_handles[i], err);
}
#endif
int iGeomDim = instance->pGB->getGeometricDimension();
if ((storage_order != iBase_BLOCKED && storage_order != iBase_INTERLEAVED)
|| new_coords_size != iGeomDim * vertex_handles_size) {
*err = iBase_INVALID_ARGUMENT;
instance->pGB->setError(*err);
return;
}

int iBase = 0;
for (int i = 0; i < vertex_handles_size; i++) {
double adNewCoords[] =
{ 0, 0, 0 };
if (storage_order == iBase_BLOCKED) {
adNewCoords[0] = new_coords[iBase];
adNewCoords[1] = new_coords[iBase + vertex_handles_size];
if (iGeomDim == 3)
adNewCoords[2] = new_coords[iBase + 2 * vertex_handles_size];
iBase++;
}
else {
adNewCoords[0] = new_coords[iBase];
adNewCoords[1] = new_coords[iBase + 1];
if (iGeomDim == 3)
adNewCoords[2] = new_coords[iBase + 2];
iBase += iGeomDim;
}

Entity *pEVertex = static_cast<Entity*>(vertex_handles[i]);
Vert *pV = dynamic_cast<Vert*>(pEVertex);
if (pV == pVInvalidVert) {
*err = iBase_INVALID_ENTITY_HANDLE;
instance->pGB->setError(*err);
continue;
}
pV->setCoords(iGeomDim, adNewCoords);
} // Done with loop over verts
}

void iMesh_createVtxArr(iMesh_Instance instance,
/*in*/const int num_verts,
			/*in*/const int storage_order,
			/*in*/const double* new_coords,
			/*in*/const int new_coords_size,
			/*inout*/iBase_EntityHandle** new_vertex_handles,
			/*inout*/int* new_vertex_handles_allocated,
			/*inout*/int* new_vertex_handles_size,
			int *err)
{
*err = iBase_SUCCESS;

if (instance->pGB == NULL || instance->pGB->pMesh() == NULL) {
if (new_coords_size == 3 * num_verts) {
instance->vInitiMesh(iMesh_Base::eVolume);
}
else {
assert(new_coords_size == 2 * num_verts);
instance->vInitiMesh(iMesh_Base::ePlane);
}
 // There's just no way to detect a surface mesh here.
}
 // If you just initialized, of course these will be okay.
int iGeomDim = instance->pGB->getGeometricDimension();
assert(new_coords_size == iGeomDim * num_verts);

TRY_ARRAY_SIZE2(new_vertex_handles, new_vertex_handles_allocated, num_verts,
		iBase_EntityHandle);
*new_vertex_handles_size = num_verts;

if (instance->pGB == NULL || instance->pGB->pMesh() == NULL) {
instance->vInitiMesh(iMesh_Base::eMesh_Type(iGeomDim));
}

int iBase = 0;
for (int i = 0; i < num_verts; i++) {
double adCoords[] =
{ 0, 0, 0 };
switch (storage_order)
{
case iBase_INTERLEAVED:
{
adCoords[0] = new_coords[iBase];
adCoords[1] = new_coords[iBase + 1];
if (iGeomDim == 3)
adCoords[2] = new_coords[iBase + 2];
iBase += iGeomDim;
}
break;
case iBase_BLOCKED:
{
adCoords[0] = new_coords[iBase];
adCoords[1] = new_coords[iBase + num_verts];
if (iGeomDim == 3)
adCoords[2] = new_coords[iBase + 2 * num_verts];
iBase++;
}
break;
default:
*err = iBase_INVALID_ARGUMENT;
instance->pGB->setError(*err);
break;
}
Vert *pV = instance->pGB->createVert(adCoords[0], adCoords[1], adCoords[2]);
(*new_vertex_handles)[i] = pV;
}
}

void iMesh_createEntArr(iMesh_Instance instance,
/*in*/const int new_ent_topology,
			/*in*/const iBase_EntityHandle* lower_order_ent_handles,
			/*in*/const int lower_order_ent_handles_size,
			/*out*/iBase_EntityHandle** new_ent_handles,
			/*out*/int* new_ent_handles_allocated,
			/*out*/int* new_ent_handles_size,
			/*inout*/int** status,
			/*inout*/int* status_allocated,
			/*out*/int* status_size,
			int *err)
{
*err = iBase_SUCCESS;
#ifndef NDEBUG
for (int i = 0; i < lower_order_ent_handles_size; i++) {
CHECK_ENTITY(instance, lower_order_ent_handles[i], err);
}
#endif

instance->pGB->createEntArr(new_ent_topology, lower_order_ent_handles,
			    lower_order_ent_handles_size, new_ent_handles,
			    new_ent_handles_allocated, new_ent_handles_size,
			    status, status_allocated, status_size, err);
}

void iMesh_deleteEntArr(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			/*in*/const int entity_handles_size, int *err)
{
*err = iBase_SUCCESS;

#ifndef NDEBUG
for (int i = 0; i < entity_handles_size; i++) {
CHECK_ENTITY(instance, entity_handles[i], err);
}
#endif

int iLast = entity_handles_size - 1;
for (int i = 0; i <= iLast; i++) {
iMesh_deleteEnt(instance, entity_handles[i], err);
}
}

void iMesh_createTag(iMesh_Instance instance,
/*in*/const char* tag_name_in,
		     /*in*/const int tag_size,
		     /*in*/const int tag_type,
		     /*out*/iBase_TagHandle* tag_handle, int *err,
		     /*in*/const int tag_name_len)
{
*err = iBase_SUCCESS;
if (instance->pGB == NULL || instance->pGB->pMesh() == NULL) {
instance->vInitiMesh();
}
char tag_name[1024];
NULL_TERMINATE(tag_name_in, tag_name_len, tag_name);
TRY(
instance->pGB->createTag(tag_name, tag_name_len, tag_size, tag_type, tag_handle)
;
);
}

void iMesh_destroyTag(iMesh_Instance instance,
/*in*/iBase_TagHandle tag_handle,
		      /*in*/const int forced, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->destroyTag(tag_handle, forced)
;
);
}

void iMesh_getTagName(iMesh_Instance instance,
/*in*/const iBase_TagHandle tag_handle,
		      /*out*/char *name, int *err, int name_len)
{
*err = iBase_SUCCESS;
 // No need to null terminate name here, as it's got intent out.
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->getTagName(tag_handle, name, name_len)
;
);
}

void iMesh_getTagSizeValues(iMesh_Instance instance,
/*in*/const iBase_TagHandle tag_handle,
			    int *tag_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(*tag_size = instance->pGB->getTagSizeValues(tag_handle)
;
);
}

void iMesh_getTagSizeBytes(iMesh_Instance instance,
/*in*/const iBase_TagHandle tag_handle,
			   int *tag_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(*tag_size = instance->pGB->getTagSizeBytes(tag_handle)
;
);
}

void iMesh_getTagHandle(iMesh_Instance instance,
/*in*/const char* tag_name_in,
			iBase_TagHandle *tag_handle, int *err, int tag_name_len)
{
*err = iBase_SUCCESS;
char tag_name[1024];
NULL_TERMINATE(tag_name_in, tag_name_len, tag_name);
TRY(*tag_handle = instance->pGB->getTagHandle(tag_name, tag_name_len)
;
);
}

void iMesh_getTagType(iMesh_Instance instance,
/*in*/const iBase_TagHandle tag_handle,
		      int *tag_type, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(*tag_type = instance->pGB->getTagType(tag_handle)
;
);
}

void iMesh_setEntSetData(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set_handle,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*in*/const void* tag_value,
			 /*in*/const int tag_value_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
CHECK_TAG(instance, tag_handle, err);

TRY(
instance->pGB->setEntSetData(entity_set_handle, tag_handle,
			     reinterpret_cast<const char*>(tag_value),
			     tag_value_size)
;
);
}

void iMesh_setEntSetIntData(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    /*in*/const int tag_value,
			    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->setEntSetIntData(entity_set, tag_handle, tag_value)
;
);
}

void iMesh_setEntSetDblData(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    /*in*/const double tag_value,
			    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->setEntSetDblData(entity_set, tag_handle, tag_value)
;
);
}

void iMesh_setEntSetEHData(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set,
			   /*in*/const iBase_TagHandle tag_handle,
			   /*in*/const iBase_EntityHandle tag_value,
			   int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_ENTITY(instance, tag_value, err);
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->setEntSetEHData(entity_set, tag_handle, tag_value)
;
);
}

void iMesh_setEntSetESHData(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    /*in*/const iBase_EntitySetHandle tag_value,
			    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_ENTITY_SET(instance, tag_value, err);
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->setEntSetESHData(entity_set, tag_handle, tag_value)
;
);
}

void iMesh_getEntSetData(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_handle,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*inout*/void* tag_value,
			 /*inout*/int* tag_value_allocated,
			 /*inout*/int* tag_value_size,
			 int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
CHECK_TAG(instance, tag_handle, err);

TRY(
instance->pGB->getEntSetData(entity_set_handle, tag_handle,
			     reinterpret_cast<char**>(tag_value),
			     tag_value_allocated, tag_value_size)
;
);
}

void iMesh_getEntSetIntData(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    int *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(*out_data = instance->pGB->getEntSetIntData(entity_set, tag_handle)
;
);
}

void iMesh_getEntSetDblData(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    double *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(*out_data = instance->pGB->getEntSetDblData(entity_set, tag_handle)
;
);
}

void iMesh_getEntSetEHData(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
			   /*in*/const iBase_TagHandle tag_handle,
			   iBase_EntityHandle *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(*out_data = instance->pGB->getEntSetEHData(entity_set, tag_handle)
;
);
}

void iMesh_getEntSetESHData(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set,
			    /*in*/const iBase_TagHandle tag_handle,
			    iBase_EntitySetHandle *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set, err);
CHECK_TAG(instance, tag_handle, err);

TRY(*out_data = instance->pGB->getEntSetESHData(entity_set, tag_handle)
;
);
}

void iMesh_getAllEntSetTags(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_handle,
			    /*out*/iBase_TagHandle** tag_handles,
			    /*out*/int* tag_handles_allocated,
			    /*out*/int* tag_handles_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);

TRY(
instance->pGB->getAllEntSetTags(entity_set_handle, tag_handles,
				tag_handles_allocated, tag_handles_size)
;
);
}

void iMesh_rmvEntSetTag(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set_handle,
			/*in*/const iBase_TagHandle tag_handle, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->rmvEntSetTag(entity_set_handle, tag_handle)
;
);
}

void iMesh_setVtxCoord(iMesh_Instance instance,
/*in*/iBase_EntityHandle vertex_handle,
		       /*in*/const double x, /*in*/const double y,
		       /*in*/const double z,
		       int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, vertex_handle, err);

int iDim;
iMesh_getGeometricDimension(instance, &iDim, err);
double adNewCoords[] =
{ x, y, z };
if (iDim == 2)
adNewCoords[2] = 0;

Entity *pEVertex = static_cast<Entity*>(vertex_handle);
Vert *pV = dynamic_cast<Vert*>(pEVertex);
if (pV == pVInvalidVert) {
*err = iBase_INVALID_ENTITY_HANDLE;
instance->pGB->setError(*err);
}
else
pV->setCoords(iDim, adNewCoords);
}

void iMesh_createVtx(iMesh_Instance instance,
/*in*/const double x, /*in*/
		     const double y,
		     /*in*/const double z,
		     /*out*/iBase_EntityHandle* new_vertex_handle, int *err)
{
*err = iBase_SUCCESS;

if (instance->pGB == NULL || instance->pGB->pMesh() == NULL) {
 // There's no way to assume anything useful here about mesh type.
instance->vInitiMesh(iMesh_Base::eVolume);
}

Vert *pV = instance->pGB->createVert(x, y, z);
*new_vertex_handle = pV;
}

void iMesh_createEnt(iMesh_Instance instance,
/*in*/const int new_entity_topology,
		     /*in*/const iBase_EntityHandle* lower_order_entity_handles,
		     /*in*/const int lower_order_entity_handles_size,
		     /*out*/iBase_EntityHandle* new_entity_handle,
		     /*out*/int* status,
		     int *err)
{
*err = iBase_SUCCESS;
#ifndef NDEBUG
for (int i = 0; i < lower_order_entity_handles_size; i++) {
CHECK_ENTITY(instance, lower_order_entity_handles[i], err);
}
#endif

instance->pGB->createEnt(new_entity_topology, lower_order_entity_handles,
			 lower_order_entity_handles_size, new_entity_handle,
			 status, err);
}

void iMesh_deleteEnt(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		     int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);

instance->pGB->deleteEnt(entity_handle, err);
}

void iMesh_getArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
		      /*in*/const int entity_handles_size,
		      /*in*/const iBase_TagHandle tag_handle,
		      /*inout*/void* tag_values,
		      /*inout*/int* tag_values_allocated,
		      /*out*/int* tag_values_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getArrData(entity_handles, entity_handles_size, tag_handle,
			  reinterpret_cast<char**>(tag_values),
			  tag_values_allocated, tag_values_size)
;
);
}

void iMesh_getIntArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*inout*/int** tag_values,
			 /*inout*/int* tag_values_allocated,
			 /*out*/int* tag_values_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getIntArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_allocated, tag_values_size)
;
);
}

void iMesh_getDblArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*inout*/double** tag_values,
			 /*inout*/int* tag_values_allocated,
			 /*out*/int* tag_values_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getDblArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_allocated, tag_values_size)
;
);
}

void iMesh_getEHArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			/*in*/const int entity_handles_size,
			/*in*/const iBase_TagHandle tag_handle,
			/*inout*/iBase_EntityHandle** tag_values,
			/*inout*/int* tag_values_allocated,
			/*out*/int* tag_values_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getEHArrData(entity_handles, entity_handles_size, tag_handle,
			    tag_values, tag_values_allocated, tag_values_size)
;
);
}

void iMesh_getESHArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*inout*/iBase_EntitySetHandle** tag_values,
			 /*inout*/int* tag_values_allocated,
			 /*out*/int* tag_values_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getESHArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_allocated, tag_values_size)
;
);
}

void iMesh_setArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
		      /*in*/const int entity_handles_size,
		      /*in*/const iBase_TagHandle tag_handle,
		      /*in*/const void* tag_values,
		      /*in*/const int tag_values_size,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->setArrData(entity_handles, entity_handles_size, tag_handle,
			  reinterpret_cast<const char*>(tag_values),
			  tag_values_size)
;
);
}

void iMesh_setIntArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*in*/const int* tag_values,
			 /*in*/const int tag_values_size,
			 int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);

TRY(
instance->pGB->setIntArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_size)
;
);
}

void iMesh_setDblArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*in*/const double* tag_values,
			 /*in*/const int tag_values_size,
			 int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->setDblArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_size)
;
);
}

void iMesh_setEHArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			/*in*/const int entity_handles_size,
			/*in*/const iBase_TagHandle tag_handle,
			/*in*/const iBase_EntityHandle* tag_values,
			/*in*/const int tag_values_size,
			int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);

TRY(
instance->pGB->setEHArrData(entity_handles, entity_handles_size, tag_handle,
			    tag_values, tag_values_size)
;
);
}

void iMesh_setESHArrData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
			 /*in*/const int entity_handles_size,
			 /*in*/const iBase_TagHandle tag_handle,
			 /*in*/const iBase_EntitySetHandle* tag_values,
			 /*in*/const int tag_values_size,
			 int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);

TRY(
instance->pGB->setESHArrData(entity_handles, entity_handles_size, tag_handle,
			     tag_values, tag_values_size)
;
);
}

void iMesh_rmvArrTag(iMesh_Instance instance,
/*in*/const iBase_EntityHandle* entity_handles,
		     /*in*/const int entity_handles_size,
		     /*in*/const iBase_TagHandle tag_handle,
		     int *err)
{
*err = iBase_SUCCESS;
CHECK_TAG(instance, tag_handle, err);

TRY(instance->pGB->rmvArrTag(entity_handles, entity_handles_size, tag_handle)
;
);
}

void iMesh_getData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		   /*in*/const iBase_TagHandle tag_handle,
		   /*inout*/void* tag_value,
		   /*inout*/int *tag_value_allocated,
		   /*out*/int *tag_value_size,
		   int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->getData(entity_handle, tag_handle,
		       reinterpret_cast<char**>(tag_value), tag_value_allocated,
		       tag_value_size)
;
);
}

void iMesh_getIntData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle, int *out_data,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(*out_data = instance->pGB->getIntData(entity_handle, tag_handle)
;
);
}

void iMesh_getDblData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle, double *out_data,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(*out_data = instance->pGB->getDblData(entity_handle, tag_handle)
;
);
}

void iMesh_getEHData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		     /*in*/const iBase_TagHandle tag_handle,
		     iBase_EntityHandle *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(*out_data = instance->pGB->getEHData(entity_handle, tag_handle)
;
);
}

void iMesh_getESHData(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle,
		      iBase_EntitySetHandle *out_data, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(*out_data = instance->pGB->getESHData(entity_handle, tag_handle)
;
);
}

void iMesh_setData(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		   /*in*/const iBase_TagHandle tag_handle,
		   /*in*/const void* tag_value,
		   /*in*/const int tag_value_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(
instance->pGB->setData(entity_handle, tag_handle,
		       reinterpret_cast<const char*>(tag_value), tag_value_size)
;
);
}

void iMesh_setIntData(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle,
		      /*in*/const int tag_value,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->setIntData(entity_handle, tag_handle, tag_value)
;
);
}

void iMesh_setDblData(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle,
		      /*in*/const double tag_value,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->setDblData(entity_handle, tag_handle, tag_value)
;
);
}

void iMesh_setEHData(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		     /*in*/const iBase_TagHandle tag_handle,
		     /*in*/const iBase_EntityHandle tag_value,
		     int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_ENTITY(instance, tag_value, err);
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->setEHData(entity_handle, tag_handle, tag_value)
;
);
}

void iMesh_setESHData(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		      /*in*/const iBase_TagHandle tag_handle,
		      /*in*/const iBase_EntitySetHandle tag_value,
		      int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_ENTITY(instance, tag_value, err);
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->setESHData(entity_handle, tag_handle, tag_value)
;
);
}

void iMesh_getAllTags(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      /*inout*/iBase_TagHandle** tag_handles,
		      /*inout*/int* tag_handles_allocated,
		      /*out*/int* tag_handles_size, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
TRY(
instance->pGB->getAllTags(entity_handle, tag_handles, tag_handles_allocated,
			  tag_handles_size)
;
);
}

void iMesh_rmvTag(iMesh_Instance instance,
/*in*/iBase_EntityHandle entity_handle,
		  /*in*/const iBase_TagHandle tag_handle, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
CHECK_TAG(instance, tag_handle, err);
TRY(instance->pGB->rmvTag(entity_handle, tag_handle)
;
);
}

void iMesh_initEntIter(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_handle,
		       /*in*/const int requested_entity_type,
		       /*in*/const int requested_entity_topology,
		       /*in*/const int /*resilient*/,
		       /*out*/iBase_EntityIterator* entity_iterator,
		       int *err)
{
 // GRUMMP's iterators are (nearly) resilient now, and there's no support
 // for treating non-resilient operators differently, so ignore the arg.
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_handle, err);
CHECK_TYPE_AND_TOPO(requested_entity_type, requested_entity_topology, err);
EntIter* pEINew = NULL;

try {
if (isRootSet(instance, entity_set_handle)) {
pEINew = new MeshIter(instance->pGB->pMesh(), requested_entity_type,
		      requested_entity_topology);
}
else {
EntitySetBase *pESB = static_cast<EntitySetBase*>(entity_set_handle);
if (pESB->qIsOrdered()) {
assert(dynamic_cast<EntityList*>(pESB));
pEINew = new EntityListIter(dynamic_cast<EntityList*>(pESB),
			    requested_entity_type, requested_entity_topology);
}
else {
assert(dynamic_cast<EntitySet*>(pESB));
pEINew = new EntitySetIter(dynamic_cast<EntitySet*>(pESB),
			   requested_entity_type, requested_entity_topology);
}
}
}
catch (iBase_ErrorType& eErr) {
*err = eErr;
instance->pGB->setError(*err);
if (pEINew)
delete pEINew;
return;
}
instance->pGB->addEntIter(pEINew);
*entity_iterator = pEINew;
}

void iMesh_getNextEntIter(iMesh_Instance instance,
/*in*/iBase_EntityIterator entity_iterator,
			  /*out*/iBase_EntityHandle* entity_handle,
			  int *has_data, int *err)
{
*err = iBase_SUCCESS;
EntIter * pEI = static_cast<EntIter*>(entity_iterator);
CHECK_ENT_ITER(instance, pEI, err);
bool qRetVal = pEI->qIncrement();
*entity_handle = pEI->pvData();
*has_data = qRetVal;
}

void iMesh_resetEntIter(iMesh_Instance instance,
/*in*/iBase_EntityIterator entity_iterator,
			int *err)
{
*err = iBase_SUCCESS;
EntIter * pEI = static_cast<EntIter*>(entity_iterator);
CHECK_ENT_ITER(instance, pEI, err);
TRY(pEI->qReset());
}

void iMesh_endEntIter(iMesh_Instance instance,
/*in*/iBase_EntityIterator entity_iterator,
		      int *err)
{
*err = iBase_SUCCESS;
EntIter * pEI = static_cast<EntIter*>(entity_iterator);
CHECK_ENT_ITER(instance, pEI, err);
instance->pGB->removeEntIter(pEI, err);
delete pEI;
}

void iMesh_getEntTopo(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      int *out_topo, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
Entity *pE = static_cast<Entity*>(entity_handle);
*out_topo = pE->getEntTopology();
}

void iMesh_getEntType(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		      int *out_type, int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);
Entity *pE = static_cast<Entity*>(entity_handle);
*out_type = pE->getEntType();
}

void iMesh_getVtxCoord(iMesh_Instance instance,
/*in*/const iBase_EntityHandle vertex_handle,
		       /*out*/double *x, /*out*/double *y, /*out*/double *z,
		       int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, vertex_handle, err);

Entity *pE = static_cast<Entity*>(vertex_handle);
if (pE->getEntType() != iBase_VERTEX) {
*err = iBase_INVALID_ENTITY_TYPE;
instance->pGB->setError(*err);
return;
}
const double * const adCoords = pE->getVert(0)->getCoords();
*x = adCoords[0];
*y = adCoords[1];
*z = adCoords[2];
}

void iMesh_getEntAdj(iMesh_Instance instance,
/*in*/const iBase_EntityHandle entity_handle,
		     /*in*/const int entity_type_requested,
		     /*inout*/iBase_EntityHandle** adj_entity_handles,
		     /*inout*/int* adj_entity_handles_allocated,
		     /*out*/int* adj_entity_handles_size,
		     int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);

instance->pGB->getEntAdj(entity_handle, entity_type_requested,
			 adj_entity_handles, adj_entity_handles_allocated,
			 adj_entity_handles_size, err);
}

void iMesh_getEnt2ndAdj(iMesh_Instance instance,
			iBase_EntityHandle entity_handle,
			int bridge_entity_type, int requested_entity_type,
			iBase_EntityHandle** adj_entity_handles,
			int* adj_entity_handles_allocated,
			int* adj_entity_handles_size, int* err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY(instance, entity_handle, err);

instance->pGB->getEnt2ndAdj(entity_handle, bridge_entity_type,
			    requested_entity_type, adj_entity_handles,
			    adj_entity_handles_allocated,
			    adj_entity_handles_size, err);
}

void iMesh_getEntArr2ndAdj(iMesh_Instance instance,
			   iBase_EntityHandle const* entity_handles,
			   int entity_handles_size, int bridge_entity_type,
			   int requested_entity_type,
			   iBase_EntityHandle** adj_entity_handles,
			   int* adj_entity_handles_allocated,
			   int* adj_entity_handles_size, int** offset,
			   int* offset_allocated, int* offset_size, int* err)
{
*err = iBase_SUCCESS;
#ifndef NDEBUG
for (int i = 0; i < entity_handles_size; i++) {
CHECK_ENTITY(instance, entity_handles[i], err);
}
#endif

instance->pGB->getEntArr2ndAdj(entity_handles, entity_handles_size,
			       bridge_entity_type, requested_entity_type,
			       adj_entity_handles, adj_entity_handles_allocated,
			       adj_entity_handles_size, offset,
			       offset_allocated, offset_size, err);
}

void iMesh_subtract(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_1,
		    /*in*/const iBase_EntitySetHandle entity_set_2,
		    /*out*/iBase_EntitySetHandle* result_entity_set,
		    int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_1, err);
CHECK_ENTITY_SET(instance, entity_set_2, err);
TRY(instance->pGB->subtract(entity_set_1, entity_set_2, result_entity_set)
;
);
}

void iMesh_intersect(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_1,
		     /*in*/const iBase_EntitySetHandle entity_set_2,
		     /*out*/iBase_EntitySetHandle* result_entity_set,
		     int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_1, err);
CHECK_ENTITY_SET(instance, entity_set_2, err);
TRY(instance->pGB->intersect(entity_set_1, entity_set_2, result_entity_set)
;
);
}

void iMesh_unite(iMesh_Instance instance,
/*in*/const iBase_EntitySetHandle entity_set_1,
		 /*in*/const iBase_EntitySetHandle entity_set_2,
		 /*out*/iBase_EntitySetHandle* result_entity_set,
		 int *err)
{
*err = iBase_SUCCESS;
CHECK_ENTITY_SET(instance, entity_set_1, err);
CHECK_ENTITY_SET(instance, entity_set_2, err);
TRY(instance->pGB->unite(entity_set_1, entity_set_2, result_entity_set)
;
);
}

 // The following implementation is courtesy of Jason Kraftcheck.

void iMesh_getAdjEntIndices(iMesh_Instance instance,
/*in*/iBase_EntitySetHandle entity_set_handle,
			    /*in*/int entity_type_requestor,
			    /*in*/int entity_topology_requestor,
			    /*in*/int entity_type_requested,
			    /*inout*/iBase_EntityHandle** entity_handles,
			    /*inout*/int* entity_handles_allocated,
			    /*out*/int* entity_handles_size,
			    /*inout*/iBase_EntityHandle** adj_entity_handles,
			    /*inout*/int* adj_entity_handles_allocated,
			    /*out*/int* adj_entity_handles_size,
			    /*inout*/int** adj_entity_indices,
			    /*inout*/int* adj_entity_indices_allocated,
			    /*out*/int* adj_entity_indices_size,
			    /*inout*/int** offset,
			    /*inout*/int* offset_allocated,
			    /*out*/int* offset_size,
			    /*out*/int *err)
{
const int allocated_entity_handles = (*entity_handles_allocated == 0);
const int allocated_indices = (*adj_entity_indices_allocated == 0);
const int allocated_offset = (*offset_allocated == 0);

 // get source entities
iMesh_getEntities(instance, entity_set_handle, entity_type_requestor,
		  entity_topology_requestor, entity_handles,
		  entity_handles_allocated, entity_handles_size, err);
if (iBase_SUCCESS != *err)
return;

 // get adjacencies
iBase_EntityHandle* all_adj_handles = 0;
int size = 0, alloc = 0;
iMesh_getEntArrAdj(instance, *entity_handles, *entity_handles_size,
		   entity_type_requested, &all_adj_handles, &alloc, &size,
		   offset, offset_allocated, offset_size, err);
if (*err != iBase_SUCCESS) {
if (allocated_entity_handles) {
free(*entity_handles);
*entity_handles = 0;
*entity_handles_allocated = 0;
}
return;
}

 // allocate or check size of adj_entity_indices
*adj_entity_indices_size = size;
if (allocated_indices) {
*adj_entity_indices = reinterpret_cast<int*>(malloc(
sizeof(iBase_EntityHandle) * size));
if (!*adj_entity_indices) {
*err = iBase_MEMORY_ALLOCATION_FAILED;
instance->pGB->setError(*err);
}
else
*adj_entity_indices_allocated = size;
}
else if (*adj_entity_indices_allocated < size) {
*err = iBase_BAD_ARRAY_DIMENSION;
instance->pGB->setError(*err);
}
if (iBase_SUCCESS != *err) {
free(all_adj_handles);
if (allocated_entity_handles) {
free(*entity_handles);
*entity_handles = 0;
*entity_handles_allocated = 0;
}
if (allocated_offset) {
free(*offset);
*offset = 0;
*offset_allocated = 0;
}
return;
}

 // Now create an array of unique sorted handles from all_adj_handles.
 // We need to create a copy because we still need all_adj_handles.  We
 // will eventually need to copy the resulting unique list into
 // adj_entity_handles, so if adj_entity_handles is already allocated and
 // of sufficient size, use it rather than allocating another temporary.
iBase_EntityHandle* unique_adj = 0;
if (*adj_entity_handles_allocated >= size) {
unique_adj = *adj_entity_handles;
}
else {
unique_adj = reinterpret_cast<iBase_EntityHandle*>(malloc(
sizeof(iBase_EntityHandle) * size));
}
std::copy(all_adj_handles, all_adj_handles + size, unique_adj);
std::sort(unique_adj, unique_adj + size);
*adj_entity_handles_size = std::unique(unique_adj, unique_adj + size)
- unique_adj;

 // If we created a temporary array for unique_adj rather than using
 // already allocated space in adj_entity_handles, allocate
 // adj_entity_handles and copy the unique handle list into it
if (*adj_entity_handles != unique_adj) {
if (!*adj_entity_handles_allocated) {
*adj_entity_handles = reinterpret_cast<iBase_EntityHandle*>(malloc(
sizeof(iBase_EntityHandle) * (*adj_entity_handles_size)));
if (!*adj_entity_handles) {
*err = iBase_MEMORY_ALLOCATION_FAILED;
instance->pGB->setError(*err);
}
else
*adj_entity_handles_allocated = *adj_entity_handles_size;
}
else if (*adj_entity_handles_allocated < *adj_entity_handles_size) {
*err = iBase_BAD_ARRAY_DIMENSION;
instance->pGB->setError(*err);
}
if (iBase_SUCCESS != *err) {
free(unique_adj);
free(all_adj_handles);
if (allocated_entity_handles) {
free(*entity_handles);
*entity_handles = 0;
*entity_handles_allocated = 0;
}
if (allocated_offset) {
free(*offset);
*offset = 0;
*offset_allocated = 0;
}
if (allocated_indices) {
free(*adj_entity_indices);
*adj_entity_indices = 0;
*adj_entity_indices_allocated = 0;
}
return;
}

std::copy(unique_adj, unique_adj + *adj_entity_handles_size,
	  *adj_entity_handles);
free(unique_adj);
unique_adj = *adj_entity_handles;
}

 // convert from adjacency list to indices into unique_adj
for (int i = 0; i < *adj_entity_indices_size; ++i)
(*adj_entity_indices)[i] = std::lower_bound(
unique_adj, unique_adj + *adj_entity_handles_size, all_adj_handles[i])
- unique_adj;
free(all_adj_handles);
}
} // end extern C

