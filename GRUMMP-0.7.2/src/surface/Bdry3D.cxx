#include <ctype.h>
#include <string.h>
#include "GR_Bdry3D.h"
#include "GR_SurfMesh.h"

Bdry3D::Bdry3D(const EntContainer<Vert>& ECV) :
    BdryRep(), apBPatches(NULL), ECVerts()
{
  for (GR_index_t i = 0; i < ECV.lastEntry(); i++) {
    Vert *pVOld = ECV.getEntry(i);
    Vert *pVNew = ECVerts.getNewEntry();
    pVNew->setCoords(3, pVOld->getCoords());
    pVNew->copyAllFlagsFrom(pVOld);
    // Don't copy vertex data any deeper than this.
  }
}

bool
Bdry3D::readBdryFile(const char strBaseFileName[])
{
  char strInFileName[FILE_NAME_LEN];
  makeFileName(strInFileName, "%s.bdry", strBaseFileName,
	       "Bdry3D::readBdryFile(const string)");

  FILE* pFileInput = fopen(strInFileName, "r");
  if (pFileInput == NULL) {
    return false;
  }
  else {
    logMessage(MSG_GLOBAL, "Reading 3D boundary geometry file %s\n",
	       strInFileName);

    GR_index_t iNVerts, iNBdryPatches;
    const int iBufSize = 2048;
    char acBuffer[iBufSize];

    //@@ Read size data
    vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
    if (2 != sscanf(acBuffer, "%u%u", &iNVerts, &iNBdryPatches)) {
      vFatalError("Problem reading total number of verts and boundary patches",
		  "3D boundary geometry input");
    }
    if (iNBdryPatches >= MAX_ENTITIES)
      vFatalError("Number of boundary patches exceeds hard-coded limit",
		  "3D boundary geometry input");

    apBPatches = new BdryPatch*[iNBdryPatches];
    iMaxNPatches = iNBdryPatches;

    //@@ Read coordinate data
    ECVerts.reserve(iNVerts);
    GR_index_t i;
    for (i = 0; i < iNVerts; i++) {
      vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
      double adCoord[3];
      if (3
	  != sscanf(acBuffer, "%lf %lf %lf", &(adCoord[0]), &(adCoord[1]),
		    &(adCoord[2]))) {
	logMessage(0, "Problem with input line: %s\n", acBuffer);
	vFatalError(
	    "Insufficient data on line while reading vertex coordinates",
	    "3D boundary geometry input");
      }
      ECVerts[i].setCoords(3, adCoord);
      ECVerts[i].setType(Vert::eBdryApex);
    }

    // Switch to reading the file directly, so that line wrap will happen
    // automatically.
    int iBdryPatch = 0;
    enum BdryPatch3D::eGeom3D eGeomType = BdryPatch3D::eUnknown;
    enum {
      eContains, eSubtract, eNormal
    } eOp;
    BdryPatch3D* pBParent = NULL;
    //@@ Read boundary geometry data.
    do {
      char strLineType[33], strLowCase[33];
      if (1 != fscanf(pFileInput, "%32s", strLineType)) {
	vFatalError("Error reading boundary geometry patch type",
		    "3D boundary geometry input");
      }
      if (strLineType[0] == '#') {
	// This line is a comment; skip it.
	char* pcRes = fgets(acBuffer, iBufSize, pFileInput);
	if (!pcRes)
	  vFatalError("Comment too long?", "3D boundary geometry input");
	continue;
      }
      for (i = 0; i < 33 && strLineType[i] != '\000'; i++)
	strLowCase[i] = tolower(strLineType[i]);
      //@@@ Identify the operation type, if any
      eOp = eNormal;
      if (strncmp(strLowCase, "subtract", 8) == 0) {
	if (pBParent == NULL) {
	  vFatalError("Error: subtracting a hole with nothing to subtract from",
		      "3D boundary geometry input");
	}
	eOp = eSubtract;
	// Now read the next word to find out what kind of bdry entity
	// this is.
	if (1 != fscanf(pFileInput, "%32s", strLineType)) {
	  vFatalError("Error reading boundary geometry patch type",
		      "3D boundary geometry input");
	}
	for (i = 0; i < 33 && strLineType[i] != '\000'; i++)
	  strLowCase[i] = tolower(strLineType[i]);
      }
      else if (strncmp(strLowCase, "contains", 8) == 0) {
	if (pBParent == NULL) {
	  vFatalError("Error: patch contained within nothing",
		      "3D boundary geometry input");
	}
	eOp = eContains;
	// Now read the next word to find out what kind of bdry entity
	// this is.
	if (1 != fscanf(pFileInput, "%32s", strLineType)) {
	  vFatalError("Error reading boundary geometry patch type",
		      "3D boundary geometry input");
	}
	for (i = 0; i < 33 && strLineType[i] != '\000'; i++)
	  strLowCase[i] = tolower(strLineType[i]);
      }

      //@@@ Identify the type of geometry entity
      if (strncmp(strLowCase, "polygon", 7) == 0) {
	eGeomType = BdryPatch3D::ePolygon;
      }
      else {
	logMessage(0, "Unknown geometric entity name: %s\n", strLineType);
	vFatalError("Error reading boundary geometry data",
		    "3D boundary geometry input");
      }

      //@@@ Need to grab BC/region data unless this is a hole in the surface.
      int iBCL, iBCR, iRegL, iRegR;
      /* Initialize for picky compilers. */
      iBCL = iBCR = iInvalidBC;
      iRegL = iRegR = iInvalidRegion;

      if (eOp != eSubtract) {
	char cTypeL, cTypeR;
	int iValueL, iValueR;
	if (4
	    != fscanf(pFileInput, " %1c%d %1c%d", &cTypeL, &iValueL, &cTypeR,
		      &iValueR)) {
	  logMessage(0,
		     "Bad boundary condition or region specification in:\n%s\n",
		     acBuffer);
	  vFatalError("Error reading boundary geometry data",
		      "3D boundary geometry input");
	}

	if ((cTypeL == 'b') || (cTypeL == 'B')) {
	  if (iValueL <= 0)
	    vFatalError("Illegal boundary condition",
			"3D boundary geometry input");
	  iBCL = iValueL;
	  iRegL = iOutsideRegion;
	}
	else {
	  if (cTypeL != 'r' && cTypeL != 'R')
	    vFatalError("Must specify either a boundary condition or a region",
			"3D boundary geometry input");
	  if (iValueL <= iInvalidRegion || iValueL > iMaxRegion)
	    vFatalError("Illegal region specification",
			"3D boundary geometry input");
	  iBCL = iInvalidBC;
	  iRegL = iValueL;
	}
	if ((cTypeR == 'b') || (cTypeR == 'B')) {
	  if (iValueR <= 0)
	    vFatalError("Illegal boundary condition",
			"3D boundary geometry input");
	  iBCR = iValueR;
	  iRegR = iOutsideRegion;
	}
	else {
	  if (cTypeR != 'r' && cTypeR != 'R')
	    vFatalError("Must specify either a boundary condition or a region",
			"3D boundary geometry input");
	  if (iValueR <= iInvalidRegion || iValueR > iMaxRegion)
	    vFatalError("Illegal region specification",
			"3D boundary geometry input");
	  iBCR = iInvalidBC;
	  iRegR = iValueR;
	}
      }

      switch (eGeomType)
	{
	case BdryPatch3D::ePolygon:
	  {
	    int iSeg, iNumSegs;
	    if (1 != fscanf(pFileInput, "%d", &iNumSegs)) {
	      vFatalError("Error reading boundary geometry data",
			  "3D boundary geometry input");
	    }
	    if (iNumSegs < 0 || iNumSegs > 2048) {
	      vFatalError("Invalid number of segments for a bdry polygon",
			  "3D boundary geometry input");
	    }
	    Vert* *apV = new Vert*[iNumSegs];
	    for (iSeg = 0; iSeg < iNumSegs; iSeg++) {
	      GR_index_t iV;
	      if (1 != fscanf(pFileInput, "%u", &iV))
		vFatalError("Error reading boundary geometry data",
			    "3D boundary geometry input");
	      if (iV > ECVerts.lastEntry())
		vFatalError("Vertex index out of range",
			    "3D boundary geometry input");
	      apV[iSeg] = ECVerts.getEntry(iV);
	    }
	    BdryPolygon *pBPoly = new BdryPolygon(iBCL, iBCR, iRegL, iRegR,
						  iNumSegs, apV);
	    //@@@ Now hook it into the geometry heirarchy
	    switch (eOp)
	      {
	      case eNormal:
		// Just an ordinary surface
		vAddPatch(pBPoly);
		iBdryPatch++;
		pBParent = pBPoly;
		break;
	      case eSubtract:
		// A hole in a surface
		pBParent->vAddChild(pBPoly);
		break;
	      case eContains:
		// A subset of a surface with a different BC or different
		// region than its parent surface.
		pBParent->vAddChild(pBPoly);
		vAddPatch(pBPoly);
		iBdryPatch++;
		break;
	      } // Switch on operation
	    delete[] apV;
	    break;
	  } // End case for polygons
	case BdryPatch3D::eUnknown:
	  {
	    assert(0);
	    break;
	  }
	} // End switch to read all geometry types.
    }
    while (iBdryPatch < iMaxNPatches); // End loop to read all boundary patches
    assert(iBdryPatch == iMaxNPatches);
    fclose(pFileInput);
    return true;
  } // Done reading .bdry file
}

bool
Bdry3D::readSMeshFile(const char strBaseFileName[])
{
  char strInFileName[FILE_NAME_LEN];
  makeFileName(strInFileName, "%s.smesh", strBaseFileName,
	       "Bdry3D::Bdry3D(const string)");

  FILE* pFileInput = fopen(strInFileName, "r");
  if (pFileInput == NULL) {
    return false;
  }
  else {
    logMessage(0, "Reading 3D boundary geometry file %s\n", strInFileName);
    const int iBufSize = 2048;
    char acBuffer[iBufSize];

    //@@ Determine file format

    // The first line of the input indicates file format.
    // Strings starting with "poly" are polygon files
    // Strings starting with "face" are face-based triangle files
    // Other strings are errors
    // The rest of the line includes size data.
    vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
    enum {
      ePolygons, eFaceBased, eUnknown
    } eFileType;

    GR_index_t iNVerts, iNEdges, iNCells;
    if (strncmp(acBuffer, "poly", 4) == 0) {
      eFileType = ePolygons;
      if (2 != sscanf(acBuffer, "%*s %d %d", &iNVerts, &iNCells))
	vFatalError("Couldn't read surface mesh size info",
		    "3D boundary geometry input");
    }
    else if (strncmp(acBuffer, "face", 4) == 0
	|| strncmp(acBuffer, "edge", 4) == 0) {
      eFileType = eFaceBased;
      GR_index_t iNSurfs;
      if (4
	  != sscanf(acBuffer, "%*s %u %u %u %u", &iNCells, &iNEdges, &iNSurfs,
		    &iNVerts))
	vFatalError("Couldn't read surface mesh size info",
		    "3D boundary geometry input");
      // Number of quads had better be zero
      assert(3 * iNCells == 2 * iNEdges);
      // The Euler formula for the surface mesh has to hold unless there
      // are handles, but the only easy way to check for handles is to
      // see if the following is an equality.  Not worth the trouble.
      // assert(iNCells + iNVerts == iNEdges + 2*iNSurfs);
    }
    else {
      eFileType = eUnknown;
    }

    //@@ Read coordinates
    ECVerts.reserve(iNVerts);
    GR_index_t i;
    for (i = 0; i < iNVerts; i++) {
      vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
      Vert *pV = pVVert(i);
      pV->setSpaceDimen(3);
      pV->setCoords(acBuffer);
      pV->setType(Vert::eBdryApex);
    }

    //@@ Read connectivity
    switch (eFileType)
      {
      case eUnknown:
	vFatalError("Unrecognized file type reading surface mesh",
		    "3D boundary geometry input");
	break;
      case eFaceBased:
	{
	  // Face-cell data.  Just read it as a SurfMesh, and use the
	  // built-in functions of faces and cells to get back the data
	  // that you need.

	  SurfMesh SM(ECVerts);

	  if (iNCells > MAX_ENTITIES) {
	    vFatalError("too many triangles in surface mesh",
			"3D boundary geometry input");
	  }
	  int *aiBC = new int[iNCells];
	  for (i = 0; i < iNCells; i++) {
	    aiBC[i] = iDefaultBC;
	  }
	  SM.setupFaces(iNEdges);       // Ensure the existence of enough
	  SM.setupTriCells(iNCells);    // blank edges and cells.
	  int iC0, iC1, iV0, iV1;
	  for (i = 0; i < iNEdges; i++) {
	    vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
	    if (4 != sscanf(acBuffer, "%d %d %d %d", &iC0, &iC1, &iV0, &iV1))
	      vFatalError("Error reading mesh connectivity",
			  "3D boundary input");
	    // The indices we just read don't get validated until -inside-
	    // the calls to pCCell and pVVert.  Klocwork doesn't like
	    // this, and it could be a problem with optimized code, but
	    // not in debug mode.  In the end, a sufficiently data file
	    // can crash the code.
	    SM.getFace(i)->assign(SM.getCell(iC0), SM.getCell(iC1),
				  SM.getVert(iV0), SM.getVert(iV1));
	    // Create cell-face pointers on the fly
	    SM.getCell(iC0)->addFace(SM.getFace(i));
	    SM.getCell(iC1)->addFace(SM.getFace(i));
	  }

	  while (!feof(pFileInput)) {
	    int iCell, iBCTmp, iRes;
	    iRes = fscanf(pFileInput, "%d %d", &iCell, &iBCTmp);
	    if (iRes == 2) {
	      assert(size_t(iCell) < iNCells && iCell >= 0);
	      assert(iBCTmp > 0);
	      aiBC[iCell] = iBCTmp;
	    }
	    else if (!feof(pFileInput))
	      vFatalError("Error reading boundary conditions",
			  "3D boundary input");
	  }
	  // Now go through and set up one BdryPolygon for each cell.
	  iMaxNPatches = iNCells;
	  apBPatches = new BdryPatch*[iMaxNPatches];
	  for (i = 0; i < iNCells; i++) {
	    Vert *apV[3];
	    apV[0] = pVVert(SM.iVertIndex(SM.getCell(i)->getVert(0)));
	    apV[1] = pVVert(SM.iVertIndex(SM.getCell(i)->getVert(1)));
	    apV[2] = pVVert(SM.iVertIndex(SM.getCell(i)->getVert(2)));

	    BdryPolygon *pBPoly = new BdryPolygon(aiBC[i], iInvalidBC,
	    iOutsideRegion,
						  iDefaultRegion, 3, apV);
	    vAddPatch(pBPoly);
	  }
	  delete[] aiBC;
	  break;
	} // End block for case: reading edge-cell data
      case ePolygons:
	{
	  //@@@ Polygonal data
	  int iNPoly = iNCells;
	  int iPoly;
	  iMaxNPatches = iNPoly;
	  apBPatches = new BdryPatch*[iNPoly];

	  logMessage(MSG_MANAGER, "Reading data for %d polygons...\n", iNPoly);
	  int iPrintInt = iNPoly / 10 + 1;
	  if (iPrintInt < 25)
	    iPrintInt = 25;
	  for (iPoly = 0; iPoly < iNPoly; iPoly++) {
	    //@@@@ Read polygonal surface, with possible holes
	    // Find the file position so that we can rewind over the line if
	    // this is a normal polygon
	    long lFPos = ftell(pFileInput);

	    vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
	    const char acMP[] = "multipoly";
	    int iNHoles = 0;
	    // Read perimeter polygon (or only polygon)
	    if (strncmp(acBuffer, acMP, 9) == 0) {
	      //@@@@@ Polygon with holes
	      // Determine the number of holes
	      int iNNested;
	      if (1 != sscanf(acBuffer, "%*s %d", &iNNested))
		vFatalError("Number of co-planar polygons not specified",
			    "Nested polygon input");
	      // The number read is actually the number of total polygons;
	      // the first is the perimeter, the rest are holes in the
	      // interior.
	      iNHoles = iNNested - 1;
	    }
	    else {
	      // Rewind the stream to put you back at the beginning of
	      // this data
	      fseek(pFileInput, lFPos, 0);
	      // Now skip comment lines and blank lines.
	      vSkipCommentLines(pFileInput);
	    }

	    // Read the size of the polygon
	    int iNumSegs;
	    if (1 != fscanf(pFileInput, "%d", &iNumSegs)) {
	      vFatalError("Error reading number of polygon sides",
			  "3D boundary geometry input");
	    }
	    if (iNumSegs <= 2)
	      vFatalError("Polygon must have at least three sides",
			  "3D boundary geometry input");
	    // Now read the polygon, and make a list of pointers to its
	    // vertices.
	    Vert* *apV = new Vert*[iNumSegs];
	    int iSeg;
	    for (iSeg = 0; iSeg < iNumSegs; iSeg++) {
	      GR_index_t iV;
	      if (1 != fscanf(pFileInput, "%u", &iV))
		vFatalError("Error reading boundary geometry data",
			    "3D boundary geometry input");
	      if (iV > ECVerts.lastEntry())
		vFatalError("Vertex index out of range",
			    "3D boundary geometry input");
	      apV[iSeg] = ECVerts.getEntry(iV);
	    }
	    // Read the optional BC
	    char *pcRes = fgets(acBuffer, iBufSize, pFileInput);
	    assert(pcRes == acBuffer);
	    int iBC = iDefaultBC;
	    {
	      int iNRead = sscanf(acBuffer, "%d", &iBC);
	      if (iNRead < 1) { // No optional BC present
		iBC = iDefaultBC;
	      }
	    }

	    // Set up the BdryPoly for the perimeter.  This is
	    // deliberately -not- deleted, because this is durable
	    // boundary data for the meshing process.

	    BdryPolygon *pBPoly = new BdryPolygon(iBC, iInvalidBC,
	    iOutsideRegion,
						  iDefaultRegion, iNumSegs,
						  apV);
	    delete[] apV;
	    vAddPatch(pBPoly);
	    BdryPolygon *pBParent = pBPoly;

	    BdryPolygon *pBPolyHole = NULL;
	    // Read each of the holes in turn.
	    for (int iHole = 0; iHole < iNHoles; iHole++) {
	      // Read the hole data, as above.
	      if (1 != fscanf(pFileInput, "%d", &iNumSegs)) {
		vFatalError("Error reading number of polygon sides",
			    "3D boundary geometry input");
	      }
	      if (iNumSegs <= 2)
		vFatalError("Polygon must have at least three sides",
			    "3D boundary geometry input");
	      // FIXME: Need to set this up so that polygon orientation
	      // gets reversed when needed to get its inside to really be
	      // on the inside.

	      // Now read the polygon, and make a list of pointers to its
	      // vertices.
	      apV = new Vert*[iNumSegs];
	      for (iSeg = 0; iSeg < iNumSegs; iSeg++) {
		GR_index_t iV;
		if (1 != fscanf(pFileInput, "%u", &iV))
		  vFatalError("Error reading boundary geometry data",
			      "3D boundary geometry input");
		if (iV > ECVerts.lastEntry())
		  vFatalError("Vertex index out of range",
			      "3D boundary geometry input");
		apV[iSeg] = ECVerts.getEntry(iV);
	      }

	      // Read the optional BC
	      pcRes = fgets(acBuffer, iBufSize, pFileInput);
	      assert(pcRes == acBuffer);
	      {
		int iNRead = sscanf(acBuffer, "%d", &iBC);
		if (iNRead < 1) { // No optional BC present
		  iBC = iInvalidBC;
		}
	      }

	      // Set up the BdryPoly for the hole.  These are deliberately
	      // -not- deleted, because this is durable boundary data for
	      // the meshing process.
	      pBPolyHole = new BdryPolygon(iBC, iInvalidBC,
	      iOutsideRegion,
					   iDefaultRegion, iNumSegs, apV);
	      delete[] apV;
	      // Whether this is a hole (Subtract) or a part of the
	      // surface with a different BC or region (Contain), it's
	      // still a child:
	      pBParent->vAddChild(pBPolyHole);

	      // If it's a part of the surface with a different BC or
	      // region (Contain), it needs to be part of the boundary
	      // description in its own right:
	      if (iBC != iInvalidBC)
		vAddPatch(pBPolyHole);
	    } // Done reading and storing holes
	    if ((iPoly + 1) % iPrintInt == 0)
	      logMessage(MSG_SERVICE, "Read polygon %d\n", iPoly);
	  } // Loop over all polygons
	  if ((iNPoly + 1) % iPrintInt != 0)
	    logMessage(MSG_SERVICE, "Read polygon %d\n", iNPoly);
	  break;
	} // end of case for reading polygonal data
      } // end of switch statement for reading connectivity data
    fclose(pFileInput);
    return true;
  }
}

bool
Bdry3D::readHINFile(const char strBaseFileName[])
{
  char strInFileName[FILE_NAME_LEN];
  if (!strstr(strBaseFileName, ".hin")) {
    makeFileName(strInFileName, "%s.hin", strBaseFileName,
		 "Bdry3D::Bdry3D(const string)");
  }
  else {
    strcpy(strInFileName, strBaseFileName);
  }
  FILE* pFileInput = fopen(strInFileName, "r");
  if (pFileInput == NULL) {
    logMessage(MSG_GLOBAL, "Unsuccessfully read %s HIN File\n",
	       strBaseFileName);

    return false;
  }
  else {
    logMessage(MSG_GLOBAL, "Reading 3D boundary geometry file %s\n",
	       strInFileName);
    // Read number of vertices.
    int iNVerts;
    int result = fscanf(pFileInput, "Vertices[%d]{", &iNVerts);
    if (result != 1) {
      vFatalError("Couldn't read number of vertices.",
		  "Reading Intel polyhedron file");
    }

    //@@ Read coordinates
    ECVerts.reserve(iNVerts);
    int i;
    for (i = 0; i < iNVerts; i++) {
      const int iBufSize = 2048;
      char acBuffer[iBufSize];
      vGetLineOrAbort(acBuffer, iBufSize, pFileInput);
      Vert *pV = pVVert(i);
      pV->setSpaceDimen(3);
      pV->setCoords(acBuffer);
      pV->setType(Vert::eBdryApex);
    }

    //@@ Polygon data
    int iNPoly;
    result = fscanf(pFileInput, "}\nPolygons[%d]{\n", &iNPoly);
    if (result != 1) {
      vFatalError("Couldn't read number of polygons.",
		  "Reading Intel polyhedron file");
    }

    iMaxNPatches = iNPoly;
    apBPatches = new BdryPatch*[iNPoly];
    Vert ***a2pV = new Vert**[iNPoly];
    int *aiNPolyVerts = new int[iNPoly];

    logMessage(MSG_MANAGER, "Reading data for %d polygons...\n", iNPoly);
    int iPrintInt = std::min(iNPoly / 10 + 1, 25);

    for (int iPoly = 0; iPoly < iNPoly; iPoly++) {
      //@@@ Read polygons
      int iNPolyVerts;
      result = fscanf(pFileInput, "Polygon[%d]{", &iNPolyVerts);
      if (result != 1) {
	vFatalError("Couldn't read number of vertices in a polygon.",
		    "Reading Intel polyhedron file");
      }

      if (iNPolyVerts <= 2)
	vFatalError("Polygon must have at least three sides",
		    "3D boundary geometry input");

      // Now read the polygon, and make a list of pointers to its
      // vertices.
      a2pV[iPoly] = new Vert*[iNPolyVerts];
      aiNPolyVerts[iPoly] = iNPolyVerts;
      for (int iVert = 0; iVert < iNPolyVerts; iVert++) {
	int iV;
	if (1 != fscanf(pFileInput, "%d", &iV))
	  vFatalError("Error reading boundary geometry data",
		      "3D boundary geometry input");
	if (unsigned(iV) > ECVerts.lastEntry() || iV < 0)
	  vFatalError("Vertex index out of range",
		      "3D boundary geometry input");
	a2pV[iPoly][iVert] = ECVerts.getEntry(iV);
      }
      int iNRead = fscanf(pFileInput, "\n}\n");
      assert(iNRead == 0);
      if ((iPoly + 1) % iPrintInt != 0)
	logMessage(MSG_SERVICE, "Read polygon %d\n", iPoly);
    }

    // Now read polyhedra (regions) and tag polygons.
    int *aiRightReg = new int[iNPoly];
    int *aiLeftReg = new int[iNPoly];
    for (int iPoly = 0; iPoly < iNPoly; iPoly++) {
      aiRightReg[iPoly] = aiLeftReg[iPoly] = iOutsideRegion;
    }

    // check for facemarkers, don't do anything with them for now
    int iNFaceMarkers;
    fpos_t iPos;
    fgetpos(pFileInput, &iPos);
    result = fscanf(pFileInput, "}\nFacemarkers[%d]{\n", &iNFaceMarkers);
    if (result != 1) {
      fsetpos(pFileInput, &iPos);
    }
    else {
      for (int iFaceMarker = 0; iFaceMarker < iNFaceMarkers; iFaceMarker++) {
	int iNFaceMarks, iFirstFaceMark;
	result = fscanf(pFileInput, "Facemark[%d]{\n%d", &iNFaceMarks,
			&iFirstFaceMark);
	if (result != 2) {
	  vFatalError("Couldn't read facemarkers.",
		      "Reading Intel polyhedron file");
	}
	for (int iFaceMark = 0; iFaceMark < iNFaceMarks; iFaceMark++) {
	  int iFM;
	  if (1 != fscanf(pFileInput, "%d", &iFM))
	    vFatalError("Error reading facemarkers",
			"Reading Intel polyhedron file");
	}
	int iNRead = fscanf(pFileInput, "\n}\n");
	assert(iNRead == 0);
      }
    }
    int iNHedra;
    result = fscanf(pFileInput, "}\nPolyhedra[%d]{\n", &iNHedra);
    if (result != 1) {
      vFatalError("Couldn't read number of polyhedra.",
		  "Reading Intel polyhedron file");
    }
    for (int iHedra = 0; iHedra < iNHedra; iHedra++) {
      //@@@ Read polygons
      int iNHedraFaces;
      result = fscanf(pFileInput, "Polyhedron[%d]{", &iNHedraFaces);
      if (result != 1) {
	vFatalError("Couldn't read number of faces in a polyhedron.",
		    "Reading Intel polyhedron file");
      }

      if (iNHedraFaces <= 3)
	vFatalError("Polyhedron must have at least four faces",
		    "3D boundary geometry input");

      // Now read the polygons, and tag them with regions.
      for (int iFace = 0; iFace < iNHedraFaces; iFace++) {
	int iP;
	if (1 != fscanf(pFileInput, "%d", &iP))
	  vFatalError("Error reading boundary geometry data",
		      "3D boundary geometry input");
	if (iP >= 0)
	  aiLeftReg[iP] = iHedra;
	else
	  aiRightReg[-iP - 1] = iHedra;
      }
      int iNRead = fscanf(pFileInput, "\n}\n");
      assert(iNRead == 0);
    }

    // check again for facemarkers, don't do anything with them for now
    fgetpos(pFileInput, &iPos);
    result = fscanf(pFileInput, "}\nFacemarkers[%d]{\n", &iNFaceMarkers);
    if (result != 1) {
      fsetpos(pFileInput, &iPos);
    }
    else {
      for (int iFaceMarker = 0; iFaceMarker < iNFaceMarkers; iFaceMarker++) {
	int iNFaceMarks, iFirstFaceMark;
	result = fscanf(pFileInput, "Facemark[%d]{\n%d", &iNFaceMarks,
			&iFirstFaceMark);
	if (result != 2) {
	  vFatalError("Couldn't read facemarkers.",
		      "Reading Intel polyhedron file");
	}
	for (int iFaceMark = 0; iFaceMark < iNFaceMarks; iFaceMark++) {
	  int iFM;
	  if (1 != fscanf(pFileInput, "%d", &iFM))
	    vFatalError("Error reading facemarkers",
			"Reading Intel polyhedron file");
	}
	int iNRead = fscanf(pFileInput, "\n}\n");
	assert(iNRead == 0);
      }
    }
    int iNMaterials;
    result = fscanf(pFileInput, "}\nMaterial[%d]{\n", &iNMaterials);
    if (result != 1) {
      vFatalError("Couldn't read materials", "Reading Intel polyhedron file");
    }
    assert(iNMaterials == iNHedra);
    // intentionally off by 1 to avoid 0 region number
    int *aiMaterials = new int[iNMaterials];
    for (int iMaterial = 0; iMaterial < iNMaterials; iMaterial++) {
      int iM;
      if (1 != fscanf(pFileInput, "%d", &iM))
	vFatalError("Error reading facemarkers",
		    "Reading Intel polyhedron file");
      aiMaterials[iMaterial] = iM + 1;
    }

    for (int iPoly = 0; iPoly < iNPoly; iPoly++) {
      if (aiLeftReg[iPoly] != iOutsideRegion)
	aiLeftReg[iPoly] = aiMaterials[aiLeftReg[iPoly]];
      if (aiRightReg[iPoly] != iOutsideRegion)
	aiRightReg[iPoly] = aiMaterials[aiRightReg[iPoly]];
    }

    for (int iPoly = 0; iPoly < iNPoly; iPoly++) {
      int iRightBC = (aiRightReg[iPoly] == iOutsideRegion ?
      iDefaultBC :
							    iInvalidBC);
      int iLeftBC = (aiLeftReg[iPoly] == iOutsideRegion ?
      iDefaultBC :
							  iInvalidBC);
      BdryPolygon *pBPoly = new BdryPolygon(iRightBC, iLeftBC,
					    aiRightReg[iPoly], aiLeftReg[iPoly],
					    aiNPolyVerts[iPoly], a2pV[iPoly]);
      logMessage(MSG_DETAIL,
		 "Poly: %d (%d verts).  Regions: %d, %d.  BCs: %d %d\n", iPoly,
		 aiNPolyVerts[iPoly], aiLeftReg[iPoly], aiRightReg[iPoly],
		 iLeftBC, iRightBC);
      logMessage(MSG_DETAIL, "   Verts: ");
      for (int i__ = 0; i__ < aiNPolyVerts[iPoly]; i__++) {
	logMessage(MSG_DETAIL, "%5d ", iVertIndex(a2pV[iPoly][i__]));
      }
      logMessage(MSG_DETAIL, "\n");
      vAddPatch(pBPoly);
      delete[] a2pV[iPoly];
    }

    delete[] aiMaterials;
    delete[] a2pV;
    delete[] aiNPolyVerts;
    delete[] aiRightReg;
    delete[] aiLeftReg;

    logMessage(MSG_MANAGER, "Successfully read %s HIN File\n", strBaseFileName);
    return true;
  }
}

Bdry3D::Bdry3D(const char strBaseFileName[], const int iOptLevel) :
    BdryRep(), apBPatches(NULL), ECVerts()
{
  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //  Try to find and read a .bdry file to create the Bdry3D object.  //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////
  bool success = readBdryFile(strBaseFileName);

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  // Following code tries (in the absence of a .bdry file) to read    //
  // an old-style .smesh file and convert it to a Bdry3D object. The  //
  // task is complicated by the fact that there are -two- different   //
  // formats for .smesh files: edge-based and polygon.  They are both //
  // relatively easy to convert, at least in principle, but the task  //
  // is somewhat tedious.                                             //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  if (!success) {
    success = readSMeshFile(strBaseFileName);
  }

  if (!success) {
    success = readHINFile(strBaseFileName);
  }

  if (!success) {
    vFatalError("Couldn't find a file with that base name for reading.",
		"3D bdry data input");
  }

  // At some point, we need to be checking the Bdry3D object for
  // orientation.  The easiest way to do this, in principle, is to
  // integrate x \hat{i} (or equivalent) over each BdryPatch that is a
  // true boundary to determine whether the volume is positive or
  // negative.  In practice, this is easy for polyhedra and much more
  // difficult for arbitrarily curved surfaces.

  vOptimizePatches(iOptLevel);
}

void
Bdry3D::vAddPatch(BdryPatch3D* const pBCh)
{
  if (iNPatches == iMaxNPatches) {
    iMaxNPatches = max(10, 2 * iMaxNPatches); // Effectively gets us off zero.
    BdryPatch **apNewB = new BdryPatch*[iMaxNPatches];
    for (int iC = 0; iC < iNPatches; iC++) {
      apNewB[iC] = apBPatches[iC];
    }
    delete[] apBPatches;
    apBPatches = apNewB;
  }
  apBPatches[iNPatches] = pBCh;
  iNPatches++;
}

Vert*
Bdry3D::pVAddVertex(const Vert* const pV0, const Vert* const pV1,
		    const double adCent[])
{
  // Create a vertex object with the appropriate coordinates
  Vert *pV = pVNewVert();
  pV->setCoords(3, adCent);
  pV->setType(Vert::eBdryCurve);
  // For every piece of the defining boundary, add this vertex where it
  // belongs, if it does.
  for (int i = iNumBdryPatches() - 1; i >= 0; i--)
    pBdry(i)->vSplitSegmentIfFound(pV0, pV1, pV);
  return pV;
}

void
BdryPatch3D::vAddChild(BdryPatch3D* const pB3D)
{
  assert(pB3D != NULL);
  if (pB3DChild == NULL) {
    // No children yet
    pB3DChild = pB3D;
  }
  else {
    // Add another child when one already exists
    BdryPatch3D* pB3DSib = pB3DChild;
    while (pB3DSib->pB3DNextSib != NULL)
      pB3DSib = pB3DSib->pB3DNextSib;
    pB3DSib->pB3DNextSib = pB3D;
  }
}

int
BdryPatch3D::iNumChildren() const
{
  int iRetVal = 0;
  if (pB3DChild != NULL) {
    BdryPatch3D* pB3DSib = pB3DChild;
    iRetVal = 1 + pB3DSib->iNumChildren();
    while (pB3DSib->pB3DNextSib != NULL) {
      pB3DSib = pB3DSib->pB3DNextSib;
      iRetVal += 1 + pB3DSib->iNumChildren();
    }
  }
  return (iRetVal);
}

int
BdryPatch3D::iNumSegments() const
{
  // For a given BdryPatch, the number of segments is the same as the
  // number of its vertices.  Also, have to add in children.
  int iRetVal = iNVerts;
  if (pB3DChild != NULL) {
    BdryPatch3D* pB3DSib = pB3DChild;
    iRetVal += pB3DSib->iNumSegments();
    while (pB3DSib->pB3DNextSib != NULL) {
      pB3DSib = pB3DSib->pB3DNextSib;
      iRetVal += pB3DSib->iNumSegments();
    }
  }
  return (iRetVal);
}

void
BdryPatch3D::vGetSegmentVerts(const int iSegIn, Vert *& pV0, Vert *& pV1) const
{
  // Segments are numbered sequentially beginning with the segment
  // connecting vert 0 and vert 1 in the parent; running through the
  // segment connecting iNVerts-1 to 0 in the parent; then proceeding
  // depth-first through all children.
  int iSeg = iSegIn;
  assert(iSeg < iNumSegments() && iSeg >= 0);
  if (iSeg < iNVerts) {
    // It's in this Patch; set up the data for return
    pV0 = apVVerts[iSeg];
    pV1 = apVVerts[(iSeg + 1) % iNVerts]; // Careful about wraparound!
    return;
  }
  // Decrement iSeg so that the alignment is right for children.
  iSeg -= iNVerts;
  // If we've gotten here, then there -must- be a child for this segment
  // to belong to.
  BdryPatch3D* pB3DSib = pB3DChild;
  if (pB3DSib->iNVerts > iSeg) {
    // The child owns this segment.
    pB3DSib->vGetSegmentVerts(iSeg, pV0, pV1);
    Vert *pVTmp = pV0;
    pV0 = pV1;
    pV1 = pVTmp;
    return;
  }
  iSeg -= pB3DSib->iNVerts;

  while (pB3DSib->pB3DNextSib != NULL) {
    pB3DSib = pB3DSib->pB3DNextSib;
    if (pB3DSib->iNVerts > iSeg) {
      // The child owns this segment.
      pB3DSib->vGetSegmentVerts(iSeg, pV0, pV1);
      Vert *pVTmp = pV0;
      pV0 = pV1;
      pV1 = pVTmp;
      return;
    }
    iSeg -= pB3DSib->iNVerts;
  }
}

void
BdryPatch3D::vSplitSegmentIfFound(const Vert* const pV0, const Vert* const pV1,
				  Vert *pV)
{
  // Does this BdryPatch have the first of the two required end verts?
  if (!qHasVert(pV0))
    return;

  // Now see if we can find the other and insert.

  for (int iV = iNVerts - 1; iV >= 0; iV--) {
    if (apVVerts[iV] == pV0) {
      // It's in this Patch; check neighbors
      if (apVVerts[(iV + 1) % iNVerts] == pV1) {
	// Found it!  Insert.
	vResize(iNVerts + 1);
	for (int i = iNVerts; i > iV + 1; i--)
	  apVVerts[i] = apVVerts[i - 1];
	apVVerts[iV + 1] = pV;
	// Update the size of this BdryPatch
	iNVerts++;
      }
      else if (apVVerts[(iV == 0 ? iNVerts - 1 : iV - 1)] == pV1) {
	// Found it!  Insert.
	vResize(iNVerts + 1);
	for (int i = iNVerts; i > iV; i--)
	  apVVerts[i] = apVVerts[i - 1];
	apVVerts[iV] = pV;
	// Update the size of this BdryPatch
	iNVerts++;
      }
      return;
    }
  }
  // If we've gotten here, then there -must- be a child that this vert
  // belongs to.
  BdryPatch3D* pB3DSib = pB3DChild;
  if (pB3DSib->qHasVert(pV0)) {
    // The child owns this segment.
    pB3DSib->vSplitSegmentIfFound(pV0, pV1, pV);
    return;
  }
  while (pB3DSib->pB3DNextSib != NULL) {
    pB3DSib = pB3DSib->pB3DNextSib;
    if (pB3DSib->qHasVert(pV0)) {
      // The child owns this segment.
      pB3DSib->vSplitSegmentIfFound(pV0, pV1, pV);
      return;
    }
  }
}

bool
BdryPatch3D::qHasVert(const Vert* const pV) const
{
  // Check to see whether the vert appears directly here, then search
  // depth-first through children.
  for (int i = iNVerts - 1; i >= 0; i--)
    if (apVVerts[i] == pV)
      return (true);
  if (pB3DChild != NULL) {
    BdryPatch3D* pB3DSib = pB3DChild;
    if (pB3DSib->qHasVert(pV))
      return (true);
    while (pB3DSib->pB3DNextSib != NULL) {
      pB3DSib = pB3DSib->pB3DNextSib;
      if (pB3DSib->qHasVert(pV))
	return (true);
    }
  }
  // If we got here, nobody knows about this vertex.
  return (false);
}

void
BdryPatch3D::vGetVertNeighbors(const Vert * const pV, Vert* &pVN0,
			       Vert* &pVN1) const
{
  assert(qHasVert(pV));
  for (int i = iNVerts - 1; i >= 0; i--) {
    if (apVVerts[i] == pV) {
      // Found it!
      // Set the two neighboring verts
      if (i == 0)
	pVN0 = apVVerts[iNVerts - 1];
      else
	pVN0 = apVVerts[i - 1];
      pVN1 = apVVerts[(i + 1) % iNVerts];
      return;
    }
  }
  if (pB3DChild != NULL) {
    BdryPatch3D* pB3DSib = pB3DChild;
    if (pB3DSib->qHasVert(pV)) {
      // It's in here.
      pB3DSib->vGetVertNeighbors(pV, pVN0, pVN1);
      return;
    }
    while (pB3DSib->pB3DNextSib != NULL) {
      pB3DSib = pB3DSib->pB3DNextSib;
      if (pB3DSib->qHasVert(pV)) {
	// It's in here.
	pB3DSib->vGetVertNeighbors(pV, pVN0, pVN1);
	return;
      }
    }
  }
  // Should never get here.
  assert2(0, "Something very fishy in BdryPatch3D::vGetVertNeighbors");
}

bool
BdryPatch3D::qDisjointFrom(const Vert* const pV) const
{
  return (iFuzzyComp(0., dDistFromBdry(pV->getCoords())) != 0);
}

bool
BdryPatch3D::qDisjointFrom(const BdryPatch* const pBP) const
{
  bool qRetVal = true;
  for (int iV = 0; iV < iNVerts; iV++) {
    qRetVal = qRetVal && pBP->qDisjointFrom(pVVert(iV));
  }
  // Now have children check themselves recursively to be sure -they-
  // are disjoint from pBP.
  if (pB3DChild != NULL) {
    BdryPatch3D* pB3DSib = pB3DChild;
    qRetVal = qRetVal && pB3DSib->qDisjointFrom(pBP);
    while (pB3DSib->pB3DNextSib != NULL) {
      pB3DSib = pB3DSib->pB3DNextSib;
      qRetVal = qRetVal && pB3DSib->qDisjointFrom(pBP);
    }
  }
  return qRetVal;
}

Vert*
BdryPatch3D::pVVert(const int i) const
{
  assert(i >= 0 && i < iNumSegments());
  if (i < iNVerts)
    return apVVerts[i];
  else {
    assert(pB3DChild != NULL);
    BdryPatch3D* pB3DSib = pB3DChild;
    int iIndex = i - iNVerts;
    while (iIndex >= pB3DSib->iNumSegments()) {
      iIndex -= pB3DSib->iNVerts;
      pB3DSib = pB3DSib->pB3DNextSib;
    }
    assert(iIndex >= 0 && iIndex < pB3DSib->iNVerts);
    return pB3DSib->pVVert(iIndex);
  }
}

void
Bdry3D::vOptimizePatches(const int iLevel)
{
  logMessage(MSG_MANAGER, "Optimizing 3D bdry representation...\n");
  logMessage(MSG_SERVICE, "  Initially requiring %u triangles.\n",
	     iNumTriangles());
  // Merge patches
  if (iLevel >= 1) {
    logMessage(MSG_DETAIL, "Merging polygonal surface patches...\n");
    logMessage(MSG_DETAIL, " Began with %d patches.\n", iNPatches);
    // Create sets of patches incident on bdry verts.
    std::set<BdryPatch3D*> *aspBPIncident =
	new std::set<BdryPatch3D*>[iNumVerts()];
    for (int iP = 0; iP < iNumPatches(); iP++) {
      BdryPatch3D* pBP = pBdry(iP);
      for (int iV = 0; iV < pBP->iNumSegments(); iV++) {
	Vert *pV = pBP->pVVert(iV);
	int iV_ = iVertIndex(pV);
	aspBPIncident[iV_].insert(pBP);
      }
    }

    // Merge polygons ONLY at this point
    bool qGlobalChange;
    do {
      qGlobalChange = false;
      for (GR_index_t iV = 0; iV < iNumVerts(); iV++) {
	Vert *pV = pVVert(iV);
	// Create an alias for the set of patches incident on this vert.
	std::set<BdryPatch3D*> &spBP = aspBPIncident[iV];
	// Try to merge patches incident on this vert.  Keep going until
	// no more merges work.
	bool qLocalChange;
	do {
	  qLocalChange = false;
	  std::set<BdryPatch3D*>::iterator iter1, iter2;
	  for (iter1 = spBP.begin(); !spBP.empty() && iter1 != spBP.end();
	      iter1++) {
	    for (iter2 = iter1, iter2++; !spBP.empty() && iter2 != spBP.end();
		iter2++) {
	      BdryPatch3D *pBPatch1 = *iter1;
	      BdryPatch3D *pBPatch2 = *iter2;
	      // Are these both polygons?
	      BdryPolygon *pBPoly1 = dynamic_cast<BdryPolygon*>(pBPatch1);
	      BdryPolygon *pBPoly2 = dynamic_cast<BdryPolygon*>(pBPatch2);
	      if (!pBPoly1->qValid() || !pBPoly2->qValid()) {
		// One or both isn't a polygon.
		continue;
	      }

	      // Try the merge.
	      int iNewIntVerts = 0;
	      logMessage(MSG_DEBUG, "About to try a merge...\n");
	      logMessage(MSG_DEBUG, "   Polygon 1 verts:");
	      for (int ii = 0; ii < pBPoly1->iNumSegments(); ii++)
		logMessage(MSG_DEBUG, " %u", iVertIndex(pBPoly1->pVVert(ii)));
	      logMessage(MSG_DEBUG, "\n   Polygon 2 verts:");
	      for (int ii = 0; ii < pBPoly2->iNumSegments(); ii++)
		logMessage(MSG_DEBUG, " %u", iVertIndex(pBPoly2->pVVert(ii)));
	      logMessage(MSG_DEBUG, "\n");
	      if (pBPoly1->qMerge(*pBPoly2, pV, iNewIntVerts)) {
		qLocalChange = qGlobalChange = true;
		pBPoly2->vMarkDeleted();
		logMessage(MSG_DEBUG, "   Result verts:");
		for (int ii = 0; ii < pBPoly1->iNumSegments(); ii++)
		  logMessage(MSG_DEBUG, " %u", iVertIndex(pBPoly1->pVVert(ii)));
		logMessage(MSG_DEBUG, "\n");

		// Delete pBP2 from the BP sets of all of its verts.
		// Decrement iter2 -first-, because it's going to get
		// invalidated during deletion otherwise.  This will
		// have iter2 pointing at the last patch that -was-
		// checked already, which is perfect.
		iter2--;
		for (int iV2 = 0; iV2 < pBPoly2->iNumSegments(); iV2++) {
		  Vert *pVSet = pBPoly2->pVVert(iV2);
		  int iVSet = iVertIndex(pVSet);
		  int iRes = aspBPIncident[iVSet].erase(pBPoly2);
		  if (iRes != 1)
		    logMessage(
			MSG_DEBUG,
			"Odd... didn't seem to need to delete that one...\n");
		} // Done deleting merged patch from all incident patch
		  // sets

		  // Now make sure that pBPoly1 -is- in the lists for all
		  // of its verts.  This will mostly be a useless
		  // undertaking, but that's okay.
		for (int iV1 = 0; iV1 < pBPoly1->iNumSegments(); iV1++) {
		  Vert *pVSet = pBPoly1->pVVert(iV1);
		  int iVSet = iVertIndex(pVSet);
		  aspBPIncident[iVSet].insert(pBPoly1);
		} // Done adding new patch for all verts it contains.

		// Now make sure that verts newly interior to patches
		// have those (and all) patches removed from their patch
		// sets.
		for (int iV2 = pBPoly1->iNumInteriorVerts() - 1;
		    iV2 >= pBPoly1->iNumInteriorVerts() - iNewIntVerts; iV2--)
		  aspBPIncident[iVertIndex(pBPoly1->pVInteriorVert(iV2))].clear();
	      } // Done with post-merge stuff.
	    } // Done with inner loop over all patches incident on this vert
	  } // Done with outer loop over all patches incident on this vert
	}
	while (qLocalChange); // Done working on this vert
      } // Done with loop over all verts
    }
    while (qGlobalChange); // Done with all merging
    logMessage(MSG_SERVICE, "  Requiring %u triangles after merging.\n",
	       iNumTriangles());
    delete[] aspBPIncident;

    // Now remove all the deleted patches from the bdry rep.
    int iF = 0, iB = iNumPatches();
    while (iF < iB) {
      while ((iF < iB) && !apBPatches[iF]->qDeleted()) {
	iF++;
      }
      do {
	iB--;
      }
      while ((iF < iB) && apBPatches[iB]->qDeleted());
      if (iF < iB) {
	BdryPatch *pBPTmp = apBPatches[iF];
	apBPatches[iF] = apBPatches[iB];
	apBPatches[iB] = pBPTmp;
      }
    }
    iNPatches = iF;
    assert((iF == iB) || (iF == iB + 1));
    logMessage(MSG_SERVICE, " Finished with %d patches.\n", iNPatches);
    logMessage(MSG_SERVICE, "  Requiring %u triangles after purging.\n",
	       iNumTriangles());
  } // Done with level 1 optimizations
}

GR_index_t
Bdry3D::iNumTriangles() const
{
  int iRetVal = 0;
  for (int i = 0; i < iNPatches; i++) {
    BdryPatch3D *pBP = pBdry(i);
    if (!pBP->qDeleted())
      iRetVal += pBP->iNumTriangles();
  }
  return iRetVal;
}

int
Bdry3D::iNumPatchesContaining(const Vert* const pVS) const
{
  int iNumHits = 0;
  for (int i = 0; i < iNumPatches(); i++) {
    BdryPatch3D *pBP = pBdry(i);
    if (!pBP->qDisjointFrom(pVS))
      iNumHits++;
  }
  return iNumHits;
}

int
BdryPolygon::iNumTriangles() const
{
  return iNumSegments() - 2 + 2 * iNumChildren();
}

bool
BdryPolygon::qMerge(const BdryPolygon& BPoly, const Vert * const pV,
		    int &iNumNewInteriorVerts)
{
  // Do they have the same BC and region data?
  if ((qIsBoundary() != BPoly.qIsBoundary()) || (qIsBoundary() && // If one's a bdry, they both are.
      (iBdryCond() != BPoly.iBdryCond())) || qHasChildren()
      || BPoly.qHasChildren()) // Don't mess with polys
      // with holes
      {
    logMessage(MSG_DEBUG, "  FAILED: BC mismatch or poly with hole.\n");
    return false;
  }

  // Checking regions is slightly trickier, because the two
  // patches could be flipped relative to each other.
  bool qRev = false;
  if ((iLeftRegion() == BPoly.iLeftRegion())
      && (iRightRegion() == BPoly.iRightRegion())) {
    // This is okay, and the orientations are the same.  Do
    // nothing.
  }
  else if ((iLeftRegion() == BPoly.iRightRegion())
      && (iRightRegion() == BPoly.iLeftRegion())) {
    // This is okay, but the orientations are different.
    qRev = true;
  }
  else {
    // Otherwise, the regions don't match up; can't merge.
    logMessage(MSG_DEBUG, "  FAILED: region mismatch.\n");
    return false;
  }

  // Do they have the same normal?
  double adNorm1[3], adNorm2[3];
  vUnitNormal(pV->getCoords(), adNorm1);
  BPoly.vUnitNormal(pV->getCoords(), adNorm2);
  double dDot = dDOT3D(adNorm1, adNorm2) * (qRev ? -1 : 1);
  if (dDot < 0.999999) { //if (iFuzzyComp(dDot, 1.) != 0)
    logMessage(MSG_DEBUG, "  FAILED (dot): not coplanar: %f. %f %f\n", dDot,
	       dMAG3D(adNorm1), dMAG3D(adNorm2));
    return false;
  }

  // Now for the sensitive test of parallelism: the cross product should
  // be nearly zero.  When it isn't, bad things can happen.  Test case:
  // aircraftsidewall, which is only flat to within about 10^-6.
  double adCross[3];
  vCROSS3D(adNorm1, adNorm2, adCross);
  double dCrossMag = dMAG3D(adCross);
  if (iFuzzyComp(dCrossMag, 0.) != 0) {
    logMessage(MSG_DEBUG, "  FAILED (cross): not coplanar: %f. %20.15f\n",
	       dCrossMag, dDot);
    return false;
  }

  Vert **apVThis = new Vert*[iNumSegments() + BPoly.iNumSegments() - 2];
  Vert **apVMerge = new Vert*[BPoly.iNumSegments()];
  int iV;
  int iNSThis = iNumSegments();
  int iNSMerge = BPoly.iNumSegments();
  int iIndThis = -1, iIndMerge = -1;
  for (iV = 0; iV < iNSThis; iV++) {
    apVThis[iV] = pVVert(iV);
    if (apVThis[iV] == pV)
      iIndThis = iV;
  }
  // When copying the second set, flip them over as required.
  if (qRev) {
    for (iV = 0; iV < iNSMerge; iV++) {
      apVMerge[iV] = BPoly.pVVert(iNSMerge - iV - 1);
      if (apVMerge[iV] == pV)
	iIndMerge = iV;
    }
  }
  else {
    for (iV = 0; iV < iNSMerge; iV++) {
      apVMerge[iV] = BPoly.pVVert(iV);
      if (apVMerge[iV] == pV)
	iIndMerge = iV;
    }
  }
  assert(iIndThis != -1);
  assert(iIndMerge != -1);

  // Identify the string of common vertices in the two patches
  // (there has to be at least one, after all...).  If the
  // patches share a vert outside that run, forget it.
  int iFirstThis = iIndThis, iLastThis = iIndThis;
  int iFirstMerge = iIndMerge, iLastMerge = iIndMerge;
  // First check upwards in *this.
  while (apVThis[iLastThis] == apVMerge[iFirstMerge]) {
    iLastThis = (iLastThis + 1) % iNSThis;
    assert(iLastThis >= 0 && iLastThis < iNSThis);
    iFirstMerge = (iFirstMerge + iNSMerge - 1) % iNSMerge;
    assert(iFirstMerge >= 0 && iFirstMerge < iNSMerge);
  }
  iLastThis = (iLastThis + iNSThis - 1) % iNSThis;
  iFirstMerge = (iFirstMerge + 1) % iNSMerge;
  // Now check downwards in *this.
  while (apVThis[iFirstThis] == apVMerge[iLastMerge]) {
    iFirstThis = (iFirstThis + iNSThis - 1) % iNSThis;
    assert(iFirstThis >= 0 && iFirstThis < iNSThis);
    iLastMerge = (iLastMerge + 1) % iNSMerge;
    assert(iLastMerge >= 0 && iLastMerge < iNSMerge);
  }
  iLastMerge = (iLastMerge + iNSMerge - 1) % iNSMerge;
  iFirstThis = (iFirstThis + 1) % iNSThis;

  assert(
      (iLastThis - iFirstThis + iNSThis) % iNSThis
	  == (iLastMerge - iFirstMerge + iNSMerge) % iNSMerge);

  // Only one common vert; merge impossible
  if (iLastThis == iFirstThis) {
    logMessage(MSG_DEBUG, "  FAILED: only one common vert.\n");
    delete[] apVThis;
    delete[] apVMerge;
    return false;
  }

  for (iV = (iLastMerge + 1) % iNSMerge; iV != iFirstMerge;
      iV = (iV + 1) % iNSMerge) {
    Vert *pVCand = apVMerge[iV];
    // Check whether any of these are in the other list (anywhere)
    for (int iV2 = 0; iV2 < iNSThis; iV2++) {
      if (pVCand == apVThis[iV2]) {
	logMessage(MSG_DEBUG,
		   "  FAILED: won't create multiply connected patch.\n");
	delete[] apVThis;
	delete[] apVMerge;
	return false;
      }
    }
  }

  // Tag the verts that are in the interior of the new patch (those
  // between iFirstThis and iLastThis, not inclusive of the ends) for
  // later insertion interior to the patch.
  iNumNewInteriorVerts = 0;
  for (iV = (iFirstThis + 1) % iNSThis; iV != iLastThis;
      iV = (iV + 1) % iNSThis) {
    vAddInteriorVert(apVThis[iV]);
    iNumNewInteriorVerts++;
  }

  // The vertex list for the patch may get larger or smaller; determine
  // which, and by how much.  The first term is the number added from
  // BPoly, the second is the number of newly-minted interior verts.
  int iNumCommon = (iLastThis - iFirstThis + 1 + iNSThis) % iNSThis;
  assert(iNumCommon >= 2);
  int iNewSize = (iNSMerge - iNumCommon) + (iNSThis - iNumCommon) + 2;
  assert(iNewSize >= 3);
  assert(iNewSize <= iNumSegments() + BPoly.iNumSegments() - 2);
  logMessage(MSG_DEBUG,
	     " Polys w/ (%d, %d) verts; %d in common; result has %d verts.\n",
	     iNSThis, iNSMerge, iNumCommon, iNewSize);

  // Now update the list of perimeter verts in *this.
  if (iNewSize > iNSThis)
    vResize(iNewSize);
  // Start with the old vertices, running from iLastThis to iFirstThis,
  // inclusive of iLastThis but -not- iFirstThis.
  int iVStored = 0;
  for (iV = iLastThis; iV != iFirstThis; iV = (iV + 1) % iNSThis) {
    apVVerts[iVStored++] = apVThis[iV];
  }
  // Now get the verts from the neighbor patch, running from iLastMerge
  // to iFirstMerge, inclusive of iLastMerge but -not- iFirstMerge.
  for (iV = iLastMerge; iV != iFirstMerge; iV = (iV + 1) % iNSMerge) {
    apVVerts[iVStored++] = apVMerge[iV];
  }
  assert(apVMerge[iV] == apVVerts[0]);

  iNVerts = iNewSize;
  delete[] apVThis;
  delete[] apVMerge;
  return true;
}

void
BdryPatch3D::vAddInteriorVert(Vert* const pV)
{
  if (iNVertsToReinsert == iVRSize) {
    int iNewSize = max(10, iVRSize * 2);
    vResizeInteriorVerts(iNewSize);
  }
  apVVertsToReinsert[iNVertsToReinsert++] = pV;
}

