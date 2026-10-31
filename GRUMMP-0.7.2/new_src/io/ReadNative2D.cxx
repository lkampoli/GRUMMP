#include <cstdio>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include "GR_Mesh.h"
#include "GR_Mesh2D.h"

#define CHECK_READ(b, a) do {						\
		int _res_ = (a);							\
		lineNum++;								\
		if (_res_ != b) {							\
			logMessage(0, "File %s, line %d: \nSource line %d",			\
					strFileName, lineNum, __LINE__);				\
					vFatalError("File format error", __func__);			\
		}									\
} while(0)

#define SET_CELL_TO_FACE(cell_) do {					\
		if (cell_ >= 0) {							\
			int ind_ = cellFaces[cell_]++;					\
			if (ind_ >= 3 && nQuads == 0) {						\
				logMessage(0, "At input line %d: ", lineNum);			\
				vFatalError("Tri with four faces", __func__);	\
			}								\
			if (ind_ >= 4 && nQuads > 0) {						\
				logMessage(0, "At input line %d: ", lineNum);			\
				vFatalError("Quad with five faces", __func__);	\
			}								\
			cellToFace[cell_][ind_] = pF;			\
		}									\
} while(0)

#define CHECK_RANGE(val_, min_, max_) do {				\
		if (val_ < min_ || val_ > GR_sindex_t(max_)) {				\
			logMessage(0, "File %s, line %d: \n", strFileName, lineNum);	\
			logMessage(0, "Source line %d\n", __LINE__);			\
			vFatalError("Index range error", __func__);			\
		}									\
} while(0)

#define CHECK_POS_RANGE(val_, max_) do {					\
		if (val_ > max_) {							\
			logMessage(0, "File %s, line %d: \n", strFileName, lineNum);	\
			logMessage(0, "Source line %s:%d\n", __FILE__,__LINE__);		\
			logMessage(0, "  Index: %d, Max: %d\n", val_, max_);		\
			vFatalError("Index range error", __func__);			\
		}									\
} while(0)

#define CHECK_VERTEX_RANGE(val_) CHECK_POS_RANGE(val_, nVerts)
#define CHECK_FACE_RANGE(val_) CHECK_POS_RANGE(val_, nFaces)
#define CHECK_CELL_RANGE(val_) CHECK_POS_RANGE(val_, int(nCells))

void
readNativeGRMESH(Mesh2D& M2D, const char strFileName[])
{
  //@@ Read data from a GRUMMP native file and set up the data in the
  //Mesh2D data structure.
  FILE *pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile) {
    vFatalError("Couldn't open input file for reading", __func__);
  }

  // This file format does not support internal bdry faces.
  GR_index_t nVerts, nFaces, nCells, nBFaces, nIBFaces = 0;
  GR_index_t lineNum = 1;

  CHECK_READ(
      5,
      (fscanf(pFInFile, "%u %u %u %u %u\n", &nCells, &nFaces, &nBFaces,
	      &nIBFaces, &nVerts)));
  //IBFaces number can be zero.

  if (nVerts == 0 || nCells == 0 || nFaces == 0 || nBFaces == 0)
    vFatalError("Number of entities is specified to be 0!", __func__);

  if (2 * nFaces < nBFaces + 3 * nCells || 4 * nCells + nBFaces < 2 * nFaces)
    vFatalError("Impossible combination of entity counts!", __func__);

  GR_index_t nQuads = 2 * nFaces - 2 * nIBFaces - nBFaces - 3 * nCells;
  GR_index_t nTris = nCells - nQuads;
  assert(nQuads <= nCells);
  assert(nTris <= nCells);
  if (nQuads > 0)
    M2D.allowNonSimplicial();
  for (GR_index_t v = 0; v < nVerts; v++) {
    double dXX, dYY;
    CHECK_READ(2, fscanf(pFInFile, "%lf%lf\n", &dXX, &dYY));
    (void) M2D.createVert(dXX, dYY);
  }

  // Set up some temporary storage for tri and quad data
  Face *(**cellToFace) = new Face**[nCells];

//	Face* (*cellToFace)[4] = NULL;
  int *cellFaces = NULL;

  // For Mixed Element Meshes, just store the extra int
  // Could use an alternative structure like an array of lists
  // if memory read in 2D is a problem it'll be here.
  cellFaces = new int[nCells];

  for (GR_index_t iC = 0; iC < nCells; iC++)
    cellToFace[iC] = new Face*[(nQuads == 0) ? 3 : 4];

  for (GR_index_t t = 0; t < nCells; t++) {
    cellFaces[t] = 0;
    cellToFace[t][0] = cellToFace[t][1] = cellToFace[t][2] = pFInvalidFace;
  }
  if (nQuads > 0)
    for (GR_index_t t = 0; t < nCells; t++)
      cellToFace[t][3] = pFInvalidFace;

  // Read face-vert and face-cell connectivity; store cell-face
  // connectivity for later.
  for (GR_index_t f = 0; f < nFaces; f++) {
    GR_sindex_t cellA, cellB;
    GR_index_t vertA, vertB;
    CHECK_READ(4,
	       fscanf(pFInFile, "%d%d %u%u\n", &cellA, &cellB, &vertA, &vertB));
    CHECK_VERTEX_RANGE(vertA);
    CHECK_VERTEX_RANGE(vertB);
    CHECK_CELL_RANGE(cellA);
    CHECK_CELL_RANGE(cellB);

    bool qExist;
    Face *pF = M2D.createFace(qExist, M2D.getVert(vertA), M2D.getVert(vertB));
    if (qExist) {
      M2D.createFace(qExist, M2D.getVert(vertA), M2D.getVert(vertB), true);
    }
    if (nIBFaces == 0) {
      assert(!qExist);
    }
    if (!qExist) {
      // These next two are macros
      SET_CELL_TO_FACE(cellA);
      SET_CELL_TO_FACE(cellB);
    }
  }

  // We'll create cells in a moment; right now, we don't have region
  // data yet.

  // First, read bdry face data.  Don't create bdry faces until after
  // cells already exist, though...

  int *bfaceToFace = new int[nBFaces];
  int *bdryConds = new int[nBFaces];

  for (GR_index_t bf = 0; bf < nBFaces; bf++) {
    /*bdryfaces: face bc verts*/
    GR_index_t face, vertA, vertB;
    int iBC;
    // Read and check data
    CHECK_READ(4,
	       fscanf(pFInFile, "%u %d %u%u\n", &face, &iBC, &vertA, &vertB));
    CHECK_VERTEX_RANGE(vertA);
    CHECK_VERTEX_RANGE(vertB);
    CHECK_FACE_RANGE(face);

    bfaceToFace[bf] = face;
    bdryConds[bf] = iBC;
  }
  int * ibfaceToFaceA;
  int * ibfaceToFaceB;
  if (nIBFaces > 0) {
    ibfaceToFaceA = new int[nIBFaces];
    ibfaceToFaceB = new int[nIBFaces];

    for (GR_index_t ibf = 0; ibf < nIBFaces; ibf++) {
      /*intbdryfaces: faces verts*/
      GR_index_t faceA, faceB, vertA, vertB;
      // Read and check data
      CHECK_READ(
	  4, fscanf(pFInFile, "%u %u %u %u\n", &faceA, &faceB, &vertA, &vertB));
      CHECK_VERTEX_RANGE(vertA);
      CHECK_VERTEX_RANGE(vertB);
      CHECK_FACE_RANGE(faceA);
      CHECK_FACE_RANGE(faceB);
      ibfaceToFaceA[ibf] = faceA;
      ibfaceToFaceB[ibf] = faceB;
    }
  }
  // Now read region data and create cells
  for (GR_index_t c = 0; c < nCells; c++) {
    int iReg = 1;
    /*cells: region*/
    CHECK_READ(1, fscanf(pFInFile, "%d\n", &iReg));
    if (cellFaces[c] == 3)
      (void) M2D.createTriCell(cellToFace[c][0], cellToFace[c][1],
			       cellToFace[c][2], iReg);
    if (cellFaces[c] == 4)
      (void) M2D.createQuadCell(cellToFace[c][0], cellToFace[c][1],
				cellToFace[c][2], cellToFace[c][3], iReg);
    if (cellFaces[c] < 3 || cellFaces[c] > 4)
      vFatalError("Cell has invalid number faces", __func__);

  }

  for (GR_index_t bf = 0; bf < nBFaces; bf++) {
    (void) M2D.createBFace(M2D.getFace(bfaceToFace[bf]), bdryConds[bf]);
  }

  if (nIBFaces > 0) {
    for (GR_index_t ibf = 0; ibf < nIBFaces; ibf++) {
      Face * fA = M2D.getFace(ibfaceToFaceA[ibf]);
      Face * fB = M2D.getFace(ibfaceToFaceB[ibf]);
      if (fB->getNumCells() == 0) {

	Cell * RightCell = fA->getRightCell();
	RightCell->replaceFace(fA, fB);
	fB->setRightCell(fA->getRightCell());
	fA->setRightCell(NULL);
      }
      else {

	Cell * RightCell = fB->getRightCell();
	RightCell->replaceFace(fB, fA);
	fA->setRightCell(fB->getRightCell());
	fB->setRightCell(NULL);
      }
      (void) M2D.createIntBdryEdge(fA, fB, pBFInvalidBFace);
    }
  }

  assert(M2D.iNumVerts() == nVerts);
  assert(M2D.getNumFaces() == nFaces);
  assert(M2D.getNumCells() == nCells);
  assert(M2D.getNumBdryFaces() == nBFaces);
  assert(M2D.getNumIntBdryFaces() == nIBFaces);
  delete[] bfaceToFace;
  delete[] bdryConds;
  if (nIBFaces > 0) {
    delete[] ibfaceToFaceA;
    delete[] ibfaceToFaceB;
  }

  assert(cellFaces && cellToFace);
  delete[] cellFaces;
  for (GR_index_t iC = 0; iC < nCells; iC++)
    delete[] cellToFace[iC];
  delete[] cellToFace;

  fclose(pFInFile);
  if (!M2D.isValid())
    vFatalError("Mesh was invalid.", "2D native mesh read");
}
void readNativeMESH(Mesh2D& M2D, const char strFileName[]) {
	//@@ Read data from a GRUMMP native file and set up the data in the
	//Mesh2D data structure.
	FILE *pFInFile = fopen(strFileName, "r");
	if (NULL == pFInFile) {
		vFatalError("Couldn't open input file for reading", __func__);
	}

	GR_index_t nVerts, nFaces, nCells, nBFaces;
	GR_index_t lineNum = 1;

	CHECK_READ(4,
			(fscanf(pFInFile, "%u %u %u %u\n", &nCells, &nFaces, &nBFaces,
					&nVerts)));
	//IBFaces number can be zero.

	if (nVerts == 0 || nCells == 0 || nFaces == 0 || nBFaces == 0)
		vFatalError("Number of entities is specified to be 0!", __func__);

	if (2 * nFaces < nBFaces + 3 * nCells || 4 * nCells + nBFaces < 2 * nFaces)
		vFatalError("Impossible combination of entity counts!", __func__);

	GR_index_t nQuads = 2 * nFaces - nBFaces - 3 * nCells;
	GR_index_t nTris = nCells - nQuads;
	assert(nQuads <= nCells);
	assert(nTris <= nCells);
	if (nQuads > 0)
		M2D.allowNonSimplicial();

	int res = -1;
	size_t bufSize = 1024;
	char* buffer = reinterpret_cast<char*>(calloc(bufSize, 1));
	Vert *pV;
	for (GR_index_t v = 0; v < nVerts; v++) {
		double dXX, dYY;
		getline(&buffer, &bufSize, pFInFile);
		int vertType = Vert::eUnknown;
		res = sscanf(buffer, "%lf%lf%d", &dXX, &dYY, &vertType);

		assert(res == 2 || res == 3);
		pV = M2D.createVert(dXX, dYY);
		if (res == 3) {
		//	cout << pV->getVertType() << "  ";
			pV->setType(vertType);
		//	cout << pV->getVertType() << endl;
		}
	}

  free (buffer);

  // Set up some temporary storage for tri and quad data
  Face *(**cellToFace) = new Face**[nCells];

//	Face* (*cellToFace)[4] = NULL;
  int *cellFaces = NULL;

  // For Mixed Element Meshes, just store the extra int
  // Could use an alternative structure like an array of lists
  // if memory read in 2D is a problem it'll be here.
  cellFaces = new int[nCells];

  for (GR_index_t iC = 0; iC < nCells; iC++)
    cellToFace[iC] = new Face*[(nQuads == 0) ? 3 : 4];

  for (GR_index_t t = 0; t < nCells; t++) {
    cellFaces[t] = 0;
    cellToFace[t][0] = cellToFace[t][1] = cellToFace[t][2] = pFInvalidFace;
  }
  if (nQuads > 0)
    for (GR_index_t t = 0; t < nCells; t++)
      cellToFace[t][3] = pFInvalidFace;

  // Read face-vert and face-cell connectivity; store cell-face
  // connectivity for later.
  for (GR_index_t f = 0; f < nFaces; f++) {
    GR_sindex_t cellA, cellB;
    GR_index_t vertA, vertB;
    CHECK_READ(4,
	       fscanf(pFInFile, "%d%d %u%u\n", &cellA, &cellB, &vertA, &vertB));
    CHECK_VERTEX_RANGE(vertA);
    CHECK_VERTEX_RANGE(vertB);
    CHECK_CELL_RANGE(cellA);
    CHECK_CELL_RANGE(cellB);

    bool qExist;
    Face *pF = M2D.createFace(qExist, M2D.getVert(vertA), M2D.getVert(vertB));
    if (qExist) {
      M2D.createFace(qExist, M2D.getVert(vertA), M2D.getVert(vertB), true);
    }
    assert(!qExist);
    if (!qExist) {
      // These next two are macros
      SET_CELL_TO_FACE(cellA);
      SET_CELL_TO_FACE(cellB);
    }
  }

  // We'll create cells in a moment; right now, we don't have region
  // data yet.

  // First, read bdry face data.  Don't create bdry faces until after
  // cells already exist, though...

  int *bfaceToFace = new int[nBFaces];
  int *bdryConds = new int[nBFaces];

  for (GR_index_t bf = 0; bf < nBFaces; bf++) {
    /*bdryfaces: face bc verts*/
    GR_index_t face, vertA, vertB;
    int iBC;
    // Read and check data
    CHECK_READ(4,
	       fscanf(pFInFile, "%u %d %u%u\n", &face, &iBC, &vertA, &vertB));
    CHECK_VERTEX_RANGE(vertA);
    CHECK_VERTEX_RANGE(vertB);
    CHECK_FACE_RANGE(face);

    bfaceToFace[bf] = face;
    bdryConds[bf] = iBC;
  }

  // Now read region data and create cells
  for (GR_index_t c = 0; c < nCells; c++) {
    int iReg = 1;
    /*cells: region*/
    CHECK_READ(1, fscanf(pFInFile, "%d\n", &iReg));
    if (cellFaces[c] == 3) {
//			(void) M2D.createTriCell(cellToFace[c][0], cellToFace[c][1], cellToFace[c][2], iReg);

      assert(cellToFace[c][0]->isValid());
      assert(cellToFace[c][1]->isValid());
      assert(cellToFace[c][2]->isValid());
      Vert *apV[3];
      apV[0] = findCommonVert(cellToFace[c][0], cellToFace[c][2]);
      apV[1] = findCommonVert(cellToFace[c][0], cellToFace[c][1]);
      apV[2] = findCommonVert(cellToFace[c][1], cellToFace[c][2]);
      (void) M2D.createTriCell(apV[0], apV[1], apV[2], iReg, true);
    }
    if (cellFaces[c] == 4)
      (void) M2D.createQuadCell(cellToFace[c][0], cellToFace[c][1],
				cellToFace[c][2], cellToFace[c][3], iReg);
    if (cellFaces[c] < 3 || cellFaces[c] > 4)
      vFatalError("Cell has invalid number faces", __func__);

  }

  for (GR_index_t bf = 0; bf < nBFaces; bf++) {
    (void) M2D.createBFace(M2D.getFace(bfaceToFace[bf]), bdryConds[bf]);
  }

  assert(M2D.iNumVerts() == nVerts);
  assert(M2D.getNumFaces() == nFaces);
  assert(M2D.getNumCells() == nCells);
  assert(M2D.getNumBdryFaces() == nBFaces);
  delete[] bfaceToFace;
  delete[] bdryConds;

  assert(cellFaces && cellToFace);
  delete[] cellFaces;
  for (GR_index_t iC = 0; iC < nCells; iC++)
    delete[] cellToFace[iC];
  delete[] cellToFace;

  fclose(pFInFile);
  if (!M2D.isValid())
    vFatalError("Mesh was invalid.", "2D native mesh read");
}

void
readTriangle(Mesh2D& M2D, const char strFileName[])
{
  char strBaseFileName[FILE_NAME_LEN], strNodeFileName[FILE_NAME_LEN],
      strEleFileName[FILE_NAME_LEN], strEdgeFileName[FILE_NAME_LEN];
//  strNeighFileName[FILE_NAME_LEN];

  const char *fileNameEnd;
  fileNameEnd = strstr(strFileName, ".poly");
  if (fileNameEnd != NULL) {
    int len = fileNameEnd - strFileName;

    if (len > 1023)
      len = 1023;
    strncpy(strBaseFileName, strFileName, len);
    strBaseFileName[len] = 0;
  }
  // there are a lot of files in here
  // strFileName contains .poly
  // rest are obvious
  sprintf(strNodeFileName, "%s.node", strBaseFileName);
  sprintf(strEleFileName, "%s.ele", strBaseFileName);
  sprintf(strEdgeFileName, "%s.edge", strBaseFileName);

  // read poly file first, then just open them all
  FILE *pFInFile = fopen(strFileName, "r");
  FILE *pFNodeFile = fopen(strNodeFileName, "r");
  FILE *pFEleFile = fopen(strEleFileName, "r");
  FILE *pFEdgeFile = fopen(strEdgeFileName, "r");

  if (NULL == pFInFile) {
    vFatalError("Couldn't open poly file for reading", __func__);
  }

  GR_index_t nVerts, nFaces, nTris, nBFaces, nEleAttr, nNodeAttr, dim,
      bdryMarkers;
  GR_index_t lineNum = 1;

  CHECK_READ(
      4,
      (fscanf(pFInFile, "%u %u %u %u\n", &nVerts, &dim, &nNodeAttr,
	      &bdryMarkers)));

  CHECK_READ(2, (fscanf(pFInFile, "%u %u\n", &nBFaces, &bdryMarkers)));

  if (dim != 2)
    vFatalError("invalid dimension", __func__);

  if (nVerts == 0) { // we must have a node file
    if (NULL == pFNodeFile) {
      vFatalError("Couldn't open node file for reading", __func__);
    }
    CHECK_READ(
	4,
	(fscanf(pFNodeFile, "%u %u %u %u\n", &nVerts, &dim, &nNodeAttr,
		&bdryMarkers)));
    if (dim != 2)
      vFatalError("invalid dimension", __func__);
  }

  // get number of triangles
  CHECK_READ(3, (fscanf(pFEleFile, "%u %u %u\n", &nTris, &dim, &nEleAttr)));
  // get number of edges
  CHECK_READ(2, (fscanf(pFEdgeFile, "%u %u\n", &nFaces, &dim)));

  if (nVerts == 0 || nTris == 0 || nFaces == 0 || nBFaces == 0)
    vFatalError("Number of entities is specified to be 0!", __func__);

  // storage for Verts
  Vert ** newVerts = new Vert*[nVerts];
  GR_index_t vertNum, vertBC;
  for (GR_index_t v = 0; v < nVerts; v++) {
    double dXX, dYY, tempAtt;
    CHECK_READ(3, fscanf(pFNodeFile, "%u %lf %lf ", &vertNum, &dXX, &dYY));
    for (GR_index_t iAtt = 0; iAtt < nNodeAttr; iAtt++)
      CHECK_READ(1, fscanf(pFNodeFile, "%lf ", &tempAtt));
    fscanf(pFNodeFile, "%d\n", &vertBC);

    dXX = round(dXX * 1.e7) / 1.e7;
    dYY = round(dYY * 1.e7) / 1.e7;
    newVerts[v] = M2D.createVert(dXX, dYY);
  }
  if (vertNum != nVerts)
    vFatalError("incorrect number of verts", __func__);

  std::map<Face*, int> faceBC;
  GR_index_t faceNum;
  GR_index_t numBdrys = 0;
  // Read face-vert connectivity, collect BCs
  for (GR_index_t f = 0; f < nFaces; f++) {
    GR_index_t vertA, vertB, BC;
    CHECK_READ(
	4, fscanf(pFEdgeFile, "%u %u %u %u\n", &faceNum, &vertA, &vertB, &BC));
    CHECK_VERTEX_RANGE(vertA - 1);
    CHECK_VERTEX_RANGE(vertB - 1);

    bool qExist;
    Face *pF = M2D.createFace(qExist, newVerts[vertA - 1], newVerts[vertB - 1]);
    if (BC > 0) {
      numBdrys++;
      faceBC.insert(std::make_pair(pF, BC));
    }
  }
//	printf("%d %d %d %d\n",nVerts,nTris,nFaces,nBFaces);

  if (numBdrys != nBFaces)
    nBFaces = numBdrys;
//	printf("%d %d %d %d\n",nVerts,nTris,nFaces,nBFaces);
//	printf("%d Verts %d Tris %d Faces %d BFaces\n",nVerts,nTris,nFaces,nBFaces);
//printf("%d < %d || %d < %d\n",2*nFaces, nBFaces + 3*nTris, 4*nTris + nBFaces, 2*nFaces);
  if (2 * nFaces < nBFaces + 3 * nTris || 3 * nTris + nBFaces < 2 * nFaces)
    vFatalError("Impossible combination of entity counts!", __func__);

  // We'll create cells in a moment; right now, we don't have region
  // data yet.

  // Now read region data and create cells
  int iReg = 1;
  for (GR_index_t c = 0; c < nTris; c++) {
    GR_index_t vertA, vertB, vertC;
    CHECK_READ(
	4,
	fscanf(pFEleFile, "%u %u %u %u\n", &vertNum, &vertA, &vertB, &vertC));
    (void) M2D.createTriCell(newVerts[vertA - 1], newVerts[vertB - 1],
			     newVerts[vertC - 1], iReg);
  }
// now deal with boundary faces
  for (std::map<Face*, int>::iterator it = faceBC.begin(); it != faceBC.end();
      ++it) {
    M2D.createBFace(it->first, it->second);
  }

  assert(M2D.iNumVerts() == nVerts);
  assert(M2D.getNumFaces() == nFaces);
  assert(M2D.getNumCells() == nTris);
  assert(M2D.getNumBdryFaces() == nBFaces);
  assert(M2D.getNumIntBdryFaces() == 0);

  fclose(pFInFile);
  if (pFNodeFile)
    fclose(pFNodeFile);
  fclose(pFEleFile);
  fclose(pFEdgeFile);

  if (!M2D.isValid())
    vFatalError("Mesh was invalid.", "readTriangle");
}

