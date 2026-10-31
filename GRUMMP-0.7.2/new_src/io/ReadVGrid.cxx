#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "GR_misc.h"
#include "GR_VolMesh.h"

// The following functions are to convert from little-endian to
// big-endian format.

static bool qFlipEndian = false;

static int
convertToInt(const unsigned char raw[4])
{
  // This implementation should be okay even with 64 bit ints, whereas a
  // union might be a little tricky.
//   return ( (reinterpret_cast<int>(raw[0]) << 24)
// 	   + (reinterpret_cast<int>(raw[1]) << 16)
// 	   + (reinterpret_cast<int>(raw[2]) << 8)
// 	   + (raw[3]) );
  static const int is = sizeof(int);
  union {
    unsigned char raw2[is];
    int dummy;
  } a;
  a.raw2[is - 1] = raw[0];
  a.raw2[is - 2] = raw[1];
  a.raw2[is - 3] = raw[2];
  a.raw2[is - 4] = raw[3];
  return a.dummy;
}

static double
convertToDouble(const unsigned char raw[8])
{
  union {
    unsigned char raw2[8];
    double dummy;
  } a;
  a.raw2[7] = raw[0];
  a.raw2[6] = raw[1];
  a.raw2[5] = raw[2];
  a.raw2[4] = raw[3];
  a.raw2[3] = raw[4];
  a.raw2[2] = raw[5];
  a.raw2[1] = raw[6];
  a.raw2[0] = raw[7];
  /*   printf("%02x%02x%02x%02x%02x%02x%02x%02x\n", */
  /* 	 raw[0], raw[1], raw[2], raw[3], */
  /* 	 raw[4], raw[5], raw[6], raw[7]); */
  return a.dummy;
}

static int
readInt(FILE * const fp)
{
  int input = -1;
  if (qFlipEndian) {
    // Just read it
    GR_index_t nRead = fread(&input, sizeof(int), 1, fp);
    assert(nRead == 1);
  }
  else {
    // Read raw data and convert
    unsigned char raw[4];
    GR_index_t nRead = fread(raw, 1, 4, fp);
    input = convertToInt(raw);
    assert(nRead == 4);
  }
  return input;
}

static double
readDouble(FILE * const fp)
{
  double input = 0;
  if (qFlipEndian) {
    // Just read it
    GR_index_t nRead = fread(&input, sizeof(double), 1, fp);
    assert(nRead == 1);
  }
  else {
    // Read raw data and convert
    unsigned char raw[8];
    GR_index_t nRead = fread(raw, 1, 8, fp);
    input = convertToDouble(raw);
    assert(nRead == 8);
  }
  return input;
}

// VGrid writes unformatted Fortran files, containing only tets.  The
// info in the .cogsg file is very limited; BC info, for instance, is in
// a separate file, which we don't read at present.
void
readVGridBinary(const char * const fileName, VolMesh * const pVM)
{
  GR_index_t i, j;
//   GR_index_t minVert = 1000000, maxVert = 0;

  GR_index_t (*tetToVert)[4] = NULL;
//   GR_index_t (*triToVert)[3] = NULL;

  double (*x)[3];

  FILE *fp = fopen(fileName, "r");
  if (fp == NULL)
    vFatalError("Could not open file", "reading VGrid file");

  GR_index_t lineLength = readInt(fp);
  GR_index_t junk = readInt(fp);
  assert(++junk); // Quiet warnings.
  GR_index_t nTet = readInt(fp);
  GR_index_t nVert = readInt(fp);
  GR_index_t nBdryPts = readInt(fp);
  GR_index_t nLayerPts = readInt(fp);
  GR_index_t nLayerTets = readInt(fp);
  junk = readInt(fp);
  junk = readInt(fp);

  GR_index_t nPredictedLength = (8 + 4 * nTet) * 4; // That last 4 is sizeof(int)...
  if (lineLength != nPredictedLength) {
    if (qFlipEndian == true) {
      vFatalError("File is neither big nor little endian",
		  "reading VGrid file");
    }
    qFlipEndian = true;
    fclose(fp);
    readVGridBinary(fileName, pVM);
    return;
  }

  logMessage(2, "%8u verts\n", nVert);
  logMessage(2, "%8u bdry verts\n", nBdryPts);
  logMessage(2, "%8u tets\n", nTet);
  logMessage(2, "%8u verts in the bdry layer\n", nLayerPts);
  logMessage(2, "%8u tets in the bdry layer\n", nLayerTets);

  {
    // Estimate GRUMMP data struct sizes; the extra word in cell sizes
    // is a min for the allocated face ptr array overhead.
    GR_index_t iTetSize = sizeof(TetCell) + 5 * sizeof(void*);
    GR_index_t iVertSize = sizeof(Vert) + 44 * sizeof(void*); // Face connectivity
    GR_index_t iFaceSize = sizeof(Face);

    GR_index_t iVertMem = nVert * iVertSize / 1024;
    GR_index_t iTetMem = nTet * iTetSize / 1024;

    GR_index_t nTriFace = 2 * nTet;
    GR_index_t iTriFaceMem = nTriFace * iFaceSize / 1024;

    logMessage(2, "Expecting %8u kB of memory for verts\n", iVertMem);
    logMessage(2, "Expecting %8u kB of memory for tets\n", iTetMem);
    logMessage(2, "Expecting %8u kB of memory for tri faces\n", iTriFaceMem);

    logMessage(2, "Total memory approx %u MB + overhead\n",
	       (iVertMem + iTriFaceMem) / 1024);
  }

  /* Now read connectivity. */

  logMessage(2, "Reading tet connectivity\n");

  tetToVert = new GR_index_t[nTet][4];

  for (j = 0; j < 4; j++) {
    for (i = 0; i < nTet; i++) {
      tetToVert[i][j] = readInt(fp) - 1;
//       minVert = (tetToVert[i][j] < minVert) ? tetToVert[i][j] : minVert;
//       maxVert = (tetToVert[i][j] > maxVert) ? tetToVert[i][j] : maxVert;
    }
  }

  // minVert should always be 1; this is a Fortran file we're reading.
//   vMessage(2, "Min vert index = %d; max = %d\n", minVert, maxVert);

  GR_index_t lineCheck = readInt(fp);
  assert(lineCheck == lineLength);

  logMessage(2, "\nReading coordinates\n");

  x = new double[nVert][3];

  lineLength = readInt(fp);
  for (j = 0; j < 3; j++) {
    for (i = 0; i < nVert; i++) {
      x[i][j] = readDouble(fp);
    }
  }
  lineCheck = readInt(fp);
  assert(lineCheck == lineLength);

  for (GR_index_t iV = 0; iV < nVert; iV++) {
    (void) pVM->createVert(x[iV]);
  }
  logMessage(2, "Set coords of %u vertices\n", nVert);
  delete[] x;

#ifdef BUILD_TETS_SLOWLY
  for (GR_index_t iT = 0; iT < nTet; iT++) {
    Vert *pV0 = pVM->getVert(tetToVert[iT][0]);
    Vert *pV1 = pVM->getVert(tetToVert[iT][1]);
    Vert *pV2 = pVM->getVert(tetToVert[iT][2]);
    Vert *pV3 = pVM->getVert(tetToVert[iT][3]);

    pVM->createTetCell(pV0, pV1, pV2, pV3, 1);
    if ((iT + 1) % (nTet / 10 + 1) == 0) {
      logMessage(2, "Create %8u tets out of %u\n", iT + 1, nTet);
    }
  }
  logMessage(2, "Created %8u tets\n", nTet);
  delete [] tetToVert;

  logMessage(2, "%8d tri faces.\n",
      pVM->iNumTriFaces());

  // Iterate over all triangles; those with only one tet are bdry faces,
  // so add them.

  GR_index_t nTri = pVM->iNumTriFaces();
  GR_index_t nBdryTri = 2*nBdryPts - 4;
  for (GR_index_t iF = 0; iF < nTri; iF++) {
    Face *pF = pVM->getFace(iF);
    if (pF->getLeftCell()->isValid() && pF->getRightCell()->isValid())
    continue;

    BFace *pBF = pVM->createBFace(pF);

    if ((iF + 1) % (nTri / 10 + 1) == 0) {
      logMessage(2, "Creating bdry tris: %8d tris of %d; %6d bdry tris of %d\n",
	  iF + 1, nTri, pVM->iNumBdryFaces(), nBdryTri);
    }
    delete [] triToVert;
  }
  logMessage(2, "Created bdry tris: %6d bdry tris of %d\n",
      pVM->iNumBdryFaces(), nBdryTri);
#else
  // Use vConvertFromCellVert instead.
  GR_index_t bdryTriToVert[1][3];
  GR_index_t *regions = new GR_index_t[nTet];
  for (i = 0; i < nTet; i++)
    regions[i] = iDefaultRegion;
  pVM->convertFromCellVert(nTet, 0, tetToVert, regions, bdryTriToVert);
#endif

  GR_index_t iNumOK = 0;
  for (GR_index_t iF = 0; iF < pVM->getNumFaces(); iF++) {
    Face *pF = pVM->getFace(iF);
    GR_index_t iOK = pF->getLeftCell()->isValid()
	&& pF->getRightCell()->isValid();
    iNumOK += iOK;
    if (iOK == 0) {
      logMessage(3, "Problem with face: verts %8u %8u %8u (face %u)\n",
		 pVM->getVertIndex(pF->getVert(0)),
		 pVM->getVertIndex(pF->getVert(1)),
		 pVM->getVertIndex(pF->getVert(2)), iF);
    }
  }
  logMessage(2, "%u faces okay out of %u\n", iNumOK, pVM->getNumFaces());

  assert(pVM->isValid());
  logMessage(2, "Done reading VGrid file\n");
}
