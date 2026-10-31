#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "GR_misc.h"
#include "GR_VolMesh.h"

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

void
readUGridBinary(const char * const fileName, VolMesh * const pVM,
		int*& bdryCond)
{
  GR_index_t nTet, nPyr, nPrism, nHex;
  GR_index_t nVert, nBdryTri, nBdryQuad;
  GR_index_t i, j;
  GR_index_t minVert = 1000000, maxVert = 0;

  GR_index_t (*tetToVert)[4] = NULL;
  GR_index_t (*pyrToVert)[5] = NULL;
  GR_index_t (*prismToVert)[6] = NULL;
  GR_index_t (*hexToVert)[8] = NULL;

  GR_index_t (*triToVert)[3] = NULL;
  GR_index_t (*quadToVert)[4] = NULL;

  unsigned char raw[36];

  double (*x)[3];

  FILE *fp = fopen(fileName, "r");
  if (!fp) {
    vFatalError("Missing ugrid file!", __func__);
  }

  GR_index_t nRead = fread(raw, 1, 28, fp);
  assert(nRead == 28);

  nVert = convertToInt(raw);
  nBdryTri = convertToInt(raw + 4);
  nBdryQuad = convertToInt(raw + 8);
  nTet = convertToInt(raw + 12);
  nPyr = convertToInt(raw + 16);
  nPrism = convertToInt(raw + 20);
  nHex = convertToInt(raw + 24);

  logMessage(2, "%8u verts\n", nVert);
  logMessage(2, "%8u bdry tris\n", nBdryTri);
  logMessage(2, "%8u bdry quads\n", nBdryQuad);
  logMessage(2, "%8u tets\n", nTet);
  logMessage(2, "%8u pyrs\n", nPyr);
  logMessage(2, "%8u prisms\n", nPrism);
  logMessage(2, "%8u hexes\n\n", nHex);

  {
    // Estimate GRUMMP data struct sizes; the extra word in cell sizes
    // is a min for the allocated face ptr array overhead.
    GR_index_t iTetSize = sizeof(TetCell) + 5 * sizeof(void*);
    GR_index_t iPyrSize = sizeof(Cell) + 6 * sizeof(void*);
    GR_index_t iPrismSize = sizeof(Cell) + 7 * sizeof(void*);
    GR_index_t iHexSize = sizeof(Cell) + 9 * sizeof(void*);
    GR_index_t iVertSize = sizeof(Vert) + 35 * sizeof(void*); // Face connectivity
    GR_index_t iFaceSize = sizeof(Face);

    GR_index_t iVertMem = nVert * iVertSize / 1024;
    GR_index_t iTetMem = nTet * iTetSize / 1024;
    GR_index_t iPyrMem = nPyr * iPyrSize / 1024;
    GR_index_t iPrismMem = nPrism * iPrismSize / 1024;
    GR_index_t iHexMem = nHex * iHexSize / 1024;

    GR_index_t nTriFace = (4 * (nTet + nPyr) + 2 * nPrism) / 2 - nBdryTri;
    GR_index_t nQuadFace = (3 * nPrism + nPyr + 6 * nHex) / 2 - nBdryQuad;

    GR_index_t iTriFaceMem = nTriFace * iFaceSize / 1024;
    GR_index_t iQuadFaceMem = nQuadFace * iFaceSize / 1024;

    logMessage(2, "Expecting %8u kB of memory for verts\n", iVertMem);
    logMessage(2, "Expecting %8u kB of memory for tets\n", iTetMem);
    logMessage(2, "Expecting %8u kB of memory for pyrs\n", iPyrMem);
    logMessage(2, "Expecting %8u kB of memory for prisms\n", iPrismMem);
    logMessage(2, "Expecting %8u kB of memory for hexes\n", iHexMem);
    logMessage(2, "Expecting %8u kB of memory for tri faces\n", iTriFaceMem);
    logMessage(2, "Expecting %8u kB of memory for quad faces\n", iQuadFaceMem);

    logMessage(
	2,
	"Total memory approx %u MB + overhead\n",
	(iVertMem + iTetMem + iPyrMem + iPrismMem + iTriFaceMem + iQuadFaceMem)
	    / 1024);
  }

  logMessage(2, "\nReading coordinates\n");

  x = new double[nVert][3];

  for (i = 0; i < nVert; i++) {
    nRead = fread(raw, 1, 24, fp);
    x[i][0] = convertToDouble(raw);
    x[i][1] = convertToDouble(raw + 8);
    x[i][2] = convertToDouble(raw + 16);
  }

  /* Now read connectivity. */

  if (nBdryTri > 0) {
    logMessage(2, "Reading bdry tri connectivity\n");
    triToVert = new GR_index_t[nBdryTri][3];
    for (i = 0; i < nBdryTri; i++) {
      for (j = 0; j < 3; j++) {
	nRead = fread(raw, 1, 4, fp);
	triToVert[i][j] = convertToInt(raw);
      }
    }
  }

  if (nBdryQuad > 0) {
    logMessage(2, "Reading bdry quad connectivity\n");
    quadToVert = new GR_index_t[nBdryQuad][4];
    for (i = 0; i < nBdryQuad; i++) {
      for (j = 0; j < 4; j++) {
	nRead = fread(raw, 1, 4, fp);
	quadToVert[i][j] = convertToInt(raw);
      }
    }
  }

  // Read bdry condition info and store it.
  if (bdryCond) {
    delete[] bdryCond;
  }
  assert(nBdryTri + nBdryQuad > 0);
  bdryCond = new int[nBdryTri + nBdryQuad];
  for (i = 0; i < nBdryTri + nBdryQuad; i++) {
    nRead = fread(raw, 1, 4, fp);
    bdryCond[i] = convertToInt(raw);
  }

  if (nTet > 0) {
    logMessage(2, "Reading tet connectivity\n");

    tetToVert = new GR_index_t[nTet][4];

    for (i = 0; i < nTet; i++) {
      for (j = 0; j < 4; j++) {
	nRead = fread(raw, 1, 4, fp);
	tetToVert[i][j] = convertToInt(raw);
	minVert = (tetToVert[i][j] < minVert) ? tetToVert[i][j] : minVert;
	maxVert = (tetToVert[i][j] > maxVert) ? tetToVert[i][j] : maxVert;
      }
    }
  }

  if (nPyr > 0) {
    logMessage(2, "Reading pyramid connectivity\n");

    pyrToVert = new GR_index_t[nPyr][5];

    for (i = 0; i < nPyr; i++) {
      for (j = 0; j < 5; j++) {
	nRead = fread(raw, 1, 4, fp);
	pyrToVert[i][j] = convertToInt(raw);
	minVert = (pyrToVert[i][j] < minVert) ? pyrToVert[i][j] : minVert;
	maxVert = (pyrToVert[i][j] > maxVert) ? pyrToVert[i][j] : maxVert;
      }
    }
  }

  if (nPrism > 0) {
    logMessage(2, "Reading prism connectivity\n");

    prismToVert = new GR_index_t[nPrism][6];

    for (i = 0; i < nPrism; i++) {
      for (j = 0; j < 6; j++) {
	nRead = fread(raw, 1, 4, fp);
	prismToVert[i][j] = convertToInt(raw);
	minVert = (prismToVert[i][j] < minVert) ? prismToVert[i][j] : minVert;
	maxVert = (prismToVert[i][j] > maxVert) ? prismToVert[i][j] : maxVert;
      }
    }
  }

  if (nHex > 0) {
    logMessage(2, "Reading hex connectivity\n");

    hexToVert = new GR_index_t[nHex][8];

    for (i = 0; i < nHex; i++) {
      for (j = 0; j < 8; j++) {
	nRead = fread(raw, 1, 4, fp);
	hexToVert[i][j] = convertToInt(raw);
	minVert = (hexToVert[i][j] < minVert) ? hexToVert[i][j] : minVert;
	maxVert = (hexToVert[i][j] > maxVert) ? hexToVert[i][j] : maxVert;
      }
    }
  }

  logMessage(2, "Min vert index = %u; max = %u\n", minVert, maxVert);

  for (GR_index_t iV = 0; iV < nVert; iV++) {
    (void) pVM->createVert(x[iV]);
  }
  logMessage(2, "Set coords of %u vertices\n", nVert);
  delete[] x;

  // Now, if minVert = 1 (Fortran numbering), decrement the indices in
  // all the connectivity tables.
  assert(minVert == 0 || minVert == 1);
  if (minVert == 1) {
    for (GR_index_t iT = 0; iT < nBdryTri; iT++) {
      triToVert[iT][0]--;
      triToVert[iT][1]--;
      triToVert[iT][2]--;
    }

    for (GR_index_t iQ = 0; iQ < nBdryQuad; iQ++) {
      quadToVert[iQ][0]--;
      quadToVert[iQ][1]--;
      quadToVert[iQ][2]--;
      quadToVert[iQ][3]--;
    }

    for (GR_index_t iT = 0; iT < nTet; iT++) {
      tetToVert[iT][0]--;
      tetToVert[iT][1]--;
      tetToVert[iT][2]--;
      tetToVert[iT][3]--;
    }

    for (GR_index_t iP = 0; iP < nPyr; iP++) {
      pyrToVert[iP][0]--;
      pyrToVert[iP][1]--;
      pyrToVert[iP][2]--;
      pyrToVert[iP][3]--;
      pyrToVert[iP][4]--;
    }

    for (GR_index_t iP = 0; iP < nPrism; iP++) {
      prismToVert[iP][0]--;
      prismToVert[iP][1]--;
      prismToVert[iP][2]--;
      prismToVert[iP][3]--;
      prismToVert[iP][4]--;
      prismToVert[iP][5]--;
    }

    for (GR_index_t iH = 0; iH < nHex; iH++) {
      hexToVert[iH][0]--;
      hexToVert[iH][1]--;
      hexToVert[iH][2]--;
      hexToVert[iH][3]--;
      hexToVert[iH][4]--;
      hexToVert[iH][5]--;
      hexToVert[iH][6]--;
      hexToVert[iH][7]--;
    }
  }
  // Done decrementing all the connectivity tables.

  // Now do the conversion
  pVM->convertFromCellVert(nTet, nPyr, nPrism, nHex, nBdryTri, nBdryQuad,
			   tetToVert, pyrToVert, prismToVert, hexToVert,
			   triToVert, quadToVert);
  if (tetToVert)
    delete[] tetToVert;
  if (pyrToVert)
    delete[] pyrToVert;
  if (prismToVert)
    delete[] prismToVert;
  if (hexToVert)
    delete[] hexToVert;
  if (triToVert)
    delete[] triToVert;
  if (quadToVert)
    delete[] quadToVert;
}
