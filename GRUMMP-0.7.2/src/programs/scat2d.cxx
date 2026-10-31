#include "GR_assert.h"
#ifndef WIN32
#include <unistd.h>
#else /* WIN32 defined */
int getopt(int argc, char *argv[], char *opstring);
#endif

#include <vector>

#include <stdio.h>
#include <string.h>
#include "GR_Mesh2D.h"
#include "GR_misc.h"
#include "GR_ADT.h"

static void
vUsage()
{
  logMessage(0, "Usage: scat2d -i basefilename [options] \n");
  logMessage(0, "   -m        Output 2D mesh [not done by default] \n\n");
  exit(1);
}

int
main(int iNArg, char *apcArgs[])
{
  char strBaseFileName[FILE_NAME_LEN], strDataFileName[FILE_NAME_LEN],
      strPointFileName[FILE_NAME_LEN], strOutputFileName[FILE_NAME_LEN];
  int iPChar;

  GRUMMPInit("scat2d");

  bool qError = true, qSaveMesh = false;
  while ((iPChar = getopt(iNArg, apcArgs, "i:m")) != EOF) {
    switch (iPChar)
      {
      case 'i':
	// Grab data file name in a way that prevents buffer overruns
	strncpy(strBaseFileName, optarg, FILE_NAME_LEN - 1);
	{
	  int iDum;
	  for (iDum = 0; iDum < FILE_NAME_LEN && strBaseFileName[iDum];
	      iDum++) {
	  }
	  if (iDum == FILE_NAME_LEN) // No null; file name too long
	    vFatalError("File name too long", "command line");
	}
	qError = false;
	break;
      case 'm':
	qSaveMesh = true;
	break;
      default:
	vUsage();
	break;
      }
  }
  if (qError)
    vUsage();

  // Open log file
  openMessageFile(strBaseFileName);

  // Open input file
  makeFileName(strDataFileName, "%s.data", strBaseFileName,
	       "main() [scat2d.C]");
  FILE *pFileIn = fopen(strDataFileName, "r");
  if (pFileIn == NULL)
    vFatalError("Scattered data file does not exist or could not be opened",
		"reading 2D scattered data");

  // Read the scattered data points and the data
  int iNVerts, iDataPerPoint;
  const int iBufSize = 2048;
  char *acBuffer = new char[iBufSize];

  vGetLineOrAbort(acBuffer, iBufSize, pFileIn);
  if (2 != sscanf(acBuffer, "%d %d", &iNVerts, &iDataPerPoint))
    vFatalError("Problem with mesh size info", "reading 2D scattered data");

  EntContainer<Vert> ECInput;
  ECInput.reserve(iNVerts);

  double **a2dData = new double*[iNVerts];
  int i;
  for (i = 0; i < iNVerts; i++) {
    a2dData[i] = new double[iDataPerPoint];
  }
  for (i = 0; i < iNVerts; i++) {
    double adCoord[3];
    char *pcData = acBuffer;
    char strErrorString[iBufSize];
    sprintf(strErrorString, "Insufficient data on line %d of file %s.", i,
	    strDataFileName);

    vGetLineOrAbort(acBuffer, iBufSize, pFileIn);
    vGetDoubleFromBuffer(&pcData, &(adCoord[0]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[1]), strErrorString, "main()");
    ECInput.getEntry(i)->setCoords(2, adCoord);

    for (int iData = 0; iData < iDataPerPoint; iData++)
      vGetDoubleFromBuffer(&pcData, &(a2dData[i][iData]), strErrorString,
			   "main()");
  }

  // Create a 2D mesh from the data points
  Mesh2D M2D(ECInput);
  assert(M2D.isSimplicial());

  // Create a tree of the cells from that mesh to make interpolation
  // efficient
  ADT* BBoxTree;
  createADTFromAllCells(&M2D, BBoxTree);

  // Read the target data points, interpolate on the fly, and spit out
  // the results.
  makeFileName(strPointFileName, "%s.points", strBaseFileName,
	       "main() [scat2d.C]");
  FILE *pFilePts = fopen(strPointFileName, "r");
  if (pFilePts == NULL)
    vFatalError("Target point file does not exist or could not be opened",
		"reading 2D target points for scattered data interp");

  int iNTargetVerts;
  vGetLineOrAbort(acBuffer, iBufSize, pFilePts);
  if (1 != sscanf(acBuffer, "%d", &iNTargetVerts))
    vFatalError("Problem with mesh size info",
		"reading 2D target points for scattered data interp");

  FILE *pF = NULL;
  makeFileName(strOutputFileName, "%s.out", strBaseFileName,
	       "main() [scat2d.C]");

  if (NULL == (pF = fopen(strOutputFileName, "w")))
    vFatalError("Couldn't open output file", "2D scattered data interp");

  for (i = 0; i < iNTargetVerts; i++) {
    double adCoord[3];
    char *pcData = acBuffer;
    char strErrorString[iBufSize];
    sprintf(strErrorString, "Insufficient data on line %d of file %s.", i,
	    strPointFileName);

    // Read the target data points
    vGetLineOrAbort(acBuffer, iBufSize, pFilePts);
    vGetDoubleFromBuffer(&pcData, &(adCoord[0]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[1]), strErrorString, "main()");

    double adBestBary[4];
    Cell *pC = findCellContainingPoint(M2D, BBoxTree, adCoord, adBestBary);
    // Spit out the points...
    fprintf(pF, "%g %g", adCoord[0], adCoord[1]);

    // Interpolate data to the target points
    for (int iData = 0; iData < iDataPerPoint; iData++) {
      double dTargetData = 0;
      for (int ii = 0; ii < 3; ii++) {
	int iVert = M2D.iVertIndex(pC->getVert(ii));
	dTargetData += a2dData[iVert][iData] * adBestBary[ii];
      }
      fprintf(pF, " %g", dTargetData);
    }
    fprintf(pF, "\n");
  }

  // Now delete a bunch of arrays.
  for (i = 0; i < iNVerts; i++) {
    delete[] a2dData[i];
  }
  delete[] a2dData;
  delete[] acBuffer;

  fclose(pF);
  fclose(pFilePts);
  fclose(pFileIn);

  if (qSaveMesh) {
    logMessage(1, "Writing 2D mesh file...\n");
    writeVTK(M2D, strBaseFileName);
    logMessage(1, "Done\n");
  }

  return (0);
}

