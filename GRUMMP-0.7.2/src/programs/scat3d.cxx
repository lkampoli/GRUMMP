#include "GR_assert.h"
#ifndef WIN32
#include <unistd.h>
#else /* WIN32 defined */
int getopt(int argc, char *argv[], char *opstring);
#endif

#include <stdio.h>
#include <string.h>
#include "GR_VolMesh.h"
#include "GR_misc.h"
#include "GR_ADT.h"

static void
vUsage()
{
  logMessage(0, "Usage: scat3d -i basefilename [options] \n");
  logMessage(0, "   -m        Output volume mesh [not done by default] \n\n");
  exit(1);
}

int
main(int iNArg, char *apcArgs[])
{
  char strBaseFileName[FILE_NAME_LEN], strDataFileName[FILE_NAME_LEN],
      strPointFileName[FILE_NAME_LEN], strOutputFileName[FILE_NAME_LEN];
  int iPChar;

  GRUMMPInit("scat3d");

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
	       "main() [scat3d.C]");
  FILE *pFileIn = fopen(strDataFileName, "r");
  if (pFileIn == NULL)
    vFatalError("Scattered data file does not exist or could not be opened",
		"reading 3D scattered data");

  // Read the scattered data points and the data
  GR_index_t iNVerts;
  int iDataPerPoint;
  const int iBufSize = 2048;
  char *acBuffer = new char[iBufSize];

  vGetLineOrAbort(acBuffer, iBufSize, pFileIn);
  if (2 != sscanf(acBuffer, "%u %d", &iNVerts, &iDataPerPoint))
    vFatalError("Problem with mesh size info", "reading 3D scattered data");

  EntContainer<Vert> ECInput;
  ECInput.reserve(iNVerts);

  double **a2dData = new double*[iNVerts];
  GR_index_t i;
  for (i = 0; i < iNVerts; i++) {
    a2dData[i] = new double[iDataPerPoint];
  }
  for (i = 0; i < iNVerts; i++) {
    double adCoord[3];
    char *pcData = acBuffer;
    char strErrorString[iBufSize];
    sprintf(strErrorString, "Insufficient data on line %u of file %s.", i,
	    strDataFileName);

    vGetLineOrAbort(acBuffer, iBufSize, pFileIn);
    vGetDoubleFromBuffer(&pcData, &(adCoord[0]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[1]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[2]), strErrorString, "main()");
    ECInput.getEntry(i)->setCoords(3, adCoord);

    for (int iData = 0; iData < iDataPerPoint; iData++)
      vGetDoubleFromBuffer(&pcData, &(a2dData[i][iData]), strErrorString,
			   "main()");
  }

  // Create a volume mesh from the data points
  VolMesh VM(ECInput);
  assert(VM.isSimplicial());

  // Create a tree of the cells from that mesh to make interpolation
  // efficient
  ADT* BBoxTree;
  createADTFromAllCells(&VM, BBoxTree);

  // Read the target data points, interpolate on the fly, and spit out
  // the results.
  makeFileName(strPointFileName, "%s.points", strBaseFileName,
	       "main() [scat3d.C]");
  FILE *pFilePts = fopen(strPointFileName, "r");
  if (pFilePts == NULL)
    vFatalError("Target point file does not exist or could not be opened",
		"reading 3D target points for scattered data interp");

  GR_index_t iNTargetVerts;
  vGetLineOrAbort(acBuffer, iBufSize, pFilePts);
  if (1 != sscanf(acBuffer, "%u", &iNTargetVerts))
    vFatalError("Problem with mesh size info",
		"reading 3D target points for scattered data interp");

  FILE *pF = NULL;
  makeFileName(strOutputFileName, "%s.out", strBaseFileName,
	       "main() [scat3d.C]");

  if (NULL == (pF = fopen(strOutputFileName, "w")))
    vFatalError("Couldn't open output file", "3D scattered data interp");

  std::vector<GR_index_t> veciQueryResult;
  for (i = 0; i < iNTargetVerts; i++) {
    double adCoord[3];
    char *pcData = acBuffer;
    char strErrorString[iBufSize];
    sprintf(strErrorString, "Insufficient data on line %u of file %s.", i,
	    strPointFileName);

    // Read the target data points
    vGetLineOrAbort(acBuffer, iBufSize, pFilePts);
    vGetDoubleFromBuffer(&pcData, &(adCoord[0]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[1]), strErrorString, "main()");
    vGetDoubleFromBuffer(&pcData, &(adCoord[2]), strErrorString, "main()");

    // Find cell (nearly?) containing the target point
    double adBestBary[4];
    Cell *pC = findCellContainingPoint(VM, BBoxTree, adCoord, adBestBary);

    // Spit out the points...
    fprintf(pF, "%g %g %g", adCoord[0], adCoord[1], adCoord[2]);

    // Interpolate data to the target points
    for (int iData = 0; iData < iDataPerPoint; iData++) {
      double dTargetData = 0;
      for (int ii = 0; ii < 4; ii++) {
	int iVert = VM.getVertIndex(pC->getVert(ii));
	dTargetData += a2dData[iVert][iData] * adBestBary[ii];
      }
      fprintf(pF, " %g", dTargetData);
    }
    fprintf(pF, "\n");
  }

  for (i = 0; i < iNVerts; i++) {
    delete[] a2dData[i];
  }
  delete[] a2dData;
  delete[] acBuffer;

  fclose(pF);
  fclose(pFilePts);
  fclose(pFileIn);

  if (qSaveMesh) {
    logMessage(1, "Writing volume mesh file...\n");
    //writeFile_VolMesh(VM, strBaseFileName);
    logMessage(1, "Done\n");
  }

  return (0);
}

