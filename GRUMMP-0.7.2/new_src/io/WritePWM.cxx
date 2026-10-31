#include <stdio.h>
#include "GR_misc.h"
#include "GR_VolMesh.h"
#include "GR_Geometry.h"

//-----------------------------------------------------------------------------
void
writeFile_PWM(VolMesh& OutMesh, const char strBaseFileName[],
	      const char strExtraSuffix[])
//-----------------------------------------------------------------------------
{
  FILE *pFOutFile = NULL;
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  if ((strlen(strBaseFileName) >= 4)
      && strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".pwm")) {
    // File name doesn't end in .pwm
    sprintf(strFileName, "%s.pwm.%s", strBaseFileName, strExtraSuffix);
  }
  else {
    // File name ends in .pwm
    sprintf(strFileName, "%s.%s", strBaseFileName, strExtraSuffix);
  }

  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing",
		"3d Pointwise mesh output");

  //-------------------------------------------------------
  // Write the PWM header details
  //-------------------------------------------------------

  // File magic
  fprintf(pFOutFile, "1994\n");

  // File authorship info
  fprintf(pFOutFile, "PWM written by GRUMMP\n");

  // Skip a couple of integers (size of int and real, I believe)
  fprintf(pFOutFile, "4\n8\n");

  // Now read the mesh size data.
  int nVerts = OutMesh.getNumVerts();
  int nTets = OutMesh.getNumTetCells();
  fprintf(pFOutFile, "%d\n 0 0 0 %d 0 0 0 0\n", nVerts, nTets);

  //-------------------------------------
  // write 3d vertex data
  //-------------------------------------
  for (int i = 0; i < nVerts; ++i) {
    Vert *pV = OutMesh.getVert(i);
    fprintf(pFOutFile, "%20.16f %20.16f %20.16f\n", pV->x(), pV->y(), pV->z());
  }

  for (int i = 0; i < nTets; i++) {
    Cell *pC = OutMesh.getCell(i);
    fprintf(pFOutFile, "%d %d %d %d\n", OutMesh.getVertIndex(pC->getVert(0)),
	    OutMesh.getVertIndex(pC->getVert(1)),
	    OutMesh.getVertIndex(pC->getVert(2)),
	    OutMesh.getVertIndex(pC->getVert(3)));
    assert(
	checkOrient3D(pC->getVert(0), pC->getVert(1), pC->getVert(2),
		      pC->getVert(3)) == 1);
  }

  fclose(pFOutFile);
}
