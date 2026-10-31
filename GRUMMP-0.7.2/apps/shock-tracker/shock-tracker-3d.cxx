/*
 * shock-tracker-3d.cxx
 *
 *  Created on: Apr 13, 2017
 *      Author: cfog
 */

#include <string.h>
#include <unistd.h>

#include "GR_InterfaceTracker.h"

static void
printUsage()
{
  logMessage(0,
	     "Usage: shock-tracker-3d -i mesh_filename -s surface_filename\n");
  logMessage(0,
	     "  mesh_filename must be the base name for a tetgen file set\n");
  logMessage(0, "  surface_filename must specify an STL file\n");
  exit(1);
}

int
main(int nArgs, char *args[])
{
  GRUMMP_Apps::InterfaceTracker IT;
  GRUMMPInit("shock-tracker-3d");
  openMessageFile("shock-track");

  char strInputMeshFileName[FILE_NAME_LEN], strSurfFileName[FILE_NAME_LEN];
  bool gotInputFileName = false, gotSurfFileName = false, experimental = false;
  int charRead = EOF;

  while ((charRead = getopt(nArgs, args, "i:s:x:")) != EOF) {
    switch (charRead)
      {
      case 'i':
	strncpy(strInputMeshFileName, optarg, FILE_NAME_LEN - 1);
	gotInputFileName = true;
	break;
      case 's':
	strncpy(strSurfFileName, optarg, FILE_NAME_LEN - 1);
	gotSurfFileName = true;
	break;
      case 'x':
	strncpy(strSurfFileName, optarg, FILE_NAME_LEN - 1);
	experimental = true;
	break;
      default:
	printUsage();
	break;
      }
  }
  if (!(gotInputFileName)) {
	printUsage();
      }

  IT.readVolumeMesh(strInputMeshFileName);

  if (gotSurfFileName) {
    IT.addInterfaceToMesh(strSurfFileName);
  }
  else if (experimental) {
    GRUMMP_Apps::readAndManipulateSurface(strSurfFileName);
  }
  else {
    IT.updateMesh();
  }

  // Output the volume mesh to VTK and Tetgen formats
  IT.writeVolumeMesh();

  exit(0);
}

