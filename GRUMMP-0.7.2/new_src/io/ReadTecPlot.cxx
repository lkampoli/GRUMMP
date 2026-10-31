/*
 * ReadTecPlot.cxx
 *
 *  Created on: Apr 4, 2017
 *      Author: cfog
 */

#include <stdio.h>
#include <sstream>

#include "GR_VolMesh.h"

void
readTecPlotASCII(const char* const fileName, VolMesh * const pVM)
{
  GR_index_t nTets, nNodes;

  FILE* pFInFile = fopen(fileName, "r");
  if (!pFInFile) {
    vFatalError("File not found", "readTecPlotASCII");
  }
  logMessage(1, "Reading TecPlot file %s\n", fileName);

  // Read title line and discard
  char *line = new char[32 * 1024 + 1]; // TecPlot line limit is 32,000 chars
  size_t lineLen = 0;
  getline(&line, &lineLen, pFInFile);

  // Read variables line and parse for the amount of data on each line
  getline(&line, &lineLen, pFInFile);
  const char* where = strstr(line, "="); // This way the space or lack thereof before = doesn't matter.

  std::istringstream iss(where + 1);
  int numVars = -1;

  while (!iss.eof()) {
    string name;
    iss >> name;
    numVars += 1;
  }
  numVars -= 3; // Coordinates aren't counted here.
  assert(numVars >= 0);

  // Read ZONE line and parse
  getline(&line, &lineLen, pFInFile);
  // Check to be sure it's a tetrahedral mesh.
  where = strcasestr(line, "tetrahedron");
  if (where == nullptr) {
    vFatalError("Not a tetrahedral mesh", "readTecPlotASCII");
  }
  where = strcasestr(line, "N=");
  if (where == nullptr) {
    vFatalError("Number of nodes not given", "readTecPlotASCII");
  }
  else {
    sscanf(where + 2, "%ud", &nNodes);
  }

  where = strcasestr(line, "E=");
  if (where == nullptr) {
    vFatalError("Number of tetrahedra not given", "readTecPlotASCII");
  }
  else {
    sscanf(where + 2, "%ud", &nTets);
  }
  logMessage(1, "  File has %d tets, %d verts\n", nTets, nNodes);

  // Read and store coordinates
  assert(pVM->getNumVerts() == 0);
  for (GR_index_t iV = 0; iV < nNodes; iV++) {
    double x, y, z;
    int result = fscanf(pFInFile, "%lg %lg %lg", &x, &y, &z);
    assert(result == 3);
    // Discard solution data
    for (int ii = 0; ii < numVars; ii++) {
      fscanf(pFInFile, "%*g");
    }
    pVM->createVert(x, y, z);
    if (iV * 10 % nNodes == 0) {
      logMessage(1, "Processed %d verts\n", iV);
    }
  }

  // Read connectivity and store it.
  GR_index_t (*connect)[4] = new GR_index_t[nTets][4];
  GR_index_t minIndex = INT_MAX, maxIndex = 0;
  for (GR_index_t iT = 0; iT < nTets; iT++) {
    int result = fscanf(pFInFile, "%u %u %u %u", &connect[iT][0],
			&connect[iT][1], &connect[iT][2], &connect[iT][3]);
    assert(result == 4);
    GR_index_t thisMin = std::min(std::min(connect[iT][0], connect[iT][1]),
				  std::min(connect[iT][2], connect[iT][3]));
    GR_index_t thisMax = std::max(std::max(connect[iT][0], connect[iT][1]),
				  std::max(connect[iT][2], connect[iT][3]));
    if (thisMin < minIndex)
      minIndex = thisMin;
    if (thisMax > maxIndex)
      maxIndex = thisMax;
  }

  if (maxIndex - minIndex + 1 > nNodes) {
    vFatalError("Some node indices are out of range", "readTecPlotASCII");
  }

  // Now create tets from verts.
  for (GR_index_t iT = 0; iT < nTets; iT++) {
    bool qAlreadyExisted = false;
    // Default region is fine here.
    pVM->createTetCell(qAlreadyExisted, pVM->getVert(connect[iT][0] - minIndex),
		       pVM->getVert(connect[iT][1] - minIndex),
		       pVM->getVert(connect[iT][2] - minIndex),
		       pVM->getVert(connect[iT][3] - minIndex));
    if (qAlreadyExisted) {
      vFatalError("Duplicate tet!", "readTecPlotASCII");
    }
  }

  // Create BFace's where needed for completeness.
  for (int i = pVM->getNumFaces() - 1; i >= 0; i--) {
    Face *pF = pVM->getFace(i);
    if (!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid()) {
      pVM->createBFace(pF);
    }
  }

  fclose(pFInFile);

  assert(pVM->isValid());
  delete[] connect;
  delete[] line;
}
