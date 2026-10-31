/*
 * ShockTopo.cxx
 *
 *  Created on: Nov 3, 2017
 *      Author: cfog
 */

#include "GR_ShockTopo.h"
#include "GR_misc.h"

namespace GRUMMP
{
  // This function is a little silly, but it does centralize this process,
  // which is likely to be helpful if I ever have to do something more
  // sophisticated here.
  static bool
  isValidID(const int ID)
  {
    return (ID >= 0);
  }
  ShockEndPoint::ShockEndPoint(const int ID, const int nCurves,
			       const int curves[]) :
      m_ID(ID), m_numCurves(nCurves)
  {
    assert(isValidID(ID));
    assert(nCurves >= 2 && nCurves <= maxIncCurves);
    for (int ii = 0; ii < nCurves; ii++) {
      assert(isValidID(curves[ii]));
      m_curveID[ii] = curves[ii];
    }
  }

  ShockCurve::ShockCurve(const int ID, const InteractionType type,
			 const int nSurfs, const int surfs[],
			 const int endPoints[]) :
      m_ID(ID), m_numSurfs(nSurfs), m_type(type)
  {
    assert(isValidID(ID));
    assert(nSurfs >= 1 && nSurfs <= maxIncSurfs);
    assert(nSurfs == getNumReqSurfs(type));
    for (int ii = 0; ii < nSurfs; ii++) {
      assert(isValidID(surfs[ii]));
      m_surfID[ii] = surfs[ii];
    }
    assert(isValidID(endPoints[0]));
    assert(isValidID(endPoints[1]));
    m_endpointID[0] = endPoints[0];
    m_endpointID[1] = endPoints[1];
  }

  ShockSurface::ShockSurface(const int ID, const int nCurves,
			     const int curves[]) :
      m_ID(ID), m_numCurves(nCurves)
  {
    assert(isValidID(ID));
    assert(nCurves >= 1 && nCurves <= maxBoundingCurves);
    for (int ii = 0; ii < nCurves; ii++) {
      assert(isValidID(curves[ii]));
      m_curveID[ii] = curves[ii];
    }
  }

  static void
  getNoncommentLine(char **lineBuffer, size_t* lineLength, FILE* stream,
		    size_t& lineNumber)
  {
    bool lineOkay = false;
    do {
      lineNumber++;
      ssize_t result = getline(lineBuffer, lineLength, stream);
      if (result == EOF) {
	vFatalError("Premature end of file", "reading shock topology file");
      }
      lineOkay = ((*lineBuffer)[0] != '#');
    }
    while (!lineOkay);
  }

  void
  ShockTopo::checkReadResult(const int result, const int expected,
			     const size_t lineNumber) const
  {
    if (result != expected) {
      logMessage(MSG_WARNING, "Error in line %lu\n", lineNumber);
      vFatalError("File format error", "reading shock topology");
    }
  }

  ShockTopo::ShockTopo(const char topoFileName[])
  {
    FILE* topoFile = fopen(topoFileName, "r");

    char* lineBuffer;
      size_t lineLength = 1024, lineNumber = 0;
    lineBuffer = reinterpret_cast<char*>(calloc(lineLength, sizeof(char)));

    int nPoints(0), nCurves(0), nSurfs(0);
    getNoncommentLine(&lineBuffer, &lineLength, topoFile, lineNumber);

    int result = sscanf(lineBuffer, "%d%d%d", &nPoints, &nCurves, &nSurfs);
    checkReadResult(result, 3, lineNumber);

    m_endpoints.reserve(nPoints);
    m_curves.reserve(nCurves);
    m_surfs.reserve(nSurfs);

    // Read and store the end point data.
    for (int ii = 0; ii < nPoints; ii++) {
      int nPtRedundant, nIncCurves, curves[ShockEndPoint::getMaxIncCurves()];
      getNoncommentLine(&lineBuffer, &lineLength, topoFile, lineNumber);

      char *readPoint = lineBuffer;
      int charsRead;
      result = sscanf(readPoint, "EP %d %d%n", &nPtRedundant, &nIncCurves,
		      &charsRead);
      checkReadResult(result, 2, lineNumber);
      readPoint += charsRead;
      assert(nPtRedundant == ii);
      assert(nIncCurves >= 2 && nIncCurves <= ShockEndPoint::getMaxIncCurves());
      // Now read the curve IDs
      for (int jj = 0; jj < nIncCurves; jj++) {
	result = sscanf(readPoint, "%d%n", curves + jj, &charsRead);
	checkReadResult(result, 1, lineNumber);
	readPoint += charsRead;
      }
      // Create a new ShockEndPoint object and add it to the shock topo data.
      ShockEndPoint SEP(ii, nIncCurves, curves);
      addShockEndPoint(SEP);
    }

    // Read and store the curve data.
    for (int ii = 0; ii < nCurves; ii++) {
      int nCurveRedundant, surfs[ShockCurve::getMaxIncSurfaces()];
      int endpoints[2], typeInt;
      getNoncommentLine(&lineBuffer, &lineLength, topoFile, lineNumber);

      char *readPoint = lineBuffer;
      int charsRead;
      result = sscanf(readPoint, "CU %d %d %d %d%n", &nCurveRedundant,
		      &typeInt, endpoints, endpoints + 1, &charsRead);
      checkReadResult(result, 4, lineNumber);
      enum ShockCurve::InteractionType type = ShockCurve::InteractionType(
	  typeInt);
      readPoint += charsRead;
      assert(nCurveRedundant == ii);
      int nIncSurfs = ShockCurve::getNumReqSurfs(type);
      // Now read the surf IDs.
      for (int jj = 0; jj < nIncSurfs; jj++) {
	result = sscanf(readPoint, "%d%n", surfs + jj, &charsRead);
	checkReadResult(result, 1, lineNumber);
	readPoint += charsRead;
      }
      // Create a new ShockCurve object and add it to the shock topo data.
      ShockCurve SC(ii, type, nIncSurfs, surfs, endpoints);
      addShockCurve(SC);
    }

    // Read and store the surface data.  Not reading a surface type at this point.
    for (int ii = 0; ii < nSurfs; ii++) {
      int nSurfRedundant, nIncCurves,
	  curves[ShockSurface::getMaxBoundingCurves()];
      getNoncommentLine(&lineBuffer, &lineLength, topoFile, lineNumber);

      char *readPoint = lineBuffer;
      int charsRead;
      result = sscanf(readPoint, "SU %d %d%n", &nSurfRedundant, &nIncCurves,
		      &charsRead);
      checkReadResult(result, 2, lineNumber);
      readPoint += charsRead;
      assert(nSurfRedundant == ii);
      assert(nIncCurves >= 2 && nIncCurves <= ShockEndPoint::getMaxIncCurves());
      // Now read the curve IDs
      for (int jj = 0; jj < nIncCurves; jj++) {
	result = sscanf(readPoint, "%d%n", curves + jj, &charsRead);
	checkReadResult(result, 1, lineNumber);
	readPoint += charsRead;
      }
      // Create a new ShockSurface object and add it to the shock topo data.
      ShockSurface SS(ii, nIncCurves, curves);
      addShockSurface(SS);
    }
    free(lineBuffer);
    assert(isValid());
  }

  ShockEndPoint
  ShockTopo::getShockEndPoint(const int ID) const
  {
    return m_endpoints[ID];
  }

  ShockCurve
  ShockTopo::getShockCurve(const int ID) const
  {
    return m_curves[ID];
  }

  ShockSurface
  ShockTopo::getShockSurface(const int ID) const
  {
    return m_surfs[ID];
  }

  void
  ShockTopo::addShockEndPoint(const ShockEndPoint SEP)
  {
    m_endpoints.push_back(SEP);
  }

  void
  ShockTopo::addShockCurve(const ShockCurve SC)
  {
    m_curves.push_back(SC);
  }

  void
  ShockTopo::addShockSurface(const ShockSurface SS)
  {
    m_surfs.push_back(SS);
  }

  int
  ShockCurve::getNumReqSurfs(const InteractionType type)
  {
    switch (type)
      {
      case ShockPerimeter:
	return 1;
      case ShockBdry:
	return 1;
      case Reflection:
	return 2;
      case TripleLine:
	return 4;
      case QuadrupleLine:
	return 5;
      case Invalid:
      default:  // There are no other valid values, but just in case.
	assert(0);
	return 0;
      }
  }

  bool
  ShockTopo::isValid() const
  {
    // Every end point should be mentioned by each of its curve
    for (int ii = getNumEndPoints() - 1; ii >= 0; ii--) {
      ShockEndPoint SEP = getShockEndPoint(ii);
      int myID = SEP.getID();
      for (int jj = SEP.getNumCurves() - 1; jj >= 0; jj--) {
	int curveID = SEP.getCurveID(jj);
	if ((getShockCurve(curveID).getEndPointID(0) != myID)
	    && (getShockCurve(curveID).getEndPointID(1) != myID)) {
	  return false;
	}
      }
    }

    // It's possible that there could be extra curves that refer to endpoints
    // that don't exist.  I'm not checking for that.

    // Every curve should be mentioned by each of its surfaces
    for (int ii = getNumCurves() - 1; ii >= 0; ii--) {
      ShockCurve SC = getShockCurve(ii);
      int myID = SC.getID();
      for (int jj = SC.getNumSurfaces() - 1; jj >= 0; jj--) {
	int surfID = SC.getSurfaceID(jj);
	ShockSurface SS = getShockSurface(surfID);
	bool foundIt = false;
	for (int kk = SS.getNumCurves() - 1; kk >= 0; kk--) {
	  if (SS.getCurveID(kk) == myID) {
	    foundIt = true;
	    break;
	  }
	} // Checking every curve for this surface
	if (!foundIt)
	  return false;
      } // Checking every surface incident on this curve
    } // Checking every curve in the topology
    return true;
  }

} // Closing out namespace GRUMMP
