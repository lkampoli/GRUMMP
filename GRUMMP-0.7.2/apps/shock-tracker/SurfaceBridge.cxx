/*
 * SurfaceBridge.cxx
 *
 *  Created on: Nov 5, 2017
 *      Author: cfog
 */

#include <math.h>

#include "GR_InterfaceTracker.h"
#include "GR_SurfMesh.h"

namespace GRUMMP_Apps
{
  void
  readAndManipulateSurface(const char fileName[])
  {
    SurfMesh SM;
    SM.readFromFile(fileName);
    std::map<Vert*, Vert*> quillPairs;

    for (GR_index_t ii = 0; ii < SM.getNumBdryFaces(); ii++) {
      BFace* pBF = SM.getBFace(ii);
      assert(pBF->getNumVerts() == 2);
      for (int jj = 0; jj < 2; jj++) {
	Vert *perimVert = pBF->getVert(jj);
	if (quillPairs.count(perimVert) != 0)
	  continue; // Already did this vert.
	double stepSize = -1;
	double stepDir[] =
	  { 0, 0, 0 };
	// Grab the neighborhood of this vert so we can use it for geometry info.
	{
	  std::set<Cell*> spCInc;
	  std::set<Vert*> spVNeigh;
	  std::set<BFace*> spBFInc;
	  bool isBdry = true;
	  findNeighborhoodInfo(perimVert, spCInc, spVNeigh, &spBFInc, &isBdry);
	  assert(isBdry);
	  assert(spBFInc.size() == 2);
	  // Find a surface normal by summing the (area-scaled) normals of all the incident
	  // triangles.
	  double normal[] =
	    { 0, 0, 0 };
	  for (Cell* thisCell : spCInc) {
	    double thisNormal[3];
	    TriCell* tri = dynamic_cast<TriCell*>(thisCell);
	    assert(tri);
	    tri->calcVecSize(thisNormal);
	    double whichWay = copysign(1, dDOT3D(thisNormal, normal));
	    for (int comp = 0; comp < 3; comp++) {
	      normal[comp] += whichWay * thisNormal[comp];
	    }
	  }
	  // Estimate lengthscale based on total area of triangles incident.
	  stepSize = 2 * sqrt(dMAG3D(normal));
	  NORMALIZE3D(normal);

	  // Now find a direction out (tangentially) from this surface.
	  for (BFace* thisBFace : spBFInc) {
	    const double *coords0 = thisBFace->getVert(0)->getCoords();
	    const double *coords1 = thisBFace->getVert(1)->getCoords();
	    Face* faceAtBdry = thisBFace->getFace();
	    Cell* triAtBdry = faceAtBdry->getOppositeCell(thisBFace);
	    assert(triAtBdry != thisBFace && triAtBdry->getNumVerts() == 3);
	    Vert* oppVert = dynamic_cast<TriCell*>(triAtBdry)->getOppositeVert(
		faceAtBdry);
	    const double *coords2 = oppVert->getCoords();
	    double edgeVec01[] = adDIFF3D(coords1, coords0);
	    double edgeVec02[] = adDIFF3D(coords2, coords0);
	    double cross[3];
	    vCROSS3D(edgeVec01, normal, cross);
	    double whichWay = copysign(1, -dDOT3D(edgeVec02, cross));
	    for (int comp = 0; comp < 3; comp++) {
	      stepDir[comp] += whichWay * cross[comp];
	    }
	  }
	  NORMALIZE3D(stepDir);
	}

	// Now create a vert out along that line, running a prescribed number
	// of local length scales.
	double newLoc[] =
	  { perimVert->x() + stepSize * stepDir[0], perimVert->y()
	      + stepSize * stepDir[1], perimVert->z() + stepSize * stepDir[2] };
	Vert *newVert = SM.createVert(newLoc);
	quillPairs.insert(make_pair(perimVert, newVert));
      } // Did both verts for this edge.

	// Now build a quad for that bdry edge
      Vert *pVA = pBF->getVert(0);
      Vert *pVD = pBF->getVert(1);
      Vert *pVB = quillPairs.find(pVA)->second;
      Vert *pVC = quillPairs.find(pVD)->second;

      Face* edgeOnBdry = pBF->getFace();
      Cell* cellNextToBdry = edgeOnBdry->getOppositeCell(pBF);
      int reg = cellNextToBdry->getRegion();
      SM.deleteBFace(pBF);
      edgeOnBdry->setFaceLoc(Face::eInterior);
      SM.createTriCell(pVA, pVB, pVC, reg);
      SM.createTriCell(pVA, pVC, pVD, reg);
      edgeOnBdry = findCommonFace(pVB, pVC);
      SM.createBFace(edgeOnBdry);

    } // Done with loop over bdry faces
    assert(SM.isValid());
      // Write a VTK file to show off the result.
    writeVTK(SM, "with-quills.vtk");
  }
}
