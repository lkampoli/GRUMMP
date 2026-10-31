/*
 * Length3D.cxx
 *
 *  Created on: Jun 27, 2017
 *      Author: cfog
 */

#include "GR_ADT.h"
#include "GR_InsertionManager.h"
#include "GR_Length.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Vec.h"
#include "GR_VolMesh.h"

namespace GRUMMP
{
  Length3D::Length3D(const Mesh* const & pM, const double& refine,
		     const double& grading, const double& minLength) :
      Length(refine, grading, minLength)
  {
    logMessage(MSG_MANAGER, "Constructing 3D length scale\n");
    const VolMesh* const pVM = dynamic_cast<const VolMesh* const >(pM);
    assert(pVM);
    m_lengthMesh = std::unique_ptr < Mesh > (new VolMesh());
    VolMesh* const pVM_int = dynamic_cast<VolMesh*>(m_lengthMesh.get());
    assert(pVM_int);
    for (GR_index_t iV = 0; iV < pM->getNumVerts(); iV++) {
      Vert *pV = pM->getVert(iV);
      pVM_int->createVert(pV->getCoords());
    }
    GR_index_t iC = 0;
    bool existsAlready = false;
    for (; iC < pVM->getNumTetCells(); iC++) {
      Cell* pC = pVM->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pVM->getVertIndex(pC->getVert(0));
      int iV1 = pVM->getVertIndex(pC->getVert(1));
      int iV2 = pVM->getVertIndex(pC->getVert(2));
      int iV3 = pVM->getVertIndex(pC->getVert(3));
      pVM_int->createTetCell(existsAlready, pVM_int->getVert(iV0),
			     pVM_int->getVert(iV1), pVM_int->getVert(iV2),
			     pVM_int->getVert(iV3));
      assert(!existsAlready);
    }
    for (; iC < pVM->getNumTetCells() + pVM->getNumPyrCells(); iC++) {
      Cell* pC = pVM->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pVM->getVertIndex(pC->getVert(0));
      int iV1 = pVM->getVertIndex(pC->getVert(1));
      int iV2 = pVM->getVertIndex(pC->getVert(2));
      int iV3 = pVM->getVertIndex(pC->getVert(3));
      int iV4 = pVM->getVertIndex(pC->getVert(4));
      pVM_int->createPyrCell(existsAlready, pVM_int->getVert(iV0),
			     pVM_int->getVert(iV1), pVM_int->getVert(iV2),
			     pVM_int->getVert(iV3), pVM_int->getVert(iV4));
      assert(!existsAlready);
    }
    for (; iC < pVM->getNumCells() - pVM->getNumPrismCells(); iC++) {
      Cell* pC = pVM->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pVM->getVertIndex(pC->getVert(0));
      int iV1 = pVM->getVertIndex(pC->getVert(1));
      int iV2 = pVM->getVertIndex(pC->getVert(2));
      int iV3 = pVM->getVertIndex(pC->getVert(3));
      int iV4 = pVM->getVertIndex(pC->getVert(4));
      int iV5 = pVM->getVertIndex(pC->getVert(5));
      pVM_int->createPrismCell(existsAlready, pVM_int->getVert(iV0),
			       pVM_int->getVert(iV1), pVM_int->getVert(iV2),
			       pVM_int->getVert(iV3), pVM_int->getVert(iV4),
			       pVM_int->getVert(iV5));
      assert(!existsAlready);
    }
    for (; iC < pVM->getNumCells(); iC++) {
      Cell* pC = pVM->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pVM->getVertIndex(pC->getVert(0));
      int iV1 = pVM->getVertIndex(pC->getVert(1));
      int iV2 = pVM->getVertIndex(pC->getVert(2));
      int iV3 = pVM->getVertIndex(pC->getVert(3));
      int iV4 = pVM->getVertIndex(pC->getVert(4));
      int iV5 = pVM->getVertIndex(pC->getVert(5));
      int iV6 = pVM->getVertIndex(pC->getVert(6));
      int iV7 = pVM->getVertIndex(pC->getVert(7));
      pVM_int->createHexCell(existsAlready, pVM_int->getVert(iV0),
			     pVM_int->getVert(iV1), pVM_int->getVert(iV2),
			     pVM_int->getVert(iV3), pVM_int->getVert(iV4),
			     pVM_int->getVert(iV5), pVM_int->getVert(iV6),
			     pVM_int->getVert(iV7));
      assert(!existsAlready);
    }
    pVM_int->purgeAllEntities();
    // Create a tree of the cells from that mesh to make interpolation
    // efficient
    createADTFromAllCells(pVM_int, m_BBoxTree);
    initLength();
    const Mesh* const tmp = getMesh();
    // Copy this to the originating mesh.
    for (GR_index_t ii = 0; ii < pM->getNumVerts(); ii++) {
      pM->getVert(ii)->setLengthScale(tmp->getVert(ii)->getLengthScale());
    }
    logMessage(MSG_MANAGER, "Done constructing 3D length scale\n");
  }

  double
  Length3D::distToNearestBdryEntity(const Vert* const vert,
				    const std::set<Cell*>& cells) const
  {
    double min_dist_squared = LARGE_DBL;
    for (auto pC : cells) {
      for (int ii = 0; ii < pC->getNumFaces(); ii++) {
	Face *pF = pC->getFace(ii);
	if (pF->hasVert(vert))
	  continue;
	if (pF->getFaceLoc() == Face::eBdryFace) {
	  Cell *bdryCell =
	      (pF->getLeftCell() == pC) ?
		  pF->getRightCell() : pF->getLeftCell();
	  TriBFaceBase* bface = dynamic_cast<TriBFaceBase*>(bdryCell);
	  assert(bface);
	  double nearest[] =
	    { 0, 0, 0 };
	  bface->findClosestPointOnBFace(vert->getCoords(), nearest);
	  double dist_squared = dDIST_SQ_3D(vert->getCoords(), nearest);

	  min_dist_squared = std::min(dist_squared, min_dist_squared);
	} // Yes, it's really a bdry face
      } // Done looping over edges of this cell
    } // Done looping over cells
    return sqrt(min_dist_squared);
  }

  double
  Length3D::queryLengthScale(const double coords[]) const
  {
    const Mesh* const mesh = m_lengthMesh.get();
    double length_scale = LARGE_DBL;
    std::vector<GR_index_t> veciQueryResult;

    // Find cell (nearly?) containing the target point
    double dEps = 1.e-8;
    do {
      double adBB[] =
	{ coords[0] - dEps, coords[0] + dEps, coords[1] - dEps, coords[1]
	    + dEps, coords[2] - dEps, coords[2] + dEps };

      veciQueryResult = m_BBoxTree->veciRangeQuery(adBB);
      dEps *= 10;
    }
    while (veciQueryResult.empty());

    double dBestValue = 1.e50;
    int iBestCell = -1;
    double adBestBary[4];
    Cell *pC;
    int iQ;
    GR_index_t iCell;
    for (iQ = veciQueryResult.size() - 1; iQ >= 0 && dBestValue > 0.5; iQ--) {
      // Get barycentric coords.  If best_value <(=) 0.5 on completion
      // then the vert lies within (on an edge of) the cell in list
      iCell = veciQueryResult[iQ];
      pC = mesh->getCell(iCell);

      double adBary[4];
      dynamic_cast<TetCell*>(pC)->calcBarycentricCoords(coords, adBary);

      double dThisValue = max(
	  max(max(fabs(adBary[0] - 0.5), fabs(adBary[1] - 0.5)),
	      fabs(adBary[2] - 0.5)),
	  fabs(adBary[3] - 0.5));

      if (dThisValue < dBestValue) {
	iBestCell = iCell;
	dBestValue = dThisValue;
	adBestBary[0] = adBary[0];
	adBestBary[1] = adBary[1];
	adBestBary[2] = adBary[2];
	adBestBary[3] = adBary[3];
      }
    } // Done checking all cells which have bboxes overlapping with the
      // vertex
    assert(iBestCell >= 0);
    pC = mesh->getCell(iBestCell);

    dEps = 1e-8;
    if (isInterpolatedLS() && dBestValue < (0.5 + dEps)) {
      //If a cell contains the vertex point, interpolate length scales from its vertices using barycentric coordinates. Otherwise, return LARGE_DBL.
      length_scale = 0;

      for (int ii = 0; ii < 4; ii++) {
	Vert* pV = pC->getVert(ii);
	length_scale += pV->getLengthScale() * adBestBary[ii];
      }
    }
    else {
      //Otherwise, use grading from the vertices of the closest cell ("closest" by barycentric coordinates, possibly not by euclidean distance).
      length_scale = LARGE_DBL;

      for (int ii = 0; ii < 4; ii++) {
	Vert* pV = pC->getVert(ii);
	length_scale = std::min(
	    length_scale,
	    pV->getLengthScale() + dDIST3D(coords, pV->getCoords()) / m_grade);
      }
    }

    return length_scale;
  }

  void
  Length3D::addToLengthMesh(double length, double x_loc, double y_loc,
			    double z_loc)
  {
    // Add a new point to the mesh.
    VolMesh* const pVM = dynamic_cast<VolMesh* const >(m_lengthMesh.get());
    // TODO:  This is okay if you only add one or two points, but if it's
    // a lot, this will kill you!
    // TODO:  A factory for inserters would clean this up immensely.
    GRUMMP::DelaunaySwapDecider3D SD3D(false);
    GRUMMP::SwapManager3D SM3D(&SD3D, pVM);
    GRUMMP::SwappingInserter3D SI3D(pVM, &SM3D);
    double coords[] =
      { x_loc, y_loc, z_loc };
    Vert *newVert = SI3D.insertPoint(coords, pVM->getCell(0));

    // Set the length scale, including grading effects
    newVert->setLengthScale(length);
    m_checkLength.insert(newVert);
    setGradedLength();
  }
} // end namespace GRUMMP

