/*
 * Length2D.cxx
 *
 *  Created on: Jun 27, 2017
 *      Author: cfog
 */

#include "GR_ADT.h"
#include "GR_InsertionManager.h"
#include "GR_Length.h"
#include "GR_Mesh2D.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Vec.h"

namespace GRUMMP
{

  Length2D::Length2D(const Mesh* const & pM, const double& refine,
		     const double& grading, const double& minLength) :
      Length(refine, grading, minLength)
  {
    const Mesh2D* const pM2D = dynamic_cast<const Mesh2D* const >(pM);
    assert(pM2D);
    m_lengthMesh = std::unique_ptr < Mesh > (new Mesh2D());
    Mesh2D* const pM2D_int = dynamic_cast<Mesh2D*>(m_lengthMesh.get());

    for (GR_index_t iV = 0; iV < pM->getNumVerts(); iV++) {
      Vert *pV = pM->getVert(iV);
      pM2D_int->createVert(pV->getCoords());
    }
    GR_index_t iC = 0;
    for (; iC < pM2D->getNumTriCells(); iC++) {
      Cell* pC = pM2D->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pM2D->getVertIndex(pC->getVert(0));
      int iV1 = pM2D->getVertIndex(pC->getVert(1));
      int iV2 = pM2D->getVertIndex(pC->getVert(2));
      pM2D_int->createTriCell(pM2D_int->getVert(iV0), pM2D_int->getVert(iV1),
			      pM2D_int->getVert(iV2));
    }
    for (; iC < pM2D->getNumCells(); iC++) {
      Cell* pC = pM2D->getCell(iC);
      if (pC->isDeleted())
	continue;
      int iV0 = pM2D->getVertIndex(pC->getVert(0));
      int iV1 = pM2D->getVertIndex(pC->getVert(1));
      int iV2 = pM2D->getVertIndex(pC->getVert(2));
      int iV3 = pM2D->getVertIndex(pC->getVert(3));
      pM2D_int->createQuadCell(pM2D_int->getVert(iV0), pM2D_int->getVert(iV1),
			       pM2D_int->getVert(iV2), pM2D_int->getVert(iV3));
    }
    for (GR_index_t f = 0; f < pM2D_int->getNumFaces(); f++) {
      Face *pF = pM2D_int->getFace(f);
      if (pF->isDeleted() || pF->getEntStatus() == Entity::ePartBdry)
	continue;
      if (pF->getLeftCell() == pCInvalidCell
	  || pF->getRightCell() == pCInvalidCell) {
	(void) pM2D_int->createBFace(pF, iDefaultBC);
      }
    }
    pM2D_int->identifyVertexTypesGeometrically();

    pM2D_int->purgeAllEntities();
    // Create a tree of the cells from that mesh to make interpolation
    // efficient
    createADTFromAllCells(pM2D_int, m_BBoxTree);
    initLength();
    const Mesh* const tmp = getMesh();
    // Copy this to the originating mesh.
    for (GR_index_t ii = 0; ii < pM->getNumVerts(); ii++) {
      pM->getVert(ii)->setLengthScale(tmp->getVert(ii)->getLengthScale());
    }
  }

  double
  Length2D::distToNearestBdryEntity(const Vert* const vert,
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
	  BdryEdgeBase* bedge = dynamic_cast<BdryEdgeBase*>(bdryCell);
	  assert(bedge);
	  double closest[] =
	    { 0, 0, 0 };
	  bedge->findClosestPointOnBFace(vert->getCoords(), closest);
	  double dist_squared = dDIST_SQ_3D(closest, vert->getCoords());
	  min_dist_squared = std::min(min_dist_squared, dist_squared);
	} // Yes, it's really a bdry face
      } // Done looping over edges of this cell
    } // Done looping over cells
    return sqrt(min_dist_squared);
  }

  double
  Length2D::queryLengthScale(const double coords[]) const
  {
    const Mesh* const mesh = m_lengthMesh.get();
    double length_scale = LARGE_DBL;
    std::vector<GR_index_t> veciQueryResult;

    // Find cell (nearly?) containing the target point
    double dEps = 1.e-8;
    do {
      double adBB[] =
	    { coords[0] - dEps, coords[0] + dEps, coords[1] - dEps, coords[1]
		+ dEps };

      veciQueryResult = m_BBoxTree->veciRangeQuery(adBB);
      dEps *= 10;
    }
    while (veciQueryResult.empty());

    double dBestValue = 1.e50;
    int iBestCell = -1;
    double adBestBary[3];
    Cell *pC;
    int iQ;
    GR_index_t iCell;
    for (iQ = veciQueryResult.size() - 1; iQ >= 0 && dBestValue > 0.5; iQ--) {
      // Get barycentric coords.  If best_value <(=) 0.5 on completion
      // then the vert lies within (on an edge of) the cell in list
      iCell = veciQueryResult[iQ];
      pC = mesh->getCell(iCell);

      double adBary[3];
      dynamic_cast<TriCell*>(pC)->calcBarycentricCoords(coords, adBary);

      double dThisValue = max(max(fabs(adBary[0] - 0.5), fabs(adBary[1] - 0.5)),
			      fabs(adBary[2] - 0.5));

      if (dThisValue < dBestValue) {
	iBestCell = iCell;
	dBestValue = dThisValue;
	adBestBary[0] = adBary[0];
	adBestBary[1] = adBary[1];
	adBestBary[2] = adBary[2];
      }
    } // Done checking all cells which have bboxes overlapping with the
      // vertex
    assert(iBestCell >= 0);
    pC = mesh->getCell(iBestCell);

    dEps = 1e-8;
    if (isInterpolatedLS() && dBestValue < (0.5 + dEps)) {
      //If a cell contains the vertex point, interpolate length scales from its vertices using barycentric coordinates. Otherwise, return LARGE_DBL.
      length_scale = 0;

      for (int ii = 0; ii < 3; ii++) {
	Vert* pV = pC->getVert(ii);
	length_scale += pV->getLengthScale() * adBestBary[ii];
      }
      assert(length_scale >= 0);
    }
    else {
      //Otherwise, use grading from the vertices of the closest cell ("closest" by barycentric coordinates, possibly not by euclidean distance).
      length_scale = LARGE_DBL;
      for (int ii = 0; ii < 3; ii++) {
	Vert* pV = pC->getVert(ii);
	length_scale = std::min(
	    length_scale,
	    pV->getLengthScale() + dDIST2D(coords, pV->getCoords()) / m_grade);
      }
      assert(length_scale >= 0);
    }

    return length_scale;
  }

  void
  Length2D::addToLengthMesh(double length, double x_loc, double y_loc,
			    double /*z_loc*/)
  {
    // Add a new point to the mesh.
    Mesh2D* const pM2D = dynamic_cast<Mesh2D* const >(m_lengthMesh.get());
    // TODO:  This is okay if you only add one or two points, but if it's
    // a lot, this will kill you!
    // TODO:  A factory for inserters would clean this up immensely.
    GRUMMP::DelaunaySwapDecider2D SD2D;
    GRUMMP::SwapManager2D SM2D(&SD2D, pM2D);
    GRUMMP::SwappingInserter2D SI2D(pM2D, &SM2D);
    SI2D.setForcedInsertion(true);
    double coords[] =
      { x_loc, y_loc };
    Vert *newVert = SI2D.insertPoint(coords, pM2D->getCell(0));

    // Set the length scale, including grading effects
    newVert->setLengthScale(length);
    m_checkLength.insert(newVert);
    setGradedLength();
  }
} // end namespace GRUMMP

