/*
 * CavityInserter.cxx
 *
 *  Created on: 2012-10-19
 *      Author: cfog
 */

#include <assert.h>
#include "GR_VolMesh.h"
#include "GR_Mesh2D.h"
#include "GR_InsertionManager.h"
namespace GRUMMP
{
  CavityInserter::CavityInserter(Mesh* mesh) :
      m_cellsToRemove(), m_bfacesToRemove(), m_facesOnHull(), m_hullInteriorFaces(), m_mesh(
	  mesh)
  {
  }

  void
  CavityInserter::recordRegions()
  {
    // Start over without any old data cluttering things up.
    m_hullFaceData.clear();
    faceIter_t iter, iterEnd = m_facesOnHull.end();
    int count = 0;
    for (iter = m_facesOnHull.begin(); iter != iterEnd; ++count, ++iter) {
      Face* face = *iter;
      if (face->isBdryFace() && face->getLeftCell()->isBdryCell()) {
	m_hullFaceData.push_back(
	    std::pair<Face*, int>(face, face->getRightCell()->getRegion()));
      }
      else {
	// This case applies both to interior faces and bdry faces with the
	// interior on their left side.
	m_hullFaceData.push_back(
	    std::pair<Face*, int>(face, face->getLeftCell()->getRegion()));
      }
    }

  }

  void
  CavityInserter::deleteHullInterior() const
  {
    // Remove all cells inside the hull

    {
      cellIter_t iter, iterEnd = m_cellsToRemove.end();
      for (iter = m_cellsToRemove.begin(); iter != iterEnd; ++iter) {
	m_mesh->deleteCell(*iter);
	assert((*iter)->isDeleted());
      }
    }

    // Faces are autoremoved.  This assertion loop is here so that
    // the code will crash if this behavior is ever changed.
#ifndef NDEBUG
    {
      faceIter_t iter, iterEnd = m_hullInteriorFaces.end();
      for (iter = m_hullInteriorFaces.begin(); iter != iterEnd; ++iter) {
	assert((*iter)->isDeleted());
      }
    }
#endif
  }

  double
  CavityInserter::dMinLengthOfNewFaces(double vert[])
  {
    //return the minimum length of all new created faces, used to reject points
    double minDist = 1.e10;
    double distance;
    recordRegions();
    for (int count = m_facesOnHull.size() - 1; count >= 0; --count) {
      Face* face = m_hullFaceData[count].first;
      for (int i = 0; i < 2; i++) {
	distance = sqrt(
	    pow((vert[0] - face->getVert(i)->x()), 2)
		+ pow((vert[1] - face->getVert(i)->y()), 2));
	minDist = std::min(minDist, distance);
      }
    }
    assert(minDist > 0.);
    return minDist;
  }

  void
  CavityInserter::connectHullToVert(Vert* const vert)
  {
    // Create new cells.

    {
      int minlayerindex = 1e9;
      int closefaceindex = 0;
      double distVtoF = 1e10;
      for (int count = m_facesOnHull.size() - 1; count >= 0; --count) {
	Face* face = m_hullFaceData[count].first;

	minlayerindex = std::min(minlayerindex, face->layerindex); //minlayerindex < face->layerindex ? minlayerindex: face->layerindex;

	const double *v = vert->getCoords();
	// TODO This clause is 2D centric.  Also, there's a non-trivial calculation
	// here that's just distance squared, which we have macros for.  CFOG 7/18.
	if (distVtoF
	    > pow(
		(v[0] - 0.5 * (face->getVert(0)->x() + face->getVert(1)->x())),
		2)
		+ pow(
		    (v[1]
			- 0.5 * (face->getVert(0)->y() + face->getVert(1)->y())),
		    2)) {
	  distVtoF =
	      pow((v[0] - 0.5 * (face->getVert(0)->x() + face->getVert(1)->x())),
		  2)
		  + pow(
		      (v[1]
			  - 0.5
			      * (face->getVert(0)->y() + face->getVert(1)->y())),
		      2);
	  closefaceindex = face->layerindex;
	}

      }
      if (closefaceindex > 1e4) {
	minlayerindex = 1e5 - 1;
      }

      for (int count = m_facesOnHull.size() - 1; count >= 0; --count) {
	Face* face = m_hullFaceData[count].first;
	if (face->hasVert(vert))
	  continue;
	int region = m_hullFaceData[count].second;
	m_mesh->createCell(minlayerindex, vert, face, region);
      }
    }
  }

  Vert*
  CavityInserter::insertPointInHull(const double candPtLoc[],
				    const Vert::VertType vertType,
				    const double splitParam)
  {

    // Create the new vertex.
    Vert* newVert = m_mesh->createVert(candPtLoc);
    assert(newVert->isValid());

    return insertPointInHull(newVert, vertType, splitParam);
  }

  Vert*
  CavityInserter::insertPointInHull(Vert* targetVert,
				    const Vert::VertType vertType,
				    const double splitParam)
  {
    assert(targetVert);
    recordRegions();
    deleteHullInterior();
    connectHullToVert(targetVert);
    updateMeshBdry(targetVert, splitParam);
    targetVert->setType(vertType);
    return targetVert;
  }
  bool
  CavityInserter::addCellToCavity(Cell* cellToAdd)
  {
    std::pair<cellIter_t, bool> ret = m_cellsToRemove.insert(cellToAdd);
    if (!ret.second)
      return false;

    for (int ii = cellToAdd->getNumFaces() - 1; ii >= 0; --ii) {
      Face* face = cellToAdd->getFace(ii);
      if (m_facesOnHull.count(face) == 0) {
	// If the face isn't on the surface of the cavity, then add it.
	m_facesOnHull.insert(face);
      }
      else {
	// If the face is on the surface of the cavity, then the fact
	// that we just found it again means it's in the interior.
	m_facesOnHull.erase(face);
	m_hullInteriorFaces.insert(face);
      }
    }
    return true;
  }

  void
  CavityInserter::addBFaceToCavity(BFace* bfaceToAdd)
  {
    m_bfacesToRemove.insert(bfaceToAdd);

    // Adding a bdry face doesn't make any changes to what's interior
    // to the hull, but we do need to change the definition of the
    // exterior, or we end up creating a cell for each bdry face that's
    // being replaced.  That's bad.
    for (int ii = bfaceToAdd->getNumFaces() - 1; ii >= 0; --ii) {
      Face* face = bfaceToAdd->getFace(ii);
      if (m_facesOnHull.count(face) != 0) {
	m_facesOnHull.erase(face);
      }
    }
  }
  void
  CavityInserter::fixIntBFacesInCavity()
  {
    bfaceIter_t iterBF, iterBFEnd = m_bfacesToRemove.end();
    for (iterBF = m_bfacesToRemove.begin(); iterBF != iterBFEnd; iterBF++) {
      BFace * bface = *iterBF;
      if (bface->getType() == Cell::eIntTriBFace) {
	if (isCellInCavity(bface->getFace(0)->getOppositeCell(bface))) {
	  Cell * oppCell = bface->getFace(1)->getOppositeCell(bface);
	  if (!isCellInCavity(oppCell)) {
	    addCellToCavity(oppCell);
	    m_facesOnHull.erase(bface->getFace(1));
	  }
	}
	if (isCellInCavity(bface->getFace(1)->getOppositeCell(bface))) {
	  Cell * oppCell = bface->getFace(0)->getOppositeCell(bface);
	  if (!isCellInCavity(oppCell)) {
	    addCellToCavity(oppCell);
	    m_facesOnHull.erase(bface->getFace(0));
	  }
	}

	for (int iF = 0; iF < bface->getNumFaces(); iF++)
	  if (m_facesOnHull.count(bface->getFace(iF)) == 1)
	    m_facesOnHull.erase(bface->getFace(iF));

      }
    }

  }

  void
  CavityInserter::makeCavityConvex(const double adPoint[])
  {
    // This has to be done a little carefully, so that facesOnHull isn't
    // changing as we iterate over it.  That gives undefined behavior: new
    // faces could be added to the list we're iterating over either before
    // or after the current iterator location, and if we delete the one we're
    // currently at in the iteration, the standard says the behavior is
    // undefined.  So do this in batch.
    size_t cellsAdded = 0;
    do {
      std::set<Cell*> cellsToAdd;
      for (Face * face : m_facesOnHull) {
	int iOrient = 0;
	if (face->getNumVerts() == 3)
	  iOrient = checkOrient3D(adPoint, face->getVert(0)->getCoords(),
				  face->getVert(1)->getCoords(),
				  face->getVert(2)->getCoords());
	else
	  iOrient = checkOrient2D(adPoint, face->getVert(0)->getCoords(),
				  face->getVert(1)->getCoords());
	switch (iOrient)
	  {
	  case 1:

	    // The left cell should be in the cavity; add if it isn't.
	    assert(face->getLeftCell()->getType() == CellSkel::eTet);
	    cellsToAdd.insert(face->getLeftCell());
	    break;
	  case -1:
	    // The right cell should be in the cavity; add if it isn't.
	    assert(face->getRightCell()->getType() == CellSkel::eTet);
	    cellsToAdd.insert(face->getRightCell());
	    break;
	  case 0:
	    assert(face->getLeftCell()->getType() == CellSkel::eTet);
	    assert(face->getRightCell()->getType() == CellSkel::eTet);

	    cellsToAdd.insert(face->getLeftCell());
	    cellsToAdd.insert(face->getRightCell());
	    break;
	  default:
	    // Can never get here.
	    assert(0);
	    break;
	  }

      }
      cellsAdded = 0;
      for (Cell* cell : cellsToAdd) {
	if (addCellToCavity(cell))
	  cellsAdded++;
      }
    }
    while (cellsAdded != 0);
    assert(isHullValid(adPoint));
  }

  bool
  CavityInserter::isHullValid(const double adPoint[])
  {
    faceIter_t iter;

    for (iter = m_facesOnHull.begin(); iter != m_facesOnHull.end(); ++iter) {
      Face * face = *iter;
      int iOrient = 0;
      if (face->getNumVerts() == 3)
	iOrient = checkOrient3D(adPoint, face->getVert(0)->getCoords(),
				face->getVert(1)->getCoords(),
				face->getVert(2)->getCoords());
      else
	iOrient = checkOrient2D(adPoint, face->getVert(0)->getCoords(),
				face->getVert(1)->getCoords());

      if (iOrient == 1) {
	if (!(!m_cellsToRemove.count(face->getRightCell())
	    && m_cellsToRemove.count(face->getLeftCell())))
	  return false;
      }
      else if (iOrient == -1) {
	if (!(m_cellsToRemove.count(face->getRightCell())
	    && !m_cellsToRemove.count(face->getLeftCell())))
	  return false;
      }
      else if (iOrient == 0)
	return false;
    }

    return true;
  }

  void
  CavityInserter2D::updateMeshBdry(Vert* const newVert,
				   const double splitParam) const
  {
    Mesh2D *mesh2D = dynamic_cast<Mesh2D*>(m_mesh);

    if (m_bfacesToRemove.empty())
      return;
    assert(m_bfacesToRemove.size() == 1);
    BdryEdgeBase *oldBFace =
	dynamic_cast<BdryEdgeBase*>(*(m_bfacesToRemove.begin()));
    assert(oldBFace);

    if (oldBFace->doesCurveExist()) {
      GRCurve* curve = oldBFace->getCurve();
      bool forward = oldBFace->isForward();

      Vert *vertA = forward ? oldBFace->getVert(0) : oldBFace->getVert(1);
      Vert *vertB = forward ? oldBFace->getVert(1) : oldBFace->getVert(0);
      double beginParam =
	  forward ? oldBFace->getVert0Param() : oldBFace->getVert1Param();
      double endParam =
	  forward ? oldBFace->getVert1Param() : oldBFace->getVert0Param();

      assert(iFuzzyComp(beginParam, endParam) == -1);
      assert(iFuzzyComp(beginParam, splitParam) == -1);
      assert(iFuzzyComp(splitParam, endParam) == -1);

      if (oldBFace->getType() == Cell::eIntBdryEdge) {

	mesh2D->createIntBdryEdge(vertA, newVert, curve, beginParam,
				  splitParam);
	mesh2D->createIntBdryEdge(newVert, vertB, curve, splitParam, endParam);
      }
      else {
	// need to be careful here
	// This handles the case of single boundary, and point that doesn't
	// lie inside any cell's circumcenter,
	if (m_cellsToRemove.empty()) {

	  Face * face = oldBFace->getFace(0);
	  int iReg = face->getOppositeCell(oldBFace)->getRegion();
	  // ugly hack, but...
	  face->removeCell(oldBFace);

	  face->setFaceLoc(Face::eInterior);
	  mesh2D->createCell(newVert, oldBFace->getFace(0), iReg);
	  oldBFace->removeFace(face);
	}

	mesh2D->createBdryEdge(vertA, newVert, curve, beginParam, splitParam);
	mesh2D->createBdryEdge(newVert, vertB, curve, splitParam, endParam);
      }

    }
    else {
      Vert *vertA = oldBFace->getVert(0);
      Vert *vertB = oldBFace->getVert(1);
      bool qIBE = false;
      if (oldBFace->getType() == Cell::eIntBdryEdge)
	qIBE = true;
      Face* faceA = findCommonFace(vertA, newVert, qIBE);
      Face* faceB = findCommonFace(vertB, newVert, qIBE);
      m_mesh->createBFace(faceA, oldBFace);
      m_mesh->createBFace(faceB, oldBFace);
    }
    mesh2D->deleteBFace(oldBFace);
  }

  void
  CavityInserter3D::updateMeshBdry(Vert* const newVert,
				   const double /*param*/) const
  {
    std::multiset<EdgeToBFace> sEBFAllEdges;

    bfaceIter_t iterBF, iterBFEnd = m_bfacesToRemove.end();
    for (iterBF = m_bfacesToRemove.begin(); iterBF != iterBFEnd; iterBF++) {
      BFace *pBF = *iterBF;
      // This assertion will have to be accompanied by one for the other
      // side for internal boundaries.
      Face *pF = pBF->getFace(0);
      assert(m_facesOnHull.count(pF) == 0);

      EdgeToBFace aEBF[3];
      aEBF[0] = EdgeToBFace(pF->getVert(0), pF->getVert(1), pBF);
      aEBF[1] = EdgeToBFace(pF->getVert(1), pF->getVert(2), pBF);
      aEBF[2] = EdgeToBFace(pF->getVert(2), pF->getVert(0), pBF);
      for (int ii = 0; ii < 3; ii++) {
	sEBFAllEdges.insert(aEBF[ii]);
      }
    }

    std::set<EdgeToBFace> sEBFSingleEdges;
    {
      std::multiset<EdgeToBFace>::iterator iterEBF, endEBF = sEBFAllEdges.end();
      for (iterEBF = sEBFAllEdges.begin(); iterEBF != endEBF; ++iterEBF) {
	if (sEBFAllEdges.count(*iterEBF) == 1)
	  sEBFSingleEdges.insert(*iterEBF);
      }
    }
    std::set<EdgeToBFace>::iterator iterEBF, endEBF = sEBFSingleEdges.end();
    for (iterEBF = sEBFSingleEdges.begin(); iterEBF != endEBF; ++iterEBF) {
      EdgeToBFace EBF = *iterEBF;
      Face *pF = findCommonFace(EBF.pVVert(0), EBF.pVVert(1), newVert, NULL,
				true);
      assert(pF->isValid());
      BFace *pBF = m_mesh->createBFace(pF, EBF.pBFace());
      assert(pBF->doFullCheck());
      assert(pF->doFullCheck());
    }

    // Now it's finally safe to delete the old BFace's.
    for (iterBF = m_bfacesToRemove.begin(); iterBF != iterBFEnd; iterBF++) {
      BFace* pBF = *iterBF;
      m_mesh->deleteBFace(pBF);
    }
  }
  void
  CavityInserter::writeCavityVTK(const char strBaseFileName[])
  {
    FILE *pFOutFile = NULL;
    char strFileName[1024];
    if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
	|| strlen(strBaseFileName) < 4) {
      // File name doesn't end in .vtk
      sprintf(strFileName, "%s.vtk", strBaseFileName);
    }
    else {
      // File name ends in .vtk
      sprintf(strFileName, "%s", strBaseFileName);
    }

    int m_dim = m_mesh->getType() == Mesh::eVolMesh ? 3 : 2;
    pFOutFile = fopen(strFileName, "w");
    if (NULL == pFOutFile)
      vFatalError("Couldn't open output file for writing",
		  "3d Vtk mesh output");
    GR_index_t nCells = m_cellsToRemove.size();
    GR_index_t nBFaces = m_bfacesToRemove.size();
    GR_index_t nno = nCells * (m_dim + 1) + nBFaces * m_dim;  // number of nodes
    //-------------------------------------------------------
    // Write the VTK header details
    //-------------------------------------------------------

    fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
    fprintf(pFOutFile, "GRUMMP Tetra example\n");
    fprintf(pFOutFile, "ASCII\n");
    fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(pFOutFile, "POINTS %d float\n", nno);

    //-------------------------------------
    // write 3d vertex data
    //-------------------------------------
    for (cellIter_t iter = m_cellsToRemove.begin();
	iter != m_cellsToRemove.end(); ++iter) {
      for (int i = 0; i < m_dim + 1; i++) {
	Vert *pV = (*iter)->getVert(i);
	fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
      }
    }
    for (bfaceIter_t iter = m_bfacesToRemove.begin();
	iter != m_bfacesToRemove.end(); ++iter) {
      for (int i = 0; i < m_dim; i++) {
	Vert *pV = (*iter)->getFace(0)->getVert(i);
	fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
      }
    }

    int counter = -1;
    fprintf(pFOutFile, "CELLS %d %d\n", nCells + nBFaces,
	    (m_dim + 2) * nCells + (m_dim + 1) * nBFaces);

    for (cellIter_t iter = m_cellsToRemove.begin();
	iter != m_cellsToRemove.end(); ++iter) {
      if (m_dim == 3) {
	fprintf(pFOutFile, "4 %d %d %d %d\n", counter + 1, counter + 2,
		counter + 3, counter + 4);
	counter = counter + 4;
      }
      else {
	fprintf(pFOutFile, "3 %d %d %d\n", counter + 1, counter + 2,
		counter + 3);
	counter = counter + 3;
      }
    }
    for (bfaceIter_t iter = m_bfacesToRemove.begin();
	iter != m_bfacesToRemove.end(); ++iter) {
      if (m_dim == 3) {
	fprintf(pFOutFile, "3 %d %d %d\n", counter + 1, counter + 2,
		counter + 3);
	counter = counter + 3;
      }
      else {
	fprintf(pFOutFile, "2 %d %d\n", counter + 1, counter + 2);
	counter = counter + 2;
      }
    }

    //-------------------------------------
    // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
    // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
    // VTK_PYRAMID = 14)
    //-------------------------------------
    fprintf(pFOutFile, "CELL_TYPES %d\n", nCells + nBFaces);
    for (cellIter_t iter = m_cellsToRemove.begin();
	iter != m_cellsToRemove.end(); ++iter) {
      if (m_dim == 3)
	fprintf(pFOutFile, "10\n");
      else
	fprintf(pFOutFile, "5\n");

    }
    for (bfaceIter_t iter = m_bfacesToRemove.begin();
	iter != m_bfacesToRemove.end(); ++iter) {
      if (m_dim == 3)
	fprintf(pFOutFile, "5\n");
      else
	fprintf(pFOutFile, "3\n");
    }
    fclose(pFOutFile);
  }
  void
  CavityInserter::writeHullVTK(const char strBaseFileName[])
  {
    FILE *pFOutFile = NULL;
    char strFileName[1024];
    if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
	|| strlen(strBaseFileName) < 4) {
      // File name doesn't end in .vtk
      sprintf(strFileName, "%s.vtk", strBaseFileName);
    }
    else {
      // File name ends in .vtk
      sprintf(strFileName, "%s", strBaseFileName);
    }
    int m_dim = m_mesh->getType() == Mesh::eVolMesh ? 3 : 2;
    pFOutFile = fopen(strFileName, "w");
    if (NULL == pFOutFile)
      vFatalError("Couldn't open output file for writing",
		  "3d Vtk mesh output");
    GR_index_t nFaces = m_facesOnHull.size();
    GR_index_t nno = nFaces * m_dim; // number of nodes
    //-------------------------------------------------------
    // Write the VTK header details
    //-------------------------------------------------------

    fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
    fprintf(pFOutFile, "GRUMMP Tetra example\n");
    fprintf(pFOutFile, "ASCII\n");
    fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(pFOutFile, "POINTS %d float\n", nno);

    //-------------------------------------
    // write 3d vertex data
    //-------------------------------------
    for (faceIter_t iter = m_facesOnHull.begin(); iter != m_facesOnHull.end();
	++iter) {
      for (int i = 0; i < m_dim; i++) {
	Vert *pV = (*iter)->getVert(i);
	fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
      }
    }

    int counter = -1;
    fprintf(pFOutFile, "CELLS %d %d\n", nFaces, (1 + m_dim) * nFaces);

    for (faceIter_t iter = m_facesOnHull.begin(); iter != m_facesOnHull.end();
	++iter) {
      if (m_dim == 3) {
	fprintf(pFOutFile, "3 %d %d %d\n", counter + 1, counter + 2,
		counter + 3);
	counter = counter + 3;
      }
      else {
	fprintf(pFOutFile, "2 %d %d\n", counter + 1, counter + 2);
	counter = counter + 2;
      }
    }

    //-------------------------------------
    // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
    // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
    // VTK_PYRAMID = 14)
    //-------------------------------------
    fprintf(pFOutFile, "CELL_TYPES %d\n", nFaces);

    for (faceIter_t iter = m_facesOnHull.begin(); iter != m_facesOnHull.end();
	++iter) {
      if (m_dim == 3)
	fprintf(pFOutFile, "5\n");
      else
	fprintf(pFOutFile, "3\n");
    }
    fclose(pFOutFile);
  }

}
