#include "GR_SurfaceInsertion.h"

#include "GR_Mesh2D.h"
#include "GR_Geometry.h"
#include "GR_GRCurve.h"
#include "GR_GRCurveGeom.h"

namespace GRUMMP
{
  SurfaceInserter::SurfaceInserter(std::shared_ptr<Length> len, double grading,
				   double scale, const char strSurfParams[]) :
      m_isoLength(len), m_IPC(new CircumcenterIPC()), m_scale(scale), m_grading(
	  grading), m_anisoLength(nullptr), m_isAnisotropic(false), m_isStructured(
	  false), m_useOffsets(false), m_offsetByNormals(false), m_offsetCoordinate(
	  std::pair<int, CubitVector>(-1, CubitVector())), m_offsetPairs(), m_numberOfOffsets(
	  0), numInitVerts(0), numInsertedVerts(0), numRemovedVerts(0), numSmoothedVerts(
	  0), numRefinedVerts(0), numInitCells(0), m_leftoverVertsToRemove()
  {
    sprintf(m_params, "%s", strSurfParams);
  }
  ;
  void
  SurfaceInserter::setRegionFromCell(Cell * cell, const int region,
				     Mesh * /*mesh*/) const
  {
    // queue of cells to test, cells we've tested
    std::queue<Cell*> cellsToTest;

    // initialize the while loop
    cellsToTest.push(cell);
    while (!cellsToTest.empty()) {

      cell = cellsToTest.front();
      cellsToTest.pop();
      // Don't bother if this cell has already been tested.
      if (cell->getRegion() == region)
	continue;

      // ensures we only check tris/tets, not boundary objects
      if (cell->getType() == Cell::eTriCell
	  || cell->getType() == Cell::eQuadCell
	  || cell->getType() == Cell::eTet) {
	cell->setRegion(region);
	// solely for visualization/debugging
//#ifndef NDEBUG
//			if (mesh)
//				mesh->writeTempMesh();
//#endif

	for (int i = 0; i < cell->getNumFaces(); i++) {
	  Face *face = cell->getFace(i);
	  cellsToTest.push(face->getOppositeCell(cell));
	}
      }
    } // end while

  }

// counts moved vertices, only used for stats
// need to determine unique ones so count is correct
  void
  SurfaceInserter::receiveMovedVerts(std::vector<Vert*>& movedVerts)
  {
    logMessage(3, "SurfaceInsertion received word about %zu moved verts.\n",
	       movedVerts.size());
    std::vector<Vert*>::iterator iter = movedVerts.begin(), iterEnd =
	movedVerts.end();

    std::pair<std::set<Vert*>::iterator, bool> insertResult;

    for (; iter != iterEnd; ++iter) {
      Vert *pV = *iter;
      pV->markInactive();
      insertResult = vertQueue.insert(pV);
      if (insertResult.second)
	numSmoothedVerts++;

    }
  }

  void
  SurfaceInserter::initializeStats(Mesh* mesh)
  {

    numInitVerts = mesh->getNumVerts();
    numInitCells = mesh->getNumCells();
    /// To determine if a vert is changed, mark all initially as active
    /// Whatever is left active has not been deleted or smoothed
    /// not ideal, but active boolean already existed, and not used elsewhere
    for (GR_index_t i = 0; i < mesh->getNumVerts(); i++) {
      mesh->getVert(i)->markActive();
    }

#ifndef NDEBUG
    mesh->writeTempMesh();
#endif

    mesh->evaluateQuality();

  }
  void
  SurfaceInserter::finalizeStats(Mesh & mesh)
  {

    mesh.evaluateQuality(); //Final

    int numUnchangedVerts = 0;

    /// count all the active verts
    for (GR_index_t iV = 0; iV < mesh.getNumVerts(); iV++)
      numUnchangedVerts += mesh.getVert(iV)->isActive();

    // Any cell that has one vert changed is considered changed
    // Potentially wrong if swapping when all four verts unchanged, but
    // should never swap, since the initial mesh is of adequate quality
    int numUnchangedCells = 0;
    for (GR_index_t iC = 0; iC < mesh.getNumCells(); iC++) {
      bool qUnchanged = true;
      Cell* cell = mesh.getCell(iC);
      if (cell->isDeleted())
	continue;
      for (int iV = 0; iV < cell->getNumVerts(); iV++)
	if (!cell->getVert(iV)->isActive()) {
	  qUnchanged = false;
	  break;
	}
      numUnchangedCells += qUnchanged;
    }

    /// write all stats to screen and log file
    logMessage(1, "\nFinal Stats -------------\n");
    logMessage(1, "%d Initial Verts\n", numInitVerts);
    logMessage(1, "%d Removed Verts\n", numRemovedVerts);
    logMessage(1, "%d Inserted Verts\n", numInsertedVerts);
    logMessage(1, "%d Refined Verts\n", numRefinedVerts);
    logMessage(1, "%d Smoothed Verts\n", numSmoothedVerts);
    logMessage(1, "%d Final Verts\n", mesh.getNumVerts());
    logMessage(1, "%d Unchanged Verts\n", numUnchangedVerts);
    logMessage(1, "%d Initial Cells\n", numInitCells);
    logMessage(1, "%d Final Cells\n", mesh.getNumCells());
    logMessage(1, "%d Unchanged Cells\n", numUnchangedCells);
    logMessage(
	1, "Percentage of Initial Verts Left unchanged (approximate) %.2f%%\n",
	(numUnchangedVerts) * 100.0 / numInitVerts);
    logMessage(
	1, "Percentage of Initial Cells Left unchanged (approximate) %.2f%%\n",
	(numUnchangedCells) * 100.0 / numInitCells);

  }

  void
  SurfaceInserter::updateVertsToRemoveAfterInsertion(
      const std::map<Vert*, Vert*>& vertMap)
  {
    std::set<Vert*>::iterator itv;
    std::set<Vert*> newVertsToRemoveAfterInsertion, newLeftOverVerts;
    std::map<Vert*, Vert*>::const_iterator itm;
    for (itv = m_vertsToRemoveAfterInsertion.begin();
	itv != m_vertsToRemoveAfterInsertion.end(); ++itv) {
      itm = vertMap.find(*itv);
      if (itm != vertMap.end()) {
	if (!(itm->second->isDeleted()))
	  newVertsToRemoveAfterInsertion.insert((*itm).second);
      }
      else {
	if (!(itm->second->isDeleted()))
	  newVertsToRemoveAfterInsertion.insert(*itv);
      }
    }

    for (itv = m_leftoverVertsToRemove.begin();
	itv != m_leftoverVertsToRemove.end(); ++itv) {
      itm = vertMap.find(*itv);
      if (itm != vertMap.end()) {
	if (!(itm->second->isDeleted()))
	  newLeftOverVerts.insert((*itm).second);
      }
      else {
	newLeftOverVerts.insert(*itv);
      }
    }

    m_vertsToRemoveAfterInsertion.clear();
    m_vertsToRemoveAfterInsertion.insert(newVertsToRemoveAfterInsertion.begin(),
					 newVertsToRemoveAfterInsertion.end());
    m_leftoverVertsToRemove.clear();
    m_leftoverVertsToRemove.insert(newLeftOverVerts.begin(),
				   newLeftOverVerts.end());
    newVertsToRemoveAfterInsertion.clear();
    newLeftOverVerts.clear();
  }

}

// contains BaryLengthScale Class

Cell *
BaryLengthScale2D::findCell(const double * coords, Cell * guess)
{

  int numViolated = 0, numTies = 0;
  Face *face, *facesViolated[3], *facesTied[3];
  Cell * cellGuess = pCInvalidCell;
  if (guess)
    m_cellguess = guess;
  m_cellguess = m_mesh->findCell(coords, m_cellguess, numViolated,
				 facesViolated, numTies, facesTied);
  Cell * intBdryEdge = pCInvalidCell;
  /// Used to allow it to step over internal boundary
  while (isInTri(coords, m_cellguess->getVert(0)->getCoords(),
		 m_cellguess->getVert(1)->getCoords(),
		 m_cellguess->getVert(2)->getCoords()) < 0) {
    cellGuess = m_cellguess;
    for (int iF = 0; iF < 3; iF++) {
      Cell * oppCell = m_cellguess->getFace(iF)->getOppositeCell(m_cellguess);
      if (oppCell->getType() == Cell::eIntBdryEdge && oppCell != intBdryEdge) {
	intBdryEdge = oppCell;
	face =
	    (oppCell->getFace(0) == m_cellguess->getFace(iF)) ?
		oppCell->getFace(1) : oppCell->getFace(0);
	oppCell = face->getOppositeCell(oppCell);
	m_cellguess = m_mesh->findCell(coords, oppCell, numViolated,
				       facesViolated, numTies, facesTied);
	break;
      }
    }
    if (m_cellguess == cellGuess)
      break;
  }
  assert(m_cellguess);
  assert(
      isInTri(coords, m_cellguess->getVert(0)->getCoords(),
	      m_cellguess->getVert(1)->getCoords(),
	      m_cellguess->getVert(2)->getCoords()) > 0);
  return m_cellguess;
}

Cell *
BaryLengthScale2D::findCell(const CubitVector & point, Cell * guess)
{
  double coords[3];
  point.get_xyz(coords);
  return findCell(coords, guess);
}

double
BaryLengthScale2D::getLengthScale(const double & param, Cell * guess)
{

  if (m_aniso)
    return getAnisoLengthScale(param, guess);
  double baryCoords[3];
  double lengthscale = 0.0;
  double coords[3];

  CubitVector cv;
  m_curve->get_curve_geom()->coord_at_param(param, cv);
  cv.get_xyz(coords);

  Cell * containingCell = findCell(cv, guess);
  dynamic_cast<TriCell*>(containingCell)->calcBarycentricCoords(coords,
								baryCoords);

  for (int iV = 0; iV < 3; iV++) {
    // negatives allowed, due to boundary issues
    assert(iFuzzyComp(baryCoords[iV], -1e-6) > -1);
    lengthscale += baryCoords[iV]
	* containingCell->getVert(iV)->getLengthScale();
  }

  return lengthscale;
}
double
BaryLengthScale2D::getAnisoLengthScale(const double & param, Cell * guess,
				       CubitVector * dir)
{
  double baryCoords[3];
  double lengthscale = 0.0;
  double coords[3];

  CubitVector cv;
  m_curve->get_curve_geom()->coord_at_param(param, cv);
  cv.get_xyz(coords);

  Cell * containingCell = findCell(cv, guess);
  dynamic_cast<TriCell*>(containingCell)->calcBarycentricCoords(coords,
								baryCoords);

  // barycentric interpolation of the metric, which is probably stupid
  // and flawed given that metric space doesn't work that way when it
  // comes to measuring distances
  double metric[3] =
    { 0., 0., 0. };

  for (int iV = 0; iV < 3; iV++) {
    assert(iFuzzyComp(baryCoords[iV], -1e-6) > -1);
    for (int iM = 0; iM < 3; iM++)
      metric[iM] += baryCoords[iV] * containingCell->getVert(iV)->getMetric(iM);
  }
  CubitVector direction;
  if (dir) {
    direction.set(*dir);
  }
  else {
    m_curve->get_curve_geom()->unit_tangent(param, direction);
  }
  direction.normalize();

  // double disc = sqrt(
  // 		metric[0] * metric[0] - 2.0 * metric[0] * metric[2]
  // 				+ 4. * metric[1] * metric[1] + metric[2] * metric[2]) / 2.0;
//	double eig1 = metric[0]/2.0 + metric[2]/2.0 - disc;
//	double eig2 = metric[0]/2.0 + metric[2]/2.0 + disc;
//	double angle = atan(metric[1]/(metric[0]/2.0 - metric[2]/2.0 - disc));
//		double tvec = atan2(direction.y(),direction.x());
//			double tvec2 = atan2(sqrt(eig1)*sin(tvec-angle),sqrt(eig2)*cos(tvec-angle));
//			double lengthscale2 = sqrt(cos(tvec2)*eig1*cos(tvec2) + sin(tvec2)*eig2*sin(tvec2));

  /* No idea what was going on here. Removing it fixed a problem that popped up in surface insertion. */
//	direction.rotate(-angle,0);
//
//	if(fabs(direction.x())*sqrt(eig2) < fabs(direction.y())*sqrt(eig1) ){
//		direction.set(sqrt(eig2*direction.x()*direction.x()/(eig1*direction.y()*direction.y())),1.0,0.0);
//	} else {
//		direction.set(1.0,sqrt(eig1*direction.y()*direction.y()/(eig2*direction.x()*direction.x())),0.0);
//	}
//	direction.rotate(angle,0);
  /* Changed to reciprocal, which is the ideal length (in euclidean space) of an edge in the given direction. */
  lengthscale = 1.
      / sqrt(
	  metric[0] * direction.x() * direction.x()
	      + 2.0 * metric[1] * direction.x() * direction.y()
	      + metric[2] * direction.y() * direction.y());

//	printf("%f %f | %f %f ls\n",direction.x(),direction.y(),lengthscale,lengthscale2);
  return lengthscale;
}
void
BaryLengthScale2D::getMetric(const CubitVector & point, double metric[],
			     Cell * guess)
{
  assert(m_aniso);
  double baryCoords[3], coords[3];
  point.get_xyz(coords);
  Cell * containingCell = findCell(point, guess);
  dynamic_cast<TriCell*>(containingCell)->calcBarycentricCoords(coords,
								baryCoords);
  for (int iM = 0; iM < 3; iM++)
    metric[iM] = 0.0;
  for (int iV = 0; iV < 3; iV++) {
    assert(iFuzzyComp(baryCoords[iV], -1e-6) > -1);
    for (int iM = 0; iM < 3; iM++) {
      metric[iM] += baryCoords[iV] * containingCell->getVert(iV)->getMetric(iM);
    }
  }

}
double
BaryLengthScale2D::getAnisoLength(const CubitVector & pointA,
				  const CubitVector & pointB, Cell * guessA,
				  Cell * guessB)
{

  double coordsA[3], coordsB[3];
  pointA.get_xyz(coordsA);
  pointB.get_xyz(coordsB);
  double metricA[3] =
    { 0.0, 0.0, 0.0 };
  double metricB[3] =
    { 0.0, 0.0, 0.0 };

  getMetric(coordsA, metricA, guessA);
  getMetric(coordsB, metricB, guessB);

  double direction[2] =
    { coordsB[0] - coordsA[0], coordsB[1] - coordsA[1] };
  return sqrt(
      0.5 * (metricA[0] + metricB[0]) * direction[0] * direction[0]
	  + 0.5 * (metricA[2] + metricB[2]) * direction[1] * direction[1]
	  + (metricA[1] + metricB[1]) * direction[0] * direction[1]);
}

double
BaryLengthScale3D::getLengthScale(const double * coords)
{
  bool status;
  double baryCoords[4];
  double lengthscale = 0.0;
  m_cellguess = m_mesh->findCell(coords, m_cellguess, status);
  dynamic_cast<TetCell*>(m_cellguess)->calcBarycentricCoords(coords,
							     baryCoords);
  for (int iV = 0; iV < 4; iV++) {
    lengthscale += baryCoords[iV] * m_cellguess->getVert(iV)->getLengthScale();
  }
  return lengthscale;
}

double
BaryLengthScale3D::getLengthScale(const Vert * const vert)
{
  double result = -1;
  std::map<const Vert* const, double>::iterator iter = m_lengthLookup.find(
      vert);
  if (iter != m_lengthLookup.end()) {
    result = iter->second;
  }
  else {
    result = getLengthScale(vert->getCoords());
    m_lengthLookup.insert(std::make_pair(vert, result));
  }
  return result;
}
