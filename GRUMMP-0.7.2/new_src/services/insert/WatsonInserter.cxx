#include <queue>
#include <set>

#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_InsertionManager.h"
#include "GR_SubsegManager.h"
#include "GR_Vertex.h"

namespace GRUMMP
{
  void
  WatsonInserter::computeHull(const double adPoint[], Cell* const seedCell,
			      BFace* const onBdryFace, const double splitParam)
  {
    m_CI->clear();
    m_bfacesEncroached.clear();

    m_candPtLoc[0] = adPoint[0];
    m_candPtLoc[1] = adPoint[1];
    m_candPtLoc[2] = adPoint[2];
    m_splitParam = splitParam;
    SimplexCell * SC = dynamic_cast<SimplexCell*>(seedCell);
    if (!seedCell->isValid()) {
      SC = findSeedNaively(adPoint);
      assert(SC->isValid());
    }

    std::queue<Cell*> cellsToTest;
    std::set<const Cell*> cellsTested;

    if (SC->isValid()) {
      assert(SC->circumscribesPoint(adPoint) || onBdryFace);

      cellsToTest.push(SC);
    }
    else {
      // If we get here, then seedCell must be valid, but not a simplex cell.
      assert(seedCell->isValid() && !SC->isValid());
      //For a boundary face, we must be extra cautious. The presence of
      //curves can cause the new vertex not to be included in the circumball
      //of the cell bounded by the boundary face. This is due to the fact that
      //we are actually redefining the convex hull of the set of points making
      //our triangulation. If this case happens, we may obtain an empty hull.
      assert(seedCell->getType() == Cell::eBdryEdge|| seedCell->getType() == Cell::eIntBdryEdge
	      || seedCell->getType() == Cell::eTriBFace
	      || seedCell->getType() == Cell::eIntTriBFace);

      cellsTested.insert(seedCell);

      for (int i = 0; i < seedCell->getNumFaces(); i++) {

	Face *face = seedCell->getFace(i);
	Cell *cell = face->getOppositeCell(seedCell);
	SC = dynamic_cast<SimplexCell*>(cell);
	// just test everything here
	if (SC) {
	  cellsToTest.push(cell);
	}
      }
    }
    // makes sure both sides of a boundary face get added to the queue
    if (onBdryFace)
      for (int i = 0; i < onBdryFace->getNumFaces(); i++)
	cellsToTest.push(onBdryFace->getFace(i)->getOppositeCell(onBdryFace));
    //The method performs a march from the top entry in the queue.
    //Marching out from seed_cell, it finds all the cells
    //having vertex's in its circumball and boundary
    //edges encroached by vertex.

    while (!cellsToTest.empty()) {

      Cell* cell = cellsToTest.front();
      cellsToTest.pop();
      // Don't bother if this cell has already been tested.
      if (!cellsTested.insert(cell).second)
	continue;

      SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);

      if (simplexCell->isValid()) {
	if (simplexCell->circumscribesPoint(adPoint)) {
	  m_CI->addCellToCavity(simplexCell);
	  for (int i = 0; i < simplexCell->getNumFaces(); i++) {
	    Face *face = simplexCell->getFace(i);
	    cellsToTest.push(face->getOppositeCell(simplexCell));
	  }
	}
      }
      else if (cell->isValid()) {
	BFace* bdryFace = dynamic_cast<BFace*>(cell);
	assert(bdryFace);

	if (m_encroachType != eNoEncroach) {
	  if (bdryFace->isPointEncroaching(adPoint, m_encroachType, true)) {
	    m_bfacesEncroached.insert(bdryFace);
	    if (bdryFace->isCoplanarWith(adPoint))
	      m_CI->addBFaceToCavity(bdryFace);

	  }
	}
      }
    }
    if (onBdryFace->isValid()) {
      // This call has to be at the end to remove whichever of its faces was
      // automatically added to the hull.

      // it is possible for the cells next to the boundary to not be in the Hull,
      // if the curves are such.
      for (int iF = 0; iF < onBdryFace->getNumFaces(); iF++) {
	SimplexCell* oppCell = dynamic_cast<SimplexCell*>((onBdryFace)->getFace(
	    iF)->getOppositeCell((onBdryFace)));
	if (!m_CI->isCellInCavity(oppCell)
	    && oppCell->circumscribesPoint(adPoint)) {
	  m_CI->addCellToCavity(oppCell);
	}
      }
      m_CI->addBFaceToCavity(onBdryFace);
    }
  }
/////////////////
/////////////////
  void
  WatsonInserter::FrontalComputeHull(const double adPoint[],
				     Cell* const seedCell,
				     BFace* const onBdryFace,
				     const double splitParam)
  {
    m_CI->clear();
    m_bfacesEncroached.clear();

    m_candPtLoc[0] = adPoint[0];
    m_candPtLoc[1] = adPoint[1];
    m_candPtLoc[2] = adPoint[2];
    m_splitParam = splitParam;
    SimplexCell * SC = dynamic_cast<SimplexCell*>(seedCell);
    if (!seedCell->isValid()) {
      SC = findSeedNaively(adPoint);
      assert(SC->isValid());
    }

    std::queue<Cell*> cellsToTest;
    std::set<const Cell*> cellsTested;

    if (SC->isValid()) {
//		assert( SC->circumscribesPoint(adPoint) || onBdryFace);

      cellsToTest.push(SC);
    }
    else {
      // If we get here, then seedCell must be valid, but not a simplex cell.
      assert(seedCell->isValid() && !SC->isValid());
      //For a boundary face, we must be extra cautious. The presence of
      //curves can cause the new vertex not to be included in the circumball
      //of the cell bounded by the boundary face. This is due to the fact that
      //we are actually redefining the convex hull of the set of points making
      //our triangulation. If this case happens, we may obtain an empty hull.
      assert(
	  seedCell->getType() == Cell::eBdryEdge
	      || seedCell->getType() == Cell::eIntBdryEdge
	      || seedCell->getType() == Cell::eTriBFace
	      || seedCell->getType() == Cell::eIntTriBFace);

      cellsTested.insert(seedCell);

      for (int i = 0; i < seedCell->getNumFaces(); i++) {

	Face *face = seedCell->getFace(i);
	Cell *cell = face->getOppositeCell(seedCell);
	SC = dynamic_cast<SimplexCell*>(cell);
	// just test everything here
	if (SC) {
	  cellsToTest.push(cell);
	}
      }
    }
    // makes sure both sides of a boundary face get added to the queue
    if (onBdryFace)
      for (int i = 0; i < onBdryFace->getNumFaces(); i++)
	cellsToTest.push(onBdryFace->getFace(i)->getOppositeCell(onBdryFace));
    //The method performs a march from the top entry in the queue.
    //Marching out from seed_cell, it finds all the cells
    //having vertex's in its circumball and boundary
    //edges encroached by vertex.

    while (!cellsToTest.empty()) {

      Cell* cell = cellsToTest.front();
      cellsToTest.pop();
      // Don't bother if this cell has already been tested.
      if (!cellsTested.insert(cell).second)
	continue;

      SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);

      if (simplexCell->isValid()) {
	if (simplexCell->circumscribesPoint(adPoint)) {
	  m_CI->addCellToCavity(simplexCell);
	  for (int i = 0; i < simplexCell->getNumFaces(); i++) {
	    Face *face = simplexCell->getFace(i);
	    cellsToTest.push(face->getOppositeCell(simplexCell));
	  }
	}
      }
      else if (cell->isValid()) {
	BFace* bdryFace = dynamic_cast<BFace*>(cell);
	assert(bdryFace);

	if (m_encroachType != eNoEncroach) {
	  if (bdryFace->isPointEncroaching(adPoint, m_encroachType, true)) {
	    m_bfacesEncroached.insert(bdryFace);
	    if (bdryFace->isCoplanarWith(adPoint))
	      m_CI->addBFaceToCavity(bdryFace);

	  }
	}
      }
    }
    if (onBdryFace->isValid()) {
      // This call has to be at the end to remove whichever of its faces was
      // automatically added to the hull.

      // it is possible for the cells next to the boundary to not be in the Hull,
      // if the curves are such.
      for (int iF = 0; iF < onBdryFace->getNumFaces(); iF++) {
	SimplexCell* oppCell = dynamic_cast<SimplexCell*>((onBdryFace)->getFace(
	    iF)->getOppositeCell((onBdryFace)));
	if (!m_CI->isCellInCavity(oppCell)
	    && oppCell->circumscribesPoint(adPoint)) {
	  m_CI->addCellToCavity(oppCell);
	}
      }
      m_CI->addBFaceToCavity(onBdryFace);
    }

  }

  double
  WatsonInserter::minLengthOfNewFaces(double newVert[])
  {
    return m_CI->dMinLengthOfNewFaces(newVert);
  }

/////////////////

  Vert*
  GRUMMP::WatsonInserter::insertPointInHull(Vert* newVert)
  {
    if (newVert) {
      // If a vertex is specified, it had better be in the right location.
      assert(newVert->x() == m_candPtLoc[0]);
      assert(newVert->y() == m_candPtLoc[1]);
      assert(newVert->z() == m_candPtLoc[2]);
      m_vertType = Vert::VertType(newVert->getVertType());
      return m_CI->insertPointInHull(newVert, m_vertType, m_splitParam);
    }
    else {
      // This is the case where a vertex isn't specified in the call.
      return m_CI->insertPointInHull(m_candPtLoc, m_vertType, m_splitParam);
    }
  }

  SimplexCell*
  WatsonInserter::findSeedNaively(const double newLoc[]) const
  {
    for (GR_index_t i = 0; i < getMesh()->getNumCells(); i++) {

      SimplexCell *cell = dynamic_cast<SimplexCell*>(getMesh()->getCell(i));

      if (cell->isDeleted())
	continue;
      if (getMesh()->getType() == Mesh::eMesh2D) {
	assert(cell->getType() == Cell::eTriCell);
	if (isInCircle(cell->getVert(0)->getCoords(),
		       cell->getVert(1)->getCoords(),
		       cell->getVert(2)->getCoords(), newLoc) > 0)
	  return cell;
      }
      else {
	assert(cell->getType() == Cell::eTet);
	if (isInsphere(cell->getVert(0)->getCoords(),
		       cell->getVert(1)->getCoords(),
		       cell->getVert(2)->getCoords(),
		       cell->getVert(3)->getCoords(), newLoc) > 0)
	  return cell;
      }
    }

    vFatalError("Unable to find a seed cell naively",
		"WatsonInserter::find_seed_naively");
    return NULL;
  }

  SimplexCell*
  WatsonInserter::findSeedFromGuesses(
      const double newLoc[], const std::set<SimplexCell*>* const guesses) const
  {
    if (guesses && !guesses->empty()) {
      std::set<SimplexCell*>::const_iterator it = guesses->begin(), it_end =
	  guesses->end();

      for (; it != it_end; ++it) {
	SimplexCell *cell = *it;
	if (cell->circumscribesPoint(newLoc))
	  return cell;
      }
    }
    return findSeedNaively(newLoc);
  }
  bool
  WatsonInserter::isHullValid()
  {
    return m_CI->isHullValid(m_candPtLoc);
  }
// WatsonInserter3D

  void
  WatsonInserter3D::compute3DHull(const double adPoint[], Cell* const seedCell,
				  BFace* const onBdryFace)
  {
    m_CI->clear();
    m_bfacesEncroached.clear();
    m_subsegsEncroached.clear();
    m_bfacesHit.clear();
    std::set<BFace*> bfacesToRemove;
    std::set<TetRefinerEdge> edges_visited;
    TopologyEntity * topEnt = NULL;
    if (onBdryFace)
      topEnt = onBdryFace->getTopologicalParent();
    m_candPtLoc[0] = adPoint[0];
    m_candPtLoc[1] = adPoint[1];
    m_candPtLoc[2] = adPoint[2];

    SimplexCell *SC = dynamic_cast<SimplexCell*>(seedCell);
    if (!seedCell->isValid()) {
      SC = findSeedNaively(adPoint);
      assert(SC->isValid());
    }

    std::queue<Cell*> cellsToTest;
    std::set<const Cell*> cellsTested;

    if (SC->isValid()) {
      assert(SC->circumscribesPoint(adPoint));
      cellsToTest.push(SC);
    }
    else {
      // If we get here, then seedCell must be valid, but not a simplex cell.
      assert(seedCell->isValid() && !SC->isValid());
      //For a boundary face, we must be extra cautious. The presence of
      //curves can cause the new vertex not to be included in the circumball
      //of the cell bounded by the boundary face. This is due to the fact that
      //we are actually redefining the convex hull of the set of points making
      //our triangulation. If this case happens, we may obtain an empty hull.
      assert(
	  seedCell->getType() == Cell::eBdryEdge
	      || seedCell->getType() == Cell::eIntBdryEdge
	      || seedCell->getType() == Cell::eTriBFace
	      || seedCell->getType() == Cell::eIntTriBFace);

      cellsTested.insert(seedCell);

      for (int i = 0; i < seedCell->getNumFaces(); i++) {

	Face *face = seedCell->getFace(i);
	Cell *cell = face->getOppositeCell(seedCell);
	SC = dynamic_cast<SimplexCell*>(cell);
	// just test everything here
	if (SC) {
	  cellsToTest.push(cell);
	}
      }
    }

    //The method performs a march from the top entry in the queue.
    //Marching out from seed_cell, it finds all the cells
    //having vertex's in its circumball and boundary
    //edges encroached by vertex.

    while (!cellsToTest.empty()) {

      Cell* cell = cellsToTest.front();
      cellsToTest.pop();
      // Don't bother if this cell has already been tested.
      if (!cellsTested.insert(cell).second)
	continue;

      SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);

      if (simplexCell->isValid()) {

	if (simplexCell->circumscribesPoint(adPoint)) {

	  m_CI->addCellToCavity(simplexCell);
	  for (int i = 0; i < simplexCell->getNumFaces(); i++) {
	    Face *face = simplexCell->getFace(i);
	    cellsToTest.push(face->getOppositeCell(simplexCell));
	  }
	}
	// Not sure what do here, but check for an encroached subseg(s), and store it
	hasEncroachedSubseg(cell, adPoint, &edges_visited);

      }
      else if (cell->isValid()) {
	BFace* bdryFace = dynamic_cast<BFace*>(cell);
	assert(bdryFace);

	m_bfacesHit.insert(bdryFace);
	if (m_encroachType != eNoEncroach) {

	  if (bdryFace->isPointEncroaching(adPoint, m_encroachType, true)) {
	    m_bfacesEncroached.insert(bdryFace);
	    if (bdryFace->isCoplanarWith(adPoint))
	      bfacesToRemove.insert(bdryFace);

	    if (onBdryFace)
	      if (topEnt == bdryFace->getTopologicalParent() && topEnt)
		bfacesToRemove.insert(bdryFace);

	  }
	  else if (onBdryFace)
	    if (topEnt == bdryFace->getTopologicalParent() && topEnt)
	      for (int iF = 0; iF < bdryFace->getNumFaces(); iF++)
		if (m_CI->isCellInCavity(
		    bdryFace->getFace(iF)->getOppositeCell(bdryFace)))
		  bfacesToRemove.insert(bdryFace);
	}
      }
    }
    if (onBdryFace->isValid()) {
      // This call has to be at the end to remove whichever of its faces was
      // automatically added to the hull.
      for (int iF = 0; iF < onBdryFace->getNumFaces(); iF++) {
	Cell* oppCell = (onBdryFace)->getFace(iF)->getOppositeCell(
	    (onBdryFace));
	if (!m_CI->isCellInCavity(oppCell))
	  m_CI->addCellToCavity(oppCell);

      }
      bfacesToRemove.insert(onBdryFace);
    }
    /// insert all the BFaces last. This is a safety step for internal boundaries..
    for (std::set<BFace*>::iterator iter = bfacesToRemove.begin();
	iter != bfacesToRemove.end(); ++iter)
      m_CI->addBFaceToCavity(*iter);

    m_CI->fixIntBFacesInCavity();

  }

  void
  WatsonInserter3D::computeSubsegHull(const double adPoint[],
				      std::vector<Cell*>& seed_cells,
				      std::vector<BFace*>& seed_bfaces)
  {

    m_CI->clear();
    m_bfacesEncroached.clear();
    m_subsegsEncroached.clear();
    std::set<BFace*> bfacesToRemove;

    std::set<TetRefinerEdge> edges_visited;

    m_candPtLoc[0] = adPoint[0];
    m_candPtLoc[1] = adPoint[1];
    m_candPtLoc[2] = adPoint[2];

    m_splitParam = -1;

    std::queue<Cell*> cellsToTest;
    std::set<const Cell*> cellsTested;

    std::vector<Cell*>::iterator it = seed_cells.begin(), it_end =
	seed_cells.end();

    for (; it != it_end; ++it) {
      cellsToTest.push(*it);
      m_CI->addCellToCavity(*it);
    }

    std::vector<BFace*>::iterator itb = seed_bfaces.begin(), itb_end =
	seed_bfaces.end();
    for (; itb != itb_end; ++itb) {
      for (int iF = 0; iF < (*itb)->getNumFaces(); iF++) {
	Cell* oppCell = (*itb)->getFace(iF)->getOppositeCell((*itb));
	cellsToTest.push(oppCell);
      }
    }
    //The method performs a march from the top entry in the queue.
    //Marching out from seed_cell, it finds all the cells
    //having vertex's in its circumball and boundary
    //edges encroached by vertex.
    // makes sure both sides of a boundary face get added to the queue

    while (!cellsToTest.empty()) {

      Cell* cell = cellsToTest.front();
      cellsToTest.pop();
      // Don't bother if this cell has already been tested.
      if (!cellsTested.insert(cell).second || cell->isDeleted())
	continue;

      SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);

      if (simplexCell->isValid()) {
	if (simplexCell->circumscribesPoint(adPoint)) {
	  m_CI->addCellToCavity(simplexCell);
	  for (int i = 0; i < simplexCell->getNumFaces(); i++) {
	    Face *face = simplexCell->getFace(i);
	    cellsToTest.push(face->getOppositeCell(simplexCell));
	  }
	}
	hasEncroachedSubseg(cell, adPoint, &edges_visited);

      }
      else if (cell->isValid()) {
	if (m_encroachType != eNoEncroach) {
	  BFace* bdryFace = dynamic_cast<BFace*>(cell);
	  assert(bdryFace);

	  if (bdryFace->isPointEncroaching(adPoint, m_encroachType, true)) {
	    m_bfacesEncroached.insert(bdryFace);
	    if (bdryFace->isCoplanarWith(adPoint))
	      bfacesToRemove.insert(bdryFace);

	    for (itb = seed_bfaces.begin(); itb != itb_end; ++itb)
	      if (((*itb)->getTopologicalParent()
		  == bdryFace->getTopologicalParent()
		  && bdryFace->getTopologicalParent() != NULL))
		bfacesToRemove.insert(bdryFace);

	  }
	  else {

	    for (itb = seed_bfaces.begin(); itb != itb_end; ++itb)
	      if (((*itb)->getTopologicalParent()
		  == bdryFace->getTopologicalParent()
		  && bdryFace->getTopologicalParent() != NULL))
		for (int iF = 0; iF < bdryFace->getNumFaces(); iF++)
		  if (m_CI->isCellInCavity(
		      bdryFace->getFace(iF)->getOppositeCell(bdryFace)))
		    bfacesToRemove.insert(bdryFace);

	  }
	}
      }
    }
    for (itb = seed_bfaces.begin(); itb != itb_end; ++itb) {
      for (int iF = 0; iF < (*itb)->getNumFaces(); iF++) {
	Cell* oppCell = (*itb)->getFace(iF)->getOppositeCell((*itb));
	if (!m_CI->isCellInCavity(oppCell))
	  m_CI->addCellToCavity(oppCell);
      }
      bfacesToRemove.insert(*itb);
    }
    for (std::set<BFace*>::iterator iter = bfacesToRemove.begin();
	iter != bfacesToRemove.end(); ++iter)
      m_CI->addBFaceToCavity(*iter);

    m_CI->fixIntBFacesInCavity();
  }
  bool
  WatsonInserter3D::hasEncroachedSubseg(Cell* const cell, const double newLoc[],
					std::set<TetRefinerEdge>* edges_visited)
  {

    assert(cell && cell->isValid() && !cell->isDeleted());
    assert(cell->getType() == Cell::eTet);
    assert(cell->getNumVerts() == 4);
    Vert* vert[2];
    Subseg* subseg;

    bool encroached_subseg = false;

    //Check the six edges for the presence of a subsegment.
    for (int i = 0; i < 6; ++i) {

      switch (i)
	{
	case 0:
	  vert[0] = cell->getVert(0);
	  vert[1] = cell->getVert(1);
	  break;
	case 1:
	  vert[0] = cell->getVert(0);
	  vert[1] = cell->getVert(2);
	  break;
	case 2:
	  vert[0] = cell->getVert(0);
	  vert[1] = cell->getVert(3);
	  break;
	case 3:
	  vert[0] = cell->getVert(1);
	  vert[1] = cell->getVert(2);
	  break;
	case 4:
	  vert[0] = cell->getVert(1);
	  vert[1] = cell->getVert(3);
	  break;
	case 5:
	  vert[0] = cell->getVert(2);
	  vert[1] = cell->getVert(3);
	  break;
	default:
	  assert(0);
	  break;
	}

      if (edges_visited->insert(TetRefinerEdge(vert[0], vert[1])).second) {

	subseg = dynamic_cast<VolMesh*>(getMesh())->areSubsegVertsConnected(
	    vert[0], vert[1]);

#ifndef NDEBUG
	if (subseg)
	  assert(
	      (vert[0]->getVertType() == Vert::eBdryApex||vert[0]->getVertType() == Vert::ePseudoSurface|| vert[0]->getVertType() == Vert::eBdryCurve)
		  && (vert[1]->getVertType() == Vert::eBdryApex || vert[1]->getVertType() == Vert::ePseudoSurface
		      || vert[1]->getVertType() == Vert::eBdryCurve));
#endif

	if (subseg
	    && (iFuzzyComp((CubitVector(newLoc) - subseg->mid_point()).length(),
			   subseg->ball_radius()) == -1)) {
	  m_subsegsEncroached.insert(subseg);
	  encroached_subseg = true;
	}
      }
    }
    return encroached_subseg;
  }

} // End of namespace GRUMMP
