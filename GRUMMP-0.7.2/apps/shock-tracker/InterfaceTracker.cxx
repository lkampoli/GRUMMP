/*
 * InterfaceTracker.cxx
 *
 *  Created on: Apr 12, 2017
 *      Author: cfog
 */

#include <algorithm>
#include <numeric>

#include "CGMApp.hpp"

#include "GR_ADT.h"
#include "GR_InterfaceTracker.h"
#include "GR_SurfaceInsertion.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_VolMesh.h"

namespace GRUMMP_Apps
{

  InterfaceTracker::InterfaceTracker() :
      m_compMesh(nullptr), m_surface(nullptr), m_length(nullptr), m_surfToVolMap(), m_numVars(
	  0), m_solnData(nullptr)
  {
    int dummyInt = 0;
    char *dummyCharArrays;
    CGMApp::instance()->startup(dummyInt, &dummyCharArrays);

    // Quality: all dihedral angles
    // Max surface angle for swapping: 5 degrees from flat
    m_compMesh = new VolMesh(2, 5);

    m_pSI3D = new GRUMMP::SwappingInserter3D(m_compMesh, nullptr);

    // Tell the mesh not to actually delete from the database
    // verts that have been removed from the mesh.
    m_compMesh->keepPhantomVerts(true);
  }

  InterfaceTracker::~InterfaceTracker()
  {
    if (m_compMesh)
      delete m_compMesh;
    if (m_solnData)
      delete[] m_solnData;
  }

  void
  InterfaceTracker::readVolumeMesh(const char volumeFilename[])
  {
    m_compMesh->readFromFile(volumeFilename, &m_numVars, &m_solnData);
    assert(m_compMesh->isValid());
    m_length = GRUMMP::Length::Create(m_compMesh, 1, 1);
  }

  void
  InterfaceTracker::writeVolumeMesh()
  {
    writeVTKLegacy(*m_compMesh, "with_surf");
    writeFileTetGen(*m_compMesh, "with_surf", "", m_numVars, m_solnData);
  }

  void
  InterfaceTracker::splitInterface()
  {
    // For now, let's stick with only simplicial meshes.
    assert(m_compMesh->isSimplicial());

    std::map<Vert*, Vert*> vertMap;

    // Identify all verts on the internal bdry faces
    // Assumption: all internal bdry faces are on shocks, not
    // on other types of internal bdrys.  This could be relaxed
    // in the future without too much trouble.
    std::set<Vert*> shockVerts;
    GR_index_t lastIBF = m_compMesh->getNumIntBdryFaces();
    for (GR_index_t iIntBF = 0; iIntBF < lastIBF; iIntBF++) {
      BFace* pIBF = m_compMesh->getIntBFace(iIntBF);
      for (int ii = 0; ii < pIBF->getNumVerts(); ii++) {
	shockVerts.insert(pIBF->getVert(ii));
      }
    }

    // Create a place to put the re-ordering for the vertices that
    // gets all the shock points put together in groups properly.
    // This is actually slightly oversized.
    GR_index_t numShockVerts = shockVerts.size();
    GR_index_t numOrigVerts = m_compMesh->getNumVerts() - numShockVerts;
    GR_index_t numShockPerimVerts = 0, numShockRegularVerts = 0;
    GR_index_t numShockInteractionVerts = 0;
    GR_index_t numTotalVerts = numOrigVerts;
    std::set<Vert*> setPerimVerts;
    std::map<Vert*, Vert*> shockRegVertMap;
    std::vector<Vert*> shockInteractionVerts;
    // Now we've got all the verts in the set.  For each one, split the
    // neighborhood into distinct traversable pieces.  Each piece after the first
    // must get its own copy of the vertex (and faces must have verts replaced).
    std::set<Vert*>::iterator iter, iterEnd = shockVerts.end();
    std::set<Face*> allIncFaces;
    std::vector<std::set<Face*> > faceGroups;
    for (iter = shockVerts.begin(); iter != iterEnd; ++iter) {
      Vert* currVert = *iter;
      faceGroups.clear();
      allIncFaces.clear();
      GR_index_t facesRemaining = currVert->getNumFaces();
      for (GR_index_t ii = 0; ii < facesRemaining; ii++) {
	allIncFaces.insert(currVert->getFace(ii));
      }
      GR_index_t facesHitSoFar = 0;
      // Now pick a face and traverse the mesh, without crossing bdry faces,
      // finding all faces incident on the current vertex.
      while (facesRemaining) {
	assert(allIncFaces.size() == facesRemaining);
	std::set<Face*> facesToCheck;
	std::set<Face*> newFaceGroup;
	Face* firstFace = *(allIncFaces.begin());
	facesToCheck.insert(firstFace);
	newFaceGroup.insert(firstFace);
	allIncFaces.erase(firstFace);
	while (!facesToCheck.empty()) {
	  std::set<Face*>::iterator iter2 = facesToCheck.begin();
	  Face *tmpFace = *iter2;
	  facesToCheck.erase(iter2); // iter2 now invalid.
	  for (int iC = 0; iC < 2; iC++) {
	    Cell *pC = tmpFace->getCell(iC);
	    // For interior cells, check all faces.
	    if (pC->isValid() && !pC->isBdryCell()) {
	      for (int ii = 0; ii < pC->getNumFaces(); ii++) {
		Face *candFace = pC->getFace(ii);
		// Don't add faces that aren't incident on this vertex.
		if (!candFace->hasVert(currVert))
		  continue;
		bool wasInserted = (newFaceGroup.insert(candFace)).second;
		// If the face hasn't been visited before, we need to
		// check its neighborhood.  Inefficient, because it
		// re-checks the cell we just visited.
		if (wasInserted) {
		  // This face shouldn't already be cued to check, I
		  // don't think.  If the assertion holds, then we
		  // can change this to a queue instead of a set. TODO
		  bool result = (facesToCheck.insert(candFace)).second;
		  assert(result);
		  allIncFaces.erase(candFace);
		}
	      }
	    } // Checked faces of this (interior) cell
	  } // Checked both cells incident on the face
	} // Done checking all faces in this traversable piece of the neighborhood
	GR_index_t facesHere = newFaceGroup.size();
	assert(facesHere <= facesRemaining);
	facesHitSoFar += facesHere;
	facesRemaining -= facesHere;
	faceGroups.push_back(newFaceGroup);
      } // Done checking all faces in the neighborhood
      if (faceGroups.size() == 1) {
	// This point is on the perimeter of the shock interface, so there's
	// no duplication required; just tag it.
	currVert->setType(Vert::eShockPerimeter);
	setPerimVerts.insert(currVert);
	numShockPerimVerts++;
	numTotalVerts++;
      }
      else {
	// Mark all the others, and the duplicates, as shock upstream for
	// now, unless there are more than two regions; in that case, we've
	// got some sort of triple (or more) point.
	// Ordinary downstream shock points will be relabeled later.
	size_t numGroups = faceGroups.size();
	bool isRegular = true;
	if (numGroups == 2) {
	  currVert->setType(Vert::eShockUpstream);
	  numShockRegularVerts++;
	  numTotalVerts += 2;
	}
	else {
	  isRegular = false;
	  currVert->setType(Vert::eShockInteraction);
	  shockInteractionVerts.push_back(currVert);
	  numTotalVerts += numGroups;
	  numShockInteractionVerts++;
	}

	for (size_t ii = 1; ii < numGroups; ii++) {
	  Vert* newVert = m_compMesh->duplicateVert(currVert);
	  std::set<Face*> thisFaceGroup = faceGroups[ii];
	  std::set<Face*>::iterator iter2, tFG_end = thisFaceGroup.end();
	  for (iter2 = thisFaceGroup.begin(); iter2 != tFG_end; ++iter2) {
	    Face *tmpFace = *iter2;
	    tmpFace->replaceVert(currVert, newVert);
	    currVert->removeFace(tmpFace);
	    newVert->addFace(tmpFace);
	  } // Replaced the vertex for all faces in this group.
	  if (isRegular) {
	    shockRegVertMap.insert(make_pair(currVert, newVert));
	  }
	  else {
	    shockInteractionVerts.push_back(newVert);
	    numShockInteractionVerts++;
	  }
	} // Done with all groups.
      } // Done with points that are in the interior of a shock surface
    } // Done with all processing for this shock point.

    logMessage(MSG_DETAIL, "Splitting %u shock faces\n", lastIBF);
    for (GR_index_t iIntBF = 0; iIntBF < lastIBF; iIntBF++) {
      BFace* pIBF = m_compMesh->getIntBFace(iIntBF);
      Face *face0 = pIBF->getFace(0);
      if (face0->getRightCell() == pIBF) {
	face0->interchangeCellsAndVerts();
      }
      assert(face0->getLeftCell() == pIBF);
      double normalFace0[3];
      face0->calcUnitNormal(normalFace0);

      Face *face1 = pIBF->getFace(1);
      if (face1->getRightCell() == pIBF) {
	face1->interchangeCellsAndVerts();
      }
      assert(face1->getLeftCell() == pIBF);
      double normalFace1[3];
      face1->calcUnitNormal(normalFace1);
      assert(dDOT3D(normalFace0, normalFace1) < -0.9999);

      // Find the normal of the underlying surface at the BFace centroid.
      RefFace* pRF = dynamic_cast<RefFace*>(pIBF->getTopologicalParent());
      assert(pRF);
      double centroid[] =
	{
	    (pIBF->getVert(0)->x() + pIBF->getVert(1)->x()
		+ pIBF->getVert(2)->x()) / 3., (pIBF->getVert(0)->y()
		+ pIBF->getVert(1)->y() + pIBF->getVert(2)->y()) / 3.,
	    (pIBF->getVert(0)->z() + pIBF->getVert(1)->z()
		+ pIBF->getVert(2)->z()) / 3. };
      CubitVector location(centroid);
      CubitVector normalSurfVec = pRF->normal_at(location, NULL);
      double normalSurface[] =
	{ normalSurfVec.x(), normalSurfVec.y(), normalSurfVec.z() };
      // Dot product: positive is downstream, negative is upstream.
      double dot0 = dDOT3D(normalFace0, normalSurface);
      logMessage(MSG_DEBUG, "dot: %f\n", dot0);
//		assert(fabs(dot0) > 0.99); // Better be pretty closely aligned.
      int iBC0 = Vert::eShockUpstream, iBC1 = Vert::eShockUpstream; // The downstream one will have to be changed.
      if (dot0 > 0) {
	// face 0 is the one on the downstream side
	iBC0 = Vert::eShockDownstream;
	for (int ii = 0; ii < face0->getNumVerts(); ii++) {
	  Vert *pV = face0->getVert(ii);
	  if (pV->getVertType() == Vert::eShockUpstream) {
	    pV->setType(Vert::eShockDownstream);
	  }
	}
      }
      else {
	// face 1 is the one on the downstream side
	iBC1 = Vert::eShockDownstream;
	for (int ii = 0; ii < face1->getNumVerts(); ii++) {
	  Vert *pV = face1->getVert(ii);
	  if (pV->getVertType() == Vert::eShockUpstream) {
	    pV->setType(Vert::eShockDownstream);
	  }
	}
      }
      // Now kill the int bdry face and create two new regular bdry faces.
      m_compMesh->deleteBFace(pIBF);
      BFace* bdryFace0 = m_compMesh->createBFace(face0);
      bdryFace0->setBdryCondition(iBC0);

      BFace* bdryFace1 = m_compMesh->createBFace(face1);
      bdryFace1->setBdryCondition(iBC1);
      logMessage(MSG_DEBUG, "iBF: %d.  BC's: %d %d  Pointers: %p %p\n", iIntBF,
		 bdryFace0->getBdryCondition(), bdryFace1->getBdryCondition(),
		 bdryFace0, bdryFace1);

      // It's not clear whether I want these to be associated with the curved
      // surface; for now, I'm not going to.
    }
    logMessage(MSG_DEBUG, "After initial split: ");
    outputNumberOfShockFaces();

    // Now set up for re-numbering.  Original verts are left alone.
    assert(numShockPerimVerts == setPerimVerts.size());
    assert(numShockRegularVerts == shockRegVertMap.size());
    assert(numShockInteractionVerts == shockInteractionVerts.size());
    assert(
	numOrigVerts + 2 * numShockRegularVerts + numShockInteractionVerts
	    + numShockPerimVerts
	    == m_compMesh->getNumVerts());
    GR_index_t *order = new GR_index_t[numTotalVerts];
    std::iota(order, order + numOrigVerts, 0);
    GR_index_t count;
    for (std::pair<Vert*, Vert*> SVMpair : shockRegVertMap) {
      // Get both vertex indices; keep it straight which is upstream
      Vert *pVUp = SVMpair.first, *pVDown = nullptr;
      if (pVUp->getVertType() == Vert::eShockUpstream) {
	pVDown = SVMpair.second;
      }
      else {
	pVDown = pVUp;
	pVUp = SVMpair.second;
      }
      assert(
	  (pVUp->getVertType() == Vert::eShockUpstream
	      && pVDown->getVertType() == Vert::eShockDownstream)
	      || (pVUp->getVertType() == Vert::eShockInteraction
		  && pVDown->getVertType() == Vert::eShockInteraction));
      GR_index_t upIndex = m_compMesh->getVertIndex(pVUp);
      GR_index_t downIndex = m_compMesh->getVertIndex(pVDown);
      order[numOrigVerts + count] = upIndex;
      order[numOrigVerts + numShockRegularVerts + count] = downIndex;
      count++;
    }

    count = 0;
    for (Vert* pV : setPerimVerts) {
      GR_index_t perimIndex = m_compMesh->getVertIndex(pV);
      order[numOrigVerts + 2 * numShockRegularVerts + count] = perimIndex;
      count++;
    }

    count = 0;
    for (Vert* pV : shockInteractionVerts) {
      GR_index_t interIndex = m_compMesh->getVertIndex(pV);
      order[numOrigVerts + 2 * numShockRegularVerts + numShockPerimVerts + count] =
	  interIndex;
      count++;
    }

    logMessage(MSG_DEBUG, "About to re-order\n");
    m_compMesh->reorderVertices(order);
    delete[] order;

    // Need to expand data set.
    if (m_solnData && m_numVars > 0) {
      GR_index_t numVerts = m_compMesh->getNumVerts();
      double *newData = new double[numVerts * m_numVars];
      std::fill(newData, newData + numVerts * m_numVars, -1.e100);
      std::copy(m_solnData, m_solnData + numOrigVerts * m_numVars, newData);
      delete[] m_solnData;
      m_solnData = newData;
      for (GR_index_t ii = numOrigVerts; ii < numVerts; ii++) {
	m_solnData[(ii + 1) * m_numVars - 1] = ii + 1;
      }
    }
    assert(m_compMesh->isValid());
  }

  void
  InterfaceTracker::swapNearShocks(GRUMMP::SwapManager3D& SM3D)
  {
    // Now fill the queue
    for (GR_index_t iV = 0; iV < m_compMesh->getNumVerts(); iV++) {
      Vert* pV = m_compMesh->getVert(iV);
      if (pV->isDeleted() || pV->getVertType() < Vert::eShockUpstream) {
	continue;
      }
      m_compMesh->moveVertEvent(pV);
    }
    m_compMesh->sendEvents();
    int numSwaps = SM3D.swapAllQueuedFacesRecursive();
    logMessage(MSG_MANAGER, "Swapped %d times for mesh quality\n", numSwaps);
  }

  void
  InterfaceTracker::improveMesh()
  {
    // Now swap near the new surfaces to improve mesh quality
    GRUMMP::DelaunaySwapDecider3D SDDel(false);
    GRUMMP::MaxMinSineSwapDecider3D SDSine(false);
    GRUMMP::SwapManager3D SM3D(&SDDel, m_compMesh);
    swapNearShocks(SM3D);
    SM3D.changeSwapDecider(&SDSine);
    swapNearShocks(SM3D);
  }

  void
  InterfaceTracker::addInterfaceToMesh(const char surfaceFilename[])
  {
    // Set up length scale

    // Need some options to control behavior of surface insertion, most
    // notably not actually deleting phantom nodes.
    GRUMMP::SurfaceInserter3D SurfIns3D(m_compMesh, m_length, 1, 1, "");
    SurfIns3D.insertSurfaceInMesh(surfaceFilename, "");

    splitInterface();
    logMessage(MSG_DEBUG, "After splitInterface: ");

    // Now swap near the new surfaces to improve mesh quality
    improveMesh();

    outputNumberOfShockFaces();
    writeVTKLegacyWithoutPurge(*m_compMesh, "raw_after");
    logMessage(MSG_DEBUG, "After writeWithoutPurge: ");
    outputNumberOfShockFaces();

    assert(m_compMesh->isValid());
  }

  GR_index_t
  InterfaceTracker::outputNumberOfShockFaces()
  {
    GR_index_t nUp = 0, nDown = 0;
    for (GR_index_t iBF = 0; iBF < m_compMesh->getNumBdryFaces(); iBF++) {
      int BC = m_compMesh->getBFace(iBF)->getBdryCondition();
      switch (BC)
	{
	case Vert::eShockUpstream:
	  nUp++;
	  break;
	case Vert::eShockDownstream:
	  nDown++;
	  break;
	default:
	  break;
	}
    }
    assert(nUp == nDown);
    logMessage(MSG_SERVICE, " Shock faces: %u up, %u down\n", nUp, nDown);
    return (nUp);
  }

  void
  InterfaceTracker::updateMesh()
  {
    logMessage(MSG_MANAGER,
	       "Updating mesh connectivity after shock movement.\n");
    // Create a search tree for all shock faces.  Well, actually, only the
    // upstream shock faces, because having the downstream ones as well is
    // redundant for this purpose.

    // First, count the number of upstream shock faces.
    GR_index_t nShockFaces = outputNumberOfShockFaces();

    // Set up an array for face bounding box data.  This needs to be
    // 6 x nShockFaces (min/max in x, y, z).
    double (*boundingBoxes)[6] = new double[nShockFaces][6];
    BFace **bdryFaces = new BFace*[nShockFaces];
    GR_index_t sf = 0;
    for (GR_index_t iBF = 0; iBF < m_compMesh->getNumBdryFaces(); iBF++) {
      BFace* pBF = m_compMesh->getBFace(iBF);
      // Only bother if it's an upstream shock face
      if (pBF->getBdryCondition() == Vert::eShockUpstream) {
	boundingBoxes[sf][0] = boundingBoxes[sf][2] = boundingBoxes[sf][4] =
	LARGE_DBL; // Mins
	boundingBoxes[sf][1] = boundingBoxes[sf][3] = boundingBoxes[sf][5] =
	    -LARGE_DBL; // Maxes
	for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
	  const double *coords = pBF->getVert(iV)->getCoords();
	  boundingBoxes[sf][0] = std::min(boundingBoxes[sf][0], coords[0]);
	  boundingBoxes[sf][1] = std::max(boundingBoxes[sf][1], coords[0]);
	  boundingBoxes[sf][2] = std::min(boundingBoxes[sf][2], coords[1]);
	  boundingBoxes[sf][3] = std::max(boundingBoxes[sf][3], coords[1]);
	  boundingBoxes[sf][4] = std::min(boundingBoxes[sf][4], coords[2]);
	  boundingBoxes[sf][5] = std::max(boundingBoxes[sf][5], coords[2]);
	}
	bdryFaces[sf] = pBF;
	sf++;
      }
    }
    assert(sf == nShockFaces);
    logMessage(MSG_SERVICE, "Building face search tree...");
    ADT shockFaceTree(nShockFaces, 3, ADT::eBBoxes,
		 reinterpret_cast<double*>(boundingBoxes));
    logMessage(MSG_SERVICE, "done.\n");
    delete[] boundingBoxes;

    // Find all (non-phantom) verts that are closer to a front face than
    // m_minDist times length scale at the vertex; these are going to
    // become phantom verts.

    // Find all (phantom) verts that are farther from a front face than
    // m_maxDist times length scale at the vertex; these are going to be
    // re-added as real computational vertices.  Also, solution will be
    // iterpolated onto these vertices.

    int numNewPhantoms = 0, numAttemptedPhantoms = 0;
    int numNewReals = 0, numAttemptedReals = 0;
    std::set<Vert*> attemptedReals, attemptedPhantoms;
    ADT *CellTree;
    createADTFromAllCells(m_compMesh, CellTree);
    int totalSwaps = 0;
    for (GR_index_t iV = 0; iV < m_compMesh->getNumVerts(); iV++) {
      if ((iV + 1) % 10000 == 0) {
	logMessage(MSG_MANAGER, "Checking vert %d of %d\r", iV + 1,
		   m_compMesh->getNumVerts());
      }
      Vert *pV = m_compMesh->getVert(iV);
      unsigned int type = pV->getVertType();
      switch (type)
	{
	case Vert::eShockUpstream:
	case Vert::eShockDownstream:
	case Vert::eShockInteraction:
	case Vert::eShockPerimeter:
	  // These don't need to be checked.
	  break;
	case Vert::eBdryTwoSide:
	case Vert::eBdryApex:
	case Vert::eBdryCurve:
	case Vert::eInteriorFixed:
	case Vert::ePseudoSurface:
	  // These may be extremely difficult to handle.
	  break;
	case Vert::eBdry:
	  // These are feasible to handle, but don't try yet.
	  break;
	case Vert::eInterior:
	  if (isCloseToShock(pV, shockFaceTree, bdryFaces, m_minDist)) {
	    attemptedPhantoms.insert(pV);
	  }
	  break;
	case Vert::ePhantom:
	  if (!isCloseToShock(pV, shockFaceTree, bdryFaces, m_maxDist)) {
	    attemptedReals.insert(pV);
	  }
	  break;
	default:
	  // These should be impossible.
	  assert(0);
	  break;
	}
    }

    // Before making any changes to the mesh, interpolate data to the phantom
    // points that are about to be re-added.
    for (Vert* pV : attemptedReals) {
      double adBary[4];
      Cell *pC = findCellContainingPoint(*m_compMesh, CellTree, pV->getCoords(),
					adBary);
      // Zero out the data for this vertex.
      int iVert = m_compMesh->getVertIndex(pV);
      for (int ii = 0; ii < m_numVars; ii++) {
	m_solnData[iVert*m_numVars + ii] = 0;
	for (int iii = 0; iii < 4; iii++) {
	  int iThisVert = m_compMesh->getVertIndex(pC->getVert(iii));
	  m_solnData[iVert*m_numVars + ii] += adBary[iii] * m_solnData[iThisVert*m_numVars + ii];
	}
      }
    }

    // Actually remove verts that should be phantoms
    for (Vert* pV : attemptedPhantoms) {
      int newSwaps = 0;
      bool success = makeVertexPhantom(pV, newSwaps);
      totalSwaps += newSwaps;
      if (success) {
	numNewPhantoms++;
      }
      numAttemptedPhantoms++;
    }
    // Actually re-add verts that should no longer be phantoms
    for (Vert* pV : attemptedReals) {
      bool success = makePhantomReal(pV);
      if (success) {
	numNewReals++;
      }
      numAttemptedReals++;

    }

    logMessage(MSG_MANAGER, "%d new phantoms out of %d attempts; %d  swaps\n",
	       numNewPhantoms, numAttemptedPhantoms, totalSwaps);
    logMessage(MSG_MANAGER, "%d phantoms now real out of %d attempts\n",
	       numNewReals, numAttemptedReals);

    delete[] bdryFaces;
    assert(m_compMesh->isValid());

    // Now swap all the cells near the shock-relevant points again.
    improveMesh();

    assert(m_compMesh->isValid());
    logMessage(MSG_MANAGER, "Done updating mesh\n");
  }

  bool
  InterfaceTracker::isCloseToShock(Vert* const pV, const ADT& BBoxTree,
				   BFace* bdryFaces[], const double ratio)
  {
//	m_compMesh->setLengthScale(pV);
    if (pV->getVertType() == Vert::ePhantom) {
      m_length->setLengthScale(pV);
    }
    double dist = ratio * pV->getLengthScale();
    double range[] =
      { pV->x() - dist, pV->x() + dist, pV->y() - dist, pV->y() + dist, pV->z()
	  - dist, pV->z() + dist };
    std::vector<GR_index_t> queryResult = BBoxTree.veciRangeQuery(range);

    // Now we have a list of all the shock faces whose bounding boxes overlap
    // with a bounding box of the specified distance around the current vertex.

    for (size_t ii = 0; ii < queryResult.size(); ii++) {
      BFace *pBF = bdryFaces[queryResult[ii]];
      double projPt[3];
      pBF->findClosestPointOnBFace(pV->getCoords(), projPt);
      double thisDist = dDIST3D(pV->getCoords(), projPt);
      if (thisDist < dist) {
//			logMessage(1, "Vert: (%15f %15f %15f) Dist: %10f  Target:  %10f\n",
//						pV->x(), pV->y(), pV->z(), thisDist, dist);
//			logMessage(1, "   Near: %p (%15f %15f %15f) \n", pBF, projPt[0], projPt[1], projPt[2]);
	return true;
      }
    }

    return false;
  }

  bool
  InterfaceTracker::makeVertexPhantom(Vert* const pV, int& numSwaps)
  {
    assert(pV->getVertType() == Vert::eInterior);
    bool success = m_compMesh->removeVert(pV, numSwaps);
    return success;
  }

  bool
  InterfaceTracker::makePhantomReal(Vert* const pV)
  {
    Vert *pVRet = m_pSI3D->insertPoint(pV->getCoords(), m_compMesh->getCell(0),
				       pV);
    return (pVRet->isValid());
  }

}

