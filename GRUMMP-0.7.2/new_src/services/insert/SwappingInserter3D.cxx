#include "GR_assert.h"
#include <math.h>
#include <stdlib.h>
#include "GR_misc.h"

#include "GR_BFace.h"
#include "GR_Geometry.h"
#include "GR_InsertionManager.h"
#include "GR_Vec.h"
#include "GR_VolMesh.h"

namespace GRUMMP
{

  static Vert*
  getNewVert(const double newPtLoc[], VolMesh* const pVM,
	     Vert* const existingVert)
  {
    Vert* newVert;
    if (existingVert == pVInvalidVert) {
      newVert = pVM->createVert(newPtLoc);
    }
    else {
      newVert = existingVert;
      newVert->setCoords(3, newPtLoc);
    }
    return newVert;
  }
/// Insert Point on subseg. Difference between insert point on edge is that
/// a sight check is needed.

  Vert*
  SwappingInserter3D::insertPointOnSubseg(const double newPtLoc[],
					  Vert* const pV0, Vert* const pV1,
					  Vert* const existingVert)
  {
    // Need to clear the queue in advance to avoid a whole bunch of
    // unnecessary swap checks when calling this occasionally during
    // predominantly Watson insertion.
    m_swapper->clearQueue();
    // This implementation works for all edges: interior, bdry, internal
    // bdry, non-manifold, etc.
    bool qIsIntBdryFace = false;
    // Get all faces common to pV0 and pV1.
    int maxNumFaces = pV0->getNumFaces();
    std::vector<Face*> vecFacesIncOnBoth;
    {
      for (int i = 0; i < maxNumFaces; i++) {
	Face* pF = pV0->getFace(i);
	assert(pF->isValid() && !pF->isDeleted() && pF->hasVert(pV0));
	if (pF->hasVert(pV1)) {
	  vecFacesIncOnBoth.push_back(pF);
	}
      }
    }

    // Identify all the incident cells, faces not incident on the edge,
    // and internal and regular bdry faces.
    int numIncOnBoth = vecFacesIncOnBoth.size();
    if (numIncOnBoth == 0)
      return pVInvalidVert;

    std::set<Cell*> setIncCells;	//, setInvalidCells;
    std::set<BFace*> setIncBFaces;
    std::map<Face*, int> mapStarFaces;
    bool isHullValid = true;
    for (int i = 0; i < numIncOnBoth; i++) {
      Face *pF = vecFacesIncOnBoth[i];
      for (int ii = 0; ii < 2; ii++) {
	Cell *pC = (ii == 0 ? pF->getLeftCell() : pF->getRightCell());
	switch (pC->getType())
	  {
	  case CellSkel::eTet:
	    {
	      for (int iF = 0; iF < 4; iF++) {
		Face *pFLoc = pC->getFace(iF);
		if (!pFLoc->hasVert(pV0) || !pFLoc->hasVert(pV1)) {

		  bool qLeft = (pFLoc->getLeftCell() == pC ? true : false);

		  int iOrient = checkOrient3D(pFLoc, newPtLoc);
		  if ((iOrient == 1 && !qLeft) || (iOrient == -1 && qLeft)) {
		    mapStarFaces.insert(std::make_pair(pFLoc, pC->getRegion()));
		    setIncCells.insert(pC);

		  }
		  else {
//							setInvalidCells.insert(pC);
//							printf("%p pC %d orient %p %p false\n",pC, checkOrient3D(pF,newPtLoc), pF,pFLoc);
		    isHullValid = false;
		  }

		}
	      }

	    }
	    break;
	  case CellSkel::eIntTriBFace:
	    qIsIntBdryFace = true;
	    // fallthrough
	    // no break
	  case CellSkel::eTriBFace:
	    setIncBFaces.insert(dynamic_cast<BFace*>(pC));
	    break;
	  default:
	    assert(0);
	    break;
	  }
      } // Loop over cells incident on the face.
    }
//	printf("qq\n");
    if (!isHullValid && qIsIntBdryFace) {
      return pVInvalidVert;
    }
//	//find the cell
//	Cell * testCell;
//	for (std::set<Cell*>::iterator iter = setIncCells.begin(); iter
//	!= setIncCells.end(); ++iter) {
//		Cell *pC = *iter;
//		if(isInTet(newPtLoc,pC))
//			testCell = pC;
//	}
//	assert(testCell->isValid());
//	insertPointInCell()
//}

    VolMesh *pVM = dynamic_cast<VolMesh*>(m_pMesh);
    assert(pVM);

    Vert* newVert = getNewVert(newPtLoc, pVM, existingVert);
    if (setIncBFaces.empty()) {
      m_nOnEdge++;
      newVert->setType(Vert::eInterior);
    }
    else {
      m_nOnBdryEdge++;
      if (qIsIntBdryFace)
	newVert->setType(Vert::eBdryTwoSide);
      else
	newVert->setType(Vert::eBdry);

    }

    // For every cell:  delete the cell.
    for (std::set<Cell*>::iterator iter = setIncCells.begin();
	iter != setIncCells.end(); ++iter) {
      Cell *pC = *iter;
      pVM->deleteCell(pC);
    }

    std::set<Face*> setFacesToSwap;
    // For every face:  create a tet with the face and pVNew.
    for (std::map<Face*, int>::const_iterator iter = mapStarFaces.begin();
	iter != mapStarFaces.end(); ++iter) {
      Face *pF = iter->first;
      int iReg = iter->second;
      bool qExist = false;
      Cell *pCNew = pVM->createTetCell(qExist, newVert, pF, iReg);
      for (int i = 0; i < 4; i++) {
	setFacesToSwap.insert(pCNew->getFace(i));
      }
      assert(!qExist);
    }

    // For every bdry face: find the two faces covering it, create new
    // BFaces for these and delete the old BFace.  There may not be any of
    // these, in which case the loop never executes.
    for (std::set<BFace*>::const_iterator iter = setIncBFaces.begin();
	iter != setIncBFaces.end(); ++iter) {
      BFace *pBF = *iter;
      assert(pBF->isValid());
      Face *pF = pBF->getFace();
      assert(pF->isValid()); // Shouldn't have been detached yet.
      assert(!pBF->isDeleted()); // Neither one should be deleted.
      assert(!pF->isDeleted());
      Vert *pVOther = pF->getVert(0);
      if (pVOther == pV0 || pVOther == pV1) {
	pVOther = pF->getVert(1);
	if (pVOther == pV0 || pVOther == pV1) {
	  pVOther = pF->getVert(2);
	}
      }
      assert(pVOther != pV0 && pVOther != pV1);
      // Grab one of the new faces.
      Face *pFNew0 = findCommonFace(newVert, pV0, pVOther, pVInvalidVert, true);
      assert(pFNew0->isValid());
      BFace *pBFNew0 = pVM->createBFace(pFNew0, pBF);
      assert(pBFNew0->doFullCheck());
      assert(pFNew0->doFullCheck());

      // Grab the other new face.
      Face *pFNew1 = findCommonFace(newVert, pV1, pVOther, pVInvalidVert, true);
      assert(pFNew1->isValid());
      BFace *pBFNew1 = pVM->createBFace(pFNew1, pBF);
      assert(pBFNew1->doFullCheck());
      assert(pFNew1->doFullCheck());

      // Delete the old BFace; the Face will auto-delete.
      pVM->deleteBFace(pBF);
    }

    // For every int bdry face: find the two faces covering it, create new
    // BFaces for these and delete the old BFace.  There may not be any of
    // these, in which case the loop never executes.

    // This loop turns out to be exactly the same as the previous one; the
    // only likely difference up front was in the creation call, but
    // internal BFaces use the same call as regular BFaces.
    //@@ Do face swapping for all the faces

    if (omp_in_parallel() == 0) {
      if (m_swapper)
	m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    else if (m_swapper)
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true, true);

    //   assert(qValid());
    return newVert;
  }
/// Insert a point on an edge
/// Preconditions:
///   1.  New tets to be created should be right-handed; this isn't checked
///       for early enough (only as an assertion at tet creation time).
///   2.  The two verts must share an edge.
  Vert*
  SwappingInserter3D::insertPointOnEdge(const double newPtLoc[],
					Vert* const pV0, Vert* const pV1,
					Vert* const existingVert)
  {
    // This implementation works for all edges: interior, bdry, internal
    // bdry, non-manifold, etc.
    bool qIsIntBdryFace = false;
    // Get all faces common to pV0 and pV1.
    int maxNumFaces = pV0->getNumFaces();
    std::vector<Face*> vecFacesIncOnBoth;
    {
      for (int i = 0; i < maxNumFaces; i++) {
	Face* pF = pV0->getFace(i);
	assert(pF->isValid() && !pF->isDeleted() && pF->hasVert(pV0));
	if (pF->hasVert(pV1)) {
	  vecFacesIncOnBoth.push_back(pF);
	}
      }
    }

    // Identify all the incident cells, faces not incident on the edge,
    // and internal and regular bdry faces.
    int numIncOnBoth = vecFacesIncOnBoth.size();
    if (numIncOnBoth == 0)
      return pVInvalidVert;

    std::set<Cell*> setIncCells;
    std::set<BFace*> setIncBFaces;
    std::map<Face*, int> mapStarFaces;
    for (int i = 0; i < numIncOnBoth; i++) {
      Face *pF = vecFacesIncOnBoth[i];
      for (int ii = 0; ii < 2; ii++) {
	Cell *pC = (ii == 0 ? pF->getLeftCell() : pF->getRightCell());
	switch (pC->getType())
	  {
	  case CellSkel::eTet:
	    setIncCells.insert(pC);
	    for (int iF = 0; iF < 4; iF++) {
	      Face *pFLoc = pC->getFace(iF);
	      if (!pFLoc->hasVert(pV0) || !pFLoc->hasVert(pV1)) {
		mapStarFaces.insert(std::make_pair(pFLoc, pC->getRegion()));
	      }
	    }
	    break;
	  case CellSkel::eIntTriBFace:
	    qIsIntBdryFace = true;
	    // fallthrough
	    // no break
	  case CellSkel::eTriBFace:
	    setIncBFaces.insert(dynamic_cast<BFace*>(pC));
	    break;
	  default:
	    assert(0);
	    break;
	  }
      } // Loop over cells incident on the face.
    }

    VolMesh *pVM = dynamic_cast<VolMesh*>(m_pMesh);
    assert(pVM);
    Vert* newVert = getNewVert(newPtLoc, pVM, existingVert);
    if (setIncBFaces.empty()) {
      m_nOnEdge++;
      newVert->setType(Vert::eInterior);
    }
    else {
      m_nOnBdryEdge++;
      if (qIsIntBdryFace)
	newVert->setType(Vert::eBdryTwoSide);
      else
	newVert->setType(Vert::eBdry);

    }

    // For every cell:  delete the cell.
    for (std::set<Cell*>::iterator iter = setIncCells.begin();
	iter != setIncCells.end(); ++iter) {
      Cell *pC = *iter;
      pVM->deleteCell(pC);
    }

    std::set<Face*> setFacesToSwap;
    // For every face:  create a tet with the face and pVNew.
    for (std::map<Face*, int>::const_iterator iter = mapStarFaces.begin();
	iter != mapStarFaces.end(); ++iter) {
      Face *pF = iter->first;
      int iReg = iter->second;
      bool qExist = false;
      Cell *pCNew = pVM->createTetCell(qExist, newVert, pF, iReg);
      for (int i = 0; i < 4; i++) {
	setFacesToSwap.insert(pCNew->getFace(i));
      }
      assert(!qExist);
    }

    // For every bdry face: find the two faces covering it, create new
    // BFaces for these and delete the old BFace.  There may not be any of
    // these, in which case the loop never executes.
    for (std::set<BFace*>::const_iterator iter = setIncBFaces.begin();
	iter != setIncBFaces.end(); ++iter) {
      BFace *pBF = *iter;
      assert(pBF->isValid());
      Face *pF = pBF->getFace();
      assert(pF->isValid()); // Shouldn't have been detached yet.
      assert(!pBF->isDeleted()); // Neither one should be deleted.
      assert(!pF->isDeleted());
      Vert *pVOther = pF->getVert(0);
      if (pVOther == pV0 || pVOther == pV1) {
	pVOther = pF->getVert(1);
	if (pVOther == pV0 || pVOther == pV1) {
	  pVOther = pF->getVert(2);
	}
      }
      assert(pVOther != pV0 && pVOther != pV1);
      // Grab one of the new faces.
      Face *pFNew0 = findCommonFace(newVert, pV0, pVOther, pVInvalidVert, true);
      assert(pFNew0->isValid());
      BFace *pBFNew0 = pVM->createBFace(pFNew0, pBF);
      assert(pBFNew0->doFullCheck());
      assert(pFNew0->doFullCheck());

      // Grab the other new face.
      Face *pFNew1 = findCommonFace(newVert, pV1, pVOther, pVInvalidVert, true);
      assert(pFNew1->isValid());
      BFace *pBFNew1 = pVM->createBFace(pFNew1, pBF);
      assert(pBFNew1->doFullCheck());
      assert(pFNew1->doFullCheck());

#ifndef OMIT_VERTEX_HINTS
      // Make sure that bdry verts have bdry hints
      for (int ii = pFNew0->getNumVerts() - 1; ii >= 0; ii--) {
	pFNew0->getVert(ii)->setHintFace(pFNew0);
      }
#endif

      // Delete the old BFace; the Face will auto-delete.
      pVM->deleteBFace(pBF);
    }

    // For every int bdry face: find the two faces covering it, create new
    // BFaces for these and delete the old BFace.  There may not be any of
    // these, in which case the loop never executes.

    // This loop turns out to be exactly the same as the previous one; the
    // only likely difference up front was in the creation call, but
    // internal BFaces use the same call as regular BFaces.
    //@@ Do face swapping for all the faces
    if (omp_in_parallel() == 0) {
      if (m_swapper)
	m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    else if (m_swapper)
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true, true);

    //   assert(qValid());
    return newVert;
  }

//Vert* SwappingInserter3D::insertPointOnCurvedIntBFace(const double newPtLoc[],
//		Face* const face,
//	  Vert* const existingVert){
//	// we need to walk from the face to the point, and find all cells it crosses.
//	// Can assume we don't hit another BFace, because otherwise
//	// we would have encroached it
//	// going to assume this only happens for internal boundaries for now
//
//	// going to create a hull consisting of the outside intbface (that can't see
//	// the new point location, and should not be deleted (then delete the bface),
//	// deleting the inside face (that can see the point), and connect the new point
//	// to those three points from the original bface to form a tet, with three new
//	// bfaces created. Then create the hull as simply as possible, walking until we
//	// have all the cells in it
//
//	// first, figure out what side we want to be on
//	Cell *apC[] =
//	{face->getLeftCell(), face->getRightCell()};
//	if (apC[0]->getType() != CellSkel::eTet) {
//		if (apC[1]->getType() != CellSkel::eTet) {
//			// At this point, we're doomed: both cells are non-tets.
//			return pVInvalidVert;
//		}
//		Cell *pCTmp = apC[0];
//		apC[0] = apC[1];
//		apC[1] = pCTmp;
//		assert(apC[0]->getType() == CellSkel::eTet);
//	}
//	BFace * pBF = dynamic_cast<BFace*> (apC[1]);
//	assert(pBF->isValid());
//	// Get the tet across the internal bdry.
//	Face * pFOther = pBF->getFace(0);
//	if (face == pFOther) pFOther = pBF->getFace(1);
//	assert(face != pFOther);
//	apC[1] = pFOther->getOppositeCell(pBF); // This tet's in the other region.
//	assert(apC[0]->getType() == CellSkel::eTet);
//	assert(apC[1]->getType() == CellSkel::eTet);
//	assert(apC[0] != apC[1]);
//
//	// now we have two tets, figure out which one we are on the same side of.
//	Vert * v0 = apC[0]->getOppositeVert(face);
//	int iOrientFace = checkOrient3D(face,newPtLoc);
//	assert(checkOrient3D(pFOther,newPtLoc) != 0);
//
//	int iOrient0 = checkOrient3D(face,v0);
//	assert(iOrientFace != 0); // if it were coplanar, we wouldn't need to be here
//	Face * pFInside;
//	Face * pFOutside;
//	if(iOrientFace == iOrient0){
//		// on same side as v0
//		pFInside = face;
//		pFOutside = pFOther;
//	} else {
//		// on other side
//		pFInside = pFOther;
//		pFOutside = face;
//
//		// swap to ensure that the cell we care about is apC[0]
//		Cell *pCTmp = apC[0];
//		apC[0] = apC[1];
//		apC[1] = pCTmp;
//	}
//	// now we know which faces.. now to start our walk from apC[0]
//	std::queue<Cell*> cellsToTest;
//	std::set<const Cell*> cellsTested;
//	cellsToTest.push(apC[0]);
//	cellsToTest.push(apC[1]);
//	CavityInserter3D CI(m_pMesh);
//
//	CI.addCellToCavity(apC[0]);
//
//	// go hunting on ONE side for cells that
//	// have a face that crosses the line of sight between the new point
//	// and one of the boundaries
//	while (!cellsToTest.empty()) {
//
//		Cell* cell = cellsToTest.front();
//		cellsToTest.pop();
//		// Don't bother if this cell has already been tested.
//		if (!cellsTested.insert(cell).second) continue;
//
//		SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);
//		if (simplexCell->isValid()) {
//			// do a sight check between a vert of bdry face-> new point loc
//			// and faces of test cells that do not touch the bdry face
//			// this SHOULD flag cells that can see it
//
//			// make sure we get all tets that contain the point in anyway.
//			if(isInTet(newPtLoc,simplexCell) > 0)
//				CI.addCellToCavity(simplexCell);
//
//
//			for (int i = 0; i < simplexCell->getNumFaces(); i++) {
//				Face * tmpFace = simplexCell->getFace(i);
//				for(int iV = 0; iV < 3; iV++){
//					// only check faces that don't contain the bface we are looking at
//					if(! tmpFace->hasVert(pBF->getVert(iV))){
//						if( calcSegmentTriangleIntersection(existingVert->getCoords(),pBF->getVert(iV)->getCoords(), tmpFace->getVert(0),
//								tmpFace->getVert(1), tmpFace->getVert(2), NULL) == 1){
//							CI.addCellToCavity(cell);
//						}
//					}
//				}
//				cellsToTest.push( tmpFace->getOppositeCell(simplexCell));
//
//			}
//		}
//	}
//	printf("Inserting on Curved BFace\n");
//	writeVTK_SingleFace(face,"tempface.vtk");
//CI.writeCavityVTK("cavity.vtk");
//if(!CI.isHullValid(newPtLoc))
//writeVTK_Legacy(*dynamic_cast<VolMesh*>(m_pMesh),"tempmesh");
////	WatsonInserter WI(dynamic_cast<VolMesh*>(m_pMesh),m_pMesh->getEncroachmentType());
////	WI.computeHull(newPtLoc,NULL);
////	WI.writeCavityVTK("Cavity.vtk");
//	Vert * newVert = CI.insertPointInHull(existingVert);
//
//	// at this point, it appears as an internal insertion. This is much easier to do.
//	// now to 'pop' the boundary face up to the new vertex.
//
//	// now we have the old intbdry face, but need to move it up to the three with the newVert
//	// easiest way is through pFInside
//
//	Cell * newCell = pFInside->getOppositeCell(pBF);
//	for(int iF = 0; iF < 4; iF++)
//		if(newCell->getFace(iF)!= pFInside)
//			m_pMesh->createBFace(newCell->getFace(iF),pBF);
//
//	//get rid of one of the duplicate faces present for the old BFace;
//	apC[1]->replaceFace(pFOutside,pFInside);
//	pFOutside->removeCell(apC[1]);
//	m_pMesh->deleteBFace(pBF);
//	pFInside->addCell(apC[1]);
//
//	// new cell is on other side, need to reset region
//	newCell->setRegion(apC[1]->getRegion());
//	newVert->setType(Vert::eBdryTwoSide);
//
//	if (omp_in_parallel() == 0){
//		if (m_swapper) m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();}
//	else
//		if (m_swapper) m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true,true);
//
//	return newVert;
//}
//Vert* SwappingInserter3D::insertPointOnCurvedIntBFace(const double newPtLoc[],
//		Face* const face,
//	  Vert* const existingVert){
//	// we need to walk from the face to the point, and find all cells it crosses.
//	// Can assume we don't hit another BFace, because otherwise
//	// we would have encroached it
//	// going to assume this only happens for internal boundaries for now
//
//	// going to create a hull consisting of the outside intbface (that can't see
//	// the new point location, and should not be deleted (then delete the bface),
//	// deleting the inside face (that can see the point), and connect the new point
//	// to those three points from the original bface to form a tet, with three new
//	// bfaces created. Then create the hull as simply as possible, walking until we
//	// have all the cells in it
//
//	// first, figure out what side we want to be on
//	Cell *apC[] =
//	{face->getLeftCell(), face->getRightCell()};
//	if (apC[0]->getType() != CellSkel::eTet) {
//		if (apC[1]->getType() != CellSkel::eTet) {
//			// At this point, we're doomed: both cells are non-tets.
//			return pVInvalidVert;
//		}
//		Cell *pCTmp = apC[0];
//		apC[0] = apC[1];
//		apC[1] = pCTmp;
//		assert(apC[0]->getType() == CellSkel::eTet);
//	}
//	BFace * pBF = dynamic_cast<BFace*> (apC[1]);
//	assert(pBF->isValid());
//	// Get the tet across the internal bdry.
//	Face * pFOther = pBF->getFace(0);
//	if (face == pFOther) pFOther = pBF->getFace(1);
//	assert(face != pFOther);
//	apC[1] = pFOther->getOppositeCell(pBF); // This tet's in the other region.
//	assert(apC[0]->getType() == CellSkel::eTet);
//	assert(apC[1]->getType() == CellSkel::eTet);
//	assert(apC[0] != apC[1]);
//
//	// now we have two tets, figure out which one we are on the same side of.
//	Vert * v0 = apC[0]->getOppositeVert(face);
//	int iOrientFace = checkOrient3D(face,newPtLoc);
//	assert(checkOrient3D(pFOther,newPtLoc) != 0);
//
//	int iOrient0 = checkOrient3D(face,v0);
//	assert(iOrientFace != 0); // if it were coplanar, we wouldn't need to be here
//	Face * pFInside;
//	Face * pFOutside;
//	if(iOrientFace == iOrient0){
//		// on same side as v0
//		pFInside = face;
//		pFOutside = pFOther;
//	} else {
//		// on other side
//		pFInside = pFOther;
//		pFOutside = face;
//
//		// swap to ensure that the cell we care about is apC[0]
//		Cell *pCTmp = apC[0];
//		apC[0] = apC[1];
//		apC[1] = pCTmp;
//	}
//	// now we know which faces.. now to start our walk from apC[0]
//	std::queue<Cell*> cellsToTest;
//	std::set<const Cell*> cellsTested;
//	cellsToTest.push(apC[0]);
//	cellsToTest.push(apC[1]);
//	CavityInserter3D CI(m_pMesh);
//
//	CI.addCellToCavity(apC[0]);
//
//	// go hunting on ONE side for cells that
//	// have a face that crosses the line of sight between the new point
//	// and one of the boundaries
//	while (!cellsToTest.empty()) {
//
//		Cell* cell = cellsToTest.front();
//		cellsToTest.pop();
//		// Don't bother if this cell has already been tested.
//		if (!cellsTested.insert(cell).second) continue;
//
//		SimplexCell* simplexCell = dynamic_cast<SimplexCell*>(cell);
//		if (simplexCell->isValid()) {
//			// do a sight check between a vert of bdry face-> new point loc
//			// and faces of test cells that do not touch the bdry face
//			// this SHOULD flag cells that can see it
//
//			// make sure we get all tets that contain the point in anyway.
//			if(isInTet(newPtLoc,simplexCell) > 0)
//				CI.addCellToCavity(simplexCell);
//
//
//			for (int i = 0; i < simplexCell->getNumFaces(); i++) {
//				Face * tmpFace = simplexCell->getFace(i);
//				for(int iV = 0; iV < 3; iV++){
//					// only check faces that don't contain the bface we are looking at
//					if(! tmpFace->hasVert(pBF->getVert(iV))){
//						if( calcSegmentTriangleIntersection(existingVert->getCoords(),pBF->getVert(iV)->getCoords(), tmpFace->getVert(0),
//								tmpFace->getVert(1), tmpFace->getVert(2),0., NULL) == 1){
//							CI.addCellToCavity(cell);
//						}
//					}
//				}
//				cellsToTest.push( tmpFace->getOppositeCell(simplexCell));
//
//			}
//		}
//	}
//	printf("Inserting on Curved BFace\n");
//	writeVTK_SingleFace(face,"tempface.vtk");
//CI.writeCavityVTK("cavity.vtk");
//if(!CI.isHullValid(newPtLoc))
//writeVTK_Legacy(*dynamic_cast<VolMesh*>(m_pMesh),"tempmesh");
////	WatsonInserter WI(dynamic_cast<VolMesh*>(m_pMesh),m_pMesh->getEncroachmentType());
////	WI.computeHull(newPtLoc,NULL);
////	WI.writeCavityVTK("Cavity.vtk");
//	Vert * newVert = CI.insertPointInHull(existingVert);
//
//	// at this point, it appears as an internal insertion. This is much easier to do.
//	// now to 'pop' the boundary face up to the new vertex.
//
//	// now we have the old intbdry face, but need to move it up to the three with the newVert
//	// easiest way is through pFInside
//
//	Cell * newCell = pFInside->getOppositeCell(pBF);
//	for(int iF = 0; iF < 4; iF++)
//		if(newCell->getFace(iF)!= pFInside)
//			m_pMesh->createBFace(newCell->getFace(iF),pBF);
//
//	//get rid of one of the duplicate faces present for the old BFace;
//	apC[1]->replaceFace(pFOutside,pFInside);
//	pFOutside->removeCell(apC[1]);
//	m_pMesh->deleteBFace(pBF);
//	pFInside->addCell(apC[1]);
//
//	// new cell is on other side, need to reset region
//	newCell->setRegion(apC[1]->getRegion());
//	newVert->setType(Vert::eBdryTwoSide);
//
//		if (omp_in_parallel() == 0){
//			if (m_swapper) m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();}
//		else
//			if (m_swapper) m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true,true);
//		}
//
//	return newVert;
//}

/// Insert a point on a tri face.
/// Preconditions:
///   1. The face must be a tri.
///   2. The new point must actually fall on the face (as measured by orient
///      predicates); this applies only for interior faces, as bdry faces may
///      be curved, in which case the new vert isn't on the triangle.
///   3. The new tets to be created must be right-handed (this can conceivably
///      fail for very flat tets); checked automatically during tet creation.
///   4. The cells on either side of the face must be tets.
  Vert*
  SwappingInserter3D::insertPointOnFace(const double newPtLoc[],
					Face* const face,
					Vert* const existingVert)
  {
    assert(face->doFullCheck());
    assert(face->getType() == Face::eTriFace);
    // This implementation works for all faces: interior, bdry, internal
    // bdry.

    // Get cells and bdry faces incident on this face.  For internal bdry
    // faces, must also get the tet on the other side of the internal
    // bdry.
    Cell *apC[] =
      { face->getLeftCell(), face->getRightCell() };
    if (apC[0]->getType() != CellSkel::eTet) {
      if (apC[1]->getType() != CellSkel::eTet) {
	// At this point, we're doomed: both cells are non-tets.
	return pVInvalidVert;
      }
      Cell *pCTmp = apC[0];
      apC[0] = apC[1];
      apC[1] = pCTmp;
      assert(apC[0]->getType() == CellSkel::eTet);
    }

    VolMesh *pVM = dynamic_cast<VolMesh*>(m_pMesh);
    assert(pVM);
    Vert *newVert = NULL;
    // used for awkward case of point not on face, and not in cell.
    BFace *pBF = pBFInvalidBFace;
    Face *pFOther = pFInvalidFace;
    if (apC[1]->getType() == CellSkel::eTriBFace) {
      pBF = dynamic_cast<BFace*>(apC[1]);
      assert(pBF->isValid());
      apC[1] = pCInvalidCell; // No tet on the other side to worry about.

      // checks if surface is curved first
      if (pBF->getTopologicalParent()) {
	bool qInvalidCell = false;
	// special set of checks for awkward case

	Cell *pC = apC[0];
	assert(apC[0]->getType() == CellSkel::eTet);

	if (pC->isValid()) {
	  for (int iF = 0; iF < 4; iF++) {
	    Face *pFLoc = pC->getFace(iF);
	    if (pFLoc != face) {
	      int iOrient = checkOrient3D(pFLoc, newPtLoc);
	      if (iOrient != -1)
		qInvalidCell = true;

	    }
	  }

	}

	// if the point can't see either adjacent cell (likely outside of both)

	if (qInvalidCell) {
	  return pVInvalidVert;
	}
      }

      m_nOnBdryTri++;
      newVert = getNewVert(newPtLoc, pVM, existingVert);
      newVert->setType(Vert::eBdry);
    }
    else if (apC[1]->getType() == CellSkel::eIntTriBFace) {
      pBF = dynamic_cast<BFace*>(apC[1]);
      assert(pBF->isValid());
      // Get the tet across the internal bdry.
      pFOther = pBF->getFace(0);
      if (face == pFOther)
	pFOther = pBF->getFace(1);
      assert(face != pFOther);
      apC[1] = pFOther->getOppositeCell(pBF); // This tet's in the other region.
      assert(apC[0]->getType() == CellSkel::eTet);
      assert(apC[1]->getType() == CellSkel::eTet);
      assert(apC[0] != apC[1]);
      // checks if surface is curved first
      if (pBF->getTopologicalParent()) {
	bool qInvalidCell[2] =
	  { false, false };
	// special set of checks for awkward case
	for (int i = 0; i < 2; i++) {
	  Cell *pC = apC[i];
	  if (pC->isValid()) {
	    for (int iF = 0; iF < 4; iF++) {
	      Face *pFLoc = pC->getFace(iF);
	      if (pFLoc != face && pFLoc != pFOther) {
		bool qLeft = (pFLoc->getLeftCell() == pC ? true : false);
		int iOrient = checkOrient3D(pFLoc, newPtLoc);
		if (!(iOrient == 1 && !qLeft) || (iOrient == -1 && qLeft))
		  qInvalidCell[i] = true;

	      }
	    }
	  }
	}

	// if the point can't see either adjacent cell (likely outside of both)
	if (qInvalidCell[0] || qInvalidCell[1]) {
	  return pVInvalidVert;
	}
      }
      m_nOnIntBdryTri++;
      newVert = getNewVert(newPtLoc, pVM, existingVert);
      newVert->setType(Vert::eBdryTwoSide);
    }
    else {
      newVert = getNewVert(newPtLoc, pVM, existingVert);
      // Confirm that the vertex is in fact on the face of the tet(s) for
      // interior face insertions.
      assert(!pBF && isInTet(newPtLoc, apC[0]) == 3);
      assert(!pBF && (!apC[1]->isValid() || isInTet(newPtLoc, apC[1]) == 3));

    }

    std::map<Face*, int> mapHullFaces;
    for (int i = 0; i < 2; i++) {
      Cell *pC = apC[i];
      if (pC->isValid()) {
	for (int iF = 0; iF < 4; iF++) {
	  Face *pFLoc = pC->getFace(iF);
	  if (pFLoc != face && pFLoc != pFOther) {

	    mapHullFaces.insert(std::make_pair(pFLoc, pC->getRegion()));
	  }
	}
	pVM->deleteCell(pC);
      }
    }
    // Create new tets (one for each hull face).
    for (std::map<Face*, int>::iterator iter = mapHullFaces.begin();
	iter != mapHullFaces.end(); ++iter) {
      Face *pFCand = iter->first;
      int iReg = iter->second;
      bool qExist = false;
      pVM->createTetCell(qExist, newVert, pFCand, iReg);
      assert(qExist == false);
    }

    // Create new bdry faces, if any.
    if (pBF) {
      for (int iV = 0; iV < 3; iV++) {
	Vert* pVA = pBF->getVert(iV);
	Vert* pVB = pBF->getVert((iV + 1) % 3);
	Face *pFNew = findCommonFace(pVA, pVB, newVert, pVInvalidVert, true);
	BFace *pBFNew = pVM->createBFace(pFNew, pBF);
	assert(pBFNew->doFullCheck());
	assert(pFNew->doFullCheck());
#ifndef OMIT_VERTEX_HINTS
	// Make sure that bdry verts have bdry hints
	for (int ii = pFNew->getNumVerts() - 1; ii >= 0; ii--) {
	  pFNew->getVert(ii)->setHintFace(pFNew);
	}
#endif

      }
      pVM->deleteBFace(pBF); // Old face(s) attached to it auto-delete.
    }
    else {
      m_nOnTri++;
      newVert->setType(Vert::eInterior);
    }

    if (omp_in_parallel() == 0) {
      if (m_swapper)
	m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    else if (m_swapper)
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true, true);

    return newVert;
  }

/// Insert a point inside a tet cell.
/// Preconditions:
///   1. The cell must be a tet.
///   2. The new point must actually fall inside the tet.
///   3. The new tets to be created must be right-handed (this can conceivably
///      fail for very flat tets).
  Vert*
  SwappingInserter3D::insertPointInCell(const double newPtLoc[],
					Cell * const cell,
					Vert * const existingVert)
  {
    assert(cell->doFullCheck());
    assert(cell->getType() == Cell::eTet);

    VolMesh *pVM = dynamic_cast<VolMesh*>(m_pMesh);
    assert(pVM);

    // The caller is responsible for being sure the point is really inside
    // the tet.  This assertion tests that, as well as testing the orientation
    // of created tets.
    assert(isInTet(newPtLoc, cell) == 4);

    m_nInCell++;

    // The order of the faces is unimportant; the create calls automatically
    // connect things properly anyway.
    Face* faces[4];
    for (int i = 0; i < 4; i++)
      faces[i] = cell->getFace(i);

    // With the new topology handling, this is trivial: delete one tet,
    // create four, and move on with your life.

    // Need to know which region first...
    int iReg = cell->getRegion();
    pVM->deleteCell(cell);

    Vert* newVert = getNewVert(newPtLoc, pVM, existingVert);
    newVert->setType(Vert::eInterior);
    for (int i = 0; i < 4; i++) {
      bool existed = true;
      Cell *newTet = pVM->createTetCell(existed, newVert, faces[i], iReg);
      assert(!existed);
      assert(newTet->doFullCheck());
    }
    if (omp_in_parallel() == 0) {
      if (m_swapper)
	m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    else if (m_swapper)
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive(true, true);
    return newVert;
  }

//@ Add a site (or some variant if it lies outside the mesh) to the mesh
// Return value indicates success
  Vert*
  SwappingInserter3D::insertPoint(const double newPtLoc_in[3],
				  Cell* const cellGuess,
				  Vert* const existingVert)
// pVNew is non-pVInvalidVert if and only if the point being
// inserted already belongs to the point set of the mesh; in this
// case, pVNew->adCoords() == adPointIn.
  {
    bool status = false;
    Cell* cell = dynamic_cast<VolMesh*>(m_pMesh)->findCell(newPtLoc_in,
							   cellGuess, status);
    double newPtLoc[] =
      { newPtLoc_in[0], newPtLoc_in[1], newPtLoc_in[2] };
    if (!status && cell->getType() != Cell::eTriBFace) {
      return pVInvalidVert;
    }

    //@@ Compute barycentric coordinates
    double adBary[4];
    TetCell* tetCell = dynamic_cast<TetCell*>(cell);
    assert(tetCell);
    tetCell->calcBarycentricCoords(newPtLoc, adBary);
    //      int iNegBary = ((adBary[0] < -1.e-10 ? 1 : 0) + (adBary[1] < -1.e-10 ? 1
    //          : 0) + (adBary[2] < -1.e-10 ? 1 : 0) + (adBary[3] < -1.e-10 ? 1 : 0));

    //@@ Snap to face/edge

    // Define the threshold of barycentric coordinate that will induce
    // snap to face/edge
    double dBaryThresh = 0.1;
    if (forceInsertion())
      dBaryThresh = 1.e-10;

    double dSum = 0.;
    int i, iLowBary = 0;
    for (i = 0; i <= 3; i++)
      if (adBary[i] < dBaryThresh) {
	adBary[i] = 0;
	iLowBary++;
      }
      else
	dSum += adBary[i];
    assert(iFuzzyComp(dSum, 0) != 0);
    assert(iLowBary <= 3);
    switch (iLowBary)
      {
      case 0:
	{ // Insertion squarely inside a cell
	  return insertPointInCell(newPtLoc, cell, existingVert);
	}
      case 1:
	{ // Insertion on a face (possibly a boundary face; let the
	  // callee deal with that question)
	  // Find the face opposite the zero barycentric coordinate
	  const Vert *pVZero;
	  if (adBary[0] < dBaryThresh)
	    pVZero = cell->getVert(0);
	  else if (adBary[1] < dBaryThresh)
	    pVZero = cell->getVert(1);
	  else if (adBary[2] < dBaryThresh)
	    pVZero = cell->getVert(2);
	  else
	    pVZero = cell->getVert(3); // (adBary[3] == 0)
	  Face *pF = cell->getOppositeFace(pVZero);
	  // Project the proposed point onto this face; this will give the
	  // face circumcenter
	  assert(pF->getType() == Face::eTriFace);
	  if (!forceInsertion())
	    dynamic_cast<TriFace*>(pF)->calcCircumcenter(newPtLoc);

	  // Check barycentrics again.
	  dynamic_cast<TetCell*>(cell)->calcBarycentricCoords(newPtLoc, adBary);
	  double dFaceBaryThresh = 0.2;
	  if (forceInsertion())
	    dFaceBaryThresh = 1.e-10;

	  iLowBary = ((adBary[0] < dFaceBaryThresh ? 1 : 0)
	      + (adBary[1] < dFaceBaryThresh ? 1 : 0)
	      + (adBary[2] < dFaceBaryThresh ? 1 : 0)
	      + (adBary[3] < dFaceBaryThresh ? 1 : 0));
	  if (iLowBary == 1) { // In the interior of the face
	    return insertPointOnFace(newPtLoc, pF, existingVert);
	  }
	  for (int iBary = 0; iBary <= 3; iBary++)
	    if (adBary[iBary] < dFaceBaryThresh)
	      adBary[iBary] = 0;

	  assert(iLowBary == 2);
	}
	// fall through
	// no break if the face case becomes an edge case
      case 2:
	{
	  // Must be insertion on an edge (possibly a boundary edge)
	  // As with face insertion, let the callee handle the boundary case

	  // Find the two verts w/ non-zero barycentrics; the new vert falls
	  // on the edge between them
	  Vert *pV0, *pV1;
	  Face *pF;

	  // Not a very elegant way to identify these things, but so what?
	  if (adBary[0] < dBaryThresh) {
	    if (adBary[1] < dBaryThresh) {
	      pV0 = cell->getVert(2);
	      pV1 = cell->getVert(3);
	    }
	    else if (adBary[2] < dBaryThresh) {
	      pV0 = cell->getVert(1);
	      pV1 = cell->getVert(3);
	    }
	    else {
	      pV0 = cell->getVert(2);
	      pV1 = cell->getVert(1);
	    }
	    pF = cell->getOppositeFace(cell->getVert(0));
	  }
	  else if (adBary[1] < dBaryThresh) {
	    if (adBary[2] < dBaryThresh) {
	      pV0 = cell->getVert(0);
	      pV1 = cell->getVert(3);
	    }
	    else {
	      pV0 = cell->getVert(2);
	      pV1 = cell->getVert(0);
	    }
	    pF = cell->getOppositeFace(cell->getVert(1));
	  }
	  else {
	    pV0 = cell->getVert(0);
	    pV1 = cell->getVert(1);
	    pF = cell->getOppositeFace(cell->getVert(2));
	  }
	  assert(pV0->isValid());
	  assert(pV1->isValid());
	  assert(pF->isValid());

	  double adActual[3];
	  if (forceInsertion()) {
	    adActual[0] = newPtLoc[0];
	    adActual[1] = newPtLoc[1];
	    adActual[2] = newPtLoc[2];
	  }
	  else {
	    adActual[0] = 0.5 * (pV0->x() + pV1->x());
	    adActual[1] = 0.5 * (pV0->y() + pV1->y());
	    adActual[2] = 0.5 * (pV0->z() + pV1->z());
	  }
	  return insertPointOnEdge(adActual, pV0, pV1, existingVert);
	} // End edge insertion case
      case 3:
	{
	  // Inserting a point that's already in the mesh; this can occur when
	  // re-inserting interior structured verts if one wasn't deleted in
	  // the first place.  Return the vert in the mesh already inserted
	  for (int q = 0; q < 4; q++)
	    if (adBary[q] > (1.0 - dBaryThresh)) {
	      return tetCell->getVert(q);
	    }

	  return pVInvalidVert;
	}
      default:
	assert2(0, "Invalid number of small barycentrics");
	return pVInvalidVert;
      } // End switch on the number of small barycentrics
  } // End of point insertion
}
