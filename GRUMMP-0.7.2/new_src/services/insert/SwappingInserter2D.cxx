#include "CubitVector.hpp"

#include "GR_assert.h"
#include "GR_Geometry.h"
#include "GR_Mesh2D.h"
#include "GR_Vec.h"
#include "GR_BFace.h"
#include "GR_Aniso.h"
#include "GR_CellQueue.h"
#include "GR_GRCurve.h"
#include "GR_InsertionManager.h"
#include "GR_QualMeasure.h"

namespace GRUMMP
{
  static Vert*
  getNewVert(const double newPtLoc[], Mesh2D* const pM2D,
	     Vert* const existingVert)
  {
    Vert* newVert;
    if (existingVert == pVInvalidVert) {
      newVert = pM2D->createVert(newPtLoc);
    }
    else {
      newVert = existingVert;
      newVert->setCoords(2, newPtLoc);
    }
    return newVert;
  }

  Vert*
  SwappingInserter2D::insertPoint(const double newPtLoc[],
				  Cell * const cellGuess,
				  Vert* const existingVert)
  {
    // Logging state on exit should be INSERT_VERTEX if it starts out that way.
    //@@ Find which cell the site falls within (or near, if outside the mesh)
    assert(!cellGuess->isValid() || cellGuess->doFullCheck());
    int numViolated = 0, numTies = 0;
    Face *facesViolated[3], *facesTied[3];
    Cell* cell = dynamic_cast<Mesh2D*>(m_pMesh)->findCell(newPtLoc, cellGuess,
							  numViolated,
							  facesViolated,
							  numTies, facesTied);
    //@@ Decide which insertion case and insert
    switch (numViolated)
      {
      case 0:
	// In this case, the point definitely lies inside the mesh.
	{
	  switch (numTies)
	    {
	    case 0:
	      return insertPointInCell(newPtLoc, cell, existingVert);
	      break;
	    case 1:
	      // Grab the correct face for this cell; only one tie.
	      return insertPointOnFace(newPtLoc, facesTied[0], existingVert);
	      break;
	    case 2:
	      // Insertion attempt where there's already a vertex. Return Existing Vertex
	      for (GR_index_t iV = 0; iV < 3; iV++)
		if (dDIST2D(cell->getVert(iV)->getCoords(),newPtLoc) < 1e-10)
		  return cell->getVert(iV);

	      return pVInvalidVert;
	      break;
	    default:
	      // Other cases are impossible.
	      assert2(0, "Unreachable; should have returned");
	      return pVInvalidVert;
	      break;
	    }
	}
	break;
      case 1:
	//@@@ Outside:  Insert at edge midside
	// Only one face has the wrong orientation wrt the proposed vert
	{
	  Face *face = facesViolated[0];
	  if (m_pMesh->areBdryChangesAllowed()) {
	    if (face->getFaceLoc() == Face::eBdryFace) {
	      return insertPointOnBdryFace(newPtLoc, face, existingVert);
	    }
	    else {
	      assert(face->getFaceLoc() == Face::eBdryTwoSide);
	      return insertPointOnInternalBdryFace(newPtLoc, face, existingVert);
	    }
	  } // Done with allowed bdry changes
	}
	break;
      case 2:
	for (GR_index_t iV = 0; iV < 3; iV++)
	  if (dDIST2D(cell->getVert(iV)->getCoords(),newPtLoc) < 1e-10)
	    return cell->getVert(iV);
	//@@@ Outside:  Logical insertion site is at existing vert; abort
	return pVInvalidVert;
	break;
      default:
	assert2(0, "Unreachable case");
	return pVInvalidVert;
	break;
      }
    return pVInvalidVert;
  }

//@ Insert a point on a boundary face
  Vert*
  SwappingInserter2D::insertPointOnBdryFace(const double newPtLoc_arg[],
					    Face* const face,
					    Vert* const existingVert)
  {
    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(m_pMesh);
    assert(pM2D);

    BdryEdgeBase* bdryFace = dynamic_cast<BdryEdgeBase*>(face->getLeftCell());
    if (!bdryFace)
      bdryFace = dynamic_cast<BdryEdgeBase*>(face->getRightCell());
    assert(bdryFace->isValid());

    if (bdryFace->getType() == Cell::eIntBdryEdge)
      return insertPointOnInternalBdryFace(newPtLoc_arg, face, existingVert);

    TriCell *triCell = dynamic_cast<TriCell *>(face->getOppositeCell(bdryFace));
    if (!triCell->isValid())
      return pVInvalidVert; // Interior cell was a quad.

    m_nOnBdryEdge++;
    assert(face->doFullCheck());
    assert(triCell->doFullCheck());
    assert(face->getLeftCell() == triCell || face->getRightCell() == triCell);
    assert(triCell->getType() == Cell::eTriCell);
    assert(bdryFace->getType() == Cell::eBdryEdge);

    double newPtLoc[3];
    if (forceInsertion()) {
      newPtLoc[0] = newPtLoc_arg[0];
      newPtLoc[1] = newPtLoc_arg[1];
    }
    else {
      CubitVector newLoc;
      bool isConcentric;
      bdryFace->getSplitData(m_splitType, newLoc, isConcentric);
      newPtLoc[0] = newLoc[0];
      newPtLoc[1] = newLoc[1];
    }
    Vert *newVert = getNewVert(newPtLoc, pM2D, existingVert);
    assert(newVert->isValid());

    if (pM2D->pCellQueue != NULL) {
      pM2D->pCellQueue->vRemoveCell(triCell);
    }

    Vert *vertA, *vertB, *vertC;
    if (face->getLeftCell() == triCell) {
      vertA = face->getVert(0);
      vertB = face->getVert(1);
    }
    else {
      vertB = face->getVert(0);
      vertA = face->getVert(1);
    }
    vertC = triCell->getOppositeVert(face);
    Face *faceA = triCell->getOppositeFace(vertA);
    Face *faceB = triCell->getOppositeFace(vertB);
    if ((checkOrient2D(vertA, newVert, vertC) != 1)
	|| (checkOrient2D(vertC, newVert, vertB) != 1)) {
      assert2(false, "Problem with orientation when inserting bdry vertex.\n");
    }
    // TODO: The aniso insertion queue needs to be an observer.
    //   if(pCellQueue!=NULL) pCellQueue->vRemoveCell(cell);
    int iReg = triCell->getRegion();
    pM2D->deleteCell(triCell);
    newVert->setType(Vert::eBdry);

    bool qExist;
    Face *faceC1 = pM2D->createFace(qExist, vertB, newVert);
    assert(!qExist);
    Face *faceC2 = pM2D->createFace(qExist, vertA, newVert);
    assert(!qExist);
    Face *faceD = pM2D->createFace(qExist, vertC, newVert);
    assert(!qExist);

    Cell *cell1 = pM2D->createTriCell(faceA, faceC1, faceD, iReg);
    Cell *cell2 = pM2D->createTriCell(faceB, faceC2, faceD, iReg);

    pM2D->createBFace(faceC1, bdryFace);
    pM2D->createBFace(faceC2, bdryFace);

    pM2D->deleteBFace(bdryFace);
    assert(cell1->doFullCheck());
    assert(cell2->doFullCheck());
    assert(faceA->doFullCheck());
    assert(faceB->doFullCheck());
    assert(faceC1->doFullCheck());
    assert(faceC2->doFullCheck());
    assert(faceD->doFullCheck());

    if (pM2D->pAnisoRef != NULL) {
      pM2D->pAnisoRef->vInterpolateMetric(newVert);
//		m_length->setLengthScale(newVert);
    }

    if (pM2D->pCellQueue != NULL) {
      pM2D->pCellQueue->vAddCell(cell1);
      pM2D->pCellQueue->vAddCell(cell2);
    }

    if (pM2D->pAnisoRef == NULL) {
      if (m_swapper != NULL)
	m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    else {

      m_swapsDone += pM2D->pAnisoRef->iAniSwap(faceA)
	  + pM2D->pAnisoRef->iAniSwap(faceB) + pM2D->pAnisoRef->iAniSwap(faceD);
    }

    // TODO: Need a general hook to set length scale here.
    //   if(pAnisoRef!=NULL){
    //     pAnisoRef->vInterpolateMetric(newVert);
    //   }
    //   spVUpdateLS.insert(newVert);

    // TODO: Observer fixes this, too.
    //   if(pCellQueue!=NULL){
    //     pCellQueue->vAddCell(cell1);
    //     pCellQueue->vAddCell(cell2);
    //   }

    //  m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
//	pM2D->markForLengthScaleCheck(newVert);
    return newVert;
  }

//@ Insert a point on a face separating two regions in the mesh.
  Vert*
  SwappingInserter2D::insertPointOnInternalBdryFace(const double newPtLoc_arg[],
						    Face* const face,
						    Vert* const existingVert)
// Note: Geometrically, pVNew may not fall exactly on pF because
// of round-off.  Once things are done, though, the connectivity will
// be the same as if it had.
  {
    //  Before:               After:
    //              B                    B
    //             /|\                  /|\         .
    //         A1 / | \ A2          A1 / | \ A2
    //           /  |  \              /  MB \       .
    //          /   |   \            /1  |  4\      .
    //         /    |    \          /_NC_N_ND_\     .
    //        C\  1 | 2  /D        C\    |    /D
    //          \   |   /            \2  |  3/
    //           \  |  /              \  MA /
    //         B1 \ | / B2          B1 \ | / B2
    //             \|/                  \|/
    //              A                    A
    //

    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(m_pMesh);
    assert(pM2D);

    BdryEdgeBase* bdryFace = dynamic_cast<BdryEdgeBase*>(face->getLeftCell());
    if (!bdryFace)
      bdryFace = dynamic_cast<BdryEdgeBase*>(face->getRightCell());
    assert(bdryFace->isValid());
    TriCell *triCell1 = dynamic_cast<TriCell *>(face->getOppositeCell(bdryFace));
    if (!triCell1->isValid())
      return pVInvalidVert; // Interior cell was a quad.

    assert(face->doFullCheck());
    assert(face->getFaceLoc() == Face::eBdryTwoSide);

    assert(triCell1->doFullCheck());
    assert(face->getLeftCell() == triCell1 || face->getRightCell() == triCell1);
    assert(triCell1->getType() == Cell::eTriCell);
    assert(bdryFace->getType() == Cell::eIntBdryEdge);

    Face* otherFace = bdryFace->getFace(0);
    if (otherFace == face)
      otherFace = bdryFace->getFace(1);

    TriCell *triCell2 = dynamic_cast<TriCell*>(otherFace->getOppositeCell(
	bdryFace));
    if (!triCell2->isValid())
      return pVInvalidVert; // Interior cell was a quad.
    assert(otherFace->doFullCheck());
    assert(otherFace->getFaceLoc() == Face::eBdryTwoSide);

    assert(triCell2->doFullCheck());
    assert(
	otherFace->getLeftCell() == triCell2
	    || otherFace->getRightCell() == triCell2);
    assert(triCell2->getType() == Cell::eTriCell);
    //    assert (triCell1->getRegion () != triCell2->getRegion ());

    m_nOnIntBdryEdge++;

    double newPtLoc[3];
    if (forceInsertion()) {
      newPtLoc[0] = newPtLoc_arg[0];
      newPtLoc[1] = newPtLoc_arg[1];
    }
    else {
      CubitVector newLoc;
      bool isConcentric;
      bdryFace->getSplitData(m_splitType, newLoc, isConcentric);
      newPtLoc[0] = newLoc[0];
      newPtLoc[1] = newLoc[1];
    }
    Vert *newVert = getNewVert(newPtLoc, pM2D, existingVert);
    assert(newVert->isValid());

    Vert *pVA = face->getVert(0);
    Vert *pVB = face->getVert(1);
    Vert *pVC = triCell1->getOppositeVert(face);
    Vert *pVD = triCell2->getOppositeVert(otherFace);

    if (checkOrient2D(pVA, pVB, pVC) == -1) {
      assert(checkOrient2D(pVB, pVA, pVD) == -1);
      Vert *pVtmp = pVA;
      pVA = pVB;
      pVB = pVtmp;
    }

    assert(checkOrient2D(pVA, pVB, pVC) == 1);
    assert(checkOrient2D(pVB, pVA, pVD) == 1);

    assert(checkOrient2D(pVA, pVD, newVert) == 1);
    assert(checkOrient2D(pVD, pVB, newVert) == 1);
    assert(checkOrient2D(pVB, pVC, newVert) == 1);
    assert(checkOrient2D(pVC, pVA, newVert) == 1);

    // Need to identify faces A1, A2, B1, B2.
    Face *pFA1 = triCell1->getOppositeFace(pVA);
    Face *pFB1 = triCell1->getOppositeFace(pVB);
    Face *pFA2 = triCell2->getOppositeFace(pVA);
    Face *pFB2 = triCell2->getOppositeFace(pVB);

    int iReg1 = triCell1->getRegion();
    int iReg2 = triCell2->getRegion();

    pM2D->deleteCell(triCell1);
    pM2D->deleteCell(triCell2);

    // Now create new faces on the internal bdry, in duplicate.
    bool qExist;
    Face *pFMA1 = pM2D->createFace(qExist, pVA, newVert);
    assert(!qExist);
    Face *pFMA2 = pM2D->createFace(qExist, pVA, newVert, true);
    assert(qExist);

    Face *pFMB1 = pM2D->createFace(qExist, pVB, newVert);
    assert(!qExist);
    Face *pFMB2 = pM2D->createFace(qExist, pVB, newVert, true);
    assert(qExist);

    // Now the new "normal" faces.
    Face *pFNC = pM2D->createFace(qExist, pVC, newVert);
    assert(!qExist);
    Face *pFND = pM2D->createFace(qExist, pVD, newVert);
    assert(!qExist);

    //@@ Create new cells.
    triCell1 = pM2D->createTriCell(pFMB1, pFA1, pFNC, iReg1);
    triCell2 = pM2D->createTriCell(pFNC, pFB1, pFMA1, iReg1);
    TriCell *triCell3 = pM2D->createTriCell(pFMA2, pFB2, pFND, iReg2);
    TriCell *triCell4 = pM2D->createTriCell(pFND, pFA2, pFMB2, iReg2);

    // The internal bdry faces must be created after the cells.
    pM2D->createIntBdryEdge(pFMA1, pFMA2, bdryFace);
    pM2D->createIntBdryEdge(pFMB1, pFMB2, bdryFace);

    // Now the old internal bdry face can be deleted.  This auto-deletes
    // the old faces connected to it.
    pM2D->deleteBFace(bdryFace);

    assert(triCell1->doFullCheck());
    assert(triCell2->doFullCheck());
    assert(triCell3->doFullCheck());
    assert(triCell4->doFullCheck());
    assert(pFMA1->doFullCheck());
    assert(pFMB1->doFullCheck());
    assert(pFMA2->doFullCheck());
    assert(pFMB2->doFullCheck());
    assert(pFA1->doFullCheck());
    assert(pFB1->doFullCheck());
    assert(pFA2->doFullCheck());
    assert(pFB2->doFullCheck());
    assert(pFNC->doFullCheck());
    assert(pFND->doFullCheck());
    if (m_swapper != NULL)
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    return newVert;
  }

//@ Insert a point on an interior face
  Vert*
  SwappingInserter2D::insertPointOnFace(const double newPtLoc_arg[],
					Face* const face,
					Vert* const existingVert)
// Note: Geometrically, the new vert may not fall exactly on the face
// because of round-off.  Once things are done, though, the
// connectivity will be the same as if it had.
  {
    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(m_pMesh);
    assert(pM2D);

    assert(face->doFullCheck());
    Cell *cell1 = face->getLeftCell();
    Cell *cell2 = face->getRightCell();

    assert(cell1->doFullCheck());
    assert(cell2->doFullCheck());

    if (cell1->getType() == Cell::eBdryEdge
	|| cell2->getType() == Cell::eBdryEdge) {
      assert(face->getFaceLoc() == Face::eBdryFace);
      return insertPointOnBdryFace(newPtLoc_arg, face, existingVert);
    }

    if (cell1->getType() == Cell::eIntBdryEdge
	|| cell2->getType() == Cell::eIntBdryEdge) {
      assert(face->getFaceLoc() == Face::eBdryTwoSide);
      return insertPointOnBdryFace(newPtLoc_arg, face, existingVert);
    }

    assert(cell1->getRegion() == cell2->getRegion());
    if ((cell1->getType() != Cell::eTriCell)
	|| (cell2->getType() != Cell::eTriCell))
      return pVInvalidVert;

    m_nOnEdge++;

    // Create the new vertex.
    double newPtLoc[2];
    if (forceInsertion()) {
      newPtLoc[0] = newPtLoc_arg[0];
      newPtLoc[1] = newPtLoc_arg[1];
    }

    else {
      newPtLoc[0] = (face->getVert(0)->x() + face->getVert(1)->x()) * 0.5;
      newPtLoc[1] = (face->getVert(0)->y() + face->getVert(1)->y()) * 0.5;
    }
    Vert *newVert = getNewVert(newPtLoc, pM2D, existingVert);
    newVert->setType(Vert::eInterior);
    assert(newVert->isValid());

    TriCell * const triCell1 = dynamic_cast<TriCell *>(cell1);
    TriCell * const triCell2 = dynamic_cast<TriCell *>(cell2);

    Vert *vertA = face->getVert(0);
    Vert *vertB = face->getVert(1);
    Vert *vertC = triCell1->getOppositeVert(face);
    Vert *vertD = triCell2->getOppositeVert(face);

    assert(checkOrient2D(vertA, vertB, vertC) == 1);
    assert(checkOrient2D(vertB, vertA, vertD) == 1);

    assert(checkOrient2D(vertA, vertD, newVert) == 1);
    assert(checkOrient2D(vertD, vertB, newVert) == 1);
    assert(checkOrient2D(vertB, vertC, newVert) == 1);
    assert(checkOrient2D(vertC, vertA, newVert) == 1);

    //@@ Create new cells
    int iReg = cell1->getRegion();
    pM2D->deleteCell(cell1);
    pM2D->deleteCell(cell2);

    cell1 = pM2D->createTriCell(vertA, vertD, newVert, iReg);
    cell2 = pM2D->createTriCell(vertD, vertB, newVert, iReg);
    Cell* cell3 = pM2D->createTriCell(vertB, vertC, newVert, iReg);
    Cell* cell4 = pM2D->createTriCell(vertC, vertA, newVert, iReg);

    assert(cell1->doFullCheck());
    assert(cell2->doFullCheck());
    assert(cell3->doFullCheck());
    assert(cell4->doFullCheck());

    // Length scale and insertion queue handling goes elsewhere.
    if (m_swapper != NULL) {
      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
    }
    return newVert;
  }

//@ Insert a point inside a cell
  Vert*
  SwappingInserter2D::insertPointInCell(const double* newPtLoc, Cell * cell,
					Vert* const existingVert)
  {
    assert(cell->doFullCheck());
    assert(cell->getType() == Cell::eTriCell);

    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(m_pMesh);
    assert(pM2D);

    int iIsInTriCheck = isInTri(newPtLoc, cell->getVert(0)->getCoords(),
				cell->getVert(1)->getCoords(),
				cell->getVert(2)->getCoords());
    switch (iIsInTriCheck)
      {
      case 2: // On an edge.
	      //   printf("Finding Face...");
	for (int iF = 0; iF < 3; iF++) {
	  Face *pF = cell->getFace(iF);
	  if (checkOrient2D(pF->getVert(0)->getCoords(),
			    pF->getVert(1)->getCoords(), newPtLoc) == 0) {
	    //	 printf("found\n");
	    return insertPointOnFace(newPtLoc, pF, existingVert);
	  }
	  assert(iF != 2);
	}
	break;
      case 3: // Interior of the triangle
	{
	  m_nInCell++;
	  if (pM2D->pCellQueue != NULL) {
	    //  printf("Removing cell-%X from queue\n",pC);
	    pM2D->pCellQueue->vRemoveCell(cell);
	  }

	  // Verts ABCD are ordered in a right-handed way.  The faces are
	  // defined relative to them so that there is no possible ambiguity.
	  // Note that the face definitions here DO NOT correspond to
	  // iFace[ABCD](iCell).
	  Vert *pVA = cell->getVert(0);
	  Vert *pVB = cell->getVert(1);
	  Vert *pVC = cell->getVert(2);

	  // Get new cells
	  int iReg = cell->getRegion();
	  pM2D->deleteCell(cell);

	  Vert *newVert = getNewVert(newPtLoc, pM2D, existingVert);

	  // Create new cells
	  Cell *cellA = pM2D->createTriCell(pVA, pVB, newVert, iReg);
	  Cell *cellB = pM2D->createTriCell(pVB, pVC, newVert, iReg);
	  Cell *cellC = pM2D->createTriCell(pVC, pVA, newVert, iReg);
	  // Connectivity for old faces gets updated automatically.

	  assert(cellA->doFullCheck());
	  assert(cellB->doFullCheck());
	  assert(cellC->doFullCheck());
	  newVert->setType(Vert::eInterior);

	  if (pM2D->pAnisoRef != NULL) {
	    pM2D->pAnisoRef->vInterpolateMetric(newVert);
//	    m_length->setLengthScale(newVert);

	    //      printf("Setting 3 new Lengths\n");
	  }

	  if (pM2D->pCellQueue != NULL) {

	    pM2D->pCellQueue->vAddCell(cellA);
	    pM2D->pCellQueue->vAddCell(cellB);
	    pM2D->pCellQueue->vAddCell(cellC);
	    //    printf("Added to the queue\n");
	  }

	  if (pM2D->pAnisoRef == NULL) {
	    if (m_swapper != NULL)
	      m_swapsDone += m_swapper->swapAllQueuedFacesRecursive();
	  }
	  else {

	    Face *FA = cellA->getOppositeFace(pVA);
	    Face *FB = cellB->getOppositeFace(pVB);
	    Face *FC = cellC->getOppositeFace(pVC);

	    Face *FD = cellA->getOppositeFace(newVert);
	    Face *FE = cellB->getOppositeFace(newVert);
	    Face *FF = cellC->getOppositeFace(newVert);
	    m_swapsDone += pM2D->pAnisoRef->iAniSwap(FA)
		+ pM2D->pAnisoRef->iAniSwap(FB) + pM2D->pAnisoRef->iAniSwap(FC)
		+ pM2D->pAnisoRef->iAniSwap(FD) + pM2D->pAnisoRef->iAniSwap(FE)
		+ pM2D->pAnisoRef->iAniSwap(FF);
	  }
	  return newVert;
	}
	break;
      case 1:
	{
	  // Case 1: New vert corresponds with an old one.
	  for (int iV = 0; iV < 3; iV++) {
	    Vert *pV = cell->getVert(iV);
	    if (dDIST2D(pV->getCoords(),newPtLoc) < 1.e-10)
	      return pV;

	  }
	  assert(0);
	}
	break;
      default:
	// Case 0: Triangle is collinear.
	// Case -1:  New vert is outside.
	// In none of these cases does insertion make any sense.

	assert(0);
	break;
      }
    return pVInvalidVert;
  }

}
