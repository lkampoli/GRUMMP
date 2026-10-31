#include "GR_config.h"
#include "GR_Classes.h"

#include "GR_BFace.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"
#include "GR_Geometry.h"

// Raw creation of new entities, with default data.

//@ Return an empty face with a specified number of verts.
Face*
VolMesh::getNewFace(const int iNV, int ID)
{
  assert(iNV == 3 || iNV == 4);
  assert(ID == omp_get_thread_num());
  switch (iNV)
    {
    case 3:
      {
	TriFace *pTF = m_ECTriF.getNewEntry(ID);
	pTF->setDefaultFlags();
	return pTF;
      }
    case 4:
      {
	QuadFace *pQF = m_ECQuadF.getNewEntry(ID);
	pQF->setDefaultFlags();
	return pQF;
      }
    default:
      assert(0);
      return (pFInvalidFace);
    }
}

//@ Return an empty bdry face with a specified number of verts.
BFace*
VolMesh::getNewBFace(const int iNBV, const bool qIsInternal, int ID)
{
  assert(iNBV == 3 || iNBV == 4);
  assert(ID == omp_get_thread_num());
  BFace *pBFRetVal = pBFInvalidBFace;
  switch (iNBV)
    {
    case 3:
      if (qIsInternal) {
	pBFRetVal = m_ECIntTriBF.getNewEntry(ID);
	pBFRetVal->setType(Cell::eIntTriBFace);
      }
      else {
	pBFRetVal = m_ECTriBF.getNewEntry(ID);
	pBFRetVal->setType(Cell::eTriBFace);
      }
      dynamic_cast<TriBFaceBase*>(pBFRetVal)->clearGeom();
      break;
    case 4:
      if (qIsInternal) {
	pBFRetVal = m_ECIntQuadBF.getNewEntry(ID);
	pBFRetVal->setType(Cell::eIntQuadBFace);
      }
      else {
	pBFRetVal = m_ECQuadBF.getNewEntry(ID);
	pBFRetVal->setType(Cell::eQuadBFace);
      }
      dynamic_cast<QuadBFaceBase*>(pBFRetVal)->clearGeom();
      break;
    default:
      assert(0);
      return (pBFInvalidBFace);
    }
  assert(pBFRetVal->isValid());
  pBFRetVal->resetAllData();
  return pBFRetVal;
}

//@ Return an empty cell with the specified number of faces and verts.
Cell*
VolMesh::getNewCell(const int iNF, const int iNV, int ID)
{
  assert((iNF == 4 && iNV == 4) || // Tetrahedron
      (iNF == 5 && (iNV == 5 || iNV == 6)) || // Pyramid and prism
      (iNF == 3 && iNV == 4) || // Flake
      (iNF == 6 && iNV == 8)); // Hexahedron
  assert(ID == omp_get_thread_num());
  switch (iNV)
    {
    case 4:
      {
	if (iNF == 4) {
	  TetCell *pTet = m_ECTet.getNewEntry(ID);
	  pTet->setDefaultFlags();
	  pTet->setType(Cell::eTet);
	  return pTet;
	}
	else {
	  FlakeCell* pFlake = m_ECFlake.getNewEntry(ID);
	  pFlake->setDefaultFlags();
	  pFlake->setType(Cell::eFlake);
	  return pFlake;
	}
      }
    case 5:
      {
	PyrCell *pPyr = m_ECPyr.getNewEntry(ID);
	pPyr->setDefaultFlags();
	pPyr->setType(Cell::ePyr);
	return pPyr;
      }
    case 6:
      {
	PrismCell *pPrism = m_ECPrism.getNewEntry(ID);
	pPrism->setDefaultFlags();
	pPrism->setType(Cell::ePrism);
	return pPrism;
      }
    case 8:
      {
	HexCell *pHex = m_ECHex.getNewEntry(ID);
	pHex->setDefaultFlags();
	pHex->setType(Cell::eHex);
	return pHex;
      }
    default:
      assert(0);
      return (pCInvalidCell);
    }
}

// The semantics of these calls are essentially identical to the ITAPS
// createEnt calls.
//   createFace returns a face with verts pV0 and pV1; if the face
//     already existed it might be opposite in sense, but that's OK.
//   createCell returns a cell with given faces.  In GRUMMP (though not
//     in ITAPS) that cell can't have existed already.
//   createCell returns a cell with given verts.  Faces are created as
//     needed.
//
// At present, canonical ordering of the input data is assumed for faces
// and vert->tri.  For vert->quad, it's mandatory, to prevent
// awkwardness like introducing edges on diagonals of a quad.

Face*
VolMesh::createFace(bool &qAlreadyExisted, Vert * const pV0, Vert * const pV1,
		    Vert * const pV2, Vert * const pV3, bool qForceDuplicate,
		    int ID)
{
// Edge3D: For now, this happens directly.  Once edges are added,
// Edge3D: it'll look a lot more like Mesh2D::createCell.

// iMesh merge: need to return both the face ptr and a flag indicating
// whether the face already existed.
  assert(pV0->isValid());
  assert(pV1->isValid());
  assert(pV2->isValid());
  Face *pF = pFInvalidFace;
  if (!pV3->isValid()) {
    // Creating a tri.
    // Check for prior existence.
    pF = findCommonFace(pV0, pV1, pV2, NULL, true);
    qAlreadyExisted = (pF != pFInvalidFace);
    // If we're forcing a duplicate, we'd better have a copy already!
    assert((!qForceDuplicate) || (qAlreadyExisted));
    if (qForceDuplicate == qAlreadyExisted) {
      // Only do this if we've got a correct duplication case or a
      // correct creation case.
      pF = getNewFace(3, ID);
      pF->setVerts(pV0, pV1, pV2);
      pF->setLeftCell(NULL);
      pF->setRightCell(NULL);
//	#pragma omp critical(Event3D)
//      {
      createFaceEvent(pF);
//      }
    }
  }
  else {
    // Creating a quad.
    m_simplexMesh = false;
    // Check for prior existence.
    pF = findCommonFace(pV0, pV1, pV2, pV3, true);
    qAlreadyExisted = (pF != pFInvalidFace);
    // If we're forcing a duplicate, we'd better have a copy already!
    assert((!qForceDuplicate) || (qAlreadyExisted));
    if (qForceDuplicate == qAlreadyExisted) {
      // Only do this if we've got a correct duplication case or a
      // correct creation case.
      pF = getNewFace(4, ID);
      pF->setVerts(pV0, pV1, pV2, pV3);
      pF->setLeftCell(NULL);
      pF->setRightCell(NULL);
//	#pragma omp critical(Event3D)
//      {
      createFaceEvent(pF);
//      }
    }
  }
// No checking possible here, since the face is created before the cell...
  return pF;
}

BFace*
VolMesh::createBFace(Face * const pF, BFace * const pBFOld)
{
  int ID = 0;
  if (omp_in_parallel() != 0)
    ID = omp_get_thread_num();
// Face must already exist
  assert(pF->isValid());
  BFace *pBF = pBFInvalidBFace;
  if (pF->getNumCells() == 1) {
    // Must have exactly one open connectivity slot for the face!
    assert(
	(pF->getLeftCell()->isValid() || pF->getRightCell()->isValid())
	    && !(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid()));

    pBF = getNewBFace(pF->getNumVerts(), false, ID);

    pBF->assignFaces(pF);
    pF->addCell(pBF);
    pF->setFaceLoc(Face::eBdryFace);
    assert(pBF->doFullCheck());
    assert(pBF->getFace() == pF);
    assert(pF->getLeftCell() == pBF || pF->getRightCell() == pBF);
//    for (int ii = 0; ii < pF->getNumVerts(); ++ii) {
//        Vert *pV = pF->getVert(ii);
//        switch (pV->getVertType()) {
//            case Vert::eUnknown:
//            case Vert::ePseudoSurface:
//            case Vert::eInterior:
//            case Vert::eInteriorFixed:
//                pV->setType(Vert::eBdry);
//                break;
//            case Vert::eBBox:
//                break;
//            default:
//                break;
//        }
//    }
  }
  else {
    // This is an internal bdry face.  The face must have cells on both sides.
//    assert(pBFOld->getType() == Cell::eIntTriBFace ||
//	   pBFOld->getType() == Cell::eIntQuadBFace);
    assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());

    bool qExists;
    Face *pFDup;
    if (pF->getNumVerts() == 3) {
      pFDup = createFace(qExists, pF->getVert(0), pF->getVert(1),
			 pF->getVert(2),
			 pVInvalidVert,
			 true, ID);
    }
    else {
      pFDup = createFace(qExists, pF->getVert(0), pF->getVert(1),
			 pF->getVert(2), pF->getVert(3), true, ID);
    }
    assert(qExists);
    pBF = getNewBFace(pF->getNumVerts(), true, ID);
    pBF->assignFaces(pF, pFDup);

    Cell *pCR = pF->getRightCell();

    pFDup->setLeftCell(pBF);
    pFDup->setRightCell(pCR);
    pF->setRightCell(pBF);
    pCR->replaceFace(pF, pFDup);

    pF->setFaceLoc(Face::eBdryTwoSide);
    pFDup->setFaceLoc(Face::eBdryTwoSide);
    assert(pBF->doFullCheck());
    assert(pCR->doFullCheck());
    assert(pF->getLeftCell()->doFullCheck());
    assert(pFDup->getRightCell()->doFullCheck());

//    for (int ii = 0; ii < pF->getNumVerts(); ++ii) {
//        Vert *pV = pF->getVert(ii);
//        switch (pV->getVertType()) {
//            case Vert::eUnknown:
//            case Vert::ePseudoSurface:
//            case Vert::eInterior:
//            case Vert::eInteriorFixed:
//                pV->setType(Vert::eBdryTwoSide);
//                break;
//            case Vert::eBBox:
//                assert(0); // This should never happen.
//                break;
//            case Vert::eBdry:
//                pV->setType(Vert::eBdryCurve);
//                break;
//            default:
//                break;
//        }
//    }
  }

#ifndef OMIT_VERTEX_HINTS
// Make sure that bdry verts have bdry hints
  for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
    pF->getVert(ii)->setHintFace(pF);
  }
#endif

// This is a really stinky way to copy bdry data, but that's what
// works currently.

  if (pBFOld->isValid()) {
    if (pBFOld->getPatchPointer()) {
      pBF->setPatch(pBFOld->getPatchPointer());
    }
    else if (pBFOld->getTopologicalParent()) {
      assert(pBFOld->getTopologicalParent());
      pBF->setTopologicalParent(pBFOld->getTopologicalParent());

    }
  }

//#pragma omp critical(Event3D)
//  {
  createBFaceEvent(pBF);
//  }
  return pBF;
}
BFace*
VolMesh::createBFaceParallel(Face * const pF, BFace * const pBFOld, int ID)
{
// Face must already exist
  assert(pF->isValid());
  BFace *pBF = pBFInvalidBFace;
  if (pF->getNumCells() == 1) {
    // Must have exactly one open connectivity slot for the face!
    assert(
	(pF->getLeftCell()->isValid() || pF->getRightCell()->isValid())
	    && !(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid()));

    pBF = getNewBFace(pF->getNumVerts(), false, ID);

    pBF->assignFaces(pF);
    pF->addCell(pBF);
    pF->setFaceLoc(Face::eBdryFace);
    assert(pBF->doFullCheck());
    assert(pBF->getFace() == pF);
    assert(pF->getLeftCell() == pBF || pF->getRightCell() == pBF);
  }
  else {
    // This is an internal bdry face.  The face must have cells on both sides.
    //    assert(pBFOld->getType() == Cell::eIntTriBFace ||
    //	   pBFOld->getType() == Cell::eIntQuadBFace);
    assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());

    bool qExists;
    Face *pFDup;
    if (pF->getNumVerts() == 3) {
      pFDup = createFace(qExists, pF->getVert(0), pF->getVert(1),
			 pF->getVert(2),
			 pVInvalidVert,
			 true, ID);
    }
    else {
      pFDup = createFace(qExists, pF->getVert(0), pF->getVert(1),
			 pF->getVert(2), pF->getVert(3), true, ID);
    }
    assert(qExists);
    pBF = getNewBFace(pF->getNumVerts(), true, ID);
    pBF->assignFaces(pF, pFDup);

    Cell *pCR = pF->getRightCell();

    pFDup->setLeftCell(pBF);
    pFDup->setRightCell(pCR);
    pF->setRightCell(pBF);
    pCR->replaceFace(pF, pFDup);

    pF->setFaceLoc(Face::eBdryTwoSide);
    pFDup->setFaceLoc(Face::eBdryTwoSide);
    assert(pBF->doFullCheck());
    assert(pCR->doFullCheck());
    assert(pF->getLeftCell()->doFullCheck());
    assert(pFDup->getRightCell()->doFullCheck());
  }

#ifndef OMIT_VERTEX_HINTS
// Make sure that bdry verts have bdry hints
  for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
    pF->getVert(ii)->setHintFace(pF);
  }
#endif

// This is a really stinky way to copy bdry data, but that's what
// works currently.

  if (pBFOld->isValid()) {
    if (pBFOld->getPatchPointer()) {
      pBF->setPatch(pBFOld->getPatchPointer());
    }
    else if (pBFOld->getTopologicalParent()) {
      pBF->setTopologicalParent(pBFOld->getTopologicalParent());
      //  		printf("old face %p %d\n",pBFOld,pBFOld->isDeleted());
      //  			  pBFOld->getFace(0)->printFaceInfo();
      //  		    BasicTopologyEntity * te = pBFOld->getTopologicalParent();
      //  		    CubitVector cv = te->center_point();
      //  		    printf("te-cp %p %f %f %f\n",te,cv.x(),cv.y(),cv.z());
      //  				printf("new face %p %d\n",pBF,pBF->isDeleted());
      //  			  pBF->getFace(0)->printFaceInfo();
      //  		     te = pBF->getTopologicalParent();
      //  		     cv = te->center_point();
      //  		    printf("te-cp %p %f %f %f\n",te,cv.x(),cv.y(),cv.z());
    }
    //     else {
    //       vFatalError("No surface geometry stored for this boundary face",
    // 		  "VolMesh::createBFace");
    //     }
  }
//#pragma omp critical(Event3D)
//  {
  createBFaceEvent(pBF);
//  }
  return pBF;
}

TetCell*
VolMesh::createTetCell(bool &qExistedAlready, Face * const pF0,
		       Face * const pF1, Face * pF2, Face * pF3, const int iReg,
		       int ID)
// pF2 and pF3 may be changed internally, but this doesn't affect data
// outside the routine, since the pointers are passed by value.
//
// This call creates a tet.
//
// Faces are stored in canonical order, even if they aren't passed that
// way.  Also, face->cell connectivity is set up as part of this
// function.
{
// Faces must already exist
  assert(pF0->isValid());
  assert(pF1->isValid());
  assert(pF2->isValid());
  assert(pF3->isValid());

// Set up the face->cell connectivity properly.  The canonical tet
// looks like this, where face 3 is the bottom (see ITAPS docs):
//
//                3
//               /|\         .
//              / | \        .
//             /  |2 \       .
//            /   |   \      .
//           0 - - - - 2
//            \ 0 | 1 /
//             \ (3) /
//              \ | /
//               \|/
//                1

  assert(pF0->getNumVerts() == 3);
  assert(pF1->getNumVerts() == 3);
  assert(pF2->getNumVerts() == 3);
  assert(pF3->getType() == Face::eTriFace);

  Cell *pCCom01 = findCommonCell(pF0, pF1);
  if (pCCom01->isValid()) {
    Cell *pCCom23 = findCommonCell(pF2, pF3);
    if (pCCom23->isValid() && pCCom23 == pCCom01) {
      qExistedAlready = true;
      TetCell* pTC = dynamic_cast<TetCell*>(pCCom01);
      assert(pTC->isValid()); // It'd better be a tet!
      return pTC;
    }
  }
  qExistedAlready = false;

// Each trio of faces  must have a common vert.
  Vert *apV[4];
  apV[0] = findCommonVert(pF0, pF2, pF3);
  apV[1] = findCommonVert(pF0, pF1, pF3);
  apV[2] = findCommonVert(pF1, pF2, pF3);
  apV[3] = findCommonVert(pF0, pF1, pF2);
  assert(apV[0]->isValid());
  assert(apV[1]->isValid());
  assert(apV[2]->isValid());
  assert(apV[3]->isValid());

  TetCell *pTC = dynamic_cast<TetCell*>(getNewCell(4, 4, ID));

  pTC->setRegion(iReg);
  int iOrient = checkOrient3D(apV[0], apV[1], apV[2], apV[3]);
//  assert(iOrient != 0);
  if (iOrient == -1) {
    Vert *tempV = apV[3];
    apV[3] = apV[1];
    apV[1] = tempV;
    Face *tempF = pF2;
    pF2 = pF3;
    pF3 = tempF;
    iOrient = 1;
    assert(-1 != checkOrient3D(apV[0], apV[1], apV[2], apV[3]));
  } 	// Now right-handed (as sketched above)
  pTC->assignFaces(pF0, pF1, pF2, pF3);

  if (dynamic_cast<TriFace*>(pF0)->isRightHanded(apV[0], apV[3], apV[1])) {
    assert(!pF0->getRightCell()->isValid());
    pF0->setRightCell(pTC);
  }
  else {
    assert(!pF0->getLeftCell()->isValid());
    pF0->setLeftCell(pTC);
  }
  if (dynamic_cast<TriFace*>(pF1)->isRightHanded(apV[1], apV[3], apV[2])) {
    assert(!pF1->getRightCell()->isValid());
    pF1->setRightCell(pTC);
  }
  else {
    assert(!pF1->getLeftCell()->isValid());
    pF1->setLeftCell(pTC);
  }
  if (dynamic_cast<TriFace*>(pF2)->isRightHanded(apV[2], apV[3], apV[0])) {
    assert(!pF2->getRightCell()->isValid());
    pF2->setRightCell(pTC);
  }
  else {
    assert(!pF2->getLeftCell()->isValid());
    pF2->setLeftCell(pTC);
  }
  if (dynamic_cast<TriFace*>(pF3)->isRightHanded(apV[0], apV[1], apV[2])) {
    assert(!pF3->getRightCell()->isValid());
    pF3->setRightCell(pTC);
  }
  else {
    assert(!pF3->getLeftCell()->isValid());
    pF3->setLeftCell(pTC);
  }

  assert(pTC->doFullCheck());
  assert(pTC->hasFace(pF0));
  assert(pTC->hasFace(pF1));
  assert(pTC->hasFace(pF2));
  assert(pTC->hasFace(pF3));
  assert(pF0->getLeftCell() == pTC || pF0->getRightCell() == pTC);
  assert(pF1->getLeftCell() == pTC || pF1->getRightCell() == pTC);
  assert(pF2->getLeftCell() == pTC || pF2->getRightCell() == pTC);
  assert(pF3->getLeftCell() == pTC || pF3->getRightCell() == pTC);
  assert(
      checkOrient3D(pTC->getVert(0), pTC->getVert(1), pTC->getVert(2),
		    pTC->getVert(3)) == 1);
//	#pragma omp critical(Event3D)
//  {
  createCellEvent(pTC);
//  }
  return pTC;
}

// Cell* VolMesh::createCell(bool &qExistedAlready,
// 			  Face * const pF0, Face * const pF1,
// 			  Face * const pF2, Face * const pF3,
// 			  Face * const pF4, const int iReg)
// // This one would be for pyramids / prisms.  Faces must be canonically
// // ordered.
// {
//   // Faces must already exist
//   assert(pF0->qValid());
//   assert(pF1->qValid());
//   assert(pF2->qValid());
//   assert(pF3->qValid());
//   assert(pF4->qValid());

//   Cell *pC;
//   if (pF0->eType() == Face::eTriFace) {
//     // This must be a pyramid
//     assert(pF0->eType() == Face::eTriFace);
//     assert(pF1->eType() == Face::eTriFace);
//     assert(pF2->eType() == Face::eTriFace);
//     assert(pF3->eType() == Face::eTriFace);
//     assert(pF4->eType() == Face::eQuadFace);

//     pC = pCNewCell(5, 5);

//     // The pyramid looks like this:
//     //
//     //                  4
//     //                 /|\               .
//     //                / | \              .
//     //               / |\  \             .
//     //              /  / |  \            .
//     //             /  |  |   \           .
//     //            /   /  \    \          .
//     //           /   |_---3-_  \         .
//     //          0--- /       ---2        .
//     //           \_ |    ---
//     //             \1---
//     //
//     // Okay, so it's a crappy pyramid; so what?  Faces are:
//     //  0:  014
//     //  1:  124
//     //  2:  234
//     //  3:  304
//     //  4:  0123

//     Vert *apV[5];
//     pF4->vAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apV));
//     apV[4] = pF0->pVVert(0);
//     if (pF4->qHasVert(apV[4])) {
//       apV[4] = pF0->pVVert(1);
//       if (pF4->qHasVert(apV[4])) {
// 	apV[4] = pF0->pVVert(2);
//       }
//     }
//     assert(!pF4->qHasVert(apV[4]));

//     // If this orientation is negative, then surely it means that the
//     // bottom face is oriented out.
//     int iOrient = checkOrient3D(apV[0], apV[1], apV[2], apV[4]);
//     assert(checkOrient3D(apV[0], apV[2], apV[3], apV[4]) == iOrient);

//     if (iOrient == 1) {
//       pF4->vSetRightCell(pC);
//     }
//     else {
//       pF4->vSetLeftCell(pC);
//       Vert *pVTmp = apV[3];
//       apV[3] = apV[0];
//       apV[0] = pVTmp;
//       pVTmp = apV[2];
//       apV[2] = apV[1];
//       apV[1] = pVTmp;
//     }

//     assert(pF0->qHasVert(apV[0]) && pF0->qHasVert(apV[1]));
//     assert(pF1->qHasVert(apV[1]) && pF1->qHasVert(apV[2]));
//     assert(pF2->qHasVert(apV[2]) && pF2->qHasVert(apV[3]));
//     assert(pF3->qHasVert(apV[3]) && pF3->qHasVert(apV[0]));

//     // Assignment for the cell is invariant w/ the orientation; this is
//     // different from createTetCell, which doesn't require faces in
//     // canonical order.

//     if (dynamic_cast<TriFace*>(pF0)->qRightHanded(apV[0], apV[4], apV[1]))
//       { pF0->vSetRightCell(pC); }
//     else
//       { pF0->vSetLeftCell(pC); }

//     if (dynamic_cast<TriFace*>(pF1)->qRightHanded(apV[1], apV[4], apV[2]))
//       { pF1->vSetRightCell(pC); }
//     else
//       { pF1->vSetLeftCell(pC); }

//     if (dynamic_cast<TriFace*>(pF2)->qRightHanded(apV[2], apV[4], apV[3]))
//       { pF2->vSetRightCell(pC); }
//     else
//       { pF2->vSetLeftCell(pC); }

//     if (dynamic_cast<TriFace*>(pF3)->qRightHanded(apV[3], apV[4], apV[0]))
//       { pF3->vSetRightCell(pC); }
//     else
//       { pF3->vSetLeftCell(pC); }

//   }

//   else {
//     // This must be a prism
//     assert(pF0->eType() == Face::eQuadFace);
//     assert(pF1->eType() == Face::eQuadFace);
//     assert(pF2->eType() == Face::eQuadFace);
//     assert(pF3->eType() == Face::eTriFace);
//     assert(pF4->eType() == Face::eTriFace);

//     pC = pCNewCell(5, 6);
//   }

//   pC->vSetRegion(iReg);

//   // Be sure to assign these with canonical ordering.
//   pC->vAssign(pF0, pF1, pF2, pF3, pF4);
//   assert(pC->iFullCheck());
//   assert(pC->qHasFace(pF0));
//   assert(pC->qHasFace(pF1));
//   assert(pC->qHasFace(pF2));
//   assert(pC->qHasFace(pF3));
//   assert(pC->qHasFace(pF4));
//   assert(pF0->getLeftCell() == pC ||
// 	 pF0->getRightCell() == pC);
//   assert(pF1->getLeftCell() == pC ||
// 	 pF1->getRightCell() == pC);
//   assert(pF2->getLeftCell() == pC ||
// 	 pF2->getRightCell() == pC);
//   assert(pF3->getLeftCell() == pC ||
// 	 pF3->getRightCell() == pC);
//   assert(pF4->getLeftCell() == pC ||
// 	 pF4->getRightCell() == pC);
//   return pC;
// }

TetCell*
VolMesh::createTetCell(bool &qAlreadyExisted, Vert * const pV, Face * const pF,
		       const int iReg)
{
// This function first creates or finds existing faces connecting the
// vert to the three edges of the face, then uses those faces to build
// a tet.  Orientation tests are done in assignment anyway, so doing
// it here is a waste of time.  The only existing orientation check is
// there to ensure that we produce canonical ordering.
  int ID = 0;
  if (omp_in_parallel() != 0)
    ID = omp_get_thread_num();
  Vert *pV0 = pF->getVert(0);
  Vert *pV1 = pF->getVert(1);
  Vert *pV2 = pF->getVert(2);

  bool qExist02, qExist10, qExist21; // These values are ignored.
  Face *pF02, *pF10, *pF21;
  TetCell *pC;
  int iOrient = checkOrient3D(pV, pV0, pV1, pV2);
//assert(iOrient != 0);
  if (iOrient == 1) {
    if (!(pF->getRightCell()->isValid() && !pF->getLeftCell()->isValid())) {
      logMessage(1, "This vert ");
      pV->printVertInfo();
      logMessage(1, "Is badly oriented with this face \n");
      pV0->printVertInfo();
      pV1->printVertInfo();
      pV2->printVertInfo();
      assert(pF->getRightCell()->isValid() && !pF->getLeftCell()->isValid());
    }
    assert(pF->getRightCell()->isValid() && !pF->getLeftCell()->isValid());
    if (omp_in_parallel() == 0) {
      pF02 = createFace(qExist02, pV0, pV2, pV);
      pF10 = createFace(qExist10, pV1, pV0, pV);
      pF21 = createFace(qExist21, pV2, pV1, pV);
      pC = createTetCell(qAlreadyExisted, pF, pF02, pF10, pF21, iReg);
    }
    else {
      pF02 = createFace(qExist02, pV0, pV2, pV, NULL, false, ID);
      pF10 = createFace(qExist10, pV1, pV0, pV, NULL, false, ID);
      pF21 = createFace(qExist21, pV2, pV1, pV, NULL, false, ID);
      pC = createTetCell(qAlreadyExisted, pF, pF02, pF10, pF21, iReg, ID);
    }
  }
  else {
    if (!(!pF->getRightCell()->isValid() && pF->getLeftCell()->isValid())) {
      logMessage(1, "This vert ");
      pV->printVertInfo();
      logMessage(1, "Is badly oriented with this face \n");
      pV0->printVertInfo();
      pV1->printVertInfo();
      pV2->printVertInfo();
      assert(!pF->getRightCell()->isValid() && pF->getLeftCell()->isValid());
    }
    assert(!pF->getRightCell()->isValid() && pF->getLeftCell()->isValid());
    if (omp_in_parallel() == 0) {
      pF02 = createFace(qExist02, pV2, pV0, pV);
      pF10 = createFace(qExist10, pV0, pV1, pV);
      pF21 = createFace(qExist21, pV1, pV2, pV);
      pC = createTetCell(qAlreadyExisted, pF, pF21, pF10, pF02, iReg);
    }
    else {
      pF02 = createFace(qExist02, pV2, pV0, pV, NULL, false, ID);
      pF10 = createFace(qExist10, pV0, pV1, pV, NULL, false, ID);
      pF21 = createFace(qExist21, pV1, pV2, pV, NULL, false, ID);
      pC = createTetCell(qAlreadyExisted, pF, pF21, pF10, pF02, iReg, ID);
    }
  }

// Orientation is checked in createTetCell(faces).
  assert(
      !qExist02
	  || (pF02->getLeftCell()->isValid() && pF02->getRightCell()->isValid()));
  assert(
      !qExist10
	  || (pF10->getLeftCell()->isValid() && pF10->getRightCell()->isValid()));
  assert(
      !qExist21
	  || (pF21->getLeftCell()->isValid() && pF21->getRightCell()->isValid()));
  assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());

// Already registered the event in the low-level create tet call.
  return pC;
}

TetCell*
VolMesh::createTetCell(bool &qAlreadyExisted, Vert * const pV0,
		       Vert * const pV1, Vert * const pV2, Vert * const pV3,
		       const int iReg, bool OMP_Flag, int ID)
// This one creates tets from verts via faces.
{
  bool qExist103, qExist213, qExist023, qExist012;
  Face *pF023, *pF103, *pF213, *pF012;
  if (!OMP_Flag) {
    assert(omp_in_parallel() == 0);
    pF023 = createFace(qExist023, pV0, pV2, pV3);
    pF103 = createFace(qExist103, pV1, pV0, pV3);
    pF213 = createFace(qExist213, pV2, pV1, pV3);
    pF012 = createFace(qExist012, pV0, pV1, pV2);
  }
  else {
    if (NUM_PROCS > 1) {
      assert(omp_in_parallel() != 0);
      assert(ID == omp_get_thread_num());
    }
    pF023 = createFace(qExist023, pV0, pV2, pV3, NULL, false, ID);
    pF103 = createFace(qExist103, pV1, pV0, pV3, NULL, false, ID);
    pF213 = createFace(qExist213, pV2, pV1, pV3, NULL, false, ID);
    pF012 = createFace(qExist012, pV0, pV1, pV2, NULL, false, ID);
  }

// The following commented code was put in for aniso adaptation in 3D
// circa 2015.  Since that approach for 3D aniso adaptation didn't pan
// out, I'm taking this out.  CFOG, 7 April 2017.
/////////////////////////////////////////////////////////////////////////////
// Trying a similiar approach for orientation for Aniso adapt as with 2D
//	double dSize023 = pF023->calcSize();
//	double dSize103 = pF103->calcSize();
//	double dSize213 = pF213->calcSize();
//	double dSize012 = pF012->calcSize();
//
//	TetCell *pTC;
//
//	if (dSize023 >= dSize103 && dSize023 >= dSize213 && dSize023 >= dSize012) {
//		pTC = createTetCell(qAlreadyExisted, pF023, pF103, pF213, pF012, iReg,
//				ID);
//	} else if (dSize103 >= dSize023 && dSize103 >= dSize213
//			&& dSize103 >= dSize012) {
//		pTC = createTetCell(qAlreadyExisted, pF103, pF213, pF023, pF012, iReg,
//				ID);
//	} else if (dSize213 >= dSize023 && dSize213 >= dSize103
//			&& dSize213 >= dSize012) {
//		pTC = createTetCell(qAlreadyExisted, pF213, pF012, pF023, pF103, iReg,
//				ID);
//	} else {
//		pTC = createTetCell(qAlreadyExisted, pF012, pF023, pF213, pF103, iReg,
//				ID);
//	}
/////////////////////////////////////////////////////////////////////////////

  TetCell *pTC = createTetCell(qAlreadyExisted, pF023, pF103, pF213, pF012,
			       iReg, ID);

#ifndef NDEBUG
  assert(pTC->hasVert(pV0));
  assert(pTC->hasVert(pV1));
  assert(pTC->hasVert(pV2));
  assert(pTC->hasVert(pV3));

// Orientation is checked in createTetCell(faces).
#endif
// Already registered the event in the low-level create tet call.
  return (pTC);
}

static void
vQuadSetCellWithOrientation(Cell * const pC, Face * const pF,
			    const Vert * const pVRef)
{
  assert(pF->getType() == Face::eQuadFace);
// The quad, if it already exists, may be in either orientation.
// This could also be checked topologically, and maybe should be.

  Cell* leftCell = pF->getLeftCell();
  Cell* rightCell = pF->getRightCell();
  assert(!leftCell->isValid() || leftCell != pC);
  assert(!rightCell->isValid() || rightCell != pC);
  if (leftCell->isValid() || rightCell->isValid()) {
    if (leftCell->isValid()) {
      // Add to the right side
      assert(!rightCell->isValid());
      pF->setRightCell(pC);
    }
    else {
      // Add to the left side
      pF->setLeftCell(pC);
    }
  }
  else {
    // No topological hints; have to use geometry instead.  Can fail for
    // badly warped quads.  TODO: Fix this.
    Vert *apV[4];
    pF->getAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apV));

    // If this orientation is negative, then surely it means that the
    // bottom face is oriented out.
    int iOrient = checkOrient3D(apV[0], apV[1], apV[2], pVRef);
    if (!(checkOrient3D(apV[0], apV[2], apV[3], pVRef) == iOrient)) {
      cout << "Bad Quad Orientation: " << endl;
      cout << apV[0]->x() << ":" << apV[0]->y() << ":" << apV[0]->z() << endl;
      cout << apV[1]->x() << ":" << apV[1]->y() << ":" << apV[1]->z() << endl;
      cout << apV[2]->x() << ":" << apV[2]->y() << ":" << apV[2]->z() << endl;
      cout << apV[3]->x() << ":" << apV[3]->y() << ":" << apV[3]->z() << endl;
      assert(checkOrient3D(apV[0], apV[2], apV[3], pVRef) == iOrient);
    }

    if (iOrient == 1) {
      pF->setRightCell(pC);
    }
    else {
      pF->setLeftCell(pC);
    }
  }
}

static void
vTriSetCellWithOrientation(Cell * const pC, Face * const pF,
			   const Vert * const pVRef)
{
  if (pF->getType() != Face::eTriFace) {
    cout << "Wrong face type...?" << endl;
    for (int iV = 0; iV < pF->getNumVerts(); iV++) {
      Vert* pV = pF->getVert(iV);
      cout << pV->x() << ":" << pV->y() << ":" << pV->z() << endl;
    }
    assert(pF->getType() == Face::eTriFace);
  }
// The tri, if it already exists, may be in either orientation.
// This could also be checked topologically, and maybe should be.
  Vert *apV[3];
  pF->getAllVertHandles(reinterpret_cast<GRUMMP_Entity**>(apV));

// If this orientation is negative, then surely it means that the
// bottom face is oriented out.
  int iOrient = checkOrient3D(apV[0], apV[1], apV[2], pVRef);
  assert(iOrient != 0);

  if (iOrient == 1) {
    assert(!pF->getRightCell()->isValid());
    pF->setRightCell(pC);
  }
  else {
    assert(!pF->getLeftCell()->isValid());
    pF->setLeftCell(pC);
  }

// This failed to work in some cases.
//   // If the verts are in cyclic order for the tri face, then pC is the
//   // right cell, otherwise the left cell.
//   if (dynamic_cast<TriFace*>(pF)->qRightHanded(pVA, pVB, pVC))
//     { pF->vSetRightCell(pC); }
//   else
//     { pF->vSetLeftCell(pC); }
}

FlakeCell*
VolMesh::createFlakeCell(bool& qExistedAlready, Vert* const pV0,
			 Vert * const pV1, Vert * const pV2, Vert * const pV3,
			 const int iReg)
{
// A flake is a cell that consists of four vertices, arranged so that there
// is one quad face 0123 and two tri faces (012, 023).  These cells are
// temporary glue when merging tets to form more complex elements, so they
// should all be deleted before the merging process is complete.
  bool qExist0123, qExist012, qExist023;
  Face* pF0123 = createFace(qExist0123, pV0, pV1, pV2, pV3);
  Face* pF012 = createFace(qExist012, pV0, pV1, pV2);
  Face* pF023 = createFace(qExist023, pV0, pV2, pV3);
  assert(pF0123->isValid());
  assert(pF012->isValid());
  assert(pF023->isValid());

// It's possible the cell already exists, so check it.
  Cell *pCCand = findCommonCell(pF0123, pF012);
  if (pCCand->isValid()) {
    if (findCommonCell(pF0123, pF023) == pCCand) {
      FlakeCell* pFC = dynamic_cast<FlakeCell*>(pCCand);
      assert(pFC); // No way it can correctly be anything else
      qExistedAlready = true;
      return pFC;
    }
  }

  qExistedAlready = false;
  FlakeCell *pFC = dynamic_cast<FlakeCell*>(getNewCell(3, 4));
  assert(pFC);
  pFC->setRegion(iReg);
// Be sure to assign these with canonical ordering.
  pFC->assignFaces(pF0123, pF012, pF023);

// There's no reliable way to get geometric orientation for these faces
// (especially the quad), so do it strictly topologically.

  pF0123->setLeftCell(pFC);
  pF012->addCell(pFC);
  pF023->addCell(pFC);

#ifndef NDEBUG
  assert(pFC->hasVert(pV0));
  assert(pFC->hasVert(pV1));
  assert(pFC->hasVert(pV2));
  assert(pFC->hasVert(pV3));

  assert(pF012->getLeftCell() == pFC || pF012->getRightCell() == pFC);
  assert(pF023->getLeftCell() == pFC || pF023->getRightCell() == pFC);
  assert(pF0123->getLeftCell() == pFC || pF0123->getRightCell() == pFC);
#endif
  return pFC;
}

PyrCell*
VolMesh::createPyrCell(bool& qExistedAlready, Vert * const pV0,
		       Vert * const pV1, Vert * const pV2, Vert * const pV3,
		       Vert * const pV4, const int iReg)
// This one creates pyramids from verts via faces.
{
  bool qExist014, qExist124, qExist234, qExist304, qExist0123;
  Face *pF014 = createFace(qExist014, pV0, pV1, pV4);
  Face *pF124 = createFace(qExist124, pV1, pV2, pV4);
  Face *pF234 = createFace(qExist234, pV2, pV3, pV4);
  Face *pF304 = createFace(qExist304, pV3, pV0, pV4);
  Face *pF0123 = createFace(qExist0123, pV0, pV1, pV2, pV3);
  assert(pF014->isValid());
  assert(pF124->isValid());
  assert(pF234->isValid());
  assert(pF304->isValid());
  assert(pF0123->isValid());

// The pyramid looks like this:
//
//                  4
//                 /|\               .
//                / | \              .
//               / |\  \             .
//              /  / |  \            .
//             /  |  |   \           .
//            /   /  \    \          .
//           /   |_---3-_  \         .
//          0--- /       ---2        .
//           \_ |    ---
//             \1---
//
// Okay, so it's a crappy pyramid; so what?  Faces are:
//  0:  014
//  1:  124
//  2:  234
//  3:  304
//  4:  0123
  return createPyrCell(qExistedAlready, pF014, pF124, pF234, pF304, pF0123,
		       iReg, pV0, pV1, pV2, pV3, pV4);
}

PyrCell*
VolMesh::createPyrCell(bool& qExistedAlready, Face* const pF014,
		       Face* const pF124, Face* const pF234, Face* const pF304,
		       Face* const pF0123, const int iReg, Vert* pV0, Vert* pV1,
		       Vert* pV2, Vert* pV3, Vert* pV4)
{
  qExistedAlready = false;

// This check is a bit awkwardly written, but it tries to minimize the
// number of calls to pCCommonCell that are actually made.

// Cell might exist already
  Cell *pCCom04 = findCommonCell(pF014, pF304);
  if (pCCom04->isValid()) {
    Cell *pCCom24 = findCommonCell(pF124, pF234);
    if (pCCom24->isValid() && pCCom04 == pCCom24) {
      Cell *pCCom01 = findCommonCell(pF014, pF0123);
      if (pCCom01->isValid() && pCCom04 == pCCom01) {
	PyrCell *pPC = dynamic_cast<PyrCell*>(pCCom01);
	if (pPC->isValid()) {
	  qExistedAlready = true;
	  return pPC;
	}
      }
    }
  }

  if (pV0->isValid()) {
    assert(
	pV1->isValid() && pV2->isValid() && pV3->isValid() && pV4->isValid());
  }
  else {
    pV0 = findCommonVert(pF014, pF304, pF0123);
    pV1 = findCommonVert(pF124, pF014, pF0123);
    pV2 = findCommonVert(pF234, pF124, pF0123);
    pV3 = findCommonVert(pF304, pF234, pF0123);
    pV4 = findCommonVert(pF014, pF124, pF234);
  }

  PyrCell *pPC = dynamic_cast<PyrCell*>(getNewCell(5, 5));
  assert(pPC);
  pPC->setRegion(iReg);
// Be sure to assign these with canonical ordering.
  pPC->assignFaces(pF014, pF124, pF234, pF304, pF0123);

  vQuadSetCellWithOrientation(pPC, pF0123, pV4);

  vTriSetCellWithOrientation(pPC, pF014, pV2);
  vTriSetCellWithOrientation(pPC, pF124, pV3);
  vTriSetCellWithOrientation(pPC, pF234, pV0);
  vTriSetCellWithOrientation(pPC, pF304, pV1);

#ifndef NDEBUG
  assert(pPC->hasVert(pV0));
  assert(pPC->hasVert(pV1));
  assert(pPC->hasVert(pV2));
  assert(pPC->hasVert(pV3));
  assert(pPC->hasVert(pV4));

// The orientation checks below can give false positives for warped quads.
//  // Now confirm that we've got orientation right
//  Vert *apV[] =
//    { pPC->getVert(0), pPC->getVert(1), pPC->getVert(2), pPC->getVert(3),
//	pPC->getVert(4) };
//  assert(checkOrient3D(apV[0], apV[1], apV[2], apV[4]) == 1);
//  assert(checkOrient3D(apV[0], apV[2], apV[3], apV[4]) == 1);
//  assert(checkOrient3D(apV[0], apV[4], apV[1], apV[2]) == 1);
//  assert(checkOrient3D(apV[1], apV[4], apV[2], apV[3]) == 1);
//  assert(checkOrient3D(apV[2], apV[4], apV[3], apV[0]) == 1);
//  assert(checkOrient3D(apV[3], apV[4], apV[0], apV[1]) == 1);
// Use a volume check for cell validity, since orientation can be wrong.
  assert(iFuzzyComp(pPC->calcSize(), 0) >= 0);

  assert(pF014->getLeftCell() == pPC || pF014->getRightCell() == pPC);
  assert(pF124->getLeftCell() == pPC || pF124->getRightCell() == pPC);
  assert(pF234->getLeftCell() == pPC || pF234->getRightCell() == pPC);
  assert(pF304->getLeftCell() == pPC || pF304->getRightCell() == pPC);
  assert(pF0123->getLeftCell() == pPC || pF0123->getRightCell() == pPC);
#endif

  createCellEvent(pPC);
  return (pPC);
}

PrismCell*
VolMesh::createPrismCell(bool &qExistedAlready, Vert * const pV0,
			 Vert * const pV1, Vert * const pV2, Vert * const pV3,
			 Vert * const pV4, Vert * const pV5, const int iReg)
// This one creates prisms from verts via faces.
{
// Verts 0, 1, 2 form a ring around the bottom.
// Verts 3, 4, 5 form a ring around the top.

  bool qExist0143, qExist1254, qExist2035, qExist012, qExist543;
  Face *pF0143 = createFace(qExist0143, pV0, pV1, pV4, pV3);
  Face *pF1254 = createFace(qExist1254, pV1, pV2, pV5, pV4);
  Face *pF2035 = createFace(qExist2035, pV2, pV0, pV3, pV5);
  Face *pF012 = createFace(qExist012, pV0, pV1, pV2);
  Face *pF543 = createFace(qExist543, pV5, pV4, pV3);

  return createPrismCell(qExistedAlready, pF0143, pF1254, pF2035, pF012, pF543,
			 iReg, pV0, pV1, pV2, pV3, pV4, pV5);
}

PrismCell*
VolMesh::createPrismCell(bool &qExistedAlready, Face * const pF0143,
			 Face * const pF1254, Face * const pF2035,
			 Face * const pF012, Face * const pF543, const int iReg,
			 Vert * pV0, Vert * pV1, Vert * pV2, Vert * pV3,
			 Vert * pV4, Vert * pV5)
{
  qExistedAlready = false;

// Cell might exist already
  Cell *pCCom03 = findCommonCell(pF0143, pF2035);
  if (pCCom03->isValid()) {
    Cell *pCCom12 = findCommonCell(pF1254, pF012);
    if (pCCom12->isValid() && pCCom03 == pCCom12) {
      Cell *pCCom14 = findCommonCell(pF1254, pF543);
      if (pCCom14->isValid() && pCCom03 == pCCom14) {
	PrismCell *pPC = dynamic_cast<PrismCell*>(pCCom03);
	if (pPC->isValid()) {
	  qExistedAlready = true;
	  return pPC;
	}
      }
    }
  }

// Creating a new one.
  if (pV0->isValid()) {
    assert(
	pV1->isValid() && pV2->isValid() && pV3->isValid() && pV4->isValid()
	    && pV5->isValid());
  }
  else {
    pV0 = findCommonVert(pF0143, pF2035, pF012);
    pV1 = findCommonVert(pF0143, pF1254, pF012);
    pV2 = findCommonVert(pF2035, pF1254, pF012);
    pV3 = findCommonVert(pF0143, pF2035, pF543);
    pV4 = findCommonVert(pF0143, pF1254, pF543);
    pV5 = findCommonVert(pF2035, pF1254, pF543);
  }
  PrismCell *pPC = dynamic_cast<PrismCell*>(getNewCell(5, 6));
  assert(pPC);
  pPC->setRegion(iReg);
// Be sure to assign these with canonical ordering.
  pPC->assignFaces(pF0143, pF1254, pF2035, pF012, pF543);

  vQuadSetCellWithOrientation(pPC, pF0143, pV2);
  vQuadSetCellWithOrientation(pPC, pF1254, pV0);
  vQuadSetCellWithOrientation(pPC, pF2035, pV1);

  vTriSetCellWithOrientation(pPC, pF012, pV3);
  vTriSetCellWithOrientation(pPC, pF543, pV0);

  assert(pF0143->getLeftCell() == pPC || pF0143->getRightCell() == pPC);
  assert(pF1254->getLeftCell() == pPC || pF1254->getRightCell() == pPC);
  assert(pF2035->getLeftCell() == pPC || pF2035->getRightCell() == pPC);
  assert(pF012->getLeftCell() == pPC || pF012->getRightCell() == pPC);
  assert(pF543->getLeftCell() == pPC || pF543->getRightCell() == pPC);

  assert(pPC->hasVert(pV0));
  assert(pPC->hasVert(pV1));
  assert(pPC->hasVert(pV2));
  assert(pPC->hasVert(pV3));
  assert(pPC->hasVert(pV4));
  assert(pPC->hasVert(pV5));

  createCellEvent(pPC);
  return pPC;
}

HexCell*
VolMesh::createHexCell(bool &qExistedAlready, Vert * const pV0,
		       Vert * const pV1, Vert * const pV2, Vert * const pV3,
		       Vert * const pV4, Vert * const pV5, Vert * const pV6,
		       Vert * const pV7, const int iReg)
// This one creates prisms from verts via faces.
{
// Verts 0, 1, 2, 3 form a ring around the bottom.
// Verts 4, 5, 6, 7 form a ring around the top.

  bool qExist0123, qExist7654, qExist0154, qExist1265, qExist2376, qExist3047;
  Face *pF0154 = createFace(qExist0154, pV0, pV1, pV5, pV4);
  Face *pF1265 = createFace(qExist1265, pV1, pV2, pV6, pV5);
  Face *pF2376 = createFace(qExist2376, pV2, pV3, pV7, pV6);
  Face *pF3047 = createFace(qExist3047, pV3, pV0, pV4, pV7);
  Face *pF0123 = createFace(qExist0123, pV0, pV1, pV2, pV3);
  Face *pF7654 = createFace(qExist7654, pV7, pV6, pV5, pV4);

  return createHexCell(qExistedAlready, pF0154, pF1265, pF2376, pF3047, pF0123,
		       pF7654, iReg, pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7);
}

HexCell*
VolMesh::createHexCell(bool &qExistedAlready, Face* const pF0154,
		       Face* const pF1265, Face* const pF2376,
		       Face* const pF3047, Face* const pF0123,
		       Face* const pF7654, const int iRegion, Vert * pV0,
		       Vert * pV1, Vert * pV2, Vert * pV3, Vert * pV4,
		       Vert * pV5, Vert * pV6, Vert * pV7)
{
// Cell might exist already.
  Cell *pCCom01 = findCommonCell(pF0123, pF0154);
  if (pCCom01->isValid()) {
    Cell *pCCom26 = findCommonCell(pF1265, pF2376);
    if (pCCom26->isValid() && pCCom01 == pCCom26) {
      Cell *pCCom47 = findCommonCell(pF7654, pF3047);
      if (pCCom47->isValid() && pCCom01 == pCCom47) {
	HexCell *pHC = dynamic_cast<HexCell*>(pCCom01);
	if (pHC->isValid()) {
	  qExistedAlready = true;
	  return pHC;
	}
      }
    }
  }

  qExistedAlready = false;

  if (pV0->isValid()) {
    assert(
	pV1->isValid() && pV2->isValid() && pV3->isValid() && pV4->isValid()
	    && pV5->isValid() && pV6->isValid() && pV7->isValid());
  }
  else {
    pV0 = findCommonVert(pF0154, pF3047, pF0123);
    pV1 = findCommonVert(pF1265, pF0154, pF0123);
    pV2 = findCommonVert(pF2376, pF1265, pF0123);
    pV3 = findCommonVert(pF3047, pF2376, pF0123);
    pV4 = findCommonVert(pF0154, pF3047, pF7654);
    pV5 = findCommonVert(pF1265, pF0154, pF7654);
    pV6 = findCommonVert(pF2376, pF1265, pF7654);
    pV7 = findCommonVert(pF3047, pF2376, pF7654);
  }

  HexCell *pHC = dynamic_cast<HexCell*>(getNewCell(6, 8));
  assert(pHC);
  pHC->setRegion(iRegion);
// Be sure to assign these with canonical ordering.
  pHC->assignFaces(pF0123, pF0154, pF1265, pF2376, pF3047, pF7654);

  vQuadSetCellWithOrientation(pHC, pF0123, pV4);
  vQuadSetCellWithOrientation(pHC, pF7654, pV0);
  vQuadSetCellWithOrientation(pHC, pF0154, pV2);
  vQuadSetCellWithOrientation(pHC, pF1265, pV3);
  vQuadSetCellWithOrientation(pHC, pF2376, pV4);
  vQuadSetCellWithOrientation(pHC, pF3047, pV1);

  assert(pF0123->hasCell(pHC));
  assert(pF7654->hasCell(pHC));
  assert(pF0154->hasCell(pHC));
  assert(pF1265->hasCell(pHC));
  assert(pF2376->hasCell(pHC));
  assert(pF3047->hasCell(pHC));

  assert(pHC->hasVert(pV0));
  assert(pHC->hasVert(pV1));
  assert(pHC->hasVert(pV2));
  assert(pHC->hasVert(pV3));
  assert(pHC->hasVert(pV4));
  assert(pHC->hasVert(pV5));
  assert(pHC->hasVert(pV6));
  assert(pHC->hasVert(pV7));

  createCellEvent(pHC);
  return pHC;
}

bool
VolMesh::deleteCell(Cell * const pC)
{
  if (omp_in_parallel() == 0) {
    if (pC->isDeleted())
      return true;
    isOKToSend(false);
    // Free up connectivity for faces.
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      Face *pF = pC->getFace(i);
      if (pF->isValid())
	pF->removeCell(pC);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
	deleteFace(pF);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTet:
	m_ECTet.deleteEntry(static_cast<TetCell*>(pC));
	break;
      case CellSkel::ePyr:
	m_ECPyr.deleteEntry(static_cast<PyrCell*>(pC));
	break;
      case CellSkel::ePrism:
	m_ECPrism.deleteEntry(static_cast<PrismCell*>(pC));
	break;
      case CellSkel::eHex:
	m_ECHex.deleteEntry(static_cast<HexCell*>(pC));
	break;
      case CellSkel::eFlake:
	m_ECFlake.deleteEntry(static_cast<FlakeCell*>(pC));
	break;
      default:
	assert(0);
	break;
      }
    isOKToSend(true);
    deleteCellEvent(pC);
    return true;
  }
  else {
    int ID = omp_get_thread_num();
    if (pC->isDeleted())
      return true;
    isOKToSend(false);
    // Free up connectivity for faces.
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      Face *pF = pC->getFace(i);
      if (pF->isValid())
	pF->removeCell(pC);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
	deleteFaceParallel(pF, ID);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTet:
	m_ECTet.deleteEntry(static_cast<TetCell*>(pC), ID);
	break;
      case CellSkel::ePyr:
	m_ECPyr.deleteEntry(static_cast<PyrCell*>(pC), ID);
	break;
      case CellSkel::ePrism:
	m_ECPrism.deleteEntry(static_cast<PrismCell*>(pC), ID);
	break;
      case CellSkel::eHex:
	m_ECHex.deleteEntry(static_cast<HexCell*>(pC), ID);
	break;
      default:
	assert(0);
	break;
      }
    isOKToSend(true);
//		#pragma omp critical(Event3D)
//		  {
    deleteCellEvent(pC);
//		  }
    return true;
  }
}

bool
VolMesh::deleteCellParallel(Cell * const pC, int ID)
{
  if (NUM_PROCS > 1) {
    assert(omp_in_parallel() != 0);
    assert(ID == omp_get_thread_num());
  }
  if (pC->isDeleted())
    return true;
  isOKToSend(false);
// Free up connectivity for faces.
  for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
    Face *pF = pC->getFace(i);
    if (pF->isValid())
      pF->removeCell(pC);
    if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
      deleteFaceParallel(pF, ID);
    }
  }
// Under the old data structures, just mark it.
  pC->markAsDeleted();
// With the new (2010) containers, also register the deletion with the
// container.
  switch (pC->getType())
    {
    case CellSkel::eTet:
      m_ECTet.deleteEntry(static_cast<TetCell*>(pC), ID);
      break;
    case CellSkel::ePyr:
      m_ECPyr.deleteEntry(static_cast<PyrCell*>(pC), ID);
      break;
    case CellSkel::ePrism:
      m_ECPrism.deleteEntry(static_cast<PrismCell*>(pC), ID);
      break;
    case CellSkel::eHex:
      m_ECHex.deleteEntry(static_cast<HexCell*>(pC), ID);
      break;
    default:
      assert(0);
      break;
    }
  isOKToSend(true);
//	#pragma omp critical(Event3D)
//  {
  deleteCellEvent(pC);
//  }
  return true;
}

// This function is not technically necessary, under the old data
// structures; deleteCell would work fine, too.
bool
VolMesh::deleteBFace(BFace * const pBF)
{
  if (omp_in_parallel() != 0) {
    int ID = omp_get_thread_num();
    return deleteBFaceParallel(pBF, ID);
  }
  if (pBF->isDeleted())
    return true;
// Under the old data structures, just mark it.
  Face *pF = pBF->getFace();
  if (pF->isValid()) {
    pF->removeCell(pBF);
    if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
      deleteFace(pF);
    }
  }

  if (pBF->getNumFaces() == 2 && pBF->getFace(1)->isValid()) {
    pF = pBF->getFace(1);
    pF->removeCell(pBF);
    if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
      deleteFace(pF);
    }
  }

  pBF->markAsDeleted();
// With the new (2010) containers, also register the container
  switch (pBF->getType())
    {
    case CellSkel::eTriBFace:
      m_ECTriBF.deleteEntry(static_cast<TriBFace*>(pBF));
      break;
    case CellSkel::eQuadBFace:
      m_ECQuadBF.deleteEntry(static_cast<QuadBFace*>(pBF));
      break;
    case CellSkel::eIntTriBFace:
      m_ECIntTriBF.deleteEntry(static_cast<IntTriBFace*>(pBF));
      break;
    case CellSkel::eIntQuadBFace:
      m_ECIntQuadBF.deleteEntry(static_cast<IntQuadBFace*>(pBF));
      break;
    default:
      assert(0);
      break;
    }
  deleteBFaceEvent(pBF);
  return true;
}

bool
VolMesh::deleteBFaceParallel(BFace * const pBF, int ID)
{
  if (pBF->isDeleted())
    return true;
// Under the old data structures, just mark it.
  Face *pF = pBF->getFace();
  if (pF->isValid()) {
    pF->removeCell(pBF);
    if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
      deleteFaceParallel(pF, ID);
    }
  }

  if (pBF->getNumFaces() == 2 && pBF->getFace(1)->isValid()) {
    pF = pBF->getFace(1);
    pF->removeCell(pBF);
    if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
      deleteFaceParallel(pF, ID);
    }
  }

  pBF->markAsDeleted();
// With the new (2010) containers, also register the container
  switch (pBF->getType())
    {
    case CellSkel::eTriBFace:
      m_ECTriBF.deleteEntry(static_cast<TriBFace*>(pBF), ID);
      break;
    case CellSkel::eQuadBFace:
      m_ECQuadBF.deleteEntry(static_cast<QuadBFace*>(pBF), ID);
      break;
    case CellSkel::eIntTriBFace:
      m_ECIntTriBF.deleteEntry(static_cast<IntTriBFace*>(pBF), ID);
      break;
    case CellSkel::eIntQuadBFace:
      m_ECIntQuadBF.deleteEntry(static_cast<IntQuadBFace*>(pBF), ID);
      break;
    default:
      assert(0);
      break;
    }
//#pragma omp critical(Event3D)
//	{
  deleteBFaceEvent(pBF);
//	}
  return true;
}

// Like ITAPS, faces that are in use can't be deleted.
bool
VolMesh::deleteFace(Face * const pF)
{
  if (pF->isDeleted())
    return true;
// Make sure that verts whose hint is going away realize it...
  for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
    Vert *pV = pF->getVert(iV);
    if (pV->isDeleted())
      continue;
//		if(pV->hasFace(pF))
    pV->removeFace(pF);
#ifndef OMIT_VERTEX_HINTS
    if (pV->getHintFace() == pF) {
      pV->updateHintFace();
    }
#endif
  }

// Under the old data structures, just mark it.
  if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
    pF->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pF->getType())
      {
      case Face::eTriFace:
	m_ECTriF.deleteEntry(static_cast<TriFace*>(pF));
	break;
      case Face::eQuadFace:
	m_ECQuadF.deleteEntry(static_cast<QuadFace*>(pF));
	break;
      default:
	assert(0);
	break;
      }
    deleteFaceEvent(pF);
    return true;
  }
  else {
    return false;
  }
}

bool
VolMesh::deleteFaceParallel(Face * const pF, int ID)
{
  if (pF->isDeleted())
    return true;
// Make sure that verts whose hint is going away realize it...
  for (int iV = pF->getNumVerts() - 1; iV >= 0; iV--) {
    Vert *pV = pF->getVert(iV);
    if (pV->isDeleted())
      continue;
    pV->removeFace(pF);
#ifndef OMIT_VERTEX_HINTS
    if (pV->getHintFace() == pF) {
      pV->updateHintFace();
    }
#endif
  }

// Under the old data structures, just mark it.
  if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
    pF->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pF->getType())
      {
      case Face::eTriFace:
	m_ECTriF.deleteEntry(static_cast<TriFace*>(pF), ID);
	break;
      case Face::eQuadFace:
	m_ECQuadF.deleteEntry(static_cast<QuadFace*>(pF), ID);
	break;
      default:
	assert(0);
	break;
      }
//#pragma omp critical(Event3D)
//{
    deleteFaceEvent(pF);
//}
    return true;
  }
  else {
    return false;
  }
}

Vert*
VolMesh::createVert(const double dX, const double dY, const double dZ)
{
  double adData[] =
    { dX, dY, dZ };
  return createVert(adData);
}

Vert*
VolMesh::createVert(const double adCoords[])
{
  int ID = 0;
  if (omp_in_parallel() != 0)
    ID = omp_get_thread_num();
  Vert *pV = getNewVert(ID);
  pV->clearHintFace();
  pV->clearFaceConnectivity();
  pV->setDefaultFlags();
  pV->setCoords(3, adCoords);
//#pragma omp critical(Event3D)
//	{
  createVertEvent(pV);
//	}
  return pV;
}

Vert*
VolMesh::createVert(const Vert* pVIn)
{
  Vert *pV = getNewVert();

// Can't use Vert::operator= here because that changes the
// face-vert connectivity.
  *pV = *pVIn;
  pV->clearHintFace();
  pV->clearFaceConnectivity();
  createVertEvent(pV);
  return pV;
}

Vert*
VolMesh::duplicateVert(const Vert* pVIn)
{
  Vert *pV = getNewVert();

// Can't use Vert::operator= here because that changes the
// face-vert connectivity.
  pV->setCoords(3, pVIn->getCoords());

  pV->copyAllFlagsFrom(pVIn);
//	pV->copyInfluence(pVIn);
// TODO Technically, should copy metric here.  But the metric may or may not
// stay part of the Vert class.
  pV->setLengthScale(pVIn->getLengthScale());
  pV->setParentEntity(pVIn->getParentEntity());
  pV->setParentTopology(pVIn->getParentTopology());

  if (pVIn->isMetricUsed()) {
    pV->setMetric(pVIn->getMetricArray());
  }

  pV->clearHintFace();
  pV->clearFaceConnectivity();
  createVertEvent(pV);
  return pV;
}

bool
VolMesh::deleteVert(Vert * const pV)
// Again, verts that are in use can't be deleted without forcing that.
{
  if (pV->isDeleted())
    return true;
  if (pV->getNumFaces() == 0) {
    // Under the old data structures, just mark it.
    if (m_keepPhantomVerts) {
      pV->markAsPhantom(true);
    }
    else {
      pV->markAsDeleted();
      m_ECVerts.deleteEntry(pV);
    }
    deleteVertEvent(pV);
    return true;
  }
  else {
    return false;
  }
}

std::pair<Subseg*, Subseg*>
VolMesh::createNewSubsegsBySplit(Subseg * oldSubseg,
				 const double splitCoords[3],
				 const double splitParam, Vert* const newVert)
{
  std::pair<Subseg*, Subseg*> newSubsegs(new Subseg(), new Subseg);

  oldSubseg->refinement_split(splitCoords, splitParam, newVert,
			      newSubsegs.first, newSubsegs.second);
  assert(newSubsegs.first->is_not_deleted());
  assert(newSubsegs.second->is_not_deleted());
  assert(oldSubseg->is_deleted());
  oldSubseg->mark_deleted();
  m_SubsegMap->add_subseg(newSubsegs.first);
  m_SubsegMap->add_subseg(newSubsegs.second);
  m_SubsegMap->remove_subseg(oldSubseg);
  return newSubsegs;
}

static Cell*
identifyCellInsidePrism(Face* pFABC, Face* pFDEF)
{
  // pFABC and pFDEF are the faces that are endcaps of a prism (in no
  // particular order).  Find the cell adjacent to pFABC whose other
  // vertex or vertices are on pFDEF.  If there isn't such a cell,
  // return nullptr.

  // Find the cell that's inside the prism and attached to tri ABC.
  Cell* pCCand = pFABC->getLeftCell();
  bool gotIt = true;
  if (pCCand->getNumVerts() == 3) {
    gotIt = false;
  }
  else {
  for (int ii = 0; ii < pCCand->getNumVerts(); ii++) {
    Vert* pVCand = pCCand->getVert(ii);
    if (pFABC->hasVert(pVCand))
      continue;

    gotIt = gotIt && pFDEF->hasVert(pVCand);
  }
  }

  if (!gotIt) {
    gotIt = true;
    pCCand = pFABC->getRightCell();
    if (pCCand->getNumVerts() == 3)
      return nullptr;
    
    for (int ii = 0; ii < pCCand->getNumVerts(); ii++) {
      Vert* pVCand = pCCand->getVert(ii);
      if (pFABC->hasVert(pVCand))
	continue;

      gotIt = gotIt && pFDEF->hasVert(pVCand);
    }
  }
  if (gotIt) {
    return pCCand;
  }
  else {
    return nullptr;
  }
}

static Cell*
identifyFlake(Vert* pVA, Vert* pVB, Vert* pVC, Vert* pVD)
{
  Face *pFABCD = findCommonFace(pVA, pVB, pVC, pVD, true);
  if (!pFABCD) {
    // No quad, so there can't be a flake.
    return nullptr;
  }
  Cell *pCCand = pFABCD->getLeftCell();
  if (pCCand->isValid() && pCCand->getNumFaces() == 3) {
    return pCCand;
  }
  else {
    pCCand = pFABCD->getRightCell();
    if (pCCand->isValid() && pCCand->getNumFaces() == 3) {
      return pCCand;
    }
  }
  return nullptr;
}

void
VolMesh::handleFlakeCell(Cell* pFCABED, Vert* pVA, Vert* pVB, Vert* pVE,
			 Vert* pVD)
{
  if (pFCABED) {
    deleteCell(pFCABED);
  }
  else {
    if (findCommonFace(pVA, pVB, pVE)->isValid()) {
      bool qExists = false;
      FlakeCell* pFC = createFlakeCell(qExists, pVA, pVB, pVE, pVD);
      assert(!qExists);
      assert(pFC);
    }
    else {
      bool qExists = false;
      FlakeCell* pFC = createFlakeCell(qExists, pVB, pVE, pVD, pVA);
      assert(!qExists);
      assert(pFC);
    }
  }
}

PrismCell*
VolMesh::replaceCellsWithPrism(Vert* pVA, Vert* pVB, Vert* pVC, Vert* pVD,
			       Vert* pVE, Vert* pVF)
{
  // ABC is the bottom face, and DEF is the top face.

  // The top and bottom faces had better exist.
  Face *pFABC = findCommonFace(pVA, pVB, pVC, nullptr, true);
  Face *pFDEF = findCommonFace(pVD, pVE, pVF, nullptr, true);
  assert(pFABC);
  assert(pFDEF);

// Find the cell that's inside the prism and attached to tri ABC.
  Cell* cells[3];
  cells[0] = identifyCellInsidePrism(pFABC, pFDEF);
  assert(cells[0]);
  if (cells[0]->getNumVerts() == 6)
    return (dynamic_cast<PrismCell*>(cells[0]));

  cells[1] = identifyCellInsidePrism(pFDEF, pFABC);

  Cell* pFCABED = identifyFlake(pVA, pVB, pVE, pVD);
  Cell* pFCBCFE = identifyFlake(pVB, pVC, pVF, pVE);
  Cell* pFCCADF = identifyFlake(pVC, pVA, pVD, pVF);

  if (cells[0]->getNumVerts() + cells[1]->getNumVerts() == 9) {
    // One of these two cells is a pyramid.  No need to find a third cell.
    assert(0);
    return nullptr;
  }
  else {
    // Need to find the other tet.  It has to share a face with cells[0].
    for (int iF = 0; iF < 4; iF++) {
      Face *pF = cells[0]->getFace(iF);
      TetCell *pTC = dynamic_cast<TetCell*>(pF->getOppositeCell(cells[0]));
      if (!pTC->isValid()) {
	// This tri face could be on the prism perimeter and have a neighboring
	// cell that's a bdry tri, pyramid, prism, or flake.  Granted, some of
	// these are possibilities we wouldn't design for, but in principle
	// they could occur.  In those cases, pTC would be invalid.
	continue;
      }
      Vert *pVOpp = pTC->getOppositeVert(pF);
      if (pFDEF->hasVert(pVOpp)) {
	cells[2] = pTC;
	break;
      }
    }

    // Now delete all the cells that are inside the prism.  That's cells[0..2] and
    // any of the flakes that are non-null.
    deleteCell(cells[0]);
    deleteCell(cells[1]);
    deleteCell(cells[2]);

    handleFlakeCell(pFCABED, pVA, pVB, pVE, pVD);
    handleFlakeCell(pFCBCFE, pVB, pVC, pVF, pVE);
    handleFlakeCell(pFCCADF, pVC, pVA, pVD, pVF);

    bool qExists = false;
    PrismCell *retVal = createPrismCell(qExists, pVA, pVB, pVC, pVD, pVE, pVF);
    return retVal;
  }

}
