#include "GR_config.h"
#include "GR_Classes.h"

#include "GR_BFace.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Geometry.h"
#include "GR_GRCurve.h"
#include "GR_Mesh2D.h"
#include "GR_InsertionManager.h"
#include "GR_Vertex.h"
#include <omp.h>

// Get new entities from the database.

Face*
Mesh2D::getNewFace(const int iNV, int ID)
{
  assert(iNV == 2);
  Face *pF;
//	assert(ID == omp_get_thread_num());
  pF = m_ECEdgeF.getNewEntry(ID);
  pF->resetAllData();
  return pF;
}

BFace*
Mesh2D::getNewBFace(const int iNBV, int ID)
{
  assert(iNBV == 2);
  BFace *pBF = m_ECEdgeBF.getNewEntry(ID);
  pBF->resetAllData();
  return pBF;
}

IntBdryEdge*
Mesh2D::getNewIntBdryEdge(const int iNBV, int ID)
{
  assert(iNBV == 2);
  IntBdryEdge *pIBE = m_ECIntEdgeBF.getNewEntry(ID);
  pIBE->resetAllData();
  return pIBE;
}

//GRUMMP::PartBdryEdge* Mesh2D::
//getNewPartBdryEdge (const int iNBV, int ID)
//{
//    assert(iNBV == 2);
//    GRUMMP::PartBdryEdge *pPBE = m_ECPartBdryEdge.getNewEntry(ID);
//    pPBE->resetAllData();
//    return pPBE;
//}

Cell*
Mesh2D::getNewCell(const int iNF, int ID)
{
//	assert(ID == omp_get_thread_num());
  TriCell *pTriC;
  QuadCell *pQuadC;

  assert(iNF == 3 || (iNF == 4 && !m_simplexMesh));
  if (iNF == 3) {
    pTriC = m_ECTri.getNewEntry(ID);
    pTriC->resetAllData();
    pTriC->setType(Cell::eTriCell);
    return pTriC;
  }
  else {
    pQuadC = m_ECQuad.getNewEntry(ID);
    pQuadC->resetAllData();
    pQuadC->setType(Cell::eQuadCell);
    return pQuadC;
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

// Face* Mesh2D::createFace(Vert * const apV[], const int iNVerts)
// {
//   assert(iNVerts == 2);
//   Vert *pV0 = apV[0];
//   Vert *pV1 = apV[1];
//   assert(pV0->isValid());
//   assert(pV1->isValid());
//   return (createFace(pV0, pV1));
// }

Face*
Mesh2D::createFace(bool& qAlreadyExisted, Vert * const pV0, Vert * const pV1,
		   const bool qForceDuplicate, int ID)
{
//	if (omp_in_parallel() != 0) assert(omp_get_num_threads() == NUM_PROCS);
//	assert(ID == omp_get_thread_num());
  assert(pV0->isValid());
  assert(pV1->isValid());
  if (pV0 == pV1) {
    logMessage(1, "%f:%f\n", pV0->x(), pV0->y());
    assert(pV0 != pV1);
  }
  // Check for prior existence.
  Face *pF;
  pF = findCommonFace(pV0, pV1, qForceDuplicate);
  qAlreadyExisted = (pF != pFInvalidFace);
  // If we're forcing a duplicate, we'd better have a copy already!
  assert((!qForceDuplicate) || (qAlreadyExisted));
  if (qForceDuplicate == qAlreadyExisted) {
    // Only do this if we've got a correct duplication case or a
    // correct creation case.

    pF = getNewFace(2, ID);
    pF->setVerts(pV0, pV1);
    pF->setLeftCell(NULL);
    pF->setRightCell(NULL);
//#pragma omp critical(Event)
//		{
    createFaceEvent(pF);
//		}
  }
//	if (qAlreadyExisted != true){
//		pF->layerindex = minlayerindex + 1;
//	}

#ifndef OMIT_VERTEX_HINTS
  // Make sure verts always have a valid hint face...

  if (!pV0->getHintFace()->isValid())
    pV0->setHintFace(pF);
  if (!pV1->getHintFace()->isValid())
    pV1->setHintFace(pF);

#endif

  assert(pF->doFullCheck());

  return pF;
}

Face*
Mesh2D::createFace(int minlayerindex, bool& qAlreadyExisted, Vert * const pV0,
		   Vert * const pV1, const bool qForceDuplicate, int ID)
{
//	if (omp_in_parallel() != 0) assert(omp_get_num_threads() == NUM_PROCS);
//	assert(ID == omp_get_thread_num());
  assert(pV0->isValid());
  assert(pV1->isValid());
  assert(pV0 != pV1);
  // Check for prior existence.
  Face *pF;
  pF = findCommonFace(pV0, pV1, qForceDuplicate);
  qAlreadyExisted = (pF != pFInvalidFace);
  // If we're forcing a duplicate, we'd better have a copy already!
  assert((!qForceDuplicate) || (qAlreadyExisted));
  if (qForceDuplicate == qAlreadyExisted) {
    // Only do this if we've got a correct duplication case or a
    // correct creation case.

    pF = getNewFace(2, ID);
    pF->setVerts(pV0, pV1);
    pF->setLeftCell(NULL);
    pF->setRightCell(NULL);
//#pragma omp critical(Event)
//		{
    createFaceEvent(pF);
//		}

    if (qAlreadyExisted != true) {
      pF->layerindex = minlayerindex + 1;
    }
  }

#ifndef OMIT_VERTEX_HINTS
  // Make sure verts always have a valid hint face...

  if (!pV0->getHintFace()->isValid())
    pV0->setHintFace(pF);
  if (!pV1->getHintFace()->isValid())
    pV1->setHintFace(pF);

#endif

  assert(pF->doFullCheck());

  return pF;
}

BFace*
Mesh2D::createBFace(Face * const pF, BFace * const pBFOld)
{
  int ID = 0;
//	if (omp_in_parallel() != 0)
//	{
  ID = omp_get_thread_num();
//		assert(NUM_PROCS == omp_get_num_threads());
//	}
  // Faces must already exist
  assert(pF->isValid());
  if (pF->getNumCells() == 1) {
    BFace *pBF = getNewBFace(2, ID);

    BdryEdge * pBE = dynamic_cast<BdryEdge *>(pBF);
    BdryEdge * pBEOld = dynamic_cast<BdryEdge *>(pBFOld);
    pBF->assignFaces(pF);
    pF->addCell(pBF);
    pF->setFaceLoc(Face::eBdryFace);

    assert(pBF->doFullCheck());
    assert(pBF->getFace() == pF);
    assert(pF->getLeftCell() == pBF || pF->getRightCell() == pBF);

#ifndef OMIT_VERTEX_HINTS
    // Make sure that bdry verts have bdry hints
    pF->getVert(0)->setHintFace(pF);
    pF->getVert(1)->setHintFace(pF);
#endif

    // This is a really stinky way to copy bdry data, but that's what
    // works currently in 2D.
    if (pBFOld->isValid()) {
      if (!(pBEOld->doesCurveExist())) {
	dynamic_cast<BdryEdgeBase *>(pBF)->setPatch(
	    reinterpret_cast<BdryPatch*>(0x1)); // NULL would trigger an assertion...
      }
      else {
	GRCurve * curve = pBEOld->getCurve();

	pBE->setCurve(curve);
	pBE->setGeometry(pBEOld->getGeometry());
	CubitVector split_coord0 = CubitVector(pBF->getVert(0)->x(),
					       pBF->getVert(0)->y(), 0);
	CubitVector split_coord1 = CubitVector(pBF->getVert(1)->x(),
					       pBF->getVert(1)->y(), 0);

	if (pBE->isForward()) {
	  pBE->setVert0Param(curve->u_from_position(split_coord0));
	  pBE->setVert1Param(curve->u_from_position(split_coord1));
	}
	else {
	  pBE->setVert1Param(curve->u_from_position(split_coord0));
	  pBE->setVert0Param(curve->u_from_position(split_coord1));
	  pF->interchangeCellsAndVerts();

	}
	pBE->setBdryCondition(pBEOld->getBdryCondition());

      }
    }
    createBFaceEvent(pBF);
    return pBF;
  } // Done with the true bdry case.
  else if (pF->getNumCells() == 2) {
    // This is an internal bdry face.  The face must have existing cells
    // on both sides so we can connect it up properly.
    assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());

    bool qExists;
    Face *pFDup = createFace(qExists, pF->getVert(0), pF->getVert(1), true, ID);
    assert(qExists);
    IntBdryEdge * pBEOld = dynamic_cast<IntBdryEdge *>(pBFOld);
    IntBdryEdge * pBF = getNewIntBdryEdge(2, ID);
    IntBdryEdge * pBE = dynamic_cast<IntBdryEdge *>(pBF);
    pBF->assignFaces(pF, pFDup);

    Cell *pCR = pF->getRightCell();
    pFDup->setLeftCell(pBF);
    pFDup->setRightCell(pCR);
    pCR->replaceFace(pF, pFDup);
    pF->setRightCell(pBF);

    pF->setFaceLoc(Face::eBdryTwoSide);
    pFDup->setFaceLoc(Face::eBdryTwoSide);
    assert(pBF->doFullCheck());
    assert(pCR->doFullCheck());
    assert(pF->getLeftCell()->doFullCheck());
    assert(pFDup->getRightCell()->doFullCheck());
#ifndef OMIT_VERTEX_HINTS
    // Make sure that bdry verts have bdry hints
    pF->getVert(0)->setHintFace(pF);
    pF->getVert(1)->setHintFace(pF);
#endif

    if (pBFOld->isValid()) {
      assert(pBFOld->getType() == Cell::eIntBdryEdge);

      if (!(pBEOld->doesCurveExist())) {
	dynamic_cast<BdryEdgeBase *>(pBF)->setPatch(
	    reinterpret_cast<BdryPatch*>(0x1)); // NULL would trigger an assertion...
      }
      else {
	GRCurve * curve = pBEOld->getCurve();

	pBE->setCurve(curve);
	pBE->setGeometry(pBEOld->getGeometry());
	CubitVector split_coord0 = CubitVector(pBF->getVert(0)->x(),
					       pBF->getVert(0)->y(), 0);
	CubitVector split_coord1 = CubitVector(pBF->getVert(1)->x(),
					       pBF->getVert(1)->y(), 0);

	double param0 = curve->u_from_position(split_coord0);
	double param1 = curve->u_from_position(split_coord1);

	if (pBE->isForward()) {
	  pBE->setVert0Param(param0);
	  pBE->setVert1Param(param1);
	}
	else {
	  pBE->setVert1Param(param1);
	  pBE->setVert0Param(param0);
	}
	pBE->setBdryCondition(pBEOld->getBdryCondition());
      }
    }
    return pBF;
  }
  else {
    assert(pF->getNumCells() == 0);
    logMessage(MSG_DETAIL, "Face with no adjacent cells...\n");
    logMessage(MSG_DETAIL, "%f: %f\n", pF->getVert(0)->x(),
	       pF->getVert(0)->y());
    logMessage(MSG_DETAIL, "%f: %f\n", pF->getVert(1)->x(),
	       pF->getVert(1)->y());
    assert(0);
    return nullptr;
  }
}
//TODO  Reza added this. I don't know when and where I run into this function instead of the one above.
//      If so, I'll need to modify it to handle parallel too.
BFace*
Mesh2D::createBFace(Face * const pF, const int i)
{
  // Faces must already exist
  assert(pF->isValid());
  // Must have exactly one open connectivity slot for the face!
  assert(XOR(pF->getLeftCell()->isValid(), pF->getRightCell()->isValid()));
  BFace *pBF = getNewBFace();

  pBF->assignFaces(pF);
  pF->addCell(pBF);
  pF->setFaceLoc(Face::eBdryFace);
  assert(pBF->doFullCheck());
  assert(pBF->getFace() == pF);
  assert(pF->getLeftCell() == pBF || pF->getRightCell() == pBF);

#ifndef OMIT_VERTEX_HINTS
  // Make sure that bdry verts have bdry hints
  pF->getVert(0)->setHintFace(pF);
  pF->getVert(1)->setHintFace(pF);
#endif

  dynamic_cast<BdryEdge*>(pBF)->setBdryCondition(abs(i));

  pBF->getVert(0)->setType(Vert::eBdry);
  pBF->getVert(1)->setType(Vert::eBdry);
  return pBF;
}

IntBdryEdge*
Mesh2D::createIntBdryEdge(Face* const pFA, Face* const pFB,
			  BFace* const pIBEOld)
{
  IntBdryEdge * pBEOld = dynamic_cast<IntBdryEdge *>(pIBEOld);
  // Faces must already exist
  assert(pFA->isValid());
  assert(pFB->isValid());
  // Must have exactly one open connectivity slot for the face!
  assert(XOR(pFA->getLeftCell()->isValid(), pFA->getRightCell()->isValid()));
  assert(XOR(pFB->getLeftCell()->isValid(), pFB->getRightCell()->isValid()));

  IntBdryEdge* pIBE = getNewIntBdryEdge();
  BdryEdgeBase * pBE = dynamic_cast<BdryEdgeBase *>(pIBE);
  pIBE->assignFaces(pFA, pFB);
  pFA->addCell(pIBE);
  pFB->addCell(pIBE);
  pFA->setFaceLoc(Face::eBdryTwoSide);
  pFB->setFaceLoc(Face::eBdryTwoSide);

  assert(pIBE->doFullCheck());
  assert(pIBE->getFace(0) == pFA);
  assert(pIBE->getFace(1) == pFB);
  assert(pFA->getLeftCell() == pIBE || pFA->getRightCell() == pIBE);
  assert(pFB->getLeftCell() == pIBE || pFB->getRightCell() == pIBE);

#ifndef OMIT_VERTEX_HINTS
  // Make sure that bdry verts have bdry hints
  pFA->getVert(0)->setHintFace(pFA);
  pFA->getVert(1)->setHintFace(pFA);
#endif
  // This is a really stinky way to copy bdry data, but that's what
  // works currently in 2D.
  if (pIBEOld->isValid()) {

    if (!(pBEOld->doesCurveExist())) {
      dynamic_cast<BdryEdgeBase *>(pIBE)->setPatch(
	  reinterpret_cast<BdryPatch*>(0x1)); // NULL would trigger an assertion...

    }
    else {
      GRCurve * curve = pBEOld->getCurve();

      // dynamic_cast<BdryEdgeBase *>(pBF)->
      //  setPatch (reinterpret_cast<BdryPatch*>(0x1)); // NULL would trigger an assertion...
      pBE->setCurve(curve);
      pBE->setGeometry(pBEOld->getGeometry());
      CubitVector split_coord0 = CubitVector(pIBE->getVert(0)->x(),
					     pIBE->getVert(0)->y(), 0);
      CubitVector split_coord1 = CubitVector(pIBE->getVert(1)->x(),
					     pIBE->getVert(1)->y(), 0);

      if (pBE->isForward()) {
	pBE->setVert0Param(curve->u_from_position(split_coord0));
	//  pBE->set_vert1_param(pBE->get_curve()->u_from_position(split_coord1));
	pBE->setVert1Param(curve->u_from_position(split_coord1));
      }
      else {
	pBE->setVert1Param(curve->u_from_position(split_coord0));
	pBE->setVert0Param(curve->u_from_position(split_coord1));
      }

      pBE->setBdryCondition(pBEOld->getBdryCondition());
    }

  }
  return pIBE;
}

//GRUMMP::PartBdryEdge* Mesh2D::createPartBdryEdge(Face* pF, GR_index_t otherPart,
//                                                bool owned)
//{
//    // TODO: Not thread-able.
//    assert(pF->isValid());
//    Vert *pV0 = pF->getVert(0);
//    Vert *pV1 = pF->getVert(1);
//
//    // Create a new edge with part bdry characteristics.
//    GRUMMP::PartBdryEdge *pPBE = getNewPartBdryEdge();
//    pPBE->setVerts(pV0, pV1);
//    Cell *pCL = pF->getLeftCell();
//    Cell *pCR = pF->getRightCell();
//    pPBE->setLeftCell(pCL);
//    pPBE->setRightCell(pCR);
//
//    // TODO: Assumes only one part per process.
//    if (owned) {
//        pPBE->setOwnerID(getPartID());
//        pPBE->setCopyID(GRUMMP::EntityGID(GRUMMP::PartGID(otherPart, 0)));
//    }
//    else {
//        pPBE->setCopyID(getPartID());
//        pPBE->setOwnerID(GRUMMP::PartGID(otherPart, 0));
//    }
//
//    // Now switch the old edge for the new one in cell connectivity.
//    if (pCL != pCInvalidCell) {
//        assert(pCR == pCInvalidCell);
//        pCL->replaceFace(pF, pPBE);
//        pF->removeCell(pCL);
//    }
//    else {
//        assert(pCR != pCInvalidCell);
//        pCR->replaceFace(pF, pPBE);
//        pF->removeCell(pCR);
//    }
//
//    // Destroy the old face
//    deleteFace(pF);
//
//    assert(pPBE->isValid());
//
//    return pPBE;
//}

// Cell* Mesh2D::createCell(Face * const apF[], const int iNFaces, const int iReg)
// // Faces are stored in canonical order.  Also, face->cell connectivity
// // is set up as part of this function.
// {
//   assert(iNFaces == 3 || iNFaces == 4);
//   if (iNFaces == 3) {
//     // Faces must already exist
//     return createTriCell(apF[0], apF[1], apF[2], iReg);
//   }
//   else {
//     return createQuadCell(apF[0], apF[1], apF[2], apF[3], iReg);
//   }
// }

TriCell*
Mesh2D::createTriCell(Face * const pF0, Face * const pF1, Face * const pF2,
		      const int iReg, int ID)
{
//	if (omp_in_parallel() != 0) assert(omp_get_num_threads() == NUM_PROCS);
//	assert(ID == omp_get_thread_num());
  // Faces must already exist
  assert(pF0->isValid());
  assert(pF1->isValid());
  assert(pF2->isValid());

  // Set up the face->cell connectivity properly.  The canonical tri
  // looks like this:
  //
  //                2
  //               / \         .
  //              /   \        .
  //             2     1       .
  //            /       \      .
  //           0 - -0- - 1

  // Each pair of faces  must have a common vert.
  Vert *apV[3];
  TriCell *pC;
  apV[0] = findCommonVert(pF0, pF2);
  apV[1] = findCommonVert(pF0, pF1);
  apV[2] = findCommonVert(pF1, pF2);
  assert(apV[0]->isValid());
  assert(apV[1]->isValid());
  assert(apV[2]->isValid());
  //#pragma omp critical(GNC1)
  //  {
  pC = dynamic_cast<TriCell*>(getNewCell(3, ID));
  //  }
  pC->setRegion(iReg);
  int iOrient = checkOrient2D(apV[0], apV[1], apV[2]);
  if (iOrient == 1) {
    // Right-handed (as sketched above)
    pC->assignFaces(pF0, pF1, pF2);
  }
  else {
    // Left-handed (faces 1 and 2 reversed above)
    pC->assignFaces(pF0, pF2, pF1);
  }
  if (XOR(iOrient == 1, (pF0->getVert(0) == apV[0]))) {
    pF0->setRightCell(pC);
  }
  else {
    pF0->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF1->getVert(0) == apV[1]))) {
    pF1->setRightCell(pC);
  }
  else {
    pF1->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF2->getVert(0) == apV[2]))) {
    pF2->setRightCell(pC);
  }
  else {
    pF2->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF0->getVert(0) == apV[0]))) {
    pF0->setRightCell(pC);
  }
  else {
    pF0->setLeftCell(pC);
  }

  createCellEvent(pC);
//	}
  return pC;
}

TriCell*
Mesh2D::createTriCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
		      const int iReg, const bool autoRotate, int ID)
{
//	if (omp_in_parallel() != 0) assert(omp_get_num_threads() == NUM_PROCS);
//	assert(ID == omp_get_thread_num());
  bool qExist;
  Face *pF01, *pF12, *pF20;
  if (!autoRotate) {
    pF01 = createFace(qExist, pV0, pV1, false, ID);
    pF12 = createFace(qExist, pV1, pV2, false, ID);
    pF20 = createFace(qExist, pV2, pV0, false, ID);
    // Some faces will exist; others won't.
  }
  else {
    double vector0[] =
      { pV0->x() - pV1->x(), pV0->y() - pV1->y() };
    double vector1[] =
      { pV1->x() - pV2->x(), pV1->y() - pV2->y() };
    double vector2[] =
      { pV2->x() - pV0->x(), pV2->y() - pV0->y() };
    double e0 = dMAG2D(vector0);
    double e1 = dMAG2D(vector1);
    double e2 = dMAG2D(vector2);

    if (e0 <= e1 && e0 <= e2) {
      pF01 = createFace(qExist, pV1, pV2, false, ID);
      pF12 = createFace(qExist, pV2, pV0, false, ID);
      pF20 = createFace(qExist, pV0, pV1, false, ID);
    }
    else if (e1 <= e0 && e1 <= e2) {
      pF01 = createFace(qExist, pV2, pV0, false, ID);
      pF12 = createFace(qExist, pV0, pV1, false, ID);
      pF20 = createFace(qExist, pV1, pV2, false, ID);
    }
    else if (e2 <= e0 && e2 <= e1) {
      pF01 = createFace(qExist, pV0, pV1, false, ID);
      pF12 = createFace(qExist, pV1, pV2, false, ID);
      pF20 = createFace(qExist, pV2, pV0, false, ID);
    }
  }
  TriCell *pC = createTriCell(pF01, pF12, pF20, iReg, ID);
  assert(pC->hasVert(pV0));
  assert(pC->hasVert(pV1));
  assert(pC->hasVert(pV2));

  return (pC);
}

TriCell*
Mesh2D::createTriCell(int minlayerindex, Vert * const pV0, Vert * const pV1,
		      Vert * const pV2, const int iReg, const bool autoRotate,
		      int ID)
{
//	if (omp_in_parallel() != 0) assert(omp_get_num_threads() == NUM_PROCS);
//	assert(ID == omp_get_thread_num());
  bool qExist;
  Face *pF01, *pF12, *pF20;
  if (!autoRotate) {
    pF01 = createFace(minlayerindex, qExist, pV0, pV1, false, ID);
    pF12 = createFace(minlayerindex, qExist, pV1, pV2, false, ID);
    pF20 = createFace(minlayerindex, qExist, pV2, pV0, false, ID);
    // Some faces will exist; others won't.
  }
  else {
    double vector0[] =
      { pV0->x() - pV1->x(), pV0->y() - pV1->y() };
    double vector1[] =
      { pV1->x() - pV2->x(), pV1->y() - pV2->y() };
    double vector2[] =
      { pV2->x() - pV0->x(), pV2->y() - pV0->y() };
    double e0 = dMAG2D(vector0);
    double e1 = dMAG2D(vector1);
    double e2 = dMAG2D(vector2);

    if (e0 <= e1 && e0 <= e2) {
      pF01 = createFace(minlayerindex, qExist, pV1, pV2, false, ID);
      pF12 = createFace(minlayerindex, qExist, pV2, pV0, false, ID);
      pF20 = createFace(minlayerindex, qExist, pV0, pV1, false, ID);
    }
    else if (e1 <= e0 && e1 <= e2) {
      pF01 = createFace(minlayerindex, qExist, pV2, pV0, false, ID);
      pF12 = createFace(minlayerindex, qExist, pV0, pV1, false, ID);
      pF20 = createFace(minlayerindex, qExist, pV1, pV2, false, ID);
    }
    else if (e2 <= e0 && e2 <= e1) {
      pF01 = createFace(minlayerindex, qExist, pV0, pV1, false, ID);
      pF12 = createFace(minlayerindex, qExist, pV1, pV2, false, ID);
      pF20 = createFace(minlayerindex, qExist, pV2, pV0, false, ID);
    }
  }
  TriCell *pC = createTriCell(pF01, pF12, pF20, iReg, ID);

  assert(pC->hasVert(pV0));
  assert(pC->hasVert(pV1));
  assert(pC->hasVert(pV2));

  return (pC);
}

QuadCell*
Mesh2D::createQuadCell(Vert * const pVc0, Vert * const pVc1, Vert * const pVc2,
		       Vert * const pVc3, const int iReg)
{

  Vert * pV0 = pVc0;
  Vert * pV1 = pVc1;
  Vert * pV2 = pVc2;
  Vert * pV3 = pVc3;

  Vert * pVtemp1;
  Vert * pVtemp2;
  // Change order of vertices if required to handle an obtuse angle. Not currently sure why we need to do this.
  int sumOrients = 0;
  sumOrients += checkOrient2D(pV0, pV1, pV2);
  sumOrients += checkOrient2D(pV1, pV2, pV3);
  sumOrients += checkOrient2D(pV2, pV3, pV0);
  sumOrients += checkOrient2D(pV3, pV0, pV1);
  if (abs(sumOrients) == 0) {
    logMessage(MSG_DETAIL, "Flipped quad...\n");
    logMessage(MSG_DETAIL, "%f: %f\n", pV0->x(), pV0->y());
    logMessage(MSG_DETAIL, "%f: %f\n", pV1->x(), pV1->y());
    logMessage(MSG_DETAIL, "%f: %f\n", pV2->x(), pV2->y());
    logMessage(MSG_DETAIL, "%f: %f\n", pV3->x(), pV3->y());
  }
  int signOrients = sumOrients / abs(sumOrients);
  if (abs(sumOrients) != 4) {
    if (checkOrient2D(pV0, pV1, pV2) != signOrients) {
      pVtemp1 = pV0;
      pV0 = pV3;
      pV3 = pV2;
      pV2 = pV1;
      pV1 = pVtemp1;
    }
    else if (checkOrient2D(pV3, pV0, pV1) != signOrients) {
      pVtemp1 = pV0;
      pVtemp2 = pV1;
      pV0 = pV2;
      pV1 = pV3;
      pV2 = pVtemp1;
      pV3 = pVtemp2;
    }
    else if (checkOrient2D(pV2, pV3, pV0) != signOrients) {
      pVtemp1 = pV0;
      pV0 = pV1;
      pV1 = pV2;
      pV2 = pV3;
      pV3 = pVtemp1;
    }
  }
  /*
   if(abs(sumOrients) != 4){
   if(checkOrient2D(pV0, pV1, pV2) != signOrients){
   pVtemp1 = pV0;
   pV0 = pV1;
   pV1 = pV2;
   pV2 = pV3;
   pV3 = pVtemp1;
   }
   else if(checkOrient2D(pV1, pV2, pV3) != signOrients){
   pVtemp1 = pV0;
   pVtemp2 = pV1;
   pV0 = pV2;
   pV1 = pV3;
   pV2 = pVtemp1;
   pV3 = pVtemp2;
   }
   else if(checkOrient2D(pV2, pV3, pV0) != signOrients){
   pVtemp1 = pV0;
   pV0 = pV3;
   pV3 = pV2;
   pV2 = pV1;
   pV1 = pVtemp1;
   }
   }
   */

  bool qExist;
  Face *pF01 = createFace(qExist, pV0, pV1);
  Face *pF12 = createFace(qExist, pV1, pV2);
  Face *pF23 = createFace(qExist, pV2, pV3);
  Face *pF30 = createFace(qExist, pV3, pV0);
  // Some faces will exist; others won't.

  QuadCell *pC = createQuadCell(pF01, pF12, pF23, pF30, iReg);
  assert(pC->hasVert(pV0));
  assert(pC->hasVert(pV1));
  assert(pC->hasVert(pV2));
  assert(pC->hasVert(pV3));
  return (pC);
}

QuadCell*
Mesh2D::createQuadCell(Face * const pF0_in, Face * const pF1_in,
		       Face * const pF2_in, Face * const pF3_in, const int iReg)
{
  // Faces must already exist
  assert(pF0_in->isValid());
  assert(pF1_in->isValid());
  assert(pF2_in->isValid());
  assert(pF3_in->isValid());

  // Set up the face->cell connectivity properly.  The canonical quad
  // looks like this (see ITAPS docs):
  //
  //      3 - -2- - 2
  //       \         \         .
  //        \         \        .
  //         3         1       .
  //          \         \      .
  //           0 - -0- - 1

  QuadCell *pC = dynamic_cast<QuadCell*>(getNewCell(4));
  pC->setRegion(iReg);

  // The given order is supposed to be canonical; this call makes sure
  // internally.
  pC->assignFaces(pF0_in, pF1_in, pF2_in, pF3_in);

  Face *pF0 = pC->getFace(0);
  Face *pF1 = pC->getFace(1);
  Face *pF2 = pC->getFace(2);
  Face *pF3 = pC->getFace(3);

  Vert* apV[4];
  apV[0] = findCommonVert(pF3, pF0);
  apV[1] = findCommonVert(pF0, pF1);
  apV[2] = findCommonVert(pF1, pF2);
  apV[3] = findCommonVert(pF2, pF3);

  int iOrient = checkOrient2D(apV[0], apV[1], apV[2]);
  if (checkOrient2D(apV[0], apV[2], apV[3]) != iOrient) {
    printf("%d %d\n", iOrient, checkOrient2D(apV[0], apV[2], apV[3]));
  }
  assert(checkOrient2D(apV[0], apV[2], apV[3]) == iOrient);

  if (XOR(iOrient == 1, (pF0->getVert(0) == apV[0]))) {
    pF0->setRightCell(pC);
  }
  else {
    pF0->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF1->getVert(0) == apV[1]))) {
    pF1->setRightCell(pC);
  }
  else {
    pF1->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF2->getVert(0) == apV[2]))) {
    pF2->setRightCell(pC);
  }
  else {
    pF2->setLeftCell(pC);
  }

  if (XOR(iOrient == 1, (pF3->getVert(0) == apV[3]))) {
    pF3->setRightCell(pC);
  }
  else {
    pF3->setLeftCell(pC);
  }

  if (iOrient == 1) {
    // Right-handed (as done above by canonicalization)
    // pC->vAssign(pF0, pF1, pF2, pF3);
  }
  else {
    // Left-handed (faces 1 and 3 reversed from above)
    pC->assignFaces(pF0, pF3, pF2, pF1);
  }

  assert(pC->doFullCheck());
  assert(pC->hasFace(pF0));
  assert(pC->hasFace(pF1));
  assert(pC->hasFace(pF2));
  assert(pC->hasFace(pF3));
  assert(pF0->doFullCheck());
  assert(pF1->doFullCheck());
  assert(pF2->doFullCheck());
  assert(pF3->doFullCheck());

  createCellEvent(pC);
  return pC;
}

// Cell* Mesh2D::createCell(Vert* const apV[], const int iNVerts, const int iReg)
// {
//   if (iNVerts == 3) {
//     return (createTriCell(apV[0], apV[1], apV[2], iReg));
//   }
//   else {
//     assert(iNVerts == 4);
//     return (createQuadCell(apV[0], apV[1], apV[2], apV[3], iReg));
//   }
// }

//bool Mesh2D::deleteCell(Cell *pC)
//{
//  if (pC->isDeleted()) return true;
//  // Free up connectivity for faces.
//   Face* pF;
////#pragma omp threadprivate (pF)
//  for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
//	  pF = pC->getFace(i);
//    if (pF->isValid()) {
//      pF->removeCell(pC);
//      if (!pF->getLeftCell()->isValid() &&
//	  !pF->getRightCell()->isValid()) deleteFace(pF);
//    }
//  }
//  // Under the old data structures, just mark it.
//  pC->markAsDeleted();
//  // With the new (2010) containers, also register the deletion with the
//  // container.
//  switch (pC->getType()) {
//  case CellSkel::eTriCell:
//    m_ECTri.deleteEntry(static_cast<TriCell*>(pC));
//    break;
//  case CellSkel::eQuadCell:
//    m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC));
//    break;
//  default:
//    assert(0);
//    break;
//  }
//  deleteCellEvent(pC);
//  return true;
//}

bool
Mesh2D::deleteCellRe(Cell * pC, int ID)
{
//	assert(ID == omp_get_thread_num());
//	assert(omp_get_num_threads() == NUM_PROCS);
  if (!pC->isValid() || pC->isDeleted())
    return true;
  // Free up connectivity for faces.
  Face* pF;
  //#pragma omp critical(deleteCellRe1234567)
  //{
  for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
    pF = pC->getFace(i);
    if (pF->isValid()) {
      pF->removeCell(pC);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid())
	deleteFaceRe(pF, ID);
    }
    //}
  }
  // Under the old data structures, just mark it.
  pC->markAsDeleted();
  // With the new (2010) containers, also register the deletion with the
  // container.
  switch (pC->getType())
    {
    case CellSkel::eTriCell:
      //#pragma omp critical(deleteCellRedeleteEntry1)
      m_ECTri.deleteEntry(static_cast<TriCell*>(pC), ID);
      break;
    case CellSkel::eQuadCell:
      m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC), ID);
      break;
    default:
      assert(0);
      break;
    }
//#pragma omp critical(Event)
//	{
  deleteCellEvent(pC);
//	}
  return true;
}

bool
Mesh2D::deleteCell(Cell * const pC)
{
  if (!pC->isValid())
    return true;
  if (omp_in_parallel() == 0) {
    if (pC->isDeleted())
      return true;
    // Free up connectivity for faces.
    Face* pF;
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      pF = pC->getFace(i);
      if (pF->isValid()) {
	pF->removeCell(pC);
	if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid())
	  deleteFace(pF);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTriCell:

	m_ECTri.deleteEntry(static_cast<TriCell*>(pC));
	break;
      case CellSkel::eQuadCell:
	m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC));
	break;
      default:
	assert(0);
	break;
      }
    deleteCellEvent(pC);
    return true;
  }
  else {
    int ID = omp_get_thread_num();
    if (pC->isDeleted())
      return true;
    // Free up connectivity for faces.
    Face* pF;
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      pF = pC->getFace(i);
      if (pF->isValid()) {
	pF->removeCell(pC);
	if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid())
	  deleteFaceRe(pF, ID);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTriCell:
	m_ECTri.deleteEntry(static_cast<TriCell*>(pC), ID);
	break;
      case CellSkel::eQuadCell:
	m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC), ID);
	break;
      default:
	assert(0);
	break;
      }
//		#pragma omp critical(Event)
//		  {
    deleteCellEvent(pC);
//		  }
    return true;
  }
  assert(0);
  return false;
}

bool
Mesh2D::deleteCellWithoutDeletingFaces(Cell * const pC)
{
  if (!pC->isValid())
    return true;
  if (omp_in_parallel() == 0) {
    if (pC->isDeleted())
      return true;
    // Free up connectivity for faces.
    Face* pF;
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      pF = pC->getFace(i);
      if (pF->isValid()) {
	pF->removeCell(pC);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTriCell:

	m_ECTri.deleteEntry(static_cast<TriCell*>(pC));
	break;
      case CellSkel::eQuadCell:
	m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC));
	break;
      default:
	assert(0);
	break;
      }
    deleteCellEvent(pC);
    return true;
  }
  else {
    int ID = omp_get_thread_num();
    if (pC->isDeleted())
      return true;
    // Free up connectivity for faces.
    Face* pF;
    for (int i = pC->getNumFaces() - 1; i >= 0; i--) {
      pF = pC->getFace(i);
      if (pF->isValid()) {
	pF->removeCell(pC);
      }
    }
    // Under the old data structures, just mark it.
    pC->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pC->getType())
      {
      case CellSkel::eTriCell:
	m_ECTri.deleteEntry(static_cast<TriCell*>(pC), ID);
	break;
      case CellSkel::eQuadCell:
	m_ECQuad.deleteEntry(static_cast<QuadCell*>(pC), ID);
	break;
      default:
	assert(0);
	break;
      }
//		#pragma omp critical(Event)
//		  {
    deleteCellEvent(pC);
//		  }
    return true;
  }
  assert(0);
  return false;
}

// This function is not technically necessary, under the old data
// structures.
bool
Mesh2D::deleteBFace(BFace * const pBF)
{
  if (!pBF->isValid())
    return true;
  if (omp_in_parallel() == 0) {
    if (pBF->isDeleted())
      return true;
    Face *pF = pBF->getFace();
    if (pF->isValid()) {
      pF->removeCell(pBF);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid())
	deleteFace(pF);
    }

    if (pBF->getNumFaces() == 2 && pBF->getFace(1)->isValid()) {
      pF = pBF->getFace(1);
      pF->removeCell(pBF);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
	deleteFace(pF);
      }
    }

    // Under the old data structures, just mark it.
    pBF->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pBF->getType())
      {
      case CellSkel::eBdryEdge:
	m_ECEdgeBF.deleteEntry(static_cast<BdryEdge*>(pBF));
	break;
      case CellSkel::eIntBdryEdge:
	m_ECIntEdgeBF.deleteEntry(static_cast<IntBdryEdge*>(pBF));
	break;
      default:
	assert(0);
	break;
      }
    deleteBFaceEvent(pBF);
    return true;
  }
  else {
    int ID = omp_get_thread_num();
    if (pBF->isDeleted())
      return true;
    Face *pF = pBF->getFace();
    if (pF->isValid()) {
      pF->removeCell(pBF);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid())
	deleteFaceRe(pF, ID);	// ********* Parallel Changes ***********
    }

    if (pBF->getNumFaces() == 2 && pBF->getFace(1)->isValid()) {
      pF = pBF->getFace(1);
      pF->removeCell(pBF);
      if (!pF->getLeftCell()->isValid() && !pF->getRightCell()->isValid()) {
	deleteFaceRe(pF, ID);	// ********* Parallel Changes ***********
      }
    }

    // Under the old data structures, just mark it.
    pBF->markAsDeleted();
    // With the new (2010) containers, also register the deletion with the
    // container.
    switch (pBF->getType())
      {
      case CellSkel::eBdryEdge:
	m_ECEdgeBF.deleteEntry(static_cast<BdryEdge*>(pBF), ID); // ********* Parallel Changes ***********
	break;
      case CellSkel::eIntBdryEdge:
	m_ECIntEdgeBF.deleteEntry(static_cast<IntBdryEdge*>(pBF), ID); // ********* Parallel Changes ***********
	break;
      default:
	assert(0);
	break;
      }
//#pragma omp critical(Event)														// ********* Parallel Changes ***********
//		{
    deleteBFaceEvent(pBF);
//		}
    return true;
  }
}

bool
Mesh2D::deleteFaceRe(Face * pF, int ID)
{
//	assert(omp_get_num_threads() == NUM_PROCS);
  if (!pF->isValid() || pF->isDeleted())
    return true;
  if ((pF->getRightCell()->isValid() || pF->getLeftCell()->isValid()))
    return false;
  //#pragma omp critical(deleteFaceReremoveFace)
  //{
  pF->getVert(1)->removeFace(pF);
  pF->getVert(0)->removeFace(pF);
  //}

#ifndef OMIT_VERTEX_HINTS
  // Make sure that verts whose hint is going away realize it...
  //#pragma omp critical(HintFace1234)
  //{
  if (pF->getVert(0)->getHintFace() == pF) {
    pF->getVert(0)->updateHintFace();
  }
  if (pF->getVert(1)->getHintFace() == pF) {
    pF->getVert(1)->updateHintFace();
  }
  //}
#endif

  // Under the old data structures, just mark it.
  pF->markAsDeleted();
  // With the new (2010) containers, also register the deletion with the
  // container.
//	assert(ID == omp_get_thread_num());
  //#pragma omp critical(deleteFaceRedeleteEntry1)
  m_ECEdgeF.deleteEntry(static_cast<EdgeFace*>(pF), ID);
//#pragma omp critical(Event)
//	{
  deleteFaceEvent(pF);
//	}
  return true;
}

// Like ITAPS, faces that are in use can't be deleted.
bool
Mesh2D::deleteFace(Face * const pF)
{
  assert(pF->getVert(0) != pF->getVert(1));

  if (!pF->isValid() || pF->isDeleted())
    return true;
  if ((pF->getRightCell()->isValid() || pF->getLeftCell()->isValid()) || ((pF->getFaceLoc() == Face::ePseudoSurface)) /*|| pF->isLocked()*/)
    return false;

  assert(pF->getFaceLoc() != Face::ePseudoSurface);
 // assert(!pF->isLocked());

  pF->getVert(1)->removeFace(pF);
  pF->getVert(0)->removeFace(pF);
#ifndef OMIT_VERTEX_HINTS
  // Make sure that verts whose hint is going away realize it...
  if (pF->getVert(0)->getHintFace() == pF) {
    pF->getVert(0)->updateHintFace();
  }
  if (pF->getVert(1)->getHintFace() == pF) {
    pF->getVert(1)->updateHintFace();
  }
#endif

  // Under the old data structures, just mark it.
  pF->markAsDeleted();
  // With the new (2010) containers, also register the deletion with the
  // container.
  m_ECEdgeF.deleteEntry(static_cast<EdgeFace*>(pF));
  deleteFaceEvent(pF);
  return true;
}

Vert*
Mesh2D::createVert(const double dX, const double dY, const double)
{
  double adData[] =
    { dX, dY };
  return createVert(adData);
}

Vert*
Mesh2D::createVert(const double adCoords[2])
{
  if (omp_in_parallel() == 0) {
    Vert *pV = getNewVert();
    pV->clearHintFace();
    pV->clearFaceConnectivity();
    pV->setDefaultFlags();
    pV->setCoords(2, adCoords);
    createVertEvent(pV);
    return pV;
  }
  else {
    int ID = omp_get_thread_num();
    Vert *pV = getNewVert(ID);
    pV->clearHintFace();
    pV->clearFaceConnectivity();
    pV->setDefaultFlags();
    pV->setCoords(2, adCoords);
//#pragma omp critical(Event)
//		{
    createVertEvent(pV);
//		}
    return pV;
  }
}

// For verts, it's non-trivial to keep the hint face up to date as faces
// are deleted, so there's no reliable test for when these are safe to
// delete.  So just mark it, and caveat emptor.
bool
Mesh2D::deleteVert(Vert * const pV)
{
  if (!pV->isValid())
    return true;
  if (omp_in_parallel() == 0) {
    if (pV->isDeleted())
      return true;
    if (pV->getNumFaces() == 0) {
      pV->clearFaceConnectivity();
      // Under the old data structures, just mark it.
      pV->markAsDeleted();
      // With the new (2010) containers, also register the deletion with the
      // container.
      m_ECVerts.deleteEntry(pV);
      deleteVertEvent(pV);
      return true;
    }
    else {
      return false;
    }
  }
  else {
    int ID = omp_get_thread_num();
    if (pV->isDeleted())
      return true;
    if (pV->getNumFaces() == 0) {
      pV->clearFaceConnectivity();
      // Under the old data structures, just mark it.
      pV->markAsDeleted();
      // With the new (2010) containers, also register the deletion with the
      // container.
      m_ECVerts.deleteEntry(pV, ID);
//			#pragma omp critical(Event)
//			{
      deleteVertEvent(pV);
//			}
      return true;
    }
    else {
      return false;
    }
  }
}

void
Mesh2D::resetVertexConnectivity()
{
  for (int i = getNumVerts() - 1; i >= 0; i--) {
    Vert *pV = getVert(i);
    if (pV->isDeleted())
      continue;
    pV->clearFaceConnectivity();
  }

  for (int i = getNumFaces() - 1; i >= 0; i--) {
    Face *pF = getFace(i);
    if (pF->isDeleted())
      continue;
    for (int ii = pF->getNumVerts() - 1; ii >= 0; ii--) {
      pF->getVert(ii)->addFace(pF);
    }
  }
}

//Private call used by TriMeshBuilder and Watson boundary insertion
IntBdryEdge*
Mesh2D::createIntBdryEdge(Vert* const pVBeg, Vert* const pVEnd,
			  GRCurve* const curve, double dBegParam,
			  double dEndParam)
{
  int ID = 0;
//	if (omp_in_parallel() != 0)
//	{
//		assert(NUM_PROCS == omp_get_num_threads());
  ID = omp_get_thread_num();
//	}
  int iRegLeft;
  int iRegRight;
#pragma omp critical(CurveBdryEdge)
  {
    iRegLeft = curve->region_left();
    iRegRight = curve->region_right();
  }
  assert(iRegRight != 0 && iRegLeft != 0);

  // Find the face connecting the two verts, if it exists
  Face *pF = findCommonFace(pVBeg, pVEnd, true);
  if (pF->getVert(0) == pVEnd) {
    pF->interchangeCellsAndVerts();
  }

  assert(pF->getVert(0) == pVBeg);
  assert(pF->getVert(1) == pVEnd);
//	Cell * pCLeft  = pF->getLeftCell();
  Cell * pCRight = pF->getRightCell();

  //  assert(pCLeft->getRegion() == iRegLeft);
  //  assert(pCRight->getRegion() == iRegRight);
  //  pCLeft->setRegion(iRegLeft);
  //  pCRight->setRegion(iRegRight);

  IntBdryEdge* pIBENew = getNewIntBdryEdge(2, ID); // ******** changes in Parallel **********
  Face* pFNew = getNewFace(2, ID);  // ******** changes in Parallel **********

  pF->setFaceLoc(Face::eBdryTwoSide);
  pFNew->setFaceLoc(Face::eBdryTwoSide);

  pF->replaceCell(pCRight, pIBENew);

  pFNew->assign(pIBENew, pCRight, pVBeg, pVEnd);

  pCRight->replaceFace(pF, pFNew);

  pIBENew->assignFaces(pF, pFNew);
  pIBENew->setCurve(curve);
  pIBENew->setVert0Param(dBegParam);
  pIBENew->setVert1Param(dEndParam);
  pVBeg->setHintFace(pF);
  pVEnd->setHintFace(pF);
  if (pVBeg->getVertType() != Vert::eBdryApex)
    pVBeg->setType(Vert::eBdryTwoSide);
  if (pVEnd->getVertType() != Vert::eBdryApex)
    pVEnd->setType(Vert::eBdryTwoSide);
  assert(pIBENew->getVert(0) == pVBeg && pIBENew->getVert(1) == pVEnd);

//#pragma omp critical(Event)									// ******** changes in Parallel **********
//		  {
  createBFaceEvent(pIBENew);
//		  }
  return pIBENew;
}

BdryEdge*
Mesh2D::createBdryEdge(Vert* const pVBeg, Vert* const pVEnd,
		       GRCurve* const curve, double dBegParam, double dEndParam)
{
  int ID = 0;
//	if (omp_in_parallel() != 0)                              // ******** changes in Parallel **********
//	{
//		assert(NUM_PROCS == omp_get_num_threads());
  ID = omp_get_thread_num();
//	}
  int iRegLeft;
  int iRegRight;
#pragma omp critical(CurveBdryEdge)
  {
    iRegLeft = curve->region_left();
    iRegRight = curve->region_right();
  }

  assert(XOR(iRegRight != iInvalidRegion, iRegLeft != iInvalidRegion));
  //"Regular" boundary edge.

  Face* pF = findCommonFace(pVBeg, pVEnd, true);
  assert(pF->isValid());
  assert(pF->hasVert(pVBeg) && pF->hasVert(pVEnd));
  pF->setFaceLoc(Face::eBdryFace);
  BdryEdge* pBENew = dynamic_cast<BdryEdge*>(getNewBFace(2, ID)); // ******** changes in Parallel **********
  assert(pBENew);
  pBENew->setCurve(curve);
  pBENew->assignFaces(pF);

#ifndef OMIT_VERTEX_HINTS
  // Make sure that bdry verts have bdry hints
  pF->getVert(0)->setHintFace(pF);
  pF->getVert(1)->setHintFace(pF);
#endif

  // Ensure that the edge and curve have the same orientation
  if ((dBegParam < dEndParam && pVBeg == pF->getVert(1))
      || (dBegParam > dEndParam && pVBeg == pF->getVert(0))) {
    pF->interchangeCellsAndVerts();
  }

  if (iRegRight == iInvalidRegion) {
    // Domain is on the left.
    Cell *pC = pF->getRightCell();
    if (pC != pCInvalidCell) {
      pC->removeFace(pF);
      pC->setRegion(iOutsideRegion);
    }
    pF->setRightCell(pBENew);
    //  pF->setLeftCell(pC);

  }
  else {
    // Domain is on the right.
    Cell *pC = pF->getLeftCell();
    if (pC != pCInvalidCell) {
      pC->removeFace(pF);
      pC->setRegion(iOutsideRegion);
    }
    pF->setLeftCell(pBENew);
    //   pF->setRightCell(pC);
  }
  if (pBENew->getVert(0) == pF->getVert(0)) {
    pBENew->setVert0Param(dBegParam);
    pBENew->setVert1Param(dEndParam);
  }
  else {
    pBENew->setVert0Param(dEndParam);
    pBENew->setVert1Param(dBegParam);
  }

  assert(
      (pF->getVert(0) == pBENew->getVert(0)
	  && pBENew->getVert0Param() < pBENew->getVert1Param())
	  || (pF->getVert(0) == pBENew->getVert(1)
	      && pBENew->getVert0Param() > pBENew->getVert1Param()));

  assert(pF->doFullCheck() == 1);
  assert(pBENew->doFullCheck() == 1);
  assert(pBENew->getFace() == pF);
  assert(pF->getLeftCell() == pBENew || pF->getRightCell() == pBENew);

  if (pVBeg->getVertType() != Vert::eBdryApex)
    pVBeg->setType(Vert::eBdry);
  if (pVEnd->getVertType() != Vert::eBdryApex)
    pVEnd->setType(Vert::eBdry);

//#pragma omp critical(Event)										// ******** changes in Parallel **********
//	{
  createBFaceEvent(pBENew);
//	}
  return pBENew;

}

bool
Mesh2D::Orientation(TriCell* Tri1, TriCell* Tri2)
{
  Face* commenface = findCommonFace(Tri1, Tri2);
  const double *r = Tri1->getOppositeVert(commenface)->getCoords();
  bool ori;
  int minlayerorder0;
  int minlayerorder1;
  if (Tri1->getLargestLayernumber() < 1e5
      and Tri2->getLargestLayernumber() < 1e5) {
    minlayerorder0 = INT_MAX;
    for (int i = 0; i < commenface->getVert(0)->getNumFaces(); i++) {
      minlayerorder0 = fmin(minlayerorder0,
			    commenface->getVert(0)->getFace(i)->layerindex);
    }
    minlayerorder1 = INT_MAX;
    for (int i = 0; i < commenface->getVert(1)->getNumFaces(); i++) {
      minlayerorder1 = fmin(minlayerorder1,
			    commenface->getVert(1)->getFace(i)->layerindex);
    }
    assert(minlayerorder0 != minlayerorder1);
  }
  else {
    minlayerorder0 = 0;
    minlayerorder1 = 0;
    for (int i = 0; i < commenface->getVert(0)->getNumFaces(); i++) {
      minlayerorder0 = fmax(minlayerorder0,
			    commenface->getVert(0)->getFace(i)->layerindex);
    }
    for (int i = 0; i < commenface->getVert(1)->getNumFaces(); i++) {
      minlayerorder1 = fmax(minlayerorder1,
			    commenface->getVert(1)->getFace(i)->layerindex);
    }
    assert(minlayerorder0 != minlayerorder1);
  }
  const double *p, *q;
  if (minlayerorder0 < minlayerorder1) {
    p = commenface->getVert(0)->getCoords();
    q = commenface->getVert(1)->getCoords();

  }
  else {
    p = commenface->getVert(1)->getCoords();
    q = commenface->getVert(0)->getCoords();
  }

  double delta = (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]);
  ori = (delta > 0.);
  if ((Tri1->getSmallestLayernumber()) % 2 == 0)
    return ori;
  else
    return !ori;
}

void
Mesh2D::TriMergetoQuad()
{
  ///
  int mergenumber = 0;
  std::cout << "Merge Game on" << std::endl;
  ///Queue All cells in tri mesh
  ///Queue For Merge Triangle to Quad;
  InsertionQueue MergeQ;
  //MergeQ(this, INT_MAX, true);
  for (GR_index_t i = 0; i < getNumTriCells(); i++) {
    TriCell* triangle = (dynamic_cast<TriCell*>(getCell(i)));
    double layernumber = triangle->getLargestLayernumber();
    layernumber += 1e-6 * triangle->getSmallestLayernumber();
    MergeQ.qAddEntrywithPriority(triangle, layernumber);
  }
//	MergeQ.print_queue_to_Screen();
  //Queue End
  while (!MergeQ.qIsQueueEmpty()) {
    bool Tri1 = false;
    bool Tri2 = false;
    ///Grab top cells to start merge
    InsertionQueueEntry IQE = MergeQ.IQETopValidEntry();
    assert(IQE.eType() == InsertionQueueEntry::eTriCell);
    TriCell* startingTri = dynamic_cast<TriCell*>(IQE.pCCell());
    if (IQE.pCCell()->isDeleted()) {
      std::cout << "this judge works" << std::endl;
      MergeQ.vPopTopEntry();
      continue;
    }
    else {
//			assert(startingTri->getSmallestLayernumber() == 0);
      Tri1 = true;
    }
    assert(startingTri->isValid());
    //Naive merge
    TriCell* secondTri = NULL;
    int layernumber = startingTri->getLargestLayernumber();
    int secondlaynumber = 0;
    ///start with first layer
//		if (layernumber > 1){
//			break;
//		}
    for (int i = 0; i < 3; i++) {
      if (startingTri->getFace(i)->isBdryFace()) {
	continue;
      }
      secondTri =
	  dynamic_cast<TriCell*>(startingTri->getFace(i)->getOppositeCell(
	      startingTri));
      if (secondTri->isValid()) {
	secondlaynumber =
	    dynamic_cast<TriCell*>(startingTri->getFace(i)->getOppositeCell(
		startingTri))->getLargestLayernumber();
      }
      else {
	continue;
      }
      if (secondlaynumber == layernumber or secondlaynumber > 1e4) {
	///Check orientation

	if (Orientation(startingTri, secondTri) == true) {
	  if (!secondTri->isDeleted()) {
	    Tri2 = true;
	    break;
	  }
	}
      }
    }
    if (Tri1 and Tri2) {
      createQuadFromTris(startingTri, secondTri);
      mergenumber++;
    }
    else {
    }
    MergeQ.vPopTopEntry();
    /////write temp mesh
//		char strName[FILE_NAME_LEN];
//		sprintf(strName, "%s_%d",m_tempMeshFileName,m_tempMeshIndex);
//		writeVTKLegacy(*this, strName);
//		m_tempMeshIndex++;
    /////

  }
  //
  std::cout << "Merge Game ends(hopefully it worked) :)" << std::endl;
  std::cout << "Merged " << mergenumber << " Quads" << std::endl;
}

void
Mesh2D::TriMergetoQuad2()
{
  ///
  double refvalue = 0.7;
  int mergenumber = 0;
  std::cout << "Merge Game2 on" << std::endl;
  ///Queue All cells in tri mesh
  ///Queue For Merge Triangle to Quad;
  InsertionQueue MergeQ;
  //MergeQ(this, INT_MAX, true);
  /////Queue all the edges with only one tri adjacent to it.
  for (GR_index_t i = 0; i < getNumFaces(); i++) {
    getFace(i)->markAsMergeable();
    Face* baseedge = getFace(i);
    if (baseedge->isBdryFace()) {
      baseedge->markAsNotMergeable();
      double edgelength = baseedge->calcSize();
      MergeQ.qAddEntrywithPriority(baseedge, edgelength - 50);
    }
  }
//	MergeQ.print_queue_to_Screen();
  //Queue End
  while (!MergeQ.qIsQueueEmpty()) {
    Face* QuadFaces[4] =
      { NULL };
    bool Tri1 = false;
    bool Tri2 = false;
    ///Grab top cells to start merge
    InsertionQueueEntry IQE = MergeQ.IQETopValidEntry();
    assert(IQE.eType() == InsertionQueueEntry::eFace);
    Face* baseedge = dynamic_cast<Face*>(IQE.pFFace());
    if (!baseedge->isValid() or baseedge->isDeleted()) {
      MergeQ.vPopTopEntry();
      continue;
    }
    TriCell* startingTri = dynamic_cast<TriCell*>(
	baseedge->getCell(0)->getType() == Cell::eTriCell ?
	    baseedge->getCell(0) : baseedge->getCell(1));
    if (!startingTri->isValid() or startingTri->getType() != Cell::eTriCell) {
//			std::cout << "this judge works" << std::endl;
      MergeQ.vPopTopEntry();
      continue;
    }
    else {

      Tri1 = true;
    }
    Face* mergeface = NULL;
    //Naive merge
    int nummergeable = 0;
    TriCell* secondTri = NULL;
    for (int i = 0; i < 3; i++) {
      if (startingTri->getFace(i) == baseedge) {
	continue;
      }
      else if (startingTri->getFace(i)->isBdryFace()) {
	continue;
      }
      else if (startingTri->getFace(i)->isMergeable()) {
	mergeface = startingTri->getFace(i);
	nummergeable++;
      }
      else {
      }
    }
//		printf("nummergeable = %d\n", nummergeable);
    /////if nothing to be merged, pop the top entry in the queue
    if (nummergeable == 0) {
      MergeQ.vPopTopEntry();
      continue;
    }
    assert(mergeface);
    secondTri =
	mergeface->getCell(0) == startingTri ?
	    dynamic_cast<TriCell*>(mergeface->getCell(1)) :
	    dynamic_cast<TriCell*>(mergeface->getCell(0));
    ///reject the merge if it generates bad shape quad.
    Vert* v0 = startingTri->getOppositeVert(mergeface);
    Vert* v2 = secondTri->getOppositeVert(mergeface);
    Vert* v1 = mergeface->getVert(0);
    Vert* v3 = mergeface->getVert(1);
    double mincos = calcsquareness(v0, v1, v2, v3);
    if (nummergeable == 1) {
      if (mincos < refvalue) {
	if (checkOrient2D(v0, v1, v2) == checkOrient2D(v0, v2, v3)) {
	  Tri2 = true;
	}
	else {
	  MergeQ.vPopTopEntry();
	  continue;
	}
      }
      else {
	MergeQ.vPopTopEntry();
	continue;

      }
    }
    if (nummergeable == 2) {
      Face* altermergeface = NULL;
      for (int i = 0; i < 3; i++) {
	if ((startingTri->getFace(i)->isBdryFace())
	    or (startingTri->getFace(i) == mergeface)
	    or (startingTri->getFace(i) == baseedge)) {
	  continue;
	}
	else {
	  altermergeface = startingTri->getFace(i);
	}
      }
      assert(altermergeface);
      TriCell *alterTri =
	  dynamic_cast<TriCell*>(altermergeface->getOppositeCell(startingTri));
      Vert* v01 = startingTri->getOppositeVert(altermergeface);
      Vert* v21 = alterTri->getOppositeVert(altermergeface);
      Vert* v11 = altermergeface->getVert(0);
      Vert* v31 = altermergeface->getVert(1);
      double mincos2 = calcsquareness(v01, v11, v21, v31);
      /////New way to check chevron pattern///
      /*
       * VC                                VD
       * ----------------------------------
       *  \              / \              /
       *    \  second   /    \  alter    /
       *      \        /       \        /
       *        \     /  starting\     /
       *          \  /             \  /
       *          VA\--------------- \ VB
       *           /                /
       *          /                /
       *         /   quad         /
       *        /                /
       *       /                /
       *      /----------------/
       *      VE               VF
       * Check the product of CE*AB and DF*AB
       */
      Vert* VA = NULL;
      Vert* VB = NULL;
      if (mergeface->hasVert(baseedge->getVert(0))) {
	VA = baseedge->getVert(0);
	VB = baseedge->getVert(1);
      }
      else {
	VA = baseedge->getVert(1);
	VB = baseedge->getVert(0);
      }
      assert(mergeface->hasVert(VA));
      assert(altermergeface->hasVert(VB));
//			Vert* VC = secondTri->getOppositeVert(mergeface);
//			Vert* VD = alterTri->getOppositeVert(mergeface);
//			if (!baseedge->isBdryFace()){
//				QuadCell *quad = dynamic_cast<QuadCell*>(baseedge->getOppositeCell(startingTri));
//				assert(quad->getType() == Cell::eQuadCell);
//				Vert* VE = NULL;
//				Vert* VF = NULL;
//				for (int i = 0; i < 4; i++){
//					if (quad->getFace(i) == baseedge){
//						continue;
//					}
//					else if(quad->getFace(i)->hasVert(VA)){
//						VE = quad->getFace(i)->getOppositeVert(VA);
//					}
//					else if(quad->getFace(i)->hasVert(VB)){
//						VF = quad->getFace(i)->getOppositeVert(VB);
//					}
//				}
////				double normAB[2] = {((VA->getCoords())[XDIR] - (VB->getCoords())[XDIR])/baseedge->calcSize(),((VA->getCoords())[YDIR] - (VB->getCoords())[YDIR])/baseedge->calcSize()};
////				double distCE = (VC->getCoords()[XDIR]-VE->getCoords()[XDIR])*(VC->getCoords()[XDIR]-VE->getCoords()[XDIR]) + (VC->getCoords()[YDIR]-VE->getCoords()[YDIR])*(VC->getCoords()[YDIR]-VE->getCoords()[YDIR]);
////				distCE = sqrt(distCE);
////				double normCE[2] = {((VC->getCoords())[XDIR] - (VE->getCoords())[XDIR])/distCE,((VC->getCoords())[YDIR] - (VE->getCoords())[YDIR])/distCE};
////				double distDF = (VD->getCoords()[XDIR]-VF->getCoords()[XDIR])*(VD->getCoords()[XDIR]-VF->getCoords()[XDIR]) + (VD->getCoords()[YDIR]-VF->getCoords()[YDIR])*(VD->getCoords()[YDIR]-VF->getCoords()[YDIR]);
////				distDF = sqrt(distDF);
////				double normDF[2] = {((VD->getCoords())[XDIR] - (VF->getCoords())[XDIR])/distDF,((VD->getCoords())[YDIR] - (VF->getCoords())[YDIR])/distDF};
////				mincos += abs(normCE[0]*normAB[0]+normCE[1]*normAB[1])*0.5;
////				mincos2 += abs(normDF[0]*normAB[0]+normDF[1]*normAB[1])*0.5;
//			}
      /////check which tri is more close to front
      int nummer = 0;
      int nummeralter = 0;
      for (int ii = 0; ii < 3; ii++) {
	if ((secondTri->getFace(ii)->isMergeable())) {
	  nummer++;
	}
	if ((alterTri->getFace(ii)->isMergeable())) {
	  nummeralter++;
	}
      }

      if (nummer == nummeralter) {
	if (mincos2 < mincos) {
	  if (mincos2 < refvalue) {
	    secondTri = alterTri;
	    mergeface = altermergeface;
	    if (checkOrient2D(v01, v11, v21) == checkOrient2D(v01, v21, v31)) {
	      Tri2 = true;
	    }
	  }
	  else {
	    MergeQ.vPopTopEntry();
	    continue;
	  }
	}
	else {
	  if (mincos < refvalue) {
	    if (checkOrient2D(v0, v1, v2) == checkOrient2D(v0, v2, v3)) {
	      Tri2 = true;
	    }
	  }
	  else {
	    MergeQ.vPopTopEntry();
	    continue;
	  }
	}
      }
      else {
	if (nummeralter < nummer) {
	  if (mincos2 < refvalue) {
	    secondTri = alterTri;
	    mergeface = altermergeface;
	    if (checkOrient2D(v01, v11, v21) == checkOrient2D(v01, v21, v31)) {
	      Tri2 = true;
	    }
	  }
	  else {
	    MergeQ.vPopTopEntry();
	    continue;
	  }
	}
	else {
	  if (mincos < refvalue) {
	    if (checkOrient2D(v0, v1, v2) == checkOrient2D(v0, v2, v3)) {
	      Tri2 = true;
	    }
	  }
	  else {
	    MergeQ.vPopTopEntry();
	    continue;
	  }
	}
      }

    }
    MergeQ.vPopTopEntry();
    if (Tri1 and Tri2) {
      int index = 0;

      for (int ii = 0; ii < 3; ii++) {
	if (startingTri->getFace(ii) != mergeface) {
	  startingTri->getFace(ii)->markAsNotMergeable();
	  QuadFaces[index] = startingTri->getFace(ii);
	  index++;
	}
	if (secondTri->getFace(ii) != mergeface) {
	  secondTri->getFace(ii)->markAsNotMergeable();
	  QuadFaces[index] = secondTri->getFace(ii);
	  index++;
	}
      }
      assert(index == 4);
      createQuadFromTris(startingTri, secondTri);
      mergenumber++;

    }
    else {
      continue;
    }
//		purgeAllEntities();
    for (int ii = 0; ii < 4; ii++) {
      assert(!QuadFaces[ii]->isDeleted());
      assert(!QuadFaces[ii]->isMergeable());
    }

//		/////
    if (mergenumber % 50 == 0) {
      std::cout << "Merged " << mergenumber << " Quads" << std::endl;
    }
//		MergeQ.clearQueue();
    Face* updatebaseedge = NULL;
    for (int ii = 0; ii < 4; ii++) {
      updatebaseedge = QuadFaces[ii];
      if (updatebaseedge->isDeleted() or (!updatebaseedge->isValid())) {
	continue;
      }
      if (updatebaseedge->isBdryFace()) {
	updatebaseedge->markAsNotMergeable();
	continue;
      }
      assert(
	  (updatebaseedge->getCell(0)->getType() == Cell::eQuadCell)
	      or (updatebaseedge->getCell(1)->getType() == Cell::eQuadCell));
      double edgelength = updatebaseedge->calcSize();
      bool base = (updatebaseedge->getCell(0)->getType() == Cell::eTriCell)
	  xor (updatebaseedge->getCell(1)->getType() == Cell::eTriCell);
      if (base) {
	TriCell* Tri = dynamic_cast<TriCell*>(
	    updatebaseedge->getCell(0)->getType() == Cell::eTriCell ?
		updatebaseedge->getCell(0) : updatebaseedge->getCell(1));
	int updatenummergeable = 0;
	for (int iii = 0; iii < 3; iii++) {
	  if (Tri->getFace(iii)->isMergeable()) {
	    assert(Tri->getFace(iii) != updatebaseedge);
	    updatenummergeable++;
	  }
	}
	assert(updatenummergeable < 3);
	int priority = 0;
	if (updatenummergeable == 0) {
	  priority = 0.;
	}
	else {
	  priority = -100. / updatenummergeable;
	}
//					if (updatenummergeable == 1){
//						for (int ii = 0; ii < 3; ii++){
//							if ((!Tri->getFace(ii)->isMergeable()) and Tri->getFace(ii) != updatebaseedge){
//								updatebaseedge = Tri->getFace(ii);
//							}
//						}
//					}
	MergeQ.qAddEntrywithPriority(
	    updatebaseedge, edgelength + priority + mergenumber * 1e-5);
      }

    }
//		std::cout << "Merged " << mergenumber << " Quads" << std::endl;
    if (mergenumber > 60000) {
      printf("exceed 60000\n");
      break;
    }
//		char strName[FILE_NAME_LEN];
//				sprintf(strName, "%s_%d",m_tempMeshFileName,m_tempMeshIndex);
//				writeVTKLegacyWithoutPurge(*this, strName);
//				m_tempMeshIndex++;
  }
  //
  std::cout << "Merge Game ends(hopefully it worked) :)" << std::endl;
  std::cout << "Merged " << mergenumber << " Quads" << std::endl;
}

double
Mesh2D::calcsquareness(Vert* const v0, Vert* const v1, Vert* const v2,
		       Vert* const v3)
{
//	double squareness;
  const double *cv0;
  const double *cv1;
  const double *cv2;
  const double *cv3;
  cv0 = v0->getCoords();
  cv1 = v1->getCoords();
  cv2 = v2->getCoords();
  cv3 = v3->getCoords();
  double l01 = (cv0[XDIR] - cv1[XDIR]) * (cv0[XDIR] - cv1[XDIR])
      + (cv0[YDIR] - cv1[YDIR]) * (cv0[YDIR] - cv1[YDIR]);
  double l12 = (cv2[XDIR] - cv1[XDIR]) * (cv2[XDIR] - cv1[XDIR])
      + (cv2[YDIR] - cv1[YDIR]) * (cv2[YDIR] - cv1[YDIR]);
  double l23 = (cv2[XDIR] - cv3[XDIR]) * (cv2[XDIR] - cv3[XDIR])
      + (cv2[YDIR] - cv3[YDIR]) * (cv2[YDIR] - cv3[YDIR]);
  double l30 = (cv0[XDIR] - cv3[XDIR]) * (cv0[XDIR] - cv3[XDIR])
      + (cv0[YDIR] - cv3[YDIR]) * (cv0[YDIR] - cv3[YDIR]);
  double l13 = (cv1[XDIR] - cv3[XDIR]) * (cv1[XDIR] - cv3[XDIR])
      + (cv1[YDIR] - cv3[YDIR]) * (cv1[YDIR] - cv3[YDIR]);
  double l02 = (cv2[XDIR] - cv0[XDIR]) * (cv2[XDIR] - cv0[XDIR])
      + (cv2[YDIR] - cv0[YDIR]) * (cv2[YDIR] - cv0[YDIR]);
  l01 = sqrt(l01);
  l12 = sqrt(l12);
  l23 = sqrt(l23);
  l30 = sqrt(l30);
  l13 = sqrt(l13);
  l02 = sqrt(l02);
  double cos0 = (l01 * l01 + l30 * l30 - l13 * l13) / (2 * l01 * l30);
  cos0 =
      abs(cos0) > abs((l01 * l01 + l12 * l12 - l02 * l02) / (2 * l01 * l12)) ?
	  abs(cos0) :
	  abs((l01 * l01 + l12 * l12 - l02 * l02) / (2 * l01 * l12));
  cos0 =
      abs(cos0) > abs((l23 * l23 + l12 * l12 - l13 * l13) / (2 * l23 * l12)) ?
	  abs(cos0) :
	  abs((l23 * l23 + l12 * l12 - l13 * l13) / (2 * l23 * l12));
  cos0 =
      abs(cos0) > abs((l23 * l23 + l30 * l30 - l02 * l02) / (2 * l23 * l30)) ?
	  abs(cos0) :
	  abs((l23 * l23 + l30 * l30 - l02 * l02) / (2 * l23 * l30));
  /////
//	///diagonal ratio
//	if (l13 > l02){
//		cos0 = l13/l02;
//	}
//	else{
//		cos0 = l02/l13;
//	}
  return cos0;
}

bool
Mesh2D::createQuadFromTris(Cell* const pCA, Cell* const pCB)
{
  assert(pCA->getType() == Cell::eTriCell);
  assert(pCB->getType() == Cell::eTriCell);

  Face * pF = findCommonFace(pCA, pCB);
  if (!pF) {
    logMessage(1, "Trying to create invalid quad cell, no common face\n");
    return false;
  }
  if (pCA->getRegion() != pCB->getRegion()) {
    logMessage(
	1,
	"Trying to create invalid quad cell, tri's have different regions\n");
    return false;
  }
  Vert * pVA = pF->getVert(0);
  Vert * pVC = pF->getVert(1);
  Vert * pVB = pCA->getOppositeVert(pF);
  Vert * pVD = pCB->getOppositeVert(pF);
  int iReg = pCA->getRegion();

  // clear out the cavity
  deleteCell(pCA);
  deleteCell(pCB);

  if (createQuadCell(pVA, pVB, pVC, pVD, iReg))
    return true;
  return false;
}

// tries to merge the cells
bool
Mesh2D::createQuadFromVerts(Vert* const pVA, Vert* const pVB, Vert* const pVC,
			    Vert* const pVD)
{
  Face * faces[4];
  Cell * cells[] =
    { NULL, NULL };
  // before even checking edges, check that these would form a convex cell
  double intersection[2];
  if (!calcLineIntersection(pVA, pVC, pVB, pVD, intersection)) {
    logMessage(1, "Trying to create invalid quad cell, would be non-convex\n");
    return false;
  }

  // make sure the edges exist
  recoverEdge(pVA, pVB, faces[0]);
  recoverEdge(pVB, pVC, faces[1]);
  recoverEdge(pVC, pVD, faces[2]);
  recoverEdge(pVD, pVA, faces[3]);
  GR_index_t numCommonCells = 0;

  for (GR_index_t i = 0; i < 4; i++) {
    Cell * cell = findCommonCell(faces[i], faces[(i + 1) % 4]);
    if (cell) {
      // lazy way to get it to put into the next cell
      cells[0 + (cells[0] != NULL)] = cell;
      numCommonCells++;
    }
  }

  // nothing we can do
  if (numCommonCells != 2) {
    logMessage(
	1,
	"Trying to create invalid quad cell, too many common cells between verts\n");
    return false;
  }
  if (cells[0]->getRegion() != cells[1]->getRegion()) {
    logMessage(
	1,
	"Trying to create invalid quad cell, tri's have different regions\n");
    return false;
  }
  int iReg = cells[0]->getRegion();
  //delete the cells
  deleteCell(cells[0]);
  deleteCell(cells[1]);
  // now that a cavity exists, create the cell
  if (createQuadCell(faces[0], faces[1], faces[2], faces[3], iReg))
    return true;
  return false;
}
bool
Mesh2D::createTrisFromQuad(Cell * const pC)
{
  double adDihed[4];
  int nDihed;
  assert(pC->getType() == Cell::eQuadCell);

  // figure out the dihedrals
  // don't split the smallest pair
  pC->calcAllDihed(adDihed, &nDihed);
  assert(nDihed == 4);

  Vert * pVA = pC->getVert(0);
  Vert * pVB = pC->getVert(1);
  Vert * pVC = pC->getVert(2);
  Vert * pVD = pC->getVert(3);

  int iReg = pC->getRegion();
  deleteCell(pC);
  if (adDihed[0] + adDihed[2] < adDihed[1] + adDihed[3]) {
    Cell * pCA = createTriCell(pVA, pVB, pVD, iReg);
    Cell * pCB = createTriCell(pVC, pVD, pVB, iReg);
    assert(pCA->doFullCheck());
    assert(pCB->doFullCheck());
  }
  else {
    Cell * pCA = createTriCell(pVB, pVC, pVA, iReg);
    Cell * pCB = createTriCell(pVD, pVA, pVC, iReg);
    assert(pCA->doFullCheck());
    assert(pCB->doFullCheck());
  }
  return true;
}
