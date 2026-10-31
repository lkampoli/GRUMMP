#include "GR_SurfMesh.h"

Vert*
SurfMesh::createVert(const double dx, const double dy, const double dz)
{
  double adData[] =
    { dx, dy, dz };
  return createVert(adData);
}

Vert*
SurfMesh::createVert(const double adCoords[])
{
  Vert *pV = getNewVert();
  pV->resetAllData();
  pV->setCoords(3, adCoords);
  return pV;
}

MultiEdge*
SurfMesh::createMultiEdge(Vert* const pV0, Vert* const pV1)
{
  assert(pV0->isValid());
  assert(pV1->isValid());
  // This thing had better not have existed before.
  assert(!findCommonFace(pV0, pV1)->isValid());
  MultiEdge *pME = getNewMultiEdge();
  pME->setVerts(pV0, pV1);
  pME->setLeftCell(NULL);
  pME->setRightCell(NULL);

#ifndef OMIT_VERTEX_HINTS
  // Make sure verts always have a valid hint face...
  if (!pV0->getHintFace()->isValid())
    pV0->setHintFace(pME);
  if (!pV1->getHintFace()->isValid())
    pV1->setHintFace(pME);
#endif

  assert(pME->doFullCheck());

  return pME;
}

MultiEdge*
SurfMesh::createMultiEdge(Face* const pF)
{
  assert(pF->isValid());
  assert(pF->getLeftCell()->isValid());
  assert(pF->getRightCell()->isValid());
  // This thing had better not have existed before.
  MultiEdge *pME = getNewMultiEdge();

  Vert *pV0 = pF->getVert(0);
  Vert *pV1 = pF->getVert(1);

  Cell *pCR = pF->getRightCell();
  Cell *pCL = pF->getLeftCell();

  pME->setVerts(pV0, pV1);
  pME->addCell(pCL);
  pME->addCell(pCR);

  pCL->replaceFace(pF, pME);
  pCR->replaceFace(pF, pME);

  pF->removeCell(pCL);
  pF->removeCell(pCR);
  deleteFace(pF);

#ifndef OMIT_VERTEX_HINTS
  // Make sure verts always have a valid hint face...
  if (!pV0->getHintFace()->isValid())
    pV0->setHintFace(pME);
  if (!pV1->getHintFace()->isValid())
    pV1->setHintFace(pME);
#endif

  assert(pME->doFullCheck());

  return pME;
}

TriCell*
SurfMesh::createTriCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
			const int iReg, const bool /*autoRotate*/, int /*ID*/)
{

  bool qExist;
  Face *pF01 = createFace(qExist, pV0, pV1, false);
  Face *pF12 = createFace(qExist, pV1, pV2, false);
  Face *pF20 = createFace(qExist, pV2, pV0, false);
  // These faces may or may not already exist.

  TriCell *pC = createTriCell(pF01, pF12, pF20, iReg);
  assert(pC->hasVert(pV0));
  assert(pC->hasVert(pV1));
  assert(pC->hasVert(pV2));
  return (pC);
}

static void
fixUpFaceCell(Face * pF, Cell *pC, Vert *pV, SurfMesh* pSM)
{
  if (pF->getRightCell()->isValid() && pF->getLeftCell()->isValid()) {
    if (pF->getType() != Face::eMultiEdge) {
      // Need to create a multi-edge.
      assert(pF->getRightCell() != pC && pF->getLeftCell() != pC);
      Face *pFOld = pF;
      pF = pSM->createMultiEdge(pF);
      pC->replaceFace(pFOld, pF);
    }
    pF->addCell(pC);
  }
  else if (pF->getRightCell()->isValid()) {
    assert(!pF->getLeftCell()->isValid());
    pF->setLeftCell(pC);
  }
  else if (pF->getLeftCell()->isValid()) {
    pF->setRightCell(pC);
  }
  else {
    if (pF->getVert(1) == pV)
      pF->setRightCell(pC);
    else
      pF->setLeftCell(pC);
  }
}

TriCell*
SurfMesh::createTriCell(Face * const pF0, Face * const pF1, Face * const pF2,
			const int iReg, int /*ID*/)
{
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
  apV[0] = findCommonVert(pF0, pF2);
  apV[1] = findCommonVert(pF0, pF1);
  apV[2] = findCommonVert(pF1, pF2);
  assert(apV[0]->isValid());
  assert(apV[1]->isValid());
  assert(apV[2]->isValid());

  TriCell *pC = dynamic_cast<TriCell*>(getNewCell(3));
  pC->resetAllData();
  pC->assignFaces(pF0, pF1, pF2);
  pC->setRegion(iReg);

  fixUpFaceCell(pF0, pC, apV[0], this);
  fixUpFaceCell(pF1, pC, apV[1], this);
  fixUpFaceCell(pF2, pC, apV[2], this);

  assert(pC->doFullCheck());
  assert(pC->hasFace(pF0));
  assert(pC->hasFace(pF1));
  assert(pC->hasFace(pF2));
  assert(pF0->doFullCheck());
  assert(pF1->doFullCheck());
  assert(pF2->doFullCheck());

  return pC;
}

QuadCell*
SurfMesh::createQuadCell(Face * const pF0, Face * const pF1, Face * const pF2,
			 Face * const pF3, const int iReg)
{
  // Faces must already exist
  assert(pF0->isValid());
  assert(pF1->isValid());
  assert(pF2->isValid());
  assert(pF3->isValid());

  // Set up the face->cell connectivity properly.  The canonical tri
  // looks like this:
  //
  //                3- - 2- - 2       
  //               /         /   .
  //              /         1    .
  //             3         /     .
  //            /         /      .
  //           0 - -0- - 1

  // Each pair of faces  must have a common vert.
  Vert *apV[4];
  apV[0] = findCommonVert(pF0, pF3);
  apV[1] = findCommonVert(pF0, pF1);
  apV[2] = findCommonVert(pF1, pF2);
  apV[3] = findCommonVert(pF2, pF3);
  assert(apV[0]->isValid());
  assert(apV[1]->isValid());
  assert(apV[2]->isValid());
  assert(apV[3]->isValid());

  QuadCell *pC = dynamic_cast<QuadCell*>(getNewCell(4));
  pC->resetAllData();
  pC->assignFaces(pF0, pF1, pF2, pF3);
  pC->setRegion(iReg);

  fixUpFaceCell(pF0, pC, apV[0], this);
  fixUpFaceCell(pF1, pC, apV[1], this);
  fixUpFaceCell(pF2, pC, apV[2], this);
  fixUpFaceCell(pF3, pC, apV[3], this);

  assert(pC->doFullCheck());
  assert(pC->hasFace(pF0));
  assert(pC->hasFace(pF1));
  assert(pC->hasFace(pF2));
  assert(pC->hasFace(pF3));
  assert(pF0->doFullCheck());
  assert(pF1->doFullCheck());
  assert(pF2->doFullCheck());
  assert(pF3->doFullCheck());

  return pC;
}

QuadCell*
SurfMesh::createQuadCell(Vert * const pV0, Vert * const pV1, Vert * const pV2,
			 Vert * const pV3, const int iReg)
{
  // Points had better be in cyclic order or bad things are gonna happen.
  bool existed;
  Face* pF01 = createFace(existed, pV0, pV1);
  Face* pF12 = createFace(existed, pV1, pV2);
  Face* pF23 = createFace(existed, pV2, pV3);
  Face* pF30 = createFace(existed, pV3, pV0);
  return createQuadCell(pF01, pF12, pF23, pF30, iReg);
}
