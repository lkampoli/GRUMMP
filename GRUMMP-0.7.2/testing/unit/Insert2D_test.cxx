#include <boost/test/unit_test.hpp>

#include "Mesh2DFixtureBase.h"

#include "GR_misc.h"
#include "GR_InsertionManager.h"
#include "GR_SwapDecider.h"

#include "CubitVector.hpp"

struct Refine2DFixture : public Mesh2DFixtureBase {
  Refine2DFixture() :
      Mesh2DFixtureBase()
  {
    // First, create some verts, forming an equilateral hex.
    BOOST_CHECK(M2D.isValid());
    static const double sqrt3over2 = sqrt(3.) / 2;
    Vert *pV0 = M2D.createVert(1, 0);
    Vert *pV1 = M2D.createVert(0.5, sqrt3over2);
    Vert *pV2 = M2D.createVert(-0.5, sqrt3over2);
    Vert *pV3 = M2D.createVert(-1, 0);
    Vert *pV4 = M2D.createVert(-0.5, -sqrt3over2);
    Vert *pV5 = M2D.createVert(0.5, -sqrt3over2);

    // This triangulation leaves a big equilateral triangle in the middle
    // of the mesh.
    (void) M2D.createTriCell(pV0, pV1, pV2, 1);
    (void) M2D.createTriCell(pV2, pV3, pV4, 1);
    (void) M2D.createTriCell(pV4, pV5, pV0, 1);
    (void) M2D.createTriCell(pV0, pV2, pV4, 1);

    Face *pF01 = findCommonFace(pV0, pV1);
    Face *pF12 = findCommonFace(pV1, pV2);
    Face *pF23 = findCommonFace(pV2, pV3);
    Face *pF34 = findCommonFace(pV3, pV4);
    Face *pF45 = findCommonFace(pV4, pV5);
    Face *pF50 = findCommonFace(pV5, pV0);

    M2D.createBFace(pF01);
    M2D.createBFace(pF12);
    M2D.createBFace(pF23);
    M2D.createBFace(pF34);
    M2D.createBFace(pF45);
    M2D.createBFace(pF50);

  }
  virtual
  ~Refine2DFixture()
  {
  }

  void
  createInternalBdryEdge()
  {
    Face* face02 = findCommonFace(M2D.getVert(0), M2D.getVert(2), true);
    Face* face24 = findCommonFace(M2D.getVert(2), M2D.getVert(4), true);
    Face* face40 = findCommonFace(M2D.getVert(4), M2D.getVert(0), true);
    BOOST_CHECK(face02 != pFInvalidFace);
    BOOST_CHECK(face24 != pFInvalidFace);
    BOOST_CHECK(face40 != pFInvalidFace);

    Cell* cell012 = M2D.getCell(0);
    cell012->setRegion(2);
    Cell* cell024 = findCommonCell(face02, face24);
    int iReg = cell024->getRegion();
    BOOST_CHECK(cell024 != pCInvalidCell);

    M2D.deleteCell(cell024);
    bool qExist;
    Face* face02Dup = M2D.createFace(qExist, face02->getVert(0),
				     face02->getVert(1), true);
    BOOST_CHECK(qExist);
    BOOST_CHECK(face02Dup != pFInvalidFace);
    BOOST_CHECK_NE(face02Dup, face02);

    cell024 = M2D.createTriCell(face02Dup, face24, face40, iReg);
    IntBdryEdge* pIBE = M2D.createIntBdryEdge(face02, face02Dup);
    BOOST_CHECK_NE(pIBE, pBFInvalidBFace);
  }

  void
  postCheckSizesInterior()
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(7, 12, 6, 6, 6, 0);
  }
  void
  postCheckSizesBoundary()
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(7, 11, 7, 5, 5, 0);
  }
  void
  postCheckSizesIntBoundary()
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(7, 14, 6, 6, 6, 0);
    BOOST_CHECK_EQUAL(2, M2D.getNumIntBdryEdgesInUse());
  }

  void
  checkAllConnectedToNewVert(Vert *pV)
  {
    for (int i = 0; i < 6; i++) {
      Face *face = findCommonFace(pV, M2D.getVert(i), true);
      BOOST_CHECK(face != pFInvalidFace);
    }
  }
};

BOOST_FIXTURE_TEST_SUITE(Refine2D, Refine2DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    // Testing: original mesh setup is correct.
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(6, 9, 6, 4, 4, 0);

    BOOST_CHECK_EQUAL(4, M2D.getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(4, M2D.getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(4, M2D.getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(5)->getNumFaces());
  }

//////////////////////////////////////////////////////////////////////
// The old insert-and-swap code is now obsolete, along with all these tests.
// Begin obsolete tests

// BOOST_AUTO_TEST_CASE(DirectInsertInCellNoSwap)
// {
//   // Testing: insert a point into a known cell, no swapping.
//   double point[] = {0.375,0.375*sqrt(3)};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell spans the vertex.
//   int swaps = M2D.iInsertInInterior(pV, pC, false);

//   BOOST_CHECK_EQUAL(0, swaps);

//   Face* pF02 = findCommonFace(M2D.pVVert(0), M2D.pVVert(2), true);
//   BOOST_CHECK(pF02 != pFInvalidFace);

//   Face* pFN0 = findCommonFace(pV, M2D.pVVert(0), true);
//   BOOST_CHECK(pFN0 != pFInvalidFace);
//   Face* pFN1 = findCommonFace(pV, M2D.pVVert(1), true);
//   BOOST_CHECK(pFN1 != pFInvalidFace);
//   Face* pFN2 = findCommonFace(pV, M2D.pVVert(2), true);
//   BOOST_CHECK(pFN2 != pFInvalidFace);

//   Cell* pCN01 = pCCommonCell(pFN0, pFN1);
//   BOOST_CHECK(pCN01 != pCInvalidCell);
//   Cell* pCN12 = pCCommonCell(pFN1, pFN2);
//   BOOST_CHECK(pCN12 != pCInvalidCell);
//   Cell* pCN20 = pCCommonCell(pFN2, pFN0);
//   BOOST_CHECK(pCN20 != pCInvalidCell);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertInCellNoSwap)
// {
//   // Testing: insert a point into a cell, no swapping.  Must search for
//   // the cell.
//   double point[] = {0.375,0.375*sqrt(3)};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(3); // This cell doesn't span the vertex.
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, false, true, pV);
//   BOOST_CHECK(result);
//   BOOST_CHECK_EQUAL(swaps, 0);

//   postCheckSizesInterior();

//   Face* pF02 = findCommonFace(M2D.pVVert(0), M2D.pVVert(2), true);
//   BOOST_CHECK(pF02 != pFInvalidFace);

//   Face* pFN0 = findCommonFace(pV, M2D.pVVert(0), true);
//   BOOST_CHECK(pFN0 != pFInvalidFace);
//   Face* pFN1 = findCommonFace(pV, M2D.pVVert(1), true);
//   BOOST_CHECK(pFN1 != pFInvalidFace);
//   Face* pFN2 = findCommonFace(pV, M2D.pVVert(2), true);
//   BOOST_CHECK(pFN2 != pFInvalidFace);

//   Cell* pCN01 = pCCommonCell(pFN0, pFN1);
//   BOOST_CHECK(pCN01 != pCInvalidCell);
//   Cell* pCN12 = pCCommonCell(pFN1, pFN2);
//   BOOST_CHECK(pCN12 != pCInvalidCell);
//   Cell* pCN20 = pCCommonCell(pFN2, pFN0);
//   BOOST_CHECK(pCN20 != pCInvalidCell);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertInCellSwap)
// {
//   // Testing: insert a point into a known cell, with swapping.
//   double point[] = {0,0};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell doesn't span the vertex.
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, true, true, pV);

//   BOOST_CHECK(result);
//   BOOST_CHECK_EQUAL(3, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertInCellSwap)
// {
//   // Testing: insert a point into a cell, with swapping.  Must search for
//   // the cell.
//   double point[] = {0,0};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(3); // This cell spans the vertex.
//   int swaps = M2D.iInsertInInterior(pV, pC, true);

//   BOOST_CHECK_EQUAL(3, swaps);
//   postCheckSizesInterior();
//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnEdgeNoSwap)
// {
//   // Testing: insert a point on a known edge, no swapping.
//   double point[] = {0.25, sqrt(3)/4};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell has the point on an edge.
//   Face *pF = findCommonFace(M2D.pVVert(0), M2D.pVVert(2), true);
//   int swaps = M2D.iInsertOnFace(pV, pF, pC, false, true);

//   BOOST_CHECK_EQUAL(0, swaps);
//   postCheckSizesInterior();

//   Face* pFN0 = findCommonFace(pV, M2D.pVVert(0), true);
//   BOOST_CHECK(pFN0 != pFInvalidFace);
//   Face* pFN1 = findCommonFace(pV, M2D.pVVert(1), true);
//   BOOST_CHECK(pFN1 != pFInvalidFace);
//   Face* pFN2 = findCommonFace(pV, M2D.pVVert(2), true);
//   BOOST_CHECK(pFN2 != pFInvalidFace);
//   Face* pFN4 = findCommonFace(pV, M2D.pVVert(4), true);
//   BOOST_CHECK(pFN4 != pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnEdgeNoSwap)
// {
//   // Testing: insert a point on an edge, no swapping.  Must search for
//   // the edge (and a containing cell.
//   double point[] = {0.25, sqrt(3)/4};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(1); // This cell does not have the point on an edge.
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, false, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(0, swaps);
//   postCheckSizesInterior();

//   Face* pFN0 = findCommonFace(pV, M2D.pVVert(0), true);
//   BOOST_CHECK(pFN0 != pFInvalidFace);
//   Face* pFN1 = findCommonFace(pV, M2D.pVVert(1), true);
//   BOOST_CHECK(pFN1 != pFInvalidFace);
//   Face* pFN2 = findCommonFace(pV, M2D.pVVert(2), true);
//   BOOST_CHECK(pFN2 != pFInvalidFace);
//   Face* pFN4 = findCommonFace(pV, M2D.pVVert(4), true);
//   BOOST_CHECK(pFN4 != pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnEdgeSwap)
// {
//   // Testing: insert a point on a known edge, with swapping.
//   double point[] = {0.25, sqrt(3)/4};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell has the point on an edge.
//   Face *pF = findCommonFace(M2D.pVVert(0), M2D.pVVert(2), true);
//   int swaps = M2D.iInsertOnFace(pV, pF, pC, true, true);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnEdgeSwap)
// {
//   // Testing: insert a point on an edge, with swapping.  Must search for
//   // the edge (and a containing cell.
//   double point[] = {0.25, sqrt(3)/4};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(1); // This cell does not have the point on an edge.
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, true, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryEdgeNoSwap)
// {
//   // Testing: insert a point on a known bdry edge, no swapping.
//   double point[] = {0, sqrt(3)/2};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell has the point on an edge.
//   Face *pF = findCommonFace(M2D.pVVert(1), M2D.pVVert(2), true);
//   int swaps = M2D.iInsertOnBdryFace(pV, pF, pC, false, false);

//   BOOST_CHECK_EQUAL(0, swaps);
//   postCheckSizesBoundary();

//   Face* pFN0 = findCommonFace(pV, M2D.pVVert(0), true);
//   BOOST_CHECK(pFN0 != pFInvalidFace);
//   Face* pFN1 = findCommonFace(pV, M2D.pVVert(1), true);
//   BOOST_CHECK(pFN1 != pFInvalidFace);
//   Face* pFN2 = findCommonFace(pV, M2D.pVVert(2), true);
//   BOOST_CHECK(pFN2 != pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryEdgeSwap)
// {
//   // Testing: insert a point on a known bdry edge, with swapping.
//   double point[] = {0, sqrt(3)/2};

//   Vert* pV = M2D.pVNewVert();
//   pV->vSetCoords(2, point);
//   Cell* pC = M2D.pCCell(0); // This cell has the point on an edge.
//   Face *pF = findCommonFace(M2D.pVVert(1), M2D.pVVert(2), true);
//   int swaps = M2D.iInsertOnBdryFace(pV, pF, pC, true, false);

//   BOOST_CHECK_EQUAL(3, swaps);
//   postCheckSizesBoundary();

//   checkAllConnectedToNewVert(pV);
// }

// The old insert-and-swap code is now obsolete, along with all these tests.
// End obsolete tests
//////////////////////////////////////////////////////////////////////

// The old code actually doesn't do snap to edge until the point fails
// the orientation test.  So these next two tests don't do what I had 
// expected.

// BOOST_AUTO_TEST_CASE(SnapInsertOnEdgeSwap)
// {
//   double point[] = {-0.5005, 0};
//   Cell* pC = M2D.pCCell(0);
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, true, false);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   Vert* pV = M2D.pVVert(6);
//   checkAllConnectedToNewVert(pV);
//   BOOST_CHECK_CLOSE(-0.5, pV->adCoords()[0], 1.e-8);
//   BOOST_CHECK_CLOSE(0, pV->adCoords()[1], 1.e-8);
// }

// BOOST_AUTO_TEST_CASE(SnapInsertOnEdgeSwapOffcenter)
// {
//   double point[] = {-0.500y5, 0};
//   Cell* pC = M2D.pCCell(0);
//   int swaps;
//   bool result = M2D.qInsertPoint(point, pC, &swaps, true, false);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   Vert* pV = M2D.pVVert(6);
//   checkAllConnectedToNewVert(pV);
//   BOOST_CHECK_CLOSE(-0.5, pV->adCoords()[0], 1.e-8);
//   BOOST_CHECK_CLOSE(0, pV->adCoords()[1], 1.e-8);
// }

  BOOST_AUTO_TEST_CASE(CavityInsertAtHullPoint)
  {
    GRUMMP::CavityInserter2D CI(&M2D);
    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      CI.addCellToCavity(M2D.getCell(i));
    }
    Vert* targetVert = M2D.getVert(0);
    Vert* newVert = CI.insertPointInHull(targetVert);
    BOOST_CHECK_EQUAL(targetVert, newVert);
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(6, 9, 6, 4, 4, 0);

    BOOST_CHECK_EQUAL(5, M2D.getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(5)->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(CavityInsertInterior)
  {
    double newPtLoc[] =
      { 0, 0 };
    GRUMMP::CavityInserter2D CI(&M2D);
    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      CI.addCellToCavity(M2D.getCell(i));
    }
    Vert* newVert = CI.insertPointInHull(newPtLoc, Vert::eInterior, -1);
    BOOST_CHECK_CLOSE(newPtLoc[0], newVert->x(), 1.e-8);
    BOOST_CHECK_CLOSE(newPtLoc[1], newVert->y(), 1.e-8);

    checkMeshSize(7, 12, 6, 6, 6, 0);

    BOOST_CHECK(M2D.isValid());

    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      Cell *pC = M2D.getCell(i);
      if (pC->isDeleted() || !pC->isValid())
	continue;
      BOOST_CHECK(pC->hasVert(newVert));
      BOOST_CHECK_EQUAL(pC->getRegion(), 1);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);
      checkCanonicalTri(pC, pVA, pVB, pVC);
    }

    BOOST_CHECK_EQUAL(3, M2D.getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(5)->getNumFaces());
    BOOST_CHECK_EQUAL(6, M2D.getVert(6)->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(CavityInsertOnBdry)
  {
    double newPtLoc[] =
      { 0, sqrt(3) / 2 };
    GRUMMP::CavityInserter2D CI(&M2D);
    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      CI.addCellToCavity(M2D.getCell(i));
    }
    CI.addBFaceToCavity(M2D.getBFace(1));
    Vert* newVert = CI.insertPointInHull(newPtLoc, Vert::eBdry, 0.5);
    BOOST_CHECK_CLOSE(newPtLoc[0], newVert->x(), 1.e-8);
    BOOST_CHECK_CLOSE(newPtLoc[1], newVert->y(), 1.e-8);

    checkMeshSize(7, 11, 7, 5, 5, 0);

    BOOST_CHECK(M2D.isValid());

    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      Cell *pC = M2D.getCell(i);
      if (pC->isDeleted() || !pC->isValid())
	continue;
      BOOST_CHECK(pC->hasVert(newVert));
      BOOST_CHECK_EQUAL(pC->getRegion(), 1);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);
      checkCanonicalTri(pC, pVA, pVB, pVC);
    }

    BOOST_CHECK_EQUAL(3, M2D.getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(5)->getNumFaces());
    BOOST_CHECK_EQUAL(6, M2D.getVert(6)->getNumFaces());

  }

  BOOST_AUTO_TEST_CASE(CavityInsertOnIntBdryEdge)
  {
    createInternalBdryEdge();
    // Testing: insert a point on a known internal bdry edge, with swapping.
    double newPtLoc[] =
      { 0.25, sqrt(3) / 4 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(0), true);

    GRUMMP::CavityInserter2D CI(&M2D);
    for (GR_index_t i = 0; i < M2D.getNumCells(); i++) {
      CI.addCellToCavity(M2D.getCell(i));
    }
    BFace* intBdryFace = dynamic_cast<IntBdryEdge*>(face->getLeftCell());
    if (!intBdryFace->isValid())
      intBdryFace = dynamic_cast<IntBdryEdge*>(face->getRightCell());
    CI.addBFaceToCavity(intBdryFace);
    Vert* newVert = CI.insertPointInHull(newPtLoc, Vert::eBdryTwoSide, 0.5);
    BOOST_CHECK_CLOSE(newPtLoc[0], newVert->x(), 1.e-8);
    BOOST_CHECK_CLOSE(newPtLoc[1], newVert->y(), 1.e-8);

    BOOST_CHECK(newVert->isValid());
    postCheckSizesIntBoundary();
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0.25, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(sqrt(3) / 4, pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_SUITE_END()

struct WatsonInsert2DFixture : public Refine2DFixture {
  GRUMMP::WatsonInserter m_WI;
  WatsonInsert2DFixture() :
      Refine2DFixture(), m_WI(&M2D)
  {
  }
  virtual
  ~WatsonInsert2DFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(WatsonInsert2D, WatsonInsert2DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    // Did we survive construction?
    BOOST_CHECK(true);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInCell)
  {
    // Testing: insert a point into a known cell via Watson
    double point[] =
      { 0, 0, 0 };

    Cell* pC = M2D.getCell(3); // This cell spans the vertex.
    m_WI.computeHull(point, pC);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInCell_NoSeed)
  {
    // Testing: insert a point via Watson without a seed cell given
    double point[] =
      { 0, 0, 0 };

    m_WI.computeHull(point, NULL);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInBall)
  {
    // Testing: insert a point into a known cell via Watson
    double point[] =
      { 0, 0, 0 };

    Cell* pC = M2D.getCell(2); // This cell circumscribes the vertex.
    m_WI.computeHull(point, pC);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnEdge)
  {
    // Testing: insert a point on an edge of a known cell, via Watson.
    double point[] =
      { -0.6, 0, 0 };

    Cell* pC = M2D.getCell(3); // This cell has the vertex on an edge.
    m_WI.computeHull(point, pC);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> iterPair =
	m_WI.getEncroachedBFaces();
    // There should be two faces encroached here.
    std::set<BFace*>::iterator begin = iterPair.first, end = iterPair.second;
    begin++;
    begin++;
    bool result = (begin == end);
    BOOST_CHECK(result);
    // BOOST_CHECK_EQUAL(begin, end);

    BOOST_CHECK(newVert->isValid());
    postCheckSizesInterior();
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnBdryEdge)
  {
    // Testing: insert a point on a known bdry edge, via Watson.
    double point[] =
      { 0, sqrt(3) / 2, 0 };

    Cell* cell = M2D.getCell(0);
    m_WI.computeHull(point, cell);

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> iterPair =
	m_WI.getEncroachedBFaces();
    // There should be one face encroached here.
    std::set<BFace*>::iterator begin = iterPair.first, end = iterPair.second;
    begin++;
    bool result = (begin == end);
    BOOST_CHECK(result);

    Vert* newVert = m_WI.insertPointInHull();
    BOOST_CHECK(newVert->isValid());
    postCheckSizesBoundary();
    Vert *pV = M2D.getVert(6);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_SUITE_END()

struct SwapInsert2DFixture : public Refine2DFixture {
  GRUMMP::DelaunaySwapDecider2D DSD;
  GRUMMP::SwapManager2D SwapMan;
  GRUMMP::SwappingInserter2D SI;

  SwapInsert2DFixture() :
      Refine2DFixture(), DSD(), SwapMan(&DSD, &M2D), SI(&M2D, &SwapMan)
  {
  }
  virtual
  ~SwapInsert2DFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(SwapInsert2D, SwapInsert2DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    // Did we survive construction?
    BOOST_CHECK(true);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertInCell)
  {
    // Testing: insert a point into a known cell, with swapping.
    double point[] =
      { 0, 0 };

    Cell* pC = M2D.getCell(3); // This cell spans the vertex.
    Vert* newVert = SI.insertPointInCell(point, pC);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_EQUAL(1, SI.getNumInCell());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdge)
  {
    // Testing: insert a point on a known edge, with swapping.
    double point[] =
      { -0.5, 0 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(4), true);
    Vert* newVert = SI.insertPointOnFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdge_Snap)
  {
    // Testing: insert a point on a known edge, with swapping.  Point gets
    // snapped to the edge midside.
    double point[] =
      { -0.5, 0.1 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(4), true);
    Vert* newVert = SI.insertPointOnFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(-0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0, pV->y(), 1.e-8);

    checkAllConnectedToNewVert(pV);
    // Because of snap-to-midside, this isn't the same as the requested point.
    BOOST_CHECK_CLOSE(-0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.0, pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdge)
  {
    // Testing: insert a point on a known edge, with swapping.
    double point[] =
      { 0, sqrt(3) / 2 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(1), true);
    Vert* newVert = SI.insertPointOnBdryFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesBoundary();
    Vert *pV = M2D.getVert(6);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(0, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(sqrt(3) / 2, pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdge_Snap)
  {
    // Testing: insert a point on a known edge, with swapping.
    double point[] =
      { 0.1, 1 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(1), true);
    Vert* newVert = SI.insertPointOnBdryFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesBoundary();
    Vert *pV = M2D.getVert(6);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(0, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(sqrt(3) / 2, pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertInCellSearch)
  {
    // Testing: insert a point into a cell, with swapping.  Must search
    // for the cell
    double point[] =
      { 0, 0 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumInCell());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesInterior();

    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdgeSearch)
  {
    // Testing: insert a point on an edge, with swapping.  Must search for
    // the cell.
    double point[] =
      { -0.5, 0 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdgeSearch_Snap)
  {
    // Testing: insert a point on an edge, with swapping.  Point gets
    // snapped to the edge midside.  Must search for
    // the cell.
    double point[] =
      { -0.5, 0.1 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(-0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0, pV->y(), 1.e-8);

    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch)
  {
    // Testing: insert a point on a bdry edge, with swapping. Must search
    // for the cell.
    double point[] =
      { 0, sqrt(3) / 2 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(3));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesBoundary();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch_Snap)
  {
    // Testing: insert a point on a bdry edge, with swapping.    Point
    // gets snapped to the edge midside.  Must search for the cell.
    double point[] =
      { 0.1, 1 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(3, SI.getNumSwapsDone());
    postCheckSizesBoundary();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(sqrt(3) / 2, pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(CreateInternalBdryEdge)
  {
    createInternalBdryEdge();
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnIntBdryEdge)
  {
    createInternalBdryEdge();
    // Testing: insert a point on a known internal bdry edge, with swapping.
    double point[] =
      { 0.25, sqrt(3) / 4 };

    Face* face = findCommonFace(M2D.getVert(2), M2D.getVert(0), true);
    Vert* newVert = SI.insertPointOnInternalBdryFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnIntBdryEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesIntBoundary();
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_SearchInsertOnIntBdryEdge)
  {
    createInternalBdryEdge();
    // Testing: insert a point on an internal bdry edge, with swapping.
    // Must search for the cell.
    double point[] =
      { 0.25, sqrt(3) / 4 };

    Vert* newVert = SI.insertPoint(point, M2D.getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = M2D.getVert(6);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnIntBdryEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesIntBoundary();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
  }
  BOOST_AUTO_TEST_CASE(IPC_GetOffCenterPoint)
  {
    GRUMMP::CircumcenterIPC IPC_CC;
    GRUMMP::OffcenterIPC2D IPC_OC;
    GRUMMP::InsertionPointCalculator::ePointType pointType;

    // testing, compute offcenter location.
    // equilateral triangle
    Vert *pVA = M2D.createVert(0, sqrt(3.) / 2.);
    Vert *pVB = M2D.createVert(-0.5, 0.);
    Vert *pVC = M2D.createVert(0.5, 0.);
    TriCell * tricell = M2D.createTriCell(pVA, pVB, pVC, 1);
    double adCC[2];
    double adOC[2];
    pointType = IPC_OC.calcInsertionPoint(tricell, adOC);
    BOOST_CHECK_EQUAL(pointType,
		      GRUMMP::InsertionPointCalculator::eCircumcenter);
    pointType = IPC_CC.calcInsertionPoint(tricell, adCC);
    BOOST_CHECK_EQUAL(pointType,
		      GRUMMP::InsertionPointCalculator::eCircumcenter);

    BOOST_CHECK_CLOSE(adOC[0], adCC[0], 1.e-8);
    BOOST_CHECK_CLOSE(adOC[1], adCC[1], 1.e-8);
    // shift top point up, making triangle much longer
    double adNewCoords[2] =
      { 0., 5.0 };
    pVA->setCoords(2, adNewCoords);

    IPC_OC.calcInsertionPoint(tricell, adOC);
    IPC_CC.calcInsertionPoint(tricell, adCC);
    BOOST_CHECK_CLOSE(adOC[0], 0.0, 1.e-8);
    BOOST_CHECK_CLOSE(adOC[1], 2.195533538112318, 1.e-8);
    // shift main point down
    adNewCoords[0] = 1.51;
    adNewCoords[1] = 0.1;
    pVA->setCoords(2, adNewCoords);
//	IPC_OC.calcInsertionPoint(tricell,adOC);
    IPC_CC.calcInsertionPoint(tricell, adCC);
    BOOST_CHECK_CLOSE(adOC[0], 0.0, 1.e-8);
    BOOST_CHECK_CLOSE(adOC[1], 2.195533538112318, 1.e-8);
    M2D.deleteCell(tricell);
    M2D.deleteVert(pVA);
    M2D.deleteVert(pVB);
    M2D.deleteVert(pVC);
  }

  BOOST_AUTO_TEST_CASE(OptimalSteinerInsertionPointCalculation)
  {
    // Testing: original mesh setup is correct.
    checkMeshSize(6, 9, 6, 4, 4, 0);

    GRUMMP::OptimalSteinerIPC2D IPC_OS;
    GRUMMP::OffcenterIPC2D IPC_OC;
    GRUMMP::InsertionPointCalculator::ePointType pointType;
    Cell * cell0 = M2D.getCell(0);
    double adPointLoc[] =
      { 0, 0 };
    pointType = IPC_OC.calcInsertionPoint(dynamic_cast<SimplexCell*>(cell0),
					  adPointLoc);
    assert(pointType == GRUMMP::InsertionPointCalculator::eCircumcenter);
    pointType = IPC_OS.calcInsertionPoint(dynamic_cast<SimplexCell*>(cell0),
					  adPointLoc);
    assert(pointType == GRUMMP::InsertionPointCalculator::eCircumcenter);

  }
  BOOST_AUTO_TEST_SUITE_END()
