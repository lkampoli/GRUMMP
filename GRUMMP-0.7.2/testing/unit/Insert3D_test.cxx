#include <boost/test/unit_test.hpp>

#include "VolMeshFixtureBase.h"

#include "GR_InsertionManager.h"
#include "GR_SwapDecider.h"
#include "GR_BaseQueueEntry.h"

#include "CubitVector.hpp"

struct Refine3DFixture : public VolMeshFixtureBase {
  Refine3DFixture() :
      VolMeshFixtureBase()
  {
    // First, create some verts, forming an equilateral hex.
    checkMeshSize(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

    BOOST_CHECK(isValid());
    Vert *pV0 = createVert(1, 0, 0);
    Vert *pV1 = createVert(1, 1, 0);
    Vert *pV2 = createVert(0, 1, 0);
    Vert *pV3 = createVert(0, 0, 0);
    Vert *pV4 = createVert(1, 0, 1);
    Vert *pV5 = createVert(1, 1, 1);
    Vert *pV6 = createVert(0, 1, 1);
    Vert *pV7 = createVert(0, 0, 1);
//    checkMeshSize(8, 0, 0, 0, 0, 0, 0, 0, 0, 0);

// The six tet configuration is more convenient for testing, because it's
// easy to put a point onto the plane in the middle of the hex or on the
// long diagonal interior edge.

// First the prism at the x=1, y=1 corner.
    bool exists;
    (void) createTetCell(exists, pV0, pV1, pV2, pV5, 1);
//    checkMeshSize(8, 4, 0, 4, 0, 1, 1, 0, 0, 0);

    (void) createTetCell(exists, pV4, pV6, pV5, pV2, 1);
//    checkMeshSize(8, 8, 0, 8, 0, 2, 2, 0, 0, 0);

    (void) createTetCell(exists, pV0, pV2, pV4, pV5, 1);
//    checkMeshSize(8, 10, 0, 10, 0, 3, 3, 0, 0, 0);

// Now the other prism.
    (void) createTetCell(exists, pV0, pV2, pV3, pV7, 1);
//    checkMeshSize(8, 14, 0, 14, 0, 4, 4, 0, 0, 0);

    (void) createTetCell(exists, pV4, pV7, pV6, pV2, 1);
//    checkMeshSize(8, 17, 0, 17, 0, 5, 5, 0, 0, 0);

    (void) createTetCell(exists, pV0, pV4, pV2, pV7, 1);
//    checkMeshSize(8, 18, 0, 18, 0, 6, 6, 0, 0, 0);

// Faces at the top and bottom.
    Face *pF012 = findCommonFace(pV0, pV1, pV2);
    Face *pF230 = findCommonFace(pV2, pV3, pV0);
    Face *pF456 = findCommonFace(pV4, pV5, pV6);
    Face *pF467 = findCommonFace(pV4, pV6, pV7);

    // Faces at the sides
    Face *pF015 = findCommonFace(pV0, pV1, pV5);
    Face *pF054 = findCommonFace(pV0, pV5, pV4);
    Face *pF237 = findCommonFace(pV2, pV3, pV7);
    Face *pF267 = findCommonFace(pV2, pV6, pV7);

    // Faces at the front and back
    Face *pF265 = findCommonFace(pV2, pV5, pV6);
    Face *pF125 = findCommonFace(pV1, pV2, pV5);
    Face *pF037 = findCommonFace(pV0, pV3, pV7);
    Face *pF047 = findCommonFace(pV0, pV4, pV7);
//    checkMeshSize(8, 18, 0, 18, 0, 6, 6, 0, 0, 0);

    assert(pF012->isValid());
    assert(pF230->isValid());
    assert(pF456->isValid());
    assert(pF467->isValid());
    assert(pF015->isValid());
    assert(pF054->isValid());
    assert(pF237->isValid());
    assert(pF267->isValid());
    assert(pF265->isValid());
    assert(pF125->isValid());
    assert(pF037->isValid());
    assert(pF047->isValid());

    createBFace(pF012);
    createBFace(pF230);
    createBFace(pF456);
    createBFace(pF467);
    createBFace(pF015);
    createBFace(pF054);
    createBFace(pF237);
    createBFace(pF267);
    createBFace(pF265);
    createBFace(pF125);
    createBFace(pF037);
    createBFace(pF047);
    checkMeshSize(8, 18, 12, 18, 0, 6, 6, 0, 0, 0);
  }
  virtual
  ~Refine3DFixture()
  {
  }

//  void createInternalBdryEdge() {
//    Face* face02 = findCommonFace(pVVert(0), pVVert(2), true);
//    Face* face24 = findCommonFace(pVVert(2), pVVert(4), true);
//    Face* face40 = findCommonFace(pVVert(4), pVVert(0), true);
//    BOOST_CHECK(face02 != pFInvalidFace);
//    BOOST_CHECK(face24 != pFInvalidFace);
//    BOOST_CHECK(face40 != pFInvalidFace);
//
//    Cell* cell012 = pCCell(0);
//    cell012->vSetRegion(2);
//    Cell* cell024 = pCCommonCell(face02, face24);
//    int iReg = cell024->iRegion();
//    BOOST_CHECK(cell024 != pCInvalidCell);
//
//    deleteCell(cell024);
//    bool qExist;
//    Face* face02Dup = createFace(qExist, face02->pVVert(0),
//				     face02->pVVert(1), true);
//    BOOST_CHECK(qExist);
//    BOOST_CHECK(face02Dup != pFInvalidFace);
//    BOOST_CHECK_NE(face02Dup, face02);
//
//    cell024 = createTriCell(face02Dup, face24, face40, iReg);
//    IntBdryEdge* pIBE = createIntBdryEdge(face02, face02Dup);
//    BOOST_CHECK_NE(pIBE, pBFInvalidBFace);
//  }

  void
  postCheckSizesInterior()
  {
    BOOST_CHECK(isValid());
    checkMeshSize(9, 30, 12, 30, 0, 12, 12, 0, 0, 0);
  }

  void
  postCheckSizesBoundaryCurve()
  {
    BOOST_CHECK(isValid());
    checkMeshSize(9, 23, 14, 23, 0, 8, 8, 0, 0, 0);
  }
  void
  postCheckSizesBoundarySurface()
  {
    BOOST_CHECK(isValid());
    checkMeshSize(9, 27, 14, 27, 0, 10, 10, 0, 0, 0);
  }

//  void postCheckSizesIntBoundary() {
//    BOOST_CHECK(qValid());
//    checkMeshSize(7, 14, 6, 6, 6, 0);
//    BOOST_CHECK_EQUAL(2, iNumIntBdryEdgesInUse());
//  }

  void
  checkAllConnectedToNewVert(Vert *pV)
  {
    for (GR_index_t i = 0; i < getNumVerts() - 1; i++) {
      if (getVert(i)->getNumFaces() == 0)
	continue;
      Face *face = findCommonFace(pV, getVert(i), true);
      BOOST_CHECK_NE(face, pFInvalidFace);
    }
  }
};

//////////////////////////////////////////////////////////////////////
// The old insert-and-swap code is now obsolete, along with all these tests.
// Begin obsolete tests

// BOOST_FIXTURE_TEST_SUITE(Refine3D, Refine3DFixture)

// BOOST_AUTO_TEST_CASE(DirectInsertInCellNoSwap)
// {
//   // Testing: insert a point into a known cell, no swapping.
//   double point[] = {0.75, 0.75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell spans the vertex.
//   int swaps = iInsertInInterior(pV, pC, false);
//   checkMeshSize(9, 24, 12, 24, 0, 9, 9, 0, 0, 0);

//   BOOST_CHECK_EQUAL(0, swaps);

//   Face* pF012 = findCommonFace(pVVert(0), pVVert(1), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF012, pFInvalidFace);
//   Face* pF015 = findCommonFace(pVVert(0), pVVert(1), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF015, pFInvalidFace);
//   Face* pF025 = findCommonFace(pVVert(0), pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF025, pFInvalidFace);
//   Face* pF125 = findCommonFace(pVVert(1), pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF125, pFInvalidFace);

//   Face* pFN01 = findCommonFace(pV, pVVert(0), pVVert(1), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN01, pFInvalidFace);
//   Face* pFN02 = findCommonFace(pV, pVVert(0), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN02, pFInvalidFace);
//   Face* pFN05 = findCommonFace(pV, pVVert(0), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN05, pFInvalidFace);
//   Face* pFN12 = findCommonFace(pV, pVVert(1), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN12, pFInvalidFace);
//   Face* pFN15 = findCommonFace(pV, pVVert(1), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN15, pFInvalidFace);
//   Face* pFN25 = findCommonFace(pV, pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN25, pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertInCellNoSwap)
// {
//   // Testing: insert a point into a known cell, no swapping.
//   double point[] = {0.75, 0.75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(5); // This cell does not span the vertex.
//   int swaps = 0;
//   bool result = qInsertPoint(point, pC, &swaps, false, true, pV);
//   BOOST_CHECK(result);
//   checkMeshSize(9, 24, 12, 24, 0, 9, 9, 0, 0, 0);

//   BOOST_CHECK_EQUAL(0, swaps);

//   Face* pF012 = findCommonFace(pVVert(0), pVVert(1), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF012, pFInvalidFace);
//   Face* pF015 = findCommonFace(pVVert(0), pVVert(1), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF015, pFInvalidFace);
//   Face* pF025 = findCommonFace(pVVert(0), pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF025, pFInvalidFace);
//   Face* pF125 = findCommonFace(pVVert(1), pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pF125, pFInvalidFace);

//   Face* pFN01 = findCommonFace(pV, pVVert(0), pVVert(1), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN01, pFInvalidFace);
//   Face* pFN02 = findCommonFace(pV, pVVert(0), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN02, pFInvalidFace);
//   Face* pFN05 = findCommonFace(pV, pVVert(0), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN05, pFInvalidFace);
//   Face* pFN12 = findCommonFace(pV, pVVert(1), pVVert(2), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN12, pFInvalidFace);
//   Face* pFN15 = findCommonFace(pV, pVVert(1), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN15, pFInvalidFace);
//   Face* pFN25 = findCommonFace(pV, pVVert(2), pVVert(5), pVInvalidVert, true);
//   BOOST_CHECK_NE(pFN25, pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertInCellSwap)
// {
//   // Testing: insert a point into a known cell, with swapping.
//   double point[] = {0.75,0.75,0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell doesn't span the vertex.
//   int swaps;
//   bool result = qInsertPoint(point, pC, &swaps, true, true, pV);

//   BOOST_CHECK(result);
//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertInCellSwap)
// {
//   // Testing: insert a point into a cell, with swapping.  Must search for
//   // the cell.
//   double point[] = {0.75,0.75,0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell spans the vertex.
//   int swaps = iInsertInInterior(pV, pC, true);

//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesInterior();
//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnEdgeNoSwap)
// {
//   double point[] = {0.5, 0.5, 0.5};
//   Vert *pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Vert *pVA = pVVert(2);
//   Vert *pVB = pVVert(4);
//   Face *pF = findCommonFace(pVA, pVB, true);
//   BOOST_CHECK_NE(pF, pFInvalidFace);
//   int swaps = iInsertOnEdge(pV, pVA, pVB, pF, false);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 26, 12, 26, 0, 10, 10, 0, 0, 0);

//   for (int i = 0; i < 8; i++) {
//     if (i == 1 || i == 3) continue;
//     Face* face = findCommonFace(pV, pVVert(i), true);
//     BOOST_CHECK_NE(face, pFInvalidFace);
//   }
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnEdgeNoSwap)
// {
//   // Testing: insert a point on an edge, no swapping.  Must find the right cell.
//   double point[] = {0.5, 0.5, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell does not contain the point.
//   int swaps = 0;
//   bool result = qInsertPoint(point, pC, &swaps, false, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 26, 12, 26, 0, 10, 10, 0, 0, 0);

//   for (int i = 0; i < 8; i++) {
//     if (i == 1 || i == 3) continue;
//     Face* face = findCommonFace(pV, pVVert(i), true);
//     BOOST_CHECK_NE(face, pFInvalidFace);
//   }
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnEdgeSwap)
// {
//   // Testing: insert a point on a known edge, with swapping.
//   double point[] = {0.5, 0.5, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Vert *pVA = pVVert(2);
//   Vert *pVB = pVVert(4);
//   Face *pF = findCommonFace(pVA, pVB, true);
//   BOOST_CHECK_NE(pF, pFInvalidFace);
//   int swaps = iInsertOnEdge(pV, pVA, pVB, pF, true);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnEdgeSwap)
// {
//   // Testing: insert a point on an edge, with swapping.  Must search for
//   // the edge (and a containing cell.
//   double point[] = {0.5, 0.5, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell does not have the point on an edge.
//   int swaps;
//   bool result = qInsertPoint(point, pC, &swaps, true, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(2, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnFaceNoSwap)
// {
//   // Testing: insert a point on a known Face, no swapping.
//   double point[] = {0.75, 0.25, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(0), pVVert(2), pVVert(4),
//                             pVInvalidVert, true);
//   Cell *pC = pF->pCCellLeft(); // It's an interior face, so either one would do.
//   int swaps = iInsertOnFace(pV, pF, pC, false);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 26, 12, 26, 0, 10, 10, 0, 0, 0);

//   for (int i = 0; i < 8; i++) {
//     if (i == 1 || i == 3 || i == 6) continue;
//     Face* face = findCommonFace(pV, pVVert(i), true);
//     BOOST_CHECK_NE(face, pFInvalidFace);
//   }
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnFaceNoSwap)
// {
//   // Testing: insert a point on an Face, no swapping.  Must find the right cell.
//   double point[] = {0.75, 0.25, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell does not contain the point.
//   int swaps = 0;
//   bool result = qInsertPoint(point, pC, &swaps, false, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 26, 12, 26, 0, 10, 10, 0, 0, 0);

//   for (int i = 0; i < 8; i++) {
//     if (i == 1 || i == 3 || i == 6) continue;
//     Face* face = findCommonFace(pV, pVVert(i), true);
//     BOOST_CHECK_NE(face, pFInvalidFace);
//   }
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnFaceSwap)
// {
//   // Testing: insert a point on a known Face, with swapping.
//   double point[] = {0.75, 0.25, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(0), pVVert(2), pVVert(4),
//                             pVInvalidVert, true);
//   Cell *pC = pF->pCCellLeft(); // It's an interior face, so either one would do.
//   int swaps = iInsertOnFace(pV, pF, pC, true);

//   BOOST_CHECK_EQUAL(4, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnFaceSwap)
// {
//   // Testing: insert a point on an Face, with swapping.  Must search for
//   // the Face (and a containing cell.
//   // Testing: insert a point on an Face, no swapping.  Must find the right cell.
//   double point[] = {0.75, 0.25, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Cell* pC = pCCell(0); // This cell does not contain the point.
//   int swaps = 0;
//   bool result = qInsertPoint(point, pC, &swaps, true, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(4, swaps);
//   postCheckSizesInterior();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryEdgeNoSwap)
// {
//   // Testing: insert a point on a known bdry edge, no swapping.
//   double point[] = {1, 1, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(1), pVVert(5), true);
//   int swaps = iInsertOnEdge(pV, pVVert(1), pVVert(5), pF, false);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 21, 14, 21, 0, 7, 7, 0, 0, 0);

//   Face* pF025 = findCommonFace(pVVert(0), pVVert(2), pVVert(5),
//                                pVInvalidVert, true);
//   BOOST_CHECK_NE(pF025, pFInvalidFace);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryEdgeSwap)
// {
//   // Testing: insert a point on a known bdry edge, with swapping.
//   double point[] = {1, 1, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(1), pVVert(5), true);
//   int swaps = iInsertOnEdge(pV, pVVert(1), pVVert(5), pF, true);

//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesBoundary();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnBdryEdgeSwap)
// {
//   // Testing: insert a point on a known bdry edge, with swapping.
//   double point[] = {1, 1, 0.5};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   int swaps = 0;
//   bool result = qInsertPoint(point, pCCell(2), &swaps, true, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesBoundary();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryFaceSwap)
// {
//   // Testing: insert a point on a known bdry face, with swapping.
//   double point[] = {1, .75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(1), pVVert(5), pVVert(0),
//                             pVInvalidVert, true);
//   int swaps = iInsertOnFace(pV, pF, pCInvalidCell, true);

//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesBoundaryFace();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnBdryFaceSwap)
// {
//   // Testing: insert a point on a known bdry face, with swapping.
//   double point[] = {1, .75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   int swaps = 0;
//   bool result = qInsertPoint(point, pCCell(5), &swaps, true, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(5, swaps);
//   postCheckSizesBoundaryFace();

//   checkAllConnectedToNewVert(pV);
// }

// BOOST_AUTO_TEST_CASE(DirectInsertOnBdryFaceNoSwap)
// {
//   // Testing: insert a point on a known bdry face, with swapping.
//   double point[] = {1, .75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   Face *pF = findCommonFace(pVVert(1), pVVert(5), pVVert(0),
//                             pVInvalidVert, true);
//   int swaps = iInsertOnFace(pV, pF, pCInvalidCell, false);
//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 23, 14, 23, 0, 8, 8, 0, 0, 0);
// }

// BOOST_AUTO_TEST_CASE(SearchInsertOnBdryFaceNoSwap)
// {
//   // Testing: insert a point on a known bdry face, with swapping.
//   double point[] = {1, .75, 0.25};

//   Vert* pV = pVNewVert();
//   pV->vSetCoords(3, point);
//   int swaps = 0;
//   bool result = qInsertPoint(point, pCCell(5), &swaps, false, true, pV);
//   BOOST_CHECK(result);

//   BOOST_CHECK_EQUAL(0, swaps);
//   checkMeshSize(9, 23, 14, 23, 0, 8, 8, 0, 0, 0);
// }

// The old insert-and-swap code is now obsolete, along with all these tests.
// End obsolete tests
//////////////////////////////////////////////////////////////////////

// BOOST_AUTO_TEST_SUITE_END()

//BOOST_AUTO_TEST_CASE(ForceInsertWatson)
//{
//  CubitVector point(0, 0, 0);
//
//  VolMesh::WatsonData watson_data;
//  compute_watson_data(point, pCCell(0), watson_data);
//
//  BOOST_CHECK(watson_data.encroached_bdry_edges.empty());
//
//  Vert* new_vert = pVNewVert();
//  new_vert->vSetCoords(2, point);
//  new_vert->vSetType(Vert::eInterior);
//  new_vert->set_parent_entity(static_cast<GeometryEntity*>(NULL));
//
//  std::vector<Cell*> new_cells;
//
//  insert_watson_interior(new_vert, watson_data, &new_cells);
//
//  checkMeshSize(7, 12, 6, 6, 6, 0);
//
//  BOOST_CHECK(qValid());
//
//  for (GR_index_t i = 0; i < iNumCells(); i++) {
//    Cell *pC = pCCell(i);
//    if (pC->qDeleted() || !pC->qValid()) continue;
//    Vert *pVA = pC->pVVert(0);
//    Vert *pVB = pC->pVVert(1);
//    Vert *pVC = pC->pVVert(2);
//    checkCanonicalTri(pC, pVA, pVB, pVC);
//  }
//
//  BOOST_CHECK_EQUAL(3, pVVert(0)->iNumFaces());
//  BOOST_CHECK_EQUAL(3, pVVert(1)->iNumFaces());
//  BOOST_CHECK_EQUAL(3, pVVert(2)->iNumFaces());
//  BOOST_CHECK_EQUAL(3, pVVert(3)->iNumFaces());
//  BOOST_CHECK_EQUAL(3, pVVert(4)->iNumFaces());
//  BOOST_CHECK_EQUAL(3, pVVert(5)->iNumFaces());
//  BOOST_CHECK_EQUAL(6, pVVert(6)->iNumFaces());
//}
//
//// BOOST_AUTO_TEST_CASE(ForceBdryInsertWatson)
//// {
////   CubitVector point(0, sqrt(3)/2, 0);
//
////   VolMesh::WatsonData watson_data;
////   M2D.compute_watson_data(point, M2D.pCCell(0), watson_data, true, false);
//
////   BOOST_CHECK(!watson_data.encroached_bdry_edges.empty());
////   Face *pF = M2D.pFFaceBetween(M2D.pVVert(1), M2D.pVVert(2));
////   BdryEdgeBase *pBEB = dynamic_cast<BdryEdgeBase*>(pF->pCCellLeft());
////   if (!pBEB) pBEB = dynamic_cast<BdryEdgeBase*>(pF->pCCellRight());
//
////   Vert *new_vert = M2D.createVert(0, sqrt(3)/2);
////   new_vert->vSetType(Vert::eBdryCurve);
////   new_vert->set_parent_entity(NULL); // Should use the curve for pBEB,
//// 				     // but there isn't one.
//
////   std::vector<Cell*> new_cells;
//
////   M2D.insert_watson_boundary(new_vert, pBEB, watson_data,
//// 			     0.5, &new_cells);
//
////   checkMeshSize(7, 15, 8, 9, 9, 0);
//
////   BOOST_CHECK(M2D.qValid());
//
////   for (int i = 0; i < 9; i++) {
////     Cell *pC = M2D.pCCell(i);
////     if (pC->qDeleted() || !pC->qValid()) continue;
////     Vert *pVA = pC->pVVert(0);
////     Vert *pVB = pC->pVVert(1);
////     Vert *pVC = pC->pVVert(2);
////     checkCanonicalTri(pC, pVA, pVB, pVC);
////   }
//
////   BOOST_CHECK_EQUAL(3, M2D.pVVert(0)->iNumFaces());
////   BOOST_CHECK_EQUAL(2, M2D.pVVert(1)->iNumFaces());
////   BOOST_CHECK_EQUAL(2, M2D.pVVert(2)->iNumFaces());
////   BOOST_CHECK_EQUAL(3, M2D.pVVert(3)->iNumFaces());
////   BOOST_CHECK_EQUAL(3, M2D.pVVert(4)->iNumFaces());
////   BOOST_CHECK_EQUAL(3, M2D.pVVert(5)->iNumFaces());
////   BOOST_CHECK_EQUAL(6, M2D.pVVert(6)->iNumFaces());
//// }
//
//BOOST_AUTO_TEST_SUITE_END()
//
struct SwapInsert3DFixture : public Refine3DFixture {
  GRUMMP::DelaunaySwapDecider3D DSD;
  GRUMMP::SwapManager3D SwapMan;
  GRUMMP::SwappingInserter3D SI;
  SwapInsert3DFixture() :
      Refine3DFixture(), DSD(false), SwapMan(&DSD, this), SI(this, &SwapMan)
  {
  }
  virtual
  ~SwapInsert3DFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(SwapInsert3D, SwapInsert3DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    // Did we survive construction?
    BOOST_CHECK(true);
    // Testing: original mesh setup is correct.
    BOOST_CHECK(isValid());
    checkMeshSize(8, 18, 12, 18, 0, 6, 6, 0, 0, 0);

    BOOST_CHECK_EQUAL(9, getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(3, getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(12, getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(8, getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(7, getVert(5)->getNumFaces());
    BOOST_CHECK_EQUAL(5, getVert(6)->getNumFaces());
    BOOST_CHECK_EQUAL(7, getVert(7)->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(SI_InsertInCell)
  {
    // Testing: insert a point into a known cell, with swapping.
    double point[] =
      { 0.75, 0.75, 0.25 };

    Cell* pC = getCell(0); // This cell spans the vertex.
    Vert* newVert = SI.insertPointInCell(point, pC);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_EQUAL(1, SI.getNumInCell());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesInterior();

    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdge)
  {
    // Testing: insert a point on a known edge, with swapping.
    double point[] =
      { 0.5, 0.5, 0.5 };

    Vert* newVert = SI.insertPointOnEdge(point, getVert(2), getVert(4));
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdge)
  {
    // Testing: insert a point on a known edge, with swapping.
    double point[] =
      { 1, 1, 0.5 };

    Vert* newVert = SI.insertPointOnEdge(point, getVert(5), getVert(1));
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesBoundaryCurve();
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertInCellSearch)
  {
    // Testing: insert a point into a cell, with swapping.  Must search
    // for the cell
    double point[] =
      { 0.75, 0.75, 0.25 };

    Vert* newVert = SI.insertPoint(point, getCell(5));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumInCell());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdgeSearch)
  {
    // Testing: insert a point on an edge, with swapping.  Must search for
    // the cell.
    double point[] =
      { 0.5, 0.5, 0.5 };

    Vert* newVert = SI.insertPoint(point, getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnEdgeSearch_Snap)
  {
    // Testing: insert a point on an edge, with swapping.  Point gets
    // snapped to the edge midside.  Must search for
    // the cell.
    double point[] =
      { 0.5, 0.49, 0.53 };

    Vert* newVert = SI.insertPoint(point, getCell(0));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);

    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch)
  {
    // Testing: insert a point on a bdry edge, with swapping. Must search
    // for the cell.
    double point[] =
      { 1, 1, 0.5 };

    Vert* newVert = SI.insertPoint(point, getCell(3));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesBoundaryCurve();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch_Snap1)
  {
    // Testing: insert a point on a bdry edge, with swapping. Must search
    // for the cell.
    double point[] =
      { 1, 1, 0.48 };

    Vert* newVert = SI.insertPoint(point, getCell(3));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesBoundaryCurve();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch_Snap2)
  {
    // Testing: insert a point on a bdry edge, with swapping. Must search
    // for the cell.
    double point[] =
      { 1, 0.99, 0.48 };

    Vert* newVert = SI.insertPoint(point, getCell(3));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesBoundaryCurve();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnBdryEdgeSearch_Snap3)
  {
    // Testing: insert a point on a bdry edge, with swapping. Must search
    // for the cell.
    double point[] =
      { 0.99, 0.99, 0.48 };

    Vert* newVert = SI.insertPoint(point, getCell(3));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnBdryEdge());
    BOOST_CHECK_EQUAL(5, SI.getNumSwapsDone());
    postCheckSizesBoundaryCurve();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnFace)
  {
    // Testing: insert a point on a known Face, with swapping.
    bool prevVal = SI.setForcedInsertion(true);
    BOOST_CHECK(!prevVal);
    prevVal = SI.setForcedInsertion(true);
    BOOST_CHECK(prevVal);
    double point[] =
      { 0.75, 0.25, 0.22 };

    Face *face = findCommonFace(getVert(0), getVert(2), getVert(4),
    pVInvalidVert,
				true);
    Vert* newVert = SI.insertPointOnFace(point, face);
    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK_EQUAL(1, SI.getNumOnFace());
    BOOST_CHECK_EQUAL(4, SI.getNumSwapsDone());
    postCheckSizesInterior();
    Vert *pV = getVert(8);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(0.75, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.25, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.22, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnFaceSearch)
  {
    // Testing: insert a point on a Face, with swapping.  Must find the face
    // first.
    (void) SI.setForcedInsertion(true);
    double point[] =
      { 0.75, 0.25, 0.22 };
    Vert* newVert = SI.insertPoint(point, getCell(5));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnFace());
    BOOST_CHECK_EQUAL(4, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0.75, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.25, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.22, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnFaceSearch_Snap1)
  {
    // Testing: insert a point on a Face, with swapping.  Must find the face
    // first.
    double point[] =
      { 0.75, 0.25, 0.25 };
    Vert* newVert = SI.insertPoint(point, getCell(5));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnFaceSearch_Snap2)
  {
    // Testing: insert a point on a Face, with swapping.  Must find the face
    // first.
    double point[] =
      { 0.75, 0.23, 0.25 };
    Vert* newVert = SI.insertPoint(point, getCell(5));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(SI_InsertOnFaceSearch_Snap3)
  {
    // Testing: insert a point on a Face, with swapping.  Must find the face
    // first.
    double point[] =
      { 0.75, 0.25, 0.22 };
    Vert* newVert = SI.insertPoint(point, getCell(5));
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert, pV);
    BOOST_CHECK_EQUAL(1, SI.getNumOnEdge());
    BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
    postCheckSizesInterior();
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(0.5, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

//BOOST_AUTO_TEST_CASE(CreateInternalBdryEdge)
//{
//  createInternalBdryEdge();
//}
//
//BOOST_AUTO_TEST_CASE(SI_InsertOnIntBdryEdge)
//{
//  createInternalBdryEdge();
//  // Testing: insert a point on a known internal bdry edge, with swapping.
//  double point[] = {0.25,sqrt(3)/4};
//
//  Face* face = findCommonFace(pVVert(2), pVVert(0), true);
//  bool result = SI.insertPointOnInternalBdryFace(point, face);
//  BOOST_CHECK(result);
//  BOOST_CHECK_EQUAL(1, SI.getNumOnIntBdryEdge());
//  BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
//  postCheckSizesIntBoundary();
//  Vert *pV = pVVert(6);
//  checkAllConnectedToNewVert(pV);
//  BOOST_CHECK_CLOSE(0.25, pV->dX(), 1.e-8);
//  BOOST_CHECK_CLOSE(sqrt(3)/4, pV->dY(), 1.e-8);
//}
//
//BOOST_AUTO_TEST_CASE(SI_SearchInsertOnIntBdryEdge)
//{
//  createInternalBdryEdge();
//  // Testing: insert a point on an internal bdry edge, with swapping.
//  // Must search for the cell.
//  double point[] = {0.25,sqrt(3)/4};
//
//  bool result = SI.insertPoint(point, pCCell(0));
//  BOOST_CHECK(result);
//  BOOST_CHECK_EQUAL(1, SI.getNumOnIntBdryEdge());
//  BOOST_CHECK_EQUAL(2, SI.getNumSwapsDone());
//  postCheckSizesIntBoundary();
//  Vert *pV = pVVert(6);
//  checkAllConnectedToNewVert(pV);
//  BOOST_CHECK_CLOSE(0.25, pV->dX(), 1.e-8);
//  BOOST_CHECK_CLOSE(sqrt(3)/4, pV->dY(), 1.e-8);
//}
  BOOST_AUTO_TEST_CASE(Face_Queue_Entry)
  {
    std::map<FaceQueueEntry*, int, FaceQueueEntryComp> fqemap;
    std::map<FaceQueueEntry*, int, FaceQueueEntryComp>::iterator fqemap_it,
	fqemap_end;
// Populate the map
    FaceQueueEntry* fqe = nullptr;
    for (GR_index_t i = 0; i < getNumFaces(); i++) {
      fqe = new FaceQueueEntry(getFace(i));
//	printf("fqe %p %p\n",fqe,getFace(i));
      fqemap.insert(std::make_pair(fqe, i));
    }

// checking if the Face is in present, even if the fqe pointers are different
    Face * face = getFace(0);
    FaceQueueEntry * new_fqe1 = new FaceQueueEntry(face);
    BOOST_CHECK_EQUAL(1, fqemap.count(new_fqe1));
// Swap the first face's vert out with one it doesn't have
    Vert * replacement_vert;
    for (GR_index_t i = 0; i < getNumVerts(); i++) {
      replacement_vert = getVert(i);
      if (!face->hasVert(replacement_vert))
	break;
    }
    face->replaceVert(face->getVert(0), replacement_vert);
// create new entry based on this face,
    FaceQueueEntry * new_fqe2 = new FaceQueueEntry(face);
// the original entry is still present, though it is deleted
    BOOST_CHECK_EQUAL(1, fqemap.count(new_fqe1));
    BOOST_CHECK_EQUAL(1, new_fqe1->isDeleted());
// the new entry is not there, but is not deleted
//BOOST_CHECK_EQUAL(0,fqemap.count(new_fqe2));
    BOOST_CHECK_EQUAL(0, new_fqe2->isDeleted());
// I can insert it, though i don't like this
    fqemap.insert(std::make_pair(new_fqe2, 10));
    BOOST_CHECK_EQUAL(1, fqemap.count(new_fqe2));

    BOOST_CHECK(fqemap.count(fqe));
    delete new_fqe1;
    // delete new_fqe2;

    while (! (fqemap.empty())) {
    	fqemap_it = fqemap.begin();
    	FaceQueueEntry *fqe = fqemap_it->first;
    	delete fqe;
    	fqemap.erase(fqemap_it);
    }
  }
  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(CavityRefine3D, Refine3DFixture)

  BOOST_AUTO_TEST_CASE(CavityInsertAtHullPoint)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    Vert *targetVert = getVert(0);
    Vert *newVert = CI.insertPointInHull(targetVert);
    BOOST_CHECK_EQUAL(targetVert, newVert);
    BOOST_CHECK(isValid());
    checkMeshSize(8, 18, 12, 18, 0, 6, 6, 0, 0, 0);

    BOOST_CHECK_EQUAL(12, getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(3, getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(9, getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(5, getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(7, getVert(5)->getNumFaces());
    BOOST_CHECK_EQUAL(8, getVert(6)->getNumFaces());
    BOOST_CHECK_EQUAL(7, getVert(7)->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(CavityInsertInterior_InCell)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    // Testing: insert a point into a known cell using cavity insertion.
    double point[] =
      { 0.75, 0.75, 0.25 };

    Vert* newVert = CI.insertPointInHull(point, Vert::eInterior);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();

    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(CavityInsertInterior_OnEdge)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    // Testing: insert a point on a known edge using cavity inserter
    double point[] =
      { 0.5, 0.5, 0.5 };

    Vert* newVert = CI.insertPointInHull(point, Vert::eInterior);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_EQUAL(newVert->getVertType(), Vert::eInterior);
    postCheckSizesInterior();

    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(CavityInsertBdry_OnBdryEdge)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    Face *pF015 = findCommonFace(getVert(0), getVert(1), getVert(5),
    NULL,
				 true);
    Face *pF054 = findCommonFace(getVert(0), getVert(5), getVert(4),
    NULL,
				 true);

    BFace *pBF015 = dynamic_cast<BFace*>(pF015->getRightCell());
    if (pBF015 == pBFInvalidBFace) {
      pBF015 = dynamic_cast<BFace*>(pF015->getLeftCell());
    }
    BFace *pBF054 = dynamic_cast<BFace*>(pF054->getRightCell());
    if (pBF054 == pBFInvalidBFace) {
      pBF054 = dynamic_cast<BFace*>(pF054->getLeftCell());
    }

    CI.addBFaceToCavity(pBF015);
    CI.addBFaceToCavity(pBF054);

    // Testing: insert a point on a known edge using cavity inserter
    double point[] =
      { 1, 0.5, 0.5 };

    Vert* newVert = CI.insertPointInHull(point, Vert::eBdry);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_EQUAL(newVert->getVertType(), Vert::eBdry);
    postCheckSizesBoundarySurface();

    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(CavityInsertBdry_OnBdryFace)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    // Testing: insert a point on a known edge using cavity inserter
    double point[] =
      { 1, 0.5, 0.75 };

    Face *pF015 = findCommonFace(getVert(0), getVert(1), getVert(5),
    NULL,
				 true);
    Face *pF054 = findCommonFace(getVert(0), getVert(5), getVert(4),
    NULL,
				 true);

    BFace *pBF015 = dynamic_cast<BFace*>(pF015->getRightCell());
    if (pBF015 == pBFInvalidBFace) {
      pBF015 = dynamic_cast<BFace*>(pF015->getLeftCell());
    }
    BFace *pBF054 = dynamic_cast<BFace*>(pF054->getRightCell());
    if (pBF054 == pBFInvalidBFace) {
      pBF054 = dynamic_cast<BFace*>(pF054->getLeftCell());
    }

    CI.addBFaceToCavity(pBF015);
    CI.addBFaceToCavity(pBF054);
    Vert* newVert = CI.insertPointInHull(point, Vert::eBdry);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert->getVertType(), Vert::eBdry);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesBoundarySurface();

    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(0.5, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.75, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(CavityInsertBdry_OnBdryCurve)
  {
    GRUMMP::CavityInserter3D CI(this);
    for (GR_index_t i = 0; i < getNumCells(); i++) {
      CI.addCellToCavity(getCell(i));
    }
    // Testing: insert a point on a known edge using cavity inserter
    double point[] =
      { 1, 1, 0.75 };

    Face *pF015 = findCommonFace(getVert(0), getVert(1), getVert(5),
    NULL,
				 true);
    Face *pF054 = findCommonFace(getVert(0), getVert(5), getVert(4),
    NULL,
				 true);
    Face *pF265 = findCommonFace(getVert(2), getVert(6), getVert(5),
    NULL,
				 true);
    Face *pF125 = findCommonFace(getVert(1), getVert(2), getVert(5),
    NULL,
				 true);

    BFace *pBF015 = dynamic_cast<BFace*>(pF015->getRightCell());
    if (pBF015 == pBFInvalidBFace) {
      pBF015 = dynamic_cast<BFace*>(pF015->getLeftCell());
    }
    BFace *pBF054 = dynamic_cast<BFace*>(pF054->getRightCell());
    if (pBF054 == pBFInvalidBFace) {
      pBF054 = dynamic_cast<BFace*>(pF054->getLeftCell());
    }
    BFace *pBF265 = dynamic_cast<BFace*>(pF265->getRightCell());
    if (pBF265 == pBFInvalidBFace) {
      pBF265 = dynamic_cast<BFace*>(pF265->getLeftCell());
    }
    BFace *pBF125 = dynamic_cast<BFace*>(pF125->getRightCell());
    if (pBF125 == pBFInvalidBFace) {
      pBF125 = dynamic_cast<BFace*>(pF125->getLeftCell());
    }

    CI.addBFaceToCavity(pBF015);
    CI.addBFaceToCavity(pBF054);
    CI.addBFaceToCavity(pBF265);
    CI.addBFaceToCavity(pBF125);

    Vert* newVert = CI.insertPointInHull(point, Vert::eBdryCurve);
    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(newVert->getVertType(), Vert::eBdryCurve);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesBoundaryCurve();

    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(1, pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(1, pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(0.75, pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_SUITE_END()

struct WatsonInsert3DFixture : public Refine3DFixture {
  GRUMMP::WatsonInserter m_WI;
  WatsonInsert3DFixture() :
      Refine3DFixture(), m_WI(this)
  {
  }
  virtual
  ~WatsonInsert3DFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(WatsonINsert3D, WatsonInsert3DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    // Did we survive construction?
    BOOST_CHECK(true);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInCell)
  {
    // Testing: insert a point into a known cell via Watson
    double point[] =
      { 0.75, 0.75, 0.25 };

    Cell* pC = getCell(0); // This cell spans the vertex.
    m_WI.computeHull(point, pC);
    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(6, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInCell_NoSeed)
  {
    // Testing: insert a point via Watson without a seed cell given
    double point[] =
      { 0.75, 0.75, 0.25 };

    m_WI.computeHull(point, NULL);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertInBall)
  {
    // Testing: insert a point into a known cell via Watson
    double point[] =
      { 0.75, 0.75, 0.25 };

    Cell* pC = getCell(1); // This cell circumscribes the vertex.
    m_WI.computeHull(point, pC);
    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    postCheckSizesInterior();
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnEdge)
  {
    // Testing: insert a point on an edge of a known cell, via Watson.
    double point[] =
      { 0.5, 0.5, 0.5 };

    Cell* pC = getCell(2); // This cell has the vertex on an edge.
    m_WI.computeHull(point, pC);
    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(12, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();

    BOOST_CHECK(newVert->isValid());
    postCheckSizesInterior();
    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
    checkAllConnectedToNewVert(pV);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnBdryFace)
  {
    // Testing: insert a point on a known bdry edge, via Watson.
    double point[] =
      { 1, 0.75, 0.3 };

    // All cells have the same circumsphere, so any one will do.
    Cell* cell = getCell(0);
    m_WI.computeHull(point, cell);

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(6, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();
    BOOST_CHECK(newVert->isValid());
    postCheckSizesBoundarySurface();

    Vert *pV = getVert(8);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnBdryEdge)
  {
    // Testing: insert a point on a known bdry edge, via Watson.
    double point[] =
      { 1, 0.5, 0.5 };

    // All cells have the same circumsphere, so any one will do.
    Cell* cell = getCell(0);
    m_WI.computeHull(point, cell);

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(10, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();
    BOOST_CHECK(newVert->isValid());
    postCheckSizesBoundarySurface();

    Vert *pV = getVert(8);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnBdryCurve)
  {
    // Testing: insert a point on a known bdry edge, via Watson.
    double point[] =
      { 1, 1, 0.5 };

    // All cells have the same circumsphere, so any one will do.
    Cell* cell = getCell(0);
    m_WI.computeHull(point, cell);

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(4, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();
    BOOST_CHECK(newVert->isValid());
    postCheckSizesBoundaryCurve();

    Vert *pV = getVert(8);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_EQUAL(pV, newVert);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(WI_InsertOnBdryApex)
  {
    // Testing: insert a point on a known bdry edge, via Watson.
    double point[] =
      { 1, 1, 0 };

    // All cells have the same circumsphere, so any one will do.
    Cell* cell = getCell(0);
    m_WI.computeHull(point, cell);

    BOOST_CHECK(m_WI.areBFacesEncroached());
    std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> encBFaces =
	m_WI.getEncroachedBFaces();
    BOOST_CHECK_EQUAL(6, distance(encBFaces.first, encBFaces.second));

    Vert* newVert = m_WI.insertPointInHull();

    // This is the vert that was at that location.  It's now disconnected from
    // the mesh but not deleted.
    Vert* vert1 = getVert(1);
    BOOST_CHECK(!vert1->isDeleted());
    BOOST_CHECK(vert1->getNumFaces() == 0);

    BOOST_CHECK(newVert->isValid());
    BOOST_CHECK(newVert->getNumFaces() == 12);

    Vert *pV = getVert(8);
    BOOST_CHECK_EQUAL(pV, newVert);
    checkAllConnectedToNewVert(pV);
    BOOST_CHECK_CLOSE(point[0], pV->x(), 1.e-8);
    BOOST_CHECK_CLOSE(point[1], pV->y(), 1.e-8);
    BOOST_CHECK_CLOSE(point[2], pV->z(), 1.e-8);
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(FlakeTest, Refine3DFixture)

  BOOST_AUTO_TEST_CASE(CreateFlake)
  {
    Vert *pV0 = getVert(0);
    Vert *pV1 = getVert(1);
    Vert *pV2 = getVert(2);
    Vert *pV3 = getVert(3);
    Vert *pV4 = getVert(4);
    Vert *pV5 = getVert(5);
    Vert *pV6 = getVert(6);
    Vert *pV7 = getVert(7);

    Cell *pC0 = getCell(0);
    Cell *pC1 = getCell(1);
    Cell *pC2 = getCell(2);

    deleteCell(pC0);
    deleteCell(pC1);
    deleteCell(pC2);

    bool qExists = false;
    // The order of verts matters:  must be cyclic, and verts ABC and ACD must
    // be triangles that are still in the mesh.  In this case, 024 and 462 are
    // the tris we're after.  Cyclically, we need 0264, or a rotation or
    // reversal.  With this order, we've got 4026 and 402, 426.  So, check!
    FlakeCell *pFC4026 = createFlakeCell(qExists, pV4, pV0, pV2, pV6);
    BOOST_CHECK(pFC4026->isValid());

    // The next flake I need is 0154, with the diagonal running from 0 to 5,
    // so tris 015 and 054
    FlakeCell *pFC0154 = createFlakeCell(qExists, pV0, pV1, pV5, pV4);
    BOOST_CHECK(pFC0154->isValid());
    Cell* pCCheck = getCell(7);
    BOOST_CHECK(pCCheck = pFC0154);

    // Finally, create flake 1265, with the diagonal running from 2 to 5, so
    // tris 125 and 256
    FlakeCell *pFC2651 = createFlakeCell(qExists, pV2, pV6, pV5, pV1);
    BOOST_CHECK(pFC2651->isValid());

    PrismCell *pPrism012456 = createPrismCell(qExists, pV0, pV1, pV2, pV4,
						 pV5, pV6);
    BOOST_CHECK(pPrism012456->isValid());

    pCCheck = getCell(6);
    BOOST_CHECK(pCCheck == pPrism012456);
    pCCheck = getCell(7);
    BOOST_CHECK(pCCheck = pFC0154);

    PrismCell *pPrism023467 = replaceCellsWithPrism(pV0, pV2, pV3, pV4, pV6,
						       pV7);
    BOOST_CHECK(pPrism023467->isValid());
    BOOST_CHECK(isValid());
    checkMeshSize(8, 17, 12, 12, 5, 6, 0, 0, 2, 0, 4);
    purgeAllEntities();
    checkMeshSize(8, 17, 12, 12, 5, 6, 0, 0, 2, 0, 4);
  }

  BOOST_AUTO_TEST_SUITE_END()
