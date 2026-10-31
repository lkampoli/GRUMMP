#ifndef Mesh2DFixtureBase_h_
#define Mesh2DFixtureBase_h_

#include "GR_Geometry.h"
#include "GR_Mesh2D.h"

struct Mesh2DFixtureBase {
  Mesh2D M2D;
  Mesh2DFixtureBase() :
      M2D()
  {
    M2D.setEncroachmentType(eBall);
    iMessageStdoutLevel = 0;
  }
  ~Mesh2DFixtureBase()
  {
  }
protected:
  void
  checkCanonicalTri(Cell* const pCA, Vert* const pV0, Vert* const pV1,
		    Vert* const pV2)
  {
    Vert *pVA = pCA->getVert(0);
    Vert *pVB = pCA->getVert(1);
    Vert *pVC = pCA->getVert(2);
    BOOST_CHECK_EQUAL(1, checkOrient2D(pVA, pVB, pVC));

    Face *pF0 = pCA->getFace(0);
    Face *pF1 = pCA->getFace(1);
    Face *pF2 = pCA->getFace(2);
    BOOST_CHECK(!pF0->hasVert(pVC));
    BOOST_CHECK(!pF1->hasVert(pVA));
    BOOST_CHECK(!pF2->hasVert(pVB));

    if (pVA == pV0) {
      BOOST_CHECK_EQUAL(pV1, pVB);
      BOOST_CHECK_EQUAL(pV2, pVC);
    }
    else if (pVA == pV1) {
      BOOST_CHECK_EQUAL(pV2, pVB);
      BOOST_CHECK_EQUAL(pV0, pVC);
    }
    else {
      BOOST_CHECK_EQUAL(pV2, pVA);
      BOOST_CHECK_EQUAL(pV0, pVB);
      BOOST_CHECK_EQUAL(pV1, pVC);
    }

    BOOST_CHECK(
	(pF0->getLeftCell() == pCA
	    && (pF0->getVert(0) == pVA && pF0->getVert(1) == pVB))
	    || (pF0->getRightCell() == pCA
		&& (pF0->getVert(0) == pVB && pF0->getVert(1) == pVA)));
    BOOST_CHECK(
	(pF1->getLeftCell() == pCA
	    && (pF1->getVert(0) == pVB && pF1->getVert(1) == pVC))
	    || (pF1->getRightCell() == pCA
		&& (pF1->getVert(0) == pVC && pF1->getVert(1) == pVB)));
    BOOST_CHECK(
	(pF2->getLeftCell() == pCA
	    && (pF2->getVert(0) == pVC && pF2->getVert(1) == pVA))
	    || (pF2->getRightCell() == pCA
		&& (pF2->getVert(0) == pVA && pF2->getVert(1) == pVC)));
  }

  void
  checkCanonicalQuad(Cell* const pCA, Vert* const pV0, Vert* const pV1,
		     Vert* const pV2, Vert* const pV3)
  {
    Vert *pVA = pCA->getVert(0);
    Vert *pVB = pCA->getVert(1);
    Vert *pVC = pCA->getVert(2);
    Vert *pVD = pCA->getVert(3);
    // These tests can fail for a non-convex quad
    BOOST_CHECK_EQUAL(1, checkOrient2D(pVA, pVB, pVC));
    BOOST_CHECK_EQUAL(1, checkOrient2D(pVA, pVC, pVD));

    Face *pF0 = pCA->getFace(0);
    Face *pF1 = pCA->getFace(1);
    Face *pF2 = pCA->getFace(2);
    Face *pF3 = pCA->getFace(3);
    BOOST_CHECK(pF0->hasVert(pVA) && pF0->hasVert(pVB));
    BOOST_CHECK(pF1->hasVert(pVB) && pF1->hasVert(pVC));
    BOOST_CHECK(pF2->hasVert(pVC) && pF2->hasVert(pVD));
    BOOST_CHECK(pF3->hasVert(pVD) && pF3->hasVert(pVA));

    if (pVA == pV0) {
      BOOST_CHECK_EQUAL(pV1, pVB);
      BOOST_CHECK_EQUAL(pV2, pVC);
      BOOST_CHECK_EQUAL(pV3, pVD);
    }
    else if (pVA == pV1) {
      BOOST_CHECK_EQUAL(pV2, pVB);
      BOOST_CHECK_EQUAL(pV3, pVC);
      BOOST_CHECK_EQUAL(pV0, pVD);
    }
    else if (pVA == pV2) {
      BOOST_CHECK_EQUAL(pV3, pVB);
      BOOST_CHECK_EQUAL(pV0, pVC);
      BOOST_CHECK_EQUAL(pV1, pVD);
    }
    else {
      BOOST_CHECK_EQUAL(pV3, pVA);
      BOOST_CHECK_EQUAL(pV0, pVB);
      BOOST_CHECK_EQUAL(pV1, pVC);
      BOOST_CHECK_EQUAL(pV2, pVD);
    }

    BOOST_CHECK(
	(pF0->getLeftCell() == pCA
	    && (pF0->getVert(0) == pVA && pF0->getVert(1) == pVB))
	    || (pF0->getRightCell() == pCA
		&& (pF0->getVert(0) == pVB && pF0->getVert(1) == pVA)));
    BOOST_CHECK(
	(pF1->getLeftCell() == pCA
	    && (pF1->getVert(0) == pVB && pF1->getVert(1) == pVC))
	    || (pF1->getRightCell() == pCA
		&& (pF1->getVert(0) == pVC && pF1->getVert(1) == pVB)));
    BOOST_CHECK(
	(pF2->getLeftCell() == pCA
	    && (pF2->getVert(0) == pVC && pF2->getVert(1) == pVD))
	    || (pF2->getRightCell() == pCA
		&& (pF2->getVert(0) == pVD && pF2->getVert(1) == pVC)));
    BOOST_CHECK(
	(pF3->getLeftCell() == pCA
	    && (pF3->getVert(0) == pVD && pF3->getVert(1) == pVA))
	    || (pF3->getRightCell() == pCA
		&& (pF3->getVert(0) == pVA && pF3->getVert(1) == pVD)));
  }

  void
  checkMeshSize(const int nVerts, const int nFaces, const int nBdryFaces,
		const int nCells, const int nTris, const int nQuads)
  {
    BOOST_CHECK_EQUAL(nVerts, M2D.iNumVerts());
    BOOST_CHECK_EQUAL(nFaces, M2D.getNumEdgesInUse());
    BOOST_CHECK_EQUAL(nTris, M2D.getNumTrisInUse());
    BOOST_CHECK_EQUAL(nQuads, M2D.getNumQuadsInUse());
    BOOST_CHECK_EQUAL(nCells, M2D.getNumQuadsInUse() + M2D.getNumTrisInUse());
    BOOST_CHECK_EQUAL(nBdryFaces, M2D.getNumBdryEdgesInUse());
  }
};

#endif
