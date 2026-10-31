#include "boost/test/unit_test.hpp"

#include <memory>

#include "GR_Length.h"
#include "GR_Mesh2D.h"

struct Length2DFixture {
  Mesh2D M2D;
  double eps = 0.01;
  Length2DFixture() :
      M2D()
  {
    BOOST_CHECK(M2D.isValid());
    Vert *pVA = M2D.createVert(0, 0);
    Vert *pVB = M2D.createVert(0, eps);
    Vert *pVC = M2D.createVert(0, 1);
    Vert *pVD = M2D.createVert(1, 1);
    Vert *pVE = M2D.createVert(1, 0);

    M2D.createTriCell(pVA, pVD, pVB);
    M2D.createTriCell(pVB, pVD, pVC);
    M2D.createTriCell(pVA, pVE, pVD);

    Face *pFAB = findCommonFace(pVA, pVB);
    Face *pFBC = findCommonFace(pVB, pVC);
    Face *pFCD = findCommonFace(pVC, pVD);
    Face *pFDE = findCommonFace(pVD, pVE);
    Face *pFEA = findCommonFace(pVE, pVA);

    M2D.createBFace(pFAB, 1);
    M2D.createBFace(pFBC, 1);
    M2D.createBFace(pFCD, 1);
    M2D.createBFace(pFDE, 1);
    M2D.createBFace(pFEA, 1);
    BOOST_CHECK(M2D.isValid());
  }
};

BOOST_FIXTURE_TEST_SUITE(Length2D, Length2DFixture)

  BOOST_AUTO_TEST_CASE(R1G1)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 1,
								       1);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), eps);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), eps);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), 1 - eps);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(), 1);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), 1);
  }

  BOOST_AUTO_TEST_CASE(R2G1)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 2,
								       1);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), (1 - eps) / 2);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(), 1. / 2);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), 1. / 2);
  }

  BOOST_AUTO_TEST_CASE(R1G2)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 1,
								       2);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), eps);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), eps);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), 0.99 / 2 + eps);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(),
		      eps + sqrt(1.9801) / 2);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), eps + 0.5);
  }

  BOOST_AUTO_TEST_CASE(R2G2)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 2,
								       2);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), (1 - eps) / 2);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(), 0.5);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), 0.5);
  }

  BOOST_AUTO_TEST_CASE(R2G4)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 2,
								       4);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), 0.99 / 4 + eps / 2);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(),
		      eps / 2 + sqrt(1.9801) / 4);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), eps / 2 + 0.25);
  }

  BOOST_AUTO_TEST_CASE(Floor)
  {
    std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 1,
								       1, 0.1);
    Mesh *temp = lenObject->getMesh();
    BOOST_CHECK_EQUAL(temp->getVert(0)->getLengthScale(), 0.1);
    BOOST_CHECK_EQUAL(temp->getVert(1)->getLengthScale(), 0.1);
    BOOST_CHECK_EQUAL(temp->getVert(2)->getLengthScale(), 0.99);
    BOOST_CHECK_EQUAL(temp->getVert(3)->getLengthScale(), 1);
    BOOST_CHECK_EQUAL(temp->getVert(4)->getLengthScale(), 1);
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(IMR6) {
  Mesh2D M2D;
  BOOST_CHECK(M2D.isValid());
  Vert *pV0 = M2D.createVert(0, 0);
  Vert *pV1 = M2D.createVert(0.08, 0);
  Vert *pV2 = M2D.createVert(0.25, 0.5);
  Vert *pV3 = M2D.createVert(0.75, 0.5);
  Vert *pV4 = M2D.createVert(1, 0);
  Vert *pV5 = M2D.createVert(2, 0);
  Vert *pV6 = M2D.createVert(2, 1);
  Vert *pV7 = M2D.createVert(1, 1);
  Vert *pV8 = M2D.createVert(0.75, 0.6);
  Vert *pV9 = M2D.createVert(0.5, 0.56);
  Vert *pV10 = M2D.createVert(0, 1);

  M2D.createTriCell(pV0, pV1, pV10);
  M2D.createTriCell(pV1, pV2, pV10);
  M2D.createTriCell(pV2, pV9, pV10);
  M2D.createTriCell(pV2, pV3, pV9);
  M2D.createTriCell(pV3, pV8, pV9);
  M2D.createTriCell(pV3, pV7, pV8);
  M2D.createTriCell(pV3, pV4, pV7);
  M2D.createTriCell(pV4, pV6, pV7);
  M2D.createTriCell(pV4, pV5, pV6);

  // Named for their originating verts
  Face *pF0 = findCommonFace(pV0, pV1);
  Face *pF1 = findCommonFace(pV1, pV2);
  Face *pF2 = findCommonFace(pV2, pV3);
  Face *pF3 = findCommonFace(pV3, pV4);
  Face *pF4 = findCommonFace(pV4, pV5);
  Face *pF5 = findCommonFace(pV5, pV6);
  Face *pF6 = findCommonFace(pV6, pV7);
  Face *pF7 = findCommonFace(pV7, pV8);
  Face *pF8 = findCommonFace(pV8, pV9);
  Face *pF9 = findCommonFace(pV9, pV10);
  Face *pF10 = findCommonFace(pV10, pV0);

  M2D.createBFace(pF0, 1);
  M2D.createBFace(pF1, 1);
  M2D.createBFace(pF2, 1);
  M2D.createBFace(pF3, 1);
  M2D.createBFace(pF4, 1);
  M2D.createBFace(pF5, 1);
  M2D.createBFace(pF6, 1);
  M2D.createBFace(pF7, 1);
  M2D.createBFace(pF8, 1);
  M2D.createBFace(pF9, 1);
  M2D.createBFace(pF10, 1);
  BOOST_CHECK(M2D.isValid());

  std::shared_ptr<GRUMMP::Length> lenObject = GRUMMP::Length::Create(&M2D, 1,
								     1);

  Mesh *temp = lenObject->getMesh();
  const double tol = 1.e-12;
  BOOST_CHECK_CLOSE(temp->getVert(0)->getLengthScale(), 0.08, tol);
  BOOST_CHECK_CLOSE(temp->getVert(1)->getLengthScale(), 0.08, tol);
  BOOST_CHECK_CLOSE(temp->getVert(2)->getLengthScale(), 7 / sqrt(1109), tol);
  BOOST_CHECK_CLOSE(temp->getVert(3)->getLengthScale(), 5. / 1282 * sqrt(641),
		    tol);
  BOOST_CHECK_CLOSE(temp->getVert(4)->getLengthScale(), sqrt(5) / 4, tol);
  BOOST_CHECK_CLOSE(temp->getVert(5)->getLengthScale(), 1, tol);
  BOOST_CHECK_CLOSE(temp->getVert(6)->getLengthScale(), 1, tol);
  BOOST_CHECK_CLOSE(temp->getVert(7)->getLengthScale(), sqrt(89) / 20, tol);
  BOOST_CHECK_CLOSE(temp->getVert(8)->getLengthScale(), 0.1, tol);
  BOOST_CHECK_CLOSE(temp->getVert(9)->getLengthScale(), 0.06, tol);
  BOOST_CHECK_CLOSE(temp->getVert(10)->getLengthScale(), sqrt(5) / 4, tol);
}
