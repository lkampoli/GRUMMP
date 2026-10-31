#define BOOST_TEST_MODULE BaseTests

#include "boost/test/unit_test.hpp"

#include "GR_Vertex.h"

BOOST_AUTO_TEST_SUITE(FaceListTests)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    FaceList FL;
    BOOST_CHECK(FL.getFace(0) == pFInvalidFace);
    BOOST_CHECK(FL.size() == 0);
  }

  BOOST_AUTO_TEST_CASE(Add)
  {
    FaceList FL;
    Face *pF0 = reinterpret_cast<Face*>(0x100);
    Face *pF1 = reinterpret_cast<Face*>(0x200);
    BOOST_CHECK(FL.addFace(pF0));
    BOOST_CHECK(FL.size() == 1);
    BOOST_CHECK(FL.hasFace(pF0));
    BOOST_CHECK(!FL.hasFace(pF1));

    BOOST_CHECK(!FL.addFace(pF0));
    BOOST_CHECK_EQUAL(FL.size(), 1);
    BOOST_CHECK(FL.hasFace(pF0));
    BOOST_CHECK(!FL.hasFace(pF1));

    BOOST_CHECK(FL.addFace(pF1));
    BOOST_CHECK_EQUAL(FL.size(), 2);
    BOOST_CHECK(FL.hasFace(pF0));
    BOOST_CHECK(FL.hasFace(pF1));
    BOOST_CHECK(!FL.addFace(pF1));

    BOOST_CHECK_EQUAL(FL.getFace(0), pF0);
    BOOST_CHECK_EQUAL(FL.getFace(1), pF1);
  }

  BOOST_AUTO_TEST_CASE(RemoveFirst)
  {
    FaceList FL;
    Face *pF0 = reinterpret_cast<Face*>(0x100);
    Face *pF1 = reinterpret_cast<Face*>(0x200);
    FL.addFace(pF0);
    FL.addFace(pF1);

    BOOST_CHECK(FL.removeFace(pF0));
    BOOST_CHECK_EQUAL(FL.size(), 1);
    BOOST_CHECK(!FL.hasFace(pF0));
    BOOST_CHECK(FL.hasFace(pF1));
    BOOST_CHECK_EQUAL(FL.getFace(0), pF1);
    BOOST_CHECK(!FL.removeFace(pF0));
  }

  BOOST_AUTO_TEST_CASE(RemoveSecond)
  {
    FaceList FL;
    Face *pF0 = reinterpret_cast<Face*>(0x100);
    Face *pF1 = reinterpret_cast<Face*>(0x200);
    FL.addFace(pF0);
    FL.addFace(pF1);

    BOOST_CHECK(FL.removeFace(pF1));
    BOOST_CHECK_EQUAL(FL.size(), 1);
    BOOST_CHECK(FL.hasFace(pF0));
    BOOST_CHECK(!FL.hasFace(pF1));
    BOOST_CHECK_EQUAL(FL.getFace(0), pF0);
    BOOST_CHECK(!FL.removeFace(pF1));
  }

  BOOST_AUTO_TEST_CASE(Clear)
  {
    FaceList FL;
    Face *pF0 = reinterpret_cast<Face*>(0x100);
    Face *pF1 = reinterpret_cast<Face*>(0x200);
    FL.addFace(pF0);
    FL.addFace(pF1);

    FL.clear();
    BOOST_CHECK_EQUAL(FL.size(), 0);
    BOOST_CHECK(!FL.hasFace(pF0));
    BOOST_CHECK(!FL.hasFace(pF1));
  }

  BOOST_AUTO_TEST_CASE(OpEquals)
  {
    FaceList FL, FL1;
    Face *pF0 = reinterpret_cast<Face*>(0x100);
    Face *pF1 = reinterpret_cast<Face*>(0x200);
    FL.addFace(pF0);
    FL.addFace(pF1);

    FL1 = FL;
    BOOST_CHECK_EQUAL(FL1.size(), 2);
    BOOST_CHECK(FL1.hasFace(pF0));
    BOOST_CHECK(FL1.hasFace(pF1));

    BOOST_CHECK_EQUAL(FL1.getFace(0), pF0);
    BOOST_CHECK_EQUAL(FL1.getFace(1), pF1);
  }

  BOOST_AUTO_TEST_SUITE_END()