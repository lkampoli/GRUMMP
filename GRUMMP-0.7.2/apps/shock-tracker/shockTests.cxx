/*
 * shockTests.cxx
 *
 *  Created on: Nov 3, 2017
 *      Author: cfog
 */

#define BOOST_TEST_MODULE ShockTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "GR_ShockTopo.h"

using GRUMMP::ShockTopo;
using GRUMMP::ShockEndPoint;
using GRUMMP::ShockCurve;
using GRUMMP::ShockSurface;

BOOST_AUTO_TEST_SUITE(ShockTopoTests)

  BOOST_AUTO_TEST_CASE(EndPoints)
  {
    int curves[] =
      { 3, 7, 4, 12 };
    ShockEndPoint SEP(33, 4, curves);

    BOOST_CHECK_EQUAL(33, SEP.getID());
    BOOST_CHECK_EQUAL(4, SEP.getNumCurves());
    BOOST_CHECK_EQUAL(3, SEP.getCurveID(0));
    BOOST_CHECK_EQUAL(7, SEP.getCurveID(1));
    BOOST_CHECK_EQUAL(4, SEP.getCurveID(2));
    BOOST_CHECK_EQUAL(12, SEP.getCurveID(3));
    BOOST_CHECK_EQUAL(8, ShockEndPoint::getMaxIncCurves());
  }

  BOOST_AUTO_TEST_CASE(Curves)
  {
    int surfs[] =
      { 1, 7, 42, 15, 13 };
    int points[] =
      { 22, 14 };
    // These constructors are here to confirm that the correct number of
    // surfaces actually works for each type of curve.
    ShockCurve SC1(24, ShockCurve::ShockPerimeter, 1, surfs, points);
    ShockCurve SC2(24, ShockCurve::ShockBdry, 1, surfs, points);
    ShockCurve SC3(24, ShockCurve::Reflection, 2, surfs, points);
    ShockCurve SC4(24, ShockCurve::TripleLine, 4, surfs, points);
    ShockCurve SC5(24, ShockCurve::QuadrupleLine, 5, surfs, points);

    // Now confirm data
    BOOST_CHECK_EQUAL(24, SC5.getID());
    BOOST_CHECK_EQUAL(22, SC5.getEndPointID(0));
    BOOST_CHECK_EQUAL(14, SC5.getEndPointID(1));
    BOOST_CHECK_EQUAL(5, SC5.getNumSurfaces());
    for (int ii = 0; ii < 5; ii++) {
      BOOST_CHECK_EQUAL(surfs[ii], SC5.getSurfaceID(ii));
    }
    BOOST_CHECK_EQUAL(8, ShockCurve::getMaxIncSurfaces());
    BOOST_CHECK_EQUAL(1,
		      ShockCurve::getNumReqSurfs(ShockCurve::ShockPerimeter));
    BOOST_CHECK_EQUAL(1, ShockCurve::getNumReqSurfs(ShockCurve::ShockBdry));
    BOOST_CHECK_EQUAL(2, ShockCurve::getNumReqSurfs(ShockCurve::Reflection));
    BOOST_CHECK_EQUAL(4, ShockCurve::getNumReqSurfs(ShockCurve::TripleLine));
    BOOST_CHECK_EQUAL(5, ShockCurve::getNumReqSurfs(ShockCurve::QuadrupleLine));
  }

  BOOST_AUTO_TEST_CASE(Surfaces)
  {
    int curves[] =
      { 5, 24, 48, 88, 4, 10, 14, 41 };
    ShockSurface SS(78, 8, curves);
    BOOST_CHECK_EQUAL(78, SS.getID());
    BOOST_CHECK_EQUAL(8, SS.getNumCurves());
    for (int ii = 0; ii < 8; ii++) {
      BOOST_CHECK_EQUAL(curves[ii], SS.getCurveID(ii));
    }
    BOOST_CHECK_EQUAL(20, ShockSurface::getMaxBoundingCurves());
  }

  BOOST_AUTO_TEST_CASE(ReadTopo)
  {
    ShockTopo ST("test-data/reflect.topo");

    BOOST_CHECK_EQUAL(6, ST.getNumEndPoints());
    BOOST_CHECK_EQUAL(7, ST.getNumCurves());
    BOOST_CHECK_EQUAL(2, ST.getNumSurfaces());
  }
  BOOST_AUTO_TEST_SUITE_END()

