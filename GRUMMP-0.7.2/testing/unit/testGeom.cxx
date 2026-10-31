#define BOOST_TEST_MODULE Geom_Tests

#include "CubitBox.hpp"
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "GR_config.h"
#include "GR_Geometry.h"
#include "GR_GRCurve.h"
#include "GR_GRCurveGeom.h"
#include "GR_GRPoint.h"
#include "GR_misc.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#define check_float_equal(a,b) {BOOST_CHECK(iFuzzyComp(a, b) == 0);}

/* Notes
 * none of the set/get z-coordinate are tested. Theres no need.. yet..
 * invalid input is not tested either, leading to a large number
 * of untested branches.
 *
 * Daniel W. Zaide, July 2013
 */

BOOST_AUTO_TEST_SUITE(Geom)

  BOOST_AUTO_TEST_CASE(TriDiag)
  {
    double** a2dLHS = new double*[5];
    for (int i = 0; i < 5; i++) {
      a2dLHS[i] = new double[3];
    }
    double* adRHS = new double[5];

    for (int i = 0; i < 5; i++) {
      a2dLHS[i][0] = static_cast<double>(i);
      a2dLHS[i][1] = static_cast<double>(i + 2);
      a2dLHS[i][2] = static_cast<double>(i + 1);
      adRHS[i] = static_cast<double>(i);
    }

    solveTriDiag(a2dLHS, adRHS, 5);
    check_float_equal(adRHS[0], -0.41176470588235);
    check_float_equal(adRHS[1], 0.82352941176471);
    check_float_equal(adRHS[2], -0.52941176470588);
    check_float_equal(adRHS[3], 0.82352941176471);
    check_float_equal(adRHS[4], 0.11764705882353);

    delete[] adRHS;
    for (int i = 0; i < 5; i++) {
      delete[] a2dLHS[i];
    }
    delete[] a2dLHS;

  }
  BOOST_AUTO_TEST_CASE(SolveQuadratic)
  {

    double dA, dB, dC, dRoot1, dRoot2;
    int iNRoots;

    //dA and dB equal to 0. No root should be found.
    dA = dB = 0.;
    dC = 10.;
    solveQuadratic(dA, dB, dC, &iNRoots, &dRoot1, &dRoot2);
    BOOST_CHECK_EQUAL(iNRoots, 0);

    //dA = 0, solution to a linear equation.
    dA = 0.;
    dB = 1.;
    dC = 1.5;
    solveQuadratic(dA, dB, dC, &iNRoots, &dRoot1, &dRoot2);
    BOOST_CHECK_EQUAL(iNRoots, 1);
    check_float_equal(dRoot1, -1.5);

    //Coefficients are such that quadratic equation does not have any real roots.
    dA = 4.;
    dB = 3.;
    dC = 1.;
    solveQuadratic(dA, dB, dC, &iNRoots, &dRoot1, &dRoot2);
    BOOST_CHECK_EQUAL(iNRoots, 0);

    //Coefficients are such that quadratic equation has only one root.
    dA = 4.;
    dB = 4.;
    dC = 1.;
    solveQuadratic(dA, dB, dC, &iNRoots, &dRoot1, &dRoot2);
    BOOST_CHECK_EQUAL(iNRoots, 1);
    check_float_equal(dRoot1, -0.5);

    //Coefficients are such that quadratic equation has two real roots.
    dA = 1.;
    dB = 4.;
    dC = -2.25;
    solveQuadratic(dA, dB, dC, &iNRoots, &dRoot1, &dRoot2);
    BOOST_CHECK_EQUAL(iNRoots, 2);
    check_float_equal(dRoot1, 0.5);
    check_float_equal(dRoot2, -9. / 2.);

  }
  BOOST_AUTO_TEST_CASE(SolveCubic)
  {

    /////////////////////////////////////
    //Testing the cubic equation solver//
    /////////////////////////////////////

    double dA, dB, dC, dD;
    double adRoots[3];
    int iNRoots;

    //dA, dB and dC are equal to 0. No root should be found.
    dA = dB = dC = 0.;
    dD = 10.;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 0);

    //dA and dB are equal to 0. Solution to linear equation.
    dA = dB = 0.;
    dC = 5.;
    dD = 10.;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 1);
    check_float_equal(adRoots[0], -2.);

    //dA = 0. Solution to a quadratic equation.
    dA = 0.;
    dB = 1.;
    dC = 4.;
    dD = -2.25;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 2);
    check_float_equal(adRoots[0], 0.5);
    check_float_equal(adRoots[1], -9. / 2.);

    //Coefficients are such that cubic equation has three different real roots.
    dA = 2.;
    dB = -4.;
    dC = -22.;
    dD = 24.;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 3);
    check_float_equal(adRoots[0], 4.);
    check_float_equal(adRoots[1], -3.);
    check_float_equal(adRoots[2], 1.);

    //Coefficients are such that cubic equation has three equal real roots.
    dA = 1.;
    dB = 6.;
    dC = 12.;
    dD = 8.;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 3);
    check_float_equal(adRoots[0], -2.);
    check_float_equal(adRoots[1], -2.);
    check_float_equal(adRoots[2], -2.);

    //Coefficients are such that cubic equation has only one real roots.
    dA = 3.;
    dB = -10.;
    dC = 14.;
    dD = 27.;
    solveCubic(dA, dB, dC, dD, &iNRoots, adRoots);
    BOOST_CHECK_EQUAL(iNRoots, 1);
    check_float_equal(adRoots[0], -1.);

  }

// As of Revision 1001, tests important things.
  BOOST_AUTO_TEST_CASE(TestLine)
  {

    CubitVector V1(0., 0., 0.);
    CubitVector V2(1., -2., 0.);

    GRLine* line = new GRLine(V1, V2);

    //Testing basic access function.
    BOOST_CHECK_EQUAL(line->topo_dim(), 1);
    BOOST_CHECK_EQUAL(line->geom_dim(), 2);
    BOOST_CHECK(line->get_curve_type() == GRCurveGeom::LINE);

    //Testing min and max param
    check_float_equal(line->min_param(), 0.);
    check_float_equal(line->max_param(), 1.);

    //Testing the bounding box
    CubitBox box = line->bounding_box();
    CubitVector mini = box.minimum();
    CubitVector maxi = box.maximum();
    check_float_equal(mini.x(), 0.);
    check_float_equal(mini.y(), -2.);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 1.);
    check_float_equal(maxi.y(), 0.);
    check_float_equal(maxi.z(), 0.);

    //Testing function returning parameter at coordinate.
    CubitVector coord(0.5, -1., 0.);
    double length = sqrt(0.5 * 0.5 + 1.);
    check_float_equal(line->param_at_coord(coord), 0.5);
    check_float_equal(line->param_at_arc_length(0., length), 0.5);

    //Testing function returning coord at parameter.
    line->coord_at_param(0.75, coord);
    check_float_equal(coord.x(), 0.75);
    check_float_equal(coord.y(), -2. * 0.75);
    check_float_equal(coord.z(), 0.);

    //Testing function returning coord at distance.
    line->coord_at_dist(length, true, coord);
    check_float_equal(coord.x(), 0.5);
    check_float_equal(coord.y(), -1.);
    check_float_equal(coord.z(), 0.);
    line->coord_at_dist(length, false, coord);
    check_float_equal(coord.x(), 0.5);
    check_float_equal(coord.y(), -1.);
    check_float_equal(coord.z(), 0.);
    line->coord_at_dist(sqrt(5.) * 2.0, false, coord);
    check_float_equal(coord.x(), 0.);
    check_float_equal(coord.y(), 0.);
    check_float_equal(coord.z(), 0.);
    line->coord_at_dist(sqrt(5.) * 2.0, true, coord);
    check_float_equal(coord.x(), 1.);
    check_float_equal(coord.y(), -2.0);
    check_float_equal(coord.z(), 0.);
    //Testing funtion returning coord at middle distance between points.
    line->coord_at_mid_dist(line->min_param(), line->max_param(), coord);
    check_float_equal(coord.x(), 0.5);
    check_float_equal(coord.y(), -1.);
    check_float_equal(coord.z(), 0.);
    CubitVector coord1;
    CubitVector coord2;
    line->coord_at_param(line->min_param(), coord1);
    line->coord_at_param(line->max_param(), coord2);
    //line->coord_at_mid_dist(coord1, coord2, coord);
    check_float_equal(coord.x(), 0.5);
    check_float_equal(coord.y(), -1.);
    check_float_equal(coord.z(), 0.);

    //Testing function returning closest point on curve.
    coord.set(5., -10., 0.);
    CubitVector closest;
    line->closest_coord_on_curve(coord, closest);
    check_float_equal(closest.x(), 1.);
    check_float_equal(closest.y(), -2.);
    check_float_equal(closest.z(), 0.);
    line->closest_coord_on_curve(coord, closest, line->min_param(), 0.9);
    line->coord_at_param(0.9, coord1);
    check_float_equal(closest.x(), coord1.x());
    check_float_equal(closest.y(), coord1.y());
    check_float_equal(closest.z(), 0.);

    coord.set(-5., 10., 0.);
    line->closest_coord_on_curve(coord, closest);
    check_float_equal(closest.x(), 0.);
    check_float_equal(closest.y(), 0.);
    check_float_equal(closest.z(), 0.);
    line->closest_coord_on_curve(coord, closest, 0.1, line->max_param());
    line->coord_at_param(0.1, coord1);
    check_float_equal(closest.x(), coord1.x());
    check_float_equal(closest.y(), coord1.y());
    check_float_equal(closest.z(), 0.);

    coord.set(0.5, -1., 0.);
    line->closest_coord_on_curve(coord, closest);
    check_float_equal(closest.x(), 0.5);
    check_float_equal(closest.y(), -1.);
    check_float_equal(closest.z(), 0.);

    coord.set(-0.39442719099992, -1.44721359549996, 0.);
    line->closest_coord_on_curve(coord, closest);
    check_float_equal(closest.x(), 0.5);
    check_float_equal(closest.y(), -1.);
    check_float_equal(closest.z(), 0.);
    check_float_equal(line->closest_param_on_curve(coord), 0.5);

    //Testing function that verify if a coord is on the curve.
    check_float_equal(line->coord_on_curve(closest), 1.);

    //Testing function that calculate first and second derivative with respect to param.
    CubitVector FD, SD;

    //	line->first_deriv(closest, FD);
    //	check_float_equal(FD.x(), 1.);
    //	check_float_equal(FD.y(), -2.);
    //	check_float_equal(FD.z(), 0.);
    line->first_deriv(0.5, FD);
    check_float_equal(FD.x(), 1.);
    check_float_equal(FD.y(), -2.);
    check_float_equal(FD.z(), 0.);

    //	line->second_deriv(closest, SD);
    //	check_float_equal(SD.x(), 0.);
    //	check_float_equal(SD.y(), 0.);
    //	check_float_equal(SD.z(), 0.);
    line->second_deriv(0.5, SD);
    check_float_equal(SD.x(), 0.);
    check_float_equal(SD.y(), 0.);
    check_float_equal(SD.z(), 0.);

    //Testing curvature query function
    check_float_equal(line->curvature(0.5), 0.);
    //	check_float_equal(line->curvature(closest), 0.);

    //Testing unit normal
    //	line->unit_normal(closest, coord);
    //	check_float_equal(coord.x(), 0.89442719099992);
    //	check_float_equal(coord.y(), 0.44721359549996);
    //	check_float_equal(coord.z(), 0.);
    line->unit_normal(0.5, coord);
    check_float_equal(coord.x(), 0.89442719099992);
    check_float_equal(coord.y(), 0.44721359549996);
    check_float_equal(coord.z(), 0.);

    //Testing unit tangent
    //	line->unit_tangent(closest, coord);
    //	check_float_equal(coord.x(), 0.44721359549996);
    //	check_float_equal(coord.y(), -0.89442719099992);
    //	check_float_equal(coord.z(), 0.);
    line->unit_tangent(0.5, coord);
    check_float_equal(coord.x(), 0.44721359549996);
    check_float_equal(coord.y(), -0.89442719099992);
    check_float_equal(coord.z(), 0.);

    //Testing arc length measurement.
    check_float_equal(line->arc_length(), 2.23606797749979);
    coord.set(0.25, -0.5, 0.);
    closest.set(0.75, -1.5, 0.);
    //	check_float_equal(line->arc_length(coord, closest), 1.11803398874989);
    check_float_equal(line->arc_length(0.25, 0.75), 1.11803398874989);

    //	line->center_point(coord1);
    //	check_float_equal(coord1.x(), 0.5);
    //	check_float_equal(coord1.y(), -1);
    //	check_float_equal(coord1.z(), 0.);
    coord1.set(-10, -10, 10);
    line->center_point(line->min_param(), line->max_param(), coord1);
    check_float_equal(coord1.x(), 0.5);
    check_float_equal(coord1.y(), -1);
    check_float_equal(coord1.z(), 0.);
    //	coord1.set(-10, -10, 10);
    //	line->center_point(V1, V2, coord1);
    //	check_float_equal(coord1.x(), 0.5);
    //	check_float_equal(coord1.y(), -1);
    //	check_float_equal(coord1.z(), 0.);

    //Testing total variation of the tangent
    //	check_float_equal(line->TVT(), 0.);
    //	check_float_equal(line->TVT(coord, closest), 0.);
    check_float_equal(line->TVT(0.25, 0.5), 0.);

    check_float_equal(line->mid_TVT(0.2, 0.4), 0.3);
    //	check_float_equal(line->mid_TVT(V1, V2), 0.5);

    //Testing coord access functions.
    check_float_equal(line->x(0), 0.);
    check_float_equal(line->x(1), 1.);
    check_float_equal(line->y(0), 0.);
    check_float_equal(line->y(1), -2.);
    check_float_equal(line->z(0), 0.);
    check_float_equal(line->z(1), 0.);

    line->get_coord(0, coord);
    check_float_equal(coord.x(), 0.);
    check_float_equal(coord.y(), 0.);
    check_float_equal(coord.z(), 0.);
    line->get_coord(1, coord);
    check_float_equal(coord.x(), 1.);
    check_float_equal(coord.y(), -2.);
    check_float_equal(coord.z(), 0.);

    check_float_equal(line->get_coord(0).x(), 0.);
    check_float_equal(line->get_coord(0).y(), 0.);
    check_float_equal(line->get_coord(0).z(), 0.);
    check_float_equal(line->get_coord(1).x(), 1.);
    check_float_equal(line->get_coord(1).y(), -2.);
    check_float_equal(line->get_coord(1).z(), 0.);

    //Testing functions to set points.
    line->set_x(0, 1.);
    line->set_y(0, 2.);
    line->set_x(1, 3.);
    line->set_y(1, 4.);
    check_float_equal(line->x(0), 1.);
    check_float_equal(line->x(1), 3.);
    check_float_equal(line->y(0), 2.);
    check_float_equal(line->y(1), 4.);
    check_float_equal(line->z(0), 0.);
    check_float_equal(line->z(1), 0.);

    line->set_coord(0, V1);
    line->set_coord(1, V2);
    check_float_equal(line->x(0), 0.);
    check_float_equal(line->x(1), 1.);
    check_float_equal(line->y(0), 0.);
    check_float_equal(line->y(1), -2.);
    check_float_equal(line->z(0), 0.);
    check_float_equal(line->z(1), 0.);

    //Testing copy constructor and operator=
    GRLine* new_line = new GRLine(*line);
    GRLine other_line;
    other_line = *line;

    check_float_equal(new_line->x(0), 0.);
    check_float_equal(new_line->x(1), 1.);
    check_float_equal(new_line->y(0), 0.);
    check_float_equal(new_line->y(1), -2.);
    check_float_equal(new_line->z(0), 0.);
    check_float_equal(new_line->z(1), 0.);
    check_float_equal(other_line.x(0), 0.);
    check_float_equal(other_line.x(1), 1.);
    check_float_equal(other_line.y(0), 0.);
    check_float_equal(other_line.y(1), -2.);
    check_float_equal(other_line.z(0), 0.);
    check_float_equal(other_line.z(1), 0.);

    delete line;
    delete new_line;

  }

  BOOST_AUTO_TEST_CASE(TestArc)
  {
    CubitVector v1(1., 1., 0.);
    CubitVector v2(1. + 0.5 * sqrt(3.), 1.5, 0.);
    CubitVector v3(1.5, 1. + 0.5 * sqrt(3.), 0.);

    GRArc* arc1 = new GRArc(v1, v2, v3);
    GRArc* arc1r = new GRArc(v2, v3, 1., true);

    v2.set(1., 0., 0.);

    GRArc* arc2 = new GRArc(v1, v2, v3);
    GRArc* arc2r = new GRArc(v2, v3, 1., true);

    v2.set(1 + 0.5 * sqrt(3.), 1.5, 0.);
    v3.set(1., 0., 0.);

    GRArc* arc3 = new GRArc(v1, v2, v3);
    GRArc* arc3r = new GRArc(v2, v3, 1., false);

    v2.set(1., 0., 0.);
    v3.set(0., 1., 0.);

    GRArc* arc4 = new GRArc(v1, v2, v3);
    GRArc* arc4r = new GRArc(v2, v3, 1., false);

    //Basic queries.
    BOOST_CHECK_EQUAL(arc1->topo_dim(), 1);
    BOOST_CHECK_EQUAL(arc1r->topo_dim(), 1);
    BOOST_CHECK_EQUAL(arc1->geom_dim(), 2);
    BOOST_CHECK_EQUAL(arc1r->geom_dim(), 2);
    BOOST_CHECK(arc1->get_curve_type() == GRCurveGeom::ARC);
    BOOST_CHECK(arc1r->get_curve_type() == GRCurveGeom::ARC);

    //Testing the bounding box.
    CubitBox box = arc1->bounding_box();
    CubitVector mini = box.minimum();
    CubitVector maxi = box.maximum();
    check_float_equal(mini.x(), 1.5);
    check_float_equal(mini.y(), 1.5);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 1. + 0.5 * sqrt(3.));
    check_float_equal(maxi.y(), 1. + 0.5 * sqrt(3.));
    check_float_equal(maxi.z(), 0.);

    box = arc2->bounding_box();
    mini = box.minimum();
    maxi = box.maximum();
    check_float_equal(mini.x(), 1.);
    check_float_equal(mini.y(), 0.);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 2.);
    check_float_equal(maxi.y(), 1. + 0.5 * sqrt(3.));
    check_float_equal(maxi.z(), 0.);

    box = arc3->bounding_box();
    mini = box.minimum();
    maxi = box.maximum();
    check_float_equal(mini.x(), 0.);
    check_float_equal(mini.y(), 0.);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 1. + 0.5 * sqrt(3.));
    check_float_equal(maxi.y(), 2.);
    check_float_equal(maxi.z(), 0.);

    //Testing dParamAtCoord
    v1.set(1. + cos(0.25 * M_PI), 1. + sin(0.25 * M_PI), 0.);
    v2.set(1. + 0.5 * sqrt(3.), 1.5, 0.);
    v3.set(1.5, 1. + 0.5 * sqrt(3), 0.);

    check_float_equal(arc1->param_at_coord(v1), 0.5);
    check_float_equal(arc1r->param_at_coord(v1), 0.5);
    check_float_equal(arc1->param_at_coord(v2), 0.);
    check_float_equal(arc1r->param_at_coord(v2), 0.);
    check_float_equal(arc1->param_at_coord(v3), 1.);
    check_float_equal(arc1r->param_at_coord(v3), 1.);

    v1.set(1. + cos(0.25 * 7. * M_PI), 1. + sin(0.25 * 7. * M_PI), 0.);
    v2.set(1. + cos(0.25 * M_PI), 1. + sin(0.25 * M_PI), 0.);

    check_float_equal(arc2->param_at_coord(v1), 0.3);
    check_float_equal(arc2r->param_at_coord(v1), 0.3);
    check_float_equal(arc2->param_at_coord(v2), 0.9);
    check_float_equal(arc2r->param_at_coord(v2), 0.9);

    check_float_equal(arc3->param_at_coord(v2), 0.0625);
    check_float_equal(arc3r->param_at_coord(v2), 0.0625);

    check_float_equal(arc4->param_at_coord(v1), 1. / 6.);
    check_float_equal(arc4r->param_at_coord(v1), 1. / 6.);
    check_float_equal(arc4->param_at_coord(v2), 0.5);
    check_float_equal(arc4r->param_at_coord(v2), 0.5);

    //Testing param at given arc length
    check_float_equal(arc2->param_at_arc_length(arc2->min_param(), M_PI/6.),
		      0.2);
    check_float_equal(arc2->param_at_arc_length(0.2, M_PI/6.), 0.4);

    //Testing coord_at_param
    arc1->coord_at_param(0.5, v1);
    arc1r->coord_at_param(0.5, v2);
    check_float_equal(v1.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v1.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(v1.z(), 0.);
    check_float_equal(v2.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v2.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(v2.z(), 0.);

    arc2->coord_at_param(0.3, v1);
    arc2r->coord_at_param(0.9, v2);
    check_float_equal(v1.x(), 1. + cos(0.25 * 7. * M_PI));
    check_float_equal(v1.y(), 1. + sin(0.25 * 7. * M_PI));
    check_float_equal(v2.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v2.y(), 1. + sin(0.25 * M_PI));

    arc3->coord_at_param(0.0625, v1);
    arc3r->coord_at_param(0.0625, v2);
    check_float_equal(v1.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v1.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(v2.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v2.y(), 1. + sin(0.25 * M_PI));

    arc4->coord_at_param(1. / 6., v1);
    arc4r->coord_at_param(0.5, v2);
    check_float_equal(v1.x(), 1. + cos(0.25 * 7. * M_PI));
    check_float_equal(v1.y(), 1. + sin(0.25 * 7. * M_PI));
    check_float_equal(v2.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v2.y(), 1. + sin(0.25 * M_PI));

    //Testing coord_at_dist
    arc2->coord_at_dist(0.2, true, v1);
    arc2r->coord_at_dist(0.2, true, v2);
    check_float_equal(v1.x(), 1.198997487421);
    check_float_equal(v1.y(), 0.02);
    check_float_equal(v2.x(), 1.198997487421);
    check_float_equal(v2.y(), 0.02);
    v1.set(0., 0., 0.);
    v2.set(0., 0., 0.);

    arc2->coord_at_dist(0.2, false, v1);
    arc2r->coord_at_dist(0.2, false, v2);
    check_float_equal(v1.x(), 1.662336879396);
    check_float_equal(v1.y(), 1.749206151998);
    check_float_equal(v2.x(), 1.662336879396);
    check_float_equal(v2.y(), 1.749206151998);

    //Testing coord_at_mid_dist
    arc2->coord_at_mid_dist(0.1, 0.3, v1);
    arc2->coord_at_param(0.2, v2);
    check_float_equal(v1.x(), v2.x());
    check_float_equal(v1.y(), v2.y());
    arc2->coord_at_param(0.2, v1);
    arc2->coord_at_param(0.6, v2);
    //arc2->coord_at_mid_dist(v1, v2, v3);
    //	arc2->coord_at_param(0.4, v1);
    //	check_float_equal(v3.x(), v1.x());
    //	check_float_equal(v3.y(), v1.y());

    //Testing closest_coord_on_curve and closest_param_on_curve.
    CubitVector v(0., 0., 0.);

    v1.set(1. + 2. * cos(0.25 * M_PI), 1. + 2. * sin(0.25 * M_PI), 0.);
    v2.set(2., 1., 0.);
    v3.set(1., 2., 0.);
    arc1->closest_coord_on_curve(v1, v);
    check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(arc1->closest_param_on_curve(v1), 0.5);
    arc1->closest_coord_on_curve(v1, v, 0., 0.4);
    arc1->coord_at_param(0.4, v1);
    check_float_equal(v.x(), v1.x());
    check_float_equal(v.y(), v1.y());
    arc1->closest_coord_on_curve(v1, v, 0.6, 1.);
    arc1->coord_at_param(0.6, v1);
    check_float_equal(v.x(), v1.x());
    check_float_equal(v.y(), v1.y());
    arc1->closest_coord_on_curve(v2, v);
    check_float_equal(v.x(), 1. + 0.5 * sqrt(3));
    check_float_equal(v.y(), 1.5);
    check_float_equal(arc1->closest_param_on_curve(v2), 0.);
    arc1->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 1.5);
    check_float_equal(v.y(), 1. + 0.5 * sqrt(3));
    check_float_equal(arc1->closest_param_on_curve(v3), 1.);

    v1.set(1. + 2. * cos(0.25 * 7. * M_PI), 1. + 2. * sin(0.25 * 7. * M_PI),
	   0.);
    v2.set(1. + 2. * cos(0.25 * M_PI), 1. + 2. * sin(0.25 * M_PI), 0.);
    v3.set(0., 1., 0.);
    arc2->closest_coord_on_curve(v1, v);
    check_float_equal(v.x(), 1. + cos(0.25 * 7. * M_PI));
    check_float_equal(v.y(), 1. + sin(0.25 * 7. * M_PI));
    check_float_equal(arc2->closest_param_on_curve(v1), 0.3);
    arc2->closest_coord_on_curve(v2, v);
    check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v.y(), 1. + cos(0.25 * M_PI));
    check_float_equal(arc2->closest_param_on_curve(v2), 0.9);
    arc2->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 1.);
    check_float_equal(v.y(), 0.);
    check_float_equal(arc2->closest_param_on_curve(v3), 0.);
    v3.set(1. + 2. * cos(0.25 * 3. * M_PI), 1. + 2. * sin(0.25 * 3. * M_PI),
	   0.);
    arc2->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 1.5);
    check_float_equal(v.y(), 1. + 0.5 * sqrt(3.));
    check_float_equal(arc2->closest_param_on_curve(v3), 1.);

    v1.set(1. + 2. * cos(0.25 * M_PI), 1. + 2. * sin(0.25 * M_PI), 0.);
    v2.set(3., 1., 0.);
    v3.set(1. + 2. * cos(0.25 * 7. * M_PI), 1. + 2. * sin(0.25 * 7. * M_PI),
	   0.);
    arc3->closest_coord_on_curve(v1, v);
    check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(arc3->closest_param_on_curve(v1), 0.0625);
    arc3->closest_coord_on_curve(v2, v);
    check_float_equal(v.x(), 1. + 0.5 * sqrt(3));
    check_float_equal(v.y(), 1.5);
    check_float_equal(arc3->closest_param_on_curve(v2), 0.);
    arc3->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 1.);
    check_float_equal(v.y(), 0.);
    check_float_equal(arc3->closest_param_on_curve(v3), 1.);

    v1.set(1. + 2. * cos(0.25 * 7. * M_PI), 1. + 2. * sin(0.25 * 7. * M_PI),
	   0.);
    v2.set(1. + 2. * cos(0.25 * M_PI), 1. + 2. * sin(0.25 * M_PI), 0.);
    v3.set(1. + 2. * cos(1.1 * M_PI), 1. + 2. * sin(1.1 * M_PI), 0.);
    arc4->closest_coord_on_curve(v1, v);
    check_float_equal(v.x(), 1. + cos(0.25 * 7. * M_PI));
    check_float_equal(v.y(), 1. + sin(0.25 * 7. * M_PI));
    check_float_equal(arc4->closest_param_on_curve(v1), 1. / 6.);
    arc4->closest_coord_on_curve(v2, v);
    check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v.y(), 1. + cos(0.25 * M_PI));
    check_float_equal(arc4->closest_param_on_curve(v2), 0.5);
    arc4->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 0.);
    check_float_equal(v.y(), 1.);
    check_float_equal(arc4->closest_param_on_curve(v3), 1.);
    v3.set(1. + 2. * cos(1.4 * M_PI), 1. + 2. * sin(1.4 * M_PI), 0.);
    arc4->closest_coord_on_curve(v3, v);
    check_float_equal(v.x(), 1.);
    check_float_equal(v.y(), 0.);
    check_float_equal(arc4->closest_param_on_curve(v3), 0.);

    //Testing rest of coord_on_curve (other parts tested before).
    v1.set(1.5, 2., 0.);
    v2.set(1. + cos(0.25 * 7. * M_PI), 1. + sin(0.25 * 7. * M_PI), 0.);
    v3.set(1. + cos(0.25 * 3. * M_PI), 1. + sin(0.25 * 3. * M_PI), 0.);
    BOOST_CHECK(arc1->coord_on_curve(v1) == false);
    BOOST_CHECK(arc1->coord_on_curve(v2) == false);
    BOOST_CHECK(arc2->coord_on_curve(v3) == false);
    BOOST_CHECK(arc3->coord_on_curve(v2) == false);
    v1.set(1. + cos(0.25 * 5. * M_PI), 1. + sin(0.25 * 5. * M_PI), 0.);
    BOOST_CHECK(arc4->coord_on_curve(v1) == false);

    delete arc2;
    delete arc3;
    delete arc1r;
    delete arc2r;
    delete arc3r;
    delete arc4r;

    //Testing first_deriv and second_deriv
    CubitVector D1, D2;
    v.set(1. + cos(0.25 * M_PI), 1. + sin(0.25 * M_PI), 0.);

    arc1->first_deriv(arc1->param_at_coord(v), D1);
    //			arc1->first_deriv(v, D2);
    check_float_equal(D1.x(), -0.37024024484);
    check_float_equal(D1.y(), 0.37024024484);
    check_float_equal(D1.z(), 0.);
    //			check_float_equal(D2.x(), -0.37024024484);
    //			check_float_equal(D2.y(), 0.37024024484);
    //			check_float_equal(D2.z(), 0.);
    arc1->second_deriv(arc1->param_at_coord(v), D1);
    //			arc1->second_deriv(v, D2);
    check_float_equal(D1.x(), -0.1938573388);
    check_float_equal(D1.y(), -0.1938573388);
    check_float_equal(D1.z(), 0.);
    //			check_float_equal(D2.x(), -0.1938573388);
    //			check_float_equal(D2.y(), -0.1938573388);
    //			check_float_equal(D2.z(), 0.);

    //Testing curvature.
    v1.set(1.5, 1. + 0.5 * sqrt(3), 0.);
    check_float_equal(arc1->curvature(0.5), 1.);
    //			check_float_equal(arc1->curvature(v1), 1.);

    //Testing normal and tangent
    arc1->unit_tangent(arc1->param_at_coord(v), D1);
    //	arc1->unit_tangent(v, D2);
    check_float_equal(D1.x(), -0.5 * sqrt(2));
    check_float_equal(D1.y(), 0.5 * sqrt(2));
    //	check_float_equal(D2.x(), -0.5 * sqrt(2));
    //	check_float_equal(D2.y(), 0.5 * sqrt(2));
    arc1->unit_normal(arc1->param_at_coord(v), D1);
    //	arc1->unit_normal(v, D2);
    check_float_equal(D1.x(), -0.5 * sqrt(2));
    check_float_equal(D1.y(), -0.5 * sqrt(2));
    //	check_float_equal(D2.x(), -0.5 * sqrt(2));
    //	check_float_equal(D2.y(), -0.5 * sqrt(2));

    //Testing arc_length.
    check_float_equal(arc1->arc_length(), 1. / 6. * M_PI);
    check_float_equal(arc4->arc_length(), 1.5 * M_PI);
    check_float_equal(arc1->arc_length(0.2, 0.6), 1. / 6. * 0.4 * M_PI);
    check_float_equal(arc4->arc_length(0.2, 0.6), 1.5 * 0.4 * M_PI);
    v1.set(1. + 0.5 * sqrt(3), 1.5, 0.);
    v2.set(1.5, 1. + 0.5 * sqrt(3), 0.);
    //	check_float_equal(arc1->arc_length(v1, v2), 1. / 6. * M_PI);
    //	check_float_equal(arc4->arc_length(v1, v2), 1. / 6. * M_PI);

    //	arc1->center_point(v);
    //	check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    //	check_float_equal(v.y(), 1. + sin(0.25 * M_PI));
    //	check_float_equal(v.z(), 0.);
    v.set(-10, -10, 10);
    arc1->center_point(arc1->min_param(), arc1->max_param(), v);
    check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    check_float_equal(v.y(), 1. + sin(0.25 * M_PI));
    check_float_equal(v.z(), 0.);
    //	v.set(-10, -10, 10);
    //	arc1->center_point(v1, v2, v);
    //	check_float_equal(v.x(), 1. + cos(0.25 * M_PI));
    //	check_float_equal(v.y(), 1. + sin(0.25 * M_PI));
    //	check_float_equal(v.z(), 0.);

    //Testing TVT and mid_TVT.

    //	check_float_equal(arc1->TVT(), 1. / 6. * M_PI);
    //	check_float_equal(arc4->TVT(), 1.5 * M_PI);
    check_float_equal(arc1->TVT(0.2, 0.6), 1. / 6. * 0.4 * M_PI);
    check_float_equal(arc4->TVT(0.2, 0.6), 1.5 * 0.4 * M_PI);
    //	check_float_equal(arc1->arc_length(v1, v2), 1. / 6. * M_PI);
    //	check_float_equal(arc4->arc_length(v1, v2), 1. / 6. * M_PI);
    check_float_equal(arc1->mid_TVT(0.1, 0.7), 0.4);
    check_float_equal(arc4->mid_TVT(0.1, 0.7), 0.4);
    arc1->coord_at_param(0.1, v3);
    arc1->coord_at_param(0.7, v);
    //	check_float_equal(arc1->mid_TVT(v3, v), 0.4);

    //Testing coord query and modification functions.
    check_float_equal(arc1->x(0), 1.);
    check_float_equal(arc1->y(0), 1.);
    arc1->get_coord(0, v);
    check_float_equal(v.x(), 1.);
    check_float_equal(v.y(), 1.);
    check_float_equal(arc1->get_coord(0).x(), 1.);
    check_float_equal(arc1->get_coord(0).y(), 1.);
    check_float_equal(arc1->get_coord(0).z(), 0.);

    arc1->set_x(0, 2.);
    arc1->set_y(0, 2.);
    check_float_equal(arc1->x(0), 2.);
    check_float_equal(arc1->y(0), 2.);
    check_float_equal(arc1->x(1), 2. + 0.5 * sqrt(3));
    check_float_equal(arc1->y(1), 2.5);
    check_float_equal(arc1->x(2), 2.5);
    check_float_equal(arc1->y(2), 2. + 0.5 * sqrt(3));

    arc1->set_x(1, 2.);
    arc1->set_y(1, 2.);
    check_float_equal(arc1->x(0), 2. - 0.5 * sqrt(3));
    check_float_equal(arc1->y(0), 1.5);
    check_float_equal(arc1->x(1), 2.);
    check_float_equal(arc1->y(1), 2.);
    check_float_equal(arc1->x(2), 2.5 - 0.5 * sqrt(3));
    check_float_equal(arc1->y(2), 1.5 + 0.5 * sqrt(3));

    arc1->set_x(2, 2.);
    arc1->set_y(2, 2.);
    check_float_equal(arc1->x(0), 1.5);
    check_float_equal(arc1->y(0), 2. - 0.5 * sqrt(3));
    check_float_equal(arc1->x(1), 1.5 + 0.5 * sqrt(3));
    check_float_equal(arc1->y(1), 2.5 - 0.5 * sqrt(3));
    check_float_equal(arc1->x(2), 2.);
    check_float_equal(arc1->y(2), 2.);

    arc1->set_coord(0, v);
    check_float_equal(arc1->x(0), 1.);
    check_float_equal(arc1->y(0), 1.);
    check_float_equal(arc1->x(1), 1. + 0.5 * sqrt(3));
    check_float_equal(arc1->y(1), 1.5);
    check_float_equal(arc1->x(2), 1.5);
    check_float_equal(arc1->y(2), 1. + 0.5 * sqrt(3));

    //Testing copy constructor and operator=
    GRArc* arc_new = new GRArc(*arc1);
    GRArc other_arc;
    other_arc = (*arc1);

    check_float_equal(arc_new->x(0), 1.);
    check_float_equal(arc_new->y(0), 1.);
    check_float_equal(arc_new->x(1), 1. + 0.5 * sqrt(3));
    check_float_equal(arc_new->y(1), 1.5);
    check_float_equal(arc_new->x(2), 1.5);
    check_float_equal(arc_new->y(2), 1. + 0.5 * sqrt(3));
    check_float_equal(arc_new->arc_length(), 1. / 6. * M_PI);
    check_float_equal(other_arc.x(0), 1.);
    check_float_equal(other_arc.y(0), 1.);
    check_float_equal(other_arc.x(1), 1. + 0.5 * sqrt(3));
    check_float_equal(other_arc.y(1), 1.5);
    check_float_equal(other_arc.x(2), 1.5);
    check_float_equal(other_arc.y(2), 1. + 0.5 * sqrt(3));
    check_float_equal(other_arc.arc_length(), 1. / 6. * M_PI);

    delete arc1;
    delete arc4;
    delete arc_new;

  }

  BOOST_AUTO_TEST_CASE(TestCubic)
  {

    double* x_coeff = new double[4];
    double* y_coeff = new double[4];

    x_coeff[0] = 2.;
    x_coeff[1] = -2.;
    x_coeff[2] = 6.;
    x_coeff[3] = -5.;
    y_coeff[0] = -1.;
    y_coeff[1] = 4.;
    y_coeff[2] = -3.;
    y_coeff[3] = 4.;

    GRCubic* cubic = new GRCubic(x_coeff, y_coeff);

    // check the coefficients
    const double * x_coeff_new = cubic->get_x_coeff();
    const double * y_coeff_new = cubic->get_y_coeff();
    check_float_equal(x_coeff_new[0], x_coeff[0]);
    check_float_equal(x_coeff_new[1], x_coeff[1]);
    check_float_equal(x_coeff_new[2], x_coeff[2]);
    check_float_equal(x_coeff_new[3], x_coeff[3]);
    check_float_equal(y_coeff_new[0], y_coeff[0]);
    check_float_equal(y_coeff_new[1], y_coeff[1]);
    check_float_equal(y_coeff_new[2], y_coeff[2]);
    check_float_equal(y_coeff_new[3], y_coeff[3]);
    BOOST_CHECK_EQUAL(cubic->topo_dim(), 1);
    BOOST_CHECK_EQUAL(cubic->geom_dim(), 2);
    BOOST_CHECK(cubic->get_curve_type() == GRCurveGeom::CUBIC);
    check_float_equal(cubic->min_param(), 0.);
    check_float_equal(cubic->max_param(), 1.);

    //Testing the bounding box
    CubitBox box = cubic->bounding_box();
    CubitVector mini = box.minimum();
    CubitVector maxi = box.maximum();
    check_float_equal(mini.x(), 1.);
    check_float_equal(mini.y(), -1.);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 2.);
    check_float_equal(maxi.y(), 4.);
    check_float_equal(maxi.z(), 0.);

    //Testing mapping between param and cartesian space.
    CubitVector v(1.875, 0.75, 0.);
    double param = cubic->param_at_coord(v);
    check_float_equal(param, 0.5);

    v.set(0., 0., 0.);
    cubic->coord_at_param(param, v);
    check_float_equal(v.x(), 1.875);
    check_float_equal(v.y(), 0.75);
    check_float_equal(v.z(), 0.);

    //Testing closest parameter from a point located near the curve.
    //Testing qCoordOnCurve from the returned coord.

    v.set(0., 5., 0.);
    check_float_equal(cubic->closest_param_on_curve(v), 1.);
    v.set(3., -2., 0.);
    check_float_equal(cubic->closest_param_on_curve(v), 0.);

    CubitVector normal;
    int num_tests = 20;

    for (int i = 0; i < num_tests; i++) {
      param = static_cast<double>(i) / static_cast<double>(num_tests);
      double move = -0.5 * pow(-1., i + 1);
      cubic->coord_at_param(param, v);
      cubic->unit_normal(param, normal);

      double x = v.x() + move * normal.x();
      double y = v.y() + move * normal.y();
      v.set(x, y, 0.);

      double tmp_param = cubic->closest_param_on_curve(v);
      //     printf("param = %e, closest = %e\n", param, tmp_param);

      check_float_equal(tmp_param, param);
      cubic->coord_at_param(param, v);

      BOOST_CHECK(cubic->coord_on_curve(v) == true);
    }

    v.set(1.3, 2.5, 0.);
    CubitVector new_v;
    cubic->closest_coord_on_curve(v, new_v);
    check_float_equal(new_v.x(), 1.595181970135);
    check_float_equal(new_v.y(), 2.598830521724);
    check_float_equal(new_v.z(), 0.);

    //Testing validity of first and second derivatives.

    param = 0.25;
    cubic->coord_at_param(param, v);
    CubitVector FD, SD;

    cubic->first_deriv(0.25, FD);
    cubic->second_deriv(0.25, SD);
    check_float_equal(FD.x(), 0.0625);
    check_float_equal(FD.y(), 3.25);
    check_float_equal(FD.z(), 0.);
    check_float_equal(SD.x(), 4.5);
    check_float_equal(SD.y(), 0.);
    check_float_equal(SD.z(), 0.);

    FD.set(0., 0., 0.);
    SD.set(0., 0., 0.);
    //	cubic->first_deriv(v, FD);
    //	cubic->second_deriv(v, SD);
    //	check_float_equal(FD.x(), 0.0625);
    //	check_float_equal(FD.y(), 3.25);
    //	check_float_equal(FD.z(), 0.);
    //	check_float_equal(SD.x(), 4.5);
    //	check_float_equal(SD.y(), 0.);
    //	check_float_equal(SD.z(), 0.);

    //Testing curvature
    check_float_equal(cubic->curvature(0.), 0.40249223594996);
    check_float_equal(cubic->curvature(0.34), 0.13850566115484);
    check_float_equal(cubic->curvature(1.), 0.06439875775199);

    cubic->coord_at_param(0.34, v);
    //	check_float_equal(cubic->curvature(v),0.13850566115484);

    //Testing tangent and normal vectors.
    CubitVector tangent;

    param = 0.75;
    cubic->coord_at_param(param, v);
    cubic->unit_tangent(param, tangent);
    cubic->unit_normal(param, normal);

    check_float_equal(tangent.x(), -0.22414769292143);
    check_float_equal(tangent.y(), 0.9745551866149);
    check_float_equal(tangent.z(), 0.);
    check_float_equal(normal.x(), -0.9745551866149);
    check_float_equal(normal.y(), -0.22414769292143);
    check_float_equal(normal.z(), 0.);

    tangent.set(0., 0., 0.);
    normal.set(0., 0., 0.);
    //	cubic->unit_tangent(v, tangent);
    //	cubic->unit_normal(v, normal);
    //
    //	check_float_equal(tangent.x(), -0.22414769292143);
    //	check_float_equal(tangent.y(), 0.9745551866149);
    //	check_float_equal(tangent.z(), 0.);
    //	check_float_equal(normal.x(), -0.9745551866149);
    //	check_float_equal(normal.y(), -0.22414769292143);
    //	check_float_equal(normal.z(), 0.);

    //Testing arc length calculation.
    cubic->coord_at_param(0.2, v);
    cubic->coord_at_param(0.6, new_v);

    check_float_equal(cubic->arc_length(), 5.19245698666);
    check_float_equal(cubic->arc_length(0.2, 0.6), 1.47603494216);
    //	check_float_equal(cubic->arc_length(v, new_v), 1.47603494216);

    //Testing center point
    //	cubic->center_point(v);
    //	check_float_equal(v.x(), 1.8467664119);
    //	check_float_equal(v.y(), 1.5577471815);
    //	check_float_equal(v.z(), 0.);
    //	cubic->coord_at_param(cubic->min_param(), v);
    //	cubic->coord_at_param(, new_v);
    cubic->center_point(cubic->min_param(), cubic->max_param(), v);
    check_float_equal(v.x(), 1.8467664119);
    check_float_equal(v.y(), 1.5577471815);
    check_float_equal(v.z(), 0.);

    //Testing total variation of the tangent angle.
    //  printf("%.16e %.16e\n", cubic->TVT(), cubic->TVT() - 1.15504911858 );
    //	check_float_equal(cubic->TVT(), 1.1550491335);
    check_float_equal(cubic->TVT(0.2, 0.6), 0.3310016578);
    //	cubic->coord_at_param(0.2, v);
    //	cubic->coord_at_param(0.6, new_v);
    //	check_float_equal(cubic->TVT(v, new_v), 0.3310016578);

    //Testing copy constructor and operator= (use arc length and TVT)

    GRCubic* new_cubic = new GRCubic(*cubic);
    check_float_equal(
	new_cubic->TVT(new_cubic->min_param(), new_cubic->max_param()),
	1.1550491335);
    check_float_equal(new_cubic->arc_length(), 5.19245698666);

    GRCubic other_cubic;
    other_cubic = *cubic;
    check_float_equal(
	other_cubic.TVT(other_cubic.min_param(), other_cubic.max_param()),
	1.1550491335);
    check_float_equal(other_cubic.arc_length(), 5.19245698666);
    // Testing set coeffs, default constructor

    GRCubic another_cubic;
    another_cubic.set_xcoeffs(x_coeff);
    another_cubic.set_ycoeffs(y_coeff);
    check_float_equal(
	another_cubic.TVT(another_cubic.min_param(), another_cubic.max_param()),
	1.1550491335);
    check_float_equal(another_cubic.arc_length(), 5.19245698666);

    delete cubic;
    delete new_cubic;
    delete [] x_coeff;
    delete [] y_coeff;
  }
  BOOST_AUTO_TEST_CASE(TestBezier)
  {

    CubitVector coord[4];
    coord[0].set(3., 2., 0.);
    coord[1].set(2., 4., 0.);
    coord[2].set(6., 3., 0.);
    coord[3].set(5., 4., 0.);

    GRBezier* bezier = new GRBezier(coord);

    BOOST_CHECK_EQUAL(bezier->topo_dim(), 1);
    BOOST_CHECK_EQUAL(bezier->geom_dim(), 2);
    BOOST_CHECK(bezier->get_curve_type() == GRCurveGeom::BEZIER);

    check_float_equal(bezier->min_param(), 0.);
    check_float_equal(bezier->max_param(), 1.);

    //Testing the bounding box
    CubitBox box = bezier->bounding_box();
    CubitVector mini = box.minimum();
    CubitVector maxi = box.maximum();
    check_float_equal(mini.x(), 2.8381049961);
    check_float_equal(mini.y(), 2.);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 5.1618950039);
    check_float_equal(maxi.y(), 4.);
    check_float_equal(maxi.z(), 0.);

    //Testing param_at_coord and param_at_arc_length
    CubitVector test_coord(4., 3.375, 0.);
    check_float_equal(bezier->param_at_coord(test_coord), 0.5);
    check_float_equal(bezier->param_at_arc_length(0.2, 1.1), 0.47307077932);

    //Testing coord_at_param
    test_coord.set(0., 0., 0.);
    bezier->coord_at_param(0.5, test_coord);
    check_float_equal(test_coord.x(), 4.);
    check_float_equal(test_coord.y(), 3.375);
    check_float_equal(test_coord.z(), 0.);

    //Testing coord_at_dist
    bezier->coord_at_param(bezier->min_param(), coord[0]);
    bezier->coord_at_param(0.2, coord[1]);
    bezier->coord_at_dist((coord[0] - coord[1]).length(), true, test_coord);
    check_float_equal(test_coord.x(), coord[1].x());
    check_float_equal(test_coord.y(), coord[1].y());
    check_float_equal(test_coord.z(), coord[1].z());

    bezier->coord_at_param(bezier->max_param(), coord[0]);
    bezier->coord_at_param(0.8, coord[1]);
    bezier->coord_at_dist((coord[0] - coord[1]).length(), false, test_coord);
    check_float_equal(test_coord.x(), coord[1].x());
    check_float_equal(test_coord.y(), coord[1].y());
    check_float_equal(test_coord.z(), coord[1].z());

    //Testing coord_at_mid_dist.
    //	bezier->coord_at_param(bezier->min_param(), coord[0]);
    //	bezier->coord_at_param(bezier->max_param(), coord[1]);
    //	//bezier->coord_at_mid_dist(coord[0], coord[1], coord[2]);
    //	check_float_equal((coord[0]-coord[2]).length(), (coord[1]-coord[2]).length());

    //Testing function to get closest point on curve.

    int n_tests = 20;
    for (int i = 0; i < n_tests; i++) {
      double param = static_cast<double>(i) / static_cast<double>(n_tests);
      double move = -0.1 * pow(-1., i + 1);
      CubitVector normal;
      bezier->coord_at_param(param, coord[0]);
      bezier->unit_normal(param, normal);

      coord[0].x(coord[0].x() + move * normal.x());
      coord[0].y(coord[0].y() + move * normal.y());

      check_float_equal(bezier->closest_param_on_curve(coord[0]), param);

      bezier->coord_at_param(param, coord[0]);
      BOOST_CHECK(bezier->coord_on_curve(coord[0]) == true);
    }

    coord[0].set(4., 3., 0.);
    bezier->closest_coord_on_curve(coord[0], coord[1]);
    check_float_equal(coord[1].x(), 3.93553046773);
    check_float_equal(coord[1].y(), 3.36392729063);
    check_float_equal(coord[1].z(), 0.);

    coord[0].set(3.5, 3.5, 0.);
    check_float_equal(bezier->closest_param_on_curve(coord[0], 0.8, 0.9), 0.8);
    check_float_equal(bezier->closest_param_on_curve(coord[0], 0.1, 0.2), 0.2);

    //Testing first and second derivatives.
    bezier->coord_at_param(0.25, coord[0]);

    bezier->first_deriv(0.25, coord[1]);
    bezier->second_deriv(0.25, coord[2]);
    check_float_equal(coord[1].x(), 2.625);
    check_float_equal(coord[1].y(), 2.4375);
    check_float_equal(coord[1].z(), 0.);
    check_float_equal(coord[2].x(), 15.);
    check_float_equal(coord[2].y(), -10.5);
    check_float_equal(coord[2].z(), 0.);

    coord[1].set(0., 0., 0.);
    coord[2].set(0., 0., 0.);
    //	bezier->first_deriv(coord[0], coord[1]);
    //	bezier->second_deriv(coord[0], coord[2]);
    //	check_float_equal(coord[1].x(), 2.625);
    //	check_float_equal(coord[1].y(), 2.4375);
    //	check_float_equal(coord[1].z(), 0.);
    //	check_float_equal(coord[2].x(), 15.);
    //	check_float_equal(coord[2].y(), -10.5);
    //	check_float_equal(coord[2].z(), 0.);

    //Testing tangent and normal vectors.
    bezier->coord_at_param(0.75, coord[0]);
    bezier->unit_tangent(0.75, coord[1]);
    bezier->unit_normal(0.75, coord[2]);
    check_float_equal(coord[1].x(), 0.94174191159484);
    check_float_equal(coord[1].y(), 0.33633639699816);
    check_float_equal(coord[1].z(), 0.);
    check_float_equal(coord[2].x(), -0.33633639699816);
    check_float_equal(coord[2].y(), 0.94174191159484);
    check_float_equal(coord[2].z(), 0.);

    coord[1].set(0., 0., 0.);
    coord[2].set(0., 0., 0.);
    //	bezier->unit_tangent(coord[0], coord[1]);
    //	bezier->unit_normal(coord[0], coord[2]);
    //	check_float_equal(coord[1].x(), 0.94174191159484);
    //	check_float_equal(coord[1].y(), 0.33633639699816);
    //	check_float_equal(coord[1].z(), 0.);
    //	check_float_equal(coord[2].x(), -0.33633639699816);
    //	check_float_equal(coord[2].y(), 0.94174191159484);
    //	check_float_equal(coord[2].z(), 0.);

    //Testing curvature
    check_float_equal(bezier->curvature(0.75), 1.19478051615393);
    bezier->coord_at_param(0.75, coord[0]);
    //	check_float_equal(bezier->curvature(coord[0]), 1.19478051615393);

    //Testing arc length
    check_float_equal(bezier->arc_length(), 3.7363838928543686);
    check_float_equal(bezier->arc_length(0.2, 0.6), 1.66765099753);
    check_float_equal(bezier->arc_length(0.2, 0.6), 1.66765099753);

    //Testing total variation tangent
    //	check_float_equal(bezier->TVT(), 4.1111687558);
    check_float_equal(bezier->TVT(0.2, 0.6), 0.89280419598);
    check_float_equal(bezier->arc_length(0.2, 0.6), 1.66765099753);

    //Testing center_point;
    bezier->center_point(bezier->min_param(), bezier->max_param(), coord[0]);
    double center_param = bezier->param_at_coord(coord[0]);
    double l1 = bezier->arc_length(bezier->min_param(), center_param);
    double l2 = bezier->arc_length(center_param, bezier->max_param());
    BOOST_CHECK(fabs(l1 - l2) < 1.e-6);

    bezier->center_point(0.2, 0.6, coord[0]);
    center_param = bezier->param_at_coord(coord[0]);
    l1 = bezier->arc_length(0.2, center_param);
    l2 = bezier->arc_length(center_param, 0.6);
    BOOST_CHECK(fabs(l1 - l2) < 1.e-6);

    //Testing mid_TVT
    // double param = bezier->mid_TVT(bezier->min_param(), bezier->max_param());
    // These two are ridiculously close to testing as equal.  Their
    // difference is something like 2% too large (and the tolerance on
    // relative difference is set at 5e-9).
    //   check_float_equal(bezier->TVT(bezier->min_param(), param),
    // 		   bezier->TVT(param, bezier->max_param()));

    //	bezier->coord_at_param(bezier->min_param(), coord[1]);
    //	bezier->coord_at_param(bezier->max_param(), coord[2]);
    //	param = bezier->mid_TVT(coord[1], coord[2]);
    //	check_float_equal(bezier->TVT(bezier->min_param(), param),
    //			bezier->TVT(param, bezier->max_param()));

    //Testing coord queries
    check_float_equal(bezier->x(0), 3.);
    check_float_equal(bezier->x(1), 2.);
    check_float_equal(bezier->x(2), 6.);
    check_float_equal(bezier->x(3), 5.);
    check_float_equal(bezier->y(0), 2.);
    check_float_equal(bezier->y(1), 4.);
    check_float_equal(bezier->y(2), 3.);
    check_float_equal(bezier->y(3), 4.);

    bezier->get_coord(1, coord[0]);
    check_float_equal(coord[0].x(), 2.);
    check_float_equal(coord[0].y(), 4.);
    check_float_equal(coord[0].z(), 0.);

    check_float_equal(bezier->get_coord(2).x(), 6.);
    check_float_equal(bezier->get_coord(2).y(), 3.);
    check_float_equal(bezier->get_coord(2).z(), 0.);

    //Testing copy constructor and operator=.
    GRBezier* new_bezier = new GRBezier(*bezier);
    check_float_equal(new_bezier->arc_length(), 3.7363838928543686);
    check_float_equal(
	new_bezier->TVT(new_bezier->min_param(), new_bezier->max_param()),
	4.1111687558);

    GRBezier other_bezier;
    other_bezier = *bezier;
    check_float_equal(other_bezier.arc_length(), 3.7363838928543686);
    check_float_equal(
	other_bezier.TVT(other_bezier.min_param(), other_bezier.max_param()),
	4.1111687558);

    //Testing functions to set coords.
    bezier->set_x(0, 4.);
    bezier->set_y(0, 1.);
    check_float_equal(bezier->arc_length(), 4.2174104931702931);
    check_float_equal(bezier->TVT(bezier->min_param(), bezier->max_param()),
		      3.949883536669);

    coord[0].set(3., 2., 0.);
    bezier->set_coord(0, coord[0]);
    check_float_equal(bezier->arc_length(), 3.7363838928543686);
    check_float_equal(bezier->TVT(bezier->min_param(), bezier->max_param()),
		      4.1111687558);

    delete bezier;
    delete new_bezier;

  }

  BOOST_AUTO_TEST_CASE(TestInterp)
  {
    // interp approximates a sine with an increasing number of points. Using a random
    // sample of points, compare the reduction in L1, L2 and LInf norms. With piecewise
    // cubic polynomials, doubling the number of interpolating points should decrease the
    // norm by around 8.

    CubitVector coord[320];

    int num_points = 5;
    std::vector<double> L1, L2, LInf;

    while (num_points <= 320) {

      coord[0].set(0., 0., 0.);
      coord[num_points - 1].set(2 * M_PI, 0., 0.);

      int i;
      for (i = 1; i < num_points - 1; i++) {
	double theta = 2 * M_PI * static_cast<double>(i)
	    / static_cast<double>(num_points - 1);
	coord[i].set(theta, sin(theta), 0.);
      }

      GRInterp* interp = new GRInterp(coord, num_points);

      double dL1 = 0.;
      double dL2 = 0.;
      double dLInf = 0.;

      //Using 1000 random sample points.
      for (i = 0; i < 1000; i++) {
	double param = drand48() * (num_points - 1);
	CubitVector location;
	interp->coord_at_param(param, location);
	double error = location.y() - sin(location.x());
	dL1 += fabs(error);
	dL2 += error * error;
	dLInf = std::max(dLInf, fabs(error));
      }

      L1.push_back(dL1 / 1000.);
      L2.push_back(sqrt(dL2) / 1000.);
      LInf.push_back(dLInf);

      num_points *= 2;
      delete interp;

    }

    //Testing the reduction in the norms: This test validates
    //the function computing the cubic polynomials.
    for (unsigned int i = 1; i < L1.size(); i++) {
      BOOST_CHECK(iFuzzyComp((L1[i - 1] / L1[i]), 6.) >= 0);
      BOOST_CHECK(iFuzzyComp((L2[i - 1] / L2[i]), 6.) >= 0);
      BOOST_CHECK(iFuzzyComp((LInf[i - 1] / LInf[i]), 6.) >= 0);
    }

    //Testing basic query functions. All tests will be conducted
    //with both a self-connected and a non self-connected spline.
    //Testing for a closed spline (periodic parameterization) will
    //also have to be conducted once the constructor is correctly
    //implemented

    num_points = 5;
    coord[0].set(0., 0., 0.);
    coord[1].set(0., 2., 0.);
    coord[2].set(2., 2., 0.);
    coord[3].set(2., 0., 0.);
    coord[4].set(4., 0., 0.);

    GRInterp* open = new GRInterp(coord, num_points);

    coord[4].set(0., 0., 0.);
    GRInterp* closed = new GRInterp(coord, num_points);

    BOOST_CHECK_EQUAL(open->topo_dim(), 1);
    BOOST_CHECK_EQUAL(closed->topo_dim(), 1);
    BOOST_CHECK_EQUAL(open->geom_dim(), 2);
    BOOST_CHECK_EQUAL(closed->geom_dim(), 2);
    BOOST_CHECK(open->get_curve_type() == GRCurveGeom::INTERP);
    BOOST_CHECK(closed->get_curve_type() == GRCurveGeom::INTERP);

    //Testing point queries. Duplicating a lot of the tests,
    //just to be on the safe side.

    //Testing a bunch of query functions for at least one param on each cubic curve.
    open->coord_at_param(0., coord[0]);
    open->unit_normal(0., coord[1]);
    open->unit_tangent(0., coord[2]);
    open->first_deriv(0., coord[3]);
    open->second_deriv(0., coord[4]);
    check_float_equal(coord[0].x(), 0.);
    check_float_equal(coord[0].y(), 0.);
    check_float_equal(coord[1].x(), -0.66436383882991978);
    check_float_equal(coord[1].y(), -0.74740931868365978);
    check_float_equal(coord[2].x(), -0.74740931868365978);
    check_float_equal(coord[2].y(), 0.66436383882991978);
    check_float_equal(coord[3].x(), -3.);
    check_float_equal(coord[3].y(), 2.6666666666666665);
    check_float_equal(coord[4].x(), 8.);
    check_float_equal(coord[4].y(), -0.99999999999999956);
    check_float_equal(open->curvature(0.), 0.28350008639725033);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 0.);

    open->coord_at_param(0.3, coord[0]);
    open->unit_normal(0.3, coord[1]);
    open->unit_tangent(0.3, coord[2]);
    open->first_deriv(0.3, coord[3]);
    open->second_deriv(0.3, coord[4]);
    check_float_equal(coord[0].x(), -0.56699999999999995);
    check_float_equal(coord[0].y(), 0.75050000000000006);
    check_float_equal(coord[1].x(), -0.93641202877174212);
    check_float_equal(coord[1].y(), -0.35090242571345975);
    check_float_equal(coord[2].x(), -0.35090242571345975);
    check_float_equal(coord[2].y(), 0.93641202877174212);
    check_float_equal(coord[3].x(), -0.87000000000000011);
    check_float_equal(coord[3].y(), 2.3216666666666668);
    check_float_equal(coord[4].x(), 6.2000000000000002);
    check_float_equal(coord[4].y(), -1.2999999999999996);
    check_float_equal(open->curvature(0.3), 0.87027005019339909);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 0.3);

    open->coord_at_param(1.4, coord[0]);
    open->unit_normal(1.4, coord[1]);
    open->unit_tangent(1.4, coord[2]);
    open->first_deriv(1.4, coord[3]);
    open->second_deriv(1.4, coord[4]);
    check_float_equal(coord[0].x(), 0.8959999999999998);
    check_float_equal(coord[0].y(), 2.2959999999999998);
    check_float_equal(coord[1].x(), -0.12263061059432923);
    check_float_equal(coord[1].y(), 0.99245238341457076);
    check_float_equal(coord[2].x(), 0.99245238341457076);
    check_float_equal(coord[2].y(), 0.12263061059432923);
    check_float_equal(coord[3].x(), 2.3200000000000003);
    check_float_equal(coord[3].y(), 0.2866666666666669);
    check_float_equal(coord[4].x(), -0.39999999999999947);
    check_float_equal(coord[4].y(), -2.3999999999999999);
    check_float_equal(open->curvature(1.4), 0.42690095572322628);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 1.4);

    open->coord_at_param(2.6, coord[0]);
    open->unit_normal(2.6, coord[1]);
    open->unit_tangent(2.6, coord[2]);
    open->first_deriv(2.6, coord[3]);
    open->second_deriv(2.6, coord[4]);
    check_float_equal(coord[0].x(), 2.0960000000000001);
    check_float_equal(coord[0].y(), 0.83999999999999964);
    check_float_equal(coord[1].x(), 0.9898903092956709);
    check_float_equal(coord[1].y(), -0.14183502939161854);
    check_float_equal(coord[2].x(), -0.14183502939161854);
    check_float_equal(coord[2].y(), -0.9898903092956709);
    check_float_equal(coord[3].x(), -0.32000000000000006);
    check_float_equal(coord[3].y(), -2.2333333333333334);
    check_float_equal(coord[4].x(), -0.39999999999999947);
    check_float_equal(coord[4].y(), 4.4408920985006262e-16);
    check_float_equal(open->curvature(2.6), 0.077788270077107383);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 2.6);

    open->coord_at_param(3.7, coord[0]);
    open->unit_normal(3.7, coord[1]);
    open->unit_tangent(3.7, coord[2]);
    open->first_deriv(3.7, coord[3]);
    open->second_deriv(3.7, coord[4]);
    check_float_equal(coord[0].x(), 2.8330000000000006);
    check_float_equal(coord[0].y(), -0.50749999999999984);
    check_float_equal(coord[1].x(), -0.26591105932845716);
    check_float_equal(coord[1].y(), 0.96399756666021619);
    check_float_equal(coord[2].x(), 0.96399756666021619);
    check_float_equal(coord[2].y(), 0.26591105932845716);
    check_float_equal(coord[3].x(), 2.870000000000001);
    check_float_equal(coord[3].y(), 0.79166666666666785);
    check_float_equal(coord[4].x(), 6.2000000000000011);
    check_float_equal(coord[4].y(), 5.5000000000000009);
    check_float_equal(open->curvature(3.7), 0.41217148391449304);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 3.7);

    open->coord_at_param(4., coord[0]);
    open->unit_normal(4., coord[1]);
    open->unit_tangent(4., coord[2]);
    open->first_deriv(4., coord[3]);
    open->second_deriv(4., coord[4]);
    check_float_equal(coord[0].x(), 4.);
    check_float_equal(coord[0].y(), 2.2204460492503131e-16);
    check_float_equal(coord[1].x(), -0.4705882352941177);
    check_float_equal(coord[1].y(), 0.88235294117647056);
    check_float_equal(coord[2].x(), 0.88235294117647056);
    check_float_equal(coord[2].y(), 0.4705882352941177);
    check_float_equal(coord[3].x(), 5.);
    check_float_equal(coord[3].y(), 2.666666666666667);
    check_float_equal(coord[4].x(), 8.);
    check_float_equal(coord[4].y(), 7.);
    check_float_equal(open->curvature(4.), 0.075106859352737615);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(open->closest_param_on_curve(coord[0]), 4.);

    closed->coord_at_param(0., coord[0]);
    closed->unit_normal(0., coord[1]);
    closed->unit_tangent(0., coord[2]);
    closed->first_deriv(0., coord[3]);
    closed->second_deriv(0., coord[4]);
    check_float_equal(coord[0].x(), 0.);
    check_float_equal(coord[0].y(), 0.);
    check_float_equal(coord[1].x(), -0.70710678118654746);
    check_float_equal(coord[1].y(), -0.70710678118654757);
    check_float_equal(coord[2].x(), -0.70710678118654757);
    check_float_equal(coord[2].y(), 0.70710678118654746);
    check_float_equal(coord[3].x(), -2.666666666666667);
    check_float_equal(coord[3].y(), 2.6666666666666665);
    check_float_equal(coord[4].x(), 7.0000000000000009);
    check_float_equal(coord[4].y(), -0.99999999999999956);
    check_float_equal(closed->curvature(0.), 0.29831067331307481);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(closed->closest_param_on_curve(coord[0]), 0.);

    closed->coord_at_param(0.3, coord[0]);
    closed->unit_normal(0.3, coord[1]);
    closed->unit_tangent(0.3, coord[2]);
    closed->first_deriv(0.3, coord[3]);
    closed->second_deriv(0.3, coord[4]);
    check_float_equal(coord[0].x(), -0.50750000000000006);
    check_float_equal(coord[0].y(), 0.75050000000000006);
    check_float_equal(coord[1].x(), -0.94648660423646402);
    check_float_equal(coord[1].y(), -0.32274309907560694);
    check_float_equal(coord[2].x(), -0.32274309907560694);
    check_float_equal(coord[2].y(), 0.94648660423646402);
    check_float_equal(coord[3].x(), -0.79166666666666685);
    check_float_equal(coord[3].y(), 2.3216666666666668);
    check_float_equal(coord[4].x(), 5.5000000000000009);
    check_float_equal(coord[4].y(), -1.2999999999999996);
    check_float_equal(closed->curvature(0.3), 0.79544821923018993);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(closed->closest_param_on_curve(coord[0]), 0.3);

    closed->coord_at_param(2.6, coord[0]);
    closed->unit_normal(2.6, coord[1]);
    closed->unit_tangent(2.6, coord[2]);
    closed->first_deriv(2.6, coord[3]);
    closed->second_deriv(2.6, coord[4]);
    check_float_equal(coord[0].x(), 2.2959999999999998);
    check_float_equal(coord[0].y(), 0.83999999999999964);
    check_float_equal(coord[1].x(), 0.99186250218518413);
    check_float_equal(coord[1].y(), -0.12731369431033715);
    check_float_equal(coord[2].x(), -0.12731369431033715);
    check_float_equal(coord[2].y(), -0.99186250218518413);
    check_float_equal(coord[3].x(), -0.28666666666666685);
    check_float_equal(coord[3].y(), -2.2333333333333334);
    check_float_equal(coord[4].x(), -2.3999999999999995);
    check_float_equal(coord[4].y(), 4.4408920985006262e-16);
    check_float_equal(closed->curvature(2.6), 0.46952482702456244);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(closed->closest_param_on_curve(coord[0]), 2.6);

    closed->coord_at_param(3.7, coord[0]);
    closed->unit_normal(3.7, coord[1]);
    closed->unit_tangent(3.7, coord[2]);
    closed->first_deriv(3.7, coord[3]);
    closed->second_deriv(3.7, coord[4]);
    check_float_equal(coord[0].x(), 0.75049999999999961);
    check_float_equal(coord[0].y(), -0.50749999999999984);
    check_float_equal(coord[1].x(), -0.32274309907560744);
    check_float_equal(coord[1].y(), -0.94648660423646402);
    check_float_equal(coord[2].x(), -0.94648660423646402);
    check_float_equal(coord[2].y(), 0.32274309907560744);
    check_float_equal(coord[3].x(), -2.3216666666666663);
    check_float_equal(coord[3].y(), 0.79166666666666785);
    check_float_equal(coord[4].x(), -1.2999999999999989);
    check_float_equal(coord[4].y(), 5.5000000000000009);
    check_float_equal(closed->curvature(3.7), 0.79544821923019038);
    coord[0].x(coord[0].x() + 0.1 * coord[1].x());
    coord[0].y(coord[0].y() + 0.1 * coord[1].y());
    check_float_equal(closed->closest_param_on_curve(coord[0]), 3.7);

    closed->coord_at_param(4., coord[0]);
    closed->unit_normal(4., coord[1]);
    closed->unit_tangent(4., coord[2]);
    closed->first_deriv(4., coord[3]);
    closed->second_deriv(4., coord[4]);
    check_float_equal(coord[0].x(), 4.4408920985006262e-16);
    check_float_equal(coord[0].y(), 2.2204460492503131e-16);
    check_float_equal(coord[1].x(), -0.70710678118654757);
    check_float_equal(coord[1].y(), -0.70710678118654735);
    check_float_equal(coord[2].x(), -0.70710678118654735);
    check_float_equal(coord[2].y(), 0.70710678118654757);
    check_float_equal(coord[3].x(), -2.6666666666666661);
    check_float_equal(coord[3].y(), 2.666666666666667);
    check_float_equal(coord[4].x(), -0.99999999999999911);
    check_float_equal(coord[4].y(), 7.);
    check_float_equal(closed->curvature(4.), 0.29831067331307476);

    //Same tests, this time with coords.
    coord[0].set(-0.195, 0.6975, 0.);
    BOOST_CHECK(open->coord_on_curve(coord[0]) == false);
    //	open->unit_normal(coord[0], coord[1]);
    //	open->unit_tangent(coord[0], coord[2]);
    //	//	open->first_deriv(coord[0], coord[3]);
    //	//	open->second_deriv(coord[0], coord[4]);
    //	check_float_equal(open->param_at_coord(coord[0]), 0.3);
    //	check_float_equal(coord[1].y(), -0.22474707977348);
    //	check_float_equal(coord[2].x(), -0.22474707977348);
    //	check_float_equal(coord[2].y(), 0.97441713353845);
    //	check_float_equal(coord[3].x(), -0.52142857142857);
    //	check_float_equal(coord[3].y(), 2.26071428571429);
    //	check_float_equal(coord[4].x(), 1.28571428571429);
    //	check_float_equal(coord[4].y(), -0.64285714285714);
    //	check_float_equal(open->curvature(coord[0]), 0.20590749929525);
    //	coord[0].x(coord[0].x() + coord[1].x() * 0.1);
    //	coord[0].y(coord[0].y() + coord[1].y() * 0.1);
    open->closest_coord_on_curve(coord[0], coord[1]);
    check_float_equal(open->param_at_coord(coord[1]), 0.21049655174688153);

    //Testing arc length.
    check_float_equal(open->arc_length(), 9.13274277552539);
    check_float_equal(closed->arc_length(), 9.0331723647651501);
    check_float_equal(open->arc_length(1.2, 3.4), 4.3052164228476411);
    check_float_equal(closed->arc_length(1.2, 3.4), 4.648416768398719);
    check_float_equal(open->arc_length(2.1, 2.6), 1.0230166802354641);
    check_float_equal(closed->arc_length(2.1, 2.6), 1.0484190883876627);

    //	open->coord_at_param(2.1, coord[0]);
    //	open->coord_at_param(2.6, coord[1]);
    //	check_float_equal(open->arc_length(coord[0], coord[1]), 1.0997069818);
    //	closed->coord_at_param(2.1, coord[0]);
    //	closed->coord_at_param(2.6, coord[1]);
    //	check_float_equal(closed->arc_length(coord[0], coord[1]), 1.12917271033);

    //Testing TVT
    check_float_equal(open->TVT(open->min_param(), open->max_param()),
		      6.3444900764484364);
    check_float_equal(closed->TVT(closed->min_param(), closed->max_param()),
		      6.2831853071795862);

    check_float_equal(open->TVT(1.2, 3.4), 3.2965883648841574);
    check_float_equal(closed->TVT(1.2, 3.4), 3.1543982776820241);
    check_float_equal(open->TVT(2.1, 2.6), 0.51565426556834204);
    check_float_equal(closed->TVT(2.1, 2.6), 0.70093308816358535);

    //	open->coord_at_param(2.1, coord[0]);
    //	open->coord_at_param(2.6, coord[1]);
    //	check_float_equal(open->TVT(coord[0], coord[1]), 0.464601524786);
    //	closed->coord_at_param(2.1, coord[0]);
    //	closed->coord_at_param(2.6, coord[1]);
    //	check_float_equal(closed->TVT(coord[0], coord[1]), 0.716270078694);

    //Testing the bounding box
    CubitBox box = open->bounding_box();
    CubitVector mini = box.minimum();
    CubitVector maxi = box.maximum();
    check_float_equal(mini.x(), -0.63113030944089876);
    check_float_equal(mini.y(), -0.56724287537299878);
    check_float_equal(mini.z(), 0.);
    check_float_equal(maxi.x(), 4.);
    check_float_equal(maxi.y(), 2.3128464544497871);
    check_float_equal(maxi.z(), 0.);

    //Testing param_at_arc_length
    double param = open->param_at_arc_length(0.2, 4.);
    check_float_equal(open->arc_length(0.2, param), 4.);

    //Testing coord_at_dist.
    open->coord_at_dist(0.5, true, coord[0]);
    open->coord_at_dist(0.5, false, coord[1]);
    open->coord_at_param(open->min_param(), coord[2]);
    open->coord_at_param(open->max_param(), coord[3]);
    double dist1 = (coord[0] - coord[2]).length();
    double dist2 = (coord[1] - coord[3]).length();
    check_float_equal(dist1, 0.5);
    check_float_equal(dist2, 0.5);

    closed->coord_at_dist(0.5, true, coord[0]);
    closed->coord_at_dist(0.5, false, coord[1]);
    closed->coord_at_param(closed->min_param(), coord[2]);
    closed->coord_at_param(closed->max_param(), coord[3]);
    dist1 = (coord[0] - coord[2]).length();
    dist2 = (coord[1] - coord[3]).length();
    check_float_equal(dist1, 0.5);
    check_float_equal(dist2, 0.5);

    //Testing coord_at_mid_dist
    open->coord_at_mid_dist(open->min_param(), open->max_param(), coord[0]);
    open->coord_at_param(open->min_param(), coord[1]);
    open->coord_at_param(open->max_param(), coord[2]);
    dist1 = (coord[0] - coord[1]).length();
    dist2 = (coord[0] - coord[2]).length();
    check_float_equal(dist1, dist2);

    //open->coord_at_mid_dist(coord[1], coord[2], coord[0]);
    dist1 = (coord[0] - coord[1]).length();
    dist2 = (coord[0] - coord[2]).length();
    check_float_equal(dist1, dist2);

    closed->coord_at_mid_dist(closed->min_param(), closed->max_param(),
			      coord[0]);
    closed->coord_at_param(closed->min_param(), coord[1]);
    closed->coord_at_param(closed->max_param(), coord[2]);
    dist1 = (coord[0] - coord[1]).length();
    dist2 = (coord[0] - coord[2]).length();
    check_float_equal(dist1, dist2);

    //Testing center_point
    open->center_point(open->min_param(), open->max_param(), coord[0]);
    param = open->param_at_coord(coord[0]);
    dist1 = open->arc_length(open->min_param(), param);
    dist2 = open->arc_length(param, open->max_param());
    check_float_equal(dist1, dist2);

    closed->center_point(closed->min_param(), closed->max_param(), coord[0]);
    param = closed->param_at_coord(coord[0]);
    dist1 = closed->arc_length(closed->min_param(), param);
    dist2 = closed->arc_length(param, closed->max_param());
    check_float_equal(dist1, dist2);

    //Testing mid_TVT
    param = open->mid_TVT(open->min_param(), open->max_param());
    double tvt1 = open->TVT(open->min_param(), param);
    double tvt2 = open->TVT(param, open->max_param());
    check_float_equal(tvt1, tvt2);

    param = closed->mid_TVT(closed->min_param(), closed->max_param());
    tvt1 = closed->TVT(closed->min_param(), param);
    tvt2 = closed->TVT(param, closed->max_param());
    BOOST_CHECK(fabs(tvt1 - tvt2) < 1.e-7);

    //Testing point coordinates query / setting functions

    check_float_equal(open->x(0), 0.);
    check_float_equal(open->x(1), 0.);
    check_float_equal(open->x(2), 2.);
    check_float_equal(open->x(3), 2.);
    check_float_equal(open->x(4), 4.);
    check_float_equal(open->x(0), 0.);
    check_float_equal(open->y(1), 2.);
    check_float_equal(open->y(2), 2.);
    check_float_equal(open->y(3), 0.);
    check_float_equal(open->y(4), 0.);

    open->get_coord(2, coord[0]);
    check_float_equal(coord[0].x(), 2.);
    check_float_equal(coord[0].y(), 2.);

    check_float_equal(open->get_coord(3).x(), 2.);
    check_float_equal(open->get_coord(3).y(), 0.);

    closed->set_x(0, 0.5);
    closed->set_y(0, 0.5);
    check_float_equal(closed->x(0), 0.5);
    check_float_equal(closed->y(0), 0.5);
    closed->set_x(0, 0.);
    closed->set_y(0, 0.);
    check_float_equal(closed->arc_length(), 8.52015078201);
    check_float_equal(closed->TVT(closed->min_param(), closed->max_param()),
		      5.24489309358);

    closed->set_x(0, -0.5);
    closed->set_y(0, -0.5);
    coord[0].set(0., 0., 0.);
    closed->set_coord(0, coord[0]);
    check_float_equal(closed->x(0), 0.);
    check_float_equal(closed->y(0), 0.);
    check_float_equal(closed->arc_length(), 8.52015078201);
    check_float_equal(closed->TVT(closed->min_param(), closed->max_param()),
		      5.24489309358);

    //Testing copy constructor and operator=

    GRInterp* new_interp = new GRInterp(*open);
    GRInterp other_interp;
    other_interp = *open;

    check_float_equal(new_interp->arc_length(), 9.13274277552539);
    check_float_equal(
	new_interp->TVT(new_interp->min_param(), new_interp->max_param()),
	6.3444900764484364);
    check_float_equal(other_interp.arc_length(), 9.13274277552539);
    check_float_equal(
	other_interp.TVT(other_interp.min_param(), other_interp.max_param()),
	6.3444900764484364);

    delete new_interp;
    delete open;
    delete closed;

  }
// Also does any of the point testing here.
  BOOST_AUTO_TEST_CASE(TestCurve)
  {

    CubitVector V1(0., 1., 0.);
    CubitVector V2(0., 2., 0.);
    GRPoint * p1 = new GRPoint(V1);
    GRPoint * p2 = new GRPoint(V2);

    GRLine* line = new GRLine(V1, V2);
    BOOST_CHECK(line->isValid());
    GRCurve* curve = new GRCurve(line, p1, p2);

    BOOST_CHECK(curve->get_curve_geom() == line);

    TBPoint * p3 = curve->get_start_point();
    TBPoint * p4 = curve->get_end_point();
    BOOST_CHECK(p3 == p1);
    BOOST_CHECK(p4 == p2);
    DLIList<CubitString*> string_list;
    DLIList<double> double_list;
    DLIList<int> int_list;

    CubitString cond_string("CURVE_CONDITION");
    CubitString reg_left("REG_LEFT");
    CubitString bdry_left("BDRY_LEFT");
    CubitString reg_right("REG_RIGHT");
    CubitString bdry_right("BDRY_RIGHT");

    string_list.append(&cond_string);
    double_list.append(0.);
    int_list.append(0);

    // both regions set to 1, no need to try anything else funny.
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(1);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(2);
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(2);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(1);

    CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
    curve->append_simple_attribute_virt(&csa);

    BOOST_CHECK(curve->qValid());

    // checks basic functionality
    BOOST_CHECK_EQUAL(2, curve->boundary_left());
    BOOST_CHECK_EQUAL(1, curve->boundary_right());
    BOOST_CHECK_EQUAL(1, curve->region_left());
    BOOST_CHECK_EQUAL(2, curve->region_right());
    CubitVector V3(0., 1.5, 0.);
    CubitVector CP = curve->center_point();
    check_float_equal(CP.x(), V3.x());
    check_float_equal(CP.y(), V3.y());

    delete curve;
    delete p1;
    delete p2;
  }
  BOOST_AUTO_TEST_CASE(segmentTriIntersection)
  {
    // lets test this new function
    // triangle is in the z-plane
    double p0[] =
      { 0.5, 0.0, 1.0 };
    double p1[] =
      { 0.0, 0.5, 0.0 };

    double t0[] =
      { 1.0, 0.0, 0.5 };
    double t1[] =
      { 0.0, 0.0, 0.5 };
    double t2[] =
      { 0.0, 1.0, 0.5 };

    double ans[] =
      { 0.25, 0.25, 0.5 };
    Vert v0;
    v0.setCoords(3, p0);
    Vert v1;
    v1.setCoords(3, p1);

    Vert f0;
    f0.setCoords(3, t0);
    Vert f1;
    f1.setCoords(3, t1);
    Vert f2;
    f2.setCoords(3, t2);
    double iSec[3] =
      { 0., 0., 0. };
    int qIntersection;
    qIntersection = calcSegmentTriangleIntersection(&v0, &v1, &f0, &f1, &f2, 0.,
						    iSec);
    BOOST_CHECK(qIntersection);
    check_float_equal(dMAG3D_SQ(ans), dMAG3D_SQ(iSec));
    // no intersection
    p1[2] = 1.0;
    v1.setCoords(3, p1);
    qIntersection = calcSegmentTriangleIntersection(&v0, &v1, &f0, &f1, &f2, 0.,
						    iSec);
    BOOST_CHECK(qIntersection == -1);
    // co planar case
    p0[2] = 0.5;
    p1[2] = 0.5;
    v0.setCoords(3, p0);
    v1.setCoords(3, p1);
    qIntersection = calcSegmentTriangleIntersection(&v0, &v1, &f0, &f1, &f2, 0.,
						    iSec);
    BOOST_CHECK(qIntersection == 0);
  }

  BOOST_AUTO_TEST_CASE(calcAreaMoments)
  {
    double adMoment[3];

    double locA[] =
      { 0, 0, 0 };
    double locB[] =
      { 1, 0, 0 };
    double locC[] =
      { 0, 2, 0 };
    double locD[] =
      { 3, 0, 0 };
    double locE[] =
      { 2, sqrt(3), 0 };
    Vert v0, v1, v2;
    v0.setCoords(2, locA);
    v1.setCoords(2, locB);
    v2.setCoords(2, locC);
    calcCovarianceMatrix(&v0, &v1, &v2, adMoment);
    check_float_equal(adMoment[0], 1. / 18.);
    check_float_equal(adMoment[1], -1. / 18.);
    check_float_equal(adMoment[2], 2. / 9.);

    v0.setCoords(2, locB);
    v1.setCoords(2, locD);
    v2.setCoords(2, locE);
    calcCovarianceMatrix(&v0, &v1, &v2, adMoment);
    check_float_equal(adMoment[0], sqrt(3.) / 6.);
    check_float_equal(adMoment[1], 0.);
    check_float_equal(adMoment[2], sqrt(3.) / 6.);
  }

  BOOST_AUTO_TEST_CASE(segmentCircleIntersection)
  {
    // lets test circle/edge intersection
    // circle at origin, radius 1
    // edge goes from (0,0) to a range of points
    double origin[] =
      { 0.0, 0.0 };
    double intersect[2];
    // first point on x-axis
    double p1[] =
      { 3.0, 0.0 };
    bool qIntersect;
    qIntersect = calcSegmentCircleIntersections(origin, p1, origin, 1.0,
						intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at (1,0)
    check_float_equal(intersect[0], 1.0);
    check_float_equal(intersect[1], 0.0);
    // same check, reverse point order
    qIntersect = calcSegmentCircleIntersections(p1, origin, origin, 1.0,
						intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at (1,0)
    check_float_equal(intersect[0], 1.0);
    check_float_equal(intersect[1], 0.0);

    double p2[] =
      { 5.0, 5.0 };
    qIntersect = calcSegmentCircleIntersections(p2, origin, origin, 1.0,
						intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at (sqrt(2)/2,sqrt(2)/2)
    check_float_equal(intersect[0], 1. / sqrt(2.));
    check_float_equal(intersect[1], 1. / sqrt(2.));

    qIntersect = calcSegmentCircleIntersections(origin, p2, origin, 1.0,
						intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at (sqrt(2)/2,sqrt(2)/2)
    check_float_equal(intersect[0], 1. / sqrt(2.));
    check_float_equal(intersect[1], 1. / sqrt(2.));

    // lets check circles in all four quadrants
    double c1[] =
      { -sqrt(3.), -sqrt(3.) };
    qIntersect = calcSegmentCircleIntersections(origin, p2, c1, 6, intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at origin
    check_float_equal(intersect[0], 0.0);
    check_float_equal(intersect[1], 0.0);

    double c2[] =
      { sqrt(3.), -sqrt(3.) };
    qIntersect = calcSegmentCircleIntersections(origin, p2, c2, 6, intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at origin
    check_float_equal(intersect[0], 0.0);
    check_float_equal(intersect[1], 0.0);

    double c3[] =
      { -sqrt(3.), sqrt(3.) };
    qIntersect = calcSegmentCircleIntersections(origin, p2, c3, 6, intersect);
    BOOST_CHECK(qIntersect == true);
    // intersection should be at origin
    check_float_equal(intersect[0], 0.0);
    check_float_equal(intersect[1], 0.0);

    double p3[] =
      { 0.3, 0.4 };
    double p4[] =
      { -0.2, -0.3 };

    qIntersect = calcSegmentCircleIntersections(p3, p4, origin, 1.0, intersect);
    BOOST_CHECK(qIntersect == false);
  }
  BOOST_AUTO_TEST_SUITE_END()
