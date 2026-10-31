/*
 * unitTestQuality.cxx
 *
 *  Created on: Jun 29, 2017
 *      Author: cfog
 */

#define BOOST_TEST_MODULE Quality
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>

#include "GR_Mesh2D.h"
#include "GR_QualMeasure.h"
#include "GR_Vec.h"
#include "GR_VolMesh.h"

struct Quality2DFixture {
  Mesh2D M2D;
  Quality2DFixture() :
      M2D()
  {
    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(0, 0.5);
    Vert *pV3 = M2D.createVert(0.5, 1.5);

    M2D.createTriCell(pV0, pV1, pV2);
    M2D.createTriCell(pV2, pV1, pV3);
  }
  ~Quality2DFixture()
  {
  }
};

struct Quality3DFixture {
  VolMesh VM;
  Quality3DFixture() :
      VM()
  {
    Vert *pV0 = VM.createVert(0, 0, 0);
    Vert *pV1 = VM.createVert(1, 0, 0);
    Vert *pV2 = VM.createVert(0.9, 1, 0);
    Vert *pV3 = VM.createVert(0, 0.9, 0);
    Vert *pV4 = VM.createVert(0, 0, 0.9);
    Vert *pV5 = VM.createVert(0.9, 0, 1);
    Vert *pV6 = VM.createVert(1, 0.9, 0.9);
    Vert *pV7 = VM.createVert(0, 1, 0.9);

    bool exists;
    VM.createTetCell(exists, pV0, pV1, pV2, pV4);
    VM.createTetCell(exists, pV1, pV2, pV4, pV6);
    VM.createTetCell(exists, pV1, pV6, pV4, pV5);
    VM.createTetCell(exists, pV0, pV2, pV3, pV4);
    VM.createTetCell(exists, pV2, pV3, pV4, pV6);
    VM.createTetCell(exists, pV3, pV4, pV6, pV7);
  }
  ~Quality3DFixture()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(Qual2DTests, Quality2DFixture)

  BOOST_AUTO_TEST_CASE(Angle2DTest)
  {
    GRUMMP::Angle2D Angle2;
    BOOST_CHECK_EQUAL(Angle2.maxNumValues(), 4);
    BOOST_CHECK_EQUAL(Angle2.getName(), "angles (2D)");

    // Check first tri
    int nVal = -1;
    double values[4];
    Angle2.evalAllFromCell(nVal, values, M2D.getCell(0));
    double correct0[] =
      { 90 * M_PI / 180, atan(0.5), atan(2) };
    BOOST_CHECK_EQUAL(nVal, 3);
    BOOST_CHECK_CLOSE(values[0], correct0[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[1], correct0[1], 1.e-8);
    BOOST_CHECK_CLOSE(values[2], correct0[2], 1.e-8);

    double result = Angle2.eval(M2D.getCell(0)->getVert(0),
				M2D.getCell(0)->getVert(1),
				M2D.getCell(0)->getVert(2));
    BOOST_CHECK_CLOSE(result, -correct0[0], 1.e-8);

    // Check second tri
    Cell *pC = M2D.getCell(1);
    Angle2.evalAllFromCell(nVal, values, pC);
    double len[] =
      { dDIST2D(pC->getVert(0)->getCoords(), pC->getVert(1)->getCoords()),
	  dDIST2D(pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords()),
	  dDIST2D(pC->getVert(2)->getCoords(), pC->getVert(0)->getCoords()) };
    double cosines[] =
      { (len[0] * len[0] + len[2] * len[2] - len[1] * len[1])
	  / (2 * len[0] * len[2]), (len[1] * len[1] + len[0] * len[0]
	  - len[2] * len[2]) / (2 * len[1] * len[0]), (len[2] * len[2]
	  + len[1] * len[1] - len[0] * len[0]) / (2 * len[2] * len[1]) };
    double correct1[] =
      { GR_acos(cosines[0]), GR_acos(cosines[1]), GR_acos(cosines[2]) };
    BOOST_CHECK_CLOSE(values[0], correct1[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[1], correct1[1], 1.e-8);
    BOOST_CHECK_CLOSE(values[2], correct1[2], 1.e-8);
    result = Angle2.eval(M2D.getCell(1)->getVert(0), M2D.getCell(1)->getVert(1),
			 M2D.getCell(1)->getVert(2));
    BOOST_CHECK_CLOSE(result, -correct1[2], 1.e-8);

    // Quad check
    Vert* pV0 = M2D.getVert(0);
    Vert* pV1 = M2D.getVert(1);
    Vert* pV2 = M2D.getVert(2);
    Vert* pV3 = M2D.getVert(3);

    M2D.deleteCell(M2D.getCell(0));
    M2D.deleteCell(M2D.getCell(1));
    M2D.allowNonSimplicial();
    Cell *quad = M2D.createQuadCell(pV0, pV1, pV3, pV2);
    Angle2.evalAllFromCell(nVal, values, quad);
    BOOST_CHECK_EQUAL(nVal, 4);
    BOOST_CHECK_CLOSE(values[0], correct0[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[1], correct0[1] + correct1[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[2], correct1[1], 1.e-8);
    BOOST_CHECK_CLOSE(values[3], correct0[2] + correct1[2], 1.e-8);

    result = Angle2.eval(pV0, pV1, pV2, pV3);
    BOOST_CHECK_CLOSE(result, -M_PI/2, 1.e-8);
  }

  BOOST_AUTO_TEST_CASE(EdgeRatio2DTest)
  {
    GRUMMP::EdgeLengthToRadius2D ER2;
    BOOST_CHECK_EQUAL(ER2.maxNumValues(), 3);
    BOOST_CHECK_EQUAL(ER2.getName(), "EdgeLength to Radius (2D)");

    Vert* pV0 = M2D.getVert(0);
    Vert* pV1 = M2D.getVert(1);
    Vert* pV2 = M2D.getVert(2);
    Vert* pV3 = M2D.getVert(3);

    double len01 = dDIST2D(pV0->getCoords(), pV1->getCoords());
    double len12 = dDIST2D(pV1->getCoords(), pV2->getCoords());
    double len20 = dDIST2D(pV2->getCoords(), pV0->getCoords());
    double len13 = dDIST2D(pV1->getCoords(), pV3->getCoords());
    double len32 = dDIST2D(pV3->getCoords(), pV2->getCoords());

    // Circumradius is trivial for right triangles.
    double circradius0 = 0.5 * len12;
    double circradius1 = 0.5 * len13;
    // Check first tri
    int nVal = -1;
    double values[] =
      { -10, -20, -30, -40 };
    ER2.evalAllFromCell(nVal, values, M2D.getCell(0));
    double correct0[] =
      { len12 / circradius0, len20 / circradius0, len01 / circradius0 };
    BOOST_CHECK_EQUAL(nVal, 3);
    BOOST_CHECK_CLOSE(values[0], correct0[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[1], correct0[1], 1.e-8);
    BOOST_CHECK_CLOSE(values[2], correct0[2], 1.e-8);

    double result = ER2.eval(M2D.getCell(0)->getVert(0),
			     M2D.getCell(0)->getVert(1),
			     M2D.getCell(0)->getVert(2));
    BOOST_CHECK_CLOSE(result, correct0[1], 1.e-8);

    // Check second tri
    Cell *pC = M2D.getCell(1);
    ER2.evalAllFromCell(nVal, values, pC);
    double correct1[] =
      { len12 / circradius1, len32 / circradius1, len13 / circradius1 };
    BOOST_CHECK_CLOSE(values[0], correct1[0], 1.e-8);
    BOOST_CHECK_CLOSE(values[1], correct1[1], 1.e-8);
    BOOST_CHECK_CLOSE(values[2], correct1[2], 1.e-8);
    result = ER2.eval(M2D.getCell(1)->getVert(0), M2D.getCell(1)->getVert(1),
		      M2D.getCell(1)->getVert(2));
    BOOST_CHECK_CLOSE(result, correct1[0], 1.e-8);

    // Can't compute this for quads
  }

  BOOST_AUTO_TEST_CASE(RadiusRatio2DTest)
  {
    GRUMMP::RadiusRatio2D RR2;
    BOOST_CHECK_EQUAL(RR2.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(RR2.getName(), "InToCircumRadius2D");

    Vert* pV0 = M2D.getVert(0);
    Vert* pV1 = M2D.getVert(1);
    Vert* pV2 = M2D.getVert(2);
    Vert* pV3 = M2D.getVert(3);

    double correct0 = 0.68328157299975, correct1 = 0.82842712474619;
    int nVal = -1;
    double values[] =
      { -10, -20, -30, -40 };
    RR2.evalAllFromCell(nVal, values, M2D.getCell(0));
    BOOST_CHECK_EQUAL(nVal, 1);
    BOOST_CHECK_CLOSE(values[0], correct0, 1.e-8);

    double result = RR2.eval(pV0, pV1, pV2);
    BOOST_CHECK_CLOSE(result, correct0, 1.e-8);

    // Check second tri
    Cell *pC = M2D.getCell(1);
    RR2.evalAllFromCell(nVal, values, pC);
    BOOST_CHECK_CLOSE(values[0], correct1, 1.e-8);
    result = RR2.eval(pV2, pV1, pV3);
    BOOST_CHECK_CLOSE(result, correct1, 1.e-8);

    // Can't compute this for quads
  }

  BOOST_AUTO_TEST_CASE(AreaPerim2DTest)
  {
    GRUMMP::AreaPerim2D AP2;
    BOOST_CHECK_EQUAL(AP2.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(AP2.getName(), "AreaPerim2D");

    Vert* pV0 = M2D.getVert(0);
    Vert* pV1 = M2D.getVert(1);
    Vert* pV2 = M2D.getVert(2);
    Vert* pV3 = M2D.getVert(3);

    double correct0 = 12 * sqrt(3) / (14 + 6 * sqrt(5)), correct1 = 6 * sqrt(3)
	/ (6 + 4 * sqrt(2));
    int nVal = -1;
    double values[] =
      { -10, -20, -30, -40 };
    AP2.evalAllFromCell(nVal, values, M2D.getCell(0));
    BOOST_CHECK_EQUAL(nVal, 1);
    BOOST_CHECK_CLOSE(values[0], correct0, 1.e-8);

    double result = AP2.eval(pV0, pV1, pV2);
    BOOST_CHECK_CLOSE(result, correct0, 1.e-8);

    // Check second tri
    Cell *pC = M2D.getCell(1);
    AP2.evalAllFromCell(nVal, values, pC);
    BOOST_CHECK_CLOSE(values[0], correct1, 1.e-8);
    result = AP2.eval(pV2, pV1, pV3);
    BOOST_CHECK_CLOSE(result, correct1, 1.e-8);

    // Can't compute this for quads
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(Quality3D, Quality3DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    BOOST_CHECK(VM.getCell(0)->calcSize() > 0);
    BOOST_CHECK(VM.getCell(1)->calcSize() > 0);
    BOOST_CHECK(VM.getCell(2)->calcSize() > 0);
    BOOST_CHECK(VM.getCell(3)->calcSize() > 0);
    BOOST_CHECK(VM.getCell(4)->calcSize() > 0);
    BOOST_CHECK(VM.getCell(5)->calcSize() > 0);
  }

  BOOST_AUTO_TEST_CASE(Dihedrals)
  {
    GRUMMP::DihedralAngles DA;
    BOOST_CHECK_EQUAL(DA.maxNumValues(), 12);
    BOOST_CHECK_EQUAL(DA.getName(), "dihedral angles");

    // Actual edge order is:
    // {01,12,20,04,14,24},
    // {12,24,41,16,26,46},
    // {16,64,41,51,56,54},
    double correctTet[6][6] =
      {
	{ M_PI / 2, 0.837981225, M_PI / 2, 1.1023762578, 0.735289577,
	    1.503999377 },
	{ 0.7059924800, 1.9970589508, 0.9346840109, 0.79484394698,
	    1.85160824327, 0.93102944265 },
	{ 0.9316007967, 1.6699735878, 0.8635542992, 1.44870491451,
	    1.09619235397, 1.54788252535 },
	{ M_PI / 2, 0.732815101786, M_PI / 2, 1.64920268168, 0.78846566730,
	    1.012436845623 },
	{ 0.48354711316367621, 2.3167639604733306, 0.48834597094593613,
	    0.61685308269162786, 2.374583414110516, 0.66023924720871774 },
	{ 1.0799867567747086, 1.4717332817279143, 1.004044000957008,
	    1.4606838371708948, 0.93160079672223628, 1.57079632679489660 } };
    // Actual edge order is:
    // {02,26,64,40,01,21,41,61},
    // {...
//    double correctPyr[2][8] =
//      {
//	{ correctTet[0][2], correctTet[1][4], correctTet[1][5],
//	    correctTet[0][3], correctTet[0][0], correctTet[0][1]
//		+ correctTet[1][0], correctTet[0][4] + correctTet[1][3],
//	    correctTet[1][3] },
//	{ 0, 0, 0, 0, 0, 0, 0, 0 } };

    double correctPrism[1][9] =
      { M_PI / 2, M_PI / 2, 1.4923899719029439, 1.4487049145113224,
	  1.4609661680267003, 1.813561639797258, 0.7884656673026047, M_PI / 2,
	  0.7884656673026047 };
    
    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      DA.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 6);
      for (int ii = 0; ii < 6; ii++) {
	BOOST_CHECK_CLOSE(values[ii], correctTet[iC][ii], 1.e-8);
      }
    }
    double result = DA.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, -M_PI/2, 1.e-8);

    // A prism check
    evalAllDihedsFromCoords(nVals, values, VM.getVert(0)->getCoords(),
		      VM.getVert(1)->getCoords(), VM.getVert(2)->getCoords(),
		      VM.getVert(4)->getCoords(), VM.getVert(5)->getCoords(),
		      VM.getVert(6)->getCoords());
    BOOST_CHECK_EQUAL(nVals, 9);
    for (int ii = 0; ii < 9; ii++) {
      BOOST_CHECK_CLOSE(values[ii], correctPrism[0][ii], 1.e-8);
    }

    // Waiting for pyramid dihedral angle calculation.
//  DA.evalAll(nVals, values, VM.getVert(0), VM.getVert(2), VM.getVert(6), VM.getVert(4), VM.getVert(1));
//  BOOST_CHECK_EQUAL(nVals, 8);
//  for (int ii = 0; ii < 8; ii++) {
//    BOOST_CHECK_CLOSE(values[ii], correctPyr[0][ii], 1.e-8);
//  }
  }

  BOOST_AUTO_TEST_CASE(Solids)
  {
    GRUMMP::SolidAngles SA;
    BOOST_CHECK_EQUAL(SA.maxNumValues(), 8);
    BOOST_CHECK_EQUAL(SA.getName(), "solid angles");

    // Actual edge order is:
    // {01,12,20,04,14,24},
    // {12,24,41,16,26,46},
    // {16,64,41,51,56,54},
    double dihedTet[6][6] =
      {
	{ M_PI / 2, 0.837981225, M_PI / 2, 1.1023762578, 0.735289577,
	    1.503999377 },
	{ 0.7059924800, 1.9970589508, 0.9346840109, 0.79484394698,
	    1.85160824327, 0.93102944265 },
	{ 0.9316007967, 1.6699735878, 0.8635542992, 1.44870491451,
	    1.09619235397, 1.54788252535 },
	{ M_PI / 2, 0.732815101786, M_PI / 2, 1.64920268168, 0.78846566730,
	    1.012436845623 },
	{ 0.48354711316367621, 2.3167639604733306, 0.48834597094593613,
	    0.61685308269162786, 2.374583414110516, 0.66023924720871774 },
	{ 1.0799867567747086, 1.4717332817279143, 1.004044000957008,
	    1.4606838371708948, 0.93160079672223628, 1.57079632679489660 } };
    double correctTet[6][4];
    for (int ii = 0; ii < 6; ii++) {
      correctTet[ii][0] = dihedTet[ii][0] + dihedTet[ii][1] + dihedTet[ii][2]
	  - M_PI;
      correctTet[ii][1] = dihedTet[ii][0] + dihedTet[ii][3] + dihedTet[ii][4]
	  - M_PI;
      correctTet[ii][2] = dihedTet[ii][1] + dihedTet[ii][3] + dihedTet[ii][5]
	  - M_PI;
      correctTet[ii][3] = dihedTet[ii][2] + dihedTet[ii][4] + dihedTet[ii][5]
	  - M_PI;
    }

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      SA.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 4);
      for (int ii = 0; ii < 4; ii++) {
	BOOST_CHECK_CLOSE(values[ii], correctTet[iC][ii], 1.e-7);
      }
    }
    double result = SA.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correctTet[0][1], 1.e-7);

    // Waiting for pyramid dihedral angle calculation.
//  SA.evalAll(nVals, values, VM.getVert(0), VM.getVert(2), VM.getVert(6), VM.getVert(4), VM.getVert(1));
//  BOOST_CHECK_EQUAL(nVals, 8);
//  for (int ii = 0; ii < 8; ii++) {
//    BOOST_CHECK_CLOSE(values[ii], correctPyr[0][ii], 1.e-8);
//  }
  }

  BOOST_AUTO_TEST_CASE(SineDihedrals)
  {
    GRUMMP::SineDihedralAngles SD;
    BOOST_CHECK_EQUAL(SD.maxNumValues(), 12);
    BOOST_CHECK_EQUAL(SD.getName(), "sine of dihedral angles");

    // Actual edge order is:
    // {01,12,20,04,14,24},
    // {12,24,41,16,26,46},
    // {16,64,41,51,56,54},
    double correctTetDihedrals[6][6] =
      {
	{ M_PI / 2, 0.837981225, M_PI / 2, 1.1023762578, 0.735289577,
	    1.503999377 },
	{ 0.7059924800, 1.9970589508, 0.9346840109, 0.79484394698,
	    1.85160824327, 0.93102944265 },
	{ 0.9316007967, 1.6699735878, 0.8635542992, 1.44870491451,
	    1.09619235397, 1.54788252535 },
	{ M_PI / 2, 0.732815101786, M_PI / 2, 1.64920268168, 0.78846566730,
	    1.012436845623 },
	{ 0.48354711316367621, 2.3167639604733306, 0.48834597094593613,
	    0.61685308269162786, 2.374583414110516, 0.66023924720871774 },
	{ 1.0799867567747086, 1.4717332817279143, 1.004044000957008,
	    1.4606838371708948, 0.93160079672223628, 1.57079632679489660 } };
    // Actual edge order is:
    // {02,26,64,40,01,21,41,61},
    // {...
//    double correctPyr[2][8] =
//      {
//	{ correctTetDihedrals[0][2], correctTetDihedrals[1][4],
//	    correctTetDihedrals[1][5], correctTetDihedrals[0][3],
//	    correctTetDihedrals[0][0], correctTetDihedrals[0][1]
//		+ correctTetDihedrals[1][0], correctTetDihedrals[0][4]
//		+ correctTetDihedrals[1][3], correctTetDihedrals[1][3] },
//	{ 0, 0, 0, 0, 0, 0, 0, 0 } };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      SD.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 6);
      for (int ii = 0; ii < 6; ii++) {
	BOOST_CHECK_CLOSE(values[ii], sin(correctTetDihedrals[iC][5 - ii]),
			  1.e-8);
      }
    }
    double result = SD.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, sin(correctTetDihedrals[0][4]), 1.e-8);

    // Waiting for pyramid dihedral angle calculation.
//  SD.evalAll(nVals, values, VM.getVert(0), VM.getVert(2), VM.getVert(6), VM.getVert(4), VM.getVert(1));
//  BOOST_CHECK_EQUAL(nVals, 8);
//  for (int ii = 0; ii < 8; ii++) {
//    BOOST_CHECK_CLOSE(values[ii], correctPyr[0][ii], 1.e-8);
//  }
  }

  BOOST_AUTO_TEST_CASE(VolumeLength)
  {
    GRUMMP::VolumeLength3D VL;
    BOOST_CHECK_EQUAL(VL.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(VL.getName(), "VolumeLength3D");

    double correct[6] =
      { 0.31405059088265896, 0.3250784746161226, 0.23337118616950436,
	  0.38749749901838848, 0.54713166425639648, 0.18403020852829255 };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      VL.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 1);
      BOOST_CHECK_CLOSE(values[0], correct[iC], 1.e-8);
    }
    double result = VL.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correct[0], 1.e-8);
    // Only applies to tets
  }

  BOOST_AUTO_TEST_CASE(VolumeArea)
  {
    GRUMMP::VolumeArea3D VA;
    BOOST_CHECK_EQUAL(VA.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(VA.getName(), "VolumeArea3D");

    double correct[6] =
      { 0.87148248854320309, 0.71573528517099227, 0.98158660716477097,
	  0.83651488997745305, 0.44674477811229102, 0.98095845755521627 };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      VA.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 1);
      BOOST_CHECK_CLOSE(values[0], correct[iC], 1.e-8);
    }
    double result = VA.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correct[0], 1.e-8);
    // Only applies to tets
  }

  BOOST_AUTO_TEST_CASE(MarcumSkewness)
  {
    GRUMMP::MarcumSkewness3D Skewness;
    BOOST_CHECK_EQUAL(Skewness.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(Skewness.getName(), "MarcumSkewness3D");

    double correct[6] =
      { 0.54570447473151784, 0.37749235054488989, 0.77232504403827884,
	  0.43953012811852188, 0.14453676704752622, 0.79888092010587752 };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      Skewness.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 1);
      BOOST_CHECK_CLOSE(values[0], correct[iC], 1.e-8);
    }
    double result = Skewness.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
				  VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correct[0], 1.e-8);
    // Only applies to tets
  }

  BOOST_AUTO_TEST_CASE(RadiusRatio3D)
  {
    GRUMMP::RadiusRatio3D RR3D;
    BOOST_CHECK_EQUAL(RR3D.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(RR3D.getName(), "RadiusRatio3D");

    double correct[6] =
      { 0.74557547677691527, 0.57827933189590119, 0.90618970111208574,
	  0.67501278782823537, 0.30668963737290167, 0.9160681003329183 };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      RR3D.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 1);
      BOOST_CHECK_CLOSE(values[0], correct[iC], 1.e-8);
    }
    RR3D.evalAllFromVerts(nVals, values, VM.getCell(0)->getVert(1),
			  VM.getCell(0)->getVert(2), VM.getCell(0)->getVert(3),
			  VM.getCell(0)->getVert(0));
    BOOST_CHECK_EQUAL(nVals, 1);
    BOOST_CHECK_CLOSE(values[0], correct[0], 1.e-8);

    double result = RR3D.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			      VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correct[0], 1.e-8);
    // Only applies to tets
  }

  BOOST_AUTO_TEST_CASE(ShortestToRadius3D)
  {
    GRUMMP::ShortestToRadius3D SR;
    BOOST_CHECK_EQUAL(SR.maxNumValues(), 1);
    BOOST_CHECK_EQUAL(SR.getName(), "ShortestToRadius3D");

    double correct[6] =
      { 0.67864477855027083, 0.68714782447315104, 0.68773193409492628,
	  0.67809934154391371, 0.68181265760533505, 0.70594550896393704 };

    int nVals = -1;
    double values[12];
    for (int iC = 0; iC < 6; iC++) {
      SR.evalAllFromCell(nVals, values, VM.getCell(iC));
      BOOST_CHECK_EQUAL(nVals, 1);
      BOOST_CHECK_CLOSE(values[0], correct[iC], 1.e-8);
    }
    double result = SR.eval(VM.getVert(0), VM.getVert(1), VM.getVert(2),
			    VM.getVert(4));
    BOOST_CHECK_CLOSE(result, correct[0], 1.e-8);
    // Only applies to tets
  }

  BOOST_AUTO_TEST_SUITE_END()

