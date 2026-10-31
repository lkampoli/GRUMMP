#include <vector>
#define BOOST_TEST_MODULE Mesh2D
#include <boost/test/unit_test.hpp>

#include "Mesh2DFixtureBase.h"
#include "GR_ADT.h"
#include "GR_Geometry.h"
#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "GR_Mesh2D.h"
#include "GR_SwapDecider.h"
#include "GR_QualMeasure.h"
#include "GR_InsertionManager.h"

#include "CubitVector.hpp"

struct MeshMod2DFixture : public Mesh2DFixtureBase {
  MeshMod2DFixture() :
      Mesh2DFixtureBase()
  {
    M2D.disallowSwapRecursion();
    M2D.allowNonSimplicial();
    // First, create some verts.  The first six can be used to create
    // two quads; verts 1, 2, 6, 7 to create two independent triangles.
    // Sketch:
    //
    //
    //
    //
    //       5 --5-- 4
    //       |       |
    //       6       4
    //       |       |
    //       3 --2-- 2 --9-- 7
    //       |       |     / |
    //       3       1   10  8
    //       |       | /     |
    //       0 --0-- 1 --7-- 6

    // lets perturb point 5 to make this slightly more interesting.
    // should form two 60 degree angles, a 150 one, and a 90
    // D. Zaide, Oct 2014

    BOOST_CHECK(M2D.isValid());
    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(1, 1);
    Vert *pV3 = M2D.createVert(0, 1);
    Vert *pV4 = M2D.createVert(1, 2);
    Vert *pV5 = M2D.createVert(1. / (1. + sqrt(3.)), 2 - 1. / (1. + sqrt(3.)));
    Vert *pV6 = M2D.createVert(2, 0);
    Vert *pV7 = M2D.createVert(2.1, 1.1); // So that the swap will happen

    (void) M2D.createQuadCell(pV0, pV1, pV2, pV3, 1);
    (void) M2D.createTriCell(pV1, pV6, pV7, 1);

    // The first face already exists.
    Face *pF2 = findCommonFace(pV2, pV3);
    // These three are new.
    bool isNew = false;
    Face *pF4 = M2D.createFace(isNew, pV2, pV4);
    Face *pF5 = M2D.createFace(isNew, pV4, pV5);
    Face *pF6 = M2D.createFace(isNew, pV5, pV3);

    // Creation from faces out of canonical order
    (void) M2D.createQuadCell(pF2, pF4, pF5, pF6, 1);

    // These two exist
    Face *pF10 = findCommonFace(pV1, pV7);
    Face *pF1 = findCommonFace(pV1, pV2);
    // This one is new
    Face *pF9 = M2D.createFace(isNew, pV2, pV7);

    (void) M2D.createTriCell(pF1, pF10, pF9, 1);

    Face *pF0 = findCommonFace(pV1, pV0);
    Face *pF3 = findCommonFace(pV0, pV3);
    Face *pF7 = findCommonFace(pV1, pV6);
    Face *pF8 = findCommonFace(pV6, pV7);
    M2D.createBFace(pF0);
    M2D.createBFace(pF3);
    M2D.createBFace(pF4);
    M2D.createBFace(pF5);
    // allows it to cover an additional function in creating the BFace
    M2D.createBFace(pF6, 1);
    M2D.createBFace(pF7, 2);
    M2D.createBFace(pF8, 3);
    M2D.createBFace(pF9, 4);
  }
  ~MeshMod2DFixture()
  {
  }

  void
  deleteTris()
  {
    // Deletes the triangles, but not the bdry faces around them.
    Face *pF = findCommonFace(M2D.getCell(0), M2D.getCell(1));

    M2D.deleteCell(M2D.getCell(0));
    M2D.deleteCell(M2D.getCell(1));

    BOOST_CHECK(M2D.getCell(0)->isDeleted());
    BOOST_CHECK(M2D.getCell(1)->isDeleted());
    BOOST_CHECK_EQUAL(pCInvalidCell, pF->getLeftCell());
    BOOST_CHECK_EQUAL(pCInvalidCell, pF->getRightCell());
    BOOST_CHECK(pF->isDeleted());
    checkMeshSize(8, 10, 8, 2, 0, 2);

  }
};

BOOST_FIXTURE_TEST_SUITE(MeshMod2D, MeshMod2DFixture)

  BOOST_AUTO_TEST_CASE(Setup)
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(8, 11, 8, 4, 2, 2);
  }

  BOOST_AUTO_TEST_CASE(VertexDegree)
  {
    BOOST_CHECK_EQUAL(2, M2D.getVert(0)->getNumFaces());
    BOOST_CHECK_EQUAL(4, M2D.getVert(1)->getNumFaces());
    BOOST_CHECK_EQUAL(4, M2D.getVert(2)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(3)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(4)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(5)->getNumFaces());
    BOOST_CHECK_EQUAL(2, M2D.getVert(6)->getNumFaces());
    BOOST_CHECK_EQUAL(3, M2D.getVert(7)->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(Orient2D)
  {
    Vert *pV0 = M2D.getVert(0);
    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV4 = M2D.getVert(4);

    BOOST_CHECK_EQUAL(1, checkOrient2D(pV0, pV1, pV2));
    BOOST_CHECK_EQUAL(-1, checkOrient2D(pV0, pV2, pV1));
    BOOST_CHECK_EQUAL(0, checkOrient2D(pV1, pV2, pV4));
  }

  BOOST_AUTO_TEST_CASE(Connect)
  {
    Vert *pV0 = M2D.getVert(0);
    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV3 = M2D.getVert(3);
    Vert *pV4 = M2D.getVert(4);
    Vert *pV5 = M2D.getVert(5);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);

    Cell *pCA = M2D.getCell(0);
    Cell *pCB = M2D.getCell(1);
    Cell *pCC = M2D.getCell(2);
    Cell *pCD = M2D.getCell(3);

    // Tris always come before quads, regardless of creation order.
    BOOST_CHECK_EQUAL(pCA->getType(), CellSkel::eTriCell);
    BOOST_CHECK_EQUAL(pCB->getType(), CellSkel::eTriCell);
    BOOST_CHECK_EQUAL(pCC->getType(), CellSkel::eQuadCell);
    BOOST_CHECK_EQUAL(pCD->getType(), CellSkel::eQuadCell);

    // Check orientation and canonical ordering for triangles
    checkCanonicalTri(pCA, pV1, pV6, pV7);
    checkCanonicalTri(pCB, pV1, pV7, pV2);

    // Check orientation and canonical ordering for quads
    checkCanonicalQuad(pCC, pV0, pV1, pV2, pV3);
    checkCanonicalQuad(pCD, pV3, pV2, pV4, pV5);
  }

  BOOST_AUTO_TEST_CASE(Delete)
  {
    deleteTris();

    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);

    BOOST_CHECK_EQUAL(3, pV1->getNumFaces());
    BOOST_CHECK_EQUAL(4, pV2->getNumFaces());
    BOOST_CHECK_EQUAL(2, pV6->getNumFaces());
    BOOST_CHECK_EQUAL(2, pV7->getNumFaces());
  }

  BOOST_AUTO_TEST_CASE(FakeSwap)
  {
    deleteTris();

    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);

    Cell *pCA = M2D.createTriCell(pV1, pV6, pV2, 1);
    Cell *pCB = M2D.createTriCell(pV2, pV6, pV7, 1);

    checkMeshSize(8, 11, 8, 4, 2, 2);

    BOOST_CHECK_EQUAL(3, pV1->getNumFaces());
    BOOST_CHECK_EQUAL(5, pV2->getNumFaces());
    BOOST_CHECK_EQUAL(3, pV6->getNumFaces());
    BOOST_CHECK_EQUAL(2, pV7->getNumFaces());

    Face *pF = findCommonFace(pCA, pCB);
    BOOST_CHECK(pF->hasVert(pV2));
    BOOST_CHECK(pF->hasVert(pV6));
  }

//BOOST_AUTO_TEST_CASE(Purge)
//{
//  deleteTris();
//
//  Vert *pV1 = M2D.getVert(1);
//  Vert *pV2 = M2D.getVert(2);
//  Vert *pV6 = M2D.getVert(6);
//  Vert *pV7 = M2D.getVert(7);
//
//  (void) M2D.createTriCell(pV1, pV6, pV2, 1);
//  (void) M2D.createTriCell(pV2, pV6, pV7, 1);
//
//  M2D.purgeAllEntities();
//  BOOST_CHECK(M2D.isValid());
//}

  BOOST_AUTO_TEST_CASE(Reconfig)
  {
    Cell *pCA = M2D.getCell(0);
    Cell *pCB = M2D.getCell(1);

    Face *pF = findCommonFace(pCA, pCB);
    int iRes = M2D.reconfigure(pF);
    BOOST_CHECK_EQUAL(1, iRes);

    checkMeshSize(8, 11, 8, 4, 2, 2);

    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);

    BOOST_CHECK_EQUAL(3, pV1->getNumFaces());
    BOOST_CHECK_EQUAL(5, pV2->getNumFaces());
    BOOST_CHECK_EQUAL(3, pV6->getNumFaces());
    BOOST_CHECK_EQUAL(2, pV7->getNumFaces());

    // The old tri's were deleted and replaced with two others in the same
    // memory slots, so use indices 0 and 1 here.
    pCA = M2D.getCell(0);
    pCB = M2D.getCell(1);
    pF = findCommonFace(pCA, pCB);
    BOOST_CHECK(pF->hasVert(pV2));
    BOOST_CHECK(pF->hasVert(pV6));
  }

  BOOST_AUTO_TEST_CASE(Swap)
  {
    Cell *pCA = M2D.getCell(0);
    Cell *pCB = M2D.getCell(1);

    Face *pF = findCommonFace(pCA, pCB);

    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);

    //  BOOST_CHECK(M2D.qDoSwap(pV1, pV7, pV6, pV2));

    int iRes = M2D.iFaceSwap_deprecated(pF);
    BOOST_CHECK_EQUAL(1, iRes);

    checkMeshSize(8, 11, 8, 4, 2, 2);

    BOOST_CHECK_EQUAL(3, pV1->getNumFaces());
    BOOST_CHECK_EQUAL(5, pV2->getNumFaces());
    BOOST_CHECK_EQUAL(3, pV6->getNumFaces());
    BOOST_CHECK_EQUAL(2, pV7->getNumFaces());

    // The old tri's were deleted and replaced with two others in the same
    // memory slots, so use indices 0 and 1 here.
    pCA = M2D.getCell(0);
    pCB = M2D.getCell(1);
    pF = findCommonFace(pCA, pCB);
    BOOST_CHECK(pF->hasVert(pV2));
    BOOST_CHECK(pF->hasVert(pV6));
  }

  BOOST_AUTO_TEST_CASE(FakeInsert)
  {
    deleteTris();

    Vert *pV1 = M2D.getVert(1);
    Vert *pV2 = M2D.getVert(2);
    Vert *pV6 = M2D.getVert(6);
    Vert *pV7 = M2D.getVert(7);
    Vert *pV8 = M2D.createVert(1.5, 0.5);

    Cell *pT0 = M2D.createTriCell(pV8, pV1, pV6, 1);
    Cell *pT1 = M2D.createTriCell(pV8, pV6, pV7, 1);
    Cell *pT2 = M2D.createTriCell(pV8, pV7, pV2, 1);
    Cell *pT3 = M2D.createTriCell(pV8, pV2, pV1, 1);

    BOOST_CHECK(M2D.isValid());

    checkCanonicalTri(pT0, pV8, pV1, pV6);
    checkCanonicalTri(pT1, pV8, pV6, pV7);
    checkCanonicalTri(pT2, pV8, pV7, pV2);
    checkCanonicalTri(pT3, pV8, pV2, pV1);

    checkMeshSize(9, 14, 8, 6, 4, 2);

    BOOST_CHECK_EQUAL(4, pV1->getNumFaces());
    BOOST_CHECK_EQUAL(5, pV2->getNumFaces());
    BOOST_CHECK_EQUAL(3, pV6->getNumFaces());
    BOOST_CHECK_EQUAL(3, pV7->getNumFaces());
    BOOST_CHECK_EQUAL(4, pV8->getNumFaces());

  }
  BOOST_AUTO_TEST_CASE(TriAndQuadAngles)
  {
    // gets the quad cell with the perturbed angle
    Cell * quad = M2D.getCell(M2D.getNumTriCells() + 1);
    BOOST_CHECK(quad->doFullCheck());
    double adCent[2], adDihed[4];
    int nDihed;
    // make sure we have the correct quad
    quad->calcCentroid(adCent);
    BOOST_CHECK(iFuzzyComp(adCent[0], 0.5915063509461097) == 0);
    BOOST_CHECK(iFuzzyComp(adCent[1], 1.4084936490538904) == 0);

    quad->calcAllDihed(adDihed, &nDihed);
    BOOST_CHECK(nDihed == 4);
    BOOST_CHECK(iFuzzyComp(adDihed[0], 60.) == 0);
    BOOST_CHECK(iFuzzyComp(adDihed[1], 90.) == 0);
    BOOST_CHECK(iFuzzyComp(adDihed[2], 60.) == 0);
    BOOST_CHECK(iFuzzyComp(adDihed[3], 150.) == 0);

  }
  BOOST_AUTO_TEST_CASE(MergingTrisIntoQuads)
  {
    int nQuads = M2D.getNumQuadCells();

    Cell * cell0 = M2D.getCell(0);
    Cell * cell1 = M2D.getCell(1);
    BOOST_CHECK(cell0->getType() == Cell::eTriCell);
    BOOST_CHECK(cell1->getType() == Cell::eTriCell);

    // lets attempt to merge these into a quadCell
    bool qQuad = M2D.createQuadFromTris(cell0, cell1);
    BOOST_CHECK(qQuad);
    BOOST_CHECK(M2D.isValid());
    BOOST_CHECK_EQUAL(nQuads + 1, M2D.getNumQuadCells());

    // now lets attempt to undo the merge
    bool qTris = M2D.createTrisFromQuad(M2D.getCell(2 + M2D.getNumTriCells()));
    BOOST_CHECK(qTris);
    BOOST_CHECK_EQUAL(nQuads, M2D.getNumQuadsInUse());
    BOOST_CHECK(M2D.isValid());

  }

//BOOST_AUTO_TEST_CASE(CurveGeomInternalBoundaryEdge)
//{
//// lets create an internal boundary edge between vertex 1 and 2
//	  Vert *pV1 = M2D.getVert(1);
//	  Vert *pV2 = M2D.getVert(2);
//	  CubitVector V1(pV1->getCoords());
//	  CubitVector V2(pV2->getCoords());
//	  GRPoint * p1 = new GRPoint(V1);
//	  GRPoint * p2 = new GRPoint(V2);
//	  GRLine* line = new GRLine(V1, V2);
//
//	  GRCurve* curve = new GRCurve(line, p1,p2);
//	  DLIList<CubitString*> string_list;
//	  DLIList<double> double_list;
//	  DLIList<int> int_list;
//
//	  CubitString cond_string("CURVE_CONDITION");
//	  CubitString reg_left( "REG_LEFT" ); CubitString bdry_left( "BDRY_LEFT" );
//	  CubitString reg_right("REG_RIGHT"); CubitString bdry_right("BDRY_RIGHT");
//
//	  string_list.append(&cond_string); double_list.append(0.); int_list.append(0);
//
//	  // both regions set to 1, no need to try anything else funny.
//	  string_list.append(&reg_left); double_list.append(0.); int_list.append(1);
//	  string_list.append(&bdry_left); double_list.append(0.); int_list.append(0);
//	  string_list.append(&reg_right); double_list.append(0.); int_list.append(1);
//	  string_list.append(&bdry_right); double_list.append(0.); int_list.append(0);
//
//	  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
//	  curve->append_simple_attribute_virt(&csa);
//
//	  BOOST_CHECK(curve->qValid());
//	  GRCurveGeom * cg = curve->get_geometry();
//
//	  IntBdryEdge * IBE = M2D.createIntBdryEdge(pV1,pV2,curve,cg->min_param(),cg->max_param());
//	  BOOST_CHECK(IBE->isValid());
//	  BOOST_CHECK_EQUAL(2,IBE->getNumFaces());
//	  BOOST_CHECK(IBE->doFullCheck());
//	  BOOST_CHECK(IBE->hasVert(pV1));
//	  BOOST_CHECK(IBE->hasVert(pV2));
//}
  BOOST_AUTO_TEST_SUITE_END()

#ifdef HAVE_MESQUITE

struct TMOPAnisoQual : public Mesh2DFixtureBase {
  TMOPAnisoQual() :
      Mesh2DFixtureBase()
  {

    M2D.disallowSwapRecursion();
    M2D.allowNonSimplicial();
    /*
     *Sketch: Mesh with 8 cells.
     *
     *	v6 ---- f11 --- v7 ---- f9 ---- v8
     *	|   \        c5  |  c4       /   |
     *	|     \          |         /     |
     *   f13     f12     f10     f8      f7
     *   |          \     |     /         |
     *	| c6         \   |   /     c3    |
     *   v3 ---- f14 --- v4 ---- f6 ---- v5
     *   |  c7         /  |   \     c2    |
     *   |          /     |     \         |
     *  f15      f2      f1      f4      f5
     *   |      /         |         \     |
     *	|   /       c0   |  c1       \   |
     *   v0 ---- f0 ---- v1 ---- f3 ---- v2
     */

    BOOST_CHECK(M2D.isValid());

    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(2, 0);
    Vert *pV3 = M2D.createVert(0, 0.1);
    Vert *pV4 = M2D.createVert(1, 0.1);
    Vert *pV5 = M2D.createVert(2, 0.1);
    Vert *pV6 = M2D.createVert(0, 0.2);
    Vert *pV7 = M2D.createVert(1, 0.2);
    Vert *pV8 = M2D.createVert(2, 0.2);

    (void) M2D.createTriCell(pV0, pV1, pV4, 1);

    Face *pF0 = findCommonFace(pV0, pV1);
    Face *pF1 = findCommonFace(pV1, pV4);
    Face *pF2 = findCommonFace(pV4, pV0);

    bool bAlreadyExisted = false;
    Face *pF3 = M2D.createFace(bAlreadyExisted, pV1, pV2);
    Face *pF4 = M2D.createFace(bAlreadyExisted, pV2, pV4);
    (void) M2D.createTriCell(pF1, pF3, pF4, 1);

    Face *pF5 = M2D.createFace(bAlreadyExisted, pV2, pV5);
    Face *pF6 = M2D.createFace(bAlreadyExisted, pV5, pV4);
    (void) M2D.createTriCell(pF4, pF5, pF6, 1);

    Face *pF7 = M2D.createFace(bAlreadyExisted, pV5, pV8);
    Face *pF8 = M2D.createFace(bAlreadyExisted, pV8, pV4);
    (void) M2D.createTriCell(pF6, pF7, pF8, 1);

    Face *pF9 = M2D.createFace(bAlreadyExisted, pV8, pV7);
    Face *pF10 = M2D.createFace(bAlreadyExisted, pV7, pV4);
    (void) M2D.createTriCell(pF8, pF9, pF10, 1);

    Face *pF11 = M2D.createFace(bAlreadyExisted, pV7, pV6);
    Face *pF12 = M2D.createFace(bAlreadyExisted, pV6, pV4);
    (void) M2D.createTriCell(pF10, pF11, pF12, 1);

    Face *pF13 = M2D.createFace(bAlreadyExisted, pV6, pV3);
    Face *pF14 = M2D.createFace(bAlreadyExisted, pV3, pV4);
    (void) M2D.createTriCell(pF12, pF13, pF14, 1);

    Face *pF15 = M2D.createFace(bAlreadyExisted, pV3, pV0);
    (void) M2D.createTriCell(pF14, pF15, pF2, 1);

    M2D.createBFace(pF0);
    M2D.createBFace(pF3);
    M2D.createBFace(pF5);
    M2D.createBFace(pF7);
    M2D.createBFace(pF9);
    M2D.createBFace(pF11);
    M2D.createBFace(pF13);
    M2D.createBFace(pF15);
  }
  ~TMOPAnisoQual()
  {
  }
};

struct MetricAnisoSwap : public Mesh2DFixtureBase {
  MetricAnisoSwap() :
      Mesh2DFixtureBase()
  {

    M2D.disallowSwapRecursion();
    M2D.allowNonSimplicial();
    /*
     *Sketch: Mesh with 18 cells.  ---->
     *
     *	v12 --- f11 --- v13 --- f30---- v14 --- f32---- v15
     *	|            / //|  \ \ __   c16   \     c17     |
     *	|  c5     /   // |    \    \ ___      \          |
     *   f12    f9   / /  f5     f27     f28___  f29      f31
     *   |    /     / /   |        \            \ __ \    |
     *	| /       / / c2 |  c14      \      c15     \ \  |
     *   v8      f7 /    v9 ---- f26---- v10 --- f24---- v11
     *   |  c4  /  /    / |            /            ___ //|
     *   |     /  f6   /  | c13     /    c12 ____ /    /  |
     *   f10  /   /   /   f3     f25  ____f23    f22     f21
     *   |  /    /  /     |   / ___ /           /         |
     *	| /    /  /   c1 |/___/     c11    /        c10  |
     *   v4    /  f4     v5 ---- f18---- v6 ---- f20---- v7
     *   | c3 / /      /  |  \ \ __   c8    \     c9      |
     *   |   / /     /    |    \    \ ___      \          |
     *   f8 / /   f2      f1     f14     f16___ f17       f19
     *   | //   /         |        \            \ __ \    |
     *	|// /       c0   |  c6       \      c7      \ \  |
     *   v0 ---- f0 ---- v1 ---- f13---- v2 ---- f15---- v3
     */

    BOOST_CHECK(M2D.isValid());

    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(2, 0);
    Vert *pV3 = M2D.createVert(3, 0);
    Vert *pV4 = M2D.createVert(0, 0.1);
    Vert *pV5 = M2D.createVert(1, 0.1);
    Vert *pV6 = M2D.createVert(2, 0.1);
    Vert *pV7 = M2D.createVert(3, 0.1);
    Vert *pV8 = M2D.createVert(0, 0.2);
    Vert *pV9 = M2D.createVert(1, 0.2);
    Vert *pV10 = M2D.createVert(2, 0.2);
    Vert *pV11 = M2D.createVert(3, 0.2);
    Vert *pV12 = M2D.createVert(0, 0.3);
    Vert *pV13 = M2D.createVert(1, 0.3);
    Vert *pV14 = M2D.createVert(2, 0.3);
    Vert *pV15 = M2D.createVert(3, 0.3);

    (void) M2D.createTriCell(pV0, pV1, pV5, 1);

    Face *pF0 = findCommonFace(pV0, pV1);
    Face *pF1 = findCommonFace(pV1, pV5);
    Face *pF2 = findCommonFace(pV5, pV0);

    bool bAlreadyExisted = false;
    Face *pF3 = M2D.createFace(bAlreadyExisted, pV5, pV9);
    Face *pF4 = M2D.createFace(bAlreadyExisted, pV9, pV0);
    (void) M2D.createTriCell(pF2, pF3, pF4, 1);

    Face *pF5 = M2D.createFace(bAlreadyExisted, pV9, pV13);
    Face *pF6 = M2D.createFace(bAlreadyExisted, pV13, pV0);
    (void) M2D.createTriCell(pF6, pF4, pF5, 1);

    Face *pF7 = M2D.createFace(bAlreadyExisted, pV13, pV4);
    Face *pF8 = M2D.createFace(bAlreadyExisted, pV4, pV0);
    (void) M2D.createTriCell(pF6, pF7, pF8, 1);

    Face *pF9 = M2D.createFace(bAlreadyExisted, pV13, pV8);
    Face *pF10 = M2D.createFace(bAlreadyExisted, pV8, pV4);
    (void) M2D.createTriCell(pF7, pF9, pF10, 1);

    Face *pF11 = M2D.createFace(bAlreadyExisted, pV13, pV12);
    Face *pF12 = M2D.createFace(bAlreadyExisted, pV12, pV8);
    (void) M2D.createTriCell(pF9, pF11, pF12, 1);

    Face *pF13 = M2D.createFace(bAlreadyExisted, pV1, pV2);
    Face *pF14 = M2D.createFace(bAlreadyExisted, pV2, pV5);
    (void) M2D.createTriCell(pF13, pF14, pF1, 1);

    Face *pF15 = M2D.createFace(bAlreadyExisted, pV2, pV3);
    Face *pF16 = M2D.createFace(bAlreadyExisted, pV3, pV5);
    (void) M2D.createTriCell(pF14, pF15, pF16, 1);

    Face *pF17 = M2D.createFace(bAlreadyExisted, pV3, pV6);
    Face *pF18 = M2D.createFace(bAlreadyExisted, pV6, pV5);
    (void) M2D.createTriCell(pF17, pF18, pF16, 1);

    Face *pF19 = M2D.createFace(bAlreadyExisted, pV3, pV7);
    Face *pF20 = M2D.createFace(bAlreadyExisted, pV7, pV6);
    (void) M2D.createTriCell(pF17, pF19, pF20, 1);

    Face *pF21 = M2D.createFace(bAlreadyExisted, pV7, pV11);
    Face *pF22 = M2D.createFace(bAlreadyExisted, pV11, pV6);
    (void) M2D.createTriCell(pF21, pF22, pF20, 1);

    Face *pF23 = M2D.createFace(bAlreadyExisted, pV11, pV5);
    (void) M2D.createTriCell(pF23, pF18, pF22, 1);

    Face *pF24 = M2D.createFace(bAlreadyExisted, pV11, pV10);
    Face *pF25 = M2D.createFace(bAlreadyExisted, pV10, pV5);
    (void) M2D.createTriCell(pF25, pF23, pF24, 1);

    Face *pF26 = M2D.createFace(bAlreadyExisted, pV10, pV9);
    (void) M2D.createTriCell(pF3, pF25, pF26, 1);

    Face *pF27 = M2D.createFace(bAlreadyExisted, pV10, pV13);
    (void) M2D.createTriCell(pF26, pF27, pF5, 1);

    Face *pF28 = M2D.createFace(bAlreadyExisted, pV11, pV13);
    (void) M2D.createTriCell(pF28, pF27, pF24, 1);

    Face *pF29 = M2D.createFace(bAlreadyExisted, pV11, pV14);
    Face *pF30 = M2D.createFace(bAlreadyExisted, pV14, pV13);
    (void) M2D.createTriCell(pF29, pF30, pF28, 1);

    Face *pF31 = M2D.createFace(bAlreadyExisted, pV11, pV15);
    Face *pF32 = M2D.createFace(bAlreadyExisted, pV15, pV14);
    (void) M2D.createTriCell(pF31, pF32, pF29, 1);

    M2D.createBFace(pF0);
    M2D.createBFace(pF13);
    M2D.createBFace(pF15);
    M2D.createBFace(pF19);
    M2D.createBFace(pF21);
    M2D.createBFace(pF31);
    M2D.createBFace(pF32);
    M2D.createBFace(pF30);
    M2D.createBFace(pF11);
    M2D.createBFace(pF12);
    M2D.createBFace(pF10);
    M2D.createBFace(pF8);
  }
  ~MetricAnisoSwap()
  {
  }
};

struct TMOPAnisoCoarsening : public Mesh2DFixtureBase {
  TMOPAnisoCoarsening() :
      Mesh2DFixtureBase()
  {

    M2D.disallowSwapRecursion();
    M2D.allowNonSimplicial();

    /*
     *
     *Sketch: Mesh with 18 cells.  ---->
     *
     *	v12 --- f26 --- v13 --- f28 ---- v14 --- f32---- v15
     *   |  \             |            /  |  c16        / |
     *   |     \    c13   |   c14    /    |          /    |
     *   f23     f24     f25     f27     f29     f30     f31
     *   |   c12     \    |     /     c15 |     /         |
     *	|             \  |  /            |  /       c17  |
     *   v8  --- f15 --- v9 ---- f18 ---- v10 --- f22---- v11
     *   |              / |            /  |  \       c11  |
     *   |    c6     /    |   c8     /    |     \         |
     *   f13     f14     f16     f17     f19     f20     f21
     *   |     /    c7    |     /         |         \     |
     *	|  /             |  /       c9   |   c10      \  |
     *   v4 ---- f1 ---- v5 ---- f8 ---- v6 ---- f10---- v7
     *   |            /   |  \           |  \       c5    |
     *   |  c0      /     |     \   c3   |     \          |
     *   f2      f0       f4      f6    f7       f9      f12
     *   |     /     c1   |          \   |           \    |
     *	|  /             |  c2        \ |   c4        \  |
     *   v0 ---- f3 ---- v1 ---- f5 ---- v2 ---- f11---- v3
     */

    BOOST_CHECK(M2D.isValid());

    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(2, 0);
    Vert *pV3 = M2D.createVert(3, 0);
    Vert *pV4 = M2D.createVert(0, 0.04);
    Vert *pV5 = M2D.createVert(1, 0.1);
    Vert *pV6 = M2D.createVert(2.6, 0.03);
    Vert *pV7 = M2D.createVert(3, 0.1);
    Vert *pV8 = M2D.createVert(0, 0.1);
    Vert *pV9 = M2D.createVert(0.5, 0.15);
    Vert *pV10 = M2D.createVert(2, 0.1);
    Vert *pV11 = M2D.createVert(2.8, 0.16);
    Vert *pV12 = M2D.createVert(0, 0.2);
    Vert *pV13 = M2D.createVert(1, 0.2);
    Vert *pV14 = M2D.createVert(2, 0.2);
    Vert *pV15 = M2D.createVert(3, 0.2);

    (void) M2D.createTriCell(pV0, pV5, pV4, 1);

    Face *pF0 = findCommonFace(pV0, pV5);
    Face *pF1 = findCommonFace(pV5, pV4);
    Face *pF2 = findCommonFace(pV4, pV0);

    bool bAlreadyExisted = false;
    Face *pF3 = M2D.createFace(bAlreadyExisted, pV0, pV1);
    Face *pF4 = M2D.createFace(bAlreadyExisted, pV1, pV5);
    (void) M2D.createTriCell(pF0, pF3, pF4, 1);

    Face *pF5 = M2D.createFace(bAlreadyExisted, pV1, pV2);
    Face *pF6 = M2D.createFace(bAlreadyExisted, pV2, pV5);
    (void) M2D.createTriCell(pF4, pF5, pF6, 1);

    Face *pF7 = M2D.createFace(bAlreadyExisted, pV2, pV6);
    Face *pF8 = M2D.createFace(bAlreadyExisted, pV6, pV5);
    (void) M2D.createTriCell(pF6, pF7, pF8, 1);

    Face *pF9 = M2D.createFace(bAlreadyExisted, pV3, pV6);
    Face *pF10 = M2D.createFace(bAlreadyExisted, pV7, pV6);
    Face *pF11 = M2D.createFace(bAlreadyExisted, pV2, pV3);
    (void) M2D.createTriCell(pF7, pF11, pF9, 1);

    Face *pF12 = M2D.createFace(bAlreadyExisted, pV3, pV7);
    (void) M2D.createTriCell(pF9, pF12, pF10, 1);

    Face *pF13 = M2D.createFace(bAlreadyExisted, pV8, pV4);
    Face *pF14 = M2D.createFace(bAlreadyExisted, pV4, pV9);
    Face *pF15 = M2D.createFace(bAlreadyExisted, pV9, pV8);
    (void) M2D.createTriCell(pF13, pF14, pF15, 1);

    Face *pF16 = M2D.createFace(bAlreadyExisted, pV5, pV9);
    (void) M2D.createTriCell(pF1, pF16, pF14, 1);

    Face *pF17 = M2D.createFace(bAlreadyExisted, pV5, pV10);
    Face *pF18 = M2D.createFace(bAlreadyExisted, pV10, pV9);
    (void) M2D.createTriCell(pF16, pF17, pF18, 1);

    Face *pF19 = M2D.createFace(bAlreadyExisted, pV6, pV10);
    (void) M2D.createTriCell(pF8, pF19, pF17, 1);

    Face *pF20 = M2D.createFace(bAlreadyExisted, pV7, pV10);
    (void) M2D.createTriCell(pF10, pF20, pF19, 1);

    Face *pF21 = M2D.createFace(bAlreadyExisted, pV7, pV11);
    Face *pF22 = M2D.createFace(bAlreadyExisted, pV11, pV10);
    (void) M2D.createTriCell(pF20, pF21, pF22, 1);

    Face *pF23 = M2D.createFace(bAlreadyExisted, pV12, pV8);
    Face *pF24 = M2D.createFace(bAlreadyExisted, pV9, pV12);
    (void) M2D.createTriCell(pF23, pF15, pF24, 1);

    Face *pF25 = M2D.createFace(bAlreadyExisted, pV9, pV13);
    Face *pF26 = M2D.createFace(bAlreadyExisted, pV13, pV12);
    (void) M2D.createTriCell(pF24, pF25, pF26, 1);

    Face *pF27 = M2D.createFace(bAlreadyExisted, pV9, pV14);
    Face *pF28 = M2D.createFace(bAlreadyExisted, pV14, pV13);
    (void) M2D.createTriCell(pF25, pF27, pF28, 1);

    Face *pF29 = M2D.createFace(bAlreadyExisted, pV10, pV14);
    (void) M2D.createTriCell(pF18, pF29, pF27, 1);

    Face *pF30 = M2D.createFace(bAlreadyExisted, pV10, pV15);
    Face *pF31 = M2D.createFace(bAlreadyExisted, pV11, pV15);
    (void) M2D.createTriCell(pF22, pF31, pF30, 1);

    Face *pF32 = M2D.createFace(bAlreadyExisted, pV15, pV14);
    (void) M2D.createTriCell(pF30, pF32, pF29, 1);

    M2D.createBFace(pF3);
    M2D.createBFace(pF5);
    M2D.createBFace(pF11);
    M2D.createBFace(pF12);
    M2D.createBFace(pF21);
    M2D.createBFace(pF31);
    M2D.createBFace(pF32);
    M2D.createBFace(pF28);
    M2D.createBFace(pF26);
    M2D.createBFace(pF23);
    M2D.createBFace(pF13);
    M2D.createBFace(pF2);
  }
  ~TMOPAnisoCoarsening()
  {
  }
};

struct TMOPAnisoEdgeInsertion : public Mesh2DFixtureBase {
  TMOPAnisoEdgeInsertion() :
      Mesh2DFixtureBase()
  {

    M2D.disallowSwapRecursion();
    M2D.allowNonSimplicial();

    /*
     *Sketch: Mesh with 6 cells.
     *
     *  v5 ---f9--- v6 ---f7--- v7
     *  |  c3      / | \    c2   |
     *  |        /   |   \       |
     *  f10   f11    |    f8    f6
     *  |    /       |       \   |
     *  |  /         |         \ |
     *  v3    c4    f12    c5   v4
     *  |  \         |         / |
     *  |    \       |       /   |
     *  f2    f1     |    f5    f4
     *  |        \   |   /       |
     *  |  c0      \ | /     c1  |
     *  v0 ---f0--- v1 ---f3--- v2
     */

    BOOST_CHECK(M2D.isValid());

    Vert *pV0 = M2D.createVert(0, 0);
    Vert *pV1 = M2D.createVert(1, 0);
    Vert *pV2 = M2D.createVert(2, 0);
    Vert *pV3 = M2D.createVert(0, 1);
    Vert *pV4 = M2D.createVert(2, 1);
    Vert *pV5 = M2D.createVert(0, 2);
    Vert *pV6 = M2D.createVert(1, 2);
    Vert *pV7 = M2D.createVert(2, 2);

    pV0->setType(Vert::eBdryApex);
    pV1->setType(Vert::eBdry);
    pV2->setType(Vert::eBdryApex);
    pV3->setType(Vert::eBdry);
    pV4->setType(Vert::eBdry);
    pV5->setType(Vert::eBdryApex);
    pV6->setType(Vert::eBdry);
    pV7->setType(Vert::eBdryApex);

    (void) M2D.createTriCell(pV0, pV1, pV3, 1);

    Face *pF0 = findCommonFace(pV0, pV1);
    Face *pF1 = findCommonFace(pV1, pV3);
    Face *pF2 = findCommonFace(pV3, pV0);

    bool bAlreadyExisted = false;
    Face *pF3 = M2D.createFace(bAlreadyExisted, pV1, pV2);
    Face *pF4 = M2D.createFace(bAlreadyExisted, pV2, pV4);
    Face *pF5 = M2D.createFace(bAlreadyExisted, pV4, pV1);
    (void) M2D.createTriCell(pF3, pF4, pF5, 1);

    Face *pF6 = M2D.createFace(bAlreadyExisted, pV4, pV7);
    Face *pF7 = M2D.createFace(bAlreadyExisted, pV7, pV6);
    Face *pF8 = M2D.createFace(bAlreadyExisted, pV6, pV4);
    (void) M2D.createTriCell(pF6, pF7, pF8, 1);

    Face *pF9 = M2D.createFace(bAlreadyExisted, pV6, pV5);
    Face *pF10 = M2D.createFace(bAlreadyExisted, pV5, pV3);
    Face *pF11 = M2D.createFace(bAlreadyExisted, pV3, pV6);
    (void) M2D.createTriCell(pF9, pF10, pF11, 1);

    Face *pF12 = M2D.createFace(bAlreadyExisted, pV1, pV6);
    (void) M2D.createTriCell(pF1, pF12, pF11, 1);
    (void) M2D.createTriCell(pF5, pF8, pF12, 1);

    M2D.createBFace(pF0, 1);
    M2D.createBFace(pF3, 1);
    M2D.createBFace(pF4, 1);
    M2D.createBFace(pF6, 1);
    M2D.createBFace(pF7, 1);
    M2D.createBFace(pF9, 1);
    M2D.createBFace(pF10, 1);
    M2D.createBFace(pF2, 1);

  }
  ~TMOPAnisoEdgeInsertion()
  {
  }
};

BOOST_FIXTURE_TEST_SUITE(TestTMOPAnisoQual, TMOPAnisoQual)

  BOOST_AUTO_TEST_CASE(Setup_AnisoQual)
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(9, 16, 8, 8, 8, 0);
  }

  BOOST_AUTO_TEST_CASE(AnisoQualityShape)
  {

    GRUMMP::TMOPQual TQ(2);
    GRUMMP::TMOPQual::eQMetric QM;
    double dMetric[3] =
      { 1, 0, 100 };

    QM = GRUMMP::TMOPQual::eShape;
    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {

      Cell *pC = M2D.getCell(i);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);

      double dT[9];
      double dQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      dQ = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		    pVB->getCoords(), pVC->getCoords());

      BOOST_CHECK(dQ < 1e-15);
    }
  }

  BOOST_AUTO_TEST_CASE(AnisoQualitySize)
  {

    GRUMMP::TMOPQual TQ(2);
    GRUMMP::TMOPQual::eQMetric QM;
    double dMetric[3] =
      { 1, 0, 100 };

    QM = GRUMMP::TMOPQual::eSize;
    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {

      Cell *pC = M2D.getCell(i);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);

      double dT[9];
      double dQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      dQ = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		    pVB->getCoords(), pVC->getCoords());

      BOOST_CHECK(dQ < 1e-15);
    }
  }

  BOOST_AUTO_TEST_CASE(AnisoQualityShapeSize)
  {

    GRUMMP::TMOPQual TQ(2);
    GRUMMP::TMOPQual::eQMetric QM;
    double dMetric[3] =
      { 1, 0, 100 };

    QM = GRUMMP::TMOPQual::eShapeSize;
    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {

      Cell *pC = M2D.getCell(i);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);

      double dT[9];
      double dQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      dQ = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		    pVB->getCoords(), pVC->getCoords());

      BOOST_CHECK(dQ < 1e-15);
    }
  }

  BOOST_AUTO_TEST_CASE(AnisoQualityShapeOrient)
  {

    GRUMMP::TMOPQual TQ(2);
    GRUMMP::TMOPQual::eQMetric QM;
    double dMetric[3] =
      { 1, 0, 100 };

    QM = GRUMMP::TMOPQual::eShapeOrient;
    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {

      Cell *pC = M2D.getCell(i);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);

      double dT[9];
      double dQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      dQ = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		    pVB->getCoords(), pVC->getCoords());

      BOOST_CHECK(dQ < 1e-15);
    }
  }

  BOOST_AUTO_TEST_CASE(AnisoQualityShapeSizeOrient)
  {

    GRUMMP::TMOPQual TQ(2);
    GRUMMP::TMOPQual::eQMetric QM;
    double dMetric[3] =
      { 1, 0, 100 };

    QM = GRUMMP::TMOPQual::eShapeSizeOrient;
    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {

      Cell *pC = M2D.getCell(i);
      Vert *pVA = pC->getVert(0);
      Vert *pVB = pC->getVert(1);
      Vert *pVC = pC->getVert(2);

      double dT[9];
      double dQ;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      dQ = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		    pVB->getCoords(), pVC->getCoords());

      BOOST_CHECK(dQ < 1e-15);
    }
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(TestAnisoSwap, MetricAnisoSwap)

  BOOST_AUTO_TEST_CASE(Setup_AnisoSwap)
  {

    BOOST_CHECK(M2D.isValid());
    checkMeshSize(16, 33, 12, 18, 18, 0);

    double dMetric[3] =
      { 1, 0, 100 };

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    GRUMMP::QualMeasure *pQM = new GRUMMP::Angle2D;
    GRUMMP::SwapDecider2D* pSD2D = new GRUMMP::MetricSwapDecider2D(pQM);
    GRUMMP::SwapManager2D *pSM2D = new GRUMMP::SwapManager2D(pSD2D, &M2D);

    pSM2D->swapAllFaces();
    M2D.purgeAllEntities();

    BOOST_CHECK(M2D.isValid());
    checkMeshSize(16, 33, 12, 18, 18, 0);

    Vert *apV[16];
    Cell *apC[18];

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      apV[i] = M2D.getVert(i);
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {
      apC[i] = M2D.getCell(i);
    }

    checkCanonicalTri(apC[0], apV[0], apV[1], apV[5]);
    checkCanonicalTri(apC[1], apV[0], apV[5], apV[4]);
    checkCanonicalTri(apC[2], apV[4], apV[5], apV[9]);
    checkCanonicalTri(apC[3], apV[4], apV[9], apV[8]);
    checkCanonicalTri(apC[4], apV[8], apV[9], apV[13]);
    checkCanonicalTri(apC[5], apV[8], apV[13], apV[12]);
    checkCanonicalTri(apC[6], apV[1], apV[2], apV[5]);
    checkCanonicalTri(apC[7], apV[5], apV[2], apV[6]);
    checkCanonicalTri(apC[8], apV[2], apV[3], apV[6]);
    checkCanonicalTri(apC[9], apV[6], apV[3], apV[7]);
    checkCanonicalTri(apC[10], apV[6], apV[7], apV[11]);
    checkCanonicalTri(apC[11], apV[5], apV[6], apV[10]);
    checkCanonicalTri(apC[12], apV[6], apV[11], apV[10]);
    checkCanonicalTri(apC[13], apV[5], apV[10], apV[9]);
    checkCanonicalTri(apC[14], apV[9], apV[10], apV[13]);
    checkCanonicalTri(apC[15], apV[13], apV[10], apV[14]);
    checkCanonicalTri(apC[16], apV[10], apV[11], apV[14]);
    checkCanonicalTri(apC[17], apV[14], apV[11], apV[15]);

    delete pSM2D;
    delete pSD2D;
    delete pQM;
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(TestTMOPTMOPAnisoCoarseningQuality, TMOPAnisoCoarsening)

  BOOST_AUTO_TEST_CASE(Setup_AnisoCoarsenQuality)
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(16, 33, 12, 18, 18, 0);
  }

  BOOST_AUTO_TEST_CASE(AnisoCoarsenQuality)
  {

    double dMetric[3] =
      { 1, 0, 100 };
    double dMaxMLength = sqrt(2);
    double dMinMLength = sqrt(2) / 2.;

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    AnisoRefinement *pAR = new AnisoRefinement(&M2D, 0);
    GRUMMP::QualMeasure *pQM = new GRUMMP::Angle2D;
    GRUMMP::SwapDecider2D* pSD2D = new GRUMMP::MetricSwapDecider2D(pQM);
    GRUMMP::InsertionPointCalculator * pIPC = new GRUMMP::AnisoEdgeIPC2D(
	dMinMLength, dMaxMLength);
    GRUMMP::SmoothingManagerMesquite2D *pSMM2D =
	new GRUMMP::SmoothingManagerMesquite2D(&M2D);
    GRUMMP::AnisoEdgeRefinementManager2D ARM2D(&M2D, dMaxMLength, dMinMLength,
					       false, pIPC, pSMM2D, pAR, pSD2D);

    M2D.getVert(0)->setType(Vert::eBdryApex);
    M2D.getVert(3)->setType(Vert::eBdryApex);
    M2D.getVert(12)->setType(Vert::eBdryApex);
    M2D.getVert(15)->setType(Vert::eBdryApex);
    M2D.getVert(9)->setType(Vert::eInterior);

    ARM2D.AnisoMetricEdgeCoarsening();
    ARM2D.AnisoMetricEdgeCoarsening();
    ARM2D.AnisoMetricEdgeCoarsening();

    M2D.purgeAllEntities();
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(12, 23, 10, 12, 12, 0);

    Vert *apV[12];
    Cell *apC[12];

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      apV[i] = M2D.getVert(i);
      logMessage(MSG_DEBUG, "Vert %d: %p (%f, %f)\n", i, apV[i], apV[i]->x(),
		 apV[i]->y());
    }

    for (unsigned i = 0; i < M2D.getNumCells(); i++) {
      apC[i] = M2D.getCell(i);
    }

    logMessage(MSG_DEBUG, "Checking tri 0\n");
    checkCanonicalTri(apC[0], apV[7], apV[11], apV[10]);
    logMessage(MSG_DEBUG, "Checking tri 1\n");
    checkCanonicalTri(apC[1], apV[0], apV[1], apV[4]);
    logMessage(MSG_DEBUG, "Checking tri 2\n");
    checkCanonicalTri(apC[2], apV[4], apV[1], apV[2]);
    logMessage(MSG_DEBUG, "Checking tri 3\n");
    checkCanonicalTri(apC[3], apV[7], apV[5], apV[11]);
    logMessage(MSG_DEBUG, "Checking tri 4\n");
    checkCanonicalTri(apC[4], apV[0], apV[4], apV[6]);
    logMessage(MSG_DEBUG, "Checking tri 5\n");
    checkCanonicalTri(apC[5], apV[2], apV[5], apV[7]);
    logMessage(MSG_DEBUG, "Checking tri 6\n");
    checkCanonicalTri(apC[6], apV[6], apV[4], apV[8]);
    logMessage(MSG_DEBUG, "Checking tri 7\n");
    checkCanonicalTri(apC[7], apV[4], apV[9], apV[8]);
    logMessage(MSG_DEBUG, "Checking tri 8\n");
    checkCanonicalTri(apC[8], apV[4], apV[7], apV[10]);
    logMessage(MSG_DEBUG, "Checking tri 9\n");
    checkCanonicalTri(apC[9], apV[4], apV[2], apV[7]);
    logMessage(MSG_DEBUG, "Checking tri 10\n");
    checkCanonicalTri(apC[10], apV[2], apV[3], apV[5]);
    logMessage(MSG_DEBUG, "Checking tri 11\n");
    checkCanonicalTri(apC[11], apV[4], apV[10], apV[9]);
    logMessage(MSG_DEBUG, "Done\n");

    delete pSMM2D;
    delete pIPC;
    delete pSD2D;
    delete pQM;
    delete pAR;
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(TestTMOPAnisoEdgeInsertion, TMOPAnisoEdgeInsertion)

  BOOST_AUTO_TEST_CASE(Setup_AnisoEdgeInsertion)
  {
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(8, 13, 8, 6, 6, 0);

  }

  BOOST_AUTO_TEST_CASE(AnisoEdgeInsertion)
  {
    BOOST_CHECK(M2D.isValid());

    double dMetric[3] =
      { 1, 0, 100 };
    double dMaxMLength = sqrt(2);
    double dMinMLength = sqrt(2) / 2.;

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    AnisoRefinement *pAR = new AnisoRefinement(&M2D, 0);
    GRUMMP::QualMeasure *pQM = new GRUMMP::Angle2D;
    GRUMMP::SwapDecider2D* pSD2D = new GRUMMP::MetricSwapDecider2D(pQM);
    GRUMMP::InsertionPointCalculator * pIPC = new GRUMMP::AnisoEdgeIPC2D(
	dMinMLength, dMaxMLength);
    GRUMMP::SmoothingManagerMesquite2D *pSMM2D =
	new GRUMMP::SmoothingManagerMesquite2D(&M2D);
    GRUMMP::AnisoEdgeRefinementManager2D ARM2D(&M2D, dMaxMLength, dMinMLength,
					       false, pIPC, pSMM2D, pAR, pSD2D);

    ARM2D.AnisoQualityEdgeInsertion();
    ARM2D.AnisoQualityEdgeInsertion();

    M2D.purgeAllEntities();
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(51, 114, 36, 64, 64, 0);

    dMetric[0] = 4;

    for (unsigned i = 0; i < M2D.getNumVerts(); i++) {
      M2D.getVert(i)->setMetric(dMetric);
    }

    ARM2D.AnisoMetricLengthEdgeInsertion();

    M2D.purgeAllEntities();
    BOOST_CHECK(M2D.isValid());
    checkMeshSize(113, 296, 40, 184, 184, 0);

    delete pSMM2D;
    delete pIPC;
    delete pSD2D;
    delete pQM;
    delete pAR;

  }

  BOOST_AUTO_TEST_SUITE_END()

#endif

BOOST_FIXTURE_TEST_SUITE(ADT_Tests, MeshMod2DFixture)
  ;

  BOOST_AUTO_TEST_CASE(Build_ADT)
  {
    ADT tree(&M2D);

    BOOST_CHECK_EQUAL(tree.iNumVars(), 2);
    BOOST_CHECK_EQUAL(tree.iNumDimensions(), 2);
    BOOST_CHECK_EQUAL(tree.iTreeSize(), M2D.getNumVerts());
  }

  BOOST_AUTO_TEST_CASE(QueryADT)
  {
    ADT tree(&M2D);
    double range[] =
      { -0.5, 1.5, -0.5, 1.5 };
    std::vector<GR_index_t> result(tree.veciRangeQuery(range));
    BOOST_CHECK_EQUAL(result.size(), 4);
    logMessage(0, "%d %d %d %d\n", result[0], result[1], result[2], result[3]);
  }

  BOOST_AUTO_TEST_CASE(MatchVert)
  {
    ADT tree(&M2D);
    Vert* target = M2D.getVert(0);
    Vert* cand = findVertex(&M2D, tree, target->getCoords());
    BOOST_CHECK_EQUAL(target, cand);
  }
  BOOST_AUTO_TEST_SUITE_END()

