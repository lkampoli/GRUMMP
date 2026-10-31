#define BOOST_TEST_MODULE VolMesh_Tests
#include "boost/test/unit_test.hpp"

#include "VolMeshFixtureBase.h"

#include "GR_VolMesh.h"

BOOST_AUTO_TEST_SUITE(VolMesh_ReadTests)

  BOOST_AUTO_TEST_CASE(Read_VMesh)
  {
    VolMesh VM;
    VM.readFromFile("meshes/sphere.vmesh");
    BOOST_CHECK(VM.isValid());

    BOOST_CHECK_MESSAGE(VM.getNumVerts() == 768,
			"VM.iNumVerts() result: " << VM.getNumVerts());
    BOOST_CHECK_MESSAGE(VM.getNumFaces() == 6163,
			"VM.iNumFaces() result: " << VM.getNumFaces());
    BOOST_CHECK_MESSAGE(VM.getNumCells() == 2821,
			"VM.iNumCells() result: " << VM.getNumCells());
    BOOST_CHECK_MESSAGE(VM.getNumBdryFaces() == 1042,
			"VM.iNumBdryFaces() result: " << VM.getNumBdryFaces());

    // Grab cells and confirm their connectivity is reasonable.
    GR_index_t nCells = VM.getNumCells();

    GR_index_t iC = 0;

    for (iC = 0; iC < nCells; iC++) {
      Cell* pC = VM.getCell(iC);
      if (pC->isDeleted())
	continue;
      BOOST_CHECK(pC->isClosed());
    }
  }

  BOOST_AUTO_TEST_CASE(Read_VTK)
  {
    VolMesh VM;
    VM.readFromFile("meshes/sphere.vtk");
    BOOST_CHECK(VM.isValid());

    BOOST_CHECK_MESSAGE(VM.getNumVerts() == 768,
			"VM.iNumVerts() result: " << VM.getNumVerts());
    BOOST_CHECK_MESSAGE(VM.getNumFaces() == 6163,
			"VM.iNumFaces() result: " << VM.getNumFaces());
    BOOST_CHECK_MESSAGE(VM.getNumCells() == 2821,
			"VM.iNumCells() result: " << VM.getNumCells());
    BOOST_CHECK_MESSAGE(VM.getNumBdryFaces() == 1042,
			"VM.iNumBdryFaces() result: " << VM.getNumBdryFaces());
    // Grab cells and confirm their connectivity is reasonable.
    GR_index_t nCells = VM.getNumCells();

    GR_index_t iC = 0;

    for (iC = 0; iC < nCells; iC++) {
      Cell* pC = VM.getCell(iC);
      if (pC->isDeleted())
	continue;
      BOOST_CHECK(pC->isClosed());
    }
  }

  BOOST_AUTO_TEST_CASE(Read_Mixed_Element_VTK)
  {
    VolMesh VM;
    VM.readFromFile("meshes/globe_mats.vtk");
    BOOST_CHECK(VM.isValid());

    BOOST_CHECK_MESSAGE(VM.getNumVerts() == 1093,
			"VM.iNumVerts() result: " << VM.getNumVerts());
    BOOST_CHECK_MESSAGE(VM.getNumFaces() == 3480,
			"VM.iNumFaces() result: " << VM.getNumFaces());
    BOOST_CHECK_MESSAGE(VM.getNumCells() == 1200,
			"VM.iNumCells() result: " << VM.getNumCells());
    BOOST_CHECK_MESSAGE(VM.getNumBdryFaces() == 200,
			"VM.iNumBdryFaces() result: " << VM.getNumBdryFaces());

    // Grab cells and confirm their connectivity is reasonable.
    GR_index_t nCells = VM.getNumCells();

    GR_index_t iC = 0;

    for (iC = 0; iC < nCells; iC++) {
      Cell* pC = VM.getCell(iC);
      if (pC->isDeleted())
	continue;
      BOOST_CHECK(pC->isClosed());
    }
  }

  BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(VolMesh_Constructors)

  BOOST_AUTO_TEST_CASE(FromCubitBox)
  {
    CubitBox box(
      { 0, 0, 0 },
		   { 1, 1, 1 });
    VolMesh VM(box, 1, 3.5);

    BOOST_CHECK(VM.isWatertight());
    BOOST_CHECK_CLOSE(VM.calcInteriorSize(), 8, 1.e-8);
    BOOST_CHECK_CLOSE(VM.calcBoundarySize(), 24, 1.e-8);
    BOOST_CHECK_EQUAL(VM.getVertIndex(VM.getVert(3)), 3);
    BOOST_CHECK_EQUAL(VM.getFaceIndex(VM.getFace(3)), 3);
    BOOST_CHECK_EQUAL(VM.getCellIndex(VM.getCell(3)), 3);
    BOOST_CHECK_EQUAL(VM.getBFaceIndex(VM.getBFace(3)), 3);
  }

  BOOST_AUTO_TEST_SUITE_END()
