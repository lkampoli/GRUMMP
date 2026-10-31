#include <vector>
#define BOOST_TEST_MODULE GRUMMP_MPI
#include <boost/test/unit_test.hpp>

#include "GR_Mesh2D.h"

BOOST_AUTO_TEST_CASE(Dummy)
{
  BOOST_CHECK(true);
}

#ifdef HAVE_MPI
#include <mpi.h>
#include "GR_Ghost.h"
#include "GR_Partition.h"

using namespace GRUMMP;

struct GRUMMPTestConfig {
  GRUMMPTestConfig()
  {
    int flag = MPI::Is_initialized();
    if (!flag) {
      MPI::Init(boost::unit_test::framework::master_test_suite().argc,
		boost::unit_test::framework::master_test_suite().argv);
      GRUMMPInit("MPI tests");
    }
  }
  ~GRUMMPTestConfig()
  {
    // GRUMMP signoff gets called automagically.
    MPI::Finalize();
  }
};

BOOST_GLOBAL_FIXTURE(GRUMMPTestConfig);

BOOST_AUTO_TEST_SUITE(MPI_GRUMMP_Tests)

  BOOST_AUTO_TEST_CASE(Basic)
  {
    int isInit = MPI::Is_initialized();
    BOOST_REQUIRE(isInit);

    BOOST_CHECK_EQUAL(PART_BITS, 12);
    BOOST_CHECK_EQUAL(PROC_BITS, 20);

    BOOST_CHECK_EQUAL(PART_UNDEFINED, 0);
//    BOOST_CHECK_EQUAL(ENTITY_UNDEFINED, 2147483647);

    int numProcs = MPI::COMM_WORLD.Get_size();
    {
      GRUMMP::Partition testPartn(Mesh::eMesh2D);
      BOOST_CHECK_EQUAL(testPartn.getNumLocalParts(), 1);
      BOOST_CHECK_EQUAL(testPartn.getNumGlobalParts(), numProcs);
    }
    {
      GRUMMP::Partition testPartn(Mesh::eMesh2D, 3);
      BOOST_CHECK_EQUAL(testPartn.getNumLocalParts(), 3);
      BOOST_CHECK_EQUAL(testPartn.getNumGlobalParts(), 3 * numProcs);

      Mesh* pM2D = testPartn.getPart(-1);
      BOOST_CHECK(pM2D == NULL);
      pM2D = testPartn.getPart(3);
      BOOST_CHECK(pM2D == NULL);

      pM2D = testPartn.getPart(0);
      BOOST_CHECK_EQUAL(pM2D->getType(), Mesh::eMesh2D);
      BOOST_CHECK_EQUAL(pM2D->getNumCells(), 0);
      BOOST_CHECK_EQUAL(pM2D->getNumVerts(), 0);

      // Only read into processor 0, so that we only have one copy to partition.
      int myRank = MPI::COMM_WORLD.Get_rank();
      if (myRank == 0) {
	pM2D->readFromFile("meshes/0012-spline.mesh");
	BOOST_CHECK_EQUAL(pM2D->getNumCells(), 173);
	BOOST_CHECK_EQUAL(pM2D->getNumVerts(), 98);
	BOOST_CHECK_EQUAL(pM2D->getNumFaces(), 271);
	BOOST_CHECK_EQUAL(pM2D->getNumBdryFaces(), 23);
	BOOST_CHECK_EQUAL(pM2D->getNumIntBdryFaces(), 0);
      }
      for (int ii = 0; ii < testPartn.getNumLocalParts(); ii++) {
	PartGID thisPGID = testPartn.getPartGID(ii);
	BOOST_CHECK_EQUAL(thisPGID.getRank(), myRank);
	BOOST_CHECK_EQUAL(thisPGID.getPartLocalID(), ii);
      }
    }
  }

  static void
  checkNonLocalEnts(const Mesh* const pM, const int expGhostedVerts,
		    const int expPartBdryVerts, const int expForeignVerts,
		    const int expGhostedCells, const int expPartBdryCells,
		    const int expForeignCells, const int expGhostedFaces,
		    const int expPartBdryFaces, const int expForeignFaces)
  {
    logMessage(3, "Counting non-interior entities.\n");
    unsigned int uniqueID = pM->getPartID().getUniqueID();
    // Count and check the number of ghosted and foreign verts and cells.
    int numForeignVerts = 0, numGhostedVerts = 0, numPartBdryVerts = 0;
    for (GR_index_t ii = 0; ii < pM->getNumVerts(); ++ii) {
      Vert* pV = pM->getVert(ii);
      if (pV->isDeleted())
	continue;
      switch (pV->getEntStatus())
	{
	case Entity::eGhosted:
	  ++numGhostedVerts;
	  break;
	case Entity::ePartBdry:
	  ++numPartBdryVerts;
	  break;
	case Entity::eForeign:
	  ++numForeignVerts;
	  logMessage(4, "Part %u: (%f %f %f)\n", uniqueID, pV->x(), pV->y(),
		     pV->z());
	  break;
	default:
	  break;
	}
    }
    BOOST_CHECK_EQUAL(numGhostedVerts, expGhostedVerts);
    BOOST_CHECK_EQUAL(numPartBdryVerts, expPartBdryVerts);
    BOOST_CHECK_EQUAL(numForeignVerts, expForeignVerts);

    int numForeignCells = 0, numGhostedCells = 0, numPartBdryCells = 0;
    for (GR_index_t ii = 0; ii < pM->getNumCells(); ++ii) {
      Cell* pC = pM->getCell(ii);
      if (pC->isDeleted())
	continue;
      switch (pC->getEntStatus())
	{
	case Entity::eGhosted:
	  ++numGhostedCells;
	  break;
	case Entity::ePartBdry:
	  ++numPartBdryCells;
	  break;
	case Entity::eForeign:
	  ++numForeignCells;
	  break;
	default:
	  break;
	}
    }

    BOOST_CHECK_EQUAL(numGhostedCells, expGhostedCells);
    BOOST_CHECK_EQUAL(numPartBdryCells, expPartBdryCells);
    BOOST_CHECK_EQUAL(numForeignCells, expForeignCells);

    int numForeignFaces = 0, numGhostedFaces = 0, numPartBdryFaces = 0;
    for (GR_index_t ii = 0; ii < pM->getNumFaces(); ++ii) {
      Face* pF = pM->getFace(ii);
      if (pF->isDeleted())
	continue;
      switch (pF->getEntStatus())
	{
	case Entity::eGhosted:
	  ++numGhostedFaces;
	  break;
	case Entity::ePartBdry:
	  ++numPartBdryFaces;
	  break;
	case Entity::eForeign:
	  ++numForeignFaces;
	  break;
	default:
	  break;
	}
    }
    BOOST_CHECK_EQUAL(numGhostedFaces, expGhostedFaces);
    BOOST_CHECK_EQUAL(numPartBdryFaces, expPartBdryFaces);
    BOOST_CHECK_EQUAL(numForeignFaces, expForeignFaces);
  }

  BOOST_AUTO_TEST_CASE(ReadParallel)
  {
    int isInit = MPI::Is_initialized();
    BOOST_REQUIRE(isInit);

    int numProcs = MPI::COMM_WORLD.Get_size();
    BOOST_CHECK(numProcs == 4); // Because of the number of files to read.

    GRUMMP::Partition testPartn(Mesh::eMesh2D);
    BOOST_CHECK_EQUAL(testPartn.getNumLocalParts(), 1);
    BOOST_CHECK_EQUAL(testPartn.getNumGlobalParts(), numProcs);

    int myRank = MPI::COMM_WORLD.Get_rank();
    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(testPartn.getPart(0));
    BOOST_CHECK(pM2D != nullptr);

    char filename[1024];
    snprintf(filename, 1024, "meshes/simple-parallel-%d.vtk", myRank);
    pM2D->readFromFile(filename);
    BOOST_CHECK_EQUAL(pM2D->getNumCells(), 5);
    BOOST_CHECK_EQUAL(pM2D->getNumVerts(), 9);
    BOOST_CHECK_EQUAL(pM2D->getNumFaces(), 13);
    BOOST_CHECK_EQUAL(pM2D->getNumBdryFaces(), 4);
    BOOST_CHECK_EQUAL(pM2D->getNumIntBdryFaces(), 0);
//    BOOST_CHECK_EQUAL(pM2D->getNumPartBdryFaces(), 4);

    checkNonLocalEnts(pM2D, 0, 5, 0, 0, 0, 0, 0, 4, 0);

//    for (GR_index_t ii = 0; ii < pM2D->getNumVerts(); ++ii) {
//        Vert *pV = pM2D->getVert(ii);
//        if (pV->isDeleted()) continue;
//        bool success = false;
//        Ghost g = pM2D->getGhostEntry(pV, success);
//        if (success) {
//            logMessage(4, "Part %u, vert %d (%f %f %f) has %d copies (%d orig)\n",
//                   pM2D->getPartID().getUniqueID(), ii, pV->x(), pV->y(), pV->z(),
//                   g.getNumCopies(), g.getNumPartBdryCopies());
//        }
//    }

    pM2D->matchPartBdrys();
    BOOST_CHECK(pM2D->isValid());
    pM2D->setupGhostData();
    BOOST_CHECK(pM2D->isValid());

    checkNonLocalEnts(pM2D, 3, 5, 7, (myRank == 3 ? 3 : 4), 0,
		      5 + myRank - (myRank == 2 ? 1 : 0), 0, 4,
		      12 + myRank - (myRank == 2 ? 1 : 0));

    pM2D->clearGhostData();
    BOOST_CHECK(pM2D->isValid());
    BOOST_CHECK_EQUAL(pM2D->getNumCells(), 5);
    BOOST_CHECK_EQUAL(pM2D->getNumVerts(), 9);
    BOOST_CHECK_EQUAL(pM2D->getNumFaces(), 13);
    BOOST_CHECK_EQUAL(pM2D->getNumBdryFaces(), 4);
    BOOST_CHECK_EQUAL(pM2D->getNumIntBdryFaces(), 0);

    checkNonLocalEnts(pM2D, 0, 5, 0, 0, 0, 0, 0, 4, 0);

    pM2D->setupGhostData();
    BOOST_CHECK(pM2D->isValid());

    checkNonLocalEnts(pM2D, 3, 5, 7, (myRank == 3 ? 3 : 4), 0,
		      5 + myRank - (myRank == 2 ? 1 : 0), 0, 4,
		      12 + myRank - (myRank == 2 ? 1 : 0));
  }

  BOOST_AUTO_TEST_CASE(Ghost)
  {
    int myRank = MPI::COMM_WORLD.Get_rank();
    PartGID myID(myRank);
    int size = MPI::COMM_WORLD.Get_size();

    ::GRUMMP::Ghost g(myID, 3);

    g.setOwnerID(myID);
    PartGID returnedID = g.getOwnerID();
    bool test = (myID == returnedID);
    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(g.getMaxNumCopies(), 3);
    BOOST_CHECK_EQUAL(g.getNumCopies(), 0);

    int otherRank = (myRank + 1) % size;
    int thirdRank = (myRank + 2) % size;
    PartGID otherID(otherRank);
    PartGID thirdID(thirdRank);
    ::Entity* fakeEnt0 = reinterpret_cast<::Entity*>(0xabcdef);
    ::GRUMMP::EntityGID fakeEGID0(otherID, fakeEnt0);

    ::Entity* fakeEnt1 = reinterpret_cast<::Entity*>(0x123456);
    ::GRUMMP::EntityGID fakeEGID1(thirdID, fakeEnt1);

    ::EntityGID fakeEGID1a(myID, fakeEnt1);

    BOOST_CHECK(g.addCopy(fakeEGID0));
    BOOST_CHECK(g.addCopy(fakeEGID1));
    BOOST_CHECK(!g.addCopy(fakeEGID1a));
    BOOST_CHECK_EQUAL(g.getNumCopies(), 2);

    test = (g.getCopyID(0) == fakeEGID0);
    BOOST_CHECK(test);

    g.resetMaxNumCopies(5);
    BOOST_CHECK_EQUAL(g.getMaxNumCopies(), 5);
    BOOST_CHECK_EQUAL(g.getNumCopies(), 2);

    ::GRUMMP::Ghost g2(g);
    BOOST_CHECK_EQUAL(g2.getMaxNumCopies(), 5);
    BOOST_CHECK_EQUAL(g2.getNumCopies(), 2);
    test = (g2.getCopyID(1) == fakeEGID1);
    BOOST_CHECK(test);

    ::Entity* fakeEnt2 = reinterpret_cast<::Entity*>(0x654321);
    ::GRUMMP::Ghost g3(g, otherID, myID, fakeEnt2);
    ::GRUMMP::EntityGID fakeEGID2(myID, fakeEnt2);
    BOOST_CHECK_EQUAL(g3.getMaxNumCopies(), 5);
    BOOST_CHECK_EQUAL(g3.getNumCopies(), 2);
    returnedID = g3.getOwnerID();
    test = (returnedID == otherID);
    BOOST_CHECK(test);
    test = (g3.getCopyID(0) == fakeEGID2);
    BOOST_CHECK(test);

    char buffer[10240]; // Way bigger than needed.
    int offsetWrite = 0, offsetRead = 0;
    g.serialize(buffer, offsetWrite);
    ::GRUMMP::Ghost g4(myID);
    g4.setOwnerID(myID);
    BOOST_CHECK_EQUAL(g4.getMaxNumCopies(), 1);
    g4.deserialize(buffer, offsetRead);
    BOOST_CHECK_EQUAL(offsetWrite, offsetRead);

    returnedID = g4.getOwnerID();
    BOOST_CHECK(returnedID == myID);
    BOOST_CHECK_EQUAL(g4.getMaxNumCopies(), 2);
    BOOST_CHECK_EQUAL(g4.getNumCopies(), 2);
    test = (g4.getCopyID(0) == fakeEGID0);
    BOOST_CHECK(test);
    test = (g4.getCopyID(1) == fakeEGID1);
    BOOST_CHECK(test);

    BOOST_CHECK_EQUAL(g.MPIdataSize(), 40);

  }

  BOOST_AUTO_TEST_SUITE_END()
#endif

