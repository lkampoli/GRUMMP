#ifndef VolMeshFixtureBase_h_
#define VolMeshFixtureBase_h_

#include "GR_Geometry.h"
#include "GR_VolMesh.h"

struct VolMeshFixtureBase : public VolMesh {
  VolMeshFixtureBase() :
      VolMesh()
  {
    setEncroachmentType(eBall);
    iMessageStdoutLevel = 0;
  }
  ~VolMeshFixtureBase()
  {
  }
protected:
  void
  checkMeshSize(const int nVerts, const int nFaces, const int nBdryFaces,
		const int nTris, const int nQuads, const int nCells,
		const int nTets, const int nPyrs, const int nPrisms,
		const int nHexes, const int nFlakes = 0)
  {
    BOOST_CHECK_EQUAL(nVerts, getNumVerts());
    BOOST_CHECK_EQUAL(nTris, getNumTrisInUse());
    BOOST_CHECK_EQUAL(nQuads, getNumQuadsInUse());
    BOOST_CHECK_EQUAL(nFaces, getNumQuadsInUse() + getNumTrisInUse());
    BOOST_CHECK_EQUAL(nBdryFaces, getNumBdryFacesInUse());
    BOOST_CHECK_EQUAL(nTets, getNumTetsInUse());
    BOOST_CHECK_EQUAL(nPyrs, getNumPyrsInUse());
    BOOST_CHECK_EQUAL(nPrisms, getNumPrismsInUse());
    BOOST_CHECK_EQUAL(nHexes, getNumHexesInUse());
    BOOST_CHECK_EQUAL(nFlakes, getNumFlakesInUse());
    BOOST_CHECK_EQUAL(
	nCells,
	getNumTetsInUse() + getNumPyrsInUse() + getNumPrismsInUse()
	    + getNumHexesInUse() + getNumFlakesInUse());
  }
};

#endif
