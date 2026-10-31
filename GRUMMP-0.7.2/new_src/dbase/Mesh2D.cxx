#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <deque>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <iterator>
#include <utility>
using std::deque;
using std::list;

using std::map;
using std::multimap;
using std::pair;
using std::set;

#include "GR_GRGeom2D.h"
#include "GR_GRCoEdge.h"
#include "GR_GRCurve.h"
#include "GR_GRPoint.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"

#include "GR_Mesh2D.h"
#include "GR_CellCV.h"
#include "GR_Geometry.h"
#include "GR_InsertionQueue.h"
#include "GR_InsertionManager.h"
#include "GR_Vec.h"
#include "GR_misc.h"
#include "GR_BFace.h"
#include "GR_BFace2D.h"
#include "GR_Sort.h"

#ifndef WIN_VC32
#include <unistd.h>
#else
#pragma warning(disable:4786) // DEBUG symbol > 255 chars
#endif

// This is here so that configure can find the library.
extern "C"
{
  void
  vGRUMMP_2D_Lib(void)
  {
  }
}

//@ Default constructor
Mesh2D::Mesh2D(const int iQualMeas) :
    Mesh(), m_ECEdgeF(), m_ECEdgeBF(), m_ECTri(), m_ECQuad(), m_ECIntEdgeBF()
{
  m_qual = new Quality(this, iQualMeas);
}

void
createMeshInBox(Mesh2D& M2D, double dXLo, double dYLo, double dXHi, double dYHi)
{
  Vert* vertLL = M2D.createVert(dXLo, dYLo);
  Vert* vertLR = M2D.createVert(dXHi, dYLo);
  Vert* vertHL = M2D.createVert(dXLo, dYHi);
  Vert* vertHR = M2D.createVert(dXHi, dYHi);

  // Set up initial mesh
  //
  //    HL-------0--------HR
  //     |              _/|
  //     |     0     __/  |
  //     |         _/     |
  //     1      _4/       3
  //     |    _/     1    |
  //     | __/            |
  //     |/               |
  //    LL-------2--------LR
  //
  // Setup sizes for initial connectivity arrays
  M2D.createTriCell(vertLL, vertLR, vertHR);
  M2D.createTriCell(vertLL, vertHR, vertHL);
  M2D.createBFace(findCommonFace(vertLL, vertLR));
  M2D.createBFace(findCommonFace(vertLR, vertHR));
  M2D.createBFace(findCommonFace(vertHR, vertHL));
  M2D.createBFace(findCommonFace(vertHL, vertLL));

  assert(M2D.isValid());
}

void
Mesh2D::markCleanNeighbors(Cell * pC)
{
  int iReg = pC->getRegion();
  assert(iReg != iInvalidRegion);
  // Cells in the exterior do not get marked with a region.  Instead,
  // they retain their invalid region tag and no traversal is done.
  for (int iF = 0; iF < pC->getNumFaces(); iF++) {
    Face *pF = pC->getFace(iF);
    if (pF->isValid()) {
      Cell *pCCand = pF->getOppositeCell(pC);
      assert(pCCand->isValid());
      if (pCCand->getType() == Cell::eBdryEdge)
	continue;
      if (pCCand->getRegion() == iInvalidRegion) {
	// If the opposite cell hasn't been marked, mark it the same as
	// this one.
	pCCand->setRegion(iReg);
	markCleanNeighbors(pCCand);
      }				// Done tagging candidate cell
    }				// Face is OK
  }				// Loop over faces
}

//@ Retrieval by index
Face*
Mesh2D::getFace(const GR_index_t i) const
{
  assert(i < getNumFaces());
  return m_ECEdgeF.getEntry(i);
}

BFace*
Mesh2D::getBFace(const GR_index_t i) const
{
  assert(i < getNumTotalBdryFaces());
  if (i < getNumBdryFaces())
    return m_ECEdgeBF.getEntry(i);
  else
    return m_ECIntEdgeBF.getEntry(i - getNumBdryFaces());
}

IntBdryEdge*
Mesh2D::getIntBFace(const GR_index_t i) const
{
  assert(i < getNumIntBdryFaces());
  return m_ECIntEdgeBF.getEntry(i);
}

Cell*
Mesh2D::getCell(const GR_index_t i) const
{
  assert(i < getNumCells());
  if (i < m_ECTri.lastEntry())
    return m_ECTri.getEntry(i);
  else
    return m_ECQuad.getEntry(i - m_ECTri.lastEntry());
}

QuadCell*
Mesh2D::getQuadCell(const GR_index_t i) const
{
  assert(i < getNumQuadCells());
  return m_ECQuad.getEntry(i);
}

/// If the face pointer is -not- a valid pointer to a face in this mesh,
/// then the function returns -1.
/// @param pF A pointer to a face
/// @seealso iBFaceIndex
/// @seealso iCellIndex
/// @seealso iVertIndex
/// @return Integer index of the face
GR_index_t
Mesh2D::getFaceIndex(const Face * const pF) const
{
  assert(pF->getType() == Face::eEdgeFace);
  return (m_ECEdgeF.getIndex(dynamic_cast<const EdgeFace *>(pF)));
}

GR_index_t
Mesh2D::getBFaceIndex(const BFace * const pBF) const
{
  switch (pBF->getType())
    {
    case Cell::eBdryEdge:
      return m_ECEdgeBF.getIndex(dynamic_cast<const BdryEdge*>(pBF));
    case Cell::eIntBdryEdge:
      return getNumBdryFaces()
	  + m_ECIntEdgeBF.getIndex(dynamic_cast<const IntBdryEdge*>(pBF));
    default:
      assert(0);
      return INT_MAX;
    }
}

GR_index_t
Mesh2D::getCellIndex(const Cell * const pC) const
{
  assert(
      (pC->getType() == Cell::eQuadCell && !m_simplexMesh)
	  || pC->getType() == Cell::eTriCell);
  if (pC->getType() == Cell::eTriCell)
    return m_ECTri.getIndex(dynamic_cast<const TriCell *>(pC));
  else {
    return m_ECQuad.getIndex(dynamic_cast<const QuadCell *>(pC))
	+ m_ECTri.lastEntry();
  }
}				//@ Creation of new entities

#ifdef ITAPS
bool
Mesh2D::isValidEntHandle(void* pvEnt) const
// Determines if the memory location pointed to by the pointer
// given is valid (within the range) for any of the entity types
// for this mesh
{
  return (m_ECEdgeF.contains(static_cast<EdgeFace*>(pvEnt))
      || m_ECTri.contains(static_cast<TriCell*>(pvEnt))
      || m_ECQuad.contains(static_cast<QuadCell*>(pvEnt))
      || m_ECVerts.contains(static_cast<Vert*>(pvEnt)));
}
#endif

//@ Validity checking of entire mesh
bool
Mesh2D::isValid() const
{

  bool qRetVal = true;
  int i;

  for (i = iNumVerts() - 1; i >= 0 && qRetVal; i--) {
    Vert *pV = getVert(i);
    qRetVal = qRetVal && pV->isValid();
    // If the vertex still belongs to the mesh...
    if (!pV->isDeleted()) {
      // and the vertex has a hint face defined...
      Face *pF = pV->getHintFace();
      if (pF->isValid())
	// that hint has to be correct.
	qRetVal = qRetVal && pF->hasVert(pV);
    }
  }
  for (i = getNumFaces() - 1; i >= 0 && qRetVal; i--) {
    Face *pF = getFace(i);
    if (!(pF->isValid())) {
      qRetVal = false;
      break;
    }

    if (pF->isDeleted())
      continue;
    // Faces on the part bdry or that are actually owned elsewhere
    // have more relaxed rules about the required topology.
    if ((pF->getEntStatus() == Entity::ePartBdry
	|| pF->getEntStatus() == Entity::eForeign)
	&& (!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid()))
      continue;
    if (!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid()) {
      if (pF->getFaceLoc() == Face::eInterior) {
	qRetVal = false;
	break;
      }
      else {
	continue;
      }
    }
    int iLoc = pF->getFaceLoc();
    switch (iLoc)
      {
      case (Face::eInterior):
	if (pF->getLeftCell()->getRegion() != pF->getRightCell()->getRegion()) {
	  qRetVal = false;
	}
	break;
      case (Face::eBdryTwoSide):
	//if (pF->pCCellLeft()->iRegion() == pF->pCCellRight()->iRegion())
	//	qRetVal = false;
	//Must make sure that at least one side is an IntBdryEdge...
	if ((pF->getLeftCell()->getType() != Cell::eIntBdryEdge)
	    && (pF->getRightCell()->getType() != Cell::eIntBdryEdge)) {
	  qRetVal = false;
	}

	break;
      default:
	break;
      }
    if (pF->doFullCheck() == 0)
      qRetVal = false;

  }
  for (i = getNumBdryFaces() - 1; i >= 0 && qRetVal; i--)
    qRetVal = qRetVal && getBFace(i)->isValid();
  for (i = getNumCells() - 1; i >= 0 && qRetVal; i--) {
    Cell *pC = getCell(i);
    qRetVal = qRetVal && pC->isValid();

    if (pC->isDeleted())
      continue;
    bool qRegOK = (pC->getRegion() != iInvalidRegion);
    qRetVal = qRetVal && qRegOK;

    // Check that the cell is actually closed
    qRetVal = qRetVal && pC->isClosed();
  }

  //   qRetVal = qRetVal && qWatertight();
  return qRetVal;
}

// TO DO: Make these non-member functions

//@ Computation of perimeter, total surface area, and watertightness of surface
double
Mesh2D::calcBoundarySize() const
{
  double dRetVal = 0;
  for (int i = getNumBdryFaces() - 1; i >= 0; i--)
    dRetVal += getBFace(i)->calcSize();
  return dRetVal;
}

double
Mesh2D::calcInteriorSize() const
{
  double dRetVal = 0;
  for (int i = getNumCells() - 1; i >= 0; i--)
    dRetVal += getCell(i)->calcSize();
  return dRetVal;
}

bool
Mesh2D::isWatertight() const
{
  double adSum[] =
    { 0., 0., 0. };
  double dArea = 0;
  bool qAnyCells = false;
  for (int i = getNumCells() - 1; i >= 0; i--) {
    Cell *pC = getCell(i);
    // Only compute the directed surface area for cells which have valid
    // connectivity.  This way, surface meshes with deleted cells (ie,
    // no valid faces) can be checked for watertightness.
    if (pC->doFullCheck() == 1) {
      qAnyCells = true;
      double adTemp[3];
      pC->calcVecSize(adTemp);
      adSum[0] += adTemp[0];
      adSum[1] += adTemp[1];
      adSum[2] += adTemp[2];
      dArea += pC->calcSize();
    }
  }
  if (!qAnyCells)
    return true;
  double dInvArea = 1. / dArea;
  vSCALE3D(adSum, dInvArea);
  if (iFuzzyComp(adSum[0], 0) == 0 && iFuzzyComp(adSum[1], 0) == 0
      && iFuzzyComp(adSum[2], 0) == 0)
    return (true);
  else {
    logMessage(2, "Directed surface areas are non-zero:");
    logMessage(2, "  x: %12.5g\n  y: %12.5g\n  z: %12.5g\n", adSum[0], adSum[1],
	       adSum[2]);
    return (false);
  }
}

void
Mesh2D::identifyVertexTypesGeometrically() const
// In addition to some topological deductions, this routine also
// detects verts at which the boundary has a sharp corner.  This is
// critical for keeping boundaries from degenerating with repeated
// coarsening.
//
// Note that no effort at all is made to define an apex as a point
// where two patches meet, because a mesh input from file will have
// one patch per bdry edge.  This would make every bdry vert a
// eBdryApex, and that's no good.  However, boundary condition changes
// are detected, and verts where the BC changes are tagged as
// eBdryApex.

{
  // Tentatively set the type of every vertex to eInterior.  Then modify
  // that to accomodate boundary verts (both internal and regular).
	for (GR_index_t iV = 0; iV < iNumVerts(); iV++) {
		Vert *pV = getVert(iV);
		if (pV->getVertType() != Vert::eBdryApex && pV->getVertType() != Vert::ePseudoSurface)
			pV->setType(Vert::eInterior);
	}

  GR_index_t iBF;
  // Mark regular boundaries.  In the process, tag intersections of
  // regular and internal boundaries as type eBdryApex.  Finally, make a
  // sorted list of bdry normals (sorted by vertex) for both regular and
  // internal boundaries.
  std::multiset<BdryNorm> msBN;
  for (iBF = 0; iBF < getNumBdryFaces(); iBF++) {
    BFace *pBF = getBFace(iBF);
    Face *pF = pBF->getFace();

    Vert *pV0 = pF->getVert(0);
    Vert *pV1 = pF->getVert(1);
    pV0->setType(Vert::eBdry);
    pV1->setType(Vert::eBdry);
    double adNorm[3];
    pBF->calcUnitNormal(adNorm);
    msBN.insert(BdryNorm(pV0, adNorm));
    msBN.insert(BdryNorm(pV1, adNorm));
  }

  for (iBF = 0; iBF < getNumIntBdryFaces(); iBF++) {
    IntBdryEdge *pIBE = getIntBFace(iBF);
    Face *pF = pIBE->getFace(0);
    Vert *pV0 = pF->getVert(0);
    Vert *pV1 = pF->getVert(1);
    pV0->setType(Vert::eBdryTwoSide);
    pV1->setType(Vert::eBdryTwoSide);
    double adNorm[3];
    pIBE->calcUnitNormal(adNorm);
    msBN.insert(BdryNorm(pV0, adNorm));
    msBN.insert(BdryNorm(pV1, adNorm));
  }
  assert(msBN.size() == getNumBdryFaces() * 2 + getNumIntBdryFaces() * 2);

  // At this point, msBN is a sorted set of BdryNorm's, with multiple
  // entries having the same vertex.
  GR_index_t iRecordsChecked = 0;
  std::multiset<BdryNorm>::iterator itFirst, itLast;

  itFirst = itLast = msBN.begin();
  do {
    Vert *pV = itFirst->pVVert();
    // Should already be tagged as a bdry vert.
    assert(pV->getVertType() != Vert::eInterior);
    int numNormals = 1;
    do {
      itLast++;
      numNormals++;
    }
    while (itLast != msBN.end() && itLast->pVVert() == pV);
    itLast--;
    numNormals--;
    assert(itLast->pVVert() == pV);

    // A vertex that's on a part bdry and a domain / internal bdry
    // will have iDist == 1.
    // TODO: Add a check that, for this case, the vertex is really
    // on a part bdry.  Currently non-trivial, because verts can't
    // be tagged yet.

    switch (numNormals)
      {
      case 0:
	// Should be impossible
	assert(0);
	break;
      case 1:
	// Do nothing
	break;
      case 2:  // Could be an apex; depends on the normals
	{
	  // Now check how close their normals are to being parallel.
	  // Deviations of as much as 20 degrees are permitted.
	  double dFirstMag = dMAG2D(itFirst->adNormal());
	  double dLastMag = dMAG2D(itLast->adNormal());
	  double dDot = dDOT2D(itFirst->adNormal(), itLast->adNormal())
	      / dFirstMag / dLastMag;
	  if (dDot < 0.9999999) {
	    // This one is now an apex!
	    pV->setType(Vert::eBdryApex);
	  }
	  break;
	}
      default:
	pV->setType(Vert::eBdryApex);
	assert(pV->getVertType() == Vert::eBdryApex);
	break;
      }
    iRecordsChecked += numNormals;
    itFirst = ++itLast;
  }
  while (iRecordsChecked < 2 * getNumBdryFaces() + 2 * getNumIntBdryFaces());

   // Below for loop with tag faces as PseudoSurface if the verts of those faces are PseudoSurface
   for (GR_index_t i = 0; i < getNumFaces(); i++) {
		Face *face = getFace(i);
		Vert *v0 = face->getVert(0);
		Vert *v1 = face->getVert(1);
		if ((v0->getVertType() == 5 && v1->getVertType() == 5)|| (v0->getVertType() == 1 && v1->getVertType() == 5)) {
			face->setFaceLoc(Face::ePseudoSurface); //Tags the faces as ePseudosurface...prevents faces from swapping
			assert(face->isValid());
			face->lock();   //prevents faces from swapping
		}
		if (face->getFaceLoc() == 3) {
			//ARM2D.m_SI2D->insertPointOnFace(insertionPoint, FaceToSplit);
		//	face->printFaceInfo(); //  Print to check if faces are tagged as ePseudosurface
		}
	}
}

void
Mesh2D::buildBdryPatches()
// Given a Mesh2D (but no knowledge of the underlying geometry), make
// a new set of patches, one per bdry face.  This is used, in
// particular, to re-create bdry patch info after mesh coarsening.
{

  vFatalError("The function is no longer defined",
	      "Mesh2D::vBuildBdryPatches()");

  //   int iBF;
  //   int iNBFaces = iNumTotalBdryFaces();
  //   int iNVerts = iNumVerts();

  //   // Record the boundary condition and region info.
  //   int *aiRightReg = new int[iNBFaces];
  //   int *aiLeftReg  = new int[iNBFaces];
  //   int *aiRightBC  = new int[iNBFaces];
  //   int *aiLeftBC   = new int[iNBFaces];

  //   for (iBF = 0; iBF < iNBFaces; iBF++) {
  //     BFace *pBF = pBFBFace(iBF);
  //     Face *pF = pBF->pFFace(0);
  //     Cell *pC = pF->pCCellOpposite(pBF);
  //     assert(pC->eType() == Cell::eTriCell);

  //     Vert *pV0 = pBF->pVVert(0);
  //     Vert *pV1 = pBF->pVVert(1);
  //     Vert *pVOpp = dynamic_cast<TriCell*>(pC)->pVVertOpposite(pF);

  //     int iOrient = iOrient2D(pV0, pV1, pVOpp);
  //     bool qIntBFace = (pBF->iNumFaces() == 2);

  //     if (qIntBFace) {
  //       // Internal bdry face!

  //       // FIX ME!  Untested code, because coarsening with internal bdrys
  //       // hasn't been written yet.

  //       assert(pBF->eType() == Cell::eIntBdryEdge);
  //       // Identify the other cell
  //       Face *pFAlt = pBF->pFFace(1);
  //       Cell *pCAlt = pFAlt->pCCellOpposite(pBF);

  //       // Get regions on both sides.
  //       switch (iOrient) {
  //       case 1:
  // 	aiLeftReg[iBF]  = pC   ->iRegion();
  // 	aiRightReg[iBF] = pCAlt->iRegion();
  // 	break;
  //       case -1:
  // 	aiLeftReg[iBF]  = pCAlt->iRegion();
  // 	aiRightReg[iBF] = pC   ->iRegion();
  // 	break;
  //       default:
  // 	// Unreachable for valid meshes
  // 	assert(0);
  // 	break;
  //       }
  //       aiRightBC[iBF] = aiLeftBC[iBF] = iInvalidBC;
  //     }
  //     else { // Regular bdry face
  //       aiRightReg[iBF] = aiLeftReg[iBF] = iInvalidRegion;
  //       aiRightBC[iBF]  = aiLeftBC[iBF]  = iInvalidBC;

  //       // Get regions on both sides.
  //       switch (iOrient) {
  //       case 1:
  // 	aiLeftReg[iBF]  = pC ->iRegion();
  // 	aiRightBC[iBF]  = pBF->iBdryCond();
  // 	break;
  //       case -1:
  // 	aiRightReg[iBF] = pC ->iRegion();
  // 	aiLeftBC[iBF]   = pBF->iBdryCond();
  // 	break;
  //       default:
  // 	// Unreachable for valid meshes
  // 	assert(0);
  // 	break;
  //       }
  //     } // Done with real bdry faces
  //     // Do some post-checking to be sure this matches qualitatively with
  //     // the data for the patch.
  // #ifndef NDEBUG
  //     BdryPatch *pBP = pBF->pPatchPointer();
  //     assert( ( aiRightReg[iBF] == pBP->iRightRegion() ||
  // 	      aiRightReg[iBF] == pBP->iLeftRegion() )
  // 	    &&
  // 	    ( aiLeftReg[iBF]  == pBP->iRightRegion() ||
  // 	      aiLeftReg[iBF]  == pBP->iLeftRegion() ) );
  //     assert( ( aiRightBC[iBF] == pBP->iBdryCond() ||
  // 	      aiRightBC[iBF] == iInvalidBC )
  // 	    &&
  // 	    ( aiLeftBC[iBF]  == pBP->iBdryCond() ||
  // 	      aiLeftBC[iBF]  == iInvalidBC ) );
  // #endif
  //   } // End of loop over BFaces

  //   // Throw away all the existing bdry patch info.
  //   B2DInfo.vSetSize(iNBFaces, iNVerts);
  //   B2DInfo.vClearPoints();

  //   // Now begin again, creating a new bdry patch for every BdryEdge.

  //   // Need to map vertex indices to point indices in the Bdry2D
  //   // struct, since Bdry2D no longer uses vertices directly.  This
  //   // mapping is stored in the following array.
  //   int *aiVertMap = new int[iNVerts];
  //   int iV;
  //   for (iV = 0; iV < iNVerts; aiVertMap[iV++] = -1);

  //   int iBPt = 0;
  //   for (iBF = 0; iBF < iNBFaces; iBF++) {
  //     BFace *pBF = pBFBFace(iBF);
  //     int iV0 = iVertIndex(pBF->pVVert(0));
  //     int iV1 = iVertIndex(pBF->pVVert(1));

  //     if (aiVertMap[iV0] == -1) {
  //       aiVertMap[iV0] = iBPt;
  //       B2DInfo.vAddPoint(pVVert(iV0)->adCoords());
  //       B2DInfo.vSetBdryPoint(iBPt, false);
  //       iBPt++;
  //     }

  //     if (aiVertMap[iV1] == -1) {
  //       aiVertMap[iV1] = iBPt;
  //       B2DInfo.vAddPoint(pVVert(iV1)->adCoords());
  //       B2DInfo.vSetBdryPoint(iBPt, false);
  //       iBPt++;
  //     }

  //     int iBP0 = aiVertMap[iV0];
  //     int iBP1 = aiVertMap[iV1];
  //     assert(iBP0 >= 0 && iBP1 >= 0);

  //     BdryPatch2D *pBP =
  //       new BdrySeg(&B2DInfo, aiLeftBC[iBF], aiRightBC[iBF], aiLeftReg[iBF],
  // 		  aiRightReg[iBF], iBP0, iBP1);
  //     B2DInfo.vAddPatchToList(pBP);
  //     // vSetPatch can handle getting the patch for an internal bdry
  //     // fixed up, too.
  //     (dynamic_cast<BdryEdgeBase*>(pBF))->vSetPatch(iBF, &B2DInfo);
  //   } // Done with boundary cases
  //   delete [] aiVertMap;

  //   if (!(qValid ()))
  //     vFoundBug ("building 2D boundary data from mesh");
}

Cell*
Mesh2D::findCell(const double newPtLoc[], Cell * const cellGuess,
		 int &numViolated, Face *facesViolated[3], int &numTies,
		 Face *facesTied[3]) const
{
  Cell *cell = cellGuess, *cellNew = cellGuess, *cellOld = cellGuess;
  if (!cell->isValid()) {
    for (GR_index_t ii = 0; ii < getNumCells(); ii++) {
      cell = getCell(ii);
      if (!cell->isDeleted())
	break;
    }
  }
  assert(cell->doFullCheck());

  int iOrient = 0;
  Face *pF = pFInvalidFace;
  Vert VTemp;
  VTemp.setCoords(2, newPtLoc);

  do {
    int i;
    numViolated = 0;
    numTies = 0;
    for (i = 0; i < 3; i++) {
      pF = cell->getFace(i);
      // If the face is pointing away from this cell, then signs need to
      // be reversed to walk in the proper direction
      int iSign = (pF->getLeftCell() == cell) ? 1 : -1;
      iOrient = checkOrient2D(pF->getVert(0), pF->getVert(1), &VTemp) * iSign;
      // If the vertex lies outside the current cell and behind the
      // current face, walk in that direction.  If the walk would
      // involve going outside the boundary, continue checking to try to
      // find another face to walk across; this improves behavior of
      // walks which try to insert a boundary point.
      if (iOrient == -1) {
	cellNew = pF->getOppositeCell(cell);
	if ((cellNew->getType() == Cell::eBdryEdge)   || (pF->getFaceLoc() == Face::eBdryTwoSide))
	  facesViolated[numViolated++] = pF;
	else
	  break;
      }
      if (iOrient == 0)
	facesTied[numTies++] = pF;
      assert(numTies >= 0 && numTies <= 2);
      assert(numViolated >= 0 && numViolated <= 2);
    }
    cellOld = cell;
    cell = cellNew;
  }
  while (iOrient == -1 && cell->getType() != Cell::eBdryEdge  && pF->getFaceLoc() != Face::eBdryTwoSide);
  // Restore identity of last interior cell visited, in case it's needed
  cell = cellOld;

  return cell;
}
void
Mesh2D::writeTempMesh()
{
  if (m_writeTempMeshes) {
    char strName[FILE_NAME_LEN];
    sprintf(strName, "%.1010s_%d", m_tempMeshFileName, m_tempMeshIndex);
//		writeVTKLegacyWithoutPurge(*this, strName);
    writeVTK(*this, strName);
    m_tempMeshIndex++;
  }
}

bool
Mesh2D::DeterGoodTri(TriCell* cell, Face* face)
{
  double longestedge;
  double m_g = gradeparameter;
  double termination = determineparameter;
  double m_angle = minangle;
  double refbdrylength = face->getVert(0)->getLengthScale()
      + 1. / m_g * 0.5 * face->calcSize();
  refbdrylength =
      refbdrylength
	  < (face->getVert(1)->getLengthScale()
	      + 1. / m_g * 0.5 * face->calcSize()) ?
	  refbdrylength :
	  (face->getVert(1)->getLengthScale()
	      + 1. / m_g * 0.5 * face->calcSize());
  //Find the longest edge
  longestedge = cell->getFace(0)->calcSize();
  longestedge =
      cell->getFace(1)->calcSize() > longestedge ?
	  cell->getFace(1)->calcSize() : longestedge;
  longestedge =
      cell->getFace(2)->calcSize() > longestedge ?
	  cell->getFace(2)->calcSize() : longestedge;
  ///Average Lengthscale
//  	double lengthscale = (cell->getVert(0)->getLengthScale() + cell->getVert(1)->getLengthScale() + cell->getVert(2)->getLengthScale())/3.;
  double hcell = cell->calcCircumradius();
  double hcenter = 1.e9;
  double largestls = 0.;
  ///Find the largest lengthscale of three verts and the length scale at circumcenter.
  for (int ii = 0; ii < cell->getNumVerts(); ii++) {
    largestls = std::max(largestls, cell->getVert(ii)->getLengthScale());
    hcenter = std::min(hcenter,
		       hcell / m_g + cell->getVert(ii)->getLengthScale());
  }
  double refratio = cell->calcShortestEdgeLength()
      / (sqrt(3.) * cell->calcCircumradius());
  bool sizeflag = true;
  for (int i = 0; i < 3; i++) {
    if (iFuzzyComp(
	cell->getFace(i)->calcSize(),
	termination * 0.5
	    * (cell->getFace(i)->getVert(0)->getLengthScale()
		+ cell->getFace(i)->getVert(0)->getLengthScale())) < 1.) {
      sizeflag = sizeflag and true;
    }
    else
      sizeflag = sizeflag and false;
  }
//  	if (face->layerindex < 20){
//  		termination = determineparameter * pow((1.+1./m_g),face->layerindex);
//  	}
//  	if (face->isBdryFace()){
//  		if ((iFuzzyComp(refratio , 2.0*sin(m_angle)/sqrt(3.)) < 1.) or (iFuzzyComp(hcell , termination*hcenter) > 0.)){
  if ((iFuzzyComp(refratio, 2.0 * sin(m_angle) / sqrt(3.)) < 1.)
      or (sizeflag == false)) {
    return false;
  }
  else
    return true;
//  	}
//  	else {
//  		if ((iFuzzyComp(refratio , 2.0*sin(m_angle)/sqrt(3.)) < 1.) or (iFuzzyComp(hcell , termination*hcenter) > 0.)){
//  		if ((iFuzzyComp(refratio , 2.0*sin(m_angle)/sqrt(3.)) < 1.) or (sizeflag == false)){
//  			return false;
//  		}
//  		else
//  			return true;
//  	}
}

void
Mesh2D::makeSimplicial()
{
  if (m_simplexMesh)
    return;
  if (getNumQuadCells() == 0) {
    m_simplexMesh = 0;
    return;
  }
  else {
    // Split every quad into two triangles.
    GR_index_t numQuads = getNumQuadCells();
    for (GR_index_t ii = 0; ii < numQuads; ii++) {
      // Grab a quad and split it into two tris.
      QuadCell* qc = dynamic_cast<QuadCell*>(getQuadCell(ii));
      assert(qc);
      if (qc->isDeleted())
	continue;
      // These verts are guaranteed to be in cyclic order.
      Vert* verts[] =
	{ qc->getVert(0), qc->getVert(1), qc->getVert(2), qc->getVert(3) };
      deleteCell(qc);
      createTriCell(verts[0], verts[1], verts[2]);
      createTriCell(verts[0], verts[2], verts[3]);
    }
  }
}
