#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "GR_Geometry.h"
#include "GR_InsertionManager.h"
#include "GR_Sort.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_VolMesh.h"

#include <algorithm>
#include <map>

using std::map;
using std::multimap;
using std::pair;

/// This file contains, for 3D meshes: (most) constructors, basic queries,

// This is here so that configure can find the library.
extern "C"
{
  void
  vGRUMMP_3D_Lib(void)
  {
  }
}

/// VolMesh constructors

VolMesh::VolMesh(const int iQualMeas, const double dMaxAngle) :
    Mesh(), m_OKToSwapSurfaceEdges(true), m_doEdgeSwaps(true), m_strictPatchChecking(
	true), m_maxAngleForSurfSwap(dMaxAngle), m_ECTriF(), m_ECQuadF(), m_ECTriBF(), m_ECQuadBF(), m_ECTet(), m_ECPyr(), m_ECPrism(), m_ECHex(), m_ECIntTriBF(), m_ECIntQuadBF(), m_Bdry3D(), m_bdryFromMesh(
	true), m_SubsegMap(new SubsegMap())
{
  m_qual = new Quality(this, iQualMeas);
  m_allowBdryChanges = true;
}

//VolMesh::VolMesh(const VolMesh& VM)
//  : Mesh (VM),
//    m_OKToSwapSurfaceEdges(VM.m_OKToSwapSurfaceEdges),
//    m_doEdgeSwaps(VM.m_doEdgeSwaps),
//    m_strictPatchChecking(VM.m_strictPatchChecking),
//    m_maxAngleForSurfSwap(VM.m_maxAngleForSurfSwap),
//    m_ECTriF(), m_ECQuadF(), m_ECTriBF(), m_ECQuadBF(),
//    m_ECTet(), m_ECPyr(), m_ECPrism(), m_ECHex(),
//    m_ECIntTriBF(), m_ECIntQuadBF(), m_Bdry3D(VM.m_Bdry3D),
//    m_bdryFromMesh(VM.m_bdryFromMesh)
//{
//  m_allowBdryChanges = VM.m_allowBdryChanges;
//  m_isLengthScaleFromCellSizes = VM.m_isLengthScaleFromCellSizes;
//  m_fillBdry = VM.m_fillBdry;
//  m_interiorSkip = VM.m_interiorSkip;
//
//  m_qual = new Quality(*VM.m_qual);
//
//  if (m_simplexMesh) {
//    assert(VM.getNumQuadFaces() == 0);
//    assert(VM.getNumPyrCells() == 0);
//    assert(VM.getNumPrismCells() == 0);
//    assert(VM.getNumHexCells() == 0);
//    assert(VM.getNumQuadBdryFaces() == 0);
//    assert(VM.getNumIntQuadBdryFaces() == 0);
//  }
//
//  // Create faces first, and then cells.  This should preserve entity ordering.
//  for (GR_index_t iF = 0; iF < VM.getNumFaces(); iF++) {
//    Vert *apV[4];
//    apV[3] = pVInvalidVert;
//    Face *pFOld = VM.getFace(iF);
//    int iNV = pFOld->getNumVerts();
//    for (int iV = 0; iV < iNV; iV++) {
//      apV[iV] = getVert(VM.getVertIndex(pFOld->getVert(iV)));
//    }
//    bool qExist;
//    Face *pFNew = createFace(qExist, apV[0], apV[1], apV[2], apV[3]);
//    assert(!qExist && pFNew->isValid());
//  }
//  assert(getNumFaces() == VM.getNumFaces() &&
//	 getNumTriFaces() == VM.getNumTriFaces());
//
//  for (GR_index_t iT = 0; iT < VM.getNumTetCells(); iT++) {
//    Vert *apV[4];
//    Cell *pCOld = VM.getCell(iT);
//    for (int iV = 0; iV < 4; iV++) {
//      apV[iV] = getVert(VM.getVertIndex(pCOld->getVert(iV)));
//    }
//    bool qExist;
//    TetCell *pTC = createTetCell(qExist, apV[0], apV[1], apV[2], apV[3]);
//    assert(!qExist && pTC->isValid());
//  }
//
//  int iOffset = VM.getNumTetCells();
//  for (GR_index_t iP = 0; iP < VM.getNumPyrCells(); iP++) {
//    Vert *apV[5];
//    Cell *pCOld = VM.getCell(iP + iOffset);
//    for (int iV = 0; iV < 5; iV++) {
//      apV[iV] = getVert(VM.getVertIndex(pCOld->getVert(iV)));
//    }
//    bool qExist;
//    PyrCell *pPC = createPyrCell(qExist, apV[0], apV[1], apV[2], apV[3],
//				 apV[4]);
//    assert(!qExist && pPC->isValid());
//  }
//
//  iOffset += VM.getNumPyrCells();
//  for (GR_index_t iP = 0; iP < VM.getNumPrismCells(); iP++) {
//    Vert *apV[6];
//    Cell *pCOld = VM.getCell(iP + iOffset);
//    for (int iV = 0; iV < 6; iV++) {
//      apV[iV] = getVert(VM.getVertIndex(pCOld->getVert(iV)));
//    }
//    bool qExist;
//    PrismCell *pPC = createPrismCell(qExist, apV[0], apV[1], apV[2], apV[3],
//				     apV[4], apV[5]);
//    assert(!qExist && pPC->isValid());
//  }
//
//  iOffset += VM.getNumPrismCells();
//  for (GR_index_t iH = 0; iH < VM.getNumHexCells(); iH++) {
//    Vert *apV[8];
//    Cell *pCOld = VM.getCell(iH + iOffset);
//    for (int iV = 0; iV < 8; iV++) {
//      apV[iV] = getVert(VM.getVertIndex(pCOld->getVert(iV)));
//    }
//    bool qExist;
//    HexCell *pHC = createHexCell(qExist, apV[0], apV[1], apV[2], apV[3],
//				 apV[4], apV[5], apV[6], apV[7]);
//    assert(!qExist && pHC->isValid());
//  }
//
//  for (GR_index_t iBF = 0; iBF < VM.getNumBdryFaces(); iBF++) {
//    Face *pF;
//    BFace *pBFOld = VM.getBFace(iBF);
//    pF = getFace(VM.getFaceIndex(pBFOld->getFace()));
//
//    // This call may not work properly in this context: somehow the
//    // -new- Bdry3D data has to get in here, and I'm not sure that's
//    // even supported... depends on how copy construction of Bdry3D data
//    // is done.
//    createBFace(pF, pBFOld);
//  }
//
//  setAllHintFaces();
//
//  assert(isValid());
//  if (!(isValid()))
//    vFoundBug("VolMesh copying");
//}

VolMesh::~VolMesh()
{
  if (m_SubsegMap) {
    delete m_SubsegMap;
    m_SubsegMap = NULL;

  }
  if (m_bdryFromMesh) {
    delete m_Bdry3D;       /// Bdry info
    m_Bdry3D = NULL;
  }
}

void
VolMesh::setBdryData(Bdry3D* bdryData)
{
  if (bdryData == NULL)
    return;
  if (m_bdryFromMesh) {
    delete m_Bdry3D;
    m_bdryFromMesh = false;
  }
  m_Bdry3D = bdryData;
}

//@ Retrieval of entities by index
Face*
VolMesh::getFace(const GR_index_t i) const
{
  assert(i < getNumFaces());
  if (i < getNumTriFaces())
    return m_ECTriF.getEntry(i);
  else
    return m_ECQuadF.getEntry(i - getNumTriFaces());
}

BFace*
VolMesh::getBFace(const GR_index_t i) const
{
  assert(i < getNumTotalBdryFaces());
  GR_index_t iNTriBF = getNumTriBdryFaces();
  if (i < iNTriBF)
    return m_ECTriBF.getEntry(i);

  GR_index_t iNRegBF = iNTriBF + getNumQuadBdryFaces();
  assert(iNRegBF == getNumBdryFaces());
  if (i < iNRegBF)
    return m_ECQuadBF.getEntry(i - iNTriBF);

  GR_index_t iNTQIntTriBF = iNRegBF + getNumIntTriBdryFaces();
  if (i < iNTQIntTriBF)
    return m_ECIntTriBF.getEntry(i - iNRegBF);

  else
    return m_ECIntQuadBF.getEntry(i - iNTQIntTriBF);
}

BFace*
VolMesh::getIntBFace(const GR_index_t i) const
{
  GR_index_t iNIntTriBF = getNumIntTriBdryFaces();
  if (i < iNIntTriBF)
    return m_ECIntTriBF.getEntry(i);

  else
    return m_ECIntQuadBF.getEntry(i - iNIntTriBF);
}

Cell*
VolMesh::getCell(const GR_index_t i) const
{
  assert(i < getNumCells());
  GR_index_t iNTet = getNumTetCells();

  if (i < iNTet)
    return m_ECTet.getEntry(i);

  GR_index_t iNTPyr = iNTet + getNumPyrCells();
  if (i < iNTPyr)
    return m_ECPyr.getEntry(i - iNTet);

  GR_index_t iNTPPrism = iNTPyr + getNumPrismCells();
  if (i < iNTPPrism)
    return m_ECPrism.getEntry(i - iNTPyr);

  GR_index_t iNTPPHex = iNTPPrism + getNumHexCells();
  if (i < iNTPPHex)
    return m_ECHex.getEntry(i - iNTPPrism);

  else
    return m_ECFlake.getEntry(i - iNTPPHex);
}

GR_index_t
VolMesh::getBFaceIndex(const BFace* const pBF) const
{
  switch (pBF->getType())
    {
    case Cell::eTriBFace:
      return m_ECTriBF.getIndex(dynamic_cast<const TriBFace*>(pBF));
    case Cell::eQuadBFace:
      assert(!m_simplexMesh);
      return m_ECQuadBF.getIndex(dynamic_cast<const QuadBFace*>(pBF))
	  + m_ECTriBF.lastEntry();
    case Cell::eIntTriBFace:
      return m_ECIntTriBF.getIndex(dynamic_cast<const IntTriBFace*>(pBF))
	  + getNumBdryFaces();
    case Cell::eIntQuadBFace:
      assert(!m_simplexMesh);
      return m_ECIntQuadBF.getIndex(dynamic_cast<const IntQuadBFace*>(pBF))
	  + getNumBdryFaces() + m_ECIntTriBF.lastEntry();
    default:
      assert(0);
      return (INT_MAX);
    }
}

GR_index_t
VolMesh::getCellIndex(const Cell* const pC) const
{
  assert(
      pC->getType() == Cell::eTet
	  || (!m_simplexMesh
	      && (pC->getType() == Cell::ePyr || pC->getType() == Cell::ePrism
		  || pC->getType() == Cell::eHex)));
  if (pC->getType() == Cell::eTet)
    return m_ECTet.getIndex(dynamic_cast<const TetCell*>(pC));
  else if (pC->getType() == Cell::ePyr)
    return m_ECPyr.getIndex(dynamic_cast<const PyrCell*>(pC))
	+ m_ECTet.lastEntry();
  else if (pC->getType() == Cell::ePrism)
    return m_ECPrism.getIndex(dynamic_cast<const PrismCell*>(pC))
	+ m_ECPyr.lastEntry() + m_ECTet.lastEntry();
  else {
    return m_ECHex.getIndex(dynamic_cast<const HexCell*>(pC))
	+ m_ECPrism.lastEntry() + m_ECPyr.lastEntry() + m_ECTet.lastEntry();
  }
}

GR_index_t
VolMesh::getFaceIndex(const Face* const pF) const
{
  assert(
      pF->getType() == Face::eTriFace
	  || (pF->getType() == Face::eQuadFace && !m_simplexMesh));
  if (pF->getType() == Face::eTriFace)
    return m_ECTriF.getIndex(dynamic_cast<const TriFace*>(pF));
  else {
    return m_ECQuadF.getIndex(dynamic_cast<const QuadFace*>(pF))
	+ m_ECTriF.lastEntry();
  }
}

#ifdef ITAPS
bool
VolMesh::isValidEntHandle(void* pvEnt) const
// Determines if the memory location pointed to by the pointer
// given is valid (within the range) for any of the entity types
// for this mesh
{
  bool qRetVal = false;
  switch (static_cast<Entity*>(pvEnt)->getEntTopology())
    {
    case iMesh_POINT:
      qRetVal = m_ECVerts.contains(static_cast<Vert*>(pvEnt));
      break;
    case iMesh_TRIANGLE:
      qRetVal = m_ECTriF.contains(static_cast<TriFace*>(pvEnt));
      break;
    case iMesh_QUADRILATERAL:
      qRetVal = m_ECQuadF.contains(static_cast<QuadFace*>(pvEnt));
      break;
    case iMesh_TETRAHEDRON:
      qRetVal = m_ECTet.contains(static_cast<TetCell*>(pvEnt));
      break;
    case iMesh_PYRAMID:
      qRetVal = m_ECPyr.contains(static_cast<PyrCell*>(pvEnt));
      break;
    case iMesh_PRISM:
      qRetVal = m_ECPrism.contains(static_cast<PrismCell*>(pvEnt));
      break;
    case iMesh_HEXAHEDRON:
      qRetVal = m_ECHex.contains(static_cast<HexCell*>(pvEnt));
      break;
    case iMesh_FLAKE:
      qRetVal = m_ECFlake.contains(static_cast<FlakeCell*>(pvEnt));
      break;
    default:
      // All others are definitely invalid.
      break;
    }
  return qRetVal;
}
#endif

//@ Validity checking of entire mesh
bool
VolMesh::isValid() const
{
  bool qRetVal = true;
  int i;
  for (i = getNumVerts() - 1; i >= 0 && qRetVal; i--) {
    Vert *pV = getVert(i);
    qRetVal = qRetVal && pV->isValid();
    if (!pV->isDeleted()) {
      Face *pF = pV->getHintFace();
      if (pF->isValid())
	qRetVal = qRetVal && pF->hasVert(pV);
    }
  }

  for (i = getNumFaces() - 1; i >= 0 && qRetVal; i--) {
    Face *pF = getFace(i);
    if (pF->isValid() || pF->isDeleted())
      continue;
    if (pF->doFullCheck() == 0)
      qRetVal = false;
  }

  for (i = getNumBdryFaces() - 1; i >= 0 && qRetVal; i--)
    qRetVal = qRetVal && getBFace(i)->isValid();
  for (i = getNumCells() - 1; i >= 0 && qRetVal; i--) {
    Cell *pC = getCell(i);
    if (!pC->isValid() || pC->isDeleted())
      continue;
    int iFull = pC->doFullCheck();
    bool qRegOK = (pC->getRegion() != iInvalidRegion);
    qRetVal = qRetVal && (iFull && qRegOK);
    // Check that the cell is actually closed
    qRetVal = qRetVal && pC->isClosed();
    // Provided that the cell actually exists, be sure that it has a
    // positive volume.
    // Test temporarily disabled to confirm topo validity without geom check.
    if (iFull == 1 && pC->getType() == CellSkel::eTet)
      qRetVal = qRetVal && (pC->calcSize() >= 0);
  }

  //    qRetVal = qRetVal && qWatertight();
  return qRetVal;
}

bool
VolMesh::isMeshOrientationOK(const bool qFix)
// Check and correct some bad orientation situations.  Specifically,
// if a mesh is specified with all tets left-handed, this will be
// detected and corrected.  More general cases, including isolated
// inverted tetrahedra in the interior of the mesh, are not yet
// handled.  Also, meshes with degenerate tetrahedra are not fixed.
{
  // The status of the two tets on either side of a face is summarized
  // by the formula 5*(left-status) + right-status, where *-status is -1
  // for reverse-oriented tets, 0 for flat tets, 1 for
  // correctly-oriented tets, 2 for boundaries, and 3 for non-tets.  The
  // following enum gives symbolic values to the results.
  enum {
    eInvInv = -6,
    eInvFlat = -5,
    eInvGood = -4,
    eInvBdry = -3,
    eInvNontet = -2,
    eFlatInv = -1,
    eFlatFlat = 0,
    eFlatGood = 1,
    eFlatBdry = 2,
    eFlatNontet = 3,
    eGoodInv = 4,
    eGoodFlat = 5,
    eGoodGood = 6,
    eGoodBdry = 7,
    eGoodNontet = 8,
    eBdryInv = 9,
    eBdryFlat = 10,
    eBdryGood = 11,
    eBdryBdry = 12,
    eBdryNontet = 13,
    eNontetInv = 14,
    eNontetFlat = 15,
    eNontetGood = 16,
    eNontetBdry = 17,
    eNontetNontet = 18
  };
  GR_index_t i, iGood, iBad, iFixed, iUnfixed;
  // Check is performed on a face-by-face basis.  Performing this check
  // on a cell-by-cell basis for tets is not necessarily thorough
  // enough, as some faces will never be used to check an orientation.
  // A more general solution to this problem may have to be found for
  // non-tetrahedral cells.
  iGood = iBad = iFixed = iUnfixed = 0;
  for (i = 0; i < getNumFaces(); i++) {
    Face *pF = getFace(i);
    Cell *pCL = pF->getLeftCell();
    Cell *pCR = pF->getRightCell();
    Vert *pV0 = pF->getVert(0);
    Vert *pV1 = pF->getVert(1);
    Vert *pV2 = pF->getVert(2);
    Vert *pVL = pVInvalidVert;
    Vert *pVR = pVInvalidVert;
    int iStatL, iStatR, iStatus;

    // Check orientation on the left.
    switch (pCL->getType())
      {
      case Cell::eTet:
	pVL = dynamic_cast<TetCell*>(pCL)->getOppositeVert(pF);
	iStatL = checkOrient3D(pV0, pV2, pV1, pVL);
	break;
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	iStatL = 2;
	break;
      default:
	iStatL = 3;
	break;
      }

    // Check orientation on the right.
    switch (pCR->getType())
      {
      case Cell::eTet:
	pVR = dynamic_cast<TetCell*>(pCR)->getOppositeVert(pF);
	iStatR = checkOrient3D(pV0, pV1, pV2, pVR);
	break;
      case Cell::eTriBFace:
      case Cell::eQuadBFace:
      case Cell::eIntTriBFace:
      case Cell::eIntQuadBFace:
	iStatR = 2;
	break;
      default:
	iStatR = 3;
	break;
      }

    assert(!pF->hasVert(pVR));
    assert(!pF->hasVert(pVL));

    iStatus = iStatL * 5 + iStatR;

    switch (iStatus)
      {
      case eGoodGood:
      case eBdryGood:
      case eGoodBdry:
      case eGoodNontet:
      case eNontetGood:
      case eNontetBdry:
      case eBdryNontet:
      case eNontetNontet:
	iGood++;
	break;
      case eInvInv:
      case eInvBdry:
      case eBdryInv:
      case eInvNontet:
      case eNontetInv:
	iBad++;
	if (qFix) {
	  pF->interchangeCells();
	  assert(!pF->getRightCell()->hasVert(pVR));
	  assert(!pF->getLeftCell()->hasVert(pVL));
	  iFixed++;
	}
	else
	  iUnfixed++;
	break;
      case eFlatInv:
	if (qFix) {
	  pF->interchangeCells();
	  Vert *pVTemp = pVL;
	  pVL = pVR;
	  pVR = pVTemp;
	  assert(pF->getRightCell()->hasVert(pVR));
	  assert(pF->getLeftCell()->hasVert(pVL));
	}
	iUnfixed++;
	iBad++;
	break;
      case eFlatGood:
	iUnfixed++;
	iBad++;
	break;
      case eInvFlat:
	if (qFix) {
	  pF->interchangeCells();
	  Vert *pVTemp = pVL;
	  pVL = pVR;
	  pVR = pVTemp;
	  assert(pF->getRightCell()->hasVert(pVR));
	  assert(pF->getLeftCell()->hasVert(pVL));
	}
	iUnfixed++;
	iBad++;
	break;
      case eGoodFlat:
	iUnfixed++;
	iBad++;
	break;
      case eFlatFlat:
      case eFlatBdry:
      case eBdryFlat:
      case eNontetFlat:
      case eFlatNontet:
	iUnfixed++;
	iBad++;
	break;
      case eInvGood:
	if (qFix) {
	  pF->interchangeCells();
	  Vert *pVTemp = pVL;
	  pVL = pVR;
	  pVR = pVTemp;
	  assert(pF->getRightCell()->hasVert(pVR));
	  assert(pF->getLeftCell()->hasVert(pVL));
	}
	iUnfixed++;
	iBad++;
	break;
      case eGoodInv:
	iUnfixed++;
	iBad++;
	break;
      case eBdryBdry:
	// This is a really weird case that should never come up.
	// The only fix is to remove the face and the corresponding bdry
	// faces.
	iUnfixed++;
	iBad++;
	break;
      default:
	// Should never be able to get here at all.
	assert(0);
	break;
      }
  }
  assert(iGood + iBad == getNumFaces());
  assert(iBad == iFixed + iUnfixed);
  logMessage(3, "Checked %u faces in mesh; %s\n", getNumFaces(),
	     "bounding tets classified as follows:\n");
  logMessage(3, "  Correct for a right-handed mesh:   %9u\n", iGood);
  logMessage(3, "  Incorrect for a right-handed mesh: %9u\n", iBad);
  if (iBad != 0) {
    if (qFix) {
      logMessage(3, "    Of these, %5u were corrected.\n", iFixed);
      logMessage(3, "              %5u could not be corrected.\n", iUnfixed);
    }
    else {
      logMessage(3, "No attempt was made to correct these errors.\n");
    }
  }
  if (iUnfixed != 0)
    vFatalError("Mesh can not be easily made into a right-handed mesh",
		"processing volume mesh data file");
  else {
    // Make sure nothing got screwed up.
    logMessage(3, "Volume mesh is properly oriented.\n");
  }
  return (isValid());
}

void
VolMesh::classifyTets() const
// Classify tets according to the scheme of Bern, Chew, Eppstein,
// and Ruppert.
{
  static const double dLilDihed = 20., dBigDihed = 150.;
  static const double dLilSolid = 3., dBigSolid = 240.;
  typedef enum {
    eNonTet = 1,
    eRound = 2,
    eNeedle = 4,
    eWedge = 8,
    eSpindle = 16,
    eSliver = 32,
    eCap = 64
  } eTetType;

  GR_index_t iRound = 0, iNeedle = 0, iWedge = 0, iSpindle = 0, iSliver = 0,
      iCap = 0;

  GR_index_t i;
  int *aiFaceType = new int[getNumFaces()];
  for (i = 0; i < getNumFaces(); i++)
    aiFaceType[i] = 0;

  logMessage(
      2, "\nClassifying tets according to which pathologies they have...\n");
  for (int iC = getNumCells() - 1; iC >= 0; iC--) {
    Cell *pC = getCell(iC);
    // Don't do anything for cells which have been deleted or are
    // otherwise bogus.
    if (pC->doFullCheck() != 1)
      continue;
    if (pC->getType() == Cell::eTet) {
      double adDihed[6];
      int iNDihed;
      dynamic_cast<TetCell*>(pC)->calcAllDihed(adDihed, &iNDihed);

      int iSmall = 0, iLarge = 0;
      for (int ii = 0; ii < iNDihed; ii++) {
	if (adDihed[ii] < dLilDihed)
	  iSmall++;
	else if (adDihed[ii] > dBigDihed)
	  iLarge++;
      }

      double dS1 = adDihed[0] + adDihed[1] + adDihed[2] - 180;
      double dS2 = adDihed[0] + adDihed[3] + adDihed[4] - 180;
      double dS3 = adDihed[1] + adDihed[3] + adDihed[5] - 180;
      double dS4 = adDihed[2] + adDihed[4] + adDihed[5] - 180;
      int iSmallSolid = ((dS1 < dLilSolid ? 1 : 0) + (dS2 < dLilSolid ? 1 : 0)
	  + (dS3 < dLilSolid ? 1 : 0) + (dS4 < dLilSolid ? 1 : 0));
      int iLargeSolid = ((dS1 > dBigSolid ? 1 : 0) + (dS2 > dBigSolid ? 1 : 0)
	  + (dS3 > dBigSolid ? 1 : 0) + (dS4 > dBigSolid ? 1 : 0));
      assert(iLarge <= 3);
      assert(iSmall <= 4);
      assert(iLargeSolid <= 1);

      int iThisType;
      if (iLargeSolid == 1) {
	iThisType = eCap;
	iCap++;
      }
      else if (iLarge > 0 && iSmall > 0) {
	iThisType = eSliver;
	iSliver++;
      }
      else if (iLarge > 0) {
	iThisType = eSpindle;
	iSpindle++;
      }
      else if (iSmall > 0) {
	iThisType = eWedge;
	iWedge++;
      }
      else if (iSmallSolid > 0) {
	iThisType = eNeedle;
	iNeedle++;
      }
      else {
	iThisType = eRound;
	iRound++;
      }

      double adRes[1];
      int iNumRes;
      calcShortestToRadius(pC, &iNumRes, adRes);
      assert(iNumRes == 1);
      logMessage(4, "[DETAIL] Tet: %7d  Type: %1d  B: %.5f\n", iC, iThisType,
		 adRes[0]);

      for (int ii = 0; ii < 4; ii++) {
	int iF = getFaceIndex(pC->getFace(ii));
	aiFaceType[iF] += iThisType;
      }
    } // if (tet)
    else {
      for (int ii = 0; ii < pC->getNumFaces(); ii++) {
	int iF = getFaceIndex(pC->getFace(ii));
	aiFaceType[iF] += eNonTet;
      }
    }
  } // End of loop over cells

    // Loop over the boundary faces and make due note of them
  for (i = 0; i < getNumBdryFaces(); i++) {
    BFace *pBF = getBFace(i);
    if (!pBF->isDeleted()) {
      int iF = getFaceIndex(pBF->getFace());
      aiFaceType[iF] += eNonTet;
    }
  }

  GR_index_t iTotal = iRound + iNeedle + iWedge + iSpindle + iSliver + iCap;
  logMessage(2, "Analyzed a total of %u tetrahedra.\n", iTotal);
  logMessage(2, "Number of well-shaped tets:%8u (%5.2f%%)\n", iRound,
	     100 * double(iRound) / double(iTotal));
  logMessage(2, "                   needles:%8u (%5.2f%%)\n", iNeedle,
	     100 * double(iNeedle) / double(iTotal));
  logMessage(2, "                    wedges:%8u (%5.2f%%)\n", iWedge,
	     100 * double(iWedge) / double(iTotal));
  logMessage(2, "                  spindles:%8u (%5.2f%%)\n", iSpindle,
	     100 * double(iSpindle) / double(iTotal));
  logMessage(2, "                   slivers:%8u (%5.2f%%)\n", iSliver,
	     100 * double(iSliver) / double(iTotal));
  logMessage(2, "                      caps:%8u (%5.2f%%)\n", iCap,
	     100 * double(iCap) / double(iTotal));

  // Now tabulate, for each face, the type of its bounding cells
  int a2iFaces[7][7];
  for (i = 0; i < 7; i++)
    for (int ii = 0; ii < 7; ii++)
      a2iFaces[i][ii] = 0;

  for (i = 0; i < getNumFaces(); i++) {
    if (getFace(i)->isDeleted())
      continue;
    switch (aiFaceType[i])
      {
      case eNonTet + eNonTet: // Neither cell is a tet of any kind (bdry,
	// pyr, prism, hex)
	a2iFaces[0][0] += 1;
	break;
      case eNonTet + eRound: // One round, one not a tet
	a2iFaces[0][1] += 1;
	break;
      case eNonTet + eNeedle: // One needle, one not a tet
	a2iFaces[0][2] += 1;
	break;
      case eNonTet + eWedge: // One wedge, one not a tet
	a2iFaces[0][3] += 1;
	break;
      case eNonTet + eSpindle: // One spindle, one not a tet
	a2iFaces[0][4] += 1;
	break;
      case eNonTet + eSliver: // One sliver, one not a tet
	a2iFaces[0][5] += 1;
	break;
      case eNonTet + eCap: // One cap, one not a tet
	a2iFaces[0][6] += 1;
	break;
      case eRound + eRound: // Both round
	a2iFaces[1][1] += 1;
	break;
      case eRound + eNeedle: // One needle, one round
	a2iFaces[1][2] += 1;
	break;
      case eRound + eWedge: // One wedge, one round
	a2iFaces[1][3] += 1;
	break;
      case eRound + eSpindle: // One spindle, one round
	a2iFaces[1][4] += 1;
	break;
      case eRound + eSliver: // One sliver, one round
	a2iFaces[1][5] += 1;
	break;
      case eRound + eCap: // One cap, one round
	a2iFaces[1][6] += 1;
	break;
      case eNeedle + eNeedle: // Both needles
	a2iFaces[2][2] += 1;
	break;
      case eNeedle + eWedge: // One wedge, one needle
	a2iFaces[2][3] += 1;
	break;
      case eNeedle + eSpindle: // One spindle, one needle
	a2iFaces[2][4] += 1;
	break;
      case eNeedle + eSliver: // One sliver, one needle
	a2iFaces[2][5] += 1;
	break;
      case eNeedle + eCap: // One cap, one needle
	a2iFaces[2][6] += 1;
	break;
      case eWedge + eWedge: // Both wedges
	a2iFaces[3][3] += 1;
	break;
      case eWedge + eSpindle: // One spindle, one wedge
	a2iFaces[3][4] += 1;
	break;
      case eWedge + eSliver: // One sliver, one wedge
	a2iFaces[3][5] += 1;
	break;
      case eWedge + eCap: // One cap, one wedge
	a2iFaces[3][6] += 1;
	break;
      case eSpindle + eSpindle: // Both spindles
	a2iFaces[4][4] += 1;
	break;
      case eSpindle + eSliver: // One sliver, one spindle
	a2iFaces[4][5] += 1;
	break;
      case eSpindle + eCap: // One cap, one spindle
	a2iFaces[4][6] += 1;
	break;
      case eSliver + eSliver: // Both slivers
	a2iFaces[5][5] += 1;
	break;
      case eSliver + eCap: // One cap, one sliver
	a2iFaces[5][6] += 1;
	break;
      case eCap + eCap: // Both caps
	a2iFaces[6][6] += 1;
	break;
      default: // Should never get here
	assert(0);
	break;
      }
  }

  logMessage(2,
	     "\nBreakdown of faces by the types of cells on either side\n\n");
  logMessage(
      2, "         Non-tet   Round  Needle   Wedge Spindle  Sliver     Cap\n");
  logMessage(2, "Non-tet %8d %7d %7d %7d %7d %7d %7d\n", a2iFaces[0][0],
	     a2iFaces[0][1], a2iFaces[0][2], a2iFaces[0][3], a2iFaces[0][4],
	     a2iFaces[0][5], a2iFaces[0][6]);
  logMessage(2, "Round   %16d %7d %7d %7d %7d %7d\n", a2iFaces[1][1],
	     a2iFaces[1][2], a2iFaces[1][3], a2iFaces[1][4], a2iFaces[1][5],
	     a2iFaces[1][6]);
  logMessage(2, "Needle  %24d %7d %7d %7d %7d\n", a2iFaces[2][2],
	     a2iFaces[2][3], a2iFaces[2][4], a2iFaces[2][5], a2iFaces[2][6]);
  logMessage(2, "Wedge   %32d %7d %7d %7d\n", a2iFaces[3][3], a2iFaces[3][4],
	     a2iFaces[3][5], a2iFaces[3][6]);
  logMessage(2, "Spindle %40d %7d %7d\n", a2iFaces[4][4], a2iFaces[4][5],
	     a2iFaces[4][6]);
  logMessage(2, "Sliver  %48d %7d\n", a2iFaces[5][5], a2iFaces[5][6]);
  logMessage(2, "Cap     %56d\n\n", a2iFaces[6][6]);

  delete[] aiFaceType;
}

void
VolMesh::initMeshInBrick(const double dXMin, const double dXMax,
			 const double dYMin, const double dYMax,
			 const double dZMin, const double dZMax)
{
  // x runs left to right
  // y runs bottom to top
  // z runs back to front
  Vert *pVLoLeftBack = createVert(dXMin, dYMin, dZMin);
  Vert *pVLoLeftFront = createVert(dXMin, dYMin, dZMax);
  Vert *pVHiLeftBack = createVert(dXMin, dYMax, dZMin);
  Vert *pVHiLeftFront = createVert(dXMin, dYMax, dZMax);
  Vert *pVLoRightBack = createVert(dXMax, dYMin, dZMin);
  Vert *pVLoRightFront = createVert(dXMax, dYMin, dZMax);
  Vert *pVHiRightBack = createVert(dXMax, dYMax, dZMin);
  Vert *pVHiRightFront = createVert(dXMax, dYMax, dZMax);

  pVLoLeftBack->setType(Vert::eBBox);
  pVLoLeftFront->setType(Vert::eBBox);
  pVHiLeftFront->setType(Vert::eBBox);
  pVHiLeftBack->setType(Vert::eBBox);
  pVLoRightBack->setType(Vert::eBBox);
  pVLoRightFront->setType(Vert::eBBox);
  pVHiRightFront->setType(Vert::eBBox);
  pVHiRightBack->setType(Vert::eBBox);

  // Set up the cells
  bool qExist;
  createTetCell(qExist, pVLoLeftBack, pVLoLeftFront, pVLoRightFront,
		pVHiLeftFront);
  assert(!qExist);
  createTetCell(qExist, pVLoRightFront, pVLoRightBack, pVLoLeftBack,
		pVHiRightBack);
  assert(!qExist);
  createTetCell(qExist, pVHiRightBack, pVHiLeftFront, pVHiLeftBack,
		pVLoLeftBack);
  assert(!qExist);
  createTetCell(qExist, pVHiLeftFront, pVHiRightBack, pVHiRightFront,
		pVLoRightFront);
  assert(!qExist);
  createTetCell(qExist, pVHiLeftFront, pVHiRightBack, pVLoRightFront,
		pVLoLeftBack);
  assert(!qExist);

  assert(getNumCells() == 5);
  assert(getNumFaces() == 16);

  // Set up the bdry faces (twelve of these, three from each of the
  // first four tets)
  createBFace(findCommonFace(pVLoLeftBack, pVLoLeftFront, pVLoRightFront));
  createBFace(findCommonFace(pVLoLeftBack, pVLoLeftFront, pVHiLeftFront));
  createBFace(findCommonFace(pVLoLeftFront, pVLoRightFront, pVHiLeftFront));

  createBFace(findCommonFace(pVLoRightFront, pVLoRightBack, pVLoLeftBack));
  createBFace(findCommonFace(pVLoRightFront, pVLoRightBack, pVHiRightBack));
  createBFace(findCommonFace(pVLoRightBack, pVLoLeftBack, pVHiRightBack));

  createBFace(findCommonFace(pVHiRightBack, pVHiLeftFront, pVHiLeftBack));
  createBFace(findCommonFace(pVHiRightBack, pVHiLeftBack, pVLoLeftBack));
  createBFace(findCommonFace(pVHiLeftFront, pVHiLeftBack, pVLoLeftBack));

  createBFace(findCommonFace(pVHiLeftFront, pVHiRightBack, pVHiRightFront));
  createBFace(findCommonFace(pVHiLeftFront, pVHiRightFront, pVLoRightFront));
  createBFace(findCommonFace(pVHiRightBack, pVHiRightFront, pVLoRightFront));

  assert(getNumBdryFaces() == 12);

  // The BC you get, by default, is iDefaultBC.  That's perfectly fine here.

  setAllHintFaces();
  setVertFaceNeighbors();

  assert(isValid());
}

VolMesh::VolMesh(const CubitBox& box, const int iQualMeas,
		 const double dMaxAngle) :
    Mesh(), m_OKToSwapSurfaceEdges(true), m_doEdgeSwaps(true), m_strictPatchChecking(
	true), m_maxAngleForSurfSwap(dMaxAngle), m_ECTriF(), m_ECQuadF(), m_ECTriBF(), m_ECQuadBF(), m_ECTet(), m_ECPyr(), m_ECPrism(), m_ECHex(), m_ECIntTriBF(), m_ECIntQuadBF(), m_Bdry3D(), m_bdryFromMesh(
	true), m_SubsegMap(new SubsegMap())
{
  m_encroachType = eBall;

  logMessage(2, "Generating the volume mesh of the model's bounding box.\n");
  m_qual = new Quality(this, iQualMeas);

  CubitVector min_vec = box.minimum();
  CubitVector max_vec = box.maximum();
  double x_min = min_vec.x();
  double x_max = max_vec.x();
  double y_min = min_vec.y();
  double y_max = max_vec.y();
  double z_min = min_vec.z();
  double z_max = max_vec.z();

  double x_mean = 0.5 * (x_min + x_max);
  double y_mean = 0.5 * (y_min + y_max);
  double z_mean = 0.5 * (z_min + z_max);
  double dx = x_max - x_min;
  double dy = y_max - y_min;
  double dz = z_max - z_min;

  // Make the region roughly cubical
  double max_aspect = 0.1;
  if (dx < max_aspect * dy)
    dx = max_aspect * dy;
  if (dx < max_aspect * dz)
    dx = max_aspect * dz;

  if (dy < max_aspect * dx)
    dy = max_aspect * dx;
  if (dy < max_aspect * dz)
    dy = max_aspect * dz;

  if (dz < max_aspect * dx)
    dz = max_aspect * dx;
  if (dz < max_aspect * dy)
    dz = max_aspect * dy;

  double x_lo = x_mean - dx;
  double x_hi = x_mean + dx;
  double y_lo = y_mean - dy;
  double y_hi = y_mean + dy;
  double z_lo = z_mean - dz;
  double z_hi = z_mean + dz;

  initMeshInBrick(x_lo, x_hi, y_lo, y_hi, z_lo, z_hi);

  // Once upon a time, we swapped at this point.  But why bother, with
  // mesh-in-a box?  All the points are co-spherical, so anything's
  // Delaunay.

  if (!isValid())
    vFoundBug("creation of a tetrahedral mesh of model's bounding box.");
  logMessage(2, "Done creating volume mesh.\n");

}

VolMesh::VolMesh(const EntContainer<Vert>& EC, const int iQualMeas,
		 const double dMaxAngle, const CubitBox* const BBox) :
    Mesh(EC), m_OKToSwapSurfaceEdges(true), m_doEdgeSwaps(true), m_strictPatchChecking(
	true), m_maxAngleForSurfSwap(dMaxAngle), m_ECTriF(), m_ECQuadF(), m_ECTriBF(), m_ECQuadBF(), m_ECTet(), m_ECPyr(), m_ECPrism(), m_ECHex(), m_ECIntTriBF(), m_ECIntQuadBF(), m_Bdry3D(), m_bdryFromMesh(
	true), m_SubsegMap(new SubsegMap())
// Create a volume mesh of a point cloud, placed inside a big box.
{
  logMessage(2, "Generating volume mesh from point cloud.\n");
  m_qual = new Quality(this, iQualMeas);

  // Set up an outer triangulation of a big brick
  double dXMin, dXMax, dYMin, dYMax, dZMin, dZMax;
  dXMin = dYMin = dZMin = 1.e100;
  dXMax = dYMax = dZMax = -1.e100;

  if (BBox) {
    CubitVector min_vec = BBox->minimum();
    CubitVector max_vec = BBox->maximum();
    dXMin = min_vec.x();
    dXMax = max_vec.x();
    dYMin = min_vec.y();
    dYMax = max_vec.y();
    dZMin = min_vec.z();
    dZMax = max_vec.z();
  }
  else {
    int iVert, iNVerts = EC.lastEntry();
    for (iVert = 0; iVert < iNVerts; iVert++) {
      double dX = EC[iVert].x();
      double dY = EC[iVert].y();
      double dZ = EC[iVert].z();
      dXMin = min(dXMin, dX);
      dXMax = max(dXMax, dX);
      dYMin = min(dYMin, dY);
      dYMax = max(dYMax, dY);
      dZMin = min(dZMin, dZ);
      dZMax = max(dZMax, dZ);
    }
  }

  double dXMean = 0.5 * (dXMin + dXMax);
  double dYMean = 0.5 * (dYMin + dYMax);
  double dZMean = 0.5 * (dZMin + dZMax);
  double dDX = dXMax - dXMin;
  double dDY = dYMax - dYMin;
  double dDZ = dZMax - dZMin;

  // Make the region roughly cubical
  double dMaxAspect = 0.1;
  if (dDX < dMaxAspect * dDY)
    dDX = dMaxAspect * dDY;
  if (dDX < dMaxAspect * dDZ)
    dDX = dMaxAspect * dDZ;

  if (dDY < dMaxAspect * dDX)
    dDY = dMaxAspect * dDX;
  if (dDY < dMaxAspect * dDZ)
    dDY = dMaxAspect * dDZ;

  if (dDZ < dMaxAspect * dDX)
    dDZ = dMaxAspect * dDX;
  if (dDZ < dMaxAspect * dDY)
    dDZ = dMaxAspect * dDY;

  double dXLo = dXMean - dDX;
  double dXHi = dXMean + dDX;
  double dYLo = dYMean - dDY;
  double dYHi = dYMean + dDY;
  double dZLo = dZMean - dDZ;
  double dZHi = dZMean + dDZ;

  initMeshInBrick(dXLo, dXHi, dYLo, dYHi, dZLo, dZHi);

  allowSwapRecursion();
  allowSurfaceEdgeChanges();
  allowEdgeSwapping();
  setSwapType(eDelaunay);

  //@@ Insert all points from the point set
  // Currently these are inserted lexicographically and swapping is
  // performed after each insertion.  Don't insert the verts forming the
  // big box.
  logMessage(1, "Adding remaining vertices to the mesh...\n");
  GR_index_t iSwaps = 0;
  GR_index_t iTenth = 1 + getNumVerts() / 10;
  if (iTenth < 25)
    iTenth = 25;
  GR_index_t iV;
  GRUMMP::DelaunaySwapDecider3D SDDel3D(true);
  GRUMMP::SwapManager3D SMDel3D(&SDDel3D, this);
  GRUMMP::SwappingInserter3D SI3D(this, &SMDel3D);
  SI3D.setForcedInsertion(true);
  for (iV = 0; iV < getNumVerts() - 8; iV++) {
    // Don't insert vertices that are already in the mesh
    Vert *pVInsert = getVert(iV);

    GR_index_t iC = getNumCells() - 1;
    Cell *pCGuess;
    while ((pCGuess = getCell(iC))->doFullCheck() != 1)
      iC--;
    int iSwapsThisTime = 0;

    // Allow swapping after insertion; force insertion to take place
    // without snap-to-edge/face.  NOTE:  This vertex already exists in
    // the point set.
    SI3D.insertPoint(pVInsert->getCoords(), pCGuess, pVInsert);
    iSwaps += iSwapsThisTime;
    if ((iV + 1) % iTenth == 0)
      logMessage(1, "%5u %s  %6u swaps performed.\n", iV,
		 "points inserted from the point cloud.", iSwaps);
  }

  logMessage(1, "%5u %s  %6u swaps performed.\n", iV,
	     "points inserted from the point cloud.", iSwaps);

  logMessage(1, "All verts added; improving connectivity.\n");

  // That's it.  The entire mesh is relevant and no exterior
  // cleanup is required.

  // Once upon a time, swapping was done at this point, but the mesh was
  // swapped during insertion, so why bother?  If the caller wants
  // something other than Delaunay, they can do it themselves.

  if (!isValid())
    vFoundBug("creation of constrained 3D mesh from point cloud");
  logMessage(2, "Done creating constrained volume mesh.\n");
}

void
VolMesh::createPatchesFromBdryFaces(const int aiBCList[])
// The array aiBCList tells what the bdry conditions are for each
// BFace.  This data makes it trivial to assemble the Bdry3D object.
{
  if (m_Bdry3D)
    delete m_Bdry3D;
  m_bdryFromMesh = true;
  m_Bdry3D = new Bdry3D(m_ECVerts);

  for (GR_index_t iBF = 0; iBF < getNumTotalBdryFaces(); iBF++) {
    Vert *apV[3];

    BFace *pBF = getBFace(iBF);
    if (pBF->isDeleted())
      continue;
    Face *pF = pBF->getFace();
    Cell *pC0 = pF->getLeftCell();
    Cell *pC1 = pF->getRightCell();

    // These verts actually need to refer to the vert array
    // in the Bdry 3D object.
    apV[0] = m_Bdry3D->pVVert(getVertIndex(pF->getVert(0)));
    apV[1] = m_Bdry3D->pVVert(getVertIndex(pF->getVert(1)));
    apV[2] = m_Bdry3D->pVVert(getVertIndex(pF->getVert(2)));

    int iBCL = iInvalidBC, iBCR = iInvalidBC;
    int iRegR = iInvalidRegion, iRegL = iInvalidRegion;
    if (pBF->getNumFaces() == 1) {
      if (pC0 == pBF) {
	// Boundary on the left!
	iBCL = aiBCList[iBF];
	iRegR = pC1->getRegion();
      }
      else {
	assert(pC1 == pBF);
	// Boundary on the right!
	iBCR = aiBCList[iBF];
	iRegL = pC0->getRegion();
      }
    }
    else {
      iRegL = pBF->getFace(0)->getOppositeCell(pBF)->getRegion();
      iRegR = pBF->getFace(1)->getOppositeCell(pBF)->getRegion();
      iBCL = aiBCList[iBF];
      iBCR = aiBCList[iBF];
    }
    BdryPolygon *pBPoly = new BdryPolygon(iBCL, iBCR, iRegL, iRegR, 3, apV);
    m_Bdry3D->vAddPatch(pBPoly);
    pBF->setPatch(pBPoly);
  } // Done with loop over all BFaces to create Bdry3D object.
  assert(m_Bdry3D->iNumBdryPatches() <= this->getNumTotalBdryFaces());
}

void
VolMesh::identifyVertexTypesGeometrically() const
{
  // TODO:  It's a hack that I'm preserving shock markings
  // here.  There has to be a better general mechanism here.
  // ITAPS-style tags, maybe?
  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
    if (getVert(iV)->getVertType() < Vert::eShockUpstream) {
      getVert(iV)->setType(Vert::eInterior);
    }
  }

  // Identify and count all boundary vertices and make a list of all the
  // normals for boundary faces incident on a given vertex.  Once this
  // list is sorted, it will be used to divide boundary verts into
  // eBdryApex, eBdryCurve, and eBdry verts.

  GR_index_t iBF;
  assert(isSimplicial());
  BdryNorm *aBN = new BdryNorm[getNumBdryFaces() * 3 + getNumIntBdryFaces() * 3];
  for (iBF = 0; iBF < getNumBdryFaces(); iBF++) {
    BFace *pBF = getBFace(iBF);
    double adNorm[3];
    pBF->calcUnitNormal(adNorm);
    for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
      aBN[3 * iBF + iV] = BdryNorm(pBF->getVert(iV), adNorm);
    }
  }
  for (iBF = 0; iBF < getNumIntBdryFaces(); iBF++) {
    BFace *pBF = getIntBFace(iBF);
    if (pBF->isDeleted())
      continue;
    double adNorm[3];
    pBF->calcUnitNormal(adNorm);
    for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
      aBN[3 * getNumBdryFaces() + 3 * iBF + iV] = BdryNorm(pBF->getVert(iV),
							   adNorm);
    }
  }
  // For each run of identical vertices, identify the number of
  // distinct normals and tag the vertex accordingly:
  //  1 -> eBdry
  //  2 -> eBdryCurve
  //  3 or more -> eBdryApex
  // Mark boundary apexes to NOT be deleted; that is, clear
  // DELETE_REQUEST for these vertices.

  sort(aBN, aBN + 3 * getNumBdryFaces() + 3 * getNumIntBdryFaces());
  GR_index_t iFirst, iLast;
  GR_index_t iBdry = 0, iBdryCurve = 0, iBdryApex = 0, iTotal = 0;
  for (iFirst = 0; iFirst < 3 * getNumBdryFaces() + 3 * getNumIntBdryFaces();
      iFirst = iLast + 1) {
    Vert *pV = aBN[iFirst].pVVert();
    if (!pV->isValid() || pV->getVertType() >= Vert::eShockUpstream) {
      iLast = iFirst;
      continue;
    }
    iTotal++;

    // Make sure this vert hasn't been touched before.
    assert(pV->getVertType() == Vert::eInterior);

    // Identify the first and last record with this vertex in it.
    iLast = iFirst;
    while (iLast < 3 * getNumBdryFaces() + 3 * getNumIntBdryFaces()
	&& aBN[iLast].pVVert() == pV) {
      iLast++;
    };
    iLast--;
    assert(aBN[iLast].pVVert() == pV);

    // For each boundary vertex, find the number of distinct normals
    // for bdry faces incident on it.  If there are more than two,
    // this is a eBdryApex.  If there are zero, this is just an
    // ordinary eBdry vert.  Otherwise, there should be exactly two.
    //   In this case, find the edges that divide these normals into
    // two groups.  If the edges are roughly collinear, then the
    // vertex is on a eBdryCurve.  If the edges aren't collinear, then
    // the vertex is a eBdryApex.  This is a knife case, which may be
    // uncommon in practice but should still be handled correctly.
    GR_index_t iNumNorms = 0;
    int aiNormTag[1000]; // Should be far fewer bdry faces incident on
    // a single bdry vert.
    GR_index_t iNorm, iTag = 0;

    assert(iFirst != iLast);
    GR_index_t iDist = iLast - iFirst;

    logMessage(4, "Vertex %u has %u incident bdry faces.\n", getVertIndex(pV),
	       iDist + 1);

    for (iNorm = 0; iNorm <= iDist; iNorm++)
      aiNormTag[iNorm] = -1;

    for (iNorm = 0; iNorm <= iDist; iNorm++) {
      if (aiNormTag[iNorm] == -1) {
	aiNormTag[iNorm] = iTag;

	const double *adNormN = aBN[iFirst + iNorm].adNormal();
	double dMagN = dMAG3D(adNormN);
	//				printf("norm %.15f %.15f %.15f\n",adNormN[0],adNormN[1],adNormN[2]);

	for (GR_index_t iiF = iNorm + 1; iiF <= iDist; iiF++) {
	  if (aiNormTag[iiF] == -1) {
	    const double *adNormF = aBN[iFirst + iiF].adNormal();
	    //						printf("inner norm %.15f %.15f %.15f\n",adNormF[0],adNormF[1],adNormF[2]);
	    double dMagF = dMAG3D(adNormF);
	    double dDot = dDOT3D(adNormF, adNormN) / dMagF / dMagN;
	    //						printf("dDot %.15f %.15f\n",dDot, (cos(m_maxAngleForSurfSwap * M_PI/180.)));
	    // All normals are directed into the domain, so there's no
	    // need to take the absolute value of dDot.  And taking the
	    // absolute value can cause faces that have normals nearly
	    // opposite each other to test wrong here.
	    if (dDot >= fabs(cos(m_maxAngleForSurfSwap * M_PI / 180.))) {
	      aiNormTag[iiF] = iTag;
	    }
	  }
	} // Loop over rest of bdry faces for this vert
	iTag++;
      } // If this bdry faces is untagged, process it
    } // Check all bdry faces incident on the vertex.
    iNumNorms = iTag;
    assert(iNumNorms >= 1);
    // Tag the vertex appropriately.
    // Fix for 0.2.1:  There is a knife case in which there are only
    // two distinct normals at a eBdryApex.
    //		printf("iNumNorms %d\n",iNumNorms);
    switch (iNumNorms)
      {
      case 1:
	pV->setType(Vert::eBdry);
	iBdry++;
	break;
      case 2:
	pV->setType(Vert::eBdryCurve);
	iBdryCurve++;
	break;
      default:
	pV->setType(Vert::eBdryApex);
	iBdryApex++;
	pV->markToKeep();
	break;
      } // end switch
  }
  //	 Lets check bdry verts for internal boundaries.
  for (iBF = 0; iBF < getNumIntBdryFaces(); iBF++) {
    BFace *pBF = getIntBFace(iBF);
    if (pBF->isDeleted())
      continue;
    for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
      if (pBF->getVert(iV)->getVertType() == Vert::eBdry)
	pBF->getVert(iV)->setType(Vert::eBdryTwoSide);
    }
  }
  logMessage(2, "Out of %u verts on the bdry, %u were apexes, %u on curves.\n",
	     iTotal, iBdryApex, iBdryCurve);
  assert(iTotal == iBdry + iBdryApex + iBdryCurve);
  delete[] aBN;
}

bool
VolMesh::areBdryFacesOKToSwap(const TriBFaceBase * const pBF0,
			      const TriBFaceBase * const pBF1) const
{
  if (m_strictPatchChecking) {
    return (pBF0->getPatchPointer() == pBF1->getPatchPointer());
  }
  else {
    // Check that the BC's are the same and that the normals are more or
    // less aligned (within the tolerances of dMaxAngleForSurfSwap).
    if (pBF0->getBdryCondition() == pBF1->getBdryCondition()) {
      double adNorm0[3], adNorm1[3];
      pBF0->calcUnitNormal(adNorm0);
      pBF1->calcUnitNormal(adNorm1);
      double dDot = dDOT3D(adNorm0, adNorm1);
      return (dDot >= cos(m_maxAngleForSurfSwap * M_PI / 180));
    }
    else {
      return false;
    }
  }
}

void
VolMesh::initializeSubsegs(const set<Subseg*>* const given_subsegs)
{
  if (!isSubsegMapEmpty())
    return;
  logMessage(2, "Initializing subsegs\n");

  assert(isSubsegMapEmpty());

  if (given_subsegs) {
    //The case where subsegments are given. This is easy...
    //Note that only the pointer gets added to m_subseg_map (not a copy of the object).
    //The memory gets de-allocated only when the subsegment gets deleted from m_mesh.
    //Be careful with that... (should maybe create a local copy if the subsegs
    //in the map must be used somewhere else, in another object).

    Subseg* subseg;
    set<Subseg*>::const_iterator it = given_subsegs->begin(), it_end =
	given_subsegs->end();

    for (; it != it_end; ++it) {

      subseg = *it;

#ifndef NDEBUG
      assert(subseg->get_beg_vert()->isValid());
      assert(subseg->get_end_vert()->isValid());
      assert(subseg->get_beg_vert() != subseg->get_end_vert());
      assert(!subseg->get_beg_vert()->isDeleted());
      assert(!subseg->get_end_vert()->isDeleted());
      hasVert(subseg->get_beg_vert());
      hasVert(subseg->get_end_vert());
#endif

      addSubseg(subseg);

    }

  }

  else {

    //Subsegments are not given. We must infer them from
    //the mesh and the underlying geometry.

    //This case should only happen with GRUMMP's native
    //geometry interface (and the following code is
    //designed as such). The code does not store geometric
    //entity pointers for vertices (makes it impossible to
    //identify orphaned subsegments inside a boundary surface).

    assert(!given_subsegs);

    //To be a subsegment, an edge must have both its vertices
    //either on a curve, or at an apex point.

    //Go through all the mesh boundary faces and identify all
    //the edges satisfying this condition. Store them in a set
    //as TetRefinerEdge's, this way each edge is stored
    //once and only once.
    BFace* bface;
    Vert* vert[2];

    TetRefinerEdgeBFace edge;
    set<TetRefinerEdgeBFace> all_edges;
    set<TetRefinerEdgeBFace>::iterator it, it_end;

    GR_index_t i, num_bfaces = getNumTotalBdryFaces();
    int j;
    bool qBdryPatchesInferred = isBdryFromMesh();

#ifndef NDEBUG
    int num_bdry_verts = 0;
    for (i = 0; i < getNumVerts(); ++i) {
      Vert* my_vert = getVert(i);
      if (my_vert->isDeleted())
	continue;
      if (my_vert->getVertType() == Vert::eBdryApex
	  || my_vert->getVertType() == Vert::eBdryCurve)
	++num_bdry_verts;
    }
#endif

    for (i = 0; i < num_bfaces; ++i) {

      bface = getBFace(i);
      if (bface->isDeleted())
	continue;

//			assert(
//					bface->getType() == Cell::eTriBFace
//							|| bface->getType() == Cell::eIntTriBFace);
//			assert(bface->getNumVerts() == 3);

      for (j = 0; j < bface->getNumVerts(); ++j) {
	vert[0] = bface->getVert(j);
	vert[1] = bface->getVert((j + 1) % bface->getNumVerts());

	if ((vert[0]->getVertType() == Vert::eBdryApex
	    || vert[0]->getVertType() == Vert::eBdryCurve)
	    && (vert[1]->getVertType() == Vert::eBdryApex
		|| vert[1]->getVertType() == Vert::eBdryCurve)) {

	  edge = TetRefinerEdgeBFace(vert[0], vert[1]);
	  it = all_edges.find(edge);

	  if (it == all_edges.end()) { //edge not in the set.
	    edge.add_bface(bface);
	    all_edges.insert(edge);
	  }
	  else { //edge already in set.
	    edge = *it;
	    assert(!edge.get_bfaces().empty());
	    edge.add_bface(bface);
	    all_edges.erase(it);
	    all_edges.insert(edge);
	  }

	}

	else
	  continue; //definitely not a subseg, move on!

      }

    }

    //Now we have a set of unique edges that are potentially
    //subsegments. Go through this set and identify which edges
    //are truly subsegments.

    Subseg* subseg;
    it = all_edges.begin();
    it_end = all_edges.end();

    for (; it != it_end; ++it) {

      subseg = NULL;
      vert[0] = it->vert(0);
      vert[1] = it->vert(1);

      //This piece of code should only be used when using GRUMMP
      //without an external geometry library (i.e. CGM). As such,
      //there should not be geometry entity stored with the vertices.
      assert(!vert[0]->getParentTopology());
      assert(!vert[1]->getParentTopology());

      switch (it->get_bfaces().size())
	{

	case 0: //These two cases should never happen with GRUMMP's native geometry,
	case 1: //otherwise something is wrong.
	  assert(0);
	  break;

	case 2:
	  {
	    BFace* bface0 = *((it->get_bfaces().begin()));
	    BFace* bface1 = *(++(it->get_bfaces().begin()));
	    assert(bface0 && bface1);
	    assert(bface0 != bface1);
	    assert(!bface0->isDeleted() && !bface1->isDeleted());
	    assert(bface0->getPatchPointer());
	    assert(bface1->getPatchPointer());
	    double dDot = 0.;
	    double norm0[3], norm1[3];
	    bface0->calcUnitNormal(norm0);
	    bface1->calcUnitNormal(norm1);
	    dDot = dDOT3D(norm0, norm1);

	    //				logMessage(2,"qq\n");
	    //				bface0->getFace(0)->printFaceInfo();
	    //				bface1->getFace(0)->printFaceInfo();
	    if ((bface0->getPatchPointer() != bface1->getPatchPointer()
		&& !qBdryPatchesInferred)
		|| (fabs(dDot) < 0.99999999 && qBdryPatchesInferred))
	      subseg = new Subseg(vert[0], vert[1]);
	    break;
	  }

	default: //If there are more than 2 distinct bfaces, we definitely have a subseg.
	  subseg = new Subseg(vert[0], vert[1]);
	  break;

	}

      if (subseg)
	addSubseg(subseg);

    }

    //Must make sure that every vertex on a curve or at an apex
    //has an entry in the subseg map.
    // not necessarily true;
    if (qBdryPatchesInferred) {
      assert(num_bdry_verts == subsegMapSize());
    }
  } //End case where subsegments are not given.

  //Assert that the entries in the subseg map are consistent.
  assert(subsegMapIsValid());
  printSubsegMap();
}

// TO DO: Make these non-member functions

double
VolMesh::calcBoundarySize() const
{
  double dRetVal = 0;
  for (int i = getNumBdryFaces() - 1; i >= 0; i--)
    dRetVal += getBFace(i)->calcSize();
  return dRetVal;
}

double
VolMesh::calcInteriorSize() const
{
  double dRetVal = 0;
  for (int i = getNumCells() - 1; i >= 0; i--)
    dRetVal += getCell(i)->calcSize();
  return dRetVal;
}

bool
VolMesh::isWatertight() const
{
  double adSum[] =
    { 0., 0., 0. };
  double adTemp[3];
  for (int i = getNumBdryFaces() - 1; i >= 0; i--) {
    BFace *pBF = getBFace(i);
    if (pBF->doFullCheck()) {
      pBF->calcVecSize(adTemp);
      adSum[0] += adTemp[0];
      adSum[1] += adTemp[1];
      adSum[2] += adTemp[2];
    }
  }
  double dCompSize = calcBoundarySize();
  double dInvCompSize = 0.001 / dCompSize;
  vSCALE3D(adSum, dInvCompSize);

  if (iFuzzyComp(adSum[0], 0.) == 0 && iFuzzyComp(adSum[1], 0.) == 0
      && iFuzzyComp(adSum[2], 0.) == 0)
    return true;
  else {
    logMessage(2, "Directed surface areas are non-zero:\n");
    logMessage(2, "  x: %12.5g\n  y: %12.5g\n  z: %12.5g\n", adSum[0], adSum[1],
	       adSum[2]);
    return false;
  }
}

Cell*
VolMesh::findCell(const double newPtLoc[], Cell* const cellGuess,
		  bool& status) const
{
  // TODO Does this routine work when there are internal bdrys?  Also, I feel like
  // this sort of code appears a lot more than once; let's clean that up.

  //@@ Find which cell the site falls within (or near, if outside the mesh)
  Cell *cell = cellGuess;
  if (cell->isValid() && !cell->isDeleted()) {
    if (cell->getType() == Cell::eTriBFace
	|| cell->getType() == Cell::eIntTriBFace) {
      Face *pFTmp = cell->getFace(0);
      cell = pFTmp->getOppositeCell(cell);
    }
  }
  if (!cell->isValid() || cell->isDeleted()
      || !(cell->getType() == CellSkel::eTet)) {
    GR_index_t i = 0;
    do {
      cell = getCell(i);
      i++;
    }
    while (cell->isDeleted() && i < getNumCells());
  }
  assert(
      cell->isValid() && !cell->isDeleted()
	  && cell->getType() == CellSkel::eTet);
  // At this point, the cell is definitely a live tet.
  Cell *cellNew = cell;

  int iOrient = 0;
  int iSteps = 0;
  Face *pF;
  do {
    status = true;
    iSteps++;
    int i;
    double adCent[3];
    cell->calcCentroid(adCent);
    for (i = 0; i < 4; i++) {
      pF = cell->getFace(i);
      // If the face is pointing away from this cell, then signs need to
      // be reversed to walk in the proper direction
      int iSign = (pF->getLeftCell() == cell) ? -1 : 1;
      iOrient = checkOrient3D(pF, newPtLoc) * iSign;
      // If the vertex lies outside the current cell and behind the
      // current face, walk in that direction.  If the walk would
      // involve going outside the boundary, continue checking to try to
      // find another face to walk across; this improves behavior of
      // walks which try to insert a boundary point.
      if (iOrient == -1) {
	if (iSteps < 2000) {
	  cellNew = pF->getOppositeCell(cell);
	  status = false;
	  if (cellNew->getType() != Cell::eTriBFace)
	    break;
	}
	else {
	  // Assume that we may be entering orbit...
	  int iOrient01 = checkOrient3D(adCent, pF->getVert(0)->getCoords(),
					pF->getVert(1)->getCoords(), newPtLoc);
	  int iOrient12 = checkOrient3D(adCent, pF->getVert(1)->getCoords(),
					pF->getVert(2)->getCoords(), newPtLoc);
	  int iOrient20 = checkOrient3D(adCent, pF->getVert(2)->getCoords(),
					pF->getVert(0)->getCoords(), newPtLoc);
	  // If the target point lies on or inside the infinite
	  // triangular pyramid defined by the cell centroid and the
	  // vertices of pF, then stepping through this face
	  // definitely moves us in the direction of the target.
	  if (iOrient01 * iOrient12 != -1 && iOrient12 * iOrient20 != -1
	      && iOrient20 * iOrient01 != -1) {
	    cellNew = pF->getOppositeCell(cell);
	    status = false;
	    if (cellNew->getType() != Cell::eTriBFace)
	      break;
	  }
	}
      }
    }
    cell = cellNew;
  }
  while (iOrient == -1 && cell->getType() != Cell::eTriBFace && iSteps < 5000);
  if (cell->getType() == Cell::eTriBFace) {
    // Nothing for it but to try old-fashioned brute-force search.
    status = false;
    GR_index_t iC = 0;
    for (; (iC < getNumTetCells()) && (status == false); iC++) {
      cell = getCell(iC);
      if (cell->isDeleted())
	continue;
      double bary[] =
	{ -1, -1, -1, -1 };
      dynamic_cast<TetCell*>(cell)->calcBarycentricCoords(newPtLoc, bary);
//			logMessage(1, "Bary: (%12f %12f %12f %12f\n",
//					bary[0], bary[1], bary[2], bary[3]);
      status = true;
      for (int ii = 0; ii < 4; ii++) {
	if (iFuzzyComp(bary[ii], 0) == -1) {
	  status = false;
	  break;
	}
      }
    }
    if (status == false)
      return nullptr; // This should cause hysteristics in the caller.
    logMessage(3, "Brute force: %6u\n", iC);
  }
  else {
    logMessage(3, "Walk did it: %6d\n", iSteps);
  }
  assert(cell);
  return cell;
}

void
VolMesh::writeTempMesh()
{
  if (m_writeTempMeshes) {
    char strName[FILE_NAME_LEN];
    sprintf(strName, "%.1000s_%d", m_tempMeshFileName, m_tempMeshIndex);
    writeVTKLegacyWithoutPurge(*this, strName);
    m_tempMeshIndex++;
  }
}
