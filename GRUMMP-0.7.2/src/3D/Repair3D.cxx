#include "GR_VolMesh.h"
#include "GR_Geometry.h"
#include "GR_Vec.h"
#include "GR_BFace.h"
#include "GR_FaceSwapInfo.h"
#include "GR_InsertionManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_SmoothingManager.h"

extern int aiEdgeReq[], aiEdgeDone[];

static double dSmallDihed = 0.;
static double dLargeDihed = 180.;
static double dSmallSolid = 0.;
static double dLargeSolid = 360.;
static const double dSmallDihedMargin = 10.;
static const double dLargeDihedMargin = 20.;
static const double dSmallSolidMargin = 3.;
static const double dLargeSolidMargin = 60.;

// The ends of the edges associated with each dihedral angle.
static const int aiVert0[] =
  { 0, 0, 0, 1, 1, 2 };
static const int aiVert1[] =
  { 1, 2, 3, 2, 3, 3 };
// These other two are used to generate calls to iEdgeSwap.
static const int aiVertOther[] =
  { 2, 3, 1, 0, 2, 0 };
static const int aiVertOpp[] =
  { 3, 1, 2, 3, 0, 1 };

static bool
qIsBdryEdge(const Face * const pFStart, const Vert * const pV0,
	    const Vert * const pV1)
// Answers the question: does one hit a boundary while traversing
// around this edge?
{
  assert(pFStart->hasVert(pV0) && pFStart->hasVert(pV1));
  const Face * pF = pFStart;
  const Cell * pC = pF->getLeftCell();
  const Vert * pV;

  if (pC->getType() == Cell::eTriBFace || pC->getType() == Cell::eIntTriBFace)
    return true;
  assert(pC->getType() == Cell::eTet);

  do {
    pV = pF->getVert(0);
    if (pV == pV0 || pV == pV1) {
      pV = pF->getVert(1);
      if (pV == pV0 || pV == pV1)
	pV = pF->getVert(2);
    }
    assert(pV != pV0 && pV != pV1 && pF->hasVert(pV));

    pC = pF->getOppositeCell(pC);
    if (pC->getType() == Cell::eTriBFace || pC->getType() == Cell::eIntTriBFace)
      return true;
    assert(pC->getType() == Cell::eTet);

    pF = dynamic_cast<const TetCell*>(pC)->getOppositeFace(pV);
  }
  while (pF != pFStart);
  return false;
}

int
VolMesh::splitLargeAngle(Cell *pC, const double dSplitLimit,
			 const double * const adDihedIn,
			 const bool qOnlyInterior)
{
  // For tets with large dihedral angles:
  //   1. Insert a point on the intersection of the angle bisector
  //      of that angle and the opposite edge.
  //   2. Smooth the location of that point.
  //   3. Swap in the vicinity of the point.
  //
  GRUMMP::OptMSSmoothingManager3D Smooth(this);
  GRUMMP::SwappingInserter3D SI3D(this, NULL);
  SI3D.setForcedInsertion(true);

  int i;
  assert(pC->getType() == Cell::eTet);
  if (pC->doFullCheck() != 1)
    return 0;
  double adDihed[6];

  Vert *apV[4];
  for (i = 0; i < 4; i++)
    apV[i] = pC->getVert(i);

  if (adDihedIn == 0) {
    dynamic_cast<TetCell*>(pC)->calcAllDihed(adDihed, &i);
    assert(i == 6);
  }
  else {
    for (i = 0; i < 6; i++) {
      adDihed[i] = adDihedIn[i];
      assert(adDihed[i] >= 0 && adDihed[i] <= 180);
    }
  }

  int iLargest = -1, iLargestInt = -1;
  double dLargest = dSplitLimit;
  double dLargestInt = dSplitLimit;
  // Determine which angles are large and which edges are boundary
  // edges.  Edges which don't fit the definition of "large" are marked
  // as interior edges by default; this does no harm, because we never
  // try to insert a point on them anyway.
//   bool qLarge[] = {false, false, false, false, false, false};
  int iNLarge = 0;
//   int iNBdry = 0;
//   bool qBEdge[] = {false, false, false, false, false, false};
  for (i = 0; i < 6; i++) {
    if (adDihed[i] > dSplitLimit) {
//       qLarge[i] = true;
      iNLarge++;
      if (adDihed[i] > dLargest) {
	dLargest = adDihed[i];
	iLargest = i;
      }
      Face *pF = pC->getOppositeFace(apV[aiVert0[i]]);
      if (qIsBdryEdge(pF, apV[aiVertOpp[i]], apV[aiVertOther[i]])) {
//	qBEdge[i] = true;
//	iNBdry ++;
      }
      else {
	if (adDihed[i] > dLargestInt) {
	  dLargestInt = adDihed[i];
	  iLargestInt = i;
	}
      }
    }
  }
  if (iNLarge == 0)
    return 0;
  assert(iLargest != -1);

  int iAngleToSplit = -1;
  // First choice: pick the edge whose opposite edge is internal that
  // has the largest dihedral angle.  Second choice: pick the edge that
  // has the largest dihedral angle.
  if (iLargestInt != -1)
    iAngleToSplit = iLargestInt;
  else if (qOnlyInterior)
    return 0;
  else
    iAngleToSplit = iLargest;
  assert(iAngleToSplit != -1);

  Vert *pVA = apV[aiVert0[iAngleToSplit]];
  Vert *pVB = apV[aiVert1[iAngleToSplit]];
  Vert *pVC = apV[aiVertOpp[iAngleToSplit]];
  Vert *pVD = apV[aiVertOther[iAngleToSplit]];

  // Don't split the edge if it is already too short.
  {
    double adDiff[] =
      { pVC->x() - pVD->x(), pVC->y() - pVD->y(), pVC->z() - pVD->z() };
    double dLen = dMAG3D(adDiff);
    double dCorrectLen = 0.5 * (pVC->getLengthScale() + pVD->getLengthScale());
    if (dLen < dCorrectLen * 0.8) {
      return 0;
    }
  }

  // Find the normal to the angle bisector plane. (Sign is irrelevant.)
  //  Find the unit normals of the faces on either side
  Face *pFC = pC->getOppositeFace(pVC);
  Face *pFD = pC->getOppositeFace(pVD);
  double adNormA[3], adNormB[3];
  pFC->calcUnitNormal(adNormA);
  pFD->calcUnitNormal(adNormB);

  //  Adjust both vectors to point into the cell.
  if (pFC->getLeftCell() == pC)
    vSCALE3D(adNormA, -1);
  if (pFD->getLeftCell() == pC)
    vSCALE3D(adNormB, -1);

  //  Add to get a vector parallel to the bisector plane.
  double adParallel[] =
	{ adNormA[0] + adNormB[0], adNormA[1] + adNormB[1], adNormA[2]
	    + adNormB[2] };

  //  Cross that vector with the edge vector; this gives the normal.
  double adEdge[3];
  adEdge[0] = (pVA->x() - pVB->x());
  adEdge[1] = (pVA->y() - pVB->y());
  adEdge[2] = (pVA->z() - pVB->z());

  double adBisectorNorm[3];
  vCROSS3D(adParallel, adEdge, adBisectorNorm);
  NORMALIZE3D(adBisectorNorm);

  assert(qFuzzyPerp3D(adBisectorNorm, adEdge));

  // Find the distance from the angle bisector plane of each of
  // the opposite points.
  adEdge[0] = (pVA->x() - pVC->x());
  adEdge[1] = (pVA->y() - pVC->y());
  adEdge[2] = (pVA->z() - pVC->z());
  double dDotC = dDOT3D(adEdge, adBisectorNorm);

  adEdge[0] = (pVA->x() - pVD->x());
  adEdge[1] = (pVA->y() - pVD->y());
  adEdge[2] = (pVA->z() - pVD->z());
  double dDotD = dDOT3D(adEdge, adBisectorNorm);

  // Use these distances to find where the opposite edge
  // intersects the angle bisector plane.
  double dSum = fabs(dDotC) + fabs(dDotD);
  double dFrac = fabs(dDotC) / dSum;

  double adNewPoint[3];
  adNewPoint[0] = pVC->x() * (1. - dFrac) + pVD->x() * dFrac;
  adNewPoint[1] = pVC->y() * (1. - dFrac) + pVD->y() * dFrac;
  adNewPoint[2] = pVC->z() * (1. - dFrac) + pVD->z() * dFrac;

#ifndef NDEBUG
  {
    double adTempEdge[3];
    adTempEdge[0] = adNewPoint[0] - pVA->x();
    adTempEdge[1] = adNewPoint[1] - pVA->y();
    adTempEdge[2] = adNewPoint[2] - pVA->z();
    assert(qFuzzyPerp3D(adBisectorNorm, adTempEdge));
  }
#endif

  // Insert this point.
  Subseg * subseg = areSubsegVertsConnected(pVC, pVD);

  Vert *pVNew = SI3D.insertPointOnEdge(adNewPoint, pVC, pVD);

  if (subseg && pVNew->isValid()) {
    double splitParam;
    subseg->compute_refinement_split_data(adNewPoint, splitParam, true);
    createNewSubsegsBySplit(subseg, adNewPoint, splitParam, pVNew);
    pVNew->setType(Vert::eBdryCurve);
  }
  if (!pVNew->isValid()) {
    assert(!areBdryChangesAllowed());
    // Site not inserted because it was on a bdry edge and bdry's are
    // immutable.
    return 0;
  }

  // Smooth nearby.
  int iSmooths = Smooth.smoothVert(pVNew);

  logMessage(4,
	     "Point inserted %s to bisect dihedral angle of %6.2f degrees.\n",
	     (iLargestInt != -1 ? "in interior" : "on boundary"),
	     (iLargestInt != -1 ? dLargestInt : dLargest));
  logMessage(4, "  %3d swaps were made; vertex was %smoved by smoothing.\n",
	     SI3D.getNumSwapsDone(), (iSmooths == 1) ? "" : "not ");

  // Swap nearby.
  // Not yet implemented; swapping will be taken care of elsewhere.
  sendEvents();
  return 1;
}

int
VolMesh::repairBadCells(bool qAllowInsertion)
{
  int iPtsInserted = 0;

  GRUMMP::MaxMinSineSwapDecider3D SDMMS(false);
  GRUMMP::SwapManager3D SineSwap(&SDMMS, this);
  GRUMMP::OptMSSmoothingManager3D Smooth(this);

#define NON_TET 0
#define ROUND 1
#define NEEDLE 2
#define WEDGE 3
#define SPINDLE 4
#define SLIVER 5
#define CAP 6

  int iDum;
  for (iDum = 0; iDum < 10; iDum++)
    aiEdgeReq[iDum] = aiEdgeDone[iDum] = 0;

  logMessage(1, "Attempting to repair bad cells by construction...\n");
  int iMaxPasses = 3, iPasses = 0;
  int iRecentSuccesses, iRecentAttempts;
  // int iNDeleted = 0;
  do {
    // Set the limiting values for cells we'll try to improve.  These
    // get updated every pass because the global actual values do, too.
    double dLocSmallSolid = dSmallSolid + dSmallSolidMargin;
    double dLocLargeSolid = dLargeSolid - dLargeSolidMargin;
    double dLocSmallDihed = dSmallDihed + dSmallDihedMargin;
    double dLocLargeDihed = dLargeDihed - dLargeDihedMargin;
    // Limit the cells that we will choose to pay attention to
    dLocSmallSolid = min(5., dLocSmallSolid);
    dLocLargeSolid = max(210., dLocLargeSolid);
    dLocSmallDihed = min(20., dLocSmallDihed);
    dLocLargeDihed = max(140., dLocLargeDihed);
    double dSplitAngle = max(dLocLargeDihed, 150.);

    double dMaxDihed = 0, dMinDihed = 180;
    double dMaxSolid = 0, dMinSolid = 360;
    iPasses++;
    logMessage(2, "Pass %d\n", iPasses);
    iRecentSuccesses = 0;
    iRecentAttempts = 0;

    GR_index_t aiAttempts[] =
      { 0, 0, 0, 0, 0, 0, 0 };
    GR_index_t aiSuccesses[] =
      { 0, 0, 0, 0, 0, 0, 0 };
    for (int iC = getNumCells() - 1; iC >= 0; iC--) {
      int iType;
      Cell *pC = getCell(iC);
      if (pC->doFullCheck() != 1)
	continue;
      // Only try to fix tets that are not known to be well-shaped
      if (pC->getType() == Cell::eTet) {
	double adDihed[6];
	int iNDihed;
	dynamic_cast<TetCell*>(pC)->calcAllDihed(adDihed, &iNDihed);

	int iSmall = 0, iLarge = 0, i;
	for (i = 0; i < iNDihed; i++) {
	  if (adDihed[i] < dLocSmallDihed)
	    iSmall++;
	  else if (adDihed[i] > dLocLargeDihed)
	    iLarge++;
	}

	double dS1 = adDihed[0] + adDihed[1] + adDihed[2] - 180;
	double dS2 = adDihed[0] + adDihed[3] + adDihed[4] - 180;
	double dS3 = adDihed[1] + adDihed[3] + adDihed[5] - 180;
	double dS4 = adDihed[2] + adDihed[4] + adDihed[5] - 180;
	int iSmallSolid = ((dS1 < dLocSmallSolid ? 1 : 0)
	    + (dS2 < dLocSmallSolid ? 1 : 0) + (dS3 < dLocSmallSolid ? 1 : 0)
	    + (dS4 < dLocSmallSolid ? 1 : 0));
	int iLargeSolid = ((dS1 > dLocLargeSolid ? 1 : 0)
	    + (dS2 > dLocLargeSolid ? 1 : 0) + (dS3 > dLocLargeSolid ? 1 : 0)
	    + (dS4 > dLocLargeSolid ? 1 : 0));
	assert(iLarge <= 3);
	assert(iSmall <= 4);
	// By the definitions used here, there can actually be more than
	// one large solid angle if the largest solid angle in the mesh
	// is sufficiently small.
	//	assert(iLargeSolid <= 1);

	if (iLargeSolid >= 1)
	  iType = CAP;
	else if (iLarge > 0 && iSmall > 0)
	  iType = SLIVER;
	else if (iLarge > 0)
	  iType = SPINDLE;
	else if (iSmall > 0)
	  iType = WEDGE;
	else if (iSmallSolid > 0) {
	  iType = NEEDLE;
	  // pC->vSetFlag();
	  // continue; // Skip to next cell; there's no work to be done here
	}
	else {
	  iType = ROUND;
	  continue; // Skip to next cell; there's no work to be done here
	}
	assert(
	    iType == CAP || iType == SLIVER || iType == WEDGE || iType == SPINDLE || iType == NEEDLE);

	int aqTried[] =
	  { 0, 0, 0, 0, 0, 0 }, iSuccess = 0, iDihed;
	int aiPrec[] =
	  { 0, 0, 0, 0 };
	Vert *apV[4];
	for (i = 0; i < 4; i++)
	  apV[i] = pC->getVert(i);

	// Make note that we found and tried to fix this cell
	aiAttempts[iType]++;
	iRecentAttempts++;

	// Start by trying to remove edges w/ small dihedrals or
	// remove one of the verts opposite those edges
	//
	// Spindles don't have small dihedral angles, and caps may not,
	// but the check is absurdly cheap, so who cares?
	for (iDihed = 0; iDihed < 6 && !iSuccess; iDihed++)
	  if (adDihed[iDihed] < dLocSmallDihed) {
	    aqTried[iDihed] = 1;
	    // Mark these verts as less likely candidates for removal
	    aiPrec[aiVert0[iDihed]]++;
	    aiPrec[aiVert1[iDihed]]++;
	    Face *pF = pC->getOppositeFace(apV[aiVertOpp[iDihed]]);
	    iSuccess = edgeSwap3D(pF, apV[aiVert0[iDihed]],
				  apV[aiVert1[iDihed]],
				  apV[aiVertOther[iDihed]]);
	  }

	// Failing that, try to remove edges w/ large dihedrals
	//
	// Wedges don't have these, and caps actually might not either.
	// Again, who cares?
	for (iDihed = 0; iDihed < 6 && !iSuccess; iDihed++)
	  if (adDihed[iDihed] > dLocLargeDihed) {
	    aqTried[iDihed] = 1;
	    // Mark these verts as less likely candidates for removal
	    aiPrec[aiVert0[iDihed]]++;
	    aiPrec[aiVert1[iDihed]]++;
	    Face *pF = pC->getOppositeFace(apV[aiVertOpp[iDihed]]);
	    iSuccess = edgeSwap3D(pF, apV[aiVert0[iDihed]],
				  apV[aiVert1[iDihed]],
				  apV[aiVertOther[iDihed]]);
	  }

	// Failing that, try to remove any edge
	for (iDihed = 0; iDihed < 6 && !iSuccess; iDihed++)
	  if (aqTried[iDihed] == 0) {
	    aqTried[iDihed] = 1;
	    Face *pF = pC->getOppositeFace(apV[aiVertOpp[iDihed]]);
	    iSuccess = edgeSwap3D(pF, apV[aiVert0[iDihed]],
				  apV[aiVert1[iDihed]],
				  apV[aiVertOther[iDihed]]);
	  }

	// If _that_ still doesn't do anything, try face swapping.
	for (int iFace = 0; iFace < 4 && !iSuccess; iFace++) {
	  Face *pF = pC->getFace(iFace);
	  iSuccess = SineSwap.swapFace(pF, false);//iFaceSwap_deprecated(pF);
	}

	// If swapping has failed, break the largest angle of slivers
	// and spidles, which have large dihedral angles.
	if (qAllowInsertion && !iSuccess
	    && (iType == SLIVER || iType == SPINDLE || iType == CAP)) {
	  iSuccess = splitLargeAngle(pC, dSplitAngle, adDihed);
	  if (iSuccess) {
	    iPtsInserted++;
//	    // Now try to remove the new point.  This may or may not
//	    // result in different connectivity than before.
//	    int iNewSwaps;
//	    (void) qRemoveVert(pVVert(iNumVerts()-1), &iNewSwaps);
	  }
	}

	// If no edge or face was swapped away, try to remove a vertex,
	// beginning with vertices that are most often opposite bad
	// edges.
//	if (!iSuccess) {
//	  for (int iPrec = 0; iPrec <= 0 && !iSuccess; iPrec++)
//	    for (int iVert = 0; iVert <= 3 && !iSuccess; iVert++)
//	      if (aiPrec[iVert] == iPrec) {
//		int iDum;
//		iSuccess = qRemoveVert(apV[iVert], &iDum);
//		if (iSuccess) {
//		  iNDeleted++;
//		  vMessage(1, "Vertex deleted\n");
//		}
//	      }
//	}

	if (iSuccess) {
	  aiSuccesses[iType]++;
	  iRecentSuccesses++;
	  sendEvents();
	}
	else {
	  // The limits on bad angles for the next pass are based only
	  // on tets that were not improved.
	  dMaxSolid = max(dMaxSolid, max(max(dS1, dS2), max(dS3, dS4)));
	  dMinSolid = min(dMinSolid, min(min(dS1, dS2), min(dS3, dS4)));

	  for (i = 0; i < iNDihed; i++) {
	    dMaxDihed = max(dMaxDihed, adDihed[i]);
	    dMinDihed = min(dMinDihed, adDihed[i]);
	  }
	  logMessage(4, "%d\n", iType);
	  logMessage(4, "%12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n",
		     adDihed[0], adDihed[1], adDihed[2], adDihed[3], adDihed[4],
		     adDihed[5]);
	  Vert *pV = pC->getVert(0);
	  logMessage(4, "%6u %12.6g %12.6g %12.6g\n", getVertIndex(pV), pV->x(),
		     pV->y(), pV->z());
	  pV = pC->getVert(1);
	  logMessage(4, "%6u %12.6g %12.6g %12.6g\n", getVertIndex(pV), pV->x(),
		     pV->y(), pV->z());
	  pV = pC->getVert(2);
	  logMessage(4, "%6u %12.6g %12.6g %12.6g\n", getVertIndex(pV), pV->x(),
		     pV->y(), pV->z());
	  pV = pC->getVert(3);
	  logMessage(4, "%6u %12.6g %12.6g %12.6g\n\n", getVertIndex(pV),
		     pV->x(), pV->y(), pV->z());
	} // End of messages for failure
      } // End of code to fix tets
    } // End of loop over cells
    logMessage(2, "Fixed%5u out of%6u needles\n", aiSuccesses[NEEDLE],
	       aiAttempts[NEEDLE]);
    logMessage(2, "Fixed%5u out of%6u wedges\n", aiSuccesses[WEDGE],
	       aiAttempts[WEDGE]);
    logMessage(2, "Fixed%5u out of%6u spindles\n", aiSuccesses[SPINDLE],
	       aiAttempts[SPINDLE]);
    logMessage(2, "Fixed%5u out of%6u slivers\n", aiSuccesses[SLIVER],
	       aiAttempts[SLIVER]);
    logMessage(2, "Fixed%5u out of%6u caps\n\n", aiSuccesses[CAP],
	       aiAttempts[CAP]);

    dSmallSolid = dMinSolid;
    dLargeSolid = dMaxSolid;
    dSmallDihed = dMinDihed;
    dLargeDihed = dMaxDihed;
    logMessage(2, "Smallest dihedral angle now %6.2f degrees\n", dSmallDihed);
    logMessage(2, "Largest  dihedral angle now %6.2f degrees\n", dLargeDihed);
    logMessage(2, "Smallest solid    angle now %6.2f degrees\n", dSmallSolid);
    logMessage(2, "Largest  solid    angle now %6.2f degrees\n\n", dLargeSolid);

    sendEvents();
    // Keep going as long as the success rate is reaonable
  }
  while ((iRecentSuccesses > int(floor(0.05 * iRecentAttempts))
      && iPasses < iMaxPasses) || (iRecentSuccesses == 0 && iPasses == 1));

  purgeCells();
  purgeFaces();
  purgeBdryFaces();
  setAllHintFaces();
  setVertFaceNeighbors();

  logMessage(1, "Done trying to fix up bad cells\n\n");
//   vMessage(2, "Removed %d vertices.\n", iNDeleted);
  logMessage(2, "Edge swapping statistics.\n");
  logMessage(2, "Swap A tets for B   Considered   Done\n");
  for (iDum = 3; iDum < 10; iDum++)
    logMessage(2, "%7d for %2d %12d %8d\n", iDum, 2 * iDum - 4, aiEdgeReq[iDum],
	       aiEdgeDone[iDum]);

//   if (qAllowInsertion)
//     iPtsInserted += iBreakBadBdryCells();
  if (!(isValid()))
    vFoundBug("repair of bad tetrahedra");
  return (iPtsInserted);
}

int
VolMesh::breakBadBdryCells(bool qAllowInsertion)
// A major obstacle to 3D mesh improvement is tetrahedra which have
// all their vertices on the boundary of the domain.  Smoothing can't
// do anything to improve these tets, and even surface smoothing would
// probably be limited in effectiveness in many cases.  This routine
// determines which tets with all their points on the boundary are of
// poor quality and attempts to break them up.
{
  int iPtsInserted = 0;
  if (qAllowInsertion) {
    // Some accounting.
    int iNBdryTets = 0, iNBTetsFixed = 0;
    int aiBTetsFound[] =
      { 0, 0, 0, 0, 0 };
    int aiBTetsFixed[] =
      { 0, 0, 0, 0, 0 };
    double dSplitAngle = max(dLargeDihed - dLargeDihedMargin, 160.);

    // Flag all vertices as boundary or interior.
//    for (int iVert = getNumVerts() - 1; iVert >= 0; iVert--)
//      getVert(iVert)->setType(Vert::eInterior);
//
//    for (int iBFace = getNumBdryFaces() - 1; iBFace >= 0; iBFace--) {
//      BFace *pBF = getBFace(iBFace);
//      if (pBF->isDeleted()) continue;
//      for (int iV = 0; iV < pBF->getNumVerts(); iV++)
//	pBF->getVert(iV)->setType(Vert::eBdry);
//    }

    for (int iCell = getNumCells() - 1; iCell >= 0; iCell--) {
      Cell *pC = getCell(iCell);
      // Only interested in valid tets with all four verts on the boundary
      if (pC->getType() != Cell::eTet || pC->doFullCheck() != 1)
	continue;

      const Vert *pV0 = pC->getVert(0);
      const Vert *pV1 = pC->getVert(1);
      const Vert *pV2 = pC->getVert(2);
      const Vert *pV3 = pC->getVert(3);

      Face *pF0 = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV0);
      Face *pF1 = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV1);
      Face *pF2 = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV2);
      Face *pF3 = dynamic_cast<TetCell*>(pC)->getOppositeFace(pV3);

      int iNBVerts = (pV0->isBdryVert() ? 1 : 0) + (pV1->isBdryVert() ? 1 : 0)
	  + (pV2->isBdryVert() ? 1 : 0) + (pV3->isBdryVert() ? 1 : 0);

      if (iNBVerts != 4)
	continue;
      iNBdryTets++;
      // For each tet that passed the above disqualifiers, count the
      // number of boundary faces.
      bool aqIsBFace[4];
      aqIsBFace[0] = pF0->getOppositeCell(pC)->getType() == Cell::eTriBFace
	  || pF0->getOppositeCell(pC)->getType() == Cell::eIntTriBFace;
      aqIsBFace[1] = pF1->getOppositeCell(pC)->getType() == Cell::eTriBFace
	  || pF1->getOppositeCell(pC)->getType() == Cell::eIntTriBFace;
      aqIsBFace[2] = pF2->getOppositeCell(pC)->getType() == Cell::eTriBFace
	  || pF2->getOppositeCell(pC)->getType() == Cell::eIntTriBFace;
      aqIsBFace[3] = pF3->getOppositeCell(pC)->getType() == Cell::eTriBFace
	  || pF3->getOppositeCell(pC)->getType() == Cell::eIntTriBFace;
      int iNBFaces = (aqIsBFace[0] ? 1 : 0) + (aqIsBFace[1] ? 1 : 0)
	  + (aqIsBFace[2] ? 1 : 0) + (aqIsBFace[3] ? 1 : 0);
      aiBTetsFound[iNBFaces]++;

      double adDihed[6];
      int iNDihed;
      pC->calcAllDihed(adDihed, &iNDihed);
      assert(iNDihed == 6);
      double dMaxDihed = 0, dMinDihed = 180;
      for (int ii = 0; ii < 6; ii++) {
	dMaxDihed = max(dMaxDihed, adDihed[ii]);
	dMinDihed = min(dMinDihed, adDihed[ii]);
      }

      // Tets that are well-shaped and connect only contiguous parts of
      // the boundary can be ignored here.
      if ((iNBFaces == 2 || iNBFaces == 3)
	  && dMaxDihed < dLargeDihed - dLargeDihedMargin
	  && dMinDihed > dSmallDihed + dSmallDihedMargin)
	continue;

      // Cell splitting is done according to the following rationale:
      // 1.  Any cell with a large angle needs to have that angle split.
      // 2.  Cells with two or three boundary faces are allowed to have
      //     small dihedral angles only at edges common to multiple
      //     boundary faces in the cell; other cells in this category with
      //     a small dihedral angle should have their LARGEST dihedral
      //     angle opposite an interior edge split to allow some wiggle
      //     room for reconnection and smoothing, provided that edge has a
      //     sufficiently large dihedral angle itself.
      // 3.  Any cell with no boundary faces is bridging among two
      //     or more boundaries and should be broken up if possible.
      //     Again this is done by splitting the largest dihedral angle
      //     available.

      if (dMaxDihed > dLargeDihed - dLargeDihedMargin) {
	logMessage(2, "%s%d bdry faces. Max = %6.2f. Min = %6.2f\n",
		   "Breaking tet w/ all pts on bdry. ", iNBFaces, dMaxDihed,
		   dMinDihed);
	if (splitLargeAngle(pC, dSplitAngle, adDihed)) {
	  iPtsInserted++;
	  iNBTetsFixed++;
	  aiBTetsFixed[iNBFaces]++;
	}
	continue;
      }

      switch (iNBFaces)
	{
	case 3:
	  logMessage(2, "Three bdry faces");
	  break; // Might eventually split the one interior face, but
		 // there are no interior edges.
	case 2:
	  if (dMinDihed < dSmallDihed + dSmallDihedMargin) {
	    // Find the common edge
	    int iComEdge = -1;
	    if (aqIsBFace[0] && aqIsBFace[1])
	      iComEdge = 5;
	    else if (aqIsBFace[0] && aqIsBFace[2])
	      iComEdge = 4;
	    else if (aqIsBFace[0] && aqIsBFace[3])
	      iComEdge = 3;
	    else if (aqIsBFace[1] && aqIsBFace[2])
	      iComEdge = 2;
	    else if (aqIsBFace[1] && aqIsBFace[3])
	      iComEdge = 1;
	    else if (aqIsBFace[2] && aqIsBFace[3])
	      iComEdge = 0;
	    assert(iComEdge != -1);

	    // Do any _other_ edges have small dihedral angles?
	    bool qOtherSmall = false;
	    for (int i = 0; i < 6 && !qOtherSmall; i++) {
	      if (i == iComEdge)
		continue;
	      if (adDihed[i] < dSmallDihed + dSmallDihedMargin)
		qOtherSmall = true;
	    }

	    // If so, split the largest dihedral opposite an interior edge,
	    // with the caveat that the angle split be larger than 60 degrees.
	    logMessage(2, "%s%d bdry faces. Max = %6.2f. Min = %6.2f\n",
		       "Breaking tet w/ all pts on bdry. ", iNBFaces, dMaxDihed,
		       dMinDihed);
	    if (splitLargeAngle(pC, 60, adDihed, true)) {
	      iNBTetsFixed++;
	      aiBTetsFixed[2]++;
	      iPtsInserted++;
	    }
	  }
	  break;
	case 1:
	  {
	    // Check out whether the fourth vertex is legitimately connected.
	    bool qOK = false;
	    if (aqIsBFace[0])
	      qOK = qIsBdryEdge(pF2, pV0, pV1) || qIsBdryEdge(pF3, pV0, pV2)
		  || qIsBdryEdge(pF1, pV0, pV3);
	    else if (aqIsBFace[1])
	      qOK = qIsBdryEdge(pF2, pV1, pV0) || qIsBdryEdge(pF3, pV1, pV2)
		  || qIsBdryEdge(pF0, pV1, pV3);
	    else if (aqIsBFace[2])
	      qOK = qIsBdryEdge(pF1, pV2, pV0) || qIsBdryEdge(pF3, pV2, pV1)
		  || qIsBdryEdge(pF0, pV2, pV3);
	    else if (aqIsBFace[3])
	      qOK = qIsBdryEdge(pF1, pV3, pV0) || qIsBdryEdge(pF2, pV3, pV1)
		  || qIsBdryEdge(pF0, pV3, pV2);

	    // If legitimately connected, treat like the two-face case.
	    if (qOK) {
	      if (dMinDihed < dSmallDihed + dSmallDihedMargin) {
		// Split the largest dihedral opposite an interior edge, with
		// the caveat that the angle split be larger than 60 degrees.
		logMessage(2, "%s%d bdry faces. Max = %6.2f. Min = %6.2f\n",
			   "Breaking tet w/ all pts on bdry. ", iNBFaces,
			   dMaxDihed, dMinDihed);
		if (splitLargeAngle(pC, 60, adDihed, true)) {
		  iNBTetsFixed++;
		  aiBTetsFixed[1]++;
		  iPtsInserted++;
		}
	      }
	      break;
	    }
	  }
	  // If we were breaking the zero-face case, deliberate fall through
	  // might be a good choice here.
	  break;
	case 0:
	  // Unconditionally split the largest angle in the tet to separate
	  // the boundaries.
	  // Currently not done. Re-examine for 0.1.8.
	  //       iNBTetsFixed++;
	  //       aiBTetsFixed[iNBFaces]++;
	  //       vMessage(2, "%s%d bdry faces. Max = %6.2f. Min = %6.2f\n",
	  //	       "Breaking tet w/ all pts on bdry. ",
	  //	       iNBFaces, dMaxDihed, dMinDihed);
	  //       iSplitLargeAngle(pC, 0, adDihed);
	  //       iPtsInserted++;
	  break;
	default: // Can't ever get here legitimately, except for meshes with
		 // only a single tet.
	  if (getNumCells() != 1)
	    vWarning("Your mesh appears to have an isolated tetrahedron.\n");
	  break;
	} // end switch on number of boundary faces for the tet
    } // end looping over tets
    assert(
	iNBdryTets
	    == aiBTetsFound[0] + aiBTetsFound[1] + aiBTetsFound[2]
		+ aiBTetsFound[3] + aiBTetsFound[4]);
    assert(
	iNBTetsFixed
	    == aiBTetsFixed[0] + aiBTetsFixed[1] + aiBTetsFixed[2]
		+ aiBTetsFixed[3] + aiBTetsFixed[4]);
    logMessage(1, "Done trying to fix tets with all points on the boundary.\n");
    logMessage(
	2, " Found %5d tets with 0 bdry faces, of which %5d needed fixing\n",
	aiBTetsFound[0], aiBTetsFixed[0]);
    logMessage(2, "       %5d           1                      %5d\n",
	       aiBTetsFound[1], aiBTetsFixed[1]);
    logMessage(2, "       %5d           2                      %5d\n",
	       aiBTetsFound[2], aiBTetsFixed[2]);
    logMessage(2, "       %5d           3                      %5d\n",
	       aiBTetsFound[3], aiBTetsFixed[3]);
    logMessage(2, "       %5d           4                      %5d\n",
	       aiBTetsFound[4], aiBTetsFixed[4]);
    logMessage(
	1, " Found %5d boundary tets overall,  of which %5d needed fixing\n",
	iNBdryTets, iNBTetsFixed);
    purgeCells();
    for (GR_index_t i = 0; i < getNumCells(); i++)
      assert(getCell(i)->doFullCheck() == 1);
    purgeFaces();
    setAllHintFaces();
    setVertFaceNeighbors();

    if (!(isValid()))
      vFoundBug("repair of bad boundary tetrahedra");
  }
  return (iPtsInserted);
}

void
repairBadCells(VolMesh* pVM)
{
  // Create a smoothing manager and a swapping manager.
  GRUMMP::MaxMinSineSwapDecider3D SDMMS(false);
  GRUMMP::SwapManager3D SineSwap(&SDMMS, pVM);

  GRUMMP::OptMSSmoothingManager3D Smooth(pVM);

  int iSMSwaps, iSmooths;
  iSMSwaps = SineSwap.swapAllFaces();
  logMessage(2, "%d swaps before global edge swapping.\n", iSMSwaps);

  // Try swapping all the edges in the mesh.  This is expensive,
  // especially since edges are checked multiple times, but it's just a
  // test...
  EntContainer<TriFace>::iterator ECiter = pVM->triFace_begin(), ECiterEnd =
      pVM->triFace_end();
  int numSwaps = 0;
  for (; ECiter != ECiterEnd; ++ECiter) {
    int swapsThisFace = 0;
    TriFace *pF = &(*ECiter);
    assert(pF->isValid() && !pF->isDeleted());
    Vert *pV0 = pF->getVert(0);
    Vert *pV1 = pF->getVert(1);
    Vert *pV2 = pF->getVert(2);
    {
      GRUMMP::FaceSwapInfo3D FC(pF, pV0, pV1, pV2);
      swapsThisFace = SineSwap.swapFace(FC);
    }
    if (swapsThisFace == 0) {
      GRUMMP::FaceSwapInfo3D FC(pF, pV1, pV2, pV0);
      swapsThisFace += SineSwap.swapFace(FC);
    }
    if (swapsThisFace == 0) {
      GRUMMP::FaceSwapInfo3D FC(pF, pV2, pV0, pV1);
      swapsThisFace += SineSwap.swapFace(FC);
    }
    numSwaps += swapsThisFace;
  }

  logMessage(2, "Did %d swaps via swapping every edge.\n", numSwaps);

  iSMSwaps = SineSwap.swapAllQueuedFaces(false);
  iSmooths = Smooth.smoothAllQueuedVerts();
  logMessage(2, "Swapped and smoothed the whole mesh; %d swaps, %d smoothed\n",
	     iSMSwaps, iSmooths);
  pVM->sendEvents();
  int numPasses = 0;
  while (iSMSwaps + iSmooths > 0 && numPasses < 10) {
    numPasses++;
    iSMSwaps = SineSwap.swapAllQueuedFaces();
    iSmooths = Smooth.smoothAllQueuedVerts();
//     vMessage(2, "Pass %d: Swapped %d, smoothed %d\n", numPasses,
// 	     iSMSwaps, iSmooths);
    pVM->sendEvents();
  }
  logMessage(2, "Finished swap/smooth propagation after edge swapping.\n");

  int fixed = 0, tried = 0;
  for (int iC = pVM->getNumTetCells() - 1; iC >= 0; iC--) {

    Cell *pC = pVM->getCell(iC);
    if (pC->isDeleted())
      continue;
    assert(pC->doFullCheck());
    double qual = SDMMS.evalTetQual(pC);
    if (qual < 0.6) {
      tried++;
      // Take extreme measures to fix it.
      // First, smooth all vertices.
      Smooth.smoothVert(pC->getVert(0));
      Smooth.smoothVert(pC->getVert(1));
      Smooth.smoothVert(pC->getVert(2));
      Smooth.smoothVert(pC->getVert(3));

      // Next, try to swap all faces.  At most one of these will succeed.
      int l_numSwaps = (SineSwap.swapFace(pC->getFace(0))
	  || SineSwap.swapFace(pC->getFace(1))
	  || SineSwap.swapFace(pC->getFace(2))
	  || SineSwap.swapFace(pC->getFace(3)));

      // If the tet still exists, try all edges.  This code actually
      // tries them all twice in the end.
      for (int ii = 0; ii < 4 && l_numSwaps == 0; ii++) {
	assert(!pC->isDeleted());
	Face *pF = pC->getFace(ii);
	assert(!pF->isDeleted());
	Vert *pV0 = pF->getVert(0);
	Vert *pV1 = pF->getVert(1);
	Vert *pV2 = pF->getVert(2);
	{
	  GRUMMP::FaceSwapInfo3D FC(pF, pV0, pV1, pV2);
	  l_numSwaps = SineSwap.swapFace(FC);
	}
	if (l_numSwaps == 0) {
	  GRUMMP::FaceSwapInfo3D FC(pF, pV1, pV2, pV0);
	  l_numSwaps += SineSwap.swapFace(FC);
	}
	if (l_numSwaps == 0) {
	  GRUMMP::FaceSwapInfo3D FC(pF, pV2, pV0, pV1);
	  l_numSwaps += SineSwap.swapFace(FC);
	}
      }

      // If successful, propagate swapping and smoothing.
      if (l_numSwaps > 0) {
	fixed++;
	numPasses = 0;
	do {
	  pVM->sendEvents();
	  numPasses++;
	  iSMSwaps = SineSwap.swapAllQueuedFaces();
	  iSmooths = Smooth.smoothAllQueuedVerts();
// 	  vMessage(2, "Pass %d: Swapped %d, smoothed %d\n", numPasses,
// 		   iSMSwaps, iSmooths);
	}
	while (iSMSwaps + iSmooths > 0 && numPasses < 10);
      } // Done with post-swap / smooth
    } // Done checking bad tets
  } // Done looping over tets
  logMessage(2, "Fixed %d tets out of %d tried.\n", fixed, tried);
}
