#include "GR_assert.h"
#include "GR_Mesh2D.h"
#include "GR_Geometry.h"
#include "GR_BFace.h"
#include "GR_Vec.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
static int iDelReq, iBdryDelReq, iBdryApexDelReq;
static int iDelDone, iBdryDelDone;

struct RemCell {
  Cell *pC;
  Face *pFOnlyV, *pFOnlyV1;
  RemCell() :
      pC(pCInvalidCell), pFOnlyV(pFInvalidFace), pFOnlyV1(pFInvalidFace)
  {
  }
  RemCell(Cell * const pCIn, Face * const pF0, Face * const pF1) :
      pC(pCIn), pFOnlyV(pF0), pFOnlyV1(pF1)
  {
  }
};

// To remove a given vertex V from a mesh, there must be a vertex V1
// adjacent to V which is visible from all other vertices adjacent to
// V.  (Note that this specifically does NOT require that the
// neighborhood be convex.)  If this requirement is met, then all faces
// and all cells incident on both V and V1 can be removed from the mesh.
// The faces in the deleted cells which are incident on only one of V
// and V1 occur in pairs and can be matched up and reconnected to ensure
// a valid mesh on output.

bool
Mesh2D::removeVert(Vert * const pV, int& iNewSwaps,
				Vert * apVPreferred[], const int iNPref)
// Remove vertex *pV from the mesh.  Decisions about whether *pV
// SHOULD be removed are made elsewhere.  The vertices listed in
// apVPreferred are given a weak preference in determining to which
// vertex pV should be moved.  This preference takes the form of
// checking these vertices first in making removal choices.
{
//	 Only works for triangular meshes at this point.
// 	attempt by Dan Zaide, to make this work, July 2014, with quads present.
//	assert (isSimplicial ());

  iNewSwaps = 0;
  iDelReq++;

  //  if (pV->getVertType () == Vert::eBdryTwoSide) {
  //    logMessage (1, "Deletion of verts on internal boundaries not yet supported.\n");
  //    return false;
  //  }


//TODO Shahzaib Added this to preserve verts on wake centerline...
    if (pV->getVertType () == Vert::ePseudoSurface) {
    //  logMessage (1, "Deletion of verts tagged as Pseudosurface is forbidden.\n");
      return false;
    }

  if (pV->getVertType() == Vert::eBdryApex) {
    logMessage(2, "Deletion of boundary apexes is forbidden.\n");
    iBdryApexDelReq++;

    assert(0);
    return false;
  }
  logMessage(4, "Deleting vertex at (%10f, %10f)\n", pV->x(), pV->y());

  // List all incident cells on *pV and all vertices in those cells.  As
  // a side effect, find out whether this is a boundary vertex.
  std::set<Cell*> spCIncident;
  std::set<Vert*> spVCand;
  std::set<BFace*> spBFIncident;
  std::set<Face*> spFUpdate;

  std::set<Cell*>::iterator iterC;
  std::set<Vert*>::iterator iterV;
  std::set<BFace*>::iterator iterBF;
  std::set<Face*>::iterator iterF;

  // Put the preferred verts at the head of the list.
  for (int ii = 0; ii < iNPref; spVCand.insert(apVPreferred[ii++])) {
  }
  bool qBdryVert = false;

  findNeighborhoodInfo(pV, spCIncident, spVCand, &spBFIncident, &qBdryVert,
		       &spFUpdate);

  if (pV->getVertType() == Vert::eBdry) {
    qBdryVert = true;
  }

  iterC = spCIncident.begin();
  int iReg = (*iterC)->getRegion();
  int iOtherReg = (*iterC)->getRegion();
  for (; iterC != spCIncident.end(); iterC++) {
    if ((*iterC)->getRegion() != iReg)
      iOtherReg = (*iterC)->getRegion();
  }

  assert(pV->getVertType() == Vert::eBdryTwoSide || iReg == iOtherReg);
  if (qBdryVert) {
    // Abort if boundaries are precious
    if (!areBdryChangesAllowed())
      return false;

    spVCand.clear();
    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
	iterBF++) {
      BFace *pBF = *iterBF;
      spVCand.insert(pBF->getVert(0));
      spVCand.insert(pBF->getVert(1));
    }
    spVCand.erase(pV);
    assert(spVCand.size() == spBFIncident.size());
    iBdryDelReq++;
    if (spVCand.size() != 2) {
      // This must be a bdry apex; don't try to remove it.
      return false;
    }

  }
  else {
    assert(spVCand.size() == spCIncident.size());
    assert(2 * spVCand.size() == spFUpdate.size());
  }

  // List all faces in incident cells that are opposite *pV
  std::vector<Face*> vecpFOpposite;
  std::vector<int> vecpFOppositeRegion;
  for (iterC = spCIncident.begin(); iterC != spCIncident.end(); iterC++) {
    Cell *pC = *iterC;
    if (pC->getType() == Cell::eQuadCell)
      return false;
    Face *pFOpp = dynamic_cast<TriCell *>(pC)->getOppositeFace(pV);
    vecpFOpposite.push_back(pFOpp);
    vecpFOppositeRegion.push_back(pC->getRegion());
    spFUpdate.erase(pFOpp);
  }

  // For each neighboring vertex, determine whether it is a valid target
  // for edge contraction.  Keep track of which faces affect which
  // vertices.  If no vertices are okay, start with the one with the
  // fewest interfering faces and try to remove that face from the hull
  // by using edge removal to break the connection from *pV to any of
  // the vertices of the face.  If that succeeds, call this routine
  // recursively to do the contraction; otherwise continue until an
  // attempt has been made to remove each problematic face.

  bool qVertOK;
  bool qAnyVertOK = false;
  double dMinDihed;
  double dMaxMinDihed = 0;
  Vert *pV1 = pVInvalidVert, *pVBest = pVInvalidVert;
  for (iterV = spVCand.begin(); iterV != spVCand.end(); ++iterV) {
    qVertOK = true;
    dMinDihed = 180;
    pV1 = *iterV;
    for (GR_index_t iF = 0; iF < vecpFOpposite.size() && qVertOK; iF++) {
      Face *pF = vecpFOpposite[iF];
      if (!pF->hasVert(pV1)) {
	// Determine whether a given face has the correct orientation
	// for the candidate contraction target.  If the orientation
	// test fails, pV1 can't see the correct side of the face pF.
	bool qResult;
	if ((pF->getLeftCell())->hasVert(pV)) {
	  // For this case, the right-hand rule applies for orientation.
	  assert(checkOrient2D(pF->getVert(0), pF->getVert(1), pV) != -1);
	  qResult = (checkOrient2D(pF->getVert(0), pF->getVert(1), pV1) == 1);
	  if (qResult) {
	    Vert* apV[] =
	      { pF->getVert(0), pF->getVert(1), pV1 };
	    TriCellCV FT(apV);
	    double dMin = FT.calcMinDihed();
	    double dMax = FT.calcMaxDihed();
	    dMinDihed = min(dMinDihed, min(dMin, 180 - dMax));
	    logMessage(4, "Min dihed (face): %10.6f  (vert): %10.6f\n",
		       FT.calcMinDihed(), dMinDihed);
	  }
	}
	else {
	  // Otherwise, the left-hand rule applies for orientation.
	  assert(checkOrient2D(pF->getVert(1), pF->getVert(0), pV) != -1);
	  qResult = (checkOrient2D(pF->getVert(1), pF->getVert(0), pV1) == 1);
	  if (qResult) {
	    Vert *apV[] =
	      { pF->getVert(1), pF->getVert(0), pV1 };
	    TriCellCV FT(apV);
	    double dMin = FT.calcMinDihed();
	    double dMax = FT.calcMaxDihed();
	    dMinDihed = min(dMinDihed, min(dMin, 180 - dMax));
	    logMessage(4, "Min dihed (face): %10.6f  (vert): %10.6f\n",
		       FT.calcMinDihed(), dMinDihed);
	  }
	}
	if (!qResult) {
	  qVertOK = false;
	}
      }
    }
    if (qVertOK) {
      qAnyVertOK = true;
      if (dMinDihed > dMaxMinDihed) {
	logMessage(4, "Best min dihedral yet.  Old: %10.6f  New: %10.6f\n",
		   dMaxMinDihed, dMinDihed);
	dMaxMinDihed = dMinDihed;
	pVBest = pV1;
      }
    }
  }

  if (!qAnyVertOK) {
    return false;
  }

  // If we get to here, we know we CAN remove *pV by contracting edge
  // (*pV, *pVBest), and that this gives the largest minimum dihedral
  // angle of any possible contraction.
  assert(qAnyVertOK);
  assert(pVBest->isValid());
  // If the best you can do is really lousy, don't do it.  This will
  // need careful modification for structured meshes.
  // if (dMaxMinDihed < 3.) return false;
  logMessage(4, "%s %u, with min dihedral of %5.2f degrees after.\n",
	     "Contracting to remove vertex", iVertIndex(pV), dMaxMinDihed);

  // At this point:
  //   spCIncident is a list of cells incident on *pV.
  //   spFUpdate contains all the faces incident on pV
  //   vecpFOpposite contains all the faces opposite pV

  // So:
  //   Delete all the cells in spCIncident
  //   For faces in vecpFOpposite that aren't incident on pVBest, create
  //     a new cell connecting them to pVBest
  //   For bdry verts, delete both BFace's incident on pV, identify the
  //     new face on the bdry incident on pVBest, and create a new BFace
  //     for it.
  //   Note that faces get deleted automagically.

  for (iterC = spCIncident.begin(); iterC != spCIncident.end(); iterC++) {
    Cell *pC = *iterC;
    deleteCell(pC);
  }

#ifndef NDEBUG
  std::set<Cell*> spCCheck;
#endif
  for (GR_index_t iF = 0; iF < vecpFOpposite.size(); iF++) {
    Face *pF = vecpFOpposite[iF];
    if (pF->hasVert(pVBest))
      continue;
#ifndef NDEBUG
    Cell *pCNew =
#endif
	createTriCell(pF->getVert(0), pF->getVert(1), pVBest,
		      vecpFOppositeRegion[iF]);
#ifndef NDEBUG
    spCCheck.insert(pCNew);
#endif
  }

  if (qBdryVert) {
    // Remove from the list of bfaces incident on *pV those not also
    // incident on *pVBest, to get a list of bfaces which will be
    // deleted.
    Vert *pVOtherBdry = pVInvalidVert;
    BFace *pBFOld = pBFInvalidBFace;
    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
	iterBF++) {
      BFace *pBF = *iterBF;
      if (pBF->hasVert(pV)) {
	if (!pBF->hasVert(pVBest)) {
	  pVOtherBdry =
	      pBF->getVert(0) == pV ? pBF->getVert(1) : pBF->getVert(0);
	  pBFOld = pBF;
	}
      }
    }
    Face * pF;
    if (pV->getVertType() == Vert::eBdryTwoSide) {
      pF = findCommonFace(pVBest, pVOtherBdry, true);
    }

    else
      pF = findCommonFace(pVBest, pVOtherBdry);

    BdryEdgeBase * pBEOld = dynamic_cast<BdryEdgeBase*>(pBFOld);
    assert(pBEOld);
    if (pBEOld->doesCurveExist()) {
      GRCurve* curve = pBEOld->getCurve();

//			Vert *vertA = forward ?  pF->getVert(0) :  pF->getVert(1);
//			Vert *vertB = forward ?  pF->getVert(1) :  pF->getVert(0);
      Vert *vertA = pF->getVert(0);
      Vert *vertB = pF->getVert(1);
      // now we need to know the parameters
      assert(vertA != vertB);

      double beginParam = curve->u_from_position(

      CubitVector(vertA->getCoords()));
      double endParam = curve->u_from_position(CubitVector(vertB->getCoords()));
      if (beginParam > endParam) {
	vertA = pF->getVert(1);
	vertB = pF->getVert(0);
	beginParam = curve->u_from_position(CubitVector(vertA->getCoords()));
	endParam = curve->u_from_position(CubitVector(vertB->getCoords()));
      }

      if (pBEOld->getType() == Cell::eIntBdryEdge) {

	createIntBdryEdge(vertA, vertB, curve, beginParam, endParam);

      }
      else {

	createBdryEdge(vertA, vertB, curve, beginParam, endParam);

      }
    }
    else

      createBFace(pF, pBFOld);
    // Fixes error with boundary orientation.
    // Basically checks if the orientation is correct w.r.t to the curve
    // after removal, interchanges if so.
//		BdryEdgeBase * pBEB = dynamic_cast<BdryEdgeBase*>(bface);
//		if(pBEB){
//			double param0 = pBEB->getVert0Param(false);
//			double param1 = pBEB->getVert1Param(false);
//			GRCurve * curve = pBEB->getCurve();
//			int iRegLeft  = curve->region_left();
//			int iRegRight = curve->region_right();
//			// okay so we want that if its forward,
//
//			printf("bface %f %f %d\n",pBEB->getVert0Param(false),pBEB->getVert1Param(false),pBEB->isForward());
//				pBEB->getVert(0)->printVertInfo();
//			pBEB->getVert(1)->printVertInfo();
//			if(pBEB->isForward() && param0 > param1){
//
//				for(int iF = 0; iF < pBEB->getNumFaces(); iF++){
//					pBEB->getFace(iF)->interchangeCellsAndVerts();
//					if (iRegRight == iInvalidRegion) {
//						if(pBEB->getFace(iF)->getRightCell()->getRegion() != iInvalidRegion){
//							pBEB->setVert0Param(param1);
//							pBEB->setVert1Param(param0);
//							pBEB->getFace(iF)->interchangeCellsAndVerts();
//						}
//					} else if(iRegLeft == iInvalidRegion) {
//						if(pBEB->getFace(iF)->getLeftCell()->getRegion() != iInvalidRegion){
//							pBEB->setVert0Param(param1);
//							pBEB->setVert1Param(param0);
//							pBEB->getFace(iF)->interchangeCellsAndVerts();
//						}
//					}
//				}
//
//			}
//		}

    for (iterBF = spBFIncident.begin(); iterBF != spBFIncident.end();
	iterBF++) {
      BFace *pBF = *iterBF;
      if (pBF->hasVert(pV)) {
	deleteBFace(pBF);
      }
    }
  }

#ifndef NDEBUG
  for (iterC = spCCheck.begin(); iterC != spCCheck.end(); iterC++) {
    Cell *pC = *iterC;
    assert(pC->hasVert(pVBest));
    assert(!pC->hasVert(pV));
    assert(pC->doFullCheck() == 1);
    assert(checkOrient2D(pC->getVert(0), pC->getVert(1), pC->getVert(2)) == 1);
  }
#endif

  // Delete *pV.
  deleteVert(pV);

  int iNumSwaps = 0;

  GRUMMP::SwapDecider2D *SwapDec2D;
  GRUMMP::SwapManager2D *SwapMan2D;
  GRUMMP::QualMeasure *pQM;

// Swap every face that was incident to the deleted vertex that survived
  pQM = new GRUMMP::Angle2D;
  SwapDec2D = new GRUMMP::DelaunaySwapDecider2D();
  SwapMan2D = new GRUMMP::SwapManager2D(SwapDec2D, this);

  GRUMMP::FaceSwapInfo2D FSI2D;
  for (iterF = spFUpdate.begin(); iterF != spFUpdate.end(); ++iterF) {
    Face *pF = *iterF;
    if (!pF->isDeleted()) {
      FSI2D.calculateFaceInfo(pF);
      iNumSwaps += SwapMan2D->swapFace(FSI2D);
    }
  }
  delete SwapMan2D;
  delete SwapDec2D;
  delete pQM;

  // Success!
  iDelDone++;
  if (qBdryVert)
    iBdryDelDone++;

  //updateLengthScale();

  return true;
}

// Yet another removal method: this one deleted everything in the
// neighborhood of the vertex being deleted, then remeshes the resulting
// star-shaped region.  (For 2D, this is guaranteed to work, so there's
// no problem with this.)  The remeshing uses a front (implemented as a
// set of faces), with each face in the front used to create a Delaunay
// triangle.  When creating triangles, new faces are created iff they
// don't already exist in the front.

// This code fails for boundary verts.

// class FrontFace2D {
//   Vert *pV0, *pV1;
//   Face *pF;
// public:
//   FrontFace2D(Vert * const pVA, Vert * const pVB,
// 	    Face * const pFIn = pFInvalidFace)
//     : pV0(pVA), pV1(pVB), pF(pFIn) {
//     assert(pV0 != pV1);
//   }
//   Vert *pVVert(const int i) const {
//     assert(i == 0 || i == 1);
//     return (i == 0 ? pV0 : pV1);
//   }
//   Face *pFFace() const {return pF;}
//   bool operator==(const FrontFace2D& FF) const {
//     return ((pV0 == FF.pV0 && pV1 == FF.pV1) ||
// 	    (pV0 == FF.pV1 && pV1 == FF.pV0));
//   }
//   bool operator<(const FrontFace2D& FF) const {
//     Vert *pVMinThis, *pVMaxThis, *pVMinOther, *pVMaxOther;
//     if (pV0 < pV1) {
//       pVMinThis = pV0; pVMaxThis = pV1;
//     }
//     else {
//       pVMinThis = pV1; pVMaxThis = pV0;
//     }
//     if (FF.pV0 < FF.pV1) {
//       pVMinOther = FF.pV0; pVMaxOther = FF.pV1;
//     }
//     else {
//       pVMinOther = FF.pV1; pVMaxOther = FF.pV0;
//     }
//     return ((pVMinThis < pVMinOther) ||
// 	    ((pVMinThis == pVMinOther) && (pVMaxThis < pVMaxOther)));
//   }
// };

bool
Mesh2D::removeVertByRemeshing(Vert * const pV, int& iNSwaps)
{
  logMessage(3, "In qRemoveVert for (%.3e,%.3e)\n", pV->x(), pV->y());
  iNSwaps = 0;
  // Find the neighborhood of pV.
  //  assert(0);
  std::set<Cell*> spCInc;
  std::set<Vert*> spVNeigh;
  std::set<Face*> spFNearby;

  {
    logMessage(3, "Finding the neighbourhood\n");
    std::set<BFace*> spBFDummy;
    bool qBdryVert;

    findNeighborhoodInfo(pV, spCInc, spVNeigh, &spBFDummy, &qBdryVert,
			 &spFNearby);

    if (qBdryVert) {
      assert(pV->isBdryVert());
      assert(spBFDummy.size() >= 2); // Could be greater for a corner
      // where regions meet.
      return false;
    }
    else {
      assert(!pV->isBdryVert());
      assert(spBFDummy.empty());
    }
  }

  std::set<Cell*>::iterator iterC = spCInc.begin();
  int iReg = (*iterC)->getRegion();
  for (; iterC != spCInc.end(); iterC++) {
    if ((*iterC)->getRegion() != iReg)
      return false;
  }

  // Rip out all the cells in the neighborhood...
  for (iterC = spCInc.begin(); iterC != spCInc.end(); iterC++) {
    Cell *pC = *iterC;

    deleteCell(pC);
  }

  // For faces, rip out those incident, and add the other to the set of
  // front faces.
  std::set<Face*> sFF;
  for (std::set<Face*>::iterator iterF = spFNearby.begin();
      iterF != spFNearby.end();) {
    Face *pF = *iterF;
    // The increment is done here, because if the face is deleted, iterF
    // would be invalidated.  So it has to be incremented -before- the
    // deletion might happen.
    std::set<Face*>::iterator iterTmp = iterF;
    ++iterF;
    if (pF->hasVert(pV)) {
      assert(pF->isDeleted());
      spFNearby.erase(iterTmp);
    }
    else {
      Vert *pV0 = pF->getVert(0);
      Vert *pV1 = pF->getVert(1);
      // Update the vertex hints for all vertices in the hull, since
      // their hints may currently be faces that are being deleted.
      pV0->setHintFace(pF);
      pV1->setHintFace(pF);

      sFF.insert(pF);
    }
  } // Done with loop over faces.

    // Mark the vertex as deleted.
  deleteVert(pV);
  int iNumNewCells = 0;
  logMessage(3, "Building the delaunay triangle\n");
  // For each face in the hull, build the Delaunay triangle that belongs
  // on it.
  do {
    assert(sFF.size() >= 3);
    // Maintaining this list correctly is easier if you just rebuild it
    // every time.
    std::set<Vert*> spVCand;
    {
      std::set<Face*>::iterator iter;
      for (iter = sFF.begin(); iter != sFF.end(); ++iter) {
	spVCand.insert((*iter)->getVert(0));
	spVCand.insert((*iter)->getVert(1));
      }
      assert(spVCand.size() <= sFF.size());
    }

    Face *pF = *sFF.begin();
    Vert *pV0, *pV1;
    if (pF->getLeftCell() == pCInvalidCell) {
      pV0 = pF->getVert(0);
      pV1 = pF->getVert(1);
    }
    else {
      pV0 = pF->getVert(1);
      pV1 = pF->getVert(0);
    }

    // Find the right vertex: the one on the correct side of the face
    // with no point in the resulting circumcircle.  In case of ties
    // between candidate verts, pick the one nearest the face midpoint.
    std::set<Vert*>::iterator iterV;
    for (iterV = spVCand.begin();
	iterV != spVCand.end() && checkOrient2D(pV0, pV1, *iterV) != 1;
	++iterV) {
    }
    assert(iterV != spVCand.end());
    Vert *pVBest = *iterV;
    ++iterV;
    for (; iterV != spVCand.end(); ++iterV) {
      Vert *pVCand = *iterV;
      if (pVCand == pV0 || pVCand == pV1)
	continue;
      if (checkOrient2D(pV0, pV1, pVCand) == 1) {
	int iIncircleResult = isInCircle(pV0, pV1, pVBest, pVCand);
	if (iIncircleResult == 1) {
	  pVBest = pVCand;
	}
	else if (iIncircleResult == 0) {
	  double adMid[] =
	    { 0.5 * (pV0->x() + pV1->x()), 0.5 * (pV0->y() + pV1->y()) };
	  if (dDIST2D(adMid, pVCand->getCoords()) <
	  dDIST2D(adMid, pVBest->getCoords()))
	    pVBest = pVCand;
	}
      }
    }
    // At this point, pVBest is definitely the one.

    // Create the new cell
    iNumNewCells++;
    Cell *pCNew = createTriCell(pV0, pV1, pVBest, iReg);
    for (int i = 0; i < 3; i++) {
      pF = pCNew->getFace(i);
      if (sFF.find(pF) != sFF.end()) {
	sFF.erase(pF);
      }
      else {
	sFF.insert(pF);
      }
    }
  }
  while (!sFF.empty());

  assert(iNumNewCells == int(spFNearby.size()) - 2);

  // Always succeeds.
  return true;
}

bool
Mesh2D::removeVertCheaply(Vert * const pV, int& iSwaps)
{
  // Simple:  remove edges until only three remain.  Then the vertex to
  // be removed is surrounded by three triangles which can be combined
  // to form one.  Afterwards, try to swap each edge that once shared a
  // triangle with the vertex.
  //
  // For boundary edges, when three edges remain, two are boundary
  // edges.  Those two edges are merged, the third disappears.
  //  printf("start ");
  if (pV->getVertType() == Vert::eBdryApex|| pV->getVertType() == Vert::ePseudoSurface || pV->getVertType() == Vert::eInteriorFixed) {
    return false;
  }
  // Can't remove bdry verts when the boundary is precious.
  if (!m_allowBdryChanges
      && (pV->getVertType() == Vert::eBdryApex
	  || pV->getVertType() == Vert::eBdryCurve
	  || pV->getVertType() == Vert::eBdryTwoSide
	  || pV->getVertType() == Vert::eBdry)) {
    assert(0);
    printf("Bdry Changes not allowed\n");
    return false;
  }
  bool qRetVal = true;
  iSwaps = 0;

  std::set<Face*> spFSwapEdges;
  {
    std::set<Vert*> spVTmp;
    std::set<Cell*> spCTmp;
    findNeighborhoodInfo(pV, spCTmp, spVTmp, NULL, NULL, &spFSwapEdges);

    std::set<Cell*>::iterator iter, iterEnd = spCTmp.end();
    for (iter = spCTmp.begin(); iter != iterEnd; ++iter) {
      if ((*iter)->getType() == Cell::eQuadCell)
	return false;
    }
  }
  std::set<Face*> spFIncEdges;
  std::set<Face*>::iterator iterF;
  for (iterF = spFSwapEdges.begin(); iterF != spFSwapEdges.end(); ++iterF) {
    Face *pF = *iterF;
    if (pF->hasVert(pV))
      spFIncEdges.insert(pF);
  }
  // Now have a list of every edge incident on the vertex.  Cycle
  // through the list, removing edges, until only three remain.
  int iLastLength = -1;
  while (spFIncEdges.size() > 3 && int(spFIncEdges.size()) != iLastLength) {
    iLastLength = spFIncEdges.size();
    for (iterF = spFIncEdges.begin(); iterF != spFIncEdges.end();) {
      Face *pF = *iterF;
      std::set<Face*>::iterator iterTmp = iterF;
      ++iterF;
      if (removeEdge(pF)) {
	spFIncEdges.erase(iterTmp);
      }
    }
  }
  int iNEdges = spFIncEdges.size();
  if (iNEdges <= 2 || iNEdges > 4) {
    // printf("iNEdges = %d, which is bad\n",iNEdges);
    // assert(0);
    //Something needs to be done here
    return false;
  }
  std::vector<Face*> vecpFIncEdges(iNEdges);
  std::copy(spFIncEdges.begin(), spFIncEdges.end(), vecpFIncEdges.begin());

  // List the cells incident on these edges
  std::set<Cell*> spCIncCells;
  for (int i = 0; i < iNEdges; i++) {
    spCIncCells.insert(vecpFIncEdges[i]->getLeftCell());
    spCIncCells.insert(vecpFIncEdges[i]->getRightCell());
  }
  for (int i = 0; i < iNEdges; i++)
    assert(vecpFIncEdges[i]->getVert(0) != vecpFIncEdges[i]->getVert(1));

  if (pV->getVertType() == Vert::eInterior) {
    // printf("Deleting interior\n");
    if (iNEdges == 4) {
      // printf("4 Edges crap\n");
      //Just take the 4 triangles and turn it into 2
      //First circle around to get the required info
      Face *pFA1 = vecpFIncEdges[0];
      Cell *pCD = pFA1->getRightCell();
      Vert *pVA = pFA1->getVert(0);
      if (pVA == pV)
	pVA = pFA1->getVert(1);

      Vert *pVB = pCD->getOppositeVert(pFA1);
      Face *pFB1 = pCD->getOppositeFace(pVA);
      Cell *pCA = pFB1->getOppositeCell(pCD);
      Vert* pVC = pCA->getOppositeVert(pFB1);
      Face *pFC1 = pCA->getOppositeFace(pVB);
      Cell* pCB = pFC1->getOppositeCell(pCA);
      Vert* pVD = pCB->getOppositeVert(pFC1);
      Face* pFD1 = pCB->getOppositeFace(pVC);
      Cell* pCC = pFD1->getOppositeCell(pCB);
      if (pCellQueue != NULL) {
	pCellQueue->vRemoveCell(pCA);
	pCellQueue->vRemoveCell(pCB);
	pCellQueue->vRemoveCell(pCC);
	pCellQueue->vRemoveCell(pCD);
	// 	pCellQueue->vRemoveFace(pFB1);  //just trying to find a bug
	// 	pCellQueue->vRemoveFace(pFC1);
	// 	pCellQueue->vRemoveFace(pFD1);
      }

      assert(spCIncCells.count(pCA) == 1);
      assert(spCIncCells.count(pCB) == 1);
      assert(spCIncCells.count(pCC) == 1);
      assert(spCIncCells.count(pCD) == 1);

      Face *pFA2 = dynamic_cast<TriCell*>(pCA)->getOppositeFace(pV);
      Face *pFB2 = dynamic_cast<TriCell*>(pCB)->getOppositeFace(pV);
      Face *pFC2 = dynamic_cast<TriCell*>(pCC)->getOppositeFace(pV);
      Face *pFD2 = dynamic_cast<TriCell*>(pCD)->getOppositeFace(pV);
      //Now I have all the required info
      int iRegA = pCA->getRegion();
      int iRegC = pCC->getRegion();
      deleteCell(pCA);
      deleteCell(pCB);
      deleteCell(pCC);
      deleteCell(pCD);
      deleteFace(pFA1);
      deleteFace(pFB1);
      deleteFace(pFC1);
      deleteFace(pFD1);

      deleteVert(pV);
      bool qExist;
      assert(pFA2->doFullCheck());
      assert(pFB2->doFullCheck());
      assert(pFC2->doFullCheck());
      assert(pFD2->doFullCheck());
      double ABDarea = dAreaTri(pVA, pVB, pVD);
      double BCDarea = dAreaTri(pVC, pVB, pVD);
      double dR1, dR2;
      if (ABDarea < BCDarea)
	dR1 = ABDarea / BCDarea;
      else
	dR1 = BCDarea / ABDarea;
      double ABCarea = dAreaTri(pVA, pVB, pVC);
      double ACDarea = dAreaTri(pVC, pVA, pVD);
      if (ABCarea < ACDarea)
	dR2 = ABCarea / ACDarea;
      else
	dR2 = ACDarea / ABCarea;

      if (dR1 > dR2) {
	Face *pFA = createFace(qExist, pVB, pVD);
	pCC = createTriCell(pFC2, pFA, pFD2, iRegC);
	assert(pCC->doFullCheck());
	pCA = createTriCell(pFB2, pFA2, pFA, iRegA);
	assert(pCA->doFullCheck());
      }
      else {
	Face *pFA = createFace(qExist, pVA, pVC);
	pCC = createTriCell(pFC2, pFB2, pFA, iRegC);
	assert(pCC->doFullCheck());
	pCA = createTriCell(pFA2, pFD2, pFA, iRegA);
	assert(pCA->doFullCheck());
      }
      if (pCellQueue != NULL) {
	pCellQueue->vAddCell(pCA);
	pCellQueue->vAddCell(pCC);
      }
      assert(pFA2->doFullCheck());
      assert(pFB2->doFullCheck());
      assert(pFC2->doFullCheck());
      assert(pFD2->doFullCheck());
      assert(pFA1->doFullCheck());
      for (iterF = spFSwapEdges.begin(); iterF != spFSwapEdges.end(); ++iterF) {
	Face *pF = *iterF;
	if (pF->isValid() && !pF->isDeleted()) {
	  // 	  pF->pVVert (0)->vMarkIllShaped ();
	  // 	  pF->pVVert (1)->vMarkIllShaped ();
	  iSwaps = iSwaps + iFaceSwap_deprecated(pF);
	}
      }

      return qRetVal;
    }

    if (iNEdges != 3 && spCIncCells.size() != 3) {
      printf("iNEdges = %d and spCInc has %zd\n", iNEdges, spCIncCells.size());
      assert(0);
    }
    // Identify all the necessary parts
    Face *pFA1 = vecpFIncEdges[0];
    Face *pFB1 = vecpFIncEdges[1];
    Face *pFC1 = vecpFIncEdges[2];

    Cell *pCA = findCommonCell(pFB1, pFC1);
    Cell *pCB = findCommonCell(pFC1, pFA1);
    Cell *pCC = findCommonCell(pFA1, pFB1);
    if (pCellQueue != NULL) {
      pCellQueue->vRemoveCell(pCA);
      pCellQueue->vRemoveCell(pCB);
      pCellQueue->vRemoveCell(pCC);
    }

    assert(spCIncCells.count(pCA) == 1);
    assert(spCIncCells.count(pCB) == 1);
    assert(spCIncCells.count(pCC) == 1);

    Face *pFA2 = dynamic_cast<TriCell*>(pCA)->getOppositeFace(pV);
    Face *pFB2 = dynamic_cast<TriCell*>(pCB)->getOppositeFace(pV);
    Face *pFC2 = dynamic_cast<TriCell*>(pCC)->getOppositeFace(pV);

    int iReg = pCA->getRegion();
    deleteCell(pCA);
    deleteCell(pCB);
    deleteCell(pCC);
    // deleteFace(pFA1);
    // deleteFace(pFB1);
    // deleteFace(pFC1);

    deleteVert(pV);
    pCA = createTriCell(pFA2, pFB2, pFC2, iReg);

    assert(pCA->doFullCheck());
    if (pCellQueue != NULL) {
      pCellQueue->vAddCell(pCA);
    }
    assert(pFA2->doFullCheck());
    assert(pFB2->doFullCheck());
    assert(pFC2->doFullCheck());
  }				// Done deleting interior vert
  else if (pV->getVertType() == Vert::eBdryTwoSide) {
    // Delete internal boundary vert
    // Version 0.1.8: Implement this
    printf("Internal boundaries not supported\n");
    assert(0);
    return false;
  }				// Done deleting internal boundary vert

  else {
    // printf("Deleting boundary vert\n");		// Delete boundary vert
    if (iNEdges == 3) {
      // Identify all the parts

      //Need to check if deleting this vertex will invert a triangle
      Vert* pVB1 = NULL;
      Vert* pVB2 = NULL;
      Vert* pVInt = NULL;
      for (int iF = 0; iF < 3; iF++) {
	if (vecpFIncEdges[iF]->isBdryFace() && pVB1 == NULL) {
	  pVB1 = vecpFIncEdges[iF]->getOppositeVert(pV);
	}
	else if (vecpFIncEdges[iF]->isBdryFace()) {
	  pVB2 = vecpFIncEdges[iF]->getOppositeVert(pV);
	}
	else {
	  pVInt = vecpFIncEdges[iF]->getOppositeVert(pV);
	}
      }
      assert(pVB1!=NULL && pVB2!=NULL && pVInt!=NULL);

      if (isInTri(pVInt, pV, pVB1, pVB2) != (-1))
	return false;

      Cell* apCIncCells[4];
      std::copy(spCIncCells.begin(), spCIncCells.end(), apCIncCells);
      Cell *pCA, *pCB;
      if (apCIncCells[0]->getType() == Cell::eTriCell) {
	pCA = apCIncCells[0];
	if (apCIncCells[1]->getType() == Cell::eTriCell)
	  pCB = apCIncCells[1];
	else if (apCIncCells[2]->getType() == Cell::eTriCell)
	  pCB = apCIncCells[2];
	else
	  pCB = apCIncCells[3];
      }
      else if (apCIncCells[1]->getType() == Cell::eTriCell) {
	pCA = apCIncCells[1];
	if (apCIncCells[2]->getType() == Cell::eTriCell)
	  pCB = apCIncCells[2];
	else
	  pCB = apCIncCells[3];
      }
      else {
	pCA = apCIncCells[2];
	pCB = apCIncCells[3];
      }
      assert(pCA->getType() == Cell::eTriCell);
      assert(pCB->getType() == Cell::eTriCell);
      assert(pCA != pCB);

      if (pCellQueue != NULL) {
	pCellQueue->vRemoveCell(pCA);
	pCellQueue->vRemoveCell(pCB);
      }
      Face *pFD = pFInvalidFace;
      for (int ii = 0; ii < 3; ii++) {
	pFD = vecpFIncEdges[ii];
	if (pFD->hasCell(pCA) && pFD->hasCell(pCB))
	  break;
      }
      assert(pFD->hasCell(pCA) && pFD->hasCell(pCB));

      Vert *pVA = dynamic_cast<TriCell*>(pCB)->getOppositeVert(pFD);
      Vert *pVB = dynamic_cast<TriCell*>(pCA)->getOppositeVert(pFD);
      Vert *pVC = (pFD->getVert(0) == pV) ? pFD->getVert(1) : pFD->getVert(0);

      Face *pFA = dynamic_cast<TriCell*>(pCA)->getOppositeFace(pV);
      Face *pFB = dynamic_cast<TriCell*>(pCB)->getOppositeFace(pV);
      Face *pFC1 = dynamic_cast<TriCell*>(pCB)->getOppositeFace(pVC);
      Face *pFC2 = dynamic_cast<TriCell*>(pCA)->getOppositeFace(pVC);

      assert(pFC1->getOppositeCell(pCB)->getType() == Cell::eBdryEdge);
      assert(pFC2->getOppositeCell(pCA)->getType() == Cell::eBdryEdge);

      BdryEdge *pBE1 = dynamic_cast<BdryEdge *>(pFC1->getOppositeCell(pCB));
      BdryEdge *pBE2 = dynamic_cast<BdryEdge *>(pFC2->getOppositeCell(pCA));

      // Originally:
      //  Cell pCA: pV, pVC, pVB; pFD, pFA, pFC2
      //  Cell pCB: pV, pVC, pVA; pFD, pFB, pFC1
      //  BFace pBE1: pFC1
      //  BFace pBE2: pFC2
      // After:
      //  Delete pFD, pFC1, pFC2, pCA, pCB, pBE1, pBE2
      //  Create pFC(pVA, pVB); pC(pFA, pFB, pFC); pBE(pFC)

      // Only delete the vertex if the boundary conditions are the same
      if (pBE1->getBdryCondition() == pBE2->getBdryCondition()) {
	int iReg = pCA->getRegion();
	bool qExist;
	deleteCell(pCA);
	deleteCell(pCB);

	// Faces are automatically deleted when their last cell is deleted.

	Face *pFC = createFace(qExist, pVA, pVB);
	assert(!qExist);
	Cell *pC = createTriCell(pFA, pFB, pFC, iReg);
	// Create a new BFace based on the old one
	BFace *pBF = createBFace(pFC, pBE1);

	// Now it's safe to delete the old BFaces.
	deleteBFace(pBE1);
	deleteBFace(pBE2);

	// Delete this last, so that its faces are all gone already, so
	// that there won't be some silent refusal to delete an entity
	// in use.
	deleteVert(pV);
#ifndef OMIT_VERTEX_HINTS
	pVC->setHintFace(pFB);
#endif
	assert(pBF->doFullCheck());
	assert(pC->doFullCheck());
	assert(pFA->doFullCheck());
	assert(pFB->doFullCheck());
	assert(pFC->doFullCheck());

	if (pCellQueue != NULL) {
	  pCellQueue->vAddCell(pC);
	}

      }
      else
	qRetVal = false;
    } // Done with three edges incident on a bdry vert
    else if (iNEdges == 2) {
      // Two edges incident on a bdry vert
      // Identify all the parts
      assert(spCIncCells.size() == 3);
      Cell* apCIncCells[3];
      std::copy(spCIncCells.begin(), spCIncCells.end(), apCIncCells);
      Cell *pCA;
      if (apCIncCells[0]->getType() == Cell::eTriCell) {
	pCA = apCIncCells[0];
      }
      else if (apCIncCells[1]->getType() == Cell::eTriCell) {
	pCA = apCIncCells[1];
      }
      else {
	pCA = apCIncCells[2];
      }
      assert(pCA->getType() == Cell::eTriCell);

      if (pCellQueue != NULL) {
	pCellQueue->vRemoveCell(pCA);
      }

      Face *pF0 = vecpFIncEdges[0];
      Face *pF1 = vecpFIncEdges[1];

      Face *pFA = dynamic_cast<TriCell*>(pCA)->getOppositeFace(pV);

      assert(pF0->getOppositeCell(pCA)->getType() == Cell::eBdryEdge);
      assert(pF1->getOppositeCell(pCA)->getType() == Cell::eBdryEdge);

      BdryEdge *pBE0 = dynamic_cast<BdryEdge *>(pF0->getOppositeCell(pCA));
      BdryEdge *pBE1 = dynamic_cast<BdryEdge *>(pF1->getOppositeCell(pCA));

      // Only delete the vertex if the boundary conditions are the same
      if (pBE0->getBdryCondition() == pBE1->getBdryCondition()) {
	deleteCell(pCA);
	deleteVert(pV);

	BFace* pBFNew = createBFace(pFA, pBE0);
	deleteBFace(pBE0);
	deleteBFace(pBE1);
	assert(pBFNew->doFullCheck());
	assert(pFA->doFullCheck());
	if (pCellQueue != NULL) {
	  pCellQueue->vAddCell(pCA);
	}
      }
      else {
	printf("Boundary conditions differ\n");
	//	assert(0);
	qRetVal = false;
      }
    }
    else
      qRetVal = false; // Too many edges incident on a bdry vert? not sure what to do!
  }				// Done deleting boundary vert
  for (iterF = spFSwapEdges.begin(); iterF != spFSwapEdges.end(); ++iterF) {
    Face *pF = *iterF;
    if (pF->isValid() && !pF->isDeleted()) {
      if (pAnisoRef) {
	iSwaps = iSwaps + (pAnisoRef->iAniSwap(pF));
      }
      else {
	iSwaps = iSwaps + iFaceSwap_deprecated(pF);
      }
    }
  }

  return qRetVal;
}
