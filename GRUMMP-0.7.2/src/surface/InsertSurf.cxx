#include "GR_assert.h"
#include "GR_Geometry.h"
#include "GR_SurfMesh.h"
#include "GR_Vec.h"
#include "GR_Bdry3D.h"

//@ Insert a point on an interior face
int
SurfMesh::insertOnFace(Vert* pVNew, Face* const pF, Cell* const pC,
		       const bool qSwap, bool* qInsertOnEdge)
// Note: Geometrically, pVNew may not fall exactly on pF because
// of round-off.  Once things are done, though, the connectivity will
// be the same as if it had.
//
// Before:                 After:
//        B                      B
//       /|\		          /|\       .
//	  / | \		         / F \      .
//  FA1/  |  \FA2	     FA1/  M  \FA2
//    /   |	\	       /N3 B N2\    .
//   /    |	 \	      /    |    \   .
//  C	 C1 | C2  D	     C-FN1---FN2-D
//   \    |	 /	      \    |    /
//    \   |	/	       \N4 F N1/
//  FB1\  |  /FB2	     FB1\  M  /FB2
//	  \ | /		         \ A /
//       \|/		          \|/
//        A	           A
//
{
  assert(pC->doFullCheck());
  assert(pF->doFullCheck());

  if (pF->getType() == Face::eMultiEdge) {
    if (qInsertOnEdge)
      *qInsertOnEdge = true;
    return insertOnMultiEdge(pVNew, dynamic_cast<MultiEdge*>(pF), qSwap);
  }

  Cell* pC1 = pF->getLeftCell();
  Cell* pC2 = pF->getRightCell();
  assert(pC1->doFullCheck());
  assert(pC2->doFullCheck());
  assert(pC1 == pC || pC2 == pC);
  BdryPatch3D *pBP1 = getBdryPatch(getCellIndex(pC1));
  BdryPatch3D *pBP2 = getBdryPatch(getCellIndex(pC2));
  if (qInsertOnEdge)
    *qInsertOnEdge = pBP1 != pBP2;

  assert(pC1->getType() == Cell::eTriCell);
  assert(pC2->getType() == Cell::eTriCell);

  TriCell* const pTC1 = dynamic_cast<TriCell*>(pC1);
  TriCell* const pTC2 = dynamic_cast<TriCell*>(pC2);

  Vert *pVA = pF->getVert(0);
  Vert *pVB = pF->getVert(1);
  Vert *pVC = pTC1->getOppositeVert(pF);
  Vert *pVD = pTC2->getOppositeVert(pF);

  // Orientation for surface is a bit tricky
//   assert(iOrient2D(pVA, pVB, pVC) == 1);
//   assert(iOrient2D(pVB, pVA, pVD) == 1);

// Delete the old cells and face
  deleteCell(pTC1);
  deleteCell(pTC2);
  deleteFace(pF);

  Cell *pCNew1 = createTriCell(pVNew, pVA, pVD);
  Cell *pCNew2 = createTriCell(pVNew, pVD, pVB);
  Cell *pCNew3 = createTriCell(pVNew, pVB, pVC);
  Cell *pCNew4 = createTriCell(pVNew, pVC, pVA);

  // Orientation for surface is a bit tricky
//   assert(iOrient2D(pVA, pVD, pVNew) == 1);
//   assert(iOrient2D(pVD, pVB, pVNew) == 1);
//   assert(iOrient2D(pVB, pVC, pVNew) == 1);
//   assert(iOrient2D(pVC, pVA, pVNew) == 1);

  setBdryPatch(pCNew1, pBP2);
  setBdryPatch(pCNew2, pBP2);
  setBdryPatch(pCNew3, pBP1);
  setBdryPatch(pCNew4, pBP1);

  //@@ Swap locally
  // Swapping is not yet defined for surface meshes
  int iRetVal = 0;
  if (qSwap) {
    Face *pFA1 = pCNew3->getOppositeFace(pVNew);
    iRetVal += iFaceSwap_deprecated(pFA1);
    Face *pFA2 = pCNew2->getOppositeFace(pVNew);
    iRetVal += iFaceSwap_deprecated(pFA2);
    Face *pFB1 = pCNew4->getOppositeFace(pVNew);
    iRetVal += iFaceSwap_deprecated(pFB1);
    Face *pFB2 = pCNew1->getOppositeFace(pVNew);
    iRetVal += iFaceSwap_deprecated(pFB2);
  }
  return (iRetVal);
}

//@ Insert a point on a MultiEdge
int
SurfMesh::insertOnMultiEdge(Vert* const pVNew, MultiEdge* const pME,
			    const bool /*qSwap*/)
// Note: Geometrically, pVNew may not fall exactly on pF because
// of round-off.  Once things are done, though, the connectivity will
// be the same as if it had. Each old cell incident on the MultiEdge
// gets split in the same way.
//
// Before:                 After:
//        B                      B
//       /|		          /|
//	  / |		         / M
//  FA1/  |	             FA1/  E
//    /   |		       /C2 2
//   /    |		      /    |
//  C	 C  M	             C-FN---
//   \    E		      \    |
//    \   |		       \C1 M
//  FB1\  |	             FB1\  E
//	  \ |		         \ 1
//       \|		          \|
//        A	           A
//
{
  assert(pME->doFullCheck());

  Vert *pVA = pME->getVert(0);
  Vert *pVB = pME->getVert(1);

  MultiEdge *pME1 = createMultiEdge(pVA, pVNew);
  MultiEdge *pME2 = createMultiEdge(pVNew, pVB);

  int iNCells = pME->getNumCells();

  // All the faces that will need to be checked for swapping.
  Face **apFCheck = new Face*[iNCells * 2];
  int iC;
  for (iC = iNCells - 1; iC >= 0; iC--) {
    Cell *pC = pME->getCell(iC);
    assert(pC->doFullCheck());
    BdryPatch3D *pBP = getBdryPatch(getCellIndex(pC));

    assert(pC->getType() == Cell::eTriCell);
    TriCell* const pTC = dynamic_cast<TriCell*>(pC);

    Vert *pVC = pTC->getOppositeVert(pME);

    Face* pFA1 = pTC->getOppositeFace(pVA);
    Face* pFB1 = pTC->getOppositeFace(pVB);

    apFCheck[iC * 2] = pFA1;
    apFCheck[iC * 2 + 1] = pFB1;

    //@@ Create new cells
    int iReg = pC->getRegion();
    deleteCell(pC);

    Cell *pC1 = createTriCell(pVA, pVNew, pVC, iReg);
    Cell *pC2 = createTriCell(pVNew, pVB, pVC, iReg);

    setBdryPatch(pC1, pBP);
    setBdryPatch(pC2, pBP);

    assert(pFA1->doFullCheck() == 1);
    assert(pFB1->doFullCheck() == 1);
    assert(pC1->doFullCheck() == 1);
    assert(pC2->doFullCheck() == 1);

  }
  // Update the two MultiEdge's
  assert(pME1->getNumCells() == pME2->getNumCells());
  assert(pME1->getNumCells() == iNCells);
  assert(pME1->doFullCheck());
  assert(pME2->doFullCheck());

  //@@ Swap locally
  // Swapping is not yet defined for surface meshes
  int iRetVal = 0;
//  if (qSwap) {
//    for (int i = 0; i < iNCells*2; i++)
//      iRetVal += iFaceSwap_deprecated(apFCheck[i]);
//  }
  delete[] apFCheck;
  return (iRetVal);
}

// //@ Insert a point inside a cell
// int SurfMesh::iInsertInInterior(Vert* const pVNew, Cell* const pC,
// 				const bool qSwap)
// {
//   assert(pC->iFullCheck());
//   assert(pC->eType() == Cell::eTriCell);
//   TriCell* pTC = dynamic_cast<TriCell*>(pC);

//   // Verts ABC are ordered in a right-handed way.  The faces are
//   // defined relative to them so that there is no possible ambiguity.
//   // Note that the face definitions here DO NOT correspond to
//   // iFace[ABCD](iCell).
//   Vert* pVA = pC->pVVert(0);
//   Vert* pVB = pC->pVVert(1);
//   Vert* pVC = pC->pVVert(2);

//   Face* pFA = pTC->pFFaceOpposite(pVA);
//   Face* pFB = pTC->pFFaceOpposite(pVB);
//   Face* pFC = pTC->pFFaceOpposite(pVC);

//   // Identify the bdry patch and region
//   BdryPatch3D *pBP = pBdryPatch(iCellIndex(pC));
//   int iReg = pC->iRegion();

//   // Get rid of the old cell.
//   deleteCell(pC);

//   // Get new cells
//   Cell* pCA = createTriCell(pVA, pVB, pVNew, iReg);
//   Cell* pCB = createTriCell(pVB, pVC, pVNew, iReg);
//   Cell* pCC = createTriCell(pVC, pVA, pVNew, iReg);

//   assert(pCA ->iFullCheck());
//   assert(pCB ->iFullCheck());
//   assert(pFC ->iFullCheck());
//   assert(pFA ->iFullCheck());
//   assert(pFB ->iFullCheck());
//   assert(pFC ->iFullCheck());

//   vSetBdryPatch(pCA, pBP);
//   vSetBdryPatch(pCB, pBP);
//   vSetBdryPatch(pCC, pBP);

//   // Do some local swapping to clean up
//   // Swapping not yet defined for surface meshes
//   int iRetVal = 0;
//   if (qSwap) {
//     // Swapping for bdry faces not yet implemented
//     iRetVal += iFaceSwap(pFA);
//     iRetVal += iFaceSwap(pFB);
//     iRetVal += iFaceSwap(pFC);
//   }
//   return (iRetVal);
// }

// //@ Add a site (or some variant if it lies outside the mesh) to the mesh
// // Return value indicates success
// bool SurfMesh::qInsertPoint(const double adPoint[3], Cell* const pC,
// 			    int * const piSwaps, const bool qSwap,
// 			    const bool qForce, Vert* pVNew)
// {
//   *piSwaps = 0;
//   //@@ The given cell had better contain the site
//   // Project the cell and site onto a plane using the cell normal.
//   double adNorm[3];
//   pC->vVecSize(adNorm);
//   NORMALIZE3D(adNorm);
// #ifndef NDEBUG
//   {
//     double dFaceLoc = dDOT3D(adNorm, pC->pVVert(0)->adCoords());
//     double dSiteLoc = dDOT3D(adNorm, adPoint);
//     assert(iFuzzyComp(dFaceLoc, dSiteLoc) == 0);
//   }
// #endif
//   if (pVNew->qValid()) {
//     assert(qForce); // Can't twiddle the point location if it's an
// 		    // existing vertex.
//     assert(1.e-12 > dDIST3D(pVNew->adCoords(), adPoint));
//   }

//   // Coordinates in the plane of the triangular face
//   Vert *pV0 = pC->pVVert(0);
//   Vert *pV1 = pC->pVVert(1);
//   Vert *pV2 = pC->pVVert(2);
//   double adX[] = adDIFF3D(pV0->adCoords(), pV1->adCoords());
//   NORMALIZE3D(adX);
//   double adY[3];
//   vCROSS3D(adNorm, adX, adY);
//   NORMALIZE3D(adY);

//   // Now project onto a plane
//   double adLoc0[] = {dDOT3D(pV0->adCoords(), adX),
// 		     dDOT3D(pV0->adCoords(), adY)};
//   double adLoc1[] = {dDOT3D(pV1->adCoords(), adX),
// 		     dDOT3D(pV1->adCoords(), adY)};
//   double adLoc2[] = {dDOT3D(pV2->adCoords(), adX),
// 		     dDOT3D(pV2->adCoords(), adY)};
//   double adLocSite[] = {dDOT3D(adPoint, adX), dDOT3D(adPoint, adY)};

//   int iOrient01 = iOrient2D(adLoc0, adLoc1, adLocSite);
//   int iOrient12 = iOrient2D(adLoc1, adLoc2, adLocSite);
//   int iOrient20 = iOrient2D(adLoc2, adLoc0, adLocSite);
//   // Make sure the inserted site lies inside or on the boundary of the
//   // cell
//   if (iOrient01*iOrient12 == -1 ||
//       iOrient01*iOrient20 == -1 ||
//       iOrient12*iOrient20 == -1) {
//     vWarning("Point not inserted in surface mesh.  Lies outside cell.");
//     return false;
//   }

//   // If the site lies strictly inside, insert it inside
//   if (iOrient01*iOrient12*iOrient20 != 0) {
//     iInterior++;
//     if (!pVNew->qValid())
//       pVNew = createVert(adPoint);
//     else
//       pVNew->vSetCoords(3, adPoint);
//     *piSwaps = iInsertInInterior(pVNew, pC, qSwap);
//     return true;
//   }
//   else {
//     // Find the face on which the site lies
//     Face *pF = pFInvalidFace;
//     for (int iF = 0; iF < 3; iF++) {
//       pF = pC->pFFace(iF);
//       int qHas0 = pF->qHasVert(pV0);
//       int qHas1 = pF->qHasVert(pV1);
//       if (iOrient01 == 0 && qHas0 && qHas1) break;
//       int qHas2 = pF->qHasVert(pV2);
//       if (iOrient12 == 0 && qHas1 && qHas2) break;
//       if (iOrient20 == 0 && qHas2 && qHas0) break;
//     }
//     if (!pVNew->qValid())
//       pVNew = createVert(adPoint);
//     else
//       pVNew->vSetCoords(3, adPoint);
//     *piSwaps = iInsertOnFace(pVNew, pF, pC, qSwap);
//     return true;
//   }
// }

//@@ Determines whether swapping is needed for local Delaunay-ness.
static bool
qDoSwapDelaunay(const Vert * pVVertA, const Vert * pVVertB,
		const Vert * pVVertC, const Vert * pVVertD)
// Is pVVertD inside the diametral sphere of pVVertA, pVVertB, pVVertC?
{
  // Find the circumcenter (in-plane) of ABC by finding intersection of
  // the perpendicular bisectors of AB and BC with the plane ABC.
  double adVecAB[] = adDIFF3D(pVVertB->getCoords(), pVVertA->getCoords());
  double adVecBC[] = adDIFF3D(pVVertC->getCoords(), pVVertB->getCoords());

  NORMALIZE3D(adVecAB);
  NORMALIZE3D(adVecBC);

  double adNorm[3];
  vCROSS3D(adVecAB, adVecBC, adNorm);
  NORMALIZE3D(adNorm);

  // The components of each of these vectors, dotted with coordinate
  // vectors from the origin, give level sets of distance parallel to
  // the vector.  Calibrate for the level set we need.
  double adRHS[3];
  {
    double adMidPtAB[] =
      { 0.5 * (pVVertA->x() + pVVertB->x()), 0.5
	  * (pVVertA->y() + pVVertB->y()), 0.5 * (pVVertA->z() + pVVertB->z()) };
    double adMidPtBC[] =
      { 0.5 * (pVVertB->x() + pVVertC->x()), 0.5
	  * (pVVertB->y() + pVVertC->y()), 0.5 * (pVVertB->z() + pVVertC->z()) };
    adRHS[0] = dDOT3D(adVecAB, adMidPtAB);
    adRHS[1] = dDOT3D(adVecBC, adMidPtBC);
    adRHS[2] = dDOT3D(adNorm, pVVertA->getCoords());
  }
  double adCircCent[3];
  solve3By3(adVecAB, adVecBC, adNorm, adRHS[0], adRHS[1], adRHS[2], adCircCent);
#ifndef NDEBUG
  {
    // Make sure we found the right point.
    double dDistA = dDIST3D(adCircCent, pVVertA->getCoords());
    double dDistB = dDIST3D(adCircCent, pVVertB->getCoords());
    double dDistC = dDIST3D(adCircCent, pVVertC->getCoords());
    assert(iFuzzyComp(dDistA, dDistB) == 0);
    assert(iFuzzyComp(dDistA, dDistC) == 0);
  }
#endif
  double dRadius = dDIST3D(adCircCent, pVVertA->getCoords());
  double dDistance = dDIST3D(adCircCent, pVVertD->getCoords());

  switch (iFuzzyComp(dRadius, dDistance))
    {
    case -1:
      return false;		// pVVertD is outside
    case 1:
      return true;		// pVVertD is inside
    case 0:
      // Gives a consistent answer, regardless of which order the
      // configuration is presented in.  (ABCD, BADC, CDAB, DCBA all
      // result in the same triangulation.)
      return (min(pVVertA, pVVertB) > min(pVVertC, pVVertD));
    default:
      assert(0);
      return false;
    }
}

//@@ Determines whether swapping will improve the maximum face angle.
static bool
qDoSwapMaxDihed(const Vert * pVVertA, const Vert * pVVertB,
		const Vert * pVVertC, const Vert * pVVertD)
{
  // This algorithm finds face angles between adjacent faces by dotting
  // unit vectors along the faces.  Small magnitude (close to -1)
  // indicates edges pointing in opposite directions, and this is bad,
  // so we choose the configuration that gives the largest min dot
  // product.
  //
  // To prevent pairs from flopping back and forth, the tie-breaker is
  // invoked if the inequality is small.

  // All of these normals are intended to point inwards.
  double adVecDA[] = adDIFF3D(pVVertD->getCoords(), pVVertA->getCoords());
  double adVecBD[] = adDIFF3D(pVVertB->getCoords(), pVVertD->getCoords());
  double adVecCB[] = adDIFF3D(pVVertC->getCoords(), pVVertB->getCoords());
  double adVecAC[] = adDIFF3D(pVVertA->getCoords(), pVVertC->getCoords());

  NORMALIZE3D(adVecDA);
  NORMALIZE3D(adVecBD);
  NORMALIZE3D(adVecCB);
  NORMALIZE3D(adVecAC);

  double dDotBDA = -dDOT3D(adVecBD, adVecDA);
  double dDotACB = -dDOT3D(adVecAC, adVecCB);
  double dDotDAC = -dDOT3D(adVecDA, adVecAC);
  double dDotCBD = -dDOT3D(adVecCB, adVecBD);

  double dMinThis = min(dDotBDA, dDotACB);
  double dMinOther = min(dDotDAC, dDotCBD);

  switch (iFuzzyComp(dMinThis, dMinOther))
    {
    case -1:
      return true;
    case 0:
      // Gives a consistent answer, regardless of which order the
      // configuration is presented in.
      return (min(pVVertA, pVVertB) > min(pVVertC, pVVertD));
    case 1:
      return false;
    default:
      assert(0);
      return false;
    }
}

bool
SurfMesh::doSwap(double /*dLenDiag*/, double /*dLenCL1*/, double /*dLenCL2*/,
		 double /*dLenCR1*/, double /*dLenCR2*/)
{
  assert(0);
  return false;
}

bool
SurfMesh::doSwap(const Vert * const pVVertA, const Vert * const pVVertB,
		 const Vert * const pVVertC, const Vert * const pVVertD,
		 const Vert * const pVVertE) const
{
  assert(pVVertA->isValid());
  assert(pVVertB->isValid());
  assert(pVVertC->isValid());
  assert(pVVertD->isValid());
  assert(!pVVertE->isValid());

  switch (getSwapType())
    {
    case eDelaunay:
      return qDoSwapDelaunay(pVVertA, pVVertB, pVVertC, pVVertD);
    case eMaxDihed:
      return qDoSwapMaxDihed(pVVertA, pVVertB, pVVertC, pVVertD);
    default:
      return (false);
    }
}

int
SurfMesh::iFaceSwap_deprecated(Face*& pF)
{
  assert(getSwapType() != eNoSwap);
  if (!pF->doFullCheck())
    return (0);

  // Never swap a face that divides two regions in a multi-region mesh.
  if (pF->getFaceLoc() == Face::eBdryTwoSide)
    return (0);

  Cell *pCLeft = pF->getLeftCell();
  Cell *pCRight = pF->getRightCell();
  if (!pCLeft->isValid() || !pCRight->isValid()
      || pCLeft->getType() != Cell::eTriCell
      || pCRight->getType() != Cell::eTriCell
      || (getBdryPatch(getCellIndex(pCLeft))
	  != getBdryPatch(getCellIndex(pCRight)))) {
    return (0);
  }

  Vert *pVVertA = pF->getVert(0);
  Vert *pVVertB = pF->getVert(1);
  Vert *pVVertC = pCLeft->getOppositeVert(pF);
  Vert *pVVertD = pCRight->getOppositeVert(pF);

  // The orientation check that the 2D code does, can't be done here.

  if (doSwap(pVVertA, pVVertB, pVVertC, pVVertD)) {
    return reconfigure(pF);
  }
  else
    return (0);
}

int
SurfMesh::reconfigure(Face*& pFEdge)
{
  if (!pFEdge->isSwapAllowed())
    return 0;
//   vMessage(2, "Swapping surface edge %d %d\n",
// 	   iVertIndex(pFEdge->pVVert(0)),
// 	   iVertIndex(pFEdge->pVVert(1)));
  assert(pFEdge->doFullCheck());
  assert(!((m_swapMeasure == eNoSwap) && m_swapRecursively));
  if (pFEdge->getFaceLoc() == Face::eBdryTwoSide)
    return (0);

  // Before:     1	After:      1
  //            / \ 	           /|\  	.
  //           C 1 D 	          C | D
  //          /     \	         /  |  \        .
  //         B-------A	        B 1 | 0 A
  //          \     /	         \  |  /
  //           B 0 A	          B | A
  //            \ /	           \|/
  //             0                  0
  //
  //
  //

  // Identify all the parts
  TriCell *pTri0 = dynamic_cast<TriCell*>(pFEdge->getLeftCell());
  TriCell *pTri1 = dynamic_cast<TriCell*>(pFEdge->getRightCell());
  Vert *pVA = pFEdge->getVert(0);
  Vert *pVB = pFEdge->getVert(1);
  Vert *pV0 = pTri0->getOppositeVert(pFEdge);
  Vert *pV1 = pTri1->getOppositeVert(pFEdge);

  Face *pFA = pTri0->getOppositeFace(pVB);
  Face *pFB = pTri0->getOppositeFace(pVA);
  Face *pFC = pTri1->getOppositeFace(pVA);
  Face *pFD = pTri1->getOppositeFace(pVB);

  BdryPatch3D* pBP = getBdryPatch(pTri0);
  assert(getBdryPatch(pTri1) == pBP);

  deleteCell(pTri0);
  deleteCell(pTri1);
  deleteFace(pFEdge);

  // Create the edge first so we can prevent non-manifold edges being
  // first in a cell's list.
  bool qExist;
  pFEdge = createFace(qExist, pV0, pV1);
  assert(!qExist);

  Cell *pC0 = createTriCell(pFEdge, pFA, pFD);
  Cell *pC1 = createTriCell(pFEdge, pFC, pFB);
  assert(pC0->getFace(0) == pFEdge);
  assert(pC1->getFace(0) == pFEdge);

  setBdryPatch(pC0, pBP);
  setBdryPatch(pC1, pBP);

  pFEdge = findCommonFace(pC0, pC1);

#ifndef OMIT_VERTEX_HINTS
  pVA->setHintFace(pFA);
  pVB->setHintFace(pFB);
#endif
  assert(pC0->doFullCheck());
  assert(pC1->doFullCheck());

  assert(pFA->doFullCheck());
  assert(pFB->doFullCheck());
  assert(pFC->doFullCheck());
  assert(pFD->doFullCheck());

  int iRetVal = 1;
  if (m_swapRecursively) {
    iRetVal += iFaceSwap_deprecated(pFA);
    iRetVal += iFaceSwap_deprecated(pFB);
    iRetVal += iFaceSwap_deprecated(pFC);
    iRetVal += iFaceSwap_deprecated(pFD);
  }
  return iRetVal;

  // This clause should be irrelevant, because all tris should have
  // non-manifold edges as their first edge.
//   if (iRetVal != 0) {
//     // Check to make sure that neither cell incident on pF now has a
//     // MultiEdge as its first face, as this can cause orientation
//     // problems.
//     Cell *pCL = pFEdge->pCCellLeft();
//     Face *pF0 = pCL->pFFace(0);
//     if (pF0->eType() == Face::eMultiEdge) {
//       vMessage(2, "Cleaning up after swap (verts %d %d)\n",
// 	       iVertIndex(pF0->pVVert(0)),
// 	       iVertIndex(pF0->pVVert(1)));
//       Face *pF1 = pCL->pFFace(1);
//       Face *pF2 = pCL->pFFace(2);
//       if (pF1->eType() == Face::eEdgeFace)
// 	pCL->vAssign(pF1, pF2, pF0);
//       else {
// 	assert(pF2->eType() == Face::eEdgeFace);
// 	pCL->vAssign(pF2, pF0, pF1);
//       }
//     }
//     Cell *pCR = pFEdge->pCCellRight();
//     pF0 = pCR->pFFace(0);
//     if (pF0->eType() == Face::eMultiEdge) {
//       vMessage(2, "Cleaning up after swap (verts %d %d)\n",
// 	       iVertIndex(pF0->pVVert(0)),
// 	       iVertIndex(pF0->pVVert(1)));
//       Face *pF1 = pCR->pFFace(1);
//       Face *pF2 = pCR->pFFace(2);
//       if (pF1->eType() == Face::eEdgeFace)
// 	pCR->vAssign(pF1, pF2, pF0);
//       else {
// 	assert(pF2->eType() == Face::eEdgeFace);
// 	pCR->vAssign(pF2, pF0, pF1);
//       }
//     }
//     return (iRetVal);
//   }
//   else {
//     assert(iRetVal == 0);
//     return 0;
//   }
}
