#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <stack>

#include "GR_config.h"
#include "GR_Geometry.h"
#include "GR_InsertionManager.h"
#include "GR_SurfMesh.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_VolMesh.h"
#include "GR_Constrain.h"
#include "GR_Bdry3D.h"

#define SMALL_ANGLE_RECOVERY

static GRUMMP::SwappingInserter3D* pSI3D = NULL;

typedef struct _SplitEdge_ {
  Vert *pV0, *pV1, *pVVolNew, *pVSurfNew;
  bool qOnBdryCurve;
  _SplitEdge_(Vert *pV0_arg, Vert *pV1_arg, Vert *pVVolNew_arg,
	      Vert *pVSurfNew_arg, bool qOnCurve) :
      pV0(pV0_arg), pV1(pV1_arg), pVVolNew(pVVolNew_arg), pVSurfNew(
	  pVSurfNew_arg), qOnBdryCurve(qOnCurve)
  {
  }
} SplitEdge;

static std::stack<SplitEdge> stackSplitEdges;

//static const int iMaxSplittings = 12;

static bool
qRemoveBdrySteinerPoint(const SplitEdge& SE, VolMesh& VM, int& iSwaps,
			bool& qDelInterior)
{
  Vert *pV = SE.pVVolNew;
  assert(pV->isValid() && pV->isDeletionRequested() && !pV->isDeleted());
  // Grab the vertex and its neighborhood.
  std::set<Vert*> spVNear;
  std::set<Face*> spFNear;
  std::set<Cell*> spCInc;
  std::set<BFace*> spBFaceInc;
  bool qBdry;
  findNeighborhoodInfo(pV, spCInc, spVNear, &spBFaceInc, &qBdry, &spFNear);

  if (spBFaceInc.size() > 4) {
    logMessage(
	3, "Possible panic!  Hit %zd incident bfaces; unsure how to proceed.\n",
	spBFaceInc.size());
  }

  logMessage(1, "Trying to remove bdry Steiner point %d at (%f, %f, %f)\n",
	     VM.getVertIndex(pV), pV->x(), pV->y(), pV->z());
  std::set<BFace*>::iterator bfaceIter = spBFaceInc.begin(), bfaceIterEnd =
      spBFaceInc.end();
  double adDir[] =
    { 0, 0, 0 };
  for (; bfaceIter != bfaceIterEnd; bfaceIter++) {
    BFace *pBF = *bfaceIter;
    double adNorm[3];
    // This provides a non-normalized normal vector, pointing into the
    // domain.  Things could be tricky for interior bdry faces.
    pBF->calcVecSize(adNorm);
    adDir[0] += adNorm[0];
    adDir[1] += adNorm[1];
    adDir[2] += adNorm[2];
  }
  NORMALIZE3D(adDir);
  logMessage(3, "  Normal: (%f, %f, %f)\n", adDir[0], adDir[1], adDir[2]);

  std::set<Face*> spFHull;
  std::set<Face*>::iterator faceIter = spFNear.begin(), faceIterEnd =
      spFNear.end();
  double dDistToIntersection = 1.e300;

  for (; faceIter != faceIterEnd; faceIter++) {
    Face *pF = *faceIter;
    if (!pF->hasVert(pV)) {
      spFHull.insert(pF);
      Vert *pV0 = pF->getVert(0);
      Vert *pV1 = pF->getVert(1);
      Vert *pV2 = pF->getVert(2);
      // Outward normal of this face.
      double adNorm[3];
      calcNormal3D(pV0->getCoords(), pV1->getCoords(), pV2->getCoords(),
		   adNorm);
      NORMALIZE3D(adNorm);
      logMessage(3, "  Hull face %d; verts %d, %d, %d\n", VM.getFaceIndex(pF),
		 VM.getVertIndex(pV0), VM.getVertIndex(pV1),
		 VM.getVertIndex(pV2));
      logMessage(3, "    Vert 0:  (%f, %f, %f)\n", pV0->x(), pV0->y(),
		 pV0->z());
      logMessage(3, "    Vert 1:  (%f, %f, %f)\n", pV1->x(), pV1->y(),
		 pV1->z());
      logMessage(3, "    Vert 2:  (%f, %f, %f)\n", pV2->x(), pV2->y(),
		 pV2->z());

      // Distance from face to this point.
      double adDiff[] = adDIFF3D(pV0->getCoords(), pV->getCoords());
      double dDistToPlane = dDOT3D(adDiff, adNorm);
      if (dDistToPlane < 0) {
	dDistToPlane *= -1;
	vSCALE3D(adNorm, -1);
      }
      logMessage(3, "  Face normal: (%f, %f, %f); dist: %f\n", adNorm[0],
		 adNorm[1], adNorm[2], dDistToPlane);

      // Get the cos of the angle between the normal from the surface
      // and normal from this face.
      double dCos = dDOT3D(adDir, adNorm) / dMAG3D(adNorm);
      logMessage(3, "  Cos of angle: %f; angle %f degrees\n", dCos,
		 GR_acos(dCos) * 180 / M_PI);

      // Find distance to this intersection.
      double dThisDist = 1.e300;
      if (dCos > 0)
	dThisDist = dDistToPlane / dCos;
      if (dThisDist < dDistToIntersection)
	dDistToIntersection = dThisDist;

      logMessage(4, "  This dist: %f; smallest so far: %f\n", dThisDist,
		 dDistToIntersection);
    } // Done processing outer hull verts.
  } // Done checking all nearby faces.

  logMessage(4, "  Dist to interior Steiner: %f\n", dDistToIntersection);
  if (dDistToIntersection < 1.e-12)
    return false;
  double adNewLoc[] =
    { pV->x() + adDir[0] * dDistToIntersection / 2, pV->y()
	+ adDir[1] * dDistToIntersection / 2, pV->z()
	+ adDir[2] * dDistToIntersection / 2 };
  logMessage(1, "  inserting an interior Steiner point at (%f, %f, %f)\n",
	     adNewLoc[0], adNewLoc[1], adNewLoc[2]);

  // Check to see whether this point lies inside one of the tets.
  GRUMMP::MinMaxDihedSwapDecider3D Del3D(true);
  GRUMMP::SwapManager3D DelSwap(&Del3D, &VM);
  pSI3D = new GRUMMP::SwappingInserter3D(&VM, &DelSwap);
  pSI3D->setForcedInsertion(true);
  Vert *pVIntSteiner = pSI3D->insertPoint(adNewLoc, *spCInc.begin());
  //  Vert *pVIntSteiner = VM.createVert(adNewLoc);
  pVIntSteiner->markToDelete();
  delete pSI3D;
  // Now need to delete all the incident tets; their shared faces go
  // with them with automagically.
  //  std::set<Cell*>::iterator cellIter;
  //  const std::set<Cell*>::iterator cellIterEnd = spCInc.end();
  //  for (cellIter = spCInc.begin(); cellIter != cellIterEnd; cellIter++) {
  //    Cell *pC = *cellIter;
  //    VM.deleteCell(pC);
  //  }
  //
  //  // Create new tets for each hull face
  //  faceIterEnd = spFHull.end();
  //  for (faceIter = spFHull.begin(); faceIter != faceIterEnd; faceIter++) {
  //    bool qExist;
  //    VM.createTetCell(qExist, pVIntSteiner, *faceIter);
  //    assert(!qExist);
  //  }
  //
  //  // And now one for each of the faces on the bdry
  //  bfaceIterEnd = spBFaceInc.end();
  //  for (bfaceIter = spBFaceInc.begin(); bfaceIter != bfaceIterEnd;
  //       bfaceIter++) {
  //    Face *pF = (*bfaceIter)->getFace();
  //    bool qExist;
  //    VM.createTetCell(qExist, pVIntSteiner, pF);
  //    assert(!qExist);
  //  }

  bool qRes;
  int iNewSwaps = 0;
  Vert* apVTargets[] =
    { SE.pV0, SE.pV1 };
  qRes = VM.removeVert(pV, iNewSwaps, apVTargets, 2);
  iSwaps += iNewSwaps;
  if (!qRes) {
    logMessage(1, "  Bdry removal failed!\n");
    //    return false;
  }
  //  assert(qRes);
  if (qRes)
    logMessage(1, "  Bdry removal succeeded with %d swaps.\n", iNewSwaps);

  qDelInterior = VM.removeVert(pVIntSteiner, iNewSwaps);
  iSwaps += iNewSwaps;
  logMessage(1, "  Interior removal %s with %d swaps.\n",
	     qDelInterior ? "succeeded" : "failed", iNewSwaps);
  return true;
}

static void
vComputeEdgeAndTriStatus(const SurfMesh& SM, std::set<ConstrainEdge> sCESubSegs,
			 std::set<ConstrainEdge> sCEInFacets,
			 std::set<ConstrainFace> sCFSubFacets,
			 int aiTriStatus[], int aiEdgeStatus[])
{
  for (GR_index_t iC = 0; iC < SM.getNumCells(); iC++) {
    aiTriStatus[iC] = 0;
  }

  for (GR_index_t iF = 0; iF < SM.getNumFaces(); iF++) {
    aiEdgeStatus[iF] = 0;
  }

  logMessage(2, "Missing %zd (sub)facets\n", sCFSubFacets.size());
  for (std::set<ConstrainFace>::iterator iter = sCFSubFacets.begin();
      iter != sCFSubFacets.end(); ++iter) {
    ConstrainFace CF = *iter;
    Vert *pV0 = SM.getVert(CF.iV0);
    Vert *pV1 = SM.getVert(CF.iV1);
    Vert *pV2 = SM.getVert(CF.iV2);
    Face *pF01 = findCommonFace(pV0, pV1, true);
    Face *pF12 = findCommonFace(pV1, pV2, true);
    assert(pF01 && pF12);
    Cell *pC = findCommonCell(pF01, pF12);
    assert(pC);
    int iC = SM.getCellIndex(pC);
    aiTriStatus[iC] = 1; // Face is missing
  }

  logMessage(2, "Missing %zd segments\n", sCESubSegs.size());
  for (std::set<ConstrainEdge>::iterator iter = sCESubSegs.begin();
      iter != sCESubSegs.end(); ++iter) {
    // Missing input edges.
    ConstrainEdge CE = *iter;
    Vert *pV0 = SM.getVert(CE.iV0);
    Vert *pV1 = SM.getVert(CE.iV1);
    Face *pF01 = findCommonFace(pV0, pV1, true);
    assert(pF01);
    int iF = SM.getFaceIndex(pF01);
    aiEdgeStatus[iF] = 2;
    int iCL = SM.getCellIndex(pF01->getLeftCell());
    int iCR = SM.getCellIndex(pF01->getRightCell());
    aiTriStatus[iCL] += 4;
    aiTriStatus[iCR] += 4;
  }

  logMessage(2, "Missing %zd in-facet edges\n", sCEInFacets.size());
  for (std::set<ConstrainEdge>::iterator iter = sCEInFacets.begin();
      iter != sCEInFacets.end(); ++iter) {
    // Missing edges created in input n-gons, n>3.
    ConstrainEdge CE = *iter;
    Vert *pV0 = SM.getVert(CE.iV0);
    Vert *pV1 = SM.getVert(CE.iV1);
    Face *pF01 = findCommonFace(pV0, pV1, true);
    assert(pF01);
    int iF = SM.getFaceIndex(pF01);
    aiEdgeStatus[iF] = 1;
    int iCL = SM.getCellIndex(pF01->getLeftCell());
    int iCR = SM.getCellIndex(pF01->getRightCell());
    aiTriStatus[iCL] += 1;
    aiTriStatus[iCR] += 1;
  }
}

static void
vWriteFile_AnnotatedSurface(SurfMesh& OutMesh, const char strBaseFileName[],
			    const char strExtraFileSuffix[],
			    std::set<ConstrainEdge> sCESubSegs,
			    std::set<ConstrainEdge> sCEInFacets,
			    std::set<ConstrainFace> sCFSubFacets, VolMesh *pVM =
			    NULL)
{
  assert(OutMesh.isSimplicial());
  FILE *pFOutFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  sprintf(strFileName, "%s.vtk%s", strBaseFileName, strExtraFileSuffix);

  /*newfile vtk*/
  printf("Opening output file %s\n", strFileName);
  if (pFOutFile)
    fclose(pFOutFile);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "surface mesh output");

  /*"# vtk DataFile Version" discard_float*/
  fprintf(pFOutFile, "# vtk DataFile Version 1.1.1\n");
  /*discard_line*/
  fprintf(pFOutFile, "GRUMMP-generated surface mesh file\n");
  /*ASCII*/
  fprintf(pFOutFile, "ASCII\n");
  /*DATASET UNSTRUCTURED_GRID*/
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  /*POINTS nverts float*/
  fprintf(pFOutFile, "POINTS %u float\n", OutMesh.iNumVerts());

  for (GR_index_t iV = 0; iV < OutMesh.iNumVerts(); iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.14g %.14g %.14g\n", OutMesh.getVert(iV)->x(),
	    OutMesh.getVert(iV)->y(), OutMesh.getVert(iV)->z());

  /*"CELLS" ncells discard_integer*/
  fprintf(pFOutFile, "CELLS %u %u\n", OutMesh.getNumCells(),
	  OutMesh.getNumCells() * 4);

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    Cell* pC = OutMesh.getCell(iC);
    fprintf(pFOutFile, "3 ");
    for (int iV = 0; iV < pC->getNumVerts(); iV++) {
      fprintf(pFOutFile, "%u ", OutMesh.iVertIndex(pC->getVert(iV)));
    }
    /*cells: 3 verts*/
    fprintf(pFOutFile, "\n");
  }

  /*CELL_TYPES ncells*/
  fprintf(pFOutFile, "CELL_TYPES %u\n", OutMesh.getNumCells());

  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    /*cells: 5*/
    fprintf(pFOutFile, "5\n");
  }

  int *aiTriStatus = new int[OutMesh.getNumCells()];
  int *aiEdgeStatus = new int[OutMesh.getNumFaces()];

  vComputeEdgeAndTriStatus(OutMesh, sCESubSegs, sCEInFacets, sCFSubFacets,
			   aiTriStatus, aiEdgeStatus);

  // More diagnostics: this time, identifying blisters (small parts of
  // the surface mesh that have a different triangulation, but the right
  // number of triangles covering them).  Blisters should, I think,
  // prove easy to remove with a surface edge swap after ditching the
  // exterior.
  if (pVM) {
    for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); iF++) {
      if (aiEdgeStatus[iF] != 0) {
	Face *pF = OutMesh.getFace(iF);
	Cell *pCL = pF->getLeftCell();
	Cell *pCR = pF->getRightCell();
	int iCL = OutMesh.getCellIndex(pCL);
	int iCR = OutMesh.getCellIndex(pCR);
	int iStatL = aiTriStatus[iCL];
	int iStatR = aiTriStatus[iCR];
	if ((iStatL == 2 || iStatL == 5) && (iStatR == 2 || iStatR == 5)) {
	  // Each has only one missing edge: the easiest case.
	  // So only four verts:
	  Vert *pVA = pF->getVert(0);
	  Vert *pVB = pF->getVert(1);
	  Vert *pVL = pCL->getOppositeVert(pF);
	  Vert *pVR = pCR->getOppositeVert(pF);
	  // Those are the surface verts.  Now get the volume verts.
	  Vert *pVA_vol = pVM->getVert(OutMesh.iVertIndex(pVA));
	  Vert *pVB_vol = pVM->getVert(OutMesh.iVertIndex(pVB));
	  Vert *pVL_vol = pVM->getVert(OutMesh.iVertIndex(pVL));
	  Vert *pVR_vol = pVM->getVert(OutMesh.iVertIndex(pVR));
	  // The two faces you have to have here are known:
	  Face *pFCand1 = findCommonFace(pVA_vol, pVL_vol, pVR_vol,
	  NULL,
					 true);
	  Face *pFCand2 = findCommonFace(pVB_vol, pVL_vol, pVR_vol,
	  NULL,
					 true);
	  if (pFCand1->isValid() && pFCand2->isValid()) {
	    aiTriStatus[iCL] = aiTriStatus[iCR] = -2;

	    // How many faces are incident on the bad edge here?
	    int iNFaces_Left = pVL_vol->getNumFaces();
	    int iCommon = 0;
	    for (int ii = 0; ii < iNFaces_Left; ii++) {
	      Face *pFTmp = pVL_vol->getFace(ii);
	      if (pFTmp->hasVert(pVR_vol))
		iCommon++;
	    }
	    // 	    vMessage(2, "Two-tri blister; %d common faces\n", iCommon);
	  }
	  else if (pFCand1->isValid() || pFCand2->isValid()) {
	    aiTriStatus[iCL] = aiTriStatus[iCR] = -1;
	  }
	}
      }
    }
  }

  fprintf(pFOutFile, "CELL_DATA %u\n", OutMesh.getNumCells());
  fprintf(pFOutFile, "SCALARS tri_status int 1\n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); iC++) {
    fprintf(pFOutFile, "%d\n", aiTriStatus[iC]);
  }

  delete[] aiTriStatus;
  delete[] aiEdgeStatus;

  /* Make sure compilers won't complain that these variables */
  /* aren't used. */
  fclose(pFOutFile);
}

void
VolMesh::splitEdge(const Vert* const pV0, const Vert* const pV1, SurfMesh& SM)
{
  // Force a split of the given surface edge within both the surface and
  // volume meshes.
  int iV0 = getVertIndex(pV0);
  int iV1 = getVertIndex(pV1);
  Vert *pVS0 = SM.getVert(iV0);
  Vert *pVS1 = SM.getVert(iV1);
  assert(dDIST3D(pVS0->getCoords(), pV0->getCoords()) < 1.e-12);
  assert(dDIST3D(pVS1->getCoords(), pV1->getCoords()) < 1.e-12);
  Face *pFS = pFInvalidFace;
  {
    std::set<Vert*> spVNear;
    std::set<Face*> spFNear;
    std::set<Cell*> spCTmp;
    std::set<BFace*> spBFTmp;
    bool qTmp;
    findNeighborhoodInfo(pVS0, spCTmp, spVNear, &spBFTmp, &qTmp, &spFNear);
    if (spVNear.count(pVS1) == 0)
      return;
    std::set<Face*>::iterator iterF;
    for (iterF = spFNear.begin(); iterF != spFNear.end(); ++iterF) {
      pFS = *iterF;
      if (pFS->hasVert(pVS0) && pFS->hasVert(pVS1))
	break;
    }
    assert(pFS->hasVert(pVS0));
    assert(pFS->hasVert(pVS1));
  }

  // Edge swap failed; have to solve this by splitting.
  double adNewLoc[3];
#ifdef SMALL_ANGLE_RECOVERY
  bool qSmall0 = pVS0->isSmallAngleVert();
  bool qSmall1 = pVS1->isSmallAngleVert();
  if ((qSmall0 && qSmall1) || (!qSmall0 && !qSmall1)) {
#endif
    // Bisect
    adNewLoc[0] = 0.5 * (pV0->x() + pV1->x());
    adNewLoc[1] = 0.5 * (pV0->y() + pV1->y());
    adNewLoc[2] = 0.5 * (pV0->z() + pV1->z());
    logMessage(3, "Splitting edge at midpoint: %12.6f %12.6f %12.6f\n",
	       adNewLoc[0], adNewLoc[1], adNewLoc[2]);
#ifdef SMALL_ANGLE_RECOVERY
  }
  else if (qSmall0) {
    // Split near pVS0.
    double dSplitDist = pVS0->getLengthScale();
    double adVec[] = adDIFF3D(pV1->getCoords(), pV0->getCoords());
    double dEdgeLen = dMAG3D(adVec);
    while (dSplitDist > dEdgeLen / 1.5)
      dSplitDist /= 2;
    assert(dSplitDist < 2. / 3. * dEdgeLen);

    NORMALIZE3D(adVec);

    adNewLoc[0] = pV0->x() + adVec[0] * dSplitDist;
    adNewLoc[1] = pV0->y() + adVec[1] * dSplitDist;
    adNewLoc[2] = pV0->z() + adVec[2] * dSplitDist;
    logMessage(3, "Splitting edge near vert %4d: (%12.6f, %12.6f, %12.6f) %f\n",
	       iV0, adNewLoc[0], adNewLoc[1], adNewLoc[2], dSplitDist);
  }
  else {
    assert(qSmall1);
    // Split near pVS1.
    double dSplitDist = pVS1->getLengthScale();
    double adVec[] = adDIFF3D(pV0->getCoords(), pV1->getCoords());
    double dEdgeLen = dMAG3D(adVec);
    while (dSplitDist > dEdgeLen / 1.5)
      dSplitDist /= 2;
    assert(dSplitDist < 2. / 3. * dEdgeLen);

    NORMALIZE3D(adVec);

    adNewLoc[0] = pV1->x() + adVec[0] * dSplitDist;
    adNewLoc[1] = pV1->y() + adVec[1] * dSplitDist;
    adNewLoc[2] = pV1->z() + adVec[2] * dSplitDist;
    logMessage(3, "Splitting edge near vert %4d: (%12.6f, %12.6f, %12.6f) %f\n",
	       iV0, adNewLoc[0], adNewLoc[1], adNewLoc[2], dSplitDist);
  }
#endif

  Vert *pVNewSurf = SM.createVert(adNewLoc[0], adNewLoc[1], adNewLoc[2]);
  bool qInsertOnCurve;
  // FIX ME  What about non-manifold edges?  Subsequent recovery of the
  // original surface is significantly harder then...
  SM.insertOnFace(pVNewSurf, pFS, pFS->getLeftCell(), true, &qInsertOnCurve);

  Vert *pVNew = pSI3D->insertPoint(adNewLoc, beginPipe(pV0, pV1));
  // This vertex was added during bdry recovery, along a curve.  Mark it
  // as a bdry curve point, and tag for possible deletion, in the case
  // that we're interested in removing as many constraining verts as possible.
  pVNew->setType(qInsertOnCurve ? Vert::eBdryCurve : Vert::eBdry);
  pVNew->markToDelete();
  stackSplitEdges.push(
      SplitEdge(const_cast<Vert*>(pV0), const_cast<Vert*>(pV1), pVNew,
		pVNewSurf, qInsertOnCurve));
}

int
VolMesh::recoverEdge(const Vert* const pV0, const Vert* const pV1,
		     const bool /*qInsertionOK*/)
{
  assert(isSimplicial());
  // commented out for now
  //  if (qInsertionOK) {
  // TODO: May want to look at this again..
  //    // FIX ME: Shouldn't bother with this unless the segment is
  //    // encroached.  Why?  Because if no segments are encroached, then
  //    // we can always recover the bdry segs.  (If segments are present
  //    // but encroached, that can always be fixed later.)
  //
  //    // Curved bdry: replace with code that splits the underlying patch
  //    // boundary properly, instead of simply taking the midpoint of a
  //    // segment.
  //    double adNewCoords[] = {0.5 * (pV0->x() + pV1->x()),
  //			    0.5 * (pV0->y() + pV1->y()),
  //			    0.5 * (pV0->z() + pV1->z())};
  //    Vert *pVBdry0 = m_Bdry3D->pVVert(getVertIndex(pV0));
  //    Vert *pVBdry1 = m_Bdry3D->pVVert(getVertIndex(pV1));
  //    m_Bdry3D->pVAddVertex(pVBdry0, pVBdry1, adNewCoords);
  //    {
  //      int iCell = 0;
  //      Cell *pC = getCell(iCell);
  //      while (!pC->doFullCheck()) {
  //	pC = getCell(++iCell);
  //      }
  //      // Any live tet will do, since the domain is concave.
  //      Vert* newVert = pSI3D->insertPoint(adNewCoords, pC);
  //      // Mark this vert correctly
  //      newVert->setType(Vert::eBdryCurve);
  //      assert(newVert->isValid());
  //    }
  //    // That particular edge certainly doesn't need to be recovered; a
  //    // return value of 1 makes sure we never try to again.
  //    return 1;
  //  }

  int iChange, iNSkip, iNPasses = 0, iRetVal = 0;
  Face* apFSkipped[20];
  do {
    iChange = iNSkip = 0;
    bool qHitEdge = false;
    iNPasses++;
    logMessage(4, "Pass %d to remove this edge\n", iNPasses);

    // Walk from pV0 towards pV1 until either (a) you hit an edge
    // or (b) you get to pV1.  Keep track of all faces you had to skip
    // over.
    int iRes = clearPipe(pV0, pV1, qHitEdge, apFSkipped, iNSkip);
    if (iRes == 1) // Pipe is unclogged
      return 1;
    if (iRes == 0) { // Changes made
      iChange = 1;
      iRetVal = -1;
    }
    if (qHitEdge) { // Try from the other end
      qHitEdge = false;
      iNSkip = 0;
      iRes = clearPipe(pV1, pV0, qHitEdge, apFSkipped, iNSkip);
      if (iRes == 1) // Pipe is unclogged; not bloody likely
	return 1;
      if (iRes == 0) { // Something has changed
	iChange = 1;
	iRetVal = -1;
      }
    }

    // If an edge was hit, try to remove it
    if (qHitEdge) {
      logMessage(3, "Unresolved edge hit!  Edge not recovered.\n");
    }
  }
  while (iChange && iNPasses < 5);

  return iRetVal;
}

int
VolMesh::recoverEdge(const Vert* const pV0In, const Vert* const pV1In,
		     SurfMesh& SM, const bool qInsertionOK)
{
  int iSwaps = 0;
  assert(isSimplicial());

  int iNSkip = 0, iRetVal = 0;
  Face* apFSkipped[20];
  bool qHitEdge = false;

  logMessage(3, "In RecoverEdge %p %p\n", pV0In, pV1In);
  // A modified approach to edge recovery:  volume mesh edges
  // Walk from pV0 towards pV1 until either (a) you hit an edge
  // or (b) you get to pV1.  Keep track of all faces you had to skip
  // over.
  int iRes, iLaps = 0;
  const Vert *pV0 = pV0In;
  const Vert *pV1 = pV1In;
  do {
    iRes = clearPipe(pV0, pV1, qHitEdge, apFSkipped, iNSkip);
    if (iRes == 1) // Pipe is unclogged
      return 1;
    const Vert *pVTmp = pV0;
    pV0 = pV1;
    pV1 = pVTmp;
    iLaps++;
  }
  while (iRes == 0 && iLaps <= 6); // Changes made; reverse and try again

  // If an edge was hit, this (almost certainly) means that an edge in
  // the surface mesh and an edge in the volume mesh both lie in the
  // same coplanar facet and that the four points in question are
  // co-circular.  In this situation, the easiest fix is to swap away
  // the offending edge in the -surface- mesh, leaving the volume mesh
  // edge alone.  (While in sufficiently complex degenerate cases, a
  // single swap isn't enough to fix this, multiple swaps (perhaps for
  // different recover-edge targets) will in fact suffice.)

  // After making this forced swap, return with value -1, to indicate
  // that something useful was done, even though an edge is not
  // certain to have been recovered.  Actually, you can do even
  // better: when the edge -has- been recovered, you know that, so you
  // can return 1.

  if (qHitEdge) {
    logMessage(3, "  Hit an edge! Trying to remove it...\n");
    int iV0 = getVertIndex(pV0);
    int iV1 = getVertIndex(pV1);
    Vert *pVS0 = SM.getVert(iV0);
    Vert *pVS1 = SM.getVert(iV1);
    assert(dDIST3D(pVS0->getCoords(), pV0->getCoords()) < 1.e-12);
    assert(dDIST3D(pVS1->getCoords(), pV1->getCoords()) < 1.e-12);
    Face *pFS = pFInvalidFace;
    {
      std::set<Vert*> spVNear;
      std::set<Face*> spFNear;
      std::set<Cell*> spCTmp;
      std::set<BFace*> spBFTmp;
      bool qTmp;
      findNeighborhoodInfo(pVS0, spCTmp, spVNear, &spBFTmp, &qTmp, &spFNear);
      assert(spVNear.count(pVS1) == 1);
      std::set<Face*>::iterator iterF;
      for (iterF = spFNear.begin(); iterF != spFNear.end(); ++iterF) {
	pFS = *iterF;
	if (pFS->hasVert(pVS0) && pFS->hasVert(pVS1))
	  break;
      }
      assert(pFS->hasVert(pVS0));
      assert(pFS->hasVert(pVS1));
    }

    // Identify the cells and patches on either side (patches should
    // be the same).
    Cell *pCLeft = pFS->getLeftCell();
    Cell *pCRight = pFS->getRightCell();
    BdryPatch3D *pBPLeft = SM.getBdryPatch(SM.getCellIndex(pCLeft));
    BdryPatch3D *pBPRight = SM.getBdryPatch(SM.getCellIndex(pCRight));
    if (pBPLeft != pBPRight) {
      // A surface edge not internal to a facet.

      // Try to swap away this edge in the volume mesh.  First
      // identify the edge.
      Vert *pVA = apFSkipped[0]->getVert(0);
      Vert *pVB = apFSkipped[0]->getVert(1);
      Vert *pVC = apFSkipped[0]->getVert(2);

      Vert *pVNorth, *pVSouth, *pVOther;
      Vert *pVV0 = getVert(iV0);
      Vert *pVV1 = getVert(iV1);

      if (checkOrient3D(pVA, pVB, pVS0, pVS1) == 0) {
	pVNorth = pVA;
	pVSouth = pVB;
	pVOther = pVC;
      }
      else if (checkOrient3D(pVB, pVC, pVS0, pVS1) == 0) {
	pVNorth = pVB;
	pVSouth = pVC;
	pVOther = pVA;
      }
      else {
	assert(checkOrient3D(pVC, pVA, pVS0, pVS1) == 0);
	pVNorth = pVC;
	pVSouth = pVA;
	pVOther = pVB;
      }
      iSwaps = edgeSwap3D(apFSkipped[0], pVNorth, pVSouth, pVOther, true, pVV0,
			  pVV1);
      if (iSwaps > 0) {
	// Success!  Well, at least we got rid of the offending edge.
	// So return -1; if the correct edge showed up in the result,
	// we'll figure it out.
	logMessage(3, "   Swapped away the edge in the volume mesh.\n");
	return -1;
      }
      else if (qInsertionOK) {
	// Edge swap failed; have to solve this by splitting.
	double adNewLoc[] =
	  { 0.5 * (pVS0->x() + pVS1->x()), 0.5 * (pVS0->y() + pVS1->y()), 0.5
	      * (pVS0->z() + pVS1->z()) };

	Vert *pVNewSurf = SM.createVert(adNewLoc[0], adNewLoc[1], adNewLoc[2]);
	iSwaps = SM.insertOnFace(pVNewSurf, pFS, pCLeft, true);

	Vert *pVNew = pSI3D->insertPoint(adNewLoc,
					 apFSkipped[0]->getLeftCell());
	pVNew->setType(Vert::eBdryCurve);
	// We may eventually work at trying to remove verts that were
	// added to constrain the boundary.  If/when that is done, we'll
	// be glad that this vert is marked for deletion.
	pVNew->markToDelete();

	logMessage(4, "   Split the edge in both surface and volume meshes.\n");
	// This edge definitely won't need to be recovered now...
	return 1;
      }
      else {
	logMessage(4, "   Not splitting; did nothing.\n");
      }
    } // Different patches on the two sides of the surface face.
    else {
      // Internal facet edge (most common case; four or more
      // cocircular points in a facet boundary)

      // Are the four verts forming those two cells (a) coplanar and (b)
      // cocircular?
      assert(pCLeft->getType() == Cell::eTriCell);
      assert(pCRight->getType() == Cell::eTriCell);
      Vert *pVLeft = dynamic_cast<TriCell*>(pCLeft)->getOppositeVert(pFS);
      Vert *pVRight = dynamic_cast<TriCell*>(pCRight)->getOppositeVert(pFS);

      if (checkOrient3D(pV0, pV1, pVLeft, pVRight) != 0) {
	logMessage(3, "    Peculiar.  I was expecting coplanar triangles.\n");
	logMessage(3, "    Nothing I can do here.\n");
	return 0;
      }

      Vert VDummy;
      double adNorm[3];
      calcNormal3D(pV0->getCoords(), pV1->getCoords(), pVLeft->getCoords(),
		   adNorm);
      NORMALIZE3D(adNorm);
      double adDummyLoc[] =
	{ pV0->x() + adNorm[0], pV0->y() + adNorm[1], pV0->z() + adNorm[2] };
      VDummy.setCoords(3, adDummyLoc);

      int iInside = isInsphere(pV0, pV1, pVLeft, &VDummy, pVRight);
      if (iInside != 0) {
	logMessage(3, "    Peculiar.  I was expecting cocircular verts.\n");
	logMessage(3, "    Nothing I can do here.\n");
	return 0;
      }

      // Find the offending edge in the VolMesh.  (apFSkipped[0] is a
      // Face that contains that edge, but which edge is it?)
      Vert *pVA = apFSkipped[0]->getVert(0);
      Vert *pVB = apFSkipped[0]->getVert(1);
      Vert *pVC = apFSkipped[0]->getVert(2);
      int iVA = getVertIndex(pVA);
      int iVB = getVertIndex(pVB);
      int iVC = getVertIndex(pVC);

      int iVLeft = SM.iVertIndex(pVLeft);
      int iVRight = SM.iVertIndex(pVRight);

      // First, the easy way: is one of the edges of ABC the same as
      // edge LR in the surface mesh?

      bool qChangeSurf = false, qExactRecovery;
      if (((iVLeft == iVA && iVRight == iVB)
	  || (iVLeft == iVB && iVRight == iVA))
	  || ((iVLeft == iVB && iVRight == iVC)
	      || (iVLeft == iVC && iVRight == iVB))
	  || ((iVLeft == iVC && iVRight == iVA)
	      || (iVLeft == iVA && iVRight == iVC))) {
	// This swap will produce in the surface mesh the edge that
	// exists in the volume mesh.
	qExactRecovery = true;
	qChangeSurf = true;
      }
      else {
	// This swap will eliminate an edge of the surface mesh that
	// intersects the desired edge of the volume mesh, but won't
	// recover the volume edge.
	qExactRecovery = false;

	// Check each possible edge in turn
	bool qABGuilty = ((checkOrient3D(pV0, pV1, pVA, pVB) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVA) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVB) == 0));
	bool qBCGuilty = ((checkOrient3D(pV0, pV1, pVB, pVC) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVB) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVC) == 0));
	bool qCAGuilty = ((checkOrient3D(pV0, pV1, pVC, pVA) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVC) == 0)
	    && (isInsphere(pV0, pV1, pVLeft, &VDummy, pVA) == 0));
	if (qABGuilty || qBCGuilty || qCAGuilty) {
	  qChangeSurf = true;
	}
      }
      if (qChangeSurf) {
	// Reconfigure the SurfMesh to remove the edge and return.
	SM.disallowSwapRecursion();
	logMessage(3, "  Trying a surface edge swap.\n");
	iSwaps = SM.reconfigure(pFS);
	assert(iSwaps == 1);
	if (qExactRecovery) {
	  logMessage(
	      3,
	      "   Made a co-circular surface swap; should have gotten the edge I wanted.\n");
	  return 1;
	}
	else {
	  logMessage(
	      3,
	      "   Made a co-circular surface swap; didn't get the edge I wanted.\n");
	  return -1;
	}
      }
      else {
	logMessage(4, "   Failed to do anything to the edge.\n");
      }
    } // Edge internal to facet
  } // Hit an edge square on

  return iRetVal;
}

static void
vRecoverEdges(VolMesh * const pVM, SurfMesh * const pSM,
	      std::set<ConstrainEdge>& sCE, int& iRecovered, int& iOtherChanges)
{
  iRecovered = iOtherChanges = 0;
  std::set<ConstrainEdge>::iterator iter, iterEnd = sCE.end();
  for (iter = sCE.begin(); iter != iterEnd;) {
    Vert *pV0 = pVM->getVert(iter->iV0);
    Vert *pV1 = pVM->getVert(iter->iV1);
    int iGotIt = pVM->recoverEdge(pV0, pV1, *pSM, false);
    std::set<ConstrainEdge>::iterator iterTmp = iter;
    ++iter;
    switch (iGotIt)
      {
      case 1: // Edge recovered
	sCE.erase(iterTmp);
	iRecovered++;
	break;
      case -1:
	iOtherChanges++;
	break;
      }
  } // Done checking non-constrained segments
}

static bool
qSplitEdges(VolMesh * const pVM, SurfMesh * const pSM,
	    std::set<ConstrainEdge>& sCESubSegs,
	    std::set<ConstrainEdge>& sCEInFacets, std::set<ConstrainFace>& sCF,
	    int& iSplittings)
{
#ifndef NDEBUG
  //  vWriteFile_AnnotatedSurface(*pSM, "presplit", "",
  //			      sCESubSegs, sCEInFacets, sCF,
  //			      pVM);
#endif
  //   // Change the swapping type.  Insertion into the volume mesh sometimes
  //   // fails with flat tets otherwise.  Since faces that have been
  //   // recovered are locked, this shouldn't hurt things too much.
  //   eSwap eSwapTmp = pVM->eSwapType();
  //   pVM->vSetSwapType(eMinSine);

  logMessage(2, "Splitting all unconstrained edges (%dth time).\n",
	     iSplittings);
  logMessage(2, "  Constraint stats before:\n");
  isBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCF);
  // Always recover subsegments first, then work on edges
  // within subfacets.
  std::set<ConstrainEdge> sCE(sCESubSegs);
  //   if (sCE.empty()) {
  sCE.insert(sCEInFacets.begin(), sCEInFacets.end());
  //   }
  // Since allowing insertion wasn't enough, force a split of
  // all edges not yet constrained.  Enable swap recursion in
  // the volume mesh while splitting.
  pVM->allowSwapRecursion();
  int iSmallVerts = 0, iSmallEdges = 0, iTotalEdges = sCE.size();
  std::set<ConstrainEdge>::iterator iter;
  int iE = 0;
  for (iter = sCE.begin(); iter != sCE.end(); ++iter, iE++) {
    ConstrainEdge CE = *iter;
    Vert *pV0 = pVM->getVert(CE.iV0);
    Vert *pV1 = pVM->getVert(CE.iV1);
    logMessage(3, "   Edge %d.  Vert %d: %s.  Vert %d: %s.\n", iE, CE.iV0,
	       (pV0->isSmallAngleVert() ? "small" : " okay"), CE.iV1,
	       (pV1->isSmallAngleVert() ? "small" : " okay"));
    iSmallVerts += (pV0->isSmallAngleVert() ? 1 : 0)
	+ (pV1->isSmallAngleVert() ? 1 : 0);
    iSmallEdges += (pV0->isSmallAngleVert() || pV1->isSmallAngleVert() ? 1 : 0);
    pVM->splitEdge(pV0, pV1, *pSM);
  }
  logMessage(2, "  Out of %d edges, %d had a small angle vert (%d verts)\n",
	     iTotalEdges, iSmallEdges, iSmallVerts);
  pVM->disallowSwapRecursion();

  logMessage(2, "  Constraint stats after:\n");
  // Reset the list of stuff to recover
  bool retVal = isBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCF);
#ifndef NDEBUG
  //  vWriteFile_AnnotatedSurface(*pSM, "postsplit","",
  //			      sCESubSegs, sCEInFacets, sCF,
  //			      pVM);
#endif
  //   exit(110);

  //   // Reset the swap type.
  //   pVM->vSetSwapType(eSwapTmp);

  return retVal;
}

static void
vRecoverFaces(VolMesh * const pVM, SurfMesh * const pSM,
	      std::set<ConstrainFace> sCF, int& iFacesRecovered,
	      int& iFaceRec3To2, int& iFaceRecEdge, const bool qSplitBadFaces)
{
  iFacesRecovered = 0;
  std::set<ConstrainFace>::iterator iter, iterEnd = sCF.end();
  for (iter = sCF.begin(); iter != iterEnd;) {
    ConstrainFace CF = *iter;
    int aiVerts[] =
      { CF.iV0, CF.iV1, CF.iV2 };
    int ii = pVM->recoverFace(aiVerts);
    std::set<ConstrainFace>::iterator iterTmp = iter;
    ++iter;
    if (ii) {
      iFacesRecovered++;
      sCF.erase(iterTmp);
      // This is an ugly hack, but it works.
      iFaceRec3To2 += ii % 1000;
      iFaceRecEdge += ii / 1000;
    }
    else if (qSplitBadFaces) {
      logMessage(3, "Splitting all edges of bad face.\n");
      // Didn't get it, so split its edges.  Yes, it's overkill to split
      // all three probably...
      pVM->splitEdge(pVM->getVert(CF.iV0), pVM->getVert(CF.iV1), *pSM);
      pVM->splitEdge(pVM->getVert(CF.iV1), pVM->getVert(CF.iV2), *pSM);
      pVM->splitEdge(pVM->getVert(CF.iV2), pVM->getVert(CF.iV0), *pSM);
    }
  }
}

// static void vRecoverAllEdges(VolMesh * const pVM,
// 			     SurfMesh * const pSM,
// 			     std::set<ConstrainEdge>& sCESubSegs,
// 			     std::set<ConstrainEdge>& sCEInFacets,
// 			     std::set<ConstrainFace>& sCF)
// {
//   std::set<ConstrainEdge> sCE(sCESubSegs);
//   if (sCE.empty()) sCE.insert(sCEInFacets.begin(), sCEInFacets.end());

//   int iSplittings = 0;
//   if (!sCE.empty()) {
//     //@@ Recover edges
//     // Verts added to the volume mesh to constrain it; must be added to
//     // surface mesh on return for constraint checking to work properly.

//     int iChange, iRecovered;
//     bool qDoInsertion = false;
//     int iTotalPasses;
//     do { // Begin recovering edges
//       iTotalPasses = 0;
//       int iNFruitlessPasses = 0;
//       do { // Begin a pass to check all unconstrained edges
// 	vMessage(3, "Recovery pass %d (%d fruitless); insertion %s\n",
// 		 iTotalPasses, iNFruitlessPasses,
// 		 (qDoInsertion ? "allowed" : "not allowed"));
// 	iTotalPasses ++;
// 	iNFruitlessPasses ++;
// 	int iVertsBefore = pVM->iNumVerts();
// 	vRecoverEdges(pVM, pSM, sCE, iRecovered, iChange);
// 	std::set<ConstrainEdge>::iterator iter;
// 	if (iRecovered) iNFruitlessPasses = 0;
// 	int iVertsAfter = pVM->iNumVerts();
// 	vMessage(2, "Status:%2d passes (%sinserting; last %d fruitless). Changes: %3d Left: %3d\n",
// 		 iTotalPasses, qDoInsertion ? "" : "not ",
// 		 iNFruitlessPasses, iRecovered, static_cast<int>(sCE.size()));
// 	if (sCE.empty()) {
// 	  qIsBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCF);
// 	  // Always recover subsegments first, then work on edges
// 	  // within subfacets.
// 	  assert(sCE.empty());
// 	  if (!sCESubSegs.empty()) {
// 	    sCE.insert(sCESubSegs.begin(), sCESubSegs.end());
// 	  }
// 	  else {
// 	    sCE.insert(sCEInFacets.begin(), sCEInFacets.end());
// 	  }
// 	  // If swapping fixed all the edges, don't turn off insertion,
// 	  // or you could loop forever trying to do surface recovery.
// 	  if (qDoInsertion && (iVertsAfter != iVertsBefore)) {
// 	    qDoInsertion = false;
// 	    assert(iNFruitlessPasses == 0);
// 	    iTotalPasses = 0;
// 	  }
// 	}
// 	for (iter = sCE.begin(); iter != sCE.end(); ++iter) {
// 	  ConstrainEdge CE = *iter;
// 	  Vert *pV0 = pVM->pVVert(CE.iV0);
// 	  Vert *pV1 = pVM->pVVert(CE.iV1);
// 	  vMessage(4, "Bad edge:\n");
// 	  vMessage(4, "%d %12.6f %12.6f %12.6f\n",
// 		   CE.iV0, pV0->dX(), pV0->dY(),
// 		   pV0->dZ());
// 	  vMessage(4, "%d %12.6f %12.6f %12.6f\n\n",
// 		   CE.iV1, pV1->dX(), pV1->dY(),
// 		   pV1->dZ());
// 	}
//       } while (iNFruitlessPasses < 2 && iChange && (!sCE.empty())
// 	       && iTotalPasses < 10); // Done with edge recovery, at
// 				      // least for now

//       // We might even be totally done...
//       if (sCE.empty()) return;
//       // Don't insert verts onto boundaries that are immutable.
//       if (pVM->qBdryChangesAllowed()) {
// 	if (qDoInsertion) {
// 	  iSplittings++;
// 	  if (iSplittings <= iMaxSplittings) {
// 	    qSplitEdges(pVM, pSM, sCESubSegs, sCEInFacets, sCF, iSplittings);
// 	    qDoInsertion = false;
// 	    // Always recover subsegments first, then work on edges
// 	    // within subfacets.
// 	    sCE.clear();
// 	    if (!sCESubSegs.empty()) {
// 	      sCE.insert(sCESubSegs.begin(), sCESubSegs.end());
// 	    }
// 	    else {
// 	      sCE.insert(sCEInFacets.begin(), sCEInFacets.end());
// 	    }
// 	  }
// 	} // Done splitting all missing edges
// 	else {
// 	  // Some diagnostic output, to see what sorts of patterns there
// 	  // are in unrecovered edges.
//  	  qIsBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCF);
// 	  // Now find all the SurfMesh Face's that are missing.
// 	  std::set<Face*> spFMissingFaces;
// 	  {
// 	    std::set<ConstrainEdge>::iterator iter, end = sCESubSegs.end();
// 	    for (iter = sCESubSegs.begin(); iter != end; ++iter) {
// 	      ConstrainEdge CE = *iter;
// 	      Vert *pV0 = pSM->pVVert(CE.iV0);
// 	      Vert *pV1 = pSM->pVVert(CE.iV1);
// 	      Face* pF = findCommonFace(pV0, pV1, true);
// 	      assert(pF->qValid());
// 	      spFMissingFaces.insert(pF);
// 	    }
// 	    // Now again for edges in facets.
// 	    end = sCEInFacets.end();
// 	    for (iter = sCEInFacets.begin(); iter != end; ++iter) {
// 	      ConstrainEdge CE = *iter;
// 	      Vert *pV0 = pSM->pVVert(CE.iV0);
// 	      Vert *pV1 = pSM->pVVert(CE.iV1);
// 	      Face* pF = findCommonFace(pV0, pV1, true);
// 	      assert(pF->qValid());
// 	      spFMissingFaces.insert(pF);
// 	    }
// 	  }
// 	  std::set<Cell*> spCMissingCells;
// 	  {
// 	    std::set<ConstrainFace>::iterator iter, end = sCF.end();
// 	    for (iter = sCF.begin(); iter != end; ++iter) {
// 	      ConstrainFace CF = *iter;
// 	      Vert *pV0 = pSM->pVVert(CF.iV0);
// 	      Vert *pV1 = pSM->pVVert(CF.iV1);
// 	      Vert *pV2 = pSM->pVVert(CF.iV2);
// 	      Face *pF01 = findCommonFace(pV0, pV1, true);
// 	      Face *pF12 = findCommonFace(pV1, pV2, true);
// 	      Face *pF20 = findCommonFace(pV2, pV0, true);
// 	      assert(pF01->qValid() && pF12->qValid() && pF20->qValid());
// 	      Cell *pC = pCCommonCell(pF01, pF12);
// 	      assert(pC->qHasFace(pF20));
// 	      spCMissingCells.insert(pC);
// 	    }
// 	  }
// 	  int *aiMissingData = new int[pSM->iNumCells()];
// 	  for (int ii = pSM->iNumCells() - 1; ii >= 0; ii--)
// 	    aiMissingData[ii] = 0;
// 	  {
// 	    std::set<Face*>::iterator iter, end = spFMissingFaces.end();
// 	    for (iter = spFMissingFaces.begin(); iter != end; ++iter) {
// 	      Face *pF = *iter;
// 	      int iCL = pSM->iCellIndex(pF->pCCellLeft());
// 	      int iCR = pSM->iCellIndex(pF->pCCellRight());
// 	      aiMissingData[iCL]++;
// 	      aiMissingData[iCR]++;
// 	    }
// 	  }
// 	  {
// 	    std::set<Cell*>::iterator iter, end = spCMissingCells.end();
// 	    for (iter = spCMissingCells.begin(); iter != end; ++iter) {
// 	      Cell *pC = *iter;
// 	      int iC = pSM->iCellIndex(pC);
// 	      if (aiMissingData[iC] == 0) aiMissingData[iC] = 4;
// 	    }
// 	  }
// 	  int aiCounts[] = {0,0,0,0,0};
// 	  for (GR_index_t i = 0; i < pSM->iNumCells(); i++) {
// 	    aiCounts[aiMissingData[i]] ++;
// // 	    if (aiMissingData[i] != 0) {
// // 	      vMessage(2, "Cell %5d:  %d\n", i, aiMissingData[i]);
// // 	    }
// 	  }
// 	  vMessage(2, "Tris fully present: %d\n", aiCounts[0]);
// 	  vMessage(2, "Tris w/ one missing edge: %d\n", aiCounts[1]);
// 	  vMessage(2, "Tris w/ two missing edge: %d\n", aiCounts[2]);
// 	  vMessage(2, "Tris w/ three missing edge: %d\n", aiCounts[3]);
// 	  vMessage(2, "Tris w/ all edges, but missing: %d\n", aiCounts[4]);
// 	  // Try insertion if you haven't recently.
// 	  qDoInsertion = true;
// 	}
//       }
//     } while (iSplittings <= iMaxSplittings && (!sCE.empty()));
//     // Done recovering edges
//   } // Edges only need to be recovered if they are any missing in the first
//     // place.
// }

void
recoverSurface(VolMesh * const pVM, SurfMesh * const pSM)
{
  //@@ Take a census to see which edges and faces need work to recover
  std::set<ConstrainEdge> sCESubSegs;
  std::set<ConstrainEdge> sCEInFacets;
  std::set<ConstrainFace> sCFSubFacets;

  // Find out how many bdry (sub)segments are not currently in the mesh.
  bool qOK = isBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCFSubFacets);
  //#ifndef NDEBUG
  //  vWriteFile_AnnotatedSurface(*pSM, "init", "",
  //			      sCESubSegs, sCEInFacets, sCFSubFacets,
  //			      pVM);
  //#endif

  int iOuterPasses = 0;
  int iFaceRec3To2 = 0, iFaceRecEdge = 0;
  int iSplittings = 0;
  bool qRecoveredAllEdges = (sCESubSegs.empty() && sCEInFacets.empty());
  bool qTriedFacesBefore = false;
  while (!qOK && iOuterPasses < 10) {
    iOuterPasses++;
    logMessage(2, "Recovering surface edges (outer pass %d)...\n",
	       iOuterPasses);
    int iInnerPasses = 0, iTotalRec;
    do {
      int iRecovered = 0, iOtherChanges;

      int iInitSubSegs = sCESubSegs.size();
      int iInitInFacets = sCEInFacets.size();
      int iInitSubFacets = sCFSubFacets.size();

      if (!sCESubSegs.empty()) {
	vRecoverEdges(pVM, pSM, sCESubSegs, iRecovered, iOtherChanges);
	logMessage(2, "Recovered %d missing subsegments; "
		   "made %d other changes.\n",
		   iRecovered, iOtherChanges);
	if (iRecovered)
	  qTriedFacesBefore = false;
      }

      //       if (iRecovered == 0 && !sCEInFacets.empty()) {
      vRecoverEdges(pVM, pSM, sCEInFacets, iRecovered, iOtherChanges);
      logMessage(2, "Recovered %d missing in-facet edges; "
		 "made %d other changes.\n",
		 iRecovered, iOtherChanges);
      if (iRecovered)
	qTriedFacesBefore = false;
      //       }

      //       if (iRecovered == 0) {
      vRecoverFaces(pVM, pSM, sCFSubFacets, iRecovered, iFaceRec3To2,
		    iFaceRecEdge, qRecoveredAllEdges && qTriedFacesBefore);
      qTriedFacesBefore = true;
      logMessage(2, "Recovered %d missing faces; "
		 "%d 3-to-2 swaps, %d full edge swaps total.\n",
		 iRecovered, iFaceRec3To2, iFaceRecEdge);
      if (iRecovered)
	qTriedFacesBefore = false;
      //       }

      qOK = isBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets, sCFSubFacets);
      qRecoveredAllEdges = (sCESubSegs.empty() && sCEInFacets.empty());

      iTotalRec = ((iInitSubSegs - sCESubSegs.size())
	  + (iInitInFacets - sCEInFacets.size())
	  + (iInitSubFacets - sCFSubFacets.size()));
      logMessage(2, "Total recoveries this pass: %d\n", iTotalRec);
      iInnerPasses++;
    }
    while (iInnerPasses < 5 && iTotalRec > 0);

    if (!qOK) {
      // Split edges.
      iSplittings++;
      qOK = qSplitEdges(pVM, pSM, sCESubSegs, sCEInFacets, sCFSubFacets,
			iSplittings);
    }
  }
  if (!qOK) {
#ifndef NDEBUG
    vWriteFile_AnnotatedSurface(*pSM, "failed", "", sCESubSegs, sCEInFacets,
				sCFSubFacets, pVM);
#endif
    vFoundBug("creation of constrained 3D mesh from boundary data");
    exit(111);
  }
}

// static void vOldSurfaceRecovery(VolMesh * const pVM,
// 				SurfMesh * const pSM)
// {
//   //@@ Take a census to see which edges and faces need work to recover
//   std::set<ConstrainEdge> sCESubSegs;
//   std::set<ConstrainEdge> sCEInFacets;
//   std::set<ConstrainFace> sCFSubFacets;

//   // Find out how many bdry (sub)segments are not currently in the mesh.
//   bool qOK = qIsBdryConstrained(pVM, pSM, sCESubSegs,
// 				sCEInFacets, sCFSubFacets);

//   int iOuterPasses = 0;
//   int iFaceRec3To2 = 0, iFaceRecEdge = 0;
//   while (!qOK && iOuterPasses < 1) {
//     iOuterPasses++;
//     vMessage(2, "Recovering surface edges...\n");
//     vRecoverAllEdges(pVM, pSM, sCESubSegs, sCEInFacets, sCFSubFacets);

//     if (!sCESubSegs.empty()) {
//       vWriteFile_AnnotatedSurface(*pSM, "failed", ".surf",
// 				  sCESubSegs, sCEInFacets, sCFSubFacets,
// 				  pVM);
//       vMessage(0, "Failed to recover an input bdry edge.\n");
//       vFatalError("Increasing iMaxSplittings in " __FILE__ " -may- help.",
//  		  "3D constrained Delaunay tetrahedralization.");
//     }

//     if (!sCEInFacets.empty()) {
//       vMessage(2, "%d %s\n", static_cast<int>(sCEInFacets.size()),
// 	       "generated surface edges within input facets are not present.");
//       vMessage(2, "If all goes well, facet recovery will recover these.\n");
//     }

//     if (!sCFSubFacets.empty()) {
//       vMessage(2, "Recovering surface facets...\n");
//       // If all edges weren't recovered, there's a definite problem
//       assert(sCESubSegs.empty());
//       //@@ Recover faces
//       pVM->vDisallowEdgeSwapping();
//       int iTotalPasses = 0, iChange;
// #ifndef NDEBUG
//       vMessage(4, "Faces before attempt to recover: \n");
//       for (GR_index_t i_ = 0; i_ < pVM->iNumFaces(); i_++) {
// 	Face *pF = pVM->pFFace(i_);
// 	GR_index_t iVA = pVM->iVertIndex(pF->pVVert(0));
// 	GR_index_t iVB = pVM->iVertIndex(pF->pVVert(1));
// 	GR_index_t iVC = pVM->iVertIndex(pF->pVVert(2));
// 	vMessage(4, "%4u %4u %4u %4u\n", i_, iVA, iVB, iVC);
//       }
// #endif

//       do {
// 	iTotalPasses++;
// 	do {
// 	  iChange = 0;
// 	  vRecoverFaces(pVM, pSM, sCFSubFacets, iChange, iFaceRec3To2,
// 			iFaceRecEdge, false);
// 	} while (iChange && !sCFSubFacets.empty());
// 	if (iChange == 0) pVM->vAllowEdgeSwapping();
// 	qOK = qIsBdryConstrained(pVM, pSM, sCESubSegs, sCEInFacets,
// 				 sCFSubFacets);
//       } while ((!sCFSubFacets.empty()) &&
// 	       (sCESubSegs.empty() && sCEInFacets.empty()) &&
// 	       iTotalPasses < 10);
//     } // Had to recover faces as a separate pass
//   }
//   vMessage(2, "Face recovery required %d 3-to-2 swaps and %d full edge swaps.\n",
// 	   iFaceRec3To2, iFaceRecEdge);
//   if (!sCFSubFacets.empty() || !sCEInFacets.empty() || !sCESubSegs.empty()) {
//     vWriteFile_AnnotatedSurface(*pSM, "failed", ".surf",
// 				sCESubSegs, sCEInFacets, sCFSubFacets,
// 				pVM);
//     vFoundBug("creation of constrained 3D mesh from boundary data");
//   }
// }

void
VolMesh::recoverConstrainedSurface(SurfMesh& SM, const int iNumOrigVerts,
				   const int iNumBBoxVerts)
{
  // You can't clear a std::stack, so this is the alternative.
  while (!stackSplitEdges.empty())
    stackSplitEdges.pop();

  recoverSurface(this, &SM);

  //@@ Strip stuff outside the domain
  logMessage(1, "Removing extraneous tets and cleaning up...\n");

  deleteExternal(&SM, iNumOrigVerts, iNumBBoxVerts);
  if (!isValid())
    vFoundBug("creation of constrained 3D mesh from boundary data");

  // Mark vertices on internal boundaries.
  // Changed by SG 09/2008. Previous version commented out below.
  for (GR_index_t iIntBFace = 0; iIntBFace < getNumIntBdryFaces();
      ++iIntBFace) {
    BFace* pBF = getIntBFace(iIntBFace);
    if (pBF->isDeleted())
      continue;
    assert(pBF->getType() == Cell::eIntTriBFace);
    assert(pBF->getNumFaces() == 2 && pBF->getNumVerts() == 3);
    pBF->getFace(0)->setFaceLoc(Face::eBdryTwoSide);
    pBF->getFace(1)->setFaceLoc(Face::eBdryTwoSide);
    for (int iVert = 0; iVert < 3; ++iVert) {
      Vert* pV = pBF->getVert(iVert);
      if (pV->getVertType() == Vert::eBdry)
	pV->setType(Vert::eBdryTwoSide);
    }
  }

  // The old piece of code which messes up vertices of type eBdryApex and eBdryCurve
  //   for (int iFace = 0; iFace < iNumFaces(); iFace++) {
  //     Face *pF = pFFace(iFace);
  //     if (pF->iFaceLoc() == Face::eBdryTwoSide) {
  //       pF->pVVert(0)->vSetType(Vert::eBdryTwoSide);
  //       pF->pVVert(1)->vSetType(Vert::eBdryTwoSide);
  //       pF->pVVert(2)->vSetType(Vert::eBdryTwoSide);
  //     }
  //   }

  // Mark vertices on external boundaries.  Must be done after internal
  // boundaries.
  for (GR_index_t iBFace = 0; iBFace < getNumBdryFaces(); iBFace++) {
    BFace *pBF = getBFace(iBFace);
    if (pBF->isDeleted())
      continue;
    if (pBF->getType() == Cell::eTriBFace
	|| pBF->getType() == Cell::eQuadBFace) {
      Face *pF = pBF->getFace();
      pF->setFaceLoc(Face::eBdryFace);
      assert(pF->isLocked());
      for (int iV = 0; iV < 3; iV++) {
	Vert *pV = pF->getVert(iV);
	if (pV->getVertType() == Vert::eBdryTwoSide) {
	  pV->setType(Vert::eBdryCurve);
	}
	int iVertInd = getVertIndex(pV);
	if (iVertInd < iNumOrigVerts) {
	  // Original verts can be tagged with the type that they have
	  // in the surface mesh.
	  pV->setType(SM.getVert(iVertInd)->getVertType());
	}
	assert(
	    pV->getVertType() != Vert::eInterior
		&& pV->getVertType() != Vert::eUnknown
		&& pV->getVertType() != Vert::eInvalid);
      } // end loop over this face
    } // Done dealing with external BFace's
  } // end loop over BFace's

  if (areBdryChangesAllowed()) {
    // Unlock all the bdry faces, which were presumably locked before this.
    for (int iFace = getNumFaces() - 1; iFace >= 0; iFace--) {
      getFace(iFace)->unlock();
    }
  }

  logMessage(1, "Surface recovery added %zd bdry Steiner points.\n",
	     stackSplitEdges.size());
  int iRemoved = -1, iRemoved2 = -1, iPass = 0, iRemainingInt = 0;
  while (!stackSplitEdges.empty() && (iRemoved != 0 || iRemoved2 != 0)) {
    iPass++;
    iRemoved = iRemoved2 = 0;
    int iAttempted = 0, iTotalSwaps = 0, iRemovedInt = 0;
    std::stack<SplitEdge> stackSplitEdgesTemp;
    while (!stackSplitEdges.empty()) {
      SplitEdge SE = stackSplitEdges.top();
      stackSplitEdges.pop();
      Vert *pV = SE.pVVolNew;
      assert(!pV->isDeleted() && pV->isDeletionRequested());

      // Remove as many vertices as possible by simple contraction.
      iAttempted++;
      int iNewSwaps = 0;
      Vert *apVTargets[] =
	{ SE.pV0, SE.pV1 };
      bool qSucc = removeVert(pV, iNewSwaps, apVTargets, 2);
      iTotalSwaps += iNewSwaps;
      if (qSucc) {
	// This would be the place to remove the surface point if we cared.
	// #ifndef NDEBUG
	//        // This assertion is incorrect occasionally, depending on the
	//        // order in which vertices are removed.
	// 	Face *pF = findCommonFace(SE.pV0, SE.pV1, true);
	// 	assert(pF->qValid());
	// #endif
	iRemoved++;
      }
      else {
	// Since the simple thing failed, try something more
	// sophisticated.
	bool qDelInt;
	qSucc = qRemoveBdrySteinerPoint(SE, *this, iNewSwaps, qDelInt);
	iTotalSwaps += iNewSwaps;
	if (qSucc) {
	  // This would be the place to remove the surface point if we cared.
	  // #ifndef NDEBUG
	  //        // This assertion is incorrect occasionally, depending on the
	  //        // order in which vertices are removed.
	  // 	  Face *pF = findCommonFace(SE.pV0, SE.pV1, true);
	  // 	  assert(pF->qValid());
	  // #endif
	  iRemoved2++;
	  if (qDelInt)
	    iRemovedInt++;
	  else
	    iRemainingInt++;
	}
	else {
	  // Re-queue this point so we don't lose it.
	  stackSplitEdgesTemp.push(SE);
	}
      }
    }
    stackSplitEdges = stackSplitEdgesTemp;
    if (iAttempted > 0) {
      logMessage(2, "Removed %4d verts by contraction\n", iRemoved);
      logMessage(2, "        %4d verts using interior Steiner points\n",
		 iRemoved2);
      logMessage(2, "        %4zd boundary Steiner points left\n",
		 stackSplitEdges.size());
      logMessage(2, "        %4d of the resulting interior Steiner points\n",
		 iRemovedInt);
      logMessage(2, "        %4d interior Steiner points left\n",
		 iRemainingInt);
      logMessage(2, "Performed %d swaps in the process.\n", iTotalSwaps);
    }
  }

  //   int iTotalSwaps = 0;
  //   int iAddlRemovals = 0;
  //   for (GR_index_t i = 0; i < iNumVerts(); i++) {
  //     Vert *pV = pVVert(i);
  //     if (pV->qDeletionRequested() && !pV->qDeleted()) {
  //       // Found an interior Steiner point to remove.
  //       int iRecentSwaps;
  //       if (qRemoveVert(pV, iRecentSwaps)) {
  // 	iTotalSwaps += iRecentSwaps;
  // 	iAddlRemovals++;
  // 	iRemainingInt--;
  //       }
  //     }
  //   }
  //   vMessage(2, "Removed %d addl interior Steiner points (%d swaps); %d remaining\n",
  // 	   iAddlRemovals, iTotalSwaps, iRemainingInt);
}

VolMesh* createVolMeshFromBdryFile(char strBaseFileName[FILE_NAME_LEN],
		int iBdryOptLevel, int iQualMeasure, double dResolution,
		double dGrading, Bdry3D*& pB3D) {
	// Read, validate, and optimize boundary rep
	pB3D = new Bdry3D(strBaseFileName, iBdryOptLevel);
	// Generate a surface mesh from the boundary rep
	SurfMesh SM(pB3D, iQualMeasure);
	SM.markSmallAngles();
	int iSurfaceVerts = SM.iNumVerts();
	int iSurfaceTris = SM.getNumTriCells();
	int iNonManifoldEdges = SM.getNumMultiEdges();
	int iHangingEdges = SM.getNumBdryFaces();
	// Create an initial tetrahedralization of that surface mesh
	VolMesh* pVMesh = new VolMesh(SM, iQualMeasure, 1. / dResolution, 1. / dGrading);
	int iBoundarySteiner = SM.iNumVerts() - iSurfaceVerts - 8;
	int iInteriorSteiner = pVMesh->getNumVerts() - iSurfaceVerts;
	assert(pVMesh->isValid());
    logMessage(1, "Mesh initialization summary:\n");
    logMessage(1, "  Initial surface mesh had %d verts, %d triangles.\n",
	       iSurfaceVerts, iSurfaceTris);
	logMessage(1, "Number of tets generated initially:  %u\n",
			pVMesh->getNumCells());
    logMessage(1, "  Non-manifold edges: %d.  Hanging edges: %d.\n",
	       iNonManifoldEdges, iHangingEdges);
    logMessage(1,
	       "  Boundary Steiner points: %d.  Interior Steiner points: %d.\n",
	       iBoundarySteiner, iInteriorSteiner);	logMessage(1, "Number of interior Steiner points:  %u\n", iInteriorSteiner);
	logMessage(1, "Number of boundary Steiner points:  %u\n", iBoundarySteiner);
	return pVMesh;
}

void VolMesh::initVolMeshFromSurfMesh(SurfMesh& SM) {
	copyVertData(*(SM.pECGetVerts()));
	m_Bdry3D = SM.getBoundaryData();
	// Set up an outer triangulation of a big brick
	double dXMin, dXMax, dYMin, dYMax, dZMin, dZMax;
	dXMin = dYMin = dZMin = 1.e100;
	dXMax = dYMax = dZMax = -1.e100;
	int iVert, iNVerts = getNumVerts();
	for (iVert = 0; iVert < iNVerts; iVert++) {
		double dX = getVert(iVert)->x();
		double dY = getVert(iVert)->y();
		double dZ = getVert(iVert)->z();
		dXMin = min(dXMin, dX);
		dXMax = max(dXMax, dX);
		dYMin = min(dYMin, dY);
		dYMax = max(dYMax, dY);
		dZMin = min(dZMin, dZ);
		dZMax = max(dZMax, dZ);
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

	double dXLo = dXMean - 3 * dDX;
	double dXHi = dXMean + 3 * dDX;
	double dYLo = dYMean - 3 * dDY;
	double dYHi = dYMean + 3 * dDY;
	double dZLo = dZMean - 3 * dDZ;
	double dZHi = dZMean + 3 * dDZ;
	int iNumOrigVerts = getNumVerts();
	initMeshInBrick(dXLo, dXHi, dYLo, dYHi, dZLo, dZHi);
	logMessage(1,
			"Inserting surface points to create unconstrained volume mesh.\n");
	// Do Delaunay insertion of all the existing verts in the mesh.  This
	// is a little weird, because the vert data is already -stored- in the
	// mesh, it just isn't currently connected to anything.
	GRUMMP::DelaunaySwapDecider3D Del3D(true);
	GRUMMP::SwapManager3D DelSwap(&Del3D, this);
	pSI3D = new GRUMMP::SwappingInserter3D(this, &DelSwap);
	pSI3D->setForcedInsertion(true);
	setSwapType(eDelaunay);
	int iStatusInterval = max(25, iNumOrigVerts / 10);
	std::set<Vert*> spVInsertLater;
	for (iVert = 0; iVert < iNumOrigVerts; iVert++) {
		Vert* pV = getVert(iVert);
		// Don't insert this point yet if it isn't connected to anything in
		// the surface mesh.
		Vert* pVS = SM.getVert(iVert);
		if (pVS->getNumFaces() == 0) {
			spVInsertLater.insert(pVS);
			continue;
		}
		// Make sure that the small angle flag is copied over from the
		// surface mesh.
		pV->markAsSmallAngleVert(pVS->isSmallAngleVert());
		// Insert, then swap.  Force insertion at the specified point and
		// use the existing vertex location.
		int iC = 0;
		while (getCell(iC)->isDeleted())
			iC++;
		const double* adCoords = pV->getCoords();
		Vert* newVert = pSI3D->insertPoint(adCoords, getCell(iC), pV);
		assert(newVert->isValid());
		if ((iVert + 1) % iStatusInterval == 0) {
			logMessage(2, "Inserted %5d verts out of %5d; made %5d swaps.\n",
					iVert + 1, iNumOrigVerts, pSI3D->getNumSwapsDone());
		}
		//TODO: fix this hack, which lets it find set up internal boundaries correctly
		pV->setType(Vert::eBdryApex);
	}
	setAllHintFaces();
	setVertFaceNeighbors();
	assert(isValid());
	bool qRealRecur = isSwapRecursionAllowed();
	disallowSwapRecursion();
	logMessage(1, "Recovering surface mesh...\n");
	// Add some dummy verts to surface mesh so that identical surface and
	// volume mesh indices are identical points.  Mark them for cleanup
	// before we exit.  These can be marked as already deleted, since
	// nothing will ever connect to them anyway.
	for (int ii = 0; ii < 8; ii++) {
		(void) (SM.createVert(0, 0, 0));
	}
	{
		GRUMMP::SwapManager3D* oldSwapper = pSI3D->changeSwapManager(NULL);
		pSI3D->setForcedInsertion(true);
		recoverConstrainedSurface(SM, iNumOrigVerts);
		pSI3D->changeSwapManager(oldSwapper);
	}
	// No purge has taken place at all yet so that the previous recovery
	// stuff would work.  Here it is.
	purgeAllEntities();
	allowSwapRecursion();
	// Swap to ensure (?) a Delaunay mesh!  But first empty out the queue, so
	// that it won't contain any faces that have already been deleted.
	DelSwap.clearQueue();
	logMessage(3, "Swapping towards a Delaunay mesh.\n");
	DelSwap.swapAllFaces(false);
	//   // Now try again to remove interior vertices.
	//   vRemoveTaggedVerts();
	// Now insert all the vertices that were removed from the bdry rep
	// before meshing, as well as verts that were extras to begin with.
	iStatusInterval = max(25, int(spVInsertLater.size() / 10));
	int iStatusValue = iStatusInterval;
	int iVertsInserted = 0;
	//TODO: FIGURE THIS OUT
	std::multimap<double, Vert*> mapVertInsertOrder;
	{
		std::set<Vert*>::iterator iter = spVInsertLater.begin();
		std::set<Vert*>::iterator iterEnd = spVInsertLater.end();
		int aiCounts[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		while (iter != iterEnd) {
			Vert* pVS = *iter;
			++iter;
			int iOffset = m_Bdry3D->iNumPatchesContaining(pVS);
			aiCounts[iOffset]++;
			double dPriority = -iOffset + drand48();
			mapVertInsertOrder.insert(std::make_pair(dPriority, pVS));
		}
		logMessage(1,
				"Out of %zd verts, %d interior, %d bdry, %d curve, %d non-manifold / vertex\n",
				spVInsertLater.size(), aiCounts[0], aiCounts[1], aiCounts[2],
				aiCounts[3]);
	}
	int iC = 0;
	while (getCell(iC)->isDeleted())
		iC++;
	Cell* pC = getCell(iC);
	std::multimap<double, Vert*>::iterator iter = mapVertInsertOrder.begin();
	std::multimap<double, Vert*>::iterator iterEnd = mapVertInsertOrder.end();
	while (iter != iterEnd) {
		Vert* pVS = iter->second;
		++iter;
		//     vMessage(2, "Inserting vert %d; priority %f (%f, %f, %f)\n",
		// 	     iVertsInserted+1, iter->first, pVS->dX(), pVS->dY(), pVS->dZ());
		assert(pVS->getNumFaces() > 0 && !pVS->isDeleted());
		const double* adCoords = pVS->getCoords();
		logMessage(4, "Attempting to insert a vertex at (%10G %10G %10G)\n",
				adCoords[0], adCoords[1], adCoords[2]);
		// This is a really slow way to do the re-insertion, but it -will- work.
		bool qHit = false;
		double adBary[4], dBestBary = -1.e20;
		//		for (GR_index_t iC = 0; !qHit && iC < getNumTetCells(); iC++) {
		//			pC = getCell(iC);
		bool status;
		iC = 0;
		while (pC->isDeleted())
			pC = getCell(iC++);
		pC = findCell(adCoords, pC, status);
		//			if (pC->isDeleted()) continue;
		TetCell* pTC = dynamic_cast<TetCell*>(pC);
		assert(pTC != NULL);
		pTC->calcBarycentricCoords(adCoords, adBary);
		double dMinBary = min(min(adBary[0], adBary[1]),
				min(adBary[2], adBary[3]));
		if (dMinBary > dBestBary)
			dBestBary = dMinBary;

		qHit = (iFuzzyComp(dMinBary, 0) != -1);
		//		}
		if (!qHit) {
			logMessage(0, "Best bary: %.16f\n", dBestBary);
			logMessage(0, "Attempting to insert a vertex at (%10G %10G %10G)\n",
					adCoords[0], adCoords[1], adCoords[2]);
			vFatalError("Tried to insert a point outside the mesh",
					"VolMesh constructor");
		}
		qHit = false;
		// This vertex may already have been inserted during surface
		// recovery, in which case we don't need to worry about it.
		double dMaxBary = max(max(adBary[0], adBary[1]),
				max(adBary[2], adBary[3]));
		if (iFuzzyComp(dMaxBary, 1) == 0) {
			continue;
		}
		Vert* newVert = pSI3D->insertPoint(adCoords, pC);
		assert(newVert->isValid());
		iVertsInserted++;
		if (iVertsInserted > iStatusValue) {
			logMessage(1,
					"Inserted %5d extra verts from the surface mesh, and made %5d swaps.\n",
					iVertsInserted, pSI3D->getNumSwapsDone());
			iStatusValue += iStatusInterval;
		}
		//		if(iVertsInserted < 2000){
		//			char tempname[FILE_NAME_LEN];
		//					sprintf(tempname,"tempMesh%d",iVertsInserted);
		//					writeVTKLegacyWithoutPurge(*this,tempname);
		//		}
	}
	setSwapType(eDelaunay);
	if (iVertsInserted > 0)
		DelSwap.swapAllFaces();

	setAllHintFaces();
	setVertFaceNeighbors();
	assert(isValid());
	setSwapType(eMinSine);
	if (qRealRecur) {
		allowSwapRecursion();
	} else {
		disallowSwapRecursion();
	}
	// Set up length scale --- shouldn't need to do this, we should be initializing length scale later with a Length3D object before using it.
	// initializeSubsegs();
	// GRUMMP::Length3D* length3D = new GRUMMP::Length3D(this, m_resolution,
	// 		m_LipschitzAlpha);
	// length3D->set_subseg_map(this->getSubsegMap());
	//initLengthScale(m_resolution, m_LipschitzAlpha);
	// Clean up
	purgeAllEntities();
	SM.purgeAllEntities();
}

//@ Triangulation of a generalized boundary geometry
VolMesh::VolMesh(SurfMesh& SM, const int iQualMeas, const double dRes,
		 const double dGrade, const double dMaxAngle) :
    Mesh(), m_OKToSwapSurfaceEdges(true), m_doEdgeSwaps(
	true), m_strictPatchChecking(true), m_maxAngleForSurfSwap(dMaxAngle), m_ECTriF(), m_ECQuadF(), m_ECTriBF(), m_ECQuadBF(), m_ECTet(), m_ECPyr(), m_ECPrism(), m_ECHex(), m_ECIntTriBF(), m_ECIntQuadBF(), m_Bdry3D(), m_bdryFromMesh(
	false), m_SubsegMap(new SubsegMap())
{
  m_resolution = dRes;
  m_LipschitzAlpha = dGrade;
  m_qual = new Quality(this, iQualMeas);
  m_Bdry3D = nullptr;

	initVolMeshFromSurfMesh(SM);
}
