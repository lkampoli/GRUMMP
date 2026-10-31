#include <set>
#include <map>
#include "GR_Bdry3D.h"
#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Sort.h"
#include "GR_SurfMesh.h"
#include "GR_Bdry2D.h"
#include "GR_Bdry3D.h"
#include "GR_InsertionManager.h"
#include "GR_SwapManager.h"
#include "GR_SwapDecider.h"
#include "GR_Geometry.h"

using std::map;

// This is here so that configure can find the library.
extern "C"
{
  void
  vGRUMMP_Surf_Lib(void)
  {
  }
}

// Obsolete, apparently....
void
SurfMesh::setVertFaceNeighbors()
{
  Mesh::setVertFaceNeighbors();
  //   vMessage(2, "Adding multi-edges as neighbors...");

  //   int i;
  //   for (i = 0; i < iNumMultiEdges(); i++) {
  //     MultiEdge* pME = pMEMultiEdge(i);
  //     if (pME->qDeleted()) continue;
  //     for (int ii = pME->iNumVerts() - 1; ii >= 0; ii --) {
  //       Vert* pV = pME->pVVert(ii);
  //       if (pV != pVInvalidVert) {
  // 	pV->vAddFace(pME);
  //       }
  //     }
  //   }
  //   vMessage(2, "done\n");
}

// //@ Copy constructor
// SurfMesh::SurfMesh(const SurfMesh& SM)
//      // Attach to the point set and connectivity, and compute the vertex
//      // hint data.  NOTE:  THE CONNECTIVITY MUST NOT BE DESTROYED BY THE
//      // ROUTINE WHICH CALLS THIS CONSTRUCTOR.  The assumption is that
//      // the connectivity is created with 'new' and not 'delete'-ed until
//      // the SurfMesh dies.
//   : Mesh2D(SM), ECME(), iOnMultiEdge(0), pB3D(), apBPatch(NULL),
//     iBPatchSize(0)
// {
//   assert(0);
// }

// SurfMesh& SurfMesh::operator=(const SurfMesh& SM)
// {
//   assert(0);
//   // Do nothing for self-assignment
//   if (this != &SM) {
//     ECVerts  = SM.ECVerts;
//     pQ       = SM.pQ;
//     ECEdgeF  = SM.ECEdgeF;
//     ECEdgeBF = SM.ECEdgeBF;
//     ECTri    = SM.ECTri;
//     ECQuad   = SM.ECQuad;
//     pB3D     = SM.pB3D;
//     ECME     = SM.ECME;

//     vSetAllHintFaces();
//   }
//   if (!(qValid()))
//     vFoundBug("surface mesh copy");
//   return (*this);
// }

//@ Create an initial surface triangulation from a Bdry3D object.

// For now, this code is set up strictly for getting initial
// triangulations from polyhedral input data.  As such, each surface
// patch is triangulated (constrained Delaunay) individually, and the
// SurfMesh constructor combines the topology info into a mesh
// (probably by using a cell-vertex to face-cell data rep conversion).

// In the longer term, this constructor will have to negotiate
// consistent discretizations of patch bounding curves and pass
// appropriate subsets of this info to the patches so they can
// discretize themselves properly.  Each patch (or possibly composite of
// patches, a la the 2D interpolated spline stuff) should triangulate
// itself on demand, returning a set of points and a connectivity.  In
// this scenario, verts may well be added by patch triangulation, and
// the constructor will have to handle this.
SurfMesh::SurfMesh(Bdry3D * const pB3DIn, const int iQualMeas) :
    Mesh2D(iQualMeas), m_ECME(), m_Bdry3D(pB3DIn), m_bdryPatches(NULL), m_numBdryPatches(
	0), m_SubsegMap(NULL), m_surfaces(NULL), m_numSurfaces(0)
{
  EntContainer<Vert> *pEC = m_Bdry3D->pECGetVerts();
  for (GR_index_t i = 0; i < pEC->lastEntry(); i++) {
    Vert *pVOld = pEC->getEntry(i);
    Vert *pVNew = m_ECVerts.getNewEntry();
    pVNew->setCoords(3, pVOld->getCoords());
    pVNew->copyAllFlagsFrom(pVOld);
    // Don't copy vertex data any deeper than this.
  }
  // Count up, in advance, how many cells there will be in the final
  // surface mesh.  Can't just use the Euler formula, because we don't
  // know how many closed surfaces and/or handles there are in the
  // input.
  GR_index_t iP, iTotalCells = 0;
  for (iP = 0; iP < m_Bdry3D->iNumBdryPatches(); iP++) {
    int iNumSegs = m_Bdry3D->pBdry(iP)->iNumSegments();
    int iNumKids = m_Bdry3D->pBdry(iP)->iNumChildren();
    iTotalCells += iNumSegs + 2 * (iNumKids - 1);
  }
  m_bdryPatches = new BdryPatch3D*[iTotalCells];
  m_numBdryPatches = iTotalCells;

  // Set up a place to put surface triangulation data as it comes
  // through.
  EntContainer<TriCellCV> ECTriCV;

  logMessage(2, "Setting up surface mesh based on geometric data...\n");
  // Turn messaging down a notch during surface triangulation.
  //  iMessageStdoutLevel++;
  //  iMessageFileLevel++;
  // For every patch, triangulate it, and add its triangles to the list
  // for future conversion.
  GR_index_t iCumCell = 0;
  for (iP = 0; iP < m_Bdry3D->iNumBdryPatches(); iP++) {
    BdryPatch3D *pBPatch = m_Bdry3D->pBdry(iP);
    // Can only handle polygons without negotiating curve discretization
    // in advance, so check for that.
    assert(pBPatch->eType() == BdryPatch3D::ePolygon);
    BdryPolygon *pBP = dynamic_cast<BdryPolygon*>(pBPatch);
    //		if (pBP->iNumSegments() == 3 && pBP->iNumChildren() == 0) {
    //			// We have a triangle with nothing inside it.  This is an easy case.
    //			Vert *pV0 = getVert(m_Bdry3D->iVertIndex(pBP->pVVert(0)));
    //			Vert *pV1 = getVert(m_Bdry3D->iVertIndex(pBP->pVVert(1)));
    //			Vert *pV2 = getVert(m_Bdry3D->iVertIndex(pBP->pVVert(2)));
    //			TriCellCV *pCCV = ECTriCV.getNewEntry();
    //			if (pBPatch->iLeftRegion() == iOutsideRegion) {
    //				pCCV->vAssign(pV0, pV2, pV1);
    //			}
    //			else {
    //				pCCV->vAssign(pV0, pV1, pV2);
    //			}
    //#ifndef NDEBUG
    //			double adNormPatch[3], adNormVerts[3];
    //			calcUnitNormal(pV0->getCoords(), pV1->getCoords(), pV2->getCoords(),
    //					adNormVerts);
    //			// First arg is a dummy anyway for polygons...
    //			pBP->vUnitNormal(adNormVerts, adNormPatch);
    //			double dDot = dDOT3D(adNormPatch, adNormVerts);
    //			//       if (pBPatch->iLeftRegion() == iOutsideRegion)
    //			// 	vSCALE3D(adNormVerts, -1);
    //			logMessage(4, "Cell: %3u  Verts: %3u, %3u, %3u.\n\t(%7f %7f %7f) Dot: %6f  Regs: %3d %3d\n",
    //					iCumCell, iVertIndex(pV0), iVertIndex(pV1), iVertIndex(pV2),
    //					adNormVerts[0], adNormVerts[1], adNormVerts[2],
    //					dDot, pBPatch->iLeftRegion(), pBPatch->iRightRegion());
    //			assert((dDot > 0 && pBPatch->iRightRegion() == iOutsideRegion) ||
    //					(dDot < 0 && pBPatch->iLeftRegion()  == iOutsideRegion) ||
    //					(pBPatch->iRightRegion() != iOutsideRegion &&
    //							pBPatch->iLeftRegion()  != iOutsideRegion));
    //#endif
    //			setBdryPatch(iCumCell, pBPatch);
    //			iCumCell++;
    //		}
    //
    //		else
    {
      // Something that we actually have to triangulate...
      Bdry2D B2D;
      pBP->B2DProjection(*m_Bdry3D, B2D);
      Mesh2D M2D(B2D);
      double adOffset[3], adBasis1[3], adBasis2[3];
      pBP->vFindProjectionBasis(adOffset, adBasis1, adBasis2);
      GRUMMP::WatsonInserter WI(&M2D);
      for (GR_index_t iVert = 0; iVert < getNumVerts(); iVert++) {
	Vert *pVS = getVert(iVert);

	if (pBP->qDisjointFrom(pVS) == 0) {
	  // need to project to the boundary;

	  const double * const adLoc3D = pVS->getCoords();

	  // To project, subtract the offset distance (adOffset).  Then take
	  // dot products between that vector (which lies in a plane through
	  // the origin) and the two basis vectors (both of which lie in that
	  // plane).
	  double adInPlane[] = adDIFF3D(adLoc3D, adOffset);

	  double adLoc2D[] =
	    { dDOT3D(adInPlane, adBasis1), dDOT3D(adInPlane, adBasis2) };

	  Cell * cellGuess = pCInvalidCell;
	  int inTri = 0;
	  for (GR_index_t iC = 0; iC < M2D.getNumCells(); iC++) {
	    inTri = isInTri(adLoc2D, M2D.getCell(iC)->getVert(0)->getCoords(),
			    M2D.getCell(iC)->getVert(1)->getCoords(),
			    M2D.getCell(iC)->getVert(2)->getCoords());
	    if (inTri > 0) {
	      cellGuess = M2D.getCell(iC);
	      break;
	    }
	  }
	  if (inTri == 1) {
	    // done, vert already exists.
	  }
	  else {
	    if (inTri == 2) {
	      if (cellGuess->getFace(0)->isBdryFace()
		  && checkOrient2D(
		      cellGuess->getFace(0)->getVert(0)->getCoords(),
		      cellGuess->getFace(0)->getVert(1)->getCoords(), adLoc2D)
		      == 0) {
		BFace * bface =
		    dynamic_cast<BFace*>(cellGuess->getFace(0)->getOppositeCell(
			cellGuess));
		assert(bface);
		WI.computeHull(adLoc2D, cellGuess, bface);
		Vert * newVert = WI.insertPointInHull();
		newVert->setType(Vert::eBdryCurve);
	      }
	      else if (cellGuess->getFace(1)->isBdryFace()
		  && checkOrient2D(
		      cellGuess->getFace(1)->getVert(0)->getCoords(),
		      cellGuess->getFace(1)->getVert(1)->getCoords(), adLoc2D)
		      == 0) {
		BFace * bface =
		    dynamic_cast<BFace*>(cellGuess->getFace(1)->getOppositeCell(
			cellGuess));
		assert(bface);
		WI.computeHull(adLoc2D, cellGuess, bface);
		Vert * newVert = WI.insertPointInHull();
		newVert->setType(Vert::eBdryCurve);
	      }
	      else if (cellGuess->getFace(2)->isBdryFace()
		  && checkOrient2D(
		      cellGuess->getFace(2)->getVert(0)->getCoords(),
		      cellGuess->getFace(2)->getVert(1)->getCoords(), adLoc2D)
		      == 0) {
		BFace * bface =
		    dynamic_cast<BFace*>(cellGuess->getFace(2)->getOppositeCell(
			cellGuess));
		assert(bface);
		WI.computeHull(adLoc2D, cellGuess, bface);
		Vert * newVert = WI.insertPointInHull();
		newVert->setType(Vert::eBdryCurve);
	      }
	      else {
		WI.computeHull(adLoc2D, cellGuess);
		Vert * newVert = WI.insertPointInHull();
		newVert->setType(Vert::eBdry);
	      }
	    }
	    else {
	      WI.computeHull(adLoc2D, cellGuess);
	      Vert * newVert = WI.insertPointInHull();
	      newVert->setType(Vert::eBdry);
	    }
	  }
	}
      }

//						char tempname[FILE_NAME_LEN];
//						sprintf(tempname,"tempMesh2D%d",iP);
//						writeVTKLegacyWithoutPurge(M2D,tempname);

      // Now the Mesh2D has the triangulation of the surface polygon, but
      // the indices of the verts do not match indices in the surface
      // mesh!  So we have to find the real vertices
      // So copy that data over.
      for (GR_index_t iC = 0; iC < M2D.getNumCells(); iC++) {
	Cell *pC = M2D.getCell(iC);
	Vert *pV2D0 = pC->getVert(0);
	Vert *pV2D1 = pC->getVert(1);
	Vert *pV2D2 = pC->getVert(2);

	double adLoc3D0[3], adLoc3D1[3], adLoc3D2[3];
	adLoc3D0[0] = (+pV2D0->x() * adBasis1[0] + pV2D0->y() * adBasis2[0]
	    + adOffset[0]);
	adLoc3D0[1] = (+pV2D0->x() * adBasis1[1] + pV2D0->y() * adBasis2[1]
	    + adOffset[1]);
	adLoc3D0[2] = (+pV2D0->x() * adBasis1[2] + pV2D0->y() * adBasis2[2]
	    + adOffset[2]);

	adLoc3D1[0] = (+pV2D1->x() * adBasis1[0] + pV2D1->y() * adBasis2[0]
	    + adOffset[0]);
	adLoc3D1[1] = (+pV2D1->x() * adBasis1[1] + pV2D1->y() * adBasis2[1]
	    + adOffset[1]);
	adLoc3D1[2] = (+pV2D1->x() * adBasis1[2] + pV2D1->y() * adBasis2[2]
	    + adOffset[2]);

	adLoc3D2[0] = (+pV2D2->x() * adBasis1[0] + pV2D2->y() * adBasis2[0]
	    + adOffset[0]);
	adLoc3D2[1] = (+pV2D2->x() * adBasis1[1] + pV2D2->y() * adBasis2[1]
	    + adOffset[1]);
	adLoc3D2[2] = (+pV2D2->x() * adBasis1[2] + pV2D2->y() * adBasis2[2]
	    + adOffset[2]);

	// These points may only be very close approximations to
	// actually points in the patch if the patch is not quite
	// precisely coplanar.  So find the nearest point.

	int iLocalInd0 = pBP->iLocalVertIndex(adLoc3D0);
	int iLocalInd1 = pBP->iLocalVertIndex(adLoc3D1);
	int iLocalInd2 = pBP->iLocalVertIndex(adLoc3D2);
	Vert * pV0, *pV1, *pV2;
	if (iLocalInd0 == -1) {
	  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
	    Vert * pV = getVert(iV);
	    if (dDIST3D(pV->getCoords(),adLoc3D0) < 1.e-6) {
	      iLocalInd0 = iV;
	      break;
	    }
	  }
	  pV0 = getVert(iLocalInd0);
	  pV0->setType(pV2D0->getVertType());

	}
	else
	  pV0 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd0));

	if (iLocalInd1 == -1) {
	  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
	    Vert * pV = getVert(iV);
	    if (dDIST3D(pV->getCoords(),adLoc3D1) < 1.e-6) {
	      iLocalInd1 = iV;
	      break;
	    }
	  }
	  pV1 = getVert(iLocalInd1);
	  pV1->setType(pV2D1->getVertType());

	}
	else
	  pV1 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd1));
	if (iLocalInd2 == -1) {
	  for (GR_index_t iV = 0; iV < getNumVerts(); iV++) {
	    Vert * pV = getVert(iV);
	    if (dDIST3D(pV->getCoords(),adLoc3D2) < 1.e-6) {
	      iLocalInd2 = iV;
	      break;
	    }
	  }
	  pV2 = getVert(iLocalInd2);
	  pV2->setType(pV2D2->getVertType());
	}
	else
	  pV2 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd2));

	assert(iLocalInd0 != -1 && iLocalInd1 != -1 && iLocalInd2 != -1);
	//Vert *pV0 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd0));

	//				Vert *pV1 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd1));
	//				Vert *pV2 = getVert(m_Bdry3D->iVertIndex(iP, iLocalInd2));

	TriCellCV *pCCV = ECTriCV.getNewEntry();
	if (pBPatch->iLeftRegion() != iOutsideRegion) {
	  pCCV->vAssign(pV0, pV2, pV1);
	}
	else {
	  pCCV->vAssign(pV0, pV1, pV2);
	}
	setBdryPatch(iCumCell, pBPatch);
	iCumCell++;
      }
    }
  }
  assert(iCumCell == ECTriCV.lastEntry());
  // Turn messaging back up a notch.
  //  iMessageStdoutLevel--;
  //  iMessageFileLevel--;

  // Convert from cell-vert to face-cell data structure and store here
  // in the SurfMesh.
  convertFromCellVert(ECTriCV);
  setAllHintFaces();
  setVertFaceNeighbors();

  for (int iC = getNumCells() - 1; iC >= 0; iC--) {
    Cell *pC = getCell(iC);
    pC->setRegion(iDefaultRegion);
  }

  // Is the thing right side out?
  double dSize = calcInteriorSize();
  if (dSize < 0) {
    logMessage(3, "Turning surface mesh right side out.\n");
    // Need to flip the surface orientation, both for the SurfMesh and
    // the Bdry3D.
    for (GR_index_t iF = 0; iF < getNumFaces(); iF++) {
      getFace(iF)->interchangeCells();
    }
    for (int iBP = 0; iBP < m_Bdry3D->iNumPatches(); iBP++) {
      m_Bdry3D->pBdry(iBP)->vSwitchOrientation();
      m_Bdry3D->pBdry(iBP)->vFlipNormal();
    }
    dSize = calcInteriorSize();
  }
  assert(dSize > 0);
  logMessage(1, "Total volume enclosed in surface mesh: %10.5g\n", dSize);

  if (!(isValid()))
    vFoundBug("surface mesh creation from Bdry3D data");
}

MultiEdge*
SurfMesh::getNewMultiEdge()
{
  MultiEdge *pME = m_ECME.getNewEntry();
  pME->setDefaultFlags();
  return pME;
}
// lets build this from a list of facevecs
SurfMesh::SurfMesh(FaceVec faces, std::set<Subseg*> subsegs,
		   const int iQualMeas) :
    Mesh2D(iQualMeas), m_ECME(), m_Bdry3D(NULL), m_bdryPatches(NULL), m_numBdryPatches(
	0), m_SubsegMap(NULL), m_surfaces(NULL), m_numSurfaces(faces.size())
{

  m_surfaces = new RefFace*[m_numSurfaces];
  std::map<Vert*, Vert*> vertMap; // mapping from input verts to output verts
  std::map<Vert*, Vert*>::iterator vertIter;

  m_SubsegMap = new SubsegMap();
  for (FaceVec::iterator it = faces.begin(); it != faces.end(); ++it) {
    Face * face = (*it).first->getFace();
    RefFace * refFace = (*it).second;
    Vert * newVerts[3];
    for (GR_index_t iV = 0; iV < 3; iV++) {
      Vert * vert = face->getVert(iV);
      vertIter = vertMap.find(vert);
      if (vertIter == vertMap.end()) {
	newVerts[iV] = createVert(vert->getCoords());
	vertMap[vert] = newVerts[iV];
      }
      else
	newVerts[iV] = (*vertIter).second;
    }
    Cell * cell = createTriCell(newVerts[0], newVerts[1], newVerts[2], 1);
    setSurface(cell, refFace);

  }
  for (std::set<Subseg*>::iterator it = subsegs.begin(); it != subsegs.end();
      ++it) {
    Subseg * subseg = *it;

    vertIter = vertMap.find(subseg->get_beg_vert());
    assert(vertIter != vertMap.end());
    Vert * v0 = vertIter->second;
    v0->setType(vertIter->first->getVertType());

    vertIter = vertMap.find(subseg->get_end_vert());
    assert(vertIter != vertMap.end());

    Vert * v1 = vertIter->second;
    v1->setType(vertIter->first->getVertType());
    if (subseg->get_ref_edge())
      m_SubsegMap->add_subseg(
	  new Subseg(v0, v1, subseg->get_beg_param(), subseg->get_end_param(),
		     subseg->get_ref_edge()));
    else if (subseg->get_ref_face())
      m_SubsegMap->add_subseg(new Subseg(v0, v1, subseg->get_ref_face()));
    else
      m_SubsegMap->add_subseg(new Subseg(v0, v1));

    createBFace(findCommonFace(v0, v1, true));
  }
}

double
SurfMesh::calcBoundarySize() const
{
  double dRetVal = 0;
  for (int iC = getNumCells() - 1; iC >= 0; iC--)
    dRetVal += getCell(iC)->calcSize();
  return dRetVal;
}

double
SurfMesh::calcInteriorSize() const
{
  double dVolX = 0;
  // The checks for size consistency can fail when there are faces that
  // should be coplanar but just fail to be.  So don't do the checks.
  // #ifndef NDEBUG
  //   double dVolY = 0, dVolZ = 0;
  // #endif
  double adNorm[3], adCent[3];
  for (int iC = getNumCells() - 1; iC >= 0; iC--) {
    if (getCell(iC)->isDeleted())
      continue;
    // If this cell of the SurfMesh lies on an interior patch, it
    // shouldn't contribute to the integrals for interior size.
    BdryPatch *pBP = getBdryPatch(iC);
    if (pBP->qIsBoundary()) {
      Cell *pC = getCell(iC);
      pC->calcCentroid(adCent);
      double dSize = pC->calcSize();
      pBP->vUnitNormal(adCent, adNorm);
      if (pBP->iLeftRegion() != iOutsideRegion) {
	vSCALE3D(adNorm, -1);
      }

      dVolX += adCent[0] * adNorm[0] * dSize;
      // See comment above.
      // #ifndef NDEBUG
      //       dVolY += adCent[1]*adNorm[1]*dSize;
      //       dVolZ += adCent[2]*adNorm[2]*dSize;
      //       double adNorm2[3];
      //       vNormal(pC->pVVert(0)->adCoords(), pC->pVVert(1)->adCoords(),
      //	      pC->pVVert(2)->adCoords(), adNorm2);
      //       vMessage(2, "Cell: %3d  Norm: (%7f %7f %7f)  Cent: (%7f %7f %7f)\n",
      //	       iC, adNorm[0], adNorm[1], adNorm[2],
      //	       adCent[0], adCent[1], adCent[2]);
      //       vMessage(2, "   Scaled cell norm: (%7f %7f %7f)\n",
      //	       adNorm2[0]/dSize, adNorm2[1]/dSize, adNorm2[2]/dSize);
      // #endif
    }
  }
  // See comment above.
  // #ifndef NDEBUG
  //   assert(iFuzzyComp(dVolX, dVolY) == 0);
  //   assert(iFuzzyComp(dVolX, dVolZ) == 0);
  // #endif
  return dVolX;
}

bool
SurfMesh::isWatertight() const
{
  double dAreaX = 0, dAreaY = 0, dAreaZ = 0;
  double adNorm[3];
  for (int iC = getNumCells() - 1; iC >= 0; iC--) {
    getCell(iC)->calcVecSize(adNorm);
    dAreaX += adNorm[0];
    dAreaY += adNorm[1];
    dAreaZ += adNorm[2];
  }
  double dCompSize = calcBoundarySize();
  double dInvCompSize = 0.001 / dCompSize;
  dAreaX *= dInvCompSize;
  dAreaY *= dInvCompSize;
  dAreaZ *= dInvCompSize;

  return ((iFuzzyComp(dAreaX, 0.) == 0) && (iFuzzyComp(dAreaY, 0.) == 0)
      && (iFuzzyComp(dAreaZ, 0.) == 0));
}

void
SurfMesh::markSmallAngles()
// Find all small angles between bdry segments in the surface mesh and
// check whether they could lead to infinite refinement because of
// mutual encroachment.  Because of the possibility of uneven edge
// splits, the edge lengths aren't considered as a factor here.
// Instead, only the angle is used (it must be greater than 45 degrees
// to prevent problems).  Vertices where we could get infinite
// refinement (let's call this implosion) are marked as having a small
// angle; this marking gets carried over into the volume mesh.  Once
// the initial length scale is established, each edge incident on a
// small angle apex is split at a small distance from the corner.
// Later on, tets incident on the small angle corner are protected
// from quality-based refinement.
{
  logMessage(2, "Indentifying small angles in the surface mesh.\n");
  GR_index_t iNumToProtect = 0;

  for (GR_index_t iV = 0; iV < iNumVerts(); iV++) {
    // Find all incident edges.
    Vert *pV = getVert(iV);
    if (!pV->getHintFace()->isValid())
      continue;
    std::set<Face*> spFNearby;
    {
      std::set<Vert*> spVNeigh;
      std::set<Cell*> spCTmp;
      std::set<BFace*> spBFTmp;
      bool qBdryVert;
      findNeighborhoodInfo(pV, spCTmp, spVNeigh, &spBFTmp, &qBdryVert,
			   &spFNearby);
    }

    // Eliminate faces that aren't incident on pV or don't lie along
    // bdry segments.  These are easily identifiable, because they have
    // incident cells from the same bdry patch.
    std::set<Face*>::iterator iterF;
    for (iterF = spFNearby.begin(); iterF != spFNearby.end();) {
      Face *pF = *iterF;
      if (!pF->hasVert(pV)) {
	std::set<Face*>::iterator iterFTmp = iterF;
	++iterF;
	spFNearby.erase(iterFTmp);
	continue;
      }

      int iCLeft = getCellIndex(pF->getLeftCell());
      int iCRight = getCellIndex(pF->getRightCell());
      BdryPatch3D* pBPLeft = getBdryPatch(iCLeft);
      BdryPatch3D* pBPRight = getBdryPatch(iCRight);
      if (pBPLeft == pBPRight) {
	std::set<Face*>::iterator iterFTmp = iterF;
	++iterF;
	spFNearby.erase(iterFTmp);
      }
      else
	++iterF;
    }

    // For each vert, find whether any pair of edges can implode.
    const int iNumEdges = int(spFNearby.size());
    Vert **apVNeigh = new Vert*[iNumEdges];
    double *adEdgeLengths = new double[iNumEdges];
    bool *aqImplodable = new bool[iNumEdges];
    bool qAnyImplodable = false;

    int iEdge = 0;
    for (iterF = spFNearby.begin(); iterF != spFNearby.end();
	++iterF, iEdge++) {
      aqImplodable[iEdge] = false;
      Face *pF = *iterF;
      apVNeigh[iEdge] =
	  (pF->getVert(0) == pV) ? pF->getVert(1) : pF->getVert(0);
    }

    for (int iOuter = 0; iOuter < iNumEdges; iOuter++) {
      double adVecOut[] = adDIFF3D(apVNeigh[iOuter]->getCoords(),
	  pV->getCoords());
      double dLenOut = dMAG3D(adVecOut);
      adEdgeLengths[iOuter] = dLenOut;
      for (int iInner = iOuter + 1; iInner < iNumEdges; iInner++) {
	double adVecIn[] = adDIFF3D(apVNeigh[iInner]->getCoords(),
	    pV->getCoords());
	double dLenIn = dMAG3D(adVecIn);
	double dCosTheta = dDOT3D(adVecOut, adVecIn) / (dLenOut * dLenIn);
	// If the angle between the edges is less than 45 degrees, then
	// there is a chance that this pair could implode.
	if (dCosTheta >= M_SQRT1_2) {
	  double dRatio = dLenOut / dLenIn;
	  // Now adjust dRatio to be in the range [0.5,1].  This works
	  // whether dRatio > 1 or < 1.
	  double dK = floor(-log(dRatio) / M_LN2);
	  double dRatio_Star = dRatio * pow(2., dK);
	  logMessage(
	      4,
	      " Cos(theta) is %g.  Unscaled ratio is %g; scaled is %g.  This is",
	      dCosTheta, dRatio, dRatio_Star);
	  // Using iFuzzyComp here to avoid possible roundoff weirdness
	  // with lengths that differ by powers of two.
	  assert(
	      iFuzzyComp(dRatio_Star, 0.5) >= 0
		  && iFuzzyComp(dRatio_Star, 1) <= 0);
	  if (iFuzzyComp(dRatio_Star, 0.5 / dCosTheta) >= 0
	      && iFuzzyComp(dRatio_Star, dCosTheta) <= 0) {
	    // So we're in the evil range.  This pair of edges could
	    // implode.
	    logMessage(4, " evil\n");
	    aqImplodable[iOuter] = aqImplodable[iInner] = true;
	    qAnyImplodable = true;
	  }
	  else {
	    logMessage(4, " okay\n");
	  }
	}
      } // Done with inner loop for edge-to-edge comparison
    } // Done with outer loop for edge-to-edge comparison

    if (qAnyImplodable) {
      pV->markAsSmallAngleVert(true);
      iNumToProtect++;

      // Find the shortest and 2nd shortest edges in the bundle.  If the
      // shortest and 2nd shortest have a length ratio > 1.5, split all
      // edges at the shortest edge length from the common vert.  Otherwise,
      // split at half the shortest edge length.  This will maximize the
      // minimum edge length after insertion.
      double dShortestLen = LARGE_DBL, dSecondLen = LARGE_DBL;
      //       int iShortest = -1;
      //       int iSecond = -1;
      for (iEdge = 0; iEdge < iNumEdges; iEdge++) {
	if (aqImplodable[iEdge]) {
	  double dThisLen = adEdgeLengths[iEdge];
	  if (dThisLen < dShortestLen) {
	    //	    iSecond = iShortest;
	    dSecondLen = dShortestLen;
	    //	    iShortest = iEdge;
	    dShortestLen = dThisLen;
	  }
	  else if (dThisLen < dSecondLen) {
	    //	    iSecond = iEdge;
	    dSecondLen = dThisLen;
	  }
	}
      }

      assert(dShortestLen <= dSecondLen);
      //       assert(iSecond != -1 && iShortest != -1);

      double dInsertionLen = dShortestLen / 2;
      //       bool qSplitShortest = true;
      if (dSecondLen / dShortestLen > 1.5) {
	dInsertionLen = dShortestLen;
	//	qSplitShortest = false;
      }
      pV->setLengthScale(dInsertionLen);
    }
    delete[] aqImplodable;
    delete[] apVNeigh;
    delete[] adEdgeLengths;
  } // Done with this vertex
  logMessage(
      2, "Out of %u vertices, found %u that may need small angle protection.\n",
      iNumVerts(), iNumToProtect);
}

// Following is the body of SurfMesh::vProtectSmallAngles, which was
// designed to split edges so that infinite encroachment could not
// occur.  That attempt to prevent implosion failed, but I'm not yet
// ready to write it off completely by eliminating the code.

//void SurfMesh::protectSmallAngles() {
//	logMessage(2, "Protecting small angles in the surface mesh.\n");
//
//	bool qChanged = false;
//	// Loop from the last initial vertex down, so that we don't end up
//	// working on verts inserted during the process.  Mind you,
//	// newly-inserted verts should be no problem, but this is more
//	// efficient, regardless.
//	do {
//		int iNumChanges = 0;
//		qChanged = false;
//		for (int iV = iNumVerts() - 1; iV >= 0; iV--) {
//			// Find all incident edges.
//			Vert *pV = getVert(iV);
//			if (!pV->getHintFace()->isValid())
//				continue;
//			std::set<Face*> spFNearby;
//			{
//				std::set<Vert*> spVNeigh;
//				std::set<Cell*> spCTmp;
//				std::set<BFace*> spBFTmp;
//				bool qBdryVert;
//				findNeighborhoodInfo(pV, spCTmp, spVNeigh, &spBFTmp, &qBdryVert,
//						&spFNearby);
//			}
//
//			// Eliminate faces that aren't incident on pV or don't lie along
//			// bdry segments.  These are easily identifiable, because they have
//			// incident cells from the same bdry patch.
//			std::set<Face*>::iterator iterF;
//			for (iterF = spFNearby.begin(); iterF != spFNearby.end(); ++iterF) {
//				Face *pF = *iterF;
//				if (!pF->hasVert(pV)) {
//					std::set<Face*>::iterator iterFTmp = iterF--;
//					spFNearby.erase(iterFTmp);
//					continue;
//				}
//				int iCLeft = getCellIndex(pF->getLeftCell());
//				int iCRight = getCellIndex(pF->getRightCell());
//				BdryPatch3D* pBPLeft = getBdryPatch(iCLeft);
//				BdryPatch3D* pBPRight = getBdryPatch(iCRight);
//				if (pBPLeft == pBPRight) {
//					std::set<Face*>::iterator iterFTmp = iterF--;
//					spFNearby.erase(iterFTmp);
//				}
//			}
//
//			// For each vert, find all pairs of edges that can implode.  For
//			// now, assume that all such pairs are going to spiral down as a
//			// single group.  This saves trying to seperate them into mutual
//			// encroaching groups.
//			const int iNumEdges = int(spFNearby.size());
//			bool *aqImplodable = new bool[iNumEdges];
//			bool qAnyImplodable = false;
//			double *adEdgeLengths = new double[iNumEdges];
//			Vert **apVNeigh = new Vert*[iNumEdges];
//
//			int iEdge = 0;
//			for (iterF = spFNearby.begin(); iterF != spFNearby.end();
//					++iterF, iEdge++) {
//				Face *pF = *iterF;
//				aqImplodable[iEdge] = false;
//				apVNeigh[iEdge] =
//						(pF->getVert(0) == pV) ?
//								pF->getVert(1) : pF->getVert(0);
//			}
//
//			// Outer loop goes all the way to iNumEdges-1 to set edge length for
//			// [iNumEdges-1].
//			for (int iOuter = 0; iOuter < iNumEdges; iOuter++) {
//				double adVecOut[] = adDIFF3D(apVNeigh[iOuter]->getCoords(),
//						pV->getCoords());
//				double dLenOut = dMAG3D(adVecOut);
//				adEdgeLengths[iOuter] = dLenOut;
//				for (int iInner = iOuter + 1; iInner < iNumEdges; iInner++) {
//					double adVecIn[] = adDIFF3D(apVNeigh[iInner]->getCoords(),
//							pV->getCoords());
//					double dLenIn = dMAG3D(adVecIn);
//					double dCosTheta = dDOT3D(adVecOut, adVecIn)
//							/ (dLenOut * dLenIn);
//					// If the angle between the edges exceeds 45 degrees, then there
//					// is no way that this pair can implode.
//					if (dCosTheta >= M_SQRT1_2) {
//						double dRatio = dLenOut / dLenIn;
//						// Now adjust dRatio to be in the range [0.5,1].  This works
//						// whether dRatio > 1 or < 1.
//						double dK = floor(-log(dRatio) / M_LN2);
//						double dRatio_Star = dRatio * pow(2., dK);
//						logMessage(4,
//								" Cos(theta) is %g.  Unscaled ratio is %g; scaled is %g.  This is",
//								dCosTheta, dRatio, dRatio_Star);
//						// Using iFuzzyComp here to avoid possible roundoff weirdness
//						// with lengths that differ by powers of two.
//						assert(
//								iFuzzyComp(dRatio_Star, 0.5) >= 0
//										&& iFuzzyComp(dRatio_Star, 1) <= 0);
//						if (iFuzzyComp(dRatio_Star, 0.5 / dCosTheta) >= 0
//								&& iFuzzyComp(dRatio_Star, dCosTheta) <= 0) {
//							// So we're in the evil range.  This pair of edges could
//							// implode.
//							logMessage(4, " evil\n");
//							aqImplodable[iOuter] = aqImplodable[iInner] = true;
//							qAnyImplodable = true;
//						} else {
//							logMessage(4, " okay\n");
//						}
//					}
//				} // Done with inner loop for edge-to-edge comparison
//			} // Done with outer loop for edge-to-edge comparison
//
//			if (qAnyImplodable) {
//				pV->markAsSmallAngleVert(true);
//				// Find the shortest and 2nd shortest edges in the bundle.  If the
//				// shortest and 2nd shortest have a length ratio > 1.5, split all
//				// edges at the shortest edge length from the common vert.  Otherwise,
//				// split at half the shortest edge length.  This will maximize the
//				// minimum edge length after insertion.
//				double dShortestLen = LARGE_DBL, dSecondLen = LARGE_DBL;
//				int iShortest = -1;
//				//	int iSecond = -1;
//				for (iEdge = 0; iEdge < iNumEdges; iEdge++) {
//					if (				// true ||
//					aqImplodable[iEdge]) {
//						double dThisLen = adEdgeLengths[iEdge];
//						if (dThisLen < dShortestLen) {
//							//	      iSecond = iShortest;
//							dSecondLen = dShortestLen;
//							iShortest = iEdge;
//							dShortestLen = dThisLen;
//						} else if (dThisLen < dSecondLen) {
//							//	      iSecond = iEdge;
//							dSecondLen = dThisLen;
//						}
//					}
//				}
//
//				assert(dShortestLen <= dSecondLen);
//				//	assert(iSecond != -1 && iShortest != -1);
//
//				// Now split these edges.
//				double dInsertionLen = dShortestLen / 2;
//				//	bool qSplitShortest = true;
//				if (dSecondLen / dShortestLen > 1.5) {
//					dInsertionLen = dShortestLen;
//					//	  qSplitShortest = false;
//					logMessage(3,
//							"Splitting at a distance of %f from the common vert.\n",
//							dInsertionLen);
//				}
//
//				for (iEdge = 0; iEdge < iNumEdges; iEdge++) {
//					if ((				// true ||
//					aqImplodable[iEdge]) && ((iEdge != iShortest)
//					//	       || (qSplitShortest)
//							)) {
//						// Split this edge.
//						Vert *pVNeigh = apVNeigh[iEdge];
//						double adVec[] = adDIFF3D(pVNeigh->getCoords(),
//								pV->getCoords());
//						NORMALIZE3D(adVec);
//						adVec[XDIR] = pV->x() + adVec[XDIR] * dInsertionLen;
//						adVec[YDIR] = pV->y() + adVec[YDIR] * dInsertionLen;
//						adVec[ZDIR] = pV->z() + adVec[ZDIR] * dInsertionLen;
//
//						// Since cell ID's are changing constantly with insertion, find a
//						// place to insert the new vertex the hard way.
//						Face *pFCand = pV->getHintFace();
//						assert(pFCand->doFullCheck() && pFCand->hasVert(pV));
//						Face *pFStart = pFCand;
//						Cell *pCCand = pFCand->getLeftCell();
//						while (!pFCand->hasVert(pVNeigh)) {
//							for (int i = 0; i < pCCand->getNumFaces(); i++) {
//								Face *pF = pCCand->getFace(i);
//								if (pF != pFCand && pF->hasVert(pV)) {
//									pFCand = pF;
//									break;
//								}
//							}
//							pCCand = pFCand->getOppositeCell(pCCand);
//							if (pFCand == pFStart)
//								break; // A way out in case something
//							// goes wrong and the loop goes infinite.
//						}
//						assert(pFCand->hasVert(pV));
//						assert(pFCand->hasVert(pVNeigh));
//						Vert *pVNew = createVert(adVec);
//						insertOnFace(pVNew, pFCand, pCCand, true);
//						qChanged = true;
//						iNumChanges++;
//						pVNew->setType(Vert::eBdryCurve);
//						logMessage(4,
//								"Inserted between verts at (%12g, %12g, %12g)",
//								pV->x(), pV->y(), pV->z());
//						logMessage(4, "and (%12g, %12g, %12g)\n", pVNeigh->x(),
//								pVNeigh->y(), pVNeigh->z());
//						logMessage(4,
//								"  New point located at (%12g, %12g, %12g)\n",
//								pVNew->x(), pVNew->y(), pVNew->z());
//					} // Done splitting an implodable edge
//				} // Done checking all edges for splitting
//			} // Done with if (any edges are implodable)
//			delete[] aqImplodable;
//			delete[] adEdgeLengths;
//			delete[] apVNeigh;
//		}
//		logMessage(2, "Done with this pass. Made %d changes.\n", iNumChanges);
//	} while (qChanged);
//}

static void
vBuildSortedEdgeList(const EntContainer<TriCellCV> & ECTriCV,
		     const GR_index_t iNEdge, FaceSort *&aFSEdges)
{
  if (aFSEdges != NULL)
    delete[] aFSEdges;
  aFSEdges = new FaceSort[iNEdge];
  GR_index_t iEdge = 0, iCell = 0, i;
  // Add edges of triangular cells

  for (i = 0; i < ECTriCV.lastEntry(); i++) {
    TriCellCV *pTCCV = ECTriCV.getEntry(i);
    if (pTCCV->isDeleted())
      continue;
    aFSEdges[iEdge++] = FaceSort(pTCCV->getVert(0), pTCCV->getVert(1), iCell);
    aFSEdges[iEdge++] = FaceSort(pTCCV->getVert(1), pTCCV->getVert(2), iCell);
    aFSEdges[iEdge++] = FaceSort(pTCCV->getVert(2), pTCCV->getVert(0), iCell);
    iCell++;
  }

  assert(iEdge == iNEdge);
  std::sort(aFSEdges, aFSEdges + iEdge);
}

//@@ This is a general conversion, for [mixed] 2D and surface meshes
// Meshes with only triangles will simply have no quads.
// Meshes with no boundary edges (surface meshes) have an empty ECBFCV.
void
SurfMesh::convertFromCellVert(EntContainer<TriCellCV> & ECTriCV)
{
  int iNTri = ECTriCV.lastEntry();
  int iNCell = iNTri;
  int iNEdge = 3 * iNTri;
  int iNQuad = 0;
  int iNBEdge = 0;

  // Set up enough blank entries in the connectivity tables to hold all
  // the data.  Don't set up edges yet, because for a non-manifold
  // surface mesh we can't tell in advance how many there will be.
  m_ECTri.reserve(iNTri);
  m_ECQuad.reserve(iNQuad);
  m_ECEdgeBF.reserve(iNBEdge);

  int i;
  // This code is needed if the SurfMesh already had data before.
  for (i = 0; i < iNTri; i++) {
    m_ECTri.getEntry(i)->resetAllData();
  }
  for (i = 0; i < iNQuad; i++) {
    m_ECQuad.getEntry(i)->resetAllData();
  }
  for (i = 0; i < iNBEdge; i++) {
    m_ECEdgeBF.getEntry(i)->resetAllData();
  }

  // If there are no quads in the mesh, tag it as simplicial so that
  // swapping, insertion, etc, will be globally legal on the mesh.
  m_simplexMesh = (iNQuad == 0);

  // This struct is used to hold edge descriptions for sorting.

  // Set up the individual edges and sort them
  FaceSort *aFSEdges = 0;

  vBuildSortedEdgeList(ECTriCV, iNEdge, aFSEdges);

  int *aiOrientErrors = new int[iNCell];
  int iNumErrors = 0, iPasses = 0;

  // Attempt to fix orientation errors.

  std::set<BdryPatch3D*> sBPFlip;
  int iNumUnmatched, iNumNonMan, iNumMisorient;
  do {
    iNumUnmatched = 0, iNumNonMan = 0, iNumMisorient = 0;
    iPasses++;
    for (i = 0; i < iNCell; aiOrientErrors[i++] = 0) {
    }
    for (i = 0; i < iNEdge;) {
      Vert *pV0 = aFSEdges[i].getVert(0);
      Vert *pV1 = aFSEdges[i].getVert(1);

      // How many faces are incident on this edge?
      int iNumInc = 1;
      while (i + iNumInc < iNEdge
	  && qIdenticalVerts(aFSEdges[i], aFSEdges[i + iNumInc]))
	iNumInc++;
      assert(i + iNumInc <= iNEdge);

      switch (iNumInc)
	{
	case 1:
	  // An unmatched edge.
	  logMessage(2, "Unmatched along edge from %u to %u\n", iVertIndex(pV0),
		     iVertIndex(pV1));
	  aiOrientErrors[aFSEdges[i].iLeftCell()] += 4;
	  iNumUnmatched++;
	  break;
	case 2:
	  if ((pV0 != aFSEdges[i + 1].getVert(1))
	      || (pV1 != aFSEdges[i + 1].getVert(0))) {
	    logMessage(2, "Orientation error along edge from %u to %u\n",
		       iVertIndex(pV0), iVertIndex(pV1));
	    iNumMisorient++;
	    aiOrientErrors[aFSEdges[i].iLeftCell()]++;
	    aiOrientErrors[aFSEdges[i + 1].iLeftCell()]++;
	  }
	  break;
	default:
	  // Could have more than three triangles even, so let's not
	  // restrict ourselves.
	  logMessage(2, "Non-manifold along edge from %u to %u\n",
		     iVertIndex(pV0), iVertIndex(pV1));
	  for (int ii = i; ii < i + iNumInc; ii++) {
	    aiOrientErrors[aFSEdges[ii].iLeftCell()] += 4;
	  }
	  iNumNonMan++;
	  break;
	}
      i += iNumInc;
    }
    // Non-manifold edges themselves aren't an error.
    iNumErrors = iNumUnmatched + iNumMisorient;

    // Now try to flip all the faces with more than one badly-oriented edge.
    int aiErrorFreq[] =
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    if (iNumErrors + iNumNonMan > 0) {
      logMessage(2, "Out of %u edges, %u had errors.\n", iNEdge, iNumErrors);
      for (i = 0; i < iNCell; i++) {
	aiErrorFreq[aiOrientErrors[i]]++;
      }
      logMessage(2, "Out of %d cells, %d had no orientation errors.\n", iNCell,
		 aiErrorFreq[0]);
      logMessage(2, "%d had one bad edge, %d had two, %d had three.\n",
		 aiErrorFreq[1], aiErrorFreq[2], aiErrorFreq[3]);
      logMessage(2, "%d had one unmatched / non-manifold (UNM) edges.\n",
		 aiErrorFreq[4]);
      logMessage(
	  2, "%d had one UNM, one bad orient; %d had one UNM, two bad orient\n",
	  aiErrorFreq[5], aiErrorFreq[6]);
      logMessage(2, "%d had two unmatched / non-manifold edges.\n",
		 aiErrorFreq[8]);
      logMessage(2, "%d had two UNM, one bad orient.\n", aiErrorFreq[9]);
      logMessage(2, "%d had three unmatched / non-manifold edges.\n",
		 aiErrorFreq[12]);
    }
    if (iNumErrors > 0) {
      // Start out with non-manifold/unmatched edges.
      if (aiErrorFreq[12] != 0) {
	logMessage(1, "Surface has triangles with three non-manifold edges; "
		   "deleting these.\n");
	for (i = 0; i < iNCell; i++) {
	  if (aiOrientErrors[i] == 12) {
	    TriCellCV* pTCCV = ECTriCV.getEntry(i);
	    logMessage(2, "Deleting cell %d (%u, %u, %u)\n", i,
		       iVertIndex(pTCCV->getVert(0)),
		       iVertIndex(pTCCV->getVert(1)),
		       iVertIndex(pTCCV->getVert(2)));
	    pTCCV->markAsDeleted();
	    iNTri--;
	  }
	}
	ECTriCV.purge();
	logMessage(2, "ECTri: %u, %u, %u\n", ECTriCV.lastEntry(),
		   ECTriCV.capacity(), ECTriCV.max_size());
	iNCell = iNTri;
	iNEdge = 3 * iNTri;
      }
      else if (aiErrorFreq[2] + aiErrorFreq[3] + aiErrorFreq[5] + aiErrorFreq[6]
	  > 0) {
	logMessage(1, "Surface mesh has inconsistent orientation!\n");
	logMessage(1, "You're in real trouble here, but I'll try to fix it.\n");
	logMessage(2, "Flipping all cells with 2 or 3 bad edges.\n");
	for (i = 0; i < iNCell; i++) {
	  if (aiOrientErrors[i] >= 2 && aiOrientErrors[i] % 4 != 0) {
	    ECTriCV.getEntry(i)->vSwitchOrientation();
	    sBPFlip.insert(getBdryPatch(i));
	  }
	}
      }
      logMessage(2, "Rebuilding edge list.\n");
      vBuildSortedEdgeList(ECTriCV, iNEdge, aFSEdges);
    }
    // Now flip all the patches that had a cell flipped.
    std::set<BdryPatch3D*>::iterator iterBP;
    for (iterBP = sBPFlip.begin(); iterBP != sBPFlip.end(); iterBP++) {
      (*iterBP)->vSwitchOrientation();
      // Flipping normals here fixes tet2 at the expense of breaking a
      // bunch of other cases.
    }
  }
  while (iNumErrors != 0 && iPasses < 25);

  // Read the data off
  int iBFace = 0, iCell;
  bool qOrientationOK = true;
  for (i = 0; i < iNEdge;) {
    Vert *pV0 = aFSEdges[i].getVert(0);
    Vert *pV1 = aFSEdges[i].getVert(1);

    // How many faces are incident on this edge?
    int iNumInc = 1;
    while (i + iNumInc < iNEdge
	&& qIdenticalVerts(aFSEdges[i], aFSEdges[i + iNumInc]))
      iNumInc++;
    assert(i + iNumInc <= iNEdge);

    switch (iNumInc)
      {
      case 0:
      case 1:
	// These cases are impossible.
	logMessage(0, "Warning: Orphan edge (no mate) along edge\n");
	switch (pV0->getSpaceDimen())
	  {
	  case 2:
	    logMessage(0, "         between vert %5u (%10g, %10g)\n",
		       iVertIndex(pV0), pV0->x(), pV0->y());
	    logMessage(0, "         and          %5u (%10g, %10g)\n",
		       iVertIndex(pV1), pV1->x(), pV1->y());
	    break;
	  case 3:
	    logMessage(0, "         between vert %5u (%10g, %10g, %10g)\n",
		       iVertIndex(pV0), pV0->x(), pV0->y(), pV0->z());
	    logMessage(0, "         and          %5u (%10g, %10g, %10g)\n",
		       iVertIndex(pV1), pV1->x(), pV1->y(), pV1->z());
	    break;
	  default:
	    assert(0);
	    break;
	  }
	break;
      case 2:
	// The standard case.  This always happens in 2D, and for surface
	// meshes only fails at edges where the mesh is locally non-manifold.
	{
	  if ((pV0 != aFSEdges[i + 1].getVert(1))
	      || (pV1 != aFSEdges[i + 1].getVert(0))) {
	    qOrientationOK = false;
	    iNumErrors++;
	    aiOrientErrors[aFSEdges[i].iLeftCell()]++;
	    aiOrientErrors[aFSEdges[i + 1].iLeftCell()]++;
	    logMessage(
		2,
		"Warning: Inconsistent orientation of input polygons along edge\n");
	    switch (pV0->getSpaceDimen())
	      {
	      case 2:
		logMessage(2, "         between vert %5u (%10g, %10g)\n",
			   iVertIndex(pV0), pV0->x(), pV0->y());
		logMessage(2, "         and          %5u (%10g, %10g)\n",
			   iVertIndex(pV1), pV1->x(), pV1->y());
		break;
	      case 3:
		logMessage(2, "         between vert %5u (%10g, %10g, %10g)\n",
			   iVertIndex(pV0), pV0->x(), pV0->y(), pV0->z());
		logMessage(2, "         and          %5u (%10g, %10g, %10g)\n",
			   iVertIndex(pV1), pV1->x(), pV1->y(), pV1->z());
		break;
	      default:
		assert(0);
	      }
	  }
	  if (qOrientationOK) {
	    // Only do this if orientation is good, to prevent potential
	    // problems with too many edges per cell or vice versa
	    int iC0 = aFSEdges[i].iLeftCell();
	    int iC1 = aFSEdges[i + 1].iLeftCell();
	    Cell *pC0, *pC1;
	    if (iC0 >= 0)
	      pC0 = getCell(iC0);
	    else {
	      BFace *pBF = getBFace(iBFace++);
	      pC0 = pBF;
	    }
	    if (iC1 >= 0)
	      pC1 = getCell(iC1);
	    else {
	      BFace *pBF = getBFace(iBFace++);
	      pC1 = pBF;
	    }

	    Face *pF = getNewFace();

	    pF->assign(pC0, pC1, pV0, pV1);
	    pC0->addFace(pF);
	    pC1->addFace(pF);
	  }
	  break;
	}
      default:
	// Handle non-manifold cases.
	logMessage(2, "Non-manifold surface along edge\n");
	logMessage(2, "         between vert %5u (%10g, %10g, %10g)\n",
		   iVertIndex(pV0), pV0->x(), pV0->y(), pV0->z());
	logMessage(2, "         and          %5u (%10g, %10g, %10g)\n",
		   iVertIndex(pV1), pV1->x(), pV1->y(), pV1->z());
	// Make an array with all the cells incident on this face.
	Cell **apCInc = new Cell*[iNumInc];
	MultiEdge *pME = getNewMultiEdge();
	for (int ii = 0; ii < iNumInc; ii++) {
	  iCell = aFSEdges[i + ii].iLeftCell();
	  if (iCell < 0)
	    vFatalError("Illegal boundary connectivity (hanging edge)\n",
			"SurfMesh::vConvertFromCellVert(...)");
	  apCInc[ii] = getCell(iCell);
	  apCInc[ii]->addFace(pME);
	}
	pME->vAssign(pV0, pV1, apCInc, iNumInc);
	delete[] apCInc;
	break;
      }
    i += iNumInc;
  }
  // Done with this, so nuke it.
  delete[] aFSEdges;

  // Check whether all faces have the same orientation that they're
  // supposed to from the CellCV's.  Those that don't will have this
  // problem because of MultiEdge's.
  for (iCell = 0; iCell < iNCell; iCell++) {
    Vert *pV0CV = ECTriCV[iCell].getVert(0);
    Vert *pV1CV = ECTriCV[iCell].getVert(1);
    Vert *pV2CV = ECTriCV[iCell].getVert(2);

    Cell *pC = getCell(iCell);
    pC->setRegion(ECTriCV[iCell].getRegion());
    Vert *pV0 = pC->getVert(0);
    Vert *pV1 = pC->getVert(1);
    Vert *pV2 = pC->getVert(2);

    if (!((pV0CV == pV0 && pV1CV == pV1 && pV2CV == pV2)
	|| (pV0CV == pV1 && pV1CV == pV2 && pV2CV == pV0)
	|| (pV0CV == pV2 && pV1CV == pV0 && pV2CV == pV1))) {
      Face *pF0 = pC->getFace(0);
      Face *pF1 = pC->getFace(1);
      Face *pF2 = pC->getFace(2);
      assert(pF0->getType() == Face::eMultiEdge);
      if (pF1->getType() == Face::eEdgeFace) {
	pC->assignFaces(pF1, pF2, pF0);
      }
      else if (pF2->getType() == Face::eEdgeFace) {
	pC->assignFaces(pF2, pF0, pF1);
      }
      else {
	logMessage(2, "%s %s\n",
		   "You have a face bounded by three multi-edges;",
		   "you're probably in trouble here.");
	logMessage(2, "  CellCV says verts are: %u %u %u\n", iVertIndex(pV0CV),
		   iVertIndex(pV1CV), iVertIndex(pV2CV));
	logMessage(2, "  SurfMesh says verts are: %u %u %u\n", iVertIndex(pV0),
		   iVertIndex(pV1), iVertIndex(pV2));
      }
    }
  }

#ifndef NDEBUG
  for (iCell = 0; iCell < iNCell; iCell++) {
    Vert *pV0CV = ECTriCV[iCell].getVert(0);
    Vert *pV1CV = ECTriCV[iCell].getVert(1);
    Vert *pV2CV = ECTriCV[iCell].getVert(2);

    Cell *pC = getCell(iCell);
    Vert *pV0 = pC->getVert(0);
    Vert *pV1 = pC->getVert(1);
    Vert *pV2 = pC->getVert(2);

    if (!((pV0CV == pV0 && pV1CV == pV1 && pV2CV == pV2)
	|| (pV0CV == pV1 && pV1CV == pV2 && pV2CV == pV0)
	|| (pV0CV == pV2 && pV1CV == pV0 && pV2CV == pV1))) {
      logMessage(2, "Possible problem with triangle orientation: \n");
      logMessage(2, "  CellCV says verts are: %d %d %d\n", iVertIndex(pV0CV),
		 iVertIndex(pV1CV), iVertIndex(pV2CV));
      logMessage(2, "  SurfMesh says verts are: %d %d %d\n", iVertIndex(pV0),
		 iVertIndex(pV1), iVertIndex(pV2));
    }
  }
#endif

  if (!qOrientationOK) {
    int aiErrorFreq[] =
      { 0, 0, 0, 0 };
    for (i = 0; i < iNCell; i++) {
      aiErrorFreq[aiOrientErrors[i]]++;
    }
    logMessage(2, "Out of %d cells, %d had no orientation errors.\n", iNCell,
	       aiErrorFreq[0]);
    logMessage(2, "%d had one bad edge, %d had two, %d had three.\n",
	       aiErrorFreq[1], aiErrorFreq[2], aiErrorFreq[3]);
    //      if (aiErrorFreq[2] + aiErrorFreq[3] > 0 &&
    //	aiErrorFreq[2] != iNCell) {
    //        // Re-orient anything with 2 or 3 bad edges.
    //        vMessage
    //	(1, "Trying heuristics to fix surface triangle orientation...\n");
    //        vMessage(2, "Flipping all cells with 2 or 3 bad edges.\n");
    //        int iNTriCV = ECTriCV.lastEntry();
    //        for (i = 0; i < iNCell; i++) {
    //	if (aiOrientErrors[i] >= 2) {
    //	  if (i < iNTriCV) {
    //	    ECTriCV.getEntry(i)->vSwitchOrientation();
    //	  }
    //	  else {
    //	    ECQuadCV.getEntry(i-iNTriCV)->vSwitchOrientation();
    //	  }
    //	}
    //        }
    //        vMessage(2, "Trying conversion again.\n");
    //        vConvertFromCellVert(ECTriCV, ECQuadCV, ECBFCV);
    //      }
    //      else {
    vFatalError("Polygons do not all have the same orientation",
		"SurfMesh::vConvertFromCellVert(...)");
    //      }
  }
  delete[] aiOrientErrors;
  assert(iBFace == iNBEdge);
}

void
SurfMesh::purgeVerts(map<Vert*, Vert*>* vert_map)
// Renumber vertices, removing those that have been marked as deleted
{
  GR_index_t *aiVertInd = new GR_index_t[iNumVerts()];
  GR_index_t iV, iVActive;

  setAllHintFaces();
  setVertFaceNeighbors();

  iVActive = 0;
  // Copy the vertex data down within the list and establish a lookup so
  // that the faces can figure out where their verts are now stored.
  for (iV = 0; iV < iNumVerts(); iV++) {
    Vert *pV = getVert(iV);
    if (!pV->isDeleted()) {
      aiVertInd[iV] = iVActive;
      Vert *pVA = getVert(iVActive);
      if (vert_map)
	vert_map->insert(std::make_pair(pV, pVA));
      // Copy the data from its old to its new location
      (*pVA) = (*pV);
      iVActive++;
    }
    else
      aiVertInd[iV] = INT_MAX;
  }

  clearVertFaceNeighbors();
  clearAllHintFaces();
  for (int iF = getNumFaces() - 1; iF >= 0; iF--) {
    Face *pF = getFace(iF);
    if (!pF->doFullCheck())
      continue;
    Vert *pVNew0 = getVert(aiVertInd[iVertIndex(pF->getVert(0))]);
    Vert *pVNew1 = getVert(aiVertInd[iVertIndex(pF->getVert(1))]);

    (dynamic_cast<EdgeFace*>(pF))->setVerts(pVNew0, pVNew1);
  }

  for (int iME = getNumMultiEdges() - 1; iME >= 0; iME--) {
    MultiEdge *pME = getMultiEdge(iME);
    if (!pME->doFullCheck())
      continue;

    Vert *pVNew0 = getVert(aiVertInd[iVertIndex(pME->getVert(0))]);
    Vert *pVNew1 = getVert(aiVertInd[iVertIndex(pME->getVert(1))]);

    (dynamic_cast<EdgeFace*>(pME))->setVerts(pVNew0, pVNew1);
  }

  delete[] aiVertInd;

  setupVerts(iVActive);
  setAllHintFaces();
  setVertFaceNeighbors();
}
void
SurfMesh::purgeAllEntities(std::map<Vert*, Vert*>* vert_map,
			   std::map<Face*, Face*>* face_map,
			   std::map<Cell*, Cell*>* cell_map,
			   std::map<BFace*, BFace*>* bface_map)
{
  // need to reassign cell data
  if (m_numSurfaces > 0) {
    std::map<Cell*, RefFace*> cellToSurf;

    for (GR_index_t i = 0; i < getNumCells(); i++)
      cellToSurf.insert(std::make_pair(getCell(i), m_surfaces[i]));
    std::map<Cell*, Cell*> cellMap;
    Mesh2D::purgeAllEntities(vert_map, face_map, &cellMap, bface_map);
    for (std::map<Cell*, Cell*>::iterator it = cellMap.begin();
	it != cellMap.end(); ++it) {
      Cell * oldCell = (*it).first;
      Cell * newCell = (*it).second;
      RefFace * oldSurface = cellToSurf.find(oldCell)->second;
      m_surfaces[getCellIndex(newCell)] = oldSurface;
    }
    if (cell_map) {
      cell_map->swap(cellMap);

    }
  }
  else
    Mesh2D::purgeAllEntities(vert_map, face_map, cell_map, bface_map);

}
SurfMesh::SurfMesh(const char strBaseFileName[], const int iQualMeas) :
    Mesh2D(iQualMeas), m_ECME(), m_Bdry3D(NULL), m_bdryPatches(NULL), m_numBdryPatches(
	0)
{
  m_isLengthScaleFromCellSizes = true;

  //@@ Read data from a user-defined format input file and convert to Mesh2D
  GR_index_t iNVerts, iNFaces, iNCells, iNIntBFaces;
  bool qFaceVert, qCellVert, qFaceCell, qCellFace, qCellRegion;
  bool qIntBFaceFace, qIntBFaceVert, qIntBFaceBC;
  GR_index_t (*a2iFaceVert)[2];
  GR_sindex_t (*a2iFaceCell)[2];
  GR_index_t (*a2iCellVert)[3], (*a2iCellFace)[3];
  int *aiCellRegion;
  int *aiIntBFaceBC;
  GR_index_t (*a2iIntBFaceFace)[2], (*a2iIntBFaceVert)[2];

  readFile_Surface(strBaseFileName, iNVerts, iNFaces, iNCells, iNIntBFaces,
		   qFaceVert, qCellVert, qFaceCell, qCellFace, qCellRegion,
		   qIntBFaceFace, qIntBFaceVert, qIntBFaceBC, m_ECVerts,
		   a2iFaceVert, a2iFaceCell, a2iCellVert, a2iCellFace,
		   aiCellRegion, a2iIntBFaceFace, aiIntBFaceBC,
		   a2iIntBFaceVert);

  if (m_ECVerts.lastEntry() == 0)
    vFatalError("No vertex coordinate data read", "surface mesh input");

  bool qConnectOK = false;
  if (qFaceVert && qFaceCell)
    qConnectOK = true;	// face-based structure

  else if (qCellVert)
    qConnectOK = true;	// FE structure

  if (!qConnectOK)
    vFatalError("Inadequate connectivity data read; consult user's manual.",
		"surface mesh input");

  //@@@ Convert face-based data to Mesh2D
  if (qFaceVert && qFaceCell) {
    // If the number of cells, faces, and boundary faces have been
    // provided, we can check for mesh validity (and the presence of
    // quads in the mesh; quads are not yet allowed).
    if (iNFaces == 0)
      vFatalError("Mesh size unknown.", "surface mesh input");
    //@@@@ Verify mesh size and that numbering starts from 0
    {
      GR_index_t iTmpNCells = 0;
      GR_index_t iTmpNBFaces = 0, iTmpNVerts = 0;
      for (GR_index_t iF = 0; iF < iNFaces; iF++) {
	GR_sindex_t iC0 = a2iFaceCell[iF][0];
	GR_sindex_t iC1 = a2iFaceCell[iF][1];
	GR_sindex_t maxTmp = std::max(iC0, iC1);
	assert(maxTmp > 0);
	iTmpNCells = std::max(iTmpNCells, GR_index_t(maxTmp));
	if (iC0 < 0 || iC1 < 0)
	  iTmpNBFaces++;
	GR_index_t iV0 = a2iFaceVert[iF][0];
	GR_index_t iV1 = a2iFaceVert[iF][1];
	iTmpNVerts = std::max((iTmpNVerts), std::max((iV0), (iV1)));
      }
      if (iNCells != 0 && iNCells != iTmpNCells + 1)
	vFatalError("Number of cells in mesh doesn't match number given.",
		    "surface mesh input");
      if ((iNVerts != 0 && iNVerts != iTmpNVerts + 1)
	  || iTmpNVerts + 1 != m_ECVerts.lastEntry())
	vFatalError("Number of verts in mesh doesn't match number given.",
		    "surface mesh input");
      iNCells = iTmpNCells + 1;
      iNVerts = iTmpNVerts + 1;
    }

    if (3 * iNCells != 2 * iNFaces)
      vFatalError("Mesh is invalid or contains quadrilaterals.",
		  "surface mesh input");

    setupFaces(iNFaces);	// Ensure the existence of enough
    setupTriCells(iNCells);	// blank faces, bdry faces, and cells.

    //@@@ Assign cells to regions
    // The default region, set when the cell was created, is 1.
    if (qCellRegion) {
      assert(aiCellRegion != NULL);
      for (GR_index_t iCell = 0; iCell < getNumCells(); iCell++)
	getCell(iCell)->setRegion(aiCellRegion[iCell]);
      // Done with aiCellRegion, so get rid of it.  Not allocated unless
      // the cell region was actually read.
      delete[] aiCellRegion;
    }

    for (GR_index_t i = 0; i < iNFaces; i++) {
      int iV0 = a2iFaceVert[i][0];
      int iV1 = a2iFaceVert[i][1];

      // Add entry to cell-face table or bdry face list on the fly
      int iC0 = a2iFaceCell[i][0];
      int iC1 = a2iFaceCell[i][1];

      Cell *pC0 = getCell(iC0);
      Cell *pC1 = getCell(iC1);

      // This needs to be done whether pC? is a real cell or a bdry face
      Face *pF = getFace(i);
      pC0->addFace(pF);
      pC1->addFace(pF);

      // Tell the face everything it needs to know.
      pF->assign(pC0, pC1, getVert(iV0), getVert(iV1));
    } // Done assigning for this face.
    setAllHintFaces();
    setVertFaceNeighbors();

    if (!(isValid()))
      vFoundBug("surface mesh input from user-format face->cell-style file");
  }

  else if (qCellVert) {
    if (iNCells == 0)
      vFatalError("Mesh size unknown.", "surface mesh input");
    //@@@@ Verify that number of verts is correct
    {
      GR_index_t iTmpNVerts = 0;
      for (GR_index_t iC = 0; iC < iNCells; iC++) {
	GR_index_t iV0 = a2iCellVert[iC][0];
	GR_index_t iV1 = a2iCellVert[iC][1];
	GR_index_t iV2 = a2iCellVert[iC][2];
	iTmpNVerts = max(max(iTmpNVerts, iV0), max(iV1, iV2));
      }
      if (iNVerts != 0 && iNVerts != iTmpNVerts + 1)
	vFatalError("Number of verts in mesh doesn't match number given.",
		    "surface mesh input");
    }

    // Set up data for all triangles
    for (GR_index_t iC = 0; iC < iNCells; iC++) {
      Vert *pV0 = getVert(a2iCellVert[iC][0]);
      Vert *pV1 = getVert(a2iCellVert[iC][1]);
      Vert *pV2 = getVert(a2iCellVert[iC][2]);
      createTriCell(pV0, pV1, pV2);
    }
    setAllHintFaces();
    setVertFaceNeighbors();

    if (!(isValid()))
      vFoundBug("surface mesh input from user-format FE-style file");
  } // Done translating cell-vert/bface-vert/bface-BC data

  //@@@ Tag faces and verts
  for (GR_index_t iVert = 0; iVert < iNumVerts(); iVert++)
    getVert(iVert)->setType(Vert::eInterior);
  for (GR_index_t iFace = 0; iFace < getNumFaces(); iFace++) {
    Face *pF = getFace(iFace);
    Cell *pCL = pF->getLeftCell();
    Cell *pCR = pF->getRightCell();
    if (pCL->getType() == Cell::eBdryEdge || pCR->getType() == Cell::eBdryEdge)
      continue; // Do boundaries separately.
    // Mark internal boundaries

    if (pCL->getRegion() != pCR->getRegion()) {
      pF->setFaceLoc(Face::eBdryTwoSide);
      pF->getVert(0)->setType(Vert::eBdryTwoSide);
      pF->getVert(1)->setType(Vert::eBdryTwoSide);
    }
    else
      // Not an internal boundary
      pF->setFaceLoc(Face::eInterior);
  }

  // Tag ordinary boundary verts, and tag verts that lie at the
  // intersection of internal and regular boundary faces as boundary
  // apexes, so that they can never be smoothed or removed.
  for (GR_index_t iBFace = 0; iBFace < getNumBdryFaces(); iBFace++) {
    BFace *pBF = getBFace(iBFace);
    pBF->getFace()->setFaceLoc(Face::eBdryFace);
    for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
      Vert *pV = pBF->getVert(iV);
      switch (pV->getVertType())
	{
	case Vert::eBdryTwoSide:
	  pV->setType(Vert::eBdryApex);
	  break;
	case Vert::eInterior:
	  pV->setType(Vert::eBdry);
	  break;
	default:
	  // Nothing needs to be done
	  break;
	}
    }
  }

  m_isLengthScaleFromCellSizes = true;

  if (qFaceVert)
    delete[] a2iFaceVert;
  if (qCellVert)
    delete[] a2iCellVert;
  if (qFaceCell)
    delete[] a2iFaceCell;
  if (qCellFace)
    delete[] a2iCellFace;
}

//void
//SurfMesh::createPatches()
//{
//  if (m_Bdry3D)
//    delete m_Bdry3D;
//  m_Bdry3D = new Bdry3D(m_ECVerts);
//
//  for (GR_index_t iC = 0; iC < getNumCells(); iC++) {
//    Vert *apV[3];
//
//    Cell *pC = getCell(iC);
//
//    apV[0] = pC->getVert(0);
//    apV[1] = pC->getVert(2);
//    apV[2] = pC->getVert(1);
//
//    int iBCL = iDefaultBC;
//    int iBCR = iInvalidBC;
//
//    int iRegL = iOutsideRegion;
//    int iRegR = iDefaultRegion;
//
//    BdryPolygon *pBPoly = new BdryPolygon(iBCL, iBCR, iRegL, iRegR, 3, apV);
//    m_Bdry3D->vAddPatch(pBPoly);
//    setBdryPatch(iC, pBPoly);
//  } // Done with loop over all BFaces to create Bdry3D object.
//  assert(m_Bdry3D->iNumBdryPatches() == this->getNumCells());
//}

GR_index_t
SurfMesh::getFaceIndex(const Face* const pF) const
{
  assert(pF->getType() == Face::eEdgeFace || pF->getType() == Face::eMultiEdge);
  if (pF->getType() == Face::eEdgeFace)
    return m_ECEdgeF.getIndex(dynamic_cast<const EdgeFace*>(pF));
  else {
    return m_ECME.getIndex(dynamic_cast<const MultiEdge*>(pF))
	+ m_ECME.lastEntry();
  }
}

Face*
SurfMesh::getFace(const GR_index_t i) const
{
  assert(i < getNumFaces());
  if (i < getNumEdgeFaces())
    return m_ECEdgeF.getEntry(i);
  else
    return m_ECME.getEntry(i - getNumEdgeFaces());
}

