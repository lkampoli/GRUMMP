#include <map>

using std::multimap;
using std::pair;

#include <stdio.h>
#include <stdlib.h>

#include "GR_VolMesh.h"

void
VolMesh::readFromFile(const char strBaseFileName[], int* numVars, double** data)
{
  m_isLengthScaleFromCellSizes = true;
  m_strictPatchChecking = false;

  logMessage(1, "Finding a reader for file %s\n", strBaseFileName);
  bool qIsUGridFile = strstr(strBaseFileName, ".ugrid");
  bool qIsVGridFile = strstr(strBaseFileName, ".cogsg");
  bool qIsVTKFile = strstr(strBaseFileName, ".vtk");
  bool qIsVMeshFile = strstr(strBaseFileName, ".vmesh")
      || strstr(strBaseFileName, ".grvmesh");
  ;
  bool qIsMEditFile = strstr(strBaseFileName, ".mesh");
  bool qIsTecPlotFile = strstr(strBaseFileName, ".tec");

  int* aiBCList = nullptr;
  if (qIsUGridFile) {
    assert(!numVars && !data);
    readUGridBinary(strBaseFileName, this, aiBCList);
  }
  else if (qIsVGridFile) {
    assert(!numVars && !data);
    readVGridBinary(strBaseFileName, this);
  }
  else if (qIsVTKFile) {
    assert(!numVars && !data);
    readVTK_ASCII(strBaseFileName, this);
  }
  else if (qIsMEditFile) {
    assert(!numVars && !data);
    readMEdit3D(strBaseFileName, this);
  }
  else if (qIsVMeshFile) {
    assert(!numVars && !data);
    readVMeshFile(strBaseFileName, this, aiBCList);
  }
  else if (qIsTecPlotFile) {
    assert(!numVars && !data);
    readTecPlotASCII(strBaseFileName, this);
  }
  else {
    readTetGenFile(strBaseFileName, this, aiBCList, numVars, data);
  }

  if (!aiBCList) {
    aiBCList = new int[getNumTotalBdryFaces()];
    GR_index_t iBF = 0;
    for (; iBF < getNumBdryFaces(); iBF++) {
      aiBCList[iBF] = iDefaultBC;
    }
    for (; iBF < getNumTotalBdryFaces(); iBF++) {
      aiBCList[iBF] = iInvalidBC;
    }
  }

  if (m_Bdry3D) {
    // Now infer which bdry faces are on which patches.  Start by
    // identifying all bdry verts, and finding which patches, if any,
    // they fall on.  While it's not very efficient to do this by
    // checking all bdry verts against all patches, this approach has
    // the virtue of working, and for meshes that aren't absolutely
    // huge, it should be no problem from the runtime point of view.
    multimap<Vert*, BdryPatch3D*> mmVertPatch;
    for (int iV = getNumVerts() - 1; iV >= 0; iV--) {
      Vert *pV = getVert(iV);
      for (int iP = m_Bdry3D->iNumPatches() - 1; iP >= 0; iP--) {
	BdryPatch3D *pBP3D = (*m_Bdry3D)[iP];
	if (!pBP3D->qDisjointFrom(pV)) {
	  logMessage(4, "Vert at (%10.6f, %10.6f, %10.6f) lies on patch %10p\n",
		     pV->x(), pV->y(), pV->z(), pBP3D);
	  mmVertPatch.insert(pair<Vert*, BdryPatch3D*>(pV, pBP3D));
	}
      }
    }

    // Now for every bdry face, identify a patch common to its verts.
    for (int iBF = getNumBdryFaces() - 1; iBF >= 0; iBF--) {
      BFace *pBF = getBFace(iBF);
      Vert *pV0 = pBF->getVert(0);
      Vert *pV1 = pBF->getVert(1);
      Vert *pV2 = pBF->getVert(2);
      pV0->setType(Vert::eBdry);
      pV1->setType(Vert::eBdry);
      pV2->setType(Vert::eBdry);

      // How many patches for each one?
      pair<multimap<Vert*, BdryPatch3D*>::iterator,
	  multimap<Vert*, BdryPatch3D*>::iterator> pRange0, pRange1, pRange2;
      pRange0 = mmVertPatch.equal_range(pV0);
      pRange1 = mmVertPatch.equal_range(pV1);
      pRange2 = mmVertPatch.equal_range(pV2);
      multimap<Vert*, BdryPatch3D*>::iterator iter0, iter1, iter2;
      BdryPatch3D* pBP3D = NULL;
      for (iter0 = pRange0.first; pBP3D == NULL && iter0 != pRange0.second;
	  iter0++) {
	for (iter1 = pRange1.first; pBP3D == NULL && iter1 != pRange1.second;
	    iter1++) {
	  if (iter0->second == iter1->second) {
	    for (iter2 = pRange2.first;
		pBP3D == NULL && iter2 != pRange2.second; iter2++) {
	      if (iter0->second == iter2->second) {
		pBP3D = iter0->second;
	      } // Inner if
	    } // Loop for vert 3
	  } // 1 and 2 have a face in common
	} // Loop for vert 2
      } // Loop for vert 1
      assert(pBP3D != NULL);
      logMessage(4, "BFace %3d lies on patch %10p\n", iBF, pBP3D);
      logMessage(4, "  coords: (%10.6f, %10.6f, %10.6f),\n", pV0->x(), pV0->y(),
		 pV0->z());
      logMessage(4, "  coords: (%10.6f, %10.6f, %10.6f),\n", pV1->x(), pV1->y(),
		 pV1->z());
      logMessage(4, "  coords: (%10.6f, %10.6f, %10.6f),\n", pV2->x(), pV2->y(),
		 pV2->z());
      pBF->setPatch(pBP3D);
    }
    delete[] aiBCList;
    aiBCList = NULL;

    // Definitively identify bdry apexes and bdry curve pts: any vert
    // with faceson
    // more than two patches is an apex; those on exactly two are bdry
    // curve pts.  This assumes specifically that there aren't any
    // internal bdrys.
    for (int iV = getNumVerts() - 1; iV >= 0; iV--) {
      Vert *pV = getVert(iV);
      if (pV->getVertType() == Vert::eBdry) {
	int iCount = mmVertPatch.count(pV);
	if (iCount == 2) {
	  // FIX ME:
	  // Could also be a bizarre apex, where the curve separating
	  // two patches has a kink in it.
	  pV->setType(Vert::eBdryCurve);
	}
	else if (iCount > 2) {
	  // FIX ME:
	  // Could be a curve with one or more internal faces.
	  pV->setType(Vert::eBdryApex);
	}
	logMessage(4, "Patch count: %d.  Vert type: %d.\n", iCount,
		   pV->getVertType());
      }
    }
  }
  else {
    logMessage(2, "Creating bdry patches directly from mesh data.\n");
    logMessage(2,
	       "You should specify a bdry file if you have one available.\n");
    createPatchesFromBdryFaces(aiBCList);
    initializeSubsegs();
    if (getNumQuadBdryFaces() == 0) {
      // Can't handle quads on the bdry.
      bool canChangeSurf = areSurfaceEdgeChangesAllowed();
      double temp = getMaxAngleForSurfaceEdgeChanges();
      allowSurfaceEdgeChanges(45);
      identifyVertexTypesGeometrically();
      if (canChangeSurf)
	allowSurfaceEdgeChanges(temp);
      else
	disallowSurfaceEdgeChanges();
    }
  }
  if (aiBCList)
    delete[] aiBCList;

  setVertFaceNeighbors();
  setAllHintFaces();
  // Check to make sure the correct number of bdry faces exist, etc.
  if (isSimplicial()) {
    assert(
	4 * getNumCells() + getNumBdryFaces() + 2 * getNumIntBdryFaces()
	    == 2 * getNumFaces());
    assert(getNumQuadFaces() == 0);
    assert(getNumQuadBdryFaces() == 0);
    assert(getNumPyrCells() == 0);
    assert(getNumPrismCells() == 0);
    assert(getNumHexCells() == 0);
  }
  if (!(isMeshOrientationOK(true)))
    vFatalError("Invalid or improperly-oriented mesh", "3D mesh input");
}

static void
addVert(GR_index_t vertList[4], const GR_index_t vert)
{
  for (int ii = 0; ii < 4; ii++) {
    if (vertList[ii] == vert)
      return;
    if (vertList[ii] == INT_MAX) {
      vertList[ii] = vert;
      return;
    }
  }
  assert(0); // Should never get here.
}

void
readVMeshFile(const char strBaseFileName[], VolMesh* const pVM, int*& aiBCList)
{
  // Now read all of the mesh data from file
  GR_index_t iNVerts, iNFaces, iNCells, iNBFaces, iNIntBFaces;
  bool qFaceVert, qFaceCell, qCellRegion;
  bool qBFaceFace, qBFaceVert, qBFaceBC;
  bool qIntBFaceFace, qIntBFaceVert, qIntBFaceBC;
  GR_index_t (*a2iFaceVert)[3] = nullptr;
  GR_sindex_t (*a2iFaceCell)[2] = nullptr;
  int* aiCellRegion = nullptr;
  GR_index_t *aiBFaceFace = nullptr, (*a2iBFaceVert)[3] = nullptr;
  GR_index_t (*a2iIntBFaceFace)[2] = nullptr, (*a2iIntBFaceVert)[3] = nullptr;
  int *aiBFaceBC = nullptr, *aiIntBFaceBC = nullptr;
  readFile_Native(strBaseFileName, iNVerts, iNFaces, iNCells, iNBFaces,
		  iNIntBFaces, qFaceVert, qFaceCell, qCellRegion, qBFaceFace,
		  qBFaceVert, qBFaceBC, qIntBFaceFace, qIntBFaceVert,
		  qIntBFaceBC, *(pVM->pECGetVerts()), a2iFaceVert, a2iFaceCell,
		  aiCellRegion, aiBFaceFace, aiBFaceBC, a2iBFaceVert,
		  a2iIntBFaceFace, aiIntBFaceBC, a2iIntBFaceVert);
  if (pVM->getNumVerts() == 0 || iNVerts == 0)
    vFatalError("No vertex coordinate data read", "volume mesh input");

  bool qConnectOK = false;
  if (qFaceVert && qFaceCell)
    qConnectOK = true; // face-based structure

  if (!qConnectOK)
    vFatalError("Inadequate connectivity data read; consult user's manual.",
		"volume mesh input");

  //@@@ Convert face-based data to VolMesh
  // This is the easiest conversion case
  if (qFaceVert && qFaceCell) {
    // If the number of cells, faces, and boundary faces have been
    // provided, we can check for mesh validity.
    if (iNFaces == 0)
      vFatalError("Mesh size unknown.", "volume mesh input");

    //@@@@ Verify mesh size
    {
      GR_sindex_t iTmpNCells = 0;
      GR_index_t iTmpNBFaces = 0, iTmpNVerts = 0;
      for (GR_index_t iF = 0; iF < iNFaces; iF++) {
	GR_sindex_t iC0 = a2iFaceCell[iF][0];
	GR_sindex_t iC1 = a2iFaceCell[iF][1];
	iTmpNCells = MAX3(iTmpNCells, iC0, iC1);
	if (iC0 < 0 || iC1 < 0)
	  iTmpNBFaces++;

	GR_index_t iV0 = a2iFaceVert[iF][0];
	GR_index_t iV1 = a2iFaceVert[iF][1];
	GR_index_t iV2 = a2iFaceVert[iF][2];
	iTmpNVerts = max(max(iTmpNVerts, iV0), max(iV1, iV2));
      }
      if (iNCells == 0)
	iNCells = iTmpNCells + 1;
      else if (iNCells != GR_index_t(iTmpNCells + 1))
	vFatalError("Number of cells in mesh doesn't match number given.",
		    "volume mesh input");

      if (iNBFaces == 0)
	iNBFaces = iTmpNBFaces;
      else if (iNBFaces != iTmpNBFaces)
	vFatalError("Number of bdry faces in mesh doesn't match number given.",
		    "volume mesh input");

      if (iNVerts == 0)
	iNVerts = iTmpNVerts + 1;
      else if (iNVerts != iTmpNVerts + 1)
	vFatalError("Number of verts in mesh doesn't match number given.",
		    "volume mesh input");
    }

    // For now, at least, we're going to continue reading only all-tet meshes
    // this way, so all the cells are tets.

    GR_index_t cellToVert[iNCells][4];
    for (GR_index_t iC = 0; iC < iNCells; iC++) {
      cellToVert[iC][0] = cellToVert[iC][1] = cellToVert[iC][2] =
	  cellToVert[iC][3] = INT_MAX;
    }
    for (GR_index_t iF = 0; iF < iNFaces; iF++) {
      // Add verts to both cells
      GR_sindex_t leftCell = a2iFaceCell[iF][0];
      GR_sindex_t rightCell = a2iFaceCell[iF][1];
      for (int ii = 0; ii < 3; ii++) {
	GR_index_t vert = a2iFaceVert[iF][ii];
	if (leftCell >= 0)
	  addVert(cellToVert[leftCell], vert);
	if (rightCell >= 0)
	  addVert(cellToVert[rightCell], vert);
      }
    }

    // Now create cells; faces are created automatically.
    for (GR_index_t iC = 0; iC < iNCells; iC++) {
      // Create each, one at a time
      bool qExisted = false;
      int region = qCellRegion ? aiCellRegion[iC] : iDefaultRegion;
      pVM->createTetCell(qExisted, pVM->getVert(cellToVert[iC][0]),
			 pVM->getVert(cellToVert[iC][1]),
			 pVM->getVert(cellToVert[iC][2]),
			 pVM->getVert(cellToVert[iC][3]), region);
      assert(!qExisted);
    }

    aiBCList = new int[iNBFaces + iNIntBFaces];
    // Now create bdry faces
    for (GR_index_t iBF = 0; iBF < iNBFaces; iBF++) {
      Vert *pV0 = pVM->getVert(a2iBFaceVert[iBF][0]);
      Vert *pV1 = pVM->getVert(a2iBFaceVert[iBF][1]);
      Vert *pV2 = pVM->getVert(a2iBFaceVert[iBF][2]);
      pV0->setType(Vert::eBdry);
      pV1->setType(Vert::eBdry);
      pV2->setType(Vert::eBdry);
      // These can't be connected both sides.
      Face *pF = findCommonFace(pV0, pV1, pV2, nullptr, false);
      assert(pF);
      // Can't specify bdry condition directly
      BFace* pBF = pVM->createBFace(pF);
      assert(pBF && pBF->getType() == Cell::eTriBFace);
      aiBCList[iBF] = aiBFaceBC[iBF];
    }

    for (GR_index_t iIntBF = 0; iIntBF < iNIntBFaces; iIntBF++) {
      Vert *pV0 = pVM->getVert(a2iIntBFaceVert[iIntBF][0]);
      Vert *pV1 = pVM->getVert(a2iIntBFaceVert[iIntBF][1]);
      Vert *pV2 = pVM->getVert(a2iIntBFaceVert[iIntBF][2]);
      if (pV0->getVertType() == Vert::eBdry)
	pV0->setType(Vert::eBdryCurve);
      else
	pV0->setType(Vert::eBdryCurve);
      if (pV1->getVertType() == Vert::eBdry)
	pV1->setType(Vert::eBdryCurve);
      else
	pV1->setType(Vert::eBdryCurve);
      if (pV2->getVertType() == Vert::eBdry)
	pV2->setType(Vert::eBdryCurve);
      else
	pV2->setType(Vert::eBdryCurve);

      // These must be connected both sides.
      Face *pF = findCommonFace(pV0, pV1, pV2, nullptr, true);
      assert(pF);
      // Can't specify bdry condition directly
      BFace* pBF = pVM->createBFace(pF);
      assert(pBF && pBF->getType() == Cell::eIntTriBFace);
      aiBCList[iNBFaces + iIntBF] = aiIntBFaceBC[iIntBF];
    }
//		// FIX ME
//		// To do Euler's formula, we would need to identify all the edges,
//		// which isn't worth doing yet.
//		pVM->reserveTriFaces(iNFaces + iNIntBFaces); // Ensure the existence of enough
//		pVM->reserveTets(iNCells); // blank faces, bdry faces, and cells.  IntBFaces aren't reserved.
//		pVM->reserveTriBFaces(iNBFaces);
//		//			m_ECIntTriBF.reserve(iNIntBFaces);
//		aiBCList = new int[iNBFaces];
//		std::multimap<int, int> CellIndexMap;
//		std::multimap<int, int>::iterator mapit;
//		GR_index_t iBFace = 0;
//		GR_index_t iIntBFace = 0;
//		GR_index_t iFaceCounter = 0;
//		for (GR_index_t i = 0; i < iNFaces; i++) {
//			// Add entry to cell-face table or bdry face list on the fly
//			GR_sindex_t iC0 = a2iFaceCell[i][0];
//			GR_sindex_t iC1 = a2iFaceCell[i][1];
//			// Tell the face everything it needs to know.
//			GR_index_t iV0 = a2iFaceVert[i][0];
//			GR_index_t iV1 = a2iFaceVert[i][1];
//			GR_index_t iV2 = a2iFaceVert[i][2];
//			bool qIntBdry = false;
//			Cell *pC0, *pC1;
//			if (findCommonFace(pVM->getVert(iV0), pVM->getVert(iV1),
//					pVM->getVert(iV2), NULL, true))
//				qIntBdry = true;
//
//			if (!qIntBdry) {
//				Face* pF = pVM->getFace(iFaceCounter);
//				iFaceCounter++;
//				if (iC0 < 0) {
//					BFace* pBF = pVM->getBFace(iBFace);
//					pBF->assignFaces(pF);
//					pF->setFaceLoc(Face::eBdryFace);
//					aiBCList[iBFace] = abs(iC0);
//					pC0 = pBF;
//					iBFace++;
//				} else {
//					pC0 = pVM->getCell(iC0);
//					pC0->addFace(pF);
//				}
//				if (iC1 < 0) {
//					BFace* pBF = pVM->getBFace(iBFace);
//					pBF->assignFaces(pF);
//					pF->setFaceLoc(Face::eBdryFace);
//					aiBCList[iBFace] = abs(iC1);
//					pC1 = pBF;
//					iBFace++;
//				} else {
//					pC1 = pVM->getCell(iC1);
//					pC1->addFace(pF);
//				}
//				pF->assign(pC0, pC1, pVM->getVert(iV0), pVM->getVert(iV1),
//						pVM->getVert(iV2));
//			} else {
//				iIntBFace++;
//			}
//		}
//		assert(iBFace == iNBFaces);
//		assert(iIntBFace == iNIntBFaces);
//		// lets handle IntBFaces;
//		//@@@ Assign BC's and check boundary consistency
//		// If the BC's were read on a boundary-face by boundary-face
//		// basis, then these BC's silently override anything read
//		// elsewhere. Also, verify that the boundary face and its
//		// corresponding face have the same vertices, if data exists to do
//		// this check.
//		//@@@@ The easy case:  real face ID read
//		if (qBFaceFace && (qBFaceBC || qBFaceVert)) {
//			for (GR_index_t iBF = 0; iBF < iNBFaces; iBF++) {
//				GR_index_t iFace = aiBFaceFace[iBF];
//				Face* pF = pVM->getFace(iFace);
//				//@@@@@ Assign BC
//				if (qBFaceBC) {
//					int iBC = aiBFaceBC[iBF];
//					Cell* pC = pF->getLeftCell();
//					if (pC->getType() == Cell::eTet) {
//						pC = pF->getRightCell();
//					}
//					assert(pC->getType() == Cell::eTriBFace);
//					aiBCList[iBF] = abs(iBC);
//				} // Done assigning BC
//				  //@@@@@ Verify vertex consistency
//				if (qBFaceVert) {
//					GR_index_t iV0 = pVM->getVertIndex(pF->getVert(0));
//					GR_index_t iV1 = pVM->getVertIndex(pF->getVert(1));
//					GR_index_t iV2 = pVM->getVertIndex(pF->getVert(2));
//					if ((iV0 != a2iBFaceVert[iBF][0]
//							&& iV0 != a2iBFaceVert[iBF][1]
//							&& iV0 != a2iBFaceVert[iBF][2])
//							|| (iV1 != a2iBFaceVert[iBF][0]
//									&& iV1 != a2iBFaceVert[iBF][1]
//									&& iV1 != a2iBFaceVert[iBF][2])
//							|| (iV2 != a2iBFaceVert[iBF][0]
//									&& iV2 != a2iBFaceVert[iBF][1]
//									&& iV2 != a2iBFaceVert[iBF][2])) {
//						vFatalError(
//								"Verts of bdry face don't match verts of actual face",
//								"volume mesh input");
//					}
//				} // Done consistency checking for vertices on boundary
//			} // Done looping over all bdry faces
//		} // Done with cases where face ID was read for each bdry faceelse
//		  //@@@@ Assign BC's if bface-vert and bdry cond data have been read
//		if (qBFaceVert && qBFaceBC) {
//			for (GR_index_t iBF = 0; iBF < iNBFaces+iNIntBFaces; iBF++) {
////	An elaborate consistency check that we can skip.  If there isn't such a face,
////  we'll find out when constructing the BFace anyway.
////				GR_index_t iV0 = a2iBFaceVert[iBF][0];
////				GR_index_t iV1 = a2iBFaceVert[iBF][1];
////				GR_index_t iV2 = a2iBFaceVert[iBF][2];
////				Vert* pV0 = pVM->getVert(iV0);
////				Vert* pV1 = pVM->getVert(iV1);
////				Vert* pV2 = pVM->getVert(iV2);
////				// The bface must already exist, but which one is it?
////				// Find all incident faces for pV0.
////				std::set<Cell*> spCIncident;
////				std::set<Vert*> spVNeigh;
////				std::set<BFace*> spBFIncident;
////				BFace* pBF = pBFInvalidBFace; /* Init for picky compilers. */
////				findNeighborhoodInfo(pV0, spCIncident, spVNeigh, &spBFIncident);
////				assert(spVNeigh.count(pV1) == 1);
////				assert(spVNeigh.count(pV2) == 1);
////				assert(spBFIncident.size() > 0);
////				for (std::set<BFace*>::iterator iter = spBFIncident.begin();
////						iter != spBFIncident.end(); ++iter) {
////					pBF = *iter;
////					assert(pBF->hasVert(pV0));
////					if (pBF->hasVert(pV1) && pBF->hasVert(pV2))
////						break;
////				}
////				assert(pBF->hasVert(pV1) && pBF->hasVert(pV2));
//				aiBCList[iBF] = abs(aiBFaceBC[iBF]);
//			} // Done looping over all bdry faces to assign BC's
//		} // Done assigning BC's based on bdry verts aloneelse
//		if (qBFaceFace || qBFaceVert || qBFaceBC) {
//			vWarning("Insufficient boundary data read to be useful; ignoring.");
//		}

  } // Done converting face-based data

  // FIX ME: Add code to set up bdry data after 3D mesh input conversion.
  //@@@ Assign cells to regions
  // The default region, set when the cell was created, is 1.
  if (qCellRegion) {
    for (GR_index_t iCell = 0; iCell < pVM->getNumCells(); iCell++)
      pVM->getCell(iCell)->setRegion(aiCellRegion[iCell]);
    // Done with aiCellRegion, so get rid of it.
    delete[] aiCellRegion;
  }
  //@@@ Tag faces and verts
  for (GR_index_t iVert = 0; iVert < pVM->getNumVerts(); iVert++)
    pVM->getVert(iVert)->setType(Vert::eInterior);
//	for (GR_index_t iFace = 0; iFace < pVM->getNumFaces(); iFace++) {
//		Face* pF = pVM->getFace(iFace);
//		Cell* pCL = pF->getLeftCell();
//		Cell* pCR = pF->getRightCell();
//		if (pF->isBdryFace())
//			continue; // Do boundaries separately.
//
//		// Mark internal boundaries
//		if (pCL->getRegion() != pCR->getRegion()) {
//			pF->setFaceLoc(Face::eBdryTwoSide);
//			//				pF->getVert(0)->setType(Vert::eBdryTwoSide);
//			//				pF->getVert(1)->setType(Vert::eBdryTwoSide);
//			//				pF->getVert(2)->setType(Vert::eBdryTwoSide);
//			pVM->createBFace(pF);
//		} else
//			// Not an internal boundary
//			pF->setFaceLoc(Face::eInterior);
//	}

//	// Tag ordinary boundary verts, and tag verts that lie at the
//	// intersection of internal and regular boundary faces as boundary
//	// apexes, so that they can never be smoothed or removed.
//	// At this point, geometric tagging of BdryApex and BdryCurve
//	// vertices isn't done.
//	for (GR_index_t iBFace = 0; iBFace < pVM->getNumBdryFaces(); iBFace++) {
//		BFace* pBF = pVM->getBFace(iBFace);
//		for (int iV = 0; iV < pBF->getNumVerts(); iV++) {
//			Vert* pV = pBF->getVert(iV);
//			switch (pV->getVertType()) {
//			case Vert::eBdryTwoSide:
//				pV->setType(Vert::eBdryApex);
//				break;
//			case Vert::eInterior:
//				pV->setType(Vert::eBdry);
//				break;
//			default:
//				// Nothing needs to be done
//				break;
//			}
//		}
//	}
  if (qFaceVert)
    delete[] a2iFaceVert;

  if (qFaceCell)
    delete[] a2iFaceCell;

  if (qBFaceFace)
    delete[] aiBFaceFace;

  if (qBFaceVert)
    delete[] a2iBFaceVert;

  if (qBFaceBC)
    delete[] aiBFaceBC;

  if (qIntBFaceFace)
    delete[] a2iIntBFaceFace;

  if (qIntBFaceVert)
    delete[] a2iIntBFaceVert;

  if (qIntBFaceBC)
    delete[] aiIntBFaceBC;
}

void
VolMesh::convertFromCellVert(const GR_index_t iNCells,
			     const GR_index_t iNBdryFaces,
			     const GR_index_t a2iCellVert[][4],
			     const GR_index_t aiRegions[],
			     const GR_index_t a2iBFaceVert[][3])
{
  m_simplexMesh = true;
  allowEdgeSwapping();
  allowSurfaceEdgeChanges(5.0);

  m_ECTriF.clear();
  m_ECQuadF.clear();
  m_ECTriBF.clear();
  m_ECQuadBF.clear();
  m_ECIntTriBF.clear();
  m_ECIntQuadBF.clear();
  m_ECTet.clear();
  m_ECPyr.clear();
  m_ECPrism.clear();
  m_ECHex.clear();

  // Create a bunch of tets.
  for (GR_index_t iC = 0; iC < iNCells; iC++) {
    bool qExist;
    (void) createTetCell(qExist, getVert(a2iCellVert[iC][0]),
			 getVert(a2iCellVert[iC][1]),
			 getVert(a2iCellVert[iC][2]),
			 getVert(a2iCellVert[iC][3]), aiRegions[iC]);
  }

  if (iNBdryFaces > 0) {
    // Now create bdry faces.
    for (GR_index_t iBF = 0; iBF < iNBdryFaces; iBF++) {
      Vert *pV0 = getVert(a2iBFaceVert[iBF][0]);
      Vert *pV1 = getVert(a2iBFaceVert[iBF][1]);
      Vert *pV2 = getVert(a2iBFaceVert[iBF][2]);
      Face *pF = findCommonFace(pV0, pV1, pV2, NULL, true);
      // This will create a regular or an internal bdry face as needed, though
      // there's no tie to geometry at this point.
      // Also, don't do this if this is already a bdry face.  This idiot-proofs
      // us against duplicates in the list of bdry faces.
      if (!pF->isBdryFace()) {
	(void) createBFace(pF);
      }
    }
  }
  else {
    // Identify places where bdry faces are needed.  Then create them, though
    // without any explicit bdry condition information.
    for (GR_index_t iF = 0; iF < getNumFaces(); iF++) {
      Face *pF = getFace(iF);
      assert(pF->isValid() && !pF->isDeleted());
      if (!pF->getLeftCell()->isValid() || !pF->getRightCell()->isValid()
	  || (pF->getLeftCell()->getRegion() != pF->getRightCell()->getRegion())) {
	// This handles both the interior and regular bdry face cases.
	(void) createBFace(pF);
      }
    }
  }

  // This will be true if the correct number of bdry faces have been provided.
  assert(
      4 * getNumTetCells() + getNumBdryFaces() + 2 * getNumIntBdryFaces()
	  == getNumFaces() * 2);

  assert(isValid());
  setAllHintFaces();
  setVertFaceNeighbors();
}

// This is the more general version, that will handle all
// GRUMMP-supported 3D cell types, and both tri and quad bdry faces.
void
VolMesh::convertFromCellVert(const GR_index_t nTet, const GR_index_t nPyr,
			     const GR_index_t nPrism, const GR_index_t nHex,
			     const GR_index_t nBdryTri,
			     const GR_index_t nBdryQuad,
			     const GR_index_t tetToVert[][4],
			     const GR_index_t pyrToVert[][5],
			     const GR_index_t prismToVert[][6],
			     const GR_index_t hexToVert[][8],
			     const GR_index_t bdryTriToVert[][3],
			     const GR_index_t bdryQuadToVert[][4])
{
  m_simplexMesh = (nPyr + nPrism + nHex + nBdryQuad == 0);
  allowEdgeSwapping();
  allowSurfaceEdgeChanges();

  assert(nBdryTri + nBdryQuad > 0);

  // Now create list of all faces and sort to pair the two sides of the face
  GR_index_t iNTriFaces = (4 * nTet + 4 * nPyr + 2 * nPrism + nBdryTri) / 2;
  GR_index_t iNQuadFaces = (nPyr + 3 * nPrism + 6 * nHex + nBdryQuad) / 2;

  m_ECTriF.clear();
  m_ECQuadF.clear();
  m_ECTriBF.clear();
  m_ECQuadBF.clear();
  m_ECIntTriBF.clear();
  m_ECIntQuadBF.clear();
  m_ECTet.clear();
  m_ECPyr.clear();
  m_ECPrism.clear();
  m_ECHex.clear();

  // Create tets.
  for (GR_index_t iC = 0; iC < nTet; iC++) {
    bool qExist;
    (void) createTetCell(qExist, getVert(tetToVert[iC][0]),
			 getVert(tetToVert[iC][1]), getVert(tetToVert[iC][2]),
			 getVert(tetToVert[iC][3]));
  }

  // Create pyrs.
  for (GR_index_t iC = 0; iC < nPyr; iC++) {
    bool qExist;
    Vert* pV0 = getVert(pyrToVert[iC][0]);
    Vert* pV1 = getVert(pyrToVert[iC][1]);
    Vert* pV2 = getVert(pyrToVert[iC][4]);
    Vert* pV3 = getVert(pyrToVert[iC][3]);
    Vert* pV4 = getVert(pyrToVert[iC][2]);
    (void) createPyrCell(qExist, pV0, pV1, pV2, pV3, pV4);
  }

  // Create prisms.
  for (GR_index_t iC = 0; iC < nPrism; iC++) {
    bool qExist;
    Vert* pV0 = getVert(prismToVert[iC][0]);
    Vert* pV1 = getVert(prismToVert[iC][1]);
    Vert* pV2 = getVert(prismToVert[iC][4]);
    Vert* pV3 = getVert(prismToVert[iC][3]);
    Vert* pV4 = getVert(prismToVert[iC][2]);
    Vert* pV5 = getVert(prismToVert[iC][5]);
    (void) createPrismCell(qExist, pV0, pV1, pV2, pV3, pV4, pV5);
  }

  // Create hexs.
  for (GR_index_t iC = 0; iC < nHex; iC++) {
    bool qExist;
    Vert* pV0 = getVert(hexToVert[iC][0]);
    Vert* pV1 = getVert(hexToVert[iC][1]);
    Vert* pV2 = getVert(hexToVert[iC][4]);
    Vert* pV3 = getVert(hexToVert[iC][3]);
    Vert* pV4 = getVert(hexToVert[iC][2]);
    Vert* pV5 = getVert(hexToVert[iC][5]);
    Vert* pV6 = getVert(hexToVert[iC][6]);
    Vert* pV7 = getVert(hexToVert[iC][7]);
    (void) createHexCell(qExist, pV0, pV1, pV2, pV3, pV4, pV5, pV6, pV7);
  }

  // Create bdry tris
  for (GR_index_t iBF = 0; iBF < nBdryTri; iBF++) {
    Vert *pV0 = getVert(bdryTriToVert[iBF][0]);
    Vert *pV1 = getVert(bdryTriToVert[iBF][1]);
    Vert *pV2 = getVert(bdryTriToVert[iBF][2]);
    Face *pF = findCommonFace(pV0, pV1, pV2);
    assert(pF->isValid());
    createBFace(pF);
  }

  // Create bdry quads
  for (GR_index_t iBF = 0; iBF < nBdryQuad; iBF++) {
    Vert *pV0 = getVert(bdryQuadToVert[iBF][0]);
    Vert *pV1 = getVert(bdryQuadToVert[iBF][1]);
    Vert *pV2 = getVert(bdryQuadToVert[iBF][2]);
    Vert *pV3 = getVert(bdryQuadToVert[iBF][3]);
    Face *pF = findCommonFace(pV0, pV1, pV2, pV3);
    assert(pF->isValid());
    createBFace(pF);
  }

  assert(getNumTriFaces() == iNTriFaces);
  assert(getNumQuadFaces() == iNQuadFaces);
  assert(isValid());
  setAllHintFaces();
  setVertFaceNeighbors();
}

void
readFile_Native(const char * const strBaseFileName, GR_index_t& iNumVerts,
		GR_index_t& iNumFaces, GR_index_t& iNumCells,
		GR_index_t& iNumBdryFaces, GR_index_t& iNumIntBdryFaces,
		bool& qFaceVert, bool& qFaceCell, bool& qCellRegion,
		bool& qBFaceFace, bool& qBFaceVert, bool& qBFaceBC,
		bool& qIntBFaceFace, bool& qIntBFaceVert, bool& qIntBFaceBC,
		EntContainer<Vert>& ECVerts, GR_index_t (*&a2iFaceVert)[3],
		GR_sindex_t (*&a2iFaceCell)[2], int *&aiCellRegion,
		GR_index_t *&aiBFaceFace, int *&aiBFaceBC,
		GR_index_t (*&a2iBFaceVert)[3],
		GR_index_t (*&a2iIntBFaceFace)[2], int *&aiIntBFaceBC,
		GR_index_t (*&a2iIntBFaceVert)[3])
{
  iNumVerts = 0;
  iNumFaces = 0;
  iNumCells = 0;
  iNumIntBdryFaces = 0;
  iNumBdryFaces = 0;
  qFaceVert = false;
  qFaceCell = false;
  qCellRegion = false;
  a2iFaceCell = reinterpret_cast<GR_sindex_t (*)[2]>(NULL);
  a2iFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);
  aiCellRegion = reinterpret_cast<int *>(NULL);
  qBFaceVert = false;
  qBFaceFace = false;
  qBFaceBC = false;
  aiBFaceFace = reinterpret_cast<GR_index_t *>(NULL);
  aiBFaceBC = reinterpret_cast<int *>(NULL);
  a2iBFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);
  qIntBFaceFace = false;
  qIntBFaceVert = false;
  qIntBFaceBC = false;
  a2iIntBFaceFace = reinterpret_cast<GR_index_t (*)[2]>(NULL);
  aiIntBFaceBC = reinterpret_cast<int *>(NULL);
  a2iIntBFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);

  int iVertOffset = 0, iFaceOffset = 0;
  int iCellOffset = 0, iBFaceOffset = 1;
  int res = -1;
  FILE *pFInFile = fopen("/dev/zero", "r");
  char strFileName[1024];
  char strTemp[1024];
  snprintf(strTemp, 1024, "vmesh");
  const char *pcWhere = strstr(strBaseFileName, strTemp);
  if (pcWhere) {
    snprintf(strFileName, 1024, "%s", strBaseFileName);
  }
  else {
    snprintf(strFileName, 1024, "%s.vmesh", strBaseFileName);
  }

  /*newfile vmesh*/
//  printf("Opening input file %s\n", strFileName);
  if (pFInFile)
    fclose(pFInFile);
  pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile)
    vFatalError("Couldn't open input file for reading", "3d mesh input");

  /*ncells nfaces nbfaces nibfaces nverts*/
  {
    char *line = NULL;
    size_t lineLen;
    getline(&line, &lineLen, pFInFile);
    res = sscanf(line, "%u %u %u %u %u\n", &iNumCells, &iNumFaces,
		 &iNumBdryFaces, &iNumIntBdryFaces, &iNumVerts);
    assert(res >= 0);
    if (res == 4) {
      // This is the normal case, in which no internal bdry faces are present.
      iNumVerts = iNumIntBdryFaces;
      iNumIntBdryFaces = 0;
    }
    free(line);
  }
  if (iNumVerts == 0)
    vFatalError(
	"Tried to read vertex data without specifying number of vertices",
	"3d mesh input");
  ECVerts.reserve(iNumVerts);
  {
    size_t bufSize = 1024;
    char* buffer = reinterpret_cast<char*>(calloc(bufSize, 1));
  for (GR_index_t iV = 0; iV < iNumVerts; iV++) {
    GR_index_t iRealVert = iV;
      double dXX, dYY, dZZ;
      getline(&buffer, &bufSize, pFInFile);
    /*verts: coords*/
      int vertType = Vert::eUnknown;
      res = sscanf(buffer, "%lf%lf%lf%d", &dXX, &dYY, &dZZ, &vertType);
      assert(res == 3 || res == 4);

      double adCoord[] =
	{ dXX, dYY, dZZ };
      Vert *pV = ECVerts.getEntry(iRealVert);
      pV->setCoords(3, adCoord);
      if (res == 4) {
	pV->setType(vertType);
      }
    }
    free(buffer);
  }


  if (iNumFaces == 0)
    vFatalError("Tried to read face data without specifying number of faces",
		"3d mesh input");
  /* Allocate both the face-vert and face-cell arrays if they haven't been.
   One or both may have to be deleted later. */
  if (!a2iFaceVert)
    a2iFaceVert = new GR_index_t[iNumFaces][3];
  if (!a2iFaceCell)
    a2iFaceCell = new GR_sindex_t[iNumFaces][2];
  if (iNumVerts == 0)
    vWarning(
	"Number of verts unknown while reading faces; sanity checking may be crippled.");
  if (iNumCells == 0)
    vWarning(
	"Number of cells unknown while reading faces; sanity checking may be crippled.");
  for (GR_index_t iF = 0; iF < iNumFaces; iF++) {
    GR_index_t iRealFace = iF;
    GR_sindex_t iCellA = iNumCells, iCellB = iNumCells;
    GR_index_t iVertA = 2 * iNumVerts, iVertB = 2 * iNumVerts;
    GR_index_t iVertC = 2 * iNumVerts;
    /*faces: cells verts*/
    res = fscanf(pFInFile, "%u%u %u%u%u\n", &iCellA, &iCellB, &iVertA, &iVertB,
		 &iVertC);
    assert(res >= 0); /* Remove an offset from the vertex indices if needed. */
    iVertA -= iVertOffset;
    iVertB -= iVertOffset;
    iVertC -= iVertOffset;
    if (((iNumVerts > 0)
	&& (iVertA >= iNumVerts || iVertC >= iNumVerts || iVertB >= iNumVerts)))
      vFatalError("Vertex index out of range", "3d mesh input");
    (a2iFaceVert)[iRealFace][0] = iVertA;
    (a2iFaceVert)[iRealFace][1] = iVertB;
    (a2iFaceVert)[iRealFace][2] = iVertC;

    /* Remove an offset from the vertex indices if needed. */
    if (iCellA >= 0)
      iCellA -= iCellOffset;
    if (iCellB >= 0)
      iCellB -= iCellOffset;
    /* Can not check cell index for out-of-range low, */
    /* because BC's have a negative index. */
    if (iNumCells > 0
	&& ((iCellA > 0 && GR_index_t(iCellA) >= iNumCells)
	    || (iCellB > 0 && GR_index_t(iCellB) >= iNumCells)))
      vFatalError("Cell index out of range", "3d mesh input");
    (a2iFaceCell)[iRealFace][0] = iCellA;
    (a2iFaceCell)[iRealFace][1] = iCellB;
  }

  qFaceVert = true;
  qFaceCell = true;

  if (iNumBdryFaces == 0)
    vFatalError(
	"Tried to read bdry face data without specifying number of bdry faces",
	"3d mesh input");
  /* Allocate the bface-face, bface-vert, and bface-BC arrays. */
  /* One or more may have to be deleted later. */
  if (!aiBFaceFace)
    aiBFaceFace = new GR_index_t[iNumBdryFaces];
  if (!aiBFaceBC)
    aiBFaceBC = new int[iNumBdryFaces];
  if (!a2iBFaceVert)
    a2iBFaceVert = new GR_index_t[iNumBdryFaces][3];
  qBFaceBC = true;
  qBFaceFace = true;
  qBFaceVert = true;
  if (iNumVerts == 0)
    vWarning(
	"# of verts unknown while reading bdry faces; sanity checking may be crippled.");
  if (iNumFaces == 0)
    vWarning(
	"# of faces unknown while reading bdry faces; sanity checking may be crippled.");
  for (GR_index_t iBF = 0; iBF < iNumBdryFaces; iBF++) {
    GR_index_t iRealBFace = iBF;
    GR_index_t iFace = iNumFaces, iBC = 0;
    GR_index_t iVertA = 2 * iNumVerts, iVertB = 2 * iNumVerts;
    GR_index_t iVertC = 2 * iNumVerts;
    /*bdryfaces: face bc verts*/
    res = fscanf(pFInFile, "%u %u %u%u%u\n", &iFace, &iBC, &iVertA, &iVertB,
		 &iVertC);
    assert(res >= 0); /* Remove offset from the vertex indices */
    iVertA -= iVertOffset;
    iVertB -= iVertOffset;
    iVertC -= iVertOffset;
    if (((iNumVerts > 0)
	&& (iVertA >= iNumVerts || iVertC >= iNumVerts || iVertB >= iNumVerts)))
      vFatalError("Vertex index out of range", "3d mesh input");
    (a2iBFaceVert)[iRealBFace][0] = iVertA;
    (a2iBFaceVert)[iRealBFace][1] = iVertB;
    (a2iBFaceVert)[iRealBFace][2] = iVertC;

    /* Remove offset from face index */
    iFace -= iFaceOffset;
    if (((iNumFaces > 0) && (iFace >= iNumFaces)))
      vFatalError("Face index out of range", "3d mesh input");
    (aiBFaceFace)[iRealBFace] = iFace;
    if (iBC <= 0)
      vFatalError("Boundary condition out of range", "3d mesh input");
    (aiBFaceBC)[iRealBFace] = iBC;
  }

  if (iNumIntBdryFaces != 0) {
    // Internal Boundaries
    if (!a2iIntBFaceFace)
      a2iIntBFaceFace = new GR_index_t[iNumIntBdryFaces][2];
    if (!a2iIntBFaceVert)
      a2iIntBFaceVert = new GR_index_t[iNumIntBdryFaces][3];
    if (!aiIntBFaceBC)
      aiIntBFaceBC = new int[iNumIntBdryFaces];
    for (GR_index_t iBF = 0; iBF < iNumIntBdryFaces; iBF++) {
      GR_index_t iRealIntBFace = iBF;
      GR_index_t iFace1 = iNumFaces, iFace2 = iNumFaces;
      GR_index_t iVertA = 2 * iNumVerts, iVertB = 2 * iNumVerts;
      GR_index_t iVertC = 2 * iNumVerts;
      /*bdryfaces: face bc verts*/
      res = fscanf(pFInFile, "%u %u %u%u%u\n", &iFace1, &iFace2, &iVertA,
		   &iVertB, &iVertC);
      assert(res >= 0); /* Remove offset from the vertex indices */
      iVertA -= iVertOffset;
      iVertB -= iVertOffset;
      iVertC -= iVertOffset;
      if (((iNumVerts > 0)
	  && (iVertA >= iNumVerts || iVertC >= iNumVerts || iVertB >= iNumVerts)))
	vFatalError("Vertex index out of range", "3d mesh input");
      (a2iIntBFaceVert)[iRealIntBFace][0] = iVertA;
      (a2iIntBFaceVert)[iRealIntBFace][1] = iVertB;
      (a2iIntBFaceVert)[iRealIntBFace][2] = iVertC;

      /* Remove offset from face index */
      iFace1 -= iFaceOffset;
      iFace2 -= iFaceOffset;
      if (((iNumFaces > 0) && (iFace1 >= iNumFaces || iFace2 >= iNumFaces)))
	vFatalError("Face index out of range", "3d mesh input");
      (a2iIntBFaceFace)[iRealIntBFace][0] = iFace1;
      (a2iIntBFaceFace)[iRealIntBFace][1] = iFace2;
      // No BC data is read, so all we can do is set it as invalid.
      aiIntBFaceBC[iRealIntBFace] = iInvalidBC;
    }

    qIntBFaceVert = true;
    qIntBFaceFace = true;
    qIntBFaceBC = true;
  }

  if (!aiCellRegion)
    aiCellRegion = new int[iNumCells];
  GR_index_t iC = 0;
  for (iC = 0; iC < iNumCells; iC++) {
    int iReg = 0;
    res = fscanf(pFInFile, "%u\n", &iReg);
    if (res == EOF)
      break;
    aiCellRegion[iC] = iReg;
    assert(iReg > 0);
  }
  // Fallback: if there's not data for all cells, put the rest in the default
  // region and hope for the best.
  for (; iC < iNumCells; iC++) {
    aiCellRegion[iC] = iDefaultRegion;
  }
  qCellRegion = true;
  /* Get rid of the face-vert array if it isn't used. */
  if (!(qFaceVert) && a2iFaceVert) {
    delete[] a2iFaceVert;
    a2iFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);
  }

  /* Get rid of the face-cell array if it isn't used. */
  if (!(qFaceCell) && (a2iFaceCell)) {
    delete[] a2iFaceCell;
    a2iFaceCell = reinterpret_cast<GR_sindex_t (*)[2]>(NULL);
  }

  /* Get rid of the cell-region array if it isn't used. */
  if (!(qCellRegion) && (aiCellRegion)) {
    delete[] aiCellRegion;
    aiCellRegion = reinterpret_cast<int *>(NULL);
  }

  /* Get rid of the bface-face array if it isn't used. */
  if (!(qBFaceFace) && (aiBFaceFace)) {
    delete[] aiBFaceFace;
    aiBFaceFace = reinterpret_cast<GR_index_t *>(NULL);
  }

  /* Get rid of the bface-BC array if it isn't used. */
  if (!(qBFaceBC) && (aiBFaceBC)) {
    delete[] aiBFaceBC;
    aiBFaceBC = reinterpret_cast<int *>(NULL);
  }

  /* Get rid of the bface-vert array if it isn't used. */
  if (!(qBFaceVert) && (a2iBFaceVert)) {
    delete[] a2iBFaceVert;
    a2iBFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);
  }

  /* Get rid of the int bface-face array if it isn't used. */
  if (!(qIntBFaceFace) && (a2iIntBFaceFace)) {
    delete[] a2iIntBFaceFace;
    a2iIntBFaceFace = reinterpret_cast<GR_index_t (*)[2]>(NULL);
  }

  /* Get rid of the int bface-vert array if it isn't used. */
  if (!(qIntBFaceVert) && (a2iIntBFaceVert)) {
    delete[] a2iIntBFaceVert;
    a2iIntBFaceVert = reinterpret_cast<GR_index_t (*)[3]>(NULL);
  }

  /* The following lines prevent compiler complaints. */
  iBFaceOffset++;
  iFaceOffset++;
  iCellOffset++;
  iVertOffset++;
  fclose(pFInFile);
}
