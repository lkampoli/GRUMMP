/*
 * WriteTriangle.cxx
 *
 *  Created on: 2014-07-09
 *      Author: zaide
 */
#include <stdio.h>
#include "GR_Mesh2D.h"

void
writeEdge(FILE * pFEdgeFile, FILE * pFPolyFile, const Mesh2D& OutMesh,
	  const Vert * pivotVert, const Vert * pivotVertOld,
	  const GR_index_t pivotVertNum, const Face *pF, GR_index_t iBC,
	  GR_index_t & iFaceCount, GR_index_t & iBFCount)
{

  if (pF->getFaceLoc() == Face::eBdryTwoSide)
    iBC = 10;
  else if (pF->getFaceLoc() == Face::eInterior)
    iBC = 0;
  else if (pivotVert && pivotVertOld && pF->getFaceLoc() == Face::eBdryFace)
    pivotVert = NULL;

  // need to change a vert.
  const Vert * v0 = pF->getVert(0);
  const Vert * v1 = pF->getVert(1);
  GR_index_t iV0 = OutMesh.getVertIndex(v0);
  GR_index_t iV1 = OutMesh.getVertIndex(v1);
  if (v0 == pivotVert)
    iV0 = pivotVertNum;
  else if (v0 == pivotVertOld)
    iV0 = pivotVertNum - 1;

  if (v1 == pivotVert)
    iV1 = pivotVertNum;
  else if (v1 == pivotVertOld)
    iV1 = pivotVertNum - 1;
  //	if(pivotVert && pivotVertOld){
  //		printf("%d %d\n", iV0+1,iV1+1);
  //		const_cast<Vert*>(v0)->printVertInfo();
  //		const_cast<Vert*>(v1)->printVertInfo();
  //		const_cast<Vert*>(pivotVert)->printVertInfo();
  //		const_cast<Vert*>(pivotVertOld)->printVertInfo();
  //	}
  fprintf(pFEdgeFile, "%u %u %u %u\n", ++iFaceCount, iV0 + 1, iV1 + 1, iBC);
  if (iBC > 0)
    fprintf(pFPolyFile, "%u %u %u %u\n", ++iBFCount, iV0 + 1, iV1 + 1, iBC);

}
void
writeEle(FILE * pFEleFile, const Mesh2D& OutMesh, const Vert * pivotVert,
	 const Vert * pivotVertOld, GR_index_t pivotVertNum, const Cell * pC,
	 GR_index_t & iCellCount)
{

  // need to change a vert.
  const Vert * v0 = pC->getVert(0);
  const Vert * v1 = pC->getVert(1);
  const Vert * v2 = pC->getVert(2);
  GR_index_t iV0 = OutMesh.getVertIndex(v0);
  GR_index_t iV1 = OutMesh.getVertIndex(v1);
  GR_index_t iV2 = OutMesh.getVertIndex(v2);
  if (v0 == pivotVert)
    iV0 = pivotVertNum;
  else if (v1 == pivotVert)
    iV1 = pivotVertNum;
  else if (v2 == pivotVert)
    iV2 = pivotVertNum;
  if (v0 == pivotVertOld)
    iV0 = pivotVertNum - 1;
  else if (v1 == pivotVertOld)
    iV1 = pivotVertNum - 1;
  else if (v2 == pivotVertOld)
    iV2 = pivotVertNum - 1;
  fprintf(pFEleFile, "%u %u %u %u\n", ++iCellCount, iV0 + 1, iV1 + 1, iV2 + 1);

}
void
writeTriangle(Mesh2D& OutMesh, const char strBaseFileName[])
{
  if (OutMesh.getNumQuadCells() != 0)
    return;

  // Triangle files all start with 1-indexing
  char strPolyFileName[FILE_NAME_LEN], strNodeFileName[FILE_NAME_LEN],
      strEleFileName[FILE_NAME_LEN], strEdgeFileName[FILE_NAME_LEN];

  sprintf(strPolyFileName, "%s.poly", strBaseFileName);
  sprintf(strNodeFileName, "%s.node", strBaseFileName);
  sprintf(strEleFileName, "%s.ele", strBaseFileName);
  sprintf(strEdgeFileName, "%s.edge", strBaseFileName);

  FILE *pFNodeFile = fopen(strNodeFileName, "w");
  FILE *pFPolyFile = fopen(strPolyFileName, "w");
  FILE *pFEleFile = fopen(strEleFileName, "w");
  FILE *pFEdgeFile = fopen(strEdgeFileName, "w");

  // first, collect all apex verts connected to an intbdry face, and store a hint
  GR_index_t numNewVerts = 0;
  // checked verts contains verts' we've looked at + the number of extra times they should appear
  std::map<Vert*, int> checkedVerts;
  std::map<Vert*, BFace*> apexVerts;
  for (GR_index_t iBF = 0; iBF < OutMesh.getNumIntBdryFaces(); iBF++) {
    BFace * bface = OutMesh.getIntBFace(iBF);
    for (GR_index_t iV = 0; iV < 2; iV++) {
      Vert * pV = bface->getVert(iV);
      if (checkedVerts.find(pV) == checkedVerts.end()) {

	std::set<Cell*> spCInc;
	std::set<Vert*> spVTmp;
	std::set<BFace*> spBFaceInc;
	bool isBdry = false;
	findNeighborhoodInfo(pV, spCInc, spVTmp, &spBFaceInc, &isBdry);
	assert(isBdry == true);
	GR_index_t nIBE = 0;
	bool qBF = false;
	BFace * bfaceHint;
	for (std::set<BFace*>::iterator itb = spBFaceInc.begin();
	    itb != spBFaceInc.end(); itb++) {
	  BFace * neighBFace = *itb;
	  if (neighBFace->getType() == Cell::eIntBdryEdge) {
	    nIBE++;
	    if (!qBF)
	      bfaceHint = neighBFace;
	  }
	  if (neighBFace->getType() == Cell::eBdryEdge) {
	    bfaceHint = neighBFace;
	    qBF = true;
	  }
	}
	assert(nIBE > 0);
	assert(bfaceHint->isValid());
	assert(
	    (qBF && bfaceHint->getType() == Cell::eBdryEdge)
		|| (!qBF && bfaceHint->getType() == Cell::eIntBdryEdge));
	if (nIBE == 1 && !qBF) {
	  pV = bfaceHint->getFace(0)->getOppositeVert(pV);
	  apexVerts.insert(std::pair<Vert*, BFace*>(pV, bfaceHint));
	  checkedVerts.insert(std::pair<Vert*, int>(pV, qBF + (nIBE)));

	}
	else if (qBF) {
	  apexVerts.insert(std::pair<Vert*, BFace*>(pV, bfaceHint));
	  checkedVerts.insert(std::pair<Vert*, int>(pV, qBF + (nIBE - 1)));

	}
	else
	  checkedVerts.insert(std::pair<Vert*, int>(pV, qBF + (nIBE - 1)));

	numNewVerts += qBF + (nIBE - 1);
      }

    }
  }
  printf("%lu apex verts, %u new Verts, %u checked verts\n", apexVerts.size(),
	 numNewVerts, checkedVerts.size());
  //	for(std::map<Vert*,BFace*>::iterator it = apexVerts.begin(); it != apexVerts.end(); ++it)
  //		(*it).first->printVertInfo();
  //	for(std::map<Vert*,int>::iterator it = checkedVerts.begin(); it != checkedVerts.end(); ++it){
  //		printf("num %d ",(*it).second);
  //		(*it).first->printVertInfo();
  //	}
  GR_index_t iBFCount = 0, iCellCount = 0, iFaceCount = 0;
  // now we can start writing the new files, its easier to check everything that will definitely
  // not be affected by doing this, otherwise we have to check against a set, rather than verttypes
  //	 I have a list of apexVerts, and a hint.

  // #Nodes, dimension of 2, 0 attributes, yes to boundary markers
  fprintf(pFNodeFile, "%u 2 0 1\n", OutMesh.getNumVerts() + numNewVerts);
  for (GR_index_t iV = 0; iV < OutMesh.getNumVerts(); ++iV) {
    Vert * pV = OutMesh.getVert(iV);
    fprintf(pFNodeFile, "%u %.16f %.16f %u\n", iV + 1, pV->x(), pV->y(),
	    pV->isBdryVert());
  }
  std::set<Cell*> cellsWritten;
  std::set<Face*> facesWritten;
  std::set<BFace*> bfacesWritten;
  // Poly File
  // no node files, dimension of 2, no attributes, yes to boundary markers
  fprintf(pFPolyFile, "0 2 0 1\n");
  // num boundary faces
  fprintf(pFPolyFile, "%u 1\n",
	  OutMesh.getNumBdryFaces() + 2 * OutMesh.getNumIntBdryFaces());
  for (GR_index_t iBF = 0; iBF < OutMesh.getNumBdryFaces(); ++iBF) {
    BFace * pBF = OutMesh.getBFace(iBF);
    Face * pF = pBF->getFace(0);
    bool qOkay = true;
    for (GR_index_t iV = 0; iV < 2; iV++) {
      Vert * pV = pF->getVert(iV);
      if (checkedVerts.count(pV)) {
	qOkay = false;
	break;
      }
    }

    if (qOkay)
      fprintf(pFPolyFile, "%u %u %u %u\n", ++iBFCount,
	      OutMesh.getVertIndex(pF->getVert(0)) + 1,
	      OutMesh.getVertIndex(pF->getVert(1)) + 1,
	      pBF->getBdryCondition());
  }

  // Ele File
  fprintf(pFEleFile, "%u 3 0\n", OutMesh.getNumCells());
  for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); ++iC) {
    Cell * pC = OutMesh.getCell(iC);
    bool qOkay = true;
    for (GR_index_t iV = 0; iV < 3; iV++) {
      Vert * pV = pC->getVert(iV);
      if (checkedVerts.count(pV)) {
	qOkay = false;
	break;
      }
    }
    if (qOkay) {
//			cellsWritten.insert(pC);
      fprintf(pFEleFile, "%u %u %u %u\n", ++iCellCount,
	      OutMesh.getVertIndex(pC->getVert(0)) + 1,
	      OutMesh.getVertIndex(pC->getVert(1)) + 1,
	      OutMesh.getVertIndex(pC->getVert(2)) + 1);
    }
  }
  // Edge File
  fprintf(pFEdgeFile, "%u 1\n", OutMesh.getNumFaces());
  for (GR_index_t iF = 0; iF < OutMesh.getNumFaces(); ++iF) {
    Face * pF = OutMesh.getFace(iF);
    bool qOkay = true;

    for (GR_index_t iV = 0; iV < 2; iV++) {
      Vert * pV = pF->getVert(iV);
      if (checkedVerts.count(pV)) {
	qOkay = false;
	break;
      }
    }

    if (qOkay) {
//			facesWritten.insert(pF);
      GR_index_t iBC = 0;
      if (pF->getCell(0)->getType() == Cell::eBdryEdge)
	iBC = dynamic_cast<BFace*>(pF->getCell(0))->getBdryCondition();
      else if (pF->getCell(1)->getType() == Cell::eBdryEdge)
	iBC = dynamic_cast<BFace*>(pF->getCell(1))->getBdryCondition();
      fprintf(pFEdgeFile, "%u %u %u %u\n", ++iFaceCount,
	      OutMesh.getVertIndex(pF->getVert(0)) + 1,
	      OutMesh.getVertIndex(pF->getVert(1)) + 1, iBC);
    }
  }

  // now all the files are ready to go

  GR_index_t pivotVertNum = OutMesh.getNumVerts();
  while (!apexVerts.empty()) {
    Vert * vert = apexVerts.begin()->first;
    BFace * bfaceHint = apexVerts.begin()->second;
    BFace * bfaceHintOld = bfaceHint;
    // need a pivot vert, which is the first vert to split
    Vert * pivotVert = vert;

    //		if(bfaceHint->getType() == Cell::eIntBdryEdge){
    //			pivotVert = bfaceHint->getFace(0)->getOppositeVert(vert);
    //		}
    Vert * pivotVertOld = NULL;
    Face * pF = bfaceHint->getFace(0);
    printf("first pivoting around ");
    pivotVert->printVertInfo();
    apexVerts.erase(pivotVert);
    std::map<Vert*, int>::iterator pivot_it = checkedVerts.find(pivotVert);
    int pivotVertCount = (*pivot_it).second;
    while (pivotVertCount > 0) {
      ((*pivot_it).second)--;

//			if(!apexVerts.count(pivotVert)) checkedVerts.erase(pivotVert);
      printf("pivoting around ");
      pivotVert->printVertInfo();
      fprintf(pFNodeFile, "%u %.16f %.16f 1\n", pivotVertNum + 1,
	      pivotVert->x(), pivotVert->y());
//			bfaceHint->printCellInfo();
      assert(bfaceHint->isValid());
      Face * pFStart = pF;
      Cell * pC = pF->getOppositeCell(bfaceHint);
      Cell * pCOld = pC;

      while (true) {

	if (!facesWritten.count(pF)) {
	  writeEdge(pFEdgeFile, pFPolyFile, OutMesh, pivotVert, NULL,
		    pivotVertNum, pF, bfaceHint->getBdryCondition(), iFaceCount,
		    iBFCount);
	  facesWritten.insert(pF);
	}
	if (pCOld->getType() == Cell::eTriCell) {
	  for (int iF = 0; iF < 3; iF++) {
	    Face *pFEdge = pC->getFace(iF);
	    assert(pFEdge->isValid());
	    if (pF == pFEdge)
	      continue;
	    if (pFEdge->getVert(0) == pivotVert
		|| pFEdge->getVert(1) == pivotVert) {
	      pF = pFEdge;
	      break;
	    }
	  }
	}

	pC = pF->getOppositeCell(pCOld);
	if (pC->getType() == Cell::eIntBdryEdge
	    || pC->getType() == Cell::eBdryEdge)
	  break;

	if (!cellsWritten.count(pCOld)) {
	  writeEle(pFEleFile, OutMesh, pivotVert, pivotVertOld, pivotVertNum,
		   pCOld, iCellCount);
	  cellsWritten.insert(pCOld);
	}
	pCOld = pC;
      } // this covers ONE side, need to fill in the LAST cell, which touches two pivot verts
	// finally do the last two, which should BOTH touch two pivotVerts;
      pivotVertOld = pivotVert;

      Face * pFInit = pF;
      pivotVert = pF->getOppositeVert(pivotVertOld);
      pivot_it = checkedVerts.find(pivotVert);
      if (pivot_it == checkedVerts.end())
	pivotVertCount = 0;
      else
	pivotVertCount = pivot_it->second;

      pivotVertNum++;
      bfaceHintOld = bfaceHint;
      bfaceHint = dynamic_cast<BFace*>(pC);

      assert(pivotVert->isBdryVert());
      if (!facesWritten.count(pF)) {
	writeEdge(pFEdgeFile, pFPolyFile, OutMesh, pivotVert, pivotVertOld,
		  pivotVertNum, pF, bfaceHint->getBdryCondition(), iFaceCount,
		  iBFCount);
	facesWritten.insert(pF);
      }
      if (!cellsWritten.count(pCOld)) {
	writeEle(pFEleFile, OutMesh, pivotVert, pivotVertOld, pivotVertNum,
		 pCOld, iCellCount);
	cellsWritten.insert(pCOld);
      }

      // now we need to fill in the other side, if we are at an int boundary edge, hop over it, and continue
      if (bfaceHint->getType() == Cell::eIntBdryEdge) {
	pF =
	    (bfaceHint->getFace(0) == pF) ?
		bfaceHint->getFace(1) : bfaceHint->getFace(0);
	pC = pF->getOppositeCell(bfaceHint);

      }
      else {
	// we need to go back and start again
//				bfaceHint->printCellInfo();
	assert(bfaceHintOld->getType() == Cell::eIntBdryEdge);
	assert(bfaceHintOld->isValid());;
//				bfaceHintOld->printCellInfo();
	pF =
	    (bfaceHintOld->getFace(0) == pFStart) ?
		bfaceHintOld->getFace(1) : bfaceHintOld->getFace(0);
	pC = pF->getOppositeCell(bfaceHintOld);
      }
      pCOld = pC;

      while (true) {

	if (!facesWritten.count(pF)) {
	  writeEdge(pFEdgeFile, pFPolyFile, OutMesh, NULL, NULL, pivotVertNum,
		    pF, bfaceHint->getBdryCondition(), iFaceCount, iBFCount);
	  facesWritten.insert(pF);
	}
	if (pCOld->getType() == Cell::eTriCell) {
	  for (int iF = 0; iF < 3; iF++) {
	    Face *pFEdge = pC->getFace(iF);
	    assert(pFEdge->isValid());
	    if (pF == pFEdge)
	      continue;
	    if (pFEdge->getVert(0) == pivotVertOld
		|| pFEdge->getVert(1) == pivotVertOld) {
	      pF = pFEdge;
	      break;
	    }
	  }
	}
//				pC->printCellInfo();

	pC = pF->getOppositeCell(pCOld);

	if (pC->getType() == Cell::eIntBdryEdge
	    || pC->getType() == Cell::eBdryEdge)
	  break;
	if (!cellsWritten.count(pCOld)) {
	  writeEle(pFEleFile, OutMesh, NULL, NULL, pivotVertNum, pCOld,
		   iCellCount);
	  cellsWritten.insert(pCOld);
	}

	pCOld = pC;
      }
      if (!facesWritten.count(pF)) {
	writeEdge(pFEdgeFile, pFPolyFile, OutMesh, NULL, NULL, pivotVertNum, pF,
		  dynamic_cast<BFace*>(pC)->getBdryCondition(), iFaceCount,
		  iBFCount);
	facesWritten.insert(pF);
      }
      if (!cellsWritten.count(pCOld)) {
	writeEle(pFEleFile, OutMesh, NULL, NULL, pivotVertNum, pCOld,
		 iCellCount);

	cellsWritten.insert(pCOld);
      }
      pF = pFInit;

      // now we are at a boundary edge, a new one.
    }			// done pivoting

  }
//	printf("Cells I'm missing\n");
//	for (GR_index_t iC = 0; iC < OutMesh.getNumCells(); ++iC){
//		Cell *pC = OutMesh.getCell(iC);
//		if(cellsWritten.count(pC) ==0 ) pC->printCellInfo();
//	}
//
//	for (GR_index_t iC = 0; iC < OutMesh.getNumFaces(); ++iC){
//		Face *pC = OutMesh.getFace(iC);
//		if(facesWritten.count(pC) ==0 ) pC->printFaceInfo();
//	}
  fprintf(pFNodeFile, "# triangle file produced by GRUMMP\n");
  fprintf(pFPolyFile, "0\n# triangle file produced by GRUMMP\n");
  fprintf(pFEleFile, "# triangle file produced by GRUMMP\n");
  fprintf(pFEdgeFile, "# triangle file produced by GRUMMP\n");

  fclose(pFNodeFile);
  fclose(pFPolyFile);
  fclose(pFEleFile);
  fclose(pFEdgeFile);

}

