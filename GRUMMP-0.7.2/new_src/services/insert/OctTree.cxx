/*
 * OctTree.cxx
 *
 *  Created on: 2014-04-13
 *      Author: zangeneh
 */
#include "GR_RefinementManager3D.h"
#include "GR_TreeConstructor.h"
#include <math.h>
#include <fstream>

namespace GRUMMP
{
  bool
  OctTree::CalcBoundingBox()
  {
    bool success = false;
    double Xminmax[2] =
      { 0 };
    double Yminmax[2] =
      { 0 };
    double Zminmax[2] =
      { 0 };
    for (GR_index_t i = 0; i < Mesh3D->getNumBdryFacesInUse(); i++) {
      BFace* BF = dynamic_cast<BFace*>(Mesh3D->getBFace(i));
      if (BF->isDeleted())
	continue;
      for (int j = 0; j < BF->getNumVerts(); j++) {
	Vert* pV = BF->getVert(j);
	assert(pV);
	assert(!pV->isDeleted());
	const double *m_loc = pV->getCoords();
	if (*m_loc <= Xminmax[0])
	  Xminmax[0] = *m_loc;
	if (*m_loc >= Xminmax[1])
	  Xminmax[1] = *m_loc;
	if (*(m_loc + 1) <= Yminmax[0])
	  Yminmax[0] = *(m_loc + 1);
	if (*(m_loc + 1) >= Yminmax[1])
	  Yminmax[1] = *(m_loc + 1);
	if (*(m_loc + 2) <= Zminmax[0])
	  Zminmax[0] = *(m_loc + 2);
	if (*(m_loc + 2) >= Zminmax[1])
	  Zminmax[1] = *(m_loc + 2);
      }
    }
    Xbound[0] = Xminmax[0];
    Xbound[1] = Xminmax[1];
    Ybound[0] = Yminmax[0];
    Ybound[1] = Yminmax[1];
    Zbound[0] = Zminmax[0];
    Zbound[1] = Zminmax[1];
    if ((Xbound[0] != 0 || Xbound[1] != 0) && (Ybound[0] != 0 || Ybound[1] != 0)
	&& (Zbound[0] != 0 || Zbound[1] != 0))
      success = true;
    return success;
  }
  int
  OctTree::Initial_OctTree_Creator(bool BoundIsGiven, int NOLayers)
  {
    if (!OctLeafList.empty())
      OctLeafList.clear();
    int Xnum = 0, Ynum = 0, Znum = 0;
    int idnum = 0;
    bool BBOK = true;
    if (!BoundIsGiven)
      BBOK = CalcBoundingBox();
    if (!BBOK)
      return 0;
    assert(MaxCircumradius != 0);
    InitialSideSize = 2 * MaxCircumradius;
    Xbound[0] -= MaxCircumradius;
    Ybound[0] -= MaxCircumradius;
    Zbound[0] -= MaxCircumradius;
    Xnum = floor(fabs(Xbound[1] - Xbound[0]) / InitialSideSize) + 1;
    Ynum = floor(fabs(Ybound[1] - Ybound[0]) / InitialSideSize) + 1;
    Znum = floor(fabs(Zbound[1] - Zbound[0]) / InitialSideSize) + 1;
    if (InitialSideSize == fabs(Xbound[1] - Xbound[0]))
      Xnum = 1;
    if (InitialSideSize == fabs(Ybound[1] - Ybound[0]))
      Ynum = 1;
    if (InitialSideSize == fabs(Zbound[1] - Zbound[0]))
      Znum = 1;
    double XYZ[3] =
      { 0 };
    double WLH[3];
    for (int ii = 0; ii < 3; ii++)
      WLH[ii] = InitialSideSize;
    for (int k = 0; k < Znum; k++) {
      XYZ[2] = Zbound[0] + k * InitialSideSize;
      for (int j = 0; j < Ynum; j++) {
	XYZ[1] = Ybound[0] + j * InitialSideSize;
	for (int i = 0; i < Xnum; i++) {
	  XYZ[0] = Xbound[0] + i * InitialSideSize;
	  Rect3D leaf(XYZ, WLH, idnum, 1);
	  OctLeafList.push_back(leaf);
	  idnum++;
	}
      }
    }
    assert(idnum == Xnum * Ynum * Znum);
    NumRects3D[0] = Xnum;
    NumRects3D[1] = Ynum;
    NumRects3D[2] = Znum;

    std::vector<Rect3D>::iterator itILS = OctLeafList.begin(), itILE =
	OctLeafList.end();
    for (; itILS != itILE; ++itILS)
      CalcBuffZone(&(*itILS));

    logMessage(2, "Number of leaves in x-direction: %d\n", Xnum);
    logMessage(2, "Number of leaves in y-direction: %d\n", Ynum);
    logMessage(2, "Number of leaves in z-direction: %d\n", Znum);

    logMessage(2, "xbounds: %f to %f\n", Xbound[0], Xbound[1]);
    logMessage(2, "Ybounds: %f to %f\n", Ybound[0], Ybound[1]);
    logMessage(2, "Zbounds: %f to %f\n", Zbound[0], Zbound[1]);
    SetNOctLayer(NOLayers);
    return (Xnum * Ynum * Znum);
  }
  int
  OctTree::OctTreeCreatorSecondPasses(int FLevel)
  {
    if (FLevel > NOctLayer || POAlter) // This number should be greater than FilterLevel defined in RefinementManager2D.
	// I know it's a naive way to do it but here it is.
	{
      WholeBufferLeaves.Clear();
      POAlter = false;
      return (NumRects3D[0] * NumRects3D[1] * NumRects3D[2]);
    }
    if (FLevel >= NOctLayer - 2)
      POAlter = true;

    std::list<Rect3D> QTCBuffer;
    std::list<Rect3D> ListBuffer;
    std::list<Rect3D> ListBufferKX[NumRects3D[1]];
    std::list<Rect3D> ListBufferKY[NumRects3D[1]];
    std::vector<Rect3D>::iterator itIL = OctLeafList.begin();
    for (int k = 0; k < NumRects3D[2]; k++) {
      for (int j = 0; j < NumRects3D[1]; j++) {
	for (int i = 0; i < NumRects3D[0]; i++) {
	  Rect3D * MLeaf = &(*itIL);
	  double XYZ[3];
	  double LWH[3];
	  double XYZDif[8][3];
	  XYZ[0] = MLeaf->x;
	  XYZ[1] = MLeaf->y;
	  XYZ[2] = MLeaf->z;
	  LWH[0] = MLeaf->length / 2;
	  LWH[1] = MLeaf->width / 2;
	  LWH[2] = MLeaf->height / 2;
	  XYZDif[0][0] = XYZ[0];
	  XYZDif[0][1] = XYZ[1];
	  XYZDif[0][2] = XYZ[2];
	  XYZDif[1][0] = XYZ[0] + LWH[0];
	  XYZDif[1][1] = XYZ[1];
	  XYZDif[1][2] = XYZ[2];
	  XYZDif[2][0] = XYZ[0];
	  XYZDif[2][1] = XYZ[1] + LWH[1];
	  XYZDif[2][2] = XYZ[2];
	  XYZDif[3][0] = XYZ[0] + LWH[0];
	  XYZDif[3][1] = XYZ[1] + LWH[1];
	  XYZDif[3][2] = XYZ[2];
	  XYZDif[4][0] = XYZ[0];
	  XYZDif[4][1] = XYZ[1];
	  XYZDif[4][2] = XYZ[2] + LWH[2];
	  XYZDif[5][0] = XYZ[0] + LWH[0];
	  XYZDif[5][1] = XYZ[1];
	  XYZDif[5][2] = XYZ[2] + LWH[2];
	  XYZDif[6][0] = XYZ[0];
	  XYZDif[6][1] = XYZ[1] + LWH[1];
	  XYZDif[6][2] = XYZ[2] + LWH[2];
	  XYZDif[7][0] = XYZ[0] + LWH[0];
	  XYZDif[7][1] = XYZ[1] + LWH[1];
	  XYZDif[7][2] = XYZ[2] + LWH[2];
	  Rect3D Leaf0(
	      XYZDif[0],
	      LWH,
	      (2 * k * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * (2 * j * NumRects3D[0] + i));
	  Rect3D Leaf1(
	      XYZDif[1],
	      LWH,
	      (2 * k * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * (2 * j * NumRects3D[0] + i) + 1);
	  Rect3D Leaf2(
	      XYZDif[2],
	      LWH,
	      (2 * k * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * ((2 * j + 1) * NumRects3D[0] + i));
	  Rect3D Leaf3(
	      XYZDif[3],
	      LWH,
	      (2 * k * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * ((2 * j + 1) * NumRects3D[0] + i) + 1);
	  Rect3D Leaf4(
	      XYZDif[4],
	      LWH,
	      ((2 * k + 1) * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * (2 * j * NumRects3D[0] + i));
	  Rect3D Leaf5(
	      XYZDif[5],
	      LWH,
	      ((2 * k + 1) * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * (2 * j * NumRects3D[0] + i) + 1);
	  Rect3D Leaf6(
	      XYZDif[6],
	      LWH,
	      ((2 * k + 1) * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * ((2 * j + 1) * NumRects3D[0] + i));
	  Rect3D Leaf7(
	      XYZDif[7],
	      LWH,
	      ((2 * k + 1) * 4 * NumRects3D[1] * NumRects3D[0])
		  + 2 * ((2 * j + 1) * NumRects3D[0] + i) + 1);
	  Rect3D *Leaves[8];
	  Leaves[0] = &Leaf0;
	  Leaves[1] = &Leaf1;
	  Leaves[2] = &Leaf2;
	  Leaves[3] = &Leaf3;
	  Leaves[4] = &Leaf4;
	  Leaves[5] = &Leaf5;
	  Leaves[6] = &Leaf6;
	  Leaves[7] = &Leaf7;
	  std::list<TetCell*>::iterator itTC = MLeaf->TetsInOct.begin();
	  // Updating the children's container with tris:
	  for (; itTC != MLeaf->TetsInOct.end(); ++itTC) {
	    int I = 0;
	    TetCell* pTC = *itTC;
	    if (!pTC->isValid())
	      continue;
	    if (pTC->isDeleted())
	      continue;
	    if (!QContainsCell(pTC, MLeaf))
	      continue;
	    bool IsQC = false;
	    for (int ii = 0; ii < 8; ii++) {
	      IsQC = QContainsCell(pTC, Leaves[ii]);
	      if (IsQC) {
		I = ii;
		Leaves[I]->TetsInOct.push_back(pTC);
		break;
	      }
	    }
	    if (!IsQC) {
	      logMessage(2, "XPosLeaf: %f\n", MLeaf->x);
	      logMessage(2, "YPosLeaf: %f\n", MLeaf->y);
	      logMessage(2, "ZPosLeaf: %f\n", MLeaf->z);
	      logMessage(2, "SideSize: %f\n", MLeaf->width);
	      logMessage(2, "SideSize: %f\n", InitialSideSize);
	      double ccc[3] =
		{ 0 };
	      pTC->calcCircumcenter(ccc);
	      logMessage(2, "XYZPos: %f %f %f\n", ccc[0], ccc[1], ccc[2]);
	      logMessage(2, "xbounds: %f to %f\n", Xbound[0], Xbound[1]);
	      logMessage(2, "Ybounds: %f to %f\n", Ybound[0], Ybound[1]);
	      logMessage(2, "Zbounds: %f to %f\n", Zbound[0], Zbound[1]);
	      assert(QContainsCell(pTC, MLeaf));
	      assert(0);
	    }
	  }
	  MLeaf->TetsInOct.clear();
	  QTCBuffer.push_back(Leaf0);
	  QTCBuffer.push_back(Leaf1);
	  ListBuffer.push_back(Leaf2);
	  ListBuffer.push_back(Leaf3);
	  ListBufferKX[j].push_back(Leaf4);
	  ListBufferKX[j].push_back(Leaf5);
	  ListBufferKY[j].push_back(Leaf6);
	  ListBufferKY[j].push_back(Leaf7);
	  ++itIL;
	}
	QTCBuffer.splice(QTCBuffer.end(), ListBuffer);
      }
      for (int counterY = 0; counterY < NumRects3D[1]; counterY++) {
	QTCBuffer.splice(QTCBuffer.end(), ListBufferKX[counterY]);
	QTCBuffer.splice(QTCBuffer.end(), ListBufferKY[counterY]);
      }

    }
    OctLeafList.clear();
    OctLeafList.insert(OctLeafList.end(), QTCBuffer.begin(), QTCBuffer.end());
    MaxCircumradius /= 4;
    InitialSideSize /= 2;
    NumRects3D[0] *= 2;
    NumRects3D[1] *= 2;
    NumRects3D[2] *= 2;
    WholeBufferLeaves.Clear();
    assert(
	OctLeafList.size()
	    == GR_index_t(NumRects3D[0] * NumRects3D[1] * NumRects3D[2]));

    std::vector<Rect3D>::iterator itILS = OctLeafList.begin(), itILE =
	OctLeafList.end();
    for (; itILS != itILE; ++itILS)
      CalcBuffZone(&(*itILS));

    return (NumRects3D[0] * NumRects3D[1] * NumRects3D[2]);
  }

  Rect3D*
  OctTree::WhichRect(const TetCell * tet)
  {
    int whichrect[3] =
      { 0 };
    double CirPos[3] =
      { 0 };
    tet->calcCircumcenter(CirPos);
    whichrect[0] = int(floor((CirPos[0] - Xbound[0]) / InitialSideSize));
    whichrect[1] = int(floor((CirPos[1] - Ybound[0]) / InitialSideSize));
    whichrect[2] = int(floor((CirPos[2] - Zbound[0]) / InitialSideSize));

    Rect3D* buffer = &(OctLeafList[whichrect[2] * NumRects3D[1] * NumRects3D[0]
	+ (whichrect[1] * NumRects3D[0] + whichrect[0])]);

    return buffer;
  }
  void
  OctTree::AddTetToRect(TetCell* tet)
  {
    assert(tet);
    Rect3D* buffer1 = WhichRect(tet);
    bool IsQT = QContainsCell(tet, buffer1);
    assert(IsQT);
    assert(buffer1);
    buffer1->TetsInOct.push_back(tet);
  }
  bool
  OctTree::QContainsCell(const TetCell * const Cell, const Rect3D * const Leaf)
  {
    assert(Cell != NULL);
    assert(!Cell->isDeleted());
    assert(Leaf);
    double PosC[3] =
      { 0 };
    bool IsIn = false;
    Cell->calcCircumcenter(PosC);
    if ((PosC[0] <= (Leaf->x + Leaf->length) && PosC[0] >= Leaf->x)
	&& (PosC[1] <= (Leaf->y + Leaf->width) && PosC[1] >= Leaf->y)
	&& (PosC[2] <= (Leaf->z + Leaf->height) && PosC[2] >= Leaf->z))
      IsIn = true;
    return IsIn;
  }
  void
  OctTree::CalcLeafIndex(const Rect3D* leaf, int &i, int &j, int &k)
  {
    assert(leaf);
    assert(leaf->LeafLayer == 1);
    k = int(double(leaf->Identifier) / double(NumRects3D[0] * NumRects3D[1]));
    assert((k + 1) * NumRects3D[0] * NumRects3D[1] >= leaf->Identifier);
    assert(k < NumRects3D[2]);
    int Remaining = leaf->Identifier - k * NumRects3D[0] * NumRects3D[1];
    j = int(double(Remaining) / double(NumRects3D[0]));
    assert(j < NumRects3D[1]);
    i = leaf->Identifier
	- (k * NumRects3D[0] * NumRects3D[1] + j * NumRects3D[0]);
    assert(i >= 0 && i < NumRects3D[0]);
  }
  int
  OctTree::FilterTet(const TetCell* tet, const Rect3D* leaf)
  {
    double sidesize = leaf->width; // with the assumption of having square rectangles
    double radius = tet->calcCircumradius();
    int ratio = floor(log2(sidesize / radius));
    return ratio;
  }
  int
  OctTree::FindIndSets(std::set<Rect3D*>& IndependentSet, bool IsNewLayer)
  {
    static GR_index_t pass = 0;
    if (!IsNewLayer) {
      IndependentSet.clear();
      if (pass >= KTInLeaves.size()) {
	pass = 0;
	WholeBufferLeaves.Clear();
	return IndependentSet.size();
      }
      std::vector<Rect3D*> IndVecBuf1 = KTInLeaves[pass];
      std::vector<Rect3D*>::iterator itInS = IndVecBuf1.begin(), itInE =
	  IndVecBuf1.end();

      IndependentSet.insert(itInS, itInE);
      pass++;
      return IndependentSet.size();
    }
    else {
      IndependentSet.clear();
      std::vector<Rect3D*> BufInds;
      Rect3D * MaxLeaf;
      // WholeBufferLeaves should be first initialised with all the leaves!
      std::map<int, Rect3D*> WBLB = WholeBufferLeaves.Queue;
      if (WBLB.empty())
	return IndependentSet.size();
      while (!WBLB.empty()) {
	std::map<int, Rect3D*>::reverse_iterator itW = WBLB.rbegin();
	MaxLeaf = itW->second;
	IndependentSet.insert(MaxLeaf);
	WholeBufferLeaves.Queue.erase(itW->first);
	WBLB.erase(itW->first);
	std::set<Rect3D*>::iterator itS = MaxLeaf->BufferZone.begin(), itE =
	    MaxLeaf->BufferZone.end();
	for (; itS != itE; ++itS) {
	  Rect3D * LeafB = *itS;
	  WBLB.erase(LeafB->SizeOrder);
	  std::set<Rect3D*>::iterator itSS = LeafB->BufferZone.begin();
	  for (; itSS != LeafB->BufferZone.end(); ++itSS) {
	    Rect3D * LeafBS = *itSS;
	    WBLB.erase(LeafBS->SizeOrder);
	    std::set<Rect3D*>::iterator itSSS = LeafBS->BufferZone.begin();
	    for (; itSSS != LeafBS->BufferZone.end(); ++itSSS) {
	      Rect3D * LeafBSS = *itSSS;
	      WBLB.erase(LeafBSS->SizeOrder);
	    }
	  }
	}
      }
      BufInds.insert(BufInds.end(), IndependentSet.begin(),
		     IndependentSet.end());
      KTInLeaves.push_back(BufInds);
      BufInds.clear();
      return IndependentSet.size();
    }
  }
  int
  OctTree::CalcBuffZone(Rect3D * const leaf)
  {
    assert(leaf);
    leaf->BufferZone.clear();
    std::vector<Rect3D*> NB;
    std::vector<Rect3D*> SB;
    std::vector<Rect3D*> EB;
    std::vector<Rect3D*> WB;
    std::vector<Rect3D*> UB;
    std::vector<Rect3D*> DB;
    // calculating 4 main directions
    NorthNeighbor(leaf, NB);
    SouthNeighbor(leaf, SB);
    WestNeighbor(leaf, WB);
    EastNeighbor(leaf, EB);
    UpNeighbor(leaf, UB);
    DownNeighbor(leaf, DB);

    // calculating orthogonal directions
    std::vector<Rect3D*> buffer;
    while (!NB.empty()) {
      Rect3D * OrthN = NB.back();
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      UpNeighbor(OrthN, buffer);
      DownNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      NB.pop_back();
    }
    while (!SB.empty()) {
      Rect3D * OrthN = SB.back();
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      UpNeighbor(OrthN, buffer);
      DownNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      SB.pop_back();
    }
    while (!EB.empty()) {
      Rect3D * OrthN = EB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      UpNeighbor(OrthN, buffer);
      DownNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      EB.pop_back();
    }
    while (!WB.empty()) {
      Rect3D * OrthN = WB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      UpNeighbor(OrthN, buffer);
      DownNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      WB.pop_back();
    }
    while (!UB.empty()) {
      Rect3D * OrthN = UB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      UB.pop_back();
    }
    while (!DB.empty()) {
      Rect3D * OrthN = DB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      DB.pop_back();
    }
    for (std::vector<Rect3D*>::iterator itt = buffer.begin();
	itt != buffer.end(); ++itt) {
      Rect3D * Orth = *itt;
      leaf->BufferZone.insert(Orth);
    }
    return leaf->BufferZone.size();
  }
  bool
  OctTree::IsOctinBufferZone(Rect3D * const leaf, Rect3D * const LeafB)
  {
//	assert(leaf);
//	assert(LeafB);
    std::set<Rect3D*>::iterator itE = leaf->BufferZone.end();
    std::set<Rect3D*>::iterator itF = leaf->BufferZone.find(LeafB);
    if (itE != itF)
      return true;
    if (leaf == LeafB)
      return true;
    return false;
  }
  int
  OctTree::NorthNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& NNeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (jj == NumRects3D[1] - 1)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < kk * NumRects3D[1] * NumRects3D[0] + (jj + 1) * NumRects3D[0] + ii;
	  i++)
	++it;
      Rect3D * FN1 = &(*it);
      if (!FN1->Wchildren) {
	NNeihbor.push_back(FN1);
	return NNeihbor.size(); // return 1;
      }
      else
	assert(0); //TODO complete this section for the cases when leaf has children

    }
    else
      assert(0); // to be done later!
    return 0;
  }
  int
  OctTree::SouthNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& SNeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (jj == 0)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < kk * NumRects3D[1] * NumRects3D[0] + (jj - 1) * NumRects3D[0] + ii;
	  i++)
	++it;
      Rect3D * FS1 = &(*it);
      if (!FS1->Wchildren) {
	SNeihbor.push_back(FS1);
	return 1;
      }
      else
	assert(0); //TODO complete this section for the cases when leaf has children
    }
    else
      assert(0); //TODO complete this section for the cases when leaf layer is not 1
    return 0;
  }
  int
  OctTree::WestNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& WNeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (ii == 0)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < kk * NumRects3D[1] * NumRects3D[0] + jj * NumRects3D[0] + ii - 1;
	  i++)
	++it;
      Rect3D * FW1 = &(*it);
      if (!FW1->Wchildren) {
	WNeihbor.push_back(FW1);
	return 1;
      }
      else
	assert(0); //TODO complete this section for the cases when leaf has children
    }
    else
      assert(0); //TODO complete this section for the cases when leaf layer is not 1
    return 0;
  }
  int
  OctTree::EastNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& ENeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (ii == NumRects3D[0] - 1)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < kk * NumRects3D[1] * NumRects3D[0] + jj * NumRects3D[0] + ii + 1;
	  i++)
	++it;
      Rect3D * FE1 = &(*it);
      if (!FE1->Wchildren) {
	ENeihbor.push_back(FE1);
	return 1;
      }
      else
	assert(0);
    }
    else
      assert(0);
    return 0;
  }
  int
  OctTree::UpNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& UNeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (kk == NumRects3D[2] - 1)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < (kk + 1) * NumRects3D[1] * NumRects3D[0] + jj * NumRects3D[0] + ii;
	  i++)
	++it;
      Rect3D * FU1 = &(*it);
      if (!FU1->Wchildren) {
	UNeihbor.push_back(FU1);
	return 1;
      }
      else
	assert(0);
    }
    else
      assert(0);
    return 0;
  }
  int
  OctTree::DownNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& DNeihbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      int ii = 0, jj = 0, kk = 0;
      CalcLeafIndex(leaf, ii, jj, kk);
      if (kk == 0)
	return 0;
      std::vector<Rect3D>::iterator it = OctLeafList.begin();
      for (int i = 0;
	  i < (kk - 1) * NumRects3D[1] * NumRects3D[0] + jj * NumRects3D[0] + ii;
	  i++)
	++it;
      Rect3D * FD1 = &(*it);
      if (!FD1->Wchildren) {
	DNeihbor.push_back(FD1);
	return 1;
      }
      else
	assert(0);
    }
    else
      assert(0);
    return 0;
  }
  void
  OctTree::writeVTKOcttree()
  {
    std::ofstream OutFileOct;
    OutFileOct.open("OctTreeMesh.vtk");
    OutFileOct << "# vtk DataFile Version 1.0" << std::endl;
    OutFileOct << "OctTree grid" << std::endl;
    OutFileOct << "ASCII" << std::endl;
    OutFileOct << "DATASET RECTILINEAR_GRID" << std::endl;
    OutFileOct << "DIMENSIONS" << " " << NumRects3D[0] + 1 << " "
	<< NumRects3D[1] + 1 << " " << NumRects3D[2] + 1 << std::endl;
    OutFileOct << "X_COORDINATES" << " " << NumRects3D[0] + 1 << " " << "double"
	<< std::endl;
    for (int i = 0; i < NumRects3D[0] + 1; i++)
      OutFileOct << Xbound[0] + i * InitialSideSize << " ";
    OutFileOct << std::endl;
    OutFileOct << "Y_COORDINATES" << " " << NumRects3D[1] + 1 << " " << "double"
	<< std::endl;
    for (int j = 0; j < NumRects3D[1] + 1; j++)
      OutFileOct << Ybound[0] + j * InitialSideSize << " ";
    OutFileOct << std::endl;
    OutFileOct << "Z_COORDINATES" << " " << NumRects3D[2] + 1 << " " << "double"
	<< std::endl;
    for (int j = 0; j < NumRects3D[1] + 1; j++)
      OutFileOct << Ybound[0] + j * InitialSideSize << " ";
  }
  bool
  OctTree::CheckEveryLeafForTets()
  {
    bool OK = false;
    for (std::vector<Rect3D>::iterator ittt = OctLeafList.begin();
	ittt != OctLeafList.end(); ++ittt) {
      for (std::list<TetCell*>::iterator itTC = (*ittt).TetsInOct.begin();
	  itTC != (*ittt).TetsInOct.end(); ++itTC) {
	OK = QContainsCell(*itTC, &(*ittt));
	if (!OK)
	  return OK;
      }
    }
    return OK;
  }

} // End of namespace
