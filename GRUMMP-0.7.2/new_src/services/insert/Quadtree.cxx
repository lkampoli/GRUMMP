/*
 * Quadtree.cxx
 *
 *  Created on: 2014-01-25
 *      Author: zangeneh
 */

#include "GR_InsertionManager.h"
#include "GR_TreeConstructor.h"
#include <math.h>
#include <fstream>

namespace GRUMMP
{
  bool
  Quadtree::CalcBoundingBox()
  {
    bool success = false;
    double Xminmax[2] =
      { 0 };
    double Yminmax[2] =
      { 0 };
    for (GR_index_t i = 0; i < Mesh->getNumBdryEdgesInUse(); i++) {
      BFace* BF = dynamic_cast<BFace*>(Mesh->getBFace(i));
      if (BF->isDeleted())
	continue;
      //assert(!BF->isDeleted());
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
      }
    }
    Xbound[0] = Xminmax[0];
    Xbound[1] = Xminmax[1];
    Ybound[0] = Yminmax[0];
    Ybound[1] = Yminmax[1];
    if ((Xbound[0] != 0 || Xbound[1] != 0)
	&& (Ybound[0] != 0 || Ybound[1] != 0))
      success = true;
    return success;
  }

  int
  Quadtree::Initial_Quadtree_Creator(bool BoundIsGiven, int NQLayers)
  {
    if (!initial_quads.empty())
      initial_quads.clear();
    int Xnum = 0;
    int Ynum = 0;
    int idnum = 0;
    bool BBOK = true;
    if (!BoundIsGiven)
      BBOK = CalcBoundingBox();
    if (!BBOK)
      return 0;
    assert(MaxCircumradius != 0);
    InitialSideSize = 2 * MaxCircumradius;
    Xnum = floor(fabs(Xbound[1] - Xbound[0]) / InitialSideSize) + 1;
    Ynum = floor(fabs(Ybound[1] - Ybound[0]) / InitialSideSize) + 1;
    if (InitialSideSize == fabs(Xbound[1] - Xbound[0]))
      Xnum = 1;
    if (InitialSideSize == fabs(Ybound[1] - Ybound[0]))
      Ynum = 1;
    for (int j = 0; j < Ynum; j++) {
      double y = Ybound[0] + j * InitialSideSize;
      for (int i = 0; i < Xnum; i++) {
	double x = Xbound[0] + i * InitialSideSize;
	Rect leaf(x, y, InitialSideSize, InitialSideSize, idnum, 1);
	initial_quads.push_back(leaf);
	idnum++;
      }
    }
    assert(idnum == Xnum * Ynum);
    NumRects[0] = Xnum;
    NumRects[1] = Ynum;

    std::vector<Rect>::iterator itILS = initial_quads.begin(), itILE =
	initial_quads.end();
    for (; itILS != itILE; itILS++)
      CalcBuffZone(&(*itILS));

    logMessage(2, "Number of leaves in x-direction: %d\n", Xnum);
    logMessage(2, "Number of leaves in y-direction: %d\n", Ynum);

    SetNQuadLayer(NQLayers);

    return (Xnum * Ynum);
  }
  Rect*
  Quadtree::WhichRect(const TriCell * tri)
  {
    int whichrect[2] =
      { 0 };
    double CirPos[2] =
      { 0 };
    tri->calcCircumcenter(CirPos);
    whichrect[0] = int(floor((CirPos[0] - Xbound[0]) / InitialSideSize));
    whichrect[1] = int(floor((CirPos[1] - Ybound[0]) / InitialSideSize));

    Rect* buffer = &(initial_quads[whichrect[1] * NumRects[0] + whichrect[0]]);

    return buffer;
  }

  void
  Quadtree::AddTriToRect(TriCell* tri)
  {
    assert(tri);
    Rect * buffer1 = WhichRect(tri);
    bool IsQT = QContainsCell(tri, buffer1);
    assert(IsQT);
    //TODO this way of accessing the iterator of initial quads needs to be changed
    //as reallocating the initial_quads vector on the fly could invalidate
    //the pointers to the specific Tris storing in each Rect
//	std::list<Rect>::iterator it = initial_quads.begin();
//	for (int i=0 ; i<whichrect[1]*NumRects[0]+whichrect[0] ; i++) it++;
//	Rect* buffer2;
//	buffer2 = &(*it);
//	assert(buffer1 == buffer2);
    assert(buffer1);
    buffer1->TrisInQuad.push_back(tri);
  }

  bool
  Quadtree::CirIntersectLeaf(const TriCell* tri, Rect* leaf)
  {
    assert(leaf);
    assert(tri);
    double CirPos[2];
    double CirRadius;
    double distance = 0, dx = 0, dy = 0;
    tri->calcCircumcenter(CirPos);
    CirRadius = tri->calcCircumradius();
    if (CirPos[0] >= leaf->x && CirPos[0] <= leaf->x)
      dx = 0;
    else if (CirPos[0] < leaf->x)
      dx = leaf->x - CirPos[0];
    else if (CirPos[0] > leaf->x + leaf->width)
      dx = CirPos[0] - (leaf->x + leaf->width);
    else
      assert(0);

    if (CirPos[1] >= leaf->y && CirPos[1] <= leaf->y)
      dy = 0;
    else if (CirPos[1] < leaf->y)
      dy = leaf->y - CirPos[1];
    else if (CirPos[1] > leaf->y + leaf->height)
      dy = CirPos[1] - (leaf->y + leaf->height);
    else
      assert(0);

    distance = sqrt(dx * dx + dy * dy);
    if (distance <= CirRadius)
      return true;
    else
      return false;
  }

  bool
  Quadtree::SubDivide(Rect* MLeaf)
  {
    assert(MLeaf);
    assert(!MLeaf->Wchildren);
    if (MLeaf->LeafLayer >= MaxLeafLayer)
      return false;
    MLeaf->Wchildren = true;
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
	double x = MLeaf->x + i * MLeaf->width;
	double y = MLeaf->y + j * MLeaf->height;
	Rect child(x, y, MLeaf->width / 2, MLeaf->height / 2,
		   MLeaf->Identifier * 10 + 2 * j + i, MLeaf->LeafLayer + 1);
	child.Mom = MLeaf;
	child.Wmother = true;
      }
    }
    return true;
  }
/// This subroutine recreates the QT structure by subdividing each leaf to 4 of its children
/// it then updates the children's TrisInQuad containers by their associated parents.
/// It also takes care of the "MaxCircumradius" and "InitialSideSize": ...
/// MaxCircumradius /= 4;     InitialSideSize /= 2;
/// the bounding coordinates would stay the same.
  int
  Quadtree::QuadTreeCreatorSecondPasses(int FLevel)
  {
    if (FLevel > NQuadLayer || PQAlter) // This number should be greater than FilterLevel defined in RefinementManager2D.
	// I know it's a naive way to do it but here it is.
	{
      WholeBufferLeaves.Clear();
      PQAlter = false;
      return (NumRects[0] * NumRects[1]);
    }
    if (FLevel >= NQuadLayer - 2)
      PQAlter = true;
    std::list<Rect> QTCBuffer;
    std::list<Rect> ListBuffer;
    std::vector<Rect>::iterator itIL = initial_quads.begin();
    for (int j = 0; j < NumRects[1]; j++) {
      for (int i = 0; i < NumRects[0]; i++) {
	Rect * MLeaf = &(*itIL);
	double XPos = MLeaf->x;
	double YPos = MLeaf->y;
	double height = MLeaf->height / 2;
	double width = MLeaf->width / 2;
	Rect Leaf0(XPos, YPos, width, height, 2 * (2 * j * NumRects[0] + i));
	Rect Leaf1(XPos + width, YPos, width, height,
		   2 * (2 * j * NumRects[0] + i) + 1);
	Rect Leaf2(XPos, YPos + height, width, height,
		   2 * ((2 * j + 1) * NumRects[0] + i));
	Rect Leaf3(XPos + width, YPos + height, width, height,
		   2 * ((2 * j + 1) * NumRects[0] + i) + 1);
	Rect *Leaves[4];
	Leaves[0] = &Leaf0;
	Leaves[1] = &Leaf1;
	Leaves[2] = &Leaf2;
	Leaves[3] = &Leaf3;
	std::list<TriCell*>::iterator itTC = MLeaf->TrisInQuad.begin();
	// Updating the children's container with tris:
	for (; itTC != MLeaf->TrisInQuad.end(); itTC++) {
	  int I = 0;
	  TriCell * pTC = *itTC;
	  if (!pTC->isValid())
	    continue;
	  if (pTC->isDeleted())
	    continue;
	  if (!QContainsCell(pTC, MLeaf))
	    continue;
	  bool IsQC = false;
	  for (int ii = 0; ii < 4; ii++) {
	    IsQC = QContainsCell(pTC, Leaves[ii]);
	    if (IsQC) {
	      I = ii;
	      Leaves[I]->TrisInQuad.push_back(pTC);
	      break;
	    }
	  }
	  if (!IsQC) {
	    logMessage(2, "XPosLeaf: %f\n", MLeaf->x);
	    logMessage(2, "YPosLeaf: %f\n", MLeaf->y);
	    logMessage(2, "SideSize: %f\n", MLeaf->width);
	    logMessage(2, "SideSize: %f\n", InitialSideSize);
	    double ccc[2] =
	      { 0 };
	    pTC->calcCircumcenter(ccc);
	    logMessage(2, "XYPos: %f %f\n", ccc[0], ccc[1]);
	    logMessage(2, "xbounds: %f to %f\n", Xbound[0], Xbound[1]);
	    logMessage(2, "Ybounds: %f to %f\n", Ybound[0], Ybound[1]);
	    assert(QContainsCell(pTC, MLeaf));
	  }
	  assert(IsQC);
	}
//			std::vector<BFace*>::iterator ItBS = MLeaf->m_newBdryFaces.begin(),
//					ItBE = MLeaf->m_newBdryFaces.end();
//			for ( ; ItBS != ItBE ; ItBS++)
//			{
//				int I = 0;
//				BFace* BF = *ItBS;
//				if (!BF->isValid() || BF->isDeleted()) continue;
//				Face * pF = BF->getFace();
//				TriCell* pTC = dynamic_cast<TriCell*>(pF->getOppositeCell(BF));
//				if (!pTC->isValid() || pTC->isDeleted()) continue;
//				bool IsQC = false;
//				for (int ii=0 ; ii<4 ; ii++)
//				{
//					IsQC = QContainsCell(pTC, Leaves[ii]);
//					if (IsQC)
//					{
//						I = ii;
//						Leaves[I]->TrisInQuad.push_back(pTC);
//						break;
//					}
//				}
//				if (!IsQC)
//					assert(0);
//			}
	MLeaf->TrisInQuad.clear();
	MLeaf->m_newBdryFaces.clear();
	QTCBuffer.push_back(Leaf0);
	QTCBuffer.push_back(Leaf1);
	ListBuffer.push_back(Leaf2);
	ListBuffer.push_back(Leaf3);
	itIL++;
      }
      QTCBuffer.splice(QTCBuffer.end(), ListBuffer);
    }
    initial_quads.clear();
    for (std::list<Rect>::iterator itC = QTCBuffer.begin();
	itC != QTCBuffer.end(); itC++)
      initial_quads.push_back(*itC);
    MaxCircumradius /= 4;
    InitialSideSize /= 2;
    NumRects[0] *= 2;
    NumRects[1] *= 2;
    WholeBufferLeaves.Clear();

    std::vector<Rect>::iterator itILS = initial_quads.begin(), itILE =
	initial_quads.end();
    for (; itILS != itILE; itILS++)
      CalcBuffZone(&(*itILS));

    return (NumRects[0] * NumRects[1]);
  }
  bool
  Quadtree::QContainsCell(const TriCell * const Cell, const Rect * const Leaf)
  {
    double PosC[2] =
      { 0 };
    bool IsIn = false;
    Cell->calcCircumcenter(PosC);
    if ((PosC[0] <= (Leaf->x + Leaf->width) && PosC[0] >= Leaf->x)
	&& (PosC[1] <= (Leaf->y + Leaf->width) && PosC[1] >= Leaf->y))
      IsIn = true;

    return IsIn;
  }
  void
  Quadtree::CalcLeafIndex(const Rect* leaf, int &i, int &j)
  {
    j = int(double(leaf->Identifier) / double(NumRects[0]));
    assert((j + 1) * NumRects[0] >= leaf->Identifier);
    assert(j < NumRects[1]);
    i = leaf->Identifier - j * NumRects[0];
    assert(i >= 0 && i < NumRects[0]);
  }

  int
  Quadtree::FilterTriangle(const TriCell* tri, const Rect* leaf)
  {
    double sidesize = leaf->width; // with the assumption of having square rectangles
    double radius = tri->calcCircumradius();
    int ratio = floor(log2(sidesize / radius));
    return ratio;
  }

  int
  Quadtree::FindIndSets(std::set<Rect*>& IndependentSet, bool IsNewLayer)
  {
    static GR_index_t pass = 0;
    if (!IsNewLayer) {
      IndependentSet.clear();
      if (pass >= KTInLeaves.size()) {
	pass = 0;
	WholeBufferLeaves.Clear();
	return 0;
      }
      std::vector<Rect*> IndVecBuf1 = KTInLeaves[pass];
      std::vector<Rect*>::iterator itInS = IndVecBuf1.begin(), itInE =
	  IndVecBuf1.end();
      IndependentSet.insert(itInS, itInE);
      pass++;
      return 0;
    }
    else {
      IndependentSet.clear();
      std::vector<Rect*> BufInds;
      Rect * MaxLeaf;
      // WholeBufferLeaves should be first initialised with all the leaves!
      std::map<int, Rect*> WBLB = WholeBufferLeaves.Queue;
      if (WBLB.empty())
	return IndependentSet.size();
      while (!WBLB.empty()) {
	std::map<int, Rect*>::reverse_iterator itW = WBLB.rbegin();
	MaxLeaf = itW->second;
	IndependentSet.insert(MaxLeaf);
	WholeBufferLeaves.Queue.erase(itW->first);
	WBLB.erase(itW->first);
	for (std::set<Rect*>::iterator itS = MaxLeaf->BufferZone.begin();
	    itS != MaxLeaf->BufferZone.end(); itS++) {
	  Rect * LeafB = *itS;
	  WBLB.erase(LeafB->SizeOrder);
	  for (std::set<Rect*>::iterator itSS = (*itS)->BufferZone.begin();
	      itSS != (*itS)->BufferZone.end(); itSS++) {
	    Rect * LeafBB = *itSS;
	    WBLB.erase(LeafBB->SizeOrder);
//					for (std::set<Rect*>::iterator itSSS = (*itSS)->BufferZone.begin();
//							itSSS != (*itSS)->BufferZone.end(); itSSS++){
//						Rect * LeafBBB = *itSSS;
//						WBLB.erase(LeafBBB->SizeOrder);
//					}
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
  Quadtree::CalcBuffZone(Rect * const leaf)
  {
    assert(leaf);
    leaf->BufferZone.clear();
    std::vector<Rect*> NB;
    std::vector<Rect*> SB;
    std::vector<Rect*> EB;
    std::vector<Rect*> WB;
    // calculating 4 main directions
    NorthNeighbor(leaf, NB);
    SouthNeighbor(leaf, SB);
    WestNeighbor(leaf, WB);
    EastNeighbor(leaf, EB);

    // calculating orthogonal directions
    std::vector<Rect*> buffer;
    while (!NB.empty()) {
      Rect * OrthN = NB.back();
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      NB.pop_back();
    }
    while (!SB.empty()) {
      Rect * OrthN = SB.back();
      EastNeighbor(OrthN, buffer);
      WestNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      SB.pop_back();
    }
    while (!EB.empty()) {
      Rect * OrthN = EB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      EB.pop_back();
    }
    while (!WB.empty()) {
      Rect * OrthN = WB.back();
      NorthNeighbor(OrthN, buffer);
      SouthNeighbor(OrthN, buffer);
      leaf->BufferZone.insert(OrthN);
      WB.pop_back();
    }

    for (std::vector<Rect*>::iterator itt = buffer.begin(); itt != buffer.end();
	itt++) {
      Rect * Orth = *itt;
      leaf->BufferZone.insert(Orth);
    }
    return leaf->BufferZone.size();
  }

  bool
  Quadtree::IsQuadinBufferZone(const Rect * const leaf, Rect * const LeafB)
  {
    assert(!leaf->BufferZone.empty());
    std::set<Rect*>::iterator itE = leaf->BufferZone.end();
    std::set<Rect*>::iterator itF = leaf->BufferZone.find(LeafB);
    if (itE != itF)
      return true;
    if (leaf == LeafB)
      return true;
    return false;
  }

  int
  Quadtree::EastNeighbor(const Rect* leaf, std::vector<Rect*>& ENeighbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0;
      CalcLeafIndex(leaf, ii, jj);
      if (ii == NumRects[0] - 1)
	return 0;
      std::vector<Rect>::iterator it = initial_quads.begin();
      for (int i = 0; i < jj * NumRects[0] + ii + 1; i++)
	it++;
      Rect * FE1 = &(*it);
      if (!FE1->Wchildren) {
	ENeighbor.push_back(FE1);
	return 1;
      }
      std::list<Rect*> buff;
      std::list<Rect*>::iterator itb = buff.begin();
      buff.push_back(FE1->Children[0]);
      buff.push_back(FE1->Children[2]);
      while (!buff.empty()) {
	itb = buff.begin();
	Rect * FNB = *itb;
	if (!FNB->Wchildren) {
	  ENeighbor.push_back(FNB);
	  buff.pop_front();
	}
	else {
	  buff.pop_front();
	  buff.push_back(FNB->Children[0]);
	  buff.push_back(FNB->Children[2]);
	}
      }
      return ENeighbor.size();
    }
    assert(0);
    return 0;
  }

  int
  Quadtree::WestNeighbor(const Rect* leaf, std::vector<Rect*>& WNeighbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0;
      CalcLeafIndex(leaf, ii, jj);
      if (ii == 0)
	return 0;
      std::vector<Rect>::iterator it = initial_quads.begin();
      for (int i = 0; i < jj * NumRects[0] + ii - 1; i++)
	it++;
      Rect * FW1 = &(*it);
      if (!FW1->Wchildren) {
	WNeighbor.push_back(FW1);
	return 1;
      }
      std::list<Rect*> buff;
      std::list<Rect*>::iterator itb = buff.begin();
      buff.push_back(FW1->Children[1]);
      buff.push_back(FW1->Children[3]);
      while (!buff.empty()) {
	itb = buff.begin();
	Rect * FNB = *itb;
	if (!FNB->Wchildren) {
	  WNeighbor.push_back(FNB);
	  buff.pop_front();
	}
	else {
	  buff.pop_front();
	  buff.push_back(FNB->Children[1]);
	  buff.push_back(FNB->Children[3]);
	}
      }
      return WNeighbor.size();
    }
    else
      assert(0);
    assert(0);
    return (0);
  }

  int
  Quadtree::SouthNeighbor(const Rect* leaf, std::vector<Rect*>& SNeighbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      assert(!leaf->Wmother);
      int ii = 0, jj = 0;
      CalcLeafIndex(leaf, ii, jj);
      if (jj == 0)
	return 0;
      std::vector<Rect>::iterator it = initial_quads.begin();
      for (int i = 0; i < (jj - 1) * NumRects[0] + ii; i++)
	it++;
      Rect * FN1 = &(*it);
      if (!FN1->Wchildren) {
	SNeighbor.push_back(FN1);
	return 1;
      }
      std::list<Rect*> buff;
      std::list<Rect*>::iterator itb = buff.begin();
      buff.push_back(FN1->Children[2]);
      buff.push_back(FN1->Children[3]);
      while (!buff.empty()) {
	itb = buff.begin();
	Rect * FNB = *itb;
	if (!FNB->Wchildren) {
	  SNeighbor.push_back(FNB);
	  buff.pop_front();
	}
	else {
	  buff.pop_front();
	  buff.push_back(FNB->Children[2]);
	  buff.push_back(FNB->Children[3]);
	}
      }
      return SNeighbor.size();

    }
    else
      assert(0); // TODO to be done later. this section is responsible when we have more than 1 layer of leaves.
    assert(0);
    return 0;
  }

  int
  Quadtree::NorthNeighbor(const Rect* leaf, std::vector<Rect*>& NNeighbor)
  {
    assert(leaf);
    if (leaf->LeafLayer == 1) {
      int jj = 0;
      int ii = 0;
      while (jj * NumRects[0] <= leaf->Identifier)
	jj++;
      jj--;
      ii = leaf->Identifier - jj * NumRects[0];
      if (jj == NumRects[1] - 1)
	return 0;
      std::vector<Rect>::iterator it = initial_quads.begin();
      for (int i = 0; i < (jj + 1) * NumRects[0] + ii; i++)
	it++;
      Rect * FN1 = &(*it);
      if (!FN1->Wchildren) {
	NNeighbor.push_back(FN1);
	return 1;
      }
      std::list<Rect*> buff;
      std::list<Rect*>::iterator itb = buff.begin();
      buff.push_back(FN1->Children[0]);
      buff.push_back(FN1->Children[1]);
      while (!buff.empty()) {
	itb = buff.begin();
	Rect * FNB = *itb;
	if (!FNB->Wchildren) {
	  NNeighbor.push_back(FNB);
	  buff.pop_front();
	}
	else {
	  buff.push_back(FNB->Children[0]);
	  buff.push_back(FNB->Children[1]);
	  buff.pop_front();
	}
      }
      return NNeighbor.size();
    }
    else {
//		assert(0); // TODO The coding is done for this part but it needs thorough debugging and testing. It got really lengthy!
//		assert(leaf->Wmother);
      int dec = leaf->Identifier % 10;
      switch (dec)
	{
	case 0:
	  Rect * FN1;
	  FN1 = leaf->Mom->Children[2];
	  if (!FN1->Wchildren) {
	    NNeighbor.push_back(FN1);
	    return 1;
	  }
	  else {
	    std::list<Rect*> buff;
	    std::list<Rect*>::iterator it = buff.begin();
	    buff.push_back(FN1->Children[0]);
	    buff.push_back(FN1->Children[1]);
	    while (!buff.empty()) {
	      it = buff.begin();
	      Rect * FNB = *it;
	      if (!FNB->Wchildren) {
		NNeighbor.push_back(FNB);
		buff.pop_front();
	      }
	      else {
		buff.push_back(FNB->Children[0]);
		buff.push_back(FNB->Children[1]);
		buff.pop_front();
	      }
	    }
	    return NNeighbor.size();
	  }
	  break;
	case 1:
	  Rect * FN2;
	  FN2 = leaf->Mom->Children[3];
	  if (!FN2->Wchildren) {
	    NNeighbor.push_back(FN2);
	    return 1;
	  }
	  else {
	    std::list<Rect*> buff;
	    std::list<Rect*>::iterator it = buff.begin();
	    buff.push_back(FN2->Children[0]);
	    buff.push_back(FN2->Children[1]);
	    while (!buff.empty()) {
	      it = buff.begin();
	      Rect * FNB = *it;
	      if (!FNB->Wchildren) {
		NNeighbor.push_back(FNB);
		buff.pop_front();
	      }
	      else {
		buff.push_back(FNB->Children[0]);
		buff.push_back(FNB->Children[1]);
		buff.pop_front();
	      }
	    }
	    return NNeighbor.size();
	  }
	  break;
	case 2:
	case 3:
	  {
	    assert(0);  // TODO
	    Rect * GrandMom = leaf->Mom;
	    std::vector<Rect*> Mothers;
	    Mothers.push_back(GrandMom);
	    while (GrandMom->LeafLayer > 1
		&& GrandMom != GrandMom->Mom->Children[0]
		&& GrandMom != GrandMom->Mom->Children[1]) {
	      GrandMom = GrandMom->Mom;
	      Mothers.push_back(GrandMom);
	    }
	    if (GrandMom->LeafLayer == 1) {
	      int jj = 0;
	      int ii = 0;
	      while (jj * NumRects[0] < leaf->Identifier)
		jj++;
	      assert((jj + 1) * NumRects[0] >= leaf->Identifier);
	      assert(jj < NumRects[1]);
	      if (jj == NumRects[1] - 1)
		return 0;
	      ii = leaf->Identifier - jj * NumRects[0];
	      assert(ii >= 0 && ii < NumRects[0]);
	      std::vector<Rect>::iterator it = initial_quads.begin();
	      for (int i = 0; i < (jj + 1) * NumRects[0] + ii; i++)
		it++;
	      Rect * NGrandMom = &(*it);
	      Mothers.pop_back();
	      std::list<Rect*> buff;
	      buff.push_back(NGrandMom);
	      std::list<Rect*>::iterator itb = buff.begin();
	      while (!buff.empty()) {
		itb = buff.begin();
		Rect * FNB = *itb;
		if (!FNB->Wchildren) {
		  NNeighbor.push_back(FNB);
		  buff.pop_front();
		}
		else {
		  int BID = -1; // some arbitrary number to show we are going down further than the leaf size.
		  if (!Mothers.empty()) {
		    BID = Mothers.back()->Identifier % 10;
		    Mothers.pop_back();
		  }
		  if (BID == 2)
		    buff.push_back(FNB->Children[0]);
		  else if (BID == 3)
		    buff.push_back(FNB->Children[1]);
		  else if (BID == -1) {
		    buff.push_back(FNB->Children[0]);
		    buff.push_back(FNB->Children[1]);
		  }
		  else
		    assert(0); // we should not end up here. if we reach here something's wrong.
		}
	      } // end of while
	    }
	    else { // this else is for the case that the root mother is not in the first layer.
	      assert(GrandMom->Wmother);
	      Rect * NGrandMom;
	      if (GrandMom == GrandMom->Mom->Children[0])
		NGrandMom = GrandMom->Mom->Children[2];
	      else if (GrandMom == GrandMom->Mom->Children[1])
		NGrandMom = GrandMom->Mom->Children[3];
	      else
		assert(0); // we should not end up here!
	      Mothers.pop_back();
	      std::list<Rect*> buff;
	      buff.push_back(NGrandMom);
	      std::list<Rect*>::iterator itb = buff.begin();
	      while (!buff.empty()) {
		itb = buff.begin();
		Rect * FNB = *itb;
		if (!FNB->Wchildren) {
		  NNeighbor.push_back(FNB);
		  buff.pop_front();
		}
		else {
		  int BID = -1; // some arbitrary number to show we are going down further than the leaf size.
		  if (!Mothers.empty()) {
		    BID = Mothers.back()->Identifier % 10;
		    Mothers.pop_back();
		  }
		  if (BID == 2)
		    buff.push_back(FNB->Children[0]);
		  else if (BID == 3)
		    buff.push_back(FNB->Children[1]);
		  else if (BID == -1) {
		    buff.push_back(FNB->Children[0]);
		    buff.push_back(FNB->Children[1]);
		  }
		  else
		    assert(0); // we should not end up here. if we reach here something's wrong.
		}
	      } // end of while
	    }
	    return NNeighbor.size();
	    // no break needed
	  }
	default:
	  logMessage(1, "leaf identifier is faulty. check subdivide function");
	  assert(0);
	  break;
	}
    }
    assert(0);
    return 0;
  }

  void
  Quadtree::writeVTKquadtree()
  {
    std::ofstream outfile2;
    outfile2.open("QuadtreeMesh.vtk");
    outfile2 << "# vtk DataFile Version 1.0" << std::endl;
    outfile2 << "Quadtree grid" << std::endl;
    outfile2 << "ASCII" << std::endl;
    outfile2 << "DATASET RECTILINEAR_GRID" << std::endl;
    outfile2 << "DIMENSIONS" << " " << NumRects[0] + 1 << " " << NumRects[1] + 1
	<< " " << 1 << std::endl;
    outfile2 << "X_COORDINATES" << " " << NumRects[0] + 1 << " " << "double"
	<< std::endl;
    for (int i = 0; i < NumRects[0] + 1; i++)
      outfile2 << Xbound[0] + i * InitialSideSize << " ";
    outfile2 << std::endl;
    outfile2 << "Y_COORDINATES" << " " << NumRects[1] + 1 << " " << "double"
	<< std::endl;
    for (int j = 0; j < NumRects[1] + 1; j++)
      outfile2 << Ybound[0] + j * InitialSideSize << " ";
    outfile2 << std::endl;
    outfile2 << "Z_COORDINATES 1 double" << std::endl;
    outfile2 << "0.0";
//	outfile2<<"Ybound[0] is:"<<Ybound[0]<<std::endl;
//	outfile2<<"Ybound[1] is:"<<Ybound[1]<<std::endl;
  }

} //end of namespace

