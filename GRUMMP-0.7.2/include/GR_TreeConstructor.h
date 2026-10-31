/*
 * GR_TreeConstructor.h
 *
 *  Created on: 2014-01-24
 *      Author: zangeneh
 */

#ifndef GR_TREECONSTRUCTOR_H_
#define GR_TREECONSTRUCTOR_H_

#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"
#include <omp.h>
namespace GRUMMP
{

#define MaxLeafLayer 3

  struct XY {
    double x, y, r;
    XY(double X, double Y, double R = 0) :
	x(X), y(Y), r(R)
    {
    }
    XY() :
	x(0), y(0), r(0)
    {
    }
  };

  struct Rect {
    double x, y, width, height;
    std::list<TriCell*> TrisInQuad;
    bool IsLocked;
    bool Wmother;
    bool Wchildren;
    int SizeOrder;
    int LeafLayer;
    int Identifier;
    Rect * Mom;
    Rect * Children[4];
    std::set<Rect*> BufferZone;
    std::vector<BFace*> m_newBdryFaces;
    Rect() :
	x(0), y(0), width(0), height(0), TrisInQuad(), IsLocked(false), Wmother(
	    false), Wchildren(false), SizeOrder(0)
    {
      Mom = NULL;
      Children[0] = NULL;
      Children[1] = NULL;
      Children[2] = NULL;
      Children[3] = NULL;
    }
    Rect(double X, double Y, double W, double H, long int id, int layer = 1) :
	x(X), y(Y), width(W), height(H), TrisInQuad(), IsLocked(false), Wmother(
	    false), Wchildren(false), SizeOrder(0), LeafLayer(layer), Identifier(
	    id)
    {
      Mom = NULL;
      Children[0] = NULL;
      Children[1] = NULL;
      Children[2] = NULL;
      Children[3] = NULL;
    }
    inline void
    Push_Back_Cell(TriCell * tri)
    {
#pragma omp critical(PushBackCells)
      {
	this->TrisInQuad.push_back(tri);
      }
    }

    inline void
    SetLock(bool lock)
    {
      //#pragma omp atomic
      IsLocked = lock;
    }
    inline bool
    isLocked()
    {
      bool buff;
      //#pragma omp atomic
      buff = IsLocked;
      return buff;
    }
    void
    UpdateCellLists(TriCell* pTC)
    {
#pragma omp critical(RecieveCell)
      {
	TrisInQuad.push_back(pTC);
      }
    }
  };

  class Quadtree {
  private:
    Mesh2D* Mesh;
    double InitialSideSize;
    bool PQAlter;
    double MaxCircumradius;
    double Xbound[2];
    double Ybound[2];
    /// NumRects[0] is the number of rects in x direction.
    /// NumRects[1] is the number of rects in y direction.
    int NumRects[2];
//	std::set<Rect*> WholeBufferLeaves;
    std::vector<std::vector<Rect*> > KTInLeaves;
    //Quadtree( const Quadtree&);
    //Quadtree& operator=( const Quadtree&);
    //friend class RefinementManager2D;
  public:
    std::vector<Rect> initial_quads;
    int NQuadLayer;
    Quadtree(Mesh2D* meshToQuad) :
	Mesh(meshToQuad), InitialSideSize(0), PQAlter(false), MaxCircumradius(
	    0), initial_quads(), NQuadLayer(5)
    {
      Xbound[0] = 0;
      Xbound[1] = 0;
      Ybound[0] = 0;
      Ybound[1] = 0;
      NumRects[0] = 0;
      NumRects[1] = 0;
    }
    struct OrderedQueue {
//		OrderedQueue& OrderedQueue::operator=(const OrderedQueue& X){
//			Queue = X.Queue;
//			return *this;
//		}
      std::map<int, Rect*> Queue;
      std::multimap<int, Rect*> MMQueue;
      void
      Insert(Rect* leaf)
      {
	int dkey = leaf->TrisInQuad.size();
	std::pair<int, Rect*> input = std::make_pair(dkey, leaf);
	MMQueue.insert(input);
      }
      void
      Update()
      {
	std::multimap<int, Rect*>::iterator itS = MMQueue.begin(), itE =
	    MMQueue.end();
	Queue.clear();
	int num = 0;
	for (; itS != itE; itS++) {
	  num++;
	  itS->second->SizeOrder = num;
	  std::pair<int, Rect*> input = std::make_pair(num, itS->second);
	  Queue.insert(input);
	}
	MMQueue.clear();
      }
//		void Delete(Rect* leaf){
////			Queue.remove(leaf);
//		}
      void
      Clear()
      {
	Queue.clear();
	MMQueue.clear();
      }
    };
    OrderedQueue WholeBufferLeaves;
    inline void
    SetNQuadLayer(int NL)
    {
      NQuadLayer = NL;
    }
    bool
    CalcBoundingBox();
    inline void
    UpdateMaxCirR(double CircR)
    {
      MaxCircumradius = CircR;
    }
    inline void
    UpdateBBound(double XB[2], double YB[2])
    {
      if (XB[0] < Xbound[0])
	Xbound[0] = XB[0];
      if (XB[1] > Xbound[1])
	Xbound[1] = XB[1];
      if (YB[0] < Ybound[0])
	Ybound[0] = YB[0];
      if (YB[1] > Ybound[1])
	Ybound[1] = YB[1];
    }
    int
    Initial_Quadtree_Creator(bool BoundIsGiven = true, int NQLayers = 6);
    int
    QuadTreeCreatorSecondPasses(int FLevel);
    Rect*
    WhichRect(const TriCell* tri);
    void
    AddTriToRect(TriCell* tri);
    bool
    SubDivide(Rect* MLeaf);
    void
    writeVTKquadtree();
    inline void
    clearKTInLeaves()
    {
      KTInLeaves.clear();
    }
    bool
    CirIntersectLeaf(const TriCell* tri, Rect* leaf);
    int
    CalcBuffZone(Rect * const leaf);
    int
    NorthNeighbor(const Rect* leaf, std::vector<Rect*>& NNeihbor);
    int
    SouthNeighbor(const Rect* leaf, std::vector<Rect*>& SNeihbor);
    int
    WestNeighbor(const Rect* leaf, std::vector<Rect*>& WNeihbor);
    int
    EastNeighbor(const Rect* leaf, std::vector<Rect*>& ENeihbor);
    int
    FindIndSets(std::set<Rect*>& IndependentSet, bool IsNewLayer = true);
    void
    CalcLeafIndex(const Rect* leaf, int &i, int &j);
    int
    FilterTriangle(const TriCell* tri, const Rect* leaf);
    bool
    IsQuadinBufferZone(const Rect * const leaf, Rect * const LeafB);
    inline bool
    IsWBLEmpty()
    {
      return WholeBufferLeaves.Queue.empty();
    }
    inline void
    FillWBZ(bool NewQuad = false)
    {
      if (!NewQuad)
	WholeBufferLeaves.Queue.insert(
	    std::pair<int, Rect*>(1, &(*initial_quads.begin())));
      WholeBufferLeaves.MMQueue.clear();
      std::vector<Rect>::iterator itW = initial_quads.begin();
      for (; itW != initial_quads.end(); itW++)
	WholeBufferLeaves.Insert(&(*itW));
      WholeBufferLeaves.Update();
    }
    bool
    QContainsCell(const TriCell * const Cell, const Rect * const Leaf);
    bool
    CheckEveryLeafForTriangles()
    {
      bool OK = false;
      for (std::vector<Rect>::iterator ittt = initial_quads.begin();
	  ittt != initial_quads.end(); ittt++) {
	for (std::list<TriCell*>::iterator itTC = (*ittt).TrisInQuad.begin();
	    itTC != (*ittt).TrisInQuad.end(); itTC++) {
	  OK = QContainsCell(*itTC, &(*ittt));
	  if (!OK)
	    return OK;
	}
      }
      return OK;
    }
    //TODO change the vector to some other stl container that does not the reallocation.
  };

  struct Rect3D {
    double x, y, z, length, width, height;
    std::list<TetCell*> TetsInOct;
    bool IsLocked;
    bool Wmother;
    bool Wchildren;
    omp_lock_t Llck;
    int LeafLayer;
    int Identifier;
    int SizeOrder;
    Rect3D * Mom;
    Rect3D * Children[8];
    std::set<Rect3D*> BufferZone;
    std::vector<BFace*> m_newBdryFaces;
    Rect3D() :
	x(0), y(0), z(0), length(0), width(0), height(0), TetsInOct(), IsLocked(
	    false), Wmother(false), Wchildren(false), Llck()
    {
      Mom = NULL;
      for (int i = 0; i < 8; i++)
	Children[i] = NULL; // i = (Z*Ny*Nx+(Y*Nx+X))
    }
    Rect3D(double XYZ[], double LWH[], long int id, int layer = 1) :
	x(XYZ[0]), y(XYZ[1]), z(XYZ[2]), length(LWH[0]), width(LWH[1]), height(
	    LWH[2]), TetsInOct(), IsLocked(false), Wmother(false), Wchildren(
	    false), Llck(), LeafLayer(layer), Identifier(id), SizeOrder(0)
    {
      Mom = NULL;
      for (int i = 0; i < 8; i++)
	Children[i] = NULL; // i = (Z*Ny*Nx+(Y*Nx+X))
    }
    inline void
    Push_Back_Cell(TetCell * tet)
    {
#pragma omp critical(PushBackCells)
      {
	this->TetsInOct.push_back(tet);
      }
    }
    void
    UpdateCellLists(TetCell* pTC)
    {
//		omp_set_lock(&Llck);
      TetsInOct.push_back(pTC);
//		omp_unset_lock(&Llck);
    }
  };

  class OctTree {
  private:
    VolMesh* Mesh3D;
    double InitialSideSize;
    bool POAlter;
    double MaxCircumradius;
    double Xbound[2];
    double Ybound[2];
    double Zbound[2];
    /// NumRects[0] is the number of rects in x direction.
    /// NumRects[1] is the number of rects in y direction.
    /// NumRects[2] is the number of rects in z direction.
    int NumRects3D[3];
//	std::set<Rect3D*> WholeBufferLeaves;
    std::vector<std::vector<Rect3D*> > KTInLeaves;
  public:
    std::vector<Rect3D> OctLeafList;
    int NOctLayer;
    OctTree(VolMesh* meshToOct) :
	Mesh3D(meshToOct), InitialSideSize(0), POAlter(false), MaxCircumradius(
	    0), OctLeafList(), NOctLayer(4)
    {
      Xbound[0] = 0;
      Xbound[1] = 0;
      Ybound[0] = 0;
      Ybound[1] = 0;
      Zbound[0] = 0;
      Zbound[1] = 0;
      NumRects3D[0] = 0;
      NumRects3D[1] = 0;
      NumRects3D[3] = 0;
    }
    struct OrderedQueue {
//		OrderedQueue& OrderedQueue::operator=(const OrderedQueue& X){
//			Queue = X.Queue;
//			return *this;
//		}
      std::map<int, Rect3D*> Queue;
      std::multimap<int, Rect3D*> MMQueue;
      void
      Insert(Rect3D* leaf)
      {
	int dkey = leaf->TetsInOct.size();
	std::pair<int, Rect3D*> input = std::make_pair(dkey, leaf);
	MMQueue.insert(input);
      }
      void
      Update()
      {
	std::multimap<int, Rect3D*>::iterator itS = MMQueue.begin(), itE =
	    MMQueue.end();
	Queue.clear();
	int num = 0;
	for (; itS != itE; itS++) {
	  num++;
	  itS->second->SizeOrder = num;
	  std::pair<int, Rect3D*> input = std::make_pair(num, itS->second);
	  Queue.insert(input);
	}
	MMQueue.clear();
      }
//		void Delete(Rect* leaf){
////			Queue.remove(leaf);
//		}
      void
      Clear()
      {
	Queue.clear();
	MMQueue.clear();
      }
    };
    OrderedQueue WholeBufferLeaves;
    inline void
    SetNOctLayer(int NL)
    {
      NOctLayer = NL;
    }
    bool
    CalcBoundingBox();
    inline void
    UpdateMaxCirR(double CircR)
    {
      MaxCircumradius = CircR;
    }
    inline void
    UpdateBBound(double XB[], double YB[], double ZB[])
    {
      if (XB[0] < Xbound[0])
	Xbound[0] = XB[0];
      if (XB[1] > Xbound[1])
	Xbound[1] = XB[1];
      if (YB[0] < Ybound[0])
	Ybound[0] = YB[0];
      if (YB[1] > Ybound[1])
	Ybound[1] = YB[1];
      if (ZB[0] < Zbound[0])
	Zbound[0] = ZB[0];
      if (ZB[1] > Zbound[1])
	Zbound[1] = ZB[1];
    }
    int
    Initial_OctTree_Creator(bool BoundIsGiven = true, int NOLayers = 5);
    int
    OctTreeCreatorSecondPasses(int FLevel);
    Rect3D*
    WhichRect(const TetCell* tet);
    void
    AddTetToRect(TetCell* tet);
    bool
    QContainsCell(const TetCell * const Cell, const Rect3D * const Leaf);
    void
    CalcLeafIndex(const Rect3D* leaf, int &i, int &j, int &k);
    int
    FilterTet(const TetCell* tet, const Rect3D* leaf);
    int
    FindIndSets(std::set<Rect3D*>& IndependentSet, bool IsNewLayer = true);
    int
    CalcBuffZone(Rect3D * const leaf);
    int
    NorthNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& NNeihbor);
    int
    SouthNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& SNeihbor);
    int
    WestNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& WNeihbor);
    int
    EastNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& ENeihbor);
    int
    UpNeighbor(const Rect3D* leaf, std::vector<Rect3D*>& UNeihbor);
    int
    DownNeighbor(const Rect3D* Leaf, std::vector<Rect3D*>& DNeihbor);
    bool
    IsOctinBufferZone(Rect3D * const leaf, Rect3D * const LeafB);
    inline void
    clearKTInLeaves()
    {
      KTInLeaves.clear();
    }
    inline bool
    IsWBLEmpty()
    {
      return WholeBufferLeaves.Queue.empty();
    }
    inline void
    FillWBZ(bool NewOct = false)
    {
      if (!NewOct)
	WholeBufferLeaves.Queue.insert(
	    std::pair<int, Rect3D*>(1, &(*OctLeafList.begin())));
      WholeBufferLeaves.MMQueue.clear();
      std::vector<Rect3D>::iterator itW = OctLeafList.begin();
      for (; itW != OctLeafList.end(); itW++)
	WholeBufferLeaves.Insert(&(*itW));
      WholeBufferLeaves.Update();
    }

//	void init_all_locks()
//	{
//		for (GR_index_t i =0 ; i<OctLeafList.size() ; i++)
//			omp_init_lock(&(OctLeafList[i].Llck));
//	}
//	void dest_all_locks()
//	{
//		for (GR_index_t i =0 ; i<OctLeafList.size() ; i++)
//			omp_destroy_lock(&(OctLeafList[i].Llck));
//	}
    // ***** Debugging functions *********
    void
    writeVTKOcttree();
    bool
    CheckEveryLeafForTets();

  };
}	//End of namespace

#endif /* GR_TREECONSTRUCTOR_H_ */
