#ifndef GR_CellQueue
#define GR_CellQueue 47

#include <map>

#include <GR_Mesh.h>
#include <GR_Vertex.h>
#include <GR_Cell.h>
#include <GR_Face.h>
#include <GR_QualMeasure.h>

//Compares cells for insertion priority
struct biggerthenCell {
  bool
  operator()(const Cell* pC1, const Cell* pC2) const
  {
    if (pC1 == pC2)
      return false;    //So that equality is dependant only on the pointer

    double dPri1, dPri2;
    dPri1 = pC1->getAnisoCircumradius();
    dPri2 = pC2->getAnisoCircumradius();

    if (iFuzzyComp(dPri1, dPri2) == 0) {
      //If the cells have the same priority, they are sorted by pointer.
      return (pC1 > pC2);
    }
    return (dPri1 > dPri2);
  }
};

struct smallerthenFace {
  bool
  operator()(const Face* pF1, const Face* pF2) const
  {
    if (pF1 == pF2)
      return false;    //So that equality is dependant only on the pointer
    double dPri1, dPri2;
    //  dPri1=pC1->pVVert(0)->dX();
    //   dPri2=pC2->pVVert(0)->dX();

    dPri1 = -(pF1->getLength());
    dPri2 = -(pF2->getLength());

    if (iFuzzyComp(dPri1, dPri2) == 0) {
      //If the cells have the same priority, they are sorted by pointer.
      return (pF1 > pF2);
    }
    return (dPri1 > dPri2);
  }
};

//struct InsertionCriteria
//{
//  bool operator()(const Face* pF1,const Face* pF2) const
//  {
//	if(pF1==pF2) return false;
//
//	double dPF1,dPF2;
//    dPF1 = (pF1->getdQualInc());
//    dPF2 = (pF2->getdQualInc());
//
//    return (dPF1<dPF2);
//  }
//};

class CellQueue {
  CellQueue(const CellQueue&) :
      pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
  {
    assert(0);
  }
  CellQueue&
  operator=(const CellQueue&)
  {
    assert(0);
    return (*this);
  }
public:
  CellQueue() :
      pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
  {
    assert(0);
  }

  CellQueue(Mesh* pM_In, double dMinSize_, double dMaxAngle_,
	    bool qAttachToMesh);

  ~CellQueue(); //need to detach cell queue (TODO)

  Mesh* pM;
  std::set<Cell*, biggerthenCell> setpCBig;
  std::set<Face*, smallerthenFace> setpFSmall;

  double dMaxBig;
  double dMinSmall;

  size_t
  iBig()
  {
    return setpCBig.size();
  }
  size_t
  iSmall()
  {
    return setpFSmall.size();
  }
  size_t
  iSize()
  {
    return iBig() + iSmall();
  }

  Cell*
  pCTopCell();
  Face*
  pFSmallestFace();

  void
  vReset();
  void
  vRemoveCell(Cell* pC);
  void
  vRemoveFace(Face* pF);
  void
  vAddCell(Cell* pC);
  bool
  qHasCell(Cell* pC);
  bool
  qNeverAdd(const Cell* const pC);
  void
  vTest();
};

/*
 * To do i queue:
 * 1.- check for quality in edge after splitting
 * 2.- if improves mark and save how much it improves
 * 3.- reorder from best to worst
 * 4.- if 2 faces are within the same cell, delete the one with the worst improvement
 * 5.- reevaluate the faces from the cell with the deleted face
 */

//
//class InsertFaceQueue {
//  InsertFaceQueue(const InsertFaceQueue&) :
//    pM(NULL), setpFSmall(), dMaxBig(), dMinSmall()
//  {assert(0);}
//  InsertFaceQueue& operator=(const InsertFaceQueue&) {assert(0); return(*this);}
// public:
//  //InsertFaceQueue(): pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
//  //{assert(0);}
//
//  InsertFaceQueue(Mesh* pM2D, GRUMMP::TMOPQual2D* pTQ);
//
//  ~InsertFaceQueue();
//
//  Mesh* pM;
//  std::set<Face*,smallerthenFace> setpFSmall;
//
//  double  dMaxBig;
//  double  dMinSmall;
//
//  size_t iSmall(){return setpFSmall.size();}
//  size_t iSize(){return iSmall();}
//
//  Cell* pCTopCell();
//  Face* pFSmallestFace();
//
//  void vReset();
//  void vRemoveCell( Cell* pC);
//  void vRemoveFace(Face* pF);
//  void vAddCell(Cell* pC);
//  bool qHasCell(Cell* pC);
//  bool qNeverAdd(const Cell* const pC);
//  void vTest();
//};
//
//
///*
// * Same process but extension to 3D:
// * 1.- run 1st for faces
// * 2.- run second for edges
// */
//
//class InsertFaceQueue3D {
//  InsertFaceQueue3D(const InsertFaceQueue3D&) :
//    pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
//  {assert(0);}
//  InsertFaceQueue3D& operator=(const InsertFaceQueue3D&) {assert(0); return(*this);}
// public:
//  InsertFaceQueue3D(): pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
//  {assert(0);}
//
//  InsertFaceQueue3D(Mesh* pM_In,double dMinSize_,double dMaxAngle_,bool qAttachToMesh);
//
//  ~InsertFaceQueue3D(); //need to detach cell queue (TODO)
//
//  Mesh* pM;
//  std::set<Cell*,biggerthenCell> setpCBig;
//  std::set<Face*,smallerthenFace> setpFSmall;
//
//  double  dMaxBig;
//  double  dMinSmall;
//
//  size_t iBig(){return setpCBig.size();}
//  size_t iSmall(){return setpFSmall.size();}
//  size_t iSize(){return iBig()+iSmall();}
//
//  Cell* pCTopCell();
//  Face* pFSmallestFace();
//
//  void vReset();
//  void vRemoveCell( Cell* pC);
//  void vRemoveFace(Face* pF);
//  void vAddCell(Cell* pC);
//  bool qHasCell(Cell* pC);
//  bool qNeverAdd(const Cell* const pC);
//  void vTest();
//};
//
//class InsertEdgeQueue3D {
//  InsertEdgeQueue3D(const InsertEdgeQueue3D&) :
//    pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
//  {assert(0);}
//  InsertEdgeQueue3D& operator=(const InsertEdgeQueue3D&) {assert(0); return(*this);}
// public:
//  InsertEdgeQueue3D(): pM(NULL), setpCBig(), setpFSmall(), dMaxBig(), dMinSmall()
//  {assert(0);}
//
//  InsertEdgeQueue3D(Mesh* pM_In,double dMinSize_,double dMaxAngle_,bool qAttachToMesh);
//
//  ~InsertEdgeQueue3D(); //need to detach cell queue (TODO)
//
//  Mesh* pM;
//  std::set<Cell*,biggerthenCell> setpCBig;
//  std::set<Face*,smallerthenFace> setpFSmall;
//
//  double  dMaxBig;
//  double  dMinSmall;
//
//  size_t iBig(){return setpCBig.size();}
//  size_t iSmall(){return setpFSmall.size();}
//  size_t iSize(){return iBig()+iSmall();}
//
//  Cell* pCTopCell();
//  Face* pFSmallestFace();
//
//  void vReset();
//  void vRemoveCell( Cell* pC);
//  void vRemoveFace(Face* pF);
//  void vAddCell(Cell* pC);
//  bool qHasCell(Cell* pC);
//  bool qNeverAdd(const Cell* const pC);
//  void vTest();
//};
//
//
//class InsertionQueue2 {
//	InsertionQueue2(const InsertionQueue2&) :
//	    pM(NULL), dMaxBig(), dMinSmall()
//	  {assert(0);}
//	InsertionQueue2& operator=(const InsertionQueue2&) {assert(0); return(*this);}
//	 public:
//	  InsertionQueue2(): pM(NULL), dMaxBig(), dMinSmall()
//	  {assert(0);}
//
//	  InsertionQueue2(Mesh* pMesh);
//	  virtual~InsertionQueue2();
//
//	  Mesh* pM;
//	  std::set<Face*,InsertionCriteria > QualPriority;
//	  std::set<EdgeFace*,InsertionCriteria > EdgeQualPriority;
//
//	  double  dMaxBig;
//	  double  dMinSmall;
//
//	  virtual void ResetQueue();
//	  virtual void RemoveFromQueue();
//	  virtual void AddToQueue();
//	  virtual bool HasQueue();
//	  virtual void GetTopQueue();
//	  virtual void SizeQueue();
//
//};
#endif 
