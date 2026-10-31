#include "GR_CellQueue.h"
#include "GR_Aniso.h"

CellQueue::CellQueue(Mesh* pM_In, double dMaxBig_ = LARGE_DBL,
		     double dMinSmall_ = LARGE_DBL, bool qAttachToMesh = true) :
    pM(pM_In), setpCBig(), setpFSmall(), dMaxBig(dMaxBig_), dMinSmall(
	dMinSmall_)
{

  if (qAttachToMesh)
    pM->attachCellQueue(this);

  printf(
      "Constructing Cell Queue containing all cells with CircumRadius>%.3e or edges smaller than %.3e\n",
      dMaxBig, dMinSmall);
  for (GR_index_t iC = 0; iC < pM->getNumCells(); iC++) {
    vAddCell(pM->getCell(iC));
  }
  printf(
      "Cell Queue constructed with %zd Big cells and %zd small faces out of %d cells, %d faces\n",
      iBig(), iSmall(), pM->getNumCells(), pM->getNumFaces());
}

CellQueue::~CellQueue()
{
  pM->attachCellQueue(NULL);
}

void
CellQueue::vReset()
{
  setpCBig.clear();
  setpFSmall.clear();

  for (GR_index_t iC = 0; iC < pM->getNumCells(); iC++) {
    vAddCell(pM->getCell(iC));
  }

}

bool
CellQueue::qNeverAdd(const Cell* pC)
{
  assert(pC->isValid());
  //  bool qRetVal;
  //a function to put limits on what cells can be added (ie nothing smaller than... or nothing in a particular spot etc...
  return false;
}

void
CellQueue::vAddCell(Cell* pC)
{
  //  vTest();
  int i;
  double dLength;
  assert(pC->isValid());
  assert(!pC->isDeleted());
  if (qNeverAdd(pC))
    return;

  //Attach values to cells and faces
  double dBiggness = (pM->pAnisoRef->dAnisoCircumRadius(pC));

  pC->setAnisoCircumradius(dBiggness);
  for (i = 0; i < pC->getNumFaces(); i++) {
    dLength = (pM->pAnisoRef->dMetricLength(pC->getFace(i)));
    pC->getFace(i)->setLength(dLength);
  }
  // assert(!pC->qIsBdryCell());
  //  vTest();
  if (dBiggness > dMaxBig) {
    setpCBig.insert(pC);
  }

  bool qDoNotAdd;
  for (i = 0; i < pC->getNumFaces(); i++) {
    Face* pF = pC->getFace(i);
    qDoNotAdd = pF->qDoNotDelete();
    dLength = pF->getLength();
    assert(!pF->isDeleted());
    if (dLength < dMinSmall && !qDoNotAdd)
      setpFSmall.insert(pF);
  }

  return;
}

bool
CellQueue::qHasCell(Cell *pC)
{
  bool qRetVal;
  qRetVal = (setpCBig.count(pC) != 0);
  return qRetVal;
}
void
CellQueue::vRemoveFace(Face *pF)
{
  if (setpFSmall.count(pF) != 0) {
    setpFSmall.erase(pF);
  }
}

void
CellQueue::vRemoveCell(Cell *pC)
{
  assert(pC->isValid());

  if (pC->isBdryCell())
    assert(0);

  //  vTest();
  //  printf("Testing for presence of Cell -%X\n",pC);
  if (setpCBig.count(pC) != 0) {
    //    printf("Big Cell Found...");
    setpCBig.erase(pC);
    //   printf("Erased\n");
  }
  for (int i = 0; i < pC->getNumFaces(); i++) {
    Face* pF = pC->getFace(i);
    if (setpFSmall.count(pF) != 0) {
      //  printf("Small Face Found...");
      setpFSmall.erase(pF);
      //     printf("Erased\n");
    }
  }

  return;
}

Cell*
CellQueue::pCTopCell()
{

  Cell* pC = NULL;

  if (!setpCBig.empty()) {
    pC = *(setpCBig.begin());
    assert(pC->isValid());
    assert(!pC->isDeleted());
    //  assert(qHasCell(pC));
  }

  //  printf("pCTopCell returning %X\n",pC);
  // assert(!pC->qDeleted());
  return pC;
}
Face*
CellQueue::pFSmallestFace()
{
  Face* pF = NULL;
  if (!setpFSmall.empty()) {
    pF = *(setpFSmall.begin());
    if (pF->isDeleted()) {   //This shouldn't be required, but it is
      setpFSmall.erase(pF);
      pF = pFSmallestFace();
      // assert(0);
    }
  }
//  assert(!pF->isDeleted());
  return pF;
}

void
CellQueue::vTest()
{


  printf("Checking %zd,%zd\n", setpCBig.size(), setpFSmall.size());

  for (std::set<Cell*>::iterator iterC = setpCBig.begin();
      iterC != setpCBig.end(); ++iterC) {
    assert(!(*iterC)->isDeleted());
  }
  printf("Big okay\n");
  for (std::set<Face*>::iterator iterF = setpFSmall.begin();
      iterF != setpFSmall.end(); ++iterF) {
    assert(!(*iterF)->isDeleted());
  }
  printf("Small okay\n");
  return;
}
