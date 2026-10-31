#include <algorithm>
#include <sstream>
#include <vector>
#include <cstdio>

using std::vector;

#include "GR_config.h"
#include "GR_ADT.h"
#include "GR_Length.h"
#include "GR_misc.h"
#include "GR_InsertionManager.h"
#include "GR_RefinementManager3D.h"
#include "GR_SwapManager.h"
#include "GR_TriMeshBuilder.h"
#include "GR_VolMesh.h"

#include "edamVert.h"
#include "./edamMesh.h"

struct clipV {
  double coords[3];
  edamVert* pEV;

  clipV(double x, double y, double z) :
      pEV(NULL)
  {
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
  }
  clipV(edamVert* pEVArg) :
      pEV(pEVArg)
  {
  }
};

struct clipE {
  clipV v0;
  clipV v1;

  clipE(clipV v0Arg, clipV v1Arg) :
      v0(v0Arg), v1(v1Arg)
  {
  }
  clipE(edamVert* pEV0, edamVert* pEV1) :
      v0(pEV0), v1(pEV1)
  {
  }
};

class edamFace {
  Face* pF;
  edamVert* m_verts[4];
  vector<edamVert*> kids;
  int m_nVerts;

  edamFace(); //never defined
  edamFace&
  operator=(const edamFace&); // never defined
public:
  edamFace(Face* pFArg, edamVert* pEV0, edamVert* pEV1, edamVert* pEV2,
	   edamVert* pEV3 = NULL) :
      pF(pFArg), kids(3), m_nVerts(3)
  {
    m_verts[0] = pEV0;
    m_verts[1] = pEV1;
    m_verts[2] = pEV2;
    if (pEV3) {
      m_nVerts = 4;
      m_verts[3] = pEV3;
    }
  }

  Face*
  getFace() const
  {
    return pF;
  }
  int
  getNumVerts() const
  {
    return m_nVerts;
  }
  edamVert*
  getVert(int iVert) const
  {
    return m_verts[iVert];
  }
  edamVert*
  getKid(int iVert) const
  {
    return kids[iVert];
  }
  edamVert*
  getKid(edamVert* pEV) const
  {
    for (int iV = 0; iV < getNumVerts(); iV++) {
      if (getVert(iV) == pEV) {
	return getKid(iV);
      }
    }
    //assert(0);
    return NULL;
  }
  void
  setKid(edamVert* pEVpar, edamVert* pEVkid)
  {
    bool qGood = false;
    for (int iV = 0; iV < getNumVerts(); iV++) {
      if (getVert(iV) == pEVpar) {
	qGood = true;
	kids[iV] = pEVkid;
      }
    }
    assert(qGood);
  }
  edamVert*
  getNextVert(const edamVert* pEV) const
  {
    int nv = getNumVerts();
    for (int iV = 0; iV < nv; iV++) {
      if (getVert((iV - 1 + nv) % nv) == pEV) {
	return getVert(iV);
      }
    }
    assert(0);
    return NULL;
  }
  edamVert*
  getPrevVert(const edamVert* pEV)
  {
    int nv = getNumVerts();
    for (int iV = 0; iV < nv; iV++) {
      if (getVert((iV + 1) % nv) == pEV) {
	return getVert(iV);
      }
    }
    assert(0);
    return NULL;
  }
  bool
  hasVert(edamVert* pEVother)
  {
    for (int ii = getNumVerts() - 1; ii >= 0; ii--) {
      if (m_verts[ii] == pEVother)
	return true;
    }
    return false;
  }
  void
  interchangeVerts()
  {
    assert(getNumVerts() == 3);
    edamVert* pEVtemp = m_verts[1];
    m_verts[1] = m_verts[2];
    m_verts[2] = pEVtemp;
    pF->interchangeCellsAndVerts();
  }

  int
  getNumCells()
  {
    return pF->getNumCells();
  }
  Cell*
  getLeftCell()
  {
    return pF->getLeftCell();
  }
  Cell*
  getRightCell()
  {
    return pF->getRightCell();
  }
  void
  calcUnitNormal(double adUnitNorm[]) const
  {
    ::calcNormal(m_verts[0]->getVert()->getCoords(),
		 m_verts[1]->getVert()->getCoords(),
		 m_verts[2]->getVert()->getCoords(), adUnitNorm);
    NORMALIZE3D(adUnitNorm);
  }
  double
  calcAngle(edamVert* pEV)
  {
    int iIndex = -1;
    for (int i = 0; i < 3; i++) {
      if (m_verts[i] == pEV)
	iIndex = i;
    }
    if (iIndex < 0) {
      cout << pEV->getLayer() << ":" << m_verts[0]->getLayer() << ":"
	  << m_verts[1]->getLayer() << ":" << m_verts[2]->getLayer() << endl;
      cout << pEV->x() << ":" << pEV->y() << ":" << pEV->z() << ":" << endl;
      cout << m_verts[0]->x() << ":" << m_verts[0]->y() << ":"
	  << m_verts[0]->z() << endl;
      cout << m_verts[1]->x() << ":" << m_verts[1]->y() << ":"
	  << m_verts[1]->z() << endl;
      cout << m_verts[2]->x() << ":" << m_verts[2]->y() << ":"
	  << m_verts[2]->z() << endl;
      if (m_verts[3])
	cout << m_verts[3]->x() << ":" << m_verts[3]->y() << ":"
	    << m_verts[3]->z() << endl;
    }
    assert(iIndex >= 0);

    int iNext = (iIndex + 1) % 3;
    int iPrev = (iIndex + 3 - 1) % 3;
    const double* adCurr = pEV->getVert()->getCoords();
    const double* adNext = m_verts[iNext]->getVert()->getCoords();
    const double* adPrev = m_verts[iPrev]->getVert()->getCoords();

    double adA[3] = adDIFF3D(adNext,adCurr);
    double adB[3] = adDIFF3D(adPrev,adCurr);
    NORMALIZE3D(adA);
    NORMALIZE3D(adB);

    return acos(dDOT3D(adA, adB));
  }

  static bool
  qPointerCompFaceAddress(const edamFace* pEF0, const edamFace* pEF1)
  {
    return pEF0->pF < pEF1->pF;
  }
};

edamFace*
edamVert::getFace(Face* pF) const
{
  for (edamFace* pEF : faceList) {
    if (pEF->getFace() == pF) {
      return pEF;
    }
  }
  return nullptr;
}
edamFace*
edamVert::getCommonFace(edamVert* pEV0, edamVert* pEV1)
{
  for (edamFace* pEF : faceList) {
    if (pEF->hasVert(pEV0) && pEF->hasVert(pEV1)) {
      return pEF;
    }
  }
  return NULL;
}
set<edamFace*>
edamVert::getCommonFaces(edamVert* pEV0)
{
  set<edamFace*> ret;
  for (edamFace* pEF : faceList) {
    if (pEF->hasVert(pEV0)) {
      ret.insert(pEF);
    }
  }
  return ret;
}
// Get the face connected to pEV0 and "this", but not to pEV1. Used to traverse the mesh locally.
edamFace*
edamVert::getOtherCommonFace(edamVert* pEV0, edamVert* pEV1)
{
  for (edamFace* pEF : faceList) {
    if (pEF->hasVert(pEV0) && !pEF->hasVert(pEV1)) {
      return pEF;
    }
  }
  cout << qIntersectClip << ":" << pEV0->qIntersectClip << ":"
      << pEV1->qIntersectClip << endl;
  cout << x() << ":" << y() << ":" << z() << endl;
  cout << pEV0->x() << ":" << pEV0->y() << ":" << pEV0->z() << endl;
  cout << pEV1->x() << ":" << pEV1->y() << ":" << pEV1->z() << endl;
  assert(0);
  return NULL;
}
bool
edamVert::hasCommonEdge(edamVert* pEV0)
{
  for (edamFace* pEF : faceList) {
    if (pEF->hasVert(pEV0)) {
      return true;
    }
  }
  return false;
}

/*
 * Identify a cell adjacent to the given vert, for use as a guess for where
 * to insert a new vertex.
 */

Cell*
findGuessCell(Vert* const pV)
{
  Face* incFace = pV->getFace(0);
  Cell* pC = incFace->getCell(0);
  if (pC->isBdryCell())
    pC = incFace->getCell(1);
  assert(!pC->isBdryCell());
  assert(pC->hasVert(pV));
  return pC;
}

/*
 * Wrappers for GRUMMP mesh modifiers, use these to keep edamVertList
 * consistent with the mesh.
 */

edamVert*
edamMesh3D_Closed::createEdamVert(double x, double y, double z, int layer)
{
  double coords[] =
    { x, y, z };
  Vert* pV = createVert(coords);
  pV->markStructured();
  m_eVL.push_back(edamVert(pV, m_eVL.size(), layer));
  return m_eVL.getLast();
}

edamFace*
edamMesh3D_Closed::createTriFace(edamVert* pEV0, edamVert* pEV1, edamVert* pEV2)
{
  Face* pF = findCommonFace(pEV0->getVert(), pEV1->getVert(), pEV2->getVert(),
  NULL,
			    true);
  assert(pF);
  assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());
//  bool qExist = false;
//  Face* pF = createFace(qExist, pEV0->getVert(), pEV1->getVert(),
//			pEV2->getVert());
//  assert(qExist);
  // If there is already an edamFace from this mesh face, find out
  // what it is.
  edamFace* pEF = pEV0->getFace(pF);
  if (pEF == nullptr) {
    pEF = new edamFace(pF, pEV0, pEV1, pEV2);
    pEV0->addFace(pEF);
    pEV1->addFace(pEF);
    pEV2->addFace(pEF);
    m_all_edamFaces.push_back(pEF);
  }
  return pEF;
}

void edamMesh3D_Closed::addNewFaces(const vector<vector<edamVert*> >& newTris,
		std::set<edamFace*>* pFaces) {
	for (vector<edamVert*> tri : newTris) {
		edamVert* vert0 = tri[0];
		edamVert* vert1 = tri[1];
		edamVert* vert2 = tri[2];
		edamFace* pEF = createTriFace(vert0, vert1, vert2);
		if (pFaces) pFaces->insert(pEF);
	}
}


edamFace*
createQuadFace(VolMesh* pM3D, edamVert* pEV0, edamVert* pEV1, edamVert* pEV2,
	       edamVert* pEV3)
{
  bool qExist;

  Face* pF = pM3D->createFace(qExist, pEV0->getVert(), pEV1->getVert(),
			      pEV2->getVert(), pEV3->getVert());
  edamFace* pEF = new edamFace(pF, pEV0, pEV1, pEV2, pEV3);
  if (!qExist) {
    pEV0->addFace(pEF);
    pEV1->addFace(pEF);
    pEV2->addFace(pEF);
    pEV3->addFace(pEF);
  }

  return pEF;
}
void
edamMesh3D_Closed::deleteTriFace(edamFace* pEF)
{
  assert(pEF->getNumCells() == 0);
  for (int iV = 0; iV < 3; iV++) {
    edamVert* pEV = pEF->getVert(iV);
    pEV->removeFace(pEF);
  }
  deleteFace(pEF->getFace());
}
void
edamMesh3D_Closed::deleteQuadFace(edamFace* pEF)
{
  assert(pEF->getNumCells() == 0);
  for (int iV = 0; iV < 4; iV++) {
    edamVert* pEV = pEF->getVert(iV);
    pEV->removeFace(pEF);
  }
  deleteFace(pEF->getFace());
}

void
edamMesh3D_Closed::deleteEdamVert(edamVert* pEV)
{
  for (edamFace* pEF : pEV->getFaceList()) {
    if (pEF->getNumVerts() == 3) {
      assert(0);
      deleteTriFace(pEF);
    }
    else {
      //deleteQuadFace(pM3D,m_eVL,pEF);
    }
  }
  edamVert *parent = pEV->getParent();
  if (parent) {
    for (edamFace* pEF : parent->getAdjFaces()) {
      pEF->setKid(parent, NULL);
    }
  }
  deleteVert(pEV->getVert());
  m_eVL.erase(pEV->getIndex());
  pEV->clearVert();
}

edamMesh3D_Closed::edamMesh3D_Closed(const double extrudeLen,
				     const double growthRatio,
				     const int targetNumLayers,
				     const eMarchDir eMD) :
    VolMesh(), m_eVL(), m_firstLayerHeight(extrudeLen), m_growthRatio(
	growthRatio), m_targetLayers(targetNumLayers), m_eMD(eMD)
{
}

edamMesh3D_Closed::~edamMesh3D_Closed()
{
  for (edamFace* EF : m_all_edamFaces) {
    delete EF;
  }
}

void
edamMesh3D_Closed::initFront(std::set<int> const& BCsToMarchFrom)
{
  std::map<int, int> indexMap;

  int iIndex = 0;
  for (size_t iV = 0; iV < getNumVerts(); iV++) {
//		assert(m_eVL.size() == iV);
//		m_eVL.push_back(terribleVert(VM.getVert(iV), iV, 0));

    Vert* pV = getVert(iV);
    if (pV->isBdryVert() || pV->getVertType() == Vert::ePseudoSurface) {
      m_eVL.push_back(edamVert(pV, iIndex, 0));
      indexMap.insert(pair<int, int>(iV, iIndex));
      iIndex++;
    }
  }

  for (int iBF = getNumBdryFaces() - 1; iBF >= 0; iBF--) {
    BFace* pBF = getBFace(iBF);
    if (pBF->isDeleted()) continue;
    int thisBC = pBF->getBdryCondition();
    if (BCsToMarchFrom.find(thisBC) == BCsToMarchFrom.end()) {
      continue;
    }
    assert(pBF->getNumVerts() == 3);

    edamVert* pEV0 = m_eVL.get(indexMap.at(getVertIndex(pBF->getVert(0))));
    edamVert* pEV1 = m_eVL.get(indexMap.at(getVertIndex(pBF->getVert(1))));
    edamVert* pEV2 = m_eVL.get(indexMap.at(getVertIndex(pBF->getVert(2))));

    switch (m_eMD)
      {
      case eOut:
	createTriFace(pEV0, pEV2, pEV1);
	break;
      case eIn:
	createTriFace(pEV0, pEV1, pEV2);
	break;
      case eBoth:
	createTriFace(pEV0, pEV1, pEV2);
	createTriFace(pEV0, pEV2, pEV1);
	break;
      }
  }
}

/*
 * Run GRUMMP's "tetra" to generate a triangular mesh of the geometry.
 */
void
meshFileFromBdryFile(const char* bdryFile, double dScale, double dGrading,
		     const char* lengthFile = NULL)
{

  char command[FILE_NAME_LEN];
  if (lengthFile) {
    sprintf(
	command,
	"/usr/anslab/src/GRUMMP/bin/linux-gnu/tetra -i %s -r %f -g %f -l %s",
	bdryFile, dScale, dGrading, lengthFile);
  }
  else {
    sprintf(command,
	    "/usr/anslab/src/GRUMMP/bin/linux-gnu/tetra -i %s -r %f -g %f",
	    bdryFile, dScale, dGrading);
  }

  if (system(command) != 0)
    assert(0);
}

///*
// * Reads a .obj file to obtain the initial surface.
// */
//void
//edamMesh3D_Closed::initFromObjFile(const char* filename,
//			    const std::set<int>& /*BCsToMarchFrom*/)
//{
//  // Marching out from all faces.
//  ifstream fObj;
//  char infilename[256];
//  sprintf(infilename, "%s.obj", filename);
//  fObj.open(infilename);
//
//  vector<vector<double>> verts;
//  vector<vector<int>> faces;
//
//  while (!fObj.eof()) {
//    string strLine;
//
//    std::getline(fObj, strLine);
//    istringstream iss(strLine);
//
//    string strType;
//    iss >> strType;
//    if (strType == "v") {
//      vector<double> vert(3);
//      iss >> vert[0] >> vert[1] >> vert[2];
//      verts.push_back(vert);
//    }
//    else if (strType == "f") {
//      vector<int> face(3);
//      iss >> face[0] >> face[1] >> face[2];
//      faces.push_back(face);
//      if (!iss.eof() && 0) {
//	vector<int> face2(3);
//	iss >> face2[1];
//	face2[0] = face[2];
//	face2[2] = face[0];
//	faces.push_back(face2);
//      }
//    }
//  }
//
//  for (size_t iV = 0; iV < verts.size(); iV++) {
//    vector<double> vert = verts.at(iV);
//    createEdamVert(vert[0], vert[1], vert[2], 0);
//  }
//
//  for (size_t iV = 0; iV < m_eVL.size(); iV++) {
//    assert(m_eVL.get(iV)->getIndex() == static_cast<int>(iV));
//    assert(getVertIndex(m_eVL.get(iV)->getVert()) == iV);
//  }
//
//  for (size_t iF = 0; iF < faces.size(); iF++) {
//    vector<int> face = faces.at(iF);
//
//    edamVert* pEV0 = m_eVL.get(face[0] - 1);
//    edamVert* pEV1 = m_eVL.get(face[1] - 1);
//    edamVert* pEV2 = m_eVL.get(face[2] - 1);
//
//    createTriFace(pEV0, pEV1, pEV2);
//  }
//}

// Assumes the face list comprises a closed manifold.
/*
 * Write a list of faces to a .off mesh file.
 *
 * Used to generate an input file for the ReebHanTun software which locates tunnels.
 */
void
OFFFaceList(vector<edamFace*> &faces, const char* fileName,
	    vector<edamVert*> &vertList, map<int, int> &tempToReal)
{

  set<edamVert*> verts;
  for (edamFace* pEF : faces) {
    for (int iV = 0; iV < 3; iV++) {
    	edamVert* pEV = pEF->getVert(iV);
    	verts.insert(pEV);
    }
  }

  vertList.assign(verts.begin(), verts.end());
  map<int, int> realToTemp;

  ofstream fOFF;
  fOFF.open(fileName);
  fOFF.precision(16);

  fOFF << "OFF" << endl;
  fOFF << verts.size() << " " << faces.size() << " 0" << endl;
  for (size_t iV = 0; iV < vertList.size(); iV++) {
    edamVert* pEV = vertList.at(iV);
    fOFF << pEV->x() << " " << pEV->y() << " " << pEV->z() << endl;

    realToTemp.insert(pair<int, int>(pEV->getIndex(), iV));
    tempToReal.insert(pair<int, int>(iV, pEV->getIndex()));
  }
  for (edamFace* pEF : faces) {
    fOFF << 3;
    for (int iV = 0; iV < 3; iV++) {
      fOFF << " " << realToTemp.at(pEF->getVert(iV)->getIndex());
    }
    fOFF << endl;
  }
}

/*
 * Write a list of faces (triangular, haven't had reason to visualize quad faces yet) to a VTK file.
 */
void
writeFaceListToVTK(const set<edamFace*> &faces, const char* fileName)
{

  set<edamVert*> verts;
  for (edamFace* pEF : faces) {
    for (int iV = 0; iV < 3; iV++) {
    	edamVert* pEV = pEF->getVert(iV);
    	verts.insert(pEV);
    }
  }

  vector<edamVert*> vertList(verts.begin(), verts.end());
  map<int, int> vertIndexMap;

  ofstream fVTK;
  fVTK.open(fileName);
  fVTK.precision(16);

  fVTK << "# vtk DataFile Version 1.0" << endl;
  fVTK << "OVMesh Tri example" << endl;
  fVTK << "ASCII" << endl;
  fVTK << endl;
  fVTK << "DATASET UNSTRUCTURED_GRID" << endl;

  fVTK << "POINTS " << verts.size() << " float" << endl;
  for (size_t iV = 0; iV < vertList.size(); iV++) {
    edamVert* pEV = vertList.at(iV);
    fVTK << pEV->x() << " " << pEV->y() << " " << pEV->z() << endl;

    vertIndexMap.insert(pair<int, int>(pEV->getIndex(), iV));
  }
  fVTK << endl;

  fVTK << "CELLS " << faces.size() << " " << faces.size() * 4 << endl;
  for (edamFace* pEF : faces) {
    fVTK << 3;
    for (int iV = 0; iV < 3; iV++) {
      fVTK << " " << vertIndexMap.at(pEF->getVert(iV)->getIndex());
    }
    fVTK << endl;
  }
  fVTK << endl;

  fVTK << "CELL_TYPES " << faces.size() << endl;
  for (size_t ii = 0; ii < faces.size(); ii++) {
    fVTK << "5" << endl;
  }

  fVTK << "CELL_DATA " << faces.size() << endl;
  fVTK << "SCALARS Cell-Type float 1" << endl;
  fVTK << "LOOKUP_TABLE default" << endl;
  for (edamFace* pEF : faces) {
    bool qClip = false;
    bool qClip2 = false;
    bool qBad = false;
    bool qGood = false;
    for (int ii = pEF->getNumVerts() - 1; ii >= 0; ii--) {
      edamVert* pEV = pEF->getVert(ii);
      if (pEV->isIntersectClip()) {
	qClip = true;
      }
      else if (pEV->isIntersectClip2()) {
	qClip2 = true;
      }
    }
    int iType = 0;
    if (qClip) {
      iType = 1;
    }
    else if (qClip2) {
      iType = 4;
    }
    else if (qBad && qGood) {
      iType = -1;
    }
    else if (qGood) {
      iType = 3;
    }
    else if (qBad) {
      iType = 2;
    }

//		int iType = -1;
//		for (terribleVert* pEV : pEF->verts) {
//			if(iType == -1){
//				iType = pEV->iLoop;
//			}
//			else if(iType != pEV->iLoop){
//				iType = -2;
//			}
//		}
    fVTK << iType << endl;
  }
}

set<edamFace*>
createFrontFaceList(int iLayer, edamVertList& m_eVL)
{
  set<edamFace*> topFaceSet;
  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
      iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner() || pEV->getNumFaces() == 0)
      continue;

    set<edamFace*> thisFaceSet = pEV->getAdjFaces();
    for (edamFace* pEF : thisFaceSet) {
      if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer + 1
	  && pEF->getVert(1)->getLayer() == iLayer + 1
	  && pEF->getVert(2)->getLayer() == iLayer + 1) {
	topFaceSet.insert(pEF);
      }
    }
  }
  return topFaceSet;
}

/*
 * Write the front faces to a VTK file. Generally a useful debugging tool.
 */
void
writeFrontToVTK(edamVertList &m_eVL, int iLayer, const char* fileName)
{
  set<edamFace*> topFaceSet = createFrontFaceList(iLayer, m_eVL);
  writeFaceListToVTK(topFaceSet, fileName);
}

/*
 * Creates a .bdry file from a list of faces (they must define a closed cavity).
 *
 * Not used anywhere in this code, but leaving it here in case it's ever useful for something.
 */
void
bdryFileFromFaceList(vector<edamFace*> &faces, map<int, int> &realToTemp)
{
  set<edamVert*> verts;
  for (edamFace* pEF : faces) {
    for (int iV = 0; iV < 3; iV++) {
      verts.insert(pEF->getVert(iV));
    }
  }

  vector<edamVert*> vertList(verts.begin(), verts.end());
  map<int, int> vertIndexMap;

  ofstream fBdry;
  fBdry.open("tempbdry.bdry");
  fBdry.precision(16);

  fBdry << verts.size() << " " << faces.size() << endl;
  for (size_t iV = 0; iV < vertList.size(); iV++) {
    edamVert* pEV = vertList.at(iV);
    fBdry << pEV->x() << " " << pEV->y() << " " << pEV->z() << endl;

    vertIndexMap.insert(pair<int, int>(pEV->getIndex(), iV));
  }

  for (edamFace* pEF : faces) {
    fBdry << "polygon b 1 r 1 3";
    for (int iV = 0; iV < 3; iV++) {
      fBdry << " " << vertIndexMap.at(pEF->getVert(iV)->getIndex());
    }
    fBdry << endl;
  }

  realToTemp = vertIndexMap;
}

//void readObj(const char* filename){
//	ifstream fObj;
//	char infilename[256];
//	sprintf(infilename,"%s.obj",filename);
//	fObj.open(infilename);
//
//	vector<vector<double>> verts;
//	vector<vector<int>> faces;
//
//	while(!fObj.eof()){
//		string strLine;
//
//		getline(fObj,strLine);
//		istringstream iss(strLine);
//
//		string strType;
//		iss >> strType;
//		if(strType == "v"){
//			vector<double> vert(3);
//			iss >> vert[0] >> vert[1] >> vert[2];
//			verts.push_back(vert);
//		}
//		else if(strType == "f"){
//			vector<int> face(3);
//			iss >> face[0] >> face[1] >> face[2];
//			faces.push_back(face);
//			if(!iss.eof()){
//				vector<int> face2(3);
//				iss >> face2[1];
//				face2[0] = face[2];
//				face2[2] = face[0];
//				faces.push_back(face2);
//			}
//		}
//	}
//
//	ofstream fBdry;
//	char outfilename[256];
//	sprintf(outfilename,"%s.bdry",filename);
//	fBdry.open(outfilename);
//	fBdry.precision(16);
//	SurfMesh* pSM = new SurfMesh();
//
//	fBdry << verts.size() << " " << faces.size() << endl;
//	for(int iV = 0; iV < verts.size(); iV++){
//		vector<double> vert = verts.at(iV);
//		fBdry << vert[0] << " " << vert[1] << " " << vert[2] << endl;
//
//		pSM->createEdamVert(vert[0], vert[1], vert[2]);
//		pSM->getVert(iV)->setType(Vert::eBdryApex);
//	}
//
//	for(int iF = 0; iF < faces.size(); iF++){
//		vector<int> face = faces.at(iF);
//		fBdry << "polygon b 1 r 1 3 " << face[0]-1 << " " << face[1]-1 << " " << face[2]-1 << endl;
//
//		pSM->createTriCell(pSM->getVert(face[0]-1), pSM->getVert(face[1]-1),
//						pSM->getVert(face[2]-1));
//	}
//
//	pSM->createPatches();
//
//	VolMesh* pMtemp = new VolMesh(true, *pSM, 0, 1, 1);
//	writeVTKLegacy(*pMtemp,filename);
//}

/*
 * Helper function for clipIntersectFaceFace. Also tags each vertex as valid or invalid (invalid ones are to be removed).
 */
int
clipIntersectFaceLine(edamFace* pEF, edamVert* pEV0, edamVert* pEV1,
		      double* adPos)
{
  double adNorm[3];
  pEF->calcUnitNormal(adNorm);
  const double* pFace = pEF->getVert(0)->getVert()->getCoords();
  const double* p0 = pEV0->getVert()->getCoords();
  const double* p1 = pEV1->getVert()->getCoords();

  double ad01[3] = adDIFF3D(p1,p0);
  double dDot01 = dDOT3D(ad01, adNorm);
  // false if line is parallel to face (hopeifully...)
  if (fabs(dDot01) < 1e-10)
    return 0;

  // false if line segment and face share a vertex
  if (pEF->hasVert(pEV0) || pEF->hasVert(pEV1))
    return 0;

  // false if line segment does not intersect plane of face
  double adFace0[3] = adDIFF3D(p0,pFace);
  double dDotFace0 = dDOT3D(adFace0, adNorm);
  double dPara = -dDotFace0 / dDot01;
  if (dPara < 1e-10 || dPara > 1 - 1e-10)
    return 0;

  double pInt[3] =
    { 0, 0, 0 };
  double pOffset[3] =
    { 0, 0, 0 };
  for (int ii = 0; ii < 3; ii++) {
    pInt[ii] = p0[ii] + dPara * ad01[ii];
    pOffset[ii] = pFace[ii] + 1e-3 * adNorm[ii];
  }

  const double* pV0 = pEF->getVert(0)->getVert()->getCoords();
  const double* pV1 = pEF->getVert(1)->getVert()->getCoords();
  const double* pV2 = pEF->getVert(2)->getVert()->getCoords();

  int iOrient01 = checkOrient3D(pV0, pV1, pInt, pOffset);
  int iOrient12 = checkOrient3D(pV1, pV2, pInt, pOffset);
  int iOrient20 = checkOrient3D(pV2, pV0, pInt, pOffset);
  if (iOrient01 == -1)
    return 0;
  if (iOrient12 == -1)
    return 0;
  if (iOrient20 == -1)
    return 0;
  if (iOrient01 + iOrient12 + iOrient20 != 3) {
    cout << pV0[0] << ":" << pV0[1] << ":" << pV0[2] << endl;
    cout << pV1[0] << ":" << pV1[1] << ":" << pV1[2] << endl;
    cout << pV2[0] << ":" << pV2[1] << ":" << pV2[2] << endl;
    cout << pEV0->x() << ":" << pEV0->y() << ":" << pEV0->z() << endl;
    cout << pEV1->x() << ":" << pEV1->y() << ":" << pEV1->z() << endl;
    cout << dDotFace0 << ":" << dPara << ":" << dDot01 << endl;
    assert(0);
  }

  for (int ii = 0; ii < 3; ii++) {
    adPos[ii] = pInt[ii];
  }
  return iOrient01 + iOrient12 + iOrient20;
}

/*
 * Returns true iff the two triangular faces intersect.
 *
 * Currently does NOT handle degenerate cases, assumes coplanar triangles do not intersect.
 */
bool
clipIntersectFaceFace(map<edamFace*, vector<clipE>> &segmentClips,
		      edamFace* pEF0, edamFace* pEF1,
		      std::set<edamVert*> &garbVerts)
{

  // Make sure we didn't give this a quad face.
  assert(pEF0->getNumVerts() == 3);
  assert(pEF1->getNumVerts() == 3);

  edamVert* pEV0 = pEF0->getVert(0);
  edamVert* pEV1 = pEF0->getVert(1);
  edamVert* pEV2 = pEF0->getVert(2);
  edamVert* pEV3 = pEF1->getVert(0);
  edamVert* pEV4 = pEF1->getVert(1);
  edamVert* pEV5 = pEF1->getVert(2);

  const double* p0 = pEV0->getVert()->getCoords();
  const double* p1 = pEV1->getVert()->getCoords();
  const double* p2 = pEV2->getVert()->getCoords();
  const double* p3 = pEV3->getVert()->getCoords();
  const double* p4 = pEV4->getVert()->getCoords();
  const double* p5 = pEV5->getVert()->getCoords();

  // False if bounding boxes do not overlap.
  if (MAX3(p0[0],p1[0],p2[0]) < MIN3(p3[0], p4[0], p5[0])
      || MIN3(p0[0],p1[0],p2[0]) > MAX3(p3[0], p4[0], p5[0]))
    return false;
  if (MAX3(p0[1],p1[1],p2[1]) < MIN3(p3[1], p4[1], p5[1])
      || MIN3(p0[1],p1[1],p2[1]) > MAX3(p3[1], p4[1], p5[1]))
    return false;
  if (MAX3(p0[2],p1[2],p2[2]) < MIN3(p3[2], p4[2], p5[2])
      || MIN3(p0[2],p1[2],p2[2]) > MAX3(p3[2], p4[2], p5[2]))
    return false;

  // False if they share an edge.
  if (pEF1->hasVert(pEF0->getVert(0)) + pEF1->hasVert(pEF0->getVert(1))
      + pEF1->hasVert(pEF0->getVert(2)) == 2)
    return false;

  double adNorm0[3];
  pEF0->calcUnitNormal(adNorm0);
  double adNorm1[3];
  pEF1->calcUnitNormal(adNorm1);
  // False if the planes are parallel (hopefully...).
  if (fabs(dDOT3D(adNorm0, adNorm1)) > 1 - 1e-10)
    return false;

  // Check for edges of face 1 crossing face 0.
  double adPos34[3];
  int i34 = clipIntersectFaceLine(pEF0, pEV3, pEV4, adPos34);
  double adPos45[3];
  int i45 = clipIntersectFaceLine(pEF0, pEV4, pEV5, adPos45);
  double adPos53[3];
  int i53 = clipIntersectFaceLine(pEF0, pEV5, pEV3, adPos53);

  // Check for edges of face 0 crossing face 1.
  double adPos01[3];
  int i01 = clipIntersectFaceLine(pEF1, pEV0, pEV1, adPos01);
  double adPos12[3];
  int i12 = clipIntersectFaceLine(pEF1, pEV1, pEV2, adPos12);
  double adPos20[3];
  int i20 = clipIntersectFaceLine(pEF1, pEV2, pEV0, adPos20);

  if (!(i34 || i45 || i53 || i01 || i12 || i20)) {
    // No intersection found.
    return false;
  }

  // Store positions of edge-face intersections.
  vector<clipV> clipVerts;
  if (i34) {
    clipVerts.push_back(clipV(adPos34[0], adPos34[1], adPos34[2]));
  }
  if (i45) {
    clipVerts.push_back(clipV(adPos45[0], adPos45[1], adPos45[2]));
  }
  if (i53) {
    clipVerts.push_back(clipV(adPos53[0], adPos53[1], adPos53[2]));
  }
  if (i01) {
    clipVerts.push_back(clipV(adPos01[0], adPos01[1], adPos01[2]));
  }
  if (i12) {
    clipVerts.push_back(clipV(adPos12[0], adPos12[1], adPos12[2]));
  }
  if (i20) {
    clipVerts.push_back(clipV(adPos20[0], adPos20[1], adPos20[2]));
  }

  // If the intersecting faces share a vertex, that vertex is an endpoint of the intersection line.
  // It must be treated as a clip vertex in the triangulation routine, but it already exists, so store it in garbVerts to prevent it from being re-created.
  if (!(i34 + i45 + i53 + i01 + i12 + i20 == 6)) {

    assert(i34 + i45 + i53 + i01 + i12 + i20 == 3);
    edamVert* pEVc = NULL;
    if (pEF1->hasVert(pEF0->getVert(0)))
      pEVc = pEF0->getVert(0);
    if (pEF1->hasVert(pEF0->getVert(1)))
      pEVc = pEF0->getVert(1);
    if (pEF1->hasVert(pEF0->getVert(2)))
      pEVc = pEF0->getVert(2);
    assert(pEVc);

    const clipV tempCV(pEVc->x(), pEVc->y(), pEVc->z());
    clipVerts.push_back(tempCV);
    garbVerts.insert(pEVc);
  }

  // Add the intersection line to the list of intersection lines for both faces.
  assert(clipVerts.size() == 2);
  segmentClips.at(pEF0).push_back(clipE(clipVerts[0], clipVerts[1]));
  segmentClips.at(pEF1).push_back(clipE(clipVerts[0], clipVerts[1]));

  // Tag the vertices of these faces. They may be decimated later to smooth things out a bit.
  pEV0->setIntersectAdj(true);
  pEV1->setIntersectAdj(true);
  pEV2->setIntersectAdj(true);
  pEV3->setIntersectAdj(true);
  pEV4->setIntersectAdj(true);
  pEV5->setIntersectAdj(true);

  return true;
}

/*
 * Insert lines of intersection into each intersected face and retriangulate them. This involves the creation
 */
//void
//edamMesh3D_Closed::clipTriangulateFaces(
//    const int iLayer,
//				 vector<edamFace*> &segmentList,
//				 map<edamFace*, vector<clipE>> &segmentClips,
//				 vector<bool> segmentIntersected,
//				 set<edamVert*> &garbVerts)
//{
//
//  // Create bounding box tree to identify clip vertices that already exist.
//  // Initializing a tree with a dummy element, don't know how else to do it.
//  double (*a2dBBclip)[6] = new double[1][6];
//  a2dBBclip[0][0] = 0;
//  a2dBBclip[0][1] = 0;
//  a2dBBclip[0][2] = 0;
//  a2dBBclip[0][3] = 0;
//  a2dBBclip[0][4] = 0;
//  a2dBBclip[0][5] = 0;
//  ADT* ADTclip = new ADT(1, 3, ADT::eBBoxes,
//			 reinterpret_cast<double*>(a2dBBclip));
//  delete a2dBBclip;
//
//  // Corresponding dummy element in vert list.
//  vector<edamVert*> clipVerts;
//  clipVerts.push_back(NULL);
//
//  vector<edamFace*> VTKfaces;
//
//  // Add existing vertices that were endpoints of intersection lines (intersection between two faces that share a vertex).
//  for (edamVert* pEVgarb : garbVerts) {
//
//    double adBB[6];
//    adBB[0] = pEVgarb->x() - 1.e-12;
//    adBB[1] = pEVgarb->x() + 1.e-12;
//    adBB[2] = pEVgarb->y() - 1.e-12;
//    adBB[3] = pEVgarb->y() + 1.e-12;
//    adBB[4] = pEVgarb->z() - 1.e-12;
//    adBB[5] = pEVgarb->z() + 1.e-12;
//
//    vector<GR_index_t> result = ADTclip->veciRangeQuery(adBB);
//    assert(result.empty());
//
//    ADTclip->vAddNode(adBB, ADTclip->iTreeSize());
//    clipVerts.push_back(pEVgarb);
//  }
//
//  // Iterate over all intersected faces.
//  for (size_t iSeg = 0; iSeg < segmentList.size(); iSeg++) {
//    if (!segmentIntersected[iSeg])
//      continue;
//
//    edamFace* pEF = segmentList.at(iSeg);
//
//    // Will project into 2D space to triangulate, find best coordinate plane to project to.
//    // This prevents degeneracies where a triangle could project to a line in 2D.
//    double adNorm[3];
//    pEF->calcUnitNormal(adNorm);
//    int iProject = 0;
//
//    if (fabs(adNorm[1]) > fabs(adNorm[0]))
//      iProject = 1;
//    if (fabs(adNorm[2]) > fabs(adNorm[iProject]))
//      iProject = 2;
//    // Need to know if we're projecting the back or front of the face.
//
//    int iFlip = (adNorm[iProject] < 0 ? -1 : 1);
//
//    // Make sure there are actually intersection lines (clipEdges) through this face.
//    vector<clipE> clipEdges = segmentClips.at(pEF);
//    assert(!clipEdges.empty());
//
//    // Split intersection lines that intersect other intersection lines.
//    // This occurs at points where 3 faces intersected (currently not handling if more than 3 faces intersect at a point).
//    bool qDone = true;
//    do {
//      qDone = true;
//      for (size_t i = 0; i < clipEdges.size() - 1; i++) {
//	for (size_t j = i + 1; j < clipEdges.size(); j++) {
//	  clipE e0 = clipEdges[i];
//	  clipE e1 = clipEdges[j];
//
//	  clipV v00 = e0.v0;
//	  clipV v01 = e0.v1;
//	  clipV v10 = e1.v0;
//	  clipV v11 = e1.v1;
//
//	  double d0[2] =
//		{ v01.coords[(iProject + 1) % 3]
//		    - v00.coords[(iProject + 1) % 3], v01.coords[(iProject + 2)
//		    % 3] - v00.coords[(iProject + 2) % 3] };
//	  double d1[2] =
//		{ v11.coords[(iProject + 1) % 3]
//		    - v10.coords[(iProject + 1) % 3], v11.coords[(iProject + 2)
//		    % 3] - v10.coords[(iProject + 2) % 3] };
//
//	  // Assume the lines do not overlap if they are parallel.
//	  if (fabs(dDOT2D(d0, d1)) > dMAG2D(d0) * dMAG2D(d1) - 1e-10) {
//	    continue;
//	  }
//
//	  double t0 =
//	      ((v00.coords[(iProject + 1) % 3] - v10.coords[(iProject + 1) % 3])
//		  * d1[1]
//		  - (v00.coords[(iProject + 2) % 3]
//		      - v10.coords[(iProject + 2) % 3]) * d1[0])
//		  / (d1[0] * d0[1] - d1[1] * d0[0]);
//	  double t1 =
//	      ((v10.coords[(iProject + 2) % 3] - v00.coords[(iProject + 2) % 3])
//		  * d0[0]
//		  - (v10.coords[(iProject + 1) % 3]
//		      - v00.coords[(iProject + 1) % 3]) * d0[1])
//		  / (d1[0] * d0[1] - d1[1] * d0[0]);
//
//	  // If the line segments intersect, split each one and add a new clip vertex at the point of intersection.
//	  if (t0 > 1e-10 && t0 < 1 - 1e-10 && t1 > 1e-10 && t1 < 1 - 1e-10) {
//	    clipEdges.erase(clipEdges.begin() + j);
//	    clipEdges.erase(clipEdges.begin() + i);
//	    clipV vNew((1 - t0) * v00.coords[0] + t0 * v01.coords[0],
//		       (1 - t0) * v00.coords[1] + t0 * v01.coords[1],
//		       (1 - t0) * v00.coords[2] + t0 * v01.coords[2]);
//	    clipEdges.push_back(clipE(v00, vNew));
//	    clipEdges.push_back(clipE(v01, vNew));
//	    clipEdges.push_back(clipE(v10, vNew));
//	    clipEdges.push_back(clipE(v11, vNew));
//
//	    qDone = false;
//	    break;
//	  }
//	}
//	if (!qDone)
//	  break;
//      }
//    }
//    while (!qDone);
//
//    // Map edge endpoints to actual (likely new) vertices on the mesh.
//    set<edamVert*> clipSet;
//    for (size_t iE = 0; iE < clipEdges.size(); iE++) {
//      clipE e = clipEdges.at(iE);
//
//      // Check if the vertex already exists.
//      double adBB0[6];
//      adBB0[0] = e.v0.coords[0] - 1.e-12;
//      adBB0[1] = e.v0.coords[0] + 1.e-12;
//      adBB0[2] = e.v0.coords[1] - 1.e-12;
//      adBB0[3] = e.v0.coords[1] + 1.e-12;
//      adBB0[4] = e.v0.coords[2] - 1.e-12;
//      adBB0[5] = e.v0.coords[2] + 1.e-12;
//      vector<GR_index_t> result = ADTclip->veciRangeQuery(adBB0);
//
//      if (result.empty()) {
//	// If not, create it and add it to the search tree.
//	e.v0.pEV = createEdamVert(e.v0.coords[0], e.v0.coords[1],
//				  e.v0.coords[2], iLayer + 1);
//	ADTclip->vAddNode(adBB0, ADTclip->iTreeSize());
//
//	clipVerts.push_back(e.v0.pEV);
//	e.v0.pEV->setIntersectClip(true);
//      }
//      else {
//	assert(result.size() == 1);
//	e.v0.pEV = clipVerts.at(result.at(0));
//      }
//      clipSet.insert(e.v0.pEV);
//
//      double adBB1[6];
//      adBB1[0] = e.v1.coords[0] - 1.e-12;
//      adBB1[1] = e.v1.coords[0] + 1.e-12;
//      adBB1[2] = e.v1.coords[1] - 1.e-12;
//      adBB1[3] = e.v1.coords[1] + 1.e-12;
//      adBB1[4] = e.v1.coords[2] - 1.e-12;
//      adBB1[5] = e.v1.coords[2] + 1.e-12;
//      result = ADTclip->veciRangeQuery(adBB1);
//      if (result.empty()) {
//	e.v1.pEV = createEdamVert(e.v1.coords[0], e.v1.coords[1],
//				  e.v1.coords[2], iLayer + 1);
//	ADTclip->vAddNode(adBB1, ADTclip->iTreeSize());
//	clipVerts.push_back(e.v1.pEV);
//	e.v1.pEV->setIntersectClip(true);
//      }
//      else {
//	assert(result.size() == 1);
//	e.v1.pEV = clipVerts[result[0]];
//      }
//      clipSet.insert(e.v1.pEV);
//
//      clipEdges.at(iE) = e;
//    }
//
//    // Make sure not to doubly include the vertices on the triangle boundary.
//    clipSet.erase(pEF->getVert(0));
//    clipSet.erase(pEF->getVert(1));
//    clipSet.erase(pEF->getVert(2));
//
//    // List the vertices to triangulate.
//    vector<edamVert*> vertList;
//    vertList.push_back(pEF->getVert(0));
//    vertList.push_back(pEF->getVert(1));
//    vertList.push_back(pEF->getVert(2));
//    vertList.insert(vertList.end(), clipSet.begin(), clipSet.end());
//
//    map<int, int> realToTemp;
//    map<int, int> tempToReal;
//
//    ofstream fPoly;
//    fPoly.open("tri.poly");
//    fPoly.precision(16);
//    fPoly << vertList.size() << " 2 0 0" << endl << endl;
//
//    vector<vector<double>> coordList;
//    // Two-way mapping between indices in temp mesh and real indices.
//    for (size_t iV = 0; iV < vertList.size(); iV++) {
//      edamVert* pEV = vertList.at(iV);
//      realToTemp.insert(pair<int, int>(pEV->getIndex(), iV));
//      tempToReal.insert(pair<int, int>(iV, pEV->getIndex()));
//      const double* coords = pEV->getVert()->getCoords();
//      fPoly << iV << " " << coords[(iProject + 1) % 3] * iFlip << " "
//	  << coords[(iProject + 2) % 3] << endl;
//
//      vector<double> coord(
//	{ coords[(iProject + 1) % 3] * iFlip, coords[(iProject + 2) % 3] });
//      coordList.push_back(coord);
//    }
//
//    // Find segments on edges of triangle...
//    double coord0[2] =
//      { coordList[0][0], coordList[0][1] };
//    double coord1[2] =
//      { coordList[1][0], coordList[1][1] };
//    double coord2[2] =
//      { coordList[2][0], coordList[2][1] };
//
//    // Find vertices on edge between verts 0 and 1. Order them by distance from vert 0.
//    set<pair<double, int>> edge01;
//    for (size_t iV = 3; iV < vertList.size(); iV++) {
//      edamVert* pEVedge = vertList.at(iV);
//      const double* coords = pEVedge->getVert()->getCoords();
//      double coord[2] =
//	{ coords[(iProject + 1) % 3] * iFlip, coords[(iProject + 2) % 3] };
//
//      if (checkOrient2D(coord0, coord, coord1) == 0) {
//	edge01.insert(pair<double, int>(dDIST2D(coord0, coord), iV));
//      }
//    }
//
//    vector<pair<double, int>> list01(edge01.begin(), edge01.end());
//    // Connect them by new edges.
//    if (!list01.empty()) {
//      clipEdges.push_back(
//	  clipE(vertList.at(0), vertList.at(list01.front().second)));
//      clipEdges.push_back(
//	  clipE(vertList.at(1), vertList.at(list01.back().second)));
//      for (size_t iE = 0; iE < list01.size() - 1; iE++) {
//	clipEdges.push_back(
//	    clipE(vertList.at(list01.at(iE).second),
//		  vertList.at(list01.at(iE + 1).second)));
//      }
//    }
//    else {
//      clipEdges.push_back(clipE(vertList.at(0), vertList.at(1)));
//    }
//
//    // Verts 1 and 2.
//    set<pair<double, int>> edge12;
//    for (size_t iV = 3; iV < vertList.size(); iV++) {
//      edamVert* pEVedge = vertList.at(iV);
//      const double* coords = pEVedge->getVert()->getCoords();
//      double coord[2] =
//	{ coords[(iProject + 1) % 3] * iFlip, coords[(iProject + 2) % 3] };
//
//      if (checkOrient2D(coord1, coord, coord2) == 0) {
//	edge12.insert(pair<double, int>(dDIST2D(coord1, coord), iV));
//      }
//    }
//    vector<pair<double, int>> list12(edge12.begin(), edge12.end());
//    if (!list12.empty()) {
//      clipEdges.push_back(
//	  clipE(vertList.at(1), vertList.at(list12.front().second)));
//      clipEdges.push_back(
//	  clipE(vertList.at(2), vertList.at(list12.back().second)));
//      for (size_t iE = 0; iE < list12.size() - 1; iE++) {
//	clipEdges.push_back(
//	    clipE(vertList.at(list12.at(iE).second),
//		  vertList.at(list12.at(iE + 1).second)));
//      }
//    }
//    else {
//      clipEdges.push_back(clipE(vertList.at(1), vertList.at(2)));
//    }
//
//    // Verts 2 and 0.
//    set<pair<double, int>> edge20;
//    for (size_t iV = 3; iV < vertList.size(); iV++) {
//      edamVert* pEVedge = vertList.at(iV);
//      const double* coords = pEVedge->getVert()->getCoords();
//      double coord[2] =
//	{ coords[(iProject + 1) % 3] * iFlip, coords[(iProject + 2) % 3] };
//
//      if (checkOrient2D(coord2, coord, coord0) == 0) {
//	edge20.insert(pair<double, int>(dDIST2D(coord2, coord), iV));
//      }
//    }
//    vector<pair<double, int>> list20(edge20.begin(), edge20.end());
//    if (!list20.empty()) {
//      clipEdges.push_back(
//	  clipE(vertList.at(2), vertList.at(list20.front().second)));
//      clipEdges.push_back(
//	  clipE(vertList.at(0), vertList.at(list20.back().second)));
//      for (size_t iE = 0; iE < list20.size() - 1; iE++) {
//	clipEdges.push_back(
//	    clipE(vertList.at(list20.at(iE).second),
//		  vertList.at(list20.at(iE + 1).second)));
//      }
//    }
//    else {
//      clipEdges.push_back(clipE(vertList.at(2), vertList.at(0)));
//    }
//
//    fPoly << endl << clipEdges.size() << " 0" << endl << endl;
//
//    for (size_t iE = 0; iE < clipEdges.size(); iE++) {
//      clipE e = clipEdges.at(iE);
//      int i0 = realToTemp.at(e.v0.pEV->getIndex());
//      int i1 = realToTemp.at(e.v1.pEV->getIndex());
//
//      fPoly << iE << " " << i0 << " " << i1 << endl;
//    }
//
//    fPoly << endl << "0";
//
//    fPoly.close();
//
//    // Use JRS triangle to triangulate the set of input vertices and edges inside the triangular face.
//    if (system("/home/cfog/Downloads/triangle/triangle -p tri.poly > /dev/null")
//	!= 0)
//      assert(0);
//
//    ifstream fEle;
//    fEle.open("tri.1.ele");
//    string strLine;
//    std::getline(fEle, strLine);
//
//    // Insert the triangles from the triangulation back into the original mesh.
//    int tri[4];
//    while (1) {
//      string line;
//      std::getline(fEle, line);
//      istringstream ISS(line);
//
//      string first;
//      ISS >> first;
//      if (first == "#")
//	break;
//
//      ISS >> tri[0] >> tri[1] >> tri[2];
//      assert(ISS.eof());
//
//      edamVert* pEV0 = m_eVL.get(tempToReal.at(tri[0]));
//      edamVert* pEV1 = m_eVL.get(tempToReal.at(tri[1]));
//      edamVert* pEV2 = m_eVL.get(tempToReal.at(tri[2]));
//
//      edamFace* pEFvtk = createTriFace(pEV0, pEV1, pEV2);
//      VTKfaces.push_back(pEFvtk);
//    }
//
//    // Delete the split triangle.
//    deleteTriFace(pEF);
//  }
//  delete ADTclip;
//}
/*
 *  Find pairs of front triangles that intersect, and clip them at lines of intersection.
 */
//bool
//edamMesh3D_Closed::clipIntersections(const int iLayer,
//				     vector<edamFace*> &segmentList)
//{
//  bool qRet = false;
//
//  // Initialization.
//  vector<bool> segmentIntersected;
//  map<edamFace*, vector<clipE>> segmentClips;
//  vector<clipE> emptyVec;
//  for (edamFace* pEF : segmentList) {
//    segmentIntersected.push_back(false);
//    segmentClips.insert(pair<edamFace*, vector<clipE>>(pEF, emptyVec));
//  }
//  set<edamVert*> garbVerts;
//
//  // Construct a bounding box tree for efficient intersection checks.
//  double (*a2dBB)[6] = new double[segmentList.size()][6];
//  for (size_t iSeg = 0; iSeg < segmentList.size(); iSeg++) {
//    edamFace* pEF = segmentList.at(iSeg);
//    edamVert* pEV0 = pEF->getVert(0);
//    edamVert* pEV1 = pEF->getVert(1);
//    edamVert* pEV2 = pEF->getVert(2);
//    a2dBB[iSeg][0] = MIN3(pEV0->x(), pEV1->x(), pEV2->x());
//    a2dBB[iSeg][1] = MAX3(pEV0->x(), pEV1->x(), pEV2->x());
//    a2dBB[iSeg][2] = MIN3(pEV0->y(), pEV1->y(), pEV2->y());
//    a2dBB[iSeg][3] = MAX3(pEV0->y(), pEV1->y(), pEV2->y());
//    a2dBB[iSeg][4] = MIN3(pEV0->z(), pEV1->z(), pEV2->z());
//    a2dBB[iSeg][5] = MAX3(pEV0->z(), pEV1->z(), pEV2->z());
//  }
//  ADT* ADTtree = new ADT(segmentList.size(), 3, ADT::eBBoxes,
//			 reinterpret_cast<double*>(a2dBB));
//  delete[] a2dBB;
//
//  // Iterate over all faces to detect intersections.
//  for (size_t iSeg = 0; iSeg < segmentList.size(); iSeg++) {
//    edamFace* pEF = segmentList.at(iSeg);
//    edamVert* pEV0 = pEF->getVert(0);
//    edamVert* pEV1 = pEF->getVert(1);
//    edamVert* pEV2 = pEF->getVert(2);
//    double adBB[6];
//    adBB[0] = MIN3(pEV0->x(), pEV1->x(), pEV2->x());
//    adBB[1] = MAX3(pEV0->x(), pEV1->x(), pEV2->x());
//    adBB[2] = MIN3(pEV0->y(), pEV1->y(), pEV2->y());
//    adBB[3] = MAX3(pEV0->y(), pEV1->y(), pEV2->y());
//    adBB[4] = MIN3(pEV0->z(), pEV1->z(), pEV2->z());
//    adBB[5] = MAX3(pEV0->z(), pEV1->z(), pEV2->z());
//
//    // Query the bounding box tree for nearby faces.
//    vector<GR_index_t> result = ADTtree->veciRangeQuery(adBB);
//
//    // Check if each nearby face actually intersects this one.
//    for (size_t ii : result) {
//      if (ii > iSeg) {
//	if (clipIntersectFaceFace(segmentClips, segmentList.at(iSeg),
//				  segmentList.at(ii), garbVerts)) {
//	  // At least one intersection has been found, so we will return true to trigger the intersection cleanup routines.
//	  qRet = true;
//
//	  // Mark intersected faces.
//	  segmentIntersected[iSeg] = true;
//	  segmentIntersected[ii] = true;
//	}
//      }
//    }
//  }
//
//  // Triangulate each intersected triangle with its intersection lines.
//  clipTriangulateFaces(iLayer, segmentList, segmentClips, segmentIntersected,
//		       garbVerts);
//
//  delete ADTtree;
//  return qRet;
//}
/*
 * Find non-manifold verts and split each one into multiple separate manifold verts.
 * Resulting front is composed of multiple disconnected, closed, manifold surfaces.
 *
 * Tag each vertex with an index iLoop indicating which of these surfaces it belongs to.
 */
//void
//edamMesh3D_Closed::reconnectIntersections(const int iLayer, int &iNumLoops)
//{
//
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//    edamVert* pEV0 = m_eVL.get(iVert);
//
//    vector<edamFace*> adjFaces = pEV0->getFaceList();
//    set<edamFace*> adjFaceSet(adjFaces.begin(), adjFaces.end());
//
//    bool qNonManifold = false;
//
//    while (!adjFaceSet.empty()) {
//
//      set<edamFace*> adjSubset;
//
//      // Starting from an arbitrary start face, move to the next adjacent face until we loop back to the start face.
//      //
//      // For manifold vertices, this traverses all adjacent faces, so this while(!adjFaceSet.empty()) loop
//      // executes only once. Otherwise, pick a new start face and repeat until no faces are left.
//      //
//      // For non-manifold vertices, faces traversed in each execution of the loop
//      // are split off by duplicating the vertex.
//      edamFace* pEFstart = *(adjFaceSet.begin());
//      edamFace* pEF = pEFstart;
//
//      do {
//	// Mark the face as visited.
//	adjSubset.insert(pEF);
//
//	// Find the next adjacent face (Clockwise around the vertex).
//	edamVert* pEV1 = pEF->getNextVert(pEV0);
//	edamVert* pEV2 = pEF->getNextVert(pEV1);
//	vector<edamFace*> cands = pEV0->getCommonFaces(pEV1);
//
//	// There should be at most 4 faces adjacent to this edge, except in a degenerate
//	// case where 3 triangles intersected along a line, which is NOT handled.
//	if (cands.size() == 4) {
//
//	  for (edamFace* pEFcand : cands) {
//	    if (pEFcand->getNextVert(pEV0) == pEV1)
//	      continue;
//	    if (checkOrient3D(pEV0->getVert(), pEV1->getVert(), pEV2->getVert(),
//			      pEFcand->getNextVert(pEV0)->getVert()) != 0) {
//	      pEF = pEFcand;
//	    }
//	  }
//	}
//	else if (cands.size() == 2) {
//	  for (edamFace* pEFcand : cands) {
//	    if (pEFcand->getNextVert(pEV0) == pEV1)
//	      continue;
//	    pEF = pEFcand;
//	  }
//	}
//	else {
//	  assert(0);
//	}
//
//      }
//      while (pEF != pEFstart);
//
//      // Break out of the loop if the vertex is manifold. No work to be done here.
//      if (adjSubset.size() != adjFaceSet.size())
//	qNonManifold = true;
//      if (!qNonManifold)
//	break;
//
//      // Create a "new" vertex in the same place, and reconnect the traversed faces to it.
//      edamVert* pEVnew = createEdamVert(pEV0->x(), pEV0->y(), pEV0->z(),
//					pEV0->getLayer());
//
//      pEVnew->setIntersectClip2(true);
//      for (edamFace* pEFsub : adjSubset) {
//	adjFaceSet.erase(pEFsub);
//
//	edamVert* pEV1 = pEFsub->getNextVert(pEV0);
//	edamVert* pEV2 = pEFsub->getNextVert(pEV1);
//	deleteTriFace(pEFsub);
//	createTriFace(pEVnew, pEV1, pEV2);
//      }
//    }
//
//    if (qNonManifold) {
//      // Delete the old vertex - it has been replaced by new ones.
//      deleteEdamVert(pEV0);
//    }
//  }
//  //The front should consist of disconnected manifolds now.
//
//  //Traverse each closed manifold in the front (by flood-fill) and assign loop tags.
//  iNumLoops = 0;
//  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//      iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->getLoop() > 0)
//      continue;
//
//    iNumLoops++;
//    pEV->setLoop(iNumLoops);
//
//    vector<edamVert*> currVerts;
//    currVerts.push_back(pEV);
//
//    while (!currVerts.empty()) {
//      vector<edamVert*> nextVerts;
//
//      for (edamVert* pEVcurr : currVerts) {
//	vector<edamFace*> adjFaces = pEVcurr->getFaceList();
//	for (edamFace* pEF : adjFaces) {
//
//	  edamVert* pEVadj = pEF->getNextVert(pEVcurr);
//	  if (pEVadj->getLoop() > 0)
//	    continue;
//
//	  pEVadj->setLoop(iNumLoops);
//	  nextVerts.push_back(pEVadj);
//	}
//      }
//
//      currVerts = nextVerts;
//    }
//  }
//}
/*
 * Delete the loops that are not valid parts of the new front (identified using winding numbers, see comment in code).
 */
//void
//edamMesh3D_Closed::findInvertedLoops(const int iLayer, const int iNumLoops)
//{
//
//  // Identify and label separate closed surfaces ("loops", using 2D terminology). Store faces belonging to each loop.
//
//  vector<set<edamFace*>> loops(iNumLoops);
//  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//      iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//    int iLoop = pEV->getLoop();
//
//    vector<edamFace*> adjFaces = pEV->getFaceList();
//    for (edamFace* pEF : adjFaces) {
//      loops.at(iLoop - 1).insert(pEF);
//
//    }
//  }
//
//  // Store the set of all faces (for the winding number computation).
//  set<edamFace*> allFaces;
//  for (set<edamFace*> loop : loops) {
//    allFaces.insert(loop.begin(), loop.end());
//  }
//
//  /*
//   * For each loop, find a point immediately outside (a.k.a. offset in the marching direction from) the loop and compute the winding number
//   * of the whole front around that point.
//   * Every loop which is a valid new front gives a winding number of...
//   * 		0, if there is no outer boundary
//   * 		4*PI, otherwise
//   */
//  vector<bool> loopsBad;
//  for (set<edamFace*> loop : loops) {
//    if (loop.empty())
//      continue;
//
//    edamFace* pEFstart = *(loop.begin());
//    double adOffset[3];
//
//    do {
//      double adCent[3];
//      double adNorm[3];
//      pEFstart->getFace()->calcCentroid(adCent);
//      pEFstart->calcUnitNormal(adNorm);
//      for (int ii = 0; ii < 3; ii++) {
//	adOffset[ii] = adCent[ii] - 1.e-7 * adNorm[ii];
//      }
//
//      cout << "Centroid: " << adCent[0] << " " << adCent[1] << " " << adCent[2]
//	  << endl;
//      cout << "Offset: " << adOffset[0] << " " << adOffset[1] << " "
//	  << adOffset[2] << endl;
//    }
//    while (0);
//
//    // Compute winding number as sum of contributions from individual faces.
//    double dWind = 0;
//    for (edamFace* pEF : allFaces) {
//      const double* ad0 = pEF->getVert(0)->getVert()->getCoords();
//      const double* ad1 = pEF->getVert(1)->getVert()->getCoords();
//      const double* ad2 = pEF->getVert(2)->getVert()->getCoords();
//      double adA[3] = adDIFF3D(ad0,adOffset);
//      double adB[3] = adDIFF3D(ad1,adOffset);
//      double adC[3] = adDIFF3D(ad2,adOffset);
//
//      double dNumer = 0;
//      dNumer += adA[0] * (adB[1] * adC[2] - adB[2] * adC[1]);
//      dNumer += adA[1] * (adB[2] * adC[0] - adB[0] * adC[2]);
//      dNumer += adA[2] * (adB[0] * adC[1] - adB[1] * adC[0]);
//      double dDenom = dMAG3D(adA) * dMAG3D(adB) * dMAG3D(adC);
//      dDenom += dDOT3D(adA,adB) * dMAG3D(adC);
//      dDenom += dDOT3D(adB,adC) * dMAG3D(adA);
//      dDenom += dDOT3D(adC,adA) * dMAG3D(adB);
//
//      double dSolidAngle = 2 * atan(dNumer / dDenom);
//      double adNorm[3];
//      pEF->calcUnitNormal(adNorm);
//      if (dSolidAngle > 0 && dDOT3D(adNorm,adA) < 0) {
//	dSolidAngle -= 2 * M_PI;
//      }
//      else if (dSolidAngle < 0 && dDOT3D(adNorm,adA) > 0) {
//	dSolidAngle += 2 * M_PI;
//      }
//      dWind += dSolidAngle;
//    }
//
//    cout << "Winding Number: " << dWind << endl;
//
//    //if(fabs(dWind-4*M_PI) > M_PI){	// TODO: automatically figure out the right winding number (or at least take it as a user input...)
//    if (fabs(dWind) > 1.e-12) {
//
//      cout << "loop is bad" << endl;
//      loopsBad.push_back(true);
//    }
//    else {
//      cout << "loop is good" << endl;
//      loopsBad.push_back(false);
//    }
//  }
//
//  // Now that we've identified bad loops, delete them.
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//
//    edamVert* pEV = m_eVL.get(iVert);
//    int iLoop = pEV->getLoop();
//    if (!loopsBad[iLoop - 1])
//      continue;
//
//    vector<edamFace*> adjFaces = pEV->getFaceList();
//    for (edamFace* pEF : adjFaces) {
//      deleteTriFace(pEF);
//
//    }
//
//    deleteEdamVert(pEV);
//  }
//}
// From "On the 'most normal' normal - Aubry - 2007", guarantees a normal visible to all faces as long as one exists.
// Much slower than the averaging approach, so we use this only where necessary.
void
bestNorm(vector<vector<double>> faceNorms, double adBestNorm[], edamVert* pEV)
{

  double scalmin = -1;
  double scal;
  double scalt;

  bool qReasonable = false;

  // Aubry 2007 Figure 4: Compute midpoints between pairs of normals.
  int n = faceNorms.size();
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      vector<double> Ni = faceNorms[i];
      vector<double> Nj = faceNorms[j];
      double Nb[3] = adSUM3D(Ni,Nj);
      NORMALIZE3D(Nb);

      scal = dDOT3D(Ni, Nb);
      if (scal < scalmin)
	continue;

      bool qBad = false;
      for (int p = 0; p < n; p++) {
	if (p == i || p == j)
	  continue;
	vector<double> Np = faceNorms[p];
	scalt = dDOT3D(Np, Nb);
	if (iFuzzyComp(scalt, scal) == -1) {
	  qBad = true;
	  break;
	}
      }
      if (qBad)
	continue;

      scalmin = scal;
      for (int ii = 0; ii < 3; ii++) {
	adBestNorm[ii] = Nb[ii];
      }
      qReasonable = true;
    }
  }

  // Aubry 2007 Figure 5: Compute centres of circles defined by triples of normals.
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      for (int k = j + 1; k < n; k++) {

	vector<double> Ni = faceNorms[i];
	vector<double> Nj = faceNorms[j];
	vector<double> Nk = faceNorms[k];
	double Nc[3];

	if (dDOT3D(Ni,Nj) > 1 - 1.e-14)
	  continue;
	if (dDOT3D(Ni,Nk) > 1 - 1.e-14)
	  continue;
	if (dDOT3D(Nj,Nk) > 1 - 1.e-14)
	  continue;

	double denom = (Ni[0] - Nj[0]) * (Ni[1] - Nk[1])
	    - (Ni[0] - Nk[0]) * (Ni[1] - Nj[1]);
	Nc[0] = ((Ni[2] - Nk[2]) * (Ni[1] - Nj[1])
	    - (Ni[2] - Nj[2]) * (Ni[1] - Nk[1])) / denom;
	Nc[1] = ((Ni[0] - Nk[0]) * (Ni[2] - Nj[2])
	    - (Ni[0] - Nj[0]) * (Ni[2] - Nk[2])) / denom;
	Nc[2] = 1. / sqrt(1 + Nc[0] * Nc[0] + Nc[1] * Nc[1]);
	Nc[0] = Nc[0] * Nc[2];
	Nc[1] = Nc[1] * Nc[2];
	scal = dDOT3D(Nc, Ni);

	// Nine permutations of the same calculation, to handle degenerate cases.
	if (!finite(scal)) {
	  denom = (Ni[1] - Nj[1]) * (Ni[2] - Nk[2])
	      - (Ni[1] - Nk[1]) * (Ni[2] - Nj[2]);
	  Nc[1] = ((Ni[0] - Nk[0]) * (Ni[2] - Nj[2])
	      - (Ni[0] - Nj[0]) * (Ni[2] - Nk[2])) / denom;
	  Nc[2] = ((Ni[1] - Nk[1]) * (Ni[0] - Nj[0])
	      - (Ni[1] - Nj[1]) * (Ni[0] - Nk[0])) / denom;
	  Nc[0] = 1. / sqrt(1 + Nc[1] * Nc[1] + Nc[2] * Nc[2]);
	  Nc[1] = Nc[1] * Nc[0];
	  Nc[2] = Nc[2] * Nc[0];
	  scal = dDOT3D(Nc, Ni);

	}
	if (!finite(scal)) {
	  denom = (Ni[2] - Nj[2]) * (Ni[0] - Nk[0])
	      - (Ni[2] - Nk[2]) * (Ni[0] - Nj[0]);
	  Nc[2] = ((Ni[1] - Nk[1]) * (Ni[0] - Nj[0])
	      - (Ni[1] - Nj[1]) * (Ni[0] - Nk[0])) / denom;
	  Nc[0] = ((Ni[2] - Nk[2]) * (Ni[1] - Nj[1])
	      - (Ni[2] - Nj[2]) * (Ni[1] - Nk[1])) / denom;
	  Nc[1] = 1. / sqrt(1 + Nc[2] * Nc[2] + Nc[0] * Nc[0]);
	  Nc[2] = Nc[2] * Nc[1];
	  Nc[0] = Nc[0] * Nc[1];
	  scal = dDOT3D(Nc, Ni);
	}
	if (!finite(scal)) {
	  denom = (Nj[0] - Nk[0]) * (Nj[1] - Ni[1])
	      - (Nj[0] - Ni[0]) * (Nj[1] - Nk[1]);
	  Nc[0] = ((Nj[2] - Ni[2]) * (Nj[1] - Nk[1])
	      - (Nj[2] - Nk[2]) * (Nj[1] - Ni[1])) / denom;
	  Nc[1] = ((Nj[0] - Ni[0]) * (Nj[2] - Nk[2])
	      - (Nj[0] - Nk[0]) * (Nj[2] - Ni[2])) / denom;
	  Nc[2] = 1. / sqrt(1 + Nc[0] * Nc[0] + Nc[1] * Nc[1]);
	  Nc[0] = Nc[0] * Nc[2];
	  Nc[1] = Nc[1] * Nc[2];
	  scal = dDOT3D(Nc, Nj);
	}
	if (!finite(scal)) {
	  denom = (Nj[1] - Nk[1]) * (Nj[2] - Ni[2])
	      - (Nj[1] - Ni[1]) * (Nj[2] - Nk[2]);
	  Nc[1] = ((Nj[0] - Ni[0]) * (Nj[2] - Nk[2])
	      - (Nj[0] - Nk[0]) * (Nj[2] - Ni[2])) / denom;
	  Nc[2] = ((Nj[1] - Ni[1]) * (Nj[0] - Nk[0])
	      - (Nj[1] - Nk[1]) * (Nj[0] - Ni[0])) / denom;
	  Nc[0] = 1. / sqrt(1 + Nc[1] * Nc[1] + Nc[2] * Nc[2]);
	  Nc[1] = Nc[1] * Nc[0];
	  Nc[2] = Nc[2] * Nc[0];
	  scal = dDOT3D(Nc, Nj);
	}
	if (!finite(scal)) {
	  denom = (Nj[2] - Nk[2]) * (Nj[0] - Ni[0])
	      - (Nj[2] - Ni[2]) * (Nj[0] - Nk[0]);
	  Nc[2] = ((Nj[1] - Ni[1]) * (Nj[0] - Nk[0])
	      - (Nj[1] - Nk[1]) * (Nj[0] - Ni[0])) / denom;
	  Nc[0] = ((Nj[2] - Ni[2]) * (Nj[1] - Nk[1])
	      - (Nj[2] - Nk[2]) * (Nj[1] - Ni[1])) / denom;
	  Nc[1] = 1. / sqrt(1 + Nc[2] * Nc[2] + Nc[0] * Nc[0]);
	  Nc[2] = Nc[2] * Nc[1];
	  Nc[0] = Nc[0] * Nc[1];
	  scal = dDOT3D(Nc, Nj);
	}
	if (!finite(scal)) {
	  denom = (Nk[0] - Ni[0]) * (Nk[1] - Nj[1])
	      - (Nk[0] - Nj[0]) * (Nk[1] - Ni[1]);
	  Nc[0] = ((Nk[2] - Nj[2]) * (Nk[1] - Ni[1])
	      - (Nk[2] - Ni[2]) * (Nk[1] - Nj[1])) / denom;
	  Nc[1] = ((Nk[0] - Nj[0]) * (Nk[2] - Ni[2])
	      - (Nk[0] - Ni[0]) * (Nk[2] - Nj[2])) / denom;
	  Nc[2] = 1. / sqrt(1 + Nc[0] * Nc[0] + Nc[1] * Nc[1]);
	  Nc[0] = Nc[0] * Nc[2];
	  Nc[1] = Nc[1] * Nc[2];
	  scal = dDOT3D(Nc, Nk);
	}
	if (!finite(scal)) {
	  denom = (Nk[1] - Ni[1]) * (Nk[2] - Nj[2])
	      - (Nk[1] - Nj[1]) * (Nk[2] - Ni[2]);
	  Nc[1] = ((Nk[0] - Nj[0]) * (Nk[2] - Ni[2])
	      - (Nk[0] - Ni[0]) * (Nk[2] - Nj[2])) / denom;
	  Nc[2] = ((Nk[1] - Nj[1]) * (Nk[0] - Ni[0])
	      - (Nk[1] - Ni[1]) * (Nk[0] - Nj[0])) / denom;
	  Nc[0] = 1. / sqrt(1 + Nc[1] * Nc[1] + Nc[2] * Nc[2]);
	  Nc[1] = Nc[1] * Nc[0];
	  Nc[2] = Nc[2] * Nc[0];
	  scal = dDOT3D(Nc, Nk);
	}
	if (!finite(scal)) {
	  denom = (Nk[2] - Ni[2]) * (Nk[0] - Nj[0])
	      - (Nk[2] - Nj[2]) * (Nk[0] - Ni[0]);
	  Nc[2] = ((Nk[1] - Nj[1]) * (Nk[0] - Ni[0])
	      - (Nk[1] - Ni[1]) * (Nk[0] - Nj[0])) / denom;
	  Nc[0] = ((Nk[2] - Nj[2]) * (Nk[1] - Ni[1])
	      - (Nk[2] - Ni[2]) * (Nk[1] - Nj[1])) / denom;
	  Nc[1] = 1. / sqrt(1 + Nc[2] * Nc[2] + Nc[0] * Nc[0]);
	  Nc[2] = Nc[2] * Nc[1];
	  Nc[0] = Nc[0] * Nc[1];
	  scal = dDOT3D(Nc, Nk);
	}

	if (!finite(scal)) {
	  cout << "Vert: " << pEV->x() << " " << pEV->y() << " " << pEV->z()
	      << endl;
	  cout << "Ni: " << Ni[0] << " " << Ni[1] << " " << Ni[2] << endl;
	  cout << "Nj: " << Nj[0] << " " << Nj[1] << " " << Nj[2] << endl;
	  cout << "Nk: " << Nk[0] << " " << Nk[1] << " " << Nk[2] << endl;
	  assert(finite(scal));
	}

	scal = dDOT3D(Nc, Ni);
	if (scal < 0) {
	  Nc[0] = -Nc[0];
	  Nc[1] = -Nc[1];
	  Nc[2] = -Nc[2];
	  scal = -scal;
	}
	if (scal < scalmin)
	  continue;

	bool qBad = false;
	for (int p = 0; p < n; p++) {
	  if (p == i || p == j || p == k)
	    continue;
	  vector<double> Np = faceNorms[p];
	  scalt = dDOT3D(Np, Nc);
	  if (iFuzzyComp(scalt, scal) == -1) {
	    qBad = true;
	    break;
	  }
	}
	if (qBad)
	  continue;

	scalmin = scal;
	for (int ii = 0; ii < 3; ii++) {
	  adBestNorm[ii] = Nc[ii];
	}
	qReasonable = true;
      }
    }
  }

  assert(qReasonable && "No normal direction visible to all adjacent faces.");
  assert(scalmin > 0);
}

/*
 * Get a bunch of information about a vertex and its in-layer neighborhood, most importantly the surface normal at the vertex.
 *
 * Return values are:
 * adUnitNorm - Unit vertex normal.
 * adAvgDir - Average of normalized directions to neighbors (not itself normalized).
 * adAvgDist - Average of non-normalized vectors to neighbors.
 * dMinDist - Minimum distance from the vertex to the polygonal boundary of its neighborhood (edges too, not just vertices).
 */
void
neighborhoodInfo(edamVert* pEV, set<edamFace*> &adjFaces, double* adUnitNorm,
		 double* adAvgDir, double* adAvgDist, double &dMinDist,
		 bool qSmoothing = false)
{

  dMinDist = LARGE_DBL;

  int iNumFaces = adjFaces.size();

  double adDirSum[3] =
    { 0, 0, 0 };
  double adDistSum[3] =
    { 0, 0, 0 };

  // Try to compute the normal as an angle-weighted average of adjacent face normals.
  double adBestNorm[] =
    { 0, 0, 0 };

  for (edamFace* pEF : adjFaces) {
    double adNorm[3];
    pEF->calcUnitNormal(adNorm);
    double dAngle = pEF->calcAngle(pEV);
    for (int ii = 0; ii < 3; ii++) {
      adBestNorm[ii] -= adNorm[ii] * dAngle;
    }
  }
  NORMALIZE3D(adBestNorm);

  // If the computed normal is not visible to all adjacent faces (likely at sharp or complex corners),
  // compute a better (but more expensive) normal.
  if (!qSmoothing) {
    for (edamFace* pEFbad : adjFaces) {
      double adNormBad[3];
      pEFbad->calcUnitNormal(adNormBad);

      if (dDOT3D(adNormBad,adBestNorm) > 0) {

	vector<vector<double>> faceNorms;
	for (edamFace* pEF : adjFaces) {
	  double adNorm[3];
	  vector<double> vecNorm(3);
	  pEF->calcUnitNormal(adNorm);
	  for (int ii = 0; ii < 3; ii++) {
	    vecNorm[ii] = -adNorm[ii];
	  }
	  faceNorms.push_back(vecNorm);
	}
	bestNorm(faceNorms, adBestNorm, pEV);

	break;
      }
    }
  }
  if (pEV->isNoDecimate()) {
    vector<vector<double>> faceNorms;
    for (edamFace* pEF : adjFaces) {
      double adNorm[3];
      vector<double> vecNorm(3);
      pEF->calcUnitNormal(adNorm);
      for (int ii = 0; ii < 3; ii++) {
	vecNorm[ii] = -adNorm[ii];
      }
      faceNorms.push_back(vecNorm);
    }
    bestNorm(faceNorms, adBestNorm, pEV);
  }

  // Find other information about the vertex neighborhood.
  double xx = pEV->x(), yy = pEV->y(), zz = pEV->z();

  for (edamFace* pEF : adjFaces) {
    edamVert* pEV1 = pEF->getNextVert(pEV);
    edamVert* pEV2 = pEF->getNextVert(pEV1);

    double adDist1[3] =
      { pEV1->x() - xx, pEV1->y() - yy, pEV1->z() - zz };
    double adDist2[3] =
      { pEV2->x() - xx, pEV2->y() - yy, pEV2->z() - zz };

    double adDist21[3] = adDIFF3D(adDist1,adDist2);
    double dProj = dDOT3D(adDist21, adDist1)
	/ (dMAG3D(adDist21) * dMAG3D(adDist1));
    double adMinDist[3];
    for (int ii = 0; ii < 3; ii++) {
      adMinDist[ii] = adDist1[ii] - adDist21[ii] * dProj;
    }
    double dMinDistProj = sqrt(
    dMAG3D_SQ(adMinDist) - pow(dDOT3D(adBestNorm, adMinDist), 2.0));

    dMinDist = fmin(dMinDist, dMinDistProj);

//    for (int ii = 0; ii < 3; ii++) {
//      adDistSum[ii] += adDist1[ii];
//    }
    NORMALIZE3D(adDist1);
    for (int ii = 0; ii < 3; ii++) {
      adDirSum[ii] += adDist1[ii];
    }
  }

  // For the averages, divide by the number of neighbors.
  for (int ii = 0; ii < 3; ii++) {
    adDirSum[ii] /= iNumFaces;
    adDistSum[ii] /= iNumFaces;
  }

  for (int ii = 0; ii < 3; ii++) {
    adUnitNorm[ii] = adBestNorm[ii];
    adAvgDir[ii] = adDirSum[ii];
    adAvgDist[ii] = adDistSum[ii];
  }
}

/*
 * Doing this with TetGen now for several reasons:
 *
 * 1) Avoid wasting time with all the tetra preprocessing.
 * 2) Don't have to mess up other GRUMMP apps while trying to make this work.
 * 3) Less mysterious boundary recovery failures (hopefully none).
 * 4) Guarantees no boundary Steiner points if requested.
 */
//void
//edamMesh3D_Closed::TetGenFill(const int iLayer, vector<edamFace*> &bottomFaces,
//		       vector<edamFace*> &topFaces)
//{
//  map<int, int> realToTemp;
//  map<int, int> tempToReal;
//
//  // Get list of verts on the boundary of the region.
//  std::set<edamVert*> pEVset;
//  for (edamFace* pEF : bottomFaces) {
//    pEVset.insert(pEF->getVert(0));
//    pEVset.insert(pEF->getVert(1));
//    pEVset.insert(pEF->getVert(2));
//  }
//  for (edamFace* pEF : topFaces) {
//    pEVset.insert(pEF->getVert(0));
//    pEVset.insert(pEF->getVert(1));
//    pEVset.insert(pEF->getVert(2));
//  }
//  vector<edamVert*> verts(pEVset.begin(), pEVset.end());
//
//  ofstream fPoly;
//  fPoly.open("fill.poly");
//  fPoly.precision(16);
//  fPoly << verts.size() << " 3 0 0" << endl << endl;
//
//  // Two-way mapping between indices in temp mesh and real indices.
//  for (size_t iV = 0; iV < verts.size(); iV++) {
//    edamVert* pEV = verts.at(iV);
//    realToTemp.insert(pair<int, int>(pEV->getIndex(), iV));
//    tempToReal.insert(pair<int, int>(iV, pEV->getIndex()));
//    fPoly << iV << " " << pEV->x() << " " << pEV->y() << " " << pEV->z()
//	<< endl;
//  }
//
//  fPoly << endl << bottomFaces.size() + topFaces.size() << " 0" << endl << endl;
//
//  for (edamFace* pEF : bottomFaces) {
//    int i0 = realToTemp.at(pEF->getVert(0)->getIndex());
//    int i1 = realToTemp.at(pEF->getVert(1)->getIndex());
//    int i2 = realToTemp.at(pEF->getVert(2)->getIndex());
//
//    fPoly << "1" << endl;
//    fPoly << "3 " << i0 << " " << i1 << " " << i2 << endl;
//  }
//  for (edamFace* pEF : topFaces) {
//    int i0 = realToTemp.at(pEF->getVert(0)->getIndex());
//    int i1 = realToTemp.at(pEF->getVert(1)->getIndex());
//    int i2 = realToTemp.at(pEF->getVert(2)->getIndex());
//
//    fPoly << "1" << endl;
//    fPoly << "3 " << i0 << " " << i2 << " " << i1 << endl;
//  }
//
//  fPoly.close();
//  if (system("/home/cfog/Downloads/tetgen1.5.0/tetgen -pY fill.poly") != 0)
//    assert(0);
//
//  string strLine;
//
//  // Check that no Steiner points were added. We'll figure out how to handle those later.
//  ifstream fNode;
//  fNode.open("fill.1.node");
//  getline(fNode, strLine);
//  istringstream iss(strLine);
//  size_t iNumNodes;
//  iss >> iNumNodes;
//
//  if (iNumNodes != verts.size()) {
//    for (size_t ii = 0; ii < verts.size(); ii++) {
//      getline(fNode, strLine);
//    }
//
//    double vert[3];
//    int iIndex;
//    while (1) {
//      getline(fNode, strLine);
//      istringstream ISS(strLine);
//
//      string first;
//      ISS >> first;
//      if (first == "#" || fNode.eof())
//	break;
//
//      istringstream ISS2(strLine);
//
//      ISS2 >> iIndex >> vert[0] >> vert[1] >> vert[2];
//      assert(ISS.eof());
//
//      edamVert* pEV = createEdamVert(vert[0], vert[1], vert[2], iLayer + 1);
//      pEV->setSteiner(true);
//      //cout << "Inserted Steiner Point!!! " << iIndex <<": "<< pEV->x()<< " " <<pEV->y()<<" "<<pEV->z()  << endl;
//
//      realToTemp.insert(pair<int, int>(pEV->getIndex(), iIndex));
//      tempToReal.insert(pair<int, int>(iIndex, pEV->getIndex()));
//      verts.push_back(pEV);
//    }
//  }
//
//  fNode.close();
//  assert(
//      iNumNodes == verts.size()
//	  && "Steiner point created! Need to add support for this!");
//
//  // Insert new tetrahedra into the actual mesh.
//  ifstream fEle;
//  fEle.open("fill.1.ele");
//  std::getline(fEle, strLine);
//
//  int tet[4];
//  while (1) {
//    string line;
//    std::getline(fEle, line);
//    istringstream ISS(line);
//
//    string first;
//    ISS >> first;
//    if (first == "#")
//      break;
//
//    ISS >> tet[0] >> tet[1] >> tet[2] >> tet[3];
//    assert(ISS.eof());
//
//    if (tempToReal.find(tet[0]) == tempToReal.end()) {
//      cout << "index " << tet[0] << endl;
//    }
//    if (tempToReal.find(tet[1]) == tempToReal.end()) {
//      cout << "index " << tet[1] << endl;
//    }
//    if (tempToReal.find(tet[2]) == tempToReal.end()) {
//      cout << "index " << tet[2] << endl;
//    }
//    if (tempToReal.find(tet[3]) == tempToReal.end()) {
//      cout << "index " << tet[3] << endl;
//    }
//    edamVert* pEV0 = m_eVL.get(tempToReal.at(tet[0]));
//    edamVert* pEV1 = m_eVL.get(tempToReal.at(tet[1]));
//    edamVert* pEV2 = m_eVL.get(tempToReal.at(tet[2]));
//    edamVert* pEV3 = m_eVL.get(tempToReal.at(tet[3]));
//
//    bool existed = false;
//    createTetCell(existed, pEV0->getVert(), pEV1->getVert(), pEV2->getVert(),
//		  pEV3->getVert());
//    assert(!existed);
//  }
//}
/*
 * Triangulate a polygonal hole in the front. Used mostly for decimation.
 *
 * This finds the weight-minimizing triangulation (currently weight = area) by searching through all the combinatorial possibilities,
 * so runtime is exponential in the number of edges...
 *
 * If no triangulation exists that preserves manifoldness of the front, return a very large weight (> 1e200) to indicate this.
 */
double
triangulateHole(vector<pair<edamVert*, edamVert*>> edges,
		vector<vector<edamVert*>> &newTris)
{

  // Weight of triangulation.
  double dWeight = LARGE_DBL;

  // Create a triangle connecting the first edge (arbitrary) to some other vertex (try all vertices to find best).
  // This splits the domain into at most two new polygons.
  // Pass those polygons recursively to this function to obtain the total triangulation weight.
  for (size_t ii = 2; ii < edges.size(); ii++) {
    double dWeightCurr = 0;

    do {
      edamVert* pEV0 = edges.at(0).first;
      edamVert* pEV1 = edges.at(0).second;
      edamVert* pEV2 = edges.at(ii).first;

      // Calculate the area of the triangle.
      double ad01[3] =
	{ pEV1->x() - pEV0->x(), pEV1->y() - pEV0->y(), pEV1->z() - pEV0->z() };
      double ad02[3] =
	{ pEV2->x() - pEV0->x(), pEV2->y() - pEV0->y(), pEV2->z() - pEV0->z() };
      double adCross[3];
      vCROSS3D(ad01, ad02, adCross);
      dWeightCurr += dMAG3D(adCross);

      // Duplicate edge produces non-manifold.
      if (ii != 2 && pEV1->hasCommonEdge(pEV2))
	dWeightCurr += 1.e200;
      if (ii != edges.size() - 1 && pEV0->hasCommonEdge(pEV2))
	dWeightCurr += 1.e200;
    }
    while (0);

    // Find left and right subdomains and triangulate those.
    vector<vector<edamVert*>> newTrisLCurr;
    vector<vector<edamVert*>> newTrisRCurr;

    if (ii != edges.size() - 1) {
      pair<edamVert*, edamVert*> newEdgeL(edges.at(0).first,
					  edges.at(ii).first);
      vector<pair<edamVert*, edamVert*>> edgesL;
      edgesL.push_back(newEdgeL);
      for (size_t iL = ii; iL < edges.size(); iL++) {
	edgesL.push_back(edges.at(iL));
      }
      dWeightCurr += triangulateHole(edgesL, newTrisLCurr);
    }

    if (ii != 2) {
      pair<edamVert*, edamVert*> newEdgeR(edges.at(ii).first,
					  edges.at(0).second);
      vector<pair<edamVert*, edamVert*>> edgesR;
      edgesR.push_back(newEdgeR);
      for (size_t iR = 1; iR < ii; iR++) {
	edgesR.push_back(edges.at(iR));
      }
      dWeightCurr += triangulateHole(edgesR, newTrisRCurr);
    }

    // Keep track of the triangulation with smallest weight so far.
    if (dWeightCurr < dWeight) {
      dWeight = dWeightCurr;
      newTris = newTrisLCurr;
      newTris.insert(newTris.end(), newTrisRCurr.begin(), newTrisRCurr.end());
      newTris.push_back(vector<edamVert*>(
	{ edges.at(0).first, edges.at(0).second, edges.at(ii).first }));
    }
  }

  return dWeight;
}

/*
 * Remove a vertex from the front and triangulate the resulting hole.
 * Return:
 * 		0, if decimation fails.
 * 		1, if decimation succeeds as normal.
 * 		2, if decimating a tetrahedral (4-vertex) front.
 */
int
edamMesh3D_Closed::decimateVert(edamVert* pEV)
{

  // Find adjacent faces and order them CCW around vertex.
  set<edamFace*> adjFaces;
  edamFace* pEFstart = pEV->getFace(0);
  edamFace* pEF = pEFstart;
  do {
    adjFaces.insert(pEF);

    edamVert* pEV0 = pEF->getNextVert(pEV);
    edamVert* pEV1 = pEF->getNextVert(pEV0);

    pEF = pEV->getOtherCommonFace(pEV1, pEV0);
  }
  while (pEF != pEFstart);

  // Make sure we iterated over all adjacent faces. Otherwise, we have a non-manifold or some other problem.
  if (static_cast<int>(adjFaces.size()) != pEV->getNumFaces()) {
    writeFaceListToVTK(adjFaces, "assertAdj.vtk");
    writeFaceListToVTK(pEV->getFaceList(), "assertFaceList.vtk");

    assert(static_cast<int>(adjFaces.size()) == pEV->getNumFaces());
  }

  // Check if the front being decimated is made up of only 4 vertices (so it's a tetrahedron). If so, delete the whole thing.
  if (adjFaces.size() == 3) {

    vector<edamVert*> verts;
    for (edamFace *pTF : adjFaces) {
      verts.push_back(pTF->getNextVert(pEV));
    }
    verts.push_back(pEV);

    edamFace* pEFbase = verts[0]->getCommonFace(verts[1], verts[2]);
    if (pEFbase) {
      cout << "Decimating a tetrahedron." << endl;

      deleteTriFace(pEFbase);
      for (edamFace* pEFdelete : adjFaces) {
	deleteTriFace(pEFdelete);
      }
      for (edamVert* pEVadjdelete : verts) {
	deleteEdamVert(pEVadjdelete);
      }
      return 2;
    }
  }

  // In the usual case, find the polygonal boundary of the vertex neighborhood and retriangulate without the decimated vertex.
  vector<pair<edamVert*, edamVert*>> edges;
  for (edamFace* pEFadj : adjFaces) {
    edamVert* pEV0 = pEFadj->getNextVert(pEV);
    edamVert* pEV1 = pEFadj->getNextVert(pEV0);
    edges.push_back(pair<edamVert*, edamVert*>(pEV0, pEV1));
  }
  vector<vector<edamVert*>> newTris;

  double dWeight = triangulateHole(edges, newTris);
  if (dWeight > 1.e199) {
    // No valid triangulation after decimating.
    return 0;
  }

  // Delete the vertex and adjacent faces.
  for (edamFace* pEFadj : adjFaces) {
    deleteTriFace(pEFadj);
  }
  deleteEdamVert(pEV);

  // Create the new faces in the triangulation.
  for (vector<edamVert*> tri : newTris) {
	  edamVert* vert0 = tri[0];
	  edamVert* vert1 = tri[1];
	  edamVert* vert2 = tri[2];
	  createTriFace(vert0, vert1, vert2);
  }
  return 1;
}

/*
 * Decimate a group of vertices.
 *
 * If it fails to decimate any one vertex, it will insert the vertex at the back of the queue and try again later,
 * so this is more robust than decimating individually. It fails only if none of the remaining vertices can be decimated.
 */
void
edamMesh3D_Closed::decimateVertList(const int iLayer, vector<edamVert*> verts)
{

  vector<edamVert*> prevVerts = verts;
  size_t iSize;
  do {
    iSize = prevVerts.size();

    vector<edamVert*> nextVerts;

    for (edamVert* pEV : prevVerts) {
      if (!pEV->getVert())
	continue;

      int iDeci = decimateVert(pEV);

      // Re-queue the vertex if decimation fails.
      if (iDeci == 0) {

	nextVerts.push_back(pEV);
      }
    }

    prevVerts = nextVerts;

    // If any vertices failed to decimate, loop back and try again.
  }
  while (prevVerts.size() < iSize && prevVerts.size() != 0);

  // Crash if decimation fails.
  if (prevVerts.size() != 0) {
    cout << "Decimation failed. Remaining m_verts not decimated:"
	<< prevVerts.size() << endl;
    for (edamVert* pEV : prevVerts) {
      cout << pEV->x() << ":" << pEV->y() << ":" << pEV->z() << endl;
    }
    writeFrontToVTK(m_eVL, iLayer, "DecimationFailed.vtk");

    assert(prevVerts.size() == 0);
  }
}

/*
 * If a vertex has a parent edge outward between its two in-layer neighbors,
 * that parent edge necessarily intersects some other new front edge.
 * "Decimate" the vertex with the illegal parent edge.
 */
/*
 * Uses incremental vertex decimation to remove flipped faces from the front.
 *
 * A face is "flipped" if its normal points in a direction opposite the edge from the parent vertex to this one.
 * This is not necessarily a problem, so this can decimate vertices that don't need to be decimated.
 * TODO: Would be better to implement the analog to 2D "deleteBadVerts" if possible.
 */
//void
//edamMesh3D_Closed::removeTempVertsFlipped(const int iLayer)
//{
//
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isNoDecimate())
//      continue;
//
//    edamVert* pEVpar = pEV->getParent();
//    assert(pEVpar);
//    double adDistPar[3] =
//      { pEV->x() - pEVpar->x(), pEV->y() - pEVpar->y(), pEV->z() - pEVpar->z() };
//
//    vector<edamFace*> adjFaces = pEV->getFaceList();
//    for (edamFace* pEF : adjFaces) {
//      double adNorm[3] =
//	{ 0, 0, 0 };
//      pEF->calcUnitNormal(adNorm);
//      if (dDOT3D(adNorm,adDistPar) > 0) {
//	int iDeci = decimateVert(pEV);
//	assert(iDeci == 1);
//
//	break;
//      }
//    }
//  }
//}
//static double
//calcMinDistToLineOpposite(const edamVert* const pEV, const edamFace* const pEF)
//{
//  const edamVert* pEV1 = pEF->getNextVert(pEV);
//  const edamVert* pEV2 = pEF->getNextVert(pEV1);
//
//  double vecTo1[3] =
//    { pEV1->x() - pEV->x(), pEV1->y() - pEV->y(), pEV1->z() - pEV->z() };
//  double vecTo2[3] =
//    { pEV2->x() - pEV->x(), pEV2->y() - pEV->y(), pEV2->z() - pEV->z() };
//
//  double vec1To2[3] = adDIFF3D(vecTo1,vecTo2);
//  double dProj = dDOT3D(vec1To2, vecTo1) / (dMAG3D(vec1To2) * dMAG3D(vec1To2));
//  double vecShortestToLine[3];
//  for (int ii = 0; ii < 3; ii++) {
//    vecShortestToLine[ii] = vecTo1[ii] - vec1To2[ii] * dProj;
//  }
//  return dMAG3D(vecShortestToLine);
//}
/*
 * Decimate vertices too close to others on the new front.
 */
//void
//edamMesh3D_Closed::removeTempVertsForProximity(const int iLayer)
//{
//  size_t layerStart = m_eVL.layerStartIndex(iLayer + 1);
//  size_t layerEnd = m_eVL.size() - 1;
//  size_t layerSize = layerEnd - layerStart + 1;
//  bool isTooClose[layerSize];
//  for (size_t ii = 0; ii < layerSize; ii++) {
//    isTooClose[ii] = false;
//  }
//  for (size_t iVert = layerEnd; iVert >= layerStart; iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isNoDecimate())
//      continue;
//
//    vector<edamFace*> adjFaces = pEV->getFaceList();
//    // Do a quick check to see if any of the neighboring verts are already
//    // marked for decimation.  If so, don't even check this one.
//    bool haveDeletedNeigh = false;
//    for (edamFace* pEF : adjFaces) {
//      edamVert* pENext = pEF->getNextVert(pEV);
//      int iNeighVert = pENext->getIndex();
//      if (isTooClose[iNeighVert - layerStart]) {
//	haveDeletedNeigh = true;
//	break;
//      }
//    }
//    for (edamFace* pEF : adjFaces) {
//      double dMinDist = calcMinDistToLineOpposite(pEV, pEF);
//      double scaling = haveDeletedNeigh ? 0.7 : 0.7;
//      if (dMinDist < scaling * pEV->getDist(pEV->getParent())) {
//	isTooClose[iVert - layerStart] = true;
//	break;
//      }
//    }
//  }
//  for (size_t iVert = layerEnd; iVert >= layerStart; iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (isTooClose[iVert - layerStart]) {
//      int iDeci = decimateVert(pEV);
//      assert(iDeci == 1);
//    }
//  }
//}
/*
 * Decimate vertices introduced by clipping.
 */
//void
//edamMesh3D_Closed::removeTempVertsFromClipping(const int iLayer)
//{
//
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isNoDecimate())
//      continue;
//
//    if (pEV->isIntersectClip2()) {
//      int iDeci = decimateVert(pEV);
//      assert(iDeci == 1);
//    }
//  }
//
//  // Decimate vertices of triangles that were split during clipping.
//  // This is an overkill solution to a problem that came up occasionally:
//  // clipping a triangle but not removing any of its vertices, leaving impossible-to-close-with-pyramids quad faces underneath.
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isIntersectAdj()) {
//      decimateVert(pEV);
//    }
//  }
////  VTKFront(m_eVL, iLayer, "intersectionsDecimated.vtk");
//}
//void
//edamMesh3D_Closed::removeTempVerts(int iLayer)
//{
//  vector<edamVert*> tempVerts;
//
//  for (int iVert = m_eVL.size() - 1; iVert >= m_eVL.layerStartIndex(iLayer + 1);
//      iVert--) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isNoDecimate())
//      continue;
//
//    if (pEV->isIntersectClip2()) {
//      tempVerts.push_back(pEV);
//      continue;
//    }
//
//    vector<edamFace*> adjFaces = pEV->getFaceList();
//    edamVert* pEVpar = pEV->getParent();
//    assert(pEVpar);
//    double adDistPar[3] =
//      { pEV->x() - pEVpar->x(), pEV->y() - pEVpar->y(), pEV->z() - pEVpar->z() };
//
//    double dMinDist = LARGE_DBL;
//    for (edamFace* pEF : adjFaces) {
//      edamVert* pEV0 = pEF->getNextVert(pEV);
//      if (!pEV0->isIntersectClip2()) {
//	dMinDist = fmin(dMinDist, pEV->getDist(pEV0));
//      }
//    }
//    if (dMinDist < 0.7 * dMAG3D(adDistPar)) {
//      tempVerts.push_back(pEV);
//      continue;
//    }
//
//    bool qBad = false;
//    for (edamFace* pEF : adjFaces) {
//      double adNorm[3] =
//	{ 0, 0, 0 };
//      pEF->calcUnitNormal(adNorm);
//      if (dDOT3D(adNorm,adDistPar) > 0) {
//	qBad = true;
//	break;
//      }
//    }
//    if (qBad) {
//      tempVerts.push_back(pEV);
//      continue;
//    }
//  }
//
//  decimateVertList(iLayer, tempVerts);
//
//  set<edamFace*, bool
//  (*)(const edamFace*, const edamFace*)> topFaceSet2(
//      edamFace::qPointerCompFaceAddress);
//  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//      iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//
//    for (int iF = 0; iF < pEV->getNumFaces(); iF++) {
//      edamFace* pEF = pEV->getFace(iF);
//      if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer + 1
//	  && pEF->getVert(1)->getLayer() == iLayer + 1
//	  && pEF->getVert(2)->getLayer() == iLayer + 1) {
//	topFaceSet2.insert(pEF);
//      }
//    }
//  }
//  vector<edamFace*> segmentList2;
//  segmentList2.insert(segmentList2.end(), topFaceSet2.begin(),
//		      topFaceSet2.end());
//  VTKFaceList(segmentList2, "intersectionsDecimated.vtk");
//}
//void
//edamMesh3D_Closed::closeTunnels(const int iLayer)
//{
//
//  // Generate a .off mesh file from the current layer.
//  vector<edamVert*> vertList;
//  map<int, int> tempToReal;
//  do {
//    set<edamFace*, bool
//    (*)(const edamFace*, const edamFace*)> topFaceSet(
//	edamFace::qPointerCompFaceAddress);
//    for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//	iVert++) {
//      edamVert* pEV = m_eVL.get(iVert);
//
//      for (int iF = 0; iF < pEV->getNumFaces(); iF++) {
//	edamFace* pEF = pEV->getFace(iF);
//	if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer + 1
//	    && pEF->getVert(1)->getLayer() == iLayer + 1
//	    && pEF->getVert(2)->getLayer() == iLayer + 1) {
//	  topFaceSet.insert(pEF);
//	}
//      }
//    }
//    vector<edamFace*> segmentList;
//    segmentList.insert(segmentList.end(), topFaceSet.begin(), topFaceSet.end());
//    VTKFaceList(segmentList, "OpenTunnels.vtk");
//    OFFFaceList(segmentList, "tunnels.off", vertList, tempToReal);
//  }
//  while (0);
//
//  // Use ReebHanTun to find closed loops of vertices around tunnels.
//  if (system("rm loops_tunnels*") != 0)
//    assert(0);
//  if (system(
//      "/home/cfog/Downloads/ReebHan/ReebHanTun-Linux/ReebHanTun -I tunnels.off -O tunnels > /dev/null")
//      != 0)
//    assert(0);
//
//  // File IO... process output of ReebHanTun.
//  ifstream fLoops;
//  fLoops.open("loops_tunnels.lop");
//  if (!fLoops.good())
//    return;
//
//  vector<vector<double>> tunnelLoops;
//  string strLine;
//  string strBegin;
//  string strcsv;
//  vector<vector<int>> loops;
//
//  while (!fLoops.eof()) {
//    std::getline(fLoops, strLine);
//    istringstream iss(strLine);
//
//    iss >> strBegin;
//    if (strBegin != "tunnel")
//      continue;
//
//    iss >> strBegin >> strBegin >> strcsv;
//
//    istringstream ss(strcsv);
//    vector<int> vals;
//
//    while (!ss.eof()) {
//      string sub;
//      std::getline(ss, sub, ',');
//      istringstream ssub(sub);
//      int val;
//      ssub >> val;
//      vals.push_back(val);
//    }
//    vals.pop_back();
//    loops.push_back(vals);
//  }
//
//  // The loop vertices are, *obnoxiously*, ordered by index rather than ordered around the loop. Fix that.
//  vector<vector<edamVert*>> loopsGood;
//  for (vector<int> loop : loops) {
//    set<edamVert*> loopSet;
//    for (int iV : loop) {
//      loopSet.insert(m_eVL.get(tempToReal.at(iV)));
//    }
//    vector<edamVert*> loopGood;
//
//    edamVert* pEVstart = m_eVL.get(tempToReal.at(loop.at(0)));
//    edamVert* pEV = pEVstart;
//    edamVert* pEVprev = NULL;
//
//    do {
//      loopGood.push_back(pEV);
//      vector<edamFace*> adjFaces = pEV->getFaceList();
//      for (edamFace* pEF : adjFaces) {
//
//	edamVert* pEV0 = pEF->getNextVert(pEV);
//	if (pEV0 != pEVprev && loopSet.find(pEV0) != loopSet.end()) {
//	  pEVprev = pEV;
//	  pEV = pEV0;
//	  break;
//	}
//      }
//    }
//    while (pEV != pEVstart);
//    loopsGood.push_back(loopGood);
//  }
//
//  // For each tunnel loop, identify the loop of vertices on each side of it and close them by triangulation.
//  for (vector<edamVert*> loopGood : loopsGood) {
//    if (loopGood.size() > 10)
//      continue;
//    set<edamVert*> loopSet(loopGood.begin(), loopGood.end());
//
//    set<edamVert*> adjLoopSet;
//    map<edamVert*, edamVert*> adjLoopMap;
//
//    for (edamVert* pEV : loopGood) {
//      vector<edamFace*> adjFaces = pEV->getFaceList();
//      for (edamFace* pEF : adjFaces) {
//
//	edamVert* pEV0 = pEF->getNextVert(pEV);
//	edamVert* pEV1 = pEF->getNextVert(pEV0);
//	if (loopSet.find(pEV0) == loopSet.end()
//	    && loopSet.find(pEV1) == loopSet.end()) {
//	  adjLoopSet.insert(pEV0);
//	  adjLoopMap.insert(pair<edamVert*, edamVert*>(pEV0, pEV1));
//	}
//      }
//    }
//
//    edamVert* pEVstart = NULL;
//    edamVert* pEV = NULL;
//    edamVert* pEVnext = NULL;
//    vector<pair<edamVert*, edamVert*>> edges0;
//    pEVstart = *(adjLoopSet.begin());
//    pEV = pEVstart;
//    do {
//      pEVnext = adjLoopMap.at(pEV);
//      edges0.push_back(pair<edamVert*, edamVert*>(pEV, pEVnext));
//      adjLoopSet.erase(pEV);
//      pEV = pEVnext;
//    }
//    while (pEV != pEVstart);
//    vector<pair<edamVert*, edamVert*>> edges1;
//    pEVstart = *(adjLoopSet.begin());
//    pEV = pEVstart;
//    do {
//      pEVnext = adjLoopMap.at(pEV);
//      edges1.push_back(pair<edamVert*, edamVert*>(pEV, pEVnext));
//      adjLoopSet.erase(pEV);
//      pEV = pEVnext;
//    }
//    while (pEV != pEVstart);
//
//    vector<vector<edamVert*>> newTris0;
//    double dWeight = triangulateHole(edges0, newTris0);
//    assert(dWeight < 1.e199);
//    for (vector<edamVert*> tri : newTris0) {
//      createTriFace(tri[0], tri[1], tri[2]);
//    }
//    vector<vector<edamVert*>> newTris1;
//    dWeight = triangulateHole(edges1, newTris1);
//    assert(dWeight < 1.e199);
//    for (vector<edamVert*> tri : newTris1) {
//      createTriFace(tri[0], tri[1], tri[2]);
//    }
//
//    for (edamVert* pEVdelete : loopGood) {
//      vector<edamFace*> adjFaces = pEVdelete->getFaceList();
//      for (edamFace* pEF : adjFaces) {
//	deleteTriFace(pEF);
//      }
//      deleteEdamVert(pEVdelete);
//    }
//  }
//
////  VTKFront(m_eVL, iLayer, "ClosedTunnels.vtk");
//}
/*
 * Find gaps between previous layer and next layer where no prisms were placed. Fill these gaps with pyramids and tets.
 */
//void
//edamMesh3D_Closed::fillHoles(const int iLayer)
//{
//
//  // Store faces that make up the hull of a region we want to tetrahedralize.
//  // These faces will be fed as input to TetGen.
//  vector<edamFace*> bottomFaces;
//  vector<edamFace*> topFaces;
//
//  set<edamFace*, bool
//  (*)(const edamFace*, const edamFace*)> bottomFaceSet(
//      edamFace::qPointerCompFaceAddress);
//  for (int iVert = m_eVL.layerStartIndex(iLayer);
//      iVert < m_eVL.layerStartIndex(iLayer + 1); iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isSteiner())
//      continue;
//
//    vector<edamFace*> adjFaces = pEV->getAdjFaces();
//
//    for (edamFace* pEF : adjFaces) {
//      if (pEF->getNumCells() == 2)
//	continue;
//      if (pEF->getKid(pEV) && !pEV->isNoQuadFaces()) {
//	continue;
//      }
//      edamVert* pEV0 = pEF->getNextVert(pEV);
//      edamVert* pEV1 = pEF->getNextVert(pEV0);
//
//      edamVert* pEV2 = NULL;
//      edamVert* pEV3 = NULL;
//      bool q2 = ((pEV2 = pEF->getKid(pEV0)) && !pEV0->isNoQuadFaces());
//      bool q3 = ((pEV3 = pEF->getKid(pEV1)) && !pEV1->isNoQuadFaces());
//
//      if (q2 && q3) {
//	// Create a pyramid to close off a quad face on the edge of the hole.
//	assert(!pEV->getCommonFace(pEV2, pEV3));
//	bottomFaceSet.insert(createTriFace(pEV, pEV2, pEV3));
//	bool qExist;
//	createPyrCell(qExist, pEV0->getVert(), pEV2->getVert(), pEV3->getVert(),
//		      pEV1->getVert(), pEV->getVert());
//	assert(!qExist);
//      }
//      else if (q2) {
//	// Place tets between the pyramids where appropriate.
//	// Not actually necessary and perhaps detrimental, since the Delaunay routine would
//	// mesh these areas anyway.
//	assert(!pEV->getCommonFace(pEV2, pEV1));
//	bottomFaceSet.insert(createTriFace(pEV, pEV2, pEV1));
//	bool existed = false;
//	createTetCell(existed, pEV->getVert(), pEV1->getVert(), pEV0->getVert(),
//		      pEV2->getVert());
//	assert(!existed);
//      }
//      else if (q3) {
//	// Do nothing. The tet we could create here will be created anyway in the if(q2) block.
//      }
//      else {
//	bottomFaceSet.insert(pEF);
//      }
//    }
//  }
//  bottomFaces.insert(bottomFaces.end(), bottomFaceSet.begin(),
//		     bottomFaceSet.end());
//
//  // Top faces are all faces on the new front that have no adjacent cells.
//  set<edamFace*, bool
//  (*)(const edamFace*, const edamFace*)> topFaceSet(
//      edamFace::qPointerCompFaceAddress);
//  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//      iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//
//    for (int iF = 0; iF < pEV->getNumFaces(); iF++) {
//      edamFace* pEF = pEV->getFace(iF);
//      if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer + 1
//	  && pEF->getVert(1)->getLayer() == iLayer + 1
//	  && pEF->getVert(2)->getLayer() == iLayer + 1) {
//	if (pEF->getNumCells() == 0) {
//	  topFaceSet.insert(pEF);
//	}
//      }
//    }
//  }
//  topFaces.insert(topFaces.end(), topFaceSet.begin(), topFaceSet.end());
//
//  // Pass the hull to TetGen.
//  if (!bottomFaces.empty()) {
//    vector<edamFace*> VTKFaces;
//    VTKFaces.insert(VTKFaces.end(), bottomFaces.begin(), bottomFaces.end());
//    //VTKFaceList(VTKFaces, "delaney.vtk");
//    VTKFaces.insert(VTKFaces.end(), topFaces.begin(), topFaces.end());
////    VTKFaceList(VTKFaces, "delaney.vtk");
//    assert(!topFaces.empty());
//
//    TetGenFill(iLayer, bottomFaces, topFaces);
//  }
//}
/*
 * Remove self-intersecting regions in the front.
 */
//void
//edamMesh3D_Closed::fixIntersections(const int iLayer)
//{
//
//  // Find all front faces.
//  set<edamFace*, bool
//  (*)(const edamFace*, const edamFace*)> topFaceSet(
//      edamFace::qPointerCompFaceAddress);
//  for (size_t iVert = m_eVL.layerStartIndex(iLayer + 1); iVert < m_eVL.size();
//      iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//
//    for (int iF = 0; iF < pEV->getNumFaces(); iF++) {
//      edamFace* pEF = pEV->getFace(iF);
//      if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer + 1
//	  && pEF->getVert(1)->getLayer() == iLayer + 1
//	  && pEF->getVert(2)->getLayer() == iLayer + 1) {
//	topFaceSet.insert(pEF);
//      }
//    }
//  }
//  vector<edamFace*> segmentList;
//  segmentList.insert(segmentList.end(), topFaceSet.begin(), topFaceSet.end());
//
//  // Check for intersections of front faces. If none exist, we're done.
//  if (!clipIntersections(iLayer, segmentList))
//    return;
//  //VTKFront(m_eVL,iLayer,"intersectionsClipped.vtk");
//
//  // Duplicate non-manifold vertices after clipping, leaving several disconnected closed surfaces.
//  int iNumLoops;
//  reconnectIntersections(iLayer, iNumLoops);
//  //VTKFront(m_eVL,iLayer,"intersectionsLabeled.vtk");
//
//  // Some of these closed surfaces will be inverted, find and remove them using winding number computation.
//  findInvertedLoops(iLayer, iNumLoops);
//  //VTKFront(m_eVL,iLayer,"intersectionsInverted.vtk");
//}
/*
 * Create prisms between parent triangles and still-intact kid triangles.
 */
//void
//edamMesh3D_Closed::createGoodCells(const int iLayer)
//{
//
//  for (int iVert = m_eVL.layerStartIndex(iLayer);
//      iVert < m_eVL.layerStartIndex(iLayer + 1); iVert++) {
//    edamVert* pEV = m_eVL.get(iVert);
//    if (pEV->isSteiner())
//      continue;
//
//    vector<edamFace*> adjFaces = pEV->getAdjFaces();
//
//    int iNumFaces = adjFaces.size();
//    if (iNumFaces < 1)
//      continue;
//
//    edamVert* pEV0 = NULL;
//    edamVert* pEV1 = NULL;
//    edamVert* pEV2 = NULL;
//    edamVert* pEV3 = NULL;
//    edamVert* pEV4 = NULL;
//    edamVert* pEV5 = NULL;
//    for (edamFace* pEF : adjFaces) {
//      if (pEF->getNumCells() == 2)
//	continue;
//
//      pEV0 = pEF->getVert(0);
//      pEV1 = pEF->getVert(1);
//      pEV2 = pEF->getVert(2);
//      if (pEV0->isNoQuadFaces() || pEV1->isNoQuadFaces()
//	  || pEV2->isNoQuadFaces())
//	continue;
//
//      // Check for the kid face. If it exists, create the prism.
//      if ((pEV3 = pEF->getKid(pEV0)) && (pEV4 = pEF->getKid(pEV1)) && (pEV5 =
//	  pEF->getKid(pEV2)) && pEV3->getCommonFace(pEV4, pEV5)
//	  && pEV3->getCommonFace(pEV4, pEV5)->getNumCells() == 0) {
//	bool qExist = false;
//	createPrismCell(qExist, pEV0->getVert(), pEV1->getVert(),
//			pEV2->getVert(), pEV3->getVert(), pEV4->getVert(),
//			pEV5->getVert());
//	assert(!qExist);
//      }
//    }
//  }
//}
/*
 * Identify upcoming collisions for the front with existing mesh entities.
 *
 * For now, identify isotropic mesh points that are within marching distance
 * and remove them.  Later, this mechanism will also be used to handle front
 * collisions.
 */

void
edamMesh3D_Closed::checkForCollisions(const int layer)
{
  // For every face in the front, check that the vert opposite that face and
  // outside the front is sufficiently far away (normal distance at least
  // 3x marching step size: this will leave a cushion between fronts that
  // are about to collide.

  // If that point is closer and it's an isotropic point, then remove it.

  // First, construct a list of all of the faces in the current front, by
  // concatentating all of the vertex neighbor lists.
  std::set<edamFace*> frontFaces;
  for (size_t ii = m_eVL.layerStartIndex(layer); ii < m_eVL.size(); ii++) {
    edamVert *pEV = m_eVL.get(ii);
    if (pEV->getNumFaces() == 0)
      continue;
    std::set<edamFace*> myFrontFaces = pEV->getAdjFaces();
    frontFaces.insert(myFrontFaces.begin(), myFrontFaces.end());
  }
  for (edamFace* pEF : frontFaces) {
    Face *pF = pEF->getFace();
    Vert *pV = nullptr;
    Cell *pC = pF->getLeftCell();
    assert(pC->isValid() && !pC->isDeleted());
    if (pC->getType() == CellSkel::eTet) {
      pV = dynamic_cast<TetCell*>(pC)->getOppositeVert(pF);
    }
    if (!pV->isValid() || pV->isStructured()) {
      pC = pF->getRightCell();
      assert(
	  pC->isValid() && !pC->isDeleted() && pC->getType() == CellSkel::eTet);
      pV = dynamic_cast<TetCell*>(pC)->getOppositeVert(pF);
    }
    assert(pV->isValid() && !pV->isStructured());
    if (pV->isBdryVert())
      continue;
    double norm[3] =
      { 0 };
    pF->calcUnitNormal(norm);
    double disp[] = adDIFF3D(pV->getCoords(), pF->getVert(0)->getCoords());
    double dot = fabs(dDOT3D(disp, norm));
    double marchDist = calcMarchingDistance(layer);
    if (dot < marchDist * 3) {
      logMessage(MSG_MANAGER, "Deleting an iso vertex.\n");
      int swaps;
      removeVert(pV, swaps);
    }
  }
}

bool
edamMesh3D_Closed::addCellAdjacentToThisFace(const Face* const pF,
					     const edamVert* const pEVCand,
					     GRUMMP::CavityInserter3D& CI) const
{
  assert(pF->getLeftCell()->isValid() && pF->getRightCell()->isValid());
  int orient = checkOrient3D(pF->getVert(0), pF->getVert(1), pF->getVert(2),
			     pEVCand->getVert());
  Cell* pC = pF->getLeftCell();
  if (pC->getType() == CellSkel::eTet) {
    double centroid[3];
    pC->calcCentroid(centroid);
    int orientTmp = checkOrient3D(pF->getVert(0)->getCoords(),
				  pF->getVert(1)->getCoords(),
				  pF->getVert(2)->getCoords(), centroid);
    if (orientTmp == orient) {
      CI.addCellToCavity(pC);
      return true;
    }
  }

  // Left cell didn't do it; try the right one.
  pC = pF->getRightCell();
  if (pC->getType() == CellSkel::eTet) {
    double centroid[3];
    pC->calcCentroid(centroid);
    int orientTmp = checkOrient3D(pF->getVert(0)->getCoords(),
				  pF->getVert(1)->getCoords(),
				  pF->getVert(2)->getCoords(), centroid);
    if (orientTmp == orient) {
      CI.addCellToCavity(pC);
      return true;
    }
  }
  return false;
}

/*
 *    Need to build a hull around the proposed vertex.  This hull must include
 *    all the faces around the original vert pEV.  Also, any lateral faces that
 *    connect edges of that star to their descendant edges have to be on the
 *    surface of the hull.  None of those faces can be crossed in building out
 *    the hull.  Also, going full Delaunay is going to be a bad idea, so including
 *    any other front faces in the hull is probably going to be a bad idea.
 */
void
edamMesh3D_Closed::computeInsertionCavity(GRUMMP::CavityInserter3D& CI,
					  edamVert* pEV, edamVert* pEVCand)
{
  set<edamFace*> adjFaces = pEV->getAdjFaces();

  CI.clear();

  logMessage(MSG_DEBUG, "New cavity\n");
  logMessage(MSG_DEBUG, "EV index: %d\n", pEV->getIndex());
  logMessage(MSG_DEBUG, "Cand index: %d\n", pEVCand->getIndex());

// Add faces built on the triangles in the previous level of the front and
// incident on the vert being marched from (pEV).
  for (edamFace* pEF : adjFaces) {
    Face* pF = pEF->getFace();

    bool ok = addCellAdjacentToThisFace(pF, pEVCand, CI);
    assert(ok);

    // Now find faces that are extruded from perimeter edges around pEV.
    // For each perimeter edge (that is, each of the edges in an adjFace that's
    // opposite to pEV) there could be zero, one, or two of these, depending on
    // the number of child verts present.

    // Identify the two other verts.
    edamVert *pEVA = nullptr, *pEVB = nullptr;
    if (pEF->getVert(0) == pEV) {
      pEVA = pEF->getVert(1);
      pEVB = pEF->getVert(2);
    }
    else {
      pEVA = pEF->getVert(0);
      if (pEF->getVert(1) == pEV) {
	pEVB = pEF->getVert(2);
      }
      else {
	pEVB = pEF->getVert(1);
      }
    }
    assert(pEVA && pEV != pEVA);
    assert(pEVB && pEV != pEVB);
    assert(pEVA != pEVB);
    assert(pEVA->getNumKids() <= 1);
    assert(pEVB->getNumKids() <= 1);
    logMessage(MSG_DEBUG, "\nEVA index: %d  Kids: %d\n", pEVA->getIndex(),
	       pEVA->getNumKids());
    logMessage(MSG_DEBUG, "EVB index: %d  Kids: %d\n", pEVB->getIndex(),
	       pEVB->getNumKids());

    if (pEVA->getNumKids() == 1) {
      edamVert *pEVAkid = pEVA->getFirstKid();
      logMessage(MSG_DEBUG, "EVAkid index: %d\n", pEVAkid->getIndex());
      if (pEVB->getNumKids() == 1) {
	// Both verts have kids; looking for two faces.
	edamVert *pEVBkid = pEVB->getFirstKid();
	logMessage(MSG_DEBUG, "EVBkid index: %d\n", pEVBkid->getIndex());
	Vert *pVA = pEVA->getVert();
	Vert *pVB = pEVB->getVert();
	Vert *pVAkid = pEVAkid->getVert();
	Vert *pVBkid = pEVBkid->getVert();
	Face *pFABAk = findCommonFace(pVA, pVB, pVAkid, nullptr, true);
	if (pFABAk->isValid()) {
	  // Have a face containing A, B, and Akid.  Must also have one with
	  // B, Akid, and Bkid.
	  Face *pFBAkBk = findCommonFace(pVB, pVAkid, pVBkid, nullptr, true);
	  assert(pFBAkBk);
	  ok = addCellAdjacentToThisFace(pFABAk, pEVCand, CI);
	  assert(ok);
	  ok = addCellAdjacentToThisFace(pFBAkBk, pEVCand, CI);
	  assert(ok);
	}
	else {
	  // No face containing A, B, and Akid.  Must have instead faces with
	  // A, B, Bkid and A, Bkid, Akid.
	  Face *pFABBk = findCommonFace(pVA, pVB, pVBkid, nullptr, true);
	  Face *pFABkAk = findCommonFace(pVA, pVBkid, pVAkid, nullptr, true);
	  assert(pFABBk);
	  assert(pFABkAk);
	  ok = addCellAdjacentToThisFace(pFABBk, pEVCand, CI);
	  assert(ok);
	  ok = addCellAdjacentToThisFace(pFABkAk, pEVCand, CI);
	  assert(ok);
	}
      }
      else {
	// Only A has a kid; looking for one face.
	Vert *pVA = pEVA->getVert();
	Vert *pVB = pEVB->getVert();
	Vert *pVAkid = pEVAkid->getVert();
	Face *pFABAk = findCommonFace(pVA, pVB, pVAkid, nullptr, true);
	assert(pFABAk);
	ok = addCellAdjacentToThisFace(pFABAk, pEVCand, CI);
	assert(ok);
      }
    }
    else if (pEVB->getNumKids() == 1) {
      // Already know that A has no kids; looking for one face.
      edamVert *pEVBkid = pEVB->getFirstKid();
      logMessage(MSG_DEBUG, "EVBkid index: %d\n", pEVBkid->getIndex());
      Vert *pVA = pEVA->getVert();
      Vert *pVB = pEVB->getVert();
      Vert *pVBkid = pEVBkid->getVert();
      Face *pFABBk = findCommonFace(pVA, pVB, pVBkid, nullptr, true);
      assert(pFABBk);
      ok = addCellAdjacentToThisFace(pFABBk, pEVCand, CI);
      assert(ok);
    }
  }

// Finish off the cavity by making it so that all its hull faces are properly
// oriented with respect to the new vert.
  CI.makeCavityConvex(pEVCand->getVert()->getCoords());
}

double
edamMesh3D_Closed::calcMarchingDistance(const int iLayer, double dExtFactor)
{
//double dExtFactor = 1;
  double dExt = m_firstLayerHeight * pow(m_growthRatio, iLayer) * dExtFactor;
  return dExt;
}

/*
 * Extrude an offset surface from each existing front. Create the surface triangulation for the next layer
 * but leave out the volume elements for now.
 *
 * This version does not attempt multiple normals from convex edges/corners.
 */
void
edamMesh3D_Closed::extrudeLayer(const int iLayer)
{
  GRUMMP::CavityInserter3D CI(this);

  int iNumVertsPrev = m_eVL.size();
  set<edamFace*> VTKfaces;

  for (int iVert = m_eVL.layerStartIndex(iLayer); iVert < iNumVertsPrev;
      iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner() || pEV->getNumFaces() == 0)
      continue;

    set<edamFace*> adjFaces = pEV->getAdjFaces();

    int iNumFaces = adjFaces.size();
    if (iNumFaces < 1)
      continue;

// Compute vertex normal and some additional information.
    double adUnitNorm[3];
    double adAvgDir[3];
    double adAvgDist[3];
    double dMinDist;
    neighborhoodInfo(pEV, adjFaces, adUnitNorm, adAvgDir, adAvgDist, dMinDist);

// Compute length of extrusion edge. Increases with vertex concavity, decreases with convexity
// to smooth out sharp corners.
    double dExtFactor = (1 + dDOT3D(adAvgDir, adUnitNorm));
    dExtFactor *= dExtFactor;
//double dExtFactor = 1;
    double dExt = calcMarchingDistance(iLayer, dExtFactor);

    double coords[] =
      { pEV->x() + adUnitNorm[0] * dExt, pEV->y() + adUnitNorm[1] * dExt,
	  pEV->z() + adUnitNorm[2] * dExt };
    edamVert* pEVkid = createEdamVert(coords[0], coords[1], coords[2],
				      iLayer + 1);


// Need to build a hull around the proposed vertex.  This hull must include
// all the faces around the original vert pEV.  Also, any lateral faces that
// connect edges of that star to their descendant edges have to be on the
// surface of the hull.  None of those faces can be crossed in building out
// the hull.  Also, going full Delaunay is going to be a bad idea, so including
// any other front faces in the hull is probably going to be a bad idea.
    computeInsertionCavity(CI, pEV, pEVkid);
    CI.writeCavityVTK("cavity-check");
    Vert* pVRet = CI.insertPointInHull(pEVkid->getVert());

    assert(pVRet == pEVkid->getVert());
    assert(pVRet->x() == coords[0]);
    assert(pVRet->y() == coords[1]);
    assert(pVRet->z() == coords[2]);

    pEVkid->setParent(pEV);
    pEV->addKid(pEVkid);
    pEVkid->setAspectRatio(pEV->getAspectRatio() / m_growthRatio);

// Create kid faces for adjacent faces if all necessary kid vertices have been extruded.
    for (edamFace* pEF : adjFaces) {
      pEF->setKid(pEV, pEVkid);
      edamVert* pEV3 = pEF->getKid(0);
      edamVert* pEV4 = pEF->getKid(1);
      edamVert* pEV5 = pEF->getKid(2);

      if (pEV3 && pEV4 && pEV5) {

	edamFace* pEFnew = createTriFace(pEV3, pEV4, pEV5);
	if (pEFnew->getFace()->getLeftCell() == nullptr
	    || pEFnew->getFace()->getRightCell() == nullptr) {
	  logMessage(MSG_MANAGER, "Fishy face: (%d %d %d)\n", pEV3->getIndex(),
		     pEV4->getIndex(), pEV5->getIndex());
	}

	VTKfaces.insert(pEFnew);
      }
//      else {
//	logMessage(MSG_MANAGER,
//		   "Failed to create new front face: %d, (%d %d %d)\n",
//		   pEV->getIndex(), pEF->getVert(0)->getIndex(),
//		   pEF->getVert(1)->getIndex(), pEF->getVert(2)->getIndex());
//      }
    }
    assert(isValid());
  }
#ifndef NDEBUG
  writeFaceListToVTK(VTKfaces, "nextLayer.vtk");
#endif
}

/*
 * Extrude an offset surface from each existing front. Create the surface triangulation for the next layer
 * but leave out the volume elements for now.
 *
 * This version attempts multiple normals at convex edges, does not currently work well.
 */
void
edamMesh3D_Closed::extrudeLayerMultiNormal(int iLayer)
{

  int iNumVertsPrev = m_eVL.size();
  set<edamFace*> VTKfaces;

  set<pair<edamVert*, edamVert*>> convexEdgeSet;
  vector<vector<edamVert*>> convexCorners;

  for (int iVert = m_eVL.layerStartIndex(iLayer); iVert < iNumVertsPrev;
      iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner())
      continue;

    vector<edamFace*> adjFaces;
    edamFace* pEFstart = *(pEV->getAdjFaces().begin());
    edamFace* pEF = pEFstart;
    do {
      adjFaces.push_back(pEF);

      edamVert* pEV0 = pEF->getNextVert(pEV);
      edamVert* pEV1 = pEF->getNextVert(pEV0);
      pEF = pEV->getOtherCommonFace(pEV1, pEV0);
    }
    while (pEF != pEFstart);

    assert(adjFaces.size() > 2);

// Compute vertex normal and some additional information.
    double adUnitNorm[3];
    double adAvgDir[3];
    double adAvgDist[3];
    double dMinDist;
    set<edamFace*> setAdjFaces;
    setAdjFaces.insert(adjFaces.begin(), adjFaces.end());
    neighborhoodInfo(pEV, setAdjFaces, adUnitNorm, adAvgDir, adAvgDist,
		     dMinDist);

// Compute length of extrusion edge. Increases with vertex concavity, decreases with convexity
// to smooth out sharp corners.
//double dExtFactor = (1 + dDOT3D(adAvgDir, adUnitNorm));
    double dExtFactor = 1;
    double dExt = m_firstLayerHeight * pow(m_growthRatio, iLayer) * dExtFactor;

// Find the dihedral angle formed at each adjacent edge.
    vector<double> adjAngles;
    for (size_t ii = 0; ii < adjFaces.size(); ii++) {
      edamFace* pEF0 = adjFaces[(ii - 1 + adjFaces.size()) % adjFaces.size()];
      edamFace* pEF1 = adjFaces[ii];

      double adNorm0[3];
      pEF0->calcUnitNormal(adNorm0);
      double adNorm1[3];
      pEF1->calcUnitNormal(adNorm1);
      double dAngle = acos(dDOT3D(adNorm0, adNorm1));

      const double* ad0 = pEF0->getNextVert(pEV)->getVert()->getCoords();
      const double* ad1 =
	  pEF1->getNextVert(pEF1->getNextVert(pEV))->getVert()->getCoords();
      double adDist01[3] = adDIFF3D(ad1,ad0);

      if (dDOT3D(adNorm0,adDist01) > 0) {
	adjAngles.push_back(M_PI + dAngle);
      }
      else {
	adjAngles.push_back(M_PI - dAngle);
      }
    }

// Mark edges for multiple normals if they exceed an angle threshold.
    double dAngleThreshold = M_PI * 4.0 / 3.0;
    vector<int> convexEdges;
    for (size_t iA = 0; iA < adjAngles.size(); iA++) {
      double dAngle = adjAngles.at(iA);
      if (dAngle > dAngleThreshold) {
	convexEdges.push_back(iA);

	edamVert* pEV0 = adjFaces.at(iA)->getNextVert(pEV);
	convexEdgeSet.insert(
	    pair<edamVert*, edamVert*>(min(pEV, pEV0), max(pEV, pEV0)));
      }
    }

// If any edges have been marked, create multiple kids at different normal directions.
    if (!convexEdges.empty()) {
      cout << convexEdges.size() << " convex edges at " << pEV->x() << ":"
	  << pEV->y() << ":" << pEV->z() << endl;

      // Can't create prisms, etc. adjacent to multiple-normal edges, because the two skinny triangles
      // extruded from the multiple-normal edge cannot be connected to neighboring quad faces properly.
      pEV->setNoQuadFaces(true);

      // If only one convex edge, need to treat an opposite edge as convex also.
      if (convexEdges.size() == 1) {
	int iA = convexEdges[0];
	int iA2 = (iA + adjAngles.size() / 2) % adjAngles.size();

	convexEdges.push_back(iA2);
	edamVert* pEV0 = adjFaces.at(iA2)->getNextVert(pEV);
	pEV0->setNoQuadFaces(true);
	convexEdgeSet.insert(
	    pair<edamVert*, edamVert*>(min(pEV, pEV0), max(pEV, pEV0)));
      }

      vector<edamVert*> kidVerts;
      // Extrude a new kid vertex for every adjacent pair of convex edges.
      for (size_t iCE = 0; iCE < convexEdges.size(); iCE++) {

	// Identify two adjacent convex edges.
	int iEdge0 = convexEdges.at(iCE);
	int iEdge1 = convexEdges.at((iCE + 1) % convexEdges.size());

	// Calculate a normal direction based on the adjacent faces between the two convex edges.
	edamFace* pEF0 = adjFaces.at(iEdge0);
	edamFace* pEF1 = adjFaces.at(
	    (iEdge1 - 1 + adjFaces.size()) % adjFaces.size());
	double adNorm0[3];
	pEF0->calcUnitNormal(adNorm0);
	double adNorm1[3];
	pEF1->calcUnitNormal(adNorm1);
	double adNormCE[3];
	for (int ii = 0; ii < 3; ii++) {
	  adNormCE[ii] = 0.7 * adUnitNorm[ii]
	      - 0.15 * (adNorm0[ii] + adNorm1[ii]);
	}
	NORMALIZE3D(adNormCE);

	// Create the kid vertex.
	edamVert* pEVkid = createEdamVert(pEV->x() + adNormCE[0] * dExt,
					  pEV->y() + adNormCE[1] * dExt,
					  pEV->z() + adNormCE[2] * dExt,
					  iLayer + 1);
	pEVkid->setParent(pEV);
	pEV->addKid(pEVkid);

	// Prevent decimation and smoothing of this vertex (not acceptable as a long-term solution...)
	pEVkid->setNoDecimate(true);
	pEVkid->setNoSmooth(true);
	kidVerts.push_back(pEVkid);

	// Assign parent-kid relations.
	for (int iFace = iEdge0; iFace != iEdge1;
	    iFace = (iFace + 1) % adjFaces.size()) {
	  edamFace* pEFadj = adjFaces.at(iFace);
	  pEFadj->setKid(pEV, pEVkid);
	}
      }

      // If this is a convex corner (more than 2 kid verts), we need to triangulate the hole between the kid
      // verts. Store them for later.
      if (kidVerts.size() > 2) {
	convexCorners.push_back(kidVerts);
      }
    }

// In the regular case, just extrude a single kid vertex.
    else {
      edamVert* pEVkid = createEdamVert(pEV->x() + adUnitNorm[0] * dExt,
					pEV->y() + adUnitNorm[1] * dExt,
					pEV->z() + adUnitNorm[2] * dExt,
					iLayer + 1);
      pEVkid->setParent(pEV);
      pEV->addKid(pEVkid);

      // Currently, never allow any children of a multiple-normal vertex to decimate or smooth.
      // Obviously not appropriate as a long-term solution.
      if (pEV->isNoDecimate()) {
	pEVkid->setNoDecimate(true);
      }
      if (pEV->isNoSmooth()) {
	pEVkid->setNoSmooth(true);
      }

      for (edamFace* pEFadj : adjFaces) {
	pEFadj->setKid(pEV, pEVkid);
      }
    }
  }

// Create new kid faces corresponding to parent faces.
  for (int iVert = m_eVL.layerStartIndex(iLayer); iVert < iNumVertsPrev;
      iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner())
      continue;
    set<edamFace*> adjFaces = pEV->getAdjFaces();

    edamVert* pEV3 = NULL;
    edamVert* pEV4 = NULL;
    edamVert* pEV5 = NULL;

    for (edamFace* pEF : adjFaces) {
      if (pEF->getVert(0) != pEV)
	continue;

      pEV3 = pEF->getKid(0);
      pEV4 = pEF->getKid(1);
      pEV5 = pEF->getKid(2);

      edamFace* pEFnew = createTriFace(pEV3, pEV4, pEV5);

      VTKfaces.insert(pEFnew);
    }
  }

// Create new faces connecting kid faces across multiple-normal convex edges.
  for (pair<edamVert*, edamVert*> edge : convexEdgeSet) {
    edamVert* pEV0 = edge.first;
    edamVert* pEV1 = edge.second;
    vector<edamFace*> faces;
    set<edamFace*> adjFaces = pEV0->getAdjFaces();
    for (edamFace* pEF : adjFaces) {
      if (pEF->hasVert(pEV1))
	faces.push_back(pEF);
    }

    assert(faces.size() == 2);
    edamFace* pEF0 = faces[0];
    edamFace* pEF1 = faces[1];
    if (pEF0->getNextVert(pEV1) == pEV0) {
      edamFace* pEFtemp = pEF0;
      pEF0 = pEF1;
      pEF1 = pEFtemp;
    }

    edamFace* pEFnew = NULL;
    edamVert* pEV3 = NULL;
    edamVert* pEV4 = NULL;
    edamVert* pEV5 = NULL;

    if ((pEV3 = pEF0->getKid(pEV0)) != (pEV4 = pEF1->getKid(pEV0))) {
      pEV5 = pEF1->getKid(pEV1);
      pEFnew = createTriFace(pEV3, pEV4, pEV5);
      VTKfaces.insert(pEFnew);
    }
    if ((pEV3 = pEF1->getKid(pEV1)) != (pEV4 = pEF0->getKid(pEV1))) {
      pEV5 = pEF0->getKid(pEV0);
      pEFnew = createTriFace(pEV3, pEV4, pEV5);
      VTKfaces.insert(pEFnew);
    }
  }

// Triangulate holes at convex corners.
  for (vector<edamVert*> corner : convexCorners) {
    vector<pair<edamVert*, edamVert*>> edges;
    for (size_t iV = 0; iV < corner.size(); iV++) {
      edges.push_back(
	  pair<edamVert*, edamVert*>(corner[iV],
				     corner[(iV + 1) % corner.size()]));
    }
    vector<vector<edamVert*>> newTris;
    double dWeight = triangulateHole(edges, newTris);
    assert(dWeight < 1.e199);
    for (vector<edamVert*> tri : newTris) {
    	edamVert* vert0 = tri[0];
    	edamVert* vert1 = tri[1];
    	edamVert* vert2 = tri[2];
      edamFace* pEFnew = createTriFace(vert0, vert1, vert2);
      VTKfaces.insert(pEFnew);
    }
  }

#ifndef NDEBUG
  writeFaceListToVTK(VTKfaces, "nextLayer.vtk");
#endif
}

/*
 * Compute "forces" acting on a vertex based on neighborhood vertices.
 */
void
edamMesh3D_Closed::forceModelSmooth(int iStartLayer, int iEndLayer,
				    int iNumPasses, double k_extru,
				    double k_neigh)
{
  for (int iPass = 0; iPass < iNumPasses; iPass++) {
    for (int iVert = m_eVL.layerStartIndex(iEndLayer + 1) - 1;
	iVert >= m_eVL.layerStartIndex(iStartLayer); iVert--) {

      edamVert* pEV = m_eVL.get(iVert);
      if (pEV->isSteiner())
	continue;
      if (pEV->isNoSmooth())
	continue;

      double adForce[3] =
	{ 0, 0, 0 };

      // Move toward the average of neighboring vertex positions. Constrained to move in tangent plane.
      do {
	set<edamFace*> adjFaces = pEV->getAdjFaces();
	int iLayer = pEV->getLayer();

	int iNumFaces = adjFaces.size();
	if (iNumFaces < 1)
	  continue;

	double adUnitNorm[3];
	double adAvgDir[3];
	double adAvgDist[3];
	double dMinDist;
	neighborhoodInfo(pEV, adjFaces, adUnitNorm, adAvgDir, adAvgDist,
			 dMinDist, true);

	//double dExtFactor = (1 + dDOT3D(adAvgDir, adUnitNorm));
	double dExtFactor = 1;
	double dExt = m_firstLayerHeight * pow(m_growthRatio, iLayer)
	    * dExtFactor;

	// Get component of adAvgDist in the tangent plane.
	double dDotAvgDist = dDOT3D(adUnitNorm, adAvgDist);
	for (int ii = 0; ii < 3; ii++) {
	  // Reduce this force at high aspect ratios, to avoid highly skew elements in boundary layer.
	  adForce[ii] += k_neigh * fmin(3 * dExt / dMinDist, 1.0)
	      * (adAvgDist[ii] - adUnitNorm[ii] * dDotAvgDist);
	}
      }
      while (0);

      // Move toward the position at which parent would currently extrude.
      do {
	edamVert* pEVp = pEV->getParent();
	if (!pEVp) {
	  cout << "Vertex has no parent: " << endl;
	  cout << pEV->x() << ":" << pEV->y() << ":" << pEV->z() << endl;
	  assert(pEVp);
	}

	set<edamFace*> adjFaces = pEVp->getAdjFaces();
	int iLayer = pEVp->getLayer();

	int iNumFaces = adjFaces.size();
	if (iNumFaces < 1)
	  continue;

	double adUnitNorm[3];
	double adAvgDir[3];
	double adAvgDist[3];
	double dMinDist;
	neighborhoodInfo(pEVp, adjFaces, adUnitNorm, adAvgDir, adAvgDist,
			 dMinDist, true);

	//double dExtFactor = (1 + dDOT3D(adAvgDir, adUnitNorm));
	double dExtFactor = 1;
	double dExt = m_firstLayerHeight * pow(m_growthRatio, iLayer)
	    * dExtFactor;

	double adParDist[3] =
	  { pEV->x() - pEVp->x(), pEV->y() - pEVp->y(), pEV->z() - pEVp->z() };
	double dDotParDist = dDOT3D(adUnitNorm, adParDist);

	// If the parent vertex appears to be on a convex corner, increase this force to prevent inverted elements.
	double k_extru_curr = (
	    1 + dDOT3D(adAvgDir, adUnitNorm) < 0.8 ? 1 : k_extru);
	for (int ii = 0; ii < 3; ii++) {
	  adForce[ii] += k_extru_curr
	      * (adUnitNorm[ii] * dDotParDist - adParDist[ii]);
	  adForce[ii] += k_extru_curr * (adUnitNorm[ii] * (dExt - dDotParDist));
	}
      }
      while (0);

      //cout << adForce[0] <<":"<<adForce[1] <<":"<<adForce[2]<<endl;

      // Determine a new position for the vertex based on the net force.
      double adNewPos[3] =
	{ pEV->x() + adForce[0], pEV->y() + adForce[1], pEV->z() + adForce[2] };
      pEV->getVert()->setCoords(3, adNewPos);
    }
  }
}

/*
 * Compute "forces" acting on a vertex based on neighborhood vertices.
 * These are adapted from the forces in terrible2D.
 */
void
edamMesh3D_Closed::forceModelSmooth2(int iStartLayer, int iEndLayer,
				     int iNumPasses)
{

  double k_parent = 0.2;
  double k_kid = 0.2;
  double k_neigh = 0.5;
  double k_org = 1;

  for (int iPass = 0; iPass < iNumPasses; iPass++) {
    for (int iVert = m_eVL.layerStartIndex(iEndLayer + 1) - 1;
	iVert >= m_eVL.layerStartIndex(iStartLayer); iVert--) {

      edamVert* pEV = m_eVL.get(iVert);
      if (pEV->isSteiner())
	continue;
      if (pEV->isNoSmooth())
	continue;
      double adForce[3] =
	{ 0, 0, 0 };

      // Some of the initial bdry may not be marched from; these points
      // have no neighborhood.  Should only occur when marching from
      // level 0 to level 1.
      if (pEV->getNumFaces() == 0) {
	assert(pEV->getLayer() == 0);
	continue;
      }
      set<edamFace*> adjFaces = pEV->getAdjFaces();
      std::set<edamVert*> adjVerts;
      for (edamFace* pEF : adjFaces) {
	adjVerts.insert(pEF->getNextVert(pEV));
      }

      double adNorm[3];
      double adAvgDir[3];
      double adAvgDist[3];
      double dMinDist;
      neighborhoodInfo(pEV, adjFaces, adNorm, adAvgDir, adAvgDist, dMinDist);

      double dExtFactor = (1 + dDOT3D(adAvgDir, adNorm));
      dExtFactor *= dExtFactor;
      //double dExtFactor = 1;
      double dExtLength = m_firstLayerHeight
	  * pow(m_growthRatio, pEV->getLayer()) * dExtFactor;

      // Force node towards its original position.
      double dAspectRatio = pEV->getAspectRatio();
      do {
	double dAspectFactor = (
	    dAspectRatio > 1 ? (dAspectRatio < 6 ? 0.5 * (1 - cos(
	    M_PI * (dAspectRatio - 1) / 5)) :
						   1) :
			       0);
	adForce[0] += k_org * dAspectFactor * (pEV->xOrig() - pEV->x())
	    / dExtLength;
	adForce[1] += k_org * dAspectFactor * (pEV->yOrig() - pEV->y())
	    / dExtLength;
	adForce[2] += k_org * dAspectFactor * (pEV->zOrig() - pEV->z())
	    / dExtLength;
      }
      while (0);

      // Move toward the average of neighboring vertex positions.
      // Direction of movement is in tangent plane.
      do {
	double adAvgPos[3] =
	  { 0, 0, 0 };
	double dAvgDist = 0;
	for (edamFace* pEF : adjFaces) {
	  edamVert* pEV0 = pEF->getNextVert(pEV);
	  if (pEV0->isSteiner())
	    cout << "adj steiner" << endl;
	  adAvgPos[0] += pEV0->x() - pEV->x();
	  adAvgPos[1] += pEV0->y() - pEV->y();
	  adAvgPos[2] += pEV0->z() - pEV->z();
	  dAvgDist += pEV->getDist(pEV0);
	}
	for (int ii = 0; ii < 3; ii++) {
	  adAvgPos[ii] /= adjFaces.size();
	}
	dAvgDist /= adjFaces.size();

	// Get component of adAvgDist in tangent plane.
	double dDotAvgPos = dDOT3D(adNorm, adAvgPos);
	double adForceDir[3];
	for (int ii = 0; ii < 3; ii++) {
	  adForceDir[ii] = adAvgPos[ii] - adNorm[ii] * dDotAvgPos;
	}

	double dAspectFactor = (
	    dAspectRatio > 1 ? (dAspectRatio < 4 ? 0.8 + 0.1 * (1 - cos(
	    M_PI * (dAspectRatio - 1) / 3)) :
						   1) :
			       0.8);
	double dDistToTargetPos = dMAG3D(adForceDir);
	double dForceMagTangent = k_neigh * dAspectFactor * dDistToTargetPos
	    / dAvgDist;
	if (dDistToTargetPos > 1.e-10) {
	  NORMALIZE3D(adForceDir);
	  for (int ii = 0; ii < 3; ii++) {
	    adForce[ii] += dForceMagTangent * adForceDir[ii];
	  }
	}
	double dDistToAvgPos = dMAG3D(adAvgPos);
	double dForceMagFree = k_neigh * (1.0 - dAspectFactor) * dDistToAvgPos
	    / dAvgDist;
	if (dDistToAvgPos > 1.e-10) {
	  NORMALIZE3D(adAvgPos);
	  for (int ii = 0; ii < 3; ii++) {
	    adForce[ii] += dForceMagFree * adAvgPos[ii];
	  }
	}
      }
      while (0);

      // Force edge length to kid node closer to ideal.
      do {
	double adForceKid[3] =
	  { 0, 0, 0 };
	int iNumKids = 0;
	for (edamFace* pEFkid : pEV->getFaceList()) {

	  edamVert* pEVkid0 = pEFkid->getNextVert(pEV);
	  edamVert* pEVkid1 = pEFkid->getPrevVert(pEV);
	  if (pEVkid0->getLayer() < pEV->getLayer()
	      || pEVkid1->getLayer() < pEV->getLayer())
	    continue;
	  if (pEVkid0->isSteiner() || pEVkid1->isSteiner())
	    continue;

	  if (adjVerts.find(pEVkid0) == adjVerts.end()) {
	    iNumKids++;
	    double dLengthIdeal = dExtLength;
	    double adDistKid[3] =
	      { pEVkid0->x() - pEV->x(), pEVkid0->y() - pEV->y(), pEVkid0->z()
		  - pEV->z() };
	    //double dLengthReal = dMAG3D(adDistKid);
	    double dLengthReal = dDOT3D(adDistKid, adNorm);
	    double dForceMag = k_kid * (dLengthReal - dLengthIdeal)
		/ (dLengthReal + dLengthIdeal);
	    if (dLengthReal < 0.01 * dLengthIdeal)
	      dForceMag = 0;

	    if (pEVkid0->getLayer() == pEV->getLayer()) {
	      NORMALIZE3D(adDistKid);
	      adForceKid[0] += dForceMag * adDistKid[0];
	      adForceKid[1] += dForceMag * adDistKid[1];
	      adForceKid[2] += dForceMag * adDistKid[2];
	    }
	    else {
	      adForceKid[0] += dForceMag * adNorm[0];
	      adForceKid[1] += dForceMag * adNorm[1];
	      adForceKid[2] += dForceMag * adNorm[2];
	    }
	  }
	  if (adjVerts.find(pEVkid1) == adjVerts.end()) {
	    iNumKids++;
	    double dLengthIdeal = dExtLength;
	    double adDistKid[3] =
	      { pEVkid1->x() - pEV->x(), pEVkid1->y() - pEV->y(), pEVkid1->z()
		  - pEV->z() };
	    //double dLengthReal = dMAG3D(adDistKid);
	    double dLengthReal = dDOT3D(adDistKid, adNorm);
	    double dForceMag = k_kid * (dLengthReal - dLengthIdeal)
		/ (dLengthReal + dLengthIdeal);
	    if (dLengthReal < 0.01 * dLengthIdeal)
	      dForceMag = 0;

	    if (pEVkid1->getLayer() == pEV->getLayer()) {
	      NORMALIZE3D(adDistKid);
	      adForceKid[0] += dForceMag * adDistKid[0];
	      adForceKid[1] += dForceMag * adDistKid[1];
	      adForceKid[2] += dForceMag * adDistKid[2];
	    }
	    else {
	      adForceKid[0] += dForceMag * adNorm[0];
	      adForceKid[1] += dForceMag * adNorm[1];
	      adForceKid[2] += dForceMag * adNorm[2];
	    }
	  }
	}
	if (iNumKids > 0) {
	  for (int ii = 0; ii < 3; ii++) {
	    adForceKid[ii] /= iNumKids;
	    adForce[ii] += adForceKid[ii];
	  }
	}
      }
      while (0);

      // Force edge length to parent node closer to ideal.
      do {
	edamVert* pEVparent = pEV->getParent();
	if (pEVparent->isSteiner())
	  cout << "parent steiner" << endl;

	// Try to compute the normal as an angle-weighted average of adjacent face normals.
	vector<double> adParNorm(
	  { 0, 0, 0 });
	set<edamFace*> adjParFaces = pEVparent->getAdjFaces();
	for (edamFace* pEF : adjParFaces) {
	  double adFaceNorm[3];
	  pEF->calcUnitNormal(adFaceNorm);
	  double dAngle = pEF->calcAngle(pEVparent);
	  for (int ii = 0; ii < 3; ii++) {
	    adParNorm[ii] -= adFaceNorm[ii] * dAngle;
	  }
	}
	NORMALIZE3D(adParNorm);

	double dLengthIdeal = dExtLength / m_growthRatio;
	double adDistPar[3] =
	  { pEVparent->x() - pEV->x(), pEVparent->y() - pEV->y(), pEVparent->z()
	      - pEV->z() };
	double dLengthReal = dMAG3D(adDistPar);
	//double dLengthReal = -dDOT3D(adDistPar,adParNorm);
	double dForceMag = k_parent * (dLengthReal - dLengthIdeal)
	    / (dLengthReal + dLengthIdeal);

	adForce[0] -= dForceMag * adParNorm[0];
	adForce[1] -= dForceMag * adParNorm[1];
	adForce[2] -= dForceMag * adParNorm[2];
      }
      while (0);

      //cout << adForce[0] <<":"<<adForce[1] <<":"<<adForce[2]<<endl;

      // Determine a new position for the vertex based on the net force.
      double dSmoothScale = dExtLength;
      adForce[0] *= dSmoothScale;
      adForce[1] *= dSmoothScale;
      adForce[2] *= dSmoothScale;
      double adNewPos[3] =
	{ pEV->x() + adForce[0], pEV->y() + adForce[1], pEV->z() + adForce[2] };
      pEV->getVert()->setCoords(3, adNewPos);
      assert(isValid());
    }
  }
}

/*
 * Stores a representative aspect ratio for each initial vertex, used to control smoothing in the boundary layers.
 */
void
edamMesh3D_Closed::initAspectRatio()
{
  for (size_t iVert = 0; iVert < m_eVL.size(); iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner())
      continue;

    set<edamFace*> adjFaces = pEV->getAdjFaces();
    double dMinDist = LARGE_DBL;

    for (edamFace* pEF : adjFaces) {
      edamVert* pEV0 = pEF->getNextVert(pEV);
      dMinDist = fmin(dMinDist, pEV->getDist(pEV0));
    }

    pEV->setAspectRatio(dMinDist / m_firstLayerHeight);
  }
}

/*
 * TODO: Get in-layer neighboring faces and store them.
 * Need to do this after each layer is complete to distinguish these faces from others during smoothing.
 */
void
edamMesh3D_Closed::setAdjFaces(const int iLayer)
{
  for (size_t iVert = m_eVL.layerStartIndex(iLayer); iVert < m_eVL.size();
      iVert++) {
    edamVert* pEV = m_eVL.get(iVert);
    if (pEV->isSteiner())
      continue;

    pEV->clearAdjFaces();
    for (edamFace* pEF : pEV->getFaceList()) {
      if (pEF->getNumVerts() == 3 && pEF->getVert(0)->getLayer() == iLayer
	  && pEF->getVert(1)->getLayer() == iLayer
	  && pEF->getVert(2)->getLayer() == iLayer && pEF->getNumCells() < 2) {
	pEV->addFace(pEF);
      }
    }
  }
}

void
edamMesh3D_Closed::createPrismUnderFrontFace(edamFace* pEF)
{
  // Get the verts for this tri.
  Vert* pVA = pEF->getVert(0)->getVert();
  Vert* pVB = pEF->getVert(1)->getVert();
  Vert* pVC = pEF->getVert(2)->getVert();

  // And their parents...
  Vert* pVApar = pEF->getVert(0)->getParent()->getVert();
  Vert* pVBpar = pEF->getVert(1)->getParent()->getVert();
  Vert* pVCpar = pEF->getVert(2)->getParent()->getVert();

  replaceCellsWithPrism(pVA, pVB, pVC, pVApar, pVBpar, pVCpar);
}

int
edamMesh3D_Closed::createPrisms(const int layer)
{
  int retVal = 0;
  set<edamFace*> newFront = createFrontFaceList(layer, m_eVL);

  for (edamFace* pEF : newFront) {
    createPrismUnderFrontFace(pEF);
  }
  return retVal;
}

/*
 * TODO: Get in-layer neighboring faces and store them.
 * Need to do this after each layer is complete to distinguish these faces from others during smoothing.
 */
void
edamMesh3D_Closed::tagComplexCorners()
{
  for (size_t iVert = 0; iVert < m_eVL.size(); iVert++) {
    edamVert* pEV = m_eVL.get(iVert);

    vector<edamFace*> adjFaces;
    edamFace* pEFstart = pEV->getFace(0);
    edamFace* pEF = pEFstart;
    do {
      adjFaces.push_back(pEF);

      edamVert* pEV0 = pEF->getNextVert(pEV);
      edamVert* pEV1 = pEF->getNextVert(pEV0);
      pEF = pEV->getOtherCommonFace(pEV1, pEV0);
    }
    while (pEF != pEFstart);

    bool qConvex = false;
    bool qConcave = false;
    for (size_t ii = 0; ii < adjFaces.size(); ii++) {
      edamFace* pEF0 = adjFaces[(ii - 1 + adjFaces.size()) % adjFaces.size()];
      edamFace* pEF1 = adjFaces[ii];

      double adNorm0[3];
      pEF0->calcUnitNormal(adNorm0);
      double adNorm1[3];
      pEF1->calcUnitNormal(adNorm1);
      double dAngle = acos(dDOT3D(adNorm0, adNorm1));

      const double* ad0 = pEF0->getNextVert(pEV)->getVert()->getCoords();
      const double* ad1 =
	  pEF1->getNextVert(pEF1->getNextVert(pEV))->getVert()->getCoords();
      double adDist01[3] = adDIFF3D(ad1,ad0);

      if (dDOT3D(adNorm0,adDist01) > 0) {
	if (dAngle > M_PI / 4.0)
	  qConvex = true;
      }
      else {
	if (dAngle > M_PI / 4.0)
	  qConcave = true;
      }
    }

    if (qConvex && qConcave) {
      cout << "complex corner" << endl;
      pEV->setNoDecimate(true);
      pEV->setNoSmooth(true);
    }
  }
}

void
edamMesh3D_Closed::advanceMesh()
{
  writeFrontToVTK(m_eVL, -1, "initial-front.vtk");

  for (int ii = 0; ii < m_targetLayers; ii++) {
    logMessage(MSG_MANAGER, "Advancing to layer %d\n", ii + 1);
    checkForCollisions(ii);
    extrudeLayer(ii);
//		if(ii == numLayers-1){
//			removeVertsVisual(pM3D, m_eVL, ii);
//			createGoodCells(pM3D, m_eVL, ii);
//			break;
//		}
//    removeTempVertsFlipped(ii); // Decimate vertices with inverted adjacent faces.
//    removeTempVertsForProximity(ii); // Decimate vertices too close to neighbor vertices.
//    fixIntersections(ii);
//    removeTempVertsFromClipping(ii);
//    //closeTunnels(pM3D,m_eVL,ii);
//    createGoodCells(ii);

    int numWithoutPrisms = createPrisms(ii);
    assert(numWithoutPrisms == 0);
    assert(getNumFlakesInUse() == 0);


//    fillHoles(ii);
//    setAdjFaces(ii + 1);
    writeFrontToVTK(m_eVL, ii, "BeforeSmooth.vtk");
//    //if(ii == numLayers-1) break;
//    forceModelSmooth2(max(ii - 2, 1), ii + 1, 5);
//    writeFrontToVTK(m_eVL, ii, "AfterSmooth.vtk");
    writeTempMesh();
  }
}

/*
 * Real main program.
 */
int
main(int iNArg, char *apcArgs[])
{

// File names, etc. Input/output file names must be set in command line, length scale files (etc.) are optional.
  int iPChar;
  char strOutFileName[FILE_NAME_LEN];
  char strInFileName[FILE_NAME_LEN];
  char strQualFileName[FILE_NAME_LEN];

// These must be set in command line.
  double extrudeLen = 0;
  double growthRatio = 0;

  int numLayers = 50;

  edamMesh3D_Closed::eMarchDir eMD = edamMesh3D_Closed::eOut;

  std::set<int> BCset;
  BCset.insert(1);

  while ((iPChar = getopt(iNArg, apcArgs, "B:d:g:i:o:n:x:")) != EOF) {
    switch (iPChar)
      {
      case 'B':
	{
	  // The associated arg is a comma-separated list of BC indices.  Parse this into BCset.
	  char *readLoc = optarg, *commaLoc = NULL;
	  int read = 0;
	  do {
	    assert(readLoc);
	    int temp = -1;
	    read = sscanf(readLoc, "%d", &temp);
	    assert(read == 1);
	    BCset.insert(temp);
	    commaLoc = strstr(readLoc, ",");
	    readLoc = commaLoc + 1;
	  }
	  while (commaLoc != NULL);
	  break;
	}
      case 'd':
	{
	  switch (optarg[0])
	    {
	    case 'i':
	    case 'I':
	      eMD = edamMesh3D_Closed::eIn;
	      break;
	    case 'o':
	    case 'O':
	      eMD = edamMesh3D_Closed::eOut;
	      break;
	    case 'b':
	    case 'B':
	      // TODO:  This actually has no effect: the code just marches out.
	      eMD = edamMesh3D_Closed::eBoth;
	      break;
	    default:
	      // This is an error, but let it slide for now.
	      break;
	    }
	  break;
	}
      case 'g':
	{
	  sscanf(optarg, "%lf", &growthRatio); // Off-wall growth ratio per layer of cells.
	  break;
	}
      case 'i':
	{
	  strncpy(strInFileName, optarg, FILE_NAME_LEN - 1); // File name for input .bdry file.
	  break;
	}
      case 'n':
	{
	  sscanf(optarg, "%d", &numLayers); // Number of layers to extrude. Controls far-field size for a mesh with no specified outer boundary.
	  break;
	}
      case 'o':
	{
	  strncpy(strOutFileName, optarg, FILE_NAME_LEN - 1); // File name for final output mesh.
	  break;
	}
      case 'x':
	{
	  sscanf(optarg, "%lf", &extrudeLen); // Off-wall spacing.
	  break;
	}
      default:
	{
	  break;
	}
      }
  }
  openMessageFile(strInFileName);

  edamMesh3D_Closed EM3D(extrudeLen, growthRatio, numLayers, eMD);
//  logMessage(MSG_MANAGER, "Done with mesh input\n");
//
//  EM3D.tagComplexCorners();
//  EM3D.setAdjFaces(0);
//  EM3D.initAspectRatio();
//  EM3D.advanceMesh();

// Proposed new program flow

// Read surface mesh (eventually generate and/or accept one passed in) and initialize volume mesh
  Bdry3D *pB3D = new Bdry3D(strInFileName);
  SurfMesh SM(pB3D);
  EM3D.initVolMeshFromSurfMesh(SM);
  EM3D.disallowBdryChanges();

// Generate iso volume mesh with very coarse length scale (something like 4x
// iso lfs, for instance; this should be super coarse and fast to generate)
  std::shared_ptr<GRUMMP::Length> sizing_field = GRUMMP::Length::Create(&EM3D,
									0.5, 1);
  GRUMMP::CircumcenterIPC IPC;
  GRUMMP::RefinementManager3D RM3D(&EM3D, &IPC, nullptr, eBall, sizing_field, OldPriorityCalculator::SIZE_ONLY);
  	  RM3D.refineMesh();

// Set up initial EDAM front based on the surface of this volume mesh,
// including setting up initial marching distance.
  EM3D.initFront(BCset);

  EM3D.setWriteTempMesh(true);
  EM3D.advanceMesh();

// After each iteration, in addition to front decimation,
//    Identify iso verts that should be removed because they're about to
//        be overrun
//    Identify places where the front is close to collision, and stop
//        marching there; perhaps before front healing?

// Terminate when there are no active faces remaining in the front.

//  logMessage(MSG_MANAGER, "Cell count %d, of which:\n", EM3D.getNumCells());
//  logMessage(MSG_MANAGER, " %10d tetrahedra\n", EM3D.getNumTetCells());
//  logMessage(MSG_MANAGER, " %10d pyramids\n", EM3D.getNumPyrCells());
//  logMessage(MSG_MANAGER, " %10d prisms\n", EM3D.getNumPrismCells());
//  logMessage(MSG_MANAGER, " %10d hexahedra\n", EM3D.getNumHexCells());

  writeVTKLegacy(EM3D, strOutFileName);
  writeNative(EM3D, strOutFileName);
  EM3D.evaluateQuality();
  snprintf(strQualFileName, FILE_NAME_LEN, "%s.qual", strOutFileName);
  EM3D.writeQualityToFile(strQualFileName);
  delete pB3D;
}

