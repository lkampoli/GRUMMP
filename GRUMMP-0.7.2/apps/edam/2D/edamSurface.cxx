#include "GR_ParseStl.h"
#include "GR_config.h"
#include "GR_ADT.h"
#include "GR_GRGeom2D.h"
#include "GR_GRPoint.h"
#include "GR_Length.h"
#include "GR_Mesh2D.h"
#include "GR_misc.h"
#include "GR_SmoothingManager.h"
#include "GR_SurfaceInsertion.h"
#include "GR_SwapManager.h"
#include "GR_SwapDecider.h"
#include "GR_TriMeshBuilder.h"
#include <iomanip>

#include "GR_SurfMeshBuilder.h"
#include "GeometryQueryTool.hpp"
#include "FacetQueryEngine.hpp"
#include "GR_GraphicsOut.h"

//Global variables.
double m_extrudeLength;
double growthRatio;
int numLayers;
int numLoopsVTK;

/*
 * Faster approximation to atan(x).
 *
 * Calls to atan(x) account for ~10% of the program's runtime (or they did, at some point), this cut it down to ~3%.
 * Ultimately a small optimization, leaving it out to simplify code.
 */
//double dAtan(double x){
//	if(fabs(x) < 1){
//		return x*(1.0584 - ((x > 0) - (x < 0))*0.273*x);
//	}
//	else if(x > 0){
//		x = 1/x;
//		return M_PI_2 - x*(1.0584 - 0.273*x);
//	}
//	else{
//		x = 1/x;
//		return - M_PI_2 - x*(1.0584 + 0.273*x);
//	}
//}
/*
 * Wrappers for existing GRUMMP mesh entities, containing some additional data for advancing front meshing.
 */

struct edamFace;

class edamVert {
private:
	Vert* pV;
	int iLayer;
	bool qDeleted;
	double dExtrudeLength;
	// for smoothing
	double dOrgX;
	double dOrgY;
	double dOrgZ;
	double dAspectRatio;
	// for front intersection
	bool qClip;
	int iLoop;

	// for internal boundaries
	edamVert* pTVcopy;
	std::vector<edamFace*> faceList;
	edamFace* pTFL;
	edamFace* pTFR;
	std::vector<edamVert*> kids;
	std::vector<edamVert*> parents;
public:
	edamVert(); //never defined
	edamVert(Vert* pVert, int layer = -1) : // this is the constructor, not default constructor, hence arguments to be given
			pV(pVert), iLayer(layer), qDeleted(false), dExtrudeLength(0), dAspectRatio(
					0), qClip(false), iLoop(0), pTVcopy(NULL), faceList(), pTFL(
			NULL), pTFR(NULL), kids(), parents() {
		if (pV) {
			dOrgX = pV->x();
			dOrgY = pV->y();
			dOrgZ = pV->z();
		}
	}

	double x() {
		return pV->x();
	}
	double y() {
		return pV->y();
	}

	double z() {
		return pV->z();
	}

	void setCoords(double* dNewPos) {
		pV->setCoords(3, dNewPos);
	}
	void setOrgCoords(double* dNewPos) {
		setCoords(dNewPos);
		dOrgX = dNewPos[0];
		dOrgY = dNewPos[1];
		dOrgZ = dNewPos[2];
	}

	edamFace*
	getFace(int iFace) {
		return faceList.at(iFace);
	}

	double
	calcCellDihed(Cell* pC);
	edamFace*
	getCommonFace(edamVert* other);

	edamVert*
	getKid(int iKid) {
		return kids.at(iKid);
	}
	edamVert*
	getFirstKid() {
		return kids.at(0);
	}
	edamVert*
	getLastKid() {
		return kids.at(kids.size() - 1);
	}
	int getNumKids() {
		return kids.size();
	}
	edamVert*
	getParent(int iPar) {
		return parents.at(iPar);
	}
	int getNumParents() {
		return parents.size();
	}
	int getNumFaces() {
		return faceList.size();
	}

	Vert* getVert() const {
		return pV;
	}

	int getLayer() const {
		return iLayer;
	}

	bool isDeleted() const {
		return qDeleted;
	}

	void setDeleted(bool deleted) {
		qDeleted = deleted;
	}

	double getExtrudeLength() const {
		return dExtrudeLength;
	}

	void setExtrudeLength(double extrudeLength_in) {
		dExtrudeLength = extrudeLength_in;
	}

	double getOrigX() const {
		return dOrgX;
	}

	double getOrigY() const {
		return dOrgY;
	}

	double getOrigZ() const {
		return dOrgZ;
	}

	double getAspectRatio() const {
		return dAspectRatio;
	}

	void setAspectRatio(double aspectRatio) {
		dAspectRatio = aspectRatio;
	}

	bool isClip() const {
		return qClip;
	}

	void setClip(bool clip) {
		qClip = clip;
	}

	int getLoop() const {
		return iLoop;
	}

	void setLoop(int loop) {
		iLoop = loop;
	}

	const edamVert* getTVcopy() const {
		return pTVcopy;
	}

	void setTVcopy(edamVert* tVcopy) {
		pTVcopy = tVcopy;
	}

	const std::vector<edamFace*> getFaceList() const {
		return faceList;
	}

	void setFaceList(const std::vector<edamFace*> fL) {
		this->faceList = fL;
	}

	void push_backFace(edamFace* pTF) {
		faceList.push_back(pTF);
	}

	void eraseFace(int iF) {
		faceList.erase(faceList.begin() + iF);
	}

	edamFace* getTfl() const {
		return pTFL;
	}

	void setTfl(edamFace* tfl) {
		pTFL = tfl;
	}

	edamFace* getTfr() const {
		return pTFR;
	}

	void setTfr(edamFace* tfr) {
		pTFR = tfr;
	}

	void push_backKids(edamVert* pTV) {
		kids.push_back(pTV);
	}

	void eraseKids(int iK) {
		kids.erase(kids.begin() + iK);
	}

	edamVert* Kidsat(int iK) {
		return kids.at(iK);
	}

	void insertKids(int iK, edamVert* pTVkL) {
		kids.insert(kids.begin() + iK, pTVkL);
	}

	int sizeKids() {
		return kids.size();
	}

	void push_backParents(edamVert* pTV) {
		parents.push_back(pTV);
	}

	void eraseParents(int iP) {
		parents.erase(parents.begin() + iP);
	}

	edamVert* vGetParent(int i) {
		return this->parents[i];
	}

	const std::vector<edamVert*>& getParents() const {
		return parents;
	}

};

struct edamVertList {
	std::vector<std::vector<edamVert*>> vertList;

	edamVertList() :
			vertList() {
	}
	void addLayer() {
		vertList.push_back(vector<edamVert*>());
	}
	void push_back(edamVert TV, int iLayer) {
		vertList.at(iLayer).push_back(new edamVert(TV));
	}
	void cullDeleted(int iLayer) {
		vector<edamVert*> newLayer;
		for (edamVert* pTV : getLayer(iLayer)) {
			if (!pTV->isDeleted()) {
				newLayer.push_back(pTV);
			} else {
				delete pTV;
			}
		}
		vertList.at(iLayer) = newLayer;
	}

	vector<edamVert*> getLayer(int iLayer) {
		return vertList.at(iLayer);
	}
	edamVert*
	get(int iLayer, int index) {
		return vertList.at(iLayer).at(index);
	}
	edamVert*
	getLast(int iLayer) {
		return vertList.at(iLayer).back();
	}

	int size(int iLayer) {
		return vertList.at(iLayer).size();
	}
};

struct edamFace {
	Face* pF;
	edamVert* pTV0;
	edamVert* pTV1;
	bool qBdry;
	bool qRightCell;

	edamFace(); //never defined
	edamFace(Face* pFace, edamVert* pTV0a, edamVert* pTV1a, bool qBoundary =
			false) :
			pF(pFace), pTV0(pTV0a), pTV1(pTV1a), qBdry(qBoundary), qRightCell(
					false) {
	}

	edamVert*
	getVert(int iVert) {
		if (iVert == 0)
			return pTV1;
		return pTV0;
	}
	edamVert*
	getOppositeVert(edamVert* pTV) {
		if (pTV == pTV0)
			return pTV1;
		return pTV0;
	}
	edamVert*
	getOppositeVert(Cell* pC) {
		for (int iF = 0; iF < pTV0->getNumFaces(); iF++) {
			edamFace* pTF = pTV0->getFace(iF);
			if (pTF != this && pTF->pF->hasCell(pC))
				return pTF->getOppositeVert(pTV0);
		}

		assert(0);
		return NULL;
	}

	int getNumCells() {
		return pF->getNumCells();
	}
	Cell*
	getLeftCell(edamVert* pTVorient = NULL) {
		assert(pTVorient == NULL || pTVorient == pTV0 || pTVorient == pTV1);
		if (pTVorient) {
			return (pTVorient->getVert() == pF->getVert(0) ?
					pF->getLeftCell() : pF->getRightCell());
		} else {
			return (pTV0->getVert() == pF->getVert(0) ?
					pF->getLeftCell() : pF->getRightCell());
		}
	}
	Cell*
	getRightCell(edamVert* pTVorient = NULL) {
		assert(pTVorient == NULL || pTVorient == pTV0 || pTVorient == pTV1);
		if (pTVorient) {
			return (pTVorient->getVert() == pF->getVert(0) ?
					pF->getRightCell() : pF->getLeftCell());
		} else {
			return (pTV0->getVert() == pF->getVert(0) ?
					pF->getRightCell() : pF->getLeftCell());
		}
	}

	void calcUnitNorm(double* adRetVal, edamVert* pTVorient = NULL) {
		assert(pTVorient == NULL || pTVorient == pTV0 || pTVorient == pTV1);
		if (pTV1 == pTVorient) {
			adRetVal[0] = pTV0->y() - pTV1->y();
			adRetVal[1] = pTV1->x() - pTV0->x();
		} else {
			adRetVal[0] = pTV1->y() - pTV0->y();
			adRetVal[1] = pTV0->x() - pTV1->x();
		}
		NORMALIZE2D(adRetVal);
	}
	void calcUnitDir(double* adRetVal, edamVert* pTVorient = NULL) {
		assert(pTVorient == NULL || pTVorient == pTV0 || pTVorient == pTV1);
		if (pTV1 == pTVorient) {
			adRetVal[0] = pTV0->x() - pTV1->x();
			adRetVal[1] = pTV0->y() - pTV1->y();
		} else {
			adRetVal[0] = pTV1->x() - pTV0->x();
			adRetVal[1] = pTV1->y() - pTV0->y();
		}
		NORMALIZE2D(adRetVal);
	}

	static bool qPointerCompFaceAddress(const edamFace* pTF0,
			const edamFace* pTF1) {
		return pTF0 < pTF1;
	}
	static bool qPointerCompLeftEndpoint(const edamFace* pTF0,
			const edamFace* pTF1) {
		return fmin(pTF0->pTV0->x(), pTF0->pTV1->x())
				< fmin(pTF1->pTV0->x(), pTF1->pTV1->x());
	}
};

double edamVert::calcCellDihed(Cell* pC) {
	edamFace* pTFLcell = NULL;
	edamFace* pTFRcell = NULL;
	for (edamFace* pTF : faceList) {
		if (pTF->getLeftCell(this) == pC)
			pTFRcell = pTF;
		if (pTF->getRightCell(this) == pC)
			pTFLcell = pTF;
	}
	if (!pTFRcell) {
		cout << "no pTFRcell" << endl;
		cout << x() << ":" << y() << endl;
	}
	if (!pTFLcell) {
		cout << "no pTFLcell" << endl;
		cout << x() << ":" << y() << endl;
	}
	if (!pC) {
		cout << "no pC" << endl;
		cout << x() << ":" << y() << endl;
	}

	const double * const adCoord0 = pV->getCoords();
	const double * const adCoord1 =
			pTFRcell->getOppositeVert(this)->getVert()->getCoords();
	const double * const adCoord2 =
			pTFLcell->getOppositeVert(this)->getVert()->getCoords();

	double adNorm0[] =
			{ +adCoord1[1] - adCoord0[1], -adCoord1[0] + adCoord0[0] };
	double adNorm2[] =
			{ +adCoord0[1] - adCoord2[1], -adCoord0[0] + adCoord2[0] };

	NORMALIZE2D(adNorm0);
	NORMALIZE2D(adNorm2);

	double dArg0 = -dDOT2D(adNorm0, adNorm2);
	return GR_acos(dArg0);
}
edamFace*
edamVert::getCommonFace(edamVert* other) {
	for (edamFace* pTF : faceList) {
		if (pTF->getOppositeVert(this) == other) {
			return pTF;
		}
	}
	return NULL;
}

/*
 * Wrappers for GRUMMP mesh modifiers, use these to keep edamVertList consistent with the mesh.
 */

edamVert*
createVert(Mesh2D* pM2D, edamVertList &TVL, double x, double y, int layer) {
	Vert* pV = pM2D->createVert(x, y);
	TVL.push_back(edamVert(pV, layer), layer); // an example of calling the constructor of class edamVert
	return TVL.getLast(layer);
}

edamVert*
createVert(Mesh2D* pM2D, edamVertList &TVL, edamVert* pTVcopied) {
	assert(pM2D);
	TVL.push_back(edamVert(pTVcopied->getVert(), pTVcopied->getLayer()),
			pTVcopied->getLayer());
	return TVL.getLast(pTVcopied->getLayer());
}

/*
 * Create a face if it doesn't exist already. Otherwise, return the existing common face between the two vertices.
 */
edamFace*
createFace(Mesh2D* pM2D, edamVert* pTV0, edamVert* pTV1, bool qBdry = false) {
	bool qExist;

	Face* pF = pM2D->createFace(qExist, pTV0->getVert(), pTV1->getVert());
	if (!qExist) {
		edamFace* pTF = new edamFace(pF, pTV0, pTV1, qBdry);

		pTV0->push_backFace(pTF);
		pTV1->push_backFace(pTF);
		//pTV0->faceList.push_back(pTF);
		//pTV1->faceList.push_back(pTF);
		if (pTV0->getLayer() > pTV1->getLayer()) {
			pTV0->push_backParents(pTV1);
			pTV1->push_backKids(pTV0);
		} else if (pTV1->getLayer() > pTV0->getLayer()) {
			pTV1->push_backParents(pTV0);
			pTV0->push_backKids(pTV1);
		}
		return pTF;
	}

	return pTV0->getCommonFace(pTV1);
}
/*
 * Create a new face, allowing for duplicates. Duplicate faces are used to represent internal boundaries.
 */
edamFace*
createFaceDuplicate(Mesh2D* pM2D, edamVert* pTV0, edamVert* pTV1, bool qBdry =
		false) {
	bool qExist;

	Face* pF = pM2D->createFace(qExist, pTV0->getVert(), pTV1->getVert());

	edamFace* pTF = new edamFace(pF, pTV0, pTV1, qBdry);

	pTV0->push_backFace(pTF);
	pTV1->push_backFace(pTF);
	if (pTV0->getLayer() > pTV1->getLayer()) {
		pTV0->push_backParents(pTV1);
		pTV1->push_backKids(pTV0);
	} else if (pTV1->getLayer() > pTV0->getLayer()) {
		pTV1->push_backParents(pTV0);
		pTV0->push_backKids(pTV1);
	}
	return pTF;
}

void deleteFace(Mesh2D* pM2D, edamFace* pTF) {
	edamVert* pTV0 = pTF->pTV0;
	edamVert* pTV1 = pTF->pTV1;

	for (int iF = 0; iF < pTV0->getNumFaces(); iF++) {
		if (pTV0->getFace(iF) == pTF) {
			pTV0->eraseFace(iF);
			//pTV0->faceList.erase(pTV0->faceList.begin() + iF);
			break;
		}
	}
	for (int iK = 0; iK < pTV0->getNumKids(); iK++) {
		if (pTV0->getKid(iK) == pTV1) {
			pTV0->eraseKids(iK);
			break;
		}
	}
	for (int iP = 0; iP < pTV0->getNumParents(); iP++) {
		if (pTV0->getParent(iP) == pTV1) {
			pTV0->eraseParents(iP);
			break;
		}
	}
	if (pTV0->getTfl() == pTF) {
		pTV0->setTfl(NULL);
	}
	if (pTV0->getTfr() == pTF) {
		pTV0->setTfr(NULL);
	}
	for (int iF = 0; iF < pTV1->getNumFaces(); iF++) {
		if (pTV1->getFace(iF) == pTF) {
			pTV1->eraseFace(iF);
			//pTV1->faceList.erase(pTV1->faceList.begin() + iF);
			break;
		}
	}
	for (int iK = 0; iK < pTV1->getNumKids(); iK++) {
		if (pTV1->getKid(iK) == pTV0) {
			pTV1->eraseKids(iK);
			break;
		}
	}
	for (int iP = 0; iP < pTV1->getNumParents(); iP++) {
		if (pTV1->getParent(iP) == pTV0) {
			pTV1->eraseParents(iP);
			break;
		}
	}
	if (pTV1->getTfl() == pTF) {
		pTV1->setTfl(NULL);
	}
	if (pTV1->getTfr() == pTF) {
		pTV1->setTfr(NULL);
	}

	Face* pF = pTF->pF;
	if (pF) {
		pM2D->deleteFace(pF);
	}
	pTF->pF = NULL;
}
void deleteVert(Mesh2D* pM2D, edamVert* pTV) {
	for (edamFace* pTF : pTV->getFaceList()) {
		assert(0);
		deleteFace(pM2D, pTF);
	}
	pM2D->deleteVert(pTV->getVert());
	pTV->setDeleted(true);
}

Cell*
createTriCell(Mesh2D* pM2D, edamVert* TV0, edamVert* TV1, edamVert* TV2) {
	createFace(pM2D, TV0, TV1);
	createFace(pM2D, TV1, TV2);
	createFace(pM2D, TV2, TV0);
	return pM2D->createTriCell(TV0->getVert(), TV1->getVert(), TV2->getVert());
}
Cell*
createQuadCell(Mesh2D* pM2D, edamVert* TV0, edamVert* TV1, edamVert* TV2,
		edamVert* TV3) {
	createFace(pM2D, TV0, TV1);
	createFace(pM2D, TV1, TV2);
	createFace(pM2D, TV2, TV3);
	createFace(pM2D, TV3, TV0);
	return pM2D->createQuadCell(TV0->getVert(), TV1->getVert(), TV2->getVert(),
			TV3->getVert());
}
Cell*
createQuadOrTriCell(Mesh2D* pM2D, edamVert* TV0, edamVert* TV1, edamVert* TV2,
		edamVert* TV3) {

	Cell *pC;
	if (TV0 == TV1) {
		pC = createTriCell(pM2D, TV0, TV2, TV3);
	} else if (TV1 == TV2) {
		pC = createTriCell(pM2D, TV0, TV1, TV3);
	} else if (TV2 == TV3) {
		pC = createTriCell(pM2D, TV0, TV1, TV2);
	} else if (TV3 == TV0) {
		pC = createTriCell(pM2D, TV1, TV2, TV3);
	} else {
		pC = createQuadCell(pM2D, TV0, TV1, TV2, TV3);
	}
	return (pC);
}
void deleteCell(Mesh2D* pM2D, Cell* pC) {
	pM2D->deleteCellWithoutDeletingFaces(pC);
}

/*
 * Reading/writing mesh files, etc.
 */

double dAcuteAngleBdryVert(Vert* pBV) {
	assert(pBV->getNumFaces() == 2);
	Vert* pV0 = pBV->getFace(0)->getOppositeVert(pBV);
	Vert* pV1 = pBV->getFace(1)->getOppositeVert(pBV);
	double dist0[2] = { pV0->x() - pBV->x(), pV0->y() - pBV->y() };
	double dist1[2] = { pV1->x() - pBV->x(), pV1->y() - pBV->y() };

	return acos(
			fmax(
					fmin(dDOT2D(dist0,dist1) / (dMAG2D(dist0) * dMAG2D(dist1)),
							1.0), -1.0));
}

/*
 * Automatically generates a geometry file representing the boundaries of the given mesh.
 * Smooth curves (angles close to PI) of 4 or more vertices are written as splines, otherwise polylines.
 *
 * Need this for mesh curvature, etc.
 */
void meshToBdryFile(const char* fileName, Mesh2D* pMin) {

	int iNumCells = pMin->getNumCells();
	for (int iCell = 0; iCell < iNumCells; iCell++) {
		Cell* pC = pMin->getCell(iCell);
		if (pC->getNumFaces() != 1) {
			pMin->deleteCell(pC);
		}
	}

	int iNumFaces = pMin->getNumFaces();
	for (int iFace = 0; iFace < iNumFaces; iFace++) {
		Face* pF = pMin->getFace(iFace);
		if (pF->getNumCells() != 1) {
			assert(pF->getNumCells() == 0);
			pMin->deleteFace(pF);
		}
	}

	int iNumVerts = pMin->getNumVerts();
	for (int iVert = 0; iVert < iNumVerts; iVert++) {
		Vert* pV = pMin->getVert(iVert);
		if (pV->getNumFaces() != 2) {
			assert(pV->getNumFaces() == 0);
			pMin->deleteVert(pV);
		}
	}

	iNumCells = pMin->getNumCells();
	for (int iCell = 0; iCell < iNumCells; iCell++) {
		Cell* pC = pMin->getCell(iCell);
		pMin->deleteCell(pC);
	}

	pMin->purgeCells();
	pMin->purgeBdryFaces();
	pMin->purgeFaces();
	pMin->purgeVerts();

	ofstream fBdry;
	fBdry.open(fileName);
	fBdry.precision(16);

	iNumVerts = pMin->getNumVerts();
	bool qCreated[iNumVerts];
	for (int iV = 0; iV < iNumVerts; iV++) {
		qCreated[iV] = false;
	}

	vector<vector<int>> vecVecIndex;

	for (int iV = 0; iV < iNumVerts; iV++) {
		Vert* pV = pMin->getVert(iV);
		if (dAcuteAngleBdryVert(pV) < M_PI * 3.0 / 4.0) {
			assert(!qCreated[iV]);
			vector<int> vecIndex;
			Vert* pVnext;
			do {
				for (int iF = 0; iF < pV->getNumFaces(); iF++) {
					Face* pF = pV->getFace(iF);
					if ((pF->getVert(0) == pV && pF->getRightCell())
							|| (pF->getVert(1) == pV && pF->getLeftCell())) {
						pVnext = pF->getOppositeVert(pV);
					}
				}
				qCreated[pMin->getVertIndex(pV)] = true;
				vecIndex.push_back(pMin->getVertIndex(pV));
				pV = pVnext;
			} while (dAcuteAngleBdryVert(pV) > M_PI * 3.0 / 4.0);
			vecIndex.push_back(pMin->getVertIndex(pV));
			vecVecIndex.push_back(vecIndex);
		}
	}

	for (int iV = 0; iV < iNumVerts; iV++) {
		Vert* pV = pMin->getVert(iV);
		if (!qCreated[iV]) {
			vector<int> vecIndex;
			Vert* pVnext;
			do {
				assert(dAcuteAngleBdryVert(pV) > M_PI * 3.0/4.0);
				assert(!qCreated[pMin->getVertIndex(pV)]);
				for (int iF = 0; iF < 2; iF++) {
					Face* pF = pV->getFace(iF);
					if ((pF->getVert(0) == pV && pF->getRightCell())
							|| (pF->getVert(1) == pV && pF->getLeftCell())) {
						pVnext = pF->getOppositeVert(pV);
					}
				}
				qCreated[pMin->getVertIndex(pV)] = true;
				vecIndex.push_back(pMin->getVertIndex(pV));
				pV = pVnext;
			} while (!qCreated[pMin->getVertIndex(pV)]);
			vecIndex.push_back(pMin->getVertIndex(pV));
			vecVecIndex.push_back(vecIndex);
		}
	}

	fBdry << iNumVerts << " " << vecVecIndex.size() << endl;
	for (int iV = 0; iV < iNumVerts; iV++) {
		Vert* pV = pMin->getVert(iV);
		fBdry << pV->x() << " " << pV->y() << endl;
	}
	for (vector<int> vecIndex : vecVecIndex) {
		if (vecIndex.size() > 3) {
			fBdry << "spline r 1 b 1 " << vecIndex.size() << " ";
		} else {
			assert(vecIndex.size() > 1);
			fBdry << "polyline r 1 b 1 " << vecIndex.size() << " ";
		}
		for (int iIndex : vecIndex) {
			fBdry << iIndex << " ";
		}
		fBdry << endl;
	}
}

/*
 * Run GRUMMP's "tri" to generate a triangular mesh of the geometry.
 */
void meshFileFromBdryFile(const char* bdryFile, double dScale, double dGrading,
		const char* lengthFile = NULL) {

	char command[FILE_NAME_LEN];
	if (lengthFile) {
		sprintf(command,
				"/usr/anslab/src/GRUMMP/bin/linux-gnu/tri -i %s -r %f -g %f -l %s",
				bdryFile, dScale, dGrading, lengthFile);
	} else {
		sprintf(command,
				"/usr/anslab/src/GRUMMP/bin/linux-gnu/tri -i %s -r %f -g %f",
				bdryFile, dScale, dGrading);
	}

	if (system(command) != 0)
		assert(0);
}

/*
 * Strip out the interior of an existing mesh so as to use it's boundary discretization for advancing front.
 * If qIgnoreNonWalls, get rid of far-field boundaries as well (advised).
 */
void initFromExistingMesh(Mesh2D*& pMin, edamVertList &TVL,
		bool qIgnoreNonWalls) {
	TVL.addLayer();

	// Insist that all boundary faces are oriented CCW around the boundary.
	int iNumBFaces = pMin->getNumBdryFaces();
	for (int iBFace = 0; iBFace < iNumBFaces; iBFace++) {
		BFace* pBF = pMin->getBFace(iBFace);
		if (pBF->isDeleted())
			continue;

		//hack
		Face* pF = pBF->getFace(0);
		Cell* pC = pF->getOppositeCell(pBF);
		assert(
				checkOrient2D(pC->getVert(0), pC->getVert(1), pC->getVert(2))
						== 1);
		if (!(pF->getLeftCell() == pBF)) {
			pF->interchangeCellsAndVerts();
		}
		assert(
				checkOrient2D(pF->getVert(0), pF->getVert(1),
						pC->getOppositeVert(pF)) == -1);
		assert(pBF->getFace(0)->getLeftCell() == pBF);
	}

	int iNumCells = pMin->getNumCells();
	for (int iCell = 0; iCell < iNumCells; iCell++) {
		Cell* pC = pMin->getCell(iCell);
		pMin->deleteCell(pC);
	}
	int iNumFaces = pMin->getNumFaces();
	for (int iFace = 0; iFace < iNumFaces; iFace++) {
		Face* pF = pMin->getFace(iFace);
		if (!pF->isBdryFace()) {
			pMin->deleteFace(pF);
		}
	}
	if (qIgnoreNonWalls) {
		for (int iBFace = 0; iBFace < iNumBFaces; iBFace++) {
			BFace* pBF = pMin->getBFace(iBFace);
			if (pBF->isDeleted())
				continue;

			assert(pBF->getFace(0)->getLeftCell() == pBF);

			if (pBF->getBdryCondition() == 5) {
				pMin->deleteBFace(pBF);
			}
		}
	}
	int iNumVerts = pMin->getNumVerts();
	for (int iVert = 0; iVert < iNumVerts; iVert++) {
		Vert* pV = pMin->getVert(iVert);
		if (pV->getNumFaces() == 0) {
			pMin->deleteVert(pV);
		}
	}
	pMin->purgeCells();
	pMin->purgeBdryFaces();
	pMin->purgeFaces();
	pMin->purgeVerts();

	Mesh2D* pMout = new Mesh2D(2);
	pMout->allowNonSimplicial();

	map<edamVert*, int> vertIndexMap;
	iNumVerts = pMin->getNumVerts();
	for (int iV = 0; iV < iNumVerts; iV++) {
		Vert* pV = pMin->getVert(iV);
		edamVert* pTV = createVert(pMout, TVL, pV->x(), pV->y(), 0);
		vertIndexMap.insert(pair<edamVert*, int>(pTV, iV));
	}

	iNumFaces = pMin->getNumFaces();
	for (int iF = 0; iF < iNumFaces; iF++) {
		Face* pF = pMin->getFace(iF);
		int iV0 = pMin->getVertIndex(pF->getVert(0));
		int iV1 = pMin->getVertIndex(pF->getVert(1));
		if (pF->getFaceLoc() == Face::eBdryTwoSide
				|| pF->getFaceLoc() == Face::ePseudoSurface) {
			continue;
		} else if (pF->getLeftCell()) {
			edamFace* pTF = createFace(pMout, TVL.get(0, iV0), TVL.get(0, iV1),
					true);
			TVL.get(0, iV0)->setTfl(pTF);
			TVL.get(0, iV1)->setTfr(pTF);
		} else {
			assert(0);
		}
	}

	int iNumDupes = 0;
	// Write two copies of verts on internal boundaries so that the rest of the code can treat them the same as other boundaries.
	for (int iV = 0; iV < iNumVerts; iV++) {
		edamVert* pTVstart = TVL.get(0, iV);
		if (pTVstart->getTVcopy())
			continue;

		Vert* pVstart = pMin->getVert(vertIndexMap.at(pTVstart));
		if (pTVstart->getNumFaces() == 0 && pVstart->getNumFaces() == 4) {
			edamVert* pTVduped = createVert(pMout, TVL, pTVstart);
			vertIndexMap.insert(
					pair<edamVert*, int>(pTVduped, iNumVerts + iNumDupes));
			iNumDupes++;

			set<edamVert*> vertSet;
			for (int iF = 0; iF < pVstart->getNumFaces(); iF++) {
				Face* pF = pVstart->getFace(iF);
				vertSet.insert(
						TVL.get(0,
								pMin->getVertIndex(
										pF->getOppositeVert(pVstart))));
			}
			assert(vertSet.size() == 2);
			vector<edamVert*> verts(vertSet.begin(), vertSet.end());

			verts[0]->setTfr(
					createFaceDuplicate(pMout, pTVstart, verts[0], true));
			pTVstart->setTfl(verts[0]->getTfr());
			verts[1]->setTfl(
					createFaceDuplicate(pMout, verts[1], pTVstart, true));
			pTVstart->setTfr(verts[1]->getTfl());
			verts[1]->setTfr(
					createFaceDuplicate(pMout, pTVduped, verts[1], true));
			pTVduped->setTfl(verts[1]->getTfr());
			verts[0]->setTfl(
					createFaceDuplicate(pMout, verts[0], pTVduped, true));
			pTVduped->setTfr(verts[0]->getTfl());

			for (int ii = 0; ii < 2; ii++) {
				edamVert* pTV = verts[ii];
				Vert* pV = pMin->getVert(vertIndexMap.at(pTV));
				while (pTV->getNumFaces() == 2 && pV->getNumFaces() == 4) {
					edamVert* pTVdupe = createVert(pMout, TVL, pTV);
					vertIndexMap.insert(
							pair<edamVert*, int>(pTVdupe,
									iNumVerts + iNumDupes));
					iNumDupes++;

					edamVert* pTVnext = NULL;

					edamVert* pTVL = pTV->getTfl()->getOppositeVert(pTV);
					edamVert* pTVR = pTV->getTfr()->getOppositeVert(pTV);
					for (int iF = 0; iF < pV->getNumFaces(); iF++) {
						Face* pF = pV->getFace(iF);
						Vert* pVopp = pF->getOppositeVert(pV);
						if (static_cast<int>(pMin->getVertIndex(pVopp))
								!= vertIndexMap.at(pTVL)
								&& static_cast<int>(pMin->getVertIndex(pVopp))
										!= vertIndexMap.at(pTVR)) {
							pTVnext = TVL.get(0, pMin->getVertIndex(pVopp));
						}
					}
					assert(pTVnext);

					edamVert* pTVprev = pTV->getTfr()->getOppositeVert(pTV);
					for (int iF = 0; iF < pTV->getNumFaces(); iF++) {
						if (pTV->getFace(iF) == pTV->getTfr()) {
							pTV->eraseFace(iF);
							break;
						}
					}
					pTVprev->setTfl(
							createFaceDuplicate(pMout, pTVprev, pTVdupe, true));
					pTVdupe->setTfr(pTVprev->getTfl());
					pTVnext->setTfr(
							createFaceDuplicate(pMout, pTVdupe, pTVnext, true));
					pTVdupe->setTfl(pTVnext->getTfr());
					deleteFace(pMout, pTV->getTfr());
					pTVnext->setTfl(
							createFaceDuplicate(pMout, pTVnext, pTV, true));
					pTV->setTfr(pTVnext->getTfl());
					pTV->setTVcopy(pTVdupe);
					pTVdupe->setTVcopy(pTV);

					pTV = pTVnext;
					pV = pMin->getVert(vertIndexMap.at(pTV));
				}
			}
		}
	}

	// Write two copies of verts on internal boundaries so that the rest of the code can treat them the same as other boundaries.
	for (int iVert = 0; iVert < iNumVerts; iVert++) {
		edamVert* pTV = TVL.get(0, iVert);
		if (pTV->getTVcopy())
			continue;

		Vert* pV = pMin->getVert(vertIndexMap.at(pTV));
		assert(pTV->getNumFaces() > 0);
		if (pTV->getNumFaces() > 2) {
			if (pTV->getNumFaces() != pV->getNumFaces()) {
				cout << "bad: " << pV->x() << ":" << pV->y() << " "
						<< pTV->getNumFaces() << ":" << pV->getNumFaces()
						<< endl;
			}
			assert(pTV->getNumFaces() == pV->getNumFaces());
			assert(pTV->getNumFaces() % 2 == 0);
			vector<pair<double, edamVert*>> vertsR;
			vector<pair<double, edamVert*>> vertsL;
			for (int iF = pTV->getNumFaces() - 1; iF >= 0; iF--) {
				edamVert* pTVopp = pTV->getFace(iF)->getOppositeVert(pTV);

				if (pTVopp->getTfl()->getOppositeVert(pTVopp) == pTV) {
					vertsR.push_back(pair<double, edamVert*>(0, pTVopp));
				} else {
					vertsL.push_back(pair<double, edamVert*>(0, pTVopp));
				}
				deleteFace(pMout, pTV->getFace(iF));
			}
			for (int iV = 0; iV < static_cast<int>(vertsR.size()); iV++) {
				edamVert* pTVR = vertsR.at(iV).second;
				double dAngle = atan2(pTVR->y() - pTV->y(),
						pTVR->x() - pTV->x());
				vertsR.at(iV).first = dAngle;
			}
			for (int iV = 0; iV < static_cast<int>(vertsL.size()); iV++) {
				edamVert* pTVL = vertsL.at(iV).second;
				double dAngle = atan2(pTVL->y() - pTV->y(),
						pTVL->x() - pTV->x());
				vertsL.at(iV).first = dAngle;
			}
			std::sort(vertsR.begin(), vertsR.end());
			std::sort(vertsL.begin(), vertsL.end());
			if (vertsL.at(0).first > vertsR.at(0).first + 1.e-10) {
				for (int ii = 0; ii < static_cast<int>(vertsL.size() - 1);
						ii++) {
					std::iter_swap(vertsL.begin() + ii,
							vertsL.begin() + ii + 1);
				}
			}

			for (int iV = 0; iV < static_cast<int>(vertsR.size()); iV++) {
				edamVert* pTVnew = pTV;
				if (iV != 0) {
					pTVnew = createVert(pMout, TVL, pTV);
					vertIndexMap.insert(
							pair<edamVert*, int>(pTVnew,
									iNumVerts + iNumDupes));
					iNumDupes++;
				}
				edamVert* pTVR = vertsR.at(iV).second;
				edamVert* pTVL = vertsL.at((iV + 1) % vertsL.size()).second;

				if (!pTVR->getCommonFace(pTVnew)) {
					pTVnew->setTfr(
							createFaceDuplicate(pMout, pTVR, pTVnew, true));
					pTVR->setTfl(pTVnew->getTfr());
				}
				if (!pTVL->getCommonFace(pTVnew)) {
					pTVL->setTfr(
							createFaceDuplicate(pMout, pTVnew, pTVL, true));
					pTVnew->setTfl(pTVL->getTfr());
				}
			}
		}
	}

	for (int iV = 0; iV < iNumVerts; iV++) {
		edamVert* pTV = TVL.get(0, iV);
		assert(pTV->getNumFaces() == 2);
	}

	assert(pMout->getNumVerts() == pMin->getNumVerts());

	delete pMin;

	pMin = pMout;
}

Mesh2D*
initMeshFromMeshFile(const char* fileName, edamVertList &TVL,
		bool qIgnoreNonWalls) {
	Mesh2D* pMin = new Mesh2D(2);
	pMin->readFromFileTerrible(fileName);

	initFromExistingMesh(pMin, TVL, qIgnoreNonWalls);

	return pMin;
}

/*
 * Returns true if the two line segments (1-2 , 3-4) intersect.
 *
 * Currently does NOT handle degenerate cases, assumes parallel lines do not overlap.
 */
bool qIntersect(double x1, double y1, double x2, double y2, double x3,
		double y3, double x4, double y4, double* adPos = NULL) {
	double d12[2] = { x2 - x1, y2 - y1 };
	double d34[2] = { x4 - x3, y4 - y3 };

	// if parallel
	if (fabs(dDOT2D(d12, d34)) > dMAG2D(d12) * dMAG2D(d34) - 1e-10) {
		return false;
	}

	double t12 = ((x1 - x3) * d34[1] - (y1 - y3) * d34[0])
			/ (d34[0] * d12[1] - d34[1] * d12[0]);
	double t34 = ((y3 - y1) * d12[0] - (x3 - x1) * d12[1])
			/ (d34[0] * d12[1] - d34[1] * d12[0]);

	// store the intersection point location
	if (adPos) {
		adPos[0] = (1 - t12) * x1 + t12 * x2;
		adPos[1] = (1 - t12) * y1 + t12 * y2;
	}
	return t12 > 1e-10 && t12 < 1 - 1e-10 && t34 > 1e-10 && t34 < 1 - 1e-10;
}
bool qIntersect(edamFace* pTF0, edamFace* pTF1, double* adPos =
NULL) {
	double x1 = pTF0->pTV0->x();
	double y1 = pTF0->pTV0->y();
	double x2 = pTF0->pTV1->x();
	double y2 = pTF0->pTV1->y();
	double x3 = pTF1->pTV0->x();
	double y3 = pTF1->pTV0->y();
	double x4 = pTF1->pTV1->x();
	double y4 = pTF1->pTV1->y();
	return qIntersect(x1, y1, x2, y2, x3, y3, x4, y4, adPos);
}

/*
 *  Find pairs of front faces that intersect, and split ("clip") them at points of intersection.
 *  Two identical vertices are added at the intersection points to divide the front into separate, closed, non-intersecting loops (we delete duplicate vertices later).
 */
void clipIntersections(Mesh2D* pM2D, edamVertList &TVL, int iLayer) {

	// Find the faces on the new front.
	std::vector<edamFace*> segmentList;
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		if (!pTV->getTfl()) {
			cout << pTV->x() << ":" << pTV->y() << endl;
			assert(pTV->getTfl());
		}
		segmentList.push_back(pTV->getTfl());
	}
	if (segmentList.empty())
		return;

	// Construct a bounding box tree for efficient intersection checks.
	double (*a2dBB)[4] = new double[segmentList.size()][4];
	for (int iSeg = 0; iSeg < static_cast<int>(segmentList.size()); iSeg++) {
		edamFace* pTF = segmentList.at(iSeg);
		a2dBB[iSeg][0] = fmin(pTF->pTV0->x(), pTF->pTV1->x());
		a2dBB[iSeg][1] = fmax(pTF->pTV0->x(), pTF->pTV1->x());
		a2dBB[iSeg][2] = fmin(pTF->pTV0->y(), pTF->pTV1->y());
		a2dBB[iSeg][3] = fmax(pTF->pTV0->y(), pTF->pTV1->y());
	}
	ADT* ADTtree = new ADT(segmentList.size(), 2, ADT::eBBoxes,
			reinterpret_cast<double*>(a2dBB));

	// Iterate over all faces to detect intersections.
	for (int ii = 0; ii < static_cast<int>(segmentList.size()); ii++) {
		edamFace* pTF = segmentList.at(ii);
		if (!pTF->pF)
			continue; // face has been deleted already

		edamVert* pTV0 = pTF->pTV0;
		edamVert* pTV1 = pTF->pTV1;
		double adBB[4];
		adBB[0] = fmin(pTF->pTV0->x(), pTF->pTV1->x());
		adBB[1] = fmax(pTF->pTV0->x(), pTF->pTV1->x());
		adBB[2] = fmin(pTF->pTV0->y(), pTF->pTV1->y());
		adBB[3] = fmax(pTF->pTV0->y(), pTF->pTV1->y());

		// Query the bounding box tree for nearby faces.
		std::vector<GR_index_t> result = ADTtree->veciRangeQuery(adBB);

		// Check if each nearby face actually intersects this one.
		for (int jj : result) {
			if (jj > ii) { // Don't check the same pair of faces twice.
				edamFace* pTFo = segmentList.at(jj);
				if (!pTFo->pF)
					continue;

				edamVert* pTV2 = pTFo->pTV0;
				edamVert* pTV3 = pTFo->pTV1;

				double adPos[2];
				if (qIntersect(pTV0->x(), pTV0->y(), pTV1->x(), pTV1->y(),
						pTV2->x(), pTV2->y(), pTV3->x(), pTV3->y(), adPos)) {
					// Intersection found!!!
					edamFace* newFaces[4];

					// Create two identical vertices at the point of intersection, delete the two faces and create four new ones.
					edamVert* pTVnew03 = createVert(pM2D, TVL, adPos[0],
							adPos[1], iLayer + 1);
					pTVnew03->setTfr(createFace(pM2D, pTV0, pTVnew03));
					pTV0->setTfl(pTVnew03->getTfr());
					newFaces[0] = pTV0->getTfl();
					pTVnew03->setTfl(createFace(pM2D, pTVnew03, pTV3));
					pTV3->setTfr(pTVnew03->getTfl());
					newFaces[3] = pTV3->getTfr();
					pTVnew03->setClip(true);
					edamVert* pTVnew21 = createVert(pM2D, TVL, adPos[0],
							adPos[1], iLayer + 1);
					pTVnew21->setTfr(createFace(pM2D, pTV2, pTVnew21));
					pTV2->setTfl(pTVnew21->getTfr());
					newFaces[2] = pTV2->getTfl();
					pTVnew21->setTfl(createFace(pM2D, pTVnew21, pTV1));
					pTV1->setTfr(pTVnew21->getTfl());
					newFaces[1] = pTV1->getTfr();
					pTVnew21->setClip(true);

					// Add new faces to the bounding box tree, as they might still intersect others.
					for (edamFace* pTFnew : newFaces) {
						double adBBnew[4];
						adBBnew[0] = fmin(pTFnew->pTV0->x(), pTFnew->pTV1->x());
						adBBnew[1] = fmax(pTFnew->pTV0->x(), pTFnew->pTV1->x());
						adBBnew[2] = fmin(pTFnew->pTV0->y(), pTFnew->pTV1->y());
						adBBnew[3] = fmax(pTFnew->pTV0->y(), pTFnew->pTV1->y());

						ADTtree->vAddNode(adBBnew, ADTtree->iTreeSize());
						segmentList.push_back(pTFnew);
					}

					deleteFace(pM2D, pTF);
					deleteFace(pM2D, pTFo);
					break;
				}
			}
		}
	}
}

/*
 * Write a list of faces to a VTK file for visualization.
 */
void VTKFaceList(vector<edamFace*> &faces, const char* filename) {
	set<edamVert*> verts;
	for (edamFace* pTF : faces) {
		for (int iV = 0; iV < 2; iV++) {
			verts.insert(pTF->getVert(iV));
		}
	}
	vector<edamVert*> vertList(verts.begin(), verts.end());
	map<edamVert*, int> vertIndexMap;

	ofstream fVTK;
	fVTK.open(filename);
	fVTK.precision(16);

	fVTK << "# vtk DataFile Version 1.0" << endl;
	fVTK << "OVMesh Tri example" << endl;
	fVTK << "ASCII" << endl;
	fVTK << endl;
	fVTK << "DATASET UNSTRUCTURED_GRID" << endl;

	fVTK << "POINTS " << verts.size() << " float" << endl;
	for (int iV = 0; iV < static_cast<int>(vertList.size()); iV++) {
		edamVert* pTV = vertList.at(iV);
		fVTK << pTV->x() << " " << pTV->y() << " 0" << endl;

		vertIndexMap.insert(pair<edamVert*, int>(pTV, iV));
	}
	fVTK << endl;

	fVTK << "CELLS " << faces.size() << " " << faces.size() * 3 << endl;
	for (edamFace* pTF : faces) {
		fVTK << 2;
		for (int iV = 0; iV < 2; iV++) {
			fVTK << " " << vertIndexMap.at(pTF->getVert(iV));
		}
		fVTK << endl;
	}
	fVTK << endl;

	fVTK << "CELL_TYPES " << faces.size() << endl;
	for (int ii = 0; ii < static_cast<int>(faces.size()); ii++) {
		fVTK << "3" << endl;
	}
}

/*
 * Special case of VTKFaceList for visualizing the entire new front. Very useful for error diagnosis.
 */
void loopsVTK(edamVertList &TVL, int iLayer) {

	char filename[256];
	sprintf(filename, "nextFront%d.vtk", numLoopsVTK);
	numLoopsVTK++;

	vector<edamFace*> faces;
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		faces.push_back(pTV->getTfl());
	}

	VTKFaceList(faces, filename);
}

/*
 * Delete vertices who have no in-layer neighbors (an entire loop of parents extruded to a single point.
 * Shouldn't happen anymore, after reworking the edge collapse...
 */
void deleteIsolatedVerts(Mesh2D* pM2D, edamVertList &TVL, int iLayer) {
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		if (!pTV->getTfl()) {
			while (pTV->getNumFaces() != 0) {
				deleteFace(pM2D, pTV->getFace(0));
			}
			deleteVert(pM2D, pTV);
		}
	}
	TVL.cullDeleted(iLayer + 1);
}
/*
 * If a vertex has a parent edge outward between its two in-layer neighbors, that parent edge necessarily intersects some other new front edge.
 * "Decimate" the vertex with the illegal parent edge.
 */
void deleteBadVerts(Mesh2D* pM2D, edamVertList &TVL, int iLayer) {
	for (int iPass = 0; iPass < 5; iPass++) {
		for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
			if (pTV->isDeleted())
				continue;

			edamVert* pTVR = pTV->getTfr()->getOppositeVert(pTV);
			edamVert* pTVL = pTV->getTfl()->getOppositeVert(pTV);

			bool qBad = false;

			for (edamVert* pTVpar : pTV->getParents()) {
				//for(int i=0; i < size; i++){
				//  edamVert* pTVpar = pTV->vGetParent(i);

				assert(pTVpar);

				// Orientation checks determine if the parent edge is illegal.
				if (checkOrient2D(pTV->getVert(), pTVR->getVert(),
						pTVL->getVert()) == 1) {
					if (checkOrient2D(pTV->getVert(), pTVR->getVert(),
							pTVpar->getVert()) != -1
							&& checkOrient2D(pTVL->getVert(), pTV->getVert(),
									pTVpar->getVert()) != -1)
						qBad = true;
				} else {
					if (checkOrient2D(pTV->getVert(), pTVR->getVert(),
							pTVpar->getVert()) != -1
							|| checkOrient2D(pTVL->getVert(), pTV->getVert(),
									pTVpar->getVert()) != -1)
						qBad = true;
				}
			}

			if (qBad) {
				// This shouldn't happen anymore (and probably doesn't), but...
				// If this front only has two vertices, delete them both,
				if (pTVR == pTVL) {
					while (pTVR->getNumFaces() != 0) {
						deleteFace(pM2D, pTVR->getFace(0));
					}
					deleteVert(pM2D, pTVR);
				}
				// If trying to decimate one of the verts on a triangular (3-vertex) front, just delete all 3 vertices.
				else if (pTVR->getCommonFace(pTVL)) {
					while (pTVR->getNumFaces() != 0) {
						deleteFace(pM2D, pTVR->getFace(0));
					}
					deleteVert(pM2D, pTVR);

					while (pTVL->getNumFaces() != 0) {
						deleteFace(pM2D, pTVL->getFace(0));
					}
					deleteVert(pM2D, pTVL);
				}
				// In the usual case, just connect the right and left neighbors.
				else {
					pTVL->setTfr(createFace(pM2D, pTVR, pTVL));
					pTVR->setTfl(pTVL->getTfr());
				}

				// Finally, delete the offending vertex.
				while (pTV->getNumFaces() != 0) {
					deleteFace(pM2D, pTV->getFace(0));
				}
				deleteVert(pM2D, pTV);
			}
		}
	}
	TVL.cullDeleted(iLayer + 1);
}

// Compute winding number as sum of contributions from individual faces.
double calcWindingNumber(const set<edamFace*> loop,
		const set<edamFace*> &allFaces) {

	edamFace* pTFstart = *(loop.begin());
	double adOffset[2];

	/*
	 * For each loop, find a point immediately outside (a.k.a. offset in the marching direction from) the loop and compute the winding number
	 * of the whole front around that point.
	 * Every loop which is a valid new front gives a winding number of...
	 * 		0, if there is no outer boundary
	 * 		2*PI, otherwise
	 */
	do {
		double adCent[2];
		double adNorm[2];
		pTFstart->pF->calcCentroid(adCent);
		pTFstart->pF->calcUnitNormal(adNorm);
		for (int ii = 0; ii < 2; ii++) {
			adOffset[ii] = adCent[ii] - 1.e-7 * adNorm[ii];
		}
	} while (0);

	double dWind = 0;
	for (edamFace* pTF : allFaces) {
		const double* ad0 = pTF->getVert(0)->getVert()->getCoords();
		const double* ad1 = pTF->getVert(1)->getVert()->getCoords();
		double adA[2] = adDIFF2D(ad0,adOffset);
		double adB[2] = adDIFF2D(ad1,adOffset);

		double dNumer = adA[0] * adB[1] - adA[1] * adB[0];
		double dDenom = dDOT2D(adA, adB);
		double dSolidAngle = atan(dNumer / dDenom);

		double adNorm[2];
		pTF->pF->calcUnitNormal(adNorm);
		if (dSolidAngle > 0 && dDOT3D(adNorm,adA) < 0) {
			dSolidAngle -= M_PI;
		} else if (dSolidAngle < 0 && dDOT3D(adNorm,adA) > 0) {
			dSolidAngle += M_PI;
		}
		dWind += dSolidAngle;
	}

	return dWind;
}

/*
 * Delete the loops which have crossed over each other or are a result of a loop crossing over itself
 */
void deleteBadLoops(Mesh2D* pM2D, edamVertList &TVL,
		vector<set<edamFace*>> &loops, double &dWindInitial,
		const int &iLayer) {

	// Store the set of all faces (for the winding number computation).
	set<edamFace*> allFaces;
	for (set<edamFace*> loop : loops) {
		allFaces.insert(loop.begin(), loop.end());
	}

	// Set the initial winding number as the winding number of the first loop in the first layer of extruded loops
	// have to check this for robustness in complicated cases with several loops inside one another.
	if (iLayer == 0)
		dWindInitial = calcWindingNumber(*loops.begin(), allFaces);

	vector<bool> loopsBad;
	for (set<edamFace*> loop : loops) {
		if (loop.empty())
			continue;

		//Calculate the winding number of this loop, ie., a point placed towards the marching direction of this loop
		double dWind = calcWindingNumber(loop, allFaces);

		if (fabs(dWindInitial) < M_PI) {
			if (fabs(dWindInitial) >= M_PI) {// TODO: automatically figure out the right winding number (or at least take it as a user input...)
				//if ((dWind) < -M_PI)
				//if( dWind > -pow(10,-3) && dWind < pow(10, -3)){
				//if( dWind > 2*M_PI-pow(10,-3) && dWind < 2*M_PI+ pow(10,-3) ){
				loopsBad.push_back(true);
			} else {
				loopsBad.push_back(false);
			}
		} else if (fabs(dWindInitial - 2 * M_PI) < M_PI) {
			if (fabs(dWind - 2 * M_PI) > M_PI) {
				loopsBad.push_back(true);
			} else {
				loopsBad.push_back(false);
			}
		}
	}

	// Now that we've identified bad loops, delete them.
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		if (loopsBad.at(pTV->getLoop() - 1)) {
			while (pTV->getNumFaces() != 0) {
				deleteFace(pM2D, pTV->getFace(0));
			}
			deleteVert(pM2D, pTV);
		}
	}
}

// "Decimate" vertices that we've added at points of self-intersection.
void decimateExtraVerticesAtIntersections(Mesh2D* pM2D, edamVertList &TVL,
		const int &iLayer) {
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		if (pTV->isDeleted())
			continue;

		if (pTV->isClip()) {
			edamVert* pTVR = pTV->getTfr()->getOppositeVert(pTV);
			edamVert* pTVL = pTV->getTfl()->getOppositeVert(pTV);

			// Decimation, see comments in duplicate code from deleteBadVerts.
			if (pTVR == pTVL) {
				while (pTVR->getNumFaces() != 0) {
					deleteFace(pM2D, pTVR->getFace(0));
				}
				deleteVert(pM2D, pTVR);
			} else if (pTVR->getCommonFace(pTVL)) {
				while (pTVR->getNumFaces() != 0) {
					deleteFace(pM2D, pTVR->getFace(0));
				}
				deleteVert(pM2D, pTVR);

				while (pTVL->getNumFaces() != 0) {
					deleteFace(pM2D, pTVL->getFace(0));
				}
				deleteVert(pM2D, pTVL);
			} else {
				pTVL->setTfr(createFace(pM2D, pTVR, pTVL));
				pTVR->setTfl(pTVL->getTfr());
			}

			while (pTV->getNumFaces() != 0) {
				deleteFace(pM2D, pTV->getFace(0));
			}
			deleteVert(pM2D, pTV);
		}
	}
	TVL.cullDeleted(iLayer + 1);
}

// Identify and label separate loops of vertices.
int identifyLoops(edamVertList &TVL, const int &iLayer) {
	// Traverse each loop by moving left along the front until we've returned to the vertex we started at.
	int iNumLoops = 0;
	for (edamVert* pTVstart : TVL.getLayer(iLayer + 1)) {
		if (pTVstart->getLoop())
			continue;
		iNumLoops++;

		edamVert* pTV = pTVstart;
		do {
			pTV->setLoop(iNumLoops);
			pTV = pTV->getTfl()->getOppositeVert(pTV);
		} while (pTV != pTVstart);
	}
	return iNumLoops;
}

/*
 * Split the new front at points of self-intersection, creating a set of closed non-intersecting loops.
 * Delete the loops that are not valid parts of the new front (identified using winding numbers, see comment in code).
 */
void clipLoops(Mesh2D* pM2D, edamVertList &TVL, double &dWindInitial,
		int iLayer) {

	// Clip intersecting loops into two different loops (Happens when two
	// different loops intersect or a loop intersects itself).
	clipIntersections(pM2D, TVL, iLayer);

	// Identify and label separate loops of vertices.
	int iNumLoops = identifyLoops(TVL, iLayer);

	// Store the faces belonging to each loop.
	vector<set<edamFace*>> loops(iNumLoops);
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		loops.at(pTV->getLoop() - 1).insert(pTV->getTfl());
	}

	// Delete the loops which have crossed over each other or are a result of a loop crossing over itself
	deleteBadLoops(pM2D, TVL, loops, dWindInitial, iLayer);

	// "Decimate" vertices that we've added at points of self-intersection.
	decimateExtraVerticesAtIntersections(pM2D, TVL, iLayer);
}

/*
 * Create quads between corresponding parent and child edges that have not been deleted, and creates triangles where child edges have been collapsed
 * or where multiple vertices were extruded at corners.
 */
void createGoodCells(Mesh2D* pM2D, edamVertList &TVL, int iLayer) {

	edamVert* pTV1 = NULL;
	edamVert* pTV2 = NULL;
	edamVert* pTV3 = NULL;

	for (edamVert* pTV0 : TVL.getLayer(iLayer)) {
		if (pTV0->getNumKids() == 0)
			continue;

		// Create triangle between two adjacent children and their common parent.
		for (int iK = 0; iK < pTV0->getNumKids() - 1; iK++) {
			pTV1 = pTV0->getKid(iK);
			pTV2 = pTV0->getKid(iK + 1);
			if (pTV1->getTfl()->getOppositeVert(pTV1) == pTV2) {
				createTriCell(pM2D, pTV0, pTV1, pTV2);
			}
		}

		pTV3 = pTV0->getTfl()->getOppositeVert(pTV0);
		if (pTV3->getNumKids() == 0)
			continue;

		// Create quad between parent and child edges, or triangle if the child edge has collapsed (pTV1 == pTV2).
		pTV1 = pTV0->getLastKid();
		pTV2 = pTV3->getFirstKid();
		if (pTV1->getTfl()->getOppositeVert(pTV1) == pTV2 || pTV1 == pTV2) {
			createQuadOrTriCell(pM2D, pTV0, pTV1, pTV2, pTV3);

			// "Fix" a bug where internal boundary faces don't get assigned cells properly.
			// Used only in findDelaunayHull.
			pTV0->getTfl()->qRightCell = true;
		}
	}

	// Remove faces between parents and kids that have no adjacent cells. The faces will *probably* be restored after the Delaunay fill,
	// but this is not guaranteed.
	for (edamVert* pTV0 : TVL.getLayer(iLayer)) {
		if (pTV0->getNumKids() == 0)
			continue;

		for (int iK = 0; iK < pTV0->getNumKids(); iK++) {
			pTV1 = pTV0->getKid(iK);
			edamFace* pTFkid = pTV0->getCommonFace(pTV1);
			assert(pTFkid);

			if (pTFkid->getNumCells() == 0) {
				deleteFace(pM2D, pTFkid);
			}
		}
	}
}

/*
 * Finds faces bounding regions in this layer which need to be triangulated. Depending on their orientation relative to these regions,
 * faces are stored in either bottomFaces or topFaces.
 */
void findDelaunayHull(edamVertList &TVL, int iLayer,
		std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> &bottomFaces,
		std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> &topFaces) {

	for (edamVert* pTV : TVL.getLayer(iLayer)) {
		edamFace* pTFL = pTV->getTfl();
		assert(pTFL);
		if (pTFL->qRightCell)
			continue;

		// Edge in previous layer that does not yet have a second adjacent cell.
		bottomFaces.insert(pTFL);

		for (edamFace* pTF : pTV->getFaceList()) {
			if (pTF->getOppositeVert(pTV)->getLayer() == iLayer + 1) {
				if (pTF->getNumCells() <= 0) {
					cout << "Face with no cells: " << endl;
					cout << pTF->pTV0->x() << ":" << pTF->pTV0->y() << endl;
					cout << pTF->pTV1->x() << ":" << pTF->pTV1->y() << endl;
				}
				assert(pTF->getNumCells() > 0);

				// Edge between parent and kid that does not yet have a second adjacent cell.
				if (!pTF->getRightCell())
					bottomFaces.insert(pTF);
				if (!pTF->getLeftCell())
					topFaces.insert(pTF);
			}
		}

		edamVert* pTVL = pTFL->getOppositeVert(pTV);
		for (edamFace* pTF : pTVL->getFaceList()) {
			if (pTF->getOppositeVert(pTVL)->getLayer() == iLayer + 1) {
				if (pTF->getNumCells() <= 0) {
					cout << "Face with no cells: " << endl;
					cout << pTF->pTV0->x() << ":" << pTF->pTV0->y() << endl;
					cout << pTF->pTV1->x() << ":" << pTF->pTV1->y() << endl;
				}
				assert(pTF->getNumCells() > 0);

				// Edge between parent and kid that does not yet have a second adjacent cell.
				if (!pTF->getRightCell())
					bottomFaces.insert(pTF);
				if (!pTF->getLeftCell())
					topFaces.insert(pTF);
			}
		}
	}
	for (edamVert* pTV : TVL.getLayer(iLayer + 1)) {
		edamFace* pTFL = pTV->getTfl();
		assert(pTFL);
		if (pTFL->getNumCells() > 0)
			continue;

		// Edge in new layer with no adjacent cells.
		topFaces.insert(pTFL);
	}
}

/*
 * Compute a constrained Delaunay triangulation bounded by the input faces. Add the output triangles to the advancing front mesh.
 *
 * This shouldn't be a source of error anymore, now that we provide a closed hull (Sept 22, 2016).
 */
void delaunayFill(Mesh2D* pMorg,
		std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> bottomFaces,
		std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> topFaces) {

	Mesh2D* pMtemp = new Mesh2D(2);
	GRGeom2D* pG2D = new GRGeom2D();

	// Find the set of vertices connected to the input faces.
	std::set<edamVert*> vertSet;
	for (edamFace* pTF : bottomFaces) {
		vertSet.insert(pTF->pTV0);
		vertSet.insert(pTF->pTV1);
	}
	for (edamFace* pTF : topFaces) {
		vertSet.insert(pTF->pTV0);
		vertSet.insert(pTF->pTV1);
	}
	std::vector<edamVert*> vertList(vertSet.begin(), vertSet.end());

	int numVerts = vertList.size();

	// Create index map...
	map<edamVert*, int> realToTemp;
	for (int iV = 0; iV < numVerts; iV++) {
		realToTemp.insert(pair<edamVert*, int>(vertList.at(iV), iV));
	}

	GRPoint* p[numVerts];
	for (int iPoint = 0; iPoint < numVerts; iPoint++) {
		edamVert* pTV = vertList.at(iPoint);
		double adEpsNorm[2] = { 0, 0 };
		if (pTV->getTVcopy()) {
			double adNormL[2];
			double adNormR[2];
			pTV->getTfl()->calcUnitNorm(adNormL);
			pTV->getTfr()->calcUnitNorm(adNormR);
			double adNorm[2] = adSUM2D(adNormL,adNormR);
			NORMALIZE2D(adNorm);

			// Offset duplicate verts slightly to avoid "cannot find seed cell naively" error.
			adEpsNorm[0] = adNorm[0] * 1.e-10;
			adEpsNorm[1] = adNorm[1] * 1.e-10;
		}
		p[iPoint] = pG2D->build_point(pTV->x() + adEpsNorm[0],
				pTV->y() + adEpsNorm[1]);
	}

	for (edamFace* pTF : topFaces) {
		pG2D->build_line(p[realToTemp.at(pTF->pTV0)],
				p[realToTemp.at(pTF->pTV1)], 1, 0, 0, 1);
	}
	for (edamFace* pTF : bottomFaces) {
		pG2D->build_line(p[realToTemp.at(pTF->pTV0)],
				p[realToTemp.at(pTF->pTV1)], 0, 1, 1, 0);
	}

	// Call the constrained Delaunay triangulation routine.
	// pMtemp->setWriteTempMesh(true);
	TriMeshBuilder* TMB = new TriMeshBuilder(pG2D, pMtemp);
	TMB->buildTerribleConstrainedDelaunayTriangulation();
	delete TMB;

	pMtemp->writeTempMesh();

	// Copy the triangles back to the advancing front mesh.
	for (int ii = 0; ii < static_cast<int>(pMtemp->getNumCells()); ii++) {
		Cell *pC = pMtemp->getCell(ii);

		int iV0 = pMtemp->getVertIndex(pC->getVert(0));
		int iV1 = pMtemp->getVertIndex(pC->getVert(1));
		int iV2 = pMtemp->getVertIndex(pC->getVert(2));

		edamVert* pTV0 = vertList.at(iV0);
		edamVert* pTV1 = vertList.at(iV1);
		edamVert* pTV2 = vertList.at(iV2);

		createTriCell(pMorg, pTV0, pTV1, pTV2);
	}

	delete pMtemp;
	delete pG2D;
}

/*
 * Calculte the mean normal of two faces
 */
void calcAdNorm(edamFace* faceL, edamFace* faceR, double adNormL[],
		double adNormR[], double adNorm[]) {
	// calculate the normals to the faces
	faceL->calcUnitNorm(adNormL);
	faceR->calcUnitNorm(adNormR);
	// calculate mean normal which is the unit vector in the direction of the sum of the two normals

	//average of the normals of adjacent faces
	adNorm[0] = adNormL[0] + adNormR[0];
	adNorm[1] = adNormL[1] + adNormR[1];
	NORMALIZE2D(adNorm);
}

//Calculate angle between two faces, given the mean normal and one face
double calcAngleAdjFaces(edamFace* faceL, double adNorm[]) {
	// Calculate angle between adjacent faces.

	double adDirL[2];
	faceL->calcUnitDir(adDirL); // orientation or direction of one of the faces, taken as the left face here

	double dTheta = 2 * GR_acos(dDOT2D(adNorm, adDirL)); // calc dot product of avg. of normals and the left face direction, twice of it is the angle bw..
	return dTheta; // return the angle between the faces
}

void checkEndPointInternalBoundary(edamVert* pTVnL, edamVert* pTVnR,
		double& dTheta, double adNorm[], double adNormL[]) {

	if (pTVnL->getVert() == pTVnR->getVert()) {
		// Endpoint of an internal boundary, must handle as a special case.
		dTheta = 2 * M_PI;
		adNorm[0] = adNormL[1];
		adNorm[1] = -adNormL[0];
	}

}

//Calculate additional number of points to be added between faces with large angle
int calcNumAddPoints(const double& dTheta, const int& iLayer) {
	return (dTheta > M_PI * 1.33 && iLayer == 0 ? ceil((dTheta - M_PI)) - 1 : 0);
}

/*
 * Extrude a single layer of mesh cells from each existing front.
 */
void extrudeAndConnectFaces(Mesh2D* pM2D, edamVertList &TVL, edamVert* pTV,
		edamVert* pTVnR, edamVert* pTVnLkid, edamVert* pTVnRkid,
		edamVert* pTVkidL, edamVert* pTVkidR, double adNormL[],
		double adNormR[], int &iNumAddPoints, double &dTheta,
		double &minDistFromOrigin) {

	// dExtFactor was used to define extrude length as a function of the angle dTheta, trying to smooth out corners more quickly between layers.
	// Not clear that this is useful in practice, sometimes produced unwanted artifacts and instability in front shape.
	//double dExtFactor = (1 + dDOT2D(adNorm,adDirL));
	//double dExtFactor = fmax((1.0 - 0.4*atan(5.0*(dTheta - M_PI))),0.8);
	//double dExtFactor = (dTheta < M_PI ? 1./sin(dTheta/2.0) : 1);
	double dExtFactor = 1;
	double dExtLength = pTV->getExtrudeLength() * dExtFactor;

	// Extrude one vertex along the left face normal direction...
	// SURF: Need 3 coords, and need to project points back onto the curved surface.
	pTVkidL = createVert(pM2D, TVL, pTV->x() + adNormL[0] * dExtLength,
			pTV->y() + adNormL[1] * dExtLength, pTV->getLayer() + 1);
	/*logMessage(MSG_MANAGER, "Left kid:  %d (%6f %6f)\n", pM2D->getVertIndex(pTVkidL->getVert()),
	 pTVkidL->x(), pTVkidL->y());*/
	pTVkidL->setAspectRatio(pTV->getAspectRatio() / growthRatio);
	pTVkidL->setExtrudeLength(pTV->getExtrudeLength() * growthRatio);
	double dist = dMAG2D(pTVkidL->getVert()->getCoords());
	if (minDistFromOrigin > dist)
		minDistFromOrigin = dist;
	// and one along the right face normal.
	// SURF: Need 3 coords, and need to project points back onto the curved surface.
	pTVkidR = createVert(pM2D, TVL, pTV->x() + adNormR[0] * dExtLength,
			pTV->y() + adNormR[1] * dExtLength, pTV->getLayer() + 1);
	/*logMessage(MSG_MANAGER, "Right kid:  %d (%6f %6f)\n", pM2D->getVertIndex(pTVkidR->getVert()),
	 pTVkidR->x(), pTVkidR->y());*/
	pTVkidR->setAspectRatio(pTV->getAspectRatio() / growthRatio);
	pTVkidR->setExtrudeLength(pTV->getExtrudeLength() * growthRatio);

	dist = dMAG2D(pTVkidR->getVert()->getCoords());
	if (minDistFromOrigin > dist)
		minDistFromOrigin = dist;

	// Kids need to be ordered CCW around the parent vertex pTV, so create the face to right one first.
	createFace(pM2D, pTV, pTVkidR);
	/*	logMessage(MSG_MANAGER, "New edge: %d %d\n", pM2D->getVertIndex(pTV->getVert()),
	 pM2D->getVertIndex(pTVkidR->getVert()));*/

	// Create faces from the new kids to their in-layer neighbors.
	if (pTVnLkid) {
		edamFace* pTFnew = createFace(pM2D, pTVkidL, pTVnLkid);
		/*		logMessage(MSG_MANAGER, "New edge: %d %d\n", pM2D->getVertIndex(pTVkidL->getVert()),
		 pM2D->getVertIndex(pTVnLkid->getVert()));*/
		pTVkidL->setTfl(pTFnew);
		pTVnLkid->setTfr(pTFnew);
	}
	if (pTVnRkid) {
		edamFace* pTFnew = createFace(pM2D, pTVnRkid, pTVkidR);
		/*		logMessage(MSG_MANAGER, "New edge: %d %d\n", pM2D->getVertIndex(pTVnRkid->getVert()),
		 pM2D->getVertIndex(pTVkidR->getVert()));*/
		pTVkidR->setTfr(pTFnew);
		pTVnRkid->setTfl(pTFnew);
	}

	// If additional marching directions must be added, space additional kids evenly (in angle) between the left and right ones, in CCW order.
	edamVert* addPoints[iNumAddPoints];
	double dRAngle = atan2(pTVnR->y() - pTV->y(), pTVnR->x() - pTV->x());
	if (iNumAddPoints > 0) {
		for (int iP = 0; iP < iNumAddPoints; iP++) {
			double dExtrudeAngle = dRAngle + 0.5 * M_PI
					+ (dTheta - M_PI) * ((iP + 1.0) / (iNumAddPoints + 1));
			double adPosNew[2] = { pTV->x() + dExtLength * cos(dExtrudeAngle),
					pTV->y() + dExtLength * sin(dExtrudeAngle) };
			// SURF: Need 3 coords, and need to project points back onto the curved surface.

			edamVert* pTVadd = createVert(pM2D, TVL, adPosNew[0], adPosNew[1],
					pTV->getLayer() + 1);
			addPoints[iP] = pTVadd;
			pTVadd->setAspectRatio(pTV->getAspectRatio() / growthRatio);
			pTVadd->setExtrudeLength(pTV->getExtrudeLength() * growthRatio);

			edamVert* pTVaddPrev = (iP == 0 ? pTVkidR : addPoints[iP - 1]);
			edamVert* pTVaddNext = (iP == iNumAddPoints - 1 ? pTVkidL : NULL);

			createFace(pM2D, pTV, pTVadd);
			pTVadd->setTfr(createFace(pM2D, pTVaddPrev, pTVadd));
			pTVaddPrev->setTfl(pTVadd->getTfr());
			if (pTVaddNext) {
				pTVaddNext->setTfr(createFace(pM2D, pTVadd, pTVaddNext));
				pTVadd->setTfl(pTVaddNext->getTfr());
			}
		}
	}
	// Otherwise, simply connect the left and right kids.
	else {
		pTVkidR->setTfl(createFace(pM2D, pTVkidR, pTVkidL));
		pTVkidL->setTfr(pTVkidR->getTfl());
	}

	// Kids need to be ordered CCW around the parent vertex pTV, so create the face to the left one last.
	createFace(pM2D, pTV, pTVkidL);
	/*	logMessage(MSG_MANAGER, "New edge: %d %d\n", pM2D->getVertIndex(pTV->getVert()),
	 pM2D->getVertIndex(pTVkidL->getVert()));*/
}

/*
 * If faces which are bw pTVkidL and pTVnLkid, and pTVkidR and pTVnRkid intersect, delete them and replace
 * them with the faces constructed from the point of intersection and pTVnLkid and pTVnRkid
 */
void dealWithIntersections(Mesh2D* pM2D, edamVertList &TVL, int iLayer) {
	for (edamVert* pTV : TVL.getLayer(iLayer)) {
		if (pTV->getNumKids() < 2)
			continue;

		edamFace* pTFL = pTV->getTfl();
		edamVert* pTVnL = pTFL->getOppositeVert(pTV);
		edamFace* pTFR = pTV->getTfr();
		edamVert* pTVnR = pTFR->getOppositeVert(pTV);
		edamVert* pTVnLkid = NULL;
		if (pTVnL->getNumKids() > 0) {
			pTVnLkid = pTVnL->getFirstKid();
		}
		edamVert* pTVnRkid = NULL;
		if (pTVnR->getNumKids() > 0) {
			pTVnRkid = pTVnR->getLastKid();
		}

		edamVert* pTVkidL = pTV->getLastKid();
		edamVert* pTVkidR = pTV->getFirstKid();

		// Find the two faces which are likely to intersect.
		edamFace* pTFkidL = pTVkidL->getCommonFace(pTVnLkid);
		edamFace* pTFkidR = pTVnRkid->getCommonFace(pTVkidR);
		assert(pTFkidL);
		assert(pTFkidR);

		// If they actually do intersect, replace with a new vertex at the point of intersection.
		double adPos[2];
		if (qIntersect(pTFkidL, pTFkidR, adPos)) {
			while (pTVkidR->getNumFaces() != 0) {
				deleteFace(pM2D, pTVkidR->getFace(0));
			}
			deleteVert(pM2D, pTVkidR);
			while (pTVkidL->getNumFaces() != 0) {
				deleteFace(pM2D, pTVkidL->getFace(0));
			}
			deleteVert(pM2D, pTVkidL);

			edamVert* pTVkidNew = createVert(pM2D, TVL, adPos[0], adPos[1],
					pTV->getLayer() + 1);
			pTVkidNew->setAspectRatio(pTV->getAspectRatio() / growthRatio);
			pTVkidNew->setExtrudeLength(pTV->getExtrudeLength() * growthRatio);

			createFace(pM2D, pTV, pTVkidNew);
			pTVnLkid->setTfr(createFace(pM2D, pTVkidNew, pTVnLkid));
			pTVkidNew->setTfl(pTVnLkid->getTfr());
			pTVkidNew->setTfr(createFace(pM2D, pTVnRkid, pTVkidNew));
			pTVnRkid->setTfl(pTVkidNew->getTfr());

		}
	}
}

// Collapse the short edge between the two kids of a smooth but convex vertex.
void collapseShortEdgeConvexVertex(Mesh2D* pM2D, edamVertList &TVL,
		const int &iLayer) {
	for (edamVert* pTV : TVL.getLayer(iLayer)) {
		if (pTV->getNumKids() != 2)
			continue;

		edamFace* pTFL = pTV->getTfl();
		edamVert* pTVnL = pTFL->getOppositeVert(pTV);
		edamFace* pTFR = pTV->getTfr();
		edamVert* pTVnR = pTFR->getOppositeVert(pTV);
		edamVert* pTVnLkid = NULL;
		if (pTVnL->getNumKids() > 0) {
			pTVnLkid = pTVnL->getFirstKid();
		}
		edamVert* pTVnRkid = NULL;
		if (pTVnR->getNumKids() > 0) {
			pTVnRkid = pTVnR->getLastKid();
		}

		edamVert* pTVkidL = pTV->getLastKid();
		edamVert* pTVkidR = pTV->getFirstKid();

		// Calculate adjacent face normals and the vertex normal (normalized average of face normals).
		double adNormL[2];
		double adNormR[2];
		pTFL->calcUnitNorm(adNormL);
		pTFR->calcUnitNorm(adNormR);
		double adNorm[2] = adSUM2D(adNormL,adNormR);
		NORMALIZE2D(adNorm);
		// Calculate angle between adjacent faces.
		double adDirL[2];
		pTFL->calcUnitDir(adDirL);
		double dTheta = 2 * GR_acos(dDOT2D(adNorm, adDirL));
		if (pTVnL->getVert() == pTVnR->getVert()) {
			// Endpoint of an internal boundary, must handle as a special case.
			dTheta = 2 * M_PI;
			adNorm[0] = adNormL[1];
			adNorm[1] = -adNormL[0];
		}
		if (dTheta < M_PI - 1.e-1 || dTheta > M_PI * 1.33)
			continue;

		//		double dExtFactor = 1;
		//		double dExtLength = pTV->dExtrudeLength * dExtFactor;

		// Combine the two kids.
		double adPos[2] = { 0.5 * (pTVkidR->x() + pTVkidL->x()), 0.5
				* (pTVkidR->y() + pTVkidL->y()) };
		//double adPos[2] = {pTV->x()+dExtLength*adNorm[0],pTV->y()+dExtLength*adNorm[1]};

		while (pTVkidR->getNumFaces() != 0) {
			deleteFace(pM2D, pTVkidR->getFace(0));
		}
		deleteVert(pM2D, pTVkidR);
		while (pTVkidL->getNumFaces() != 0) {
			deleteFace(pM2D, pTVkidL->getFace(0));
		}
		deleteVert(pM2D, pTVkidL);

		edamVert* pTVkidNew = createVert(pM2D, TVL, adPos[0], adPos[1],
				pTV->getLayer() + 1);
		pTVkidNew->setAspectRatio(pTV->getAspectRatio() / growthRatio);
		pTVkidNew->setExtrudeLength(pTV->getExtrudeLength() * growthRatio);

		createFace(pM2D, pTV, pTVkidNew);
		pTVnLkid->setTfr(createFace(pM2D, pTVkidNew, pTVnLkid));
		pTVkidNew->setTfl(pTVnLkid->getTfr());
		pTVkidNew->setTfr(createFace(pM2D, pTVnRkid, pTVkidNew));
		pTVnRkid->setTfl(pTVkidNew->getTfr());
	}
}

// Maintain unit aspect ratio by collapsing short edges in the new front.
void maintainUnitAspectRatio(Mesh2D* pM2D, edamVertList &TVL,
		const int &iLayer) {
	for (edamVert* pTV : TVL.getLayer(iLayer)) {
		if (pTV->getNumKids() == 0)
			continue;

		edamVert* pTVkid = pTV->getLastKid();
		edamVert* pTVkL = pTVkid->getTfl()->getOppositeVert(pTVkid);

		double adDist[2] =
				{ pTVkL->x() - pTVkid->x(), pTVkL->y() - pTVkid->y() };
		double tri_threshold = 0.828; // ~ 2*tan(M_PI/8) minimize normalized max deviation of angle (scaled by 1/90 for quad, 1/60 for tri)
		if (dMAG2D(adDist) < tri_threshold * pTV->getExtrudeLength()) {
			double adPosNew[2] = { pTVkid->x() + 0.5 * adDist[0], pTVkid->y()
					+ 0.5 * adDist[1] };
			pTVkL->setOrgCoords(adPosNew);

			edamVert* pTVkR = pTVkid->getTfr()->getOppositeVert(pTVkid);
			assert(pTVkL);
			assert(pTVkR);
			assert(pTVkL != pTVkR);
			if (pTVkL->getCommonFace(pTVkR)) {
				while (pTVkL->getNumFaces() != 0) {
					deleteFace(pM2D, pTVkL->getFace(0));
				}
				deleteVert(pM2D, pTVkL);
				while (pTVkR->getNumFaces() != 0) {
					deleteFace(pM2D, pTVkR->getFace(0));
				}
				deleteVert(pM2D, pTVkR);
			} else {
				pTVkL->setTfr(createFace(pM2D, pTVkR, pTVkL));
				pTVkR->setTfl(pTVkL->getTfr());

				// Kind of a hack... this keeps the kids in CCW order.
				for (int iP = pTVkid->getNumParents() - 1; iP >= 0; iP--) {
					edamVert* pTVpar = pTVkid->getParent(iP);

					for (int iK = 0; iK < pTVpar->getNumKids(); iK++) {
						if (pTVpar->getKid(iK) == pTVkid) {
							edamFace* pTFdelete = pTVpar->getCommonFace(pTVkid);
							deleteFace(pM2D, pTFdelete);

							createFace(pM2D, pTVpar, pTVkL);
							pTVpar->eraseKids(pTVpar->sizeKids() - 1); //previously pTVpar->kids.erase( pTvpar->kids.end()-1 );
							pTVpar->insertKids(iK, pTVkL); //previously pTVpar->kids.insert( pTVpar->kids.begin() +iK, pTVkL);

							assert(pTVpar->Kidsat(iK) == pTVkL);
						}
					}
				}
			}

			while (pTVkid->getNumFaces() != 0) {
				deleteFace(pM2D, pTVkid->getFace(0));
			}
			deleteVert(pM2D, pTVkid);
		}
	}
	TVL.cullDeleted(iLayer + 1);
}

double extrudeLayer(Mesh2D* pM2D, edamVertList &TVL, double &dWindInitial,
		int iLayer) {
	double minDistFromOrigin = DBL_MAX;
	TVL.addLayer();
	for (edamVert* pTV : TVL.getLayer(iLayer)) {

		// Find neighboring faces and vertices.
		edamFace* pTFL = pTV->getTfl(); //get face on the left of the vertex
		edamVert* pTVnL = pTFL->getOppositeVert(pTV); //get vertex on the left of the vertex
		edamFace* pTFR = pTV->getTfr(); //get face on the right of the vertex
		edamVert* pTVnR = pTFR->getOppositeVert(pTV); //get vertex on the right of the vertex
		edamVert* pTVnLkid = NULL; //declare the first kid of the left vertex and assign it as NULL
		if (pTVnL->getNumKids() > 0) { //assign the first kid of the left vertex to pTVnLkid if it exists
			pTVnLkid = pTVnL->getFirstKid();
		}
		edamVert* pTVnRkid = NULL; //declare the last kid of the right vertex and assign it as NULL
		if (pTVnR->getNumKids() > 0) { //assign the last kid of the right vertex to pTVnRkid if it exists
			pTVnRkid = pTVnR->getLastKid();
		}

		// SURF: Finding marching directions will be different on surfaces.  Need a normal to the edge in the tangent plane.
		// Jasmeet: Also, the normals would have three coordinates now.

		// Calculate angle between adjacent faces.
		double adNormL[2];
		double adNormR[2];
		double adNorm[2];

		calcAdNorm(pTFL, pTFR, adNormL, adNormR, adNorm); // calc mean normal
		double dTheta = calcAngleAdjFaces(pTFL, adNorm); //calculate angle bw adj faces

		// modify dTheta and adNorm if the point (pTV) is an endpoint of an internal boundary
		checkEndPointInternalBoundary(pTVnL, pTVnR, dTheta, adNorm, adNormL);

		// Calculate number of additional points to add at convex corners.
		int iNumAddPoints = calcNumAddPoints(dTheta, iLayer);

		// Instead of extruding, by default, a single kid vertex from every parent, we now extrude a kid edge from every edge. This guarantees that,
		// if the new front is oriented correctly, it does not intersect the previous front.
		do {
			edamVert* pTVkidL = NULL;
			edamVert* pTVkidR = NULL;
			extrudeAndConnectFaces(pM2D, TVL, pTV, pTVnR, pTVnLkid, pTVnRkid,
					pTVkidL, pTVkidR, adNormL, adNormR, iNumAddPoints, dTheta,
					minDistFromOrigin);
		} while (0);
	}
	/* Structure of vertices around pTV
	 *
	 *                 pTVkidL *          * pTVkidR
	 *                          \        /
	 *                           \      /
	 *                            \    /
	 *                             \  /
	 * 	                           * pTV
	 *              pTVnLkid        / \         pTVnRkid
	 *                 *           /   \           *
	 *                  \         /     \         /
	 *                   \       /  	   \       /
	 *                    \     /    	    \     /
	 *                     \   /  already  \   /
	 *                      \ /  meshed or  \ /
	 *                       *    boundary   *
	 *                     pTVnL           pTVnR
	 */

	//loopsVTK(TVL,iLayer);
	// There will usually be a self-intersection between two faces adjacent to concave corners. In this case,
	// delete both kids and replace with one kid at the point of intersection.
	// SURF: The edges don't actually exactly intersect.  This makes detection harder, as well as changing how to fix the problem.
	dealWithIntersections(pM2D, TVL, iLayer);

	// Collapse the short edge between the two kids of a smooth but convex vertex.
	collapseShortEdgeConvexVertex(pM2D, TVL, iLayer);

	TVL.cullDeleted(iLayer + 1); //delete all the vertices from the list which are deleted from this layer

	std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> bottomFaces(
			edamFace::qPointerCompFaceAddress);
	std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> topFaces(
			edamFace::qPointerCompFaceAddress);

	//loopsVTK(TVL,iLayer);
	deleteIsolatedVerts(pM2D, TVL, iLayer);
	//loopsVTK(TVL,iLayer);
	clipLoops(pM2D, TVL, dWindInitial, iLayer);
	//loopsVTK(TVL,iLayer);

	// Maintain unit aspect ratio by collapsing short edges in the new front.
	maintainUnitAspectRatio(pM2D, TVL, iLayer);

	//loopsVTK(TVL,iLayer);

	deleteBadVerts(pM2D, TVL, iLayer);
	createGoodCells(pM2D, TVL, iLayer);

	//loopsVTK(TVL,iLayer);

	findDelaunayHull(TVL, iLayer, bottomFaces, topFaces);

	if (!bottomFaces.empty() || !topFaces.empty()) {
		vector<edamFace*> VTKfaces;
		VTKfaces.insert(VTKfaces.end(), bottomFaces.begin(), bottomFaces.end());
		VTKfaces.insert(VTKfaces.end(), topFaces.begin(), topFaces.end());

		VTKFaceList(VTKfaces, "DelaunayHull2D.vtk");
		// SURF: How we fill these bits where we stopped marching will have to be
		// different.  Probably steal this from the 3D minimum area triangulation code.
		delaunayFill(pM2D, bottomFaces, topFaces);
	}

	//pM2D->writeTempMesh();
	return minDistFromOrigin;
}

void adROT2D(double* adVec, double dAngle) {
	double dTemp = adVec[0] * cos(dAngle) - adVec[1] * sin(dAngle);
	adVec[1] = adVec[0] * sin(dAngle) + adVec[1] * cos(dAngle);
	adVec[0] = dTemp;
}

/*
 * Compute "forces" acting on a vertex based on neighborhood vertices and cells.
 *
 * This smoothing can tangle the mesh, mostly happens when two opposing fronts are very close together.
 * Will have trouble merging nearly-parallel opposing fronts until this is fixed.
 *
 * ***Issue currently suppressed by crippling smoothing of front vertices. This probably has negative effects on mesh quality.
 */

//SURF:  Smoothing will have to be done differently.  At the very least, new point locations will have to be constrained
// to lie on the surface.
void forceModelSmooth(Mesh2D* pM2D, edamVertList &TVL, int iEndLayer,
		int iStartLayer, int iNumPasses, double k_angle, double k_extru,
		double k_neigh, double k_org) {

	iStartLayer = max(iStartLayer, 1);

	for (int iPass = 0; iPass < iNumPasses; iPass++) {
		for (int iLayer = iStartLayer; iLayer <= iEndLayer; iLayer++) {
			for (edamVert* pTV : TVL.getLayer(iLayer)) {
				double adForce[2] = { 0, 0 };

				double adNormL[2];
				double adNormR[2];
				if (!pTV->getTfl()) {
					pM2D->writeTempMesh();
					cout << pTV->x() << ":" << pTV->y() << endl;
				}
				pTV->getTfl()->calcUnitNorm(adNormL);
				pTV->getTfr()->calcUnitNorm(adNormR);
				double adNorm[2] =
						adSUM2D(adNormL,adNormR);
						NORMALIZE2D(adNorm);

						double adDirL[2];
						pTV->getTfl()->calcUnitDir(adDirL);

						// dExtFactor was used to define extrude length as a function of the angle dTheta, trying to smooth out corners more quickly between layers.
						// Not clear that this is useful in practice, sometimes produced unwanted artifacts and instability in front shape.
						//double dTheta = 2*GR_acos(dDOT2D(adNorm,adDirL));
						//double dExtFactor = (1 + dDOT2D(adNorm,adDirL));
						//double dExtFactor = fmax((1.0 - 0.4*atan(5.0*(dTheta - M_PI))),0.8);
						//double dExtFactor = (dTheta < M_PI ? 1./sin(dTheta/2.0) : 1);
								double dExtFactor = 1;
								double dExtLength = pTV->getExtrudeLength() * dExtFactor;

								double dMinFLength = LARGE_DBL,
				dMaxFLength = 0;
				for (int iF = 0; iF < pTV->getNumFaces(); iF++) {
					edamVert* pTVother = pTV->getFace(iF)->getOppositeVert(pTV);
					double adDistOther[2] = { pTVother->x() - pTV->x(),
							pTVother->y() - pTV->y() };
					double dLengthFace = dMAG2D(adDistOther);
					dMinFLength = fmin(dMinFLength, dLengthFace);
					dMaxFLength = fmax(dMaxFLength, dLengthFace);
				}

				// Force node towards its original position.
				do {
					double dAspectRatio1 = dMaxFLength / dMinFLength;
					double dAspectRatio = fmax(dAspectRatio1,
							pTV->getAspectRatio());

					double dAspectFactor = (
							dAspectRatio > 1 ?
									(dAspectRatio < 21 ? 0.5 * (1 - cos(
									M_PI * (dAspectRatio - 1) / 20)) :
															1) :
									0);
					adForce[0] += k_org * dAspectFactor
							* (pTV->getOrigX() - pTV->x()) / dExtLength;
					adForce[1] += k_org * dAspectFactor
							* (pTV->getOrigY() - pTV->y()) / dExtLength;
				} while (0);

				// Force edge length to kid node closer to ideal.
				for (edamFace* pTFkid : pTV->getFaceList()) {
					edamVert* pTVkid = pTFkid->getOppositeVert(pTV);
					if (pTVkid->getLayer() < pTV->getLayer())
						continue;
					if (pTFkid == pTV->getTfl() || pTFkid == pTV->getTfr())
						continue;

					double dLengthIdeal = dExtLength;
					double adDistKid[2] = { pTVkid->x() - pTV->x(), pTVkid->y()
							- pTV->y() };
					double dLengthReal = dMAG2D(adDistKid);
					double dForceMag = k_extru * (dLengthReal - dLengthIdeal)
							/ (dLengthReal + dLengthIdeal);
					if (pTV->getNumParents() > 1
							|| pTVkid->getLayer() == pTV->getLayer()) {
						NORMALIZE2D(adDistKid);
						adForce[0] += dForceMag * adDistKid[0];
						adForce[1] += dForceMag * adDistKid[1];
					} else {
						adForce[0] += dForceMag * adNorm[0];
						adForce[1] += dForceMag * adNorm[1];
					}
				}

				// Force edge length to parent node closer to ideal.
				double dAvgParentDist[2] = { 0, 0 };
				for (edamVert* pTVparent : pTV->getParents()) {
					double dLengthIdeal = dExtLength / growthRatio;
					double adDistPar[2] = { pTVparent->x() - pTV->x(),
							pTVparent->y() - pTV->y() };
					double dLengthReal = dMAG2D(adDistPar);
					double dForceMag = k_extru * (dLengthReal - dLengthIdeal)
							/ (dLengthReal + dLengthIdeal);
					NORMALIZE2D(adDistPar);
					if (pTV->getNumParents() > 1) {
						adForce[0] += dForceMag * adDistPar[0];
						adForce[1] += dForceMag * adDistPar[1];
					} else {
						adForce[0] -= dForceMag * adNorm[0];
						adForce[1] -= dForceMag * adNorm[1];
					}
					dAvgParentDist[0] += adDistPar[0];
					dAvgParentDist[1] += adDistPar[1];
				}
				NORMALIZE2D(dAvgParentDist);

				do {
					// Force node to the middle of its two neighbors.
					edamVert* pTVL = pTV->getTfl()->getOppositeVert(pTV);
					double adDistL[2] = { pTVL->x() - pTV->x(), pTVL->y()
							- pTV->y() };
					double dLengthL = dMAG2D(adDistL);
					edamVert* pTVR = pTV->getTfr()->getOppositeVert(pTV);
					double adDistR[2] = { pTVR->x() - pTV->x(), pTVR->y()
							- pTV->y() };
					double dLengthR = dMAG2D(adDistR);
					double dForceMag = k_neigh * (dLengthL - dLengthR)
							/ (dLengthL + dLengthR);

					adForce[0] += dForceMag * dAvgParentDist[1];
					adForce[1] -= dForceMag * dAvgParentDist[0];
				} while (0);

				// Force angles to 90 and 60 degrees for quads and triangles respectively.
				for (int iF = 0; iF < pTV->getNumFaces(); iF++) {
					edamFace* pTF = pTV->getFace(iF);
					Cell* pC = pTF->getLeftCell(pTV);
					if (!pC)
						continue;

					double dCellDihed = pTV->calcCellDihed(pC);
					double dMoveDir[2];
					pTF->calcUnitDir(dMoveDir, pTV);
					adROT2D(dMoveDir, 0.5 * dCellDihed);

					double dAngleFactor;
					if (pC->getNumFaces() == 3) {
						dAngleFactor = (M_PI / 3.0 - dCellDihed) / (M_PI / 3.0);
					} else {
						dAngleFactor = (M_PI / 2.0 - dCellDihed) / (M_PI / 2.0);
					}

					for (int ii = 0; ii < 2; ii++) {
						adForce[ii] += k_angle * dMoveDir[ii] * dAngleFactor;
					}
				}

				// Move node based on computed net force.
				double dSmoothScale = fmin(dExtLength, dMinFLength);
				adForce[0] *= dSmoothScale;
				adForce[1] *= dSmoothScale;
				double adPosNew[2] = { pTV->x() + adForce[0], pTV->y()
						+ adForce[1] };

				// Partial solution to mesh invalidation - prevent movement of front vertices from extending meshed region.
				if (!pTV->getTfl()->getRightCell()) {
					bool qBadR =
							(checkOrient2D(pTV->getVert()->getCoords(),
									pTV->getTfr()->getOppositeVert(pTV)->getVert()->getCoords(),
									adPosNew) == 1);
					bool qBadL =
							(checkOrient2D(pTV->getVert()->getCoords(),
									pTV->getTfl()->getOppositeVert(pTV)->getVert()->getCoords(),
									adPosNew) == -1);
					if (qBadR && qBadL) {
						adPosNew[0] = pTV->x();
						adPosNew[1] = pTV->y();
					} else if (qBadR) {
						double dDotNormR = dDOT2D(adForce, adNormR);
						assert(dDotNormR > 0);
						for (int ii = 0; ii < 2; ii++) {
							adPosNew[ii] -= adNormR[ii] * dDotNormR;
						}
					} else if (qBadL) {
						double dDotNormL = dDOT2D(adForce, adNormL);
						assert(dDotNormL > 0);
						for (int ii = 0; ii < 2; ii++) {
							adPosNew[ii] -= adNormL[ii] * dDotNormL;
						}
					}
				}

				// Set the new position of the vertex
				pTV->setCoords(adPosNew);
			}
		}
	}
}

/*
 * Set boundary conditions for the completed mesh.
 *
 * TODO: Currently assumes walls at all initial boundaries and far-field at open fronts.
 * 		 Probably want to use BC info from input mesh instead.
 */
void setBCs(Mesh2D* pM2D, edamVertList &TVL) {

	int iFirstLayerSize = TVL.getLayer(0).size();

	int iNumBFs = pM2D->getNumBdryFaces();
	for (int iBF = 0; iBF < iNumBFs; iBF++) {
		BFace* pBF = pM2D->getBFace(iBF);

		int iV = pM2D->getVertIndex(pBF->getVert(0));
		if (iV < iFirstLayerSize) {
			dynamic_cast<BdryEdge*>(pBF)->setBdryCondition(1);
		} else {
			dynamic_cast<BdryEdge*>(pBF)->setBdryCondition(5);
		}
	}
}

/*
 * Sets the initial marching thickness for each boundary vertex.
 *
 * Specify point-wise lengths in a file (and interpolate) using '-X fileName' (not thoroughly tested...) or a global value using '-x length'.
 */
void initExtrudeLength(Mesh2D* /*pM2D*/, edamVertList &TVL,
		const char* strExtLenFileName = NULL) {

	// If a length scale file is provided, create a length2D object to interpolate the given values.
	if (strExtLenFileName) {

		Mesh2D* pMlen = new Mesh2D(2);
		GRGeom2D* pG2D = new GRGeom2D();
		vector<double> lengthList;

		FILE *pFLenFile = fopen(strExtLenFileName, "r");
		if (NULL == pFLenFile)
			vFatalError("Couldn't open file for reading", "length scale input");

		// Bounding box 5000 chords away from the origin. This is bad because it will break on geometries with sufficiently large scale (which is arbitrary).
		double dBBox = 5000;
		GRPoint* P0 = pG2D->build_point(-dBBox, -dBBox);
		GRPoint* P1 = pG2D->build_point(dBBox, -dBBox);
		GRPoint* P2 = pG2D->build_point(dBBox, dBBox);
		GRPoint* P3 = pG2D->build_point(-dBBox, dBBox);

		pG2D->build_line(P0, P1, 1, 0, 0, 1);
		pG2D->build_line(P1, P2, 1, 0, 0, 1);
		pG2D->build_line(P2, P3, 1, 0, 0, 1);
		pG2D->build_line(P3, P0, 1, 0, 0, 1);

		// Read data from file.
		double x_loc = -LARGE_DBL, y_loc = -LARGE_DBL, length = -LARGE_DBL;
		while (1) {
			int iDataRead = fscanf(pFLenFile, "%lf%lf%lf", &(x_loc), &(y_loc),
					&(length));
			if (iDataRead != 3)
				break;

			pG2D->build_point(x_loc, y_loc);
			lengthList.push_back(length);
		}

		// Construct the triangulation, needed for interpolation with barycentric coordinates.
		TriMeshBuilder* TMB = new TriMeshBuilder(pG2D, pMlen);
		TMB->buildConstrainedDelaunayTriangulation();
		delete TMB;

		//writeVTKLegacy(*pMlen,"extLenMesh");

		int jj = 0;
		for (int ii = 0; ii < static_cast<int>(pMlen->getNumVerts()); ii++) {
			Vert* pV = pMlen->getVert(ii);

			// Set length scales on the bounding box to the global extrude length specified with -x.
			//
			// Approximately, initial mesh points outside the convex hull of the specified length points will
			// get a length scale linearly interpolated between two nearby points on the hull.
			if (fabs(pV->x()) + fabs(pV->y()) >= dBBox * 0.9) {
				pV->setLengthScale(m_extrudeLength);
			}
			// Set length scales at points specified in the file.
			else {
				pV->setLengthScale(lengthList.at(jj));
				jj++;
			}
		}
		{
			// Initialize the Length2D object, including storing length scale in pMLen's verts.
			std::shared_ptr<GRUMMP::Length> len = GRUMMP::Length::Create(pMlen,
					1, 1);
			// The length scale was created only for the side effect
			// of storing LS values on all the verts in the mesh.
		}

		fclose(pFLenFile);
	}

	for (edamVert* pTV : TVL.getLayer(0)) {
		int iNumFaces = pTV->getNumFaces();

		// Calculate longest adjacent face to each vertex to assign local "aspect ratio" values.
		double dMaxFaceLength = 0, dMinFaceLength = LARGE_DBL;
		for (int iF = 0; iF < iNumFaces; iF++) {
			edamVert* pTVopp = pTV->getFace(iF)->getOppositeVert(pTV);
			double adDist[2] =
					{ pTVopp->x() - pTV->x(), pTVopp->y() - pTV->y() };
			dMaxFaceLength = fmax(dMaxFaceLength, dMAG2D(adDist));
			dMinFaceLength = std::min(dMinFaceLength, dMAG2D(adDist));
		}

		double extrudeLen = std::min(dMinFaceLength, m_extrudeLength);
		if (strExtLenFileName) {
			extrudeLen = std::min(extrudeLen, pTV->getVert()->getLengthScale());
		}
		pTV->setExtrudeLength(extrudeLen);
		pTV->setAspectRatio(dMaxFaceLength / extrudeLen);
	}
}

/*
 * Return the quality of the resulting quad after combining the pair of triangles adjacent to this face.
 * Currently: decreases with maximum absolute value of (internal angle - PI/2) for resulting quad.
 */
double qCombineQuality(edamFace* pTF) {
	assert(pTF->getLeftCell()->getNumFaces() == 3);
	assert(pTF->getRightCell()->getNumFaces() == 3);

	double maxDev = 0;
	Cell* pCR = pTF->getRightCell();
	Cell* pCL = pTF->getLeftCell();

	edamVert* pTV0 = pTF->pTV0;
	edamVert* pTV1 = pTF->getOppositeVert(pCR);
	edamVert* pTV2 = pTF->pTV1;
	edamVert* pTV3 = pTF->getOppositeVert(pCL);

	double dihed0 = pTV0->calcCellDihed(pCR) + pTV0->calcCellDihed(pCL);
	double dihed1 = pTV1->calcCellDihed(pCR);
	double dihed2 = pTV2->calcCellDihed(pCR) + pTV2->calcCellDihed(pCL);
	double dihed3 = pTV3->calcCellDihed(pCL);

	maxDev = fmax(maxDev, fabs(dihed0 - M_PI_2));
	maxDev = fmax(maxDev, fabs(dihed1 - M_PI_2));
	maxDev = fmax(maxDev, fabs(dihed2 - M_PI_2));
	maxDev = fmax(maxDev, fabs(dihed3 - M_PI_2));

	return -maxDev;
}
/*
 * Returns true if removing face pTF0 obtains a better resulting quad than removing pTF1.
 * (Yes, this is a "greater than" comparator, because we need to sort in descending order of quality).
 */
bool qCombineBadnessComp(edamFace* pTF0, edamFace* pTF1) {
	return qCombineQuality(pTF0) > qCombineQuality(pTF1);
}

/*
 * For verts between iStartLayer and iEndLayer, combine adjacent triangular faces to quads where possible.
 * Currently combining in descending order of resulting quad quality, can stop combining after a minimum quality threshold is reached.
 */
void combineTrisToQuads(Mesh2D* pM2D, edamVertList &TVL, int iStartLayer,
		int iEndLayer, double dMinQual = -LARGE_DBL) {
	std::set<edamFace*, bool (*)(const edamFace*, const edamFace*)> faceSet(
			edamFace::qPointerCompFaceAddress);

	// Construct a set of candidate pairs of triangles (defined by their common face) for merging (deletion of the face).
	for (int iLayer = iStartLayer; iLayer < iEndLayer; iLayer++) {
		for (edamVert* pTV : TVL.getLayer(iLayer)) {
			for (int iF = 0; iF < pTV->getNumFaces(); iF++) {

				edamFace* pTF = pTV->getFace(iF);
				if (pTF->getNumCells() == 2) {
					edamVert* pTV0 = pTF->getOppositeVert(pTV);

					// For now, don't delete faces which define the interface between layers.
					// Don't really need this, but it allows easy traversal of loops.
					if (pTV->getTfl() == pTF || pTV->getTfr() == pTF
							|| pTV0->getTfl() == pTF || pTV0->getTfr() == pTF)
						continue;

					Cell* pC0 = pTF->getLeftCell();
					Cell* pC1 = pTF->getRightCell();
					if (pC0->getNumFaces() == 3 && pC1->getNumFaces() == 3) {
						faceSet.insert(pTF);
					}
				}
			}
		}
	}

	// Sort in descending order of resulting quad quality.
	std::vector<edamFace*> faceList(faceSet.begin(), faceSet.end());
	std::sort(faceList.begin(), faceList.end(), qCombineBadnessComp);

	for (edamFace* pTF : faceList) {
		assert(pTF->pF);
		Cell* pC0 = pTF->getLeftCell();
		Cell* pC1 = pTF->getRightCell();

		// Check that both adjacent faces are still triangles.
		if (pC0->getNumFaces() != 3)
			continue;
		if (pC1->getNumFaces() != 3)
			continue;

		// Finish when remaining merges would produce quads below the quality threshold.
		if (qCombineQuality(pTF) < dMinQual)
			break;

		// Find cells and vertices that make up the quad:
		//
		//					    pTV2
		//					   / ||	\ .
		//					  /  ||	 \ .
		//					 /   ||	  \ .
		//					/    ||	   \ .
		//				  pTV3   pTF  pTV1
		//					\    ||	   /
		//					 \   ||	  /
		//					  \  ||	 /
		//					   \ || /
		//					    pTV0
		//
		Cell* pCR = pTF->getRightCell();
		Cell* pCL = pTF->getLeftCell();
		edamVert* pTV0 = pTF->pTV0;
		edamVert* pTV1 = NULL;
		for (int iF = 0; iF < pTV0->getNumFaces(); iF++) {
			edamFace* pTFv1 = pTV0->getFace(iF);
			if (pTFv1 != pTF && pTFv1->pF->hasCell(pCR))
				pTV1 = pTFv1->getOppositeVert(pTV0);
		}
		edamVert* pTV2 = pTF->pTV1;
		edamVert* pTV3 = NULL;
		for (int iF = 0; iF < pTV0->getNumFaces(); iF++) {
			edamFace* pTFv3 = pTV0->getFace(iF);
			if (pTFv3 != pTF && pTFv3->pF->hasCell(pCL))
				pTV3 = pTFv3->getOppositeVert(pTV0);
		}

		// Don't merge two triangles if both are on the boundary. This can easily lead to concave quads which can't be fixed by smoothing.
		if (pTV1->getCommonFace(pTV2)->qBdry
				&& pTV2->getCommonFace(pTV3)->qBdry) {
			continue;
		}

		// Update the connectivity.
		deleteCell(pM2D, pC0);
		deleteCell(pM2D, pC1);
		deleteFace(pM2D, pTF);
		createQuadCell(pM2D, pTV0, pTV1, pTV2, pTV3);
	}
}

/*
 * Remove all mesh cells at (x > xCutoff). Using this to remove the far ends of wake lines and trim the cell count.
 * Jasmeet: not being used anywhere else.
 */
void vRemoveTail(Mesh2D* pM2D, double xCutoff) {
	for (size_t iV = 0; iV < pM2D->getNumVerts(); iV++) {
		Vert* pV = pM2D->getVert(iV);

		if (pV->x() > xCutoff) {
			for (int iF = 0; iF < pV->getNumFaces(); iF++) {
				Face* pF = pV->getFace(iF);
				Cell* pC;
				if ((pC = pF->getRightCell())) { // Intentional assignment in condition.
					pM2D->deleteCell(pC);
				}
				if ((pC = pF->getLeftCell())) { // Intentional assignment in condition.
					pM2D->deleteCell(pC);
				}
				pM2D->deleteFace(pF);
			}
			pM2D->deleteVert(pV);
		}
	}
}

SurfMeshBuilder* initializeSurfaceMesh(char * strInFileName, int iNArg,
		char *apcArgs[]) {
	SurfMeshBuilder* SMB = NULL;

	// Start CGM instance and initialize the arguments for the import_facets function
	// import_facets function details the meaning of each of these arguments

	CGMApp::instance()->startup(iNArg, apcArgs);
	CubitBoolean use_feature_angle = CUBIT_TRUE;
	double feature_angle = 135;
	double tolerance = 1.e-6;
	int interp_order = 4;
	CubitBoolean smooth_non_manifold = CUBIT_TRUE;
	CubitBoolean split_surfaces = CUBIT_FALSE;
	CubitBoolean stitch = CUBIT_TRUE;
	CubitBoolean improve = CUBIT_TRUE;
	DLIList<CubitQuadFacet*> quad_facet_list;
	DLIList<CubitFacet*> tri_facet_list;
	DLIList<Surface*> surface_list;

	// Build geometry from a faceted representation by calling CGM

	GeometryQueryTool *gqt = GeometryQueryTool::instance();
	FacetQueryEngine *fqe = FacetQueryEngine::instance();
	fqe->import_facets(strInFileName, use_feature_angle, feature_angle,
			tolerance, interp_order, smooth_non_manifold, split_surfaces,
			stitch, improve, quad_facet_list, tri_facet_list, surface_list,
			STL_FILE);

	// Set names of output files carrying faces and edges
	GraphicsOut go_surf("CGMsurf.vtk");
	GraphicsOut go_curv("CGMcurv.vtk");

	// Write faces and edges in VTK format
	go_surf.output_faces();
	go_curv.output_edges();

	//Set the compare tolerance to a lower value for all FacetSurface.
	//This is needed to project points onto the surface interpolated from facets.

	// this is here for now, because its not set up to handle multiple volumes,
	// but the capability now exists

	assert(gqt->num_ref_volumes() == 1);
	RefVolume * refVolume = gqt->get_first_ref_volume();
	DLIList<RefFace*> all_surfaces;
	refVolume->ref_faces(all_surfaces);
	int num_surfaces = all_surfaces.size();
	RefFace* surface;
	FacetSurface* facet_surface;
	for (int i = 0; i < num_surfaces; ++i) {
		surface = all_surfaces.get_and_step();
		facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());
		assert(facet_surface);
		facet_surface->get_eval_tool()->compare_tol(1.e-6);
	}

	double allowed_TVT = 0.45 * M_PI;

	SMB = new SurfMeshBuilder(refVolume, SurfMeshBuilder::FACET, allowed_TVT);
	SMB->sample_all_surfaces();
	SMB->output_stitched_surface("final-sampled.vtk");

	return SMB;
}

int main(int iNArg, char *apcArgs[]) {

	// Input and output file names
	char strOutFileName[FILE_NAME_LEN];
	char strInFileName[FILE_NAME_LEN];

	// These must be set in command line.
	m_extrudeLength = 0;
	growthRatio = 0;

	// Number of layers to extrude
	numLayers = 50;

	// Next option character returned by getopt.
	// Returns '-1' when all command line options are parsed, '?' when unknown option detected
	int iPChar;

	// Read the command line input
	while ((iPChar = getopt(iNArg, apcArgs, "i:o:")) != EOF) {
		switch (iPChar) {
		case 'i': {
			strncpy(strInFileName, optarg, FILE_NAME_LEN - 1); // File name for input .bdry file.
			break;
		}
		case 'o': {
			strncpy(strOutFileName, optarg, FILE_NAME_LEN - 1); // File name for final output mesh.
			break;
		}
		default: {
			break;
		}
		}
	}

	edamVertList TVL;
	Mesh2D* pM2D = NULL;

	// parse input stl file
	if (*strInFileName) {
		auto info = stl::parse_stl(strInFileName);
		vector<stl::triangle> triangles = info.triangles;
		cout << "STL HEADER = " << info.name << std::endl;
		cout << "# triangles = " << triangles.size() << std::endl;
	} else {
		cout << "No input file name.";
	}

	SurfMeshBuilder *SMB = initializeSurfaceMesh(strInFileName, iNArg, apcArgs);

//			std::vector<CubitVector> pointList;
//			for(int iB = 0; iB < pVMesh->getNumBdryFaces(); iB++){
//				TriBFace * bface = dynamic_cast<TriBFace*>(pVMesh->getBFace(iB));
//				double split[3];
//				bface->calcSplitPoint(split);
//				pointList.push_back(split);
//
//			}
//			writeVTKPointSet(pointList,NULL,"projections","");
//			abort();

	/*numLoopsVTK = 0;

	 // File names, etc. Input/output file names must be set in command line, length scale files (etc.) are optional.
	 int iPChar;
	 char strOutFileName[FILE_NAME_LEN];
	 char strInFileName[FILE_NAME_LEN];
	 char strInMeshFileName[FILE_NAME_LEN];
	 char strQualFileName[FILE_NAME_LEN];
	 char strTriLenFileName[FILE_NAME_LEN];
	 bool qTriLenGiven = false;
	 char strExtLenFileName[FILE_NAME_LEN];
	 bool qExtLenGiven = false;
	 bool qGRmesh = false;
	 bool qGenerateBdrySpacing = false;
	 bool qReadBdryFromMeshFile = false;
	 bool qMakeMeshAllTris = false;

	 // Inputs to 'tri' executable used to generate the initial mesh. Can be optionally set in command line, otherwise use these defaults.
	 double dTriReso = 1;
	 double dTriGrad = 1;

	 // Smoothing parameters. Can be optionally set in command line, otherwise use these defaults.
	 double k_angle = 0.08;
	 double k_extru = 0.2;
	 double k_neigh = 0.1;
	 double k_org = 1.0;

	 // Initial Winding Number, critical in deciding which loops to delete and which ones to keep. Currently set for exterior mesh,
	 // but doesn't matter as it is assigned later
	 double dWindInitial = 0.0;

	 // These must be set in command line.
	 m_extrudeLength = 0;
	 growthRatio = 0;

	 numLayers = 50;

	 while ((iPChar = getopt(
	 iNArg, apcArgs,
	 "b:g:G:i:I:l:n:o:R:stx:X:0:1:2:3:"))
	 != EOF) {
	 switch (iPChar) {

	 case 'b':
	 {
	 int i;
	 sscanf(optarg, "%d", &i); // Use '-b 1' if the mesh generated by tri is a .grmesh file. Happens when internal boundaries exist in the mesh.
	 if (i == 0) {
	 qGRmesh = false;
	 }
	 else if (i == 1) {
	 qGRmesh = true;
	 }
	 else {
	 assert(0);
	 }
	 break;
	 }
	 case 'g':
	 {
	 sscanf(optarg, "%lf", &growthRatio); // Off-wall growth ratio per layer of cells.
	 break;
	 }
	 case 'G':
	 {
	 sscanf(optarg, "%lf", &dTriGrad); // '-g' input to tri.
	 break;
	 }
	 case 'i':
	 {
	 strncpy(strInFileName, optarg, FILE_NAME_LEN - 1); // File name for input .bdry file.
	 break;
	 }
	 case 'I':
	 {
	 strncpy(strInMeshFileName, optarg, FILE_NAME_LEN - 1);
	 qGenerateBdrySpacing = false;
	 qReadBdryFromMeshFile = true;
	 break;
	 }
	 case 'l':
	 {
	 qTriLenGiven = true;
	 strncpy(strTriLenFileName, optarg, FILE_NAME_LEN - 1); // Length scale file for generation of initial mesh in tri.
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
	 case 'R':
	 {
	 sscanf(optarg, "%lf", &dTriReso); // '-r' input to tri.
	 break;
	 }
	 case 's':
	 {
	 // Use an external call to tri to generate the bdry discretization.
	 qGenerateBdrySpacing = true;
	 qReadBdryFromMeshFile = true;
	 break;
	 }
	 case 't':
	 {
	 qMakeMeshAllTris = true;
	 break;
	 }
	 case 'x':
	 {
	 sscanf(optarg, "%lf", &m_extrudeLength); // Off-wall spacing.
	 break;
	 }
	 case 'X':
	 {

	 int i;
	 sscanf(optarg, "%d", &i); // Use '-b 1' if the mesh generated by tri is a .grmesh file. Happens when internal boundaries exist in the mesh.
	 if (i == 0) {            Shahzaib for .bdry file with 0012 airfoil and streamline at TE/LE...if you want to treat as internal boundary (that r 1 r 1 in .bdry file)...use -b1
	 qGRmesh = false;
	 }
	 else if (i == 1) {
	 qGRmesh = true;
	 }
	 else {
	 assert(0);
	 }

	 qExtLenGiven = true;
	 strncpy(strExtLenFileName, optarg, FILE_NAME_LEN - 1); // Off-wall spacing (specified as an interpolation mesh to allow different spacings at different parts of the boundary).

	 break;
	 }
	 case '0':
	 {
	 sscanf(optarg, "%lf", &k_angle);
	 break;
	 }
	 case '1':
	 {
	 sscanf(optarg, "%lf", &k_extru);
	 break;
	 }
	 case '2':
	 {
	 sscanf(optarg, "%lf", &k_neigh);
	 break;
	 }
	 case '3':
	 {
	 sscanf(optarg, "%lf", &k_org);
	 break;
	 }
	 default:
	 {
	 break;
	 }
	 }
	 }

	 edamVertList TVL;
	 Mesh2D* pM2D = NULL;

	 // Run GRUMMP's "tri" to generate a boundary discretization.
	 if (qGenerateBdrySpacing) {
	 meshFileFromBdryFile(strInFileName, dTriReso, dTriGrad,
	 (qTriLenGiven ? strTriLenFileName : NULL));
	 sprintf(strInMeshFileName, "%s", strInFileName);
	 }

	 if (qReadBdryFromMeshFile) {
	 if (qGRmesh) {
	 sprintf(strInMeshFileName + strlen(strInMeshFileName), ".grmesh");
	 }
	 else {
	 sprintf(strInMeshFileName + strlen(strInMeshFileName), ".mesh");
	 }
	 // Load the boundary discretization.
	 pM2D = initMeshFromMeshFile(strInMeshFileName, TVL, true);
	 }
	 else {
	 // Generate the mesh bdry spacing directly in the mesh.
	 pM2D = new Mesh2D(2);
	 GRGeom2D* pG2D = new GRGeom2D(strInFileName);
	 pM2D->setEncroachmentType(eBall);
	 {
	 TriMeshBuilder* TMB = new TriMeshBuilder(pG2D, pM2D);
	 TMB->buildConstrainedDelaunayTriangulation();
	 delete TMB;
	 }
	 char strSurfInsParams[] = "";
	 smoothBdryRefine(pM2D, pG2D, dTriGrad, dTriReso, strSurfInsParams);

	 initFromExistingMesh(pM2D, TVL, true);
	 }

	 // pM2D->setWriteTempMesh(true);
	 initExtrudeLength(pM2D, TVL, (qExtLenGiven ? strExtLenFileName : NULL));

	 time_t t = clock();
	 double minDistFromOrigin = 0, requestedMinDist = 5000;
	 for (int ii = 0; ii < numLayers && minDistFromOrigin < requestedMinDist;
	 ii++) {
	 printf("%d\n", ii);

	 // Generate layer of cells.
	 minDistFromOrigin = extrudeLayer(pM2D, TVL, dWindInitial, ii);
	 // Smooth.
	 forceModelSmooth(pM2D, TVL, ii + 1, ii - 3, 5, k_angle, k_extru,
	 k_neigh, k_org);
	 // Combine triangles to quads where possible.
	 combineTrisToQuads(pM2D, TVL, max(ii, 2), ii + 2);

	 // Smooth.
	 forceModelSmooth(pM2D, TVL, ii + 1, ii - 3, 5, k_angle, k_extru,
	 k_neigh, k_org);
	 }
	 t = clock() - t;
	 printf("Time: %f seconds\n", t / 1000000.0);

	 //combineTrisToQuads(pM2D,TVL,1,2);

	 pM2D->purgeAllEntities();
	 setBCs(pM2D, TVL);

	 if (qMakeMeshAllTris) {
	 pM2D->makeSimplicial();
	 pM2D->purgeAllEntities();
	 GRUMMP::MinMaxAngleSwapDecider2D MMASD2D;
	 GRUMMP::SwapManager2D SM2D(&MMASD2D, pM2D);
	 SM2D.swapAllFaces();
	 GRUMMP::OptMSSmoothingManager* OMS = GRUMMP::OptMSSmoothingManager::Create(
	 pM2D);
	 if (OMS) {
	 OMS->smoothAllVerts(2);
	 delete OMS;
	 }
	 }

	 printf("Cell count: %d of which\n", pM2D->getNumCells());
	 printf("  %6d quads\n", pM2D->getNumQuadCells());
	 printf("  %6d tris\n", pM2D->getNumTriCells());

	 // File output.
	 writeVTK(*pM2D, strOutFileName);
	 writeNative(*pM2D, strOutFileName, false, "");
	 pM2D->evaluateQuality();
	 makeFileName(strQualFileName, "%s.qual", strOutFileName,
	 "main() [edam2d.cxx]");
	 pM2D->writeQualityToFile(strQualFileName);
	 //	meshToBdryFile("terrible.bdry", pM2D);

	 delete pM2D;*/
}
