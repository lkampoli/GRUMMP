#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;
// Generally speaking, the functions in this file are generic topology
// change functions that are members of Face, Cell, etc.  These
// functions are intended only to be called by Mesh member functions
// that are changing topology in ways that only valid entities can ever
// exist (although topology may be incomplete at times).

void Face::assign(Cell* const pCA, Cell* const pCB, Vert* const pV0,
		Vert* const pV1, Vert* const pV2, Vert* const pV3) {
	assert(isValid());
	m_leftCell = pCA;
	m_rightCell = pCB;
	setFaceLoc();
	setVerts(pV0, pV1, pV2, pV3);
}

/// This topology change function is only used in createFace (both in 2D
/// and 3D) and in vAssign.
void Face::setVerts(Vert* pV0, Vert* pV1, Vert* pV2, Vert* pV3) {
//#pragma omp critical(FacesetVerts123)

	assert(isValid());
	m_verts[0] = pV0;
	m_verts[1] = pV1;
	m_verts[2] = pV2;
	m_verts[3] = pV3;
	int nVerts = getNumVerts();
	for (int i = 0; i < nVerts; i++) {
		assert(m_verts[i]->isValid());
		m_verts[i]->addFace(this);
#ifndef OMIT_VERTEX_HINTS
		if (!m_verts[i]->getHintFace()->isValid())
			m_verts[i]->setHintFace(this);
#endif
	}

}

void Mesh::meltIntBdryFaces(std::vector<Face*>& to_lock) {

	GR_index_t i;

	for (i = 0; i < getNumIntBdryFaces(); i++) {
		BFace *pBF = getIntBFace(i);

		Face* pFLeft = pBF->getFace(0);
		Face* pFRight = pBF->getFace(1);

		logMessage(1, "L: %p   R: %p\n", pFLeft, pFRight);
		Cell* pCRight = pFRight->getOppositeCell(pBF);

		pCRight->removeFace(pFRight);
		pFRight->removeCell(pCRight);
		deleteBFace(pBF);
		deleteFace(pFRight);

		pFLeft->addCell(pCRight);
		pCRight->addFace(pFLeft);
		assert(pFLeft->doFullCheck());
		assert(pCRight->doFullCheck());

		to_lock.push_back(pFLeft);

		// cout<<getFaceIndex(pFLeft)<<endl;

	}

}

//Tag faces and their verts as ePseudosurface(pFLeft from above) that we don't want to lose during anisotropic mesh adaptation as a result of mesh operations
//such as swapping, coarsening, etc.Working good now.....
void Mesh::lock_faces(std::vector<Face*>& faces) {

	for (GR_index_t i = 0; i < faces.size(); i++) {

		Face *face = faces[i];

		Vert *v0;
		Vert *v1;

		v0 = face->getVert(0);
		v1 = face->getVert(1);
		assert(v0->isValid());
		assert(v1->isValid());

		//eBdryApex Vert should not be tagged as PseudoSurface
		if (v0->getVertType() != Vert::eBdryApex)
			v0->setType(Vert::ePseudoSurface);
		//cout << getVertIndex(v0) << endl;
		v1->setType(Vert::ePseudoSurface);
		assert(v1->isValid());

		if ((v0->getVertType() == 5 && v1->getVertType() == 5)|| (v0->getVertType() == 1 && v1->getVertType() == 5)) {
			face->setFaceLoc(Face::ePseudoSurface); //prevents faces from swapping and face deletion
			assert(face->isValid());
			face->lock();   //prevents faces from swapping
		}

	//	face->printFaceInfo(); // Print to check if faces and verts are tagged as ePseudosurface

	}
}

