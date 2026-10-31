/*
 * SmoothingManager2D.cxx
 *
 *  Created on: Nov 18, 2011
 *      Author: cfog
 */

#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_GRCurve.h"
#include "GR_Cell.h"
#include "GR_BFace2D.h"
#include "GR_Mesh2D.h"
#include "OptMS.h"
#include "GR_misc.h"

#ifdef HAVE_MESQUITE
#include "Mesquite.hpp"
#include "MsqIMesh.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "PlanarDomain.hpp"
#endif

#include "GR_SmoothingManager.h"

#define iMaxVertNeigh MAX_NUM_PTS

namespace GRUMMP
{
  OptMSSmoothingManager2D::OptMSSmoothingManager2D(Mesh2D * const pMeshToSmooth,
						   int qualMeasure) :
      OptMSSmoothingManager(pMeshToSmooth, 2)
  {
    if (qualMeasure != OPTMS_DEFAULT)
      setSmoothingGoal(qualMeasure);
    m_neighCoords = new double*[iMaxVertNeigh];
    m_faceConn = new int*[iMaxVertNeigh];
    for (int i = 0; i < iMaxVertNeigh; i++) {
      m_neighCoords[i] = new double[2];
      m_faceConn[i] = new int[2];
    }
    setSmoothingThreshold(30.);
    setSmoothingTechnique('f');
  }

  OptMSSmoothingManager2D::~OptMSSmoothingManager2D()
  {
    assert(m_neighCoords && m_faceConn);
    for (int i = 0; i < iMaxVertNeigh; i++) {
      delete[] m_neighCoords[i];
      delete[] m_faceConn[i];
    }
    delete[] m_neighCoords;
    delete[] m_faceConn;
  }

  void
  OptMSSmoothingManager2D::setSmoothingGoal(const int goal)
  {
    SMsetSmoothFunction(m_smoothData, goal);
  }

  int
  OptMSSmoothingManager2D::smoothVert(Vert* vert)
  {
    if (vert->isDeleted() || !vert->isSmoothable()) {
      return 0;
    }
    m_smoothData->local_mesh->Quad = 0;
    double adCurrPt[] =
      { vert->x(), vert->y() };

    bool qBdryVert = false;
    bool qIntBdryVert = false;
    int iNeigh = 0;
    Face *pF = vert->getFace(0);
    Face *pFIBE;
    Cell *pC;
    Cell *pCIBE;
//	int vi = m_pMesh->getVertIndex(vert);
//	printf("vert %d %f %f\n",vi, vert->x(),vert->y());
//	Vert * v2 = m_pMesh->getVert(54);
//	printf("vert %d %f %f\n",vi, v2->x(),v2->y());
//	for(int ff = 0; ff < m_pMesh->getNumFaces(); ff++){
//		Face * face = m_pMesh->getFace(ff);
//		int vi0 = m_pMesh->getVertIndex(face->getVert(0));
//		int vi1 = m_pMesh->getVertIndex(face->getVert(1));
//		printf("face %d %d %d %f %f %f %f\n",ff,vi0,vi1,face->getVert(0)->x(),
//				face->getVert(0)->y(),face->getVert(1)->x(),face->getVert(1)->y());
//
//	}
    if (pF->getVert(0) == vert)
      pC = pF->getLeftCell();
    else
      pC = pF->getRightCell();

    Face *pFStart = pF;
    do {
//		if (pC->getType() == Cell::eQuadCell) return 0;

      if (pC->getType() == Cell::eQuadCell && iNeigh < iMaxVertNeigh) {
	m_smoothData->local_mesh->Quad = 1;
	Vert *vertOpp;
	if (pF->getVert(0) == vert)
	  vertOpp = pF->getVert(1);
	else
	  vertOpp = pF->getVert(0);
	m_neighCoords[iNeigh][0] = vertOpp->x();
	m_neighCoords[iNeigh][1] = vertOpp->y();
	iNeigh++;
      }
      if (pC->getType() == Cell::eTriCell && iNeigh < iMaxVertNeigh) {
	Vert *vertOpp;
	if (pF->getVert(0) == vert)
	  vertOpp = pF->getVert(1);
	else
	  vertOpp = pF->getVert(0);
	m_neighCoords[iNeigh][0] = vertOpp->x();
	m_neighCoords[iNeigh][1] = vertOpp->y();
	iNeigh++;
      }
      else if (pC->getType() == Cell::eBdryEdge) {
	qBdryVert = true;
	break;
      }

      pC = pF->getOppositeCell(pC);
      if (pC->getType() == Cell::eIntBdryEdge) {
	pFIBE = pF;
	pCIBE = pC;
	Face *pFEdge = pC->getFace(0);
	assert(pFEdge->isValid());
	if (pF == pFEdge)
	  pFEdge = pC->getFace(1);
	pF = pFEdge;
	qIntBdryVert = true;
      }
      else if (pC->getType() == Cell::eTriCell) {
	for (int iF = 0; iF < 3; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }
      else if (pC->getType() == Cell::eQuadCell) {
	for (int iF = 0; iF < 4; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }
      else {
	qBdryVert = true;
      }
    }
    while (pF != pFStart && !qBdryVert);

    // Boundary vertex smoothing is handled in a separate routine.
    if (vert->getVertType() == Vert::eBdryApex
	|| vert->getVertType() == Vert::eInteriorFixed) {
      // Can't smooth these
      return 0;
    }
    else if (qBdryVert) {

      if (pC->getType() == Cell::eQuadCell)
	return 0;

      // Note by Charles Boivin - 21/09/99
      // This seems to be reached for vertices on internal boundaries too
      // Made sure the assertion was not triggered for no reason
      assert(
	  (pC->getType() == Cell::eBdryEdge)
	      || (pC->getType() == Cell::eIntBdryEdge));
      if (pC->getType() == Cell::eBdryEdge) {
	BdryEdge *pBE = dynamic_cast<BdryEdge *>(pC);
	return smoothBdryVert(vert, pF, pBE);
      }
      else {
	assert((pCIBE->getType() == Cell::eIntBdryEdge));
	//printf("Smoothing Along Internal Boundary Edge\n");
	IntBdryEdge *pBE = dynamic_cast<IntBdryEdge *>(pCIBE);
	//printf("Int Vert %f %f\n",vert->x(),vert->y());
	return smoothIntBdryVert(vert, pFIBE, pBE);
      }
    }

    if (vert->getVertType() == Vert::eBdryApex|| vert->getVertType() == Vert::eInteriorFixed || vert->getVertType() == Vert::ePseudoSurface) {
      // Can't smooth these
      return 0;
    }

    // Note by Charles Boivin - 21/09/99
    // This seems to be reached for vertices on internal boundaries too
    // Made sure the assertion was not triggered for no reason
    if (qBdryVert) {
      assert((pC->getType() == Cell::eBdryEdge));
      BdryEdge *pBE = dynamic_cast<BdryEdge *>(pC);
      return smoothBdryVert(vert, pF, pBE);
    }
    else if (qIntBdryVert) {
      assert((pCIBE->getType() == Cell::eIntBdryEdge));
      IntBdryEdge *pBE = dynamic_cast<IntBdryEdge *>(pCIBE);
      return smoothIntBdryVert(vert, pFIBE, pBE);
    }

    // List all opposite faces (w/ correct orientation) using indices
    // into the location array
//	if (pC->getType() == Cell::eQuadCell){
//		for (int ii = 1; ii < iQuads; ii++){
//			m_faceConn[ii][0] = (ii)*2;
//			m_faceConn[ii][1] = 2*ii + 1;
//			m_faceConn[ii][2] = 2*ii - 2;
//		}
//		m_faceConn[0][0] = 0;
//		m_faceConn[0][1] = 1;
//		m_faceConn[0][2] = iNeigh-2;
//	}
//	else
//	{
    m_faceConn[0][0] = 0;
    m_faceConn[0][1] = iNeigh - 1;
    int ii;
    for (ii = 1; ii < iNeigh; ii++) {
      m_faceConn[ii][0] = ii;
      m_faceConn[ii][1] = ii - 1;
    }

    if (iNeigh <= iMaxVertNeigh) {
      double adLoc[] =
	{ adCurrPt[0], adCurrPt[1] };
      SMsmooth(iNeigh, iNeigh, adCurrPt, m_neighCoords, m_faceConn,
	       m_smoothData, 0);

      if (dDIST2D (adLoc, adCurrPt) > 1.e-8) {
	vert->setCoords(2, adCurrPt);
	m_pMesh->moveVertEvent(vert);
#ifndef NDEBUG
	m_pMesh->writeTempMesh();
#endif
	return 1;

      }
    }

    return 0;
  }

  int
  OptMSSmoothingManager2D::smoothBdryVert(Vert * const vert,
					  Face * const faceOnBdry,
					  BdryEdge * const bdryEdgeInit)
  {

//	assert (m_pMesh->isSimplicial ());
    assert(bdryEdgeInit->hasFace(faceOnBdry));
    assert(faceOnBdry->hasVert(vert));

    if (!m_pMesh->areBdryChangesAllowed()
	|| vert->getVertType() == Vert::eBdryApex || vert->isSmallAngleVert()
	|| vert->isShellVert()) {
      return 0;

    }

    // This code shouldn't be here if theres an internal boundary.
    if (bdryEdgeInit->getType() == Cell::eIntBdryEdge)
      return 0;

    // Also misbehaves for periodic boundaries, for the same reasons as
    // for internal boundaries.

    //   if (bdryEdgeInit->pPatchPointer()->qIsPeriodic())
    //     return false;

    int iNeigh = 0;

    Cell *pC = bdryEdgeInit;
    Face *pF = faceOnBdry;
    do {
      if (iNeigh < iMaxVertNeigh) {
	Vert *vertOpp;
	if (pF->getVert(0) == vert)
	  vertOpp = pF->getVert(1);
	else
	  vertOpp = pF->getVert(0);
	m_neighCoords[iNeigh][0] = vertOpp->x();
	m_neighCoords[iNeigh][1] = vertOpp->y();
	iNeigh++;
      }

      pC = pF->getOppositeCell(pC);
      if (pC->getType() == Cell::eIntBdryEdge) {
	iNeigh--;
	Face *pFEdge = pC->getFace(0);
	assert(pFEdge->isValid());
	if (pF == pFEdge)
	  pFEdge = pC->getFace(1);
	pF = pFEdge;
      }
      if (pC->getType() == Cell::eTriCell) {
	for (int iF = 0; iF < 3; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }
      if (pC->getType() == Cell::eQuadCell)
	return 0;
      if (pC->getType() == Cell::eQuadCell) {
	for (int iF = 0; iF < 4; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }
    }
    while (pC->getType() == Cell::eQuadCell || pC->getType() == Cell::eTriCell
	|| pC->getType() == Cell::eIntBdryEdge);

    // Modified 21/09/99 by Charles Boivin
    // Added the possibility of an internal boundary.
    assert(
	(pC->getType() == Cell::eBdryEdge)
	    || (pC->getType() == Cell::eIntBdryEdge));
    BdryEdgeBase *bdryEdgeFinal = dynamic_cast<BdryEdgeBase *>(pC);

    // Check whether it's okay to smooth this vertex.
    //   Are the boundary conditions the same?

    int bc1 = bdryEdgeInit->getBdryCondition(true), bc2 =
	bdryEdgeInit->getBdryCondition(false);
    int bc3 = bdryEdgeFinal->getBdryCondition(true), bc4 =
	bdryEdgeFinal->getBdryCondition(false);

    if (bc1 != bc3 && bc1 != bc4 && bc2 != bc3 && bc2 != bc4) {
      return 0;
    }
    //   Are the faces co-linear?  If so, the unit normals will have a
    //   dot product which has magnitude of almost precisely 1.
    {
      double adNorm_Init[3], adNorm_Final[3];
      faceOnBdry->calcUnitNormal(adNorm_Init);
      pF->calcUnitNormal(adNorm_Final);
      double dDot = fabs(dDOT2D(adNorm_Init, adNorm_Final));
      if (iFuzzyComp(dDot, 1.) != 0)
	return 0;
      // Co-linear!  Set the surface normal accordingly.  Orientation of
      // the normal (inward vs. outward) doesn't matter here.
      SMsetNormal(adNorm_Init);
    }

    // Now that we know we're going to smooth, set up the connectivity,
    // etc and get on with it.

    // List all opposite faces (w/ correct orientation) using indices
    // into the location array. There are iNeigh-1 such faces.
    if ((faceOnBdry->getRightCell() == static_cast<Cell *>(bdryEdgeInit)
	&& faceOnBdry->getVert(1) == vert)
	|| (faceOnBdry->getLeftCell() == static_cast<Cell *>(bdryEdgeInit)
	    && faceOnBdry->getVert(0) == vert)) {
      for (int ii = 0; ii < iNeigh - 1; ii++) {
	m_faceConn[ii][0] = ii + 1;
	m_faceConn[ii][1] = ii;
      }
    }
    else {
      for (int ii = 0; ii < iNeigh - 1; ii++) {
	m_faceConn[ii][0] = ii;
	m_faceConn[ii][1] = ii + 1;
      }
    }

    if (iNeigh <= iMaxVertNeigh) {
      double adCurrPt[] =
	{ vert->x(), vert->y() };
      double adLoc[] =
	{ adCurrPt[0], adCurrPt[1] };
      SMsmooth(iNeigh, iNeigh - 1, adCurrPt, m_neighCoords, m_faceConn,
	       m_smoothData, 1);

      if (dDIST2D (adLoc, adCurrPt) > 1.e-8) {
	// If we're using legit curves (class GRCurve) as the underlying
	// geometry, then we need to update the vertex parameter when
	// we're done here.
	GRCurve *initCurve = bdryEdgeInit->getCurve();
	GRCurve *finalCurve = bdryEdgeFinal->getCurve();
	assert(initCurve == finalCurve);
	if (initCurve != NULL) {
	  // Which end of the edge is the current vert for both edges?
	  bool isInitVert0 = (bdryEdgeInit->getVert(0) == vert);
	  bool isFinalVert0 = (bdryEdgeFinal->getVert(0) == vert);
	  double midParam = -1;
	  CubitVector loc(adCurrPt[0], adCurrPt[1], 0), onCurveLoc, tangent,
	      curvature;
	  initCurve->closest_point(loc, onCurveLoc, &tangent, &curvature,
				   &midParam);
	  if (isInitVert0)
	    bdryEdgeInit->setVert0Param(midParam);
	  else
	    bdryEdgeInit->setVert1Param(midParam);
	  if (isFinalVert0)
	    bdryEdgeFinal->setVert0Param(midParam);
	  else
	    bdryEdgeFinal->setVert1Param(midParam);
	}
	vert->setCoords(2, adCurrPt);

	m_pMesh->moveVertEvent(vert);
#ifndef NDEBUG
	m_pMesh->writeTempMesh();
#endif
	return 1;
      }
    }

    return 0;
  }
#ifdef HAVE_MESQUITE
  using namespace ITAPS_GRUMMP;
  using namespace Mesquite;

  SmoothingManagerMesquite2D::SmoothingManagerMesquite2D(Mesh2D* const pM2D) :
      SmoothingManager(pM2D)
  {
    Mesquite::MsqError MesqErr;
    //getting msqMesh (mesh from Mesquite) from pM2D
    iMesh_Plane *iMP = new iMesh_Plane(pM2D);
    assert(iMP != NULL);
    iMesh_Instance instance = new iMesh_Instance_Private(iMP);
    assert(instance);
    //   Create a Mesquite::MsqIMesh from iMP

    int iMeshErr;
    iBase_EntitySetHandle root_set;
    iMesh_getRootSet(instance, &root_set, &iMeshErr);

    iBase_TagHandle fixed_tag;
    iMesh_createTag(instance, "fixed", 1, iBase_INTEGER, &fixed_tag, &iMeshErr,
		    5);
    if (iMeshErr != iBase_SUCCESS) {
      vFatalError("failed to create a tag",
		  "Mesquite smoothing initialization");
    }

    // Now set data for it, by looping over verts and transcribing.
    for (unsigned int iV = 0; iV < instance->pGB->getNumVertices(); iV++) {

      Vert *pV = instance->pGB->pMesh()->getVert(iV);
      int iType = pV->getVertType();
      if (iType == Vert::eBdryApex || iType == Vert::eBdryCurve
	  || iType == Vert::eBdryTwoSide || iType == Vert::eBdry || iType == Vert::ePseudoSurface) //ePseudoSurface Vert needed here,else code stop at ePseudoSurface Vert (Shahzaib)
	iMesh_setIntData(instance, pV, fixed_tag, 1, &iMeshErr);
      else
	iMesh_setIntData(instance, pV, fixed_tag, 0, &iMeshErr);

      assert(iMeshErr == iBase_SUCCESS);
    }

    msqMesh = new Mesquite::MsqIMesh(instance, root_set, iBase_FACE, MesqErr,
				     &fixed_tag);

    if (MSQ_CHKERR(MesqErr)) {
      if (msqMesh)
	delete msqMesh;
      std::cout << MesqErr << std::endl;
      vFatalError("failed to create a Mesquite mesh",
		  "Mesquite smoothing initialization");
    }
  }

  SmoothingManagerMesquite2D::~SmoothingManagerMesquite2D()
  {
    iMesh_Instance instance = msqMesh->get_imesh_instance();
    delete msqMesh;
    delete instance;
  }
  
  int
  SmoothingManagerMesquite2D::smoothAllVerts(const int)

  {
    Mesquite::MsqError MesqErr;
    //Getting the normal plane and then smoothing
    Mesquite::Vector3D normal(0, 0, 1);
    Mesquite::Vector3D point(0, 0, 0);
    Mesquite::PlanarDomain my_mesh_plane(normal, point);
    Queue.run_instructions(msqMesh, &my_mesh_plane, MesqErr);
    if (MesqErr)
      std::cout << MesqErr << std::endl;
    return 1;			// Should be actual number of vertices.
  }

  void
  SmoothingManagerMesquite2D::vSetTargetfromMetric()
  {
    Mesquite::MsqError MesqErr;
    //calculating the target from metric. lid_tag1, lid_tag2, lid_tag3, lid_tag4 are
    //the components of the target matrix.

    std::vector<Mesquite::Mesh::ElementHandle> elements;
    msqMesh->get_all_elements(elements, MesqErr);
    std::vector<double> lid1((m_pMesh->getNumCells()));
    std::vector<double> lid2((m_pMesh->getNumCells()));
    std::vector<double> lid3((m_pMesh->getNumCells()));
    std::vector<double> lid4((m_pMesh->getNumCells()));
    double adMetric[3];
    double adEigVals[2];
    double adTest[4], adTestOri[4];
    double adEigVec1[2], adEigVec2[2];
    for (unsigned int i = 0; i < m_pMesh->getNumCells(); i++) {
      Vert *pV0 = m_pMesh->getCell(i)->getVert(0);
      Vert *pV1 = m_pMesh->getCell(i)->getVert(1);
      Vert *pV2 = m_pMesh->getCell(i)->getVert(2);

      ////////////////////////////////////////////////
       ////////////////////////////////////////////////
//       const double adLoc1[2]={pV0->x(),pV0->y()};
//       const double adLoc2[2]={pV1->x(),pV1->y()};
//       const double adLoc3[2]={pV2->x(),pV2->y()};
//       const double adNewPoint[2]={pV2->x(),pV2->y()};
//       double adVal1[3],adVal2[3],adVal3[3];
//
//       adVal1[0]=pV0->getMetric(0);
//       adVal1[1]=pV0->getMetric(1);
//       adVal1[2]=pV0->getMetric(2);
//       adVal2[0]=pV1->getMetric(0);
//       adVal2[1]=pV1->getMetric(1);
//       adVal2[2]=pV1->getMetric(2);
//       adVal3[0]=pV2->getMetric(0);
//       adVal3[1]=pV2->getMetric(1);
//       adVal3[2]=pV2->getMetric(2);
//       double dDetM;
//       double dX[3]={adLoc1[0],adLoc2[0],adLoc3[0]};
//       double dY[3]={adLoc1[1],adLoc2[1],adLoc3[1]};
//       double adNewLoc[3]={1,adNewPoint[0],adNewPoint[1]};
//       double adSol[3];
//       int i,j;
//       double a2dMi[3][3];
//
//       dDetM=dX[1]*dY[2]-dY[1]*dX[2]-dX[0]*dY[2]+dY[0]*dX[2]+dX[0]*dY[1]-dY[0]*dX[1];
//       if(fabs(dDetM)<1e-10){
//       //Just use the average
//       for(i=0;i<3;i++){
//       adMetric[i]=1.0/3.0*(adVal1[i]+adVal2[i]+adVal3[i]);
//       }
//
//       }
//
//       else{
//       a2dMi[0][0]=dX[1]*dY[2]-dY[1]*dX[2];
//       a2dMi[0][1]=dX[2]*dY[0]-dX[0]*dY[2];
//       a2dMi[0][2]=dX[0]*dY[1]-dY[0]*dX[1];
//       a2dMi[1][0]=dY[1]-dY[2];
//       a2dMi[1][1]=dY[2]-dY[0];
//       a2dMi[1][2]=dY[0]-dY[1];
//       a2dMi[2][0]=dX[2]-dX[1];
//       a2dMi[2][1]=dX[0]-dX[2];
//       a2dMi[2][2]=dX[1]-dX[0];
//       for(j=0;j<3;j++){
//       for(i=0;i<3;i++){
//       a2dMi[i][j]/=dDetM;
//       }
//       }
//
//       for(j=0;j<3;j++){
//       adSol[j]=0;
//       for(i=0;i<3;i++){
//       adSol[j]+=adNewLoc[i]*a2dMi[i][j];
//       }
//       }
//       for(j=0;j<3;j++){
//       adMetric[j]=adSol[0]*adVal1[j];
//       adMetric[j]+=adSol[1]*adVal2[j];
//       adMetric[j]+=adSol[2]*adVal3[j];
//       }
//       }
//       	  ///////////////////////////////////////////////
		  ////////////////////////////////////////////////
		  //vLinearInterpolate(adLoc0,adVal0,adLoc1,adVal1,adLoc2,adVal2,const double adNewPoint[2],adNMetric,3);
      adMetric[0] = (pV0->getMetric(0) + pV1->getMetric(0) + pV2->getMetric(0))
	  / 3.;
      adMetric[1] = (pV0->getMetric(1) + pV1->getMetric(1) + pV2->getMetric(1))
	  / 3.;
      adMetric[2] = (pV0->getMetric(2) + pV1->getMetric(2) + pV2->getMetric(2))
	  / 3.;
      //Computing the eigenvalues of the Metric
      double dRoot = sqrt(
	  (adMetric[0] - adMetric[2]) * (adMetric[0] - adMetric[2])
	      + 4 * adMetric[1] * adMetric[1]);
      adEigVals[0] = 0.5 * ((adMetric[0] + adMetric[2]) + dRoot);//Larger eigenvalue
      adEigVals[1] = 0.5 * ((adMetric[0] + adMetric[2]) - dRoot);//smaller eigenvalue

      //Computing the normalized eigenvectors of Metric

      if (adMetric[1] < 1e-6 && adMetric[1] > -1e-6) {

	adEigVec1[0] = 1;
	adEigVec1[1] = 0;
	adEigVec2[0] = 0;
	adEigVec2[1] = 1;
	double adRes[2], adDiff[2];
	adRes[0] = adMetric[0] * adEigVec1[0] + adMetric[1] * adEigVec1[1];

	adRes[1] = adMetric[1] * adEigVec1[0] + adMetric[2] * adEigVec1[1];

	adDiff[0] = fabs(adRes[0] - adEigVals[0] * adEigVec1[0]);
	adDiff[1] = fabs(adRes[1] - adEigVals[0] * adEigVec1[1]);
	if (adDiff[0] > 1e-6 || adDiff[1] > 1e-6) {
	  adEigVec1[0] = 0;
	  adEigVec1[1] = 1;
	  adEigVec2[0] = 1;
	  adEigVec2[1] = 0;
	}

      }
      else {
	// EigenVector corresponding to larger eigenvalue
	adEigVec1[0] = adMetric[0] - adEigVals[0];
	adEigVec1[1] = -adMetric[1];
	double dNorm = sqrt(
	    adEigVec1[0] * adEigVec1[0] + adEigVec1[1] * adEigVec1[1]);
	adEigVec1[0] /= dNorm;
	adEigVec1[1] /= dNorm;

	//EigenVector corresponding to smaller eigenvalue
	adEigVec2[0] = adMetric[0] - adEigVals[1];
	adEigVec2[1] = -adMetric[1];
	dNorm = sqrt(adEigVec2[0] * adEigVec2[0] + adEigVec2[1] * adEigVec2[1]);
	adEigVec2[0] /= dNorm;
	adEigVec2[1] /= dNorm;
      }

      //setting the Target based on the eigenvectors and eigenvalues

//	  		        adTest[0]=adEigVec1[0]*sqrt(adEigVals[1]);
//	  		        adTest[2]=adEigVec1[1]*sqrt(adEigVals[1]);
//
//	  		        adTest[1]=adEigVec2[0]*sqrt(adEigVals[0]);
//	  		        adTest[3]=adEigVec2[1]*sqrt(adEigVals[0]);

      adTest[0] = adEigVec1[0] * (1. / (sqrt(adEigVals[0])));
      adTest[2] = adEigVec1[1] * (1. / (sqrt(adEigVals[0])));

      adTest[1] = adEigVec2[0] * (1. / (sqrt(adEigVals[1])));
      adTest[3] = adEigVec2[1] * (1. / (sqrt(adEigVals[1])));

      //testing the orientation of cell
      if (adTest[0] * adTest[3] > adTest[1] * adTest[2]) {
	adTestOri[0] = adTest[0];
	adTestOri[1] = adTest[1];
	adTestOri[2] = adTest[2];
	adTestOri[3] = adTest[3];
      }
      else {
	adTestOri[0] = adTest[1];
	adTestOri[1] = adTest[0];
	adTestOri[2] = adTest[3];
	adTestOri[3] = adTest[2];
      }

      //testing the norm of the first column
      if (adTestOri[1] * adTestOri[1] + adTestOri[3] * adTestOri[3]
	  > adTestOri[0] * adTestOri[0] + adTestOri[2] * adTestOri[2]) {
	lid1[i] = -adTestOri[1];
	lid2[i] = adTestOri[0];
	lid3[i] = -adTestOri[3];
	lid4[i] = adTestOri[2];
      }
      else {
	lid1[i] = adTestOri[0];
	lid2[i] = adTestOri[1];
	lid3[i] = adTestOri[2];
	lid4[i] = adTestOri[3];

      }
      if (adMetric[1] < 1e-6 && adMetric[1] > -1e-6) { //zero anti-diagonal entries (approximately)
	if (fabs(adMetric[0] - adMetric[2]) < 1e-6) {
	  lid1[i] = 1;
	  lid2[i] = 0.5;
	  lid3[i] = 0;
	  lid4[i] = sqrt(3) / 2.;
	}
      }
    }
    const char LOCAL_ID_NAME1[] = "LOCAL_ID1";
    const char LOCAL_ID_NAME2[] = "LOCAL_ID2";
    const char LOCAL_ID_NAME3[] = "LOCAL_ID3";
    const char LOCAL_ID_NAME4[] = "LOCAL_ID4";
    Mesquite::TagHandle lid_tag1 = msqMesh->tag_create(LOCAL_ID_NAME1,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag2 = msqMesh->tag_create(LOCAL_ID_NAME2,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag3 = msqMesh->tag_create(LOCAL_ID_NAME3,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag4 = msqMesh->tag_create(LOCAL_ID_NAME4,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    size_t num_cells = (m_pMesh->getNumCells());
    msqMesh->tag_set_vertex_data(lid_tag1, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid1), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag2, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid2), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag3, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid3), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag4, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid4), MesqErr);

  }

  /*
   * For the local smoothing function of SmoothingManagerMesquite2D, the maxTime at the AnisWrapper.cxx for 2D smoothing must be set really low (1 or 2)
   */
  int
  SmoothingManagerMesquite2D::smoothVert(Vert* pVtoSmooth)
  {

//	  Smooths only the desired point by using Mesquite in a small patch from the neighbouring cells of the point.
//	  Only smooths interior points
    if (pVtoSmooth->getVertType() != Vert::eInterior || pVtoSmooth->getVertType() == Vert::ePseudoSurface)
      return 0;

//	  creates a small patch from the neigbouring cells
    Mesh2D M2DPatch;
    std::set<Cell*> spC;
    std::set<Vert*> spV;
    findNeighborhoodInfo(pVtoSmooth, spC, spV);

    std::set<Vert*>::iterator iterTest, iterTestE;

    iterTest = spV.begin();
    iterTestE = spV.end();

    std::set<Cell*>::iterator iterC, iterCE;
    iterC = spC.begin();
    iterCE = spC.end();

    int iVSize;
    int iFSize;
    int iCSize;

    iVSize = spV.size();
    iFSize = iVSize * 2;
    iCSize = spC.size();

//	  First we create arrays for cells, faces, and vertices. This arrays are set in a counterclockwise order
    Cell *pC[iCSize + 1];
    Vert *pV[iVSize];
    Face *pF[iFSize];

    pC[0] = (*iterC);

    assert(pC[0]->doFullCheck());
    std::set<Vert*> spVOrder;
    pF[1] = pC[0]->getOppositeFace(pVtoSmooth);

    if (checkOrient2D(pVtoSmooth, pF[1]->getVert(0), pF[1]->getVert(1)) == 1) {
      spVOrder.insert(pF[1]->getVert(0));
      spVOrder.insert(pF[1]->getVert(1));
      pV[0] = pF[1]->getVert(0);
      pV[1] = pF[1]->getVert(1);
      pF[0] = findCommonFace(pV[0], pVtoSmooth, true);
    }
    else {
      spVOrder.insert(pF[1]->getVert(1));
      spVOrder.insert(pF[1]->getVert(0));
      pV[0] = pF[1]->getVert(1);
      pV[1] = pF[1]->getVert(0);
      pF[0] = findCommonFace(pV[1], pVtoSmooth, true);
    }

    pF[2] = findCommonFace(pVtoSmooth, pV[1], true);
    pC[1] = pF[2]->getOppositeCell(pC[0]);
    pV[2] = pC[1]->getOppositeVert(pF[2]);
    spVOrder.insert(pV[2]);
    pF[3] = findCommonFace(pV[1], pV[2], true);

    for (unsigned int j = 2; j < spV.size(); j++) {
      pF[j * 2] = findCommonFace(pV[j], pVtoSmooth, true);
      pC[j] = pF[j * 2]->getOppositeCell(pC[j - 1]);
      if (pC[j] == pC[0])
	break;
      pV[j + 1] = pC[j]->getOppositeVert(pF[j * 2]);
      spVOrder.insert(pV[j + 1]);
      pF[(j * 2) + 1] = findCommonFace(pV[j], pV[j + 1], true);
    }

    std::set<Vert*>::iterator iter, iterE;
    iter = spV.begin();
    iterE = spV.end();
    Vert *pVNeigh[spV.size()];

    int iV = 0;
    double adM[3];

//	  Creating a patch, cells, faces and vertices are ordered in a counterclockwise configuration
    Vert *pVSmoothNew = M2DPatch.createVert(pVtoSmooth->x(), pVtoSmooth->y());

    adM[0] = pVtoSmooth->getMetric(0);
    adM[1] = pVtoSmooth->getMetric(1);
    adM[2] = pVtoSmooth->getMetric(2);

    pVSmoothNew->setMetric(adM);
    pVSmoothNew->setType(Vert::eInterior);

    for (int j = 0; j < iVSize; j++) {
      Vert *pViter = pV[j];
      pVNeigh[iV] = M2DPatch.createVert(pViter->x(), pViter->y());

      adM[0] = pViter->getMetric(0);
      adM[1] = pViter->getMetric(1);
      adM[2] = pViter->getMetric(2);

      pVNeigh[iV]->setMetric(adM);
      pVNeigh[iV]->setType(Vert::eBdry);
      iV++;
    }

    Face *pFNew[spVOrder.size() * 2];

    int iNF = 0;
    bool bAlreadyExisted = false;
    for (unsigned int j = 0; j < spV.size(); j++) {
      pFNew[iNF] = M2DPatch.createFace(bAlreadyExisted, pVSmoothNew,
				       pVNeigh[j]);
      assert(pFNew[iNF]->doFullCheck());
      iNF++;
      int k = j + 1;
      if (k == iVSize)
	k = 0;
      pFNew[iNF] = M2DPatch.createFace(bAlreadyExisted, pVNeigh[j], pVNeigh[k]);
      assert(pFNew[iNF]->doFullCheck());
      iNF++;
    }

    for (unsigned int j = 0; j < spV.size(); j++) {
      Cell *pCCheck;
      int k = j + 1;
      if (k == iVSize)
	k = 0;
      pCCheck = M2DPatch.createTriCell(pVSmoothNew, pVNeigh[j], pVNeigh[k], 1,
				       true);
      assert(pCCheck->doFullCheck());
    }

//	  setting boundaries for the patch
    for (unsigned int j = 1; j < spV.size() * 2; j = j + 2) {
      int iTestNF = pFNew[j]->getNumCells();
      assert(iTestNF == 1);
      M2DPatch.createBFace(pFNew[j]);
    }

    assert(M2DPatch.isValid());

//	  calling mesquite to smooth the patch, then setting new coordinates in the point to be smoothed
//	  GRUMMP::SmoothingManagerMesquite2D *pSMM2DPatch = new GRUMMP::SmoothingManagerMesquite2D(&M2DPatch);
//	  pSMM2DPatch-> vSetTargetfromMetric();
//	  pSMM2DPatch->smoothAllVerts(1);

    double dNewCoords[2];

    dNewCoords[0] = pVSmoothNew->x();
    dNewCoords[1] = pVSmoothNew->y();

    pVtoSmooth->setCoords(2, dNewCoords);

//	  delete pSMM2DPatch;

    return 1;
  }

#endif

  int
  OptMSSmoothingManager2D::smoothIntBdryVert(Vert * const vert,
					     Face * const faceOnBdry,
					     IntBdryEdge * const bdryEdgeInit)
  {

    assert(bdryEdgeInit->hasFace(faceOnBdry));
    assert(faceOnBdry->hasVert(vert));

    if (!m_pMesh->areBdryChangesAllowed()|| vert->getVertType() == Vert::eBdryApex || vert->isSmallAngleVert()|| vert->isShellVert()) {
      return 0;
    }

    int iNeigh = 0;

    Cell *pC = bdryEdgeInit;
    Face *pF = faceOnBdry;
    Face *pF2 = pF;
    Vert *vertOpp = vert;
    Vert *vertFirst = vert;

    while (vertOpp != vertFirst || (vertOpp == vertFirst && iNeigh < 2)) {
      if (pC->getType() == Cell::eQuadCell)
	return 0;
      if (iNeigh < iMaxVertNeigh) {

	if (pF->getVert(0) == vert)
	  vertOpp = pF->getVert(1);
	else
	  vertOpp = pF->getVert(0);
	if (iNeigh == 0) {
	  vertFirst = vertOpp;
	}

	if ((vertOpp != vertFirst) || (vertOpp == vertFirst && iNeigh == 0)) {
	  m_neighCoords[iNeigh][0] = vertOpp->x();
	  m_neighCoords[iNeigh][1] = vertOpp->y();
	  iNeigh++;
	}
      }

      pC = pF->getOppositeCell(pC);
      if (pC->getType() == Cell::eIntBdryEdge) {
	iNeigh--;
//			for (int iF = 0; iF < 2; iF++) {
//				Face *pFEdge = pC->getFace(iF);
//				assert (pFEdge->isValid ());
//				if (!((pFEdge->getVert(0) == faceOnBdry->getVert(0) && pFEdge->getVert(1) == faceOnBdry->getVert(1) ) ||
//						(pFEdge->getVert(0) == faceOnBdry->getVert(1) && pFEdge->getVert(1) == faceOnBdry->getVert(0) )))
//					pF2 = pFEdge;
//				if (pF == pFEdge) continue;
//				if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
//					pF = pFEdge;
//					break;
//				}
//			}
	Face *pFEdge = pC->getFace(0);
	assert(pFEdge->isValid());
	if (pF == pFEdge)
	  pFEdge = pC->getFace(1);
	pF = pFEdge;
      }
      if (pC->getType() == Cell::eTriCell) {
	for (int iF = 0; iF < 3; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }
      if (pC->getType() == Cell::eQuadCell) {
	for (int iF = 0; iF < 4; iF++) {
	  Face *pFEdge = pC->getFace(iF);
	  assert(pFEdge->isValid());
	  if (pF == pFEdge)
	    continue;
	  if (pFEdge->getVert(0) == vert || pFEdge->getVert(1) == vert) {
	    pF = pFEdge;
	    break;
	  }
	}
      }

    }
    iNeigh++;
//	printf("vertex to smooth: %f %f %d\n",vert->x(),vert->y(),vert->getNumFaces());
//
//	for(int i = 0; i < iNeigh; i++){
//		printf("Vert %f %f\n",m_neighCoords[i][0] ,m_neighCoords[i][1] );
//	}
    // Modified 21/09/99 by Charles Boivin
    // Added the possibility of an internal boundary.
    assert((pC->getType() == Cell::eIntBdryEdge));
    assert(pF2->isValid());
    BdryEdgeBase *pBE_Final = dynamic_cast<IntBdryEdge *>(pC);

    //need to make sure we have both boundary internal boundary faces, or at least faces corresponding to them

    double adNorm_Init[3], adNorm_Final[3];
    faceOnBdry->calcUnitNormal(adNorm_Init);
    pF2->calcUnitNormal(adNorm_Final);
    double dDot = fabs(dDOT2D(adNorm_Init, adNorm_Final));

    GRCurve * curve1 = bdryEdgeInit->getCurve();
    GRCurve * curve2 = pBE_Final->getCurve();

    if (curve1 != curve2)
      return 0;

    if (dDot < 0.999)
      return 0;
    // Co-linear!  Set the surface normal accordingly.  Orientation of
    // the normal (inward vs. outward) doesn't matter here.
    CubitVector vertcv, normal;
    vertcv.set(vert->getCoords());
    (curve1->get_curve_geom())->closest_coord_on_curve(vertcv, vertcv);
    curve1->get_curve_geom()->unit_normal(
	curve1->get_curve_geom()->param_at_coord(vertcv), normal);
    normal.get_xyz(adNorm_Final);
    SMsetNormal(adNorm_Final);

    // Now that we know we're going to smooth, set up the connectivity,
    // etc and get on with it.

    m_faceConn[0][1] = 0;
    m_faceConn[0][0] = iNeigh - 1;
    int ii;
    for (ii = 1; ii < iNeigh; ii++) {
      m_faceConn[ii][1] = ii;
      m_faceConn[ii][0] = ii - 1;
    }

    if (iNeigh <= iMaxVertNeigh) {
      double adCurrPt[] =
	{ vert->x(), vert->y() };
      double adLoc[] =
	{ adCurrPt[0], adCurrPt[1] };
      SMsmooth(iNeigh, iNeigh, adCurrPt, m_neighCoords, m_faceConn,
	       m_smoothData, 0);

      // project back onto the curve.
      GRCurveGeom * cg = curve1->get_curve_geom();
      CubitVector coord, closest;
      coord.set(adCurrPt[0], adCurrPt[1], 0.);
      //double adPt[] = {adCurrPt[0], adCurrPt[1]};
      cg->closest_coord_on_curve(coord, closest);
      adCurrPt[0] = closest.x();
      adCurrPt[1] = closest.y();
      // Put this in to make sure distances arent too big
//		if(dDIST2D(adPt,adCurrPt) < dDIST2D(adLoc,adPt))
//			return 0;
      //TODO put in check for whether the move is safe, now that we've
      // switched to closest point on curve - ZAIDE
      if (dDIST2D (adLoc, adCurrPt) > 1.e-8) {
	vert->setCoords(2, adCurrPt);
	m_pMesh->moveVertEvent(vert);
#ifndef NDEBUG
	m_pMesh->writeTempMesh();
#endif
	return 1;
      }
    }
    return 0;
  }

}
