/*
 * SmoothingManager3D.cxx
 *
 *  Created on: Nov 18, 2011
 *      Author: cfog
 */

#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_BFace.h"
#include "GR_VolMesh.h"
#include "OptMS.h"
#include "FacetSurface.hpp"
#include "GR_SmoothingManager.h"
#include "GR_QualMeasure.h"

#ifdef HAVE_MESQUITE
#include "Mesquite.hpp"
#include "MsqIMesh.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "PlanarDomain.hpp"
#include "MeshInterface.hpp"
#endif

#define iMaxVertNeigh MAX_NUM_PTS

namespace GRUMMP
{
  OptMSSmoothingManager3D::OptMSSmoothingManager3D(
      VolMesh * const pMeshToSmooth, int qualMeasure) :
      OptMSSmoothingManager(pMeshToSmooth, 3), m_incCells(), m_neighVerts(), m_incBdryFaces(), m_nearFaces()
  {
    if (qualMeasure != OPTMS_DEFAULT)
      setSmoothingGoal(qualMeasure);
    m_neighCoords = new double*[iMaxVertNeigh];
    m_faceConn = new int*[2 * iMaxVertNeigh];
    for (int i = 0; i < iMaxVertNeigh; i++) {
      m_neighCoords[i] = new double[3];
      m_faceConn[i] = new int[3];
      m_faceConn[i + iMaxVertNeigh] = new int[3];
    }
    setSmoothingThreshold(25.);
    setSmoothingTechnique('f');
  }

  OptMSSmoothingManager3D::~OptMSSmoothingManager3D()
  {
    assert(m_neighCoords && m_faceConn);
    for (int i = 0; i < iMaxVertNeigh; i++) {
      delete[] m_neighCoords[i];
      delete[] m_faceConn[i];
      delete[] m_faceConn[i + iMaxVertNeigh];
    }
    delete[] m_neighCoords;
    delete[] m_faceConn;
  }

  void
  OptMSSmoothingManager3D::setSmoothingGoal(const int goal)
  {
    SMsetSmoothFunction(m_smoothData, goal + 20);
  }
/// attempted switch to handle curved surfaces done by checking for a surface,
/// and then using the normal from the surface
/// and using the original smoothing, and then projecting
/// I make an assumption that if you are eBdryTwoSide, that
/// you are on one surface and one surface only. Otherwise, the vertex
/// would be qBdryApex, or qBdryCurve
  int
  OptMSSmoothingManager3D::smoothVert(Vert* vert)
  {

    bool qBdryVert, qOnSurf = false;
    // Set up coordinates of the current vertex
    double adCurrPt[3];
    // May not use this;
    FacetSurface * facetSurf = NULL;
    adCurrPt[0] = vert->x();
    adCurrPt[1] = vert->y();
    adCurrPt[2] = vert->z();

    if (vert->isDeleted() || !vert->isSmoothable()) {
      return 0;
    }

    // Set up coordinates of the current vertex
    adCurrPt[0] = vert->x();
    adCurrPt[1] = vert->y();
    adCurrPt[2] = vert->z();

    // Reset lists of nearby entities
    m_incCells.clear();
    m_neighVerts.clear();
    m_incBdryFaces.clear();
    m_nearFaces.clear();

    findNeighborhoodInfo(vert, m_incCells, m_neighVerts, &m_incBdryFaces,
			 &qBdryVert, &m_nearFaces);
    // Verify that the right number of each kind of entity exists.
//    assert(!m_pMesh->isSimplicial() ||
//        (((!qBdryVert || vert->getVertType() == Vert::eBdryTwoSide )&&
//                m_incCells.size() == m_neighVerts.size()*2 - 4)
//            ||
//            (qBdryVert &&
//                m_incCells.size() == (m_neighVerts.size()*2 -
//                    m_incBdryFaces.size() - 2))));

    if (qBdryVert && vert->getVertType() != Vert::eBdryTwoSide) {
      // Can only smooth a bdry vertex if it lies in a coplanar patch
      // of the surface and if all its surrounding bdry faces have the
      // same boundary condition.  If these conditions are met, then
      // we need to set up the surface normal prior to smoothing.
      bool qSurfSmoothOK = true;
      std::set<BFace*>::iterator iterBF = m_incBdryFaces.begin();
      TriBFaceBase* pTBF = dynamic_cast<TriBFaceBase*>(*iterBF);
      if (!pTBF->isValid()) {
	// This is a quad bface instead
	return 0;
      }
      int iBC = pTBF->getBdryCondition();
      double adNorm[3];
      pTBF->getFace()->calcUnitNormal(adNorm);
      BasicTopologyEntity * BTE = pTBF->getTopologicalParent();
      if (BTE)
	qOnSurf = true;
      for (iterBF++; iterBF != m_incBdryFaces.end() && qSurfSmoothOK;
	  iterBF++) {
	TriBFace* pTBFTest = dynamic_cast<TriBFace*>(*iterBF);
	if (!pTBFTest->isValid()) {
	  // This is a quad bface instead
	  return 0;
	}
	if (pTBFTest->getBdryCondition() != iBC)
	  qSurfSmoothOK = false;
	else {
	  double adNormTest[3];
	  pTBFTest->getFace()->calcUnitNormal(adNormTest);
	  double dDot = fabs(dDOT3D(adNorm, adNormTest));
	  if (iFuzzyComp(dDot, 1.) != 0
	      && !(qOnSurf && BTE != pTBFTest->getTopologicalParent()))
	    qSurfSmoothOK = false;
	}
      }
      if (qSurfSmoothOK && m_pMesh->areBdryChangesAllowed()) {
	if (qOnSurf) {
	  Surface * surf = dynamic_cast<RefFace*>(BTE)->get_surface_ptr();
	  facetSurf = dynamic_cast<FacetSurface*>(surf);
	  CubitVector closest, normal;
	  facetSurf->closest_point(CubitVector(vert->getCoords()), &closest,
				   &normal);
	  normal.get_xyz(adNorm);
	}
	SMsetNormal(adNorm);
      }
      else {
	return 0;
      }
    }
    else if (vert->getVertType() == Vert::eBdryTwoSide) {
      // Can only smooth an internal bdry vertex if it lies in a coplanar
      // patch of the internal surface and if its incident cells fall into
      // only two regions.
      bool qSurfSmoothOK = true;
      // First check for region membership.
      int aiIncRegs[] =
	{ iInvalidRegion, iInvalidRegion };
      std::set<Cell*>::iterator iterC = m_incCells.begin();
      aiIncRegs[0] = (*iterC)->getRegion();
      int numRegions = 1;
      for (iterC++; iterC != m_incCells.end(); iterC++) {
	int iReg = (*iterC)->getRegion();
	if (iReg != aiIncRegs[0]) {
	  if (aiIncRegs[1] == iInvalidRegion) {
	    aiIncRegs[1] = iReg;
	    numRegions++;
	  }
	  else if (iReg != aiIncRegs[1]) {
	    // More than two regions incident; bail out now.
	    numRegions++;
	    return 0;
	  }
	}
      }
      // Make sure there really are two regions in this case.
      assert(numRegions == 1 || aiIncRegs[1] != iInvalidRegion);

      // Identify all faces that separate the two regions.
      std::set<TriFace*> spTFRegionBdry;
      std::set<Face*>::iterator iterF;
      for (iterF = m_nearFaces.begin(); iterF != m_nearFaces.end(); ++iterF) {
	Face *pF = *iterF;
	if (pF->hasVert(vert) && pF->getFaceLoc() == Face::eBdryTwoSide) {
	  assert(
	      pF->getLeftCell()->getRegion()
		  != pF->getRightCell()->getRegion());
	  if (pF->getType() != Face::eTriFace)
	    return 0;
	  spTFRegionBdry.insert(dynamic_cast<TriFace*>(pF));
	}
      }

      // Now check normals of all faces that separate the two regions.

      BasicTopologyEntity * BTE =
	  (*m_incBdryFaces.begin())->getTopologicalParent();
      if (BTE)
	qOnSurf = true;
      std::set<TriFace*>::iterator iterTF = spTFRegionBdry.begin();
      TriFace* pTF = *iterTF;
      double adNorm[3];
      pTF->calcUnitNormal(adNorm);
      for (iterTF++; iterTF != spTFRegionBdry.end(); iterTF++) {
	TriFace* pTFTest = *iterTF;
	double adNormTest[3];
	pTFTest->calcUnitNormal(adNormTest);
	double dDot = fabs(dDOT3D(adNorm, adNormTest));
	if (iFuzzyComp(dDot, 1.) != 0 && !qOnSurf)
	  qSurfSmoothOK = false;
      }
      if (qSurfSmoothOK && m_pMesh->areBdryChangesAllowed()) {
	// Set stuff up so that the internal bdry point remains on the
	// internal bdry.
	if (qOnSurf) {
	  Surface * surf = dynamic_cast<RefFace*>(BTE)->get_surface_ptr();
	  facetSurf = dynamic_cast<FacetSurface*>(surf);
	  CubitVector closest, normal;
	  facetSurf->closest_point(CubitVector(vert->getCoords()), &closest,
				   &normal);
	  normal.get_xyz(adNorm);
	}
	SMsetNormal(adNorm);
	qBdryVert = true;
      }
      else {
	return 0;
      }
    }

    // If there are too many neighbors, punt
    if (m_neighVerts.size() > iMaxVertNeigh
	|| m_incCells.size() > 2 * iMaxVertNeigh) {
      return 0;
    }

    // List the coordinates of all vertices
    int ii = 0;
    std::set<Vert*>::iterator iterV;
    std::map<Vert*, int> mVertIndex;
    for (iterV = m_neighVerts.begin(); iterV != m_neighVerts.end();
	++iterV, ii++) {
      Vert *vertTemp = *iterV;
      mVertIndex.insert(std::pair<Vert*, int>(vertTemp, ii));
      m_neighCoords[ii][0] = vertTemp->x();
      m_neighCoords[ii][1] = vertTemp->y();
      m_neighCoords[ii][2] = vertTemp->z();
    }

    // List all opposite faces (w/ correct orientation) using indices
    // into the location array
    std::set<Cell*>::iterator iterC;
    for (iterC = m_incCells.begin(), ii = 0; iterC != m_incCells.end();
	iterC++, ii++) {
      Cell *pC = *iterC;
      if (pC->getType() != Cell::eTet)
	return 0;
      Face *pF = dynamic_cast<TetCell*>(pC)->getOppositeFace(vert);
      if (pF->getLeftCell() == pC) {
	m_faceConn[ii][0] = (mVertIndex.find(pF->getVert(0)))->second;
	m_faceConn[ii][1] = (mVertIndex.find(pF->getVert(2)))->second;
	m_faceConn[ii][2] = (mVertIndex.find(pF->getVert(1)))->second;
      }
      else {
	m_faceConn[ii][0] = (mVertIndex.find(pF->getVert(0)))->second;
	m_faceConn[ii][1] = (mVertIndex.find(pF->getVert(1)))->second;
	m_faceConn[ii][2] = (mVertIndex.find(pF->getVert(2)))->second;
      }
    }

    double adPrevPt[] =
      { adCurrPt[0], adCurrPt[1], adCurrPt[2] };
    SMsmooth(m_neighVerts.size(), m_incCells.size(), adCurrPt, m_neighCoords,
	     m_faceConn, m_smoothData, (qBdryVert ? 1 : 0));
    if (qBdryVert && qOnSurf) {
      CubitVector closest;
      facetSurf->closest_point(CubitVector(adCurrPt), &closest);
      closest.get_xyz(adCurrPt);
    }
    double dDistMoved = dDIST3D(adCurrPt, adPrevPt);

    if (dDistMoved > 1.e-8) {
      vert->setCoords(3, adCurrPt);
      m_pMesh->moveVertEvent(vert);
      return 1;
    }
    else
      return 0;
  }

#ifdef HAVE_MESQUITE
  using namespace ITAPS_GRUMMP;
  using namespace Mesquite;

  SmoothingManagerMesquite3D::SmoothingManagerMesquite3D(VolMesh* const pM3D,
							 bool bIsoTarget) :
      SmoothingManager(pM3D)
  {

    bITarget = bIsoTarget;
    Mesquite::MsqError MesqErr;
    //getting msqMesh (mesh from Mesquite) from pM#D
    iMesh_Volume *iMP = new iMesh_Volume(pM3D);
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
      vFatalError("failed to create a tag", "TMOP 3D smoothing initialization");
    }

    // Now set data for it, by looping over verts and transcribing.
    for (unsigned int iV = 0; iV < instance->pGB->getNumVertices(); iV++) {

      Vert *pV = instance->pGB->pMesh()->getVert(iV);
      int iType = pV->getVertType();
      if (iType == Vert::eBdryApex || iType == Vert::eBdryCurve
	  || iType == Vert::eBdryTwoSide || iType == Vert::eBdry)
	iMesh_setIntData(instance, pV, fixed_tag, 1, &iMeshErr);
      else
	iMesh_setIntData(instance, pV, fixed_tag, 0, &iMeshErr);
      assert(iMeshErr == iBase_SUCCESS);
    }

    msqMesh = new Mesquite::MsqIMesh(instance, root_set, iBase_REGION, MesqErr,
				     &fixed_tag);

    if (MSQ_CHKERR(MesqErr)) {
      if (msqMesh)
	delete msqMesh;
      std::cout << MesqErr << std::endl;
      vFatalError("failed to create a Mesquite mesh",
		  "TMOP 3D smoothing initialization");
    }
    if (MesqErr)
      assert(0);
  }

  int
  SmoothingManagerMesquite3D::smoothAllVerts(const int)
  {

    Mesquite::MsqError MesqErr;
    Queue.run_instructions(msqMesh, MesqErr);
    if (MesqErr)
      std::cout << MesqErr << std::endl;
    return 1;
  }

  /*
   * At the moment only rigth angle tetrahedra are set as targets,
   * due to the complexity of 3D, a higher range of targets is needed for this smoothing
   * approach to work accordingly.
   */
  void
  SmoothingManagerMesquite3D::vSetTargetfromMetric()
  {

    Mesquite::MsqError MesqErr;

    std::vector<Mesquite::Mesh::ElementHandle> elements;
    msqMesh->get_all_elements(elements, MesqErr);
    std::vector<double> lid1(m_pMesh->getNumCells());
    std::vector<double> lid2(m_pMesh->getNumCells());
    std::vector<double> lid3(m_pMesh->getNumCells());
    std::vector<double> lid4(m_pMesh->getNumCells());
    std::vector<double> lid5(m_pMesh->getNumCells());
    std::vector<double> lid6(m_pMesh->getNumCells());
    std::vector<double> lid7(m_pMesh->getNumCells());
    std::vector<double> lid8(m_pMesh->getNumCells());
    std::vector<double> lid9(m_pMesh->getNumCells());

    for (unsigned int i = 0; i < m_pMesh->getNumCells(); i++) {
      Cell *pC = m_pMesh->getCell(i);
      Vert *pV0 = pC->getVert(0);
      Vert *pV1 = pC->getVert(1);
      Vert *pV2 = pC->getVert(2);
      Vert *pV3 = pC->getVert(3);

      double adMetric[6];
      double adEigVals[3];
      double A[3][3];
      double V[3][3];
      double adEigVec0[3], adEigVec1[3], adEigVec2[3];

      for (unsigned ii = 0; ii < 6; ii++)
	adMetric[ii] = (pV0->getMetric(ii) + pV1->getMetric(ii)
	    + pV2->getMetric(ii) + pV3->getMetric(ii)) / 4.;

      A[0][0] = adMetric[0];
      A[0][1] = adMetric[1];
      A[0][2] = adMetric[3];

      A[1][0] = adMetric[1];
      A[1][1] = adMetric[2];
      A[1][2] = adMetric[4];

      A[2][0] = adMetric[3];
      A[2][1] = adMetric[4];
      A[2][2] = adMetric[5];

      EigenDecomposition(A, V, adEigVals);

      adEigVec0[0] = V[0][0];
      adEigVec0[1] = V[1][0];
      adEigVec0[2] = V[2][0];

      adEigVec1[0] = V[0][1];
      adEigVec1[1] = V[1][1];
      adEigVec1[2] = V[2][1];

      adEigVec2[0] = V[0][2];
      adEigVec2[1] = V[1][2];
      adEigVec2[2] = V[2][2];

      NORMALIZE3D(adEigVec0);
      NORMALIZE3D(adEigVec1);
      NORMALIZE3D(adEigVec2);
//
//		  adEigVec0[0] = adEigVec0[0]*sqrt(adEigVals[1]*adEigVals[2]);
//		  adEigVec0[1] = adEigVec0[1]*sqrt(adEigVals[1]*adEigVals[2]);
//		  adEigVec0[2] = adEigVec0[2]*sqrt(adEigVals[1]*adEigVals[2]);
//
//		  adEigVec1[0] = adEigVec1[0]*sqrt(adEigVals[0]*adEigVals[2]);
//		  adEigVec1[1] = adEigVec1[1]*sqrt(adEigVals[0]*adEigVals[2]);
//		  adEigVec1[2] = adEigVec1[2]*sqrt(adEigVals[0]*adEigVals[2]);
//
//		  adEigVec2[0] = adEigVec2[0]*sqrt(adEigVals[0]*adEigVals[1]);
//		  adEigVec2[1] = adEigVec2[1]*sqrt(adEigVals[0]*adEigVals[1]);
//		  adEigVec2[2] = adEigVec2[2]*sqrt(adEigVals[0]*adEigVals[1]);
//

      adEigVec0[0] = adEigVec0[0] * sqrt(1. / sqrt(adEigVals[0]));
      adEigVec0[1] = adEigVec0[1] * sqrt(1. / sqrt(adEigVals[0]));
      adEigVec0[2] = adEigVec0[2] * sqrt(1. / sqrt(adEigVals[0]));

      adEigVec1[0] = adEigVec1[0] * sqrt(1. / sqrt(adEigVals[1]));
      adEigVec1[1] = adEigVec1[1] * sqrt(1. / sqrt(adEigVals[1]));
      adEigVec1[2] = adEigVec1[2] * sqrt(1. / sqrt(adEigVals[1]));

      adEigVec2[0] = adEigVec2[0] * sqrt(1. / sqrt(adEigVals[2]));
      adEigVec2[1] = adEigVec2[1] * sqrt(1. / sqrt(adEigVals[2]));
      adEigVec2[2] = adEigVec2[2] * sqrt(1. / sqrt(adEigVals[2]));

      double dMagEigVec0 = dMAG3D(adEigVec0);
      double dMagEigVec1 = dMAG3D(adEigVec1);
      double dMagEigVec2 = dMAG3D(adEigVec2);

      double adT0[3] =
	{ 0, 0, 0 };
      double adT1[3] =
	{ 0, 0, 0 };
      double adT2[3] =
	{ 0, 0, 0 };
      double adT3[3] =
	{ 0, 0, 0 };

      if (dMagEigVec0 >= dMagEigVec1 && dMagEigVec0 >= dMagEigVec2) {
	for (int ii = 0; ii < 3; ii++) {
	  adT1[ii] = adEigVec0[ii];
	  adT2[ii] = adEigVec1[ii];
	  adT3[ii] = adEigVec2[ii];
	}

	if (checkOrient3D(adT0, adT1, adT2, adT3) != 1) {
	  for (int ii = 0; ii < 3; ii++) {
	    adT2[ii] = adEigVec2[ii];
	    adT3[ii] = adEigVec1[ii];
	  }
	}

      }

      if (dMagEigVec1 >= dMagEigVec0 && dMagEigVec1 >= dMagEigVec2) {
	for (int ii = 0; ii < 3; ii++) {
	  adT1[ii] = adEigVec1[ii];
	  adT2[ii] = adEigVec0[ii];
	  adT3[ii] = adEigVec2[ii];
	}
	if (checkOrient3D(adT0, adT1, adT2, adT3) != 1) {
	  for (int ii = 0; ii < 3; ii++) {
	    adT2[ii] = adEigVec2[ii];
	    adT3[ii] = adEigVec0[ii];
	  }
	}
      }

      if (dMagEigVec2 >= dMagEigVec0 && dMagEigVec2 >= dMagEigVec0) {

	for (int ii = 0; ii < 3; ii++) {
	  adT1[ii] = adEigVec2[ii];
	  adT2[ii] = adEigVec0[ii];
	  adT3[ii] = adEigVec1[ii];
	}
	if (checkOrient3D(adT0, adT1, adT2, adT3) != 1) {
	  for (int ii = 0; ii < 3; ii++) {
	    adT2[ii] = adEigVec1[ii];
	    adT3[ii] = adEigVec0[ii];
	  }
	}
      }

      lid1[i] = adT1[0];
      lid2[i] = adT2[0];
      lid3[i] = adT3[0];
      lid4[i] = adT1[1];
      lid5[i] = adT2[1];
      lid6[i] = adT3[1];
      lid7[i] = adT1[2];
      lid8[i] = adT2[2];
      lid9[i] = adT3[2];

      if (adMetric[1] < 1e-6 && adMetric[1] > -1e-6 && adMetric[3] < 1e-6
	  && adMetric[4] > -1e-6 && adMetric[4] < 1e-6 && adMetric[3] > -1e-6) {
	if (fabs(adMetric[0] - adMetric[2]) < 1e-6
	    && fabs(adMetric[2] - adMetric[5]) < 1e-6
	    && fabs(adMetric[5] - adMetric[0]) < 1e-6) {
	  lid1[i] = 0.5;
	  lid2[i] = -0.5;
	  lid3[i] = 0;
	  lid4[i] = sqrt(2) / 2.;
	  lid5[i] = sqrt(2) / 2.;
	  lid6[i] = 0;
	  lid7[i] = 0.5;
	  lid8[i] = 0.5;
	  lid9[i] = 1;
	}
      }

    }
    const char LOCAL_ID_NAME1[] = "LOCAL_ID1";
    const char LOCAL_ID_NAME2[] = "LOCAL_ID2";
    const char LOCAL_ID_NAME3[] = "LOCAL_ID3";
    const char LOCAL_ID_NAME4[] = "LOCAL_ID4";
    const char LOCAL_ID_NAME5[] = "LOCAL_ID5";
    const char LOCAL_ID_NAME6[] = "LOCAL_ID6";
    const char LOCAL_ID_NAME7[] = "LOCAL_ID7";
    const char LOCAL_ID_NAME8[] = "LOCAL_ID8";
    const char LOCAL_ID_NAME9[] = "LOCAL_ID9";
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
    Mesquite::TagHandle lid_tag5 = msqMesh->tag_create(LOCAL_ID_NAME5,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag6 = msqMesh->tag_create(LOCAL_ID_NAME6,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag7 = msqMesh->tag_create(LOCAL_ID_NAME7,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag8 = msqMesh->tag_create(LOCAL_ID_NAME8,
						       Mesquite::Mesh::DOUBLE,
						       1, NULL, MesqErr);
    Mesquite::TagHandle lid_tag9 = msqMesh->tag_create(LOCAL_ID_NAME9,
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
    msqMesh->tag_set_vertex_data(lid_tag5, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid5), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag6, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid6), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag7, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid7), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag8, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid8), MesqErr);
    msqMesh->tag_set_vertex_data(lid_tag9, num_cells,
				 Mesquite::arrptr(elements),
				 Mesquite::arrptr(lid9), MesqErr);
  }
#endif

}
