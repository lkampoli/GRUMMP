#include <cstdio>
#include <set>
#include <list>
#include <deque>

#include "GR_GRCurve.h"
#include "GR_GRGeom2D.h"
#include "GR_GRPoint.h"
#include "GR_Mesh2D.h"

///// Read 2D mesh file and coordinate bdry entities with a geometry.
//Mesh2D::Mesh2D(const char strFileName[],
//	       GRGeom2D* const geometry,
//	       const int iQualMeas)
//  : Mesh(), m_ECEdgeF(), m_ECEdgeBF(), m_ECTri(), m_ECQuad(), m_ECIntEdgeBF()
//{
//
//
//  // These will both be factored out eventually.
//  m_qual = new Quality(this, iQualMeas);
//
//  // Read the actual data
//
//  readFromFile(strFileName);
//
//  relateMeshToGeometry(geometry);
//
//
//}

void
Mesh2D::readFromFileTerrible(const char strFileName[])
{
  if (strcasestr(strFileName, ".vtk") != NULL) {
    // Read a VTK legacy file
    // TO DO: Need an actual clearMeshData function, at some point...
    //     clearMeshData();
    readVTK(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".grmesh") != NULL) {
    // Read a GRUMMP native file
    //     clearMeshData();
    readNativeGRMESH(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".mesh") != NULL) {
    // Read a GRUMMP native file
    //     clearMeshData();
    readNativeMESH(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".poly") != NULL) {
    // Read a poly file
    //     clearMeshData();
    readTriangle(*this, strFileName);
  }
  else {
    // File format not currently recognized
    vFatalError("File format not recognized", "Mesh2D::readMeshFile");
  }

//	identifyVertexTypesGeometrically();
//
//	setAllHintFaces ();
//	setVertFaceNeighbors();
}

void
Mesh2D::readFromFile(const char strFileName[], int* /*numVars*/,
		     double** /*data*/)
{
  if (strcasestr(strFileName, ".vtk") != NULL) {
    // Read a VTK legacy file
    // TO DO: Need an actual clearMeshData function, at some point...
    //     clearMeshData();
    readVTK(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".grmesh") != NULL) {
    // Read a GRUMMP native file
    //     clearMeshData();
    readNativeGRMESH(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".mesh") != NULL) {
    // Read a GRUMMP native file
    //     clearMeshData();
    readNativeMESH(*this, strFileName);
  }
  else if (strcasestr(strFileName, ".poly") != NULL) {
    // Read a poly file
    //     clearMeshData();
    readTriangle(*this, strFileName);
  }
  else {
    // File format not currently recognized
    vFatalError("File format not recognized", "Mesh2D::readMeshFile");
  }

  identifyVertexTypesGeometrically();

  setAllHintFaces();
  setVertFaceNeighbors();
}

void
Mesh2D::relateMeshToGeometry(GRGeom2D* const geometry)
{
  if (!geometry)
    return;

  if (geometry->initialized()) {
    std::list<GRCurve*> curve_list;
    geometry->get_curves(curve_list);
    std::list<GRCurve*>::iterator it_curve, it_curve_end;
    std::list<GRPoint*> point_list;
    geometry->get_points(point_list);
    std::list<GRPoint*>::iterator it_point, it_point_end;

    BdryEdgeBase* bedge;
    std::set<BdryEdgeBase*> all_bedges;

    for (size_t i = 0; i < getNumBdryFaces(); i++) {
      bedge = dynamic_cast<BdryEdgeBase*>(getBFace(i));
      assert(bedge);
      all_bedges.insert(bedge);
    }
    for (GR_index_t i = 0; i < getNumIntBdryFaces(); i++) {
      bedge = dynamic_cast<BdryEdgeBase*>(getIntBFace(i));
      assert(bedge);
      all_bedges.insert(bedge);
    }
    std::deque<BdryEdgeBase*> bedge_chain;
    std::deque<BdryEdgeBase*>::iterator it_chain;

    GRPoint* point;
    GRCurve* curve;

    Vert* vert[2];
    Face* face;
    double vert_param[2];
    CubitVector vert_coord[2];

    BdryEdgeBase* bdry_edge[2];

    while (!all_bedges.empty()) {

      //Finding the first bdry edge of the chain.

      bdry_edge[0] = bdry_edge[1] = *(all_bedges.begin());
      all_bedges.erase(all_bedges.begin());

      vert[0] = bdry_edge[0]->getVert(0);
      vert[1] = bdry_edge[1]->getVert(1);
      vert_coord[0].set(vert[0]->x(), vert[0]->y(), 0.);
      vert_coord[1].set(vert[1]->x(), vert[1]->y(), 0.);

      it_curve = curve_list.begin();
      it_curve_end = curve_list.end();
      do {
	curve = *it_curve;
	vert_coord[0].set(vert[0]->x(), vert[0]->y(), 0.);
	vert_coord[1].set(vert[1]->x(), vert[1]->y(), 0.);
	if (curve->get_curve_geom()->coord_on_curve(vert_coord[0],
						    &vert_param[0])
	    && curve->get_curve_geom()->coord_on_curve(vert_coord[1],
						       &vert_param[1]))
	  break;
      }
      while (++it_curve != it_curve_end);
      assert(it_curve != it_curve_end);
      curve_list.erase(it_curve);

      bdry_edge[0]->setCurve(curve);
      bdry_edge[0]->setVert0Param(vert_param[0]);
      bdry_edge[0]->setVert1Param(vert_param[1]);

      assert(bedge_chain.empty());
      bedge_chain.push_back(bdry_edge[0]);

      //Building the chain.
      do {
	assert(vert[0] || vert[1]);
	for (int i = 0; i < 2; i++) {
	  if (!vert[i])
	    continue;

	  if (iFuzzyComp(vert_param[i], curve->get_curve_geom()->min_param())
	      == 0
	      || iFuzzyComp(vert_param[i], curve->get_curve_geom()->max_param())
		  == 0) {
	    // A curve end point
	    if (!vert[i]->getParentEntity()) {
	      it_point = point_list.begin();
	      it_point_end = point_list.end();
	      do {
		point = *it_point;
		if (iFuzzyComp(vert[i]->x(), point->X()) == 0
		    && iFuzzyComp(vert[i]->y(), point->Y()) == 0) {
		  vert[i]->setParentEntity(point);
		  break;
		}
	      }
	      while (++it_point != it_point_end);
	      if (it_point != it_point_end)
		point_list.erase(it_point);
	    }
	    // Done getting parent entity for vert at curve end point.
	    vert[i] = NULL;
	  }
	  else {
	    // Interior to a curve
	    if (!vert[i]->getParentEntity())
	      vert[i]->setParentEntity(curve);

	    std::set<BFace*> neigh_bfaces;
	    std::set<BFace*>::iterator it_bface;
	    {
	      bool is_bdry_vert;
	      std::set<Cell*> dummy_cell_set;
	      std::set<Vert*> dummy_vert_set;

	      findNeighborhoodInfo(vert[i], dummy_cell_set, dummy_vert_set,
				   &neigh_bfaces, &is_bdry_vert);
	      assert(is_bdry_vert);
	    }
	    if (neigh_bfaces.size() != 2) {
	      // Can this ever happen?
	      vert[i] = NULL;
	      continue;
	    }

	    it_bface = neigh_bfaces.begin();

	    // I'm nervous about *++it_bface evaluating to
	    // *(++it_bface).  Also, does C++ guarantee the order in
	    // which a series of variables are initialized in a single
	    // statement?  Not sure.  So I'm splitting what used to be
	    // one statement into several here.
	    BdryEdgeBase *edge_cand1, *edge_cand2;
	    edge_cand1 = dynamic_cast<BdryEdgeBase*>(*it_bface);
	    it_bface++;
	    assert(it_bface != neigh_bfaces.end());
	    edge_cand2 = dynamic_cast<BdryEdgeBase*>(*it_bface);
	    assert(edge_cand1 && edge_cand2);
	    assert(edge_cand1 == bdry_edge[i] || edge_cand2 == bdry_edge[i]);

	    bdry_edge[i] =
		(bdry_edge[i] == edge_cand1 ? edge_cand2 : edge_cand1);
	    assert(bdry_edge[i]->hasVert(vert[i]));

	    all_bedges.erase(bdry_edge[i]);

	    if (i == 0)
	      bedge_chain.push_front(bdry_edge[i]);
	    else
	      bedge_chain.push_back(bdry_edge[i]);

	    bdry_edge[i]->setCurve(curve);

	    vert[i] = (
		vert[i] == bdry_edge[i]->getVert(0) ?
		    bdry_edge[i]->getVert(1) : bdry_edge[i]->getVert(0));

	    vert_coord[i].set(vert[i]->x(), vert[i]->y(), 0.);

	    double this_param = vert_param[i];
	    double curve_param = -10;
	    bool result = curve->get_curve_geom()->coord_on_curve(vert_coord[i],
								  &curve_param);
	    assert(result);
	    vert_param[i] = curve_param;

	    if (vert[i] == bdry_edge[i]->getVert(0)) {
	      bdry_edge[i]->setVert0Param(vert_param[i]);
	      bdry_edge[i]->setVert1Param(this_param);
	    }
	    else {
	      bdry_edge[i]->setVert0Param(this_param);
	      bdry_edge[i]->setVert1Param(vert_param[i]);
	    }
	  } // Done with the case of a vert inside the curve
	} // Done checking both ends of the current chain
      }
      while (vert[0] || vert[1]);
      // Exit from the loop when the chain covers an entire curve

      assert(!vert[0] && !vert[1]);
      int param_compare;

      if (bedge_chain.size() >= 2) {

	BdryEdgeBase *bedge1, *bedge2;
	double min1, min2, max1, max2;

	for (int i = 0; i < 2; i++) {
	  switch (i)
	    {
	    case 0:
	      it_chain = bedge_chain.begin();
	      bedge1 = *it_chain;
	      bedge2 = *(++it_chain);
	      break;
	    case 1:
	      it_chain = --bedge_chain.end();
	      bedge1 = *it_chain;
	      bedge2 = *(--it_chain);
	      break;
	    default:
	      assert(0);
	      break;
	    }

	  min1 = std::min(bedge1->getVert0Param(), bedge1->getVert1Param());
	  max1 = std::max(bedge1->getVert0Param(), bedge1->getVert1Param());
	  min2 = std::min(bedge2->getVert0Param(), bedge2->getVert1Param());
	  max2 = std::max(bedge2->getVert0Param(), bedge2->getVert1Param());

	  if (iFuzzyComp(max1, max2) == 0) {
	    assert(iFuzzyComp(min1, max1) == -1);

	    if (iFuzzyComp(max1, bedge1->getVert0Param()) == 0)
	      bedge1->setVert1Param(curve->get_curve_geom()->max_param());
	    else
	      bedge1->setVert0Param(curve->get_curve_geom()->max_param());
	  }

	  if (iFuzzyComp(min1, min2) == 0) {
	    assert(iFuzzyComp(min1, max1) == -1);
	    if (iFuzzyComp(min1, bedge1->getVert0Param()) == 0)
	      bedge1->setVert1Param(curve->get_curve_geom()->min_param());
	    else
	      bedge1->setVert0Param(curve->get_curve_geom()->min_param());
	  }
	}
      }

      for (it_chain = bedge_chain.begin(); it_chain != bedge_chain.end();
	  ++it_chain) {
	bedge = *it_chain;
	vert_param[0] = bedge->getVert0Param();
	vert_param[1] = bedge->getVert1Param();
	param_compare = iFuzzyComp(vert_param[0], vert_param[1]);

	for (int i = 0; i < bedge->getNumFaces(); i++) {
	  face = bedge->getFace(i);

	  if ((face->getVert(0) == bedge->getVert(0) && param_compare == 1)
	      || (face->getVert(0) == bedge->getVert(1) && param_compare == -1)) {

	    //face is inverted with respect to curve.
	    Cell *cell_left = face->getLeftCell(), *cell_right =
		face->getRightCell();
	    Vert *beg_vert = face->getVert(0), *end_vert = face->getVert(1);

	    face->interchangeCellsAndVerts();

#ifndef NDEBUG
	    assert(face->isValid());
	    assert(cell_left->isValid());
	    assert(cell_right->isValid());
	    assert(beg_vert->isValid());
	    assert(end_vert->isValid());
#endif
	  }

	  if (param_compare == -1)
	    assert(bedge->isForward());
	  else if (param_compare == 1)
	    assert(!bedge->isForward());
	  else
	    assert(0);

	} // Done checking orientation of faces for this bdry edge
      } // Done checking the whole chain for this curve

      logMessage(4, "curve = %p, min_param = %lf, max_param = %lf\n", curve,
		 curve->get_curve_geom()->min_param(),
		 curve->get_curve_geom()->max_param());

      for (std::deque<BdryEdgeBase*>::iterator it_deque = bedge_chain.begin();
	  it_deque != bedge_chain.end(); ++it_deque) {
	BFace* my_bface = *it_deque;
	BdryEdgeBase* my_bedge = dynamic_cast<BdryEdgeBase*>(my_bface);
	logMessage(
	    4,
	    "bdry_edge = %p, curve = %p, vert0_param = %lf, vert1_param = %lf, parent0 = %p, parent1 = %p\n",
	    my_bedge, my_bedge->getCurve(), my_bedge->getVert0Param(),
	    my_bedge->getVert1Param(), my_bedge->getVert(0)->getParentEntity(),
	    my_bedge->getVert(1)->getParentEntity());
      }
      logMessage(2, "\n");

      bedge_chain.clear();

    } // Done assigning all bdry edges to a curve

    m_isLengthScaleFromCellSizes = false;
  } // Done matching up the mesh to an existing geometry

  else {
    logMessage(3, "Recreating boundary data by inference\n");
    // Recreate boundary data by inference.
    // Potentially a bad issue on curved surfaces IF a vert is connected to two close-to colinear bdry edges,
    // but is identified as Vert::eBdry
    typedef std::map<Vert*, GRPoint*> VertToPoint;
    typedef VertToPoint::iterator itVTP;
    VertToPoint vert_to_point;
    // used to connect all colinear bdry faces together...
    // This is based on identifyVertexTypesGeometrically no longer
    // classifying everything as a bdryApex.. if a vert is
    // not a bdryApex, it must have two neighboring colinear bdryfaces
    // we can traverse them to find apexes, and then build the chain as a line

    // visitedBFaces contains any BFace we've seen (so we don't start from a BFace on a chain
    std::set<BFace*> visitedBFaces;
    // this contains the BFaces on our local chain, to join up with a single GR_LINE
    std::set<BFace*> chainBFaces;

    std::set<IntBdryEdge*> visitedIntBEdges;
    std::set<IntBdryEdge*> chainIntBEdges;

    // used for findNeighborhoodInfo
    std::set<Cell*> spCIncident;
    std::set<Vert*> spVCand;
    std::set<BFace*> spBFIncident;
    // Now the internal bdry faces, if any.

    for (GR_index_t iIBF = 0; iIBF < getNumIntBdryFaces(); iIBF++) {
      // This array has to have live data, or we wouldn't be here.
      IntBdryEdge * pIBE = getIntBFace(iIBF);
      if (visitedIntBEdges.count(pIBE) == 0) {

	visitedIntBEdges.insert(pIBE);
	chainIntBEdges.insert(pIBE);
	Face *pF0 = pIBE->getFace(0);
	Face *pF1 = pIBE->getFace(1);
	Vert *pV0 = pIBE->getVert(0), *pV1 = pIBE->getVert(1);
	IntBdryEdge * pIBE0 = pIBE, *pIBE1 = pIBE;
	bool qBdryVert = true;
	// Declare that our chain goes from pV0 to pV1, and is complete when we've hit two apexes
	while (pV0->getVertType() != Vert::eBdryApex
	    || pV1->getVertType() != Vert::eBdryApex) {
	  // If this vert isn't an apex, keep the chain going using pIBE0, pIBE1
	  if (pV0->getVertType() != Vert::eBdryApex) {
	    // get the other BFace connected to this vert
	    findNeighborhoodInfo(pV0, spCIncident, spVCand, &spBFIncident,
				 &qBdryVert);
	    assert(spBFIncident.size() == 2);

	    // continue our chain along the BFace we aren't already on
	    IntBdryEdge * incBFace =
		dynamic_cast<IntBdryEdge*>(*spBFIncident.begin());
	    if (incBFace == pIBE0)
	      incBFace = dynamic_cast<IntBdryEdge*>(*(++spBFIncident.begin()));
	    pIBE0 = incBFace;
	    pV0 = incBFace->getFace(0)->getOppositeVert(pV0);
	    visitedIntBEdges.insert(pIBE0);
	    chainIntBEdges.insert(pIBE0);
	  }
	  if (pV1->getVertType() != Vert::eBdryApex) {
	    findNeighborhoodInfo(pV1, spCIncident, spVCand, &spBFIncident,
				 &qBdryVert);
	    assert(spBFIncident.size() == 2);
	    IntBdryEdge * incBFace =
		dynamic_cast<IntBdryEdge*>(*spBFIncident.begin());
	    if (incBFace == pIBE1)
	      incBFace = dynamic_cast<IntBdryEdge*>(*(++spBFIncident.begin()));
	    pIBE1 = incBFace;
	    pV1 = incBFace->getFace(0)->getOppositeVert(pV1);
	    visitedIntBEdges.insert(pIBE1);
	    chainIntBEdges.insert(pIBE1);
	  }
	}

	// build the GR_Line from start to finish
	GRPoint* point0 = static_cast<GRPoint*>(NULL);
	GRPoint* point1 = static_cast<GRPoint*>(NULL);

	std::pair<itVTP, bool> ins_pair0 = vert_to_point.insert(
	    std::make_pair(pV0, point0));
	std::pair<itVTP, bool> ins_pair1 = vert_to_point.insert(
	    std::make_pair(pV1, point1));
	if (ins_pair0.second) {
	  point0 = geometry->build_point(pV0->x(), pV0->y());
	  ins_pair0.first->second = point0;
	}
	else {
	  point0 = ins_pair0.first->second;
	  assert(point0);
	}

	if (ins_pair1.second) {
	  point1 = geometry->build_point(pV1->x(), pV1->y());
	  ins_pair1.first->second = point1;
	}
	else {
	  point1 = ins_pair1.first->second;
	  assert(point1);
	}
	Cell *pCL = pF0->getOppositeCell(pIBE);
	int iLeftReg = pCL->getRegion();
	Cell *pCR = pF1->getOppositeCell(pIBE);
	int iRightReg = pCR->getRegion();

	//Create new polyline between point0 and point1
	GRCurve * curve = geometry->build_line(point0, point1, iLeftReg,
					       iRightReg, 0, 0);
	// for each BFace, connect it to the curve and set the parameters accordingly
	for (std::set<IntBdryEdge*>::iterator itB = chainIntBEdges.begin();
	    itB != chainIntBEdges.end(); ++itB) {
	  IntBdryEdge *tmpIntBEdge = *itB;

	  Face * F0 = tmpIntBEdge->getFace(0);
	  Face * F1 = tmpIntBEdge->getFace(1);

	  tmpIntBEdge->setCurve(curve);
	  double param0 = curve->u_from_position(
	      CubitVector(tmpIntBEdge->getVert(0)->getCoords()));
	  double param1 = curve->u_from_position(
	      CubitVector(tmpIntBEdge->getVert(1)->getCoords()));

	  if (tmpIntBEdge->isForward()) {

	  }
	  else {
	    F0->interchangeCellsAndVerts();
	    F1->interchangeCellsAndVerts();

	  }

	  tmpIntBEdge->setVert0Param(param0);
	  tmpIntBEdge->setVert1Param(param1);
	  if (tmpIntBEdge->getVert(0)->getVertType() != Vert::eBdryApex)
	    tmpIntBEdge->getVert(0)->setParentEntity(curve);
	  if (tmpIntBEdge->getVert(1)->getVertType() != Vert::eBdryApex)
	    tmpIntBEdge->getVert(1)->setParentEntity(curve);
	  tmpIntBEdge->getVert0Param();
	  tmpIntBEdge->getVert1Param();

	}
	chainIntBEdges.clear();
      }
    }
    assert(visitedIntBEdges.size() == getNumIntBdryFaces());
    visitedIntBEdges.clear();

    // This works the same as internal boundaries
    for (GR_index_t iBF = 0; iBF < getNumBdryFaces(); iBF++) {
      BFace *pBF = getBFace(iBF);
      if (visitedBFaces.count(pBF) == 0) {

	visitedBFaces.insert(pBF);
	chainBFaces.insert(pBF);
	Face *pF = pBF->getFace(0);
	Vert *pV0 = pF->getVert(0), *pV1 = pF->getVert(1);

	BFace * pBF0 = pBF, *pBF1 = pBF;
	bool qBdryVert = true;
	while (pV0->getVertType() != Vert::eBdryApex
	    || pV1->getVertType() != Vert::eBdryApex) {
	  if (pV0->getVertType() != Vert::eBdryApex) {

	    findNeighborhoodInfo(pV0, spCIncident, spVCand, &spBFIncident,
				 &qBdryVert);
	    assert(spBFIncident.size() == 2);
	    BFace * incBFace = *spBFIncident.begin();
	    if (incBFace == pBF0)
	      incBFace = *(++spBFIncident.begin());
	    pBF0 = incBFace;
	    pV0 = incBFace->getFace(0)->getOppositeVert(pV0);
	    visitedBFaces.insert(pBF0);
	    chainBFaces.insert(pBF0);
	  }
	  if (pV1->getVertType() != Vert::eBdryApex) {
	    findNeighborhoodInfo(pV1, spCIncident, spVCand, &spBFIncident,
				 &qBdryVert);
	    assert(spBFIncident.size() == 2);
	    BFace * incBFace = *spBFIncident.begin();
	    if (incBFace == pBF1)
	      incBFace = *(++spBFIncident.begin());
	    pBF1 = incBFace;
	    pV1 = incBFace->getFace(0)->getOppositeVert(pV1);
	    visitedBFaces.insert(pBF1);
	    chainBFaces.insert(pBF1);
	  }
	}
	// Now to build them using connectedVerts;

	GRPoint* point0 = static_cast<GRPoint*>(NULL);
	GRPoint* point1 = static_cast<GRPoint*>(NULL);

	std::pair<itVTP, bool> ins_pair0 = vert_to_point.insert(
	    std::make_pair(pV0, point0));
	std::pair<itVTP, bool> ins_pair1 = vert_to_point.insert(
	    std::make_pair(pV1, point1));
	if (ins_pair0.second) {
	  point0 = geometry->build_point(pV0->x(), pV0->y());
	  ins_pair0.first->second = point0;
	}
	else {
	  point0 = ins_pair0.first->second;
	  assert(point0);
	}

	if (ins_pair1.second) {
	  point1 = geometry->build_point(pV1->x(), pV1->y());
	  ins_pair1.first->second = point1;
	}
	else {
	  point1 = ins_pair1.first->second;
	  assert(point1);
	}
	int iBC = pBF->getBdryCondition();
	assert(iBC > 0);

	// Find the adjacent region.
	Cell *pC = pF->getOppositeCell(pBF);
	int iReg = pC->getRegion();
	GRCurve* curve = NULL;

	if (pF->getLeftCell() == pC)
	  // Cell is on the left, bdry on the right.
	  curve = geometry->build_line(point0, point1, iReg, 0, 0, iBC);
	else
	  // Cell is on the right, bdry on the left.
	  curve = geometry->build_line(point0, point1, 0, iReg, iBC, 0);
	for (std::set<BFace*>::iterator itB = chainBFaces.begin();
	    itB != chainBFaces.end(); ++itB) {
	  BFace *tmpBFace = *itB;
	  BdryEdge* bedge = dynamic_cast<BdryEdge*>(tmpBFace);
	  assert(bedge);
	  bedge->setCurve(curve);
//					if(bedge->isForward()){
//printf("%.16f %.16f | %.16f %.16f \n",point0->X(),point0->Y(),point1->X(),point1->Y());

	  bedge->setVert0Param(
	      curve->get_curve_geom()->param_at_coord(
		  CubitVector(tmpBFace->getVert(0)->getCoords())));
	  bedge->setVert1Param(
	      curve->get_curve_geom()->param_at_coord(
		  CubitVector(tmpBFace->getVert(1)->getCoords())));
//						printf("%f %f params %d\n",bedge->getVert0Param(),bedge->getVert1Param(),bedge->isForward());
//						bedge->getVert(0)->printVertInfo();
//						bedge->getVert(1)->printVertInfo();
//					} else {
//						bedge->setVert0Param(curve->get_curve_geom()->param_at_coord(CubitVector(tmpBFace->getVert(1)->getCoords())));
//						bedge->setVert1Param(curve->get_curve_geom()->param_at_coord(CubitVector(tmpBFace->getVert(0)->getCoords())));
//					}
	  if (tmpBFace->getVert(0)->getVertType() != Vert::eBdryApex)
	    tmpBFace->getVert(0)->setParentEntity(curve);
	  if (tmpBFace->getVert(1)->getVertType() != Vert::eBdryApex)
	    tmpBFace->getVert(1)->setParentEntity(curve);
	  bedge->getVert0Param();
	  bedge->getVert1Param();

	}
	chainBFaces.clear();
      }

    }
    assert(visitedBFaces.size() == getNumBdryFaces());
    visitedBFaces.clear();
    geometry->set_initialized();
    m_isLengthScaleFromCellSizes = true;

    //    assert(getNumIntBdryFaces() == 0);
    //		if (getNumIntBdryFaces() > 0) {
    //			for (GR_index_t iIBF = 0; iIBF < getNumIntBdryFaces(); iIBF++) {
    //				// This array has to have live data, or we wouldn't be here.
    //				IntBdryEdge * pIBE = getIntBdryEdge(iIBF);
    //				Face *pF0 = pIBE->getFace(0);
    //				Face *pF1 = pIBE->getFace(1);
    //
    //				// Need this both for checking and for setting up patches.
    //				Vert *pV0 = pIBE->getVert(0);
    //				Vert *pV1 = pIBE->getVert(1);
    //
    //				GRPoint* point0 = static_cast<GRPoint*>(NULL);
    //				GRPoint* point1 = static_cast<GRPoint*>(NULL);
    //
    //				std::pair<itVTP, bool> ins_pair0
    //				= vert_to_point.insert( std::make_pair(pV0, point0) );
    //				std::pair<itVTP, bool> ins_pair1
    //				= vert_to_point.insert( std::make_pair(pV1, point1) );
    //
    //				if( ins_pair0.second ) {
    //					point0 = geometry->build_point(pV0->x(), pV0->y());
    //					ins_pair0.first->second = point0;
    //				}
    //				else {
    //					point0 = ins_pair0.first->second;
    //					assert(point0);
    //				}
    //
    //				if( ins_pair1.second ) {
    //					point1 = geometry->build_point(pV1->x(), pV1->y());
    //					ins_pair1.first->second = point1;
    //				}
    //				else {
    //					point1 = ins_pair1.first->second;
    //					assert(point1);
    //				}
    //
    //				// Grab both region tags, working from IBF -> face -> cell.
    //				Cell *pCL = pF0->getOppositeCell(pIBE);
    //				int iLeftReg = pCL->getRegion();
    //				Cell *pCR = pF1->getOppositeCell(pIBE);
    //				int iRightReg = pCR->getRegion();
    //
    //				//Create new polyline between point0 and point1.
    //				GRCurve* curve = geometry->build_line(point0, point1,
    //						iLeftReg, iRightReg,
    //						0, 0);
    //				pIBE->setCurve(curve);
    //				pIBE->setVert0Param(curve->get_curve_geom()->min_param());
    //				pIBE->setVert1Param(curve->get_curve_geom()->max_param());
    //				pV0->setParentEntity(curve);
    //				pV1->setParentEntity(curve);
    //			}
    //
    //		} // Done with setting up internal bdry faces w/ face links

    // Now turn to regular bdry faces
    //		for (GR_index_t iBF = 0; iBF < getNumBdryFaces(); iBF++) {
    //			BFace *pBF = getBFace(iBF);
    //			Face *pF = pBF->getFace(0);
    //			Vert *pV0 = pF->getVert(0), *pV1 = pF->getVert(1);
    //
    //			GRPoint* point0 = static_cast<GRPoint*>(NULL);
    //			GRPoint* point1 = static_cast<GRPoint*>(NULL);
    //
    //			std::pair<itVTP, bool> ins_pair0
    //			= vert_to_point.insert( std::make_pair(pV0, point0) );
    //			std::pair<itVTP, bool> ins_pair1
    //			= vert_to_point.insert( std::make_pair(pV1, point1) );
    //
    //			if( ins_pair0.second ) {
    //				point0 = geometry->build_point(pV0->x(), pV0->y());
    //				ins_pair0.first->second = point0;
    //			}
    //			else {
    //				point0 = ins_pair0.first->second;
    //				assert(point0);
    //			}
    //
    //			if( ins_pair1.second ) {
    //				point1 = geometry->build_point(pV1->x(), pV1->y());
    //				ins_pair1.first->second = point1;
    //			}
    //			else {
    //				point1 = ins_pair1.first->second;
    //				assert(point1);
    //			}
    //
    //			// Grab the BC from the incoming data; this info will get
    //			// transfered to the new curve.
    //			int iBC = pBF->getBdryCondition();
    //			assert(iBC > 0);
    //
    //			// Find the adjacent region.
    //			Cell *pC = pF->getOppositeCell(pBF);
    //			int iReg = pC->getRegion();
    //
    //			GRCurve* curve = NULL;
    //
    //			if (pF->getLeftCell() == pC)
    //				// Cell is on the left, bdry on the right.
    //				curve = geometry->build_line(point0, point1, iReg, 0, 0, iBC);
    //			else
    //				// Cell is on the right, bdry on the left.
    //				curve = geometry->build_line(point0, point1, 0, iReg, iBC, 0);
    //
    //			BdryEdgeBase* bedge = dynamic_cast<BdryEdgeBase*>(pBF);
    //			assert(bedge);
    //
    //			bedge->setCurve(curve);
    //			bedge->setVert0Param(curve->get_curve_geom()->min_param());
    //			bedge->setVert1Param(curve->get_curve_geom()->max_param());
    //			pV0->setParentEntity(curve);
    //			pV1->setParentEntity(curve);
    //
    //		} // Done setting up all BFace stuff, including building patches.

    //		assert (geometry->num_curves() == getNumBdryFaces() + getNumIntBdryFaces());
    //		geometry->set_initialized();
    //		m_isLengthScaleFromCellSizes = true;
  }

} // Done relating mesh to geometry

