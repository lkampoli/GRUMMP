/*
 * InsertionPointCalculator.cxx
 *
 *  Created on: 2014-08-06
 *      Author: zaide
 */

#include "GR_InsertionManager.h"
#include "GR_SmoothingManager.h"
#include "GR_Cell.h"
#include "GR_TMOPQual.h"

namespace GRUMMP
{
// this one is obvious, and basically just another level to get to the point
  InsertionPointCalculator::ePointType
  CircumcenterIPC::calcInsertionPoint(const SimplexCell * const cell,
				      double adPointLoc[]) const
  {

    assert(cell->getType() == Cell::eTriCell || cell->getType() == Cell::eTet);
    cell->calcCircumcenter(adPointLoc);

    return eCircumcenter;
  }

// this is more complicated
  InsertionPointCalculator::ePointType
  OffcenterIPC2D::calcInsertionPoint(const SimplexCell * const cell,
				     double adPointLoc[]) const
  {

    assert(cell->getType() == Cell::eTriCell);
    cell->calcCircumcenter(adPointLoc);

    if (iFuzzyComp(cell->calcMinSolid() * M_PI / 180., m_alpha) > -1) {
      return eCircumcenter;
    }
    const Face * pFShortest = cell->getFace(0);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(1)->calcSize() ?
	    pFShortest : cell->getFace(1);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(2)->calcSize() ?
	    pFShortest : cell->getFace(2);

    // lets get the shortest edge, and preserve orientation
    const Vert * pVA = pFShortest->getVert(0);
    const Vert * pVB = pFShortest->getVert(1);

    assert(iFuzzyComp(dDist(pVA, pVB), cell->calcShortestEdgeLength()) == 0);

    // circumradius of shortestFace + circumcenter of cell
    const double * const adA = pVA->getCoords();
    const double * const adB = pVB->getCoords();
    double dAB2 = dDIST_SQ_2D(adA, adB);
//	double dR2 = dDIST_SQ_2D(adPointLoc,adA);
    // compute circumradius of A-B-Circumcenter using R^2 = (a*b*c)/4A and
    // isosceles triangle identities.
    double adAB[2] =
      { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };
    double adCircOffset[2];
//	double alpha = 55.*M_PI/180.;//2*0.447832396928932;

//	double tempCoeff = 0.5*sqrt(1./sin(alpha/2.)/sin(alpha/2.)-1.0);
    double t = m_coeff * sqrt(dAB2 / dDIST_SQ_2D(adPointLoc, adAB));
//	double t = tempCoeff*sqrt(dAB2/dDIST_SQ_2D(adPointLoc,adAB));

    adCircOffset[0] = adAB[0] + t * (adPointLoc[0] - adAB[0]);
    adCircOffset[1] = adAB[1] + t * (adPointLoc[1] - adAB[1]);

    //	if(dR2*dR2/(4.*dR2-dAB2) <= m_beta2*dAB2){
    //we are done!
    //		assert(cell->calcMinSolid()*M_PI/180. > m_alpha/2.*0.999999);

    //		assert(iFuzzyComp(t,1.0) > -1);
    //		adCircOffset[0] = adAB[0] + t*(adPointLoc[0]-adAB[0]);
    //		adCircOffset[1] = adAB[1] + t*(adPointLoc[1]-adAB[1]);
//	if(0){
    if (dDIST_SQ_2D(adAB,adCircOffset) > dDIST_SQ_2D(adAB, adPointLoc)) {
      return eCircumcenter;
    }
    else {
//		double delta2 = dDIST2D(adCircOffset,adPointLoc);
//		printf("delta2 %.16f %.16f\n",delta2,dDIST2D(adPointLoc,adA)*(1-1/sqrt(2.)));
      adPointLoc[0] = adCircOffset[0];
      adPointLoc[1] = adCircOffset[1];
      // choose offcenter as point D that lies on the bisector of AB and use
      // law of sines to determine the length of the legs based on alpha
      // and ||D-A||^2 = coeff ||A-B||^2
      // writing D = 0.5*(A+B)+t*(0.5*(A+B)-C) for parameter t, we can solve for t
      //		assert(cell->calcMinSolid()*M_PI/180. < m_alpha/2.*1.00000001);
      //
      //		double adAB[2] = {0.5*(adA[0] + adB[0]), 0.5*(adA[1] + adB[1])};
      //		double t = m_coeff*sqrt(dAB2/dDIST_SQ_2D(adPointLoc,adAB));
      //
      //		assert(iFuzzyComp(t,1.0) < 1);
      //
      //		adPointLoc[0] = adAB[0] + t*(adPointLoc[0]-adAB[0]);
      //		adPointLoc[1] = adAB[1] + t*(adPointLoc[1]-adAB[1]);

      // checks that the offcenter is the same distance from A, B
      //		assert(iFuzzyComp(dDIST_SQ_2D(adPointLoc,adA),dDIST_SQ_2D(adPointLoc,adB)) == 0);
      //		assert(iFuzzyComp(dDIST_SQ_2D(adPointLoc,adA),dAB2) > -1);
      //#ifndef NDEBUG
      //		double adRow0[] = {adA[0] - adB[0],
      //				adA[1] - adB[1]};
      //		double adRow1[] = {adA[0] - adPointLoc[0],
      //				adA[1] - adPointLoc[1]};
      //
      //		double dRHS0 = 0.5 * (dDOT2D (adA, adA) - dDOT2D (adB, adB));
      //		double dRHS1 = 0.5 * (dDOT2D (adA, adA) - dDOT2D (adPointLoc, adPointLoc));
      //
      //		double dDet = adRow0[0] * adRow1[1] - adRow1[0] * adRow0[1];
      //		double dDet0 = dRHS0 * adRow1[1] - dRHS1 * adRow0[1];
      //		double dDet1 = dRHS1 * adRow0[0] - dRHS0 * adRow1[0];
      //		double adCircSmall[2];
      //		adCircSmall[0] = dDet0 / dDet;
      //		adCircSmall[1] = dDet1 / dDet;
      //
      //
      //		double A4 = sqrt(dAB2)*sqrt(4.*dDIST_SQ_2D(adPointLoc,adA)-dAB2);
      //		double R2 = dDIST_SQ_2D(adPointLoc,adA)*dDIST_SQ_2D(adPointLoc,adA)*dAB2/(A4*A4);
      //		assert(iFuzzyComp(R2/dAB2,m_beta2) == 0);
      //		//				printf("%.16f %.16f\n",dDIST_SQ_2D(adCircSmall,adA),dDIST_SQ_2D(adCircSmall,adB));
      //		assert(iFuzzyComp(dDIST_SQ_2D(adCircSmall,adA),dDIST_SQ_2D(adCircSmall,adB)) == 0);
      //
      //		assert(iFuzzyComp(dDIST_SQ_2D(adCircSmall,adA),dDIST_SQ_2D(adCircSmall,adPointLoc)) == 0);
      //
      //#endif

      return eOffcenter;

    }
  }

// this is more complicated
  InsertionPointCalculator::ePointType
  OptimalSteinerIPC2D::calcInsertionPoint(const SimplexCell * const cell,
					  double adPointLoc[]) const
  {
    //this needs to figure out stuff. This is hard
    assert(cell->getType() == Cell::eTriCell);

    //		printf("%f min angle %f\n",cell->calcMinSolid()*M_PI/180., m_alpha );
    //	writeVTK_SingleCell(cell,"tempcell.vtk");

    cell->calcCircumcenter(adPointLoc);
    if (iFuzzyComp(cell->calcMinSolid() * M_PI / 180., m_alpha) > -1) {
      return eCircumcenter;
    }
    const Face * pFShortest = cell->getFace(0);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(1)->calcSize() ?
	    pFShortest : cell->getFace(1);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(2)->calcSize() ?
	    pFShortest : cell->getFace(2);
    // given the shortest face, find the petal of interest.
    // lets get the shortest edge, and preserve orientation
    const Vert * pVA = pFShortest->getVert(0);
    const Vert * pVB = pFShortest->getVert(1);

    assert(iFuzzyComp(dDist(pVA, pVB), cell->calcShortestEdgeLength()) == 0);

    const double * const adA = pVA->getCoords();
    const double * const adB = pVB->getCoords();

    double dCircRadius2 = dDIST_SQ_2D(adPointLoc, adA);
    double dAB2 = dDIST_SQ_2D(adA, adB);
    assert(iFuzzyComp(sqrt(dAB2), cell->calcShortestEdgeLength()) == 0);

    double adAB[2] =
      { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };

    double t = 0.95 * m_coeffRad * sqrt(dAB2 / dDIST_SQ_2D(adPointLoc, adAB));
    // this is circumcenter of the petal.
    double adPetalCirc[2];
    adPetalCirc[0] = adAB[0] + t * (adPointLoc[0] - adAB[0]);
    adPetalCirc[1] = adAB[1] + t * (adPointLoc[1] - adAB[1]);
    // makes sure this is equidistance from points A and B and equal to the prescribed radius
//	assert(iFuzzyComp(dDIST_SQ_2D(adA,adPetalCirc),m_petalRad2*dAB2) == 0);
//	assert(iFuzzyComp(dDIST_SQ_2D(adB,adPetalCirc),m_petalRad2*dAB2) == 0);

    // Type I. Need to find new point location
    double tOffcenter = 0.95 * m_coeffOff
	* sqrt(dAB2 / dDIST_SQ_2D(adPointLoc, adAB));
    //	printf("%f toff %f tnormal\n",tOffcenter,t);
    if (iFuzzyComp(tOffcenter, 1.0) < 1) {
//		assert(iFuzzyComp(tOffcenter,1.0) < 1);
      adPointLoc[0] = adAB[0] + tOffcenter * (adPointLoc[0] - adAB[0]);
      adPointLoc[1] = adAB[1] + tOffcenter * (adPointLoc[1] - adAB[1]);
      //		printf("offcenter %f %f | %f |dist %.16f\n",adPointLoc[0],adPointLoc[1], t,dDIST2D(adPetalCirc,adPointLoc));
      return eOffcenter;
    }

    //	printf("circumcenter %f %f | t %f |dist %f\n",adPointLoc[0],adPointLoc[1], t,dDIST2D(adAB,adPointLoc));
    //	printf("Petal info %f %f | rad %.16f |dist %f\n",adPetalCirc[0],adPetalCirc[1],sqrt(m_petalRad2*dAB2),dDIST2D(adAB,adPetalCirc));
    //	assert(iFuzzyComp(t,1.0) < 1);

    //	assert(isInCircle(adPetalCirc,adA,adB,adPointLoc) > -1);
    // now we have the petal, can begin to figure things out.
    // need to find the point in the petal that is furthest away from all neighboring verts.
    // if not, we need to look at voronoi in detail
    // every face is a perpendicular bisector of an edge
    // given the two neighboring cells of a face, the voronoi edge connects their circumcenters

    // Will never explicitly store voronoi edges, other than given a face,
    // if both its neighboring cells have been checked, their circumcenters are stored
    // and the edge can be reconstructed from them.

    // store the "winning" pointType and its distance away
    ePointType pointType = eCircumcenter;
    // distance away from neighboring verts, squared
    // its location is already in adPointLoc
    double furthestDist2 = dCircRadius2;

    double adIntersect[2]; // intersection point
    // circumcenters of cell0 and cell1 of the face to check
    // also gives endpoints of voronoi edge
    double adCircs[2][2];
    double distToCheck2; // new distance to check, squared
    // faces to check.
    std::queue<const Face*> facesToTest;
    std::set<const Face*> facesTested;

    // map of cells to circumcenters, to avoid recomputing circumcenters multiple times
    // also used to check if we've looked at a circumcenter as a candidate point

    // hacky way do to things. Could design simple struct, but use this for now
    std::map<const SimplexCell*, VoronoiNode> circCentMap;
    double adCirc[2] =
      { adPointLoc[0], adPointLoc[1] };
    circCentMap.insert(std::make_pair(cell, adCirc));
    //	printf("first cell %p\n",cell);
    // first two faces
    facesToTest.push(cell->getOppositeFace(pVA));
    facesToTest.push(cell->getOppositeFace(pVB));

    std::vector<VoronoiNode> ptsChecked;
    std::vector<double> dChecked;
    std::vector<ePointType> ptTypes;
//	ptsChecked.push_back(adPointLoc);
//	dChecked.push_back(furthestDist2);
//	ptTypes.push_back(CIRCUMCENTER);
    assert(
	iFuzzyComp(dDIST_SQ_2D(adCirc, adPetalCirc), (m_petalRad2 * dAB2)) < 1);
    while (!facesToTest.empty()) {

      const Face * face = facesToTest.front();
      //		assert(face != pFShortest);
      facesToTest.pop();
      // Don't bother if this face has been tested, or is a boundary face
      // Internal boundary faces complicate things, but lets worry about that
      // at another time.
      if (!facesTested.insert(face).second || face->isBdryFace())
	continue;

      // get the two cells on either side of a face
      const SimplexCell * cells[] =
	{ dynamic_cast<const SimplexCell*>(face->getCell(0)),
	    dynamic_cast<const SimplexCell*>(face->getCell(1)) };

      assert(cells[0] && cells[1]); // should always be valid, since we skip boundary faces

      // try and find them. should only be able to find one of them
      std::map<const SimplexCell*, VoronoiNode>::iterator cellIts[] =
	{ circCentMap.find(cells[0]), circCentMap.find(cells[1]) };
      // assertion to guarantee that
      //						printf("face %p | cells %p %p | %d \n",face,cells[0],cells[1],circCentMap.size());
      if (cellIts[0] != circCentMap.end() && cellIts[1] != circCentMap.end())
	continue;

      assert(
	  cellIts[0] == circCentMap.end() || cellIts[1] == circCentMap.end());
      // index in cells[] of which cell is new (hasn't been checked)
      int newCellIndex = (cellIts[0] == circCentMap.end()) ? 0 : 1;
      // for each cell, if we haven't seen it, add it
      // this fills adCircs[] with the voronoi edge
      for (int i = 0; i < 2; i++) {
	if (newCellIndex == i) {
	  // we haven't seen this point
	  double adTempCirc[2];
	  cells[i]->calcCircumcenter(adTempCirc);
	  circCentMap.insert(std::make_pair(cells[i], adTempCirc));
	  adCircs[i][0] = adTempCirc[0];
	  adCircs[i][1] = adTempCirc[1];
	}
	else {
	  adCircs[i][0] = cellIts[i]->second[0];
	  adCircs[i][1] = cellIts[i]->second[1];
	}
      }
      //		printf("circ centers %f %f | %f %f\n",adCircs[0][0],adCircs[0][1],adCircs[1][0],adCircs[1][1]);
      //		writeVTK_SingleCell(cells[0],"tempcell0.vtk");
      //		writeVTK_SingleCell(cells[1],"tempcell1.vtk");
      //		printf("%f %f %f hmm\n",adPetalCirc[0],adPetalCirc[1],sqrt(m_petalRad2*dAB2));

      //		printf("circ centers %.16f %.16f | %.16f %.16f\n",adCircs[0][0],adCircs[0][1],adCircs[1][0],adCircs[1][1]);
      //				printf("face %p | cells %p %p\n",face,cells[0],cells[1]);
      bool qIntersect = calcSegmentCircleIntersections(adCircs[0], adCircs[1],
						       adPetalCirc,
						       m_petalRad2 * dAB2,
						       adIntersect);
      //		printf("intersect %d?\n",qIntersect);
      if (qIntersect == true) {
	// we intersect, so get the distance from a vertex to the intersection with the petal
	distToCheck2 = dDIST_SQ_2D(face->getVert(0)->getCoords(), adIntersect);

	//			assert(iFuzzyComp(dDIST_SQ_2D(face->getVert(0)->getCoords(),adIntersect),
	//					dDIST_SQ_2D(face->getVert(1)->getCoords(),adIntersect)) == 0);
//			ptsChecked.push_back(adIntersect);
//			ptTypes.push_back(SINK);
//			dChecked.push_back(distToCheck2);
	// if its further away , the intersection is "winning"
	if (distToCheck2 > furthestDist2) {
	  adPointLoc[0] = adIntersect[0];
	  adPointLoc[1] = adIntersect[1];
	  furthestDist2 = distToCheck2;

	  pointType = eSink;
	}
	// this branch of the voronoi diagram is done, since following it takes us outside the petal
      }
      else {
	//			printf("circ centers %.16f %.16f | %.16f %.16f\n",adCircs[0][0],adCircs[0][1],adCircs[1][0],adCircs[1][1]);
	assert(
	    iFuzzyComp(dDIST_SQ_2D(adCircs[newCellIndex], adPetalCirc), (m_petalRad2 * dAB2)) < 1);

	// we've hit a circumcenter inside the petal but not from our initial cell
	// get the circumradius of the unchecked cells' circumcircle
	distToCheck2 = dDIST_SQ_2D(face->getVert(0)->getCoords(),
				   adCircs[newCellIndex]);
//			ptsChecked.push_back(adCircs[newCellIndex]);
//			ptTypes.push_back(OTHERCIRCUMCENTER);
//			dChecked.push_back(distToCheck2);
	//			cell->circumscribesPoint(adCircs[newCellIndex]);
	// if its further away, its "winning"
	if (distToCheck2 > furthestDist2) {
	  adPointLoc[0] = adCircs[newCellIndex][0];
	  adPointLoc[1] = adCircs[newCellIndex][1];
	  furthestDist2 = distToCheck2;

	  pointType = eOthercircumcenter;
	}
	// since we are still inside the petal, add neighboring faces
	// that aren't the face we started with
	for (int iFace = 0; iFace < 3; iFace++) {
	  const Face * newFace = cells[newCellIndex]->getFace(iFace);
	  if (newFace != face)
	    facesToTest.push(newFace);
	}
      }

    } // if empty, we are done, with the location stored in adPointLoc

//	printf("Points checked %d\n",ptsChecked.size());
//	for(int iP = 0; iP < min(int(ptsChecked.size()),10); iP++){
//		printf("Point seen %f %f %d | %f\n",ptsChecked[iP][0],ptsChecked[iP][1],ptTypes[iP], dChecked[iP]);
//	}
//
//	printf("Point chosen %f %f %d  out of %d | %f\n",adPointLoc[0],adPointLoc[1],pointType,ptsChecked.size(), furthestDist2);
//

    return pointType;
  }

  double
  OffcenterEngwirda2D::FindAltitude(const Vert* pVA, const Vert* pVB,
				    double circumcenter[],
				    const Vert* pVC) const
  {
    /// Find the position of inserted point
    /// altitude
//	double LSPA = pVA->getLengthScale()
//			+ 0.5 * (pVA->getLengthScale() + pVB->getLengthScale()) / m_g;
//	//double LSPA = pVA->getLengthScale() + pVA->getLengthScale()/m_g;
//	double LSPB = pVB->getLengthScale()
//			+ 0.5 * (pVA->getLengthScale() + pVB->getLengthScale()) / m_g;
//	//double LSPB = pVB->getLengthScale() + pVB->getLengthScale()/m_g;
//	LSPA = 0.5 * (LSPA + pVA->getLengthScale());
//	//LSPA = min(0.5*(LSPA + pVA->getLengthScale()),1.1*cell->calcShortestEdgeLength());
//	LSPB = 0.5 * (LSPB + pVB->getLengthScale());
//	//LSPB = min(0.5*(LSPB + pVB->getLengthScale()),1.1*cell->calcShortestEdgeLength());
//		const double *const adA = pVA->getCoords();
//		const double *const adB = pVB->getCoords();
//	double dAB2 = dDIST_SQ_2D(adA,adB);
//	double alti = sqrt(0.5 * (LSPA + LSPB) * 0.5 * (LSPA + LSPB) - 0.25 * dAB2);

    //////////////////////////////////////////////
    //bisection method
    const double * const adA = pVA->getCoords();
    const double * const adB = pVB->getCoords();
    const double * const adC = pVC->getCoords();

    double midAB[2] =
      { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };
    double h0 = 2.
	* sqrt(pow(adC[0] - midAB[0], 2) + pow(adC[1] - midAB[1], 2));
    double h1 = 0.;
//	printf("circumcenter: %6.4f\t%6.4f\n",circumcenter[0],circumcenter[1]);
    double t = 0.;
    double location[2] =
      { 0., 0. };
    double distA = 0.;
    double distC = 0.;
    double lsLoc = 0.;
    double answer = 1.e5;
    int b0, b1;

    t = (h0) / sqrt(dDIST_SQ_2D(circumcenter, midAB));
    location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
    location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
    distA = sqrt(dDIST_SQ_2D(location, adA));
    distC = sqrt(dDIST_SQ_2D(location, adC));
    lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
		     pVA->getLengthScale() + distA / m_g);
    lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
    answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
    b0 = (answer > 0) ? 1 : -1;
    t = (h1) / sqrt(dDIST_SQ_2D(circumcenter, midAB));
    location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
    location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
    distA = sqrt(dDIST_SQ_2D(location, adA));
    distC = sqrt(dDIST_SQ_2D(location, adC));
    lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
		     pVA->getLengthScale() + distA / m_g);
    lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
    answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
    b1 = (answer > 0) ? 1 : -1;
    if (b0 * b1 > 0) {
      return sqrt(dDIST_SQ_2D(circumcenter, midAB));
    }
    assert(b0 * b1 < 0);
    double hc = 0.5 * (h0 + h1);
    while (true) {
      hc = 0.5 * (h0 + h1);
      t = hc / sqrt(dDIST_SQ_2D(circumcenter, midAB));
      location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
      location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
      distA = sqrt(dDIST_SQ_2D(location, adA));
      distC = sqrt(dDIST_SQ_2D(location, adC));
      lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
		       pVA->getLengthScale() + distA / m_g);
      lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
      answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
      if (fabs(answer) < 1e-4 * distA)
	break;
      else {
	int bc = (answer > 0) ? 1 : -1;
	if (bc * b0 < 0) {
	  h1 = hc;
	  b1 = bc;
	}
	else {
	  h0 = hc;
	  b0 = bc;
	}
      }

    }

    double alti = t * sqrt(dDIST_SQ_2D(circumcenter, midAB));
//
////	bool B0 = iFuzzyComp(start0,)
//	double alti;
    return alti;
  }
#ifdef HAVE_MESQUITE
  InsertionPointCalculator::ePointType
  AnisoEdgeIPC2D::calcInsertionPoint(const Face * const pF,
				     double adPointLoc[]) const
  {

    Cell*pCL = pF->getLeftCell();
    Cell*pCR = pF->getRightCell();

    const Vert *pVA = pF->getVert(0);
    const Vert *pVB = pF->getVert(1);
    const Vert *pVC = pCL->getOppositeVert(pF);
    const Vert *pVD = pCR->getOppositeVert(pF);
    assert(pCL->isValid());
    assert(pCR->isValid());

    double dCoordsA[2];
    double dCoordsB[2];
    dCoordsA[0] = pVA->x();
    dCoordsA[1] = pVA->y();
    dCoordsB[0] = pVB->x();
    dCoordsB[1] = pVB->y();

    double QualBef[2];
    double QualAft[4];
    double dT[9];
    double adMetric[3];
    double dOldQual = 1e10;
    double dNewQual;
    double dFrate = 0.5;
    bool bContinue = false;
    GRUMMP::TMOPQual::eQMetric QM;
    GRUMMP::TMOPQual TQ(2);

    for (int i = 0; i < 3; i++)
      adMetric[i] = (pVA->getMetric(i) + pVB->getMetric(i)) / 2.;

    QM = Aniso_QM;

    if (adMetric[1] < dMetricTol && adMetric[1] > -dMetricTol) {
      if (fabs(adMetric[0] - adMetric[2]) < dMetricTol) {
	QM = Iso_QM;
      }
    }

    adPointLoc[0] = (pVA->x() + (pVB->x() - pVA->x()) * dFrate);
    adPointLoc[1] = (pVA->y() + (pVB->y() - pVA->y()) * dFrate);

    do {

      double dMlenghtA = 0;
      double dMlenghtB = 0;
      double dDelM;
      double dQualMaxAft;
      double adIP[2];

      for (int i = 0; i < 3; i++)
	adMetric[i] = pVA->getMetric(i);
      dMlenghtA = TQ.dLogMetricLength(dCoordsA, adMetric, dCoordsB, adMetric);

      for (int i = 0; i < 3; i++)
	adMetric[i] = pVB->getMetric(i);
      dMlenghtB = TQ.dLogMetricLength(dCoordsA, adMetric, dCoordsB, adMetric);

      dDelM = dMlenghtA / (dMlenghtA + dMlenghtB);
      adIP[0] = (pVA->x() + (pVB->x() - pVA->x()) * dFrate);
      adIP[1] = (pVA->y() + (pVB->y() - pVA->y()) * dFrate);

      for (int i = 0; i < 3; i++)
	adMetric[i] = (pVA->getMetric(i) + pVB->getMetric(i)) / 2.;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      QualBef[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			    pVB->getCoords(), pVC->getCoords());

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA,
				       pVC);
      QualAft[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			    adPointLoc, pVC->getCoords());

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVB,
				       pVC);
      QualAft[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, adPointLoc,
			    pVB->getCoords(), pVC->getCoords());

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVD, pVB);
      QualBef[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			    pVD->getCoords(), pVB->getCoords());

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA,
				       pVD);
      QualAft[2] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			    pVD->getCoords(), adPointLoc);

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVD,
				       pVB);
      QualAft[3] = TQ.dEval(QM, CellSkel::eTriCell, dT, adPointLoc,
			    pVD->getCoords(), pVB->getCoords());

      dQualMaxAft = std::max(std::max(QualAft[0], QualAft[1]),
			     std::max(QualAft[2], QualAft[3]));
      dNewQual =
	  dQualMaxAft < std::max(QualBef[0], QualBef[1]) ? dQualMaxAft : -1;

      double dML = 1e10;
      double dMinML = 1e10;
      double adM[3];
      double adC[2];
      const Vert *apV[4] =
	{ pVA, pVB, pVC, pVD };

      for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 3; j++)
	  adM[j] = apV[i]->getMetric(j);

	adC[0] = apV[i]->x();
	adC[1] = apV[i]->y();

	dML = TQ.dLogMetricLength(adPointLoc, adMetric, adC, adM);
	dMinML = std::min(dML, dMinML);
      }
      dNewQual = dMinML < dMinLength ? -1 : dNewQual;

      if (dNewQual < dOldQual && dNewQual > 0) {
	adPointLoc[0] = adIP[0];
	adPointLoc[1] = adIP[1];
	dOldQual = dNewQual;
	bContinue = true;
	if (dFrate > dDelM)
	  dFrate = (dFrate - 0.05) < dDelM ? dDelM : dFrate - 0.05;
	if (dFrate < dDelM)
	  dFrate = (dFrate + 0.05) > dDelM ? dDelM : dFrate + 0.05;
      }
      else
	bContinue = false;

    }
    while (bContinue == true);

    return eAnisoEdge;
  }

  InsertionPointCalculator::ePointType
  AnisoEdgeIPC2D::calcInsertionPoint(const BFace * const pBF,
				     double adPointLoc[]) const
  {

    const Vert *pVA = pBF->getVert(0);
    const Vert *pVB = pBF->getVert(1);
    const Vert *pVC;
    const Face *pF = pBF->getFace(0);
    if (pF->getRightCell()->isBdryCell())
      pVC = pF->getLeftCell()->getOppositeVert(pF);

    else {
      pVC = pF->getRightCell()->getOppositeVert(pF);
      assert(pF->getLeftCell()->isBdryCell());
    }

    assert(pVA->isValid());
    assert(pVB->isValid());
    assert(pVC->isValid());

    double dCoordsA[2];
    double dCoordsB[2];
    dCoordsA[0] = pVA->x();
    dCoordsA[1] = pVA->y();
    dCoordsB[0] = pVB->x();
    dCoordsB[1] = pVB->y();

    double QualBef;
    double QualAft[2];
    double dT[9];
    double adMetric[3];
    double dOldQual = 1e10;
    double dNewQual;
    double dFrate = 0.5;
    bool bContinue = false;
    GRUMMP::TMOPQual::eQMetric QM;
    GRUMMP::TMOPQual TQ(2);

    for (int i = 0; i < 3; i++)
      adMetric[i] = (pVA->getMetric(i) + pVB->getMetric(i)) / 2.;

    QM = Aniso_QM;

    if (adMetric[1] < dMetricTol && adMetric[1] > -dMetricTol) {
      if (fabs(adMetric[0] - adMetric[2]) < dMetricTol) {
	QM = Iso_QM;
      }
    }

    adPointLoc[0] = (pVA->x() + (pVB->x() - pVA->x()) * dFrate);
    adPointLoc[1] = (pVA->y() + (pVB->y() - pVA->y()) * dFrate);

    do {

      double dMlenghtA = 0;
      double dMlenghtB = 0;
      double dDelM;
      double dQualMaxAft;
      double adIP[2];

      for (int i = 0; i < 3; i++)
	adMetric[i] = pVA->getMetric(i);
      dMlenghtA = TQ.dLogMetricLength(dCoordsA, adMetric, dCoordsB, adMetric);

      for (int i = 0; i < 3; i++)
	adMetric[i] = pVB->getMetric(i);
      dMlenghtB = TQ.dLogMetricLength(dCoordsA, adMetric, dCoordsB, adMetric);

      dDelM = dMlenghtA / (dMlenghtA + dMlenghtB);

      adIP[0] = (pVA->x() + (pVB->x() - pVA->x()) * dFrate);
      adIP[1] = (pVA->y() + (pVB->y() - pVA->y()) * dFrate);

      for (int i = 0; i < 3; i++)
	adMetric[i] = (pVA->getMetric(i) + pVB->getMetric(i)) / 2.;

      TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
      QualBef = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			 pVB->getCoords(), pVC->getCoords());

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA,
				       pVC);
      QualAft[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			    adPointLoc, pVC->getCoords());

      TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVB,
				       pVC);
      QualAft[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, adPointLoc,
			    pVB->getCoords(), pVC->getCoords());

      dQualMaxAft = std::max(QualAft[0], QualAft[1]);
      dNewQual = dQualMaxAft < QualBef ? dQualMaxAft : -1;

      double dML = 1e10;
      double dMinML = 1e10;
      double adM[3];
      double adC[2];
      const Vert *apV[3] =
	{ pVA, pVB, pVC };

      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++)
	  adM[j] = apV[i]->getMetric(j);

	adC[0] = apV[i]->x();
	adC[1] = apV[i]->y();

	dML = TQ.dLogMetricLength(adPointLoc, adMetric, adC, adM);
	dMinML = std::min(dML, dMinML);
      }
      dNewQual = dMinML < dMinLength ? -1 : dNewQual;

      if (dNewQual < dOldQual && dNewQual > 0) {
	adPointLoc[0] = adIP[0];
	adPointLoc[1] = adIP[1];
	dOldQual = dNewQual;
	bContinue = true;
	if (dFrate > dDelM)
	  dFrate = (dFrate - 0.05) < dDelM ? dDelM : dFrate - 0.05;
	if (dFrate < dDelM)
	  dFrate = (dFrate + 0.05) > dDelM ? dDelM : dFrate + 0.05;
      }
      else
	bContinue = false;

    }
    while (bContinue == true);

    return eAnisoEdge;
  }

  InsertionPointCalculator::ePointType
  AnisoEdgeIPC2D::calcInsertionPoint(const SimplexCell * const cell,
				     double adPointLoc[]) const
  {
    cell->isValid();
    adPointLoc[0] = 0;
    adPointLoc[1] = 0;
    assert(0);
    return eAnisoEdge;
  }

  InsertionPointCalculator::ePointType
  AnisoEdgeIPC3D::calcInsertionPoint(const SimplexCell * const cell,
				     double adPointLoc[]) const
  {
    cell->isValid();
    adPointLoc[0] = 0;
    adPointLoc[1] = 0;
    assert(0);
    return eAnisoEdge;
  }

  InsertionPointCalculator::ePointType
  AnisoEdgeIPC3D::calcInsertionPoint(const Vert* const pV0,
				     const Vert* const pV1,
				     double adPointLoc[]) const
  {

    assert(pV0->getSpaceDimen() == 3 && pV1->getSpaceDimen() == 3);
    double adCoord1[3];
    double adCoord2[3];

    adCoord1[0] = pV0->x();
    adCoord1[1] = pV0->y();
    adCoord1[2] = pV0->z();

    adCoord2[0] = pV1->x();
    adCoord2[1] = pV1->y();
    adCoord2[2] = pV1->z();

    adPointLoc[0] = (adCoord1[0] + adCoord2[0]) / 2.;
    adPointLoc[1] = (adCoord1[1] + adCoord2[1]) / 2.;
    adPointLoc[2] = (adCoord1[2] + adCoord2[2]) / 2.;

    return eAnisoEdge;
  }
#endif

///Offcenter calculator in Engwirda thesis
  InsertionPointCalculator::ePointType
  OffcenterEngwirda2D::calcInsertionPoint(const SimplexCell * const cell,
					  double adPointLoc[]) const
  {
    assert(cell->getType() == Cell::eTriCell);
    ///adPointLoc here is the Circumcenter of cell
    cell->calcCircumcenter(adPointLoc);
    ///if the minimum angle of the cell is not less than the threshold angle.
    ///insert Circumcenter
    if (iFuzzyComp(cell->calcMinSolid() * M_PI / 180., m_alpha) > -1) {
      //logMessage(1,"eCircumcenter\n");
      return eCircumcenter;
    }
    ///Find the shortest edge of the cell
    const Face * pFShortest = cell->getFace(0);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(1)->calcSize() ?
	    pFShortest : cell->getFace(1);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(2)->calcSize() ?
	    pFShortest : cell->getFace(2);

    // lets get the shortest edge, and preserve orientation
    const Vert * pVA = pFShortest->getVert(0);
    const Vert * pVB = pFShortest->getVert(1);
    const Vert * pVC = cell->getOppositeVert(pFShortest);

    assert(iFuzzyComp(dDist(pVA, pVB), cell->calcShortestEdgeLength()) == 0);

    /// Coordinates of the end points of the all three edges
    /// The length square of the shortest edge.
    const double * const adA = pVA->getCoords();
    const double * const adB = pVB->getCoords();
//	const double *const adC = pVC->getCoords();
    double dAB2 = dDIST_SQ_2D(adA, adB);
    /// The mid point of the shortest edge
    double adAB[2] =
      { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };

    ///
    double adShapeOptimal[2];
    double adSizeOptimal[2];
    /// Find the position of inserted point
    /// altitude
    double alti1 = FindAltitude(pVA, pVB, adPointLoc, pVC);
    double alti2 = FindAltitude(pVB, pVA, adPointLoc, pVC);
    double alti = std::min(alti1, alti2);
//	printf("alti: %6.4f\n",alti);
    double t = m_coeff * sqrt(dAB2 / dDIST_SQ_2D(adPointLoc, adAB));
    //can't lower than
    double tlow = sqrt(dAB2) / 2. * tan(m_alpha)
	/ sqrt(dDIST_SQ_2D(adPointLoc, adAB));
    t = t > tlow ? t : tlow;

    adShapeOptimal[0] = adAB[0] + t * (adPointLoc[0] - adAB[0]);
    adShapeOptimal[1] = adAB[1] + t * (adPointLoc[1] - adAB[1]);

    double tsize = alti / sqrt(dDIST_SQ_2D(adPointLoc, adAB));
    adSizeOptimal[0] = adAB[0] + tsize * (adPointLoc[0] - adAB[0]);
    adSizeOptimal[1] = adAB[1] + tsize * (adPointLoc[1] - adAB[1]);

    /// find the point minimises the distance to the edge midpoint
    double dmidShape = dDIST_SQ_2D(adShapeOptimal, adAB);
    double dmidSize = dDIST_SQ_2D(adSizeOptimal, adAB);

    double adOptimal[2];
    if (iFuzzyComp(dmidShape, dmidSize) > -1) {
      adOptimal[0] = adSizeOptimal[0];
      adOptimal[1] = adSizeOptimal[1];
    }
    else {
      adOptimal[0] = adShapeOptimal[0];
      adOptimal[1] = adShapeOptimal[1];
    }
    /// Make sure the offcenter is on the dual of the shortest edge
    if (iFuzzyComp(dDIST_SQ_2D(adAB, adOptimal), dDIST_SQ_2D(adAB, adPointLoc))
	> -1) {
      //logMessage(1,"eCircumcenter\n");
      return eCircumcenter;
    }
    else {
      adPointLoc[0] = adOptimal[0];
      adPointLoc[1] = adOptimal[1];

      if (iFuzzyComp(dmidShape, dmidSize) > -1) {
	//logMessage(1,"eSizeOptimal\n");
	return eSizeOptimal;
      }
      else {
	//logMessage(1,"eShapeOptimal\n");
	return eShapeOptimal;
      }

    }
  }
////////////
//////////////
/////////////////
//Advancing Frontal
  InsertionPointCalculator::ePointType
  AF2D::calcInsertionPoint(const SimplexCell * const cell,
			   double adPointLoc[]) const
  {
    assert(cell->getType() == Cell::eTriCell);
    cell->calcCircumcenter(adPointLoc);
    //find the shortest edge
    const Face * pFShortest = cell->getFace(0);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(1)->calcSize() ?
	    pFShortest : cell->getFace(1);
    pFShortest =
	pFShortest->calcSize() < cell->getFace(2)->calcSize() ?
	    pFShortest : cell->getFace(2);
    //Choose shortest edge as base edge
    //find all three vertices
    if (cell->getFace(0)->isBdryFace()) {
      pFShortest = cell->getFace(0);
    }
    if (cell->getFace(1)->isBdryFace()
	and (cell->getFace(1)->calcSize() < pFShortest->calcSize())) {
      pFShortest = cell->getFace(1);
    }
    if (cell->getFace(2)->isBdryFace()
	and (cell->getFace(2)->calcSize() < pFShortest->calcSize())) {
      pFShortest = cell->getFace(2);
    }
    const Vert * pVA = pFShortest->getVert(0);
    const Vert * pVB = pFShortest->getVert(1);
//	const Vert * pVC = cell->getOppositeVert(pFShortest);
//	assert(iFuzzyComp(dDist(pVA,pVB), cell->calcShortestEdgeLength()) == 0);

//	const Face * indexShortest = cell->getFace(0);
//	indexShortest = indexShortest->getOppositeCell(cell)->getIndex() < cell->getFace(1)->getOppositeCell(cell)->getIndex() ? indexShortest : cell->getFace(1);
//	indexShortest = indexShortest->getOppositeCell(cell)->getIndex() < cell->getFace(2)->getOppositeCell(cell)->getIndex() ? indexShortest : cell->getFace(2);
//	if (cell->getFace(0)->isBdryFace()){
//		indexShortest = cell->getFace(0);
//	}
//	if (cell->getFace(1)->isBdryFace()){
//		indexShortest = cell->getFace(1);
//	}
//	if (cell->getFace(2)->isBdryFace()){
//		indexShortest = cell->getFace(2);
//	}
//	const Vert * pVA = indexShortest->getVert(0);
//	const Vert * pVB = indexShortest->getVert(1);
//	const Vert * pVC = cell->getOppositeVert(indexShortest);

    /// Coordinates of the end points of the all three edges
    /// The length square of the shortest edge.
    const double * const adA = pVA->getCoords();
    const double * const adB = pVB->getCoords();
    /// The mid point of the shortest edge
    double adAB[2] =
      { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };
    // find the candidate point
    double t = dDist(pVA, pVB) / 2. * sqrt(3.)
	/ sqrt(dDIST_SQ_2D(adPointLoc, adAB));

    //find the length scale of new points
//	double LSPA = pVA->getLengthScale() + dDist(pVA,pVB)/2.*sqrt(3.)/m_g;
    //double LSPA = pVA->getLengthScale() + pVA->getLengthScale()/m_g;
//	double LSPB = pVB->getLengthScale() + dDist(pVA,pVB)/2.*sqrt(3.)/m_g;
//	double reflength = (pVA->getLengthScale()+pVB->getLengthScale()+(LSPA+LSPB)/2.)/3.;
//	double reflength = (pVA->getLengthScale()+pVB->getLengthScale()+pVC->getLengthScale())/3.;
//	if(iFuzzyComp(dDist(pVA,pVB),pVC->getLengthScale()) > - 1.){
//		t = t*(1. - 1./m_g);
//	}
//	else
    t = t * (1. + 1. / m_g);
    adPointLoc[0] = adAB[0] + t * (adPointLoc[0] - adAB[0]);
    adPointLoc[1] = adAB[1] + t * (adPointLoc[1] - adAB[1]);
    return eAF;

  }

} //end namespace GRUMMP
