#include "GR_CurveDiscretizer.h"
#include "GR_GeomComp.h"
#include "GR_BaryCentricLengthScale.h"
#include "GR_misc.h"

#include <list>
#include <set>

using std::list;

void
CurveSample::init_sample(unsigned int num_points, bool compute_tangent)
{

  clear_sample();
  assert(m_curve_sample.empty());

  if (num_points < 3)
    num_points = 3;
  unsigned int num_sample = num_points - 2; //num_points - number of end points

  double min_param, max_param, param_diff, fraction;
  m_curve->get_param_range(min_param, max_param);
  param_diff = max_param - min_param;
  assert(param_diff > 0.);

  CurveSampleData csd_beg, csd_end;
  CubitVector coord, tangent;
  csd_beg.set_param(min_param);
  csd_end.set_param(max_param);
  csd_beg.set_coord(m_curve->start_coordinates());
  csd_end.set_coord(m_curve->end_coordinates());

  if (compute_tangent) {

    m_curve->position_from_u(min_param + 0.0001, coord);
    m_curve->tangent(coord, tangent);
    csd_beg.set_tangent(tangent);

    m_curve->position_from_u(max_param - 0.0001, coord);
    m_curve->tangent(coord, tangent);
    csd_end.set_tangent(tangent);

  }

  m_curve_sample.push_back(csd_beg);

  for (unsigned int i = 0; i < num_sample; i++) {

    fraction = static_cast<double>(i + 1) / static_cast<double>(num_sample + 1);

    csd_beg.set_param(min_param + (param_diff * fraction));
    m_curve->position_from_u(csd_beg.get_param(), coord);
    csd_beg.set_coord(coord);

    if (compute_tangent) {
      m_curve->tangent(*csd_beg.get_coord(), tangent);
      csd_beg.set_tangent(tangent);
    }

    m_curve_sample.push_back(csd_beg);

  }

  m_curve_sample.push_back(csd_end);

  assert(m_curve_sample.size() >= 3);
  assert(m_curve_sample.size() == num_points);

}

void
CurveSample::sample_tvt(double threshold)
{

  assert(m_curve);
  assert(iFuzzyComp(threshold, 0.) == 1);
  assert(iFuzzyComp(threshold, 0.1) == -1);

  init_sample(10, true);

  double mid_param, tvt1, tvt2, tvt_diff;
  CubitVector mid_coord, mid_tangent;
  CurveSampleData sample_data;

  list<CurveSampleData>::iterator it_prev, it;
  it_prev = it = m_curve_sample.begin();

  // Klocwork complains about this iterator being dereferenced when it
  // can be equal to end().  But it's never == end() at loop entry and
  // never changed in the loop until the very last line, so I disagree
  // with Klocwork's analysis here.  Unless, maybe, the analysis is
  // looking at cases where there are initially zero entries in
  // m_curve_sample?
  while (++it != m_curve_sample.end()) {

    assert(it != it_prev);
    assert(it_prev->get_tangent() && it->get_tangent());
    assert(it_prev->get_param() < it->get_param());

    mid_param = 0.5 * (it_prev->get_param() + it->get_param());
    m_curve->position_from_u(mid_param, mid_coord);
    m_curve->tangent(mid_coord, mid_tangent);

    tvt1 = GR_acos(*it_prev->get_tangent() % *it->get_tangent());
    tvt2 = GR_acos(*it_prev->get_tangent() % mid_tangent)
	+ GR_acos(*it->get_tangent() % mid_tangent);

    tvt_diff = fabs(tvt2 - tvt1);

    if (tvt_diff < threshold) {
      it_prev = it;
    }
    else {
      sample_data.set_param(mid_param);
      sample_data.set_coord(mid_coord);
      sample_data.set_tangent(mid_tangent);
      m_curve_sample.insert(it, sample_data);
      it = it_prev;
    }

  }

}

double
CurveSample::compute_tvt(double beg_param, double end_param)
{

  assert(!m_curve_sample.empty());

  double min_param, max_param;
  m_curve->get_param_range(min_param, max_param);

  double u1, u2;
  if (beg_param < min_param)
    u1 = min_param;
  else
    u1 = beg_param;
  if (end_param > max_param)
    u2 = max_param;
  else
    u2 = end_param;

  assert(iFuzzyComp(u1, u2) == -1);
  assert(iFuzzyComp(u1, min_param) >= 0 && iFuzzyComp(u1, max_param) <= 0);
  assert(iFuzzyComp(u2, min_param) >= 0 && iFuzzyComp(u2, max_param) <= 0);

  CubitVector beg_tangent, end_tangent, beg_coord, end_coord;
  m_curve->position_from_u(u1, beg_coord);
  m_curve->position_from_u(u2, end_coord);
  m_curve->tangent(beg_coord, beg_tangent);
  m_curve->tangent(end_coord, end_tangent);

  list<CurveSampleData>::iterator it, it_beg = std::find_if(
      m_curve_sample.begin(), m_curve_sample.end(), ParamGreaterThan(u1)),
      it_end = --std::find_if(it_beg, m_curve_sample.end(),
			      ParamGreaterThan(u2));

  assert(iFuzzyComp(u1, it_beg->get_param()) != 1);
  assert(iFuzzyComp(u2, it_end->get_param()) != -1);
  assert(it_beg->get_tangent());
  assert(it_end->get_tangent());

  //This will happen if u1 and u2 are in the same sample interval.
  if (it_beg->get_param() > it_end->get_param())
    return GR_acos(beg_tangent % end_tangent);

  double tvt = 0.;
  CubitVector *tangent1 = it_beg->get_tangent(), *tangent2;

  tvt += GR_acos(beg_tangent % *tangent1);

  for (it = it_beg; it != it_end;) {
    tangent2 = (++it)->get_tangent();
    assert(tangent1 && tangent2);
    tvt += GR_acos(*tangent1 % *tangent2);
    tangent1 = tangent2;
  }

  tangent2 = it_end->get_tangent();
  tvt += GR_acos(*tangent2 % end_tangent);

  return tvt;

}

void
TangentDiscretizer::discretize(DiscreteCurve& discrete_curve)
{

  discrete_curve.clear();

  RefEdge* curve = m_sample.get_curve();
  double lo_param, hi_param, mid_param, tvt;
  CubitVector coord;

  curve->get_param_range(lo_param, hi_param);
  curve->position_from_u(lo_param, coord);
  discrete_curve.insert(CurvePoint(lo_param, coord));
  curve->position_from_u(hi_param, coord);
  discrete_curve.insert(CurvePoint(hi_param, coord));

  DiscreteCurve::iterator it_prev, it;
  it_prev = it = discrete_curve.begin();

  while (++it != discrete_curve.end()) {

    lo_param = it_prev->get_param();
    hi_param = it->get_param();
    assert(hi_param > lo_param);

    tvt = m_sample.compute_tvt(lo_param, hi_param);
    if (tvt > m_tvt_bound) {
      mid_param = 0.5 * (lo_param + hi_param);
      curve->position_from_u(mid_param, coord);
      it = discrete_curve.insert(CurvePoint(mid_param, coord)).first;
      it_prev = --it;
    }
    else {
      it_prev = it;
    }

  }

}

//void UniformDiscretizer::discretize(DiscreteCurve& discrete_curve) {
//
//	discrete_curve.clear();
//	assert(0);
//
//}
void
LengthScaleDiscretizer::discretize(DiscreteCurve& discrete_curve)
{
  discrete_curve.clear();
  // uses 2D version, in 3D
  RefEdge* curve = m_sample.get_curve();
  double loParam, hiParam, currentParam;
  double startLS, currentLS, endLS;
  CubitVector coord;
  double length = 0., currentLength = 0.;
  double change_in_arclen = LARGE_DBL, arclen_old;
  std::vector<double> arclength_at_new_points;
  std::vector<double> param_to_insert;

  length = curve->get_arc_length();
  curve->get_param_range(loParam, hiParam);

  // initialize everything
  curve->position_from_u(loParam, coord);
  discrete_curve.insert(CurvePoint(loParam, coord));
  startLS = m_bls->getLengthScale(coord);

  curve->position_from_u(hiParam, coord);
  discrete_curve.insert(CurvePoint(hiParam, coord));
  endLS = m_bls->getLengthScale(coord);

  param_to_insert.push_back(loParam);
  arclength_at_new_points.push_back(0.0);
  // provided we have one case
  if (iFuzzyComp(startLS, length) <= 0) {
    currentParam = curve->u_from_arc_length(loParam, startLS);
  }
  else {
    currentParam = hiParam;
  }

  currentLS = startLS;

  while (currentLength + currentLS < length - endLS * 0.5) {
    // we want the point at currentLS added to the list of points
    param_to_insert.push_back(currentParam);
    currentLength = curve->length_from_u(loParam, currentParam);
    arclength_at_new_points.push_back(currentLength);
    // get next potential point's location
    curve->position_from_u(currentParam, coord);
    currentLS = m_bls->getLengthScale(coord);
    // this is the next parameter
    currentParam = curve->u_from_arc_length(currentParam, currentLS);

  }
  arclength_at_new_points.push_back(length);
  param_to_insert.push_back(hiParam);
  if (param_to_insert.size() == 2) {
    // nothing to do if just endpoints
  }
  else {

    // initialize the length scales
    GR_index_t numPoints = static_cast<int>(param_to_insert.size());
    double * ls_at_halfpoint = new double[numPoints - 1];
    double * arclen_new = new double[numPoints];

    for (GR_index_t ii = 0; ii < numPoints; ii++)
      arclen_new[ii] = arclength_at_new_points[ii];

    for (GR_index_t ii = 0; ii < numPoints - 1; ii++) {
      curve->position_from_u(
	  curve->u_from_arc_length(
	      param_to_insert[ii],
	      0.5
		  * (arclength_at_new_points[ii + 1]
		      - arclength_at_new_points[ii])),
	  coord);
      ls_at_halfpoint[ii] = m_bls->getLengthScale(coord);
    }

    // lets smooth internal points. We know how many. lets figure out where to put them
    // uses Gauss-Seidel, fast enough right now
    int gs_iter = 0;
    while (change_in_arclen > 1.e-2 && gs_iter++ < 1000) {
      change_in_arclen = 0.0;
      for (GR_index_t ii = 1; ii < numPoints - 1; ii++) {
	arclen_old = arclen_new[ii];
	arclen_new[ii] = 1. / (ls_at_halfpoint[ii] + ls_at_halfpoint[ii - 1])
	    * (ls_at_halfpoint[ii - 1] * arclen_new[ii + 1]
		+ ls_at_halfpoint[ii] * arclen_new[ii - 1]);
	param_to_insert[ii] = curve->u_from_arc_length(
	    param_to_insert[ii - 1], arclen_new[ii] - arclen_new[ii - 1]);
	change_in_arclen += fabs(arclen_old - arclen_new[ii]) / arclen_old;
	curve->position_from_u(
	    curve->u_from_arc_length(
		param_to_insert[ii - 1],
		0.5 * (arclen_new[ii] - arclen_new[ii - 1])),
	    coord);
	ls_at_halfpoint[ii - 1] = m_bls->getLengthScale(coord);

      }
      arclength_at_new_points.clear();
      param_to_insert.clear();
    }

    for (GR_index_t ii = 1; ii < numPoints - 1; ii++) {
      currentParam = curve->u_from_arc_length(0.0, arclen_new[ii]);
      curve->position_from_u(currentParam, coord);
      discrete_curve.insert(CurvePoint(currentParam, coord));
    }
    delete[] ls_at_halfpoint;
    delete[] arclen_new;
  }
  // now we need to check if there are internal boundary intersections
  if (m_backgroundMesh->getNumIntBdryFaces() == 0)
    return;

  CurveDiscretizer::DiscreteCurve::const_iterator it = discrete_curve.begin();
  CurveDiscretizer::DiscreteCurve::const_iterator it_end =
      --discrete_curve.end();
  for (; it != it_end;) {
    CurvePoint cp1 = *it;
    CurvePoint cp2 = *(++it);
    CubitVector cv1(cp1.get_coord());
    CubitVector cv2(cp2.get_coord());
    double coords1[3];
    double coords2[3];
    cv1.get_xyz(coords1);
    cv2.get_xyz(coords2);
    GR_index_t iB = 0;
    bool qFound = false;
    CubitVector intersect;
    do {
      cv1.set(cp1.get_coord());
      cv2.set(cp2.get_coord());
      cv1.get_xyz(coords1);
      cv2.get_xyz(coords2);
      BFace * bface = m_backgroundMesh->getIntBFace(iB);
      int intersections = computeRefEdgeIntersectionsWithBFace(curve, coords1,
							       coords2, bface,
							       intersect);

      if (intersections == 1) {

	// presumably we have the coordinates now
	double dist1 = (cv1 - intersect).length();
	double dist2 = (cv2 - intersect).length();
	double new_param = curve->u_from_position(intersect);
	;
	CurvePoint newCurvePoint(new_param, intersect);
	DiscreteCurve::iterator newPoint =
	    discrete_curve.insert(newCurvePoint).first;

	// figure out which existing one is closer, remove it, and insert the other
	if (dist1 < dist2)
	  discrete_curve.erase(cp1);
	else
	  discrete_curve.erase(cp2);

	intBdryIntersections.insert(&(*newPoint));
	// we found one, so move on from that edge
	++it;
	qFound = true;
      } // found intersect
    }
    while (++iB < m_backgroundMesh->getNumIntBdryFaces() && qFound == false);
  }
}

void
QuasiStructDiscretizer::discretize(DiscreteCurve& discrete_curve)
{
  discrete_curve.clear();

  RefEdge* curve = m_sample.get_curve();
  CubitVector align(1.0, 0.0, 0.0);
  RefFace* surface;
  FacetSurface * facet_surface;
  DLIList<RefFace*> ref_faces;
  curve->ref_faces(ref_faces);
  surface = ref_faces.get_and_step();
  double tol = 1.e-6;
  DLIList<RefEdge*> rEdges;
  DLIList<RefVertex*> rVerts;
  facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());

  curve->ref_vertices(rVerts);
  CubitVector edge0, edge1, edge2;
  CubitVector point0, point1, point2;
  std::set<std::pair<Vert*, Vert*> > edges;

  double loParam, hiParam, param;
  CubitVector coord;
  curve->get_param_range(loParam, hiParam);

  // initialize everything
  curve->position_from_u(loParam, coord);
  discrete_curve.insert(CurvePoint(loParam, coord));

  curve->position_from_u(hiParam, coord);
  discrete_curve.insert(CurvePoint(hiParam, coord));

  std::vector<CubitVector> points;
  std::set<CubitVector, CubitVectorCompare> uniquePoints;
//

  for (GR_index_t iF = 0; iF < m_backgroundMesh->getNumBdryFaces(); iF++) {

    DLIList<CubitVector*> intersects;
    Face * face = m_backgroundMesh->getBFace(iF)->getFace(0);

    bool qE0 = edges.insert(
	std::make_pair<Vert*, Vert*>(face->getVert(1), face->getVert(2))).second;
    bool qE1 = edges.insert(
	std::make_pair<Vert*, Vert*>(face->getVert(0), face->getVert(2))).second;
    bool qE2 = edges.insert(
	std::make_pair<Vert*, Vert*>(face->getVert(0), face->getVert(1))).second;

    point0.set(face->getVert(0)->getCoords());
    point1.set(face->getVert(1)->getCoords());
    point2.set(face->getVert(2)->getCoords());

    edge0 = ~(point1 - point2);
    edge1 = ~(point0 - point2);
    edge2 = ~(point0 - point1);

    bool q0 = (face->getVert(1)->getVertType() != Vert::eInterior
	&& face->getVert(2)->getVertType() != Vert::eInterior);
    bool q1 = (face->getVert(0)->getVertType() != Vert::eInterior
	&& face->getVert(2)->getVertType() != Vert::eInterior);
    bool q2 = (face->getVert(0)->getVertType() != Vert::eInterior
	&& face->getVert(1)->getVertType() != Vert::eInterior);

    double dot0 = fabs(edge0 % align);
    double dot1 = fabs(edge1 % align);
    double dot2 = fabs(edge2 % align);

    //			double normal[3];
    //			face->calcNormal(normal);
    //			NORMALIZE3D(normal);
    //			printf("normal %f %f %f %d\n",iF,normal[0],normal[1],normal[2]);
    //			if(fabs(normal[0]) < 0.999999 && fabs(normal[1]) < 0.999999 && fabs(normal[2]) < 0.999999) continue;
    // find edge that is tangential to either y or z
    if (dot0 > (1. - tol) && qE0) {
      //			double edgelength[] = {(point1-point2).length(),(point0-point2).length(),(point1-point0).length()};
      //			if(edgelength[0] < edgelength[1] && edgelength[0] < edgelength[2]){
      facet_surface->get_eval_tool()->get_intersections(point1, point2,
							intersects, true);
      //				printf("intersects 0 size %d\n",intersects.size());
    }
    else if (dot1 > (1. - tol) && qE1) {
      //			if(edgelength[1] < edgelength[2] && edgelength[1] < edgelength[0]){
      facet_surface->get_eval_tool()->get_intersections(point0, point2,
							intersects, true);
      //				printf("intersects 1 size %d\n",intersects.size());
    }
    else if (dot2 > (1. - tol) && qE2) {
      //			if(edgelength[2] < edgelength[1] && edgelength[2] < edgelength[0]){
      facet_surface->get_eval_tool()->get_intersections(point0, point1,
							intersects, true);
      //				printf("intersects 2 size %d\n",intersects.size());
    }
    for (int iInt = 0; iInt < intersects.size(); iInt++) {
      CubitVector * cv = intersects.get_and_step();
      //				CubitVector newpoint;
      //				surface->find_closest_point_trimmed(CubitVector(*cv),newpoint);
      CubitVector point(*cv);

      if ((q0 && qE0) || (q1 && qE1) || (q2 && qE2))
	points.push_back(point);

    }
  }
  uniquePoints.insert(points.begin(), points.end());
  points.assign(uniquePoints.begin(), uniquePoints.end());
  uniquePoints.clear();
//	printf("points size %d\n",points.size());
//	for(int iRefVert = 0; iRefVert < rVerts.size(); iRefVert++){
//		RefVertex * refVert = rVerts.get_and_step();
//		printf("ref point %f %f %f\n", refVert->coordinates()[0],refVert->coordinates()[1],refVert->coordinates()[2]);
//
//	}
  for (GR_index_t iP = 0; iP < points.size(); iP++) {
    //		printf("qwe %f %f %f\n",points[iP][0],points[iP][1],points[iP][2]);
    bool qRefVert = false;
    for (int iRefVert = 0; iRefVert < rVerts.size(); iRefVert++) {
      RefVertex * refVert = rVerts.get_and_step();
      if ((refVert->coordinates() - points[iP]).length() < 1.e-6) {
	//				uniquePoints.insert(points[iP]);
	param = curve->u_from_position(points[iP]);
	discrete_curve.insert(CurvePoint(param, points[iP]));
	qRefVert = true;
	break;
      }
    }
    if (!qRefVert) {

      CubitVector closest;
      curve->closest_point_trimmed(points[iP], closest);
      if ((closest - points[iP]).length() < 1.e-3) {
	//					points.push_back(CubitVector(closest));
	//				printf("curve point %f %f %f | %f\n", points[iP][0],points[iP][1],points[iP][2],(closest-points[iP]).length() );
	//				uniquePoints.insert(points[iP]);
	param = curve->u_from_position(points[iP]);
	discrete_curve.insert(CurvePoint(param, points[iP]));
      }

    }

  }

  //	if(m_backgroundMesh->getNumIntBdryFaces() == 0) return;
  //
  //	CurveDiscretizer::DiscreteCurve::const_iterator it     =   discrete_curve.begin();
  //	CurveDiscretizer::DiscreteCurve::const_iterator it_end = --discrete_curve.end();
  //	for( ; it != it_end; ) {
  //		CurvePoint cp1 = *it;
  //		CurvePoint cp2 = *(++it);
  //		CubitVector cv1(cp1.get_coord());
  //		CubitVector cv2(cp2.get_coord());
  //		double coords1[3];
  //		double coords2[3];
  //		cv1.get_xyz(coords1);
  //		cv2.get_xyz(coords2);
  //		GR_index_t iB = 0;
  //		bool qFound = false;
  //		while(iB < m_backgroundMesh->getNumIntBdryFaces() && qFound == false){
  //			cv1.set(cp1.get_coord());
  //			cv2.set(cp2.get_coord());
  //			cv1.get_xyz(coords1);
  //			cv2.get_xyz(coords2);
  //			BFace * bface = m_backgroundMesh->getIntBFace(iB);
  //			iB++;
  //			if(checkOrient3D(bface->getVert(0)->getCoords(),
  //					bface->getVert(1)->getCoords(),bface->getVert(2)->getCoords(),coords1)*checkOrient3D(bface->getVert(0)->getCoords(),
  //							bface->getVert(1)->getCoords(),bface->getVert(2)->getCoords(),coords2) == 1) continue;
  //
  //			double bFaceNormal[3];
  //			bface->calcUnitNormal(bFaceNormal);
  //			NORMALIZE3D(bFaceNormal);
  //			CubitVector normalCV(bFaceNormal), bfaceCV(bface->getVert(0)->getCoords());
  //
  //			CubitVector tangent1,tangent2;
  //			double dist = (cv1-cv2).length();
  //			// basic approach is to walk along the curve in both directions until we have two coordinates very close to eachother
  //			// on either side of the internal boundary, and then check for intersection.
  //			int coordIter = 0;
  //			while(dist > 1.e-6 && ++coordIter < 100){
  //				// move the points closer before checking.
  //				curve->tangent(cv1,tangent1);
  //				curve->tangent(cv2,tangent2);
  //				tangent1.normalize();
  //				tangent2.normalize();
  //				// find intersection between tangent with plane
  //				double dist1 = ((bfaceCV-cv1) % normalCV)/(tangent1 % normalCV);
  //				double dist2 = ((bfaceCV-cv2) % normalCV)/(tangent2 % normalCV);
  //				cv1 = cv1 + dist1/2.*tangent1;
  //				cv2 = cv2 + dist2/2.*tangent2;
  //				// move the points slightly closer, still on the surface
  //				curve->closest_point_trimmed(cv1,cv1);
  //				curve->closest_point_trimmed(cv2,cv2);
  //				dist = (fabs(dist2)+fabs(dist1));
  //
  //			}
  //
  //			cv1.get_xyz(coords1);
  //			cv2.get_xyz(coords2);
  //			// bring the coords close together. Check if they intersect with the triangle
  //			double intersection[3];
  //			int intersect = calcSegmentTriangleIntersection(coords1,coords2,bface->getVert(0),
  //					bface->getVert(1),bface->getVert(2),0.,intersection);
  //
  //			if(intersect == 1){
  //
  //				double projection[3];
  //				CubitVector CVprojection, CVintersect;
  //				double dDiff = LARGE_DBL;
  //				int iter = 0;
  //				while(dDiff > 1e-10 && ++iter < 50){
  //					// do your projection stuff
  //					CVintersect.set(intersection);
  //					curve->closest_point_trimmed(CVintersect,CVprojection);
  //					CVprojection.get_xyz(projection);
  //					bface->findClosestPointOnBFace(projection,intersection);
  //					//					closest_on_triangle(projection,bface->getVert(0)->getCoords(),
  //					//							bface->getVert(1)->getCoords(),bface->getVert(2)->getCoords(),intersection);
  //					//					dDiff = dDIST3D(intersection,projection);
  //
  //				}
  //				// presumably we have the coordinates now
  //				//				double dist = (cv1-cv2).length();
  //				cv1.set(cp1.get_coord());
  //				cv2.set(cp2.get_coord());
  //				double dist1 = (cv1-CVprojection).length();
  //				double dist2 = (cv2-CVprojection).length();
  //				double new_param = curve->u_from_position(CVprojection);
  //				//				printf("CV %f %f %f\n",CVprojection.x(),CVprojection.y(),CVprojection.z());
  //				//								printf("testing this %f %f %f | %f %f %f\n",cv1.x(),cv1.y(),
  //				//										cv1.z(),cv2.x(),cv2.y(),
  //				//										cv2.z());
  //				CurvePoint newCurvePoint(new_param,CVprojection);
  //				if(dist1 < dist2){
  //					discrete_curve.erase(cp1);
  //					DiscreteCurve::iterator newPoint = discrete_curve.insert(newCurvePoint).first;
  //					intBdryIntersections.insert(&(*newPoint));
  //
  //				} else {
  //					discrete_curve.erase(cp2);
  //					DiscreteCurve::iterator newPoint = discrete_curve.insert(newCurvePoint).first;
  //					intBdryIntersections.insert(&(*newPoint));
  //				}
  //				++it;
  //				qFound = true;
  //			}
  //		}
  //	}
}
//void CurvatureDiscretizer::discretize(DiscreteCurve& discrete_curve) {
//
//	discrete_curve.clear();
//	assert(0);
//
//}
