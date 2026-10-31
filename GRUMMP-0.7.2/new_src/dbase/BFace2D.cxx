#include "GR_BFace2D.h"
#include "GR_Geometry.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_Vertex.h"
#include "GR_GRCurve.h"
#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "GeometryEntity.hpp"
#include <algorithm>

BdryEdgeBase::BdryEdgeBase(const bool qIsInternal, Face* const pFLeft,
			   Face* const pFRight, GRCurve* const parent_curve,
			   double param_vert0, double param_vert1) :
    BFace(qIsInternal, pFLeft, pFRight), m_curve(parent_curve), m_vert0_param(
	param_vert0), m_vert1_param(param_vert1), m_bdryPatch(nullptr)
{

  forbidOffsetInsertion();

}

BdryEdgeBase&
BdryEdgeBase::operator=(const BdryEdgeBase& BEB)
{

  if (this != &BEB) {

    this->BFace::operator=(BEB);
    m_bdryCond = BEB.m_bdryCond;
    m_curve = BEB.m_curve;
    m_vert0_param = BEB.m_vert0_param;
    m_vert1_param = BEB.m_vert1_param;

    m_bdryPatch = BEB.m_bdryPatch;

  }

  return *this;

}

bool
BdryEdgeBase::isCoplanarWith(const double pointLoc[]) const
{
  return (checkOrient2D(getVert(0)->getCoords(), getVert(1)->getCoords(),
			pointLoc) == 0);
}

void
BdryEdgeBase::setCurve(GRCurve* const curve)
{
  m_curve = curve;
  m_bdryCond = -1;
}

GRCurve*
BdryEdgeBase::getCurve() const
{
  assert(m_curve);
  return m_curve;
}

bool
BdryEdgeBase::doesCurveExist() const
{
  if (m_curve == NULL)
    return false;
  return true;
}

void
BdryEdgeBase::setVert0Param(double param)
{

  assert(m_curve);

  assert(
      iFuzzyComp(param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(param, m_curve->get_curve_geom()->max_param()) != 1);

  m_vert0_param = param;

}

void
BdryEdgeBase::setVert1Param(double param)
{

  assert(m_curve);

  assert(
      iFuzzyComp(param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(param, m_curve->get_curve_geom()->max_param()) != 1);

  m_vert1_param = param;

}

double
BdryEdgeBase::getVert0Param()
{
  assert(m_curve);

  assert(
      iFuzzyComp(m_vert0_param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(m_vert0_param, m_curve->get_curve_geom()->max_param())
	      != 1);
  assert(
      iFuzzyComp(m_vert1_param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(m_vert1_param, m_curve->get_curve_geom()->max_param())
	      != 1);

  int compare = iFuzzyComp(m_vert0_param, m_vert1_param);
  // Need to flip edge orientation if:
  //   v0_par > v1_par && isForward
  //   v0_par < v1_par && !isForward
  //
  // So in other words things are okay if exactly one of v0_par < v1_par (compare = -1)
  // and isForward() are true.  Use exclusive or.
  if (!((compare == 1) xor isForward())) {
    // Flip the orientation of the face.
    getFace()->interchangeCellsAndVerts();
  }
  assert(
      ((compare == -1) and isForward()) or ((compare == 1) and !isForward()));

  return m_vert0_param;
}

double
BdryEdgeBase::getVert1Param()
{
  assert(m_curve);

  assert(
      iFuzzyComp(m_vert0_param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(m_vert0_param, m_curve->get_curve_geom()->max_param())
	      != 1);
  assert(
      iFuzzyComp(m_vert1_param, m_curve->get_curve_geom()->min_param()) != -1
	  && iFuzzyComp(m_vert1_param, m_curve->get_curve_geom()->max_param())
	      != 1);

  int compare = iFuzzyComp(m_vert0_param, m_vert1_param);
  // Need to flip edge orientation if:
  //   v0_par > v1_par && isForward
  //   v0_par < v1_par && !isForward
  //
  // So in other words things are okay if exactly one of v0_par < v1_par (compare = -1)
  // and isForward() are true.  Use exclusive or.
  if (!((compare == 1) xor isForward())) {
    // Flip the orientation of the face.
    getFace()->interchangeCellsAndVerts();
  }
  assert(
      ((compare == -1) and isForward()) or ((compare == 1) and !isForward()));

  return m_vert1_param;
}

int
BdryEdgeBase::getBdryCondition(const bool left) const
{

  assert(isValid());

  if (left)
    return getLeftBdryCondition();
  else
    return getRightBdryCondition();

}

int
BdryEdgeBase::getLeftBdryCondition() const
{

  if (getNumFaces() == 2) {
    assert(dynamic_cast<IntBdryEdge*>(const_cast<BdryEdgeBase*>(this)));
    return iInvalidBC;
  }

  assert(getNumFaces() == 1);

  if (isForward())
    return m_curve->boundary_left();
  else
    return m_curve->boundary_right();

}

int
BdryEdgeBase::getRightBdryCondition() const
{

  if (getNumFaces() == 2) {
    assert(dynamic_cast<IntBdryEdge*>(const_cast<BdryEdgeBase*>(this)));
    return iInvalidBC;
  }

  assert(getNumFaces() == 1);

  if (isForward())
    return m_curve->boundary_right();
  else
    return m_curve->boundary_left();

}

double
BdryEdgeBase::calcLength() const
{

  const Vert *vert0 = getVert(0), *vert1 = getVert(1);

  double vector[] =
    { vert1->x() - vert0->x(), vert1->y() - vert0->y() };

  return dMAG2D(vector);

}

CubitVector
BdryEdgeBase::calcMidpoint() const
{

  const Vert *vert0 = getVert(0), *vert1 = getVert(1);

  return CubitVector(0.5 * (vert0->x() + vert1->x()),
		     0.5 * (vert0->y() + vert1->y()), 0.);

}

double
BdryEdgeBase::getSplitData(const SplitType split_type, CubitVector& split_coord,
			   bool& concentric_shell_split)
{

  const Vert *vert0 = getVert(0), *vert1 = getVert(1);

  if (m_curve) {
    bool vert0_small = vert0->isSmallAngleVert(), vert1_small =
	vert1->isSmallAngleVert(), vert0_shell = vert0->isShellVert(),
	vert1_shell = vert1->isShellVert();

    bool forward = isForward();

    GRCurveGeom * cg = m_curve->get_curve_geom();
    if (vert0_small && vert1_shell) {
      concentric_shell_split = true;
      if (forward) {
	double splitParam = cg->param_at_arc_length(getVert0Param(),
						    0.5 * calcLength());
	cg->coord_at_param(splitParam, split_coord);
	return splitParam;

      }
      else {
	double length = cg->arc_length(getVert1Param(), getVert0Param());
	double splitParam = cg->param_at_arc_length(
	    getVert1Param(), length - 0.5 * calcLength());
	cg->coord_at_param(splitParam, split_coord);
	return splitParam;

      }
    }

    if (vert1_small && vert0_shell) {

      concentric_shell_split = true;
      if (forward) {
	double length = cg->arc_length(getVert0Param(), getVert1Param());

	double splitParam = cg->param_at_arc_length(
	    getVert0Param(), length - 0.5 * calcLength());
	cg->coord_at_param(splitParam, split_coord);
	return splitParam;

      }
      else {
	double splitParam = cg->param_at_arc_length(getVert1Param(),
						    0.5 * calcLength());
	cg->coord_at_param(splitParam, split_coord);
	return splitParam;

      }
    }

    if (vert0_shell && vert1_shell) {

      double param_vert0 = forward ? getVert0Param() : getVert1Param();
      double param_vert1 = forward ? getVert1Param() : getVert0Param();

      concentric_shell_split = true;

      return m_curve->get_curve_geom()->coord_at_mid_dist(param_vert0,
							  param_vert1,
							  split_coord);
    }

    concentric_shell_split = false;

    double split_param;
    double lo_param = min(m_vert0_param, m_vert1_param); //forward ? getVert0Param() : getVert1Param();
    double hi_param = max(m_vert0_param, m_vert1_param); //forward ? getVert1Param() : getVert0Param();

    switch (split_type)
      {

      case MID_TVT:
	split_param = m_curve->get_curve_geom()->mid_TVT(lo_param, hi_param);
	break;

      case MID_PARAM:
	split_param = 0.5 * (lo_param + hi_param);
	break;

      default:
	vFatalError("Unknown split type specified",
		    "BdryEdgeBase::get_split_data()");
	break;

      }

    m_curve->get_curve_geom()->coord_at_param(split_param, split_coord);

    return split_param;
  }
  else {
    // Treat as a linear segment and split it.
    split_coord.set(0.5 * (getVert(0)->x() + getVert(1)->x()),
		    0.5 * (getVert(0)->y() + getVert(1)->y()), 0);
    return 0.5;
  }
}
CubitBox
BdryEdgeBase::bounding_box() const
{

  CubitVector mini, maxi;

  const Vert* vert0 = getVert(0);
  const Vert* vert1 = getVert(1);

  mini.set(std::min(vert0->x(), vert1->x()), std::min(vert0->y(), vert1->y()),
	   0.);
  maxi.set(std::max(vert0->x(), vert1->x()), std::max(vert0->y(), vert1->y()),
	   0.);

  return CubitBox(mini, maxi);

}

void
BdryEdgeBase::calcCentroid(double centroid[]) const
{

  assert(isValid());
  const Vert *vert0 = getVert(0), *vert1 = getVert(1);

  centroid[0] = 0.5 * (vert0->x() + vert1->x());
  centroid[1] = 0.5 * (vert0->y() + vert1->y());

}

void
BdryEdgeBase::calcCircumcenter(double circcenter[]) const
{

  calcCentroid(circcenter);

}

void
BdryEdgeBase::setGeometry(GeometryEntity* const geom_ent)
{

  m_curve = dynamic_cast<GRCurve*>(geom_ent);
  assert(m_curve);

}

GeometryEntity*
BdryEdgeBase::getGeometry()
{

  return getCurve();

}

eEncroachResult
BdryEdgeBase::isEncroached(const eEncroachType enc_type) const
{

  // Works for both regular and internal bdry edges.
  for (int iF = 0; iF < getNumFaces(); iF++) {
    const Face *pF = getFace(iF);
    Cell *pC = pF->getOppositeCell(this);
    // This cell could be a quad and we can still tell whether the bdry
    // edge is encroached on.
    if (!pC)
      return eClean;
    if (pC->getType() == Cell::eQuadCell)
      return eClean;
    for (int iV = 0; iV < pC->getNumVerts(); iV++) {
      Vert *pV = pC->getVert(iV);
      if (!hasVert(pV)) {
	enum eEncroachResult eER = isPointEncroaching(pV->getCoords(), enc_type,
						      false);
	if (eER != eClean)
	  return eER;
      }
    }
  }
  return eClean;

}

eEncroachResult
BdryEdgeBase::isPointEncroachingBall(const double coord[],
				     const bool tie_break) const
{

  const Vert *v1 = getVert(0), *v2 = getVert(1);
  double c1[] =
    { v1->x(), v1->y() };
  double c2[] =
    { v2->x(), v2->y() };

  double center[] =
    { 0.5 * (c1[0] + c2[0]), 0.5 * (c1[1] + c2[1]) };
  double radius = 0.5 * hypot(c2[0] - c1[0], c2[1] - c1[1]);
  double dist = hypot(coord[0] - center[0], coord[1] - center[1]);

  switch (iFuzzyComp(dist, radius))
    {
    case -1:
      return eSplit;
    case 1:
      return eClean;
    case 0:
      if (tie_break)
	return eSplit;
      else
	return eClean;
    default:
      assert(0);
      return eClean;
    }

}

eEncroachResult
BdryEdgeBase::isPointEncroachingLens(const double coord[],
				     const bool tie_break) const
{
  return isEncroachingLens(coord, tie_break, eLens);
}

eEncroachResult
BdryEdgeBase::isPointEncroachingNewLens(const double coord[],
					const bool tie_break) const
{
  return isEncroachingLens(coord, tie_break, eNewLens);
}

eEncroachResult
BdryEdgeBase::isEncroachingLens(const double coord[], const bool tie_break,
				const enum eEncroachType enc_type) const
{

  const Vert *v1 = getVert(0), *v2 = getVert(1);
  assert(v1->getSpaceDimen() == 2 && v2->getSpaceDimen() == 2);
  assert(enc_type == eLens || enc_type == eNewLens);

  double mid[] =
    { 0.5 * (v1->x() + v2->x()), 0.5 * (v1->y() + v2->y()) };
  double norm[] =
    { v1->y() - v2->y(), v2->x() - v1->x() };
  double mag = dMAG2D(norm);
  double rad = mag / sqrt(3.);
  norm[0] /= mag;
  norm[1] /= mag;

  double cent1[] =
    { mid[0] + 0.5 * rad * norm[0], mid[1] + 0.5 * rad * norm[1] };
  double cent2[] =
    { mid[0] - 0.5 * rad * norm[0], mid[1] - 0.5 * rad * norm[1] };

  int in1 = iFuzzyComp(dDIST2D(coord, cent1), rad);
  int in2 = iFuzzyComp(dDIST2D(coord, cent2), rad);

  //Both negative, definitely encroaching
  if (in1 == -1 && in2 == -1)
    return eSplit;
  // If one is negative and one zero, then use the tiebreaker
  else if ((in1 == -1 && in2 == 0) || (in1 == 0 && in2 == -1))
    return tie_break ? eSplit : eClean;
  else
    return eClean;

}

void
BdryEdgeBase::findClosestPointOnBFace(const double pt[3], double close[3])
{
  // Could easily return the distance, too, but the API doesn't
  // have a slot for that.
  if (doesCurveExist()) {
    GRCurve* curve = dynamic_cast<GRCurve*>(getCurve());
    assert(curve);
    double beg_param, end_param;
    if (isForward()) {
      beg_param = getVert0Param();
      end_param = getVert1Param();
    }
    else {
      beg_param = getVert1Param();
      end_param = getVert0Param();
    }

    assert(curve);
    assert(iFuzzyComp(beg_param, end_param) == -1);

    CubitVector vert_coord(pt[0], pt[1], pt[2]), closest_coord;
    curve->get_curve_geom()->closest_coord_on_curve(vert_coord, closest_coord,
						    beg_param, end_param);

    close[0] = closest_coord.x();
    close[1] = closest_coord.y();
    close[2] = closest_coord.z();
  }
  else {
    // Work from the discrete edge directly.
    // Find the normal to the edge.
    double norm[3];
    calcUnitNormal(norm);

    // Find projection distance to the edge.
    double startVec[] =
      { getVert(0)->x() - pt[0], getVert(0)->y() - pt[1], getVert(0)->z()
	  - pt[2] };
    double projDist = -dDOT3D(norm, startVec);
    double projPt[] =
      { pt[0] - projDist * norm[0], pt[1] - projDist * norm[1], pt[2]
	  - projDist * norm[2] };
    int orientStart = checkOrient2D(pt, projPt, getVert(0)->getCoords());
    int orientEnd = checkOrient2D(pt, projPt, getVert(1)->getCoords());
    if (orientStart * orientEnd == -1) {
      // pt projects into the interior of the edge
      close[0] = projPt[0];
      close[1] = projPt[1];
      close[2] = projPt[2];
      return;
    }

    // Compute d^2 for each endpoint.
    double startDist2 = dMAG3D_SQ(startVec);
    double endDist2 = dDIST_SQ_3D(getVert(1)->getCoords(), pt);

    // Take the smaller
    if (startDist2 < endDist2) {
      close[0] = getVert(0)->x();
      close[1] = getVert(0)->y();
      close[2] = getVert(0)->z();
      return;
    }
    else {
      close[0] = getVert(1)->x();
      close[1] = getVert(1)->y();
      close[2] = getVert(1)->z();
      return;
    }
  }
}

// CubitVector BdryEdgeBase::
// get_offset_loc(const CubitVector&) const {

//   assert(0);

//   CubitVector v;

//   switch(pVVert(0)->iSpaceDimen()) {

//   case 2: {
//     Vert dummy_vert;
//     double dummy_coord[] = { trigger.x(), trigger.y() };
//     dummy_vert.vSetCoords(2, dummy_coord );

//     int left = iOrient2D(pVVert(0), pVVert(1), &dummy_vert);
//     if(left == 1)
//       v = get_offset_loc(true);
//     else if(left == -1)
//       v = get_offset_loc(false);
//     else
//       v = get_split_loc();
//     break;
//   }

//   case 3: { //not implemented yet
//     vFatalError("This function is not implemented in 3D",
// 		"BdryEdge::get_offset_loc(const double*)");
//     v.set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
//     break;
//   }
//   default: {
//     assert(0);
//     v.set(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
//     break;
//   }

//   }

//   return v;

// }

// CubitVector BdryEdgeBase::
// get_offset_loc(const bool) const {

//   assert(0);

//   CubitVector offset_loc;

//   switch(pVVert(0)->iSpaceDimen()) {

//   case 2: {
//     GRCurve* c = dynamic_cast<GRCurve*>(curve);

//     int reg_left, reg_right;
//     if(sense == CUBIT_FORWARD) {
//       reg_left  = c->region_left();
//       reg_right = c->region_right();
//     }
//     else {
//       assert(sense == CUBIT_REVERSED);
//       reg_right = c->region_left();
//       reg_left  = c->region_right();
//     }

//     if(left && reg_left == 0)
//       vFatalError("There is no region on the left side of this curve",
//  		  "BdryEdge::get_offset_loc(const bool)");
//     if(!left && reg_right == 0)
//       vFatalError("There is no region on the right side of this curve",
//  		  "BdryEdge::get_offset_loc(const bool)");

//     if(left) assert(bdry_cond(false) != 0);
//     else     assert(bdry_cond(true)  != 0);

//     const Vert *v1 = pVVert(0), *v2 = pVVert(1);
//     double mid[]   = { 0.5 * (v1->dX() + v2->dX()), 0.5 * (v1->dY() + v2->dY()) };
//     double norm[]  = { v1->dY() - v2->dY(), v2->dX() - v1->dX() };
//     double mag     = dMAG2D(norm);
//     double rad     = mag / sqrt(3.);
//     norm[0] /= mag; norm[1] /= mag;

//     if(left) offset_loc.set( mid[0] + 0.50010001 * rad * norm[0],
// 			     mid[1] + 0.50010001 * rad * norm[1],
// 			     0. );
//     else     offset_loc.set( mid[0] - 0.50010001 * rad * norm[0],
// 			     mid[1] - 0.50010001 * rad * norm[1],
// 			     0 );
//     break;
//   }

//   return offset_loc;

// }
