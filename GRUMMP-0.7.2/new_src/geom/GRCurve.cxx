#include "GR_config.h"
#include "CastTo.hpp"
#include "CoEdgeSM.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"
#include "TopologyBridge.hpp"
#include "GR_GRLoop.h"
#include "GR_GRCoEdge.h"
#include "GR_GRCurve.h"
#include "GR_GRPoint.h"
#include "GR_misc.h"

GRCurve::GRCurve() :
    start_point(NULL), end_point(NULL), coedge_list(), curve_geom(NULL), attrib_set(), sense(
	CUBIT_UNKNOWN), periodic(false)
{
  coedge_list.clean_out();
}

GRCurve::GRCurve(GRCurveGeom* cg, TBPoint* start, TBPoint* end) :
    start_point(start), end_point(end), coedge_list(), curve_geom(cg), attrib_set(), sense(
	CUBIT_FORWARD), periodic(start == end)
{
  coedge_list.clean_out();
}

GRCurve::~GRCurve()
{
  remove_all_simple_attribute_virt();
  coedge_list.clean_out();

  delete_geometry();
  curve_geom = NULL;

  start_point = NULL;
  end_point = NULL;
}

//double GRCurve::
//get_arc_length(const CubitVector& p, const int which_end) {
//
//	assert(qValid());
//	assert(curve_geom->coord_on_curve(p));
//	double param = curve_geom->param_at_coord(p);
//	if(which_end) {
//		double max_param = curve_geom->max_param();
//		return curve_geom->arc_length(param, max_param);
//	}
//	else {
//		double min_param = curve_geom->min_param();
//		return curve_geom->arc_length(min_param, param);
//	}
//
//	assert(0);
//	return LARGE_DBL;
//
//}

CubitStatus
GRCurve::position_from_fraction(const double fraction, CubitVector& position)
{
  assert(qValid());
  assert(iFuzzyComp(fraction, 0.) >= 0 && iFuzzyComp(fraction, 1.) <= 0);

  double minp = curve_geom->min_param();
  double maxp = curve_geom->max_param();

  curve_geom->coord_at_param(minp + fraction * (maxp - minp), position);
  return CUBIT_SUCCESS;
}

CubitBoolean
GRCurve::is_periodic(double& period)
{
  assert(qValid());
  if (periodic) {
    period = 1.;
    return CUBIT_TRUE;
  }
  else {
    return CUBIT_FALSE;
  }
}

CubitStatus
GRCurve::closest_point(CubitVector const& location,
		       CubitVector& closest_location, CubitVector* tangent,
		       CubitVector* curvature, double* param)
{
  assert(qValid());
  curve_geom->closest_coord_on_curve(location, closest_location);
  if (tangent != NULL) {
    curve_geom->unit_tangent(curve_geom->param_at_coord(closest_location),
			     *tangent);
  }
  if (curvature != NULL) {
    curve_geom->curvature(curve_geom->param_at_coord(closest_location),
			  *curvature);
  }
  if (param != NULL) {
    *param = curve_geom->param_at_coord(closest_location);
  }
  return CUBIT_SUCCESS;
}

CubitStatus
GRCurve::get_center_radius(CubitVector& center, double& radius)
{
  assert(qValid());
  if (!dynamic_cast<GRArc*>(curve_geom)) {
    assert(0);
    return CUBIT_FAILURE;
  }
  else {
    dynamic_cast<GRArc*>(curve_geom)->get_center(center);
    radius = dynamic_cast<GRArc*>(curve_geom)->get_radius();
    return CUBIT_SUCCESS;
  }
}

void
GRCurve::get_parents_virt(DLIList<TopologyBridge*>& parents)
{
  assert(qValid());
  CAST_LIST_TO_PARENT(coedge_list, parents);
}

void
GRCurve::get_children_virt(DLIList<TopologyBridge*>& children)
{
  assert(qValid());
  children.append(start_point);
  if (start_point != end_point) {
    children.append(end_point);
  }
}

void
GRCurve::get_loops(DLIList<GRLoop*>& loop_list)
{
  assert(qValid());
  DLIList<GRCoEdge*> ce_list;
  get_coedges(ce_list);
  ce_list.reset();

  for (int i = ce_list.size(); i--;) {
    GRCoEdge* ce = ce_list.get_and_step();
    GRLoop* loop = dynamic_cast<GRLoop*>(ce->get_loop());
    if (loop)
      loop_list.append_unique(loop);
  }
}

void
GRCurve::get_coedges(DLIList<GRCoEdge*>& ce_list)
{
  assert(qValid());
  int i;
  coedge_list.reset();
  for (i = 0; i < coedge_list.size(); i++)
    if (GRCoEdge* ce = dynamic_cast<GRCoEdge*>(coedge_list.next(i)))
      ce_list.append(ce);
}

void
GRCurve::get_curves(DLIList<GRCurve*>& curve_list)
{
  assert(qValid());
  curve_list.append(this);
}

void
GRCurve::get_points(DLIList<GRPoint*>& point_list)
{
  assert(qValid());
  GRPoint* p;
  if ((p = dynamic_cast<GRPoint*>(start_point)))
    point_list.append(p);
  if ((start_point != end_point) && (p = dynamic_cast<GRPoint*>(end_point)))
    point_list.append(p);
}

int
GRCurve::region_left()
{

  assert(qValid());
  int i, return_value;
  DLIList<CubitSimpleAttrib*> csa_list;
  CubitString cond_string("CURVE_CONDITION");

  get_simple_attribute(cond_string, csa_list);
  csa_list.reset();
  assert(csa_list.size() == 1);
  csa_list.next(0)->string_data_list()->reset();
  csa_list.next(0)->int_data_list()->reset();

  CubitString s("REG_LEFT");
  for (i = 0; i < csa_list.next(0)->string_data_list()->size(); i++) {
    if (*(csa_list.next(0)->string_data_list()->next(i)) == s) {
      return_value = *(csa_list.next(0)->int_data_list()->next(i));
      break;
    }
    if (i + 1 == csa_list.next(0)->string_data_list()->size()) {
      logMessage(
	  0, "ERROR: The curve's left region is not initialized properly.\n");
      abort();
    }
  }

  delete csa_list.next(0);

  return return_value;

}

int
GRCurve::region_right()
{

  assert(qValid());
  int i, return_value;
  DLIList<CubitSimpleAttrib*> csa_list;
  CubitString cond_string("CURVE_CONDITION");

  get_simple_attribute(cond_string, csa_list);
  csa_list.reset();
  assert(csa_list.size() == 1);
  csa_list.next(0)->string_data_list()->reset();
  csa_list.next(0)->int_data_list()->reset();

  CubitString s("REG_RIGHT");
  for (i = 0; i < csa_list.next(0)->string_data_list()->size(); i++) {
    if (*(csa_list.next(0)->string_data_list()->next(i)) == s) {
      return_value = *(csa_list.next(0)->int_data_list()->next(i));
      break;
    }
    if (i + 1 == csa_list.next(0)->string_data_list()->size()) {
      logMessage(
	  0, "ERROR: The curve's right region is not initialized properly.\n");
      abort();
    }
  }

  delete csa_list.next(0);

  return return_value;

}

int
GRCurve::boundary_left()
{

  assert(qValid());
  int i, return_value;
  DLIList<CubitSimpleAttrib*> csa_list;
  CubitString cond_string("CURVE_CONDITION");

  get_simple_attribute(cond_string, csa_list);
  csa_list.reset();
  assert(csa_list.size() == 1);
  csa_list.next(0)->string_data_list()->reset();
  csa_list.next(0)->int_data_list()->reset();

  CubitString s("BDRY_LEFT");
  for (i = 0; i < csa_list.next(0)->string_data_list()->size(); i++) {
    if (*(csa_list.next(0)->string_data_list()->next(i)) == s) {
      return_value = *(csa_list.next(0)->int_data_list()->next(i));
      break;
    }
    if (i + 1 == csa_list.next(0)->string_data_list()->size()) {
      logMessage(
	  0, "ERROR: The curve's left boundary is not initialized properly.\n");
      abort();
    }
  }

  delete csa_list.next(0);

  return return_value;

}

int
GRCurve::boundary_right()
{

  assert(qValid());
  int i, return_value;
  DLIList<CubitSimpleAttrib*> csa_list;
  CubitString cond_string("CURVE_CONDITION");

  get_simple_attribute(cond_string, csa_list);
  csa_list.reset();
  assert(csa_list.size() == 1);
  csa_list.next(0)->string_data_list()->reset();
  csa_list.next(0)->int_data_list()->reset();

  CubitString s("BDRY_RIGHT");
  for (i = 0; i < csa_list.next(0)->string_data_list()->size(); i++) {
    if (*(csa_list.next(0)->string_data_list()->next(i)) == s) {
      return_value = *(csa_list.next(0)->int_data_list()->next(i));
      break;
    }
    if (i + 1 == csa_list.next(0)->string_data_list()->size()) {
      logMessage(
	  0,
	  "ERROR: The curve's right boundary is not initialized properly.\n");
      abort();
    }
  }

  delete csa_list.next(0);

  return return_value;

}
void
GRCurve::set_left_region(int iRegion)
{
  assert(qValid());

  DLIList<CubitSimpleAttrib*> csa_list;
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);
  string_list.append(&reg_left);
  double_list.append(0.);
  int_list.append(iRegion);
  string_list.append(&bdry_left);
  double_list.append(0.);
  int_list.append(boundary_left());
  string_list.append(&reg_right);
  double_list.append(0.);
  int_list.append(region_right());
  string_list.append(&bdry_right);
  double_list.append(0.);
  int_list.append(boundary_right());

  get_simple_attribute(cond_string, csa_list);

  csa_list.reset();
  assert(csa_list.size() == 1);
  CubitSimpleAttrib * csa = csa_list.next(0);
  remove_simple_attribute_virt(csa);
  delete csa;
  CubitSimpleAttrib csaNew(&string_list, &double_list, &int_list);
  append_simple_attribute_virt(&csaNew);
}

void
GRCurve::set_right_region(int iRegion)
{
  assert(qValid());

  DLIList<CubitSimpleAttrib*> csa_list;
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);
  string_list.append(&reg_left);
  double_list.append(0.);
  int_list.append(region_left());
  string_list.append(&bdry_left);
  double_list.append(0.);
  int_list.append(boundary_left());
  string_list.append(&reg_right);
  double_list.append(0.);
  int_list.append(iRegion);
  string_list.append(&bdry_right);
  double_list.append(0.);
  int_list.append(boundary_right());

  get_simple_attribute(cond_string, csa_list);

  csa_list.reset();
  assert(csa_list.size() == 1);
  CubitSimpleAttrib * csa = csa_list.next(0);
  remove_simple_attribute_virt(csa);

  CubitSimpleAttrib csaNew(&string_list, &double_list, &int_list);
  append_simple_attribute_virt(&csaNew);

}

bool
GRCurve::qValid() const
{
  assert(this);
  if (curve_geom == NULL) {
    logMessage(
	0, "ERROR: The curve geometry object of this curve is unitialized\n");
    return false;
  }
  if (start_point == NULL) {
    logMessage(0, "ERROR: The start point of this curve is unitialized\n");
    return false;
  }
  if (end_point == NULL) {
    logMessage(0, "ERROR: The end point of this curve is unitialized\n");
    return false;
  }

  return true;

}
