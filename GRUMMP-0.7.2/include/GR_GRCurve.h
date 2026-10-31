#ifndef GR_GRCURVE_H
#define GR_GRCURVE_H 1

#include "GR_config.h"
#include "GR_misc.h"
#include "CubitDefines.h"   //Cubit definition
#include "Curve.hpp"        //CGM pure virtual curve class
#include "DLIList.hpp"      //CGM Doubly-linked list class
#include "GR_GRAttrib.h"    //GRUMMP attributes
#include "GR_GRAttribSet.h" //GRUMMP attribute set
#include "GR_GRCurveGeom.h" //GRUMMP curve geometry interface

//Forward declarations
class CoEdgeSM;
class CubitVector;
class TBPoint;
class GRPoint;
class GRCoEdge;
class GRLoop;
class TopologyBridge;

class GRCurve : public Curve {

private:

  TBPoint* start_point;               //Pointer to the start point
  TBPoint* end_point;                 //Pointer to the end point
  DLIList<CoEdgeSM*> coedge_list;   //Pointer to list of coedges of this curve

  GRCurveGeom* curve_geom;          //Curve geometry
  GRAttribSet attrib_set;           //Curve attributes
  CubitSense sense;                 //Sense w/r geometry
  bool periodic;                   //True if start and end points are the same

public:

  //Constructors
  GRCurve();
  GRCurve(GRCurveGeom* cg, TBPoint* start, TBPoint* end);
  //GRCurve(GRCurveGeom* cg, TBPoint* start, TBPoint* end, DLIList<CoEdgeSM*>& coedges);
  //GRCurve(GRCurveGeom* cg, TBPoint* start, TBPoint* end, CubitSense curve_sense);

private:

  //Copy constructor and operator= are private for now.
  // GRCurve(const GRCurve&); // Never defined
  // GRCurve& operator=(const GRCurve&); // Never defined

public:

  //Destructor
  virtual
  ~GRCurve();

  //Deletes the geometry object pointed by curve_geom
  void
  delete_geometry();
  //Returns the curve evaluation tool.
  GRCurveGeom*
  get_curve_geom();
  //Sets the geometry for this curve.
  //void set_curve_geom(GRCurveGeom* cg);
//  //Sets start point
//  void set_start_point(TBPoint* start);
//  //Sets end point
//  void set_end_point(TBPoint* end);
//  //Sets both end points.
//  void set_points(TBPoint* start, TBPoint* end);
  //Add one co-edge to the co-edge list.
  void
  add_coedge(CoEdgeSM* ce);
  //Add a list of co-edges to the co-edge list.
  void
  add_coedges(DLIList<CoEdgeSM*>& ce_list);

  //Pure virtual in base class. Attribute management.
  virtual void
  append_simple_attribute_virt(CubitSimpleAttrib* csa);

  virtual void
  remove_simple_attribute_virt(CubitSimpleAttrib* csa);

  virtual void
  remove_all_simple_attribute_virt();

  virtual CubitStatus
  get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list);

  virtual CubitStatus
  get_simple_attribute(const CubitString& name,
		       DLIList<CubitSimpleAttrib*>& csa_list);

  //Write GRAttribs out to file
  CubitStatus
  save_attribs(FILE* file);
  //Read GRAttribs from file
  CubitStatus
  restore_attribs(FILE* file, unsigned int endian);

  //Returns a box bounding the curve (useful for self-intersection check)
  virtual CubitBox
  bounding_box() const;
  //Returns the GRUMMP geometry engine (which does not exist).
  virtual GeometryQueryEngine*
  get_geometry_query_engine() const
  {
    vFatalError("This function is not supported in GRUMMP",
		"GRCurve::get_geometry_query_engine()");
    return NULL;
  }

  //Returns the length of the curve (pure virtual)
  virtual double
  measure();

  //Returns the arc length from param1 to param2
  virtual double
  length_from_u(double param1, double param2);

//  //Return the arc length along the entire curve (somewhat redundant with measure()).
  virtual double
  get_arc_length();

  // these exist to avoid warnings. Don't call them. Less to maintain.
  virtual double
  get_arc_length(const CubitVector&, const CubitVector&)
  {
    assert(0);
    return 0;
  }
  virtual double
  get_arc_length(const CubitVector&, int)
  {
    assert(0);
    return 0;
  }
  //Returns the center point of the curve based on arc length.
  virtual CubitVector
  center_point();

  //Calculates the position at a fraction of the entire curve in parameter space.
  virtual CubitStatus
  position_from_fraction(const double fraction, CubitVector& position);

  //Returns CUBIT_TRUE if the parameterization is periodic.
  //THE PERIOD RETURNED ALWAYS EQUALS 1 (not correctly implemented yet)
  virtual CubitBoolean
  is_periodic(double& period);

  //Returns parameter range of the curve. Takes orientation into account
  //so lower_bound can be greater than upper_bound.
  virtual CubitBoolean
  get_param_range(double& lower_bound, double& upper_bound);

  virtual CubitStatus
  get_ellipse_params(CubitVector &/*center*/, CubitVector &/*normal*/,
		     CubitVector &/*major_axis*/,
		     double &/*radius_ratio*/) const
  {
    assert(0);
    return CUBIT_FAILURE;
  }

  virtual CubitStatus
  get_spline_params(bool &/*rational*/, int &/*degree*/,
		    DLIList<CubitVector> &/*ctrl_pts*/,
		    DLIList<double> &/*ctrl_pt_wts*/,
		    DLIList<double> &/*knots*/) const
  {
    assert(0);
    return CUBIT_FAILURE;
  }
  //Calculates the closest location on a curve to a given location
  //Can optionally return the tangent, curvature and parameter at this point.
  // - The tangent is always along the positive direction of the curve.
  virtual CubitStatus
  closest_point(CubitVector const& location, CubitVector& closest_location,
		CubitVector* tangent = NULL, CubitVector* curvature = NULL,
		double* param = NULL);

  //Gets the location at a given parameter on the curve.
  virtual CubitStatus
  position_from_u(double param, CubitVector& position);

  //Returns the parameter at a given location.
  double
  u_from_position(const CubitVector& position);

  //Returns the parameter of the point arc_length away from start point
  //located at root_param.
  virtual double
  u_from_arc_length(double root_param, double arc_length);

  //Verifies that a point is lying directly on a curve
  virtual CubitBoolean
  is_position_on(const CubitVector& position);

  //Finds the underlying arc's center point and radius
  //ONLY VALID FOR ARCS!
  virtual CubitStatus
  get_center_radius(CubitVector& center, double& radius);

  //Returns CUBIT_PNT_OUTSIDE, CUBIT_PNT_INSIDE, CUBIT_PNT_BOUNDARY or CUBIT_PNT_UNKNOWN
  //depending if the point is inside or outside the boundary. u and v are parametric
  //coords on a trimmed surface. NOT YET IMPLEMENTED!
  virtual CubitPointContainment
  point_containment(const CubitVector&)
  {
    assert(0);
    return CUBIT_PNT_UNKNOWN;
  }
  virtual CubitPointContainment
  point_containment(double, double)
  {
    assert(0);
    return CUBIT_PNT_UNKNOWN;
  }
  virtual CubitPointContainment
  point_containment(CubitVector&, double, double)
  {
    assert(0);
    return CUBIT_PNT_UNKNOWN;
  }

  virtual CubitStatus
  get_interior_extrema(DLIList<CubitVector*>& /* interior_points */,
		       CubitSense& /* return_sense */)
  {
    return CUBIT_FAILURE;
  }
  //Returns start and end param for the curve.
  virtual double
  start_param();
  virtual double
  end_param();

  //Adjacency queries.
  /*   void get_bodies  ( DLIList<GRBody   *>& bodies   ); */
  /*   void get_lumps   ( DLIList<GRLump   *>& lumps    ); */
  /*   void get_shells  ( DLIList<GRShell  *>& shells   ); */
  /*   void get_surfaces( DLIList<GRSurface*>& surfaces ); */
  void
  get_loops(DLIList<GRLoop*>& loops);
  void
  get_coedges(DLIList<GRCoEdge*>& coedges);
  void
  get_curves(DLIList<GRCurve*>& curve_list);
  void
  get_points(DLIList<GRPoint*>& point_list);

  //Returns a list of parent and child topological entities.
  void
  get_parents_virt(DLIList<TopologyBridge*>& parents);
  void
  get_children_virt(DLIList<TopologyBridge*>& children);

  //Gets pointer to start and end points.
  TBPoint*
  get_start_point();
  TBPoint*
  get_end_point();

  //Gets the sense of the curve.
  CubitSense
  get_sense();

  //Boundary conditions and regions
  int
  region_left();
  int
  region_right();
  int
  boundary_left();
  int
  boundary_right();
  // lets you reset regions well after creation
  void
  set_left_region(int iRegion);
  void
  set_right_region(int iRegion);
  //Checks if this curve has co-edges.
  bool
  has_parent_coedge();

  bool
  qValid() const;

};

//Inline function definitions

inline void
GRCurve::delete_geometry()
{
  if (curve_geom)
    delete curve_geom;
  curve_geom = NULL;
}

inline void
GRCurve::add_coedge(CoEdgeSM* ce)
{
  coedge_list.append_unique(ce);
}
inline void
GRCurve::add_coedges(DLIList<CoEdgeSM*>& ce_list)
{
  assert(qValid());
  coedge_list += ce_list;
}

inline GRCurveGeom*
GRCurve::get_curve_geom()
{
  return curve_geom;
}

inline TBPoint*
GRCurve::get_start_point()
{
  return start_point;
}
inline TBPoint*
GRCurve::get_end_point()
{
  return end_point;
}

inline CubitSense
GRCurve::get_sense()
{
  return sense;
}

inline bool
GRCurve::has_parent_coedge()
{
  return coedge_list.size() > 0;
}

inline void
GRCurve::append_simple_attribute_virt(CubitSimpleAttrib* csa)
{
  assert(qValid());
  attrib_set.append_attribute(csa);
}

inline void
GRCurve::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
{
  assert(qValid());
  attrib_set.remove_attribute(csa);
}

inline void
GRCurve::remove_all_simple_attribute_virt()
{
  assert(qValid());
  attrib_set.remove_all_attributes();
}

inline CubitStatus
GRCurve::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
{
  assert(qValid());
  return attrib_set.get_attributes(csa_list);
}

inline CubitStatus
GRCurve::get_simple_attribute(const CubitString& name,
			      DLIList<CubitSimpleAttrib*>& csa_list)
{
  assert(qValid());
  return attrib_set.get_attributes(name, csa_list);
}

inline CubitStatus
GRCurve::save_attribs(FILE* file)
{
  assert(qValid());
  return attrib_set.save_attributes(file);
}

inline CubitStatus
GRCurve::restore_attribs(FILE* file, unsigned int endian)
{
  assert(qValid());
  return attrib_set.restore_attributes(file, endian);
}

inline CubitBox
GRCurve::bounding_box() const
{
  assert(qValid());
  return curve_geom->bounding_box();
}

inline double
GRCurve::measure()
{
  assert(qValid());
  return curve_geom->arc_length();
}

inline double
GRCurve::length_from_u(double param1, double param2)
{
  assert(qValid());
  return curve_geom->arc_length(param1, param2);
}

inline double
GRCurve::get_arc_length()
{
  return curve_geom->arc_length();
}

inline CubitVector
GRCurve::center_point()
{
  assert(qValid());
  CubitVector CP;
  curve_geom->center_point(curve_geom->min_param(), curve_geom->max_param(),
			   CP);
  return CP;
}

inline CubitBoolean
GRCurve::get_param_range(double& lower_bound, double& upper_bound)
{
  assert(qValid());
  lower_bound = curve_geom->min_param();
  upper_bound = curve_geom->max_param();
  return CUBIT_TRUE;
}

inline CubitStatus
GRCurve::position_from_u(double param, CubitVector& position)
{
  assert(qValid());
  curve_geom->coord_at_param(param, position);
  return CUBIT_SUCCESS;
}

inline double
GRCurve::u_from_position(const CubitVector& position)
{
  assert(qValid());
  return curve_geom->param_at_coord(position);
}

inline double
GRCurve::u_from_arc_length(double root_param, double arc_length)
{
  assert(qValid());
  return curve_geom->param_at_arc_length(root_param, arc_length);
}

inline CubitBoolean
GRCurve::is_position_on(const CubitVector& position)
{
  assert(qValid());
  if (curve_geom->coord_on_curve(position))
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

inline double
GRCurve::start_param()
{
  assert(qValid());
  return curve_geom->min_param();
}

inline double
GRCurve::end_param()
{
  assert(qValid());
  return curve_geom->max_param();
}

#endif
