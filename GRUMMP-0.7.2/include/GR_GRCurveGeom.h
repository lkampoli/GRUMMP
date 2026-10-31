#ifndef GR_GRCURVEGEOM_H
#define GR_GRCURVEGEOM_H 1

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Face.h"
#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include <set>
#include <algorithm>
#include <iterator>
#include <functional>

class CubitBox;

class GRCurveGeom {

public:

  enum curve_type {
    UNKNOWN, LINE, ARC, CUBIC, BEZIER, INTERP
  };

protected:

  int n_pts;           //Number of points defining the curve
  CubitVector* points; //Points locations
  int n_dim;           //Curve geometric dimension
  double cg_tolerance; //Used for anything iterative
  mutable double m_minRadiusOfCurvature;
  //Functor used to compare parameters on a curve. Two params
  //that are almost equal (to a small tolerance) are considered equal.
  struct param_compare {
    bool
    operator()(const double& d1, const double& d2) const
    {
      if (iFuzzyComp(d1, d2) == -1) {
	return true;
      }
      else {
	return false;
      }
    }
  };

private:

  //Disabling copy construction and operator=.
  GRCurveGeom(const GRCurveGeom&); // Never defined
  GRCurveGeom&
  operator=(const GRCurveGeom&); // Never defined

public:

  GRCurveGeom() :
      n_pts(0), points(NULL), n_dim(0), cg_tolerance(1e-12), m_minRadiusOfCurvature(0)
  {
  }
  virtual
  ~GRCurveGeom()
  {
    if (points)
      delete[] points;
  }

  //Returns the curve number of topological dimensions.
  int
  topo_dim() const
  {
    return 1;
  }
  void
  set_tolerance(double tolerance)
  {
    cg_tolerance = tolerance;
  }
  //Returns the curve number of geometric dimensions.
  virtual int
  geom_dim() const
  {
    return n_dim;
  }
  int
  get_num_points() const
  {
    return n_pts;
  }
  //Returns the type of geometric entity.
  virtual curve_type
  get_curve_type() const = 0;

  //Returns the min and max param for this curve;
  virtual double
  min_param() const = 0;
  virtual double
  max_param() const = 0;

  //Calculates the bounding box for this curve.
  virtual CubitBox
  bounding_box() const = 0;

  //Returns parameter at a given coordinate.
  virtual double
  param_at_coord(const CubitVector& coord) const = 0;

  //Gets the parameter located arc length "length" from "param_start"
  virtual double
  param_at_arc_length(const double param_start, const double length) const = 0;

  //Returns coordinate at a given parameter.
  virtual void
  coord_at_param(const double param, CubitVector& coord) const = 0;

  //Computes coordinate at a given Euclidean distance from start or end
  //Returns corresponding parameter.
  virtual double
  coord_at_dist(const double dist, const bool start,
		CubitVector& coord) const = 0;

  //Computes the coordinate at middle Euclidean distance from two coords or params on the curve.
  //Returns corresponding parameter.
  virtual double
  coord_at_mid_dist(const double param1, const double param2,
		    CubitVector& middle) const = 0;
  //  virtual double coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
  //				   CubitVector& middle) const = 0;

  //Given a coordinate near the curve, computes neareast coord on the curve.
  //Optionally, can specify a param range between which to restrict the search.
  virtual void
  closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			 const double param1 = -1.,
			 const double param2 = -1.) const = 0;

  //Given a coordinate near the curve, computes neareast parameter on the curve.
  //Optionally, can specify a param range between which to restrict the search.
  virtual double
  closest_param_on_curve(const CubitVector& coord, const double param1 = -1.,
			 const double param2 = -1.) const = 0;

  //Returns true if coord is directly on the curve.
  virtual bool
  coord_on_curve(const CubitVector& coord, double* const param_X = NULL,
		 double* const param_Y = NULL,
		 double* const param_Z = NULL) const = 0;

  //Returns the first and second derivative with respect to the parameter
  virtual void
  first_deriv(const double param, CubitVector& FD) const = 0;
  //virtual void first_deriv(const CubitVector& coord, CubitVector& FD) const = 0;

  virtual void
  second_deriv(const double param, CubitVector& SD) const = 0;
  // virtual void second_deriv(const CubitVector& coord, CubitVector& SD) const = 0;

  //Curvature at param and at coord (located on the curve).
  virtual double
  curvature(const double param) const = 0;
  //  virtual double curvature(const CubitVector& coord) const = 0;

  //Curvature vector
  virtual void
  curvature(const double param, CubitVector& curvature) const = 0;
  //  virtual void curvature(const CubitVector& coord, CubitVector& curvature) const = 0;

  double
  getMinRadiusOfCurvature() const {
	  if (m_minRadiusOfCurvature == 0) {
		  m_minRadiusOfCurvature = LARGE_DBL;
		  for (double para = min_param(); para < max_param();
				  para += (max_param() - min_param())
				/ (get_num_points()*20)) {
			  m_minRadiusOfCurvature = std::min(m_minRadiusOfCurvature,
					  1. / curvature(para));
		  }
	  }
	  return m_minRadiusOfCurvature;
  }

  //Unit normal points to the left side of the curve.
  //ONLY WORKS IN 2D!!!
  virtual void
  unit_normal(const double param, CubitVector& normal) const = 0;
  //  virtual void unit_normal(const CubitVector& coord, CubitVector& normal) const = 0;

  //Tangent vector has the same "sense" as the curve.
  virtual void
  unit_tangent(const double param, CubitVector& tangent) const = 0;
  //  virtual void unit_tangent(const CubitVector& coord, CubitVector& tangent) const = 0;

  //Total length, length between dParam1 and dParam2, between adCoord1 and adCoord2.
  virtual double
  arc_length() const = 0;
  virtual double
  arc_length(const double param1, const double param2) const = 0;
  //  virtual double arc_length(const CubitVector& coord1, const CubitVector& coord2) const = 0;

  //Rerturns the center point based on arc length. Either based on total
  //arc length or between params or points.
  //  virtual void center_point(CubitVector& mid_point) const = 0;
  virtual void
  center_point(const double param1, const double param2,
	       CubitVector& mid_point) const = 0;
  //  virtual void center_point(const CubitVector& coord1, const CubitVector& coord2,
  //			    CubitVector& mid_point) const = 0;

  //Rerturns the point located in the middle of the parameterization range.
  //  void mid_point(CubitVector& mp) const;
  //  void mid_point(const double param1, const double param2,
  //		 CubitVector& mp) const;
  //  void mid_point(const CubitVector& coord1, const CubitVector& coord2,
  //		 CubitVector& mp) const;

  //Total variation of the tangent angle
  //  virtual double TVT() const = 0;
  virtual double
  TVT(const double param1, const double param2) const = 0;
  //  virtual double TVT(const CubitVector& coord1, const CubitVector& coord2) const = 0;

  //Finds the param at which the TVT is half of what it is between param1 and param2.
  virtual double
  mid_TVT(const double param1, const double param2) const = 0;
  //  virtual double mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const = 0;

  //Returns the local interior extrema for the curve.
  virtual void
  interior_extrema(DLIList<CubitVector*>& point_list) const = 0;

  virtual bool
  line_intersect(const double coords1[3], const double coords2[3],
		 const double begParam, const double endParam,
		 std::vector<double> & params) const = 0;
  //Gets points[pt] in some form or another.
  virtual double
  x(const int pt) const = 0;
  virtual double
  y(const int pt) const = 0;
  virtual double
  z(const int pt) const = 0;
  virtual void
  get_coord(const int pt, CubitVector& coord) const = 0;
  virtual const CubitVector&
  get_coord(const int pt) const = 0;

  //Sets Points[iPt] to new location
  virtual void
  set_x(const int pt, const double x) = 0;
  virtual void
  set_y(const int pt, const double y) = 0;
  virtual void
  set_z(const int pt, const double z) = 0;
  virtual void
  set_coord(const int pt, const CubitVector& new_coord) = 0;

  virtual bool
  isValid() const = 0;
  virtual bool
  closed_param()
  {
    return false;
  }

};

class GRLine : public GRCurveGeom {

public:

  GRLine();
  GRLine(const CubitVector& coord_start, const CubitVector& coord_end);
  GRLine(const GRLine& L);
  GRLine&
  operator=(const GRLine& L);
  virtual
  ~GRLine();

  GRCurveGeom::curve_type
  get_curve_type() const
  {
    return GRCurveGeom::LINE;
  }

  double
  min_param() const
  {
    return 0.;
  }
  double
  max_param() const
  {
    return 1.;
  }

  CubitBox
  bounding_box() const;

  double
  param_at_coord(const CubitVector& coord) const;

  double
  param_at_arc_length(const double param_start, const double length) const;

  void
  coord_at_param(const double param, CubitVector& coord) const;

  double
  coord_at_dist(const double dist, const bool start, CubitVector& coord) const;

  double
  coord_at_mid_dist(const double param1, const double param2,
		    CubitVector& middle) const;
  //  double coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
  //			   CubitVector& middle) const;

  void
  closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			 const double param1 = -1.,
			 const double param2 = -1.) const;

  double
  closest_param_on_curve(const CubitVector& coord, const double param1 = -1.,
			 const double param2 = -1.) const;

  bool
  coord_on_curve(const CubitVector& coord, double* const param_X = NULL,
		 double* const param_Y = NULL,
		 double* const param_Z = NULL) const;

  void
  first_deriv(const double param, CubitVector& FD) const;
  //void first_deriv(const CubitVector& coord, CubitVector& FD) const;

  void
  second_deriv(const double param, CubitVector& SD) const;
  //  void second_deriv(const CubitVector& coord, CubitVector& SD) const;

  double
  curvature(const double param) const;
  //  double curvature(const CubitVector& coord) const;
  void
  curvature(const double param, CubitVector& curvature) const;
  //  void curvature(const CubitVector& coord, CubitVector& curvature) const;

  void
  unit_normal(const double param, CubitVector& normal) const;
  //  void unit_normal(const CubitVector& coord, CubitVector& normal) const;

  void
  unit_tangent(const double param, CubitVector& tangent) const;
  //  void unit_tangent(const CubitVector& coord, CubitVector& tangent) const;

  double
  arc_length() const;
  double
  arc_length(const double param1, const double param2) const;
  //  double arc_length(const CubitVector& coord1, const CubitVector& coord2) const;

  //void center_point(CubitVector& mid_point) const;
  void
  center_point(const double param1, const double param2,
	       CubitVector& mid_point) const;
  //  void center_point(const CubitVector& coord1, const CubitVector& coord2,
  //		    CubitVector& mid_point) const;

  //  double TVT() const;
  double
  TVT(const double param1, const double param2) const;
  //  double TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  double
  mid_TVT(const double param1, const double param2) const;
  //  double mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  void
  interior_extrema(DLIList<CubitVector*>& point_list) const;

  bool
  line_intersect(const double coords1[3], const double coords2[3],
		 const double begParam, const double endParam,
		 std::vector<double> & params) const;

  double
  x(const int pt) const;
  double
  y(const int pt) const;
  double
  z(const int pt) const;
  void
  get_coord(const int pt, CubitVector& coord) const;
  const CubitVector&
  get_coord(const int pt) const;

  //Sets Points[iPt] to new location
  void
  set_x(const int pt, const double x);
  void
  set_y(const int pt, const double y);
  void
  set_z(const int pt, const double z);
  void
  set_coord(const int pt, const CubitVector& new_coord);

  bool
  isValid() const;

};

class GRArc : public GRCurveGeom {

  double radius;
  double angle_beg;
  double angle_span;

public:

  GRArc();
  GRArc(const CubitVector& center, const CubitVector& beg_pt,
	const CubitVector& end_pt);
  GRArc(const CubitVector& beg_pt, const CubitVector& end_pt, const double rad,
	const bool short_arc = true);
  GRArc(const GRArc& arc);
  GRArc&
  operator=(const GRArc& arc);
  ~GRArc();

  int
  topo_dim() const
  {
    return 1;
  }
  int
  geom_dim() const
  {
    return n_dim;
  }

  GRCurveGeom::curve_type
  get_curve_type() const
  {
    return GRCurveGeom::ARC;
  }

  //Returns the min and max param for this curve
  double
  min_param() const
  {
    return 0.;
  }
  double
  max_param() const
  {
    return 1.;
  }

  CubitBox
  bounding_box() const;

  double
  param_at_coord(const CubitVector& coord) const;
  double
  param_at_arc_length(const double param_start, const double length) const;
  void
  coord_at_param(const double param, CubitVector& coord) const;

  double
  coord_at_dist(const double dist, const bool start, CubitVector& coord) const;

  double
  coord_at_mid_dist(const double param1, const double param2,
		    CubitVector& middle) const;
  //  double coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
  //			   CubitVector& middle) const;

  void
  closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			 const double param1 = -1.,
			 const double param2 = -1.) const;
  double
  closest_param_on_curve(const CubitVector& coord, const double param1 = -1.,
			 const double param2 = -1.) const;

  bool
  coord_on_curve(const CubitVector& coord, double* const param_X = NULL,
		 double* const param_Y = NULL,
		 double* const param_Z = NULL) const;

  void
  first_deriv(const double param, CubitVector& FD) const;
  //  void first_deriv(const CubitVector& coord, CubitVector& FD) const;

  void
  second_deriv(const double param, CubitVector& SD) const;
  //  void second_deriv(const CubitVector& coord, CubitVector& SD) const;

  double
  curvature(const double param) const;
  //  double curvature(const CubitVector& coord) const;

  void
  curvature(const double, CubitVector&) const
  {
    assert(0);
  }
  //  void curvature(const CubitVector&, CubitVector&) const
  //  { assert(0); }

  void
  unit_normal(const double param, CubitVector& normal) const;
  //  void unit_normal(const CubitVector& coord, CubitVector& normal) const;

  void
  unit_tangent(const double param, CubitVector& tangent) const;
  //  void unit_tangent(const CubitVector& coord, CubitVector& tangent) const;

  double
  arc_length() const;
  double
  arc_length(const double param1, const double param2) const;
  //  double arc_length(const CubitVector& coord1, const CubitVector& coord2) const;

  //  void center_point(CubitVector& mid_point) const;
  void
  center_point(const double param1, const double param2,
	       CubitVector& mid_point) const;
  //  void center_point(const CubitVector& coord1, const CubitVector& coord2,
  //		    CubitVector& mid_point) const;

  //  double TVT() const;
  double
  TVT(const double param1, const double param2) const;
  //  double TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  double
  mid_TVT(const double param1, const double param2) const;
  //  double mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  void
  interior_extrema(DLIList<CubitVector*>& point_list) const;

  bool
  line_intersect(const double coords1[3], const double coords2[3],
		 const double begParam, const double endParam,
		 std::vector<double> & params) const;

  double
  x(const int pt) const;
  double
  y(const int pt) const;
  double
  z(const int pt) const;
  void
  get_coord(const int pt, CubitVector& coord) const;
  const CubitVector&
  get_coord(const int pt) const;

  void
  set_x(const int pt, const double x);
  void
  set_y(const int pt, const double y);
  void
  set_z(const int pt, const double z);
  void
  set_coord(const int pt, const CubitVector& new_coord);

  void
  get_center(CubitVector& center)
  {
    center = points[0];
  }
  double
  get_radius()
  {
    return radius;
  }

  bool
  isValid() const;

};

class GRCubic : public GRCurveGeom {

  //Class for cubic parametric curves. The curves are defined as:
  //x(t) = x_coeff[0] + t * x_coeff[1] + t^2 * x_coeff[2] + t^3 * x_coeff[3]
  //Likewise for y(t).
  //This class only works for 2D.

  double* x_coeff;
  double* y_coeff;

  std::set<double, param_compare> crit_pts;

public:

  GRCubic();
  GRCubic(const double x_coefficient[4], const double y_coefficient[4]);
  GRCubic(const GRCubic& cubic);
  GRCubic&
  operator=(const GRCubic& cubic);
  ~GRCubic();

  void
  set_coeffs(const double xcoeff[4], const double ycoeff[4]);
  void
  set_xcoeffs(const double xcoeff[4]);
  void
  set_ycoeffs(const double ycoeff[4]);

  int
  topo_dim() const
  {
    return 1;
  }
  int
  geom_dim() const
  {
    return n_dim;
  }

  GRCurveGeom::curve_type
  get_curve_type() const
  {
    return GRCurveGeom::CUBIC;
  }

  //Returns the min and max param for this curve
  double
  min_param() const
  {
    return 0.;
  }
  double
  max_param() const
  {
    return 1.;
  }

  CubitBox
  bounding_box() const;

  double
  param_at_coord(const CubitVector& coord) const;
  double
  param_at_arc_length(const double param_start, const double length) const;

  void
  coord_at_param(const double param, CubitVector& coord) const;

  double
  coord_at_dist(const double dist, const bool start, CubitVector& coord) const;

  double
  coord_at_mid_dist(const double param1, const double param2,
		    CubitVector& middle) const;
  //  double coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
  //			   CubitVector& middle) const;

  void
  closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			 const double param1 = -1.,
			 const double param2 = -1.) const;
  double
  closest_param_on_curve(const CubitVector& coord, const double param1 = -1.,
			 const double param2 = -1.) const;

  bool
  coord_on_curve(const CubitVector& coord, double* const param_X = NULL,
		 double* const param_Y = NULL,
		 double* const param_Z = NULL) const;

  void
  first_deriv(const double param, CubitVector& FD) const;
  //  void first_deriv(const CubitVector& coord, CubitVector& FD) const;

  void
  second_deriv(const double param, CubitVector& SD) const;
  //  void second_deriv(const CubitVector& coord, CubitVector& SD) const;

  double
  curvature(const double param) const;
  //  double curvature(const CubitVector& coord) const;

  void
  curvature(const double, CubitVector&) const
  {
    assert(0);
  }
  //  void curvature(const CubitVector&, CubitVector&) const
  //  { assert(0); }

  void
  unit_normal(const double param, CubitVector& normal) const;
  //  void unit_normal(const CubitVector& coord, CubitVector& normal) const;

  void
  unit_tangent(const double param, CubitVector& tangent) const;
  //  void unit_tangent(const CubitVector& coord, CubitVector& tangent) const;

  double
  arc_length() const;
  double
  arc_length(const double param1, const double param2) const;
  //  double arc_length(const CubitVector& coord1, const CubitVector& coord2) const;

  //  void center_point(CubitVector& mid_point) const;
  void
  center_point(const double param1, const double param2,
	       CubitVector& mid_point) const;
  //  void center_point(const CubitVector& coord1, const CubitVector& coord2,
  //		    CubitVector& mid_point) const;

  //  double TVT() const;
  double
  TVT(const double param1, const double param2) const;
  //  double TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  double
  mid_TVT(const double param1, const double param2) const;
  //  double mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  void
  interior_extrema(DLIList<CubitVector*>& point_list) const;

  bool
  line_intersect(const double coords1[3], const double coords2[3],
		 const double begParam, const double endParam,
		 std::vector<double> & params) const;

  double
  x(const int pt) const;
  double
  y(const int pt) const;
  double
  z(const int pt) const;
  void
  get_coord(const int pt, CubitVector& coord) const;
  const CubitVector&
  get_coord(const int pt) const;

  void
  set_x(const int pt, const double x);
  void
  set_y(const int pt, const double y);
  void
  set_z(const int pt, const double z);
  void
  set_coord(const int pt, const CubitVector& new_coord);

  void
  get_center(CubitVector& center)
  {
    center.set(LARGE_DBL, LARGE_DBL, LARGE_DBL);
    vFatalError("This function is only defined for arcs.",
		"GRCubic::get_center(CubitVector&)\n");
  }
  double
  get_radius()
  {
    vFatalError("This function is only defined for arcs.",
		"GRCubic::get_radius()\n");
    return DBL_MAX;
  }

  void
  get_crit_pts(std::set<double>& critical_points,
	       const double add_param = 0.) const
  {

    std::transform(crit_pts.begin(), crit_pts.end(),
		   std::inserter(critical_points, critical_points.begin()),
		   std::bind2nd(std::plus<double>(), add_param));

  }

  const double*
  get_x_coeff() const
  {
    return x_coeff;
  }
  const double*
  get_y_coeff() const
  {
    return y_coeff;
  }

  bool
  isValid() const;

private:

  void
  compute_crit_pts();
  void
  insert_crit_pts(const int num_roots, const double root1, const double root2);

  bool
  newton(const CubitVector& coord, const double lower_param,
	 const double upper_param, double& return_param) const;
  bool
  bisection(const CubitVector& coord, const double lower_param,
	    const double upper_param, double& return_param) const;

};

class GRSpline : public GRCurveGeom {

  //This is the base class for all flavors of cubic parametric spline

public:

  GRSpline();
  virtual
  ~GRSpline();

protected:

  int num_cubics;
  GRCubic** cubics;
  bool closed;

  GRSpline(const GRSpline&); // Never defined
  GRSpline&
  operator=(const GRSpline&); // Never defined

public:

  int
  topo_dim() const
  {
    return 1;
  }
  int
  geom_dim() const
  {
    return n_dim;
  }

  virtual GRCurveGeom::curve_type
  get_curve_type() const = 0;

  //Returns the min and max param for this curve
  virtual double
  min_param() const = 0;
  virtual double
  max_param() const = 0;

  CubitBox
  bounding_box() const;

  double
  param_at_coord(const CubitVector& coord) const;
  double
  param_at_arc_length(const double param_start, const double length) const;

  void
  coord_at_param(const double param, CubitVector& coord) const;

  double
  coord_at_dist(const double dist, const bool start, CubitVector& coord) const;

  double
  coord_at_mid_dist(const double param1, const double param2,
		    CubitVector& middle) const;
  //  double coord_at_mid_dist(const CubitVector& coord1, const CubitVector& coord2,
  //			   CubitVector& middle) const;

  void
  closest_coord_on_curve(const CubitVector& coord, CubitVector& closest,
			 const double param1 = -1.,
			 const double param2 = -1.) const;
  double
  closest_param_on_curve(const CubitVector& coord, const double param1 = -1.,
			 const double param2 = -1.) const;

  bool
  coord_on_curve(const CubitVector& coord, double* const param_X = NULL,
		 double* const param_Y = NULL,
		 double* const param_Z = NULL) const;

  void
  first_deriv(const double param, CubitVector& FD) const;
  //  void first_deriv(const CubitVector& coord, CubitVector& FD) const;

  void
  second_deriv(const double param, CubitVector& SD) const;
  //  void second_deriv(const CubitVector& coord, CubitVector& SD) const;

  double
  curvature(const double param) const;
  //  double curvature(const CubitVector& coord) const;

  void
  curvature(const double, CubitVector&) const
  {
    assert(0);
  }
  //  void curvature(const CubitVector&, CubitVector&) const
  //  { assert(0); }

  void
  unit_normal(const double param, CubitVector& normal) const;
  //  void unit_normal(const CubitVector& coord, CubitVector& normal) const;

  void
  unit_tangent(const double param, CubitVector& tangent) const;
  //  void unit_tangent(const CubitVector& coord, CubitVector& tangent) const;

  double
  arc_length() const;
  double
  arc_length(const double param1, const double param2) const;
  //  double arc_length(const CubitVector& coord1, const CubitVector& coord2) const;

  //  void center_point(CubitVector& mid_point) const;
  void
  center_point(const double param1, const double param2,
	       CubitVector& mid_point) const;
  //  void center_point(const CubitVector& coord1, const CubitVector& coord2,
  //		    CubitVector& mid_point) const;

  //  double TVT() const;
  double
  TVT(const double param1, const double param2) const;
  //  double TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  double
  mid_TVT(const double param1, const double param2) const;
  //  double mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const;

  void
  interior_extrema(DLIList<CubitVector*>& point_list) const;
  bool
  line_intersect(const double coords1[3], const double coords2[3],
		 const double begParam, const double endParam,
		 std::vector<double> & params) const;
  virtual double
  x(const int pt) const = 0;
  virtual double
  y(const int pt) const = 0;
  virtual double
  z(const int pt) const = 0;
  virtual void
  get_coord(const int pt, CubitVector& coord) const = 0;
  virtual const CubitVector&
  get_coord(const int pt) const = 0;

  virtual void
  set_x(const int pt, const double x) = 0;
  virtual void
  set_y(const int pt, const double y) = 0;
  virtual void
  set_z(const int pt, const double z) = 0;
  virtual void
  set_coord(const int pt, const CubitVector& new_coord) = 0;

  void
  get_center(CubitVector& center)
  {
    center.set(LARGE_DBL, LARGE_DBL, LARGE_DBL);
    vFatalError("This function is only defined for arcs.",
		"GRSpline::get_center(CubitVector&)\n");
  }
  double
  get_radius()
  {
    vFatalError("This function is only defined for arcs.",
		"GRSpline::get_radius()\n");
    return DBL_MAX;
  }

  virtual bool
  isValid() const = 0;

  virtual bool
  closed_param()
  {
    return (iFuzzyComp(points[0].x(), points[n_pts - 1].x()) == 0
	&& iFuzzyComp(points[0].y(), points[n_pts - 1].y()) == 0);
  }

};

class GRBezier : public GRSpline {

public:

  GRBezier();
  GRBezier(const CubitVector input_points[4]);
  GRBezier(const GRBezier& B);
  GRBezier&
  operator=(const GRBezier& B);
  ~GRBezier();

  GRCurveGeom::curve_type
  get_curve_type() const
  {
    return GRCurveGeom::BEZIER;
  }
  double
  min_param() const
  {
    return 0.;
  }
  double
  max_param() const
  {
    return 1.;
  }

  double
  x(const int pt) const;
  double
  y(const int pt) const;
  double
  z(const int pt) const;
  void
  get_coord(const int pt, CubitVector& coord) const;
  const CubitVector&
  get_coord(const int pt) const;

  void
  set_x(const int pt, const double x);
  void
  set_y(const int pt, const double y);
  void
  set_z(const int pt, const double z);
  void
  set_coord(const int pt, const CubitVector& new_coord);

  bool
  isValid() const;

private:

  void
  compute_coeffs(double x_coeff[4], double y_coeff[4]);
  void
  set_coeffs();

};

class GRInterp : public GRSpline {

public:

  GRInterp();
  GRInterp(const CubitVector input_pts[], const int num_pts,
	   const bool is_closed = false);
  GRInterp(const GRInterp& interp);
  GRInterp&
  operator=(const GRInterp& interp);
  ~GRInterp();

  GRCurveGeom::curve_type
  get_curve_type() const
  {
    return GRCurveGeom::INTERP;
  }
  double
  min_param() const
  {
    return 0.;
  }
  double
  max_param() const
  {
    return num_cubics;
  }

  double
  x(const int pt) const;
  double
  y(const int pt) const;
  double
  z(const int pt) const;
  void
  get_coord(const int pt, CubitVector& coord) const;
  const CubitVector&
  get_coord(const int pt) const;

  void
  set_x(const int pt, const double x);
  void
  set_y(const int pt, const double y);
  void
  set_z(const int pt, const double z);
  void
  set_coord(const int pt, const CubitVector& new_coord);

  bool
  isValid() const;

private:
  void
  compute_coeffs_notAKnot(double (*x_coeff)[4], double (*y_coeff)[4]);
  void
  compute_coeffs_natural(double (*x_coeff)[4], double (*y_coeff)[4]);
  void
  compute_coeffs_closed(double (*x_coeff)[4], double (*y_coeff)[4]);

};

///INLINE DEFINITIONS FOR BASE CLASS///

//inline void GRCurveGeom::mid_point(CubitVector& mp) const {
//  assert(isValid()); coord_at_param(0.5 * (max_param() + min_param()), mp);
//}
//
//inline void GRCurveGeom::mid_point(const double param1,
//				   const double param2,
//				   CubitVector& mp) const {
//  assert(isValid());
//  assert(iFuzzyComp(param1, min_param()) >= 0 &&
//	 iFuzzyComp(param1, max_param()) <= 0);
//  assert(iFuzzyComp(param2, min_param()) >= 0 &&
//	 iFuzzyComp(param2, max_param()) <= 0);
//  coord_at_param(0.5 * (param1 + param2), mp);
//}
//
//inline void GRCurveGeom::mid_point(const CubitVector& coord1,
//				   const CubitVector& coord2,
//				   CubitVector& mp) const {
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//  coord_at_param(0.5 * (param_at_coord(coord1) + param_at_coord(coord2)), mp);
//}

///INLINE DEFINITIONS FOR GRLine CLASS///

inline void
GRLine::second_deriv(const double, CubitVector& SD) const
{
  SD.set(0., 0., 0.);
}
//inline void GRLine::second_deriv(const CubitVector&,
//				 CubitVector& SD) const
//{ SD.set(0., 0., 0.); }

inline double
GRLine::curvature(const double) const
{
  return 0.;
}
//inline double GRLine::curvature(const CubitVector&) const { return 0.; }

inline void
GRLine::curvature(const double, CubitVector& curv_vector) const
{
  curv_vector.set(0., 0., 0.);
}

//inline void GRLine::curvature(const CubitVector&,
//			       CubitVector& curv_vector) const
//{ curv_vector.set(0., 0., 0.); }

//inline double GRLine::TVT() const { return 0.; }
inline double
GRLine::TVT(const double param1, const double param2) const
{
  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);
  return 0.;
}
//inline double GRLine::TVT(const CubitVector& coord1,
//			  const CubitVector& coord2) const {
//  assert(isValid());
//  assert(coord_on_curve(coord1));
//  assert(coord_on_curve(coord2));
//  return 0.;
//}

inline double
GRLine::mid_TVT(const double param1, const double param2) const
{
  assert(isValid());
  assert(
      iFuzzyComp(param1, min_param()) >= 0
	  && iFuzzyComp(param1, max_param()) <= 0);
  assert(
      iFuzzyComp(param2, min_param()) >= 0
	  && iFuzzyComp(param2, max_param()) <= 0);
  return (0.5 * (param1 + param2));
}

//inline double GRLine::mid_TVT(const CubitVector& coord1, const CubitVector& coord2) const {
//  return ( 0.5 * ( param_at_coord(coord1) + param_at_coord(coord2) ) );
//}

//No interior extremum for a line.
inline void
GRLine::interior_extrema(DLIList<CubitVector*>& point_list) const
{
  point_list.clean_out();
}

inline bool
GRLine::line_intersect(const double coords1[3], const double coords2[3],
		       const double begParam, const double endParam,
		       std::vector<double> & params) const
{
  // taken from wolfram
  CubitVector p1, p2;
  coord_at_param(begParam, p1);
  coord_at_param(endParam, p2);

  CubitVector p3, p4, a, b, c;
  p3.set(coords1);
  p4.set(coords2);
  a = p2 - p1;
  b = p4 - p3;
  c = p3 - p1;

  // if not coplanar, return nothing
  CubitVector colinear = a * b;
  double coplanar = c % colinear;
  if (iFuzzyComp(coplanar, 0.) != 0) {
    return false;
  }
  else {
    // every point is an intersection point, return endpoints
    if (iFuzzyComp(colinear.x(), 0.) == 0 && iFuzzyComp(colinear.y(), 0.) == 0
	&& iFuzzyComp(colinear.z(), 0.) == 0) {
      //need to check
      CubitVector d = c * a;
      if (iFuzzyComp(d.length(), 0.) == 0) {
	params.push_back(0.);
	params.push_back(1.);
	return true;
      }

    }
    else {
      double param1 = ((c * b) % (a * b)) / colinear.length_squared();
      double param2 = ((c * a) % (a * b)) / colinear.length_squared();
      if ((iFuzzyComp(param1, 1.) < 1 && iFuzzyComp(param1, 0.) > -1)
	  && (iFuzzyComp(param2, 1.) < 1 && iFuzzyComp(param2, 0.) > -1)) {
	// we are only interested in param1, which corresponds to the GRLine
	params.push_back(param1);
	return true;
      }
    }
  }
  return false;
}
inline double
GRLine::x(const int pt) const
{
  assert(pt >= 0 && pt < n_pts);
  return points[pt].x();
}
inline double
GRLine::y(const int pt) const
{
  assert(pt >= 0 && pt < n_pts);
  return points[pt].y();
}
inline double
GRLine::z(const int pt) const
{
  assert(pt >= 0 && pt < n_pts);
  return points[pt].z();
}
inline void
GRLine::get_coord(const int pt, CubitVector& coord) const
{
  assert(pt >= 0 && pt < n_pts);
  coord.set(points[pt].x(), points[pt].y(), points[pt].z());
}
inline const CubitVector&
GRLine::get_coord(const int pt) const
{
  assert(pt >= 0 && pt < n_pts);
  return points[pt];
}

inline void
GRLine::set_x(const int pt, const double new_x)
{
  assert(n_pts == 2);
  assert(points != NULL);
  assert(pt >= 0 && pt < n_pts);

  points[pt].x(new_x);

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    vFatalError("Both points of a line are at the same coordinate.",
		"GRLine::set_x(const int, const double)");
  }
}
inline void
GRLine::set_y(const int pt, const double new_y)
{
  assert(n_pts == 2);
  assert(points != NULL);
  assert(pt >= 0 && pt < n_pts);

  points[pt].y(new_y);

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    vFatalError("Both points of a line are at the same coordinate.",
		"GRLine::set_y(const int, const double)");
  }
}
inline void
GRLine::set_z(const int pt, const double new_z)
{
  assert(n_pts == 2);
  assert(points != NULL);
  assert(pt >= 0 && pt < n_pts);

  if (n_dim == 2) {
    assert(0);
  }
  points[pt].z(new_z);

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    vFatalError("Both points of a line are at the same coordinate.",
		"GRLine::vSety(const int, const double)");
  }
}
inline void
GRLine::set_coord(const int pt, const CubitVector& new_coord)
{
  assert(n_pts == 2);
  assert(points != NULL);
  assert(pt >= 0 && pt < n_pts);
  if (n_dim == 2) {
    assert(iFuzzyComp(new_coord.z(), 0.) == 0);
  }

  points[pt].set(new_coord.x(), new_coord.y(), new_coord.z());

  if (iFuzzyComp(points[0].x(), points[1].x()) == 0
      && iFuzzyComp(points[0].y(), points[1].y()) == 0
      && iFuzzyComp(points[0].z(), points[1].z()) == 0) {
    vFatalError("Both points of a line are at the same coordinate.",
		"GRLine::set_coord(const int, const CubitVector&).");
  }
}

#endif 
