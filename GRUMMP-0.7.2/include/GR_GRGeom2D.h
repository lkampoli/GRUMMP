#ifndef GR_GRGEOM2D_H
#define GR_GRGEOM2D_H 1

/* 2D GRUMMP geometry class: */
/* Reads boundary data from a two dimensional boundary file */
/* then initializes and stores the geometry and topology accordingly. */

#include "GR_config.h"
#include "GR_misc.h"
#include "CubitString.hpp"
#include <vector>
#include <list>
#include <map>

class CubitVector;
class GRPoint;
class GRCurve;
class GRCoEdge;
class GRLoop;

class GRGeom2D {

  //To know whether the geometry is initialized of not.
  bool geom_initialized;

public:

  GRGeom2D();
  GRGeom2D(const char* const bdry_file);
  ~GRGeom2D();

  //Functions giving access to geometry entity
  void
  get_points(std::list<GRPoint*>& points) const;
  void
  get_curves(std::list<GRCurve*>& curves) const;

  GR_index_t
  num_points() const;
  GR_index_t
  num_curves() const;

  //Builds a point at coords, returns pointer to that GRPoint 
  GRPoint*
  build_point(double x_coord, double y_coord);

  //Builds a polyline between two GRPoint's, returns curve pointer.
  GRCurve*
  build_line(GRPoint* const point1, GRPoint* const point2, int region_left,
	     int region_right, int bound_left, int bound_right);

  // offset the curve identified by curveNumber by the coordinate offset
  std::pair<GRCurve *, GRCurve *>
  offset_curve_linear(GR_index_t curveNumber, const CubitVector & offset);
  std::map<GRCurve *, GRCurve *>
  offset_curves(std::list<GRCurve*> & curvesToOffset, double dist);

  bool
  initialized() const
  {
    return geom_initialized;
  }
  void
  set_initialized()
  {
    geom_initialized = true;
  }
  void
  set_uninitialized()
  {
    geom_initialized = false;
  }

private:

  //Member variables containing the number of points and curves.
  GR_index_t n_points, n_curves, n_constructed_curves;

  //Geometric entities stored in these lists
  std::list<GRPoint*> point_list;
  std::list<GRCurve*> curve_list;
  std::list<GRCoEdge*> coedge_list;
  std::list<GRLoop*> loop_list;

  //Struct containing curve info extracted from file.
  struct curve_info {
    CubitString curve_type;
    CubitString l_cond;
    int left;
    CubitString r_cond;
    int right;
    std::vector<double> curve_data;
    curve_info() :
	curve_type(), l_cond(), left(-1), r_cond(), right(-1), curve_data()
    {
    }
  };

  //Vector containing the boundary point coordinates
  std::vector<CubitVector> v_points;
  //Vector containing the string defining every curves.
  std::vector<curve_info> v_curves;

  ///////////////////////
  ///PARSING FUNCTIONS///
  ///////////////////////
  //Opens the boundary file and returns a point to this file.
  FILE*
  open_bdry_file(const char* const bdry_file) const;
  //Reads and initialized the number of points and curves.
  void
  read_num_entities(FILE* const f);
  //Reads the point coords and stores into a vector
  void
  read_points(FILE* const f);
  //Reads the curve info and stores the string into a vector
  void
  read_curves(FILE* const f);
  //Checks patch name validity by comparing strings
  bool
  valid_patch_type(const char* const str) const;
  //Special line reader similar to vGetLineOrAbort()
  void
  get_line(char buffer[], const int size, FILE *f) const;

  ////////////////////////
  ///BUILDING FUNCTIONS///
  ////////////////////////
  //Once the bdry data has been read in the constructor, init geometry
  //by calling this function.
  void
  build_geometry();

  //Builds one curve (curve type built depends on ci.curve_type)
  void
  build_curve(const GR_index_t curve_id,
	      std::map<GR_index_t, GRPoint*>& created_points,
	      std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);
  //Builds a vertex "curve" (a vertex embeded in a face)
  void
  build_vertex(const curve_info* const ci,
	       std::map<GR_index_t, GRPoint*>& created_points,
	       std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);
  //Builds a polyline (or a single line)
  void
  build_polyline(const curve_info* const ci,
		 std::map<GR_index_t, GRPoint*>& created_points,
		 std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);
  //Builds an arc
  void
  build_arc(const curve_info* const ci,
	    std::map<GR_index_t, GRPoint*>& created_points,
	    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves,
	    const bool short_arc = true);

  //Builds a circle
  void
  build_circle(const curve_info* const ci,
	       std::map<GR_index_t, GRPoint*>& created_points,
	       std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);
  //Builds a bezier curve
  void
  build_bezier(const curve_info* const ci,
	       std::map<GR_index_t, GRPoint*>& created_points,
	       std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);
  //Builds a cubic interpolated spline
  void
  build_interp(const curve_info* const ci,
	       std::map<GR_index_t, GRPoint*>& created_points,
	       std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);

  //Build a co-edge given a region id.
  void
  build_coedge(int region, GRCurve* const curve,
	       std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves);

  //Build the loops for a "surface" given the co-edges of this "surface"
  void
  build_loop(std::vector<GRCoEdge*>& coedge_vector,
	     std::vector<GRLoop*>& created_loops);

  ///////////////////
  ///VERIFICATION////
  ///////////////////

  //Verifies that the regions were correctly specified in the boundary file.
  void
  verify_region_conditions(int region, std::vector<GRLoop*>& created_loops);
public:
  void
  writeGeometryVTK(const char strBaseFileName[]);

};

inline void
GRGeom2D::get_points(std::list<GRPoint*>& points) const
{
  points = point_list;
}

inline void
GRGeom2D::get_curves(std::list<GRCurve*>& curves) const
{
  curves = curve_list;
}

inline GR_index_t
GRGeom2D::num_points() const
{
  return static_cast<int>(point_list.size());
}

inline GR_index_t
GRGeom2D::num_curves() const
{
  return static_cast<int>(curve_list.size());
}

#endif
