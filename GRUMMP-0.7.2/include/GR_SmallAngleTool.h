#ifndef GR_SMALL_ANGLE_TOOL
#define GR_SMALL_ANGLE_TOOL

#include "GR_config.h"
#include <map>
#include <set>
#include <vector>

class Vert;
class Subseg;

//This class is basically a wrapper around an STL
//multimap. It carries operationd required to handle
//small input angle in 3D. Implemented as a singleton
//class.

class SmallAngleTool {

private:

  static SmallAngleTool* _instance_;

  typedef std::multimap<Vert*, Subseg*> SmallAngleMap;
  typedef SmallAngleMap::iterator ItMap;
  typedef SmallAngleMap::const_iterator ItcMap;

  SmallAngleMap vert_map;
  double minimum_angle;   //The minimum acceptable angle (degrees)

  //No constructor can be called from the outside
  SmallAngleTool();
  SmallAngleTool(const SmallAngleTool&) :
      vert_map(), minimum_angle(90.)
  {
    assert(0);
  }
  SmallAngleTool&
  operator=(const SmallAngleTool&)
  {
    assert(0);
    return *this;
  }

  //Computes the angle between two subsegs joined at vertex.
  double
  compute_angle(const Vert* const vertex, const Subseg* const subseg1,
		const Subseg* const subseg2) const;

  //This is only called by the initialize function.
  void
  remove_acceptable_angle_entries();

public:

  //Destructor
  ~SmallAngleTool();

  static SmallAngleTool*
  instance();

  //Clears the contents of the vert_map.
  void
  clear_small_angle_data()
  {
    vert_map.clear();
  }

  //Changes the minimum acceptable angle (default = 90 deg.) 
  void
  set_minimum_angle(const double angle);

  //Initializes the small angle data. Takes a vector
  //containing all the subsegments as argument.
  void
  initialize(const std::vector<Subseg*>& all_subsegs);

  //Returns all the subsegments that can be found twice in the map
  //(subsegs bounded by two apex verts)
  void
  get_multiple_entries(std::set<Subseg*>& multiple_entries) const;

  //Returns the number of subsegments currently in a small angle complex.
  int
  size() const
  {
    return static_cast<int>(vert_map.size());
  }

  //Tries to add a subsegment to vert_map.
  //If check_angle is true, will check tp make sure the new subseg
  //is part of a small angle complex. Obviously, this slows down the
  //addition process.
  void
  add_subseg(Subseg* const subseg, const bool check_angle = false);

  //Tries to remove a subsegment from vert_map
  void
  remove_subseg(Subseg* const subseg);

  //Populates the vector subsegments with all the subsegments
  //part of the small angle complex attached to vertex.
  void
  get_subsegs_of_vert(Vert* const vertex,
		      std::vector<Subseg*>& subsegments) const;

  //Populates all_subsegs with all the subsegments in all small angle complex (complices?)
  void
  get_all_subsegs(std::vector<Subseg*>& subsegments) const;

  //Computes the split radius for a collection of subsegments to split.
  double
  compute_split_radius(const std::vector<Subseg*>& subsegs_to_split) const;

  //Checks if vertex is the base of a small angle complex
  bool
  vertex_has_small_angle(Vert* const vertex) const;

  //Checks if subseg is in a small angle complex
  bool
  subseg_has_small_angle(const Subseg* const subseg) const;

  //Prints the address in the map to screen. Just for
  //debugging purposes.
  void
  print_map_to_screen() const;

};

#endif
