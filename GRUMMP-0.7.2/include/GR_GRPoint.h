#ifndef GR_GRPOINT_H
#define GR_GRPOINT_H

#include "GR_config.h"
#include "GR_misc.h"

#include "Point.hpp"         //CGM point class
#include "DLIList.hpp"       //CGM doubly-linked list class
#include "CubitDefines.h"    //CGM definitions
#include "CubitVector.hpp"   //CGM spatial vector
#include "GR_GRAttrib.h"     //GRUMMP attributes
#include "GR_GRAttribSet.h"  //GRUMMP attribute set

class GRCurve;
class GRCoEdge;
class GRLoop;

class GRPoint : public TBPoint {

private:

  CubitVector* loc;        //Point location
  DLIList<Curve*> curves;  //Curves attached to this point
  GRAttribSet attrib_set;  //Attribute set attached to this point

public:

  //Constructors
  GRPoint();
  //GRPoint(const CubitVector& _loc_, const DLIList<Curve*> _curves_);
  GRPoint(const CubitVector& _loc_);
  //GRPoint(const CubitVector* const _loc_, const DLIList<Curve*> _curves_);
  //GRPoint(const CubitVector* const _loc_);
  GRPoint(const double x, const double y, const double z);

private:

  //Copy constructor and operator= private for now.
//  GRPoint(const GRPoint& GRP);
//  GRPoint& operator=(const GRPoint& GRP);

public:

  virtual
  ~GRPoint();

  //Adding a curve attached to this point
  void
  add_curve(Curve* const c)
  {
    curves.append_unique(c);
  }
  int
  num_curves()
  {
    return curves.size();
  }

  //Coordinates queries.
  virtual CubitVector
  coordinates() const;

  //Attribute management. Derived from pure abstract functions.
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
  //Write GRAttribs out to file.
  CubitStatus
  save_attribs(FILE* file_ptr);
  //Read FactAttribs from file. 
  CubitStatus
  restore_attribs(FILE* file_ptr, unsigned int endian);

  //Returns query engine.
  virtual GeometryQueryEngine*
  get_geometry_query_engine() const
  {
    vFatalError("This function is not defined",
		"GRPoint::get_geometry_query_engine()");
    return NULL;
  }

  //At this point, not sure if these should be defined.
  void
  get_parents_virt(DLIList<TopologyBridge*>&)
  {
    assert(0);
  }
  void
  get_children_virt(DLIList<TopologyBridge*>&)
  {
  } //Returns an empty list.

    //Coordinates queries.
  double
  X() const;
  double
  Y() const;
  double
  Z() const;
  void
  get_loc(double* const pt_loc) const;

  //Setting coordinates / moving point.
  void
  set_x(const double xx);
  void
  set_y(const double yy);
  void
  set_z(const double zz);
  void
  set_loc(const double* const _loc_);
  void
  set_loc(const CubitVector* const new_loc);
  void
  move_by(const CubitVector* const delta_loc);

  //Adjacency queries
  void
  get_loops(DLIList<GRLoop*>& l);
  void
  get_coedges(DLIList<GRCoEdge*>& ce);
  void
  get_curves(DLIList<GRCurve*>& c);
  void
  get_points(DLIList<GRPoint*>& p);

private:

  //Check for validity of the instantiated object.
  bool
  qValid() const;

};

//Inline definitions

//Appends an attribute to the geometric entity.
//The name is attached to the underlying solid model entity.
inline void
GRPoint::append_simple_attribute_virt(CubitSimpleAttrib* csa)
{
  attrib_set.append_attribute(csa);
}

//Remove a simple attribute attached to this geometry entity. The name is 
//removed from the underlying BODY this points to.
inline void
GRPoint::remove_simple_attribute_virt(CubitSimpleAttrib* csa)
{
  attrib_set.remove_attribute(csa);
}

//Removes all simple attributes attached to this geometry entity. Also
//moves lingering GTC attributes.
inline void
GRPoint::remove_all_simple_attribute_virt()
{
  attrib_set.remove_all_attributes();
}

//Gets the attributes attached to this geometry entity. The name is 
//attached to the underlying BODY this points to.
inline CubitStatus
GRPoint::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
{
  return attrib_set.get_attributes(csa_list);
}

inline CubitStatus
GRPoint::get_simple_attribute(const CubitString& name,
			      DLIList<CubitSimpleAttrib*>& csa_list)
{
  return attrib_set.get_attributes(name, csa_list);
}

//Write GRAttribs out to file
inline CubitStatus
GRPoint::save_attribs(FILE* f)
{
  return attrib_set.save_attributes(f);
}

//Read GRAttribs from file
inline CubitStatus
GRPoint::restore_attribs(FILE* f, unsigned int endian)
{
  return attrib_set.restore_attributes(f, endian);
}

inline CubitVector
GRPoint::coordinates() const
{
  return *loc;
}

inline double
GRPoint::X() const
{
  assert(qValid());
  return loc->x();
}
inline double
GRPoint::Y() const
{
  assert(qValid());
  return loc->y();
}
inline double
GRPoint::Z() const
{
  assert(qValid());
  return loc->z();
}

inline void
GRPoint::get_loc(double* const pt_loc) const
{
  assert(qValid());
  pt_loc[0] = loc->x();
  pt_loc[1] = loc->y();
  pt_loc[2] = loc->z();
}

inline void
GRPoint::set_x(const double xx)
{
  loc->x(xx);
}
inline void
GRPoint::set_y(const double yy)
{
  loc->y(yy);
}
inline void
GRPoint::set_z(const double zz)
{
  loc->z(zz);
}

inline void
GRPoint::set_loc(const double* const _loc_)
{
  loc->set(_loc_);
}
inline void
GRPoint::set_loc(const CubitVector* const new_loc)
{
  loc->set(*new_loc);
}

inline void
GRPoint::move_by(const CubitVector* const delta_loc)
{
  assert(qValid());
  loc->x(loc->x() + delta_loc->x());
  loc->y(loc->y() + delta_loc->y());
  loc->z(loc->z() + delta_loc->z());
}

#endif
