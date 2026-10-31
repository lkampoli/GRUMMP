#ifndef GR_GRCOEDGE_H
#define GR_GRCOEDGE_H 1

#include "GR_config.h"
#include "GR_misc.h"
#include "CoEdgeSM.hpp"
#include "CubitDefines.h"
#include "CubitUtil.hpp"
#include "Curve.hpp"
#include "DLIList.hpp"
#include "LoopSM.hpp"
#include "GR_GRAttribSet.h"

class CubitSimpleAttrib;
class GRPoint;
class GRCurve;
class TopologyBridge;

class GRCoEdge : public CoEdgeSM {
  GRCoEdge(const GRCoEdge&);
  GRCoEdge&
  operator=(const GRCoEdge&);

  Curve* curve;         //The curve associated to this co-edge
  LoopSM* loop;         //The loop in which this co-edge is found
  CubitSense ce_sense;  //The sense

public:

  GRCoEdge();
  GRCoEdge(Curve* _curve_, LoopSM* _loop_, CubitSense _sense_);
  GRCoEdge(Curve* _curve_, CubitSense _sense_);

  virtual
  ~GRCoEdge();

  //This is not defined for GRUMMP geometry.
  virtual GeometryQueryEngine*
  get_geometry_query_engine() const;

  //Pure virtual in base class. Attribute management.
  virtual void
  append_simple_attribute_virt(CubitSimpleAttrib* csa);

  virtual void
  remove_simple_attribute_virt(CubitSimpleAttrib* csa);

  virtual void
  remove_all_simple_attribute_virt();

  virtual CubitStatus
  get_simple_attribute(DLIList<CubitSimpleAttrib*>& attrib_list);
  virtual CubitStatus
  get_simple_attribute(const CubitString& name,
		       DLIList<CubitSimpleAttrib*>& attrib_list);

  //Returns the sense of the coedge wrt the underlying edge
  CubitSense
  sense();

  //Adjacency queries
  /*   void get_bodies  ( DLIList<FacetBody   *>& bodies   ); */
  /*   void get_lumps   ( DLIList<FacetLump   *>& lumps    ); */
  /*   void get_shells  ( DLIList<FacetShell  *>& shells   ); */
  /*   void get_surfaces( DLIList<FacetSurface*>& surfaces ); */
  void
  get_loops(DLIList<GRLoop *>& loops);
  void
  get_coedges(DLIList<GRCoEdge *>& coedges);
  void
  get_curves(DLIList<GRCurve *>& curves);
  void
  get_points(DLIList<GRPoint *>& points);

  //Returns parents and children of this topological entity
  void
  get_parents_virt(DLIList<TopologyBridge*>& parents);
  void
  get_children_virt(DLIList<TopologyBridge*>& children);

  //Sets the loop pointer associated with this co-edge.
  void
  add_loop(LoopSM* _loop_);

  //Sets the sense of this co-edge with respect to the face loop.
  void
  set_sense(CubitSense _sense_);

  //Reverses the sense of the co-edge.
  void
  reverse_sense();

  //Returns the curve associated with this co-edge.
  Curve*
  get_curve() const;

  //Returns the loop associated with this co-edge.
  LoopSM*
  get_loop() const;

  //Sets the loop pointer to NULL.
  void
  remove_loop();

  //Sets the curve pointer to NULL.
  void
  remove_curve();

private:

  bool
  qValid() const;

};

inline GeometryQueryEngine*
GRCoEdge::get_geometry_query_engine() const
{
  vFatalError("This function is not supported in GRUMMP",
	      "GRCoEdge::get_geometry_query_engine()");
  return NULL;
}

inline void
GRCoEdge::append_simple_attribute_virt(CubitSimpleAttrib*)
{
  assert(qValid());
  return;
}

inline void
GRCoEdge::remove_simple_attribute_virt(CubitSimpleAttrib*)
{
  assert(qValid());
  return;
}

inline void
GRCoEdge::remove_all_simple_attribute_virt()
{
  assert(qValid());
  return;
}

inline CubitStatus
GRCoEdge::get_simple_attribute(DLIList<CubitSimpleAttrib*>&)
{
  assert(qValid());
  return CUBIT_FAILURE;
}

inline CubitStatus
GRCoEdge::get_simple_attribute(const CubitString&, DLIList<CubitSimpleAttrib*>&)
{
  assert(qValid());
  return CUBIT_FAILURE;
}

inline void
GRCoEdge::get_parents_virt(DLIList<TopologyBridge*>& parents)
{
  assert(qValid());
  parents.append(loop);
}

inline void
GRCoEdge::get_children_virt(DLIList<TopologyBridge*>& children)
{
  assert(qValid());
  children.append(curve);
}

inline CubitSense
GRCoEdge::sense()
{
  return ce_sense;
}

inline void
GRCoEdge::add_loop(LoopSM* _loop_)
{
  loop = _loop_;
}

inline void
GRCoEdge::set_sense(CubitSense _sense_)
{
  ce_sense = _sense_;
}

inline void
GRCoEdge::reverse_sense()
{
  assert(qValid());
  CubitUtil::opposite_sense(ce_sense);
}

inline Curve*
GRCoEdge::get_curve() const
{
  return curve;
}

inline LoopSM*
GRCoEdge::get_loop() const
{
  return loop;
}

inline void
GRCoEdge::remove_loop()
{
  loop = NULL;
}

inline void
GRCoEdge::remove_curve()
{
  curve = NULL;
}

#endif
