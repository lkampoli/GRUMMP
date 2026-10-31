#ifndef GR_GRLOOP_H
#define GR_GRLOOP_H 1

#include "GR_config.h"
#include "GR_misc.h"
#include "CubitDefines.h"
#include "CubitBox.hpp"
#include "LoopSM.hpp"
#include "DLIList.hpp"

class CubitSimpleAttrib;
class GRPoint;
class GRCurve;
class GRCoEdge;
class Surface;
class TopologyBridge;

class GRLoop : public LoopSM {

private:
  GRLoop(const GRLoop&);
  GRLoop&
  operator=(const GRLoop&);
  Surface* surface;
  DLIList<CoEdgeSM*> coedge_list;

public:

  //Default constructor
  GRLoop();
  //Constructor, takes a pointer to the surface looped by this loop
  //and a list of co_edge forming this loop.
//  GRLoop(Surface* surf, DLIList<CoEdgeSM*>& ce_list);
  //Takes a list of co_edge forming this loop.
  GRLoop(DLIList<CoEdgeSM*>& ce_list);
  //Destructor
  virtual
  ~GRLoop();

  LoopType
  loop_type()
  {
    return LOOP_TYPE_UNKNOWN;
  }
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

  //Bounding box function
  virtual CubitBox
  bounding_box() const;
  //Not yet implemented in GRUMMP.
  virtual GeometryQueryEngine*
  get_geometry_query_engine() const;

  //Adjacency queries.
  //void get_bodies  (DLIList<GRBody   *>& bodies   );
  //void get_lumps   (DLIList<GRLump   *>& lumps    );
  //void get_shells  (DLIList<GRShell  *>& shells   );
  //void get_surfaces(DLIList<GRSurface*>& surfaces );
//  void get_loops   (DLIList<GRLoop   *>& loops   );
  void
  get_coedges(DLIList<GRCoEdge *>& coedges);
//  void get_curves  (DLIList<GRCurve  *>& curves  );
//  void get_points  (DLIList<GRPoint  *>& points  );

  // These two are inherited virtual functions and so must remain.
  void
  get_parents_virt(DLIList<TopologyBridge*>& parents);
  void
  get_children_virt(DLIList<TopologyBridge*>& children);

  //Assigns the surface point to _surface_
  void
  add_surface(Surface* _surface_);

  //Add co-edges to the loop.
  void
  add_coedges(DLIList<CoEdgeSM*>& ce_list);

  //Reverse the sense of loop.
  void
  reverse();

  //Returns to pointer to the surface
  Surface*
  get_surface() const;

  //Nullifies the pointer to the surface
  void
  remove_surface();

//  void disconnect_all_coedges();

  bool
  qValid() const;

};

//Inline definitions

inline void
GRLoop::append_simple_attribute_virt(CubitSimpleAttrib*)
{
  assert(0);
  return;
}

inline void
GRLoop::remove_simple_attribute_virt(CubitSimpleAttrib*)
{
  assert(0);
  return;
}

inline void
GRLoop::remove_all_simple_attribute_virt()
{
  assert(0);
  return;
}

inline CubitStatus
GRLoop::get_simple_attribute(DLIList<CubitSimpleAttrib*>&)
{
  assert(0);
  return CUBIT_FAILURE;
}

inline CubitStatus
GRLoop::get_simple_attribute(const CubitString&, DLIList<CubitSimpleAttrib*>&)
{
  assert(0);
  return CUBIT_FAILURE;
}

inline CubitBox
GRLoop::bounding_box() const
{
  assert(0);
  CubitBox box;
  return box;
}

inline GeometryQueryEngine*
GRLoop::get_geometry_query_engine() const
{
  vFatalError("This function is not supported in GRUMMP",
	      "GRLoop::get_geometry_query_engine()");
  return static_cast<GeometryQueryEngine*>(NULL);
}

inline void
GRLoop::add_surface(Surface *_surface_)
{
  surface = _surface_;
}

inline void
GRLoop::add_coedges(DLIList<CoEdgeSM*>& ce_list)
{
  if (coedge_list.size() == 0)
    coedge_list += ce_list;
  else
    assert(0);
}

inline void
GRLoop::reverse()
{
  coedge_list.reverse();
}

inline Surface*
GRLoop::get_surface() const
{
  return surface;
}

inline void
GRLoop::remove_surface()
{
  surface = NULL;
}

#endif
