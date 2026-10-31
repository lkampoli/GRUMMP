#include "CastTo.hpp"
#include "CoEdgeSM.hpp"
#include "Point.hpp"
#include "Surface.hpp"
#include "TopologyBridge.hpp"

#include "GR_GRLoop.h"
#include "GR_GRCoEdge.h"
#include "GR_GRCurve.h"
#include "GR_GRPoint.h"

GRLoop::GRLoop() :
    surface(NULL), coedge_list()
{
  coedge_list.clean_out();
}

//GRLoop::
//GRLoop(Surface* surf, DLIList<CoEdgeSM*>& ce_list) :
//surface(surf), coedge_list()
//{
//  coedge_list.clean_out();
//  coedge_list += ce_list;
//}

GRLoop::GRLoop(DLIList<CoEdgeSM*>& ce_list) :
    surface(NULL), coedge_list()
{
  coedge_list.clean_out();
  coedge_list += ce_list;
}

GRLoop::~GRLoop()
{
  surface = NULL;
  coedge_list.clean_out();
}

void
GRLoop::get_parents_virt(DLIList<TopologyBridge*>& parents)
{

  parents.append(surface);

}

void
GRLoop::get_children_virt(DLIList<TopologyBridge*>& children)
{

  CAST_LIST_TO_PARENT(coedge_list, children);

}

//void GRLoop::
//get_loops(DLIList<GRLoop*>& loop_list) {
//
//  loop_list.clean_out();
//  loop_list.append(this);
//
//}

void
GRLoop::get_coedges(DLIList<GRCoEdge*>& ce_list)
{

  ce_list.clean_out();
  coedge_list.reset();
  for (int i = 0; i < coedge_list.size(); i++)
    if (GRCoEdge* coedge = dynamic_cast<GRCoEdge*>(coedge_list.next(i)))
      ce_list.append(coedge);

}

//void GRLoop::
//get_curves(DLIList<GRCurve*>& curve_list) {
//
//  curve_list.clean_out();
//  DLIList<GRCoEdge*> ce_list;
//  get_coedges(ce_list);
//  ce_list.reset();
//
//  for (int i = coedge_list.size(); i--; ) {
//    GRCoEdge* coedge = ce_list.get_and_step();
//    GRCurve* curve = dynamic_cast<GRCurve*>(coedge->get_curve());
//    if(curve) curve_list.append_unique(curve);
//  }
//
//}

//void GRLoop::
//get_points(DLIList<GRPoint*>& point_list ) {
//
//  point_list.clean_out();
//  DLIList<GRCurve*> curve_list;
//  get_curves(curve_list);
//  curve_list.reset();
//
//  for (int i = curve_list.size(); i--; ) {
//    GRCurve* curve = curve_list.get_and_step();
//    GRPoint* point = dynamic_cast<GRPoint*>(curve->get_start_point());
//    if (point) point_list.append_unique(point);
//    point = dynamic_cast<GRPoint*>(curve->get_end_point());
//    if (point) point_list.append_unique(point);
//  }
//
//}

//void GRLoop::
//disconnect_all_coedges() {
//
//  coedge_list.reset();
//  for (int i = coedge_list.size(); i--; ) {
//    CoEdgeSM* ce = coedge_list.get_and_step();
//    GRCoEdge* coedge = dynamic_cast<GRCoEdge*>(ce);
//    if (coedge) {
//      assert(coedge->get_loop() == this);
//      coedge->remove_loop();
//    }
//  }
//
//  coedge_list.clean_out();
//
//}

bool
GRLoop::qValid() const
{
  assert(this);
  if (surface == NULL) {
    logMessage(0, "The surface pointer in GRLoop is unitialized.\n");
    return false;
  }
  if (coedge_list.size() == 0) {
    logMessage(0, "The co-edge list in GRLoop is empty.\n");
    return false;
  }

  return true;

}
