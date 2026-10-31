#include "GR_config.h"
#include "CubitSimpleAttrib.hpp"
#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "GR_GRCoEdge.h"
#include "GR_GRLoop.h"
#include "GR_misc.h"
#include "TopologyBridge.hpp"

GRCoEdge::GRCoEdge() :
    curve(NULL), loop(NULL), ce_sense(CUBIT_UNKNOWN)
{
}

GRCoEdge::GRCoEdge(Curve* _curve_, LoopSM* _loop_, CubitSense _sense_) :
    curve(_curve_), loop(_loop_), ce_sense(_sense_)
{
}

GRCoEdge::GRCoEdge(Curve* _curve_, CubitSense _sense_) :
    curve(_curve_), loop(NULL), ce_sense(_sense_)
{
}

GRCoEdge::~GRCoEdge()
{

  remove_loop();
  remove_curve();

}

void
GRCoEdge::get_loops(DLIList<GRLoop*>& loop_list)
{
  if (GRLoop* _loop_ = dynamic_cast<GRLoop*>(loop))
    loop_list.append_unique(_loop_);
}

void
GRCoEdge::get_coedges(DLIList<GRCoEdge*>& coedge_list)
{
  coedge_list.append(this);
}

void
GRCoEdge::get_curves(DLIList<GRCurve*>& curve_list)
{
  if (GRCurve* _curve_ = dynamic_cast<GRCurve*>(curve))
    curve_list.append(_curve_);
}

void
GRCoEdge::get_points(DLIList<GRPoint*>& point_list)
{
  if (GRCurve* _curve_ = dynamic_cast<GRCurve*>(curve))
    _curve_->get_points(point_list);
}

bool
GRCoEdge::qValid() const
{

  if (curve == NULL) {
    logMessage(0, "ERROR: The curve pointer is not set in this GRCoEdge.");
    return false;
  }
  if (loop == NULL) {
    logMessage(0, "ERROR: The loop pointer is not set in this GRCoEdge.");
    return false;
  }
  if (ce_sense == CUBIT_UNKNOWN) {
    logMessage(0, "ERROR: The edge sense is not set in this GRCoEdge.");
    return false;
  }

  return true;

}
