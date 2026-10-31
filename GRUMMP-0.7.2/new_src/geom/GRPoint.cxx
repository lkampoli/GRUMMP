#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "GR_GRCoEdge.h"
#include "GR_GRLoop.h"

GRPoint::GRPoint() :
    loc(NULL), curves(), attrib_set()
{
  loc = new CubitVector(-LARGE_DBL, -LARGE_DBL, -LARGE_DBL);
  curves.clean_out();
  assert(this->qValid());
}

//GRPoint::
//GRPoint(const CubitVector& _loc_, const DLIList<Curve*> _curves_) :
//  loc(NULL), curves(), attrib_set()
//{
//  loc = new CubitVector(_loc_);
//  curves.clean_out();
//  curves += _curves_;
//  assert(qValid());
//}

GRPoint::GRPoint(const CubitVector& _loc_) :
    loc(NULL), curves(), attrib_set()
{
  loc = new CubitVector(_loc_);
  curves.clean_out();
  assert(this->qValid());
}

//GRPoint::
//GRPoint(const CubitVector* const _loc_, const DLIList<Curve*> _curves_) :
//  loc(NULL), curves(), attrib_set()
//{
//  loc = new CubitVector(*_loc_);
//  curves.clean_out();
//  curves += _curves_;
//  assert(qValid());
//}
//
//GRPoint::
//GRPoint(const CubitVector* const _loc_) :
//  loc(NULL), curves(), attrib_set()
//{
//  loc = new CubitVector(*_loc_);
//  curves.clean_out();
//}
//
GRPoint::GRPoint(double x, double y, double z) :
    loc(NULL), curves(), attrib_set()
{
  loc = new CubitVector(x, y, z);
}

//GRPoint::
//GRPoint(const GRPoint& grp) :
//  Point(), loc(NULL), curves(), attrib_set()
//{
//  assert(grp.qValid());
//  loc = new CubitVector(*(grp.loc));
//  curves.clean_out();
//  curves = grp.curves;
//  assert(qValid());
//}
//
//GRPoint& GRPoint::
//operator=(const GRPoint& grp) {
//  assert(grp.qValid());
//  loc->x(grp.loc->x());
//  loc->y(grp.loc->y());
//  loc->z(grp.loc->z());
//  curves.clean_out();
//  curves = grp.curves;
//  assert(qValid());
//  return *this;
//}

GRPoint::~GRPoint()
{
  remove_all_simple_attribute_virt();
  if (loc) {
    delete loc;
  }
  curves.clean_out();
}

void
GRPoint::get_loops(DLIList<GRLoop*>& l)
{
  DLIList<GRCoEdge*> coedge_list;
  get_coedges(coedge_list);
  coedge_list.reset();

  for (int i = coedge_list.size(); i--;) {
    GRCoEdge* coedge = coedge_list.get_and_step();
    GRLoop* loop = dynamic_cast<GRLoop*>(coedge->get_loop());
    if (loop)
      l.append_unique(loop);
  }
}

void
GRPoint::get_coedges(DLIList<GRCoEdge*>& ce)
{
  DLIList<GRCurve*> curve_list;
  get_curves(curve_list);
  curve_list.reset();
  for (int i = curve_list.size(); i--;)
    curve_list.get_and_step()->get_coedges(ce);
}

void
GRPoint::get_curves(DLIList<GRCurve*>& c)
{
  curves.reset();
  for (int i = 0; i < curves.size(); i++)
    if (GRCurve* curve = dynamic_cast<GRCurve*>(curves.next(i)))
      c.append(curve);
}

void
GRPoint::get_points(DLIList<GRPoint*>& p)
{
  p.append(this);
}

bool
GRPoint::qValid() const
{
  if (!loc) {
    return false;
  }
  if (loc->x() < -1.e200) {
    return false;
  }
  if (loc->y() < -1.e200) {
    return false;
  }
  if (loc->z() < -1.e200) {
    return false;
  }
  return true;
}
