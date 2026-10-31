#ifndef GR_CURVE_DISCRETIZER_H
#define GR_CURVE_DISCRETIZER_H

#include "GR_config.h"
#include "GR_BaryCentricLengthScale.h"

#include "CubitVector.hpp"
#include "RefEdge.hpp"
#include "FacetSurface.hpp"
#include "FacetEvalTool.hpp"
#include "RefVertex.hpp"
#include <list>
#include <set>

class RefEdge;

class CurveSampleData {

  double m_param;
  CubitVector* m_coord;
  CubitVector* m_tangent;

  CurveSampleData&
  operator=(const CurveSampleData&)
  {
    assert(0);
    return *this;
  }

public:

  CurveSampleData() :
      m_param(-LARGE_DBL), m_coord(NULL), m_tangent(NULL)
  {
  }
  CurveSampleData(const CurveSampleData& CSD) :
      m_param(CSD.m_param), m_coord(
	  CSD.m_coord ? new CubitVector(*CSD.m_coord) : NULL), m_tangent(
	  CSD.m_tangent ? new CubitVector(*CSD.m_tangent) : NULL)
  {
  }
  virtual
  ~CurveSampleData()
  {
    remove_coord();
    remove_tangent();
  }

  void
  set_param(double param)
  {
    m_param = param;
  }
  void
  set_coord(const CubitVector& coord)
  {
    remove_coord();
    m_coord = new CubitVector(coord);
  }
  void
  set_tangent(const CubitVector& tangent)
  {
    remove_tangent();
    m_tangent = new CubitVector(tangent);
  }

  double
  get_param() const
  {
    return m_param;
  }
  CubitVector*
  get_coord() const
  {
    assert(m_coord);
    return m_coord;
  }
  CubitVector*
  get_tangent() const
  {
    assert(m_tangent);
    return m_tangent;
  }

  void
  remove_coord()
  {
    if (m_coord) {
      delete m_coord;
    }
    m_coord = NULL;
  }
  void
  remove_tangent()
  {
    if (m_tangent) {
      delete m_tangent;
    }
    m_tangent = NULL;
  }

};

class CurveSample {

  //curve to sample
  RefEdge* m_curve;

  //stores the sample
  std::list<CurveSampleData> m_curve_sample;

  CurveSample() :
      m_curve(NULL), m_curve_sample()
  {
    assert(0);
  }
  CurveSample(const CurveSample&) :
      m_curve(NULL), m_curve_sample()
  {
    assert(0);
  }
  CurveSample&
  operator=(const CurveSample&)
  {
    assert(0);
    return *this;
  }

public:

  CurveSample(RefEdge* const curve) :
      m_curve(curve), m_curve_sample()
  {
  }
  virtual
  ~CurveSample()
  {
  }

  RefEdge*
  get_curve() const
  {
    assert(m_curve);
    return m_curve;
  }

  void
  init_sample(unsigned int num_points, bool compute_tangent);
  void
  clear_sample()
  {
    m_curve_sample.clear();
  }

  //constructs a sample of num_points equally in PARAMETRIC space
  void
  sample(int num_points = 10)
  {
    init_sample(num_points, false);
  }
  //constructs a sample such that the tvt can be approximated fairly
  //accurately by summing the tangent difference between sampling points.
  void
  sample_tvt(double threshold = 1.e-4);

  //computes total variation of the tangent between params.
  double
  compute_tvt(double beg_param, double end_param = LARGE_DBL);

private:

  struct ParamGreaterThan {
    const double m_param;
    ParamGreaterThan(const double param) :
	m_param(param)
    {
    }

    bool
    operator()(const CurveSampleData& sample_data)
    {
      return (sample_data.get_param() > m_param);
    }
  };

};

/// This class discretizes all refEdges into subsegs,
/// breaking them into CurvePoints.
/// Different Types of Discretizer, only two, one by TVT
/// Other based on lengthscale of underlying intersecting mesh

class CurveDiscretizer {

protected:

  CurveSample m_sample;

  CurveDiscretizer() :
      m_sample(NULL)
  {
    assert(0);
  }
  CurveDiscretizer(const CurveDiscretizer&) :
      m_sample(NULL)
  {
    assert(0);
  }
  CurveDiscretizer&
  operator=(const CurveDiscretizer&)
  {
    assert(0);
    return *this;
  }

public:

  class CurvePoint {

    double m_param;
    CubitVector m_coord;

  public:

    CurvePoint() :
	m_param(-1.), m_coord()
    {
    }
    CurvePoint(double param, const CubitVector& coord) :
	m_param(param), m_coord(coord)
    {
    }
    CurvePoint(const CurvePoint& CP) :
	m_param(CP.m_param), m_coord(CP.m_coord)
    {
    }
    virtual
    ~CurvePoint()
    {
    }

    double
    get_param() const
    {
      return m_param;
    }
    const CubitVector&
    get_coord() const
    {
      return m_coord;
    }

    void
    set_param(double param)
    {
      m_param = param;
    }
    void
    set_coord(const CubitVector& coord)
    {
      m_coord.set(coord);
    }

    bool
    operator<(const CurvePoint& pt) const
    {
      return m_param < pt.m_param;
    }

  private:

    CurvePoint&
    operator=(const CurvePoint&)
    {
      assert(0);
      return *this;
    }

  };

  enum CDType {
    NONE, UNIFORM, TANGENT, CURVATURE, LENGTHSCALE, QUASISTRUCT
  };

  typedef std::set<CurvePoint> DiscreteCurve;

  CurveDiscretizer(RefEdge* const curve) :
      m_sample(curve)
  {
    m_sample.sample_tvt();
  }

  virtual
  ~CurveDiscretizer()
  {
  }

  virtual CDType
  get_type() const = 0;

  virtual void
  discretize(DiscreteCurve& discrete_curve) = 0;

};

//This is to discretize the curves based on the total vartiation of 
//the tangent angle.

class TangentDiscretizer : public CurveDiscretizer {

  double m_tvt_bound;

public:

  TangentDiscretizer(RefEdge* const curve, double tvt_bound) :
      CurveDiscretizer(curve), m_tvt_bound(tvt_bound)
  { /*assert( m_tvt_bound > 0 && m_tvt_bound <= 0.5 * M_PI );*/
  }

  virtual
  ~TangentDiscretizer()
  {
  }

  CDType
  get_type() const
  {
    return TANGENT;
  }

  void
  discretize(DiscreteCurve& discrete_curve);

private:

  TangentDiscretizer() :
      CurveDiscretizer(), m_tvt_bound(-LARGE_DBL)
  {
    assert(0);
  }
  TangentDiscretizer(const TangentDiscretizer&) :
      CurveDiscretizer(), m_tvt_bound(-LARGE_DBL)
  {
    assert(0);
  }
  TangentDiscretizer&
  operator=(const TangentDiscretizer&)
  {
    assert(0);
    return *this;
  }

};

/// Used for surface insertion, samples based on underlying mesh using
/// barycentriclengthscale object.
class LengthScaleDiscretizer : public CurveDiscretizer {

  const VolMesh * m_backgroundMesh;
  BaryLengthScale3D * m_bls;
  std::set<const CurvePoint*> intBdryIntersections;
public:

  LengthScaleDiscretizer(RefEdge* const curve, const VolMesh * backgroundMesh) :
      CurveDiscretizer(curve), m_backgroundMesh(backgroundMesh), m_bls(
	  new BaryLengthScale3D(backgroundMesh)), intBdryIntersections()
  {
  }
  virtual
  ~LengthScaleDiscretizer()
  {
  }

  CDType
  get_type() const
  {
    return LENGTHSCALE;
  }

  void
  discretize(DiscreteCurve& discrete_curve);

  bool
  isIntBdryIntersection(const CurvePoint * cp)
  {
    return intBdryIntersections.count(const_cast<CurvePoint *>(cp));
  }

private:

  LengthScaleDiscretizer() :
      CurveDiscretizer(), m_backgroundMesh(NULL), m_bls(NULL)
  {
    assert(0);
  }
  LengthScaleDiscretizer(const LengthScaleDiscretizer&) :
      CurveDiscretizer(), m_backgroundMesh(NULL), m_bls(NULL)
  {
    assert(0);
  }
  LengthScaleDiscretizer&
  operator=(const LengthScaleDiscretizer&)
  {
    assert(0);
    return *this;
  }

};

class QuasiStructDiscretizer : public CurveDiscretizer {

  const VolMesh * m_backgroundMesh;
  std::set<const CurvePoint*> intBdryIntersections;
public:

  QuasiStructDiscretizer(RefEdge* const curve, const VolMesh * backgroundMesh) :
      CurveDiscretizer(curve), m_backgroundMesh(backgroundMesh), intBdryIntersections()
  {
  }
  virtual
  ~QuasiStructDiscretizer()
  {
  }

  CDType
  get_type() const
  {
    return QUASISTRUCT;
  }

  void
  discretize(DiscreteCurve& discrete_curve);

  bool
  isIntBdryIntersection(const CurvePoint * cp)
  {
    return intBdryIntersections.count(const_cast<CurvePoint *>(cp));
  }

private:

  QuasiStructDiscretizer() :
      CurveDiscretizer(), m_backgroundMesh(NULL)
  {
    assert(0);
  }
  QuasiStructDiscretizer(const QuasiStructDiscretizer&) :
      CurveDiscretizer(), m_backgroundMesh(NULL)
  {
    assert(0);
  }
  QuasiStructDiscretizer&
  operator=(const QuasiStructDiscretizer&)
  {
    assert(0);
    return *this;
  }

};

//Inserts point uniformely on the curve.

//class UniformDiscretizer : public CurveDiscretizer {
//
//  int m_num_points;
//
// public:
//
//  UniformDiscretizer(RefEdge* const curve, int num_points)
//    : CurveDiscretizer(curve), m_num_points(num_points) { assert(m_num_points > 1); }
//  virtual ~UniformDiscretizer() { }
//
//  Type get_type() const { return UNIFORM; }
//
//  void discretize(DiscreteCurve& discrete_curve);
//
// private:
//
//  UniformDiscretizer()
//    : CurveDiscretizer(), m_num_points(-1) { assert(0); }
//  UniformDiscretizer(const UniformDiscretizer&)
//    : CurveDiscretizer(), m_num_points(-1) { assert(0); }
//  UniformDiscretizer& operator=(const UniformDiscretizer&)
//    { assert(0); return *this; }
//
//};
//
//class CurvatureDiscretizer : public CurveDiscretizer {
//
//  double m_ratio;
//
// public:
//
//  CurvatureDiscretizer(RefEdge* const curve, double ratio)
//    : CurveDiscretizer(curve), m_ratio(ratio) { assert(ratio > 0.); }
//  virtual ~CurvatureDiscretizer() {}
//
//  Type get_type() const { return CURVATURE; }
//
//  void discretize(DiscreteCurve& discrete_curve);
//
// private:
//
//  CurvatureDiscretizer()
//    : CurveDiscretizer(), m_ratio(-1.) { assert(0); }
//  CurvatureDiscretizer(const UniformDiscretizer&)
//    : CurveDiscretizer(), m_ratio(-1.) { assert(0); }
//  CurvatureDiscretizer& operator=(const UniformDiscretizer&)
//    { assert(0); return *this; }
//};

#endif
