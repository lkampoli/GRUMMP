#ifndef GR_SwapDecider_h
#define GR_SwapDecider_h

#include <string>
#include <vector>

#include "GR_BFace.h"
#include "GR_EdgeSwapInfo.h"
#include "GR_FaceSwapInfo.h"
#include "GR_QualMeasure.h"
#include "GR_SwapManager.h"
#include "GR_Vertex.h"
#include "GR_Aniso.h"
#include "GR_Aniso.h"

#ifdef HAVE_MESQUITE
#include "iMesh_GRUMMP.hh"
#include "Mesquite.hpp"
#include "TMetric.hpp"
#include "AWMetric.hpp"
#include "TMPQualityMetric.hpp"
#include "TargetCalculator.hpp"
#include "RefMeshTargetCalculator.hpp"
#include "ReferenceMesh.hpp"
#include "MeshImpl.hpp"
#endif

#undef SIM_ANNEAL_TEST

namespace GRUMMP
{
  class SwapDecider2D {
  protected:
    QualMeasure *m_pQM;
    std::string m_name;
  private:
    SwapDecider2D(const SwapDecider2D& SD) :
	m_pQM(SD.m_pQM), m_name(SD.m_name)
    {
    }
    SwapDecider2D&
    operator=(const SwapDecider2D&);
  public:
    SwapDecider2D(QualMeasure* const pQM, const std::string name) :
	m_pQM(pQM), m_name(name)
    {
    }
    virtual
    ~SwapDecider2D()
    {
    }
    // No protected data accessors at this point.  Add them later if
    // they're needed.
    const std::string&
    getName() const
    {
      return m_name;
    }
    virtual bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };

  class DelaunaySwapDecider2D : public SwapDecider2D {
  private:
  public:
    DelaunaySwapDecider2D() :
	SwapDecider2D(NULL, "Delaunay")
    {
    }
    ~DelaunaySwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };

  class MinMaxAngleSwapDecider2D : public SwapDecider2D {
  private:
  public:
    MinMaxAngleSwapDecider2D() :
	SwapDecider2D(new Angle2D(), "minmax angle")
    {
    }
    ~MinMaxAngleSwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };

  class DelaunaySwapDeciderSurf : public SwapDecider2D {
  private:
  public:
    DelaunaySwapDeciderSurf() :
	SwapDecider2D(NULL, "Delaunay")
    {
    }
    ~DelaunaySwapDeciderSurf()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };
//  class MinMaxAngleSwapDeciderSurf : public SwapDecider2D {
//  private:
//  public:
//    MinMaxAngleSwapDeciderSurf() :
//      SwapDecider2D(NULL,
//		    "minmax angle") {}
//    ~MinMaxAngleSwapDeciderSurf() {}
//    bool doFaceSwap(const FaceSwapInfo2D& FC) const;
//  };
  class UniformDegreeSwapDecider2D : public SwapDecider2D {
  private:
  public:
    UniformDegreeSwapDecider2D() :
	SwapDecider2D(NULL, "uniform degree")
    {
    }
    ~UniformDegreeSwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };

  class AnisoSwapDecider2D : public SwapDecider2D {
  private:
  public:
    AnisoSwapDecider2D() :
	SwapDecider2D(NULL, "Aniso")
    {
    }
    ~AnisoSwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
  };

  class NullSwapDecider2D : public SwapDecider2D {
  private:
  public:
    NullSwapDecider2D() :
	SwapDecider2D(NULL, "No-swap")
    {
    }
    ~NullSwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D&) const
    {
      return false;
    }
  };

  /*
   * Currently only working for Min Max Angle in metric space. Other criteria were tested in metric space,
   * but best result were obtain with Min Max angle.
   */
  class MetricSwapDecider2D : public SwapDecider2D {
  private:
    bool bMetric;
  public:
    MetricSwapDecider2D(GRUMMP::QualMeasure* pQM2D, bool bM = true) :
	SwapDecider2D(pQM2D, "Metric2D")
    {
      bMetric = bM;
    }
    ~MetricSwapDecider2D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo2D& FC) const;
    bool
    isMetricOK(Vert *pV) const;

  };

  class SwapDecider3D {
  protected:
    bool m_allowBdryChanges, m_allowEdgeSwaps;
    double m_maxBdryAngle;
    QualMeasure *m_pQM;
    std::string m_name;
  private:
    SwapDecider3D(const SwapDecider3D& SD) :
	m_allowBdryChanges(SD.m_allowBdryChanges), m_allowEdgeSwaps(
	    SD.m_allowEdgeSwaps), m_maxBdryAngle(SD.m_maxBdryAngle), m_pQM(
	    SD.m_pQM), m_name(SD.m_name)
    {
    }
    SwapDecider3D&
    operator=(const SwapDecider3D&);
    double
    calcOrigQuality(const EdgeSwapInfo& ES) const;
    bool
    areBdryFacesOKToSwap(const BFace* const pBF0,
			 const BFace* const pBF1) const;
  public:
    SwapDecider3D(const bool bdryChangeable, const bool edgesSwappable,
		  const bool maxBdryAngle, QualMeasure* const pQM,
		  const std::string name) :
	m_allowBdryChanges(bdryChangeable), m_allowEdgeSwaps(edgesSwappable), m_maxBdryAngle(
	    maxBdryAngle), m_pQM(pQM), m_name(name)
    {
    }
    virtual
    ~SwapDecider3D()
    {
    }
    // No protected data accessors at this point.  Add them later if
    // they're needed.
    const std::string&
    getName() const
    {
      return m_name;
    }
    virtual bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
    bool
    doEdgeSwap(EdgeSwapInfo& ES, Vert * const pVWant0 = pVInvalidVert,
	       Vert * const pVWant1 = pVInvalidVert, bool qWant = true) const;
    bool
    doBdryEdgeSwap(EdgeSwapInfo& ES) const;
    double
    evalTetQual(const Cell* const pC) const;
  };

  class NullSwapDecider3D : public SwapDecider3D {
  private:
  public:
    NullSwapDecider3D(const bool bdryChangeable, const double maxBdryAngle = 0) :
	SwapDecider3D(bdryChangeable, false, maxBdryAngle, NULL, "No-swap")
    {
    }
    ~NullSwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D&) const
    {
      return false;
    }
  };

  class DelaunaySwapDecider3D : public SwapDecider3D {
  private:
  public:
    DelaunaySwapDecider3D(const bool bdryChangeable, const double maxBdryAngle =
			      0) :
	SwapDecider3D(bdryChangeable, false, maxBdryAngle, NULL, "Delaunay")
    {
    }
    ~DelaunaySwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
  };

  class MaxMinSineSwapDecider3D : public SwapDecider3D {
  private:
  public:
    MaxMinSineSwapDecider3D(const bool bdryChangeable,
			    const double maxBdryAngle = 0) :
	SwapDecider3D(bdryChangeable, true, maxBdryAngle,
		      new SineDihedralAngles(), "maxmin sine dihedral")
    {
    }
    ~MaxMinSineSwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
  };

  class MinMaxDihedSwapDecider3D : public SwapDecider3D {
  private:
  public:
    MinMaxDihedSwapDecider3D(const bool bdryChangeable,
			     const double maxBdryAngle = 0) :
	SwapDecider3D(bdryChangeable, true, maxBdryAngle, new DihedralAngles(),
		      "minmax dihedral")
    {
    }
    ~MinMaxDihedSwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
  };

  class UniformDegreeSwapDecider3D : public SwapDecider3D {
  private:
  public:
    UniformDegreeSwapDecider3D(const bool bdryChangeable,
			       const double maxBdryAngle = 0) :
	SwapDecider3D(bdryChangeable, false, maxBdryAngle, NULL,
		      "uniform degree")
    {
    }
    ~UniformDegreeSwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
  };

  /*
   * Class intended to be used only for swapping in metric space, not all quality measures are able to swap
   * in metric space at the moment
   */
  class MetricSwapDecider3D : public SwapDecider3D {
  private:
  public:
    MetricSwapDecider3D(const bool bdryChangeable,
			const double maxBdryAngle = 0,
			GRUMMP::QualMeasure* pQM3D = NULL) :
	SwapDecider3D(bdryChangeable, true, maxBdryAngle, pQM3D, "Metric3D")
    {
    }
    ~MetricSwapDecider3D()
    {
    }
    bool
    doFaceSwap(const FaceSwapInfo3D& FC) const;
  };

} // namespace GRUMMP

#endif
