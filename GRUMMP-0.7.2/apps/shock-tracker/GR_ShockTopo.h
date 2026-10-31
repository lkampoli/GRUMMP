/*
 * GR_ShockTopo.h
 *
 *  Created on: Nov 3, 2017
 *      Author: cfog
 */

#ifndef APPS_SHOCK_TRACKER_GR_SHOCKTOPO_H_
#define APPS_SHOCK_TRACKER_GR_SHOCKTOPO_H_

#include <stddef.h>
#include <assert.h>
#include <vector>

namespace GRUMMP
{
  class ShockEndPoint {
    static constexpr int maxIncCurves = 8;
    int m_curveID[maxIncCurves];
    int m_ID, m_numCurves;
  public:
    ShockEndPoint(const int ID, const int nCurves, const int curves[]);

    int
    getID() const
    {
      return m_ID;
    }

    int
    getNumCurves() const
    {
      return m_numCurves;
    }
    int
    getCurveID(const int which) const
    {
      assert(which >= 0 && which < m_numCurves);
      return m_curveID[which];
    }

    static constexpr int
    getMaxIncCurves()
    {
      return maxIncCurves;
    }
  };

  class ShockCurve {
  public:
    enum InteractionType {
      ShockPerimeter, ShockBdry, Reflection, TripleLine, QuadrupleLine, Invalid
    };
  private:
    static constexpr int maxIncSurfs = 8;
    int m_surfID[maxIncSurfs];
    int m_endpointID[2];
    int m_ID, m_numSurfs;
    InteractionType m_type;
  public:
    ShockCurve(const int ID, const InteractionType type, const int nSurfs,
	       const int surfs[], const int endPoints[]);

    int
    getEndPointID(const int which) const
    {
      assert(which == 0 || which == 1);
      return m_endpointID[which];
    }

    int
    getID() const
    {
      return m_ID;
    }

    int
    getNumSurfaces() const
    {
      return m_numSurfs;
    }

    int
    getSurfaceID(const int which) const
    {
      assert(which >= 0 && which < m_numSurfs);
      return m_surfID[which];
    }

    static constexpr int
    getMaxIncSurfaces()
    {
      return maxIncSurfs;
    }

    static int
    getNumReqSurfs(const InteractionType type);
  };

  class ShockSurface {
    static constexpr int maxBoundingCurves = 20;
    int m_curveID[maxBoundingCurves];
    int m_ID, m_numCurves;
  public:
    ShockSurface(const int ID, const int nCurves, const int curves[]);

    int
    getID() const
    {
      return m_ID;
    }

    int
    getCurveID(const int which) const
    {
      assert(which >= 0 && which < m_numCurves);
      return m_curveID[which];
    }

    int
    getNumCurves() const
    {
      return m_numCurves;
    }

    static constexpr int
    getMaxBoundingCurves()
    {
      return maxBoundingCurves;
    }
  };

  class ShockTopo {
    std::vector<ShockEndPoint> m_endpoints;
    std::vector<ShockCurve> m_curves;
    std::vector<ShockSurface> m_surfs;

    void
    checkReadResult(const int result, const int expected,
		    const size_t lineNumber) const;

  public:
    ShockTopo(const char topoFileName[]);
    size_t
    getNumEndPoints() const
    {
      return m_endpoints.size();
    }
    size_t
    getNumCurves() const
    {
      return m_curves.size();
    }
    size_t
    getNumSurfaces() const
    {
      return m_surfs.size();
    }
    ShockEndPoint
    getShockEndPoint(const int ID) const;
    ShockCurve
    getShockCurve(const int ID) const;
    ShockSurface
    getShockSurface(const int ID) const;
    void
    addShockEndPoint(const ShockEndPoint SEP);
    void
    addShockCurve(const ShockCurve SC);
    void
    addShockSurface(const ShockSurface SS);
    bool
    isValid() const;
  };
}


#endif /* APPS_SHOCK_TRACKER_GR_SHOCKTOPO_H_ */
