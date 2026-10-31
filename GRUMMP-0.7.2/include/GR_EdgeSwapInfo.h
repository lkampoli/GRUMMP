#ifndef GR_EdgeSwapInfo_h
#define GR_EdgeSwapInfo_h

#include <set>

#include "GR_config.h"

namespace GRUMMP
{
  class FaceSwapInfo3D;

/// \brief A data carrier for the information needed for edge
/// swapping.
///
/// At creation time, this includes only the edge identity.  Decisions
/// about whether to reconfigure an edge are made by a SwapDecider,
/// which knows about quality measures.  During the decision process,
/// info about the "after" configuration is added.  Then this info is
/// passed to a SwapManager reconfiguration function.
  class EdgeSwapInfo {
  public:
    std::set<Vert*> m_TotalVerts;
    EdgeSwapInfo() :
	m_origBFace0(NULL), m_origBFace1(NULL), m_northVert(NULL), m_southVert(
	NULL), m_numFaces(0), m_bestCanon(0), m_bestPerm(0), m_tooManyFaces(0), m_isBdryEdge(
	    false)
    {
    }
    EdgeSwapInfo(const FaceSwapInfo3D& FC);
    unsigned
    getNumOrigCells() const
    {
      return m_numFaces - (m_isBdryEdge ? 1 : 0);
    }
    EdgeSwapInfo&
    operator=(const EdgeSwapInfo& old)
    {
      for (int i = 0; i < m_maxFaces + 1; i++) {
	m_outerVerts[i] = old.m_outerVerts[i];
	m_origCells[i] = old.m_origCells[i];
      }
      m_origBFace0 = old.m_origBFace0;
      m_origBFace1 = old.m_origBFace1;
      m_northVert = old.m_northVert;
      m_southVert = old.m_southVert;
      m_numFaces = old.m_numFaces;
      m_bestCanon = old.m_bestCanon;
      m_bestPerm = old.m_bestPerm;
      m_tooManyFaces = old.m_tooManyFaces;
      m_isBdryEdge = old.m_isBdryEdge;

      m_TotalVerts.clear();
      for (std::set<Vert*>::iterator it = old.m_TotalVerts.begin();
	  it != old.m_TotalVerts.end(); ++it)
	m_TotalVerts.insert(*it);

      return *this;
    }
    unsigned
    getNumOrigVerts() const
    {
      return m_numFaces;
    }
    Cell *
    getOrigCell(const unsigned u) const
    {
      assert(u < getNumOrigCells());
      return m_origCells[u];
    }
    BFace *
    getBFace(const int i) const
    {
      assert(i == 0 || i == 1);
      return (i == 0 ? m_origBFace0 : m_origBFace1);
    }
    Vert *
    getNorthVert() const
    {
      return m_northVert;
    }
    Vert *
    getSouthVert() const
    {
      return m_southVert;
    }
    Vert *
    getVert(const unsigned u) const
    {
      assert(u < m_numFaces);
      return m_outerVerts[u];
    }
    void
    setNewConfig(const unsigned bestCanon, const unsigned bestPerm)
    {
      m_bestCanon = bestCanon;
      m_bestPerm = bestPerm;
    }
    void
    getNewConfig(unsigned& bestCanon, unsigned& bestPerm) const
    {
      bestCanon = m_bestCanon;
      bestPerm = m_bestPerm;
    }
    bool
    isBdryEdge() const
    {
      return m_isBdryEdge;
    }
    bool
    isSwappableEdge() const
    {
      return !m_tooManyFaces && m_numFaces > 0;
    }
  private:
    // Private data
    enum {
      m_maxFaces = 10
    };
    // The "+1" in these arrays is to allow space for bdry edge swapping
    Vert *m_outerVerts[m_maxFaces + 1];
    Cell *m_origCells[m_maxFaces + 1];
    BFace *m_origBFace0, *m_origBFace1;

    //Vert *m_TotalVerts[m_maxFaces+3];  // "+3" comes from the number of outerVerts[m_maxFaces+1] and the number of two Verts on the poles.

    Vert *m_northVert, *m_southVert;
    // These variables are all small, so packing them in is no problem.
    // Number of faces should never be more than 10.
    // The number of canonical configurations is no larger than 7 for
    // cases currently coded.
    // The number of permutations is no larger than numFaces
    unsigned m_numFaces :8, m_bestCanon :6, m_bestPerm :6;
    // Some status indicators.
    bool m_tooManyFaces :1, m_isBdryEdge :1;
    //EdgeSwapInfo();
    //EdgeSwapInfo(const EdgeSwapInfo&);
//    ~EdgeSwapInfo(){m_TotalVerts.clear()};
//    EdgeSwapInfo& operator=(const EdgeSwapInfo&);
    // These are called (directly or indirectly) from the constructor.
    void
    gatherEdgeInfo(Face* const face, Vert* const otherVert);
    void
    gatherBdryEdgeInfo(Face* const face, Vert* const otherVert);
    void
    gatherIntBdryEdgeInfo(Face* const face, Vert* const otherVert);

  };
} // namespace GRUMMP
#endif
