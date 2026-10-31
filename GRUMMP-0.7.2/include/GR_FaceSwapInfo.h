#ifndef GR_FaceSwapInfo_h
#define GR_FaceSwapInfo_h

#include "GR_config.h"
#include "GR_Vertex.h"
namespace GRUMMP
{
  class FaceSwapInfo2D {
  private:
    Face *m_pF;
    TriCell* m_pTris[2];
    Vert* m_pVerts[4];
    bool m_nonSwappable, m_orientOK;
//#pragma omp threadprivate (m_pF, m_pTris, m_pVerts, m_nonSwappable, m_orienrOK)
  public:
    FaceSwapInfo2D() :
	m_pF(pFInvalidFace), m_nonSwappable(false), m_orientOK(false)
    {
      m_pTris[0] = m_pTris[1] = NULL;
      m_pVerts[0] = m_pVerts[1] = m_pVerts[2] = m_pVerts[3] = pVInvalidVert;
    }
    virtual
    ~FaceSwapInfo2D()
    {
    }
    ;
    void
    calculateFaceInfo(Face* pF);
    bool
    isSwappable() const
    {
      return (!m_nonSwappable && m_orientOK);
    }

    Face*
    getFace() const
    {
      return m_pF;
    }
    Vert*
    getVertA() const
    {
      return m_pVerts[0];
    }
    Vert*
    getVertB() const
    {
      return m_pVerts[1];
    }
    Vert*
    getVertC() const
    {
      return m_pVerts[2];
    }
    Vert*
    getVertD() const
    {
      return m_pVerts[3];
    }
    Vert*
    getVert(GR_index_t i) const
    {
      assert(i < 4);
      return m_pVerts[i];
    }
    void
    getTris(Cell* pTris[2]) const;
    inline void
    setTris(Face* pF);
  };

  class FaceSwapInfo3D {
  public:
    /// Detailed categorization of faces into swappable and unswappable.
    ///
    /// This taxonomy uses Joe's, plus adds some cases for boundaries,
    /// plus a catch-all that is used for definitely unswappable cases
    /// without finishing all the work to categorize them completely.
    /// Because who really cares if you can't swap the face anyway?
    enum faceCat {
      eT23 = 0,
      eT32 = 1,
      eT22 = 2,
      eT44 = 3,
      eN32 = 4,
      eN44 = 5,
      eN40 = 6,
      eN30 = 7,
      eN20 = 8,
      eBdry = 9,
      eBdryReflex = 10,
      eOther = 11
    };

  private:
    double m_maxAngle;
    Face *m_pF;
    TetCell* m_pTets[4];
    TriBFace* m_pBFaces[2];
    Vert* m_pVerts[6];
    faceCat m_fc;
    int m_nTets;
    // For this class, default copy construction and operator= actually
    // give sensible results.
  public:
    FaceSwapInfo3D(double maxSurfaceAngle = 0) :
	m_maxAngle(maxSurfaceAngle), m_pF(NULL), m_fc(eOther), m_nTets(0)
    {
    }

    /// A constructor used to force edge swap consideration
    FaceSwapInfo3D(Face* face, Vert* north, Vert* south, Vert* other);

    FaceSwapInfo3D &
    operator=(const FaceSwapInfo3D &old)
    {
      m_maxAngle = old.m_maxAngle;
      m_pF = old.m_pF;
      for (int i = 0; i < 4; i++)
	m_pTets[i] = old.m_pTets[i];
      m_pBFaces[0] = old.m_pBFaces[0];
      m_pBFaces[1] = old.m_pBFaces[1];
      for (int i = 0; i < 6; i++)
	m_pVerts[i] = old.m_pVerts[i];
      m_fc = old.m_fc;
      m_nTets = old.m_nTets;

      return *this;
    }

    /// Determine what kind of face this is.
    ///
    /// On return, m_pTets contains the left and right tets for the input
    /// face.  Any additional entries are also tets that can be
    /// constructed from the same five verts as the first two.
    ///
    /// Verts are ordered like this:
    ///
    /// 2->3, 3->2 swaps: verts A, B, C are always the ones on the internal
    /// triangle in the two tet state.
    ///
    /// 2->2, 4->4 swaps: verts A, B, C are always the ones on the (a, for
    /// 4->4) current internal triangle.  A and B are two of the four
    /// coplanar verts (D and E are the others).
    ///
    /// Vert D is to the left of tri 012 (orientation(ABCD) is inverted).
    /// Vert E is to the right (orientation(ABCE) is correct).
    ///
    /// Vert F (only for 4->4 swaps) is the analog to vert C, but on the
    /// other side of plane ABDE.
    ///
    /// @param[in]  pF       Face to be categorized.
    ///
    /// @return  The category the face falls into.
    faceCat
    categorizeFace(Face* pF);
    faceCat
    getFaceCat() const
    {
      return m_fc;
    }

    Face*
    getFace() const
    {
      return m_pF;
    }
    Vert*
    getVert(GR_index_t i) const
    {
      assert(i < 6);
      return m_pVerts[i];
    }

    // All the tets as described above.
    void
    getTets(int& nTets, Cell* pTets[4]) const;
    // These functions return the vertex pointers as described above.
    Vert*
    getVertA() const
    {
      return (m_fc <= eT44) ? m_pVerts[0] : NULL;
    }
    Vert*
    getVertB() const
    {
      return (m_fc <= eT44) ? m_pVerts[1] : NULL;
    }
    Vert*
    getVertC() const
    {
      return (m_fc <= eT44) ? m_pVerts[2] : NULL;
    }
    Vert*
    getVertD() const
    {
      return (m_fc <= eT44) ? m_pVerts[3] : NULL;
    }
    Vert*
    getVertE() const
    {
      return (m_fc <= eT44) ? m_pVerts[4] : NULL;
    }
    Vert*
    getVertF() const
    {
      return (m_fc == eT44) ? m_pVerts[5] : NULL;
    }

    TriBFace*
    getBFaceA() const
    {
      return m_pBFaces[0];
    }
    TriBFace*
    getBFaceB() const
    {
      return m_pBFaces[1];
    }

    // For N32 and N44 cases, North and South are the two poles, and
    // Other is on the "equator".
    Vert*
    getVertNorth() const
    {
      return (m_fc == eN32 || m_fc == eN44) ? m_pVerts[0] : NULL;
    }

    Vert*
    getVertSouth() const
    {
      return (m_fc == eN32 || m_fc == eN44) ? m_pVerts[1] : NULL;
    }

    Vert*
    getVertOther() const
    {
      return (m_fc == eN32 || m_fc == eN44) ? m_pVerts[2] : NULL;
    }
  };

} // namespace GRUMMP
#endif
