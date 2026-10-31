#ifndef GR_EdgeCanon_h
#define GR_EdgeCanon_h

#include "GR_config.h"
#include "GR_Util.h"

namespace GRUMMP
{
  class EdgeConfig {
    int m_numPerm, m_numFaces;
    int (*m_faceVert)[3];
    EdgeConfig(const EdgeConfig&); // Deliberately not defined
    EdgeConfig&
    operator=(const EdgeConfig&); // Deliberately not defined
  public:
    EdgeConfig() :
	m_numPerm(0), m_numFaces(0), m_faceVert(NULL)
    {
    }
    ~EdgeConfig()
    {
      if (m_faceVert)
	delete[] m_faceVert;
    }
    void
    setup(const int numPerm, const int numFaces)
    {
      m_numPerm = numPerm;
      m_numFaces = numFaces;
      m_faceVert = new int[numFaces][3];
    }
    int
    getNumPerm() const
    {
      return m_numPerm;
    }
    int
    getFaceVert(const int face, const int vert) const
    {
      VALIDATE_INPUT(face >= 0 && face < m_numFaces);
      VALIDATE_INPUT(vert >= 0 && vert < 3);
      return m_faceVert[face][vert];
    }
    void
    setFaceVert(const int face, const int v0, const int v1, const int v2);
  };

  class ConfSet {
    int m_numCanon;
    EdgeConfig *m_configs;
    ConfSet(const ConfSet&); // Deliberately not defined.
    ConfSet&
    operator=(const ConfSet&); // Deliberately not defined
  public:
    ConfSet() :
	m_numCanon(0), m_configs(NULL)
    {
    }
    ~ConfSet()
    {
      if (m_configs)
	delete[] m_configs;
    }
    int
    getNumCanon() const
    {
      return m_numCanon;
    }
    void
    setNumCanon(const int i)
    {
      VALIDATE_INPUT(m_numCanon == 0);
      VALIDATE_INPUT(i > 0);
      m_numCanon = i;
      m_configs = new EdgeConfig[i];
    }
    EdgeConfig&
    getConfig(const int i) const
    {
      VALIDATE_INPUT(i >= 0 && i < m_numCanon);
      return m_configs[i];
    }
  };

// A file-local set of data describing the possible outcomes of
// reconfiguration.  The set is indexed by the number of equatorial
// verts in the initial configuration.
  class EdgeConfigs {
    ConfSet *aCS;
  protected:
    EdgeConfigs();
    ~EdgeConfigs()
    {
      delete[] aCS;
    }
    EdgeConfigs(const EdgeConfigs&);              // intentionally undefined 
    EdgeConfigs&
    operator=(const EdgeConfigs&);   // intentionally undefined
  public:
    static const EdgeConfigs&
    getInstance();
    const ConfSet&
    getConfSet(const int i) const
    {
      VALIDATE_INPUT(i >= 3 && i < 8);
      return aCS[i - 3];
    }
  };
}
#endif
