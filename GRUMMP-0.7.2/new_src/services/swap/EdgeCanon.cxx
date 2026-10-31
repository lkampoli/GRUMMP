#include "GR_assert.h"
#include "GR_EdgeCanon.h"

namespace GRUMMP
{

  void
  EdgeConfig::setFaceVert(const int face, const int v0, const int v1,
			  const int v2)
  {
    VALIDATE_INPUT(face >= 0 && face < m_numFaces);
    VALIDATE_INPUT(v0 >= 0 && v0 < m_numFaces + 2);
    VALIDATE_INPUT(v1 >= 0 && v1 < m_numFaces + 2);
    VALIDATE_INPUT(v2 >= 0 && v2 < m_numFaces + 2);
    m_faceVert[face][0] = v0;
    m_faceVert[face][1] = v1;
    m_faceVert[face][2] = v2;
  }

  const EdgeConfigs&
  EdgeConfigs::getInstance()
  {
    static EdgeConfigs instance;
    return instance;
  }

  EdgeConfigs::EdgeConfigs()
  {
    // Only need to create these for 3 through 7 verts; the offset is
    // taken care of in the retrieval call.
    aCS = new ConfSet[5];

    // If the original config has two or fewer verts, there is no
    // canonical config at all.

    //@ One config for three-to-two swapping.
    {
      aCS[0].setNumCanon(1);
      EdgeConfig& Conf = aCS[0].getConfig(0);
      Conf.setup(1, 1);
      Conf.setFaceVert(0, 0, 1, 2);
    }

    //@ One config for four-to-four swapping.
    {
      aCS[1].setNumCanon(1);
      EdgeConfig& Conf = aCS[1].getConfig(0);
      Conf.setup(2, 2);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
    }

    //@ One config for five-to-six swapping
    {
      aCS[2].setNumCanon(1);
      EdgeConfig& Conf = aCS[2].getConfig(0);
      Conf.setup(5, 3);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 4);
    }

    //@ Four configs for six-to-ten swapping
    aCS[3].setNumCanon(4);

    //@@ Fan configuration
    {
      EdgeConfig& Conf = aCS[3].getConfig(0);
      Conf.setup(6, 4);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 4);
      Conf.setFaceVert(3, 0, 4, 5);
    }

    //@@ Zig-zag 1 configuration
    {
      EdgeConfig& Conf = aCS[3].getConfig(1);
      Conf.setup(3, 4);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 5);
      Conf.setFaceVert(3, 3, 4, 5);
    }

    //@@ Zig-zag 2 configuration
    {
      EdgeConfig& Conf = aCS[3].getConfig(2);
      Conf.setup(3, 4);
      Conf.setFaceVert(0, 1, 2, 3);
      Conf.setFaceVert(1, 0, 1, 3);
      Conf.setFaceVert(2, 0, 3, 4);
      Conf.setFaceVert(3, 0, 4, 5);
    }

    //@@ Triangle configuration
    {
      EdgeConfig& Conf = aCS[3].getConfig(3);
      Conf.setup(2, 4);
      Conf.setFaceVert(0, 0, 2, 4);
      Conf.setFaceVert(1, 0, 1, 2);
      Conf.setFaceVert(2, 2, 3, 4);
      Conf.setFaceVert(3, 0, 4, 5);
    }

    //@ Six configs for seven-to-twelve swapping
    aCS[4].setNumCanon(6);
    //@@ Fan configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(0);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 4);
      Conf.setFaceVert(3, 0, 4, 5);
      Conf.setFaceVert(4, 0, 5, 6);
    }

    //@@ Zig-zag configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(1);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 6);
      Conf.setFaceVert(2, 2, 3, 6);
      Conf.setFaceVert(3, 3, 5, 6);
      Conf.setFaceVert(4, 3, 4, 5);
    }

    //@@ Zig-fan configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(2);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 1, 2, 3);
      Conf.setFaceVert(1, 0, 1, 3);
      Conf.setFaceVert(2, 0, 3, 4);
      Conf.setFaceVert(3, 0, 4, 5);
      Conf.setFaceVert(4, 0, 5, 6);
    }

    //@@ Fan-zig configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(3);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 4);
      Conf.setFaceVert(3, 0, 4, 6);
      Conf.setFaceVert(4, 4, 5, 6);
    }

    //@@ Fan-tri configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(4);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 3);
      Conf.setFaceVert(2, 0, 3, 5);
      Conf.setFaceVert(3, 0, 5, 6);
      Conf.setFaceVert(4, 3, 4, 5);
    }

    //@@ Fan-tri configuration
    {
      EdgeConfig& Conf = aCS[4].getConfig(5);
      Conf.setup(7, 5);
      Conf.setFaceVert(0, 0, 1, 2);
      Conf.setFaceVert(1, 0, 2, 4);
      Conf.setFaceVert(2, 0, 4, 5);
      Conf.setFaceVert(3, 0, 5, 6);
      Conf.setFaceVert(4, 2, 3, 4);
    }

    // Check to be sure that vertex index increases around each face
#ifndef NDEBUG
    for (int iOrig = 4; iOrig <= 7; iOrig++) {
      const ConfSet& CS = getConfSet(iOrig);
      for (int iConf = 0; iConf < CS.getNumCanon(); iConf++) {
	for (int iTri = 0; iTri < iOrig - 2; iTri++) {
	  assert(
	      (CS.getConfig(iConf).getFaceVert(iTri, 0)
		  < CS.getConfig(iConf).getFaceVert(iTri, 1))
		  && (CS.getConfig(iConf).getFaceVert(iTri, 1)
		      < CS.getConfig(iConf).getFaceVert(iTri, 2)));
	}
      }
    }
#endif
  }
}
