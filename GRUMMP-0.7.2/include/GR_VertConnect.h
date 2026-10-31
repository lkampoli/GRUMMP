#ifndef GR_VertConnect
#define GR_VertConnect 1

#include "stdlib.h"
#include "GR_config.h"
#include "GR_Vertex.h"
#include "GR_Mesh.h"

// This class is designed to store information about the neighbors of a
// vertex. These neighbors may be defined topologically or by means of a
// conflict graph (in which verts are "neighbors" if they are too close
// together according to current length scale).

class VertConnect {
private:
  size_t iNumNeigh;
  Vert **apVNeigh;
  double dPri;
  VertConnect(const VertConnect&) :
      iNumNeigh(0), apVNeigh(NULL), dPri(0)
  {
    assert(0);
  }
  VertConnect&
  operator=(const VertConnect&)
  {
    assert(0);
    return (*this);
  }

public:
  VertConnect() :
      iNumNeigh(0), apVNeigh(NULL), dPri(0)
  {
  }
  ~VertConnect()
  {
    if (apVNeigh) {
      delete[] apVNeigh;
      apVNeigh = NULL;
    }
  }

  void
  vSetPri(double dP)
  {
    dPri = dP;
  }

  double
  vGetPri()
  {
    return dPri;
  }

  void
  vSetSize(const size_t iN)
  {
    iNumNeigh = iN;
    if (apVNeigh)
      delete[] apVNeigh;
    apVNeigh = new Vert*[iN];
  }
  size_t
  iSize() const
  {
    return iNumNeigh;
  }
  void
  vSetNeigh(const int i, Vert* const pV)
  {
    assert(i >= 0 && unsigned(i) < iNumNeigh);
    apVNeigh[i] = pV;
  }
  Vert *
  pVNeigh(const int i) const
  {
    assert(i >= 0 && unsigned(i) < iNumNeigh);
    return apVNeigh[i];
  }
  Vert **
  ppVVerts() const
  {
    return apVNeigh;
  }
  int
  iCountKeptNeigh() const;
  double
  dClosestLegitSurfPoint(const Vert *pV,
			 const Vert::VertType VTLegitVertType) const;
  bool
  qNormalNeighborExists(const Vert * const pVCurr, const Vert * const pVPrev,
			Vert * &pVNext) const;
};

int
iCountKept2nd(const int iV, const VertConnect aVC[], const Mesh* const pMesh);

#endif
