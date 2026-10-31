#ifndef GR_CellCV
#define GR_CellCV 1

#include "GR_config.h"

#include "GR_CellSkel.h"
#include "GR_Vertex.h"

static CellCV* const pCInvalidCellCV = static_cast<CellCV*>(NULL);

/// This class and its descendants are deprecated; don't use them in new
/// code. CFOG, 10/15/2012
class CellCV : public CellSkel {
private:
  // These don't make sense.
  const Face*
  getFace(const int) const
  {
    assert(0);
    return pFInvalidFace;
  }
  const Cell*
  getCell(const int) const
  {
    assert(0);
    return pCInvalidCell;
  }
public:
  /// Override copy constructor
  CellCV(const CellCV& CCV) :
      CellSkel(), iNV(CCV.iNV), ppVVerts(new Vert*[iNV])
  {
    for (int ii = 0; ii < iNV; ii++) {
      ppVVerts[ii] = CCV.ppVVerts[ii];
    }
  }
protected:
  CellCV&
  operator=(const CellCV& CCV)
  {
    if (this != &CCV) {
      logMessage(2, "In CellCV op=\n");
      delete[] ppVVerts;
      iNV = CCV.iNV;
      ppVVerts = new Vert*[iNV];
      for (int ii = 0; ii < iNV; ii++) {
	ppVVerts[ii] = CCV.ppVVerts[ii];
      }
      copyAllFlags(CCV);
    }
    return (*this);
  }
protected:
  ///
  int iNV;
  ///
  Vert** ppVVerts;
public:
  ///
  void
  resetAllData()
  {
    for (int i = 0; i < iNV; ppVVerts[i++] = pVInvalidVert)
      ;
  }
  ///
  inline
  CellCV(const int iNVerts, Vert** const ppVIn = 0, int iDim = 2);
  ///
  inline virtual
  ~CellCV();
  ///
  bool
  qValid() const
  {
    return (this != pCInvalidCellCV);
  }
  ///
  inline int
  iFullCheck() const;
  ///
  int
  getNumVerts() const
  {
    return iNV;
  }
  ///
  int
  getNumFaces() const
  {
    if (getVert(0)->getSpaceDimen() == 2) {
      return (getNumVerts());
    }
    else {
      switch (getNumVerts())
	{
	case 3:
	  return 3; // Tri in surface mesh
	case 4:
	  return 4; // Tet; or quad in surface mesh
	case 5:
	  return 5; // Pyr
	case 6:
	  return 5; // Prism
	case 8:
	  return 6; // Hex
	default:
	  assert(0);
	  return -1;
	}
    }
  }
  int
  getNumCells() const
  {
    return 1;
  }
  inline virtual const Vert*
  getVert(const int i) const;
  inline virtual Vert*
  getVert(const int i)
  {
    return (const_cast<Vert*>(static_cast<const CellCV*>(this)->getVert(i)));
  }
  void
  getAllVertHandles(GRUMMP_Entity* aHandles[]) const
  {
    for (int i = 0; i < iNV; i++)
      aHandles[i] = ppVVerts[i];
  }
  void
  getAllFaceHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  void
  getAllCellHandles(GRUMMP_Entity*[]) const
  {
    assert(0);
  }
  int
  getEntType() const
  {
    if (getVert(0)->getSpaceDimen() == 2 || getNumVerts() == 3)
      // Incorrect for QuadCellCV; must override there.  For QuadCV's
      // that have been declared as CellCV's, there doesn't seem to be
      // any immediate hope.
      return iBase_FACE;
    else
      return iBase_REGION;
  }
  int
  getEntTopology() const
  {
    switch (getNumVerts())
      {
      case 2:
	return iMesh_LINE_SEGMENT;
      case 3:
	return iMesh_TRIANGLE;
	// Incorrect for QuadCellCV; must override there.  For QuadCV's
	// that have been declared as CellCV's, there doesn't seem to be
	// any immediate hope.
      case 4:
	return iMesh_TETRAHEDRON;
      case 5:
	return iMesh_PYRAMID;
      case 6:
	return iMesh_PRISM;
      case 8:
	return iMesh_HEXAHEDRON;
      default:
	assert(0);
	return iMesh_ALL_TOPOLOGIES; // Bogus!
      }
  }
  ///
  inline bool
  hasVert(const Vert* pV) const;
  ///
  virtual double
  calcSize() const;
  ///
  virtual void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const;
  ///
  virtual void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const;
};

inline
CellCV::CellCV(const int iNVerts, Vert** const ppVIn, int iDim) :
    iNV(iNVerts), ppVVerts(NULL)
{
  assert(iNVerts >= 2 && iNVerts <= 8);
  ppVVerts = new Vert*[iNVerts];
  int i;
  if (ppVIn)
    for (i = 0; i < iNVerts; i++) {
      ppVVerts[i] = ppVIn[i];
      assert2(ppVVerts[i] != pVInvalidVert, "Uninitialized vert");
    }
  else
    for (i = 0; i < iNVerts; i++)
      ppVVerts[i] = pVInvalidVert;
  switch (iNVerts)
    {
    case 2:
      break;
    case 3:
      m_cellType = eTriCV;
      break;
    case 4:
      if (iDim == 2)
	m_cellType = eQuadCV;
      else
	m_cellType = eTetCV;
      break;
    case 5:
      m_cellType = ePyrCV;
      break;
    case 6:
      m_cellType = ePrismCV;
      break;
    case 8:
      m_cellType = eHexCV;
      break;
    default:
      break;
    }
}

inline
CellCV::~CellCV()
{
  if (ppVVerts) {
    delete[] ppVVerts;
    ppVVerts = NULL;
  }
}

inline const Vert*
CellCV::getVert(const int i) const
{
  assert2(i >= 0 && i < iNV, "Index out of range");
  return ppVVerts[i];
}

inline int
CellCV::iFullCheck() const
{
  int iRetVal = qValid();
  int iNVerts = getNumVerts();
  for (int i = 0; iRetVal && i < iNVerts; i++) {
    const Vert* pVThisVert = getVert(i);
    iRetVal = iRetVal && pVThisVert->isValid();
  }
  return (iRetVal);
}

inline bool
CellCV::hasVert(const Vert* pV) const
{
  assert(pV->isValid());
  for (int i = 0; i < iNV; i++)
    if (ppVVerts[i] == pV)
      return (true);
  return (false);
}

///
class TriCellCV : public CellCV {
public:
  ///
  TriCellCV(Vert** const ppVIn = 0) :
      CellCV(3, ppVIn, 2)
  {
  }
  TriCellCV(Vert * const pVA, Vert * const pVB, Vert * const pVC) :
      CellCV(3, NULL, 2)
  {
    ppVVerts[0] = pVA;
    ppVVerts[1] = pVB;
    ppVVerts[2] = pVC;
  }
  ///
  void
  vAssign(Vert* const pV0, Vert* const pV1, Vert* const pV2)
  {
    ppVVerts[0] = pV0;
    ppVVerts[1] = pV1;
    ppVVerts[2] = pV2;
  }
  void
  vSwitchOrientation()
  {
    Vert *pVTmp = ppVVerts[2];
    ppVVerts[2] = ppVVerts[1];
    ppVVerts[1] = pVTmp;
  }
};
//
/////
//class QuadCellCV : public CellCV {
//public:
//  ///
//  QuadCellCV(Vert** const ppVIn = 0) : CellCV(4, ppVIn, 2) {}
//  void vSwitchOrientation()
//    {
//      Vert *pVTmp = ppVVerts[3];
//      ppVVerts[3] = ppVVerts[0];
//      ppVVerts[0] = pVTmp;
//
//      pVTmp = ppVVerts[2];
//      ppVVerts[2] = ppVVerts[1];
//      ppVVerts[1] = pVTmp;
//    }
//  int getEntType() const {return iBase_FACE;}
//  int getEntTopology() const {return iMesh_QUADRILATERAL;}
//};

///
class TetCellCV : public CellCV {
public:
  ///
  TetCellCV(Vert** const ppVIn = 0) :
      CellCV(4, ppVIn, 3)
  {
  }
  TetCellCV(Vert * const pVA, Vert * const pVB, Vert * const pVC,
	    Vert * const pVD) :
      CellCV(4, NULL, 3)
  {
    ppVVerts[0] = pVA;
    ppVVerts[1] = pVB;
    ppVVerts[2] = pVC;
    ppVVerts[3] = pVD;
  }
  ~TetCellCV()
  {
  }
  void
  vCircumcenter(double adCircCent[3]) const;
  double
  calcSize() const;
  void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const;
  void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const;
};

///
//class PyrCellCV : public CellCV {
//public:
//  ///
//  PyrCellCV(Vert** const ppVIn = 0) : CellCV(5, ppVIn, 3) {}
//};
//
/////
//class PrismCellCV : public CellCV {
//public:
//  ///
//  PrismCellCV(Vert** const ppVIn = 0) : CellCV(6, ppVIn, 3) {}
//};
//
/////
//class HexCellCV : public CellCV {
//public:
//  ///
//  HexCellCV(Vert** const ppVIn = 0) : CellCV(8, ppVIn, 3) {}
//};

#endif
