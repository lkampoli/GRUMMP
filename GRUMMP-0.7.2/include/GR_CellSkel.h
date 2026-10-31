#ifndef GR_CellSkel
#define GR_CellSkel 1

#include "GR_config.h"
#include "GR_Entity.h"
#include "GR_misc.h"
#include "GR_Vertex.h"

class CellSkel : public Entity {
public:
  ///
  enum eCellType {
    eTriCell = 0,
    eQuadCell,
    eTet,
    ePyr,
    ePrism,
    eHex,
    eFlake,
    eBdryEdge,
    eIntBdryEdge,
    eTriBFace,
    eIntTriBFace,
    eQuadBFace,
    eIntQuadBFace,
    eTriCV,
    eQuadCV,
    eTetCV,
    ePyrCV,
    ePrismCV,
    eHexCV,
    eInvalid
  };
protected:
  // This data doesn't all get used by all CellSkel descendents, but it
  // makes sense to define it here, because that way the bitfields
  // really get packed tightly by the compiler.
  unsigned int m_region :iRegionBits, m_cellType :6;
  bool m_deleted :1;
  ///
  CellSkel&
  operator=(const CellSkel& CS)
  {
    if (&CS != this) {
      copyAllFlags(CS);
    }
    return *this;
  }
public:
  ///
  CellSkel() :
      m_region(0), m_cellType(eInvalid), m_deleted(false)
  {
  }
  ///
  virtual
  ~CellSkel()
  {
  }
  ///
  bool
  isValid() const
  {
    return (this != NULL);
  }
  ///
  // Asks if it is a bdry cell.. return false by default
  virtual bool
  isBdryCell() const
  {
    return false;
  }
  ///
  virtual int
  doFullCheck() const;
  ///
  void
  setDefaultFlags()
  {
    m_region = 0;
    // Type is not reset; this is set once and for all by the
    // derived-class constructor.
    m_deleted = false;
  }
  ///
  void
  copyAllFlags(const CellSkel &C)
  {
    m_region = C.m_region;
    m_cellType = C.m_cellType;
    m_deleted = C.m_deleted;
  }
  ///
  void
  setRegion(int iReg)
  {
    assert(iReg >= 0 && iReg < iMaxRegionLabel);
    m_region = (iReg & 0x7F);
  }
  ///
  int
  getRegion() const
  {/*assert(isValid());*/
    return m_region;
  }
  ///
  int
  getType() const
  {
    assert(isValid());
    return m_cellType;
  }
  ///
  void
  setType(eCellType eCT)
  {
    m_cellType = (eCT & 0x1F);
  }
  ///
  void
  markAsDeleted()
  {
    m_deleted = true;
  }
  ///
  bool
  isDeleted() const
  {
    assert(isValid());
    return m_deleted;
  }
  ///
  virtual int
  getNumFaces() const = 0;
  ///
  virtual int
  getNumVerts() const = 0;
  ///
  virtual const Vert*
  getVert(const int i) const = 0;
  ///
  virtual Vert*
  getVert(const int i) = 0;
  ///
  virtual double
  calcSize() const
  {
    assert(0);
    return -1.e300;
  }
  ///
  double
  calcMaxDihed() const;
  ///
  double
  calcMinDihed() const;
  ///
  double
  calcMaxSolid() const;
  ///
  double
  calcMinSolid() const;
  ///
  virtual void
  calcAllDihed(double adDihed[], int* const piNDihed, const bool in_degrees =
		   true) const = 0;
  ///
  virtual void
  calcAllSolid(double adSolid[], int* const piNSolid, const bool in_degrees =
		   true) const = 0;
  ///
  virtual double
  calcShortestEdgeLength() const
  {
    assert(0);
    return -1.e300;
  }

  void
  printCellInfo() const;
};

#endif
