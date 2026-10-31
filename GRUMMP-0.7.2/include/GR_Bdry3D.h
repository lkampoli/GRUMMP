#ifndef GR_Bdry3D
#define GR_Bdry3D 1

#include <vector>

#include "GR_config.h"

#include "GR_misc.h"
#include "GR_Vertex.h"
#include "GR_Cell.h"
#include "GR_misc.h"
#include "GR_EntContainer.h"
#include "GR_Bdry.h"

// This is the base class for all 3D boundary data (including internal
// boundaries between regions in a multi-region domain).  Eventually
// this will have a flock of derived classes hanging from it.
class BdryPatch3D : public BdryPatch {
  /// Copy construction disallowed
  BdryPatch3D(const BdryPatch3D&); // Never defined
  /// operator= disallowed
  BdryPatch3D&
  operator=(const BdryPatch3D&); // Never defined
protected:
  int iNVerts, iVSize, iNVertsToReinsert, iVRSize;
  Vert **apVVerts;
  Vert **apVVertsToReinsert;
  BdryPatch3D* pB3DChild;
  BdryPatch3D* pB3DNextSib;
public:
/// Constructor and destructor  
  BdryPatch3D() :
      BdryPatch(), iNVerts(0), iVSize(0), iNVertsToReinsert(0), iVRSize(0), apVVerts(
      NULL), apVVertsToReinsert(NULL), pB3DChild(
      NULL), pB3DNextSib(NULL)
  {
  }
  BdryPatch3D(const int iBCL, const int iBCR, const int iRL, const int iRR,
	      const int iNV, Vert* apVV[]) :
      BdryPatch(iBCL, iBCR, iRL, iRR), iNVerts(iNV), iVSize(iNV), iNVertsToReinsert(
	  0), iVRSize(0), apVVerts(NULL), apVVertsToReinsert(NULL), pB3DChild(
      NULL), pB3DNextSib(NULL)
  {
    // There had better be a finite number of vertices, which had better
    // be provided as input.
    assert(iNV > 0 && apVV);
    apVVerts = new Vert*[iNV];
    for (int i = 0; i < iNV; i++)
      apVVerts[i] = apVV[i];
  }
  virtual
  ~BdryPatch3D()
  {
    // Disconnect children and sibs
    if (pB3DChild) {
      pB3DChild->~BdryPatch3D();
      pB3DChild = NULL;
    }
    if (pB3DNextSib) {
      pB3DNextSib->~BdryPatch3D();
      pB3DNextSib = NULL;
    }
    // Get rid of vertex data
    if (apVVerts) {
      delete[] apVVerts;
      apVVerts = NULL;
    }
    if (apVVertsToReinsert) {
      delete[] apVVertsToReinsert;
      apVVertsToReinsert = NULL;
    }
  }
  virtual void
  vComputeMapping(const BdryPatch* const pBP, int& iDim, double& dScale,
		  double adTrans[3], double a2dRot[3][3]) const = 0;
  enum eGeom3D {
    ePolygon, eUnknown
  };
  virtual enum eGeom3D
  eType() const = 0;
  void
  vResize(const int iNewSize)
  {
    assert(iNewSize > iNVerts);
    if (iNewSize > iVSize) {
      Vert **apVNew = new Vert*[iNewSize];
      for (int i = 0; i < iNVerts; i++)
	apVNew[i] = apVVerts[i];
      delete[] apVVerts;
      apVVerts = apVNew;
      iVSize = iNewSize;
    }
  }
  void
  vResizeInteriorVerts(const int iNewSize)
  {
    assert(iNewSize > iNVertsToReinsert);
    if (iNewSize > iVRSize) {
      Vert **apVNew = new Vert*[iNewSize];
      for (int i = 0; i < iNVertsToReinsert; i++)
	apVNew[i] = apVVertsToReinsert[i];
      delete[] apVVertsToReinsert;
      apVVertsToReinsert = apVNew;
      iVRSize = iNewSize;
    }
  }
  Vert *
  pVVert(const int i) const;
  Vert *
  pVInteriorVert(const int i) const
  {
    assert(i >= 0 && i < iNVertsToReinsert);
    return apVVertsToReinsert[i];
  }
  void
  vAddInteriorVert(Vert* const pV);
  void
  vAddChild(BdryPatch3D* const pB3D);
  bool
  qHasChildren() const
  {
    return (pB3DChild != NULL);
  }
  int
  iNumChildren() const;
  int
  iNumSegments() const;
  int
  iNumInteriorVerts() const
  {
    return iNVertsToReinsert;
  }
  BdryPatch3D*
  pBP3DChild() const
  {
    return pB3DChild;
  }
  BdryPatch3D*
  pBP3DSibling() const
  {
    return pB3DNextSib;
  }
  void
  vGetSegmentVerts(const int iSeg, Vert *& pV0, Vert *& pV1) const;
  bool
  qHasVert(const Vert* const pV) const;
  void
  vGetVertNeighbors(const Vert* const pV, Vert* &pVN0, Vert* &pVN1) const;
  void
  vSplitSegmentIfFound(const Vert* const pV0, const Vert* const pV1, Vert *pV);
  // Stuff related to bdry patch optimization
//  void vGetVertsToReinsert(::std::vector<Vert*> vecVert);
  virtual int
  iNumTriangles() const = 0;
  // Geometric stuff
  virtual double
  dDistFromBdry(const double adLoc[], int * const iSide = NULL) const = 0;
  virtual bool
  qDisjointFrom(const Vert* const pV) const;
  virtual bool
  qDisjointFrom(const BdryPatch* const pBP) const;
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const = 0;
  virtual void
  vFlipNormal() = 0;

};

class BdryPolygon : public BdryPatch3D {
private:
  // Disable the default constructor
  BdryPolygon();
  double adNormal[3];
public:
  BdryPolygon(const int iBCL, const int iBCR, const int iRL, const int iRR,
	      const int iNV, Vert *apVV[]);
  virtual
  ~BdryPolygon()
  {
  }
  virtual enum eGeom3D
  eType() const
  {
    return ePolygon;
  }

  virtual void
  vComputeMapping(const BdryPatch* const pBP, int& iDim, double& dScale,
		  double adTrans[3], double a2dRot[3][3]) const;

  virtual double
  dDistFromBdry(const double adLoc[], int * const piSide = NULL) const;
  virtual void
  vUnitNormal(const double adLoc[], double adNorm[]) const;
  virtual void
  vFlipNormal();
  void
  B2DProjection(const Bdry3D& B3D, Bdry2D& B2D) const;
  // virtual void vCurvature(const double adLoc[], double a2dCurve[2][2]);
  int
  iLocalVertIndex(const double adLoc[3]) const;
  void
  vFindProjectionBasis(double adOffset[3], double adBasis1[3],
		       double adBasis2[3]) const;
  // Stuff related to bdry patch optimization; returns true if merge
  // succeeded.
  bool
  qMerge(const BdryPolygon& BPoly, const Vert * const pV,
	 int &iNumNewInteriorVerts);
  virtual int
  iNumTriangles() const;
};

class Bdry3D : public BdryRep {
protected:
  BdryPatch **apBPatches;
private:
  EntContainer<Vert> ECVerts;
  Bdry3D(const Bdry3D&) :
      BdryRep(), apBPatches(NULL), ECVerts()
  {
    assert(0);
  }
  Bdry3D&
  operator=(const Bdry3D&)
  {
    assert(0);
    return (*this);
  }
public:
  Bdry3D() :
      BdryRep(), apBPatches(NULL), ECVerts()
  {
  }
  Bdry3D(const EntContainer<Vert>& ECV);
  Bdry3D(const char strBaseFileName[], const int iOptLevel = 1);
  ~Bdry3D()
  {
    if (apBPatches != NULL) {
      for (int i = 0; i < iNPatches; i++) {
	delete apBPatches[i];
      }
      delete[] apBPatches;
    }
  }
  ///
  EntContainer<Vert>*
  pECGetVerts()
  {
    return (&ECVerts);
  }
  GR_index_t
  iNumBdryPatches() const
  {
    return iNPatches;
  }
  GR_index_t
  iNumVerts() const
  {
    return ECVerts.lastEntry();
  }
  GR_index_t
  iNumTriangles() const;
  BdryPatch3D*
  pBdry(const int i) const
  {
    assert(i >= 0 && i < iNPatches);
    return (dynamic_cast<BdryPatch3D*>(apBPatches[i]));
  }
  void
  vAddPatch(BdryPatch3D* const pBCh);
  BdryPatch3D*
  operator[](const int i) const
  {
    return pBdry(i);
  }
  Vert*
  pVVert(const GR_index_t i) const
  {
    assert(i < ECVerts.lastEntry());
    return ECVerts.getEntry(i);
  }
  Vert *
  pVNewVert()
  {
    Vert *pV = ECVerts.getNewEntry();
    pV->setDefaultFlags();
    return (pV);
  }
  // Add a single vertex consistently to the boundary
  Vert*
  pVAddVertex(const Vert* const pV0, const Vert* const pV1,
	      const double adCent[]);
  // Add vertices to the boundary so that no vertex encroaches on any
  // curve in the boundary.  For polygonal boundaries, this means simply
  // that no vertex can be strictly inside the equatorial sphere of any
  // boundary segment or subsegment.
  void
  vFixEncroachment();
  // This should return integer indices rather than pointers because
  // the vertex pointers will be different in the 3D mesh.  The indices
  // are what's needed to be sure to recover edges.
  GR_index_t
  iVertIndex(const int iPoly, const int iVert) const
  {
    assert(iPoly < iNPatches && iPoly >= 0);
    return ECVerts.getIndex(pBdry(iPoly)->pVVert(iVert));
  }
  GR_index_t
  iVertIndex(const Vert* const pV) const
  {
    return ECVerts.getIndex(pV);
  }
  bool
  readBdryFile(const char strBaseFileName[]);
  bool
  readSMeshFile(const char strBaseFileName[]);
  bool
  readHINFile(const char strBaseFileName[]);
  void
  vOptimizePatches(const int iLevel = 1);
  int
  iNumPatchesContaining(const Vert* const pVS) const;
};

#endif
