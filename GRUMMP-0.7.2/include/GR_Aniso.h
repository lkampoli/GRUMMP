#ifndef GR_Aniso
#define GR_Aniso 47

#include <stdio.h>

#include "GR_Cell.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Mesh2D.h"
#include "GR_LengthAniso2D.h"
#include "GR_CellQueue.h"
#include "GR_Mesh.h"
#include "GR_InsertionManager.h"
#include "SMsmooth.h"

using GRUMMP::SwappingInserter2D;
using GRUMMP::SwappingInserter3D;
class AnisoRefinement {
public:

  enum eMetricType {
    e2DTestQuarterCircle, e2DTestCross
  };

  AnisoRefinement(Mesh2D* pM2D_, char strMetricFileName_[FILE_NAME_LEN],
		  double dAlpha_);
  AnisoRefinement(Mesh2D* pM2D_, double dAlpha_);

  AnisoRefinement(Mesh2D* pM2D_, LengthAniso2D* length_, int iOrder_,
		  double dAlpha_);
  ~AnisoRefinement()
  {
    if (pvSmoothData)
      SMfinalizeSmoothing(pvSmoothData);
    if (pSI2D)
      delete pSI2D;
  }
  //Objects
  Mesh2D* pM2D;
  LengthAniso2D* m_length;
  SwappingInserter2D* pSI2D;
  void *pvSmoothData;
  char strMetricFileName[FILE_NAME_LEN];
  int iOrder;
  double dAlpha;
public:
  int iNMetricValues; //will be either 3 or 5 when metric is read

  //Public Functions

public:
  void
  vAssignMetricsToVerts(FILE *pFMetricFile);
  void
  vAssignMetricsToVerts(void
  (vFunc)(const double adLoc[], double adRetVal[]));

  /*
   * Sets a metric based on an analytic metric function (mesh must be a unit square mesh):
   *
   * - Quarter circle anisotropic features
   * - Cross anisotropic features
   */
  void
  vAssignMetricsToVertsTest(enum eMetricType Type);
  void
  vHealBdry();

  int
  iRefine(CellQueue* pCellQueue, int iMaxVertsAdded);
  int
  iAniCoarsen(CellQueue* pCellQueue, int IMaxVertsRemoved);

#ifdef HAVE_MESQUITE
  bool
  bRemoveVertexwithQual(Face* const FaceVertRemove);
  bool
  bRemoveVertex(Vert *pVert);
#endif

  void
  vImproveMesh(CellQueue* pCellQueue, int iMaxVerts, double dRate,
	       bool qVerbose);

  int
  iAniSwap(Face* pF);
  bool
  qAniInCircle(const Vert * const pVertA, const Vert * const pVertB,
	       const Vert * const pVertC, const Vert * const pVertD);
  bool
  qDoesSwapImprove(const Face* const pF);

  void
  vSmoothMetric();
  void
  vSmoothMetric(Vert* pV);
  bool
  qSmoothMetric(Vert* pV, const Vert* const pV1);

  double
  dAvgCircRadius();
  double
  dAvgArea();

  void
  vMetricLengths(const Cell* const pCell, double adLengths[3]);
  double
  dMetricLength(const double adLoc0[2], const double adLoc1[2],
		const double adMetric[5]);
  double
  dMetricLength(const Face* const pF);

  void
  vAniMidPoint(Face* pFFace, double adMidPoint[2]);

  bool
  qAniEncroachFace(BFace* pBF, double adLoc[2]);

  Cell*
  pCAnisoCircumCenter(Cell* pC, double adCircCen[2]);

  int
  iCheckBoundary(double adLoc[2]);
  int
  iCheckBoundary(double adLoc[2], Cell* pCellwLoc);

  double
  dAnisoCircumRadius(const Vert* const pVVert1, const Vert* const pVVert2,
		     const Vert* const pVVert3);
  double
  dAnisoCircumRadius(const Cell* const pC);

  double
  dAnisoArea(const Vert* const pVVert1, const Vert* const pVVert2,
	     const Vert* const pVVert3);
//	double dLogAnisoArea(const Vert* const pVVert1,const Vert* const pVVert2,const Vert* const pVVert3);
  double
  dAnisoArea(const Cell* const pC);

  void
  vInterpolateMetric(Vert* pVVert);
  void
  vInterpolateMetric(Vert* pVVert, std::set<Vert*>* setpVVerts);

  int
  iFixFace(Face* pF);
  int
  iFixFace(Cell* pC);

  void
  vWriteDistribution(char strOutName[]);

  bool
  qIsCurvedBdryOkay(Face* pF, const double *adVertLocC);

  void
  vFixForCurve();
};

#ifdef HAVE_MESQUITE
/*
 * This class is intended to be an extension of AnisoRefinement class into 3D. At them moment is only has
 * few functions for assigning metrics and removing points
 */

class AnisoRefinement3D {

public:

//	Two analytics test cases, must be set on a unitary cube mesh
  enum eMetricType {
    e3DTestQuarterCircle, e3DTestCross
  };

  AnisoRefinement3D(VolMesh* pVMesh_);
  ~AnisoRefinement3D()
  {

  }

  VolMesh* pVMesh;
  int iOrder;
  int iNMetricValues;

//		Assigns a metric for every point in the mesh, only analytic quality metrics at the moment
  void
  vAssignMetricsToVertsTest(enum eMetricType Type);

//		Decides which point to delete in an edge based on an average quality of its neighbouring cells,
//		the pioint is deleted by edge collapse
  bool
  bRemoveVertexwithQual(Vert *pV0, Vert *pV1);

//		Removes a point by edge collapse
  bool
  bRemoveVertex(Vert *pVert);
};
#endif
typedef struct virtualvert {
  Vert* pVRealVert;
  double dvX, dvY;  //virtual location of the vertex
} VirtualVert;

class VirtualTransMesh2D {
  std::vector<Cell *> vecpCCells;
  std::vector<VirtualVert> vecvvVerts;
  AnisoRefinement* pAnisoRef;
  Face* pEndOfMeshBdryFace;

  VirtualTransMesh2D() :
      vecpCCells(), vecvvVerts(), pAnisoRef(), pEndOfMeshBdryFace()
  {
    assert(0);
  }
  VirtualTransMesh2D(const VirtualTransMesh2D&);
  VirtualTransMesh2D&
  operator=(const VirtualTransMesh2D&);
public:
  VirtualTransMesh2D(AnisoRefinement* pAnisoRef, Cell* pCStart);

  VirtualVert
  vvFindVirtualVert(Vert* pVReal);

  friend class AnisoRefinement;

  void
  vAddCell(Face* pFforAdd);

  double adStart[2];
  double adDest[2];

  bool
  qDoesLineIntersectFace(double adP1[2], double adP2[2], Face* pF);

  void
  vBuildMeshToDest();
  void
  vAdaptToNonEmptyCC();
  void
  vUnTransformLoc(double adVirtLoc[2], double adRealLoc[2]);
  void
  vUnTransformDest(double adRealLoc[2])
  {
    vUnTransformLoc(adDest, adRealLoc);
  }
};

void
vTriCircumcenter(double adA[2], double adB[2], double adC[2], double adRes[2]);

double
dCircumRadius(double dL1, double dL2, double dL3);
double
dArea(double dL1, double dL2, double dL3);

double
dAniBiggness(const Cell* const pCell);
double
dAniPoorness(const Cell* const pCell);

bool
qCircumCirclesEmpty(double dLenDiag, double dLenL1, double dLenL2,
		    double dLenR1, double dLenR2); //Would a quad with these lengths be delaunay triangulated?

double
dSpaceDensity(const Vert* const pV, const double adNorm[2]);

#endif
