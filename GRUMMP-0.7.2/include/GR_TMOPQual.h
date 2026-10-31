/*
 * GR_TMOPQual.h
 *
 *  Created on: Jun 30, 2017
 *      Author: cfog
 */

#ifndef INCLUDE_GR_TMOPQUAL_H_
#define INCLUDE_GR_TMOPQUAL_H_

#include <set>

#include "GR_config.h"
#include "GR_Cell.h"

#ifdef HAVE_MESQUITE
#include "TQualityMetric.hpp"
#endif

namespace GRUMMP
{

#ifdef HAVE_MESQUITE
  /*
   * This class is set to measure quality with TMOP for the ideal target, a unit tetrahedron.
   * Quality is based on TMOP with a shape quality metric. A metric doesn't need to be
   * previously defined, all cells use the ideal target.
   */

  /*
   * This class evaluates quality cells based on TMOP.
   * By setting a target defined by an eigen decomposition of a metric, quality can be computed
   * for several different quality criteria. A metric must be previously set
   */
  class TMOPQual {

  public:
    TMOPQual(int dimension = 2) :
	iVs(0), iNM(0), iT(0), M_Qmetric(), bMaxQualperCell(true)
    {
      if (dimension == 3)
	iNMetricValues = 6;
      else
	iNMetricValues = 3;
    }
    ~TMOPQual()
    {
    }

    enum eQMetric {
      eDefault = 0, eShape, eSize, eShapeSize, eShapeOrient, eShapeSizeOrient
    };
    enum eType {
      eAverage, eLinear, eLogEuclidean
    };

// Computes a target for a given set of points
    void
    SetTargetfromMetric(double *adTarget, CellSkel::eCellType Type,
			const Vert* pV0, const Vert* pV1, const Vert* pV2,
			const Vert* pV3 = NULL, const Vert* pV4 = NULL,
			const Vert* pV5 = NULL, const Vert* pV6 = NULL,
			const Vert* pV7 = NULL);

// Computes a target for a metric and a given set of points.
// Used to compute targets for a quality of a possible cell to create
    void
    SetTargetfromMetricNewVertice(double *adTarget, CellSkel::eCellType Type,
				  double *Metric, const Vert* pV0,
				  const Vert* pV1, const Vert* pV2 = NULL,
				  const Vert* pV3 = NULL,
				  const Vert* pV4 = NULL,
				  const Vert* pV5 = NULL,
				  const Vert* pV6 = NULL,
				  const Vert* pV7 = NULL);

//	 Evaluates quality with TMOP based on a set of points
    double
    dEval(enum eQMetric eQM, CellSkel::eCellType Type, const Vert* pV0,
	  const Vert* pV1, const Vert*pV2, const Vert* pV3 = NULL,
	  const Vert* pV4 = NULL, const Vert* pV5 = NULL,
	  const Vert* pV6 = NULL, const Vert* pV7 = NULL);

//	 Evaluates quality with TMOP based on a set of coords
    double
    dEval(enum eQMetric eQM, CellSkel::eCellType Type, double* adTarget,
	  const double *adLoc0, const double *adLoc1, const double *adLoc2,
	  const double *adLoc3 = NULL, const double *adLoc4 = NULL,
	  const double *adLoc5 = NULL, const double *adLoc6 = NULL,
	  const double *adLoc7 = NULL);

//	 Evaluates which configuratation is the best target
    double
    dEvalgetBestTarget(enum eQMetric eQM, CellSkel::eCellType Type, double* dT,
		       double *dTbest, const double *adLoc0,
		       const double *adLoc1, const double *adLoc2);

// Functions for a Log-Euclidean metric interpolation in 2D
    void
    vRLamdaRt2D(double *adVec, double *adEig, double *adM);
    void
    vLog2D(double *adM, double *adMLog);
    void
    vExp2D(double *adM, double *adMExp);

    // Functions for a Log-Euclidean metric interpolation in 3D
    void
    vRLamdaRt3D(double *adVec, double *adEig, double *adM);
    void
    vLog3D(double *adM, double *adMLog);
    void
    vExp3D(double *adM, double *adMExp);

// Computes metric length in 2D and 3D
    double
    dMLength2D(const Face* const pF);
    double
    dMLength2D(const double adLoc0[2], const double adLoc1[2],
	       const double adMet[3]);
    double
    MLength3D(const double adLoc0[3], const double adLoc1[3],
	      const double adMet[6]);

//	Computes log metric length
    double
    dLogMetricLength(const Face* const pF);
    double
    dLogMetricLength(const Vert *pV0, const Vert *pV1);
    double
    dLogMetricLength(double *Coords0, double *adMet0, double *Coords1,
		     double *adMet1);

//	Computes log metric interpolation for a cell
    void
    dLogMetricCell(const Cell *pC, double *adM);
    void
    dLogMetricCell(const Vert *pV0, const Vert *pV1, const Vert *pV2,
		   double *adM);
    void
    dLogMetricCell(double *adMet0, double *adMet1, double *adMet2, double *adM);

//	Computes log metric interpolation for a new point
    void
    LogEuclideanInterpMetric(Vert* pV, std::set<Vert*> setpV);
    void
    LogEuclideanInterpMetric(double *dCoords, double *adNewM,
			     std::set<Vert*> setpV);

  protected:

//	Sets a quality metric
    Mesquite::TMetric*
    SingleQualMetric(enum eQMetric eQM = eDefault);

//	Returns true if metric is available
    bool
    isMetricAvailable(CellSkel::eCellType Type, const Vert* pV0,
		      const Vert* pV1, const Vert*pV2, const Vert* pV3 = NULL,
		      const Vert* pV4 = NULL, const Vert* pV5 = NULL,
		      const Vert* pV6 = NULL, const Vert* pV7 = NULL);

//	Checks if a metric is available
    bool
    CheckMetric(const Vert *apV[]);

//	 Computes a target from an eigen decomposition
    void
    SetTargetfromEigen(CellSkel::eCellType Type, double *dT, double *dVec,
		       double *dEV);

//	 Computes eigen decomposition in 2D and 3D
    void
    getEigenDecomposition(double *dA, double *dVec, double *dEV);
    void
    EigenDecomposition2D(double *dA, double *dVec, double *dEV);
    void
    EigenDecomposition3D(double *dA, double *dVec, double *dEV);

//	 Functions for eigen decomposition in 3D
    void
    tred(double V[3][3], double d[3], double e[3]);
    void
    tql(double V[3][3], double d[3], double e[3]);

    int iVs;
    int iNM;
    int iT;
    int iNMetricValues;

  public:

    void
    CheckTypeCell(CellSkel::eCellType Type);
    void
    SetMetricValue(int dMetricValue)
    {
      iNMetricValues = dMetricValue;
    }
    ;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    //	 For Merging Triangular cells into Quad cells

    /*
     * Merges a triangular mesh in to a quad mesh. Merge is based on quality in a greedy approach.
     *
     * bool Maxqual = if set, merges a quad if both tri are best quality for merging
     * pM2D = mesh to merge
     * dQRange = quality threshold for which merging is allowed
     */
    void
    MergeTriMeshToQuadMesh(bool MaxQual, Mesh2D *pM2D, double *dQRange = NULL);

    /*
     * Evaluates all cell and sets which cells should merge, by identifying which edges from the quadrilateral formed between two adjacent cells should be erased.
     *
     * pM2D = mesh to evaluate
     * dQRaneg = quality threshold for which merging is allowed.
     */
    void
    EvalMeshTritoQuad(Mesh2D *pM2D, double *dQRange);

    /*
     * Evaluates a cell for merge with its neighbouring cells
     *
     * pC = cell to evalaute for merge
     * dQRange = quality threshold for which merging is allowed
     */
    bool
    EvalCellToQuad(Cell *pC, double *dQRange);

    /*
     * Evaluates if a specific face is the best quality option to merge a cell
     *
     * pFToMerge = specific face to evaluate
     * pF1 = other face of cell
     * pF2 = other face of cell
     */
    bool
    EvalFaceToMerge(Face *pFToMerge, Face *pF1, Face *pF2);

    /*
     * For a triangular cell A, which best quality for merge is cell B.
     * Checks if for cell B best quality for merge is also cell A
     *
     * pF = Face between two cells to merge
     * pC = cell to merge
     */
    bool
    CheckOppositeCellToQuad(Face *pF, Cell* pC);

    /*
     * Evaluates quality with TMOP for a quad cell
     *
     * pCell = cell to evaluate quality
     * dNewCellsCoords = evaluates quality for a potential nonexistent cell by giving a set of coordinates
     * ForNewCell = if set true, must give the coords of the potential cell with dNewCellsCoords
     */
    double
    EvalQuadCell(QuadCell* pCell, double *NewCellCoords, bool ForNewCell);

    /*
     * Evaluates quality with TMOP for two triangular cells to merge into a quad
     *
     * pFace = face between to triangular cells to merge
     * dNewCellsCoords = evaluates quality for a potential nonexistent cell by giving a set of coordinates
     * ForNewCell = if set true, must give the coords of the potential cell with dNewCellsCoords
     */
    double
    EvalTriToQuad(Face *pFace, double *NewCellCoords, bool ForNewCell);

    /*
     * Evaluates quality with TMOP for four diferent points of a quad
     *
     * pV0 -pV3 = points that define a quad
     * dNewCellsCoords = evaluates quality for a potential nonexistent cell by giving a set of coordinates
     * ForNewCell = if set true, must give the coords of the potential cell with dNewCellsCoords
     */
    double
    EvalTriToQuad(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3,
		  double *NewCellCoords, bool ForNewCell);

    /*
     * Computes a Target based on an eigen decomposition of a metric for a quad cell, uses average metric from all points
     *
     * pV0 - pV3 = points that define a quad
     * Target = returns an array of size 12 with a quad target
     */
    void
    SetQuadTargetfromMetric(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3,
			    double *Target) const;

    /*
     * Checks for convexity for between four different points
     *
     * pV0 - pV3 = points that define a quad
     */
    bool
    CheckConvexQuad(Vert *pV0, Vert *pV1, Vert *pV2, Vert *pV3);

    /*
     * Merges 2 triangular cells that share a face into a quad
     *
     * pM2D = pointer to mesh
     * pF = face adjacent the cells to merge
     */
    QuadCell *
    CombineTrisToQuad(Mesh2D *pM2D, Face *pF);

    Mesquite::TMetric *M_Qmetric;
    bool bMaxQualperCell;
    //////////////////////////////////////////////////////////////////////////////////////////////////////
  };
#endif

} // end namespace GRUMMP

#endif /* INCLUDE_GR_TMOPQUAL_H_ */
