#ifndef GR_Quality
#define GR_Quality 1

#include <stdio.h>
#include "GR_config.h"
#include "GR_QualMeasure.h"

void
calcCircumcenter(const double *adLoc0, const double *adLoc1,
		 const double *adLoc2, double adRes[]);

///
void
calcShortestToRadius(const CellSkel* const pC, int * const piN, double adRes[]);
void
calcShortestToRadius(const double *adLoc0, const double *adLoc1,
		     const double *adLoc2, int * const piN, double adRes[]);
void
calcShortestToRadius(const double *adLoc0, const double *adLoc1,
		     const double *adLoc2, const double *adLoc3,
		     int * const piN, double adRes[]);

///
class Quality {
public:
  /// Should we tabulate all values, min, or max?
  enum eEvalType {
    eAll, eMin, eMax
  };
private:
  ///
  const Mesh *pM;
  ///
  GRUMMP::QualMeasure* pQM;
  ///
  void
  (*pvfQualityFunc)(const CellSkel* const pC, int * const piN,
		    double adResult[]);
  ///
  int iQualEvals, iMaxEvals, iNBins;
  ///
  double dMaxQual, dMinQual, dScale;
  ///
  int *aiNumQualVals;
  ///
  double *adSignifImp; // Statistical measure of the significance of the
		       // improvement in the mesh quality measure since
		       /// the last assessment.
  double **a2dBinFrac;
  ///
  double *adMinQual, *adMaxQual;
  ///
  eEvalType eET;
  ///
  int
  iComputeBin(const double dQual) const;
  ///
  Quality&
  operator=(const Quality &);
  ///
  Quality(const Quality &);
public:
  ///
  Quality(Mesh* const pMIn, const int iMeas);
  ///
  ~Quality();
  ///
  void
  vEvaluate();
  /// Evaluate for only one cell.
  double
  dEvaluate(const CellSkel * const pC) const;
  ///
  void
  vWriteToFile(const char strQualFileName[]) const;
  ///
  double
  dWorstQuality() const
  {
    return adMinQual[iQualEvals - 1];
  }
  ///
  bool
  qMaximizeMinValue() const
  {
    return eET == eMin;
  }
  ///
  eEvalType
  getEvalType() const
  {
    return eET;
  }
  ///
  void
  setEvalType(const eEvalType eET_in)
  {
    eET = eET_in;
  }
  ///
  void
  vSetMesh(const Mesh * const pMIn)
  {
    pM = pMIn;
  }
  ///
  void
  vImportQualityData(const double adData[]);
};

#endif
