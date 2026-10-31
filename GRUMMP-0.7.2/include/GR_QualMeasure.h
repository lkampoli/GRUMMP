#ifndef GR_QualMeasure_h
#define GR_QualMeasure_h

#include <math.h>
#include <float.h>
#include <string>
#include <set>

#include "GR_config.h"
#include "GR_Cell.h"
#include "GR_Util.h"

void
tred(double V[3][3], double d[3], double e[3]);
void
tql(double V[3][3], double d[3], double e[3]);
void
EigenDecomposition(double A[3][3], double V[3][3], double d[3]);

namespace GRUMMP
{
  class QualMeasure {
  protected:
    /// String description of the quality measure.
    std::string m_name;
    int m_dim;
    bool bEvalMetric;
  public:
    QualMeasure(std::string name, const int dim) :
	m_name(name), m_dim(dim), bEvalMetric(false)
    {
    }
    virtual
    ~QualMeasure()
    {
    }

    std::string
    getName() const
    {
      return m_name;
    }
  protected:
    virtual bool
    shouldMinimize() const = 0;

  public:
    /// \brief This is the value that optimizing code should compute, as
    /// maximizing this will maximize the min or minimize the max, as
    /// appropriate.

    /// Because this mode of use will often be focusing on scenarios where
    /// verts are known rather than cells, the function signature
    /// accomodates anything from tris to hexes. The callee has the
    /// responsibility to something reasonable with this.

    /// @param vert0   First vertex of the cell.
    /// @param vert1   Second vertex of the cell.
    double
    eval(const Vert* vert0, const Vert* vert1, const Vert* vert2,
	 const Vert* vert3 = NULL, const Vert* vert4 = NULL, const Vert* vert5 =
	 NULL,
	 const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    double
    eval(const Cell* pC) const;

//  To evaluate quality in metric space, a metric must be previously set. Not all quality measures are
//  able to be evaluated in metric space
    void
    vEvalInMetric(bool bMetric)
    {
      bEvalMetric = bMetric;
    }

    /// \brief Evaluate all values of the quality function for the cell.

    /// An alterantive for Quality in Metric Space, by vert Coords
    virtual void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const = 0;
    /// An alternate interface to support usage when verts are known
    /// rather than cells.
    ///
    /// Contrary to custom, the output variables are put first to
    /// accommodate a variable number of args.
    virtual void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;

    /// \brief Evaluate all values of the quality function for the cell.
    ///
    /// In many cases, there's only one value (for instance, aspect ratio
    /// measures.  But in other cases (for instance, dihedral angles),
    /// there's more than one.
    ///
    /// Contrary to custom, the output variables are put first for
    /// consistent ordering with the other evalAll call.
    ///
    /// @param nValues Number of quality values for this cell.
    /// @param values  Array in which quality values are returned.  The
    /// caller is responsible for ensuring there's enough space in values
    /// for all the data.
    /// @param pCell   Cell to check quality for.
    void
    evalAllFromCell(int& nValues, double* values, const Cell* pCell) const;

    /// \brief Return the number of separate values evalAll returns.
    virtual int
    maxNumValues() const = 0;
  };

  /*
   * 2D QualMeasures
   */

  class Angle2D : public QualMeasure {
  public:
    Angle2D() :
	QualMeasure("angles (2D)", 2)
    {
    }
    ~Angle2D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return true;
    }
  public:
    virtual void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    virtual void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 4;
    } // 4 for quads
  };

  class RadiusRatio2D : public QualMeasure {
  public:
    RadiusRatio2D() :
	QualMeasure("InToCircumRadius2D", 2)
    {
    }
    ~RadiusRatio2D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  class AreaPerim2D : public QualMeasure {
  public:
    AreaPerim2D() :
	QualMeasure("AreaPerim2D", 2)
    {
    }
    ~AreaPerim2D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  class EdgeLengthToRadius2D : public QualMeasure {
  public:
    EdgeLengthToRadius2D() :
	QualMeasure("EdgeLength to Radius (2D)", 2)
    {
    }
    ~EdgeLengthToRadius2D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 3;
    } // doesn't make sense for quads
  };

  /*
   * 3D QualMeasures
   */

  class DihedralAngles : public QualMeasure {
  public:
    DihedralAngles() :
	QualMeasure("dihedral angles", 3)
    {
    }
    ~DihedralAngles()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return true;
    }
  public:
    void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 12;
    }
  };

  class SineDihedralAngles : public QualMeasure {
  public:
    SineDihedralAngles() :
	QualMeasure("sine of dihedral angles", 3)
    {
    }
    ~SineDihedralAngles()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 12;
    }
  };

  class SolidAngles : public QualMeasure {
    DihedralAngles m_DA;
  public:
    SolidAngles() :
	QualMeasure("solid angles", 3), m_DA()
    {
    }
    ~SolidAngles()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 8;
    }
  };

  class ShortestToRadius3D : public QualMeasure {
  public:
    ShortestToRadius3D() :
	QualMeasure("ShortestToRadius3D", 3)
    {
    }
    ~ShortestToRadius3D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  class RadiusRatio3D : public QualMeasure {
  public:
    RadiusRatio3D() :
	QualMeasure("RadiusRatio3D", 3)
    {
    }
    ~RadiusRatio3D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 =
		      NULL,
		      const double *adLoc4 = NULL, const double *adLoc5 =
		      NULL,
		      const double *adLoc6 = NULL, const double *adLoc7 =
		      NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  /*
   * This class evaluates quality based on a volume length measure. It's an approach sugggested by
   * by Parthasarathy, Graichen, and Hathaway, but computed in the metric space. Volume is computed in metric space
   * and divided by the cube of its root means square length in metric space, then set so the highest quality is 0.
   */
  class VolumeLength3D : public QualMeasure {
  public:
    VolumeLength3D() :
	QualMeasure("VolumeLength3D", 3)
    {
    }
    ~VolumeLength3D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromVerts(int& nValues, double* values, const Vert* vert0,
		     const Vert* vert1, const Vert* vert2, const Vert* vert3 =
		     NULL,
		     const Vert* vert4 = NULL, const Vert* vert5 = NULL,
		     const Vert* vert6 = NULL, const Vert* vert7 = NULL) const;
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 = NULL, const double *adLoc4 = NULL,
		      const double *adLoc5 = NULL, const double *adLoc6 = NULL,
		      const double *adLoc7 = NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  class VolumeArea3D : public QualMeasure {
  public:
    VolumeArea3D() :
	QualMeasure("VolumeArea3D", 3)
    {
    }
    ~VolumeArea3D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 = NULL, const double *adLoc4 = NULL,
		      const double *adLoc5 = NULL, const double *adLoc6 = NULL,
		      const double *adLoc7 = NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };

  class MarcumSkewness3D : public QualMeasure {
  public:
    MarcumSkewness3D() :
	QualMeasure("MarcumSkewness3D", 3)
    {
    }
    ~MarcumSkewness3D()
    {
    }
  protected:
    virtual bool
    shouldMinimize() const
    {
      return false;
    }
  public:
    void
    evalAllFromCoords(int& nValues, double* values, const double *adLoc0,
		      const double *adLoc1, const double *adLoc2,
		      const double *adLoc3 = NULL, const double *adLoc4 = NULL,
		      const double *adLoc5 = NULL, const double *adLoc6 = NULL,
		      const double *adLoc7 = NULL) const;
    int
    maxNumValues() const
    {
      return 1;
    }
  };
}
// namespace GRUMMP

#endif
