#ifndef GR_DISCRETECURVTOOL_H
#define GR_DISCRETECURVTOOL_H

#include "GR_config.h"

#include <math.h>
#include <list>
#include <map>
#include <utility>
#include <vector>

class CubitVector;
class Face;
class SurfTri;
class TriFace;
class Vert;

class DiscreteCurvTool {

  typedef std::vector<double> StarAngles;
  typedef std::multimap<const Vert*, std::pair<bool, double> > Cotangents;

  static const double ANGLE_COEFF;

  const Vert* m_vertex;
  double m_area;
  StarAngles m_star_angles;
  Cotangents m_cotan;

private:

  DiscreteCurvTool();
  DiscreteCurvTool(const DiscreteCurvTool&);
  DiscreteCurvTool&
  operator=(const DiscreteCurvTool&);

public:

  virtual
  ~DiscreteCurvTool();

  struct VertTriplet {

    Vert* v0, *v1, *v2;

    VertTriplet() :
	v0(NULL), v1(NULL), v2(NULL)
    {
    }

    VertTriplet(Vert* vv0, Vert* vv1, Vert* vv2) :
	v0(vv0), v1(vv1), v2(vv2)
    {
    }

    VertTriplet(const VertTriplet& triplet) :
	v0(triplet.v0), v1(triplet.v1), v2(triplet.v2)
    {
    }

    VertTriplet&
    operator=(const VertTriplet& triplet)
    {
      if (this != &triplet) {
	v0 = triplet.v0;
	v1 = triplet.v1;
	v2 = triplet.v2;
      }
      return *this;
    }

    void
    set(int i, Vert* vertex)
    {
      assert(i == 0 || i == 1 || i == 2);
      switch (i)
	{
	case 0:
	  v0 = vertex;
	  break;
	case 1:
	  v1 = vertex;
	  break;
	case 2:
	  v2 = vertex;
	  break;
	}
    }

  };

  typedef std::pair<double, double> Curvatures;

  static DiscreteCurvTool*
  instance();

  static Curvatures
  principal_curvatures(const double mean_curvature,
		       const double gaussian_curvature);

  void
  init_tool(const Vert* const vertex,
	    const std::list<TriFace*>& star_triangles);
  void
  init_tool(const Vert* const vertex, const std::list<Face*>& star_triangles);
  void
  init_tool(const Vert* const vertex,
	    const std::list<SurfTri*>& star_triangles);
  void
  init_tool(const Vert* const vertex,
	    const std::list<VertTriplet>& star_triangles);

  double
  mean_curvature() const;

  double
  mean_curvature(const Vert* const vertex,
		 const std::list<TriFace*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return mean_curvature();
  }
  double
  mean_curvature(const Vert* const vertex,
		 const std::list<Face*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return mean_curvature();
  }

  double
  gaussian_curvature() const;

  double
  gaussian_curvature(const Vert* const vertex,
		     const std::list<TriFace*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return gaussian_curvature();
  }
  double
  gaussian_curvature(const Vert* const vertex,
		     const std::list<Face*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return gaussian_curvature();
  }

  Curvatures
  principal_curvatures() const
  {
    return principal_curvatures(mean_curvature(), gaussian_curvature());
  }

  Curvatures
  principal_curvatures(const Vert* const vertex,
		       const std::list<TriFace*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return principal_curvatures();
  }
  Curvatures
  principal_curvatures(const Vert* const vertex,
		       const std::list<Face*>& star_triangles)
  {
    init_tool(vertex, star_triangles);
    return principal_curvatures();
  }

private:

  void
  init_area(const Face* const face);

  void
  init_area(const Vert* const vert0, const Vert* const vert1,
	    const Vert* const vert2);

  //Not implemented. Keep private for now.
  Curvatures
  principal_curvatures(CubitVector& princ_direction1,
		       CubitVector& princ_direction2) const;

  void
  assert_valid_star() const;

};

#endif
