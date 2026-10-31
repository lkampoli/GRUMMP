#include "GR_DiscreteCurvTool.h"
#include "GR_Face.h"
#include "GR_SurfIntersect.h"
#include "GR_Vertex.h"

#include "CubitVector.hpp"

#include <algorithm>
#include <math.h>
#include <list>
#include <numeric>
#include <utility>

using std::list;

const double DiscreteCurvTool::ANGLE_COEFF = M_PI / 180;

DiscreteCurvTool::DiscreteCurvTool() :
    m_vertex(NULL), m_area(0), m_star_angles(), m_cotan()
{
}

DiscreteCurvTool*
DiscreteCurvTool::instance()
{

  static DiscreteCurvTool DCT;
  return &DCT;

}
DiscreteCurvTool::~DiscreteCurvTool()
{
//	delete m_vertex;
}
void
DiscreteCurvTool::init_tool(const Vert* const vertex,
			    const list<TriFace*>& star_triangles)
{

  m_area = 0.;
  m_star_angles.clear();
  m_cotan.clear();

  m_vertex = vertex;

  list<TriFace*>::const_iterator it = star_triangles.begin(), it_end =
      star_triangles.end();

  for (; it != it_end; ++it)
    init_area(*it);

}

void
DiscreteCurvTool::init_tool(const Vert* const vertex,
			    const list<Face*>& star_triangles)
{

  m_area = 0.;
  m_star_angles.clear();
  m_cotan.clear();

  m_vertex = vertex;

  list<Face*>::const_iterator it = star_triangles.begin(), it_end =
      star_triangles.end();

  for (; it != it_end; ++it)
    init_area(*it);

}

//Tool initializer designed to be used with the SurfTri
//class as defined in GR_SurfIntersect.

void
DiscreteCurvTool::init_tool(const Vert* const vertex,
			    const list<SurfTri*>& star_triangles)
{

  m_area = 0.;
  m_star_angles.clear();
  m_cotan.clear();

  m_vertex = vertex;

  SurfTri* triangle;
  list<SurfTri*>::const_iterator it = star_triangles.begin(), it_end =
      star_triangles.end();

  for (; it != it_end; ++it) {
    triangle = *it;
    init_area(triangle->get_vertex(0), triangle->get_vertex(1),
	      triangle->get_vertex(2));
  }

}

void
DiscreteCurvTool::init_tool(const Vert* const vertex,
			    const list<VertTriplet>& star_triangles)
{

  m_area = 0.;
  m_star_angles.clear();
  m_cotan.clear();

  m_vertex = vertex;

  VertTriplet triplet;
  list<VertTriplet>::const_iterator it = star_triangles.begin(), it_end =
      star_triangles.end();

  for (; it != it_end; ++it) {
    triplet = *it;
    init_area(triplet.v0, triplet.v1, triplet.v2);
  }

}

DiscreteCurvTool::Curvatures
DiscreteCurvTool::principal_curvatures(const double mean_curvature,
				       const double gaussian_curvature)
{

  double delta = mean_curvature * mean_curvature - gaussian_curvature;
  delta = delta < 0. ? 0. : sqrt(delta);

  return std::make_pair(mean_curvature + delta, mean_curvature - delta);

}

//Not implemented yet since there is no real need for principal
//directions with isotropic meshing. But if need be, refer to
//the paper cited in the header file (GR_SurfaceSampler.h).

DiscreteCurvTool::Curvatures
DiscreteCurvTool::principal_curvatures(CubitVector& /*princ_direction1*/,
				       CubitVector& /*princ_direction2*/) const
{
  assert(0);
  return principal_curvatures(mean_curvature(), gaussian_curvature());
}

double
DiscreteCurvTool::mean_curvature() const
{

  //Given correct initialization for m_vertex, 
  //this function returns the mean curvature at m_vertex.

  assert(!m_cotan.empty());

  const Vert* ring_vert;
  double cotan1, cotan2;

  CubitVector vertex_coord = CubitVector(m_vertex->getCoords());
  CubitVector vector_sum(0., 0., 0.);

  Cotangents::const_iterator it = m_cotan.begin(), it_end = m_cotan.end();

  do {

    assert(m_cotan.count(it->first) == 2);

    ring_vert = it->first;

    cotan1 =
	it->second.first ?
	    it->second.second : 1. / tan(ANGLE_COEFF * it->second.second);
    cotan2 =
	(++it)->second.first ?
	    it->second.second : 1. / tan(ANGLE_COEFF * it->second.second);

    vector_sum += (cotan1 + cotan2)
	* (vertex_coord - CubitVector(ring_vert->getCoords()));

  }
  while (++it != it_end);

  return (0.25 / m_area) * vector_sum.length();

}

double
DiscreteCurvTool::gaussian_curvature() const
{

  //Given correct initialization for m_vertex, 
  //this function returns the mean curvature at m_vertex.

  assert(!m_star_angles.empty());

  double angle_sum = std::accumulate(m_star_angles.begin(), m_star_angles.end(),
				     0.);

  return ((2. * M_PI) - (angle_sum * ANGLE_COEFF)) / m_area;

}

void
DiscreteCurvTool::init_area(const Face* const face)
{

  assert(!face->isDeleted());
  assert(face->getNumVerts() == 3);

  init_area(face->getVert(0), face->getVert(1), face->getVert(2));

}

void
DiscreteCurvTool::init_area(const Vert* const vert0, const Vert* const vert1,
			    const Vert* const vert2)
{

  assert(vert0 && vert1 && vert2);

  const Vert *v0, *v1, *v2;

  switch (m_vertex == vert0 ? 1 : m_vertex == vert1 ? 2 : 3)
    {
    case 1:
      v0 = vert0;
      v1 = vert1;
      v2 = vert2;
      break;
    case 2:
      v0 = vert1;
      v1 = vert0;
      v2 = vert2;
      break;
    case 3:
      v0 = vert2;
      v1 = vert0;
      v2 = vert1;
      break;
    default:
      vFatalError("", "");
      break;
    }

  assert(v0 == m_vertex);

  //Get the coords of the three vertices, compute the vectors
  //representing the triangle's edges and finally compute the 
  //three angles. 

  CubitVector coord[] =
    { CubitVector(v0->getCoords()), CubitVector(v1->getCoords()), CubitVector(
	v2->getCoords()) };

  CubitVector vec[] =
    { coord[1] - coord[0], coord[2] - coord[0], coord[2] - coord[1] };

  double angles[] =
    { vec[0].interior_angle(vec[1]), (-vec[0]).interior_angle(vec[2]),
	(-vec[1]).interior_angle(-vec[2]) };

  assert(
      iFuzzyComp(angles[0], 0.) == 1 && iFuzzyComp(angles[1], 0.) == 1
	  && iFuzzyComp(angles[2], 0.) == 1);
  assert(
      iFuzzyComp(angles[0], 180.) == -1 && iFuzzyComp(angles[1], 180.) == -1
	  && iFuzzyComp(angles[2], 180.) == -1);

  m_star_angles.push_back(angles[0]);

  //The "mixed" area (as defined in the paper) is computed.
  //If the triangle is non-obtuse, the mixed area = voronoi area.
  //Otherwise, approximate by a fraction of the triangle's area.

  if (angles[0] <= 90. && angles[1] <= 90. && angles[2] <= 90.) {

    double coeff[] =
      { 1. / tan(ANGLE_COEFF * angles[1]), 1. / tan(ANGLE_COEFF * angles[2]) };

    m_cotan.insert(std::make_pair(v1, std::make_pair(true, coeff[1])));
    m_cotan.insert(std::make_pair(v2, std::make_pair(true, coeff[0])));

    m_area += (vec[0].length_squared() * coeff[1]
	+ vec[1].length_squared() * coeff[0]) / 8.;
  }
  else {

    m_cotan.insert(std::make_pair(v1, std::make_pair(false, angles[2])));
    m_cotan.insert(std::make_pair(v2, std::make_pair(false, angles[1])));

    double vec1[] = adDIFF3D(vert0->getCoords(), vert1->getCoords()), vec2[] =
    adDIFF3D(vert0->getCoords(), vert2->getCoords()), vec3[3];
    vCROSS3D(vec1, vec2, vec3);

    double triangle_area = 0.5 * dMAG3D(vec3);

    if (angles[0] > 90.)
      m_area += 0.50 * triangle_area;
    else
      m_area += 0.25 * triangle_area;

  }

}

void
DiscreteCurvTool::assert_valid_star() const
{

  assert(0);

}
