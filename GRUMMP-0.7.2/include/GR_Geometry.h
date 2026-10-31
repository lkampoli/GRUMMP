#ifndef GR_Geometry
#define GR_Geometry 1

#include "GR_config.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"

/// Distance between two points
inline double
dDist(const Vert * const pV0, const Vert * const pV1)
{
  assert(pV1->getSpaceDimen() == pV0->getSpaceDimen());
  if (pV0->getSpaceDimen() == 2) {
    return dDIST2D(pV0->getCoords(), pV1->getCoords());
  }
  else {
    return dDIST3D(pV0->getCoords(), pV1->getCoords());
  }
}

/// Normal of an edge in 2D
void
calcNormal2D(const double adA[2], const double adB[2], double adResult[2]);
/// Normal of a triangle in 2D
void
calcNormal2D(const double adA[2], const double adB[2], const double adC[2],
	     double * const adResult);
/// Normal of a triangle in 3D
void
calcNormal3D(const double adA[3], const double adB[3], const double adC[3],
	     double * const adResult);
/// Normal at the center of a bilinear quad in 3D
void
calcQuadNormal(const double pt0[3], const double pt1[3], const double pt2[3],
	       const double pt3[3], double norm[3]);
///
int
checkOrient3D(const Face* const pF, const double ad[3]);
///
int
checkOrient3D(const Face* const pF, const Vert* const pV);
///
int
checkOrient3D(const Vert* const pVA, const Vert* const pVB,
	      const Vert* const pVC, const Vert* const pVD);
///
int
checkOrient3D(const double adA[3], const double adB[3], const double adC[3],
	      const double adD[3]);
///
int
isInsphere(const Vert* const pVA, const Vert* const pVB, const Vert* const pVC,
	   const Vert* const pVD, const Vert* const pVE);
///
int
isInsphere(const double adA[3], const double adB[3], const double adC[3],
	   const double adD[3], const double adE[3]);
/// Is testPt inside the tet formed by ptA-ptD (right-handed order).
int
isInTet(const double testPt[3], const double ptA[3], const double ptB[3],
	const double ptC[3], const double ptD[3]);

/// Is testPt inside the tet?
int
isInTet(const double testPt[3], const Cell* const cell);

////
//int isCoplanarInTet(const Vert* const pV,     const Vert* const pVertA,
//		   const Vert* const pVertB, const Vert* const pVertC,
//		   const Vert* const pVertD);
///
int
isInTri(const Vert* const pV, const Vert* const pVertA,
	const Vert* const pVertB, const Vert* const pVertC);
///
int
isInTri(const double adVert[2], const double adVertA[2],
	const double adVertB[2], const double adVertC[2]);
///
int
checkOrient2D(const Vert* const pVA, const Vert* const pVB,
	      const Vert* const pVC);
///
int
checkOrient2D(const double adA[2], const double adB[2], const double adC[2]);
/// lets me do relevant stuff with surface meshes
/// pA, pB, pC all double arrays of size 3
/// pD, optional, but orientation only makes sense if projecting all.
//void projectTo2D(const Vert* const pVertA, const Vert* const pVertB,	const Vert* const pVertC,
//		double * pA, double * pB, double * pC, double * pD = NULL);
//void projectTo2D(const double adA[3],const double adB[3], const double adC[3],
//		double * pA, double * pB, double * pC, double * pD = NULL);
/// int iOrient2D(const EdgeFace* const pFE, const double ad[3]);
int
isInCircle(const Vert* const pVA, const Vert* const pVB, const Vert* const pVC,
	   const Vert* const pVD);
int
isInCircle(const double adA[2], const double adB[2], const double adC[2],
	   const double adD[2]);
///
//int isReflex(const EdgeFace* const pF);

///If line passing by pVA and pVB intersect line passing
///by pVC and pVD, return intersection coords.
bool
calcLineIntersection(const Vert* const pVA, const Vert* const pVB,
		     const Vert* const pVC, const Vert* const pVD,
		     double intersection[2]);

/// Does ray starting from pVA towards pVB intersect
/// line going from pVC to pVD.
/// return -1 if no intersection,
///         0 if intersect at either end of the line
///         1 if intersect the line.
/// This might be somewhat unstable, use with care...
//int intersectRays(const Vert* const pVA, const Vert* pVB,
//		  const Vert* const pVC, const Vert* pVD,
//		  double* intersection = NULL);

/// Does the segment going from pVA to pVB intersect
/// the segement going from pVC to pVD?
int
intersectSegments(const Vert* const pVA, const Vert* pVB, const Vert* const pVC,
		  const Vert* pVD);
double
dAreaTriFromEdgeLengths(double len0, double len1, double len2);
double
dAreaTri(const Vert* const pV0, const Vert* pV1, const Vert* const pV2);
double
dVolTet(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	const Vert* const pV3);
double
dInradius(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	  const Vert* const pV3);
double
dCircumradius(const Vert* const pV0, const Vert* pV1, const Vert* const pV2,
	      const Vert* const pV3);
/// calculates intersection between a line segment and a triangle, all defined as "verts"
/// used for surface insertion
/// return -1 if no intersection,
///         0 if co planar
///         1 if intersect the line.

/// either way, if the line (not segment) intersects the plane of the triangle, returns intersection on plane
int
calcSegmentTriangleIntersection(const Vert * const pV0, const Vert * const pV1,
				const Vert* const pFV0, const Vert* const pFV1,
				const Vert* const pFV2, double tol,
				double * intersection);

int
calcSegmentTriangleIntersection(const double l0[3], const double l1[3],
				const Vert* const pFV0, const Vert* const pFV1,
				const Vert* const pFV2, double tol,
				double * intersection);
// assumes one point is inside or on the circle. Returns true if intersection found, false if not
// this is NOT a generic function, it is specific to calculation of petals in optimized steiner insertion (See Erten and Ungor)
// this takes in the radius squared, avoiding a bunch of unnecessary square roots.
bool
calcSegmentCircleIntersections(const double adA[2], const double adB[2],
			       const double adCent[2], const double dRadius2,
			       double adIntersect[]);

void
calcCircumcenter3D(const Vert* const pVA, const Vert* pVB,
		   const Vert* const pVC, double * adCircCent);
// calculates area moment of a triangle/tetrahedron around its centroid
// really the covariance matrix
// adMoment[3] = int [ x^2, xy; xy, y^2] dA
// adMoment[5] = int [ x^2, xz, xy; xz, y^2, yz; xy, yz, z^2] dV
void
calcCovarianceMatrix(const Vert* const pV1, const Vert* pV2,
		     const Vert* const pV3, double * adMoment);
void
calcCovarianceMatrix(const Vert* const pV1, const Vert* pV2,
		     const Vert* const pV3, const Vert * const pV4,
		     double * adMoment);
int
closestOnTriangle(const double pt[3], const double t0[3], const double t1[3],
		  const double t2[3], double close[3]);

#endif
