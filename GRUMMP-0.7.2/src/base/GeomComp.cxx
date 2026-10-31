#include "GR_GeomComp.h"

// Refactoring this to reuse for functions related to CGM,
// D. Zaide, Oct 2014

bool
computeRefEdgeIntersectionsWithBFace(RefEdge * curve, const double * coordsA,
				     const double * coordsB, BFace * bface,
				     CubitVector & intersectCubitVec)
{

  // check if the edge we are on has verts on either side first
  // obviously fails if the curve cuts through and then doubles back,
  // but I'm not worried about those cases
  if (checkOrient3D(bface->getVert(0)->getCoords(),
		    bface->getVert(1)->getCoords(),
		    bface->getVert(2)->getCoords(), coordsA)
      * checkOrient3D(bface->getVert(0)->getCoords(),
		      bface->getVert(1)->getCoords(),
		      bface->getVert(2)->getCoords(), coordsB) == 1)
    return false;

  CubitVector pointA(coordsA), pointB(coordsB);
  double faceNormal[3];
  bface->calcUnitNormal(faceNormal);
  CubitVector normalCubitVec(faceNormal), faceCubitVec(
      bface->getVert(0)->getCoords());
  CubitVector tangentA, tangentB;

  double dist = (pointA - pointB).length();
  double distA, distB;
  // basic approach is to walk along the curve in both directions until we have two coordinates very close to eachother
  // on either side of the internal boundary, and then check for intersection.
  // this gets around the case where the curve wraps around a boundary face.
  GR_index_t coordIter = 0;
  while (dist > 1.e-6 && ++coordIter < 100) {
    // move the points closer before checking.
    curve->tangent(pointA, tangentA);
    curve->tangent(pointB, tangentB);
    tangentA.normalize();
    tangentB.normalize();
    // find intersection between tangent with plane
    distA = ((faceCubitVec - pointA) % normalCubitVec)
	/ (tangentA % normalCubitVec);
    distB = ((faceCubitVec - pointB) % normalCubitVec)
	/ (tangentB % normalCubitVec);
    pointA = pointA + distA / 2. * tangentA;
    pointB = pointB + distB / 2. * tangentB;
    // move the points half the distance closer, and
    // keep iterating
    curve->closest_point_trimmed(pointA, pointA);
    curve->closest_point_trimmed(pointB, pointB);
    dist = (fabs(distA) + fabs(distB));

  }

  // bring the coords close together.
  // check the linear intersection
  double intersection[3];
  if (calcSegmentTriangleIntersection(coordsA, coordsB, bface->getVert(0),
				      bface->getVert(1), bface->getVert(2), 0.,
				      intersection) == 1) {
    // iterate a few times, projecting back and forth to get
    // a little bit closer
    double projection[3];
    CubitVector CVprojection, CVintersect;
    double dDiff = LARGE_DBL;
    int iter = 0;
    while (dDiff > 1e-10 && ++iter < 10) {
      CVintersect.set(intersection);
      curve->closest_point_trimmed(CVintersect, CVprojection);
      CVprojection.get_xyz(projection);
      bface->findClosestPointOnBFace(projection, intersection);
    }
    intersectCubitVec.set(intersection);
    return true;
  }
  return false;
}
;
bool
computeRefFaceIntersectionsWithBFace(RefFace * refFace, const double * coordsA,
				     const double * coordsB,
				     const BFace * bface,
				     CubitVector & intersectCubitVec)
{

  // if they are on the same side as the planar face, assume it does not intersect
  if (checkOrient3D(bface->getVert(0)->getCoords(),
		    bface->getVert(1)->getCoords(),
		    bface->getVert(2)->getCoords(), coordsA)
      * checkOrient3D(bface->getVert(0)->getCoords(),
		      bface->getVert(1)->getCoords(),
		      bface->getVert(2)->getCoords(), coordsB) == 1)
    return false;

  CubitVector pointA(coordsA), pointB(coordsB);
  CubitVector normalA = refFace->normal_at(pointA);
  CubitVector normalB = refFace->normal_at(pointB);

  double bfaceNormal[3];
  bface->calcUnitNormal(bfaceNormal);
  double dist = (pointA - pointB).length();
  double distA, distB;
  CubitVector normalCV(bfaceNormal), bfaceCV(bface->getVert(0)->getCoords());
  // we should hit this in 100 iterations or less
  GR_index_t ii = 0;
  while (dist > 1e-6 && ++ii < 100) {
    // now attempt to find coordinates pointA, pointB,
    // closer to the plane but still on opposite sides

    //project the vector pointA-pointB onto the plane at pointA;
    CubitVector pointAproj = ~(pointB - pointA
	- ((pointB - pointA) % normalA) * normalA);
    //project the vector pointB-pointA onto the plane at pointB;
    CubitVector pointBproj = ~(pointA - pointB
	- ((pointA - pointB) % normalB) * normalB);

    // move pointA along this vector, half of the way.
    distA = ((bfaceCV - pointA) % normalCV) / (pointAproj % normalCV); //fabs(normalCV % (cpointA - bfaceCV));
    distB = ((bfaceCV - pointB) % normalCV) / (pointBproj % normalCV); //fabs(normalCV % (cpointB - bfaceCV));
    // new edge closer to the boundary face
    pointA = pointA + distA / 2. * pointAproj;
    pointB = pointB + distB / 2. * pointBproj;
    dist = fabs(distA) + fabs(distB);
    refFace->get_point_normal(pointA, normalA);
    refFace->get_point_normal(pointB, normalB);

  }

  double intersection[3], projection[3];
  CubitVector projectionCubitVec;

  // if we intersect on the face
  if (calcSegmentTriangleIntersection(coordsA, coordsB, bface->getVert(0),
				      bface->getVert(1), bface->getVert(2), 0.,
				      intersection) == 1) {

    double dDiff = LARGE_DBL;
    GR_index_t iter = 0;
    // again, should not iterate till 50
    while (dDiff > 1.e-8 && ++iter < 50) {
      // do your projection stuff
      intersectCubitVec.set(intersection);
      refFace->find_closest_point_trimmed(intersectCubitVec,
					  projectionCubitVec);
      projectionCubitVec.get_xyz(projection);
      closestOnTriangle(projection, bface->getVert(0)->getCoords(),
			bface->getVert(1)->getCoords(),
			bface->getVert(2)->getCoords(), intersection);
      dDiff = dDIST3D(intersection, projection);
    }

    // tolerance
    if (dDiff < 1.e-6)
      return true;
  }
  return false;

}
//void projectPointOnRefFaceToBdryFace(const BFace * bdryFace, RefFace * refFace, const CubitVector & point,
//		CubitVector & projection){
//
//}
