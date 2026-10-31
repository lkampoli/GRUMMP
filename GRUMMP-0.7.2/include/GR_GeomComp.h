#ifndef GR_GEOMCOMP_H
#define GR_GEOMCOMP_H

#include "GR_BFace.h"
#include "GR_Geometry.h"
#include "GR_Face.h"

#include "CubitVector.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
// This will be re-purposed into a set of functions that tie in with CGM,
// and do things that CGM can't do geometry related

// Challenge of this is that theres a fair amount of bouncing between
// double arrays of size 3, and CubitVectors

// Struct used for sorting CubitVectors
struct CubitVectorCompare {
  bool
  operator()(const CubitVector& lhs, const CubitVector& rhs) const
  {
    for (int i = 0; i < 3; i++) {
      if (iFuzzyComp(lhs[i], rhs[i]) != 0)
	return lhs[i] < rhs[i];
    }
    return false;
  }

};
bool
computeRefEdgeIntersectionsWithBFace(RefEdge * curve, const double * coordsA,
				     const double * coordsB, BFace * bface,
				     CubitVector & intersectCubitVec);

bool
computeRefFaceIntersectionsWithBFace(RefFace * refFace, const double * coordsA,
				     const double * coordsB,
				     const BFace * bface,
				     CubitVector & intersectCubitVec);

void
projectPointOnRefFaceToBdryFace(const BFace * bdryFace, RefFace * refFace,
				const CubitVector & point,
				CubitVector & projection);
#endif
