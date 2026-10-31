#ifndef GR_LENGTHANISO2D_H
#define GR_LENGTHANISO2D_H

#include "GR_LengthAniso.h"

class LengthAniso2D : public LengthAniso {

private:

  LengthAniso2D(); // Never defined
  //Length2D(const Length2D&); // Never defined
  LengthAniso2D&
  operator=(const LengthAniso2D&); // Never defined

  virtual void
  evaluateLengthFromNeighbors(Vert* const vertex, double* adMetric,
			      const std::set<Vert*>& neigh_verts) const;

  virtual void
  evaluateLengthFromMesh(Vert* const vertex, double* adMetric) const;

public:

  //Public constructor.
  //length_type: the type of length scale used:  eAutomatic, eTree.
  //refine: the refinement constant,
  //grade: the grading constant
  LengthAniso2D(const eLengthType& lengthType, double refine, double grade,
		double min_length = LARGE_DBL) :
      LengthAniso(lengthType, refine, grade, min_length)
  {
  }

  //Provides a background mesh for length scale interpolation.
  virtual void
  provideMesh(Mesh* pM2D, char* strLenFileName = NULL);
};

#endif
