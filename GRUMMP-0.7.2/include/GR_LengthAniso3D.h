#ifndef GR_LENGTHANISO3D_H
#define GR_LENGTHANISO3D_H

#include "GR_LengthAniso.h"

class LengthAniso3D : public LengthAniso {

private:

  LengthAniso3D(); // Never defined
  //Length2D(const Length2D&); // Never defined
  LengthAniso3D&
  operator=(const LengthAniso3D&); // Never defined

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
  LengthAniso3D(const eLengthType& lengthType, double refine, double grade,
		double min_length = LARGE_DBL) :
      LengthAniso(lengthType, refine, grade, min_length)
  {
  }

  //Provides a background mesh for length scale interpolation.
  virtual void
  provideMesh(Mesh* pM3D, char* strLenFileName = NULL);
};

#endif
