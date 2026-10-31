#include "SMsmooth.h"

/* This is a very simple implementation of surface normals, for which
 the normal is constant.  This is suitable only for flat surfaces, in
 both 2D and 3D. */
static double stored_normal[MAX_DIM];

void
SMsetNormal(const double new_norm[MAX_DIM])
{
  SM_COPY_VECTOR(stored_normal, new_norm, MAX_DIM);
}

void
SMunitNormal(const double location[MAX_DIM], double normal[MAX_DIM])
{
  // This does nothing useful, except silencing a compiler warning.
  normal[0] = location[0];
  SM_COPY_VECTOR(normal, stored_normal, MAX_DIM);
}

