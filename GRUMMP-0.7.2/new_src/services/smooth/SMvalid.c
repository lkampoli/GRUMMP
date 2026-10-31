#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

int
SMvalidityCheck(SMlocal_mesh *local_mesh)
{
  /*  validity tests : 1. the mesh is still valid based on righthandedness
   2. the new mesh is indeed better than the original
   */

  int valid = 1;

  /* check that the mesh is still valid, based on right handedness. 
   Returns a 1 or a 0 */
  valid = SMvalidMesh(local_mesh);
//  assert(valid);

  /* check to see that the mesh didn't get worse */
  if (local_mesh->original_value > local_mesh->current_active_value) {
    SM_DEBUG_ACTION(
	1,
	{ printf("The local mesh got worse; initial value %f; final value %f\n", local_mesh->original_value, local_mesh->current_active_value ); });
    valid = 0;
  }

  SM_DEBUG_ACTION(
      1,
      { if (!valid) { int i; printf("The mesh is no longer valid\n"); WRITE_BINARY_ORDERED_PTS(local_mesh); for (i=0;i<local_mesh->dimension;i++) fprintf(stderr,"free_new[%d] = %f;  ",i,local_mesh->free_vtx[i]); fprintf(stderr,"\n");

      PRINT_ORDERED_PTS(local_mesh); } });
  return (valid);
}

int
SMvalidMesh(SMlocal_mesh *local_mesh)
{
  int i, valid;
  int num_tri;
  int ind1, ind2, ind3;
  double *free_vtx;

  valid = 1;
  num_tri = local_mesh->num_tri;
  free_vtx = local_mesh->free_vtx;

  if (local_mesh->dimension == 2) {
    for (i = 0; i < num_tri && valid; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      valid = SMorient2D(free_vtx, local_mesh->incident_vtx[ind1],
			 local_mesh->incident_vtx[ind2]);
    }
  }

  if (local_mesh->dimension == 3) {
    for (i = 0; i < num_tri && valid; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      ind3 = local_mesh->vtx_connectivity[i][2];
      valid = (SMorient3D(local_mesh->incident_vtx[ind1],
			  local_mesh->incident_vtx[ind2],
			  local_mesh->incident_vtx[ind3], free_vtx) == 1);
    }
  }
  return (valid);
}

int
SMorient2D(double *vtx1, double *vtx2, double *vtx3)
{
  int valid;
  double cross;

  valid = 1;
  cross = (vtx2[XDIR] - vtx1[XDIR]) * (vtx3[YDIR] - vtx2[YDIR])
      - (vtx3[XDIR] - vtx2[XDIR]) * (vtx2[YDIR] - vtx1[YDIR]);

  if (cross < 0)
    valid = 0;
  return (valid);
}

int
SMorient3D(double *vtx1, double *vtx2, double *vtx3, double *free_vtx)
/* Returns 1 if tet ABCD is right-handed, -1 if it's left-handed,
 and 0 if it's essentially a tie. */
{
  double dXa = vtx1[0];
  double dYa = vtx1[1];
  double dZa = vtx1[2];

  double dXb = vtx2[0];
  double dYb = vtx2[1];
  double dZb = vtx2[2];

  double dXc = vtx3[0];
  double dYc = vtx3[1];
  double dZc = vtx3[2];

  double dXd = free_vtx[0];
  double dYd = free_vtx[1];
  double dZd = free_vtx[2];

  double dDX2 = dXb - dXa;
  double dDX3 = dXc - dXa;
  double dDX4 = dXd - dXa;

  double dDY2 = dYb - dYa;
  double dDY3 = dYc - dYa;
  double dDY4 = dYd - dYa;

  double dDZ2 = dZb - dZa;
  double dDZ3 = dZc - dZa;
  double dDZ4 = dZd - dZa;

  /* dDet is proportional to the cell volume */
  double dDet = dDX2 * dDY3 * dDZ4 + dDX3 * dDY4 * dDZ2 + dDX4 * dDY2 * dDZ3
      - dDZ2 * dDY3 * dDX4 - dDZ3 * dDY4 * dDX2 - dDZ4 * dDY2 * dDX3;

  /* Compute a length scale based on edge lengths. */
  double dScale = (sqrt(
      (dXa - dXb) * (dXa - dXb) + (dYa - dYb) * (dYa - dYb)
	  + (dZa - dZb) * (dZa - dZb))
      + sqrt(
	  (dXa - dXc) * (dXa - dXc) + (dYa - dYc) * (dYa - dYc)
	      + (dZa - dZc) * (dZa - dZc))
      + sqrt(
	  (dXa - dXd) * (dXa - dXd) + (dYa - dYd) * (dYa - dYd)
	      + (dZa - dZd) * (dZa - dZd))
      + sqrt(
	  (dXb - dXc) * (dXb - dXc) + (dYb - dYc) * (dYb - dYc)
	      + (dZb - dZc) * (dZb - dZc))
      + sqrt(
	  (dXb - dXd) * (dXb - dXd) + (dYb - dYd) * (dYb - dYd)
	      + (dZb - dZd) * (dZb - dZd))
      + sqrt(
	  (dXc - dXd) * (dXc - dXd) + (dYc - dYd) * (dYc - dYd)
	      + (dZc - dZd) * (dZc - dZd))) / 6.;

  double dEps = 1.e-13;

  /* Use the length scale to get a better idea if the tet is flat or
   just really small. */
  dDet /= (dScale * dScale * dScale);

  if (dDet > dEps)
    return (1);
  else if (dDet < -dEps)
    return (-1);
  else
    return (0);
}
