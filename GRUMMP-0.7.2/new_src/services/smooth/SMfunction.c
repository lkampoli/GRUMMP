#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

void
SMcomputeFunction(SMlocal_mesh *local_mesh, SMparam *smooth_param,
		  double *function)
{
  int fcn_ind, ind1, ind2, ind3;
  int i, j, num_values, num_tri;
  int dimension;

  SMfunction_ptr2D ComputeFunctionValues2D;
  SMfunction_ptr3D ComputeFunctionValues3D;

  dimension = local_mesh->dimension;
  if (dimension == 2) {
    double function1[3];

    ComputeFunctionValues2D = smooth_param->ComputeFunctionValues2D;
    fcn_ind = 0;
    ind1 = 0;
    ind2 = 0;
    ind3 = 0;

    num_tri = local_mesh->num_tri;
    for (i = 0; i < num_tri; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      ComputeFunctionValues2D(local_mesh->free_vtx,
			      local_mesh->incident_vtx[ind1],
			      local_mesh->incident_vtx[ind2], function1,
			      &num_values);

      for (j = 0; j < num_values; j++)
	function[fcn_ind++] = function1[j];
    }
  }
  else if (dimension == 3) {
    double function1[6];
    ComputeFunctionValues3D = smooth_param->ComputeFunctionValues3D;
    fcn_ind = 0;
    ind1 = 0;
    ind2 = 0;
    ind3 = 0;

    num_tri = local_mesh->num_tri;
    for (i = 0; i < num_tri; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      ind3 = local_mesh->vtx_connectivity[i][2];
      ComputeFunctionValues3D(local_mesh->free_vtx,
			      local_mesh->incident_vtx[ind1],
			      local_mesh->incident_vtx[ind2],
			      local_mesh->incident_vtx[ind3], function1,
			      &num_values);
      for (j = 0; j < num_values; j++)
	function[fcn_ind++] = function1[j];
    }
  }
}

void
SMcomputeGradient(SMlocal_mesh *local_mesh, SMparam *smooth_param,
		  double **gradient)
{
  int i, j, num_values, num_tri;
  int ind1, ind2, ind3, grad_ind;
  int dimension;

  dimension = local_mesh->dimension;

  if (dimension == 2) {
    double gradient1[3][2];
    SMgradfunc_ptr2D ComputeGradientValues2D =
	smooth_param->ComputeGradientValues2D;

    ind1 = 0;
    ind2 = 0;
    ind3 = 0;
    grad_ind = 0;
    num_tri = local_mesh->num_tri;
    for (i = 0; i < num_tri; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      ComputeGradientValues2D(local_mesh->free_vtx,
			      local_mesh->incident_vtx[ind1],
			      local_mesh->incident_vtx[ind2], gradient1,
			      &num_values);
      /* First cut at surface smoothing. CFO-G, 4/98. */
      /* Project onto the tangent line to the surface at this point. */
      if (smooth_param->is_surface_vertex) {
	/*   Find the surface unit normal */
	double unit_norm[MAX_DIM];
	SMunitNormal(local_mesh->free_vtx, unit_norm);
	/*   Remove the normal component of the search direction */
	for (j = 0; j < num_values; j++) {
	  double dot_prod;
	  SM_DOT(dot_prod, unit_norm, gradient1[j], 2);
	  gradient1[j][0] -= dot_prod * unit_norm[0];
	  gradient1[j][1] -= dot_prod * unit_norm[1];
	}
      }
      /* End of changes for surface smoothing. CFO-G, 4/98. */
      for (j = 0; j < num_values; j++) {
	SM_COPY_VECTOR2(gradient[grad_ind], gradient1[j]);
	grad_ind++;
      }
    }
  }
  else if (dimension == 3) {
    SMgradfunc_ptr3D ComputeGradientValues3D =
	smooth_param->ComputeGradientValues3D;
    double gradient1[6][3];

    ind1 = 0;
    ind2 = 0;
    ind3 = 0;
    grad_ind = 0;
    num_tri = local_mesh->num_tri;
    for (i = 0; i < num_tri; i++) {
      ind1 = local_mesh->vtx_connectivity[i][0];
      ind2 = local_mesh->vtx_connectivity[i][1];
      ind3 = local_mesh->vtx_connectivity[i][2];
      ComputeGradientValues3D(local_mesh->free_vtx,
			      local_mesh->incident_vtx[ind1],
			      local_mesh->incident_vtx[ind2],
			      local_mesh->incident_vtx[ind3], gradient1,
			      &num_values);
      /* First cut at surface smoothing. CFO-G, 4/98. */
      /* Project onto the tangent line to the surface at this point. */
      if (smooth_param->is_surface_vertex) {
	/*   Find the surface unit normal */
	double unit_norm[MAX_DIM];
	SMunitNormal(local_mesh->free_vtx, unit_norm);
	/*   Remove the normal component of the search direction */
	for (j = 0; j < num_values; j++) {
	  double dot_prod;
	  SM_DOT(dot_prod, unit_norm, gradient1[j], 3);
	  gradient1[j][0] -= dot_prod * unit_norm[0];
	  gradient1[j][1] -= dot_prod * unit_norm[1];
	  gradient1[j][2] -= dot_prod * unit_norm[2];
	}
      }
      /* End of changes for surface smoothing. CFO-G, 4/98. */
      for (j = 0; j < num_values; j++) {
	SM_COPY_VECTOR3(gradient[grad_ind], gradient1[j]);
	grad_ind++;
      }
    }
  }
}

void
SMcopyFunction(int num_values, double *function1, double *function2)
{
  int i;

  for (i = 0; i < num_values; i++) {
    function2[i] = function1[i];
  }
}

void
SMcopyActive(SMactive *active1, SMactive *active2)
{
  int i;

  active2->num_active = active1->num_active;
  active2->num_equal = active1->num_equal;
  active2->true_active_value = active1->true_active_value;
  for (i = 0; i < active1->num_active; i++) {
    active2->active_ind[i] = active1->active_ind[i];
  }
}
