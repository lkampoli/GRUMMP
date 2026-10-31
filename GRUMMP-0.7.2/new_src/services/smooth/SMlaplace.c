#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

void
SMlaplaceSmooth(SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
  int i, valid = 0;
  int num_values;
  double min_value, func_diff;
  double *laplace_function;

  local_mesh->lap_done = TRUE;

  num_values = local_mesh->num_values;

  SMcentroidSmoothMesh(local_mesh->num_incident_vtx, local_mesh->incident_vtx,
		       local_mesh->free_vtx, local_mesh->dimension,
		       smooth_param->is_surface_vertex);

  /* check the validity of the new point */
  valid = SMvalidMesh(local_mesh);
  if (!valid && (local_mesh->validity == OPTMS_VALID_MESH)) {

    /* don't use this step if this step makes a previously valid mesh invalid*/
    SM_DEBUG_PRINT(2, "Did Not Accept Laplacian Smoothing\n");
    SM_COPY_VECTOR(local_mesh->free_vtx, local_mesh->original_pt,
		   local_mesh->dimension);
    /* write this out so that we can keep stats */
    local_mesh->lap_info->lap_invalid = TRUE;
    local_mesh->lap_info->lap_accepted = FALSE;
    local_mesh->lap_info->laplacian_value = local_mesh->original_value;
  }
  else {

    /* compute the new function values */
    SMcomputeFunction(local_mesh, smooth_param,
		      local_mesh->lap_info->laplacian_function);
    laplace_function = local_mesh->lap_info->laplacian_function;

    /* find the minimum value */
    min_value = 1E300;
    for (i = 0; i < num_values; i++) {
      min_value = MIN(laplace_function[i], min_value);
    }

    /* if this step improves the function, keep it */
    func_diff = min_value - local_mesh->original_value;

    SM_DEBUG_ACTION(
	2,
	{ fprintf(stdout,"The improvement from Laplacian smoothing is %f\n", func_diff);});
    if (func_diff > 0) {

      SM_DEBUG_PRINT(2, "Accepted Laplacian Smoothing \n");
      local_mesh->lap_info->lap_accepted = TRUE;

      /* if the mesh is now has gone from invalid to valid change the status
       of mesh->validity */
      if (valid)
	local_mesh->validity = OPTMS_VALID_MESH;

      /* set the iteration counter to 1 */
      smooth_param->iter_count += 1;
      local_mesh->lap_info->lap_improvement = func_diff;
      local_mesh->lap_info->laplacian_value = min_value;
      local_mesh->current_active_value = min_value;

    }
    else {
      /* don't use this step */
      SM_DEBUG_PRINT(2, "Did Not Accept Laplacian Smoothing\n");
      local_mesh->lap_info->lap_accepted = FALSE;
      local_mesh->lap_info->laplacian_value = local_mesh->original_value;
      SM_COPY_VECTOR(local_mesh->lap_info->laplacian_function,
		     local_mesh->original_function, local_mesh->num_values);
      SM_COPY_VECTOR(local_mesh->free_vtx, local_mesh->original_pt,
		     local_mesh->dimension);
    }
  }
}

void
SMcentroidSmoothMesh(int num_incident_vtx, double **incident_vtx,
		     double *free_vtx, int dimension,
		     const int is_surface_vertex)
{
  int i, j;
  double avg[MAX_DIM], movement[MAX_DIM];

  for (j = 0; j < dimension; j++) {
    avg[j] = 0.;
    for (i = 0; i < num_incident_vtx; i++)
      avg[j] += incident_vtx[i][j];
    movement[j] = (free_vtx[j] - avg[j] / num_incident_vtx);
  }
  /* First cut at surface smoothing. CFO-G, 4/98. */
  /* Project onto the tangent to the surface at this point. */
  if (is_surface_vertex) {
    /*   Find the surface unit normal */
    double unit_norm[MAX_DIM], dot_prod = 0; /* Init for picky compilers. */
    SMunitNormal(free_vtx, unit_norm);
    /*   Remove the normal component of the search direction */
    SM_DOT(dot_prod, unit_norm, movement, dimension);
    movement[0] -= dot_prod * unit_norm[0];
    movement[1] -= dot_prod * unit_norm[1];
    if (dimension == 3)
      movement[2] -= dot_prod * unit_norm[2];
  }
  for (j = 0; j < dimension; j++) {
    free_vtx[j] += movement[j];
    SM_DEBUG_ACTION(
	3, { fprintf(stdout,"final --  avg[%d] = %f\n",j, free_vtx[j]); });
  }
  /* End of changes for surface smoothing. CFO-G, 4/98. */
}

