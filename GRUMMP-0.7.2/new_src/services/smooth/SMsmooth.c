#include <stdio.h>
#include "SMsmooth.h"

/*@
 SMsmooth - This is the main routine that optimizes a local submesh.  Most of the quality
 functions available require that the local submesh is initially valid.  If the user suspects
 that the mesh contains invalid elements, a quality assessment should be done, and,
 if necessary, SMuntangle should be used to try to rectify the problem.

 Input Parameters:
 +  num_incident_vtx - the number of incident vertices in the local mesh
 .  num_tri - the number of incident triangles or tetrahedra
 .  free_vtx - the coordinates of the free vertex
 a vector of length equal to the problem dimension in x, y, z order
 .  vtx_list - a matrix of the coordinates of the incident vtx;
 matrix dimensions are num_incident_vtx by problem dimension
 .  vtx_connectivity - a matrix that gives the connectivity info for the
 incident vertices. matrix dimensions are num_incident_vtx by
 the problem dimension. Note: this assumes that the connectivity
 given is for a right handed triangle or tetrahedra the free vertex
 ordered first
 -  ext_smooth_data - data structure for mesh smoothing; created in
 SMinitSmoothing and cast to the local data structure

 Output Parameter:
 .  free_vtx -  contains the new coordinates of the free vertex

 Note:
 This function can be called only after SMinitSmoothing has been called to create
 the ext_smooth_data data structure.  Once the mesh has been optimized,
 SMfinalizeSmoothing should be called to release the memory.

 .seealso SMinitSmoothing(), SMsetSmoothTechnique(), SMsetSmoothFunction(),
 SMsetSmoothThreshold(), SMfinializeSmoothing()
 @*/
int
SMsmooth(int num_incident_vtx, int num_tri, double *free_vtx, double **vtx_list,
	 int **vtx_connectivity, void *ext_smooth_data, int is_surface_vert)
{
  int i, do_optimization = FALSE;
  int dimension;
  double min_val;
  SMsmooth_data *smooth_data;
  SMlocal_mesh *local_mesh;
  SMparam *smooth_param;
  SMstats *smooth_stats;
  SM_DEBUG_PRINT(2, "\nSmoothing a new local submesh\n");
  if (ext_smooth_data == NULL) {
    SM_ERROR(
	"The smoothing data structure has either not been created or has been corrupted\n",
	{
	});
    return (1);
  }

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  local_mesh = smooth_data->local_mesh;
  smooth_param = smooth_data->smooth_param;
  smooth_stats = smooth_data->smooth_stats;
  /* Set the flag for whether this is a surface smoothing case. */
  smooth_param->is_surface_vertex = is_surface_vert;
  /* initialize the local mesh triangles and the smoothing parameters */
  SMinitLocalMesh(num_incident_vtx, num_tri, free_vtx, vtx_list,
		  vtx_connectivity, local_mesh, smooth_param);
  dimension = local_mesh->dimension;
//  if (num_incident_vtx == (num_tri*2)){
//	  smooth_param->function_values_per_tri = 4;
//  }
//  else{
//	  smooth_param->function_values_per_tri = 3;
//  }
  /* check for some simple incorrect input problems */
  if (dimension == 2) {
    if (is_surface_vert && (num_incident_vtx != num_tri + 1)) {
      SM_ERROR("Incorrect input, num_incident_vtx != num_tri+1 on surface\n", {
      });
      return (1);
    }
    else if (!is_surface_vert && (num_incident_vtx != num_tri)) {
      SM_ERROR("Incorrect input, num_incident_vtx != num_tri\n", {
      });
      return (1);
    }
  }
  else if (dimension == 3) {
    if (is_surface_vert
	&& ((num_tri > 2 * num_incident_vtx - 3)
	    || (num_tri < num_incident_vtx - 2))) {
      SM_ERROR("Incorrect input, num_tri out of range on surface\n", {
      });
      printf("V: %2d T: %2d Comp1: %2d Comp2: %2d\n", num_incident_vtx, num_tri,
	     2 * num_incident_vtx - 3 - num_tri,
	     num_tri - (num_incident_vtx - 2));
      return (1);
    }
    if (!is_surface_vert && (num_tri != 2 * num_incident_vtx - 4)) {
      SM_ERROR("Incorrect input, 2*num_incident_vtx-4 != num_tri\n", {
      });
      return (1);
    }
  }

  /* if the technique chosen is Stupid Laplacian, do it and leave */
  if (smooth_param->smooth_technique == STUPID_LAPLACE) {
    SMcentroidSmoothMesh(local_mesh->num_incident_vtx, local_mesh->incident_vtx,
			 local_mesh->free_vtx, local_mesh->dimension,
			 is_surface_vert);
    SM_COPY_VECTOR(free_vtx, local_mesh->free_vtx, dimension);
    STATS_ON(
	{ SM_DEBUG_PRINT(2,"Creating stats \n"); SMaccumulateStats(local_mesh,smooth_param,smooth_stats); STATS_OFF});
    return (0); /* leaving after stupid laplacian smoothing */
  }

  /* check for an invalid mesh; if it's invalid return and ask the user to use SMuntangle */
  if (!SMvalidityCheck(local_mesh)) {
    printf(
	"The local mesh is invalid in SMsmooth; please check your input or use SMuntangle to create a valid mesh\n");
    assert(0);
    return (1); /* leaving after an invalid mesh was found */
  }
  /* The mesh is known to be valid at this point. */
  local_mesh->validity = OPTMS_VALID_MESH;

  /* compute the original function values */
  SMcomputeFunction(local_mesh, smooth_param, local_mesh->original_function);

  /* find the minimum value */
  min_val = 1.E300;
  for (i = 0; i < local_mesh->num_values; i++) {
    min_val = MIN(local_mesh->original_function[i], min_val);
  }
  local_mesh->original_value = local_mesh->current_active_value = min_val;

  SM_DEBUG_ACTION(
      3, {fprintf(stdout,"The initial minimum value is %f\n",min_val); });

  /**********************************************************************************************************
   Methods available in the smoothing code:

   LAPLACE_ONLY: take a laplacian step only if it improves the mesh and leave
   OPTIMIZATION_ONLY: for each local submesh do optimization based smoothing.
   This method is not recommended due to it's expense
   COMBINED 1: if the minimum value in the local submesh is too small, do
   optimization based smoothing, otherwise do Laplacian smoothing
   COMBINED 2: always do a Laplacian step; if the minimum value of the local
   submesh is too small following this step, perform optimization-based smoothing.
   This is the recommended technique as it gives general improvement to the mesh
   and eliminates extremal angles at roughly twice the cost of Laplacian smoothing
   used alone.
   COMBINED 3: if the minimum value of the local submesh exceeds the threshold value,
   do nothing, otherwise use the COMBINED 2 approach.  This is the cheapest of the
   three combined techniques and gives reasonable results (although quality metric
   distribution can clump near the threshold value.
   **********************************************************************************************************/
  switch (smooth_param->smooth_technique)
    {
    case LAPLACE_ONLY:
      SMinitLap(local_mesh->num_values, local_mesh->lap_info);
      SMlaplaceSmooth(local_mesh, smooth_param);
      do_optimization = FALSE;
      local_mesh->opt_info->opt_iter_count = 0;
      break;
    case OPTIMIZATION_ONLY:
      do_optimization = TRUE;
      break;
    case COMBINED1:
      if (smooth_param->lap_accept_value <= local_mesh->current_active_value) {
	/* do Laplacian */
	SMinitLap(local_mesh->num_values, local_mesh->lap_info);
	SMlaplaceSmooth(local_mesh, smooth_param);
	do_optimization = FALSE;
      }
      else {
	do_optimization = TRUE;
      }
      break;
    case COMBINED2:
      SMinitLap(local_mesh->num_values, local_mesh->lap_info);
      SMlaplaceSmooth(local_mesh, smooth_param);
      if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
	do_optimization = TRUE;
	SM_DEBUG_PRINT(
	    2, "Current Active value less than threshold, Doing Optimization\n");
      }
      else if (smooth_param->lap_accept_value
	  <= local_mesh->current_active_value) {
	do_optimization = FALSE;
	local_mesh->opt_info->opt_iter_count = 0;
	SM_DEBUG_PRINT(
	    2,
	    "Current Active Value exceeds threshold, Not Doing Optimization\n");
      }
      break;
    case COMBINED3:
      if (smooth_param->lap_accept_value <= local_mesh->current_active_value) {
	/* do nothing */
	do_optimization = FALSE;
      }
      else {
	/* do laplacian */
	SMinitLap(local_mesh->num_values, local_mesh->lap_info);
	SMlaplaceSmooth(local_mesh, smooth_param);
	/* decide if optimization is necessary */
	if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
	  do_optimization = TRUE;
	  SM_DEBUG_PRINT(
	      2,
	      "Current Active value less than threshold, Doing Optimization\n");
	}
	else if (smooth_param->lap_accept_value
	    <= local_mesh->current_active_value) {
	  do_optimization = FALSE;
	  local_mesh->opt_info->opt_iter_count = 0;
	  SM_DEBUG_PRINT(
	      2,
	      "Current Active Value exceeds threshold, Not Doing Optimization\n");
	}
      }
      break;
    case COMBINED:
    case FLOATING_THRESHOLD:
      SMinitLap(local_mesh->num_values, local_mesh->lap_info);
      SMlaplaceSmooth(local_mesh, smooth_param);
      if (smooth_param->lap_accept_value > local_mesh->current_active_value) {
	do_optimization = TRUE;
	SM_DEBUG_PRINT(
	    2, "Current Active value less than threshold, Doing Optimization\n");
      }
      else if (smooth_param->lap_accept_value
	  <= local_mesh->current_active_value) {
	do_optimization = FALSE;
	local_mesh->opt_info->opt_iter_count = 0;
	SM_DEBUG_PRINT(
	    2,
	    "Current Active Value exceeds threshold, Not Doing Optimization\n");
      }
      break;
    default:
      SM_DEBUG_PRINT(1, "Invalid smoothing technique\n")
      ;
      return (1);
    }

  /* start the optimization process */

  if (do_optimization == TRUE) {
    local_mesh->opt_done = TRUE;
    /* initialize the optimization structure */
    SMinitOpt(local_mesh->num_values, smooth_param->maxit,
	      local_mesh->opt_info);
    local_mesh->opt_info->prev_active_values[0] = local_mesh->original_value;
    if ((local_mesh->lap_done == 0)
	|| (local_mesh->lap_info->lap_accepted == FALSE)) {
      SM_COPY_VECTOR(local_mesh->opt_info->function,
		     local_mesh->original_function, local_mesh->num_values);
      local_mesh->opt_info->iter_count = 0;
    }
    else {
      SM_COPY_VECTOR(local_mesh->opt_info->function,
		     local_mesh->lap_info->laplacian_function,
		     local_mesh->num_values);
      local_mesh->opt_info->prev_active_values[1] =
	  local_mesh->lap_info->laplacian_value;
      local_mesh->opt_info->iter_count = 1;
    }
    SMfindActiveSet(local_mesh->num_values, local_mesh->opt_info->function,
		    smooth_param->active_eps, local_mesh->opt_info->active);
    SM_DEBUG_ACTION(
	3,
	{
	/* Print the active set */
	SMprintActiveSet(local_mesh->opt_info->active, local_mesh->opt_info->function); });

    /* terminate if we have found an equilibrium point or if the step is
     too small to be worthwhile continuing */
    while ((local_mesh->opt_info->status != EQUILIBRIUM)
	&& (local_mesh->opt_info->status != STEP_TOO_SMALL)
	&& (local_mesh->opt_info->status != IMP_TOO_SMALL)
	&& (local_mesh->opt_info->status != FLAT_NO_IMP)
	&& (local_mesh->opt_info->status != ZERO_SEARCH)
	&& (local_mesh->opt_info->status != MAX_ITER_EXCEEDED)) {

      /* increase the iteration count by one */
      /* smooth_param->iter_count += 1; */
      local_mesh->opt_info->iter_count += 1;
      local_mesh->opt_info->opt_iter_count += 1;
      if (local_mesh->opt_info->iter_count > MAX_SM_ITER)
	local_mesh->opt_info->status = MAX_ITER_EXCEEDED;

      SM_DEBUG_PRINT(3, "\n");
      SM_DEBUG_ACTION(
	  3,
	  { fprintf(stdout,"ITERATION %d \n",local_mesh->opt_info->iter_count); });

      /* compute the gradient */
      SMcomputeGradient(local_mesh, smooth_param,
			local_mesh->opt_info->gradient);

      SM_DEBUG_PRINT(3, "computing the search direction \n");
      SMsearchDirection(local_mesh);

      /* if there are viable directions to search */
      if ((local_mesh->opt_info->status != ZERO_SEARCH)
	  && (local_mesh->opt_info->status != MAX_ITER_EXCEEDED)) {

	SM_DEBUG_PRINT(3, "Computing the projections of the gradients \n");
	SMgetGradientProjections(local_mesh->opt_info);

	SM_DEBUG_PRINT(3, "Computing the initial step size \n");
	SMcomputeAlpha(local_mesh->opt_info);

	SM_DEBUG_PRINT(3, "Testing whether to accept this step \n");
	SMstepAcceptance(local_mesh, smooth_param);

	SM_DEBUG_ACTION(
	    3,
	    {
	    /* Print the active set */
	    SMprintActiveSet(local_mesh->opt_info->active, local_mesh->opt_info->function); });

	/* check for equilibrium point */
	if (local_mesh->opt_info->active->num_active >= 2) {
	  SM_DEBUG_PRINT(3, "Testing for an equilibrium point \n");
	  local_mesh->opt_info->equilibrium_pt = SMcheckEquilibrium(
	      local_mesh->opt_info);

	  SM_DEBUG_ACTION(
	      2,
	      { if (local_mesh->opt_info->equilibrium_pt) fprintf(stdout,"Optimization Exiting: An equilibrium point \n"); });
	}

	/* record the values */
	local_mesh->current_active_value =
	    local_mesh->opt_info->active->true_active_value;
	SM_RECORD_ITER_VALUE(local_mesh->opt_info);
      }
      else {
	/* decrease the iteration count by one */
	/* smooth_param->iter_count -= 1; */
	local_mesh->opt_info->iter_count -= 1;
	SM_DEBUG_ACTION(
	    2,
	    { fprintf(stdout,"Optimization Exiting: No viable directions; equilibrium point \n");
	    /* Print the old active set */
	    SMprintActiveSet(local_mesh->opt_info->active, local_mesh->opt_info->function); });
      }
    }
#ifndef NDEBUG
    SM_DEBUG_PRINT(3, "Checking the validity of the mesh\n");
    if (!SMvalidityCheck(local_mesh))
      fprintf(stdout, "The final mesh is not valid\n");
#endif
    SM_DEBUG_ACTION(
	2,
	{fprintf(stdout,"Number of optimization iterations %d\n",local_mesh->opt_info->iter_count);});
  }

  /* checking the global minimum value */
  if (smooth_param->global_min_value > local_mesh->current_active_value) {
    smooth_param->global_min_value = local_mesh->current_active_value;
  }

  SM_COPY_VECTOR(free_vtx, local_mesh->free_vtx, dimension);
  SM_DEBUG_PRINT(2, "Returning to example program\n");

  STATS_ON(
      { SM_DEBUG_PRINT(2,"Creating stats \n"); SMaccumulateStats(local_mesh,smooth_param,smooth_stats); STATS_OFF});

  return (0);
}

