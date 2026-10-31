#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

void
SMstepAcceptance(SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
  int i;
  int num_values, num_steps;
  int valid, step_status;
  int accept_alpha;
  int dimension;
  double alpha;
  double estimated_improvement;
  double current_improvement = A_BIG_NEG_NUMBER;
  double previous_improvement = A_BIG_NEG_NUMBER;
  double current_percent_diff = A_BIG_POS_NUMBER;
  double original_point[MAX_DIM];
  SMoptimal *opt_info;

  opt_info = local_mesh->opt_info;
  num_values = opt_info->num_values;
  alpha = opt_info->alpha;
  valid = 1; /* Initialize for picky compilers. */

  dimension = local_mesh->dimension;

  step_status = STEP_NOT_DONE;
  opt_info->status = 0;
  num_steps = 0;

  if (alpha < smooth_param->min_step_size) {
    opt_info->status = IMP_TOO_SMALL;
    step_status = STEP_DONE;
    SM_DEBUG_PRINT(3, "Alpha starts too small, no improvement\n");
  }
  if (alpha > 1E16) {
    SM_ERROR("Error: The initial step size\n", {
      SMwrite_ordered_points(local_mesh)
      ;
    });
    opt_info->status = IMP_TOO_SMALL;
    step_status = STEP_DONE;
  }

  /* save the original function and active set */
  SM_COPY_VECTOR(original_point, local_mesh->free_vtx, dimension);
  SM_COPY_VECTOR(opt_info->original_function, opt_info->function, num_values);
  SMcopyActive(opt_info->active, opt_info->original_active);

  while (step_status == STEP_NOT_DONE) {

    num_steps++;
    if (num_steps >= 100)
      step_status = STEP_DONE;

    accept_alpha = FALSE;
    while (!accept_alpha && alpha > smooth_param->min_step_size) {

      /* make the step */
      for (i = 0; i < dimension; i++) {
	local_mesh->free_vtx[i] += alpha * opt_info->search[i];
      }

      /* assume alpha is acceptable */
      accept_alpha = TRUE;

      /* check for valid step */
      valid = 1;
      valid = SMvalidMesh(local_mesh);

      /* never take a step that makes a valid mesh invalid */
      if (!valid && (local_mesh->validity == OPTMS_VALID_MESH)) {
	accept_alpha = FALSE;
	for (i = 0; i < dimension; i++) {
	  local_mesh->free_vtx[i] -= alpha * opt_info->search[i];
	}
	alpha = alpha / 2;
	opt_info->alpha = alpha;
	SM_DEBUG_ACTION(
	    2,
	    { fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); });

	if (alpha < smooth_param->min_step_size) {
	  opt_info->status = STEP_TOO_SMALL;
	  step_status = STEP_DONE;
	  SM_DEBUG_PRINT(2, "Step too small\n");
	  /* get back the original point, function, and active set */
	  SM_COPY_VECTOR(local_mesh->free_vtx, original_point, dimension);
	  SM_COPY_VECTOR(opt_info->function, opt_info->original_function,
			 num_values);
	  SMcopyActive(opt_info->original_active, opt_info->active);
	}
      }
    }

    if ((valid || local_mesh->validity == OPTMS_INVALID_MESH)
	&& (alpha > smooth_param->min_step_size)) {
      /* compute the new function and active set */
      SMcomputeFunction(local_mesh, smooth_param, opt_info->function);
      SMfindActiveSet(local_mesh->opt_info->num_values, opt_info->function,
		      smooth_param->active_eps, opt_info->active);

      /* estimate the minimum improvement by taking this step */
      estimated_improvement = SMgetMinEstimate(local_mesh->opt_info);
      SM_DEBUG_ACTION(
	  3,
	  { fprintf(stdout,"The estimated improvement for this step: %f\n", estimated_improvement); });

      /* calculate the actual increase */
      current_improvement = opt_info->active->true_active_value
	  - opt_info->prev_active_values[opt_info->iter_count - 1];

      SM_DEBUG_ACTION(
	  3, {fprintf(stdout,"Actual improvement %f\n",current_improvement);});

      /* calculate the percent difference from estimated increase */
      current_percent_diff = fabs(current_improvement - estimated_improvement)
	  / fabs(estimated_improvement);

      /* determine whether to accept a step */
      if ((previous_improvement > current_improvement)
	  && (previous_improvement > 0)) {
	/* accept the previous step - it was better */
	SM_DEBUG_PRINT(2, "Accepting the previous step\n");

	/* add alpha in again (previous step) */
	for (i = 0; i < dimension; i++) {
	  local_mesh->free_vtx[i] += alpha * opt_info->search[i];
	}

	/* does this make an invalid mesh valid? */
	valid = 1;
	valid = SMvalidMesh(local_mesh);
	if (valid && local_mesh->validity == OPTMS_INVALID_MESH) {
	  local_mesh->validity = OPTMS_VALID_MESH;
	}

	/* copy test function and active set */
	SM_COPY_VECTOR(opt_info->function, opt_info->test_function, num_values);
	SMcopyActive(opt_info->test_active, opt_info->active);

	opt_info->status = STEP_ACCEPTED;
	step_status = STEP_DONE;

	/* check to see that we're still making good improvements */
	if (previous_improvement < smooth_param->min_acceptable_imp) {
	  opt_info->status = IMP_TOO_SMALL;
	  step_status = STEP_DONE;
	  SM_DEBUG_PRINT(2, "Optimization Exiting: Improvement too small\n");
	}

      }
      else if (((current_improvement > estimated_improvement)
	  || (current_percent_diff < .1)) && (current_improvement > 0)) {
	/* accept this step, exceeded estimated increase or was close */
	opt_info->status = STEP_ACCEPTED;
	step_status = STEP_DONE;

	/* does this make an invalid mesh valid? */
	valid = 1;
	valid = SMvalidMesh(local_mesh);
	if (valid && local_mesh->validity == OPTMS_INVALID_MESH) {
	  local_mesh->validity = OPTMS_VALID_MESH;
	}

	/* check to see that we're still making good improvements */
	if (current_improvement < smooth_param->min_acceptable_imp) {
	  SM_DEBUG_PRINT(2, "Optimization Exiting: Improvement too small\n");
	  opt_info->status = IMP_TOO_SMALL;
	  step_status = STEP_DONE;
	}
      }
      else if ((current_improvement < 0) && (previous_improvement < 0)
	  && (fabs(current_improvement) < smooth_param->min_acceptable_imp)
	  && (fabs(previous_improvement) < smooth_param->min_acceptable_imp)) {

	/* we are making no progress, quit */
	opt_info->status = FLAT_NO_IMP;
	step_status = STEP_DONE;
	SM_DEBUG_PRINT(2, "Optimization Exiting: Flat no improvement\n");

	/* get back the original point, function, and active set */
	SM_COPY_VECTOR(local_mesh->free_vtx, original_point, dimension);
	SM_COPY_VECTOR(opt_info->function, opt_info->original_function,
		       num_values);
	SMcopyActive(opt_info->original_active, opt_info->active);

      }
      else {
	/* halve alpha and try again */
	/* subtract out the old step */
	for (i = 0; i < dimension; i++)
	  local_mesh->free_vtx[i] -= alpha * opt_info->search[i];

	/* halve step size */
	alpha = alpha / 2;
	opt_info->alpha = alpha;
	SM_DEBUG_ACTION(
	    3,
	    {fprintf(stdout,"Step not accepted, the new alpha %f\n",alpha); });

	if (alpha < smooth_param->min_step_size) {
	  /* get back the original point, function, and active set */
	  SM_DEBUG_PRINT(2, "Optimization Exiting: Step too small\n");
	  SM_COPY_VECTOR(local_mesh->free_vtx, original_point, dimension);
	  SM_COPY_VECTOR(opt_info->function, opt_info->original_function,
			 num_values);
	  SMcopyActive(opt_info->original_active, opt_info->active);
	  opt_info->status = STEP_TOO_SMALL;
	  step_status = STEP_DONE;
	}
	else {
	  SM_COPY_VECTOR(opt_info->test_function, opt_info->function,
			 num_values);
	  SMcopyActive(opt_info->active, opt_info->test_active);
	  previous_improvement = current_improvement;
	}
      }
    }
  }
}

int
SMcheckEquilibrium(SMoptimal *opt_info)
{
  int i, j;
  int equil;
  int num_active;
  int ind1, ind2;
  int dimension;
  double min;
  double **dir, **G;
  double mid_vec[MAX_DIM], mid_cos, test_cos;

  dimension = opt_info->dimension;

  num_active = opt_info->active->num_active;
  equil = 0;
  ind1 = ind2 = -1;

  if (num_active == opt_info->num_values) {
    equil = 1;
    opt_info->status = EQUILIBRIUM;
    SM_DEBUG_PRINT(3, "All the function values are in the active set\n");
    return equil;
  }

  /* set up the active gradient directions */
  dir = SMgetActiveDirections(num_active, opt_info->gradient,
			      opt_info->active->active_ind, dimension);

  /* normalize the active directions */
  for (j = 0; j < num_active; j++)
    SM_NORMALIZE(dir[j], dimension);

  if (dimension == 2) {
    /* form the grammian */
    (void) SMformGrammian(num_active, dir, opt_info->G, dimension);
    G = opt_info->G;

    /* find the minimum element in the upper triangular portion of G*/
    min = 1;
    for (i = 0; i < num_active; i++) {
      for (j = i + 1; j < num_active; j++) {
	if (fabs(-1 - G[i][j]) < MACHINE_EPS) {
	  equil = 1;
	  opt_info->status = EQUILIBRIUM;
	  SM_DEBUG_PRINT(3, "The gradients are antiparallel, eq. pt\n");
	  SMfreeActiveDirections(num_active, dir);
	  return equil;
	}
	if (G[i][j] < min) {
	  ind1 = i;
	  ind2 = j;
	  min = G[i][j];
	}
      }
    }

    if ((ind1 != -1) && (ind2 != -1)) {
      /* find the diagonal of the parallelepiped */
      for (j = 0; j < dimension; j++) {
	mid_vec[j] = .5 * (dir[ind1][j] + dir[ind2][j]);
      }
      SM_NORMALIZE(mid_vec, dimension);
      SM_DOT(mid_cos, dir[ind1], mid_vec, dimension);

      /* test the other vectors to be sure they lie in the cone */
      for (i = 0; i < num_active; i++) {
	if ((i != ind1) && (i != ind2)) {
	  SM_DOT(test_cos, dir[i], mid_vec, dimension);
	  if ((test_cos < mid_cos) && fabs(test_cos - mid_cos) > MACHINE_EPS) {
	    SM_DEBUG_PRINT(3, "An equilibrium point \n");
	    equil = 1;
	    opt_info->status = EQUILIBRIUM;
	    SMfreeActiveDirections(num_active, dir);
	    return equil;
	  }
	}
      }
    }
  }
  if (dimension == 3) {
    equil = SMconvexHullTest(dir, num_active);
    if (equil == 1)
      opt_info->status = EQUILIBRIUM;
  }
  SMfreeActiveDirections(num_active, dir);
  return equil;
}

double
SMgetMinEstimate(SMoptimal *opt_info)
{
  int i;
  int num_active;
  double alpha;
  double final_est, est_imp;

  num_active = opt_info->active->num_active;
  alpha = opt_info->alpha;

  final_est = A_BIG_POS_NUMBER;
  for (i = 0; i < num_active; i++) {
    est_imp = alpha * opt_info->gs[opt_info->active->active_ind[i]];
    if (est_imp < final_est)
      final_est = est_imp;
  }
  if (final_est == 0) {
    final_est = A_BIG_POS_NUMBER;
    for (i = 0; i < opt_info->num_values; i++) {
      est_imp = alpha * opt_info->gs[i];
      if ((est_imp < final_est) && (fabs(est_imp) > MACHINE_EPS)) {
	final_est = est_imp;
      }
    }
  }
  return final_est;
}

void
SMcomputeVerticalStep(SMlocal_mesh *local_mesh)
{
  int i;
  int *PDG_ind;
  int num_LI;
  int dimension;
  double **N, *R, *y, y_0;
  double v[MAX_DIM];
  double **gradient;
  SMoptimal *opt_info;

  opt_info = local_mesh->opt_info;
  num_LI = opt_info->num_LI;
  gradient = opt_info->gradient;
  PDG_ind = opt_info->PDG_ind;

  dimension = local_mesh->dimension;

  /* form the projection matrix */
  N = SMformVerticalMatrix(num_LI, opt_info->PDG);

  /* form the RHS */
  R = SMformVerticalRHS(num_LI, opt_info->function, opt_info->PDG_ind,
			opt_info->active->true_active_value);

  /* determine the step */
  switch (num_LI)
    {
    case 1:
      y_0 = R[0] / N[0][0];
      for (i = 0; i < dimension; i++) {
	v[i] = gradient[PDG_ind[0]][i] * y_0;
      }
      break;
    case 2:
      y = SMsolve2x2(N[0][0], N[0][1], N[1][0], N[1][1], R[0], R[1]);
      for (i = 0; i < dimension; i++) {
	v[i] = gradient[PDG_ind[0]][i] * y[0] + gradient[PDG_ind[1]][i] * y[1];
      }
      MY_FREE(y)
      ;
      break;
    case 3:
      y = SMsolveSPD3x3(N, R);
      for (i = 0; i < dimension; i++) {
	v[i] = gradient[PDG_ind[0]][i] * y[0] + gradient[PDG_ind[1]][i] * y[1]
	    + gradient[PDG_ind[2]][i] * y[2];
      }
      MY_FREE(y)
      ;
      break;
    }

  /* add in the correction */
  for (i = 0; i < dimension; i++) {
    local_mesh->free_vtx[i] += v[i];
  }

  /* free the variables */
  for (i = 0; i < num_LI; i++) {
    MY_FREE(N[i]);
  }
  MY_FREE(N);
  MY_FREE(R);

  /* update the function and gradient */
  SM_DEBUG_ACTION(
      2,
      { SMcomputeFunction(local_mesh,smooth_param,local_mesh->opt_info->test_function); SMfindActiveSet(local_mesh->opt_info->num_values, local_mesh->opt_info->test_function, smooth_param->active_eps, local_mesh->opt_info->test_active); SMprintActiveSet(local_mesh->opt_info->test_active, local_mesh->opt_info->test_function); });
}

double *
SMformVerticalRHS(int n, double *function, int ind[MAX_DIM], double max_value)
{
  int i;
  double *R;

  MY_MALLOCN(R, (double *), sizeof(double) * n, 1);
  for (i = 0; i < n; i++) {
    R[i] = max_value - function[ind[i]];
  }
  return (R);
}

void
SMstepToCusp(SMlocal_mesh *local_mesh)
{
  int ind0, ind1, i;
  double f0, f1, proj0, proj1, denom;
  double alpha_cusp, s_perp[MAX_DIM];
  SMoptimal *opt_info;
#ifdef LOCALTEST
  FILE *fp;
#endif

  opt_info = local_mesh->opt_info;
  s_perp[XDIR] = -1. * opt_info->search[YDIR];
  s_perp[YDIR] = opt_info->search[XDIR];
  s_perp[ZDIR] = 0.;

  ind0 = opt_info->active->active_ind[0];
  ind1 = opt_info->active->active_ind[1];

  f0 = opt_info->function[ind0];
  f1 = opt_info->function[ind1];
  SM_DOT(proj0, opt_info->gradient[ind0], s_perp, local_mesh->dimension);
  SM_DOT(proj1, opt_info->gradient[ind1], s_perp, local_mesh->dimension);
  denom = proj0 - proj1;
  if (!LESS_THAN_MACHINE_EPS(denom)) {
    alpha_cusp = (f1 - f0) / denom;
    for (i = 0; i < local_mesh->dimension; i++) {
      local_mesh->free_vtx[i] += alpha_cusp * s_perp[i];
    }
  }

  /*
   if (opt_info->function[ind0] == opt_info->active->true_active_value) {
   f_min      = opt_info->function[ind0];
   f_other    = opt_info->function[ind1];
   SM_DOT(proj_min,opt_info->gradient[ind0],s_perp,opt_info->dimension);
   SM_DOT(proj_other,opt_info->gradient[ind1],s_perp,opt_info->dimension);
   } else {
   f_min      = opt_info->function[ind1];
   f_other    = opt_info->function[ind0];
   SM_DOT(proj_min,opt_info->gradient[ind1],s_perp,opt_info->dimension);
   SM_DOT(proj_other,opt_info->gradient[ind0],s_perp,opt_info->dimension);
   }

   if ( fabs(proj_min - proj_other) > MACHINE_EPS) {
   alpha_cusp = (f_min - f_other)/(proj_other - proj_min);

   for (i=0;i<opt_info->dimension;i++) {
   local_mesh->free_vtx[i] += alpha_cusp*s_perp[i];
   }
   }
   */
  SM_DEBUG_ACTION(
      2,
      { SMcomputeFunction(local_mesh,smooth_param,local_mesh->opt_info->test_function); SMfindActiveSet(local_mesh->opt_info->num_values, local_mesh->opt_info->test_function, smooth_param->active_eps, local_mesh->opt_info->test_active); SMprintActiveSet(local_mesh->opt_info->test_active, local_mesh->opt_info->test_function); });

}

void
SMgetGradientProjections(SMoptimal *opt_info)
{
  int i;

  for (i = 0; i < opt_info->num_values; i++) {
    SM_DOT(opt_info->gs[i], opt_info->gradient[i], opt_info->search,
	   opt_info->dimension);
  }
}

void
SMcomputeAlpha(SMoptimal *opt_info)
{
  int i, j, steepest;
  int ind;
  int num_values;
  int *active;
  double *function;
  double steepest_function;
  double steepest_grad;
  double alpha_i;
  double alpha;

  num_values = opt_info->num_values;
assert(num_values > 0);
  			    MY_MALLOC(active, (int *), sizeof(int) * num_values, 1);

  alpha = A_BIG_POS_NUMBER;

  for (i = 0; i < opt_info->num_values; i++) {
    active[i] = 0;
  }
  for (j = 0; j < opt_info->active->num_active; j++) {
    ind = opt_info->active->active_ind[j];
    active[ind] = 1;
  }

  function = opt_info->function;
  steepest = opt_info->steepest;
  steepest_function = function[steepest];
  steepest_grad = opt_info->gs[steepest];
  for (i = 0; i < num_values; i++) {
    /* if it's not active */
    if (!active[i]) {
      alpha_i = function[i] - steepest_function;

      if (fabs(opt_info->gs[steepest] - opt_info->gs[i]) > 1E-13) {
	/* compute line intersection */
	alpha_i = alpha_i / (steepest_grad - opt_info->gs[i]);
      }
      else {
	/* the lines don't intersect - it's not under consideration*/
	alpha_i = 0;
      }
      if ((alpha_i > 0) && (fabs(alpha_i) < fabs(alpha))) {
	alpha = fabs(alpha_i);
      }
    }
  }
  opt_info->alpha = alpha;
  MY_FREE(active);
  SM_DEBUG_ACTION(3, { fprintf(stdout,"The initial step size: %f\n",alpha); });
}

