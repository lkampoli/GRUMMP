#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SMsmooth.h"

SMlocal_mesh *
SMmallocLocalMesh(int num_pts, SMparam *smooth_param)
{
  int i, num_tri, num_values;
  SMlocal_mesh *local_mesh;

  MY_MALLOCN(local_mesh, (SMlocal_mesh *), sizeof(SMlocal_mesh), 1);

  num_tri = 2 * num_pts;
  num_values = num_tri * smooth_param->function_values_per_tri;
  local_mesh->max_num_pts = num_pts;
  local_mesh->max_num_tri = num_tri;
  local_mesh->num_values = num_values;

  /* malloc more space than you think you'll ever need, so you only need to do
   it once.  Based on 3D problem */
  MY_MALLOCN(local_mesh->incident_vtx, (double **), sizeof(double *) * num_pts,
	     1);
  MY_MALLOCN(local_mesh->vtx_connectivity, (int **), sizeof(int *) * num_tri, 1);
  for (i = 0; i < num_pts; i++) {
    MY_MALLOCN(local_mesh->incident_vtx[i], (double *), sizeof(double)*MAX_DIM,
	       1);
  }
  for (i = 0; i < num_tri; i++) {
    MY_MALLOCN(local_mesh->vtx_connectivity[i], (int *), sizeof(int)*MAX_DIM, 1);
  }
  MY_MALLOCN(local_mesh->original_pt, (double *), sizeof(double )*MAX_DIM, 1);
  MY_MALLOCN(local_mesh->original_function, (double *),
	     sizeof(double) * num_values, 1);

  local_mesh->lap_info = SMmallocLap(num_values);
  local_mesh->opt_info = SMmallocOpt(num_values);
  local_mesh->lp_info = SMmallocLP(5, num_tri);
  SMinitLap(num_values, local_mesh->lap_info);
  SMinitOpt(num_values, 20, local_mesh->opt_info);

  return local_mesh;
}

SMlap_info *
SMmallocLap(int num_values)
{
  SMlap_info *lap_info;
  MY_MALLOCN(lap_info, (SMlap_info *), sizeof(SMlap_info), 1);
  MY_MALLOCN(lap_info->laplacian_function, (double *),
	     sizeof(double) * num_values, 1);
  return lap_info;
}

SMoptimal *
SMmallocOpt(int num_values)
{
  int i;
  SMoptimal *opt_info;

  MY_MALLOCN(opt_info, (SMoptimal *), sizeof(SMoptimal), 1);

  MY_MALLOCN(opt_info->function, (double *), sizeof(double) * num_values, 1);
  MY_MALLOCN(opt_info->test_function, (double *), sizeof(double) * num_values,
	     1);
  MY_MALLOCN(opt_info->original_function, (double *),
	     sizeof(double) * num_values, 1);
  MY_MALLOCN(opt_info->gradient, (double **), sizeof(double *) * num_values, 1);
  MY_MALLOCN(opt_info->gs, (double *), sizeof(double) * num_values, 1);
  MY_MALLOCN(opt_info->G, (double **), sizeof(double *)*MAX_G_NUM, 1);
  MY_MALLOCN(opt_info->PDG, (double **), sizeof(double *)*MAX_DIM, 1);

  for (i = 0; i < MAX_DIM; i++) {
    MY_MALLOCN(opt_info->PDG[i], (double *), sizeof(double)*MAX_DIM, 1);
  }
  opt_info->num_grad_slots = num_values;
  for (i = 0; i < num_values; i++) {
    MY_MALLOCN(opt_info->gradient[i], (double *), sizeof(double)*MAX_DIM, 1);
  }
  for (i = 0; i < MAX_G_NUM; i++) {
    MY_MALLOCN(opt_info->G[i], (double *), sizeof(double)*MAX_G_NUM, 1);
  }

  MY_MALLOCN(opt_info->prev_active_values, (double *), sizeof(double) * 50, 1);

  /* malloc the active information sets */
  opt_info->active = SMmallocActive(num_values);
  opt_info->test_active = SMmallocActive(num_values);
  opt_info->original_active = SMmallocActive(num_values);

  return opt_info;
}

SMlp *
SMmallocLP(int num_active, int num_constraints)
{
  int i;
  SMlp *lp_info;

  MY_MALLOCN(lp_info, (SMlp *), sizeof(SMlp), 1);

  MY_MALLOCN(lp_info->ipivot, (int *), sizeof(int) * num_active, 1);
  MY_MALLOCN(lp_info->Amat, (double **), sizeof(double *) * num_constraints, 1);
  MY_MALLOCN(lp_info->Amat_T, (double **), sizeof(double *) * num_active, 1);
  MY_MALLOCN(lp_info->Amat_T_O, (double **), sizeof(double *) * num_active, 1);
  MY_MALLOCN(lp_info->feasible_x, (double *), sizeof(double) * num_constraints,
	     1);
  MY_MALLOCN(lp_info->free_ind, (int *), sizeof(int) * num_constraints, 1);
  MY_MALLOCN(lp_info->active_ind, (int *), sizeof(int) * num_active, 1);
  MY_MALLOCN(lp_info->b, (double *), sizeof(double) * num_constraints, 1);
  MY_MALLOCN(lp_info->c, (double *), sizeof(double) * num_constraints, 1);
  MY_MALLOCN(lp_info->pi, (double *), sizeof(double) * num_active, 1);
  MY_MALLOCN(lp_info->s, (double *), sizeof(double) * num_constraints, 1);
  MY_MALLOCN(lp_info->alpha, (double *), sizeof(double) * num_active, 1);
  MY_MALLOCN(lp_info->step, (double *), sizeof(double) * num_constraints, 1);
  MY_MALLOCN(lp_info->Bmat, (double *),
	     sizeof(double) * num_active * num_active, 1);
  MY_MALLOCN(lp_info->Bmat_T, (double *),
	     sizeof(double) * num_active * num_active, 1);
  MY_MALLOCN(lp_info->AAmat, (double **), sizeof(double *) * num_active, 1);

  for (i = 0; i < num_active; i++) {
    MY_MALLOCN(lp_info->AAmat[i], (double*), sizeof(double) * num_constraints,
	       1);
    MY_MALLOCN(lp_info->Amat_T_O[i], (double *),
	       sizeof(double) * num_constraints, 1);
    MY_MALLOCN(lp_info->Amat_T[i], (double *), sizeof(double) * num_constraints,
	       1);
  }
  for (i = 0; i < num_constraints; i++)
    MY_MALLOCN(lp_info->Amat[i], (double *), sizeof(double) * num_active, 1);

  return lp_info;
}

void
SMinitLocalMesh(int num_incident_vtx, int num_tri, double *free_vtx,
		double **vtx_list, int **vtx_connectivity,
		SMlocal_mesh *local_mesh, SMparam *smooth_param)
{
  int i, j;
  int num_values;
  int dimension;
  double min[3], max[3];
  double coord;

  local_mesh->num_incident_vtx = num_incident_vtx;
  local_mesh->num_tri = num_tri;
  dimension = local_mesh->dimension;

  num_values = num_tri * smooth_param->function_values_per_tri;
  local_mesh->num_values = num_values;

  SM_COPY_VECTOR(local_mesh->free_vtx, free_vtx, dimension);
  SM_COPY_VECTOR(local_mesh->original_pt, free_vtx, dimension);

  /* find the min and max coordinates in the local mesh */
  for (i = 0; i < dimension; i++) {
    min[i] = 1E300;
    max[i] = -1E300;
  }
  for (i = 0; i < num_incident_vtx; i++) {
    for (j = 0; j < dimension; j++) {
      coord = vtx_list[i][j];
      local_mesh->incident_vtx[i][j] = coord;
      if (coord < min[j])
	min[j] = coord;
      if (coord > max[j])
	max[j] = coord;
    }
  }
  if (local_mesh->free_vtx[XDIR] < min[XDIR])
    min[XDIR] = local_mesh->free_vtx[XDIR];
  if (local_mesh->free_vtx[YDIR] < min[YDIR])
    min[YDIR] = local_mesh->free_vtx[YDIR];
  if (local_mesh->free_vtx[XDIR] > max[XDIR])
    max[XDIR] = local_mesh->free_vtx[XDIR];
  if (local_mesh->free_vtx[YDIR] > max[YDIR])
    max[YDIR] = local_mesh->free_vtx[YDIR];

  if (vtx_connectivity != NULL) {
    for (i = 0; i < num_tri; i++)
      SM_COPY_VECTOR(local_mesh->vtx_connectivity[i], vtx_connectivity[i],
		     dimension);
  }
  else {
    fprintf(stderr, "A null connectivity vector was sent\n");
  }

  SM_COPY_VECTOR(local_mesh->min, min, dimension);
  SM_COPY_VECTOR(local_mesh->max, max, dimension);
  local_mesh->original_value = 0;
  local_mesh->current_active_value = 0;
  for (i = 0; i < local_mesh->num_values; i++) {
    local_mesh->original_function[i] = 0.0;
  }

  local_mesh->lap_done = FALSE;
  local_mesh->opt_done = FALSE;

  SM_DEBUG_ACTION(3, {PRINT_ORDERED_PTS(local_mesh); });
  SM_DEBUG_ACTION(3, {WRITE_BINARY_ORDERED_PTS(local_mesh);});

}

void
SMinitLap(int num_values, SMlap_info *lap_info)
{
  int i;

  lap_info->laplacian_value = 0;
  lap_info->lap_improvement = 0;
  lap_info->lap_invalid = FALSE;
  lap_info->lap_accepted = FALSE;
  for (i = 0; i < num_values; i++) {
    lap_info->laplacian_function[i] = 0.;
  }
}

void
SMinitLP(SMlocal_mesh *local_mesh)
{
  int i, j;
  int num_active, num_constraints;
  SMlp *lp_info;

  lp_info = local_mesh->lp_info;

  num_active = 5;
  num_constraints = local_mesh->num_tri;

  for (i = 0; i < num_constraints; i++) {
    lp_info->feasible_x[i] = 0;
    lp_info->b[i] = 0;
    lp_info->c[i] = 0;
    lp_info->free_ind[i] = 0;
    lp_info->s[i] = 0;
    lp_info->step[i] = 0;
  }
  for (i = 0; i < num_active; i++) {
    lp_info->active_ind[i] = 0;
    lp_info->pi[i] = 0;
    lp_info->alpha[i] = 0;
    lp_info->ipivot[i] = 0;
    for (j = 0; j < num_constraints; j++) {
      lp_info->AAmat[i][j] = 0;
      lp_info->Amat_T[i][j] = 0;
      lp_info->Amat_T_O[i][j] = 0;
      lp_info->Amat[j][i] = 0;
    }
  }
  for (i = 0; i < num_active * num_active; i++) {
    lp_info->Bmat[i] = 0;
    lp_info->Bmat_T[i] = 0;
  }
}

void
SMinitSmoothParam(int technique, int FunctionID, double Threshold,
		  void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  SMparam *smooth_param;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_param = smooth_data->smooth_param;

  smooth_param->iter_count = 0;
  smooth_param->new_init_pt_option = SM_NONE;
  smooth_param->maxit = 20;
  smooth_param->conv_eps = 1E-10;
  smooth_param->active_eps = .00003;
  smooth_param->min_acceptable_imp = .000001;
  smooth_param->min_step_size = 1E-8;
  smooth_param->function_id = FunctionID;
  smooth_param->global_min_value = INIT_MIN_VALUE;
  SMsetSmoothTechnique(ext_smooth_data, technique);
  SMsetSmoothFunction(ext_smooth_data, FunctionID);
  SMsetSmoothThreshold(ext_smooth_data, Threshold);
}

SMprocinfo *
SMinitProcinfo(void)
{
#ifdef PARALLEL_LOG
  int nprocs, myid;
#endif
  SMprocinfo *procinfo;

#ifdef PARALLEL_LOG
  /* set up procinfo */
  procinfo = (SMprocinfo *)malloc(sizeof(SMprocinfo));
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  procinfo->nprocs = nprocs;
  procinfo->myid = myid;
  procinfo->procset = MPI_COMM_WORLD;
#else
  procinfo = (SMprocinfo *) malloc(sizeof(SMprocinfo));
  procinfo->nprocs = 1;
  procinfo->myid = 0;
#endif

  return procinfo;
}

SMquality_table *
SMmallocQualityTable(void)
{
  int i;
  char function_name[10][128];
  double target[10];
  int dimension[10];
  SMquality_table *quality_table;

  MY_MALLOCN(quality_table, (SMquality_table *), sizeof(SMquality_table), 1);

  quality_table->num_functions = 9;
  /* 2D functions */
  strcpy(function_name[0], "Angle");
  target[0] = 60;
  dimension[0] = 2;
  strcpy(function_name[1], "Deviation from Equil");
  target[1] = 0;
  dimension[1] = 2;
  strcpy(function_name[2], "Scaled Jacobian");
  target[2] = .86;
  dimension[2] = 2;
  strcpy(function_name[3], "Triangle Area");
  target[3] = 0;
  dimension[3] = 2;

  /* 3D functions */
  strcpy(function_name[4], "Dihedral Angle");
  target[4] = 72;
  dimension[4] = 3;
  strcpy(function_name[5], "Scaled Jacobian");
  target[5] = .86;
  dimension[5] = 3;
  strcpy(function_name[6], "Edge Length/Volume");
  target[6] = .1;
  dimension[6] = 3;
  strcpy(function_name[7], "Tet Area");
  target[7] = 0;
  dimension[7] = 3;
  strcpy(function_name[8], "Condition Number");
  target[8] = 1;
  dimension[8] = 3;

  MY_MALLOC(quality_table->measure, (SMqualityMeasure **),
	    sizeof(SMqualityMeasure *) * quality_table->num_functions, 1);
  for (i = 0; i < quality_table->num_functions; i++) {
    MY_MALLOC(quality_table->measure[i], (SMqualityMeasure *),
	      sizeof(SMqualityMeasure), 1);
    strcpy(quality_table->measure[i]->name, function_name[i]);
    quality_table->measure[i]->dimension = dimension[i];
    quality_table->measure[i]->target = target[i];
    quality_table->measure[i]->min_value = A_BIG_POS_NUMBER;
    quality_table->measure[i]->max_value = A_BIG_NEG_NUMBER;
    quality_table->measure[i]->avg_value = 0.0;
    quality_table->measure[i]->num_function_values = 0;
  }
  quality_table->num_tangled_elements = 0;
  quality_table->initialized = FALSE;
  return quality_table;
}

void
SMinitOpt(int num_values, int maxit, SMoptimal *opt_info)
{
  int i, j;

  /* for the purposes of initialization will be set to zero after */
  opt_info->num_values = num_values;
  opt_info->equilibrium_pt = 0;
  opt_info->step_too_small = 0;
  opt_info->step_accepted = 0;
  opt_info->status = 0;
  opt_info->iter_count = 0;
  opt_info->opt_iter_count = 0;
  opt_info->opt_improvement = 0;
  opt_info->maxit = maxit;
  opt_info->steepest = 0;
  opt_info->alpha = 0;

  for (i = 0; i < MAX_DIM; i++) {
    opt_info->search[i] = 0;
    opt_info->PDG_ind[i] = -1;
    for (j = 0; j < MAX_DIM; j++)
      opt_info->PDG[i][j] = 0;
  }
  for (i = 0; i < num_values; i++) {
    opt_info->function[i] = 0;
    opt_info->test_function[i] = 0;
    opt_info->original_function[i] = 0;
    opt_info->gs[i] = 0;
    for (j = 0; j < MAX_DIM; j++) {
      opt_info->gradient[i][j] = 0;
    }
  }
  if (num_values > MAX_G_NUM) {
    for (i = 0; i < MAX_G_NUM; i++) {
      for (j = 0; j < MAX_G_NUM; j++)
	opt_info->G[i][j] = 0;
    }
  }
  else {
    for (i = 0; i < num_values; i++) {
      for (j = 0; j < num_values; j++)
	opt_info->G[i][j] = 0;
    }
  }

  for (i = 0; i < 20; i++)
    opt_info->prev_active_values[i] = 0;
}

SMstats *
SMinitStats(SMstats *smooth_stats)
{
  smooth_stats->stats_initialized = 1;
  smooth_stats->total_verts_smoothed = 0;
  smooth_stats->num_equil = 0;
  smooth_stats->num_started_equil = 0;
  smooth_stats->num_zero_search = 0;
  smooth_stats->num_imp_too_small = 0;
  smooth_stats->num_flat_no_imp = 0;
  smooth_stats->num_step_too_small = 0;
  smooth_stats->num_max_iter_exceeded = 0;
  smooth_stats->num_surface_verts_opted = 0;
  smooth_stats->num_verts_opted = 0;
  smooth_stats->opt_count = 0.0;
  smooth_stats->opt_iter_count = 0;
  smooth_stats->num_verts_laplaced = 0;
  smooth_stats->num_lap_enough = 0;
  smooth_stats->num_lap_invalid = 0;
  smooth_stats->num_lap_worse = 0;
  smooth_stats->no_improvement = 0;
  smooth_stats->avg_improvement = 0.0;
  smooth_stats->avg_active_val = 0.0;
  smooth_stats->global_minimum_val = INIT_MIN_VALUE;

  return smooth_stats;
}

SMactive *
SMmallocActive(int num_values)
{
  SMactive *active;

  MY_MALLOCN(active, (SMactive *), sizeof(SMactive), 1);
  MY_MALLOCN(active->active_ind, (int *), sizeof(int) * num_values, 1);

  return active;
}

void
SMconvertToDegrees(int function_id, double *value)
{
  static const double dInvRad = 180 / M_PI;
  switch (function_id)
    {
    case MAX_MIN_ANGLE:
    case MIN_MAX_ANGLE:
    case MAX_MIN_INTERIOR_ANGLE:
    case MAX_MIN_DIHEDRAL:
    case MIN_MAX_DIHEDRAL:
    case MAX_MIN_SOLID:
      /* value is currently in radians */
      *value = *value * dInvRad;
      break;
    case MIN_MAX_COSINE:
    case MAX_MIN_COSINE:
    case MIN_MAX_INTERIOR_COSINE:
    case MIN_MAX_COSINE_DIHEDRAL:
    case MAX_MIN_COSINE_DIHEDRAL:
      /* value is a cosine */
      *value = SAFE_ACOS(*value) * dInvRad;
      break;
    case MAX_MIN_SINE:
    case MAX_MIN_INTERIOR_SINE:
    case MAX_MIN_SINE_DIHEDRAL:
      /* value is a sine */
      *value = SAFE_ASIN(*value) * dInvRad;
      break;
    default:
      break; /* For other cases, the conversion makes no sense. */
    }
  SM_DEBUG_ACTION(
      3, {fprintf(stdout,"Converting the min value to degrees %f\n",*value);});
}
