#include <stdlib.h>
#include <stdio.h>
#include "SMsmooth.h"
/*@
 SMfinalizeSmoothing - This routine frees all the memory allocated in SMinitSmoothing
 including the data structure smooth_data.  This routine should be called when mesh 
 optimization is complete.

 Input Parameters:
 .  smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be called
 prior to calling this routine.

 .seealso SMinitSmoothing()
 @*/
void
SMfinalizeSmoothing(void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;

  MY_FREE(smooth_data->smooth_stats);
  SMfreeLocalMesh(smooth_data->local_mesh);
  SMfreeParam(smooth_data->smooth_param);
  SMfreeProcinfo(smooth_data->smooth_procinfo);
  SMfreeQualityTable(smooth_data->quality_table);
  MY_FREE(smooth_data);
}

void
SMfreeLocalMesh(SMlocal_mesh *local_mesh)
{
  int i;

  for (i = 0; i < local_mesh->max_num_pts; i++) {
    MY_FREE(local_mesh->incident_vtx[i]);
  }
  for (i = 0; i < local_mesh->max_num_tri; i++) {
    MY_FREE(local_mesh->vtx_connectivity[i]);
  }
  MY_FREE(local_mesh->original_pt);
  MY_FREE(local_mesh->incident_vtx);
  MY_FREE(local_mesh->vtx_connectivity);
  MY_FREE(local_mesh->original_function);
  MY_FREE(local_mesh->lap_info->laplacian_function);
  MY_FREE(local_mesh->lap_info);
  SMfreeOptimal(local_mesh);
  SMfreeLP(local_mesh, 5, local_mesh->max_num_tri);
  MY_FREE(local_mesh);
}

void
SMfreeOptimal(SMlocal_mesh *local_mesh)
{
  int i;
  int num_values;
  SMoptimal *opt_info;

  opt_info = local_mesh->opt_info;

  num_values = opt_info->num_grad_slots;

  for (i = 0; i < num_values; i++)
    MY_FREE(opt_info->gradient[i]);
  for (i = 0; i < MAX_G_NUM; i++)
    MY_FREE(opt_info->G[i]);

  for (i = 0; i < MAX_DIM; i++) {
    MY_FREE(opt_info->PDG[i]);
  }

  MY_FREE(opt_info->function);
  MY_FREE(opt_info->test_function);
  MY_FREE(opt_info->original_function);
  MY_FREE(opt_info->gradient);
  MY_FREE(opt_info->gs);
  MY_FREE(opt_info->G);
  MY_FREE(opt_info->PDG);

  /* free the active information sets */
  SMfreeActive(opt_info->active);
  SMfreeActive(opt_info->test_active);
  SMfreeActive(opt_info->original_active);

  MY_FREE(opt_info->prev_active_values);

  MY_FREE(opt_info);
}

void
SMfreeLP(SMlocal_mesh *local_mesh, int num_active, int num_constraints)
{
  int i;
  SMlp *lp_info;

  lp_info = local_mesh->lp_info;
  for (i = 0; i < num_active; i++) {
    MY_FREE(lp_info->AAmat[i]);
    MY_FREE(lp_info->Amat_T[i]);
    MY_FREE(lp_info->Amat_T_O[i]);
  }
  for (i = 0; i < num_constraints; i++)
    MY_FREE(lp_info->Amat[i]);

  MY_FREE(lp_info->ipivot);
  MY_FREE(lp_info->Amat);
  MY_FREE(lp_info->Amat_T);
  MY_FREE(lp_info->Amat_T_O);
  MY_FREE(lp_info->feasible_x);
  MY_FREE(lp_info->free_ind);
  MY_FREE(lp_info->active_ind);
  MY_FREE(lp_info->b);
  MY_FREE(lp_info->c);
  MY_FREE(lp_info->pi);
  MY_FREE(lp_info->s);
  MY_FREE(lp_info->alpha);
  MY_FREE(lp_info->step);
  MY_FREE(lp_info->Bmat);
  MY_FREE(lp_info->Bmat_T);
  MY_FREE(lp_info->AAmat);

  MY_FREE(lp_info);
}

void
SMfreeActive(SMactive *active)
{
  MY_FREE(active->active_ind);
  MY_FREE(active);
}

void
SMfreeParam(SMparam *smooth_param)
{
  MY_FREE(smooth_param);
}

void
SMfreeProcinfo(SMprocinfo *procinfo)
{
  MY_FREE(procinfo);
}

void
SMfreeQualityTable(SMquality_table *quality_table)
{
  int i;
  for (i = 0; i < quality_table->num_functions; i++) {
    MY_FREE(quality_table->measure[i]);
  }
  MY_FREE(quality_table->measure);
  MY_FREE(quality_table);
}

void
SMfreeActiveDirections(const int num_active, double **dir)
{
  int i;
  for (i = 0; i < num_active; i++) {
    MY_FREE(dir[i]);
  }
  MY_FREE(dir);
}

