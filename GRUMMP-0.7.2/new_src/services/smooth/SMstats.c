#include <stdio.h>
#include <math.h>
#include "GR_misc.h"
#include "SMsmooth.h"

void
SMaccumulateStats(SMlocal_mesh *local_mesh, SMparam *smooth_param,
		  SMstats *smooth_stats)
{
  double imp;

  /* this vert was smoothed */
  smooth_stats->total_verts_smoothed += 1;

  /* Laplacian statistics */
  if (local_mesh->lap_done) {
    smooth_stats->num_verts_laplaced += 1;
    if (local_mesh->lap_info->lap_accepted == FALSE) {
      if (local_mesh->lap_info->lap_invalid == TRUE) {
	smooth_stats->num_lap_invalid += 1;
      }
      else {
	smooth_stats->num_lap_worse += 1;
      }
    }
  }

  /* if combined - how many was laplacian enough for */
  if ((smooth_param->smooth_technique != OPTIMIZATION_ONLY)
      && (smooth_param->smooth_technique != LAPLACE_ONLY)
      && (local_mesh->lap_done) && (local_mesh->lap_info->lap_accepted)
      && (local_mesh->lap_info->laplacian_value > smooth_param->lap_accept_value)) {
    smooth_stats->num_lap_enough += 1;
  }

  /* if optimization used */
  if (local_mesh->opt_done) {
    smooth_stats->num_verts_opted += 1;
    if (smooth_param->is_surface_vertex)
      smooth_stats->num_surface_verts_opted += 1;
    smooth_stats->opt_iter_count += local_mesh->opt_info->opt_iter_count;

    /* optimization termination status */
    switch (local_mesh->opt_info->status)
      {
      case EQUILIBRIUM:
	smooth_stats->num_equil += 1;
	if (local_mesh->opt_info->iter_count == 0)
	  smooth_stats->num_started_equil += 1;
	break;
      case ZERO_SEARCH:
	smooth_stats->num_zero_search += 1;
	break;
      case IMP_TOO_SMALL:
	smooth_stats->num_imp_too_small += 1;
	break;
      case FLAT_NO_IMP:
	smooth_stats->num_flat_no_imp += 1;
	break;
      case STEP_TOO_SMALL:
	smooth_stats->num_step_too_small += 1;
	break;
      case MAX_ITER_EXCEEDED:
	smooth_stats->num_max_iter_exceeded += 1;
	break;
      default:
	fprintf(stderr, "An invalid termination status\n");
	WRITE_BINARY_ORDERED_PTS(local_mesh)
	;
	break;
      }
  }

  /* active value in local submesh */
  smooth_stats->avg_active_val += local_mesh->current_active_value;

  /* global active value in the mesh */
  if (smooth_stats->global_minimum_val > local_mesh->current_active_value) {
    smooth_stats->global_minimum_val = local_mesh->current_active_value;
  }

  /* improvement in local submesh */
  imp = local_mesh->current_active_value - local_mesh->original_value;
  smooth_stats->avg_improvement += imp;

  /* was there actually any improvement */
  if (imp == 0) {
    smooth_stats->no_improvement += 1;
  }
}

void
SMprintStats(SMsmooth_data *smooth_data)
{
  SMstats *smooth_stats;
  int num_verts, num_verts_laplaced, num_verts_opted;
  double min_value;

  smooth_stats = smooth_data->smooth_stats;

  logMessage(
      2,
      "******************************************************************\n");
  num_verts = smooth_stats->total_verts_smoothed;
  num_verts_laplaced = smooth_stats->num_verts_laplaced;
  num_verts_opted = smooth_stats->num_verts_opted;
  min_value = smooth_stats->global_minimum_val;
  SMconvertToDegrees(smooth_data->smooth_param->function_id, &min_value);

  logMessage(1, "Approximate global minimum value (degrees)  %11.5f \n",
	     min_value);
  logMessage(1, "# of vertices smoothed                      %11d \n",
	     num_verts);
  logMessage(2, "# of calls to Laplacian smoothing           %11d \n",
	     num_verts_laplaced);
  if (num_verts_laplaced != 0) {
    logMessage(2, "     Invalid Laplacian steps         %6d (%4.1f%%) \n",
	       smooth_stats->num_lap_invalid,
	       (100. * smooth_stats->num_lap_invalid) / num_verts_laplaced);
    logMessage(2, "     Non-improvement Laplacian steps %6d (%4.1f%%) \n",
	       smooth_stats->num_lap_worse,
	       (100. * smooth_stats->num_lap_worse) / num_verts_laplaced);
  }
  logMessage(2, "# of calls to optimization smoothing        %11d \n",
	     num_verts_opted);
  logMessage(2, "# of calls to surface smoothing             %11d \n",
	     smooth_stats->num_surface_verts_opted);
  if (num_verts_opted != 0) {
    logMessage(3, "     Average iterations/optimization call %4.2f \n",
	       smooth_stats->opt_iter_count / ((double) num_verts_opted));
  }
  logMessage(2, "# of vertices with no improvement           %11d \n",
	     smooth_stats->no_improvement);
  if (num_verts > 0) {
    logMessage(2, "Average final active value                  %11.3f \n",
	       (smooth_stats->avg_active_val / (double) num_verts));
    logMessage(3, "Average improvement (all smoothed vertices) %11.5f \n",
	       (smooth_stats->avg_improvement / (double) num_verts));
  }
  logMessage(3, "Termination status: \n");
  logMessage(3, "  Laplacian smoothing enough      %6d (%4.1f%%) \n",
	     smooth_stats->num_lap_enough,
	     (100. * smooth_stats->num_lap_enough) / num_verts);
  logMessage(3, "  Equilibrium pt found            %6d (%4.1f%%) \n",
	     smooth_stats->num_equil,
	     (100. * smooth_stats->num_equil) / num_verts);
  logMessage(3, "     Started at equilibrium       %6d (%4.1f%%) \n",
	     smooth_stats->num_started_equil,
	     (100. * smooth_stats->num_started_equil) / num_verts);
  logMessage(3, "  Zero search direction           %6d (%4.1f%%) \n",
	     smooth_stats->num_zero_search,
	     (100. * smooth_stats->num_zero_search) / num_verts);
  logMessage(3, "  Improvement too small           %6d (%4.1f%%) \n",
	     smooth_stats->num_imp_too_small,
	     (100. * smooth_stats->num_imp_too_small) / num_verts);
  logMessage(3, "  Flat no improvement             %6d (%4.1f%%) \n",
	     smooth_stats->num_flat_no_imp,
	     (100. * smooth_stats->num_flat_no_imp) / num_verts);
  logMessage(3, "  Step too small                  %6d (%4.1f%%) \n",
	     smooth_stats->num_step_too_small,
	     (100. * smooth_stats->num_step_too_small) / num_verts);
  logMessage(3, "  Max Iter Exceeded               %6d (%4.1f%%) \n",
	     smooth_stats->num_max_iter_exceeded,
	     (100. * smooth_stats->num_max_iter_exceeded) / num_verts);
  logMessage(
      2,
      "******************************************************************\n");
}

void
SMwriteStatsFile(SMstats *smooth_stats, int smooth_count)
{
  FILE *fp;
  int num_verts;

  num_verts = smooth_stats->total_verts_smoothed;

  fp = fopen("Smooth.stats", "a");
  fprintf(fp, "\n");
  fprintf(
      fp,
      "******************************************************************\n");
  fprintf(fp,
	  "                            ITERATION  %d                        \n",
	  smooth_count);
  fprintf(
      fp,
      "******************************************************************\n");
  fprintf(fp, "  The number of verts smoothed                            %d \n",
	  num_verts);
  fprintf(fp, "  The number of equilibrium pts found                     %d \n",
	  smooth_stats->num_equil);
  fprintf(fp, "  The number of verts that started with an equil point    %d \n",
	  smooth_stats->num_started_equil);
  if (smooth_stats->num_verts_opted == 0) {
    fprintf(fp, "  The average optimization count is 0\n");
  }
  else {
    fprintf(fp,
	    "  The average optimization iteration count                %f \n",
	    smooth_stats->opt_count / ((double) smooth_stats->num_verts_opted));
  }
  fprintf(fp, "  The number of verts with no improvement                 %d \n",
	  smooth_stats->no_improvement);
  fprintf(fp, "  The average active value                                %f \n",
	  (smooth_stats->avg_active_val / (double) num_verts));
  fprintf(fp, "  The average improvement                                 %f \n",
	  (smooth_stats->avg_improvement / (double) num_verts));
  fprintf(
      fp,
      "******************************************************************\n");
  fclose(fp);
}
