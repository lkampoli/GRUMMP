#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

void
SMwriteLocalMesh(FILE *fp, SMlocal_mesh *local_mesh)
{
  SMwriteLocalTriangleList(fp, local_mesh);
}

void
SMwriteLocalAxes(FILE *fp, SMlocal_mesh *mesh)
{
  if (mesh->dimension == 2) {
    if (mesh->free_vtx[XDIR] < mesh->min[XDIR])
      mesh->min[XDIR] = mesh->free_vtx[XDIR];
    if (mesh->free_vtx[YDIR] < mesh->min[YDIR])
      mesh->min[YDIR] = mesh->free_vtx[YDIR];
    if (mesh->free_vtx[XDIR] > mesh->max[XDIR])
      mesh->max[XDIR] = mesh->free_vtx[XDIR];
    if (mesh->free_vtx[YDIR] > mesh->max[YDIR])
      mesh->max[YDIR] = mesh->free_vtx[YDIR];

    fprintf(
	fp,
	"axes('Xlim',[%16.10f,%16.10f],'Ylim',[%16.10f,%16.10f],'Box','on ',",
	mesh->min[XDIR], mesh->max[XDIR], mesh->min[YDIR], mesh->max[YDIR]);
    fprintf(fp, "'XTick',[%16.10f,%16.10f],'YTick',[%16.10f,%16.10f]);\n",
	    mesh->min[XDIR], mesh->max[XDIR], mesh->min[YDIR], mesh->max[YDIR]);
    fprintf(fp, "hold on;\n");
  }
  else if (mesh->dimension == 3) {
    printf("Can't look at three dimensions in matlab\n");
  }
}

void
SMwriteLocalTriangleList(FILE *fp, SMlocal_mesh *local_mesh)
{
  int i, num_tri;
  double *vtx1, *vtx2, *vtx3;
  double x_0, x_1, x_2;
  double y_0, y_1, y_2;

  SMwriteLocalAxes(fp, local_mesh);

  num_tri = local_mesh->num_tri;
  for (i = 0; i < num_tri; i++) {
    vtx1 = local_mesh->free_vtx;
    vtx2 = local_mesh->incident_vtx[local_mesh->vtx_connectivity[i][0]];
    vtx3 = local_mesh->incident_vtx[local_mesh->vtx_connectivity[i][1]];

    x_0 = vtx1[XDIR];
    y_0 = vtx1[YDIR];
    x_1 = vtx2[XDIR];
    y_1 = vtx2[YDIR];
    x_2 = vtx3[XDIR];
    y_2 = vtx3[YDIR];

    fprintf(fp, "line([%f, %f],[%f, %f]);\n", x_0, x_1, y_0, y_1);
    fprintf(fp, "line([%f, %f],[%f, %f]);\n", x_1, x_2, y_1, y_2);
    fprintf(fp, "line([%f, %f],[%f, %f]);\n", x_0, x_2, y_0, y_2);
  }
}

void
SMwriteActiveSet(FILE *fp, SMlocal_mesh *local_mesh)
{
  int i, ind;
  int tri_num;
  int vert_id, ind1, ind2;
  double x_pt, y_pt, x_1, x_2, y_1, y_2;
  double rise, run, slope, b_int;
  double midx, midy;
  double xdiff, x_temp, y_temp;

  /*PRINT_FUNCTION_VALUES(local_mesh->opt_info); */
  /* marks the minimum angle set with cyan * */
  if (local_mesh->opt_info->num_values == local_mesh->num_tri * 3) {
    ind1 = ind2 = 0;
    for (i = 0; i < local_mesh->opt_info->active->num_active; i++) {
      ind = local_mesh->opt_info->active->active_ind[i];

      tri_num = (int) ind / 3;
      vert_id = ind % 3;
      ind1 = local_mesh->vtx_connectivity[tri_num][0];
      ind2 = local_mesh->vtx_connectivity[tri_num][1];
      if (vert_id == 0) {
	x_pt = local_mesh->free_vtx[XDIR];
	y_pt = local_mesh->free_vtx[YDIR];
	x_1 = local_mesh->incident_vtx[ind1][XDIR];
	x_2 = local_mesh->incident_vtx[ind2][XDIR];
	y_1 = local_mesh->incident_vtx[ind1][YDIR];
	y_2 = local_mesh->incident_vtx[ind2][YDIR];
      }
      else if (vert_id == 1) {
	x_pt = local_mesh->incident_vtx[ind1][XDIR];
	y_pt = local_mesh->incident_vtx[ind1][YDIR];
	x_1 = local_mesh->free_vtx[XDIR];
	y_1 = local_mesh->free_vtx[YDIR];
	x_2 = local_mesh->incident_vtx[ind2][XDIR];
	y_2 = local_mesh->incident_vtx[ind2][YDIR];
      }
      else {
	x_pt = local_mesh->incident_vtx[ind2][XDIR];
	y_pt = local_mesh->incident_vtx[ind2][YDIR];
	x_1 = local_mesh->free_vtx[XDIR];
	y_1 = local_mesh->free_vtx[YDIR];
	x_2 = local_mesh->incident_vtx[ind1][XDIR];
	y_2 = local_mesh->incident_vtx[ind1][YDIR];
      }

      midx = (x_1 + x_2) / 2;
      midy = (y_1 + y_2) / 2;
      xdiff = x_pt - midx;

      rise = y_pt - midy;
      run = x_pt - midx;

      slope = rise / run;
      b_int = midy - slope * midx;
      x_temp = x_pt - .1 * xdiff;
      y_temp = slope * x_temp + b_int;

      fprintf(fp, "plot( %f , %f, 'c*'); \n", x_temp, y_temp);
    }
  }
  if (local_mesh->opt_info->num_values == local_mesh->num_tri) {
    ind1 = ind2 = 0;
    for (i = 0; i < local_mesh->opt_info->active->num_active; i++) {
      ind = local_mesh->opt_info->active->active_ind[i];

      tri_num = (int) ind;
      ind1 = local_mesh->vtx_connectivity[tri_num][0];
      ind2 = local_mesh->vtx_connectivity[tri_num][1];
      x_pt = local_mesh->free_vtx[XDIR];
      y_pt = local_mesh->free_vtx[YDIR];
      x_1 = local_mesh->incident_vtx[ind1][XDIR];
      x_2 = local_mesh->incident_vtx[ind2][XDIR];
      y_1 = local_mesh->incident_vtx[ind1][YDIR];
      y_2 = local_mesh->incident_vtx[ind2][YDIR];

      fprintf(fp, "line([%f, %f],[%f, %f],'Color','c');\n", x_1, x_2, y_1, y_2);
      fprintf(fp, "line([%f, %f],[%f, %f],'Color','c');\n", x_pt, x_2, y_pt,
	      y_2);
      fprintf(fp, "line([%f, %f],[%f, %f],'Color','c');\n", x_1, x_pt, y_1,
	      y_pt);
    }
  }
}

void
SMwriteSearch(FILE *fp, SMlocal_mesh *local_mesh)
{
  double x_1, y_1, x_2, y_2;

  x_1 = local_mesh->free_vtx[XDIR];
  y_1 = local_mesh->free_vtx[YDIR];
  x_2 = x_1 + local_mesh->opt_info->search[XDIR];
  y_2 = y_1 + local_mesh->opt_info->search[YDIR];

  fprintf(fp, "line([%f, %f],[%f, %f],'Color','m');\n", x_1, x_2, y_1, y_2);

}

void
SMwritePoint(FILE *fp, double x, double y)
{
  fprintf(fp, "plot(%f,%f,'r*');\n", x, y);
}

