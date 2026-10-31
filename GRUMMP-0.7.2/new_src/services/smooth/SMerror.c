#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

int
SMerror(int line, char *func, char* file, char *dir, int n, int p, char *mess)
{
  int ierr = 0 * p;

  printf("OPT-MS ERROR:  %s()  line %d  in %s%s\n", func, line, dir, file);
  switch (n)
    {
    case OPTMS_MEM_ERR:
      printf("OPT-MS ERROR:    Out of memory. \n");
      break;
      /*   default: */
    }

  if (mess) {
    printf("OPT-MS ERROR:   %s\n", mess);
  }
  return (ierr);
}

int
SMwrite_ordered_points(SMlocal_mesh *local_mesh)
{
  FILE *fp;
  int i99, j99;

  if ((fp = fopen("debug.ascii", "w")) == NULL) {
    fprintf(stderr, "Can't open debug.ascii for writing\n");
    return (1);
  }

  fprintf(fp, "%d  %d\n", local_mesh->num_incident_vtx, local_mesh->num_tri);
  for (i99 = 0; i99 < local_mesh->dimension; i99++)
    fprintf(fp, "%f  ", local_mesh->free_vtx[i99]);
  fprintf(fp, "\n");
  for (i99 = 0; i99 < local_mesh->num_incident_vtx; i99++) {
    for (j99 = 0; j99 < local_mesh->dimension; j99++)
      fprintf(fp, "%f  ", local_mesh->incident_vtx[i99][j99]);
    fprintf(fp, "\n");
  }
  for (i99 = 0; i99 < local_mesh->num_tri; i99++) {
    for (j99 = 0; j99 < local_mesh->dimension; j99++)
      fprintf(fp, "%d  ", local_mesh->vtx_connectivity[i99][j99]);
    fprintf(fp, "\n");
  }
  return (0);
}

