#include "SM_config.h"
#include "SMdefs.h"
/* EXTERNAL FUNCTIONS */
/* functions needed in the user code */
#ifdef __cplusplus
extern "C"
{
#endif
  void
  SMinitSmoothing(char technique, int FunctionID, double AcceptFunction,
		  void **smooth_data);

  void
  SMsmooth(int num_pts, int num_tet, double *free_vtx,
	   double (*vtx_list)[MAX_DIM], int (*vtx_connectivity)[MAX_DIM],
	   void *smooth_data, int is_surface_vert);
  void
  SMsetSmoothTechnique(void *ext_smooth_data, char technique);
  void
  SMinitGlobalMinValue(void *ext_smooth_data);
  void
  SMsetSmoothThreshold(void *smooth_data, double accept);
  void
  SMsetSmoothFunction(void *smooth_data, int FunctionID);

  void
  SMinitSmoothStats(void *smooth_data);
  void
  SMprintSmoothStats(void *smooth_data);

  void
  SMfinalizeSmoothing(void *smooth_data);

#ifdef __cplusplus
}
#endif

