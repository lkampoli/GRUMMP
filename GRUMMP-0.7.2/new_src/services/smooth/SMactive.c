#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

#define SMinitActiveInfo(max_pts, active) \
{  \
  int i99; \
  active->num_active = 0; \
  active->num_equal = 0; \
  active->true_active_value = 0; \
  for (i99=0;i99<max_pts;i99++) { \
    active->active_ind[i99] = 0; \
  } \
}

void
SMfindActiveSet(int num_values, double *function, double active_eps,
		SMactive *active_info)
{
  int i, ind;
  int num_active, num_equal;
  double function_val;
  double active_value0;
  double true_active;
  double temp;

  active_info->num_active = 0;
  active_info->num_equal = 0;
  active_info->true_active_value = 0;
  for (i = 0; i < num_values; i++) {
    active_info->active_ind[i] = 0;
  }
  /*    SMinitActiveInfo(num_values, active_info); */

  /* the first function value is our initial active value */
  num_active = 1;
  num_equal = 0;
  active_info->active_ind[0] = 0;
  true_active = function[0];

  /* first sort out the active set... 
   all vals within active_eps of smallest val */

  for (i = 1; i < num_values; i++) {
    function_val = function[i];
    true_active = MIN(function_val, true_active);
    active_value0 = function[active_info->active_ind[0]];
    temp = fabs(function_val - active_value0);
    if (function_val < active_value0) {
      if (temp > active_eps) {
	num_active = 1;
	num_equal = 0;
	active_info->active_ind[0] = i;
      }
      else if (temp < active_eps) {
	num_active += 1;
	ind = num_active - 1;
	active_info->active_ind[ind] = i;
	if (fabs(function_val - active_value0) < MACHINE_EPS) {
	  num_equal += 1;
	}
      }
    }
    else {
      if (temp < active_eps) {
	num_active += 1;
	ind = num_active - 1;
	active_info->active_ind[ind] = i;
	if (fabs(function_val - active_value0) < MACHINE_EPS) {
	  num_equal += 1;
	}
      }
    }
  }
  active_info->true_active_value = true_active;
  active_info->num_active = num_active;
  active_info->num_equal = num_equal;
}

void
SMprintActiveSet(SMactive *active, double *function)
{
  int i;
  /* This next line is designed to eliminate compiler warnings. */
  double d = *function;
  d++;

  /* print the active set */
  for (i = 0; i < active->num_active; i++) {
    SM_DEBUG_ACTION(
	2,
	{ fprintf(stdout,"Active value %d:   %f \n", i+1,function[active->active_ind[i]]); });
  }
}
