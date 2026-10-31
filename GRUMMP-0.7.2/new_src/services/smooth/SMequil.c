#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

#define NO_EQUIL 101
#define CHECK_TOP_DOWN 102
#define CHECK_BOTTOM_UP 103
#define TWO_PT_PLANE 104
#define THREE_PT_PLANE 105
#define CHECK_Y_COORD_DIRECTION 106
#define CHECK_X_COORD_DIRECTION 107
#define CHECK_Z_COORD_DIRECTION 108
#define EQUIL 109
#define HULL_TEST_ERROR 110

int
SMconvexHullTest(double **vec, int num_vec)
{
  int status, equil, dir_done;
  double pt1[3], pt2[3], pt3[3];
  double normal[3];

  /* tries to determine equilibrium for the 3D case */
  equil = FALSE;
  status = CHECK_Z_COORD_DIRECTION;
  dir_done = -1;

  if (num_vec <= 2)
    status = NO_EQUIL;

  while ((status != EQUIL) && (status != NO_EQUIL)
      && (status != HULL_TEST_ERROR)) {
    if (status == CHECK_Z_COORD_DIRECTION) {
      SMfindPlanePoints(ZDIR, YDIR, vec, num_vec, pt1, pt2, pt3, &status);
      dir_done = 2;
    }
    if (status == CHECK_Y_COORD_DIRECTION) {
      SMfindPlanePoints(YDIR, XDIR, vec, num_vec, pt1, pt2, pt3, &status);
      dir_done = 1;
    }
    if (status == CHECK_X_COORD_DIRECTION) {
      SMfindPlanePoints(XDIR, ZDIR, vec, num_vec, pt1, pt2, pt3, &status);
      dir_done = 0;
    }
    if (status == TWO_PT_PLANE) {
      pt3[0] = 0.;
      pt3[1] = 0.;
      pt3[2] = 0.;
    }
    if ((status == TWO_PT_PLANE) || (status == THREE_PT_PLANE)) {
      SMfindPlaneNormal(pt1, pt2, pt3, normal);
      equil = SMcheckVectorDots(vec, num_vec, normal);
      if (equil == 1) {
	switch (dir_done)
	  {
	  case 2:
	    equil = 0;
	    status = CHECK_Y_COORD_DIRECTION;
	    break;
	  case 1:
	    equil = 0;
	    status = CHECK_X_COORD_DIRECTION;
	    break;
	  case 0:
	    equil = 1;
	    status = EQUIL;
	    break;
	  }
      }
      else if (equil == 0) {
	status = NO_EQUIL;
      }
      else {
	SM_DEBUG_PRINT(2, "Equilibrium point not equal 0 or 1\n");
	status = HULL_TEST_ERROR;
      }
    }
  }
  switch (status)
    {
    case NO_EQUIL:
      SM_DEBUG_PRINT(3, "Not an equilibrium point\n")
      ;
      equil = 0;
      break;
    case EQUIL:
      SM_DEBUG_PRINT(3, "An equilibrium point\n")
      ;
      equil = 1;
      break;
    default:
      SM_DEBUG_ACTION(
	  3,
	  { fprintf(stdout,"Failed to determine equil or not; status = %d\n",status); })
      ;
      break;
    }
  return equil;
}

int
SMcheckVectorDots(double **vec, int num_vec, double *normal)
{
  int i, ind, equil;
  double test_dot, dot;

  equil = FALSE;
  SM_DOT(test_dot, vec[0], normal, 3);
  ind = 1;
  while ((fabs(test_dot) < MACHINE_EPS) && (ind < num_vec)) {
    SM_DOT(test_dot, vec[ind], normal, 3);
    ind++;
  }

  for (i = ind; i < num_vec; i++) {
    SM_DOT(dot, vec[i], normal, 3);
    if (((dot > 0 && test_dot < 0) || (dot < 0 && test_dot > 0))
	&& (fabs(dot) > MACHINE_EPS)) {
      equil = TRUE;
      return equil;
    }
  }
  return equil;
}

void
SMfindPlaneNormal(double pt1[3], double pt2[3], double pt3[3], double *cross)
{
  int i;
  double vecA[3], vecB[3];

  for (i = 0; i < 3; i++) {
    vecA[i] = pt2[i] - pt1[i];
    vecB[i] = pt3[i] - pt1[i];
  }
  cross[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
  cross[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
  cross[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
  SM_NORMALIZE(cross, 3);
}

void
SMfindPlanePoints(int dir1, int dir2, double **vec, int num_vec, double *pt1,
		  double *pt2, double *pt3, int *status)
{
  int i;
  int ind[50], num_min, num_max;
  int rotate = CW;
  int num_rotated = 0;
  double pt_1, pt_2;
  double min, inv_slope;
  double min_inv_slope = 0.;
  double max;
  double max_inv_slope = 0;
  double inv_origin_slope = 0;
  double *junk = pt3;
  junk++;

  *status = CHECK_BOTTOM_UP;
  /* find the minimum points in dir1 starting at -1 */
  num_min = 0;
  ind[0] = -1;
  ind[1] = -1;
  ind[2] = -1;
  min = 1.0;
  for (i = 0; i < num_vec; i++) {
    if (vec[i][dir1] < min) {
      min = vec[i][dir1];
      ind[0] = i;
      num_min = 1;
    }
    else if (fabs(vec[i][dir1] - min) < MACHINE_EPS) {
      ind[num_min++] = i;
    }
  }
  if (min >= 0)
    *status = NO_EQUIL;

  if (*status != NO_EQUIL) {
    switch (num_min)
      {
      case 1: /* rotate to find the next point */
	SM_COPY_VECTOR3(pt1, vec[ind[0]])
	;
	pt_1 = pt1[dir1];
	pt_2 = pt1[dir2];
	if (pt1[dir2] <= 0) {
	  rotate = CCW;
	  max_inv_slope = A_BIG_NEG_NUMBER;
	}
	if (pt1[dir2] > 0) {
	  rotate = CW;
	  min_inv_slope = A_BIG_POS_NUMBER;
	}
	switch (rotate)
	  {
	  case CCW:
	    for (i = 0; i < num_vec; i++) {
	      if (i != ind[0]) {
		inv_slope = (vec[i][dir2] - pt_2) / (vec[i][dir1] - pt_1);
		if ((inv_slope > max_inv_slope)
		    && (fabs(inv_slope - max_inv_slope) > MACHINE_EPS)) {
		  ind[1] = i;
		  max_inv_slope = inv_slope;
		  num_rotated = 1;
		}
		else if (fabs(inv_slope - max_inv_slope) < MACHINE_EPS) {
		  ind[2] = i;
		  num_rotated++;
		}
	      }
	    }
	    break;
	  case CW:
	    for (i = 0; i < num_vec; i++) {
	      if (i != ind[0]) {
		inv_slope = (vec[i][dir2] - pt_2) / (vec[i][dir1] - pt_1);
		if ((inv_slope < min_inv_slope)
		    && (fabs(inv_slope - max_inv_slope) > MACHINE_EPS)) {
		  ind[1] = i;
		  min_inv_slope = inv_slope;
		  num_rotated = 1;
		}
		else if (fabs(inv_slope - min_inv_slope) < MACHINE_EPS) {
		  ind[2] = i;
		  num_rotated++;
		}
	      }
	    }
	    break;
	  }
	switch (num_rotated)
	  {
	  case 0:
	    SM_DEBUG_PRINT(3, "No points in the rotation ... odd\n")
	    ;
	    *status = HULL_TEST_ERROR;
	    break;
	  case 1:
	    SM_DEBUG_PRINT(3, "Found a line in the convex hull\n")
	    ;
	    SM_COPY_VECTOR3(pt2, vec[ind[1]])
	    ;
	    *status = TWO_PT_PLANE;
	    break;
	  default:
	    SM_DEBUG_PRINT(3, "Found 2 or more points in the rotation\n")
	    ;
	    if (fabs(pt_1) > MACHINE_EPS)
	      inv_origin_slope = pt_2 / pt_1;
	    switch (rotate)
	      {
	      case CCW:
		if (inv_origin_slope >= max_inv_slope)
		  *status = NO_EQUIL;
		else
		  *status = CHECK_TOP_DOWN;
		break;
	      case CW:
		if (inv_origin_slope <= min_inv_slope)
		  *status = NO_EQUIL;
		else
		  *status = CHECK_TOP_DOWN;
		break;
	      }
	    break;
	  }
	break;
      case 2: /* use these two points to define the plane */
	SM_DEBUG_PRINT(3, "Found two minimum points to define the plane\n")
	;
	SM_COPY_VECTOR3(pt1, vec[ind[0]])
	;
	SM_COPY_VECTOR3(pt2, vec[ind[1]])
	;
	*status = TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	SM_DEBUG_ACTION(
	    3,
	    {fprintf(stdout,"Found 3 or more points in min plane %f\n",min);})
	if (vec[ind[0]][dir1] >= 0)
	  *status = NO_EQUIL;
	else
	  *status = CHECK_TOP_DOWN;
	break;
      }
  }

  /***************************/
  /*  failed to find any information, checking top/down this coord*/
  /***************************/

  if (*status == CHECK_TOP_DOWN) {
    /* find the maximum points in dir1 starting at 1 */
    num_max = 0;
    ind[0] = -1;
    ind[1] = -1;
    ind[2] = -1;
    max = -1.0;
    for (i = 0; i < num_vec; i++) {
      if (vec[i][dir1] > max) {
	max = vec[i][dir1];
	ind[0] = i;
	num_max = 1;
      }
      else if (fabs(vec[i][dir1] - max) < MACHINE_EPS) {
	ind[num_max++] = i;
      }
    }
    if (max <= 0)
      *status = NO_EQUIL;

    if (*status != NO_EQUIL) {
      switch (num_max)
	{
	case 1: /* rotate to find the next point */
	  SM_COPY_VECTOR3(pt1, vec[ind[0]])
	  ;
	  pt_1 = pt1[dir1];
	  pt_2 = pt1[dir2];
	  if (pt1[dir2] < 0) {
	    rotate = CW;
	    min_inv_slope = A_BIG_POS_NUMBER;
	  }
	  if (pt1[dir2] >= 0) {
	    rotate = CCW;
	    max_inv_slope = A_BIG_NEG_NUMBER;
	  }
	  switch (rotate)
	    {
	    case CCW:
	      for (i = 0; i < num_vec; i++) {
		if (i != ind[0]) {
		  inv_slope = (vec[i][dir2] - pt_2) / (vec[i][dir1] - pt_1);
		  if (inv_slope > max_inv_slope) {
		    ind[1] = i;
		    max_inv_slope = inv_slope;
		    num_rotated = 1;
		  }
		  else if (fabs(inv_slope - max_inv_slope) < MACHINE_EPS) {
		    ind[2] = i;
		    num_rotated++;
		  }
		}
	      }
	      break;
	    case CW:
	      for (i = 0; i < num_vec; i++) {
		if (i != ind[0]) {
		  inv_slope = (vec[i][dir2] - pt_2) / (vec[i][dir1] - pt_1);
		  if (inv_slope < min_inv_slope) {
		    ind[1] = i;
		    min_inv_slope = inv_slope;
		    num_rotated = 1;
		  }
		  else if (fabs(inv_slope - min_inv_slope) < MACHINE_EPS) {
		    ind[2] = i;
		    num_rotated++;
		  }
		}
	      }
	      break;
	    }
	  switch (num_rotated)
	    {
	    case 0:
	      SM_DEBUG_PRINT(3, "No points in the rotation ... odd\n")
	      ;
	      *status = HULL_TEST_ERROR;
	      break;
	    case 1:
	      SM_DEBUG_PRINT(3, "Found a line in the convex hull\n")
	      ;
	      SM_COPY_VECTOR3(pt2, vec[ind[1]])
	      ;
	      *status = TWO_PT_PLANE;
	      break;
	    default:
	      SM_DEBUG_PRINT(3, "Found 2 or more points in the rotation\n")
	      ;
	      /* check to see if rotation got past origin */
	      inv_origin_slope = pt_2 / pt_1;
	      switch (rotate)
		{
		case CCW:
		  if (inv_origin_slope >= max_inv_slope)
		    *status = NO_EQUIL;
		  else if (dir1 == 2)
		    *status = CHECK_Y_COORD_DIRECTION;
		  else if (dir1 == 1)
		    *status = CHECK_X_COORD_DIRECTION;
		  else
		    *status = EQUIL;
		  break;
		case CW:
		  if (inv_origin_slope <= min_inv_slope)
		    *status = NO_EQUIL;
		  else if (dir1 == 2)
		    *status = CHECK_Y_COORD_DIRECTION;
		  else if (dir1 == 1)
		    *status = CHECK_X_COORD_DIRECTION;
		  else
		    *status = EQUIL;
		  break;
		}
	      break;
	    }
	  break;
	case 2: /* use these two points to define the plane */
	  SM_COPY_VECTOR3(pt1, vec[ind[0]])
	  ;
	  SM_COPY_VECTOR3(pt2, vec[ind[1]])
	  ;
	  *status = TWO_PT_PLANE;
	  break;
	default: /* check to see if all > 0 */
	  SM_DEBUG_ACTION(3, {fprintf(stdout,"Found 3 in max plane %f\n",max);})
	  ;
	  if (vec[ind[0]][dir1] <= 0)
	    *status = NO_EQUIL;
	  else if (dir1 == 2)
	    *status = CHECK_Y_COORD_DIRECTION;
	  else if (dir1 == 1)
	    *status = CHECK_X_COORD_DIRECTION;
	  else
	    *status = EQUIL;
	  break;
	}
    }
  }
}

