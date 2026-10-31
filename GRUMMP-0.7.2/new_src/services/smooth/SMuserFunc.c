#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

/********************************************************************************************************
 SMinitSmoothing
 *******************************************************************************************************/
/*@ SMinitSmoothing - Initializes the smoothing data structure and sets
 values for the smoothing technique, the mesh quality function,
 and threshold usage.

 Input Parameters:
 .  dimension - an integer indicating the problem dimension, either 2 or 3

 .  technique - an integer argument to set the smoothing technique used to
 adjust grid point location.
 .vb
 Input one of:
 OPTMS_LAPLACIAN_ONLY      (or L)
 OPTMS_OPTIMIZATION_ONLY   (or O)
 OPTMS_COMBINED            (or C)
 OPTMS_COMBINED1           (or 1)
 OPTMS_COMBINED2           (or 2)
 OPTMS_COMBINED3           (or 3)
 OPTMS_FLOATING_THRESHOLD  (or F)
 OPTMS_STUPID_LAPLACE      (or S)
 OPTMS_TECHNIQUE_DEFAULT   (or C)
 OPTMS_DEFAULT             (or -1)

 Note that either the COMBINED or FLOATING_THRESHOLD technique
 is recommended.  The COMBINED approach is the default.
 .ve

 .  FunctionID - an integer argument used to set the mesh quality measure to be
 optimized.

 .vb
 In 2D input one of
 MAX_MIN_ANGLE  (or 1): maximize the minimum angle
 MIN_MAX_COSINE (or 2): minimize the maximum cosine of the angle
 MAX_MIN_COSINE (or 3): maximize the minimum cosine of the angle
 MAX_MIN_SINE   (or 4): maximize the minimum sine of the angle
 MIN_MAX_ANGLE  (or 5): minimize the maximum angle
 MIN_MAX_JACOBIAN_DIFF (or 6): minimize the maximum square of the
 difference of the current jacobian and the jacobian of an equilateral
 triangle (scaled by the jacobian of an equilateral triangle)
 MAX_MIN_SCALED_JACOBIAN (or 7): maximize the minimum scaled jacobian
 for each of the three vertices of a triangle (J/(L1*L2)) where L1 and L2 are
 the lengths of the incident edges.  Same as MAX_MIN_SINE in the feasible
 region, but returns negative angles for inverted elements
 MAX_MIN_AREA_LENGTH_RATIO (or 8):  Computes the ratio of the the area
 of the triangle and the sum of the squares of the length of the edges
 MIN_MAX_LENGTH_AREA_RATIO (or 9): Computes the negtive inverse of the
 MAX_MIN_AREA_LENGTH_RATIO
 FUNCTION_DEFAULT_2D (which is MAX_MIN_SINE) (or 4)
 OPTMS_DEFAULT (-1, which will result in a choice of MAX_MIN_SINE)

 In 3D input one of
 MAX_MIN_DIHEDRAL (or 21): maximize the minimum angle
 MIN_MAX_DIHEDRAL (or 22): minimize the maximum angle
 MAX_MIN_COSINE_DIHEDRAL (or 23): maximize the minimum cosine of the angle
 MIN_MAX_COSINE_DIHEDRAL (or 24): minimize the maximum cosine of the angle
 MAX_SINE_DIHEDRAL (or 25): maximize the minimum sine of the angle
 MAX_MIN_SOLID (or 30): maximize the minimum solid angle
 FUNCTION_DEFAULT_3D (which is MAX_SINE_DIHEDRAL) (or 25)
 OPTMS_DEFAULT (-1, which will result in a choice of MAX_SINE_DIHEDRAL)

 The default in 2D is MAX_MIN_SINE and in 3D is MAX_SINE_DIHEDRAL.
 .ve
 -  Threshold - a double argument that sets the degree value of the threshold
 used in either the COMBINED technique, which has a fixed value, or
 the FLOATING_THRESHOLD technique, which allows the threshold to vary.
 The default values for the quality measures that depend on angle measures
 are 30 and 15 degrees for the COMBINED approach for 2D and 3D, respectively,
 and 10 and 15 degrees for the FLOATING_THRESHOLD approach.  For the
 measures that depend on Jacobian ratios the default is .25.

 Output Parameters:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing

 .seealso  SMsetSmoothTechnique(), SMsetSmoothFunction(), SMsetSmoothThreshold(),
 SMfinalizeSmoothing()
 @*/
void
SMinitSmoothing(int dimension, int technique, int FunctionID, double Threshold,
		void **ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  MY_MALLOCN(smooth_data, (SMsmooth_data *), sizeof(SMsmooth_data), 1);

  MY_MALLOC(smooth_data->smooth_param, (SMparam *), sizeof(SMparam), 1);

  (*ext_smooth_data) = (SMsmooth_data *) smooth_data;
  smooth_data->dimension = dimension;

  SMinitSmoothParam(technique, FunctionID, Threshold, (*ext_smooth_data));

  smooth_data->local_mesh = SMmallocLocalMesh(MAX_NUM_PTS,
					      (smooth_data)->smooth_param);
  smooth_data->local_mesh->dimension = dimension;
  smooth_data->local_mesh->opt_info->dimension = dimension;

  MY_MALLOCN(smooth_data->smooth_stats, (SMstats *), sizeof(SMstats), 1);
  smooth_data->smooth_stats->stats_initialized = 0;
  SMinitSmoothStats(smooth_data);
  smooth_data->smooth_procinfo = SMinitProcinfo();
  smooth_data->quality_table = SMmallocQualityTable();
}

/********************************************************************************************************
 SMsetProblemDimension
 *******************************************************************************************************/
/*@
 SMsetProblemDimension - This function allows the user to set the dimension of the smoothing
 problem.  Opt-MS currently supports 2D planar smoothing (in the x-y plane)
 and 3D smoothing.  This function call be called any time after SMinitSmoothing
 to set or change the dimension of the problem.

 Input Parameters:
 +  smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be called
 prior to calling this routine.
 -  dimension - an integer argument to set the dimension of the problem, either 2 or 3.

 .seealso SMinitSmoothing()
 @*/
void
SMsetProblemDimension(void *smooth_data, int dimension)
{
  SMsmooth_data *int_smooth_data;

  int_smooth_data = (SMsmooth_data *) smooth_data;
  int_smooth_data->dimension = dimension;
  int_smooth_data->local_mesh->dimension = dimension;
  int_smooth_data->local_mesh->opt_info->dimension = dimension;
}

void
SMsetFunctionID(void *m_smoothData, int FunctionID)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) m_smoothData;
  smooth_data->smooth_param->function_id = FunctionID;
  SMsetSmoothFunction(smooth_data, FunctionID);
}

/********************************************************************************************************
 SMsetSmoothTechnique
 *******************************************************************************************************/
/*@
 SMsetSmoothTechnique - This function allows the user to change the technique used for mesh
 smoothing at any time during the mesh improvement process.

 Input Parameters:
 +  smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be called
 prior to calling this routine.
 -  technique - a character argument to set the smoothing technique used to
 adjust grid point location.
 .vb
 Input one of:
 OPTMS_LAPLACIAN           (or S)
 OPTMS_SMART_LAPLACIAN     (or L)
 OPTMS_OPTIMIZATION_ONLY   (or O)
 OPTMS_COMBINED            (or C)
 OPTMS_COMBINED1           (or 1)
 OPTMS_COMBINED2           (or 2)
 OPTMS_COMBINED3           (or 3)
 OPTMS_FLOATING_THRESHOLD  (or F)
 OPTMS_STUPID_LAPLACE      (or S)
 OPTMS_TECHNIQUE_DEFAULT   (or C)
 OPTMS_DEFAULT             (or -1)

 Note that either the COMBINED or FLOATING_THRESHOLD technique
 is recommended.  The COMBINED approach is the default.

 Note that either the COMBINED or FLOATING_THRESHOLD technique
 is recommended.  The COMBINED approach is the default.
 .ve
 .seealso SMinitSmoothing()
 @*/

void
SMsetSmoothTechnique(void *ext_smooth_data, int technique)
{
  SMsmooth_data *smooth_data;
  int SmoothTechnique;

  smooth_data = (SMsmooth_data *) ext_smooth_data;

  switch (technique)
    {
    case STUPID_LAPLACE:
    case 'S':
    case 's':
      SM_DEBUG_PRINT(1,
		     "Setting the smoothing technique to regular LAPLACIAN\n")
      ;
      SmoothTechnique = STUPID_LAPLACE;
      break;
    case LAPLACE_ONLY:
    case 'L':
    case 'l':
      SM_DEBUG_PRINT(1,
		     "Setting the smoothing technique to SMART_LAPLACE_ONLY\n")
      ;
      SmoothTechnique = LAPLACE_ONLY;
      break;
    case OPTIMIZATION_ONLY:
    case 'O':
    case 'o':
      SM_DEBUG_PRINT(1,
		     "Setting the smoothing technique to OPTIMIZATION_ONLY\n")
      ;
      SmoothTechnique = OPTIMIZATION_ONLY;
      break;
    case COMBINED:
    case 'C':
    case 'c':
      SM_DEBUG_PRINT(1, "Setting the smoothing technique to COMBINED \n")
      ;
      SmoothTechnique = COMBINED;
      break;
    case COMBINED1:
    case '1':
      SM_DEBUG_PRINT(1, "Setting the smoothing technique to COMBINED1 \n")
      ;
      SmoothTechnique = COMBINED1;
      break;
    case COMBINED2:
    case '2':
      SM_DEBUG_PRINT(1, "Setting the smoothing technique to COMBINED2 \n")
      ;
      SmoothTechnique = COMBINED2;
      break;
    case COMBINED3:
    case '3':
      SM_DEBUG_PRINT(1, "Setting the smoothing technique to COMBINED3 \n")
      ;
      SmoothTechnique = COMBINED3;
      break;
    case FLOATING_THRESHOLD:
    case '4':
    case 'F':
    case 'f':
      SM_DEBUG_PRINT(
	  1,
	  "Setting the smoothing technique to COMBINED4 (FLOATING_THRESHOLD) \n")
      ;
      SmoothTechnique = FLOATING_THRESHOLD;
      break;
    default:
      SM_DEBUG_PRINT(
	  1, "You have entered an incorrect option or the default option\n")
      ;
      SM_DEBUG_PRINT(1, "Setting the smoothing technique to COMBINED\n")
      ;
      SmoothTechnique = COMBINED;
      break;
    }
  smooth_data->smooth_param->smooth_technique = SmoothTechnique;
}

/********************************************************************************************************
 SMsetSmoothFunction
 *******************************************************************************************************/
/*@
 SMsetSmoothFunction - This function allows the user to change the mesh quality function that
 is optimized at any time during the smoothing process.

 Input Parameters:
 +   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be called
 prior to calling this routine.
 -   FunctionID - an integer argument used to set the mesh quality measure to be
 optimized.
 .vb
 In 2D input one of
 MAX_MIN_ANGLE  (or 1): maximize the minimum angle
 MIN_MAX_COSINE (or 2): minimize the maximum cosine of the angle
 MAX_MIN_COSINE (or 3): maximize the minimum cosine of the angle
 MAX_MIN_SINE   (or 4): maximize the minimum sine of the angle
 MIN_MAX_ANGLE  (or 5): minimize the maximum angle
 MIN_MAX_JACOBIAN_DIFF (or 6): minimize the maximum square of the
 difference of the current jacobian and the jacobian of an equilateral
 triangle (scaled by the jacobian of an equilateral triangle)
 MAX_MIN_SCALED_JACOBIAN (or 7): maximize the minimum scaled jacobian
 for each of the three vertices of a triangle (J/(L1*L2)) where L1 and L2 are
 the lengths of the incident edges.  Same as MAX_MIN_SINE in the feasible
 region, but returns negative angles for inverted elements
 MAX_MIN_AREA_LENGTH_RATIO (or 8):  Computes the ratio of the the area
 of the triangle and the sum of the squares of the length of the edges
 MIN_MAX_LENGTH_AREA_RATIO (or 9): Computes the negtive inverse of the
 MAX_MIN_AREA_LENGTH_RATIO
 FUNCTION_DEFAULT_2D (which is MAX_MIN_SINE) (or 4)
 OPTMS_DEFAULT (-1, which will result in a choice of MAX_MIN_SINE)

 In 3D input one of
 MAX_MIN_DIHEDRAL (or 21): maximize the minimum angle
 MIN_MAX_DIHEDRAL (or 22): minimize the maximum angle
 MAX_MIN_COSINE_DIHEDRAL (or 23): maximize the minimum cosine of the angle
 MIN_MAX_COSINE_DIHEDRAL (or 24): minimize the maximum cosine of the angle
 MAX_SINE_DIHEDRAL (or 25): maximize the minimum sine of the angle
 FUNCTION_DEFAULT_3D (which is MAX_SINE_DIHEDRAL) (or 25)
 OPTMS_DEFAULT (-1, which will result in a choice of MAX_SINE_DIHEDRAL)

 The default in 2D is MAX_MIN_SINE and in 3D is MAX_SINE_DIHEDRAL.
 .ve
 .seealso SMinitSmoothing()
 @*/
void
SMsetSmoothFunction(void *ext_smooth_data, int FunctionID)
{
  SMsmooth_data *smooth_data;
  SMparam *smooth_param;
  int dimension;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_param = smooth_data->smooth_param;
  dimension = smooth_data->dimension;
  if (FunctionID == OPTMS_DEFAULT) {
    if (dimension == 2) {
      smooth_param->function_id = MAX_MIN_SINE;
    }
    else {
      smooth_param->function_id = MAX_MIN_SINE_DIHEDRAL;
    }
  }
  else {
    smooth_param->function_id = FunctionID;
  }

  /* set the default values */
  if (dimension == 2) {
    smooth_param->function_values_per_tri = 3;
  }
  else {
    smooth_param->function_values_per_tri = 6;
  }
  smooth_param->ComputeFunctionValues2D = SMcomputeTriSines;
  smooth_param->ComputeGradientValues2D = SMcomputeSineGradients;
  smooth_param->ComputeFunctionValues3D = vSineDihedrals;
  smooth_param->ComputeGradientValues3D = vGradSineDihedrals;

  switch (smooth_param->function_id)
    {
    case MAX_MIN_ANGLE:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeTriAngles;
      smooth_param->ComputeGradientValues2D = SMcomputeAngGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MAX_MIN_ANGLE\n")
      ;
      break;
    case MIN_MAX_ANGLE:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeNegTriAngles;
      smooth_param->ComputeGradientValues2D = SMcomputeNegAngGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MIN_MAX_ANGLE\n")
      ;
      break;
    case MIN_MAX_COSINE:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeNegTriCosines;
      smooth_param->ComputeGradientValues2D = SMcomputeNegCosGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MIN_MAX_COSINE\n")
      ;
      break;
    case MAX_MIN_COSINE:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeTriCosines;
      smooth_param->ComputeGradientValues2D = SMcomputeCosGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MAX_MIN_COSINE\n")
      ;
      break;
    case MAX_MIN_SINE:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeTriSines;
      smooth_param->ComputeGradientValues2D = SMcomputeSineGradients;
////	///////Quad////////
//	smooth_param->function_values_per_tri=4;
//	smooth_param->ComputeFunctionValuesQuad2D = SMcomputeQuadSines;
//	smooth_param->ComputeGradientValuesQuad2D = SMcomputeSineGradientsQuad;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MAX_MIN_SINE\n")
      ;
      break;
    case MIN_MAX_JACOBIAN_DIFF:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeTriJacobians;
      smooth_param->ComputeGradientValues2D = SMcomputeJacobianGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MIN_MAX_JACOBIAN_DIFF\n")
      ;
      break;
    case MAX_MIN_SCALED_JACOBIAN:
      smooth_param->function_values_per_tri = 3;
      smooth_param->ComputeFunctionValues2D = SMcomputeScaledTriJacobians;
      smooth_param->ComputeGradientValues2D = SMcomputeScaledJacobianGradients;
      SM_DEBUG_PRINT(1, "Setting the 2D function to MAX_MIN_SCALED_JACOBIAN\n")
      ;
      break;
    case MAX_MIN_AREA_LENGTH_RATIO:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeAreaLengthRatio;
      smooth_param->ComputeGradientValues2D = SMcomputeAreaLengthRatioGradients;
      SM_DEBUG_PRINT(1,
		     "Setting the 2D function to MAX_MIN_AREA_LENGTH_RATIO\n")
      ;
      break;
    case MIN_MAX_LENGTH_AREA_RATIO:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeLengthAreaRatio;
      smooth_param->ComputeGradientValues2D = NULL;
      /*        smooth_param->ComputeGradientValues2D = SMcomputeLengthAreaRatioGradients; */
      SM_DEBUG_PRINT(1,
		     "Setting the 2D function to MIN_MAX_LENGTH_AREA_RATIO\n")
      ;
      break;
    case MAX_MIN_INTERIOR_ANGLE:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriAngles;
      smooth_param->ComputeGradientValues2D = SMcomputeInteriorAngGradients;
      break;
    case MAX_MIN_INTERIOR_SINE:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriSines;
      smooth_param->ComputeGradientValues2D = SMcomputeInteriorSineGradients;
      break;
    case MAX_MIN_SQUARENESS:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeInteriorSquareness;
      smooth_param->ComputeGradientValues2D =
	  SMcomputeInteriorSquarenessGradients;
      break;
    case MAX_MIN_EDGERATIO:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeInteriorEdgeratio;
      smooth_param->ComputeGradientValues2D =
	  SMcomputeInteriorEdgeratioGradients;
      break;
    case MIN_MAX_INTERIOR_COSINE:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMcomputeInteriorTriCosines;
      smooth_param->ComputeGradientValues2D = SMcomputeInteriorCosGradients;
      break;
    case MAX_MIN_INTERIOR_SCALED_JACOBIAN:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D =
	  SMcomputeInteriorScaledTriJacobians;
      smooth_param->ComputeGradientValues2D =
	  SMcomputeInteriorScaledJacobianGradients;
      break;
    case MIN_MAX_NORM_JAC_SQUARED_2D:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues2D = SMNormJacSquared2D;
      smooth_param->ComputeGradientValues2D =
	  SMcomputeNormJacSquaredGradients2D;
      break;

    case MAX_MIN_DIHEDRAL:
      smooth_param->function_values_per_tri = 6;
      smooth_param->ComputeFunctionValues3D = vDihedrals;
      smooth_param->ComputeGradientValues3D = vGradDihedrals;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MAX_MIN_DIHEDRAL\n")
      ;
      break;
    case MIN_MAX_DIHEDRAL:
      smooth_param->function_values_per_tri = 6;
      smooth_param->ComputeFunctionValues3D = vNegateDihedrals;
      smooth_param->ComputeGradientValues3D = vNegateGradDihedrals;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MIN_MAX_DIHEDRAL\n")
      ;
      break;
    case MIN_MAX_COSINE_DIHEDRAL:
      smooth_param->function_values_per_tri = 6;
      smooth_param->ComputeFunctionValues3D = vNegateCosineDihedrals;
      smooth_param->ComputeGradientValues3D = vNegateGradCosineDihedrals;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MIN_MAX_COSINE_DIHEDRAL\n")
      ;
      break;
    case MAX_MIN_COSINE_DIHEDRAL:
      smooth_param->function_values_per_tri = 6;
      smooth_param->ComputeFunctionValues3D = vCosineDihedrals;
      smooth_param->ComputeGradientValues3D = vGradCosineDihedrals;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MAX_MIN_COSINE_DIHEDRAL\n")
      ;
      break;
    case MAX_MIN_SINE_DIHEDRAL:
      smooth_param->function_values_per_tri = 6;
      smooth_param->ComputeFunctionValues3D = vSineDihedrals;
      smooth_param->ComputeGradientValues3D = vGradSineDihedrals;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MAX_SINE_DIHEDRAL\n")
      ;
      break;
    case MAX_MIN_SCALED_JACOBIAN_3D:
      smooth_param->function_values_per_tri = 4;
      smooth_param->ComputeFunctionValues3D = vScaledJacobian;
      smooth_param->ComputeGradientValues3D = vGradScaledJacobian;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MAX_MIN_SCALED_JACOBIAN\n")
      ;
      break;
    case MIN_MAX_SRMS_VOLUME_RATIO:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues3D = vSMRSVolumeRatio;
      smooth_param->ComputeGradientValues3D = vGradSMRSVolumeRatio;
      SM_DEBUG_PRINT(1,
		     "Setting the 3D function to MIN_MAX_SRMS_VOLUME_RATIO\n")
      ;
      break;
    case MIN_MAX_NORM_JAC_SQUARED_3D:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues3D = vNormJacSquared;
      smooth_param->ComputeGradientValues3D = vGradNormJacSquared;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MIN_MAX_NORM_JAC_SQUARED\n")
      ;
      break;
    case MIN_MAX_CONDITION:
      smooth_param->function_values_per_tri = 1;
      smooth_param->ComputeFunctionValues3D = vCondition;
      smooth_param->ComputeGradientValues3D = vGradCondition;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MIN_MAX_CONDITION\n")
      ;
      break;
    case MAX_MIN_SOLID:
      smooth_param->function_values_per_tri = 4;
      smooth_param->ComputeFunctionValues3D = vSolids;
      smooth_param->ComputeGradientValues3D = vGradSolids;
      SM_DEBUG_PRINT(1, "Setting the 3D function to MAX_MIN_SOLID\n")
      ;
      break;

    default:
      SM_DEBUG_PRINT(1, "Using the default functions in both 2 and 3D\n")
      ;
      break;
    }
}

void
SMset2DUserQualityFunction(void *ext_smooth_data, int values_per_tri,
			   SMfunction_ptr2D userQualityFunc,
			   SMgradfunc_ptr2D userQualityGrad)
{
  SMsmooth_data *smooth_data;
  SMparam *smooth_param;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_param = smooth_data->smooth_param;
  smooth_param->function_values_per_tri = values_per_tri;
  smooth_param->ComputeFunctionValues2D = userQualityFunc;
  smooth_param->ComputeGradientValues2D = userQualityGrad;
  SM_DEBUG_PRINT(1, "Setting the quality function to user function \n");
}

/********************************************************************************************************
 SMsetSmoothThreshold
 *******************************************************************************************************/
/*@
 SMsetSmoothThreshold - This function allows the user to change the value of the threshold
 when using the combined or the FLOATING_THRESHOLD techniques at any
 time in the smoothing process.

 Input Parameters:
 +   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be called
 prior to calling this routine.
 -  Threshold - a double argument that sets the degree value of the threshold
 used in either the COMBINED technique, which has a fixed value, or
 the FLOATING_THRESHOLD technique, which allows the threshold to vary.
 The default values for the quality measures that depend on angles are 15 and
 30 degrees for the COMBINED approach for 3D and 2D respectively, and 10 and
 15 degrees for the FLOATING_THRESHOLD approach.  For the measures that
 depend on Jacobian ratios the default is .25.

 .seealso SMinitSmoothing()
 @*/
void
SMsetSmoothThreshold(void *ext_smooth_data, double Threshold)
{
  SMsmooth_data *smooth_data;
  SMparam *smooth_param;
  int technique;
  int function_id;
  double radians = 0; /* Init for picky compilers */
  double add_value;
  int dimension;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_param = smooth_data->smooth_param;
  dimension = smooth_data->dimension;
  technique = smooth_param->smooth_technique;
  function_id = smooth_param->function_id;

  /***********************************************************
   set the thresholds for the functions that rely on angles
   for their measure of quality
   ***********************************************************/

  if (function_id == MAX_MIN_ANGLE || function_id == MIN_MAX_COSINE
      || function_id == MAX_MIN_COSINE || function_id == MAX_MIN_SINE
      || function_id == MIN_MAX_ANGLE || function_id == MAX_MIN_SCALED_JACOBIAN
      || function_id == MAX_MIN_DIHEDRAL || function_id == MIN_MAX_DIHEDRAL
      || function_id == MAX_MIN_SOLID || function_id == MAX_MIN_COSINE_DIHEDRAL
      || function_id == MIN_MAX_COSINE_DIHEDRAL
      || function_id == MAX_MIN_SINE_DIHEDRAL) {

    if (dimension == 2) {
      /* set the default values */
      if (Threshold <= OPTMS_DEFAULT) {
	if (technique == FLOATING_THRESHOLD)
	  Threshold = 15.;
	else
	  Threshold = 30.;
      }
      add_value = 5.;
    }
    else {
      if (Threshold <= OPTMS_DEFAULT) {
	if (technique == FLOATING_THRESHOLD)
	  Threshold = 10.;
	else
	  Threshold = 15.;
      }
      add_value = 5.;
      if (function_id == MAX_MIN_SOLID) {
	Threshold /= 5;
	add_value /= 5;
      }
    }

    /* set the radians to the default or passed-in value */
    radians = Threshold * 3.14159 / 180;

    /* modify if necessary for the floating technique */
    if ((technique == FLOATING_THRESHOLD)
	&& (smooth_param->global_min_value != INIT_MIN_VALUE)) {
      SMconvertToDegrees(smooth_param->function_id,
			 &(smooth_param->global_min_value));
      /* if it's valid */
      if (smooth_param->global_min_value < 180
	  && smooth_param->global_min_value > 0) {
	radians = (smooth_param->global_min_value + add_value) * 3.14159 / 180.;
	SM_DEBUG_ACTION(
	    1,
	    { fprintf(stdout,"Set floating threshold: min value %f threshold %f\n", smooth_param->global_min_value, smooth_param->global_min_value + add_value); });
      }
    }
  }

  /***********************************************************
   set the thresholds for the functions that rely on jacobians
   for their measure of quality (note that scaled jacobian is really sine)
   ***********************************************************/

  if ((function_id == MIN_MAX_JACOBIAN_DIFF)
      || (function_id == MAX_MIN_AREA_LENGTH_RATIO)) {
    if (Threshold <= OPTMS_DEFAULT) {
      Threshold = .25;
    }
    if (Threshold > 1.0) {
      Threshold = .25;
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"User defined threshold must be less than 1, setting to the default value\n"); });
    }
  }

  switch (smooth_param->function_id)
    {
    case MAX_MIN_ANGLE:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_ANGLE: %f\n",radians); })
      ;
      smooth_param->lap_accept_value = radians;
      break;
    case MIN_MAX_ANGLE:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MIN_MAX_ANGLE: %f\n",radians); })
      ;
      smooth_param->lap_accept_value = radians;
      break;
    case MIN_MAX_COSINE:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MIN_MAX_COSINE: %f\n", cos(radians)); })
      ;
      smooth_param->lap_accept_value = cos(radians);
      break;
    case MAX_MIN_COSINE:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_COSINE: %f\n", cos(radians)); })
      ;
      smooth_param->lap_accept_value = cos(radians);
      break;
    case MAX_MIN_SINE:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_SINE: %f\n", sin(radians)); })
      ;
      smooth_param->lap_accept_value = sin(radians);
      break;
    case MAX_MIN_SCALED_JACOBIAN:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_SCALED_JACOBIAN: %f\n", sin(radians)); })
      ;
      smooth_param->lap_accept_value = sin(radians);
      break;
    case MIN_MAX_JACOBIAN_DIFF:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MIN_MAX_JACOBIAN_DIFF: %f\n", Threshold); })
      ;
      smooth_param->lap_accept_value = Threshold;
      break;
    case MAX_MIN_AREA_LENGTH_RATIO:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_AREA_LENGTH_RATIO: %f\n", Threshold); })
      ;
      smooth_param->lap_accept_value = Threshold;
      break;
    case MAX_MIN_DIHEDRAL:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_DIHEDRAL: %f\n", radians); })
      ;
      smooth_param->lap_accept_value = radians;
      break;
    case MIN_MAX_DIHEDRAL:
      SM_DEBUG_ACTION(
	  1,
	  { printf("Setting the threshold value for MIN_MAX_DIHEDRAL: %f\n", radians); })
      ;
      smooth_param->lap_accept_value = radians;
      break;
    case MIN_MAX_COSINE_DIHEDRAL:
      SM_DEBUG_ACTION(
	  1,
	  { printf("Setting the threshold value for MIN_MAX_COSINE_DIHEDRAL: %f\n",cos(radians)); })
      ;
      smooth_param->lap_accept_value = cos(radians);
      break;
    case MAX_MIN_COSINE_DIHEDRAL:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_COSINE_DIHEDRAL: %f\n",cos(radians)); })
      ;
      smooth_param->lap_accept_value = cos(radians);
      break;
    case MAX_MIN_SINE_DIHEDRAL:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_SINE_DIHEDRAL: %f\n", sin(radians)); })
      ;
      smooth_param->lap_accept_value = sin(radians);
      break;
    case MAX_MIN_SOLID:
      SM_DEBUG_ACTION(
	  1,
	  { fprintf(stdout,"Setting the threshold value for MAX_MIN_SOLID: %f\n", radians); })
      ;
      smooth_param->lap_accept_value = radians;
      break;
    default:
      if (dimension == 2) {
	SM_DEBUG_PRINT(
	    1, "The smoothing function threshold has not been defined\n");
	SM_DEBUG_PRINT(1,
		       "Setting the default acceptable function value to .5\n");
	smooth_param->lap_accept_value = sin(radians);
      }
      else {
	SM_DEBUG_PRINT(
	    1, "The smoothing function threshold has not been defined\n");
	SM_DEBUG_PRINT(1,
		       "Setting the default acceptable function value to .5\n");
	smooth_param->lap_accept_value = sin(radians);
      }
      break;
    }
}

/********************************************************************************************************
 SMinitGlobalMinValue
 *******************************************************************************************************/
/*@

 SMinitGlobalMinValue - This routine initializes the global minimum value of the quality metric
 to a very large value.   This is used in the floating threshold technique in which the
 minimum value of the quality metric is tracked for the next iteration.

 Input Parameter:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing which must be
 called  prior to calling this routine.  SMinitSmoothingStats must also have been
 called prior to calling this routine.

 .seealso SMinitSmoothing(), SMsetSmoothThreshold()
 @*/
void
SMinitGlobalMinValue(void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;
  SM_DEBUG_PRINT(1, "Setting the global minimum value to 1E300\n");
  smooth_data->smooth_param->global_min_value = INIT_MIN_VALUE;
}

/********************************************************************************************************
 SMinitSmoothStats
 *******************************************************************************************************/
/*@
 SMinitSmoothStats - If statistics gathering has been enabled in the configure process,
 then this routine intializes the statistics structure that records the number
 of cells smoothed, the number of equilibrium points found, and the reason for
 algorithm termination.

 Input Parameter:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.

 Note:
 It is useful to call this routine in conjunction with SMprintSmoothStats during each
 global pass over the mesh so that the incremental improvment can be monitored.

 .seealso SMprintSmoothStats()
 @*/
void
SMinitSmoothStats(void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;
  STATS_ON( {
    SMinitStats(smooth_data->smooth_stats)
    ;
  STATS_OFF});
}

/********************************************************************************************************
 SMprintSmoothStats
 *******************************************************************************************************/
/*@
 SMprintSmoothStats - If statistics gathering has been enabled in the configure process,
 then this routine prints the statistics that have been accumulated by the smoothing
 code since the last call to SMinitSmoothStats.

 Input Parameter:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing which must be
 called  prior to calling this routine.  SMinitSmoothingStats must also have been
 called prior to calling this routine.

 Note:
 .vb
 The information printed includes
 - the total number of nodes smoothed,
 - the number for which Laplacian smoothing was used and the number
 of the those that resulted in an invalid mesh and/or no improvement to
 the mesh
 - the number of nodes for which optimization-based smoothing was used
 (including the average iteration count),
 - the number of cells with no improvement,
 - the averate active value and average improvement, and
 - the termination status for the cells that were smoothed.
 .ve

 .seealso SMinitSmoothing(), SMinitSmoothStats()
 @*/
void
SMprintSmoothStats(void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;
  STATS_ON(
      { if (smooth_data->smooth_stats->stats_initialized) { SMprintStats(smooth_data); } else { fprintf(stderr,"[SM] You did not initialize stats for printing\n"); } STATS_OFF});
}

/********************************************************************************************************
 SMinitQualityTable
 *******************************************************************************************************/
/*@
 SMinitQualityTable - This function allows the user to take advantage of the
 quality metrics implemented in the OptMS code to analyze the quality of their
 mesh.

 Input Parameter:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.   This routine should be called before each
 global pass of measuring quality or information from the previous pass will
 continue to be accumulated.

 Note:
 .n      In 2D the min, max, and average values of the triangle angles, deviation from
 an equilateral triangle, the scaled jacobians, and the triangle area are printed.
 .n
 .n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
 scaled Jacobians, the ratio of sum of the squares of the length of the edges raised
 to the 3/2 power to the volume, and the tetrahdral volume are printed.
 .n
 .n      If any triangle areas or tetrahedral volumes are negative, the mesh is considered
 to be invalid, and SMuntangle should be called to try to create a valid mesh.

 .seealso SMinitSmoothing(), SMuntangle(), SMaccumulateQualityInformation(),
 SMprintQualityInformation()
 @*/
void
SMinitQualityTable(void *ext_smooth_data)
{
  int i;
  SMsmooth_data *smooth_data;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_data->quality_table->initialized = TRUE;
  smooth_data->quality_table->num_tangled_elements = 0;
  ;
  smooth_data->quality_table->mesh_validity = 1;

  for (i = 0; i < smooth_data->quality_table->num_functions; i++) {
    smooth_data->quality_table->measure[i]->min_value = A_BIG_POS_NUMBER;
    smooth_data->quality_table->measure[i]->max_value = A_BIG_NEG_NUMBER;
    smooth_data->quality_table->measure[i]->avg_min_value = 0.0;
    smooth_data->quality_table->measure[i]->avg_max_value = 0.0;
    smooth_data->quality_table->measure[i]->avg_value = 0.0;
    smooth_data->quality_table->measure[i]->num_function_values = 0;
    smooth_data->quality_table->measure[i]->num_elements = 0;
  }
}

/********************************************************************************************************
 SMaccumulateQualityInformation
 *******************************************************************************************************/
/*@
 SMaccumulateQualityInformation - This function computes the quality information
 for a right-handed triangle or tetrahedra.

 Input Parameters:
 +   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.
 -   vtx - a matrix containing the coordinates of the nodes of the triangle or tetrahedra.
 Of dimension 3 x 2 for triangles, and 4 x 3 for tetrahedra.

 Notes:
 In 2D the min, max, and average values of the triangle angles, deviation from
 an equilateral triangle, the scaled jacobians, and the triangle area are computed.
 .n
 .n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
 scaled Jacobians, the ratio of sum of the squares of the length of the edges raised
 to the 3/2 power to the volume, and the tetrahdral volume are computed.

 .seealso SMinitSmoothing(), SMinitQualityTable(), SMprintQualityInformation()
 @*/
void
SMaccumulateQualityInformation(void *ext_smooth_data, double **vtx)
{
  int i, j;
  int num_values;
  int dimension;
  double function[6];

  SMsmooth_data *smooth_data;
  SMquality_table *quality_table;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  dimension = smooth_data->dimension;
  quality_table = smooth_data->quality_table;

  for (i = 0; i < quality_table->num_functions; i++) {
    /* 2D functions */
    if (dimension == 2) {
      switch (i)
	{
	case 0:
	  SMcomputeTriAngles(vtx[0], vtx[1], vtx[2], function, &num_values);
	  for (j = 0; j < num_values; j++)
	    function[j] = function[j] / 3.14159 * 180;
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	case 1:
	  SMcomputeTriJacobians(vtx[0], vtx[1], vtx[2], function, &num_values);
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	case 2:
	  SMcomputeScaledTriJacobians(vtx[0], vtx[1], vtx[2], function,
				      &num_values);
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	case 3:
	  SMcomputeTriArea(vtx[0], vtx[1], vtx[2], function, &num_values);
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  if (function[0] < 0) {
	    quality_table->num_tangled_elements++;
	    quality_table->mesh_validity = 0;
	  }
	  break;
	}
    }
    if (dimension == 3) {
      switch (i)
	{
	/* 3D functions */
	case 4:
	  vDihedrals(vtx[0], vtx[1], vtx[2], vtx[3], function, &num_values);
	  for (j = 0; j < num_values; j++)
	    function[j] = function[j] / 3.14159 * 180;
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	case 5:
	  vScaledJacobian(vtx[0], vtx[1], vtx[2], vtx[3], function,
			  &num_values);
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	case 6:
	  vSMRSVolumeRatio(vtx[0], vtx[1], vtx[2], vtx[3], function,
			   &num_values);
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  for (j = 0; j < num_values; j++)
	    function[j] = -function[j];
	  break;
	case 7:
	  vComputeTetVolume(vtx[0], vtx[1], vtx[2], vtx[3], function,
			    &num_values);
	  if (function[0] < -1E16) {
	    printf("Bad things in tet volume calculation %f \n", function[0]);
	  }
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  if (function[0] < 0) {
	    quality_table->num_tangled_elements++;
	    quality_table->mesh_validity = 0;
	  }
	  break;
	case 8:
	  vCondition(vtx[0], vtx[1], vtx[2], vtx[3], function, &num_values);
	  for (j = 0; j < num_values; j++)
	    function[j] = -function[j];
	  SMinsertQualityInfo(quality_table, i, function, num_values);
	  break;
	}
    }
  }
}

void
SMinsertQualityInfo(SMquality_table *quality_table, int measure_id,
		    double *function, int num_values)
{
  int i;
  double element_min = A_BIG_POS_NUMBER;
  double element_max = A_BIG_NEG_NUMBER;

  quality_table->measure[measure_id]->num_function_values += num_values;
  quality_table->measure[measure_id]->num_elements++;

  for (i = 0; i < num_values; i++) {

    if (function[i] < element_min)
      element_min = function[i];
    if (function[i] > element_max)
      element_max = function[i];

    quality_table->measure[measure_id]->avg_value += function[i];
  }
  if (element_min < quality_table->measure[measure_id]->min_value) {
    quality_table->measure[measure_id]->min_value = element_min;
  }
  if (element_max > quality_table->measure[measure_id]->max_value) {
    quality_table->measure[measure_id]->max_value = element_max;
  }
  quality_table->measure[measure_id]->avg_min_value += element_min;
  quality_table->measure[measure_id]->avg_max_value += element_max;

}

/********************************************************************************************************
 SMprintQualityInformation
 *******************************************************************************************************/
/*@
 SMprintQualityInformation - This function prints the quality information accumulated
 since the last call to SMinitQualityInformation.

 Input Parameters:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.

 Notes:
 In 2D the min, max, and average values of the triangle angles, deviation from
 an equilateral triangle, the scaled jacobians, and the triangle area are printed.
 .n
 .n      In 3D, the min, max, and average values of the tetrahedral dihedral angles,
 scaled jacobians, the ratio of sum of the squares of the length of the edges raised
 to the 3/2 power to the volume, and the tetrahdral volume are printed.
 .n
 .n      If any triangle areas or tetrahedral volumes are negative, the mesh is considered
 to be invalid, and SMuntangle should be called to try to create a valid mesh.

 .seealso SMinitSmoothing(), SMinitQualityInformation(), SMaccumulateQualityInformation()
 @*/
void
SMprintQualityInformation(void *ext_smooth_data)
{
  int i;
  int dimension;
  int measure_ind = 0; /* Init for picky compilers. */
  SMsmooth_data *smooth_data;
  SMquality_table *quality_table;
  SMqualityMeasure *measure;

  smooth_data = (SMsmooth_data *) ext_smooth_data;
  dimension = smooth_data->dimension;
  quality_table = smooth_data->quality_table;

  if (dimension == 2)
    measure_ind = 0;
  if (dimension == 3)
    measure_ind = 4;

  printf("\nMesh Quality Information for %d Elements\n",
	 quality_table->measure[measure_ind]->num_elements);
  printf(
      "-------------------------------------------------------------------------------------\n");
  printf(
      "Quality Metric (Target)           Min Value  Max Value  Avg Value  Avg Min    Avg Max\n");
  printf(
      "-------------------------------------------------------------------------------------\n");
  for (i = 0; i < quality_table->num_functions; i++) {
    if (quality_table->measure[i]->dimension == dimension) {
      measure = quality_table->measure[i];
      printf("%20s (%3.2e) %10.2e %10.2e %10.2e %10.2e %10.2e (%d)\n",
	     measure->name, measure->target, measure->min_value,
	     measure->max_value,
	     measure->avg_value / measure->num_function_values,
	     measure->avg_min_value / measure->num_elements,
	     measure->avg_max_value / measure->num_elements,
	     measure->num_function_values);
    }
  }
  printf("-----------------------------------------------------------\n");
  printf("There are %d invalid elements in the mesh \n\n",
	 quality_table->num_tangled_elements);
}

/********************************************************************************************************
 SMinvalidMesh
 *******************************************************************************************************/
/*@
 SMinvalidMesh - This function returns a Boolean value after testing to see if the
 mesh contains any invalid elements.  The test is based on information contained
 in the quality table and this information must be accumulated before this test can
 be made.  This is useful in determining if the mesh requires untangling.

 Input Parameters:
 .   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.

 Return Value:
 a Boolean value of 1 if the mesh contains invalid elements or 0 if the
 mesh contains no invalid elements.

 .seealso SMinitSmoothing, SMinitQualityInformation(), SMaccumulateQualityInformtaion,
 SMuntangle()
 @*/
int
SMinvalidMesh(void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;
  if (smooth_data->quality_table->mesh_validity == 0)
    return (1);
  return (0);
}

/********************************************************************************************************
 SMsetMeshValidity
 *******************************************************************************************************/
/*@
 SMsetMeshValidity - This function allows the user to set the mesh validity without using the
 quality functions provided by the Opt-MS package.  This is useful if the user knows the
 mesh requires untangling and doesn't want to evaluate the quality of the entire mesh to
 determine if there are inverted elements.

 Input Parameters:
 +   mesh_validity - a Boolean value of 1 if the mesh is valid and 0 if the mesh contains
 elements with negative area
 -   smooth_data - a void data structure that contains the context and data structures
 for smoothing.  This structure is created in SMinitSmoothing, which must be
 called  prior to calling this routine.

 .seealso SMinitSmoothing, SMinitQualityInformation(), SMaccumulateQualityInformtaion,
 SMuntangle()
 @*/
void
SMsetMeshValidity(int mesh_validity, void *ext_smooth_data)
{
  SMsmooth_data *smooth_data;
  smooth_data = (SMsmooth_data *) ext_smooth_data;
  smooth_data->quality_table->mesh_validity = mesh_validity;
}
