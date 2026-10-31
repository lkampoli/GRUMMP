#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SMsmooth.h"

int
SMformGrammian(int num_vecs, double **vec, double **G, int dimension)
{
  int i, j;

  if (num_vecs > MAX_G_NUM) {
    fprintf(stderr, "Exceeded maximum allowed active values\n");
    fprintf(stderr, "Maximum allowed %d You requested %d \n",
    MAX_G_NUM,
	    num_vecs);
    fprintf(stderr, "Change MAX_G_NUM in SMsmooth.h \n");
    return (1);
  }
  /* form the grammian with the dot products of the gradients */
  for (i = 0; i < num_vecs; i++) {
    for (j = i; j < num_vecs; j++) {
      G[i][j] = 0.;
      SM_DOT(G[i][j], vec[i], vec[j], dimension);
      G[j][i] = G[i][j];
    }
  }
  return (0);
}

void
SMformPDGrammian(SMoptimal *opt_info)
{
  int i, k, g_ind_1;
  int dimension;
  int num_active, *active_ind, *PDG_ind;
  double **G, **PDG;

  dimension = opt_info->dimension;

  /* this assumes the grammian has been formed */
  num_active = opt_info->active->num_active;
  active_ind = opt_info->active->active_ind;
  PDG_ind = opt_info->PDG_ind;
  G = opt_info->G;
  PDG = opt_info->PDG;

  /* use the first gradient in the active set */
  g_ind_1 = 0;
  PDG[0][0] = G[0][0];
  PDG_ind[0] = active_ind[0];

  /* test the rest and add them as appropriate */
  k = 1;
  i = 1;
  while ((k < dimension) && (i < num_active)) {
    PDG[0][k] = PDG[k][0] = G[0][i];
    PDG[k][k] = G[i][i];
    if (k == 2) { /* add the dot product of g1 and g2 */
      PDG[1][k] = PDG[k][1] = G[g_ind_1][i];
    }
    if (SMnonSingularTest(k + 1, PDG)) {
      PDG_ind[k] = active_ind[i];
      if (k == 1)
	g_ind_1 = i;
      k++;
    }
    i++;
  }
  opt_info->num_LI = k;
}

double **
SMformVerticalMatrix(int num, double **PDG)
{
  int i, j;
  double **N;

  MY_MALLOCN(N, (double **), sizeof(double *) * num, 1);
  for (i = 0; i < num; i++)
    MY_MALLOCN(N[i], (double *), sizeof(double) * num, 1);

  for (i = 0; i < num; i++) {
    N[i][i] = PDG[i][i] + 1;
    for (j = i + 1; j < num; j++) {
      N[i][j] = N[j][i] = PDG[i][j] + 1;
    }
  }
  return (N);
}

double **
SMformReducedMatrix(int num_active, double **G)
{
  int i, j;
  double **P;

assert(num_active >= 1);
  			    MY_MALLOCN(P, (double **), sizeof(double *) * (num_active - 1), 1);
  for (i = 0; i < num_active - 1; i++)
    MY_MALLOCN(P[i], (double *), sizeof(double) * (num_active - 1), 1);

  for (i = 0; i < num_active - 1; i++) {
    P[i][i] = G[0][0] - 2 * G[0][i + 1] + G[i + 1][i + 1];
    for (j = i + 1; j < num_active - 1; j++) {
      P[i][j] = G[0][0] - G[0][j + 1] - G[i + 1][0] + G[i + 1][j + 1];
      P[j][i] = P[i][j];
    }
  }
  return (P);
}

int
SMnonSingularTest(int n, double **A)
{
  int test;
  double cond;

  test = FALSE;
  switch (n)
    {
    case 1:
      if (A[0][0] > 0)
	test = TRUE;
      break;
    case 2:
      if (fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]) > MACHINE_EPS)
	test = TRUE;
      break;
    case 3:
      /* calculate the condition number */
      cond = SMcondition3x3(A);
      if (cond < 1E14)
	test = TRUE;
      /*         if (cond > 1E13 && cond <= 1E14) printf("Nearly singular 3x3 matrix: condition number %e\n",cond);*/
      break;
    default:
      printf("This PD test doesn't work for n > 3");
      break;
    }
  return (test);
}

double *
SMsolve2x2(double a11, double a12, double a21, double a22, double b1, double b2)
{
  double *x;
  double factor;

  /* if the system is not singular, solve it */
  if (fabs(a11 * a22 - a21 * a12) > MACHINE_EPS) {
    MY_MALLOCN(x, (double *), sizeof(double) * 2, 1);
    if (fabs(a11) > MACHINE_EPS) {
      factor = (a21 / a11);
      x[1] = (b2 - factor * b1) / (a22 - factor * a12);
      x[0] = (b1 - a12 * x[1]) / a11;
    }
    else if (fabs(a21) > MACHINE_EPS) {
      factor = (a11 / a21);
      x[1] = (b1 - factor * b2) / (a12 - factor * a22);
      x[0] = (b2 - a22 * x[1]) / a21;
    }
  }
  else {
    x = NULL;
  }
  return x;
}

void
SMsolve2x2_t(double a11, double a12, double a21, double a22, double b1,
	     double b2, double *x)
{
  double factor;

  /* if the system is not singular, solve it */
  if (fabs(a11 * a22 - a21 * a12) > MACHINE_EPS) {
    if (fabs(a11) > MACHINE_EPS) {
      factor = (a21 / a11);
      x[1] = (b2 - factor * b1) / (a22 - factor * a12);
      x[0] = (b1 - a12 * x[1]) / a11;
    }
    else if (fabs(a21) > MACHINE_EPS) {
      factor = (a11 / a21);
      x[1] = (b1 - factor * b2) / (a12 - factor * a22);
      x[0] = (b2 - a22 * x[1]) / a21;
    }
  }
  else {
    printf("Singular 2x2 system!\n");
    x = NULL;
  }
}

void
SMsolve3x3(double **A, double *B, double *x)
{
  int pivot, row1, row2;
  double L10, L20, a22_bar, a23_bar, a32_bar, a33_bar;
  double a33_dbar, b2_bar, b3_bar, b3_dbar;
  double a23_dbar, b2_dbar;

  /* if nonsingular, solve it */
  if (SMnonSingularTest(3, A)) {
    if (A[0][0] != 0) {
      pivot = 0;
      row1 = 1;
      row2 = 2;
    }
    else if (A[1][0] != 0) {
      pivot = 1;
      row1 = 0;
      row2 = 2;
    }
    else {
      pivot = 2;
      row1 = 0;
      row2 = 1;
    }

    L10 = A[row1][pivot] / A[pivot][pivot];
    L20 = A[row2][pivot] / A[pivot][pivot];
    a22_bar = A[row1][row1] - L10 * A[pivot][row1];
    a23_bar = A[row1][row2] - L10 * A[pivot][row2];
    a32_bar = A[row2][row1] - L20 * A[pivot][row1];
    a33_bar = A[row2][row2] - L20 * A[pivot][row2];

    b2_bar = B[row1] - L10 * B[pivot];
    b3_bar = B[row2] - L20 * B[pivot];

    if (a22_bar != 0) {
      a33_dbar = a33_bar - (a32_bar / a22_bar) * a23_bar;
      b3_dbar = b3_bar - a32_bar / a22_bar * b2_bar;
      x[row2] = b3_dbar / a33_dbar;
      x[row1] = (b2_bar - a23_bar * x[row2]) / a22_bar;
    }
    else {
      a23_dbar = a23_bar - (a22_bar / a32_bar) * a33_bar;
      b2_dbar = b2_bar - a22_bar / a32_bar * b3_bar;
      x[row1] = b2_dbar / a23_dbar;
      x[row2] = (b3_bar - a33_bar * x[row2]) / a32_bar;
    }
    x[pivot] = (B[pivot] - A[pivot][row1] * x[row1] - A[pivot][row2] * x[row2])
	/ A[pivot][pivot];
  }
  else {
    /* return an error */
    x = NULL;
    printf("Trying to solve a singular 3x3 system \n");
  }
}

double *
SMsolveSPD3x3(double **A, double *B)
{
  double *x;
  double L11, L21, L31, L22, L32, L33;
  double Y1, Y2, Y3;

  /* solves a 3x3 symmetric positive definite system 
   using cholesky factorization */

  MY_MALLOCN(x, (double *), sizeof(double) * 3, 1);

  /* Find the cholesky factor L */
  L11 = sqrt(A[0][0]);
  L21 = A[1][0] / L11;
  L22 = sqrt(A[1][1] - L21 * L21);
  L31 = A[2][0] / L11;
  L32 = (A[2][1] - L31 * L21) / L22;
  L33 = sqrt(A[2][2] - L31 * L31 - L32 * L32);

  /* Find the solution to the lower triangular system */
  Y1 = B[0] / L11;
  Y2 = (B[1] - L21 * Y1) / L22;
  Y3 = (B[2] - L31 * Y1 - L32 * Y2) / L33;

  /* Find the solution of the upper triangular system */
  x[2] = Y3 / L33;
  x[1] = (Y2 - L32 * x[2]) / L22;
  x[0] = (Y1 - L21 * x[1] - L31 * x[2]) / L11;

  return (x);
}

void
SMtransposeMatrix(double *mat, int rows, int cols, double *mat_T)
{
  int i, j, k;
  k = 0;
  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++)
      mat_T[k++] = mat[i * rows + j];
  }
}

void
SMtransposeMatrix2(double **mat, int rows, int cols, double **mat_T)
{
  int i, j;
  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++)
      mat_T[j][i] = mat[i][j];
  }
}

double
SMcondition3x3(double **A)
{
  double cond;
  double a00, a01, a02, a10, a11, a12, a20, a21, a22;
  double t2, t4, t6, t8, t10, t12, t14, t16;
  double t18, t20, t25, t30, t35, t40, t42, t45, t46;
  double t47, t50, t55, t57, t60, t61, t62, t65, t67;

  a00 = A[0][0];
  a01 = A[0][1];
  a02 = A[0][2];
  a10 = A[1][0];
  a11 = A[1][1];
  a12 = A[1][2];
  a20 = A[2][0];
  a21 = A[2][1];
  a22 = A[2][2];

  t2 = a00 * a00;
  t4 = a01 * a01;
  t6 = a02 * a02;
  t8 = a10 * a10;
  t10 = a11 * a11;
  t12 = a12 * a12;
  t14 = a20 * a20;
  t16 = a21 * a21;
  t18 = a22 * a22;
  t20 = sqrt(t2 + t4 + t6 + t8 + t10 + t12 + t14 + t16 + t18);
  t25 = (a11 * a22 - a12 * a21) * (a11 * a22 - a12 * a21);
  t30 = (a01 * a22 - a02 * a21) * (a01 * a22 - a02 * a21);
  t35 = (a01 * a12 - a02 * a11) * (a01 * a12 - a02 * a11);
  t40 = (a10 * a22 - a12 * a20) * (a10 * a22 - a12 * a20);
  t42 = a02 * a20;
  t45 = (a00 * a22 - t42) * (a00 * a22 - t42);
  t46 = a00 * a12;
  t47 = a02 * a10;
  t50 = (t46 - t47) * (t46 - t47);
  t55 = (a10 * a21 - a11 * a20) * (a10 * a21 - a11 * a20);
  t57 = a01 * a20;
  t60 = (a00 * a21 - t57) * (a00 * a21 - t57);
  t61 = a00 * a11;
  t62 = a01 * a10;
  t65 = (t61 - t62) * (t61 - t62);
  t67 = sqrt(t25 + t30 + t35 + t40 + t45 + t50 + t55 + t60 + t65);
  cond = t20 * t67
      / fabs(
	  t61 * a22 - t46 * a21 - t62 * a22 + t47 * a21 + t57 * a12
	      - t42 * a11);

  return (cond);
}

double
SMdeterminant3x3(double a1[3], double a2[3], double a3[3])
{
  double determinant;

  determinant = a1[0] * (a2[1] * a3[2] - a3[1] * a2[2]);
  determinant -= a2[0] * (a1[1] * a3[2] - a3[1] * a1[2]);
  determinant += a3[0] * (a1[1] * a2[2] - a2[1] * a1[2]);

  return (determinant);
}

void
SMmultiply3x3(double a1[3], double a2[3], double a3[3], double b1[3],
	      double b2[3], double b3[3], double r1[3], double r2[3],
	      double r3[3])
{
  r1[0] = a1[0] * b1[0] + a2[0] * b1[1] + a3[0] * b1[2];
  r1[1] = a1[1] * b1[0] + a2[1] * b1[1] + a3[1] * b1[2];
  r1[2] = a1[2] * b1[0] + a2[2] * b1[1] + a3[2] * b1[2];

  r2[0] = a1[0] * b2[0] + a2[0] * b2[1] + a3[0] * b2[2];
  r2[1] = a1[1] * b2[0] + a2[1] * b2[1] + a3[1] * b2[2];
  r2[2] = a1[2] * b2[0] + a2[2] * b2[1] + a3[2] * b2[2];

  r3[0] = a1[0] * b3[0] + a2[0] * b3[1] + a3[0] * b3[2];
  r3[1] = a1[1] * b3[0] + a2[1] * b3[1] + a3[1] * b3[2];
  r3[2] = a1[2] * b3[0] + a2[2] * b3[1] + a3[2] * b3[2];
}

void
SMtranspose3x3(double a1[3], double a2[3], double a3[3], double b1[3],
	       double b2[3], double b3[3])
{
  b1[0] = a1[0];
  b1[1] = a2[0];
  b1[2] = a3[0];
  b2[0] = a1[1];
  b2[1] = a2[1];
  b2[2] = a3[1];
  b3[0] = a1[2];
  b3[1] = a2[2];
  b3[2] = a3[2];
}

void
SMadjoint3x3(double a1[3], double a2[3], double a3[3], double b1[3],
	     double b2[3], double b3[3])
{
  b1[0] = a2[1] * a3[2] - a3[1] * a2[2];
  b1[1] = a3[1] * a1[2] - a1[1] * a3[2];
  b1[2] = a1[1] * a2[2] - a2[1] * a1[2];

  b2[0] = a3[0] * a2[2] - a2[0] * a3[2];
  b2[1] = a1[0] * a3[2] - a3[0] * a1[2];
  b2[2] = a2[0] * a1[2] - a1[0] * a2[2];

  b3[0] = a2[0] * a3[1] - a3[0] * a2[1];
  b3[1] = a3[0] * a1[1] - a1[0] * a3[1];
  b3[2] = a1[0] * a2[1] - a2[0] * a1[1];
}

double
SMfrobenius_norm_squared3x3(double a1[3], double a2[3], double a3[3])
{
  double norma;

  norma = a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2];
  norma += a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2];
  norma += a3[0] * a3[0] + a3[1] * a3[1] + a3[2] * a3[2];

  return (norma);
}

double
SMfrobenius_norm_squared2x2(double a1[2], double a2[2])
{
  double norma;

  norma = a1[0] * a1[0] + a1[1] * a1[1];
  norma += a2[0] * a2[0] + a2[1] * a2[1];

  return (norma);
}

double
SMfrobenius_norm_squared_adj3x3(double a1[3], double a2[3], double a3[3])
{
  double norm_adja;
  double tmp[3];

  tmp[0] = a1[1] * a2[2] - a1[2] * a2[1];
  tmp[1] = -(a1[0] * a2[1] - a1[2] * a2[0]);
  tmp[2] = a1[0] * a2[1] - a1[1] * a2[0];
  norm_adja = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];

  tmp[0] = a2[1] * a3[2] - a2[2] * a3[1];
  tmp[1] = -(a2[0] * a3[1] - a2[2] * a3[0]);
  tmp[2] = a2[0] * a3[1] - a2[1] * a3[0];
  norm_adja += tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];

  tmp[0] = a3[1] * a1[2] - a3[2] * a1[1];
  tmp[1] = -(a3[0] * a1[1] - a3[2] * a1[0]);
  tmp[2] = a3[0] * a1[1] - a3[1] * a1[0];
  norm_adja += tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];

  return (norm_adja);
}

