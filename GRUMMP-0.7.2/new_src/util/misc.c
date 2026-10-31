#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "GR_AdaptPred.h"
#include "GR_Vec.h"

FILE *pFMsg;

#include "GR_config.h"
#define misc_implementation
#include "GR_misc.h"
#undef misc_implementation

char strExecName[1024];

#ifdef WIN32
#define  RAND48_SEED_0  (0x330e)
#define  RAND48_SEED_1  (0xabcd)
#define  RAND48_SEED_2  (0x1234)
#define  RAND48_MULT_0  (0xe66d)
#define  RAND48_MULT_1  (0xdeec)
#define  RAND48_MULT_2  (0x0005)
#define  RAND48_ADD     (0x000b)

unsigned short __rand48_seed[3] = {
  RAND48_SEED_0,
  RAND48_SEED_1,
  RAND48_SEED_2
};

unsigned short __rand48_mult[3] = {
  RAND48_MULT_0,
  RAND48_MULT_1,
  RAND48_MULT_2
};
unsigned short __rand48_add = RAND48_ADD;

unsigned short *seed48(unsigned short xseed[3])
{
  static unsigned short sseed[3];

  sseed[0] = __rand48_seed[0];
  sseed[1] = __rand48_seed[1];
  sseed[2] = __rand48_seed[2];
  __rand48_seed[0] = xseed[0];
  __rand48_seed[1] = xseed[1];
  __rand48_seed[2] = xseed[2];
  __rand48_mult[0] = RAND48_MULT_0;
  __rand48_mult[1] = RAND48_MULT_1;
  __rand48_mult[2] = RAND48_MULT_2;
  __rand48_add = RAND48_ADD;
  return sseed;
}

void srand48(long seed)
{
  __rand48_seed[0] = RAND48_SEED_0;
  __rand48_seed[1] = (unsigned short) seed;
  __rand48_seed[2] = (unsigned short) (seed >> 16);
  __rand48_mult[0] = RAND48_MULT_0;
  __rand48_mult[1] = RAND48_MULT_1;
  __rand48_mult[2] = RAND48_MULT_2;
  __rand48_add = RAND48_ADD;
}

static void __dorand48(unsigned short xseed[3])
{
  unsigned long accu;
  unsigned short temp[2];

  accu = (unsigned long) __rand48_mult[0] * (unsigned long) xseed[0] + (unsigned long) __rand48_add;
  temp[0] = (unsigned short) accu; /* lower 16 bits */

  accu >>= sizeof(unsigned short) * 8;
  accu += (unsigned long) __rand48_mult[0] * (unsigned long) xseed[1] + (unsigned long) __rand48_mult[1] * (unsigned long) xseed[0];
  temp[1] = (unsigned short) accu; /* middle 16 bits */

  accu >>= sizeof(unsigned short) * 8;
  accu += __rand48_mult[0] * xseed[2] + __rand48_mult[1] * xseed[1] + __rand48_mult[2] * xseed[0];

  xseed[0] = temp[0];
  xseed[1] = temp[1];
  xseed[2] = (unsigned short) accu;
}

static double erand48(unsigned short xseed[3])
{
  __dorand48(xseed);
  return ldexp((double) xseed[0], -48) +
  ldexp((double) xseed[1], -32) +
  ldexp((double) xseed[2], -16);
}

double drand48(void)
{
  return erand48(__rand48_seed);
}

#endif /* WIN32 */

void
GRUMMPInit(const char strExecNameIn[])
{
  GRUMMP_AT_EXIT(GRUMMPSignOff);
  strExecName[1023] = 0x0;
  strncpy(strExecName, strExecNameIn, (size_t) (1023));
  /* Print banner */
  vInfo(strExecName);
  /* This call includes registration of smoothing events. */
#ifdef SUMAA_LOG
  SUMAA_LOG_EVENT_BEGIN(TOTAL_TIME);
#endif
  /* Print banner */
  vInfo(strExecName);
  /* Initialization for adaptive predicates */
  exactinit();
  srand48((long int) (12345678));
}

void
GRUMMPSignOff(void)
{
  logMessage(
      0,
      "\nExecutable %s (linked to GRUMMP version %s) exiting successfully.\n\n",
      strExecName, GRUMMP_VERSION);
  closeMessageFile();
}

/* These levels are #define'd in config.h */
int iMessageStdoutLevel = MESSAGE_LEVEL_STDOUT;
int iMessageFileLevel = MESSAGE_LEVEL_FILE;

void
openMessageFile(const char strBaseFileName[])
{
  char strMsgFileName[FILE_NAME_LEN];

  makeFileName(strMsgFileName, "%s.msg", strBaseFileName, "vOpenMessageFile");
  pFMsg = fopen(strMsgFileName, "w");
  if (NULL == pFMsg)
    vFatalError("Couldn't open message file for writing.", "vOpenMessageFile");
}

void
logMessage(const int i, const char *acFormat, ...)
{
  /*   SUMAA_LOG_EVENT_BEGIN(PROBE3); */

  va_list VArgList1, VArgList2;
  va_start(VArgList1, acFormat);
  va_copy(VArgList2, VArgList1);

  if (i <= iMessageStdoutLevel) {
    vprintf(acFormat, VArgList1);
    fflush(stdout);
  }
  if (i <= iMessageFileLevel && pFMsg) {
    vfprintf(pFMsg, acFormat, VArgList2);
    fflush(pFMsg);
  }

  va_end(VArgList1);
  va_end(VArgList2);
  /*   SUMAA_LOG_EVENT_END(PROBE3); */
}

void
closeMessageFile(void)
{
  if (pFMsg)
    fclose(pFMsg);
}

bool
qFuzzyPerp3D(const double adAIn[3], const double adBIn[3])
{
  double adCross[3], dMag;
  double adA[3], adB[3];
  adA[0] = adAIn[0];
  adA[1] = adAIn[1];
  adA[2] = adAIn[2];
  adB[0] = adBIn[0];
  adB[1] = adBIn[1];
  adB[2] = adBIn[2];
  NORMALIZE3D(adA);
  NORMALIZE3D(adB);
  vCROSS3D(adA, adB, adCross);
  dMag = dMAG3D(adCross);
  return (iFuzzyComp(dMag, 1.0) == 0);
}

int
iFuzzyComp(const double dA, const double dB)
{
  double dSum = fabs(dA) + fabs(dB);
  double dDiff = dA - dB;
  const double dTol = 5.e-9;
  const double dFloor = 4.e-5;

  dSum = (dSum > dFloor) ? dSum : dFloor;
  if (dDiff > dTol * dSum)
    return (1);
  else if (dDiff < -dTol * dSum)
    return (-1);
  else
    return (0);
}

/* Not a stunningly efficient approach for this, probably, but slow
 processing of file input is a ridiculously minor efficiency issue for
 GRUMMP; all the time is spent on mesh manipulation. */
/* IRIX cc won't take a const int iBufSize, so here's the hack: */
#define iBufSize 10240
void
vSkipCommentLines(FILE *pFInput)
{
  long lFPos;
  char acBuffer[iBufSize];
  do {
    lFPos = ftell(pFInput);
    /* Skip comments and blank lines.  EOF is always unexpected,
     because we always know whether to read another line. */
    if (NULL == fgets(acBuffer, iBufSize, pFInput))
      vFatalError("Unexpected end of file", "vSkipCommentLines");
  }
  while (acBuffer[0] == '!' || acBuffer[0] == '#' || acBuffer[0] == '\n'
      || acBuffer[0] == '\000');
  fseek(pFInput, lFPos, 0);
}
/* The other end of the iBufSize hack. */
#undef iBufSize

void
vGetLineOrAbort(char acBuffer[], const int iBufSize, FILE *pFInput)
{
  int i;
  do {
    do {
      /* Skip comments and blank lines.  EOF is always unexpected,
       because we always know whether to read another line. */
      void *res = fgets(acBuffer, iBufSize, pFInput);
      if (NULL == res)
	vFatalError("Unexpected end of file", "vGetLineOrAbort");
    }
    while (acBuffer[0] == '!' || acBuffer[0] == '#' || acBuffer[0] == '\n'
	|| acBuffer[0] == '\000');
    /* Strip whitespace off the end of the string by replacing
     whitespace with nulls. */
    for (i = (int) (strlen(acBuffer)) - 1;
	i >= 0 && isspace((int) (acBuffer[i])); i--)
      acBuffer[i] = '\000';
  }
  while (strlen(acBuffer) == 0); /* If the resulting string is empty,
   get another line. */
  /* Convert the string to lower case */
  for (i = 0; acBuffer[i] != '\000'; i++) {
    acBuffer[i] = (char) (tolower((int) (acBuffer[i])));
  }
  return;
}

void
vGetDoubleFromBuffer(char **ppcBuf, double * const pdData,
		     const char * const strErrorString,
		     const char * const strContextString)
{
  char *pcBufInit;
  pcBufInit = *ppcBuf;
  *pdData = strtod(pcBufInit, ppcBuf);

  if (pcBufInit == *ppcBuf)
    vFatalError(strErrorString, strContextString);
}

void
vGetIntFromBuffer(char **ppcBuf, int * const piData,
		  const char * const strErrorString,
		  const char * const strContextString)
{
  char *pcBufInit = *ppcBuf;
  *piData = (int) strtol(pcBufInit, ppcBuf, 10);
  if (pcBufInit == *ppcBuf)
    vFatalError(strErrorString, strContextString);
}

void
solve3By3(const double adRow1_in[3], const double adRow2_in[3],
	  const double adRow3_in[3], double dRHS1, double dRHS2, double dRHS3,
	  double adRes[3])
{
  double adRow1[3];
  double adRow2[3];
  double adRow3[3];

  adRow1[0] = adRow1_in[0];
  adRow1[1] = adRow1_in[1];
  adRow1[2] = adRow1_in[2];
  adRow2[0] = adRow2_in[0];
  adRow2[1] = adRow2_in[1];
  adRow2[2] = adRow2_in[2];
  adRow3[0] = adRow3_in[0];
  adRow3[1] = adRow3_in[1];
  adRow3[2] = adRow3_in[2];

  /* Pivot first column */
  if (fabs(adRow1[0]) >= fabs(adRow2[0])
      && fabs(adRow1[0]) >= fabs(adRow3[0])) {
  }
  else if (fabs(adRow2[0]) >= fabs(adRow3[0])) {
    double dTemp = dRHS1;
    dRHS1 = dRHS2;
    dRHS2 = dTemp;
    vSWAP3D(adRow1, adRow2);
  }
  else {
    double dTemp = dRHS1;
    dRHS1 = dRHS3;
    dRHS3 = dTemp;
    vSWAP3D(adRow1, adRow3);
  }

  /* Eliminate first column */
  dRHS2 -= dRHS1 * adRow2[0] / adRow1[0];
  dRHS3 -= dRHS1 * adRow3[0] / adRow1[0];
  adRow2[1] = adRow2[1] - adRow1[1] * (adRow2[0] / adRow1[0]);
  adRow2[2] = adRow2[2] - adRow1[2] * (adRow2[0] / adRow1[0]);
  adRow3[1] = adRow3[1] - adRow1[1] * (adRow3[0] / adRow1[0]);
  adRow3[2] = adRow3[2] - adRow1[2] * (adRow3[0] / adRow1[0]);

  /* Pivot second column */
  if (fabs(adRow2[1]) >= fabs(adRow3[1])) {
  }
  else {
    double dTemp = dRHS2;
    dRHS2 = dRHS3;
    dRHS3 = dTemp;
    vSWAP3D(adRow2, adRow3);
  }

  /* Eliminate second column */
  dRHS3 -= dRHS2 * adRow3[1] / adRow2[1];
  adRow3[2] = adRow3[2] - adRow2[2] * (adRow3[1] / adRow2[1]);

  /* Solve for dRHS3 and back-substitute */
  adRes[2] = dRHS3 /= adRow3[2];
  dRHS2 -= adRow2[2] * dRHS3;
  dRHS1 -= adRow1[2] * dRHS3;

  /* Solve for dRHS2 and back-substitute */
  adRes[1] = dRHS2 /= adRow2[1];
  dRHS1 -= adRow1[1] * dRHS2;

  /* Solve for dRHS1 */
  adRes[0] = dRHS1 /= adRow1[0];
}

void
solveTriDiag(double* const a2dLHS[3], double* const adRHS, const int iSize)
{
  // iSize should be the actual number of rows in the matrix.  This means
  // that the rows are numbered from 0 to iSize - 1.
  int i;
  a2dLHS[0][0] = a2dLHS[iSize - 1][2] = 0;

  /* Forward elimination */
  for (i = 0; i < iSize - 1; i++) {
    a2dLHS[i][2] /= a2dLHS[i][1];
    adRHS[i] /= a2dLHS[i][1];
    a2dLHS[i + 1][1] -= (a2dLHS[i][2] * a2dLHS[i + 1][0]);
    adRHS[i + 1] -= (a2dLHS[i + 1][0] * adRHS[i]);
  }

  /* Last line of elimination */
  adRHS[iSize - 1] /= a2dLHS[iSize - 1][1];

  /* Back-substitution */
  for (i = iSize - 2; i >= 0; i--) {
    adRHS[i] -= adRHS[i + 1] * a2dLHS[i][2];
  }

}

void
solveQuadratic(double dA, double dB, double dC, int* iNRoots, double* dRoot1,
	       double* dRoot2)
{

  //Finds the (possibly two) REAL roots of an equation
  //of the form dA * x^2 + dB * x + dC = 0. The number
  //of real roots is returned in iNRoots, and the two
  //roots in dRoot1 and dRoot2.

  {
    //Added this scaling function so that the biggest coeff is 1.
    //The iFuzzyComp calls are more robust this way.
    double dScale = MAX3(fabs(dA), fabs(dB), fabs(dC));
    dA /= dScale;
    dB /= dScale;
    dC /= dScale;
  }
  if (iFuzzyComp(dA, 0.) == 0) {
    if (iFuzzyComp(dB, 0.) == 0) {
      *iNRoots = 0;
      *dRoot1 = - LARGE_DBL;
      *dRoot2 = - LARGE_DBL;
      return;
    }
    else {
      *iNRoots = 1;
      *dRoot1 = -(dC / dB);
      *dRoot2 = - LARGE_DBL;
      return;
    }
  }
  else {
    double dRad = (dB * dB) - (4. * dA * dC);
    if (iFuzzyComp(dRad, 0.) == -1) {
      //No REAL roots.
      *iNRoots = 0;
      *dRoot1 = - LARGE_DBL;
      *dRoot2 = - LARGE_DBL;
      return;
    }
    else if (iFuzzyComp(dRad, 0.) == 0) {
      *iNRoots = 1;
      *dRoot1 = -(dB / (2. * dA));
      *dRoot2 = - LARGE_DBL;
      return;
    }
    else {
      assert(iFuzzyComp(dRad, 0.) == 1);
      *iNRoots = 2;
      double dRadSqrt = sqrt(dRad);
      *dRoot1 = (-dB + dRadSqrt) / (2. * dA);
      *dRoot2 = (-dB - dRadSqrt) / (2. * dA);
      return;
    }
  }

}

void
solveCubic(double dA, double dB, double dC, double dD, int* iNRoots,
	   double adRoots[3])
{

  {
    //Added this scaling function so that the biggest coeff is 1.
    //The iFuzzyComp calls are more robust this way.
    double dScale = MAX(fabs(dA), MAX(fabs(dB), MAX(fabs(dC), fabs(dD))));
    dA /= dScale;
    dB /= dScale;
    dC /= dScale;
    dD /= dScale;
  }

  if (iFuzzyComp(dA, 0.) == 0) {
    if (iFuzzyComp(dB, 0.) == 0) {
      if (iFuzzyComp(dC, 0.) == 0) {
	*iNRoots = 0;
	adRoots[0] = adRoots[1] = adRoots[2] = - LARGE_DBL;
	return;
      }
      else {
	*iNRoots = 1;
	adRoots[0] = -dD / dC;
	adRoots[1] = adRoots[2] = - LARGE_DBL;
	return;
      }
    }
    else {
      solveQuadratic(dB, dC, dD, iNRoots, &adRoots[0], &adRoots[1]);
      adRoots[2] = -LARGE_DBL;
      return;
    }
  }
  else {
    double dF = (dC / dA) - (pow(dB / dA, 2.) / 3.);
    double dG = (2. * pow(dB, 3.)) / (27. * pow(dA, 3.))
	- (dB * dC) / (3. * pow(dA, 2.)) + (dD / dA);
    double dH = (0.25 * pow(dG, 2.)) + (pow(dF, 3.) / 27.);

    switch (iFuzzyComp(dH, 0.))
      {
      case -1:
	{ //Three different real roots
	  double dI = sqrt(0.25 * pow(dG, 2.) - dH);
	  double dJ;
	  if (dI < 0.) {
	    dJ = -pow(-dI, 1. / 3.);
	  }
	  else {
	    dJ = pow(dI, 1. / 3.);
	  }

	  double dK = GR_acos(-0.5 * dG / dI);
	  double dM = cos(dK / 3.);
	  double dN = sqrt(3.) * sin(dK / 3.);
	  double dP = -dB / (3. * dA);

	  *iNRoots = 3;
	  adRoots[0] = 2 * dJ * cos(dK / 3.) - dB / (3. * dA);
	  adRoots[1] = -dJ * (dM + dN) + dP;
	  adRoots[2] = -dJ * (dM - dN) + dP;
	  break;
	}
      case 0:
	{ //All three roots are real and equal
	  *iNRoots = 3;
	  double dQ = dD / dA;
	  if (dQ < 0.) {
	    adRoots[0] = adRoots[1] = adRoots[2] = pow(-dQ, 1. / 3.);
	  }
	  else {
	    adRoots[0] = adRoots[1] = adRoots[2] = -pow(dQ, 1. / 3.);
	  }
	  break;
	}
      case 1:
	{ //One real root, two complex conjugate roots.
	  double dR = -0.5 * dG + sqrt(dH);
	  double dS = -0.5 * dG - sqrt(dH);

	  double dT, dU;
	  if (dR < 0.) {
	    dT = -pow(-dR, 1. / 3.);
	  }
	  else {
	    dT = pow(dR, 1. / 3.);
	  }
	  if (dS < 0.) {
	    dU = -pow(-dS, 1. / 3.);
	  }
	  else {
	    dU = pow(dS, 1. / 3.);
	  }

	  *iNRoots = 1;
	  adRoots[0] = dT + dU - dB / (3. * dA);
	  adRoots[1] = adRoots[2] = - LARGE_DBL;
	  break;
	}
      default:
	assert(0); //should not get here.
	break;
      }

  }

}

double
calcDet4By4(double a2dMat4[4][4])
{
  return (a2dMat4[0][0]
      * (a2dMat4[1][1] * a2dMat4[2][2] * a2dMat4[3][3]
	  + a2dMat4[1][2] * a2dMat4[2][3] * a2dMat4[3][1]
	  + a2dMat4[1][3] * a2dMat4[2][1] * a2dMat4[3][2]
	  - a2dMat4[3][1] * a2dMat4[2][2] * a2dMat4[1][3]
	  - a2dMat4[3][2] * a2dMat4[2][3] * a2dMat4[1][1]
	  - a2dMat4[3][3] * a2dMat4[2][1] * a2dMat4[1][2])
      - a2dMat4[0][1]
	  * (a2dMat4[1][0] * a2dMat4[2][2] * a2dMat4[3][3]
	      + a2dMat4[1][2] * a2dMat4[2][3] * a2dMat4[3][0]
	      + a2dMat4[1][3] * a2dMat4[2][0] * a2dMat4[3][2]
	      - a2dMat4[3][0] * a2dMat4[2][2] * a2dMat4[1][3]
	      - a2dMat4[3][2] * a2dMat4[2][3] * a2dMat4[1][0]
	      - a2dMat4[3][3] * a2dMat4[2][0] * a2dMat4[1][2])
      + a2dMat4[0][2]
	  * (a2dMat4[1][0] * a2dMat4[2][1] * a2dMat4[3][3]
	      + a2dMat4[1][1] * a2dMat4[2][3] * a2dMat4[3][0]
	      + a2dMat4[1][3] * a2dMat4[2][0] * a2dMat4[3][1]
	      - a2dMat4[3][0] * a2dMat4[2][1] * a2dMat4[1][3]
	      - a2dMat4[3][1] * a2dMat4[2][3] * a2dMat4[1][0]
	      - a2dMat4[3][3] * a2dMat4[2][0] * a2dMat4[1][1])
      - a2dMat4[0][3]
	  * (a2dMat4[1][0] * a2dMat4[2][1] * a2dMat4[3][2]
	      + a2dMat4[1][1] * a2dMat4[2][2] * a2dMat4[3][0]
	      + a2dMat4[1][2] * a2dMat4[2][0] * a2dMat4[3][1]
	      - a2dMat4[3][0] * a2dMat4[2][1] * a2dMat4[1][2]
	      - a2dMat4[3][1] * a2dMat4[2][2] * a2dMat4[1][0]
	      - a2dMat4[3][2] * a2dMat4[2][0] * a2dMat4[1][1]));
}

/* Compute the vector normal for a set of three verts */
void
calcNormal(const double adA[3], const double adB[3], const double adC[3],
	   double adRes[3])
{
  double adTempB[3];
  double adTempC[3];
  adTempB[0] = adB[0] - adA[0];
  adTempB[1] = adB[1] - adA[1];
  adTempB[2] = adB[2] - adA[2];

  adTempC[0] = adC[0] - adA[0];
  adTempC[1] = adC[1] - adA[1];
  adTempC[2] = adC[2] - adA[2];

  vCROSS3D(adTempB, adTempC, adRes);
}

/* Compute the vector normal for a set of three verts */
void
calcUnitNormal(const double adA[3], const double adB[3], const double adC[3],
	       double adRes[3])
{
  calcNormal(adA, adB, adC, adRes);
  NORMALIZE3D(adRes);
}

void
makeFileName(char *strNewName, const char* strFormat, const char* strBaseName,
	     const char* strCaller)
{
  int iLen = snprintf(strNewName, (FILE_NAME_LEN), strFormat, strBaseName);
  if (-1 == iLen || FILE_NAME_LEN == iLen) {
    logMessage(-1, "Problem caused by call from %s to makeFileName()\n",
	       strCaller);
    vFatalError("Quality file name too long", "makeFileName()");
  }
}

#ifndef HAVE_SNPRINTF
int snprintf(char *s, size_t len, const char *format, ...) {
  char *s_tmp;
  va_list args;
  size_t result, tmp_len;

  tmp_len = 10*len;
  if (tmp_len < 1024)
  tmp_len = 1024;
  s_tmp = (char*) calloc(len*10, sizeof(char));
  va_start(args, format);
  (void) vsprintf(s_tmp, format, args);
  va_end(args);
  result = strlen(s_tmp);
  if (result > tmp_len)
  vFatalError("Buffer overrun", "internal snprintf");
  if (result < len)
  (void) strncpy(s, s_tmp, (size_t) result+1);
  if (result >= len) {
    result = len;
    vWarning("String too long in internal snprintf; truncated!");
    (void) strncpy(s, s_tmp, (size_t) len);
    s[len-1] = '\0';
  }
  free(s_tmp);
  return result;
}
#endif

double
GR_acos(const double dArg)
{
  if (dArg > 1.0)
    return 0;
  else if (dArg < -1.0)
    return M_PI;
  else
    return acos(dArg);
}
