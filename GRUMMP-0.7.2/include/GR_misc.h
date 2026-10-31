#ifndef GR_misc
#define GR_misc 1

#include "GR_config.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef misc_implementation
#ifdef __cplusplus
extern "C"
{
#endif
  extern int iMessageStdoutLevel;
  extern int iMessageFileLevel;
  extern FILE* pFMsg;
#ifdef __cplusplus
}
#endif
#endif

#ifndef FILE_NAME_LEN
#define FILE_NAME_LEN 1024
#define FILE_NAME_LEN_STR "1024"
#endif

#ifndef __cplusplus
#undef MAX
#undef MIN

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MIN3(a,b,c) ((a) > (b) ? MIN(b,c) : MIN(a,c))
#define MID3(a,b,c) ((a) > (b) ? ((b) > (c) ? (b) : MIN(a,c)) : ((a) > (c) ? (a) : MIN(b,c)))
#define MAX3(a,b,c) ((a) > (b) ? MAX(a,c) : MAX(b,c))

#else

/* Include std::min and std::max. */
#include <algorithm>
using std::min;
using std::max;
#define MIN3(a,b,c) ((a) > (b) ? min(b,c) : min(a,c))
#define MID3(a,b,c) ((a) > (b) ? ((b) > (c) ? (b) : min(a,c)) : ((a) > (c) ? (a) : min(b,c)))
#define MAX3(a,b,c) ((a) > (b) ? max(a,c) : max(b,c))

#endif

#define XOR(a,b) ((a) ? (!b) : (b))

#undef XDIR
#undef YDIR
#undef ZDIR

#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define iInvalidBC INT_MAX
#define iDefaultBC INT_MAX - 1
/* Region numbers are unsigned N-bit integers when attached to cells. */
/* Legal region numbers run from 1 to 2^N - 1, with 0 marking cells with */
/* unknown or invalid region and 2^N - 1 cells outside the domain. */
#define iInvalidRegion 0
#define iRegionBits 7
#define iMaxRegionLabel (1 << iRegionBits)
#define iOutsideRegion iMaxRegionLabel - 1
#define iMaxRegion iMaxRegionLabel - 2
#define iInsideRegion 1
#define iDefaultRegion 1

/* Define a macro so that gcc can use function and variable attributes
 without messing up other compilers. */
#ifdef __GNUC__
#define ATTRIBUTE(a) __attribute__(a)
#else
#define ATTRIBUTE(a)
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define MSG_GLOBAL 0
#define MSG_MANAGER 1
#define MSG_WARNING 1
#define MSG_SERVICE 2
#define MSG_DETAIL 3
#define MSG_DEBUG 4
#define MSG_DEBUG_INACTIVE 5
  /* Writes messages to stdout and to a log file.  The global variables
   * iMessageStdoutLevel and iMessageFileLevel control what level messages
   * will appear where (everything at or below those variable values.
   *
   * In general, messages in production code should be tagged as:
   * 0: Messages from executables, abort messages, other error messages.
   * 1: Algorithmic milestones from manager-level code
   *    Examples: finished a pass of swapping, inserted X0000 vertices,
   *    etc.
   * 2: Algorithmic milestones from lower-level code
   * 	  Examples:
   * 3: Operational / decision details from lower-level code
   * 4: Active debugging messages
   * 5: Inactive debugging message: Never appear anywhere
   */
  void
  logMessage(const int i, const char *acFormat, ...)
      ATTRIBUTE ((format (printf, 2, 3)));

  void
  openMessageFile(const char strInFileName[]);
  void
  closeMessageFile(void);
  /* The following function is used for sorting Lists.  It has to be declared
   here instead of in List.h for compilers that instantiate templates at
   link time. */
  bool
  qFuzzyPerp3D(const double adA[3], const double adB[3]);
  int
  iFuzzyComp(const double dA, const double dB);
  void
  vGetLineOrAbort(char acBuffer[], const int iBufSize, FILE *pInFile);
  void
  vSkipCommentLines(FILE *pInFile);
  void
  vGetDoubleFromBuffer(char **ppcBuf, double * const pdData,
		       const char * const strErrorString,
		       const char * const strContextString);
  void
  vGetIntFromBuffer(char **ppcBuf, int * const piData,
		    const char * const strErrorString,
		    const char * const strContextString);
  void
  solve3By3(const double adRow1_in[3], const double adRow2_in[3],
	    const double adRow3_in[3], double dRHS1, double dRHS2, double dRHS3,
	    double dResult[3]);
  void
  solveTriDiag(double* const a2dLHS[3], double* const adRHS, const int iSize);
  void
  solveQuadratic(double dA, double dB, double dC, int* iNRoots, double* dRoot1,
		 double* dRoot2);
  void
  solveCubic(double dA, double dB, double dC, double dD, int* iNRoots,
	     double adRoots[3]);
  double
  calcDet4By4(double a2dMat4[4][4]);
  void
  calcUnitNormal(const double adA[3], const double adB[3], const double adC[3],
		 double adRes[3]);
  void
  calcNormal(const double adA[3], const double adB[3], const double adC[3],
	     double adRes[3]);
  void
  makeFileName(char *strNewName, const char* strFormat, const char* strBaseName,
	       const char* strCaller);
  void
  GRUMMPInit(const char strExecName[]);
  void
  GRUMMPSignOff(void);

  void
  vGRUMMP_2D_Lib(void);
  void
  vGRUMMP_3D_Lib(void);
  void
  vGRUMMP_Surf_Lib(void);

  /* A range-safe version of the arccosine, to prevent numerical
   difficulties when dArg = 1 + machine epsilon. */
  double
  GR_acos(const double dArg);

#ifdef __cplusplus
}
#endif

#define vFatalError(message, function) \
do { \
  logMessage(0, "Fatal error in %s\n%s\n", function, message); \
  abort(); \
} while (0)

#define vPrintErrorLine(line, file) \
do { \
  logMessage(0, "Error on line %d of %s\n", line, file);\
} while (0)

#define vFoundBug(arg_message) \
do { \
  logMessage(0, "\n"); \
  logMessage(0, "Bug found in GRUMMP version %s\n", GRUMMP_VERSION); \
  logMessage(0, "  Internal error: %s\n", arg_message); \
  logMessage(0, "  In file %s, line %d.\n", __FILE__, __LINE__); \
  logMessage(0, "Please report this bug via email to GRUMMP-bug@mech.ubc.ca\n"); \
  logMessage(0, "Include this entire message, the command line used, and if\n"); \
  logMessage(0, "possible the entire input file used.\n\n"); \
  abort(); \
} while (0)

#define vWarning(message) logMessage(0, "Warning: %s\n", message)

#define vPrintVersion() do {logMessage(0, "Linked to GRUMMP libraries version %s.\n\n", GRUMMP_VERSION);} while (0)

#define vInfo(a) do {logMessage(0, "Running executable %s on %s\n", a, HOST_OS); \
vPrintVersion();} while (0)

#endif

