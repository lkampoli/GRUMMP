#ifndef VC__GRUMMMP_00300_CONFIG_H__INCLUDED
#define VC__GRUMMMP_00300_CONFIG_H__INCLUDED

#ifdef WIN32
#define WIN_VC32

#include <io.h>

/* symbol > 255 chars (DEBUG std::templates) */
#pragma warning(disable:4786) 

/* To enable logging, define "SUMAA_LOG=1" in VC IDE */

/* Define if you want to use polygonal/polyhedral input files */
/* instead of new-style boundary descriptions.  The old-style */
/* input files will no longer be supported at some near-future */
/* release, although they may not be removed for some time yet. */
#undef USE_DEPRECATED

/* bool type not known to C compiler */
#ifndef __cplusplus
# define bool  unsigned int
# define true  1
# define false 0
#endif

/* Define to unsigned if sys/types.h doesn't define it */
#define mode_t unsigned  

#define HAVE_SNPRINTF 1
#define snprintf _snprintf

#define GRUMMP_AT_EXIT(a) atexit(a)

/* Define to enable assertions */
/* NDEBUG managed by VC++ IDE */

#ifdef NDEBUG
#define MESSAGE_LEVEL_STDOUT 1
#define MESSAGE_LEVEL_FILE 1
#else
#define MESSAGE_LEVEL_STDOUT 1
#define MESSAGE_LEVEL_FILE 1
#endif

/* The following is fine for IEEE double precision, or anything else */
/* with at least ten bits of binary exponent. */
#define LARGE_DBL 1.e+300

/* drand48() on Win32 platforms */
/* See: WIN32_rand48.h, Win32_drand48.c */

/* "drand48" returns a pseudo-random number in [0.0,1.0). */
/*  Before the first call to "drand48", you should call   */
/* "srand48" or "seed48" to initialize random "seed".     */

#ifdef __cplusplus
extern "C" {
#endif
  double drand48(void);
#ifdef __cplusplus
}
#endif

#define GRUMMP_VERSION "0.3.0"
#define GRUMMP_SERIAL 00300
#define HOST_OS "Win32"

/* Be sure these two always get included along with config info. */
#include "GR_assert.h"

#ifdef __cplusplus
#include "GR_Classes.h"
#endif

/* Add missing constants */
#ifndef M_PI
#define M_PI       3.14159265358979323846   /* pi */
#define M_PI_2     1.57079632679489661923   /* pi/2 */
#define M_LOG2E    1.4426950408889634074    /* log 2e */
#define M_SQRT2    1.41421356237309504880   /* sqrt(2) */
#define M_SQRT1_2  0.70710678118654752440   /* sqrt(1/2) */
#define M_LN2      0.69314718055994530942   /* ln(2) */
#define M_LN10     2.30258509299404568402   /* ln(10) */
#endif

#include <float.h>          /* DBL_MAX,  _finite*/
#define snprintf _snprintf
#define finite   _finite

int getopt(int argc, char *argv[], char *opstring);
extern char *optarg; /* defined in base/GetOpts.cxx */

#endif /* WIN32 */
#endif /* VC__GRUMMMP_00300_CONFIG_H__INCLUDED */
