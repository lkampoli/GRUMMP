/* include/GR_config.h.  Generated from GR_config.h.in by configure.  */
#ifndef CONFIG_H_GRUMMP
#define CONFIG_H_GRUMMP

/* include/GR_config.h.in.  Generated from configure.in by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Unique, monotonically increasing serial number for this version */
#define GRUMMP_SERIAL 00702

/* Version number as a string. */
#define GRUMMP_VERSION "0.7.2"

/* Define to 1 if compiling with MPI */
#define HAVE_MPI 1

/* Define to 1 if compiling with OpenMP */
#define HAVE_OPENMP 1

/* Define to 1 if you have the Mesquite headers available. */
/* #undef HAVE_MESQUITE */

/* Define to 1 if you have the Mesquite library available. */
/* #undef HAVE_MESQUITE_LIB */

/* Define to 1 if you have the Mesquite iMesh library available. */
/* #undef HAVE_MSQIMESH_LIB */

/* Define to 1 if you have the <unordered_map> header file. */
/* #undef HAVE_CXX0X */

/* Define to 1 if you have the `atexit' function. */
#define HAVE_ATEXIT 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Define to 1 if you have the `drand48' function. */
#define HAVE_DRAND48 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `dl' library (-ldl). */
#define HAVE_LIBDL 1

/* Define to 1 if you have the `fl' library (-lfl). */
/* #undef HAVE_LIBFL */

/* Define to 1 if you have the `l' library (-ll). */
/* #undef HAVE_LIBL */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `TKernel' library (-lTKernel). */
/* #undef HAVE_LIBTKERNEL */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Is OpenCascade present? */
/* #undef HAVE_OCC */

/* Define to 1 if you have the `on_exit' function. */
#define HAVE_ON_EXIT 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strcasecmp' function. */
#define HAVE_STRCASECMP 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strtod' function. */
#define HAVE_STRTOD 1

/* Define to 1 if you have the `strtol' function. */
#define HAVE_STRTOL 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/timeb.h> header file. */
#define HAVE_SYS_TIMEB_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <tr1/unordered_map> header file. */
#define HAVE_TR1 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <values.h> header file. */
#define HAVE_VALUES_H 1

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* Does the C compiler support a bool type? */
/* #undef HAVE_bool_c */

/* What machine type is this? */
#define HOST_OS "LINUX_GNU"

/* Build with support for ITAPS */
#ifndef ITAPS
#define ITAPS 1
#endif

/* Compile without assertions and extra checks in this is set. */
#define NDEBUG 1

/* Turn debugging output off in the smoothing code. */
#define SM_DEBUG_LEVEL 0

/* Output statistics for smoothing results? */
#define SM_STATS 1 

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
/* #undef YYTEXT_POINTER */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the type of a signed integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef int32_t */

/* Define to the type of an unsigned integer type of width exactly 32 bits if
   such a type exists and the standard includes do not define it. */
/* #undef uint32_t */

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

#ifndef __cplusplus
/* Define to empty if the const keyword is broken in the C compiler */
/* #undef const */
#endif

/* Define if the bool type is known to your compiler */
/* #undef HAVE_bool_c */

/* Define if vprintf exists */
#define HAVE_VPRINTF 1

/* #undef SUNOS4 */
/* #undef SOLARIS2 */
#define LINUX_GNU 1
/* #undef LINUX */
/* #undef HPUX10 */
/* #undef OSF4 */

#ifdef NDEBUG
#define MESSAGE_LEVEL_STDOUT 1
#define MESSAGE_LEVEL_FILE 2
#else
#define MESSAGE_LEVEL_STDOUT 1
#define MESSAGE_LEVEL_FILE 3
#endif

/* The following is fine for IEEE double precision, or anything else */
/* with at least ten bits of binary exponent. */
#define LARGE_DBL 1.e+300

#ifdef SUNOS4
/* Sun OS 4.1.x has some missing declarations */
#include <stdio.h>
#ifdef __cplusplus
extern "C" {
#endif
int fflush( FILE *stream);
int fclose( FILE *stream);
int tolower (int c);
long int strtol(const char *nptr, char **endptr, int base);
int printf( const char *format, ...);
int fprintf( FILE *stream, const char *format, ...);
#ifdef __cplusplus
}
#else
#include <varargs.h>
int vprintf(const char *format, va_list ap);
int vfprintf(FILE *stream, const char *format, va_list ap);
char *vsprintf(char *s, const char *format, va_list ap);
#endif
#endif

#if ((!defined(HAVE_DRAND48)) || (defined(SUNOS)))
#ifdef __cplusplus
extern "C" {
#endif
double drand48(void);
void srand48(long seed);
#ifdef __cplusplus
}
#endif
#endif

#ifndef __cplusplus

#ifndef HAVE_bool_c
#define bool unsigned int
#define true 1
#define false 0
#endif

#endif

#if ((!defined(HAVE_SNPRINTF)) || (defined(_AIX) && !defined(__GNUG__)))
#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
int snprintf ( char *str, size_t n, const char *format, ... );
#ifdef __cplusplus
}
#endif
#endif

#if defined(HAVE_ATEXIT)
/* atexit() is defined */
#define GRUMMP_AT_EXIT(a) atexit(a)
#else
#if defined(HAVE_ON_EXIT)
/* atexit() is not defined but on_exit() is.  Sun 4, for example, has */
/* no prototype for atexit(), although the function appears to exist. */
#define GRUMMP_AT_EXIT(a) on_exit(a, (void*)(0))
#else
/* neither atexit() nor on_exit() is defined */
#define GRUMMP_AT_EXIT(a) 
#endif /* HAVE_ON_EXIT */
#endif /* HAVE_ATEXIT */

#if ((defined(HPUX10)) && (!defined(__GNUG__)))
/* HPUX10 native compiler pukes on inlines with loops instead of
   converting them to non-inline functions. */
#define BROKEN_INLINE 
#endif

#ifdef WIN32
#ifdef __cplusplus
extern "C" {
#endif
int getopt(int argc, char * const argv[],
                  const char *optstring);
#ifdef __cplusplus
}
#endif
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
typedef uint32_t GR_index_t;
typedef int32_t  GR_sindex_t;

#ifdef __cplusplus
#ifdef HAVE_TR1
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#define GR_map std::tr1::unordered_map
#define GR_set std::tr1::unordered_set
#else
#ifdef HAVE_CXX0X 
#include <unordered_set>
#include <unordered_map>
#define GR_map std::unordered_map
#define GR_set std::unordered_set
#else
#include <set>
#include <map>
#define GR_map std::map
#define GR_set std::set
#endif /* HAVE_CXX0X (else) */
#endif /* HAVE_TRI */
#endif /* __cplusplus */

#define TEMPLATE_DEFS_INCLUDED 

/* Be sure these two always get included along with config info. */
#include "GR_assert.h"
#ifdef __cplusplus
#include "GR_Classes.h"
#endif

#endif
 
