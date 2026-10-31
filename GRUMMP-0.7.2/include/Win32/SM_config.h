#ifndef VC__SM_CONFIG_H__INCLUDED
#define VC__SM_CONFIG_H__INCLUDED

#ifdef WIN32
#define WIN_VC32

/* symbol > 255 chars (DEBUG std::templates) */
#pragma warning(disable:4786) 

/* To enable logging, define "SUMAA_LOG=1" in VC IDE*/

/* Define to turn on statistic acculumation and output */
#define SM_STATS

/* Level of internal debugging output [0,1,2,3] */
#define SM_DEBUG_LEVEL (1)

/* Notify the smoothing code that it is -not- the */
/* highest level of logging, so that it will not */
/* add its own time to global time. */
#define SM_NOT_GLOBAL

/* Define to spew out matlab data on optimization */
#if (SM_DEBUG_LEVEL >= 3)
#define LOCALTEST 1
#else
#undef LOCALTEST
#endif

/* This is always defined for the auto-differentiation code. */
#define AD_DERIV_TYPE_ARRAY

/* Add missing constants */
/* #include "Win32_GR_values.h" */
#ifndef M_PI
#define M_PI  3.14159265358979323846  /* pi */
#endif

#endif /* WIN32 */
#endif /* VC__SM_CONFIG_H__INCLUDED */
