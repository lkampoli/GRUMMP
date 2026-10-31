#ifndef VC__SUMMA_CONFIG_H__INCLUDED
#define VC__SUMMA_CONFIG_H__INCLUDED

#ifdef WIN32
#define WIN_VC32

/* symbol > 255 chars (DEBUG std::templates) */
#pragma warning(disable:4786) 

/* To enable logging, define "SUMAA_LOG=1" in VC IDE */

/* Lite logging is used whenever logging is done. */
#ifdef SUMAA_LOG
#ifndef SUMAA_LOG_LITE
#define SUMAA_LOG_LITE
#endif
#endif

#endif /* WIN32 */
#endif /* VC__SUMMA_CONFIG_H__INCLUDED */
