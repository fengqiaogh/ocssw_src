/*
 * $Header: /disk01/seadas/Repository/seadas/inc/swfinc/config.h,v 1.5 2006/10/06 04:23:23 seadas Exp $
 * $Log: config.h,v $
 * Revision 1.5  2006/10/06 04:23:23  seadas
 * MDM: cleaned up comments
 *
 * Revision 1.4  2006/10/06 03:46:38  seadas
 *  ifdef check for __APPLE__
 *
 * Revision 1.3  2006/07/14 17:03:53  seadas
 * MAR: svn update. 7/14/06.
 *
 * Revision 4.16  1995/12/07 17:36:11  seawifsd
 * updated support prototype for IRIX_5
 *
 * Revision 4.15  1995/05/31 17:41:39  seawifsd
 * typo NEED_BLKCKR -> NEED_BLKCLR.
 *
 * Revision 4.14  1995/05/12 13:39:49  seawifsd
 * added '#include <limits.h>' to get SHRT_MAX, SHRT_MIN, ...
 *
 * Revision 4.13  1995/05/11 15:57:08  seawifsd
 * added support that even when the upgraded OS is not been tested, hopefully
 * code can still be successfully compiled.
 *
 * Revision 4.12  1995/05/11 15:51:29  seawifsd
 * added support for IRIX64_6(64 bit IRIX running IRIX 6.x)
 *
 * Revision 4.11  1995/04/03 14:23:32  seawifsd
 * added support for OSF/1.
 *
 * Revision 4.10  1995/01/17 19:59:00  seawifsd
 * Jan. 17, 1994, V4.10
 *
 * Revision 4.1  1995/01/17 14:15:14  seawifsd
 * Jan. 9, 1994, 4.0
 *
 * Revision 3.4  1994/11/29 16:53:09  seawifsd
 * fixed prototyping of gethostname routine when some macro definitions are
 * in certain states(_SGI_SOURCE, _POSIX_SOURCE, _XOPEN_SOURCE).
 *
 * Revision 3.3  1994/11/08 18:47:15  seawifsd
 * Nov. 8, 1994, 3.3a3
 *
 * Revision 3.3  1994/11/08 15:05:12  seawifsd
 * Nov. 8, 1994, 3.3a2
 *
 */


#ifndef CONFIG_H_
#define CONFIG_H_

#include <stdio.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>

typedef uint32_t uaddr;

#ifdef IRIX_4

#ifdef PROTOTYPE
extern char *basename(char *path);
extern char *dirname(char *path);
extern int gethostname(char *name, int namelen);
#else
extern char *basename();
extern char *dirname();
extern int gethostname();
#endif /* PROTOTYPE */

#endif /* IRIX_4 -> IRIX 4.0.5C */


#if (defined(IRIX_5) || defined(IRIX64_6) || (defined(__sgi) && !defined(IRIX_4)))

/* basename(),dirname() */
#include <libgen.h>
#if defined(_SGI_SOURCE) && !defined(_POSIX_SOURCE) && !defined(_XOPEN_SOURCE)
/* gethostname() */
#include <unistd.h>
/* strcasecmp() */
#include <string.h>
#else
extern int gethostname(char *name, int namelen);
extern int strcasecmp(const char *, const char *);
extern int strncasecmp(const char *, const char *);
#endif /* _SGI_SOURCE && !_POSIX_SOURCE && !_XOPEN_SOURCE */

#endif /* IRIX_5 || IRIX64_6 || (!IRIX_4 && __sgi))*/


#if (defined(IRIX_4) || defined(IRIX_5) || defined(IRIX64_6) || defined(__sgi))
/* bzero(),bcopy(),bcmp(),blkclr() */
#include <bstring.h>
#endif /* IRIX_4 || IRIX_5 || IRIX64_6 || __sgi */

#ifdef SunOS_4
#endif /* SunOS 4.1.3 */


#ifdef SunOS_5
/* bzero */
extern void bzero(register char *sp, int len);
/* bzero only exist in BSD compatible mode */
#define NEED_BZERO
#define NEED_BLKCLR
#define NEED_SETLINEBUF
#define NEED_TRUNC
#endif /* Solaris -> SunOS 5.3 */

#ifdef OSF1_V3
#define NEED_CFTIME
#define NEED_BLKCLR
#define NEED__TIMEZONE
#define NEED__ALTZONE
#define NEED_ALTZONE
#endif /* OSF1_V3 */


/* define non-ANSI related */
#if (!defined(_POSIX_SOURCE) && !defined(__EXTENSIONS__))

/*
 * ANSI does not define this in limits.h but POSIX do
 */
#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

/*
 * In 4.0.5C, this is defined in string.h but said that is not in the
 * standard. In 5.1, this is not defined when you use -ansi
 */
extern char * strdup(const char *);

extern void tzset(void);
extern int cftime(char *, char *, const time_t *);

extern FILE * popen(const char *, const char *);
extern int pclose(FILE *);

#endif /* !_POSIX_SOURCE && !__EXTENSIONS__ */

/*
 * in SGI 5.2, setlinebuf is defined only when 
 * _SGI_SOURCE && !_POSIX_SOURCE && !_XOPEN_SOURCE
 *
 * while in SGI 4.0.5C, setlinebuf is defined when 
 * !_POSIX_SOURCE && __EXTENSIONS__
 */
#if defined(_POSIX_SOURCE) || !defined(__EXTENSION__)

/* FILE requires the inclusion of stdio.h */
#ifndef OSF1_V3
#if defined(LINUX) && !(defined(__APPLE__))
extern void setlinebuf(FILE *stream);
#else
extern int setlinebuf(FILE *stream);
#endif
#endif /* !OSF1_V3 */

/* trunc.h */
extern double trunc(double);
#ifndef OSF1_V3
/* #pragma no side effects (trunc) */
#endif /* !OSF1_V3 */

#endif /* _POSIX_SOURCE || !__EXTENSION__ */

#ifdef NEED_BZERO
#define bzero(p,l) {int i; register char *sp; sp = p ; for (i=0; i < l ; i++) {*sp = 0;sp++;}}
#endif /* NEED_BZERO */
#ifdef NEED_BLKCLR
#define blkclr(p,l) bzero(p,l)
#endif /* NEED_BLKCLR */
#ifdef NEED_SETLINEBUF
/*
#define setlinebuf(lun)	{char buf[BUFSIZ + 8]; setvbuf(lun,buf,_IOLBF,BUFSIZ);}
 */
#define setlinebuf(lun) setvbuf(lun,NULL,_IOLBF,0)
#endif /* NEED_SETLINEBUF */
#ifdef NEED_TRUNC
#define trunc(x) (((x) > 0)?floor(x):ceil(x))
#endif /* NEED_TRUNC */
#ifdef NEED_CFTIME
/*
   this is not complete correct if s is dynamically declared than sizeof(s)
   is always 4.
 */
#define cftime(s,fmt,tick) strftime(s,sizeof(s),fmt,localtime(tick))
#endif /* NEED_CFTIME */

#ifdef NEED__TIMEZONE
#define _timezone timezone
#endif /* NEED__TIMEZONE */

#ifdef NEED__ALTZONE
#define _altzone altzone
#endif /* NEED__ALTZONE */

#ifdef NEED_ALTZONE
/* this is not really a fix. But altzone is not used anyway. */
#define altzone  timezone
#endif /* NEED_ALTZONE */

#endif /* CONFIG_H_ */
