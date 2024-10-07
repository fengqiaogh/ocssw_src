/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.5/L012_Util/util/usr/usrmac.h,v 4.13 1995/08/03 20:15:51
$Log: usrmac.h,v $
Revision 4.14  1996/01/17 16:48:10  seawifsd
added TIME_REPORT feature that when macro ON_TIME() and
OFF_TIME() are used in the code, unix time will be printed
and difference will be printed also.


Revision 4.13  1995/08/03 20:15:51  seawifsd
added macros MSEC2SEC and SEC2MSEC.

Revision 4.12  1995/04/03 14:21:04  seawifsd
added support to find out duplicate memory free.

Revision 4.11  1995/02/15 20:34:30  seawifsd
defined AX_PLUS_B(X,A,B), R_AX_PLUS_B(Y,A,B) and INRANGE(v,lo,hi)

Revision 4.10  1995/01/17 19:59:06  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:15:20  seawifsd
Jan. 9, 1994, 4.0

Revision 3.7  1995/01/06 18:52:10  seawifsd
re-defined CHKPRT1/CHKPRT2 macros for ANSI/non-ANSI situation.

Revision 3.6  1994/12/15 16:50:45  seawifsd
added macro definitions for MALLOCFAILEXIT,CHKMALLOC,and CHKCALLOC.

Revision 3.5  1994/11/28 18:23:11  seawifsd
fixed macro NOPRINT definition caused by changing the default from PRINT to
NOPRINT.

Revision 3.4  1994/11/09 18:36:08  seawifsd
changed default from 'debuging if NOPRINT is not defined' to
'debuging only if PRINT is defined'.

Revision 3.3  1994/11/08 18:47:21  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:05:18  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.4  1994/11/03 20:16:46  frank
To prevent conflicts in the definition of MIN/MAX/MOD macros, those
macros are renamed to MINV,MAXV and MODV respectively.

Revision 1.1.1.3  1994/10/04 15:38:39  frank
made sure macros TRUE and FALSE are properly defined.
added '#include "wifs_conf.h"' to deal with application specific features.

Revision 1.1.1.2  1994/05/23 14:08:20  frank
include the platform/compiler dependent include file config.h.

Revision 1.2  1994/05/11 15:54:39  seawifst
May 6, 1994, version 1.2
        Summary 1.1.1.1 to 1.1.1.1(total 1 revisions):
 * Make sure PROTOTYPE is set when __STDC__ is defined.

Revision 1.1.1.1  1994/05/06 14:23:11  frank
Make sure PROTOTYPE is set when __STDC__ is defined.

Revision 1.1  1994/04/19 13:52:17  seawifst
Initial revision


 */


#ifndef USRMAC_H_
#define USRMAC_H_

/*
#define FAKETILT
 */
/*
#define	FAKEGAIN
 */
/*
#define FIXSETDIMS
 */
#ifndef PRIVATE
#define PRIVATE static
#endif
#ifndef bool
#define bool int
#endif
#if !defined(TRUE) || ((TRUE) != 1)
#define TRUE (1)
#endif
#if !defined(FALSE) || ((FALSE) != 0)
#define FALSE (0)
#endif

#define S_FILE_FOUND  -3
#define S_FILE_NOT_FOUND -2
#define S_ERROR   -1
#define S_SUCCESS  0
#define S_END_OF_FILE  1


#ifdef PRINT
#undef NOPRINT
#define PRINTF printf("(%s,%d) ",__FILE__,__LINE__); printf
#define FPRINTF fprintf(stderr,"(%s,%d) ",__FILE__,__LINE__); fprintf
#else
#define NOPRINT 1
#define PRINTF if (!NOPRINT) printf
#define FPRINTF if (!NOPRINT) fprintf
#endif /* PRINT */

/*
  Example:
   CHKPRT1("%d",SDS->rank) -> PRINTF("%s=""%d","SDS->rank",SDS->rank)
   CHKPRT2("%d,%d",SDS->rec,i) -> PRINTF("%s=""%d,%d","SDS->rec,i",SDS->rec,i)
 */
#ifdef __STDC__
#define CHKPRT1(FMT,ARG) PRINTF("%s="# FMT,# ARG,ARG)
#define CHKPRT2(FMT,ARG1,ARG2) PRINTF("%s="# FMT,# ARG1 #ARG2,ARG1,ARG2)
#define CHKIPRT(EX) PRINTF("%s=%d\n", # EX ,EX)
#define CHKFPRT(EX) PRINTF("%s=%f\n", # EX ,EX)
#define CHKSPRT(EX) PRINTF("%s=%s\n", # EX ,EX)
#else
#define CHKPRT1(FMT,ARG) PRINTF("%s="FMT,"ARG",ARG)
#define CHKPRT2(FMT,ARG1,ARG2) PRINTF("%s="FMT,"ARG1,ARG2",ARG1,ARG2)
#define CHKIPRT(EX) PRINTF("%s=%d\n","EX",EX)
#define CHKFPRT(EX) PRINTF("%s=%f\n","EX",EX)
#define CHKSPRT(EX) PRINTF("%s=%s\n","EX",EX)
#endif /* ANSI */

#define MAXV(a,b) (((a) > (b)) ? (a) : (b))
#define MINV(a,b) (((a) > (b)) ? (b) : (a))
#define MODV(a,b) ((a) - ((a)/(b))*(b))

#define AX_PLUS_B(X,A,B) (A) * (X) + (B)
#define R_AX_PLUS_B(Y,A,B) ((Y) - (B)) / (A)
#define INRANGE(v,lo,hi) (((v - lo) >= 0) && ((hi - v) >= 0))

#ifdef CHECK_DUPLICATE_FREE
#ifndef MM_PTR
#define MM_PTR
static void *mm_ptr = (void *) 0x01;
#endif /* !MM_PTR */
#define free(x) if ((x != NULL) && (x != mm_ptr)) { free(x) ; x = mm_ptr;} else {fprintf(stderr,"%s %d: '%s' already free'd\n",__FILE__,__LINE__,#x);}
#endif /* CHECK_DUPLICATE_FREE */

#ifdef TIME_REPORT
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

static struct timeval stime_struct;
static struct timeval etime_struct;
static clock_t sclock;
static clock_t eclock;

#define ON_TIME() gettimeofday(&stime_struct);sclock=clock();fprintf(stderr,"%s :  ON : %ld.%ld %ld\n",FUNC,stime_struct.tv_sec,stime_struct.tv_usec,sclock)
#define OFF_TIME() gettimeofday(&etime_struct);eclock=clock();fprintf(stderr,"%s : OFF : %ld.%ld %f %ld %ld\n",FUNC,etime_struct.tv_sec,etime_struct.tv_usec,(etime_struct.tv_sec - stime_struct.tv_sec)+(etime_struct.tv_usec - stime_struct.tv_usec)/1.0E6,eclock,(eclock - sclock))
#else
#define ON_TIME()
#define OFF_TIME()
#endif /* TIME_REPORT */

#define TRACE  trace

#ifdef __STDC__
#ifndef PROTOTYPE
#define PROTOTYPE
#endif /* PROTOTYPE */
#endif /* __STDC__ */

#ifdef PROTO
#undef PROTO
#endif /* PROTO */
#ifdef PROTOTYPE
#define  PROTO(x) x
#else
#define  PROTO(x) ()
#endif

#ifndef byte
#define byte unsigned char
#endif

/*
   define size of different data type
 */
#define SC sizeof(char)
#define SB sizeof(byte)
#define SS sizeof(short int)
#define SI sizeof(int)
#define SL sizeof(int32_t int)
#define SF sizeof(float)
#define SD sizeof(double)

/*
   Following macros are used to indicate where in the program, we use
   fake routines to supply data or we don't process the data yet, or
   the definition of some constants is not confirmed, if
   you see compiler complaining about 'NODATA','FAKEDATA', or 'FAKEDEF'
   are been redefined. 
 */
#define NODATA
#define FAKEDATA
#define FAKEDEF

/*
   for code need to be resolve
 */
#define ISSUE


/*
   Following macros are used to turn on/off addition information, or do
   some additional operations.
   NOTE: you must use OPTIONS_START and OPTIONS_STOP in pairs.
   default is turn off the OPTIONS if not defined. i.e. all the codes
   in between OPTIONS_START and OPTIONS_STOP will not be executed.
   You can specify it in the compile time also by
   cc -DOPTIONS=OPTIONS_OFF		; turn off
   cc -DOPTIONS=OPTIONS_ON		; turn on
 */
#define OPTIONS_ON 1
#define OPTIONS_OFF 0
#ifndef OPTIONS
#define OPTIONS OPTIONS_OFF
#endif
#define OPTIONS_START(); if (OPTIONS) {
#define OPTIONS_STOP();  }

#define MSEC2SEC(x) ((x + 500)/1000)
#define SEC2MSEC(x) (x * 1000)

#include "config.h"
#include "wifs_conf.h"
#endif /* USRMAC_H_ */
