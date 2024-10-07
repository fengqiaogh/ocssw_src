/*
$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.5/L012_Util/util/usr/wifs_conf.h,v 4.21 1996/05/07 14:28:03 seawifsd Exp seawifsd $
$Log: wifs_conf.h,v $
Revision 4.21  1996/05/07 14:28:03  seawifsd
defined REAL_TILT.  (NAK)

Revision 4.20  1996/01/16 22:26:12  seawifsd
fixed a bug that will make compiling to a wrong image when both
IO_SPEC_V44 and IO_SPEC_V43 are NOT defined.

Revision 4.19  1995/12/07 19:48:06  seawifsd
updated NONPROD_SPEC version

Revision 4.18  1995/12/07 17:32:13  seawifsd
updated to put IO_SPEC_V44 as the default.

Revision 4.17  1995/08/03 20:18:13  seawifsd
defined DONT_BREAK_TILT.

Revision 4.16  1995/05/31 17:43:00  seawifsd
added NO_BROWSE_REC_L2 and NO_GET_FLAG_NAMES to support IO_SPEC_V43

Revision 4.15  1995/05/12 13:28:56  seawifsd
1. upgraded to IO_SPEC_V43, PROD_SPEC_V27, and NONPROD_SPEC_V12.
2. eliminated macros related to version before IO_SPEC_V42 which includes
NO_SMOOTH_SCAN_TEMP_SIDE, QC_LA865, NO_GET_SCALE_OFFSET, NO_ADD_CSOL_Z,
NO_ENTRY_YEAR_DAY, and DONT_OVERWRITE_SLOPE_INTERCEPT

Revision 4.14  1995/05/08 17:33:38  seawifsd
made NO_CAL_OFFSET as default.

Revision 4.13  1995/02/21 16:47:35  seawifsd
set NO_SMOOTH_SCAN_TEMP_SIDE and NO_CAL_OFFSET for version before IO_SPEC_V42.

Revision 4.12  1995/02/15 20:37:34  seawifsd
elliminated IO_SPEC_V33 related code and rearrange to set IO_SPEC_V42
as the default.

Revision 4.11  1995/01/18 22:05:01  seawifsd
added new features REAL_CTIME and REAL_TILT(for future).

Revision 4.10  1995/01/17 19:59:11  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:15:23  seawifsd
Jan. 9, 1994, 4.0

Revision 3.9  1994/12/23 19:02:40  seawifsd
added check to make sure this file be included only once.

Revision 3.8  1994/12/15 20:12:52  seawifsd
made PROD_SPEC_V25 as the default.

Revision 3.7  1994/12/15 16:27:47  seawifsd
made IO_SPEC_V41 as current default and replace all V40 occurrences by V41.
added macro definitions for NONPROD_SPEC.

Revision 3.6  1994/12/06 19:26:41  seawifsd
made IO_SPEC_V40 is the current version and IO_SPEC_V33 is the last version.
made TRUE_ATTR_TYPE as default for IO_SPEC_V40.
made TRANSPOSE_IMAGE as default for IO_SPEC_V40 so that L1A image data are
pixel interlaced instead of band interlaced.

Revision 3.5  1994/11/28 18:22:32  seawifsd
made 'TRUE_ATTR_TYPE' defined as default.

Revision 3.4  1994/11/08 19:07:39  seawifsd
made ONCE_PER_PROCESS as default.

Revision 3.3  1994/11/08 18:47:24  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:05:20  seawifsd
Nov. 8, 1994, 3.3a2

 */


#ifndef WIFS_CONF_H_
#define WIFS_CONF_H_

/*
   This include file is used to attempt to make sure multiple versions of
   software can be selectablly compiled by specifying the right macros
   in the wifs_make_pre.inc file.
 */

/*
   versions for different I/O spec.
 */

#define LATEST_IO_SPEC IO_SPEC_V44
#define LAST_IO_SPEC IO_SPEC_V44
#define LATEST_PROD_SPEC PROD_SPEC_V28
#define LAST_PROD_SPEC  PROD_SPEC_V28
#define LATEST_NONPROD_SPEC NONPROD_SPEC_V13
#define LAST_NONPROD_SPEC NONPROD_SPEC_V13

/* if both undefined, use latest version */
#if (!defined(IO_SPEC_V44) && !defined(IO_SPEC_V43))
#define IO_SPEC_V44
#endif /* !IO_SPEC_V44 && !IO_SPEC_V43 */
/* if version is older than previous version, use previous version instead */
#if defined(IO_SPEC_V10) || defined(IO_SPEC_V33) || defined(IO_SPEC_V41) || defined(IO_SPEC_V42) || defined(IO_SPEC_V43)
#define IO_SPEC_V44
#ifdef IO_SPEC_V10
#define SHOULD_NOT_USE_IO_SPEC_V10_USE_IO_SPEC_V44_INSTEAD
#define SHOULD_NOT_USE_IO_SPEC_V10_USE_IO_SPEC_V44_INSTEAD
#undef IO_SPEC_V10
#endif /* IO_SPEC_V10 */
#ifdef IO_SPEC_V33
#define SHOULD_NOT_USE_IO_SPEC_V33_USE_IO_SPEC_V44_INSTEAD
#define SHOULD_NOT_USE_IO_SPEC_V33_USE_IO_SPEC_V44_INSTEAD
#undef IO_SPEC_V33
#endif /* IO_SPEC_V33 */
#ifdef IO_SPEC_V41
#define SHOULD_NOT_USE_IO_SPEC_V41_USE_IO_SPEC_V44_INSTEAD
#define SHOULD_NOT_USE_IO_SPEC_V41_USE_IO_SPEC_V44_INSTEAD
#undef IO_SPEC_V41
#endif /* IO_SPEC_V41 */
#ifdef IO_SPEC_V42
#define SHOULD_NOT_USE_IO_SPEC_V42_USE_IO_SPEC_V44_INSTEAD
#define SHOULD_NOT_USE_IO_SPEC_V42_USE_IO_SPEC_V44_INSTEAD
#undef IO_SPEC_V42
#endif /* IO_SPEC_V42 */
#ifdef IO_SPEC_V43
#define SHOULD_NOT_USE_IO_SPEC_V42_USE_IO_SPEC_V44_INSTEAD
#define SHOULD_NOT_USE_IO_SPEC_V42_USE_IO_SPEC_V44_INSTEAD
#undef IO_SPEC_V43
#endif /* IO_SPEC_V43 */
#endif /* IO_SPEC_V10 || IO_SPEC_V33 || IO_SPEC_V41 || IO_SPEC_V42 || IO_SPEC_V43 */


#if defined(IO_SPEC_V43) || defined(IO_SPEC_V44)

/* define this only when a new set of simulation data are available */
/* Fred Patt tells me it is available today (7-May-1996) Norman Kuring */
#define REAL_TILT

/* set dbopen once for every process as default */
#define ONCE_PER_PROCESS


#define LOCAL_PTIME

/* add valid_min,valid_max,min,max in MFSDObj */
#define MFSD_MINMAX
/* 3-code letter for different HRPT stations. Example: NSG -> _HNSG */
#define MULTI_HRPTTYPE
/* set MFSD->rec_flag to TRUE if the variable is a record variable */
#define REC_FLAG
/* store information about which dimension contain 'nsamp' */
#define NSAMP_DIMS
/* check boundary of image values */
#define CHK_IMAGE
/* calculate pixel min/max values */
#define CAL_PIXEL_MINMAX
/* to access every object */
/* #define EVERYTHING */


#endif /* IO_SPEC_V43 || IO_SPEC_V44 */



/* when not defined, additional attributes will be created to map the	*/
/* SeaWiFS slope/intercept values to NCSA HDF's convention.		*/
/* This will be defined indefinitely until there is a need.		*/
#define NO_CAL_OFFSET

/* when not defined, L1A scene will be created when the tilt changed	*/
#define DONT_BREAK_TILT

#endif /* WIFS_CONF_H_ */
