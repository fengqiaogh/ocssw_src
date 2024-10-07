/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/osc/ecv.h,v 4.10 1995/01/17 19:58:40 seawifsd Exp seawifsd $
$Log: ecv.h,v $
Revision 4.10  1995/01/17 19:58:40  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:53  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:56  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:54  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.2  1994/05/10 18:51:33  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:49:00  seawifst
Initial revision


 */


/* Engineering ConVersion */
#ifndef _ECV_FILE
#define _ECV_FILE
/*----------------------------------------------------------------------------*/
/* Engineering Conversion Type for all the values			*/
/* There should be no need for discrete field, but ...			*/
/* Use negative for discrete field and positive for analog field	*/
/* 0 indicate no conversion						*/
#define NECT0     0

/* Discrete Engineering Conversion Type 1 - for discrete field		*/
#define DECT1     -1

/* Analog Engineering Conversion Type 1 - (Y = mX+b)			*/
/* two floating point values (m and b) are following the conversion	*/
/* type field								*/
#define AECT1     1
/* Analog Engineering Conversion Type 2 - for 4 byte GPS		*/
#define AECT2     2

/* define a filled parameters for NECT0(no conversion)			*/
#define FILLAECT1 NECT0,0.0,0.0

#ifdef __ANSI_CPP__
#define ECVM(n)  n##_ECV
#else
#define ECVM(n)  n##_ECV
#endif /* ANSI */

#endif /* _ECV_FILE */
