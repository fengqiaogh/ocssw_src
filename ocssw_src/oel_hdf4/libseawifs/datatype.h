/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/data/datatype.h,v 4.10 1995/01/17 19:57:54 seawifsd Exp seawifsd $
$Log: datatype.h,v $
Revision 4.10  1995/01/17 19:57:54  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:12  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:19  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:17  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/05/23 16:15:45  frank
included 'datatype_proto.h' for prototyping.

Revision 1.2  1994/05/10 18:47:18  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:29:50  seawifst
Initial revision


 */


/*
 * define data type -2(unknown), 0(normal LAC), 1(lunar cal), 2(solar cal),
 * 3(intergain cal), 4(detector check), 15(gac)
 */
#define UNKNOWN_TYPE -2
#define HRPTTYPE 0
#define LACTYPE  0
#define LUNTYPE  1
#define SOLTYPE  2
#define IGCTYPE  3
#define TDITYPE  4
#define GACTYPE  15
#define UNKNOWNTYPESTR "XXX"
#define LACSTR "LAC"
#define LUNSTR "LUN"
#define SOLSTR "SOL"
#define IGCSTR "IGC"
#define TDISTR "TDI"
#define GACSTR "GAC"
#define HRPTSTR "HRPT"
/* this is actually not a real data type but is used to identify 	*/
/* converted instrument telemetry file for CAL/VAL			*/
#define TLMSTR "TLM"

#define DATATYPE_SET(x) ((x == GACTYPE) || (x == HRPTTYPE) || (x == LACTYPE) || (x == LUNTYPE) || (x == SOLTYPE) || (x == IGCTYPE) || (x == TDITYPE))
#define NEW_DATATYPE(prev,next) (DATATYPE_SET(prev) && DATATYPE_SET(next) && (prev != next))

//#include	"usrhdr.h"
#include "usrmac.h"
#include "datatype_proto.h"
