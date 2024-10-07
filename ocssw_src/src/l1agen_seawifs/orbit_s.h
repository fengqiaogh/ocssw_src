/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/mops/prepnav/orbit_s.h,v 4.11 1995/04/03 14:51:20 seawifsd Exp seawifsd $
$Log: orbit_s.h,v $
Revision 4.11  1995/04/03 14:51:20  seawifsd
changed references of 'long' to 'int' to support OSF/1.

Revision 4.10  1995/01/17 19:58:38  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:51  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:55  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:53  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/10/04 16:58:42  frank
made sure only be included once.
changed arrange size using macro name instead of the constant.

Revision 1.2  1994/05/10 18:51:06  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:38:02  seawifst
Initial revision


 */

#ifndef ORBIT_S_H_
#define ORBIT_S_H_

#include "swl0_parms.h"

typedef struct orbit_struct {
    double torb[MAXFRAMES];
    double pos[MAXFRAMES][3];
    double vel[MAXFRAMES][3];
    int nvec;
    int iyr;
    int iday;
} orbit_sType;

#endif /* ORBIT_S_H_ */
