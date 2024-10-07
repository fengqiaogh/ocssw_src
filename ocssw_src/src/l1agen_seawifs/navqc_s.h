/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/mops/prepnav/navqc_s.h,v 4.11 1995/04/03 14:51:16 seawifsd Exp seawifsd $
$Log: navqc_s.h,v $
Revision 4.11  1995/04/03 14:51:16  seawifsd
changed references of 'long' to 'int' to support OSF/1.

Revision 4.10  1995/01/17 19:58:37  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:50  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:54  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:52  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/10/04 16:48:48  frank
made sure only be included once.

Revision 1.2  1994/05/10 18:51:03  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:37:57  seawifst
Initial revision


 */

#ifndef NAVQC_S_H_
#define NAVQC_S_H_

typedef struct navqc_struct {
    float sun_tol_1[2];
    float sun_tol_2[2];
    float sun_del_1;
    float sun_del_2;
    float ear_tol_wd[2];
    float ear_tol_ph[2];
    float ear_del_wd;
    float ear_del_ph;
    int yearmin;
    int yearmax;
    float sectol1;
    float sectol2;
    float sc_att[3][2];
    float att_del[3];
} navqc_sType;

#endif /* NAVQC_S_H_ */
