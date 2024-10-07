/*

$Header: /app/shared/RCS/irix-5.2/seawifsd/src/hdfio/Shared.V4.2/L012_Util/util/mops/prepnav/navctl_s.h,v 4.11 1995/04/03 14:51:13 seawifsd Exp seawifsd $
$Log: navctl_s.h,v $
Revision 4.11  1995/04/03 14:51:13  seawifsd
changed references of 'long' to 'int' to support OSF/1.

Revision 4.10  1995/01/17 19:58:37  seawifsd
Jan. 17, 1994, V4.10

Revision 4.1  1995/01/17 14:14:50  seawifsd
Jan. 9, 1994, 4.0

Revision 3.3  1994/11/08 18:46:54  seawifsd
Nov. 8, 1994, 3.3a3

Revision 3.3  1994/11/08 15:04:52  seawifsd
Nov. 8, 1994, 3.3a2

Revision 1.1.1.1  1994/10/04 16:56:29  frank
made sure only be included once.
added a flag to turn on/off generating some debug files.

Revision 1.2  1994/05/10 18:50:59  seawifst
May 6, 1994 version 1.2

Revision 1.1  1994/04/19 13:37:52  seawifst
Initial revision


 */

#ifndef NAVCTL_S_H_
#define NAVCTL_S_H_

typedef struct navctl_struct {
    int procatt;
    int redoyaw;
    float yawtol;
    int nefpts;
    int neskip;
    int nsfpts;
    int nsskip;
    float msenoff[3][3];
    float tiltcos[3];
    float tiltcos2[3];
    float tiltfor;
    float tiltaft;
    float sun_mat[3][3][3];
    float sun_scal[3][2];
    float sun_bias[3][2];
    float ear_mat[2][3][3];
    float ear1sca;
    float ear2sca;
    float e1biasic;
    float e2biasic;
    float e1biasoc;
    float e2biasoc;
    int lvdbug;
} navctl_sType;

#endif /* NAVCTL_S_H_ */
