      subroutine get_xfm( smat, navctl, tilt, att, attxfm)
c
c  get_xfm( smat, msenoff, tilt, att, attxfm)
c
c  Purpose:  extract the ECEF-to-orbital transformation matrix from the
c  sensor transformation matrix.
c
c  Calling Arguments:
c     
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     smat(3,3)    R*4      I      size 3 by 3 sensor transformation matrix
c     navctl       struct   I      navigation control parameter structure
c     tilt         R*4      I      tilt angle
c     att(3)       R*4      I      size 3 array of attitude angles
c     attxfm(3,3)  R*4      O      size 3 by 3 ECEF-to-orbital frame matrix
c
c  By: F. S. Patt, SAIC GSC, 24 Sep 98
c
c  Notes:  This routine essentially performs the reverse of the matrix 
c  operations used by ellxfm, in order to obtain the ECEF-to-orbital 
c  transformation.
c
c  Modification History:
c
      implicit none
c
#include "navctl_s.fin"
c
      type(navctl_struct) :: navctl
      real*4 smat(3,3), att(3), attxfm(3,3), tilt
      real*4 sm1(3,3), sm2(3,3), sm3(3,3)

c  Compute rotation matrix for tilt angle, transpose and apply
      call eaxis( navctl%tiltcos, tilt, sm1)
      call xpose( sm1, sm2)
      
      call matmpy( sm2, smat, sm3)

c  Apply transpose of sensor offset matrix
      call xpose( navctl%msenoff, sm1)
      call matmpy( sm1, sm3, sm2)
     
c  Convert Euler angles to matrix
      call euler(att,sm3)
      call xpose( sm3, sm1)

c   Apply attitude offset matrix 
      call matmpy( sm1, sm2, attxfm)

      return
      end
