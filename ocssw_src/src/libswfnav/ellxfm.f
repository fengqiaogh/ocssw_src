        subroutine ellxfm(attxfm, att_ang, tilt, p, navctl, smat, coef)
c
c  ellxfm(attxfm, att_ang, p, smat, coef)
c
c  Purpose: get sensor orientation matrix and scan ellipse 
c           coefficients for navigation of SeaWiFS data
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  attxfm       R*4      I      size 3 by 3 attitude transform matrix
c  att_ang      R*4      I      size 3 spacecraft yaw, roll and pitch
c  tilt         R*4      I      tilt angle 
c  p            R*4      I      size 3 spacecraft position from GPS
c  navctl       struct   I      controls for processing data
c  smat(3,3)    R*4      O      Sensor orientation matrix
c  coef(6)      R*4      O      Scan path coefficients
c
c  By: W. Robinson, GSC, 29 Mar 93
c
c  Notes: This is taken from the ORIENT routine by F Patt.  The 
c       attitude transform matrix replaces the computation of
c       the same made from the position and velocity (this is done 
c       elsewhere).  What follows is the description of ORIENT
c
c  This subroutine performs a simple calculation of the sensor
c  orientation from the orbit position vector and input values of the
c  attitude offset angles.  The calculations assume that the angles
c  represent the yaw, roll and pitch offsets between the local vertical
c  reference frame (at the spacecraft position) and the sensor frame.
c  Sensor tilt angles are assumed to be included in the pitch angle.
c  The outputs are the matrix which represents the transformation from
c  the geocentric rotating to sensor frame, and the coefficients which
c  represent the Earth scan track in the sensor frame. 
c
c  The reference ellipsoid uses an equatorial radius of 6378.137 km and
c  a flattening factor of 1/298.257 (WGS 1984). 
c
c  Modification History:
c
c  Added sensor offset matrix navctl.msenoff to calculation of sensor
c  transformation matrix - F. S. Patt, March 11, 1994.
c
c       Subprograms Called:
c
c       CROSSP          Compute cross product of two vectors
c       EULER           Compute matrix from Euler angles
c       MATMPY          Multiply two 3x3 matrices
c
      implicit none
#include "nav_cnst.fin"
#include "navctl_s.fin"
      type(navctl_struct) :: navctl

        real att_ang(3), attxfm(3,3)
        real smat(3,3),coef(6), p(3)
c
        real*8 rd
        real*4 sm1(3,3), sm2(3,3), sm3(3,3), tilt
        integer*4 i, j

        rd = 1.d0/omf2

c  Convert Euler angles to matrix
        call euler(att_ang,sm1)

c   Apply attitude offset matrix 
        call matmpy(sm1,attxfm,sm2)

c   Apply sensor offset matrix
        call matmpy(navctl%msenoff,sm2,sm1)

c  Compute rotation matrix for tilt angle and apply
        call eaxis(navctl%tiltcos,tilt,sm2)
        call matmpy(sm2,sm1,sm3)

c  Compute coefficients of intersection ellipse in scan plane
        coef(1) = 1.d0+(rd-1.d0)*sm3(1,3)*sm3(1,3)
        coef(2) = (rd-1.d0)*sm3(1,3)*sm3(3,3)*2.d0
        coef(3) = 1.d0+(rd-1.d0)*sm3(3,3)*sm3(3,3)
        coef(4) = (sm3(1,1)*p(1)+sm3(1,2)*p(2)+sm3(1,3)*p(3)*rd)*2.d0
        coef(5) = (sm3(3,1)*p(1)+sm3(3,2)*p(2)+sm3(3,3)*p(3)*rd)*2.d0
        coef(6) = p(1)*p(1)+p(2)*p(2)+p(3)*p(3)*rd-re*re

c  Transfer sensor orientation matrix to output array
        do i=1,3
            do j=1,3
                smat(i,j) = sm3(i,j)
            end do
        end do

        return
        end
