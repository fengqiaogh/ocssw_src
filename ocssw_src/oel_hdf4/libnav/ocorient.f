        subroutine ocorient(pos,vel,att,rm,coef)

c  This subroutine performs a simple calculation of the sensor
c  orientation from the orbit position vector and input values of the
c  attitude offset angles.  The calculations assume that the angles
c  represent the roll, pitch and yaw offsets between the local vertical
c  reference frame (at the spacecraft position) and the sensor frame. 
c  Sensor tilt angles are assumed to be included in the pitch angle. 
c  The outputs are the matrix which represents the transformation from
c  the geocentric rotating to sensor frame, and the coefficients which
c  represent the Earth scan track in the sensor frame.

c  The reference ellipsoid uses an equatorial radius of 6378.137 km and
c  a flattening factor of 1/298.257 (WGS 1984). 


c  Calling Arguments

c  Name         Type    I/O     Description
c
c  pos(3)       R*4      I      Orbit Position Vector ECEF (km)
c  vel(3)       R*4      I      Orbit Velocity Vector ECEF (km/sec)
c  att(3)       R*4      I      Attitude Offsets (Euler angles in 
c                                degrees); referenced to local vertical 
c                                coordinates; order is roll, pitch, yaw 
c                                (X, Y, Z)
c  rm(3,3)      R*4      O      Sensor orientation matrix
c  coef(10)     R*4      O      Scan path coefficients
c
c       Subprograms Called (attached):

c       CROSSP          Compute cross product of two vectors
c       EULER           Compute matrix from Euler angles
c       MATMPY          Multiply two 3x3 matrices
c
c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               July 22, 1992
c
c       Modification History:
c
c       Added improved calculation of local vertical reference frame and
c        modified calling argument names.  
c        F. S. Patt, September 30, 1992
c
c       Expanded vector normalization in-line to eliminate subroutine 
c        call.  F. S. Patt, October 19, 1992
c
c       Eliminated redundant calculations, changed to correspond to
c       paper, "Exact closed-form geolocation algorithm for earth
c       survey sensors", International Journal of Remote Sensing, Patt
c       and Gregg, 1993.  W. Gregg, 4/5/93.
c
c       Modified to support three-dimensional view vectors by computing
c       coefficients array with all 10 ellipsoid terms.  
c       F. S. Patt, November 25, 1996
c
      real*8 pi,radeg,re,rem,f,omf2,omegae
      real vc(3),xpri(3),ypri(3),zpri(3),yrp(3,3),rn(3,3),rm(3,3)
      real pos(3),vel(3),att(3),coef(10)
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
c
c  Compute correction to orbit velocity vector in Earth-centered
c   Earth-fixed (ECEF) frame; this involves subtracting effect of Earth
c   rotation rate on velocity to get correct scan plane orientation in
c   ECEF frame.
      vc(1) = vel(1) - omegae*pos(2)
      vc(2) = vel(2) + omegae*pos(1)
      vc(3) = vel(3)

c  Determine nadir frame reference axes
c  Uses method of local ellipsoid approximation good to 0.3 arcsecond
c   Compute Z axis as local nadir vector
        pm = dsqrt(DBLE(pos(1)*pos(1)+pos(2)*pos(2)+pos(3)*pos(3)))
        omf2p = (omf2*rem + pm - rem)/pm
        pxy = pos(1)*pos(1)+pos(2)*pos(2)
        temp = dsqrt(DBLE(pos(3)*pos(3) + omf2p*omf2p*pxy))
        zpri(1) = -omf2p*pos(1)/temp
        zpri(2) = -omf2p*pos(2)/temp
        zpri(3) = -pos(3)/temp
c   Compute Y axis along negative orbit normal
        call crossp(vc,zpri,ypri)
        yprim = dsqrt(DBLE(ypri(1)*ypri(1) + ypri(2)*ypri(2) +
     +                     ypri(3)*ypri(3)))
        ypri(1) = -ypri(1)/yprim
        ypri(2) = -ypri(2)/yprim
        ypri(3) = -ypri(3)/yprim
c   Compute X axis to complete orthonormal triad
        call crossp(ypri,zpri,xpri)
c  Store in matrix 
        do i=1,3
            rn(1,i)=xpri(i)
            rn(2,i)=ypri(i)
            rn(3,i)=zpri(i)
        end do
c
c  Convert attitude (Euler) angles to YRP matrix
        call oceuler(att,yrp)

c  Compute sensor orientation matrix 
        call matmpy(yrp,rn,rm)

c  Compute coefficients of intersection ellipse in scan plane
        rd=1.d0/omf2
        coef(1) = 1.d0+(rd-1.d0)*rm(1,3)*rm(1,3)
        coef(2) = 1.d0+(rd-1.d0)*rm(2,3)*rm(2,3)
        coef(3) = 1.d0+(rd-1.d0)*rm(3,3)*rm(3,3)
        coef(4) = (rd-1.d0)*rm(1,3)*rm(2,3)*2.d0
        coef(5) = (rd-1.d0)*rm(1,3)*rm(3,3)*2.d0
        coef(6) = (rd-1.d0)*rm(2,3)*rm(3,3)*2.d0
        coef(7) = (rm(1,1)*pos(1)+rm(1,2)*pos(2)+rm(1,3)*
     *            pos(3)*rd)*2.d0
        coef(8) = (rm(2,1)*pos(1)+rm(2,2)*pos(2)+rm(2,3)*
     *            pos(3)*rd)*2.d0
        coef(9) = (rm(3,1)*pos(1)+rm(3,2)*pos(2)+rm(3,3)*
     *            pos(3)*rd)*2.d0
        coef(10) = pos(1)*pos(1)+pos(2)*pos(2)+pos(3)*pos(3)*rd-re*re
c
        return
        end

c       subroutine crossp(v1,v2,v3)
c  originally here but removed as it is a independent routine in library

        subroutine oceuler(a,xm)
c  Computes coordinate transformation matrix corresponding to Euler 
c  sequence.  The order of angles in the input array is yaw, roll, 
c  pitch; according to OSC, the order of the rotations is the reverse
c  of this; the roll and pitch angles are about the negative Y and Z
c  axes, respectively, while the yaw angle is about the positive X axis. 

c  Reference:  Wertz, Appendix E; OSC, personal communication 

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  a(3)         R*4      I      Input Array of Euler Angles (degrees)
c  xm(3,3)      R*4      O      Output Transformation Matrix
c
c  Frederick S. Patt, GSC, sometime in 1992.
c
c  Modification history:
c
c
c  Modified to change order of rotations to pitch, roll, yaw (-Z, -Y, -X)
c  F. S. Patt, GSC, September 29, 1996.
c
c  Removed negative signs on Y and Z rotations for OCTS; order of rotations
c  is yaw, pitch, roll (Z, Y, X)

        real xm1(3,3),xm2(3,3),xm3(3,3),xm(3,3),xmm(3,3),a(3)
        real*8 radeg
        radeg = 180.d0/3.14159265359d0

c  Initialize all matrix elements to zero.      
        do i=1,3
          do j=1,3
            xm1(i,j) = 0.d0
            xm2(i,j) = 0.d0
            xm3(i,j) = 0.d0
          end do
        end do

c  Compute sines and cosines
        c1=dcos(a(1)/radeg)
        s1=dsin(a(1)/radeg)
        c2=dcos(a(2)/radeg)
        s2=dsin(a(2)/radeg)
        c3=dcos(a(3)/radeg)
        s3=dsin(a(3)/radeg)

c  Convert individual rotations to matrices
        xm1(1,1)=1.d0
        xm1(2,2)=c1
        xm1(3,3)=c1
        xm1(2,3)=s1
        xm1(3,2)=-s1
        xm2(2,2)=1.d0
        xm2(1,1)=c2
        xm2(3,3)=c2
        xm2(3,1)=s2
        xm2(1,3)=-s2
        xm3(3,3)=1.d0
        xm3(2,2)=c3
        xm3(1,1)=c3
        xm3(1,2)=s3
        xm3(2,1)=-s3

c  Compute total rotation as xm1*xm2*xm3 
        call matmpy(xm2,xm3,xmm)
        call matmpy(xm1,xmm,xm)
        return
        end

        subroutine matmpy(xm1,xm2,xm3)
c  Computes matrix product of 3x3 matrices xm1 and xm2

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  xm1(3,3)     R*4      I      Input Matrix
c  xm2(3,3)     R*4      I      Input Matrix
c  xm3(3,3)     R*4      O      Output Matrix

        real xm1(3,3),xm2(3,3),xm3(3,3)
        do i=1,3
            do j=1,3
                xm3(i,j) = 0.d0
                do k=1,3
                    xm3(i,j) = xm3(i,j) + xm1(i,k)*xm2(k,j)
                end do
            end do
        end do
        return
        end
