        subroutine euler(a,xm)
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

        real xm1(3,3),xm2(3,3),xm3(3,3),xm(3,3),xmm(3,3),a(3)
        radeg = 180.d0/3.14159265359d0

c  Initialize all matrix elements to zero.      
        do i=1,3
          do j=1,3
            xm1(i,j) = 0.d0
            xm2(i,j) = 0.d0
            xm3(i,j) = 0.d0
          end do
        end do

c  Compute sines and cosines; use negative of Y and Z rotations
        c1=cos(a(1)/radeg)
        s1=sin(a(1)/radeg)
        c2=cos(a(2)/radeg)
        s2=-sin(a(2)/radeg)
        c3=cos(a(3)/radeg)
        s3=-sin(a(3)/radeg)

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
