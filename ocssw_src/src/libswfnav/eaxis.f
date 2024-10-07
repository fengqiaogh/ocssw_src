        subroutine eaxis(e,phi,xm)
c  Computes coordinate transformation matrix corresponding to Euler axis e
c   and rotation angle phi

c  Reference:  Wertz, Appendix E

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  e(3)         R*4      I      Input Euler axis (must be a unit vector)
c  phi          R*4      I      Input rotation angle (degrees)
c  xm(3,3)      R*4      O      Output Transformation Matrix

        real xm(3,3),e(3)
        real*8 pi,radeg,re,rem,f,omf2,omegae
        common /gconst/pi,radeg,re,rem,f,omf2,omegae

        cp = cos(phi/radeg)
        sp = sin(phi/radeg)
        omcp = 1.0 - cp

        xm(1,1) = cp + e(1)*e(1)*omcp  
        xm(1,2) = e(1)*e(2)*omcp + e(3)*sp
        xm(1,3) = e(1)*e(3)*omcp - e(2)*sp
        xm(2,1) = e(1)*e(2)*omcp - e(3)*sp
        xm(2,2) = cp + e(2)*e(2)*omcp  
        xm(2,3) = e(2)*e(3)*omcp + e(1)*sp
        xm(3,1) = e(1)*e(3)*omcp + e(2)*sp
        xm(3,2) = e(2)*e(3)*omcp - e(1)*sp
        xm(3,3) = cp + e(3)*e(3)*omcp  

        return
        end
