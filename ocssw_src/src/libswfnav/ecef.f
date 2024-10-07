      subroutine ecef(gha,posi,veli,pose,vele)
c
c  Converts inertial Earth-centered coordinates of position and
c  velocity into Earth-centered, Earth-fixed (ECEF) coordinates,
c  using the Greenwich Hour Angle as the rotation angle.
c
      implicit real*8 (a-h,o-z)
      dimension posi(3),veli(3)
      dimension pose(3),vele(3)
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
c
      rgha = gha/radeg
      cosgha = dcos(rgha)
      singha = dsin(rgha)
c
c  Position rotation
      pose(1) = posi(1)*cosgha + posi(2)*singha
      pose(2) = -posi(1)*singha + posi(2)*cosgha
      pose(3) = posi(3)
c
c  Velocity rotation
      term1 = veli(1)*cosgha + veli(2)*singha
      term2 = omegae*(-posi(1)*singha + posi(2)*cosgha)
      vele(1) = term1+term2
      term1 = -veli(1)*singha + veli(2)*cosgha
      term2 = omegae*(-posi(1)*cosgha - posi(2)*singha)
      vele(2) = term1+term2
      vele(3) = veli(3)
c
      return
      end
