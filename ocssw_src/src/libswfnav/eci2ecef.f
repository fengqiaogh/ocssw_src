      subroutine eci2ecef(time,ecivec,ecrvec)
c
c  Converts inertial Earth-centered coordinates of position and
c  velocity into Earth-centered, Earth-fixed (ECEF) coordinates,
c  using the Greenwich Hour Angle as the rotation angle.
c
      implicit none
c
      real*8    time            ! julian vector time
      real*8    ecivec(6)       ! eci vector (x,y,z,vx,vy,vz)
      real*8    ecrvec(6)       ! eci vector (x,y,z,vx,vy,vz)
c
      real*8    gha
      integer*4 year
      integer*4 day
      real*8    sec
c
      logical*4 firstCall
      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
      data firstCall / .true. /
c
      if (firstCall) then
          firstCall = .false.
          call cdata
      endif
c
      call jul2yds(time,year,day,sec)
      call gha2000(year,day+sec/86400.D0,gha)
c
      call ecef(gha,ecivec(1),ecivec(4),ecrvec(1),ecrvec(4))

c      rgha   = gha/radeg
c      cosgha = dcos(rgha)
c      singha = dsin(rgha)

c
c  Position rotation
c      ecrvec(1) =  ecivec(1)*cosgha + ecivec(2)*singha
c      ecrvec(2) = -ecivec(1)*singha + ecivec(2)*cosgha
c      ecrvec(3) =  ecivec(3)

c
c  Velocity rotation
c
c      term1 = ecivec(4)*cosgha + ecivec(5)*singha
c      term2 = omegae*(-ecivec(4)*singha + ecivec(5)*cosgha)
c      ecrvec(4) = term1+term2
c
c      term1 = -ecivec(4)*singha + ecivec(5)*cosgha
c      term2 = omegae*(-ecivec(4)*cosgha - ecivec(5)*singha)
c      ecrvec(5) = term1+term2
c
c      ecrvec(6) = ecivec(6)
c
      return
      end
