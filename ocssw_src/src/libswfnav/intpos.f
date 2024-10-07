c -------------------------------------------------------------------
c subroutine intpos
c
c Orbit vector interpolation.  
c
c Algorithm by Fred Patt.
c -------------------------------------------------------------------
      subroutine intpos(sec1,pos1,vel1,sec2,pos2,vel2,sec,pos,vel)
c
      implicit none
c
      real*8    sec1,sec2,sec
      real*4    pos1(3),pos2(3),pos(3)
      real*4    vel1(3),vel2(3),vel(3)
      real*4    a0(3),a1(3),a2(3),a3(3)
      real*8    dt
      real*8    x,x2,x3
      integer*4 j
c
      dt = sec2 - sec1
      do j=1,3
          a0(j) = pos1(j)
          a1(j) = vel1(j)*dt
          a2(j) = 3.D0*pos2(j)    - 3.D0*pos1(j) 
     .          - 2.D0*vel1(j)*dt - vel2(j)*dt
          a3(j) = 2.D0*pos1(j)    - 2.D0*pos2(j) 
     .          + vel1(j)*dt      + vel2(j)*dt
      end do
c
      x  = (sec - sec1)/dt
      x2 = x*x
      x3 = x2*x
c
      do j=1,3
          pos(j) = a0(j) + a1(j)*x + a2(j)*x2 + a3(j)*x3
          vel(j) = (a1(j) + 2.*a2(j)*x + 3.*a3(j)*x2)/dt
      end do
c
      return
      end






