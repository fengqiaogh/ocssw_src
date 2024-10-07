c --------------------------------------------------------------
c Subroutine EarthSunDist
c
c Returns Earth-Sun distance in AU.
c
c Modification History:
c
c - More accurrate formulation. B Franz, 26 June 1998.
c
c --------------------------------------------------------------
      real*8 function earthsundist(day,dutc)
c
      implicit none
c
      integer*4 day
      real*8    dutc 
      real*8    pi2
      real*8    rjd
c
      pi2    = 2.D0 * 3.141592654
      rjd    = day + dutc/24.0
c
c      EarthSunDist = 1.D0 + 0.0167 * cos(pi2*(rjd-3.0)/365.0)
c
      EarthSunDist = (1.00014 
     .      - 0.01671*cos(1.0*pi2*(0.9856002831*rjd-3.4532868)/360.0)
     .      - 0.00014*cos(2.0*pi2*(0.9856002831*rjd-3.4532868)/360.0))

c
      return
      end

