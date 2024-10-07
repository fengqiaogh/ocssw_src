        subroutine ydsadd(iy, id, sec, deld, delsec)
c
c  ydsadd(iy, id, sec, deld, delsec)
c
c  Purpose: add a delta day, second to a yesr, day, second
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  iy           I*4     I/O     year value to be changed
c  id           I*4     I/O     day of year to be changed
c  sec          R*8     I/O     second of day to be changed
c  deld         I*4      I      delta day to add
c  delsec       R*8      I      delta second to add
c
c  By: W. Robinson, GSC, 19 Mar 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
c
      integer*4 iy, id, deld
      real*8 sec, delsec
c
      integer*4 dnew, i4y, julday, jd, mon, day, refday
      real*8 snew, secinday
      parameter (secinday = 60 * 60 * 24 )
c
c       find new day with addition
c
      dnew = id + deld
c
c       convert this to julian day
c
      i4y = iy
      julday = jd( i4y, 1, dnew )
c
c       add the seconds and adjust seconds and julian day if required
c
      snew = sec + delsec
      julday = julday + int( snew / secinday )
      sec = dmod(snew, secinday)
c
      if( sec .lt. 0 ) then
         sec = sec + secinday
         julday = julday - 1
      end if
c
c       convert the julian day to a date
c
      call jddate( julday, i4y, mon, day )
      refday = jd( i4y, 1, 1)
      dnew = julday - refday + 1
c
      id = dnew
      iy = i4y
c
c       and exit
c
      return
      end
