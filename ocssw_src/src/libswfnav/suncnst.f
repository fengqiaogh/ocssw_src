      subroutine suncnst(gaclac,navqc,sunrng,sun)
c
c  suncnst(gaclac,navqc,sunrng,sun)
c
c  Purpose: check consistency of the sun sensor data for the 3 sensors
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  gaclac       I*4      I      flag for GAC or LAC data.  If LAC
c                               data, there is 1 TLM for 3 lines
c                               else only once every 5 lines
c                               and the lines are 4scan lines apart
c  navqc        struct   I      navigation quality control info
c  sunrng       I*4     I/O     size 2 by 3 array of active range for
c                               the 3 sun sensors
c  sun          struct  I/O     sun sensor data structure
c
c  By: W. Robinson, GSC, 1 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Corrected logic of difference comparison to use absolute value.
c  F. S. Patt, GSC, December 3, 1997.
c
c  Modified tolerance multiplier for GAC data.
c  F. S. Patt, SAIC GSC, June 25, 1998.
c
c  Added logic to change sunrng according to first and last unflagged values.
c  F. S. Patt, SAIC GSC, August 14, 1998.
c
c  Added a check for maximum gap between samples, since the usefulness of this
c  check degrades with gap size.
c  F. S. Patt, SAIC GSC, August 25, 1998.


      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
      type(sun_struct) :: sun(3)
      type(navqc_struct) :: navqc
c
      integer*4 gaclac, sunrng(2,3)
c
      real*4 toldif(2)
c
      integer*4 i1, i2, j1, j2, nper, nrng, isens, maxgap
      logical found, end, gottwo
      real*4 tolmult, diff(2)
      data maxgap/6/

c
c
c       set up some controls
c
      nper = 1
      tolmult = 1.
      if( gaclac .eq. 1 ) then
        nper = 5       ! # actual lines per tlm line
c        tolmult = 4.   ! second tolerence multiplier
      end if
c
c       loop over the 3 sensors ( if they are active)
c
      do isens = 1,3
        if( sunrng(1,isens) .ne. -1 ) then
c
c           check the consistency of angles
c
          i2 = 0          ! i1, i2 are pointers to consecutive good values
          end = .FALSE.   ! true if the last unflagged value was found
          found = .FALSE. ! to signal that a consistent pair was found
          gottwo = .FALSE.! to signal that previous pair was consistent
          nrng = sunrng(2,isens) ! only go searching flags to end
c                                    of active range
c
c           start out by finding the next unflagged sun angle set
c
          call fndflg(sun(isens)%flag, nrng, sunrng(1,isens), i1 )
          if( i1 .le. 0 ) then
c
c              no unflagged values found at all, return with error
c
             write( 6, 100 ) isens
  100        format(' SUNCNST: no unflagged sun angle values found',
     1        /,' for sensor:',i7)
             sunrng(1,isens) = -1
             go to 980
          end if
c
c           place the next unflagged location in i2 and compare values
c
          do while( .not. end )
             call fndflg(sun(isens)%flag, nrng, (i1 + 1), i2 )
             if( i2 .le. 0 ) then
                end = .TRUE.
             else
c
c                 do the actual consistency checks
c
                diff(1) = sun(isens)%ang(1,i2) - sun(isens)%ang(1,i1)
                diff(2) = sun(isens)%ang(2,i2) - sun(isens)%ang(2,i1)
                toldif(1) = navqc%sun_del_1 * ( i2 - i1 ) *
     1            nper * tolmult
                toldif(2) = navqc%sun_del_2 * ( i2 - i1 ) *
     1            nper * tolmult
c
                if( ( abs(diff(1)) .gt. toldif(1) ) .or.
     1               ( abs(diff(2)) .gt. toldif(2) ) .or.
     2               ( (i2-i1) .gt. maxgap ) )  then
                   if ( .not. gottwo) sun(isens)%flag(i1) = 1
                   gottwo = .FALSE.
c                   sun(isens)%flag(i2) = 1
                else 
                   if (.not.found) then
                      j1 = i1
                      found = .TRUE.
                   end if
                   gottwo = .TRUE.
                   j2 = i2
                end if
c
c                 for next pair, move secind pointer to the first
                i1 = i2
             end if
          end do
c
c           make sure a consistent pair was found
c
          if( .not. found ) then
             write( 6, 500 ) isens
  500        format(' SUNCNST: no consistent pairs of sun sensor',/,
     1       ' angles were found for sensor:',i7)
             sunrng(1,isens) = -1
             go to 980
          else

c       check first and last unflagged values vs sunrng
             if (j1 .gt. sunrng(1,isens)) sunrng(1,isens) = j1
             if (j2 .lt. sunrng(2,isens)) sunrng(2,isens) = j2
          end if
c
  980     continue
        end if
      end do
c
c       and end
c
  990 continue
      return
      end
