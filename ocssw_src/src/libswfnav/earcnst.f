      subroutine earcnst(gaclac,navqc,earrng,earth)
c
c  earcnst(gaclac,navqc,earrng,earth)
c
c  Purpose: check consistency of the earth sensor data for the 2 sensors
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  gaclac       I*4      I      flag for GAC or LAC data.  If LAC
c                               data, there is 1 TLM for 3 lines
c                               else only once every 15 lines
c                               and the lines are 14 scan lines apart
c  navqc        struct   I      navigation quality control info
c  earrng       I*4     I/O     size 2 (low, hi) by 2 (sensor 1, 2) 
c                               array of active range for
c                               the 2 sun sensors
c  earth        struct  I/O     earth sensor data structure
c
c  By: W. Robinson, GSC, 13 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Corrected logic of difference comparison to use absolute value.
c  F. S. Patt, GSC, December 3, 1997.
c
c  Fixed bug which caused width tolerance to be used for both angles.
c  F. S. Patt, SAIC GSC, April 26. 1998
c
c  Added a check for maximum gap between samples, since the usefulness of this
c  check degrades with gap size.
c  F. S. Patt, SAIC GSC, September 14, 1998.
c
c  Added logic to change earrng according to first and last unflagged values.
c  F. S. Patt, SAIC GSC, February 3, 1999.

      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
      type(earth_struct) :: earth(2)
      type(navqc_struct) :: navqc
c
      integer*4 gaclac, earrng(2,2)
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
c       loop over the 2 sensors ( if they are active)
c
      do isens = 1,2
        if( earrng(1,isens) .ne. -1 ) then
c
c           check the consistency of the width and phase
c
          i2 = 0          ! i1, i2 are pointers to consecutive good values
          end = .FALSE.   ! true if the last unflagged value was found
          found = .FALSE. ! to signal that a consistent pair was found
          gottwo = .FALSE.! to signal that previous pair was consistent
          nrng = earrng(2,isens)   ! only go searching flags to end
c                                    of active range
c
c           start out by finding the next unflagged earth angle  set
c
          call fndflg(earth(isens)%flag, nrng, earrng(1,isens), i1 )
          if( i1 .le. 0 ) then
c
c              no unflagged values found at all, return with error
c
             write( 6, 100 ) isens
  100        format(' EARCNST: no unflagged earth angle values found',
     1        /,' for sensor:',i7)
             earrng(1,isens) = -1
             go to 980
          end if
c
c           place the next unflagged location in i2 and compare values
c
          do while( .not. end )
             call fndflg(earth(isens)%flag, nrng, (i1 + 1), i2 )
             if( i2 .le. 0 ) then
                end = .TRUE.
             else
c
c                 do the actual consistency checks
c
                diff(1) = earth(isens)%widphse(1,i2) - 
     1                    earth(isens)%widphse(1,i1)
                diff(2) = earth(isens)%widphse(2,i2) - 
     1                    earth(isens)%widphse(2,i1)
                toldif(1) = navqc%ear_del_wd* ( i2 - i1 ) *
     1            nper * tolmult
                toldif(2) = navqc%ear_del_ph* ( i2 - i1 ) *
     1            nper * tolmult
c
                if( ( abs(diff(1)) .gt. toldif(1) ) .or.
     1               ( abs(diff(2)) .gt. toldif(2) ) .or.
     2               ( (i2-i1) .gt. maxgap )     ) then
                   if ( .not. gottwo) earth(isens)%flag(i1) = 1
                   gottwo = .FALSE.
c                   earth(isens)%flag(i2) = 1
                else
                   if (.not.found) then
                      j1 = i1
                      found = .TRUE.
                   end if
                   gottwo = .TRUE.
                   j2 = i2
                end if

c
c                 for next pair, move second pointer to the first
                i1 = i2
             end if
          end do

c         check for last point passing check   
          if ( .not. gottwo) earth(isens)%flag(i1) = 1
c
c           make sure a consistent pair was found
c

          if( .not. found ) then
             write( 6, 500 ) isens
  500        format(' EARCNST: no consistent pairs of earth sensor',/,
     1       ' width, phase were found for sensor:',i7)
             earrng(1,isens) = -1
             go to 980
          else

c       check first and last unflagged values vs sunrng
             if (j1 .gt. earrng(1,isens)) earrng(1,isens) = j1
             if (j2 .lt. earrng(2,isens)) earrng(2,isens) = j2

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
