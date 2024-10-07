      subroutine proctim2(input,nframes,navqc,time,timref,
     1     nlines,iret)
c
c  proctim2(input,nframes,navqc,time,timref,nlines,iret)
c
c  Purpose: process SeaStar time tag data into times for every scan line
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  input        struct   I      input data structure containing 
c                               SeaStar ID and time tag
c  nframes      I*4      I      number of lines of time data
c  navqc        struct   I      navigation quality control info
c  time         R*8      O      array of time in seconds relative to 
c                               timref for every scan line
c  timref       R*8      O      size 3 reference time at start line
c                               of data: year, day, sec
c  nlines       I*4      O      number of scan lines covered
c  iret         I*4      O      return code, 0 - good
c
c  By: W. Robinson, GSC, 18 Mar 93
c
c  Notes:  
c
c  Modification History:
c
c  10 August 1993 - Modified by Frederick S. Patt, GSC to process input
c    SeaStar time tags
c
c  Commented out redundant call to fndflg.  F.S. Patt, GSC, August 15, 1996.
c
c  Added time offset of 4 LAC lines to time tags (to remove navigation
c  residuals) and removed check for consecutive scan times (does not work 
c  unless scans are filled and is otherwise not needed).  
c  F.S. Patt, SAIC GSC, February 24, 1998.
c
c  Modified time offset to be 6 LAC lines (1 second).
c  F.S. Patt, SAIC GSC, January 19, 1999.
c  
c  Moved addition of time offset to end of routine (added to timref(3)) 
c  to fix occasional problem with flagged frames at the end of the day.
c  F.S. Patt, SAIC GSC, May 31, 2001.
c
c  Fixed bug in logic which determined if first unflagged time is after 
c  a day crossing (changed .gt. to .ge.).  F.S. Patt, SAIC GSC, May 31, 2001.
c
c  Added capability to apply fixed time shifts for specified periods to 
c  correct extended spacecraft time error periods.
c  F.S. Patt, SAIC, October 1, 2002.

      implicit none
#include "nav_cnst.fin"
#include "input_s.fin"
#include "timtag_s.fin"
#include "navqc_s.fin"
      type(input_struct) :: input(maxlin)
      type(timtag_struct) :: timtag(maxlin)
      type(navqc_struct) :: navqc
c
      integer*4 gaclac, ntim, itim, iret
      real*8 time(maxlin), timref(3)
      real*8 secinday
c
      integer*4 nframes, flag(maxlin), flagsec(maxlin), daycross
      integer*4 i1, i2, crosspt, nper, ilin, nlines, jd
      integer*4 iy1, iy2, id1, id2
      logical found, endsec
      real*8 delsec, tolmult, delfac
      real*8 r8sec, r8dels

      parameter (secinday = 60.d0 * 60.d0 * 24.d0 )
      parameter (delsec = 1.d0 / 6.d0)
      parameter (delfac = 6.d0)
c
c
c       initialize flag arrays for seconds and for day/year
c

      if (input(1)%sc_id(2).eq.15) then
        gaclac = 1
      else
        gaclac = 0
      end if

      ntim = nframes
      iret = 0
      nper = 1
      tolmult = 1.
      if( gaclac .eq. 1 ) then
         nper = 5       ! # actual lines per tlm line
         tolmult = 4.d0   ! second tolerence multiplier
      end if
c
      do itim = 1, ntim
         flag(itim) = 1
         flagsec(itim) = 1
         timtag(itim)%sec = input(itim)%msec/1.d3
      end do
      nlines = ntim * nper
c
c       flag second, day or year values outside tolerence
c
      do itim = 1, ntim
         if( ( input(itim)%iyear .ge. navqc%yearmin )  .and.
     1       ( input(itim)%iyear .le. navqc%yearmax )  .and.
     1       ( input(itim)%iday .ge. 1 )               .and. 
     1       ( input(itim)%iday .le. 366 )  ) flag(itim) = 0
c
         if( ( timtag(itim)%sec .ge. 0 )  .and.
     1       ( timtag(itim)%sec .le. secinday )  ) flagsec(itim) = 0
      end do
c
c       check the consistency of the second values
c
      daycross = 0   ! day crossing flag
      crosspt = 0    ! index of crossing point
      i2 = 0          ! i1, i2 are pointers to consecutive good sec values
      endsec = .FALSE.   ! true if the last unflagged sec was found
      found = .FALSE. ! to signal that a consistent pair was found
c
c       start out by finding the next unflagged second value
c
      call fndflg(flagsec, ntim, 1, i1 )
      if( i1 .le. 0 ) then
c
c          no unflagged values found at all, return with error
c
         iret = -1
         write( 6, 100 ) 
  100    format(' PROCTIM: no unflagged tag second values found')
         go to 990
      end if
c
c       place the next unflagged second location in i2 and compare values
c
      do while( .not. endsec )
         call fndflg(flagsec, ntim, (i1 + 1), i2 )
         if( i2 .le. 0 ) then
            endsec = .TRUE.
         else
            found = .TRUE.  
c
c             adjust the second value past any day crossing detected
            if( daycross .ne. 0 ) 
     1          timtag(i2)%sec = timtag(i2)%sec + secinday
c
c                   in this success case, the day boundary was crossed
            iy1 = input(i1)%iyear
            id1 = input(i1)%iday
            iy2 = input(i2)%iyear
            id2 = input(i2)%iday
            if ( (jd(iy2,0,id2) - jd(iy1,0,id1)) .eq. 1) then
                  daycross = 1
                  crosspt = i2
                  timtag(i2)%sec = timtag(i2)%sec + secinday
c               end if
            end if
c
c             for next pair, move secind pointer to the first
            i1 = i2
         end if
      end do
c
c       make sure a consistent second pair was found
c
      if( .not. found ) then
         iret = -1
         write( 6, 500 )
  500    format(' PROCTIM: no consistent pairs of tag seconds were 
     1       found')
         go to 990
      end if
c
c       interpolate the second values to the output time array
c
c       first, move the second values from the timetag array to 
c       the output array
c
      do ilin = 1,ntim
         time( ( ilin - 1 ) * nper + 1 ) = timtag( ilin )%sec
      end do
c
c       extrapolate any lines required at the start of the segment
c
c          use first 2 good times to extrapolate or interpolate below
c
      call fndflg(flagsec, ntim, 1, i1 )
c
      if( flagsec(1) .eq. 1) then
c
         do ilin = 1, (i1 - 1) * nper
            time( ilin ) = timtag(i1)%sec - delsec * tolmult *
     1                    ( (i1 - 1) * nper + 1 - ilin )
         end do
      end if
c
c       interpolate through the available times, start with first 2 times
c       found above
c
      endsec = .FALSE.
      do while( .not. endsec )
         call fndflg(flagsec, ntim, (i1 + 1), i2 )
         if( i2 .le. 0 ) then
            endsec = .TRUE.
         else
            if( ( i2 * nper - i1 * nper ) .gt. 1 ) then
c                there are spaces to fill in output time array
               do ilin = (i1 - 1) * nper + 2, (i2 - 1) * nper
                  time(ilin) = timtag(i1)%sec + delsec * tolmult *
     1                         ( ilin - ( i1 - 1 ) * nper - 1 )
               end do
            end if
c
c             find next pair to interpolate
            i1 = i2
         end if
      end do
c
c       extrapolate times to the end of the segment
c
      if( ( flagsec(ntim) .eq. 1 )  .or.  ( gaclac .eq. 1 ) ) then
         do ilin = (i1 - 1 ) * nper + 2, ntim * nper
            time(ilin) = timtag(i1)%sec + delsec * tolmult *
     1                   ( ilin - ( i1 - 1 ) * nper - 1 )
         end  do
      end if
c
c       locate consistent day and year with a good seconds value for 
c       the first value in a pair
c
      call fndflg(flag, ntim, 1, i1 )
      if( i1 .le. 0 ) then
         iret = -1
         write( 6, 300 ) 
  300    format(' PROCTIM: no unflagged tag year/day values found')
         go to 990
      end if
c
      found = .FALSE.
c
      do while( .not. found )
         call fndflg(flag, ntim, (i1 + 1), i2)
         if( i2 .le. 0 ) then
            write( 6, 400 )
  400       format(' PROCTIM: unable to find a good year/ day sequence ',
     1        'in the segment')
            go to 990
         else
            if( ( input(i1)%iyear .ne. input(i2)%iyear )  .or.
     1      ( input(i1)%iday .ne. input(i2)%iday )  .or.
     1      ( flagsec(i1) .eq. 1 ) ) then
c
c             a start was not found, set the next pair and search again
               i1 = i2
            else
               found = .TRUE.
            end if
         end if
      end do
c
c       adjust the year, day for any day boundary crossings
c       by subtracting 1 day from the date
c
      r8sec = 0.
      r8dels = 0.
      iy1 = input(i1)%iyear
      id1 = input(i1)%iday
c     
      if( ( crosspt .gt. 0 ) .and. ( i1 .ge. crosspt ) ) 
     1     call ydsadd( iy1, id1, r8sec, -1, r8dels )
c
      timref(1) = iy1
      timref(2) = id1 
c
c       reconcile the times to start at the reference time
c
      timref(3) = time(1)
      nlines = ntim * nper
      do ilin = 1, nlines
         time(ilin) = time(ilin) - timref(3)
      end do

c       add timetag offset to start time
      timref(3) = timref(3) + delfac*delsec
c
c       check for time shift period
      call checkshift (timref)

c       and end
c
  990 continue
      return
      end

      subroutine checkshift (timref)
c
c  checkshift(timref)
c
c  Purpose:  apply large time shifts for specified periods
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  timref       R*4     I/O     size 3 reference time at start line
c                               of data: year, day, sec
c
c  By: F.S. Patt, SAIC, 1 October 2002
c
c  Notes:  
c
c  Modification History:
c
c  Added period in 2006 to correct time shift during GPS outage.
c  F.S. Patt, 17 April 2006
c
c  Added special case for 2010 with linear shift from days 57 to 88
c  F.S. Patt, 5 April 2010

      implicit none

      real*8 timref(3)
      real*8 r8jd, secinday
      integer*4 i, iy, id, jd, nshift

      parameter (secinday = 60.d0 * 60.d0 * 24.d0 )
      parameter (nshift = 5)
      
      real*8 sjd1(nshift), sjd2(nshift), shft(nshift)

c     Defined time shift periods and shifts
c       Start date/time
      data sjd1/2450995.700d0,         ! Day 181, 1998, 1648 UT
     *          2452544.625d0,         ! Day 269, 2002, 1500 UT
     *          2453807.74375d0,       ! Day 71, 2006, 1751 UT        
     *          2453985.850d0,         ! Day 249, 2006, 20:24 UT
     *          2455255.d0/           ! Day 57, 2010, 12:00 UT
c       End date/time
      data sjd2/2450996.650d0,         ! Day 182, 1998, 1536 UT
     *          2452546.125d0,         ! Day 271, 2002, 0300 UT
     *          2453809.075d0,         ! Day 73, 2006, 0148 UT
     *          2453986.200d0,         ! Day 250, 2006, 0448 UT 
     *          2455285.8d0/           ! Day 88, 2010, 19:12 UT
c       Shift in seconds
      data shft/8.d0,8.d0,-2.d0,-3.d0, 0.d0/
      
c     Convert reference time to Julian day
      iy = timref(1)
      id = timref(2)
      r8jd = jd(iy, 1, id) + timref(3)/secinday

c     Compute 2010 shift period as linear function of time
      shft(5) = (r8jd - sjd1(5))*0.41d0

c     Check against defined shift periods
      do i = 1,nshift

c       If reference time is within shift time range, apply shift
         if ((r8jd.gt.sjd1(i)).and.(r8jd.lt.sjd2(i))) then
            timref(3) = timref(3) + shft(i)
            write(*,*) 'Time shift applied ', shft(i)
         end if
      end do

      return
      end
