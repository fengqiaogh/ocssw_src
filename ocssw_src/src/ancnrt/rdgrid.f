      subroutine rdgrid(filename, type, ozone, julval, year,
     .   start, end, acntime, rval)
c ********************************************************************
c     Routine reads the daily gridded data from real-time
c      Meteor-3/TOMS, TOVS, EPTOMS or ADTOMS.
c
c     Program derived from TOMS project program RDGRID.FOR
c  filename  character  I  text file with ozone values
c  type      integer    I  type 0 = TOMS, 1 = TOVS, and 2 = EPTOMS
c                          3 = ADTOMS
c  ozone     integer*2  O  ozone array/ozone

c
c     Brian Schieber, SAIC/GSC 3/93
c     Mods:
c       BDS, 4/23/93 - ozone array set at integer*2
c     BDS, 9/29/95 - support reading start/end times from file.
c       TOVS files have 1 extra header line for Start/End times:
c       Start_Time: 1995123010101001 End_Time: 1995321010101001
c            *** Therefore TOMS format not currently supported *****
c     BDS, 9/3/96 - support reading new EPTOMS file
c                    - add "TYPE" (integer) argument for "TOMS" or "TOVS"
c     BDS, 9/12/96 - support three types: TOMS, TOVS, and EPTOMS
c                    - add Ascention time Metadata for EPTOMS data type.
c  Example header for each of the types:
c   TOMS:
c     [ Day: 100         April 10, 1993  Real Time Nimbus-7 TOMS]
c   TOVS:
c     [ Day: 274      October   1, 1995  Real Time NOAA TOVS]
c     [ Start_Time: 1995273183828000  End_Time: 1995275071153000]
c   EPTOMS:
c     [ Day: 247 Sep  3, 1996 Near Realtime     EP/TOMS    OZONE    Asc LECT: 11:16 AM]
c   ADTOMS:
c     [ Day: 255 Sep 11, 1996 Near Realtime    ADEOS/TOMS  OZONE    Asc LECT: 10:41 PM]
c ********************************************************************

c incoming args and local variables

      character   filename*(*)
      character   header*80
      integer     type
      integer*2   ozone
      integer   rval
      integer   julval
      integer   year
      character month*8
      character acntime*8
      character day*2
      character start*(*)
      character end*(*)
      dimension   ozone(360,180)
c
c     open the input file
c
      open(1,file=filename,status='old', err=998)
c
c     read 3 header lines, skipping last 2
c
c     type = 0, TOMS   type = 1, TOVS, type = 2, EPTOMS, type = 3,ADTOMS
      if (type .EQ. 0) then                   ! TOMS
        read(1, 10) julval, month, day, year
   10   format(6x, i3, 1x, a3, 1x, a2, 2x, i4)
        acntime = ""
      else if (type .EQ. 1) then                ! TOVS
        read(1, 15) julval, month, day, year
15      format(6x, i3, 6x, a8, 1x, a2, 2x, i4)
        read(1, 20) start, end
20      format(13x, a16, 12x, a16)
        acntime = ""
      else if ((type .EQ. 2) .OR. (type .EQ. 3)) then ! EPTOMS or ADTOMS
        read(1, 40) julval, month, day, year, acntime
40      format(6x, i3, 1x, a3, 1x, a2, 2x, i4, 49x, a8)
      else
        print *,
     1'Error in rdgrid.f, TYPE not "TOMS", "TOVS", "EPTOMS" or "ADTOMS"'
        rval = 1
        return
      end if

c      print *, 'julian ', julval
c      print *, 'month [', month,']'
c      print *, 'day [', day,']'
c      print *, 'year ', year

c   if (type .EQ. 1) then
c     print *, 'start ', start
c     print *, 'end   ', end
c   endif

      if ((type .EQ. 2) .OR. (type .EQ. 3)) then
c        print *,'acntime [',acntime,']'
        if( ( acntime(8:8) .NE. 'M' ) .AND. ( acntime(8:8) .NE. 'm' ) ) then
          rval = 1
          close(1)
          return
        endif
      endif

      read(1,'(a80)') header
      read(1,'(a80)') header
c     
c     read in the data into the array ozone
c     store in West to East, North to South order
c     (ASCII files in South to North so swap)
c     
      do 30 i=180,1,-1
        read(1,'(1x,25i3)') (ozone(j,i),j=1,360)
30    continue

      close(1)
      goto 999
  998 continue

      print *,'Error reading file ',filename
      rval = 1
      return

  999 continue
      rval = 0
      return
      end
