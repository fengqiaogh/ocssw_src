      subroutine read_gps(input,nframes,scal_p,xmaglm,gpsvec,nsig,
     *  igyr,igday,gpsec,secst,secend,ngps)
c
c  Purpose:  This subroutine extracts the valid GPS data from the input
c               telemetry structure.  Only data points with four or more
c               GPS signals are included.  Also, the time tags for all
c               input minor frames are checked to determine the time 
c               range of the data.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  input        struct   I      input data structure
c  nframes      I*4      I      Number of frames in input structure
c  scal_p       R*8      I      Scale factor for position vectors
c  xmaglm(2)    R*8      I      Limits on position magnitude
c  gpsvec(6,*)  R*8      O      Output GPS orbit vectors
c  nsig(*)      I*4      O      Number of GPS signals for each vector
c  igyr         I*4      O      Time tag year
c  igday        I*4      O      Time tag day-of-year
c  gpsec(*)     R*8      O      Time tag seconds of igday
c  secend       R*8      O      End of interval time in seconds
c  ngps         I*4      O      Number of GPS vectors
c
c  By: Frederick S. Patt, GSC, December 21, 1993
c
c  Notes:  
c
c  Modification History:
c
c  Modified to unpack GPS time tag from telemetry (previously used
c  minor frame time tag).  F. S. Patt, GSC, September 14, 1996
c
c  Modified to perform sanity check on time tags (e.g., all tags 
c  fall within +/- 1 day of first valid time tag.
c  F. S. Patt, GSC, October 4, 1996.
c
c  Fixed a bug in validating the first GPS time code.  
c  F. S. Patt, SAIC GSC, September 7, 1997.
c
c  Added additional range checking to the seconds-of-day.
c  F. S. Patt, SAIC GSC, September 10, 1997
c
c  Changed integer*1 to byte for Sun OS compatibility, B. A. Franz,
c  GSC, November 14, 1997.
c
c  Modified to reject GPS samples with time tags later than minor frame.
c  F. S. Patt, SAIC GSC, July 22, 1998.
c
c  Modified to reject GPS samples prior to current time range.
c  F. S. Patt, SAIC GSC, March 3, 2000
c
c  Modified to remove checks on minor frame times (SWl01 takes care of this).
c  F. S. Patt, SAIC GSC, May 20, 2000
c
c  Fixed bug which caused error message to be output for duplicate GPS sample.
c  F. S. Patt, SAIC GSC, June 8, 2000 
c
c  Fixed bug introduced in indexing of input array during generation of 
c  single-source code for SGI and Linux. 
c  F. S. Patt, SAIC GSC, January 19, 2001 
c
c  Modified duplicate time check to use a tolerance instead of equivalence.
c  F. S. Patt, SAIC, February 17, 2004

      implicit none
#include "nav_cnst.fin"
#include "input_s.fin"

      type(input_struct) :: input(maxlin)

      real*8 gpsvec(6,maxlin), gpsec(maxlin), scal_p, secst, secend
      real*8 sec, xmaglm(2), pmag
      integer*4 nsig(maxlin), nframes, igyr, igday, ngps
      integer*4 iframe, mframe, if1, j, iy2, im2, id2, ms, mp 
      integer*4 jd, jd1, jd2, jdtol
      integer*2 ib2
      byte      b2(2)
      logical gotone,gotime
      data jdtol/1/
      equivalence (ib2,b2)

c  Find two good minor frames with the same day in the time tag
      if1 = 1
      ngps = 0
          gotime = .true.
          igyr = input(if1)%iyear
          igday = input(if1)%iday
          jd1 = jd(igyr,1,igday)
          secst = input(if1)%msec/1000.d0
          secend = secst

c  Find first good vector
      gotone = .false. 
      dowhile (.not.gotone)
        mframe = input(if1)%sc_id(1)/128

c  Check minor frame time tag against time range
        iy2 = input(if1)%iyear
        id2 = input(if1)%iday
        jd2 = jd(iy2,1,id2)
        sec = input(if1)%msec/1000.d0 
          sec = sec +(jd2-jd1)*864.d2
          if (sec.lt.secst) secst = sec
          if (sec.gt.secend) secend = sec
        if((input(if1)%flag.eq.0).and.(mframe.eq.1)) then

c  Unpack number of GPS signals from telemetry data and check for minimum
          ms = input(if1)%sc_dis(12)
          if (ms.ge.4) then

c  Get GPS time tag
c   Year is two bytes
#ifdef LINUX
            b2(1) = input(if1)%sc_dis(15)
            b2(2) = input(if1)%sc_dis(14)
#else
            b2(1) = input(if1)%sc_dis(14)
            b2(2) = input(if1)%sc_dis(15)
#endif
            iy2 = ib2
c   Month and day 
            im2 = input(if1)%sc_dis(16)
            id2 = input(if1)%sc_dis(17)
            jd2 = jd(iy2,im2,id2)
c  Convert hours, minutes, seconds, fraction to seconds of day
            sec = input(if1)%sc_dis(18)*36.d2 
     *           + input(if1)%sc_dis(19)*6.d1 + input(if1)%sc_dis(20)
     *           + input(if1)%sc_ana(21) 

            sec = sec + (jd2-jd1)*864.d2

            
            if ( (sec.ge.(secst-10.)) .and. (sec.lt.secend) ) then

              ngps = ngps + 1
              gotone = .true.
              if (sec.lt.secst) secst = sec

c  Store data in output array
              gpsec(ngps) = sec
              nsig(ngps) = ms
              do j=1,3
                gpsvec(j,ngps) = input(if1)%sc_ana(j)*scal_p
                gpsvec(j+3,ngps) = input(if1)%sc_ana(j+3)
              end do
c Perform limit check on position magnitude
              pmag = sqrt(gpsvec(1,ngps)**2 + gpsvec(2,ngps)**2
     *          + gpsvec(3,ngps)**2)
              if ((pmag.lt.xmaglm(1)).or.(pmag.gt.xmaglm(2))) then
                gotone = .false.
                ngps = ngps - 1  
              end if
            end if
          end if
        end if
        if1 = if1 + 1
        if (if1.gt.nframes) go to 999
      end do
      mp = ms

c  Now unpack remaining GPS vectors; check for >4 signals at previous point.
      do iframe = if1,nframes
        iy2 = input(iframe)%iyear
        id2 = input(iframe)%iday
        jd2 = jd(iy2,1,id2)
        sec = input(iframe)%msec/1000.d0 
        sec = sec + (jd2-jd1)*864.d2

c  Check time tag against range
        if (sec.lt.secst) secst = sec
        if (sec.gt.secend) secend = sec

c  Check for minor frame 1 and valid data
        mframe = input(iframe)%sc_id(1)/128
        if ((input(iframe)%flag.eq.0).and.(mframe.eq.1)) then

c  Unpack number of GPS signals from telemetry data
          ms = input(iframe)%sc_dis(12)

c  Get GPS time tag
c   Year is two bytes
#ifdef LINUX
            b2(1) = input(iframe)%sc_dis(15)
            b2(2) = input(iframe)%sc_dis(14)
#else
            b2(1) = input(iframe)%sc_dis(14)
            b2(2) = input(iframe)%sc_dis(15)
#endif
          iy2 = ib2
c   Month and day 
          im2 = input(iframe)%sc_dis(16)
          id2 = input(iframe)%sc_dis(17)
          jd2 = jd(iy2,im2,id2)
          sec = input(iframe)%sc_dis(18)*36.d2 
     *         + input(iframe)%sc_dis(19)*6.d1 
     *         + input(iframe)%sc_dis(20)
     *         + input(iframe)%sc_ana(21) 

c  Convert hours, minutes, seconds, fraction to seconds of day
          sec = sec + (jd2-jd1)*864.d2

c  Check for time tag out of range
          if ((sec.ge.secst).and.(sec.le.secend)) then

c  Check for < minimum GPS signals or duplicated time tag
             if ((ms.ge.4).and.(mp.ge.4).and.
     *            (dabs(sec-gpsec(ngps)).gt.1.d-9)) then
                ngps = ngps + 1
             
c     Store data in output array
                gpsec(ngps) = sec
                nsig(ngps) = ms
                do j=1,3
                   gpsvec(j,ngps) = input(iframe)%sc_ana(j)*scal_p
                   gpsvec(j+3,ngps) = input(iframe)%sc_ana(j+3)
                end do
c     Perform limit check on position magnitude
                pmag = sqrt(gpsvec(1,ngps)**2 + gpsvec(2,ngps)**2
     *               + gpsvec(3,ngps)**2)
                if ((pmag.lt.xmaglm(1)).or.(pmag.gt.xmaglm(2))) then
                   ngps = ngps - 1  
                end if
             
c  Save number of GPS signals
                mp = ms
             end if
          else
             print *,'READ_GPS:  Invalid date',iy2,im2,id2,
     *            ' record',iframe
          end if
        end if
      end do
      
  999 print *,'READ_GPS:',ngps,'  valid GPS vectors'
      print *,' starting at ',igyr,igday,secst
      print *,' ending   at ',igyr,igday,secend

      return
      end
