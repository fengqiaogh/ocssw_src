      subroutine tiltcomp( nlines, tlm, timref, time, gaclac, navqc, 
     1  tiltpr, tiltfl )

c $Header$
c $Log$
c
c  tiltcomp( nlines, tlm, gaclac, navqc, tiltpr, tiltfl )
c
c  Purpose: process tilt angle data to a line-by-line array
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nlines       I*4      I      number of scan lines covered
c  tlm          struct   I      telemetry structure containing tilt info
c  timref       R*8      I      size 3 reference time at start line
c                               of data: year, day, sec
c  time         R*8      I      array of time in seconds relative to
c                               timref for every scan line
c  gaclac       I*4      I      flag for GAC or LAC data.  If LAC
c                               data, a time exists for each scan
c                               line, else only once every 5 lines
c  navqc        struct   I      navigation quality control info
c  tiltpr       R*4      O      size nlines processed tilt data 
c                               for each line and for two tilt measurements
c  tiltfl       I*4      O      flags for goodness of tilts 0 - good,
c                               1 - bad
c
c  By: W. Robinson, GSC, 13 Apr 93
c
c  Notes:  
c
c  Modification History:
c
c  Fixed error in last call to runfit3t.  F. S. Patt, January 12, 1995.
c
c  Fixed error in indexing of smoothed tilt angles.  F.S. Patt, Oct. 31, 1997 
c
c  Fixed bug in initializing counter for end-of-scene processing.  F. S. Patt,
c  SAIC GSC, May 21, 1998.
c  
c  Added a consistency check on the tilt angles during tilt changes; include 
c  a final limit check on tilt angles to catch processing problems.
c  F. S. Patt, SAIC GSC, August 4, 1998.

 
      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
      type(tlm_struct) :: tlm
      type(navqc_struct) :: navqc
c
      integer*4 gaclac, nlines, tiltfl(nlines)
c
      integer*4 nper, ilin, j, k, i1, in(4), ntilt, n1, n2, tfl(100)
      integer*4 nskip, iret, nfit, jp, sgn
      real*8 timref(3),time(nlines)
      real*4 tiltpr(nlines)
      real*4 tout(100),tlcl(100),tdif(100),tdlim(2),tldif,tlim(2)
      logical gottwo
      data nfit/9/,nskip/2/,tdif/100*0.0/,tdlim/-1.,10./,tlim/-20.,20./
c
c
c       initialize output arrays for processed tilt
c
      nper = 1
      if( gaclac .eq. 1 ) nper = 5
c
      do ilin = 1, nlines
         tiltpr(ilin) = 0
         tiltfl(ilin) = 1
      end do
c
c       First transfer the tilt info
c
      do ilin = 1, tlm%ntlm
          
         do j = 1, nper
          
c
            tiltfl((ilin - 1) * nper + j ) = tlm%tilt(1)%flag(j,ilin)
c            tiltfl(2,(ilin - 1) * nper + j ) = tlm%tilt(2)%flag(j,ilin)
c
            tiltpr((ilin - 1) * nper + j ) = tlm%tilt(1)%ang(j,ilin)
c            tiltpr(2,(ilin - 1) * nper + j ) = tlm%tilt(2)%ang(j,ilin)
         end do
      end do

c       Add consistency check at some future date
c     call tiltcnst(

c       Now fill in missing values

c      do itilt=1,2
      ilin = 1
      i1 = 1

c       Initialize by finding first 4 static tilt times
      do k=1,4
        call fndflg(tiltfl(1),nlines,i1,in(k))
c       If end of range encountered, process data and exit
        if (in(k).eq.-1) then

c       If 0 or 1 static values, attempt smoothing
          if (k.lt.3) then
            print *,'Insufficient static tilt data samples'
            print *,' Smoothing attempted'
            do j=1,nlines
              if (tiltfl(j).eq.2) tiltfl(j) = 0
            end do
            ntilt = nlines
            call runfit3t(nfit,nskip,tlcl,tiltpr(1),ntilt,1,tiltfl(1),
     *        1,time(1),tdif,tout,tfl,iret)
            do j=1,nlines
               tiltpr(j) = tout(j)
               tiltfl(j) = 2 - tfl(j)
            end do
            go to 990

c       Else if 2 static values, process data
          else if (k.eq.3) then
            in(3) = in(2)
          end if
        end if
      i1 = in(k) + 1
      end do

c       Now handle cases of non-static tilt at start of interval
c        Change tilt change data flags to be used in smoothing
      n2 = 0
      n1 = 0
      do j=1,in(3)-1
        if (tiltfl(j).eq.2) then
          tiltfl(j) = 0
          n2 = n2 + 1
          if (j.lt.in(1)) n1 = n1 + 1
        end if
      end do

c       Try to perform smoothing  
      ntilt = in(3)
        call runfit3t(nfit,nskip,tlcl,tiltpr(1),ntilt,1,tiltfl(1),
     *    1,time(1),tdif,tout,tfl,iret)

c       Fill in values before first pointer if necessary
      if (in(1).gt.1) then
         if ((n1.gt.0).and.(iret.eq.0)) then
            do j=1,in(1)-1
               tiltpr(j) = tout(j)
               tiltfl(j) = 2 - tfl(j)
            end do
         else
            do j=1,in(1)-1
               tiltpr(j) = tiltpr(in(1))
               tiltfl(j) = tiltfl(in(1))
            end do
         end if
      end if

c       Fill in values between first and second pointers
      if (in(2).gt.(in(1)+1)) then
        if (tiltpr(in(2)).eq.tiltpr(in(1))) then
          do j=in(1)+1,in(2)-1
            tiltpr(j) = tiltpr(in(1))
            tiltfl(j) = tiltfl(in(1))
          end do
        else
          do j=in(1)+1,in(2)-1
            tiltpr(j) = tout(j)
            tiltfl(j) = 2 - tfl(j)
          end do
        end if
      end if


c       Now fill remaining array locations with good values
      dowhile (in(4).ne.-1)

c       Check for tilt change between second and third pointers
        if (tiltpr(in(3)).eq.tiltpr(in(2))) then
          if (in(3).gt.(in(2)+1)) then
            do j=in(2)+1,in(3)-1
              tiltpr(j) = tiltpr(in(2))
              tiltfl(j) = tiltfl(in(2))
            end do
          end if
          
c       Else smooth and load tilt change data
        else

c       Perform consistency check on tilt angles
           gottwo = .false.
           jp = in(2)
           sgn = sign(1.0,(tiltpr(in(3)) - tiltpr(in(2))))
           do j=in(2)+1,in(3)-1
              if (tiltfl(j).eq.2) then
                 if (time(j).ne.time(jp)) then
                    tldif = sgn*(tiltpr(j) - tiltpr(jp)) / 
     *                   (time(j)-time(jp))
                    if ((tldif.lt.tdlim(1)).or.(tldif.gt.tdlim(2))) then
                       if (.not.gottwo) tiltfl(jp) = 1
                       gottwo = .false.
                    else                       
                       tiltfl(j) = 0
                       gottwo = .true.
                    end if
                    jp = j
                 else
                    tiltfl(j) = 1
                 end if
              end if
           end do
           if (.not.gottwo) tiltfl(jp) = 1
           ntilt = in(3) - in(2) + 1
           call runfit3t(nfit,nskip,tlcl,tiltpr(in(2)),ntilt,1,
     1          tiltfl(in(2)),1,time(in(2)),tdif,tout,tfl,iret)
           do j=in(2)+1,in(3)-1
              tiltpr(j) = tout(j-in(2)+1)
              tiltfl(j) = 2 - tfl(j-in(2)+1)
           end do
        end if
        
c       Find next static value and shift pointers
        do j=1,3
          in(j) = in(j+1)
        end do
        i1 = in(3) + 1
        call fndflg(tiltfl(1),nlines,i1,in(4))

c       End of main processing loop
      end do
        
c       Now process data at end of interval
c        Change tilt change data flags to be used in smoothing
      n1 = 0
      n2 = 0
      do j=in(2),nlines
        if (tiltfl(j).eq.2) then
          tiltfl(j) = 0
          n2 = n2 + 1
          if (j.gt.in(3)) n1 = n1 + 1
        end if
      end do

c       Try to perform smoothing  
      ntilt = nlines - in(1) + 1
        call runfit3t(nfit,nskip,tlcl,tiltpr(in(1)),ntilt,1,
     1    tiltfl(in(1)),1,time(in(1)),tdif,tout,tfl,iret)

c       Fill in values between second and third pointers
      if (in(3).gt.(in(2)+1)) then
        if (tiltpr(in(3)).eq.tiltpr(in(2))) then
          do j=in(2)+1,in(3)-1
            tiltpr(j) = tiltpr(in(2))
            tiltfl(j) = tiltfl(in(2))
          end do
        else
          do j=in(2)+1,in(3)-1
            tiltpr(j) = tout(j-in(1)+1)
            tiltfl(j) = 2 - tfl(j-in(1)+1)
          end do
        end if
      end if

c       Fill in values after third pointer if necessary
      if (in(3).lt.nlines) then
        if ((n1.gt.0).and.(iret.eq.0)) then
          do j=in(3)+1,nlines
            tiltpr(j) = tout(j-in(1)+1)
            tiltfl(j) = 2 - tfl(j-in(1)+1)
          end do
        else
          do j=in(3)+1,nlines
            tiltpr(j) = tiltpr(in(3))
            tiltfl(j) = tiltfl(in(3))
          end do
        end if
      end if

c  
c       and end
c
  990 continue

      do j=1,nlines
c       Perform final check to limit tilt angles to +/- 20 degrees
         if (tiltpr(j).gt.tlim(2)) then
            tiltpr(j) = tlim(2)
            tiltfl(j) = 1
         else if (tiltpr(j).lt.tlim(1)) then
            tiltpr(j) = tlim(1)
            tiltfl(j) = 1
         end if
         
c       print *,tiltpr(j),tiltfl(j)
      end do

      return
      end
