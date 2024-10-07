      subroutine navtlm(input,nframes,navctl,gaclac,tlm,tilts)

c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/navtlm.f,v 1.1 1995/01/17 23:02:29 seawifsd Exp seawifsd $
c $Log: navtlm.f,v $
c Revision 1.1  1995/01/17 23:02:29  seawifsd
c Initial revision
c
c
c  navtim(input,nframes,navctl,tlm)
c
c  Purpose: Fill navigation telemetry structure from input data structure
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  input        struct   I      input data structure containing SeaStar ID 
c                               and time tag, converted spacecraft and
c                               instrument telemetry 
c  nframes      I*4      I      number of minor frames of input data
c  gaclac       I*4      I      data type indicator (1=GAC)
c  tlm          struct   O      output navigation telemetry structure
c
c  By: Frederick S. Patt, GSC, 10 August 93
c
c  Notes:       
c
c  Instrument telemetry locations are as defined in the OSC L-band ICD.
c   Instrument telemetry occurs in minor frames 1 and 3 for LAC and 
c   for GAC lines 1, 3, 4, 6, 7, 9, 10, 12, 13 and 15 in major frame.

c  Spacecraft telemetry locations are currently assumed as follows:
c   Analog words 7 - 9:    Attitude angles yaw, roll, pitch
c   Analog words 10 - 15:  Sun sensor angles
c   Analog words 16 - 19:  Earth scanner angles (phase and width)
c   Discrete words 1 - 3:  Sun sensor status (active sensor flag)
c
c Spacecraft telemetry is contained in minor frame 1 for LAC and GAC
c       
c
c  Modification History:
c
c  June 24, 1994:  Added tilt structure to pass tilt flags and record 
c  ranges to L1A software.
c  Frederick S. Patt, GSC
c
c  January 5, 1995:  Corrected code to set tilt flags for fill frames.
c  Frederick S. Patt, GSC
c
c  March 30, 1995:  Modified code to accept GAC minor frame sequence
c  of 1, 3, 2.  Frederick S. Patt, GSC 
c
c  May 2, 1995:  Modified code to use actual tilt telemetry
c  Frederick S. Patt, GSC
c
c  Sept. 24, 1997:  Modified to unpack Sun sensor time tags.
c  Frederick S. Patt, GSC
c
c  October 30, 1997:  Increased tilt angle tolerance to allow test with
c  new static tilt values.  Frederick S. Patt, GSC
c
c  December 1, 1997:  Corrected tilt state logic to check for missing tilt
c  telemetry at the start of a tilt change.  Frederick S. Patt, GSC
c
c  December 4, 1997:  Limited tilt states to 20, to avoid overflow of
c  tilts structure when processing corrupted telemetry.  Also changed
c  logic to ignore indeterminate tilt states, and therefor assume the
c  actual state was equal to the previous state.
c  B. A. Franz, GSC
c
c  February 25, 1998: Modified tilt-state logic to incorporate state 
c  filtering for errors in tilt telemetry. B. A. Franz, GSC
c 
c  July 14, 1998: Modified tilt consistency checking logic to extend tilt
c  change state through any adjacent missing or flagged values.
c  F. S. Patt, SAIC GSC
c
c  August 2, 1998:  Included flagging of tilt telemetry for single-line tilt
c  states; added a limit check for tilt motor angles.  F. S. Patt, SAIC GSC
c
c  September 22, 1998:  Modified logic for checking attitude sensor delta times
c  to reject positive values (since sensor times can't be later than frame 
c  times).  F. S. Patt, SAIC GSC

         implicit none
#include "input_s.fin"
#include "navctl_s.fin"
#include "tlm_str.fin"
#include "tilt_s.fin"
      type(input_struct) :: input(maxlin)
      type(navctl_struct) :: navctl
      type(tlm_struct) :: tlm
      type(tilt_states) :: tilts

      real*4 ttilt, tilttol, tanglim(2)
      integer*4 nframes, gaclac, nper, i, j, k, ip, mf, ilin
      integer*4 itilt, tflag, ntlts, istime
      byte stime(4)

      integer*4 tilt_flags(maxlin)
      integer*4 nlines
      logical gottwo
      equivalence (istime,stime)
      data tilttol/0.25/, tanglim/48.,225./

      tlm%ntlm = nframes
      ntlts = 0

c  Determine if LAC or GAC
      if (input(1)%sc_id(2).eq.15) then
        gaclac = 1
        nper = 5
      else
        gaclac = 0
        nper = 1
      end if

c     !
c     ! Initialize tilt flag array
c     !
      nlines = nframes*nper
      do i=1,nlines
          tilt_flags(i) = -1
      enddo

c     !
c     ! Extract telemetry from each frame
c     !
      do i=1, nframes

c  Get minor frame number
        mf = input(i)%sc_id(1)/128

c  If not fill frame    
        if (input(i)%flag.eq.0) then
        
c  Unpack tilt telemetry 
          do j=1,nper
            ilin = mf - 1 + j
c  If second of 3 scan lines then no tilt telemetry
            if (mod(ilin,3).eq.2) then
              tlm%tilt(1)%flag(j,i) = 1
              tlm%tilt(2)%flag(j,i) = 1
              tlm%tilt(1)%ang(j,i) = 0.
              tlm%tilt(2)%ang(j,i) = 0.
              tilt_flags((i-1)*nper+j) = -1
            else
c  Else get tilt telemetry 
c  Convert tilt motor angles to tilt angle
              call conv_tilt(input(i)%inst_ana(23,j),
     *            input(i)%inst_ana(22,j),ttilt)
c   Check tilt bits
              itilt = input(i)%inst_dis(16,j)*4 
     *          + input(i)%inst_dis(17,j)*2 + input(i)%inst_dis(18,j)
              if (itilt.eq.7) then
c  Tilt change in progress; get tilt from analog telemetry 
                 if ((input(i)%inst_ana(22,j).gt.tanglim(1)).and.
     *                (input(i)%inst_ana(22,j).lt.tanglim(2)).and.
     *                (input(i)%inst_ana(23,j).gt.tanglim(1)).and.
     *                (input(i)%inst_ana(23,j).lt.tanglim(2))) then
                    tlm%tilt(1)%flag(j,i) = 2
                    tlm%tilt(2)%flag(j,i) = 2
c               tlm%tilt(1)%ang(j,i) = input(i)%inst_ana(22,j)
c               tlm%tilt(2)%ang(j,i) = input(i)%inst_ana(23,j)
                    tlm%tilt(1)%ang(j,i) = ttilt
                    tlm%tilt(2)%ang(j,i) = ttilt
                    tflag = 3
                 else
                    tlm%tilt(1)%flag(j,i) = 1
                    tlm%tilt(2)%flag(j,i) = 1
                    tlm%tilt(1)%ang(j,i) = 0.
                    tlm%tilt(2)%ang(j,i) = 0.
                    tflag = -1
                 end if
 

              else if (itilt.eq.3) then
c  Tilt nadir aligned
                tlm%tilt(1)%flag(j,i) = 0
                tlm%tilt(2)%flag(j,i) = 0
                tlm%tilt(1)%ang(j,i) = 0.
                tlm%tilt(2)%ang(j,i) = 0.
                tflag = 0
                
              else if (itilt.eq.5) then
c  Tilt aft aligned
                tlm%tilt(1)%flag(j,i) = 0
                tlm%tilt(2)%flag(j,i) = 0
                tlm%tilt(1)%ang(j,i) = navctl%tiltaft
                tlm%tilt(2)%ang(j,i) = navctl%tiltaft
                tflag = 2

              else if (itilt.eq.6) then
c  Tilt forward aligned
                tlm%tilt(1)%flag(j,i) = 0
                tlm%tilt(2)%flag(j,i) = 0
                tlm%tilt(1)%ang(j,i) = navctl%tiltfor
                tlm%tilt(2)%ang(j,i) = navctl%tiltfor
                tflag = 1

              else
c  Tilt status unknown
                tlm%tilt(1)%flag(j,i) = 1
                tlm%tilt(2)%flag(j,i) = 1
                tlm%tilt(1)%ang(j,i) = 0.
                tlm%tilt(2)%ang(j,i) = 0.
                tflag = -1
              end if

c  Check for inconsistency between tilt bits and angles
              if (abs(ttilt-tlm%tilt(1)%ang(j,i)) .gt. tilttol) then
                tlm%tilt(1)%flag(j,i) = 1
                tlm%tilt(2)%flag(j,i) = 1
                tlm%tilt(1)%ang(j,i) = 0.
                tlm%tilt(2)%ang(j,i) = 0.
                tflag = -1
              end if

              tilt_flags((i-1)*nper+j) = tflag

            end if    ! inst telemetry frame

          end do

c  Else if fill record, set flags
        else
          do j=1,nper
              tlm%tilt(1)%flag(j,i) = 1
              tlm%tilt(2)%flag(j,i) = 1
              tlm%tilt(1)%ang(j,i) = 0.
              tlm%tilt(2)%ang(j,i) = 0.
          end do
        end if

c  If minor frame 1, unpack ACS telemetry
        if ((input(i)%flag.eq.0).and.(mf.eq.1)) then

c  Unpack attitude angles
            do j=1,3
              tlm%sc_att%att(j,i) = input(i)%sc_ana(j+6)
            end do
            tlm%sc_att%flag(i) = 0

c  Unpack Sun sensor angles and status 
            do j=1,3 
              if (input(i)%sc_dis(j).eq.3) then
                tlm%sun(j)%active(i) = 0
              else
                tlm%sun(j)%active(i) = 1
              end if            
              tlm%sun(j)%flag(i) = 0
c             tlm%sun(j)%deltim(i) = 0.0
              tlm%sun(j)%ang(1,i) = input(i)%sc_ana(2*j+8)
              tlm%sun(j)%ang(2,i) = input(i)%sc_ana(2*j+9)
c  Get delta time
#ifdef LINUX
              do k=1,4
                 stime(k)=input(i)%sc_dis(16+4*j+5-k)
              end do
#else
              do k=1,4
                 stime(k)=input(i)%sc_dis(16+4*j+k)
              end do
#endif
              istime = mod(istime,86400000)
              tlm%sun(j)%deltim(i) = (istime - input(i)%msec)/1.d3
              if ((tlm%sun(j)%deltim(i).gt.0.).or.
     *             (tlm%sun(j)%deltim(i).lt.-10.))
     *             tlm%sun(j)%flag(i) = 1
            end do

c  Unpack Earth scanner angles
            do j=1,2
              tlm%earth(j)%active(i) = 0
              tlm%earth(j)%flag(i) = 0
              tlm%earth(j)%deltim(i) = 0.0
              tlm%earth(j)%widphse(1,i) = input(i)%sc_ana(2*j+15)
              tlm%earth(j)%widphse(2,i) = input(i)%sc_ana(2*j+14)
              if (tlm%earth(j)%widphse(2,i).gt.180.) 
     *  tlm%earth(j)%widphse(2,i) = tlm%earth(j)%widphse(2,i) - 360.0
#ifdef LINUX
              do k=1,4
                 stime(k)=input(i)%sc_dis(28+4*j+5-k)
              end do
#else
              do k=1,4
                 stime(k)=input(i)%sc_dis(28+4*j+k)
              end do
#endif
              istime = mod(istime,86400000)
              tlm%earth(j)%deltim(i) = (istime - input(i)%msec)/1.d3
              if ((tlm%earth(j)%deltim(i).gt.0.).or.
     *             (tlm%earth(j)%deltim(i).lt.-10.))
     *             tlm%earth(j)%flag(i) = 1
           end do

        else

c  Set Attitude flags
            do j=1,3
              tlm%sc_att%att(j,i) = 0.
            end do
            tlm%sc_att%flag(i) = 1

c  Set Sun sensor angles, flags and status 
            do j=1,3 
              tlm%sun(j)%active(i) = 1
              tlm%sun(j)%flag(i) = 1
              tlm%sun(j)%deltim(i) = 0.
              tlm%sun(j)%ang(1,i) = 0.
              tlm%sun(j)%ang(2,i) = 0.
            end do

c  Set Earth scanner angles
            do j=1,2
              tlm%earth(j)%active(i) = 1
              tlm%earth(j)%flag(i) = 1
              tlm%earth(j)%deltim(i) = 0.0
              tlm%earth(j)%widphse(1,i) = 0.
              tlm%earth(j)%widphse(2,i) = 0.
            end do  
          
        end if
      end do

c     !
c     ! Filter-out single line tilt states
c     !

c     Find first unflagged value
      ip = 1
      do while (tilt_flags(ip) .eq. -1)
         ip = ip + 1
      end do

      gottwo = .false.

      do i = ip+1, nlines

c     Find next unflagged value
         if (tilt_flags(i) .ne. -1) then
            if (tilt_flags(i) .eq. tilt_flags(ip)) then
               gottwo = .true.
            else
               if (.not.gottwo) then
                  tilt_flags(ip) = -1
c        Also flag tilt telemetry for this line
                  j = mod( (ip-1), nper) + 1
                  k = (ip-1)/nper + 1
                  tlm%tilt(1)%flag(j,k) = 1
                  tlm%tilt(2)%flag(j,k) = 1
               end if
               gottwo = .false.
            end if
            ip = i
         end if
      end do
      if (.not.gottwo) then 
         tilt_flags(ip) = -1
         j = mod( (ip-1), nper) + 1
         k = (ip-1)/nper + 1
         tlm%tilt(1)%flag(j,k) = 1
         tlm%tilt(2)%flag(j,k) = 1
      end if
c     !
c     ! Now fill in missing tilt flags
c     !

      ip = 1
      do while (tilt_flags(ip) .eq. -1)
         ip = ip + 1
      end do

      do i=1,ip
         tilt_flags(i) = tilt_flags(ip)
      end do

      do i = ip+1, nlines

c     Find next unflagged value
         if (tilt_flags(i) .ne. -1) then

c       Fill missing values if needed
            if ((i - ip) .gt. 1) then

c             Use max value for flag to default to tilt change state
               tflag = max(tilt_flags(ip),tilt_flags(i))
               do j=ip+1, i-1
                  tilt_flags(j) = tflag
               end do
            end if
            ip = i
         end if
      end do

      do i=ip, nlines
         tilt_flags(i) = tilt_flags(ip)
      end do

c     !
c     ! Locate and store tilt-change periods
c     !
      ntlts = 1
      tilts%tilt_ranges(1,ntlts) = 1
      tilts%tilt_flags(ntlts)    = tilt_flags(1)
c
      do i=2,nlines

          if (tilt_flags(i) .ne. tilt_flags(i-1)) then 
              tilts%tilt_ranges(2,ntlts) = i-1 
              if (ntlts .eq. 20) then
                  write(*,*)
                  write(*,*) "navtlm: Error: excessive tilt changes."
                  write(*,*)
                  goto 10
              endif
              ntlts = ntlts+1
              tilts%tilt_ranges(1,ntlts) = i
              tilts%tilt_flags(ntlts)    = tilt_flags(i)
          endif

      enddo
c
  10  continue 
      tilts%tilt_ranges(2,ntlts) = nper*nframes
      tilts%ntilts = ntlts
      
              
  990 continue
      return
      end
