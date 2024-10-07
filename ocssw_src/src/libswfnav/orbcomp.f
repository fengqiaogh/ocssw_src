        subroutine orbcomp(input, nframes, orbit, ier)

c $Header$
c $Log$
c
c  Purpose:  This routine performs filtering of the GPS data for 
c       subsequent navigation processing.  It unpacks the GPS data from 
c       the converted spacecraft telemetry and fits the valid data to an 
c       orbit model.  The filtered vectors are stored at 1-minute 
c       intervals in the structure  ORBIT.  Filtering methods are 
c       described in the document TBD.
c
c       The orbit filtering routines information from three external
c       files:  elements.dat, a direct-access binary file which stores 
c                the mean element sets used to initialize the orbit 
c                integrator and save the results of the orbit filtering;
c               asap_parms.dat, which contains the gravitational model 
c                terms and other parameters needed by ASAP to integrate 
c                the orbit; and
c               orbctl.nl, which contains user-specifiable parameters 
c                in FORTRAN namelist format (see below).
c       All files are assumed to be located in the directory specified 
c       by the environmental variable ORBCTL.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  input        struct   I      Input data structure
c  nframes      I*4      I      Number of frames in input structure
c  orbit        struct   I      Output orbit data structure
c  ier          I*4      I      Error code:     =0, success
c
c
c  Parameters in namelist /orbctl/:
c
c  Name         Type   Default          Description
c  --------     ----   -------          -----------
c  updtol(6)    R*8     0.002, 1.d-5,   Tolerance for iteration of orbital      
c                       0.001, 0.001,    element updates
c                       0.1, 0.001
c  s0(6)        R*8     1.d4,2*4.d6,    State weights for least-squares
c                       3*4.d5           determination of element updates
c  gps_scal_p   R*8     1.0             Scale factor for GPS position vectors
c  fit          L*4     .true.          Use fitted elements if available
c  iplt         I*4     0               If >0, write ASAP output to file
c  pmaglm(2)    R*8     7070., 7110     Validation limits for checking
c                                        orbit position magnitude
c  li, mi       I*4     21, 12          Order and degree of gravity model
c  ndmax        I*4     7               Maximum days to propagate previous
c                                        orbital elements
c
c  By: Frederick S. Patt, GSC, December 21, 1993
c
c  Notes:  
c
c  Modification History:
c
c  Modified to allow for processing of non-contiguous data:  added calls
c  to add_elements and vec2mean, and changes logic of calls to put_elements.
c  F.S. Patt, GSC, October 26, 1994.
c
c  Modified to limit maximum updates to elements per iteration.
c  F. S. Patt, GSC, January 28, 1996.
c
c  Modified to add limit checks to orbit position magnitude
c  F. S. Patt, GSC, May 8, 1996.
c
c  Modified to limit number of iterations on GPS fitting.
c  F. S. Patt, SAIC GSC, July 22, 1998.
c
c  Modified to return error if initial elements are more than N days 
c  prior to data, where N is a namelist parameter defaulting to 7.  
c  F. S. Patt, SAIC GSC, January 11, 2001.
c
c  Modified to force GPS fitting for multi-orbit GAC files.
c  F. S. Patt, SAIC, April 2, 2003.
c
c  Modified to update error message for stale elements.dat file to 
c  reflect end of NASA data collection on December 23, 2004.
c  F. S. Patt, SAIC, January 5, 2005.

      implicit none
#include "nav_cnst.fin"
#include "input_s.fin"
#include "orbit_s.fin"

      type(input_struct) :: input(maxlin)
      type(orbit_struct) :: orbit

      real*8 gpsvec(6,maxlin),gpsec(maxlin),asap(6,maxlin),tsap(maxlin)
      real*8 vecs(6,maxlin),driv(6,3,maxlin)
      real*8 orbinit(6),orbupd(6),orbend(6),updorb(6),maxupd(6)
      real*8 secinit,secst,secend,cdrg,ge,aj2,tdif,oneorb
      real*8 updtol(6),s0(6),gps_scal_p,tdifmax,wdifac,pmaglm(2)
      integer*4 nsig(maxlin),nframes,igyr,igday,ngps,jd,i,j,ilast,ir
      integer*4 iyinit,idinit,irec,iddif,nstp,nstr,niter,maxit
      integer*4 lun,li,mi,iplt,ier,nmin,ndmax
      logical iter,fit,write,init
      character*80 filnm
C  Namelist parameters
      namelist /orbctl/updtol,s0,gps_scal_p,fit,nmin,iplt,pmaglm,li,mi,
     *     ndmax
      data ge/3.9860050D5/,aj2/0.10826270E-02/
      data li/21/,mi/21/,iplt/0/,nmin/3/,lun/17/
      data updtol/0.002d0,1.d-5,1.d-3,1.d-3,0.1d0,1.d-3/
      data s0/1.d4,2*4.d6,3*4.d5/
      data gps_scal_p/1.d0/,tdifmax/86400.d0/,ndmax/7/
      data pmaglm/7070.d0, 7110.d0/
      data fit/.true./
      data maxit/10/
      data oneorb/5940.d0/

C  Open file and read namelist
      filnm = '$ORBCTL/orbctl.nl'
      call filenv(filnm,filnm)
      open (18,file=filnm,status='old')
      read (18,nml=orbctl)       
      write (*,nml=orbctl)
      close (18)

C  Set maximum update per iteration to 500 times tolerance 
C   (a SWAG, may get more sophisticated in the future)
      do i=1,6
        maxupd(i) = 500.* updtol(i)
      end do    
     
C  Unpack GPS data 
      call read_gps(input,nframes,gps_scal_p,pmaglm,gpsvec,nsig,
     *  igyr,igday,gpsec,secst,secend,ngps)

C  If multi-orbit GAC data, force GPS fitting
      if ((input(1)%sc_id(2).eq.15).and.((secend-secst).gt.oneorb)) then
         fit = .false.
      end if

C  If less than minimum number of GPS, use fitted elements
      if (ngps.lt.nmin) fit = .true.

C  Get initial elements from last processing
      call get_elements(igyr,igday,secst,secend,fit,orbinit,cdrg,
     *  iyinit,idinit,secinit,irec,ier)
      if (ier.ne.0) go to 999   

C  If elements are not recent, propagate orbit and insert placeholders in 
C    element file
      iddif = jd(igyr,1,igday) - jd(iyinit,1,idinit)
      tdif = iddif*864.d2 + secst - secinit 
      if (tdif.gt.tdifmax) then

C   Check for stale elements.dat tfile
         if (tdif.gt.(ndmax*tdifmax)) then
            ier = 1
            orbit%nvec = 0
            print *,' '
            print *,'**************************************************'
            print *,' '
            print *,' ELEMENTS MORE THAN',ndmax,' DAYS PRIOR TO DATA'
            print *,' '
            print *,'         FOR DATA PRIOR TO DECEMBER 24, 2004'
            print *,'     PLEASE DOWNLOAD THE LATEST ELEMENTS.DAT FILE'
            print *,'         FROM THE NASA/GSFC SEAWIFS PROJECT'
            print *,' '
            print *,'             FOR ALL SUBSEQUENT DATA'
            print *,'     PLEASE CONTACT ORBIMAGE (www.orbimage.com)'
            print *,' '
            print *,'**************************************************'
            print *,' '
            go to 999

         else
            
            call add_elements(orbinit,cdrg,iyinit,idinit,secinit,
     *           igyr,igday,secst,tdifmax,irec)

C  Reduce state weights to account for uncertainty in propagation
            wdifac = (tdif/tdifmax)**2
            do i=1,6
               s0(i) = s0(i)/wdifac
            end do

         end if
      end if

      write = .not.fit

C  If less than minimum number of data points or we are using fitted elements
C   do not fit elements to GPS.
      if (ngps.lt.nmin) then
        iter = .false.
        write = .false.
      else
        iter = (.not.fit)
      end if

      do i=1,6
        orbupd(i) = orbinit(i)
      end do

C  Compute number of integrated vectors required at 1 minute intervals
      iddif = jd(igyr,1,igday) - jd(iyinit,1,idinit)
      nstp = iddif*1440.d0 + (secend - secinit)/60.d0 + 90

C  Call ASAP to integrate orbit
      call asaps(li,mi,iplt,orbupd,iyinit,idinit,secinit,nstp,cdrg,
     *    tsap,asap)

      niter = 0

      do while (iter.and.(niter.lt.maxit))

C  Now start fitting algorithm

C  Rotate ASAP vectors to ECEF and interpolate to GPS times
         call asap_rot_int(nstp,iyinit,idinit,tsap,asap,ngps,
     *        igyr,igday,gpsec,vecs)
      
C     Compute partial derivatives for orbit vectors with respect to elements
         call pderiv(ngps,igyr,igday,gpsec,vecs,orbupd,iyinit,idinit,
     *        secinit,driv)

C  Compute updates to orbital elements 
         call fitgps(ngps,gpsvec,nsig,vecs,driv,s0,updorb)

C  Update orbital elements 
         iter = .false.

c  Do some checks on updates

c  If semimajor axis or mean anomaly change by large amounts, hold other
c   elements constant
         if ((abs(updorb(1)).gt.maxupd(1)).or.
     *        (abs(updorb(6)).gt.maxupd(6))) then
            do i=2,5 
               updorb(i) = 0.d0
            end do

         else

c  Else check for update in other elements greater than maximum iteration
            do i=2,5
               if (abs(updorb(i)).gt.maxupd(i))
     *              updorb(i) = sign(maxupd(i),updorb(i))
            end do
         end if
         
c  If eccentricity would be negative, limit correction
         if ((updorb(2)+orbupd(2)).lt.0.d0) updorb(2) = -orbupd(2)

c  Apply updates
         do i=1,6
            orbupd(i) = orbupd(i) + updorb(i)

C  Check if elements have converged within tolerance; if not,
C    another iteration is required
            if (abs(updorb(i)).gt.updtol(i)) iter = .true.
         end do
         orbupd(6) = orbupd(6) - updorb(5)
         print *,updorb

         call asaps(li,mi,iplt,orbupd,iyinit,idinit,secinit,nstp,cdrg,
     *        tsap,asap)

         niter = niter + 1

      end do

      if (niter.ge.maxit) then
         ier = 1
         write = .false.
         print *, ' ORBCOMP:  MAX ITERATIONS PERFORMED'
      end if

C  Now perform final orbit processing
      call asap_rots(iyinit,idinit,tsap,asap,nstp,vecs)

      nstr = iddif*1440.d0 + (secst - secinit)/60.d0
      orbit%nvec = nstp - nstr + 1
      orbit%iyr = iyinit
      orbit%iday = idinit
      do i=nstr,nstp 
        orbit%torb(i-nstr+1) = tsap(i)
        do j=1,3
          orbit%pos(j,i-nstr+1) = vecs(j,i)
          orbit%vel(j,i-nstr+1) = vecs(j+3,i)
        end do
      end do

C  If elements were computed for this interval, write elements for start
C   and end of interval to file
      if (write) then

c  Open elements file
        filnm = '$ELEMENTS/elements.dat'
        call filenv(filnm,filnm)
        open(lun,file=filnm,status='old',access='direct',err=990,
     *    recl=128)

c  First write any placeholder records at ascending node crossings
        ir = irec
        init = .true.
        do i=2,nstr
          if ((asap(3,i).gt.0.).and.(asap(3,i-1)).lt.0.) then
            ir = ir + 1
            call vec2mean(asap(1,i),ge,aj2,orbupd,orbend,ier)
            call put_elements(lun,orbend,cdrg,iyinit,idinit,tsap(i),
     *        ir,init)
            ilast = i
          end if
        end do

c  Write fitted elements to file
        init = .false.
        if (ir.eq.irec) then
          call put_elements(lun,orbupd,cdrg,iyinit,idinit,secinit,
     *        irec,init)
        else
          call put_elements(lun,orbend,cdrg,iyinit,idinit,tsap(ilast),
     *        ir,init)
        end if

c  Write elements for remainder of interval at ascending node crossings
        do i=nstr+1,nstp
          if ((asap(3,i).gt.0.).and.(asap(3,i-1)).lt.0.) then
            init = .true.
            ir = ir + 1
            call vec2mean(asap(1,i),ge,aj2,orbupd,orbend,ier)
            call put_elements(lun,orbend,cdrg,iyinit,idinit,tsap(i),
     *        ir,init)

c  If not end of interval, write as fitted elements
            if ((nstp-i).gt.90) then
              init = .false.
              call put_elements(lun,orbend,cdrg,iyinit,idinit,tsap(i),
     *        ir,init)
            end if
          end if
        end do
          
        close (lun)
      end if
      return

  990 print *, 'PUT_ELEMENTS:  Error opening elements file'
      close(lun)

 999  return
      end
