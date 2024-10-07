      subroutine get_elements(iyr,iday,sec,secend,fit,orbout,cdrg,
     *  iyorb,idorb,secorb,irec,ier)

c $Header$
c $Log$                                                                        
c
c  Purpose:  This routine reads orbital elements for use by the ASAP orbit
c               integrator from a file, according to the specified time.  
c               The first element set is read which precedes the specified 
c               time.  The file contains both fitted and unfitted (to GPS 
c               data) element sets.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  iyr          I*4      I      Year of date for requested elements
c  iday         I*4      I      Day of date for requested elements
c  sec          R*8      I      Seconds of day for requested elements
c  secend       R*8      I      Seconds of day at end of interval
c  fit          L*4     I/O     Fitting flag for elements
c                                =.true., use fitted elements if available
c                                =.false., use unfitted elements only
c                                output value reflects elements used.
c  orbout(6)    R*8      O      Orbital elements
c  cdrg         R*8      O      Drag coefficient for ASAP model
c  iyorb        I*4      O      Year of element epoch
c  idorb        I*4      O      Day of element epoch
c  secorb       R*8      O      Seconds of day of element epoch
c  irec         I*4      O      Record number of elements in file
c  ier          I*4      O      Error code:     =0, success
c                                               =1, error
c
c  By: Frederick S. Patt, GSC, December 23, 1993
c
c  Notes:  
c
c  Modification History:
c
c  Added check for non-fitted elements within time span of data, which will
c  cause fit flag to be set to false.  F.S. Patt, GSC, October 26, 1994.
c
c  Fixed a typo in the comparison of orbital element epochs to input times.
c  F.S. Patt, GSC, August 15, 1996.
c
c  Added conditional compile for Sun OS record length differences. B. A.
c  Franz, GSC, November 14, 1997.
c
c  Added record number in direct access read for Sun OS compatibility.
c  B. A. Franz, November 14, 1997.
c
c  Added check for last record in file more than one orbit before end
c  of data.  F.S. Patt, SAIC GSC, February 12, 1998
c
c  Fixed a bug introduced in the previous change, which caused the fit
c  flag in the elements.dat record not to be checked if only the last record
c  was read.  F.S. Patt, SAIC GSC, March 10, 1998.

      implicit none
#include "nav_cnst.fin"

      real*8 orbout(6),orb1(6),orb2(6),cdrg,sec,secend
      real*8 secorb,tdif,s,cd,tdifen,oneorb
      integer*4 iyr,iday,iyorb,idorb,irec,ier,lun,spare
      integer*4 jd,jdr,jdo,i,nrecs,iy,id
      character*256 filnm
      logical fit,forb
      data oneorb/5940.d0/
      data lun/17/

c  Open elements file
      filnm = '$ELEMENTS/elements.dat'
      call filenv(filnm,filnm)
      open(lun,file=filnm,status='old',access='direct',err=990,
     *  recl=128, convert='big_endian')

c  Read header record
      read(lun,rec=1,err=990) nrecs
        
      jdr = jd(iyr,1,iday)

c  Read last record and check if more than one orbit before end of data
      irec = nrecs 
      write(*,*) irec
      read (lun,rec=irec,err=990) iy,id,s,orb1,orb2,cd,forb,spare
      write(*,*) iy,id,s,orb1,orb2,cd,forb
      jdo = jd(iy,1,id)
      tdif = (jdo-jdr)*864.d2 + s - sec
      tdifen = tdif + sec - secend
      if (tdifen.lt.(-oneorb)) fit = .false. 

c  Read backwards through file to find most recent elements set which
c   precedes start of data
      dowhile (tdif.gt.0.)

        irec = irec - 1
        if (irec.le.1) then
          print *,'GET_ELEMENTS:',
     *'  No elements in file preceding requested time'
          go to 990
        end if
        read (lun,rec=irec,err=990) iy,id,s,orb1,orb2,cd,forb,spare
        jdo = jd(iy,1,id)
        tdif = (jdo-jdr)*864.d2 + s - sec

c  Check if unfitted elements in interval
        tdifen = tdif + sec - secend
        if (tdifen.lt.0.) fit = (fit.and.forb)
      end do

c  Set fit flag using input and record flags
      fit = (fit.and.forb)
      
c  Select element set from record according to logical flags
      if (fit) then
        do i=1,6
          orbout(i) = orb2(i)
        end do
      else
        do i=1,6
          orbout(i) = orb1(i)
        end do
      end if
      iyorb = iy
      idorb = id
      secorb = s
      cdrg = cd
      print *,'GET_ELEMENTS:',iyorb,idorb,secorb,fit
      print *,orbout
      ier = 0
      close(lun)

      return

  990 print *, 'Error reading elements file'
      ier = 1
      close(lun)

      return
      end
