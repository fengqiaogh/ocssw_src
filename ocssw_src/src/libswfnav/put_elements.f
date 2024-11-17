      subroutine put_elements(lun,orbupd,cdrg,iyr,iday,sec,irec,init)

c $Header$
c $Log$
c
c  Purpose:  This routine writes orbital elements to a file.  The updated
c               elements for the current record and the initial elements 
c               for the following record are written.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  lun          I*4      I      Unit number for elements file
c  orbupd(6)    R*8      I      Updated orbital elements
c  cdrg         R*8      I      Drag coefficient for ASAP model
c  iyr          I*4      I      Year of element epoch
c  iday         I*4      I      Day of element epoch
c  sec          R*8      I      Seconds of day of element epoch
c  irec         I*4      I      Record number of elements in file
c  init         L*4      I      Initial or updated elements flag
c
c  By: Frederick S. Patt, GSC, December 23, 1993
c
c  Notes:  
c
c  Modification History:
c
c  Modified to write one set of elements per call, with argument specifying
c  whether this is first or second set in record - F.S. Patt, 10/21/94.
c
c  Modified to not overwrite updated (fitted) elements with initial elements.
c  F. S. Patt, 1/28/96
c
      implicit none
#include "nav_cnst.fin"

      real*8 orbupd(6),cdrg,sec
      real*8 orb1(6),orb2(6),secorb,secin
      integer*4 iyr,iday,iyorb,idorb,idin,irec,lun,spare
      integer*4 jd,jdr,jdo,i,nrecs
      logical forb,init

c  Read header record
      read(lun,err=990,rec=1) nrecs
        
c  If initial elements for record, write integrated elements to file.
      if (init) then

c  First check for existing record with fitted elements
        forb = .false.
        if (irec.le.nrecs) read (lun,rec=irec,err=990) iyorb,idorb,
     *      secorb,orb1,orb2,cdrg,forb,spare

c  If record doesn't exist or not previously fitted, write elements
        if (.not.forb) then
          do i=1,6
            orb1(i) = orbupd(i)
            orb2(i) = 0.d0
          end do
          idorb = iday + sec/864.d2
          iyorb = iyr       
        
c  Check for year end rollover
          if (idorb.ge.365) then
            jdr = jd(iyorb,1,idorb)
            call jdate(jdr,iyorb,idorb)
          end if
          secorb = mod(sec,864.d2)
          print *,'PUT_ELEMENTS:',irec,iyorb,idorb,secorb
          print *,orb1,cdrg,forb
          write (lun,rec=irec,err=991) iyorb,idorb,secorb,orb1,orb2,
     *      cdrg,forb,spare

c  Check to see if number of records in file needs to be updated
          if (irec.gt.nrecs) then
            nrecs = irec
            write(lun,rec=1,err=991) nrecs
          end if
        end if

c  Else these are updated elements
      else

c  Read record irec from file
        read (lun,rec=irec,err=990) iyorb,idorb,secorb,orb1,orb2,
     *    cdrg,forb,spare
        idin = iday + sec/864.d2
        secin = mod(sec,864.d2)
        jdo = jd(iyorb,1,idorb)
        jdr = jd(iyr,1,idin)
        if ((jdo.ne.jdr).or.(secorb.ne.secin)) then
          print *,'PUT_ELEMENTS:  Input time does not match record',
     *     irec,' time'
        else

c  Load updated elements and write to file
          do i=1,6
            orb2(i) = orbupd(i)
          end do
          forb = .true.
          print *,'PUT_ELEMENTS:',irec,iyorb,idorb,secorb
          print *,orb2,cdrg,forb
          write (lun,rec=irec,err=991) iyorb,idorb,secorb,orb1,orb2,
     *      cdrg,forb,spare
        end if
      end if

      return

  990 print *, 'PUT_ELEMENTS:  Error reading elements file'
      close(lun)

      return
  991 print *, 'PUT_ELEMENTS:  Error writing elements file'
      close(lun)

      return
      end
