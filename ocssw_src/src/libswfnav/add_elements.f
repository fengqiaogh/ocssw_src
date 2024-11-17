      subroutine add_elements(orbinit,cdrg,iyinit,idinit,secinit,
     *  igyr,igday,secst,tdifmax,irec)

c $Header: /app/shared/RCS/irix-5.2/seawifsd/src/mops/mopsV4.1/swfnav/add_elements.f,v 1.1 1995 /01/17 23:01:58 seawifsd Exp seawifsd $
c $Log: add_elements.f,v $
c Revision 1.1  1995/01/17 23:01:58  seawifsd
c Initial revision
c

c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  orbinit(6)   R*8     I/O     Orbital elements
c  cdrg         R*8      O      Drag coefficient for ASAP model
c  iyinit       I*4     I/O     Year of element epoch
c  idinit       I*4     I/O     Day of element epoch
c  secinit      R*8     I/O     Seconds of day of element epoch
c  igyr         I*4      I      Year of data start time
c  igday        I*4      I      Day of data start time
c  secst        R*8      I      Seconds of data start time
c  tdifmax      R*8      I      Maximum allowed time between element epoch
c                                and data start time
c  irec         I*4     I/O     Record number for elements in file
c
c  By: Frederick S. Patt, GSC, October 20, 1994
c
c  Notes:  
c
c  Modification History:
c
c  Added conditional compile for Sun OS record length differences. B. A.
c  Franz, GSC, November 14, 1997.
c

      real*8 orbinit(6),secinit,secst,cdrg,ge,aj2,tdifmax
      real*8 asap(6,1500),tsap(1500),tdif,orbupd(6)
      integer*4 igyr,igday,jd,i,j,ilast,lun,jdl
      integer*4 iyinit,idinit,irec,iddif,nstp,li,mi,iplt,ier
      character*80 filnm
      logical*4 init
      data ge/3.9860050D5/,aj2/0.10826270E-02/
      data li/30/,mi/15/,iplt/0/,lun/17/

      iddif = jd(igyr,1,igday) - jd(iyinit,1,idinit)
      tdif = iddif*864.d2 + secst - secinit 
      print *,'ADD_ELEMENTS: time difference',tdif

c  Open elements file
      filnm = '$ELEMENTS/elements.dat'
      call filenv(filnm,filnm)
      open(lun,file=filnm,status='old',access='direct',err=990,
     *  recl=128)

c  Continue while time difference between data start time and element epoch 
c   exceeds maximum allowed
      dowhile (tdif.gt.tdifmax)
        
c  Integrate orbit using current elements

        nstp = tdif/60.d0 + 1 
        if (nstp.gt.1500) nstp = 1500
        call asaps(li,mi,iplt,orbinit,iyinit,idinit,secinit,nstp,cdrg,
     *    tsap,asap)

c  Write any placeholder records at ascending node crossings
        init = .true.
        i = 2
        dowhile (i.le.nstp)
          if ((asap(3,i).gt.0.).and.(asap(3,i-1)).lt.0.) then
            irec = irec + 1
            ier = 1
            i = i - 1
            dowhile ((ier.ne.0.).and.(i.le.nstp))
              i = i + 1
              call vec2mean(asap(1,i),ge,aj2,orbinit,orbupd,ier)
            end do
            call put_elements(lun,orbupd,cdrg,iyinit,idinit,tsap(i),
     *        irec,init)
            ilast = i
          end if
          i = i + 1
        end do

c  Update current elements and epoch (keep eccentricity from initial elements
        orbinit(1) = orbupd(1)
        do j=3,6
          orbinit(j) = orbupd(j)
        end do
        idinit = idinit + tsap(ilast)/864.d2

c  Check for year end rollover 
        if (idinit.ge.365) then
          jdl = jd(iyinit,1,idinit)
          call jdate(jdl,iyinit,idinit)
        end if
        secinit = mod(tsap(ilast),864.d2)
        iddif = jd(igyr,1,igday) - jd(iyinit,1,idinit)
        tdif = iddif*864.d2 + secst - secinit 

      end do
      close(lun)
      return

  990 print *, 'PUT_ELEMENTS:  Error opening elements file'
      close(lun)
      return

      end
