c ------------------------------------------------------------------
c function fget_ancillary
c
c F77 wrapper for the SeaWiFS ancillary hdf file reader, 
c get_ancillary. Calling sequence is identical.  Relies
c on the companion C-wrapper get_ancillary_, which was 
c added to getanc.c.
c
c BA Franz, GSC/NASA SeaWiFS, 1/97
c ------------------------------------------------------------------

      integer*4 function fget_ancillary(
     .          lat,lon,nsamp,syear,sday,eday,msec,
     .          file1,file2,file3,parmID,
     .          outval, qcflag)
c
      implicit none
c
#define FLEN    255
c
      integer*2   nsamp
      real*4      lat(nsamp)
      real*4      lon(nsamp)
      integer*2   syear
      integer*2   sday
      integer*2   eday
      integer*4   msec
      character*FLEN file1
      character*FLEN file2
      character*FLEN file3
      integer*2   parmID
      real*4      outval(nsamp)
      integer*2   qcflag(nsamp)
      integer*4   status
c
      integer*2   i
      integer*4   n1, n2, n3
      integer*4   lenstr
      byte        bfile1(FLEN+1)
      byte        bfile2(FLEN+1)
      byte        bfile3(FLEN+1)
c
      integer*4   get_ancillary
      external    get_ancillary 
c
c     !
c     ! Make C-style string from f77 string
c     !
      n1 = lenstr(file1)
      do i=1,n1
          bfile1(i) = ICHAR(file1(i:i))
      enddo
      bfile1(n1+1) = 0
c
      n2 = lenstr(file2)
      do i=1,n2
          bfile2(i) = ICHAR(file2(i:i))
      enddo
      bfile2(n2+1) = 0
c
      n3 = lenstr(file3)
      do i=1,n3
          bfile3(i) = ICHAR(file3(i:i))
      enddo
      bfile3(n3+1) = 0
c
      status = get_ancillary(lat,lon,%val(nsamp),
     .        %val(syear),%val(sday),%val(eday),%val(msec),
     .        bfile1,bfile2,bfile3,
     .        %val(parmID),outval,qcflag)
c
      fget_ancillary = status
      return
c
      end
