      subroutine rdfiles(infile,gdasc,nfiles,status)
c*****************************************************************
c *** program to read input file names from an ASCII file
c
c     Jim Firestone, GSC/SAIC, Feb. 1993
c
c     Inputs:
c
c     INFILE: Char*80 - name of file containing ASCII files to read
c       This file contains as its first line the number of files to 
c       follow, with the value written in columns 2-5. The second 
c                 and subsequent lines contain the file names, one per line.
c
c     Outputs:
c
c     GDASC: Char(10000)*255 - names of up to 10000 files 
c                       containing gridded ASCII data.
c
c     NFILES: Integer - the number of files returned in the GDASC array.
c
c     STATUS: Integer - status return code for the routine, 0 for normal
c       termination or 1 if INFILE could not be read.
c*****************************************************************
c
      implicit none
      character infile*80, gdasc(10000)*255
      integer nfiles, status, i, lnstrg, len
c
      status = 0
      open(9,err=998,file=infile,recl=80)
      read(9,10) nfiles
10    format(i3)
c10    format(1x,i4)
      do i = 1, nfiles
        read(9,20) gdasc(i)
20      format(a255)
        len = lnstrg(gdasc(i))+1
        gdasc(i)(len:len) = achar(0)
      end do
      close(9)
      goto 999
998   write(*,*) ' error opening file ', infile, ' terminating...'
      status = 1
999   continue
      return
      end
