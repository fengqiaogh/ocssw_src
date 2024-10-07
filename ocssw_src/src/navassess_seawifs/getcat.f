      subroutine getcat(lucat,eptime,maxcat,smaglm,iqlimt,numcat,
     *  idncat,datcat,ierr)

c  This subroutine reads the island catalog and loads the data into
c  output arrays.  Some arguments are included for compatibility with the
c  FDF pattern matching code but not actually used.  This routine may be 
c  augmented at some later date to filter the islands according to the data
c  limits but for now it passes all islands in the catalog.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  lucat        I*4      I      Unit number (not used)
c  eptime       R*8      I      Epoch time (not used)
c  maxcat       I*4      I      Maximum number of catalog entries
c  smaglm       R*4      I      Magnitude limit (not used)
c  iqlimt       I*4      I      Catalog quality limits (1=gaclac flag,
c                                 2=west flag)
c  numcat       I*4      O      Number of catalog entries
c  idncat(*)    I*4      O      ID numbers of catalog entries
c  datcat(7,*)  R*4      O      Catalog entries (first three are unit vector 
c                                 components, fourth is size, last three are 
c                                 not used)
c  ierr         I*4      O      Error code

c       Subprograms Called:

c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               June 13, 1994
c
c       Modification History:   
     

      real*8 eptime
      real*4 datcat(7,16000)
      integer*4 idncat(16000),iqlimt(6)
      character*80 filnm

      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae

      ierr = 0
      numcat = 0

c  Open catalog file
      if (iqlimt(1).eq.0) then
        filnm = '$CATALOG/islands.lac'
      else
        filnm = '$CATALOG/islands.gac'
      end if
      call filenv(filnm,filnm)
      open(file=filnm,unit=52,status='old',iostat=istat)

c  Read and convert data
      dowhile (.true.)
        read(52,*,end=600) xlon,xlat,wlon,wlat
        numcat = numcat + 1
        idncat(numcat) = numcat
        datcat(1,numcat) = cos(xlon/radeg)*cos(xlat/radeg)
        datcat(2,numcat) = sin(xlon/radeg)*cos(xlat/radeg)
        datcat(3,numcat) = sin(xlat/radeg)
        wlon = wlon*cos(xlat/radeg)
        datcat(4,numcat) = amax1(wlat,wlon)
      end do    

  600 close(52)
      return

  999 write(*,*) 'Error opening catalog file'
      return
      end

