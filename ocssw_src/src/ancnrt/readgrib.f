C***************************************************************** 
C* FILE:                        readgrib
C*
C* PURPOSE:             read NMC produced data GRIBs.
C*
C* DESCRIPTION: 
C*      read 1 data record (360,181) and flip it.
C*
C* INPUT PARAMETERS:            
C*      char *(*) filename      -       file to read
C*
C* OUTPUT PARAMETERS:           
C*      real datrtn                             - array of return values
C*      integer year, month, day, hour - temporal info of read data
C*
C* COMMON AREAS:        none
C*
C* LOCAL VARIABLES: 
C*              pds     - header string
C*              data    - full globe of data
C*
C* SUBROUTINES CALLED:
C*              none
C*
C* HISTORY:
C*
C* AUTHOR:      Brian D. Schieber, GSC, 6/92
C*
C* MODIFICATION HISTORY:
C*  BDS, 9/18/96        Support NEW 1x1 deg GDAS1 data files
C*
C***************************************************************** 

      subroutine readgrib(filename, datrtn, year, month, day, hour)
c
c     Mods:
c
      character*(*)    filename
      real             datarr(65160)
      real             datrtn(65160)
      real             global(360,181)
      integer          i, j, jj, k
      integer          ipdsl, lenkgds, nwords
      integer          year, month, day, hour
                integer                   kgds(200)

      character * 1    pds(50)
c
c     equivalence nhem and shem so equator overlaps
c
c      equivalence      (shem(1),datarr(1))
c      equivalence      (nhem(1),datarr(5221))

c     equivalence datarr to global array to get lat/lon indices for image flip
c      equivalence      (datarr(1), global(360,181))
c
c     open file
c   
      open(9, err=998, file=filename, type='unknown', 
     +     form='unformatted')
c
c     read grid
c
      read (9) ipdsl, lenkgds, nwords
      read (9) (pds(jj),jj=1,ipdsl)
      IF (lenkgds.ne.0) read (9) (kgds(jj),jj=1,lenkgds)
      read (9) (datarr(j),j=1,nwords)

                year  = ichar (pds(13))
c               print *,'year ', year
                month = ichar (pds(14))
c               print *,'month ', month
                day   = ichar (pds(15))
c               print *,'day ', day
                hour  = ichar (pds(16))
c               print *,'hour ', hour
c
c               k = 0
c      do j = 181, 1, -1
c       do i = 1, 360
c                               k = k+1
c                       datrtn(k) = global(j,i)
c       end do
c      end do

        do i = 1, 65160
                datrtn(i) = datarr(i)
      end do
c
      goto 1000

 998  print *, 'error opening file'
      status = -1    

1000  continue
      return 
      end

