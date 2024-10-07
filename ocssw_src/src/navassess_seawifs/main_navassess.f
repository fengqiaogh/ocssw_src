c
      program nav_assess

      real*8 dangtl,dmagtl,pangtl,tangtl,tminco
      real*8 dmagtg,pangtg,tangtg,tmincg
      integer*4 prod_ID(2), irec, ind, gaclac
      integer sffattr, sfrattr
      real*4 xlat(1000),xlon(1000)
      real*4 wlat(1000),wlon(1000)
      integer*4 nump(1000),nill(1000),nilp(1000)
      character*80 infile,parmfile
      character*32 statfile
      character*1 bytefile(32), isl(3)
      logical west

      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae

      integer*4 levdbg(8),ludbug
      common /cmdebg/levdbg,ludbug
      common /idparm/dangtl,dmagtl,pangtl,tangtl,tminco,
     *  dmagtg,pangtg,tangtg,tmincg,imatch
      namelist /idnl/imatch,dangtl,dmagtl,pangtl,tangtl,tminco,
     *  dmagtg,pangtg,tangtg,tmincg,levdbg,ludbug
      data imatch/3/,dangtl/0.2/
      data dmagtl/0.025/,pangtl/0.01/,tangtl/0.01/,tminco/30.0/
      data dmagtg/0.1/,pangtg/0.04/,tangtg/0.03/,tmincg/15.0/
      data levdbg/8*0/,ludbug/99/,nllun/50/
      data isl/'I','S','L'/
      equivalence (statfile,bytefile)


      parmfile = '$IDPARMS/idparms.nl'
      call filenv(parmfile,parmfile)
      open(file=parmfile,unit=nllun,status='old',iostat=istat)
      read(nllun,idnl)
      write(*,idnl)
      close(nllun)

c      write( 6, 100 ) 
c  100 format( ' Enter the name of the HDF file to open' )
c      read( 5, 200 ) infile
c  200 format( a )
        call getarg(1,infile)
      write( 6, 300 ) infile
  300 format( ' nav_assess.f: input file =',/,a,/ )

c       call the open routine
      call get_l1a_open( infile, prod_ID, npix, nlin, iret )
c
      write( 6, 400 ) prod_ID(1), npix, nlin, iret
  400 format( ' test, open: prod_ID = ', i7, '   npix = ',i7,/,
     1     ' nlin = ',i7, ' iret = ',i7 )
c
      if (iret.eq.-1) go to 999

      gaclac = 0
      if (npix.eq.248) gaclac = 1

c  Get file name to use for output stats file

      ind = sffattr(prod_ID, 'Product Name')
      iret = sfrattr(prod_ID, ind, statfile)
      do i=1,3
         bytefile(i+15) = isl(i)
      end do

      call cdata

      call find_islands(prod_ID,npix,west,numi,xlon,xlat,wlon,wlat,
     *  nump,nill,nilp)

      call get_l1a_close( prod_ID, iret )
      write( 6,600) iret
 600  format(' test, close: iret = ', i7)

      if (numi.gt.0) then
        open(55, file=statfile)

        call id_drv(gaclac,west,numi,xlon,xlat,wlon,wlat,nilp,nump)
      end if    

      stop

 999  write(*,*) 'Invalid file name'

      stop

      end
