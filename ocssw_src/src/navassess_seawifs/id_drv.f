      subroutine id_drv(gaclac,west,numi,xlon,xlat,wlon,wlat,nilp,nump)

c  This subroutine is a driver for the FDF pattern matching routine IDENTY.
c  It converts the input island latitudes and longitudes to unit vectors and
c  sets up some additional variables and parameters for the identification 
c  algorithm.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  gaclac       I*4      I      GAC/LAC flag (1 = GAC)
c  west         L*4      I      Flag to indicate scene starts in Western
c                               hemisphere (sets longitude range at 0 to 360
c                               instead of -180 to 180 degrees)
c  numi         I*4      I      Number of islands found
c  xlon(1000)   R*4      I      Longitude centroid of each island
c  xlat(1000)   R*4      I      Latitude centroid of each island
c  wlon(1000)   R*4      I      Longitude extent of each island
c  wlat(1000)   R*4      I      Latitude extent of each island
c  nilp(1000)   I*4      I      Pixel number in scan line for each island
c  nump(1000)   I*4      I      Number of pixels for each island

c       Subprograms Called:     identy

c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               June 13, 1994
c
c       Modification History:   
     

      real*4 xlon(1000),xlat(1000),wlon(1000),wlat(1000)
      integer*4 nilp(1000), nump(1000), gaclac
      real*8 eptime,dangtl,dmagtl,pangtl,tangtl,tminco
      real*8 timclm(1000),dmagtg,pangtg,tangtg,tmincg
      real*4 gcrclm(3,1000),briclm(1000),vrmclm(1000),vrpclm(1000)
      real*4 datcat(7,16000),skyclm(10,3,1000)
      real*4 dv(3),dlon(1000),dlat(1000),dlom,dloq,dlam,dlaq
      integer*4 nobclm(1000),mrkclm(1000),klmstr(1000)
      integer*4 mrkstr(1000),lblclm(1000),idfclm(1000),nrfclm(1000)
      integer*4 idncat(16000),mapclm(10,1000),idfhst(1000),iqlimt(6)
      logical west

      real*8 pi,radeg,re,rem,f,omf2,omegae
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
      common /idparm/dangtl,dmagtl,pangtl,tangtl,tminco,
     *  dmagtg,pangtg,tangtg,tmincg,imatch

      data eptime/0.d0/,timclm/1000*0.d0/,ifxclm/10/
      data vrmclm/1000*0.d0/,vrpclm/1000*0.d0/,iqlimt/6*0/
      data nobclm/1000*1/,mrkclm/1000*0/,mrkstr/1000*0/
      data sinc/0.00159/
      data dmismn/0.02/,dnotmx/0.015/,dmismg/0.07/,dnotmg/0.05/
      

      open (54,file='nav_assess.out')
c      open (55,file='nav_stats.out')

c  Define position matching parameters

      sind = sinc
      inpix = 643
c  If GAC then
      if (gaclac.eq.1) then
        dmagtl = dmagtg
        pangtl = pangtg
        tangtl = tangtg
        tminco = tmincg
        sind = sinc*4.0
        inpix = 125
        dmismn = dmismg
        dnotmx = dnotmg
      end if

      rporm = 7083./rem
      rpor2 = rporm*rporm
   
c  Convert input positions to unit vectors and
c  load island parameters into arrays
      do i=1,numi
        gcrclm(1,i) = cos(xlon(i)/radeg)*cos(xlat(i)/radeg)
        gcrclm(2,i) = sin(xlon(i)/radeg)*cos(xlat(i)/radeg)
        gcrclm(3,i) = sin(xlat(i)/radeg)

c  Store island "magnitude" as largest dimension
        zlon = wlon(i)*cos(xlat(i)/radeg)
        briclm(i) = amax1(zlon,wlat(i)) + 0.01

c  Estimate position uncertainty from pixel number and number of pixels
        th = sind*(nilp(i)-inpix)
        vrpclm(i) = rporm*cos(th)/sqrt(1.0-rpor2*(sin(th))**2) - 1
        vrpclm(i) = vrpclm(i)/sqrt(float(nump(i)))
        klmstr(i) = i
        lblclm(i) = i
      end do
      numstr = numi
      numclm = numi

c  Store gaclac and west flags in catalog quality array
      iqlimt(1) = gaclac
      if (west) iqlimt(2) = 1


c  Call identification routine
      call identy (eptime, imatch, dangtl, dmagtl, pangtl, tangtl,
     *           tminco, smaglm, iqlimt, maxcat, idncat, datcat,
     *           ifxclm, numclm, timclm, gcrclm, briclm, vrmclm,
     *           vrpclm, nobclm, mrkclm, numstr, klmstr, mrkstr, 
     *           lblclm, idfclm, nrfclm, mapclm, skyclm, idfhst,
     *           ircode)                                                        

c  Print results
      write(*,*)'IDENTY return code',ircode
      nid = 0
      mis = 0
      nnt = 0
      do i=1,numi
        dp = 0.0

c   If at least one candidate, compute position difference
        if (nrfclm(i).ge.1) then
          do j=1,3
            dv(j) = gcrclm(j,i) - skyclm(1,j,i)
            dp = dp + dv(j)*dv(j)
          end do
          dp = radeg*sqrt(dp)

c   If this is a positive ID
          if (mrkclm(i).eq.0) then 
            nid = nid + 1

c   Compute difference in longitude and latitude and write stats to file
            clat2 = gcrclm(1,i)**2 + gcrclm(2,i)**2
            clat = sqrt(clat2)
            dlat(nid) = radeg*dv(3)/clat
            dlon(nid) = radeg*(dv(2)*gcrclm(1,i)-dv(1)*gcrclm(2,i))/clat2
            l = lblclm(i)
            write (55,1200) lblclm(i),xlon(l),xlat(l),nump(l),nilp(l),
     *          mapclm(1,i),dlon(nid),dlat(nid)
 1200       format (i6,2f12.5,3i8,2f10.5)

c    If position difference exceeds tolerance, this may be a mis-ID
            if (dp.gt.dmismn) mis = mis + 1

c   Else if not ID but position difference is small
          else if (dp.lt.dnotmx) then
            nnt = nnt + 1
          end if
          write(*,*)i,lblclm(i),mrkclm(i),nrfclm(i),mapclm(1,i),dp
        end if
        write (54,*) i,lblclm(i),mrkclm(i),nrfclm(i),mapclm(1,i),dp
      end do
      write(*,*)'IDs',nid,'  Mis',mis,'  Not',nnt

c  If at least 5 IDs, compute summary statistics
      if (nid.ge.5) then
        call mediqr(dlon,nid,dlom,dloq)
        write(*,*)'Longitude median, IQR =',dlom,dloq
        call mediqr(dlat,nid,dlam,dlaq)
        write(*,*)' Latitude median, IQR =',dlam,dlaq
      end if

      return
      end

