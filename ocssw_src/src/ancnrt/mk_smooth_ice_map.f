!  this routine averages the ncep daily ice fields onto a 0.5 deg lat/lon grid
!  the averaging window is a circle with a radius of 50 km.

!  input:  filename1 is the file name for the ncep ice map.
!  the file contains a character (1-byte) array that is 4320 by 2160

      subroutine mk_smooth_ice_map(char_ice, frac_ice_smoothed)
      implicit none

      real(8), parameter :: re=6378.137d3 !earth equatorial radius (meters)
      real(8), parameter :: rp=6356.752d3 !earth      polar radius (meters)
      real(8), parameter :: ffac=(rp/re)*(rp/re)

      integer(4), parameter :: nlat=2160, nlon=4320

      character(1), intent(in) :: char_ice(nlon,nlat)
      real(4),     intent(out) :: frac_ice_smoothed(720,360)

      integer(4) istart
      integer(4) ilat,ilon,jlat,jlon,klat,klon,kklon

      real(4) xlat,xlatg,xlon
      real(4) sinlat0,coslat0
      real(4) cel0(3),dif(3),rsq,rsq_max,wt

      real(4) coslat(nlat),cel(3,nlon,nlat)
      integer(4) nstep(nlat)

      real(4) frac_ice(nlon,nlat), dummy_lat(nlon)
      real(8) tsum(0:1)

      real(4) cosd, sind, tand
      real(8) datand

      data istart/1/

      if(istart.eq.1) then
         istart=0

         rsq_max=(2*sind(0.5*0.45d0))**2 !0.45 degree = 0.45 * 111.2 km/deg = 50 km radius

         do ilat=1,nlat
            xlat=(ilat-1)/12.-89.95833
            xlatg=datand(ffac*tand(xlat)) !convert from geodetic to geocentric
            sinlat0=sind(xlatg)
            coslat0=cosd(xlatg)

            coslat(ilat)=coslat0

            nstep(ilat)=nint(6./coslat0) + 1
            if(nstep(ilat).ge.nlon/2) nstep(ilat)=nlon/2 -1
            if(ilat.le.12 .or.ilat.ge.nlat-11)  nstep(ilat)=nlon/2 -1

            do ilon=1,nlon
               xlon=(ilon-1)/12.+ 0.04167
               cel(1,ilon,ilat)=cosd(xlon)*coslat0
               cel(2,ilon,ilat)=sind(xlon)*coslat0
               cel(3,ilon,ilat)=sinlat0
            enddo               !ilon
         enddo                  !ilat
      endif

      frac_ice=ichar(char_ice)*0.01

!      write(*,*) 'Smoothing ice file'
      frac_ice_smoothed = 0
      do ilat=1,nlat

         xlat=(ilat-1)/12.-89.95833
         xlatg=datand(ffac*tand(xlat)) !convert from geodetic to geocentric
         sinlat0=sind(xlatg)
         coslat0=cosd(xlatg)

!         write(*,*) ilat

         do ilon=1,nlon

            xlon=(ilon-1)/12.+ 0.04167

            cel0(1)=cosd(xlon)*coslat0
            cel0(2)=sind(xlon)*coslat0
            cel0(3)=sinlat0

            tsum=0

            do klat=ilat-6,ilat+6
               if(klat.lt.1 .or. klat.gt.nlat) cycle
               wt=coslat(klat)

               do klon=ilon-nstep(klat),ilon+nstep(klat)+1
                  kklon=klon
                  if(kklon.lt.   1) kklon=kklon+nlon
                  if(kklon.gt.nlon) kklon=kklon-nlon

                  dif=cel(:,kklon,klat)-cel0

                  rsq=dif(1)*dif(1) + dif(2)*dif(2) + dif(3)*dif(3)

                  if(rsq.gt.rsq_max) cycle

                  tsum(0)=tsum(0) + wt
                  tsum(1)=tsum(1) + wt*frac_ice(kklon,klat)

               enddo            !klon
            enddo               !klat

            jlat=1 + int(2*(xlat+90))
            jlon=1 + int(2*xlon)
            if(jlat.lt.1 .or. jlat.gt.360) stop 'error1'
            if(jlon.lt.1 .or. jlon.gt.720) stop 'error2'

            frac_ice_smoothed(jlon,jlat)=tsum(1)/tsum(0)

         enddo                  !ilon
      enddo                     !ilat

      if(minval(frac_ice_smoothed).lt.0 .or. minval(frac_ice_smoothed).gt.1)  stop 'frac_ice_smooth oob, pgm stopped'
      
      return
      end
