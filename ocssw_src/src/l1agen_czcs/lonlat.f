      subroutine lonlat(alon,alat,xlon,ylat)
c
c  Interpolates lineraly across anchor points to fill up lons and
c  lats.
c
c  Modified to offset scan by a fixed number of pixels
c  F.S. Patt, SAIC GSC, December 21, 2000.
      parameter(npix=1968,nanc=77)
      integer ianc(nanc)
      real alon(nanc),alat(nanc)
      real xlon(npix),ylat(npix)
      data ianc /1,16,31,46,61,76,91,106,121,136,151,166,181,196,216,
     +236,256,276,296,316,341,366,391,416,441,466,496,526,556,591,626,
     +666,706,751,796,841,886,931,984,1037,1082,1127,1172,1217,1262,
     +1302,1342,1377,1412,1442,1472,1502,1527,1552,1577,1602,1627,1652,
     +1672,1692,1712,1732,1752,1772,1787,1802,1817,1832,1847,1862,1877,
     +1892,1907,1922,1937,1952,1968/
      data poff/5.5/
c
c  Interpolate across scan lines
      do ia = 1,nanc-1
       if (alon(ia) .ge. 90.0 .and. alon(ia+1) .le. -90.0)then
         rm = ((alon(ia+1)-alon(ia)+360.0))/
     $        (float(ianc(ia+1))-float(ianc(ia)))
        else if( alon(ia) .lt. -90.0 .and. alon(ia+1) .gt. 90.0)then
         rm = ((alon(ia+1)-alon(ia)-360.0))/
c     $        (float(ianc(ia))-float(ianc(ia+1)))
     $        (float(ianc(ia+1))-float(ianc(ia)))
        else
         rm = (alon(ia+1)-alon(ia))/(float(ianc(ia+1))-float(ianc(ia)))
       end if
c
       b = alon(ia) - rm*float(ianc(ia))
       rm2 = (alat(ia+1)-alat(ia))/(float(ianc(ia+1))-float(ianc(ia)))
       b2 = alat(ia) - rm2*float(ianc(ia))
       i1 = ianc(ia) + poff + 1
       i2 = ianc(ia+1) + poff
       if (ia.eq.1) i1 = 1
       if (ia.eq.(nanc-1)) i2 = npix
       do i = i1,i2
        xlon(i) = rm*(float(i)-poff)+b
        if (xlon(i) .gt. 180.0)then
         xlon(i)=((xlon(i) - 180.0) - 180.0)
        end if
        ylat(i) = rm2*(float(i)-poff)+b2
       enddo
      enddo
c      xlon(npix) = alon(nanc)
c      ylat(npix) = alat(nanc)
c
      return
      end
