      subroutine put_l1a_metadata( prod_ID, nlin, pos, smat, coef, iret)

c  put_l1a_metadata( prod_ID, nlin, pos, smat, coef, iret)
c
c  Purpose:  computes the granule-level and scan-line geographic metadata 
c     fields and writes them to a SeaWiFS L1A file 
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     prod_ID      I*4      I      HDF file ID for open file 
c     nlin         I*4      I      number of scan lines in the L1A file
c     pos(3,*)     R*4      I      size 3 x nlin array of orbit vectors
c     smat(3,3,*)  R*4      I      size 3 x 3 x nlin array of sensor matrices
c     coef(6,*)    R*4      I      size 3 x nlin array of ellipse coefficients
c     iret         I*4      O      return code
c                                    =0, success
c                                    =-1, failure    
c
c  By: F. S. Patt, SAIC GSC, 5 Feb 99
c
c  Notes:  
c
c  Modification History:
c

#include "nav_cnst.fin"

      real*4 smat(3,3,*), pos(3,*), coef(6,*)
      real*4 xlon(maxlin,3), xlat(maxlin,3), xla(3) 
      real*4 xlo(3), solz(3), sola(3), senz(3), sena(3)
      real*4 csolz(maxlin),  sunr(3,maxlin)
      integer*4 prod_ID, iret
      integer*4 sr_id, sl_id, at_id, ind, dfflt
      integer*4 istart(3), istr(3), idims(3)
      integer*4 nflag(maxlin)
      integer sfn2index, sfselect, sfwdata, sfrdata, sfendacc
      integer sfsattr
      logical cross
      data istart/3*0/, istr/3*1/, dfflt/5/

      iret = 0

c  First read Sun reference vectors from file
      ind = sfn2index(prod_ID, 'sun_ref')
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting index for sun_ref'
         return
      end if
      sr_id = sfselect(prod_ID, ind)
      if (sr_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting sun_ref'
         return
      end if
      idims(1) = 3
      idims(2) = nlin
      iret = sfrdata(sr_id, istart, istr, idims, sunr)
      if (iret.eq.-1) then
        write(*,*) 'Error reading sun_ref'
        return
      end if
      iret = sfendacc(sr_id)


c  Next compute scan line end and center coordinates for each line
c   Also find min and max coordinates

      nsta = 1
      ninc = 642
      npix = 3
      cross = .false.

      do ilin = 1, nlin
         call geonav( pos(1,ilin), smat(1,1,ilin), coef(1,ilin),
     1        sunr(1,ilin), nsta, ninc, npix, xla, xlo, 
     2        solz, sola, senz, sena )
         csolz(ilin) = solz(2)
         do j=1,3
            xlon(ilin,j) = xlo(j)
            xlat(ilin,j) = xla(j)
         end do
         

c  Check for date line crossing
         if ( xlo(1) .gt. xlo(3) ) cross = .true.

      end do

c  Find coordinate extrema
c   Not worried about exceptional cases here
      xlanth = max(xlat(1,1),xlat(1,2))
      xlasth = min(xlat(nlin,2),xlat(nlin,3) )
      xlowst = xlon(nlin,1)
      xloest = xlon(1,3)
      icntr = nlin/2
      xlactr = xlat(icntr,2)
      xloctr = xlon(icntr,2)
      solzct = csolz(icntr)

c  Write scan line metadata

c  Scan start latitude      
      call put_scan_metadata (prod_ID, 'slat', nlin, xlat(1,1), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing slat'
        return
      end if
      
c  Scan center latitude      
      call put_scan_metadata (prod_ID, 'clat', nlin, xlat(1,2), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing clat'
        return
      end if
      
c  Scan end latitude      
      call put_scan_metadata (prod_ID, 'elat', nlin, xlat(1,3), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing elat'
        return
      end if
      
c  Scan start longitude      
      call put_scan_metadata (prod_ID, 'slon', nlin, xlon(1,1), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing slon'
        return
      end if
      
c  Scan center longitude      
      call put_scan_metadata (prod_ID, 'clon', nlin, xlon(1,2), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing clon'
        return
      end if
      
c  Scan end longitude      
      call put_scan_metadata (prod_ID, 'elon', nlin, xlon(1,3), 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing elon'
        return
      end if
      
c  Scan center solar zenith      
      call put_scan_metadata (prod_ID, 'csol_z', nlin, csolz, 
     1     iret )
      if (iret.eq.-1) then
        write(*,*) 'Error writing csol_z'
        return
      end if
      
c  End of writing scan line metadat

c  Write scene coordinate metadata

      iret = sfsattr(prod_ID, 'Scene Center Solar Zenith', dfflt, 
     1     1, solzct )
      if (iret.eq.-1) then
        write(*,*)'Error writing Scene Center Solar Zenith'
        return
      endif

      iret = sfsattr(prod_ID, 'Scene Center Latitude', dfflt, 
     1     1, xlactr)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Scene Center Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Scene Center Longitude', dfflt, 
     1     1, xloctr)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Scene Center Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Upper Right Latitude', dfflt, 
     1     1, xlat(1,3))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Upper Right Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Upper Right Longitude', dfflt, 
     1     1, xlon(1,3))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Upper Right Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Upper Left Latitude', dfflt, 
     1     1, xlat(1,1))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Upper Left Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Upper Left Longitude', dfflt, 
     1     1, xlon(1,1))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Upper Left Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Lower Left Latitude', dfflt, 
     1     1, xlat(nlin,1))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Lower Left Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Lower Left Longitude', dfflt, 
     1     1, xlon(nlin,1))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Lower Left Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Lower Right Latitude', dfflt, 
     1     1, xlat(nlin,3))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Lower Right Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Lower Right Longitude', dfflt, 
     1     1, xlon(nlin,3))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Lower Right Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Start Center Latitude', dfflt, 
     1     1, xlat(1,2))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Start Center Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Start Center Longitude', dfflt, 
     1     1, xlon(1,2))      
      if (iret.eq.-1) then
        write(*,*)'Error writing Start Center Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'End Center Latitude', dfflt, 
     1     1, xlat(nlin,2))      
      if (iret.eq.-1) then
        write(*,*)'Error writing End Center Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'End Center Longitude', dfflt, 
     1     1, xlon(nlin,2))      
      if (iret.eq.-1) then
        write(*,*)'Error writing End Center Longitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Northernmost Latitude', dfflt, 
     1     1, xlanth)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Northernmost Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Southernmost Latitude', dfflt, 
     1     1, xlasth)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Southernmost Latitude'
        return
      endif
 
      iret = sfsattr(prod_ID, 'Westernmost Longitude', dfflt, 
     1     1,  xlowst)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Westernmost Longitude'
        return
      endif
 
      at_id = sfsattr(prod_ID, 'Easternmost Longitude', dfflt, 
     1     1,  xloest)      
      if (iret.eq.-1) then
        write(*,*)'Error writing Easternmost Longitude'
        return
      endif
 
      return
      end

      subroutine put_scan_metadata (prod_ID, sdname, idim, sddata, 
     1     iret )

      real*4 sddata(*)
      integer*4 prod_ID, idim, iret
      integer*4 istart(2), istr(2), idims(2), ind, sd_id
      integer sfn2index, sfselect, sfwdata, sfendacc
      character*(*) sdname
      data istart/2*0/, istr/2*1/

      ind = sfn2index(prod_ID, sdname)
      if (ind.eq.-1) then
         iret = -1
         write(*,*) 'Error getting sds index'
         return
      end if

      sd_id = sfselect(prod_ID, ind)
      if (sd_id.eq.-1) then
         iret = -1
         write(*,*) 'Error selecting sds'
         return
      end if

      idims(1) = idim
      iret = sfwdata(sd_id, istart, istr, idims, sddata)
      if (iret.eq.-1) then
        write(*,*) 'Error writing sds'
        return
      end if

      iret = sfendacc(sd_id)
      
      return
      end
