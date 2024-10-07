      subroutine get_l1a_open( infile, prod_ID, npix, nlin, iret )

c  This routine opens the Level 1A HDF file and reads the attributes

      character*80 infile,bname
      integer*4 prod_ID(2), at_id
      integer*4 npix, nlin 
      integer*4 amode, jd, jdr
      integer*2 iyr, idy
      integer sfstart, sffattr, sfrattr
      data amode/1/, jdr/2450696/
                        
      prod_ID(1) = sfstart(infile, amode)
      if (prod_ID(1).eq.-1) then
        write(*,*)'Error opening HDF file'
        iret = -1
        return
      endif     

      at_id = sffattr(prod_ID, 'Number of Scan Lines')
      if (at_id.eq.-1) then
        write(*,*)'Error getting index for Number of Scan Lines'
        iret = -1
        return
      endif

      iret = sfrattr(prod_ID, at_id, nlin)      
      if (iret.eq.-1) then
        write(*,*)'Error getting value for Number of Scan Lines'
        return
      endif

      at_id = sffattr(prod_ID, 'Pixels per Scan Line')
      if (at_id.eq.-1) then
        write(*,*)'Error getting index for Pixels per Scan Line'
        iret = -1
        return
      endif

      iret = sfrattr(prod_ID, at_id, npix)
      if (iret.eq.-1) then
        write(*,*)'Error getting value for Pixels per Scan Line'
        return
      endif

      at_id = sffattr(prod_ID, 'Start Year')
      if (at_id.eq.-1) then
        write(*,*)'Error getting index for Start Year'
        iret = -1
        return
      endif

      iret = sfrattr(prod_ID, at_id, iyr)
      if (iret.eq.-1) then
        write(*,*)'Error getting value for Start Year'
        return
      endif

      at_id = sffattr(prod_ID, 'Start Day')
      if (at_id.eq.-1) then
        write(*,*)'Error getting index for Start Day'
        iret = -1
        return
      endif

      iret = sfrattr(prod_ID, at_id, idy)
      if (iret.eq.-1) then
        write(*,*)'Error getting value for Start Day'
        return
      endif

      lyr = iyr
      ldy = idy

      prod_ID(2) = jd(lyr, 1, ldy) - jdr

      return
      end


