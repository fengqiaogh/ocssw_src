      subroutine get_l1a_open( infile, amode, prod_ID, isday, 
     *     nlin, iret )

c  get_l1a_open( infile, amode, prod_ID, nlin, iret )
c
c  Purpose: opens a L1A HDF file and reads the number of scan lines attribute
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     infile       char     I      size 80 L1A file name string
c     amode        I*4      I      mode for HDF file open
c                                    =1, read-only
c                                    =3, read/write
c     prod_ID      I*4      O      HDF file ID
c     isday        I*4      O      Start day of file 
c     nlin         I*4      O      number of scan lines in the L1A file
c     iret         I*4      O      return code
c                                    =0, success
c                                    =-1, failure    
c
c  By: F. S. Patt, SAIC GSC, 24 Sep 98
c
c  Notes:  
c
c  Modification History:
c
c  Modified to read and return 'Start Day' attribute from file.
c  F. S. Patt, SAIC GSC, April 22, 1999.
c
c  Modified to read start year and convert start day to Julian day, 
c  to solve year-rollover problem.
c  F. S. Patt, SAIC, November 5, 2002.


      character*80 infile
      integer*4 prod_ID, at_id
      integer*4 nlin
      integer*4 amode
      integer*4 isday,isyr
      integer*2 iday,iyr
      integer sfstart, sffattr, sfrattr

      iret = 0

      prod_ID = sfstart(infile, amode)
      if (prod_ID.eq.-1) then
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

      at_id = sffattr(prod_ID, 'Start Day')
      if (at_id.eq.-1) then
        write(*,*)'Error getting index for Start Day'
        iret = -1
        return
      endif

      iret = sfrattr(prod_ID, at_id, iday)      
      if (iret.eq.-1) then
        write(*,*)'Error getting value for Start Day'
        return
      endif

      isday = iday

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

      isyr = iyr

      isday = jd(isyr,1,isday)

      return
      end


