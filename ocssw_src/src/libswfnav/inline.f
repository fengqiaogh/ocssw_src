        subroutine inline(lun,icomm,line,iret)
c
c  inline(lun,icomm,line,iret)
c
c  Purpose: read the next line from a file that is not a comment
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  lun          I*4      I      Unit # of file to read (already open) 
c  icomm        char*1   I      Comment character '!' or otherwise
c  line         char*80  O      returned line
c  iret         I*4      O      return code, 0 = good
c
c  By: W. Robinson, GSC, 22 Mar 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
c
      integer*4 lun, iret, istat
      character icomm*1, line*80
      logical end
c
c
c       start, read in the line
c
      iret = 0
      end = .false.
      do while( .not. end )
         read(lun,100,iostat=istat)line
  100    format(a)
         if( istat .ne. 0 ) then
            iret = -1
            end = .true.
         else
            if( line(1:1) .ne. icomm(1:1) ) end = .true.
         end if
      end do
c
c       and end
c
      return
      end
