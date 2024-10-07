      subroutine actrng(nitem, active, irange)
c
c  actrng(nitem, active, irange)
c
c  Purpose: find the active range of a set of flags
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nitem        I*4      I      number of items in active array
c  active       I*4      I      array of active flage, 0 = active
c  irange       I*4      O      size 2 start, end of active range
c
c  By: W. Robinson, GSC,  1 Apr 93
c
c  Notes: this will find the extremes of the active range, ie, from
c         the first active to the last active values 
c
c  Modification History:
c
      implicit none
c
      integer*4 nitem, active(nitem), irange(2), ilin
c
c
c       start, set range to bad, loop and find first and last
c
      irange(1) = -1
      irange(2) = -1
c
      do ilin = 1, nitem
c
         if( irange(1) .eq. -1  .and.  
     1     active(ilin) .eq. 0 ) irange(1) = ilin
c
         if( active(ilin) .eq. 0 ) irange(2) = ilin
c
      end do
c
      return
      end
