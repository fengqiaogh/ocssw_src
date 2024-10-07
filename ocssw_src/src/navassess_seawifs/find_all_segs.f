      subroutine find_all_segs(ill,ilp,ilf,ist,nsegs,jll,jlp,nj)
    
c  This subroutine finds all contiguous segments in a list of segments
c  from a specified starting point.  For each contiguous segment found, 
c  the list is searched for additional contiguous segments until no more
c  are found or the list is exhausted.

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  ill(*)       I*4      I      Array of segment line numbers
c  ilp(2,*)     I*4      I      Array of segment endpoint pixel numbers
c  ilf(*)       I*4      I      Array of flags (1=segment already checked)
c  ist          I*4      I      Starting point in list
c  nsegs        I*4      I      Number of segments in list
c  jll(*)       I*4     I/O     Array of contiguous segment line numbers
c  jlp(2,*)     I*4     I/O     Array of contiguous segment end pixel numbers
c  nj           I*4      O      Number of contiguous segments found

c       Subprograms Called:

c       find_segs       Find segments in list adjacent to one segment


c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               April 17, 1994
c
c       Modification History:   
     
      integer*4 ill(*),ilp(2,*),ilf(*)
      integer*4 jll(*),jlp(2,*)
      integer*4 iln(10)
      logical end

      ij = 1
      nj = 1

      dowhile (ij.le.nj)

c  Find segments contiguous to each segment
      call find_segs(ill,ilp,ilf,ist,nsegs,jll(ij),jlp(1,ij),iln,nln)

c  If more continguous segments found, add to list
      if(nln.gt.0) then
        do i=1,nln
          ilf(iln(i)) = 1
          nj = nj + 1
          jll(nj) = ill(iln(i))
          jlp(1,nj) = ilp(1,iln(i))
          jlp(2,nj) = ilp(2,iln(i))
        end do
      end if
      ij = ij + 1

      end do

      return
      end       
