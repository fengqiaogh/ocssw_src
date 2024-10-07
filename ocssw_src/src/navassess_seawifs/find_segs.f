      subroutine find_segs(ill,ilp,ilf,i1,imax,jll,jlp,iln,nln)

c  This subroutine finds additional segments in the list 
c  adjacent to the current segment

c  Calling Arguments

c  Name         Type    I/O     Description
c
c  ill(*)       I*4      I      Array of segment line numbers
c  ilp(2,*)     I*4      I      Array of segment endpoint pixel numbers
c  ilf(*)       I*4      I      Array of flags (1=segment already checked)
c  i1           I*4      I      First segment to check in list
c  imax         I*4      I      Number of segments in list
c  jll          I*4      I      Line number of current segment
c  jlp(2)       I*4      I      End pixel numbers of current segment
c  iln(10)      I*4      O      Indices of continguous segments
c  nln          I*4      O      Number of contiguous segments found

c       Subprograms Called:

c       Program written by:     Frederick S. Patt
c                               General Sciences Corporation
c                               April 17, 1994
c
c       Modification History:   
     


      integer*4 ill(*),ilp(2,*),ilf(*),jll,jlp(2),iln(10)


      nln = 0
      do i=i1,imax
        if (ilf(i).eq.0) then
          if ((ill(i)-jll).gt.1) then
            return
          else if ((ilp(1,i).le.(jlp(2)+1)).and.
     *          (ilp(2,i).ge.(jlp(1)-1)).and.
     *          (abs((ill(i)-jll)).eq.1)) then
            nln = nln + 1
            iln(nln) = i 
          end if
        end if
      end do

      return
      end
