        subroutine xpose( in, out )
c
c  xpose( in, out )
c
c  Purpose: form the transpose of matrix in to matrix out
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  in           R*4      I      size 3 by 3 input matrix
c  out          R*4      O      size 3 by 3 output matrix
c
c  By: W. Robinson, GSC, 15 Apr 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
      real*4 in(3,3), out(3,3)
c
      integer*4 i, j
c
c
      do i = 1, 3
        do j = 1, 3
          out(j,i) = in(i,j)
        end do
      end do
c
      return
      end
