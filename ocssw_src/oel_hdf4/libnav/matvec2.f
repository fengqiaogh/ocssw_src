        subroutine matvec2(xm, nxm, v1, v2)
c
c  matvec2(xm, nxm, v1, v2)
c
c  Purpose: create v2 as the product of matrix xm and vector v1
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  xm           R*8      I      size nxm by nxm matrix
c  nxm          I*4      I      size of matrix and vectors
c  v1           R*8      I      size nxm Input Vector
c  v2           R*8      O      size nxm Output Vector
c
c  By: W. Robinson, GSC,  2 Apr 93 from matvec by F Patt
c
c  Notes:  
c
c  Modification History:
c
      implicit none
      integer*4 nxm
      real*8 v1(nxm), v2(nxm), xm(nxm,nxm)
c
      integer*4 i, j
c
      do i=1, nxm
         v2(i) = 0.0
         do j=1, nxm
            v2(i) = v2(i) + xm(i,j)*v1(j) 
         end do
      end do
      return
      end
