       subroutine mediqr(dat,ndat,xmed,xiqr)
c
c  Purpose:  This routine computes the median and interquartile range
c               for an input data array.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  dat          R*8      I      Input data array
c  ndat         I*4      I      Number of items in input data array
c  xmed         R*8      O      Median of data
c  xiqr         R*8      O      Interquartile range of data
c
c  By: Frederick S. Patt, GSC, January 4, 1994
c
c  Notes:  Uses Select routine from Numerical Recipes for FORTRAN, p334.
c
c  Modification History:
c
      implicit none

#include "nav_cnst.fin"
      real*8 dat(3*maxlin),xmed,xiqr
      real*8 da1(3*maxlin),select
      integer*4 ndat,i,iq,ir,im

c  Initialize index array
      do i=1,ndat
        da1(i) = dat(i)
      end do

c  Compute median and IQR by calling function Select

      im = ndat/2 + 1
      iq = (ndat-1)/4 + 1
      ir = ndat - iq + 1
      xmed = select(im,ndat,da1)
      xiqr = select(ir,ndat,da1) - select(iq,ndat,da1)

      return
      end
