      subroutine fitgps(ngps,gps,nsig,vecs,driv,s0,updorb)
c
c  Purpose:  This routine uses the GPS, the interpolated ASAP vectors and
c               the partial derivatives to compute updates to the orbital
c               elements.  It uses a weighted batch least-squares algorithm
c               to minimize the orbit position differences and also performs 
c               data rejection based on an initial fit to the data.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  ngps         I*4      I      Number of data samples
c  gps(6,ngps)  R*8      I      GPS orbit vectors
c  nsig(ngps)   R*8      I      Number of GPS signals per sample
c  vecs(6,ngps) R*8      I      ASAP orbit vectors
c  driv(6,3,ngps) R*8    I      Partial derivatives of orbit position with
c                                respect to initial orbital elements
c  s0(6)        R*8      I      State weights to use in least-squares   
c  updorb(6)    R*8      O      Updates to orbital elements
c
c  By: Frederick S. Patt, GSC, December 23, 1993
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "nav_cnst.fin"

      real*8 gps(6,maxlin), vecs(6,maxlin), driv(6,3,maxlin), updorb(6)
      real*8 pd(3,maxlin), pe(3,maxlin), up(6), w(6,6), v(6), s0(6)
      real*8 xtol, xmed, xiqr
      integer*4 nsig(maxlin), ngps, nobs, n, min, i, j, k, l, ier

      nobs = 3*ngps

c  Compute position differences
      do i=1,ngps
        do j=1,3
          pd(j,i) = gps(j,i) - vecs(j,i)

c  Scale derivatives for eccentricity and argument of perigee to 
c   avoid numerical roundoff errors in matrix inversion
          driv(2,j,i) = driv(2,j,i)/1.d2
          driv(5,j,i) = driv(5,j,i)*1.d2
        end do
      end do

c  Set tolerance for points with large differences
      call mediqr(pd,nobs,xmed,xiqr)
      xtol = 2.d0*xiqr
      min = 5

c  Zero out state matrix and vector
      do k=1,6
        v(k) = 0.d0
        do l=1,6
          w(k,l) = 0.d0      
        end do
      end do

c  Perform initial estimate of correction

c  Generate state equations
      n = 0
      do i=1,ngps
        if((nsig(i).ge.min).and.
     *    (abs(pd(1,i)-xmed).lt.xtol).and.
     *    (abs(pd(2,i)-xmed).lt.xtol).and.
     *    (abs(pd(3,i)-xmed).lt.xtol)) then
          n = n + 1
          do j=1,3
            do k=1,6
              v(k) = v(k) + pd(j,i)*driv(k,j,i)
              do l=1,6
                w(k,l) = w(k,l) + driv(k,j,i)*driv(l,j,i)
              end do
            end do
          end do
        end if
      end do
      print *,'FITGPS:',n,xmed,xiqr
    
c  Add state weights
      do k=1,6
        w(k,k) = w(k,k) + s0(k)
      end do
  
c  Solve for state update
      call invert(w,v,6,6,up,ier)
      print *,up

c  Recompute fitting residuals
      do i=1,ngps
        do j=1,3
          pe(j,i)=pd(j,i)
          do k=1,6
            pe(j,i) = pe(j,i) - up(k)*driv(k,j,i)
          end do
        end do
      end do
        
c  Compute updated tolerance
      call mediqr(pe,nobs,xmed,xiqr)
      xtol = 2.d0*xiqr
      min = 4

c  Zero out state matrix and vector
      do k=1,6
        v(k) = 0.d0
        do l=1,6
          w(k,l) = 0.d0      
        end do
      end do

c  Perform final estimate of correction

c  Generate state equations
      n = 0
      do i=1,ngps
        if((nsig(i).ge.min).and.
     *    (abs(pe(1,i)-xmed).lt.xtol).and.
     *    (abs(pe(2,i)-xmed).lt.xtol).and.
     *    (abs(pe(3,i)-xmed).lt.xtol)) then
          n = n + 1
          do j=1,3
            do k=1,6
              v(k) = v(k) + pd(j,i)*driv(k,j,i)
              do l=1,6
                w(k,l) = w(k,l) + driv(k,j,i)*driv(l,j,i)
              end do
            end do
          end do
        end if
      end do
      print *,'FITGPS:',n,xmed,xiqr
            
c  Add state weights
      do k=1,6
        w(k,k) = w(k,k) + s0(k)
      end do
  
c  Solve for state update
      call invert(w,v,6,6,up,ier)
      print *,up

c  Scale corrections for eccentricity and argument of perigee
      up(2) = up(2)/1.d2
      up(5) = up(5)*1.d2

c  Load corrections into output array
      
      do k=1,6
        updorb(k) = up(k)
      end do

      return
      end
