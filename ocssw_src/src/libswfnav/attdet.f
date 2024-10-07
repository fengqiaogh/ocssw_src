        subroutine attdet( navblk, navctl, nad_bod, attxfm, 
     1     pos, vel, widphse1, widphfl1, widphse2, widphfl2, 
     1     sun_bod, sunbodfl, att_prev, delta_t, covar)
c
c  attdet( navblk, navctl, nad_bod, attxfm, 
c       pos, vel, widphse1, widphfl1, widphse2, 
c       widphfl2, sun_bod, sunbodfl, att_prev, delta_t, covar )
c
c  Purpose: perform the actual attitude determination for one line
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  navblk       struct  I/O     navigation block including:
c                               in: sun_ref - sun ref in geo coord sys
c                                   lcl_vert - local vertical
c                               out: sen_frame - attitude matrix
c                                    att_ang - attitude angles in the
c                                       order yaw, roll, pitch
c  navctl       struct   I      navigation controls. Used are:
c                               redoyaw - # iterations to redo yaw 
c                               computation
c                               yawtol - difference between yaw guess
c                               and computed yaw 
c  nad_bod      R*4      I      size 3 Unit Nadir Vector in Spacecraft Frame
c  attxfm       R*4      I      size 3 by 3 attitude transform matrix - 
c                               the transform from geocentric to 
c                               orbit coordinates
c  pos          R*4      I      size 3 orbit Position Vector (km)
c  vel          R*4      I      size 3 orbit Velocity Vector (km/sec)
c  widphse1     R*4      I      size 2 processed earth width 
c                               and phase for sensor 1
c  widphfl1     I*4      I      flag for earth sensor 1 - 0 = good
c  widphse2     R*4      I      size 2 processed earth width 
c                               and phase for sensor 2
c  widphfl2     I*4      I      flag for earth sensor 2 - 0 = good
c  sun_bod      R*4      I      size 3 weighted sun vectors in Spacecraft 
c                               frame
c  sunbodfl     I*4      I      flag for Sun sensor data - 0 = good
c  att_prev     R*4      I      size 3 attitude angles from previous scan line
c  delta_t      R*4      I      time difference from previous scan line
c                               (set to 0 if first line in frame, negative if
c                               propagating backwards in time)
c  covar        R*4      I      size 3 array of covariance matrix diagonal 
c                               elements
c
c
c  By: W. Robinson, GSC, 14 Apr 93
c
c  Notes: the initial yaw used to create the nadir vector from the
c       earth sensors is assumed to be the previous value.  If the 
c       resultant yaw computed here id different from the guess by 
c       the value yawtol, and the redoyaw is set, the nadir vector 
c       will be re-computed by calling routine earth again.  This will
c       be repeated up to redoyaw number of times 
c
c  Modification History:
c
c  Modified extraction of Euler angles from attitude matrix to reflect
c  actual order and sign of rotations implemented by OSC in the ACS.
c  The order is 2-1-3 in the ACS (not the Spacecraft) frame.  The ACS
c  axes defined as:  X toward the velocity (spacecraft -Y axis); 
c  Y toward the negative orbit normal (spacecraft -Z axis); Z toward 
c  the nadir (spacecraft X axis).  This is equivalent to a 
c  pitch-roll-yaw transformation, with the pitch and roll angles being 
c  negative about the spacecraft Z and Y axes.
c  F. S. Patt, GSC, October 4, 1996. 
c
c  Modified the order of operations to first transform the matrix constructed
c  from the reference vectors to the orbital frame, then use that matrix to
c  compute the orbit-to-spacecraft transformation.  This will give the same
c  results as before, but allows the attitude matrix to be a small rotation.
c  F. S. Patt, SAIC GSC, February 3, 1998. 
c
c  Modified attitude angle calculation to weight the yaw angle according to
c  the Sun/Nadir angle to prevent blowing up at subsolar point.
c  F. S. Patt, SAIC GSC, March 9, 1998.
c  
c  Rewrote this routine to implement attitude determination using a Kalman
c  filter instead of the deterministic algorithm used previously.  This was
c  done to resolve multiple issues:  reduce jitter at sensor transitions,
c  improve yaw accuracy at the subsolar point and handle sensor gaps.  The
c  initial version uses a simple dynamics model which assumes an inertial
c  spin axis (not necessarily parallel to the orbit normal).
c  All of the equations are taken from Wertz, Section 13.5.2, pp. 462-465.
c  F. S. Patt, SAIC GSC, June 1, 1998.
c
c  Reduced the state noise value (from 1.e-8 to 0.25e-8) to reduce the yaw
c  variation in the middle of the orbit.
c  F. S. Patt, SAIC GSC, June 22, 1998.
c
c  Corrected a (minor) error in the calculation of the partial derivatives.
c  F. S. Patt, SAIC GSC, October 19, 1998.
c 
c  Reduced the state noise value for roll and yaw to 6.e-10, to constrain the 
c  yaw angle near the subsolar point.  F. S. Patt, SAIC GSC, November 25, 1998.
c
c  Added data initialization statements for p and q to ensure they get saved 
c  between calls.   F. S. Patt, SAIC GSC, April 14, 1999
c
c  Modifed calculation of Sun vector measurement covariance to avoid 
c  numerical round-off errors for Sun vector very close to -X axis.
c  F. S. Patt, SAIC GSC, August 27, 1999
c
c  Increased state noise for roll and yaw to 2.5e-9 based on results of testing
c  with revised horizon scanner calibration and Sun sensor alignments.
c  F. S. Patt, SAIC GSC, September 20, 1999
c
c  Modified state noise calculation to use delta_t**2 (instead of abs(delta_t))
c  to reduce noise for HRPT/LAC, and changed qvar (from 2.5e-9 to 4.0e-9) to 
c  offset this for GAC; this is based on initial results of testing for the 
c  mergered 1-km (MLAC) files.
c  F. S. Patt, SAIC, November 22, 2002

      implicit none
#include "nav_cnst.fin"
#include "navblk_s.fin"
#include "navctl_s.fin"
c
      type(navblk_struct) :: navblk   ! note we only want the block for 1 line
      type(navctl_struct) :: navctl
      integer*4 widphfl1, widphfl2, sunbodfl
      real*4 nad_bod(3), sun_bod(3), attxfm(3,3), widphse1(2), 
     1     widphse2(2), pos(3), vel(3), att_prev(3), delta_t, covar(3)
c
      real*8 xk(3,6), p(3,3), r(6,6), g(6,3), xid(3,3), d(3,3), q(3,3)
      real*8 t33(3,3), t66(6,6), t36(3,6), yg(6), gt(3,6), dt(3,3) 
      real*8 t6(6), t3(3), r33(3,3), s33(3,3), xx(6,3)
      real*4 dmda(3,3,3), dvda(3,3), cpit, spit, crol, srol, cyaw, syaw
      real*4 att_ap(3), apm(3,3), sunvar, earvar, qvar(3), vmag, tv(3)
      integer*4 nyaw, ier, i, j, nadbodfl
      real*4 pmag, yawest, sr(3), nr(3)
      logical recomp, update

c   Set variances for Sun and Earth vectors and state noise
      data sunvar/1.e-6/, earvar/3.e-6/, qvar/2*4.0e-9,4.0e-9/
c   Identity matrix
      data xid/1.d0,3*0.d0,1.d0,3*0.d0,1.d0/
      data d/1.d0,3*0.d0,1.d0,3*0.d0,1.d0/
      data p/1.d0,3*0.d0,1.d0,3*0.d0,1.d0/
      data q/9*0.d0/
      data r/36*0.d0/
    
c   Compute state propagation and noise covariance matrices
      d(2,1) = delta_t*.00106d0
      d(1,2) = -d(2,1)
      do i=1,3
c         q(i,i) = abs(delta_t)*qvar(i)
         q(i,i) = delta_t*delta_t*qvar(i)
      end do

c   If first scan of scene (delta_t = 0), set state covariance and attitude
      if (delta_t.eq.0.0) then
         do i=1,3
            att_ap(i) = 0.0
            do j=1,3
               p(i,j) = xid(i,j)
               apm(i,j) = xid(i,j)
            end do
         enddo
      else

c   Compute a priori attitude and matrix based on previous attitude

c      Propagate attitude using simple dynamics model (Wertz, eq. 13-86)
c      May replace this with a more sophisticated model at a later date
c      Need to do this in-line to handle mixed precision
         do i=1,3
            att_ap(i) = d(i,1)*att_prev(1) + d(i,2)*att_prev(2)
     *           + d(i,3)*att_prev(3)
         end do
         call euler(att_ap, apm)

c      Propagate state covariance matrix and add state noise (Wertz, eq. 13-88)
         call dmatmp(d,p,t33,3,3,3)
         call dxpose(d,dt,3,3)
         call dmatmp(t33,dt,p,3,3,3)
         do i=1,3
            do j=1,3
               p(i,j) = p(i,j) + q(i,j)
            end do
         end do
      end if
      
c      call matmpy(apm, attxfm, m)

c   Compute reference vectors in a prior spacecraft frame
      pmag = vmag(pos)
c      call matvec(m, pos, nr)
c      do i=1,3 
c         nr(i) = -nr(i)/pmag
c      end do
c      call matvec(m, navblk%sun_ref(1), sr)

c   Compute partial derivatives of attitude matrix with respect to angles
      cyaw = cos(att_ap(1)/radeg)
      syaw = sin(att_ap(1)/radeg)
      crol = cos(att_ap(2)/radeg)
      srol = sin(att_ap(2)/radeg)
      cpit = cos(att_ap(3)/radeg)
      spit = sin(att_ap(3)/radeg)
      do i=1,3
         dmda(1,i,1) = 0.0
         dmda(2,i,1) = apm(3,i)
         dmda(3,i,1) = -apm(2,i)
         dmda(i,1,2) = -apm(i,3)*cpit
         dmda(i,2,2) = apm(i,3)*spit
         dmda(i,1,3) = apm(i,2)
         dmda(i,2,3) = -apm(i,1)
         dmda(i,3,3) = 0.0
      end do
      dmda(1,3,2) = crol
      dmda(2,3,2) = -syaw*srol
      dmda(3,3,2) = -cyaw*srol

c   Initialize for yaw iteration if necessary

      nyaw = 0
      recomp = .true.
      update = .false.
      yawest = att_ap(1)
c
      dowhile( recomp )
         call earth( pos, vel, widphse1, widphfl1, widphse2, widphfl2,
     1        yawest, navctl, nad_bod, nadbodfl, navblk%sun_ref(3) )

c          Initialize arrays
         do i=1,6
            yg(i) = 0.d0
            r(i,i) = 1.d0
            do j=1,3
               g(i,j) = 0.d0
               xk(j,i) = 0.d0
            end do
         end do   

c          Compute partial derivative matrix and difference vector

c          If there is valid Sun data
         if (sunbodfl.eq.0) then

c       Compute partials of Sun vector with respect to attitude angles 
            call matvec(attxfm, navblk%sun_ref(1), tv)

            do i=1,3
               call matvec(dmda(1,1,i),tv,dvda(1,i))
            end do
            call matvec(apm,tv,sr)

c          Load into G matrix and compute difference vector
            do i=1,3
               yg(i) = sun_bod(i) - sr(i)
               r(i,i) = sunvar*(1.0001d0 - sr(i)*sr(i))
c               r(i,i) = sunvar
               do j=1,3
                  g(i,j) = dvda(i,j)
               end do  
            end do   
         end if

c          If there is valid Earth data
         if (nadbodfl.eq.0) then

c          Compute partials of Earth vector with respect to attitude angles 
            call matvec(attxfm, pos,tv)
            do i=1,3 
               tv(i) = -tv(i)/pmag
            end do
            do i=1,3
               call matvec(dmda(1,1,i),tv,dvda(1,i))
            end do
            call matvec(apm,tv,nr)

c          Load into G matrix and compute difference vector
            do i=1,3
               yg(i+3) = nad_bod(i) - nr(i)
               r(i+3,i+3) = earvar
               do j=1,3
                  g(i+3,j) = dvda(i,j)
               end do  
            end do               
         end if

c      If there are observations, correct attitude and update covariance
         if ( (sunbodfl.eq.0) .or. (nadbodfl.eq.0) ) then

c          Compute gain matrix (Wertz, eq. 13-74)
            call dxpose(g,gt,6,3)
            call dmatmp(p,gt,t36,3,3,6)
            call dmatmp(g,t36,t66,6,3,6)
            do i=1,6
               t66(i,i) = t66(i,i) + r(i,i)
            end do
            call invert(t66,yg,6,6,t6,ier)
            call dmatmp(gt,t66,t36,3,6,6)
            call dmatmp(p,t36,xk,3,3,6)

c          Compute update to attitude angles
            call dmatmp(xk,yg,t3,3,6,1)
            do i=1,3
               navblk%att_ang(i) = att_ap(i) + radeg*t3(i)
            end do
            update = .true.

c      Else, use propagated attitude
         else
            do i=1,3
               navblk%att_ang(i) = att_ap(i)
            end do

         end if
c
c          if the re-compute yaw value is greater than the 
c          the tolerence and 
c
         nyaw = nyaw + 1
         if( navctl%redoyaw .ge. nyaw .and. 
     1      abs( navblk%att_ang(1) - yawest) .gt. navctl%yawtol ) then
            recomp = .true.
            yawest = navblk%att_ang(1)
         else
            recomp = .false.
         end if
      end do
      
c       Update covariance matrix (Wertz, eq. 13-76)
      call dmatmp(xk,g,t33,3,6,3)
      do i=1,3
         do j=1,3
            t33(i,j) = xid(i,j) - t33(i,j)
         end do
      end do
      call dmatmp(t33,p,r33,3,3,3)
      call dxpose(t33,s33,3,3)
      call dmatmp(r33,s33,p,3,3,3)
      call dmatmp(xk,r,t36,3,6,6)
      call dxpose(xk,xx,3,6)
      call dmatmp(t36,xx,t33,3,6,3)
      do i=1,3
         do j=1,3
            p(i,j) = p(i,j) + t33(i,j)
         end do
      end do

      if (navctl%lvdbug.gt.2) then
         write (66,*) ((att_ap(i) - att_prev(i)),i=1,3)
         write (66,*) ((navblk%att_ang(i) - att_ap(i)),i=1,3)
         write (66,*) yg
         do i=1,6
            t6(i) = 0.d0
            do j=1,3
               t6(i) = t6(i) + g(i,j)*t3(j)
            end do
         end do
         write (66,*) t6
      end if

c      Add diagonal elements of covariance matrix to output array
      do i=1,3
         covar(i) = p(i,i)
      end do

c
c       and end
c
      return
      end

        subroutine dmatmp(dm1, dm2, dm3, l, m, n)
c
c  dmatmp( dm1, dm2, dm3, l, m, n)
c
c  Purpose: multiply two general double-precision matrices
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     dm1          R*8      I       Input matrix of size l-by-m
c     dm2          R*8      I       Input matrix of size m-by-n
c     dm3          R*8      O       Output matrix of size l-by-n
c     l            I*4      I       Number of rows in dm1 and dm3
c     m            I*4      I       Number of columns in dm1 and rows in dm2
c     n            I*4      I       Number of columns in dm2 and dm3
c
c  By: F. S. Patt, SAIC GSC, June 2, 1998
c
c
c  Modification History:
c
      implicit none
      integer*4 l,m,n
      real*8 dm1(l,m),dm2(m,n),dm3(l,n)
      integer*4 i,j,k

c  Compute dm3 as dm1 X dm2
      do i=1,l
         do j=1,n
            dm3(i,j) = 0.d0
            do k=1,m
               dm3(i,j) = dm3(i,j) + dm1(i,k)*dm2(k,j)
            end do
         end do
      end do

      return
      end

        subroutine dxpose(din, dout, n, m)
c
c  dxpose( din, dout, n, m )
c
c  Purpose: form the transpose of matrix din to matrix dout
c
c  Calling Arguments:
c
c     Name         Type    I/O     Description
c     --------     ----    ---     -----------
c     din          R*8      I      size n by m input matrix
c     dout         R*8      O      size m by n output matrix
c     n            I*4      I      number of rows in din
c     m            I*4      I      number of columns in din 
c
c  By: F. S. Patt, SAIC GSC, June 2, 1998
c
c  Notes:  
c
c  Modification History:
c
      implicit none
      integer*4 n, m
      real*8 din(n,m), dout(m,n)
c
      integer*4 i, j
c
c
      do i = 1, n
        do j = 1, m
          dout(j,i) = din(i,j)
        end do
      end do
c
      return
      end
