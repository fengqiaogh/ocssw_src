      subroutine asap_rot_int(nstp,iyinit,idinit,tsap,asap,ngps,
     * igyr,igday,gpsec,vecs)
c
c  Purpose: Rotate ASAP orbit vectors to ECEF coordinates and interplolate
c   		to GPS sample times
c  
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nstp		I*4	 I	Number of ASAP vectors
c  iyinit	I*4	 I	ASAP time tag year
c  idinit	I*4	 I	ASAP time tag day
c  tsap(nstp)	R*8	 I	ASAP time tag seconds-of-day
c  asap(6,nstp)	R*8	 I	ASAP orbit vectors (inertial)
c  ngps		I*4	 I	Number of GPS sample times
c  igyr		I*4	 I	GPS time tag year
c  igday	I*4	 I	GPS time tag year
c  gpsec(nstp)	R*8	 I	GPS time tag seconds-of-day
c  vecs(6,nstp)	R*8	 I	Interpolated ASAP vectors
c
c
c  By: Frederick S. Patt, GSC, December 22, 1993
c
c  Notes:  Method uses cubic polynomial to match positions and velocities
c   at input data points.
c
c  Modification History:
c
c  Corrected logic to allow for times out of ascending order
c
c  Added error return code
c  F. S. Patt, July 19, 2011

      implicit none
c
      real*8 asap(6,*),tsap(*),vecs(6,*),gpsec(*)
      real*8 a0(3),a1(3),a2(3),a3(3),vt(6,2)
      real*8 tsap1,tsap2,t,dift,dt,x,x2,x3
      integer*4 nstp,iyinit,idinit,ngps,igyr,igday,jd
      integer*4 ind, i, j, i1, nr
      data nr/2/

      ind = 1
      tsap2 = -1.d6
      tsap1 = 1.d6
      
      dift = (jd(igyr,1,igday) - jd(iyinit,1,idinit))*864.d2

c  Start main interpolation loop      
      do i=1,ngps
         t = dift + gpsec(i)

c  Check if GPS time is outside range of current ASAP vectors
         if ((t.gt.tsap2).or.(t.lt.tsap1)) then
            if (t.gt.tsap(ind+1)) then
               do while (t.gt.tsap(ind+1))
                  ind = ind + 1
                  if (ind.ge.nstp) then
                     print *, 'GPS times after available ASAP data'
                     do j=1,6
                        vecs(j,i) = 0.d0
                     end do
                     ind = 1
                     go to 990
                  end if
               end do
            else
               do while (t.lt.tsap(ind))
                  ind = ind - 1
                  if (ind.lt.1) then
                     print *, 'GPS times before available ASAP data'
                     do j=1,6
                        vecs(j,i) = 0.d0
                     end do
                     ind = 1
                     go to 990
                  end if
               end do
            end if
   
c  Set up cubic interpolation
            dt = tsap(ind+1) - tsap(ind)
           call asap_rots(iyinit,idinit,tsap(ind),asap(1,ind),nr,vt)
           do j=1,3
              a0(j) = vt(j,1)
              a1(j) = vt(j+3,1)*dt
              a2(j) = 3.d0*vt(j,2) - 3.d0*vt(j,1) 
     *             - 2.d0*vt(j+3,1)*dt - vt(j+3,2)*dt
              a3(j) = 2.d0*vt(j,1) - 2.d0*vt(j,2) 
     *             + vt(j+3,1)*dt + vt(j+3,2)*dt
           end do
           tsap1 = tsap(ind)
           tsap2 = tsap(ind+1)
        end if
      
c     Interpolate orbit position and velocity components to GPS time
        x = (t - tsap1)/dt
        x2 = x*x
        x3 = x2*x
        do j=1,3
           vecs(j,i) = a0(j) + a1(j)*x + a2(j)*x2 + a3(j)*x3
           vecs(j+3,i) = (a1(j) + 2.*a2(j)*x + 3.*a3(j)*x2)/dt
        end do
      
      end do
      return

 990  call EXIT(1)

      end
