      subroutine asap_int2( nstp, tsap, asap, pos_erri, ngps,
     * gpsec, vecs, pos_erro )
c
c  Purpose: Rotate ASAP orbit vectors to ECEF coordinates and interplolate
c               to GPS sample times
c  
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  nstp         I*4      I      Number of ASAP vectors
c  tsap(nstp)   R*8      I      ASAP time tag seconds-of-day
c  asap(6,nstp) R*8      I      ASAP orbit vectors (inertial)
c  pos_erri(nstp) R*8    I      position error , meters 
c  ngps         I*4      I      Number of GPS sample times
c  gpsec(ngps)  R*8      I      GPS time tag seconds-of-day
c  vecs(6,ngps) R*8      O      Interpolated ASAP vectors
c  pos_erro(ngps) R*4    O      interpolated position errors (linearly)
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
c  W. Robinson, SAIC, 19 Dec 2005  adapted for use with CZCS daily orbit data
c    and to linearly interpolate the position error information
c
      implicit none
c
      real*8 asap(6,*),tsap(*),vecs(6,*),gpsec(*)
      real*8 a0(3),a1(3),a2(3),a3(3),vt(6,2)
      real*8 tsap1,tsap2,t,dift,dt,x,x2,x3
      real*8 pos_erri(*)
      real*4 pos_erro(*)
      integer*4 nstp,iyinit,idinit,ngps,igyr,igday,jd
      integer*4 ind, i, j, i1, nr
      data nr/2/

      ind = 1
      tsap2 = -1.d6
      tsap1 = 1.d6
      
      dift = 0.

c  Start main interpolation loop      
      do i=1,ngps
        t = dift + gpsec(i)

c  Check if GPS time is outside range of current ASAP vectors
        if ((t.gt.tsap2).or.(t.lt.tsap1)) then
          if (t.gt.tsap(ind+1)) then
            dowhile (t.gt.tsap(ind+1))
              ind = ind + 1
              if (ind.ge.nstp) then
                write(*,*) 'GPS times after available ASAP data'
                do j=1,6
                  vecs(j,i) = 0.d0
                end do
                ind = 1
                pos_erro(i) = -1
                go to 990
              end if
            end do
          else
            dowhile (t.lt.tsap(ind))
              ind = ind - 1
              if (ind.lt.1) then
                write(*,*) 'GPS times before available ASAP data'
                do j=1,6
                  vecs(j,i) = 0.d0
                end do
                ind = 1
                pos_erro(i) = -1
                go to 990
              end if
            end do
          end if
            

c  Set up cubic interpolation
          dt = tsap(ind+1) - tsap(ind)
c         call asap_rots(iyinit,idinit,tsap(ind),asap(1,ind),nr,vt)
          do j=1,6
             vt(j,1) = asap(j,ind)
             vt(j,2) = asap(j,ind+1)
          end do
          do j=1,3
            a0(j) = vt(j,1)
            a1(j) = vt(j+3,1)*dt
            a2(j) = 3.d0*vt(j,2) - 3.d0*vt(j,1) 
     *              - 2.d0*vt(j+3,1)*dt - vt(j+3,2)*dt
            a3(j) = 2.d0*vt(j,1) - 2.d0*vt(j,2) 
     *              + vt(j+3,1)*dt + vt(j+3,2)*dt
          end do
          tsap1 = tsap(ind)
          tsap2 = tsap(ind+1)
        end if

c  Interpolate orbit position and velocity components to GPS time
c  also, just linearly interpolate the position error values
        x = (t - tsap1)/dt
        x2 = x*x
        x3 = x2*x
        do j=1,3
          vecs(j,i) = a0(j) + a1(j)*x + a2(j)*x2 + a3(j)*x3
          vecs(j+3,i) = (a1(j) + 2.*a2(j)*x + 3.*a3(j)*x2)/dt
        end do

        pos_erro(i) = pos_erri(ind) * ( 1. - x ) + pos_erri(ind+1) * x

  990   continue
      end do
      return

      end       
