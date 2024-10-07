      subroutine suntol(sun,navqc,sunrng)
c
c  suntol(sun,navqc,sunrng)
c
c  Purpose: Check sun sensor values to see that they are within 
c           the required tolerence
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  sun          struct  I/O     sun sensor data atructure
c  navqc        struct   I      navigation quality control info
c  sunrng       I*4      I      size 2 by 3 array of active range for
c                               the 3 sun sensors
c
c  By: W. Robinson, GSC,  1 Apr 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
      type(sun_struct) :: sun(3)
      type(navqc_struct) :: navqc
c
      integer*4 sunrng(2,3)
c
      integer*4 isens, ilin
c
c
c       flag any unflagged sun sensor data in active range with
c       angles out of tolerence
c
      do isens = 1,3
        if( sunrng(1,isens) .ne. -1 ) then
c
          do ilin = sunrng(1,isens), sunrng(2,isens)
            if( sun(isens)%flag(ilin) .eq. 0 ) then
              if(( sun(isens)%ang(1,ilin) .lt. navqc%sun_tol_1(1) ) .or.
     1           ( sun(isens)%ang(1,ilin) .gt. navqc%sun_tol_1(2) ) .or.
     1           ( sun(isens)%ang(2,ilin) .lt. navqc%sun_tol_2(1) ) .or.
     1           ( sun(isens)%ang(2,ilin) .gt. navqc%sun_tol_2(2) )    )
     1             sun(isens)%flag(ilin) = 1
            end if
          end do
c
        end if
      end do
c
c       and end
c
  990 continue
      return
      end
