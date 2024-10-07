      subroutine eartol(earth,navqc,earrng)
c
c  eartol(earth,navqc,earrng)
c
c  Purpose: Check earth sensor values to see that they are within 
c           the required tolerences
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  earth        struct  I/O     earth sensor data structure
c  navqc        struct   I      navigation quality control info
c  earrng       I*4      I      size 2 (start-end) by 2 ( sens 1, 2 ) 
c                               array of active range for
c                               the 2 earth sensors
c
c  By: W. Robinson, GSC,  13 Apr 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
c
      type(earth_struct) :: earth(2)
      type(navqc_struct) :: navqc
c
      integer*4 earrng(2,2)
c
      integer*4 isens, ilin
c
c
c       flag any unflagged earth sensor data in active range with
c       angles out of tolerence
c
      do isens = 1,2
        if( earrng(1,isens) .ne. -1 ) then
c
          do ilin = earrng(1,isens), earrng(2,isens)
            if( earth(isens)%flag(ilin) .eq. 0 ) then
              if(( earth(isens)%widphse(1,ilin) .lt. 
     1                           navqc%ear_tol_wd(1) ) .or.
     1           ( earth(isens)%widphse(1,ilin) .gt. 
     1                           navqc%ear_tol_wd(2) ) .or.
     1           ( earth(isens)%widphse(2,ilin) .lt. 
     1                           navqc%ear_tol_ph(1) ) .or.
     1           ( earth(isens)%widphse(2,ilin) .gt. 
     1                           navqc%ear_tol_ph(2) )       )
     1             earth(isens)%flag(ilin) = 1
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
