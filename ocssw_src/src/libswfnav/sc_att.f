      subroutine sc_att(gaclac,tlm,navqc,att_ang,attangfl,
     1     iret)
c
c  sc_att(gaclac,tlm,navqc,att_ang,attangfl,iret)
c
c  Purpose: process spacecraft provided attitude to every scan line
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  gaclac       I*4      I      flag for GAC or LAC data.  If LAC
c                               data, a time exists for each scan 
c                               line, else only once every 5 lines
c                               and the lines are 4scan lines apart
c  tlm          struct   I      telemetry data structure containing
c                               the S/C attitude
c  navqc        struct   I      navigation quality control info
c  att_ang      R*4      O      3 by nlines array of spacecraft yaw,
c                               roll and pitch
c  attangfl     I*4      O      array of flags for the attitude values
c  iret         I*4      O      return code, 0 - good
c
c  By: W. Robinson, GSC, 25 Mar 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
#include "tlm_str.fin"
#include "navqc_s.fin"
      type(tlm_struct) :: tlm
      type(navqc_struct) :: navqc
c
      integer*4 gaclac, iret, attangfl(maxlin)
      real*4 att_ang(3,maxlin)
c
      real*4 toldif(3)
c
      integer*4 i1, i2, nper, ilin, nlines, i
      logical found, end
      real*4 tolmult, diff(3)

c
c
c       set up some controls
c
      nper = 1
      tolmult = 1.
      if( gaclac .eq. 1 ) then
         nper = 5       ! # actual lines per tlm line
         tolmult = 4.   ! second tolerence multiplier
      end if
c
      nlines = nper * tlm%ntlm
c
c       flag any unflagged attitude with pitch, roll, or 
c       yaw out of tolerence
c
      do ilin = 1, tlm%ntlm
        if( tlm%sc_att%flag(ilin) .eq. 0 ) then
          if( ( tlm%sc_att%att(1,ilin) .lt. navqc%sc_att(1,1) ) .or.
     1        ( tlm%sc_att%att(1,ilin) .gt. navqc%sc_att(2,1) ) .or.
     1        ( tlm%sc_att%att(2,ilin) .lt. navqc%sc_att(1,2) ) .or.
     1        ( tlm%sc_att%att(2,ilin) .gt. navqc%sc_att(2,2) ) .or.
     1        ( tlm%sc_att%att(3,ilin) .lt. navqc%sc_att(1,3) ) .or.
     1        ( tlm%sc_att%att(3,ilin) .gt. navqc%sc_att(2,3) )     )
     1          tlm%sc_att%flag(ilin) = 1
        end if
      end do
c
c       check the consistency of angles
c
      i2 = 0          ! i1, i2 are pointers to consecutive good values
      end = .FALSE.   ! true if the last unflagged value was found
      found = .FALSE. ! to signal that a consistent pair was found
c
c       start out by finding the next unflagged S/C attitude value
c
      call fndflg(tlm%sc_att%flag, tlm%ntlm, 1, i1 )
      if( i1 .le. 0 ) then
c
c          no unflagged values found at all, return with error
c
         iret = -1
         write( 6, 100 ) 
  100    format(' SC_ATT: no unflagged S/C based attitude values found')
         go to 990
      end if
c
c       place the next unflagged location in i2 and compare values
c
      do while( .not. end )
         call fndflg(tlm%sc_att%flag, tlm%ntlm, (i1 + 1), i2 )
         if( i2 .le. 0 ) then
            end = .TRUE.
         else
            found = .TRUE.  
c
c             do the actual consistency checks
c
            do i = 1,3
               diff(i) = tlm%sc_att%att(i,i2) - tlm%sc_att%att(i,i1)
               toldif(i) = navqc%att_del(i) * ( i2 - i1 ) * 
     1            nper * tolmult
            end do
c
            if( ( diff(1) .gt. toldif(1) ) .or.
     1          ( diff(2) .gt. toldif(2) ) .or.
     1          ( diff(3) .gt. toldif(3) )      ) then
               tlm%sc_att%flag(i1) = 1
               tlm%sc_att%flag(i2) = 1
            end if
c
c             for next pair, move secind pointer to the first
            i1 = i2
         end if
      end do
c
c       make sure a consistent pair was found
c
      if( .not. found ) then
         iret = -1
         write( 6, 500 )
  500    format(' SC_ATT: no consistent pairs of S/C based attitude were 
     1       found')
         go to 990
      end if
c
c       fit the attitude values to the line-by-line array
c
      call fitrng(tlm%sc_att%att, tlm%ntlm, 3, tlm%sc_att%flag, 
     1     nper, att_ang, attangfl)
c
c       and end
c
  990 continue
      return
      end
