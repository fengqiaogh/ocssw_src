        subroutine fndflg(flag,nflag,istndx,next)
c
c  fndflg(flag,nflag,istndx,next)
c
c  Purpose: find girst good flag in a flag array
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  flag         I*4      I      Flag array (0 - good, 1 - bad)
c  nflag        I*4      I      Number of valid flag entries
c  istndx       I*4      I      index to start checking on
c  next         I*4      O      Index of first good flag, -1 if
c                                none found before end
c
c  By: W. Robinson, GSC, 18 Mar 93
c
c  Modification History:
c
      implicit none
      integer*4 flag(1), nflag, istndx, next, iptr
c
c
      next = -1
c
c       just start at istndx and find next good flag
c
      if( istndx .le. nflag ) then
c
        do iptr = istndx, nflag
          if( flag(iptr) .eq. 0 ) then
            next = iptr
            go to 990
          end if
        end do
c
      end if
c
  990 return
      end
