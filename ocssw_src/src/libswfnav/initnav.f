        subroutine initnav(input,nframes,navctl,navqc,orbit,ierr)
c
c  initnav(input,nframes,navctl,navqc,orbit)
c
c  Purpose:  Initilizes navigation constants, reads control and QC parameters
c    into structures and processes GPS orbit data
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  input        struct   I      Structure for input data 
c  nframes       I*4     I      Number of input minor frames
c  navctl       struct   O      Structure for navigation control parameters
c  navqc        struct   O      Structure for navigation QC parameters
c  orbit        struct   O      Structure for processed orbit data
c
c  By: Frederick S. Patt, GSC, 5 July 93
c
c  Notes: 
c
c  Modification History:
c
c  Added ierr return parameter. B. A. Franz, GSC, 24 December 1997.
c

        implicit none
c
#include "nav_cnst.fin"
#include "navqc_s.fin"
#include "navctl_s.fin"
#include "input_s.fin"
#include "orbit_s.fin"
c
c
      type(navqc_struct) :: navqc
      type(navctl_struct) :: navctl
      type(input_struct) :: input(maxlin)
      type(orbit_struct) :: orbit
      integer*4 nframes, ierr

c  Initilize navigation constants
      call cdata

c  Read control parameter file
      call readctl(navctl, ierr)
      if (ierr .ne. 0) return

c  Read QC parameter file
      call readqc(navqc, ierr)
      if (ierr .ne. 0) return

c  Process orbit data
      call orbcomp(input, nframes, orbit, ierr)
      if (ierr .ne. 0) return

      return
      end
