        subroutine readqc(navqc,ierr)
c
c  readqc(navqc,ierr)
c
c  Purpose: read in the navigation quality control parameters
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  navqc        struct   O      navigation quality control structure
c  ierr         I*4      O      return code, 0 is good
c
c  By: W. Robinson, GSC, 22 Mar 93
c
c  Notes:  
c
c  Modification History:
c
c  Changed environmental variable from NAVCTL to NAVQC.  F.S. Patt, 
c  May 3, 1994.

      implicit none
c
#include "navqc_s.fin"
c
      type(navqc_struct) :: navqc
      integer*4 iret, istat, i, j, ierr
      character*80 line,filnm
c
c
c       start, open the file of nav qc parameters
c
      ierr = -1
      filnm = '$NAVQC/navqc.dat'
      call filenv(filnm,filnm)
      open(file=filnm,unit=7,status='old',iostat=istat)
      if( istat .ne. 0 )then
         ierr = -1
         write(6,100)
  100    format(' READQC: unable to open the file navqc.dat')
      else
c
c          read in sun sensor qa values
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_tol_1(1)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_tol_1(2)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_tol_2(1)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_tol_2(2)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_del_1
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sun_del_2
c
c          earth sensor QA info
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_tol_wd(1)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_tol_wd(2)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_tol_ph(1)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_tol_ph(2)
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_del_wd
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%ear_del_ph
c
c          magnetometer (future)
c
c          time tag qa
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navqc%yearmin
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navqc%yearmax
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sectol1
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navqc%sectol2
c
c          S/C provided attitude qa
c
         do i = 1,3
            do j = 1,2
               call inline( 7, '!', line, iret )
               if( iret .ne. 0 ) go to 990
               read(line,200,err=990)navqc%sc_att(j,i)
            end do
         end do
c
         do i = 1,3
            call inline( 7, '!', line, iret )
            if( iret .ne. 0 ) go to 990
            read(line,200,err=990)navqc%att_del(i)
         end do
c
c       close the file
c
         ierr = 0
  990    continue
         close(unit=7)
c
         if( ierr .ne. 0 ) write(6,400)
  400    format(' READQC: not all nav parameters provided in the QC ',
     1     'file')
      end if
c
c       and end
c
      return
  200 format(f12.4)
  300 format(i12)
      end
