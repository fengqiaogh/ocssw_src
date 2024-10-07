        subroutine readctl(navctl,ierr)
c
c  readctl(navctl,ierr)
c
c  Purpose: read in the navigation processing control parameters
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  navctl       struct   O      navigation processing control structure
c  ierr         I*4      O      return code, 0 is good
c
c  By: W. Robinson, GSC, 22 Mar 93
c
c  Notes:  
c
c  Modification History:
c
      implicit none
c
#include "navctl_s.fin"
c
      type(navctl_struct) :: navctl
      integer*4 iret, istat, i, j, k, ierr
      character*80 line,filnm
c
c
c       start, open the file of nav control parameters
c
      ierr = -1
      filnm = '$NAVCTL/navctl.dat'
      call filenv(filnm,filnm)
      open(file=filnm,unit=7,status='old',iostat=istat)
      if( istat .ne. 0 )then
         ierr = -1
         write(6,100)
  100    format(' READCTL: unable to open the file navctl.dat')
      else
c
c          read in attitude processing flag
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%procatt
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%redoyaw
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%yawtol
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%nefpts
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%neskip
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%nsfpts
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%nsskip
c
         do i = 1,3
            do j = 1,3
               call inline( 7, '!', line, iret )
               if( iret .ne. 0 ) go to 990
               read(line,200,err=990)navctl%msenoff(j,i)
            end do
         end do
c
         do i = 1,3
            call inline( 7, '!', line, iret )
            if( iret .ne. 0 ) go to 990
            read(line,200,err=990)navctl%tiltcos(i)
         end do
c
         do i = 1,3
            call inline( 7, '!', line, iret )
            if( iret .ne. 0 ) go to 990
            read(line,200,err=990)navctl%tiltcos2(i)
         end do
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%tiltfor
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%tiltaft
c
         do i = 1,3
            do j = 1,3
               do k = 1,3
                  call inline( 7, '!', line, iret )
                  if( iret .ne. 0 ) go to 990
                  read(line,200,err=990)navctl%sun_mat(k,j,i)
               end do
            end do
         end do
c
         do i = 1,3
            do j = 1,2
               call inline( 7, '!', line, iret )
               if( iret .ne. 0 ) go to 990
               read(line,200,err=990)navctl%sun_scal(j,i)
            end do
         end do
c
         do i = 1,3
            do j = 1,2
               call inline( 7, '!', line, iret )
               if( iret .ne. 0 ) go to 990
               read(line,200,err=990)navctl%sun_bias(j,i)
            end do
         end do
c
         do i = 1,2
            do j = 1,3
               do k = 1,3
                  call inline( 7, '!', line, iret )
                  if( iret .ne. 0 ) go to 990
                  read(line,200,err=990)navctl%ear_mat(k,j,i)
               end do
            end do
         end do
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%ear1sca
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%ear2sca
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%e1biasic
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%e2biasic
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%e1biasoc
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,200,err=990)navctl%e2biasoc
c
         call inline( 7, '!', line, iret )
         if( iret .ne. 0 ) go to 990
         read(line,300,err=990)navctl%lvdbug
c
c          close the file
c
         ierr = 0
  990    continue
         close(unit=7)
c
         if( ierr .ne. 0 ) write(6,400)
  400    format(' READCTL: not all nav parameters provided.')
      end if
c
c       and end
c
      return
  200 format(f12.4)
  300 format(i12)
      end
