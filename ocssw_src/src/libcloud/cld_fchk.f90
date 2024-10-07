SUBROUTINE cld_fchk( status, file, line )
!
!  to check the status returns from the nf_ netcdf calls
!
integer, intent (in) :: status, line
character(*), intent (in) :: file
  if( status .ne. 0 ) then
    print*, "Code failure, file: ", file, ' Line: ', line, ' Status: ', status
    stop 10
  endif
return
end
