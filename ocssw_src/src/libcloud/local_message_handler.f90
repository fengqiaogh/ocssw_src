 subroutine local_message_handler(message, severity,functionname)
   
   use  nonscience_parameters
   implicit none

!   include 'PGS_SMF.f'
!   include 'PGS_IO.f'
!   include 'PGS_MODIS_39500.f'

   character*(*), intent(in) :: message, functionname
   integer,       intent(in) :: severity
   select case(severity)
      case(success)  
!        do nothing for operational algorithm
      case(warning)      
		  print*, "WARNING: ", message
!         call modis_smf_setdynamicmsg(MODIS_W_GENERIC,message,functionname)
      case(error)
		  print*, "ERROR: ", message
		  stop
!         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,functionname)
      case(failure)
		  print*, "FAILURE, CRITICAL ERROR: ", message
		  stop
!         call modis_smf_setdynamicmsg(MODIS_F_GENERIC,message,functionname)
   end select
 end subroutine local_message_handler
