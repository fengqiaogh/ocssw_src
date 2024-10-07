      subroutine Driver_MOD_PR06OD( )

!-------------------------------------------------------------------
!
!!F77
!
!!Description:
! Driver program for MODIS cloud property algorithm.
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
!
!!Revision History:
! Revision 3.0  2004/08/24 mag
! the new driver
!
!!Team-unique Header:
! Cloud Retrieval Group, NASA/GSFC, Code 913, Greenbelt, Maryland, USA
!
!!Credit and Reference:
! Programmed by:
! Mark Gray (mag)
! L-3 GSI
! NASA/GSFC, code 913
! Greenbelt, Maryland, USA
!
!  W. Robinson, SAIC -re-worked to remove work file
!
!!End
!---------------------------------------------------------------------
  use names
  use mod06_run_settings
  use ch_xfr, only: cm_from_l2, c2_g_year, c2_g_day, c2_sensor_id
  use cld_tbl_names ! virtually use them all

  implicit none
 
  !.....Parameter Declarations


  !*****************************************************
  !     Variable Declarations
  !*****************************************************

  real :: statistics(7)
  integer :: ErrorFlag
  integer :: i

  integer :: status

  !
  !
  !--------------------------------------------------------------------------
  !   now, this mainly sets the table names wheras the rest of the information 
  !   comes through l2gen
  !---------------------------------------------------------------------------
  MYTIME = 0
  ! WDR we have made time the input so...
  MYYEAR = c2_g_year
  MYDAY = c2_g_day
  !  derive the month - it is used - just get a rough value
  MYMONTH = c2_g_day /30 + 1
  if ( MYMONTH > 12 ) MYMONTH = 12
  !
  !  WDR set the cloud mask to be from the l2
  cm_from_l2 = 1
	
  ! water and ice (over land) no-wind tables
  call get_cld_tbl( c2_sensor_id, SCAT_ICE, Aice_library, status )

  call get_cld_tbl( c2_sensor_id, SCAT_WATER, Awater_library, status )
	
  ! ice wind-dependent libs
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_03_M_S, Alibnames_ice(1), status )
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_SD_03_M_S, &
    Alibnames_ice_sdev(1), status )
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_07_M_S, Alibnames_ice(2), status )
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_SD_07_M_S, &
    Alibnames_ice_sdev(2), status )
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_15_M_S, Alibnames_ice(3), status )
  call get_cld_tbl( c2_sensor_id, SCAT_ICE_SD_15_M_S, &
    Alibnames_ice_sdev(3), status )

  ! water wind-dependent libs
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_03_M_S, Alibnames_water(1), &
    status )
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_SD_03_M_S, &
   Alibnames_water_sdev(1), status )
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_07_M_S, Alibnames_water(2), &
   status )
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_SD_07_M_S, &
   Alibnames_water_sdev(2), status )
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_15_M_S, Alibnames_water(3), &
    status )
  call get_cld_tbl( c2_sensor_id, SCAT_WATER_SD_15_M_S, &
    Alibnames_water_sdev(3), status )

  ! single scattering / phase functions
  call get_cld_tbl( c2_sensor_id, SINGLE_SCAT, Aphase_library, status )
  
  ! transmittance coefficients
  call get_cld_tbl( c2_sensor_id, TRANSMITTANCE_COEFFS, &
    Atransmittance_library, status )
	
  ! sfc, snow albedo map
  call get_cld_tbl( c2_sensor_id, ECOLOGY_MAP, Aecosystem_data_name, status )
  call get_cld_tbl( c2_sensor_id, SNOW_ALBEDO, Asnowicealbedo_data_name, &
    status )
	
  ! emissivity map
  call get_cld_tbl( c2_sensor_id, GLOBAL_EMISS, Aemissivity_name, status )

  !...........Perform cloud optical properties retrieval
  call MOD_PR06OD(statistics,ErrorFlag)



End subroutine Driver_MOD_PR06OD
