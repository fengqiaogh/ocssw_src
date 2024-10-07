module specific_ancillary

   implicit none

   INTEGER :: NumberOfShortWavelengths  ! to change for OCI
   INTEGER, PARAMETER :: NumberOfLongWavelengths  = 2  
   INTEGER, DIMENSION(8) :: bandIndexMapSW=(/ 1,2,3,4,5,6,7,8 /)
   INTEGER, DIMENSION(NumberOfLongWavelengths)  :: bandIndexMapLW=(/ 7,8 /)

   integer, parameter :: trans_nband = 7
   integer, parameter :: index_2way = trans_nband
   
! Ackerman 1971
   real, parameter ::    ozone_absorp_coeff = 5.56209e-5  !   2.07e-21 * 2.687e16
   real, parameter ::    o3_coef_470 = 4.06e-22 * 2.687e16
   real, parameter ::    o3_coef_550 = 3.36e-21 * 2.687e16


   logical :: do_37, have_fascode

   logical :: have_band(9)


contains

   subroutine setup_atm_corr
!  WDR it is of interest that have_band is only used in specific_ancillary.f90
!  These are not in the radiance measurement order, so...?  
   use ch_xfr, only : c2_cmp_there, c2_sensor_id, OCI_ID, OCIS_ID
   use mod06_run_settings, only : band_0370

   do_37 = .true.
   if( c2_cmp_there(band_0370) == 0 ) do_37 = .false.
   if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) &
     do_37 = .false.
   have_fascode = .true.
   !  make have_band include the 2.2 um band
   !             470,     550,    650,    0.86,   0.94,   1.2,  1.6,  2.1, 2.2
   have_band = (/.true., .true., .true., .true., .true., .true., .true., &
     .true., .true. /)
   if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) ) then
     NumberOfShortWavelengths = 8
   else
     have_band(9) = .false.
     NumberOfShortWavelengths = 7
   endif

   end subroutine setup_atm_corr
!   
!  *** remap_bands()
!
   subroutine remap_bands(tau2way, temp_trans, sdev2way, temp_sdev)
   use mod06_run_settings
   use ch_xfr, only : c2_sensor_id, OCI_ID, OCIS_ID
      
   real, dimension(:), intent(in) :: temp_trans, temp_sdev
   real, dimension(:), intent(inout) :: tau2way, sdev2way
   
   !  for the non OCI branch
   if ( (c2_sensor_id /= OCI_ID) .and. (c2_sensor_id /= OCIS_ID) ) then
     tau2way(band_0065) = temp_trans(1)
     tau2way(band_0086) = temp_trans(2)
     tau2way(band_0124) = temp_trans(4)
     tau2way(band_0163) = temp_trans(5)
     tau2way(band_0213) = temp_trans(6)
     tau2way(band_0935) = temp_trans(3)
     tau2way(band_0370) = temp_trans(7)

     if (temp_sdev(1) /= -10000.) then 
       sdev2way(band_0065) = temp_sdev(1)
       sdev2way(band_0086) = temp_sdev(2)
       sdev2way(band_0124) = temp_sdev(4)
       sdev2way(band_0163) = temp_sdev(5)
       sdev2way(band_0213) = temp_sdev(6)
       sdev2way(band_0935) = temp_sdev(3)
       sdev2way(band_0370) = temp_sdev(7)
     endif
   else  ! for the OCI and OCIS
     tau2way(band_0065) = temp_trans(1)
     tau2way(band_0086) = temp_trans(2)
     tau2way(band_0124) = temp_trans(4)
     tau2way(band_0163) = temp_trans(6)
     tau2way(band_0213) = temp_trans(7)
     tau2way(band_0935) = temp_trans(3)
     tau2way(band_0226) = temp_trans(8)

     if (temp_sdev(1) /= -10000.) then
       sdev2way(band_0065) = temp_sdev(1)
       sdev2way(band_0086) = temp_sdev(2)
       sdev2way(band_0124) = temp_sdev(4)
       sdev2way(band_0163) = temp_sdev(6)
       sdev2way(band_0213) = temp_sdev(7)
       sdev2way(band_0935) = temp_sdev(3)
       sdev2way(band_0226) = temp_sdev(8)
     endif
   endif
   
   end subroutine remap_bands


! This subroutine is basically a stub as far as MODIS is concerned
! subroutine get_specific_ancillary(   mod06input_filedata,   &
!                       pressure_name, &
!                       temperature_name, &
!                       ctm_name, &
!                       sfc_temp_name, &
!                       cpi_name, &
!                                mod06_start,       &
!                                mod06_stride,      &
!                                mod06_edge, EXECUTE_STANDARD, REPLACE_ALBEDO)
!                  
!   integer,     intent(inout),dimension(:)  :: mod06_start, mod06_stride, mod06_edge, mod06input_filedata
!  logical, intent(inout) :: EXECUTE_STANDARD, REPLACE_ALBEDO
! 
!  character(*), intent(inout) :: pressure_name, temperature_name, ctm_name, sfc_temp_name, cpi_name
! 
!#ifdef CT_1KM
!
!  pressure_name = "cloud_top_pressure_1km"
!  temperature_name = "cloud_top_temperature_1km"
!  ctm_name = "cloud_top_method_1km"
!  sfc_temp_name = "surface_temperature_1km"
!  cpi_name = "Cloud_Phase_Infrared_1km"
!
!    mod06_start  = mod06_start
!    mod06_edge   = mod06_edge
!
!#else
!
!  pressure_name = "Cloud_Top_Pressure"
!  temperature_name = "Cloud_Top_Temperature"
!  ctm_name = "Cloud_Height_Method"
!  sfc_temp_name = "Surface_Temperature"
!  cpi_name = "Cloud_Phase_Infrared"
!
!    mod06_start  = floor(real(mod06_start)/5.)
!    mod06_edge   = floor(real(mod06_edge)/5.)
!
!
!#endif
!
!    mod06_stride = 1
!
!  EXECUTE_STANDARD = .true.
!  REPLACE_ALBEDO = .false.
!
! 
! end subroutine get_specific_ancillary


! subroutine get_cloud_top_properties(mapi_filedata, &
!                          pressure_sds, &
!                          temperature_sds, &
!                          ctm_sds, &
!                          sfc_temp_sds, &
!                          cpi_sds, &
!                                    start,    &
!                                    stride,   &
!                                    edge,     &
!                          start_1km, &
!                          edge_1km, &
!                                    cloud_top_pressure, &
!                                    cloud_top_temperature, &
!                          cloud_height_method, &
!                          surface_temperature, &
!                          cloud_phase_infrared, &
!                          EXECUTE_STANDARD, &
!                                    status)
!
!  use GeneralAuxType
!   use nonscience_parameters
!   use mod06_run_settings
!  use general_array_io,only: read_int_array, read_byte_array
!  implicit none
!
!  character(*), intent(in)                  :: pressure_sds, temperature_sds, ctm_sds, sfc_temp_sds, cpi_sds
!  integer,      intent(in)                  :: mapi_filedata(:)
!  integer,      intent(inout), dimension(:) :: start, stride, edge, start_1km, edge_1km
!  real(single), dimension(:,:), intent(out) :: cloud_top_pressure, cloud_top_temperature, surface_temperature
!  integer*1, dimension(:,:), intent(out ) :: cloud_height_method, cloud_phase_infrared
!  integer,      intent(out)                 :: status
!  logical, intent(in) :: EXECUTE_STANDARD
!
!  integer(integer_twobyte)                  :: fill_value
!  integer                                   :: i, j, local_i, local_j, dimsizes(2)
!  real(single), dimension(:,:), allocatable :: ctp_temp_array, ctt_temp_array, sfc_temp_array
!  integer*1, dimension(:,:), allocatable :: ctm_temp_array, cpi_temp_array
!  logical :: offset
!  integer :: localstride(2)
!  
!
!  cloud_top_pressure = 0.
!  cloud_top_temperature = 0.
!  cloud_height_method = 0
!  cloud_phase_infrared = 0.
!
!  offset = .true.
!
!#ifdef CT_1KM
!
!  call read_int_array(mapi_filedata, pressure_sds, start, stride, edge, offset, cloud_top_pressure, status)
!  call read_int_array(mapi_filedata, temperature_sds, start, stride, edge, offset, cloud_top_temperature, status)
!  call read_int_array(mapi_filedata, sfc_temp_sds, start, stride, edge, offset, surface_temperature, status)   
!  call read_byte_array(mapi_filedata, ctm_sds, start, stride, edge, cloud_height_method, status)
!  call read_byte_array(mapi_filedata, cpi_sds, start, stride, edge, cloud_phase_infrared, status)
!
!
!#else
!
!  allocate(ctp_temp_array(edge(1), edge(2)))
!  allocate(ctt_temp_array(edge(1), edge(2)))
!  allocate(ctm_temp_array(edge(1), edge(2)))
!  allocate(sfc_temp_array(edge(1), edge(2)))
!  allocate(cpi_temp_array(edge(1), edge(2)))
!
!  call read_int_array(mapi_filedata, pressure_sds, start, stride, edge, offset, ctp_temp_array, status)
!  call read_int_array(mapi_filedata, temperature_sds, start, stride, edge, offset, ctt_temp_array, status)
!  call read_int_array(mapi_filedata, sfc_temp_sds, start, stride, edge, offset, sfc_temp_array, status)
!  call read_byte_array(mapi_filedata, ctm_sds, start, stride, edge, ctm_temp_array, status)
!  call read_byte_array(mapi_filedata, cpi_sds, start, stride, edge, cpi_temp_array, status)
!
!  dimsizes(1) = size(cloud_top_pressure, 1)
!  dimsizes(2) = size(cloud_top_pressure, 2)
!
!  if (status /= success) then
!    call local_message_handler('Problem reported, see earlier message/s',status,'get_cloud_top_properties') 
!  endif 
!
!
!  do i = 1, dimsizes(1)
!    local_i = (i - 1) / 5 + 1
!    if (local_i > edge(1) ) local_i = edge(1)
!    do j = 1, dimsizes(2)
!       local_j = (j - 1) / 5 + 1
!     if (local_j > edge(2)) local_j = edge(2)
!       cloud_top_pressure(i,j) =  ctp_temp_array(local_i, local_j)
!       cloud_top_temperature(i,j) = ctt_temp_array(local_i, local_j)
!       cloud_height_method(i,j) = ctm_temp_array(local_i, local_j)
!       surface_temperature(i,j) = sfc_temp_array(local_i, local_j)
!       cloud_phase_infrared(i,j) = cpi_temp_array(local_i, local_j)
!           
!    enddo
!  enddo
!  deallocate( ctt_temp_array, ctp_temp_array, ctm_temp_array, sfc_temp_array, cpi_temp_array )
!
!#endif
!
!
! end subroutine get_cloud_top_properties


end module specific_ancillary
