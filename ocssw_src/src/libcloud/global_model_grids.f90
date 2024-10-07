module global_model_grids
  
   implicit none
  
  real ::  ice_grid( 0:719, 0:359, 1)
  real :: ozn_grid(0:359, 0:180, 1)


  logical  :: grids_are_read
  logical :: icec_is_read
  logical :: NISE_is_read
  logical :: snow_stats_are_read = .false.
  
  type ancillary_type

    real, dimension(:), allocatable :: mixr_profile
	real, dimension(:), allocatable :: temp_profile 
	real, dimension(:), allocatable :: height_profile
	real, dimension(:), allocatable :: pressure_profile
	real, dimension(:), allocatable :: o3_profile
  
    
	real :: wind_speed
	real :: Ts
	real :: Ps
	real :: col_o3
		
	real :: seaice_fraction
	real :: snow_fraction
	
	integer*1 :: surface_level
    integer*1 :: trop_level
	integer*1 :: LSM
	
  end type ancillary_type
  
contains

	subroutine get_model_idx(lat, lon, i, j) 


		real, intent(in) :: lat, lon
		integer, intent(inout ) :: i,j

		real :: x,y, x0, dx, y0, dy
	
      	x = min( max( lon,  -179.99 ), 179.99 )
      	if( x > -999. .and. x < 0.0 ) x = lon+ 360.0
      	x0 = 0.0
      	dx = 1.0
      	i = int( ( x - x0 + 0.5*dx ) / dx ) 

      	if( i .eq. 360 ) i = 0

      	y = min( max( lat, -89.99 ), 89.99 )
      	y0 = 90.0
      	dy = -1.0
      	j = int( ( y - y0 + 0.5*dy ) / dy )

	  	i = i+1
	  	j = j+1


	end subroutine get_model_idx
  

	subroutine get_model_idx_geos5(grid_xstart, grid_ystart, lat, lon, model_i, model_j)

		integer, intent(in) :: grid_xstart, grid_ystart
		real, intent(in) :: lat, lon
		integer, intent(inout) :: model_i, model_j

		real :: minLat, maxLat, minLon, maxLon
		real, parameter :: dlon = 0.625 ! 5/8
		real, parameter :: dlat = 0.5

		model_i = int ((lon + 180) / dlon ) - grid_xstart 
		model_j = int ((lat + 90) / dlat + 1 ) - grid_ystart 
		
		if (model_i < 1) model_i = 1
		if (model_j < 1) model_j = 1
		
	end subroutine get_model_idx_geos5

	subroutine find_geos5_bounds(istart, iend, jstart, jend, lat, lon)

		integer, intent(inout) :: istart, iend, jstart, jend
		real, dimension(:,:), intent(in) :: lat, lon

		real :: minLat, maxLat, minLon, maxLon
		real, parameter :: dlon = 5./8.
		real, parameter :: dlat = 0.5


		minLat  = minval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
		maxLat  = maxval( Lat,  mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )

		minLon = minval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
		maxLon = maxval( Lon, mask = ( (Lon >= -180.000 .and. Lon <= 180.0000) ) )
				
		istart = int ((minlon + 180) / dlon ) 
		if (istart < 0) istart = 0
		iend = int ((maxlon + 180) / dlon ) 
		if (iend > (360/dlon-1)) iend = 360/dlon - 1
		
		jstart = int ((minlat + 90) / dlat )  
		jend = int ((maxlat + 90) / dlat ) 
				
				
	end subroutine find_geos5_bounds

	real function get_W(RH, T, P)

! After Rogers and Yau. 

		real*8, intent(in) :: RH, T, P

		real, parameter :: c = 4187.	
		real, parameter :: cpv = 1870.
		real, parameter :: L0 = 2.501e6
		real, parameter :: T0 = 273.
	
		real :: L	
		real :: esat
		real :: ws
		real :: W
	
		if (T < 235.) then 
			L = 2.83e6   ! latent heat of ice, doesn't change very much
		else
			L = L0 - (c - cpv) * (T - T0) ! latent heat of water
		endif
	
		esat = 6.11657*exp((L/461.51)*(1/273. - 1/T) ) ! saturation vapor pressure in mb

		ws = esat * 0.622 / (P - esat) ! saturation mixing ratio in kg/kg	

		W = RH / 100. * ws ! mixing ratio in kg/kg
	
		get_W = W * 1000. ! mixing ratio in g/kg
	


	end function get_W


  
end module global_model_grids
