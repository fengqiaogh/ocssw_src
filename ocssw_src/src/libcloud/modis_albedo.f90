module modis_albedo

! HI
!****************************************************************************
  ! !F90
  !
  ! !Description:
  !    This module contains the subroutines used to provide albedo and ecosystem
  !     maps for processing a MOD06 granule.
  !    There is only one callable routine, getAlbedoEco, which returns albedo 
  !     values and IGBP ecosystem classification for each pixel of the 
  !     MOD06 granule,
  !     and for the wavelengths specified.  Also returned are values of 
  !     snow albedos
  !     by ecosystem class for the wavelengths specified.
  !  
  ! !Callable routines:
  !    getAlbedoEco()
  ! 
  ! !Revision History:
  !
  ! Revision 1.0  2003/12/18  12:43:43  EGMoody
  ! Initial revision.
  !
  ! !Team-Unique Header:
  !   Cloud Retrieval Group, NASA Goddard Space Flight Center
  !
  ! !References and Credits:
  !   Written by
  !    Eric Moody
  !    Climate and Radiation Branch, Code 913
  !    NASA/GSFC
  !    Greenbelt MD 20771
  !    Eric.Moody@gsfc.nasa.gov
  !
  ! !Design Notes:
  !
  ! !END
  !*****************************************************************************


  !Dependencies:
  use GeneralAuxType
  use nonscience_parameters
  use core_arrays, only : cloudmask
  use MOD06AlbedoEcoModule
  use mod06_run_settings, only : set_albedo_bands
  implicit none

  include 'netcdf.inc'

  private

  public :: getAlbedoEco, init_snow_stats, lat_start, lat_end

  integer, parameter :: MAX_VAR_DIMS = 10
  integer :: lat_start, lat_end

  !local variables:
  !counters:
  integer  (kind = integer_fourbyte)  :: i,j,k,l,m,n

  !HDF variables:
  integer  (kind = integer_fourbyte)  :: HDFstatus
  integer  (kind = integer_fourbyte  ),  &
             dimension(MAX_VAR_DIMS)     :: hdfStart, hdfStride, hdfEdge
  integer  (kind = integer_fourbyte)     :: sds_id, sds_index

  integer  (kind = integer_fourbyte), parameter  :: NumSnowAlbStatWaves = 10
  integer  (kind = integer_fourbyte), parameter  :: NumEcosystems = 18
  integer  (kind = integer_fourbyte), parameter  ::  NumSnowTypes = 2

  real :: SnowAlbedoStats (1:NumSnowAlbStatWaves, 0:NumEcosystems-1, &
          1:NumSnowTypes)
  real :: ProcessSnowAlbStats (1:set_albedo_bands,  0:NumEcosystems-1, &
          1:NumSnowTypes)

  integer  (kind = integer_fourbyte), parameter  :: NumSeaIceAlbWaves = 8
  real     (kind = single), dimension(NumSeaIceAlbWaves) :: meltSeaIceAlbedos, &
            drySeaIceAlbedos
  real (kind = single), dimension(set_albedo_bands) :: &
    ProcessMeltSeaIceAlbedos, ProcessDrySeaIceAlbedos,  TransitionAlbedo, & 
    Slope, Intercept

  real   (kind = single), parameter :: PercentageOfDry = 0.800d0
  real :: Wavelength(10)

contains


  subroutine init_snow_stats(SnowAlbedoFN, WavelengthText) 
  use ch_xfr, only: c2_sensor_id, OCI_ID, OCIS_ID
  character (len = *), intent(in ) :: SnowAlbedoFN
  character (len = *),  dimension(:), intent(in ) :: WavelengthText
		
  real (kind = single), dimension(NumSnowAlbStatWaves) :: PossStatWavelengths
  real (kind = single), dimension(NumSeaIceAlbWaves) :: PossSeaIceWavelengths

  integer :: errorLevel, i, NumAlbWavesProcess
  logical :: FoundWave

  NumAlbWavesProcess = set_albedo_bands

  !*****************************************************************************
  ! Read in the snow albedo statistics from the corresponding file:
  !*****************************************************************************
  call ReadSnowAlbStats  ( SnowAlbedoFN, NumSnowTypes, NumSnowAlbStatWaves,  &
    NumEcosystems, SnowAlbedoStats, errorLevel )

  !  WDR note that errorLevel is not dealt with - we'll just kill down below
  !*****************************************************************************
  ! For each wavelength to be processed, store the albedo snow statistic in the
  !  cooresponding index in a temp array. Note that 8-10 are broad band, and 
  !  not included.
  !*****************************************************************************
  !Load in the order of the statistical wavelengths from the file:
  PossStatWavelengths(1)  = 0.470000
  PossStatWavelengths(2)  = 0.555000
  PossStatWavelengths(3)  = 0.659000
  PossStatWavelengths(4)  = 0.858000
  PossStatWavelengths(5)  = 1.240000
  PossStatWavelengths(6)  = 1.640000
  PossStatWavelengths(7)  = 2.130000
  PossStatWavelengths(8)  = 0.0 
  PossStatWavelengths(9)  = 0.0
  PossStatWavelengths(10) = 0.0
    
  !Loop over the wavelengths to be processed, storing the cooresponding
  ! data from the complete snowAlbedoStats array.  Do not process 3.7,
  ! the last index as there are no cooresponding values, instead compute
  ! from the 2.1um at the end:
    
  do i = 1, NumAlbWavesProcess
    read( WavelengthText(i), * )  Wavelength(i) 
  end do
   
  ProcessSnowAlbStats= 0.
  do i = 1, NumAlbWavesProcess-1
    FoundWave = .false.

    do j = 1, NumSnowAlbStatWaves-3
      if ( mod(Wavelength(i),PossStatWavelengths(j)) < 0.05 ) then
        ProcessSnowAlbStats(i,:,:) = SnowAlbedoStats(j,:,:)
        FoundWave = .true.
        goto 100
      end if
    end do
 100 continue
    if (.not. FoundWave) then
      return
    end if
  end do
    
  if( ( c2_sensor_id == OCI_ID ) .OR. ( c2_sensor_id == OCIS_ID ) )then   
              !  6th spot is 2.2 for OCI dupe the 2.1
    ProcessSnowAlbStats(6,:,:) = ProcessSnowAlbStats(5,:,:)
    print*,'WDR NOTE both albedo and snow stats are cloned from 2.1 to 2.2 um'
  else   ! 6th spot is 3.7 for modis
    !compute the 3.7um (assume 6 = 3.6 and 5 = 2.1um):
    ProcessSnowAlbStats(6,:,:) = ProcessSnowAlbStats(5,:,:) * 0.5000
  endif
                        
  !*****************************************************************************
  ! For each wavelength to be processed, store the albedo sea ica albedo in the
  !  cooresponding index in a temp array. 
  !*****************************************************************************
  !Load in the order of the statistical wavelengths from the file:
  PossSeaIceWavelengths(1) = 0.470000
  PossSeaIceWavelengths(2) = 0.555000
  PossSeaIceWavelengths(3) = 0.659000
  PossSeaIceWavelengths(4) = 0.858000
  PossSeaIceWavelengths(5) = 1.240000
  PossSeaIceWavelengths(6) = 1.640000
  PossSeaIceWavelengths(7) = 2.130000
  PossSeaIceWavelengths(8) = 3.700000
    
  ! Set the sea ice values to be equal to the dry perm snow/ice 
  ! values or percentage of them:
  drySeaIceAlbedos(1:7)  = SnowAlbedoStats(1:7,15,1)
  meltSeaIceAlbedos(1:7) = SnowAlbedoStats(1:7,15,1) * PercentageOfDry
  !3.7 is half of 2.1:
  drySeaIceAlbedos(8)  =  drySeaIceAlbedos(7) * 0.5000d0
  meltSeaIceAlbedos(8) = meltSeaIceAlbedos(7) * 0.5000d0
                         
  !Loop over the wavelengths to be processed, storing the cooresponding
  ! data from the complete snowAlbedoStats array.  Process 3.7um for this
  ! sea ice case
  do i = 1, NumAlbWavesProcess
    FoundWave = .false.
    do j = 1, NumSeaIceAlbWaves
       if ( mod(Wavelength(i),PossSeaIceWavelengths(j)) < 0.05 ) then
          ProcessDrySeaIceAlbedos(i)  = drySeaIceAlbedos(j)
          ProcessMeltSeaIceAlbedos(i) = meltSeaIceAlbedos(j)
          FoundWave = .true.
          goto 101
       end if
    end do
  101 continue
    if (.not. FoundWave) then
!       status = 19
!       return
    end if
  end do

  end subroutine init_snow_stats

! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
  subroutine getAlbedoEco (Lat, Long, JulianDay, EcosystemFN,  &
    emissivity_name, Debug, WavelengthText, ProcessLandAlbedo,  &
    icec, Albedo, emissivity, Ecosystem, Status, sfc_info )

  ! WDR access c2_alt_snowice and albedo
  use ch_xfr, only: c2_alt_snowice, c2_sfc_albedo, c2_sensor_id, OCI_ID, &
    OCIS_ID

  ! EM -- June 17, 2005
  ! Apparently ncdump does not handle the metadata very well, and there 
  ! was some binary conversion 
  ! issue such that the metadata values are mixed up.  If you use something 
  ! like HDFlook, they are fine.
  ! Here is how ncdump prints them out:
  !  byte IGBP_Land_Cover_Type(NumLatPoints, NumLongPoints) ;
  !    IGBP_Land_Cover_Type:long_name = "IGBP Land Cover Classification Type" ;
  !    IGBP_Land_Cover_Type:units = "class number" ;
  !    IGBP_Land_Cover_Type:valid_range = '\0', '\376' ;
  !    IGBP_Land_Cover_Type:_FillValue = '\377' ;
  !    IGBP_Land_Cover_Type:water = '\0' ;
  !    IGBP_Land_Cover_Type:evergreen needleleaf forest = '\1' ;
  !    IGBP_Land_Cover_Type:evergreen broadleaf forest = '\2' ;
  !    IGBP_Land_Cover_Type:deciduous needleleaf forest = '\3' ;
  !    IGBP_Land_Cover_Type:deciduous broadleaf forest = '\4' ;
  !    IGBP_Land_Cover_Type:mixed forests = '\5' ;
  !    IGBP_Land_Cover_Type:closed shrubland = '\6' ;
  !    IGBP_Land_Cover_Type:open shrublands = '\7' ;
  !    IGBP_Land_Cover_Type:woody savannas = '\10' ;
  !    IGBP_Land_Cover_Type:savannas = '\11' ;
  !    IGBP_Land_Cover_Type:grasslands = '\12' ;
  !    IGBP_Land_Cover_Type:permanent wetlands = '\13' ;
  !    IGBP_Land_Cover_Type:croplands = '\14' ;
  !    IGBP_Land_Cover_Type:urban and built-up = '\15' ;
  !    IGBP_Land_Cover_Type:cropland/natural vegetation mosaic = '\16' ;
  !    IGBP_Land_Cover_Type:snow and ice = '\17' ;
  !    IGBP_Land_Cover_Type:barren or sparsely vegetated = '\20' ;
  !    IGBP_Land_Cover_Type:unclassified = '\376' ;
  ! If you notice it skips 7 and goes to 10 (open shrub to woody savanna) 
  ! and then 17 to 20 (snow/ice to barren).
  !
  ! Here is how my code actually assigns them:
  !
  !*****************************************************************************
  ! Create the Ecosystem_Order array, used in HDF initialization:
  !*****************************************************************************
  !    Ecosystem_Order(0)  = 'Ocean / Water'
  !    Ecosystem_Order(1)  = 'Evergreen Needle Forest'
  !    Ecosystem_Order(2)  = 'Evergreen Broad Forest'
  !    Ecosystem_Order(3)  = 'Deciduous Needle Forest'
  !    Ecosystem_Order(4)  = 'Deciduous Broad Forest'
  !    Ecosystem_Order(5)  = 'Mixed Forest'
  !    Ecosystem_Order(6)  = 'Closed Shrubs'
  !    Ecosystem_Order(7)  = 'Open / Shrubs'
  !    Ecosystem_Order(8)  = 'Woody Savanna'
  !    Ecosystem_Order(9)  = 'Savanna'
  !    Ecosystem_Order(10) = 'Grassland'
  !    Ecosystem_Order(11) = 'Wetland'
  !    Ecosystem_Order(12) = 'Cropland'
  !    Ecosystem_Order(13) = 'Urban'
  !    Ecosystem_Order(14) = 'Crop Mosaic'
  !    Ecosystem_Order(15) = 'Antarctic / Permanant Snow'
  !    Ecosystem_Order(16) = 'Barren / Desert'
  !    Ecosystem_Order(17) = 'Tundra'   
  !
  !Ncdump is quite deficient in this respect and I get asked this all the time.
  !But it is not a coding issue. -- E. Moody

  real      (kind = single ), dimension(:,:), intent(in ) :: Lat, Long
  character (len = *), intent(in ) :: EcosystemFN
  character (len = *), intent(in ) :: emissivity_name
  integer   (kind = integer_fourbyte ),  intent(in ) :: JulianDay
  logical,  intent(in ) :: Debug
  character (len = *), dimension(:),  intent(in ) :: WavelengthText
  logical,  dimension(:), intent(in ) :: ProcessLandAlbedo
  real      (kind = single), dimension(:,:),  intent(in ) :: icec 
  real (kind = single ), dimension(:,:,:),  intent(out) :: emissivity
  integer*2, dimension(:,:,:), intent(out) :: Albedo 
  integer*1, dimension(:,:), intent(out) :: Ecosystem, sfc_info
  integer   (kind = integer_fourbyte), intent(out) :: Status

  ! !Description:
  !    This returns albedo values and IGBP ecosystem classification for each 
  !    pixel 
  !    of the MOD06 granule, and for the wavelengths specified.  Also 
  !    returned are 
  !    values of snow albedos by ecosystem class for the wavelengths specified.
  !  
  ! !Input Parameters:
  !    Lat, Long       : Latitude and Longitude arrays from granule.
  !    WavelengthFN    : The albedo filenames for each specified wavelength 
  !                      to process.
  !    EcosystemFN     : The filename for the IGBP Ecosystem Classification map.
  !    Julian Day      : Julian Day
  !    Debug           : Logical Flag for printing during debugging.
  !    WavelengthText  : Text array of wavelengths to be processed.
  !
  ! !Output Parameters:
  !    Wavelength      : Real values of the Text array.
  !    Albedo          : 3D array containing the albedo values for each 
  !                     lat/lon point
  !                      and for each wavelength specified.
  !    SnowAlbedo      : Snow albedos by ecosystem class for each 
  !                       wavelength specified.
  !    Ecosystem       : IGBP ecosystem classification for each lat/lon point.
  !    Status          : Flag specifying if this routine worked correctly.
  !
  ! !Revision History:
  !   See Module revision history at the beginning of the file.
  !
  ! !Team-Unique Header:
  !   Cloud Retrieval Group, NASA Goddard Space Flight Center
  !
  ! !References and Credits:
  !   Written by
  !    Eric Moody
  !    Climate and Radiation Branch, Code 913
  !    NASA/GSFC
  !    Greenbelt MD 20771
  !    Eric.Moody@gsfc.nasa.gov
  !
  ! !Design Notes:
  !  Albedos values for each lat/lon point are provided by finding 
  !  the nearest neighbor
  !    from the ancillary albedo maps.
  !  Wavelengths to be processed are specified through the WavelengthText 
  !  array and
  !    corresponding WavelengthFN array.
  !  There is no ancillary files for 3.7 micron, so this is assumed to be 
  !  half the 2.1 value.
  !  The 3.7 micron always has to be last in the WavelengthText array and 
  !   WavelengthFN array,
  !    however the WavelengthFN for 3.7 is not used, so it can be fill.
  !  As the albedo ancillary data is for land only, water values are 
  !  assumed to be 5% for
  !    all wavelengths.
  !  Sea ice is approximated with perm snow/ice values for a winter time, 
  !  and a percentage 
  !    during a melt season.  A transition period is approximated between 
  !  winter and melt.
  !  
  !  The following status error values are used:
  !  1 = Array Sizes not valid
  !  2 = Invalid min/max long/lat
  !  3 = Invalid global map start/stop x/y
  !  4 = Could not gain access to Albedo HDF file
  !  5 = Could not gain access to Albedo SDS
  !  6 = Could not read Albedo SDS
  !  7 = Could not end access to Albedo SDS
  !  8 = Could not end access to Albedo HDF file
  !  9 = Could not gain access to Ecosystem HDF file
  !  10 = Could not gain access to Ecosystem SDS
  !  11 = Could not read Ecosystem SDS
  !  12 = Could not end access to Ecosystem SDS
  !  13 = Could not end access to Ecosystem HDF file
  !  14 = Trouble with Albedo Snow Statistics
  !  15 = Trouble with Reading NISE data
  !  16 = Pixel's computed global map index invalid
  !  17 = Pixel's geolocation invalid
  !  18 = Wavelengths do not match for snow statistics
  !  19 = Wavelengths do not match for sea ice values
  !
  !
  ! !END

  !Local Variables

  !
  !Size variables:
  !
  !Number of Albedo Wavelengths to be processed.
  integer  (kind = integer_fourbyte) :: NumAlbWavesProcess
  !Number of Columns and Rows of the granule to be processed
  integer  (kind = integer_fourbyte) :: NumberOfCols, NumberOfRows

  !
  !Min/Max Lat/Lon:
  !
  real     (kind = single) :: minLat, maxLat, minLong, maxLong

  ! 
  ! Albedo/Ecosystem map variables:
  !
  ! Global map dims and maximum number of albedo wavelengths:
  integer  (kind = integer_fourbyte)  :: NumMapCols, NumMapRows
  integer  (kind = integer_fourbyte)  :: NumMapCols_alb, NumMapRows_alb
                                                                 
  ! Stores the resolution (in degrees) of the global map
  real     (kind = single)    :: Map_Resolution
  real     (kind = single)  :: Map_Resolution_alb
  ! Stores the upper right corner global map coords, used to compute rest 
  ! of map coords.
  real  (kind = single)   :: westernMostLongitude, northernMostLatitude
  real  (kind = single)   :: westernMostLongitude_alb, northernMostLatitude_alb
  ! Indices for section of global map that is needed for supplying values 
  ! for granule:
  integer  (kind = integer_fourbyte)  :: startX, stopX, startY, stopY, &
    NumberOfXPoints,  NumberOfYPoints  
  integer  (kind = integer_fourbyte) :: startX_alb, stopX_alb,  &
    startY_alb, stopY_alb
  real (kind = single) :: startXval, stopXval, startYval, stopYval 
  ! Stores the real and integer values of the section of the global albedo map.
  ! Used to scale the albedo values read in as integers from the HDF files.
  real  (kind = single ), parameter :: ScaleFactor = 0.00100
  ! Store the NISE snow type:
  INTEGER(integer_fourbyte), allocatable, DIMENSION(:,:)    :: snowIceType
  ! Store the section of the global ecosystem map needed for supplying 
  ! `values for granule:
  integer (kind = integer_onebyte), allocatable, dimension(:,:) :: EcosystemMap 
  ! Dummy index variables:
  integer  (kind = integer_fourbyte) :: xIndex, yIndex 
  integer  (kind = integer_fourbyte) :: xIndex_alb, yIndex_alb 
  ! HDF variables:
  ! File IDs
  integer  (kind = integer_fourbyte) :: EcoMapFID, AlbMapFID
  ! SDS IDs
  integer  (kind = integer_fourbyte) :: EcoMapSDSID, AlbMapSDSID    
  ! Stores the SDS name
  character(len =  MAX_SDS_NAME_LEN) :: SDSName

  ! Error handling:
  INTEGER(integer_fourbyte)    :: errorLevel

  ! Wavelength descriptions used in snow statistics:
  !Variables used for sea ice albedos:
  LOGICAL :: arcticMeltingSeason, WinterToMeltTransition, MeltToWinterTransition
  integer  (kind = integer_fourbyte), parameter :: TransitionDays = 10, &
    NHMeltBeg = 152, NHMeltEnd = 244, SHMeltBeg = 334, SHMeltEnd = 61
    
  !Water albedo values:
  real  (kind = single), parameter  :: WaterAlbedo = 0.050000
  real, dimension(:), allocatable :: albedo_real4
  real, parameter :: emissivity_fill = 0.985 ! according to Wisconsin
	
  ! emissivity map variables
  real, parameter :: emissivity_res = 0.05
  integer, parameter :: num_emiss_cols = 7200
  integer, parameter :: num_emiss_rows = 3600
  real, dimension(:,:,:), allocatable :: emissivity_map
  integer*2, dimension(:,:), allocatable :: emissivity_map_int
  integer :: emiss_x, emiss_y
  real :: emiss_west, emiss_north
	
  integer :: file_id, var_id, err_code, start(2), stride(2), edge(2)
	
  integer :: FAIL = -1

  ! WDR this will support the re-read only when lat, long range falls outside 
  ! the previous range of min, max (lat, lon)
  real, dimension(2) :: min_geo_sav = (/ -9000., -9000. /), &
                         max_geo_sav = (/ -9000., -9000. /)
  real :: geo_margin = 5.  ! 5 degree extra margin to allocate for

  !*****************************************************************************
  ! Determine the size of the arrays:  
  !*****************************************************************************
  !Determine dimensions:
  NumberOfCols      = size( Albedo(:,:,:),   dim = 1 )
  NumberOfRows      = size( Albedo(:,:,:),   dim = 2 )
  NumAlbWavesProcess = size( Albedo(:,:,:),   dim = 3 )  
	
  if( .not. allocated( albedo_real4) ) &
    allocate(albedo_real4(NumAlbWavesProcess))
	
    !print*,'c,R ',NumberOfCols,NumberOfRows
    !Error checking:
    if ( (NumberOfCols < 1)      .or. &
      (NumberOfRows < 1)      .or. &
      (NumAlbWavesProcess < 2) .or. &
      (NumEcosystems < 1)          ) then
      status = error
      call local_message_handler( &
        'Problem with albedo array size, contact SDST',status,'getAlbedoEco')
      return     
    end if

    !************************************************************************
    ! Determine the min/max of lat/lon, which will be used to determine what 
    ! portion
    !  of the albedo/eco maps need to be read in. 
    !************************************************************************
    ! Determine min and max:
    if (maxval(Lat) == -999. .or. maxval(Long) == -999.) then
      ! there is no gelocation data in the scan.
      Albedo(:,:,:) = -999.  
	   emissivity(:,:,:) = -999.
      return
    else
      minLat  = minval( Lat,  &
        mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
      maxLat  = maxval( Lat,  &
        mask = ( (Lat  >= -90.0000 .and. Lat  <= 90.00000) ) )
      minLong = minval( Long, &
        mask = ( (Long >= -180.000 .and. Long <= 180.0000) ) )
      maxLong = maxval( Long, &
        mask = ( (Long >= -180.000 .and. Long <= 180.0000) ) )
    endif

    ! WDR see if we need to re-read the grids for geo bounds outside current
    !PRINT*, 'Decide whether to read the albedo/emiss tables'
    if( ( minLat < min_geo_sav(1) ) .or. ( maxLat > max_geo_sav(1) ) .or. &
        ( minLong < min_geo_sav(2) ) .or. ( maxLong > max_geo_sav(2) ) ) then
      !PRINT*, 'Reading albedo/emiss tables'
      !PRINT*, 'minLat, maxLat, minLong, maxLong ', minLat, maxLat, minLong, maxLong
      !PRINT*, 'min_geo_sav, max_geo_sav ', min_geo_sav, max_geo_sav

      ! set the new margin-enhanced range and re-read from files
      minLat = minLat - geo_margin
      maxLat = maxLat + geo_margin
      minLong = minLong - geo_margin
      maxLong = maxLong + geo_margin

      if( minLat < -90. ) minLat = -90.
      if( maxLat > 90. ) maxLat = 90.
      if( minLong < -180. ) minLong = -180.
      if( maxLong > 180. ) maxLong = 180.

      min_geo_sav(1) = minLat
      min_geo_sav(2) = minLong
      max_geo_sav(1) = maxLat
      max_geo_sav(2) = maxLong
      !PRINT*, 'NEW min_geo_sav, max_geo_sav ', min_geo_sav, max_geo_sav
 
  
      !**********************************************************************
      ! Set up the global albedo/ecosystem map resolution and determine the 
      ! start/stop
      !  points for the portion of the global map to be read in.  This 
      ! is done by
      !  computing the map indices from the min/max lat/lon values.
      !**********************************************************************
      ! 9-30-11 For Collection 4 (Moody) surface albedo and ecosystem maps 
      !  were identical resolution,
      !  the resolution of the C5 SA dataset however is double C4, thus 
      ! appropriate modifications  
      !   have been made below to read the higher resolution Albedo map 
      ! data. - GTA
      
      NumMapCols = 21600
      NumMapRows = 10800
        
      Map_Resolution = 360.000d0 / float(NumMapCols)
      
      ! Compute the WesternMostLongitude and NothernMostLatitude of the maps:
      westernMostLongitude = Map_Resolution / 2.0 - 180.0
      northernMostLatitude = 90.0 - Map_Resolution / 2.0
      !print*,'mapres=',Map_Resolution,westernMostLongitude,northernMostLatitude
      !print*,'minLong=',minLong,maxLong,minLat,maxLat
     
      ! Because the albedo map res. is different than the ecosystem map, 
      ! save the albedo map
      !  parameters for later determination surface albedo for each 
      ! granule pixel  - GTA
        
      ! Compute the start/stop x/y indices based on min/max lat/lon values:
      ! use anint so that it rounds to nearest index.
      startx = anint( (minLong-westernMostLongitude) / (Map_Resolution) ) + 1
      stopx  = anint( (maxLong-westernMostLongitude) / (Map_Resolution) ) + 1
      starty = anint( (northernMostLatitude-maxLat ) / (Map_Resolution) ) + 1
      stopy  = anint( (northernMostLatitude-minLat ) / (Map_Resolution) ) + 1
      
      ! Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
      ! indices of one beyond the map bounds, so correct that here:
      if (startx == 0) startx = 1
      if (stopx  == NumMapCols+1 ) stopx = NumMapCols
      if (starty == 0) starty = 1
      if (stopy  == NumMapRows+1 ) stopy = NumMapRows
      
      ! Compute the number of points:
      NumberOfXPoints = stopx - startx + 1
      NumberOfYPoints = stopy - starty + 1
  
      ! Error Checking:
      if ( (startx < 1) .or. (stopx > NumMapCols) .or. &
        (starty < 1) .or. (stopy > NumMapRows) .or. &
        (startx > stopx) .or. (starty > stopy)      ) then
        ! status = error
        ! call local_message_handler('Problem with loop index, contact SDST', &
        !  status,'getAlbedoEco')
        ! return  
  	     if (stopx > NumMapCols) stopx = NumMapCols
        if (stopy > NumMapRows) stopy = NumMapRows
        if (startx > stopx) startx = stopx-1
        if (starty > stopy) starty = stopy-1
        if (startx < 1) startx = 1
        if (starty < 1) starty = 1
  		
      end if
  
      !***********************************************************************
      ! Allocate the global albedo and ecosytem map storage arrays based 
      ! upon these
      !  start and stop values.  Note that the arrays will be indexed using 
      ! start/stop
      !  variables, instead of 1-NumberOfZPoints, for ease of selecting 
      ! the nearest
      !  neighbor to the inputted lat/lon grid.
      !************************************************************************
      ! WDR re-allocate here
      if( allocated(EcosystemMap) ) deallocate( EcosystemMap )
      allocate( EcosystemMap (startx:stopx, starty:stopy ) )
    			    
      !*************************************************************************
      ! Set up the HDF read variables, making note that HDF arrays start at 0
      !  instead of F90's 1:
      !************************************************************************
      ! Define the starting point
      hdfStart (:) = 1
      hdfStart (1) = startX
      hdfStart (2) = startY
        
      ! Define the stride = 1 so that it doesn't skip any data.
      hdfStride(:) = 1
        
      ! Define the Number of Points to read:
      hdfEdge  (:) = 1
      hdfEdge  (1) = NumberOfXPoints
      hdfEdge  (2) = NumberOfYPoints
          
      !************************************************************************
      ! Read in the portion of the ecosystem map from the corresponding file:
      !***********************************************************************
      ! Define the SDS Name:
      SDSName = "IGBP_Land_Cover_Type"
      
      ! Open the HDF file for reading:
      status = nf_open( Trim(EcosystemFN), NF_NOWRITE, EcoMapFID )
      if (EcoMapFID == Fail) then
        status = failure
        call local_message_handler('Problem detected ecosysyem file sfstart', &
          status,'getAlbedoEco')
      end if
      call cld_fchk( status, __FILE__, __LINE__ )
        
      ! Obtain the data set ID's:
      status = nf_inq_varid( EcoMapFID, trim(SDSName), EcoMapSDSID )
      if (EcoMapSDSID == Fail) then
        status = failure
        call local_message_handler('Problem detected ecosystem file sfselect', &
          status,'getAlbedoEco')
      end if
      call cld_fchk( status, __FILE__, __LINE__ )
      
      ! Read in the data:
      HDFstatus = nf_get_vara_int1( EcoMapFID, EcoMapSDSID, hdfStart, &
        hdfEdge, EcosystemMap )
      !Error checking:
      if (HDFstatus == FAIL) then
        status = failure
        call local_message_handler('Problem detected ecosystem file sfrdata', &
           status,'getAlbedoEco')
      end if
      call cld_fchk( HDFstatus, __FILE__, __LINE__ )
          
      HDFstatus = nf_close( EcoMapFID   )
      if (HDFstatus == FAIL) then
        status = failure
        call local_message_handler('Problem detected ecossystem file sfend', &
          status,'getAlbedoEco')
      end if
      call cld_fchk( HDFstatus, __FILE__, __LINE__ )

      !************************************************************************
      ! Read in the portion of the emissivity map from the corresponding file:
      !***********************************************************************
    	    
      ! Compute the WesternMostLongitude and NothernMostLatitude of the maps:
      emiss_west = emissivity_res / 2.0 - 180.0
      emiss_north = 90.0 - emissivity_res / 2.0
        
      ! Compute the start/stop x/y indices based on min/max lat/lon values:
      !use anint so that it rounds to nearest index.
      startx = anint( (minLong-emiss_west) / (emissivity_res) ) + 1
      stopx  = anint( (maxLong-emiss_west) / (emissivity_res) ) + 1
      starty = anint( (emiss_north-maxLat ) / (emissivity_res) ) + 1
      stopy  = anint( (emiss_north-minLat ) / (emissivity_res) ) + 1
        
      ! Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
      ! indices of one beyond the map bounds, so correct that here:
      if (startx == 0) startx = 1
      if (stopx  == num_emiss_cols+1 ) stopx = num_emiss_cols
      if (starty == 0) starty = 1
      if (stopy  == num_emiss_rows+1 ) stopy = num_emiss_rows
        
      ! Compute the number of points:
      NumberOfXPoints = stopx - startx + 1
      NumberOfYPoints = stopy - starty + 1
  
      ! Error Checking:
      if ( (startx < 1) .or. (stopx > num_emiss_cols) .or. &
        (starty < 1) .or. (stopy > num_emiss_rows) .or. &
        (startx > stopx) .or. (starty > stopy)      ) then
    
    		if (stopx > NumMapCols) stopx = NumMapCols
    		if (stopy > NumMapRows) stopy = NumMapRows
    		if (startx > stopx) startx = stopx-1
    		if (starty > stopy) starty = stopy-1
    		if (startx < 1) startx = 1
    		if (starty < 1) starty = 1
    
        ! status = error
        ! call local_message_handler('Problem with loop index, contact SDST', &
        ! status,'getAlbedoEco')
        ! return  
      end if
    	
      ! Define the starting point
      hdfStart (:) = 1
      hdfStart (1) = startX
      hdfStart (2) = startY
        
      ! Define the stride = 1 so that it doesn't skip any data.
      hdfStride(:) = 1
        
      ! Define the Number of Points to read:
      hdfEdge  (:) = 1
      hdfEdge  (1) = NumberOfXPoints
      hdfEdge  (2) = NumberOfYPoints	
    	
      ! WDR re-allocate here instead
      if( allocated( emissivity_map ) ) then
        deallocate( emissivity_map )
      end if
      allocate( emissivity_map(startx:stopx, starty:stopy, 2) )
      allocate( emissivity_map_int (startx:stopx, starty:stopy ) )
    
      status = nf_open( trim(emissivity_name), NF_NOWRITE, EcoMapFID )
      call cld_fchk( status, __FILE__, __LINE__ )
    
      do i=1, 2
  	
        if (i==1) SDSName = "emiss7"
        if (i==2) SDSName = "emiss14"
    	
        status = nf_inq_varid( EcoMapFID, trim(SDSName), EcoMapSDSID )
        call cld_fchk( status, __FILE__, __LINE__ )
    		
        HDFstatus = nf_get_vara_int2( EcoMapFID, EcoMapSDSID, hdfStart, &
          hdfEdge, emissivity_map_int )
        call cld_fchk( status, __FILE__, __LINE__ )

        emissivity_map(:,:,i) = emissivity_map_int * ScaleFactor
    		
      end do
    
      status = nf_close( EcoMapFID )
      call cld_fchk( status, __FILE__, __LINE__ )
    
      deallocate(emissivity_map_int)
      ! print*,'to snowalb'
  
    end if
  
    !*****************************************************************************
    ! Read in the NISE snow type from the corresponding file:
    !*****************************************************************************
    allocate( snowIceType   (1:NumberOfCols, 1:NumberOfRows   ) ) 
    snowIceType = 0  ! WDR-UIV
  
    if (errorLevel == FAIL) then
      status = failure
      call local_message_handler('Problem reported from GetNISEType', &
        status,'getAlbedoEco')
      return
    end if
                                       
    ! print*,'tp latlon loop',NumberOfCols,NumberOfRows
  
    ! WDR we insert the l2 snow-ice array here
    snowIceType = c2_alt_snowice
    !*****************************************************************************
    ! Loop over each pixel in the lat/lon grid, compute the global coordinates and
    !  set the pixel's albedo and ecosystem values to the corresponding global
    !  albedo and ecosystem map values.  This is essentially a nearest neighbor 
    !  calculation.
    !*****************************************************************************
    do i = 1, NumberOfCols
      do j = 1, NumberOfRows
        
        sfc_info(i,j) = 0
        ! Check to see if the geolocation is valid, also check we have a cloud. 
        if ( (Lat (i,j) >= -90.0000) .and. &
          (Lat (i,j) <=  90.0000) .and. &
          (Long(i,j) >= -180.000) .and. &
          (Long(i,j) <=  180.000) ) then 
  
        ! put the emissivity values where they belong 
  		emiss_x = anint( (Long(i,j)-emiss_west) / (emissivity_res) ) + 1
  	   emiss_y = anint( (emiss_north-Lat (i,j)) / (emissivity_res) ) + 1
  
        if (emiss_x == 0) emiss_x = 1
        if (emiss_x  == num_emiss_cols+1 ) emiss_x = num_emiss_cols
        if (emiss_y == 0) emiss_y = 1
        if (emiss_y  == num_emiss_rows+1 ) emiss_y = num_emiss_rows
  
  
  	   emissivity(i,j, 1) = emissivity_map(emiss_x, emiss_y, 1)
  		emissivity(i,j, 2) = emissivity_map(emiss_x, emiss_y, 2)
  			
        if (emissivity(i,j,1) < 0.) emissivity(i,j,1) = emissivity_fill
        if (emissivity(i,j,2) < 0.) emissivity(i,j,2) = emissivity_fill
  			
        ! end emissivity additions
  
        ! print*,'rows=',j,i,Map_Resolution,Long(i,j),westernMostLongitude, &
        !  northernMostLatitude,Lat (i,j)
           
        ! Valid geolocation, so proceed.
            
        ! Compute the ecosystem map x,y indices from this pixel's lat/long values.
        xIndex = anint( (Long(i,j)-westernMostLongitude) / (Map_Resolution) ) + 1
        yIndex = anint( (northernMostLatitude-Lat (i,j)) / (Map_Resolution) ) + 1

        ! Due to roundoff, geocoordinates of -90, 90, -180, 180 can result in 
        ! indices of one beyond the map bounds, so correct that here:
        if (xIndex == 0) xIndex = 1
          if (xIndex  == NumMapCols+1 ) xIndex = NumMapCols
          if (yIndex == 0) yIndex = 1
          if (yIndex  == NumMapRows+1 ) yIndex = NumMapRows
          ! print*,'xIndex, NumMapCols,yIndex, NumMapRows=',xIndex, &
          ! NumMapCols,yIndex, NumMapRows
          ! Checks to make sure indices are within bounds of the map:
          if (      (xIndex < 1) .or. (xIndex > NumMapCols) &
            .or. (yIndex < 1) .or. (yIndex > NumMapRows) ) then
            status = error
            call local_message_handler( &
       'Problem detected, bad array indexes (main loop) ',status,'getAlbedoEco')
              return
          end if
  
          ! Implement Rich Frey (UWisc) interpolation scheme in which he 
          ! sets the map
          ! to 200 if it should be set to snow/ice.  For sea ice, he also 
          ! uses the icec data:
#ifndef GEOS5_SFC
          if (Lat(i,j) >= 40.0 .or. Lat(i,j) <= -40.0) then
            if (EcosystemMap(xIndex,yIndex) == 0) then
              ! This is water, so check if the interpolation flagged as 200, 
              ! or the icec data
              ! says this is ice, or the NISE has misidentified it as Perm Snow:
              if (     (icec(i,j) > 0.5 .and. icec(i,j) <= 1.0) &
                .or. (snowIceType(i,j) == 101) &
                .or. (snowIceType(i,j) == 200) ) then
                ! This is water, so this will be sea ice.  Set snowIceType to 95%
                ! so that the sea ice flag is triggered
                snowIceType(i,j) = 95
              end if
            else
              if (     (snowIceType(i,j) >= 1 .and. snowIceType(i,j) <= 100) &
                .or. (snowIceType(i,j) == 200) ) then
                ! This is snow covered land, so set to dry snow:
                snowIceType(i,j) = 103
              end if
            end if
          end if
#else
          if (EcosystemMap(xIndex,yIndex) == 0) then 
            if ( snowIceType(i,j) >= 30 .and. snowIceType(i,j) <= 100) then 
              snowIceType(i,j) = snowIceType(i,j) 
          	 else
          		snowIceType(i,j) = 0
          	 endif
#ifdef SSFR_INST
  		    if ( icec(i,j) > 0. .and. icec(i,j) <= 1.0 ) &
  				snowIceType(i,j) = nint(icec(i,j)*100)
#else
  		    if ( icec(i,j) > 0.5 .and. icec(i,j) <= 1.0 ) &
  				snowIceType(i,j) = nint(icec(i,j)*100)
#endif
  		  endif
  
#endif
  
          ! Store the ecosystem value:
          Ecosystem(i,j) = EcosystemMap(xIndex,yIndex)
  		  ! print*,xIndex,yIndex,EcosystemMap(xIndex,yIndex)
          ! WDR we use the sfc albedo from l2gen instead
          albedo_real4(:) = c2_sfc_albedo( i, j, :)
  		  		  
          ! Set the albedo to either land values, water, or snow/ice, depending
          ! on the scenario:
          if (      (snowIceType(i,j) >= 1)   &
            .and. (snowIceType(i,j) <= 100) &
            .and. (Ecosystem(i,j) == 0)     ) then
  
            Ecosystem(i,j) = 15
  
            !This is sea ice, so set the either wet or dry scenario:
                
  			 sfc_info(i,j) = 1      
                
            ! Determine arctic melting season and set albedo values:
            if (lat(i,j) >= 0.0) then
            ! Arctic melting season: June 1st (day 152) to September 1st (day 244):
            SELECT CASE (JulianDay)
              CASE ( (NHMeltBeg-TransitionDays):NHMeltBeg-1)
                ! transition period between winter and melt, so linear 
                ! fit between them
                ! to get albedo value:
                Slope (:) = (ProcessMeltSeaIceAlbedos(:) - &
                  ProcessDrySeaIceAlbedos(:)) / real(TransitionDays) 
                Intercept(:) = ProcessDrySeaIceAlbedos(:) - Slope(:)  &
                  * real(NHMeltBeg-TransitionDays) 
                TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)
                albedo_real4(:) = (TransitionAlbedo(:) *  &
                  real(snowIceType(i,j))/100.000)         &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
              CASE ( NHMeltBeg:NHMeltEnd )
                ! Melt Sea Ice Values:
                albedo_real4(:) = (ProcessMeltSeaIceAlbedos(:) * &
                  real(snowIceType(i,j))/100.000) &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
              CASE ( NHMeltEnd+1:(NHMeltEnd+TransitionDays) )
                ! transition period between melt and winter, so linear 
                ! fit between them
                ! to get albedo value:
                Slope (:) = (ProcessDrySeaIceAlbedos(:)- &
                  ProcessMeltSeaIceAlbedos(:)) / real(TransitionDays) 
                Intercept(:) = ProcessMeltSeaIceAlbedos(:) - Slope(:)  &
                  * real(NHMeltEnd) 
                TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                albedo_real4(:) = (TransitionAlbedo(:) * &
                  real(snowIceType(i,j))/100.000)         &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
              CASE DEFAULT
                ! Dry Sea Ice values::
                albedo_real4(:) = (ProcessDrySeaIceAlbedos(:) *  &
                  real(snowIceType(i,j))/100.000) &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
            END SELECT
  	
          else           
            ! Determine southern hemisphere melting season:
            ! Antarctic melting season: December 1st (day 334) to 
            ! March 1st (day 61):
            SELECT CASE (JulianDay)
              CASE ( (SHMeltBeg-TransitionDays):SHMeltBeg-1)
                ! transition period between winter and melt, so linear 
                ! fit between them
                ! to get albedo value:
                Slope (:) = (ProcessMeltSeaIceAlbedos(:)- &
                  ProcessDrySeaIceAlbedos(:))   &
                  / real(TransitionDays) 
                Intercept(:) = ProcessDrySeaIceAlbedos(:) - Slope(:)  &
                  * real(SHMeltBeg-TransitionDays) 
                TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:)    
                albedo_real4(:) = (TransitionAlbedo(:) * &
                  real(snowIceType(i,j))/100.000)         &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
                      
              CASE ( SHMeltBeg:366 )
                ! Melt Sea Ice Values:
                albedo_real4(:) = (ProcessMeltSeaIceAlbedos(:) * &
                  real(snowIceType(i,j))/100.000) &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
              CASE ( 1:SHMeltEnd )
                ! Melt Sea Ice Values:
                albedo_real4(:) = (ProcessMeltSeaIceAlbedos(:) * &
                  real(snowIceType(i,j))/100.000) &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
              CASE (SHMeltEnd+1:(SHMeltEnd+TransitionDays) )
                ! transition period between melt and winter, so linear fit 
                !  between them
                ! to get albedo value:
                Slope (:) = (ProcessDrySeaIceAlbedos(:)- &
                  ProcessMeltSeaIceAlbedos(:)) / real(TransitionDays) 
                Intercept(:) = ProcessMeltSeaIceAlbedos(:) - Slope(:) * &
                  real(SHMeltEnd) 
              
                TransitionAlbedo(:) = Slope(:) * real(JulianDay) + Intercept(:) 
                albedo_real4(:) = (TransitionAlbedo(:) * &
                  real(snowIceType(i,j))/100.000)         &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
	
              CASE DEFAULT
                ! Dry Sea Ice values::
                albedo_real4(:) = (ProcessDrySeaIceAlbedos(:) * &
                  real(snowIceType(i,j))/100.000) &
                  + (WaterAlbedo * real(100-snowIceType(i,j))/100.000)
            END SELECT
				
          end if

          Albedo(i,j,:) = nint(albedo_real4(:)*1000.)

      else if ( Ecosystem(i,j) == 15 .or. snowIceType(i,j) == 101 ) then
        ! This is permanent snow/ice, so set to either map or stat value:
        ! WDR new logic with c2_sfc_albedo
        if ((Ecosystem(i,j) == 15) .and. &
          (      c2_sfc_Albedo(i,j,1) > 0 &
          .and. c2_sfc_Albedo(i,j,1) < 1.) ) then

          ! This pixel is perm snow class and has a valid map value,
          ! so set to the map's value:
          Albedo(i,j,:) = nint( (1.0 - emissivity(i,j,1))*1000.)

          sfc_info(i,j) = 3

        else
          ! This pixel isn't flag as perm ice eco, so set it to the
          ! statistical value for perm ice:
          Ecosystem(i,j) = 15
          Albedo(i,j,:) = nint (ProcessSnowAlbStats(:, Ecosystem(i,j), 1)*1000.)
          if( ( c2_sensor_id /= OCI_ID ) .AND. ( c2_sensor_id /= OCIS_ID ) ) &
            Albedo(i,j,6) = nint( (1.0 - emissivity(i,j,1))*1000.)
          sfc_info(i,j) = 3
        end if
#ifndef GEOS5_SFC
        else if ( snowIceType(i,j) == 103 ) then
          !This is dry snow, so set the pixel to dry snow stats:
          ! if (Ecosystem(i,j) == 0) Ecosystem(i,j) = 15
          Albedo(i,j,:) = nint( ProcessSnowAlbStats(:, Ecosystem(i,j), 1) * 1000.)
          if( ( c2_sensor_id /= OCI_ID ) .AND. &
              ( c2_sensor_id /= OCIS_ID ) ) then &
            Albedo(i,j,6) = nint( (1.0 - emissivity(i,j,1))*1000.)
          Ecosystem(i,j)= 15
  
          sfc_info(i,j) = 3

#else
        ! use the snow fraction as true fraction 
        else if (snowIceType(i,j) >= 30 .and. snowIceType(i,j) <= 100 .and. &
          Ecosystem(i,j) /= 0) then 
          Albedo(i,j,:) = nint( (c2_sfc_albedo(i,j,:)* &
            (1.-snowIceType(i,j)*0.01) + &
            snowIceType(i,j)*0.01* &
            ProcessSnowAlbStats(:, Ecosystem(i,j), 1)) * 1000.)
          if( ( c2_sensor_id /= OCI_ID ) .AND. ( c2_sensor_id /= OCIS_ID ) ) &
            Albedo(i,j,6) = nint( (1.0 - emissivity(i,j,1))*1000.)
          Ecosystem(i,j)= 15
  
          sfc_info(i,j) = 3

#endif

        else if ( snowIceType(i,j) == 104 ) then
          ! This is wet snow, so set the pixel to wet snow stats:
          ! if (Ecosystem(i,j) == 0) Ecosystem(i,j) = 15
          Albedo(i,j,:) = nint (ProcessSnowAlbStats(:, Ecosystem(i,j), 2) * 1000.)
          if( ( c2_sensor_id /= OCI_ID ) .AND. ( c2_sensor_id /= OCIS_ID ) ) &
            Albedo(i,j,6) = nint ((1.0 - emissivity(i,j,1))*1000.)
          Ecosystem(i,j)= 15
  
          sfc_info(i,j) = 3
  
  
        else if (Ecosystem(i,j) == 0) then
          ! Water, so set all wavelengths to the water albedo value, 5%:
          Albedo(i,j,:) = nint(WaterAlbedo*1000.)
  
          sfc_info(i,j) = 0
        else
  
          ! This is a snow-free, ocean free, and not perm snow/ice pixel,
          ! so just store the value from the albedo map:
          if( ( c2_sensor_id == OCI_ID ) .or. ( c2_sensor_id == OCIS_ID ) )then
            ! for added 2.2 in OCI
            Albedo(i,j,:) = c2_sfc_albedo(i,j,:) * 1000.
          else
            Albedo(i,j,1:5) = c2_sfc_albedo(i,j,:) * 1000.
            ! use the 1-emissivity for the 3.7um albedo map
            Albedo(i,j,6) = nint((1.0 - emissivity(i,j,1))*1000.)
          endif
          sfc_info(i,j) = 2
        end if

      else
    
        ! Bad geolocation, so set albedos and emissivity to fill value
        Albedo(i,j,:) = -999
        emissivity(i,j,:) = -999.
   
      end if
    end do
  end do
  !*****************************************************************************
  ! Deallocate arrays:
  !*****************************************************************************
  ! WDR we'll be keeping some of these arrays
  !WDRdeallocate( AlbedoMap        )
  !WDRdeallocate( EcosystemMap     ) 
  !WDRdeallocate( emissivity_map)
  deallocate( snowIceType      )
  !WDRdeallocate(albedo_real4)

  if (status > 0 )then
    status = error
    call local_message_handler('Undetermined problem detected in getAlbedoEco',&
      status,'getAlbedoEco')
  endif	
	
  end subroutine getAlbedoEco

end module modis_albedo
