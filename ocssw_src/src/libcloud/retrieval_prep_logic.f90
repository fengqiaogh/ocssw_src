module retrieval_prep_logic

  use science_parameters

  implicit none
  
  
  integer, private :: TOTAL_POINTS, total_taus, max_allowable_tau
  real, dimension(:), allocatable, private :: taux, rfx
  
  logical :: GO_PRINT
  
contains

subroutine init_retrieval(library_taus)

  real, dimension(:), intent(in) :: library_taus
  integer :: i

  total_taus = size(library_taus)
  TOTAL_POINTS = total_taus+1
  max_allowable_tau = library_taus(total_taus)
  
  allocate(taux(TOTAL_POINTS))
  allocate(rfx(TOTAL_POINTS))

  taux(1) = 0
  do i = 1, total_taus
    taux(i+1) = library_taus(i)
  enddo

end subroutine init_retrieval

subroutine cleanup_retrieval

  deallocate(taux)
  deallocate(rfx)
  
end subroutine cleanup_retrieval

subroutine compute_water_path(tau, re, density, library_re, &
  extinction_efficiency, water_path)
  ! added by G. Wind 11.07.05 to compute the water path properly

  use modis_numerical_module
  use GeneralAuxType
  implicit none

  real, intent(in) :: tau, re, density
  real, dimension(:), intent(in) :: extinction_efficiency, library_re
  real, intent(out) :: water_path

  integer :: ilo, ihi
  real :: x(2), y(2), Qe

  ilo = 0 
  ihi = 0 

  ! determine the bounds for interpolation
  call bisectionsearch(library_re, re, ilo, ihi)
  x(1) = library_re(ilo)
  x(2) = library_re(ihi)
  y(1) = extinction_efficiency(ilo)
  y(2) = extinction_efficiency(ihi)
  
  Qe = linearinterpolation(x,y,re)
  
  water_path = ( 4. * tau * re * density ) / (3. * Qe)

end subroutine compute_water_path

subroutine vis_nonabsorbing_science(reflectance_nonabsorbing, &
  nonabsorbing_index, nonabsorbing_albedo, library_taus, &
  library_radii, sfr,fti1,fti0, rfi, theta, theta0, phi, &
  cloudtop_pressure, process, optical_thickness_vector)
  !
  ! WDR some info
  !  I believe that for each radius in the table, this finds the optical
  !  thickness values that give the closest reflectance to the non-absorbing
  !  band's reflectance.  Before doing that, it corrects the table refl for 
  !  the sfc albedo.  It also (over land) makes a (iterative for 0.66 
  !  non-abs band for land) correction to the input point reflectance for 
  !  rayleigh effects from the atmosphere, cloud and both.
  !
  !  reflectance_nonabsorbing  IN  the reflectance to get COT vector for
  !  nonabsorbing_index  IN  index of non abs band (to get right table look at)
  !  nonabsorbing_albedo  IN  the sfc albedo in non abs band
  !  library_taus  IN  tau axis values in refl table
  !  library_radii  IN  re axis values (ice or water or...)
  !  sfr  IN (re x Lam x tau) spherical albedo (ice or water or...)
  !  fti1  IN  (  x  x  ) fluxdownwater_sensor?
  !  fti0  IN  (  x  x  ) fluxdownwater_solar?
  !  rfi  IN  ( tau,lambda, re) reflectance table for pixel view geometry
  !  theta  IN  sen zen
  !  theta0  IN  sol zen
  !  phi  IN  rel azimuth
  !  cloudtop_pressure  IN  cloud top pressure
  !  process  IN  structure with phase info
  !  optical_thickness_vector  OUT  output set of COT candidates for the radii
  !      in the table

  ! sep modification w/Rayleigh scattering, 24 May 2005, etc.

  use GeneralAuxType
  use modis_cloudstructure
  use modis_numerical_module
  use mod06_run_settings
  use libraryarrays, only : reflibA, reflibB,rayleigh_liq, rayleigh_ice

  implicit none

  type(cloudphase), intent(in) :: process
  integer, intent(in)  ::  nonabsorbing_index

  real, intent(in)     ::  reflectance_nonabsorbing, nonabsorbing_albedo,          &
                           theta, theta0, phi, cloudtop_pressure
  real(single), dimension (:),     intent(in)  ::   library_taus, library_radii
  real(single), dimension (:,:,:), intent(in)  ::   sfr,  fti1, fti0, rfi                      
  
  real, intent(out)    ::   optical_thickness_vector(:)

  
  
  !local variables
  integer :: radii, i, ii
  integer :: correction_iteration, maximum_correction_iterations
  real     :: reflectance_corrected,  local_reflectance_nonabsorbing, &
    optical_thickness
  !  real :: max_allowable_tau
  !  integer :: TOTAL_POINTS
  !  real, dimension(:), allocatable  :: taux, rfx
  real, dimension(:), allocatable ::  rfi1
  
  !  parameter (max_allowable_tau = 200.)
   
  !   TOTAL_POINTS = size(library_taus)+1

   radii = size(library_radii)

   allocate( rfi1(radii) )
  !   allocate( taux(TOTAL_POINTS) )
  !   allocate( rfx(TOTAL_POINTS))


  ! add the tau=0 to the working list of taus
  !   taux(1) = 0
  !   do i = 1, total_taus
  !     taux(i+1) = library_taus(i)
  !   enddo
   
  !  calculated refl value for tauc=inf at non absorbing wavelength 
  !  for all radii
  !  now this is simply the maximum library reflectance. 
  if (.not. COX_MUNK) then 
    rfi1(1:radii) = rfi(total_taus,nonabsorbing_index,1:radii)  
  else
    rfi1(1:radii) = rfi(TOTAL_POINTS,nonabsorbing_index,1:radii)  
  endif
  
  DO i = 1, radii   

    !  sep, 25 May: this has to be inside the re loop, otherwise 
    !  local_reflectance_nonabsorbing decreases
    !  with each subsequent re, in addition to each subsequent iterations
    local_reflectance_nonabsorbing = reflectance_nonabsorbing
  
    !  Perform Rayleigh correction iteration loop only for 0.65 um band
#if ASTER_INST
    if (nonabsorbing_index==1 .or. nonabsorbing_index == 2) then!{
#else
    if (nonabsorbing_index==1) then!{
#endif
      maximum_correction_iterations = 3
    else!}{
      maximum_correction_iterations = 1   
    end if
  
    !  rfx index=1 is set to the surface albedo
    if (.not. COX_MUNK) then 
      rfx(1) = nonabsorbing_albedo
  
      !  Modify the reflectances by non-absorbing sfc albedo
      rfx(2:TOTAL_POINTS) = rfi(1:total_taus, nonabsorbing_index, i) + &
        (nonabsorbing_albedo * fti1(1:total_taus, nonabsorbing_index, i) * &
        fti0(1:total_taus, nonabsorbing_index, i)) / &
        ( 1.0 - nonabsorbing_albedo * sfr(1:total_taus, nonabsorbing_index, i))
    else
      rfx(1:TOTAL_POINTS) = rfi(1:TOTAL_POINTS, nonabsorbing_index, i)
    endif 
          
    ! save the nonabsorbing library off to the side so we can use it 
    ! if necessary during the alternate retrieval
    reflibA(1:TOTAL_POINTS, i) = rfx(:)   
            
    do correction_iteration=1, maximum_correction_iterations
    
      ! wind fix 6.24.05 correct for divide-by-zero and also indexing 
      !  changed to (i) from (1:radii)
      if ( (rfi1(i) - local_reflectance_nonabsorbing) < 0.0001) then!{
        optical_thickness_vector(i) = max_allowable_tau
        exit
      endif!}
      
      call interpolate_refl_cot (local_reflectance_nonabsorbing, rfx, taux, &
        optical_thickness)
           
      optical_thickness_vector(i) = optical_thickness
        
#if ASTER_INST
      if ( (nonabsorbing_index == 1 .or. nonabsorbing_index == 2) .and. &
        optical_thickness > 0. .and. maximum_correction_iterations /= 1 ) then!{
#else
      if (nonabsorbing_index == 1 .and. optical_thickness > 0. .and. &
        maximum_correction_iterations /= 1 ) then!{
#endif    
        call rayleighcorrection (reflectance_nonabsorbing, cloudtop_pressure, &
          process, optical_thickness, nonabsorbing_albedo, &
          fti1(:,nonabsorbing_index,i),fti0(:,nonabsorbing_index,i), &
          sfr(:,nonabsorbing_index,i), nonabsorbing_index,i,theta0,theta,phi, &
          reflectance_corrected)
  
        local_reflectance_nonabsorbing = reflectance_corrected
  
        ! collect the rayleigh-corrected reflectance for later storage to disk
        if (process%icecloud == 1) then 
          rayleigh_ice(i) = local_reflectance_nonabsorbing
        else if (process%watercloud == 1) then 
          rayleigh_liq(i) = local_reflectance_nonabsorbing
        endif
        
        if (local_reflectance_nonabsorbing > rfi1(i)) then!{ 
          optical_thickness_vector(i) = max_allowable_tau 
          exit
        endif!}
  
      else!}{
        exit
      endif!}
  
    end do   ! Rayleigh correction loop
  
    if ( (rfi1(i) - local_reflectance_nonabsorbing) < 0.0001) &
      optical_thickness_vector(i) = max_allowable_tau
      
    call interpolate_refl_cot (local_reflectance_nonabsorbing, rfx, taux, &
      optical_thickness) 
    optical_thickness_vector(i) = optical_thickness

  END DO      ! radii loop for nonabsorbing band optical thickness

  do i = 1, radii
    ! if (optical_thickness_vector(i) > 100.) optical_thickness_vector(i) = 100.
    if(rfi1(i) < reflectance_nonabsorbing )  &
      optical_thickness_vector(i) = max_allowable_tau 
  end do      

  deallocate (rfi1)
  !   deallocate(taux, rfx)

end subroutine vis_nonabsorbing_science
!

subroutine vis_absorbing_science(optical_thickness_vector,  &
  reflectance_absorbing, absorbing_index, absorbing_albedo, &
  library_taus, library_radii, sfr,fti1,fti0, rfi, &
  residual, maxradii, debug )
  !
  !  Use the tau(Re) found for the non-absorb refl and look in the 
  !  abs table at the same locations and get the fraction diff (relative 
  !  to the pixel reflectance) between
  !  the table refl and the pixel refl, output as residual.  Again, the table
  !  refl are corrected for the sfc albedo before doing this 
  !  
  !  optical_thickness_vector  IN  the best COT found for each table radius
  !  reflectance_absorbing  IN  the target abs refl to match
  !  absorbing_index   IN  which abs band being used
  !  absorbing_albedo  IN  sfc albedo
  !  library_taus  IN  tau axis in table
  !  library_radii  IN  reff axis in table
  !  sfr  IN  spherical_albedo_water
  !  fti1  IN  fluxdownwater_sensor
  !  fti0  IN  fluxdownwater_solar
  !  rfi  IN  reflectance table
  !  residual  OUT  residual for each radii again, this is the 
  !     abs refl(point) - refl(radius).  so this goes to solveretrieval 
  !     to declare the winner.
  !  maxradii  OUT  
  !  debug  IN  debug switch
  !
  ! sep, 4 May: example call
  !  WDR nore args, so ??
  !
  !        call vis_absorbing_science(optical_thickness_vector,          &
  !          band_measurements(xpoint, absorbingband_index, ypoint), &
  !          absorbingband_index,                 &
  !          surface_albedo(xpoint, ypoint,absorbingband_index), &
  !          library_wavelengths,                  &
  !          library_taus,                         &
  !          ice_radii,                            &
  !          asymmetry_ice,                        &
  !          singlescattering_ice,                 &
  !          extinction_ice,                       &
  !          ak_ice,                               &
  !          el_ice,                               &
  !          em_ice,                               &
  !          en_ice,                               &
  !          spherical_albedo_ice,                 &
  !          int_escapefuncice_solar,              &
  !          int_fluxdownice_sensor,               &
  !          int_fluxdownice_solar,                &
  !          int_escapefuncice_sensor,             &
  !          int_reflectance_ice,                  &
  !          residual,maxradii,debug)

  ! with 
  !   1) a known vector of optical thicknesses 
  !   2) absorbing band reflectance for the scene of interest
  !   3) library values for the scene of interest 
  !
  ! cloud top effective particle radius is calculated

  use GeneralAuxType
  use modis_numerical_module
  use mod06_run_settings
  use libraryarrays, only: reflibB

  implicit none

  real, intent(in) :: optical_thickness_vector(:), reflectance_absorbing,  &
                      absorbing_albedo       

  integer, intent(in)  :: absorbing_index

  real(single), intent(in) :: library_taus(:), library_radii(:), &
    sfr(:,:,:), fti1(:,:,:), fti0(:,:,:), rfi(:,:,:)

  real, intent(out)    ::   residual(:)
  logical, intent(in)  ::   debug
  integer, intent(out) ::   maxradii

  !  local variables
  real     ::  rf1

  real, allocatable    :: rfi1(:)
  real, allocatable    ::  localoptical_thickness_vector(:)
  !  real, allocatable :: taux(:)

  integer  ::  radii, i, ii
  
  !   real, dimension(:), allocatable ::  rfx
  !   integer :: TOTAL_POINTS
  integer  :: lowbound, highbound, taubounds(2)
  real     :: iftau(2), refvals(2)
  real     :: interp_tau

  lowbound = 0 ! WDR-UIV
  highbound = 0 ! WDR-UIV
 
  ! TOTAL_POINTS = size(library_taus)+1
 
  !  total_taus = size(library_taus)
  radii = size(library_radii)

  !  allocate(rfx(TOTAL_POINTS))

  allocate(rfi1(radii))
  !  allocate(taux(TOTAL_POINTS))
  allocate(localoptical_thickness_vector(radii))

  ! science fixer, for when the semi-infinite layer library reflectance 
  !  goes wiggly
  !  apparently, a (-) tau can occur?
  do i = radii, 2, -1
    maxradii = i
    if (optical_thickness_vector(i) > 0. ) exit
  enddo
  
  ! insert the zero tau into the tau vector
  !  taux(1) = 0.
  !  do i = 1, total_taus
  !     taux(i+1) = library_taus(i)
  !  enddo

  !  get reflectance residual (library - measurement)  for the absorbing band

  localoptical_thickness_vector = optical_thickness_vector

  !  calculated refl value for tauc=inf at absorbing wavelength for all radii

  if (.not. COX_MUNK) then 
    rfi1(1:radii) = rfi(total_taus,absorbing_index,1:radii)
  else
    rfi1(1:radii) = rfi(TOTAL_POINTS,absorbing_index,1:radii)
  endif

  !  reflection functions for channel 2 and 3 with fixed tc
  !  all variables must be redefined because this part is used
  !  after line-scanning for ch-1.

  !  optical thickness  for the appropriate wavelength = 
  !    qext(wavelength)/qext(0.65 um) * tauc(0.65 um)
  !  determine the reflectance at calculated optical thickness

  do i = 1,maxradii

    if (.not. COX_MUNK) then 
      rfx(1) =  absorbing_albedo

      rfx(2:TOTAL_POINTS) = rfi(1:(TOTAL_POINTS-1),absorbing_index,i)  + &
        fti1(1:(TOTAL_POINTS-1),absorbing_index,i) * &
        fti0(1:(TOTAL_POINTS-1),absorbing_index,i) * &
        absorbing_albedo /            &
        ( 1.-absorbing_albedo*sfr(1:(TOTAL_POINTS-1),absorbing_index,i))
    else
      rfx(1:TOTAL_POINTS) = rfi(1:TOTAL_POINTS,absorbing_index,i) 
    endif

    ! save the absorbing library off to the side so we can use it if necessary 
    ! during the alternate retrieval
    reflibB(1:TOTAL_POINTS, i) = rfx(:)   

    call bisectionsearch(taux, localoptical_thickness_vector(i), &
      lowbound, highbound)
    taubounds(1) = lowbound
    taubounds(2) = highbound
    iftau(1)    = taux(taubounds(1))
    iftau(2)    = taux(taubounds(2))

    refvals(1) = rfx(lowbound)
    refvals(2) = rfx(highbound)

    rf1 = linearinterpolation(iftau,refvals,localoptical_thickness_vector(i))

    residual(i) = rf1/reflectance_absorbing- 1.0

  enddo


  deallocate(rfi1, localoptical_thickness_vector)
  !  deallocate(taux, rfx)

end subroutine vis_absorbing_science

subroutine interpolate_refl_cot(reflectance, reflectance_vector, &
  optical_thickness_vector, optical_thickness)
  !
  !  reflectance
  !  reflectance_vector
  !  optical_thickness_vector
  !  optical_thickness
  ! 
  !  We need to interpolate the reflectance vector to get the optical 
  !  thickness that corresponds to the reflectance given.  We cannot 
  !  assume however that the reflectance
  !  is uniformly increasing as a function of the optical thickness.  
  !  There are cases
  !  where there maybe multiple solutions of the optical thickness 
  !  for a given reflectance
  !  In this case we take the larger solution and move on.

  use GeneralAuxType
  use modis_numerical_module       

  implicit none

  real, intent(in) :: reflectance, reflectance_vector(:), &
    optical_thickness_vector(:)
  real, intent(out):: optical_thickness

  integer  :: refl_size, tau_size
  integer  :: lowbound, highbound, refbounds(2)
  real     :: ifref(2), tauvals(2)
  real     :: interp_tau
   
  lowbound = 0  ! WDR-UIV
  highbound = 0 ! WDR-UIV
   
  refl_size = size(reflectance_vector)

  if (reflectance < reflectance_vector(1) ) then!{
    !     reflectance is out of range. no interpolation is possible
    optical_thickness = - 1.
    return
  endif!}
      
  ! reflectance is super-bright, thick cloud, so just set to max 
  !  tau and be done with it, 
  ! no interpolation possible, so be it. 
  if ( reflectance > reflectance_vector(refl_size) ) then 
    tau_size = size(optical_thickness_vector)
    optical_thickness = optical_thickness_vector(tau_size)
    return
  endif
    
  ! doesn't matter nothing. We will use the standard linear 
  ! interpolation with the library that big.   
  ! we can probably get away with it. 

  call bisectionsearch(reflectance_vector, reflectance, lowbound, highbound)
  refbounds(1) = lowbound
  refbounds(2) = highbound
  ifref(1)    = reflectance_vector(refbounds(1))
  ifref(2)    = reflectance_vector(refbounds(2))

  tauvals(1) = optical_thickness_vector(lowbound)
  tauvals(2) = optical_thickness_vector(highbound)
  optical_thickness = linearinterpolation(ifref,tauvals,reflectance)

end subroutine interpolate_refl_cot

subroutine rayleighcorrection (reflectance, cloudtoppressure, process,  &
  optical_thickness, nonabsorbing_galbedo, fti1, fti0, sfr, iw, ir,     &
  solarzenith, sensorzenith, azimuth, reflectance_corrected)
  !
  !  reflectance  IN  non-absorbing reflectance, I believe only the 0.66 
  !       um band value is used
  !  cloudtoppressure  IN  pressure of cloud top
  !  process  IN  processing info, has cloud phase
  !  optical_thickness  IN  optical thickness at current point
  !  nonabsorbing_galbedo  IN  surface albedo
  !  fti1  IN  fluxdownwater_sensor at the current table location
  !  fti0  IN  fluxdownwater_solar at the current table location
  !  sfr  IN spherical albedo (ice or water or...)
  !  iw  IN  non absorbing index, probably 1 for .66 um
  !  ir  IN  table radii index
  !  solarzenith, sensorzenith, azimuth  point 's view geometry
  !  reflectance_corrected  OUT  final reflectance, rayleigh corrected
  !
  ! Versioning history: 
  !
  !  Initial, but incomplete, portions written by M. Gray
  !  Modified by B. Wind, G. Wind to add cloud albedo
  !  5-16-05, etc.: Modified by S. Platnick to include non-zero 
  !  surface albedo and reflectances in the asymptotic regime,
  !  as well as changes to loop logic/structure in vis_nonabsorbing_science
  !  This algorithm is based on the Rayleigh correction algorithm 
  !  described in Wang and King (1997). The 
  !  Rayleigh correction algorithm described in the paper is derived 
  !  for a black surface and describes
  !  results for water clouds only.  Here the algorithm has been 
  !  extended to a non-zero surface albedo 
  !  and uses forward library water and ice cloud fluxes based on 
  !  the decision tree cloud phase.

  ! VARIABLE TRANSLATION
  !
  ! nonabsorbing_galbedo = sfc albedo in nonabsorbing band; scalar
  ! fti0, fti1     = transmitted fluxes interpolated to mu0, mu; 
  !     indexing: (tau index)
  ! sfr            = spherical albedo; indexing: (tau index)
  ! iw, ir         = band and effective radius indices, respectively
 
   use GeneralAuxType
   use libraryarrays
   use modis_cloudstructure, only: cloudphase

   implicit none

  !intent in/out variables
  type(cloudphase), intent(in)  :: process
  integer, intent(in) :: iw, ir
  real, intent(in)    ::  reflectance, cloudtoppressure, &
    optical_thickness, nonabsorbing_galbedo, &
    solarzenith, sensorzenith, azimuth
  real, intent(in)    :: sfr(:), fti0(:), fti1(:)

  real, intent(out)   :: reflectance_corrected

  !local variables
  real  :: taur, rd, cloud_albedo
  real  :: mu0, mu, flux0, flux1, pray, x, refl_s
  real  :: local_int_fluxup_solar,local_int_fluxup_sensor

  real, parameter :: surfacepressure = 1013.  ! Fixed sfc P!
  real, parameter :: Cm = 0.84
  real, parameter :: taur0(2) = (/0.044, 0.025/)
  !  the second position of taur0 is for ASTER only, none of 
  !  the other instruments would use it, 
  !  because none of the other instruments do rayleigh correction 
  !  for non-0.65um band

  local_int_fluxup_solar=0.; local_int_fluxup_sensor=0.
  reflectance_corrected=0.

  !  Rayleigh optical depth from TOA to cloud-top
  !  WDR TOA-TOC thickness, Wang&King Eq 12
  taur = taur0(iw)*cloudtoppressure/surfacepressure

  mu0 = cos(d2r*solarzenith)
  mu  = cos(d2r*sensorzenith)

  if(process%icecloud == 1) then!{
    call interp_lib_reflFlux_cloudAlbedo( mu0, mu, optical_thickness, &
      nonabsorbing_galbedo, sfr, fti1, fti0, iw, ir, &
      library_fluxsolarzenith, library_fluxsensorzenith, flux_up_ice_solar, &
      flux_up_ice_sensor, local_int_fluxup_solar, local_int_fluxup_sensor)
  else!}{
    call interp_lib_reflFlux_cloudAlbedo(mu0, mu, optical_thickness, &
      nonabsorbing_galbedo, sfr, fti1, fti0, iw, ir, &
      library_fluxsolarzenith, library_fluxsensorzenith, &
      flux_up_water_solar, flux_up_water_sensor, local_int_fluxup_solar, &
      local_int_fluxup_sensor)
  endif!}


  !  Rayleigh scattering phase function, where x=scattering angle 
  !  (Eq. 3, Wang and King)
  !  sep, 13 May: confirmed that pi is not needed in pray formula 
  !  since it is for reflectance, not intensity.
  !  sep, 16 May: w/KIng, confirmed that this formula is correct 
  !  for the usual definition of relative azimuth
  !  and mu, mo0 defined to be always positive.

  x = -mu0*mu + sqrt((1.-mu0**2)*(1.-mu**2))*cos(d2r*azimuth)
  pray = 3.*(1.+x**2)/4.0

  !  compute single scattering reflectance from molecules 
  !  bounded by cloud layer (portion of Eq. 11, Wang and King)

  refl_s = taur*pray/(4.*mu0*mu) + &
    0.5*taur* local_int_fluxup_sensor * exp(-taur/mu )/mu0 + &
    0.5*taur* local_int_fluxup_solar  * exp(-taur/mu0)/mu

  !  remove rayleigh scattering effects (Eq. 11, Wang and King)

  reflectance_corrected = (reflectance - refl_s)*exp(Cm*taur*(1./mu0 + 1./mu))

end subroutine rayleighcorrection

subroutine interp_lib_reflFlux_cloudAlbedo(miu0, miu, &
  optical_thickness, nonabsorbing_galbedo, sfr, fti1, fti0, iw, ir, &
  fluxsolarzenithangles, fluxsensorzenithangles,         &
  fluxup_solar, fluxup_sensor, interp_fluxup_solar, interp_fluxup_sensor)
  !
  !  miu0  IN  cos( solz )
  !  miu   IN  cos( senz )
  !  optical_thickness  IN  optical thickness at current point
  !  nonabsorbing_galbedo  IN  surface albedo
  !  sfr   IN spherical albedo (ice or water or...)
  !  fti1  IN  fluxdownwater_sensor at the current table location
  !  fti0  IN  fluxdownwater_solar at the current table location
  !  iw    IN  non absorbing index, probably 1 for .66 um
  !  ir    IN  table radii index
  !  fluxsolarzenithangles  IN list of the solz angles in the table
  !  fluxsensorzenithangles  IN  list of the senz angles in the table
  !  fluxup_solar  IN  array of table up sol flux (n_solz, n_tau, n_wav, 
  !      n_ice_rad )
  !  fluxup_sensor  IN  array of table up sen flux (n_solz, n_tau, n_wav, 
  !      n_ice_rad )
  !  interp_fluxup_solar  OUT  sol flux at the solz, tau, wave, and radius
  !  interp_fluxup_sensor OUT  sen flux at the solz, tau, wave, and radius

  use GeneralAuxType
  use modis_numerical_module
  use libraryarrays
  implicit none

  !intent in/out variables
  integer, intent(in) :: iw, ir
  real, intent(in)    ::  miu0, miu, optical_thickness, nonabsorbing_galbedo
  real, intent(in)    ::  sfr(:), fti0(:), fti1(:)
  real, dimension(:), intent(in) :: fluxsolarzenithangles, &
    fluxsensorzenithangles
  real, dimension(:,:,:,:), intent(in) :: fluxup_solar, fluxup_sensor

  real, intent(out)  :: interp_fluxup_solar, interp_fluxup_sensor

  !local variables
  integer  :: lowbound, highbound, its, nts, taubounds(2)
  real     :: scaledTau, iftaus(2), angles(2), functionvalues(2)
  real     :: q1, interp_fti0, interp_fti1, interp_sfr
  real, dimension(:), allocatable :: fluxup_solar_tau,fluxup_sensor_tau

  !  ---------------------------------------------------------------
  !  Subroutine to interpolate .65 micron band cloud albedo library 
  !  in angle and tau space
  !  --------------------------------------------------------------

  !  sep NOTES: 
  !
  ! 'bisectionsearch' subroutine (looks like numerical recipe code) 
  ! and function 'linearinterpolation'
  !   both exist in modis_numerical_module.f90
  !       
  !  Get the cloud-top albedo, including the effect of a non-zero surface albedo
  !
  !   For non-asymptotic regime:
  !
  !    With a non-zero surface albedo, the cloud-top albedo is:
  !
  !      A(tau, mu, galb) = A(tau, mu, galb=0) + 
  !         T(tau, mu, galb=0)*galb*Tsph(tau)/(1 - Rsph(tau)*galb)
  !
  !       where,
  !
  !       A(tau, mu0, galb=0) = local_int_fluxup_solar
  !       A(tau, mu,  galb=0) = local_int_fluxup_sensor
  !       T(tau, mu,  galb=0) = transmitted flux = fti1/0 
  !            = interp_fluxdown_sensor/solar
  !       Rsph(tau)           = "spherical" albedo for galb=0 
  !       Tsph(tau)           = "spherical" transmittance for galb=0 = 1 
  !              - Rsph for conservative scattering
  !
  !   For asymptotic regime, the cloud-surface albedo for conservative 
  !   scattering can be simply written as (from M. King notes):
  !
  !     A(tau, mu, galb) = 
  !        1 - 4*(1-galb)*K(mu0)/[3*(1-galb)*(1-g)*(tau+2q0) + 4*galb]
  
  !    NOTE: from modis_frontend_module.f90:1522:  
  !     allocate(flux_up_ice_solar (number_fluxsolarzenith, number_taus, 
  !     number_wavelengths, number_iceradii))

   
  interp_fluxup_solar=0.; interp_fluxup_sensor=0.

  !  allocations and initializations
  nts = total_taus !size(library_taus)
  allocate(fluxup_solar_tau(nts),fluxup_sensor_tau(nts))
  interp_fluxup_solar = 0.; interp_fluxup_sensor = 0.
  taubounds = 0; lowbound=0; highbound=0; iftaus=0. 

  !  -------------------------------------------------------
  !  NON-ASYMPTOTIC REGIME: CLOUD-TOP ALBEDO w/BLACK SURFACE
  !  -------------------------------------------------------

  !  sep, 24May: scaledTau<0.5 was not being explicitly counted 
  !  for in non-asymptotic calculations. 
  !  Fortunately (or by design) the values returned for 
  !  library_taus(0) and fluxup_solar_tau(0) = 0,
  !  despite the arrays not being declared for a zero index. 
  !  Following modified to treat scaledTau<0.5 explicitly.

  !
  !  Tau index bounds for the scaled optical thickness
  !
  if (optical_thickness < library_taus(1)) then!{
    taubounds(1) = 1
    iftaus(1)    = 0.
    iftaus(2)    = library_taus(taubounds(1))
  else!}{
    call bisectionsearch(library_taus, optical_thickness, lowbound, highbound)
    taubounds(1) = lowbound
    taubounds(2) = highbound
    iftaus(1)    = library_taus(taubounds(1))
    iftaus(2)    = library_taus(taubounds(2))
  end if

  !
  !  Interpolation of fti0, fti1, and srf to the scaled optical 
  !  thickness (for incorporating surface reflectance)
  !
  if (nonabsorbing_galbedo > 0.) then 
   
    if (optical_thickness < library_taus(1)) then!{
      functionvalues(1) = 1.
      functionvalues(2) = fti0(1)
      interp_fti0 = linearinterpolation(iftaus,functionvalues,optical_thickness)
      functionvalues(1) = 1.
      functionvalues(2) = fti1(1)
      interp_fti1 = linearinterpolation(iftaus,functionvalues,optical_thickness)
      functionvalues(1) = 0.
      functionvalues(2) = sfr(1)
      interp_sfr  = linearinterpolation(iftaus,functionvalues,optical_thickness)
    else!}{
      functionvalues(1) = fti0(lowbound)
      functionvalues(2) = fti0(highbound)
      interp_fti0 = linearinterpolation(iftaus,functionvalues,optical_thickness)
      functionvalues(1) = fti1(lowbound)
      functionvalues(2) = fti1(highbound)
      interp_fti1 = linearinterpolation(iftaus,functionvalues,optical_thickness)
      functionvalues(1) = sfr(lowbound)
      functionvalues(2) = sfr(highbound)
      interp_sfr  = linearinterpolation(iftaus,functionvalues,optical_thickness)
    end if
  
  else
   
  endif

  !
  !  Interpolation for solar zenith angles and scaled optical thickness
  !
  call bisectionsearch(library_fluxsolarzenith, miu0, lowbound, highbound)
  angles(1) = fluxsolarzenithangles(lowbound)
  angles(2) = fluxsolarzenithangles(highbound)

  do its = 1,nts
    functionvalues(1) = fluxup_solar(lowbound,its, iw,ir)
    functionvalues(2) = fluxup_solar(highbound,its,iw,ir)
    fluxup_solar_tau(its) = linearinterpolation(angles,functionvalues,miu0)
  enddo

  if (optical_thickness < library_taus(1)) then!{
    functionvalues(1) = 0.
    functionvalues(2) = fluxup_solar_tau(taubounds(1))
  else!}{
    functionvalues(1) = fluxup_solar_tau(taubounds(1))
    functionvalues(2) = fluxup_solar_tau(taubounds(2))
  end if
  interp_fluxup_solar = linearinterpolation(iftaus,functionvalues, &
    optical_thickness)

  !
  !  Interpolation for sensor zenith angles and scaled optical thickness
  !
  call bisectionsearch(library_fluxsensorzenith, miu, lowbound, highbound)
  angles(1) = fluxsensorzenithangles(lowbound)
  angles(2) = fluxsensorzenithangles(highbound)

  do its = 1,nts
    functionvalues(1) = fluxup_sensor(lowbound, its,iw,ir)
    functionvalues(2) = fluxup_sensor(highbound,its,iw,ir)
    fluxup_sensor_tau(its) = linearinterpolation(angles,functionvalues,miu)
  enddo

  if (optical_thickness < library_taus(1)) then!{
    functionvalues(1) = 0.
    functionvalues(2) = fluxup_sensor_tau(taubounds(1))
  else!}{
    functionvalues(1) = fluxup_sensor_tau(taubounds(1))
    functionvalues(2) = fluxup_sensor_tau(taubounds(2))
  end if
  interp_fluxup_sensor = linearinterpolation(iftaus,functionvalues, &
    optical_thickness)


  !  --------------------------------
  !  NON-ASYMPTOTIC REGIME: w/SURFACE
  !  --------------------------------

  !  A(tau, mu, galb) = A(tau, mu, galb=0) + 
  !     T(tau, mu, galb=0)*galb*Tsph(tau)/(1 - Rsph(tau)*galb)
  if (nonabsorbing_galbedo > 0.) then 

    interp_fluxup_solar  = interp_fluxup_solar + &
      interp_fti0*nonabsorbing_galbedo*(1-interp_sfr)/ &
      (1 - interp_sfr*nonabsorbing_galbedo)
    interp_fluxup_sensor = interp_fluxup_sensor + &
      interp_fti1*nonabsorbing_galbedo*(1-interp_sfr)/ &
      (1 - interp_sfr*nonabsorbing_galbedo)
  else
  !  cox-munk  
  !  we do nothing, the surface assumed not to contribute 
  !  anything in this case. This would be the case 
  !  when 0.86um band saturated. 
  endif

  deallocate(fluxup_solar_tau,fluxup_sensor_tau)

end subroutine interp_lib_reflFlux_cloudAlbedo

subroutine nir_absorbing_science(platform_name, optical_thickness_vector,  &
  reflectance_absorbing, absorbing_index, absorbing_albedo, &
  xpoint, ypoint, CTT, thermal_trans_1way, thermal_trans_2way, &
  library_taus, library_radii, sfr,fti1,fti0,fri1, rfi,  cl_emis, sf_emis, &
  residual, maxradii, channel_number_37, emission_uncertainty_pw, &
  emission_uncertainty_Tc, sigma_R37_PW, debug)

  ! sep, 4 May: example call
  !
  !   call nir_absorbing_science(platform_name, &
  !     optical_thickness_vector,          &
  !     band_measurements(xpoint, absorbingband_index, ypoint), &
  !     absorbingband_index - 1,                 &
  !     surface_albedo(xpoint, ypoint,absorbingband_index), &
  !     cloud_top_temperature(xpoint, ypoint), &
  !     surface_temperature(xpoint, ypoint), &
  !     thermal_trans_1way,            &
  !     thermal_trans_2way,            &
  !     solar_zenith_angle(xpoint,ypoint),    &
  !     library_taus,                         &
  !     ice_radii,                            &
  !     spherical_albedo_ice,       sfr        &
  !     int_fluxdownice_sensor,     fti1          &
  !     int_fluxdownice_solar,      fti0       &
  !     int_fluxupice_sensor,        fr1          &
  !     int_reflectance_ice,        rfi          &
  !     int_cloud_emissivity_water, &
  !     int_surface_emissivity_water, &
  !     residual, maxradii, channel_37, debug)


  use GeneralAuxType
  use mod06_run_settings
  use science_parameters
  use libraryarrays, only : reflibB
  use planck_functions
  use core_arrays, only : thermal_correction_twoway_low, &
    thermal_correction_twoway_high, thermal_correction_oneway_low, &
    thermal_correction_oneway_high, Transprime_2way, Transprime_1way, &
    Tc_low_for_delta, &
    Tc_high_for_delta, surface_temperature, solar_zenith_angle, const_C, &
    Bprime_Ts, Bprime_Tc, abovecloud_watervapor
     
  implicit none

  character*(*),intent(in)::  platform_name

  real, intent(in) :: optical_thickness_vector(:), reflectance_absorbing,  &
                      absorbing_albedo

  integer, intent(in) ::  absorbing_index, channel_number_37, xpoint, ypoint

  real, intent(in)    :: thermal_trans_2way,thermal_trans_1way, CTT

  real(single), intent(in) :: library_taus(:), library_radii(:), &
    sfr(:,:,:), fti1(:,:,:), fti0(:,:,:), fri1(:,:,:), &
    rfi(:,:,:), cl_emis(:,:,:), sf_emis(:,:,:)
    
  real, intent(inout) ::   emission_uncertainty_pw(:), &
    emission_uncertainty_Tc(:), sigma_R37_PW(:)
  logical, intent(in)  ::  debug
  real, intent(out)    ::   residual(:)
  integer, intent(out) ::   maxradii

  integer              :: radii, wavelengths,i
  real                 :: ReflectanceSolarMeasured,rf1, rtherm37
  real                 :: tc, local_trans_1way, local_trans_2way
  real, allocatable    :: localoptical_thickness_vector(:), rfi1(:)
  !  real, allocatable ::  taux(:)

  real :: rtherm37_high, rtherm37_low, rsm_low, rsm_high
  real :: emission_low, emission_high, emission_reg
  real :: solar_zenith

  real :: Es, Ec, Es_gnd, B_Tc, sigmaZ1, sigmaZ2
  real :: B_Tg, B_Tc_low, B_Tc_high, Z2_term, Trans_ratio, sigmaPW_squared

  !  integer :: TOTAL_POINTS
  ! TOTAL_POINTS = size(library_taus) + 1

  solar_zenith = solar_zenith_angle(xpoint, ypoint)

  !  total_taus = size(library_taus)
  radii = size(library_radii)

  ! Early exit criteria.  If we have no thermal information a retrieval is
  ! impossible

  if (CTT < 0. ) then!{
    residual  =-999.
    return
  endif!}    

  allocate(localoptical_thickness_vector(radii))
  allocate(rfi1(radii))
  !  allocate(taux(TOTAL_POINTS))

  !  taux(1) = 0.
  !  do i = 1, total_taus 
  !     taux(i+1) = library_taus(i)
  !  enddo
 

  ! science fixer, for when the semi-infinite layer library 
  ! reflectance goes wiggly
  do i = radii, 2, -1
    maxradii = i
    if (optical_thickness_vector(i) > 0. ) exit
  enddo

  if (.not. COX_MUNK) then 
  rfi1(1:radii) = rfi(total_taus,absorbing_index,1:radii)
  else
  rfi1(1:radii) = rfi(TOTAL_POINTS,absorbing_index,1:radii)  
  endif

  localoptical_thickness_vector = optical_thickness_vector

  if (thermal_trans_1way < 0. .or. thermal_trans_1way > 1.) then!{
    local_trans_1way = 1.
  else!}{
    local_trans_1way = thermal_trans_1way
  endif!}

  if (thermal_trans_2way < 0. .or. thermal_trans_2way > 1.) then!{
    local_trans_2way = 1.
  else!}{
    local_trans_2way = thermal_trans_2way
  endif!}

  ! why call modis_planck ten billion times for each Re when this 
  !  number never changes
  ! as long as we're processing a single pixel. 

  B_Tc = modis_planck(platform_name, CTT, channel_number_37, 1)
  B_Tg = modis_planck(platform_name, surface_temperature(xpoint, ypoint), &
    channel_number_37, 1)
  B_Tc_low = modis_planck(platform_name, Tc_low_for_delta, &
    channel_number_37, 1)
  B_Tc_high = modis_planck(platform_name, Tc_high_for_delta,  &
    channel_number_37, 1)
  sigmaZ1 = const_C * (reflectance_absorbing / (local_trans_2way**2)) &
    * Transprime_2way * &
    ((2.*watervapor_error)*abovecloud_watervapor(xpoint, ypoint))

  do i = 1, maxradii
    
  if (.not. COX_MUNK) then 
    call toa_radiance37 (platform_name, taux, &
      localoptical_thickness_vector(i), sfr(:,absorbing_index,i), &
      rfi1(i), fti0(:,absorbing_index,i), fti1(:,absorbing_index,i),  &
      fri1(:,absorbing_index,i), rfi(:,absorbing_index,i), &
      absorbing_albedo, B_Tg, B_Tc, rf1, &
      rtherm37, channel_number_37, reflibB(:,i), Es, Ec)

    Es_gnd = 1. - absorbing_albedo
  else
    call toa_radiance37_cox_munk(platform_name, taux, &
      localoptical_thickness_vector(i), B_Tg, B_Tc, &
      rfi(:, absorbing_index, i), cl_emis(:,1,i), sf_emis(:,1,i), rf1, & 
      rtherm37, channel_number_37, Es, Ec)

    Es_gnd = Es
    reflibB(:,i) = rfi(:, absorbing_index, i)
  endif
  !   thermal component of the top of the atmosphere measurement 
  !  in radiance units
    
  emission_reg = rtherm37
  rtherm37 = rtherm37 * local_trans_1way

  Trans_ratio = Transprime_1way / local_trans_1way - Transprime_2way /  &
    local_trans_2way
  sigmaPW_squared = ((2.*watervapor_error)* &
    abovecloud_watervapor(xpoint, ypoint))**2
  Z2_term = (Es * Bprime_Ts)**2 * (2*delta_Ts)**2 + &
    (    Ec*Bprime_Tc + ( Es_gnd + Ec*B_Tc) * Trans_ratio )**2 *  &
    sigmaPW_squared 

  sigmaZ2 = const_C * (local_trans_1way / local_trans_2way) * sqrt(Z2_term)  
  

  ! keeping it squared as we always need it squared later. 
  sigma_R37_pw(i) = sigmaZ1**2 + sigmaZ2**2 - 2*abs(sigmaZ1*sigmaZ2)
  

  ! solar component of the top of the atmosphere measurement in 
  ! reflectance units

  ! Fri Sep 27, 2002: 3.7 um solar iradiance changed to 11.297 from 10.715 (sep)
  ! Th, June 2, 2005: 3.7 um solar iradiance  changed from 11.297 to 
  !    11.739 W/m2-um based on Juan Fontenla (LASP, CU) solar model

  !   above cloud atmosphere corrected  solar component
  !   of the top of the atmosphere measurement in reflectance units
#ifdef SIM_NORAD
  ReflectanceSolarMeasured = const_C * reflectance_absorbing !11.739)
#else
  ReflectanceSolarMeasured = const_C * (reflectance_absorbing-rtherm37) !11.739)
#endif
  ReflectanceSolarMeasured = ReflectanceSolarMeasured / local_trans_2way
  
  
  !   Compare this with the library value rf1
  residual(i) = ( rf1 / ReflectanceSolarMeasured ) - 1.

  !   rtherm37 > reflectance_absorbing is not physical. Strictly speaking
  !   the answer here would be infinity.

#ifndef SIM_NORAD
  if (rtherm37 .gt. reflectance_absorbing) then!{
    residual(i) = 10000.
  endif!}
#endif
  
  enddo

  deallocate( rfi1, localoptical_thickness_vector)
!  deallocate(taux)
  
end subroutine nir_absorbing_science

subroutine toa_radiance37(platform_name, taux, tc, sfr, rfi1, &
  fti0, fti1, fri1, rfi, galbedo, B_Tg, B_Tc, rf1,  &
  rtherm37, channel_number_37, reflib, Es, Ec )
  !
  !  platform_name
  !  taux
  !  tc
  !  sfr
  !  rfi1    
  !  fti0
  !  fti1
  !  fri1
  !  rfi
  !  galbedo
  !  B_Tg
  !  B_Tc
  !  rf1
  !  rtherm37
  !  channel_number_37
  !  reflib
  !  Es
  !  Ec
  !call toa_radiance37 (platform_name, &
  !                         taux, &
  !            local_optical_thickness_vector(i), &
  !                         sfr(:,absorbing_index,i), &
  !            rfi1(i),                    &
  !                         fti0(:,absorbing_index,i), &
  !            fti1(:,absorbing_index,i),  &
  !                         fri1(:,absorbing_index,i), &
  !            rfi(:,absorbing_index,i), &
  !                         absorbing_albedo,                   &
  !                         surface_temperature,&
  !            cloud_top_temperature, &
  !                         rf1,                                        &
  !                         rtherm37)

  use GeneralAuxType
  use modis_numerical_module
  use mod06_run_settings
  implicit none

  character*(*),intent(in)        ::  platform_name

  real, intent(in) :: taux(:)

  real(single), intent(in) :: sfr(:),   &
                              fti1(:), fri1(:), fti0(:),rfi(:)
  real, dimension(:), intent(inout) :: reflib(:)
  real, intent(in)  :: tc, galbedo, rfi1, B_Tg, B_Tc
  integer, intent(in) :: channel_number_37

  real, intent(out) :: rtherm37, rf1, Es, Ec
  
  integer   ::  i, nts
  real      :: sfr1,rthermc, rthermg, frefl1,  trans
  
  real, dimension(:), allocatable :: tranx, frefl,  local_sfr
  !  real, dimension(:), allocatable :: rfx
  
  integer  :: lowbound, highbound, taubounds(2)
  real     :: iftaus(2), functionvalues(2)  
  
  real intensity, intensity_g
  !  integer :: TOTAL_POINTS

  lowbound = 0  ! WDR-UIV
  highbound = 0  ! WDR-UIV
  
  
  !  TOTAL_POINTS = size(taux)
  nts = TOTAL_POINTS - 1

  allocate(tranx(TOTAL_POINTS), frefl(TOTAL_POINTS), local_sfr(TOTAL_POINTS))
  !  allocate(rfx(TOTAL_POINTS))
  
  !note: previous developer suggests finding a 
  !better way to pass these variables. no reason was given
  
  !compute reflectance, transmittance, cloud albedo, 
  !and spherical albedo of cloud layer bounded by 
  !lambertian sea surface at the top for band 3.7 micron.
  local_sfr(1) = 0.
  local_sfr(2:TOTAL_POINTS) = sfr(1:nts)

  ! set reflectance, transmissivity and upward flux for optical thickness == 0.
  !   sep, 3 May: this rfx zero tau value is for a black surface, 
  !  but rfx(2),... include the surface albedo as they should
  !   rfx(1)   = 0.
  rfx(1)   = galbedo

  tranx(1) = 1.0 
  frefl(1) = 0.
 
  do i = 1, nts
  
    rfx(i+1) = rfi(i)+fti1(i)*fti0(i)*galbedo/ (1.-galbedo*sfr(i))
  
    ! downward flux
    tranx(i+1) = fti1(i)                  
  
    !  upward flux
    frefl(i+1) = fri1(i)                  
  
  enddo

  reflib = rfx

  !
  !  Tau index bounds for the scaled optical thickness
  !
  call bisectionsearch(taux, tc, lowbound, highbound)
  taubounds(1) = lowbound
  taubounds(2) = highbound
  iftaus(1)    = taux(taubounds(1))
  iftaus(2)    = taux(taubounds(2))
  !
  !  Interpolation of fti0, fti1, and srf to the scaled optical 
  !  thickness (for incorporating surface reflectance)
  !
   
  functionvalues(1) = rfx(lowbound)
  functionvalues(2) = rfx(highbound)
  rf1 = linearinterpolation(iftaus,functionvalues,tc)
  functionvalues(1) = frefl(lowbound)
  functionvalues(2) = frefl(highbound)
  frefl1 = linearinterpolation(iftaus,functionvalues,tc)
  functionvalues(1) = tranx(lowbound)
  functionvalues(2) = tranx(highbound)
  trans  = linearinterpolation(iftaus,functionvalues,tc)
  functionvalues(1) = local_sfr(lowbound)
  functionvalues(2) = local_sfr(highbound)
  sfr1  = linearinterpolation(iftaus,functionvalues,tc)


  ! emission from the cloud (radiance units) at 3.7um
  intensity = B_Tc
  rthermc = (1.- trans - frefl1) * intensity

  Ec = (1.- trans - frefl1)

  ! surface emission from the ground (radiance units) at 3.7um
  !  intensity = modis_planck(platform_name, ttop, 20, 1)
  ! *** sep, 2 May: should use gtemp for sfc instead of ttop. Fixed on 3 May.
  intensity_g = B_Tg
  rthermg = trans*(1.-galbedo) * intensity_g / (1.-sfr1*galbedo)
  
  Es = trans*(1.-galbedo) / (1.-sfr1*galbedo)
  
  ! total thermal radiance contribution at band 3.7 micron
  rtherm37 = rthermc + rthermg


  deallocate(tranx, frefl, local_sfr)
  ! deallocate(rfx)

end subroutine toa_radiance37

subroutine toa_radiance37_cox_munk(platform_name, taux, tc, B_Tg, &
  B_Tc, rfi, cl_emis, sf_emis, rf1, rtherm37, channel_number_37, Es, Ec)
  !
  !  platform_name
  ! taux
  !  tc
  !  B_Tg
  !  B_Tc
  !  rfi
  !  cl_emis
  !  sf_emis\
  !  rf1
  !  rtherm37
  !  channel_number_37
  !  Es
  !  Ec
  !
  use GeneralAuxType
  use modis_numerical_module
  use mod06_run_settings
  implicit none

  character*(*),intent(in)        ::  platform_name
  real, intent(in) :: taux(:)
  real(single), intent(in) :: cl_emis(:), sf_emis(:), rfi(:)
  real, intent(in)  :: tc, B_Tg, B_Tc
  integer, intent(in) :: channel_number_37
  real, intent(out) :: rtherm37, rf1, Es, Ec
  
  integer   ::  i, nts
  real      :: surface_emissivity, cloud_emissivity, rthermc, rthermg
  
  integer  :: lowbound, highbound, taubounds(2)
  real     :: iftaus(2), functionvalues(2)  
  
  real intensity, intensity_g
  !  integer :: TOTAL_POINTS

  lowbound = 0  ! WDR-UIV
  highbound = 0  ! WDR-UIV
  
  !  TOTAL_POINTS = size(taux)
  nts = TOTAL_POINTS
  !
  !  Tau index bounds for the scaled optical thickness
  !
  call bisectionsearch(taux, tc, lowbound, highbound)
  taubounds(1) = lowbound
  taubounds(2) = highbound
  iftaus(1)    = taux(taubounds(1))
  iftaus(2)    = taux(taubounds(2))
  !
  !  Interpolation of fti0, fti1, and srf to the scaled optical 
  !  thickness (for incorporating surface reflectance)
  !
  functionvalues(1) = rfi(lowbound)
  functionvalues(2) = rfi(highbound)
  rf1 = linearinterpolation(iftaus,functionvalues,tc)
  functionvalues(1) = cl_emis(lowbound)
  functionvalues(2) = cl_emis(highbound)
  cloud_emissivity = linearinterpolation(iftaus,functionvalues,tc)
  functionvalues(1) = sf_emis(lowbound)
  functionvalues(2) = sf_emis(highbound)
  surface_emissivity = linearinterpolation(iftaus,functionvalues,tc)

  ! emission from the cloud (radiance units) at 3.7um
  intensity = B_Tc
  rthermc = cloud_emissivity * intensity

  Ec = cloud_emissivity

  ! surface emission from the ground (radiance units) at 3.7um
  !  intensity = modis_planck(platform_name, ttop, 20, 1)
  ! *** sep, 2 May: should use gtemp for sfc instead of ttop. Fixed on 3 May.
  intensity_g = B_Tg
  rthermg = surface_emissivity * intensity_g
  
  Es = surface_emissivity
  
  ! total thermal radiance contribution at band 3.7 micron
  rtherm37 = rthermc + rthermg


end subroutine toa_radiance37_cox_munk

subroutine calculate_new_Tc(platform_name, Tc, Tg, galbedo, wlen, tau, re, &
  lib_taus, lib_res, sph_albedo, down_flux_sensor, up_flux_sensor, &
  cloud_emiss, surface_emiss, newTc, PRN)
  !
  !  platform_name
  !  Tc
  !  Tg
  !  galbedo
  !  wlen
  !  tau
  !  re
  !  lib_taus
  !  lib_res
  !  sph_albedo
  !  down_flux_sensor
  !  up_flux_sensor
  !  cloud_emiss
  !  surface_emiss\
  !  newTc
  !  PRN
  !
  use modis_numerical_module
  use mod06_run_settings
  use core_arrays, only: Tc_low_for_delta, Tc_high_for_delta
  use planck_functions
  use nonscience_parameters
  implicit none
  
  character(*), intent(in) :: platform_name
  real, intent(in) :: Tc, Tg, tau, re, galbedo
  integer, intent(in) :: wlen
  real, dimension(:), intent(in) :: lib_taus, lib_res
  real, dimension(:,:,:), intent(in) :: sph_albedo, down_flux_sensor, &
    up_flux_sensor, cloud_emiss, surface_emiss
  real, intent(inout) :: newTc
  logical, intent(in) :: PRN
  
  integer  :: lowbound_tau, highbound_tau 
  real :: taubounds(2), rebounds(2)
  integer  :: lowbound_re, highbound_re
  real     :: iftaus(2), forint(2,2)  

  !  integer :: TOTAL_POINTS, total_taus
  !   real, dimension(:), allocatable :: taux

  real :: ec, eg
  real :: BTg, BTc, BTcr, temp_T_low, temp_T_high
  integer :: i
  real :: frefl1, trans, sfr1

  lowbound_tau = 0  ! WDR-UIV
  highbound_tau = 0  ! WDR-UIV
  lowbound_re = 0  ! WDR-UIV
  highbound_re = 0  ! WDR-UIV

  ! total_taus = size(lib_taus)
  !   TOTAL_POINTS = total_taus + 1

  !  allocate(taux(TOTAL_POINTS))

  !  taux(1) = 0.
  !  do i = 1, total_taus
  !     taux(i+1) = lib_taus(i)
  !  enddo
  
  call bisectionsearch(taux, tau, lowbound_tau, highbound_tau)
  call bisectionsearch(lib_res, re, lowbound_re, highbound_re)
  
  taubounds(1) = taux(lowbound_tau)
  taubounds(2) = taux(highbound_tau)
  rebounds(1) = lib_res(lowbound_re)
  rebounds(2) = lib_res(highbound_re)


  if (COX_MUNK) then 
  
    forint(1,1) = cloud_emiss(lowbound_tau, 2, lowbound_re)
    forint(1,2) = cloud_emiss(highbound_tau, 2, lowbound_re)
    forint(2,1) = cloud_emiss(lowbound_tau, 2, highbound_re)
    forint(2,2) = cloud_emiss(highbound_tau, 2, highbound_re)


    ec = bilinear_interpolation(  taubounds, rebounds, tau, re, forint, 0)
              
    forint(1,1) = surface_emiss(lowbound_tau, 2, lowbound_re)
    forint(1,2) = surface_emiss(highbound_tau, 2, lowbound_re)
    forint(2,1) = surface_emiss(lowbound_tau, 2, highbound_re)
    forint(2,2) = surface_emiss(highbound_tau, 2, highbound_re)

    eg = bilinear_interpolation( taubounds, rebounds, tau, re, forint, 0)

  else
    ! in Lambertian libraries we are missing
    lowbound_tau = lowbound_tau - 1
    highbound_tau = highbound_tau - 1

    if(highbound_re < 1) highbound_re = 1
    if(lowbound_re < 1) lowbound_re = 1
    if (highbound_tau < 1) highbound_tau = 1
    if (lowbound_tau < 1) lowbound_tau = 1

    if (lowbound_tau == 0) then 
      forint(1,1) = 0.
      forint(2,1) = 0.
    else
      forint(1,1) = sph_albedo(lowbound_tau, wlen, lowbound_re)
      forint(2,1) = sph_albedo(lowbound_tau, wlen, highbound_re)
    endif

    forint(2,2) = sph_albedo(highbound_tau, wlen, highbound_re)
    forint(1,2) = sph_albedo(highbound_tau, wlen, lowbound_re)


    sfr1 = bilinear_interpolation(  taubounds, rebounds, tau, re, forint, 0)
    
    if (lowbound_tau == 0) then 
      forint(1,1) = 1.
      forint(2,1) = 1.
    else
      forint(1,1) = down_flux_sensor(lowbound_tau, wlen, lowbound_re)
      forint(2,1) = down_flux_sensor(lowbound_tau, wlen, highbound_re)
    endif
    
    forint(2,2) = down_flux_sensor(highbound_tau, wlen, highbound_re)
    forint(1,2) = down_flux_sensor(highbound_tau, wlen, lowbound_re)
    
        
    trans = bilinear_interpolation(  taubounds, rebounds, tau, re, forint, 0)
    
    if (lowbound_tau == 0) then 
      forint(1,1) = 0.
      forint(2,1) = 0.
    else
      forint(1,1) = up_flux_sensor(lowbound_tau, wlen, lowbound_re)
      forint(2,1) = up_flux_sensor(lowbound_tau, wlen, highbound_re)
    endif
    
    forint(2,2) = up_flux_sensor(highbound_tau, wlen, highbound_re)
    forint(1,2) = up_flux_sensor(highbound_tau, wlen, lowbound_re)
        
    frefl1 = bilinear_interpolation(  taubounds, rebounds, tau, re, forint, 0)

    ec = 1.- trans - frefl1
    eg = trans*(1.-galbedo) / (1.-sfr1*galbedo)

  endif
  
  BTcr = modis_planck(platform_name, Tc, channel_11um, 1)
  BTg = modis_planck(platform_name, Tg, channel_11um, 1)

  BTc = (BTcr - eg*BTg) / ec
  if (BTc < 0.) then 
    newTc = fillvalue_real
    return
  endif

  newTc = modis_bright(platform_name, BTc, channel_11um, 1)

  ! recalculate Tc_low and Tc_high
  BTcr = modis_planck(platform_name, Tc_low_for_delta, channel_11um, 1) 
  BTc = (BTcr - eg*BTg) / ec
  if (BTc < 0.) then 
    newTc = fillvalue_real
    return
  endif
  temp_T_low = modis_bright(platform_name, BTc, channel_11um, 1)
        
  BTcr = modis_planck(platform_name, Tc_high_for_delta, channel_11um, 1)  
  BTc = (BTcr - eg*BTg) / ec
  if (BTc < 0.) then 
    newTc = fillvalue_real
    return
  endif
  temp_T_high = modis_bright(platform_name, BTc, channel_11um, 1)
  
  if (temp_T_high < newTc .or. temp_T_low > newTc) then
    newTc = fillvalue_real
    return
  endif
  ! overwrite these variables only when nothing broke down
  Tc_low_for_delta = temp_T_low
  Tc_high_for_delta = temp_T_high
    
end subroutine calculate_new_Tc



end module retrieval_prep_logic
