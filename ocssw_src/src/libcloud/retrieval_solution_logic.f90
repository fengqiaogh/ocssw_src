module retrieval_solution_logic

  implicit none

  integer, parameter :: TOP_NK_ASL_CONTOUR_ICE = 1, &
    TOP_NK_ASL_CONTOUR_WATER = 2
  
contains

  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine find_zero_crossings (residual, crossing_vector, crossing_num)
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  !   *** S. Platnick, 17-Oct-2005
  !   *** G. Wind 20-Oct-2005 -- removed the dynamic memory 
  !  allocation to speed things up a bit
  !
  !  residual
  !  crossing_vector
  !  crossing_num
  !
  implicit none

  REAL, intent(in)     :: residual(:)
  INTEGER, intent(out) :: crossing_num, crossing_vector(:)

  !   local variables
  ! INTEGER, allocatable, dimension(:) :: k
  integer :: k(18)

  INTEGER :: i,length

  length = size(residual)
  ! allocate (k(length))

  ! G.Wind 12.19.05: Replaced the where statement below with a do loop.  This 
  ! eliminates ambiguities arising from the use of nonconformable arrays. 
  k =   1
  !   where(residual<0) k=-1

  do i=1, length
    if ( residual(i) < 0 ) k(i) = -1
  end do  
  ! crossing index in vector is defined for the smaller index 
  !  between which the zero lies
  crossing_num=0; crossing_vector=0.
  Do i=1,length-1
    if( abs(k(i+1)-k(i)) > 1) then!{
      crossing_num = crossing_num + 1
      crossing_vector(crossing_num) = i
    end if
  End do
    
  ! deallocate (k)
    
end subroutine find_zero_crossings


  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine find_local_minima (residual, local_min_vector, local_min_num)
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  !   *** S. Platnick, 17-Oct-2005
  !   This give the local minima for the absolute residual, not the residual**2.
  !   10-18-05: modify so that minima at end points are not 
  !  considered local minima
  !    to the extent that they do not indicate a change in curvature.
  !  residual
  !  local_min_vector
  !  local_min_num

  implicit none

  REAL, intent(in)     :: residual(:)
  INTEGER, intent(out) :: local_min_num, local_min_vector(:)

  !   local variables
  INTEGER i,length

  length = size(residual)

  local_min_num=0; local_min_vector=0.

  if (length > 2) then!{
    Do i=2,length-1
      if( residual(i) < residual(i-1) .AND. residual(i) < residual(i+1) ) then!{
        local_min_num = local_min_num + 1
        local_min_vector(local_min_num) = i
      end if
    End do
  end if

end subroutine find_local_minima

  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine find_local_maxima (residual, local_max_vector, local_max_num)
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  !   *** S. Platnick, 17-Oct-2005
  !   This give the local maxima for the absolute residual, not the residual**2.
  !   10-18-05: modify so that minima at end points are not 
  !   considered local maxima
  !    to the extent that they do not indicate a change in curvature.
  !
  !  residual
  !  local_max_vector
  !  local_max_num
  !
  implicit none

  REAL, intent(in)     :: residual(:)
  INTEGER, intent(out) :: local_max_num, local_max_vector(:)

  !   local variables
  INTEGER i,length

  length = size(residual)

  local_max_num=0; local_max_vector=0.

  if (length > 2) then!{
    Do i=2,length-1
      if( residual(i) > residual(i-1) .AND. residual(i) > residual(i+1) ) then!{
        local_max_num = local_max_num + 1
        local_max_vector(local_max_num) = i
      end if
    End do
  end if

end subroutine find_local_maxima

  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine solution_re (re, residual, crossing_vector, crossing_num,    &
      local_min_vector, local_min_num, local_max_vector, local_max_num, &
      retrieval, quality )
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  use science_parameters 

  !  *** S. Platnick, 17-Oct-2005
  !      Using input from sep routines 'find_zero_crossings', 
  !          'find_local_minima', and 'find_local_maxima'
  !
  !  *** S. Platnick, 7-March-2006
  !
  !   Modification of collection 5 solution logic:
  !   1. process_zero_crossing_residual = .FALSE., i.e., re for 
  !      pixels with no zero crossings are set to fill
  !   2. solve_for_zero_crossing modified to use spline interpolation 
  !      and iterative root finder
  !   3. No longer need 'find_local_minima' and 'find_local_maxima' 
  !      routines, though left intact for now since
  !      they provide useful information that might still be needed sometime.
  !   4. No retrieval when crossing_num > 2

  implicit none

  integer, intent(in)    :: crossing_vector(:), crossing_num, &
    local_min_vector(:), local_min_num, local_max_vector(:), local_max_num
  real, intent(in)       :: residual(:), re(:)

  real, intent(out)      :: retrieval
  integer*1, intent(out) :: quality

  integer                :: residual_crossing_index, index1(1)
  real                   :: fillval_r4, resid_lo, resid_hi, slope_1, slope_2
  logical                :: process_zero_crossing_residual


  fillval_r4 = INVALID_ATTEMPTED_BUT_FAILED_RE; residual_crossing_index = 0

  ! ===================================
  ! SOLUTION LOGIC QUALITY ASSIGNMENTS: 
  ! ===================================
  ! Can be used for diagnostic purposes and setting retrieval 
  ! confidence QA in core science (not presently being
  !  done as of this writing)
  !
  !   quality =  -6  => spline interpolation, iteration /< iternum-5
  !   quality =  -5  => crossing_num < 0 (should never occur, but 
  !         included for completeness)
  !   quality =  -4  => crossing_num = 0 and process_zero_crossing_residual
  !        =.FALSE.
  !   quality =  -3  => crossing_num = 0 and process_zero_crossing_residual
  !        =.TRUE.
  !   quality =  -2  => crossing_num > 2
  !   quality =  -1  => DEFAULT: this quality index should never occur, 
  !         and if it does something went awry
  !   quality =   0  => occur for truncated residual vector having fill value(s)
  !   quality =   1  => linear interpolation w/maxval(abs(residual)) > 
  !         resid_max_for_spline
  !   quality =   2  => linear interpolation when the root lies in the 
  !         last interval of the residual vector
  !   quality =   3  => spline interpolation (or no interpolation), 
  !         and iteration < iternum-5
  !
  ! If ever used in MOD06 for assigning confidence QA then possible 
  !         negative assignments are:
  !    -6->3 or 2;   -5->0;   -4->0;   -3->1 or 0;   
  !     -2->1 or 0 depending on solution logic
  
  ! --------------------------------------------------------------
  ! .TRUE. corresponds to the logic that had apparently been in 
  !  place since the beginning of the c5 deliveries
  ! 2 March 06: Changed to .FALSE., i.e., for a non-zero crossing 
  !  residual, return a fill value. Do not return with
  !   the library radius corresponding to smallest residual as the solution.
  ! -------------------------------------------------------------

  process_zero_crossing_residual = .false.   

  ! -------------------------
  ! Default retrieval output:
  ! -------------------------
  retrieval = fillval_r4; quality = -1

  ! ------------------------------------------
  ! Check for a proper value for crossing_num:
  ! ------------------------------------------

  if (crossing_num < 0) then
    quality = -5
    return
  end if
  ! ---------------------------------------------------------------------------
  ! If there are more than 2 roots for the residual vector return a fill value:
  ! ---------------------------------------------------------------------------

  if (crossing_num > 2) then
    quality = -2
    return
  end if

  ! -----------------------------------------------------------------
  ! The calling routine in C5 truncates the residual array to get 
  !  rid of contiguous fill values starting with the largest radii, and
  !  then passes the useful max radii to the re solution routines. 
  !  However, the calling routine does not look for remaining fill values
  !  that are elsewhere in the residual vector but with non-fill 
  !  neighbors. If such a residual vector remains, then return with a
  !  fill value for effective radius.
  !
  !  However, I can't determine if resdiual vector is initialized 
  !  to fillvalue_real !?
  !
  ! ----------------------------------------------------------------
  
  !  if (min(residual) == fillvalue_real) then
  !     quality = 0
  !     return
  !  end if
  
  ! ------------------------------------------------------------------------
  ! If there are NO roots for the residual vector (i.e., no zero crossings):
  ! ------------------------------------------------------------------------
  IF (crossing_num == 0) THEN
    if (process_zero_crossing_residual) then

      index1 = MINLOC(abs(residual))
      retrieval = re(index1(1))
      quality = -3
    else
      quality = -4
    end if
    return
  END IF
  
  ! --------------------------------------------------
  ! If there are 1 or 2 roots for the residual vector:
  ! --------------------------------------------------

  IF (crossing_num == 1 .OR. crossing_num == 2) THEN

  ! If there are 2 zero crossings for the residual vector 
  !  then choose the crossing having the negative slope
  !  regardless of the cloud phase. A derivative of residual 
  !  w.r.t. re < 0 corresponds to the physical situation
  !  where single scattering albedo decreases with re (and no 
  !  significant change in effective radius or extinction).

  If (crossing_num == 1) Then
    residual_crossing_index = crossing_vector(crossing_num)

  Else
    slope_1 = residual(crossing_vector(1)+1) - residual(crossing_vector(1))
    slope_2 = residual(crossing_vector(2)+1) - residual(crossing_vector(2))
    !   resid_lo = 0.5 * ( residual(crossing_vector(1)) + &
    !      residual(crossing_vector(1) + 1) )
    !   resid_hi = 0.5 * ( residual(crossing_vector(2)) + &
    !     residual(crossing_vector(2) + 1) )

    !  For crossing_num =2, both slopes can't be of the same sign. 
    !  If they are something's not correct, but choose larger residual.
    if ( slope_1 <= 0. ) then
      residual_crossing_index = crossing_vector(1)
    end if
    if ( slope_2 <= 0. ) then 
      residual_crossing_index = crossing_vector(2)
    end if


  End If
  
  ! it could potentially happen that slopes are > 0, then this 
  !  would be a segfault because residual_crossing_index 
  ! defaults to 0. --- G.Wind 12.15.11
  if (residual_crossing_index > 0) &
    call solve_for_zero_crossing (re, residual, residual_crossing_index, &
      fillval_r4, local_min_vector, local_min_num, local_max_vector, &
      local_max_num, retrieval, quality)
  END IF    

end subroutine solution_re

  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine solve_for_zero_crossing (re, residual, &
    residual_crossing_index, fillval_r4, &
    local_min_vector, local_min_num, local_max_vector, local_max_num, &
    retrieval, quality)
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  use GeneralAuxType
  use spline_module

  !  *** S. Platnick, 7-March-2006
  !
  !   Modification of collection 5 solution logic:
  !
  !   1. Use natural spline interpolation instead of quadratic(s) fits.
  !      The piece-wise quadratic approach was producing discontinuities 
  !      in the re histograms. Epecially for ice clouds where
  !      the residual vector often has changes in curvature which exacerbates 
  !      the piece-wise approach. The advantage of the 
  !      quadratic fits is that they are analytic and computationally 
  !      fast; the spline fit must use an iterative procedure for
  !      for finding the root of the residual vector.
  !
  !   2. Also, see comments in solution_re routine.

  implicit none

  real,    intent (in)   :: re(:), residual(:), fillval_r4
  integer, intent (in)   :: residual_crossing_index,  &
    local_min_vector(:), local_min_num, local_max_vector(:), local_max_num

  real,      intent(inout) :: retrieval
  integer*1, intent(inout) :: quality
   
  integer i, icount, iternum
  real    y2(size(re)), z, delta_resid, delta_resid_fraction, &
    resid_min_allowable, resid_max_for_spline
  real    yl(2), xl(2), re_lower, re_upper, re_iter
   
  character*15 numerical_type

  retrieval=fillval_r4

  !  -------------------------------------------------------------
  !  If the root is close enough to a library re, don't bother 
  !   interpolating for the solution
  !   (also, spline bisection search was found to fail for a residual 
  !   index of zero which did occur for test granules)
  !  --------------------------------------------------------------

  resid_min_allowable = 0.0001
  if (abs(residual(residual_crossing_index)) < resid_min_allowable) then
    retrieval = re(residual_crossing_index)
    quality = 3
    return
  end if
  if (abs(residual(residual_crossing_index+1)) < resid_min_allowable) then
    retrieval = re(residual_crossing_index+1)
    quality = 3
    return
  end if

  !  -----------------------------------------------------------------
  !  If bounding library points are so near to one another as to be 
  !    representationally equal, no need to interpolate
  !  ----------------------------------------------------------------

  if (real_s_equal(residual(residual_crossing_index), &
    residual(residual_crossing_index+1))) then
    retrieval = re(residual_crossing_index)
    quality = 3
    return
  end if

  !  -----------------------------------------------------------------
  !  Use either a linear or spline interpolation to determine the 
  !  residual vector root:
  !  ------------------------------------------------------------------
  numerical_type = 'spline'

  !  The following 2 Do Loops set linear interpolation for residual 
  !  roots near a local minimum or maximum. Commented out
  !   to use spline interpolation in these instances. Use of linear 
  !   interpoltion gave some re histogram discontinuities in
  !   test granules, at least for ice clouds where local min and max 
  !   seemed to occur more often.

  !   Do i=1,local_max_num
  !     if (residual_crossing_index == local_max_vector(i) .OR. &
  !       residual_crossing_index+1 == local_max_vector(i)) then
  !       numerical_type = 'linear'
  !       exit
  !     end if
  !   End do
  !   IF (numerical_type /= 'linear') THEN
  !     Do i=1,local_min_num
  !       if (residual_crossing_index == local_min_vector(i) .OR. &
  !         residual_crossing_index+1 == local_min_vector(i)) then
  !         numerical_type = 'linear'
  !         exit
  !       end if
  !     End do
  !   End If


  !  Use linear interpolation if the root is before end of re vector. 
  !   Especially for ice clouds, the spline appears to give undue 
  !   curvature to fit the large distance from 60 to 90 microns. For 
  !   water clouds, the accuracy in establishing the root between
  !   28 and 30 microns isn't relevant. Use of this logic did not 
  !   appear to cause any undue features in the re histogram.

  ! The following innocuous check on 'linear' remains for historical 
  !  reasons in case the block immediately above ever gains favor
 IF (numerical_type /= 'linear') THEN
   if (residual_crossing_index == size(re)-1) then
     numerical_type = 'linear'
     quality = 2
   end if
  End If

  !  Use linear interpolation if there are large discontinuities in 
  !  the residual which would give rise to spline oscillations
  !
  !   Example print output for 3.7 um residual, MYD Katrina granule: 
  !   iterational doesn't converge because the spline does not provide
  !   a root between 30 and 35 um points due to large oscillation set 
  !   up by large residual at larger radii.
  !
  !   *** number of spline iterations: icount = 30  residual_crossing_index= 6
  !   re:        5.0   10.0   15.0   20.0   25.0   30.0   35.0   
  !   40.0   45.0   50.0   55.0   60.0   90.0
  !   residual: 3.9147  2.1260  0.6555  0.0979  -0.0182  -0.0686  0.0729  
  !    0.4781  2.1382  10000.0  383.20  10000.0  10000.0

  resid_max_for_spline = 100.
  IF (numerical_type /= 'linear') THEN
    if ( maxval(abs(residual)) > resid_max_for_spline) then
      numerical_type = 'linear'
      quality = 1
    end if
  End If

  !  ------------------------------------------------------------------
  !  Implement interpolation scheme and solve for residual vector root:
  !  ------------------------------------------------------------------

  IF (numerical_type == 'linear') THEN
    yl(1) = residual(residual_crossing_index); yl(2) = &
      residual(residual_crossing_index+1)
    if (sign(1.,yl(1))*sign(1.,yl(2)) <= 0.) then
      call linear_interpolate_for_root (yl, &
        re(residual_crossing_index:residual_crossing_index+1), retrieval)
    end if
  END IF

  IF (numerical_type == 'spline') THEN
    call spline (size(re), re, residual, y2)

    iternum = 30; icount = 1; delta_resid_fraction=0.001

    delta_resid = abs( delta_resid_fraction * &
      max(residual(residual_crossing_index), &
      residual(residual_crossing_index+1)) )
    if(delta_resid < resid_min_allowable) delta_resid = resid_min_allowable
    re_lower = re(residual_crossing_index)
    re_upper = re(residual_crossing_index+1)
    re_iter = 0.5*(re_lower + re_upper)
  
    !  The interpolation finds the y (residual) for the input x (re), so 
    !  they iterate till the y they find is close to 0.  Now, if the 
    !  relation is reversed, and y is re and x is residual, you could do 
    !  the spline 1 time for x = 0.  All the min, max points could be used to 
    !  isolate the correct segment so a function wuld be analyzed.`
    DO WHILE (icount < iternum)
      z= splint (size(re), re_iter, re, residual, y2)

      if(abs(z) < delta_resid) then
        retrieval = re_iter
        exit
      end if
    
      if(sign(1.,z)*sign(1.,residual(residual_crossing_index)) > 0.) then
        re_lower = re_iter
      else
        re_upper = re_iter
      end if  
      re_iter = 0.5*(re_lower + re_upper)
      icount = icount + 1
    END DO

    ! If consistently near the max iteration number, it is possible a 
    !  higher max iteration number might be more suitable.
    if(icount < iternum-5) then
      quality =  3
    else
      quality = -6
    end if

  END IF

end subroutine solve_for_zero_crossing


  logical function is_water_phase(library_radii)

  use libraryarrays

  real(single), intent(in) :: library_radii(:)

  ! if it's not ice phase, it's water
  is_water_phase = (size(library_radii) .ne. size(ice_radii)) 

  ! NOTE:  if at some point ice & water libs have equal sizes, 
  !  this needs reconsideration

  end function is_water_phase

subroutine ray_corr_nearest(refl_source, As, iw, tau, re, Pc,  &
    sfr, fti1, fti0, fluxup_solar, fluxup_sensor, &
    theta, theta0, phi, location, crefl)
  !    don't know what this does specifically
  !  refl_source
  !  As
  !  iw
  !  tau
  !  re
  !  Pc
  !  sfr
  !  fti1
  !  fti0
  !  fluxup_solar
  !  fluxup_sensor
  !  theta
  !  theta0
  !  phi
  !  location  I  radius location
  !  crefl  O returned reflectance?

  use science_parameters
  use libraryarrays
  use mod06_run_settings
  use modis_numerical_module

  real, intent(in) :: tau, re, Pc, sfr, fti1, fti0, theta, theta0, phi, &
    As, refl_source
  integer, dimension(2), intent(in) :: location
  integer, intent(in) :: iw
  real, dimension(:,:,:,:), intent(in) :: fluxup_solar, fluxup_sensor
  real, intent(inout) :: crefl
  
  real, parameter :: Ps = 1013.
  real, parameter :: Cm = 0.84
  real, parameter :: taur0(2) = (/0.044, 0.025/)

  real :: taur, mu0, mu, x, pray, refl_s
  real :: fti1_as, fti0_as
  real :: angles(2), functionvalues(2)
  real :: fluxup_solar_tau, fluxup_sensor_tau
  integer :: lowbound, highbound

  lowbound = 0  ! WDR-UIV
  highbound = 0  ! WDR-UIV

  taur = taur0(iw) * Pc / Ps
  mu0 = cos(d2r*theta0)
  mu = cos(d2r*theta)

  x = -mu0*mu + sqrt((1.-mu0**2)*(1.-mu**2))*cos(d2r*phi)
    pray = 3.*(1.+x**2)/4.0

  !
  !  Interpolation for solar zenith angles and scaled optical thickness
  !
  call bisectionsearch(library_fluxsolarzenith, mu0, lowbound, highbound)
  angles(1) = library_fluxsolarzenith(lowbound)
  angles(2) = library_fluxsolarzenith(highbound)

  if (location(1) == 1) then 
    functionvalues = 0.
  else 
    functionvalues(1) = fluxup_solar(lowbound,location(1)-1,1,location(2))
    functionvalues(2) = fluxup_solar(highbound,location(1)-1,1,location(2))
  endif

  fluxup_solar_tau = linearinterpolation(angles,functionvalues,mu0)

  !
  !  Interpolation for sensor zenith angles and scaled optical thickness
  !
  call bisectionsearch(library_fluxsensorzenith, mu, lowbound, highbound)
  angles(1) = library_fluxsensorzenith(lowbound)
  angles(2) = library_fluxsensorzenith(highbound)

  if (location(1) == 1) then 
    functionvalues = 0.
  else 
    functionvalues(1) = fluxup_sensor(lowbound,location(1)-1,1,location(2))
    functionvalues(2) = fluxup_sensor(lowbound,location(1)-1,1,location(2))
  endif
  fluxup_sensor_tau = linearinterpolation(angles,functionvalues,mu)

  fti1_as = fluxup_sensor_tau
  fti0_as = fluxup_solar_tau

  if (As > 0.) then 
    fti0_as = fti0_as + fti0*As*(1-sfr)/(1 - sfr*As)
    fti1_as = fti1_as + fti1*As*(1-sfr)/(1 - sfr*As)
  endif

  refl_s = taur*pray/(4.*mu0*mu) + 0.5*taur* fti1_as * exp(-taur/mu )/mu0 + &
    0.5*taur* fti0_as  * exp(-taur/mu0)/mu
  crefl = (refl_source - refl_s)*exp(Cm*taur*(1./mu0 + 1./mu))

end subroutine ray_corr_nearest

subroutine solveretrieval_nearest(xx_pt,yy_pt,Ram, Rbm, twobands, &
    radii, tau, re, lib_dist, phase_liquid, Ram_corr,quality_in,CH37_IDX, &
    CTopT,CH37_NUM, platFormName)
  !  xx_pt, IN  x, loc of point in the chunk
  !  yy_pt  IN  y loc of point in the chunk
  !  Ram  IN  non-abs reflectance
  !  Rbm   IN  absorbing band reflectance
  !  twobands  IN  size 2 indicies of non-abs and absorbing indiciew
  !  radii   IN  array of table re values
  !  tau    OUT  retrieved optical thickness
  !  re   OUT  retrieved final radius
  !  lib_dist  OUT  I believe a uncertainty measure (root sum square) for the 
  !       derived re / tau?
  !  phase_liquid   IN  switch for liquid phase if true
  !  Ram_corr    
  !  quality_in
  !  --- extra args? ---
  !  CH37_IDX
  !  CTopT,CH37_NUM
  !  platFormName
  !
  !  alternate solution exterior, just the boundary (for water, re=4, 
  !  re=30 and tau = 158.8)
  !  for ice re=5,re=60 and tau = 158.8 (top, bottom and rightmost

  use GeneralAuxType
  use core_arrays
  use science_parameters
  use nonscience_parameters
  use libraryinterpolates
  use libraryarrays
  use planck_functions

  integer, intent(in) ::xx_pt,yy_pt
  real, intent(in) :: Ram, Rbm
  real, dimension(:), intent(in) :: radii
  integer, dimension(:), intent(in) :: twobands
  real, intent(inout) :: tau, re, lib_dist, Ram_corr
  logical, intent(in) :: phase_liquid
    integer, intent(in) :: quality_in
    CHARACTER(len=*),INTENT(IN),OPTIONAL::platFormName
  INTEGER,INTENT(IN),OPTIONAL::CH37_NUM,CH37_IDX
  REAL, INTENT(IN), OPTIONAL::CTopT
  
  REAL :: Pc, theta, theta0, phi
  real::CTopT_lcl, intCTop,intSurf
  integer :: size_tau, size_re

  size_tau = number_taus + 1
  size_re = size(reflibA, 2)
  
  Pc = cloud_top_pressure(xx_pt,yy_pt)
  theta = sensor_zenith_angle(xx_pt,yy_pt)
  theta0= solar_zenith_angle(xx_pt,yy_pt)
  phi=relative_azimuth_angle(xx_pt,yy_pt)
  !
  !  WDR have the CTopT go to a local variable and handle presence then
  IF(PRESENT(CTopT)) THEN
    CTopT_lcl = CTopT
  else
    CTopT_lcl = 0.
  endif
  !  3.7 not in PACE, so return later
  IF(PRESENT(CH37_IDX) .AND. PRESENT(CTopT) .and. PRESENT(platFormName) &
    .and. PRESENT(CH37_NUM) .and. CTopT_lcl > 0.0)THEN     
    intCTop= modis_planck(platFormName, CTopT_lcl, CH37_NUM, 1)
    intSurf = modis_planck(platFormName, surface_temperature(xx_pt,yy_pt), &
      CH37_NUM, 1)
       
    if( (.NOT. COX_MUNK))THEN
      if(phase_liquid)then
        call calc37RadianceLibLamb(intCTop,intSurf, &
          thermal_correction_oneway(1), &
          thermal_correction_twoway(1), theta0, &
          spherical_albedo_water(:,CH37_IDX,:), &
          int_fluxdownwater_sensor(:,CH37_IDX,:), &
          int_fluxdownwater_solar(:,CH37_IDX,:), &
          int_fluxupwater_sensor(:,CH37_IDX,:))
      else
        call calc37RadianceLibLamb(intCTop,intSurf,  &
          thermal_correction_oneway(1), thermal_correction_twoway(1), theta0, &
          spherical_albedo_ice(:,CH37_IDX,:), &
          int_fluxdownice_sensor(:,CH37_IDX,:), &
          int_fluxdownice_solar(:,CH37_IDX,:), &
          int_fluxupice_sensor(:,CH37_IDX,:))
      endif
            
    elseif(COX_MUNK)THEN
      if(phase_liquid)then
        call calc37RadianceLibCM(intCTop,intSurf, &
          thermal_correction_oneway(1), &
          thermal_correction_twoway(1), theta0, &
          int_cloud_emissivity_water(:,1,:), &
          int_surface_emissivity_water(:,1,:))
      else
        call calc37RadianceLibCM(intCTop,intSurf,  &
          thermal_correction_oneway(1), thermal_correction_twoway(1), &
          theta0, int_cloud_emissivity_ice(:,1,:), &
          int_surface_emissivity_ice(:,1,:))
      endif
    endif
  ENDIF   ! End 3.7 related logic
 
  if(quality_in == -4 .OR. quality_in == -3)then
    CALL ASL_boundary(Ram, Rbm, twobands, radii, tau, re, lib_dist, Pc, & 
      theta, theta0, phi, phase_liquid, Ram_corr)                   
  else
    CALL ASL_interior(Ram, Rbm, twobands, radii, tau, re, lib_dist, Pc, & 
      theta, theta0, phi, phase_liquid, Ram_corr)
  endif
               
  return

end subroutine solveretrieval_nearest

subroutine calc37RadianceLibLamb(intensity,intensity_g, &
    thermal_trans_1way, thermal_trans_2way, solar_zenith, &
    sfr,fti1,fti0,fri1)

  use GeneralAuxType
  use science_parameters
  use libraryarrays, only : reflibB,rad37lib

  implicit none

  real, intent(in)     ::   intensity,intensity_g,solar_zenith

  real, intent(in)      :: thermal_trans_2way,thermal_trans_1way
  real(single), intent(in)  :: sfr(:,:),  fti1(:,:), fti0(:,:), fri1(:,:)
!  rfi(:,:)

  integer              :: total_taus, radii, wavelengths,i,j
  real                 :: absorbing_albedo
  real                 :: tc, local_trans_1way, local_trans_2way

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


  !populate rad37lib here, in reflectance units
  absorbing_albedo = reflibB(1,1)  
  
  rad37lib(1,:) = (1.0 - absorbing_albedo) * intensity_g 
  rad37lib(2:,:) = (1.0- fti1(:,:) - fri1(:,:)) * intensity + &
    fti1(:,:) * (1.0 - absorbing_albedo) * intensity_g / &
    (1.0 - absorbing_albedo * sfr(:,:))
  rad37lib(:,:) = rad37lib(:,:) * local_trans_1way * pi/ &
    ( cos(solar_zenith*d2r)*solar_constant_37)                                 
        rad37lib(:,:) = rad37lib(:,:) + local_trans_2way * reflibB(:,:)
          
end subroutine calc37RadianceLibLamb

subroutine calc37RadianceLibCM(intensity,intensity_g,     &
    thermal_trans_1way, thermal_trans_2way, solar_zenith,  &
    cl_emis,sf_emis)

  use GeneralAuxType
  use science_parameters
  use libraryarrays, only : reflibB,rad37lib

  implicit none

  real, intent(in) :: intensity,intensity_g,solar_zenith
  real, intent(in)      :: thermal_trans_2way,thermal_trans_1way
  real(single), intent(in)     ::  cl_emis(:,:), sf_emis(:,:)
                            
  !   real,intent(out)::rad037lib(:,:),rtherm037lib(:,:)
                           

  integer              :: total_taus, radii, wavelengths,i,j
  real                 :: absorbing_albedo
  real                 :: tc, local_trans_1way, local_trans_2way

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

  !populate rad37lib here, in reflectance units

  rad37lib(:,:) = cl_emis(:,:) * intensity + sf_emis(:,:) * intensity_g
  rad37lib(:,:) = rad37lib(:,:) * local_trans_1way * pi/ &
    ( cos(solar_zenith*d2r)*solar_constant_37)
  rad37lib(:,:)  = rad37lib(:,:) + local_trans_2way * reflibB(:,:)
                      
end subroutine calc37RadianceLibCM

subroutine ASL_interior(Ram, Rbm_in, bands, radii, tau, re, lib_dist, Pc, & 
                   theta, theta0, phi, phase_liquid, Ram_corr)
  !
  !  For solutions inside the LUT with multiple solutions - derive one and
  !  essentialy an uncertainty in the form of a distance - This is for
  !alternate solution interior                   
  !
  !  Ram
  !  Rbm_in
  !  bands
  !  radii
  !  tau
  !  re
  !  lib_dist
  !  Pc
  !  theta
  !  theta0
  !  phi
  !  phase_liquid
  !  Ram_corr
  !
  use GeneralAuxType
  use libraryarrays
  use libraryinterpolates
  use nonscience_parameters
  use mod06_run_settings
  use science_parameters

  real, intent(in) :: Ram, Rbm_in, Pc, theta, theta0, phi
  real, dimension(:), intent(in) :: radii
  integer, dimension(:), intent(in) :: bands
  real, intent(inout) :: tau, re, lib_dist, Ram_corr
  logical, intent(in) :: phase_liquid

  
  real, dimension(:,:), allocatable ::  rel_residual,reflibBL
  integer :: size_tau, size_re, loc(2), rel_loc(2)
  logical :: re_altered
  
  real :: Ramc, crefl,Rbm
  integer :: iter, max_iter

  real :: albedoA, albedoB,normConst
  real :: sfr, fti1, fti0

  size_tau = number_taus + 1
  size_re = size(reflibA, 2)
  allocate(rel_residual(size_tau, size_re),reflibBL(size_tau, size_re))
  
  reflibBL(:,:) = reflibB(:,:)
  Rbm = Rbm_in
      
  if( bands(2) == band_0370)then 
    normConst = pi/(cos(theta0*d2r)*solar_constant_37)
    Rbm = Rbm_in * normConst
    reflibBL(:,:) = rad37lib(:,:)
  endif

  albedoA = reflibA(1,1)
  albedoB = reflibBL(1,1)


#if ASTER_INST
  if (bands(1) == band_0065 .or. bands(1) == band_0086) then 
#else
  if (bands(1) == band_0065) then 
#endif  
    max_iter = 3
  else 
    max_iter = 1
  endif
  
  Ramc = Ram
  iter = 1
  
  do while (iter <= max_iter) 

    ! the point has no chance, it is darker than the surface albedo
    if (Ramc < albedoA) then
      tau = fillvalue_real
      re  = fillvalue_real
      lib_dist = MAX_COST_FUNCTION
      Ram_corr = Ramc
      deallocate (rel_residual, reflibBL )  ! WDR add the reflibBL deallocate
      return
    endif
    ! collect the value of cost function
    rel_residual = sqrt ( (reflibA - Ramc)**2 + (reflibBL-Rbm)**2 ) / &
      sqrt (Ramc**2 + Rbm**2)
    loc = minloc(rel_residual)  
    lib_dist = rel_residual(loc(1), loc(2))
    re = radii(loc(2))
  
    ! don't want any 2's and such 
    if (re < 4. ) re = 4.
    ! the library taus don't have 0., so we need to explicitly count for 
    !  it. But we can't put in a 0. straight because
    !  L3 code can't handle an explicit 0. value. 
    if (loc(1) == 1) then 
      tau = 0.01
    else
      tau = library_taus(loc(1)-1)
    endif
  
    ! call rayleigh here
#if ASTER_INST
    if (bands(1) == band_0065 .or. bands(1) == band_0086) then 
#else
    if (bands(1) == band_0065) then 
#endif  
      if (phase_liquid) then 
        if (loc(1) == 1) then 
          sfr = 0.
          fti1 = 1.
          fti0 = 1.
        else
          sfr = spherical_albedo_water(loc(1)-1,bands(1),loc(2))
          fti1 = int_fluxdownwater_sensor(loc(1)-1,bands(1),loc(2))
          fti0 = int_fluxdownwater_solar(loc(1)-1,bands(1),loc(2))
        endif
          
        call ray_corr_nearest (Ram, albedoA, bands(1), tau, re, Pc, &
          sfr, fti1, fti0, flux_up_water_solar, flux_up_water_sensor, &
          theta, theta0, phi, loc, crefl)
      else 
        if (loc(1) == 1) then 
          sfr = 0.
          fti1 = 1.
          fti0 = 1.
        else
          sfr = spherical_albedo_ice(loc(1)-1,bands(1),loc(2))
          fti1 = int_fluxdownice_sensor(loc(1)-1,bands(1),loc(2))
          fti0 = int_fluxdownice_solar(loc(1)-1,bands(1),loc(2))
        endif
        
        call ray_corr_nearest (Ram, albedoA, bands(1), tau, re, Pc, &
          sfr, fti1, fti0, flux_up_ice_solar, flux_up_ice_sensor, &
          theta, theta0, phi, loc, crefl)
      endif             
      Ramc = crefl
  
    endif
    
    iter = iter + 1

  end do
  
  Ram_corr = Ramc

  if (Ramc < albedoA) then
    tau = fillvalue_real
    re  = fillvalue_real
    lib_dist = MAX_COST_FUNCTION
    deallocate (rel_residual, reflibBL )  !  WDR add this to complete deallocs
    return
  endif
  ! collect the value of cost function
  rel_residual = sqrt ( (reflibA - Ramc)**2 + (reflibBL-Rbm)**2 ) / &
    sqrt (Ramc**2 + Rbm**2)
  loc = minloc(rel_residual)  
  lib_dist = rel_residual(loc(1), loc(2)) 
  
  ! do the alternate retrieval  
  re = radii(loc(2))
  ! don't want any 2's and such 
  if (re < 4. ) re = 4.
  ! the library taus don't have 0., so we need to explicitly count 
  !  for it. But we can't put in a 0. straight because
  ! L3 code can't handle an explicit 0. value. 
  if (loc(1) == 1) then 
    tau = 0.01
  else
    tau = library_taus(loc(1)-1)
  endif

  deallocate (rel_residual,reflibBL)

end subroutine ASL_interior


subroutine ASL_boundary(Ram, Rbm_in, bands, radii, tau, re, lib_dist, Pc, & 
    theta, theta0, phi, phase_liquid, Ram_corr)

  !  For solutions outside the LUT - derive one and 
  !  essentialy an uncertainty in the form of a distance - This is for 
  !alternate solution exterior, just the boundary (for water, re=4, 
  !  re=30 and tau = 158.8)
  !  for ice re=5,re=60 and tau = 158.8 (top, bottom and rightmost)
  !
  !  Ram  IN  non-abs reflectance
  !  Rbm_in  IN  absorbing band reflectance
  !  bands  IN  size 2 indicies of non-abs and absorbing indicies
  !  radii  IN  array of table re values
  !  tau  IN/OUT  retrieved optical thickness
  !  re  IN/OUT  retrieved radius
  !  lib_dist  a distance outside the LUT, a measure of amount off
  !  Pc  IN  cloud top pressure
  !  theta  IN  I believe the sat zen angle
  ! theta0  IN  I believe the sol zen angle
  ! phi   IN  I believe the delta azimuth angle
  ! phase_liquid  IN  switch for liquid phase if true
  !  Ram_corr  OUT  may be abs band refl

  use GeneralAuxType
  use libraryarrays
  use libraryinterpolates
  use nonscience_parameters
  use mod06_run_settings
  use science_parameters

  real, intent(in) :: Ram, Rbm_in, Pc, theta, theta0, phi
  real, dimension(:), intent(in) :: radii
  integer, dimension(:), intent(in) :: bands
  real, intent(inout) :: tau, re, lib_dist, Ram_corr
  logical, intent(in) :: phase_liquid
  
  real, dimension(:,:), allocatable :: temp_refl,reflibBL
  real,dimension(:),allocatable::rel_residual
  integer :: size_tau, size_re, loc, rel_loc(1),ire_ct_1, &
    temp_array_size,s_rt,iL,iR
  logical :: re_altered
  
  real :: Ramc, crefl, Rbm
  integer :: iter, max_iter,iiiloc

  real :: albedoA, albedoB,normConst
  real :: sfr, fti1, fti0

  size_tau = number_taus + 1
  size_re = size(reflibA, 2)

  Ramc = Ram
  Rbm = Rbm_in

  allocate(reflibBL(size_tau, size_re))
  
  reflibBL(:,:) = reflibB(:,:)
    
  if( bands(2) == band_0370)then 
    normConst = pi/(cos(theta0*d2r)*solar_constant_37)
    Rbm = Rbm_in * normConst
    reflibBL(:,:) = rad37lib(:,:)
  endif
  
  !  set the top contour (CT 1) for ice =1 (re =5) for water = 2(re = 4)
  !  can do it in mod06_run_settings

  ire_ct_1 = TOP_NK_ASL_CONTOUR_ICE 
  if(phase_liquid)ire_ct_1 = TOP_NK_ASL_CONTOUR_WATER
  
  !  reflibA is the modified copied non-abs reflectance table made back in 
  !  vis non absorbing work and reflibBL is the abs table
  albedoA = reflibA(1,ire_ct_1)
  albedoB = reflibBL(1,ire_ct_1)
  
  s_rt = size_tau+size_re
  temp_array_size = size_tau + s_rt - (ire_ct_1 + 1)
  
  allocate(rel_residual(temp_array_size))
  allocate(temp_refl(2,temp_array_size))
  ! EXPLAIN YOUR WORK!  
  temp_refl(1,1:size_tau) = reflibA(1:size_tau,ire_ct_1)
  
  temp_refl(1, size_tau+1 : s_rt-ire_ct_1) = &
    reflibA(size_tau,ire_ct_1+1:size_re)  !
  
  temp_refl(1, s_rt-ire_ct_1+1 : temp_array_size) = &
    reflibA(1:size_tau-1,size_re)

  !   assign the radiance+emission of the 3.7 LUT
  temp_refl(2, 1 : size_tau) = reflibBL(1:size_tau,ire_ct_1)
  temp_refl(2, size_tau+1 : s_rt-ire_ct_1) = &
    reflibBL(size_tau,ire_ct_1+1:size_re)  !
  temp_refl(2, s_rt-ire_ct_1+1 : temp_array_size) = &
    reflibBL(1:size_tau-1,size_re)
  
  !  abs refl is below sfc albedo: region I
  ! the point has no chance, it is darker than the surface albedo
  if (Ramc < albedoA  .OR.  & 
    ( Rbm < minval(reflibBL(size_tau,ire_ct_1:size_re)) .and. &
    Ramc > minval(reflibA(size_tau,ire_ct_1:size_re)) ) .OR. & 
    ( Rbm > maxval(reflibBL(size_tau,ire_ct_1:size_re)) .and. &
    Ramc > maxval(reflibA(size_tau,ire_ct_1:size_re)) )) then

    tau = fillvalue_real
    re  = fillvalue_real
    lib_dist = MAX_COST_FUNCTION
    Ram_corr = Ramc
    deallocate (rel_residual, temp_refl, reflibBL)
    return
    
  endif
  ! collect the value of cost function
  rel_residual = sqrt ( (temp_refl(1,:) - Ramc)**2 + &
    (temp_refl(2,:) -Rbm)**2 ) / sqrt (Ramc**2 + Rbm**2)
  rel_loc = minloc(rel_residual)  
  lib_dist = rel_residual(rel_loc(1))
  
  !  reclaculate lib_dist with both x and y in terms of reflectances, 
  !  if band(2)=3.7

  ! do the alternate retrieval  
  if(rel_loc(1) <= size_tau)then
    ! contour 1
    re = radii(ire_ct_1)
    loc = rel_loc(1)
       
    if(bands(1) == band_0065)then
      if(phase_liquid) then       
        Ramc = rayleigh_liq(ire_ct_1)
      else
        Ramc = rayleigh_ice(ire_ct_1)
      endif
            
      CALL calcDistanceAndMinLoc(Ramc,Rbm,temp_refl(1,1:size_tau), &
        temp_refl(2,1:size_tau),lib_dist,loc)
      Ram_corr = Ramc
    endif
             
    if (loc == 1) then 
      tau = 0.01
    else
      tau = library_taus(loc-1)
    endif
      
    if(Ramc > maxval(reflibA(size_tau,ire_ct_1:size_re)))then 
      if(bands(1) == band_0163  .AND. bands(2) == band_0213 )then
        tau = fillvalue_real
      else
        if((Rbm - reflibBL(size_tau,ire_ct_1)) < 0.0)then
          do iiiloc = ire_ct_1+1,size_re
            iL = iiiloc-1
            iR = iiiloc
            if((Rbm - reflibBL(size_tau,iiiloc)) >=0.0)exit
          enddo     
          tau = MAX_TAU_RTRIEVED
          re=radii(iL)

          if(abs(reflibBL(size_tau,iR)-reflibBL(size_tau,iL)) > 0.000001)then
            re = (Rbm-reflibBL(size_tau,iL))*(radii(iR) - &
              radii(iL))/(reflibBL(size_tau,iR)-reflibBL(size_tau,iL)) + &
              radii(iL)
            if(re > radii(iR) .OR. re < radii(iL))re = radii(iL)
          endif
        endif
      endif         
    endif
      
  elseif(rel_loc(1) >  size_tau .and. rel_loc(1) <= s_rt - ire_ct_1)then
    !contour 4
    loc = rel_loc(1) - size_tau + 1
    re = radii(loc)

    if(loc==1 .OR. &
      ((Rbm - reflibBL(size_tau,loc)) < 0.0 .AND. loc <= size_re))then
      do iiiloc = loc , size_re
        if((Rbm - reflibBL(size_tau,iiiloc)) >=0.0)exit
        loc = iiiloc
      enddo
          
      iL = loc
      iR = loc+1
    elseif((Rbm - reflibBL(size_tau,loc)) >= 0.0 .AND. loc <= size_re)then    
      do iiiloc = loc-1 , 1,-1
        if((Rbm - reflibBL(size_tau,iiiloc)) <= 0.0)exit
        loc = iiiloc
      enddo
      iL = loc-1
      iR = loc
    endif
       
    if(IR <= size_re .and. IL >= 1)then
      re = radii(loc)
      if(abs(reflibBL(size_tau,iR)-reflibBL(size_tau,iL)) > 0.000001 &
        .and. abs(Rbm - reflibBL(size_tau,loc)) > 0.0)then
        re = radii(loc) + (Rbm-reflibBL(size_tau,loc))*(radii(iR) - &
          radii(iL))/(reflibBL(size_tau,iR)-reflibBL(size_tau,iL))  
        if(re > radii(iR) .OR. re < radii(iL))re = radii(loc)            
      endif
    endif
        
    if(bands(1) == band_0065)then
      if(phase_liquid)then       
        Ramc = rayleigh_liq(loc)
      else
        Ramc = rayleigh_ice(loc)
      endif
      Ram_corr = Ramc
    endif
        
    tau = fillvalue_real
    if(bands(1) == band_0163  .AND. bands(2) == band_0213 )then
      tau = fillvalue_real
    else
      tau = MAX_TAU_RTRIEVED
    endif
        
    if(Ramc > 0 .and. Ramc < maxval(reflibA(1:size_tau,size_re)) .and. &
      iR > size_re)THEN
      tau = fillvalue_real
      re = fillvalue_real
      lib_dist = MAX_COST_FUNCTION          
    endif
  else
    !contour 2
    re = radii(size_re)
    loc = rel_loc(1) - (s_rt-ire_ct_1) + 1

    if(bands(1) == band_0065)then
      if(phase_liquid)then       
        Ramc = rayleigh_liq(size_re)
      else
        Ramc = rayleigh_ice(size_re)
      endif
            
      CALL calcDistanceAndMinLoc(Ramc,Rbm, &
        temp_refl(1,s_rt-ire_ct_1+1: temp_array_size), &
        temp_refl(2,s_rt-ire_ct_1+1: temp_array_size),lib_dist,loc)
        Ram_corr = Ramc            
    endif
      
    if ( loc == 1 ) then 
      tau = 0.01
    else
      tau = library_taus(loc-1 )
    endif
    
    if(Ramc > maxval(reflibA(1:size_tau,size_re)))then
      tau =  fillvalue_real
      re = fillvalue_real
      lib_dist = MAX_COST_FUNCTION

    endif
              
  endif     
    
  if(tau > library_taus(size_tau-1)-1.0 )then
    tau =  fillvalue_real
    re = fillvalue_real
    lib_dist = MAX_COST_FUNCTION
  endif
    
  ! don't want any 2's and such 
  ! if (re < 4. ) re = 4.
  ! the library taus don't have 0., so we need to explicitly count 
  !  for it. But we can't put in a 0. straight because
  ! L3 code can't handle an explicit 0. value. 

  deallocate (rel_residual,temp_refl,reflibBL)

end subroutine ASL_boundary

subroutine calcDistanceAndMinLoc(R1,R2,R1vec, R2vec,mindist,loc_index)
  !
  !  R1
  !  R2
  !  R1vec
  !  R2vec
  !  mindist
  !  loc_index

  REAL,INTENT(IN)::R1,R2
  REAL,DIMENSION(:),INTENT(IN)::R1vec,R2vec
  REAL,INTENT(OUT)::mindist
  INTEGER,INTENT(OUT)::loc_index

  real,dimension(:),allocatable::rel_residual
  INTEGER::rel_loc(1)

  allocate(rel_residual(size(R1vec)))

  rel_residual = sqrt ( (R1vec - R1)**2 + (R2vec -R2)**2 ) / &
    sqrt (R1**2 + R2**2)
  rel_loc = minloc(rel_residual)  
  mindist = rel_residual(rel_loc(1))
  loc_index = rel_loc(1)
    
  deallocate(rel_residual)
    
  return
end subroutine calcDistanceAndMinLoc

subroutine solve_retrieval_noSWIR(optical_thickness_vector, &
    library_radii, effective_radius, optical_thickness) 
  !
  !  optical_thickness_vector
  !  library_radii
  !  effective_radius
  !  optical_thickness
    
  use nonscience_parameters
  use science_parameters
  use libraryarrays, only:library_taus                 

  real, intent(in)   :: optical_thickness_vector(:)
  real, intent(in) :: library_radii(:)
  real, intent(in) :: effective_radius
  real, intent(out)  :: optical_thickness


  integer :: ilo, ihi

  ! we have no SWIR (or knocked it out on purpose) 
  ! we are fixing the value of effective radius on some value that 
  !  is generally agreed on 
  ! as being reasonable. We will return an optical thickness retrieval 
  !  that's all nice
  ! and interpolated, but the re is going to be a constant
  ! There will be no ASL retrieval because we are just staying on the 
  !  re=const curve in the library

  call findinterpolationrange(3, effective_radius, size(library_radii), &
    real(library_radii), ilo,  ihi)

  optical_thickness = lagrangeinterp(effective_radius, &
    real(library_radii(ilo:ilo+2)), optical_thickness_vector(ilo:ilo+2) )
                                      
  if (optical_thickness > 150.) optical_thickness = 150.
  if (optical_thickness < 0.) then 
    optical_thickness = fillvalue_real
  endif

end subroutine solve_retrieval_noSWIR

subroutine solveretrieval( residual, optical_thickness_vector, &
  library_radii, effective_radius, optical_thickness, debug, &
  use_nearest, quality_out)
  !
  !  residual IN  vector of fractional diff of rho(abs) measured - 
  !      rho(abs) of the rho library for the set of solved tau (by the 
  !      non-abs algorithm) at each re (make sense?)
  !  optical_thickness_vector IN  tau solved for at each re
  !  library_radii  IN  the re axis values
  !  effective_radius  OUT  solved re
  !  optical_thickness  OUT  solved tau
  !  debug  IN  switch for debugging
  !  use_nearest  OUT  apparently, this failed and solveretrieval_nearest needs
  !       to be called
  !  quality_out  OUT  many states of the re retrieval

  use GeneralAuxType
  use nonscience_parameters
  use science_parameters
  use libraryarrays, only:library_taus                    

  real, intent(in)   :: residual(:), optical_thickness_vector(:)
  real(single), intent(in) :: library_radii(:)
  logical, intent(in)  :: debug
  logical, intent(inout) :: use_nearest
  real, intent(out)  :: effective_radius, optical_thickness
  integer,intent(out)::quality_out

  real :: max_allowable_tau,local_max_allowable           !!added
  integer:: total_taus                                    !!added

  integer              :: number_local_minima, ilo, ihi, zero_count
  integer  :: local_minima(8)

  INTEGER, parameter   :: nre=18
  INTEGER :: crossing_num, crossing_vector(nre), local_max_num, &
    local_max_vector(nre), local_min_num, local_min_vector(nre), k
  integer(integer_onebyte) :: quality

  integer:: res_size

  real  re_water_outofbounds_high, re_water_outofbounds_low
  logical :: re_altered, correction_made

  local_minima(:) = 0
  re_water_outofbounds_high = 30.; re_water_outofbounds_low = 4.
  
  correction_made = .FALSE.
  
  ! get effective radius solution                        !added/changed
  ! changed for water residual going in from re = 4
  !  This decides if we are doing ice or water retr by checking the 
  !  library array size  The find_zero_crossings, find_local_minima,
  !  find_local_maxima gets the indicies of the tops, bottoms and 
  !  (value before) the zero crossing of the residual and the # of these
  !  Then, solution_re
  if(is_water_phase(library_radii))then
    call find_zero_crossings (residual(TOP_NK_ASL_CONTOUR_WATER:), &
      crossing_vector,  crossing_num)
    call find_local_minima   (residual(TOP_NK_ASL_CONTOUR_WATER:), &
      local_min_vector, local_min_num)
    call find_local_maxima   (residual(TOP_NK_ASL_CONTOUR_WATER:), &
      local_max_vector, local_max_num)

    call solution_re (library_radii(TOP_NK_ASL_CONTOUR_WATER:), &
      residual(TOP_NK_ASL_CONTOUR_WATER:) , crossing_vector, crossing_num, &
      local_min_vector, local_min_num, local_max_vector, local_max_num, &
      effective_radius, quality )
  else
    call find_zero_crossings (residual, crossing_vector,  crossing_num)
    call find_local_minima   (residual, local_min_vector, local_min_num)
    call find_local_maxima   (residual, local_max_vector, local_max_num)

    call solution_re (library_radii, residual, crossing_vector, crossing_num, &
      local_min_vector, local_min_num, local_max_vector, local_max_num, &
                        effective_radius, quality )
  endif

  quality_out = quality   !added/changed
    
  use_nearest = .false.
  if ( quality == -2 .or. quality == -3 .or. quality == -4 &
    .OR. quality == -6) then    !added/changed

    use_nearest = .true.
    RETURN

  endif
  ! get optical thickness solution

  ! PARTIAL RETRIEVAL LOGIC ADDENDUM FOR RETENTION OF TAU SOLUTION 
  !  FOR FAILED RE
  !  order is important here--the following block needs to be done 
  !  here before the call to findinterpolationrange
  ! for performing a partial tau retrieval after having retrieved the 
  !  maximum or minimum desireable, or failed to retrieve liquid library re

  ! sep: 3=8-06 
  !  For water clouds there are 4 possibilities at this point: re>=30, 
  !  re<=4, re=INVALID, or re valid
  !  For ice clouds there are 2 possibilities at this point: re=INVALID 
  !  or re valid
  !  It was decided to perform "partial" optical thickness retrievals 
  !  for any non-valid case by temporarily assigning re a 
  !  value of 10 or 30 microns for ice and water, respectively. 
  !  The logic to account for this is simpler than the previous version.

  re_altered = .false.
  if (is_water_phase(library_radii) .and. &  ! changed/added <= to <
    ( (effective_radius < re_water_outofbounds_low) .or. & 
    (effective_radius == INVALID_ATTEMPTED_BUT_FAILED_RE) ) ) then!{
    effective_radius = 10.
    re_altered = .true.
  else!}{
  ! for performing a partial tau retrieval after having failed to 
  !  retrieve ice library re
    if (effective_radius == INVALID_ATTEMPTED_BUT_FAILED_RE) then!{
      effective_radius = 30.
      re_altered = .true.
    endif!}
  endif!}
  ! END PARTIAL RETRIEVAL LOGIC ADDENDUM FOR RETENTION OF TAU 
  !   SOLUTION FOR FAILED RE 

  call findinterpolationrange(3, effective_radius, size(library_radii), &
    real(library_radii), ilo, ihi)
                          
  !!! added !!!  Cloud folks emphasis
  !!! correction starts here !!!!!!!!
  total_taus = size(library_taus)
  max_allowable_tau = library_taus(total_taus)
  local_max_allowable = max_allowable_tau - 1.0

  !  we go in if tau_vector(ilo:ilo+2) has a 158.78 
  !  or -1 (i.e, lagrange interpolation points)
  if(.not.(use_nearest) .and. &
    (maxval(optical_thickness_vector(ilo:ilo+2)) > local_max_allowable  &
    .OR. minval(optical_thickness_vector(ilo:ilo+2)) < 0))then

    correction_made = .TRUE.                                    
    if(ilo <= 2)then
      if(maxval(optical_thickness_vector(1:2)) > local_max_allowable  & 
        .OR. minval(optical_thickness_vector(1:2)) < 0)then 
        optical_thickness = fillvalue_real
        use_nearest = .true.
        quality_out = -4
      else
        ilo =1
        optical_thickness = optical_thickness_vector(ilo) + &
          ( ( effective_radius - real(library_radii(ilo))) /  &
          (real(library_radii(ilo+1)-library_radii(ilo))) ) * &
          (optical_thickness_vector(ilo+1) - optical_thickness_vector(ilo) )  
      endif

    elseif(maxval(optical_thickness_vector(ilo-1:ilo)) > local_max_allowable  &
      .OR. minval(optical_thickness_vector(ilo-1:ilo)) < 0)then
      optical_thickness = fillvalue_real
      use_nearest = .true.
      quality_out = -4
            
    elseif(maxval(optical_thickness_vector(ilo-1:ilo)) < local_max_allowable & 
      .AND. minval(optical_thickness_vector(ilo-1:ilo)) > 0)then
      ilo = ilo-2
      if(optical_thickness_vector(ilo+3) < local_max_allowable & 
        .AND. optical_thickness_vector(ilo+3) > 0)ilo = ilo+1

      optical_thickness = optical_thickness_vector(ilo+1) + &
        ( ( effective_radius - real(library_radii(ilo+1))) /  &
        (real(library_radii(ilo+2)-library_radii(ilo+1))) ) * &
        (optical_thickness_vector(ilo+2) - optical_thickness_vector(ilo+1) )   
    else
      !  this else can be removed. If all above fails it should come here. 
      !  Tested 5 granuels and didn't come here
      !  but to be safe make it filled
      !  optical_thickness = fillvalue_real
      use_nearest = .true.
      quality_out = -4
    endif    
  endif

  !!! correction done !!!! The end of their correction

  if(.NOT.(correction_made))then 
    optical_thickness = lagrangeinterp(effective_radius, &
      real(library_radii(ilo:ilo+2)), &
      optical_thickness_vector(ilo:ilo+2) )
  endif

  if( optical_thickness < optical_thickness_vector(ilo) .AND. &
    (.NOT.(correction_made))) then!{    !!!!!!added /changed

    call findinterpolationrange(2, effective_radius, size(library_radii), &
      real(library_radii), ilo,  ihi)
      optical_thickness = optical_thickness_vector(ilo) + &
      ( ( effective_radius - real(library_radii(ilo))) / &
      (real(library_radii(ilo+1)-library_radii(ilo))) ) * &
      (optical_thickness_vector(ilo+1) - optical_thickness_vector(ilo) ) 
  endif!}
  
  ! KILL RE FOR PARTIAL RETRIEVAL
  ! Now that the optical thickness has been determined, set the 
  !  effective radius to fill if re_altered = .TRUE.
  if (re_altered) effective_radius = fillvalue_real

  ! RANGE CHECKING of optical thickness and effective radius

  ! OPTICAL THICKNESS 

  ! tau low boundary of range
  if (optical_thickness /= fillvalue_real .and. optical_thickness < 0.) then 
    optical_thickness = fillvalue_real
    effective_radius = fillvalue_real
    use_nearest = .true.
    quality_out = -4
  endif 
  
  if (optical_thickness /= fillvalue_real .and. &
    optical_thickness < 0.01 .and. optical_thickness > 0.) then!{ 
    optical_thickness = 0.01
  endif!}
  ! end tau low boundary of range

  ! tau high boundary of range
  !  QA for large optical thickness
  !   if (optical_thickness > 150.) then!{
  !     optical_thickness_outofbounds = 2
  !   elseif (optical_thickness > 100. .and. optical_thickness <= 150.) then!}{
  !     optical_thickness_outofbounds = 1
  !   else!}{
  ! optical_thickness_outofbounds = 0
  !   endif!}

  !   if ( optical_thickness > 100.) then!{
  !     optical_thickness = 100.
  !   endif!}
  ! end tau high boundary of range

  ! EFFECTIVE RADIUS 

  ! NOTE: Ideally, it would be more satisfying to assign re retrievals 
  !  outside of the library space to fill value,
  !  even here after other checks have done the same however, retrievals 
  !  can be within epsilon of boundary 
  !  (e.g., Liam's pixel with ice re = 90.00000763).  That is the intent 
  !  of this strict range enforcement here.

  ! re low boundary of range
  if (effective_radius /= fillvalue_real .and. &
    effective_radius < library_radii(1)) then!{ 
    effective_radius = library_radii(1)
  endif!}
  ! end re low boundary of range

  ! re high boundary of range
  if (effective_radius > library_radii(size(library_radii))) then!{ 
     effective_radius = library_radii(size(library_radii))
  endif!}
  ! end re high boundary of range

end subroutine solveretrieval

subroutine linear_interpolate_for_root(residual, rad, re)
  !
  !  residual
  !  rad
  !  re
  !
  implicit none

  real, intent(in)     :: residual(2), rad(2)
  real, intent(out)    :: re

  re = -1.0*( rad(1)*residual(2)/(residual(1)-residual(2)) + &
    rad(2)*residual(1)/(residual(2)-residual(1)) )

end subroutine linear_interpolate_for_root

subroutine findinterpolationrange( n1, xx, n,x, k1, k2)
  !
  !  n1
  !  xx
  !  n
  !  x
  !  k1
  !  k2
  !
  implicit none

  integer, intent(in)   :: n1, n
  real, intent(in)      :: xx, x(n)

  integer,intent(out)   :: k1, k2
  
  integer               :: i, n2

  if (n > n1) then!{
    n2 = 0.5*n1
    i = 1
    do while (xx > x(i) .and. i < n)
      i = i + 1
    enddo
    k1 = max(1,i-n2)
    k2 = k1 + n1 - 1
    if (k2 > n) then!{
      k2 = n
      k1 = k2 - n1 + 1
    end if
  else!}{
  !    if # of exist data set n <= n1, i
  !    use n1 = n for interp_mas_geoolation.
    k1 = 1
    k2 = n
  end if

end subroutine findinterpolationrange

real function lagrangeinterp( xx, x, y)
  ! in the previous interation of the COT retrieval code this was 
  !  known as 'BINTPB'
  !
  !  xx
  !  x
  !  y
  !
  implicit none
  real, intent(in)  :: x(3), y(3), xx 
  lagrangeinterp = (xx-x(2))* (xx-x(3))/ (x(1)-x(2))/ (x(1)-x(3))*y(1) + &
                   (xx-x(3))* (xx-x(1))/ (x(2)-x(3))/ (x(2)-x(1))*y(2) + &
                   (xx-x(1))* (xx-x(2))/ (x(3)-x(1))/ (x(3)-x(2))*y(3)
end function lagrangeinterp

  ! sssssssssssssssssssssssssssssssssssssssssssssss
subroutine check_for_signchange (x, signchange)
  ! sssssssssssssssssssssssssssssssssssssssssssssss
  !
  !  x
  !  signchange
  !
  !  *** Copied from collect 5 core_science code

  ! this routine reports on sign change for a set of 3 values
  ! signchnage == 0 -> no sign change
  ! signchange == 1 -> sign changes between x1 and x2
  ! signchange == 2 -> sign changes bweteen x2 and x3

  implicit none

  real, intent(in)      :: x(3)
  integer, intent(out)  :: signchange

  if ( (x(1)==0 .and. x(2)==0) .or. (x(1)==0 .and. x(3)==0) .or. &
    (x(2)==0 .and. x(3)==0) ) Then
    signchange = 0
  else!}{
    if ( sign(1.,x(1))*sign(1.,x(2))*sign(1.,x(3)) /= -1.) Then
      signchange =1
      if ( x(1) > 0 .AND. x(2) > 0 .AND. x(3) > 0) then!{
        signchange =0
      end if

      if ( x(1) < 0 .AND. x(2) < 0 .AND. x(3) < 0) then!{
        signchange =0
      end if

    else!}{
      signchange =1
    endif!}
  endif!}

  if (signchange == 1) then!{
    if ( sign(1.,x(1))*sign(1.,x(2)) == -1. ) signchange = 1
    if ( sign(1.,x(2))*sign(1.,x(3)) == -1. ) signchange = 2
  endif!}

end subroutine check_for_signchange

subroutine quad_interpolate_for_root ( radii, residual, radius_solution, &
    status)
  !
  !  radii
  !   residual
  !  radius_solution
  !  status
  !
  use GeneralAuxType
  implicit none

  real(single), intent(in)  :: radii(:), residual(:)

  real(single), intent(out) :: radius_solution
  integer,      intent(out) :: status

  real(double) :: d1, d2, d3, a0, a1, a2, root1, root2, z

  d1 = residual(1) / ((radii(1)-radii(2))*(radii(1)-radii(3)))
  d2 = residual(2) / ((radii(2)-radii(1))*(radii(2)-radii(3)))
  d3 = residual(3) / ((radii(3)-radii(1))*(radii(3)-radii(2)))

  a0 = d1*radii(2)*radii(3) + d2*radii(1)*radii(3) + d3*radii(1)*radii(2)
  a1 = -d1*(radii(2)+radii(3)) - d2*(radii(1)+radii(3)) - d3*(radii(1)+radii(2))
  a2 = d1 + d2 + d3


  z = a1**2 - 4*a0*a2

  if(z >= 0. .and. a2 /= 0.) then!{
    root1 = (-a1 + sqrt(z))/(2*a2)
    root2 = (-a1 - sqrt(z))/(2*a2)
  else!}{
    root1 = -999; root2 = -999
    status = 1
    return
  end if
  if (root1 < radii(1) .and. root2 < radii(1)) then!{
    radius_solution =-999.
    status = 2
  elseif (root1 > radii(3) .and. root2 > radii(3)) then!}{
    radius_solution =-999.
    status = 3
  elseif (root1 >= radii(1) .and. root1 <= radii(3)) then!}{
    radius_solution = root1
    status = 4
  elseif(root2 >= radii(1) .and. root2 <= radii(3)) then!}{
    radius_solution = root2
    status = 5
  else!}{
    radius_solution = -999.
    status = 1
  endif!}

end subroutine quad_interpolate_for_root

end module retrieval_solution_logic
