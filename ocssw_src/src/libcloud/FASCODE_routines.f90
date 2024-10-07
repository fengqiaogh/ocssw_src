module FASCODE_routines

	implicit none

! This is a rewrite / translation to Fortran 90 of a set of subroutines that do the forward
! model calculations for CT and OD codes using UW-Madison coefficient files
!		-- Gala Wind 11.30.2011

!c * PLOD/PFAAST subprograms for infrared transmittance
!C 	at 101-level AIRS/MODIS pressure coordinate
!c .... version of 18.03.03

!contents:

!c     block data reference_atmosphere
!c     subroutine calpir
!c     subroutine conpir
!c     subroutine gphite
!c     subroutine taudoc
!c     subroutine tauwtr

!c *********************
#ifndef AIRBORNE 

	integer, parameter :: nl = 101

	real, parameter :: pstd(nl) = (/0.0050,    0.0161,    0.0384,    0.0769,    0.1370, &
         0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,&
         1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,&
         5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,&
        14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,&
        29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,&
        51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,&
        83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,&
       125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,&
       180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,&
       247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,&
       328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,&
       424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,&
       535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,&
       661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,&
       802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,&
       958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/)
	
	
	real, parameter :: tstd(nl) = (/190.19, 203.65, 215.30, 226.87, 237.83, &
       247.50, 256.03, 263.48, 267.09, 270.37, 266.42, 261.56, 256.40, &
       251.69, 247.32, 243.27, 239.56, 236.07, 232.76, 230.67, 228.71, &
       227.35, 226.29, 225.28, 224.41, 223.61, 222.85, 222.12, 221.42, &
       220.73, 220.07, 219.44, 218.82, 218.23, 217.65, 217.18, 216.91, &
       216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, &
       216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.71, &
       216.71, 216.72, 216.81, 217.80, 218.77, 219.72, 220.66, 222.51, &
       224.57, 226.59, 228.58, 230.61, 232.61, 234.57, 236.53, 238.48, &
       240.40, 242.31, 244.21, 246.09, 247.94, 249.78, 251.62, 253.45, &
       255.26, 257.04, 258.80, 260.55, 262.28, 264.02, 265.73, 267.42, &
       269.09, 270.77, 272.43, 274.06, 275.70, 277.32, 278.92, 280.51, &
       282.08, 283.64, 285.20, 286.74, 288.25, 289.75, 291.22, 292.68/)
	
	real, parameter :: wstd(nl) = (/  0.001,  0.001,  0.002,  0.003,  0.003, &
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
        0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, &
        0.003,  0.003,  0.004,  0.004,  0.005,  0.005,  0.007,  0.009, &
        0.011,  0.012,  0.014,  0.020,  0.025,  0.030,  0.035,  0.047, &
        0.061,  0.075,  0.089,  0.126,  0.162,  0.197,  0.235,  0.273, &
        0.310,  0.356,  0.410,  0.471,  0.535,  0.601,  0.684,  0.784, &
        0.886,  0.987,  1.094,  1.225,  1.353,  1.519,  1.686,  1.852, &
        2.036,  2.267,  2.496,  2.721,  2.947,  3.170,  3.391,  3.621, &
        3.848,  4.084,  4.333,  4.579,  4.822,  5.061,  5.296,  5.528/)
	
	real, parameter :: ostd(nl) = (/0.47330,0.27695,0.28678,0.51816,0.83229, &
      1.18466,1.69647,2.16633,3.00338,3.76287,4.75054,5.61330,6.33914, &
      7.03675,7.50525,7.75612,7.81607,7.69626,7.56605,7.28440,7.01002, &
      6.72722,6.44629,6.17714,5.92914,5.69481,5.47387,5.26813,5.01252, &
      4.68941,4.35141,4.01425,3.68771,3.37116,3.06407,2.77294,2.50321, &
      2.24098,1.98592,1.74840,1.54451,1.34582,1.17824,1.02513,0.89358, &
      0.78844,0.69683,0.62654,0.55781,0.50380,0.45515,0.42037,0.38632, &
      0.35297,0.32029,0.28832,0.25756,0.22739,0.19780,0.16877,0.14901, &
      0.13190,0.11511,0.09861,0.08818,0.07793,0.06786,0.06146,0.05768, &
      0.05396,0.05071,0.04803,0.04548,0.04301,0.04081,0.03983,0.03883, &
      0.03783,0.03685,0.03588,0.03491,0.03395,0.03368,0.03349,0.03331, &
      0.03313,0.03292,0.03271,0.03251,0.03190,0.03126,0.03062,0.02990, &
      0.02918,0.02850,0.02785,0.02721,0.02658,0.02596,0.02579,0.02579/)


	  integer, parameter :: Fnm = 101-1
	  
#else
	  integer, parameter :: Fnm = 66-1

#endif

	  integer, parameter :: Fncd = 8+1, Fnco = 9+1, Fncl = 2+1, Fncs=11+1, Fncc=4+1 

#ifdef AHI_INST
	integer, parameter :: Fnr = 10
      real ::	coefd(Fncd,Fnm,Fnr),&
				coefo(Fnco,Fnm,Fnr), &
                coefl(Fncl,Fnm,Fnr),&
				coefs(Fncs,Fnm,Fnr), &
                coefc(Fncc,Fnm,Fnr)
	
#endif

#ifdef SEVIRI_INST
	integer, parameter :: Fnr = 8
      real ::	coefd(Fncd,Fnm,Fnr),&
				coefo(Fnco,Fnm,Fnr), &
                coefl(Fncl,Fnm,Fnr),&
				coefs(Fncs,Fnm,Fnr), &
                coefc(Fncc,Fnm,Fnr)
	
#endif

#ifdef MAS_INST
	integer, parameter :: Fnr = 22
      real ::	coefd(Fncd,Fnm,Fnr),&
				coefo(Fnco,Fnm,Fnr), &
                coefl(Fncl,Fnm,Fnr),&
				coefs(Fncs,Fnm,Fnr), &
                coefc(Fncc,Fnm,Fnr)
#endif

#ifdef EMAS_INST
	integer, parameter ::  Fnr = 13
      real ::	coefd(Fncd,Fnm,Fnr),&
				coefo(Fnco,Fnm,Fnr), &
                coefl(Fncl,Fnm,Fnr),&
				coefs(Fncs,Fnm,Fnr), &
                coefc(Fncc,Fnm,Fnr)
#endif

	real :: pavg(Fnm),tref(Fnm),wref(Fnm),oref(Fnm)
	real ::  tavg(Fnm),wamt(Fnm),oamt(Fnm),secz(Fnm), secz_2way(Fnm)


	integer, parameter :: Fnxd = 8, Fnxo=9, Fnxw = 2 + 11, Fnxc =4

	real ::  xdry(Fnxd,Fnm),xozo(Fnxo,Fnm),xwet(Fnxw,Fnm),xcon(Fnxc,Fnm)
	real ::  xdry_2way(Fnxd,Fnm),xozo_2way(Fnxo,Fnm),xwet_2way(Fnxw,Fnm),xcon_2way(Fnxc,Fnm)
	
	character*6 cinit
	integer :: init 

contains

	subroutine init_FASCODE
	
#ifdef AIRBORNE
		cinit = 'zz'
#else
		cinit = 'zzzzzz'
#endif	
	
		init = 0
		
		xdry = 0.
		xozo = 0.
		xwet = 0.
		xcon = 0.
		
		xdry_2way = 0.
		xozo_2way = 0.
		xwet_2way = 0.
		xcon_2way = 0.
		
	
	end subroutine init_FASCODE

	real function secant (z)
		
		real, intent(in) :: z
			
		secant = 1./cos(0.01745329*z)
		
	end function secant


      subroutine calpir(t_avg_ref,  amt_wet_ref, amt_ozo_ref, &
                             t_avg,      amt_wet,     amt_ozo, &
                             p_avg,      sec_theta,   n_layers, &
                             n_dry_pred, n_wet_pred,  n_ozo_pred, &
                             n_con_pred, &
                             pred_dry,   pred_wet,    pred_ozo, &
                             pred_con, do_init) 
! ... version of 18.03.03

!  PURPOSE:

!    Routine to calculate the predictors for the dry (temperature), 
!      wet and ozone components of a fast transmittance model for a
!      scanning satellite based instrument.

!  REFERENCES:

!    AIRS FTC package science notes and software, S. Hannon and L. Strow,
!      Uni. of Maryland, Baltimore County (UMBC)

!  CREATED:

!    19-Sep-1996 HMW

!  ARGUMENTS:

!      Input
!    -----------
!     t_avg_ref  - REAL*4 reference layer average temperature array (K)

!    amt_wet_ref - REAL*4 reference water vapour amount array (k.mol)/cm^2

!    amt_ozo_ref - REAL*4 reference ozone amount array (k.mol)/cm^2

!      t_avg     - REAL*4 layer average temperature array (K)

!     amt_wet    - REAL*4 water vapour amount array (k.mol)/cm^2

!     amt_ozo    - REAL*4 ozone amount array (k.mol)/cm^2

!      p_avg     - REAL*4 layer average pressure array (mb)

!    sec_theta   - REAL*4 secant of the zenith angle array

!     n_layers   - INT*4 Number of atmospheric layers

!    n_dry_pred  - INT*4 number of dry (temperature) predictors

!    n_wet_pred  - INT*4 number of water vapour predictors

!    n_ozo_pred  - INT*4 number of ozone predictors

!    n_con_pred  - INT*4 number of water vapour continuum predictors

!      Output
!    -----------
!     pred_dry   - REAL*4 dry gas (temperature) predictor matrix

!     pred_wet   - REAL*4 water vapour predictor matrix

!     pred_ozo   - REAL*4 ozone predictor matrix

!     pred_con   - REAL*4 water vapour continuum predictor matrix

!  COMMENTS:

!    Levels or Layers?
!    -----------------
!      Profile data is input at a number of *LAYERS*.

!    Layer Numbering pt. A
!    ---------------------
!      Layer 1   => Atmosphere between LEVELs 1 & 2
!      Layer 2   => Atmosphere between LEVELs 2 & 3
!                        .
!                        .
!                        .
!      Layer L-1 => Atmosphere between LEVELs L-1 & L

!    Layer Numbering pt. B
!    ---------------------
!      For the HIS instrument, Layer 1 is at the top of the atmosphere
!        and Layer L-1 is at the surface.    

!    Layer Numbering pt. C
!    ---------------------
!      In this routine the number of *LAYERS* is passed in the argument
!        list, _not_ the number of LEVELS.  This was done to improve
!        the readability of this code, i.e. loop from 1->L(ayers) 
!        rather than from 1->L(evels)-1.

!=======================================================================

!-----------------------------------------------------------------------
!                 Turn off implicit type declaration
!-----------------------------------------------------------------------

       implicit none

!------------------------------------------------------------------------
!                             Arguments
!------------------------------------------------------------------------

! -- Input

      integer*4 :: n_layers, &
               n_dry_pred, n_wet_pred, n_ozo_pred, n_con_pred

      real*4  ::  t_avg_ref(*), amt_wet_ref(*), amt_ozo_ref(*),&
               t_avg(*),     amt_wet(*),     amt_ozo(*), &
               p_avg(*),     sec_theta(*)
! WDR in/output
      integer ::  do_init

! -- Output

      real*4 ::   pred_dry(n_dry_pred, *), &
               pred_wet(n_wet_pred, *),&
               pred_ozo(n_ozo_pred, *),&
               pred_con(n_con_pred, *)

!------------------------------------------------------------------------
!                           Local variables
!------------------------------------------------------------------------

! -- Parameters
#ifdef AIRBORNE
     integer, parameter ::  max_layers = 65
#else
      integer, parameter ::  max_layers = 100
#endif
      integer, parameter ::  max_dry_pred = 8, &
							max_wet_pred = 13,&
							 max_ozo_pred = 9,&
							  max_con_pred = 4

! -- Scalars

      integer*4 :: l

! -- Arrays

!     ....Pressure
      real*4  ::  p_dp(max_layers),&
               p_norm(max_layers)

!     ....Temperature
      real*4   :: delta_t(max_layers), &
               t_ratio(max_layers), &
               pw_t_ratio(max_layers)      ! Pressure weighted

!     ....Water vapour
      real*4  ::  wet_ratio(max_layers),&
               pw_wet(max_layers),   &      ! Pressure weighted
               pw_wet_ref(max_layers),  &    ! Pressure weighted
               pw_wet_ratio(max_layers)    ! Pressure weighted

!     ....Ozone
      real*4   :: ozo_ratio(max_layers), &
               pw_ozo_ratio(max_layers), &  ! Pressure weighted
               pow_t_ratio(max_layers)     ! Pressure/ozone weighted

!************************************************************************
!                         ** Executable code **
!************************************************************************

!------------------------------------------------------------------------
!                   -- Check that n_layers is o.k. --
!------------------------------------------------------------------------

! WDR to avoid init parts that are the same
    if( do_init == 1 ) then
      if( n_layers > max_layers )then
        write(*,'(/10x,''*** calpir : n_layers > max_layers'')')
        stop
      end if 

!------------------------------------------------------------------------
!         -- Check that numbers of predictors is consistent --
!------------------------------------------------------------------------

!     ---------------------------------
!     # of dry (temperature) predictors
!     ---------------------------------

      if( n_dry_pred /= max_dry_pred )then
        write(*,'(/10x,''*** calpir : invalid n_dry_pred'')')
        stop
      end if 

!     ----------------------------
!     # of water vapour predictors
!     ----------------------------

      if( n_wet_pred /= max_wet_pred )then
        write(*,'(/10x,''*** calpir : invalid n_wet_pred'')')
        stop
      end if 

!     ---------------------
!     # of ozone predictors
!     ---------------------

      if( n_ozo_pred /= max_ozo_pred )then
        write(*,'(/10x,''*** calpir : invalid n_ozo_pred'')')
        stop
      end if 

!     --------------------------------------
!     # of water vapour continuum predictors
!     --------------------------------------

      if( n_con_pred /= max_con_pred )then
        write(*,'(/10x,''*** calpir : invalid n_con_pred'')')
        stop
      end if 

!------------------------------------------------------------------------
!         -- Calculate ratios, offsets, etc, for top layer --
!------------------------------------------------------------------------

!     ------------------
!     Pressure variables
!     ------------------

      p_dp(1)   = p_avg(1) * ( p_avg(2) - p_avg(1) )
      p_norm(1) = 0.0

!     ---------------------
!     Temperature variables
!     ---------------------

      delta_t(1)    = t_avg(1) - t_avg_ref(1)
      t_ratio(1)    = t_avg(1) / t_avg_ref(1)
      pw_t_ratio(1) = 0.0

!     ----------------
!     Amount variables
!     ----------------

!     ....Water vapour
 
      wet_ratio(1)    = amt_wet(1) / amt_wet_ref(1)
      pw_wet(1)       = p_dp(1) * amt_wet(1)
      pw_wet_ref(1)   = p_dp(1) * amt_wet_ref(1)
      pw_wet_ratio(1) = wet_ratio(1)

!     ....Ozone

      ozo_ratio(1)    = amt_ozo(1) / amt_ozo_ref(1)
      pw_ozo_ratio(1) = 0.0
      pow_t_ratio(1)  = 0.0

!------------------------------------------------------------------------
!         -- Calculate ratios, offsets, etc, for all layers --
!------------------------------------------------------------------------

      do l = 2, n_layers

!       ------------------
!       Pressure variables
!       ------------------

        p_dp(l) = p_avg(l) * ( p_avg(l) - p_avg(l-1) )
        p_norm(l) = p_norm(l-1) + p_dp(l)

!       ---------------------
!       Temperature variables
!       ---------------------

        delta_t(l)    = t_avg(l) - t_avg_ref(l)
        t_ratio(l)    = t_avg(l) / t_avg_ref(l)
        pw_t_ratio(l) = pw_t_ratio(l-1) + ( p_dp(l) * t_ratio(l-1) )

!       ----------------
!       Amount variables
!       ----------------

!       ..Water vapour

        wet_ratio(l)  = amt_wet(l) / amt_wet_ref(l)
        pw_wet(l)     = pw_wet(l-1) + ( p_dp(l) * amt_wet(l) )
        pw_wet_ref(l) = pw_wet_ref(l-1) + ( p_dp(l) * amt_wet_ref(l) )
        
!       ..Ozone

        ozo_ratio(l)    = amt_ozo(l) / amt_ozo_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l-1) + &
                           ( p_dp(l) * ozo_ratio(l-1) )
        pow_t_ratio(l)  = pow_t_ratio(l-1) + &
                           ( p_dp(l) * ozo_ratio(l-1) * delta_t(l-1) )

      end do

!------------------------------------------------------------------------
!              -- Scale the pressure dependent variables --
!------------------------------------------------------------------------

      do l = 2, n_layers

        pw_t_ratio(l)   = pw_t_ratio(l) / p_norm(l)
        pw_wet_ratio(l) = pw_wet(l) / pw_wet_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l) / p_norm(l)
        pow_t_ratio(l)  = pow_t_ratio(l) / p_norm(l)
 
      end do

   ! WDR end init portion
     do_init = 0
   end if
!------------------------------------------------------------------------
!                     -- Load up predictor arrays --
!------------------------------------------------------------------------

      do l = 1, n_layers

!       ----------------------
!       Temperature predictors
!       ----------------------

        pred_dry(1,l) = sec_theta(l)
        pred_dry(2,l) = sec_theta(l) * sec_theta(l)
        pred_dry(3,l) = sec_theta(l) * t_ratio(l)
        pred_dry(4,l) = pred_dry(3,l) * t_ratio(l)
        pred_dry(5,l) = t_ratio(l)
        pred_dry(6,l) = t_ratio(l) * t_ratio(l)
        pred_dry(7,l) = sec_theta(l) * pw_t_ratio(l)
        pred_dry(8,l) = pred_dry(7,l) / t_ratio(l) 

!       -----------------------
!       Water vapour predictors
!       -----------------------

        pred_wet(1,l)  = sec_theta(l) * wet_ratio(l)
        pred_wet(2,l)  = sqrt( pred_wet(1,l) )
        pred_wet(3,l)  = pred_wet(1,l) * delta_t(l)
        pred_wet(4,l)  = pred_wet(1,l) * pred_wet(1,l)
        pred_wet(5,l)  = abs( delta_t(l) ) * delta_t(l) * pred_wet(1,l)
        pred_wet(6,l)  = pred_wet(1,l) * pred_wet(4,l)
        pred_wet(7,l)  = sec_theta(l) * pw_wet_ratio(l)
        pred_wet(8,l)  = pred_wet(2,l) * delta_t(l)
        pred_wet(9,l)  = sqrt( pred_wet(2,l) )
        pred_wet(10,l) = pred_wet(7,l) * pred_wet(7,l)
        pred_wet(11,l) = sqrt( pred_wet(7,l) )
        pred_wet(12,l) = pred_wet(1,l)
! +++ old
!       pred_wet(13,l) = pred_wet(2,l)
! +++ new
        pred_wet(13,l) = pred_wet(1,l) / pred_wet(10,l)

!       ----------------
!       Ozone predictors
!       ----------------

        pred_ozo(1,l) = sec_theta(l) * ozo_ratio(l)
        pred_ozo(2,l) = sqrt( pred_ozo(1,l) )
        pred_ozo(3,l) = pred_ozo(1,l) * delta_t(l)
        pred_ozo(4,l) = pred_ozo(1,l) * pred_ozo(1,l)
        pred_ozo(5,l) = pred_ozo(2,l) * delta_t(l)
        pred_ozo(6,l) = sec_theta(l) * pw_ozo_ratio(l)
        pred_ozo(7,l) = sqrt( pred_ozo(6,l) ) * pred_ozo(1,l)
        pred_ozo(8,l) = pred_ozo(1,l) * pred_wet(1,l)
        pred_ozo(9,l) = sec_theta(l) * pow_t_ratio(l) * pred_ozo(1,l)

!       ---------------------------------
!       Water vapour continuum predictors
!       ---------------------------------

        pred_con(1,l) = sec_theta(l) * wet_ratio(l) / ( t_ratio(l) * t_ratio(l) )  
        pred_con(2,l) = pred_con(1,l) * pred_con(1,l) / sec_theta(l)
        pred_con(3,l) = sec_theta(l) * wet_ratio(l) / t_ratio(l) 
        pred_con(4,l) = pred_con(3,l) * wet_ratio(l)

      end do
         
      end subroutine calpir
	  

      subroutine conpir( p, t, w, o, n_levels, i_dir, &
                              p_avg, t_avg, w_amt, o_amt)
! ... version of 19.09.96

!  PURPOSE:

!    Function to convert atmospheric water vapour (g/kg) and ozone (ppmv)
!      profiles specified at n_levels layer BOUNDARIES to n_levels-1
!      integrated layer amounts of units (k.moles)/cm^2.  The average
!      LAYER pressure and temperature are also returned.

!  REFERENCES:

!    AIRS LAYERS package science notes, S. Hannon and L. Strow, Uni. of
!      Maryland, Baltimore County (UMBC)

!  CREATED:

!    19-Sep-1996 HMW

!  ARGUMENTS:

!     Input
!    --------
!       p     - REAL*4 pressure array (mb)

!       t     - REAL*4 temperature profile array (K)

!       w     - REAL*4 water vapour profile array (g/kg)

!       o     - REAL*4 ozone profile array (ppmv)

!    n_levels - INT*4 number of elements used in passed arrays

!     i_dir   - INT*4 direction of increasing layer number

!                 i_dir = +1, Level(1) == p(top)         } satellite/AC
!                             Level(n_levels) == p(sfc)  }    case

!                 i_dir = -1, Level(1) == p(sfc)         } ground-based
!                             Level(n_levels) == p(top)  }    case

!     Output
!    --------
!     p_avg   - REAL*4 average LAYER pressure array (mb)

!     t_avg   - REAL*4 average LAYER temperature (K)

!     w_amt   - REAL*4 integrated LAYER water vapour amount array (k.moles)/cm^2

!     o_amt   - REAL*4 integrated LAYER ozone amount array (k.moles)/cm^2

!  ROUTINES:

!    Subroutines:
!    ------------
!      gphite      - calculates geopotential height given profile data.

!    Functions:
!    ----------
!      NONE

!  COMMENTS:

!    Levels or Layers?
!    -----------------
!      Profile data is input at a number of *LEVELS*.  Number densitites
!        are calculated for *LAYERS* that are bounded by these levels.
!        So, for L levels there are L-1 layers.

!    Layer Numbering
!    ---------------
!      Layer 1   => Atmosphere between LEVELs 1 & 2
!      Layer 2   => Atmosphere between LEVELs 2 & 3
!                        .
!                        .
!                        .
!      Layer L-1 => Atmosphere between LEVELs L-1 & L

!=======================================================================

!-----------------------------------------------------------------------
!              -- Prevent implicit typing of variables --
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!                           -- Arguments --
!-----------------------------------------------------------------------

! -- Arrays

      real*4  ::  p(*), t(*), w(*), o(*), &
               p_avg(*), t_avg(*), w_amt(*), o_amt(*)

! -- Scalars

      integer*4 :: n_levels, i_dir

!-----------------------------------------------------------------------
!                         -- Local variables --
!-----------------------------------------------------------------------

! -- Parameters

#ifdef AIRBORNE
      integer*4, parameter :: max_levels = 66 ! Maximum number of layers
#else
      integer*4, parameter :: max_levels = 101 ! Maximum number of layers
#endif
      real, parameter :: r_equator = 6.378388e+06,& ! Earth radius at equator
                         r_polar   = 6.356911e+06, & ! Earth radius at pole
                         r_avg     = 0.5*(r_equator+r_polar) 

      real, parameter ::  g_sfc = 9.80665          ! Gravity at surface
! 
      real*4, parameter :: rho_ref = 1.2027e-12     ! Reference air "density"

      real, parameter :: mw_dryair = 28.97, &       ! Molec. wgt. of dry air (g/mol)
                 mw_h2o    = 18.0152,    &  ! Molec. wgt. of water
                 mw_o3     = 47.9982      ! Molec. wgt. of ozone

      real, parameter ::  R_gas = 8.3143,  &         ! Ideal gas constant (J/mole/K)
                 R_air = 0.9975*R_gas     ! Gas constant for air (worst case) 

! -- Scalars

      integer*4 l, l_start, l_end, l_indx

      real*4    rho1, rho2, p1, p2, w1, w2, o1, o2, z1, z2, &
               c_avg, g_avg, z_avg, w_avg, o_avg, &
               dz, dp, r_hgt, A, B

! -- Arrays

      real*4 ::   z(max_levels),  &            ! Pressure heights (m)
               g(max_levels),     &         ! Acc. due to gravity (m/s/s)
               mw_air(max_levels), &        ! Molec. wgt. of air (g/mol)
               rho_air(max_levels), &       ! air mass density (kg.mol)/m^3
               c(max_levels),        &      ! (kg.mol.K)/(N.m)
               w_ppmv(max_levels)          ! h2o LEVEL amount (ppmv)

!***********************************************************************
!                         ** Executable code **
!***********************************************************************

!-----------------------------------------------------------------------
!           -- Calculate initial values of pressure heights --
!-----------------------------------------------------------------------

      call gphite( p, t, w, 0.0, n_levels, i_dir, z)

!-----------------------------------------------------------------------
!      -- Set loop bounds for direction sensitive calculations --
!      -- so loop iterates from surface to the top             --
!-----------------------------------------------------------------------

      if( i_dir > 0 )then

!       --------------------
!       Data stored top down
!       --------------------

        l_start = n_levels
        l_end   = 1

      else

!       ---------------------
!       Data stored bottom up
!       ---------------------

        l_start = 1
        l_end   = n_levels

      end if

!-----------------------------------------------------------------------
!          -- Air molecular mass and density, and gravity --
!          -- as a function of LEVEL                      --
!-----------------------------------------------------------------------

!     -----------------------
!     Loop from bottom to top
!     -----------------------

      do l = l_start, l_end, -1*i_dir

!       ---------------------------------
!       Convert water vapour g/kg -> ppmv
!       ---------------------------------

        w_ppmv(l) = 1.0e+03 * w(l) * mw_dryair / mw_h2o

!       -----------------------------------------
!       Calculate molecular weight of air (g/mol)
!       ----------------------------------------

        mw_air(l) = ( ( 1.0 - (w_ppmv(l)/1.0e+6) ) * mw_dryair ) + &
                   ( ( w_ppmv(l)/1.0e+06 ) * mw_h2o )

!       ----------------
!       Air mass density
!       ----------------

        c(l) = 0.001 * mw_air(l) / R_air    ! 0.001 factor for g -> kg
        rho_air(l) = c(l) * p(l) / t(l)

!       -------
!       Gravity
!       -------

        r_hgt = r_avg + z(l)                !  m
        g(l) = g_sfc -    &                  !  m/s^2
              g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

      end do
 
!-----------------------------------------------------------------------
!                        -- LAYER quantities --
!-----------------------------------------------------------------------

!     -----------------------
!     Loop from bottom to top
!     -----------------------

      do l = l_start, l_end+i_dir, -1*i_dir

!       -------------------------------------------------------
!       Determine output array index.  This is done so that the
!       output data is always ordered from 1 -> L-1 regardless
!       of the orientation of the input data.  This is true by
!       default only for the bottom-up case.  For the top down
!       case no correction would give output layers from 2 -> L
!       -------------------------------------------------------

        if( i_dir > 0 )then

          l_indx = l - 1

        else

          l_indx = l

        end if

!       ---------------------------------------
!       Assign current layer boundary densities
!       ---------------------------------------
 
        rho1 = rho_air(l)
        rho2 = rho_air(l-i_dir)
 
!       ---------
!       Average c
!       ---------

        c_avg = ( (rho1*c(l)) + (rho2*c(l-i_dir)) ) / ( rho1 + rho2 )

!       ---------
!       Average t
!       ---------

        t_avg(l_indx) = ( (rho1*t(l)) + (rho2*t(l-i_dir)) ) / ( rho1 + rho2 )

!       ---------
!       Average p
!       ---------

        p1 = p(l)
        p2 = p(l-i_dir)

        z1 = z(l)
        z2 = z(l-i_dir)

        dp = p2 - p1

        A = log(p2/p1) / (z2-z1)
        B = p1 / exp(A*z1)

        p_avg(l_indx) = dp / log(p2/p1)

!       ------------------------------------------------
!       LAYER thickness (rather long-winded as it is not
!       assumed the layers are thin) in m. Includes
!       correction for altitude/gravity.
!       ------------------------------------------------

!       ...Initial values
        g_avg = g(l)
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

!       ...Calculate z_avg
        z_avg = z(l) + ( 0.5*dz )

!       ...Calculate new g_avg
        r_hgt = r_avg + z_avg 
        g_avg = g_sfc - g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

!       ...Calculate new dz
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

!       ----------------------------------------
!       Calculate LAYER amounts for water vapour
!       ----------------------------------------

        w1 = w_ppmv(l)
        w2 = w_ppmv(l-i_dir)

        w_avg =  ( (rho1*w1) + (rho2*w2) ) / ( rho1+rho2 )

        w_amt(l_indx) = rho_ref * w_avg * dz * p_avg(l_indx) / t_avg(l_indx)

!       ---------------------------------
!       Calculate LAYER amounts for ozone
!       ---------------------------------

        o1 = o(l)
        o2 = o(l-i_dir)

        o_avg =  ( (rho1*o1) + (rho2*o2) ) / ( rho1+rho2 )

        o_amt(l_indx) = rho_ref * o_avg * dz * p_avg(l_indx) / t_avg(l_indx)

      end do

      end subroutine conpir
	  


      subroutine gphite( p, t, w, z_sfc, n_levels, i_dir, z)
! .... version of 18.05.00

! PURPOSE:

!  Routine to compute geopotential height given the atmospheric state.
!    Includes virtual temperature adjustment.

! CREATED:

! 19-Sep-1996 Received from Hal Woolf, recoded by Paul van Delst
! 18-May-2000 Logic error related to z_sfc corrected by Hal Woolf

! ARGUMENTS:

!     Input
!    --------
!       p     - REAL*4 pressure array (mb)

!       t     - REAL*4 temperature profile array (K)

!       w     - REAL*4 water vapour profile array (g/kg)

!     z_sfc   - REAL*4 surface height (m).  0.0 if not known.

!    n_levels - INT*4 number of elements used in passed arrays

!     i_dir   - INT*4 direction of increasing layer number

!                 i_dir = +1, Level(1) == p(top)         } satellite/AC
!                             Level(n_levels) == p(sfc)  }    case

!                 i_dir = -1, Level(1) == p(sfc)         } ground-based
!                             Level(n_levels) == p(top)  }    case

!     Output
!     --------
!        z     - REAL*4 pressure level height array (m)

!  COMMENTS:

!    Dimension of height array may not not be the same as that of the
!      input profile data.

!=======================================================================

!-----------------------------------------------------------------------
!              -- Prevent implicit typing of variables --
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!                           -- Arguments --
!-----------------------------------------------------------------------

! -- Arrays

      real*4  ::  p(*), t(*), w(*),  z(*)

! -- Scalars

      integer*4 :: n_levels, i_dir

      real*4  ::  z_sfc

!-----------------------------------------------------------------------
!                         -- Local variables --
!-----------------------------------------------------------------------

! -- Parameters

      real*4, parameter :: rog = 29.2898, &
                 fac = 0.5 * rog 

! -- Scalars

      integer*4 :: i_start, i_end, l

      real*4 ::   v_lower, v_upper, algp_lower, algp_upper, hgt

!***********************************************************************
!                         ** Executable code **
!***********************************************************************

!-----------------------------------------------------------------------
!  -- Calculate virtual temperature adjustment and exponential       --
!  -- pressure height for level above surface.  Also set integration --
!  -- loop bounds                                                    --
!-----------------------------------------------------------------------

      if( i_dir > 0 )then

!       --------------------
!       Data stored top down
!       --------------------

        v_lower = t(n_levels) * ( 1.0 + ( 0.00061 * w(n_levels) ) )

        algp_lower = alog( p(n_levels) )

        i_start = n_levels-1
        i_end   = 1

      else

!       ---------------------
!       Data stored bottom up
!       ---------------------

        v_lower = t(1) * ( 1.0 + ( 0.00061 * w(1) ) )

        algp_lower = alog( p(1) )

        i_start = 2
        i_end   = n_levels

      end if

!-----------------------------------------------------------------------
!                     -- Assign surface height --
!-----------------------------------------------------------------------

      hgt = z_sfc

! .. Following added 18 May 2000 ... previously, z(n_levels) for downward
!       (usual) case was not defined!

      if(i_dir.gt.0) then
         z(n_levels) = z_sfc
      else
         z(1) = z_sfc
      endif

! .. End of addition

!-----------------------------------------------------------------------
!             -- Loop over layers always from sfc -> top --
!-----------------------------------------------------------------------

      do l = i_start, i_end, -1*i_dir

!       ----------------------------------------------------
!       Apply virtual temperature adjustment for upper level
!       ----------------------------------------------------

        v_upper = t(l)
        if( p(l) >= 300.0 ) v_upper = v_upper * ( 1.0 + ( 0.00061 * w(l) ) )

!       ----------------------------------------------------- 
!       Calculate exponential pressure height for upper layer
!       ----------------------------------------------------- 

        algp_upper = alog( p(l) )

!       ----------------
!       Calculate height
!       ----------------

        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )

!       -------------------------------
!       Overwrite values for next layer
!       -------------------------------

        v_lower = v_upper
        algp_lower = algp_upper

!       ---------------------------------------------
!       Store heights in same direction as other data
!       ---------------------------------------------

        z(l) = hgt

      end do


      end subroutine gphite


	subroutine taudoc(nc,nx,ny,cc,xx,tau)
! * Strow-Woolf model ... for dry, ozo(ne), and wco (water-vapor continuum)
! .... version of 05.09.02
	integer, intent(in) :: nc, ny, nx
	real, intent(in) :: cc(:,:), xx(:,:)
	real, intent(inout) :: tau(*)

!   cc(11,100)
!	dimension cc(nc,ny),xx(nx,ny),tau(*)

	real, parameter :: trap = -999.99
	
	integer :: last, j, i
	real :: taulyr, taulev, yy
	
	last=ny
	taulyr=1.
	taulev=1.
	
	tau(1)=taulev

	do j=1,ny
		
		if (taulev == 0.) then 
			tau(j+1) = taulev
			exit
		endif
	
		if(j > last) then 
			taulev=taulev*taulyr
			tau(j+1)=taulev
			exit
		endif
			
		yy=cc(nc,j)
	
		if(yy == trap) then

			last=j
			taulev=taulev*taulyr
			tau(j+1)=taulev
			exit 
			
		endif

		do i=1,nx
			yy=yy+cc(i,j)*xx(i,j)
		enddo
	
		yy=amax1(yy,0.)
		if(yy > 0.) taulyr=exp(-yy)
			
		taulev=taulev*taulyr
		tau(j+1)=taulev
		
	end do
	
	end subroutine taudoc


	subroutine taudry(nc,nx,ny,cc,xx,tau)
! * Strow-Woolf model ... for dry, ozo(ne), and wco (water-vapor continuum)
! .... version of 05.09.02
	integer, intent(in) :: nc, ny, nx
	real, intent(in) :: cc(:,:), xx(:,:)
	real, intent(inout) :: tau(:)

!   cc(11,100)
!	dimension cc(nc,ny),xx(nx,ny),tau(*)

	real, parameter :: trap = -999.99
	
	integer :: j, i
	real :: taulyr, taulev, yy
	
	taulyr=1.
	taulev=1.	
	tau(1)=taulev

	do j=1,ny
		
		if (taulev == 0.) then 
			tau(j+1) = taulev
			exit
		endif
	
		yy=cc(nc,j)
	
		if(yy == trap) then

			taulev=0.
			tau(j+1)=taulev
			exit 
			
		endif

		do i=1,nx
			yy=yy+cc(i,j)*xx(i,j)
		enddo
	
		yy=amax1(yy,0.)
		if(yy > 0.) taulyr=exp(-yy)
			
		taulev=taulev*taulyr
		tau(j+1)=taulev
		
	end do
	
	end subroutine taudry




	subroutine tauwtr(ncs,ncl,nxs,nxl,nxw,ny,ccs,ccl,xx,tau)
! * Strow-Woolf model ... for 'wet' (water-vapor other than continuum)
! .... version of 05.09.02

		integer, intent(in) :: ncs, ncl, nxs, nxl, nxw, ny
		real, intent(in) :: ccs(:,:), ccl(:,:), xx(:,:)
		real, intent(inout) :: tau(:)

!	dimension ccs(ncs,ny),ccl(ncl,ny),xx(nxw,ny),tau(*)

		real :: odsum, taulyr, taulev, yy
		integer :: j, i 

	odsum=0.
	taulyr=1.
	taulev=1.
	tau(1)=taulev

	do  j=1,ny
		if(odsum < 5.) then
			yy=ccs(ncs,j)
			do i=1,nxs
				yy=yy+ccs(i,j)*xx(i,j)
			enddo
		else
			yy=ccl(ncl,j)
			do i=1,nxl
				yy=yy+ccl(i,j)*xx(i+11,j)
			enddo
		endif
		yy=amax1(yy,0.)
		odsum=odsum+yy
		if(yy > 0.) taulyr=exp(-yy)
		taulev=taulev*taulyr
		tau(j+1)=taulev
	
	end do

	end subroutine tauwtr







end module FASCODE_routines

