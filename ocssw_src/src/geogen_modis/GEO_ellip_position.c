#include "PGS_CSC.h"
#include "PGS_MODIS_35251.h"
#include "GEO_earth.h"
#include "GEO_global_arrays.h"
#include "GEO_util.h"
#include "smfio.h"

int GEO_ellip_position(
	int const       scan_number,
	int const       sample_number,
	int const       num_detectors, 
	double          u_inst[][3], 
        double          ecr_sc_sample_position[][3],
        double          ecr_sc_velocity[][3],
        double          T_inst2ecr[][3][3],
        double          ecr_sample_position[][MAX_PADDED][3],
        double          ellip_sample_position[][MAX_PADDED][3],
        unsigned char   sample_flags[][MAX_PADDED],
        double          sample_view_vec[][MAX_PADDED][3]
	)
/*
!C*****************************************************************************
!Description:   
        Routine in Earth Location group of the Level-1A geolocation software to
        calculate the ellipsoidal position of a pixel. It transforms each view
        vector to the ECR frame and solves for the magnitude of the vector
        (See Geolocation ATBD section 3.1.3.2).  Finally, it computes the
        geodetic latitude and longitude for the pixel and stores them in the
        global array.  If the view vector for a detector does not intersect the
        earth ellipsoid, no ecr_sample_position nor ellip_sample_position will
        be determined, and the sample_view_vec will not be rescaled for that
        detector.
 
        GEO_ellip_position initializes the terrain sample position validity
        sample_flags, setting NO_ELLIPSE_INTERSECT for detectors with view
        vectors that do not intersect the earth ellipsoid.  The
        INVALID_INPUT_DATA flag is not set if a routine failure occurs, it
        being set for each detector in the frame in GEO_locate_one_scan. 

!Input Parameters:
        scan_number             The scan number
        sample_number           The sample number within the scan
        num_detectors           Number of detectors for the sample
        u_inst                  View vectors for each detector in the
                                instrument coordinate system (meters).
        ecr_sc_sample_position  The spacecraft position.
        ecr_sc_velocity         The spacecraft velocity.
        T_inst2ecr              Matrix for conversion from instrument to ECR
                                coordinates.
 
!Output Parameters:
        ecr_sample_position     Ground positions in ECR coordinates.
        ellip_sample_position   Ground positions without terrain correction.
        sample_flags            Error flags for each sample.
        sample_view_vec         View vector from satellite to ground position.
       
Return Values:
        FAIL            If arguments are invalid or PGS_CSC_GetEarthFigure()
                        fails.
        SUCCESS         Otherwise 

Externally Defined: 
        DETECTORS_QKM           "GEO_geo.h" 
        FAIL                    "GEO_basic.h"
        INVALID_INPUT_DATA      "GEO_geo.h"
        MAX_SCAN_NUMBER         "GEO_geo.h"
        MODIS_E_BAD_INPUT_ARG   "PGS_MODIS_35251.h"
        MODIS_E_BAD_VIEW_VEC    "PGS_MODIS_35251.h"
        MODIS_E_GEO             "PGS_MODIS_35251.h"
        MODIS_W_SCAN_OFF        "PGS_MODIS_35251.h"
        NO_ELLIPSE_INTERSECT    "GEO_geo.h" 
        MAX_PADDED              "GEO_geo.h"
        PGSd_GEO_ERROR_VALUE    "PGS_TD.h"
        PGS_S_SUCCESS           "PGS_SMF.h"
        SUCCESS                 "GEO_basic.h"    

Called by: GEO_earth_location
 
Routines called:
        PGS_CSC_ECRtoGEO        "PGS_CSC.h"
        PGS_CSC_GetEarthFigure  "PGS_CSC.h"
        modsmf                  "smfio.h"

!Revision History:
		$Log: GEO_ellip_position.c,v $
		Revision 6.4  2010/06/29 18:45:31  kuyper
		Replaced floating point equality test with an inequality.

		Revision 6.3  2009/05/31 19:02:17  kuyper
		Corrected to initialize p_prime_sq, up_prod, u_prime_sq where needed.
		Improved error messages.

		Revision 6.2  2009/05/15 17:00:18  xgeng
		Minor correction on error messages.

                Xu Geng (xu.geng@saic.com)
		Revision 6.1  2009/05/06 19:34:43  xgeng
		Code updated to match the PDL change made in revision 6.1, 6.2 and 6.3.

		Revision 4.4  2003/09/03 17:17:04  vlin
		corrected a typo.

		Revision 4.3  2003/08/28 15:56:51  kuyper
		Corrected prolog.

		Revision 4.2  2003/08/22 14:24:21  kuyper
		Removed T_inst2ecr and sc_attitude globals.

		Revision 4.1  2003/08/08 19:36:28  vlin
		Calls to GEO_interp_ephemeris_attitude() and GEO_get_T_inst2ecr() are moved
		  up to a higher level. As a result, ecr_sc_position, ecr_sc_velocity,
		  scan_start_time, and T_inst2ecr are now inputs rather than outputs.
		Removed intitialization of sc_attitude.
		Shifted messages to MODIS_E_GEO, where applicable.

		Revision 2.4  1999/01/29 20:14:01  kuyper
		Removed max_extrap from GEO_interp_ephemeris() argument list.

                4/25/95
                Ruiming Chen 
                Finished coding.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

Requirements:
		PR03-F-3.2.2-1
		PR03-F-3.2.3-1
		PR03-F-3.2.3-2
		PR03-F-3.2.4-1
		PR03-I-1
		PR03-I-2

!END*************************************************************************
*/

{
  static PGSt_double equat=-1.0; /* equatorial radius of the earth 	*/
  static PGSt_double polar=-1.0; /* polar radius of the earth.		*/
  PGSt_double p_prime[3];	/* rescaled sc position */
  double p_prime_sq = 0.0;	/* square of p_prime_length */
  PGSt_SMF_status PGS_error_code;	/* SDP toolkit function return value */
  int det;	                /* detector index */
  int ecr;			/* iteration parameter */
  char msgbuf[128] = "";
  static char filefunc[] = __FILE__ ", GEO_ellip_position";

/****************************************************************************
calculate from here
****************************************************************************/

  if( scan_number < 0 || scan_number >= MAX_SCAN_NUMBER ||
      num_detectors < 0 || num_detectors > DETECTORS_QKM ||
      sample_number < 0 || sample_number >= MAX_PADDED ||
      u_inst == NULL || ecr_sc_sample_position == NULL || 
      ecr_sc_velocity == NULL || T_inst2ecr == NULL ||
      ecr_sample_position == NULL || ellip_sample_position == NULL ||
      sample_flags == NULL || sample_view_vec == NULL )
  {
    sprintf(msgbuf, "scan:%d detectors:%d sample:%d u_inst:%p " 
                    "ecr_sc_sample_position:%p ecr_sc_velocity:%p " 
                    "T_inst2ecr:%p ecr_sample_position:%p " 
                    "ellip_sample_position:%p sample_flags:%p " 
                    "sample_view_vec:%p",
	scan_number, num_detectors, sample_number, (void *)u_inst, 
        (void *)ecr_sc_sample_position, (void *)ecr_sc_velocity, 
        (void *)T_inst2ecr,  (void *)ecr_sample_position, 
        (void *)ellip_sample_position,  (void *)sample_flags, 
        (void *)sample_view_vec   
    );

    modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
    return FAIL;
  }

  if (ecr_sc_sample_position[sample_number][0] >= PGSd_GEO_ERROR_VALUE) {
      for (det = 0; det < num_detectors; det++)
           sample_flags[det][sample_number] = INVALID_INPUT_DATA;
      return SUCCESS;
  }

  for (det=0; det<num_detectors ; det++)
       sample_flags[det][sample_number] = 0;
  
  if (equat < 0.0) { 
  /* the equatorial and polar earth radii have not been set */
      PGS_error_code = PGS_CSC_GetEarthFigure(EARTH_MODEL, &equat, &polar);
      if (PGS_error_code!=PGS_S_SUCCESS)
      {
          modsmf(MODIS_E_GEO, "PGS_CSC_GetEarthFigure(\"" EARTH_MODEL "\")",
	      filefunc);
          return FAIL;
      }
  }

  /* Rescale orbit position vector */
  p_prime[0] = ecr_sc_sample_position[sample_number][0] / equat;	
  p_prime[1] = ecr_sc_sample_position[sample_number][1] / equat;
  p_prime[2] = ecr_sc_sample_position[sample_number][2] / polar;

  for (ecr = 0; ecr < 3; ecr ++) 
    p_prime_sq += p_prime[ecr]*p_prime[ecr];

  /* Loop through detectors in sample */
  for (det = 0; det < num_detectors; det++)
  {
    PGSt_double u_prime[3];	/* rescaled view vec */
    double up_prod_sq;		/* square of up_prod */
    double sqrt_arg;	 /* argument of square root in the calculation of d */
    double d;			/* scaling factor of the viewing vector */
    double up_prod = 0.0;	/* product between u_prime and p_prime */
    double u_prime_sq = 0.0;	/* square of u_prime_length */

    /* Transform the viewing vector to ECR */ 
    for (ecr = 0; ecr < 3; ecr ++) {  
      int inst;			/* iteration parameter */
      sample_view_vec[det][sample_number][ecr] = 0.0;

      for (inst = 0; inst < 3; inst ++) {
        sample_view_vec[det][sample_number][ecr] += 
          T_inst2ecr[sample_number][ecr][inst]*u_inst[det][inst]; 
      }
    } 
    
    /* Rescale viewing vector */
    u_prime[0] = sample_view_vec[det][sample_number][0] / equat;
    u_prime[1] = sample_view_vec[det][sample_number][1] / equat;
    u_prime[2] = sample_view_vec[det][sample_number][2] / polar;

    /* Compute dot product of u_prime and p_prime  
       and squared magnitudes of u_prime and p_prime */
    for (ecr = 0; ecr < 3; ecr ++) { 
	up_prod += u_prime[ecr]*p_prime[ecr]; 
	u_prime_sq += u_prime[ecr]*u_prime[ecr];
    } 

    up_prod_sq = up_prod * up_prod; 

    /* Check for negative under square root (scan off of Earth) */
    sqrt_arg = up_prod_sq - u_prime_sq*(p_prime_sq - 1.0);
    if (sqrt_arg < 0.0 || up_prod>=0.0) { /* No real ellipsoid intersection*/
      /* call SDP function to report event */
      modsmf(MODIS_W_SCAN_OFF, "", filefunc);
      sample_flags[det][sample_number] |= NO_ELLIPSE_INTERSECT;
    }
    else if ( (-up_prod - sqrt(sqrt_arg))/DBL_MAX >=  u_prime_sq)
    {	/* Calculation of d will produce overflow.	*/
      modsmf(MODIS_E_BAD_VIEW_VEC, "", filefunc);
      sample_flags[det][sample_number] |= INVALID_INPUT_DATA;
    }
    else
    {
      PGSt_double latitude;	/* latitude of the ellipsoid position	*/
      PGSt_double longitude;	/* longitude of the ellipsoid position	*/
      PGSt_double altitude;	/* altitude of the ellipsoid position	*/

      d = (-up_prod - sqrt(sqrt_arg)) / u_prime_sq;

      /* calculate the geocentric position vector */
      for (ecr = 0; ecr < 3; ++ecr) {
	sample_view_vec[det][sample_number][ecr] *= d; 
	ecr_sample_position[det][sample_number][ecr] = 
	  ecr_sc_sample_position[sample_number][ecr] +
	  sample_view_vec[det][sample_number][ecr];
      }

      /* calculate the geodetic lat and lon */

      if (PGS_CSC_ECRtoGEO( ecr_sample_position[det][sample_number],
	  EARTH_MODEL, &longitude, &latitude, &altitude) != PGS_S_SUCCESS)
      {
        modsmf(MODIS_E_GEO, "PGS_CSC_ECRtoGEO(\"" EARTH_MODEL "\")", filefunc);
	sample_flags[det][sample_number] |= INVALID_INPUT_DATA;
	continue;
      } 
      
      /* store latitude and longitude in global variable */
      ellip_sample_position[det][sample_number][LAT_IDX] = latitude;
      ellip_sample_position[det][sample_number][LON_IDX] = longitude;
      ellip_sample_position[det][sample_number][HT_IDX] = 0.0;
    }
  }

  return SUCCESS;
} 

