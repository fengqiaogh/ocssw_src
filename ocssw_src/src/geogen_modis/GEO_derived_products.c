#include "PGS_CBP.h" 
#include "PGS_MODIS_35251.h"
#include "GEO_output.h"
#include "GEO_global_arrays.h" 
#include "GEO_util.h" 
#include "GEO_geo.h"

PGSt_SMF_status GEO_derived_products(
        int		const num_frames,
        int		const num_detectors,
        double		ecr_sample_position[][MAX_SCAN_SAMPLE][3],
        double		ecr_sc_position[][3],
        unsigned char	pixel_flags[][MAX_SCAN_SAMPLE],
        double		terrain_sample_position[][MAX_SCAN_SAMPLE][3],
        char		utc[],
        PGSt_double	toff[],
        double		sample_solar_angles[][MAX_SCAN_SAMPLE][2],
        double		sample_to_sensor[][MAX_SCAN_SAMPLE][3]
)
/*
!C*****************************************************************************
!Description:   
        Subroutine in output group of the Level-1A geolocation software to
        compute the derived products, consisting of the sensor azimuth and
        zenith angles, the pixel-to-sensor range, and the solar azimuth and
        zenith angles.

!Input Parameters:
        num_frames              The number of frames to process
        num_detectors           The number of detectors in the sample
        ecr_sample_position     The ellipsoid intersection points in ECR
                                    coordinates.
        ecr_sc_position         The spacecraft positions in ECR coordinates.
        terrain_sample_position Terrain_corrected ground positions
        utc                     The time of the start of a scan, as a character
                                    string in UTC format.
        toff                    Offsets from 'utc'

!Output Parameters:
        sample_solar_angles     Solar azimuth and zenith angles
        sample_to_sensor        Sensor azimuth and zenith angles and range

!Input/Output Parameters:
        pixel_flags             Pixel flags.

Return value:
        MODIS_E_GEO             If any subroutine fails.
        MODIS_E_BAD_INPUT_ARG   Invalid input argument
        PGS_S_SUCCESS           Otherwise

Externally Defined:
        azimuth_index                   "GEO_global_arrays.h"
        MAX_DETECTORS                   "GEO_geo.h"
        MAX_SCAN_SAMPLE                 "GEO_geo.h"
        MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
        MODIS_E_GEO                     "PGS_MODIS_35251.h"
        NO_ELLIPSE_INTERSECT            "GEO_geo.h"
        PGS_S_SUCCESS                   "PGS_SMF.h"
        PGSCBP_W_BAD_CB_VECTOR          "PGS_CBP_6.h"
        PGSCSC_W_PREDICTED_UT1          "PGS_CSC.h"
        range_index                     "GEO_global_arrays.h"
        SENSORAZIM_FVALUE               "GEO_geo.h"
        SENSORZEN_FVALUE                "GEO_geo.h"
        SOLARAZIM_FVALUE                "GEO_geo.h"
        SOLARZEN_FVALUE                 "GEO_geo.h"
        z_angle_index                   "GEO_global_arrays.h"

Called by:
        GEO_locate_one_scan

Routines Called:
        GEO_vec_mul3                    "GEO_utils.h"
        GEO_vec_length3                 "GEO_utils.h"
        GEO_vec_prod3                   "GEO_utils.h"
        modsmf                          "smfio.h"
        PGS_CBP_Earth_CB_Vector         "PGS_CBP.h"
        PGS_CSC_ECItoECR                "PGS_CSC.h"

!Revision History:
	$Log: GEO_derived_products.c,v $
	Revision 6.2  2010/06/29 18:44:37  kuyper
	Replaced equality test involving floating point numbers with an inequality.

	Revision 6.1  2010/06/18 20:08:07  kuyper
	Removed const qualifiers from pointer to array arguments.

	Revision 4.5  2003/10/27 01:01:57  vlin
	expanded message buffer size

	Revision 4.4  2003/09/24 18:46:49  vlin
	Added some lines to check input arguments

	Revision 4.3  2003/08/26 21:06:25  kuyper
	Corrected obsolete setting of toff.

	Revision 4.2  2003/08/18 20:33:52  kuyper
	Move two globals to GEO_locate_one_scan.c.

	Revision 4.1  2003/08/12 16:09:40  vlin
	Changed to use parameters rather than global variables,
	process entire scan in one call, and allowing greater efficiency
	in the SDP Toolkit function calls.

	Revision 1.11  1997/10/22 20:21:50  jjb
	Added initialization of (global) outputs.

	10/3/95
	Frederick S. Patt
	Revised calculation of range to account for terrain correction

!Team-unique Header:

        This software is developed by the MODIS Science Data Support
        Team for the National Aeronautics and Space Administration,
        Goddard Space Flight Center, under contract NAS5-32373.

!END
******************************************************************************/

{
  int    sample_number, idet = 0, j = 0; 
  double clon = 0.0;   /* cosine of longitude */
  double clat = 0.0;   /* cosine of latitude */
  double slon = 0.0;   /* sine of longitude */
  double slat = 0.0;   /* sine of latitude */
  double n[3] = {0.0}; /* normal vector */
  double e[3] = {0.0}; /* East vector */
  double a[3] = {0.0}; /* vector perpendicular to n and e */
  double v[3] = {0.0}; /* pixel-to-sensor vector */
  double vdotn = 0.0;  /* Dot product of vector with n */
  double vdote = 0.0;  /* Dot product of vector with e */
  double vdota = 0.0;  /* Dot product of vector with a */
  double sunecr[3] = {0.0}; /* Sun vector in ECR frame */
  PGSt_double sun[MAX_SCAN_SAMPLE][3] = {0.0};  /* Sun vector in ECI frame */
  PGSt_double tveci[MAX_SCAN_SAMPLE][6] = {0.0};
  PGSt_double tvecr[MAX_SCAN_SAMPLE][6] = {0.0};
  PGSt_SMF_status PGS_error_code;
  char msgbuf[256] = "";
  char filefunc[] = __FILE__ ", GEO_derived_products";

  if (num_frames < 0    || num_frames > MAX_SCAN_SAMPLE || 
      num_detectors < 0 || num_detectors > MAX_DETECTORS ||
      utc == NULL || toff == NULL) {
      sprintf(msgbuf, "num_frames: %d, num_detectors: %d,\nutc: %s, toff:%p",
              num_frames, num_detectors, utc, (void *)toff);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  if (ecr_sample_position == NULL || ecr_sc_position == NULL ||
      pixel_flags == NULL || terrain_sample_position == NULL ||
      sample_solar_angles == NULL || sample_to_sensor == NULL) {
      sprintf(msgbuf, "ecr_sample_position: %p, ecr_sc_position:%p, \n "
        "pixel_flags: %p, terrain_sample_position: %p, \n"
        "sample_solar_angles: %p, sample_to_sensor: %p", 
        (void *)ecr_sample_position,(void *)ecr_sc_position,(void *)pixel_flags,
        (void *)terrain_sample_position,(void *)sample_solar_angles,
        (void *)sample_to_sensor);
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      return MODIS_E_BAD_INPUT_ARG;
  }

  for (idet = 0; idet < num_detectors; idet++) { 
  for (sample_number = 0; sample_number < num_frames; sample_number++) {
    sample_to_sensor[idet][sample_number][z_angle_index] = 
		(double) SENSORZEN_FVALUE;
    sample_to_sensor[idet][sample_number][azimuth_index] = 
		(double) SENSORAZIM_FVALUE;
    sample_to_sensor[idet][sample_number][range_index] = (double) RANGE_FVALUE;
    sample_solar_angles[idet][sample_number][z_angle_index] = 
		(double) SOLARZEN_FVALUE;
    sample_solar_angles[idet][sample_number][azimuth_index] = 
		(double) SOLARAZIM_FVALUE;
  }
  }

  /* Get Sun vector */
  PGS_error_code = PGS_CBP_Earth_CB_Vector(num_frames,utc,toff,PGSd_SUN,sun);
  if (PGS_error_code != PGS_S_SUCCESS && 
      PGS_error_code != PGSCBP_W_BAD_CB_VECTOR) {
      sprintf(msgbuf,"PGS_CBP_Earth_CB_Vector(%s)", utc);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
  }
  
  for (sample_number = 0; sample_number < num_frames; sample_number++)
  for (j=0; j<3; j++)
      tveci[sample_number][j] = sun[sample_number][j];

  PGS_error_code = PGS_CSC_ECItoECR(num_frames,utc,toff,tveci,tvecr);
  if (PGS_error_code != PGS_S_SUCCESS &&
      PGS_error_code != PGSCSC_W_PREDICTED_UT1) {
      sprintf(msgbuf,"PGS_CSC_ECItoECR(%s)", utc);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);
      return MODIS_E_GEO;
  }

  for (sample_number = 0; sample_number < num_frames; sample_number++) {
    for (j=0; j<3; j++)
        sunecr[j] = (double)tvecr[sample_number][j];

    for (idet = 0; idet < num_detectors; idet++) { 
      if (pixel_flags[idet][sample_number] < NO_ELLIPSE_INTERSECT) {
      /* Compute the normal and East vectors */
      clon = cos(terrain_sample_position[idet][sample_number][1]);
      clat = cos(terrain_sample_position[idet][sample_number][0]);
      slon = sin(terrain_sample_position[idet][sample_number][1]);
      slat = sin(terrain_sample_position[idet][sample_number][0]);
      n[0] = clat*clon;
      n[1] = clat*slon;
      n[2] = slat;
      e[0] = -slon;
      e[1] = clon;
      
      /* Compute the perpendicular vector */
      GEO_vec_mul3(n, e, a);

      /* Compute the sample-to-sensor vector and store in v */
      for (j = 0; j < 3; j++)
	v[j] = ecr_sc_position[sample_number][j] 
	  - ecr_sample_position[idet][sample_number][j];

      /* Compute the sensor azimuth and zenith angles and range */
      GEO_vec_prod3(v, n, &vdotn);
      GEO_vec_prod3(v, e, &vdote);
      GEO_vec_prod3(v, a, &vdota);
      sample_to_sensor[idet][sample_number][z_angle_index] =
      atan2(sqrt(vdote*vdote + vdota*vdota),vdotn); 
      sample_to_sensor[idet][sample_number][azimuth_index] = atan2(vdote,vdota);
      GEO_vec_length3(v,&sample_to_sensor[idet][sample_number][range_index]);
    
      if (sunecr[0] >= PGSd_GEO_ERROR_VALUE)
          pixel_flags[idet][sample_number] |= INVALID_SOLAR_ANGLES;
      else { 
          GEO_vec_prod3(sunecr, n, &vdotn);
          GEO_vec_prod3(sunecr, e, &vdote);
          GEO_vec_prod3(sunecr, a, &vdota);
          sample_solar_angles[idet][sample_number][1] = 
                 atan2(sqrt(vdote*vdote + vdota*vdota),vdotn);
          sample_solar_angles[idet][sample_number][0] = atan2(vdote,vdota);
      }
      }
    } /* End of detector loop */
  } 
  return PGS_S_SUCCESS;
}
