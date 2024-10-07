#include "PGS_SMF.h"
#include "PGS_MODIS_35251.h"
#include "GEO_output.h"
#include "GEO_earth.h"
#include "GEO_util.h"
#include "GEO_validation.h"

PGSt_SMF_status GEO_write_one_scan( 
	PGSt_double const	frame_time[],
	int			const scan_number,
	l1a_data_struct const	* const l1a_data,
	GEO_param_struct const	* const geo_params,
	double			frame_to_sensor[][MAX_FRAMES][3],
	double			frame_solar_angles[][MAX_FRAMES][2],
	double			terrain_frame_position[][MAX_FRAMES][3],
	uint8			frame_flags[][MAX_FRAMES],
	int8			hires_offsets[][DETECTORS_QKM][SAMPLES_QKM],
	PGSt_double const	rpy[],
	uint32			frame_quality[][MAX_FRAMES],
	uint8			frame_landsea[][MAX_FRAMES],
	uint8			frame_waterpresent[][MAX_FRAMES],
	uint8			land_seamask_qaflag,
	MODFILE			* const geo_file
)
/*
!C*****************************************************************************
!Description:   
		subroutine in output group of the Level-1A geolocation
                software to write scan-level data to the output product.
		It calls routines to write the scan data and the scan
                          l1a_data->frame_data[scan_number].EV_frames,
		metadata to the geolocation file.

!Input Parameters:
	scan_number		the scan number
        sc_ev_frame_state	array of spacecraft kinematic states for each
				    Earth View frame.
        l1a_data		data read from L1A file
        geo_params		geolocation parameters read from input file
        maneuver_list		A list of spacecraft maneuvers
	frame_to_sensor		pixel-to-satellite angles and range
	frame_solar_angles	solar azimuth and zenith angles 
        terrain_frame_position	terrain-corrected geodetic coordinated for
				    pixels
        frame_flags		pixel flags
        hires_offsets		High resolutions offsets from bilinearly
				    interpolated low resolution data.
	rpy			Roll, pitch, and yaw components of thermal
				    attitude correction.
	frame_quality		Ephemeris and Attitude quality flags.
	frame_landsea		Land/Water mask.
	frame_waterpresent	Sum of weights of high-resolution pixels with
				    water aggregated into low-resolution mask.
	land_seamask_qaflag	Indicates whether there were problems with the
				Land/Water mask.
	geo_file		M-API structure for the geolocation file

!Output Parameters:
        None

Return parameter:
        MODIS_E_BAD_INPUT_ARG   If input values are found to be invalid.
        MODIS_E_GEO             If any subroutine fails
        PGS_S_SUCCESS           Otherwise

Externally Defined:
        GOOD_DATA                       "GEO_validation.h"
        MAX_FRAMES                      "GEO_geo.h"
        MAX_SCAN_NUMBER                 "GEO_geo.h"
        MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
        MODIS_E_GEO                     "PGS_MODIS_35251.h"
        PGS_S_SUCCESS                   "PGS_SMF.h"
        SAMPLES_QKM                     "GEO_geo.h"

Called by:
	GEO_locate_one_scan

Routines called:
        GEO_solar_and_lunar_vectors - spacecraft, solar, and lunar
                                      vector metadata.
	GEO_write_scan_data         - write scan data
	GEO_write_scan_metadata     - write scan-level data
        modsmf                      - writes error status messages to log
		
!Revision History:
 * $Log: GEO_write_one_scan.c,v $
 * Revision 6.3  2011/02/18 21:58:21  kuyper
 * In order to resolve feature-request bug 3446, changed to receive landsea
 *    mask, rather than calling function to fill it in.
 *
 * Revision 6.2  2010/05/27 15:28:49  kuyper
 * Helped resolve Bug 1969 by adding rpy as an input parameters to be passed to
 *   GEO_write_scan_metadata().
 * Helped resolve Bug 2470 by removing a parameter.
 * Helped resolve Bug 2472 by adding frame_quality as an input parameters to be
 *   passed to GEO_write_scan_metadata().
 * Changed to expect status codes from GEO_write_scan_metadata().
 *
 * Revision 6.1  2009/05/27 14:47:05  xgeng
 * Change to return a status code, and to expect one from GEO_write_scan_data().
 * Pass frame_time on to GEO_solar_and_lunar_vectors().
 * Move code for rearranging arrays to GEO_write_scan_data, where it can be
 *   combined with code that changes the value and formate of the data, avoiding
 *   need to keep a copy in scan_data.
 * Rename pixel_flags as frame flags.
 * Change to use DETECTORS_1KM rather than num_detectors for all loop limits.
 * Changed MAX_SCAN_SAMPLE to MAX_FRAMES.
 * Added frame_time, hires_offsets and terrain_frame_position as input parameters.
 * Corrected frame_time's data type, removed sc_ev_frame_state.
 * The above changes match with PDL revision 6.1 to 6.3.
 * 
 * Xu Geng (xu.geng@saic.com)
 * 
 * Revision 5.2  2004/10/27 19:55:00  vlin
 * call to function GEO_solar_and_lunar_vectors updated.
 * vlin@saicmodis.com
 *
 * Revision 5.1  2004/08/25 19:06:59  vlin
 * Input parameters maneuver_list, sample_to_sensor, and
 *       sample_solar_angles added.
 * Change to check for  !=SUCCESS instead of  ==FAIL.
 *
 * Revision 4.3  2004/04/09 21:01:16  kuyper
 * Corrected typo.
 *
 * Revision 4.2  2003/12/30 21:21:43  kuyper
 * Changed message for failure of GEO_solar_and_lunar_vectors.
 * Added file information to message for GEO_write_scan_data call.
 *
 * Revision 4.1  2003/02/21 22:56:11  kuyper
 * Corrected to use void* pointers with %p format code.

Requirements:
                PR03-F-3.5.1-1
                PR03-F-3.5.2-1
                PR03-F-3.5.3-1
                PR03-F-3.5.3-2
                PR03-F-3.6.1-1
                PR03-F-3.6.2-1
                PR03-F-3.6.3-2
                PR03-F-4.3-1
                PR03-F-4.3-2
                PR03-F-4.4-1
                PR03-F-4.4-2
                CCR-319

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
  swath_elem_struct	scan_data;
  celestial_bodies_struct cb_vectors;
  /* offsets into sensor, solar arrays. */
  PGSt_SMF_status retval=PGS_S_SUCCESS;	/* Return value for this routine. */
  unsigned char flag;	/* qaflag value for scan.	*/
  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE] = "";
  char *filefunc = __FILE__ ", GEO_write_one_scan";

  /* Begin program logic */
  /* Only need to validate pointers dereferenced in this module. */
  if (geo_file==NULL || l1a_data==NULL ||
      scan_number < 0 || scan_number >= MAX_SCAN_NUMBER) 
  {
      sprintf(msgbuf, "geo_file = %p, l1a_data = %p, scan_number = %d",
	  (void*)geo_file, (void*)l1a_data, scan_number); 
      modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
      
      return MODIS_E_BAD_INPUT_ARG;
  }

  /* compute the celestial body vector metadata of spacecraft, sun, and moon */
  if(GEO_solar_and_lunar_vectors(frame_time, &l1a_data->frame_data[scan_number],
        geo_params, &l1a_data->fill_values, &cb_vectors) != PGS_S_SUCCESS) {
    modsmf(MODIS_E_GEO, "GEO_solar_and_lunar_vectors()", filefunc);
    retval = MODIS_E_GEO;
  }

  /* Write scan metadata */
  if (GEO_write_scan_metadata(scan_number, l1a_data, &cb_vectors, geo_file,
      rpy, frame_quality) != PGS_S_SUCCESS)
  { 
      sprintf(msgbuf, "GEO_write_scan_metadata(%d, \"%s\")", scan_number,
	  geo_file->filename);
      modsmf(MODIS_E_GEO, msgbuf, filefunc);

      retval = MODIS_E_GEO;
  }

  if(l1a_data->frame_data[scan_number].L1A_scan_quality[0] == 0)
    return retval;	/* nothing to write out and no functional error */

  scan_data.num_detectors = DETECTORS_1KM;
  scan_data.num_frames = l1a_data->frame_data[scan_number].EV_frames;

  flag = (l1a_data->mirr_data[scan_number].impulse_flag!=GOOD_DATA ||
      scan_data.num_frames==0);

  /* We have 9 seperate flags, anticipating that in the future they will
    be set in different locations for different reasons. */
  scan_data.lat_qaflag = flag;
  scan_data.lon_qaflag = flag;
  scan_data.height_qaflag = flag;
  scan_data.sensorazimuth_qaflag = flag;
  scan_data.sensorzenith_qaflag = flag;
  scan_data.solarazimuth_qaflag = flag;
  scan_data.solarzenith_qaflag = flag;
  scan_data.range_qaflag = flag;
  scan_data.land_seamask_qaflag = land_seamask_qaflag;

  if(flag==0)
  {
    /* Write geolocation scan data */
    if (GEO_write_scan_data(geo_file, scan_number, &scan_data,
	(GEO_param_struct  *)geo_params, terrain_frame_position,
	frame_to_sensor, frame_solar_angles, frame_flags, frame_landsea,
	frame_waterpresent, hires_offsets)
      != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_write_scan_data(\"%s\",%d)",
	    geo_file->filename, scan_number);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }
  }	/* impulse_flag == GOOD_DATA */

  return retval;
}
