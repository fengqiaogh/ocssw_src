#include "GEO_earth.h"
#include "GEO_global_arrays.h"
#include "GEO_input.h"
#include "GEO_inst.h"
#include "GEO_output.h"
#include "GEO_validation.h"
#include "PGS_MODIS_35251.h"
#include <strings.h>

/* sample time for all the pixels in a scan (relative to sector start time */
double *sample_time;

PGSt_SMF_status GEO_locate_one_scan(
        GEO_param_struct const	* const geo_params,
        l1a_data_struct		* const l1a_data,
        int			const scan_number,
        qa_metadata_struct	* const qa_metadata,
        GEO_bcoord_struct	* const bounding_coords,
        MODFILE			* const geo_file
)
/*****************************************************************************
!C

!Description:
	Routine in Output group of the Level-1A geolocation software to
	geolocate one scan and write to output file.  For each sample in a
	scan, it call routines to:  compute and validate the Earth locations
	for the pixels, and compute and validate the sensor and solar angles
	for the pixels.  It then calls the routine to write the data for the
	entire scan to the geolocation product file.

!Input Parameters:

	geo_params	structure containing geolocation process parameters
			from the parameter file.
	scan_number	the scan number
	geo_file	structure used to reference the geolocation file in 
                        the MAPI library

!Output Parameters:

	qa_metadata	structure for accumulating data quality measures
	bounding_coords	structure of the granules bounding coordinates

Input/Output Parameters:
	l1a_data	data read from L1A file.

Return values:
	MODIS_E_BAD_INPUT_ARG	One of the input arguments had an invalid value
	MODIS_E_GEO		A subroutine failed.
	PGS_S_SUCCESS		Otherwise.

Externally Defined:
	BAD_DATA                        "GEO_validation.h"
	DBL_MIN                         "float.h"
	DETECTORS_1KM                   "GEO_geo.h"
	DETECTORS_QKM                   "GEO_geo.h"
	GEO_DOUBLE_FILLVALUE            "GEO_output.h"
	GOOD_DATA                       "GEO_validation.h"
	INVALID_INPUT_DATA              "GEO_geo.h"
	INVALID_SOLAR_ANGLES            "GEO_geo.h"
	MAX_FRAMES                      "GEO_geo.h"
	MAX_PADDED                      "GEO_geo.h"
	MAX_SCAN_NUMBER                 "GEO_geo.h"
	MODIS_E_BAD_INPUT_ARG           "PGS_MODIS_35251.h"
	MODIS_E_GEO                     "PGS_MODIS_35251.h"
	NO_ELLIPSE_INTERSECT            "GEO_geo.h"
	SAMPLES_QKM                     "GEO_geo.h"
	sample_time                     "GEO_global_arrays.h"
	SUCCESS                         "GEO_basic.h"
	TIMECODEASIZE			"smfio.h"

	Note: sample_time has external linkage, but is to be defined in this
	module.

Called by:
	GEO_locate_one_granule

Routines called:
	GEO_aggregate                   "GEO_earth.h"
	GEO_cumulate_GRing		"GEO_input.h"
	GEO_derived_products		"GEO_output.h"
	GEO_earth_location		"GEO_earth.h"
	GEO_get_bounding_coords		"GEO_output.h"
	GEO_get_sample_time             "GEO_inst.h"
        GEO_interp_ECR                  "GEO_earth.h"
	GEO_validate_derived_products	"GEO_validation.h"
	GEO_validate_earth_location	"GEO_validation.h"
	GEO_write_one_scan		"GEO_output.h"
	modsmf				"smfio.h"

!Revision History:
$Log: GEO_locate_one_scan.c,v $
Revision 6.7  2011/02/18 21:59:34  kuyper
In order to resolve feature-request Bug 3446, changed to retrieve landsea
  mask at high resolution, and then aggregate to low resolution, filling in
  new waterpresent array.

Revision 6.6  2010/06/29 20:11:22  kuyper
Declared many arrays static, solely to reduce strain on stack.

Revision 6.5  2010/06/18 20:11:25  kuyper
Corrected handling of int32 arguments passed to sprintf().

Revision 6.4  2010/05/27 17:58:07  kuyper
Helped resolve Bug 1969 by creating rpy, filling it in by calling
  GEO_interp_ECR(), and passing it to GEO_write_one_scan().
Helped resolve Bug 2470 by removing a parameter.
Helped resolve Bug 2472 by filling in, aggregating, and writing out
  the ephemeris and attitude quality flags, and passing the entire Geolocation
  parameters struct to GEO_interp_ECR().
Resolved Bug 2949 by removing pointless and dangerous code.
Stopped passing num_detectors as an argument to GEO_aggregate().
Changed to return a status code.

James Kuyper		James.R.Kuyper@NASA.gov

Revision 6.3  2009/05/31 22:08:40  kuyper
Corrected to initialize sample_flags.
Corrected to use DETECTORS_1KM as loop limit.

Revision 6.2  2009/05/30 23:32:12  kuyper
Corrected definition of sample_time to be at file scope.
Removed macros no longer referenced.
Corrected to validate geo_file.
Corrected frame->samp.
Changed to pass l1b_data to GEO_earth_location(), rather than
  scan_start_time[scan_number]
Corrected expected return value for GEO_validate_derived_products().

Revision 6.1  2009/05/28 22:51:25  kuyper
Changed to support both high-resolution and low-resolution data, connected
  by a call to GEO_aggregate.
Changed most global arrays to local arrays.
Moved definition and setting of sample_time global into this function.
Moved definition and setting of parameter globals into GEO_earth_location().
Corrected to cumulate gflags only after validation of derived products,
  resolving Bug 2199.
Changed interface of most subroutines.
Dropped qa_flag member of frame_state_struct.
Improved error messaging.
Removed sc_ev_frame_state argument from GEO_write_one_scan; long obsolete..

James Kuyper Jr.	James.R.Kuyper@nasa.gov

Revision 5.4  2006/01/27 17:19:06  kuyper
Corrected to set impulse_flag to GOOD_DATA only if GEO_interp_ECR() call
  succeeds.

Revision 5.3  2005/12/02 22:42:05  kuyper
Corrected to write scan data that does not depend upon earth view start time,
  even if earth view start time is a fill value.

Revision 5.2  2004/10/14 20:22:53  kuyper
Changed to pass geo_params->sol_elev_cor to GEO_interp_ECR().

Revision 5.1  2004/08/30 17:02:27  vlin
1. input parameter maneuver_list added.
2. Check for invalid EV start time.
3. call to GEO_write_one_scan updated.

Revision 4.8  2004/05/12 19:45:26  kuyper
Corrected to initialize impulse flag with GOOD_DATA, where appropriate,
rather than relying upon default initialization (which is not actually
guaranteed to occur).

Revision 4.7  2004/04/09 18:09:08  kuyper
Corrected to make MODIS_E_GEO error from GEO_interp_ECR() non-fatal.

Revision 4.6  2003/12/26 19:08:49  vlin
Two functions added in "Routines called" section.

Revision 4.5  2003/10/22 19:31:06  vlin
Beef up message for invalid input argument(s).

Revision 4.4  2003/08/26 19:31:53  kuyper
Corrected test for INVALID_SOLAR_ANGLES bit.

Revision 4.3  2003/08/22 22:08:49  kuyper
Changed to give a T_sc2ecr work area to GEO_interp_ECR().

Revision 4.2  2003/08/22 14:24:50  kuyper
Moved in definition of T_inst2ecr.
Corrected transfer of eulerAngles.

Revision 4.1  2003/08/13 20:40:15  kuyper
Added calls to GEO_get_sample_time() and GEO_interp_ECR().
Changed interface of GEO_derived_products to process an entire scan in one
  call, requiring splitting of the frame loop.
Moved the definitions of sample_to_sensor and sample_solar_angles into this
  module.
Added code to move eulerangles to sc_attitude.

5/23/95
Frederick S. Patt (patt@modis-xl.gsfc.nasa.gov)
Finished coding.

Requirements:
	PR03-F-3.1.1-1
	PR03-F-3.1.1-2
	PR03-F-3.1.2-1
	PR03-F-3.1.2-2
	PR03-F-3.1.3-1
	PR03-F-3.1.4-1
	PR03-F-3.2.1-1
	PR03-F-3.2.1-2
	PR03-F-3.2.2-1
	PR03-F-3.2.3-1
	PR03-F-3.2.3-2
	PR03-F-3.2.4-1
	PR03-F-3.3.1-1
	PR03-F-3.3.2-1
	PR03-F-3.3.3-1
	PR03-F-3.3.3-2
	PR03-F-3.4.1-1
	PR03-F-3.4.2-1
	PR03-F-3.4.3-1
	PR03-F-3.4.4-1
	PR03-F-3.5.1-1
	PR03-F-3.5.2-1
	PR03-F-3.5.3-1
	PR03-F-3.5.3-2
	PR03-F-3.5.3-3
	PR03-F-3.6.1-1
	PR03-F-3.6.2-1
	PR03-F-3.6.3-1
	PR03-F-3.6.3-2
	PR03-F_4.2-1
	PR03-F_4.2-2


!Team-unique Header:

	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

!END
**************************************************************************/

{
    frame_state_struct	sc_ev_frame_state[MAX_FRAMES] = {{0.0}};
    int			samp, frame, det;
    int                 missing_count, outofbounds_count; /* flag counters */
    uint32 mask;

    /* sample time relative to sector start time*/
    static double *frame_time;

    /* The largest of the following arrays are declared static, not so their
     * values can be retained from one call to the next, but to minimize
     * strain on stack.
     */
    /* pixel-to-satellite angles (azimuth, zenith), and range */
    typedef double frame_to_sensor_t[MAX_FRAMES][3];
    static frame_to_sensor_t *frame_to_sensor;

    /* solar azimuth and zenith angles for all the pixels in a scan */
    /* azimuth, zenith */
    typedef double frame_solar_angles_t[MAX_FRAMES][2];
    static frame_solar_angles_t *frame_solar_angles;
   
    /* instrument to ECR transformation */
    typedef double T_inst2ecr_t[3][3];
    static T_inst2ecr_t *T_inst2ecr;

    /* ECR ground location */
    typedef double ecr_sample_position_t[MAX_PADDED][3];
    static ecr_sample_position_t *ecr_sample_position;

    typedef double ecr_frame_position_t[MAX_FRAMES][3];
    static ecr_frame_position_t *ecr_frame_position;

    /* terrain corrected position in geodetic coordinates:
     * latitude, longitude, height
     */
    typedef double terrain_sample_position_t[MAX_PADDED][3];
    static terrain_sample_position_t *terrain_sample_position;

    typedef double terrain_frame_position_t[MAX_FRAMES][3];
    static terrain_frame_position_t *terrain_frame_position;

    /* Land/Water mask, and WaterPresent  */
    typedef uint8 sample_landsea_t[MAX_PADDED];
    static sample_landsea_t *sample_landsea;

    typedef uint8 frame_landsea_t[MAX_FRAMES];
    static frame_landsea_t *frame_landsea;

    typedef uint8 frame_waterpresent_t[MAX_FRAMES];
    static frame_waterpresent_t *frame_waterpresent;

    uint8 land_seamask_qaflag;

    /* geolocation flags */
    typedef uint8 sample_flags_t[MAX_PADDED];
    static sample_flags_t *sample_flags;

    typedef uint8 frame_flags_t[MAX_FRAMES];
    static frame_flags_t *frame_flags;

    /* spacecraft ECR position */
    typedef double ecr_sc_t[3];
    static ecr_sc_t *ecr_sc_sample_position;
    static ecr_sc_t *ecr_sc_frame_position;

    /* spacecraft ECR velocity */
    static ecr_sc_t *ecr_sc_velocity;
    
    static sc_state_struct *sc_state;

    char		asciiUTC[TIMECODEASIZE];
    PGSt_SMF_status	interp_ecr_status;

    typedef int8 hires_offsets_t[DETECTORS_QKM][SAMPLES_QKM];
    static hires_offsets_t *hires_offsets;

    static PGSt_double *offsets;
    
    typedef double T_sc2ecr_t[3][3];
    static T_sc2ecr_t *T_sc2ecr;

    typedef uint32 *sample_quality_t[2];
    static sample_quality_t *sample_quality;

    typedef uint32 *frame_quality_t[MAX_FRAMES];
    static frame_quality_t *frame_quality;
    
    static int firstRun = 1;
    if(firstRun) {
        firstRun = 0;
        
        sample_time = (double*) malloc(MAX_PADDED * sizeof(double));
        frame_time = (double*) malloc(MAX_FRAMES * sizeof(double));
        frame_to_sensor = (frame_to_sensor_t*) calloc(DETECTORS_1KM, sizeof(frame_to_sensor_t));
        frame_solar_angles = (frame_solar_angles_t*) calloc(DETECTORS_1KM, sizeof(frame_solar_angles_t));
        T_inst2ecr = (T_inst2ecr_t*) calloc(MAX_PADDED, sizeof(T_inst2ecr_t));
        ecr_sample_position = (ecr_sample_position_t*) calloc(DETECTORS_QKM, sizeof(ecr_sample_position_t));
        ecr_frame_position = (ecr_frame_position_t*) calloc(DETECTORS_1KM, sizeof(ecr_frame_position_t));
        terrain_sample_position = (terrain_sample_position_t*) calloc(DETECTORS_QKM, sizeof(terrain_sample_position_t));
        terrain_frame_position = (terrain_frame_position_t*) calloc(DETECTORS_1KM, sizeof(terrain_frame_position_t));
        sample_landsea = (sample_landsea_t*) calloc(DETECTORS_QKM, sizeof(sample_landsea_t));
        frame_landsea = (frame_landsea_t*) calloc(DETECTORS_1KM, sizeof(frame_landsea_t));
        frame_waterpresent = (frame_waterpresent_t*) calloc(DETECTORS_1KM, sizeof(frame_waterpresent_t));
        sample_flags = (sample_flags_t*) calloc(DETECTORS_QKM, sizeof(sample_flags_t));
        frame_flags = (frame_flags_t*) calloc(DETECTORS_1KM, sizeof(frame_flags_t));
        ecr_sc_sample_position = (ecr_sc_t*) calloc(MAX_PADDED, sizeof(ecr_sc_t));
        ecr_sc_frame_position = (ecr_sc_t*) calloc(MAX_FRAMES, sizeof(ecr_sc_t));
        ecr_sc_velocity = (ecr_sc_t*) calloc(MAX_PADDED, sizeof(ecr_sc_t));
    	sc_state = (sc_state_struct*) calloc(MAX_PADDED, sizeof(sc_state_struct));
        hires_offsets = (hires_offsets_t*) calloc(3, sizeof(hires_offsets_t));
        offsets = calloc(MAX_PADDED, sizeof(PGSt_double));
        T_sc2ecr = (T_sc2ecr_t*) calloc(MAX_PADDED, sizeof(T_sc2ecr_t));
        sample_quality = (sample_quality_t*) calloc(MAX_PADDED, sizeof(sample_quality_t));
        frame_quality = (frame_quality_t*) calloc(2, sizeof(frame_quality_t));
        
    }
    
    bzero(sample_flags, DETECTORS_QKM * sizeof(sample_flags_t));
    
    
    PGSt_double		rpy[3] =
	{THERMCORR_FVALUE, THERMCORR_FVALUE, THERMCORR_FVALUE};
    int			N_samp;
    int			padded_samples;
    char		filefunc[] = __FILE__ ", GEO_locate_one_scan";
    char		msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
    PGSt_SMF_status	retval = PGS_S_SUCCESS;

    if(geo_params == NULL || l1a_data == NULL || qa_metadata == NULL ||
	geo_file == NULL || scan_number < 0 || scan_number >= MAX_SCAN_NUMBER)
    {
        sprintf(msgbuf, "geo_params:%p, l1a_data:%p, \n"
               "qa_metadata:%p, geo_file:%p, scan_number:%d",
	       (void *)geo_params, (void *)l1a_data, (void *)qa_metadata,
	       (void*)geo_file, scan_number);
	modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);

	return MODIS_E_BAD_INPUT_ARG;
    }

    l1a_data->mirr_data[scan_number].impulse_flag = BAD_DATA;

    if ((fabs(scan_start_time[scan_number]-
                 l1a_data->fill_values.EV_start_time) < DBL_MIN) 
        || (l1a_data->frame_data[scan_number].EV_frames < 0)
	|| (l1a_data->frame_data[scan_number].EV_frames > MAX_FRAMES))
    {
       sprintf(msgbuf, "EV_start_time:%f EV_frames:%ld",
	   scan_start_time[scan_number],
	   (long)(l1a_data->frame_data[scan_number].EV_frames));
       modsmf(MODIS_E_BAD_INPUT_ARG, msgbuf, filefunc);
       l1a_data->frame_data[scan_number].EV_frames = 0;
       padded_samples = 0;
       interp_ecr_status = MODIS_E_BAD_INPUT_ARG;
    }
    else
    {
	N_samp = geo_params->geometry_params.N_samp
	    [geo_params->geometry_params.band_number];
	padded_samples =
	    (l1a_data->frame_data[scan_number].EV_frames + 1) * N_samp - 1;

	for (samp = 0; samp < padded_samples; samp++)
	{
	    sample_time[samp] =
		GEO_get_sample_time(&geo_params->geometry_params, samp);
	    if(sample_time[samp] >= GEO_DOUBLE_FILLVALUE)
	    {
		sprintf(msgbuf, "GEO_get_sample_time(%d)", samp);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		return MODIS_E_GEO;
	    }
	    offsets[samp] = (PGSt_double)sample_time[samp];
	}

	for(frame=0; frame<l1a_data->frame_data[scan_number].EV_frames; frame++)
	    frame_time[frame] = sample_time[N_samp * (frame + 1) - 1];

	if(padded_samples < MAX_PADDED)	/* Mark the end of the data. */
	    sample_time[padded_samples] = GEO_DOUBLE_FILLVALUE;

	interp_ecr_status = GEO_interp_ECR(
	    l1a_data->frame_data[scan_number].EV_start, padded_samples, offsets,
	    geo_params, asciiUTC, sc_state, T_sc2ecr, T_inst2ecr,
	    ecr_sc_sample_position, ecr_sc_velocity, rpy, sample_quality);
	if(interp_ecr_status !=  PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "GEO_interp_ECR(%f, %d)",
		l1a_data->frame_data[scan_number].EV_start, padded_samples);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    /* Most MODIS_E_GEO returns are for reasons that should be
	     * non-fatal.
	     */
	    if(interp_ecr_status != MODIS_E_GEO)
		retval = MODIS_E_GEO;
	}
	else
	    l1a_data->mirr_data[scan_number].impulse_flag = GOOD_DATA;
    }

    /* Loop on samples in scan */

    for (samp=0; samp<padded_samples; samp++)
    {
	if(interp_ecr_status != PGS_S_SUCCESS)
	{
	    for (det = 0; det < geo_params->num_detectors; det++)
		sample_flags[det][samp] |= INVALID_INPUT_DATA;
	}
	else if (l1a_data->mirr_data[scan_number].impulse_flag == GOOD_DATA)
	{
	     /* Perform Earth location */

	    if(GEO_earth_location(scan_number, samp,
		(GEO_param_struct*)geo_params, sample_time[samp], l1a_data,
		ecr_sc_sample_position, ecr_sc_velocity, T_inst2ecr,
		sample_flags, ecr_sample_position, terrain_sample_position)
		!= PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "GEO_earth_location(%d, %d)",
		    scan_number, frame);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);
		for (det = 0; det < geo_params->num_detectors; det++)
		    sample_flags[det][samp] |= INVALID_INPUT_DATA;

		retval = MODIS_E_GEO;
	    }

	    /* Validate Earth location */
	    if (GEO_validate_earth_location(
		l1a_data->mirr_data[scan_number].impulse_flag,
		terrain_sample_position, samp, geo_params->num_detectors,
		sample_flags) != SUCCESS)
	    {
		sprintf(msgbuf, "GEO_validate_earth_location(%d, %d)",
		    scan_number, samp);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		retval = MODIS_E_GEO;
	    }
	}

    }  /* End of sample loop */

    if(interp_ecr_status == PGS_S_SUCCESS)
    {
	if(l1a_data->mirr_data[scan_number].impulse_flag!=GOOD_DATA)
	    land_seamask_qaflag = 1;
	else
	{
	    if(GEO_landsea_mask(padded_samples, geo_params->num_detectors,
		terrain_sample_position, sample_flags, &land_seamask_qaflag,
		sample_landsea) != PGS_S_SUCCESS)
	    {
		sprintf(msgbuf, "GEO_landsea_mask(%d, %d) failed.\n",
		    padded_samples, geo_params->num_detectors);
		modsmf(MODIS_E_GEO, msgbuf, filefunc);

		retval = MODIS_E_GEO;
	    }
	}

	if(GEO_aggregate(l1a_data->frame_data[scan_number].EV_frames, N_samp,
	    geo_params->hires_scale, sample_flags, ecr_sample_position,
	    ecr_sc_sample_position, terrain_sample_position, sample_quality,
	    sample_landsea, ecr_frame_position, terrain_frame_position,
	    frame_flags, ecr_sc_frame_position, hires_offsets, frame_quality,
	    frame_landsea, frame_waterpresent) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "GEO_aggregate(%d, %ld, %d, %f)",
		geo_params->num_detectors,
		(long)l1a_data->frame_data[scan_number].EV_frames, N_samp,
		geo_params->hires_scale);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	if(GEO_derived_products(l1a_data->frame_data[scan_number].EV_frames,
	    DETECTORS_1KM, ecr_frame_position, ecr_sc_frame_position,
	    frame_flags, terrain_frame_position, asciiUTC, frame_time,
	    frame_solar_angles, frame_to_sensor) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "GEO_derived_products(%ld, %d, \"%s\")",
		(long)l1a_data->frame_data[scan_number].EV_frames,
		DETECTORS_1KM, asciiUTC);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	for(frame=0; frame<l1a_data->frame_data[scan_number].EV_frames; frame++)
	{
	    int ecr;

	    if((frame_flags[0][frame] & INVALID_SOLAR_ANGLES) == 0)
	    {
		/* Validate derived products */
		if (GEO_validate_derived_products(frame, DETECTORS_1KM,
		    geo_params->range_scale, frame_to_sensor, frame_flags)
		    != PGS_S_SUCCESS)
		{
		    sprintf(msgbuf, "GEO_validate_derived_products(%d, %d, %f)",
			frame, DETECTORS_1KM, geo_params->range_scale);
		    modsmf(MODIS_E_GEO, msgbuf, filefunc);

		    retval = MODIS_E_GEO;
		}
	    }

	    /* Update data quality measurments */
	    missing_count = 0;
	    outofbounds_count = 0;
	    for (det = 0; det < DETECTORS_1KM; det++)
	    {
		int bit;

		if ((frame_flags[det][frame] & INVALID_INPUT_DATA) != 0)
		    missing_count++;

		if ((frame_flags[det][frame] & INVALID_INPUT_DATA) == 0 &&
		    (frame_flags[det][frame] & NO_ELLIPSE_INTERSECT) != 0 )
		    outofbounds_count++;

		for (mask=1, bit=0; bit<8; bit++, mask*=2)
		{
		    if ((unsigned)frame_flags[det][frame] & mask)
			qa_metadata->cumulated_gflags[bit]++;
		}
	    }
	    qa_metadata->no_of_pixels += DETECTORS_1KM;
	    qa_metadata->missingdata += missing_count;
	    qa_metadata->outofboundsdata += outofbounds_count;

	    /* Subsampled frame data. */
	    samp = N_samp*(frame + 1) - 1;
    
	    for (ecr = 0; ecr < 3; ecr++)
	    {
		sc_ev_frame_state[frame].positionECR[ecr] =
		    ecr_sc_sample_position[samp][ecr];
		sc_ev_frame_state[frame].velocityECR[ecr] =
		    ecr_sc_velocity[samp][ecr];
		/* This array is ued only by GEO_cumulate_GRing(), and those
		 * are the only two fields used by that function.
		 */
	    }
	}
    }

    /* Write one scan of geolocation data to file */
    if (GEO_write_one_scan(frame_time, scan_number, l1a_data, geo_params,
	frame_to_sensor, frame_solar_angles, terrain_frame_position,
	frame_flags, hires_offsets, rpy, frame_quality, frame_landsea,
	frame_waterpresent, land_seamask_qaflag, geo_file) != PGS_S_SUCCESS)
    {
	sprintf(msgbuf, "GEO_write_one_scan(%d, \"%s\")",
	    scan_number, geo_file->filename);
	modsmf(MODIS_E_GEO, msgbuf, filefunc);

	retval = MODIS_E_GEO;
    }

    if(interp_ecr_status == PGS_S_SUCCESS)
    {
	if (GEO_cumulate_GRing(geo_params,
	    l1a_data->frame_data[scan_number].EV_frames, sc_ev_frame_state,
	    frame_flags, ecr_frame_position) != PGS_S_SUCCESS)
	{
	    sprintf(msgbuf, "GEO_cumulate_GRing(%ld)",
		(long)l1a_data->frame_data[scan_number].EV_frames);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}

	if (GEO_get_bounding_coords(terrain_frame_position, frame_flags,
	    DETECTORS_1KM,
	    l1a_data->frame_data[scan_number].EV_frames, bounding_coords)
	    != SUCCESS)
	{
	    sprintf(msgbuf, "GEO_get_bounding_coords(%d, %ld)", DETECTORS_1KM,
		(long)l1a_data->frame_data[scan_number].EV_frames);
	    modsmf(MODIS_E_GEO, msgbuf, filefunc);

	    retval = MODIS_E_GEO;
	}
    }
  
    return retval;
}

