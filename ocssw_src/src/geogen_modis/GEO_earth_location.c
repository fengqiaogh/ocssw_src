/* file: GEO_earth_location.c */

#include "PGS_MODIS_35251.h"
#include "smfio.h"
#include "GEO_earth.h"
#include "GEO_inst.h"

/* polynomial coefficients for mirror encoder-to-angle conversion */
double const * poly_coef;
int poly_degree;

/* mirr errors alpha, beta, gamma */
double alpha;
double beta;
double gammaa;

/* telescope to instrument frame transformation matrix */
double T_tel2inst[3][3];

/* mirror to instrument frame transformation matrix */
double T_mirr2inst[3][3];

/* mirror side 1 angle range (radians) */
double mirr_side1_range[2];

/* viewing vector in telescope coord */
/*  This is derived from instrument geometry parameters */
double u_tel[DETECTORS_QKM][3];

/* Earth sector start time for the scans in a granule */
double scan_start_time[MAX_SCAN_NUMBER];

/*
 * An array of pointers, one for each scan, to an array of mirror impulse
 * encoder (encoder counts) values for each impulse of that scan.
 */
double * mirr_impulse_enc[MAX_SCAN_NUMBER];

/*
 * An array of pointers, one for each scan, to an array fo mirror impulse
 * time (seconds) values for each impulse of that scan.
 */
double * mirr_impulse_time[MAX_SCAN_NUMBER];


PGSt_SMF_status GEO_earth_location(
	int			const scan_number,
	int			const sample_number,
        GEO_param_struct const	* geo_params,
        PGSt_double		sample_time,
        l1a_data_struct		* l1a_data,
        double			ecr_sc_sample_position[MAX_PADDED][3],
        double			ecr_sc_velocity[MAX_PADDED][3],
        double			T_inst2ecr[MAX_PADDED][3][3],
        unsigned char		sample_flags[DETECTORS_QKM][MAX_PADDED],
        double		ecr_sample_position[DETECTORS_QKM][MAX_PADDED][3],
        double		terrain_sample_position[DETECTORS_QKM][MAX_PADDED][3]
	)
/*
!C******************************************************************************
!Description:   
	Top-level routine in Earth Location group of the Level-1A geolocation
	software to calculate earth location for one sample in a scan.  This
	routine calls the routines which get the instrument view vectors,
	compute the ellipsoidal position and perform the terrain correction

!Input Parameters:
	scan_number             the scan number
	sample_number           the sample number for the ideal band
	geo_params              Geolocation parameters struct
	sample_time             Time when sample was collected, as an offset
				from the scan_start_time. Not used in this
				version of the code, but will be passed on to
				GEO_get_view_vec() in a future version.
	l1a_data		data read from the L1A file.
	ecr_sc_sample_position  interpolated sample ECR spacecraft postion
	ecr_sc_velocity         interpolated sample ECR spacecraft velocity

!Output Parameters:
	sample_flags            terrain sample position validity flags
	sample_time             pixel sample time offset from sector start time
	T_inst2ecr              telescope to instrument frame transformation
				matrices for each frame in a scan [INPUT]
	ecr_sample_position     ECR ground position for each sample.
	terrain_sample_position Geodetic ground positions for each sample.

Return Values:
	MODIS_E_GEO_BAD_INPUT_ARG       If geo_params is null.
	MODIS_E_GEO                     If any subroutine fails.
	PGS_S_SUCCESS                   Otherwise.

Externally Defined:
	From "GEO_global_arrays.h":
	Note: all of these global variables should be defined in this module;
	they are "externally defined" only in the sense that they have external
	linkage.

	poly_coef               polynomial coefficients for mirror
				encoder-to-angle conversion [OUTPUT]
	poly_degree             degree of the mirror polynomial [OUTPUT]

	alpha                   mirror wedge angle (radians) alpha [OUTPUT]
	beta                    mirror wedge angle (radians) beta [OUTPUT]
	gammaa                  mirror wedge angle (radians) gammaa [OUTPUT]
	T_tel2inst              telescope to instrument frame transformation
				matrix [OUTPUT]
	T_mirr2inst             mirror to instrument frame transformation
				matrix [OUTPUT]
	mirr_side1_range        mirror side 1 angle range (radians) [OUTPUT]
	scan_start_time         Earth View scan start times [OUTPUT]
	mirr_impulse_enc        an array of pointers, one for each scan, to an
				array of mirror impulse encoder (encoder
				counts) values for each impulse of that
				scan [OUTPUT]
	mirr_impulse_time       an array of pointers, one for each scan, to an
				array of mirror impulse time (seconds)
				values for each impulse of that scan [OUTPUT]

	CHAN_A				"GEO_parameters.h"
	CHAN_B				"GEO_parameters.h"
	DETECTORS_QKM                   "GEO_geo.h"
	MODIS_E_GEO                     "PGS_MODIS_35251.h"
	MODIS_E_GEO_BAD_INPUT_ARG       "PGS_MODIS_35251.h"
	MAX_PADDED                      "GEO_geo.h"
	PGS_S_SUCCESS                   "PGS_SMF.h"
	SUCCESS                         "GEO_basic.h"

Called by:
	GEO_locate_one_scan             "GEO_output.h"

Routines Called:
	GEO_get_view_vec                "GEO_inst.h"
	GEO_ellip_position              "GEO_earth.h"
	GEO_terrain_correct             "GEO_earth.h"
	modsmf                          "smfio.h"

!Revision History:
   $Log: GEO_earth_location.c,v $
   Revision 6.4  2013/05/23 18:30:55  jkuyper
   Change l1a_data argument of GEO_earth_location() to be a pointer to
     non-const, to allow initialization of mirr_impulse_enc.

   Revision 6.3  2011/02/14 21:24:44  kuyper
   Corrected const-qualification of several pointers.

   Revision 6.2  2009/05/30 23:50:38  kuyper
   Added global variables that were supposed to be moved to this module
     from GEO_locate_one_scan.c.
   Changed to take l1a_data as a parameter, and to calculate
     scan_start_time, mirr_impulse_ang, and mirr_impulse_time from it.
   Changed to set poly_coef.
   Made error messages more informative.

   Revision 6.2  2009/05/30 18:29:32  kuyper
   Added global variables that were supposed to be moved to this module
     from GEO_locate_one_scan.c.
   Changed to take l1a_data as a parameter, and to calculate
     scan_start_time, mirr_impulse_ang, and mirr_impulse_time from it.
   Made error messages more informative.

   Revision 6.1  2009/05/29 13:57:24  xgeng
   Move definition, initialization of some global variables into this module
     from GEO_locate_one_scan().
   Changed sample_view_vec and ellip_sample_position from global arrays to local.
   Changed to pass what used to be global arrays as arguments to and from
     GEO_ellip_position(), and GEO_terrain_correct().
   Changed to expect status code values from GEO_ellip_position() and
     GEO_terrain_correct(), and to return a status code.
   Cleaned up error detection and error messaging.
   Changed MAX_SCAN_SAMPLES to MAX_PADDED.
   The above changes match with PDL revision 6.1 and 6.2.
   Please note that 'sample_time' and 'scan_start_time' are not currently used.
     They are added at this time, so that in some future revision,
     we can pass them through GEO_get_view_vec() and GEO_get_inst_mirr_normal()
     all the way down to GOE_interp_mirr_enc().

   Xu Geng (xu.geng@saic.com)

   Revision 2.1  1997/10/21 18:16:22  kuyper
   Returned from ClearCase

 * Revision 1.6.1.1  1997/07/21  22:20:17  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.6  1997/07/21  16:24:34  kuyper
 * Baselined Version 1
 *
 *Parallel Development:
 * Revision 1.7  1997/04/21  22:28:51  fhliang
 * removed unused argument 'scan_number' (LL.112-115).
 *
 * Revision 1.6  1997/03/26  18:02:57  fhliang
 * Initial revision of GEO_earth_location.c.
 *
				Revision 1.5  1997/02/13 19:28:56  kuyper
				Merged seed files.

				Revision 1.4  1996/07/24 20:59:10  kuyper
				Standardized order of #include files.
				Declared arguments const.

				Revision 1.3  1996/07/23 22:59:14  kuyper
				Inserted required '!' in comments.
				Changed constants to double to avoid conversions.
				Removed ret_val.

				Revision 1.2  1996/07/18 15:58:00  kuyper
				Included self-checking header file, and other needed header files.


                4/25/95
                Ruiming Chen (rchen@ltpmail.gsfc.nasa.gov)
                Finished coding.
   
                7/3/95
                Tracey W. Holmes (holmes@modis-xl.gsfc.nasa.gov)
                Added SDP error messages.

                10/10/95
                Tracey W. Holmes
                Added debug option.

!Team-unique Header:
                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END*************************************************************************
*/
{
  int i, j, det; 
  double u_inst[DETECTORS_QKM][3] = {0.0};

  typedef double padded3_t[MAX_PADDED][3];
  static padded3_t *sample_view_vec;
  static padded3_t *ellip_sample_position;
  static int firstRun = 1;
  if(firstRun) {
    firstRun = 0;
    sample_view_vec = (padded3_t*) calloc(DETECTORS_QKM, sizeof(padded3_t));
    ellip_sample_position = (padded3_t*) calloc(DETECTORS_QKM, sizeof(padded3_t));
  }

  char msgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  char filefunc[] = __FILE__ ", GEO_earth_location";

  if (geo_params == NULL) {
    modsmf(MODIS_E_BAD_INPUT_ARG, "geo_params is null", filefunc);
    return MODIS_E_BAD_INPUT_ARG;
  } 

  if (poly_degree < geo_params->poly_degree) {
    alpha = geo_params->mirror_model.alpha; 
    beta = geo_params->mirror_model.beta;
    gammaa = geo_params->mirror_model.gammaa;
    mirr_side1_range[0]  = geo_params->mirror_model.mirr_side1_range[0];
    mirr_side1_range[1]  = geo_params->mirror_model.mirr_side1_range[1];
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
      {
        T_mirr2inst[i][j] = geo_params->coord_trans.T_mirr2inst[i][j];
        T_tel2inst[i][j]  = geo_params->coord_trans.T_tel2inst[i][j];
      }
    for (det = 0; det < geo_params->num_detectors; det++)
    { 
      for (i = 0; i < 3; i++)
        u_tel[det][i] = geo_params->u_tel[det][i]; 
    }
    poly_degree = geo_params->poly_degree;
  } 
  scan_start_time[scan_number] = l1a_data->frame_data[scan_number].EV_start;
  mirr_impulse_enc[scan_number] =
      l1a_data->mirr_data[scan_number].mirr_impulse_enc;
  mirr_impulse_time[scan_number] =
      l1a_data->mirr_data[scan_number].mirr_impulse_time;
  if(l1a_data->mirr_data[scan_number].mirr_chan == CHAN_A)
      poly_coef = geo_params->poly_coef[CHAN_A];
  else
      poly_coef = geo_params->poly_coef[CHAN_B];

  /* get view vector */
  if (GEO_get_view_vec(scan_number, sample_number, geo_params->num_detectors,
      u_inst) != SUCCESS)
  {
    /* call SDP function to report error */
    sprintf(msgbuf, "GEO_get_view_vec(%d, %d, %d)",
	scan_number, sample_number, geo_params->num_detectors);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);
    return MODIS_E_GEO;
  }

  /* get ellipsoid position */
  if (GEO_ellip_position(scan_number, sample_number, geo_params->num_detectors,
      u_inst, ecr_sc_sample_position, ecr_sc_velocity, T_inst2ecr,
      ecr_sample_position, ellip_sample_position, sample_flags, sample_view_vec)
      != PGS_S_SUCCESS)
  {
    /* call SDP function to report error */
    sprintf(msgbuf, "GEO_ellip_position(%d, %d, %d)",
	scan_number, sample_number, geo_params->num_detectors);
    modsmf(MODIS_E_GEO, msgbuf, filefunc);
    return MODIS_E_GEO;
  }

  /* get terrain corrected position */
  if (GEO_terrain_correct(sample_number, geo_params->num_detectors,
      sample_view_vec, ecr_sample_position, ellip_sample_position, sample_flags,
      terrain_sample_position ) != PGS_S_SUCCESS)
  {
    /* call SDP function to report error */
    sprintf(msgbuf, "GEO_terrain_correct(%d, %d)",
	sample_number, geo_params->num_detectors);
    modsmf(MODIS_E_GEO, msgbuf, filefunc); 
    return MODIS_E_GEO;
  }

  return PGS_S_SUCCESS;
}

