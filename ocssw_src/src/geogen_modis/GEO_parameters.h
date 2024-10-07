#ifndef GEO_PARAMETERS_H
#define GEO_PARAMETERS_H
#include "PGS_PC.h"
#include "GEO_geo.h"
#include "GEO_product.h"
#include "hdfi.h"

/*
!C-INC**************************************************************************
!Description:	The declaration of the global parameter variables for the
		Level-1A geolocation software

!Input Parameters: N/A

!Output Parameters: N/A

!Revision History:
 * $Log: GEO_parameters.h,v $
 * Revision 6.9  2011/02/14 21:28:31  kuyper
 * Corrected to const-qualify *poly_coef.
 *
 * Revision 6.8  2010/07/01 19:38:43  kuyper
 * Corrected definition of INPUTPOINTERS.
 *
 * Revision 6.7  2010/04/23 14:55:09  kuyper
 * Converted new numbers controlling setting of quality flags into parameters
 *   retrieved from the parameter file.
 *
 * Revision 6.6  2010/04/08 18:01:21  kuyper
 * Dropped orbit and attitude globals that are no longer needed.
 * Helped resolve Bug 17 by adding field to fill_values_struct.
 *
 * Revision 6.5  2009/05/30 18:09:48  kuyper
 * Corrected size of u_tel.
 *
 * Revision 6.4  2009/05/22 16:33:43  xgeng
 * reinstated member qa_flag for frame_state_struct.
 *
 * Revision 6.3  2009/05/19 14:44:49  kuyper
 * Reinstated T_mirr2inst, u_tel globals; cannot be dropped until after
 * GEO_get_inst_mirr_normal() and GEO_get_view_vec() have been opened.
 *
 * Revision 6.2  2009/04/15 21:18:29  kuyper
 * Removed member from frame_state_struct.
 *
 * Revision 6.1  2009/03/12 17:59:47  kuyper
 * Changed MAX_DETECTORS to DETECTORS_QKM.
 * Added hires_scale.
 * Dropped globals no longer in use.
 *
 * Revision 5.4  2005/09/27 15:12:00  vlin
 * Four fieldnames added to orbit_validation_params_struct
 * vlin@saicmodis.com
 *
 * Revision 5.3  2004/08/30 21:15:23  kuyper
 * Entered first non-dummy value for NUM_SOL_ELEV.
 *
 * Revision 5.2  2004/08/23 18:56:31  vlin
 * fill_values_struct added, one member "fill_values" added to l1a_data_struct.
 *
 * Revision 5.1  2004/08/02 19:36:29  vlin
 * NCSA acknowledgement removed from prologue,
 * solar elevation correction matrix added.
 *
 * Revision 4.5  2003/08/22 22:18:52  kuyper
 * Changed T_inst2SD to T_sc2SD.
 *
 * Revision 4.4  2003/08/05 21:12:19  kuyper
 * Corrected location of processing metadata fields.
 *
 * Revision 4.3  2003/01/29 20:17:30  vlin
 * Corrected declaration for encoder_gap
 *
 * Revision 4.2  2002/11/21 22:54:59  kuyper
 * Added in packet_interval.
 *
 * Revision 4.1  2002/11/15 21:24:29  kuyper
 * Added items needed to support temperture retrieval, expanded input pointers,
 *   and reprocessing and processingenvironment metadata.

!Team-unique Header:

                This software is developed by the MODIS Science Data Support
                Team for the National Aeronautics and Space Administration,
                Goddard Space Flight Center, under contract NAS5-32373.

!END****************************************************************************
*/

/* Identifying tag required on first line of geolocation parameter file. */
#define PARAM_TAG "MOD03_L1A_GEOLOCATION"
#define SCAN_TYPE_LEN 10

/**********************************************************
Variables from the parameter file
**********************************************************/

/* polynomial coefficients for mirror encoder-to-angle conversion */
extern double const * poly_coef;
extern int poly_degree;

/* mirr errors alpha, beta, gamma */
extern double alpha;
extern double beta;
extern double gammaa;

/* mirror side 1 angle range (radians) */
extern double mirr_side1_range[2];

/* telescope to instrument frame transformation matrix */
extern double T_tel2inst[3][3];

/* mirror to instrument frame transformation matrix */
extern double T_mirr2inst[3][3];

/* viewing vector in telescope coord */
/*  This is derived from instrument geometry parameters */
extern double u_tel[DETECTORS_QKM][3];

/* scale factors for converting ancillary input data */
typedef struct{
	double S_position;      /* orbit position (1/8 meters/count) */
	double S_velocity;      /* orbit velocity (1/4096 meters/sec/count) */
	double S_attitude;      /* attitude angles (1 arcsec/count) */
	double S_angvel;        /* angular velocity (0.5 arcsec/sec/count) */
} ancil_scale_struct;

/* location of data words in l1a data segments */
extern int sector_word[2];      /* Earth sector encoder word [side 1/side 2] */
extern int vernier_word[2];     /* Vernier count word [side 1/side 2] */
typedef struct{
	int time_words[2];	/* ancillary time start word [current/prior] */
	int orbit_words[2];	/* orbit data start word [current/prior] */
	int attit_words[2];	/* attitude data start word [current/prior] */
} ancil_word_struct;

typedef struct{
	ancil_word_struct ancil_words;
	ancil_scale_struct ancil_scale_factors;
} ancil_data_param_struct;

typedef struct {
           uint8     num_impulse;   /* number of good impulse data */
           int       impulse_flag;  /* encoder data flag    */
           uint16    mirr_side;                       /* mirror side*/
           double    mirr_impulse_enc[MAX_IMPULSE_NUMBER];  /* mirror impulse encoder */
           double    mirr_impulse_time[MAX_IMPULSE_NUMBER]; /* times (TAI seconds)  */
           uint16    mirr_chan;
        } mirr_data_struct;

typedef struct {
           int32    L1A_scan_id;
           int32    L1A_scan_quality[4]; /* Scan quality array              */
           float64  EV_start;            /* Earth View scan start time (TAI seconds) */
           int32    EV_frames;           /* Earth View frames in each scan  */
           float64  SD_start;            /* Solar Diffusor scan start times (TAI seconds) */
           int32    SD_frames;           /* Solar Diffusor frames in each scan */
           float64  SV_start;            /* Space View scan start times (TAI seconds)  */
           int32    SV_frames;           /* Space View frames in each */
           int8     SCI_ABNORM;          /* spacecraft in unusual state
                                            (i.e. maneuver) flag
                                            0 = maneuver
                                            1 = normal                     */
           int8     SCI_STATE;           /* instrument test mode flag
                                            0 = MODIS test mode
                                            1 = normal                     */
           char8    Scan_type[SCAN_TYPE_LEN]; /*'Day', 'Night', or 'Other'  */
        } frame_data_struct;

typedef struct{
           float64  time;                /* frame (mid-)sample time, TAI seconds   */
           float64  positionECR[3];      /* ECR postion, meters                    */
           float64  velocityECR[3];      /* ECR velocity, meters/sec               */
           float64  eulerAngles[3];      /* attitude euler angles (roll, pitch, yaw) */
                                        /* wrt the orbital reference frame (?), radians. */
           float64  T_inst2ecr[3][3];    /* instrument to ecr coordinate transformation matrix. */
           int      qa_flag;          /*     kinematic state quality (FAIL, SUCCESS)  */
        } frame_state_struct;

typedef struct {
	   frame_state_struct sc;
           float64  sun_unit_vector[3];
           float64 SD_sun_zenith;
           float64 SD_sun_azimuth;
           float64 moon_unit_vector[3];
        } celestial_bodies_struct;
/* extern celestial_bodies_struct cb_vectors; */

typedef struct {
        int32   max_earth_frames;
        int32   max_sd_frames;
        int32   max_sv_frames;
        int32   incomplete_scans;
        int32   missing_packets;
        int32   bad_CRC_packets;
        int32   discarded_packets;
        int32   Extract_Pixel_Offset;
        int32   Extract_Pixel_Count;
        int32   Extract_Line_Offset;
        int32   Extract_Line_Count;
       } l1a_metadata_struct;

typedef struct {
        char    localversionid[17];
        char    pgeversion[11];
        char    processversion[11];
        char    productionhistory[32];
	char	processingenvironment[PGSd_PC_VALUE_LENGTH_MAX];
	char	reprocessingplanned[PGSd_PC_VALUE_LENGTH_MAX];
	char	reprocessingactual[PGSd_PC_VALUE_LENGTH_MAX];
      } version_metadata_struct;

typedef struct {
        char    rangebeginningdate[11];
        char    rangebeginningtime[16];
        char    rangeendingdate[11];
        char    rangeendingtime[16];
        char    operationmode[8];
        char    localinputgranuleid[PGSd_PC_FILE_PATH_MAX];
        char    granulenumber[8];
        char    platformshortname[20];
        version_metadata_struct version_metadata;
      } ECS_metadata_struct;

typedef struct {
        PGSt_integer    orbitnumber;
        PGSt_double     equatorcrossinglongitude;
        char            equatorcrossingdate[11];
        char            equatorcrossingtime[16];
      } EPH_metadata_struct;

typedef struct {
           int16     scan_number;
           float64   EV_start_time;
           float64   SD_start_time;
           float64   SV_start_time;
           float32   SD_sun_zenith;
           float32   SD_sun_azimuth;
           float32   moon_vector;
           float32   sun_ref;
           uint16    mirr_side;
	   uint16    raw_mir_enc;
           float64   impulse_enc;
           float64   impulse_time;
           int32     L1_scan_quality;
           int8      geo_scan_quality;
           float64   EV_center_time;
           float64   orb_pos;
           float64   orb_vel;
           float64   T_inst2ECR;
           float64   attitude_angels;
           int32     EV_frames;
           int32     SD_frames;
           int32     SV_frames;
           uint8     num_impulse;
           char8     Scan_type[SCAN_TYPE_LEN];
        } fill_values_struct;

#define NUM_TEMPS	6
#define NUM_SOL_ELEV    36

typedef struct {
            int			num_scans;
            frame_data_struct	frame_data[MAX_SCAN_NUMBER];
            mirr_data_struct	mirr_data[MAX_SCAN_NUMBER];
            l1a_metadata_struct	granule_metadata;
            ECS_metadata_struct	ECS_metadata;
	    float32		temperatures[NUM_TEMPS];
            fill_values_struct  fill_values;
        } l1a_data_struct;
/* extern l1a_data_struct l1a_data; */

typedef struct {
            int     band_number; /* band number to geolocate (0 is ideal band) */
            double  t_reset;     /* time to reset sample at beginning of frame */
            double  t_frame;     /* sample frame time */
            /* offset in IFOV units to get time of first sample for a band */
            double  F_offset[MAX_BAND_NUMBER+1];
            /*  number of samples per frame for each band */
            uint16  N_samp[MAX_BAND_NUMBER+1];
            float64 focal_length[MAX_BAND_NUMBER+1]; /* focal length */
            /* detector position for each band */
            float64 det_space[MAX_BAND_NUMBER+1];
            float64 det_position[MAX_BAND_NUMBER+1][2];
            float64 band_position[MAX_BAND_NUMBER+3]; /* band center position */
        } focal_plane_geometry_struct;
/* extern focal_plane_geometry_struct geometry_params; */

typedef struct {
    /* Number of sample periods from sector start to start of data collection  */
            double  N_reset;
            double  mirr_abs_limit[2]; /* mirror encoder time absolute limits */
            double  mirr_del_limit; /* mirror encoder time delta limits */
            int     sample_impulse; /* encoder sample period in impulses */
            int     sector_word[2]; /* Earth sector encoder word [side 1/side 2] */
            int     vernier_word[2]; /* Vernier count word [side 1/side 2] */
            double  t_encoder; /* encoder time scale factor */
            double  t_vernier; /* vernier offset count time */
	    double  encoder_gap[2]; /* Range affected by encoder_adjustment. */
	    double  encoder_adjustment;	/* Size of encoder correction.	*/
	    double  packet_interval;	/* Interval between ancillary packets.*/
        }  mirror_preparation_struct;
/* extern  mirror_preparation_struct       mirror_prep_params; */
/* offsets into mirr_abs_limit */
enum {lower_limit, upper_limit};

typedef struct {
            double  mirr_side1_range[2]; /* mirror side 1 angle range (radians) */
            float64 alpha;
            float64 beta;
            float64 gammaa;
        }  mirror_model_struct;
/* extern  mirror_model_struct  mirror_model; */

typedef struct {
            double position_abs_limit[2]; /* orbit position absolute limits */
            double position_mag_limit[2]; /* orbit position magnitude limits */
            double velocity_abs_limit[2]; /* orbit velocity absolute limits */
            double velocity_mag_limit[2]; /* orbit velocity magnitude limits */
            double ang_mom_limit[2];      /* angular momentum magnitude limits */
            double ang_mom_z_limit[2];    /* angular momentum Z component */
            /* orbit position/velocity consistency limit */
            double orbit_consistency;
            double descend_time_0[2];
	              /* the descend time corresponding to orbit number 0 */
            double orbit_tolerance[2];
            double period[2];
            int    transition_orbit;
	    double eph_max_short_gap;	/* Maximum length of a short gap */
        } orbit_validation_params_struct;
/* extern orbit_validation_params_struct  orbit_valid_params; */

typedef struct {
            double  angvel_abs_limit[2]; /* angular velocity absolute limits */
            double  angvel_del_limit; /* angular velocity delta limits */
            double  attitude_abs_limit[2]; /* attitude angle absolute limits */
            double  attitude_del_limit; /* attitude angle delta limits */
            double  attit_consistency;     /* angle/angular velocity consistency limit */
	    double att_max_short_gap;	/* Maximum length of a short gap. */
	    short  att_valid_range[2];	/* Limits for valid entrained attitude */
        } attitude_valid_struct;
/* extern attitude_valid_struct       attit_valid_params; */

typedef struct {
            float64 T_inst2sc[3][3];	/* instrument to spacecraft */
            float64 T_mirr2inst[3][3];	/* mirror to instrument frame */
            float64 T_tel2inst[3][3];	/* telescope to instrument frame */
            float64 T_sc2SD[3][3];	/* spacecraft to solar diffuser */
	    int	    rpy_count;	     	/* count of rows in rpy_inst2sc */
	    double  (*rpy_inst2sc)[3];	/* Roll/Pitch/Yaw corrections for
					 * T_inst2sc */
	    double  *rpy_times;		/* The corresponding TAI times */
        } internal_coord_trans_struct;
/* extern internal_coord_trans_struct coord_trans; */

/*
 * ELEC_SIDES has two possible values, which are defined below.
 */
#define CHAN_A 0
#define CHAN_B 1

#define NUM_TEMP_EQNS	7
#define TEMP_ORDER	6

typedef struct {
	    char			revision[32];
            focal_plane_geometry_struct geometry_params;
            mirror_preparation_struct   mirror_prep_params;
            mirror_model_struct         mirror_model;
            ancil_data_param_struct     ancil_params;
	    double			max_non_gap;
	    /* Minimum interval between packets that is considered a gap. */
            orbit_validation_params_struct    orbit_valid_params;
            attitude_valid_struct       attit_valid_params;
            internal_coord_trans_struct coord_trans;
            /* number of detectors for band_number */
            int                         num_detectors;
            /* telescope view vector in meters. */
            double                      u_tel[DETECTORS_QKM][3];
            /* maximum orbit/attitude extrapolation limit (seconds) */
            double                      max_extrap;
            /* polynomial coefficients for mirror encoder-to-angle conversion */
            double                      poly_coef[ELEC_SIDES][MAX_POLY_DEGREE+1];
            int                         poly_degree;
            /* scale factors for scaling geolocation output values:  */
            double                      angle_scale; /* counts/radian   */
            double                      hires_scale; /* kilometers IFOV */
            double                      range_scale; /* counts/meter    */
            float32                     RMS_error; /* estimated geo RMS error (meters). */
	    char			spacecraft_ID[8];

	    float32			temp_range[NUM_TEMP_EQNS][2];
	    /* The valid output range of the temperature conversion functions.*/
	    float32			temp_coeff[NUM_TEMP_EQNS][TEMP_ORDER];
            double                      sol_elev_cor[NUM_SOL_ELEV][3];
	    /* The polynomial coefficients of the temperature conversion
	     * functions.
	     */
        } GEO_param_struct;

typedef struct {
            char  exclusionflag[2]; /*  Flag indicating whether points
                                        are on an inner (exclusion) G-ring.*/
            int   sequenceno[4]; /*     Sequence numbers corresponding to
                                        perimeter latitudes and longitudes. */
            double  latitude[4]; /*     Latitudes of a series of points
                                        representing the perimeter of the
                                        granule spatial coverage. */
            double  longitude[4]; /*    Longitudes of a series of points
                                        representing the perimeter of the
                                        granule spatial coverage. */
} GEO_GRing_struct;

#define INPUTPOINTERS (2 + 2*MAX_EA_FILES)
typedef struct {
	char	inputpointer[INPUTPOINTERS][PGSd_PC_UREF_LENGTH_MAX];
} pointer_metadata_struct;

typedef struct {
	int	no_of_pixels;	/* total number of pixels in granule */
	int	missingdata;	/* no. of pixels for which there are
				  insufficient data to perform geolocation */
	int	outofboundsdata;/* no. of pixels that have un-
				  geolocatable views */
        int     retval;         /* FAIL/SUCCESS, depending upon whether a
                                  fatal error was detected */
        char    rms_error[9];   /*      string value of RMS geolocation error*/
        uint32  cumulated_gflags[8]; /* Cumulated count of pixels for which the
                                        corresponding bit of gflags was set.
                                        LSB corresponds to first element.*/

} qa_metadata_struct;

#endif

