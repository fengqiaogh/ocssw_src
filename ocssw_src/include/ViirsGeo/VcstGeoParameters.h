/**************************************************************************
 *
 * NAME: VcstGeoParameters
 *
 * DESCRIPTION: Contains variables and structures used by the ProSdrViirs
 *  geolocation routines.
 *
 *
 **************************************************************************/
/* file: VcstGeoParameters.h */

#ifndef GEO_PARAMETERS_H
#define GEO_PARAMETERS_H

#include <VcstCmnConsts.h>
#include <VcstCalLutStructures.h>

#define DM_BADADDRESS ((void *) -1L)

const int EV_FRAMES_375M = 6400;
const int EV_FRAMES_750M_UNAGG = 6304;
const int EV_FRAMES_750M = 3200;
const int EV_FRAMES_DNB = 4064;

const int M_DETECTORS = 16;
const int I_DETECTORS = 32;
//const int NUM_DETECTORS_DNB = 16;
const int DNB_DETECTORS = 16;

const int I_VIIRS_COLS = 6400;
const int I_VIIRS_ROWS = Number_of_Scans*I_DETECTORS;
const int M_VIIRS_COLS = 3200;
const int M_VIIRS_ROWS = Number_of_Scans*M_DETECTORS;
const int M_A_VIIRS_COLS = 5024;
const int DNB_VIIRS_COLS = 4064;
const int DNB_VIIRS_ROWS = Number_of_Scans*M_DETECTORS;

// VIIRS Moderate dual gain number of Columns
const int M_D_VIIRS_COLS = 6304;

//const int VEC_SIZE = 3;
const int VEC_RPY_SIZE = 3;

const int NUM_BITS_PER_BYTE = 8;
//const int NUM_DETECTORS = 432;
//const int NUM_BANDS = 22;

//const int NUM_THERMISTORS = 26;
const int MAX_POLY_DEGREE_PLUS_ONE = 5;
const int MAX_BAND_NUMBER_PLUS_ONE = 37;
//const int NUM_ELECTRONICS_SIDE = 2;
const int MAX_THERMISTOR_ID_LEN = 40;
const int NUM_AGG_ZONES = 5;
//const int NUM_MOON_OFFSET_LIMITS = 4;

const int MAX_BAND_NUMBER = 22;
const int MAX_SCAN_NUMBER = Number_of_Scans;
const int MAX_DETECTORS = 32;
const int MAX_POLY_DEGREE = 4;
const int CHEBY_ORDER = 4;

/* telescope encoder impulses/Earth view */
const int MAX_IMPULSE_NUMBER_TEL = 1290;
const int ENCODER_LENGTH_TEL = 1290;

/* mirror encoder impulses/Earth view    */
const int MAX_IMPULSE_NUMBER_HAM = 1290;
const int ENCODER_LENGTH_HAM = 1290;

const int SECTOR_LENGTH = 32; /* no. of view sector def, words/scan    */
const int ANCIL_LENGTH = 64; /* no. of S/C ancillary data words/scan  */

/* used for flags to indicate data quality */
const int GOOD_DATA = 0;
const int BAD_DATA = -1;
const int SERVO_CNTRL_BAD_DATA = -2;
const int DEGRADED_DATA = 2;

const int GEO_FAIL = -1; /* GEO_FAIL represents a default index if no index
                                of the next unflagged value is found */

const int mirr_lower_limit = 0;
const int mirr_upper_limit = 1;

const int tel_lower_limit = 0;
const int tel_upper_limit = 1;

// const double MAX_ANGLE = PIO2;

// Desciption of interpolation rectangles
const int MOD_VIIRSI_NRCTNGL = 400;
const int MOD_VIIRSI_TWIDTH = 8;
const int MOD_MIDDLE_ROW = 7;

const int DNB_VIIRSI_NRCTNGL = 508;
const int DNB_VIIRSI_TWIDTH = 8;
const int DNB_MIDDLE_ROW = 7;

const int IMG_VIIRSI_NRCTNGL = 400;
const int IMG_VIIRSI_TWIDTH = 16;
const int IMG_MIDDLE_ROW = 15;

const int VIIRSI_MAX_NRCTNGL = DNB_VIIRSI_NRCTNGL;

/* used for RPY time-based inst2sc adjustments */
const int MAX_RPY_COUNT = 1000;
const int RPY_TIMECODE_SIZE = 20;

/************************************************************************
 * Constants for the VIIRS resolution aggregation zones at the SDR level
 **********************************************************************/
const int NumModImgAggZones = 5; // 5 zones because zone 3,4 are the same

const int ViirsModPixPerZone[NumModImgAggZones] ={640, 368, 1184, 368, 640};

const int ViirsImgPixPerZone[NumModImgAggZones] ={1280, 736, 2368, 736, 1280};

// For DNB the array begins at the edge of scan and ends at nadir.  The
// last number is the 2 nadir zones combined
const int NumDnbAggZones = 63; // 63 zones because zones 32,33 are the same
const int ViirsDnbPixPerZone[NumDnbAggZones] ={80, 16, 64, 64, 64, 32, 24, 72, 40, 56, 40, 48, 32, 48, 32, 72,
    72, 72, 80, 56, 80, 64, 64, 64, 64, 64, 72, 80, 72, 88, 72, 368,
    72, 88, 72, 80, 72, 64, 64, 64, 64, 64, 80, 56, 80, 72, 72,
    72, 32, 48, 32, 48, 40, 56, 40, 72, 24, 32, 64, 64, 64, 16, 80};


// enumerated indexes for the temperature correction values

enum AttitudeType {
    ROLL_IDX = 0, // roll index
    PITCH_IDX, // pitch index
    YAW_IDX // yaw index
};

enum GeoFillActionType {
    FILL_POINT = 0, // fill a single point
    FILL_DETECTOR_ARRAY, // fill all the detectors for a frame
    FILL_SCAN_ARRAY // fill values for an entire scan level array
};

const int MAX_SCAN_SAMPLE = 6400;
const int NUM_DNB_ZONES = 64;

const int NUM_TEMPERATURE_INPUTS = 9 + NUM_THERMISTORS;
const int VECTOR_SIZE = 3; /* thermal correction vector */

/* flag conditions */
const unsigned char INVALID_INPUT_DATA = 128;
const unsigned char NO_ELLIPSE_INTERSECT = 64;
const unsigned char BAD_TERRAIN = 32;
const unsigned char INVALID_SENSOR_ANGLES = 8;


const double MIN_TERRAIN_HEIGHT = -450.0; /* units in meters */
const double MAX_TERRAIN_HEIGHT = 9600.0; /* units in meters */

const int ORDER_CHEBY = 4;

enum LocationType {
    LAT_IDX = 0,
    LON_IDX,
    HT_IDX
};

enum DirectionType {
    SAT_AZM_IDX = 0,
    SAT_ZEN_IDX,
    SOL_AZM_IDX,
    SOL_ZEN_IDX,
    LUN_AZM_IDX,
    LUN_ZEN_IDX
};


const unsigned short MAX_UINT16_VAL = 65535;

enum DirectionDistanceType {
    azimuth_index = 0,
    z_angle_index,
    range_index
};

/* Constants for temperature correction LUT */
const int MAX_TEMPERATURE_BRKTS = 11; /*this is a guess */
const int MAX_TEMPERATURE_COEFS = 22; /*this is a guess */
const int MAX_TEMPERATURE_IND_VARS = 22; /*this is a guess */
const int IDXARRAY_SIZE = 9; /* set to match input file */

/* constants used for interpolating moderate unaggregated from aggregated */
const double ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0;
const double HALF = ONE / TWO;
const double THIRD = ONE / THREE;
const double TWOTHIRDS = TWO / THREE;
const double FOURTHIRDS = THIRD + ONE;
const double FIVETHIRDS = TWOTHIRDS + ONE;
const double FOURTH = ONE / FOUR;
const double THREEFOURTHS = THREE / FOUR;
const double FIVEFOURTHS = FOURTH + ONE;
const double SEVENFOURTHS = THREEFOURTHS + ONE;

const int NUM_UNAGG_PRODUCTS = 8; // number of products interpolated

// indices into arrays holding the interpolated products

enum InterpArrayIndexType {
    I_LON_IDX = 0,
    I_LAT_IDX,
    I_HT_IDX,
    I_SATA_IDX,
    I_SATZ_IDX,
    I_SOLA_IDX,
    I_SOLZ_IDX,
    I_RANGE_IDX
};

enum ViirsGeoLimitCheck {
    LOWEREND = 0,
    UPPEREND,
    MAX_LIMIT_CHECK
};

// Start of aggregation zones (Moderate aggregated)
const int START_AGG_ZONE1 = 0;
const int START_AGG_ZONE2 = 640;
const int START_AGG_ZONE3 = 1008;
const int START_AGG_ZONE4 = 1600;
const int START_AGG_ZONE5 = 2192;
const int START_AGG_ZONE6 = 2560;

// End of aggregation zones (Moderate aggregated)
const int END_AGG_ZONE1 = START_AGG_ZONE2 - 1;
const int END_AGG_ZONE2 = START_AGG_ZONE3 - 1;
const int END_AGG_ZONE3 = START_AGG_ZONE4 - 1;
const int END_AGG_ZONE4 = START_AGG_ZONE5 - 1;
const int END_AGG_ZONE5 = START_AGG_ZONE6 - 1;
const int END_AGG_ZONE6 = 3200;

// Start of aggregation zones (Moderate unaggregated)
const int START_UNAGG_ZONE2 = 0;
const int START_UNAGG_ZONE3 = 736;
const int START_UNAGG_ZONE5 = 4288;

// End of aggregation zones (Moderate unaggregated)
const int END_UNAGG_ZONE2 = START_UNAGG_ZONE3 - 1;
const int END_UNAGG_ZONE4 = START_UNAGG_ZONE5 - 1;
const int END_UNAGG_ZONE5 = 5023;

typedef struct {
    int band_number; //  band no. to geolocate (0=ideal band)
    char pad0[4];
    double latch_to_center;
    double t_reset; //  time to reset sample at beginning of frame

    unsigned short N_samp[MAX_BAND_NUMBER_PLUS_ONE];
    char pad1[6];
    double focal_length[MAX_BAND_NUMBER_PLUS_ONE];
    double det_space_track[MAX_BAND_NUMBER_PLUS_ONE];
    double det_space_scan[MAX_BAND_NUMBER_PLUS_ONE];
    double DNB_space_track[32];
    double DNB_space_scan[32];
    double det_position[MAX_BAND_NUMBER_PLUS_ONE][2];

    //  band center position 
    double band_position[MAX_BAND_NUMBER_PLUS_ONE];

    double earth_view_delay;
    double detector_sampling_rate;
    double scan_length;
    int agg_zone_bounds[NUM_AGG_ZONES];
    int DNB_aggregation[32][2];
    int DNB_ag_zone_bounds[64][3];
    char implicit_pad1[4];
    double scan_ang_coef_tel;
    double scan_ang_coef_mirr;
    double scan_ang_offsets[2];
} focal_plane_geometry_struct;

typedef struct {
    double enc_scale; // Scale to convert from 14-bit to 16-bit
    double mirr_abs_limit[MAX_LIMIT_CHECK]; // mirror encoder time absolute limits
    double mirr_del_limit; // mirror encoder time delta limits
    double tel_abs_limit[MAX_LIMIT_CHECK]; // telescope encoder time absolute limits
    double tel_del_limit; // telescope encoder time delta limits
    int sample_impulse_mirr; // encoder sample period in impulses
    int sample_impulse_tel; // encoder sample period in impulses
    int A_bit_adj[2]; // Encoder adjust [even SOS/odd SOS]
    int B_HAM_adj[2]; // HAM enc adjust [side A/side B]
    double t_encoder; // encoder time scale factor
} mirror_preparation_struct;

typedef struct {
    double mirr_side1_range[2]; //  mirror side 1 angle range (radians)
    double alpha;
    double beta;
    double gammaa;
} mirror_model_struct;

typedef struct {
    double position_abs_limit[2]; //  orbit position absolute limits
    double position_mag_limit[2]; //  orbit position magnitude limits
    double velocity_abs_limit[2]; //  orbit velocity absolute limits
    double velocity_mag_limit[2]; //  orbit velocity magnitude limits
    double ang_mom_limit[2]; //  angular momentum magnitude limits
    double ang_mom_z_limit[2]; //  angular momentum Z component
    double orbit_consistency; //  orbit pos./vel. consistency limit
} orbit_validation_params_struct;

typedef struct {
    double attitude_abs_limit[2]; //  attitude angle absolute limits
} attitude_valid_struct;

typedef struct {
    double T_inst2sc[3][3]; //  instrument to spacecraft
    double T_mirr2inst[3][3]; //  mirror to instrument frame
    double T_aft2inst[3][3]; //  aft optics to instrument frame
    double T_inst2SD[3][3]; //  instrument to solar diffuse
    double T_tel2inst[3][3]; //  telescope to instrument
    int rpy_count; // number of rpy adjustments that are included
    char rpy_times[MAX_RPY_COUNT][RPY_TIMECODE_SIZE]; // timestamp for each rpy adjustment
    double rpy_tai[MAX_RPY_COUNT]; // timestamp converted to TAI double
    double rpy_values[MAX_RPY_COUNT][3]; // rpy value for each adjustment
} internal_coord_trans_struct;

typedef struct {
    //  Number of thermistors currently used in thermal correction. 
    int num_thermistor;

    //  Identifier for each thermistor. 
    unsigned char thermistor_id[MAX_THERMISTOR_ID_LEN][NUM_THERMISTORS];

    // Coefficient array for determination of instrument
    // temperature from thermistor readings. 
    char implicit_pad1[4];
    double thermistor_coeffs[NUM_THERMISTORS][6];
} thermistor_parameter_struct;

// Filename,VcstGeoParameters
// Type,GEO_param_struct

typedef struct {
    unsigned char revision[10];
    char pad1[6];
    focal_plane_geometry_struct geometry_params;
    mirror_preparation_struct mirror_prep_params;
    mirror_model_struct mirror_model;
    internal_coord_trans_struct coord_trans;
    thermistor_parameter_struct thermistor_params;

    double Mag[3][3];
    double basis_in[3][3];
    double basis_out[3][3];


    // polynomial coefficients for HAM encoder-to-angle conversion
    double poly_coef_mirr[MAX_POLY_DEGREE_PLUS_ONE];

    // polynomial coefficients for telescope encoder-to-angle conversion
    double poly_coef_tel[NUM_ELECTRONICS_SIDE][MAX_POLY_DEGREE_PLUS_ONE];

    double tel_ref;
    double min_cos_view; //  minimum cos(SensorZenAngle)
    //  for near limb test

    int band_type;

    int num_detectors; //  no. of band_number detectors

    // degree of polynomial used to calculate tel and mirr pos   
    int poly_degree;

    unsigned short N_frame; //  no. frames per scan

    char pad2[2];


} GEO_param_struct;


#endif
