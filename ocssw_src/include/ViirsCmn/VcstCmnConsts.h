/*****************************************************************************
 *
 * NAME:  VcstCmnConsts.h
 *
 * DESCRIPTION: Defines constants used for Common Geolocation and Calibration.
 *
 * Adapted as a consolidation of several IDPS header files published by
 * Raytheon Company.
 *
 * REFERENCES:
 * Meeus, Jean, "Astronomical Algorithms, 2nd Edition," Willman-Bell Inc,
 * Richmond VA, 1998, 477 pp.
 *
 *
 *****************************************************************************/

#ifndef VcstCmnGeoConsts_h
#define VcstCmnGeoConsts_h

#include <math.h>
#include <string>
#include <VcstViirsStructs.h>

const int TEL_START_ENC_NOMINAL_NPP = 30940;
const int HAM_START_ENC_NOMINAL_0_NPP = 10514;
const int HAM_START_ENC_NOMINAL_1_NPP = 10513;
const int SCAN_ENCDR_START_MAX_NPP = 32767;

const int TEL_START_ENC_NOMINAL_J1 = 31002;
const int HAM_START_ENC_NOMINAL_0_J1 = 10579;
const int HAM_START_ENC_NOMINAL_1_J1 = 10579;
const int SCAN_ENCDR_START_MAX_J1 = 32767;

const int TEL_START_ENC_NOMINAL_J2 = 30814;
const int HAM_START_ENC_NOMINAL_0_J2 = 10465;
const int HAM_START_ENC_NOMINAL_1_J2 = 10465;
const int SCAN_ENCDR_START_MAX_J2 = 32767;

const signed char SCE_INVALID = -1;
const char SCE_A_SIDE_ON = 0; //DR4759 scan control electronics side
const char SCE_B_SIDE_ON = 1;
const int SCAN_SYNC_FAILURE = 1; //DR4795 tel/ham sync failure
const int SCAN_SYNC_NORMAL = 0;
const int SCAN_SYNC_UNDETERMINED = -1;

const int VCST_SUCCESS = 0;
const int VCST_FAIL = 1;
const int VCST_WARNING = 2;
const int VCST_STOP = 3;
const int VCST_INIT_FAIL = 4;
const int VCST_CROSSGRAN_FAIL = 5;
const int VCST_TILE_FAIL = 6;
const int VCST_NO_DATA = 7;
const int VCST_THREAD_LAUNCHED = 8;
const int VCST_GEO_WARNING = 12;
const int VCST_GEO_OLD_TLE = 14;

const int VIIRS_SCANS = Number_of_Scans;
const int SC_DIARY_RECORDS = SC_Diary_Records;
const int SC_DIARY_RECORDS_10HZ = SC_Diary_Records_10Hz;
const int SC_DIARY_RECORDS_1HZ = SC_Diary_Records_1Hz;
const int VCST_BANDS = 22;

// L1A fills
const double VCST_DOUBLE_FILL = -999.9e0;
const float VCST_FLOAT_FILL = -999.9;
const long long VCST_INT64_FILL = -999;
const int VCST_INT32_FILL = -999;
const short VCST_SHORT_FILL = -999;
const short VCST_USHORT_FILL = 32767;
const signed char VCST_BYTE_FILL = -1;
const unsigned char VCST_UBYTE_FILL = 255;

// Ellipsoid Intersection Failed
const double ELLIPSOID_FLOAT64_FILL = -999.9e0;
const float ELLIPSOID_FLOAT32_FILL = -999.9;

// Algorithm could not compute the pixel because of a software problem
// ie. couldn't converge
const double ERR_FLOAT64_FILL = -999.9e0;
const float ERR_FLOAT32_FILL = -999.9;

//Missing - C3S provided a fill value or AP missing
const double MISS_FLOAT64_FILL = -999.9e0;
const float MISS_FLOAT32_FILL = -999.9;
const long long MISS_INT64_FILL = -999;
const int MISS_INT32_FILL = -999;
const unsigned char MISS_UINT8_FILL = 255;

// The pixel was not to be computed because it is not applicable to situation
const double NA_FLOAT64_FILL = -999.9e0;
const long long NA_INT64_FILL = -999;

// Value Does Not Exist/Missing Scan
const double VDNE_FLOAT64_FILL = -999.9e0;

// Fill Tests
const long long INT64_FILL_TEST = -990;
const double FLOAT64_FILL_TEST = -999.0e0;
const float FLOAT32_FILL_TEST = -999.0;

// Error codes for Cmn Geo

const int CMNGEO_WARNING = 12;
const int CMNGEO_OLD_TLE = 14;

// Band groups
const int NUM_DNB_BANDS = 1;
const int NUM_REFL_375M_BANDS = 3;
const int NUM_REFL_750M_BANDS = 11;
const int NUM_REFL_750M_SG_BANDS = 5;
const int NUM_EMISS_375M_BANDS = 2;
const int NUM_EMISS_750M_BANDS = 5;
const int NUM_EMISS_750M_SG_BANDS = 4;
const int NUM_EMISS_750M_DG_BANDS = 1;
const int NUM_REFLECTIVE_BANDS = 15;
const int NUM_REFLECTIVE_BANDS_WITHOUT_DNB = NUM_REFLECTIVE_BANDS
        - NUM_DNB_BANDS;
const int NUM_EMISSIVE_BANDS = 7; // m12,m13,m14,m15,m16,i4,i5
const int NUM_REFLECTIVE_BT = 3; // M7,M8,M10
const int NUM_VIS_BANDS = 9;
const int NUM_SMIR_BANDS = 8;
const int NUM_LW_BANDS = 4;

const int EMISSIVE_BANDS_INDEX[NUM_EMISSIVE_BANDS] ={3, 4, 16, 17, 18, 19, 20};
const int VIS_BAND_INDEX[NUM_VIS_BANDS] = {0, 1, 5, 6, 7, 8, 9, 10, 11}; /* I1,2 M1-7 */
const int SMIR_BAND_INDEX[NUM_SMIR_BANDS] = {2, 3, 12, 13, 14, 15, 16, 17}; /* I3,4 M8-13*/
const int LWIR_BAND_INDEX[NUM_LW_BANDS] = {4, 18, 19, 20}; /* I5 M14-16 */

const int CMNGEO_MIN_MAX_DIM = 2;

const int JPL_EPHEM_ROW = 460;
const int JPL_EPHEM_COL = 826;
const int JPL_IPT_ROW = 13;
const int JPL_IPT_COL = 3;

const int SD_MATRIX_COL = 3;
const int SD_MATRIX_ROW = 3;

const int SDSM_COEF = 6;
const int SDSM_DETECTORS = 8;
const int SDSM_SAMPLES = 5;
const int BB_THERMISTORS = 6;

const unsigned char SDSM_POS_HOME = 0;
const unsigned char SDSM_POS_HOME_RAW = 0;
const unsigned char SDSM_POS_SD_VIEW = 1;
const unsigned char SDSM_POS_SD_VIEW_RAW = 28;
const unsigned char SDSM_POS_SUN_VIEW = 2;
const unsigned char SDSM_POS_SUN_VIEW_RAW = 67;

const double TAI93_TAI58_SEC = 1104537627.0;
const double AU_TO_METERS = 149597870700.0;

const int MOON_INSIDE_SV_KOB = 1;
const int MOON_OUTSIDE_SV_KOB = 0;
const int UNDETERMINED_MOON_IN_SV_KOB = -1;
const int TRK_UPPER = 0; // track upper limit index
const int TRK_LOWER = 1; // track lower limit index
const int SCN_UPPER = 2; // scan upper limit index
const int SCN_LOWER = 3; // scan lower limit index

const int MAX_STRING_LENGTH = 255;

const int HAM_ENCDR_START_NOT_NOMINAL = 1; //DR4767 non nominal ham encoder start
const int HAM_ENCDR_START_IS_NOMINAL = 0;
const int HAM_ENCDR_START_UNDETERMINED = -1;
const int TEL_ENCDR_START_NOT_NOMINAL = 1;
const int TEL_ENCDR_START_IS_NOMINAL = 0;
const int TEL_ENCDR_START_UNDETERMINED = -1;

//DR4777 bit 5 for non nominal telescope encoder start
const unsigned char SDR_SCAN_QUALITY_TEL_START_IS_NOMINAL = 0x00; //00000000
const unsigned char SDR_SCAN_QUALITY_TEL_START_NOT_NOMINAL = 0x10; //000x0000
//DR4795 bit 4 for scan sync failure
const unsigned char SDR_SCAN_QUALITY_SCAN_SYNC_NORMAL = 0x00; //00000000
const unsigned char SDR_SCAN_QUALITY_SCAN_SYNC_FAILURE = 0x08; //0000x000

const unsigned char VRDR_ENG_ENCODER_VALID_DATA = 0;
const unsigned char VRDR_ENG_ENCODER_NO_DATA = 1;
const unsigned char VRDR_ENG_ENCODER_MASK = 0x1;
const unsigned char VRDR_ENG_ENCODER_SHIFT = 0;

// size of a quaternion
const int QUAT_SIZE = 4;

// size of a vector
const int VEC_SIZE = 3;

// time in microseconds of a full viirs scan
const double VIIRS_SCAN_RATE = 1.786400;

// Angular values are from Y3261
const double START_OF_EV_DEG = -56.04;
const double MIDDLE_OF_SV_DEG = -65.40;
const double MIDDLE_OF_SD_DEG = 157.42;
const double MIDDLE_OF_BB_DEG = 100.0;

// the solar diffuser center view time offset from the start of the
// scan time.  start of scan time is the beginning of the earth view
// (VIIRS_SCAN_RATE / 360.0) * (MIDDLE_OF_SD_DEG - START_OF_EV_DEG);

const double SD_VIEW_CENTER_OFFSET_FROM_START_OF_SCAN = 1.059236;

// the space view center view time offset from the start of the
// scan time.  start of scan time is the begining of the earthview
// Note the correct space view to consider is actually just prior to
// to the earth view and the offset becomes negative
// (VIIRS_SCAN_RATE / 360.0) * (MIDDLE_OF_SV_DEG - START_OF_EV_DEG);

const double SV_VIEW_CENTER_OFFSET_FROM_START_OF_SCAN = -0.046446;

// the black body center view time offset from the start of the
// scan time.  start of scan time is the beginning of the earth view
// (VIIRS_SCAN_RATE / 360.0) * (MIDDLE_OF_BB_DEG - START_OF_EV_DEG);

const double BB_VIEW_CENTER_OFFSET_FROM_START_OF_SCAN = 0.774305;

// Number of seconds in 30 minutes
const double SEC_IN_30MIN = 1800.0e+0;

// Fraction of a second to check if a leap second
// has occurred in UT1mUTC data.
//TODO: Randy's code sets this to .5
// and does the interpolation differently for leap seconds!
const double LEAPSEC_DETECT = 0.65e0;

// Difference between estimated lat and calculated lat
const double MAX_LAT_DIFF = 1.e-10;

// Instrument frame to ECR transformation matrix
const double tMtrxInst2SC[VEC_SIZE][VEC_SIZE] = {
    { 1.0, 0.0, 0.0},
    { 0.0, 1.0, 0.0},
    { 0.0, 0.0, 1.0}
};

// Constant for the value 1000.0
const double KILO = 1000.0e0;

// Constant for the value 1/1000.0
// (used for meters to kilometers conversion)
const double INV_KILO = 1.0e-3;

// Maximum number of seconds of difference between 2 packets.
// If difference exceeds this amount, interpolate for missing packets
const double PKT_TIME_TOL = 1.8e0;

// Minimum number of points before and after granule start/stop time
// that CmnGeo would like to have
const int MIN_POINTS = 2;

// Approximately 1cm when muliplied by Equitoral radius
const double APPROX_1CM = 1.5e-09;

// The following constants are Flag definitions

// Flag for converting Eci to Ecr
const int ECI2ECR = 0;

// Flag for converting Ecr to Eci
const int ECR2ECI = 1;

// Flag for calculating Sun values
const int SUN_FLAG = 0;

// Flag for calculating Moon values
const int MOON_FLAG = 1;

// Constant for NOVAS-C designation of Earth
const int NOVASC_EARTH = 3;

// Constant for NOVAS-C designation of Sun
const int NOVASC_SUN = 10;

// Constant for NOVAS-C designation of Moon
const int NOVASC_MOON = 11;

// Constant for NOVAS-C function equ2hor.  This selects the
// standard atmosphere refraction model
const short NOVASC_REFRACT = 0;

// Max number of eclipses to cache
const short PRO_GEO_MAX_ECLIPSES = 3;

// Time of the New Moon for Jan 6, 2000, in TJD/JDE (Terrestrial Dynamic 
// Julian Day = Julian Day Ephemeris, in Meeus' terminology).  From
// Meeus, pg. 349 (eq 49.1).
const double TJD_JAN2000_NEWMOON = 2451550.09766;

// Length of the mean synodic period of the Moon, i.e., the mean interval 
// between two consecutive New Moons, in days.  From Meeus, pg. 349
// (eq. 49.1).
const double MOON_PERIOD_DAYS = 29.530588861;

// Angular semidiameter of the Sun at a distance of 1 AU, in radians
// (equiv. to Meeus value of 959".63).  See Meeus, pg. 389.
const double SUN_SEMIDIAM = 4.6524e-3;

// Radius of the Moon, derived from the Meeus' value of "k" = 0.272481, 
// which is the ratio of the Moon's mean radius to equatorial radius of
// the Earth, times the equatorial radius of the Earth from ProCmnPhysConst.h
// Units are meters.  See Meeus, pg. 390
const double LUNAR_RADIUS = 1.73792e6;

// The following are from ProCmnMathConsts.h

// index for X component in an array
const int X_COM = 0;

// index for Y component in an array
const int Y_COM = 1;

// index for Z component in an array
const int Z_COM = 2;

// index for Roll component in an attitude array
const int ROLL_COM = 0;

// index for Pitch component in an attitude array
const int PITCH_COM = 1;

// index for Yaw component in an attitude array
const int YAW_COM = 2;

// index for I component in quaternion array
const int I_COM = 0;

// index for J component in quaternion array
const int J_COM = 1;

// index for K component in quaternion array
const int K_COM = 2;

// index for W component in quaternion array
const int W_COM = 3;

// Minimum quaternion magnitude
const double QUAT_MIN_MAG = 9.999999e-1;

// Maximum quaternion magnitude
const double QUAT_MAX_MAG = 1.000001e0;

// Value of PI/2 in Float32 format
const float FLOAT32_PIO2 = 1.570795;

// Value of PI in float format
const float FLOAT32_PI = 3.141591;

const double PI = M_PI; // math constant pi
const double PIO2 = M_PI_2; // pi/2
const double PIO4 = M_PI_4; // pi/4
const double TREPIO2 = 3.0L * M_PI / 2.0L; // 3*pi/2
///const double  TWOPI =    2.0L*M_PI;           // 2*pi

// used to convert degrees to radians
//const double  DEG2RAD =  M_PI/180.0L;         // (pi/180)

// used to convert radians to degrees
//const double  RAD2DEG =  180.0L/M_PI;         // (180/pi)

// used to convert degrees to arcsec
const double DEG2ARCSEC = 3600.0L; // seconds in a degree

// Universal Gas Constant - used to compute virtual temp.  J/(kg*K)
const float DRYGAS = 287.05;

// Average Atmospheric Pressure at Sea Level
const double AVG_PRESS_SEALVL = 1013.25e0;

// Gravity (m/s^2)
const double GRAVITY = 9.80665;

// Radius of the Area Weighting Reference Sphere (m)
// This is NOT the radius of the earth and should not be used for any
// calculations other than by Area Weighting to determine if a pixel
// is near the pole.
const double EARTH_RADIUS_METERS = 6371007.181;

// The following constant definitions are the WGS84 Earth Ellipsoid
// constants.  The main reference for the WGS84 ellipsoid is
// NIMA Physical Geodesy web page: 164.214.2.59/GandG/wgs-84/egm96.htm
// See the updated page with the link to the NIMA TR8350.2 document:
// http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
// The Flattening Factor(f) is computed from f = 1/298.257223563.
const double EQUAT_RAD = 6.37813700000000e+6; // Equatoral rad., meters WGS84
const double POLAR_RAD = 6.35675231424518e+6; // Polar radius, meters WGS84
const double ECCEN_SQ = 6.69437999014132e-3; // Eccentricity Squared WGS84
const double FLATFAC = 3.35281066474748071e-3; // Flattening Factor WGS84

// Earth Gravitational Parameter mu (G*M) in m^3 per s^2 from WGS84.
// The central term in the Earth's gravitational field (GM) is known with
// much greater accuracy than either 'G', the universal gravitational
// constant, or 'M', the mass of the Earth.  The refined value accounting
// for the mass of the atmosphere is (3986004.418 +/- 0.008) e+8 m^3/s^2.
// GPS OCS applications take advantage of the improved value, however,
// Section 3.2.3.2 and the ICD-GPS-200 recommends the original WGS 84 GM
// value to avoid introduction of error to the GPS user.
const double EARTH_GRAV_mu = 3.986005000e+14;

// USAF Orbit Analyst Manuals, circa 1978.  1 - eccen_sqr
const double DETIC2CENTRIC = 9.93305620009859e-1;
const double CENTRIC2DETIC = 1.00673949674228e+0;

// The following constant definitions are for time conversions

const double TAI2IET = 1.0e+06; // Conversion factor
const double MIN_IN_HOUR = 60.0e+0; // Number of minutes in an hour
const double SEC_IN_HOUR = 3600.0e+0; // Number of seconds in an hour
const double MJD_CONV_FAC = 2.4000005e+6; // Factor to convert AJD to MJD
const double SEC_IN_DAY = 8.64e+04; // Number of seconds in a day
const double UJD58 = 2.43620450e+06; // Jan 1 1958  UJD format
const double JAN012030 = 2.272147232e+09; // Jan 1 2030  TAI format
const double TJD_CONV_FAC = 32.184e+0; // Factor to convert TAI to TJD
const double DEG_IN_HOUR = 15.0e+0; // Number of degrees in an hour

// The following constant definitions are for polarstereographic dataset
const double MINUS30 = -0.523598775598299e0; // -30 degrees in radians
const double PLUS30 = 0.523598775598299e0; // 30 degrees in radians

// The following constant definitions are for nwp ancillary granulation
// declare constant for calculation of water vapor mixing ratio (r)
const double GAS = 621.97; //-- ratio of the molecular weight
//-- of water vapor to dry air
//-- in units of grams/kilogram

// declare just a few of the more popular of the twenty SI prefixes
const double MICRO = 0.000001; //-- scale by 1/1000000th
const double MILLI = 0.001; //-- scale by 1/1000th
const double CENTI = 0.01; //-- scale by 1/100th
const double DECI = 0.1; //-- scale by 1/10th
const double DEKA = 10.0; //-- scale by 10x
const double HECTO = 100.0; //-- scale by 100x

// Kelvin/Celsius conversion factor
const double TCOEFF = 273.15;

// Constant used to generate surface reflectance
// multiplier to convert pascal to atmospheres (1 atm = 101325 pascal)
const float PRESS_CONV = 1.0 / 101325.0;

// Standard Atmosphere Surface Pressure
const double STDPSL = 1013.0;

// Moist air adiabatic lapse rate is 6.5 K/Km (equivalent to 6.5 C/Km)
// Converted value would be .0065 C/m
const double MOIST_AIR_LAPSE_RATE = 6.5 / 1000;

// Constant used to convert atm-cm to Dobson units
const double ATM_CM2DOBSON = 1000.0;

/*******************************************************************************
 *
 * NAME: maptrig_main
 *
 * DESCRIPTION: Several constants always needed for map conversions are defined
 * in a way that makes them global.  Each function that uses these constants
 * does an include on "maptrig_func.h"
 *
 * REFERENCES: none
 *
 * LIMITATIONS: none
 *
 * REVISION/EVENT HISTORY:
 * DATE        PR#      AUTHOR            Build    DESCRIPTION
 * ---------   ---      ------            -----    -----------
 * 01MAR2010   021991   J. Gibbs          SenChar  Added defines
 *
 * NOTES (MISCELLANEOUS) SECTION:
 *
 *******************************************************************************/

const double pi = 3.14159265358979e+0; /*  widely known math constant pi */
const double pio2 = 1.57079632679490e+0; /*  pi/2    */
const double trepio2 = 4.71238898038470e+0; /*  3*pi/2  */
const double pio4 = 7.85398163397448e-1; /*  pi/4    */
const double twopi = 6.28318530717960e+0; /*  2*pi    */
const double deg2rad = 1.74532925199433e-2; /*  converts degrees to radians */
const double rad2deg = 5.72957795130822e+1; /*  converts radians to degrees */

/***************************************************************************
 WGS84 Earth Ellipsoid constants needed by earth_radius and other functions
 The main reference for the WGS84 ellipsoid is NIMA Physical Geodesy
 source document on internet:
 http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
 Equatorial radius and Inverse Flattening factor on page 3-2.

 NOTE: These constants can have 15 significant digits, or as many digits
 as desired, because NIMA specifically defined the WGS84 Equatorial radius
 to exactly 6,378,137 meters, with no fraction after that. They also
 defined the number 1/f (see flatinv below) with as many zeroes behind it
 as you want to add. So all the numbers associated with the WGS84 Ellipsoid
 can be calculated to as many significant digits as desired.

 NOTE: Inverse of the Earth Flattening factor is from the definition
 of WGS84, and is 1/f. (It is inside a comment block because it is not
 directly used, but the other constants shown here are calculated from
 it.)
 double flatinv=298.257223563000e+0;

 NOTE: WGS84 Earth flattening factor is f. (It is inside a comment block
 because it is not directly used.)
 double flatten=3.35281066474748e-3;

 NOTE: The following constants are directly used by earth_radius, and/or
 other functions.

 equatorial radius, in kilometers, is from the definition of WGS84
 equatorial radius in meters
 polar radius in meters
 eccentricity squared is from e^2 = f(2 - f)
 detic2centric is explained in earth_radius.cpp
 centric2detic is the inverse of the above
 delta is explained in the comments in earth_radius.cpp
 **********************/

const double eq_rad_km = 6.37813700000000e+3; /* equatorial radius, KM */
const double eq_radm = 6.37813700000000e+6; /* equatorial radius, meters */
const double pole_radm = 6.35675231424518e+6; /* polar_radius, meters */
const double eccen_sqr = 6.69437999014132e-3; /* e^2 = f(2 - f) */
const double detic2centric = 9.93305620009859e-1; /* 1 - e^2 */
const double centric2detic = 1.00673949674228e+0; /* 1 / (1 - e^2) */
const double delta = 6.73949674227643e-3; /* ( 1/(1-f)^2 ) - 1 */

#endif
