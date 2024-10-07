/******************************************************************************
 *  NAME: VcstCmnGeo.h
 *
 *  DESCRIPTION: Object class that serves as the intermediary between the NASA
 *  L1A file format and the IDPS Common Geolocation interfaces. It reads L1A
 *  navigation data and stores it into lists of spacecraft diary records
 *  compatible with the VcstCmnGeo object class.  It initializes common
 *  geolocation for use by other functions.
 *
 *  Created on: December 15, 2014
 *      Author: Sam Anderson, VCST
 *
 ******************************************************************************/

#ifndef VcstCmnGeo_H_
#define VcstCmnGeo_H_

#include <cmath>
#include <novas.h>
#include <VcstTime.h>
#include <VcstCmnGeoStructs.h>
#include <VcstCmnLutInputItem.h>

#include <hdf5.h>
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

struct viirsCmnLutPtrs {
    // Pointers to LUT data structures
    InstrumentLutType* scDataLut;
    CmnGeoPolarUT1Type* pwUT1Ptr;
    JPLEphemLutType* jplEphemPtr;
    SAALutType* saaCoeff;
    CmnGeoParamLutType* paramLut;
    proSdrViirsCalQALUT* QALUT;
};

// Structure of LUT pointers used by Common Geolocation

class VcstCmnGeo {
public:

    // numeric constants
    static constexpr double MSECPERDAY = 86400000000.0; //24*60*60=86400 million
    static constexpr double TAI93_TAI58_SEC = 1104537627.0;
    static constexpr double RAD_TO_DEG = 180.0 / M_PI;

    static const std::string SN_APID_0; // = "CRITICAL";
    static const std::string SN_APID_8; // = "ADCS_HKH";
    static constexpr int ADSTATE_BYTE = 169;
    static constexpr int ADSTATE_BITS = 3;
    static constexpr int ADSTATE_NUM_OF_BITS = 3;
    static constexpr int ADMANDONE_BYTE = 255;
    static constexpr int ADMANDONE_BIT = 7;
    static constexpr int ADFFTID_BYTE = 339;
    static constexpr int SCAN_ALL = -1;

    static constexpr bool MEAN_SID = true; // Flag for getting Mean Sidereal Time.
    static constexpr bool APP_SID = false; // Flag for getting Apparent Sidereal Time.

    /**
     *  VIIRS platform identification string
     */

    std::string platform_;

    /**
     *  Global data used by geolocation
     */

    int act_scans_;
    int sc_diary_records_;

    int scan_number_[VIIRS_SCANS];
    unsigned char ham_side_[VIIRS_SCANS];
    unsigned char sensor_mode_[VIIRS_SCANS];
    unsigned char se_a_anlg_pwr_on_[VIIRS_SCANS];
    unsigned char se_b_anlg_pwr_on_[VIIRS_SCANS];
    unsigned char dp_servo_in_use_[VIIRS_SCANS];
    unsigned char se_a_tele_pos_known_[VIIRS_SCANS];
    unsigned char se_b_tele_pos_known_[VIIRS_SCANS];
    unsigned char se_a_mtrs_stopped_[VIIRS_SCANS];
    unsigned char se_b_mtrs_stopped_[VIIRS_SCANS];
    unsigned char encoderFlag_[VIIRS_SCANS];
    unsigned char es_se_a_teleham_scansyn_[VIIRS_SCANS];
    unsigned char es_se_b_teleham_scansyn_[VIIRS_SCANS];
    unsigned char cp_blk_pwr_sel_[VIIRS_SCANS];
    unsigned char dp_nonrdt_fpie_pwr_[VIIRS_SCANS];
    unsigned short ham_start_enc_[VIIRS_SCANS];
    unsigned short tel_start_enc_[VIIRS_SCANS];

    unsigned short ham_encoder_[VIIRS_SCANS][Encoder_Reading];
    unsigned short tel_encoder_[VIIRS_SCANS][Encoder_Reading];

    short electronics_side_[VIIRS_SCANS];
    bool scan_sync_failure_[VIIRS_SCANS];
    bool tel_start_not_nominal_[VIIRS_SCANS];
    bool sensor_not_nominal_[VIIRS_SCANS];
    bool scan_data_missing_[VIIRS_SCANS];
    short scan_sync_failure_cnt_;
    int scanSyncFailCt_;

    double start_TAI93sec_[VIIRS_SCANS];
    double end_TAI93sec_[VIIRS_SCANS];
    double ev_mid_TAI93sec_[VIIRS_SCANS];
    double sv_mid_TAI93sec_[VIIRS_SCANS];
    double bb_mid_TAI93sec_[VIIRS_SCANS];
    double sd_mid_TAI93sec_[VIIRS_SCANS];

    double start_TAI58sec_[VIIRS_SCANS];
    double end_TAI58sec_[VIIRS_SCANS];
    double mid_TAI58sec_[VIIRS_SCANS];
    double midGranTAI58sec_;

    float sc_solar_zenith_[VIIRS_SCANS];
    float sc_solar_azimuth_[VIIRS_SCANS];

    float solar_zenith_[VIIRS_SCANS];
    float solar_azimuth_[VIIRS_SCANS];
    float earth_sun_distance_[VIIRS_SCANS];
    float solar_j2000_[VIIRS_SCANS][VEC_SIZE];
    float solar_inst_[VIIRS_SCANS][VEC_SIZE];

    float att_quat_ev_[VIIRS_SCANS][QUAT_SIZE];
    float att_ang_[VIIRS_SCANS][VEC_SIZE];
    float orb_pos_ev_[VIIRS_SCANS][VEC_SIZE];
    float orb_vel_ev_[VIIRS_SCANS][VEC_SIZE];

    float att_quat_sd_[VIIRS_SCANS][QUAT_SIZE];
    float orb_pos_sd_[VIIRS_SCANS][VEC_SIZE];
    float orb_vel_sd_[VIIRS_SCANS][VEC_SIZE];

    float earth_moon_distance_[VIIRS_SCANS];
    float lunar_j2000_[VIIRS_SCANS][VEC_SIZE];
    float lunar_inst_[VIIRS_SCANS][VEC_SIZE];
    bool moon_in_sv_kob_[VCST_BANDS][VIIRS_SCANS];

    double earth_sun_distance_avg_;

    // Pointer to structure of LUT pointer

    viirsCmnLutPtrs pLut_;
    /*
     *  Constructor
     */

    VcstCmnGeo();

    /*
     *  Destructor
     */

    ~VcstCmnGeo();

    /**
     *  Initialize navigation based on granule sequence enum
     */

    int initialize(GRAN_SEQ_ENUM gran_seq);

    /**
     *  Determine spacecraft maneuver state based on the information
     *  Bus Critical Telemetry data and ADCS Housekeeping Telemetry.
     */

    static const std::string getManeuverStr(const BusCritTlmType * busCritTlm,
            const AdcsHskTlmType * adcsHskTlm);

    /**
     * This function will find a point in the ephatt_ points map
     * field for a specified time.  The point that is found is
     * used to determine the location on the earth.  The Roll, Pitch,
     * and Yaw produced by this method are in arcseconds.
     */

    int satPosAtt(const double inTAI, CmnGeoEAPointType& outPt,
            bool performOrbitCheck = false);

    /**
     * This function calculates the x, y, z vector to the Sun or Moon
     * in spacecraft frame at the location of the satellite.
     */

    int vectorAtSat(const CmnGeoEAPointType& inPt, const int vecFlag,
            double vecSat[VEC_SIZE]);

    /**
     * This function will calculate the satellite azimuth and
     * zenith angles from the sample position
     */

    void calcSatAzmZen(const double inDLat, const double inLon,
            const double inSamp2Sat[VEC_SIZE], double& outSatAzm,
            double& outSatZen);

    /**
     * This function will calculate the geodetic latitude/longitude and
     * satellite azimuth/zenith angles.  It will combine the sensor
     * exit vector with the roll, pitch, & yaw of the point.
     */

    int ellipIntersect(const CmnGeoEAPointType& inPt,
            const double inRotMat[VEC_SIZE][VEC_SIZE],
            const double viewVec[VEC_SIZE], double& dLat, double& lon,
            double& satAzm, double& satZen, double& range);

    /**
     *  Compute solar and lunar angles at lat and lon
     */

    int sunMoonAngles(const double inTAI, const double inLat,
            const double inLon, const int smFlag, double& sunAzm,
            double& sunZen, double& moonAzm, double& moonZen, double& moonPhase,
            double& moonIFrac);

    /**
     *  This function will calculate the azimuth and zenith angles
     *  to the Sun from the intersection on the ellipsoid.
     */

    int sunAngles(const double inTAI, const float inLat, const float inLon,
            float& sunAzm, float& sunZen);

    int sunAngles(const double inTAI, const double inLat, const double inLon,
            double& sunAzm, double& sunZen);

    /**
     * This function will calculate the azimuth and zenith angles
     * to the Moon from the intersection on the ellipsoid.  It will
     * also calculate the Moon phase and Moon illumination.
     */

    int moonAngles(const double inTAI, const float inLat, const float inLon,
            float& moonAzm, float& moonZen, float& moonPhase, float& moonIFrac);

    /**
     * For a given location, estimate the number of single event upsets
     * per year.  A higher number indicates a greater danger of
     * increased radiation within the South Atlantic Anomaly (SAA) causing
     * earth location errors within optical encoders and false hits within
     * detector arrays.
     */

    int getSAAIntensity(const float inLat, const float inLon,
            float& outIntensity);

    /**
     * Check if a solar eclipse is occurring at the given time and
     * ground location, using detailed NOVAS routines to verify local
     * conditions.
     */

    int checkSolarEclipse(const double inTAI, const double inLat,
            const double inLon, bool& outEclipseOccurring);

    /**
     *  Write navigation data to netcdf file
     */

    int writeNavigation();

    /**
     *  Write geolocation data to netcdf file
     */

    int writeCmnGeo();


private:

    CmnGeoScEAType ephemAtt_;
    CmnGeoEAPointType sdPoint_[VIIRS_SCANS];
    CmnGeoEAPointType evPoint_[VIIRS_SCANS];
    CmnGeoEphAttType ephatt_;

    site_info locPixel_;
    body sun_;
    body earth_;
    body moon_;

    VcstTime* timeAPI_;
    polarWanderDataType *polarWanderDataPtr_;

    // Map to store LUT input item objects.

    std::map<std::string, VcstLutInputItem*> cmnGeoLutItems_;

    /**
     * State holder for the DEBUG message throttle.  Protect against
     * flooding the logs with DEBUG messages during a maneuver.
     */

    int dbgMsgThrottle_;

    /**
     *  Initialize NOVAS library
     */

    int initialize_novas();

    /**
     *  Read LUT data required for navigation
     */

    int initialize_LUT_data();

    /**
     *  Read L1A data and initialize navigation data based on granule sequence enum
     */

    int initialize_L1A_data(GRAN_SEQ_ENUM gran_seq);

    /**
     *  Check telencoder for abnormal status 
     */
    
    int Check_Tel_Start_Not_Nominal();

    /**
     *  Check scan sync status for all scans
     */

    int check_scan_sync();

    /**
     * This function forces the re-calculation of the Polar Wander values
     */

    int setPolarWander(const double a_tai);

    /**
     * Fills the ephatt_ structure with information specific to
     * the current granule.
     */

    int setupEphAtt(const CmnGeoScEAType * const sceaPtr);

    /**
     * This function will initialize the points map with byte-align
     * S/C E&A data.
     */

    int initPointsMap(const CmnGeoScEAType* inPtr);

    /**
     * This function will fill one point in the map of points in
     * the ephatt_ structure.
     */

    int fillEphAttPoint(CmnGeoEAPointType& pntNode);

    /**
     * This function will perform quadratic interpolation of the
     * three points.  It only fills in the ephemeris portion of
     * the data in outputPt.  Interpolation of the attitude
     * data is done elsewhere.
     */

    int quadraticInterpEph(const CmnGeoEAPointType* firstPt,
            const CmnGeoEAPointType* secondPt, const CmnGeoEAPointType* thirdPt,
            CmnGeoEAPointType* outputPt);

    /**
     * This function will convert position, velocity, or other vectors
     * in either ECR or ECI coords to ECI or ECR coords.
     */

    int convCoordSys(const double inTAI, const double inPos[VEC_SIZE],
            const double inVel[VEC_SIZE], const int convFlag,
            double outPos[VEC_SIZE], double outVel[VEC_SIZE]);

    /**
     * This function will convert a vector into geodetic latitude,
     * geocentric latitude, and longitude
     */

    void convVec2LatLon(const double inVec[VEC_SIZE], double& dLat,
            double& cLat, double& lon);

    /**
     * This function will calculate the ECR & ECI acceleration.
     * It uses the position information and the Earth Gravitational
     * parameter to determine the acceleration
     */

    void calculateAccel(CmnGeoEAPointType& pntNode);

    /**
     * This function will build the quaternion matrix
     */

    void buildQuatMatrix(const double Qw, const double Qi, const double Qj,
            const double Qk, double qMat2eci[VEC_SIZE][VEC_SIZE]);

    /**
     * This function will build the orbit frame matrix using the
     * ECR position and ECI velocity vectors.  The Z axis of the orbit
     * frame is determined using the ECR position vector and the geodetic
     * latitude.  The Z axis is then converted to ECI coordinates.
     */

    int buildECIOrbFrame(const CmnGeoEAPointType& inPt,
            double outOrb[VEC_SIZE][VEC_SIZE]);

    /**
     * This function will transform the ECI position, velocity and
     * quaternions to roll, pitch, and yaw in the ECI system (in radians).
     */

    int calcGDRollPitchYaw(CmnGeoEAPointType& pntNode);

    /**
     * This function will find the polar wander / UT1-UTC data
     * based on an input time.  This function may also update
     * the variable holding the time used to calculate Polar Wander
     * parameters.
     */

    int getPolarUT1UTC(const double inTAI, double& outTDTmUT1,
            double& outUT1mUTC, double& outXWander, double& outYWander,
            double& outTJD, bool bForceCalc = false) const;

    /**
     * This function will calculate the Apparent Sidereal Time, in hours,
     * for the input time.  Also returns Polar Wander values and the Solar
     * System Dynamic Julian Date for use by convCoordSys().
     */

    int appSiderealTime(const double inTAI, const bool getmean, double& gast,
            double& xWander, double& yWander, double& tjd);

    /**
     * This function computes earth radius when input latitude is geocentric.
     */

    double earth_radius_C(double rlat);

    /**
     * This function computes earth radius
     */

    double earth_radius_D(double rlat);

    /**
     * This function will perform the inverse of the NOVAS-C function
     * pnsw().  It will convert ECI coordinates to ECF coordinates.
     */

    void invPNSW(const double inTDT, const double gast, const double Xpole,
            const double Ypole, double eciVect[VEC_SIZE],
            double ecfVect[VEC_SIZE]);

    /**
     * This function will perform the inverse of the NOVAS-C function
     * wobble().
     */

    void invWobble(const double x, const double y, const double rotEq[VEC_SIZE],
            double geoEq[VEC_SIZE]);

    /**
     *  Compute solar vectors
     */

    int computeSolarVectors();

    /**
     *  Compute lunar vectors
     */

    int computeLunarVectors();

    /**
     *  Compute earth view vectors
     */

    int computeEarthViewVectors();

    /**
     * Checks lunar data for moon-in-SV conditions
     */

    int check_moon_in_sv();

    /**
     *  Validate lat and lon
     */

    inline bool validLat(const double inLat) {
        return ((inLat < -PIO2) || (inLat > PIO2) || std::isnan(inLat)) ?
                false : true;
    }

    inline bool validLon(const double inLon) {
        return ((inLon < (-PI)) || (inLon > PI) || std::isnan(inLon)) ? false : true;
    }

    /**
     * Convert signed 14 bit short to 16 bit
     */

    short convert14To16bit(short value14bit);

    /**
     * Determine if platform is little endian
     */

    bool isPlatformLittleEndian();

};

#endif /* VcstCmnGeo_H_ */
