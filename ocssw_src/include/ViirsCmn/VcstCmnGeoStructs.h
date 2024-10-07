/*****************************************************************************
 *
 * NAME:  VcstCmnGeoStructs.h
 *
 * DESCRIPTION: Defines the data structures used in VcstCmnGeo.
 *
 * Adapted directly from ProSdrCmnGeoStruct.h developed by Raytheon Company.
 *
 *****************************************************************************/

#ifndef VcstCmnGeoStructs_h
#define VcstCmnGeoStructs_h

#include <map>
#include <string>
#include <vector>
#include <VcstCmnConsts.h>
#include <VcstCalLutStructures.h>

// Type definition for the Cmn Geo Anomaly Stage.

typedef enum {
    GEO_ANOMALY_STAGE_0 = 0,
    GEO_ANOMALY_STAGE_1 = 1,
    GEO_ANOMALY_STAGE_2 = 2,
    GEO_ANOMALY_STAGE_3 = 3
} CmnGeoAnomalyStage;

// Type definition for the Cmn Geo E&A Check Error.

typedef enum {
    GEO_EA_CHK_NO_ERROR = 0, // No Errors
    GEO_EA_CHK_ANG_MOMENTUM = 1, // Angular Momentum
    GEO_EA_CHK_ANG_MOMENTUM_Z = 2, // Angular Momentum Z Component
    GEO_EA_CHK_ORIBIT_CONSISTENCY = 4, // Orbit Consitency
    GEO_EA_CHK_POS_ABS = 8, // Position Absolute
    GEO_EA_CHK_POS_MAG = 16, // Position Magnitude
    GEO_EA_CHK_VEL_ABS = 32, // Velocity Absolute
    GEO_EA_CHK_VEL_MAG = 64, // Velocity Magnitude
    GEO_EA_CHK_ATTITUDE = 128, // Attitude RPY
} CmnGeoEACheckError;

/*****************************************************************
 * Definitions for the VIIRS resolution types
 ****************************************************************/

enum BandType {
    BAND_TYPE_MOD = 0, BAND_TYPE_IMG, BAND_TYPE_DNB, BAND_TYPE_MOD_UNAGG
};

// Contains E&A information for a specific time

struct CmnGeoEAPointType {
    double tai;
    long long iet;
    double cLat;
    double dLat;
    double lon;
    double alt;
    double direc;
    double direcR;
    double roll;
    double pitch;
    double yaw;
    double iPos[VEC_SIZE];
    double iVel[VEC_SIZE];
    double iAcc[VEC_SIZE];
    double rPos[VEC_SIZE];
    double rVel[VEC_SIZE];
    double rAcc[VEC_SIZE];
    double quat[QUAT_SIZE];
    double matQ[VEC_SIZE][VEC_SIZE];
    CmnGeoAnomalyStage eaStage;
    unsigned short eaError;
    std::string maneuverType;

    CmnGeoEAPointType() {
        InitPoint();
    }
    ;

    void InitPoint(void) {
        tai = 0.0;
        iet = 0;
        cLat = 0.0;
        dLat = 0.0;
        lon = 0.0;
        alt = 0.0;
        direc = 0.0;
        direcR = 0.0;
        roll = 0.0;
        pitch = 0.0;
        yaw = 0.0;
        for (int i = 0; i < VEC_SIZE; i++) {
            iPos[i] = 0.0;
            iVel[i] = 0.0;
            iAcc[i] = 0.0;
            rPos[i] = 0.0;
            rVel[i] = 0.0;
            rAcc[i] = 0.0;
            matQ[i][0] = 0.0;
            matQ[i][1] = 0.0;
            matQ[i][2] = 0.0;
        }
        for (int j = 0; j < QUAT_SIZE; j++) {
            quat[j] = 0.0;
        }
        eaStage = GEO_ANOMALY_STAGE_0;
        eaError = 0;
        maneuverType = "(empty)";
    }
    ;
};

// Type definition for the list of E&A points.
// The key for this map is the TAI time of the point.
typedef std::map<double, CmnGeoEAPointType> CmnGeoEphAttPointsMap;

// Contains E&A information for the entire granule

struct CmnGeoEphAttType {
    CmnGeoEphAttPointsMap points;
    std::string satId;
    std::string version;
    long long ietCreated;
    long long ietBgnTime;
    long long ietEndTime;
    double taiBgnTime;
    double taiEndTime;
    double deltaTAI;

};

// Contains subset of ADCS_HSK (APID 8) data

struct AdcsHskTlmType {
    double taiTime;
    unsigned char adManDone;
    unsigned short adFftId;
};

// Contains subset of Bus Critical (APID 0) data

struct BusCritTlmType {
    double taiTime;
    unsigned char adState;
};

// Contains S/C Position/Velocity data

struct CmnGeoEphemType {
    double taiTime;
    float pos[VEC_SIZE];
    float vel[VEC_SIZE];
};

// Contains S/C Quaternion data

struct CmnGeoAttType {
    double taiTime;
    float quat[QUAT_SIZE];
};

// Contains S/C Ephem/Attitude byte aligned data

struct CmnGeoScEAType {
    std::vector<CmnGeoEphemType> ephem;
    std::vector<CmnGeoAttType> attitude;
    std::vector<AdcsHskTlmType> adcsHskTlm;
    std::vector<BusCritTlmType> busCritTlm;
    double startTime; // Start time of granule in TAI
    double stopTime; // Stop time of granule in TAI
};

/*******************************************************************
 *
 * orbital_elements - structure which contains the orbital elements
 *    tle_string1 & 2      - Two-Line-Element strings
 *    sat_ID               - the satellite identifier string
 *    epochtime            - epoch time of the element set
 *    epochTAI             - epoch time in TAI
 *    epoch_rev            - rev number of ascending node
 *    satrec               - satellite record used by SGP4 routines.
 *
 *  These elements have been removed from the structure since they
 *  are also in the "satrec" structure used by SGP4.  The comments
 *  are left here so that the names used by SGP4 can be more easily
 *  translated and understood.
 *    firtimedermeanmot    - first time derivative of mean motion
 *    sectimedermeanmot    - second time derivative of mean motion
 *    sectimedermeanmotexp - exponent of sec time derivative of mean mot
 *    bstardragterm        - BSTAR drag term
 *    bstardragtermexp     - exponent of BSTAR drag term
 *    inclination          - inclination
 *    rightascenascnode    - right acension of the ascending node
 *    eccentricity         - eccentricity
 *    argperigee           - argument of perigee
 *    meananomaly          - mean anomaly
 *    meanmotion           - mean motion
 *
 *******************************************************************/

// Provide a forward reference to elsetrec to avoid conflicts between
// SGP4's definition of PI and ones used elsewhere in PRO.
struct elsetrec;
// elsetrec is defined in sgp4unit.h

struct OrbitalElementsType {
    char tle_string1[80];
    char tle_string2[80];

    char sat_ID[32]; // NORAD names are up to 24 characters
    double epochtime; // TLE year-days.fractionOfDay
    double epochTAI; // TLE epoch time in TAI seconds
    int epoch_rev; // revolution number at epoch (not stored
    //   in the SGP4 satrec)

    double TDTmUT1; // diff between Terrestrial Dynamic Time and
    //    UT1 (units of seconds)
    double UT1mUTC; // diff between UT1 - UTC
    double xWander; // X polar wander
    double yWander; // Y polar wander

    double UT1SatEpoch; // TLE epoch time in UT1.  Note satrec
    //    stores the epoch in UTC.
    elsetrec *satrec; // pointer to a record of SGP4 element sets

};

//  Variables to hold the values of the polar wander data within
//  an algorithm run, will be reset between taskings.

struct polarWanderDataType {
    double TDTmUT1; // diff Terrestrial Dynamic - UT1, seconds
    double UT1mUTC; // difference between UT1 and UTC, seconds
    double prevTAI; // TAI on a previous call
    double xWander; // x coordinate of polar wander, arc seconds
    double yWander; // y coordinate of polar wander, arc seconds

};

// USNO Polar Wander / UT1-UTC data

struct CmnGeoPolarUT1Type {
    double txyTAI[13000];
    double xWander[13000];
    double yWander[13000];
    double UT1mUTC[13000];
    double dateIERSBullA; // TAI format date
    double xwCon[5];
    double ywCon[5];
    double xyConA[2];
    double xyConC[2];
    double bessel[3];
    double conUT2mUT1[4];
    double conUT1mUTC[3];
    int numPts;
    // XLC's "power" alignment implicitly puts a pad in this location.
    unsigned char pad[4];
};

struct InstrumentLutType {
    int tel_start_enc_nominal;
    int ham_start_enc_nominal[2];
    int scan_encdr_start_max;
};

struct CmnGeoParamLutType {
    double AngularMomentumLimit[CMNGEO_MIN_MAX_DIM];
    double AngularMomentumZLimit[CMNGEO_MIN_MAX_DIM];
    double OrbitConsistency;
    double PositionAbsLimit[CMNGEO_MIN_MAX_DIM];
    double PositionMagLimit[CMNGEO_MIN_MAX_DIM];
    double VelocityAbsLimit[CMNGEO_MIN_MAX_DIM];
    double VelocityMagLimit[CMNGEO_MIN_MAX_DIM];
    double AttitudeAbsLimit[CMNGEO_MIN_MAX_DIM];
};

struct JPLEphemLutType {
    double tjdStart;
    double tjdEnd;
    double span;
    double au;
    double emrat;
    int ipt[JPL_IPT_ROW][JPL_IPT_COL];
    char pad[4];
    double coeffs[JPL_EPHEM_ROW][JPL_EPHEM_COL];
};

struct SAALutType {
    double centerLat;
    double centerLon;
    double maxIndex;
    double latHeight;
    double lonWidth;
};

struct SolarDiffVoltLutType {
    float voltCoef[SDSM_SAMPLES][SDSM_DETECTORS][SDSM_COEF];
    int voltLowerLimit[SDSM_SAMPLES][SDSM_DETECTORS];
    int voltUpperLimit[SDSM_SAMPLES][SDSM_DETECTORS];
};

struct SolarDiffRotationMatrixLutType {
    float sdMatrix[SD_MATRIX_ROW][SD_MATRIX_COL];
    float sdsmMatrix[SD_MATRIX_ROW][SD_MATRIX_COL];
};

struct proSdrViirsCalQALUT {
    unsigned char Detector_Quality_Flag_Values[NUM_DETECTORS][8];
    float moon_offset_limits[NUM_BANDS][NUM_MOON_OFFSET_LIMITS];
    float saa_threshold;
};

struct L1ANavDataType {
    double adcs_time[SC_DIARY_RECORDS_10HZ];
    short adfftid[SC_DIARY_RECORDS_10HZ];
    unsigned char admandone[SC_DIARY_RECORDS_10HZ];
    double bus_time[SC_DIARY_RECORDS];
    unsigned char adstate[SC_DIARY_RECORDS];
    double att_time[SC_DIARY_RECORDS_10HZ];
    float att_quat[SC_DIARY_RECORDS_10HZ][QUAT_SIZE];
    double orb_time[SC_DIARY_RECORDS_10HZ];
    float orb_pos[SC_DIARY_RECORDS_10HZ][VEC_SIZE];
    float orb_vel[SC_DIARY_RECORDS_10HZ][VEC_SIZE];
};

#endif
