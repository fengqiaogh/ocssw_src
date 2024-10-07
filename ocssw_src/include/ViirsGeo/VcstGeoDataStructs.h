/**************************************************************************
 *
 * NAME: VcstGeoDataStructs.h
 *
 * DESCRIPTION: contains structure definitions used in the ProSdrViirs 
 *  geolocation routines.  This file differs from GEO_parameters in that
 *  it defines input and output structures utilized by the I portion of the 
 *  IPO.
 *
 *
 **************************************************************************/

#ifndef VcstGeoDataStructs_H
#define VcstGeoDataStructs_H

#include <vector>
#include <VcstGeoParameters.h>
#include <VcstPolarStereographicDataSet.h>

/* Local definitions */

typedef struct ViirsGeoEncoderType {
    /* total number of good impulse data */
    int num_impulse_HAM[MAX_SCAN_NUMBER];
    int num_impulse_tel[MAX_SCAN_NUMBER];

    /* validation flags for telescope and mirror data (each scan) */
    int HAM_impulse_flag[MAX_SCAN_NUMBER];
    int tel_impulse_flag[MAX_SCAN_NUMBER];

    /* telescope and mirror impulse encoder for the scans in a granule */
    double HAM_impulse_enc[MAX_SCAN_NUMBER][MAX_IMPULSE_NUMBER_HAM];
    double tel_impulse_enc[MAX_SCAN_NUMBER][MAX_IMPULSE_NUMBER_TEL];

    /* telescope and mirror impulse time for scans in a granule */
    double HAM_impulse_time[MAX_SCAN_NUMBER][MAX_IMPULSE_NUMBER_HAM];
    double tel_impulse_time[MAX_SCAN_NUMBER][MAX_IMPULSE_NUMBER_TEL];

    /*DR4767 non-nominal ham encoder start*/
    signed char ham_start_not_nominal[VIIRS_SCANS];

} ViirsGeoEncoderType;


// Nadir values on a per scan basis.

typedef struct ViirsGeoNadirType {
    double lat[MAX_SCAN_NUMBER];
    double lon[MAX_SCAN_NUMBER];
} ViirsGeoNadirType;

typedef struct ViirsMoonArrays {
    double moon_vectors[VIIRS_SCANS][VEC_SIZE];
} ViirsMoonArrays;

typedef struct ViirsGeoProcType {
    double t_sync[MAX_SCAN_NUMBER];

    double sample_time[MAX_SCAN_NUMBER][MAX_SCAN_SAMPLE];

    double matTempCorr[3][3]; // Temperature correction matrix

    double DNB_aft_const;

    double u_aft[MAX_DETECTORS][3];

    double n_normal_side_1[3]; // Normal vector to HAM side 1
    double n_normal_side_2[3]; // Normal vector to HAM side 2

    //wwang J1 code change
    //Holds EV nadir frame #, #of agg zones, and # of pixels per agg zone
    //derived from geoParams LUT.
    int evNadirFrame;  //nadir frame number
    int numAggZones;   //number of agg zones
    int numPixelsPerAggZone[NUM_DNB_ZONES]; //number of samples per agg zone

    // Structure to hold Nadir latitude and longitude values for each scan.
    ViirsGeoNadirType nadir_values;

    // Contains decoded Telescope and Ham encoder times and values.
    ViirsGeoEncoderType telHamEncData;

    // Pointer to the Polarstereographic Data Set.
    VcstPolarStereographicDataSet* psds;
    VcstPolarStereographicDataSet* psds_i;
    VcstPolarStereographicDataSet* psds_m;
    VcstPolarStereographicDataSet* psds_d;

    // terrain_correction flag
    int TC;

    // State holder for the DEBUG message throttle.  Protect against
    // flooding the logs with DEBUG messages during a maneuver.
    int dbgMsgThrottle;

    /*store Scan Control Electronics side, DR4759*/
    signed char SCE_side[VIIRS_SCANS];

} ViirsGeoProcType;

typedef struct {
    unsigned char scan_mode[VIIRS_SCANS]; /* each scan day-night     */
    unsigned char mode; /* granule mode night-day-mixed         */
    unsigned char pad[3]; /* 3 Byte fill to fit to a 32 bit boundry */
} viirs_SDRhdr_type;


// A structure to contain the geolocation output parameters from
// the geolocatePixel() method.
// NOTE:  The GEO field height is not in the Geo Pixel structure because
// it is set from the Terrain Info in the storeGranule() method for IMG data
// and copied directly into the DMS buffer in calcModFromImg() for MOD data.

typedef struct ViirsGeoPixelType {
    double lat;
    double lon;
    double satazm;
    double satzen;
    double sunazm;
    double sunzen;
    double moonazm;
    double moonzen;
    double moonPhase;
    double moonFrac;
    double range;
    unsigned char qFlag;
} ViirsGeoPixelType;

// A structure to contain the Full Sized geolocation output arrays from
// the various VIIRS geolocation methods.  This structure is sized such
// that it can contain the largest output arrays and be reused within
// the VIIRS Geolocation routines.
// NOTE:  The GEO field height is not in the Full Geo structure because
// it is set from the Terrain Info in the storeGranule() method for IMG data
// and copied directly into the DMS buffer in calcModFromImg() for MOD data.

typedef struct ViirsGeoFullType {
    double lat[I_DETECTORS][I_VIIRS_COLS];
    double lon[I_DETECTORS][I_VIIRS_COLS];
    double satazm[I_DETECTORS][I_VIIRS_COLS];
    double satzen[I_DETECTORS][I_VIIRS_COLS];
    double sunazm[I_DETECTORS][I_VIIRS_COLS];
    double sunzen[I_DETECTORS][I_VIIRS_COLS];
    double moonazm[I_DETECTORS][I_VIIRS_COLS];
    double moonzen[I_DETECTORS][I_VIIRS_COLS];
    double range[I_DETECTORS][I_VIIRS_COLS];
    double grow[I_DETECTORS][I_VIIRS_COLS];
    double gcol[I_DETECTORS][I_VIIRS_COLS];
    unsigned char qFlags[I_DETECTORS][I_VIIRS_COLS];
    //Add pixel phase and fraction.
    double moonFrac[I_DETECTORS][I_VIIRS_COLS];
    double moonPhase[I_DETECTORS][I_VIIRS_COLS];
} ViirsGeoFullType;

typedef struct {
    long long scanStartTime[VIIRS_SCANS];
    long long scanMidTime[VIIRS_SCANS];
    float lat[M_DETECTORS][M_VIIRS_COLS];
    float lon[M_DETECTORS][M_VIIRS_COLS];
    float sunzen[M_DETECTORS][M_VIIRS_COLS];
    float sunazm[M_DETECTORS][M_VIIRS_COLS];
    float satzen[M_DETECTORS][M_VIIRS_COLS];
    float satazm[M_DETECTORS][M_VIIRS_COLS];
    float height[M_DETECTORS][M_VIIRS_COLS];
    float range[M_DETECTORS][M_VIIRS_COLS];
    float scPosition[VIIRS_SCANS][VEC_SIZE];
    float scVelocity[VIIRS_SCANS][VEC_SIZE];
    float scAttitude[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionStart[VIIRS_SCANS][VEC_SIZE];
    float scVelocityStart[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeStart[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionEnd[VIIRS_SCANS][VEC_SIZE];
    float scVelocityEnd[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeEnd[VIIRS_SCANS][VEC_RPY_SIZE];
    float scSunZen[VIIRS_SCANS];
    float scSunAzm[VIIRS_SCANS];
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
    unsigned char scanQuality[VIIRS_SCANS];
    unsigned char scanQuality2[VIIRS_SCANS];
    unsigned char pixelQuality[M_DETECTORS][M_VIIRS_COLS];
    unsigned char landwater[M_DETECTORS][M_VIIRS_COLS];
} viirs_SDR_MOD_Fgeoloc_type;

typedef struct {
    long long scanStartTime[VIIRS_SCANS];
    long long scanMidTime[VIIRS_SCANS];
    float lat[M_DETECTORS][M_A_VIIRS_COLS];
    float lon[M_DETECTORS][M_A_VIIRS_COLS];
    float sunzen[M_DETECTORS][M_A_VIIRS_COLS];
    float sunazm[M_DETECTORS][M_A_VIIRS_COLS];
    float satzen[M_DETECTORS][M_A_VIIRS_COLS];
    float satazm[M_DETECTORS][M_A_VIIRS_COLS];
    float height[M_DETECTORS][M_A_VIIRS_COLS];
    float range[M_DETECTORS][M_A_VIIRS_COLS];
    float scPosition[VIIRS_SCANS][VEC_SIZE];
    float scVelocity[VIIRS_SCANS][VEC_SIZE];
    float scAttitude[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionStart[VIIRS_SCANS][VEC_SIZE];
    float scVelocityStart[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeStart[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionEnd[VIIRS_SCANS][VEC_SIZE];
    float scVelocityEnd[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeEnd[VIIRS_SCANS][VEC_RPY_SIZE];
    float scSunZen[VIIRS_SCANS];
    float scSunAzm[VIIRS_SCANS];
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
    unsigned char scanQuality[VIIRS_SCANS];
    unsigned char scanQuality2[VIIRS_SCANS];
    unsigned char pixelQuality[M_DETECTORS][M_A_VIIRS_COLS];
    unsigned char landwater[M_DETECTORS][M_A_VIIRS_COLS];
} viirs_SDR_MOD_Unagg_Fgeoloc_type;

typedef struct {
    long long scanStartTime[VIIRS_SCANS];
    long long scanMidTime[VIIRS_SCANS];
    float lat[DNB_DETECTORS][DNB_VIIRS_COLS];
    float lat_tc[DNB_DETECTORS][DNB_VIIRS_COLS];
    float lon[DNB_DETECTORS][DNB_VIIRS_COLS];
    float lon_tc[DNB_DETECTORS][DNB_VIIRS_COLS];
    float sunzen[DNB_DETECTORS][DNB_VIIRS_COLS];
    float sunazm[DNB_DETECTORS][DNB_VIIRS_COLS];
    float satzen[DNB_DETECTORS][DNB_VIIRS_COLS];
    float satazm[DNB_DETECTORS][DNB_VIIRS_COLS];
    float moonzen[DNB_DETECTORS][DNB_VIIRS_COLS];
    float moonazm[DNB_DETECTORS][DNB_VIIRS_COLS];
    float height[DNB_DETECTORS][DNB_VIIRS_COLS];
    float height_tc[DNB_DETECTORS][DNB_VIIRS_COLS];
    float range[DNB_DETECTORS][DNB_VIIRS_COLS];
    float scPosition[VIIRS_SCANS][VEC_SIZE];
    float scVelocity[VIIRS_SCANS][VEC_SIZE];
    float scAttitude[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionStart[VIIRS_SCANS][VEC_SIZE];
    float scVelocityStart[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeStart[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionEnd[VIIRS_SCANS][VEC_SIZE];
    float scVelocityEnd[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeEnd[VIIRS_SCANS][VEC_RPY_SIZE];
    float scSunZen[VIIRS_SCANS];
    float scSunAzm[VIIRS_SCANS];
    //Add pixel phase and fraction.
    float moonPhase[DNB_DETECTORS][DNB_VIIRS_COLS];
    float moonFrac[DNB_DETECTORS][DNB_VIIRS_COLS];
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
    unsigned char scanQuality[VIIRS_SCANS];
    unsigned char scanQuality2[VIIRS_SCANS];
    unsigned char pixelQuality[DNB_DETECTORS][DNB_VIIRS_COLS];
    unsigned char pixelQuality_tc[DNB_DETECTORS][DNB_VIIRS_COLS];
    unsigned char landwater[DNB_DETECTORS][DNB_VIIRS_COLS];
} viirs_SDR_DNB_Fgeoloc_type;

typedef struct {
    long long scanStartTime[VIIRS_SCANS];
    long long scanMidTime[VIIRS_SCANS];
    float lat[I_DETECTORS][I_VIIRS_COLS];
    float lon[I_DETECTORS][I_VIIRS_COLS];
    float sunzen[I_DETECTORS][I_VIIRS_COLS];
    float sunazm[I_DETECTORS][I_VIIRS_COLS];
    float satzen[I_DETECTORS][I_VIIRS_COLS];
    float satazm[I_DETECTORS][I_VIIRS_COLS];
    float height[I_DETECTORS][I_VIIRS_COLS];
    float range[I_DETECTORS][I_VIIRS_COLS];
    float scPosition[VIIRS_SCANS][VEC_SIZE];
    float scVelocity[VIIRS_SCANS][VEC_SIZE];
    float scAttitude[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionStart[VIIRS_SCANS][VEC_SIZE];
    float scVelocityStart[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeStart[VIIRS_SCANS][VEC_RPY_SIZE];
    float scPositionEnd[VIIRS_SCANS][VEC_SIZE];
    float scVelocityEnd[VIIRS_SCANS][VEC_SIZE];
    float scAttitudeEnd[VIIRS_SCANS][VEC_RPY_SIZE];
    float scSunZen[VIIRS_SCANS];
    float scSunAzm[VIIRS_SCANS];
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
    unsigned char scanQuality[VIIRS_SCANS];
    unsigned char scanQuality2[VIIRS_SCANS];
    unsigned char pixelQuality[I_DETECTORS][I_VIIRS_COLS];
    unsigned char landwater[I_DETECTORS][I_VIIRS_COLS];
} viirs_SDR_IMG_Fgeoloc_type;

/***************************************************************************
NOTE: These structures provide a grid row and column for
every pixel in the granule.  The MDS used for the granule is also
copied here. These structures make granulating the Ancillary Data
a very fast process.  All granule grids are polar stereographic.
 **********************/


typedef struct {
    double grow[M_DETECTORS][M_VIIRS_COLS];
    double gcol[M_DETECTORS][M_VIIRS_COLS];
    double ctr_grow[VIIRS_SCANS][MOD_VIIRSI_NRCTNGL];
    double ctr_gcol[VIIRS_SCANS][MOD_VIIRSI_NRCTNGL];
    mds_type gmds;
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
} viirs_SDR_MOD_growcol_type;

typedef struct {
    double grow[DNB_DETECTORS][DNB_VIIRS_COLS];
    double gcol[DNB_DETECTORS][DNB_VIIRS_COLS];
    double ctr_grow[VIIRS_SCANS][DNB_VIIRSI_NRCTNGL];
    double ctr_gcol[VIIRS_SCANS][DNB_VIIRSI_NRCTNGL];
    mds_type gmds;
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
} viirs_SDR_DNB_growcol_type;

typedef struct {
    double grow[I_DETECTORS][I_VIIRS_COLS];
    double gcol[I_DETECTORS][I_VIIRS_COLS];
    double ctr_grow[VIIRS_SCANS][IMG_VIIRSI_NRCTNGL];
    double ctr_gcol[VIIRS_SCANS][IMG_VIIRSI_NRCTNGL];
    mds_type gmds;
    viirs_SDRhdr_type SDRhdr;
    int act_scans; /* SDR structure may be less than full */
} viirs_SDR_IMG_growcol_type;

// A structure to contain the pointers to the DMS output buffers.

typedef struct ViirsGeoOutputType {
    viirs_SDR_MOD_Fgeoloc_type* mod;
    viirs_SDR_MOD_Fgeoloc_type* mod_tc;
    //  viirs_SDR_MOD_Unagg_Fgeoloc_type* mod_unagg;
    viirs_SDR_IMG_Fgeoloc_type* img;
    viirs_SDR_IMG_Fgeoloc_type* img_tc;
    viirs_SDR_DNB_Fgeoloc_type* dnb;
    viirs_SDR_MOD_growcol_type* mod_grid;
    viirs_SDR_MOD_growcol_type* mod_tc_grid;
    viirs_SDR_IMG_growcol_type* img_grid;
    viirs_SDR_DNB_growcol_type* dnb_grid;
} ViirsGeoOutputType;

// Interpolation rectangle information.  Contains row and column information
// for the interpolation data.  Interpolation rectangle numbers for a scan
// are pushed onto the badRec vector if a point cannot be geolocated.

typedef struct {
    std::vector<int> badRec[VIIRS_SCANS];
    int numRctngl;
    int bgnRow;
    int midRow;
    int endRow;
    int lcol[VIIRSI_MAX_NRCTNGL];
    int mcol[VIIRSI_MAX_NRCTNGL];
    int rcol[VIIRSI_MAX_NRCTNGL];
} ViirsGeoRctnglType;


#endif

