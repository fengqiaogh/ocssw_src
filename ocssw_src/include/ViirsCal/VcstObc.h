/******************************************************************************
 *  NAME: VcstObc.h
 *
 *  DESCRIPTION: Object class that generates an OBC file from a VIIRS L1A file.
 *
 *  Created on: January 2, 2015
 *      Author: Sam Anderson, VCST
 *
 ******************************************************************************/

#ifndef VcstObc_H_
#define VcstObc_H_

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <netinet/in.h>

#include <VcstCmnGeo.h>
#include <VcstViirsStructs.h>
#include <VcstCmnGeoStructs.h>
#include <VcstViirsThermal.h>
#include <netcdf>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


// Structure of LUT pointers used for temperatures

struct obcLutPtrs {
    SolarDiffRotationMatrixLutType* SdRotMatPtr_;
    SolarDiffVoltLutType* SdsmVoltPtr_;
    proSdrViirsCalRVSLUT* rvsLUT;
    proSdrViirsCalDgAnDnLmtLUT* dgAnDnLimits;
    proSdrViirsCalFPredictedTableLUT* TableFPredicted;
    ProSdrViirsCalSolarSpectralIradLUT* solarSpectralIradLUT;
};

class ViirsThermal;

class VcstObc {
public:

    int orbit_number_;
    int format_version_;
    int instrument_number_;
    string platform_;
    string startDirection_;
    string endDirection_;
    string day_night_flag_;
    string time_coverage_start_;
    string time_coverage_end_;
    string pge_start_time_;
    string pge_end_time_;
    string versionid_;
    string history_;
    string source_files_;

    // FIR filter configuration

    static constexpr int NUM_TAPS = 201;

    // constants

    static constexpr int RSB_BANDS = 14;
    static constexpr int TEB_BANDS = 7;
    static constexpr int MIRROR_SIDE = 2;

    static constexpr int SATURATED_DN = 4095;
    static constexpr int MAX_NUM_HISTORY_ORBITS = 480;

    // Static constants and arrays

    static const int processAtNight[Viirs_Bands];

    static constexpr double DEGtoRAD = 0.017453292519943296e0;
    static constexpr double RADtoDEG = 57.295779513082321e0;

    static constexpr int C_COEFFS = 4;
    static constexpr int A_COEFFS = 3;
    static constexpr int T_FP_LEVELS = 5;
    static constexpr int T_ELEC_LEVELS = 5;

    static constexpr unsigned char NIGHT_GRAN = 0;
    static constexpr unsigned char DAY_GRAN = 1;
    static constexpr unsigned char MIXED_GRAN = 2;

    // global constants
    static constexpr char SCE_A_SIDE_ON = 0;
    static constexpr char SCE_B_SIDE_ON = 1;
    static constexpr char SCE_INVALID = -1;

    static constexpr unsigned char L1B_SCAN_STATE_HA_MIRROR_B_SIDE = 0x01; //0000000x
    static constexpr unsigned char L1B_SCAN_STATE_ELECTRONICS_SIDE = 0x02; //000000x0
    static constexpr unsigned char L1B_SCAN_STATE_SENSOR_NIGHT_MODE = 0x04; //00000x00
    static constexpr unsigned char L1B_SCAN_STATE_SENSOR_NOT_OPERATIONAL = 0x08; //0000x000
    static constexpr unsigned char L1B_SCAN_STATE_NO_CALIBRATION = 0xff; //xxxxxxxx

    static constexpr unsigned char L1B_SCAN_QUALITY_MOON_IN_SPACE_VIEW = 0x01; //0000000x
    static constexpr unsigned char L1B_SCAN_QUALITY_NO_EV_DATA = 0x02; //000000x0
    static constexpr unsigned char L1B_SCAN_QUALITY_SENSOR_NOT_NOMINAL = 0x04; //00000x00
    static constexpr unsigned char L1B_SCAN_QUALITY_SCAN_SYNC_FAILURE = 0x08; //0000x000
    static constexpr unsigned char L1B_SCAN_QUALITY_TEL_START_NOT_NOMINAL = 0x10; //000x0000
    static constexpr unsigned char L1B_SCAN_QUALITY_BB_TEMP_NOT_NOMINAL = 0x20; //00x00000
    static constexpr unsigned char L1B_SCAN_QUALITY_LWIR_TEMP_NOT_NOMINAL = 0x40; //0x000000

    // numeric constants
    static constexpr long long MSECPERDAY = 86400000000ll; //24*60*60=86400 million
    static constexpr double TAI93_TAI58_SEC = 1104537627.0;
    static constexpr double RAD_TO_DEG = 180.0 / M_PI;
    static constexpr int RSR_SIZE = 1000;
    static constexpr int SOLAR_IRAD_SIZE = 49951;

    // General Pixel quality Flags
    static constexpr short NUM_PIXEL_QUALITY_FLAGS = 4;
    static constexpr unsigned char L1B_PIXEL_QUALITY_GOOD = 0x00; //00000000
    static constexpr unsigned char L1B_PIXEL_QUALITY_SUBSTITUTE_CAL = 0x01; //0000000x
    static constexpr unsigned char L1B_PIXEL_QUALITY_OUTOFRANGE_RAD = 0x02; //000000x0
    static constexpr unsigned char L1B_PIXEL_QUALITY_OUTOFRANGE_RFLBT = 0x02; //000000x0
    static constexpr unsigned char L1B_PIXEL_QUALITY_ALLSATURATION = 0x04; //00000x00
    static constexpr unsigned char L1B_PIXEL_QUALITY_TEMP_NOT_NOMINAL = 0x08; //0000x000
    static constexpr unsigned char L1B_PIXEL_QUALITY_LOWER_MASK = 0x0f; //0000xxxx
    static constexpr unsigned char L1B_PIXEL_QUALITY_UPPER_MASK = 0xf0; //xxxx0000
    static constexpr unsigned char L1B_PIXEL_QUALITY_NOT_ALLSATURATION = 0xfb; //xxxxx0xx
    static constexpr unsigned char L1B_PIXEL_QUALITY_NOT_OUTOFRANGE = 0xfd; //xxxxxx0x

    // DG Pixel quality Flags
    static constexpr short NUM_DG_PIXEL_QUALITY_FLAGS = 8;
    static constexpr short NUM_UADG_PIXEL_QUALITY_FLAGS = 6;
    static constexpr unsigned char L1B_PIXEL_QUALITY_LOW_GAIN_STATE = 0x10; //000x0000
    static constexpr unsigned char L1B_PIXEL_QUALITY_MIXED_GAIN_STATE = 0x20; //00x00000
    static constexpr unsigned char L1B_PIXEL_QUALITY_DG_ANOMALY = 0x40; //0x000000
    static constexpr unsigned char L1B_PIXEL_QUALITY_SOMESATURATION = 0x80; //00000x00

    // DNB Pixel quality Flags
    static constexpr short NUM_DNB_PIXEL_QUALITY_FLAGS = 5;
    static constexpr unsigned char L1B_PIXEL_QUALITY_DNB_STRAY_LIGHT = 0x10; //000x0000

    // Bad Pixel quality Flags
    static constexpr short NUM_BAD_PIXEL_QUALITY_FLAGS = 5;
    static constexpr short BAD_PIXEL_QUAL_SHIFT = 8;
    static constexpr unsigned short L1B_PIXEL_QUALITY_BOWTIE = 0x0100; //0000000x00000000
    static constexpr unsigned short L1B_PIXEL_QUALITY_MISSING = 0x0200; //000000x000000000
    static constexpr unsigned short L1B_PIXEL_QUALITY_NOCALIBRATE = 0x0400; //00000x0000000000
    static constexpr unsigned short L1B_PIXEL_QUALITY_DEAD_DETECTOR = 0x0800; //0000x00000000000
    static constexpr unsigned short L1B_PIXEL_QUALITY_NOISY_DETECTOR = 0x1000; //000x000000000000

    static constexpr short TOTAL_SG_PIXEL_QUALITY_FLAGS = NUM_PIXEL_QUALITY_FLAGS + NUM_BAD_PIXEL_QUALITY_FLAGS;
    static constexpr short TOTAL_DG_PIXEL_QUALITY_FLAGS = NUM_DG_PIXEL_QUALITY_FLAGS + NUM_BAD_PIXEL_QUALITY_FLAGS;
    static constexpr short TOTAL_UADG_PIXEL_QUALITY_FLAGS = NUM_UADG_PIXEL_QUALITY_FLAGS + NUM_BAD_PIXEL_QUALITY_FLAGS;
    static constexpr short TOTAL_DNB_PIXEL_QUALITY_FLAGS = NUM_DNB_PIXEL_QUALITY_FLAGS + NUM_BAD_PIXEL_QUALITY_FLAGS;

    // Row quality flags
    static constexpr unsigned char L1B_ROW_QUALITY_SUBSTITUTE_CAL = 0x01; //00000000
    static constexpr unsigned char L1B_ROW_QUALITY_NO_CALIBRATION = 0x02; //000000x0
    static constexpr unsigned char L1B_ROW_QUALITY_MOON_IN_SPACE_VIEW = 0x04; //000000x0
    static constexpr unsigned char L1B_ROW_QUALITY_DG_ANOMALY = 0x08; //0000x000
    static constexpr unsigned char L1B_ROW_QUALITY_DEAD_DETECTOR = 0x10; //000x0000
    static constexpr unsigned char L1B_ROW_QUALITY_BB_TEMP_NOT_NOMINAL = 0x20; //00x00000
    static constexpr unsigned char L1B_ROW_QUALITY_LWIR_TEMP_NOT_NOMINAL = 0x40; //0x000000
    static constexpr unsigned char L1B_ROW_QUALITY_DNB_STRAY_LIGHT = 0x80; //0x000000

    // Moon-in-SV processing parameters
    static constexpr int MOON_IN_SV_SCANS_TO_SKIP_FORWARD = 6;
    static constexpr int MOON_IN_SV_SCANS_TO_SKIP_BACK = 4;

    // Fill values
    static constexpr float L1B_FILL_FLOAT = -999.9f;
    static constexpr double L1B_FILL_DOUBLE = -999.9L;
    static constexpr short L1B_FILL_SHORT = -999;
    static constexpr unsigned short L1B_FILL_USHORT = 65535;
    static constexpr char L1B_FILL_BYTE = -1;
    static constexpr unsigned char L1B_FILL_UBYTE = 255;

    static constexpr float VDNE_FLOAT32_FILL = -999.9f;
    static constexpr float ELLIPSOID_FLOAT32_FILL = -999.9f;
    static constexpr float ERR_FLOAT32_FILL = -999.9f;
    static constexpr float ONBOARD_PT_FLOAT32_FILL = -999.9f;
    static constexpr float MISS_FLOAT32_FILL = -999.9f;
    static constexpr float NA_FLOAT32_FILL = -999.9f;

    // Fill value test thresholds
    static constexpr float FLOAT_FILL_TEST = -999.0f;
    static constexpr short SHORT_FILL_TEST = -990;
    static constexpr char BYTE_FILL_TEST = 0;
    static constexpr unsigned char UBYTE_FILL_TEST = 247;

    // Operational constants
    static constexpr short GRANULES = 3;
    static const float var_limit_[Viirs_Bands - 1];
    static const float asp_factor_[Viirs_Bands - 1];
    static constexpr float sv_filter_bw_ = 0.01;
    static constexpr float f_filter_bw_ = 0.005;

    // Trim Table and ranges
    static const int TrimTable_I_[Iband_detectors][2];
    static const int TrimTable_M_[Mband_detectors][2];

    static const float radOffset_[Viirs_Bands - 1];
    static const float radRange_[Viirs_Bands - 1];
    static const float reflOffset_[Viirs_Bands - 1];
    static const float reflRange_[Viirs_Bands - 1];
    static const unsigned short destRange_;
    static const unsigned short out_of_range_dn_[Viirs_Bands - 1];
    static const unsigned int destRange_M13_;

    // Radiance/reflectance scale factors
    float  refl_scale_factor_[Viirs_Bands - 1];
    float  refl_add_offset_[Viirs_Bands - 1];
    float  rad_scale_factor_[Viirs_Bands - 1];
    float  rad_add_offset_[Viirs_Bands - 1];

    // special processing
    bool processModByScan_;
    int scanToProcess_;
    int extract_pixel_start_;
    int extract_pixel_stop_;

    static const string band_prefix_[Viirs_Bands];

    GRAN_SEQ_ENUM granSeq_;
    bool bValidSeq_;

    /**
     *  Adjacent granules, before and after
     */

    VcstObc* obcSeq_[GRANULES];

    // Pointer to instance of ViirsCmnGeo object
    VcstCmnGeo* pCmnGeo_;

    // Pointer to instance of ViirsThermal object
    ViirsThermal* pTherm_;

    // Pointer to structure of LUT pointer defined above
    obcLutPtrs pLut_;

    // Granule data
    int act_scans_;

   // Engineering Status Data
    unsigned char mode_;
//    double startTAI93sec_[VIIRS_SCANS]; /* TAI scan start time */
//    double endTAI93sec_[VIIRS_SCANS]; /* TAI scan end time */
//    double sv_midTAI93sec_[VIIRS_SCANS]; /* TAI scan ev mid time */
//    double bb_midTAI93sec_[VIIRS_SCANS]; /* TAI scan ev mid time */
//    double sd_midTAI93sec_[VIIRS_SCANS]; /* TAI scan ev mid time */
//    int scan_number_[VIIRS_SCANS];
//    unsigned char ham_side_[VIIRS_SCANS];
//    unsigned char sensor_mode_[VIIRS_SCANS];
    unsigned char cal_metadata_[VIIRS_SCANS][Cal_Metadata];

    // Derived status

//    int electronics_side_[VIIRS_SCANS];
//    bool scan_sync_failure_[VIIRS_SCANS];
//    char tel_start_not_nominal_[VIIRS_SCANS];
//    int scan_sync_failure_cnt_;
    bool scan_no_calibration_[VIIRS_SCANS];
    unsigned char scan_state_[VIIRS_SCANS];
    unsigned char scan_quality_[VIIRS_SCANS];

    // asp offset data

    short ASP_I_[Number_of_Ibands][Number_of_Scans][Mirror_Sides][Iband_detectors];
    short ASP_M_[Number_of_Mbands + 1][Number_of_Scans][Mirror_Sides][Mband_detectors];
    short ASP_DIFF_I_[Number_of_Ibands][Number_of_Scans][Mirror_Sides][Iband_detectors];
    short ASP_DIFF_M_[Number_of_Mbands][Number_of_Scans][Mirror_Sides][Mband_detectors];

    // Calibration data

    float SV_DN_M_[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    float SV_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    float SV_DN_DNB_[Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    float SD_DN_M_[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    float SD_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    float SD_DN_DNB_[Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    float BB_DN_M_[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    float BB_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    float BB_DN_DNB_[Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    unsigned char GAIN_STATE_[Number_of_DG_bands][Number_of_Scans][Mband_detectors];
    unsigned char DNB_sequence_[Number_of_Scans];
    unsigned char DNB_sequence_valid_min_, DNB_sequence_valid_max_;

    unsigned char SDSM_active_[Number_of_Scans];
    unsigned char SDSM_position_[Number_of_Scans];
    float SDSM_sample_[Number_of_Scans][SDSM_DETECTORS][SDSM_SAMPLES];

    // Processed digital counts data

    float AVG_SV_DN_M_[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors];
    float AVG_SV_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Parity];
    float STD_SV_DN_M_[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors];
    float STD_SV_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Parity];
    float AVG_BB_DN_M_[Number_of_Mbands][Number_of_Scans][Mband_detectors];
    float AVG_BB_DN_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Parity];

    // Validated calibration data

    float SV_M_[Number_of_Mbands][Number_of_Scans][Mband_detectors][Gain_States];
    float SV_I_[Number_of_Ibands][Number_of_Scans][Iband_detectors][Agg_Zones][Parity];
    float BB_M_[Number_of_MTbands][Number_of_Scans][Mband_detectors][Gain_States];
    float BB_I_[Number_of_ITbands][Number_of_Scans][Iband_detectors][Agg_Zones][Parity];

    // Granule averages

    float AVG_SV_M_[Number_of_Mbands][Mirror_Sides][Mband_detectors][Gain_States];
    float AVG_SV_I_[Number_of_Ibands][Mirror_Sides][Iband_detectors][Agg_Zones][Parity];
    float AVG_BB_M_[Number_of_MTbands][Mirror_Sides][Mband_detectors][Gain_States];
    float AVG_BB_I_[Number_of_ITbands][Mirror_Sides][Iband_detectors][Agg_Zones][Parity];

    // Solar calibration

    double CP_M_[Number_of_MRbands][Number_of_Scans][Mband_detectors][Gain_States][C_Coefs];
    double CP_I_[Number_of_IRbands][Number_of_Scans][Iband_detectors][C_Coefs];

    // Quality

    unsigned char rowQuality_M_[Number_of_Mbands][Number_of_Scans*Mband_detectors];
    unsigned char rowQuality_I_[Number_of_Ibands][Number_of_Scans*Iband_detectors];

    // Diagnostics

    bool bLunar_;
    bool bFilters_;
    double fir_coeff_[NUM_TAPS];

    float diag_sv_I_[5][4][GRANULES * VIIRS_SCANS][32][3][2];
    float diag_sv_M_[16][4][GRANULES * VIIRS_SCANS][16][2];
    float diag_f_I_[2][2][GRANULES * VIIRS_SCANS][32][3][2];
    float diag_f_M_[5][2][GRANULES * VIIRS_SCANS][16][2];
                              /**
      *  Class constructor
      */

    VcstObc();

    /**
     *  Class destructor
     */

    ~VcstObc();

    /**
     *  Initialize all data
     */

    int initialize(GRAN_SEQ_ENUM gran_seq);

    /**
     *  Initialize calibration data for a single band
     */

    int initialize(VIIRS_BAND_ENUM band);

    /**
     * process_bands()
     *
     * Compute calibration data for all bands
     */

    int process_bands();

    /**
     * Write OBC file
     */

    int write_obc_data();

    void setHistory(std::string history) {
        history_ = history;
    }

    std::string getHistory() {
        return history_;
    }

    void setSource(std::string source) {
        source_files_ = source;
    }

    std::string getSource() {
        return source_files_;
    }

    inline const bool doFilters() {
        return bFilters_;
    };

    inline const void setFilters(const bool bFilters) {
        bFilters_ = bFilters;
    };

private:

    //Map to store LUT input item objects.

    std::map<std::string, VcstLutInputItem*> obcLutItems_;

    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data(GRAN_SEQ_ENUM sequence);

    /**
     *  Initialize LUT data
     */

    int initialize_LUT_data();

    /**
     *  Initialize L1B data
     */

    int initialize_L1B_data();

    /**
     * read_L1A_attributes()
     *
     * Read L1A attributes
     */

    int read_L1A_attributes(NcFile* nc_input);

    /**
     * read_L1A_cal_data()
     *
     * Read L1A calibration data for all bands
     */

    int read_L1A_cal_data_img(NcFile* nc_input);
    int read_L1A_cal_data_mod(NcFile* nc_input);
    int read_L1A_cal_data_dnb(NcFile* nc_input);

    /**
     * check_scan_state()
     *
     * Check all scans state and set state flags
     */

    int check_scan_state();

    /**
     * check_scan_quality()
     *
     * Check all scans quality and set quality flags
     */

    int check_scan_quality();

    /**
     * compute_refl_scale_factor()
     *
     * Compute reflectance scale factors
     */

    int compute_refl_scale_factor(VIIRS_BAND_ENUM band);

    /**
     * compute_cal_stats()
     *
     * Compute calibration statistics
     */

    int compute_cal_stats(VIIRS_BAND_ENUM band);

    /**
     * swap_even_odd
     *
     * Swap even and odd calibration frames
     */

    int swap_even_odd(VIIRS_BAND_ENUM band);

    /**
     * validate_cal()
     *
     * Validate SV and BB
     */

    int validate_cal(VIIRS_BAND_ENUM band);

    /**
     * compute_Cprime()
     *
     * Apply F factors to calibration coefficients
     */

    int compute_Cprime(VIIRS_BAND_ENUM band);

    /**
     * NAME: filter_SV()
     *
     * Filter SV data
     */

    int filter_sv(VIIRS_BAND_ENUM band);

    /**
     * NAME: filter_F()
     *
     * Filter TEB F factor data
     */

    int filter_F(VIIRS_BAND_ENUM band);

    /**
     * write_global_attributes()
     *
     * Write global attributes to file
     */

    int write_global_attributes(NcFile* nc_output);

    /**
     * write_obc_eng_data()
     *
     * Write OBC engineering data to file
     */

    int write_obc_eng_data(NcFile* nc_output);

    /**
     * write_obc_nav_data()
     *
     * Write OBC navigation data to file
     */

    int write_obc_nav_data(NcFile* nc_output);

    /**
     * write_obc_cal_data()
     *
     * Write OBC calibration data to file
     */

    int write_obc_cal_data(NcFile* nc_output);

    /**
     * write_diagnostics()
     *
     * Write OBC diagnostics data to file
     */

    int write_diagnostics(NcFile* nc_output);

    /**
     * NAME: generate_filt_coeff( )
     *
     * DESCRIPTION: Compute FIR_low_pass filter taps.
     */

    int generate_filt_coeff(float cutoff_frequency);

    /**
     *  isValidGranuleSequence()
     *
     * Determine whether a granule sequence is valid.  Return true
     * if valid, otherwise return false.
     */

    bool isValidGranuleSequence();

    /**
     * convertSigned14BitTo16bit()
     *
     * Convert signed 14 bit short to 16 bit
     */

    short convert14To16bit(short value14bit);

    /**
     * isPlatformLittleEndian()
     *
     * Determine if platform is little endian
     */

    bool isPlatformLittleEndian();
};

#endif /* VcstObc_H_ */
