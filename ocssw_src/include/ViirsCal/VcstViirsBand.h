/*******************************************************************************
 *
 * NAME: ViirsBand.h
 *
 * DESCRIPTION: Base class for all ViirsBand object classes. Support generic data
 * and functions that pertain to all VIIRS bands.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSBAND_H_
#define VIIRSBAND_H_

#include <VcstCmnConsts.h>
#include <VcstObc.h>

class ViirsBand {
public:

    /**
      *  Class constructor
      */

     ViirsBand(VcstObc* pObc, VIIRS_BAND_ENUM band);

     /**
      *  Class destructor
      */

     virtual ~ViirsBand();

     /**
      *  Initialize the ViirsBand class
      */

     virtual int initialize();

     /**
      *  Calibrate all scans in a granule.
      */

     int calibrate();

     /**
      * Precheck_scan to determine if it is good for calibration.
      * Returns true or false.
      */

     virtual bool precheck_scan(int scan);

     /**
      *  Calibrate a single scan.
      */

     virtual int calibrate_scan(int scan);

     /**
      *  Get L1b data pointer to pass to l2gen
      */

     virtual int get_l1bdata(float *l1bdata);

    /**
     *  Initialize L1A data
     */

    virtual int initialize_L1A_data();

     /**
      *  Write L1B output data
      */

     virtual int write_data(const NcFile* nc_output);
     virtual int write_rsb(const NcFile* nc_output, float *obs_val,
             unsigned short *short_val, unsigned short *pixelQuality,
             int8_t *unc_idx, bool bRad);
     virtual int write_teb(const NcFile* nc_output, float *obs_val,
             unsigned short *short_val, unsigned int *long_val,
             unsigned short *pixelQuality, int8_t *unc_idx, bool bRad);
     virtual int write_cdg(const NcFile* nc_output);

     // Configuration data
     // Static constants and arrays
     static const int processAtNight[Viirs_Bands];

    // FIR filter configuration
    static constexpr int NUM_TAPS = 201;

    // constants
    static constexpr int RSB_BANDS = 14;
    static constexpr int TEB_BANDS = 7;
    static constexpr int MIRROR_SIDE = 2;

    static constexpr int SATURATED_DN = 4095;
    static constexpr int MAX_NUM_HISTORY_ORBITS = 480;

    static constexpr double DEGtoRAD = 0.017453292519943296e0;
    static constexpr double RADtoDEG = 57.295779513082321e0;

    static constexpr int C_COEFFS = 4;
    static constexpr int A_COEFFS = 3;
    static constexpr int T_FP_LEVELS = 5;
    static constexpr int T_ELEC_LEVELS = 5;

    static constexpr unsigned char NIGHT_GRAN = 0;
    static constexpr unsigned char DAY_GRAN = 1;
    static constexpr unsigned char MIXED_GRAN = 2;

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
    static constexpr unsigned char L1B_PIXEL_QUALITY_SOMESATURATION = 0x80; //x0000000

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
    static constexpr unsigned int L1B_FILL_UINT_M13 = 327680;
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

    string band_str_;
    string ev_var_str_;

    // float attributes:
    float obs_fill_value_;
    float obs_valid_min_;
    float obs_valid_max_;

    // Conversion factors
    float  flt_valid_min_;
    float  flt_valid_max_;
    float  flt_scale_factor_;
    float  flt_add_offset_;
    float  shrt_rad_scale_factor_;
    float  shrt_rad_add_offset_;

    // short attributes:
    unsigned short shrt_fill_value_;
    unsigned short shrt_valid_min_;
    unsigned short shrt_valid_max_;
    unsigned int M13_fill_value_;
    unsigned int M13_valid_min_;
    unsigned int M13_valid_max_;
    float shrt_scale_factor_;
    float shrt_add_offset_;

    // brightness temperature attributes
    float bt_fill_value_;
    float bt_valid_min_;
    float bt_valid_max_;

    // quality attributes:
    unsigned short pixel_quality_fill_value_;
    vector<unsigned short> pixel_quality_flag_masks_;

    //detector quality flag attributes:
    vector<unsigned char> detector_quality_flag_;
    unsigned char detector_quality_flag_masks_[8];

    // uncertainty index attributes
    float uncidx_scale_factor_;
    unsigned short uncidx_fill_value_;
    unsigned short uncidx_valid_min_;
    unsigned short uncidx_valid_max_;
    vector<unsigned int> uncertainty_index_flag_values_;
    vector<unsigned int> m13_uncertainty_index_flag_values_;

    // pointers used for writing output
    unsigned short *shrtPtr_;
    unsigned int *longPtr_;
    unsigned short *qualPtr_;
    int8_t  *uncPtr_;
    float *radPtr_;
    float *cdgPtr_;

protected:

    bool bLWIR_;

    // OBC pointer

    VcstObc* pObc_;

    //  Band identifiers

    VIIRS_BAND_ENUM band_; // Viirs band type
    VIIRS_TEB_BAND_ENUM teb_band_; // Viirs TEB band type
    VIIRS_RSB_BAND_ENUM rsb_band_; // Viirs TEB band type

    short det_432_; // Viirs-wide detector index

    /**
     *  Initialize L1B output data
     */

    virtual int initialize_L1B_data();

    /**
     * check_pixel_quality()
     *
     * Check EV and set quality flags
     */

    virtual int check_pixel_quality(int band, int scan, int det,
            short& DN_EV, unsigned short& qual);

};

#endif /* VIIRSBAND_H_ */
