/*******************************************************************************
 *
 * NAME: ViirsBandDnb.h
 *
 * DESCRIPTION: Subclass of ViirsBand object class to provide data structures and
 * processes specific to day-night band (DNB) calibration.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSBANDDNB_H_
#define VIIRSBANDDNB_H_

#include <VcstViirsBand.h>

struct dnbLutPtrs {
    // Pointers to LUT data structures
    proSdrViirsCalDnbFrameToZoneLUT* ftz;
    ProSdrViirsCalDnbLgsGainsLUT* TableDnbLgsGains;
    ProSdrViirsCalDnbGainRatiosLUT* TableDnbGainRatios;
    proSdrViirsDnbDn0Type* dnbDn0;
    proSdrViirsCaldnbRVSLUT* RVS_DNB;
    proSdrViirsCalDnbStrayLightLUT* strayLight;
};

class ViirsBandDnb : public ViirsBand {
public:

    /**
     *  Class constructor
     */

    ViirsBandDnb(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandDnb();

    /**
     *  Initialize all data required for calibration
     */

    int initialize();

    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data();

    /**
     *  Initialize LUT data
     */

    int initialize_LUT_data();

    /**
     *  Initialize L1B data
     */

    int initialize_L1B_data();

    /**
     *  Write L1B output data
     */

    int write_data(const NcFile* nc_output);

    // Get L1b data to pass to l2gen

    int get_l1bdata(float *l1bdata) {
        return 0;
    }

protected:

    // observation attributes:

    float obs_fill_value_;
    float obs_valid_min_;
    float obs_valid_max_;

    // quality attributes:

    unsigned short pixel_quality_fill_value_;
    vector<unsigned short> pixel_quality_flag_masks_;

    // scan quality flags

    static constexpr unsigned char L1B_SCAN_QUALITY_DNB_STRAY_LIGHT = 0x80; //x0000000

    // constants

    static constexpr int DETECTORS = 16;
    static constexpr int CAL_SAMPLES = 64;
    static constexpr int EV_PIXELS = 4064;
    static constexpr int C_COEFFS = 3;
    static constexpr int GAIN_STATES = 3;
    static constexpr int GAIN_RATIOS = 2;
    static constexpr int MIRROR = 2;
    static constexpr int AGG_ZONES = 32;
    static constexpr int SZA_BINS = 469;
    static constexpr int HEMISPHERES = 2;

    static constexpr int DNB_HIGH = 1;
    static constexpr int DNB_MED = 2;
    static constexpr int DNB_LOW = 3;

    static constexpr float radOffset = -0.0000000015;
    static constexpr float radRange = 0.04;

    static constexpr unsigned short SATURATED_DNB_DN = 8191;
    static constexpr unsigned short out_of_range_dn = 7860;
    //	float max_rad_;

    // DPP Configuration Data

    unsigned char dp_dnb_1a_1b_stage_[VIIRS_SCANS];
    unsigned char dp_dnb_dark_sub_eth_[VIIRS_SCANS];
    unsigned char dp_dnb_tmg_mode_[VIIRS_SCANS];

    // calibration and earth view data

    float DN_SV_[VIIRS_SCANS][DETECTORS][CAL_SAMPLES];
    float DN_BB_[VIIRS_SCANS][DETECTORS][CAL_SAMPLES];

    short DN_EV_[VIIRS_SCANS][DETECTORS][EV_PIXELS];
    unsigned char EV_GAIN_[VIIRS_SCANS][DETECTORS][EV_PIXELS];

    double CRVS_[DETECTORS][EV_PIXELS][GAIN_STATES][MIRROR][C_COEFFS];

    bool strayLightRegion_[VIIRS_SCANS];

    // output arrays

    float radiance_[VIIRS_SCANS * DETECTORS][EV_PIXELS];
    unsigned short pixelQuality_[VIIRS_SCANS * DETECTORS][EV_PIXELS];
    unsigned char rowQuality_[VIIRS_SCANS * DETECTORS];
    unsigned char badDetector_[DETECTORS];

    //  Pointer to structure of LUT pointer defined above

    dnbLutPtrs pLut_;

    //Map to store LUT input item objects.

    std::map<std::string, VcstLutInputItem*> dnbLutItems_;

    /**
     * compute_CRVS()
     *
     * Compute CRVS coefficients
     */

    int compute_CRVS();

    /**
     * calibrate_scan()
     *
     * Perform calibration for a single DNB scan
     */

    int calibrate_scan(int scan);

    /**
     * precheck_scan()
     *
     * Precheck_scan to determine if it is good for calibration.
     * Returns true or false.
     */

    bool precheck_scan(int scan);

    /**
     * check_pixel_quality()
     *
     * Check a pixel's quality and set quality flags as indicated.
     */
    using ViirsBand::check_pixel_quality;
    virtual int check_pixel_quality(int scan, int detector, int frame);

    /**
     * check_limits()
     *
     * Check whether a pixel is within its appropriate range, and set quality
     * flags if it is not.
     */

    int check_limits(int row, int frame);

};

#endif /* VIIRSBANDIMG_H_ */
