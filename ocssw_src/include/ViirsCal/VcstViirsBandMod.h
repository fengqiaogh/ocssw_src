/*******************************************************************************
 *
 * NAME: ViirsBandMod.h
 *
 * DESCRIPTION: Subclass of ViirsBand object class that provides data structures
 * and processes common to all moderate resolution (MOD) bands required to support
 * calibration processing. It is the base class for the the ViirsBandModDG and
 * ViirsBandModSG object classes.
 *
 *  Created on: Aug 25, 2014
 *      Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSBANDMOD_H_
#define VIIRSBANDMOD_H_

#include <VcstViirsBand.h>

class ViirsBandMod : public ViirsBand {

public:
    // constants

    static const int M_BANDS = 16;

    static const int EV_PIXELS = 3200;
    static const int DETECTORS = 16;
    static const int CAL_SAMPLES = 48;
    static const int MIRROR = 2;
    static const int GAIN_STATE = 2;
    static const int C_COEFFS = 4;

    // dual-gain sample attributes

    vector<unsigned short> sample_quality_flag_masks_;

    /**
     *  Band identifiers
     */
    VIIRS_M_BAND_ENUM mod_band_; // Viirs MOD band index
    VIIRS_MOD_SG_BAND_ENUM sg_band_; // Viirs SG band type
    VIIRS_MOD_DG_BAND_ENUM dg_band_; // Viirs SG band type
    VIIRS_MOD_RSB_BAND_ENUM mrsb_band_; // Viirs M TEB band type
    VIIRS_MOD_TEB_BAND_ENUM mteb_band_; // Viirs M TEB band type

    /**
     *  Output arrays
     */

    typedef short dn_ev_array[DETECTORS][EV_PIXELS];
    typedef float flt_pix_array[EV_PIXELS];
    typedef unsigned short shrt_pix_array[EV_PIXELS];
    typedef unsigned int int_pix_array[EV_PIXELS];
    typedef int8_t byte_pix_array[EV_PIXELS];
    dn_ev_array *DN_EV_;
    flt_pix_array *radiance_;
    shrt_pix_array *short_val_;
    int_pix_array *long_val_;
    shrt_pix_array *pixelQuality_;
    byte_pix_array *unc_idx_;

    unsigned char badDetector_[DETECTORS];

    /**
     *  Class constructor
     */

    ViirsBandMod(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandMod();

    /**
     *  Initialize class
     */

    virtual int initialize();

    /**
     *  Get a pointer to L1B data
     */

    int get_l1bdata(float *l1bdata);


    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data();

protected:

    /**
     *  Initialize L1B data
     */

    virtual int initialize_L1B_data();

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
     * Check EV and set quality flags
     */

    virtual int check_pixel_quality(int band, int scan, int det,
            short& DN_EV, unsigned short& qual);

    /**
     * calibrate_scan()
     *
     * Calibrate TEB scan
     */

    virtual int calibrate_scan(int scan);

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single RSB pixel and set quality flags as indicated.
     */

    virtual int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& refl_bt, unsigned int& shrt,
            unsigned short& qual);

    /**
     * NAME: getSVFromEV()
     *
     * DESCRIPTION: Obtain background reference from EV rather than SV.
     */

    virtual int getSVFromEV();
};


class ViirsBandModSGRsb : public ViirsBandMod {

public:

    /**
     *  Class constructor
     */

    ViirsBandModSGRsb(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandModSGRsb();

protected:

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single RSB pixel and set quality flags as indicated.
     */

    int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& refl, unsigned int& shrt,
            unsigned short& qual);

};


class ViirsBandModSGTeb : public ViirsBandMod {

public:

     /**
     *  Class constructor
     */

    ViirsBandModSGTeb(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandModSGTeb();


protected:

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single RSB pixel and set quality flags as indicated.
     */

    int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& refl, unsigned int& shrt,
            unsigned short& qual);

};

#endif /* VIIRSBANDMOD_H_ */
