/*******************************************************************************
 *
 * NAME: ViirsBandImg.h
 *
 * DESCRIPTION: Subclass of ViirsBand object class that provides data structures
 * and processes common to all imagery (IMG) bands required to support calibration
 * processing. It is the base class for the the ViirsBandImgRsb and ViirsBandImgTeb
 * object classes.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSBANDIMG_H_
#define VIIRSBANDIMG_H_

#include <VcstViirsCal.h>
#include <VcstViirsBand.h>

class ViirsBandImg : public ViirsBand {

public:
    /**
     *  Class constructor
     */

    ViirsBandImg(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandImg();

    /**
     *  Initialize class
     */

    virtual int initialize();

    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data();

    // Get L1b data to pass to l2gen

    virtual int get_l1bdata(float *l1bdata) {
        return 0;
    }

    // constants
    static const int I_BANDS = 5;
    static const int DETECTORS = 32;
    static const int C_COEFFS = 4;
    static const int CAL_SAMPLES = 96;
    static const int EV_PIXELS = 6400;
    static const int MIRROR = 2;
    static const int AGG_ZONES = 3;
    static const int PARITY = 2;

    static const int MIN_1_1_AGG_AREA_375M = 1280;
    static const int MAX_1_1_AGG_AREA_375M = 5119;
    static const int MIN_1_2_1_AGG_AREA_375M = 1280;
    static const int MAX_1_2_1_AGG_AREA_375M = 2016;
    static const int MIN_1_2_2_AGG_AREA_375M = 4384;
    static const int MAX_1_2_2_AGG_AREA_375M = 5120;
    static const int MIN_1_3_AGG_AREA_375M = 2016;
    static const int MAX_1_3_AGG_AREA_375M = 4384;

    short DN_EV_[VIIRS_SCANS][DETECTORS][EV_PIXELS];

    // Output arrays
    float radiance_[VIIRS_SCANS*DETECTORS][EV_PIXELS];
    unsigned short short_val_[VIIRS_SCANS*DETECTORS][EV_PIXELS];
    unsigned short pixelQuality_[VIIRS_SCANS*DETECTORS][EV_PIXELS];
    int8_t unc_idx_[VIIRS_SCANS*DETECTORS][EV_PIXELS];

    unsigned char badDetector_[DETECTORS];

    // Band identifiers
    VIIRS_I_BAND_ENUM img_band_; // Viirs IMG band index
    VIIRS_IMG_TEB_BAND_ENUM iteb_band_; // Viirs IMG TEB band type

protected:

    /**
     *  Initialize L1B data
     */

    virtual int initialize_L1B_data();

   /**
     * precheck_scan
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

    int check_pixel_quality(int band, int scan, int det,
            short& DN_EV, unsigned short& qual);

    /**
     * calibrate_scan()
     *
     * Calibrate TEB scan
     */

    int calibrate_scan(int scan);

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single RSB pixel and set quality flags as indicated.
     */

    virtual int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& refl_bt, unsigned short& shrt,
            unsigned short& qual);


};


class ViirsBandImgRsb : public ViirsBandImg {
public:

    /**
     *  Class constructor
     */

    ViirsBandImgRsb(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandImgRsb();

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single RSB pixel and set quality flags as indicated.
     */

    int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& refl, unsigned short& shrt,
            unsigned short& qual);

};

class ViirsBandImgTeb : public ViirsBandImg {
public:
    /**
     *  Class constructor
     */

    ViirsBandImgTeb(VcstObc* pObc, VIIRS_BAND_ENUM band);

    /**
     *  Class destructor
     */

    ~ViirsBandImgTeb();

     /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single TEB pixel and set quality flags as indicated.
     */

    int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, float& bt, unsigned short& shrt,
            unsigned short& qual);

};

#endif /* VIIRSBANDIMG_H_ */
