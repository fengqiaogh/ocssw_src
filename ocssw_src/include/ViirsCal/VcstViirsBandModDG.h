/*******************************************************************************
 *
 * NAME: ViirsBandModDG.h
 *
 * DESCRIPTION: Subclass of ViirsBandMod object class that provides data structures
 * and processes common to all dual-gain moderate resolution bands required to
 * support calibration processing. It is the base class for the the
 * ViirsBandModDGRsb and ViirsBandModDGTeb object classes.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSBANDMODDG_H_
#define VIIRSBANDMODDG_H_

#include <VcstViirsBandMod.h>

class ViirsBandModDG : public ViirsBandMod {
protected:

    static const int EV_SAMPLES = 6304;

    // Aggregation zones for moderate resolution dual-gain bands

    static const int FIRST_AG_START = 0; //  FIRST AGG ZONE USE 1 Pixel
    static const int FIRST_AG_END = 639; //  END OF FIRST ZONE
    static const int SECOND_AG_START = 640; //  START OF SECOND AGG ZONE 2 Pixels Avg
    static const int SECOND_AG_END = 1375; //  END OF SECOND AGG ZONE
    static const int FOURTH_AG_START = 4928; //  START OF FOURTH AGG ZONE 2 PIXELS Avg
    static const int FOURTH_AG_END = 5663; //  END OF FOURTH AGG ZONE
    static const int FIFTH_AG_START = 5664; //  START OF LAST AGG ZONE
    static const int FIFTH_AG_END = 6303; //  END OF LAST AGG ZONE

    static const int HIGH_GAIN = 0;
    static const int LOW_GAIN = 1;

    // Trim Table for unaggregated moderate resolution dual gain bands

    static const int TrimTable_DG_[DETECTORS][2];

    typedef short dn_ev_array[DETECTORS][EV_SAMPLES];
    typedef unsigned char ev_gain_array[DETECTORS][EV_SAMPLES];
    typedef float ev_rad_array[EV_SAMPLES];
    typedef unsigned short ev_flag_array[EV_SAMPLES];

    dn_ev_array *DN_EV_;
    ev_gain_array *EV_GAIN_;
    ev_rad_array *EV_RAD_;
    ev_flag_array *EV_FLAG_;

    bool bCdg_;

public:

    /**
     *  Class constructor
     */

    ViirsBandModDG(VcstObc* pObc, VIIRS_BAND_ENUM band, bool bCdg);

    /**
     *  Class destructor
     */

    ~ViirsBandModDG();

    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data();


protected:

    /**
     * calibrate_scan
     *
     * Calibrate dual-gain scan
     */

    int calibrate_scan(int scan);

    /**
     * NAME: calibrate_sample()
     *
     * DESCRIPTION: Calibrate single DG frame and set quality flags as indicated.
     */

    virtual int calibrate_sample(int row, int scan, int det, int& frm, short MS,
         short DN_EV, float& rad, unsigned short& qual, short& gain);

    /**
     * check_pixel_quality()
     *
     * Check EV and set quality flags
     */
    using ViirsBand::check_pixel_quality;
    virtual int check_pixel_quality(int band, int scan, int det,
            short& DN_EV, unsigned short& qual, short gain);

    /**
     * NAME: calibrate_pixel()
     *
     * DESCRIPTION: Calibrate single aggregated DG pixel and set quality
     * flags as indicated.
     */
    using ViirsBandMod::calibrate_pixel;
    virtual int calibrate_pixel(int row, int scan, int det, int& frm, short MS,
            short* DN_EV, float& rad, float& refl,
            unsigned int& shrt, unsigned short& qual);

    /**
     *  Write calibrated dual-gain output data as radiance with float data type
     */

    int write_cdg( const NcFile* nc_output );

    /**
     * NAME: getSVFromEV()
     *
     * DESCRIPTION: Obtain background reference from EV rather than SV.
     */

    int getSVFromEV();
};


class ViirsBandModDGRsb : public ViirsBandModDG {
public:

    /**
     *  Class constructor
     */

    ViirsBandModDGRsb(VcstObc* pObc, VIIRS_BAND_ENUM band, bool bCdg);

    /**
     *  Class destructor
     */

    ~ViirsBandModDGRsb();

protected:

    /**
     * NAME: calibrate_sample()
     *
     * DESCRIPTION: Calibrate single DG frame and set quality flags as indicated.
     */

    int calibrate_sample(int row, int scan, int det, int& frm, short MS,
            short DN_EV, float& rad, unsigned short& qual, short& gain);

};


class ViirsBandModDGTeb : public ViirsBandModDG {
public:

    /**
     *  Class constructor
     */

    ViirsBandModDGTeb(VcstObc* pObc, VIIRS_BAND_ENUM band, bool bCdg);

    /**
     *  Class destructor
     */

    ~ViirsBandModDGTeb();

 protected:

    /**
     * NAME: calibrate_sample()
     *
     * DESCRIPTION: Calibrate single TEB sample and set quality flags as indicated.
     */

    int calibrate_sample(int row, int scan, int det, int& frm,
        short MS, short DN_EV, float& rad, unsigned short& qual, short& gain);

};

#endif /* VIIRSBANDMODDG_H_ */
