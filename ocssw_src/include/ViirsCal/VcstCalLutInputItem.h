/*******************************************************************************
 *
 * NAME: VcstCalLutInputItem.h
 *
 * DESCRIPTION: Subclass of the the VcstLutInputItem object class to provide input
 * file processing specific to calibration LUTs--mainly endian conversion.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VcstCalLutInputItem_h
#define VcstCalLutInputItem_h

#include <VcstLutInputItem.h>

class VcstCalLutInputItem : public VcstLutInputItem {
public:

    /**
     * Constructor
     *
     */

    VcstCalLutInputItem(const std::string& groupName, size_t size);

    /**
     * Destructor
     *
     */

    virtual ~VcstCalLutInputItem();

protected:

    /**
     * Convert endianness of a calibration LUT file.
     *
     */

    int convertEndianness();

    int convert_DG_ANOMALY_DN_LIMITS_LUT();
    int convert_DNB_STRAY_LIGHT_CORRECTION_LUT();
    int convert_SDR_DNB_DN0_LUT();
    int convert_SDR_DNB_RVS();
    int convert_SDR_DNB_FRAME_TO_ZONE();
    int convert_SDR_DNB_LGS_GAINS_LUT();
    int convert_SDR_DNB_GAIN_RATIOS_LUT();
    int convert_SDR_F_PREDICTED_LUT();
    int convert_SDR_GAIN_LUT();
    int convert_SDR_HAM_ER_TABLE();
    int convert_SDR_RTA_ER_TABLE();
    int convert_SDR_OBC_ER_TABLE();
    int convert_SDR_OBC_RR_TABLE();
    int convert_SDR_EBBT_TABLE();
    int convert_SDR_TELE_COEFFS();
    int convert_SDR_SOLAR_IRAD_LUT();
    int convert_SDR_RSR_LUT();
    int convert_SDR_OBS_TO_PIXELS();
    int convert_SDR_RADIOMETRIC_PARAMETERS();
    int convert_SDR_EMISSIVE_LUT();
    int convert_SDR_REFLECTIVE_LUT();
    int convert_SDR_RVS_LUT();
    int convert_SDR_BB_TEMP_COEFFS();
    int convert_SDR_DELTA_C_LUT();
    int convert_SDR_COEFF_A_LUT();
    int convert_SDR_COEFF_B_LUT();
    int convert_SDR_RELATIVE_SPECTRAL_RESPONSE_LUT();

};



#endif
