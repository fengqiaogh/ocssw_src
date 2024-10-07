/*******************************************************************************
 *
 * NAME: ViirsThermal.h
 *
 * DESCRIPTION: Object class that provides data structures and processes specific
 * to the processing and management of temperature data required for calibration.
 * The ViirsGranule object class maintains a pointer to a ViirsThermal object
 * class, through which it is able to gain access to all available temperature
 * information.
 *
 * Created on: Aug 25, 2014
 *     Author: Sam Anderson, VCST
 *
 *******************************************************************************/

#ifndef VIIRSTHERMAL_H_
#define VIIRSTHERMAL_H_

#include <VcstObc.h>
#include <VcstCmnGeo.h>
#include "VcstCalLutInputItem.h"

// Structure of LUT pointers used for temperatures

struct tempLutPtrs {
    proSdrViirsCalEmissiveLUT* EmissiveLUT;
    proSdrViirsCalRMParametersLUT* rmParameters;
    proSdrViirsCalHAMERLUT* HAMERLUT;
    proSdrViirsCalRTAERLUT* RTAERLUT;
    proSdrViirsCalOBCERLUT* OBCERLUT;
    proSdrViirsCalOBCRRLUT* OBCRRLUT;
    proSdrViirsCalLtoEBBTLUT* EBBTLUT;
    proSdrViirsCalTeleCoeffLUT* tCoeffs;
    proSdrViirsCalBBTempCoeffs* bBTempCoeffs;
    proSdrViirsCalDetectorResponseLUT* TableA;
    proSdrViirsCalInstrumentResponseLUT* TableB;
    proSdrViirsCalDeltaCTempLUT* DeltaCTable;
    proSdrViirsCalGainTableLUT* GainTable;
};

// Define static arrays

class ViirsGranule;
class VcstObc;

class ViirsThermal {
public:

    static constexpr int FOCAL_PLANES = 3;
    static constexpr int REF_RESISTANCES = 3;
    static constexpr int TEMP_COEFFS = 6;
    static constexpr int TEMP_MEASUREMENTS = 2;
    static constexpr int BB_THERMISTORS = 6;
    static constexpr int NUM_THERM = 20;

    static constexpr int TEMP_DN_LOW = -9000;
    static constexpr int THERM_LOW_VALUE = -8192;
    static constexpr int THERM_HI_VALUE = 8192;

    static constexpr int TEB_LUT_SIZE = 65536;
    static constexpr int TEB_LUT_VALUES = 65528;
    static constexpr int TEB_LUT_SIZE_M13 = 327681;
    static constexpr int TEB_LUT_VALUES_M13 = 327671;

    static constexpr double TOLLERANCE = 1.0E-20; // defined zero
    static constexpr double LOWBTEMP = 0.0;

    static constexpr int C_COEFFS = 4;
    static constexpr int A_COEFFS = 3;
    static constexpr int T_FP_LEVELS = 5;
    static constexpr int T_ELEC_LEVELS = 5;
    static constexpr int DGSTATES = 2;

    // DCR Coefficients

    static constexpr float C0_dcr = -5.381061f;
    static constexpr float C1_dcr = 1.357921e-003f;
    static constexpr float C2_dcr = 8.193453e-009f;
    static constexpr float C3_dcr = -1.679609e-012f;

    static constexpr int DCR_LIMIT = 8192;

    // Focal plane assembly types

    static constexpr unsigned char VISNIR = 0;
    static constexpr unsigned char SMWIR = 1;
    static constexpr unsigned char LWIR = 2;
    static constexpr unsigned char BG_DNB = 3;
    static const unsigned char band_to_fp_[Viirs_Bands];

    /**
     *  Class constructor
     */

    ViirsThermal();
    ViirsThermal(VcstObc* pObc);

    /**
     *  Class destructor
     */

    ~ViirsThermal();

    /**
     *  Initialize class data
     */

    int initialize(GRAN_SEQ_ENUM gran_seq);

    /**
     * interpolate_Temp_to_Rad
     *
     * @param bandIndex     band of interest for the interpolation
     * @param *array_tp     array of temperatures
     * @param *array_rad    radiances corresponding to temperatures in array_tp
     * @param indices[][2]  number of bands for which there is data, offset
     * @param tp            temperature for which radiance is derived
     *
     * @return rad_interp derived radiance
     */

    double interpolate_Temp_to_Rad(int bandIndex,
            double *array_tp,
            double *array_rad,
            int indices[][2],
            double tp);

    /**
     * interpolate_L_to_EBBT
     *
     * @param bandIndex     band of interest for the interpolation
     * @param *array_tp     array of temperatures
     * @param *array_rad    radiances corresponding to temps in array_tp
     * @param indices[][2]  number of bands with data
     * @param rad
     *
     * @return tp_interp  derived temperature
     */

    double interpolate_L_to_EBBT(int bandIndex,
            double *array_tp,
            double *array_rad,
            int indices[][2],
            double rad);

    /**
     * compute_C
     *
     * Generate calibration coefficients
     */

    int compute_C(VIIRS_BAND_ENUM band);

    /**
     * compute_F()
     *
     * Compute F factors for TEB calibration
     */

    int compute_F(VIIRS_BAND_ENUM band);

    /**
     *  Temperature data in the form used by the calibration algorithms
     */

    float Telec_[VIIRS_SCANS]; // Calculate C and C'
    float fp_[FOCAL_PLANES][VIIRS_SCANS]; // Calculate C and C' (Tdet)
    float Tcav_[VIIRS_SCANS]; // Calibrate TEB
    float Trta_[VIIRS_SCANS]; // Calibrate TEB
    float Ttele_[VIIRS_SCANS]; // Calibrate TEB
    float Tsh_[VIIRS_SCANS]; // Calibrate TEB
    float mir_avg_[VIIRS_SCANS]; // Calibrate TEB (Tham)
    float bb_avg_[VIIRS_SCANS]; // Calibrate TEB (Tobc)
    float Tomm_[VIIRS_SCANS];
    float RrtaRhamDiff_[Number_of_Tbands][VIIRS_SCANS];

    bool bbTempNotNominal_[VIIRS_SCANS]; //DR4710 bb not nominal (warm up/cool down)
    bool lwirFPAtempNonNominal_[VIIRS_SCANS]; //DR4501 lwir fpa test
    bool badThermValue_[VIIRS_SCANS];
    bool badInterpolateValue_[VIIRS_SCANS];

    float bt_lut_[7][TEB_LUT_SIZE];
    float bt_lut_M13_[TEB_LUT_SIZE_M13];

    double C_I_[Number_of_Ibands][VIIRS_SCANS][Iband_detectors][C_Coefs];
    double C_M_[Number_of_Mbands][VIIRS_SCANS][Mband_detectors][Gain_States][C_Coefs];

    float F_I_[Number_of_ITbands][VIIRS_SCANS][Iband_detectors][Agg_Zones][Parity];
    float F_M_[Number_of_MTbands][VIIRS_SCANS][Mband_detectors][Gain_States];

    float AVG_F_I_[Number_of_ITbands][Mirror_Sides][Iband_detectors][Agg_Zones][Parity];
    float AVG_F_M_[Number_of_MTbands][Mirror_Sides][Mband_detectors][Gain_States];

    //  Pointer to VcstObc object class

    VcstObc* pObc_;

    //  Pointer to structure of LUT pointer defined above

    tempLutPtrs pLut_;

    // L1A raw data
    short bbTemp_[VIIRS_SCANS][BB_THERMISTORS];
    short mf_tel_blkhd_py_[VIIRS_SCANS];
    short mf_scan_cavity_nxp_[VIIRS_SCANS];
    short mf_scan_cavity_baf_pz_[VIIRS_SCANS];
    short mf_scan_cavity_baf_nz_[VIIRS_SCANS];
    short ap_lw_cca_[VIIRS_SCANS];
    short ap_vn_cca_[VIIRS_SCANS];
    short dp_dnb_cca_[VIIRS_SCANS];
    short ap_sm_cca_[VIIRS_SCANS];
    short mf_ao_blkhd_px_nz_[VIIRS_SCANS];
    short mf_ao_blkhd_nx_pz_[VIIRS_SCANS];
    short mf_stopassy_baff_nz_[VIIRS_SCANS];
    short mf_fold_mirror_blkhd_[VIIRS_SCANS];
    short mf_ham_blkhd_[VIIRS_SCANS];
    short ft_lw_cfpa_hi_rsl_[VIIRS_SCANS];
    short ft_lw_cfpa_lo_rsl_[VIIRS_SCANS];
    short ft_sm_cfpa_hi_rsl_[VIIRS_SCANS];
    short ft_sm_cfpa_lo_rsl_[VIIRS_SCANS];
    short ft_vis_nir_fpa_[VIIRS_SCANS];
    short ham_tmp1_[VIIRS_SCANS];
    short ham_tmp2_[VIIRS_SCANS];
    short ct_prec_tref_mux1ca1_[VIIRS_SCANS];
    short ct_prec_tref_mux1ca2_[VIIRS_SCANS];
    short ct_prec_tref_mux1ca3_[VIIRS_SCANS];
    short bb_htr_temp_[VIIRS_SCANS];


protected:

    /**
     *  Initialize L1A data
     */

    int initialize_L1A_data(GRAN_SEQ_ENUM gran_seq);

    /**
     *  Initialize LUT data
     */

    int initialize_LUT_data();

    // Temperature Data

    unsigned char ft_lw_80k_[VIIRS_SCANS];
    unsigned char ft_lw_82k_[VIIRS_SCANS];

private:

    // Actual number of scans in granule
    int act_scans_;

    //Map to store LUT input item objects.

    std::map<std::string, VcstLutInputItem*> tempLutItems_;

    static const float bbSetTempCoefs_[6];

    struct thermistor_struct {
        short buffer[VIIRS_SCANS];
        short max_valid_value;
    };

    thermistor_struct thermistors_[NUM_THERM];

    int quality_flag_[Viirs_Bands][VIIRS_SCANS];

    /**
     * convert_temp_poly()
     *
     * Converts digital number reading of thermistor to temperature in Kelvins
     */

    float convert_temp_poly(short tdn, float c[6], float toffset);

    /**
     * coeffMult()
     *
     * Compute value of polynomial and return as float
     */

    float coeffMult(const float coeffs[], int numCoef, float floatDN);

    /**
     * black_body()
     *
     * Compute black body temperatures
     */

    int black_body();

    /**
     * focal_plane()
     *
     * Compute focal plane temperature
     */

    int focal_plane();

    /**
     * opto_mechanical()
     *
     * Compute TOMM weighted average of opto-mechanical module thermistors
     */

    int opto_mechanical();

    /**
     * scan_mirror()
     *
     * Compute scan mirror temperature
     */

    int scan_mirror();

    /**
     * shield()
     *
     * Compute shield temperature
     */

    int shield();

    /**
     * rta()
     *
     * Compute rotating telescope assembly (RTA) temperature
     */

    int rta();

    /**
     * telescope()
     *
     * Compute telescope temperature
     */

    int telescope();

    /**
     * cavity()
     *
     * Compute cavity temperature
     */

    int cavity();

    /**
     * electronics()
     *
     * Compute electronics temperature
     */

    int electronics();

    /**
     * generate_rad_bt_lut()
     *
     * Generate 100,000 element LUTs from EBBT LUT
     */

    int generate_rad_bt_lut();

    /**
     * convertSigned14BitTo16bit()
     *
     * Convert signed 14 bit short to 16 bit, and return the converted short.
     */

    short convert14To16bit(short value14bit);


};

#endif /* VIIRSTEMPERATURES_H_ */
