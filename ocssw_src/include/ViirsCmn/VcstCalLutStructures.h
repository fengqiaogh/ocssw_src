/*******************************************************************************
 *
 * NAME: ViirsLutStructures.h
 *
 * DESCRIPTION: Header file with definitions of constants and data structures for
 * all calibration LUTs .
 *
 * Created on: Feb 13, 2015
 *     Author: Sam Anderson, VCST
 *
 *********************************************************************************/

#ifndef VIIRSLUTSTRUCTURES_H_
#define VIIRSLUTSTRUCTURES_H_

const int NUM_ELECTRONICS_SIDE = 2;
const int NUM_VIIRS_BAND = 22;
const int MAX_NUM_DETECTOR = 32;
const int NUM_M_I_GAIN_STATES = 2;
const int NUM_MIRROR_SIDES = 2;
const int NUM_TDET_LEVEL = 5;
const int DIM_NUM4 = 4;
const int NUM_TELE_LEVEL = 5;
const int MAX_NUM_C_COEF = 4;
const int NUM_FOCAL_PLANES = 3;
const int telecThermistorNum = 25;
const int NUM_REFL_750M_DG_BANDS = 6;
const int NUM_DETECTORS_750M = 16;
const int PRO_VIIRS_MIN_MAX_DIM = 2;
const int NUM_AGG_SEQ = 36;
const int MAX_NUM_GAIN = 3;
const int NUM_REFL_PLUS_DNB_BANDS = 15;
const int NUM_TCAV = 2;
const int NUM_TELEC_THERM = 4;
const int NUM_TEMP_COEFFS = 6;
const int MAX_NUM_VIIRS_SCANS = 48;
const int NUM_OF_EBBT_FILES = 7;
const int NUM_OF_HAMER_FILES = 7;
const int NUM_OF_OBCER_FILES = 7;
const int NUM_OF_OBCRR_FILES = 7;
const int NUM_OF_RTAER_FILES = 7;
const int MAX_HAM_ER_INDEX = 7498;
const int NUM_TFPLW = 2;
const int NUM_TFPSM = 2;
const int NUM_THERMISTORS = 26;
const int NUM_TMIR = 2;
const int NUM_TOMM_THERM = 5;
const int NUM_TRTA = 2;
const int NUM_TSH = 2;
const int NUM_TTELE = 2;
const int NUM_TTEL_TSH_TCAV = 3;
const int NUM_DETECTORS_DNB = 16;
const int MAX_RSR_VALUES = 1000;
const int OneTempTwoTolerances = 3;
const int MAX_RTA_ER_INDEX = 7498;
const int NUM_375M_BANDS = 5;
const int EV_375M_FRAMES = 6400;
const int NUM_DETECTORS_375M = 32;
const int NUM_750M_SG_BANDS = 9;
const int EV_750M_SG_FRAMES = 3200;
const int NUM_MAX_THERMISTORS = 26;
const int EV_DNB_FRAMES = 4064;
const int NUM_SZA_BINS = 469;
const int NUM_750M_DG_BANDS = 7;
const int X_Y_SIZE = 49951;
const int MAX_OBC_ER_INDEX = 7498;
const int NUM_MAX_EV_CT_PREC_TREF_MUX1CA = 3;
const int EV_FRAMES_750m_DG = 6304;
const int EV_750M_DG_FRAMES = 6304;
const int NUM_HEMISPHERES = 2;
const int MAX_EBBT_INDEX = 1050000;
const int NUM_BB_THERMISTORS = 6;
const int MAX_OBC_RR_INDEX = 7498;
const int NUM_T_MIR_THERMISTORS = 2;
const int DNB_VIIRS_RDR_COLS = 4064;
const int N_DNB_BANDS = 4;
const int NUM_BANDS = 22;
const int NUM_DETECTORS = 432;
const int NUM_MOON_OFFSET_LIMITS = 4;
const int NUM_ZONES = 32;
const int NUM_DNB_GAIN_RATIO_COEFFS = 3;
const int NUM_GAIN_RATIOS = 2;
const int SHORT_DNB_FRAMES = 127;

const int RSR_SIZE = 1000;
const int SOLAR_IRAD_SIZE = 49951;

struct proSdrViirsCalDetectorResponseLUT {
    double data[NUM_ELECTRONICS_SIDE][NUM_VIIRS_BAND][MAX_NUM_DETECTOR][NUM_M_I_GAIN_STATES][NUM_MIRROR_SIDES][NUM_TDET_LEVEL][DIM_NUM4];
};

struct proSdrViirsCalInstrumentResponseLUT {
    double data[NUM_ELECTRONICS_SIDE][NUM_VIIRS_BAND][MAX_NUM_DETECTOR][NUM_MIRROR_SIDES][NUM_TELE_LEVEL][DIM_NUM4];
};

// Delta C Temperature LUT

struct proSdrViirsCalDeltaCTempLUT {
    double DeltaC[NUM_ELECTRONICS_SIDE][NUM_VIIRS_BAND][MAX_NUM_DETECTOR][NUM_M_I_GAIN_STATES][NUM_MIRROR_SIDES][MAX_NUM_C_COEF][telecThermistorNum];
    double Tele[NUM_ELECTRONICS_SIDE][NUM_TELE_LEVEL];
    double Tdet[NUM_ELECTRONICS_SIDE][NUM_FOCAL_PLANES][NUM_TDET_LEVEL];
};

struct proSdrViirsCalDgAnDnLmtLUT {
    float data[NUM_REFL_750M_DG_BANDS][NUM_DETECTORS_750M][PRO_VIIRS_MIN_MAX_DIM];
};

// C Coefficients LUT

struct proSdrViirsCalDnbCCoeffLUT {
    double DNBCoeffs[NUM_ELECTRONICS_SIDE][NUM_ZONES][NUM_DETECTORS_DNB][MAX_NUM_GAIN][3];
};

struct proSdrViirsDnbDn0Type {
    float data[EV_DNB_FRAMES][NUM_DETECTORS_DNB][MAX_NUM_GAIN][NUM_MIRROR_SIDES];
};

// Frame to Zone Table

struct proSdrViirsCalDnbFrameToZoneLUT {
    int FrameToZone[EV_DNB_FRAMES];
};

struct ProSdrViirsCalDnbGainRatiosLUT {
    double DnbGainRatios[NUM_AGG_SEQ][NUM_DETECTORS_DNB][NUM_GAIN_RATIOS][NUM_DNB_GAIN_RATIO_COEFFS];
};

struct ProSdrViirsCalDnbLgsGainsLUT {
    signed char Fit_type;
    unsigned char implicit_pad0[7];
    long long T_ref;
    double F_ref[NUM_AGG_SEQ][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
    double F_param_1[NUM_AGG_SEQ][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
    double F_param_2[NUM_AGG_SEQ][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
};

struct proSdrViirsCaldnbRVSLUT {
    float dnbrvs[EV_DNB_FRAMES][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
};

struct proSdrViirsCalDnbStrayLightLUT {
    float MAX_RADIANCE_STRAY;
    float VIIRS_STRAY_SZA_GRID[NUM_HEMISPHERES][NUM_SZA_BINS];
    float VIIRS_DNB_SDR_STRAY_OFFSET[NUM_HEMISPHERES][NUM_SZA_BINS][EV_DNB_FRAMES][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
};

struct shortViirsCalDnbStrayLightLUT {
    float MAX_RADIANCE_STRAY;
    float VIIRS_STRAY_SZA_GRID[NUM_HEMISPHERES][NUM_SZA_BINS];
    float VIIRS_DNB_SDR_STRAY_OFFSET[NUM_HEMISPHERES][NUM_SZA_BINS][SHORT_DNB_FRAMES][NUM_DETECTORS_DNB][NUM_MIRROR_SIDES];
};

struct proSdrViirsCalLtoEBBTLUT {
    int EBBT_indices[NUM_OF_EBBT_FILES][2];
    double L_to_EBBT_tp[MAX_EBBT_INDEX];
    double L_to_EBBT_rad[MAX_EBBT_INDEX];
};

struct proSdrViirsCalEmissiveLUT {
    short wordBoundaryPad;
    short SV_DN_first_frame_to_use;
    short SV_DN_number_of_frames_to_use;
    short SV_DN_moon_include_frames;
    short BB_DN_first_frame_to_use;
    short BB_DN_number_of_frames_to_use;
    int T_mir_function_flag[NUM_T_MIR_THERMISTORS];
    float t_mir_default;
    float BB_Weight[NUM_BB_THERMISTORS];
};

struct proSdrViirsCalFPredictedTableLUT {
    signed char Fit_type[NUM_VIIRS_BAND];
    unsigned char implicit_pad0[2];
    long long T_ref;
    double F_ref[NUM_VIIRS_BAND][MAX_NUM_DETECTOR][MAX_NUM_GAIN][NUM_MIRROR_SIDES];
    double F_param_1[NUM_VIIRS_BAND][MAX_NUM_DETECTOR][MAX_NUM_GAIN][NUM_MIRROR_SIDES];
    double F_param_2[NUM_VIIRS_BAND][MAX_NUM_DETECTOR][MAX_NUM_GAIN][NUM_MIRROR_SIDES];
};

// Gain Table

struct proSdrViirsCalGainTableLUT {
    double data[NUM_VIIRS_BAND];
};

struct proSdrViirsCalHAMERLUT {
    int HAM_ER_indices[NUM_OF_HAMER_FILES][2];
    double HAM_ER_tp[MAX_HAM_ER_INDEX];
    double HAM_ER_rad[MAX_HAM_ER_INDEX];
};

struct proSdrViirsCalOBCERLUT {
    int OBC_ER_indices[NUM_OF_OBCER_FILES][2];
    double OBC_ER_tp[MAX_OBC_ER_INDEX];
    double OBC_ER_rad[MAX_OBC_ER_INDEX];
};

struct proSdrViirsCalOBCRRLUT {
    int OBC_RR_indices[NUM_OF_OBCRR_FILES][2];
    double OBC_RR_tp[MAX_OBC_RR_INDEX];
    double OBC_RR_rad[MAX_OBC_RR_INDEX];
};

// Reflective Lookup Table Structure

struct proSdrViirsCalObsToPixelsLUT {
    int pixels[EV_FRAMES_750m_DG];
};

struct ProSdrViirsCalOnboardDnbOffsetsLUT {
    short offsets[DNB_VIIRS_RDR_COLS][NUM_DETECTORS_DNB][N_DNB_BANDS];
};

struct proSdrViirsCalRMParametersLUT {
    int Telec_Therm_Indexes[NUM_TELEC_THERM];
    float Telec_Therm_Weights[NUM_TELEC_THERM];
    int Tsh_Indexes[NUM_TSH];
    float Tsh_Weights[NUM_TSH];
    int Ttele_Indexes[NUM_TTELE];
    float Ttele_Weights[NUM_TTELE];
    float Ttele_Offset;
    int Trta_Indexes[NUM_TRTA];
    float Trta_Weights[NUM_TRTA];
    float Trta_Offset;
    int Tcav_Indexes[NUM_TCAV];
    float Tcav_Weights[NUM_TCAV];
    int Tmir_Indexes[NUM_TMIR];
    float Tmir_Weights[NUM_TMIR];
    int Tfpsm_Indexes[NUM_TFPSM];
    float Tfpsm_Weights[NUM_TFPSM];
    int Tfplw_Indexes[NUM_TFPLW];
    float Tfplw_Weights[NUM_TFPLW];
    int Tomm_Thermister_Indexes[NUM_TOMM_THERM];
    float Tomm_Thermister_Weights[NUM_TOMM_THERM];
    float Ttel_Tsh_Tcav_Weights[NUM_TTEL_TSH_TCAV];
    float Tmax;
    float Tmin;
    int Tomm_for_Tfpsm_Switch;
    int Tomm_for_Tfplw_Switch;
    float bbNominalTemp_Tolerances[OneTempTwoTolerances];
    float lwirNominalTolerance;
};

struct proSdrViirsCalReflectiveLUT {
    short DN_obc_avg_first_frame;
    short DN_obc_avg_num_frames;
    short RSB_SV_DN_moon_include_frames;
    short wordBoundaryPad;
};

struct ProSdrViirsCalRelativeSpectralResponseLUT {
    double wavelength[NUM_REFL_PLUS_DNB_BANDS][MAX_RSR_VALUES];
    double rsr[NUM_REFL_PLUS_DNB_BANDS][MAX_RSR_VALUES];
    int numRsr[NUM_REFL_PLUS_DNB_BANDS];
    unsigned char implicit_pad0[4];
};

struct ProSdrViirsCalSolarSpectralIradLUT {
    double solar_spectral_irad[NUM_REFL_PLUS_DNB_BANDS];
};

struct proSdrViirsCalRTAERLUT {
    int RTA_ER_indices[NUM_OF_RTAER_FILES][2];
    double RTA_ER_tp[MAX_RTA_ER_INDEX];
    double RTA_ER_rad[MAX_RTA_ER_INDEX];
};

struct proSdrViirsCalRVSLUT {
    float RVS_375m[NUM_375M_BANDS][NUM_DETECTORS_375M][EV_375M_FRAMES][NUM_MIRROR_SIDES];
    float RVS_750m_SG[NUM_750M_SG_BANDS][NUM_DETECTORS_750M][EV_750M_SG_FRAMES][NUM_MIRROR_SIDES];
    float RVS_750m_DG[NUM_750M_DG_BANDS][NUM_DETECTORS_750M][EV_750M_DG_FRAMES][NUM_MIRROR_SIDES];
    float RVS_375m_SV[NUM_375M_BANDS][NUM_DETECTORS_375M][NUM_MIRROR_SIDES];
    float RVS_750m_SV_SG[NUM_750M_SG_BANDS][NUM_DETECTORS_750M][NUM_MIRROR_SIDES];
    float RVS_750m_SV_DG[NUM_750M_DG_BANDS][NUM_DETECTORS_750M][NUM_MIRROR_SIDES];
    float RVS_375m_BB[NUM_375M_BANDS][NUM_DETECTORS_375M][NUM_MIRROR_SIDES];
    float RVS_750m_BB_SG[NUM_750M_SG_BANDS][NUM_DETECTORS_750M][NUM_MIRROR_SIDES];
    float RVS_750m_BB_DG[NUM_750M_DG_BANDS][NUM_DETECTORS_750M][NUM_MIRROR_SIDES];
};

// Solar Irad Table

struct proSdrViirsCalSolarIradLUT {
    double x[X_Y_SIZE];
    double y[X_Y_SIZE];
};

struct proSdrViirsCalTeleCoeffLUT {
    float teleCoeffs[NUM_MAX_THERMISTORS][NUM_TEMP_COEFFS];
    float bbTempOffsetDef;
    float bbTempGainDef;
    float bbTempAdditive;
    float bbTempSpare[3];
    float Spare1[12];
    float filterThreshold[NUM_MAX_THERMISTORS];
    float Spare2[3];
    float defaultValue[NUM_MAX_THERMISTORS];
    float evCtPrecTrefMux1Ca[NUM_MAX_EV_CT_PREC_TREF_MUX1CA];
};

typedef struct proSdrViirsCalBBTempCoeffStruct {
    double I0;
    double V0;
    double Rp;
    double G;
    double Const1;
} proSdrViirsCalBBTempCoeffStruct;

typedef struct proSdrViirsCalBBTempCoeffs {
    proSdrViirsCalBBTempCoeffStruct thermister[NUM_BB_THERMISTORS];
} proSdrViirsCalBBTempCoeffs;

#endif /* VIIRSLUTSTRUCTURES_H_ */
