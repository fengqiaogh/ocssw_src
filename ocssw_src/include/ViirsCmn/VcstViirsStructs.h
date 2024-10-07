/******************************************************************************
 * NAME: ViirsStructures.h
 *
 * DESCRIPTION: Header file that consolidates structure and constant definitions
 * used by the VIIRS geolocation, OBC, and calibration processes.
 *
 *
 * Created on: Nov 13, 2014
 *     Author: Sam Anderson, VCST
 *
 ******************************************************************************/

#ifndef VIIRSSTRUCTURES_H_
#define VIIRSSTRUCTURES_H_

#define Viirs_Bands     22

// L1A dimensions

#define Number_of_Scans_5min 169
#define SC_Diary_Records_5min 320
#define Number_of_Scans_6min 203
#define SC_Diary_Records_6min 381
// don's short L1A file
//#define Number_of_Scans 		34

#define Number_of_Scans   203
#define SC_Diary_Records   381
#define SC_Diary_Records_1Hz   380
#define SC_Diary_Records_10Hz  3801

#define Mband_detectors   16
#define Iband_detectors   32
#define Mband_Pixels   3200
#define Iband_Pixels   6400
#define DNB_Pixels    4064
#define Mband_Samples   6304
#define Mband_Cal_Samples  48
#define Iband_Cal_Samples  96
#define DNB_Cal_Samples   64
#define Number_of_Mbands  16
#define Number_of_Ibands  5
#define Number_of_Tbands  7
#define Number_of_MTbands  5
#define Number_of_ITbands  2
#define Number_of_Rbands  14
#define Number_of_MRbands  11
#define Number_of_IRbands  3
#define Number_of_DNBs   1   // old
#define Agg_Zones     3
#define Parity      2
#define High_Gain         0
#define Low_Gain     1
#define Gain_States    2
#define C_Coefs     4
//#define Number_of_DNBs			3		// new
#define Quaternion_Elements  4
#define Vector_Elements   3
#define EV_APIDs    24
#define HR_Metadata    146
#define Cal_Metadata   134
#define Eng_Status    8
#define Eng_Block    128
#define ASP_Offsets    3072
#define SDSM_Data    256
#define Encoder_Reading   1290
#define Mirror_Sides    2
#define Number_of_DG_bands   7
#define Number_of_BB_temps  6
#define SDSM_Samples   5
#define SDSM_Detectors   8

#define M13_LUT_Radiance_Values  327681
#define LUT_Radiance_Values  65536
#define Number_of_Iband_Lines  Iband_detectors*Number_of_Scans   //5408
#define Number_of_Iband_Pixels  6400
#define Number_of_Mband_Lines Mband_detectors*Number_of_Scans   //2704
#define Number_of_Mband_Pixels 3200
#define Number_of_DNB_Lines  Mband_detectors*Number_of_Scans   //2704
#define Number_of_DNB_Pixels 4064

// Compression

const bool bShuffleFilter = true;
const bool bDeflateFilter = true;
const int deflateLevel = 5;

// Enumerations

enum VIIRS_PLATFORM_ENUM {
    NPP,
    J1,
    J2,
    PLATFORM_MAXIMUM
};

enum GRAN_SEQ_ENUM {
    BEFORE, CURRENT, AFTER, GRAN_SEQ_MAXIMUM
};

enum VIIRS_CATEGORY_ENUM {
    ALL_BANDS,
    IMG_BANDS,
    MOD_BANDS,
    DNB_BAND,
    RSB_BANDS,
    TEB_BANDS,
    CATEGORY_MAXIMUM
};

enum VIIRS_BAND_ENUM {
    IMG_1, IMG_2, IMG_3, IMG_4, IMG_5,
    MOD_1, MOD_2, MOD_3, MOD_4, MOD_5,
    MOD_6, MOD_7, MOD_8, MOD_9, MOD_10,
    MOD_11, MOD_12, MOD_13, MOD_14, MOD_15,
    MOD_16, DNB_, BAND_MAXIMUM
};

enum VIIRS_I_BAND_ENUM {
    I_1, I_2, I_3, I_4, I_5, I_BAND_MAXIMUM
};

enum VIIRS_M_BAND_ENUM {
    M_1, M_2, M_3, M_4, M_5, M_6, M_7, M_8,
    M_9, M_10, M_11, M_12, M_13, M_14, M_15, M_16,
    M_BAND_MAXIMUM
};

enum VIIRS_RSB_BAND_ENUM {
    RSB_I1, RSB_I2, RSB_I3, RSB_M1, RSB_M2, RSB_M3,
    RSB_M4, RSB_M5, RSB_M6, RSB_M7, RSB_M8, RSB_M9,
    RSB_M10, RSB_M11, RSB_BAND_MAXIMUM
};

enum VIIRS_TEB_BAND_ENUM {
    TEB_I4, TEB_I5, TEB_M12, TEB_M13, TEB_M14, TEB_M15,
    TEB_M16, TEB_BAND_MAXIMUM
};

enum VIIRS_MOD_SG_BAND_ENUM {
    SG_M6, SG_M8, SG_M9, SG_M10, SG_M11,
    SG_M12, SG_M14, SG_M15, SG_M16, SG_BAND_MAXIMUM
};

enum VIIRS_MOD_DG_BAND_ENUM {
    DG_M1, DG_M2, DG_M3, DG_M4, DG_M5, DG_M7, DG_M13, DG_BAND_MAXIMUM
};

enum VIIRS_IMG_TEB_BAND_ENUM {
    I_TEB_I4, I_TEB_I5, I_TEB_BAND_MAXIMUM
};

enum VIIRS_MOD_RSB_BAND_ENUM {
    M_RSB_M1, M_RSB_M2, M_RSB_M3,
    M_RSB_M4, M_RSB_M5, M_RSB_M6, M_RSB_M7, M_RSB_M8, M_RSB_M9,
    M_RSB_M10, M_RSB_M11, M_RSB_BAND_MAXIMUM
};

enum VIIRS_MOD_TEB_BAND_ENUM {
    M_TEB_M12, M_TEB_M13, M_TEB_M14, M_TEB_M15,
    M_TEB_M16, M_TEB_BAND_MAXIMUM
};

enum VIIRS_SENSOR_MODES {
    VIIRS_SENSOR_LAUNCH = 0,
    VIIRS_SENSOR_ACTIVATION,
    VIIRS_SENSOR_OG,
    VIIRS_SENSOR_DIAGNOSTIC,
    VIIRS_SENSOR_OP_DAY,
    VIIRS_SENSOR_OP_NIGHT,
    VIIRS_SENSOR_SAFE,
    VIIRS_SENSOR_INVALID
};

// L1A Earth View and Calibration data structures

struct ImgL1ADataType {
    short SV_I[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    short BB_I[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    short SD_I[Number_of_Ibands][Number_of_Scans][Iband_detectors][Iband_Cal_Samples];
    short EV_I[Number_of_Scans][Iband_detectors][Iband_Pixels];
};

struct ModL1ADataType {
    short SV_M[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    short BB_M[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    short SD_M[Number_of_Mbands + 1][Number_of_Scans][Mband_detectors][Mband_Cal_Samples];
    short EV_M_SG[Number_of_Scans][Mband_detectors][Mband_Pixels];
    short EV_M_DG[Number_of_Scans][Mband_detectors][Mband_Samples];
};

struct DnbL1ADataType {
    short SV_DNB[Number_of_DNBs][Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    short BB_DNB[Number_of_DNBs][Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    short SD_DNB[Number_of_DNBs][Number_of_Scans][Mband_detectors][DNB_Cal_Samples];
    short EV_DNB[Number_of_Scans][Mband_detectors][Number_of_DNB_Pixels];
};

// L1A Engineering data structures

// analog signal processor data

struct ASPAnalogDataType {
    short p30vPin;
    unsigned char aspBiases[60];
    unsigned char other[25];
    unsigned char reserved[41];
};

// analog signal processor hardware configuration data .. should be 128 bytes

struct ASPConfigL1ADataType {
    // revision for VISNIR regiser table
    unsigned char ap_visnir_reg_tbl_rev;

    unsigned char ap_smwir_reg_tbl_rev; // revision for SMWIR regiser table
    unsigned char ap_lwir_reg_tbl_rev; // revision for LWIR regiser table
    unsigned char ap_fp_ideal_tbl_rev; // revision for FPA ideal table

    // LW imaging integration ticks by imaging frame sync width
    unsigned char ap_lw_ifs_width : 6;
    unsigned char ap_det_connected : 1; // ASP Detector state
    unsigned char ap_dc_restore : 1; // ASP global DC restore

    // SM radiometric integration ticks by radiometric frame sync width
    unsigned char ap_lw_rfs_width : 6;
    unsigned char ap_dual_gain_3_ops : 2; // ASP dual gain ops

    // SM radiometric integration ticks by radiometric frame sync width
    unsigned char ap_sm_ifs_width : 6;
    unsigned char ap_fpa_st_3_ops : 2; // ASP FPA self test ops

    // SM radiometric integration ticks by radiometric frame sync width
    unsigned char ap_sm_rfs_width : 6;
    unsigned char fill1 : 2;

    // VN imaging integration ticks by imaging frame sync width
    unsigned char ap_vn_ifs_width : 6;
    unsigned char fill2 : 2;

    // VN imaging integration ticks by imaging frame sync width
    unsigned char ap_vn_rfs_width : 6;
    unsigned char fill3 : 2;

    // Focal Plane Temperature LW CFPA high resolution temperature
    short ft_lw_cfpa_hi_rsl;

    // Focal Plane Temperature LW CFPA wide range temperature
    short ft_lw_cfpa_lo_rsl;

    // Focal Plane Temperature SM CFPA high resolution temperature
    short ft_sm_cfpa_hi_rsl;

    // Focal Plane Temperature SM CFPA wide range temperature
    short ft_sm_cfpa_lo_rsl;

    // Focal Plane Temperature VIS/NIR FPA temperature
    short ft_vis_nir_fpa;

    unsigned char fill4[108]; // 864 bits of fill data
};

// analog signal processor offset data

struct ASPOffsetL1ADataType {
    // for 0 ham side
    short m1_ham0[16];
    short m2_ham0[16];
    short m3_ham0[16];
    short m4_ham0[16];
    short m5_ham0[16];
    short m6_ham0[16];
    short m7_ham0[16];
    short m8_ham0[16];
    short m9_ham0[16];
    short m10_ham0[16];
    short m11_ham0[16];
    short m12_ham0[16];
    short m13_ham0[16];
    short m14_ham0[16];
    short m15_ham0[16];
    short m16a_ham0[16];
    short m16b_ham0[16];
    short i1e_ham0[16];
    short i1o_ham0[16];
    short i2e_ham0[16];
    short i2o_ham0[16];
    short i3e_ham0[16];
    short i3o_ham0[16];
    short i4e_ham0[16];
    short i4o_ham0[16];
    short i5e_ham0[16];
    short i5o_ham0[16];

    // for 1 ham side
    short m1_ham1[16];
    short m2_ham1[16];
    short m3_ham1[16];
    short m4_ham1[16];
    short m5_ham1[16];
    short m6_ham1[16];
    short m7_ham1[16];
    short m8_ham1[16];
    short m9_ham1[16];
    short m10_ham1[16];
    short m11_ham1[16];
    short m12_ham1[16];
    short m13_ham1[16];
    short m14_ham1[16];
    short m15_ham1[16];
    short m16a_ham1[16];
    short m16b_ham1[16];
    short i1e_ham1[16];
    short i1o_ham1[16];
    short i2e_ham1[16];
    short i2o_ham1[16];
    short i3e_ham1[16];
    short i3o_ham1[16];
    short i4e_ham1[16];
    short i4o_ham1[16];
    short i5e_ham1[16];
    short i5o_ham1[16];

    char fill[1344];
};

// DNB configuration data .. should be 128 bytes

struct DNBConfigL1ADataType {
    unsigned char dp_aggreg_mode_rev;
    unsigned char dp_threshold_tbl_rev;
    unsigned char dp_1a_offsets_tbl_rev;
    short dp_dnb_ccd;
    unsigned char dp_1b_offsets_tbl_rev;
    unsigned char dp_2_offsets_tbl_rev;
    unsigned char dp_3_offsets_tbl_rev;
    unsigned char fill[120];
};

// digital pre-processor configuration data .. should be 128 bytes

const int dpp_config_index_1 = 0;

struct DPPConfigL1ADataType1 {
    unsigned char dp_reg_tbl_rev;
    unsigned char dp_state_tran_tbl_rev;
    unsigned char dp_band_proc_tbl_rev;
    unsigned char dp_heat_ctrl_tbl_rev;
    unsigned char dp_macro_cmd_tbl_rev;
    unsigned char dp_crit_tele_tbl_rev;
    unsigned short dp_stor_cmd_tbl_rev;

    unsigned char ps_sec_e_isog_on : 1;
    unsigned char ps_sec_d_csog_on : 1;
    unsigned char ps_sec_c_se_on : 1;
    unsigned char ps_sec_b_apfp_on : 1;
    unsigned char dp_servo_in_use : 1;
    unsigned char dp_nonrdt_fpie_pwr : 1;
    unsigned char dp_hrd_pkt_norm_test : 1;
    unsigned char dp_dn_m_l_gain_pkt : 1;

    unsigned char se_b_tele_pos_known : 1;
    unsigned char se_b_mtrs_stopped : 1;
    unsigned char se_b_mtr_coil_driver : 1;
    unsigned char se_b_anlg_pwr_on : 1;
    unsigned char se_a_tele_pos_known : 1;
    unsigned char se_a_mtrs_stopped : 1;
    unsigned char se_a_mtr_coil_driver : 1;
    unsigned char se_a_anlg_pwr_on : 1;

    unsigned char dp_ap_m16_tdi_on : 1;
    unsigned char cp_blk_pwr_sel : 1;
    unsigned char dp_dn_aggreg_mod : 6;

    unsigned char dp_scan_encdr_delta;
    unsigned char dummy;

    unsigned char se_servo_pwr_sel : 1;
    unsigned char dp_ap_self_test : 1;
    unsigned char spare1 : 6;

    unsigned char se_b_ham_mir_side : 1;
    unsigned char se_a_ham_mir_side : 1;
    unsigned char ap_dc_fast_restor : 1;
    unsigned char dp_dnb_dark_sub_eth : 1;
    unsigned char dp_dnb_dark_sub_cal : 1;
    unsigned char dp_dnb_tmg_mode : 1;
    unsigned char dp_dnb_1a_1b_stage : 2;
};

const int dpp_config_index_2 = 15;

struct DPPConfigL1ADataType2 {
    short se_a_ham_mtr_curr;
    short se_a_tele_mtr_curr;
    short se_b_ham_mtr_curr;
    short se_b_tele_mtr_curr;
    short ct_prec_tref_mux1ca1;
    short ct_prec_tref_mux1ca2;
    short ct_prec_tref_mux1ca3;
    short ft_adc_ref;
    short ct_adc_ref_lw_stpt;
    short ft_ckt_gnd;
    short ft_lw_cfpa_htr_pwr;
    short ft_lw_setpt_ref;
    short ft_sm_cfpa_htr_pwr;
    short ft_sm_setpt_ref;
    short se_a_ham_rate_error;
    short se_a_tele_rate_error;
    short se_b_ham_rate_error;
    short se_b_tele_rate_error;
};

// solar diffuser ...
//struct SDSML1ADataType
//{
//   unsigned short   position : 8;
//   short   			samples[5][8]; // Data is sample by AMP, 8 SDSM detector, 5 samples
//   // the following combines the "other" and "reserved" field
//   // described in the MDFCB for a total of 1400 bits
//   short    		preamp;
//   unsigned char    other[173];
//};

// temperature data

struct EngTempL1ADataType {
    short bbTemp[6];
    short ham_tmp1; // half angle mirror t1
    short ham_tmp2; // half angle mirror t2
    short hm_cr_cs_prt; // CR cold stage PRT temperature
    short hm_cr_is_prt; // CR intermediate stage PRT temperature
    short hm_cr_os_prt; // CR outer stage PRT temperature

    unsigned short ap_vn_vref_inhibit : 12;
    unsigned short fill1 : 4;

    unsigned int fill2;
    short hm_ham_db_op_htrA;
    short hm_ham_db_op_htrB;
    short hm_tele_db_op_htrA;
    short hm_tele_db_op_htrB;

    short mf_ao_blkhd_nx_pz; // Aft optics
    short mf_ao_blkhd_px_nz; // Aft optics
    short mf_ao_km2_nx; // Aft optics KM2_NX
    short mf_ao_km2_ny;
    short mf_ao_km3_px;
    short mf_fold_mirror_blkhd;
    short mf_ham_blkhd; // half angle mirror bulkhead
    short mf_km2_nxny; // NXNY kinematic mount, 2 deg contraint
    short mf_km2_nxpy; // NXPY kinematic mount, 2 deg contraint
    short mf_km1_pxny; // PXNY kinematic mount, 1 deg contraint
    short mf_km2_pxpy; // PXPY kinematic mount, 2 deg contraint
    short mf_nadir_rad_nxp; // Nadir radiator NXPY
    short mf_nadir_rad_pxn; // Nadir radiator PXNY
    short mf_scan_cavity_nxn;
    short mf_scan_cavity_nxp;
    short mf_scan_cavity_baf_nz; // scan cavity baffle NZ
    short mf_scan_cavity_baf_pz; // scan cavity baffle PZ
    short mf_scan_cavity_bknd_n;
    short mf_space_rad_nz; // space radiator NZ
    short mf_space_rad_pz; // space radiator PZ
    short mf_stopassy_baff_nz;
    short mf_tel_blkhd_nypz;
    short mf_tel_blkhd_py;

    short se_a_hammtr_dfbear; // HAM DF Bearing Temp MTR_A
    short se_a_telemtr_dfbear; // Tele DF Bearing Temp MTR_A
    short se_b_hammtr_dfbear; // HAM DF Bearing Temp MTR_B
    short se_b_telemtr_dfbear; // Tele DF Bearing Temp MTR_B
    short ap_lw_cca; // Long wave IR CCA
    short ap_sm_cca; // Short/Med wave IR CCA
    short ap_vn_cca; // Visible/Near IR CCA

    short dp_dnb_cca; // DNB CCA
    short dp_dpp_cca; // Digital preprocessor CCA
    short dp_fpie_ad_cca; // FPIE A/D CCA
    short dp_fpie_clk_cca; // FPIE CLK CCA
    short ft_cca;
    short power_supply1;
    short power_supply2;
    short se_a_cca;
    short se_b_cca;

    unsigned char hm_ham_op_htr_cntl : 1; // FSW echo C_HM09 for HAM heater
    unsigned char hm_dnb_htr_cntl : 1; // FSW echo C_HM06A for DNB heater
    unsigned char ft_sm_heater : 1; // short/med wave IR heater, 0=off, 1=on
    unsigned char ft_lw_heater : 1; // Long wave IR heater, 0=off, 1=on
    unsigned char ft_sm_82k : 1; // short/med wave 80K setpoint, 0=off, 1=on
    unsigned char ft_sm_80k : 1; // short/med wave 78K setpoint, 0=off, 1=on
    unsigned char ft_lw_82k : 1; // long wave 80K setpoint, 0=off, 1=on
    unsigned char ft_lw_80k : 1; // long wave 78K setpoint, 0=off, 1=on

    unsigned char hm_tel_op_htr_cntl : 1; // FSW echo C_HM10 for Tele Op heater
    unsigned char filler1 : 7;

    unsigned short hm_ham_ophtr_temp : 14; // indicates heater setpoint
    unsigned char filler2 : 2;

    unsigned short hm_tel_ophtr_temp : 14; // indicates heater setpoint
    unsigned char filler3 : 2;

    unsigned short ap_bb_avg_temp : 14; // average of the BB thermistor         : 2;
    unsigned char filler4 : 2;

    unsigned short bb_htr_temp : 14; // FSW echo C_BB01
    unsigned char filler5 : 2;

    unsigned short bb_select; // indicates which BB thermistor used to get avg
    unsigned short filler6;
};

const int bb_htr_temp_index = 122;

struct EncoderReadingL1ADataType {
    unsigned short teleEncoding[1290];
    unsigned short hamEncoding[1290];
};

struct EngReservedL1ADataType {
    unsigned short telAngleStart;
    unsigned short hamAngScanStart;
    unsigned short hrdValues[16];
    unsigned char other[3];
    unsigned char reserved[89];
};

// engineering status data

struct EngStatusL1ADataType {
    unsigned char fpaCalDataGain : 1;
    unsigned char sdsmActive : 1;
    unsigned char hamSide : 1;
    unsigned char servoLocked : 1;
    unsigned char sensorMode : 4;

    unsigned char es_se_b_teleham_scansyn : 1;
    unsigned char es_se_a_teleham_scansyn : 1;
    unsigned char testDataPattern : 4;
    unsigned char spare : 2;

    unsigned char es_sd_sdsm_mtr_step_cnt : 7;
    unsigned char other : 1;

    unsigned char reserved[5];
};

// L1B program data structures

struct ScanLineL1ADataType {
    double scan_start_time[Number_of_Scans];
    double scan_end_time[Number_of_Scans];
    short scan_start_CCSDS_day[Number_of_Scans];
    int scan_start_CCSDS_msec[Number_of_Scans];
    short scan_start_CCSDS_usec[Number_of_Scans];
    short scan_end_CCSDS_day[Number_of_Scans];
    int scan_end_CCSDS_msec[Number_of_Scans];
    short scan_end_CCSDS_usec[Number_of_Scans];
    int VIIRS_scan_number[Number_of_Scans];
    unsigned char ham_side[Number_of_Scans]; // new version
    unsigned char sensor_mode[Number_of_Scans]; // new version
    //	short 			ham_side[Number_of_Scans];         	  // old version
    //	short			sensor_mode[Number_of_Scans];         // old version
    unsigned char HR_metadata[Number_of_Scans][EV_APIDs][HR_Metadata];
    unsigned char cal_metadata[Number_of_Scans][Cal_Metadata];
};

struct EngineeringL1ADataType {
    unsigned char eng_status[Number_of_Scans][Eng_Status];
    unsigned char eng_temp[Number_of_Scans][Eng_Block];
    unsigned char eng_reserved[Number_of_Scans][Eng_Block];
    unsigned char dnb_config[Number_of_Scans][Eng_Block];
    unsigned char dpp_config[Number_of_Scans][Eng_Block];
    unsigned char asp_config[Number_of_Scans][Eng_Block];
    unsigned char asp_offsets[Number_of_Scans][ASP_Offsets];
    unsigned char asp_analog[Number_of_Scans][Eng_Block];
    unsigned char sdsm_data[Number_of_Scans][SDSM_Data];
    unsigned short ham_encoder[Number_of_Scans][Encoder_Reading];
    unsigned short tel_encoder[Number_of_Scans][Encoder_Reading];
};

/**
 *  isImg(), isTeb(), isRSB(), isDG()
 *
 * Band categories identification
 */

inline bool isImg(VIIRS_BAND_ENUM band) {
    return ( band < MOD_1);
}

inline bool isMod(VIIRS_BAND_ENUM band) {
    return ((band > IMG_5) &&
            (band < DNB_));
}

inline bool isDG(VIIRS_BAND_ENUM band) {
    return ((band == MOD_1) ||
            (band == MOD_2) || (band == MOD_3) || (band == MOD_4) ||
            (band == MOD_5) || (band == MOD_7) || (band == MOD_13));
}

inline bool isDG(VIIRS_M_BAND_ENUM band) {
    return ((band == M_1) ||
            (band == M_2) || (band == M_3) || (band == M_4) ||
            (band == M_5) || (band == M_7) || (band == M_13));
}

inline bool isTeb(VIIRS_BAND_ENUM band) {
    return ((band == IMG_4) ||
            (band == IMG_5) || (band == MOD_12) || (band == MOD_13) ||
            (band == MOD_14) || (band == MOD_15) || (band == MOD_16));
}

/**
 *  ib_index(), mb_index(), dg_index(), teb_index(), rsb_index(),
 *  iteb_index(), mteb_index(), irsb_index(), mrsb_index
 *
 * Band categories identification
 */

inline int ib_index(VIIRS_BAND_ENUM band) {
    if (isImg(band)) return band;
    else return -1;
}

inline int mb_index(VIIRS_BAND_ENUM band) {
    if (isMod(band)) return (band - MOD_1);
    else return -1;
}

inline int sg_index(VIIRS_BAND_ENUM band) {
    if (isMod(band)) {
        if (band == MOD_6) return 0;
        else if (band == MOD_8) return 1;
        else if (band == MOD_9) return 2;
        else if (band == MOD_10) return 3;
        else if (band == MOD_11) return 4;
        else if (band == MOD_12) return 5;
        else if (band == MOD_14) return 6;
        else if (band == MOD_15) return 7;
        else if (band == MOD_16) return 8;
        else return -1;
    } else return -1;
}

inline int dg_index(VIIRS_BAND_ENUM band) {
    if (isMod(band)) {
        if (band <= MOD_5) return mb_index(band);
        else if (band == MOD_7) return 5;
        else if (band == MOD_13) return 6;
        else return -1;
    } else return -1;
}

inline int teb_index(VIIRS_BAND_ENUM band) {
    if (band == IMG_4) return 0;
    else if (band == IMG_5) return 1;
    else if (band == MOD_12) return 2;
    else if (band == MOD_13) return 3;
    else if (band == MOD_14) return 4;
    else if (band == MOD_15) return 5;
    else if (band == MOD_16) return 6;
    else return -1;
}

inline int rsb_index(VIIRS_BAND_ENUM band) {
    if (band == IMG_1) return 0;
    else if (band == IMG_2) return 1;
    else if (band == IMG_3) return 2;
    else if (band == MOD_1) return 3;
    else if (band == MOD_2) return 4;
    else if (band == MOD_3) return 5;
    else if (band == MOD_4) return 6;
    else if (band == MOD_5) return 7;
    else if (band == MOD_6) return 8;
    else if (band == MOD_7) return 9;
    else if (band == MOD_8) return 10;
    else if (band == MOD_9) return 11;
    else if (band == MOD_10) return 12;
    else if (band == MOD_11) return 13;
    else return -1;
}

inline int iteb_index(VIIRS_BAND_ENUM band) {
    if (band == IMG_4) return 0;
    else if (band == IMG_5) return 1;
    else return -1;
}

inline int irsb_index(VIIRS_BAND_ENUM band) {
    if (band == IMG_1) return 0;
    else if (band == IMG_2) return 1;
    else if (band == IMG_3) return 2;
    else return -1;
}

inline int mteb_index(VIIRS_BAND_ENUM band) {
    if (band == MOD_12) return 0;
    else if (band == MOD_13) return 1;
    else if (band == MOD_14) return 2;
    else if (band == MOD_15) return 3;
    else if (band == MOD_16) return 4;
    else return -1;
}

inline int mrsb_index(VIIRS_BAND_ENUM band) {
    if (band == MOD_1) return 0;
    else if (band == MOD_2) return 1;
    else if (band == MOD_3) return 2;
    else if (band == MOD_4) return 3;
    else if (band == MOD_5) return 4;
    else if (band == MOD_6) return 5;
    else if (band == MOD_7) return 6;
    else if (band == MOD_8) return 7;
    else if (band == MOD_9) return 8;
    else if (band == MOD_10) return 9;
    else if (band == MOD_11) return 10;
    else return -1;
}

#endif /* VIIRSSTRUCTURES_H_ */
