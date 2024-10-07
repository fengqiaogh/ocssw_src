/**************************************************************************
 *
 * NAME: VcstParamsReader
 *
 * DESCRIPTION: Singleton object class. Reads file path data from program
 * configuration file (pcf) and makes it available to applications via
 * various get commands.
 *
 *  Created on: Aug 25, 2014
 *      Author: Sam Anderson, VCST
 *
 **************************************************************************/

#ifndef _VcstParamsReader_h_
#define _VcstParamsReader_h_

#include <string>
#include <iostream>
#include <vector>

#include <clo.h>

extern const std::string VIIRS_PLATFORM;
extern const std::string PLATFORM_NPP;
extern const std::string PLATFORM_J1;
extern const std::string PLATFORM_J2;

//
// SDR CmnGeolocation short name constants
//
extern const std::string CMNGEO_JPL_EPHEM;
extern const std::string CMNGEO_USNO_PW_UT1;
extern const std::string CMNGEO_SAA_COEFF;
extern const std::string CMNGEO_ANOMALY_TLE;
extern const std::string CMNGEO_PARAM_LUT;
extern const std::string CMNGEO_PLATFORM_LUT;

//
// Geolocation Resolution specific parameter luts
//
extern const std::string GEO_DNB_PARAM;
extern const std::string GEO_IMG_PARAM;
extern const std::string GEO_MOD_PARAM;

//
// SDR short name constants for Terrain Correction TERECO tiles
//
extern const std::string ECO_TILE;
extern const std::string LWM_PATH;
extern const std::string TILE_ID_METADATA;

//
// Thresholds and LUT's needed for processing VIIRS SDR
//
extern const std::string VIIRS_DG_ANOMALY_DN_LIMITS_LUT;
extern const std::string VIIRS_DNB_STRAY_LIGHT_LUT;
extern const std::string VIIRS_DNB_STRAY_LIGHT_CORRECTION_LUT;
extern const std::string VIIRS_SDR_DNB_DN0_LUT;
extern const std::string VIIRS_SDR_DNB_RVS;
extern const std::string VIIRS_SDR_DNB_FRAME_TO_ZONE;
extern const std::string VIIRS_SDR_F_PREDICTED_LUT;
extern const std::string VIIRS_SDR_DNB_F_PREDICTED_LUT;
extern const std::string VIIRS_SDR_GAIN_LUT;
extern const std::string VIIRS_SDR_HAM_ER_TABLE;
extern const std::string VIIRS_SDR_RTA_ER_TABLE;
extern const std::string VIIRS_SDR_OBC_ER_TABLE;
extern const std::string VIIRS_SDR_OBC_RR_TABLE;
extern const std::string VIIRS_SDR_EBBT_TABLE;
extern const std::string VIIRS_SDR_TELE_COEFFS;
extern const std::string VIIRS_SDR_SOLAR_IRAD_LUT;
extern const std::string VIIRS_SDR_RSR_LUT;
extern const std::string VIIRS_SDR_OBS_TO_PIXELS;
extern const std::string VIIRS_SDR_RADIOMETRIC_PARAMETERS;
extern const std::string VIIRS_SDR_QA_LUT;
extern const std::string VIIRS_SDR_EMISSIVE_LUT;
extern const std::string VIIRS_SDR_REFLECTIVE_LUT;
extern const std::string VIIRS_SDR_RVS_LUT;
extern const std::string VIIRS_SDR_BB_TEMP_COEFFS;
extern const std::string VIIRS_SDR_DNB_C_COEFFS;
extern const std::string VIIRS_SDR_DELTA_C_LUT;
extern const std::string VIIRS_SDR_COEFF_A_LUT;
extern const std::string VIIRS_SDR_COEFF_B_LUT;
extern const std::string VIIRS_SDR_DNB_GAIN_RATIOS_LUT;
extern const std::string VIIRS_SDR_DNB_LGS_GAINS_LUT;
extern const std::string VIIRS_SDR_RELATIVE_SPECTRAL_RESPONSE_LUT;
extern const std::string VIIRS_SDR_SOLAR_SPECTRAL_IRAD_LUT;

//
// VIIRS SDR Solar Diffuser View Processing Short Names.
//
extern const std::string VIIRS_SOLAR_DIFF_OBC_IP;
extern const std::string VIIRS_SOLAR_DIFF_GEO_IP;
extern const std::string VIIRS_SOLAR_DIFF_REFL_LUT;
extern const std::string VIIRS_SOLAR_DIFF_RVS_LUT;
extern const std::string VIIRS_SOLAR_DIFF_PROC_COEFFS;
extern const std::string VIIRS_SOLAR_DIFF_HIS_AGG;
extern const std::string VIIRS_SOLAR_DIFF_AGG_LUT;
extern const std::string VIIRS_SOLAR_DIFF_HIS_AGG_AC;
extern const std::string VIIRS_SOLAR_DIFF_VOLT_LUT;
extern const std::string VIIRS_SOLAR_DIFF_ROT_MATRIX_LUT;
extern const std::string VIIRS_SOLAR_DIFF_SDSM_TIME_LUT;
extern const std::string VIIRS_SOLAR_DIFF_SDSM_BRDF_LUT;
extern const std::string VIIRS_SOLAR_DIFF_SDSM_TRANS_SCREEN_LUT;
extern const std::string VIIRS_SOLAR_DIFF_TRANS_SCREEN_LUT;

//
// VIIRS netCDF4 LUT Path
//
extern const std::string VIIRS_XML;
extern const std::string VIIRS_NETCDF_LUT_PATH;
extern const std::string VIIRS_L1A;
extern const std::string VIIRS_L1A_BEFORE;
extern const std::string VIIRS_L1A_AFTER;
extern const std::string VIIRS_LEAP_SEC_PATH;
extern const std::string VIIRS_L1B_IMG;
extern const std::string VIIRS_L1B_MOD;
extern const std::string VIIRS_L1B_DNB;
extern const std::string VIIRS_L1B_CDG;
extern const std::string VIIRS_L1B_OBC;
extern const std::string VIIRS_L1B_NAV;
extern const std::string VIIRS_GEO_IMG;
extern const std::string VIIRS_GEO_MOD;
extern const std::string VIIRS_GEO_DNB;

void VL1_set_optionList(clo_optionList_t* list);
clo_optionList_t* VL1_get_optionList();
std::string VL1_get_option(const std::string& name);
std::string VL1_get_group(const std::string& group);
std::string VL1_get_netcdf_group(const std::string& group);

void VL1_add_options(clo_optionList_t* list);
void VL1_copy_options();
std::string VL1_get_source(std::vector<std::string> sourcesList);
std::string VL1_get_history(int argc, char* argv[]);

#endif
