#ifndef  L1B_TABLES_H
#define  L1B_TABLES_H

#include  "Granule.h" /* contain hdf.h, mfhdf.h, etc, and the definition of macros */

/*
!C-INC**********************************************************************
!Description:  Header file L1B_Tables.h to be included in files where the 
               lookup table routines are defined.

!Revision History:
 $Log: L1B_Tables.h,v $
 Revision 2.29  2017-04-25 10:18:32-04  xgeng
 version number updated to 6.2.2

 Revision 2.27  2017-03-09 15:40:04-05  xgeng
 updated the version number to 6.2.0 for 6.1 reprocessing

 Revision 2.25  2016-07-22 13:23:31-04  xgeng
 Applied xtalk correction to Terra LWIR PV bands (27-30)

 Revision 2.23  2015-03-13 10:41:18-04  xgeng
 Fix the anomaly before and after the sector rotation due to the timing mismatch between the command performed and the status updated.

 Revision 2.22  2014-10-02 13:39:49-04  xgeng
 Added doi attributes to NRT product (V6.1.18).

 Revision 2.21  2014-01-29 09:11:24-05  xgeng
 Migrate to SDP Toolkit 5.2.17

 Revision 2.20  2012-06-27 15:57:04-04  xgeng
 add the doi metadata

 Revision 2.19  2011-04-29 12:26:36-04  ltan
 PGE02_VERSION updated to V6.1.12

 Revision 2.18  2011-04-07 14:41:30-04  xgeng
 1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

 Revision 2.16  2010-11-15 10:15:56-05  xgeng
 Code version is updated to 6.1.8

 Revision 2.13  2010/02/04 22:20:17  xgeng
 PGE02_VERSION updated to v6.1.6

 Revision 2.9  2009/07/24 21:03:12  xgeng
 PGE02_VERSION updated to v6.1.2

 Revision 2.7  2008/12/24 21:39:32  xgeng
 L1B Code update to V6.1.0/V6.1.1

 Revision 2.5  2008/11/18 20:20:00  xgeng
 merge branch for V6.0.0

 Revision 2.4.1.12  2008/06/05 21:53:06  xgeng
 Added a new QA table name "DET_QUAL_FLAG2_VALS_LUT_NAME"; Added "Detector_Quality_Flag2_Values" and "Detector_Quality_Flag2" in structure common_QA_tables_t; updated PGE02_VERSION

 Revision 2.4  2008/02/01 22:24:51  ltan
 PGE02_VERSION update for v5.0.40

 Revision 2.1  2006/10/30 15:05:15  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log. Comments updated.

 Revision 03.21, October 15, 2004  Razor Issue #201
 Added a new dimension of Mirror Side to array "Band_21_b1".
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 03.21, October 15, 2004  Razor Issue #199
 Added a new macro definition of "SWIR_CORR_SENDING_DETECTOR_LUT_NAME".
 Added a new one-dimension array of "SWIR_corr_sending_detector" in
 structure SWIR_correction_tables_t.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 03.20, November 7, 2003
 Increased buffer sizes for ALGORITHMPACKAGEMATURITYCODE and MISSION_PHASE
 to 15.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.13, October 31, 2003  Razor Issue #195
 Removed LUT ProcessingCenter from the QA LUT structure. Removed macros
 "MAX_PROCESSINGCENTER_BUFFER" and "PROCESSING_CENTER_LUT_NAME.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 03.12  October 24, 2003  Razor Issue #196 (formerly Issue #184)
 Added LUTs roll_threshold_angle, pitch_threshold_angle, and
 yaw_threshold_angle to QA LUT structure. Added macro
 "NUM_ATTITUDE_ANGLES".
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 03.11 March 26, 2003  Razor Issue #190
 Changed name of LUTs "B26_B5_Corr_Terra", "B26_B5_Corr_Switch_Terra",
   "B26_B5_Frame_Offset_Terra" to "B26_B5_Corr", "B26_B5_Corr_Switch",
   "B26_B5_Frame_Offset" since are now applicable to both Terra and Aqua 
   processing
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.16  March 26, 2003  Razor Issue #191
 Added SWIR OOB sending band LUT structure.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.15 Oct. 3, 2002    Razor Issue #187
 Removed R_star from refl_tables_t structure. 
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 03.14 July 2, 2002    Razor Issue #161
 Added function prototypes for BDSM_index and Expand_BDSM_LUTs.
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 03.13  June 5, 2002  Razor Issue #183
 Change type of dn_sat_ev to float64.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.12  March 25, 2002  Razor Issue #178
 Remove ADC Correction and Associated LUTs
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.11  March 14, 2002  Razor Issue #182
 Applicable to MODIS/TERRA (PFM) processing only.
 Added LUTs "B26_B5_Corr_Terra", "B26_B5_Corr_Switch_Terra",
   "B26_B5_Frame_Offset_Terra"
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.1   March 8, 2002   Razor Issue #174
 Removed LUTs "RVS_250m", "RVS_500m", "RVS_1km_RefSB", "RVS_1km_Emiss_SV", 
     "RVS_1km_Emiss_BB", "RVS_1km_Emiss_EV"
 Added LUTs "RVS_RefSB", "RVS_TEB", "RVS_BB_SV_Frame_No"
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 03.07  November 6, 2001  Razor Issue #166
 Added runtime_params to argument list for Read_Lookup_Tables 
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.06 March 26, 2001, Razor Issue 159
 Added new LUT dn_sat_ev and removed LUT DN_sat
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.05 March 5, 2001, Razor Issue 156
 Added TDLUT_PIECEWISE_LINEAR, TDLUT_ALGORITHM_ATTR_NAME, and 
    TDLUT_PIECEWISE_LINEAR_TIMES.  Added function prototype for
    TDLUT_ReadPiecewiseLinearFunction.
 Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov)
  
 
 Revision 03.04 December 19, 2000, Razor issue 148
 Changed num_overlap_scans_b1_T_bb to num_overlap_scans_b1 and
 deleted num_overlap_scans_temps.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.03 December 8, 2000, Razor issue 143
 Esun, Radiance Scales and Offsets
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.02 December 1, 2000
 Changed for new SWIR correction algorithm, issue 145.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.02 October 29, 2000
 Added common QA LUT "Control options" for Razor issue 142.
 NOTE: earlier added emissive LUT "SV_DN_moon_include_frames"
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 ... (gap) ...

 Revision 03.01 December 22, 1999
 Revised the variables and structure members.  The same variables which
 define the LUTs are also used by the LUTs generation code.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.00 November 17, 1999
 Completely revised the organization to allow for valid range checking
 and to enable ease of writing meaningful error messages.
 Add new versioning strategy.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.17 August 23, 1999
 Removed LUT INT_correction_switch
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.16 August 12, 1999
 Added macros and variables for LUTs L_Max and L_Min
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.15 April 12, 1999
 Removed A3 lookup table.
 Steven Gu(zgu@mcst.gsfc.nasa.gov)
 
 Revision 02.10 April 9, 1998
 Changed LUT "Time_gain_factor" to
 "Time_gain_factor_rad".
 Added 2 new Reflective LUTs:
 "Time_gain_factor_refl", and "L_SD".
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 April 1998
 Removed the Cleanup_Emiss_Tables(), read_SIL_All_Bands(),
 SIL_Lookup(), and SIL_Cleanup() prototypes and the SIL_t
 type definition. (they're obsolete).
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 added the SWIR_correction_tables_t
 structure and an instance of this new 
 type in the refl_tables_t structure.
 updated the refl_table dimension constants.
 updated the refl_table table names.
 got rid of DIL_t and it's supporting functions.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 updated the emiss_tables_t
 added a new structure: QA_tables_t
 updated the emiss_table dimension constants
 updated the emiss_table table names
 deleted constants from V2.0 of V vs L algo
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.00 Dec. 1996
 Summarized L1B_Tables.h, HDF_Tables.h, and part of Gen_Lib.h; 
 added new lookup tables for L1B_v2. 
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

 Revision 01.00 1996/01/29 13:54:59
 Initial development
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden(baden@highwire.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

 Currently, the maximum rank is 5, which covers almost all LUTs.

 Versioning Strategy
 -------------------
 1.  The PGE version currently represents the version of the code itself.
     When the code changes, the PGE version must also change.
     Therefore, this is set as a macro "PGE02_VERSION".  For each
     actual Level 1B delivery, of the code, we need to agree with the
     DAAC and SDST as to what the PGE version will be and then set this
     macro accordingly.  Note that this macro should be set prior to
     generating any LUTs.  The LUT generation software will also use this
     macro to set the PGE version in each LUT file (see item 2 below).
     In the code, this macro value will be checked against the value of
     PGE02_Version retrieved from the PCF file to make sure they are
     identical.
 2.  Each LUT file contains the PGE version as a LUT.  The code will
     check that the PGE version set in the LUT matches the code macro
     PGE02_VERSION.  This is to help prevent out-of-date LUT files
     from being used with a release of the code (which has happened
     numerous times outside of MCST).
 3.  Each LUT file contain the MCST version (which is set in the archive
     metadata "AlgorithmPackageVersion").  These should all agree with
     each other.  When a LUT update is delivered to SDST, only the MCST
     version changes -- the PGE version remains the same.

 Problems with Versioning Strategy
 ---------------------------------
 1.  The PGE version really includes MOD_PR02QA code as well.  If that
     code has to change, then we may have to deliver another version of L1B
     so that the PGE version is changed (???).
 2.  SDST has talked about adding a "process version" at some point.
     This could change the strategy if they ever do.

 Writing the Versions to the product
 -----------------------------------
 The value in PGE02_Version retrieved from the PCF file will be written to
 the product.

!END********************************************************************
*/

/* 
 * The following macros define look up table types, e.g. constant LUT,
 * step function time-dependent LUT, and the names of the attributes
 * of SDS LUTs.
 */

#define TDLUT_INVALID                 -1
#define TDLUT_CONSTANT                0
#define TDLUT_STEPFUNCTION            1
#define TDLUT_PIECEWISE_LINEAR        2   
#define TDLUT_ALGORITHM_ATTR_NAME     "algorithm"
#define TDLUT_STEPFUNCTION_TIMES      "times"
#define TDLUT_PIECEWISE_LINEAR_TIMES  "times"

/*
 * The following structure is used to define all the information necessary
 * to read in a LUT from an ASCII file and write it to the HDF file (for
 * LUT generation) or read the LUT from the HDF file (L1B code).
 * Some conventions that apply for all L1B LUTs:
 * 1) Strings and scalars are stored as global attributes.
 * 2) Arrays of non-character data are stored as SDSs.
 * NOTE: the structure members "ascii_file", "revdet" and "dimnames"
 * are not used in the L1B code.  They are only used in the LUT generation
 * software.
 */

typedef struct 
{
  char *name;         /* LUT name in the HDF file */
  char *ascii_file;   /* name of the ASCII file holding the LUT. */
  int32 kind;         /* attribute or SDS), using macros in L1B_Tables.h */
  int32 type;         /* data type (DFNT ...) of the LUT */
  int32 rank;         /* matrix rank (maximum = 5) */
  int32 dims[5];      /* number of values in each dimension */
  char *dimnames[5];  /* dimension names of the LUTs */
  void  *data;        /* generic pointer, dynamically assigned */
  boolean revdet;     /* flag denoting to reverse detector order */
  char  *a_lb;        /* ASCII form of the lower bound */
  char  *a_ub;        /* ASCII form of the upper bound */
  char  *a_fv;        /* ASCII form of the fill value */
} LUT_Definition_t;

/*
 * The following variables define all LUTs in each file.  The reason that
 * these are external is because the LUTs generation
 * code links to L1B_Tables and uses these variables.
 */

extern LUT_Definition_t refl_luts[];
extern LUT_Definition_t emiss_luts[];
extern LUT_Definition_t qa_luts[];


/*-----------------------------------------------------------------------
                       PGE VERSION MACRO
------------------------------------------------------------------------*/

#define PGE02_VERSION   "6.2.2_obpg"

/*-----------------------------------------------------------------------
                   DIMENSIONS AND INDEXING MACROS
           (See also "Granule.h" for other dimension macros)
------------------------------------------------------------------------*/

/*---------------------------------------------
   Common Lookup Table Dims
----------------------------------------------*/

#define MAX_SERIAL_NUMBER_BUFFER                  31
#define MAX_PGE_VERSION_BUFFER                    11
#define MAX_MCST_VERSION_BUFFER                   21
#define MAX_ASSOCIATEDPLATFORMSHORTNAME_BUFFER    21
#define MAX_ALGORITHMPACKAGEACCEPTANCEDATE_BUFFER 11
#define MAX_ALGORITHMPACKAGEMATURITYCODE_BUFFER   15
#define MAX_MISSION_PHASE_BUFFER                  15
#define NUM_DN_VALUES                 4096   /* 2^12 = 4096 (0..4095) */
#define NUM_PRI_RED_SYSTEMS                        2
#define NUM_BITS_IN_UINT8                          8

#define NUM_1ST_ORDER_COEFFS                       2
#define NUM_2ND_ORDER_COEFFS                       3
#define NUM_4TH_ORDER_COEFFS                       5
/*
 * The following is used specifically for the emissive and reflective 
 * 4th order polynomial evaluations as part of the uncertainty calculation.
 * The arrays dimensioned "NUM_4TH_ORDER_COEFFS" are the coefficients.
 */
#define EVAL_4TH_ORDER_POLYNOMIAL(p,a,x) \
  p = a[0] + x * (a[1] + x * (a[2] + x * (a[3] + x * a[4])));

/*
 * The following is used specifically for the emissive and reflective 
 * 1st order polynomial evaluations as part of the RVS calculation.
 * The arrays dimensioned "NUM_2ND_ORDER_COEFFS" are the coefficients.
 * Input "x" is the variable and input "y" is x*x, previously computed.
 */
#define EVAL_2ND_ORDER_POLYNOMIAL(p,a,x,y) \
  p = a[0] + x * a[1] + y * a[2];

/*---------------------------------------------
   SWIR-band (Reflective) Table Dims
----------------------------------------------*/

#define NUM_SWIR_BANDS                 4
#define MAX_NUM_SWIR_SUBSAMPLES        2
#define MAX_DETECTORS_PER_SWIR_BAND   20
#define MAX_SWIR_BAND_EV_FRAMES     2708
#define MAX_NUM_SWIR_RSR_WL           27     

/*---------------------------------------------
   Emissive Tables Dimensions and other macros

 Use caution when changing any of these, especially
 those related to polynomial evaluations.  For
 efficiency reasons, the code may not be general.
----------------------------------------------*/

#define NUM_PC_XT_BANDS            5
#define NUM_PC_XT_PARAMETERS       4
#define MAX_NUM_RSR_vs_LAMBDA     49

#define NUM_PVLW_XT_BANDS          4
#define NUM_PVLW_XT_COEFF          3

#define NUM_a0_vs_T_inst_COEFF     3  /* see caution above */
#define NUM_a2_vs_T_inst_COEFF     3  /* see caution above */
/* obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng */
/*
#define NUM_UI_PARAMETERS          8
#define NUM_UI_POLYNOMIAL_COEFF    2
#define NUM_FI_POLYNOMIAL_COEFF    5  
*/
#define NUM_T_MIR_THERMISTORS  2     /* number of values in LUT */
#define NUM_T_INS_THERMISTORS  4     /* number of values in LUT */
#define INDEX_TP_AO_SMIR_OBJ   0     /* LUT index */
#define INDEX_TP_AO_LWIR_OBJ   1     /* LUT index */
#define INDEX_TP_AO_SMIR_LENS  2     /* LUT index */
#define INDEX_TP_AO_LWIR_LENS  3     /* LUT index */

#define NUM_T_CAV_THERMISTORS  4     /* number of values in LUT */
#define INDEX_TP_MF_CALBKHD_SR 0     /* LUT index */
#define INDEX_TP_SR_SNOUT      1     /* LUT index */
#define INDEX_TP_MF_Z_BKHD_BB  2     /* LUT index */
#define INDEX_TP_MF_CVR_OP_SR  3     /* LUT index */

#define NUM_U2_FRAME           7     /* number of frames in uncertainty term u2 LUT */
#define NUM_RSB_RVS_COEFFS     5     /* number of coefficients used to calcuate RVS */

/*--------------------------------------------
   Dimensions for QA Tables
---------------------------------------------*/
#define NUM_MOON_OFFSET_LIMITS 4
#define NUM_ATTITUDE_ANGLES    3

/*
 * The following macros define the meaning of the index within the
 * "NUM_MOON_OFFSET_LIMITS" dimension of "moon_offset_limits"
 * in the common_QA_tables_t structure.
 */
#define TRK_UPPER 0  /* track upper limit index */
#define TRK_LOWER 1  /* track lower limit index */
#define SCN_UPPER 2  /* scan upper limit index */
#define SCN_LOWER 3  /* scan lower limit index */

/*-----------------------------------------------------------------------
                         LUT Information Macros
------------------------------------------------------------------------*/

/*
 * The following macros are used to determine the "kind" of LUT
 * as it exists inside the LUT file.
 */

#define GLOBAL_ATTRIBUTE_LUT  0      /* Must use read_attribute */
#define SDS_LUT               1      /* Must use read_sds_rankn */

/*---------------------------------------------
      Common Lookup Table names
(These tables should exist in all three files)
----------------------------------------------*/

#define PGE_VERSION_LUT_NAME             "PGE_Version"
#define MCST_VERSION_LUT_NAME            "MCST_Version"

/*---------------------------------------------
    Reflective Lookup Table names
----------------------------------------------*/

#define REFL_SERIAL_NUMBER_LUT_NAME        "Serial Number of Reflective LUT"
#define M0_LUT_NAME                        "m0"
#define M1_LUT_NAME                        "m1"
#define DN_STAR_MAX_LUT_NAME               "dn_star_Max"
#define DN_STAR_MIN_LUT_NAME               "dn_star_Min"
#define K_INST_LUT_NAME                    "K_inst"
#define K_FPA_LUT_NAME                     "K_FPA"
#define T_INST_REF_LUT_NAME                "T_inst_ref"
#define T_FPA_REF_LUT_NAME                 "T_FPA_ref"
#define RVS_RSB_LUT_NAME                   "RVS_RefSB"
/* obsolete due to RSB UI algorithm update, 2/19/2011, Xu Geng */
/*
#define SIGMA_RVS_RSB_LUT_NAME             "Sigma_RVS_RSB"
#define SIGMA_M1_LUT_NAME                  "Sigma_m1" 
#define SIGMA_K_INST_LUT_NAME              "Sigma_K_inst"
#define SIGMA_T_INST_LUT_NAME              "Sigma_T_inst"
#define SIGMA_PV_RESID_ELEC_LUT_NAME       "Sigma_PV_Resid_Elec"
#define SIGMA_R_STAR_LIN_RESID_UCOEFF_LUT_NAME  "Sigma_R_Star_Lin_Resid_Ucoeff"
#define RSB_NEDL_LUT_NAME                  "RSB_NEdL"
#define SIGMA_RSB_ADC_LUT_NAME             "Sigma_RSB_ADC"
*/
#define DN_OBC_1ST_FRAME_LUT_NAME          "DN_obc_avg_first_frame_to_use"
#define DN_OBC_NUM_FRAMES_LUT_NAME         "DN_obc_avg_number_of_frames_to_use"
#define SWIR_CORRECTION_SWITCH_LUT_NAME    "SWIR_OOB_correction_switch"
#define X_OOB_0_LUT_NAME                   "X_OOB_0"
#define X_OOB_1_LUT_NAME                   "X_OOB_1"
#define X_OOB_2_LUT_NAME                   "X_OOB_2"
#define RSB_SPECIFIED_UNCERTAINTY_LUT_NAME "RSB_specified_uncertainty"
#define RSB_UI_SCALING_FACTOR_LUT_NAME     "RSB_UI_scaling_factor"
#define E_SUN_OVER_PI_LUT_NAME             "E_sun_over_pi"
#define RSB_SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME "RSB_SV_DN_moon_include_frames"
#define DN_SAT_EV_LUT_NAME                 "dn_sat_ev"
#define B26_B5_CORR_LUT_NAME               "B26_B5_Corr"
#define B26_B5_CORR_SWITCH_LUT_NAME        "B26_B5_Corr_Switch"
#define B26_B5_FRAME_OFFSET_LUT_NAME       "B26_B5_Frame_Offset"
#define SWIR_CORR_SENDING_BAND_LUT_NAME    "SWIR_OOB_corr_sending_band"
#define SWIR_CORR_SENDING_DETECTOR_LUT_NAME "SWIR_OOB_corr_sending_detector"
/* new luts tables due to RSB uncertainty algorithm update 2/22/2011, Xu Geng */
#define U1_LUT_NAME                        "u1"
#define U2_LUT_NAME                        "u2"
#define U3_LUT_NAME                        "u3"
#define U4_LUT_NAME                        "u4"
#define U2_FRAMES_LUT_NAME                 "u2_frames"  
#define SWIR_UI_FACTOR_LUT_NAME            "swir_ui_factor"  

/*---------------------------------------------
    Emissive Lookup Table names
----------------------------------------------*/

#define EMISS_SERIAL_NUMBER_LUT_NAME       "Serial Number of Emissive LUT"
#define EPSILON_BB_LUT_NAME                "epsilon_bb"
#define EPSILON_CAV_LUT_NAME               "epsilon_cav"
#define DELTA_T_BB_BETA_LUT_NAME           "delta_T_bb_beta"
#define DELTA_T_BB_DELTA_LUT_NAME          "delta_T_bb_delta"
#define PCX_TALK_LUT_NAME                  "PC_XT"
#define RVS_TEB_LUT_NAME                   "RVS_TEB"
#define RVS_BB_SV_FRAME_NO_LUT_NAME        "RVS_BB_SV_Frame_No"
#define RSR_LUT_NAME                       "RSR"
#define WAVELENGTH_LUT_NAME                "WAVELENGTH"
#define NUM_WL_INCREMENT_LUT_NAME          "NWL"
#define CALIB_A0_LUT_NAME                  "A0"
#define CALIB_A2_LUT_NAME                  "A2"
/* obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng */
/*
#define UI_UCOEFF_LUT_NAME                 "Ucoeff"
#define SIGMA_TEB_PV_RESID_ELEC_LUT_NAME   "Sigma_TEB_PV_resid_elec"
#define SIGMA_TEB_ADC_LUT_NAME             "Sigma_TEB_ADC"
#define UCOEFF_CALIBR_RESID_LUT_NAME       "Ucoeff_Calibr_resid"
#define BAND_21_UNCERT_LSAT_LUT_NAME       "Band_21_Uncert_Lsat"
*/
#define T_INS_FUNCTION_FLAG_LUT_NAME       "T_ins_function_flag"
#define T_INS_DEFAULT_LUT_NAME             "T_ins_default"
#define T_INS_OFFSET_LUT_NAME              "T_ins_offset"
#define T_CAV_FUNCTION_FLAG_LUT_NAME       "T_cav_function_flag"
#define T_CAV_DEFAULT_LUT_NAME             "T_cav_default"
#define T_MIR_FUNCTION_FLAG_LUT_NAME       "T_mir_function_flag"
#define T_MIR_DEFAULT_LUT_NAME             "T_mir_default"
#define BB_WEIGHT_LUT_NAME                 "BB_Weight"
#define BB_DN_1ST_FRAME_LUT_NAME           "BB_DN_first_frame_to_use"
#define BB_DN_NUM_FRAMES_LUT_NAME          "BB_DN_number_of_frames_to_use"
#define SV_DN_1ST_FRAME_LUT_NAME           "SV_DN_first_frame_to_use"
#define SV_DN_NUM_FRAMES_LUT_NAME          "SV_DN_number_of_frames_to_use"
#define SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME "SV_DN_moon_include_frames"
#define PCX_CORRECTION_SWITCH_LUT_NAME     "PCX_correction_switch"
#define OVERLAP_SCANS_B1_LUT_NAME          "num_overlap_scans_b1"
#define BAND_21_B1_LUT_NAME                "Band_21_b1"
#define L_MIN_LUT_NAME                     "L_Min"
#define L_MAX_LUT_NAME                     "L_Max"
#define TEB_SPECIFIED_UNCERTAINTY_LUT_NAME "TEB_specified_uncertainty"
#define TEB_UI_SCALING_FACTOR_LUT_NAME     "TEB_UI_scaling_factor"
/* new luts tables due to TEB uncertainty algorithm update 3/22/2011, Xu Geng */
#define SIGMA_A0_LUT_NAME                  "sigma_a0"
#define SIGMA_A2_LUT_NAME                  "sigma_a2"
#define SIGMA_RVS_EV_LUT_NAME              "sigma_RVS_EV"
#define SIGMA_EPSILON_BB_LUT_NAME          "sigma_epsilon_BB"
#define SIGMA_EPSILON_CAV_LUT_NAME         "sigma_epsilon_CAV"
#define SIGMA_L_LAMBDA_LUT_NAME            "sigma_L_lambda"
#define SIGMA_L_TBB_LUT_NAME               "sigma_L_Tbb"
#define SIGMA_L_TSM_LUT_NAME               "sigma_L_Tsm"
#define SIGMA_L_TCAV_LUT_NAME              "sigma_L_Tcav"
#define SIGMA_B1_BAND21_LUT_NAME           "sigma_b1_B21"
#define XT_UI_FACTOR_LUT_NAME              "xt_ui_factor"
#define PVLW_XT_COEFF_LUT_NAME             "PVLW_XT_coeff"

/*---------------------------------------------
    QA Lookup Table information
----------------------------------------------*/

#define QA_SERIAL_NUMBER_LUT_NAME          "QA serial number"
#define PLATFORM_SHORT_NAME_LUT_NAME       "ASSOCIATEDPLATFORMSHORTNAME"
#define PACKAGE_ACCEPT_DATE_LUT_NAME       "ALGORITHMPACKAGEACCEPTANCEDATE"
#define PACKAGE_MATURITY_CODE_LUT_NAME     "ALGORITHMPACKAGEMATURITYCODE"
#define DET_QUAL_FLAG_VALS_LUT_NAME        "Detector Quality Flag Values"
#define DET_QUAL_FLAG2_VALS_LUT_NAME       "Detector Quality Flag2 Values"
#define MOON_OFFSET_LIMITS_LUT_NAME        "Moon Offset Limits"
#define MISSION_PHASE_LUT_NAME             "mission phase"
#define CONTROL_OPTIONS_LUT_NAME           "Control options"
#define BASE_VARI_VISUAL_FPA_LUT_NAME      "visual FPA base variance"
#define BASE_VARI_NIR_FPA_LUT_NAME         "NIR FPA base variance"
#define BB_TEMP_VARIANCE_LUT_NAME          "T_BB_Variance"
#define BB_AVG_TEMP_VAR_LUT_NAME           "BB Average Temperature Variance"
#define LWIR_FPA_TEMP_VAR_LUT_NAME         "LWIR FPA Temperature Variance"
#define MWIR_FPA_TEMP_VAR_LUT_NAME         "MWIR FPA Temperature Variance"
#define MIR_SIDE_1_TEMP_VAR_LUT_NAME       "MirrorSide 1 Temperature Variance"
#define MIR_SIDE_2_TEMP_VAR_LUT_NAME       "MirrorSide 2 Temperature Variance"
#define MIR_AVG_TEMP_VAR_LUT_NAME          "Mirror Average Temperature Variance"
#define INST_TEMP_VAR_LUT_NAME             "Instrument Temperature Variance"
#define CAVITY_TEMP_VAR_LUT_NAME           "Cavity Temperature Variance"
#define EMISS_NEdL_LUT_NAME                "NEdL"
#define CALIB_A1_LUT_NAME                  "a1"
#define ROLL_THRESHOLD_LUT_NAME            "Spacecraft_Roll_Threshold_Angle"
#define PITCH_THRESHOLD_LUT_NAME           "Spacecraft_Pitch_Threshold_Angle"
#define YAW_THRESHOLD_LUT_NAME             "Spacecraft_Yaw_Threshold_Angle"

/*-----------------------------------------------------------------------
                         LUT Structures
------------------------------------------------------------------------*/

/* 
 * SWIR Correction Lookup Table Structure
 */

typedef struct 
{
  int16   SWIR_correction_switch;
  int16   SWIR_corr_sending_band;
  int16   SWIR_corr_sending_detector[DETECTORS_PER_1KM_BAND];
  float32 X_OOB_0[NUM_SWIR_BANDS]
                 [MAX_DETECTORS_PER_SWIR_BAND]
                 [MAX_NUM_SWIR_SUBSAMPLES]
                 [NUM_MIRROR_SIDES];
  
  float32 X_OOB_1[NUM_SWIR_BANDS]
                 [MAX_DETECTORS_PER_SWIR_BAND]
                 [MAX_NUM_SWIR_SUBSAMPLES]
                 [NUM_MIRROR_SIDES];
  
  float32 X_OOB_2[NUM_SWIR_BANDS]
                 [MAX_DETECTORS_PER_SWIR_BAND]
                 [MAX_NUM_SWIR_SUBSAMPLES]
                 [NUM_MIRROR_SIDES];
                 
} SWIR_correction_tables_t;

/*
 * Reflective Lookup Table Structure
 */

typedef struct 
{
  char    Serial_Number[MAX_SERIAL_NUMBER_BUFFER];
  char    PGE_Version  [MAX_PGE_VERSION_BUFFER];
  char    MCST_Version [MAX_MCST_VERSION_BUFFER];

  int16   DN_obc_avg_first_frame_to_use;

  int16   DN_obc_avg_number_of_frames_to_use;

  float32 K_inst[NUM_REFLECTIVE_BANDS]
                [MAX_DETECTORS_PER_BAND]
                [MAX_SAMPLES_PER_BAND]
                [NUM_MIRROR_SIDES];

  float32 K_FPA[NUM_REFLECTIVE_BANDS]
               [MAX_DETECTORS_PER_BAND]
               [MAX_SAMPLES_PER_BAND]
               [NUM_MIRROR_SIDES];

  float32 m0[NUM_REFLECTIVE_BANDS]
            [MAX_DETECTORS_PER_BAND]
            [MAX_SAMPLES_PER_BAND]
            [NUM_MIRROR_SIDES];

  float32 m1[NUM_REFLECTIVE_BANDS]
            [MAX_DETECTORS_PER_BAND]
            [MAX_SAMPLES_PER_BAND]
            [NUM_MIRROR_SIDES];
 
  float32 dn_star_Max[NUM_REFLECTIVE_BANDS];

  float32 dn_star_Min[NUM_REFLECTIVE_BANDS];

  float32 T_inst_ref;

  float32 T_FPA_ref[NUM_FOCAL_PLANES];

  float32 RVS_RefSB[NUM_REFLECTIVE_BANDS]
                   [MAX_DETECTORS_PER_BAND]
                   [NUM_MIRROR_SIDES]
                   [NUM_RSB_RVS_COEFFS];

  /* obsolete due to RSB uncertainty algorithm update, 2/19/2011, Xu Geng */
  /*
  float32 Sigma_RVS_RSB[NUM_REFLECTIVE_BANDS]
                       [NUM_MIRROR_SIDES];

  float32 Sigma_m1[NUM_REFLECTIVE_BANDS]
                  [MAX_DETECTORS_PER_BAND]
                  [MAX_SAMPLES_PER_BAND]
                  [NUM_MIRROR_SIDES];

  float32 Sigma_K_inst[NUM_REFLECTIVE_BANDS]
                      [MAX_DETECTORS_PER_BAND]
                      [MAX_SAMPLES_PER_BAND]
                      [NUM_MIRROR_SIDES];

  float32 Sigma_T_inst;

  float32 Sigma_PV_Resid_Elec[NUM_REFLECTIVE_BANDS]
                             [MAX_DETECTORS_PER_BAND]
                             [MAX_SAMPLES_PER_BAND];

  float32 Sigma_R_Star_Lin_Resid_Ucoeff[NUM_REFLECTIVE_BANDS]
                                       [MAX_DETECTORS_PER_BAND]
                                       [MAX_SAMPLES_PER_BAND]
                                       [NUM_MIRROR_SIDES]
                                       [NUM_4TH_ORDER_COEFFS];

  float32 RSB_NEdL[NUM_REFLECTIVE_BANDS]
                  [MAX_DETECTORS_PER_BAND]
                  [MAX_SAMPLES_PER_BAND]
                  [NUM_MIRROR_SIDES];

  float32 Sigma_RSB_ADC[NUM_REFLECTIVE_BANDS]
                       [MAX_DETECTORS_PER_BAND];
  */

  float32 RSB_specified_uncertainty[NUM_REFLECTIVE_BANDS];

  float32 RSB_UI_scaling_factor[NUM_REFLECTIVE_BANDS];

  SWIR_correction_tables_t SWIR_correction_tables; /* All SWIR LUTs */
  float32 E_sun_over_pi[NUM_REFLECTIVE_DETECTORS];
  int16   RSB_SV_DN_moon_include_frames;

  float64 dn_sat_ev[NUM_REFLECTIVE_BANDS]
                   [MAX_DETECTORS_PER_BAND]
                   [MAX_SAMPLES_PER_BAND]
                   [NUM_MIRROR_SIDES];

  float32 B26_B5_Corr[DETECTORS_PER_1KM_BAND];
  int16   B26_B5_Corr_Switch;
  int16   B26_B5_Frame_Offset[DETECTORS_PER_1KM_BAND];

  /* new luts tables due to RSB uncertainty algorithm update 2/22/2011, Xu Geng */
  float32 u1[NUM_REFLECTIVE_DETECTORS];
  float32 u2_samples[NUM_REFLECTIVE_DETECTORS]
                    [NUM_MIRROR_SIDES]
                    [NUM_U2_FRAME];
  float32 u3[NUM_REFLECTIVE_DETECTORS]
            [NUM_MIRROR_SIDES];
  float32 u4_coeffs[NUM_REFLECTIVE_BANDS]
                   [MAX_DETECTORS_PER_BAND]
                   [MAX_SAMPLES_PER_BAND]
                   [NUM_MIRROR_SIDES]
                   [NUM_2ND_ORDER_COEFFS];
  int16   u2_frames[NUM_U2_FRAME]; 
  float32 swir_ui_factor[NUM_SWIR_BANDS];
 
  } refl_tables_t;


/*
 * Emissive Lookup Table Structure
 * NOTE: PC_XT values in the lookup table are in percent.
 */

typedef struct 
{
  char     Serial_Number   [MAX_SERIAL_NUMBER_BUFFER];
  char     PGE_Version     [MAX_PGE_VERSION_BUFFER];
  char     MCST_Version    [MAX_MCST_VERSION_BUFFER];
  float32  epsilon_bb      [NUM_EMISSIVE_DETECTORS];
  float32  epsilon_cav     [NUM_EMISSIVE_DETECTORS];
  float32  delta_T_bb_beta [NUM_EMISSIVE_DETECTORS];
  float32  delta_T_bb_delta[NUM_EMISSIVE_DETECTORS];
  float32  PC_XT           [NUM_PC_XT_BANDS]
                           [DETECTORS_PER_1KM_BAND]
                           [NUM_PC_XT_PARAMETERS];
  float32  RSR             [NUM_EMISSIVE_DETECTORS]
                           [MAX_NUM_RSR_vs_LAMBDA];
  float32  wavelength      [NUM_EMISSIVE_DETECTORS]
                           [MAX_NUM_RSR_vs_LAMBDA];
  float32  A0              [NUM_a0_vs_T_inst_COEFF]
                           [NUM_MIRROR_SIDES]
                           [NUM_EMISSIVE_DETECTORS];
  float32  A2              [NUM_a2_vs_T_inst_COEFF]
                           [NUM_MIRROR_SIDES]
                           [NUM_EMISSIVE_DETECTORS];
 
  int16    NUM_RSR_vs_Lambda[NUM_EMISSIVE_DETECTORS]; 
  float32  CW               [NUM_EMISSIVE_BANDS];

  /* obsolete due to TEB uncertainty algorithm update, 3/22/2011, Xu Geng */
  /*
  float32  Ucoeff      [NUM_EMISSIVE_DETECTORS]
                       [NUM_UI_PARAMETERS]
                       [NUM_UI_POLYNOMIAL_COEFF]
                       [NUM_FI_POLYNOMIAL_COEFF];
  float32  Sigma_TEB_PV_resid_elec[NUM_EMISSIVE_DETECTORS];
  float32  Sigma_TEB_ADC[NUM_EMISSIVE_DETECTORS];
  float32  Ucoeff_Calibr_resid[NUM_EMISSIVE_DETECTORS]
                              [NUM_4TH_ORDER_COEFFS];
  float32 Band_21_Uncert_Lsat;
  */
  int16 SV_DN_first_frame_to_use;
  int16 SV_DN_number_of_frames_to_use;
  int16 SV_DN_moon_include_frames;
  int16 BB_DN_first_frame_to_use;
  int16 BB_DN_number_of_frames_to_use;
  int16 num_overlap_scans_b1;
  int8  PCX_correction_switch;
  int32   T_ins_function_flag[NUM_T_INS_THERMISTORS];
  float32 T_ins_default;
  float32 T_ins_offset[NUM_T_INS_THERMISTORS];
  int32   T_cav_function_flag[NUM_T_CAV_THERMISTORS];
  int32   T_mir_function_flag[NUM_T_MIR_THERMISTORS];
  float32 T_cav_default;
  float32 T_mir_default;
  float32 BB_Weight[NUM_BB_THERMISTORS];

  float32 RVS_TEB[NUM_EMISSIVE_BANDS]
                 [DETECTORS_PER_1KM_BAND]
                 [NUM_MIRROR_SIDES]
                 [NUM_2ND_ORDER_COEFFS];
  int16   RVS_BB_SV_Frame_No[2]; 
  
  float32 Band_21_b1[DETECTORS_PER_1KM_BAND][NUM_MIRROR_SIDES];
  float32 L_Max[NUM_EMISSIVE_BANDS];
  float32 L_Min[NUM_EMISSIVE_BANDS];
  float32 TEB_specified_uncertainty[NUM_EMISSIVE_BANDS];
  float32 TEB_UI_scaling_factor[NUM_EMISSIVE_BANDS];

  /* new luts tables due to TEB uncertainty algorithm update 3/22/2011, Xu Geng */ 
  float32 sigma_a0[NUM_a0_vs_T_inst_COEFF]
                  [NUM_MIRROR_SIDES]
                  [NUM_EMISSIVE_DETECTORS]; 
  float32 sigma_a2[NUM_a2_vs_T_inst_COEFF]
                  [NUM_MIRROR_SIDES]
                  [NUM_EMISSIVE_DETECTORS]; 
  float32 sigma_RVS_EV[NUM_EMISSIVE_BANDS]
                      [DETECTORS_PER_1KM_BAND]
                      [NUM_MIRROR_SIDES]
                      [NUM_2ND_ORDER_COEFFS];
  float32 sigma_epsilon_BB[NUM_EMISSIVE_BANDS];
  float32 sigma_epsilon_CAV[NUM_EMISSIVE_BANDS]; 
  float32 sigma_L_lambda[NUM_EMISSIVE_BANDS]
                        [NUM_1ST_ORDER_COEFFS]; 
  float32 sigma_L_Tbb[NUM_EMISSIVE_BANDS];
  float32 sigma_L_Tsm[NUM_EMISSIVE_BANDS];
  float32 sigma_L_Tcav[NUM_EMISSIVE_BANDS];
  float32 sigma_b1_B21[DETECTORS_PER_1KM_BAND][NUM_MIRROR_SIDES];
  float32 xt_ui_factor[NUM_EMISSIVE_BANDS][DETECTORS_PER_1KM_BAND];
  float64 PVLW_XT_coeff    [NUM_PVLW_XT_BANDS]        /*receiving bands*/
                           [DETECTORS_PER_1KM_BAND]   /*detector in receiving bands*/ 
                           [NUM_PVLW_XT_BANDS]        /*sending bands*/
                           [DETECTORS_PER_1KM_BAND]   /*detector in sending bands*/ 
                           [NUM_PVLW_XT_COEFF];

} emiss_tables_t;

/*
 * common QA tables 
 */

/*
 * Control options within L1B:
 * 1) Split scans -- control whether or not to treat as missing if a split
 *    scan is detected.  There are only two values, ON or OFF.
 * 2) Bad scan quality -- control whether or not to treat as missing if
 *    an invalid value of Scan quality array is detected.
 *    There are only two values, ON or OFF.
 */

typedef enum {
  SPLIT_SCAN_CONTROL,
  BAD_SCAN_QUALITY_CONTROL,
  NUM_CONTROL_OPTIONS
} control_options_defs_t;

typedef struct 
{
  char    Serial_Number[MAX_SERIAL_NUMBER_BUFFER];
  char    PGE_Version[MAX_PGE_VERSION_BUFFER];
  char    MCST_Version[MAX_MCST_VERSION_BUFFER];
  char    AssociatedPlatformShortname
            [MAX_ASSOCIATEDPLATFORMSHORTNAME_BUFFER];
  char    AlgorithmPackageAcceptanceDate
            [MAX_ALGORITHMPACKAGEACCEPTANCEDATE_BUFFER];	 
  char    AlgorithmPackageMaturityCode
            [MAX_ALGORITHMPACKAGEMATURITYCODE_BUFFER];
  int8    dead_detector[NUM_DETECTORS];
  int8    noisy_detector[NUM_DETECTORS];
  int8    dead_subframe[NUM_HIGH_RESOLUTION_SUBFRAMES];
  int8    noisy_subframe[NUM_HIGH_RESOLUTION_SUBFRAMES];
  uint8   Detector_Quality_Flag_Values[NUM_DETECTORS]
                                      [NUM_BITS_IN_UINT8];
  uint8   Detector_Quality_Flag[NUM_DETECTORS];
  uint8   Detector_Quality_Flag2_Values[NUM_HIGH_RESOLUTION_DETECTORS]
                                       [NUM_BITS_IN_UINT8];
  uint8   Detector_Quality_Flag2[NUM_HIGH_RESOLUTION_DETECTORS];
  float32 moon_offset_limits[NUM_BANDS]
                            [NUM_MOON_OFFSET_LIMITS];
  char    mission_phase[MAX_MISSION_PHASE_BUFFER];
  uint8   control_options[NUM_CONTROL_OPTIONS];
  float32 roll_threshold_angle;
  float32 pitch_threshold_angle;
  float32 yaw_threshold_angle;
} common_QA_tables_t;

/*
 * refl QA tables 
 */

typedef struct 
{
  float32  var_visual_FPA;
  float32  var_NIR_FPA;
} refl_QA_tables_t;

/*
 * emiss QA tables 
 */

typedef struct 
{
  float32 var_T_bb[NUM_BB_THERMISTORS];
  float32 var_T_bb_avg;
  float32 var_T_lwir;
  float32 var_T_mwir;
  float32 var_T_mir1;
  float32 var_T_mir2;
  float32 var_T_mir_avg;
  float32 var_T_ins;
  float32 var_T_cav;
  float32 NEdL[NUM_EMISSIVE_DETECTORS];
  float32 a1  [NUM_EMISSIVE_DETECTORS];
} emiss_QA_tables_t;

/*
 * QA Lookup table structure
 */

typedef struct 
{
  common_QA_tables_t common_QA_tables;
  refl_QA_tables_t   refl_QA_tables;
  emiss_QA_tables_t  emiss_QA_tables;
} QA_tables_t;

/*
 * L1B lookup tables structure
 */

typedef struct 
{
  refl_tables_t  refl;
  emiss_tables_t emiss;
  QA_tables_t    QA;
} lookup_tables_t;

/*--------------------------------------------------------
              Function Prototypes 
---------------------------------------------------------*/

PGSt_SMF_status Read_Lookup_Tables     
                      (L1A_granule_t          *L1A_Gran,
                       lookup_tables_t        *tables,
                       Run_Time_Parameters_t  *runtime_params);

PGSt_SMF_status Read_Refl_Tables       
                      (L1A_granule_t          *,
                       refl_tables_t          *);

PGSt_SMF_status Read_QA_Tables
                      (L1A_granule_t          *,
                       QA_tables_t            *);

PGSt_SMF_status Read_Emiss_Tables
                      (L1A_granule_t          *,
                       emiss_tables_t         *);

PGSt_SMF_status Read_L1B_SDS_LUT
                      (int32                  sd_id,
                       char                   *name,
                       int32                  data_type,
                       int32                  rank,
                       int32                  *dims,
                       float64                data_collection_TAI_time,
                       void                   *data);

int32 TDLUT_GetAlgorithm
                      (int32                  sd_id,
                       char                   *name);

PGSt_SMF_status TDLUT_ReadStepFunction
                      (int32                  sd_id,
                       char                   *name,
                       int32                  data_type,
                       int32                  rank,
                       int32                  *dims,
                       float64                data_collection_TAI_time,
                       void                   *data);

PGSt_SMF_status TDLUT_ReadPiecewiseLinearFunction
                      (int32                  sd_id,
                       char                   *name,
                       int32                  data_type,
                       int32                  rank,
                       int32                  *dims,
                       float64                data_collection_TAI_time,
                       void                   *data);

PGSt_SMF_status Expand_BDSM_LUT
                      (void                   *data,
                       void                   **data_new,
                       int32                  data_type,
                       int32                  lead_dim,
                       int32                  *n_bytes);

int32 BDSM_index      (char                   *ascii_file);

#endif


