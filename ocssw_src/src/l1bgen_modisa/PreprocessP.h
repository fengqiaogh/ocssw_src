#ifndef	PREPROCESSP_H
#define	PREPROCESSP_H

#include    "Preprocess.h"
/*
!C-INC********************************************************************** 
!Description:      Private header file to be used only in Preprocess.c. 
                   Contains variable and macro definition, function prototypes.  

!Revision History:
 $Log: PreprocessP.h,v $
 Revision 1.18  2015-03-17 10:09:24-04  xgeng
 Fix the anomaly before and after the sector rotation due to the timing mismatch between the command performed and the status updated.

 Revision 1.17  2011-07-28 15:37:21-04  xgeng
 Change the function prototype to pass LWIR FPA temperature to Get_Emiss_Coeff_Per_Scan

 Revision 1.15  2011-04-07 14:40:48-04  xgeng
 1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

 Revision 1.14  2010-11-15 11:22:05-05  xgeng
 Defined a new constant ELECTRONICS_BOTH

 Revision 1.13  2009/07/24 19:46:41  xgeng
 Added Last_Valid_Scans member to structure L1A_granule_OBCEng_t and Overlap_OBCEng_t to fix the sector-rotation abnormal

 Revision 1.11  2008/01/24 15:50:37  ltan
 Changed to relax the RVS correction limit range from [0.8, 1.2] to [0.4, 2.4].

 Revision 1.10  2006/10/27 15:01:36  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 01.15   April 16, 2003
 Changed the nominal platform height from 750.259 to its correct value, 705.295. The
 platform height is used only to determine the pixel offset when the Moon in SV
 keep-out-box (KOB) algorithm is used.  The incorrect value resulted in an actual SV  KOB
 slightly smaller than intended (approximately 3 pixels in the along- scan direction and
 .9 pixels in the along-track direction).  Since the actual  desired boundaries of the SV
 KOB are at +/-15 pixels along track and +/-55 pixels along scan and the pixels where SV
 data are taken are +/-5 along track and +/-25 along scan, the impact of the error was
 minimal.
 Alice Isaacman, SAIC GSO    (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 01.14, March 27, 2003    Razor issue #173
 In the initialization of the Engineering_Coefficients_t structures 
 Engineering_Coefficients_PFM and Engineering_Coefficients_FM1, enclosed the
 rows of the 2-dimensional arrays C_FP3, C_FP4, C_RC, and BB_a with braces
 for ANSI-C compliance.
 Liqin Tan,  SAIC GSO (ltan@saicmodis.com)

 Revision 01.13, April 15, 2002    Razor issue #166
 Removed GRANULE_TIME_INTERVAL and GRANULE_TIMEDIFF_TOLERANCE
 Added INNER_GRANULE_TIMEDIFF_TOLERANCE and OUTER_GRANULE_TIMEDIFF_TOLERANCE.
 Change Read_Overlap_OBCEng prototype to allow passing of boolean pointer
   "scan_gap".
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.12  March 25, 2002  Razor Issue #178
 Strip out ADC Correction
 Alice_Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.11, March 11, 2002  Razor Issue #174
 Change function prototypes to pass 
   Emissive SV and BB RVS corrections.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.10, January 29, 2002  Razor Issue #175
 Change function prototypes as appropriate to pass satellite ID.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.09, Dec. 28, 2001, Razor issue #154
 Added temperature coefficient structure and made Terra- and Aqua- specific
 coefficient structures.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.08, Nov. 5, 2001, Razor issue 166
 Added GRANULE_TIME_INTERVAL and GRANULE_TIMEDIFF_TOLERANCE.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.07, Dec 7, 2000, Razor issue 146.
 Remove Fill_DN_Sat_Prime, add MirrorSide to Fill_xxx_DN_OBC_Avg.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.06  October 24, 2000
 Added code for granule average QA values as per issue 140.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.05 Sep 5, 2000
 Made several changes for Razor issue 130 (errors computing "b1" when
 near a sector rotation or if Ecal ends in the last 40 scans of the leading
 granule or starts in the 1st 40 scans of the trailing granule).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.04 November 23, 1999
 Read_LeadingOverlap and Read_TrailingOverlap replaced by Read_Overlap_OBCEng.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.03 September 30, 1999
 Added num_thermistor_outliers to Temperature_t structure.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.02 May 6, 1999
 Removed ADC_correction, Best_Pixel_Avg_int16(), Calculate_SV_Best_Pixel()
 and Get_L1B_SD_SRCA_BB_SV(). 
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb 17, 1999
 Moved ADC_correction into this module from Reflective_Cal as part of
 DN_to_DN_star to avoid circular include dependencies.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 Nov. 1998
 part of original Preprocess.h into this header file
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)
 
!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END***********************************************************************
*/


/* numbers of different kinds of bands used as dimensionality */

#define MAX_NUM_BANDS_PER_TYPE     11
#define NUM_BAND_TYPES             5
#define NUM_PV_VIS_BANDS           7
#define NUM_PV_NIR_BANDS           11
#define NUM_SM_BANDS               10
#define NUM_PV_LW_BANDS            4
#define NUM_PC_LW_BANDS            6

#define ELECTRONICS_PRIMARY        0  /* Electronics side "A" */
#define ELECTRONICS_REDUNDANT      1  /* Electronics side "B" */
#define ELECTRONICS_BOTH           2  /* both side of Electronics on */
#define MAX_NUM_OVERLAP_SCANS      MAX_NUM_SCANS
#define MAX_OVERLAP_TRACK_DIM      MAX_1KM_TRACK_DIM
#define MAX_TOTAL_XGRAN_SCANS      (2*MAX_NUM_OVERLAP_SCANS + MAX_NUM_SCANS)
#define GEO_NUM_SCANS_NAME         "Number of Scans"
#define GEO_SD_SUN_AZ_SDS_NAME     "SD Sun azimuth"
#define GEO_SD_SUN_ZEN_SDS_NAME    "SD Sun zenith"
#define OBC_SD_SUN_AZ_SDS_NAME     "SD Sun Azimuth"
#define OBC_SD_SUN_ZEN_SDS_NAME    "SD Sun Zenith" 
#define DN_BG_FIRST_FRAME_SDS_NAME "DN_obc_avg_first_frame_to_use"
#define DN_BG_NUM_FRAMES_SDS_NAME  "DN_obc_avg_number_of_frames_to_use"
#define VEC_DIM_NAME               "vecdim"
#define NUM_BANDS_DIM_NAME         "num_bands"
#define MOON_IN_KOB_SDS_NAME       "Moon in keep-out-box"
#define VEC_DIM                    3
#define MOON_VECTOR_SDS_NAME       "Moon Vector"
#define MAX_ENG_MEM_DATA_DIM       550
#define MAX_NUM_VDATA              100
#define MAX_VDATA_SIZE             1024
#define NOMINAL_PIXEL_WIDTH        1.0
#define NOMINAL_PLATFORM_HEIGHT    705.295
#define PIXEL_WIDTH_RADIANS        (NOMINAL_PIXEL_WIDTH/NOMINAL_PLATFORM_HEIGHT)
#define SV_CNTR_ANG_FROM_Y         (-8.425 * (PGS_PI/180.0))
#define SC_ANCIL_SIZE              128
#define MOON_X_DATA_POSITION_IN_SC_ANCIL 61
#define CONVERT_TO_UINT8(x,y) \
  if ( x < 0 )      y = 0;    \
  else if ( x >= 9) y = 255;  \
  else  y = (uint8) (255.0 * (x / 9.0));

#define  Celsius_To_Kelvin_Offset  273.15
#define DN_SAT_9BITS   511  /* saturation level for 9-bit temperature DN */
#define DN_SAT_10BITS 1023  /* saturation level for 10-bit voltage DN */
#define DN_SAT_12BITS 4095  /* saturation level for 12-bit temperature DN */

  /*
   * Emissive L_vs_DN Algorithm Planck function MACRO
   * (needs two constants, C_1 and C_2,
   * for radiance in the units of W/m^2/micron/sr)
   */
#define  C_1  1.19106e+8
#define  C_2  1.43879e+4
#define  L(lambda,T) (C_1 / pow((double) (lambda), (double) 5) / \
                      (exp((double) (C_2 / ((lambda) * (T))) ) - 1) )

   /* Upper and Lower limits on calculated RVS Correction terms */
#define RVS_CORRECTION_UPPER_LIMIT  2.4e+0
#define RVS_CORRECTION_LOWER_LIMIT  4.0e-1  

   /* Time elapsed between MODIS scans in seconds */
#define SCAN_TIME_INTERVAL  1.47717

  /* Highest acceptable number of dropped scans between granules */
#define DROPPED_SCAN_UPPER_LIMIT 5

   /* Tolerance of difference between estimated and actual time differences.  */
#define INNER_GRANULE_TIMEDIFF_TOLERANCE  .5*SCAN_TIME_INTERVAL
#define OUTER_GRANULE_TIMEDIFF_TOLERANCE                                   \
         (DROPPED_SCAN_UPPER_LIMIT + .5)*SCAN_TIME_INTERVAL

/* satellite ID */

#define TERRA 0
#define AQUA  1
#define INVALID_SATELLITE_ID -1

   /* Bands for MODIS/Aqua (FM1) which show dn BB saturation on BB warmup: 
      (indices within Emissive Bands): */
static int16 bb_sat_bands_indices_aqua[NUM_AQUA_BB_SAT_BANDS] = {12, 14, 15};

  /*  Structure which holds engineering coefficients */
  typedef struct {
    int32     Rp       ;
    float32   V0_CPA   ;
    float32   I_CPA    ;
    float32   V0_CPB   ;
    float32   I_CPB    ;   
      
    float32 C_FP1[6];  
    float32 C_FP2[6];  
    float32 C_FP3[3][6];   
    float32 C_FP4[3][6];   
    float32 C_MIR_A[6];   
    float32 C_MIR_B[6];   
    float32 C_CAV_A[6];    
    float32 C_CAV_B[6];    
    float32 C_ADC[6];   
    float32 C_INS_A[6];
    float32 C_INS_B[6];
    float32 C_RC[NUM_T_RC_VALUES][6];  
    float32 C_RC_VR_FPA[6];   
    float32 BB_a[NUM_BB_THERMISTORS][4];   
     } Engineering_Coefficients_t;
     
  /*   ENGINEERING COEFFICIENTS FOR TERRA (PFM) 
   *
   * -------------------------------------------------------
   *   DN to Eng Units equation coefficients
   *   Reference: SBRS Doc. 151840 REV C, T20-5A and T20-6A
   * -------------------------------------------------------
   */ 

   Engineering_Coefficients_t Engineering_Coefficients_PFM = {
      2.0000E+07,         /* Rp     */
     -3.216341,           /* V0_CPA */
      0.000332142,        /* I_CPA */ 
     -3.215385,           /* V0_CPB */
      0.000332376,        /* I_CPB */ 

    /*
     * T_fp1 --> TA_AO_VIS_FPA
     * T_fp2 --> TA_AO_NIR_FPA
     * unit = ^oC
     * Ref. T20-5A    SBRS Doc. 151840 REV C
     */ 
    { /* TERRA (PFM) C_FP1[6] */
      -5.8957e+01, 3.5611e-02, 2.7192e-06, -1.2548e-09, 1.4924e-13, 0. },   
    { /* TERRA (PFM) C_FP2[6]*/
      -5.7755e+01, 3.5526e-02, 2.8295e-06, -1.2865e-09, 1.5373e-13, 0. },   

    /*
     * T_fp3 --> TA_RC_SMIR_CFPA
     * T_fp4 --> TA_RC_LWIR_CFPA
     * unit = K
     * Ref. T20-5A    SBRS Doc. 151840 REV C
     */ 
    { /*[FP_SPSs][6]*/  /*TA_RC_SMIR_CFPA*/   /* TERRA (PFM) C_FP3[3][6] */
     {4.9970e+01, 1.8360e-02, -9.1174e-07, 1.7673e-10, -1.3356e-14, 0.},
     {4.9912e+01, 1.8426e-02, -9.6036e-07, 1.9249e-10, -1.5208e-14, 0.},
     {4.9846e+01, 1.8413e-02, -9.4652e-07, 1.8741e-10, -1.4575e-14, 0.}
    }, 
      
    { /*[FP_SPSs][6]*/  /*TA_RC_LWIR_CFPA*/   /* TERRA (PFM) C_FP4[3][6] */
     {8.0062e+01, 1.4493e-03, -6.8684e-09, 1.9367e-12, -2.5410e-16, 0.},
     {8.2050e+01, 1.4442e-03, -5.6085e-09, 1.0498e-12, -1.0151e-16, 0.},
     {8.5048e+01, 1.4388e-03, -8.8412e-09, 2.5727e-12, -2.9433e-16, 0.}
    }, 
          
    /*
     * T_mir --> TP_SA_RCT1_MIR & TP_SA_RCT2_MIR, 
     * they use the same coefficients
     * unit = ^oC 
     * Ref. T20-6A    SBRS Doc. 151840 REV C
     */             

    /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    { /* TERRA (PFM) C_MIR_A[6] */
      8.6953e+01, -7.0286e-02, 4.1191e-05, -1.6989e-08, 
                  3.6372e-12, -3.2037e-16 },    
    { /* TERRA (PFM) C_MIR_B[6] */
      8.6953e+01, -7.0286e-02, 4.1191e-05, -1.6989e-08, 
                  3.6372e-12, -3.2037e-16 },    

    /*
     *   T_cav --> TP_MF_CALBKHD_SR, TP_SA_A_MTR
     *   unit = ^oC
     *   Ref. T20-6
     */
 
    /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    {  /* TERRA (PFM) C_CAV_A[6] */
      8.6954e+01, -5.6241e-01, 2.6375e-03, -8.7030e-06, 
              1.4904e-08, -1.0501e-11 },   
    {  /* TERRA (PFM) C_CAV_B[6] */
      8.6954e+01, -5.6241e-01, 2.6375e-03, -8.7030e-06, 
              1.4904e-08, -1.0501e-11 },    
    
    /*
     *  T_ADC --> TA_PVLW_PWB4_10
     *              TA_PVSM_PWB6_12
     *   unit = ^oC
     *   Ref. T20-5A    SBRS Doc. 151840 REV C
     */
     {     /* TERRA (PFM) C_ADC[6] */
       -1.8888e+01, 4.9724e-01, -2.3834e-03, 8.1515e-06, 
                   -1.3848e-08, 9.2793e-12},  
    /*
     *   T_INS --> TP_AO_SMIR_OBJ, TP_AO_SMIR_LENS, 
     *             TP_AO_LWIR_OBJ, TP_AO_LWIR_LENS
     *   unit = ^oC
     *   Ref. T20-6A    SBRS Doc. 151840 REV C
     */ 

    /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    {  /* TERRA (PFM) C_INS_A[6] */
      8.6953e+01, -7.0286e-02, 4.1191e-05, -1.6989e-08, 
              3.6372e-12, -3.2037e-16},   
    {  /* TERRA (PFM) C_INS_B[6] */
      8.6953e+01, -7.0286e-02, 4.1191e-05, -1.6989e-08, 
              3.6372e-12, -3.2037e-16},  
    
    /* RC temperature conversion
     * num bits = 12 (max DN value = 4095)
     * units = Kelvin
     * ref. T20-5A    SBRS Doc. 151840 REV C
     */
    {  /* TERRA (PFM) C_RC[NUM_T_RC_VALUES][6] */
     {-2.6607e+00, 2.9239e-02, -1.2789e-06, 1.4882e-10, 0., 0.}, /* TA_RC_CS */
     {-1.3726e+02, 7.5134e-02, -2.7944e-06, 3.5353e-09, 0., 0.}, /* TA_RC_CS_OG */
     { 3.1238e+01, 2.6812e-02,  1.8002e-06, 2.8576e-11, 0., 0.}, /* TA_RC_IS */
     {-1.4316e+02, 8.0885e-02, -4.5854e-06, 3.7057e-09, 0., 0.}, /* TA_RC_IS_OG */
     {-1.4102e+02, 7.9490e-02, -4.3164e-06, 3.6498e-09, 0., 0.}  /* TA_RC_OS_OG */
    }, 
      
    /* RC FPA voltage heater conversion (LWIR or SMIR)
     * num bits = 10 (max DN value = 1023)
     * units = volts DC
     *    VR_RC_LW_FPA_HTR
     *    VR_RC_SM_FPA_HTR
     * ref T20-5    SBRS Doc. 151840 REV C
     */
    { /* TERRA (PFM) C_RC_VR_FPA[6] */
      -1.4970e+01, 2.9238e-02, 0., 0., 0., 0. }, 

    /*
     * All engineering coefficients for Terra OBC BB are described in 
     * "Preliminary Investigation of the On-Board Calibrator (OBC) Blackbody
     * Temperature Calibration Algorithm" by J. Xiong and T. Dorman,
     * 03/24/1998.
     */
    { /* TERRA (PFM) BB_a[NUM_BB_THERMISTORS][4] */
     {9.298819E-04,	3.028024E-04,	-4.973001E-06,	3.284554E-07},
     {9.521969E-04,	2.972846E-04,	-4.348669E-06,	3.039369E-07},
     {9.614388E-04,	2.905839E-04,	-3.482982E-06,	2.769668E-07},
     {9.533260E-04,	2.949576E-04,	-4.035258E-06,	2.947707E-07},
     {9.946773E-04,	2.814390E-04,	-2.494828E-06,	2.345659E-07},
     {9.808668E-04,	2.872997E-04,	-3.192201E-06,	2.592090E-07},
     {9.699586E-04,	2.888108E-04,	-3.350490E-06,	2.674701E-07},
     {9.550755E-04,	2.944409E-04,	-3.996653E-06,	2.917877E-07},
     {9.491089E-04,	2.977076E-04,	-4.445135E-06,	3.096761E-07},
     {1.011649E-03,	2.766499E-04,	-1.983597E-06,	2.138617E-07},
     {1.028646E-03,	2.700306E-04,	-1.171597E-06,	1.830096E-07},
     {1.001408E-03,	2.788388E-04,	-2.210990E-06,	2.223562E-07}         
    }
   };


  /*   ENGINEERING COEFFICIENTS FOR AQUA (FM1) 
   *
   * -------------------------------------------------------
   *   DN to Eng Units equation coefficients
   *   Reference: SBRS Doc. 151840 REV C, T20-5B and T20-6B 
   * -------------------------------------------------------
   */ 
          
   Engineering_Coefficients_t Engineering_Coefficients_FM1 = {
        2.0000E+07 ,     /* Rp     */
       -3.225313   ,     /* V0_CPA */
        0.000333193,     /* I_CPA */
       -3.228853   ,     /* V0_CPB */
        0.000334205,     /* I_CPB */
      
    /*
     * T_fp1 --> TA_AO_VIS_FPA
     * T_fp2 --> TA_AO_NIR_FPA
     * unit = ^oC
     * Ref. T20-5B    SBRS Doc. 151840 REV C
     */ 
    { /* AQUA (FM1) C_FP1[6] */ 
      -5.8639e+01, 3.6188e-02, 2.6495e-06, -1.3749e-09, 1.7933e-13, 0. }, 
    { /* AQUA (FM1) C_FP2[6] */ 
      -5.8907e+01, 3.5310e-02, 2.8666e-06, -1.3021e-09, 1.5744e-13, 0. }, 
    
    /*
     * T_fp3 --> TA_RC_SMIR_CFPA
     * T_fp4 --> TA_RC_LWIR_CFPA
     * unit = K
     * Ref. T20-5B    SBRS Doc. 151840 REV C
     */ 
    { /*[FP_SPSs][6]*/  /*TA_RC_SMIR_CFPA*/  /* AQUA (FM1) C_FP3[3][6] */
     {4.9904e+01, 1.8157e-02, -6.8273e-07, 1.0882e-10, -7.0948e-15, 0.},
     {4.9832e+01, 1.8200e-02, -7.0521e-07, 1.1395e-10, -7.5264e-15, 0.},
     {4.9762e+01, 1.8208e-02, -7.0835e-07, 1.1460e-10, -7.5779e-15, 0 }
    }, 
      
    { /*[FP_SPSs][6]*/  /*TA_RC_LWIR_CFPA*/  /* AQUA (FM1) C_FP4[3][6] */
     {8.0049e+01, 1.4429e-03, -2.0639e-09, 2.0749e-14, 0., 0.},
     {8.2079e+01, 1.4374e-03, -2.0972e-09, 3.9308e-14, 0., 0.},
     {8.5062e+01, 1.4287e-03, -1.8224e-09, 9.6394e-14, 0., 0}
    }, 
      
    /*
     * T_mir --> TP_SA_RCT1_MIR & TP_SA_RCT2_MIR, 
     * they use the same coefficients
     * unit = ^oC 
     * Ref. T20-6B    SBRS Doc. 151840 REV C
     */             

     /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    { /* AQUA (FM1) C_MIR_A[6] */
      8.6330e+01, -6.9256e-02, 4.0275e-05, -1.6550e-08, 
                  3.5331e-12, -3.1087e-16 },   
    { /* AQUA (FM1) C_MIR_B[6] */
      8.7354E+01,	-7.0526E-02,	4.1313E-05,	-1.6883E-08,	
                  3.5664E-12,	-3.0902E-16 },  
    /*
     *   T_cav --> TP_MF_CALBKHD_SR, TP_SA_A_MTR
     *   unit = ^oC
     *   Ref. T20-6
     */
 
     /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    { /* AQUA (FM1) C_CAV_A[6] */
      8.6334e+01, -5.5420e-01, 2.5791e-03, -8.4794e-06, 
                  1.4482e-08, -1.0193e-11 },    
    { /* AQUA (FM1) C_CAV_B[6] */
      8.7352e+01, -5.6419e-01, 2.6441e-03, -8.6454e-06, 
                  1.4614e-08, -1.0132e-11 },  

    /*
     *  T_ADC --> TA_PVLW_PWB4_10
     *              TA_PVSM_PWB6_12
     *   unit = ^oC
     *   Ref. T20-5B    SBRS Doc. 151840 REV C
     */

    { /* AQUA (FM1) C_ADC[6] */
      -1.8888e+01, 4.9724e-01, -2.3834e-03, 8.1515e-06, 
                   -1.3848e-08, 9.2793e-12},   
    /*
     *   T_INS --> TP_AO_SMIR_OBJ, TP_AO_SMIR_LENS, 
     *             TP_AO_LWIR_OBJ, TP_AO_LWIR_LENS
     *   unit = ^oC
     *   Ref. T20-6B    SBRS Doc. 151840 REV C
     */ 

    /* Common A/B side for Terra (PFM), Separate A/B sides for Aqua (FM1). */
    { /* AQUA (FM1) C_INS_A[6] */
      8.6330e+01, -6.9256e-02, 4.0275e-05, -1.6550e-08, 
                 3.5331e-12, -3.1087e-16},  
    { /* AQUA (FM1) C_INS_B[6] */
      8.7354e+01, -7.0526e-02, 4.1313e-05, -1.6883e-08, 
                 3.5664e-12, -3.0902e-16},  
                 
    /* RC temperature conversion
     * num bits = 12 (max DN value = 4095)
     * units = Kelvin
     *    VR_RC_LW_FPA_HTR
     *    VR_RC_SM_FPA_HTR
     * ref. T20-5B    SBRS Doc. 151840 REV C
     */
          
    { /* AQUA (FM1) C_RC[NUM_T_RC_VALUES][6] */
     {-2.6607e+00, 2.9239e-02, -1.2789e-06, 1.4882e-10, 0., 0.}, /* TA_RC_CS */
     {-2.0315e+02, 1.3710e-01, -2.1847e-05, 5.4693e-09, 0., 0.}, /* TA_RC_CS_OG */
     { 3.1238e+01, 2.6812e-02,  1.8002e-06, 2.8576e-11, 0., 0.}, /* TA_RC_IS */
     {-2.0199e+02, 1.3610e-01, -2.1620e-05, 5.4342e-09, 0., 0.}, /* TA_RC_IS_OG */
     {-2.0916e+02, 1.4514e-01, -2.5230e-05, 5.8560e-09, 0., 0.} /* TA_RC_OS_OG */
    }, 

    /* RC FPA voltage heater conversion (LWIR or SMIR)
     * num bits = 10 (max DN value = 1023)
     * units = volts DC
     * ref T20-5B    SBRS Doc. 151840 REV C
     */
    { /* AQUA (FM1) C_RC_VR_FPA[6] */ 
      -1.4970e+01, 2.9238e-02, 0., 0., 0., 0. },  
                   
    { /* AQUA (FM1) BB_a[NUM_BB_THERMISTORS][4] */
     {1.056236E-03, 2.599020E-04, -5.425454E-08, 1.415182E-07},
     {1.044243E-03, 2.644832E-04, -5.928912E-07, 1.628909E-07},
     {1.047991E-03, 2.626818E-04, -3.507229E-07, 1.533612E-07},
     {1.023606E-03, 2.710050E-04, -1.348365E-06, 1.935789E-07},
     {9.557013E-04, 2.970908E-04, -4.471715E-06, 3.073990E-07},
     {9.921567E-04, 2.816634E-04, -2.570464E-06, 2.376100E-07},
     {9.778825E-04, 2.874839E-04, -3.200168E-06, 2.605131E-07},
     {9.724649E-04, 2.883692E-04, -3.300496E-06, 2.659593E-07},
     {1.061567E-03, 2.579380E-04,  1.435982E-07, 1.358619E-07},
     {1.048034E-03, 2.628100E-04, -3.993544E-07, 1.539620E-07},
     {1.103644E-03, 2.444421E-04,  1.686106E-06, 7.536176E-08},
     {1.083254E-03, 2.507835E-04,  9.546836E-07, 1.038510E-07}
    }     
   };

  typedef struct {
    float32 moon_vector[MAX_NUM_SCANS][VEC_DIM];
    int8    moon_in_SV_KOB[MAX_NUM_SCANS][NUM_BANDS];
  } Moon_arrays_t;

  typedef struct 
  {
    float32 fp      [NUM_FOCAL_PLANES][MAX_NUM_SCANS];
    float32 mir1                      [MAX_NUM_SCANS];
    float32 mir2                      [MAX_NUM_SCANS];
    float32 mir_avg                   [MAX_NUM_SCANS];
    float32 scn_mtr                   [MAX_NUM_SCANS];
    float32 cav                       [MAX_NUM_SCANS];
    float32 ins                       [MAX_NUM_SCANS];
    float32 bb    [NUM_BB_THERMISTORS][MAX_NUM_SCANS];
    float32 bb_avg                    [MAX_NUM_SCANS];
    int8    fp_set_point_state        [MAX_NUM_SCANS];
    uint16  lwir_adc_index            [MAX_NUM_SCANS];
    uint16  mwir_adc_index            [MAX_NUM_SCANS];
    uint16  nir_adc_index             [MAX_NUM_SCANS];    
    uint16  vis_adc_index             [MAX_NUM_SCANS];
    uint16  pclw_adc_index            [MAX_NUM_SCANS];   
    uint8   num_thermistor_outliers   [MAX_NUM_SCANS];     /* range [0-12] */
    float32 ins_temp[NUM_T_INS_THERMISTORS][MAX_NUM_SCANS];
    float32 cav_temp[NUM_T_CAV_THERMISTORS][MAX_NUM_SCANS];
    float32 rc_temp       [NUM_T_RC_VALUES][MAX_NUM_SCANS];
    float32 vr_lwir                        [MAX_NUM_SCANS];
  } Temperatures_t; 

  typedef struct 
  {
    int32 num_scans;
    int16 BB_1km_night[MAX_1KM_TRACK_DIM]
                      [NUM_1000M_NIGHT_BANDS]
                      [BB_1km_FRAMES];
    int16 SV_1km_night[MAX_1KM_TRACK_DIM]
                      [NUM_1000M_NIGHT_BANDS]
                      [SV_1km_FRAMES];
    int16 MirrorSide[MAX_NUM_SCANS];
    boolean  Ecal_On[MAX_NUM_SCANS][NUM_BANDS];
    boolean  Sector_Rotated[MAX_NUM_SCANS]; 
    Temperatures_t temps;
  } L1A_granule_OBCEng_t;

  typedef struct 
  {
    int32 num_scans;
    int16 BB_1km_night[MAX_OVERLAP_TRACK_DIM]
                      [NUM_1000M_NIGHT_BANDS]
                      [BB_1km_FRAMES];
    int16 SV_1km_night[MAX_OVERLAP_TRACK_DIM]
                      [NUM_1000M_NIGHT_BANDS]
                      [SV_1km_FRAMES];
    int16  MirrorSide[MAX_NUM_SCANS];
    boolean Ecal_On[MAX_NUM_SCANS][NUM_BANDS];
    boolean Sector_Rotated[MAX_NUM_SCANS];
    Temperatures_t temps;
  } Overlap_OBCEng_t;

  void init_L1A_granule_OBCEng 
                        (L1A_granule_OBCEng_t *);

  PGSt_SMF_status Calculate_Temp_QA 
                        (int32                num_scans,
                         Temperatures_t       *temps,
			       	           QA_tables_t          *tables,
                         QA_Data_t            *QA);

  PGSt_SMF_status Process_OBCEng_Refl
                        (int32                sd_id, 
                         int32                num_scans,
                         float32              *T_ins, 
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,  
                         Moon_arrays_t        *moon_arrays,
                         L1A_granule_OBCEng_t *L1A_OBCEng, 
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl,
                         Preprocess_Refl_t    *PP_Refl);

  PGSt_SMF_status Process_OBCEng_Emiss 
                        (L1A_granule_t        *L1A_Gran,
                         L1A_granule_OBCEng_t *L1A_OBCEng,
                         Emiss_Cal_Coeff_t    *RVS_Coeff,
                         emiss_tables_t       *tables,
                         QA_tables_t          *QA_tables,
                         Moon_arrays_t        *moon_arrays,
                         QA_Data_t            *QA,
                         Preprocess_Emiss_t   *PP,
                         DN_OBC_Avg_t         *DN_OBC_Avg);

  PGSt_SMF_status Get_Emiss_Coeff_Per_Scan 
                        (int8                 moon_in_SV_KOB,
                         int16                *BB,
                         int16                *SV,
                         int16                B,
                         int16                D,
                         int16                D_emiss,
                         int16                MS,
                         Emiss_Cal_Coeff_t    *RVS_Coeff,
                         emiss_tables_t       *tables,
                         float32              T_bb,
                         float32              T_mir,
                         float32              T_cav,
                         float32              T_ins,
                         float32              T_fp_lwir,
                         float32              *Xdn_bb_31,
                         float32              *Xb1,
                         float32              *xdLbb,
                         float32              *xdnbb,
                         float32              *xLbb,
                         float32              *xLcav,
                         float32              *xdnsv,
                         float32              *xdnsv_var,
                         uint32               *sv_omask,
                         float32              *SNR,
                         int32                satellite_ID);

  PGSt_SMF_status Calculate_Planck 
                        (float32              *RSR,
                         float32              *wl,
                         int16                size,
                         float32              T,
                         float32              *planck);

  PGSt_SMF_status Read_L1A_OBCEng 
                        (emiss_tables_t       *emiss_tables,
                         L1A_granule_t        *L1A_Gran,
                         L1A_granule_OBCEng_t *L1A_OBCEng);

  PGSt_SMF_status Read_Overlap_OBCEng 
                        (L1A_granule_t        *L1A_Gran,
                         emiss_tables_t       *emiss_tables,
                         PGSt_PC_Logical      lun,
                         Overlap_OBCEng_t     *Overlap_OBCEng,
                         boolean              *scan_gap);


  PGSt_SMF_status Write_L1B_OBCEng 
                        (L1A_granule_t       *L1A_Gran,
                         lookup_tables_t     *tables,
                         Moon_arrays_t       *moon_arrays,
                         DN_OBC_Avg_t        *DN_OBC_Avg);

  PGSt_SMF_status Copy_EngMemData 
                        (int32                in_sd_id,
                         int32                out_sd_id,
                         int32                num_scans);

  PGSt_SMF_status Copy_ScanMetadata  
                        (int32                in_sd_id,
                         int32                out_sd_id,
                         int32                num_scans);

  PGSt_SMF_status Copy_PixelQualityData 
                        (int32                in_sd_id,
                         int32                out_sd_id,
                         int32                num_scans);

  PGSt_SMF_status Copy_EngVdata   
                        (int32                in_v_id,
                         int32                out_v_id);

  PGSt_SMF_status Read_Convert_Temperatures 
                        (emiss_tables_t       *tables,
                         int32                v_id,
                         int32                start_scan,
                         int32                num_scans,
                         int32                satellite_ID, 
                         Temperatures_t       *temps);

  PGSt_SMF_status Calculate_PP_Planck_Mir
                        (int32                num_scans,
                         emiss_tables_t       *tables,
                         Preprocess_Emiss_t   *PP);

  PGSt_SMF_status Get_Leading_Gran_Emiss_Coeff
                        (Overlap_OBCEng_t     *Leading_OBCEng,
                         int16                B,
                         int16                D,
                         int16                D_emiss,
                         Emiss_Cal_Coeff_t    *RVS_Coeff,
                         emiss_tables_t       *tables,
                         Moon_arrays_t        *moon_arrays,
                         float32              Xdn_bb_31[]
                                                [DETECTORS_PER_1KM_BAND],
                         float32              *Xb1,
                         int16                *XMS,
                         int32                satellite_ID);

  PGSt_SMF_status Get_Middle_Gran_Emiss_Coeff
                        (L1A_granule_OBCEng_t *L1A_OBCEng,
                         int16                B,
                         int16                D,
                         int16                D_emiss,
                         Emiss_Cal_Coeff_t    *RVS_Coeff,
                         emiss_tables_t       *tables,
                         emiss_QA_tables_t    *QA_tables,
                         Moon_arrays_t        *moon_arrays,
                         float32              Xdn_bb_31[]
                                                [DETECTORS_PER_1KM_BAND],
                         float32              *Xb1,
                         int16                *XMS,
                         QA_Emiss_t           *QA,
                         Preprocess_Emiss_t   *PP,
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         int32                satellite_ID);

  PGSt_SMF_status Get_Trailing_Gran_Emiss_Coeff
                        (int32                num_scans_middle,
                         Overlap_OBCEng_t     *Trailing_OBCEng,
                         int16                B,
                         int16                D,
                         int16                D_emiss,
                         Emiss_Cal_Coeff_t    *RVS_Coeff,
                         emiss_tables_t       *tables,
                         Moon_arrays_t        *moon_arrays,
                         float32              Xdn_bb_31[]
                                                [DETECTORS_PER_1KM_BAND],
                         float32              *Xb1,
                         int16                *XMS,
                         int32                satellite_ID);

  PGSt_SMF_status Get_All_Emiss_Coeff
                        (L1A_granule_OBCEng_t *L1A_OBCEng,
                         int16                B,
                         int16                D_emiss,
                         float32              *Xb1,
                         int16                *Xmir,
                         emiss_tables_t       *tables,
                         emiss_QA_tables_t    *QA_tables,
                         QA_Emiss_t           *QA,
                         Preprocess_Emiss_t   *PP);

  PGSt_SMF_status Compute_BB_Temperature
                        (emiss_tables_t       *tables,
                         int32                start_scan,
                         int32                num_scans,
                         uint16               *CP_index,
                         uint16               **raw_BB_temp,
                         int32                satellite_ID, 
                         Temperatures_t       *temps); 

  PGSt_SMF_status Get_Electronics_index
                        (int32                num_scans,
                         uint16               *Electronics_pri,
                         uint16               *Electronics_red,
                         uint16               *Electronics_index);

  PGSt_SMF_status Write_Geo_OBC_SDS
                        (int32                obc_sd_id,
                         int32                num_scans);

  PGSt_SMF_status Fill_250m_DN_OBC_Avg
                        (int32                sd_id,
                         int32                num_scans,
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,
                         Moon_arrays_t        *moon_arrays,
                         boolean              Ecal_On[][NUM_BANDS],
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl);

  PGSt_SMF_status Fill_500m_DN_OBC_Avg
                        (int32                sd_id,
                         int32                num_scans,
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,
                         Moon_arrays_t        *moon_arrays,
                         boolean              Ecal_On[][NUM_BANDS],
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl);

  PGSt_SMF_status Fill_1km_day_DN_OBC_Avg
                        (int32                sd_id, 
                         int32                num_scans,
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,
                         Moon_arrays_t        *moon_arrays,
                         boolean              Ecal_On[][NUM_BANDS],
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl);

  PGSt_SMF_status Fill_Band_26_DN_OBC_Avg
                        (int32                sd_id,
                         int32                num_scans,
                         int16                *MirrorSide,
                         refl_tables_t        *refl_tables,
                         Moon_arrays_t        *moon_arrays,
                         boolean              Ecal_On[][NUM_BANDS],
                         DN_OBC_Avg_t         *DN_OBC_Avg,
                         QA_Refl_t            *QA_refl);

  PGSt_SMF_status Check_For_Moon_in_SV_KOB
                        (int32                num_scans,
                         common_QA_tables_t   *common_QA_tables,
                         QA_Data_t            *QA,
                         Moon_arrays_t        *moon_arrays);

  PGSt_SMF_status Get_DN_Avg_SDev_Rejects
                        (int32                start_index,
                         int32                N,
                         int32                index_increment,
                         int16                *DN_array,
                         int16                DN_upper_valid_limit,
                         int16                DN_lower_valid_limit,
                         float32              *mean_DN,
                         float32              *sdev_DN,
                         int16                *rejects);

  PGSt_SMF_status Pack_Rejects_In_Outlier_Mask
                        (int32                N,
                         int16                *rejects,
                         uint32               *packed_rejects);

  PGSt_SMF_status Cross_Granule_Sliding_Average
                        (int32                num_overlapScans,
                         int32                num_scans_middle,
                         float32              *Xarray,
                         int16                *Xmir,
                         float32              *avg);

  PGSt_SMF_status Granule_Average_Temperature
                        (int32                num_values,
                         float32              *array,
                         float32              *avg);

  PGSt_SMF_status Get_Temp_Avg_And_Variance
                        (int32                N,
                         float32              *T,
                         float32              *avg,
                         float32              *var);

  PGSt_SMF_status Fill_Invalid_Temp_DNs
                        (int32                num_scans,
                         uint16               miss_DN,
                         uint16               sat_DN,
                         uint16               *DN);

  PGSt_SMF_status Adjust_dn_star_Min
                        (float32              *dn_star_Min,
                         int32                L1A_v_id,
                         int32                num_scans);

  PGSt_SMF_status Check_For_Ecal_On
                        (int32                lun,
                         int32                num_scans,
                         int32                v_id,
                         boolean              Ecal_On[][NUM_BANDS]);

  PGSt_SMF_status Check_If_Sector_Rotated
                        (int32                lun,
                         int32                num_scans,
                         int32                v_id,
                         boolean              Sector_Rotated[]);

  PGSt_SMF_status sort_int16_array 
                        (int32                n, 
                         int16                *a,
                         int32                *indx);

  PGSt_SMF_status Get_DN_Avg_SDev_Rejects_LowN
                        (int32                N_include,
                         int32                start_index,
                         int32                N,
                         int32                index_increment,
                         int16                *DN_array,
                         int16                DN_upper_valid_limit,
                         int16                DN_lower_valid_limit,
                         float32              *mean_DN,
                         float32              *sdev_DN,
                         int16                *rejects);

  PGSt_SMF_status Calculate_RVS_Correction
                        (lookup_tables_t      *tables,
                         L1B_granule_t        *L1B_Gran);


#endif



