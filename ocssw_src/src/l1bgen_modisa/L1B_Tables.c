#include    "L1B_Tables.h" 
#include    "HDF_Lib.h"
#include    "PGS_PC.h"
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    "HDF_Lib.h"         /* for error messages invalidinputfile
                                   and corruptinputfile */ 
#include    <math.h>
#include    <stdlib.h>

/***************************************************************************
Developer's note: (Jim Rogers)

In the functions "Read_..._Tables", there is a return statement after each
L1BErrorMsg function call.  These are actually not necessary if the last
argument to L1BErrorMsg is "True". These returns make it easier to test
error out conditions (with SMF_ERROR used in UNIT_TEST_MODE_ONLY -- it does
not actually exit).
****************************************************************************/

extern int16 RFLAG;

static char *invalidlutfile = 
               "This is most likely due to an invalid LUT file.";
static PGSt_SMF_status Read_LUT_Tables 
                           (L1A_granule_t      *L1A_Gran,
                            int32              lun,
                            LUT_Definition_t   *luts);

PGSt_SMF_status  Read_Lookup_Tables 
                           (L1A_granule_t         *L1A_Gran,
                            lookup_tables_t       *tables,
                            Run_Time_Parameters_t *runtime_params)
/*
!C***************************************************************
!Description:      Read all L1B lookup tables. 

!Input Parameters:
        L1A_granule_t         *L1A_Gran        contains satellite id and 
                                               EV start time 
        lookup_tables_t       *tables
        Run_Time_Parameters_t *runtime_params  contains LUT version for check    

!Output Parameters:
        lookup_tables_t     *tables

!Revision History:

 Revision 02.20, October 31, 2003 (Razor Issue #195)
 Delete LUT ProcessingCenter from the QA LUT array.
 Liqin Tan  SAIC/GSO (ltan@saicmodis.com)

 Revision 02.19 , October 24, 2003  Razor Issue #196 (formerly Issue #184)
 Added "ROLL_THRESHOLD_LUT_NAME", 'PITCH_THRESHOLD_LUT_NAME', and
 "YAW_THRESHOLD_LUT_NAME" in the LUT_Definition_t array "qa_luts[]".
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 02.18   April 16, 2003
 Changed lower bound on m1_table LUT values to a small positive epsilon.
 Since the m1 values are used as divisors to generate the R* LUT values
 (see revision 02.16), they should never be allowed to be 0.
 Alice Isaacman, SAIC GSO    (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.17   March 26, 2003   Razor Issue #191
 Added LUT for SWIR OOB correction sending band.
 Alice Isaacman, SAIC GSO    (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.16, Oct 3, 2002  Razor Issue #187
 Removed R_star from code and LUTs. Added zero divide error checks for m1 and
   E_sun_over_pi  
 Liqin Tan, SAIC GSO (ltan@saicmodis.com) 

 Revision 02.15, April 23, 2001  Razor Issue #167
 Changed references to "PDF" to "PCF"
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.14, March 25, 2001  Razor Issue #178
 Removed ADC Correction LUTs
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.13 November 6, 2001  Razor issue #167
 Added check on MCST Version read from runtime parameters
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

 Revision 02.12 March 5, 2001  Razor issue #156
 Added new function TDLUT_ReadPiecewiseLinearFunction 
    and appropriate call.
 Added check on reported vs. actual rank in TDLUT_ReadStepFunction
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 02.11 November 17, 1999
 Added checking that MCST versions are consistent.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Apr 1998
 Removed Malloc_Lookup_Tables().
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.00 Jan 1998
 Added the QA component of lookup tables.
 Zhidong Hao (hao@gscmail.gsfc.nasa.gov)

 Revision 01.00 Dec 1996
 Initial development.
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Read_Lookup_Tables";

  /*
   * Read tables
   */ 
  returnStatus = Read_Refl_Tables (L1A_Gran, &tables->refl);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, 
                "Read_Refl_Tables", 0, NULL, True);
 
  returnStatus = Read_Emiss_Tables (L1A_Gran, &tables->emiss);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, 
                "Read_Emiss_Tables", 0, NULL, True);

  returnStatus = Read_QA_Tables (L1A_Gran, &tables->QA);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, 
                "Read_QA_Tables", 0, NULL, True);

#ifndef NOCHECKLUT
  /*
   * The table value of PGE version has already been checked against the
   * code macro PGE02_VERSION.  Each LUT file also contains the MCST
   * version.  Check that these are all the same.
   */

  if (strcmp(tables->refl.MCST_Version, 
             tables->QA.common_QA_tables.MCST_Version)) 
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "MCST version in reflective LUT file does "
                "not match that in QA LUT file.",
                NULL, 0,
                "LUT files are invalid.  All must have the "
                "same MCST version.",
                True);
    return returnStatus;
  }
  if (strcmp(tables->emiss.MCST_Version, 
      tables->QA.common_QA_tables.MCST_Version)) 
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "MCST version in emissive LUT file does not "
                "match that in QA LUT file.",
                NULL, 0,
                "LUT files are invalid.  All must have the "
                "same MCST version.",
                True);
    return returnStatus;
  }

  /*
   * The MCST versions of the LUTs are identical; now check the MCST version 
   * against the MCST version given in the PCF file.
   */
        
  if (strcmp(runtime_params->MCST_LUT_Version, 
            tables->QA.common_QA_tables.MCST_Version)) 
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
      "MCST version in PCF file does not match that in LUT files.",
      NULL, 0,
      "LUT files must have the same MCST version as "
                "specified in the PCF file.",
      True);
    return returnStatus;
  }
#endif

  return(MODIS_S_OK);
}

/********************* MACRO ASSIGN_DATA_PTR ********************************
This macro assigns the dataptr to the luts[i].data member.
This macro is used in Read_Refl_Tables, Read_Emiss_Tables and Read_QA_Tables.
There should be the following local variables in those functions:
  i       (loop counter)
  lun     (assigned to be the logical unit number)
Macro variables:
  luts    (array of luts, either refl_luts, emiss_luts or qa_luts)
  lutname (name of particular LUT)
  dataptr (structure member to be assigned.)
****************************************************************************/

#define ASSIGN_DATA_PTR(luts,lutname,dataptr)                           \
  i = 0;                                                                \
  while (luts[i].name)                                                  \
  {                                                                     \
    if (!strcmp(lutname,luts[i].name))                                  \
      break;                                                            \
    i++;                                                                \
  }                                                                     \
  if (!luts[i].name)                                                    \
  {                                                                     \
    char errmsg[512];                                                   \
    sprintf(errmsg, "LUT name %s not found in LUTs array.", lutname);   \
    returnStatus = MODIS_F_NOK;                                         \
    L1BErrorMsg(location, returnStatus, errmsg, NULL,                   \
                lun, "*** CODE BUG ***", True);                         \
    return returnStatus;                                                \
  }                                                                     \
  luts[i].data = (VOIDP) dataptr;


PGSt_SMF_status Read_Refl_Tables (L1A_granule_t *L1A_Gran, 
                                  refl_tables_t *tables)
/*
!C***************************************************************
!Description:      Read all reflective lookup tables. 

!Input Parameters:
        L1A_granule_t     *L1A_Gran      contains satellite id and 
                                         EV start time
	refl_tables_t     *tables        address of empty tables       

!Output Parameters:
        refl_tables_t     *tables        address of filled tables 

!Revision History:
 $Log: L1B_Tables.c,v $
 Revision 1.35  2017-04-20 09:41:55-04  xgeng
 PGE02 update for MODIS processing on CentOS 7 system

 Revision 1.34  2011-09-07 09:47:18-04  xgeng
 Added 3 new luts for default b1 algorithm change

 Revision 1.31  2011-04-07 14:38:42-04  xgeng
 1. RSB &TEB uncertainty algorithm update; 2. The quadratic RSB RVS changed to 4th order.

 Revision 1.29  2008-12-16 16:37:26-05  xgeng
 Extended rvs ascii lut validation limitation

 Revision 1.27  2008/11/18 16:28:20  xgeng
 merge the branch for V6.0.1

 Revision 1.26.2.3  2008/06/02 15:39:54  xgeng
 1.Added a new table piece "DET_QUAL_FLAG2_VALS_LUT_NAME" into the LUT_Definition_t array "qa_luts[]". 2.Added the calculation of Detector_Quality_Flag2. 3.Generated noisy and dead subframe list.

 Revision 1.26  2005/01/18 21:57:37  ltan
 Added new file attributes prefixed with HDFEOS_FractionalOffset

 
 Revision 02.42  October 15, 2004  Razor Issue #199
 Added "swir_oob_sending_detector_table" to the Reflective LUTs.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 02.41  March 26, 2003   Razor Issue #191
 Check range of sending band to use for SWIR OOB correction
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.41  March 26, 2003   Razor Issue #190
 Added LUTs "B26_B5_Corr", "B26_B5_Corr_Switch",
 "B26_B5_Frame_Offset" (previously added to Terra code; 
 Razor Issue #182.)
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.39 June 28, 2002    Razor Issue #161
 Handle BDSM LUTs without fill values.
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 02.38  June 5, 2002  Razor Issue #183
 Change type of dn_sat_ev to float64.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.37, March 25, 2001  Razor Issue #178
 Removed ADC Correction LUTs
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 02.35, March 8, 2002           Razor Issue #174
 Removed LUTs "RVS_250m", "RVS_500m", "RVS_1km_RefSB"
 Added LUT "RVS_RefSB"
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 02.33 January 27, 2002 
 Razor Issue #175:  Added new LUTs BB_T_SAT_SWITCH_AQUA, BB_T_SAT_AQUA,
   B1_DEFAULT_B33_35_36 to emiss_tables_t for use with Aqua

 Revision 02.34, March 26, 2001
 Razor Issue #179: Changed maximum number of allowable wavelengths for 
   RSR tables to 66 to accomodate Aqua; changed fill value to 0.0.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.32, March 26, 2001
 Added new LUT dn_sat_ev and removed LUT DN_sat
 Alice Isaacman, SAIC GSC (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 02.31, December 22, 1999
 Redid the tables array to make it work for LUT generation also.
 Added Read_LUT_Tables to avoid a lot of code duplication.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.30, November 17, 1999
 Completely revised methodology of reading in tables.  Using a loop to
 allow for better error messages and for valid range checking (if enabled).
 Checks for PGE version consistency also added.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.12 Sept. 03, 1999
 Replaced all old SWIR correction tables with new ones due to the new 
 SWIR correction algorithm.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) 

 Revision 02.11 August 1999
 Added delta_DN_RSB
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Apr. 9, 1998
 Changed Time_gain_factor to Time_gain_factor_rad
 and added Time_gain_factor_refl. 
 Added L_SD table read.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Apr. 8, 1998
 Changed the read table calls to the new generalized interfaces.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Apr. 2, 1998
 V2.1 reflective tables added.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Nov. 22, 1996 
 Initial development.
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 band, det, sample, mirr_side;
  int32 i, lun = REFLECTIVE_TABLES_FILE;
  char *location = "Read_Refl_Tables";
  SWIR_correction_tables_t *swir_tables = &tables->SWIR_correction_tables;

  /*
   * Assign the array pointers to the data member in refl_luts.
   */ 

  ASSIGN_DATA_PTR(refl_luts, 
      PGE_VERSION_LUT_NAME,        
          tables->PGE_Version);

  ASSIGN_DATA_PTR(refl_luts, 
      MCST_VERSION_LUT_NAME,      
          tables->MCST_Version);

  ASSIGN_DATA_PTR(refl_luts, 
      REFL_SERIAL_NUMBER_LUT_NAME,
          tables->Serial_Number);

  ASSIGN_DATA_PTR(refl_luts, 
      M0_LUT_NAME,                
          tables->m0);

  ASSIGN_DATA_PTR(refl_luts, 
      M1_LUT_NAME,                
          tables->m1);

  ASSIGN_DATA_PTR(refl_luts, 
      DN_STAR_MAX_LUT_NAME,       
          tables->dn_star_Max);

  ASSIGN_DATA_PTR(refl_luts, 
      DN_STAR_MIN_LUT_NAME,       
          tables->dn_star_Min);

  ASSIGN_DATA_PTR(refl_luts, 
      K_INST_LUT_NAME,            
          tables->K_inst);

  ASSIGN_DATA_PTR(refl_luts, 
      K_FPA_LUT_NAME,             
          tables->K_FPA);

  ASSIGN_DATA_PTR(refl_luts, 
      T_INST_REF_LUT_NAME,         
         &tables->T_inst_ref);

  ASSIGN_DATA_PTR(refl_luts, 
      T_FPA_REF_LUT_NAME,         
          tables->T_FPA_ref);

  ASSIGN_DATA_PTR(refl_luts, 
      RVS_RSB_LUT_NAME,          
          tables->RVS_RefSB);

  /*
   * obsolete due to RSB UI algorithm update, 2/23/2011, Xu Geng
   */
  /*
  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_RVS_RSB_LUT_NAME,     
          tables->Sigma_RVS_RSB);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_M1_LUT_NAME,          
          tables->Sigma_m1);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_K_INST_LUT_NAME,      
          tables->Sigma_K_inst);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_T_INST_LUT_NAME,       
         &tables->Sigma_T_inst);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_PV_RESID_ELEC_LUT_NAME,
          tables->Sigma_PV_Resid_Elec);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_R_STAR_LIN_RESID_UCOEFF_LUT_NAME, 
          tables->Sigma_R_Star_Lin_Resid_Ucoeff);

  ASSIGN_DATA_PTR(refl_luts, 
      RSB_NEDL_LUT_NAME,          
          tables->RSB_NEdL);

  ASSIGN_DATA_PTR(refl_luts, 
      SIGMA_RSB_ADC_LUT_NAME,     
          tables->Sigma_RSB_ADC); 
  */

  ASSIGN_DATA_PTR(refl_luts, 
      DN_OBC_1ST_FRAME_LUT_NAME,   
         &tables->DN_obc_avg_first_frame_to_use);

  ASSIGN_DATA_PTR(refl_luts, 
      DN_OBC_NUM_FRAMES_LUT_NAME,  
         &tables->DN_obc_avg_number_of_frames_to_use);

  ASSIGN_DATA_PTR(refl_luts, 
      SWIR_CORRECTION_SWITCH_LUT_NAME, 
          &swir_tables->SWIR_correction_switch);

  ASSIGN_DATA_PTR(refl_luts, 
      SWIR_CORR_SENDING_BAND_LUT_NAME, 
          &swir_tables->SWIR_corr_sending_band);

  ASSIGN_DATA_PTR(refl_luts,
      SWIR_CORR_SENDING_DETECTOR_LUT_NAME,
          &swir_tables->SWIR_corr_sending_detector);

  ASSIGN_DATA_PTR(refl_luts, 
      X_OOB_0_LUT_NAME,            
          swir_tables->X_OOB_0);

  ASSIGN_DATA_PTR(refl_luts, 
      X_OOB_1_LUT_NAME,            
          swir_tables->X_OOB_1);

  ASSIGN_DATA_PTR(refl_luts, 
      X_OOB_2_LUT_NAME,            
          swir_tables->X_OOB_2);

  ASSIGN_DATA_PTR(refl_luts, 
      RSB_SPECIFIED_UNCERTAINTY_LUT_NAME, 
          tables->RSB_specified_uncertainty);

  ASSIGN_DATA_PTR(refl_luts, 
      RSB_UI_SCALING_FACTOR_LUT_NAME, 
          tables->RSB_UI_scaling_factor);

  ASSIGN_DATA_PTR(refl_luts, 
      E_SUN_OVER_PI_LUT_NAME,     
          tables->E_sun_over_pi);

  ASSIGN_DATA_PTR(refl_luts, 
      RSB_SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME, 
          &tables->RSB_SV_DN_moon_include_frames);

  ASSIGN_DATA_PTR(refl_luts, 
      DN_SAT_EV_LUT_NAME,         
          tables->dn_sat_ev);

  ASSIGN_DATA_PTR(refl_luts, 
      B26_B5_CORR_LUT_NAME,         
          tables->B26_B5_Corr);
  
  ASSIGN_DATA_PTR(refl_luts, 
      B26_B5_CORR_SWITCH_LUT_NAME,         
          &tables->B26_B5_Corr_Switch);
  
  ASSIGN_DATA_PTR(refl_luts, 
      B26_B5_FRAME_OFFSET_LUT_NAME,         
          tables->B26_B5_Frame_Offset);

  ASSIGN_DATA_PTR(refl_luts,
      U1_LUT_NAME,
          tables->u1);

  ASSIGN_DATA_PTR(refl_luts,
      U2_LUT_NAME,
          tables->u2_samples);

  ASSIGN_DATA_PTR(refl_luts,
      U3_LUT_NAME,
          tables->u3);

  ASSIGN_DATA_PTR(refl_luts,
      U4_LUT_NAME,
          tables->u4_coeffs);

  ASSIGN_DATA_PTR(refl_luts,
      U2_FRAMES_LUT_NAME,
          tables->u2_frames);

  ASSIGN_DATA_PTR(refl_luts,
      SWIR_UI_FACTOR_LUT_NAME,
          tables->swir_ui_factor);
  
  /*
   * Read all tables.
   */

  returnStatus = Read_LUT_Tables(L1A_Gran, 
                                 REFLECTIVE_TABLES_FILE, 
                                 refl_luts);
  if (returnStatus != MODIS_S_OK) 
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Read_LUT_Tables", REFLECTIVE_TABLES_FILE, 
                NULL, True);
    return returnStatus;
  }


  /*
   * Check some quantities that could cause code failures and normalize the
   * Sigma RVS quantities by the RVS.  If full checking has been enabled,
   * some of these checks may be redundant.  However, in operations, we do not
   * expect to turn on the full checking (that is for LUT generation only).
   */      

    /* MODIS_BAND20_INDEX = 21    NUM_BANDS = 38    MODIS_BAND26_INDEX = 27  */

  if (swir_tables->SWIR_corr_sending_band < MODIS_BAND20_INDEX - 1 ||
      swir_tables->SWIR_corr_sending_band > NUM_BANDS - 1 || 
      swir_tables->SWIR_corr_sending_band == MODIS_BAND26_INDEX - 1)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "Band to use for SWIR OOB Correction is out of range.",
                NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }
  
  if (tables->DN_obc_avg_first_frame_to_use < 0 ||
      tables->DN_obc_avg_first_frame_to_use >= MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "DN_OBC first frame to use is out of range.",
                NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if (tables->DN_obc_avg_number_of_frames_to_use <= 0 || 
      (tables->DN_obc_avg_first_frame_to_use + 
       tables->DN_obc_avg_number_of_frames_to_use) > MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "DN_OBC number of frames to use is out of range "
                "or is inconsistent\n"
                "with the value of DN_OBC first frame to use.",
                NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  /* Ensure that m1 isn't zero (divide by m1 occurs in code) */
  for ( band = 0; band < NUM_REFLECTIVE_BANDS; band++ )
    for ( det = 0; det < MAX_DETECTORS_PER_BAND; det++ )
      for ( sample = 0; sample < MAX_SAMPLES_PER_BAND; sample++ )
        for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ )
          if ( fabs((double) tables->m1[band][det][sample][mirr_side]) 
              < TOLERANCE )
          {
            returnStatus = MODIS_F_OUT_OF_RANGE;
            L1BErrorMsg(location, returnStatus,
                  "Bad m1 LUT (Zero values detected).",
                  NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
            return returnStatus;
          }

  /* Ensure that E_sun_over_pi isn't zero (divide by E_sun_over_pi occurs in code) */
  for ( det = 0; det < NUM_REFLECTIVE_DETECTORS; det++ )
    if ( fabs((double) tables->E_sun_over_pi[det])
              < TOLERANCE )
          {
            returnStatus = MODIS_F_OUT_OF_RANGE;
            L1BErrorMsg(location, returnStatus,
                  "Bad E_sun_over_pi LUT (Zero values detected).",
                  NULL, REFLECTIVE_TABLES_FILE, invalidlutfile, True);
            return returnStatus;
          }

  /* 
   * Check if the LUT RSB_specified_uncertainty has zero or negative values. 
   * The values should be all positive.
   */

  for (band = 0; band < NUM_REFLECTIVE_BANDS; band++)
    if (tables->RSB_specified_uncertainty[band] <= TOLERANCE)
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus,   
                  "Bad RSB_specified_uncertainty LUT "
                  "(Zero or negative values detected).",
                  NULL, REFLECTIVE_TABLES_FILE, NULL, True);
      return returnStatus;
    }

  /* If no 250m bands then use rescaling */
  if ((RFLAG & 1) == 1) {
    tables->dn_star_Max[0] = 32767;
    tables->dn_star_Max[1] = 32767;
  }

  /* If no 500m bands then use rescaling */
  if ((RFLAG & 2) == 2) {
    tables->dn_star_Max[2] = 32767;
    tables->dn_star_Max[3] = 32767;
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Read_Emiss_Tables (L1A_granule_t  *L1A_Gran, 
                                   emiss_tables_t *tables)
/*
!C***************************************************************
!Description:      Read in all Emissive lookup tables.

!Input Parameters:
         L1A_granule_t      *L1A_Gran      contains satellite id and 
                                           EV start time
         emiss_tables_t     *tables        address of set of empty tables

!Output Parameters:
         emiss_tables_t     *tables        address of set of filled tables

!Revision History:
 Revision 02.36, October 15, 2004  Razor Issue #201
 Added the new dimension of Mirror Side to the LUT definition of band_21_b1
 for Emissive LUTs.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 02.35, March 8, 2002
 Removed LUTs "RVS_1km_Emiss_SV", "RVS_1km_Emiss_BB", "RVS_1km_Emiss_EV"
 Added LUTs "RVS_TEB", "RVS_BB_SV_Frame_No"
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.31, December 22, 1999
 Redid the tables array to make it work for LUT generation also.
 Added Read_LUT_Tables to avoid a lot of code duplication.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.30, November 17, 1999
 Completely revised methodology of reading in tables.  Using a loop to
 allow for better error messages and for valid range checking (if enabled).
 Checks for PGE version consistency also added.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.18 August 1999
 Removed LUT INT_correction_switch.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.17 August 1999
 Added delta_DN_TEB
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.16 August 12, 1999
 Added L_Max, L_Min Luts
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.15 April 12, 1999
 Removed cubic term a3
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Replaced the V2.0 tables by V2.1 tables input
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 01.00 Nov. 1996
 Initial development.
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 i, j, band, lun = EMISSIVE_TABLES_FILE;
  char *location = "Read_Emiss_Tables";

  /*
   * Assign the array pointers to the data member in emiss_luts.
   */ 


  ASSIGN_DATA_PTR(emiss_luts,
                     PGE_VERSION_LUT_NAME,
                         tables->PGE_Version);

  ASSIGN_DATA_PTR(emiss_luts,
                     MCST_VERSION_LUT_NAME,
                         tables->MCST_Version);

  ASSIGN_DATA_PTR(emiss_luts,
                     EMISS_SERIAL_NUMBER_LUT_NAME,
                         tables->Serial_Number);

  ASSIGN_DATA_PTR(emiss_luts,
                     EPSILON_BB_LUT_NAME,
                         tables->epsilon_bb);

  ASSIGN_DATA_PTR(emiss_luts,
                     EPSILON_CAV_LUT_NAME,
                         tables->epsilon_cav);

  ASSIGN_DATA_PTR(emiss_luts,
                     DELTA_T_BB_BETA_LUT_NAME,
                         tables->delta_T_bb_beta);

  ASSIGN_DATA_PTR(emiss_luts,
                     DELTA_T_BB_DELTA_LUT_NAME,
                         tables->delta_T_bb_delta);

  ASSIGN_DATA_PTR(emiss_luts,
                     PCX_TALK_LUT_NAME,
                         tables->PC_XT);

  ASSIGN_DATA_PTR(emiss_luts,
                     RVS_TEB_LUT_NAME,
                         tables->RVS_TEB);

  ASSIGN_DATA_PTR(emiss_luts,
                     RVS_BB_SV_FRAME_NO_LUT_NAME,
                         tables->RVS_BB_SV_Frame_No);

  ASSIGN_DATA_PTR(emiss_luts,
                     RSR_LUT_NAME,
                         tables->RSR);

  ASSIGN_DATA_PTR(emiss_luts,
                     WAVELENGTH_LUT_NAME,
                         tables->wavelength);

  ASSIGN_DATA_PTR(emiss_luts,
                     NUM_WL_INCREMENT_LUT_NAME,
                         tables->NUM_RSR_vs_Lambda);

  ASSIGN_DATA_PTR(emiss_luts,
                     CALIB_A0_LUT_NAME,
                         tables->A0);

  ASSIGN_DATA_PTR(emiss_luts,
                     CALIB_A2_LUT_NAME,
                         tables->A2);
  /*
   * obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng
   */
  /*

  ASSIGN_DATA_PTR(emiss_luts,
                     UI_UCOEFF_LUT_NAME,
                         tables->Ucoeff);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_TEB_PV_RESID_ELEC_LUT_NAME,
                         tables->Sigma_TEB_PV_resid_elec);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_TEB_ADC_LUT_NAME,
                         tables->Sigma_TEB_ADC);

  ASSIGN_DATA_PTR(emiss_luts,
                     UCOEFF_CALIBR_RESID_LUT_NAME,
                         tables->Ucoeff_Calibr_resid);

  ASSIGN_DATA_PTR(emiss_luts,
                     BAND_21_UNCERT_LSAT_LUT_NAME,
                             &tables->Band_21_Uncert_Lsat);

  */

  ASSIGN_DATA_PTR(emiss_luts,
                     T_INS_FUNCTION_FLAG_LUT_NAME,
                         tables->T_ins_function_flag);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_INS_DEFAULT_LUT_NAME,
                         &tables->T_ins_default);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_INS_OFFSET_LUT_NAME,
                         tables->T_ins_offset);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_CAV_FUNCTION_FLAG_LUT_NAME,
                         tables->T_cav_function_flag);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_CAV_DEFAULT_LUT_NAME,
                         &tables->T_cav_default);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_MIR_FUNCTION_FLAG_LUT_NAME,
                         tables->T_mir_function_flag);

  ASSIGN_DATA_PTR(emiss_luts,
                     T_MIR_DEFAULT_LUT_NAME,
                         &tables->T_mir_default);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_WEIGHT_LUT_NAME,
                         tables->BB_Weight);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_DN_1ST_FRAME_LUT_NAME,
                         &tables->BB_DN_first_frame_to_use);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_DN_NUM_FRAMES_LUT_NAME,
                         &tables->BB_DN_number_of_frames_to_use);

  ASSIGN_DATA_PTR(emiss_luts,
                     SV_DN_1ST_FRAME_LUT_NAME,
                         &tables->SV_DN_first_frame_to_use);

  ASSIGN_DATA_PTR(emiss_luts,
                     SV_DN_NUM_FRAMES_LUT_NAME,
                         &tables->SV_DN_number_of_frames_to_use);

  ASSIGN_DATA_PTR(emiss_luts,
                     SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME,
                         &tables->SV_DN_moon_include_frames);

  ASSIGN_DATA_PTR(emiss_luts,
                     PCX_CORRECTION_SWITCH_LUT_NAME,
                         &tables->PCX_correction_switch);

  ASSIGN_DATA_PTR(emiss_luts,
                     OVERLAP_SCANS_B1_LUT_NAME,
                         &tables->num_overlap_scans_b1);

  ASSIGN_DATA_PTR(emiss_luts,
                     BAND_21_B1_LUT_NAME,
                         tables->Band_21_b1);

  ASSIGN_DATA_PTR(emiss_luts,
                     L_MIN_LUT_NAME,
                         tables->L_Min);

  ASSIGN_DATA_PTR(emiss_luts,
                     L_MAX_LUT_NAME,
                         tables->L_Max);

  ASSIGN_DATA_PTR(emiss_luts,
                     TEB_SPECIFIED_UNCERTAINTY_LUT_NAME,
                         tables->TEB_specified_uncertainty);

  ASSIGN_DATA_PTR(emiss_luts,
                     TEB_UI_SCALING_FACTOR_LUT_NAME,
                         tables->TEB_UI_scaling_factor);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_T_SAT_SWITCH_AQUA_LUT_NAME,
                         &tables->BB_T_sat_switch_aqua);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_T_SAT_AQUA_LUT_NAME,
                         tables->BB_T_sat_aqua);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_T_SAT_DEFAULT_B1_BASELINE_AQUA_LUT_NAME,
                         tables->BB_T_sat_default_b1_baseline_aqua);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_T_SAT_DEFAULT_B1_C1_AQUA_LUT_NAME,
                         tables->BB_T_sat_default_b1_c1_aqua);

  ASSIGN_DATA_PTR(emiss_luts,
                     BB_T_SAT_DEFAULT_B1_TLWIR_BASELINE_AQUA_LUT_NAME,
                         &tables->BB_T_sat_default_b1_Tlwir_baseline_aqua);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_A0_LUT_NAME,
                         tables->sigma_a0);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_A2_LUT_NAME,
                         tables->sigma_a2);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_RVS_EV_LUT_NAME,
                         tables->sigma_RVS_EV);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_EPSILON_BB_LUT_NAME,
                         tables->sigma_epsilon_BB);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_EPSILON_CAV_LUT_NAME,
                         tables->sigma_epsilon_CAV);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_L_LAMBDA_LUT_NAME,
                         tables->sigma_L_lambda);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_L_TBB_LUT_NAME,
                         tables->sigma_L_Tbb);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_L_TSM_LUT_NAME,
                         tables->sigma_L_Tsm);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_L_TCAV_LUT_NAME,
                         tables->sigma_L_Tcav);

  ASSIGN_DATA_PTR(emiss_luts,
                     SIGMA_B1_BAND21_LUT_NAME,
                         tables->sigma_b1_B21);

  ASSIGN_DATA_PTR(emiss_luts,
                     PCX_UI_FACTOR_LUT_NAME,
                         tables->pcx_ui_factor);

  /*
   * Read all tables.
   */

  returnStatus = Read_LUT_Tables(L1A_Gran, 
                                 EMISSIVE_TABLES_FILE, 
                                 emiss_luts);
  if (returnStatus != MODIS_S_OK) 
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Read_LUT_Tables", EMISSIVE_TABLES_FILE, NULL, True);
    return returnStatus;
  }


  /*
   * Check some quantities that could cause code failures. If full
   * checking has been enabled, some of these checks may be
   * redundant.  In operations, we do not expect to turn on the full
   * checking (that is for LUT generation only).
   * 
   */

   /*
    * BAND30 is the index of band 31 if band 26 is not included.
    */

 for (i = 0; i < 2; i++) 
  {
    if (tables->RVS_BB_SV_Frame_No[i] < 0 ||
      tables->RVS_BB_SV_Frame_No[i] >= EV_1km_FRAMES)
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus,
         "Value(s) of Frame to use for BB or SV RVS "
	       " correction is out of range.",
         NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
      return returnStatus;
    }
  }
      
      
   /*
    * Check that each wavelength is positive (this could be done in a call to
    * Check_Valid_Range if the dead elements had a non-zero fill value).
    */

  for (i = 0; i < NUM_EMISSIVE_DETECTORS; i++) 
  {
    for (j = 0; j < tables->NUM_RSR_vs_Lambda[i]; j++) 
    {
      if (tables->wavelength[i][j] < TOLERANCE) 
      {
        returnStatus = MODIS_F_OUT_OF_RANGE;
        L1BErrorMsg(location, returnStatus,
                    "Value(s) of input RSR wavelength are zero.",
                    NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
        return returnStatus;
      }
    }
  }

  /*
   * Check the consistency of the BB and SV number of frames to use with the
   * first frame to use.
   */

  if (tables->BB_DN_first_frame_to_use < 0 ||
      tables->BB_DN_first_frame_to_use >= MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "BB DN first frame to use is out of range.",
                NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if (tables->BB_DN_number_of_frames_to_use <= 0 || 
      (tables->BB_DN_first_frame_to_use + 
       tables->BB_DN_number_of_frames_to_use) > MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "BB DN number of frames to use is out of "
                "range or is inconsistent\n"
                "with the value of DN_OBC first frame to use.",
                NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if (tables->SV_DN_first_frame_to_use < 0 ||
      tables->SV_DN_first_frame_to_use >= MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "SV DN first frame to use is out of range.",
                NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if (tables->SV_DN_number_of_frames_to_use <= 0 || 
      (tables->SV_DN_first_frame_to_use + 
       tables->SV_DN_number_of_frames_to_use) > MAX_1KM_OBC_FRAME_DIM)
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "SV DN number of frames to use is out of range or "
                "is inconsistent\n"
                "with the value of DN_OBC first frame to use.",
                NULL, EMISSIVE_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  /* 
   * Check if the LUT TEB_specified_uncertainty has zero or negative values.
   * The values should be all positive.
   */
   
   for (band = 0; band < NUM_EMISSIVE_BANDS; band++)
     if (tables->TEB_specified_uncertainty[band] <= TOLERANCE)
     {
       returnStatus = MODIS_F_OUT_OF_RANGE;
       L1BErrorMsg(location, returnStatus,
                   "Bad TEB_specified_uncertainty LUT (Zero "
                   "or negative values detected).",
                   NULL, EMISSIVE_TABLES_FILE, NULL, True);
       return returnStatus;
     }

  return(MODIS_S_OK);
}


PGSt_SMF_status Read_QA_Tables (L1A_granule_t *L1A_Gran, 
                                QA_tables_t   *QA_tables)
/*
!C***************************************************************
!Description:      Read in all QA lookup tables in this module.

!Input Parameters:
   L1A_granule_t  *L1A_Gran      contains satellite id and 
                                 EV start time

!Output Parameters:
   QA_tables_t    *QA_tables     Holds all QA tables.

!Revision History:

   Revision 02.32, October 29, 2000
   Add LUT "Control Parameters", Razor issue 142.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.31, December 22, 1999
   Redid the tables array to make it work for LUT generation also.
   Added Read_LUT_Tables to avoid a lot of code duplication.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.30, November 17, 1999
   Completely revised methodology of reading in tables.  Using a loop to
   allow for better error messages and for valid range checking (if enabled).
   Checks for PGE version consistency also added.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.10 Mar 1998
   Included all QA tables
   Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

   Revision 01.00 Jan 1998
   Initial development.
   Zhidong Hao (hao@gscmail.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 i, j, lun = QA_TABLES_FILE;
  int16 D_490, R, B, D, S_520, subsamp;
  int8  mask;
  char *location = "Read_QA_Tables";

  common_QA_tables_t  *common_QA_LUT = &QA_tables->common_QA_tables;
  refl_QA_tables_t    *refl_QA_LUT   = &QA_tables->refl_QA_tables;
  emiss_QA_tables_t   *emiss_QA_LUT  = &QA_tables->emiss_QA_tables;

  /*
   * Assign the array pointers to the data member in qa_luts.
   */ 

  ASSIGN_DATA_PTR(qa_luts,
                     PGE_VERSION_LUT_NAME,
                         common_QA_LUT->PGE_Version);

  ASSIGN_DATA_PTR(qa_luts,
                     MCST_VERSION_LUT_NAME,
                         common_QA_LUT->MCST_Version);

  ASSIGN_DATA_PTR(qa_luts,
                     QA_SERIAL_NUMBER_LUT_NAME,
                         common_QA_LUT->Serial_Number);

  ASSIGN_DATA_PTR(qa_luts,
                     PLATFORM_SHORT_NAME_LUT_NAME,
                         common_QA_LUT->AssociatedPlatformShortname);

  ASSIGN_DATA_PTR(qa_luts,
                     PACKAGE_ACCEPT_DATE_LUT_NAME,
                         common_QA_LUT->AlgorithmPackageAcceptanceDate);

  ASSIGN_DATA_PTR(qa_luts,
                     PACKAGE_MATURITY_CODE_LUT_NAME,
                         common_QA_LUT->AlgorithmPackageMaturityCode);

  ASSIGN_DATA_PTR(qa_luts,
                     DET_QUAL_FLAG_VALS_LUT_NAME,
                         common_QA_LUT->Detector_Quality_Flag_Values);

  ASSIGN_DATA_PTR(qa_luts,
                     DET_QUAL_FLAG2_VALS_LUT_NAME,
                         common_QA_LUT->Detector_Quality_Flag2_Values);

  ASSIGN_DATA_PTR(qa_luts,
                     MOON_OFFSET_LIMITS_LUT_NAME,
                         common_QA_LUT->moon_offset_limits);

  ASSIGN_DATA_PTR(qa_luts,
                     MISSION_PHASE_LUT_NAME,
                         common_QA_LUT->mission_phase);

  ASSIGN_DATA_PTR(qa_luts,
                      CONTROL_OPTIONS_LUT_NAME,
                          common_QA_LUT->control_options);

  ASSIGN_DATA_PTR(qa_luts,
                     BASE_VARI_VISUAL_FPA_LUT_NAME,
                         &refl_QA_LUT->var_visual_FPA);

  ASSIGN_DATA_PTR(qa_luts,
                     BASE_VARI_NIR_FPA_LUT_NAME,
                         &refl_QA_LUT->var_NIR_FPA);

  ASSIGN_DATA_PTR(qa_luts,
                     BB_TEMP_VARIANCE_LUT_NAME,
                         emiss_QA_LUT->var_T_bb);

  ASSIGN_DATA_PTR(qa_luts,
                     BB_AVG_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_bb_avg);

  ASSIGN_DATA_PTR(qa_luts,
                     LWIR_FPA_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_lwir);

  ASSIGN_DATA_PTR(qa_luts,
                     MWIR_FPA_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_mwir);

  ASSIGN_DATA_PTR(qa_luts,
                     MIR_SIDE_1_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_mir1);

  ASSIGN_DATA_PTR(qa_luts,
                     MIR_SIDE_2_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_mir2);

  ASSIGN_DATA_PTR(qa_luts,
                     MIR_AVG_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_mir_avg);

  ASSIGN_DATA_PTR(qa_luts,
                     INST_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_ins);

  ASSIGN_DATA_PTR(qa_luts,
                     CAVITY_TEMP_VAR_LUT_NAME,
                         &emiss_QA_LUT->var_T_cav);

  ASSIGN_DATA_PTR(qa_luts,
                     EMISS_NEdL_LUT_NAME,
                         emiss_QA_LUT->NEdL);

  ASSIGN_DATA_PTR(qa_luts,
                     CALIB_A1_LUT_NAME,
                         emiss_QA_LUT->a1);

  ASSIGN_DATA_PTR(qa_luts,
                     ROLL_THRESHOLD_LUT_NAME,
                         &common_QA_LUT->roll_threshold_angle);

  ASSIGN_DATA_PTR(qa_luts,
                     PITCH_THRESHOLD_LUT_NAME,
                         &common_QA_LUT->pitch_threshold_angle);

  ASSIGN_DATA_PTR(qa_luts,
                     YAW_THRESHOLD_LUT_NAME,
                         &common_QA_LUT->yaw_threshold_angle);

  /*
   * Read all tables.
   */

  returnStatus = Read_LUT_Tables(L1A_Gran, 
                                 QA_TABLES_FILE, 
                                 qa_luts);
  if (returnStatus != MODIS_S_OK) 
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Read_LUT_Tables", QA_TABLES_FILE, NULL, True);
    return returnStatus;
  }

  
/*
 *   Calculate Detector_Quality_Flag and Determine Dead Detector List
 *   and Noisy  Detector List. The LUT "Detector Quality Flag Values"
 *   contains the array  used to fill the L1B global attribute
 *   "Detector Quality Flag". Each array  element of this LUT (having
 *   a value of 0 or 1) sets one bit of the  corresponding L1B
 *   attribute. The second dimension of this LUT cycles through the 8
 *   bits of each word of the attribute, with array element [w][0]
 *   corresponding to the least significant bit w and array element
 *   [w][7] corresponding to the most significant bit. The global
 *   attributes Dead Detector List and Noisy Detector List used to be
 *   implemented from LUTs. These two attributes don't provide enough
 *   information for detector behavior and are superseded by Detector
 *   Quality Flag. For backward compatibility, they are retained. The
 *   array element [w][0] of the LUT "Detector Quality Flag Values" is
 *   noisy detector flag for detector w and [w][1] is the dead
 *   detector flag for detector w.
 */
  

  for (i = 0; i < NUM_DETECTORS; i++)
  {
    mask = 1;
    common_QA_LUT->Detector_Quality_Flag[i] = 0;
    for (j = 0; j < NUM_BITS_IN_UINT8; j++)
    {
      if (common_QA_LUT->Detector_Quality_Flag_Values[i][j] == 1)
        common_QA_LUT->Detector_Quality_Flag[i] |= mask;

      mask = mask << 1;
    }

    common_QA_LUT->noisy_detector[i] = 
        common_QA_LUT->Detector_Quality_Flag_Values[i][0];
    common_QA_LUT->dead_detector[i]  = 
        common_QA_LUT->Detector_Quality_Flag_Values[i][1];
  }

/*
 *   Calculate Detector_Quality_Flag2.
 *   The LUT "Detector Quality Flag2 Values" contains the array
 *   used to fill the L1B global attribute "Detector Quality Flag2".
 *   Each array  element of this LUT (having a value of 0 or 1)
 *   sets one bit of the corresponding L1B attribute. The second
 *   dimension of this LUT cycles through the 8 bits of each word
 *   of the attribute, with array element [w][0] corresponding to
 *   the least significant bit w and array element [w][7]
 *   corresponding to the most significant bit.
 *   The element [w][0] of the LUT "Detector Quality Flag2 Values"
 *   is for subframe 1 and [w][3] is for subframe 4 for detector w.
 *   bits 0 to 3 are used for noisy subframes, and bits 4 to 7 are
 *   used for dead subframes. Only 250m and 500m bands are considered.
 */

  D_490 = 0; 
  S_520 = 0; 
  for (R = 0; R <= INDEX_500M; R++)
  {
    for (B = 0; B < L1A_BANDS_AT_RES[R]; B++)
    {
      for (D = 0; D < DETECT_PER_BAND_AT_RES[R]; D++, D_490++)
      {
        mask = 1;
        common_QA_LUT->Detector_Quality_Flag2[D_490] = 0;
        for (j = 0; j < NUM_BITS_IN_UINT8; j++)
        {
          if (common_QA_LUT->Detector_Quality_Flag2_Values[D_490][j] == 1)
            common_QA_LUT->Detector_Quality_Flag2[D_490] |= mask;

          mask = mask << 1;
        }
        /*generate noisy and dead subframe list*/ 
        i=DETECT_PER_BAND_AT_RES[R]/DETECTORS_PER_1KM_BAND;
        for (subsamp = 0; subsamp < i; subsamp++, S_520++)
        { 
          common_QA_LUT->noisy_subframe[S_520] =
            common_QA_LUT->Detector_Quality_Flag2_Values[D_490][subsamp];
          common_QA_LUT->dead_subframe[S_520] =
            common_QA_LUT->Detector_Quality_Flag2_Values[D_490][subsamp+4];
        }
      }
    }
  }

 /*----------------------------------------------------------------------------
    Note: all the Emissive QA lookup table value should not be zero, because
          they will be divided by the calculated values in the production code. 
          Physically, all these values should be positive.
  -----------------------------------------------------------------------------*/
  for (i = 0; i < 12; i++)
    if ( emiss_QA_LUT->var_T_bb[i] < TOLERANCE )
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, 
                  "Table \"var_T_bb\" is zero or negative.",
                  NULL, QA_TABLES_FILE, invalidlutfile, True);
      return returnStatus;
    }
  
  if ( emiss_QA_LUT->var_T_bb_avg < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_bb_avg\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }
  if ( emiss_QA_LUT->var_T_lwir < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "Table \"var_T_lwir\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_mwir < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_mwir\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_mir1 < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_mir1\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_mir2 < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_mir2\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_mir_avg < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_mir_avg\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_ins < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_ins\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  if ( emiss_QA_LUT->var_T_cav < TOLERANCE )
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
                "Table \"var_T_cav\" is zero or negative.",
                NULL, QA_TABLES_FILE, invalidlutfile, True);
    return returnStatus;
  }

  for (i = 0; i < NUM_EMISSIVE_DETECTORS; i++)
    if ( emiss_QA_LUT->NEdL[i] < TOLERANCE )
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, 
                "Table \"NEdL\" is zero or negative.",
                  NULL, QA_TABLES_FILE, invalidlutfile, True);
      return returnStatus;
    }

  for (i = 0; i < NUM_EMISSIVE_DETECTORS; i++)
    if ( emiss_QA_LUT->a1[i] < TOLERANCE )
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, 
                "Table \"a1\" is zero or negative.",
                  NULL, QA_TABLES_FILE, invalidlutfile, True);
      return returnStatus;
    }

  /*
   * Check the platform name.  Currently, the allowed name for Terra is "AM-1"
   * and possibly "Terra", and the allowed name for Aqua could be "PM-1" or
   * "Aqua". These names should be consistent with the satellite id determined
   * from the metadata SHORTNAME in L1A granule. 
   */

  if (!((strcmp(common_QA_LUT->AssociatedPlatformShortname, "AM-1") == 0 ||
       strcmp(common_QA_LUT->AssociatedPlatformShortname, "Terra") == 0) &&
       L1A_Gran->satellite_id == TERRA) &&
      !((strcmp(common_QA_LUT->AssociatedPlatformShortname, "PM-1") == 0 ||
       strcmp(common_QA_LUT->AssociatedPlatformShortname, "Aqua") == 0) &&
       L1A_Gran->satellite_id == AQUA))
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "*** WRONG LUTS ARE INSTALLED ***\n"
                "Satellite instrument in the LUTs is inconsistent "
                "with the value\n"
                "inside the middle L1A granule (lun 500001).",
                NULL, QA_TABLES_FILE, 
                "This is probably due to improper installation of "
                "LUT files.", True);
    return returnStatus;
  }

  return(MODIS_S_OK);
}

static 
PGSt_SMF_status Read_LUT_Tables (L1A_granule_t    *L1A_Gran, 
                                 int32            lun, 
                                 LUT_Definition_t *luts)
/*
!C**************************************************************************
!Description:
   Read all lookup tables.  The "data" member of the luts structure will
   contain the data read in.  This should be already allocated (or assigned
   to the appropriate buffer) prior to calling this function.

!Input Parameters:
   L1A_granule_t    *L1A_Gran      contains satellite id and EV start time
   int32 lun                The logical unit number for the file to be read
   LUT_Definition_t *luts   Array of structures, one per LUT.  All data are
                            assigned and the "data" member of the luts
                            structure will contain the data read in.  This
                            should be already allocated (or assigned to the
                            appropriate buffer) prior to calling this function.

!Output Parameters:
   LUT_Definition_t *luts   "data" member filled with the data read from file.

!Revision History:
   Revision 01.00, December 22, 1999
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  PGSt_integer    version      = 1;
  char            file_name[PGSd_PC_FILE_PATH_MAX];
  int32           file_id      = 0;
  int32 i;
  int32 nluts;
  char *location = "Read_LUT_Tables";
  int32 i_bdsm = -1;                /* Index of LUT in refl_lut structure */
  int32 rank = 0;                   /* Intrinsic LUT rank */
  int32 *dims_ptr;                  /* Pointer to dimension arrays */ 

  /*
   * Initially, only check the bounds of the luts during LUT generation.  
   * Later, when bounds are stabilized, perhaps implement checking 
   * bounds operationally
   */

#ifdef ENABLE_LUT_VALID_RANGE_CHECKING
  int32 count, j;
#endif

  /*
   * Obtain the file name from the PCF and open the file for 
   * Science Data access.
   */

  returnStatus = PGS_PC_GetReference (lun, &version, file_name);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, 
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", lun, NULL, True);
    return returnStatus;
  }
  file_id = SDstart(file_name, DFACC_RDONLY);
  if (file_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, 
                "Could not open file for SD read access.",
                "SDstart", lun,
                "The file may be missing, corrupted or "
                "not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Determine the number of luts to read.
   */

  nluts = 0;
  while (luts[nluts].name) nluts++;

  /*
   * Read and check all LUTs.
   */

  for (i = 0; i < nluts; i++) 
  {
    if (!luts[i].data) 
    {
      char errmsg[PGS_SMF_MAX_MSGBUF_SIZE];
      returnStatus = MODIS_F_NOK;
      sprintf(errmsg, "\"data\" member not assigned for LUT %s", 
                    luts[i].name);
      L1BErrorMsg(location, returnStatus, errmsg,
                  NULL, lun, "*** CODE BUG ***", True);
      return returnStatus;
    }

    if (luts[i].kind == GLOBAL_ATTRIBUTE_LUT) 
      {
	returnStatus = read_attribute(file_id, luts[i].name,
				      luts[i].type, luts[i].data);
	if (returnStatus != MODIS_S_OK) 
	  {
	    L1BErrorMsg(location, returnStatus, NULL,
			"read_attribute", lun, NULL, True);
	    return returnStatus;
	  }
      }
    else 
      {

	/*
	 * If a Reflective BDSM LUT, set dimensions expected in HDF file
	 */

	rank     = luts[i].rank;
	dims_ptr = luts[i].dims;

	if ((i_bdsm = BDSM_index (luts[i].ascii_file)) >= 0)
	  {
	    rank = 1;
	    dims_ptr[0] = NUM_REFL_INDICES;
	    dims_ptr[1] = 0;
	    dims_ptr[2] = 0;
	    dims_ptr[3] = 0;
	  }

	/*
	 * Read the LUT
	 */
	returnStatus = Read_L1B_SDS_LUT(file_id, 
					luts[i].name,
					luts[i].type, 
					rank,
					dims_ptr, 
					L1A_Gran->data_collection_time,
					luts[i].data);
	if (returnStatus != MODIS_S_OK) 
	  {
	    L1BErrorMsg(location, returnStatus, NULL,
			"Read_L1B_SDS_LUT", lun, NULL, True);
	    return returnStatus;
	  }

	/*
	 * If a Reflective BDSM LUT, expand array into full BDSM
	 */

	if (i_bdsm >= 0)
	  {
	    void *data_new = NULL;
	    int32 n_bytes;
	    returnStatus = Expand_BDSM_LUT(luts[i].data, &data_new, 
					   luts[i].type, 1, &n_bytes);
	    if (returnStatus != MODIS_S_OK) 
	      {
		L1BErrorMsg(location, returnStatus, NULL,
			    "Expand_BDSM_LUT", lun, NULL, False);
		return returnStatus;
	      }
	    memcpy(luts[i].data, data_new, n_bytes);
	    free(data_new);
    	  }
      }

#ifndef NOCHECKLUT
    /*
     * Check the PGE version against the code macro.  If it does not
     * match, then don't bother reading any of the rest of the tables.
     */

    if (!strcmp(luts[i].name, PGE_VERSION_LUT_NAME)) 
    {
      if (strcmp(luts[i].data, PGE02_VERSION)) 
      {
        char msg[PGS_SMF_MAX_MSGBUF_SIZE];
        returnStatus = MODIS_F_OUT_OF_RANGE;
        sprintf(msg, "PGE version in the file (%s) does not "
                "match the code value of %s",
                (char *) luts[i].data, PGE02_VERSION);
        L1BErrorMsg(location, returnStatus, msg, NULL, lun,
                    "The LUT file is incompatible with this "
                    "version of the L1B code.",
                    True);
        return returnStatus;
      }
    }
#endif

#ifdef ENABLE_LUT_VALID_RANGE_CHECKING
    for (count = 1, j = 0; j < luts[i].rank; j++)
      count *= luts[i].dims[j];
    returnStatus = Check_Valid_Range(luts[i].name, 
                                     luts[i].type,
                                     luts[i].a_lb, 
                                     luts[i].a_ub,
                                     luts[i].a_fv, 
                                     count, 
                                     luts[i].data);
    if (returnStatus != MODIS_S_OK) 
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Check_Valid_Range", lun, NULL, True);
      return returnStatus;
    }
#endif
  }

  /*
   * Close the file.
   */

  if (SDend(file_id) == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not close access to the HDF file.",
                "SDend", lun, NULL, True);
    return returnStatus;
  }

  return(MODIS_S_OK);
}
 
/************************************************************************

Global variables which define the LUTs:

1. refl_luts     (defines LUTs for Reflective_Lookup_Tables_file)
2. emiss_luts    (defines LUTs for Emissive_Lookup_Tables_file)
3. qa_luts       (defines LUTs for QA_Lookup_Tables_file)

*************************************************************************/

/*
 * Define the reflective LUT information as LUT_Definition_t array.
 * This array must be NULL-terminated.
 */

LUT_Definition_t refl_luts[] = 
{
  { REFL_SERIAL_NUMBER_LUT_NAME, "refl_lut_serial_number.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_SERIAL_NUMBER_BUFFER, 0, 0, 0, 0},
    {"MAX_SERIAL_NUMBER_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PGE_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_PGE_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_PGE_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { MCST_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_MCST_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_MCST_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { K_INST_LUT_NAME,  "k_inst_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "-0.1", "0.1", "-999" },

  { K_FPA_LUT_NAME,  "k_fpa_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "-0.1", "0.1", "-999" },

  { M0_LUT_NAME,  "m0_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "0.0", "1.0", "-999" },

  { M1_LUT_NAME,  "m1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "1.0E-20", "1.0", "-999" },

  { RVS_RSB_LUT_NAME,  "rvs_rsb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, NUM_MIRROR_SIDES, 
     NUM_RSB_RVS_COEFFS, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "NUM_MIRROR_SIDES", "NUM_RSB_RVS_COEFFS", NULL},
    NULL, True, "-4.0E-4", "2.4", "-999" },

  /*
   * obsolete due to RSB UI algorithm update, 2/23/2011, Xu Geng
   */
  /*
  { SIGMA_RSB_ADC_LUT_NAME,  "sigma_rsb_adc_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", NULL, NULL, NULL},
    NULL, True, NULL, NULL, "-999" },

  { SIGMA_RVS_RSB_LUT_NAME,  "sigma_rvs_rsb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_REFLECTIVE_BANDS, NUM_MIRROR_SIDES, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", "NUM_MIRROR_SIDES", NULL, NULL, NULL},
    NULL, False, "0", "0.03", NULL },

  { SIGMA_M1_LUT_NAME,  "sigma_m1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "0.0", "1.0", "-999" },

  { SIGMA_K_INST_LUT_NAME,  "sigma_k_inst_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "0", "0.01", "-999" },

  { SIGMA_T_INST_LUT_NAME, "sigma_t_inst_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_PV_RESID_ELEC_LUT_NAME,  "sigma_pv_resid_elec_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, 0, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND",
     NULL, NULL},
    NULL, True, NULL, NULL, "-999" },

  { SIGMA_R_STAR_LIN_RESID_UCOEFF_LUT_NAME,  
    "sigma_r_star_lin_resid_ucoeff_table.asc",
    SDS_LUT, DFNT_FLOAT32, 5,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, NUM_4TH_ORDER_COEFFS},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", "NUM_4TH_ORDER_COEFFS"},
    NULL, True, NULL, NULL, "-999" },

  { RSB_NEDL_LUT_NAME,  "rsb_nedl_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "0", "1", "-999" },
  */

  { T_INST_REF_LUT_NAME, "t_inst_ref_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "250", "300", NULL },

  { T_FPA_REF_LUT_NAME,  "t_fpa_ref_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_FOCAL_PLANES, 0, 0, 0, 0},
    {"NUM_FOCAL_PLANES", NULL, NULL, NULL, NULL},
    NULL, False, "80", "300", NULL },

  { SWIR_CORRECTION_SWITCH_LUT_NAME, 
    "swir_oob_correction_switch_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { SWIR_CORR_SENDING_BAND_LUT_NAME, 
    "swir_oob_sending_band_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "20", "36", NULL },

  { SWIR_CORR_SENDING_DETECTOR_LUT_NAME,
    "swir_oob_sending_detector_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {DETECTORS_PER_1KM_BAND, 0, 0, 0, 0},
    {"DETECTORS_PER_1KM_BAND", NULL, NULL, NULL, NULL},
    NULL, True, "0", "9", NULL },

  { X_OOB_0_LUT_NAME,  "x_oob_0_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_SWIR_BANDS, MAX_DETECTORS_PER_SWIR_BAND,
     MAX_NUM_SWIR_SUBSAMPLES, NUM_MIRROR_SIDES, 0},
    {"NUM_SWIR_BANDS", "MAX_DETECTORS_PER_SWIR_BAND",
     "MAX_NUM_SWIR_SUBSAMPLES", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "-100", "100", "-999" },

  { X_OOB_1_LUT_NAME,  "x_oob_1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_SWIR_BANDS, MAX_DETECTORS_PER_SWIR_BAND,
     MAX_NUM_SWIR_SUBSAMPLES, NUM_MIRROR_SIDES, 0},
    {"NUM_SWIR_BANDS", "MAX_DETECTORS_PER_SWIR_BAND",
     "MAX_NUM_SWIR_SUBSAMPLES", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "-100", "100", "-999" },

  { X_OOB_2_LUT_NAME,  "x_oob_2_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_SWIR_BANDS, MAX_DETECTORS_PER_SWIR_BAND,
     MAX_NUM_SWIR_SUBSAMPLES, NUM_MIRROR_SIDES, 0},
    {"NUM_SWIR_BANDS", "MAX_DETECTORS_PER_SWIR_BAND",
     "MAX_NUM_SWIR_SUBSAMPLES", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "-100", "100", "-999" },

  { DN_OBC_1ST_FRAME_LUT_NAME, 
    "dn_obc_avg_first_frame_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "49", NULL },

  { DN_OBC_NUM_FRAMES_LUT_NAME, 
    "dn_obc_avg_number_of_frames_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "1", "50", NULL },

  { DN_STAR_MIN_LUT_NAME,  "dn_star_min_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_BANDS, 0, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "-40", "0", NULL },

  { DN_STAR_MAX_LUT_NAME,  "dn_star_max_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_BANDS, 0, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "4095", "4095", NULL },


  { RSB_SPECIFIED_UNCERTAINTY_LUT_NAME,  
    "rsb_specified_uncertainty_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_BANDS, 0, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { RSB_UI_SCALING_FACTOR_LUT_NAME,  
    "rsb_ui_scaling_factor_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_BANDS, 0, 0, 0, 0},
    {"NUM_REFLECTIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { E_SUN_OVER_PI_LUT_NAME,  "e_sun_over_pi_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_REFLECTIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { RSB_SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME,
    "rsb_sv_dn_moon_include_frames_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "50", NULL },

  { DN_SAT_EV_LUT_NAME,  "dn_sat_ev_table.asc",
    SDS_LUT, DFNT_FLOAT64, 4,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, 
     MAX_SAMPLES_PER_BAND, NUM_MIRROR_SIDES, 0},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND", 
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", NULL},
    NULL, True, "0", "4095", "-999" },

  { B26_B5_CORR_LUT_NAME,  "b26_b5_corr_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {DETECTORS_PER_1KM_BAND, 0, 0, 0, 0},
    {"DETECTORS_PER_1KM_BAND", NULL, NULL, NULL, NULL},
    NULL, True, "0.", "1.", "-999" },

  { B26_B5_CORR_SWITCH_LUT_NAME,  "b26_b5_corr_switch_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0},
    {"1", NULL, NULL, NULL, NULL},
    NULL, FALSE, "0", "1", NULL },
    
  { B26_B5_FRAME_OFFSET_LUT_NAME,  "b26_b5_frame_offset_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {DETECTORS_PER_1KM_BAND, 0, 0, 0, 0},
    {"DETECTORS_PER_1KM_BAND", NULL, NULL, NULL, NULL},
    NULL, True, "-10", "10", "-999" },

  { U1_LUT_NAME,  "u1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_REFLECTIVE_DETECTORS, 0,
     0, 0, 0},
    {"NUM_REFLECTIVE_DETECTORS", NULL,
     NULL, NULL, NULL},
    NULL, True, "0", NULL, NULL },

  { U2_LUT_NAME,  "u2_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_REFLECTIVE_DETECTORS,
     NUM_MIRROR_SIDES, NUM_U2_FRAME, 0, 0},
    {"NUM_REFLECTIVE_DETECTORS",
     "NUM_MIRROR_SIDES", "NUM_U2_FRAME", NULL, NULL},
    NULL, True, "0", NULL, NULL },

  { U3_LUT_NAME,  "u3_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_REFLECTIVE_DETECTORS,
     NUM_MIRROR_SIDES, 0, 0, 0},
    {"NUM_REFLECTIVE_DETECTORS",
     "NUM_MIRROR_SIDES", NULL, NULL, NULL},
    NULL, True, "0", NULL, NULL },

  { U4_LUT_NAME,  "u4_table.asc",
    SDS_LUT, DFNT_FLOAT32, 5,
    {NUM_REFLECTIVE_BANDS, MAX_DETECTORS_PER_BAND, MAX_SAMPLES_PER_BAND,
     NUM_MIRROR_SIDES, NUM_2ND_ORDER_COEFFS},
    {"NUM_REFLECTIVE_BANDS", "MAX_DETECTORS_PER_BAND",
     "MAX_SAMPLES_PER_BAND", "NUM_MIRROR_SIDES", "NUM_2ND_ORDER_COEFFS"},
    NULL, True, NULL, NULL, "-999.000" },

  { U2_FRAMES_LUT_NAME,  "u2_frames_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {NUM_U2_FRAME, 0, 0, 0, 0},
    {"NUM_U2_FRAME", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1353", NULL },

  { SWIR_UI_FACTOR_LUT_NAME,  "swir_ui_factor_table.asc",
    SDS_LUT, DFNT_FLOAT, 1,
    {NUM_SWIR_BANDS, 0, 0, 0, 0},
    {"NUM_SWIR_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { NULL }    /* Must be NULL-terminated */
};

/*
 * Define the emissive LUT information as LUT_Definition_t array.
 * This array must be NULL-terminated.
 */

LUT_Definition_t emiss_luts[] = 
{
  { EMISS_SERIAL_NUMBER_LUT_NAME, "emiss_lut_serial_number.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_SERIAL_NUMBER_BUFFER, 0, 0, 0, 0},
    {"MAX_SERIAL_NUMBER_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PGE_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_PGE_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_PGE_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { MCST_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_MCST_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_MCST_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { EPSILON_BB_LUT_NAME,  "epsilon_bb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "0.9", "1.1", NULL },

  { EPSILON_CAV_LUT_NAME,  "epsilon_cav_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "0.5", "1.0", NULL },

  { DELTA_T_BB_BETA_LUT_NAME,  "delta_t_bb_beta_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "-0.5", "0.5", NULL },

  { DELTA_T_BB_DELTA_LUT_NAME,  "delta_t_bb_delta_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "-0.5", "0.5", NULL },

  { PCX_TALK_LUT_NAME,  "pc_xt_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_PC_XT_BANDS, DETECTORS_PER_1KM_BAND, 
     NUM_PC_XT_PARAMETERS, 0, 0},
    {"NUM_PC_XT_BANDS", "DETECTORS_PER_1KM_BAND", 
     "NUM_PC_XT_PARAMETERS", NULL, NULL},
    NULL, True, "-15", "15", NULL },

  { RSR_LUT_NAME,  "rsr_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_EMISSIVE_DETECTORS, MAX_NUM_RSR_vs_LAMBDA, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", "MAX_NUM_RSR_vs_LAMBDA", 
     NULL, NULL, NULL},
     NULL, True, "0.01", "1.0", "0.0" },

  { WAVELENGTH_LUT_NAME,  "wavelength_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_EMISSIVE_DETECTORS, MAX_NUM_RSR_vs_LAMBDA, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", "MAX_NUM_RSR_vs_LAMBDA", 
     NULL, NULL, NULL},
    NULL, True, "3.0", "15.5", "-999." },

  { NUM_WL_INCREMENT_LUT_NAME,  "nwl_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "24", "67", NULL },

  { CALIB_A0_LUT_NAME,  "a0_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_a0_vs_T_inst_COEFF, NUM_MIRROR_SIDES, 
     NUM_EMISSIVE_DETECTORS, 0, 0},
    {"NUM_a0_vs_T_inst_COEFF", "NUM_MIRROR_SIDES", 
     "NUM_EMISSIVE_DETECTORS", NULL, NULL},
    NULL, True, "-100.", "300.", NULL },

  { CALIB_A2_LUT_NAME,  "a2_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_a2_vs_T_inst_COEFF, NUM_MIRROR_SIDES, 
     NUM_EMISSIVE_DETECTORS, 0, 0},
    {"NUM_a2_vs_T_inst_COEFF", "NUM_MIRROR_SIDES", 
     "NUM_EMISSIVE_DETECTORS", NULL, NULL},
    NULL, True, "-1.", "1.", NULL },

  /*
   * obsolete due to TEB UI algorithm update, 3/22/2011, Xu Geng
   */
  /*
  { UI_UCOEFF_LUT_NAME,  "ucoeff_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_EMISSIVE_DETECTORS, NUM_UI_PARAMETERS, 
     NUM_UI_POLYNOMIAL_COEFF, NUM_FI_POLYNOMIAL_COEFF, 0},
    {"NUM_EMISSIVE_DETECTORS", "NUM_UI_PARAMETERS", 
     "NUM_UI_POLYNOMIAL_COEFF", "NUM_FI_POLYNOMIAL_COEFF", NULL},
    NULL, True, NULL, NULL, "-999." },

  { SIGMA_TEB_PV_RESID_ELEC_LUT_NAME,  
   "sigma_teb_pv_resid_elec_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, NULL, NULL, "-999." },

  { SIGMA_TEB_ADC_LUT_NAME,  "sigma_teb_adc_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, NULL, NULL, "-999." },

  { UCOEFF_CALIBR_RESID_LUT_NAME,  
    "ucoeff_calibr_resid_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_EMISSIVE_DETECTORS, NUM_4TH_ORDER_COEFFS, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", "NUM_4TH_ORDER_COEFFS", 
    NULL, NULL, NULL},
    NULL, True, NULL, NULL, "-999." },

  { BAND_21_UNCERT_LSAT_LUT_NAME,  
    "band_21_uncert_lsat_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0", "0.2", NULL },
*/

  { BB_DN_1ST_FRAME_LUT_NAME,  
    "bb_dn_first_frame_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "49", NULL },

  { BB_DN_NUM_FRAMES_LUT_NAME,  
    "bb_dn_number_of_frames_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "1", "50", NULL },

  { SV_DN_1ST_FRAME_LUT_NAME,  
    "sv_dn_first_frame_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "49", NULL },

  { SV_DN_NUM_FRAMES_LUT_NAME,  
    "sv_dn_number_of_frames_to_use_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "1", "50", NULL },

  { SV_DN_MOON_INCLUDE_FRAMES_LUT_NAME, 
    "sv_dn_moon_include_frames_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "50", NULL },

  { OVERLAP_SCANS_B1_LUT_NAME,  "num_overlap_scans_b1_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "100", NULL },

  { PCX_CORRECTION_SWITCH_LUT_NAME,  
    "pcx_correction_switch_table.asc",
    SDS_LUT, DFNT_INT8, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { T_INS_FUNCTION_FLAG_LUT_NAME,  
    "t_ins_function_flag_table.asc",
    SDS_LUT, DFNT_INT32, 1,
    {NUM_T_INS_THERMISTORS, 0, 0, 0, 0},
    {"NUM_T_INS_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { T_INS_DEFAULT_LUT_NAME,  "t_ins_default_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "200.", "300.", NULL },

  { T_INS_OFFSET_LUT_NAME,  "t_ins_offset_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_T_INS_THERMISTORS, 0, 0, 0, 0},
    {"NUM_T_INS_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "-15.", "15.", NULL },

  { T_CAV_FUNCTION_FLAG_LUT_NAME,  
    "t_cav_function_flag_table.asc",
    SDS_LUT, DFNT_INT32, 1,
    {NUM_T_CAV_THERMISTORS, 0, 0, 0, 0},
    {"NUM_T_CAV_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { T_MIR_FUNCTION_FLAG_LUT_NAME,  
    "t_mir_function_flag_table.asc",
    SDS_LUT, DFNT_INT32, 1,
    {NUM_T_MIR_THERMISTORS, 0, 0, 0, 0},
    {"NUM_T_MIR_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { T_CAV_DEFAULT_LUT_NAME,  
    "t_cav_default_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "230", "300", NULL },

  { T_MIR_DEFAULT_LUT_NAME,  "t_mir_default_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "230", "300", NULL },

  { BB_WEIGHT_LUT_NAME,  "bb_weight_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_BB_THERMISTORS, 0, 0, 0, 0},
    {"NUM_BB_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { RVS_TEB_LUT_NAME,  "rvs_teb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_EMISSIVE_BANDS, DETECTORS_PER_1KM_BAND,
        NUM_MIRROR_SIDES, NUM_2ND_ORDER_COEFFS, 0},
    {"NUM_EMISSIVE_BANDS", "DETECTORS_PER_1KM_BAND",
        "NUM_MIRROR_SIDES", "NUM_2ND_ORDER_COEFFS", NULL},
    NULL, True, "-2.0E-4", "1.2", "-999" },

  { RVS_BB_SV_FRAME_NO_LUT_NAME, "rvs_bb_sv_frame_no_table.asc",
    SDS_LUT, DFNT_INT16, 1,
    {2, 0, 0, 0, 0},
    {"2", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1353", NULL },
  
  { BAND_21_B1_LUT_NAME,  "band_21_b1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {DETECTORS_PER_1KM_BAND, NUM_MIRROR_SIDES, 0, 0, 0},
    {"DETECTORS_PER_1KM_BAND", "NUM_MIRROR_SIDES", NULL, NULL, NULL},
    NULL, True, "0.005", "0.1", NULL },

  { L_MAX_LUT_NAME,  "l_max_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0}, 
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "110", NULL },

  { L_MIN_LUT_NAME,  "l_min_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0}, 
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "-10.", "0", NULL },

  { TEB_SPECIFIED_UNCERTAINTY_LUT_NAME,  
    "teb_specified_uncertainty_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0}, 
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { TEB_UI_SCALING_FACTOR_LUT_NAME,  
    "teb_ui_scaling_factor_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0}, 
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { BB_T_SAT_SWITCH_AQUA_LUT_NAME,  
    "bb_t_sat_switch_aqua_table.asc",
    SDS_LUT, DFNT_INT8, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { BB_T_SAT_AQUA_LUT_NAME,  
    "bb_t_sat_aqua_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_AQUA_BB_SAT_BANDS, 0, 0, 0, 0}, 
    {"NUM_AQUA_BB_SAT_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, "270.0", "315.0", NULL },

  { BB_T_SAT_DEFAULT_B1_BASELINE_AQUA_LUT_NAME,  
    "bb_t_sat_default_b1_baseline_aqua_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_AQUA_BB_SAT_BANDS, DETECTORS_PER_1KM_BAND, 
     NUM_MIRROR_SIDES, 0, 0}, 
    {"NUM_AQUA_BB_SAT_BANDS", "DETECTORS_PER_1KM_BAND", 
     "NUM_MIRROR_SIDES", NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { BB_T_SAT_DEFAULT_B1_C1_AQUA_LUT_NAME,
    "bb_t_sat_default_b1_c1_aqua_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_AQUA_BB_SAT_BANDS, DETECTORS_PER_1KM_BAND,
     NUM_MIRROR_SIDES, 0, 0},
    {"NUM_AQUA_BB_SAT_BANDS", "DETECTORS_PER_1KM_BAND",
     "NUM_MIRROR_SIDES", NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { BB_T_SAT_DEFAULT_B1_TLWIR_BASELINE_AQUA_LUT_NAME,
    "bb_t_sat_default_b1_Tlwir_baseline_aqua_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0},
    {"LWIR FPA Temperature", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_A0_LUT_NAME,
    "sigma_a0_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_a0_vs_T_inst_COEFF, NUM_MIRROR_SIDES,
     NUM_EMISSIVE_DETECTORS, 0, 0},
    {"NUM_a0_vs_T_inst_COEFF", "NUM_MIRROR_SIDES",
     "NUM_EMISSIVE_DETECTORS", NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { SIGMA_A2_LUT_NAME,
    "sigma_a2_table.asc",
    SDS_LUT, DFNT_FLOAT32, 3,
    {NUM_a2_vs_T_inst_COEFF, NUM_MIRROR_SIDES,
     NUM_EMISSIVE_DETECTORS, 0, 0},
    {"NUM_a2_vs_T_inst_COEFF", "NUM_MIRROR_SIDES",
     "NUM_EMISSIVE_DETECTORS", NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { SIGMA_RVS_EV_LUT_NAME,
    "sigma_RVS_ev_table.asc",
    SDS_LUT, DFNT_FLOAT32, 4,
    {NUM_EMISSIVE_BANDS, DETECTORS_PER_1KM_BAND,
     NUM_MIRROR_SIDES, NUM_2ND_ORDER_COEFFS, 0},
    {"NUM_EMISSIVE_BANDS", "DETECTORS_PER_1KM_BAND",
     "NUM_MIRROR_SIDES", "NUM_2ND_ORDER_COEFFS", NULL},
    NULL, True, NULL, NULL, NULL },

  { SIGMA_EPSILON_BB_LUT_NAME,
    "sigma_epsilon_bb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_EPSILON_CAV_LUT_NAME,
    "sigma_epsilon_cav_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_L_LAMBDA_LUT_NAME,
    "sigma_L_lambda_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_EMISSIVE_BANDS, NUM_1ST_ORDER_COEFFS, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", "NUM_1ST_ORDER_COEFFS", NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_L_TBB_LUT_NAME,
    "sigma_L_Tbb_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_L_TSM_LUT_NAME,
    "sigma_L_Tsm_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_L_TCAV_LUT_NAME,
    "sigma_L_Tcav_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_BANDS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { SIGMA_B1_BAND21_LUT_NAME,
    "sigma_b1_B21_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {DETECTORS_PER_1KM_BAND, NUM_MIRROR_SIDES, 0, 0, 0},
    {"DETECTORS_PER_1KM_BAND", "NUM_MIRROR_SIDES", NULL, NULL, NULL},
    NULL, True, NULL, NULL, NULL },

  { PCX_UI_FACTOR_LUT_NAME,
    "pcx_ui_factor_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_PC_XT_BANDS, 0, 0, 0, 0},
    {"NUM_PC_XT_BANDS", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { NULL }    /* Must be NULL-terminated */
};

/*
 * Define the quality assurance LUT information as LUT_Definition_t array.
 * This array must be NULL-terminated.
 */

LUT_Definition_t qa_luts[] = {

  { QA_SERIAL_NUMBER_LUT_NAME, "qa_lut_serial_number.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_SERIAL_NUMBER_BUFFER, 0, 0, 0, 0},
    {"MAX_SERIAL_NUMBER_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PGE_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_PGE_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_PGE_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { MCST_VERSION_LUT_NAME, NULL, GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_MCST_VERSION_BUFFER, 0, 0, 0, 0},
    {"MAX_MCST_VERSION_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PLATFORM_SHORT_NAME_LUT_NAME, 
    "associatedplatformshortname_table.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_ASSOCIATEDPLATFORMSHORTNAME_BUFFER, 0, 0, 0, 0},
    {"MAX_ASSOCIATEDPLATFORMSHORTNAME_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PACKAGE_ACCEPT_DATE_LUT_NAME, 
    "algorithmpackageacceptancedate_table.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_ALGORITHMPACKAGEACCEPTANCEDATE_BUFFER, 0, 0, 0, 0},
    {"MAX_ALGORITHMPACKAGEACCEPTANCEDATE_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { PACKAGE_MATURITY_CODE_LUT_NAME, 
    "algorithmpackagematuritycode_table.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_ALGORITHMPACKAGEMATURITYCODE_BUFFER, 0, 0, 0, 0},
    {"MAX_ALGORITHMPACKAGEMATURITYCODE_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { DET_QUAL_FLAG_VALS_LUT_NAME, 
    "detector_quality_flag_values_table.asc",
    SDS_LUT, DFNT_UINT8, 2,
    {NUM_DETECTORS, NUM_BITS_IN_UINT8, 0, 0, 0},
    {"NUM_DETECTORS", "NUM_BITS_IN_UINT8", NULL, NULL, NULL},
    NULL, True, "0", "1", NULL },

  { DET_QUAL_FLAG2_VALS_LUT_NAME,
    "detector_quality_flag2_values_table.asc",
    SDS_LUT, DFNT_UINT8, 2,
    {NUM_HIGH_RESOLUTION_DETECTORS, NUM_BITS_IN_UINT8, 0, 0, 0},
    {"NUM_HIGH_RESOLUTION_DETECTORS", "NUM_BITS_IN_UINT8", NULL, NULL, NULL},
    NULL, True, "0", "1", NULL },

  { MISSION_PHASE_LUT_NAME, "mission_phase_table.asc",
    GLOBAL_ATTRIBUTE_LUT, DFNT_CHAR8, 1,
    {MAX_MISSION_PHASE_BUFFER, 0, 0, 0, 0},
    {"MAX_MISSION_PHASE_BUFFER", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { CONTROL_OPTIONS_LUT_NAME, "control_options_table.asc",
    SDS_LUT, DFNT_UINT8, 1,
    {NUM_CONTROL_OPTIONS, 0, 0, 0, 0},
    {"NUM_CONTROL_OPTIONS", NULL, NULL, NULL, NULL},
    NULL, False, "0", "1", NULL },

  { BASE_VARI_VISUAL_FPA_LUT_NAME, 
    "visual_fpa_base_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { BASE_VARI_NIR_FPA_LUT_NAME, 
    "nir_fpa_base_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, NULL, NULL, NULL },

  { BB_TEMP_VARIANCE_LUT_NAME, "t_bb_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_BB_THERMISTORS, 0, 0, 0, 0},
    {"NUM_BB_THERMISTORS", NULL, NULL, NULL, NULL},
    NULL, False, "0.0002", "0.1", NULL },

  { BB_AVG_TEMP_VAR_LUT_NAME, 
    "bb_average_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.00003", "0.1", NULL },

  { LWIR_FPA_TEMP_VAR_LUT_NAME, 
    "lwir_fpa_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.000001", "0.1", NULL },

  { MWIR_FPA_TEMP_VAR_LUT_NAME, 
    "mwir_fpa_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0005", "0.1", NULL },

  { MIR_SIDE_1_TEMP_VAR_LUT_NAME, 
    "mirrorside_1_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0001", "0.1", NULL },

  { MIR_SIDE_2_TEMP_VAR_LUT_NAME, 
    "mirrorside_2_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0001", "0.1", NULL },

  { MIR_AVG_TEMP_VAR_LUT_NAME, 
    "mirror_average_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.001", "0.1", NULL },

  { INST_TEMP_VAR_LUT_NAME, 
    "instrument_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.001", "0.1", NULL },

  { CAVITY_TEMP_VAR_LUT_NAME, 
    "cavity_temperature_variance_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0}, {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.01", "0.1", NULL },

  { EMISS_NEdL_LUT_NAME, "nedl_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "0.0002", "0.1", NULL },

  { CALIB_A1_LUT_NAME, "a1_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {NUM_EMISSIVE_DETECTORS, 0, 0, 0, 0},
    {"NUM_EMISSIVE_DETECTORS", NULL, NULL, NULL, NULL},
    NULL, True, "0.0002", "0.1", NULL },

  { MOON_OFFSET_LIMITS_LUT_NAME, "moon_offset_limits_table.asc",
    SDS_LUT, DFNT_FLOAT32, 2,
    {NUM_BANDS, NUM_MOON_OFFSET_LIMITS, 0, 0, 0},
    {"NUM_BANDS", "NUM_MOON_OFFSET_LIMITS", NULL, NULL, NULL},
    NULL, False, "-200", "200", NULL },

  { ROLL_THRESHOLD_LUT_NAME, "roll_threshold_angle_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0},
    {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0", "360.0", NULL },

  { PITCH_THRESHOLD_LUT_NAME, "pitch_threshold_angle_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0},
    {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0", "360.0", NULL },

  { YAW_THRESHOLD_LUT_NAME, "yaw_threshold_angle_table.asc",
    SDS_LUT, DFNT_FLOAT32, 1,
    {1, 0, 0, 0, 0},
    {"1", NULL, NULL, NULL, NULL},
    NULL, False, "0.0", "360.0", NULL },

  { NULL }    /* Must be NULL-terminated */
};


PGSt_SMF_status Read_L1B_SDS_LUT(int32   sd_id, 
                                 char    *name, 
                                 int32   data_type, 
                                 int32   rank,
                                 int32   *dims, 
                                 float64 data_collection_TAI_time,
                                 void    *data)
/*
!C**************************************************************************
!Description: This function reads one level 1B LUT implemented as SDS. 

!Input Parameters:
   int32  sd_id       HDF file ID of the LUT file being read.
   char   *name       Name of the LUT as implemented in the HDF file
   int32  data_type   HDF data type of the data being read.  This also
                      matches the data type of the memory addressed by
                      the "data" pointer.
   int32  rank        Array rank of the data.  This matches the intrinsic
                      array rank of the data in the HDF file.
   int32  *dims       Array containing the size of each dimension of the
                      intrinsic array.
   float64 data_collection_TAI_time  The data collection time in TAI
                                     seconds (seconds since 1/1/1993)

!Output Parameters:
   void   *data       Address of memory being filled by this function.
                      The data type and number of elements of this array
                      should be consistent with the inputs.

!Revision History:
   Revision 01.00 June 23, 2000
   Initial development.
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK; 
  int32 algorithm;
  int32 start[MAX_VAR_DIMS];
  int32 i;
  char errmsgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  char *location = "Read_L1B_SDS_LUT";

      /* Check input parameters for gross validility */

  if (sd_id == FAIL || !name || rank <= 0 || rank > MAX_VAR_DIMS || 
      !dims || !data || data_collection_TAI_time <= 0) 
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus,
                "One or more input parameters are invalid.\n",
                NULL, 0, "Code bug.", False);
    return returnStatus;
  }

      /* Determine the time dependent algorithm for this LUT */

  algorithm = TDLUT_GetAlgorithm(sd_id, name);

  if (algorithm == TDLUT_CONSTANT) 
  { /* constant LUT */
    for (i = 0; i < MAX_VAR_DIMS; i++)
      start[i] = 0;

    returnStatus = read_sds_rankn(sd_id, name, data_type, rank,
                                  start, dims, data);
    if (returnStatus != MODIS_S_OK) 
    {
      L1BErrorMsg(location, returnStatus, NULL,
                    "read_sds_rankn", 0, NULL, False);
    }
  }
  else if (algorithm == TDLUT_STEPFUNCTION) 
  { /* step function LUT */
    returnStatus = TDLUT_ReadStepFunction(sd_id, name, data_type,
                                          rank, dims, 
                                          data_collection_TAI_time,
                                          data);
    if (returnStatus != MODIS_S_OK) 
    {
      L1BErrorMsg(location, returnStatus, NULL,
                    "TDLUT_ReadStepFunction", 0, NULL, False);
    }
  }
  else if (algorithm == TDLUT_PIECEWISE_LINEAR) 
  { /* piecewise linear function LUT */
    returnStatus = TDLUT_ReadPiecewiseLinearFunction(sd_id, name, data_type,
                                          rank, dims, 
                                          data_collection_TAI_time,
                                          data);
    if (returnStatus != MODIS_S_OK) 
    {
      L1BErrorMsg(location, returnStatus, NULL,
                    "TDLUT_ReadPiecewiseLinearFunction", 0, NULL, False);
    }
  }

    else 
    { /* algorithm flag is invalid */
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "The algorithm flag is invalid for LUT \"%s\".",
            name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
  }

  return returnStatus;
}

int32 TDLUT_GetAlgorithm(int32 sd_id, 
                         char *name)
/*
!C**************************************************************************
!Description: Given a HDF file ID and a LUT name (LUT is assumed to be an SDS),
              return the time-dependence algorithm flag associated with the LUT.

!Input Parameters:
   int32  sd_id       HDF file ID of the LUT file being read.
   char   *name       Name of the LUT as implemented in the HDF file

!Output Parameters:
   (none)

!Revision History:
   Revision 01.00 June 22, 2000
   Initial development.
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   1.  The algorithm is implemented as an SDS attribute, of int32 type
       (scalar).  The name of the attribute is given by the macro
       TDLUT_ALGORITHM_ATTR_NAME.  If this attribute is missing, then
       it is assumed that the algorithm is TDLUT_CONSTANT.

   2.  Initialize the SDS identifier to FAIL.
!END********************************************************************
*/
{
  int32 algorithm_flag;
  int32 sds_id = FAIL;             /* SDS identifier */
  int32 sds_index;                 /* SDS index */ 
  int32 attr_index;                /* index of an attribute of the SDS */
  int32 nattrs;                    /* number of attributes of the SDS */
  int32 count;                     /* number of values of an attribute */
  int32 rank;                      /* rank of the SDS */
  int32 dims[MAX_VAR_DIMS];        /* dimensionality of the SDS */
  int32 data_type;                 /* data type of an attributes or the SDS */
  char buffer[MAX_NC_NAME];        /* buffer to store a name */
  char errmsgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  intn hdf_return;
  char *location = "TDLUT_GetAlgorithm";

    /* Check the gross validility of the input parameters */
 
  if (sd_id == FAIL || !name) 
  {
    L1BErrorMsg(location, MODIS_F_INVALID_ARGUMENT, 
                "Input parameter sd_id is -1 or name is NULL\n",
                NULL, 0, "Code bug.", False);     
    return TDLUT_INVALID;
  }

    /* Find the index of the SDS from its name */

  sds_index = SDnametoindex(sd_id, name);
  if (sds_index == FAIL) 
  {
    sprintf(errmsgbuf, 
            "Could not find SDS \"%s\" in the file.", name);
    L1BErrorMsg(location, MODIS_F_READ_ERROR, errmsgbuf, 
                "SDnametoindex", 0, invalidinputfile, False);
    return TDLUT_INVALID;
  }

    /* Get the identifier of the SDS */
 
  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) 
  {
    sprintf(errmsgbuf, 
            "Could not open access to SDS \"%s\".", name);
    L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    return TDLUT_INVALID;
  }

    /* Get the information of the SDS, especially the attribute information */

  hdf_return = SDgetinfo(sds_id, 
                         buffer, 
                         &rank, 
                         dims, 
                         &data_type, 
                         &nattrs);
  if (hdf_return == FAIL) 
  {
    algorithm_flag = TDLUT_INVALID;
    sprintf(errmsgbuf, 
            "Could not get information about SDS \"%s\".", name);
    L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, "SDgetinfo", 0,
                corruptinputfile, False);
    goto TDLUT_GetAlgorithm_exit; 
  }

  algorithm_flag = TDLUT_CONSTANT;

    /* If there is no attribute for the SDS, assume the algorithm is
     * TDLUT_CONSTANT.
     */
 
  if (nattrs == 0)
    goto TDLUT_GetAlgorithm_exit;

    /* Loop through the attributes to find the attribute 
     * TDLUT_ALGORITHM_ATTR_NAME. If it is missing, the algorithm is
     * assumed to be TDLUT_CONSTANT. If an error occurs, return 
     * TDLUT_INVALID.
     */

  for (attr_index = 0; attr_index < nattrs; attr_index++) 
  {
    hdf_return = SDattrinfo(sds_id, 
                            attr_index, 
                            buffer, 
                            &data_type, 
                            &count);
    if (hdf_return == FAIL) 
    {
      algorithm_flag = TDLUT_INVALID;
      sprintf(errmsgbuf, 
              "Could not get information of an attribute of SDS \"%s\""
              "with a valid index.", name);
      L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, "SDattrinfo", 0,
                  corruptinputfile, False);
      goto TDLUT_GetAlgorithm_exit;
    }
    if (strcmp(buffer, TDLUT_ALGORITHM_ATTR_NAME) == 0) 
    {
      if (data_type != DFNT_INT32 || count != 1) 
      {
        algorithm_flag = TDLUT_INVALID;
        sprintf(errmsgbuf, 
                "The data type or the count is not correct for the"
                "attribute \"%s\" of SDS \"%s\".", TDLUT_ALGORITHM_ATTR_NAME,
                name);
        L1BErrorMsg(location, MODIS_F_NOK, errmsgbuf, NULL, 0, 
                    invalidinputfile, False);
        goto TDLUT_GetAlgorithm_exit;
      } 
      else 
      {
        hdf_return = SDreadattr(sds_id, attr_index, &algorithm_flag);
        if (hdf_return == FAIL) 
        {
          algorithm_flag = TDLUT_INVALID;
          sprintf(errmsgbuf, 
                  "Failed to read attribute \"%s\" of SDS \"%s\".",
                  TDLUT_ALGORITHM_ATTR_NAME, name);
          L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, 
                      "SDreadattr", 0, corruptinputfile, False);
          goto TDLUT_GetAlgorithm_exit;
        }  
      } 
      break;   /* out of the "for" loop */
    }
  }

TDLUT_GetAlgorithm_exit:

  if (sds_id != FAIL) SDendaccess(sds_id);
  
  return algorithm_flag;
} 
     
PGSt_SMF_status TDLUT_ReadStepFunction (int32   sd_id,
                                        char    *name,
                                        int32   data_type,
                                        int32   rank,
                                        int32   *dims,
                                        float64 data_collection_TAI_time, 
                                        void    *data)
/*
!C**************************************************************************
!Description: This function reads one level 1B LUT which has a step-function 
              time dependence.
 
!Input Parameters:
   int32  sd_id       HDF file ID of the LUT file being read.
   char   *name       Name of the LUT as implemented in the HDF file
   int32  data_type   HDF data type of the data being read.  This also
                      matches the data type of the memory addressed by
                      the "data" pointer.
   int32  rank        Array rank of the data.  This matches the intrinsic
                      array rank of the data in the HDF file.
   int32  *dims       Array containing the size of each dimension of the
                      intrinsic array.
   float64 data_collection_TAI_time  The data collection time in TAI
                                     seconds (seconds since 1/1/1993)

!Output Parameters:
   void   *data       Address of memory being filled by this function.
                      The data type and number of elements of this array
                      should be consistent with the inputs.

!Revision History:
   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00 June 23, 2000
   Initial development.
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The actual LUT implemented in the HDF LUT file has an extra leading
   dimension from the intrinsic rank and dimensions indicated by the
   input variable.  The extra leading dimension corresponds to TAI
   time of the steps.

   There must be an additional SDS attribute having the name given by the
   macro TDLUT_STEPFUNCTION_TIMES.  This attribute must be of type
   float64.  The number of elements of this attribute must match the
   size of the leading dimension of the LUT as implemented in the HDF file.
   The times in this attribute are assumed monotonically increasing.
!END********************************************************************
*/
{
#define TDLUT_STEPFUNCTION_MAX_TIMES 1000
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 sds_id = FAIL;             /* SDS identifier */
  int32 sds_index;                 /* SDS index */
  int32 attr_index;                /* index of an attribute of the SDS */
  int32 nattrs;                    /* number of attributes of the SDS */
  int32 count;                     /* number of values of an attribute */
  int32 f_rank;                    /* rank of the SDS */
  int32 full_dims[MAX_VAR_DIMS];   /* leading dimension + intrinsic dims */
  int32 f_data_type;               /* data type of an attributes or the SDS */
  char buffer[MAX_NC_NAME];        /* buffer to store a name */
  char errmsgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  intn hdf_return;
  int32 start[MAX_VAR_DIMS];
  int32 edge [MAX_VAR_DIMS];
  float64 times[TDLUT_STEPFUNCTION_MAX_TIMES];
  int32 i;
  char *location = "TDLUT_ReadStepFunction";

    /* Check input parameters for gross validity */

  if (sd_id == FAIL || !name || rank <= 0 || rank >= MAX_VAR_DIMS || 
      !dims || !data || data_collection_TAI_time <= 0) 
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus,
                "One or more input parameters are invalid.\n",
                NULL, 0, "Code bug.", False);
    return returnStatus;
  }

    /* Find the index of the SDS from its name */

  sds_index = SDnametoindex(sd_id, name);
  if (sds_index == FAIL) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "Could not find SDS \"%s\" in the file.", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }
  
    /* Open SDS access to the data set */

  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not open access to SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }
                     
    /* Get information about the SDS */

  hdf_return = SDgetinfo(sds_id, buffer, &f_rank, full_dims, 
                         &f_data_type, &nattrs);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not get information about SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDgetinfo", 0,
                corruptinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }
  if (data_type != f_data_type) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "SDS \"%s\" has data type mismatch.\n"
            "File data type = %ld, expected type = %ld\n",
            name, (long)f_data_type, (long)data_type);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

  /*  Check that the HDF rank is one more than the actual rank. */  
  if ((rank + 1) != f_rank) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "SDS \"%s\" has a rank mismatch.\n"
            "HDF File rank = %ld, expected rank = %ld\n",
            name, (long)f_rank, (long)(rank + 1));
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

    /* For step function time-dependent LUT, number of attributes
     * should not be zero. The attributes include at least "algorithm"
     * and "times".
     */

  if (nattrs == 0) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "The number of attributes for the step function LUT"
            " \"%s\" is 0.", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

    /* Find the attribute TDLUT_STEPFUNCTION_TIMES */

  attr_index = SDfindattr(sds_id, TDLUT_STEPFUNCTION_TIMES);
  if (attr_index == FAIL) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "Step function lut \"%s\" does not have attribute \"%s\".",
            name, TDLUT_STEPFUNCTION_TIMES);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDfindattr", 0,
                invalidinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }
 
    /* Get the attribute info */
 
  hdf_return = SDattrinfo(sds_id, attr_index, buffer, &f_data_type, &count);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not get info of the attribute \"%s\" of SDS \"%s\".",
            TDLUT_STEPFUNCTION_TIMES, name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDattrinfo", 0,
                corruptinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

  /*
   * Check the data type and count of the attribute
   * TDLUT_STEPFUNCTION_TIMES. The data type should be float64
   * (DFNT_FLOAT64) and the count should not exceed
   * TDLUT_STEPFUNCTION_MAX_TIMES.
   */


  if (f_data_type != DFNT_FLOAT64 || 
      count > TDLUT_STEPFUNCTION_MAX_TIMES) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "Either the data type is not float64 or the number of\n"
            "values are too large (>500) for the "
            "attribute\n \"%s\" of SDS \"%s\".",
            TDLUT_STEPFUNCTION_TIMES, name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
    goto TDLUT_ReadStepFunction_exit;
  }

    /* Read the attribute TDLUT_STEPFUNCTION_TIMES */

  hdf_return = SDreadattr(sds_id, attr_index, times);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Failed to read attribute \"%s\" of SDS \"%s\".",
            TDLUT_STEPFUNCTION_TIMES, name);
    L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, "SDreadattr", 0,
                corruptinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

    /* Determine the index of the leading dimension to read the data by
     * comparing the data collection time with the times as read from the
     * attribute. 
     */

  if (data_collection_TAI_time < times[0])
    start[0] = 0;
  else if (data_collection_TAI_time >= times[count - 1])
    start[0] = count - 1;
  else 
  {
    for (i = 0; i < count - 1; i++) 
    {
      if (data_collection_TAI_time >= times[i] && 
          data_collection_TAI_time < times[i + 1]) 
      {
        start[0] = i;
        break;
      }
    }
  } 
    
  edge[0] = 1;

    /* Set values for other elements of "start" and "edge" arrays */

  for (i = 1; i < MAX_VAR_DIMS; i++) 
  {
    start[i] = 0; 
    if (i <= rank) 
      edge[i] = dims[i - 1];
    else
      edge[i] = 0;
  }

    /* Read the data from the LUT */

  hdf_return = SDreaddata(sds_id, start, NULL, edge, data);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not read data from SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    goto TDLUT_ReadStepFunction_exit;
  }

TDLUT_ReadStepFunction_exit:

  if (sds_id != FAIL) SDendaccess(sds_id);

  return returnStatus;
}   
   
/*=============================================== */


PGSt_SMF_status TDLUT_ReadPiecewiseLinearFunction 
                                       (int32   sd_id, 
                                        char    *name, 
                                        int32   data_type, 
                                        int32   rank, 
                                        int32   *dims, 
                                        float64 data_collection_TAI_time, 
                                        void    *data)
/*
!C**************************************************************************
!Description: This function reads one level 1B LUT which is piecewise   
              linear time dependent.
 
!Input Parameters:
   int32  sd_id       HDF file ID of the LUT file being read.
   char   *name       Name of the LUT as implemented in the HDF file
   int32  data_type   HDF data type of the data being read.  This also
                      matches the data type of the memory addressed by
                      the "data" pointer.
   int32  rank        Array rank of the data.  This matches the intrinsic
                      array rank of the data in the HDF file.
   int32  *dims       Array containing the size of each dimension of the
                      intrinsic array.
   float64 data_collection_TAI_time  The data collection time in TAI
                                     seconds (seconds since 1/1/1993)

!Output Parameters:
   void   *data       Address of memory being filled by this function.
                      The data type and number of elements of this array
                      should be consistent with the inputs.

!Revision History:
   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00 March 04, 2001
   Initial development.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The actual LUT implemented in the HDF LUT file has an extra leading
   dimension from the intrinsic rank and dimensions indicated by the
   input variable.  The extra leading dimension corresponds to the TAI
   times of the LUT tables which can be interpolated.
   
   The LUT which is being interpolated must have data type float32.

   There must be an additional SDS attribute having the name given by the
   macro TDLUT_PIECEWISE_LINEAR_TIMES.  This attribute must be of type
   float64.  The number of elements of this attribute must match the
   size of the leading dimension of the LUT as implemented in the HDF file.
!END*********************************************************************/

{

/********************* TDLUT_PIECEWISE_LINEAR_MAX_TIMES ******************/
/*  Maximum number of separate time dependent lookup tables allowed in a */
/*        file.                                                          */
#define TDLUT_PIECEWISE_LINEAR_MAX_TIMES 1000
   
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int32 sds_id = FAIL;             /* SDS identifier                    */
  int32 sds_index;                 /* SDS index                          */
  int32 attr_index;                /* index of an attribute of the SDS   */
  int32 nattrs;                    /* number of attributes of the SDS    */
  int32 num_lut_copies;            /* number of LUT copies in the file   */
  int32 f_rank;                    /* rank of the SDS                    */
  int32 full_dims[MAX_VAR_DIMS];   /* leading dimension + intrinsic dims */
  int32 f_data_type;               /* data type of an attribute or SDS   */
  char  buffer[MAX_NC_NAME];       /* buffer to store a name             */
  char  errmsgbuf[PGS_SMF_MAX_MSGBUF_SIZE];
  intn  hdf_return;
  int32 start[MAX_VAR_DIMS];
  int32 next_start[MAX_VAR_DIMS];
  int32 edge [MAX_VAR_DIMS];
  float64 times[TDLUT_PIECEWISE_LINEAR_MAX_TIMES];
  float64 data_coll_time;          /* Scaled data collection time         */
  float64 lut_data_time;           /* Scaled data time of first LUT       */
  float64 next_lut_data_time;      /* Scaled data time of second LUT      */
  float64 slope_denom;             /* Denominator of Linear Slope         */
  float64 xproportion;             /* Linear Proportion (x - x1)/(x2 - x1)*/
  int32 num_lut_elements;          /* Total number of elements in one LUT */
  char *location = "TDLUT_ReadPiecewiseLinearFunction";
  void *data2 = NULL;              /* Pointer to data in second LUT       */
  int32 pw_size = sizeof(float32); /* Allowed size of PW Linear LUT data  */
  int32 i;

  /* Check input parameters for gross validity                            */

  if (sd_id == FAIL || !name || rank <= 0 || rank >= MAX_VAR_DIMS || 
      !dims || !data || data_collection_TAI_time <= 0) 
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus,
                "One or more input parameters are invalid.\n",
                NULL, 0, "Code bug.", False);
    return returnStatus;
  }

 /*  Determine the input data type of the LUT.  It must be float32.       */
 if (data_type != DFNT_FLOAT32) 
 {
   returnStatus = MODIS_F_NOK;
   sprintf(errmsgbuf,
       "\"%s\" must have data type float32 to be Piecewise Linear.", name) ;
   L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
   return returnStatus;
 }
 
  /* Find the index of the SDS from its name                              */

  sds_index = SDnametoindex(sd_id, name);
  if (sds_index == FAIL) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "Could not find SDS \"%s\" in the file.", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDnametoindex", 0,
                invalidinputfile, False);
    return returnStatus;
  }
  
  /* Open SDS access to the data set                                      */

  sds_id = SDselect (sd_id, sds_index);
  if (sds_id == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not open access to SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDselect", 0,
                corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }
                     
  /* Get information about the SDS                                        */

  hdf_return = SDgetinfo(sds_id, buffer, &f_rank, full_dims, 
                         &f_data_type, &nattrs);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not get information about SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDgetinfo", 0,
                corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }
  

  /*  Check that the HDF data type is the same as that expected.          */  
  if (data_type != f_data_type) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "SDS \"%s\" has data type mismatch.\n"
            "HDF File data type = %ld, expected type = %ld\n",
            name, (long)f_data_type, (long)data_type);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /*  Check that the HDF rank is one more than the actual rank.           */  
  if ((rank + 1) != f_rank) 
  {
    returnStatus = MODIS_F_READ_ERROR;
    sprintf(errmsgbuf, 
            "SDS \"%s\" has a rank mismatch.\n"
            "HDF File rank = %ld, expected rank = %ld\n",
            name, (long)f_rank, (long)(rank + 1));
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /*  For piecewise linear time-dependent LUTs, the number of attributes  */
  /*  should be two or more. (The attributes must include at least        */
  /*  "algorithm" and "times".)                                           */
     

  if (nattrs < 2) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
       "There must be at least 2 attributes for a Piecewise Linear LUT but"
            " \"%s\" has %ld.", name, (long)nattrs);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0,
                invalidinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;  
  }

    /* Find the attribute TDLUT_PIECEWISE_LINEAR_TIMES                    */

  attr_index = SDfindattr(sds_id, TDLUT_PIECEWISE_LINEAR_TIMES);
  if (attr_index == FAIL) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "Piecewise linear function LUT \"%s\" does not have attribute \"%s\".",
            name, TDLUT_PIECEWISE_LINEAR_TIMES);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDfindattr", 0,
                invalidinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }
 
    /* Get the attribute info                                             */
 
  hdf_return = SDattrinfo(sds_id, 
                          attr_index, 
                          buffer, 
                          &f_data_type, 
                          &num_lut_copies);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not get info of the attribute \"%s\" of SDS \"%s\".",
            TDLUT_PIECEWISE_LINEAR_TIMES, name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDattrinfo", 0,
                corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /* Check the data type and count of the attribute                       */
  /* TDLUT_PIECEWISE_LINEAR_TIMES.  The data type should be float64       */
  /* (DFNT_FLOAT64) and the count should not exceed                       */
  /* TDLUT_PIECEWISE_LINEAR_MAX_TIMES.                                    */

  if (f_data_type != DFNT_FLOAT64 ||
      num_lut_copies > TDLUT_PIECEWISE_LINEAR_MAX_TIMES) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf,
            "Either the data type is not float64 or the number of\n"
            "values are too large (>500) for the "
            "attribute\n \"%s\" of SDS  \"%s\".",
             TDLUT_PIECEWISE_LINEAR_TIMES, name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /* Read the attribute TDLUT_PIECEWISE_LINEAR_TIMES                      */

  hdf_return = SDreadattr(sds_id, attr_index, times);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Failed to read attribute \"%s\" of SDS \"%s\".",
            TDLUT_PIECEWISE_LINEAR_TIMES, name);
    L1BErrorMsg(location, MODIS_F_HDF_ERROR, errmsgbuf, "SDreadattr", 0,
                corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /* Determine the index of the leading dimension to read the data by     */
  /* comparing the data collection time with the times as read from the   */
  /* attribute.                                                           */
 

  /*  First special case:  even though the LUT is designated piecewise    */
  /*    linear, there is actually only one table in it.                   */
  if (num_lut_copies == 1) 
  {
    start[0] = 0;
  }

  /*  Second special case:  granule time before time of first LUT         */
  else if (data_collection_TAI_time < times[0]) 
  {
    start[0]      = 0;
    next_start[0] = 1;
  }
  /*  Third special case:  granule time after time of last LUT            */
  else if (data_collection_TAI_time >= times[num_lut_copies-1]) 
  {
    start[0]      = num_lut_copies - 2;
    next_start[0] = num_lut_copies - 1;
  }
  /*  All other cases                                                     */
  else 
  {
    for (i = 0; i < num_lut_copies - 1; i++) 
    {
      if (data_collection_TAI_time >= times[i] &&
          data_collection_TAI_time < times[i+1]) 
      {
        start[0] = i;
        next_start[0] = i+1;
        break;
      }
    }
  }

  edge[0] = 1;

  /* Set values for other elements of "start" and "edge" arrays           */

  for (i = 1; i <= rank; i++) 
  {
    start[i] = 0;
    next_start[i] = 0;
    edge[i] = dims[i - 1];
  }

  /* Read the data from the LUT                                           */

  hdf_return = SDreaddata(sds_id, start, NULL, edge, data);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not read first set of data from"
                    " SDS \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDreaddata", 0,
               corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /*  Exit if there is only one record; no interpolation can be done.     */
 
  if (num_lut_copies == 1) return returnStatus;

  /* Get the times and their difference, which must be positive.          */
  /* Scale the TAI times.                                                 */
  data_coll_time = data_collection_TAI_time;
  lut_data_time = times[start[0]];
  next_lut_data_time = times[next_start[0]];
  slope_denom = (next_lut_data_time - lut_data_time);

  /* If the denominator of the slope is 0 then there is a problem.        */
  /* Likewise, a negative slope indicates the times are not in            */
  /* increasing order.                                                    */

  if (slope_denom <= 0) 
  {
    returnStatus = MODIS_F_NOK;
    sprintf(errmsgbuf, 
            "LUT associated with \"%s\" has times out of order.",
                    name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
         goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /* Linear Proportion (x - x1)/(x2 - x1)                                 */
  xproportion = (data_coll_time - lut_data_time)/slope_denom;
  
  /*  Get total number of LUT elements                                    */
  for (num_lut_elements = 1, i = 0; i < rank; i++)
    num_lut_elements *= dims[i];
 
  /* Read the data from the next higher LUT                               */
 
  data2 = malloc((unsigned) num_lut_elements*pw_size);
  if (!(data2)) 
  {
    returnStatus = MODIS_F_OUT_OF_MEMORY;
    sprintf(errmsgbuf,
       "Could not allocate memory for second LUT associated "
                    "with \"%s\".", name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, NULL, 0, NULL, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  hdf_return = SDreaddata(sds_id, next_start, NULL, edge, data2);
  if (hdf_return == FAIL) 
  {
    returnStatus = MODIS_F_HDF_ERROR;
    sprintf(errmsgbuf, 
            "Could not read second set of data from"
                    " SDS \"%s\"." , name);
    L1BErrorMsg(location, returnStatus, errmsgbuf, "SDreaddata", 0,
                corruptinputfile, False);
    goto TDLUT_ReadPiecewiseLinearFunction_exit;
  }

  /*  Do the interpolation.                                               */
 
  /*  Begin Interpolation Block                                           */
  {
    float32 first_val, second_val, interp_val; 
    float32 *first_data = (float32 *) data;
    float32 *next_data  = (float32 *) data2;

    for (i = 0; i < num_lut_elements; i++)  
    {
      first_val   = *first_data;
      second_val  = *next_data;
      interp_val  = first_val + xproportion*(second_val - first_val);
      *first_data = interp_val;
      first_data++;
      next_data++;
    }

  }
  /*  End Interpolation Block                                             */


TDLUT_ReadPiecewiseLinearFunction_exit:

  /*  Free up the memory used by data2                                    */
  if (data2 != NULL) 
    free(data2);

  if (sds_id != FAIL) SDendaccess(sds_id);

  return returnStatus;
}   

PGSt_SMF_status Expand_BDSM_LUT (
    void *data,      /* generic pointer to the unexpanded LUT array */
    void **data_new, /* generic pointer to the array to fill and return */
    int32 data_type, /* HDF data type of the data in the file */
    int32 lead_dim,  /* NUM_TIMES in time-dependent LUTs */
    int32 *n_bytes   /* size of filled array in bytes */
    )
/*
!C****************************************************************************
!Description:
   Given a LUT array containing values only for valid BDSM combinations,
   expand the array with fill values occupying invalid BDSM combinations.

!Input Parameters:
   void *data         Generic pointer to the unexpanded LUT array
   int32 data_type    HDF data type of the data in the file
   int32 lead_dim     NUM_TIMES in time-dependent LUTs

!Output Parameters:
   void **data_new    Generic pointer to the array to fill and return

!Revision History:

 Revision 01.02 August 9, 2002
 Removed illegal printf statement.
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 01.01 June 28, 2002    Razor Issue #161
 Handle expansion in memory, not files.  Moved to L1B_Tables.c.
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 01.00 December 31, 2001
 Initial development
 Gwyn Fireman, SAIC-GSO <fireman@mcst.gsfc.nasa.gov>

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  typedef enum {
    INT8, UINT8, INT16,
    UINT16, INT32, UINT32,
    FLOAT32, FLOAT64,
    NUM_DATA_TYPES
  } alltypes_t;
  int32 all_data_types[NUM_DATA_TYPES] = {
    DFNT_INT8, DFNT_UINT8, DFNT_INT16,
    DFNT_UINT16, DFNT_INT32, DFNT_UINT32,
    DFNT_FLOAT32, DFNT_FLOAT64
  };
  int32   data_var_size[NUM_DATA_TYPES] = {
    sizeof(int8), sizeof(uint8), sizeof(int16),
    sizeof(uint16), sizeof(int32), sizeof(uint32),
    sizeof(float32), sizeof(float64)
  };

  void * datap;
  float64 fill_64 = -999.0;
  int32 rank = 4;
  int32 dims[] = {0,	                /* override 2-d array dimensions */
		  NUM_REFLECTIVE_BANDS,
		  MAX_DETECTORS_PER_BAND,
		  MAX_SAMPLES_PER_BAND,
		  NUM_MIRROR_SIDES};
  int32 i, i_type, arrsize;
  int32 i_T, i_B, i_D, i_S, i_M;	/* loop indices for 
					   time, band, det, sf, ms */
  char *location = "Expand_BDSM_LUT";


/************************ MACRO FOR INSERTING LUT VALUES *******************/
#define INSERT_LUT_VALUES                                           \
  /* for every LUT time, */                                         \
  i = 0;                                                            \
  for (i_T = 0; i_T < dims[0]; i_T++)                               \
  {                                                                 \
    /*                                                              \
     *  For every 250m Band,                                        \
     *    For every 250-m Detector, SubFrame and Mirror Side,       \
     *      Write value to new array.                               \
     */                                                             \
    for (i_B = 0; i_B < NUM_250M_BANDS; i_B++)                      \
      for (i_D = 0; i_D < MAX_DETECTORS_PER_BAND; i_D++)            \
	for (i_S = 0; i_S < MAX_SAMPLES_PER_BAND; i_S++)            \
	  for (i_M = 0; i_M < NUM_MIRROR_SIDES; i_M++)              \
	    {data_temp[i] = *data_orig++;                           \
	    i++;}                                                   \
                                                                    \
    /*                                                              \
     *  For every 500m Band,                                        \
     *    For every 250-m Detector, SubFrame and Mirror Side,       \
     *      If Detector and SubFrame are meaningful for 500-m bands,\
     *        Write value to new array.                             \
     */                                                             \
    for (i_B = 0; i_B < NUM_500M_BANDS; i_B++)                      \
      for (i_D = 0; i_D < MAX_DETECTORS_PER_BAND; i_D++)            \
	for (i_S = 0; i_S < MAX_SAMPLES_PER_BAND; i_S++)            \
	  for (i_M = 0; i_M < NUM_MIRROR_SIDES; i_M++)              \
	    {if ((i_D < DETECTORS_PER_500M_BAND) &&                 \
		 (i_S < NUM_500M_SUBSAMP))                          \
	      data_temp[i] = *data_orig++;                          \
	    i++;}                                                   \
                                                                    \
    /*                                                              \
     *  For every 1km Reflective Solar Band,                        \
     *    For every 250-m Detector, SubFrame and Mirror Side,       \
     *      If Detector and SubFrame are meaningful for 1km bands,  \
     *        Write value to new array.                             \
     */                                                             \
    for (i_B = 0; i_B < NUM_1000M_REFL_BANDS; i_B++)                \
      for (i_D = 0; i_D < MAX_DETECTORS_PER_BAND; i_D++)            \
	for (i_S = 0; i_S < MAX_SAMPLES_PER_BAND; i_S++)            \
	  for (i_M = 0; i_M < NUM_MIRROR_SIDES; i_M++)              \
	    {if ((i_D < DETECTORS_PER_1KM_BAND) &&                  \
		 (i_S < NUM_1KM_SUBSAMP))                           \
	      data_temp[i] = *data_orig++;                          \
	    i++;}                                                   \
                                                                    \
  }                             /* i_T */

/************************ END MACRO ************************/
 

  /* Calculate total number of values in new LUT array */

  dims[0] = lead_dim;
  for (arrsize = 1, i = 0; i < rank+1; i++)
    arrsize *= dims[i];
  
  /* Determine index of the data type. */

  for (i_type = 0; i_type < NUM_DATA_TYPES; i_type++) {
    if (all_data_types[i_type] == data_type) break;
  }
  if (i_type == NUM_DATA_TYPES)
    {
      returnStatus = MODIS_F_INVALID_ARGUMENT;
      L1BErrorMsg(location, returnStatus, "Invalid data type", 
		  "Expand_BDSM_LUT", 0, NULL, False);
      return returnStatus;
    }
  *n_bytes = arrsize * data_var_size[i_type];

  /* Allocate memory for new array */

  datap = malloc((unsigned) (*n_bytes));
  if (!(datap))
    {
      returnStatus = MODIS_F_OUT_OF_MEMORY;
      L1BErrorMsg(location, returnStatus, "Could not allocate memory", 
		  "Expand_BDSM_LUT", 0, NULL, False);
      return returnStatus;
    }

  /* Cast array to proper type, initialize to fill value and expand */

  if (data_type == DFNT_INT8)
    {
      int8 *data_orig = (int8 *)data;
      int8 *data_temp = (int8 *)datap;
      int8 fill_value = (int8)fill_64;
      for (i = 0; i < arrsize; i++)
	data_temp[i] = fill_value;
      INSERT_LUT_VALUES;
      *data_new = (void *) data_temp;
    }
  else if (data_type == DFNT_INT16)
    {
      int16 *data_orig = (int16 *)data;
      int16 *data_temp = (int16 *)datap;
      int16 fill_value = (int16)fill_64;
      for (i = 0; i < arrsize; i++)
	data_temp[i] = fill_value;
      INSERT_LUT_VALUES;
      *data_new = (void *) data_temp;
    }
  else if (data_type == DFNT_INT32)
    {
      int32 *data_orig = (int32 *)data;
      int32 *data_temp = (int32 *)datap;
      int32 fill_value = (int32)fill_64;
      for (i = 0; i < arrsize; i++)
	data_temp[i] = fill_value;
      INSERT_LUT_VALUES;
      *data_new = (void *) data_temp;
    }
  else if (data_type == DFNT_FLOAT32)
    {
      float32 *data_orig = (float32 *)data;
      float32 *data_temp = (float32 *)datap;
      float32 fill_value = (float32)fill_64;
      for (i = 0; i < arrsize; i++)
	data_temp[i] = fill_value;
      INSERT_LUT_VALUES;
      *data_new = (void *) data_temp;
    }
  else if (data_type == DFNT_FLOAT64)
    {
      float64 *data_orig = (float64 *)data;
      float64 *data_temp = (float64 *)datap;
      float64 fill_value =  fill_64;
      for (i = 0; i < arrsize; i++)
	data_temp[i] = fill_value;
      INSERT_LUT_VALUES;
      *data_new = (void *) data_temp;
    }
  else returnStatus = MODIS_F_INVALID_ARGUMENT;
  
  return returnStatus;
  
}				/* end of Expand_BDSM_LUT */

int32 BDSM_index (char *ascii_file)
/*
!C****************************************************************************
!Description:
   Function returns index of LUT in refl_lut structure,
   or -1 if LUT is not a Reflective BDSM LUT.

!Input Parameters:
   char *ascii_file             name of the ASCII file holding the LUT

!Output Parameters:
 (none)

!Revision History:

 Revision 01.01 June 28, 2002    Razor Issue #161
 Moved to L1B_Tables.c.
 Gwyn Fireman, SAIC-GSO <Gwyn.Fireman@gsfc.nasa.gov>

 Revision 01.00 December 31, 2001
 Initial development
 Gwyn Fireman, SAIC-GSO <fireman@mcst.gsfc.nasa.gov>

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32 returnStatus = -1;	/* return status */
  int32 i_lut;			/* index of LUT in refl_lut structure */
  int32 nluts;			/* number of LUTs in refl_lut structure */
  char **dimnames;		/* name of each LUT dimension */

  /* 
   *  If LUT ascii name is null, then exit
   */
  if (!ascii_file)
    return returnStatus;

  /* 
   *  If LUT is not in refl_luts list, then exit
   */
  nluts = 0;				/* find # of LUTs in structure */
  while (refl_luts[nluts].name) nluts++;
  
  for (i_lut = 0; i_lut < nluts; i_lut++)
    {
      if (refl_luts[i_lut].ascii_file)	/* don't test against nulls */
	if (!strcmp (ascii_file, refl_luts[i_lut].ascii_file))
	  break;		/* found LUT name */
    }

  if (i_lut == nluts)
    return returnStatus;

  /* 
   * If LUT is not of rank 4, then exit
   */
  if (refl_luts[i_lut].rank != 4)
    return returnStatus;

  /* 
   *  If LUT does not have BDSM dimensions, then exit
   */
  dimnames = refl_luts[i_lut].dimnames;
  if (strcmp (dimnames[0], "NUM_REFLECTIVE_BANDS") ||
      strcmp (dimnames[1], "MAX_DETECTORS_PER_BAND") ||
      strcmp (dimnames[2], "MAX_SAMPLES_PER_BAND") ||
      strcmp (dimnames[3], "NUM_MIRROR_SIDES"))
    return returnStatus;

  /* 
   *  LUT passes all tests;
   *  return index of LUT in refl_lut structure
   */
  return returnStatus = i_lut;

}				/* end of BDSM_LUT_index */
