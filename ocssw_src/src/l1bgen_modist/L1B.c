#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    <math.h>
#include    <time.h>
#include    "PGS_TD.h"
#include    "Preprocess.h"      /* for Preprocess_L1A */
#include    "Metadata.h"        /* for Gran_Meta_Cal, Write_Gran_Metadata */
#include    "Emissive_Cal.h"    /* for Emissive_Cal */
#include    "Reflective_Cal.h"  /* for Reflective_Cal */
#include    "L1B_Tables.h"      /* for Read_Lookup_Tables */
#include    "L1B_Setup.h"       /* for L1B_Setup */
#include    "Granule.h"         /* for Read_L1A_EV_Scan, Aggregate_L1B,
                                 *     Write_L1B_EV_Scan, Close_L1A_Granule,
                                 *     Close_L1B_Granule */

int main(void)
/*
!C**********************************************************************
!Description:  L1B main routine.

!Input Parameters:

!Output Parameters:

!Revision History:
 $Log: L1B.c,v $
 Revision 1.15  2008-11-18 14:38:34-05  xgeng
 merge branch for V6.0.0

 Revision 1.13.2.2  2008/06/05 14:07:59  xgeng
 Removed the call of Fill_Dead_Detector_SI to remain the filled value for pixels affected by dead detector and noisy subframe

   Revision 1.13  2005/01/18 19:34:44  ltan
   MOD_PR02_TERRA update to V5.0.4
   
   Revision 02.19  October 16, 2004  Razor Issue #200
   Added "void" in the function main's parameter list to comply with ESDIS 
   guideline.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 02.18  March 26, 2003  
   Removed Check on Terra satellite and "_Terra" in Band 5 to Band 26
   crosstalk correction names because algorithm now applies to 
   Terra or Aqua.  Removed passing L1A_Gran to Band_26_Crosstalk_Correction
   because it is no longer necessary to check on which platform is being
   processed in that routine.
   Alice Isaacman  SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.17  March 19, 2002  Razor Issue #182
   Code added only applicable to MODIS/TERRA (PFM) Processing:
   Added call to Band_26_Crosstalk_Correction which corrects Band 26
   values using aggregated Band 5 values and LUT correction factors.
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.16  March 11, 2002  Razor Issue #174
   Added Passing of L1B_Gran to Preprocess.c to pass needed RVS correction
   coefficients.
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.15 November 12, 2001   Razor Issue #169
   Added boolean skip_night_hi_res to variable list and changed name 
     of boolean isdaymode to daymode_scan.
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.14  November 6, 2001
   Changed call to Read_Lookup_Tables to include runtime parameters so that
   a check can be made on MCST Version from the PCF file (Razor issue #167).
   Alice Isaacman (Alice.R.Isaacman.1@gsfc.nasa.gov), SAIC GSO

   Revision 02.13  October 29, 2000
   Razor issue 142, Add LUT to control scan checking QA for split scans
   and scan quality array.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.12  September 21, 2000
   Added counting num_rsb_at_night_scans as per Razor issue 137.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   ... (many changes not logged) ...

 Revision 02.11 Feb. 8, 1999
 Moved call to function Write_L1B_ScanMeta to L1B (main).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 
 Revision 02.10 June 1998
 Modified the Preprocess_L1A_Data() calls to include the refl_tables into
 parameters
 Zhenying Gu (zgu@ltpmail.gsfc.nasa.gov)
 
 Revision 02.10 Apr. 1998
 Removed the Cleanup_Preprocess_Data(&PP), because PP is changed to static
 allocation.  Removed the Cleanup_Emiss_Tables(), because Planck Integration
 is decided to carry on with in-line calculation. It may will be switched
 back to using SIL tables, then we may need to add it back.
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)
 
 Revision 02.10 Apr. 10, 1998
 Removed the momos declaration (moved to the L1B_granule_t in Granule.h)
 Removed the momos_t parmater from the L1B_Setup() and Reflective_Cal()
 calls - they're now included in the L1B_Gran parameter.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Apr. 1998
 Replaced L1B_Metadata() with gran_meta_cal() for the time being. It will
 be gran_meta_cal() and Write_Granule_Metadata() later.
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
 
 Revision 02.10 Apr. 1998
 Changed tables from dynamic to static allocation.
 Changed Cleanup_Lookup_Tables() to Cleanup_Emiss_Tables() 
 Added the momos_t parameter to the L1B_Setup() and
 Reflective_Cal() function calls.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Add argument MirrorSide to the Emissive_Cal argument list 
 Replace PP by PP->PP_Emiss in Emissive_Cal
 Replace PP by PP->PP_Refl in Reflective_Cal
 Shi-Yue Qiu (syqiu@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Mar. 1998
 Changed order of calibration to emissive followed by reflective
 (since reflective calibration needs band 28 radiances).
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.10 Jan. 1998
 Introduced lookup_tables_t, which includes refl_tables_t, emiss_tables_t, and 
 the new QA_tables_t.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 02.00 Jan. 1997
 Introduced refl_tables_t, emiss_tables_t, and Preprocess_Data_t,
 restructured L1A_granule_t, L1A_Scan_t, L1B_granule_t, and L1B_scan_t,
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.05 1995/12/07 14:15:00
 Changed Therm_ prefix on lookup tables to Emiss_ prefix
 Joan Baden(baden@highwire.gsfc.nasa.gov)

 Revision 01.04 1995/11/20 11:24:19
 Changed call from Setup() to L1B_Setup()
 John Hannon(hannon@highwire.gsfc.nasa.gov)

 Revision 01.03 1995/11/13 09:47:09
 Changed Derived_Housekeeping_t to Derived_Engineering_Data_t and DH to DED
 John Hannon(hannon@highwire.gsfc.nasa.gov)

 Revision 01.02 1995/10/20 11:17:25
 Implement direct PGS toolkit calls and populate comments
 Joan Baden(baden@highwire.gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)

 Revision 01.01 1995/05/17
 Initial header development
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

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
  PGSt_SMF_status      returnStatus = MODIS_S_OK;
  Run_Time_Parameters_t *runtime_params;
  lookup_tables_t      *tables;
  QA_Data_t            *QA;            /*Contains all L1B QA items */
  Preprocess_Data_t    *PP;            /*Contains OBC/Eng averages/statictcs*/
  L1A_granule_t        *L1A_Gran;
  L1A_Scan_t           *L1A_Scan;
  L1B_granule_t        *L1B_Gran;
  L1B_Scan_t           *L1B_Scan;
  L1B_Scan_Metadata_t  *L1B_Scan_Meta;
  L1B_Gran_Metadata_t  *L1B_Gran_Meta;

  char                 *location = "main";
  time_t tnow;
  
    /*Scan index (0,L1A_Gran.num_scans)*/
  int16                S            = 0;  

    /* Flag denoting whether a scan is "Day" */
  boolean              daymode_scan;      

    /* 
     * Flag denoting whether to skip writing a 250m- or 500m-resolution
     * granule containing only NIGHT mode scans 
     */
  boolean              skip_night_hi_res = FALSE;
  
  /*
   *   Flag to processing to signal whether or not to perform
   *   Band 26 correction from Band 5
   */
  boolean              do_B26_B5_corr = FALSE;
  
  // Allocate structures
  runtime_params = (Run_Time_Parameters_t*) calloc(1, sizeof(Run_Time_Parameters_t));
  tables = (lookup_tables_t*) calloc(1, sizeof(lookup_tables_t));
  QA = (QA_Data_t*) calloc(1, sizeof(QA_Data_t));
  PP = (Preprocess_Data_t*) calloc(1, sizeof(Preprocess_Data_t));
  L1A_Gran = (L1A_granule_t*) calloc(1, sizeof(L1A_granule_t));
  L1A_Scan = (L1A_Scan_t*) calloc(1, sizeof(L1A_Scan_t));
  L1B_Gran = (L1B_granule_t*) calloc(1, sizeof(L1B_granule_t));
  L1B_Scan = (L1B_Scan_t*) calloc(1, sizeof(L1B_Scan_t));
  L1B_Scan_Meta = (L1B_Scan_Metadata_t*) calloc(1, sizeof(L1B_Scan_Metadata_t));
  L1B_Gran_Meta = (L1B_Gran_Metadata_t*) calloc(1, sizeof(L1B_Gran_Metadata_t));

  setlinebuf(stdout);
  printf("MODIS Terra L1B version %s, built %s %s\n",
         PGE02_VERSION, __DATE__, __TIME__);

    /*
     * Read all run-time parameters from PCF.
     */

  returnStatus = Read_Run_Time_Parameters(runtime_params);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg(location, returnStatus, 
                "Could not read run-time parameters from PCF",
                "Read_Run_Time_Parameters", 0, NULL, True);
  }

    /* 
     * Open and read L1A files 
     */

  returnStatus = Open_and_Read_L1A(runtime_params,
                                   L1A_Gran, &skip_night_hi_res);
  if (returnStatus != MODIS_S_OK) {
    L1BErrorMsg(location, returnStatus, 
                "Could not open and read middle L1A granule",
                "Open_and_Read_L1A", FIRST_L1A_GRANULE, NULL, True);
  }
 
    /*
     * Read all lookup tables
     */

  returnStatus = Read_Lookup_Tables (L1A_Gran, tables, runtime_params);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, "Read_Lookup_Tables", 
                0, NULL, True);

  /* Determine whether Band 5 to Band 26 crosstalk correction wanted */  
  do_B26_B5_corr = tables->refl.B26_B5_Corr_Switch;

    /*
     * Determine other scans to be treated as missing.  This may adjust the
     * "missing_scan" array in L1A_Gran. (Must be called before L1B_Setup).
     */

  returnStatus = Determine_Other_Missing_Scans (tables, L1A_Gran);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, "Determine_Other_Missing_Scans",
                0, NULL, True);

    /* 
     * Read/Process/Write OBC/Eng data, compute PP
     */
  
  returnStatus = Preprocess_L1A_Data 
                         (tables,
                          L1A_Gran,
                          L1B_Gran,
                          QA,
                          PP);
                          
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"Preprocess_L1A_Data() in main(void), L1B.c");

    /*
     * Prepare for L1B EV Calibration: Open L1A EV SDSs, Open L1B EV files, 
     * Create L1B EV SDSs, & read/process ScanMetaData,
     * and calculate the MOMOs.
     */

  returnStatus = L1B_Setup (tables,
                            L1A_Gran,
                            L1B_Gran,
                            L1A_Scan,
                            L1B_Scan,
                            QA,
                            L1B_Scan_Meta,
                            skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"L1B_Setup() in main(void), L1B.c");
   
    /*
     * Read/Process/Write EV data in scan sequence
     */

  for (S = 0; S < L1A_Gran->num_scans; S++ )
  {
    if (S % 10 == 0) {
      time(&tnow);
      printf("scan: %d out of %d %s",
             S, L1A_Gran->num_scans, asctime(localtime(&tnow)));
    }

      /*
       * Skip missing scans.
       */

    if (L1A_Gran->missing_scan[S] == True) continue;
 
      /*
       * Read one scan of EV data
       */

    returnStatus = Read_L1A_EV_Scan(S, L1A_Gran, L1A_Scan);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR (returnStatus, "Read_L1A_EV_Scan() in main(void), L1B.c");

      /* 
       * Calibrate one scan of EV data with Emissive_Cal() & Reflective_Cal()
       */
    returnStatus = Emissive_Cal(S, 
                                L1A_Gran->MirrorSide[S],
                                L1B_Gran,
                                L1A_Scan,
                                L1B_Scan,
                                &(PP->PP_Emiss),
                                &(tables->emiss),
                                &(tables->QA),
                                QA,
                                L1B_Scan_Meta);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR (returnStatus, "Emissive_Cal() in main(void), L1B.c");
    
    
     /************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
      Reflective_Cal_Band_Flag = ALL_REFLECTIVE_BANDS;
#endif /* WRITE_BAND_26_SDS */
     /************************** End Band 26 Section ***************************/

      /*Only do Reflective_Cal() for Day Scans*/

    if (strcmp(L1A_Gran->ScanType[S], "Day") == SAME)
    {
      daymode_scan = True;
      returnStatus = Reflective_Cal(S, 
                                    L1A_Gran,
                                    L1B_Gran,
                                    L1A_Scan,
                                    L1B_Scan,
                                    PP,
                                    &(tables->refl),
                                    &(tables->QA.common_QA_tables),
                                    &(QA->QA_common));
      if (returnStatus != MODIS_S_OK)
        SMF_ERROR (returnStatus, "Reflective_Cal() in main(void), L1B.c");
     
        /*
         * Aggregation (needed only for Day-mode and reflective bands)
         */
    
      returnStatus = Aggregate_L1B (L1B_Scan);
      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus, "Aggregate_L1B() in main(void), L1B.c");

    /*
     * If in day mode and the Band 5 to Band 26 crosstalk
     * correction switch is set, perform the correction.
     */ 
	
      if (do_B26_B5_corr == TRUE)
	    { 
        returnStatus = Band_26_Crosstalk_Correction (
                         L1B_Scan,
                         tables->refl.B26_B5_Frame_Offset,
#ifdef USE_B5_RAD_OFFSET           
                         L1B_Gran->SO.rad_offset_RefSB[MODIS_BAND5_INDEX],
#endif /* USE_B5_RAD_OFFSET */            
                         L1B_Gran->b26_fr_b5_scaled_corr,
                         &(QA->QA_common),
                         L1B_Gran->valid_pixels,
                         L1B_Gran->negative_value_below_noise_pixels,
                         L1B_Gran->bad_data_flag,
                         daymode_scan,
                         do_B26_B5_corr);
        
      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus, 
          "Band_26_Crosstalk_Correction() in main(void), L1B.c");
	    }

    }

    else 
      
    {
      daymode_scan = False;
      QA->QA_common.num_rsb_at_night_scans++;
    }

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

    if (strcmp(L1A_Gran->ScanType[S], "Day") != SAME)
    {
      Reflective_Cal_Band_Flag = REFLECTIVE_BAND_26_ONLY;
      returnStatus = Reflective_Cal(S, 
                                    L1A_Gran,
                                    L1B_Gran,
                                    L1A_Scan,
                                    L1B_Scan,
                                    PP,
                                    &(tables->refl),
                                    &(tables->QA.common_QA_tables),
                                    &(QA->QA_common));
      if (returnStatus != MODIS_S_OK)
        SMF_ERROR(returnStatus, "Reflective_Cal(), in main(void), L1B.c");
    }
    returnStatus = Copy_Band_26_Data(L1B_Scan);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR(returnStatus, "Copy_Band_26_Data(), in main(void), L1B.c");

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

      /*
       * Fill the dead-detector SI values with an artificial value from
       * the adjacent live detectors.
       */
/*  remain filled value 
    returnStatus = Fill_Dead_Detector_SI
                             (daymode_scan,
                              tables->QA.common_QA_tables.dead_detector,
                              L1B_Scan,
                              L1B_Gran,
                              &(QA->QA_common));
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR (returnStatus,"Fill_Dead_Detector_SI() in main(void)");
*/

      /*
       * Write one scan of EV data
       */
        
    returnStatus = Write_L1B_EV_Scan (S, 
                                      L1B_Gran,
                                      L1B_Scan,
                                      daymode_scan);
    if (returnStatus != MODIS_S_OK)
      SMF_ERROR (returnStatus,"Write_L1B_EV_Scan() in main(void), L1B.c");
      
  }

    /*
     * Write scan metadata
     */

  returnStatus = Write_L1B_ScanMeta (L1B_Scan_Meta,
                                     L1A_Gran,
                                     QA,
                                     skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Write_L1B_ScanMeta() in main(void), L1B.c");
      
    /*
     * Write granule metadata
     */

  returnStatus = Gran_Meta_Cal(L1A_Gran,
                               L1B_Gran,
                               PP,
                               QA,
                               L1B_Scan_Meta,
                               L1B_Gran_Meta);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"Gran_Meta_Cal() in main(void), L1B.c");
   
  returnStatus = Write_Gran_Metadata(runtime_params,
                                     L1B_Gran_Meta,
                                     QA,
                                     PP,
                                     tables,
                                     L1A_Gran,
                                     skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"L1B_Metadata() in main(void), L1B.c");
    
    /* 
     * Done with this granule
     */
  returnStatus = Close_L1A_Granule (L1A_Gran,
                                    L1A_Scan);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"Close_L1A_Granule() in main(void), L1B.c");

  returnStatus = Close_L1B_Granule (L1B_Gran,
                                    L1B_Scan,
                                    skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus,"Close_L1B_Granule() in main(void), L1B.c");

  return(0);

}


