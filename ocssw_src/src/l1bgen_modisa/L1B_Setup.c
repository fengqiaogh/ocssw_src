/*****************************************************************************

File: L1B_Setup.c

External functions:
   L1B_Setup
   Write_L1B_ScanMeta
   Determine_Other_Missing_Scans

Other functions:
   Open_L1A_EV_SDS
   Calculate_Earth_Sun_Distance
   Calculate_RSB_Cal_Coeff
   Init_L1B_ScaleOffset
   Calculate_B26_B5_Correction
   Open_W_L1B_Granule
   Create_L1B_Swath
   Write_Swath_Band_Number
   Open_L1B_EV_SDS
   Set_L1B_EV_SDS_Attrs
   Set_Unit_Range_Fillvalue
   Get_SDS_id
   Set_SDS_Attributes
   Create_Band_Subsetting_SDS
   Copy_Geo_SDS
   Scan_Meta_Cal
   Calculate_DCR_Change
   Init_QA_Parameters
   Determine_Split_Scans
   Get_Split_Scan_Indexes
   Set_UI_ConvertToPercent_Attrs

Note:
   In the the various functions, there is a return statement after each
   L1BErrorMsg function call.  These are actually not necessary if the last
   argument to L1BErrorMsg is "True". These returns make it easier to test
   error out conditions (with SMF_ERROR used in UNIT_TEST_MODE_ONLY -- it does
   not actually exit).

Revision History:
   $Log: L1B_Setup.c,v $
   Revision 1.26  2015-03-17 10:48:11-04  xgeng
   Fix the anomaly before and after the sector rotation due to the timing mismatch between the command performed and the status updated.

   Revision 1.25  2010-11-15 11:37:49-05  xgeng
   Added a QA flag to identify the scan with both sides of PCLW electronics on.

   Revision 1.24  2008/11/18 16:25:15  xgeng
   merge the branch for V6.0.1

   Revision 1.23.1.3  2008/06/02 15:38:38  xgeng
   Initialize dead_subframe_pixels and num_dead_subframe_EV_data


   Revision 1.23  2006/10/27 14:34:39  ltan
   Removed the "Mixed" option from ScanType since L1A "Scan Type" never "Mixed".  Some comments regarding Night mode handling corrected.


 *****************************************************************************/

#include    <string.h>
#include    <math.h>
/* Replace <malloc.h> with the standard head file <stdlib.h> for ANSI-C 
   compliance (Razor Issue 173). */
#include    <stdlib.h>
#include    <time.h>
#include    "L1B_Setup.h"
#include    "L1B_Tables.h"
#include    "Metadata.h"
#include    "HDF_Lib.h"
#include    "PGS_PC.h"
#include    "PGS_TD.h"
#include    "PGS_Error_Codes.h"
#include    "FNames.h"
#include    "L1B_SetupP.h"  /* local prototypes and structures */

/*
 * L1A_EV_SDS_NAME -- names of the earth-view SDSs in the L1A granule.
 *                    Used in opening access to those SDSs.
 */
char *L1A_EV_SDS_NAME[NUM_L1A_RESOLUTIONS] = {
  "EV_250m",   "EV_500m",   "EV_1km_day",   "EV_1km_night"
};

/*-----------------------------------------
  L1B file spec
  (SI = Scaled Integer)
  (UI = Uncertainty Index)
  (SU = Samples Used)
  L1B_EV_SDS_index_t    Defines macros for the indices in the SDS name arrays.
  L1B_EV_SDS_NAME       Names of all earth-view SDSs in the L1B files.
  L1B_EV_SDS_LONG_NAME  Long names for the EV SDSs
-----------------------------------------*/
typedef enum 
{
  SI_250m,         UI_250m,
  SI_500m,         UI_500m,
  SI_1km,          UI_1km,
  SI_1km_EMISS,    UI_1km_EMISS,
  SI_250to500m,    UI_250to500m,   SU_250to500m,
  SI_250to1km,     UI_250to1km,    SU_250to1km,
  SI_500to1km,     UI_500to1km,    SU_500to1km,
  NUM_L1B_EV_SDS
} L1B_EV_SDS_index_t;

char *L1B_EV_SDS_NAME[NUM_L1B_EV_SDS] = 
{
  "EV_250_RefSB",           "EV_250_RefSB_Uncert_Indexes",
  "EV_500_RefSB",           "EV_500_RefSB_Uncert_Indexes",
  "EV_1KM_RefSB",           "EV_1KM_RefSB_Uncert_Indexes",
  "EV_1KM_Emissive",        "EV_1KM_Emissive_Uncert_Indexes",
  "EV_250_Aggr500_RefSB",   "EV_250_Aggr500_RefSB_Uncert_Indexes",   
                                "EV_250_Aggr500_RefSB_Samples_Used",
  "EV_250_Aggr1km_RefSB",   "EV_250_Aggr1km_RefSB_Uncert_Indexes",   
                                "EV_250_Aggr1km_RefSB_Samples_Used",
  "EV_500_Aggr1km_RefSB",   "EV_500_Aggr1km_RefSB_Uncert_Indexes",   
                                "EV_500_Aggr1km_RefSB_Samples_Used"
};

char *L1B_EV_SDS_LONG_NAME[NUM_L1B_EV_SDS] = 
{
  "Earth View 250M Reflective Solar Bands Scaled Integers",
  "Earth View 250M Reflective Solar Bands Uncertainty Indexes",
  "Earth View 500M Reflective Solar Bands Scaled Integers",
  "Earth View 500M Reflective Solar Bands Uncertainty Indexes",
  "Earth View 1KM Reflective Solar Bands Scaled Integers",
  "Earth View 1KM Reflective Solar Bands Uncertainty Indexes",
  "Earth View 1KM Emissive Bands Scaled Integers",
  "Earth View 1KM Emissive Bands Uncertainty Indexes",
  "Earth View 250M Aggregated 500M Reflective Solar Bands "
                                   "Scaled Integers",
  "Earth View 250M Aggregated 500M Reflective Solar Bands "
                                   "Uncertainty Indexes",
  "Earth View 250M Aggregated 500M Reflective Solar Bands "
                                   "Number of Samples Used in Aggregation",
  "Earth View 250M Aggregated 1km Reflective Solar Bands "
                                   "Scaled Integers",
  "Earth View 250M Aggregated 1km Reflective Solar Bands "
                                   "Uncertainty Indexes",
  "Earth View 250M Aggregated 1km Reflective Solar Bands "
                                   "Number of Samples Used in Aggregation",
  "Earth View 500M Aggregated 1km Reflective Solar Bands "
                                   "Scaled Integers",
  "Earth View 500M Aggregated 1km Reflective Solar Bands "
                                   "Uncertainty Indexes",
  "Earth View 500M Aggregated 1km Reflective Solar Bands "
                                   "Number of Samples Used in Aggregation"
};

/*
 * Band_subsetting_names   Names of the swath fields for bands subsetting SDSs.
 * L1B_EV_DIM_NAME         Dimension names.  1st column for band subsetting.
 */

char *Band_subsetting_names[NUM_L1A_RESOLUTIONS] = {
       "Band_250M",
       "Band_500M",
       "Band_1KM_RefSB",
       "Band_1KM_Emissive"
};

#define  L1B_EV_SDS_RANK            3

char *L1B_EV_DIM_NAME[NUM_L1A_RESOLUTIONS][L1B_EV_SDS_RANK] = {
  {"Band_250M",           "40*nscans",     "4*Max_EV_frames"},
  {"Band_500M",           "20*nscans",     "2*Max_EV_frames"},
  {"Band_1KM_RefSB",      "10*nscans",     "Max_EV_frames"},
  {"Band_1KM_Emissive",   "10*nscans",     "Max_EV_frames"}
};

/* Mapping offsets in scan and track directions:
*/
int32 L1B_EV_DIM_OFFSET[NUM_L1A_RESOLUTIONS][L1B_EV_SDS_RANK-1] = {
  {1, 0}, {0, 0}, {0, 0}, {0, 0}
};

/* Mapping fractional offsets in scan and track directions:
*/
float32 L1B_EV_DIM_FRAC_OFFSET[NUM_L1A_RESOLUTIONS][L1B_EV_SDS_RANK-1] = {
  {0.5, 0.0}, {0.5, 0.0}, {0.0, 0.0}, {0.0, 0.0}
};

  /*Neal's notation of geo_sds*/

  /* NOTE: only the name, src_name and type are currently being used
   *       in the code.
   */
/* Enclosed the rows in the initializer of array GEO_SDS with braces for
   ANSI-C compliance (Razor Issue 173). */
geolocation_sds_t GEO_SDS[NUM_GEO_SDS] = {
/* name            src_name         type          units       
     valid_range[2]  FillValue line_numbers frame_numbers scale_factor */
  {"Latitude"     , "Latitude"     , DFNT_FLOAT32, "degrees"   , 
   {-90.,90.}     , -999.          , "3,8"       , "3,8,13,...", 1.f}  ,
  {"Longitude"    , "Longitude"    , DFNT_FLOAT32, "degrees"   ,  
   {-180., 180.}  , -999.          , "3,8"       , "3,8,13,...", 1.f}  ,
  {"Height"       , "Height"       , DFNT_INT16  , "meters"    ,  
   {-400,10000}   , -32767         , "3,8"       , "3,8,13,...", 1.f}  ,
  {"SensorZenith" , "SensorZenith" , DFNT_INT16  , "degrees"   , 
   {0,18000}      , -32767         , "3,8"       , "3,8,13,...", 0.01f},
  {"SensorAzimuth", "SensorAzimuth", DFNT_INT16  , "degrees"   , 
   {-18000,18000} , -32767         , "3,8"       , "3,8,13,...", 0.01f},
  {"Range"        , "Range"        , DFNT_UINT16 , "meters"    , 
   {27000.,65535.}, 0              , "3,8"       , "3,8,13,...", 25.f} ,
  {"SolarZenith"  , "SolarZenith"  , DFNT_INT16  , "degrees"   , 
   {0,18000}      , -32767         , "3,8"       , "3,8,13,...", 0.01f},
  {"SolarAzimuth" , "SolarAzimuth" , DFNT_INT16  , "degrees"   ,   
   {-18000,18000} , -32767         , "3,8"       , "3,8,13,...", 0.01f},
  {"gflags"       , "gflags"       , DFNT_UINT8  , ""          , 
   {0.,0.}        , 255            , ""          , ""          , 1.}
};

extern int16 RFLAG;
extern int16 RSCL_FLAG;

/******************************** FUNCTIONS *********************************/

PGSt_SMF_status L1B_Setup (lookup_tables_t      *tables,
                           L1A_granule_t        *L1A_Gran,
                           L1B_granule_t        *L1B_Gran,
                           L1A_Scan_t           *L1A_Scan,
                           L1B_Scan_t           *L1B_Scan,
                           QA_Data_t            *QA,
                           L1B_Scan_Metadata_t  *L1B_Scan_Meta,
                           boolean              skip_night_hi_res)
/*
!C************************************************************************
!Description:
    This routine performs a variety of functions in preparation for EV
    calibration: (1) opens SDS access to each of the four L1A EV SDSs,
    (2) calculates radiance and reflectance coefficients, populating the
    "momos_t" structure, (3) sets the radiance, reflectance and emissive
    scales and offsets for the L1B data products, (4) creates HDFEOS Swaths
    and data fields for each L1B output file,  (5) creates band-subsetting
    SDSs, (6) opens SDS access for each of the EV SDSs in the L1B EV files
    and create attributes for all, (7) reads subsampled SDSs from geolocation
    file and write those data and their attributes into the 1km L1B file,
    (8) assigns nadir-frame lattitude and longitude to members of L1B_Scan_Meta
    to be written later in Write_L1B_ScanMeta, (9) assigns members of the
    L1B_Scan_Meta to be written to the L1B EV files later in
    Write_L1B_ScanMeta, and (10) initializes the QA values of total number of
    pixels, number of valid pixels, number of saturated pixels, number of
    missing pixels and pixels representing negative values below noise.

!Input Parameters:
    lookup_tables_t *tables
    L1A_granule_t   *L1A_Gran
    boolean         skip_night_hi_res   True if and only if all scans are 
                                     night mode scans and writing of 250m
                                     and 500m night mode granules has been
                                     turned off

!Output Parameters:
    L1B_granule_t   *L1B_Gran
    L1A_Scan_t      *L1A_Scan
    L1B_Scan_t      *L1B_Scan
    QA_Data_t       *QA

!Revision History:
 (continue at top of the file)

 Revision 02.15 March 26, 2003 Razor Issue #191
 Added assignment of SWIR OOB correction band from LUT value
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.14  March 26, 2003 Razor Issue #190
 Added Call to Calculate_B26_B5_Correction
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.13 November 15, 2001  Razor Issue #169
 Changed call to Open_W_L1B_Granule to include boolean skip_night_hi_res
 Added boolean skip_night_hi_res to input parameters to L1B_Setup
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.12  Feb 10, 1999
 Moved call to function Write_L1B_ScanMeta to L1B (main).
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 Feb. 1999
 The following L1A_Gran members are now read in Read_L1A_OBCEng:
 Number of Day mode scans, EV start time and Scan Type.  This eliminates
 the need for Get_EV_Start_Time and Read_L1A_ScanMeta.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 April 1998
 Added scan_meta_cal(). 
 Zhenying Gu(zgu@gscmail.gsfc.nasa.gov)
  
 Revision 02.10 April 1998
 Added the Get_EV_Start_Time() and Calculate_MOMOs() calls.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.00 March 1997
 Exluded everything associated with LookupTables & PreprocessData;
 changed processing metadata one scan a time to all scans at once.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Updated to match Version 1 Design Document
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir Kvaran(geir@highwire.gsfc.nasa.gov)

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
  PGSt_SMF_status  returnStatus = MODIS_S_OK;
  char *location = "L1B_Setup";   /* L1B location */
  SWIR_correction_tables_t *swir_tables = 
      &tables->refl.SWIR_correction_tables;

  /*
   * Assign L1B_Scan->band_X based on reflective LUT value
   * Assignment is as night band index for use with L1A data.
   */
   /* MODIS_BAND20_INDEX = 21    NUM_BANDS = 38    MODIS_BAND26_INDEX = 27  */


  if ((swir_tables->SWIR_corr_sending_band >= MODIS_BAND20_INDEX - 1 && 
       swir_tables->SWIR_corr_sending_band < MODIS_BAND26_INDEX - 1) ||
      (swir_tables->SWIR_corr_sending_band > MODIS_BAND26_INDEX - 1 && 
       swir_tables->SWIR_corr_sending_band < NUM_BANDS - 1))
    L1B_Scan->band_X = swir_tables->SWIR_corr_sending_band - MODIS_BAND20_INDEX + 1;    
  else
  {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "Band to use for SWIR OOB Correction is out of range.",
                NULL, REFLECTIVE_TABLES_FILE, 
                "This is most likely due to an invalid LUT file.", True);
    return returnStatus;
  }

  /*
   * Assign L1B_Gran->num_scans & L1B_Gran->num_day_scans
   */
  L1B_Gran->num_scans = L1A_Gran->num_scans;
  L1B_Gran->num_day_scans = L1A_Gran->num_day_scans;


  returnStatus = Open_L1A_EV_SDS (L1A_Gran, L1A_Scan);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,"Open_L1A_EV_SDS() in L1B_Setup");

  returnStatus = Calculate_Earth_Sun_Distance(L1A_Gran, 
                                              &L1B_Gran->Earth_Sun_Dist);
  if (returnStatus != MODIS_S_OK)
    L1BErrorMsg(location, returnStatus, NULL, 
                "Calculate_Earth_Sun_Distance",
                0, NULL, True);

  returnStatus = Calculate_RSB_Cal_Coeff(tables, 
                                         L1B_Gran->Earth_Sun_Dist,
                                         &L1B_Gran->RSB_Cal_Coeff);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,"Calculate_RSB_Cal_Coeff() in L1B_Setup");
  
  returnStatus = Init_L1B_ScaleOffset(&L1B_Gran->SO, 
                                      &L1B_Gran->RSB_Cal_Coeff,
                                      L1B_Gran->Earth_Sun_Dist, tables);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus,"Init_L1B_ScaleOffset() in L1B_Setup");

  /*
   * Scale the Band 5 to Band 26 crosstalk correction coefficients by the
   *   ratio of the Band 5 to Band 26 scales and store the result in the
   *   L1B_Gran->RSB_Cal_Coeff structure.
   */ 
  
  returnStatus = Calculate_B26_B5_Correction
                  (tables->refl.B26_B5_Corr,
                   L1B_Gran->b26_fr_b5_scaled_corr,
                   &L1B_Gran->SO);

  if (returnStatus != MODIS_S_OK)
    SMF_ERROR (returnStatus, 
             "Calculate_B26_B5_Correction() in L1B_Setup");
  
  returnStatus = Open_W_L1B_Granule (tables, L1B_Gran, L1B_Scan, 
      skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Open_W_L1B_Granule() in L1B_Setup");

  returnStatus = Copy_Geo_SDS (L1B_Gran, skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Copy_Geo_SDS() in L1B_Setup");
  
  returnStatus = Scan_Meta_Cal(tables,
                               L1A_Gran, L1B_Gran, L1B_Scan_Meta, QA);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Scan_Meta_Cal() in L1B_Setup");
  
  returnStatus = Calculate_DCR_Change(L1A_Gran, QA, L1B_Scan_Meta);
  if(returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Calculate_DCR_Change() in Gran_Meta_Cal()");

  returnStatus = Init_QA_Parameters (L1A_Gran, L1B_Gran, QA);
  if (returnStatus != MODIS_S_OK)
    SMF_ERROR(returnStatus, "Init_QA_Parameters() in L1B_Setup");

  return(MODIS_S_OK);
}

PGSt_SMF_status Open_L1A_EV_SDS (L1A_granule_t  *L1A_Gran,
                                 L1A_Scan_t     *L1A_Scan)
/*
!C****************************************************************************
!Description: Create SDS access to each of the four L1A EV SDSs. The sds_ids
              are stored in L1A_Scan. These SDS accesses are not terminated
              until after all processing has been complete.

!Input Parameters:
      L1A_Gran (->sd_id)      SD file access ID -- already opened.

!Output Parameters:
      L1A_Scan (->sds_id[R])  SDS access IDs for each EV SDS in L1A file.

!Revision History:
 (continue at top of the file)

 Revision 02.10 March 1997
 Modified to match the change in L1A_Scan_t. 
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 02.00 1996/07/02
 Version 2, initial redesign
 Neal Devine(neal.devine@gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Update to match Version 1 Design Document
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Neal Devine(neal.devine@gsfc.nasa.gov)
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

!END**************************************************************************
*/
{
  PGSt_SMF_status  returnStatus  = MODIS_S_OK;
  char *location = "Open_L1A_EV_SDS";
  int16   R          = 0;
  int32   sds_index  = 0;
  char errmsg[256];

  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
  {
    if ((RFLAG & (1 << R)) != 0) continue;

    sds_index = SDnametoindex (L1A_Gran->sd_id, L1A_EV_SDS_NAME[R]);
    if (sds_index == FAIL)
    {
      sprintf(errmsg, "Could not get SDS index of \"%s\".", 
              L1A_EV_SDS_NAME[R]);
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDnametoindex", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    L1A_Scan->sds_id[R] = SDselect (L1A_Gran->sd_id, sds_index);
    if (L1A_Scan->sds_id[R] == FAIL)
    {
      sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
              L1A_EV_SDS_NAME[R]);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDselect", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Calculate_Earth_Sun_Distance( L1A_granule_t  *L1A_Gran,
                                              float32        *Earth_Sun_Dist)
/*
!C****************************************************************************
!Description: This function calculates the Earth Sun distance based on the TAI
              time.

!Input Parameters:
  L1A_granule_t        *L1A_Gran       contains number of scans in middle 
                                       L1A granule and EV start time.

!Output Parameters:
  float32 *            Earth_Sun_Dist  to store Earth Sun distance

!Revision History:
 $Log: L1B_Setup.c,v $
 Revision 1.26  2015-03-17 10:48:11-04  xgeng
 Fix the anomaly before and after the sector rotation due to the timing mismatch between the command performed and the status updated.

 Revision 1.25  2010-11-15 11:37:49-05  xgeng
 Added a QA flag to identify the scan with both sides of PCLW electronics on.

 Revision 1.24  2008/11/18 16:25:15  xgeng
 merge the branch for V6.0.1

 Revision 1.23.1.3  2008/06/02 15:38:38  xgeng
 Initialize dead_subframe_pixels and num_dead_subframe_EV_data

 (continue at top of the file)

 Revision 1.20  2005/03/01 20:34:31  ltan
 Code changes to correct the HDFEOS dimension mapping offset for QKM data and
 1KM geolocation from 0 to 1.

 Revision 1.19  2005/01/18 21:57:36  ltan
 Added new file attributes prefixed with HDFEOS_FractionalOffset


  Revision 1.02    June 3, 2003  Razor Issue #193
  The calculation of the Earth-Sun distance was amended after it was
  discovered that the memo on which it is based has a typographical
  error; the "+" in front of the second order term should be a "-". 
  Comments were also clarified.
  Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

  Revision 01.01 November 7, 2001
  Change use of PI to PGS_PI (Razor issue #168); clean up some spelling
    and grammar
  Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

  Revision 01.00 Nov. 04, 1999
  Initial development.
  Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Calculate_Earth_Sun_Distance";
  float64 TAItime;
  float64 EMA;
  int32   middle_scan;

  if (L1A_Gran == NULL || Earth_Sun_Dist == NULL)
    L1BErrorMsg(location, MODIS_F_INVALID_ARGUMENT, 
        "Input parameter is NULL pointer.", NULL, 0, NULL, True);

  /*
   * The formula calculating earth sun distance uses one TAI time for 
   * the granule.  The memo referred to below does not specify with time 
   * to use. Logically, it would be near the middle of the granule.
   */


  middle_scan = L1A_Gran->num_scans/2;
  TAItime = L1A_Gran->EVStartTime_TAIsecond[middle_scan];
  
/*
 * Calculate EMA using the expression for m appropriate for use with secTAI93
 * (TAItime) input.   This equation appears on the last page of the memo, "The
 * Earth-Sun Distance  Correction for the MODIS Reflective Product", M0612. 
 * The additional modulo  360 equation in this reference can be ignored as it
 * will not affect the answer.
 */

  EMA = -2.15802 + 0.0000114074107 * TAItime - 
          1.56645046e-23 * TAItime * TAItime;

/*
 * Calculate Earth-Sun distance using equation (2) with e = 0.167075 from the 
 * memo "The Earth-Sun Distance Correction for the MODIS Reflective Product", 
 * M0612.
 */

  *Earth_Sun_Dist = 1.0001396 - 0.0167068 * cos((double)(EMA * PGS_PI/180.0)) 
                    - 0.0001396 * cos((double)(2 * EMA * PGS_PI/180.0));

  return returnStatus;
}

PGSt_SMF_status Calculate_RSB_Cal_Coeff(lookup_tables_t *tables,
                                        float32         E_S_Dist,
                                        RSB_Cal_Coeff_t *RSB_Cal_Coeff)
/*
!C****************************************************************************
!Description:
 This function calculates the reflective calibration coefficient m1_des_sq, 
 which is the values of LUT m1 multiplied by the square of the Earth-Sun
 distance for the current granule, and it's maximum values. It also calculates 
 minimum values of R*, from the time-dependent LUT quantities.  The calculated 
 values are stored in the RSB_Cal_Coeff_t structure for use in Reflective_Cal 
 and Init_L1B_ScaleOffset.

!Input Parameters:
 lookup_tables_t *  tables    Contains structure members for m1 and R*
                              time-dependent LUTs and dead detector list.
 float32            E_S_Dist  earth-sun distance
 
!Output Parameters:
 RSB_Cal_Coeff_t *  RSB_Cal_Coeff   Contains the structure members to be filled
                                    in this array: m1_des_sq, m1_des_sq_max, 
                                    R_star, R_star_min.

!Revision History:
 (continue at top of the file)

 Revision 01.01 Nov, 1999
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 01.00 May, 1999
 See SDF, Initial development.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/

{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16 band, det, sample, mirr_side;
  int16 det_490;  /* index for all detectors */
  int16 num_detectors[NUM_REFLECTIVE_BANDS] = 
          {40, 40, 20, 20, 20, 20, 20, 10,
           10, 10, 10, 10, 10, 10, 10, 10,
           10, 10, 10, 10, 10, 10};
  
  int16 num_subsamps [NUM_REFLECTIVE_BANDS] = 
          {4, 4, 2, 2, 2, 2, 2, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  int32 i, n;
  float32 *work;
  work = malloc(MAX_DETECTORS_PER_BAND * 
                MAX_SAMPLES_PER_BAND *
                NUM_MIRROR_SIDES *
                sizeof(float32));
  float32 m1;
  
  /* Calculate m1 and m1 max */
  
  det_490 = 0; 
  for (band = 0; band < NUM_REFLECTIVE_BANDS; band++ ) 
  {
    n = 0;

  /*
   * (NUM_REFLECTIVE_DETECTORS - 10) is the index of the first
   * detector in band 20. MODIS_BAND26_INDEX_AT_RES equals the
   * difference of band index between band 26 and band 20, which
   * is 6. Every band has 10 detectors for 1km bands.
   */

    if (band == NUM_REFLECTIVE_BANDS - 1) /* band 26 */
      det_490 = (NUM_REFLECTIVE_DETECTORS - 10) + 
                   MODIS_BAND26_INDEX_AT_RES * 10;

    for ( det = 0; det < num_detectors[band]; det++, det_490++ ) 
    {
      for ( sample = 0; sample < num_subsamps[band]; sample++ ) 
      {
        for ( mirr_side = 0; mirr_side < NUM_MIRROR_SIDES; mirr_side++ ) 
        {
          m1 = tables->refl.m1[band][det][sample][mirr_side];
          RSB_Cal_Coeff->m1_des_sq[band][det][sample][mirr_side] = 
                                  m1 * E_S_Dist * E_S_Dist;
         
          if (tables->QA.common_QA_tables.dead_detector[det_490] != 1)
          {
            work[n] = RSB_Cal_Coeff->m1_des_sq[band][det][sample][mirr_side];
            n++;
          }
        }  /* end loop through mirr_side */
      }    /* end loop through sample */
    }      /* end loop through det */


    RSB_Cal_Coeff->m1_des_sq_max[band] = work[0];
    for (i = 0; i < n; i++)
      if (RSB_Cal_Coeff->m1_des_sq_max[band] < work[i])
        RSB_Cal_Coeff->m1_des_sq_max[band] = work[i];
    
    if (RSB_Cal_Coeff->m1_des_sq_max[band] == 0 || n == 0)
      SMF_ERROR(MODIS_F_NOK, 
          "invalid m1 maximum value in Calculate_RSB_Cal_Coeff()");
  }/* end loop through band */ 
  
  free(work);
  
  return (returnStatus);
} 

PGSt_SMF_status Init_L1B_ScaleOffset (L1B_ScaleOffset_t *SO, 
                                      RSB_Cal_Coeff_t   *RSB_Cal_Coeff,
                                      float             E_S_Dist,
                                      lookup_tables_t   *tables)
/*
!C**********************************************************************
!Description: This routine sets the radiance, reflectance and emissive scales
              and offsets for the L1B data products.

!Input Parameters:
      lookup_tables_t      *tables        contains L_Min, L_Max, dn_star_Min
                                          dn_star_Max
      RSB_Cal_Coeff_t      *RSB_Cal_Coeff   contains m1_des_sq_max and 
                                          R_star_min
      float                E_S_Dist       Earth-sun distance in AU

!Output Parameters:
      L1B_ScaleOffset_t  *SO              contains all scales and offsets 
                                          for reflectance, radiance and count

!Revision History:
 (continue at top of the file)

 Revision 02.13 August 12, 1999
 Used L_Max, L_Min Luts
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.12 May 18, 1999
 Changed to match RSB calibration coefficient algorithm change.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov) and Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 
 Revision 02.11 Feb 16, 1999
 Removed outdated and unused variable DN_Star_Max_Band1_7.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.10 Apr 6, 1998
 Commented out the Reflective Scale and Offset expressions,
 awaiting the V2.1 algorithm.
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.00 Dec. 1996.
 Modified Scale and Offset expressions according to V2 spec.
 Zhidong Hao (hao@acrobat.gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Updated to match Version 1 Design Document.
 Neal Devine(neal.devine@gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)

 Revision 01.00 1995/12/05
 Initial development
 Shi_Yue Qiu(syqiu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   Reflectance does not contain the Sun_Adj factor which varies from slightly
   less than 1 to 1.03something.  See comment below.

   Replaced `SIL_t *E_Sun' by `refl_tables_t *refl_tables' and
   removed Sun_Adj as an input parameter, as it is not needed here.
     --- zh 11/21/96

!END********************************************************************
*/
{
  int16   B        = 0;
  int16   B_emiss  = 0;
  float32 E_sun_over_pi_B[NUM_REFLECTIVE_BANDS];
  int16 D, D_330;

  /*set dn_star_Max, dn_star_Min*/
  for (B = 0; B < NUM_REFLECTIVE_BANDS; B++)
  {
    SO->dn_star_Max[B] = tables->refl.dn_star_Max[B];
    SO->dn_star_Min[B] = tables->refl.dn_star_Min[B];
    if (SO->dn_star_Max[B] <= SO->dn_star_Min[B])
      SMF_ERROR(MODIS_F_NOK, 
          "LUT dn_star_Max <= dn_star_Min, Init_L1B_ScaleOffset");
  }
     
  /*set Emiss radiance scale/offset*/  
  B_emiss = 0;
  for (B = 0; B < L1A_BANDS_AT_RES[INDEX_1000M_EMISS]; B++)
  {
    if (B == MODIS_BAND26_INDEX_AT_RES)
      continue;
    SO->rad_scale_Emiss[B_emiss]    = 
                        (tables->emiss.L_Max[B_emiss] - 
                        tables->emiss.L_Min[B_emiss]) / DN15_SAT;
    SO->rad_offset_Emiss[B_emiss]   = 
                        -DN15_SAT * tables->emiss.L_Min[B_emiss]/
                        (tables->emiss.L_Max[B_emiss] - 
                        tables->emiss.L_Min[B_emiss]);
    B_emiss++;
  }

  /*
   * Calculate the mean E_sun_over_pi for each band.
   */

  D_330 = 0;
  for (B = 0; B < NUM_250M_BANDS; B++) 
  {
    E_sun_over_pi_B[B] = 0;
    for (D = 0; D < DETECTORS_PER_250M_BAND; D++, D_330++)
      E_sun_over_pi_B[B] += tables->refl.E_sun_over_pi[D_330];
    E_sun_over_pi_B[B]   /= (float32) DETECTORS_PER_250M_BAND;
  }
  for (     ; B < (NUM_250M_BANDS+NUM_500M_BANDS); B++) 
  {
    E_sun_over_pi_B[B] = 0;
    for (D = 0; D < DETECTORS_PER_500M_BAND; D++, D_330++)
      E_sun_over_pi_B[B] += tables->refl.E_sun_over_pi[D_330];
    E_sun_over_pi_B[B]   /= (float32) DETECTORS_PER_500M_BAND;
  }
  for (     ; B < NUM_REFLECTIVE_BANDS; B++) 
  {
    E_sun_over_pi_B[B] = 0;
    for (D = 0; D < DETECTORS_PER_1KM_BAND; D++, D_330++)
      E_sun_over_pi_B[B] += tables->refl.E_sun_over_pi[D_330];
    E_sun_over_pi_B[B]   /= (float32) DETECTORS_PER_1KM_BAND;
  }

  /*set RefSB radiance, reflectance and counts scale/offset*/
 
  for (B = 0; B < NUM_REFLECTIVE_BANDS; B++)
  {
    SO->counts_scale_RefSB[B]  = 
        (SO->dn_star_Max[B] - SO->dn_star_Min[B])/DN15_SAT;
    SO->counts_offset_RefSB[B] =  
        -DN15_SAT * SO->dn_star_Min[B] /
        (SO->dn_star_Max[B] - SO->dn_star_Min[B]);

    SO->refl_scale_RefSB[B] = 
        RSB_Cal_Coeff->m1_des_sq_max[B] * SO->counts_scale_RefSB[B];
    SO->refl_offset_RefSB[B] = 
        SO->counts_offset_RefSB[B];

    SO->rad_scale_RefSB[B] = 
        (E_sun_over_pi_B[B] / (E_S_Dist * E_S_Dist))
        * SO->refl_scale_RefSB[B];
    SO->rad_offset_RefSB[B] = 
        SO->counts_offset_RefSB[B];
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Open_W_L1B_Granule 
                  (lookup_tables_t       *tables,
                   L1B_granule_t         *L1B_Gran,
                   L1B_Scan_t            *L1B_Scan,
                   boolean               skip_night_hi_res)
/*
!C**********************************************************************
!Description: Create HDFEOS Swaths and data fields for each L1B output file.
              Also get SD & V HDF interface ids.  The band-subsetting data
              fields (one-dimensional arrays) are implemented as Vdata while
              the other fields are implemented as SDSs.  Separate band-
              subsetting SDSs are also created.  Open SDS access for each of
              the EV SDSs in the L1B EV files and create attributes for all.
              SDS access will remain open until after all scan-by-scan
              processing is complete.

!Input Parameters:
   lookup_tables_t  *tables
   L1B_granule_t    *L1B_Gran
   L1B_Scan_t       *L1B_Scan
   boolean          skip_night_hi_res   True if and only if all scans are 
                                     night mode scans and writing of 250m
                                     and 500m night mode granules has been
                                     turned off

!Output Parameters:
   L1B_Gran (->sw_f_id[])     Swath file ID for each of the L1B EV files.
   L1B_Scan_t       *L1B_Scan

!Revision History:
 (continue at top of the file)

   Revision 02.01, November 14, 2001, Razor Issue #169
   Added input parameter skip_night_hi_res and changed logic so that
   250m and 500m data are not written when granule has no day mode
   scans and runtime parameter Write_Night_Mode_HiRes_Data is False.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 02.01, Dec 5, 2000, Razor issue 141.
   Added "tables" to argument list, add to arg list of "Set_L1B_EV_SDS_Attrs"
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.00 Feb 16, 1999
   Split Write_Band_Subsetting_SDS and its subordinate function
   Write_Subsetting_SDS into two functions: Set_L1B_EV_SDS_Attrs and
   Create_Band_Subsetting_SDS for clarity.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 02.00 March 1997
   Modified to match version 2 design.
   Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

   Revision 01.01 1996/04/05
   Update to match Version 1 Design Document
   John Hannon(hannon@highwire.gsfc.nasa.gov)
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
   1. The PGS_PC logical reference IDs are defined in Fnames.h.  These are
      consecutive integers for the EV files and are indexed as such below.

!END********************************************************************
*/
{
  PGSt_SMF_status       returnStatus = MODIS_S_OK;
  PGSt_PC_Logical       logical_id   = 0;
  PGSt_integer          version      = 1;
  char                  fname[PGSd_PC_FILE_PATH_MAX];
  int16  f  = 0;
  char *location = "Open_W_L1B_Granule";
  int16 start_output_res = INDEX_L1B_250m;
  
  /*
   * If we are in night mode and we do not wish to write 250m and 500m data
   * sets, start with the 1KM output data sets.
   */

  if (skip_night_hi_res == True)
  {
    start_output_res = INDEX_L1B_1km;
    logical_id = L1B_EV_1000M_FILE;
  }
  else
  {
    start_output_res = INDEX_L1B_250m;
    logical_id = L1B_EV_250M_FILE;
  }

  for (f = start_output_res; f < NUM_L1B_EV_FILES; f++, logical_id++)
  {
    version = 1;

    if (PGS_PC_GetReference(logical_id, &version, fname) != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_FILE_NOT_FOUND;
      L1BErrorMsg(location, returnStatus, 
          "Could not retrieve file name from PCF.", 
          "PGS_PC_GetReference", logical_id, NULL, True);
      return returnStatus;
    }

    if ((L1B_Gran->sw_f_id[f] = SWopen(fname, DFACC_CREATE)) == FAIL)
    {
      returnStatus = MODIS_F_FILE_NOT_CREATED;
      L1BErrorMsg(location, returnStatus, 
                  "Could not create output granule.",
                  "SWopen", logical_id, NULL, True);
      return returnStatus;
    }
  }

  returnStatus = Create_L1B_Swath (L1B_Gran, skip_night_hi_res); 
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, 
                    "Create_L1B_Swath", 0, NULL, True);
    return returnStatus;
  }

  returnStatus = Open_L1B_EV_SDS (L1B_Gran, 
                                  L1B_Scan, 
                                  skip_night_hi_res); 
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, 
                    "Open_L1B_EV_SDS", 0, NULL, True);
    return returnStatus;
  }

  returnStatus = Set_L1B_EV_SDS_Attrs (tables, 
                                       L1B_Gran, 
                                       L1B_Scan,
                                       skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, 
                "Set_L1B_EV_SDS_Attrs", 0, NULL, True);
    return returnStatus;
  }

  returnStatus = Create_Band_Subsetting_SDS (L1B_Gran,
                                             skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, 
                "Create_Band_Subsetting_SDS", 0, NULL, True);
    return returnStatus;
  }

  return(MODIS_S_OK);
}


PGSt_SMF_status Create_L1B_Swath (L1B_granule_t    *L1B_Gran,
                                  boolean          skip_night_hi_res)
/*
!C****************************************************************
!Description: Create HDFEOS Swaths and data fields for each L1B output file.
              Also get SD & V HDF interface ids. The band-subsetting data
              fields (one-dimensional arrays) are implemented as Vdata while
              the other data fields (for EV data) are implemented as SDSs.
              The values for the band-subsetting fields are assigned here.
              For the EV data, after fixing fields here, accesses will be
              opened to the SDSs in another routine and used normally
              thereafter.

!Input Parameters:
    L1B_Gran (->sw_f_id[])           HDF-EOS swath file IDs for each EV file.
    boolean         skip_night_hi_res   True if and only if all scans are 
                                     night mode scans and writing of 250m
                                     and 500m night mode granules has been
                                     turned off


!Output Parameters:
      L1B_Gran (->v_id[])       Vdata interface IDs for each EV file.
      L1B_Gran (->sd_id[])      SD interface IDs for each EV file.

!Revision History:
 (continue at top of the file)
 
 Revision 02.12, January 12, 2005 (Razor Issue #202)
 Added new file attribute prefixed with HDFEOS_FractionalOffset for each
 dimension of the swath in every MODIS product file.
 Liqin Tan  SAIC/GSO (ltan@saicmodis.com)

 Revision 02.11 November 15, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic so that
 if skip_night_hi_res is True, output data sets for 250m and 500m resolution
 data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.10 May 13, 1999
 Added band 26 section
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.01 Feb 10, 1999
 In defining the data fields for band subsetting, changed to use the global
 variables Band_subsetting_names and L1B_EV_DIM_NAME defined at the top of
 the file.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.00 April 1997
 Modified to match changes in data structures.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.00 1997/02/21
 Initial development
 Neal Devine (neal.devine@gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   Attributes not yet implemented.
   The opening of target SDS's is lengthy enough to warrant a separate sub
   function.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus;
  char     *location = "Create_L1B_Swath";
  int16    file_index  = 0; /*file index*/
  int16    R           = 0; /*L1A resolution index*/
  int16    isds        = 0; /*sds_index*/
  int32    track_dim   = 0;
  int32    frame_dim   = 0;
  int16    BG_index    = 0; /*Band_Group_index*/
  int16    num_BGs     = 0; /*num_band_groups*/
  int32    sw_id[NUM_L1B_EV_FILES] = {0, 0, 0};
  intn     hdf_return  = 0;
  int32    evfile_luns[NUM_L1B_EV_FILES] = {
            L1B_EV_250M_FILE,
            L1B_EV_500M_FILE,
            L1B_EV_1000M_FILE
  };
  int16    start_output_res = INDEX_L1B_250m;
#define MAX_ATTR_NAME_SIZE2   72
  char     attr_name_buffer[MAX_ATTR_NAME_SIZE2];

/*------------------------------------------------------------------------*/
#define SW_DEF_EV_SDS(SI_name, UI_name, str)                               \
if (SWdefdatafield(sw_id[file_index], SI_name, str, DFNT_UINT16,           \
       HDFE_NOMERGE) == FAIL)                                              \
  {                                                                        \
    returnStatus = MODIS_F_HDF_ERROR;                                      \
    L1BErrorMsg(location, returnStatus, SI_name,                           \
                "SWdefdatafield", evfile_luns[file_index], NULL, True);    \
    return returnStatus;                                                   \
  }                                                                        \
if (SWdefdatafield(sw_id[file_index], UI_name, str, DFNT_UINT8,            \
       HDFE_NOMERGE) == FAIL)                                              \
  {                                                                        \
    returnStatus = MODIS_F_HDF_ERROR;                                      \
    L1BErrorMsg(location, returnStatus, UI_name,                           \
                "SWdefdatafield", evfile_luns[file_index], NULL, True);    \
    return returnStatus;                                                   \
  }
/*------------------------------------------------------------------------*/
#define SW_DEF_EV_Aggr_SDS(SI_name, UI_name, SU_name, str)                 \
SW_DEF_EV_SDS(SI_name,UI_name,str);                                        \
if (SWdefdatafield(sw_id[file_index], SU_name, str, DFNT_INT8,             \
       HDFE_NOMERGE) == FAIL)                                              \
  {                                                                        \
    returnStatus = MODIS_F_HDF_ERROR;                                      \
    L1BErrorMsg(location, returnStatus, SU_name,                           \
                "SWdefdatafield", evfile_luns[file_index], NULL, True);    \
    return returnStatus;                                                   \
  }
/*------------------------------------------------------------------------*/

  /*
   * If we are in night mode and we do not wish to write 250m and 500m data
   * sets, start with the 1KM output data sets.
   */
  
  if (skip_night_hi_res == True)
    start_output_res = INDEX_L1B_1km;
  else
    start_output_res = INDEX_L1B_250m;
  
  /* create EV swaths */
  for (file_index = start_output_res; 
          file_index < NUM_L1B_EV_FILES; file_index++)
  {
    /* The V2 Delivery Guide suggests to use "MOD_SWATH_" 
       Here we use the name in the file spec, "MODIS_SWATH_" 
    */
    sw_id[file_index] = SWcreate(L1B_Gran->sw_f_id[file_index],
                                 "MODIS_SWATH_Type_L1B");
    if (sw_id[file_index] == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL,
                  "SWcreate", evfile_luns[file_index], NULL, True);
      return returnStatus;
    }

    /* Get equivalent HDF v_id & sd_id and assign to L1B_Gran members. 
    */
    if (EHidinfo(L1B_Gran->sw_f_id[file_index],
                 &L1B_Gran->v_id[file_index],
                 &L1B_Gran->sd_id[file_index]) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL,
                  "EHidinfo", evfile_luns[file_index], NULL, True);
      return returnStatus;
    }

    /* define band_dim and data field
    */
    num_BGs = file_index + 1;
    if (file_index == INDEX_L1B_1km)  num_BGs++; 

    for (BG_index = 0; BG_index < num_BGs; BG_index++)
    {
      R = BG_index;    /* R = [0], [0, 1], [0, 1, 2, 3] */
      hdf_return = SWdefdim(sw_id[file_index], L1B_EV_DIM_NAME[R][0],
                            (int32)L1B_BANDS_AT_RES[R]);
      if (hdf_return == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, L1B_EV_DIM_NAME[R][0],
                    "SWdefdim", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
    }

    /* The Offset and Increment entries in the HDF-EOS metadata need to be
       set appropriately to preserve the geolocation accuracy when the data
       stored has a different resolution than the geolocation (see "MODIS 
       Science Software Delivery Guide" -- 9/23/2004). Because HDF-EOS does 
       not support fractional offset values, it doesn't properly handle 
       sub-pixel offsets in Swath files. In order to better document the 
       true offset, file global attributes of fractional offset are added
       to each files (the addition to the "integer" offset provided in the 
       HDF-EOS metadata).
    */
    sprintf(attr_name_buffer,
            "HDFEOS_FractionalOffset_%s_MODIS_SWATH_Type_L1B",
            L1B_EV_DIM_NAME[num_BGs-1][1]);

    if(SDsetattr(L1B_Gran->sd_id[file_index],
                 attr_name_buffer, DFNT_FLOAT32, 1,
                 (VOIDP)&L1B_EV_DIM_FRAC_OFFSET[num_BGs-1][0]) == FAIL)
    {
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write L1B EV file dimensional attribute.",
                  "SDsetattr", evfile_luns[file_index], NULL, True);
    }

    sprintf(attr_name_buffer,
            "HDFEOS_FractionalOffset_%s_MODIS_SWATH_Type_L1B",
            L1B_EV_DIM_NAME[num_BGs-1][2]);

    if(SDsetattr(L1B_Gran->sd_id[file_index],
                 attr_name_buffer, DFNT_FLOAT32, 1,
                 (VOIDP)&L1B_EV_DIM_FRAC_OFFSET[num_BGs-1][1]) == FAIL)
    {
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write L1B EV file dimensional attribute.",
                  "SDsetattr", evfile_luns[file_index], NULL, True);
    }

    /* define track & frame dims
    */
    R = file_index;

    track_dim = DETECT_PER_BAND_AT_RES[R] * L1B_Gran->num_scans;
    frame_dim = BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;

    if (SWdefdim(sw_id[file_index], L1B_EV_DIM_NAME[R][1], 
                 track_dim) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, L1B_EV_DIM_NAME[R][1],
                  "SWdefdim", evfile_luns[file_index], NULL, True);
      return returnStatus;
    }

    if (SWdefdim(sw_id[file_index], L1B_EV_DIM_NAME[R][2], 
                 frame_dim) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, L1B_EV_DIM_NAME[R][2],
                  "SWdefdim", evfile_luns[file_index], NULL, True);
      return returnStatus;
    }

    /* define geo track & frame dims, and geo/data fields, and dimmaps
    */
    if (file_index == INDEX_L1B_1km)
    { 
      if (SWdefdim(sw_id[file_index], "2*nscans",
                   (int32)(2*L1B_Gran->num_scans)) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "2*nscans",
                    "SWdefdim", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefdim(sw_id[file_index], "1KM_geo_dim",
                   (int32)(EV_1km_FRAMES/5 + 1)) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "1KM_geo_dim",
                    "SWdefdim", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefgeofield(sw_id[file_index], "Latitude",
                        "2*nscans,1KM_geo_dim", DFNT_FLOAT32,
                        HDFE_NOMERGE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "Latitude",
                    "SWdefgeofield", evfile_luns[file_index], 
                     NULL, True);
        return returnStatus;
      }

      if (SWdefgeofield(sw_id[file_index], "Longitude",
                        "2*nscans,1KM_geo_dim", DFNT_FLOAT32,
                        HDFE_NOMERGE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "Longitude",
                    "SWdefgeofield", evfile_luns[file_index], 
                    NULL, True);
        return returnStatus;
      }

      SW_DEF_EV_SDS(L1B_EV_SDS_NAME[SI_1km],
                    L1B_EV_SDS_NAME[UI_1km],
                    "Band_1KM_RefSB,10*nscans,Max_EV_frames");

      SW_DEF_EV_SDS(L1B_EV_SDS_NAME[SI_1km_EMISS],
                    L1B_EV_SDS_NAME[UI_1km_EMISS],
                    "Band_1KM_Emissive,10*nscans,Max_EV_frames");

      SW_DEF_EV_Aggr_SDS(L1B_EV_SDS_NAME[SI_250to1km],
                         L1B_EV_SDS_NAME[UI_250to1km],
                         L1B_EV_SDS_NAME[SU_250to1km],
                         "Band_250M,10*nscans,Max_EV_frames");

      SW_DEF_EV_Aggr_SDS(L1B_EV_SDS_NAME[SI_500to1km],
                         L1B_EV_SDS_NAME[UI_500to1km],
                         L1B_EV_SDS_NAME[SU_500to1km],
                         "Band_500M,10*nscans,Max_EV_frames");

      for (isds = INDEX_HEIGHT; isds < NUM_GEO_SDS; isds++)
      {
        if (SWdefdatafield(sw_id[file_index], GEO_SDS[isds].name,
                           "2*nscans,1KM_geo_dim",
                           GEO_SDS[isds].type,HDFE_NOMERGE) == FAIL)
        {
          returnStatus = MODIS_F_HDF_ERROR;
          L1BErrorMsg(location, returnStatus, GEO_SDS[isds].name,
                      "SWdefdatafield", evfile_luns[file_index], NULL, True);
          return returnStatus;
        }
      }
      
      if (SWdefdimmap(sw_id[file_index], "2*nscans", "10*nscans",
                      GEO_OFFSET, GEO_STRIDE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "2*nscans to 10*nscans",
                    "SWdefdimmap", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefdimmap(sw_id[file_index], "1KM_geo_dim", "Max_EV_frames",
                      GEO_OFFSET, GEO_STRIDE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "1KM_geo_dim to Max_EV_frames",
                    "SWdefdimmap", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      /* define data fields for band_subsetting
      */

      for (R = 0; R < 4; R++) {
        if (SWdefdatafield(sw_id[file_index],
                           Band_subsetting_names[R],
                           L1B_EV_DIM_NAME[R][0],
                           DFNT_FLOAT32, HDFE_NOMERGE) == FAIL)
        {
          returnStatus = MODIS_F_HDF_ERROR;
          L1BErrorMsg(location, returnStatus, Band_subsetting_names[R],
                      "SWdefdatafield", evfile_luns[file_index], NULL, True);
          return returnStatus;
        }
      }

    }
    else
    {
      if (SWdefdim(sw_id[file_index], "10*nscans",
                   (int32)(10*L1B_Gran->num_scans)) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "10*nscans",
                    "SWdefdim", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefdim(sw_id[file_index], "Max_EV_frames",
                   (int32)EV_1km_FRAMES) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "Max_EV_frames",
                    "SWdefdim", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefgeofield(sw_id[file_index], "Latitude",
                        "10*nscans,Max_EV_frames", DFNT_FLOAT32,
                        HDFE_NOMERGE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "Latitude",
                    "SWdefgeofield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefgeofield(sw_id[file_index], "Longitude",
                        "10*nscans,Max_EV_frames", DFNT_FLOAT32,
                        HDFE_NOMERGE) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, "Longitude",
                    "SWdefgeofield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      switch(file_index)
      {
        case INDEX_L1B_250m:
          SW_DEF_EV_SDS(L1B_EV_SDS_NAME[SI_250m],
                        L1B_EV_SDS_NAME[UI_250m],
                        "Band_250M,40*nscans,4*Max_EV_frames");

          if (SWdefdatafield(sw_id[file_index],
                             Band_subsetting_names[INDEX_250M],
                             L1B_EV_DIM_NAME[INDEX_250M][0],
                             DFNT_FLOAT32, HDFE_NOMERGE) == FAIL)
          {
            returnStatus = MODIS_F_HDF_ERROR;
            L1BErrorMsg(location, returnStatus, 
                        Band_subsetting_names[INDEX_250M],
                        "SWdefdatafield", 
                        evfile_luns[file_index], NULL, True);
            return returnStatus;
          }

          break;

        case INDEX_L1B_500m:
          SW_DEF_EV_SDS(L1B_EV_SDS_NAME[SI_500m],
                        L1B_EV_SDS_NAME[UI_500m],
                        "Band_500M,20*nscans,2*Max_EV_frames");

          SW_DEF_EV_Aggr_SDS(L1B_EV_SDS_NAME[SI_250to500m],
                             L1B_EV_SDS_NAME[UI_250to500m],
                             L1B_EV_SDS_NAME[SU_250to500m],
                             "Band_250M,20*nscans,2*Max_EV_frames");

          for (R = 0; R < 2; R++) {
            if (SWdefdatafield(sw_id[file_index],
                               Band_subsetting_names[R],
                               L1B_EV_DIM_NAME[R][0],
                               DFNT_FLOAT32, HDFE_NOMERGE) == FAIL)
            {
              returnStatus = MODIS_F_HDF_ERROR;
              L1BErrorMsg(location, returnStatus, 
                          Band_subsetting_names[R],
                          "SWdefdatafield", 
                          evfile_luns[file_index], NULL, True);
              return returnStatus;
            }
          }

      }

      R = file_index;
     
      if (SWdefdimmap(sw_id[file_index], 
                      "10*nscans", 
                      L1B_EV_DIM_NAME[R][1], 
                      L1B_EV_DIM_OFFSET[R][0], BAND_RATIO_AT_RES[R]) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, 
                    "L1B_EV_DIM_NAME[R][1] to 10*nscans",
                    "SWdefdimmap", 
                    evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWdefdimmap(sw_id[file_index], 
                      "Max_EV_frames", 
                      L1B_EV_DIM_NAME[R][2],
                      L1B_EV_DIM_OFFSET[R][1], BAND_RATIO_AT_RES[R]) == FAIL)
      {
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, 
                    "L1B_EV_DIM_NAME[R][2] to 10*nscans",
                    "SWdefdimmap", 
                    evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
    }

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
    if (file_index == INDEX_L1B_1km)
    { 
        /*
         * Create the band 26 SI and UI SDSs as part of the HDF-EOS
         * swath.  These are 2D SDSs, using the dimensions indicated.
         * These are only created for the L1B 1km EV product.
         */
      SW_DEF_EV_SDS(BAND_26_SI_SDS_NAME,
                    BAND_26_UI_SDS_NAME,
                    "10*nscans,Max_EV_frames");
    }
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  }

  for (file_index = start_output_res; 
            file_index < NUM_L1B_EV_FILES; file_index++)
  {

       /* detach from swath (this should be done prior to accessing members
        * to ensure that fields are set properly).
        */

    if (SWdetach(sw_id[file_index]) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL,
                  "SWdetach", evfile_luns[file_index], NULL, True);
      return returnStatus;
    }

      /* Write the band-subsetting fields numerical values to file.
       */

    returnStatus = Write_Swath_Band_Number(file_index, L1B_Gran);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Write_Swath_Band_Number", 
                  evfile_luns[file_index], NULL, True);
      return returnStatus;
    }
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Write_Swath_Band_Number(int32 file_index,
                                        L1B_granule_t *L1B_Gran)
/*
!C**********************************************************************
!Description:   For one EV file, this function writes the numerical values
                of the band-subsetting SDS data fields.
   
!Input Parameters:
    file_index              index of the EV file that we will write data to.
    L1B_Gran (->sw_f_id[])  array of swath file IDs for the EV files.

!Output Parameters:
      
!Revision History:
 (continue at top of the file)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

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
  PGSt_SMF_status   returnStatus = MODIS_S_OK;
  char *location = "Write_Swath_Band_Number";

  int32 swId;
  float32 data_250m[2] = {1, 2};
  float32 data_500m[5] = {3, 4, 5, 6, 7};
  float32 data_1km_refl[15] = 
            {8, 9, 10, 11, 12, 13, 13.5, 14, 14.5,
             15, 16, 17, 18, 19, 26};
  float32 data_1km_emiss[16] = 
            {20, 21, 22, 23, 24, 25, 27, 28, 29,
             30, 31, 32, 33, 34, 35, 36};
  int32 evfile_luns[NUM_L1B_EV_FILES] = 
            {
              L1B_EV_250M_FILE,
              L1B_EV_500M_FILE,
              L1B_EV_1000M_FILE
            };

  swId = SWattach(L1B_Gran->sw_f_id[file_index], "MODIS_SWATH_Type_L1B");
  if (swId == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL,
                "SWattach", evfile_luns[file_index], NULL, True);
    return returnStatus;
  }

  switch(file_index)
  {
    case INDEX_L1B_1km:
      if (SWwritefield(swId, "Band_250M", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_250m) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_250M",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWwritefield(swId, "Band_500M", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_500m) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_500M",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWwritefield(swId, "Band_1KM_RefSB", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_1km_refl)== FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_1KM_RefSB",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }

      if (SWwritefield(swId, "Band_1KM_Emissive", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_1km_emiss)== FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_1KM_Emissive",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
      break;

    case INDEX_L1B_500m: 
      if (SWwritefield(swId, "Band_250M", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_250m) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_250M",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
      
      if (SWwritefield(swId, "Band_500M", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_500m) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_500M",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
      break;

    case INDEX_L1B_250m: 
      if (SWwritefield(swId, "Band_250M", (int32*)NULL, (int32*)NULL,
                       (int32*)NULL, (VOIDP)data_250m) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, "Band_250M",
                    "SWwritefield", evfile_luns[file_index], NULL, True);
        return returnStatus;
      }
      break; 
  }

  if (SWdetach(swId) == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL,
                "SWdetach", evfile_luns[file_index], NULL, True);
    return returnStatus;
  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Open_L1B_EV_SDS (L1B_granule_t    *L1B_Gran,
                                 L1B_Scan_t       *L1B_Scan,
                                 boolean          skip_night_hi_res)
/*
!C****************************************************************
!Description: Given the L1B file swaths and their EV data fields created
              earlier (implemented as SDSs), open SDS access for each of the
              SDSs.  Additionally, set the "long name" attribute for each SDS.
              SDS access will remain open until after all scan-by-scan
              processing is complete.

!Input Parameters:
    L1B_Gran (->sd_id)               SD file IDs for each of the L1B EV files.
    boolean         skip_night_hi_res   True if and only if all scans are 
                                     night mode scans and writing of 250m
                                     and 500m night mode granules has been
                                     turned off
!Output Parameters:
      L1B_Scan              SDS access ID members:
                               SI_sds_id[]
                               EV_250m_Aggr500m_RefSB_sds_id
                               EV_250m_Aggr1km_RefSB_sds_id
                               EV_500m_Aggr1km_RefSB_sds_id
                               UI_sds_id[]
                               EV_250m_Aggr500m_RefSB_UI_sds_id
                               EV_250m_Aggr1km_RefSB_UI_sds_id
                               EV_500m_Aggr1km_RefSB_UI_sds_id
                               EV_250m_Aggr500m_RefSB_SU_sds_id
                               EV_250m_Aggr1km_RefSB_SU_sds_id
                               EV_500m_Aggr1km_RefSB_SU_sds_id

!Revision History:
 (continue at top of the file)

 Revision 03.11 November 15, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic so that
   if skip_night_hi_res is True, output data sets for 250m and 500m resolution
   data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 03.10 May   1999
 Added Band 26 SDS.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 03.00 April 1997
 Modified to match changes in data structures.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 02.0 1996/06/28
 Version 2, initial redesign
 Neal Devine(neal.devine@gsfc.nasa.gov)

 Revision 01.02 1996/04/05
 Update to match Version 1 Design Document
 Neal Devine(neal.devine@gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)

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
  PGSt_SMF_status returnStatus;
  char *location = "Open_L1B_EV_SDS";
  int32    sds_index    = 0;
  int16    file_index  = 0;
  int32 evfile_luns[NUM_L1B_EV_FILES] = {
    L1B_EV_250M_FILE,
    L1B_EV_500M_FILE,
    L1B_EV_1000M_FILE
  };

  int32 scaled_250[2] = {13,16};
  int32 scaled_500[2] = {10,12};

  /*
   * ----------------------------------------------------------------------
   * This macro opens SDS access to a pair of L1B EV SDS's (Scaled_Integer 
   * and Uncert_Index) and sets the long-name attribute for each SDS.
   *   Input:  sd_id, SI_name, UI_name, SI_longName, UI_longName
   *   Output: SI_id, UI_id
   * ----------------------------------------------------------------------
   */
#define CREATE_EV_SDS(sd_id, SI_name, UI_name, SI_longName, UI_longName,    \
        SI_id,UI_id)                                                        \
if((sds_index = SDnametoindex(sd_id, SI_name)) == FAIL)                     \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SI_name,                            \
                "SDnametoindex", evfile_luns[file_index], NULL, True);      \
    return returnStatus;                                                    \
  }                                                                         \
if((SI_id = SDselect(sd_id, sds_index)) == FAIL)                            \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SI_name,                            \
                "SDselect", evfile_luns[file_index], NULL, True);           \
    return returnStatus;                                                    \
  }                                                                         \
if(SDsetattr(SI_id, "long_name", DFNT_CHAR8,(int32)strlen(SI_longName),     \
             SI_longName) == FAIL)                                          \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SI_name,                            \
                "SDsetattr", evfile_luns[file_index], "long_name", True);   \
    return returnStatus;                                                    \
  }                                                                         \
if((sds_index = SDnametoindex(sd_id, UI_name)) == FAIL)                     \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, UI_name,                            \
                "SDnametoindex", evfile_luns[file_index], NULL, True);      \
    return returnStatus;                                                    \
  }                                                                         \
if((UI_id = SDselect(sd_id, sds_index)) == FAIL)                            \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, UI_name,                            \
                "SDselect", evfile_luns[file_index], NULL, True);           \
    return returnStatus;                                                    \
  }                                                                         \
if(SDsetattr(UI_id, "long_name", DFNT_CHAR8,(int32)strlen(UI_longName),     \
             UI_longName) == FAIL)                                          \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, UI_name,                            \
                "SDsetattr", evfile_luns[file_index], "long_name", True);   \
    return returnStatus;                                                    \
  } 

  /*
   * -----------------------------------------------------------------------
   * This macro opens SDS access to a set of three L1B EV_Aggr SDS's 
   * (SI, UI and SU) and sets the long-name attribute for each SDS
   *   Input:  sd_id, SI_name, UI_name, SU_name, SI_LName, UI_LName, SU_LName
   *   Output: SI_id, UI_id, SU_id
   * -----------------------------------------------------------------------
   */
#define CREATE_EV_Aggr_SDS(sd_id, SI_name, UI_name, SU_name, SI_LName,      \
                           UI_LName, SU_LName, SI_id, UI_id, SU_id)         \
CREATE_EV_SDS(sd_id, SI_name, UI_name, SI_LName, UI_LName, SI_id, UI_id);   \
if((sds_index = SDnametoindex(sd_id, SU_name)) == FAIL)                     \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SU_name,                            \
                "SDnametoindex", evfile_luns[file_index], NULL, True);      \
    return returnStatus;                                                    \
  }                                                                         \
if((SU_id = SDselect(sd_id, sds_index)) == FAIL)                            \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SU_name,                            \
                "SDselect", evfile_luns[file_index], NULL, True);           \
    return returnStatus;                                                    \
  }                                                                         \
if(SDsetattr(SU_id, "long_name", DFNT_CHAR8,(int32)strlen(SU_LName),        \
             SU_LName) == FAIL)                                             \
  {                                                                         \
    returnStatus = MODIS_F_HDF_ERROR;                                       \
    L1BErrorMsg(location, returnStatus, SU_name,                            \
                "SDsetattr", evfile_luns[file_index], "long_name", True);   \
    return returnStatus;                                                    \
  }
/*--------------------------------------------------------------------*/

  /* Create 250m and 500m files only when skip_night_hi_res is False. */
  
  if (skip_night_hi_res == False)
  {
    /*250m file*/
    file_index = INDEX_L1B_250m;

    CREATE_EV_SDS (L1B_Gran->sd_id[INDEX_L1B_250m],
                   L1B_EV_SDS_NAME[SI_250m],
                   L1B_EV_SDS_NAME[UI_250m],
                   L1B_EV_SDS_LONG_NAME[SI_250m],
                   L1B_EV_SDS_LONG_NAME[UI_250m],
                   L1B_Scan->SI_sds_id[INDEX_250M],
                   L1B_Scan->UI_sds_id[INDEX_250M]);

    /*500m file*/
    file_index = INDEX_L1B_500m;
    CREATE_EV_SDS (L1B_Gran->sd_id[INDEX_L1B_500m],
                   L1B_EV_SDS_NAME[SI_500m],
                   L1B_EV_SDS_NAME[UI_500m],
                   L1B_EV_SDS_LONG_NAME[SI_500m],
                   L1B_EV_SDS_LONG_NAME[UI_500m],
                   L1B_Scan->SI_sds_id[INDEX_500M],
                   L1B_Scan->UI_sds_id[INDEX_500M]);

    CREATE_EV_Aggr_SDS (L1B_Gran->sd_id[INDEX_L1B_500m],
                        L1B_EV_SDS_NAME[SI_250to500m],
                        L1B_EV_SDS_NAME[UI_250to500m],
                        L1B_EV_SDS_NAME[SU_250to500m],
                        L1B_EV_SDS_LONG_NAME[SI_250to500m],
                        L1B_EV_SDS_LONG_NAME[UI_250to500m],
                        L1B_EV_SDS_LONG_NAME[SU_250to500m],
                        L1B_Scan->EV_250m_Aggr500m_RefSB_sds_id,
                        L1B_Scan->EV_250m_Aggr500m_RefSB_UI_sds_id,
                        L1B_Scan->EV_250m_Aggr500m_RefSB_SU_sds_id);
  }
  
  /*1km file*/
  file_index = INDEX_L1B_1km;
  CREATE_EV_SDS (L1B_Gran->sd_id[INDEX_L1B_1km],
                 L1B_EV_SDS_NAME[SI_1km],
                 L1B_EV_SDS_NAME[UI_1km],
                 L1B_EV_SDS_LONG_NAME[SI_1km],
                 L1B_EV_SDS_LONG_NAME[UI_1km],
                 L1B_Scan->SI_sds_id[INDEX_1000M_REFL],
                 L1B_Scan->UI_sds_id[INDEX_1000M_REFL]);

  CREATE_EV_Aggr_SDS (L1B_Gran->sd_id[INDEX_L1B_1km],
                      L1B_EV_SDS_NAME[SI_250to1km],
                      L1B_EV_SDS_NAME[UI_250to1km],
                      L1B_EV_SDS_NAME[SU_250to1km],
                      L1B_EV_SDS_LONG_NAME[SI_250to1km],
                      L1B_EV_SDS_LONG_NAME[UI_250to1km],
                      L1B_EV_SDS_LONG_NAME[SU_250to1km],
                      L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id,
                      L1B_Scan->EV_250m_Aggr1km_RefSB_UI_sds_id,
                      L1B_Scan->EV_250m_Aggr1km_RefSB_SU_sds_id);

  CREATE_EV_Aggr_SDS (L1B_Gran->sd_id[INDEX_L1B_1km],
                      L1B_EV_SDS_NAME[SI_500to1km],
                      L1B_EV_SDS_NAME[UI_500to1km],
                      L1B_EV_SDS_NAME[SU_500to1km],
                      L1B_EV_SDS_LONG_NAME[SI_500to1km],
                      L1B_EV_SDS_LONG_NAME[UI_500to1km],
                      L1B_EV_SDS_LONG_NAME[SU_500to1km],
                      L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id,
                      L1B_Scan->EV_500m_Aggr1km_RefSB_UI_sds_id,
                      L1B_Scan->EV_500m_Aggr1km_RefSB_SU_sds_id);

  CREATE_EV_SDS (L1B_Gran->sd_id[INDEX_L1B_1km],
                 L1B_EV_SDS_NAME[SI_1km_EMISS],
                 L1B_EV_SDS_NAME[UI_1km_EMISS],
                 L1B_EV_SDS_LONG_NAME[SI_1km_EMISS],
                 L1B_EV_SDS_LONG_NAME[UI_1km_EMISS],
                 L1B_Scan->SI_sds_id[INDEX_1000M_EMISS],
                 L1B_Scan->UI_sds_id[INDEX_1000M_EMISS]);
  
  if ((RSCL_FLAG & 1) == 1)
    {
      if (SDsetattr(L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id,
                    "Rescaled Ocean R", DFNT_INT32,
                    1, &scaled_250[0]) == FAIL)
        {
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, "Rescaled Ocean R",
                      "SDsetattr", 0, NULL, False);
          return returnStatus;
        }
    }

  if ((RSCL_FLAG & 2) == 2)
    {
      if (SDsetattr(L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id,
                    "Rescaled Ocean NIR", DFNT_INT32,
                    1, &scaled_250[1]) == FAIL)
        {
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, "Rescaled Ocean NIR",
                      "SDsetattr", 0, NULL, False);
          return returnStatus;
        }
    }

  if ((RSCL_FLAG & 4) == 4)
    {
      if (SDsetattr(L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id,
                    "Rescaled Ocean B", DFNT_INT32,
                    1, &scaled_500[0]) == FAIL)
        {
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, "Rescaled Ocean B",
                      "SDsetattr", 0, NULL, False);
          return returnStatus;
        }
    }

  if ((RSCL_FLAG & 8) == 8)
    {
      if (SDsetattr(L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id,
                    "Rescaled Ocean G", DFNT_INT32,
                    1, &scaled_500[1]) == FAIL)
        {
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, "Rescaled Ocean G",
                      "SDsetattr", 0, NULL, False);
          return returnStatus;
        }
    }

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
/*
 * Open SDS access to the band 26 SDSs and save the
 * SDS ids in the L1B_Scan structure.  SDS access is closed in function
 * Close_L1B_Granule.  Note that the SDSs are not really created here
 * depite the name of the macro.  They were created in Create_L1B_Swath.
 */
  CREATE_EV_SDS (L1B_Gran->sd_id[INDEX_L1B_1km],
                 BAND_26_SI_SDS_NAME,
                 BAND_26_UI_SDS_NAME,
                 BAND_26_SI_SDS_LONG_NAME,
                 BAND_26_UI_SDS_LONG_NAME,
                 L1B_Scan->Band26.SI_sds_id,
                 L1B_Scan->Band26.UI_sds_id);
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  return(MODIS_S_OK);
}


PGSt_SMF_status Set_L1B_EV_SDS_Attrs (lookup_tables_t  *tables,
                                      L1B_granule_t    *L1B_Gran,
                                      L1B_Scan_t       *L1B_Scan,
                                      boolean          skip_night_hi_res)
/*
!C****************************************************************
!Description: This routine creates and set the units, valid range and
              fill-value attributes for each of the EV SDSs (SI, UI and SU). 
              This routine also writes the following attributes of scaled
              integer SDSs: (1) band names, (2) reflectance scales, offsets
              and units (only if reflectance SDSs are present), (3) radiance
              scales, offsets and units and (4) corrected counts scales,
              offsets and units.

!Input Parameters:
    L1B_Gran  ->sd_id[]       SD file ids for the EV files.
              ->SO            Scale-Offset struct, contains the radiance,
                              reflectance and corrected counts scales,
                                offsets.
    L1B_Scan                  SDS access ids used in Get_SDS_id and in
                              Set_Unit_Range_Fillvalue.
    boolean  skip_night_hi_res   True if and only if all scans are
                              night mode scans and writing of 250m
                              and 500m night mode granules has been
                              turned off
!Output Parameters:
      None

!Revision History:
 (continue at top of the file)

 Revision 02.21 November 15, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic so that
   if skip_night_hi_res is True, output data sets for 250m and 500m resolution
   data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.20 May 13, 1999
 Added band 26 section.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.11 Sept 14 1998
 Change dimension attribute to sds attribute based on
 file specs V2.1.1
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.10 April 10 1998
 Changed to the V2.1 L1B_Gran->momos (from the V2 
 L1B_Gran->SO.rad_scale_RefSB and L1B_Gran->SO.refl_scale_RefSB)
 David Catozzi (cato@ltpmail.gsfc.nasa.gov)

 Revision 02.00 April 1997
 Modified to match changes in data structures.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.00 1997/03/01  
 Initial development
 Neal Devine (neal.devine@gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
  Note that the band dim sds's are used for two distinct purposes:
  dim attributes, which happen to create float32 sdss, and band subsetting
  numbers.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus;
  char *location = "Set_L1B_EV_SDS_Attrs";
  int16  f    = 0; /*L1B EV file index (0,2)*/
  int16  loc  = 0; /*location in an array*/
  int16  num_sds    = 0; /*num band subsetting sds's within one file*/
  int32  sds_ids[4] = {0, 0, 0, 0};  /* added by gu */

  typedef enum {
    radiance,
    reflectance,
    counts,
    Num_Kinds
  } ScaleOffset_Kind_t;
 
  float32 *scale[Num_Kinds][NUM_L1A_RESOLUTIONS];
  float32 *offset[Num_Kinds][NUM_L1A_RESOLUTIONS];
  char    *BandNames[NUM_L1A_RESOLUTIONS] = {
           "1,2",
           "3,4,5,6,7",
           "8,9,10,11,12,13lo,13hi,14lo,14hi,15,16,17,18,19,26",
           "20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36"
  };
  char    *rad_units    = "Watts/m^2/micrometer/steradian";
  char    *refl_units   = "none";
  char    *counts_units = "counts";
  int32   evfile_luns[NUM_L1B_EV_FILES] = 
                                          {
                                            L1B_EV_250M_FILE,
                                            L1B_EV_500M_FILE,
                                            L1B_EV_1000M_FILE
                                          };
  int16 start_output_res = INDEX_L1B_250m;

  returnStatus = Set_Unit_Range_Fillvalue(L1B_Scan, skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Set_Unit_Range_Fillvalue", 0, NULL, True);
    return returnStatus;
  }

  /*
   * Add additional attributes to uncertainty index SDSs.
   */

  returnStatus = Set_UI_ConvertToPercent_Attrs(tables, 
                                               L1B_Scan, 
                                               skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Set_UI_ConvertToPercent_Attrs", 0, NULL, True);
    return returnStatus;
  }
 
    /********************************
     *   250m scales and offsets    *
     ********************************/
    scale[radiance][INDEX_250M]     = &(L1B_Gran->SO.rad_scale_RefSB[0]);
    scale[reflectance][INDEX_250M]  = &(L1B_Gran->SO.refl_scale_RefSB[0]);

    offset[radiance][INDEX_250M]    = &(L1B_Gran->SO.rad_offset_RefSB[0]);
    offset[reflectance][INDEX_250M] = &(L1B_Gran->SO.refl_offset_RefSB[0]);

    scale[counts][INDEX_250M]       = &(L1B_Gran->SO.counts_scale_RefSB[0]);
    offset[counts][INDEX_250M]      = &(L1B_Gran->SO.counts_offset_RefSB[0]);


    /********************************
     *   500m scales and offsets    *
     ********************************/
    loc = L1B_BANDS_AT_RES[INDEX_250M];

    scale[radiance][INDEX_500M]     = 
        &(L1B_Gran->SO.rad_scale_RefSB[loc]);
    scale[reflectance][INDEX_500M]  = 
        &(L1B_Gran->SO.refl_scale_RefSB[loc]);

    offset[radiance][INDEX_500M]    = 
        &(L1B_Gran->SO.rad_offset_RefSB[loc]);
    offset[reflectance][INDEX_500M] = 
        &(L1B_Gran->SO.refl_offset_RefSB[loc]);

    scale[counts][INDEX_500M]       = 
        &(L1B_Gran->SO.counts_scale_RefSB[loc]);
    offset[counts][INDEX_500M]      = 
        &(L1B_Gran->SO.counts_offset_RefSB[loc]);

  /********************************
   * 1Km Refl scales and offsets  *
   ********************************/
  loc += L1B_BANDS_AT_RES[INDEX_500M];

  scale[radiance][INDEX_1000M_REFL]     = 
        &(L1B_Gran->SO.rad_scale_RefSB[loc]);
  scale[reflectance][INDEX_1000M_REFL]  = 
        &(L1B_Gran->SO.refl_scale_RefSB[loc]);

  offset[radiance][INDEX_1000M_REFL]    = 
        &(L1B_Gran->SO.rad_offset_RefSB[loc]);
  offset[reflectance][INDEX_1000M_REFL] = 
        &(L1B_Gran->SO.refl_offset_RefSB[loc]);

  scale[counts][INDEX_1000M_REFL]       = 
        &(L1B_Gran->SO.counts_scale_RefSB[loc]);
  offset[counts][INDEX_1000M_REFL]      = 
        &(L1B_Gran->SO.counts_offset_RefSB[loc]);


  /********************************
   * 1Km Emiss scales and offsets *
   ********************************/
  scale[radiance][INDEX_1000M_EMISS]  = 
        &(L1B_Gran->SO.rad_scale_Emiss[0]);
  offset[radiance][INDEX_1000M_EMISS] = 
        &(L1B_Gran->SO.rad_offset_Emiss[0]);

  /*******************************
   * 250m holds rescaled R & NIR *
   *******************************/
  if ((RFLAG & 1) == 1) {
    scale[reflectance][INDEX_250M][0] = scale[counts][INDEX_250M][0]*
      scale[reflectance][INDEX_1000M_REFL][5]/scale[counts][INDEX_1000M_REFL][5];
    scale[reflectance][INDEX_250M][1] = scale[counts][INDEX_250M][1]*
      scale[reflectance][INDEX_1000M_REFL][10]/scale[counts][INDEX_1000M_REFL][10];
  }

  /*****************************
   * 500m holds rescaled B & G *
   *****************************/
  if ((RFLAG & 2) == 2) {
    scale[reflectance][INDEX_500M][0] = scale[counts][INDEX_500M][0]*
      scale[reflectance][INDEX_1000M_REFL][2]/scale[counts][INDEX_1000M_REFL][2];
    scale[reflectance][INDEX_500M][1] = scale[counts][INDEX_500M][1]*
      scale[reflectance][INDEX_1000M_REFL][4]/scale[counts][INDEX_1000M_REFL][4];
  }

  /*
   * If we are in night mode and we do not wish to write 250m and 500m data
   * sets, start with the 1KM output data sets.
   */

  if (skip_night_hi_res == True)
    start_output_res = INDEX_L1B_1km;
  else
    start_output_res = INDEX_L1B_250m;
 
  for (f = start_output_res; f < NUM_L1B_EV_FILES; f++)
  {
    returnStatus = Get_SDS_id(f, L1B_Scan, &num_sds, sds_ids);

    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Get_SDS_id", evfile_luns[f], NULL, True);
      return returnStatus;
    }

    returnStatus = Set_SDS_Attributes(sds_ids, BandNames, 
                                      (float32 **)scale, 
                                      (float32 **)offset, 
                                      rad_units,
                                      refl_units, 
                                      counts_units, 
                                      num_sds);

    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Set_SDS_Attributes", evfile_luns[f], NULL, True);
      return returnStatus;
    }
   
  } /*f*/

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
    /*
     * Set the radiance, reflectance and corrected counts scales for
     * the Band 26 scaled integer SDS.
     */

  f   = INDEX_L1B_1km;             /* file index, for LUN identification */
  loc = NUM_REFLECTIVE_BANDS - 1;  /* index of Band 26 in SO arrays */

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "radiance_scales",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.rad_scale_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"radiance_scales\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "radiance_offsets",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.rad_offset_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"radiance_offsets\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "radiance_units",
                DFNT_CHAR8,
                (int32)strlen(rad_units), rad_units) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"radiance_units\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "reflectance_scales",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.refl_scale_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"reflectance_scales\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "reflectance_offsets",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.refl_offset_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"reflectance_offsets\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "reflectance_units",
                DFNT_CHAR8,
                (int32)strlen(refl_units), refl_units) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"reflectance_units\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "corrected_counts_scales",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.counts_scale_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"corrected_counts_scales\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "corrected_counts_offsets",
                DFNT_FLOAT32,
                1, &(L1B_Gran->SO.counts_offset_RefSB[loc])) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"corrected_counts_offsets\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "corrected_counts_units",
                DFNT_CHAR8,
                (int32)strlen(counts_units), counts_units) == FAIL)
  {
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus,
                "Could not write Band 26 SI SDS "
                   "attribute "
                "\"corrected_counts_units\".",
                "SDsetattr", evfile_luns[f], NULL, True);
    return returnStatus;
  }

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  return(MODIS_S_OK);
}

PGSt_SMF_status Set_Unit_Range_Fillvalue(L1B_Scan_t *L1B_Scan,
                                         boolean    skip_night_hi_res)
/*
!C**********************************************************************
!Description:   This routine creates and set the units, valid range and
                fill-value attributes for each of the EV SDSs (SI, UI and SU).
   
!Input Parameters:
    L1B_Scan_t *  L1B_Scan   
    boolean  skip_night_hi_res   True if and only if all scans are
                              night mode scans and writing of 250m
                              and 500m night mode granules has been
                              turned off

!Output Parameters:
      
!Revision History:
 (continue at top of the file)

 Revision 01.11  October 16, 2004  Razor Issue #200
 Casted Int32 variables in sprintf calls to "long" with the
 format specifier "%ld" for better code portability.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 01.10  May 13, 1999
 Added band 26 section.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01  Feb 18, 1999
 Made the SU_valid_range an array of pointers, since there are different
 valid ranges for the different aggregated SDSs.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

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
  PGSt_SMF_status   returnStatus = MODIS_S_OK;
  char   *location = "Set_Unit_Range_Fillvalue";
  char   errmsg[256];
  char   *errmsgfmt = "Could not write attribute \"%s\", "
                      "type[%ld], sds_number[%ld].";
  typedef enum {
                SI,
                UI,
                SU,
                Num_Data_types
               } EV_Data_Type_t;
  typedef enum { RefSB_250_Data,
                 RefSB_500_Data,
                 RefSB_1KM_Data,
                 Emissive_Data,
                 Aggr_250_500_RefSB_Data,
                 Aggr_250_1km_RefSB_Data,
                 Aggr_500_1km_RefSB_Data,
                 NUM_EV_DATA_SETS
               } ev_Data_Sets_t;
  typedef enum { Aggr_250_500_RefSB_Samples,
                 Aggr_250_1km_RefSB_Samples,
                 Aggr_500_1km_RefSB_Samples,
                 NUM_EV_AGGR_SAMPLES
               } ev_Aggr_Samples_t;      
  int32  type;
  int32  sds_number;
  int32  num_sds_ids[Num_Data_types] = {7, 7, 3};
  int32  R;
  char   *units = "none";
  uint16 SI_valid_range[2] = {0, 32767};
  uint8  UI_valid_range[2] = {0, 15};
  int8   SU_valid_range_X2[2] = {0, 6};  /* aggregating by a factor of 2 */
  int8   SU_valid_range_X4[2] = {0, 28}; /* aggregating by a factor of 4 */
  int8   *SU_valid_range[3];
  uint16 SI_fillvalue = 65535;
  uint8  UI_fillvalue = 255;
  int8   SU_fillvalue = -1;  
  int32  SI_sds_ids[7];
  int32  UI_sds_ids[7];
  int32  SU_sds_ids[3];
  int32  *sds_ids[Num_Data_types];
  int16  start_output_res = INDEX_L1B_250m;
     /* in the loop below, these are indexed with sds_number for SU case */

                  /* 250m_Aggr500 */
  SU_valid_range[Aggr_250_500_RefSB_Samples] = SU_valid_range_X2;  
                  /* 250m_Aggr1km */
  SU_valid_range[Aggr_250_1km_RefSB_Samples] = SU_valid_range_X4;  
                  /* 500m_Aggr1km */
  SU_valid_range[Aggr_500_1km_RefSB_Samples] = SU_valid_range_X2;  
  
  sds_ids[0] = SI_sds_ids;
  sds_ids[1] = UI_sds_ids;
  sds_ids[2] = SU_sds_ids;

  for(R = 0; R < NUM_L1A_RESOLUTIONS; R++)
  {   
    SI_sds_ids[R] = L1B_Scan->SI_sds_id[R];
    UI_sds_ids[R] = L1B_Scan->UI_sds_id[R];
  }
  
  SI_sds_ids[Aggr_250_500_RefSB_Data] = 
                              L1B_Scan->EV_250m_Aggr500m_RefSB_sds_id;
  SI_sds_ids[Aggr_250_1km_RefSB_Data] = 
                              L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id;
  SI_sds_ids[Aggr_500_1km_RefSB_Data] = 
                              L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id;
  
  UI_sds_ids[Aggr_250_500_RefSB_Data] = 
                              L1B_Scan->EV_250m_Aggr500m_RefSB_UI_sds_id;
  UI_sds_ids[Aggr_250_1km_RefSB_Data] = 
                              L1B_Scan->EV_250m_Aggr1km_RefSB_UI_sds_id;
  UI_sds_ids[Aggr_500_1km_RefSB_Data] = 
                              L1B_Scan->EV_500m_Aggr1km_RefSB_UI_sds_id;

  SU_sds_ids[Aggr_250_500_RefSB_Samples] = 
                              L1B_Scan->EV_250m_Aggr500m_RefSB_SU_sds_id;
  SU_sds_ids[Aggr_250_1km_RefSB_Samples] = 
                              L1B_Scan->EV_250m_Aggr1km_RefSB_SU_sds_id;
  SU_sds_ids[Aggr_500_1km_RefSB_Samples] = 
                              L1B_Scan->EV_500m_Aggr1km_RefSB_SU_sds_id;

  /*
   * The following shows the relationship between index, type, and sds#.
   * index type sds#        good night indices  = [2, 3, 5, 6]
   *   0     0   0     SI_name = EV_250_RefSB
   *   1     0   1     SI_name = EV_500_RefSB
   *   2     0   2     SI_name = EV_1KM_RefSB
   *   2     0   3     SI_name = EV_1KM_Emissive
   *   1     0   4     SI_name = EV_250_Aggr500_RefSB
   *   2     0   5     SI_name = EV_250_Aggr1km_RefSB
   *   2     0   6     SI_name = EV_500_Aggr1km_RefSB
   *                        good night indices = [2, 3, 5, 6]
   *   0     1   0     UI_name = EV_250_RefSB_Uncert_Indexes
   *   1     1   1     UI_name = EV_500_RefSB_Uncert_Indexes
   *   2     1   2     UI_name = EV_1KM_RefSB_Uncert_Indexes
   *   2     1   3     UI_name = EV_1KM_Emissive_Uncert_Indexes
   *   1     1   4     UI_name = EV_250_Aggr500_RefSB_Uncert_Indexes
   *   2     1   5     UI_name = EV_250_Aggr1km_RefSB_Uncert_Indexes
   *   2     1   6     UI_name = EV_500_Aggr1km_RefSB_Uncert_Indexes
   *                        good night indices = none
   *   1     2   0     SU_name = EV_250_Aggr500_RefSB_Samples_Used
   *   2     2   1     SU_name = EV_250_Aggr1km_RefSB_Samples_Used
   *   2     2   2     SU_name = EV_500_Aggr1km_RefSB_Samples_Used
   */


  
  for (type = 0; type < Num_Data_types; type++)
  {
    if (type < SU && skip_night_hi_res == True) 
      start_output_res = INDEX_L1B_1km;
    else
      start_output_res = INDEX_L1B_250m;
        
    for (sds_number = start_output_res; 
           sds_number < num_sds_ids[type]; sds_number++)
               
    {
      if (skip_night_hi_res == True)
        if ((type < SU && sds_number == Aggr_250_500_RefSB_Data) || 
            (type == SU && sds_number == Aggr_250_500_RefSB_Samples)) 
          continue;
      
      if (SDsetattr(sds_ids[type][sds_number], "units",
                    DFNT_CHAR, 4, (VOIDP)units) == FAIL)
      {
        sprintf(errmsg, errmsgfmt, "units", (long)type, (long)sds_number);
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
        return returnStatus;
      }

      switch(type)
      {
       case SI:

         if (SDsetattr(sds_ids[type][sds_number], "valid_range",
                       DFNT_UINT16, 2, (VOIDP)SI_valid_range) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "valid_range", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }

         if (SDsetattr(sds_ids[type][sds_number], "_FillValue",
                       DFNT_UINT16, 1,(VOIDP)&SI_fillvalue) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "_FillValue", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }
         break;
       
       case UI:

         if (SDsetattr(sds_ids[type][sds_number], "valid_range",
                       DFNT_UINT8, 2, (VOIDP)UI_valid_range) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "valid_range", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }

         if (SDsetattr(sds_ids[type][sds_number], "_FillValue",
                       DFNT_UINT8, 1,(VOIDP)&UI_fillvalue) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "_FillValue", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }
         break;

       case SU:

         if (SDsetattr(sds_ids[type][sds_number], "valid_range",
                       DFNT_INT8, 2, 
                       (VOIDP)SU_valid_range[sds_number]) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "valid_range", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }
 
         if (SDsetattr(sds_ids[type][sds_number], "_FillValue",
                       DFNT_INT8, 1,(VOIDP)&SU_fillvalue) == FAIL)
         {
           sprintf(errmsg, errmsgfmt, "_FillValue", (long)type, (long)sds_number);
           returnStatus = MODIS_F_WRITE_ERROR;
           L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }
         break;
       default:
         {
           returnStatus = MODIS_F_NOK;
           L1BErrorMsg(location, returnStatus,
                       "Invalid \"type\" ... must be code bug.",
                       "SDsetattr", 0, NULL, True);
           return returnStatus;
         }
      }
    }
  }/* end loop through type */

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
      /*
       * Set the units, valid_range and _FillValue attributes for the
       * Band 26 SDSs.
       */
      /* Scaled Integer SDS */

  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "units",
                DFNT_CHAR, 4, (VOIDP)units) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 SI SDS "
                   "attribute \"units\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }
  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "valid_range",
                DFNT_UINT16, 2, (VOIDP)SI_valid_range) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 SI SDS "
                   "attribute \"valid_range\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }
  if (SDsetattr(L1B_Scan->Band26.SI_sds_id, "_FillValue",
                DFNT_UINT16, 1,(VOIDP)&SI_fillvalue) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 SI SDS "
                   "attribute \"_FillValue\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }

      /* Undertainty SDS */

  if (SDsetattr(L1B_Scan->Band26.UI_sds_id, "units",
                DFNT_CHAR, 4, (VOIDP)units) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 UI SDS "
                   "attribute \"units\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }
  if (SDsetattr(L1B_Scan->Band26.UI_sds_id, "valid_range",
                DFNT_UINT8, 2, (VOIDP)UI_valid_range) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 UI SDS "
                   "attribute \"valid_range\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }
  if (SDsetattr(L1B_Scan->Band26.UI_sds_id, "_FillValue",
                DFNT_UINT8, 1,(VOIDP)&UI_fillvalue) == FAIL)
  {
    strcpy(errmsg, "Could not write Band 26 UI SDS "
                   "attribute \"_FillValue\".");
    returnStatus = MODIS_F_WRITE_ERROR;
    L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, True);
    return returnStatus;
  }

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  return(MODIS_S_OK);
}

PGSt_SMF_status Get_SDS_id(int32      f, 
                           L1B_Scan_t *L1B_Scan,
                           int16      *num_sds, 
                           int32      *sds_id)
/*
!C**********************************************************************
!Description:   For one L1B EV file of index f, this routine gets the number
                of EV, scaled integer (SI), SDSs and an array of the SDS ids. 
                These include both aggregated and non-aggregated SI SDSs.

!Input Parameters:
      f                Index of one of the L1B EV files.
      L1B_Scan         Members which are SDS ids for scaled integer (SI) SDSs.

!Output Parameters:
      num_sds          Number of SDS ids assigned to array sds_id.
      sds_id           Array of SDS ids pertaining to the particular file
                       of index f.  This array should be allocated large
                       enough to hold at least 4 ids.

!Revision History:
 (continue at top of the file)

 Revision 01.01 Feb 10, 1999
 Added address of num_sds to argument list to be able to assign the variable
 in this function.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
  The ordering of the SDS ids follows the ordering in resolution_index_t
  so that we can later use the macros INDEX_250M, ... INDEX_1000M_EMISS.

!END********************************************************************
*/
{
  switch(f)
    {
      case INDEX_L1B_250m:
        *num_sds = 1;
        sds_id[INDEX_250M] = L1B_Scan->SI_sds_id[INDEX_250M];
        break;

      case INDEX_L1B_500m:
        *num_sds = 2;
        sds_id[INDEX_250M]        =
            L1B_Scan->EV_250m_Aggr500m_RefSB_sds_id;
        sds_id[INDEX_500M]        =
            L1B_Scan->SI_sds_id[INDEX_500M];
        break;

      case INDEX_L1B_1km:
        *num_sds =
            4;
        sds_id[INDEX_250M]        =
            L1B_Scan->EV_250m_Aggr1km_RefSB_sds_id;
        sds_id[INDEX_500M]        =
            L1B_Scan->EV_500m_Aggr1km_RefSB_sds_id;
        sds_id[INDEX_1000M_REFL]  =
            L1B_Scan->SI_sds_id[INDEX_1000M_REFL]; 
        sds_id[INDEX_1000M_EMISS] =
            L1B_Scan->SI_sds_id[INDEX_1000M_EMISS];
        break;
      
      default:
        break;
    }
  
   return(MODIS_S_OK);
}

PGSt_SMF_status Set_SDS_Attributes(int32   *sds_id, 
                                   char    **BandNames, 
                                   float32 **scale, 
                                   float32 **offset,
                                   char    *rad_units,
                                   char    *refl_units,
                                   char    *counts_units,
                                   int32   num_sds)
/*
!C**********************************************************************
!Description:   This routine writes the the following attributes of
                scaled integer SDSs: (1) band names, (2) reflectance
                scales, offsets and units (only if reflectance SDSs are
                present), (3) radiance scales, offsets and units and
                (4) corrected counts scales, offsets and units.

!Input Parameters:
      int32 *      sds_id        array of SDS ids for EV SDSs.  The ordering
                                 follows the order in resolution_index_t
                                 so we can use the macros defined there.
      char  **     BandNames     Array of band names for each SDS.
      float32 **   scale
      float32 **   offset
      char    *    rad_units     String containing radiance units.
      char    *    refl_units    String containing reflectance units.
      char    *    counts_units  String containing corrected counts units.
      int32        num_sds       Number of SDSs in array sds_id.
!Output Parameters:
      
!Revision History:
 (continue at top of the file)

 Revision 01.01  October 16, 2004  Razor Issue #200
 Casted Int32 variables in sprintf calls to "long" with the
 format specifier "%ld" for better code portability.
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   If there is an error, write the L1BErrorMsg but do not error-out.
   Let the parent error out so that LUN can be identified.

!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Set_SDS_Attributes";
  char *errmsgfmt = "Could not write attribute \"%s\" for R = %ld";
  char errmsg[256];
  int32 R;
  typedef enum {
    radiance,
    reflectance,
    counts,
    Num_Kinds
  } ScaleOffset_Kind_t;
 
  int16 sds_dim = 4;  
  int32 dim_size[4] = {2, 5, 15, 16};

  for (R = 0; R < num_sds; R++)
  {  
    if (SDsetattr(sds_id[R],"band_names",DFNT_CHAR8,
                  (int32)strlen(BandNames[R]), BandNames[R]) == FAIL) 
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      sprintf(errmsg, errmsgfmt, "band_names", (long)R);
      L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
      return returnStatus;
    }

    if (SDsetattr(sds_id[R],"radiance_scales",DFNT_FLOAT32,
                  dim_size[R],scale[radiance * sds_dim + R]) == FAIL)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      sprintf(errmsg, errmsgfmt, "radiance_scales", (long)R);
      L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
      return returnStatus;
    }

    if (SDsetattr(sds_id[R],"radiance_offsets",DFNT_FLOAT32,
                  dim_size[R], offset[radiance * sds_dim + R]) == FAIL)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      sprintf(errmsg, errmsgfmt, "radiance_offsets", (long)R);
      L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
      return returnStatus;
    }

    if (SDsetattr(sds_id[R],"radiance_units",DFNT_CHAR8,
                  (int32)strlen(rad_units), rad_units) == FAIL)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      sprintf(errmsg, errmsgfmt, "radiance_units", (long)R);
      L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
      return returnStatus;
    }

    if (R != INDEX_1000M_EMISS)
    {
      if (SDsetattr(sds_id[R],"reflectance_scales",DFNT_FLOAT32,
                    dim_size[R],scale[reflectance * sds_dim + R]) 
                    == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "reflectance_scales", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }

      if (SDsetattr(sds_id[R],"reflectance_offsets",DFNT_FLOAT32,
                    dim_size[R],offset[reflectance * sds_dim + R]) 
                    == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "reflectance_offsets", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }

      if (SDsetattr(sds_id[R],"reflectance_units",DFNT_CHAR8,
                    (int32)strlen(refl_units),refl_units) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "reflectance_units", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }

      if (SDsetattr(sds_id[R],"corrected_counts_scales",DFNT_FLOAT32,
                    dim_size[R],scale[counts * sds_dim + R]) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "corrected_counts_scales", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }

      if (SDsetattr(sds_id[R],"corrected_counts_offsets",DFNT_FLOAT32,
                    dim_size[R],offset[counts * sds_dim + R]) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "corrected_counts_offsets", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }

      if (SDsetattr(sds_id[R],"corrected_counts_units",DFNT_CHAR8,
                    (int32)strlen(counts_units),counts_units) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        sprintf(errmsg, errmsgfmt, "corrected_counts_units", (long)R);
        L1BErrorMsg(location, returnStatus, errmsg, 
                       "SDsetattr", 0, NULL, False);
        return returnStatus;
      }
    }

  } /*R*/ 

  return(MODIS_S_OK);
}

PGSt_SMF_status Create_Band_Subsetting_SDS(L1B_granule_t *L1B_Gran,
                                           boolean       skip_night_hi_res)
/*
!C**********************************************************************
!Description:   For each of the L1B output granule files, this routine
                creates band-subsetting SDSs, writes the
                numerical equivalents of the band values and writes
                long-name attributes for each.  Note that this creates new
                SDSs, independent of the Swath data fields (which are
                implemented as Vdata, not as SDSs).
   
!Input Parameters:
    L1B_Gran (->sd_id[])    SD file interface ID for all L1B EV files.
    boolean  skip_night_hi_res   True if and only if all scans are
                              night mode scans and writing of 250m
                              and 500m night mode granules has been
                              turned off
!Output Parameters:
      None

!Revision History:
 (continue at top of the file)

  Revision 01.04  January 17, 2002 Razor Issue #171
  Improve portability of code to 64-bit mode.
  Cast strlen returns to int32 in calls to SDSetattr.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.03 November 16, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic so that
   if skip_night_hi_res is True, output data sets for 250m and 500m resolution
   data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.02 Feb 12, 1999
 Formed this by including the loop through L1B output files in
 Write_Subsetting_SDS and changing the name to an appropriate name.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb 10, 1999
 Made the Band_subsetting_names variable a global at top of file so that the
 same names can be used in create Swath.  Used the L1B_EV_DIM_NAME variable
 for dimension name rather than re-defining it here.  Added SDendaccess.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1998
 Initial development
 Zhenying Gu(zgu@ltpmail.gsfc.nasa.gov)

!Team-unique Header:

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
     file index   res.   # bnd. sub. SDSs     
     ----------  ------  ---------------
         0        250m          1
         1        500m          2
         2        1km           4

!END********************************************************************
*/  
{
  PGSt_SMF_status returnStatus;
  char    *location = "Create_Band_Subsetting_SDS";
  int32   f;
  int32   hdf_return;
  int32   sds_index;
  int32   sds_id;
  int32   num_sds = 0;
  int32   R       = 0;
  int32   dims[NUM_L1A_RESOLUTIONS] = {2, 5, 15, 16};
  float32 data_250m[2] = {1, 2};
  float32 data_500m[5] = {3, 4, 5, 6, 7};
  float32 data_1km_refl[15] = {8, 9, 10, 11, 12, 13, 13.5, 14, 14.5, 15,
                               16, 17, 18, 19, 26};
  float32 data_1km_emiss[16] = {20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31,
                                32, 33, 34, 35, 36};
  float32 *data[NUM_L1A_RESOLUTIONS];  
  char  *SDS_LongName[NUM_L1A_RESOLUTIONS] = 
      {
        "250M Band Numbers for Subsetting",
        "500M Band Numbers for Subsetting",
        "1KM Reflective Solar Band Numbers for Subsetting",
        "1KM Emissive Band Numbers for Subsetting"
      };
  int32 evfile_luns[NUM_L1B_EV_FILES] = {
    L1B_EV_250M_FILE,
    L1B_EV_500M_FILE,
    L1B_EV_1000M_FILE
  };
  char errmsg[256];
  int16 start_output_res = INDEX_L1B_250m;
  
  data[0] = data_250m;
  data[1] = data_500m;
  data[2] = data_1km_refl;
  data[3] = data_1km_emiss;

  /*
   * If we are in night mode and we do not wish to write 250m and 
   * 500m data sets, start with the 1KM output data sets.
   */

  if (skip_night_hi_res == True)
    start_output_res = INDEX_L1B_1km;
  else
    start_output_res = INDEX_L1B_250m;

  for (f = start_output_res; f < NUM_L1B_EV_FILES; f++)
  {

    num_sds = f + 1;
    if (f == INDEX_L1B_1km)
      num_sds++;

    for (R = 0; R < num_sds; R++)
    {

      returnStatus = write_sds_rank1(L1B_Gran->sd_id[f],
                                     Band_subsetting_names[R], 
                                     L1B_EV_DIM_NAME[R][0], dims[R],
                                     "float32", data[R]);
      if(returnStatus != MODIS_S_OK)
      {
        L1BErrorMsg(location, returnStatus, NULL,
                    "write_sds_rank1", evfile_luns[f], NULL, True);
        return returnStatus;
      }
   
      sds_index = SDnametoindex(L1B_Gran->sd_id[f], 
                        Band_subsetting_names[R]);
      if(sds_index == FAIL)
      {
        sprintf(errmsg, "Could not get SDS index of \"%s\".", 
                        Band_subsetting_names[R]);
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDnametoindex", evfile_luns[f], NULL, True);
        return returnStatus;
      }

      sds_id = SDselect(L1B_Gran->sd_id[f], sds_index);
      if(sds_id == FAIL)
      {
        sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
                        Band_subsetting_names[R]);
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDselect", evfile_luns[f], NULL, True);
        return returnStatus;
      }

      if(SDsetattr(sds_id, "long_name", DFNT_CHAR8, 
                   (int32)strlen(SDS_LongName[R]),
                   (VOIDP)SDS_LongName[R]) == FAIL)
      {
        sprintf(errmsg, 
                "Could not write long_name attribute for SDS \"%s\".",
                Band_subsetting_names[R]);
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDsetattr", evfile_luns[f], NULL, True);
        return returnStatus;
      }

      hdf_return = SDendaccess(sds_id);
      if (hdf_return == FAIL)
      {
        sprintf(errmsg, "Could not end access to SDS \"%s\".", 
                        Band_subsetting_names[R]);
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDendaccess", evfile_luns[f], NULL, True);
        return returnStatus;
      }

    }   /* R */
  }   /* f */
  return(MODIS_S_OK);
}

PGSt_SMF_status Copy_Geo_SDS (L1B_granule_t    *L1B_Gran,
                              boolean          skip_night_hi_res)

/*
!C****************************************************************
!Description:
    Read subsampled SDSs from geolocation file and write those data and
    their attributes into the 1km L1B file.  For the 250m and 500m L1B
    granules, read and write the full geolocation latitude and longitude
    SDS (but not others).  Only the units, valid_range and fill-value
    attributes are written to the 250m and 500m granules.  Note that SDSs
    have been previously created (in function Create_L1B_Swath), but open
    access, read/write, and end access, are performed for both input and
    output here.

!Input Parameters:
    L1B_granule_t    *L1B_Gran    
    boolean  skip_night_hi_res   True if and only if all scans are
                              night mode scans and writing of 250m
                              and 500m night mode granules has been
                              turned off

!Output Parameters:
    L1B_granule_t    *L1B_Gran

!Revision History:
 (continue at top of the file)

 Revision 02.14  March 27, 2003  Razor Issue #173 
 Iinitialized both "attr_name" and "attr_value" to NULL for ANSI-C compliance.
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 02.13 November 16, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic 
   so that if skip_night_hi_res is True, output data sets for 250m and 
   500m resolution data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.12, Feb 22, 1999
 Moved calculations of Nadir_Frame ... to Scan_Meta_Cal, deleted 
 L1B_ScanMeta from argument list since no longer needed.
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.11, Feb 19, 1999
 Added SDendaccess for the L1B sds_ids in the loop through 250m and
 500m files when writing the lat and long SDSs.
 Also, for the geo_sds_id in the 1st loop through isds, the SDendaccess 
 was moved to be just before the end-loop curly brace (after the loop
 through L1B granules).
 Jim Rogers (rogers@msct.gsfc.nasa.gov)

 Revision 02.10 Oct. 1997
 Modified to copy SDS attributes instead of write.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 02.00 April 1997
 Modified to match some changes in data structures.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.00 1997/02/23  
 Initial development
 Neal Devine(neal.devine@gsfc.nasa.gov)

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
  PGSt_SMF_status  returnStatus   = MODIS_S_OK;
  char             *location      = "Copy_Geo_SDS";
  PGSt_integer     Version        = 1;
  char             fname[512];
  int32            geo_sd_id      = 0;
  geo_sds_index_t  isds           = INDEX_LATITUDE; /*Geo sds index*/
  int16            i              = 0; /*generic index*/
  int32            geo_sds_index  = 0; /*Geo sds index in hdf file*/
  int32            geo_sds_id     = 0;
  int32            geo_attr_index = 0;
  int32            L1B_sds_index  = 0;
  int32            L1B_sds_id     = 0;
  int32            geo_start[2]   = {GEO_OFFSET, GEO_OFFSET};
  int32            geo_stride[2]  = {GEO_STRIDE, GEO_STRIDE};
  int32            start[2]       = {0, 0};
  int32            edge[2]        = {0, 0};
  intn             hdf_return     = FAIL;
  int32            evfile_luns[NUM_L1B_EV_FILES] 
                                  = {
                                     L1B_EV_250M_FILE,
                                     L1B_EV_500M_FILE,
                                     L1B_EV_1000M_FILE
                                    };

       /* these two are used for setting gflags attributes */
  char            *attr_name = NULL; 
  char            *attr_value = NULL; 

#define MAX_ATTR_NAME_SIZE   26
#define MAX_ATTR_BUFFER_SIZE 72

         /* used in SDattrinfo() */
  char            attr_name_buffer[MAX_ATTR_NAME_SIZE];  

  int32           attr_data_type  = 0;
  int32           attr_count      = 0;
  char            attr_buffer[MAX_ATTR_BUFFER_SIZE];
  float32         *buffer;
  buffer = malloc(10 * MAX_NUM_SCANS * EV_1km_FRAMES * sizeof(float32));
  char            errmsg[256];


#define COPY_ATTR(attr_name)                                             \
  geo_attr_index = SDfindattr(geo_sds_id, attr_name);                    \
  if (geo_attr_index == FAIL)                                            \
  {                                                                      \
    sprintf(errmsg, "Could not get attribute index for \"%s\".",         \
                    attr_name);                                          \
    returnStatus = MODIS_F_READ_ERROR;                                   \
    L1BErrorMsg(location, returnStatus, errmsg,                          \
                "SDfindattr", GEOLOCATION_FILE, NULL, True);             \
    return returnStatus;                                                 \
  }                                                                      \
  hdf_return = SDattrinfo(geo_sds_id, geo_attr_index, attr_name_buffer,  \
                          &attr_data_type, &attr_count);                 \
  if (hdf_return == FAIL)                                                \
  {                                                                      \
    sprintf(errmsg, "Could not get attribute info for \"%s\".",          \
            attr_name);                                                  \
    returnStatus = MODIS_F_HDF_ERROR;                                    \
    L1BErrorMsg(location, returnStatus, errmsg,                          \
                "SDattrinfo", GEOLOCATION_FILE, NULL, True);             \
    return returnStatus;                                                 \
  }                                                                      \
  hdf_return = SDreadattr(geo_sds_id, geo_attr_index,                    \
                         (void *)attr_buffer);                           \
  if (hdf_return == FAIL)                                                \
  {                                                                      \
    sprintf(errmsg, "Could not read attribute \"%s\".", attr_name);      \
    returnStatus = MODIS_F_READ_ERROR;                                   \
    L1BErrorMsg(location, returnStatus, errmsg,                          \
                "SDreadattr", GEOLOCATION_FILE, NULL, True);             \
    return returnStatus;                                                 \
  }                                                                      \
  hdf_return = SDsetattr(L1B_sds_id, attr_name, attr_data_type,          \
                         attr_count, (void *)attr_buffer);               \
  if (hdf_return == FAIL)                                                \
  {                                                                      \
    sprintf(errmsg, "Could not write attribute \"%s\".", attr_name);     \
    returnStatus = MODIS_F_WRITE_ERROR;                                  \
    L1BErrorMsg(location, returnStatus, errmsg,                          \
                "SDsetattr", evfile_luns[i], NULL, True);                \
    return returnStatus;                                                 \
  }

  /* Get Geoloc file name */

  Version = 1;
  if ( PGS_PC_GetReference(GEOLOCATION_FILE, &Version, fname) != 
                           PGS_S_SUCCESS )
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, 
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  /* Open geolocation file */

  geo_sd_id = SDstart(fname, DFACC_RDONLY);
  if (geo_sd_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, 
                "Could not open file for SD read access.",
                "SDstart", GEOLOCATION_FILE,
                "The file may be missing, corrupted or "
                "not an HDF-4 file.", True);
    return returnStatus;
  }

  /* loop over all Geo_SDSs.  This first loop is for L1B EV 250m 
   * and EV 500m granules only, so if skip_night_hi_res is True, will not
   * be executed.
   */

  if (skip_night_hi_res == False)
  {
    edge[0] = 10 * L1B_Gran->num_scans;
    edge[1] = EV_1km_FRAMES ;

    for (isds = INDEX_LATITUDE; isds < INDEX_HEIGHT; isds++)
    {
      geo_sds_index = SDnametoindex(geo_sd_id, GEO_SDS[isds].src_name);
      if (geo_sds_index == FAIL)
      {
        sprintf(errmsg, "Could not get SDS index of \"%s\".", 
            GEO_SDS[isds].src_name);
        returnStatus = MODIS_F_READ_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDnametoindex", GEOLOCATION_FILE, NULL, True);
        return returnStatus;
      }

      geo_sds_id = SDselect(geo_sd_id, geo_sds_index);
      if (geo_sds_id == FAIL)
      {
        sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
            GEO_SDS[isds].src_name);
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDselect", GEOLOCATION_FILE, NULL, True);
        return returnStatus;
      }

      hdf_return = SDreaddata(geo_sds_id, start, NULL, edge, 
                              (void *)buffer);
      if (hdf_return == FAIL)
      {
        sprintf(errmsg, "Could not read data for SDS \"%s\".", 
            GEO_SDS[isds].src_name);
        returnStatus = MODIS_F_READ_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDreaddata", GEOLOCATION_FILE, NULL, True);
        return returnStatus;
      }

      for (i = INDEX_L1B_250m; i < INDEX_L1B_1km; i++)
      {
        L1B_sds_index = SDnametoindex(L1B_Gran->sd_id[i], 
                                      GEO_SDS[isds].name);
        if (L1B_sds_index == FAIL)
        {
          sprintf(errmsg, "Could not get SDS index of \"%s\".", 
                  GEO_SDS[isds].name);
          returnStatus = MODIS_F_HDF_ERROR;
          L1BErrorMsg(location, returnStatus, errmsg,
                      "SDnametoindex", evfile_luns[i], NULL, True);
          return returnStatus;
        }

        L1B_sds_id = SDselect(L1B_Gran->sd_id[i], L1B_sds_index);
        if (L1B_sds_id == FAIL)
        {
          sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
                  GEO_SDS[isds].name);
          returnStatus = MODIS_F_HDF_ERROR;
          L1BErrorMsg(location, returnStatus, errmsg,
                      "SDselect", evfile_luns[i], NULL, True);
          return returnStatus;
        }

        hdf_return = SDwritedata(L1B_sds_id, start, 
                                 NULL, edge, (void *)buffer);
        if (hdf_return == FAIL)
        {
          sprintf(errmsg, "Could not write SDS \"%s\".", 
                          GEO_SDS[isds].name);
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, errmsg,
                      "SDwritedata", evfile_luns[i], NULL, True);
          return returnStatus;
        }

        COPY_ATTR("units");
        COPY_ATTR("valid_range");
        COPY_ATTR("_FillValue"); 

        /* Done with this L1B_sds_id
         */
        hdf_return = SDendaccess(L1B_sds_id);
        if (hdf_return == FAIL)
        {
          sprintf(errmsg, "Could not end access for SDS \"%s\".", 
              GEO_SDS[isds].name);
          returnStatus = MODIS_F_HDF_ERROR;
          L1BErrorMsg(location, returnStatus, errmsg,
                      "SDendaccess", evfile_luns[i], NULL, True);
          return returnStatus;
        }

      }

      /* End access for this geolocation SDS.
       */
      hdf_return = SDendaccess(geo_sds_id);
      if (hdf_return == FAIL)
      {
        sprintf(errmsg, "Could not end access for SDS \"%s\".", 
            GEO_SDS[isds].src_name);
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDendaccess", GEOLOCATION_FILE, NULL, True);
        return returnStatus;
      }

    }    /* for (isds = INDEX_LATITUDE; isds < INDEX_HEIGHT; isds++) */
  }      /* if (skip_night_hi_res == False) */

  /* Loop through all geolocations SDSs.  This loop is for the 
   * L1B EV 1km granule.
   */
 
  edge[0] = 2 * L1B_Gran->num_scans;
  edge[1] = 1 + (EV_1km_FRAMES - 3) / 5;

  for (isds = INDEX_LATITUDE; isds < NUM_GEO_SDS; isds++)
  {
    /* Open a geo sds */

    geo_sds_index = SDnametoindex(geo_sd_id, GEO_SDS[isds].src_name);
    if (geo_sds_index == FAIL)
    {
      sprintf(errmsg, "Could not get SDS index of \"%s\".", 
          GEO_SDS[isds].src_name);
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDnametoindex", GEOLOCATION_FILE, NULL, True);
      return returnStatus;
    }

    geo_sds_id = SDselect(geo_sd_id, geo_sds_index);
    if (geo_sds_id == FAIL)
    {
      sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
          GEO_SDS[isds].src_name);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDselect", GEOLOCATION_FILE, NULL, True);
      return returnStatus;
    }

    /* Open an L1B sds */
    
    L1B_sds_index = SDnametoindex(L1B_Gran->sd_id[INDEX_L1B_1km], 
        GEO_SDS[isds].name);
    if (L1B_sds_index == FAIL)
    {
      sprintf(errmsg, "Could not get SDS index of \"%s\".", 
          GEO_SDS[isds].name);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDnametoindex", evfile_luns[INDEX_L1B_1km], NULL, True);
      return returnStatus;
    }

    L1B_sds_id = SDselect(L1B_Gran->sd_id[INDEX_L1B_1km], L1B_sds_index);
    if (L1B_sds_id == FAIL)
    {
      sprintf(errmsg, "Could not get SDS ID for \"%s\".", 
          GEO_SDS[isds].name);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDselect", evfile_luns[INDEX_L1B_1km], NULL, True);
      return returnStatus;
    }

    /* Read an sds from geo file into buffer */
    
    hdf_return = SDreaddata(geo_sds_id, geo_start, geo_stride, 
        edge, (void *)buffer);
    if (hdf_return == FAIL)
    {
      sprintf(errmsg, "Could not read data for SDS \"%s\".", 
          GEO_SDS[isds].src_name);
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDreaddata", GEOLOCATION_FILE, NULL, True);
      return returnStatus;
    }

    /* Write the sds stored in the buffer to L1B_1km file */
    
    hdf_return = SDwritedata(L1B_sds_id, start, NULL, edge, (void *)buffer);
    if (hdf_return == FAIL)
    {
      sprintf(errmsg, "Could not write SDS \"%s\".", GEO_SDS[isds].name);
      returnStatus = MODIS_F_WRITE_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDwritedata", evfile_luns[INDEX_L1B_1km], NULL, True);
      return returnStatus;
    }

  /*
   * Copy attributes: 
   *   copy units, valid_range, _FillValue for all sdss except gflags, 
   *   and scale_factor for all except Latitude and Longitude and gflags
   * Set other attributes:
   *   set line_numbers and frame_numbers for all sdss except gflags
   */

    if (isds != INDEX_GFLAGS)
    {
      i = INDEX_L1B_1km;             /* to set correct lun in the macro */
      COPY_ATTR("units");
      COPY_ATTR("valid_range");
      COPY_ATTR("_FillValue");

      if (isds > INDEX_HEIGHT)
      {
        COPY_ATTR("scale_factor");
      }

      hdf_return = SDsetattr(L1B_sds_id, "line_numbers", DFNT_CHAR8,
                             (int32)strlen("3,8"), (VOIDP)"3,8");
      if (hdf_return == FAIL)
      {
        sprintf(errmsg, "Could not write attribute \"line_numbers\" "
                        "to SDS \"%s\".",
                        GEO_SDS[isds].name);
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDsetattr", evfile_luns[INDEX_L1B_1km], NULL, True);
        return returnStatus;
      }

      hdf_return = SDsetattr(L1B_sds_id, "frame_numbers", DFNT_CHAR8,
                             (int32)strlen("3,8,13,..."), 
                             (VOIDP)"3,8,13,...");
      if (hdf_return == FAIL)
      {
        sprintf(errmsg, "Could not write attribute \"frame_numbers\" "
                        "to SDS \"%s\".",
                GEO_SDS[isds].name);
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "SDsetattr", evfile_luns[INDEX_L1B_1km], 
                    NULL, True);
        return returnStatus;
      }
    }
    else
    {
      i = INDEX_L1B_1km;             /* to set correct lun in the macro */
      COPY_ATTR("_FillValue");
      for (i = 0; i < 5; i++)
      {
        switch(i)
        {
          case 0: attr_name = "Bit 7(MSB)"; 
                  attr_value = "1 = invalid input data";        
                  break;
          case 1: attr_name = "Bit 6"; 
                  attr_value = "1 = no ellipsoid intersection"; 
                  break;
          case 2: attr_name = "Bit 5"; 
                  attr_value = "1 = no valid terrain data";     
                  break;
          case 3: attr_name = "Bit 4"; 
                  attr_value = "1 = DEM missing or of inferior quality"; 
                  break;
          case 4: attr_name = "Bit 3"; 
                  attr_value = "1 = invalid sensor range";
        }
        hdf_return = SDsetattr(L1B_sds_id, attr_name, DFNT_CHAR8,
                               (int32)strlen(attr_value), attr_value);
        if (hdf_return == FAIL)
        {
          sprintf(errmsg,  
                  "Could not write attribute \"%s\" to SDS \"%s\".",
                  attr_name, GEO_SDS[isds].name);
          returnStatus = MODIS_F_WRITE_ERROR;
          L1BErrorMsg(location, returnStatus, errmsg,
                      "SDsetattr", evfile_luns[INDEX_L1B_1km], NULL, True);
          return returnStatus;
        }
      }
    } 
      
    /* Done with this geo_sds */
    
    hdf_return = SDendaccess(geo_sds_id);
    if (hdf_return == FAIL)
    {
      sprintf(errmsg, "Could not end access for SDS \"%s\".", 
                      GEO_SDS[isds].src_name);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDendaccess", GEOLOCATION_FILE, NULL, True);
      return returnStatus;
    }

    /* Done with this L1B_sds */
    
    hdf_return = SDendaccess(L1B_sds_id);
    if (hdf_return == FAIL)
    {
      sprintf(errmsg, "Could not end access for SDS \"%s\".", 
              GEO_SDS[isds].name);
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, errmsg,
                  "SDendaccess", evfile_luns[INDEX_L1B_1km], NULL, True);
      return returnStatus;
    }

  }  /* for (isds = INDEX_LATITUDE; isds < NUM_GEO_SDS; isds++) */

  /* Done with the geo file */
  
  hdf_return = SDend(geo_sd_id);
  if (hdf_return == FAIL)
  {
     returnStatus = MODIS_F_HDF_ERROR;
     L1BErrorMsg(location, returnStatus, NULL, "SDend", 
                  GEOLOCATION_FILE,
                  "Memory or the disk file must have become corrupted.", 
                  True);
     return returnStatus;
  }

  free(buffer);
  
  return(MODIS_S_OK);
}

PGSt_SMF_status Scan_Meta_Cal(lookup_tables_t     *tables,
                              L1A_granule_t       *L1A_Gran,
                              L1B_granule_t       *L1B_Gran,
                              L1B_Scan_Metadata_t *L1B_Scan_Meta,
                              QA_Data_t           *QA)
/*
!C****************************************************************
!Description:
    This routine assigns members of the L1B_Scan_Meta to be written to
    the L1B EV files later in Write_L1B_ScanMeta.  Values already
    determined are simply assigned while other data, such as telemetry
    and SRCA calibration mode, are read in from the L1A granule and
    then assigned.

!Input Parameters:
    L1A_granule_t       *L1A_Gran
    L1B_granule_t       *L1B_Gran
    QA_Data_t           *QA
!Output Parameters:
    L1B_Scan_Metadata_t *L1B_Scan_Meta
    QA_Data_t           *QA
!Revision History:
 (continue at top of the file)

 Revision 01.06  October 24, 2003   Razor Issue #196 (formerly Issue #184)
 All changes are in the function "Scan_Meta_Cal". Added parameter "tables",
 array "Attitude_Angles", and constant "rtod". Moved the code which opens the
 geolocation file to just before the Scan-by-Scan loop which sets bit flags. 
 Added the code which reads "Attitude_Angles" from the geolocation file. 
 Changed the setting of QA flag bit 1 (Spacecraft Maneuver) so that the code 
 will check the attitude angles against their threshold angles in the QA LUT. 
 If any one of them is more than its threshold angle set the "maneuver" bit 
 for that scan. Added a new "Sci Abnormal" bit (bit 26) which is set when 
 SS_FR_SCIABNORM was set (the way the "maneuver" flag is originaly set).
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 01.05  April 23, 2003
  Changed the telemetry mnemonic which determines whether the nadir aperture
  door (NAD) is open from "CR_DR_NAD_OPEN" to the  equivalent "CR_DR_NAD_CLSD"
  mnemonic and altered the L1B code logic to correctly use the new value.  It
  was discovered that at two different times in the history of the MODIS/Terra
  instrument the mnemonic "CR_DR_NAD_OPEN" did not correctly reflect the state
  of the NAD, whereas the mnemonic "CR_DR_NAD_CLSD" has been consistently
  reliable.  
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)
 
 Revision 01.04 March 27, 2003  Razor Issue #173
 Enclosed the rows in the initializer of array-of-struct "temp" with braces 
 for ANSI-C compliance. 
 Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

 Revision 01.03 April 16, 2002  Razor Issue #166
 Added leading/trailing granule scan gap to bit QA flags.
 Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.02 August 26, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb 20, 1999
 Fixed some minor bugs.  Moved the assignment of nadir frame values into
 macro NADIR_FRAME_GEO_VALS.  Moved the assignment of nadir frame
 latitude and longitude into here from Create_L1B_Swath to have all
 these in one place.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1998/4/8
 Original Development
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)

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
  PGSt_SMF_status   returnStatus;
  PGSt_integer      Version       = 1;
  int32             geo_sd_id;
  char              fname[512];
  int               f             = 0;
  int               S             = 0;
  int               start_scan    = 0;
  int16             *SRCA_mode;
  SRCA_mode = (int16*) malloc(MAX_NUM_SCANS * sizeof(int16));
  int16             *i16_buffer;
  i16_buffer = (int16*) malloc(10 * MAX_NUM_SCANS * EV_1km_FRAMES * sizeof(int16));
  float32           *f32_buffer;
  f32_buffer = (float32*) malloc(10 * MAX_NUM_SCANS * EV_1km_FRAMES * sizeof(float32));
  int32             buf_index     = 0;
  float32           *SD_Sun_zenith;
  SD_Sun_zenith = (float32*) malloc(MAX_NUM_SCANS * sizeof(float32));
  float32           *SD_Sun_Azimuth;
  SD_Sun_Azimuth = (float32*) malloc(MAX_NUM_SCANS * sizeof(float32));
  float32           Z_MODIS;
  char              *location     = "Scan_Meta_Cal";
  char              upper_bound_str[20];
  typedef float64 Attitude_Angles_t[NUM_ATTITUDE_ANGLES];
  Attitude_Angles_t *Attitude_Angles;
  Attitude_Angles = (Attitude_Angles_t*) malloc(MAX_NUM_SCANS * sizeof(Attitude_Angles_t));
  float32           rtod          = 180.0/PGS_PI;


  enum {
    SCAN_SECTOR_ROTATION,
    PC_LW_A,
    PC_LW_B,
    PV_VIS_A,
    PV_VIS_B,
    PV_NIR_A,
    PV_NIR_B,
    PV_SW_A,
    PV_SW_B,
    PV_LW_A,
    PV_LW_B,
    SDSM_A,
    SDSM_B,
    RC_CSH,
    RC_ISH,
    RC_OSH,
    BB_HEATER_A,
    BB_HEATER_B,
    SD_DOOR,
    SD_SCREEN,
    SCI_ABNORM,
    NAD_DOOR,
    MACRO_ID,
    PCLW_ADC_PRI,
    PCLW_ADC_RED,
    NUM_FIELD
  }i;
  struct {
    char *vname;     /*vdata_name*/
    char *fname;     /*field_name*/
       /*data buffer; all eng_vdata are type of uint16*/
    uint16 buffer[MAX_NUM_SCANS];  
  } temp[NUM_FIELD] = {
    /* XG 3/12/2015 
     * use SET_FR_ENC_DELTA in the field Command Parameters to get the 
     * sector rotation status
     */
    /*{"Telemetry Major Cycle 6 of 7",      "CS_FR_ENC_DELTA",  {0} },*/
    {"Command Parameters",                "SET_FR_ENC_DELTA", {0} },
    {"Telemetry Major Cycle 4A of 7",     "CR_PCLWA_ECAL_ON", {0} },
    {"Telemetry Major Cycle 4A of 7",     "CR_PCLWB_ECAL_ON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVVISA_ECALON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVVISB_ECALON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVNIRA_ECALON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVNIRB_ECALON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVSMA_ECAL_ON", {0} },
    {"Telemetry Major Cycle 4B of 7",     "CR_PVSMB_ECAL_ON", {0} },
    {"Telemetry Major Cycle 4A of 7",     "CR_PVLWA_ECAL_ON", {0} },
    {"Telemetry Major Cycle 4A of 7",     "CR_PVLWB_ECAL_ON", {0} },
    {"Telemetry Major Cycle 5B of 7",     "CR_SM_SDSM_A_ON",  {0} },
    {"Telemetry Major Cycle 5B of 7",     "CR_SM_SDSM_B_ON",  {0} },
    {"Telemetry Major Cycle 5A of 7",     "CR_RC_CSHTR_ON",   {0} },
    {"Telemetry Major Cycle 5A of 7",     "CR_RC_ISHTR_ON",   {0} },
    {"Telemetry Major Cycle 5A of 7",     "CR_RC_OSHTR_ON",   {0} },
    {"Telemetry Major Cycle 0 of 7",      "CR_BB_A_PWR_ON",   {0} },
    {"Telemetry Major Cycle 0 of 7",      "CR_BB_B_PWR_ON",   {0} },
    {"Telemetry Major Cycle 3A of 7",     "CR_DR_SDD_OPEN",   {0} },
    {"Telemetry Major Cycle 3A of 7",     "CR_DR_SDS_OPEN",   {0} },
    {"Telemetry Major Cycle 1 of 7",      "SS_FR_SCIABNORM",  {0} },
    {"Telemetry Major Cycle 3A of 7",     "CR_DR_NAD_CLSD",   {0} },
    {"Telemetry Major Cycle All Part 3",  "SS_CP_MACRO_ID",   {0} },
    {"Telemetry Major Cycle 4A of 7" ,    "CR_PCLW_A_ON",     {0} }, 
    {"Telemetry Major Cycle 4A of 7" ,    "CR_PCLW_B_ON",     {0} } 
  };

/********************** MACRO NADIR_FRAME_GEO_VALS ************************
  For one geolocation SDS, this macro will read the data from the
  geolocation SDS (includes all scan lines and all frames) and assign
  one value for the nadir frame per scan to an appropriate variable
  in L1B_Scan_Meta.  "geoname" is a string that is the name of the
  geolocation SDS and "var" is the variable to assign values to. 
  "scalefactor" converts numbers in the SDS to match the units of the
  data.  The geolocation file should already be open and the variable
  geo_sd_id should thus be valid.  "buf" is the buffer (either int16
  or float 32) that we read data into.  Lat and Long use float32, the
  others use int16.
****************************************************************************/
#define NADIR_FRAME_GEO_VALS(geoname,var,scalefactor,buf)               \
  returnStatus = read_sds_rank2(geo_sd_id, geoname,                     \
                                10 * L1A_Gran->num_scans,               \
                                EV_1km_FRAMES, buf);                    \
  if (returnStatus != MODIS_S_OK)                                       \
  {                                                                     \
    L1BErrorMsg(location, returnStatus, NULL,                           \
                "read_sds_rank2", GEOLOCATION_FILE, NULL, True);        \
    return returnStatus;                                                \
  }                                                                     \
  for (S = 0; S < L1A_Gran->num_scans; S++) {                           \
    buf_index = (S * 10 + 4) * EV_1km_FRAMES + NADIR_1km_FRAME_NUM - 1; \
    var[S] = buf[buf_index] * scalefactor;                              \
  }
 
  for (i = SCAN_SECTOR_ROTATION; i < NUM_FIELD; i++)
  {
    returnStatus = read_vdata (L1A_Gran->v_id, 
                               start_scan,
                               L1A_Gran->num_scans,
                               temp[i].vname, 
                               temp[i].fname,
                               (VOIDP)temp[i].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "read_vdata", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    if (i == MACRO_ID)
      strcpy(upper_bound_str, "31");
    else if (i == SCAN_SECTOR_ROTATION)
      strcpy(upper_bound_str, "16383");
    else
      strcpy(upper_bound_str, "1");

    returnStatus = Check_Valid_Range(temp[i].fname, DFNT_UINT16,
                                     NULL, upper_bound_str, NULL,
                                     L1A_Gran->num_scans, 
                                     (void *) temp[i].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
  }


  /*
   * Set the SRCA calibration mode based on macro ID.  This will be a
   * more accurate determination than using the SDS "SRCA calibration
   * mode". If the mode is "undetermined", use the electronics
   * redundancy vector to determine if the SRCA is actually on or
   * off.  If ON, then either the macro has not turned on yet or this
   * is a special SRCA event (not one of the normal OAs).
   */


  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    if (temp[MACRO_ID].buffer[S] == 15 ||
        temp[MACRO_ID].buffer[S] == 16 ||
        temp[MACRO_ID].buffer[S] == 17)
      SRCA_mode[S] = 0;                         /* Radiometric */
    else if (temp[MACRO_ID].buffer[S] == 18 ||
             temp[MACRO_ID].buffer[S] == 19 ||
             temp[MACRO_ID].buffer[S] == 20 ||
             temp[MACRO_ID].buffer[S] == 21)
      SRCA_mode[S] = 2;                         /* Spectral */
    else if (temp[MACRO_ID].buffer[S] == 22 ||
             temp[MACRO_ID].buffer[S] == 23 ||
             temp[MACRO_ID].buffer[S] == 24)
      SRCA_mode[S] = 1;                         /* Spatial */
    else
      SRCA_mode[S] = 3;                         /* undetermined */
  }

  /* calculate Nadir frame quantities and assign to L1B_Scan_Meta */

  if (PGS_PC_GetReference(GEOLOCATION_FILE, &Version, fname) != 
                                                  PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, 
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  /* Open geolocation file
   */
  geo_sd_id = SDstart(fname, DFACC_RDONLY);
  if (geo_sd_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, 
                "Could not open file for SD read access.",
                "SDstart", GEOLOCATION_FILE,
                "The file may be missing, corrupted or "
                "not an HDF-4 file.", True);
    return returnStatus;
  }

  /* Read s/c attitude angles */
  returnStatus = read_sds_rank2(geo_sd_id, "attitude_angles",
                                L1A_Gran->num_scans,
                                NUM_ATTITUDE_ANGLES, Attitude_Angles);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank2",
                GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
  /* initialize all bits to zero */

    L1B_Scan_Meta->Bit_QA_Flags[S] = 0;

  /* bit 0: Moon within defined limits of SVP */

    if (QA->QA_refl.moon_in_SV_KOB_RSB[S] == MOON_INSIDE_SV_KOB ||
        QA->QA_emiss.moon_in_SV_KOB_TEB[S] == MOON_INSIDE_SV_KOB)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000001;

  /* bit 1: Spacecraft Maneuver */
  
    if ( fabs((double) (Attitude_Angles[S][0]))*rtod >
                tables->QA.common_QA_tables.roll_threshold_angle ||
         fabs((double) (Attitude_Angles[S][1]))*rtod >
                tables->QA.common_QA_tables.pitch_threshold_angle ||
         fabs((double) (Attitude_Angles[S][2]))*rtod >
                tables->QA.common_QA_tables.yaw_threshold_angle )
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000002;

  /* bit 2: Sector Rotation */

    if ( temp[SCAN_SECTOR_ROTATION].buffer[S])
    {
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000004;
       QA->QA_common.Sector_Rotation[S] = True;
    }
    else
       QA->QA_common.Sector_Rotation[S] = False;

  /* Set QA flag for both PCLW A and B on */ 
    if ( temp[PCLW_ADC_PRI].buffer[S] == 1 && temp[PCLW_ADC_RED].buffer[S] == 1)
    {
       QA->QA_common.Electronic_Anomaly[S] = True;
    }
    else
       QA->QA_common.Electronic_Anomaly[S] = False;


  /* bit 3: Negative Radiance Beyond Noise Level, is set in Emissive_Cal. */

  /* bit 4: PC Ecal on */

    if (temp[PC_LW_A].buffer[S] || temp[PC_LW_B].buffer[S])
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000010;

  /* bit 5: PV Ecal on */

    if (temp[PV_VIS_A].buffer[S]    || temp[PV_VIS_B].buffer[S]
        || temp[PV_NIR_A].buffer[S] || temp[PV_NIR_B].buffer[S]
        || temp[PV_SW_A].buffer[S]  || temp[PV_SW_B].buffer[S]
        || temp[PV_LW_A].buffer[S]  || temp[PV_LW_B].buffer[S])
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000020;

  /* bit 6: SD Door Open */

    if (temp[SD_DOOR].buffer[S] == 1)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000040;

  /* bit 7: SD Screen Down */

    if (temp[SD_SCREEN].buffer[S] == 0)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000080;

  /*
   * check NAD door open and this will be used in Reflective_Cal
   * and Emissive_Cal 
   */

    if (temp[NAD_DOOR].buffer[S] == 0)
       QA->QA_common.NAD_Door_Open[S] = True;
    else
    {
       QA->QA_common.NAD_Door_Open[S] = False;

   /* bit 8: NAD closed. Added by CCR-508. */

       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000100;
    }

  /* bit 9: SDSM on */

    if (temp[SDSM_A].buffer[S] || temp[SDSM_B].buffer[S])
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000200;

  /* bit 10: Radcooler Heaters On */

    if (temp[RC_CSH].buffer[S] || temp[RC_ISH].buffer[S]
        || temp[RC_OSH].buffer[S])
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000400;

  /*
   * bit 11: Day mode bands telemetered at night. Added by CCR-508.
   * It is set at the end of this function.
   */

  /*
   * bit 12: Linear Emissive Calibration.
   * The current emissive algorithm is non-linear. So it is not set.
   */
 
  /* bit 13: DC Restore Change. Set in Calculate_DCR_Change() */

  /* bit 14: Unused. */

  /* bit 15: BB Heater on */
    if (temp[BB_HEATER_A].buffer[S] || temp[BB_HEATER_B].buffer[S])
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00008000;

  /* bit 16: Missing Leading Granule */

    if (QA->QA_common.missing_leading_granule == True)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00010000;

  /* bit 17: Missing Trailing Granule */

    if (QA->QA_common.missing_trailing_granule == True)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00020000;

  /* bit 18, 19: SRCA calibration mode */

    if (SRCA_mode[S] == 2 || SRCA_mode[S] == 3 || SRCA_mode[S] == -1)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00040000;
    if (SRCA_mode[S] == 1 || SRCA_mode[S] == 3 || SRCA_mode[S] == -1)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00080000;

  /* bit 20: Moon within the SV keep-out box for RSB (any band)*/

    if (QA->QA_refl.moon_in_SV_KOB_RSB[S] == MOON_INSIDE_SV_KOB)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00100000;

  /* bit 21: Moon within the SV keep-out box for Emissive bands (any band)*/

    if (QA->QA_emiss.moon_in_SV_KOB_TEB[S] == MOON_INSIDE_SV_KOB)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00200000;
  
  /* bit 22: all SV data are bad for any subsamples, detectors and bands */

    if (QA->QA_refl.all_SV_DN_bad[S] == 1)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00400000;

  /* bit 23: all BB data are bad for any subsamples, detectors and bands */

    if (QA->QA_refl.all_BB_DN_bad[S] == 1)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00800000;

  /* bit 24: Leading Granule Scan Gap */

    if (QA->QA_common.leading_granule_scan_gap == True)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x01000000;

  /* bit 25: Trailing Granule Scan Gap */

    if (QA->QA_common.trailing_granule_scan_gap == True)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x02000000;
    
  /* bit 26: Sci Abnormal*/

    if (temp[SCI_ABNORM].buffer[S] == 0)
       L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x04000000;

  /* bits 27 - 31: Remaining bits reserved for future use */

  }

  /********************************
   *  Initialize variable values  *
   ********************************/
 
  L1B_Scan_Meta->num_scans = L1A_Gran->num_scans;
  L1B_Scan_Meta->EV_frames = EV_1km_FRAMES;
  L1B_Scan_Meta->Nadir_Frame_Number = NADIR_1km_FRAME_NUM;

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    strcpy(L1B_Scan_Meta->ScanType[S], L1A_Gran->ScanType[S]);
    L1B_Scan_Meta->MirrorSide[S]            = 
          L1A_Gran->MirrorSide[S];
    L1B_Scan_Meta->EVStartTime_TAIsecond[S] =
          L1A_Gran->EVStartTime_TAIsecond[S];
  }

  for (f = 0; f< NUM_L1B_EV_FILES; f++)
    L1B_Scan_Meta->v_id[f] = L1B_Gran->v_id[f];

  returnStatus = read_sds_rank1(geo_sd_id, "SD Sun azimuth",
                                L1A_Gran->num_scans, SD_Sun_Azimuth);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank1",
                GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  returnStatus = read_sds_rank1(geo_sd_id, "SD Sun zenith",
                                L1A_Gran->num_scans, SD_Sun_zenith);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_sds_rank1",
                GEOLOCATION_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Validation of the SWIR out-of-band correction algorithm requires
   * the collection of reflective band data at night. In this case, the
   * scan type is "Day", even it is really night. The following calculation
   * determines the MODIS nadir pixel sees day or night. If it is night and
   * the scan type is "Day", set bit 11 to 1. See MCST internal memo "How
   * To Determine Whether the MODIS Nadir Pixel Seeds Day Or Night" by
   * Joe Esposito and Bruce Berriman, Nov 1, 1999 (Draft) and CCR-508.
   */

  /*
   * The SD sun zenith and azimuth angles in the MOD03 file are in radians.
   */

  for (S = 0; S < L1A_Gran->num_scans; S++)
  {
    Z_MODIS = - 0.34582 * sin((double)SD_Sun_zenith[S]) *
                          cos((double)SD_Sun_Azimuth[S]) +
                0.923830 * cos((double)SD_Sun_zenith[S]);

    /* bit 11: Day mode bands telemetered at night */

    if (Z_MODIS > 0 && strcmp(L1A_Gran->ScanType[S], "Day") == 0)
      L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00000800;
  }

  NADIR_FRAME_GEO_VALS("SolarAzimuth",
                       L1B_Scan_Meta->Nadir_Frame_Solar_Azimuth,
                       SOLAR_AZIMUTH_ZENITH_SCALE_FACTOR,i16_buffer);

  NADIR_FRAME_GEO_VALS("SolarZenith",
                       L1B_Scan_Meta->Nadir_Frame_Solar_Zenith,
                       SOLAR_AZIMUTH_ZENITH_SCALE_FACTOR,i16_buffer);

  NADIR_FRAME_GEO_VALS("Latitude",
                       L1B_Scan_Meta->Nadir_Frame_Latitude,
                       1.0,f32_buffer);

  NADIR_FRAME_GEO_VALS("Longitude",
                       L1B_Scan_Meta->Nadir_Frame_Longitude,
                       1.0,f32_buffer);

  if (SDend(geo_sd_id) == FAIL)
  {
     returnStatus = MODIS_F_HDF_ERROR;
     L1BErrorMsg(location, returnStatus, NULL, "SDend", 
                 GEOLOCATION_FILE,
                 "Memory or the disk file must have become corrupted.", 
                 True);
     return returnStatus;
  }

  free(SRCA_mode);
  free(i16_buffer);
  free(f32_buffer);
  free(SD_Sun_zenith);
  free(SD_Sun_Azimuth);
  free(Attitude_Angles);

  return(MODIS_S_OK);   
}  

PGSt_SMF_status Calculate_DCR_Change(L1A_granule_t       *L1A_Gran, 
                                     QA_Data_t           *QA, 
                                     L1B_Scan_Metadata_t *L1B_Scan_Meta)
/*
!C****************************************************************
!Description:   Read in the DCR offset values of the last scan in 
                leading granule and all the DCR values of L1A granule.
                Compare the DCR values of the current scan with the 
                same ones in previous scan. Set DCR changes to 1, if 
                they are different. Otherwise set them to 0.

!Input Parameters: 
     L1A_granule_t *L1A_Gran
!Output Parameters:
     QA_Data_t     *QA

!Revision History:
 (continue at top of the file)

 Revision 01.02 8/26/1999
 Added checking if the dcr offset values are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 01.01 1/26/1999
 Two minor bug fixes.  Move DCR_550 = 0 into loop through scans,
 change the last "or" to an "and".
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1/26/1999 
 Original Development
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
  
!Team-unique Header:
    This software is developed by the MODIS Characterization Support
    Team (MCST)for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.
  
!References and Credits:
    The time setting code was borrowed from Paul Fishers timea.c

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.
 
!Design Notes:
  
!END********************************************************************
*/

{
  PGSt_SMF_status returnStatus;
  int32           S, B, B_emiss, D, D_1km, num_scans; 
  int32           DCR_550;
  int32           file_id;
  int16           flag         = 0;
  int8            Leading_DCR_offset[NUM_DCR_VALUE];
  typedef int8    DCR_offset_t[NUM_DCR_VALUE];
  DCR_offset_t    *DCR_offset;
  DCR_offset = (DCR_offset_t*) malloc(MAX_NUM_SCANS * sizeof(DCR_offset_t));
  typedef int8    XDCR_offset_t[NUM_DCR_VALUE]; 
  XDCR_offset_t   *XDCR_offset;
  XDCR_offset = (XDCR_offset_t*) malloc((MAX_NUM_SCANS + 1) * sizeof(XDCR_offset_t));
  PGSt_integer    Version      = 1;
  char            file_name[PGSd_PC_FILE_PATH_MAX];
  char            *location    = "Calculate_DCR_Change";
  
  /* Read DCR_offset from the last scan of the leading granule, if present */

  if (QA->QA_common.missing_leading_granule == False) {
    returnStatus = PGS_PC_GetReference (LEADING_L1A_GRANULE, 
                                        &Version, file_name);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_FILE_NOT_FOUND;
      L1BErrorMsg(location, returnStatus, 
                  "Could not retrieve file name from PCF.",
                  "PGS_PC_GetReference", LEADING_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

      /* Open file */

    file_id = SDstart(file_name, DFACC_RDONLY); /*for sds interfaces*/
    if (file_id == FAIL)
    {
      returnStatus = MODIS_F_FILE_NOT_OPENED;
      L1BErrorMsg(location, returnStatus, 
                  "Could not open file for SD read access.",
                  "SDstart", LEADING_L1A_GRANULE,
                  "The file may be missing, corrupted or "
                  "not an HDF-4 file.", True);
      return returnStatus;
    }

    returnStatus = read_attribute (file_id, "Number of Scans",
                                   DFNT_INT32, (void *)&num_scans);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, 
                  "Could not read Number of Scans.",
                  "read_attribute", LEADING_L1A_GRANULE, 
                  Invalid_MOD01_Msg, True);
      return returnStatus;
    }

    returnStatus = read_part_sds_rank2 (file_id, 
                                        "fpa_dcr_offset", 
                                        num_scans - 1,
                                        0, 
                                        1, 
                                        NUM_DCR_VALUE, 
                                        Leading_DCR_offset);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, 
                  "Could not read fpa_dcr_offset.",
                  "read_part_sds_rank2", LEADING_L1A_GRANULE, 
                  Invalid_MOD01_Msg, True);
      return returnStatus;
    }

    if (SDend(file_id) == FAIL)
    {
       returnStatus = MODIS_F_HDF_ERROR;
       L1BErrorMsg(location, returnStatus, NULL, "SDend", 
                   LEADING_L1A_GRANULE,
                   "Memory or the disk file must have become "
                   "corrupted.", True);
      return returnStatus;
    }

  }

      /* Read DCR_offset from all scans of the middle granule */

  returnStatus = read_sds_rank2(L1A_Gran->sd_id, "fpa_dcr_offset", 
                                L1A_Gran->num_scans,
                                NUM_DCR_VALUE, DCR_offset);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, 
                "Could not read fpa_dcr_offset.",
                "read_sds_rank2", FIRST_L1A_GRANULE, 
                Invalid_MOD01_Msg, True);
    return returnStatus;
  }

      /* Place data in the XDCR_offset array.  If the leading granule
       * is missing, use the values from the first scan of the
       * middle granule to begin the XDCR array.
       */

  if (QA->QA_common.missing_leading_granule == False) {
    for(DCR_550 = 0; DCR_550 < NUM_DCR_VALUE; DCR_550++)
      XDCR_offset[0][DCR_550] = Leading_DCR_offset[DCR_550];
  }
  else {
    for(DCR_550 = 0; DCR_550 < NUM_DCR_VALUE; DCR_550++)
      XDCR_offset[0][DCR_550] = DCR_offset[0][DCR_550];
  }

  for(S = 0; S < L1A_Gran->num_scans; S++) {
    for(DCR_550 = 0; DCR_550 < NUM_DCR_VALUE; DCR_550++)
      XDCR_offset[S+1][DCR_550] = DCR_offset[S][DCR_550];
  }
      
  for(S = 0; S < L1A_Gran->num_scans; S++)
  {
    DCR_550 = 0;
    flag = 0;

    for(D_1km = 0; D_1km < DETECTORS_PER_1KM_BAND; D_1km++)
    {
      for(B = 0; B < NUM_250M_BANDS; B++)
      {     
        for(D = 0; D < BAND_RATIO_AT_RES[INDEX_250M]; D++, DCR_550++)
        {
          if(XDCR_offset[S+1][DCR_550] == XDCR_offset[S][DCR_550])
            QA->QA_refl.change_dc_restore_250m
                [S][B][(DETECTORS_PER_1KM_BAND - D_1km) * 4 - D - 1] = 0;
          else
          {
            QA->QA_refl.change_dc_restore_250m
                [S][B][(DETECTORS_PER_1KM_BAND - D_1km) * 4 - D - 1] = 1;
            flag = 1;
          }
        }
      }
      
      for(B = 0; B < NUM_500M_BANDS; B++)
      { 
        for(D = 0; D < BAND_RATIO_AT_RES[INDEX_500M]; D++, DCR_550++)
        {
          if(XDCR_offset[S+1][DCR_550] == XDCR_offset[S][DCR_550])
            QA->QA_refl.change_dc_restore_500m
                [S][B][(DETECTORS_PER_1KM_BAND -D_1km) * 2 - D - 1] = 0;
          else
          {
            QA->QA_refl.change_dc_restore_500m
                [S][B][(DETECTORS_PER_1KM_BAND -D_1km) * 2 - D - 1] = 1;
            flag = 1;
          }
        }
      } 
        
      for(B = 0; B < NUM_1000M_REFL_BANDS - 1; B++, DCR_550++)
      {
        if(XDCR_offset[S+1][DCR_550] == XDCR_offset[S][DCR_550])
          QA->QA_refl.change_dc_restore_1km
              [S][B][DETECTORS_PER_1KM_BAND - D_1km - 1] = 0;
        else
        {
          QA->QA_refl.change_dc_restore_1km
              [S][B][DETECTORS_PER_1KM_BAND - D_1km - 1] = 1;
          flag = 1;
        }
      }

      B_emiss = 0;
      for(B = 0; B <= NUM_EMISSIVE_BANDS; B++, DCR_550++)
      {
        if(B < BAND31)
        {
          if(B == BAND26)
          {
            if(XDCR_offset[S+1][DCR_550] == XDCR_offset[S][DCR_550])
              QA->QA_refl.change_dc_restore_1km
                  [S][NUM_1000M_REFL_BANDS - 1]
                    [DETECTORS_PER_1KM_BAND - D_1km - 1] = 0;
            else
            {
              QA->QA_refl.change_dc_restore_1km
                  [S][NUM_1000M_REFL_BANDS - 1]
                    [DETECTORS_PER_1KM_BAND - D_1km - 1] = 1;
              flag = 1;
            }
            continue;
          }
   
          if(XDCR_offset[S+1][DCR_550] == XDCR_offset[S][DCR_550])
            QA->QA_emiss.change_dc_restore
                [S][B_emiss][DETECTORS_PER_1KM_BAND - D_1km - 1] = 0;
          else
          {
            QA->QA_emiss.change_dc_restore
                [S][B_emiss][DETECTORS_PER_1KM_BAND - D_1km - 1] = 1;
            flag = 1;
          } 
        }
        
        else
        {
          if(XDCR_offset[S+1][DCR_550]     == XDCR_offset[S][DCR_550] &&
             XDCR_offset[S+1][DCR_550 + 1] == XDCR_offset[S][DCR_550 + 1])
            QA->QA_emiss.change_dc_restore
                [S][B_emiss][DETECTORS_PER_1KM_BAND - D_1km - 1] = 0;
          else
          {
            QA->QA_emiss.change_dc_restore
                [S][B_emiss][DETECTORS_PER_1KM_BAND - D_1km - 1] = 1;
            flag = 1;
          }
          DCR_550++;
        }  

        B_emiss++;
      }
    }/*end D_1km */

    if (flag == 1)
      L1B_Scan_Meta->Bit_QA_Flags[S] |= 0x00002000;
  }/*end S */

  free(DCR_offset);
  free(XDCR_offset);
  
  return(returnStatus);
}

PGSt_SMF_status  Init_QA_Parameters (L1A_granule_t   *L1A_Gran,
                                     L1B_granule_t   *L1B_Gran,
                                     QA_Data_t       *QA)
/*
!C**********************************************************************
!Description:
    This routine initializes the QA values of total numbr of pixels, number of
    valid pixels, number of saturated pixels, number of missing pixels and
    pixels representing negative values below noise.  Initially, all pixels
    are assumed to be valid.  Then it checks missing scans. If there is missing
    scan, the function update these values and set bad data flag to be true.
    Later, in Emissive_Cal and Reflective_Cal, these values will be adjusted 
    on a pixel by pixel basis if some values are determined to be missing, 
    saturated, etc.

!Input Parameters:
    L1A_granule_t   *L1A_Gran    contains L1A granule sd_id and scan type and
                                 scan quality array
    L1B_granule_t   *L1B_Gran    contains number scans and number day mode scans 

!Output Parameters:
    L1B_granule_t   *L1B_Gran    contains number of total pixels, missing pixels,
                                 valid pixels, saturate pixels, negative value
                                 below noise level pixels and bad data flag for
                                 each band

    QA_Data_t       *QA          contains QA data such as number of missing scans, 
                                 number of dead detector EV data, number of sector 
                                 rotation EV data, etc. 

!Revision History:
 (continue at top of the file)

   Revision 01.04  February 7, 2002   Razor Issue 180
   Added counting num_b1_lt_0 
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.03  October 29, 2000
   The missing scans are now flagged in an array in L1A_Gran, Razor issue 142.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.02  September 21, 2000
   Added counting num_rsb_at_night_scans as per Razor issue 137.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   ... (many changes not logged) ...

   Revision 01.01  Nov. 16, 1999
   Added checking missing scans and consequentially update the QA parameters
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

   Revision 01.00  Oct. 1997
   Initial development
   Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

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
  int16  B    = 0;  /* band index within L1A resolution */
  int16  B_38 = 0;
  int16  D_490 = 0; /* detector index within set of MODIS bands (0-489) */
  int16  R    = 0;
  int16  S    = 0;  /* scan index */

  QA->QA_common.num_missing_scans      = 0;
  QA->QA_common.num_rsb_at_night_scans = 0;
  for(D_490 = 0; D_490 < NUM_DETECTORS; D_490++)
  {
    QA->QA_common.num_missing_data_in_scans[D_490]        = 0;
    QA->QA_common.num_dead_detector_EV_data[D_490]        = 0;
    QA->QA_common.num_sector_rotation_EV_data[D_490]      = 0;
    QA->QA_common.num_saturated_EV_data[D_490]            = 0;
    QA->QA_common.num_no_bg_DN_EV_data[D_490]             = 0;
    QA->QA_common.num_moon_in_SVP_TEB_EV_data[D_490]      = 0;
    QA->QA_common.num_bad_dn_star_star_RSB_EV_data[D_490] = 0;
    QA->QA_common.num_exceed_max_for_scaling[D_490]       = 0;
    QA->QA_common.num_nadir_door_closed_EV_data[D_490]    = 0;
    QA->QA_common.num_negative_b1[D_490]                  = 0;
    if(D_490 < NUM_HIGH_RESOLUTION_DETECTORS) 
         QA->QA_common.num_dead_subframe_EV_data[D_490]   = 0;
  }

  for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
    for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++)
    {
      if (R != INDEX_1000M_EMISS || B_38 == MODIS_BAND26_INDEX)
      {
        L1B_Gran->valid_pixels[B_38] = L1B_Gran->num_day_scans * 
                                         DETECT_PER_BAND_AT_RES[R] *
                                         BAND_RATIO_AT_RES[R] * 
                                         EV_1km_FRAMES;
      }
      else 
      {
        L1B_Gran->valid_pixels[B_38] = L1B_Gran->num_scans * 
                                         DETECT_PER_BAND_AT_RES[R] *
                                         BAND_RATIO_AT_RES[R] * 
                                         EV_1km_FRAMES;
      }

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
          /*
           * This resets the number of valid pixels for band 26 to
           * include all scans, both day and night.
           */
      if (B_38 == MODIS_BAND26_INDEX)
        L1B_Gran->valid_pixels[B_38] = L1B_Gran->num_scans * 
                                       DETECT_PER_BAND_AT_RES[R] *
                                       BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

      L1B_Gran->total_pixels[B_38]                      = 
                  L1B_Gran->num_scans * 
                    DETECT_PER_BAND_AT_RES[R] * 
                    BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
      L1B_Gran->saturated_pixels[B_38]                  = 0; 
      L1B_Gran->missing_pixels[B_38]                    = 
                  L1B_Gran->total_pixels[B_38] - 
                    L1B_Gran->valid_pixels[B_38];
      L1B_Gran->negative_value_below_noise_pixels[B_38] = 0;
      L1B_Gran->interpolated_pixels[B_38]               = 0;
      L1B_Gran->dead_detector_pixels[B_38]              = 0;
      L1B_Gran->dead_subframe_pixels[B_38]              = 0;
      L1B_Gran->bad_data_flag[B_38]                     = 0;
    }

      /*
       * Adjust qa counter variables for scans that are completely missing.
       */  

    for (S = 0; S < L1B_Gran->num_scans; S++) 
    {  
      if (L1A_Gran->missing_scan[S] == True)
      {
        B_38 = 0;
        for (R = 0; R < NUM_L1A_RESOLUTIONS; R++)
          for (B = 0; B < L1A_BANDS_AT_RES[R]; B++, B_38++)
          {
            if (strcmp(L1A_Gran->ScanType[S], "Day") == SAME ||
                (R == INDEX_1000M_EMISS && B_38 != MODIS_BAND26_INDEX))
            {
              L1B_Gran->missing_pixels[B_38] += 
                  DETECT_PER_BAND_AT_RES[R] * 
                  BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
              L1B_Gran->valid_pixels[B_38] -= 
                  DETECT_PER_BAND_AT_RES[R] * 
                  BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
            }
/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

  /* 
   * If WRITE_BAND_26_SDS is defined, the valid pixels for band 26
   * is computed based on the number of all scans. So, if the scan
   * is missing, number of the missing pixels and valid pixels should
   * be recalculated even if the scan is in night mode.
   */

            else if (B_38 == MODIS_BAND26_INDEX)
            {
              L1B_Gran->missing_pixels[B_38] += 
                  DETECT_PER_BAND_AT_RES[R] * 
                  BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
              L1B_Gran->valid_pixels[B_38] -= 
                  DETECT_PER_BAND_AT_RES[R] * 
                  BAND_RATIO_AT_RES[R] * EV_1km_FRAMES;
            }

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

            L1B_Gran->bad_data_flag[B_38] = 1;
          }
        QA->QA_common.num_missing_scans++;
      }
    }
     
  return MODIS_S_OK;
}

PGSt_SMF_status Write_L1B_ScanMeta (L1B_Scan_Metadata_t *L1B_Scan_Meta, 
                                    L1A_granule_t       *L1A_Gran,
                                    QA_Data_t           *QA,
                                    boolean             skip_night_hi_res)
/*
!C**************************************************************************
!Description:   This routine writes level 1B scan metadata into EV files and 
                also writes bit QA flags into OBC file.

!Input Parameters:
      L1B_Scan_Metadata_t  *L1B_Scan_Meta   contains L1A scan metadata
      L1A_granule_t        *L1A_Gran        contains number of scans, v_id
                                            and scan quality array
      QA_Data_t            *QA              contains number of thermistor 
                                            outliers
      boolean  skip_night_hi_res   True if and only if all scans are
                                night mode scans and writing of 250m
                                and 500m night mode granules has been
                                turned off

!Output Parameters:
      none
!Revision History:
 (continue at top of the file)

 Revision 02.20  March 27, 2003	 Razor Issue #173 		
 Enclosed the rows in the initializer of array-of-struct "vd_field" with braces
 for ANSI-C compliance. 
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 02.19  November 19, 2001  Razor Issue #169
 Added boolean skip_night_hi_res to input parameters and altered logic 
   so that if skip_night_hi_res is True, output data sets for 250m and 
   500m resolution data are not created.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.18 September 30, 1999
 The structure member "num_thermistor_outliers" is now indexed thorugh scan.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 02.17 August 1999
 Added checking if the scan quality array data are valid in L1A file.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

 Revision 02.16 May 1999
 Wrote bit QA flags into OBC file also
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov) 

 Revision 02.00 April 1997
 Combined V_Data_Create & Write_L1B_ScanMeta of version 1.
 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)

 Revision 01.01 1996/04/05
 Update to match Version 1 Design Document
 Neal Devine(devine@ltpmail.gsfc.nasa.gov)
 John Hannon(hannon@highwire.gsfc.nasa.gov)
 Joan Baden (baden@highwire.gsfc.nasa.gov)

 Revision 01.00 1993
 Initial development
 Geir E. Kvaran(geir@highwire.gsfc.nasa.gov)

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
  PGSt_SMF_status   returnStatus    = MODIS_S_OK;
  char    *location                 = "Write_L1B_ScanMeta";
  PGSt_integer       Version        = 1;
  int16   f                         = 0; /*L1B_file index (0,2)*/
  int16   i                         = 0; /*vdata field index*/
  int32   S                         = 0; /*Scan index*/
  int16   num_fields                = 14;
  int32   vd_id                     = 0;
  int32   list_len                  = 0;
  char    *list                     = NULL;
  int32   scan_number               = 0;
  int32   complete_scan_flag        = 0;
  char8   scan_type[4]              = {' ', ' ', ' ', ' '};
  int32   mirror_side               = 0;
  int32   EV_frames                 = 0;
  int32   nadir_frame               = 0;
  float32 nadir_frame_latitude      = 0;
  float32 nadir_frame_longitude     = 0;
  float32 nadir_frame_solar_azimuth = 0;
  float32 nadir_frame_solar_zenith  = 0;
  int32   BB_thermistor_outliers    = 0;
  uint32  bit_QA_flags              = 0;
  uint16  sector_rotation_dn[MAX_NUM_SCANS];
  float32 sector_rotation_angle;
  char    buffer[64];
  char    *buffer_ptr;
  int32   obc_file_id;
  intn    hdf_return;
  char    file_name[512];
  struct  { char  *name; int32 type; int32 order; } vd_field[14] = 
    {
     {"Scan Number"                      , DFNT_INT32  , 1}, 
     {"Complete Scan Flag"               , DFNT_INT32  , 1},
     {"Scan Type"                        , DFNT_CHAR8  , 4},
     {"Mirror Side"                      , DFNT_INT32  , 1},
     {"EV Sector Start Time"             , DFNT_FLOAT64, 1}, 
     {"EV_Frames"                        , DFNT_INT32  , 1}, 
     {"Nadir_Frame_Number"               , DFNT_INT32  , 1}, 
     {"Latitude of Nadir Frame"          , DFNT_FLOAT32, 1},
     {"Longitude of Nadir Frame"         , DFNT_FLOAT32, 1},
     {"Solar Azimuth of Nadir Frame"     , DFNT_FLOAT32, 1},
     {"Solar Zenith of Nadir Frame"      , DFNT_FLOAT32, 1},
     {"No. OBC BB thermistor outliers"   , DFNT_INT32  , 1},  
     {"Bit QA Flags"                     , DFNT_UINT32 , 1},
     {"Sector Rotation Angle"            , DFNT_FLOAT32, 1}
    };
  int32 evfile_luns[NUM_L1B_EV_FILES] = 
                                       {
                                         L1B_EV_250M_FILE,
                                         L1B_EV_500M_FILE,
                                         L1B_EV_1000M_FILE
                                       };
  int16 start_output_res = INDEX_L1B_250m;
     
#define PACK_MEMBER(member, buffer_ptr) \
memcpy(buffer_ptr, (void *)&member, sizeof(member)); \
buffer_ptr += sizeof(member)

  if (skip_night_hi_res == True)
    start_output_res = INDEX_L1B_1km;
  else
    start_output_res = INDEX_L1B_250m;
  
  for (f = start_output_res; f < NUM_L1B_EV_FILES; f++)
  {

    /*create a new v_data*/

    if ((vd_id = VSattach(L1B_Scan_Meta->v_id[f],-1,"w")) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not attach to vdata for writing Level "
                  "1B Swath Metadata",
                  "VSattach", evfile_luns[f], NULL, True);
      return returnStatus;
    }

    /*give it a name*/

    if (VSsetname(vd_id,"Level 1B Swath Metadata") == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not create vdata name for writing "
                  "Level 1B Swath Metadata",
                  "VSsetname", evfile_luns[f], NULL, True);
      return returnStatus;
    }

    /*define and set fields of the v_data*/

    list_len = 0;
    for (i = 0; i < num_fields; i++)
    {
      if (VSfdefine(vd_id, 
                    vd_field[i].name, 
                    vd_field[i].type, 
                    vd_field[i].order) == FAIL)
      {
        char errmsg[256];
        if (vd_field[i].name)
          sprintf(errmsg, 
                  "Could not define vdata field for \"%s\".", 
                  vd_field[i].name);
        else
          strcpy(errmsg, "Could not define vdata field -- name is NULL.");
        returnStatus = MODIS_F_HDF_ERROR;
        L1BErrorMsg(location, returnStatus, errmsg,
                    "VSfdefine", evfile_luns[f], NULL, True);
        return returnStatus;
      }

      /*total num of char's in all names*/
      list_len += strlen(vd_field[i].name);
    }

    /*allocate memory for char's & commas & \0*/
    list = (char *)malloc((list_len + num_fields + 1) * sizeof(char)); 
    if (!list)
    {
      returnStatus = MODIS_F_OUT_OF_MEMORY;
      L1BErrorMsg(location, returnStatus, NULL, "malloc", 0, NULL, True);
      return returnStatus;
    }

    /* Initialize list */

    list[0] = '\0';

    /*Append names and commas*/

    for (i = 0; i < num_fields - 1; i++)
    {
      strcat(list,vd_field[i].name);
      strcat(list,",");
    }

    /*Last name is not followed by a comma*/

    strcat(list,vd_field[i].name);
 
    /*set fields*/

    if (VSsetfields(vd_id,list) == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, NULL,
                  "VSsetfields", evfile_luns[f], NULL, True);
      return returnStatus;
    }

    /* read sector rotation dn */
    /* XG 3/12/2015 
     * use SET_FR_ENC_DELTA from "Command Parameters" field for sector rotation status 
     */
    /*
    returnStatus = read_vdata(L1A_Gran->v_id, 0, 
                              L1B_Scan_Meta->num_scans, 
                              "Telemetry Major Cycle 6 of 7", 
                              "CS_FR_ENC_DELTA", 
                              (VOIDP)sector_rotation_dn);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "read_vdata", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
    returnStatus = Check_Valid_Range("CS_FR_ENC_DELTA", 
                                     DFNT_UINT16, 
                                     NULL, 
                                     "16383", 
                                     NULL,
                                     L1A_Gran->num_scans, 
                                     (void *) sector_rotation_dn);
    */

    returnStatus = read_vdata(L1A_Gran->v_id, 0, 
                              L1B_Scan_Meta->num_scans, 
                              "Command Parameters", 
                              "SET_FR_ENC_DELTA", 
                              (VOIDP)sector_rotation_dn);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "read_vdata", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
    returnStatus = Check_Valid_Range("SET_FR_ENC_DELTA", 
                                     DFNT_UINT16, 
                                     NULL, 
                                     "16383", 
                                     NULL,
                                     L1A_Gran->num_scans, 
                                     (void *) sector_rotation_dn);

    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL, "Check_Valid_Range",
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    for (S = 0; S < L1B_Scan_Meta->num_scans; S++)
    {
      buffer_ptr = buffer;

      scan_number = S + 1;
      PACK_MEMBER(scan_number, buffer_ptr);

      if (L1A_Gran->scan_quality[S][1] != 0)
        complete_scan_flag = 0;
      else
        complete_scan_flag = 1;


      PACK_MEMBER(complete_scan_flag, buffer_ptr);

      if (strcmp(L1B_Scan_Meta->ScanType[S],"Day") == SAME) 
        scan_type[0] = 'D';
      else if (strcmp(L1B_Scan_Meta->ScanType[S],"Night") == SAME) 
        scan_type[0] = 'N';
      else
        scan_type[0] = 'O';
      memcpy(buffer_ptr, (void *)scan_type, sizeof(scan_type));
      buffer_ptr += sizeof(scan_type);

      mirror_side               = 
          L1B_Scan_Meta->MirrorSide[S];
      PACK_MEMBER(mirror_side, buffer_ptr);

      PACK_MEMBER(L1B_Scan_Meta->EVStartTime_TAIsecond[S], buffer_ptr);

      EV_frames                 = 
          L1B_Scan_Meta->EV_frames;
      PACK_MEMBER(EV_frames, buffer_ptr);

      nadir_frame               = 
          L1B_Scan_Meta->Nadir_Frame_Number;
      PACK_MEMBER(nadir_frame, buffer_ptr);

      nadir_frame_latitude      = 
          L1B_Scan_Meta->Nadir_Frame_Latitude[S];
      PACK_MEMBER(nadir_frame_latitude, buffer_ptr);

      nadir_frame_longitude     = 
          L1B_Scan_Meta->Nadir_Frame_Longitude[S];
      PACK_MEMBER(nadir_frame_longitude, buffer_ptr);

      nadir_frame_solar_azimuth = 
          L1B_Scan_Meta->Nadir_Frame_Solar_Azimuth[S];
      PACK_MEMBER(nadir_frame_solar_azimuth, buffer_ptr);

      nadir_frame_solar_zenith  = 
          L1B_Scan_Meta->Nadir_Frame_Solar_Zenith[S];
      PACK_MEMBER(nadir_frame_solar_zenith, buffer_ptr);

      BB_thermistor_outliers    = 
          (int32) QA->QA_emiss.num_thermistor_outliers[S];
      PACK_MEMBER(BB_thermistor_outliers, buffer_ptr);

      bit_QA_flags              = 
          L1B_Scan_Meta->Bit_QA_Flags[S];
      PACK_MEMBER(bit_QA_flags, buffer_ptr);
     
      sector_rotation_angle     = 
          (180.0/8192.0)*(float32)sector_rotation_dn[S];
      PACK_MEMBER(sector_rotation_angle, buffer_ptr);
 
      if (VSwrite(vd_id, 
                  (unsigned char *)buffer, 
                  1, 
                  FULL_INTERLACE) == FAIL)
      {
        returnStatus = MODIS_F_WRITE_ERROR;
        L1BErrorMsg(location, returnStatus, NULL,
                    "VSwrite", evfile_luns[f], NULL, True);
        return returnStatus;
      }

    }  /* S */

    free(list);
      
    VSdetach(vd_id);
  
  }/* f */
 
  /*
   * Convert logical ID to physical file name.
   */

  returnStatus = PGS_PC_GetReference (L1B_OBC_FILE, &Version, file_name);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, 
                "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Open file for Write.  The file should already exist since most of the
   * data was written in the Preprocess module (the file was closed there).
   */

  obc_file_id = SDstart (file_name, DFACC_RDWR);   /* for SD interfaces */
  if (obc_file_id == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, 
                "Could not open file for SD read/write access.",
                "SDstart", L1B_OBC_FILE,
                "The file may be missing, corrupted or "
                "not an HDF-4 file.", True);
    return returnStatus;
  }

  /*
   * Write the Bit_QA_Flags as a 1D SDS to the OBC file.
   */

  returnStatus = write_sds_rank1(obc_file_id, 
                                 "Bit QA Flags", 
                                 NUM_SCANS_DIM_NAME,
                                 L1B_Scan_Meta->num_scans, 
                                 "uint32", 
                                 L1B_Scan_Meta->Bit_QA_Flags);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "write_sds_rank1",
                L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*
   * Close the file.
   */

  hdf_return = SDend(obc_file_id);
  if (hdf_return == FAIL)
  {
     returnStatus = MODIS_F_HDF_ERROR;
     L1BErrorMsg(location, returnStatus, NULL, "SDend", L1B_OBC_FILE,
                 "Memory or the disk file must have become corrupted.", 
                 True);
     return returnStatus;
  }
  
  return(MODIS_S_OK);
}


PGSt_SMF_status Determine_Other_Missing_Scans (lookup_tables_t *tables,
                                               L1A_granule_t   *L1A_Gran)
/*
!C**************************************************************************
!Description:
   This function examines validity of certain L1A data and, based on options
   defined in the LUTs, determines additional scans to be treated as
   completely missing (meaning that data will not be calibrated for any
   band of the scan).

!Input Parameters:
   lookup_tables_t *tables           Pointer to all L1B LUTs.  The common
                                     QA LUT "control_options" is used.
   L1A_granule_t   *L1A_Gran         All members are filled, including
                                     structure member "missing_scan", based
                                     on normal values of L1A data.

!Output Parameters:
   L1A_granule_t   *L1A_Gran         structure member "missing_scans" may
                                     be adjusted based on the data.

!Revision History:
 (continue at top of the file)

   Revision 01.00 October 29, 2000
   Initial development.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

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
  boolean split_scan[MAX_NUM_SCANS];
  boolean bad_scan_quality[MAX_NUM_SCANS];
  int32 S;
  char *location = "Determine_Other_Missing_Scans";
  common_QA_tables_t *common_QA_LUT = &tables->QA.common_QA_tables;

  /*
   * Set all elements of the array "split_scan" to True or False.
   */

  if (common_QA_LUT->control_options[SPLIT_SCAN_CONTROL] == ON) 
  {
    returnStatus = Determine_Split_Scans(L1A_Gran, split_scan);
    if (returnStatus != MODIS_S_OK) {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Determine_Split_Scans", 0, NULL, True);
    }
  }
  else {
    for (S = 0; S < L1A_Gran->num_scans; S++)
      split_scan[S] = False;
  }

  /*
   * Set all elements of the array "bad_scan_quality" to True or False.
   */

  if (common_QA_LUT->control_options[BAD_SCAN_QUALITY_CONTROL] == ON) 
  {
    for (S = 0; S < L1A_Gran->num_scans; S++) {
      if ((L1A_Gran->scan_quality[S][0] != 0 &&
           L1A_Gran->scan_quality[S][0] != 1) ||
          L1A_Gran->scan_quality[S][1] < 0  ||
          L1A_Gran->scan_quality[S][2] < 0  ||
          L1A_Gran->scan_quality[S][3] < 0)
        bad_scan_quality[S] = True;
      else
        bad_scan_quality[S] = False;
    }
  }
  else 
  {
    for (S = 0; S < L1A_Gran->num_scans; S++)
      bad_scan_quality[S] = False;
  }

  /*
   * adjust the values of "missing_scan", if necessary.
   */

  for (S = 0; S < L1A_Gran->num_scans; S++) 
  {
    if (split_scan[S] || bad_scan_quality[S])
      L1A_Gran->missing_scan[S] = True;
  }


  /*
   * Check the mirror side flags.  So far, there has not been a
   * situation where the mirror side flag was not valid but there were
   * data on the scan.  However, just in case this happens, this will
   * allow the scan to be skipped rather than causing the code to error
   * out.
   */


  for (S = 0; S < L1A_Gran->num_scans; S++) 
  {
    if (!(L1A_Gran->MirrorSide[S] == 0 || L1A_Gran->MirrorSide[S] == 1))
      L1A_Gran->missing_scan[S] = True;
  }


  return MODIS_S_OK;
}

PGSt_SMF_status Determine_Split_Scans(L1A_granule_t *L1A_Gran,
                                      boolean *split_scan)
/*
!C**************************************************************************
!Description:
   This function determines scan indices of "split scans".  A split scan
   arises from a bit flip in the scan time of a packet.  The L1A code takes
   the out-of-order packet and starts a new scan in the granule.  Thus, the
   original scan is split into two parts.  A split scan causes the number
   of scans in the granule to exceed the normal 203 or 204 scans.

   See the design notes for comments on the algorithm for detecting a split
   scan.

!Input parameters:

   L1A_granule_t   *L1A_Gran       contains L1A variables for SD start
                                   time, scan quality and mirror side

!Output parameters:

   boolean         split_scan      array which defines if a scan

!Revision History:
 (continue at top of the file)

   Revision 01.02, October 29, 2000
   Moved from L1B.c and removed adjusting the QA counters. (issue 142)
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.01 September 20, 2000
   Correct logic errors which led to PGE02 failures (see Razor issue 136).
   Call new function, Get_Split_Scan_Indexes.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00 May 2, 2000
   Initial development
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

   The algorithm for finding the two parts of a split scan begins by
   finding an abnormally short time interval between two scans.  The
   "SD start time" is used for this because in examining many cases of
   split scans, it appeared to be the most reliable time among the
   various sector start times ("EV start time" was the least reliable).

   The two scans with the abnormally short time interval are not
   necessarily the two parts of the split scan.  If there are no
   repeated mirror side indexes or no missing packets, then the
   abnormal time interval is assumed to not be a marker for a split
   scan.

   If the short time interval does correspond to a split scan, then at
   least one of the two scans will be part of the split scan.  To find
   the other part, we expand the range of scan indexes by 1 in each
   direction and look at "Mirror side" and "Scan quality array".  A
   repeated mirror side index is the next most reliable indicator.
   However, occasionally, one part of the split scan will have a mirror
   side index of -1.  In this case, we then turn to the number of
   missing packets (recorded in the second element of the scan quality
   array).  The two adjacent scans with large numbers of missing packets
   are the next most reliable indicator.

   There is no attempt to find the two parts of a split scan if they
   reside in different granules.

!END********************************************************************
*/
{
#define SCAN_INTERVAL  1.4771
  PGSt_SMF_status returnStatus;
  int32 i;         /* scan index in loop through scans */
  int32 j;         /* used for scan interval */
  int32 S_split_1; /* scan index of 1st part of split scan */
  int32 S_split_2; /* scan index of 2nd part of split scan */
  float64 SD_start_time[MAX_NUM_SCANS];
  char *location = "Determine_Split_Scans";

  for (i = 0; i < L1A_Gran->num_scans; i++)
  {
    split_scan[i]   = False;
  }

  /* Read SD start time */

  returnStatus = read_sds_rank1 (L1A_Gran->sd_id, "SD start time",
                                 L1A_Gran->num_scans,
                                 (void *)SD_start_time);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, 
                "Could not read SD start time.",
                "read_sds_rank1", 
                FIRST_L1A_GRANULE, Invalid_MOD01_Msg, True);
    return returnStatus;
  }

  /* Loop through all scans looking for potential split scans */

  for (i = 1; i < L1A_Gran->num_scans; i++)
  {

    /*
     * Look at the time interval between scans i-1 and i.  If the time
     * interval is smaller than 1/2 the normal time interval between
     * MODIS scans, then define this to be "abnormally short".
     */

    j = (int32)((SD_start_time[i] - SD_start_time[i-1])/SCAN_INTERVAL + 0.5);
    if (j <= 0)
    {

      /*
       * An abnormally short time interval was detected between scans
       * i-1 and i.  Get the two indexes of the split scan.  If the
       * indexes return as invalid indexes, assume that the abnormally
       * short time interval does not correspond to a split scan.
       */

      returnStatus = Get_Split_Scan_Indexes(i-1, 
                                            L1A_Gran->num_scans,
                                            L1A_Gran->MirrorSide,
                                            L1A_Gran->scan_quality,
                                            &S_split_1, 
                                            &S_split_2);
      if (returnStatus != MODIS_S_OK)
      {
        L1BErrorMsg(location, returnStatus, NULL,
                    "Get_Split_Scan_Indexes", 0, NULL, False);
        return MODIS_S_OK;
      }

      if (S_split_1 >= 0 && S_split_1 < L1A_Gran->num_scans)
        split_scan[S_split_1] = True;
      if (S_split_2 >= 0 && S_split_2 < L1A_Gran->num_scans)
        split_scan[S_split_2] = True;

    }
  }

  return MODIS_S_OK;
}

PGSt_SMF_status Get_Split_Scan_Indexes 
                       (int32 S1,
                        int32 num_scans,
                        int16 mirror_side[],
                        int32 scan_quality[]
                                [SCAN_QUALITY_ARRAY_NUM_ELEMENTS],
                        int32 *S_split_1,
                        int32 *S_split_2)
/*
!C**************************************************************************
!Description:
   Given a pair of scans with an abnormally short time interval,
   examine the mirror side and scan quality array indexes to find the
   actual split scan indexes.  If the split scan indexes cannot be
   determined, return invalid values (-1) for the indexes.

!Input parameters:
   int32 S1                 1st scan index of the pair of scans with an
                            abnormally short time interval
   int32 num_scans          Number of scans in the granule
   int32 mirror_side[]      mirror side indexes
   int32 scan_quality[][]   Scan quality array

!Output parameters:
   int32 *S_split_1         Index of 1st part of split scan (or -1
                            if there is no split scan)
   int32 *S_split_2         Index of 2nd part of split scan (or -1
                            if there is no split scan)

!Revision History:
 (continue at top of the file)

   Revision 01.01  October 16, 2004  Razor Issue #200
   Casted Int32 variables in sprintf calls to "long" with the
   format specifier "%ld" for better code portability.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.00 September 20, 2000
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
    This software is developed by the MODIS Science Data Support
    Team for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END********************************************************************
*/
{
  int32 S, S_begin, S_end;
  int32 num_invalid_ms, n, max_n;
  char *location = "Get_Split_Scan_Indexes";

  /*
   * Check S1 and num_scans
   */

  if (S1 < 0 || S1 >= num_scans || num_scans <= 0)
  {
    char errmsg[256];
    sprintf (errmsg, 
             "Input error in S1 or num_scans:\nS1 = %ld\nnum_scans = %ld\n",
             (long)S1, (long)num_scans);
    L1BErrorMsg(location, MODIS_F_INVALID_ARGUMENT, errmsg,
                NULL, 0, NULL, True);
    return MODIS_F_INVALID_ARGUMENT;
  }

  /*
   * Initialize the return indexes to invalid values.
   */

  *S_split_1 = -1;
  *S_split_2 = -1;

  /*
   * Find the begin and end scan indexes of the range of scans
   * to look through.  The normal range is a set of 4 scans, starting
   * at the scan just before S1.
   */

  S_begin = S1 - 1;
  if (S_begin < 0) S_begin = 0;
  S_end = S1 + 2;
  if (S_end >= num_scans) S_end = num_scans - 1;

  /*
   * determine if any of the scans has an invalid mirror side index.
   */

  num_invalid_ms = 0;
  for (S = S_begin; S <= S_end; S++) {
    if (!(mirror_side[S] == 0 || mirror_side[S] == 1))
      num_invalid_ms++;
  }

  /*
   * If all the mirror side indexes are valid, base the determination
   * on mirror side.  Otherwise, base the determination on the pair
   * with the largest number of missing packets.
   */

  if (num_invalid_ms == 0) {
    for (S = S_begin; S < S_end; S++) {
      if (mirror_side[S] == mirror_side[S+1]) {
        *S_split_1 = S;
        *S_split_2 = S+1;
        break;
      }
    }
  }
  else {
    max_n = 0;
    for (S = S_begin; S < S_end; S++) {
      n = 0;
      if (scan_quality[S][1] > 0)
        n+= scan_quality[S][1];
      if (scan_quality[S+1][1] > 0)
        n+= scan_quality[S+1][1];
      if (n > max_n) {
        max_n = n;
        if (scan_quality[S][0] == 0) {
          *S_split_1 = S;    /* extra scan, blank */
          *S_split_2 = -1;   /* assume this one is OK */
        }
        else if (scan_quality[S+1][0] == 0) {
          *S_split_1 = -1;   /* assume this one is OK */
          *S_split_2 = S+1;  /* extra scan, blank */
        }
        else {
          *S_split_1 = S;
          *S_split_2 = S+1;
        }
      }
    }
    if (max_n < 2800) {
      *S_split_1 = -1;
      *S_split_2 = -1;
    }
  }

  return MODIS_S_OK;
}


PGSt_SMF_status Set_UI_ConvertToPercent_Attrs
                            (lookup_tables_t  *tables,
                             L1B_Scan_t       *L1B_Scan,
                             boolean          skip_night_hi_res)
/*
!C****************************************************************************
!Description:
   Add attributes to all EV UI SDSs which allow conversion of UI to
   percent uncertainty.  These attributes are: "specified_uncertainty",
   "scaling_factor" and "uncertainty_units".

!Input Parameters:
   lookup_tables_t  *tables      Contains values for the attributes in LUTs
   L1B_Scan_t       *L1B_Scan    Contains the UI SDS IDs (opened)
   boolean          skip_night_hi_res   True if and only if all scans are 
                                     night mode scans and writing of 250m
                                     and 500m night mode granules has been
                                     turned off

!Output Parameters:
   (none -- the HDF files themselves are modified)

!Revision History:
 (continue at top of the file)

  Revision 01.02  January 17, 2002 Razor Issue #171
  Improve portability of code to 64-bit mode.
  Cast strlen returns to int32 in calls to SDSetattr.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.01, November 16, 2001, Razor Issue #169
     Added input parameter skip_night_hi_res and changed logic so that
     250m and 500m data are not written when granule has no day mode
     scans and runtime parameter Write_Night_Mode_HiRes_Data is False.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.00  Nov. 29, 2000
   Initial development
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Science Data Support
   Team for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   The arrays contained in the LUTs are for all emissive bands (TEB LUTs)
   and all reflective bands (RSB LUTs).  Within these linear arrays, the
   numbers to be set in the attributes are consecutive.  Thus, use a macro
   where the beginning index and number of values define which subset of
   the arrays to write into the attributes.

!END**************************************************************************
*/
{
  PGSt_SMF_status returnStatus;
  char *units_attr_value = "percent";
  char *location = "Set_UI_ConvertToPercent_Attrs";

#define REFL_250M_START_INDEX    0
#define REFL_500M_START_INDEX    2
#define REFL_1000M_START_INDEX   7
#define REFL_BAND26_START_INDEX 21

/* Macro: SET_UI_SDS_CONV_ATTRS
 *    sds_id       = opened SDS ID for the SDS to set attributes for
 *    spec_uncert  = variable containing specified uncertainty float array
 *    scale_factor = variable containing scaling factor float array
 *    indx         = the starting element index in the arrays
 *    numval       = number of values in the attribute
 */
#define SET_UI_SDS_CONV_ATTRS(sds_id, spec_uncert, scale_factor, indx,     \
                              numval, lun)                                 \
  if (SDsetattr(sds_id, "specified_uncertainty", DFNT_FLOAT32,             \
                numval, &(spec_uncert[indx])) == FAIL) {                   \
    returnStatus = MODIS_F_WRITE_ERROR;                                    \
    L1BErrorMsg(location, returnStatus,                                    \
                "Could not write UI SDS attribute"                         \
                     "\"specified_uncertainty\".",                         \
                "SDsetattr", lun, NULL, True);                             \
    return returnStatus;                                                   \
  }                                                                        \
  if (SDsetattr(sds_id, "scaling_factor", DFNT_FLOAT32,                    \
                numval, &(scale_factor[indx])) == FAIL) {                  \
    returnStatus = MODIS_F_WRITE_ERROR;                                    \
    L1BErrorMsg(location, returnStatus,                                    \
                "Could not write UI SDS attribute \"scaling_factor\".",    \
                "SDsetattr", lun, NULL, True);                             \
    return returnStatus;                                                   \
  }                                                                        \
  if (SDsetattr(sds_id, "uncertainty_units", DFNT_CHAR8,                   \
             (int32)strlen(units_attr_value), units_attr_value) == FAIL) { \
    returnStatus = MODIS_F_WRITE_ERROR;                                    \
    L1BErrorMsg(location, returnStatus,                                    \
                "Could not write UI SDS attribute \"uncertainty_units\".", \
                "SDsetattr", lun, NULL, True);                             \
    return returnStatus;                                                   \
  }

  if (skip_night_hi_res == False)
  {
    SET_UI_SDS_CONV_ATTRS(L1B_Scan->UI_sds_id[INDEX_250M], 
                          tables->refl.RSB_specified_uncertainty, 
                          tables->refl.RSB_UI_scaling_factor, 
                          REFL_250M_START_INDEX, NUM_250M_BANDS, 
                          L1B_EV_250M_FILE);

    SET_UI_SDS_CONV_ATTRS(L1B_Scan->UI_sds_id[INDEX_500M], 
                          tables->refl.RSB_specified_uncertainty, 
                          tables->refl.RSB_UI_scaling_factor, 
                          REFL_500M_START_INDEX, NUM_500M_BANDS, 
                          L1B_EV_500M_FILE);

    SET_UI_SDS_CONV_ATTRS(L1B_Scan->EV_250m_Aggr500m_RefSB_UI_sds_id, 
                          tables->refl.RSB_specified_uncertainty, 
                          tables->refl.RSB_UI_scaling_factor, 
                          REFL_250M_START_INDEX, NUM_250M_BANDS, 
                          L1B_EV_500M_FILE);
  }

  SET_UI_SDS_CONV_ATTRS(L1B_Scan->UI_sds_id[INDEX_1000M_REFL], 
                        tables->refl.RSB_specified_uncertainty, 
                        tables->refl.RSB_UI_scaling_factor, 
                        REFL_1000M_START_INDEX, NUM_1000M_REFL_BANDS, 
                        L1B_EV_1000M_FILE);

  SET_UI_SDS_CONV_ATTRS(L1B_Scan->UI_sds_id[INDEX_1000M_EMISS], 
                        tables->emiss.TEB_specified_uncertainty, 
                        tables->emiss.TEB_UI_scaling_factor, 
                        0, NUM_EMISSIVE_BANDS, 
                        L1B_EV_1000M_FILE);

  SET_UI_SDS_CONV_ATTRS(L1B_Scan->EV_250m_Aggr1km_RefSB_UI_sds_id, 
                        tables->refl.RSB_specified_uncertainty, 
                        tables->refl.RSB_UI_scaling_factor, 
                        REFL_250M_START_INDEX, NUM_250M_BANDS, 
                        L1B_EV_1000M_FILE);

  SET_UI_SDS_CONV_ATTRS(L1B_Scan->EV_500m_Aggr1km_RefSB_UI_sds_id, 
                        tables->refl.RSB_specified_uncertainty, 
                        tables->refl.RSB_UI_scaling_factor, 
                        REFL_500M_START_INDEX, NUM_500M_BANDS, 
                        L1B_EV_1000M_FILE);

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS

  SET_UI_SDS_CONV_ATTRS(L1B_Scan->Band26.UI_sds_id, 
                        tables->refl.RSB_specified_uncertainty, 
                        tables->refl.RSB_UI_scaling_factor, 
                        REFL_BAND26_START_INDEX, 1, 
                        L1B_EV_1000M_FILE);

#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

  return MODIS_S_OK;
}

PGSt_SMF_status Calculate_B26_B5_Correction
                (float32            *original_correction,
                 float32            *scaled_correction,
                 L1B_ScaleOffset_t  *ScaleOffset)

/*
!C****************************************************************************
!Description:
  This function scales crosstalk correction factors which are used to adjust 
  Band 26 data by crosstalk from Band 5.  The scaling factor is the Band 5
  radiance scale divided by the Band 26 radiance scale.

!Input Parameters:
   float32           *original_correction  Lookup table correction values
   L1B_ScaleOffset_t *ScaleOffset          Granule-specific scales and offsets
     
!Output Parameters:
   float32           *scaled_correction    Correction values scaled by the
                                           granule-specific scales and offsets

!Revision History:
 (continue at top of the file)

 Initial development  March 19, 2002      Razor Issue #182
 Added to Aqua code   March 26, 2003      Razor Issue #190
 Alice Isaacman, SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

!Team-unique Header:
 This software is developed by the MODIS Characterization Support
 Team (MCST)for the National Aeronautics and Space Administration,
 Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
 HDF portions developed at the National Center for Supercomputing
 Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:

!END**************************************************************************
*/

{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  int16    D            = 0;     /* Detector */     
  float32  scale_factor = 1.;    /* Scaling factor */


  scale_factor = ScaleOffset->rad_scale_RefSB[MODIS_BAND5_INDEX]/
                   ScaleOffset->rad_scale_RefSB[BAND_26_REFL_INDEX];
  for (D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    scaled_correction[D] = scale_factor*original_correction[D];

  return returnStatus;
}   
