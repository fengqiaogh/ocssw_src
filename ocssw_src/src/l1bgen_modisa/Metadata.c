/*****************************************************************************

File: Metadata.c

External functions:
   Write_Gran_Metadata
   Gran_Meta_Cal

Other functions:
   Write_Global_Metadata
   Get_Electronics_Status
   Get_Elec_Config_Status_Per_Gran
   Get_Elec_Config_Status

Revision History:
   $Log: Metadata.c,v $
   Revision 1.34  2017-08-02 10:41:20-04  xgeng
   Aqua V6.2.1 update for C6.1 reprcoessing

   Revision 1.33  2014-10-02 13:36:48-04  xgeng
   Added doi attributes to NRT product

   Revision 1.32  2012-06-26 15:07:42-04  xgeng
   added the doi metadata

   Revision 1.27  2009-08-31 13:39:45-04  xgeng
   Corrected a bug in the output of dead subframe list and noisy subframe list

   Revision 1.26  2009/07/24 20:00:37  xgeng
   Added an extension "NRT" to the LOCALGRANULEID to indentify the Near-Real-Time product.

   Revision 1.21  2008/11/18 16:30:09  xgeng
   merge the branch for V6.0.1

   Revision 1.20.1.3  2008/06/02 15:44:09  xgeng
   Calculated the deadSubframeDataPercent for 250m and 500m bands; Added the effect of dead subframe on the bad_data_flag and uncalibratedDataPercent; Output Detector Quality Flag2, deadSubframeDataPercent, dead_subframe, noisy_subframe.

   Revision 1.20  2006/10/27 14:35:32  ltan
   Changed to write PGEVersion value extracted from the PCF file (LUN 800500). Some comments regarding Night mode handling corrected.


*****************************************************************************/

/*------------------------------------------------------------------
  Must include PGS_MET.h before hdf includes as a typedef of
  ATTRIBUTE in CUC/odldef.h conflicts with a #def in hdf includes.
------------------------------------------------------------------*/
#include  "MetadataP.h"
#include  "Preprocess.h"
#include  "L1B_Tables.h"
#include  "PGS_TD.h"
#include  "HDF_Lib.h"
#include  "FNames.h"
#include  "PGS_Error_Codes.h"
#include  <time.h>
#include  <libgen.h>
#include  "granule_extract_metadata.h"

#define TIMECODEASIZE        28
#define MECS_CORE    "CoreMetadata.0"
#define MECS_ARCH    "ArchiveMetadata.0"

#define CALIBRATIONQUALITY_MACRO       "marginal"
#define NADIRPOINTING_MACRO            "Y"
#define ALGORITHMPACKAGENAME_MACRO     "MODIS Level 1B Algorithm Package"
#define AUTOMATICQUALITYFLAG_MACRO 		 "Suspect" 
#define AUTOMATICQUALITYFLAGEXPLANATION_MACRO	"not being investigated" 
#define ANCILLARYINPUTTYPE_MACRO 		   "Geolocation"
#define INSTRUMENTNAME_MACRO	"Moderate Resolution Imaging SpectroRadiometer"

pgs_meta_t pgs_in;
char *pgs_out_mdHandle;

PGSt_SMF_status  Write_Gran_Metadata 
                      (Run_Time_Parameters_t *runtime_params,
                       L1B_Gran_Metadata_t   *L1B_Gran_Meta,
                       QA_Data_t             *QA,
                       Preprocess_Data_t     *PP,
                       lookup_tables_t       *tables,
                       L1A_granule_t         *L1A_Gran,
                       boolean               skip_night_hi_res)
/*
!C****************************************************************
!Description:   Read metadata from input files and write to L1B files.

!Input Parameters:
     L1B_Gran_Metadata_t *L1B_Gran_Meta   Contains metadata
     QA_Data_t *QA                        Contains metadata 
     Preprocess_Data_t *PP                Contains metadata
     lookup_tables_t *tables              Contains metadata
     L1A_granule_t *L1A_Gran              Contains satellite id 
     boolean skip_night_hi_res       True if and only if all scans are
                                       NIGHT mode and runtime parameter
                                       Write_Night_Mode_HiRes_Data is False.   
!Output Parameters:
     none

!Revision History:
 $Log: Metadata.c,v $
 Revision 1.34  2017-08-02 10:41:20-04  xgeng
 Aqua V6.2.1 update for C6.1 reprcoessing

 Revision 1.33  2014-10-02 13:36:48-04  xgeng
 Added doi attributes to NRT product

 Revision 1.32  2012-06-26 15:07:42-04  xgeng
 added the doi metadata

 Revision 1.27  2009-08-31 13:39:45-04  xgeng
 Corrected a bug in the output of dead subframe list and noisy subframe list

 Revision 1.26  2009/07/24 20:00:37  xgeng
 Added an extension "NRT" to the LOCALGRANULEID to indentify the Near-Real-Time product.

 Revision 1.21  2008/11/18 16:30:09  xgeng
 merge the branch for V6.0.1

 Revision 1.20.1.3  2008/06/02 15:44:09  xgeng
 Calculated the deadSubframeDataPercent for 250m and 500m bands; Added the effect of dead subframe on the bad_data_flag and uncalibratedDataPercent; Output Detector Quality Flag2, deadSubframeDataPercent, dead_subframe, noisy_subframe.

 (continue at top of the file)

 Revision 1.7  2005/11/21 22:04:11  seadas
 MDM: fixed minor bug in the LOCALGRANULEID field so now only the filename
      and not the path is written

 Revision 1.6  2005/10/27 20:59:38  seadas
 MAR: Updated from swdev 10/27/05.

 Revision 02.17 October 31, 2003   Razor Issue #195
 Changed the source of PROCESSINGCENTER from the Run Time Parameters
 read from PCF file instead of from QA LUTs.
 Liqin Tan   SAIC GSO   (ltan@saicmodis.com)

 Revision 02.16 September 15, 2003   Razor Issue #181
 Added PROCESSINGENVIRONMENT to Archived Metadata.
 Liqin Tan   SAIC GSO   (ltan@saicmodis.com)

 Revision 02.15  September 3, 2002    Razor Issue #188
 Change the form of ProductionHistory Archived Metadata variable.
 Liqin Tan   SAIC GSO   (ltan@saicmodis.com)

 Revision 02.14  December 13, 2002    Razor Issue #189
 Added ProductionHistory Archived Metadata variable
 generation.
 Alice Isaacman   SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.13  July 8, 2002    Razor Issue #185
 Change source of copied metadata to Geolocation file.
 Alice Isaacman   SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.11 November 12, 2001 (Razor Issue #169)
  Added input parameter skip_night_hi_res and changed logic so that
  250m and 500m data are not written when granule has no day mode
  scans and runtime parameter Write_Night_Mode_HiRes_Data is False.
  Moved output_file_indices_t to MetadataP.h.
 Alice Isaacman   SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.10 1998/4/8 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
 Added functions to set, get, copy attribute based on original macros.
 Added new granule metadata.
  
 Revision 01.03 1997/6/4 Zhidong Hao (hao@barebackride.gsfc.nasa.gov)
 
 Revision 01.02 1997/2/27 Neal Devine (neal.devine@gsfc.nasa.gov)
 Rewrote ECSmetadata code for clarity and consistent with CCR.

 Revision 01.01 1996/09/26 Neal Devine (neal.devine@gsfc.nasa.gov)
 Added multiple input and ancilinput files, taken from pcf, to metadata.

 Revision 01.00 1996/04/05 Neal Devine (neal.devine@gsfc.nasa.gov)
 Original Development

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
  PGSt_PC_Logical       in_fileID    = 0;
  intn                  hdf_return   = FAIL;
  int32                 i            = 0; /* for general indexing */
  PGSt_integer          Version      = 1;
  char                  OBC_fname[PGSd_PC_FILE_PATH_MAX];

  PGSt_SMF_status       returnStatus = MODIS_S_OK;
  char                  *location = "Write_Gran_Metadata";
  PGSt_MET_all_handles  mdhandles;
  PGSt_integer          version      = 1;
  int32                 VersionID;
  struct tm             *curtime   = NULL;
  time_t                t;
  char                  *ptr = NULL;
  char                  beginDateTimeA[30];
  char                  beginDateTimeB[30];
  char                  procDateTimeA[30];
  char                  procDateTimeB[30];
  char                  shortname[10];
  char                  sub_shortname[6];
  char                  ProductionHistory[MAX_PRODUCTIONHISTORY_SIZE];
  char                  MCSTversion[MAX_MCST_VERSION_BUFFER];
  char                  MOD01_ProductionHistory[MAX_PRODUCTIONHISTORY_SIZE];
  char                  MOD03_ProductionHistory[MAX_PRODUCTIONHISTORY_SIZE];
  const char            delimiter[] = "_";
       /* The following match the organization of "input_file_luns" */

  typedef enum 
  {
    INDEX_LEADING_GRAN,
    INDEX_MIDDLE_GRAN,
    INDEX_TRAILING_GRAN,
    INDEX_GEO_FILE,
    INDEX_REFL_TABLES,
    INDEX_EMISS_TABLES,
    INDEX_QA_TABLES,
    NUM_INPUT_FILES
  } input_file_macros_t;

       /* Logical unit numbers for all input files. */

  int32 input_file_luns[NUM_INPUT_FILES] = 
  {
    LEADING_L1A_GRANULE,
    FIRST_L1A_GRANULE,
    TRAILING_L1A_GRANULE,
    GEOLOCATION_FILE,
    REFLECTIVE_TABLES_FILE,
    EMISSIVE_TABLES_FILE,
    QA_TABLES_FILE
  };

       /* The following are macros for the output files */


       /* The following are LUNs for the MCF files and the output granules */

  int32 mcf_file_luns[NUM_OUTPUT_FILES] = 
  {
    MCF_FILE_QKM,
    MCF_FILE_HKM,
    MCF_FILE_1KM,
    MCF_FILE_OBC
  };
  int32 output_file_luns[NUM_OUTPUT_FILES] = 
  {
    L1B_EV_250M_FILE,
    L1B_EV_500M_FILE,
    L1B_EV_1000M_FILE,
    L1B_OBC_FILE
  };

  char  input_file_UR_name[NUM_INPUT_FILES][PGSd_PC_UREF_LENGTH_MAX]; 
  char  *input_filename_p[NUM_INPUT_FILES];
  int32 out_sd_id[NUM_OUTPUT_FILES] = {FAIL, FAIL, FAIL, FAIL};
  char  L1B_modis_filenames[132];

  int   R = 0;   /* file_index, etc. */

  /*
   * Starting output file index; changes to INDEX_L1B_EV_1000M_FILE
   * when HiRes night mode is not written.
   */
  int16 start_output_file_index = INDEX_L1B_EV_250M_FILE;    

  char str[80];  /* parameter value string of additional attributes */

  float64  eastboundingcoordinate[4]  ;
  float64  westboundingcoordinate[4]  ;
  float64  northboundingcoordinate[4] ;
  float64  southboundingcoordinate[4] ;
  char     daynightflag[COMMON_TEXT_SIZE];
  char     equatorcrossingdate[MAX_DATE_TIME_SIZE];
  char     equatorcrossingtime[MAX_DATE_TIME_SIZE];

  float64  equatorcrossinglongitude;

  char     exclusiongringflag[COMMON_TEXT_SIZE];
  char     rangebeginningdate[MAX_DATE_TIME_SIZE];
  char     rangebeginningtime[MAX_DATE_TIME_SIZE];
  char     rangeendingdate[MAX_DATE_TIME_SIZE];
  char     rangeendingtime[MAX_DATE_TIME_SIZE];
  float64  gringpointlongitude[4];
  float64  gringpointlatitude[4];
  int32    gringpointsequenceno[4];
 
  int orbitnumber                     = 0; 
  
  const char* extension = "NRT.hdf";
  const char* product_version = "NRT.061";
  if(strcmp(runtime_params->ReprocessingActual, "Near Real Time"))
  {
    extension += 4;
    product_version += 4;

  }

  /* This is not necessary, but to supress 
     the warning messages (UMR) when run Purify
  */
  for (i = 0; i < 30; i++)
  {
    beginDateTimeA[i] = '\0';
    beginDateTimeB[i] = '\0';
    procDateTimeA[i] = '\0';
    procDateTimeB[i] = '\0';
  }
  
  in_fileID    = FIRST_L1A_GRANULE;
  out_sd_id[INDEX_L1B_EV_250M_FILE] = 
      L1B_Gran_Meta->L1B_Gran_sd_id[INDEX_L1B_250m];
  out_sd_id[INDEX_L1B_EV_500M_FILE] = 
      L1B_Gran_Meta->L1B_Gran_sd_id[INDEX_L1B_500m];
  out_sd_id[INDEX_L1B_EV_1000M_FILE] = 
      L1B_Gran_Meta->L1B_Gran_sd_id[INDEX_L1B_1km];

  /*--------------------------------------------
    Reopen L1B_OBC file to write metadata
  --------------------------------------------*/
  returnStatus = PGS_PC_GetReference (L1B_OBC_FILE, &Version, OBC_fname);
  if (returnStatus != PGS_S_SUCCESS)
  {
    returnStatus = MODIS_F_FILE_NOT_FOUND;
    L1BErrorMsg(location, returnStatus, 
        "Could not retrieve file name from PCF.",
                "PGS_PC_GetReference", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  out_sd_id[INDEX_L1B_OBC_FILE] = SDstart (OBC_fname, DFACC_RDWR); 
  if (out_sd_id[INDEX_L1B_OBC_FILE] == FAIL)
  {
    returnStatus = MODIS_F_FILE_NOT_OPENED;
    L1BErrorMsg(location, returnStatus, 
        "Could not open file for SD read/write access.",
                "SDstart", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  /*---------------------------------- 
    Generate processing time.
  ----------------------------------*/
  if((t = time(NULL)) == -1)
    SMF_ERROR(MODIS_W_TIME_INCORRECT,
        "time() failed in Write_Gran_Metadata(), Metadata.c");

  if(!(curtime = (struct tm *)gmtime(&t)))
    SMF_ERROR(MODIS_F_NOK, 
        "localtime() failed in Write_Gran_Metadata(), Metadata.c"); 

  /* obtain CCSDS ASCII time code A */
  if((strftime(procDateTimeA, TIMECODEASIZE, 
      "%Y-%m-%dT%H:%M:%S.000000Z", curtime)) == 0)
    SMF_ERROR(MODIS_W_TIME_INCORRECT, 
        "strftime() failed in Write_Gran_Metadata(), Metadata.c");
  
  /*---------------------------------------------------------------------
    Generate MODIS Output file names for each of 3 resolutions and OBC.
  ---------------------------------------------------------------------*/
  pgs_in.fileID = in_fileID;   
  pgs_in.version = 1;   
  pgs_in.hdfAttr = "CoreMetadata.0";
  
  get_string_attr( "RANGEBEGINNINGDATE", beginDateTimeA);
  beginDateTimeA[10] = 'T';
  get_string_attr( "RANGEBEGINNINGTIME", &beginDateTimeA[11]);
  beginDateTimeA[26] = 'Z';

  if (PGS_TD_ASCIItime_AtoB(beginDateTimeA, beginDateTimeB) != PGS_S_SUCCESS)
    SMF_ERROR(MODIS_F_HDF_ERROR, 
        "PGS_TD_ASCIItime_AtoB() in Write_Gran_Metadata(), Metadata.c");

  if(PGS_TD_ASCIItime_AtoB(procDateTimeA, procDateTimeB) != PGS_S_SUCCESS )
    SMF_ERROR(MODIS_F_HDF_ERROR, 
        "PGS_TD_ASCIItime_AtoB() in Write_Gran_Metadata(), Metadata.c");

  /* 
   * Generate LOCALGRANULEID
   */

  if (L1A_Gran->satellite_id == TERRA)
    strcpy(sub_shortname, "MOD02");
  else /* must be "Aqua", otherwise the code error out at Read_QA_Tables. */
    strcpy(sub_shortname, "MYD02");

  /* initialize the doi metadata */  
  strcpy(L1B_Gran_Meta->doi_attr_name, "identifier_product_doi");
  strcpy(L1B_Gran_Meta->doi_authority_attr_name, "identifier_product_doi_authority");
  strcpy(L1B_Gran_Meta->doi_authority_para_value, "http://dx.doi.org");

  /*------------------------------------------------
    Get input & ancil file names from pcf file
  ------------------------------------------------*/
  for ( R = 0 ; R < NUM_INPUT_FILES ; R++ )
  {
    if (R == INDEX_LEADING_GRAN && 
        QA->QA_common.missing_leading_granule == True)
      ;
    else if (R == INDEX_TRAILING_GRAN && 
      QA->QA_common.missing_trailing_granule == True)
      ;
    else
    {
      version = 1;
      returnStatus = PGS_PC_GetUniversalRef(input_file_luns[R], &version,
                                            &input_file_UR_name[R][0]);
      if (returnStatus != PGS_S_SUCCESS)
      {
        returnStatus = MODIS_F_FILE_NOT_FOUND;
        L1BErrorMsg(location, returnStatus, 
            "Could not get universal reference from PCF.",
                    "PGS_PC_GetUniversalRef", input_file_luns[R], NULL, True);
        return returnStatus;
      }
    }
  }

  /*
   * Assemble the INPUTPOINTER from the universal references of the input files
   * except that the geolocation UR is placed in ANCILLARYINPUTPOINTER.
   * The logic implemented also excludes any MOD01 granule that is treated as
   * "missing".
   */

  i = 0;
  for ( R = 0; R < NUM_INPUT_FILES ; R++ )
  {
    if ( R == INDEX_GEO_FILE)
      continue;  /* Geolocation file name is put in ancillaryinputpointer */
    if (R == INDEX_LEADING_GRAN && 
        QA->QA_common.missing_leading_granule == True)
      continue;
    else if (R == INDEX_TRAILING_GRAN && 
      QA->QA_common.missing_trailing_granule == True)
      continue;

    input_filename_p[i] = &input_file_UR_name[R][0];
    for (ptr = &input_file_UR_name[R][0]; *ptr != '\0'; ptr++)
      if (*ptr == '/')                    /* skip over the path name */
        input_filename_p[i] = ptr+1;
    i++;
  }
  input_filename_p[i] = NULL; /* Indicates end of array to PGS_MET_SetAttr */

  /*------------------------------------------------
     Initialize metadata structures 
  ------------------------------------------------*/
  /*
   * If this is a NIGHT mode granule and the runtime parameter
   * Write_Night_Mode_HiRes_Data is False, do not write the 250m and 500m
   * resolution files, and so start with output file number 2, not 0.
   */

  if (skip_night_hi_res == True)
    start_output_file_index = INDEX_L1B_EV_1000M_FILE;
  else
    start_output_file_index = INDEX_L1B_EV_250M_FILE;


  for (R = start_output_file_index; R < NUM_OUTPUT_FILES; R++)
  {
    returnStatus = PGS_MET_Init(mcf_file_luns[R], mdhandles);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_FILE_NOT_FOUND;
      L1BErrorMsg(location, returnStatus,
                  "Could not initialize the PGS Metadata from MCF file.",
                  "PGS_MET_Init", mcf_file_luns[R], NULL, True);
      return returnStatus;
    }

    /* Check if the "ShortName" value agrees with the satellite id. */

    ptr = shortname;
    returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "ShortName", &ptr);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not get the value of \"ShortName\" from MCF file.",
                  "PGS_MET_GetSetAttr", mcf_file_luns[R], NULL, True);
      return returnStatus;
    }
    if (!((L1A_Gran->satellite_id == TERRA && 
        strncmp(shortname, "MOD", 3) == 0) ||
            (L1A_Gran->satellite_id == AQUA && 
        strncmp(shortname, "MYD", 3) == 0))) 
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus,
                  "ShortName in MCF file is not consistent with the"
                  " instrument name in L1A granule.",
                  NULL, mcf_file_luns[R], 
                  "The MCF file could be invalid.", True);
      return returnStatus;
    }

    /* generate the doi parameter value for each L1B product */
    ptr = &L1B_Gran_Meta->doi_para_value[R][0];
    sprintf(ptr, "10.5067/MODIS/%s.%s", shortname, product_version);

    /*
     * Read the VersionID from MCF and create the LocalGranuleID
     */

    returnStatus = PGS_MET_GetSetAttr(mdhandles[1], "VersionID", &VersionID);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_READ_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not get the value of \"VersionID\" from MCF file.",
                  "PGS_MET_GetSetAttr", mcf_file_luns[R], NULL, True);
      return returnStatus;
    }
    if (VersionID < 1 || VersionID > 999)
    {
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus,
                  "VersionID is not in range of [1-999].",
                  NULL, mcf_file_luns[R], 
                  "The MCF file could be invalid.", True);
      return returnStatus;
    }

#ifdef SEADAS
    Version = 1;
    returnStatus = PGS_PC_GetReference (output_file_luns[R], &Version,
                                        L1B_modis_filenames);
    sprintf(L1B_modis_filenames, "%s", basename(L1B_modis_filenames));
#else
    sprintf(L1B_modis_filenames,
            "%s%s.A%4.4s%3.3s.%2.2s%2.2s.%03d.%4.4s%3.3s%2.2s%2.2s%2.2s.%s",
            sub_shortname, file_sub_name[R],
            &beginDateTimeB[0], &beginDateTimeB[5],
            &beginDateTimeB[9],&beginDateTimeB[12],
            VersionID, &procDateTimeB[0],
            &procDateTimeB[5],&procDateTimeB[9],&procDateTimeB[12],
            &procDateTimeB[15], extension);
#endif

  /*------------------------------------------------------
     Set all "coremetadata.0" attributes not 
     initialized via PGS_MET_Init (not PCF/MCF).
  ------------------------------------------------------*/
  pgs_out_mdHandle = mdhandles[1];  /* [1] is for inventory (Core) metadata */

  pgs_in.fileID = in_fileID;   
  pgs_in.version = 1;   
  pgs_in.hdfAttr = "CoreMetadata.0";

  set_string_attr   ( "PGEVERSION", runtime_params->PGE02_Version);

  copy_string_attr  ( "RANGEBEGINNINGDATE"   , rangebeginningdate );
  copy_string_attr  ( "RANGEBEGINNINGTIME"   , rangebeginningtime );
  copy_string_attr  ( "RANGEENDINGDATE"      , rangeendingdate );
  copy_string_attr  ( "RANGEENDINGTIME"      , rangeendingtime );

  copy_string_attr  ( "DAYNIGHTFLAG"         , daynightflag );
  
  set_string_attr ( "ADDITIONALATTRIBUTENAME.1"  , 
                    "AveragedBlackBodyTemperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_bb_gran);
  set_string_attr ( "PARAMETERVALUE.1"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.2"  , 
                    "AveragedMirrorTemperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_mir_gran);
  set_string_attr ( "PARAMETERVALUE.2"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.3"  , 
                    "AveragedFocalPlane1Temperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_fp1);
  set_string_attr ( "PARAMETERVALUE.3"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.4"  , 
                    "AveragedFocalPlane2Temperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_fp2);
  set_string_attr ( "PARAMETERVALUE.4"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.5"  , 
                    "AveragedFocalPlane3Temperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_fp3);
  set_string_attr ( "PARAMETERVALUE.5"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.6"  , 
                    "AveragedFocalPlane4Temperature" );
  sprintf(str, "%7.2f", PP->PP_Emiss.T_fp4);
  set_string_attr ( "PARAMETERVALUE.6"           , str );
  set_string_attr ( "ADDITIONALATTRIBUTENAME.7"  , 
                    "CalibrationQuality" );

  set_string_attr ( "PARAMETERVALUE.7"           , 
                    CALIBRATIONQUALITY_MACRO);

  set_string_attr ( "ADDITIONALATTRIBUTENAME.8"  , 
                    "MissionPhase" );

  set_string_attr ( "PARAMETERVALUE.8"           , 
                    tables->QA.common_QA_tables.mission_phase );

  set_string_attr ( "ADDITIONALATTRIBUTENAME.9"  , "NadirPointing" );

  set_string_attr ( "PARAMETERVALUE.9"           , NADIRPOINTING_MACRO);

  set_string_attr ( "ADDITIONALATTRIBUTENAME.10" , L1B_Gran_Meta->doi_attr_name );

  set_string_attr ( "PARAMETERVALUE.10"          , &L1B_Gran_Meta->doi_para_value[R][0] );

  set_string_attr ( "ADDITIONALATTRIBUTENAME.11" , L1B_Gran_Meta->doi_authority_attr_name );

  set_string_attr ( "PARAMETERVALUE.11"          , L1B_Gran_Meta->doi_authority_para_value );

  switch(R)
  {
    case INDEX_L1B_EV_250M_FILE:
      set_attr        ( "QAPERCENTMISSINGDATA.1"     , 
                        &L1B_Gran_Meta->qapercent_missing_250m);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.1"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_250m);
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.1"   , 
                        &L1B_Gran_Meta->qapercent_outofbound_250m); 
      set_string_attr ( "PARAMETERNAME.1"            , "EV_250_RefSB" );
      set_string_attr   ( "AUTOMATICQUALITYFLAG.1", 
                        AUTOMATICQUALITYFLAG_MACRO);  
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.1", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO); 
      break;   
    case INDEX_L1B_EV_500M_FILE:
      set_attr        ( "QAPERCENTMISSINGDATA.1"     , 
                        &L1B_Gran_Meta->qapercent_missing_500m);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.1"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_500m);
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.1"   , 
                        &L1B_Gran_Meta->qapercent_outofbound_500m); 
      set_string_attr ( "PARAMETERNAME.1"            , "EV_500_RefSB" );
      set_string_attr   ( "AUTOMATICQUALITYFLAG.1", 
                        AUTOMATICQUALITYFLAG_MACRO);  
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.1", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO); 
      break;
    case INDEX_L1B_EV_1000M_FILE:
      set_attr        ( "QAPERCENTMISSINGDATA.1"     , 
                        &L1B_Gran_Meta->refl_1km_qapercent_missing);
      set_attr        ( "QAPERCENTMISSINGDATA.2"     , 
                        &L1B_Gran_Meta->emiss_qapercent_missing);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.1"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_refl_1km);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.2"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_emiss);
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.1"   , 
                        &L1B_Gran_Meta->refl_1km_qapercent_outofbound); 
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.2"   , 
                        &L1B_Gran_Meta->emiss_qapercent_outofbound);
      set_string_attr ( "PARAMETERNAME.1"            , "EV_1KM_RefSB" );
      set_string_attr ( "PARAMETERNAME.2"            , "EV_1KM_Emissive" );
      set_string_attr   ( "AUTOMATICQUALITYFLAG.1", 
                        AUTOMATICQUALITYFLAG_MACRO);  
      set_string_attr   ( "AUTOMATICQUALITYFLAG.2", 
                        AUTOMATICQUALITYFLAG_MACRO);     
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.1", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO); 
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.2", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO);   
      break;
    default: 
      set_attr        ( "QAPERCENTMISSINGDATA.1"     , 
                        &L1B_Gran_Meta->refl_1km_qapercent_missing);
      set_attr        ( "QAPERCENTMISSINGDATA.2"     , 
                        &L1B_Gran_Meta->emiss_qapercent_missing);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.1"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_refl_1km);
      set_attr        ( "QAPERCENTINTERPOLATEDDATA.2"  , 
                        &L1B_Gran_Meta->qapercent_interpolated_emiss);
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.1"   , 
                        &L1B_Gran_Meta->refl_1km_qapercent_outofbound); 
      set_attr        ( "QAPERCENTOUTOFBOUNDSDATA.2"   , 
                        &L1B_Gran_Meta->emiss_qapercent_outofbound);
      set_string_attr ( "PARAMETERNAME.1"            , "EV_1KM_RefSB" );
      set_string_attr ( "PARAMETERNAME.2"            , "EV_1KM_Emissive" );
      set_string_attr   ( "AUTOMATICQUALITYFLAG.1", 
                        AUTOMATICQUALITYFLAG_MACRO);  
      set_string_attr   ( "AUTOMATICQUALITYFLAG.2", 
                        AUTOMATICQUALITYFLAG_MACRO);     
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.1", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO); 
      set_string_attr   ( "AUTOMATICQUALITYFLAGEXPLANATION.2", 
                        AUTOMATICQUALITYFLAGEXPLANATION_MACRO); 
      break; 
 }  

  set_string_attr ( "REPROCESSINGPLANNED", 
                        runtime_params->ReprocessingPlanned);
  set_string_attr ( "REPROCESSINGACTUAL",  
                        runtime_params->ReprocessingActual);

  set_ptrstring_attr ( "INPUTPOINTER"            , input_filename_p);
  
  set_string_attr ( "ANCILLARYINPUTTYPE.1", 
                        ANCILLARYINPUTTYPE_MACRO);
  set_string_attr ( "ANCILLARYINPUTPOINTER.1", 
                        &input_file_UR_name[INDEX_GEO_FILE][0]);
  set_string_attr ( "LOCALGRANULEID"             , 
                        L1B_modis_filenames); 

  /*
   * Generate the new "ProductionHistory" based on the MOD01 and MOD03
   * PGEVERSIONS.  Using the PGEVERSION ensures backwards compatibility with
   * earlier version of MOD01 and MOD03 which did not write the
   * "ProductionHistory" metadata variable.
   */

  /* Get PGEVERSION from L1A in granule MOD01 */
  get_string_attr   ( "PGEVERSION", MOD01_ProductionHistory);

  /* Make the first part of the PGE02 ProductionHistory: */
  sprintf(ProductionHistory, "PGE02:");
  strcpy (MCSTversion, tables->QA.common_QA_tables.MCST_Version);
  safe_strcat( ProductionHistory,
               strtok(MCSTversion, delimiter),
               /* Eliminate the "_Aqua/Terra" designator on the LUT version.*/
               MAX_PRODUCTIONHISTORY_SIZE );
  safe_strcat(ProductionHistory, ";PGE01:", MAX_PRODUCTIONHISTORY_SIZE);

  /*---------------------------------------------
    Setup Archive metadata fields
  ---------------------------------------------*/

  pgs_in.fileID = GEOLOCATION_FILE;
  pgs_in.version = 1;   
  pgs_in.hdfAttr = "CoreMetadata.0";

  /* Get PGEVERSION from geolocation file MOD03 */
  get_string_attr   ( "PGEVERSION", MOD03_ProductionHistory);

  /* Add on the second part of the PGE02 ProductionHistory: */
  if ( strcmp (MOD03_ProductionHistory, MOD01_ProductionHistory) ) 
     /* MOD01 and MOD03 from different PGE01 versions */
  {
    safe_strcat(ProductionHistory, MOD01_ProductionHistory,
        MAX_PRODUCTIONHISTORY_SIZE);
    safe_strcat(ProductionHistory, ":mod01", MAX_PRODUCTIONHISTORY_SIZE);
    safe_strcat(ProductionHistory, ";PGE01:", MAX_PRODUCTIONHISTORY_SIZE);
    safe_strcat(ProductionHistory, MOD03_ProductionHistory,
        MAX_PRODUCTIONHISTORY_SIZE);
    safe_strcat(ProductionHistory, ":mod03", MAX_PRODUCTIONHISTORY_SIZE);
  }
  else 
     /* MOD01 and MOD03 from same PGE01 version */
  {
    safe_strcat(ProductionHistory, MOD01_ProductionHistory,
        MAX_PRODUCTIONHISTORY_SIZE);
  }

  copy_attr         ( "ORBITNUMBER.1"        , &orbitnumber ); 
  copy_attr         ( "EQUATORCROSSINGLONGITUDE.1", &equatorcrossinglongitude);

  copy_string_attr  ( "EQUATORCROSSINGDATE.1", equatorcrossingdate);
  copy_string_attr  ( "EQUATORCROSSINGTIME.1", equatorcrossingtime);
  copy_string_attr  ( "EXCLUSIONGRINGFLAG.1" , exclusiongringflag);
  
  copy_attr         ( "GRINGPOINTLATITUDE.1" ,  gringpointlatitude);
  copy_attr         ( "GRINGPOINTLONGITUDE.1" , gringpointlongitude);
  copy_attr         ( "GRINGPOINTSEQUENCENO.1", gringpointsequenceno);


  pgs_out_mdHandle = mdhandles[2];  /* [2] is for Archive metadata */

  /* Write PGE02 ProductionHistory into Archive metadata */
  set_string_attr   ( "PRODUCTIONHISTORY", ProductionHistory);

  pgs_in.hdfAttr = "ArchiveMetadata.0";
  copy_attr ( "EASTBOUNDINGCOORDINATE"  ,  eastboundingcoordinate);
  copy_attr ( "WESTBOUNDINGCOORDINATE"  ,  westboundingcoordinate);
  copy_attr ( "NORTHBOUNDINGCOORDINATE"  ,  northboundingcoordinate);
  copy_attr ( "SOUTHBOUNDINGCOORDINATE"  ,  southboundingcoordinate);

  set_string_attr   ( "ALGORITHMPACKAGENAME" , 
                        ALGORITHMPACKAGENAME_MACRO); 
  set_string_attr   ( "ALGORITHMPACKAGEACCEPTANCEDATE",
            tables->QA.common_QA_tables.AlgorithmPackageAcceptanceDate);
  set_string_attr   ( "ALGORITHMPACKAGEMATURITYCODE",
            tables->QA.common_QA_tables.AlgorithmPackageMaturityCode);
  set_string_attr   ( "ALGORITHMPACKAGEVERSION", 
                        tables->QA.common_QA_tables.MCST_Version);
  set_string_attr   ( "INSTRUMENTNAME", INSTRUMENTNAME_MACRO);
  set_string_attr   ( "PROCESSINGCENTER",
                        runtime_params->ProcessingCenter);
  set_string_attr   ( "PROCESSINGENVIRONMENT",
                        runtime_params->ProcessingEnvironment);

  /*-----------------------------------------------------
    Write the core and archieve metadata to hdf files 
  -----------------------------------------------------*/
    returnStatus = PGS_MET_Write(mdhandles[1], MECS_CORE, out_sd_id[R]);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not write ECS Core metadata to output granule.",
                  "PGS_MET_Write", output_file_luns[R], NULL, True);
      return returnStatus;
    }

    returnStatus = PGS_MET_Write(mdhandles[2], MECS_ARCH, out_sd_id[R]);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_WRITE_ERROR;
      L1BErrorMsg(location, returnStatus,
                  "Could not write ECS Archive metadata to output granule.",
                  "PGS_MET_Write", output_file_luns[R], NULL, True);
      return returnStatus;
    }

    PGS_MET_Remove();
  }
  
  returnStatus = Write_Global_Metadata(L1B_Gran_Meta, 
                                       QA, 
                                       tables,
                                       out_sd_id[INDEX_L1B_OBC_FILE],
                                       skip_night_hi_res);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL,
                "Write_Global_Metadata", 0, NULL, True);
    return returnStatus;
  }

 
  /*-------------------------
    Close L1B_OBC file
  -------------------------*/
  hdf_return = SDend(out_sd_id[INDEX_L1B_OBC_FILE]);
  if (hdf_return == FAIL)
  {
    returnStatus = MODIS_F_HDF_ERROR;
    L1BErrorMsg(location, returnStatus, NULL,
                "SDend", L1B_OBC_FILE, NULL, True);
    return returnStatus;
  }

  return(MODIS_S_OK);
}

void  get_attr( char *field, void *value)
/*
!C****************************************************************
!Description:
   Helper function to Write_Gran_Metadata().  Reads a non-string data item
   from ECS metadata.

!Input Parameters: 
   char        *field        name of data to read

!Output Parameters:
   char        *value        buffer to read data into

!Revision History:
 (continue at top of the file)

   Revision 01.01 November 24, 1999
   Error message revised.
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
   pgs_in is a global variable which holds the file LUN and other info
   needed for PGS_MET_GetPCAttr.

!END********************************************************************
*/
{  
  int ret;
  ret = PGS_MET_GetPCAttr(pgs_in.fileID, 
                          pgs_in.version, 
                          pgs_in.hdfAttr, 
                          field, 
                          value);
  if (ret != PGS_S_SUCCESS)
  {
    char errorMessage[256];
    sprintf(errorMessage, "Could not get metadata \"%s\" from %s.", 
            field, pgs_in.hdfAttr);
    L1BErrorMsg("get_attr", MODIS_F_READ_ERROR, errorMessage,
                "PGS_MET_GetPCAttr", pgs_in.fileID, NULL, True);
  }
}

void  get_string_attr( char *field, char *value)
/*
!C****************************************************************
!Description:
   Helper function to Write_Gran_Metadata().  Reads a string data item
   from ECS metadata.

!Input Parameters: 
   char        *field        name of data to read

!Output Parameters:
   char        *value        buffer to read data into

!Revision History:
 (continue at top of the file)

   Revision 01.01 November 24, 1999
   Error message revised.
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
   pgs_in is a global variable which holds the file LUN and other info
   needed for PGS_MET_GetPCAttr.

!END********************************************************************
*/
{
  int ret;
  ret = PGS_MET_GetPCAttr(pgs_in.fileID, 
                          pgs_in.version, 
                          pgs_in.hdfAttr, 
                          field, 
                          &value);
  if (ret != PGS_S_SUCCESS)
  {
    char errorMessage[256];
    sprintf(errorMessage, "Could not get metadata \"%s\" from %s.", 
            field, pgs_in.hdfAttr);
    L1BErrorMsg("get_string_attr", MODIS_F_READ_ERROR, errorMessage,
                "PGS_MET_GetPCAttr", pgs_in.fileID, NULL, True);
  }
}

void set_attr( char *field, void *value)
/*
!C*************************************************************************
!Description:
   Helper function to Write_Gran_Metadata().  Writes data (non-string)
   into the appropriate ECS metadata buffer identified by pgs_out_mdHandle.

!Input Parameters: 
   char        *field         Name of field to write.
   int         value[]        Pointer to value(s) to be written

!Output Parameters:
   none

!Revision History:
 (continue at top of the file)

   Revision 01.01 November 24, 1999
   Error message revised.
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
   pgs_out_mdHandle is a global variable identifying the buffer to write
   the data into (either ECS core metadata or ECS archive metadata).

!END************************************************************************
*/
{
  int ret;
  ret = PGS_MET_SetAttr(pgs_out_mdHandle, field, value);
  if (ret != PGS_S_SUCCESS)
  {
    char errorMessage[256];
    sprintf(errorMessage, "Toolkit error setting metadata \"%s\".", field);
    L1BErrorMsg("set_attr", MODIS_F_NOK, errorMessage,
                "PGS_MET_SetAttr", 0, NULL, True);
  }   
}

void set_string_attr( char *field, char *value)
/*
!C**************************************************************************
!Design Notes:
  
!Description:
   Helper function to Write_Gran_Metadata().  Writes data (string)
   into the appropriate ECS metadata buffer identified by pgs_out_mdHandle.

!Input Parameters: 
   char   *field         Name of field to write.
   char   *value         Pointer to string value to be written

!Output Parameters:
   none

!Revision History:
 (continue at top of the file)

   Revision 01.01 November 24, 1999
   Error message revised.
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
   pgs_out_mdHandle is a global variable identifying the buffer to write
   the data into (either ECS core metadata or ECS archive metadata).

!END********************************************************************
*/
{
  int ret;
  ret = PGS_MET_SetAttr(pgs_out_mdHandle, field, &value);
  if (ret != PGS_S_SUCCESS)
  {
    char errorMessage[256];
    sprintf(errorMessage, "Toolkit error setting metadata \"%s\".", field);
    L1BErrorMsg("set_string_attr", MODIS_F_NOK, errorMessage,
                "PGS_MET_SetAttr", 0, NULL, True);
  }
}


void set_ptrstring_attr( char *field, char **value )
/*
!C*************************************************************************
!Description:
   Helper function to Write_Gran_Metadata().  Writes data (array of strings)
   into the appropriate ECS metadata buffer identified by pgs_out_mdHandle.

!Input Parameters: 
   char  *field         Name of field to write.
   char  **value        Pointer to array of strings to be written

!Output Parameters:
   none

!Revision History:
 (continue at top of the file)

   Revision 01.01 November 24, 1999
   Error message revised.
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
   pgs_out_mdHandle is a global variable identifying the buffer to write
   the data into (either ECS core metadata or ECS archive metadata).

!END********************************************************************
*/
{
  int ret;
  ret = PGS_MET_SetAttr(pgs_out_mdHandle, field, value);
  if (ret != PGS_S_SUCCESS)
  {
    char errorMessage[256];
    sprintf(errorMessage, "Toolkit error setting metadata \"%s\".", field);
    L1BErrorMsg("set_ptrstring_attr", MODIS_F_NOK, errorMessage,
                "PGS_MET_SetAttr", 0, NULL, True);
  }
}

void copy_attr( char *field, void *value)
/*
!C****************************************************************
!Description:   Helper function to Write_Gran_Metadata()

!Input Parameters: 
     pgs_meta_t  pgs_in
     char        *field
     void        *value

!Output Parameters:
     none

!Revision History:
 (continue at top of the file)
 
 Revision 01.00 1998/4/8 
 Original Development
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
  
!Team-unique Header:
  
!References and Credits:
    This software is developed by the MODIS Characterization Support
    Team (MCST)for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.
 
    The time setting code was borrowed from Paul Fishers timea.c

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.
 
!Design Notes:
  
!END********************************************************************
*/
{
  get_attr( field, value);
  set_attr( field, value);
}

void copy_string_attr( char *field, char *value)
/*
!C****************************************************************
!Description:   Helper function to Write_Gran_Metadata()

!Input Parameters: 
     pgs_meta_t  pgs_in
     char        *field
     int         value[]

!Output Parameters:
     none

!Revision History:
 (continue at top of the file)

 Revision 01.00 1998/4/8 
 Original Development
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)
  
!Team-unique Header:
  
!References and Credits:
    This software is developed by the MODIS Characterization Support
    Team (MCST)for the National Aeronautics and Space Administration,
    Goddard Space Flight Center, under contract NAS5-32373.
 
    The time setting code was borrowed from Paul Fishers timea.c

    HDF portions developed at the National Center for Supercomputing
    Applications at the University of Illinois at Urbana-Champaign.
 
!Design Notes:
  
!END********************************************************************
*/
{
  get_string_attr( field, value);
  set_string_attr( field, value);
}

PGSt_SMF_status Gran_Meta_Cal(L1A_granule_t       *L1A_Gran,
                              L1B_granule_t       *L1B_Gran, 
                              Preprocess_Data_t   *PP,
                              QA_Data_t           *QA, 
                              L1B_Scan_Metadata_t *L1B_Scan_Meta,
                              L1B_Gran_Metadata_t *L1B_Gran_Meta)
/*
!C****************************************************************
!Description:   pass the values from L1A_Gran, L1B_Gran, PP, tables, QA, to 
                L1B_Gran_Meta and calculate some granule metadata values.

!Input Parameters: 
     L1A_granule_t *L1A_Gran   contains information about middle L1A granule
     L1B_granule_t *L1B_Gran   contains L1B granule level information
     Preprocess_Data_t *PP     contains focal plane set point
     QA_Data_t *QA             contains QA data
     L1B_Gran_Metadata_t *L1B_Gran_Meta contains bit QA flag

!Output Parameters:
     L1B_Gran_Metadata_t *L1B_Gran_Meta  contains granule level metadata

!Revision History:
 (continue at top of the file)

   Revision 01.06  June 23, 2003, Razor Issue #192
    Changed the telemetry mnemonic which determines whether the nadir aperture
    door (NAD) is open from "CR_DR_NAD_OPEN" to the  equivalent "CR_DR_NAD_CLSD"
    mnemonic and altered the L1B code logic to correctly use the new value.  It
    was discovered that at two different times in the history of the MODIS/Terra
    instrument the mnemonic "CR_DR_NAD_OPEN" did not correctly reflect the state
    of the NAD, whereas the mnemonic "CR_DR_NAD_CLSD" has been consistently
    reliable.  
   Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.05 March 27, 2003, Razor Issue #173
   Enclosed the rows in the initializer of array-of-struct "temp" with braces
   for ANSI-C compliance.
   Liqin Tan, SAIC GSO (ltan@saicmodis.com)

   Revision 01.04  April 25, 2002  Razor Issue #183
   Changed initialization of structure "temp" to match definition.
   Gwyn Fireman, SAIC GSO (fireman@mcst.gsfc.nasa.gov)

   Revision 01.03 September 21, 2000
   Corrected calculation of % RSB in night mode, razor issue 137.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   ... (many changes not logged) ...
 
 Revision 01.02 August 28, 1999
 Added checking if the data read from L1A are valid.
 Zhenying Gu (zgu@mcst.gsfc.nasa.gov) and Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.01 Feb 22, 1999
 Change dthe calculation of num_day_scans to a simple assignment from
 L1A_Gran->num_day_scans.
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
  char *location = "Gran_Meta_Cal";
  int32             B, B_38, D, S;
  int32             missing_pixels, outofbound_pixels, 
                    interpolated_pixels, total_pixels;

    /* The following three variables are used to calculate percentages
     * of some low quality data.
     */

  int32 num_pixels_per_250m_det;   
  int32 num_pixels_per_500m_det;
  int32 num_pixels_per_1km_det;

  enum {
    NAD_Door,
    SVD_Door,
    SDD_Door,
    SDD_Screen,
    NUM_Config
  }Door_Screen_Config;
  
  struct {
    char *vname;     /*vdata_name*/
    char *fname;     /*field_name*/
    uint16 buffer;  
  } temp[NUM_Config] = {
   {"Telemetry Major Cycle 3A of 7",      "CR_DR_NAD_CLSD",  1},
   {"Telemetry Major Cycle 3A of 7",      "CR_DR_SVD_OPEN",  0},
   {"Telemetry Major Cycle 3A of 7",      "CR_DR_SDD_OPEN",  0},
   {"Telemetry Major Cycle 3A of 7",      "CR_DR_SDS_OPEN",  0} 
  };
  int32 day_bands_at_night;
  uint32 day_bands_at_night_mask;

    /* Calculate the percentages of different kinds of low quality data
     * for each detector.
     */
 
   
  num_pixels_per_250m_det = L1A_Gran->num_scans * EV_250m_FRAMES;
  num_pixels_per_500m_det = L1A_Gran->num_scans * EV_500m_FRAMES;
  num_pixels_per_1km_det  = L1A_Gran->num_scans * EV_1km_FRAMES;
 
  L1B_Gran_Meta->missAllScanDataPercent = 
      100*(float32)QA->QA_common.num_missing_scans/L1A_Gran->num_scans;
 
  for (D = 0; D < NUM_250M_BANDS*DETECTORS_PER_250M_BAND; D++)
  {
    L1B_Gran_Meta->missInScanDataPercent[D] =
         100*(float32)QA->QA_common.num_missing_data_in_scans[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->nightRSBPercent[D] = 100. * (float32)
                     QA->QA_common.num_rsb_at_night_scans / 
                                    L1A_Gran->num_scans;
    L1B_Gran_Meta->deadDetectorDataPercent[D] =
         100*(float32)QA->QA_common.num_dead_detector_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->deadSubframeDataPercent[D] =
         100*(float32)QA->QA_common.num_dead_subframe_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->sectorRotateDataPercent[D] =
         100*(float32)QA->QA_common.num_sector_rotation_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->saturatedDataPercent[D] =
         100*(float32)QA->QA_common.num_saturated_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->noBGDataPercent[D] =
         100*(float32)QA->QA_common.num_no_bg_DN_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->badDNStarStarRSBDataPercent[D] =
         100*(float32)QA->QA_common.num_bad_dn_star_star_RSB_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->exceedMaxForScalingPercent[D] =
         100*(float32)QA->QA_common.num_exceed_max_for_scaling[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->NADClosedDataPercent[D] =
         100*(float32)QA->QA_common.num_nadir_door_closed_EV_data[D]/
                                    num_pixels_per_250m_det;
    L1B_Gran_Meta->moonInSVPTEBDataPercent[D] = 0; 
  }

  for (; D < NUM_250M_BANDS*DETECTORS_PER_250M_BAND + 
                  NUM_500M_BANDS*DETECTORS_PER_500M_BAND; D++)
  {
    L1B_Gran_Meta->missInScanDataPercent[D] =
         100*(float32)QA->QA_common.num_missing_data_in_scans[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->nightRSBPercent[D] = 100. * (float32)
                     QA->QA_common.num_rsb_at_night_scans / 
                                    L1A_Gran->num_scans;
    L1B_Gran_Meta->deadDetectorDataPercent[D] =
         100*(float32)QA->QA_common.num_dead_detector_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->deadSubframeDataPercent[D] =
         100*(float32)QA->QA_common.num_dead_subframe_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->sectorRotateDataPercent[D] =
         100*(float32)QA->QA_common.num_sector_rotation_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->saturatedDataPercent[D] =
         100*(float32)QA->QA_common.num_saturated_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->noBGDataPercent[D] =
         100*(float32)QA->QA_common.num_no_bg_DN_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->badDNStarStarRSBDataPercent[D] =
         100*(float32)QA->QA_common.num_bad_dn_star_star_RSB_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->exceedMaxForScalingPercent[D] =
         100*(float32)QA->QA_common.num_exceed_max_for_scaling[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->NADClosedDataPercent[D] =
         100*(float32)QA->QA_common.num_nadir_door_closed_EV_data[D]/
                                    num_pixels_per_500m_det;
    L1B_Gran_Meta->moonInSVPTEBDataPercent[D] = 0; 
  }

  for (; D < NUM_DETECTORS; D++)
  {
    L1B_Gran_Meta->missInScanDataPercent[D] =
         100*(float32)QA->QA_common.num_missing_data_in_scans[D]/
                                    num_pixels_per_1km_det;
 
    if (D < NUM_REFLECTIVE_DETECTORS - 10 || 
                   (D >= (NUM_REFLECTIVE_DETECTORS - 10)
      + BAND26 * 10 && D < NUM_REFLECTIVE_DETECTORS + BAND26 * 10))
      L1B_Gran_Meta->nightRSBPercent[D] = 100. * (float32)
                     QA->QA_common.num_rsb_at_night_scans / 
                                    L1A_Gran->num_scans;
    else
      L1B_Gran_Meta->nightRSBPercent[D] = 0;

/************************* Begin Band 26 Section **************************/
#ifdef WRITE_BAND_26_SDS
    if (D >= (NUM_REFLECTIVE_DETECTORS - 10) + BAND26 * 10 && 
        D < NUM_REFLECTIVE_DETECTORS + BAND26 * 10)
      L1B_Gran_Meta->nightRSBPercent[D] = 0;
#endif /* WRITE_BAND_26_SDS */
/************************** End Band 26 Section ***************************/

    L1B_Gran_Meta->deadDetectorDataPercent[D] =
         100*(float32)QA->QA_common.num_dead_detector_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->sectorRotateDataPercent[D] =
         100*(float32)QA->QA_common.num_sector_rotation_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->saturatedDataPercent[D] =
         100*(float32)QA->QA_common.num_saturated_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->noBGDataPercent[D] =
         100*(float32)QA->QA_common.num_no_bg_DN_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->badDNStarStarRSBDataPercent[D] =
         100*(float32)QA->QA_common.num_bad_dn_star_star_RSB_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->exceedMaxForScalingPercent[D] =
         100*(float32)QA->QA_common.num_exceed_max_for_scaling[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->NADClosedDataPercent[D] =
         100*(float32)QA->QA_common.num_nadir_door_closed_EV_data[D]/
                                    num_pixels_per_1km_det;
    L1B_Gran_Meta->moonInSVPTEBDataPercent[D] =
         100*(float32)QA->QA_common.num_moon_in_SVP_TEB_EV_data[D]/
                                    num_pixels_per_1km_det;
  }

  for (D = 0; D < NUM_DETECTORS; D++){
    L1B_Gran_Meta->uncalibratedDataPercent[D] =
                        L1B_Gran_Meta->missAllScanDataPercent         +
                        L1B_Gran_Meta->nightRSBPercent[D]             +
                        L1B_Gran_Meta->missInScanDataPercent[D]       +
                        L1B_Gran_Meta->deadDetectorDataPercent[D]     +
                        L1B_Gran_Meta->sectorRotateDataPercent[D]     +
                        L1B_Gran_Meta->saturatedDataPercent[D]        +
                        L1B_Gran_Meta->noBGDataPercent[D]             +
                        L1B_Gran_Meta->badDNStarStarRSBDataPercent[D] +
                        L1B_Gran_Meta->exceedMaxForScalingPercent[D]  +
                        L1B_Gran_Meta->NADClosedDataPercent[D]        +
                        L1B_Gran_Meta->moonInSVPTEBDataPercent[D];
    if(D < NUM_HIGH_RESOLUTION_DETECTORS) 
      L1B_Gran_Meta->uncalibratedDataPercent[D] +=
                        L1B_Gran_Meta->deadSubframeDataPercent[D] ; 
  } 

    /* Calculate last value of bit QA flag and if there is any change for each bit */

  QA->QA_common.bit_QA_flags_change = 0;
  QA->QA_common.bit_QA_flags_last_value = 
      L1B_Scan_Meta->Bit_QA_Flags[L1A_Gran->num_scans - 1];
  for (S = 0; S < L1A_Gran->num_scans - 1; S++)
  {
    QA->QA_common.bit_QA_flags_change = QA->QA_common.bit_QA_flags_change|
                                        (L1B_Scan_Meta->Bit_QA_Flags[S+1]^
                                         L1B_Scan_Meta->Bit_QA_Flags[S]);
  }

    /******** Special case for Day bands at Night bit (bit 11) **********
     * For bit 11 in the bit_QA_flags_change and bit_QA_flags_last_value,
     * the logic is different than that defined above.  If there are more
     * than 15 scans in the Bit QA Flags with bit 11 set to 1, then set
     * bit 11 to 1 in BOTH attributes.  If there are 15 or fewer scans in
     * the Bit QA Flags with bit 11 set to 1, then set bit 11 in BOTH
     * attributes to zero. This is done as a stop-gap measure to prevent
     * the granule-level flag in the MCST QA databases from being in error.
     */

  day_bands_at_night = 0;
  day_bands_at_night_mask = 0x00000800;
  for (S = 0; S < L1A_Gran->num_scans - 1; S++)
  {
    if (L1B_Scan_Meta->Bit_QA_Flags[S] & day_bands_at_night_mask)
       day_bands_at_night++;
  }
  if (day_bands_at_night > 15)
  {
    QA->QA_common.bit_QA_flags_change |= day_bands_at_night_mask;
    QA->QA_common.bit_QA_flags_last_value |= day_bands_at_night_mask;
  }
  else
  {
    QA->QA_common.bit_QA_flags_change &= ~day_bands_at_night_mask;
    QA->QA_common.bit_QA_flags_last_value &= ~day_bands_at_night_mask;
  }

    /*******************************************************************/
 
  L1B_Gran_Meta->L1A_Gran_sd_id    = L1A_Gran->sd_id;
  L1B_Gran_Meta->L1B_Gran_sd_id[0] = L1B_Gran->sd_id[0];
  L1B_Gran_Meta->L1B_Gran_sd_id[1] = L1B_Gran->sd_id[1];
  L1B_Gran_Meta->L1B_Gran_sd_id[2] = L1B_Gran->sd_id[2];
  L1B_Gran_Meta->incomplete_scans  = L1A_Gran->incomplete_scans;
  L1B_Gran_Meta->max_ev_frames     = L1A_Gran->max_ev_frames;
  L1B_Gran_Meta->num_scans         = L1A_Gran->num_scans;
  L1B_Gran_Meta->num_day_scans     = L1A_Gran->num_day_scans;
  L1B_Gran_Meta->num_night_scans   = L1A_Gran->num_night_scans;
  L1B_Gran_Meta->Extract_Pixel_Offset  = L1A_Gran->Extract_Pixel_Offset;
  L1B_Gran_Meta->Extract_Pixel_Count   = L1A_Gran->Extract_Pixel_Count;
  L1B_Gran_Meta->Extract_Line_Offset   = L1A_Gran->Extract_Line_Offset;
  L1B_Gran_Meta->Extract_Line_Count    = L1A_Gran->Extract_Line_Count;
 
  for (B = 0; B < NUM_BANDS; B++)
  {   
    L1B_Gran_Meta->total_pixels[B]     = L1B_Gran->total_pixels[B];
    L1B_Gran_Meta->valid_pixels[B]     = L1B_Gran->valid_pixels[B]; 
    L1B_Gran_Meta->saturated_pixels[B] = L1B_Gran->saturated_pixels[B];
    L1B_Gran_Meta->missing_pixels[B]   = L1B_Gran->missing_pixels[B]; 
    L1B_Gran_Meta->validEVPercent[B] = 
        100 *(float32)L1B_Gran_Meta->valid_pixels[B] / 
          L1B_Gran_Meta->total_pixels[B]; 
    L1B_Gran_Meta->satEVPercent[B]   = 
        100 *(float32)L1B_Gran_Meta->saturated_pixels[B] / 
          L1B_Gran_Meta->total_pixels[B];
    L1B_Gran_Meta->missEVPercent[B]  = 
        100 *(float32)L1B_Gran_Meta->missing_pixels[B] / 
          L1B_Gran_Meta->total_pixels[B];   
 
      /* The bad data flag is not set in Reflective_Cal and Emissive_Cal for
       * dead detectors, because the EV data for that pixel may be interpolated.
       * Check the number of interpolated pixels and dead detector pixels for
       * each band. If they are not equal, set this flag.
       * should also consider dead subframe
       */
 
    if (L1B_Gran->bad_data_flag[B] == 1 ||
          L1B_Gran->interpolated_pixels[B] != 
            (L1B_Gran->dead_detector_pixels[B]+L1B_Gran->dead_subframe_pixels[B])) 
      L1B_Gran->bad_data_flag[B] = 1;  
  }
  
  missing_pixels = 0;
  outofbound_pixels = 0;
  interpolated_pixels = 0;
  total_pixels = 0;

  for (B = 0; B < NUM_250M_BANDS; B++)
  {
    missing_pixels += L1B_Gran->missing_pixels[B];
    outofbound_pixels += L1B_Gran->saturated_pixels[B] 
                         + L1B_Gran->negative_value_below_noise_pixels[B];
    interpolated_pixels += L1B_Gran->interpolated_pixels[B]; 
    total_pixels += L1B_Gran->total_pixels[B];
  }

  L1B_Gran_Meta->qapercent_missing_250m = 
           (int32)(100 * (float32)missing_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_missing_250m < 0 || 
      L1B_Gran_Meta->qapercent_missing_250m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
        "QAPERCENTMISSINGDATA for 250m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  } 
    
  L1B_Gran_Meta->qapercent_outofbound_250m = 
           (int32)(100 * (float32)outofbound_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_outofbound_250m < 0 ||
      L1B_Gran_Meta->qapercent_outofbound_250m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
        "QAPERCENTOUTOFBOUNDSDATA for 250m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->qapercent_interpolated_250m = 
           (int32)(100 * (float32) interpolated_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_interpolated_250m < 0 ||
      L1B_Gran_Meta->qapercent_interpolated_250m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
       "QAPERCENTINTERPOLATEDDATA for 250m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  missing_pixels = 0;
  outofbound_pixels = 0;
  interpolated_pixels = 0;
  total_pixels = 0;
  
  for (B = NUM_250M_BANDS; B < NUM_250M_BANDS + NUM_500M_BANDS; B++)
  {
    missing_pixels += L1B_Gran->missing_pixels[B];
    outofbound_pixels += L1B_Gran->saturated_pixels[B] 
                         + L1B_Gran->negative_value_below_noise_pixels[B]; 
    interpolated_pixels += L1B_Gran->interpolated_pixels[B];
    total_pixels += L1B_Gran->total_pixels[B];
  } 

  L1B_Gran_Meta->qapercent_missing_500m = 
           (int32)(100 * (float32)missing_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_missing_500m < 0 ||
      L1B_Gran_Meta->qapercent_missing_500m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
       "QAPERCENTMISSINGDATA for 500m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->qapercent_outofbound_500m = 
           (int32)(100 * (float32)outofbound_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_outofbound_500m < 0 ||
      L1B_Gran_Meta->qapercent_outofbound_500m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
       "QAPERCENTOUTOFBOUNDSDATA for 500m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->qapercent_interpolated_500m = 
           (int32)(100 * (float32) interpolated_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->qapercent_interpolated_500m < 0 ||
      L1B_Gran_Meta->qapercent_interpolated_500m > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
       "QAPERCENTINTERPOLATEDDATA for 500m bands is out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  missing_pixels = 0;
  outofbound_pixels = 0;
  interpolated_pixels = 0;
  total_pixels = 0;   
  for(; B < NUM_REFLECTIVE_BANDS - 1; B++) /* band 26 */
  {
    missing_pixels += L1B_Gran->missing_pixels[B];
    outofbound_pixels += L1B_Gran->saturated_pixels[B] 
                         + L1B_Gran->negative_value_below_noise_pixels[B]; 
    interpolated_pixels += L1B_Gran->interpolated_pixels[B];
    total_pixels += L1B_Gran->total_pixels[B];
  } 

  missing_pixels += L1B_Gran->missing_pixels[MODIS_BAND26_INDEX];
  outofbound_pixels += L1B_Gran->saturated_pixels[MODIS_BAND26_INDEX] 
     + L1B_Gran->negative_value_below_noise_pixels[MODIS_BAND26_INDEX]; 
  interpolated_pixels += L1B_Gran->interpolated_pixels[MODIS_BAND26_INDEX];
  total_pixels += L1B_Gran->total_pixels[MODIS_BAND26_INDEX];
  
  L1B_Gran_Meta->refl_1km_qapercent_missing = 
           (int32)(100 * (float32)missing_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->refl_1km_qapercent_missing < 0 ||
      L1B_Gran_Meta->refl_1km_qapercent_missing > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus, 
       "QAPERCENTMISSINGDATA for 1km reflective bands is "
        "out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->refl_1km_qapercent_outofbound = 
           (int32)(100 * (float32)outofbound_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->refl_1km_qapercent_outofbound < 0 ||
      L1B_Gran_Meta->refl_1km_qapercent_outofbound > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
       "QAPERCENTOUTOFBOUNDSDATA for 1km reflective bands is "
        "out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->qapercent_interpolated_refl_1km = 
           (int32)(100 * (float32) interpolated_pixels/total_pixels + 0.5); 
  if (L1B_Gran_Meta->qapercent_interpolated_refl_1km < 0 ||
      L1B_Gran_Meta->qapercent_interpolated_refl_1km > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "QAPERCENTINTERPOLATEDDATA for 1km reflective bands is"
                " out of range [0, 100]",
                NULL, 0, NULL, True);
  }
  

  missing_pixels = 0;
  outofbound_pixels = 0;
  interpolated_pixels = 0;
  total_pixels = 0;  
 
  for(; B < NUM_BANDS; B++)
  {
    if(B == MODIS_BAND26_INDEX)
      continue;
   
    missing_pixels += L1B_Gran->missing_pixels[B];
    outofbound_pixels += L1B_Gran->saturated_pixels[B] 
                         + L1B_Gran->negative_value_below_noise_pixels[B]; 
    interpolated_pixels += L1B_Gran->interpolated_pixels[B];
    total_pixels += L1B_Gran->total_pixels[B];
  } 

  L1B_Gran_Meta->emiss_qapercent_missing = 
           (int32)(100 * (float32)missing_pixels/total_pixels + 0.5);
  if (L1B_Gran_Meta->emiss_qapercent_missing < 0 ||
      L1B_Gran_Meta->emiss_qapercent_missing > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "QAPERCENTMISSINGDATA for emissive bands is "
                "out of range [0, 100]",
                NULL, 0, NULL, True);
  }

  L1B_Gran_Meta->emiss_qapercent_outofbound = 
           (int32)(100 * (float32)outofbound_pixels/total_pixels + 0.5);  
  if (L1B_Gran_Meta->emiss_qapercent_outofbound < 0 ||
      L1B_Gran_Meta->emiss_qapercent_outofbound > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "QAPERCENTOUTOFBOUNDSDATA for emissive bands is "
                "out of range [0, 100]",
                NULL, 0, NULL, True);
  }
 
  L1B_Gran_Meta->qapercent_interpolated_emiss = 
           (int32)(100 * (float32)interpolated_pixels/total_pixels + 0.5); 
  if (L1B_Gran_Meta->qapercent_interpolated_emiss < 0 ||
      L1B_Gran_Meta->qapercent_interpolated_emiss > 100) {
    returnStatus = MODIS_F_OUT_OF_RANGE;
    L1BErrorMsg(location, returnStatus,
                "QAPERCENTINTERPOLATEDDATA for emissive bands is"
                " out of range [0, 100]", 
                NULL, 0, NULL, True);
  }
 
  L1B_Gran_Meta->FPSetPointState = PP->PP_Emiss.fp_set_point_state[0];
  for (S = 1; S < L1A_Gran->num_scans; S++) 
    if(L1B_Gran_Meta->FPSetPointState != PP->PP_Emiss.fp_set_point_state[S])
      L1B_Gran_Meta->FPSetPointState = 0;
    
  returnStatus = MODIS_S_OK;
  
  /* New QA */
 
  L1B_Gran_Meta->Door_Screen_Configuration = 0; 
  for (Door_Screen_Config = NAD_Door; 
      Door_Screen_Config < NUM_Config; Door_Screen_Config++)
  {
    returnStatus = read_vdata (L1A_Gran->v_id, 
                               0, 
                               1, 
                               temp[Door_Screen_Config].vname, 
                               temp[Door_Screen_Config].fname, 
                               (VOIDP)&temp[Door_Screen_Config].buffer);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "read_vdata", FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    if (temp[Door_Screen_Config].buffer > 1)
    {
      char msgbuf[256];
      sprintf(msgbuf, "Value[0] of \"%s\" is out of the range [0,1].",
              temp[Door_Screen_Config].fname);
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, msgbuf, NULL, FIRST_L1A_GRANULE,
                  Invalid_MOD01_Msg, True);
      return returnStatus;
    }
  } 
 
  if (temp[NAD_Door].buffer == 0)
    L1B_Gran_Meta->Door_Screen_Configuration |= 0x80;
  if (temp[SVD_Door].buffer == 1)
    L1B_Gran_Meta->Door_Screen_Configuration |= 0x40;    
  if (temp[SDD_Door].buffer == 1)
    L1B_Gran_Meta->Door_Screen_Configuration |= 0x20;
  if (temp[SDD_Screen].buffer == 1)
    L1B_Gran_Meta->Door_Screen_Configuration |= 0x10;
   
  for (B = 0; B < NUM_REFLECTIVE_BANDS - 1; B++)
  {
    if(L1B_Gran->bad_data_flag[B] == 1)
      L1B_Gran_Meta->Reflective_Band_Identification[B] = 1;
    else
      L1B_Gran_Meta->Reflective_Band_Identification[B] = 0;
  }
  /* for band 26 */
  if(L1B_Gran->bad_data_flag[MODIS_BAND26_INDEX] == 1)
    L1B_Gran_Meta->Reflective_Band_Identification[B] = 1;
  else
    L1B_Gran_Meta->Reflective_Band_Identification[B] = 0;  

  B_38 = NUM_REFLECTIVE_BANDS - 1;

  for (B = 0; B < NUM_EMISSIVE_BANDS; B++, B_38++)
  {
    if(B_38 == MODIS_BAND26_INDEX)
      B_38++;   
    if(L1B_Gran->bad_data_flag[B_38] == 1)
      L1B_Gran_Meta->Emissive_Band_Identification[B] = 1;
    else
      L1B_Gran_Meta->Emissive_Band_Identification[B] = 0;
  }

  L1B_Gran_Meta->All_L1B_Error_Flag_Off = 0;

  for (B = 0; B < NUM_EMISSIVE_BANDS; B++)
    for(D = 0; D < DETECTORS_PER_1KM_BAND; D++)
    {
      L1B_Gran_Meta->Thermal_Detector_Noise[B][D] = 
          QA->QA_emiss.NEdL[B * 10 + D] ;
      L1B_Gran_Meta->Thermal_Detector_Relative_Response_Change[B][D] =
          QA->QA_emiss.change_b1[B * 10 + D] ;
    }
 
  /* 
   * Determine the Electronics Configuration Status and the Electronics 
   * Configuration Change. 
   */ 

  returnStatus = Get_Elec_Config_Status(&QA->QA_common, 
                                        L1A_Gran->v_id, 
                                        L1A_Gran->num_scans,
                                        L1B_Gran_Meta->Elec_config_status, 
                                        L1B_Gran_Meta->Elec_config_change);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, 
        "Get_Elec_Config_Status", 0, NULL, True);
    return returnStatus;
  }

  /* Store the Earth-Sun Distance into L1B_Gran_Meta */

  L1B_Gran_Meta->Earth_Sun_Dist = L1B_Gran->Earth_Sun_Dist;

  return(returnStatus);
}   

  
PGSt_SMF_status Write_Global_Metadata
                       (L1B_Gran_Metadata_t *L1B_Gran_Meta,
                        QA_Data_t           *QA,
                        lookup_tables_t     *tables,
                        int32               OBC_sd_id,
                        boolean             skip_night_hi_res)
/*
!C****************************************************************
!Description:  This function writes the MODIS Level 1B Product Granule
               Metadata stored as global attributes, which are described
               in 1.3) of the EV file specifications, and MODIS Level 1B
               QA Granule Metadata, which are described in 1.4) of
               of the EV file specifications into EV files and OBC file. 

!Input Parameters: 
     L1B_Gran_Metadata_t  *L1B_Gran_Meta   contains granule metadata 
     QA_Data_t            *QA              contains QA granule metadata
     lookup_tables_t      *tables          contains granule metadata 
                                           implemented directly from LUTs
     int32                OBC_sd_id        OBC file id for sds interface
     boolean skip_night_hi_res       True if and only if all scans are
                                       NIGHT mode and runtime parameter
                                       Write_Night_Mode_HiRes_Data is False.   

!Output Parameters:
     none

!Revision History:
 (continue at top of the file)

  Revision 01.03  January 17, 2002   Razor Issue #172
  Improve portability of code to 64-bit mode.
  Cast strlen returns to int32 in calls to SDSetattr.
  Alice Isaacman, SAIC GSO  (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.02 November 19, 2000 (Razor Issue #169)
  Added input parameter skip_night_hi_res and changed logic so that
  250m and 500m data are not written when granule has no day mode
  scans and runtime parameter Write_Night_Mode_HiRes_Data is False.
 Alice Isaacman   SAIC GSO   (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.01 1999/8/20
 Moved Manager's quality indices to common QA
 Zhenying Gu (zgu@gscmail.gsfc.nasa.gov)

 Revision 01.00 1998/10/30 
 Original Development(This function is part of the original function 
 Write_Gran_Metadata())
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
  PGSt_SMF_status returnStatus;
  char            *location = "Write_Global_Metadata";
  int32           out_sd_id[4] = {0, 0, 0, 0};
  int32           R;
  int32           hdf_return = FAIL;
  int32           lun;
  /*
   * Starting output file index; changes to INDEX_L1B_EV_1000M_FILE
   * when HiRes night mode is not written.
   */
  int16 start_output_file_index = INDEX_L1B_EV_250M_FILE;    
  
  out_sd_id[0] = L1B_Gran_Meta->L1B_Gran_sd_id[0];
  out_sd_id[1] = L1B_Gran_Meta->L1B_Gran_sd_id[1];
  out_sd_id[2] = L1B_Gran_Meta->L1B_Gran_sd_id[2];
  out_sd_id[3] = OBC_sd_id;

  /*---------------------------------------
    Write Native HDF metadata 
  ---------------------------------------*/
  /*
   * If this is a NIGHT mode granule and the runtime parameter
   * Write_Night_Mode_HiRes_Data is False, do not write the 250m and 500m
   * resolution files, and so start with output file number 2, not 0.
   */

  if (skip_night_hi_res == True)
    start_output_file_index = INDEX_L1B_EV_1000M_FILE;
  else
    start_output_file_index = INDEX_L1B_EV_250M_FILE;

  for ( R = start_output_file_index ; R < NUM_OUTPUT_FILES ; R++)
  { 
    /* Get logical number for the file for L1BErrorMsg */
    
    switch (R)
    {
      case INDEX_L1B_EV_250M_FILE:  lun = L1B_EV_250M_FILE;
              break;
      case INDEX_L1B_EV_500M_FILE:  lun = L1B_EV_500M_FILE;
              break;
      case INDEX_L1B_EV_1000M_FILE: lun = L1B_EV_1000M_FILE;
              break;
      case INDEX_L1B_OBC_FILE:      lun = L1B_OBC_FILE;
              break;
      default: lun = 0;
    }

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Number of Scans",
                            DFNT_INT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->num_scans);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Number of Scans.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Number of Day mode scans",
                            DFNT_INT32, 1, 
                           (VOIDP)&L1B_Gran_Meta->num_day_scans);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Number of Day mode scans.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Number of Night mode scans",
                            DFNT_INT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->num_night_scans);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Number of Night mode scans.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Incomplete Scans",
                            DFNT_INT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->incomplete_scans);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Incomplete Scans.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Max Earth View Frames",
                            DFNT_INT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->max_ev_frames);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Max Earth View Frames.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "%Valid EV Observations", 
                            DFNT_FLOAT32, NUM_BANDS, 
                            (VOIDP)L1B_Gran_Meta->validEVPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write %Valid EV Observations.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                  "%Saturated EV Observations", 
                            DFNT_FLOAT32, NUM_BANDS, 
                            (VOIDP)L1B_Gran_Meta->satEVPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write %Saturated EV Observations.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% L1A EV All Scan Data are Missing",
                            DFNT_FLOAT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->missAllScanDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % L1A EV All Scan Data are Missing.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                  "% L1A EV RSB DN Not in Day Mode",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->nightRSBPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % L1A EV RSB DN Not in Day Mode.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% L1A EV DN Missing Within Scan",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->missInScanDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % L1A EV DN Missing Within Scan.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% Dead Detector EV Data",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->deadDetectorDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % Dead Detector EV Data.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R],
                            "% Dead Subframe EV Data",
                            DFNT_FLOAT32, NUM_HIGH_RESOLUTION_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->deadSubframeDataPercent);

    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % Dead Subframe EV Data.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% Sector Rotation EV Data",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->sectorRotateDataPercent);

    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % Sector Rotation EV Data.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% Saturated EV Data",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->saturatedDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % Saturated EV Data.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% TEB EV Data With Moon in SVP",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->moonInSVPTEBDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % TEB EV Data With Moon in SVP.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% EV Data Where Cannot Compute BG DN",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->noBGDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % EV Data Where Cannot Compute BG DN.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                           "% RSB EV Data With dn** Below Scale",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->badDNStarStarRSBDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % RSB EV Data With dn** Below Scale.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% EV Data Where Nadir Door Closed",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->NADClosedDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % EV Data Where Nadir Door Closed.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "% EV Data Not Calibrated",
                            DFNT_FLOAT32, NUM_DETECTORS,
                            (VOIDP)L1B_Gran_Meta->uncalibratedDataPercent);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write % EV Data Not Calibrated.",
                  "SDsetattr", lun, NULL, True);

    write_extract_metadata(out_sd_id[R],
                           L1B_Gran_Meta->Extract_Pixel_Offset,
                           L1B_Gran_Meta->Extract_Pixel_Count,
                           L1B_Gran_Meta->Extract_Line_Offset,
                           L1B_Gran_Meta->Extract_Line_Count);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Bit QA Flags Last Value",
                            DFNT_UINT32, 1,
                            (VOIDP)&QA->QA_common.bit_QA_flags_last_value);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Bit QA Flags Last Value.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Bit QA Flags Change",
                            DFNT_UINT32, 1,
                            (VOIDP)&QA->QA_common.bit_QA_flags_change);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Bit QA Flags Change.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Granule Average QA Values",
                            DFNT_FLOAT32, MAX_NUM_GRAN_AVERAGES,
                            (VOIDP) QA->QA_common.granule_averages);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Granule Average QA Values.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Electronics Redundancy Vector",
                            DFNT_UINT32, 2, 
                            (VOIDP)L1B_Gran_Meta->Elec_config_status);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Electronics Redundancy Vector.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Electronics Configuration Change",
                            DFNT_UINT32, 2, 
                            (VOIDP)L1B_Gran_Meta->Elec_config_change);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Electronics Configuration Change.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                      "Reflective LUT Serial Number and Date of Last Change",
                      DFNT_CHAR8, (int32)strlen(tables->refl.Serial_Number),
                      (VOIDP)tables->refl.Serial_Number);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Reflective LUT Serial Number.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                        "Emissive LUT Serial Number and Date of Last Change",
                        DFNT_CHAR8, (int32)strlen(tables->emiss.Serial_Number),
                        (VOIDP)tables->emiss.Serial_Number);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Emissive LUT Serial Number.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "QA LUT Serial Number and Date of Last Change",
                            DFNT_CHAR8, 
                            (int32)strlen(tables->QA.common_QA_tables.Serial_Number),
                            (VOIDP)tables->QA.common_QA_tables.Serial_Number);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write QA LUT Serial Number.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Focal Plane Set Point State", 
                            DFNT_INT8, 1, 
                            (VOIDP)&L1B_Gran_Meta->FPSetPointState);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Focal Plane Set Point State.",
                  "SDsetattr", lun, NULL, True);
     
    hdf_return = SDsetattr (out_sd_id[R], 
                            "Doors and Screens Configuration", DFNT_INT8, 1, 
                            (VOIDP)&L1B_Gran_Meta->Door_Screen_Configuration);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Doors and Screens Config.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Reflective Bands With Bad Data", 
                            DFNT_INT8,
                            NUM_REFLECTIVE_BANDS, 
                      (VOIDP)L1B_Gran_Meta->Reflective_Band_Identification);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Reflective Bands With Bad Data.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Emissive Bands With Bad Data", 
                            DFNT_INT8, 
                            NUM_EMISSIVE_BANDS, 
                      (VOIDP)L1B_Gran_Meta->Emissive_Band_Identification);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Emissive Bands With Bad Data.",
                  "SDsetattr", lun, NULL, True);
 
    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Black Body Thermistors", 
                            DFNT_UINT8, NUM_THERMISTORS, 
                            (VOIDP)QA->QA_emiss.noise_T_bb);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Black Body Thermistors.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Average BB Temperature", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_bb_avg);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Average BB Temperature.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in LWIR FPA Temperature", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_lwir);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in LWIR FPA Temperature.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in MWIR FPA Temperature", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_mwir);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in MWIR FPA Temperature.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Scan Mirror Thermistor #1", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_mir1);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Scan Mirror Thermistor #1.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Scan Mirror Thermistor #2", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_mir2);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Scan Mirror Thermistor #2.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Scan Mirror Thermistor Average", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_mir_avg);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Scan Mirror Thermistor Avg.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Instrument Temperature", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_ins);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Instrument Temperature.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Cavity Temperature", DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_emiss.noise_T_cav);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Cavity Temperature.",
                  "SDsetattr", lun, NULL, True);

    returnStatus = write_sds_rank2 (out_sd_id[R], 
                                    "Noise in Thermal Detectors",
                                    "number of emissive bands", 
                                    "detectors per 1km band", 
                                    NUM_EMISSIVE_BANDS, 
                                    DETECTORS_PER_1KM_BAND,   
                                    "uint8", 
                                    L1B_Gran_Meta->Thermal_Detector_Noise);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Noise in Thermal Detectors.",
                  "write_sds_rank2", lun, NULL, True);

    returnStatus = write_sds_rank2 
                    (out_sd_id[R], 
                     "Change in relative responses of thermal detectors", 
                     "number of emissive bands", 
                     "detectors per 1km band", 
                     NUM_EMISSIVE_BANDS, 
                     DETECTORS_PER_1KM_BAND, 
                     "uint8", 
           (VOIDP)L1B_Gran_Meta->Thermal_Detector_Relative_Response_Change);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Change in relative responses of "
                  "thermal detectors.",
                  "write_sds_rank2", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Temperature of NIR FPA", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_refl.noise_T_nir);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Noise in Temperature of NIR FPA.",
                  "SDsetattr", lun, NULL, True);
   
    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noise in Temperature of Vis FPA", 
                            DFNT_UINT8, 1, 
                            (VOIDP)&QA->QA_refl.noise_T_vis);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Noise in Temperature of Vis FPA.",
                  "SDsetattr", lun, NULL, True);

    returnStatus = write_sds_rank3 (out_sd_id[R],
                                    "DC Restore Change for Thermal Bands", 
                                    "number of scans",
                                    "number of emissive bands", 
                                    "detectors per 1km band", 
                                    L1B_Gran_Meta->num_scans, 
                                    NUM_EMISSIVE_BANDS,
                                    DETECTORS_PER_1KM_BAND, 
                                    "int8", 
                                    QA->QA_emiss.change_dc_restore);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write DC Restore Change for Thermal Bands.",
                  "write_sds_rank3", lun, NULL, True);

    returnStatus = write_sds_rank3 
                           (out_sd_id[R], 
                            "DC Restore Change for Reflective 250m Bands",
                            "number of scans", 
                            "number of 250m bands",
                            "detectors per 250m band", 
                            L1B_Gran_Meta->num_scans,
                            NUM_250M_BANDS, 
                            DETECTORS_PER_250M_BAND, 
                            "int8",
                            QA->QA_refl.change_dc_restore_250m);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
          "Could not write DC Restore Change for Reflective 250m Bands.",
                "write_sds_rank3", lun, NULL, True);
    
    returnStatus = write_sds_rank3                                          
                            (out_sd_id[R], 
                            "DC Restore Change for Reflective 500m Bands",
                            "number of scans", 
                            "number of 500m bands",
                            "detectors per 500m band", 
                            L1B_Gran_Meta->num_scans,
                            NUM_500M_BANDS, 
                            DETECTORS_PER_500M_BAND, 
                            "int8",
                            QA->QA_refl.change_dc_restore_500m);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
          "Could not write DC Restore Change for Reflective 500m Bands.",
                  "write_sds_rank3", lun, NULL, True);
   
    returnStatus = write_sds_rank3
                            (out_sd_id[R], 
                            "DC Restore Change for Reflective 1km Bands",
                            "number of scans", 
                            "number of 1km reflective bands",
                            "detectors per 1km band", 
                            L1B_Gran_Meta->num_scans,
                            NUM_1000M_REFL_BANDS, 
                            DETECTORS_PER_1KM_BAND, 
                            "int8",
                            QA->QA_refl.change_dc_restore_1km);
    if (returnStatus != MODIS_S_OK)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
            "Could not write DC Restore Change for Reflective 1km Bands.",
                  "write_sds_rank3", lun, NULL, True);
    
    hdf_return = SDsetattr (out_sd_id[R], 
                            "Dead Detector List", 
                            DFNT_INT8, 
                            NUM_DETECTORS,
                            (VOIDP)tables->QA.common_QA_tables.dead_detector);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Dead Detector List.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Noisy Detector List", 
                            DFNT_INT8, 
                            NUM_DETECTORS,
                    (VOIDP)tables->QA.common_QA_tables.noisy_detector);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Noisy Detector List.",
                  "SDsetattr", lun, NULL, True);

    /*output noisy and dead subframe list */
    hdf_return = SDsetattr (out_sd_id[R],
                            "Dead Subframe List",
                            DFNT_INT8,
                            NUM_HIGH_RESOLUTION_SUBFRAMES,
                            (VOIDP)tables->QA.common_QA_tables.dead_subframe);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Dead Subframe List.",
                  "SDsetattr", lun, NULL, True);

    hdf_return = SDsetattr (out_sd_id[R],
                            "Noisy Subframe List",
                            DFNT_INT8,
                            NUM_HIGH_RESOLUTION_SUBFRAMES,
                    (VOIDP)tables->QA.common_QA_tables.noisy_subframe);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Noisy Subframe List.",
                  "SDsetattr", lun, NULL, True);

    /* Detector Quality Flag */

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Detector Quality Flag",
                            DFNT_UINT8, 
                            NUM_DETECTORS, 
                   (VOIDP)tables->QA.common_QA_tables.Detector_Quality_Flag);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Detector Quality Flag.",
                  "SDsetattr", lun, NULL, True);

    /* Detector Quality Flag2 */

    hdf_return = SDsetattr (out_sd_id[R],
                            "Detector Quality Flag2",
                            DFNT_UINT8,
                            NUM_HIGH_RESOLUTION_DETECTORS,
                   (VOIDP)tables->QA.common_QA_tables.Detector_Quality_Flag2);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write Detector Quality Flag2.",
                  "SDsetattr", lun, NULL, True);
     
    /* Earth-Sun distance */
    
    hdf_return = SDsetattr (out_sd_id[R], 
                            "Earth-Sun Distance",
                            DFNT_FLOAT32, 1, 
                            (VOIDP)&L1B_Gran_Meta->Earth_Sun_Dist);
    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write Earth-Sun distance.",
                  "SDsetattr", lun, NULL, True);

    /* Solar irradiance */

    hdf_return = SDsetattr (out_sd_id[R], 
                            "Solar Irradiance on RSB Detectors over pi",
                            DFNT_FLOAT32, 
                            NUM_REFLECTIVE_DETECTORS,
                            (VOIDP) tables->refl.E_sun_over_pi);

    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR, 
                  "Could not write E_sun_over_pi.",
                  "SDsetattr", lun, NULL, True);

    /* set doi global attribute */
    hdf_return = SDsetattr (out_sd_id[R],
                            L1B_Gran_Meta->doi_attr_name,
                            DFNT_CHAR8, (int32)strlen(L1B_Gran_Meta->doi_para_value[R]),
                            &L1B_Gran_Meta->doi_para_value[R][0]);

    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write identifier_product_doi.",
                  "SDsetattr", lun, NULL, True);

    /* set doi authority global attribute */
    hdf_return = SDsetattr (out_sd_id[R],
                            L1B_Gran_Meta->doi_authority_attr_name,
                            DFNT_CHAR8, (int32)strlen(L1B_Gran_Meta->doi_authority_para_value),
                            L1B_Gran_Meta->doi_authority_para_value);

    if (hdf_return == FAIL)
      L1BErrorMsg(location, MODIS_F_WRITE_ERROR,
                  "Could not write identifier_product_doi_authority.",
                  "SDsetattr", lun, NULL, True);


  }

  return(MODIS_S_OK);
}

PGSt_SMF_status Get_Electronics_Status(int32   v_id, 
                                       int32   num_scans,
                                       char    *vname, 
                                       char    *fname,
                                       int16   *final_value, 
                                       int16   *is_changed,
                                       boolean *no_valid_value)
/*
!C**********************************************************************
!Description: This function computes the final status of the electronics
              and determines if any change occurred within the granule
              for the electronics.

!Input Parameters:
   int32   v_id       file id used for vdata interface
   int32   num_scans  number of scans in the granule represented by v_id
   char *  vname      vdata name of the telemetry point for the interested
                      electronics
   char *  fname      telemetry field name for the electronics

!Output Parameters:
   int16 * final_value  the address of the final value of the telemetry
                        point corresponding to the interested electronics.
                        The final value is in set [0, 1].
   int16 * is_changed   the address of the change status of the telemetry
                        point. 0 = no change, 1 = change occured within
                        granule.
   boolean *no_valid_value  True = no valid value for the telemetry field
                            within the granule, False = there is at least
                            one valid value for the telemetry field within
                            the granule.

!Revision History:
 (continue at top of the file)

 Revision 01.00 Sept. 7, 1999
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   If there is no valid scan for the telemetry field corresponding
 to the electronics, the final value of the telemetry field is set
 to 0 and the change status is set to 0.

!END********************************************************************
*/

{
  PGSt_SMF_status  returnStatus;
  char *location = "Get_Electronics_Status";
  uint16  values[MAX_NUM_SCANS];  /* values of the telemetry point for
                                     all scans in the granule */
  uint16  temp_final_value;       /* temporary variable to store the
                                     last valid value */
  uint16  last_valid_scans[MAX_NUM_SCANS]; /* the scan from which the
                               information is taken for the current scan */
  boolean first_valid_value_found; /* flag to indicate if there has been a
                               valid value found while looping through scans */

  int16   S;  /* scan index */

  /* Check the input parameters */

  if (vname == NULL || fname == NULL || final_value == NULL || 
      is_changed == NULL)
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus,
                "One of vname, fname, final_value or is_changed is invalid.",
                NULL, 0, NULL, True);
    return returnStatus;
  }

  if (num_scans <= 0 || num_scans > MAX_NUM_SCANS)
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, "num_scans is invalid.",
                NULL, 0, NULL, True);
    return returnStatus;
  }

  /* Initialization */

  *is_changed = 0;
  first_valid_value_found = False;
  temp_final_value = 0;

  /* Read the vdata field and the "LAST_VALID_SCAN" field */

  returnStatus = read_vdata(v_id, 0, num_scans, vname, fname, values);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_vdata", 0, NULL, False);
    return returnStatus;
  }

  returnStatus = read_vdata(v_id, 0, num_scans, vname, 
      "LAST_VALID_SCAN", last_valid_scans);
  if (returnStatus != MODIS_S_OK)
  {
    L1BErrorMsg(location, returnStatus, NULL, "read_vdata", 0, NULL, False);
    return returnStatus;
  }

  for (S = 0; S < num_scans; S++)
  {
    if (values[S] != 0 && values[S] != 1)
    {
      char msgbuf[256];
      sprintf(msgbuf, 
              "Value[%d] of telemetry \"%s\" is out of the range [0,1].",
              S, fname);
      returnStatus = MODIS_F_OUT_OF_RANGE;
      L1BErrorMsg(location, returnStatus, msgbuf, NULL, 0, 
                  Invalid_MOD01_Msg, False);
      return returnStatus;
    }


    /*
     * If current value is valid and different from the temporary final
     * value, it means there is a change. Update the final value.
     */

 
    if (last_valid_scans[S] != L1A_MISSING_ENG_PACKET)
    {
      if (first_valid_value_found == True && values[S] != temp_final_value)
        *is_changed = 1;

      temp_final_value = values[S];
      first_valid_value_found = True;
    }
  }

  *final_value = temp_final_value;
  *no_valid_value = !first_valid_value_found;

  return returnStatus;
}

PGSt_SMF_status Get_Elec_Config_Status_Per_Gran
                                 (int32 v_id,
                                  int32 num_scans,
                                  uint32 *Elec_config_status,
                                  uint32 *Elec_config_change,
                                  uint32 *Elec_config_invalid_flag)
/*
!C**********************************************************************
!Description: This function determines the final status of a set of 
              telemetry fields electronics and determines if any change 
              occurred within the granule for each of the fields. It also 
              indicates if there is valid data for the telemetry fields 
              within the granule.

!Input Parameters:
   int32    v_id           vdata interface for L1A file
   int32    num_scans      number of scans in the L1A file

!Output Parameters:

   uint32   Elec_config_status[E_VECTOR_SIZE]        
                 the final valid value for the corresponding
                 telemetry point within the granule with one
                 bit representing one field.
                 Values are 0=OFF and 1=ON
   uint32   Elec_config_change[E_VECTOR_SIZE]        
                 identifies if any change occurred within
                 the granule for the corresponding telemetry
                 point with one bit representing one field.
   uint32   Elec_config_invalid_flag[E_VECTOR_SIZE]  
                 identifies if there is no valid value for
                 the telemetry fields within the whole granule

!Revision History:
 (continue at top of the file)

 Revision 01.01  March 2003, Razor Issue #173
 Enclosed the rows in the initializer of array-of-struct "Elec_config_fields" 
 with braces for ANSI-C compliance. 
 Liqin Tan, SAIC GSO (ltan@saicmodis.com)

 Revision 01.00 Oct. 29, 1999
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
  If the final valid value for a telemetry field is 0, set the
  corresponding bit (see file specification for which bit represents
  which telemetry field.) in Elec_config_status to be 0. If the final
  valid value for a telemetry field  is 1, set the bit to 1. If there
  is no valid value for a telemetry within the  granule, set the
  corresponding bit in Elec_config_status to 0 and set the
  corresponding bit in Elec_config_invalid_flag to 1. If there is valid
  value change for a telemetry field within the granule, set the
  corresponding bit in Elec_config_change to 1. Otherwise set it to
  0.  
!END********************************************************************
*/
{
  PGSt_SMF_status returnStatus = MODIS_S_OK;
  char *location = "Get_Elec_Config_Status_Per_Gran";
  int16 i;
  int16 final_value;  /* the value of electronics field in the last scan */
  int16 is_changed;   /* the change status of the electronics field in */
                      /*     the granule */
  boolean no_valid_value;  /* flag indicating if there is valid value */
                      /*     for a telemetry field within the granule */
  int16 ivar;         /* index for number of word used for electronics */
                      /*     configuration status (or change) */
  int32 mask;         /* bit mask used to set the every bit of the */
                      /*     attributes */
  int32 total_fields; /* index to keep track of all the fields have */
                      /*     been looped through */
  int32 num_elec_config_elements = 0;  /* number of the telemetry fields */
                      /*     in the Electronics configuration */

  struct {
   char *vname;
   char *fname;
  }Elec_config_fields[] = {
   {"Telemetry Major Cycle 0 of 7",  "CR_BB_A_PWR_ON"},
   {"Telemetry Major Cycle 0 of 7",  "CR_BB_B_PWR_ON"},
   {"Telemetry Major Cycle 1 of 7",  "CR_CE_A_ON"},
   {"Telemetry Major Cycle 1 of 7",  "CR_CE_B_ON"},
   {"Telemetry Major Cycle 1 of 7",  "CR_CP_A_ON_M"},
   {"Telemetry Major Cycle 1 of 7",  "CR_CP_B_ON_M"},
   {"Telemetry Major Cycle 3B of 7", "CR_FI_A_ON"},
   {"Telemetry Major Cycle 3B of 7", "CR_FI_B_ON"},
   {"Telemetry Major Cycle 3B of 7", "CR_FO_BLK1_ON"},
   {"Telemetry Major Cycle 3B of 7", "CR_FO_BLK2_ON"},
   {"Telemetry Major Cycle 3B of 7", "CR_FO_BLK3_ON"},
   {"Telemetry Major Cycle 3B of 7", "CR_FO_BLK4_ON"},
   {"Telemetry Major Cycle 3C of 7", "CR_FR_A_ON"},
   {"Telemetry Major Cycle 3C of 7", "CR_FR_B_ON"},
   {"Telemetry Major Cycle 4B of 7", "CS_FR_PC_DCR_ON"},
   {"Telemetry Major Cycle 4B of 7", "CS_FR_PV_DCR_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PCLW_A_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PCLW_B_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PVLW_A_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PVLW_B_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVSM_A_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVSM_B_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVNIR_A_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVNIR_B_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVVIS_A_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVVIS_B_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PVLWA_ECAL_ON"},
   {"Telemetry Major Cycle 4A of 7", "CR_PVLWB_ECAL_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVNIRA_ECALON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVNIRB_ECALON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVSMA_ECAL_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVSMB_ECAL_ON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVVISA_ECALON"},
   {"Telemetry Major Cycle 4B of 7", "CR_PVVISB_ECALON"},
   {"Telemetry Major Cycle 5A of 7", "CR_RC_SMHTR_ON"},
   {"Telemetry Major Cycle 5A of 7", "CR_RC_LWHTR_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SA_A_SCAN_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SA_B_SCAN_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SM_SDSM_A_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SM_SDSM_B_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SR_A_ON"},
   {"Telemetry Major Cycle 5B of 7", "CR_SR_B_ON"},
   {"Telemetry Major Cycle 6 of 7",  "CR_TG_A_ON"},
   {"Telemetry Major Cycle 6 of 7",  "CR_TG_B_ON"},
   {NULL,                            NULL}
  };
  
  /* check input parameters */

  if (num_scans <= 0 || num_scans > MAX_NUM_SCANS)
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, "number of scans is out of range",
                NULL, 0, NULL, True);
    return returnStatus;
  }

  /* initialization */

  total_fields = 0;

  /* count up number of elements in Elec_config_fields */

  while (Elec_config_fields[num_elec_config_elements].vname!= NULL)
    num_elec_config_elements++;

  if (num_elec_config_elements > E_VECTOR_SIZE * 32)
  {
    returnStatus = MODIS_F_NOK;
    L1BErrorMsg(location, returnStatus, "macro E_VECTOR_SIZE is too small",
                NULL, 0, NULL, True);
    return returnStatus;
  }
 
  for (ivar = 0; ivar < E_VECTOR_SIZE; ivar++)
  {
    /* initialization */

    Elec_config_status[ivar] = 0;
    Elec_config_change[ivar] = 0;
    Elec_config_invalid_flag[ivar] = 0;
   
    /* mask starts at the least significant bit */
 
    mask = 1;

    for (i = 0; i < 32 && total_fields < num_elec_config_elements; 
                i++, total_fields++)
    {
      /* 
       * Get the final value, change status and valid value information 
       * for the field.
       */

      returnStatus = Get_Electronics_Status
                               (v_id, 
                                num_scans, 
                                Elec_config_fields[total_fields].vname,
                                Elec_config_fields[total_fields].fname, 
                                &final_value, 
                                &is_changed, 
                                &no_valid_value);
      if (returnStatus != MODIS_S_OK)
      {
        L1BErrorMsg(location, returnStatus, NULL,
                    "Get_Electronics_Status", 0, NULL, False);
        return returnStatus;
      }

      /* 
       * Pack final value, change status and valid value information
       * into three output variables, respectively. Each bit in each
       * variable represents a specific value for a telemetry field.
       */
 
      if (final_value == 1)
        Elec_config_status[ivar] |= mask;

      if (is_changed == 1)
        Elec_config_change[ivar] |= mask;

      if (no_valid_value == True)
        Elec_config_invalid_flag[ivar] |= mask;

      /* shift mask to next more significant bit */

      mask *= 2;
    }
  }

  return returnStatus;
}  

       

PGSt_SMF_status Get_Elec_Config_Status(QA_Common_t *QA_common,
                                       int32       v_id,
                                       int32       num_scans,
                                       uint32      *Elec_config_status,
                                       uint32      *Elec_config_change)
/*
!C**********************************************************************
!Description: This function computes the final status of the electronics
              and determines if any change occurred within the granule
              or between the last scan of the leading granule and the first
              scan of the middle granule for the electronics. It also 
              indicates if there is valid data for the telemetry fields 
              within the granule.

!Input Parameters:
   QA_Common_t *QA_common  containing the information if the leading
                           granule is missing
   int32    v_id           vdata interface for middle L1A granule
   int32    num_scans      number of scans in the middle L1A granule

!Output Parameters:
   uint32   Elec_config_status[E_VECTOR_SIZE]   
                        the final valid value for the corresponding
                        telemetry point within the granule with one
                        bit represents one field.
                        Values are 0=OFF and 1=ON
                                                
   uint32   Elec_config_change[E_VECTOR_SIZE]   
                        identifies if any change occurred within
                        the granule or between the last valid value of
                        the leading granule and the last valid value of
                        the middle granule for the corresponding telemetry
                        point with one bit represent one field.

!Revision History:
 (continue at top of the file)

 Revision 01.00 Nov. 1, 1999
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!Team-unique Header:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:
   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
     Check the valid values of the electronics telemetry fields for every
  scan in the leading L1A granule and middle L1A granule. See function
  Get_Elec_Config_Status_Per_Gran  for the listing and order of fields. The
  final valid value of the telemetry point  corresponds to a electronics
  within the middle L1A granule is considered to be the  status of the
  electronics in that granule. The corresponding bit of the attribute 
  "Electronics Configuration Status" is set to be 1 if the final valid
  value is 1 and  0 if that value is 0. If there is any change among the
  valid values of one of the  telemetry points within the granule or
  between the final valid value of the leading  granule and that of the
  middle granule, the corresponding bit of the attribute  "Electronics
  Configuration Change" is set to 1. If there is no valid value for some 
  of the electronics telemetry fields in the middle granule, the most
  significant bits  of the second words of both attributes are set to 1 as
  a flag for the analysts. If  there is no valid value for some of the
  electronics telemetry fields in the leading  granule, the comparision of
  the final values of those fields in the leading granule  and those in the
  middle granule is not made.
!END********************************************************************
*/

{
  PGSt_SMF_status returnStatus;
  char *location = "Get_Elec_Config_Status";
  intn            hdf_return;
  int32    i;
  int32    Leading_v_id;      /* vdata interface for leading L1A granule */
  int32    Leading_sd_id;     /* sds interface for leading L1A granule */ 
  int32    num_scans_leading; /* number of scans in leading L1A granule */
  uint32   Leading_Elec_config_status[E_VECTOR_SIZE]; 
                              /* The Electronics Configuration Status
                                 in leading L1A granule */
  uint32   Leading_Elec_config_change[E_VECTOR_SIZE]; 
                              /* The Electronics Configuration Change
                                 in leading L1A granule */
  uint32   Leading_Elec_config_invalid_flag[E_VECTOR_SIZE]; 
                              /* identifies if there are valid values
                                 for the telemetry points within 
                                 leading L1A granule */
  uint32   Elec_config_invalid_flag[E_VECTOR_SIZE];         
                              /* identifies if there are valid values
                                 for the telemetry points within 
                                 middle L1A granule */
  PGSt_integer       Version      = 1;
  char     file_name[PGSd_PC_FILE_PATH_MAX];  
                              /* leading granule file name */

  /* Check if the input parameter num_scans has valid value */
 
  if (num_scans <= 0 || num_scans > MAX_NUM_SCANS)
  {
    returnStatus = MODIS_F_INVALID_ARGUMENT;
    L1BErrorMsg(location, returnStatus, "number of scans is out of range",
                NULL, 0, NULL, True);
    return returnStatus;
  }


  /*
   * Check if the leading granule is missing. If it is missing, only the
   * Electronics  Configuration Status and the Electronics Configuration
   * Change for the middle granule need to be calculated.
   */

  
  if (QA_common->missing_leading_granule == True)
  {
    returnStatus = Get_Elec_Config_Status_Per_Gran(v_id, 
                                                   num_scans, 
                                                   Elec_config_status,
                                                   Elec_config_change, 
                                                   Elec_config_invalid_flag);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Get_Elec_Config_Status_Per_Gran", 
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }
  }

  else
  {
    /* Get the file name of leading granule */

    returnStatus = PGS_PC_GetReference (LEADING_L1A_GRANULE, 
                                        &Version, 
                                        file_name);
    if (returnStatus != PGS_S_SUCCESS)
    {
      returnStatus = MODIS_F_FILE_NOT_FOUND;
      L1BErrorMsg(location, returnStatus, 
                  "Could not retrieve file name from PCF.",
                  "PGS_PC_GetReference", LEADING_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    Leading_sd_id = SDstart(file_name, DFACC_RDONLY); /*for sds interfaces*/
    if (Leading_sd_id == FAIL)
    {
      returnStatus = MODIS_F_FILE_NOT_OPENED;
      L1BErrorMsg(location, returnStatus, 
                  "Could not open file for SD read access.",
                  "SDstart", 
                  LEADING_L1A_GRANULE,
                  "The file may be missing, corrupted or not an HDF-4 file.", 
                  True);
      return returnStatus;
    }

    Leading_v_id = Hopen(file_name,DFACC_RDONLY,0); /* for vdate interface */
    if (Leading_v_id == FAIL)                       
    {
      returnStatus = MODIS_F_FILE_NOT_OPENED;
      L1BErrorMsg(location, returnStatus, 
                  "Could not open file for Vdata read access.",
                  "Hopen", 
                  LEADING_L1A_GRANULE,
                  "The file may be corrupted or not an HDF-4 file.", 
                  True);
      return returnStatus;
    }

    hdf_return = Vstart(Leading_v_id); 
    if (hdf_return == FAIL)
    {
      returnStatus = MODIS_F_HDF_ERROR;
      L1BErrorMsg(location, returnStatus, 
                  "Could not initialize Vdata interface.",
                  "Vstart", 
                  LEADING_L1A_GRANULE,
                  "The file may be corrupted or not an HDF-4 file.", 
                  True);
      return returnStatus;
    }

    returnStatus = read_attribute(Leading_sd_id, "Number of Scans",
                                   DFNT_INT32, (void *)&num_scans_leading);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, 
                  "Could not read Number of Scans.",
                  "read_attribute", 
                  LEADING_L1A_GRANULE, 
                  Invalid_MOD01_Msg, 
                  True);
      return returnStatus;
    }

    /* 
     * Determine the Electronics Configuration Status and the Electronics Configuration
     * Change within the leading granule. Also determine if there is valid data for each
     * of the telemetry fields compose the Electronics Configuration in the leading granule.
     */
  
    returnStatus = Get_Elec_Config_Status_Per_Gran
                                          (Leading_v_id, 
                                           num_scans_leading, 
                                           Leading_Elec_config_status,
                                           Leading_Elec_config_change, 
                                           Leading_Elec_config_invalid_flag);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Get_Elec_Config_Status_Per_Gran", 
                  LEADING_L1A_GRANULE, NULL, True);
      return returnStatus;
    }

    /* Do the same thing as above for the middle granule */

    returnStatus = Get_Elec_Config_Status_Per_Gran
                                          (v_id,
                                           num_scans,
                                           Elec_config_status,
                                           Elec_config_change,
                                           Elec_config_invalid_flag);
    if (returnStatus != MODIS_S_OK)
    {
      L1BErrorMsg(location, returnStatus, NULL,
                  "Get_Elec_Config_Status_Per_Gran", 
                  FIRST_L1A_GRANULE, NULL, True);
      return returnStatus;
    }


  /*
   * If a telemetry field value changed within the middle granule, set the
   * corresponding bit to 1. If there is a change between the last valid
   * value of a telemetry field in the leading granule and that in the
   * middle granule, also set the corresponding bit  to 1.
   */

  
    for (i = 0; i < E_VECTOR_SIZE; i++)
       Elec_config_change[i] = Elec_config_change[i] |
         ((~Leading_Elec_config_invalid_flag[i])&
         (Leading_Elec_config_status[i]^Elec_config_status[i]));
   
     /* Close the leading L1A file */
  
    hdf_return = SDend(Leading_sd_id);
    if (hdf_return == FAIL)
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, 
                   "SDend", LEADING_L1A_GRANULE,
                   "Memory or the disk file must have become corrupted.", 
                   True);
  
    hdf_return = Vend(Leading_v_id);
    if (hdf_return == FAIL)
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, 
                   "Vend", LEADING_L1A_GRANULE,
                   "Memory or the disk file must have become corrupted.", 
                   True);

    hdf_return = Hclose(Leading_v_id);
    if (hdf_return == FAIL)
       L1BErrorMsg(location, MODIS_F_HDF_ERROR, NULL, 
                   "Hclose", LEADING_L1A_GRANULE,
                   "Memory or the disk file must have become corrupted.", 
                   True);
  }


  /*
   * If any of the telemetry fields has no valid value within the middle
   * granule, set the most significant bit of the second words of the
   * Electronics Configuration Status and the Electronics Configuration
   * Change to be 1.
   */
 

  if (Elec_config_invalid_flag[0] > 0 || Elec_config_invalid_flag[1] > 0)
  {
    Elec_config_change[1] |= 0x80000000;
    Elec_config_status[1] |= 0x80000000;
  }

  return returnStatus;
} 

/* end of Metadata.c */    

