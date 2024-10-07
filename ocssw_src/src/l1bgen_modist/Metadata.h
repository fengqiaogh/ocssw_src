#ifndef	METADATA_H
#define METADATA_H

#include    "Preprocess.h" /* Contain some variable definitions */
                           /*   used in Metadata.c */
#include    "L1B_Setup.h"  /* Contains L1B_Scan_Metadata_t */
#include    "PGS_MET.h"
 
/*
!C-INC**********************************************************************
!Description:  Header file Metadata.h to be included in files where to 
               populate and write granule metadata.

!Revision History:
 $Log: Metadata.h,v $
 Revision 1.9  2012-06-27 15:57:04-04  xgeng
 add the doi metadata

 Revision 1.8  2008-11-18 15:21:53-05  xgeng
 merge branch for V6.0.0

 Revision 1.7.2.3  2008/06/05 21:53:57  xgeng
 Added deadSubframeDataPercent to L1B_Gran_Metadata_t.

 Revision 1.7  2006/10/30 15:00:13  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.
 
 Revision 02.13 November 13, 2001  (Razor Issue #169)
 Added skip_night_hi_res to Write_Gran_Metadata header.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 02.12 Dec 5, 2000
 Razor issue 147
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
   
 Revision 02.11 Feb 8, 1999
 Moved declaration of Scan_Meta_Cal to L1B_SetupP.h
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)
   
 Revision 02.10 April 1998
 Changed some macros to functions, added new granule metadata,
 added new data structure  L1B_Gran_Metadata_t.
 Zhenying Gu(zgu@gscmail.gsfc.nasa.gov)
   
 Revision 01.00 1996 
 Initial development
 Zhidong Hao(hao@ltpmail.gsfc.nasa.gov)  

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   
!END********************************************************************
*/
#define E_VECTOR_SIZE 2
#define NUM_THERMISTORS    12
#define SERIAL_NUMBER_SIZE 30
#define MAX_DATE_TIME_SIZE 30
#define NUM_DCR_VALUE      550 

typedef struct {

  int32     L1A_Gran_sd_id;
  int32     L1B_Gran_sd_id[NUM_L1B_EV_FILES];

  /* For L1B metadata
  */
  int32     num_night_scans;
  int32     incomplete_scans;
  int32     max_ev_frames;
  int32     num_scans;                                         
  int32     num_day_scans;
  int32     Extract_Pixel_Offset;
  int32     Extract_Pixel_Count;
  int32     Extract_Line_Offset;
  int32     Extract_Line_Count;

  float32 validEVPercent[NUM_BANDS];
  float32 satEVPercent[NUM_BANDS];
  float32 missEVPercent[NUM_BANDS];
  uint32  elecRedVec[2];

  uint32 total_pixels[NUM_BANDS];
  uint32 valid_pixels[NUM_BANDS];
  uint32 saturated_pixels[NUM_BANDS];
  uint32 missing_pixels[NUM_BANDS];

  /*
   * Focal_Plane_Set_Point_State
   */
   
  int8    FPSetPointState;
  

  char    Reflective_LUT_Serial_Number[SERIAL_NUMBER_SIZE];                                    
  char    Emissive_LUT_Serial_Number[SERIAL_NUMBER_SIZE];               
  
  char    QA_LUT_Serial_Number[SERIAL_NUMBER_SIZE]; /* New QA */

  char    doi_attr_name[PGSd_MET_MAX_STRING_SET_L];
  char    doi_para_value[NUM_L1B_EV_FILES+1][PGSd_MET_MAX_STRING_SET_L];
  char    doi_authority_attr_name[PGSd_MET_MAX_STRING_SET_L];
  char    doi_authority_para_value[PGSd_MET_MAX_STRING_SET_L];

  /* New QA */
  int8   Door_Screen_Configuration;
  int8   Reflective_Band_Identification[NUM_REFLECTIVE_BANDS];
  int8   Emissive_Band_Identification[NUM_EMISSIVE_BANDS];
  int8   All_L1B_Error_Flag_Off;

  uint8  Thermal_Detector_Noise
             [NUM_EMISSIVE_BANDS][DETECTORS_PER_1KM_BAND];
  uint8  Thermal_Detector_Relative_Response_Change
             [NUM_EMISSIVE_BANDS][DETECTORS_PER_1KM_BAND];

  int32 qapercent_missing_250m;
  int32 qapercent_outofbound_250m;
  int32 qapercent_interpolated_250m;
  int32 qapercent_missing_500m;
  int32 qapercent_outofbound_500m;
  int32 qapercent_interpolated_500m;
  int32 refl_1km_qapercent_missing;
  int32 refl_1km_qapercent_outofbound;
  int32 qapercent_interpolated_refl_1km;
  int32 emiss_qapercent_missing;
  int32 emiss_qapercent_outofbound;
  int32 qapercent_interpolated_emiss;
  float32 Earth_Sun_Dist;
 
  /*
   * Electronics Configuration Status and 
   * Electronics Configuration Change variables 
   */
  
  uint32 Elec_config_status[E_VECTOR_SIZE];
  uint32 Elec_config_change[E_VECTOR_SIZE]; 

  float32 missAllScanDataPercent;
  float32 nightRSBPercent[NUM_DETECTORS];
  float32 missInScanDataPercent[NUM_DETECTORS];
  float32 deadDetectorDataPercent[NUM_DETECTORS];
  float32 deadSubframeDataPercent[NUM_HIGH_RESOLUTION_DETECTORS];
  float32 sectorRotateDataPercent[NUM_DETECTORS];
  float32 saturatedDataPercent[NUM_DETECTORS];
  float32 noBGDataPercent[NUM_DETECTORS];
  float32 moonInSVPTEBDataPercent[NUM_DETECTORS];
  float32 badDNStarStarRSBDataPercent[NUM_DETECTORS];
  float32 exceedMaxForScalingPercent[NUM_DETECTORS];
  float32 NADClosedDataPercent[NUM_DETECTORS];
  float32 uncalibratedDataPercent[NUM_DETECTORS];
} L1B_Gran_Metadata_t;


PGSt_SMF_status Gran_Meta_Cal (L1A_granule_t *, 
                                L1B_granule_t *, 
                                Preprocess_Data_t *, 
                                QA_Data_t *, 
                                L1B_Scan_Metadata_t *,
                                L1B_Gran_Metadata_t *);

PGSt_SMF_status Write_Gran_Metadata 
                              (Run_Time_Parameters_t *runtime_params,
                                L1B_Gran_Metadata_t   *L1B_Gran_Meta,
                                QA_Data_t             *QA,
                                Preprocess_Data_t     *PP,
                                lookup_tables_t       *tables,
                                L1A_granule_t         *L1A_Gran,
                                boolean               skip_night_hi_res);

#endif



