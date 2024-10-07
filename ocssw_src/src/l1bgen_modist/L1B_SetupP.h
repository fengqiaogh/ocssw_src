#ifndef L1B_SETUPP_H
#define L1B_SETUPP_H

/*
!C-INC***********************************************************************
!Description:  Private header file containing macros, structure definitions 
               and external function declarations.  Used only in L1B_Setup.c.

!Revision History:
   $Log: L1B_SetupP.h,v $
   Revision 1.12  2006-10-30 10:00:12-05  ltan
   Changed for ANSI-C compliance. Correction for the generation of code change log.

   NOTE:  as of April 2003, the Band 5 to Band 26 correction is also applied to
   MODIS/AQUA  (FM1) data.  This change did not involve a revision to TERRA code.
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.11, October 24, 2003  Razor Issue #196 (formerly Issue #184)
   Added parameter "tables" to "Scan_Meta_Cal" prototype.
   Liqin Tan, SAIC GSO  (ltan@saicmodis.com)

   Revision 01.10, April 15, 2002  Razor Issue #182
   Applicable to MODIS/TERRA (PFM) processing only.
   Added function prototype for Calculate_B26_B5_Correction and macro
     BAND_26_REFL_INDEX for Band 26 crosstalk correction
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.08, November 19, 2001    Razor Issue #169
   Added argument skip_night_hi_res to the following modules' argument lists:
     Copy_Geo_SDS
     Open_W_L1B_Granule
     Create_L1B_Swath 
     Open_L1B_EV_SDS
     Set_Unit_Range_Fillvalue
     Set_L1B_EV_SDS_Attrs
     Create_Band_Subsetting_SDS 
     Set_UI_ConvertToPercent_Attrs   
   Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)
   
   Revision 01.07, Dec 8, 2000, Razor issue 143
   Added argument to Init_L1B_ScaleOffset
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.06, October 29, 2000
   Moved Determine_Split_Scans and Get_Split_Scan_Indexes function
   prototypes from from L1B.c (Razor issue 142).
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.05 August 12, 1999
   Changed the prototype of Init_L1B_ScaleOffset()
   Zhenying Gu (zgu@mcst.gsfc.nasa.gov)

   Revision 01.04 Feb 22, 1999
   Changed argument list of Copy_Geo_SDS, dropping the L1B_ScanMeta.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.03 Feb 16, 1999
   Removed declarations of Write_Band_Subsetting_SDS and Write_Subsetting_SDS
   since they were superceded by Set_L1B_EV_SDS_Attrs and
   Create_Band_Subsetting_SDS (added declarations for those).
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.02 Feb 10, 1999
   Removed declaration of Write_L1B_ScanMeta since that has become and external
   function called by main(). Changed argument list of Get_SDS_id to allow
   setting of num_sds in the same function as the related array.
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.01 Feb 8, 1999
   Moved declaration of Scan_Meta_Cal into here since
   this is local in scope to L1B_Setup.c
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00 Feb 2, 1999
   Initial development
   This was formed from declarations in L1B.c and GranuleP.h as part of the
   process of separating L1B_Setup out into its own module.
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

/* Band index for Band 5 to Band 26 crosstalk correction. */
#define BAND_26_REFL_INDEX NUM_REFLECTIVE_BANDS - 1

  /* Geolocation parameters -- neal */
#define  GEO_STRIDE    5
#define  GEO_OFFSET    2

  /* Geolocation SDS ordering */

typedef enum {
  INDEX_LATITUDE,
  INDEX_LONGITUDE,
  INDEX_HEIGHT,
  INDEX_SENSORZENITH,
  INDEX_SENSORAZIMUTH,
  INDEX_RANGE,
  INDEX_SOLARZENITH,
  INDEX_SOLARAZIMUTH,
  INDEX_GFLAGS,
  NUM_GEO_SDS
} geo_sds_index_t;


  /* Geolocation SDS structure */

  /* NOTE: only the name, src_name and type are currently being used
   *       in the code.
   */

typedef struct {
  char      *name;
  char      *src_name;
  int32     type;
  char      *units;
  float32   valid_range[2]; /* not all are of type float32 */
  float32   FillValue;
  char      *line_numbers;
  char      *frame_numbers;
  float32   scale_factor;
  /* char    *long_name; */
} geolocation_sds_t;


  /*
   * The following is used in Scan_Meta_Cal:
   */
#define SOLAR_AZIMUTH_ZENITH_SCALE_FACTOR 0.01

  /*
   * The following is used in Calculate_M1O1 and Calculate_M2O2:
   */
#define  SECONDS_IN_DAY  86400

  /*
   * Function prototypes for those functions in L1B_Setup.c
   */

  PGSt_SMF_status Open_L1A_EV_SDS                                             
                  (L1A_granule_t       *L1A_Gran, 
                   L1A_Scan_t          *L1A_Scan);

  PGSt_SMF_status Calculate_Earth_Sun_Distance                                             
                  (L1A_granule_t       *L1A_Gran, 
                   float32             *Earth_Sun_Dist);

  PGSt_SMF_status Calculate_RSB_Cal_Coeff                                             
                  (lookup_tables_t     *tables, 
                   float32             E_S_Dist, 
                   RSB_Cal_Coeff_t     *RSB_Cal_Coeff);

  PGSt_SMF_status Init_L1B_ScaleOffset                                             
                  (L1B_ScaleOffset_t   *SO, 
                   RSB_Cal_Coeff_t     *RSB_Cal_Coeff, 
                   float32             E_S_Dist, 
                   lookup_tables_t     *tables);

  PGSt_SMF_status Copy_Geo_SDS                                             
                  (L1B_granule_t       *L1B_Gran, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Scan_Meta_Cal                                             
                  (lookup_tables_t     *tables,
                   L1A_granule_t       *L1A_Gran,
                   L1B_granule_t       *L1B_Gran,
                   L1B_Scan_Metadata_t *L1B_Scan_Meta,
                   QA_Data_t           *QA);

  PGSt_SMF_status Open_W_L1B_Granule                                             
                  (lookup_tables_t     *tables, 
                   L1B_granule_t       *L1B_Gran, 
                   L1B_Scan_t          *L1B_Scan, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Init_QA_Parameters                                             
                  (L1A_granule_t       *L1A_Gran, 
                   L1B_granule_t       *L1B_Gran, 
                   QA_Data_t           *QA);

  PGSt_SMF_status Create_L1B_Swath                                             
                  (L1B_granule_t       *L1B_Gran, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Open_L1B_EV_SDS                                             
                  (L1B_granule_t       *L1B_Gran, 
                   L1B_Scan_t          *L1B_Scan, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Get_SDS_id                                             
                  (int32               f,
                   L1B_Scan_t          *L1B_Scan,
                   int16               *num_sds,
                   int32               *sds_id);

  PGSt_SMF_status Set_SDS_Attributes                                             
                  (int32               *sds_id, 
                   char                **BandNames, 
                   float32             **scale, 
                   float32             **offset, 
                   char                *rad_units, 
                   char                *refl_units, 
                   char                *counts_units, 
                   int32               num_sds);

  PGSt_SMF_status Write_Swath_Band_Number                                             
                  (int32               file_index, 
                   L1B_granule_t       *L1B_Gran);

  PGSt_SMF_status Set_Unit_Range_Fillvalue                                             
                  (L1B_Scan_t          *L1B_Scan, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Set_L1B_EV_SDS_Attrs                                             
                  (lookup_tables_t     *tables, 
                   L1B_granule_t       *L1B_Gran, 
                   L1B_Scan_t          *L1B_Scan, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Create_Band_Subsetting_SDS                                             
                  (L1B_granule_t       *L1B_Gran, 
                   boolean             skip_night_hi_res);

  PGSt_SMF_status Calculate_DCR_Change                                             
                  (L1A_granule_t       *L1A_Gran, 
                   QA_Data_t           *QA, 
                   L1B_Scan_Metadata_t *L1B_Scan_Meta);

  PGSt_SMF_status Determine_Split_Scans                                             
                  (L1A_granule_t       *L1A_Gran,
                   boolean             *split_scan);

  PGSt_SMF_status Get_Split_Scan_Indexes                                             
                  (int32               S1,
                   int32               num_scans,
                   int16               mirror_side[],
                   int32               scan_quality
                                         [][SCAN_QUALITY_ARRAY_NUM_ELEMENTS],
                   int32               *S_split_1,
                   int32               *S_split_2);

  PGSt_SMF_status Set_UI_ConvertToPercent_Attrs                                             
                  (lookup_tables_t     *tables,
                   L1B_Scan_t          *L1B_Scan,
                   boolean             skip_night_hi_res);


  PGSt_SMF_status Calculate_B26_B5_Correction
                             (float32            *original_correction,
                              float32            *scaled_correction,
                              L1B_ScaleOffset_t  *ScaleOffset);

#endif
