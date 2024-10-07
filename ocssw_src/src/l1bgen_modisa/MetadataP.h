#ifndef	METADATAP_H
#define METADATAP_H

#include    "Metadata.h"
 
/*
!C-INC**********************************************************************
!Description:  Private header file to be used in Metadata.c only.

!Revision History:
 $Log: MetadataP.h,v $
 Revision 1.7  2012-06-25 14:22:25-04  xgeng
 added doi metadata

 Revision 1.6  2006-10-27 11:01:35-04  ltan
 Changed for ANSI-C compliance. Correction for the generation of code change log.

 Revision 01.03 December 13, 2002  (Razor Issue #188)
 Added MAX_PRODUCTIONHISTORY_SIZE
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov) 

 Revision 01.02 November 19, 2001  (Razor Issue #169)
 Added skip_night_hi_res to Write_Global_Metadata header.
 Moved output_file_indices_t from Metadata.c.
 Alice Isaacman, SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

 Revision 01.01 Feb 9, 1999
 Moved define SOLAR_AZIMUTH_ZENITH_SCALE_FACTOR to L1B_SetupP.h.
 Jim Rogers (rogers@mcst.gsfc.nasa.gov)

 Revision 01.00 1998 
 part of original Metadata.h
 Initial development
 Zhenying Gu(zgu@mcst.gsfc.nasa.gov)

!References and Credits:
   This software is developed by the MODIS Characterization Support
   Team (MCST)for the National Aeronautics and Space Administration,
   Goddard Space Flight Center, under contract NAS5-32373.

   HDF portions developed at the National Center for Supercomputing
   Applications at the University of Illinois at Urbana-Champaign.

!Design Notes:
   
!END********************************************************************
*/
#define L1A_MISSING_ENG_PACKET 65535
#define MAX_PRODUCTIONHISTORY_SIZE 255

typedef struct 
{
  PGSt_PC_Logical fileID;
  PGSt_integer    version;
  char            *hdfAttr;
} pgs_meta_t;

typedef enum
{
  INDEX_L1B_EV_250M_FILE,
  INDEX_L1B_EV_500M_FILE,
  INDEX_L1B_EV_1000M_FILE,
  INDEX_L1B_OBC_FILE,
  NUM_OUTPUT_FILES
} output_file_indices_t;


PGSt_SMF_status Write_Global_Metadata
                   (L1B_Gran_Metadata_t *L1B_Gran_Meta,
                    QA_Data_t           *QA,
                    lookup_tables_t     *tables,
                    int32               OBC_sd_id,
                    boolean             skip_night_hi_res);

PGSt_SMF_status Get_Electronics_Status
                   (int32              v_id,
                    int32              num_scans,
                    char               *vname,
                    char               *fname,
                    int16              *final_value,
                    int16              *is_changed,
                    boolean            *no_valid_value);

PGSt_SMF_status Get_Elec_Config_Status_Per_Gran
                    (int32             v_id,
                     int32             num_scans,
                     uint32            *Elec_config_status,
                     uint32            *Elec_config_change,
                     uint32            *Elec_config_invalid_flag);

PGSt_SMF_status Get_Elec_Config_Status
                     (QA_Common_t      *QA_common,
                      int32            v_id,
                      int32            num_scans,
                      uint32           *Elec_config_status,
                      uint32           *Elec_config_change);

void get_attr( char *, void *);

void get_string_attr( char *, char *);

void set_attr(char *, void *);

void set_string_attr( char *, char *);

void set_ptrstring_attr( char *, char **);

void copy_string_attr( char *, char *);

void copy_attr( char *, void *);

#endif












