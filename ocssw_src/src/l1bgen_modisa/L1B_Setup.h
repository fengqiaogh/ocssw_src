#ifndef L1B_SETUP_H
#define L1B_SETUP_H

/*
!C-INC***********************************************************************
!Description:  Public header file containing macros, structure definitions 
               and external function declarations.

!Revision History:
   $Log: L1B_Setup.h,v $
   Revision 1.7  2006-10-27 11:50:35-04  ltan
   Changed for ANSI-C compliance. Correction for the generation of code change log.

   Revision 01.03  November 19, 2001   Razor Issue #169
   Added skip_night_hi_res to argument lists
   Alice Isaacman  SAIC GSO (Alice.R.Isaacman.1@gsfc.nasa.gov)

   Revision 01.02  October 29, 2000
   Changes made for issue 142, Added external function
   "Determine_Other_Missing_Scans".
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.01 Feb 10, 1999
   Moved declaration of Write_L1B_ScanMeta into here to make as external
   function called by main().
   Jim Rogers (rogers@mcst.gsfc.nasa.gov)

   Revision 01.00 Feb 2, 1999
   Initial development
   This was formed from declarations in L1B.c as part of the
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

#include    "L1B_Tables.h"
#include    "Granule.h"

/*
 * L1B_Scan_Metadata_t: holds data for writing the scan metadata to the
 * L1B output granules.  Most of the members are assigned/filled in
 * L1B_Setup().  Additional member(s) are filled in Emissive_Cal during
 * the loop through scans.  After the loop through scans, these data are
 * written by "Write_L1B_ScanMeta".
 */

typedef struct {
  int32     v_id[NUM_L1B_EV_FILES];
  int32     sd_id[NUM_L1B_EV_FILES];
  int32     num_scans;
  int32     Complete_Scan_Flag;
  int32     EV_frames;
  float32   Nadir_Frame_Latitude[MAX_NUM_SCANS];
  float32   Nadir_Frame_Longitude[MAX_NUM_SCANS];
  int32     Nadir_Frame_Number;
  float32   Nadir_Frame_Solar_Azimuth[MAX_NUM_SCANS];
  float32   Nadir_Frame_Solar_Zenith[MAX_NUM_SCANS];
  int16     MirrorSide[MAX_NUM_SCANS];
  char      ScanType[MAX_NUM_SCANS][SCAN_TYPE_TEXT_SIZE];
  float64   EVStartTime_TAIsecond[MAX_NUM_SCANS];
  uint32    Bit_QA_Flags[MAX_NUM_SCANS];
} L1B_Scan_Metadata_t;

/*
 * Function prototypes.
 */

  PGSt_SMF_status Determine_Other_Missing_Scans (lookup_tables_t *tables,
                                                 L1A_granule_t *L1A_Gran);

  PGSt_SMF_status  L1B_Setup (lookup_tables_t     *lookup_tables,
                              L1A_granule_t       *L1A_Gran,
                              L1B_granule_t       *L1B_Gran,
                              L1A_Scan_t          *L1A_Scan,
                              L1B_Scan_t          *L1B_Scan,
                              QA_Data_t           *QA,
                              L1B_Scan_Metadata_t *L1B_Scan_Meta,
                              boolean             skip_night_hi_res);

  PGSt_SMF_status  Write_L1B_ScanMeta 
                             (L1B_Scan_Metadata_t *L1B_Scan_Meta,
                              L1A_granule_t       *L1A_Gran,
                              QA_Data_t           *QA,
                              boolean             skip_night_hi_res);

#endif




