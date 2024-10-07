#ifndef MD_metadata_H
#define MD_metadata_H

#include "PGS_PC.h"
/*
!C-INC************************************************************************

!Description:  This header file contains the metadata definitions and 
               structures for the ECS Granule Inventory Metadata (section 1.1.1
               of the MODIS Level 1A Data Product Format), MODIS Level 1A 
               Specific Granule Metadata (section 1.2 of the MODIS Level 1A 
               Data Product Format), and Scan Level Metadata (section 2 of the
               MODIS Level 1A Data Product Format).

!Input Parameters:
               N/A

!Output Parameters:
               N/A

Externally Defined:  
               int16                     (hdfi.h)
               int32                     (hdfi.h)
               PGSd_PC_VALUE_LENGTH_MAX  (PGS_PC.h)

!Revision History:
    $Log: MD_metadata.h,v $
    Revision 5.1  2004/09/30 18:54:46  seaton
    Updated to run Collection 5 L1A code.
    seaton@saicmodis.com

    Revision 4.2  2003/01/07 21:14:25  vlin
    MD_ECS_GRA_INV_MET_t struct field name "operation_mode" removed

    Revision 4.1  2002/10/30 14:31:16  vlin
    Added some fields to structure name "MD_ECS_GRA_INV_MET_t.
    vlin@saicmodis.com

               Revision 2.1 1998/10/07   13:00 EDT
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Modified #defines to reference mapiL1A.h #defines

               Revision 2.0  1997/07/15  16:32
               Qi Huang/RDC (qhuang@ltpmail.gsfc.nasa.gov)
               Copied, modified, and renamed from modis_IO.h

               Revision 1.0  1996/09/24  14:45 EDT
               Keith Degnan/SAIC/GSC (keith.degnan@gsfc.nasa.gov)
               Created the PDL for module

!Team-unique Header:

         This software is developed by the MODIS Science Data Support Team
         for the National Aeronautics and Space Administration, 
         Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:

Design Notes: 
               The ".h" file below was specifically written for development
               in C. Any other language choice may require reworking of the
               ".h" file before coding can begin.

!END
**************************************************************************/

#include "hdfi.h"
#include "mapiL1A.h"

#define   MD_TIMECODEALEN                       28
#define   MD_TIMECODEADATELEN                   11
#define   MD_TIMECODEATIMELEN                   16
#define   MD_ADDIATTRNAMELEN                    20
#define   MD_NUM_POINTERS                        3
#define   MD_POINTER_NAME_L                     50


#define   MD_NUM_SCANS_TXT                       M01NUMBER_OF_SCANS
#define   MD_NUM_DAY_SCANS_TXT                   M01NUMBER_DAY_SCANS 
#define   MD_NUM_NIGHT_SCANS_TXT                 M01NUMBER_NIGHT_SCANS 
#define   MD_MAX_TOTAL_FRAMES_TXT                M01MAX_TOTAL_FRAMES
#define   MD_MAX_EARTH_FRAMES_TXT                M01MAX_EARTH_FRAMES
#define   MD_MAX_SD_FRAMES_TXT                   M01MAX_SD_FRAMES
#define   MD_MAX_SRCA_FRAMES_TXT                 M01MAX_SRCA_FRAMES
#define   MD_MAX_BB_FRAMES_TXT                   M01MAX_BB_FRAMES
#define   MD_MAX_SV_FRAMES_TXT                   M01MAX_SV_FRAMES
#define   MD_SCAN_TYPES_TXT                      M01SCAN_TYPES
#define   MD_INCOMPLETE_SCANS_TXT                M01INCOMPL_SCANS
#define   MD_MISSING_PKTS_TXT                    M01MISSING_PACKETS
#define   MD_BAD_CRC_PKTS_TXT                    M01PACKTS_BAD_CRC
#define   MD_DISCARDED_PKTS_TXT                  M01DISCARD_PACKETS

#define   MECS_PRODHISTORY                       "PRODUCTIONHISTORY"
#define   MD_PROCESSVERSION                      "PROCESSVERSION" 
#define   MD_OTHER_STRING                        "Other"
#define   MD_MIXED_SCAN                          "Mixed"
#define   MD_DAY_SCAN                            "Day"
#define   MD_NIGHT_SCAN                          "Night"
#define   MD_BOTH                                "Both"
#define   MD_NA                                  "NA"  
#define   MD_MODIS_BOTH                          "MODIS_Both"
#define   MD_MODIS_DAY                           "MODIS_Day"
#define   MD_MODIS_NIGHT                         "MODIS_Night"
#define   MD_GRANULENUMBER                       "GRANULENUMBER"
#define   MD_MIDNIGHT                            "T00:00:00.000000Z"
#define   MD_NOT_PROCESSED                       "processed once"
#define   EASTBOUNDINGCOORDNIATE                 "EastBoundingCoordinate"
#define   WESTBOUNDINGCOORDNIATE                 "WestBoundingCoordinate"
#define   SOUTHBOUNDINGCOORDNIATE                "SouthBoundingCoordinate"
#define   NORTHBOUNDINGCOORDNIATE                "NorthBoundingCoordinate"
#define   EASTBOUNDVALUE                         -180.
#define   WESTBOUNDVALUE                          180.
#define   SOUTHBOUNDVALUE                          90.
#define   NORTHBOUNDVALUE                         -90.



#define   MD_NUM_FRAME_COUNT_ARRAY               6
#define   MD_NUM_CCSDS_APIDS                     3


/**************************************************************************/
/*  The following are the default values.                                 */
/**************************************************************************/

#define   MD_MAX_MISSING_PKTS_IN_SCAN         3032
#define   MD_NO_VALID_DATA_IN_SCAN               0
#define   MD_SOME_VALID_DATA_IN_SCAN             1


/**************************************************************************/
/*  The following are the indexes to the Frame count array.              */
/**************************************************************************/
#define   MD_TOTAL_FRAMES_IN_SCAN                0
#define   MD_EV_FRAMES_IN_SCAN                   1
#define   MD_SD_FRAMES_IN_SCAN                   2
#define   MD_SRCA_FRAMES_IN_SCAN                 3
#define   MD_BB_FRAMES_IN_SCAN                   4
#define   MD_SV_FRAMES_IN_SCAN                   5


/**************************************************************************/
/*  The following are the indexes to the Scan quality array.              */
/**************************************************************************/

#define   MD_SCAN_DATA_PRESENCE                  0
#define   MD_MISSING_PACKET                      1
#define   MD_BAD_CHECKSUM_PACKET                 2
#define   MD_DISCARDED_PACKET                    3


typedef   char MD_INPUT_POINTER[MD_NUM_POINTERS][MD_POINTER_NAME_L];

/**************************************************************************/
/*  The following structure defines the ECS Granule Inventory Metadata    */
/*    (section 1.1.1 of the MODIS Level 1A Data Product Format)           */
/**************************************************************************/

typedef struct
{
  char     *pgeversion;
  char      rangebeginningdate[MD_TIMECODEADATELEN];
  char      rangebeginningtime[MD_TIMECODEATIMELEN];
  char      rangeendingdate[MD_TIMECODEADATELEN];
  char      rangeendingtime[MD_TIMECODEATIMELEN];
  char      day_night_flag[6];
  int       orbit_num_1;
  double    equatorcrossinglongitude_1;
  char     *equatorcrossingdate_1;
  char     *equatorcrossingtime_1;
  char      exclusiongringflag_1[2];
  double    gringpointlatitude_1[4];
  double    gringpointlongitude_1[4];
  int       gringpointsequenceno_1[4];
  char     *additionalattributename_1;
  char     *additionalattributename_2;
  char      parametervalue_1[5];
  char     *parametervalue_2;
  char     *reprocessingactual;
  char     *reprocessingplanned;
  char     *processingenvironment;
  char     *localversionid;
  char      productionhistory[16];
} MD_ECS_GRA_INV_MET_t;


/**************************************************************************/
/*  The following structure defines the MODIS Level 1A Specific Granule   */
/*    Metadata (section 1.1.1 of the MODIS Level 1A Data Product Format)  */
/**************************************************************************/
  
typedef struct
{
  int32     num_scans;
  int32     num_day_scans;
  int32     num_night_scans;
  int32     max_total_frames;
  int32     max_earth_frames;
  int32     max_sd_frames;
  int32     max_srca_frames;
  int32     max_bb_frames;
  int32     max_sv_frames;
  char      scan_types_product[10];
  int32     incomplete_scans;
  int32     missing_packets;
  int32     packets_bad_crc;
  int32     discarded_packets;
} MD_L1A_SPECIFIC_MET_t;


/**************************************************************************/
/*  The following structure defines the Scan Level Metadata (section      */
/*    1.1.1 of the MODIS Level 1A Data Product Format)                    */
/**************************************************************************/

typedef struct
{
  int16     scan_num;
  int16     frame_count_array[6];
  char      scan_type[10];
  float64   sd_start_time;
  float64   srca_start_time;
  float64   bb_start_time;
  float64   sv_start_time;
  float64   ev_start_time;
  int16     srca_cal_mode;
  int16     packet_scan_count;
  int16     ccsds_apids[3];
  int16     packet_expedited_data_flag;
  int16     mirror_side;
  int32     scan_qual_array[4];
} MD_SCAN_MET_t;

#endif /* MD_metadata_H */
