#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_TD.h"
#include "MD_metadata.h"
#include "SC_scan.h"


void  create_missing_scans (int16           prev_scan_num,
                            PGSt_double     scan_rate,
                            PGSt_double     *SD_start_time,
                            MD_SCAN_MET_t   *scan_metadata,
                            SC_PIXEL_QUALITY_DATA_t   *pixel_qual_data)

/*
!C**********************************************************************

!Description:  This function defaults the scan metadata to create a missing
               scan.  Only the scan metadata will be defaulted as the rest
               of the scan is already set with fill values.

!Input Parameters:
               int16  prev_scan_num            **  previous scan number    **
               PGSt_double  scan_rate          **  scan rate               **

!Output Parameters:
               MD_SCAN_MET_t  scan_metadata    **  Defaulted scan metadata **

!Input/Output Parameters:
               PGSt_double SD_start_time       **  SD start time           **

Return Values: 
               None

Externally Defined:  
               MD_SCAN_META_t          (MD_metadata.h)
               PGSt_double             (PGS_TD.h)

Called By:
               handle_missing_scans

Routines Called:
               initialize_scan_metadata

!Revision History:
               Revision 1.0  1997/08/28  17:30
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               None

!END***************************************************************************
*/
{
  /******************************************************************************/
  /*                   Declare Global Variables                                 */
  /******************************************************************************/

  extern PGSt_double  global_time_offset_array[5];  
					     /*  array containing the time      */
					     /*   offsets for each sector from  */
					     /*   the beginning of the scan     */
					     /*   (ie. SD sector)               */


  /******************************************************************************/
  /*                                                                            */
  /* CALL initialize_scan_metadata to initialize the scan metadata              */
  /*   INPUT:  None                                                             */
  /*   OUTPUT: MD_SCAN_META_t                                                   */
  /*   RETURN: None                                                             */
  /*                                                                            */
  /******************************************************************************/

  initialize_scan_metadata(scan_metadata);


  /******************************************************************************/
  /*                                                                            */
  /* Set MD_SCAN_META_t.scan_num to prev_scan_num + 1                           */
  /*                                                                            */
  /* Set SD_start_time to SD_start_time + scan_rate                             */
  /*                                                                            */
  /* Set MD_SCAN_META_t.sd_start_time to SD_start_time                          */
  /* Set MD_SCAN_META_t.srca_start_time to SD_start_time +                      */
  /*    global_time_offset_array[1]                                             */
  /* Set MD_SCAN_META_t.bb_start_time to SD_start_time +                        */
  /*    global_time_offset_array[2]                                             */
  /* Set MD_SCAN_META_t.sv_start_time to SD_start_time +                        */
  /*    global_time_offset_array[3]                                             */
  /* Set MD_SCAN_META_t.ev_start_time to SD_start_time +                        */
  /*    global_time_offset_array[4]                                             */
  /*                                                                            */
  /******************************************************************************/

  scan_metadata->scan_num = prev_scan_num + 1;

  *SD_start_time = *SD_start_time + scan_rate;

  scan_metadata->sd_start_time = *SD_start_time;
  scan_metadata->srca_start_time = *SD_start_time + global_time_offset_array[1];
  scan_metadata->bb_start_time = *SD_start_time + global_time_offset_array[2];
  scan_metadata->sv_start_time = *SD_start_time + global_time_offset_array[3];
  scan_metadata->ev_start_time = *SD_start_time + global_time_offset_array[4];


  /******************************************************************************/
  /*                                                                            */
  /* CALL initialize_pixel_qual_data to initialize the pixel quality data       */
  /*   INPUT:  None                                                             */
  /*   OUTPUT: SC_PIXEL_QUALITY_DATA_t                                          */
  /*   RETURN: None                                                             */
  /*                                                                            */
  /******************************************************************************/

  initialize_pixel_qual_data(pixel_qual_data);


  /******************************************************************************/
  /*                                                                            */
  /* CALL finalize_pixel_qual_data to finalize the pixel quality data           */
  /*   INPUT:  MD_SCAN_META_t, SC_PIXEL_QUALITY_DATA_t                          */
  /*   OUTPUT: SC_PIXEL_QUALITY_DATA_t                                          */
  /*   RETURN: None                                                             */
  /*                                                                            */
  /******************************************************************************/

  finalize_pixel_qual_data(pixel_qual_data, scan_metadata);


  return;
}
