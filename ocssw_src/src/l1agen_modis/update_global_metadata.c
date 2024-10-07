#include <ctype.h>
#include "L1A_prototype.h"
#include "hdf.h"
#include "hdfi.h"
#include "MD_metadata.h"


void update_global_metadata (MD_SCAN_MET_t          *scan_meta,
                             MD_ECS_GRA_INV_MET_t   *ecs_gra_inv_met,
                             MD_L1A_SPECIFIC_MET_t  *l1a_specific_met)
/*
!C***************************************************************************

!Description:   Function update_global_metadata updates ECS Inventory Granule 
                Metadata (MD_ECS_INV_MET_t) and MODIS L1A Specific Metadata 
                (MD_L1A_SPECIFIC_MET_t).
                 
!Input Parameters:
                MD_SCAN_MET_t          scan_meta         ** Scan level metadata   **

!Output Parameters:
                None

!Input/Output Parameters:
                MD_ECS_INV_MET_t       ecs_gra_inv_met   ** ECS Inventory Granule
                                                             Metadata             **
                MD_L1A_SPECIFIC_MET_t  l1a_specific_met  ** L1A Specific Metadata **

Return Values:                
                None

Externally Defined:
                MD_SCAN_MET_t                (MD_metadata.h)
                MD_ECS_INV_MET_t             (MD_metadata.h)
                MD_L1A_SPECIFIC_MET_t        (MD_metadata.h)
                MD_INPUT_POINTER             (MD_metadata.h)
                MD_TOTAL_FRAMES_IN_SCAN      (MD_metadata.h)
                MD_EV_FRAMES_IN_SCAN         (MD_metadata.h)
                MD_SD_FRAMES_IN_SCAN         (MD_metadata.h)
                MD_SRCA_FRAMES_IN_SCAN       (MD_metadata.h)
                MD_BB_FRAMES_IN_SCAN         (MD_metadata.h)
                MD_SV_FRAMES_IN_SCAN         (MD_metadata.h)
                MD_MISSING_PACKET            (MD_metadata.h)
                MD_BAD_CHECKSUM_PACKET       (MD_metadata.h)
                MD_DISCARDED_PACKET          (MD_metadata.h)
                MD_BOTH                      (MD_metadata.h)
                MD_MODIS_BOTH                (MD_metadata.h)
                MD_MIXED_SCAN                (MD_metadata.h)
                MD_DAY_STRING                (MD_metadata.h)
                MD_NIGHT_STRING              (MD_metadata.h)
                MAX                          (hdfi.h)

Called By:
                process_a_granule
                handle_missing_scans

Routines Called:
                None

!Revision History:

                revision 1.0 1997/08/20  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support Team 
                Team (SDST) for the National Aeronautics and Space Administration
                (NASA), Goddard Space Flight Center (GSFC), under contract 
                NAS5-32373.

!References and Credits:
                None

!Design Notes:
                None

!END***************************************************************************
*/
{
  
  /****************************************************************************/
  /*                                                                          */
  /*  IF MD_SCAN_MET_t.scan_type is not equal to "Other"                      */
  /*  THEN                                                                    */
  /*     IF MD_ECS_INV_MET_t.day_night_flag is not "BOTH"                     */
  /*     THEN                                                                 */
  /*        IF MD_ECS_INV_MET_t.day_night_flag is not equal to                */
  /*           MD_SCAN_MET_t.scan_type AND MD_ECS_INV_MET_t.day_night_flag    */
  /*           is not equal to "NA"                                           */
  /*        THEN                                                              */
  /*           Set MD_ECS_INV_MET_t.day_night_flag to "Both"                  */
  /*           Set MD_ECS_INV_MET_t.operation_mode to "MODIS_Both"            */
  /*           Set MD_L1A_SPECIFIC_MET_t.scan_types_product to "Mixed"        */
  /*        ELSE                                                              */
  /*           Set MD_ECS_INV_MET_t.day_night_flag to MD_SCAN_MET_t.scan_type */
  /*           Set MD_L1A_SPECIFIC_MET_t.scan_types_product to                */
  /*              MD_SCAN_MET_t.scan_type                                     */
  /*           Set MD_ECS_INV_MET_t.operation_mode to "MODIS_" appended to    */
  /*              MD_SCAN_MET_t.scan_type                                     */
  /*        ENDIF                                                             */
  /*     ENDIF                                                                */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (strcmp(scan_meta->scan_type,MD_OTHER_STRING) != 0)
    if (strcmp(ecs_gra_inv_met->day_night_flag,MD_BOTH) != 0)
      {
        if ((strcmp(ecs_gra_inv_met->day_night_flag,scan_meta->scan_type) != 0)
            && (strcmp(ecs_gra_inv_met->day_night_flag,MD_NA) != 0))
          {
            strcpy(ecs_gra_inv_met->day_night_flag,MD_BOTH);
            memset(l1a_specific_met->scan_types_product, '\0', 
                   sizeof(l1a_specific_met->scan_types_product));
            strcpy(l1a_specific_met->scan_types_product,MD_MIXED_SCAN);
          }
        else
          {
            strcpy(ecs_gra_inv_met->day_night_flag,scan_meta->scan_type);
            memset(l1a_specific_met->scan_types_product, '\0', 
                   sizeof(l1a_specific_met->scan_types_product));
            strcpy(l1a_specific_met->scan_types_product,scan_meta->scan_type);
          }
      }


  /****************************************************************************/
  /*                                                                          */
  /*  Increment MD_L1A_SPECIFIC_MET_t.num_scans by 1                          */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->num_scans ++;


  /****************************************************************************/
  /*                                                                          */
  /*  IF MD_SCAN_MET_t.scan_type is "DAY"                                     */
  /*  THEN                                                                    */
  /*     Increment MD_L1A_SPECIFIC_MET_t.num_day_scans by 1                   */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (strcmp(scan_meta->scan_type,MD_DAY_SCAN) == 0)
    l1a_specific_met->num_day_scans ++;


  /****************************************************************************/
  /*                                                                          */
  /*  IF MD_SCAN_MET_t.scan_type is "NIGHT"                                   */
  /*  THEN                                                                    */
  /*    Increment MD_L1A_SPECIFIC_MET_t.num_night_scans by 1                  */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (strcmp(scan_meta->scan_type,MD_NIGHT_SCAN) == 0)
    l1a_specific_met->num_night_scans ++;


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_total_frames to the larger of either the  */
  /*     previous MD_L1A_SPECIFIC_MET_t.max_total_frames or                   */
  /*     MD_SCAN_MET_t.frame_count_array[MD_TOTAL_FRAMES_IN_SCAN]             */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_total_frames = MAX(l1a_specific_met->max_total_frames,
      scan_meta->frame_count_array[MD_TOTAL_FRAMES_IN_SCAN]);


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_earth_frames to the larger of either the  */
  /*     previous MD_L1A_SPECIFIC_MET_t.max_earth_frames or                   */
  /*     MD_SCAN_MET_t.frame_count_array[MD_EV_FRAMES_IN_SCAN]                */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_earth_frames = MAX(l1a_specific_met->max_earth_frames,
      scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN]);


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_sd_frames to the larger of either the     */
  /*     previous MD_L1A_SPECIFIC_MET_t.max_sd_frames or                      */
  /*     MD_SCAN_MET_t.frame_count_array[MD_SD_FRAMES_IN_SCAN]                */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_sd_frames = MAX(l1a_specific_met->max_sd_frames,
      scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN]);


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_srca_frames to the larger of either the   */ 
  /*     previous MD_L1A_SPECIFIC_MET_t.max_srca_frames or                    */
  /*     MD_SCAN_MET_t.frame_count_array[MD_SRCA_FRAMES_IN_SCAN]              */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_srca_frames = MAX(l1a_specific_met->max_srca_frames,
      scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN]);


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_bb_frames to the larger of either the     */
  /*     previous MD_L1A_SPECIFIC_MET_t.max_bb_frames or                      */
  /*     MD_SCAN_MET_t.frame_count_array[MD_BB_FRAMES_IN_SCAN]                */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_bb_frames = MAX(l1a_specific_met->max_bb_frames,
      scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN]);


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.max_sv_frames to the larger of either the     */ 
  /*     previous MD_L1A_SPECIFIC_MET_t.max_sv_frames or                      */ 
  /*     MD_SCAN_MET_t.frame_count_array[MD_SV_FRAMES_IN_SCAN]                */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->max_sv_frames = MAX(l1a_specific_met->max_sv_frames,
      scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN]);


  /*******************************************************************************/
  /*                                                                             */
  /*  IF (MD_SCAN_MET_t.scan_qual_array[MD_MISSING_PACKET] is greater than 0) OR */
  /*     (MD_SCAN_MET_t.scan_qual_array[MD_BAD_CHECKSUM] is greater than 0) OR   */
  /*     (MD_SCAN_MET_t.scan_qual_array[MD_DISCARDED_PACKET] is greater than 0)  */
  /*  THEN                                                                       */
  /*     Set MD_L1A_SPECIFIC_MET_t.incomplete_scans to                           */
  /*        MD_L1A_SPECIFIC_MET_t.incomplete_scans + 1                           */
  /*  ENDIF                                                                      */
  /*                                                                             */
  /*******************************************************************************/

  if ( (scan_meta->scan_qual_array[MD_MISSING_PACKET] > 0) ||
       (scan_meta->scan_qual_array[MD_BAD_CHECKSUM_PACKET] > 0) ||
       (scan_meta->scan_qual_array[MD_DISCARDED_PACKET] > 0) )
    l1a_specific_met->incomplete_scans ++; 


  /****************************************************************************/
  /*                                                                          */
  /*  IF MD_L1A_SPECIFIC_MET_t.num_scans is equal to 1                        */
  /*  THEN                                                                    */
  /*     Set MD_L1A_SPECIFIC_MET_t.missing_packets to                         */
  /*        MD_SCAN_MET_t.scan_qual_array[MD_MISSING_PACKET]                  */
  /*  ELSE                                                                    */
  /*     Set MD_L1A_SPECIFIC_MET_t.missing_packets to                         */
  /*        MD_L1A_SPECIFIC_MET_t.missing_packets +                           */
  /*        MD_SCAN_MET_t.scan_qual_array[MD_MISSING_PACKET]                  */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (l1a_specific_met->num_scans == 1)
    l1a_specific_met->missing_packets = scan_meta->scan_qual_array[MD_MISSING_PACKET];
  else
    l1a_specific_met->missing_packets += scan_meta->scan_qual_array[MD_MISSING_PACKET];


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.packets_bad_crc to                            */
  /*     MD_L1A_SPECIFIC_MET_t.packets_bad_crc +                              */
  /*     MD_SCAN_MET_t.scan_qual_array[MD_BAD_CHECKSUM]                       */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->packets_bad_crc += scan_meta->scan_qual_array[MD_BAD_CHECKSUM_PACKET];


  /****************************************************************************/
  /*                                                                          */
  /*  Set MD_L1A_SPECIFIC_MET_t.discarded_packets to                          */
  /*     MD_L1A_SPECIFIC_MET_t.discarded_packets +                            */
  /*     MD_SCAN_MET_t.scan_qual_array[MD_DISCARDED_PACKET]                   */
  /*                                                                          */
  /****************************************************************************/

  l1a_specific_met->discarded_packets += scan_meta->scan_qual_array[MD_DISCARDED_PACKET];


  /****************************************************************************/
  /*RETURN                                                                    */
  /****************************************************************************/

  return;

}
