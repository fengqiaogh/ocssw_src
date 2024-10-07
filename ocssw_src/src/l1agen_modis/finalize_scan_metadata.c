#include "MD_metadata.h"
#include "PH_pkt_hdr.h"
#include "hdfi.h"
#include "L1A_prototype.h"


void   finalize_scan_metadata (MD_SCAN_MET_t   *scan_meta,
                               int16           num_packets)

/*
!C************************************************************************

!Description:  This function will set the "Missing packets in scan", "Total
               Frames in scan", and "XX_start_time" fields in the scan-level 
               metadata.

!Input Parameters:
               int16               num_packets    ** Number of packets in  
                                                     this scan             **

!Output Parameters:
               None

!Input/Output Parameters:
               MD_SCAN_MET_t       *scan_meta     ** Scan metadata         **

Return Values: 
               None

Externally Defined:  
               int16                                           (hdfi.h)
               MAX(a,b)                                        (hdfi.h)
               MD_SCAN_MET_t                                   (MD_metadata.h)
               MD_TOTAL_FRAMES_IN_SCAN                         (MD_metadata.h)
               MD_EV_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_SD_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_SRCA_FRAMES_IN_SCAN                          (MD_metadata.h)
               MD_BB_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_SV_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_MISSING_PACKET                               (MD_metadata.h)
               MD_NIGHT_SCAN                                   (MD_metadata.h)
               PH_REALISTIC_NUM_SD_PACKETS                     (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_SRCA_PACKETS                   (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_BB_PACKETS                     (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_SV_PACKETS                     (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_EV_DAY_PACKETS                 (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_EV_NIGHT_PACKETS               (PH_pkt_hdr.h)
               global_time_offset_array                        (level1a)

Called By:
               process_a_scan

Routines Called:
               None

!Revision History:
               Revision 2.2  1997/09/03  15:20
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code

               Revision 2.1  1997/08/27  10:35
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 2.0  1997/07/25  10:00
               Tom Johnson/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Updated for Version 2 from Version 1

               Revision 1.0  1997/06/18  16:40 EDT  
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               This routine assumes that the packet has been validated 
               before this routine is executed.  No validation of data 
               is performed in this routine.

               The following CODE was written in C language.

!END**********************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                      Define Global Variables                           */
  /*                                                                        */
  /**************************************************************************/

  extern PGSt_double  global_time_offset_array[5];  
					 /*  array containing the time      */
					 /*   offsets for each sector from  */
					 /*   the beginning of the scan     */
					 /*   (ie. SD sector)               */


  /**************************************************************************/
  /*                                                                        */
  /*           Declare the local variables and initialize them.             */
  /*                                                                        */
  /**************************************************************************/

  int16  max_SD_packets;          /* Maximum number of SD packets in scan   */
  int16  max_SRCA_packets;        /* Maximum number of SDCA packets in scan */
  int16  max_BB_packets;          /* Maximum number of BB packets in scan   */
  int16  max_SV_packets;          /* Maximum number of SV packets in scan   */
  int16  max_EV_packets;          /* Maximum number of EV packets in scan   */


  max_SD_packets   = 0;
  max_SRCA_packets = 0;
  max_BB_packets   = 0;
  max_SV_packets   = 0;
  max_EV_packets   = 0;



  /*************************************************************************/
  /*                                                                       */
  /* Set max_SD_packets to the larger of either                            */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SD_FRAMES_IN_SCAN] * 2) or    */
  /*     PH_REALISTIC_NUM_SD_PACKETS                                       */
  /*                                                                       */
  /* Set max_SRCA_packets to the larger of either                          */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SRCA_FRAMES_IN_SCAN] * 2) or  */
  /*     PH_REALISTIC_NUM_SRCA_PACKETS                                     */
  /*                                                                       */
  /* Set max_BB_packets to the larger of either                            */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_BB_FRAMES_IN_SCAN] * 2) or    */
  /*     PH_REALISTIC_NUM_BB_PACKETS                                       */
  /*                                                                       */
  /* Set max_SV_packets to the larger of either                            */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SV_FRAMES_IN_SCAN] * 2) or    */
  /*     PH_REALISTIC_NUM_SV_PACKETS                                       */
  /*                                                                       */
  /*************************************************************************/

     max_SD_packets =
                 MAX((scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN] * 2), 
                     PH_REALISTIC_NUM_SD_PACKETS);

     max_SRCA_packets = 
                 MAX((scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN] * 2), 
                     PH_REALISTIC_NUM_SRCA_PACKETS);

     max_BB_packets = 
                 MAX((scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN] * 2), 
                     PH_REALISTIC_NUM_BB_PACKETS);

     max_SV_packets =
                 MAX((scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN] * 2), 
                     PH_REALISTIC_NUM_SV_PACKETS);


  /*************************************************************************/
  /*                                                                       */
  /* IF MD_SCAN_MET_t.scan_type is equal to "Night"                        */
  /* THEN                                                                  */
  /*    Set max_EV_packets to the larger of either                         */
  /*    MD_SCAN_MET_t.frame_count_array[MD_EV_FRAMES_IN_SCAN] or           */
  /*    PH_REALISTIC_NUM_EV_NIGHT_PACKETS                                  */
  /* ELSE                                                                  */
  /*    Set max_EV_packets to the larger of either                         */
  /*    (MD_SCAN_MET_t.frame_count_array[MD_EV_FRAMES_IN_SCAN] * 2) or     */
  /*    PH_REALISTIC_NUM_EV_DAY_PACKETS                                    */
  /* ENDIF                                                                 */
  /*                                                                       */
  /*************************************************************************/

     if (strcmp(scan_meta->scan_type, MD_NIGHT_SCAN) == 0) 
        max_EV_packets = MAX(scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN],
                             PH_REALISTIC_NUM_EV_NIGHT_PACKETS);

     else 
        max_EV_packets = 
                  MAX((scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN] * 2),
                      PH_REALISTIC_NUM_EV_DAY_PACKETS);



  /*************************************************************************/
  /*                                                                       */
  /* Set MD_SCAN_MET_t.scan_qual_array[MD_MISSING_PACKET] to               */
  /*     (max_SD_packets + max_SRCA_packets + max_BB_packets +             */
  /*     max_SV_packets + max_EV_packets + PH_REALISTIC_NUM_ENG_PACKETS) - */
  /*     num_packets                                                       */
  /*                                                                       */
  /* Set MD_SCAN_MET_t.frame_count_array[MD_TOTAL_FRAMES_IN_SCAN] to       */
  /*     MD_SCAN_MET_t.frame_count_array[MD_EV_FRAMES_IN_SCAN] +           */
  /*     MD_SCAN_MET_t.frame_count_array[MD_SD_FRAMES_IN_SCAN] +           */
  /*     MD_SCAN_MET_t.frame_count_array[MD_SRCA_FRAMES_IN_SCAN] +         */
  /*     MD_SCAN_MET_t.frame_count_array[MD_BB_FRAMES_IN_SCAN] +           */
  /*     MD_SCAN_MET_t.frame_count_array[MD_SV_FRAMES_IN_SCAN]             */
  /*                                                                       */
  /*************************************************************************/

     scan_meta->scan_qual_array[MD_MISSING_PACKET] = 
                              (max_SD_packets   +
                               max_SRCA_packets +
                               max_BB_packets   + 
                               max_SV_packets   +
                               max_EV_packets   +
                               PH_REALISTIC_NUM_ENG_PACKETS) - num_packets;

     scan_meta->frame_count_array[MD_TOTAL_FRAMES_IN_SCAN] =
               scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN] +
               scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN]    +
               scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN]  +
               scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN]    +
               scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN];



  /*************************************************************************/
  /*                                                                       */
  /*           Set SD start time if it is fill value                       */
  /*                                                                       */
  /*************************************************************************/
  /*                                                                       */
  /* IF MD_SCAN_MET_t.sd_start_time is less than 0 (time is fill value)    */
  /* THEN                                                                  */
  /*    IF MD_SCAN_MET_t.srca_start_time is not less than 0                */
  /*    THEN                                                               */
  /*       Set MD_SCAN_MET_t.sd_start_time to                              */
  /*          MD_SCAN_MET_t.srca_start_time - global_time_offset_array[1]  */
  /*    ELSE                                                               */
  /*       IF MD_SCAN_MET_t.bb_start_time is not less than 0               */
  /*       THEN                                                            */
  /*          Set MD_SCAN_MET_t.sd_start_time to                           */
  /*             MD_SCAN_MET_t.bb_start_time - global_time_offset_array[2] */
  /*       ELSE                                                            */
  /*          IF MD_SCAN_MET_t.sv_start_time is not less than 0            */
  /*          THEN                                                         */
  /*             Set MD_SCAN_MET_t.sd_start_time to                        */
  /*                MD_SCAN_MET_t.sv_start_time -                          */
  /*                global_time_offset_array[3]                            */
  /*          ELSE                                                         */
  /*             Set MD_SCAN_MET_t.sd_start_time to                        */
  /*                MD_SCAN_MET_t.ev_start_time -                          */
  /*                global_time_offset_array[4]                            */
  /*          ENDIF                                                        */
  /*       ENDIF                                                           */
  /*    ENDIF                                                              */
  /* ENDIF                                                                 */
  /*                                                                       */
  /*************************************************************************/

     if (scan_meta->sd_start_time < 0.0) { 
        if (scan_meta->srca_start_time >= 0.0)
           scan_meta->sd_start_time = scan_meta->srca_start_time -
                                                 global_time_offset_array[1];
        else if (scan_meta->bb_start_time >= 0.0)
           scan_meta->sd_start_time = scan_meta->bb_start_time -
                                                 global_time_offset_array[2];
        else if (scan_meta->sv_start_time >= 0.0)
           scan_meta->sd_start_time = scan_meta->sv_start_time -
                                                 global_time_offset_array[3];
        else if (scan_meta->ev_start_time >= 0.0)
           scan_meta->sd_start_time = scan_meta->ev_start_time -
                                                 global_time_offset_array[4];
     } 
     if (scan_meta->sd_start_time >= 0.0)
     {
     /*************************************************************************/
     /*                                                                       */
     /* IF MD_SCAN_MET_t.srca_start_time is less than 0  (time is fill value) */
     /* THEN                                                                  */
     /*    Set MD_SCAN_MET_t.srca_start_time to MD_SCAN_MET_t.sd_start_time + */
     /*        global_time_offset_array[1]                                    */
     /* ENDIF                                                                 */
     /*                                                                       */
     /* IF MD_SCAN_MET_t.bb_start_time is less than 0  (time is fill value)   */
     /* THEN                                                                  */
     /*    Set MD_SCAN_MET_t.bb_start_time to MD_SCAN_MET_t.sd_start_time +   */
     /*        global_time_offset_array[2]                                    */
     /* ENDIF                                                                 */
     /*                                                                       */
     /* IF MD_SCAN_MET_t.sv_start_time is less than 0  (time is fill value)   */
     /* THEN                                                                  */
     /*    Set MD_SCAN_MET_t.sv_start_time to MD_SCAN_MET_t.sd_start_time +   */
     /*        global_time_offset_array[3]                                    */
     /* ENDIF                                                                 */
     /*                                                                       */
     /* IF MD_SCAN_MET_t.ev_start_time is less than 0  (time is fill value)   */
     /* THEN                                                                  */
     /*    Set MD_SCAN_MET_t.ev_start_time to MD_SCAN_MET_t.sd_start_time +   */
     /*        global_time_offset_array[4]                                    */
     /* ENDIF                                                                 */
     /*                                                                       */
     /*************************************************************************/

	if (scan_meta->srca_start_time < 0.0) 
	   scan_meta->srca_start_time = scan_meta->sd_start_time + 
				       global_time_offset_array[1];

	if (scan_meta->bb_start_time < 0.0) 
	   scan_meta->bb_start_time = scan_meta->sd_start_time + 
				       global_time_offset_array[2];

	if (scan_meta->sv_start_time < 0.0) 
	   scan_meta->sv_start_time = scan_meta->sd_start_time + 
				       global_time_offset_array[3];

	if (scan_meta->ev_start_time < 0.0) 
	   scan_meta->ev_start_time = scan_meta->sd_start_time + 
				       global_time_offset_array[4];

     }
 } /* End of routine finalize_scan_metadata */
