#include "MD_metadata.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"
#include "hdfi.h"
#include "L1A_prototype.h"


void   finalize_pixel_qual_data (SC_PIXEL_QUALITY_DATA_t   *scan_pixel,
                                 MD_SCAN_MET_t             *scan_meta)

/*
!C************************************************************************

!Description:  This function will set the unused parts of the pixel quality 
               data to fill values.

!Input Parameters:
               

!Output Parameters:
               MD_SCAN_MET_t            *scan_meta   ** The MODIS scan level **
	                                             ** metadata structure   **

!Input/Output Parameters:
               SC_PIXEL_QUALITY_DATA_t  *scan_pixel  ** The MODIS scan pixel **
                                                     ** quality structure    **

Return Values: 
               None

Externally Defined:  
               SC_PIXEL_QUALITY_DATA_t                    (SC_scan.h)
               SC_FILL_VALUE                              (SC_scan.h)
               MD_SCAN_MET_t                              (MD_metadata.h)
               MD_EARTH_FRAMES_IN_SCAN                    (MD_metadata.h)
               MD_SD_FRAMES_IN_SCAN                       (MD_metadata.h)
               MD_SRCA_FRAMES_IN_SCAN                     (MD_metadata.h)
               MD_BB_FRAMES_IN_SCAN                       (MD_metadata.h)
               MD_SV_FRAMES_IN_SCAN                       (MD_metadata.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT     (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_SD_FRAMES                 (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_SRCA_FRAMES               (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_BB_FRAMES                 (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_SV_FRAMES                 (PH_pkt_hdr.h)
               PH_REALISTIC_NUM_EV_FRAMES                 (PH_pkt_hdr.h)
               MAX                                        (hdfi.h)

Called By:
               process_a_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/05  14:00
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code

               Revision 1.0  1997/08/27  11:00
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Created because of PDL walkthru comments

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The following CODE was written in C language.

!END**********************************************************************
*/

 {
  /*************************************************************************/
  /*                                                                       */
  /*          Declare the local variables and initialize them.             */
  /*                                                                       */
  /*************************************************************************/

  int  max_SD_frames;            /* Maximum number of SD packets in scan   */
  int  max_SRCA_frames;          /* Maximum number of SDCA packets in scan */
  int  max_BB_frames;            /* Maximum number of BB packets in scan   */
  int  max_SV_frames;            /* Maximum number of SV packets in scan   */
  int  max_EV_frames;            /* Maximum number of EV packets in scan   */
  int  frame;                    /* Scan Pixel frame index counter         */


  max_SD_frames   = 0;
  max_SRCA_frames = 0;
  max_BB_frames   = 0;
  max_SV_frames   = 0;
  max_EV_frames   = 0;



  /*************************************************************************/
  /*                                                                       */
  /* Set max_SD_frames to the larger of either                             */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SD_FRAMES_IN_SCAN] or         */
  /*     PH_REALISTIC_NUM_SD_FRAMES                                        */
  /*                                                                       */
  /* FOR frame equals max_SD_frames to                                     */
  /*                           (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1)    */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SD_pix_qual[frame][0] to SC_FILL_VALUE */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SD_pix_qual[frame][1] to SC_FILL_VALUE */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     max_SD_frames = MAX(scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN],
                         PH_REALISTIC_NUM_SD_FRAMES);

     for (frame = max_SD_frames;
          frame <= (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1);
          frame++) {
        scan_pixel->SD_pix_qual[frame][0] = SC_FILL_VALUE;
        scan_pixel->SD_pix_qual[frame][1] = SC_FILL_VALUE;
     }



  /*************************************************************************/
  /*                                                                       */
  /* Set max_SRCA_frames to the larger of either                           */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SRCA_FRAMES_IN_SCAN] or       */
  /*     PH_REALISTIC_NUM_SRCA_FRAMES                                      */
  /*                                                                       */
  /* FOR frame equals max_SRCA_frames to                                   */
  /*                           (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1)    */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SRCA_pix_qual[frame][0] to             */
  /*        SC_FILL_VALUE                                                  */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SRCA_pix_qual[frame][1] to             */
  /*        SC_FILL_VALUE                                                  */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     max_SRCA_frames = MAX(scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN],
                           PH_REALISTIC_NUM_SRCA_FRAMES);

     for (frame = max_SRCA_frames;
          frame <= (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1);
          frame++) {
        scan_pixel->SRCA_pix_qual[frame][0] = SC_FILL_VALUE;
        scan_pixel->SRCA_pix_qual[frame][1] = SC_FILL_VALUE;
     }



  /*************************************************************************/
  /*                                                                       */
  /* Set max_BB_frames to the larger of either                             */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_BB_FRAMES_IN_SCAN] or         */
  /*     PH_REALISTIC_NUM_BB_FRAMES                                        */
  /*                                                                       */
  /* FOR frame equals max_BB_frames to                                     */
  /*                           (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1)    */
  /*    Set SC_PIXEL_QUALITY_DATA_t.BB_pix_qual[frame][0] to SC_FILL_VALUE */
  /*    Set SC_PIXEL_QUALITY_DATA_t.BB_pix_qual[frame][1] to SC_FILL_VALUE */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     max_BB_frames = MAX(scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN],
                         PH_REALISTIC_NUM_BB_FRAMES);

     for (frame = max_BB_frames;
          frame <= (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1);
          frame++) {
        scan_pixel->BB_pix_qual[frame][0] = SC_FILL_VALUE;
        scan_pixel->BB_pix_qual[frame][1] = SC_FILL_VALUE;
     }



  /*************************************************************************/
  /*                                                                       */
  /* Set max_SV_frames to the larger of either                             */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_SV_FRAMES_IN_SCAN] or         */
  /*     PH_REALISTIC_NUM_SV_FRAMES                                        */
  /*                                                                       */
  /* FOR frame equals max_SV_frames to                                     */
  /*                           (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1)    */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SV_pix_qual[frame][0] to SC_FILL_VALUE */
  /*    Set SC_PIXEL_QUALITY_DATA_t.SV_pix_qual[frame][1] to SC_FILL_VALUE */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     max_SV_frames = MAX(scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN],
                         PH_REALISTIC_NUM_SV_FRAMES);

     for (frame = max_SV_frames;
          frame <= (PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX - 1);
          frame++) {
        scan_pixel->SV_pix_qual[frame][0] = SC_FILL_VALUE;
        scan_pixel->SV_pix_qual[frame][1] = SC_FILL_VALUE;
     }



  /*************************************************************************/
  /*                                                                       */
  /* Set max_EV_frames to the larger of either                             */
  /*     (MD_SCAN_MET_t.frame_count_array[MD_EV_FRAMES_IN_SCAN] or         */
  /*     PH_REALISTIC_NUM_EV_FRAMES                                        */
  /*                                                                       */
  /* FOR frame equals max_EV_frames to                                     */
  /*                         (PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT - 1)  */
  /*    Set SC_PIXEL_QUALITY_DATA_t.EV_pix_qual[frame][0] to SC_FILL_VALUE */
  /*    Set SC_PIXEL_QUALITY_DATA_t.EV_pix_qual[frame][1] to SC_FILL_VALUE */
  /* ENDFOR                                                                */
  /*                                                                       */
  /*************************************************************************/

     max_EV_frames = MAX(scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN],
                         PH_REALISTIC_NUM_EV_FRAMES);

     for (frame = max_EV_frames; 
          frame <= (PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT - 1);
          frame++) {
        scan_pixel->EV_pix_qual[frame][0] = SC_FILL_VALUE;
        scan_pixel->EV_pix_qual[frame][1] = SC_FILL_VALUE;
     }

 }  /* End of routine finalize_pixel_qual_data */


