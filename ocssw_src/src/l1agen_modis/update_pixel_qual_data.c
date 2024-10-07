#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"
#include "hdfi.h"
#include "L1A_prototype.h"


void   update_pixel_qual_data (PH_PACKET_HEADER_t       pkt_header,
                               int16                    qual_value,           
                               SC_PIXEL_QUALITY_DATA_t  *scan_pix)

/*
!C************************************************************************

!Description:  This function updates the scan pixel quality arrays, by 
               using the qual_value as the value to be inserted into the 
               pixel quality array locations.  

!Input Parameters:
               PH_PACKET_HEADER_t  pkt_header     ** Current L0 data packet **
                                                  ** packet header          **

               int16               qual_value     ** Pixel quality value    **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_PIXEL_QUALITY_DATA_t  *scan_pix ** The MODIS scan pixel   **
                                                  ** quality structure      **

Return Values: 
               None 

!Externally Defined:  
               int16                                      (hdfi.h)
               SC_PIXEL_QUALITY_DATA_t                    (SC_scan.h)
               PD_DN_FIRST_IFOV_DAY_PKT_1                 (PD_pkt_data.h)
               PD_DN_LAST_IFOV_DAY_PKT_2                  (PD_pkt_data.h)
               PD_DN_FIRST_IFOV_DAY_PKT_1                 (PD_pkt_data.h)
               PD_DN_LAST_IFOV_DAY_PKT_2                  (PD_pkt_data.h)
               PD_DN_NUM_IFOVS_IN_NIGHT_PKT               (PD_pkt_data.h)
               PH_PRI_LONG_PKT_LENGTH                     (PH_pkt_hdr.h) 
               PH_PRI_SHORT_PKT_LENGTH                    (PH_pkt_hdr.h) 
               PH_SEC_PKT_TYPE_DAY_GROUP                  (PH_pkt_hdr.h) 
               PH_SEC_PKT_TYPE_NIGHT_GROUP                (PH_pkt_hdr.h) 
               PH_SEC_PKT_TYPE_ENG1_GROUP                 (PH_pkt_hdr.h) 
               PH_SEC_PKT_TYPE_ENG2_GROUP                 (PH_pkt_hdr.h) 
               PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP         (PH_pkt_hdr.h) 
               PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP        (PH_pkt_hdr.h) 
               PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE  (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE (PH_pkt_hdr.h)  
               PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE     (PH_pkt_hdr.h) 
               PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH           (PH_pkt_hdr.h) 
               PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE
                                                          (PH_pkt_hdr.h)
Called By:
               process_a_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/02  13:43 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments in CODE.

               Revision 2.1  1997/08/27  10:30
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 2.0  1997/08/11  11:40 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design. 

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
               The CODE below was developed in C language.

               This routine assumes that the data contained in the
               input parameters have been validated before this
               routine is executed.  No validation of data is performed
               in this routine.

!END**********************************************************************
*/

 {

  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */ 
  /**************************************************************************/

  int16      sequence;             /* Packet sequence location              */
  int16      frame;                /* Calculated frame number; used to      */
                                   /* properly place data in science arrays */


  frame    = 0;

  
  /**************************************************************************/
  /*                                                                        */
  /* IF ((PH_PACKET_HEADER_t.pkt_type is not equal to                       */
  /*                                    PH_SEC_PKT_TYPE_ENG1_GROUP)         */
  /*    AND (PH_PACKET_HEADER_t.pkt_type is not equal to                    */
  /*                                    PH_SEC_PKT_TYPE_ENG2_GROUP) )       */
  /* THEN                                                                   */
  /*                                                                        */
  /**************************************************************************/

  if ((pkt_header.pkt_type != PH_SEC_PKT_TYPE_ENG1_GROUP) &&
      (pkt_header.pkt_type != PH_SEC_PKT_TYPE_ENG2_GROUP))
    {


  /**************************************************************************/
  /*                                                                        */
  /* Calculate the array indices for the pixel quality array, based on      */
  /* this packet's location information.                                    */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /*    IF (PH_PACKET_HEADER_t.pkt_type equals PH_SEC_PKT_TYPE_DAY_GROUP)   */
  /*    THEN                                                                */
  /*       set sequence equal to 0                                          */
  /*       IF (PH_PACKET_HEADER_t.sequence_flag equals                      */
  /*          PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP)                          */
  /*       THEN                                                             */
  /*          set sequence equal to 1                                       */
  /*       ENDIF                                                            */
  /*                                                                        */
  /**************************************************************************/

     if (pkt_header.pkt_type == PH_SEC_PKT_TYPE_DAY_GROUP) {
        sequence = 0;
        if (pkt_header.sequence_flag == PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP)
           sequence = 1;



  /**************************************************************************/
  /* Determine the current packet's frame count and write the qual_value    */
  /* to the appropriate pixel quality array                                 */
  /*                                                                        */
  /**************************************************************************/
  /*       IF (PH_PACKET_HEADER_t.source_ID_type_flag equals                */
  /*          PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH)                             */
  /*       THEN                                                             */
  /*          IF PH_PACKET_HEADER_t.earth_frame_cnt is greater than 0       */
  /*          THEN                                                          */
  /*             set frame equal to PH_PACKET_HEADER_t.earth_frame_cnt - 1  */
  /*             set SC_PIXEL_QUALITY_DATA_t.EV_pix_qual[frame][sequence]   */
  /*                to qual_value                                           */
  /*          ENDIF                                                         */
  /*                                                                        */
  /**************************************************************************/

        if (pkt_header.source_ID_type_flag == PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH){
           if (pkt_header.earth_frame_cnt > 0) {
              frame = pkt_header.earth_frame_cnt - 1 ;
              scan_pix->EV_pix_qual[frame][sequence] = qual_value;
           }
        }
        


  /**************************************************************************/
  /*       ELSE        ** Must be a calibration frame **                    */
  /*          IF PH_PACKET_HEADER_t.cal_frame_cnt is greater than 0         */
  /*          THEN                                                          */
  /*             set frame equal to PH_PACKET_HEADER_t.cal_frame_cnt - 1    */
  /*                                                                        */
  /*             SWITCH (PH_PACKET_HEADER_t.cal_type)                       */
  /*                CASE PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE    */
  /*                   set SC_PIXEL_QUALITY_DATA_t.SD_pix_qual[frame]       */
  /*                                                          [sequence]    */
  /*                      to qual_value                                     */
  /*                                                                        */
  /*                CASE PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE          */
  /*                   set SC_PIXEL_QUALITY_DATA_t.SRCA_pix_qual[frame]     */
  /*                                                            [sequence]  */
  /*                      to qual_value                                     */
  /*                                                                        */
  /*                CASE PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE         */
  /*                   set SC_PIXEL_QUALITY_DATA_t.BB_pix_qual[frame]       */
  /*                                                          [sequence]    */
  /*                      to qual_value                                     */
  /*                                                                        */
  /*                CASE PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE             */
  /*                   set SC_PIXEL_QUALITY_DATA_t.SV_pix_qual[frame]       */
  /*                                                          [sequence]    */
  /*                      to qual_value                                     */
  /*                                                                        */
  /*             END-SWITCH                                                 */
  /*          ENDIF                                                         */
  /*       ENDIF                                                            */
  /*                                                                        */
  /**************************************************************************/

        else if (pkt_header.cal_frame_cnt > 0) {
           frame = pkt_header.cal_frame_cnt - 1;
           switch (pkt_header.cal_type) {
              case PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE:
                 scan_pix->SD_pix_qual[frame][sequence] = qual_value; 
                 break;
 
              case PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE:      
                 scan_pix->SRCA_pix_qual[frame][sequence] = qual_value; 
                 break;

              case PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE:
                 scan_pix->BB_pix_qual[frame][sequence] = qual_value; 
                 break;

              case PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE:
                 scan_pix->SV_pix_qual[frame][sequence] = qual_value; 
                 break;
           }
        }
     }

  /**************************************************************************/
  /*    ELSE   ** Set frame count for nightmode packets **                  */
  /*       IF (PH_PACKET_HEADER_t.pkt_type equals                           */ 
  /*          PH_SEC_PKT_TYPE_NIGHT_GROUP)                                  */
  /*       THEN ** Night packets can only be collected at Earth sector **   */
  /*          set frame equal to PH_PACKET_HEADER_t.earth_frame_cnt - 1     */
  /*          IF frame is greater than 0                                    */
  /*          THEN                                                          */
  /*             set SC_PIXEL_QUALITY_DATA_t.EV_pix_qual[frame][0] to       */
  /*                qual_value                                              */
  /*             set SC_PIXEL_QUALITY_DATA_t.EV_pix_qual[frame][1] to       */
  /*                qual_value                                              */
  /*          ENDIF                                                         */
  /*       ENDIF                                                            */
  /*    ENDIF                                                               */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/
 
     else if (pkt_header.pkt_type == PH_SEC_PKT_TYPE_NIGHT_GROUP) {
        frame = pkt_header.earth_frame_cnt - 1;
        if (pkt_header.earth_frame_cnt > 0) {
           scan_pix->EV_pix_qual[frame][0] = qual_value; 
           scan_pix->EV_pix_qual[frame][1] = qual_value; 
        }

     }

   } /* End check for Eng pkt  */

 } /* End of routine update_pixel_qual_data */
