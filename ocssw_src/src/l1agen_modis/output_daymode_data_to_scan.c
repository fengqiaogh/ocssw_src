#include "L1A_prototype.h"
#include "hdfi.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"


void   output_daymode_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                    uint16              *pkt_contents,
                                    SC_SCAN_DATA_t      *L1A_scan)

/*
!C************************************************************************

!Description:  This routine determines whether the packet is an Earth View or
               Callibration type daymode packet and then calls the appropiate
               routines to output the packet's data to the scan's structure.

!Input Parameters:
               PH_PACKET_HEADER_t  *pkt_header   ** Buffer that contains pkt **
                                                 ** header info of the pckt  **

               uint16              *pkt_contents ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               SC_SCAN_DATA_t      *L1A_scan     ** The MODIS scan structure **
                                                 ** currently being built    **

Return Values: NONE

Externally Defined:
               uint16					        (hdfi.h)
               PH_PACKET_HEADER_t                               (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH                 (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_CAL                   (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE  (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE        (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE       (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE           (PH_pkt_hdr.h)
               SC_SCAN_DATA_t                                   (SC_scan.h)

Called By:
               put_pkt_cont_in_scan

Routines Called:
               put_cal_data_in_scan
               put_earth_data_in_scan

!Revision History:
               Revision 2.0  1997/09/09  11:33 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code for version 2.

               Revision 1.1  1997/08/27   9:35
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

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

               This routine was designed totally under the assumption that
               the packet header has been previously validated and there
               is no need to test for error conditions.

!END*************************************************************************
*/

 {
  /**************************************************************************/
  /* Figure out whether this packet is an Earth View packet or a            */
  /* Calibration packet, and then figure out which Cal Sector this packet   */
  /* is in.                                                                 */
  /*                                                                        */
  /**************************************************************************/
  /* IF PH_PACKET_HEADER_t.source_ID_type_flag equals                       */
  /*                                     PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH   */
  /* THEN                                                                   */
  /*    CALL put_earth_data_in_scan to output the packet's Earth View       */
  /*         Sector data to the scan structure                              */
  /*      INPUTS:  PH_PACKET_HEADER_t, pkt_contents, SC_SCAN_DATA_t.EV_250m,*/
  /*               SC_SCAN_DATA_t.EV_500m, SC_SCAN_DATA_t.EV_1km_day,       */
  /*               SC_SCAN_DATA_t.EV_1km_night                              */
  /*      OUTPUTS: SC_SCAN_DATA_t.EV_250m, SC_SCAN_DATA_t.EV_500m,          */
  /*               SC_SCAN_DATA_t.EV_1km_day, SC_SCAN_DATA_t.EV_1km_night   */
  /*      RETURN:  NONE                                                     */
  /*                                                                        */
  /**************************************************************************/

     if (pkt_header->source_ID_type_flag == PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH) {
        put_earth_data_in_scan (pkt_header, pkt_contents, 
                                L1A_scan->EV_250m,
                                L1A_scan->EV_500m, 
                                L1A_scan->EV_1km_day, 
                                L1A_scan->EV_1km_night);
     }

  /**************************************************************************/
  /* ELSE                                                                   */
  /*    SWITCH (PH_PACKET_HEADER_t.cal_type)                                */
  /*       CASE PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE             */
  /*          CALL put_cal_data_in_scan to output the packet's Solar        */
  /*               Diffuser Sector data to the scan structure               */
  /*            INPUTS:  PH_PACKET_HEADER_t, pkt_contents,                  */
  /*                     SC_SCAN_DATA_t.SD_250m, SC_SCAN_DATA_t.SD_500m     */
  /*                     SC_SCAN_DATA_t.SD_1km_day,                         */
  /*                     SC_SCAN_DATA_t.SD_1km_night                        */
  /*            OUTPUTS: SC_SCAN_DATA_t.SD_250m, SC_SCAN_DATA_t.SD_500m,    */
  /*                     SC_SCAN_DATA_t.SD_1km_day,                         */
  /*                     SC_SCAN_DATA_t.SD_1km_night                        */
  /*            RETURN:  NONE                                               */
  /*                                                                        */
  /**************************************************************************/

     else {
        switch (pkt_header->cal_type) {
           case PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE:
              put_cal_data_in_scan (pkt_header, pkt_contents,
                                    L1A_scan->SD_250m,
                                    L1A_scan->SD_500m,
                                    L1A_scan->SD_1km_day,
                                    L1A_scan->SD_1km_night);

              break;
                                   

  /**************************************************************************/
  /*       CASE PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE                   */
  /*          CALL put_cal_data_in_scan to output the packet's SRCA         */
  /*               Sector data to the scan structure                        */
  /*            INPUTS:  PH_PACKET_HEADER_t, pkt_contents,                  */
  /*                     SC_SCAN_DATA_t.SRCA_250m, SC_SCAN_DATA_t.SRCA_500m */
  /*                     SC_SCAN_DATA_t.SRCA_1km_day,                       */
  /*                     SC_SCAN_DATA_t.SRCA_1km_night                      */
  /*            OUTPUTS: SC_SCAN_DATA_t.SRCA_250m, SC_SCAN_DATA_t.SRCA_500m,*/
  /*                     SC_SCAN_DATA_t.SRCA_1km_day,                       */
  /*                     SC_SCAN_DATA_t.SRCA_1km_night                      */
  /*            RETURN:  NONE                                               */
  /*                                                                        */
  /**************************************************************************/

           case PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE:
              put_cal_data_in_scan (pkt_header, pkt_contents,
                                    L1A_scan->SRCA_250m,
                                    L1A_scan->SRCA_500m,
                                    L1A_scan->SRCA_1km_day,
                                    L1A_scan->SRCA_1km_night);

              break;



  /**************************************************************************/
  /*       CASE PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE                  */
  /*          CALL put_cal_data_in_scan to output the packet's BlackBody    */
  /*               Sector data to the scan structure                        */
  /*            INPUTS:  PH_PACKET_HEADER_t, pkt_contents,                  */
  /*                     SC_SCAN_DATA_t.BB_250m, SC_SCAN_DATA_t.BB_500m,    */
  /*                     SC_SCAN_DATA_t.BB_1km_day,                         */
  /*                     SC_SCAN_DATA_t.BB_1km_night                        */
  /*            OUTPUTS: SC_SCAN_DATA_t.BB_250m, SC_SCAN_DATA_t.BB_500m,    */
  /*                     SC_SCAN_DATA_t.BB_1km_day,                         */
  /*                     SC_SCAN_DATA_t.BB_1km_night                        */
  /*            RETURN:  NONE                                               */
  /*                                                                        */
  /**************************************************************************/

           case PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE:      
              put_cal_data_in_scan (pkt_header, pkt_contents,
                                    L1A_scan->BB_250m,
                                    L1A_scan->BB_500m,
                                    L1A_scan->BB_1km_day,
                                    L1A_scan->BB_1km_night);

              break;



  /**************************************************************************/
  /*       CASE PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE                      */
  /*          CALL put_cal_data_in_scan to output the packet's Space        */
  /*               View Sector data to the scan structure                   */
  /*            INPUTS:  PH_PACKET_HEADER_t, pkt_contents,                  */
  /*                     SC_SCAN_DATA_t.SV_250m, SC_SCAN_DATA_t.SV_500m,    */
  /*                     SC_SCAN_DATA_t.SV_1km_day,                         */
  /*                     SC_SCAN_DATA_t.SV_1km_night                        */
  /*            OUTPUTS: SC_SCAN_DATA_t.SV_250m, SC_SCAN_DATA_t.SV_500m,    */
  /*                     SC_SCAN_DATA_t.SV_1km_day,                         */
  /*                     SC_SCAN_DATA_t.SV_1km_night                        */
  /*            RETURN:  NONE                                               */
  /*    END_SWITCH                                                          */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

           case PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE:
              put_cal_data_in_scan (pkt_header, pkt_contents,
                                    L1A_scan->SV_250m,
                                    L1A_scan->SV_500m,
                                    L1A_scan->SV_1km_day,
                                    L1A_scan->SV_1km_night);
                                        
              break;

        }   /* end_switch */
     }      /* endif      */

 }  /* End of routine output_daymode_data_to_scan */
