#include "L1A_prototype.h"
#include "hdfi.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"


void   output_nightmode_data_to_scan (PH_PACKET_HEADER_t  *pkt_header,
                                      uint16              *pkt_contents,
                                      SC_EV_1KM_NIGHT      EV_1km_night )

/*
!C************************************************************************

!Description:  The routine takes an unpacked nightmode packet and places the
               packet's contents in the scan structure. 

!Input Parameters:
               PH_PACKET_HEADER_t  *pkt_header   ** The structure containing **
                                                 ** the current packet's     **
                                                 ** packet header            **

               uint16              *pkt_contents ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               None
                         
!Input/Output Parameters:                         
               SC_EV_1KM_NIGHT     EV_1km_night  ** array that holds earth   **
                                                 ** view sector 1km night    **
                                                 ** radiances                **

Return Values: 
               None


Externally Defined:      
               int16                                  (hdfi.h)
               uint16                                 (hdfi.h)
               PH_PACKET_HEADER_t                     (PH_pkt_hdr.h)
               PD_DN_FIRST_1KM_NIGHT_BAND             (PD_pkt_data.h)
               PD_DN_LAST_1KM_NIGHT_BAND              (PD_pkt_data.h)
               PD_DN_NUM_IFOVS_IN_NIGHT_PKT           (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_DETECTORS           (PD_pkt_data.h)
               SC_EV_1KM_NIGHT                        (SC_scan.h)

Called By:
               put_pkt_cont_in_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/09  11:50 EDT 
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Recreating module per version2 developement.

               Revision 1.1  1997/08/27   9:40
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
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/
  /*                                                                       */
  /* set pkt_cont_pos equal to PH_NUM_12BIT_WORDS_IN_HEADER                */
  /* set start_ifov equal to 0                                             */
  /* set end_ifov equal to PD_NUM_IFOVS_IN_NIGHT_PKT - 1                   */
  /*                                                                       */
  /*************************************************************************/

  int              start_ifov;     /* starting ifov position within scan     */
  int              end_ifov;       /* ending ifov position within scan       */
  int              det_num;        /* Calculated detector number; used to    */
                                   /* properly place data in science arrays  */
  int              band;           /* Band indicator (used as loop variable) */
  int              band_num;       /* Band indicator corrected for each      */
                                   /* detector group's start position        */
  int              line_num;       /* Calculated line number; used to pro-   */
                                   /* perly extract data in science arrays   */
  int              frame;          /* Calculated frame number; used to pro-  */
                                   /* perly extract data in science arrays   */
  int              pkt_cont_pos;   /* index into packet_contents             */


  det_num      = 0;
  band         = 0;
  band_num     = 0;
  line_num     = 0;
  frame        = 0;
  start_ifov   = 0;
  end_ifov     = PD_DN_NUM_IFOVS_IN_NIGHT_PKT;
  pkt_cont_pos = PH_NUM_12BIT_WORDS_IN_HEADER;


  /**************************************************************************/
  /*                                                                        */
  /* First check the Validity of the packet's frame count, then extract     */
  /* the nightmode data                                                     */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* IF PH_PACKET_HEADER_t.earth_frame_cnt is greater than zero             */
  /* THEN                                                                   */
  /*    set frame equal to PH_PACKET_HEADER_t.earth_frame_cnt               */
  /*    FOR det_num equal to start_ifov to end_ifov                         */
  /*     {The 1km night data is all that is in the packet, with one    }    */
  /*     {sample per band                                              }    */
  /*       FOR band equal to PD_DN_FIRST_1KM_NIGHT_BAND to                  */
  /*                                        PD_DN_LAST_1KM_NIGHT_BAND       */
  /*          set band_num equal to band - PD_DN_FIRST_1KM_NIGHT_BAND       */
  /*     {Correct line number to account for instrument detector       }    */
  /*     {numbering                                                    }    */
  /*          set line_num equal to (PD_DN_NUM_1KMNIGHT_DETECTORS - 1) -    */
  /*                                                              det_num   */
  /*          set EV_1km_night(line_num,band_num,frame) equal to            */
  /*                                          pkt_contents(pkt_cont_pos)    */
  /*          set pkt_cont_pos to pkt_cont_pos + 1                          */
  /*       ENDFOR (next band)                                               */
  /*    ENDFOR (next ifov)                                                  */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

     if (pkt_header->earth_frame_cnt > 0) {
        frame = pkt_header->earth_frame_cnt - 1;
        for (det_num = start_ifov; det_num < end_ifov; det_num++) {
           for (band  = PD_DN_FIRST_1KM_NIGHT_BAND; 
                band <= PD_DN_LAST_1KM_NIGHT_BAND;
                band++) {
              band_num = band - PD_DN_FIRST_1KM_NIGHT_BAND;
              line_num = (PD_DN_NUM_1KMNIGHT_DETECTORS - 1) - det_num;
              EV_1km_night[line_num][band_num][frame] =
                                                    pkt_contents[pkt_cont_pos];
              pkt_cont_pos++;
           } 
        }
     }

 } /* End of routine output_nightmode_data_to_scan */
