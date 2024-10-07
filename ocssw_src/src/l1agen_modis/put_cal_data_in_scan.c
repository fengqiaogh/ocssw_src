#include "hdfi.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "SC_scan.h"
#include "L1A_prototype.h"


void   put_cal_data_in_scan (PH_PACKET_HEADER_t  *pkt_header,
                             uint16              *pkt_cont,
                             SC_CAL_250M          SC_250m,
                             SC_CAL_500M          SC_500m,
                             SC_CAL_1KM_DAY       SC_1km_day, 
                             SC_CAL_1KM_NIGHT     SC_1km_night)

/*
!C************************************************************************

!Description:  The routine takes an unpacked daymode packet containing
               calibration data at one of the calibration targets, extracts
               it, and places the data in the proper arrays within the scan
               structure. 

!Input Parameters:
               PH_PACKET_HEADER_t  *pkt_header   ** The structure containing **
                                                 ** the current packet's     **
                                                 ** packet header            **

               int16               *pkt_cont     ** The structure containing **
                                                 ** the current packet's     **
                                                 ** unpacked contents        **

!Output Parameters:
               None

!Input/Output Parameters:
               SC_CAL_250M         SC_250m       ** array that holds earth/  **
                                                 ** calibration sector 250m  **
                                                 ** radiances                **

               SC_CAL_500M         SC_500m       ** array that holds earth/  **
                                                 ** calibration sector 500m  **
                                                 ** radiances                **

               SC_CAL_1KM_DAY      SC_1km_day    ** array that holds earth/  **
                                                 ** calibration sector 1km   **
                                                 ** day radiances            **

               SC_CAL_1KM_NIGHT    SC_1km_night  ** array that holds earth/  **
                                                 ** calibration sector 1km   **
                                                 ** night radiances          **

Return Values: 
               None

Externally Defined:
               int16                                   (hdfi.h)
               SC_CAL_250M                             (SC_scan.h)
               SC_CAL_500M                             (SC_scan.h)
               SC_CAL_1KM_DAY                          (SC_scan.h)
               SC_CAL_1KM_NIGHT                        (SC_scan.h)
               PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP      (PH_pkt_hdr.h)
               PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP     (PH_pkt_hdr.h)
               PD_DN_FIRST_250M_BAND                   (PD_pkt_data.h)
               PD_DN_LAST_250M_BAND                    (PD_pkt_data.h)
               PD_DN_BAND_RATIO_250M                   (PD_pkt_data.h)
               PD_DN_FIRST_500M_BAND                   (PD_pkt_data.h)
               PD_DN_LAST_500M_BAND                    (PD_pkt_data.h)
               PD_DN_BAND_RATIO_500M                   (PD_pkt_data.h)
               PD_DN_FIRST_1KM_DAY_BAND                (PD_pkt_data.h)
               PD_DN_LAST_1KM_DAY_BAND                 (PD_pkt_data.h)
               PD_DN_FIRST_1KM_NIGHT_BAND              (PD_pkt_data.h)
               PD_DN_LAST_1KM_NIGHT_BAND               (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_DETECTORS              (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_DETECTORS            (PD_pkt_data.h)
               PD_DN_NUM_IFOVS_IN_DAY_PKT              (PD_pkt_data.h)
               PD_DN_FIRST_IFOV_DAY_PKT_1              (PD_pkt_data.h)
               PD_DN_LAST_IFOV_DAY_PKT_1               (PD_pkt_data.h)
               PD_DN_FIRST_IFOV_DAY_PKT_2              (PD_pkt_data.h)
               PD_DN_LAST_IFOV_DAY_PKT_2               (PD_pkt_data.h)

Called By:
               output_daymode_data_to_scan

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/06/18  16:40 EDT 
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.  (original design based on version1 
                                  routine: output_cal_data_to_scan)

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The CODE below was developed in C language.

               This routine was designed totally under the asumption that 
               the packet header has been priviously validated and there 
               is no need to test for error conditions.

               The structure of "SC_*" is determined by the Level 1A 
               project specification.

!END*************************************************************************
*/

 {
  /*************************************************************************/
  /*                                                                       */
  /*          Declare the local variables and initialize them.             */
  /*                                                                       */
  /*************************************************************************/
  /*                                                                       */
  /* set start_ifov to 0                                                   */
  /* set end_ifov equal to PD_DN_NUM_IFOVS_IN_DAY_PKT                      */
  /* set pkt_cont_pos equal to PH_NUM_12BIT_WORDS_IN_HEADER                */
  /*                                                                       */
  /*************************************************************************/

  int        start_ifov;           /* starting ifov position within scan     */
  int        end_ifov;             /* ending ifov position within scan       */
  int        ifov;                 /* loop variable holding current ifov     */
  int        band;                 /* Band indicator (used as loop variable) */
  int        band_num;             /* Band indicator corrected for each      */
                                   /* detector group's start position        */
  int        sample;               /* Sample indicator for 250m and 500m     */
                                   /* data (used as loop variable)           */
  int        det;                  /* detector multiplier for 250m and 500m  */
                                   /* data (used as loop variable)           */
  int        line_num;             /* Calculated line number; used to        */
                                   /* properly place data in science arrays  */
  int        frame;                /* Calculated frame number; used to       */
                                   /* properly place data in science arrays  */
  int        frame_ifov_offset;    /* Offset used to adjust a packet's       */
                                   /* IFOVs to the proper scan location      */
  int        scan_ifov;            /* IFOV within the scan where the         */
                                   /* current packet's contents belong       */
  int        pkt_cont_pos;         /* index into packet_contents             */
  int8       frame_cnt;            /* Frame count for Calibration packets    */
  

  band              = 0;
  band_num          = 0;
  sample            = 0;
  det               = 0;
  line_num          = 0;
  frame             = 0;
  frame_cnt         = 0;
  frame_ifov_offset = 0;
  scan_ifov         = 0;
  ifov              = 0;
  start_ifov        = 0;
  end_ifov          = PD_DN_NUM_IFOVS_IN_DAY_PKT;
  pkt_cont_pos      = PH_NUM_12BIT_WORDS_IN_HEADER;



  /**************************************************************************/
  /*                                                                        */
  /* Get packet frame count. Then determine whether this packet contains    */
  /* IFOVs 1 through 5 or 6 through 10.                                     */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set frame_cnt equal to PH_PACKET_HEADER_t.cal_frame_cnt                */
  /*                                                                        */
  /* Set frame_ifov_offset equal to 0                                       */
  /* IF PH_PACKET_HEADER_t.sequence_flag equals                             */
  /*                                 PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP    */
  /* THEN                                                                   */
  /*    Set frame_ifov_offset equal to 5                                    */
  /* ENDIF                                                                  */
  /*                                                                        */
  /**************************************************************************/

     frame_cnt = pkt_header->cal_frame_cnt;

     frame_ifov_offset = 0;
     if (pkt_header->sequence_flag == PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP)
        frame_ifov_offset = 5;


  /**************************************************************************/
  /*                                                                        */
  /*          Extract the packet contents into the scan structure.          */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* FOR ifov equal to start_ifov to end_ifov                               */
  /*    set scan_ifov equal to (ifov + frame_ifov_offset)                   */
  /*                                                                        */
  /**************************************************************************/

     for (ifov = start_ifov; ifov < end_ifov; ifov++) {
        scan_ifov = ifov + frame_ifov_offset;


  /**************************************************************************/
  /*                                                                        */
  /* Within each ifov set, the 250 meter data appear first, organized by    */
  /* band_num, each band's radiances organized by the 4 samples taken by    */
  /* each detector, then by the 4 detectors within the ifov. And the line   */
  /* number is corrected to account for the instrument detector numbering.  */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /*    FOR band equal to PD_DN_FIRST_250M_BAND to PD_DN_LAST_250M_BAND     */
  /*       set band_num equal to (band - PD_DN_FIRST_250M_BAND)             */
  /*       FOR sample equal to 0 upto PD_DN_BAND_RATIO_250M                 */
  /*          set frame equal to ((frame_cnt - 1) * PD_DN_BAND_RATIO_250M) +*/
  /*                             sample                                     */
  /*          FOR det equals to 0 upto PD_DN_NUM_250M_DETECTORS_IN_IFOV     */
  /*             set line_num equal to ((PD_DN_NUM_250M_DETECTORS - 1) -    */
  /*                 ((scan_ifov * PD_DN_NUM_250M_DETECTORS_IN_IFOV) + det))*/
  /*             set SC_250m(line_num,band_num,frame) equal to              */
  /*                                                 pkt_cont(pkt_cont_pos) */
  /*             set pkt_cont_pos to  pkt_cont_pos + 1                      */
  /*          ENDFOR (next det)                                             */
  /*       ENDFOR (next sample)                                             */
  /*    ENDFOR (next band)                                                  */
  /*                                                                        */
  /**************************************************************************/

        for (band = PD_DN_FIRST_250M_BAND; band <= PD_DN_LAST_250M_BAND; band++){
           band_num = band - PD_DN_FIRST_250M_BAND;
           for (sample = 0; sample < PD_DN_BAND_RATIO_250M; sample++){
              frame = ((frame_cnt - 1) * PD_DN_BAND_RATIO_250M) + sample;
              for (det = 0; det < PD_DN_NUM_250M_DETECTORS_IN_IFOV; det++){
                 line_num = ((PD_DN_NUM_250M_DETECTORS - 1) -
                      ((scan_ifov * PD_DN_NUM_250M_DETECTORS_IN_IFOV) + det));
                 SC_250m[line_num][band_num][frame] = pkt_cont[pkt_cont_pos];
                 pkt_cont_pos++;
              }
           }
        }

 
  /**************************************************************************/
  /*                                                                        */
  /* The 500 meter data is next, with each band's radiances organized be    */
  /* band_num, then by the 2 samples taken by each detector, then by the 2  */
  /* detectors within the ifov. And the line number is corrected to         */
  /* account for the instrument detector numbering.                         */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /*    FOR band equal to PD_DN_FIRST_500M_BAND to PD_DN_LAST_500M_BAND     */
  /*       set band_num equal to (band - PD_DN_FIRST_500M_BAND)             */
  /*       FOR sample equal to 0 upto PD_DN_BAND_RATIO_500M                 */
  /*          set frame equal to ((frame_cnt - 1) * PD_DN_BAND_RATIO_500M) +*/
  /*                             sample                                     */
  /*          FOR det equals to 0 upto PD_DN_NUM_500M_DETECTORS_IN_IFOV     */
  /*             set line_num equal to ((PD_DN_NUM_500M_DETECTORS - 1) -    */
  /*             ((scan_ifov * PD_DN_NUM_500M_DETECTORS_IN_IFOV) + det))    */
  /*             set SC_500m(line_num,band_num,frame) equal to              */
  /*                                             pkt_cont(pkt_cont_pos)     */
  /*             set pkt_cont_pos to  pkt_cont_pos + 1                      */
  /*          ENDFOR (next det)                                             */
  /*       ENDFOR (next sample)                                             */
  /*    ENDFOR (next band)                                                  */
  /*                                                                        */
  /**************************************************************************/

        for (band = PD_DN_FIRST_500M_BAND; band <= PD_DN_LAST_500M_BAND; band++){
           band_num = band - PD_DN_FIRST_500M_BAND;
           for (sample = 0; sample < PD_DN_BAND_RATIO_500M; sample++){
              frame = ((frame_cnt - 1) * PD_DN_BAND_RATIO_500M) + sample;
              for (det = 0; det < PD_DN_NUM_500M_DETECTORS_IN_IFOV; det++){
                 line_num = ((PD_DN_NUM_500M_DETECTORS - 1) -
                       ((scan_ifov * PD_DN_NUM_500M_DETECTORS_IN_IFOV) + det));
                 SC_500m[line_num][band_num][frame] = pkt_cont[pkt_cont_pos];
                 pkt_cont_pos++;
              }
           }
        }


  /**************************************************************************/
  /*                                                                        */
  /* The 1km day data is next in each ifov set, with one sample per band    */
  /* And the line number is corrected to account for the instrument         */
  /* detector numbering.                                                    */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /*    set frame equal to frame_cnt - 1                                    */
  /*    FOR band = PD_DN_FIRST_1KM_DAY_BAND to PD_DN_LAST_1KM_DAY_BAND      */
  /*       set band_num equal to (band - PD_DN_FIRST_1KM_DAY_BAND)          */
  /*       set line_num equal to ((PD_DN_NUM_1KMDAY_DETECTORS - 1) -        */
  /*                             scan_ifov)                                 */
  /*       set SC_1km_day(line_num,band_num,frame) equal to                 */
  /*                                             pkt_cont(pkt_cont_pos)     */
  /*       set pkt_cont_pos to pkt_cont_pos + 1                             */
  /*    ENDFOR (next band)                                                  */
  /*                                                                        */
  /**************************************************************************/

        frame = frame_cnt - 1;
        for (band = PD_DN_FIRST_1KM_DAY_BAND;
             band <= PD_DN_LAST_1KM_DAY_BAND;
             band++) {
           band_num = band - PD_DN_FIRST_1KM_DAY_BAND;
           line_num = ((PD_DN_NUM_1KMDAY_DETECTORS - 1) - scan_ifov);
           SC_1km_day[line_num][band_num][frame] = pkt_cont[pkt_cont_pos];
           pkt_cont_pos++;
        }


  /**************************************************************************/
  /*                                                                        */
  /* The 1km night data is last in each ifov set, with one sample per band  */
  /* And the line number is corrected to account for the instrument         */
  /* detector numbering.                                                    */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /*    set frame equal to frame_cnt - 1                                    */
  /*    FOR band = PD_DN_FIRST_1KM_NIGHT_BAND to PD_DN_LAST_1KM_NIGHT_BAND  */
  /*       set band_num equal to (band - PD_DN_FIRST_1KM_NIGHT_BAND)        */
  /*       set line_num equal to ((PD_DN_NUM_1KMNIGHT_DETECTORS - 1 ) -     */
  /*                             scan_ifov)                                 */
  /*       set SC_1km_night(line_num,band_num,frame) equal to               */
  /*                                             pkt_cont(pkt_cont_pos)     */
  /*       set pkt_cont_pos to pkt_cont_pos + 1                             */
  /*    ENDFOR (next band)                                                  */
  /* ENDFOR (next ifov)                                                     */
  /*                                                                        */
  /**************************************************************************/

        frame = frame_cnt - 1;
        for (band = PD_DN_FIRST_1KM_NIGHT_BAND;
             band <= PD_DN_LAST_1KM_NIGHT_BAND;
             band++) {
           band_num = band - PD_DN_FIRST_1KM_NIGHT_BAND;
           line_num = ((PD_DN_NUM_1KMNIGHT_DETECTORS - 1 ) - scan_ifov);
           SC_1km_night[line_num][band_num][frame] = pkt_cont[pkt_cont_pos];
           pkt_cont_pos++;
        }
     }     /* next ifov */

 } /* End of routine put_cal_data_in_scan */
