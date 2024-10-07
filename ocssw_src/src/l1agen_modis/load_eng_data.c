#include "PGS_MODIS_35005.h"
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_TYPES.h"
#include "PH_pkt_hdr.h"
#include "EN_eng_data.h"
#include "L1A_prototype.h"


PGSt_SMF_status  load_eng_data (PGSt_double         scan_rate,
                                EN_VDATA_TYPE_t     *eng_data,
                                PGSt_IO_L0_Packet   *pkt,
                                PH_PACKET_HEADER_t  *pkt_header,
                                PGSt_IO_L0_VirtualDataSet  L0_file)

/*
!C************************************************************************

!Description:  This routine will populate the engineering data lists with 
               data from any engineering packets that are encountered.

!Input Parameters:
               PGSt_IO_L0_VirtualDataSet  L0_file
                                         ** Virtual data set - L0 file    **

               PGSt_double         scan_rate
                                         ** The length of each scan in    **
                                         ** seconds                       **

!Output Parameters:
               PGSt_IO_L0_Packet   *pkt  ** Packet buffer which contains  **
                                         ** packed packet                 **

               PH_PACKET_HEADER_t  *pkt_header
                                         ** Buffer that contains the pkt  **
                                         ** header info of the packet     **

!Input/Output Parameters:
               EN_VDATA_TYPE_t     *eng_data
                                         ** The structure containing the  **
                                         ** engineering data              **

Return Values: 
               MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
               MODIS_W_PRELOAD_ENG_DATA        (PGS_MODIS_35005.h)
               MODIS_F_PKT_READ_FAILED         (PGS_MODIS_35005.h)
               MODIS_W_NO_MORE_PACKETS         (PGS_MODIS_35005.h)
               MODIS_E_CORRUPT_PKT_HDR         (PGS_MODIS_35005.h)
               MODIS_E_CHECKSUM_NOT_VALID      (PGS_MODIS_35005.h)

Externally Defined:  
               PGSt_SMF_status                 (PGS_SMF.h) 
               PGS_TRUE                        (PGS_SMF.h)
               PGS_FALSE                       (PGS_SMF.h)
               PGSt_double                     (PGS_TYPES.h)
               PGSt_IO_L0_Packet               (PGS_IO.h)
               PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
               PH_SEC_PKT_TYPE_ENG1_GROUP      (PH_pkt_hdr.h) 
               PH_SEC_PKT_TYPE_ENG2_GROUP      (PH_pkt_hdr.h) 
               PH_PACKET_HEADER_t              (PH_pkt_hdr.h)
               EN_VDATA_TYPE_t                 (EN_eng_data.h)

Called By:     
               initialize_level1a

Routines Called:
               read_a_packet
               log_fmt_msg                
               unpack_packet_header       
               unpack_packet_contents     
               check_checksum  
               process_eng_packet           
               compute_SD_start_time      

!Revision History:
               Revision 2.1  1997/08/27  22:20 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 2.0  1997/07/15  11:20 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Created Module for version 2 based on load_decom_data.pdl
               of version 1.

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
               The CODE was written in C language.

!END**********************************************************************
*/

 {
  /**************************************************************************/
  /*                      Define Global Variables                           */
  /*                                                                        */
  /**************************************************************************/

  extern PGSt_double global_first_gran_start_time;
                                          /* Start time for the first L1A   */
                                          /* output granule to be created   */


  /*************************************************************************/
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/
  /* set returnStatus equal to MODIS_S_SUCCESS                             */
  /* set first_time_in_loop to PGS_TRUE                                    */
  /* set next_scan_start_time to -1.0                                      */
  /* set continue_processing equal to PGS_TRUE                             */
  /* set scan_num to 0                                                     */
  /*                                                                       */
  /*************************************************************************/

  char   *routine = "load_eng_data";
  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  PGSt_SMF_status  tempStatus;     /* SMF-style message returned by function */
  PGSt_SMF_boolean first_time_in_loop;
                                   /* loop control variable                  */
  PGSt_SMF_boolean continue_processing;
                                   /* loop control variable                  */
  PGSt_double next_scan_start_time;/* variable that holds the next scan      */
                                   /* start time                             */
  PGSt_double current_scan_start_time;
                                   /* variable that holds the current scan   */
                                   /* start time                             */
  uint16  pkt_contents[PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT +
                      PH_NUM_12BIT_WORDS_IN_HEADER + 1];
                                   /* The current packet, with the 12-bit    */
                                   /* words unpacked into 16-bit integers    */
  int               scan_num;      /* Scan counter                           */


  tempStatus              = MODIS_S_SUCCESS;
  returnStatus            = MODIS_S_SUCCESS;
  first_time_in_loop      = PGS_TRUE;
  continue_processing     = PGS_TRUE;
  next_scan_start_time    = -1.0;
  current_scan_start_time = 0.0;
  scan_num                = 0;


  /*************************************************************************/
  /* DO_WHILE (continue_processing equals PGS_TRUE)                        */
  /*    CALL read_a_packet to read a packet of data from the L0 virtual    */
  /*         data set                                                      */
  /*       INPUTS:  L0_file                                                */
  /*       OUTPUTS: L0_file, pkt                                           */
  /*       RETURN:  returnStatus                                           */
  /*                                                                       */
  /*    IF returnStatus is equal to MODIS_S_SUCCESS                        */
  /*    THEN                                                               */
  /*       CALL unpack_packet_header to unpack the packet's header         */
  /*            information to determine what kind of packet it is         */
  /*          INPUTS:  pkt, PH_PACKET_HEADER_t                             */
  /*          OUTPUTS: NONE                                                */
  /*          RETURN:  tempStatus                                          */
  /*                                                                       */
  /*       IF tempStatus is equal to MODIS_S_SUCCESS                       */
  /*       THEN                                                            */
  /*          IF (PH_PACKET_HEADER_t.pkt_TAI_time is greater than or equal */
  /*              to first_gran_start_time)                                */
  /*          THEN                                                         */
  /*                                                                       */
  /*************************************************************************/

  while (continue_processing == PGS_TRUE) {
     returnStatus = read_a_packet (&L0_file, pkt);
    
     if (returnStatus == MODIS_S_SUCCESS) {
        tempStatus = unpack_packet_header (pkt, pkt_header);

        if (tempStatus == MODIS_S_SUCCESS) {

           if (pkt_header->pkt_TAI_time >= global_first_gran_start_time) {



  /*************************************************************************/
  /*             IF (next_scan_start_time < 0) -> { this is to ensure that */
  /*                                                this block of code     */
  /*                                                is executed only once }*/
  /*             THEN                                                      */
  /*                CALL compute_SD_start_time to calculate the start time */
  /*                     for the current scan of the current packet        */
  /*                   INPUTS:  PH_PACKET_HEADER_t                         */
  /*                   OUTPUTS: current_scan_start_time                    */
  /*                   RETURN:  NONE                                       */
  /*                                                                       */
  /*                set next_scan_start_time to current_scan_start_time +  */
  /*                                            scan_rate                  */
  /*             ENDIF                                                     */
  /*                                                                       */
  /*************************************************************************/

              if (next_scan_start_time < 0.0) {
                 compute_SD_start_time (pkt_header, &current_scan_start_time);
                 next_scan_start_time = current_scan_start_time + scan_rate;
              }



  /*************************************************************************/
  /*             IF (current_scan_start_time is greater than or equal to   */
  /*                 the first_gran_start_time)                            */
  /*             THEN                                                      */
  /*                set continue_processing equal to PGS_FALSE             */
  /*                IF first_time_in_loop is equal to PGS_TRUE             */
  /*                                       -> { this is the first time in  */
  /*                                            while loop }               */
  /*                THEN                                                   */
  /*                   set returnStatus to MODIS_W_PRELOAD_ENG_DATA        */
  /*                   CALL log_fmt_msg to report warning in preloading    */
  /*                        eng data                                       */
  /*                      INPUTS: returnStatus, routine, msg               */
  /*                      OUTPUT: None                                     */
  /*                      RETURN: None                                     */
  /*                ENDIF                                                  */
  /*             ELSE                                                      */
  /*                IF (packet's time is greater than or equal to the      */
  /*                    next_scan_start_time)                              */
  /*                THEN                                                   */
  /*                   set continue_processing equal to PGS_FALSE          */
  /*                ENDIF                                                  */
  /*             ENDIF                                                     */
  /*          ENDIF                                                        */
  /*                                                                       */
  /*************************************************************************/

              if (current_scan_start_time >= global_first_gran_start_time) {
                 continue_processing = PGS_FALSE;

                 if (first_time_in_loop) {
                    returnStatus = MODIS_W_PRELOAD_ENG_DATA;
                    log_fmt_msg (MODIS_W_PRELOAD_ENG_DATA, routine, "");
                 }
              }
              else 
                 if (pkt_header->pkt_TAI_time >= next_scan_start_time) 
                    continue_processing = PGS_FALSE;
           }



  /*************************************************************************/
  /*          IF ( (continue_processing equals PGS_TRUE) AND               */
  /*               (PH_PACKET_HEADER_t.pkt_type is equal to                */
  /*                PH_SEC_PKT_TYPE_ENG1_GROUP OR                          */
  /*                PH_SEC_PKT_TYPE_ENG2_GROUP) )                          */
  /*          THEN                                                         */
  /*             CALL unpack_packet_contents to unpack the packet contents */
  /*                  of the virtual data set                              */
  /*                INPUTS:  pkt                                           */
  /*                OUTPUTS: pkt_contents                                  */
  /*                RETURN:  None                                          */
  /*                                                                       */
  /*             CALL check_checksum to check the packet checksum          */
  /*                INPUTS:  pkt_contents, PH_PACKET_HEADER_t              */
  /*                OUTPUTS: NONE                                          */
  /*                RETURN:  tempStatus                                    */
  /*                                                                       */
  /*************************************************************************/

           if ( (continue_processing == PGS_TRUE) &&
                ((pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) ||
                 (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP)) ) {
              unpack_packet_contents (pkt, pkt_header, pkt_contents);
              tempStatus = check_checksum (*pkt_header, pkt_contents);



  /*************************************************************************/
  /*             IF tempStatus is equal to MODIS_S_SUCCESS                 */
  /*             THEN                                                      */
  /*                CALL process_eng_packet to process the Eng data from   */
  /*                     the Eng packets and store it in the Eng data      */
  /*                     structure                                         */
  /*                   INPUTS:  PH_PACKET_HEADER_t, scan_number, pkt,      */
  /*                            eng_data                                   */
  /*                   OUTPUTS: eng_data                                   */
  /*                   RETURN:  NONE                                       */
  /*             ELSE                                                      */
  /*                set tempStatus to MODIS_E_CHECKSUM_NOT_VALID           */
  /*                set msg equal to "descrepancy in packet checksum"      */
  /*                CALL log_fmt_msg to report error in checksum           */
  /*                     descrepancy                                       */
  /*                   INPUTS: Status, routine, msg                        */
  /*                   OUTPUT: None                                        */
  /*                   RETURN: None                                        */
  /*             ENDIF                                                     */
  /*          ENDIF                                                        */
  /*                                                                       */
  /*************************************************************************/

              if (tempStatus == MODIS_S_SUCCESS)
                 process_eng_packet (eng_data, scan_num, pkt_header, pkt);
              else {
                 tempStatus = MODIS_E_CHECKSUM_NOT_VALID;
                 log_fmt_msg (MODIS_E_CHECKSUM_NOT_VALID, routine, "");
              }
           }
        }



  /*************************************************************************/
  /*       ELSE                                                            */
  /*          set tempStatus to MODIS_E_CORRUPT_PKT_HDR                    */
  /*          CALL log_fmt_msg to report error in unpacking packet header  */
  /*             INPUTS: Status, routine, msg                              */
  /*             OUTPUT: None                                              */
  /*             RETURN: None                                              */
  /*       ENDIF                                                           */
  /*                                                                       */
  /*************************************************************************/

        else {
           tempStatus = MODIS_E_CORRUPT_PKT_HDR;
           log_fmt_msg (MODIS_E_CORRUPT_PKT_HDR, routine, "");
        }
     }

  /*************************************************************************/
  /*    ELSE                                                               */
  /*       set continue_processing equal to PGS_FALSE                      */
  /*       set msg equal to "error received in getting packet from virtual */
  /*           data set"                                                   */
  /*       CALL log_fmt_msg to report error in getting packet from virtual */
  /*            data set"                                                  */
  /*          INPUTS: returnStatus, routine, msg                           */
  /*          OUTPUT: None                                                 */
  /*          RETURN: None                                                 */
  /*    ENDIF                                                              */
  /*                                                                       */
  /*    set first_time_in_loop to PGS_FALSE                                */
  /*                                                                       */
  /* END_WHILE                                                             */
  /*                                                                       */
  /* return returnStatus                                                   */
  /*                                                                       */
  /*************************************************************************/

     else {
        continue_processing = PGS_FALSE;
        log_fmt_msg (returnStatus, routine,
                     "error received in getting packet from virtual data set");
     }

     first_time_in_loop = PGS_FALSE;
  }
    
  return returnStatus;

 } /* End of routine load_eng_data */
