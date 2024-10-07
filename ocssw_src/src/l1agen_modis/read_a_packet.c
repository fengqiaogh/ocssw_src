#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "PH_pkt_hdr.h"
#include "PC_pcf_info.h"
#include "PD_pkt_data.h"
#include "L1A_prototype.h"

PGSt_SMF_status  read_a_packet (PGSt_IO_L0_VirtualDataSet  *L0_file,
                                PGSt_IO_L0_Packet          *pkt)

/*
!C************************************************************************

!Description:  This function reads a packet of data from the L0 file.

!Input Parameters:
               None

!Output Parameters:
               PGSt_IO_L0_Packet          *pkt      ** L0 data packet       **

!Input/Output Parameters:
               PGSt_IO_L0_VirtualDataSet  *L0_file  ** L0 virtual data set  **

Return Values: 
               MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
               MODIS_F_PKT_READ_FAILED         (PGS_MODIS_35005.h)
               MODIS_W_NO_MORE_PACKETS         (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status                 (PGS_SMF.h)
               PGSt_IO_L0_Packet               (PGS_IO.h)
               PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
               PGSt_Logical                    (PGS_TYPES)
               PGS_TRUE                        (PGS_SMF.h)
               PGS_FALSE                       (PGS_SMF.h)
               PC_PRIOR_L0_PCF_ID              (PC_pcf_info.h)
               PC_CURRENT_L0_PCF_ID            (PC_pcf_info.h)
               PD_PKT_BUF_MAX                  (PD_pkt_data.h)
               global_L0_logical               (level1a)
               PGSIO_W_L0_END_OF_VIRTUAL_DS    (PGS_IO.h)
               PGSIO_M_L0_HEADER_CHANGED       (PGS_IO.h)
               PGSIO_W_L0_PKT_BUF_TRUNCATE     (PGS_IO.h)

Called By:
               process_a_packet
               load_eng_data

Routines Called:
               PGS_SMF_TestSuccessLevel
               PGS_IO_L0_GetPacket 
               PGS_IO_L0_Close
               get_valid_L0_file
               validate_L0_header
               log_fmt_msg  

!Revision History:
               Revision 2.3  2000/07026
               John Seaton
               Added aqua compatability modifications.

               Revision 2.2 2000/03/21  11:52
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Added check for packet truncate return code from
               PGS_IO_L0_GetPacket. This will cause the PGE to continue
               processing instead of failing exit code 1.

               Revision 2.1  1997/08/27  10:50
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 2.0  1997/07/24  15:30
               Tom Johnson/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Updated from Version 1 to Version 2

               Revision 1.0  1997/06/18  16:40 EDT 
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               The details for the packet data locations were taken from 
               Santa Barbara Remote Sensing (SBRS) Contract Data Requirement 
               List (CDRL)305, MODIS Engineering Telemetry Description, 
               Tables T30-5A, T30-5B, T30-5C, T30-5D, and T30-5E.


!END*************************************************************************
*/

{
  /***************************************************************************/
  /*                                                                         */
  /*                   Declare Global Variables                              */
  /*                                                                         */
  /***************************************************************************/

  extern  PGSt_PC_Logical  global_L0_logical;
                                        /*  L0 logical file unit number      */

  /***************************************************************************/
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  static PGSt_SMF_status  lastStatus = MODIS_S_SUCCESS;

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  PGSt_SMF_status  L1A_status;     /* SMF-style message returned by function */
  PGSt_SMF_status  PGS_status;     /* SMF-style message returned by function */
  char             *routine = "read_a_packet";
  char             msg[300];
  PGSt_double      L0_start_time;
  PGSt_double      L0_stop_time;


  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  Set returnStatus equal to MODIS_S_SUCCESS                              */
  /*                                                                         */
  /*  set routine to "read_a_packet"    (done during declaration)            */
  /*                                                                         */
  /***************************************************************************/

  returnStatus = MODIS_S_SUCCESS;
  

  /***************************************************************************/
  /*                                                                         */
  /*  IF lastStatus is not equal to MODIS_S_SUCCESS                          */
  /*  THEN                                                                   */
  /*                                                                         */
  /*     Set returnStatus to lastStatus                                      */
  /*                                                                         */
  /*     Set msg to "no more packets in current L0 file"                     */
  /*     CALL log_fmt_msg to report that there are no more packets           */
  /*        in the L0 file                                                   */
  /*       INPUTS:  returnStatus, routine, msg                               */
  /*       OUTPUTS: None                                                     */
  /*       RETURN:  None                                                     */
  /*                                                                         */
  /***************************************************************************/

  if (lastStatus != MODIS_S_SUCCESS)
    {
     returnStatus = lastStatus;
     sprintf(msg, "L0 LUN: %d", global_L0_logical);
     log_fmt_msg (MODIS_W_NO_MORE_PACKETS, routine, msg);
    }


  /***************************************************************************/
  /*                                                                         */
  /*  ELSE                                                                   */
  /*                                                                         */
  /*     CALL PGS_IO_L0_GetPacket to retrieve a packet of data from the      */
  /*        L0 file                                                          */
  /*       INPUTS:  L0_file, pkt_buff_size                                   */
  /*       OUTPUTS: pkt                                                      */
  /*       RETURN:  returnStatus                                             */
  /*                                                                         */
  /***************************************************************************/

  else
    {
     PGS_status = PGS_IO_L0_GetPacket(*L0_file, PD_PKT_BUF_MAX, pkt);


     /************************************************************************/
     /*                                                                      */
     /*  SWITCH (PGS_status)                                                 */
     /*     CASE PGSIO_W_L0_END_OF_VIRTUAL_DS                                */
     /*                                                                      */
     /*        IF global_L0_logical equals PRIOR_L0_PCF_ID                   */
     /*        THEN                                                          */
     /*           CALL PGS_IO_L0_Close to close the prior L0 file            */
     /*             INPUTS:  L0_file                                         */
     /*             OUTPUTS: None                                            */
     /*             RETURN:  returnStatus                                    */
     /*                                                                      */
     /*           CALL PGS_SMF_TestSuccessLevel to determine if the prior    */
     /*              L0 file was closed                                      */
     /*             INPUTS:  returnStatus                                    */
     /*             OUTPUTS: None                                            */
     /*             RETURN:  result                                          */
     /*                                                                      */
     /*           IF result is equal to PGS_FALSE                            */
     /*           THEN                                                       */
     /*              set msg to "unable to close the prior L0 file"          */
     /*              CALL log_fmt_msg to report that the prior L0 file       */
     /*                 could not be closed                                  */
     /*                INPUTS:  returnStatus, routine, msg                   */
     /*                OUTPUTS: None                                         */
     /*                RETURN:  None                                         */
     /*           ENDIF                                                      */
     /*                                                                      */
     /*           Set global_L0_logical to CURRENT_L0_PCF_ID                 */
     /*                                                                      */
     /*         NOTE: Modis timestamps on both instruments, Terra and Aqua,  */
     /*         use segmented timestamps. Therefore, no change is needed for */
     /*         Aqua delivery to read packet times. The Terra toolkit macros */
     /*         PGSd_EOS_AM1 are used so L0 files are opened using the       */
     /*         toolkit routines which read segmented time.                  */
     /*                                                                      */
     /*           CALL get_valid_L0_file to open and validate the current    */
     /*              L0 file                                                 */
     /*             INPUTS:  spacecraft_tag                                  */
     /*             OUTPUTS: L0_file, L0_start_time, L0_stop_time            */
     /*             RETURN:  returnStatus                                    */
     /*                                                                      */
     /*           IF returnStatus is not equal to MODIS_S_SUCCESS            */
     /*           THEN                                                       */
     /*              Set msg to "unable to open the current L0 file"         */
     /*              CALL log_fmt_msg to report that the current L0 file     */
     /*                 could not be opened                                  */
     /*                INPUTS:  returnStatus, routine, msg                   */
     /*                OUTPUTS: None                                         */
     /*                RETURN:  None                                         */
     /*              Set returnStatus to MODIS_F_PKT_READ_FAILED             */
     /*           ENDIF                                                      */
     /*        ELSE                                                          */
     /*           set lastStatus to MODIS_W_NO_MORE_PACKETS                  */
     /*        ENDIF                                                         */
     /************************************************************************/

     switch (PGS_status)
       {
        case PGSIO_W_L0_END_OF_VIRTUAL_DS:
           if (global_L0_logical == PC_PRIOR_L0_PCF_ID)
             {
              PGS_status = PGS_IO_L0_Close (*L0_file);
              if (PGS_SMF_TestSuccessLevel(PGS_status) == PGS_FALSE)
                 log_fmt_msg (MODIS_E_PGS_IO_LO_CLOSE, routine, 
                              "Unable to close the prior L0 file");

              global_L0_logical = PC_CURRENT_L0_PCF_ID;
              L1A_status = get_valid_L0_file(PGSd_EOS_AM1,
                                             L0_file, 
                                             &L0_start_time,
                                             &L0_stop_time);

              if (L1A_status != MODIS_S_SUCCESS)
                {
                 sprintf(msg, "L0 LUN: %d", global_L0_logical);
                 log_fmt_msg (MODIS_E_GET_VALID_L0_FILE, routine, msg);
                 lastStatus = MODIS_F_PKT_READ_FAILED;
                }
             }
           else
             {
              lastStatus = MODIS_W_NO_MORE_PACKETS;
             }
           break;


     /************************************************************************/
     /*                                                                      */
     /*     CASE PGSIO_M_L0_HEADER_CHANGED                                   */
     /*        CALL validate_L0_header to validate the L0 header             */
     /*          INPUT:  L0_file                                             */
     /*          OUTPUT: None                                                */
     /*          RETURN: L1A_status                                          */
     /*                                                                      */
     /*        IF L1A_status is not equal to MODIS_S_SUCCESS                 */
     /*        THEN                                                          */
     /*           Set returnStatus to MODIS_F_PKT_READ_FAILED                */
     /*           Set msg to "unable to validate the L0 header               */
     /*              successfully"                                           */
     /*           CALL log_fmt_msg to report that the L0 header could not    */
     /*              be validated successfully                               */
     /*             INPUT:  L1A_status, routine, msg                         */
     /*             OUTPUT: None                                             */
     /*             RETURN: None                                             */
     /*        ENDIF                                                         */
     /*                                                                      */
     /************************************************************************/

        case PGSIO_M_L0_HEADER_CHANGED:
           L1A_status = validate_L0_header (*L0_file);
           if (L1A_status != MODIS_S_SUCCESS)
             {
              returnStatus = MODIS_F_PKT_READ_FAILED;
              sprintf(msg, "Unable to validate the L0 header successfully L0 LUN: %d",
                      global_L0_logical);
              log_fmt_msg (MODIS_F_PKT_READ_FAILED, routine, msg);
             }
           break;


     /************************************************************************/
     /*                                                                      */
     /*     CASE PGS_S_SUCCESS or PGSIO_W_L0_PKT_BUF_TRUNCATE                */
     /*      The return code PGSIO_W_L0_PKT_BUF_TRUNCATE is returned when    */
     /*      the Toolkit routine reads a packet that is too large to fit     */
     /*      in the allocated packet space. The oversize portion of the      */
     /*      packet is truncated. This prevents the PGE from exiting 1 as it */
     /*      would do if this return code went into the default case.        */
     /*        Break                                                         */
     /*                                                                      */
     /************************************************************************/

        case PGS_S_SUCCESS:
        case PGSIO_W_L0_PKT_BUF_TRUNCATE:
           break;


     /************************************************************************/
     /*                                                                      */
     /*      DEFAULT                                                         */
     /*         Set returnStatus to MODIS_F_PKT_READ_FAILED                  */
     /*   END-SWITCH                                                         */
     /*                                                                      */
     /************************************************************************/

        default:
           returnStatus = MODIS_F_PKT_READ_FAILED;
       }


  /***************************************************************************/
  /*                                                                         */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

    }

  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}
