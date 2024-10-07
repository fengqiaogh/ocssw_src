#include "L1A_prototype.h"
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdfi.h"
#include "PD_pkt_data.h"
#include "FP_failed_pkt_queue.h"

PGSt_SMF_status  accumulate_failed_packets (PGSt_IO_L0_Packet  *pkt, 
                                            FP_QUEUE_t         failed_pkts)

/*
!C************************************************************************

!Description:  This function puts the failed packet "pkt" on the failed
               packets queue "failed_pkts".

!Input Parameters:
               PGSt_IO_LO_Packet  *pkt         ** The MODIS packet     **

!Output Parameters:
               None

!Input/Output Parameters:
               FP_QUEUE_t         failed_pkts  ** The failed pkt queue ** 

Return Values:
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_ENQUEUE              (PGS_MODIS_35005.h)
               MODIS_E_MALLOC_FAILED        (PGS_MODIS_35005.h)

Externally Defined:  
               FP_QUEUE_t                   (FP_failed_pkt_queue.h)
               PGSt_IO_L0_Packet            (PGS_IO.h)
               PD_PKT_BUF_MAX               (PD_pkt_data.h)
               PGSt_SMF_status              (PGS_SMF.h)
               int8                         (hdfi.h)

Called By:
               process_a_scan

Routines Called:
               log_fmt_msg
               enqueue

!Revision History:
               Revision 2.0  1997/08/27  10:00
               Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate PDL walkthru comments

               Revision 1.0  1997/07/14  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Original design.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               None

!END************************************************************************
*/

{
  /***************************************************************************/
  /*              Declare and Initialize Local Variables                     */
  /***************************************************************************/

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  char             *routine = "accumulate_failed_packets";
  int8             *copy_of_packet;
  int              i;

  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  Set returnStatus equal to MODIS_S_SUCCESS                              */
  /*                                                                         */
  /*  set routine to "accumulate_failed_packets" (done during declaration)   */
  /*                                                                         */
  /***************************************************************************/

  returnStatus = MODIS_S_SUCCESS;

  
  /***************************************************************************/
  /*                                                                         */
  /*  dynamically allocate 650 bytes of memory for the failed packet         */
  /*                                                                         */
  /*  IF ( unable to allocate memory for the failed packet )                 */
  /*  THEN                                                                   */
  /*                                                                         */
  /*     Set returnStatus to MODIS_E_MALLOC_FAILED                           */
  /*     CALL log_fmt_msg to report that memory could not be allocated for   */
  /*        the failed packet                                                */
  /*       INPUTS:  MODIS_E_MALLOC_FAILED, msg, "memory allocation error"    */
  /*       OUTPUTS: None                                                     */
  /*       RETURN:  None                                                     */
  /*                                                                         */
  /***************************************************************************/

  copy_of_packet = (int8 *) malloc(650);

  if (copy_of_packet == NULL)
    {
     returnStatus = MODIS_E_MALLOC_FAILED;
     log_fmt_msg (MODIS_E_MALLOC_FAILED, routine, 
          "Need an additional 650 bytes to store the failed packet data");
    }


  /***************************************************************************/
  /*                                                                         */
  /*  ELSE                                                                   */
  /*     DO FOR ( i = 0 to 649 )                                             */
  /*        IF (i < PD_PKT_BUF_MAX)                                          */
  /*        THEN                                                             */
  /*           Set copy_of_packet[i] to pkt[i]                               */
  /*        ELSE                                                             */
  /*           Set copy_of_packet[i] to 0                                    */
  /*        ENDIF                                                            */
  /*     END DO                                                              */
  /*                                                                         */
  /*     CALL enqueue to put the failed packet on the queue                  */
  /*       INPUTS:  copy_of_packet, failed_pkts                              */
  /*       OUTPUTS: None                                                     */
  /*       RETURN:  returnStatus                                             */
  /*                                                                         */
  /*     IF ( returnStatus is MODIS_E_MALLOC_FAILED )                        */
  /*     THEN                                                                */
  /*        Set msg to "failed packet could not be put in the queue"         */
  /*        Set logStatus to MODIS_E_ENQUEUE                                 */
  /*        CALL log_fmt_msg to report that the failed packet could not be   */
  /*           put in the queue                                              */
  /*          INPUTS:  logStatus, routine, msg                               */
  /*          OUTPUTS: None                                                  */
  /*          RETURN:  None                                                  */
  /*     ENDIF                                                               */
  /*                                                                         */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  else
    {
     memset (copy_of_packet, 0, 650);
     for (i=0; i<PD_PKT_BUF_MAX; i++)
        copy_of_packet[i] = pkt[i];

     returnStatus = enqueue ((char *)copy_of_packet, failed_pkts);

     if (returnStatus == MODIS_E_MALLOC_FAILED)
        log_fmt_msg (MODIS_E_ENQUEUE, routine, 
                     "Could not allocate the failed packets linked list node");
    }


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}
