#include "L1A_prototype.h"
#include "FP_failed_pkt_queue.h"
#include "PGS_MODIS_35005.h"
#include "PGS_SMF.h"
#include "hdf.h"
#include "mapiL1A.h"
#include "hfile.h"


PGSt_SMF_status  write_failed_packets (FP_QUEUE_t failed_pkts)

/*
!C*****************************************************************************

!Description:  This function writes the failed packets for the scan to the
               L1A file's 'Discarded Packets' Vdata.

!Input Parameters:
               FP_QUEUE_t    failed_pkts     ** The failed_pkts queue **

!Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               FAIL                         (HDF)

Externally Defined:
               FP_QUEUE_t                   (FP_failed_pkt_queue.h)
               PGSt_SMF_status              (PGS_SMF.h)
               MODIS_E_WRITE_VDATA          (PGS_MODIS_35005.h)
               FP_VDATA_MAX_BLOCK_NUM       (FP_failed_pkt_queue.h)

Called By:
               write_scan

Routines Called:
               dequeue
               write_Vdata
               log_fmt_msg

!Revision History:
               Revision 2.0  1999/08/20  10:47
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added call to free the packet memory.

               Revision 1.1  1997/09/03  10:55
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate walkthrough comments

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

!END***************************************************************************
*/

{
  /***************************************************************************/
  /*                                                                         */
  /*              Define Local Variables                                     */
  /*                                                                         */
  /***************************************************************************/

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  PGSt_SMF_status  tempStatus;     /* SMF-style message returned by function */
  char             *next_failed_packet;
  char             *routine = "write_failed_packets";
  extern int32     global_H_ID;
  filerec_t        *file_rec;


  /***************************************************************************/
  /*                                                                         */
  /*  set routine to "write_failed_packets"  (done during declaration)       */
  /*                                                                         */
  /*  set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                         */
  /***************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /***************************************************************************/
  /*                                                                         */
  /*  CALL dequeue to take the first packet off the failed_pkts queue        */
  /*  INPUTS:  failed_pkts                                                   */
  /*  OUTPUTS: failed_pkts                                                   */
  /*  RETURN:  next_failed_packet                                            */
  /*                                                                         */
  /***************************************************************************/

  next_failed_packet = dequeue(failed_pkts);


  /***************************************************************************/
  /*                                                                         */
  /*  Call HAatom_object to get the current number of vdata blocks used      */
  /*    INPUTS: global_H_ID                                                  */
  /*    OUTPUTS: None                                                        */
  /*    RETURNS: filerec_t structure (file_rec)                              */
  /*                                                                         */
  /*  DO_WHILE (next_failed_packet is not NULL)                              */
  /*                                                                         */
  /*     IF file_rec->maxref does not exceed the max number of vdata blocks  */
  /*       to be safely written (FP_VDATA_MAX_BLOCK_NUM)                     */
  /*                                                                         */
  /*       CALL write_Vdata to write 1 record (650 bytes) for each failed    */
  /*          packet to the Discarded Packets Vdata for this scan            */
  /*         INPUTS:  "Discarded Packets", next_failed_packet, 1             */
  /*         OUTPUTS: None                                                   */
  /*         RETURN:  tempStatus                                             */
  /*                                                                         */
  /*       IF (tempStatus is FAIL)                                           */
  /*       THEN                                                              */
  /*          set Status to MODIS_E_WRITE_VDATA                              */
  /*          CALL log_fmt_msg to report that the failed packet could not be */
  /*             written to the L1A granule                                  */
  /*            INPUTS:  Status, routine, msg                                */
  /*            OUTPUTS: None                                                */
  /*            RETURN:  None                                                */
  /*          set returnStatus to FAIL                                       */
  /*       ENDIF                                                             */
  /*     ENDIF                                                               */
  /*                                                                         */
  /***************************************************************************/

  file_rec = HAatom_object(global_H_ID);

  while (next_failed_packet != NULL)
    {
     if (file_rec->maxref <= FP_VDATA_MAX_BLOCK_NUM) {
       tempStatus = write_Vdata(M01DISCARDED_PACKETS, 
                                (unsigned char *) next_failed_packet, 1);

       if (tempStatus == FAIL)
         {
          log_fmt_msg(MODIS_E_WRITE_VDATA, routine, "Vdata Name: Discarded packets");
          returnStatus = FAIL;
         }
     }


  /***************************************************************************/
  /*                                                                         */
  /*     free the space for the dequeued failed packet (next_failed_packet)  */
  /*                                                                         */
  /*      -----------------------------------------------------------        */
  /*                                                                         */
  /*     CALL dequeue to take the first packet off the failed_pkts queue     */
  /*       INPUTS:  failed_pkts                                              */
  /*       OUTPUTS: failed_pkts                                              */
  /*       RETURN:  next_failed_packet                                       */
  /*                                                                         */
  /*  END_WHILE                                                              */
  /*                                                                         */
  /***************************************************************************/

     free(next_failed_packet);

     next_failed_packet = dequeue(failed_pkts);

    }


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN returnStatus                                                    */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}

