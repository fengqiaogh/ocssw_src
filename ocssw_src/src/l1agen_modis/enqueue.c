#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "FP_failed_pkt_queue.h"


PGSt_SMF_status  enqueue ( char *item,  FP_QUEUE_t Q )

/*
!C************************************************************************

!Description:  Puts the item on the end of the queue.

!Input Parameters:
               char          *item         **  The item to put on the   **
                                           **  end of the queue         **

!Output Parameters:
               None

!Input/Output Parameters:
               FP_QUEUE_t    Q             **  The FIFO queue on which  **
                                           **  to put the item          **

Return Values:
               MODIS_S_SUCCESS             (PGS_MODIS_35005.h)
               MODIS_E_MALLOC_FAILED       (PGS_MODIS_35005.h)

Externally Defined:  
               PGSt_SMF_status             (PGS_SMF.h)
               FP_QUEUE_t                  (FP_failed_pkt_queue.h)

Called By:
               accumulate_failed_packets

Routines Called:
               log_fmt_msg

!Revision History:
               Revision 1.0  1997/09/15  10:35 
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Original code.

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
  char             *routine = "enqueue";
  node_ptr         tmp_cell; 

  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  set routine to "enqueue" (done during declaration)                     */
  /*  set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                         */
  /***************************************************************************/

   returnStatus = MODIS_S_SUCCESS;
 

  /***************************************************************************/
  /*                                                                         */
  /*  Dynamically allocate memory for the next item (tmp_cell) in the queue  */
  /*                                                                         */
  /*  IF unable to allocate memory                                           */ 
  /*  THEN                                                                   */
  /*     Set msg to "unable to allocate memory"                              */
  /*     Set returnStatus to MODIS_E_MALLOC_FAILED                           */
  /*     CALL log_fmt_msg to report the memory allocation error              */
  /*       INPUTS:  logStatus, routine, msg                                  */
  /*       OUTPUTS: None                                                     */
  /*       RETURN:  None                                                     */
  /*                                                                         */
  /*  ELSE                                                                   */
  /*     Set tmp_cell.element to item                                        */
  /*     Set tmp_cell.prev to Q.tail.prev                                    */
  /*     Set Q.tail.prev.next to tmp_cell                                    */
  /*     Set Q.tail.prev to tmp_cell                                         */
  /*     Set tmp_cell.next to Q.tail                                         */
  /*                                                                         */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  tmp_cell = (node_ptr) malloc( sizeof ( struct node ) );
  if (tmp_cell == NULL) 
    {
     returnStatus = MODIS_E_MALLOC_FAILED;
     log_fmt_msg (MODIS_E_MALLOC_FAILED, routine, 
                  "Unable to allocate the failed packets linked list node");
    }

  else
    {
     tmp_cell->element = item;
     tmp_cell->prev = Q->tail->prev;
     Q->tail->prev->next = tmp_cell;
     Q->tail->prev = tmp_cell;
     tmp_cell->next = Q->tail;
    }


  /***************************************************************************/
  /*                                                                         */
  /* RETURN returnStatus                                                     */
  /*                                                                         */
  /***************************************************************************/

  return (returnStatus);

}
