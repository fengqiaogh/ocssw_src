#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "FP_failed_pkt_queue.h"


void   free_queue (FP_QUEUE_t  Q)

/*
!C************************************************************************

!Description:  Disposes of the queue, releasing dynamically allocated 
               memory.

!Input Parameters:
               FP_QUEUE_t  Q          ** The FIFO queue to dispose of **

!Output Parameters:
               None

!Input/Output Parameters:
               None

Return Values:
               None

Externally Defined:  
               FP_QUEUE_t             (FP_failed_pkt_queue.h)
               MODIS_E_FREE_QUEUE     (PGS_MODIS_35005.h)
               MODIS_E_NULL_POINTER   (PGS_MODIS_35005.h)
 
Called By:
               level1a

Routines Called:
               log_fmt_msg

!Revision History:
               Revision 2.0  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Simplified code. Removed call to is_empty.
               Added checks for NULL input arguments.

               Revision 1.0  1997/09/15   9:55
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Original code development.

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
  /*                                                                         */
  /*              Declare and Initialize Local Variables                     */
  /*                                                                         */
  /***************************************************************************/

  char             *routine = "free_queue";

  /***************************************************************************/
  /*                                                                         */
  /* set routine to "free_queue" (done during declaration)                   */
  /*                                                                         */
  /* IF ((Q not NULL) AND (Q->head not NULL))                                */ 
  /* THEN                                                                    */
  /*   IF (NOT(Q->head->next == Q->tail))                                    */
  /*   THEN                                                                  */
  /*      CALL log_fmt_msg to report that the queue is not empty             */
  /*        INPUTS:  Status, routine,                                        */
  /*         "There are still unwritten packets in the Failed Packets Queue" */
  /*        OUTPUTS: None                                                    */
  /*        RETURNS: None                                                    */
  /*   ENDIF                                                                 */
  /*                                                                         */
  /*   Free the dynamically allocated memory for Q.head                      */
  /*   Free the dynamically allocated memory for Q.tail                      */
  /*   Free the dynamically allocated memory for Q                           */
  /* ENDIF                                                                   */
  /*                                                                         */
  /***************************************************************************/

    if ((Q != NULL) && (Q->head != NULL)) {
       if (!(Q->head->next == Q->tail))
         log_fmt_msg (MODIS_E_FREE_QUEUE, routine, 
            "There are still unwritten packets in the Failed Packets Queue");

       free(Q->head);
       free(Q->tail);
       free(Q);
    }

  /***************************************************************************/
  /*                                                                         */
  /* ELSE                                                                    */
  /*   Figure out which argument was NULL and report error in LogStatus file */
  /* ENDELSE                                                                 */
  /*                                                                         */
  /***************************************************************************/

    else {
      if (Q == NULL)
        log_fmt_msg(MODIS_E_NULL_POINTER, routine,
                "The Q structure was NULL");
      else
        log_fmt_msg(MODIS_E_NULL_POINTER, routine, 
                 "The Q->head structure was NULL");          
   }

}
