#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "FP_failed_pkt_queue.h"


FP_QUEUE_t  make_queue (void)

/*
!C************************************************************************

!Description:  Creates a new linked list FIFO queue.

!Input Parameters:
               None

!Output Parameters:
               None

Return Values:
               FP_QUEUE_t  Q               ** The newly created FIFO queue **
               MODIS_E_MALLOC_FAILED       (PGS_MODIS_35005.t)

Externally Defined:  
               FP_QUEUE_t                  (FP_failed_pkt_queue.h)

Called By:
               level1a

Routines Called:
               log_fmt_msg

!Revision History:
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

!END************************************************************************
*/

{

  /**************************************************************************/
  /*              Declare and Initialize Local Variables                    */
  /*                                                                        */
  /**************************************************************************/

  FP_QUEUE_t  Q;

  char        *routine = "make_queue";

  /**************************************************************************/

  /**************************************************************************/ 
  /*                                                                        */
  /*  set routine to "make_queue"  (done during declaration)                */
  /*                                                                        */
  /**************************************************************************/ 


  /**************************************************************************/ 
  /*                                                                        */
  /*  dynamically allocate memory for the queue (Q)                         */
  /*  IF ( Q equals NULL )                                                  */
  /*  THEN                                                                  */
  /*     CALL log_fmt_msg to report the memory allocation error             */
  /*       INPUTS:  Status, routine, msg                                    */
  /*       OUTPUTS: None                                                    */
  /*       RETURN:  None                                                    */
  /*                                                                        */
  /**************************************************************************/ 

  Q = (struct queue_record *) malloc( sizeof (struct queue_record) );
  if (Q == NULL)
     log_fmt_msg (MODIS_E_MALLOC_FAILED, routine, 
           "unable to allocate memory for the failed packets queue");
  

  /**************************************************************************/ 
  /*                                                                        */
  /*  ELSE                                                                  */
  /*     dynamically allocate memory for the queue (Q.head)                 */
  /*     IF ( Q.head equals NULL )                                          */
  /*     THEN                                                               */
  /*        set Status to MODIS_E_MALLOC_FAILED                             */
  /*        CALL log_fmt_msg to report the memory allocation error          */
  /*          INPUTS:  Status, routine, msg                                 */
  /*          OUTPUTS: None                                                 */
  /*          RETURN:  None                                                 */
  /*        set Q to NULL                                                   */
  /*                                                                        */
  /**************************************************************************/ 

  else
    {
     Q->head = (node_ptr) malloc( sizeof (struct node) );
  
     if ( Q->head == NULL) 
       {
        log_fmt_msg (MODIS_E_MALLOC_FAILED, routine, 
               "unable to allocate memory for the failed packet queue's head");
        Q = NULL;
       }


     /***********************************************************************/ 
     /*                                                                     */
     /*  ELSE                                                               */
     /*     dynamically allocate memory for the queue (Q.tail)              */
     /*     IF ( Q.head.next equals NULL )                                  */
     /*     THEN                                                            */
     /*        set Status to MODIS_E_MALLOC_FAILED                          */
     /*        CALL log_fmt_msg to report the memory allocation error       */
     /*          INPUTS:  Status, routine, msg                              */
     /*          OUTPUTS: None                                              */
     /*          RETURN:  None                                              */
     /*        set Q to NULL                                                */
     /*                                                                     */
     /***********************************************************************/ 

     else
       {
        Q->head->next = (node_ptr) malloc( sizeof (struct node));
        if (Q->head->next == NULL)
          {
           log_fmt_msg (MODIS_E_MALLOC_FAILED, routine, 
                  "unable to allocate memory for the failed packet queue's tail");
           Q = NULL;
          }


        /********************************************************************/ 
        /*                                                                  */
        /*  ELSE                                                            */
        /*     set Q.tail to Q.head.next                                    */
        /*     set Q.tail.prev to Q.head                                    */
        /*                                                                  */
        /********************************************************************/ 

        else
          {
           Q->tail = Q->head->next;
           Q->tail->prev = Q->head;


  /**************************************************************************/ 
  /*                                                                        */
  /*        ENDIF                                                           */
  /*     ENDIF                                                              */
  /*  ENDIF                                                                 */
  /*                                                                        */
  /*  return Q                                                              */
  /*                                                                        */
  /**************************************************************************/ 

          }
       }
    }

  return (Q);
}
