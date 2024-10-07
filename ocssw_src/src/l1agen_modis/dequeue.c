#include "FP_failed_pkt_queue.h"
#include "L1A_prototype.h"


char *dequeue ( FP_QUEUE_t  Q )

/*
!C************************************************************************

!Description:  Removes the first item from the front of the queue and
               returns it. If the queue is empty, NULL is returned.

!Input Parameters:
               None

!Output Parameters:
               None

!Input/Output Parameters:
               FP_QUEUE_t  Q                ** The FIFO queue from which to
                                               dequeue an item                **

Return Values: 
               NULL
               item                         ** the item that is first in
                                               line (or NULL if Q is empty)   **

Externally Defined:  
               FP_QUEUE_t                     (FP_failed_pkt_queue.h)

Called By:
               write_failed_packets

Routines Called:
               None

!Revision History:
               Revision 2.0  2001/01/04
               John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
               Simplified code. Removed call to is_empty

               Revision 1.0  1997/07/14  15:58 EDT
               David Catozzi/SAIC/GSC (cato@ltpmail.gsfc.nasa.gov)
               Originated code development.

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

  node_ptr  first_cell;
  char      *item=NULL;


  /***************************************************************************/
  /*                                                                         */
  /* IF (NOT((Q == NULL) OR (Q->head == NULL))                               */
  /* THEN                                                                    */
  /*   CALL is_empty to ascertain whether or not there are any items on the  */ 
  /*     queue.                                                              */
  /*    INPUTS:  Q                                                           */
  /*    OUTPUTS: None                                                        */
  /*    RETURNS: queue_is_empty                                              */
  /*                                                                         */
  /*   IF ( queue_is_empty equals True )                                     */
  /*   THEN                                                                  */
  /*     set item to NULL                                                    */
  /*   ELSE                                                                  */
  /*                                                                         */
  /*     set first_cell to Q.head.next                                       */
  /*     set item to first_cell.element                                      */
  /*     set Q.head.next to Q.head.next.next                                 */
  /*                                                                         */
  /*     free the dynamically allocated memory (first_cell)                  */
  /*                                                                         */
  /*     set Q.head.next.prev to Q.head                                      */
  /*   ENDIF                                                                 */
  /* ENDIF                                                                   */
  /*                                                                         */
  /***************************************************************************/

  if (!(Q == NULL) || (Q->head == NULL)) {

    if( Q->head->next == Q->tail ) /* check for empty queue */
       item = NULL;
    else {
       first_cell = Q->head->next;
       item = first_cell->element;
       Q->head->next = Q->head->next->next;

       free( first_cell );

       Q->head->next->prev = Q->head;
    }
  } 


  /***************************************************************************/
  /*                                                                         */
  /*  RETURN item                                                            */
  /*                                                                         */
  /***************************************************************************/

  return (item);

}
