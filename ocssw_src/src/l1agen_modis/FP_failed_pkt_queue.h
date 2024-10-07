#ifndef   FP_FAILED_PKT_QUEUE_H
#define   FP_FAILED_PKT_QUEUE_H

/*
!C-INC************************************************************************

!Description:  This header file contains the definitions and structures used
               for processing failed packets to be put into the Discarded
               Packets (section 5 of the MODIS Level 1A Data Product Format).

!Input Parameters:
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:  
               N/A

Called By:
               N/A

Routines Called:
               N/A

!Revision History:
               Revision 1.0  1996/09/24  14:45 EDT
               David Catozzi/GSC (cato@ltpmail.gsfc.nasa.gov)
               Created original header file

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               The ".h" file below was specifically written for development
               in C. Any other language choice may require reworking of the
               ".h" file before coding can begin.

!END**********************************************************************
*/

#include <stdio.h>
#include <stdlib.h>

/***********************************************************/
/*          Constants                                      */
/***********************************************************/

#define FP_VDATA_MAX_BLOCK_NUM  64400

/***********************************************************
 *          data  types                                    *
 ***********************************************************/

typedef struct node *node_ptr;

struct node
{
    char *element;
    node_ptr next;    
    node_ptr prev;
};


struct queue_record
{
   node_ptr head;
   node_ptr tail;
};

typedef struct queue_record *FP_QUEUE_t; 


#endif
