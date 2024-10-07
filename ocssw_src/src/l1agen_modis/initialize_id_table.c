#include "L1A_prototype.h"
#include "VU_vdata_utility.h"

void   initialize_id_table  (void)

/*
!C************************************************************************

!Description:  This function prepares the Vdata ID table for use.

!Input Parameters:
               None

!Output Parameters:
               None

Return Values: 
               None

Externally Defined:  
               VU_MAX_NUMBER_OF_VDATAS        (VU_vdata_utility.h)
               VU_EMPTY                       (VU_vdata_utility.h)
               global_VU_ID_TABLE_READY       (level1a)
               global_VU_VDATA_ID             (level1a)

Called By:
               remember

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/10/02  13:35 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

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
               This is a private function - only to be used by functions 
               internal to the Vdata_id_table.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                  Declare the global variables.                         */
  /*                                                                        */
  /**************************************************************************/

  extern  VU_ID_TABLE  global_VU_VDATA_ID[VU_MAX_NUMBER_OF_VDATAS];

  extern   int   global_VU_ID_TABLE_READY;


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set index to 0                                                         */
  /*                                                                        */
  /**************************************************************************/

  int16 index = 0;


  /**************************************************************************/
  /*                                                                        */
  /* DO_FOR ( index = 0 to VU_MAX_NUMBER_OF_VDATAS-1 )                      */
  /*    set global_VU_VDATA_ID[index].id to VU_EMPTY                        */
  /*    set global_VU_VDATA_ID[index].name to "EMPTY_SLOT"                  */
  /* END_DO                                                                 */
  /*                                                                        */
  /* set global_VU_ID_TABLE_READY to TRUE                                   */
  /*                                                                        */
  /**************************************************************************/

     for (index = 0; index < VU_MAX_NUMBER_OF_VDATAS; index++) {
        global_VU_VDATA_ID[index].id = VU_EMPTY;
        strcpy(global_VU_VDATA_ID[index].name, EN_EMPTY_SLOT);
     }

     global_VU_ID_TABLE_READY = TRUE;

 } /* End of routine initialize_id_table */
