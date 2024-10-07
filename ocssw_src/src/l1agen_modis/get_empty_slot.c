#include "L1A_prototype.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"

int16  get_empty_slot (void)

/*
!C************************************************************************

!Description:  This function searches the Vdata ID table for an empty slot in 
               which to insert a Vdata name/id pair.

!Input Parameters:
               None

!Output Parameters:
               None

Return Values: 
               FAIL                        (hdf.h)
               index                       **  the index of the empty slot  **
                                           **  if Vdata_id_table is full    **

Externally Defined:  
               int16                       (hdfi.h)
               TRUE                        (hdf.h)
               FALSE                       (hdf.h)
               VU_MAX_NUMBER_OF_VDATAS     (VU_vdata_utility.h)
               VU_EMPTY                    (VU_vdata_utility.h)
               global_VU_VDATA_ID          (level1a)

Called By:
               remember

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/10/02  12:20 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.1  1997/09/03  10:55
               Tom Johnson/GSC     (johnson@ltpmail.gsfc.nasa.gov)
               Incorporate walkthrough comments

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


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set index to 0                                                         */
  /*                                                                        */
  /* Set Found to FALSE                                                     */
  /*                                                                        */
  /**************************************************************************/

  int16 index = 0;
 
  int16 Found = FALSE;


  /**************************************************************************/
  /*                                                                        */
  /* DO-WHILE (index is less than or equal to VU_MAX_NUMBER_OF_VDATAS-1)    */
  /*          AND (Found is equal to FALSE)                                 */
  /*    IF (global_VU_VDATA_ID[index].id is equal to VU_EMPTY)              */
  /*    THEN                                                                */
  /*       Set Found to TRUE                                                */
  /*    ELSE                                                                */
  /*       Increment index                                                  */
  /*    ENDIF                                                               */
  /*                                                                        */
  /* END-WHILE                                                              */
  /*                                                                        */
  /**************************************************************************/

     while ((index < VU_MAX_NUMBER_OF_VDATAS) && (Found == FALSE)) {
        if (global_VU_VDATA_ID[index].id == VU_EMPTY)
           Found = TRUE;
        else
           index++;
     }


  /**************************************************************************/
  /*                                                                        */
  /* IF Found is equal to FALSE                                             */
  /* THEN                                                                   */
  /*    Set index to FAIL                                                   */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* RETURN index                                                           */
  /*                                                                        */
  /**************************************************************************/

     if (Found == FALSE)
        index = FAIL;

     return index;

 }  /* End of routine get_empty_slot */
