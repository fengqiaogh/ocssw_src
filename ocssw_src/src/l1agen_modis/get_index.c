#include "L1A_prototype.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"

int16  get_index (char  *Vdata_name)

/*
!C************************************************************************

!Description:  This function retrieves the Vdata_id_table index for Vdata_name.

!Input Parameters:
               char  *Vdata_name         ** The name of the Vdata for     **
                                         ** which to get the associated   **
                                         ** Vdata_id_table index          **
!Output Parameters:
               None

Return Values: 
               FAIL                        (hdf.h)
               index                       ** if found in Vdata_id_table    **
 
Externally Defined:  
               int16                       (hdfi.h)
               VU_MAX_NUMBER_OF_VDATAS     (VU_vdata_utility.h)
               global_VU_VDATA_ID          (level1a)

Called By:
               forget

Routines Called:
               equal_strings

!Revision History:
               Revision 2.0  1997/10/02  12:55 EDT
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
  /* set index to 0                                                         */
  /* set Found to FALSE                                                     */
  /*                                                                        */
  /**************************************************************************/

  int16  index;
  int16  Found;
  int16  result;


  index = 0;
  Found = FALSE;


  /**************************************************************************/
  /*                                                                        */
  /* DO-WHILE (index is less than or equal to 0 to                          */
  /*          VU_MAX_NUMBER_OF_VDATAS-1) AND (Found is equal to FALSE)      */
  /*    CALL equal_strings to test string equality                          */
  /*      INPUTS:  global_VU_VDATA_ID[index].name, Vdata_name               */
  /*      OUTPUTS: None                                                     */
  /*      RETURN:  result                                                   */
  /*                                                                        */
  /*    IF (result is equal to TRUE)                                        */
  /*    THEN                                                                */
  /*       set Found to TRUE                                                */
  /*    ELSE                                                                */
  /*       increment index                                                  */
  /*    ENDIF                                                               */
  /*                                                                        */
  /* END_DO                                                                 */
  /*                                                                        */
  /**************************************************************************/

     while (( index < VU_MAX_NUMBER_OF_VDATAS) && (Found == FALSE)) {
        result = equal_strings(global_VU_VDATA_ID[index].name, Vdata_name);
        
        if (result == TRUE)
           Found = TRUE;
        else
           index++;
     }



  /**************************************************************************/
  /*                                                                        */
  /* IF (Found is equal to FALSE)                                           */
  /* THEN                                                                   */
  /*    set index to FAIL                                                   */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* RETURN index                                                           */
  /*                                                                        */
  /**************************************************************************/

     if (Found == FALSE)
        index = FAIL;

     return index;

 }  /* End of routine get_index */
