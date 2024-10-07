#include "L1A_prototype.h"
#include "EN_eng_data.h"
#include "hdfi.h"

 void  reset_last_valid_scan ( EN_VDATA_TYPE_t  *eng_data )

/*
!C************************************************************************

!Description:  This function resets the LAST_VALID_SCAN field for each Vdata
               in the eng_data array structure. All LAST_VALID_SCAN fields 
               that aren't EN_INITIAL_LAST_VALID_SCAN_VALUE are reset to 
               EN_INITIAL_FIELD_VALUE.

!Input Parameters:
               None

!Output Parameters:
               None

!Input/Output Parameters:
               EN_VDATA_TYPE_t  eng_data   ** The Vdata array structure (All **
                          	           ** LAST_VALID_SCAN fields that    **
                          		   ** aren't EN_INITIAL_LAST_VALID_  **
                                           ** SCAN_VALUE are reset to        **
                                           ** EN_INITIAL_FIELD_VALUE)        **
Return Values: 
               None

Externally Defined:
               EN_VDATA_TYPE_t                       (EN_eng_data.h)
               EN_NUM_VDATAS                         (EN_eng_data.h)
               EN_INITIAL_FIELD_VALUE                (EN_eng_data.h)
               EN_INITIAL_LAST_VALID_SCAN_VALUE      (EN_eng_data.h)
               uint16                                (hdfi.h)

Called By:     
               process_a_granule

Routines Called:
               None

!Revision History:
               Revision 2.0  1997/09/18  13:50 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

               Revision 1.0  1997/07/16  15:58 EDT
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
               This code was developed in C language.

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/

  uint16  v;                                  /* index to access eng fields */


  /**************************************************************************/  
  /*                                                                        */
  /* DO FOR ( v = 0 to EN_NUM_VDATAS-1 )                                    */
  /*    IF ( eng_data[v].field[0].value does NOT equal                      */
  /*                  EN_INITIAL_LAST_VALID_SCAN_VALUE )                    */
  /*    THEN                                                                */
  /*       set eng_data[v].field[0].value to EN_INITIAL_FIELD_VALUE         */
  /*    ENDIF                                                               */
  /* END_DO                                                                 */
  /*                                                                        */
  /**************************************************************************/  

  for (v = 0; v < EN_NUM_VDATAS; v++)
     if (eng_data[v].field[0].value != EN_INITIAL_LAST_VALID_SCAN_VALUE)
        eng_data[v].field[0].value = EN_INITIAL_FIELD_VALUE;

 }  /* End of routine reset_last_valid_scan */
