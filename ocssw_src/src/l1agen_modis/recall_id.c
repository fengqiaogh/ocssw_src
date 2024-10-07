#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"
#include "VU_vdata_utility.h"

int32  recall_id ( char  *Vdata_name )

/*
!C************************************************************************

!Description:  This function gets the Vdata id number corresponding to the 
               specified Vdata name.

!Input Parameters:
               char     *Vdata_name        ** The name of the Vdata for which
                                              to recall the id              **

!Output Parameters:
               None

Return Values: 
               FAIL                              (hdf.h)
               Vdata_id                    ** The Vdata id corresponging to **
                                           ** the specified Vdata (or FAIL) **

Externally Defined:  
               PGS_SMF_status                    (PGS_SMF.h)
               int16                             (hdfi.h)
               int32                             (hdfi.h)
               TRUE                              (hdf.h)
               FALSE                             (hdf.h)
               global_VU_ID_TABLE_READY          (level1a)
               VU_MAX_NUMBER_OF_VDATAS           (VU_vdata_utility.h)
               MODIS_W_ID_TABLE_NOT_INIT         (PGS_MODIS_35005.h)
               global_VU_VDATA_ID                (level1a)

Called By:
               write_Vdata
               close_Vdata

Routines Called:
               log_fmt_msg
               equal_strings

!Revision History:
               Revision 2.0  1997/10/01  22:05 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.

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
  /* Set routine to "recall_id"                                             */
  /*                                                                        */
  /* set returnStatus to FAIL                                               */
  /*                                                                        */
  /**************************************************************************/

  char  *routine = "recall_id";           /* Variable to hold routine name  */

  int32  returnStatus;          /* int32-style message returned by function */

  int16  strings_are_equal;

  int ctr;

  int Found;


  returnStatus = FAIL;



  /**************************************************************************/
  /*                                                                        */
  /* {first ensure that the id table has been initialized}                  */
  /* IF ( global_VU_ID_TABLE_READY equals False )                           */
  /* THEN                                                                   */
  /*    CALL log_fmt_msg to report that the id table is empty               */
  /*      INPUTS:  Status, routine, ""                                      */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /**************************************************************************/
  
     if (!global_VU_ID_TABLE_READY)
        log_fmt_msg(MODIS_W_ID_TABLE_NOT_INIT, routine, "");
   
  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    set Found to FALSE                                                  */
  /*    set ctr to 0                                                        */
  /*    DO-WHILE (ctr is less than or equal to (VU_MAX_NUMBER_OF_VDATAS-1)) */
  /*             AND (Found is equal to FALSE)                              */
  /*       CALL equal_strings to test for string equality                   */
  /*         INPUTS:  VDATA_ID[ctr].name, Vdata_name                        */
  /*         OUTPUTS: None                                                  */
  /*         RETURN:  strings_are_equal                                     */
  /*                                                                        */
  /*       IF (strings_are_equal equals True)                               */
  /*       THEN                                                             */
  /*          set Found to TRUE                                             */
  /*       ELSE                                                             */
  /*          Increment ctr                                                 */
  /*       ENDIF                                                            */
  /*    END-WHILE                                                           */
  /*                                                                        */
  /*    IF (Found is equal to TRUE)                                         */
  /*    THEN                                                                */
  /*       set returnStatus to global_VU_VDATA_ID[ctr].id                   */
  /*    ENDIF                                                               */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* RETURN returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     else {
        Found = FALSE;
        ctr   = 0;
        while ((ctr < VU_MAX_NUMBER_OF_VDATAS) && (Found == FALSE)) {

           strings_are_equal = equal_strings(global_VU_VDATA_ID[ctr].name, 
                                             Vdata_name);

           if (strings_are_equal == TRUE) 
              Found = TRUE;
           else
              ctr++;
        }

        if (Found == TRUE)
           returnStatus = global_VU_VDATA_ID[ctr].id;
     }

     return returnStatus;

 } /* End of routine recall_id */
