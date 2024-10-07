#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "VU_vdata_utility.h"
#include "hdf.h"
#include "hdfi.h"

PGSt_SMF_status  end_Vdata_access_to_file (MODFILE  *L1A_file)

/*
!C************************************************************************

!Description:  This function checks to see that all of the Vdata have been
               detached from the file.

!Input Parameters:
               None

!Output Parameters:
               None

!Input/Output Parameters:
               MODFILE    *L1A_file         **  Level 1A file           **

Return Values:
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_ATTACHED_VDATAS      (PGS_MODIS_35005.h)
               FAIL                         (hdf.h)

Externally Defined:  
               PGSt_SMF_status              (PGS_SMF.h)
               int16                        (hdfi.h)
               global_H_ID                  (level1a)

Called By:
               end_eng_access_to_file

Routines Called:
               log_fmt_msg
               attached_Vdata_counter

!Revision History:
               Revision 2.0  1997/10/01  12:45 EDT
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

  extern   int32   global_H_ID;


  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set routine to "end_Vdata_access_to_file"                              */
  /*                                                                        */
  /* Set returnStatus to MODIS_S_SUCCESS                                    */
  /*                                                                        */
  /**************************************************************************/

  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */

  char  *routine = "end_Vdata_access_to_file";

  char  msg[300];

  int16  num_attached_Vdatas;


  returnStatus = MODIS_S_SUCCESS;

  
  /**************************************************************************/
  /*                                                                        */
  /*            CHECK TO SEE IF ALL Vdatas HAVE BEEN UNATTACHED             */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* CALL attached_Vdata_counter to get the number of attached Vdatas       */
  /*   INPUTS:  VU_REPORT                                                   */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: num_attached_Vdatas                                         */
  /*                                                                        */
  /* IF ( num_attached_Vdatas > 0 )                                         */
  /* THEN                                                                   */
  /*    set Status to MODIS_E_ATTACHED_VDATAS                               */
  /*    set msg to "Attempting to close the hdf file with Vdata(s) still    */
  /*                attached."                                              */
  /*    CALL log_fmt_msg to report that there are still Vdata attached to   */
  /*         the HDF file                                                   */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUTS: None                                                     */
  /*      RETURN:  None                                                     */
  /*                                                                        */
  /*    set returnStatus to FAIL                                            */
  /* ELSE                                                                   */
  /*    set global_H_ID to NULL                                             */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* return returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     num_attached_Vdatas = attached_Vdata_counter(VU_REPORT);

     if ( num_attached_Vdatas > 0 ) {
        sprintf(msg, "Number of Vdata(s) still attached: %d", num_attached_Vdatas);
        log_fmt_msg(MODIS_E_ATTACHED_VDATAS, routine, msg);

        returnStatus = FAIL;
     }
     else {
        L1A_file->hdf_id = global_H_ID;
	/*       global_H_ID = NULL;  */
     }

     return returnStatus;

 } /* End of routine end_Vdata_access_to_file */
