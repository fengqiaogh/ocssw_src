#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "hdf.h"
#include "hdfi.h"


PGSt_SMF_status  close_Vdata ( char *Vdata_name )

/*
!C************************************************************************

!Description:  This function detaches the specified Vdata from the 
               currently open hdf file.

!Input Parameters:
               char     *Vdata_name       ** The name of the Vdata to close **

!Output Parameters:
               None

Return Values:
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               FAIL                         (hdf.h)

Externally Defined:  
               int32                        (hdfi.h)

Called By:
               end_eng_data_access_to_file

Routines Called:
               recall_id
               VSdetach
               forget
               log_fmt_msg

!Revision History:
               Revision 2.0  1997/10/01  16:45 EDT
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
	       HDF portions developed at the National Center for 
               Supercomputing Applications at the University of Illinois 
               at Urbana-Champaign.               

!Design Notes:
               None

!END************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*          Declare the local variables and initialize them.              */
  /*                                                                        */
  /**************************************************************************/
  /*                                                                        */
  /* Set routine to "close_Vdata"                                           */
  /* Set returnStatus to MODIS_S_SUCCESS                                    */
  /* declare Vdata_id to be variables of type int32                         */
  /*                                                                        */
  /**************************************************************************/

  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */

  char  *routine = "close_Vdata";

  char  msg[300];

  int32  Vdata_id;


  returnStatus = MODIS_S_SUCCESS;


  /**************************************************************************/
  /*                                                                        */
  /* CALL recall_id to get the Vdata_id for this Vdata                      */
  /*   INPUTS:  Vdata_name                                                  */
  /*   OUTPUTS: None                                                        */
  /*   RETURNS: Vdata_id                                                    */
  /*                                                                        */
  /* IF Vdata_id is not equal to FAIL                                       */
  /* THEN                                                                   */
  /*    CALL VSdetach to detach the Vdata from the file                     */
  /*      INPUTS:  Vdata_id                                                 */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /*    CALL forget to remove the Vdata's name/id pair from the lookup table*/
  /*      INPUTS:  Vdata_name                                               */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /**************************************************************************/

     Vdata_id = recall_id(Vdata_name);
     if (Vdata_id != FAIL) {
        VSdetach(Vdata_id);
        forget(Vdata_name);
     }


  /**************************************************************************/
  /*                                                                        */
  /* ELSE                                                                   */
  /*    Set Status to MODIS_E_RECALL_ID                                     */
  /*    Set msg to "Vdata Name: <vdata name>"                               */
  /*    Set routine to "close_Vdata"                                        */
  /*    CALL log_fmt_msg to report that the Vdata id could not be retrieved */
  /*      INPUTS:  Status, routine, msg                                     */
  /*      OUTPUT:  None                                                     */
  /*      RETURN:  None                                                     */
  /*                                                                        */
  /*    Set returnStatus to FAIL                                            */
  /* ENDIF                                                                  */
  /*                                                                        */
  /* RETURN returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

     else {
        sprintf(msg, "Vdata Name: %s", Vdata_name);
        log_fmt_msg (MODIS_E_RECALL_ID, routine, msg);

        returnStatus = FAIL;
     }

     return returnStatus;

 } /* End of routine close_Vdata.c */
