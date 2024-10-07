#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "EN_eng_data.h"
#include "hdf.h"
#include "mapiL1A.h"


PGSt_SMF_status  end_eng_data_access_to_file (MODFILE          *L1A_file,
                                              EN_VDATA_TYPE_t  *eng_data)

/*
!C************************************************************************

!Description:  This function is a driver for other Vdata cleanup functions.

!Input Parameters:
               EN_VDATA_TYPE_t  *eng_data   ** The Vdata array structure **

!Output Parameters:
               None

!Input/Output Parameters:
               MODFILE          *L1A_file   ** The Level 1A file         **

Return Values:
               MODIS_S_SUCCESS                  (PGS_MODIS_35005.h)
               FAIL                             (hdf.h)

Externally Defined:
               PGSt_SMF_status                  (PGS_SMF.h)
               EN_NUM_VDATAS                    (EN_eng_data.h)
               EN_VDATA_TYPE_t                  (EN_eng_data.h)
               MODIS_E_CLOSE_VDATA              (PGS_MODIS_35005.h)
               MODIS_E_END_VDATA_FILE           (PGS_MODIS_35005.h)

Called By:
               process_a_granule

Routines Called:
               close_Vdata
               end_Vdata_access_to_file
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

  PGSt_SMF_status  tempStatus;    /* SMF-style message returned by function */

  char    *routine = "end_eng_data_access_to_file";

  char    msg[300];

  int     i;
  
  /**************************************************************************/
  /*                                                                        */
  /*  set routine to "end_eng_data_access_to_file" (done during declaration)*/
  /*                                                                        */
  /*  set returnStatus to MODIS_S_SUCCESS                                   */
  /*                                                                        */
  /**************************************************************************/

  returnStatus = MODIS_S_SUCCESS;
  

  /**************************************************************************/
  /*                                                                        */
  /*  DO FOR ( i = 0 to EN_NUM_VDATAS-1 )                                   */
  /*                                                                        */
  /*     set Vdata_name to eng_data[i].vdata_name                           */
  /*                                                                        */
  /*     CALL close_Vdata to detach it from the file                        */
  /*       INPUTS:  Vdata_name                                              */
  /*       OUTPUTS: None                                                    */
  /*       RETURN:  tempstatus                                              */
  /*                                                                        */
  /*     IF (tempstatus is FAIL)                                            */
  /*     THEN                                                               */
  /*        set Status to MODIS_E_CLOSE_VDATA                               */
  /*        set msg to "close_Vdata failed (couldn't detach an eng_data     */
  /*           Vdata)"                                                      */
  /*        CALL log_fmt_msg to report that the eng_data Vdata could not    */
  /*           be detached from the HDF file                                */
  /*          INPUTS:  Status, routine, msg                                 */
  /*          OUTPUTS: None                                                 */
  /*          RETURN:  None                                                 */
  /*        set returnStatus to FAIL                                        */
  /*     ENDIF                                                              */
  /*  END DO                                                                */
  /*                                                                        */
  /**************************************************************************/

  for (i=0; i < EN_NUM_VDATAS; i++)
    {
     tempStatus = close_Vdata(eng_data[i].vdata_name);
     if (tempStatus == FAIL)
       {
        sprintf(msg, "Vdata Name: %s", eng_data[i].vdata_name);
        log_fmt_msg(MODIS_E_CLOSE_VDATA, routine, msg);
        returnStatus = FAIL;
       }
    }


  /**************************************************************************/
  /*                                                                        */
  /*  set Vdata_name to "Discarded Packets"                                 */
  /*                                                                        */
  /*  CALL close_Vdata to detach it from the file                           */
  /*    INPUTS:  Vdata_name                                                 */
  /*    OUTPUTS: None                                                       */
  /*    RETURN:  tempstatus                                                 */
  /*                                                                        */
  /*  IF (tempstatus is FAIL)                                               */
  /*  THEN                                                                  */
  /*     set Status to MODIS_E_CLOSE_VDATA                                  */
  /*     set msg to "close_Vdata failed (couldn't detach the                */
  /*        discarded_packets Vdata)"                                       */
  /*     CALL log_fmt_msg report that the Discarded Packets Vdata could     */
  /*        not be detached from the HDF file                               */
  /*       INPUTS:  Status, routine, msg                                    */
  /*       OUTPUTS: None                                                    */
  /*       RETURN:  None                                                    */
  /*     set returnStatus to FAIL                                           */
  /*  ENDIF                                                                 */
  /*                                                                        */
  /**************************************************************************/

  tempStatus = close_Vdata(M01DISCARDED_PACKETS);
  if (tempStatus == FAIL)
    {
     log_fmt_msg(MODIS_E_CLOSE_VDATA, routine, "Vdata Name: discarded_packets");
     returnStatus = FAIL;
    }


  /**************************************************************************/
  /*                                                                        */
  /*  CALL end_Vdata_access_to_file to do Vdata cleanup                     */
  /*    INPUTS:  None                                                       */
  /*    OUTPUTS: None                                                       */
  /*    RETURN:  tempstatus                                                 */
  /*                                                                        */
  /*  IF (tempstatus is FAIL)                                               */
  /*  THEN                                                                  */
  /*     set Status to MODIS_E_END_VDATA_FILE                               */
  /*     set msg to "end_vdata_access_to_file failed (couldn't do Vdata     */
  /*        cleanup)"                                                       */
  /*     CALL log_fmt_msg to report that the Vdata cleanup failed           */
  /*       INPUTS:  Status, routine, msg                                    */
  /*       OUTPUTS: None                                                    */
  /*       RETURN:  None                                                    */
  /*     set returnStatus to FAIL                                           */
  /*  ENDIF                                                                 */
  /*                                                                        */
  /**************************************************************************/

  tempStatus = end_Vdata_access_to_file(L1A_file);

  if (tempStatus == FAIL)
    {
     log_fmt_msg(MODIS_E_END_VDATA_FILE, routine, "");
     returnStatus = FAIL;
    }


  /**************************************************************************/
  /*                                                                        */
  /*  return returnStatus                                                   */
  /*                                                                        */
  /**************************************************************************/

  return (returnStatus);

 }
