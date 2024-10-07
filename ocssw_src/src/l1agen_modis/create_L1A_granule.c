#include "L1A_prototype.h"
#include "hdf.h"
#include "PGS_MODIS_35005.h"
#include "EN_eng_data.h"
#include "PC_pcf_info.h"
#include "mapi.h"
#include "PGS_PC.h"
#include "PGS_SMF.h"

PGSt_SMF_status create_L1A_granule (EN_VDATA_TYPE_t *eng_data, 
                                    int16 nscans,
                                    MODFILE **L1A_file_ptr)
/*
!C*****************************************************************************

!Description:   Function create_L1A_granule opens an L1A file for output, and
                initializes all SDSs and Vdata structures for writing.
                 
!Input Parameters:
                EN_VDATA_TYPE_t  eng_data        **  Engineering data              **
                int16            nscans          **  Number of scans in a granule  **

!Output Parameters:
                MODFILE          **L1A_file_ptr  **  Address of L1A file pointer   **
 
Returns Values:                
                MODIS_S_SUCCESS                 (PGS_MODIS_35005.h)
                MODIS_E_CREATE_L1A_GRANULE      (PGS_MODIS_35005.h)

Externally Defined:
                MODFILE                         (mapi.h)
                PC_L1A_PROD_PCF_ID              (PC_pcf_info.h)

Called By:
                process_a_granule

Routines Called:
                PGS_SMF_TestSuccessLevel        (PGS_SMF.h)
                log_fmt_msg
                PGS_PC_GetReference             (PGS_PC.h)
                openMODISfile                   (mapi.h)
                init_L1A_HDF_vdatas
                init_L1A_HDF_sdss

!Revision History:
                revision 1.0 1997/08/25  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support Team 
                Team (SDST) for the National Aeronautics and Space Administration 
                (NASA), Goddard Space Flight Center (GSFC), under contract 
                NAS5-32373.
        
!References and Credits:
                None

!Design Notes:
                None

!END*************************************************************************
*/
{
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/

  char                 *routine = "create_L1A_granule";
  char                 msg[300];
  PGSt_SMF_status      returnStatus;
  PGSt_SMF_status      PGS_status;
  PGSt_integer         L1A_version_num;
  static PGSt_integer  version = 0;
  char                 L1A_filename[256];

  /*************************************************************************/
  /*                                                                       */
  /*  Set routine to "create_L1A_granule"                                  */
  /*                                                                       */
  /*  Set returnStatus to MODIS_S_SUCCESS                                  */
  /*                                                                       */
  /*************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /*************************************************************************/
  /*                                                                       */
  /*  CALL PGS_PC_GetReference to retrieve the L1A file name               */
  /*    INPUT:  L1A_PROD_PCF_ID, L1A_version_num                           */
  /*    OUTPUT: L1A_filename, L1A_version_num                              */
  /*    RETURN: PGS_status                                                 */
  /*                                                                       */
  /*************************************************************************/

  version ++;
  L1A_version_num = version;
  PGS_status = PGS_PC_GetReference(PC_L1A_PROD_PCF_ID, &L1A_version_num,
                                   L1A_filename);


  /***************************************************************************/
  /*                                                                         */
  /*  CALL PGS_SMF_TestSuccessLevel to determine if the retrieve is          */
  /*    successful                                                           */
  /*    INPUT:  PGS_status                                                   */
  /*    ONTPUT: None                                                         */
  /*    RETURN: TestStatus                                                   */
  /*                                                                         */
  /*  IF TestStatus is not equal to PGS_TRUE                                 */
  /*  THEN                                                                   */
  /*     Set returnStatus to MODIS_E_CREATE_L1A_GRANULE                      */
  /*     Set msg to "The L1A_filename could not be retrieved from the PCF    */
  /*       file"                                                             */
  /*     CALL log_fmt_msg to report that the L1A_filename was not retrieved  */
  /*       INPUT:  PGS_status, routine, msg                                  */
  /*       OUTPUT: None                                                      */
  /*       RETURN: None                                                      */
  /*  ELSE                                                                   */
  /*                                                                         */
  /*     Set access_mode to "w" (write only)                                 */
  /*     CALL openMODISfile to open the L1A file                             */
  /*       INPUT:  filename, access_mode                                     */
  /*       OUPUT:  None                                                      */
  /*       RETURN: L1A_file_ptr                                              */
  /*                                                                         */
  /*     IF L1A_file_ptr is equal to NULL                                    */
  /*     THEN                                                                */
  /*        Set returnStatus to MODIS_E_CREATE_L1A_GRANULE                   */
  /*        Set msg to "The L1A granule could not be opened for granule      */
  /*           (L1A_filename)"                                               */
  /*        CALL log_fmt_msg to report that the opening L1A file failed      */
  /*          INPUT:  routine, msg, filename                                 */
  /*          OUTPUT: None                                                   */
  /*          RETURN: None                                                   */
  /*     ELSE                                                                */
  /*                                                                         */
  /*        CALL init_L1A_HDF_vdatas to initialize the HDF vdatas in the     */
  /*          L1A file                                                       */
  /*          INPUT:  L1A_file_ptr, eng_data                                 */
  /*          OUTPUT: None                                                   */
  /*          RETURN: PGS_status                                             */
  /*                                                                         */
  /*        IF L1A_status is not equal to MODIS_S_SUCCESS                    */
  /*        THEN                                                             */
  /*           Set returnStatus to MODIS_E_CREATE_L1A_GRANULE                */
  /*           Set msg to "The L1A HDF vdatas could not be created for       */
  /*             granule (L1A_filename)"                                     */
  /*           CALL log_fmt_msg to report that the initializing L1A HDF      */
  /*              vdatas failed                                              */
  /*             INPUT:  L1A_status, routine, msg                            */
  /*             OUTPUT: None                                                */
  /*             RETURN: None                                                */
  /*                                                                         */
  /*           CALL completeMODISfile to close the L1A granule               */
  /*             INPUT:  L1A_file_ptr, NULL, NULL, 0                         */
  /*             OUTPUT: None                                                */
  /*             RETURN: mapi_status                                         */
  /*                                                                         */
  /*           IF mapi_status is MFAIL                                       */
  /*           THEN                                                          */
  /*              Set returnStatus to MODIS_E_CREATE_L1A_GRANULE             */
  /*              Set msg to "The L1A granule could not be closed properly"  */
  /*              CALL log_fmt_msg to report that closing the L1A file       */
  /*                 was not successful                                      */
  /*                INPUT:  returnStatus, routine, msg                       */
  /*                OUTPUT: None                                             */
  /*                RETURN: None                                             */
  /*           ENDIF                                                         */
  /*        ELSE                                                             */
  /*                                                                         */
  /*           CALL init_L1A_HDF_sdss to initialize the HDF sdss in the      */
  /*             L1A file                                                    */
  /*             INPUT:  L1A_file_ptr, nscans                                */
  /*             OUTPUT: None                                                */
  /*             RETURN: PGS_status                                          */
  /*                                                                         */
  /*           IF L1A_status is not equal to MODIS_S_SUCCESS                 */
  /*           THEN                                                          */
  /*              Set returnStatus to MODIS_E_CREATE_L1A_GRANULE             */
  /*              Set msg to "The L1A HDF sdss could not be created for      */
  /*                granule (L1A_filename)"                                  */
  /*              CALL log_fmt_msg to report that the initializing L1A       */
  /*                HDF sdss failed                                          */
  /*                INPUT:  L1A_status, routine, msg                         */
  /*                OUTPUT: None                                             */
  /*                RETURN: None                                             */
  /*                                                                         */
  /*              CALL completeMODISfile to close the L1A granule            */
  /*                INPUT:  L1A_file_ptr, NULL, NULL, 0                      */
  /*                OUTPUT: None                                             */
  /*                RETURN: mapi_status                                      */
  /*                                                                         */
  /*              IF mapi_status is MFAIL                                    */
  /*              THEN                                                       */
  /*                 Set returnStatus to MODIS_E_CREATE_L1A_GRANULE          */
  /*                 Set msg to "The L1A granule could not be closed         */
  /*                    properly"                                            */
  /*                 CALL log_fmt_msg to report that closing the L1A file    */
  /*                    was not successful                                   */
  /*                   INPUT:  returnStatus, routine, msg                    */
  /*                   OUTPUT: None                                          */
  /*                   RETURN: None                                          */
  /*              ENDIF                                                      */
  /*           ENDIF                                                         */
  /*        ENDIF                                                            */
  /*     ENDIF                                                               */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

  if (PGS_SMF_TestSuccessLevel(PGS_status) != PGS_TRUE)
    {
      returnStatus = MODIS_E_CREATE_L1A_GRANULE;
      log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, 
                   "The L1A filename could not be retrieved from the PCF file");
    }

  else
    {
      *L1A_file_ptr = openMODISfile(L1A_filename,"w"); 
      if (*L1A_file_ptr == NULL)
        {
          returnStatus = MODIS_E_CREATE_L1A_GRANULE;
          sprintf(msg,"The L1A filename %s could not be opened", L1A_filename);
          log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, msg);
        }
      else
        {
         if ((PGS_status = init_L1A_HDF_vdatas(eng_data, *L1A_file_ptr))
               != MODIS_S_SUCCESS)
           {
            returnStatus = MODIS_E_CREATE_L1A_GRANULE;
            sprintf(msg,"The L1A HDF vdatas could not be created for file: %s",
                   L1A_filename);
            log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, msg);
            if (completeMODISfile(L1A_file_ptr, NULL, NULL, 
                               0) == MFAIL)
              {
               returnStatus = MODIS_E_CREATE_L1A_GRANULE;
               sprintf(msg,"The L1A file: %s could not be closed properly", 
                      L1A_filename);
               log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, msg);
              }
           }
         else
           {
            if ((PGS_status =  init_L1A_HDF_sdss(*L1A_file_ptr, nscans))
                 != MODIS_S_SUCCESS)
              {
               returnStatus = MODIS_E_CREATE_L1A_GRANULE;
               sprintf(msg,"The L1A HDF sdss could not be created for file: %s",
                      L1A_filename);
               log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, msg);
               if (completeMODISfile(L1A_file_ptr, NULL, NULL, 
                               0) == MFAIL)
                 {
                  returnStatus = MODIS_E_CREATE_L1A_GRANULE;
                  sprintf(msg,"The L1A file: %s could not be closed properly",
                         L1A_filename);
                  log_fmt_msg(MODIS_E_CREATE_L1A_GRANULE, routine, msg);
                 }
              }
           } 
        }
     }
                
  /************************/
  /*  RETURN returnStatus */
  /************************/
  return (returnStatus);
}
