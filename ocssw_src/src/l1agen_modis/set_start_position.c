#include "PGS_IO.h"
#include "PGS_TD.h"
#include "PGS_SMF.h"
#include "PGS_TYPES.h"
#include "PGS_MODIS_35005.h"
#include "PC_pcf_info.h"
#include "L1A_prototype.h"


PGSt_SMF_status  set_start_position (PGSt_IO_L0_VirtualDataSet  L0_file,
                                     PGSt_double     *preload_start_time,
                                     PGSt_double     start_time,
                                     PGSt_double     stop_time)

/*
!C************************************************************************

!Description:  This function sets a pointer in the L0 file to start reading 
               data at a specific time (preload_start_time).


!Input Parameters:
               PGSt_IO_L0_VirtualDataSet  L0_file     
                                                ** L0 virtual data set     **

               PGSt_double  preload_start_time  ** Start time for          **
                                                ** preloading Eng data     **

!Output Parameters:
               None

Return Values: 
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_F_L0_SETSTART_FAILED         (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status                 (PGS_SMF.h)
               PGS_TRUE                        (PGS_SMF.h)
               PGS_FALSE                       (PGS_SMF.h)
               PGSt_SMF_boolean                (PGS_SMF.h)
               PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
               PGSd_EOS_AM1                    (PGS_TD.h)
               PGSt_Logical                    (PGS_TYPES.h)

Called By:     
               initialize_level1a

Routines Called:
               PGS_IO_L0_SetStart 
               PGS_IO_L0_Close
               get_valid_L0_file
               validate_L0_header
               log_fmt_msg  

!Revision History:
               Revision 2.2  2000/07/14
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Added code for Aqua instrument.

	       Revision 2.1  2000/04/06  11:54
               John Seaton/GSC (seaton@ltpmail.gsfc.nasa.gov)
	       Fix to change start/stop times when trying to process L0 file
               in LUN 599002. Call to get_valid_L0_file was using local vars
               to store start/stop times. These variables were never used in
               the comparison under case PGSIO_W_L0_TIME_NOT_FOUND. Thus the
               pre-load start time was not getting set and processing fails
               with return code 1.

               Revision 2.0  1997/08/28  02:23
               Timi Adelekan/GSC/SAIC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design

               Revision 1.0  1997/07/24  15:30
               Tom Johnson/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Original design

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes:
               The details for the packet data locations were taken from 
               Santa Barbara Remote Sensing (SBRS) Contract Data Requirement 
               List (CDRL)305, MODIS Engineering Telemetry Description, 
               Tables T30-5A, T30-5B, T30-5C, T30-5D, and T30-5E.


!END*************************************************************************
*/

 {
  /**************************************************************************/
  /*                                                                        */
  /*                      Define Global Variables                           */
  /*                                                                        */
  /**************************************************************************/

  extern PGSt_double global_first_gran_start_time;
                                          /* Start time for the first L1A   */
                                          /* output granule to be created   */

  extern  PGSt_PC_Logical  global_L0_logical;
                                          /*  L0 logical file unit number   */

  /**************************************************************************/
  /*              Define and Initialize Local Variables                     */
  /**************************************************************************/
  /*                                                                        */
  /* set returnStatus equal to MODIS_S_SUCCESS                              */
  /* Set continue_trying to PGS_TRUE                                        */
  /*                                                                        */
  /**************************************************************************/

  char             *routine = "set_start_position";
  PGSt_SMF_status  returnStatus;  /* SMF-style message returned by function */
  PGSt_SMF_status     PGSstatus;  /* SMF-style message returned by function */
  PGSt_SMF_status     L1Astatus;  /* SMF-style message returned by function */
  PGSt_SMF_boolean continue_trying; /* Boolean variable                     */


  L1Astatus       = MODIS_S_SUCCESS;
  PGSstatus       = MODIS_S_SUCCESS;
  returnStatus    = MODIS_S_SUCCESS;
  continue_trying = PGS_TRUE;


  /**************************************************************************/
  /*                                                                        */
  /* DO-WHILE continue_trying                                               */
  /*                                                                        */
  /*    CALL PGS_IO_L0_SetStart to set the virtual file pointer to the      */
  /*         first available packet after the preload_start_time            */
  /*      INPUTS:  L0_file, preload_start_time                              */
  /*      OUTPUTS: None                                                     */
  /*      RETURN:  PGSstatus                                                */
  /*                                                                        */
  /**************************************************************************/

     while (continue_trying) {
        PGSstatus = PGS_IO_L0_SetStart (L0_file, *preload_start_time);


  /**************************************************************************/
  /*                                                                        */
  /*    SWITCH (PGSstatus)                                                  */
  /*       CASE PGSIO_W_L0_TIME_NOT_FOUND                                   */
  /*          IF stop_time is greater than or equal to                      */
  /*             global_first_gran_start_time                               */
  /*          THEN                                                          */
  /*             Set preload_start_time to start_time                       */
  /*          ELSE                                                          */
  /*             IF global_L0_logical equals PRIOR_L0_PCF_ID                */
  /*             THEN                                                       */
  /*                CALL PGS_IO_L0_Close to close the prior L0 file         */
  /*                  INPUTS:  L0_file                                      */
  /*                  OUTPUTS: None                                         */
  /*                  RETURN:  PGSstatus                                    */
  /*                                                                        */
  /**************************************************************************/

        switch (PGSstatus) {
           case PGSIO_W_L0_TIME_NOT_FOUND:
              if (stop_time >= global_first_gran_start_time)
                 *preload_start_time = start_time;
              else {
                 if (global_L0_logical == PC_PRIOR_L0_PCF_ID) {
                    PGSstatus = PGS_IO_L0_Close (L0_file);


  /**************************************************************************/
  /*                                                                        */
  /*                CALL PGS_SMF_TestSuccessLevel to determine if the prior */
  /*                     L0 file was closed                                 */
  /*                  INPUTS:  PGSstatus                                    */
  /*                  OUTPUTS: None                                         */
  /*                  RETURN:  result                                       */
  /*                                                                        */
  /*                IF result is equal to PGS_FALSE                         */
  /*                THEN                                                    */
  /*                   CALL log_fmt_msg to report that the prior L0 file    */
  /*                      could not be closed                               */
  /*                     INPUTS:  PGSstatus, routine, msg                   */
  /*                     OUTPUTS: None                                      */
  /*                     RETURN:  None                                      */
  /*                ENDIF                                                   */
  /*                                                                        */
  /*                Set global_L0_logical to CURRENT_L0_PCF_ID              */
  /*                                                                        */
  /*         NOTE: Modis timestamps on both instruments, Terra and Aqua,    */
  /*          use segmented timestamps. Therefore, no change is needed for  */
  /*          Aqua delivery to read packet times. The Terra toolkit macros  */
  /*          PGSd_EOS_AM1 are used so L0 files are opened using the        */
  /*          toolkit routines which read segmented time.                   */
  /*                                                                        */
  /*                CALL get_valid_L0_file                                  */
  /*                    INPUTS: PGSd_EOS_AM1                                */
  /*                    OUTPUTS: L0_file, start_time, stop_time             */
  /*                    RETURNS: L1Astatus                                  */
  /*                ENDIF                                                   */
  /*                                                                        */
  /**************************************************************************/

                    if (!PGS_SMF_TestSuccessLevel(PGSstatus))
                       log_fmt_msg (MODIS_E_PGS_IO_LO_CLOSE, routine,
                                    "unable to close the prior L0 file");

                    global_L0_logical = PC_CURRENT_L0_PCF_ID;

                    L1Astatus = get_valid_L0_file ( PGSd_EOS_AM1, 
                                                    &L0_file,
                                                    &start_time,
                                                    &stop_time);

  /**************************************************************************/
  /*                                                                        */
  /*                IF L1Astatus is not equal to MODIS_S_SUCCESS            */
  /*                THEN                                                    */
  /*                   CALL log_fmt_msg to report that the current L0 file  */
  /*                      could not be opened                               */
  /*                     INPUTS:  L1Astatus, routine,                       */
  /*                               "unable to open the current L0 file"     */
  /*                     OUTPUTS: None                                      */
  /*                     RETURN:  None                                      */
  /*                                                                        */
  /*                   Set returnStatus to MODIS_F_L0_SETSTART_FAILED       */
  /*                   Set continue_trying to PGS_FALSE                     */
  /*                ENDIF                                                   */
  /*                                                                        */
  /**************************************************************************/

                    if (L1Astatus != MODIS_S_SUCCESS) {
                       log_fmt_msg (MODIS_E_GET_VALID_L0_FILE, routine,
                                    "unable to open the current L0 file");

                       returnStatus = MODIS_F_L0_SETSTART_FAILED;
                       continue_trying = PGS_FALSE;
                    }
                 }



  /**************************************************************************/
  /*                                                                        */
  /*             ELSE                                                       */
  /*                Set continue_trying to PGS_FALSE                        */
  /*                Set returnStatus to MODIS_F_L0_SETSTART_FAILED          */
  /*                CALL log_fmt_msg to report that the preload start time  */
  /*                   is outside of the current L0 file                    */
  /*                  INPUTS:  returnStatus, routine,                       */
  /*                        "preload start time outside of current L0 file" */
  /*                  OUTPUTS: None                                         */
  /*                  RETURN:  None                                         */
  /*             ENDIF                                                      */
  /*          ENDIF                                                         */
  /*                                                                        */
  /**************************************************************************/

                 else {
                    continue_trying = PGS_FALSE;
                    returnStatus = MODIS_F_L0_SETSTART_FAILED;
                    log_fmt_msg (MODIS_F_L0_SETSTART_FAILED, routine,
                               "preload start time outside of current L0 file");
                 }
              }

              break;

  /**************************************************************************/
  /* If call to PGS_IO_L0_SetStart returns PGSIO_E_L0_SEEK_PACKET, the      */
  /* start position could not be found in the L0 file. This may be due to   */
  /* missing data in the L0 file. When this value is returned and the       */
  /* current working L0 file is stored in LUN 599001, the process should    */
  /* try to open the 599002 LUN and see if its data is valid.               */
  /**************************************************************************/

          case PGSIO_E_L0_SEEK_PACKET:

  /**************************************************************************/
  /* If global_L0_logical == PC_PRIOR_L0_PCF_ID (599001) then               */
  /*   CALL PGS_IO_L0_Close to close the 599001 L0 file                     */
  /*     INPUTS: L0_file                                                    */
  /*     OUTPUTS: None                                                      */
  /*     RETURNS: PGSstatus                                                 */
  /**************************************************************************/

              if (global_L0_logical == PC_PRIOR_L0_PCF_ID) {
                 PGSstatus = PGS_IO_L0_Close (L0_file);

  /**************************************************************************/
  /* If the call to PGS_IO_L0_Close to close the 599001 L0 file fails, then */
  /*   CALL log_fmt_msg to report that the file could not be closed         */
  /*     INPUTS: MODIS_E_PGS_IO_LO_CLOSE, routine,                          */
  /*             "unable to close the prior L0 file"                        */
  /*     OUTPUTS: None                                                      */
  /*     RETURNS: None                                                      */
  /**************************************************************************/
  
                 if (!PGS_SMF_TestSuccessLevel(PGSstatus))
                    log_fmt_msg (MODIS_E_PGS_IO_LO_CLOSE, routine,
                                 "unable to close the prior L0 file");
  
  /**************************************************************************/
  /* set the global_L0_logical = PC_CURRENT_L0_PCF_ID (599002)              */
  /*                                                                        */
  /*         NOTE: Modis timestamps on both instruments, Terra and Aqua,    */
  /*          use segmented timestamps. Therefore, no change is needed for  */
  /*          Aqua delivery to read packet times. The Terra toolkit macros  */
  /*          PGSd_EOS_AM1 are used so L0 files are opened using the        */
  /*          toolkit routines which read segmented time.                   */
  /*                                                                        */
  /*            CALL get_valid_L0_file                                      */
  /*                INPUTS: PGSd_EOS_AM1                                    */
  /*                OUTPUTS: L0_file, start_time, stop_time                 */
  /*                RETURNS: L1Astatus                                      */
  /*            ENDIF                                                       */
  /**************************************************************************/

                 global_L0_logical = PC_CURRENT_L0_PCF_ID;
  
                 L1Astatus = get_valid_L0_file ( PGSd_EOS_AM1,
                                                 &L0_file,
                                                 &start_time,
                                                 &stop_time);

  /**************************************************************************/
  /* If the return from get_valid_L0_file is not MODIS_S_SUCCESS then       */
  /*   CALL log_fmt_msg to report that the L0 file could not be opened      */
  /*      INPUTS: MODIS_E_GET_VALID_L0_FILE, routine,                       */
  /*              "unable to open the current L0 file"                      */
  /*      OUTPUTS: None                                                     */
  /*      RETURNS: None                                                     */
  /*                                                                        */
  /*   set returnStatus = MODIS_F_L0_SETSTART_FAILED                        */
  /*   set continue_trying = PGS_FALSE                                      */
  /**************************************************************************/
                    
                 if (L1Astatus != MODIS_S_SUCCESS) {
                    log_fmt_msg (MODIS_E_GET_VALID_L0_FILE, routine,
                                 "unable to open the current L0 file");
                                                   
                    returnStatus = MODIS_F_L0_SETSTART_FAILED;
                    continue_trying = PGS_FALSE;
                 }
              }
  
  /***************************************************************************/
  /* The current L0 file is sotred in LUN 599002 already and processing      */
  /* should be halted.                                                       */
  /*  else                                                                   */
  /*     set continue_trying = PGS_FALSE                                     */
  /*     set returnStatus = MODIS_F_L0_SETSTART_FAILED                       */
  /*     CALL log_fmt_msg to report that the start time is outside of the    */
  /*         current L0 file                                                 */
  /*       INPUTS: MODIS_F_L0_SETSTART_FAILED, routine,                      */
  /*               "preload start time outside of current L0 file"           */
  /*       OUTPUTS: None                                                     */
  /*       RETURNS: None                                                     */
  /*  ENDIF                                                                  */
  /***************************************************************************/
              else {
                 continue_trying = PGS_FALSE;
                 returnStatus = MODIS_F_L0_SETSTART_FAILED;
                 log_fmt_msg (MODIS_F_L0_SETSTART_FAILED, routine,
                            "preload start time outside of current L0 file");
              }

          break;
                     




  /**************************************************************************/
  /*                                                                        */
  /*       CASE PGSIO_M_L0_HEADER_CHANGED                                   */
  /*          CALL validate_L0_header to validate the L0 header             */
  /*            INPUT:  L0_file                                             */
  /*            OUTPUT: None                                                */
  /*            RETURN: L1A_status                                          */
  /*                                                                        */
  /*          IF L1A_status is not equal to MODIS_S_SUCCESS                 */
  /*          THEN                                                          */
  /*             Set returnStatus to MODIS_F_L0_SETSTART_FAILED             */
  /*             CALL log_fmt_msg to report that the L0 file header could   */
  /*                  not be validated                                      */
  /*               INPUT:  L1A_status, routine, msg                         */
  /*               OUTPUT: None                                             */
  /*               RETURN: None                                             */
  /*          ENDIF                                                         */
  /*                                                                        */
  /*          Set continue_trying to PGS_FALSE                              */
  /*                                                                        */
  /**************************************************************************/

           case PGSIO_M_L0_HEADER_CHANGED:
             /* Comment out for now; use "inhouse" tools another time.
              L1Astatus = validate_L0_header (L0_file);
 
              if (L1Astatus != MODIS_S_SUCCESS) {
                 returnStatus = MODIS_F_L0_SETSTART_FAILED;
                 log_fmt_msg (MODIS_F_L0_HEADER_VAL_FAILED, routine, "");
              }

              continue_trying = PGS_FALSE;
             */
              break;
           
  

  /**************************************************************************/
  /*                                                                        */
  /*       CASE PGS_S_SUCCESS                                               */
  /*          Set continue_trying to PGS_FALSE                              */
  /*                                                                        */
  /*       DEFAULT                                                          */
  /*          Set continue_trying to PGS_FALSE                              */
  /*          Set returnStatus to MODIS_F_L0_SETSTART_FAILED                */
  /*    END_SWITCH                                                          */
  /* END_WHILE                                                              */
  /*                                                                        */
  /* RETURN returnStatus                                                    */
  /*                                                                        */
  /**************************************************************************/

           case PGS_S_SUCCESS:
              continue_trying = PGS_FALSE;
              break;

           default:
              continue_trying = PGS_FALSE;
              returnStatus = MODIS_F_L0_SETSTART_FAILED;
              break;
        }
     }

     return returnStatus;

 } /* End of routine set_start_position */
