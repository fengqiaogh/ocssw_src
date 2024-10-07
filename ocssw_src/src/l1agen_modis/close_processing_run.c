#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_IO_L0.h"
#include "PGS_MODIS_35005.h"

void close_processing_run (PGSt_IO_L0_VirtualDataSet L0_file)

/*
!C***************************************************************************

!Description:   Function close_processing_run handles the wrap-up processing
                for each run of the MODIS Level 1A software.
                 
!Input Parameters:
                L0_file         **  L0 file descriptor  **

!Output Parameters:
                None

Return Values:                
                None

Externally Defined:
                PGSt_SMF_status                         (PGS_SMF.h)
                PGSt_IO_L0_VirtualDataSet               (PGS_IO_L0.h)

Called By:
                level1a

Routines Called:
                PGS_SMF_TestSuccessLevel                (PGS_SMF.h)
                PGS_IO_L0_Close                         (PGS_IO_L0_Wrap.h)
                log_fmt_msg

!Revision History:
                revision 2.0 1997/09/22  14:30:00
                Tom Johnson/GSC   (johnson@ltpmail.gsfc.nasa.gov)
                Original code

                revision 1.0 1997/08/05  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original design

!Team-unique Header:
                This software is developed by the MODIS Science Data Support 
                Team (SDST) for the National Aeronautics and Space 
                Administration (NASA), Goddard Space Flight Center (GSFC), 
                under contract NAS5-32373.

!References and Credits:

!Design Notes:


!END***************************************************************************
*/

{
  /***************************************************************************/
  /*              Declare and Initialize Local Variables                     */
  /***************************************************************************/

  extern PGSt_PC_Logical  global_L0_logical;
  PGSt_SMF_status  PGS_status;   /* SMF-style message returned by function */
  char             *routine = "close_processing_run";
  char             msg[300];

  /***************************************************************************/


  /***************************************************************************/
  /*                                                                         */
  /*  set routine to "close_processing_run" (done during declaration)        */
  /*                                                                         */
  /*                                                                         */
  /*  CALL PGS_IO_L0_Close to close the L0 file                              */
  /*    INPUT:  L0_file                                                      */
  /*    OUTPUT: None                                                         */
  /*    RETURN: PGS_status                                                   */
  /*                                                                         */
  /*  CALL PGS_SMF_TestSuccessLevel to determine if the closing file was     */
  /*     successful                                                          */
  /*    INPUT:  PGS_status                                                   */
  /*    OUPUT:  None                                                         */
  /*    RETURN: TestStatus                                                   */
  /*                                                                         */
  /*  IF TestStatus is not equal to PGS_TRUE                                 */
  /*  THEN                                                                   */
  /*    Set msg to "function close_processing_run was unable to close the    */
  /*       L0 file successfully"                                             */
  /*    CALL log_fmt_msg to report that the closing L0 file was not          */
  /*       successful                                                        */
  /*      INPUT:  PGS_status, routine, msg                                   */
  /*      OUTPUT: None                                                       */
  /*      RETURN: None                                                       */
  /*  ENDIF                                                                  */
  /*                                                                         */
  /***************************************************************************/

   PGS_status = PGS_IO_L0_Close(L0_file);
   if (!PGS_SMF_TestSuccessLevel(PGS_status)) 
     {
      sprintf (msg, "Unable to close the L0 file successfully LUN: %d",
                    global_L0_logical);
      log_fmt_msg(MODIS_E_PGS_IO_LO_CLOSE, routine, msg);
     }
}

