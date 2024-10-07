#include "L1A_prototype.h"
#include "mapi.h"
#include "MD_metadata.h"
#include "PGS_MODIS_35005.h"
#include "mapiL1A.h"
#include "PGS_SMF.h"

PGSt_SMF_status   write_scan_metadata (MODFILE        *L1A_file_ptr, 
                                       MD_SCAN_MET_t  *scan_meta)

/*
!C*****************************************************************************

!Description:  This function will write out the Scan Level Metadata (section
               2 of the MODIS Level 1A Data Product Format) to the L1A file.

!Input Parameters:
               MODFILE    L1A_file_ptr     ** A pointer to the L1A file   **
               MD_SCAN_MET_t  scan_meta    ** Scan level metadata for 
                                              this scan                   **

!Output Parameters: 
               None

Return Values: 
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_ARRAY_OUTPUT_ERR     (PGS_MODIS_35005.h) 

Externally Defined:
               PGSt_SMF_status              (PGS_SMF.h)
               MD_SCAN_MET_t                (MD_metadata.h)
               MODFILE                      (mapi.h)
               M01SCAN_NUMBER               (mapiL1A.h)
               M01FRAME_COUNT_ARRAY         (mapiL1A.h)
               M01SCAN_TYPE                 (mapiL1A.h)
               M01SD_START_TIME             (mapiL1A.h)
               M01SRCA_START_TIME           (mapiL1A.h)
               M01BB_START_TIME             (mapiL1A.h)
               M01SV_START_TIME             (mapiL1A.h)
               M01EV_START_TIME             (mapiL1A.h)
               M01SRCA_CALIBRATION_MODE     (mapiL1A.h)
               M01PACKET_SCAN_COUNT         (mapiL1A.h)
               M01CCSDS_APID                (mapiL1A.h)
               M01PACKET_QL                 (mapiL1A.h)
               M01MIRROR_SIDE               (mapiL1A.h)
               M01SCAN_QUALITY_ARRAY        (mapiL1A.h)

Called By:
               write_scan
               handle_missing_scans

Routines Called:
               putMODISarray
               log_fmt_msg

!Revision History:
               revision 1.0 1997/09/04  17:30:00
               Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
               Original development

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               None

!END**********************************************************************
*/
{
  char                         *routine = "write_scan_metadata";
  char                         msg[300];
  PGSt_SMF_status              returnStatus;
  long int                     start[2];
  long int                     dimsize[2] = {1,0};
  long int                     dimsize_1[1] = {1};

  start[0] = scan_meta->scan_num -1;
  start[1] = 0;

  /******************************************************************************/
  /*                                                                            */
  /*  set returnStatus to MODIS_S_SUCCESS                                       */
  /*                                                                            */
  /******************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Scan number                                 */
  /*  Compute the start for the Scan number                                     */
  /*  CALL putMODISarray to write the Scan number to the L1A granule            */
  /*    INPUTS: MODFILE, M01SCAN_NUMBER, NULL, start, dimsize array,            */
  /*            MD_SCAN_MET_t.scan_num                                          */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set msg to "The Scan number could not be written to the L1A granule"   */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Scan number could not be written to*/
  /*        the L1A granule                                                     */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SCAN_NUMBER ,NULL, start, dimsize_1,
                    &scan_meta->scan_num) == MFAIL)
    {
      sprintf(msg,"The Scan number could not be written to the L1A granule %s ",
              L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Frame count array                           */
  /*  Compute the start for the Frame count array                               */
  /*  CALL putMODISarray to write the Frame count array to the L1A granule      */
  /*    INPUTS: MODFILE, M01FRAME_COUNT_ARRAY, NULL, start, dimsize array,      */
  /*            MD_SCAN_MET_t.frame_count_array                                 */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Frame count array could not be written to the L1A      */
  /*        granule"                                                            */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Frame count array could not be     */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/
                                                      
  dimsize[1] = 6;

  if (putMODISarray(L1A_file_ptr, M01FRAME_COUNT_ARRAY, NULL, start, dimsize,
                    scan_meta->frame_count_array) == MFAIL)
    {
      sprintf(msg,"The Frame count array could not be written to the L1A granule: %s ",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Scan Type                                   */
  /*  Compute the start for the Scan Type                                       */
  /*  CALL putMODISarray to write the Scan Type to the L1A granule              */
  /*    INPUTS: MODFILE, M01SCAN_TYPE, NULL, start, dimsize array,              */
  /*            MD_SCAN_MET_t.scan_type                                         */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Scan Type could not be written to the L1A granule"     */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Scan Type could not be written to  */
  /*        the L1A granule                                                     */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/
                                                      
  dimsize[1] = 10;

  if (putMODISarray(L1A_file_ptr, M01SCAN_TYPE, NULL, start, dimsize,
                    scan_meta->scan_type) == MFAIL)
    {
      sprintf(msg,"The Scan Type could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the SD start time                               */
  /*  Compute the start for the SD start time                                   */
  /*  CALL putMODISarray to write the SD start time to the L1A granule          */
  /*    INPUTS: MODFILE, M01SD_START_TIME, NULL, start, dimsize array,          */
  /*            MD_SCAN_MET_t.sd_start_time                                     */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The SD start time could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the SD start time could not be         */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SD_START_TIME, NULL, start, dimsize_1,
                    &scan_meta->sd_start_time) == MFAIL)
    {
      sprintf(msg,"The SD start time could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the SRCA start time                             */
  /*  Compute the start for the SRCA start time                                 */
  /*  CALL putMODISarray to write the SRCA start time to the L1A granule        */
  /*    INPUTS: MODFILE, M01SRCA_START_TIME, NULL, start, dimsize array,        */
  /*            MD_SCAN_MET_t.srca_start_time                                   */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The SRCA start time could not be written to the L1A        */
  /*        granule"                                                            */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the SRCA start time could not be       */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SRCA_START_TIME, NULL, start, dimsize_1,
                    &scan_meta->srca_start_time) == MFAIL)
    {
      sprintf(msg,"The SRCA start time could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the BB start time                               */
  /*  Compute the start for the BB start time                                   */
  /*  CALL putMODISarray to write the BB start time to the L1A granule          */
  /*    INPUTS: MODFILE, M01BB_START_TIME, NULL, start, dimsize array,          */
  /*            MD_SCAN_MET_t.bb_start_time                                     */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The BB start time could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the BB start time could not be         */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01BB_START_TIME, NULL, start, dimsize_1,
                    &scan_meta->bb_start_time) == MFAIL)
    {
      sprintf(msg,"The BB start time could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the SV start time                               */
  /*  Compute the start for the SV start time                                   */
  /*  CALL putMODISarray to write the SV start time to the L1A granule          */
  /*    INPUTS: MODFILE, M01SV_START_TIME, NULL, start, dimsize array,          */
  /*            MD_SCAN_MET_t.sv_start_time                                     */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The SV start time could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the SV start time could not be         */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SV_START_TIME, NULL, start, dimsize_1,
                    &scan_meta->sv_start_time) == MFAIL)
    {
      sprintf(msg,"The SV start time could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the EV start time                               */
  /*  Compute the start for the EV start time                                   */
  /*  CALL putMODISarray to write the EV start time to the L1A granule          */
  /*    INPUTS: MODFILE, M01EV_START_TIME, NULL, start, dimsize array,          */
  /*            MD_SCAN_MET_t.ev_start_time                                     */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The EV start time could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the EV start time could not be         */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01EV_START_TIME, NULL, start, dimsize_1,
                    &scan_meta->ev_start_time) == MFAIL)
    {
      sprintf(msg,"The EV start time could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the SRCA calibration mode                       */
  /*  Compute the start for the SRCA calibration mode                           */
  /*  CALL putMODISarray to write the SRCA calibration mode to the L1A granule  */
  /*    INPUTS: MODFILE, M01SRCA_CALIBRATION_MODE, NULL, start, dimsize array,  */
  /*            MD_SCAN_MET_t.srca_cal_mode                                     */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The SRCA calibration mode could not be written to the L1A  */
  /*        granule"                                                            */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the SRCA calibration mode could not be */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SRCA_CALIBRATION_MODE, NULL, start, dimsize_1,
                    &scan_meta->srca_cal_mode) == MFAIL)
    {
      sprintf(msg,"The SRCA calibration mode could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Packet scan count                           */
  /*  Compute the start for the Packet scan count                               */
  /*  CALL putMODISarray to write the Packet scan count to the L1A granule      */
  /*    INPUTS: MODFILE, M01PACKET_SCAN_COUNT, NULL, start, dimsize array,      */
  /*            MD_SCAN_MET_t.packet_scan_count                                 */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Packet scan count could not be written to the L1A      */
  /*        granule"                                                            */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Packet scan count could not be     */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01PACKET_SCAN_COUNT, NULL, start, dimsize_1,
                    &scan_meta->packet_scan_count) == MFAIL)
    {
      sprintf(msg,"The Packet scan count could not be written to the L1A granule %s",
              L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the CCSDS Application Identifiers               */
  /*  Compute the start for the CCSDS Application Identifiers                   */
  /*  CALL putMODISarray to write the CCSDS Application Identifiers to the      */
  /*     L1A granule                                                            */
  /*    INPUTS: MODFILE, M01CCSDS_APID, NULL, start, dimsize array,             */
  /*            MD_SCAN_MET_t.ccsds_apids                                       */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The CCSDS Application Identifiers could not be written to  */
  /*        the L1A granule"                                                    */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the CCSDS Application Identifiers could*/
  /*        not be written to the L1A granule                                   */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/
                                                      
  dimsize[1] = 3;

  if (putMODISarray(L1A_file_ptr, M01CCSDS_APID, NULL, start, dimsize,
                    scan_meta->ccsds_apids) == MFAIL)
    {
      sprintf(msg,"The CCSDS Application Identifiers could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Packet expedited data flag                  */
  /*  Compute the start for the Packet expedited data flag                      */
  /*  CALL putMODISarray to write the Packet expedited data flag to the L1A     */
  /*     granule                                                                */
  /*    INPUTS: MODFILE, M01PACKET_QL, NULL, start, dimsize array,              */
  /*            MD_SCAN_MET_t.packet_expedited_data_flag                        */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Packet expedited data flag could not be written to the */
  /*        L1A granule"                                                        */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Packet expedited data flag could   */
  /*        not be written to the L1A granule                                   */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01PACKET_QL, NULL, start, dimsize_1,
                    &scan_meta->packet_expedited_data_flag) == MFAIL)
    {
      sprintf(msg,"The Packet expedited data flag could not be written to the L1A granule %s",
              L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Mirror side                                 */
  /*  Compute the start for the Mirror side                                     */
  /*  CALL putMODISarray to write the Mirror side to the L1A granule            */
  /*    INPUTS: MODFILE, M01MIRROR_SIDE, NULL, start, dimsize array,            */
  /*            MD_SCAN_MET_t.mirror_side                                       */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Mirror side could not be written to the L1A granule"   */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Mirror side could not be           */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01MIRROR_SIDE, NULL, start, dimsize_1,
                    &scan_meta->mirror_side) == MFAIL)
    {
      sprintf(msg,"The Mirror side could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
                                                      

  /******************************************************************************/
  /*                                                                            */
  /*  Compute dimsize array for the Scan quality array                          */
  /*  Compute the start for the Scan quality array                              */
  /*  CALL putMODISarray to write the Scan quality array to the L1A granule     */
  /*    INPUTS: MODFILE, M01SCAN_QUALITY_ARRAY, NULL, start, dimsize array,     */
  /*            MD_SCAN_MET_t.scan_qual_array                                   */
  /*    OUTPUT: None                                                            */
  /*    RETURN: mapiStatus                                                      */
  /*                                                                            */
  /*  IF mapiStatus equals MFAIL                                                */
  /*  THEN                                                                      */
  /*     set routine to "write_scan_metadata"                                   */
  /*     set msg to "The Scan quality array could not be written to the L1A     */
  /*        granule"                                                            */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                           */
  /*     CALL log_fmt_msg to report that the Scan quality array could not be    */
  /*        written to the L1A granule                                          */
  /*       INPUTS: returnStatus, routine, msg                                   */
  /*       OUTPUT: None                                                         */
  /*       RETURN: None                                                         */
  /*  ENDIF                                                                     */
  /*                                                                            */
  /******************************************************************************/

  dimsize[1] = 4;

  if (putMODISarray(L1A_file_ptr, M01SCAN_QUALITY_ARRAY, NULL, start, dimsize,
                    scan_meta->scan_qual_array) == MFAIL)
    {
      sprintf(msg,"The Scan quality array could not be written to the L1A granule %s",
              L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /******************************************************************************/
  /*                                                                            */
  /*  RETURN returnStatus                                                       */
  /*                                                                            */
  /******************************************************************************/
                                                      
  return (returnStatus);
                                                      
}
