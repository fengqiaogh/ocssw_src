#include "L1A_prototype.h"
#include "mapi.h"
#include "PGS_MODIS_35005.h"
#include "SC_scan.h"
#include "MD_metadata.h"
#include "PD_pkt_data.h"
#include "PH_pkt_hdr.h"
#include "mapiL1A.h"
#include "PGS_SMF.h"

PGSt_SMF_status  write_scan_data (MODFILE         *L1A_file_ptr, 
                                  SC_SCAN_DATA_t  *L1A_scan,
                                  int16           scan_num,
                                  char            *scan_type)

/*
!C**********************************************************************

!Description:  This function will write out the Scan Data (section 4 of the
               MODIS Level 1A Data Product Format) to the L1A file.

!Input Parameters:
               MODFILE  *L1A_file_ptr       ** A pointer to the L1A file   **
               SC_SCAN_DATA_t  *L1A_scan    ** L1A scan data for this scan **
               int16           scan_num     ** scan number                 **
               char            *scan_type   ** scan type                   **

!Output Parameters: 
               None

Return Values: MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_ARRAY_OUTPUT_ERR     (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status                       (PGS_SMF.h)
               SC_SCAN_DATA_t                        (SC_scan.h)
               MODFILE                               (mapi.h)
               M01SD_250M                            (mapiL1A.h)
               M01SD_500M                            (mapiL1A.h)
               M01SD_1M_DAY                          (mapiL1A.h)
               M01SD_1M_NITE                         (mapiL1A.h)
               M01SRCA_250M                          (mapiL1A.h)
               M01SRCA_500M                          (mapiL1A.h)
               M01SRCA_1M_DAY                        (mapiL1A.h)
               M01SRCA_1M_NITE                       (mapiL1A.h)
               M01BB_250M                            (mapiL1A.h)
               M01BB_500M                            (mapiL1A.h)
               M01BB_1M_DAY                          (mapiL1A.h)
               M01BB_1M_NITE                         (mapiL1A.h)
               M01SV_250M                            (mapiL1A.h)
               M01SV_500M                            (mapiL1A.h)
               M01SV_1M_DAY                          (mapiL1A.h)
               M01SV_1M_NITE                         (mapiL1A.h)
               M01EV_250M                            (mapiL1A.h)
               M01EV_500M                            (mapiL1A.h)
               M01EV_1M_DAY                          (mapiL1A.h)
               M01EV_1M_NITE                         (mapiL1A.h)
               M01FPA_AEM_CONFIG                     (mapiL1A.h)
               M01SCIENCE_STATE                      (mapiL1A.h)
               M01SCIENCE_ABNORM                     (mapiL1A.h)
               M01FPA_DCR_OFFSET                     (mapiL1A.h)
               M01RAW_MIR_ENC                        (mapiL1A.h)
               M01RAW_VS_DEF                         (mapiL1A.h)
               M01RAW_VS_ACT                         (mapiL1A.h)
               M01RAW_SCI_ENG                        (mapiL1A.h)
               M01RAW_HK_TELEM                       (mapiL1A.h)
               M01RAW_SC_ANCIL                       (mapiL1A.h)
               M01RAW_PARAM                          (mapiL1A.h)
               M01RAW_PV_GAINS                       (mapiL1A.h)
               PD_DN_NUM_250M_DETECTORS              (PD_pkt_data.h)
               PD_DN_NUM_250M_BANDS                  (PD_pkt_data.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX    (PH_pkt_hdr.h)
               PD_DN_BAND_RATIO_250M                 (PD_pkt_data.h)
               PD_DN_NUM_500M_DETECTORS              (PD_pkt_data.h)
               PD_DN_NUM_500M_BANDS                  (PD_pkt_data.h)   
               PD_DN_BAND_RATIO_500M                 (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_DETECTORS            (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_BANDS                (PD_pkt_data.h)
               PD_DN_BAND_RATIO_1KM                  (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_DETECTORS          (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_BANDS              (PD_pkt_data.h)
               PD_DN_BAND_RATIO_1KM                  (PD_pkt_data.h)

Called By:
               write_scan

Routines Called:
               putMODISarray
               log_fmt_msg

!Revision History:
               revision 1.0 1997/09/03  17:30:00
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

!END***************************************************************************
*/
{
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/

  char               *routine = "write_scan_data";
  char               msg[300];
  PGSt_SMF_status    returnStatus;
  long int           start[3] = {0, 0, 0};
  long int           start_250[3] = {0, 0, 0};
  long int           start_500[3] = {0, 0, 0};
  long int           start_day[3] = {0, 0, 0};
  long int           start_night[3] = {0, 0, 0};
  long int           dimsize[3];
  long int           dimsize_cal_250[3] = {PD_DN_NUM_250M_DETECTORS,
                                           PD_DN_NUM_250M_BANDS,
                                           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX *
                                           PD_DN_BAND_RATIO_250M};
  long int           dimsize_cal_500[3] = {PD_DN_NUM_500M_DETECTORS,
                                           PD_DN_NUM_500M_BANDS,
                                           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX *
                                           PD_DN_BAND_RATIO_500M};
  long int           dimsize_cal_day[3] = {PD_DN_NUM_1KMDAY_DETECTORS,
                                           PD_DN_NUM_1KMDAY_BANDS,
                                           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX *
                                           PD_DN_BAND_RATIO_1KM};
  long int           dimsize_cal_night[3] = {PD_DN_NUM_1KMNIGHT_DETECTORS,
                                             PD_DN_NUM_1KMNIGHT_BANDS,
                                             PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX *
                                             PD_DN_BAND_RATIO_1KM};


  start[0] = scan_num - 1;
  start_250[0] = (scan_num - 1) * PD_DN_NUM_250M_DETECTORS;
  start_500[0] = (scan_num - 1) * PD_DN_NUM_500M_DETECTORS;
  start_day[0] = (scan_num - 1) * PD_DN_NUM_1KMDAY_DETECTORS;
  start_night[0] = (scan_num - 1) * PD_DN_NUM_1KMNIGHT_DETECTORS;

  /*****************************************************************************/
  /*                                                                           */
  /*  set returnStatus to MODIS_S_SUCCESS                                      */
  /*                                                                           */
  /*****************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SD 250m data                               */
  /*  Compute the start for the SD 250m data                                   */
  /*  CALL putMODISarray to write the SD 250m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01SD_250M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.SD_250m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SD 250m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SD 250m data could not be         */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SD_250M,NULL,start_250,dimsize_cal_250,
                    L1A_scan->SD_250m) == MFAIL)
    {
      sprintf(msg,"The SD 250m data could not be written to the L1A granule %s",
            L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SD 500m data                               */
  /*  Compute the start for the SD 500m data                                   */
  /*  CALL putMODISarray to write the SD 500m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01SD_500M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.SD_500m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SD 500m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SD 500m data could not be         */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SD_500M,NULL,start_500,dimsize_cal_500,
                    L1A_scan->SD_500m) == MFAIL)
    {
      sprintf(msg,"The SD 500m data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SD 1km day data                            */
  /*  Compute the start for the SD 1km day data                                */
  /*  CALL putMODISarray to write the SD 1km day data to the L1A granule       */
  /*    INPUTS: MODFILE, M01SD_1M_DAY, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.SD_1km_day                                      */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SD 1km day data could not be written to the L1A       */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SD 1km day data could not be      */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SD_1KM_DAY,NULL,start_day,dimsize_cal_day,
                    L1A_scan->SD_1km_day) == MFAIL)
    {
      sprintf(msg,"The SD 1km day data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SD 1km night data                          */
  /*  Compute the start for the SD 1km night data                              */
  /*  CALL putMODISarray to write the SD 1km night data to the L1A granule     */
  /*    INPUTS: MODFILE, M01SD_1M_NITE, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.SD_1km_night                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SD 1km night data could not be written to the L1A     */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SD 1km night data could not be    */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SD_1KM_NITE,NULL,start_night,dimsize_cal_night,
                    L1A_scan->SD_1km_night) == MFAIL)
    {
      sprintf(msg,"The SD 1km_night data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SRCA 250m data                             */
  /*  Compute the start for the SRCA 250m data                                 */
  /*  CALL putMODISarray to write the SRCA 250m data to the L1A granule        */
  /*    INPUTS: MODFILE, M01SRCA_250M, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.SRCA_250m                                       */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SRCA 250m data could not be written to the L1A        */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SRCA 250m data could not be       */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SRCA_250M,NULL,start_250,dimsize_cal_250,
                    L1A_scan->SRCA_250m) == MFAIL)
    {
      sprintf(msg,"The SRCA 250m data could not be written to the L1A granule %s",
            L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SRCA 500m data                             */
  /*  Compute the start for the SRCA 500m data                                 */
  /*  CALL putMODISarray to write the SRCA 500m data to the L1A granule        */
  /*    INPUTS: MODFILE, M01SRCA_500M, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.SRCA_500m                                       */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SRCA 500m data could not be written to the L1A        */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SRCA 500m data could not be       */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SRCA_500M,NULL,start_500,dimsize_cal_500,
                    L1A_scan->SRCA_500m) == MFAIL)
    {
      sprintf(msg,"The SRCA 500m data could not be written to the L1A granule %s",
            L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SRCA 1km day data                          */
  /*  Compute the start for the SRCA 1km day data                              */
  /*  CALL putMODISarray to write the SRCA 1km day data to the L1A granule     */
  /*    INPUTS: MODFILE, M01SRCA_1M_DAY, NULL, start, dimsize array,           */
  /*            SC_SCAN_DATA_t.SRCA_1km_day                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SRCA 1km day data could not be written to the L1A     */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SRCA 1km day data could not be    */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SRCA_1KM_DAY,NULL,start_day,dimsize_cal_day,
                    L1A_scan->SRCA_1km_day) == MFAIL)
    {
      sprintf(msg,"The SRCA 1km day data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }
  

  /*****************************************************************************/  
  /*                                                                           */
  /*  Compute dimsize array for the SRCA 1km night data                        */
  /*  Compute the start for the SRCA 1km night data                            */
  /*  CALL putMODISarray to write the SRCA 1km night data to the L1A granule   */
  /*    INPUTS: MODFILE, M01SRCA_1M_NITE, NULL, start, dimsize array,          */
  /*            SC_SCAN_DATA_t.SRCA_1km_night                                  */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SRCA 1km night data could not be written to the L1A   */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SRCA 1km night data could not be  */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SRCA_1KM_NITE,NULL,start_night,dimsize_cal_night,
                    L1A_scan->SRCA_1km_night) == MFAIL)
    {
      sprintf(msg,"The SRCA 1km night data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the BB 250m data                               */
  /*  Compute the start for the BB 250m data                                   */
  /*  CALL putMODISarray to write the BB 250m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01BB_250M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.BB_250m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The BB 250m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the BB 250m data could not be written */
  /*        to the L1A granule                                                 */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01BB_250M,NULL,start_250,dimsize_cal_250,
                    L1A_scan->BB_250m) == MFAIL)
    {
      sprintf(msg,"The BB 250m data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the BB 500m data                               */
  /*  Compute the start for the BB 500m data                                   */
  /*  CALL putMODISarray to write the BB 500m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01BB_500M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.BB_500m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The BB 500m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the BB 500m data could not be written */
  /*        to the L1A granule                                                 */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01BB_500M,NULL,start_500,dimsize_cal_500,
                    L1A_scan->BB_500m) == MFAIL)
    {
      sprintf(msg,"The BB 500m data could not be written to the L1A granule %s",
            L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the BB 1km day data                            */
  /*  Compute the start for the BB 1km day data                                */
  /*  CALL putMODISarray to write the BB 1km day data to the L1A granule       */ 
  /*    INPUTS: MODFILE, M01BB_1M_DAY, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.BB_1km_day                                      */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The BB 1km day data could not be written to the L1A       */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the BB 1km day data could not be      */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01BB_1KM_DAY,NULL,start_day,dimsize_cal_day,
                    L1A_scan->BB_1km_day) == MFAIL)
    {
      sprintf(msg,"The BB 1km day data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the BB 1km night data                          */
  /*  Compute the start for the BB 1km night data                              */
  /*  CALL putMODISarray to write the BB 1km night data to the L1A granule     */
  /*    INPUTS: MODFILE, M01BB_1M_NITE, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.BB_1km_night                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The BB 1km night data could not be written to the L1A     */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the BB 1km night data could not be    */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01BB_1KM_NITE,NULL,start_night,dimsize_cal_night,
                    L1A_scan->BB_1km_night) == MFAIL)
    {
      sprintf(msg,"The BB 1km night data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SV 250m data                               */
  /*  Compute the start for the SV 250m data                                   */
  /*  CALL putMODISarray to write the SV 250m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01SV_250M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.SV_250m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SV 250m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SV 250m data could not be written */
  /*        to the L1A granule                                                 */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SV_250M,NULL,start_250,dimsize_cal_250,
                    L1A_scan->SV_250m) == MFAIL)
    {
      sprintf(msg,"The SV 250m data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SV 500m data                               */
  /*  Compute the start for the SV 500m data                                   */
  /*  CALL putMODISarray to write the SV 500m data to the L1A granule          */
  /*    INPUTS: MODFILE, M01SV_500M, NULL, start, dimsize array,               */
  /*            SC_SCAN_DATA_t.SV_500m                                         */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SV 500m data could not be written to the L1A granule" */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SV 500m data could not be written */
  /*        to the L1A granule                                                 */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SV_500M,NULL,start_500,dimsize_cal_500,
                    L1A_scan->SV_500m) == MFAIL)
    {
      sprintf(msg,"The SV 500m data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SV 1km day data                            */
  /*  Compute the start for the SV 1km day data                                */
  /*  CALL putMODISarray to write the SV 1km day data to the L1A granule       */
  /*    INPUTS: MODFILE, M01SV_1M_DAY, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.SV_1km_day                                      */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SV 1km day data could not be written to the L1A       */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SV 1km day data could not be      */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SV_1KM_DAY,NULL,start_day,dimsize_cal_day,
                    L1A_scan->SV_1km_day) == MFAIL)
    {
      sprintf(msg,"The SV 1km day data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the SV 1km night data                          */
  /*  Compute the start for the SV 1km night data                              */
  /*  CALL putMODISarray to write the SV 1km night data to the L1A granule     */
  /*    INPUTS: MODFILE, M01SV_1M_NITE, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.SV_1km_night                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The SV 1km night data could not be written to the L1A     */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the SV 1km night data could not be    */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SV_1KM_NITE,NULL,start_night,dimsize_cal_night,
                    L1A_scan->SV_1km_night) == MFAIL)
    {
      sprintf(msg,"The SV 1km night data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  IF scan_type is equal to Day                                             */
  /*  THEN                                                                     */
  /*     Compute dimsize array for the EV 250m data                            */
  /*     Compute the start for the EV 250m data                                */
  /*     CALL putMODISarray to write the EV 250m data to the L1A granule       */
  /*       INPUTS: MODFILE, M01EV_250M, NULL, start, dimsize array,            */
  /*               SC_SCAN_DATA_t.EV_250m                                      */
  /*       OUTPUT: None                                                        */
  /*       RETURN: mapiStatus                                                  */
  /*                                                                           */
  /*     IF mapiStatus equals MFAIL                                            */
  /*     THEN                                                                  */
  /*        set routine to "write_scan_data"                                   */
  /*        set msg to "The EV 250m data could not be written to the L1A       */
  /*           granule"                                                        */
  /*        set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                       */
  /*        CALL log_fmt_msg to report that the EV 250m data could not be      */
  /*           written to the L1A granule                                      */
  /*          INPUTS: returnStatus, routine, msg                               */
  /*          OUTPUT: None                                                     */
  /*          RETURN: None                                                     */
  /*     ENDIF                                                                 */
  /*                                                                           */
  /*****************************************************************************/

  if (strcmp(scan_type,MD_DAY_SCAN) == 0)
    {
     dimsize[0] = PD_DN_NUM_250M_DETECTORS;
     dimsize[1] = PD_DN_NUM_250M_BANDS;
     dimsize[2] = PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_250M;

     if (putMODISarray(L1A_file_ptr, M01EV_250M, NULL, start_250, dimsize,
        L1A_scan->EV_250m) == MFAIL)
       {
        sprintf(msg,"The EV 250m data could not be written to the L1A granule %s",
               L1A_file_ptr->filename);
        returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
        log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
       }


  /*****************************************************************************/
  /*                                                                           */
  /*     Compute dimsize array for the EV 500m data                            */
  /*     Compute the start for the EV 500m data                                */
  /*     CALL putMODISarray to write the EV 500m data to the L1A granule       */
  /*       INPUTS: MODFILE, M01EV_500M, NULL, start, dimsize array,            */
  /*               SC_SCAN_DATA_t.EV_500m                                      */
  /*       OUTPUT: None                                                        */
  /*       RETURN: mapiStatus                                                  */
  /*                                                                           */
  /*     IF mapiStatus equals MFAIL                                            */
  /*     THEN                                                                  */
  /*        set routine to "write_scan_data"                                   */
  /*        set msg to "The EV 500m data could not be written to the L1A       */
  /*           granule"                                                        */
  /*        set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                       */
  /*        CALL log_fmt_msg to report that the EV 500m data could not be      */
  /*           written to the L1A granule                                      */
  /*          INPUTS: returnStatus, routine, msg                               */
  /*          OUTPUT: None                                                     */
  /*          RETURN: None                                                     */
  /*     ENDIF                                                                 */
  /*                                                                           */
  /*****************************************************************************/

     dimsize[0] = PD_DN_NUM_500M_DETECTORS;
     dimsize[1] = PD_DN_NUM_500M_BANDS;
     dimsize[2] = PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_500M;

     if (putMODISarray(L1A_file_ptr, M01EV_500M, NULL, start_500, dimsize,
        L1A_scan->EV_500m) == MFAIL)
       {
        sprintf(msg,"The EV 500m data could not be written to the L1A granule %s",
               L1A_file_ptr->filename);
        returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
        log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
       }


  /*****************************************************************************/
  /*                                                                           */
  /*     Compute dimsize array for the EV 1km day data                         */
  /*     Compute the start for the EV 1km day data                             */
  /*     CALL putMODISarray to write the EV 1km day data to the L1A granule    */
  /*       INPUTS: MODFILE, M01EV_1M_DAY, NULL, start, dimsize array,          */
  /*               SC_SCAN_DATA_t.EV_1km_day                                   */
  /*       OUTPUT: None                                                        */
  /*       RETURN: mapiStatus                                                  */
  /*                                                                           */
  /*     IF mapiStatus equals MFAIL                                            */
  /*     THEN                                                                  */
  /*        set routine to "write_scan_data"                                   */
  /*        set msg to "The EV 1km day data could not be written to the L1A    */
  /*           granule"                                                        */
  /*        set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                       */
  /*        CALL log_fmt_msg to report that the EV 1km day data could not be   */
  /*           written to the L1A granule                                      */
  /*          INPUTS: returnStatus, routine, msg                               */
  /*          OUTPUT: None                                                     */
  /*          RETURN: None                                                     */
  /*     ENDIF                                                                 */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

     dimsize[0] = PD_DN_NUM_1KMDAY_DETECTORS;
     dimsize[1] = PD_DN_NUM_1KMDAY_BANDS;
     dimsize[2] = PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_1KM;

     if (putMODISarray(L1A_file_ptr, M01EV_1KM_DAY, NULL, start_day, dimsize, 
        L1A_scan->EV_1km_day) == MFAIL)
       {
        sprintf(msg,"The EV 1km day data could not be written to the L1A granule %s",
               L1A_file_ptr->filename);
        returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
        log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
       }
    }

  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the EV 1km night data                          */
  /*  Compute the start for the EV 1km night data                              */
  /*  CALL putMODISarray to write the EV 1km night data to the L1A granule     */
  /*    INPUTS: MODFILE, M01EV_1M_NITE, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.EV_1km_night                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The EV 1km night data could not be written to the L1A     */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the EV 1km night data could not be    */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[0] = PD_DN_NUM_1KMNIGHT_DETECTORS;
  dimsize[1] = PD_DN_NUM_1KMNIGHT_BANDS;
  dimsize[2] = PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_1KM;

  if (putMODISarray(L1A_file_ptr,M01EV_1KM_NITE,NULL,start_night,dimsize,
                    L1A_scan->EV_1km_night) == MFAIL)
    {
      sprintf(msg,"The EV 1km night data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the FPA/AEM Config data                        */
  /*  Compute the start for the FPA/AEM Config data                            */
  /*  CALL putMODISarray to write the FPA/AEM Config data to the L1A granule   */
  /*    INPUTS: MODFILE, M01FPA_AEM_CONFIG, NULL, start, dimsize array,        */
  /*            SC_SCAN_DATA_t.fpa_aem_config                                  */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The FPA/AEM Config data could not be written to the L1A   */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the FPA/AEM Config data could not be  */ 
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[0] = 1;
  dimsize[1] = PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS;
  dimsize[2] = 0;

  if (putMODISarray(L1A_file_ptr,M01FPA_AEM_CONFIG,NULL,start,dimsize,
                    L1A_scan->fpa_aem_config) == MFAIL)
    {
      sprintf(msg,"The FPA/AEM Config data could not be written to the L1A granule %s",
            L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Science State data                         */
  /*  Compute the start for the Science State data                             */
  /*  CALL putMODISarray to write the Science State data to the L1A granule    */
  /*    INPUTS: MODFILE, M01SCIENCE_STATE, NULL, start, dimsize array,         */
  /*            SC_SCAN_DATA_t.science_state                                   */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Science State data could not be written to the L1A    */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Science State data could not be   */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/
  dimsize[1] = 0;

  if (putMODISarray(L1A_file_ptr,M01SCIENCE_STATE,NULL,start,dimsize
                    ,&L1A_scan->science_state) == MFAIL)
    {
      sprintf(msg,"The Science State data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Science Abnorm data                        */
  /*  Compute the start for the Science Abnorm data                            */
  /*  CALL putMODISarray to write the Science Abnorm data to the L1A granule   */
  /*    INPUTS: MODFILE, M01SCIENCE_ABNORM, NULL, start, dimsize array,        */
  /*            SC_SCAN_DATA_t.science_abnorm                                  */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Science Abnorm data could not be written to the L1A   */
  /*        granule"                                                           */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Science Abnorm data could not be  */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  if (putMODISarray(L1A_file_ptr,M01SCIENCE_ABNORM,NULL,start,dimsize,
                    &L1A_scan->science_abnormal) == MFAIL)
    {
      sprintf(msg,"The Science Abnorm data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the PV & PC DCR Offset data                    */
  /*  Compute the start for the PV & PC DCR Offset data                        */
  /*  CALL putMODISarray to write the PV & PC DCR Offset data to the L1A       */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01FPA_DCR_OFFSET, NULL, start, dimsize array,        */
  /*            SC_SCAN_DATA_t.fpa_dcr_offset                                  */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The PV & PC DCR Offset data could not be written to the   */
  /*        L1A granule"                                                       */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the PV & PC DCR Offset data could not */
  /*        be written to the L1A granule                                      */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E1P1_NUM_FPA_DCR_OFFSETS;

  if (putMODISarray(L1A_file_ptr,M01FPA_DCR_OFFSET,NULL,start,dimsize,
                    L1A_scan->fpa_dcr_offset) == MFAIL)
    {
      sprintf(msg,"The PV & PC DCR Offset data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Earth Encoder Times data                   */
  /*  Compute the start for the Earth Encoder Times data                       */
  /*  CALL putMODISarray to write the Earth Encoder Times data to the L1A      */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_MIR_ENC, NULL, start, dimsize array,           */
  /*            SC_SCAN_DATA_t.raw_mir_enc                                     */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Earth Encoder Times data could not be written to the  */
  /*        L1A granule"                                                       */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Earth Encoder Times data could    */
  /*        not be written to the L1A granule                                  */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E1P2_NUM_EARTH_ENCODER_TIMES;

  if (putMODISarray(L1A_file_ptr, M01RAW_MIR_ENC, NULL, start, dimsize,
      L1A_scan->raw_mir_enc) == MFAIL)
    {
      sprintf(msg,"The Mirror Encoder Times data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }

  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the View Sector Definitions data               */
  /*  Compute the start for the View Sector Definitions data                   */
  /*  CALL putMODISarray to write the View Sector Definitions data to the L1A  */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_VS_DEF, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.raw_vs_def                                      */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The View Sector Definitions data could not be written to  */
  /*        the L1A granule"                                                   */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the View Sector Definitions data could*/
  /*        could not be written to the L1A granule                            */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS;

  if (putMODISarray(L1A_file_ptr,M01RAW_VS_DEF,NULL,start,dimsize,L1A_scan->raw_vs_def)
      == MFAIL)
    {
      sprintf(msg,"The View Sector Definitions data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the View Sector Actuals data                   */
  /*  Compute the start for the View Sector Actuals data                       */
  /*  CALL putMODISarray to write the View Sector Actuals data to the L1A      */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_VS_ACT, NULL, start, dimsize array,            */
  /*            SC_SCAN_DATA_t.raw_vs_act                                      */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The View Sector Actuals data could not be written to the  */
  /*        L1A granule"                                                       */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the View Sector Actuals data could    */
  /*        not be written to the L1A granule                                  */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*********************************************************************************/

  dimsize[1] = PD_E1P2_NUM_VIEW_SECTOR_ACTUALS;

  if (putMODISarray(L1A_file_ptr,M01RAW_VS_ACT,NULL,start,dimsize,L1A_scan->raw_vs_act)
      == MFAIL)
    {
      sprintf(msg,"The View Sector Actuals data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Eng Data                                   */
  /*  Compute the start for the Eng Data                                       */
  /*  CALL putMODISarray to write the Eng Data to the L1A granule              */
  /*    INPUTS: MODFILE, M01RAW_SCI_ENG, NULL, start, dimsize array,           */
  /*            SC_SCAN_DATA_t.raw_sci_eng                                     */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Eng Data cuold not be written to the L1A granule"     */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Eng Data could not be             */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = SC_NUM_SCI_ENG_BYTES_IN_SDS;

  if (putMODISarray(L1A_file_ptr,M01RAW_SCI_ENG,NULL,start,dimsize,L1A_scan->raw_sci_eng)
      == MFAIL)
    {
      sprintf(msg,"The Raw Science Eng Data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Current/Prior HK Tlmy data                 */
  /*  Compute the start for the Current/Prior HK Tlmy data                     */
  /*  CALL putMODISarray to write the Current/Prior HK Tlmy data to the L1A    */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_HK_TELEM, NULL, start, dimsize array,          */
  /*            SC_SCAN_DATA_t.raw_hk_telem                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Current/Prior HK Tlmy data could not be written to    */
  /*        the L1A granule"                                                   */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Current/Prior HK Tlmy data could  */ 
  /*        not be written to the L1A granule                                  */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E2P1_NUM_HK_TELEM_BYTES;

  if (putMODISarray(L1A_file_ptr,M01RAW_HK_TELEM,NULL,start,dimsize,L1A_scan->raw_hk_telem)
      == MFAIL)
    {
      sprintf(msg,"The Current/Prior HK Tlmy data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Current/Prior S/C Ancill data              */
  /*  Compute the start for the Current/Prior S/C Ancill data                  */
  /*  CALL putMODISarray to write the Current/Prior S/C Ancill data to the L1A */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_SC_ANCIL, NULL, start, dimsize array,          */
  /*            SC_SCAN_DATA_t.raw_sc_ancil                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Current/Prior S/C Ancill data could not be written to */
  /*        the L1A granule"                                                   */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Current/Prior S/C Ancill data     */
  /*        could not be written to the L1A granule                            */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E2P1_NUM_SC_ANCIL_WORDS;

  if (putMODISarray(L1A_file_ptr,M01RAW_SC_ANCIL,NULL,start,dimsize,L1A_scan->raw_sc_ancil)
      == MFAIL)
    {
      sprintf(msg,"The Current/Prior S/C Ancill data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the Command Parameters data                    */
  /*  Compute the start for the Command Parameters data                        */
  /*  CALL putMODISarray to write the Command Parameters data to the L1A       */
  /*     granule                                                               */
  /*    INPUTS: MODFILE, M01RAW_PARAM, NULL, start, dimsize array,             */
  /*            SC_SCAN_DATA_t.raw_param                                       */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The Command Parameters data could not be written to the   */
  /*        L1A granule"                                                       */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the Command Parameters data could not */
  /*        be written to the L1A granule                                      */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E2P1_NUM_PARAM_BYTES;

  if (putMODISarray(L1A_file_ptr,M01RAW_PARAM,NULL,start,dimsize,L1A_scan->raw_param)
      == MFAIL)
    {
      sprintf(msg,"The Command Parameters data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  Compute dimsize array for the PV Gains data                              */
  /*  Compute the start for the PV Gains data                                  */
  /*  CALL putMODISarray to write the PV Gains data to the L1A granule         */
  /*    INPUTS: MODFILE, M01RAW_PV_GAINS, NULL, start, dimsize array,          */
  /*            SC_SCAN_DATA_t.raw_pv_gains                                    */
  /*    OUTPUT: None                                                           */
  /*    RETURN: mapiStatus                                                     */
  /*                                                                           */
  /*  IF mapiStatus equals MFAIL                                               */
  /*  THEN                                                                     */
  /*     set routine to "write_scan_data"                                      */
  /*     set msg to "The PV Gains data could not be written to the L1A granule"*/
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                          */
  /*     CALL log_fmt_msg to report that the PV Gains data could not be        */
  /*        written to the L1A granule                                         */
  /*       INPUTS: returnStatus, routine, msg                                  */
  /*       OUTPUT: None                                                        */
  /*       RETURN: None                                                        */
  /*  ENDIF                                                                    */
  /*                                                                           */
  /*****************************************************************************/

  dimsize[1] = PD_E2P2_NUM_PV_GAINS;

  if (putMODISarray(L1A_file_ptr,M01RAW_PV_GAINS,NULL,start,dimsize,L1A_scan->raw_pv_gains)
      == MFAIL)
    {
      sprintf(msg,"The PV Gains data could not be written to the L1A granule %s",
             L1A_file_ptr->filename);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /*****************************************************************************/
  /*                                                                           */
  /*  return returnStatus                                                      */
  /*                                                                           */
  /*****************************************************************************/

  return (returnStatus);
}
