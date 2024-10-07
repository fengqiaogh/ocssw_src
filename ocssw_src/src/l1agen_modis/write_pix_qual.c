#include "L1A_prototype.h"
#include "mapi.h"
#include "hdf.h"
#include "PGS_MODIS_35005.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"
#include "mapiL1A.h"
#include "PGS_SMF.h"

PGSt_SMF_status  write_pix_qual ( MODFILE                  *L1A_file_ptr, 
                                  SC_PIXEL_QUALITY_DATA_t  *pix_qual,
                                  int16                    scan_num )
/*
!C****************************************************************************

!Description:  This function will write out the Pixel Quality Data (section
               3 of the MODIS Level 1A Data Product Format) to the L1A file.

!Input Parameters:
               MODFILE                  L1A_file_ptr   ** A pointer to the L1A 
                                                          file               **
               SC_PIXEL_QUALITY_DATA_t  pix_qual       ** Pixel quality data
                                                          for this scan      **
               int16                    scan_num       ** scan number        **

!Output Parameters: 
               None

Return Values: MODIS_S_SUCCESS                           (PGS_MODIS_35005.h)
               MODIS_E_ARRAY_OUTPUT_ERR                  (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status                           (PGS_SMF.h)
               SC_PIXEL_QUALITY_DATA_t                   (SC_scan.h)
               MODFILE                                   (mapi.h)
               M01SD_PIX_QUAL                            (mapiL1A.h)
               M01SRCA_PIX_QUAL                          (mapiL1A.h)
               M01BB_PIX_QUAL                            (mapiL1A.h)
               M01SV_PIX_QUAL                            (mapiL1A.h)
               M01EV_PIX_QUAL                            (mapiL1A.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX        (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT    (PH_pkt_hdr.h)

Called By:
               write_scan

Routines Called:
               putMODISarray
               log_fmt_msg

!Revision History:
               revision 1.0 1997/08/29  17:30:00
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
  /****************************************************************************/
  /*                                                                          */
  /*              Define and Initialize Local Variables                       */
  /*                                                                          */
  /****************************************************************************/

  char                      *routine = "write_pix_qual";
  char                      msg[300];
  PGSt_SMF_status           returnStatus;
  long int                  start[3] = {0, 0, 0};
  long int                  dimsize[3] = {1, 0, 0};

  start[0] = scan_num - 1;

  /****************************************************************************/
  /*                                                                          */
  /*  set returnStatus to MODIS_S_SUCCESS                                     */
  /*                                                                          */
  /****************************************************************************/

  returnStatus =  MODIS_S_SUCCESS;


  /****************************************************************************/
  /*                                                                          */
  /*  Compute dimsize array for the SD pixel quality data                     */
  /*  Compute the start for the SD pixel quality data                         */
  /*  CALL putMODISarray to write the SD pixel quality data to the L1A        */
  /*     granule                                                              */
  /*    INPUTS: MODFILE, M01SD_PIX_QUAL, NULL, start, dimsize array,          */
  /*            SC_PIXEL_QUALITY_DATA_t.SD_pix_qual                           */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus equals MFAIL                                              */
  /*  THEN                                                                    */
  /*     set routine to "write_pix_qual"                                      */
  /*     set msg to "The SD pixel quality data could not be written to the    */
  /*        L1A granule"                                                      */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                         */
  /*     CALL log_fmt_msg to report that the SD pixel quality data could not  */
  /*        be written to the L1A granule                                     */
  /*       INPUTS: returnStatus, routine, msg                                 */
  /*       OUTPUT: None                                                       */
  /*       RETURN: None                                                       */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  dimsize[1] = PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX;
  dimsize[2] = PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP;

  if (putMODISarray(L1A_file_ptr, M01SD_PIX_QUAL, NULL, start, dimsize,
                    pix_qual->SD_pix_qual) == MFAIL)
    {
      sprintf(msg,"The SD pixel quality data could not be written to the L1A granule: %s  Scan Number: %d",
              L1A_file_ptr->filename, scan_num);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /****************************************************************************/
  /*                                                                          */
  /*  Compute dimsize array for the SRCA pixel quality data                   */
  /*  Compute the start for the SRCA pixel quality data                       */
  /*  CALL putMODISarray to write the SRCA pixel quality data to the L1A      */
  /*     granule                                                              */
  /*    INPUTS: MODFILE, M01SRCA_PIX_QUAL, NULL, start, dimsize array,        */
  /*            SC_PIXEL_QUALITY_DATA_t.SRCA_pix_qual                         */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus equals MFAIL                                              */
  /*  THEN                                                                    */
  /*     set routine to "write_pix_qual"                                      */
  /*     set msg to "The SRCA pixel quality data could not be written to the  */
  /*        L1A granule"                                                      */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                         */
  /*     CALL log_fmt_msg to report that the SRCA pixel quality data could not*/
  /*        be written to the L1A granule                                     */
  /*       INPUTS: returnStatus, routine, msg                                 */
  /*       OUTPUT: None                                                       */
  /*       RETURN: None                                                       */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SRCA_PIX_QUAL, NULL, start, dimsize,
                    pix_qual->SRCA_pix_qual) == MFAIL)
    {
      sprintf(msg,"The SRCA pixel quality data could not be written to the L1A granule: %s Scan Number: %d",
              L1A_file_ptr->filename, scan_num);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }

 
  /****************************************************************************/
  /*                                                                          */
  /*  Compute dimsize array for the BB pixel quality data                     */
  /*  Compute the start for the BB pixel quality data                         */
  /*  CALL putMODISarray to write the BB pixel quality data to the L1A        */
  /*     granule                                                              */
  /*    INPUTS: MODFILE, M01BB_PIX_QUAL, NULL, start, dimsize array,          */
  /*            SC_PIXEL_QUALITY_DATA_t.BB_pix_qual                           */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus equals MFAIL                                              */
  /*  THEN                                                                    */
  /*     set routine to "write_pix_qual"                                      */
  /*     set msg to "The BB pixel quality data could not be written to the    */
  /*        L1A granule"                                                      */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                         */
  /*     CALL log_fmt_msg to report that the BB pixel quality data could not  */
  /*        be written to the L1A granule                                     */
  /*       INPUTS: returnStatus, routine, msg                                 */
  /*       OUTPUT: None                                                       */
  /*       RETURN: None                                                       */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01BB_PIX_QUAL, NULL, start, dimsize,
                    pix_qual->BB_pix_qual) == MFAIL)
    {
      sprintf(msg,"The BB pixel quality data could not be written to the L1A granule: %s  Scan Number: %d ",
             L1A_file_ptr->filename, scan_num);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /****************************************************************************/
  /*                                                                          */
  /*  Compute dimsize array for the SV pixel quality data                     */
  /*  Compute the start for the SV pixel quality data                         */
  /*  CALL putMODISarray to write the SV pixel quality data to the L1A        */
  /*     granule                                                              */
  /*    INPUTS: MODFILE, M01SV_PIX_QUAL, NULL, start, dimsize array,          */
  /*            SC_PIXEL_QUALITY_DATA_t.SV_pix_qual                           */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus equals MFAIL                                              */
  /*  THEN                                                                    */
  /*     set routine to "write_pix_qual"                                      */
  /*     set msg to "The SV pixel quality data could not be written to the    */
  /*        L1A granule"                                                      */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                         */
  /*     CALL log_fmt_msg to report that the SV pixel quality data could not  */
  /*        be written to the L1A granule                                     */
  /*       INPUTS: returnStatus, routine, msg                                 */
  /*       OUTPUT: None                                                       */
  /*       RETURN: None                                                       */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if (putMODISarray(L1A_file_ptr, M01SV_PIX_QUAL, NULL, start, dimsize,
                    pix_qual->SV_pix_qual) == MFAIL)
    {
      sprintf(msg,"The SV pixel quality data could not be written to the L1A granule: %s  Scan Number: %d",
             L1A_file_ptr->filename, scan_num);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /****************************************************************************/
  /*                                                                          */
  /*  Compute dimsize array for the EV pixel quality data                     */
  /*  Compute the start for the EV pixel quality data                         */
  /*  CALL putMODISarray to write the EV pixel quality data to the L1A        */
  /*     granule                                                              */
  /*    INPUTS: MODFILE, M01EV_PIX_QUAL, NULL, start, dimsize array,          */
  /*            SC_PIXEL_QUALITY_DATA_t.EV_pix_qual                           */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus equals MFAIL                                              */
  /*  THEN                                                                    */
  /*     set routine to "write_pix_qual"                                      */
  /*     set msg to "The EV pixel quality data could not be written to the    */
  /*        L1A granule"                                                      */
  /*     set returnStatus to MODIS_E_ARRAY_OUTPUT_ERR                         */
  /*     CALL log_fmt_msg to report that the EV pixel quality data could not  */
  /*        be written to the L1A granule                                     */
  /*       INPUTS: returnStatus, routine, msg                                 */
  /*       OUTPUT: None                                                       */
  /*       RETURN: None                                                       */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  dimsize[1] = PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT;

  if (putMODISarray(L1A_file_ptr, M01EV_PIX_QUAL, NULL, start, dimsize,
                    pix_qual->EV_pix_qual) == MFAIL)
    {
      sprintf(msg,"The EV pixel quality data could not be written to the L1A granule: %s  Scan Number: %d",
             L1A_file_ptr->filename, scan_num);
      returnStatus = MODIS_E_ARRAY_OUTPUT_ERR;
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
    }


  /****************************************************************************/
  /*                                                                          */
  /*  return returnStatus                                                     */
  /*                                                                          */
  /****************************************************************************/

  return (returnStatus);
}
