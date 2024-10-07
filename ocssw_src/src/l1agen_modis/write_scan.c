#include "L1A_prototype.h"
#include "mapi.h"
#include "PGS_MODIS_35005.h"
#include "SC_scan.h"
#include "hdf.h"
#include "MD_metadata.h"
#include "EN_eng_data.h"
#include "FP_failed_pkt_queue.h"
#include "PGS_SMF.h"
#include "hfile.h"

PGSt_SMF_status  write_scan ( MODFILE                  *L1A_file_ptr, 
                              SC_SCAN_DATA_t           *L1A_scan, 
                              MD_SCAN_MET_t            *scan_meta, 
                              SC_PIXEL_QUALITY_DATA_t  *pix_qual, 
                              FP_QUEUE_t               failed_pkts, 
                              EN_VDATA_TYPE_t          *eng_data )

/*
!C**********************************************************************

!Description:  This function is a driver for other 'write' functions.

!Input Parameters:
               MODFILE                  *L1A_file_ptr   ** A pointer to the L1A 
                                                           file                **
               SC_SCAN_DATA_t           *L1A_scan       ** L1A scan data for 
                                                           this scan           **
               MD_SCAN_MET_t            scan_meta       ** Scan metadata for 
                                                           this scan           **
               SC_PIXEL_QUALITY_DATA_t  pix_qual        ** The pixel quality 
                                                           data for this scan  **
               FP_QUEUE_t               failed_pkts     ** The discarded packets 
                                                           for this scan       **
               EN_VDATA_TYPE_t          eng_data        ** Engineering data for 
                                                           this scan           **

!Output Parameters: 
               None

Return Values: 
               MODIS_S_SUCCESS              (PGS_MODIS_35005.h)
               MODIS_E_WRITE_SCAN_FAIL      (PGS_MODIS_35005.h)
               MODIS_F_WRITE_ENG_DATA_F     (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status              (PGS_SMF.h)
               MODFILE                      (mapi.h)
               SC_SCAN_DATA_t               (SC_scan.h)
               SC_PIXEL_QUALITY_DATA_t      (SC_scan.h)
               MD_SCAN_MET_t                (MD_metadata.h)
               EN_VDATA_TYPE_t              (EN_eng_data.h)
               FP_QUEUE_t                   (FP_failed_pkt_queue.h)

Called By:
               process_a_granule

Routines Called:
               write_scan_metadata
               write_pix_qual
               write_scan_data
               write_failed_packets
               write_eng_data
               log_fmt_msg

!Revision History:
               Revision 2.0 1999/08/20  11:33:00
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added code to not write past the vdata block limit for
               the discarded packet vdata, and free memory used by
               non-written discarded packets.

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

!END************************************************************************
*/
{
  PGSt_SMF_status        returnStatus;
  char                   *routine = "write_scan";
  char                   msg[300];
  int                    temp_status;
  int                    scan_num;
    
  /**********************************************************************************/
  /*                                                                                */
  /*  set returnStatus to MODIS_S_SUCCESS                                           */
  /*                                                                                */
  /**********************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /**********************************************************************************/
  /*                                                                                */
  /*  CALL write_scan_metadata to write out the Scan Level Metadata (section 2 of   */
  /*     the MODIS Level 1A Data Product Format)                                    */
  /*    INPUTS:  L1A_file_ptr, MD_SCAN_MET_t                                        */
  /*    OUTPUTS: None                                                               */
  /*    RETURN:  temp_status                                                        */
  /*                                                                                */
  /**********************************************************************************/

  temp_status =  write_scan_metadata(L1A_file_ptr, scan_meta);


  /**********************************************************************************/
  /*                                                                                */
  /*  IF temp_status is not equal to MODIS_S_SUCCESS                                */
  /*  THEN                                                                          */
  /*     set msg to "The Scan Level Metadata could not be written to the L1A        */
  /*        granule for scan number (number)"                                       */
  /*     CALL log_fmt_msg to report that the Scan Level Metadata could not be       */
  /*        written to the L1A file                                                 */
  /*       INPUTS: temp_status, routine, msg                                        */
  /*       OUTPUT: None                                                             */
  /*       RETURN: None                                                             */
  /*     set returnStatus to MODIS_E_WRITE_SCAN_FAIL                                */
  /*  ENDIF                                                                         */
  /*                                                                                */
  /**********************************************************************************/

  if (temp_status != MODIS_S_SUCCESS)
    {
      sprintf(msg,"The Scan Level Metadata could not be written to the L1A granule: %s scan number: %d",
              L1A_file_ptr->filename, scan_meta->scan_num);
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
      returnStatus = MODIS_E_WRITE_SCAN_FAIL;
    }


  /**********************************************************************************/
  /*                                                                                */
  /*  set scan_num to MD_SCAN_MET_t.scan_num                                        */
  /*                                                                                */
  /**********************************************************************************/

  scan_num = scan_meta->scan_num;


  /**********************************************************************************/
  /*                                                                                */
  /*  CALL write_pix_qual to write out the Pixel Quality Data (section 3 of the     */
  /*     MODIS Level 1A Data Product Format)                                        */
  /*    INPUTS:  L1A_file_ptr, SC_PIXEL_QUALITY_DATA_t, scan_num                    */
  /*    OUTPUTS: None                                                               */
  /*    RETURN:  temp_status                                                        */
  /*                                                                                */
  /*  IF temp_status is not equal to MODIS_S_SUCCESS                                */
  /*  THEN                                                                          */
  /*     set msg to "The Pixel Quality Data could not be written to the L1A granule */
  /*        for scan number (number)"                                               */
  /*     CALL log_fmt_msg to report that the Pixel Quality Data could not           */
  /*        be written to the L1A file                                              */
  /*       INPUTS: temp_status, routine, msg                                        */
  /*       OUTPUT: None                                                             */
  /*       RETURN: None                                                             */
  /*     set returnStatus to MODIS_E_WRITE_SCAN_FAIL                                */
  /*  ENDIF                                                                         */
  /*                                                                                */
  /**********************************************************************************/

  temp_status = write_pix_qual(L1A_file_ptr, pix_qual, scan_num);

  if (temp_status != MODIS_S_SUCCESS)
    {
      sprintf(msg,"The Pixel Quality Data could not be written to the L1A granule: %s  scan number %d",
             L1A_file_ptr->filename, scan_num);
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
       returnStatus = MODIS_E_WRITE_SCAN_FAIL;
    }


  /**********************************************************************************/
  /*                                                                                */
  /*  CALL write_scan_data to write out the Scan Data (section 4 of the             */
  /*     MODIS Level 1A Data Product Format)                                        */
  /*    INPUTS:  L1A_file_ptr, SC_SCAN_DATA_t, scan_num                             */
  /*    OUTPUTS: None                                                               */
  /*    RETURN:  temp_status                                                        */
  /*                                                                                */
  /*  IF temp_status is not equal to MODIS_S_SUCCESS                                */
  /*  THEN                                                                          */
  /*     set msg to "The Scan Data could not be written to the L1A granule for      */
  /*        scan number (number)"                                                   */
  /*     CALL log_fmt_msg to report that the Scan Data could not be written to the  */
  /*        L1A file                                                                */
  /*       INPUTS: temp_status, routine, msg                                        */
  /*       OUTPUT: None                                                             */
  /*       RETURN: None                                                             */
  /*     set returnStatus to MODIS_E_WRITE_SCAN_FAIL                                */
  /*  ENDIF                                                                         */
  /*                                                                                */
  /**********************************************************************************/

  temp_status = write_scan_data(L1A_file_ptr, L1A_scan, scan_num, scan_meta->scan_type);

  if (temp_status != MODIS_S_SUCCESS)
    {
      sprintf(msg,"The Scan Data could not be written to the L1A granule: %s  scan number: %d",
             L1A_file_ptr->filename, scan_num);
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
      returnStatus = MODIS_E_WRITE_SCAN_FAIL;
    }

  /**********************************************************************************/
  /*                                                                                */
  /* CALL write_failed_packets to write out the Discarded Packets (section 5 of the */
  /*    MODIS Level 1A Data Product Format)                                         */
  /*   INPUTS:  L1A_file_ptr, failed_pkts                                           */
  /*   OUTPUTS: None                                                                */
  /*   RETURN:  temp_status                                                         */
  /*                                                                                */
  /**********************************************************************************/

  temp_status = write_failed_packets(failed_pkts);

  /**********************************************************************************/
  /*                                                                                */
  /*  IF temp_status is not equal to MODIS_S_SUCCESS                                */
  /*  THEN                                                                          */
  /*     set msg to "The Discarded Packets could not be written to the L1A granule  */
  /*        for scan number (number)"                                               */
  /*     CALL log_fmt_msg to report that the Discarded Packets could not be written */
  /*        to the L1A file                                                         */
  /*       INPUTS: temp_status, routine, msg                                        */
  /*       OUTPUT: None                                                             */
  /*       RETURN: None                                                             */
  /*     set returnStatus to MODIS_E_WRITE_SCAN_FAIL                                */
  /*  ENDIF                                                                         */
  /*                                                                                */
  /**********************************************************************************/

  if (temp_status !=  MODIS_S_SUCCESS)
    {
      sprintf(msg,"The Discarded Packets could not be written to the L1A granule: %s  scan number: %d",
             L1A_file_ptr->filename, scan_num);
      log_fmt_msg(MODIS_E_ARRAY_OUTPUT_ERR, routine, msg);
      returnStatus = MODIS_E_WRITE_SCAN_FAIL;
    }
    

  /**********************************************************************************/
  /*                                                                                */
  /*  CALL write_eng_data to write out the Engineering Data (section 6 of the MODIS */
  /*     Level 1A Data Product Format)                                              */
  /*    INPUTS:  L1A_file_ptr, eng_data                                             */
  /*    OUTPUTS: None                                                               */
  /*    RETURN:  temp_status                                                        */
  /*                                                                                */
  /*  IF temp_status is not equal to MODIS_S_SUCCESS                                */
  /*  THEN                                                                          */
  /*     set routine to "write_scan"                                                */
  /*     set msg to "The Engineering Data could not be written to the L1A granule   */
  /*        for scan number (number)"                                               */
  /*     CALL log_fmt_msg to report that the Engineering Data could not be written  */
  /*        to the L1A file                                                         */
  /*       INPUTS: temp_status, routine, msg                                        */
  /*       OUTPUT: None                                                             */
  /*       RETURN: None                                                             */
  /*     set returnStatus to temp_status                                            */
  /*  ENDIF                                                                         */
  /*                                                                                */
  /**********************************************************************************/
  
  temp_status = write_eng_data(eng_data);

  if (temp_status != MODIS_S_SUCCESS)
    {
      sprintf(msg,"The Engineering Data could not be written to the L1A granule: %s  scan number %d",
             L1A_file_ptr->filename, scan_num);
      log_fmt_msg(MODIS_F_WRITE_ENG_DATA_FAIL, routine, msg);
      returnStatus = temp_status;
    }
    

  /**********************************************************************************/
  /*                                                                                */
  /*  return returnStatus                                                           */
  /*                                                                                */
  /**********************************************************************************/

  return (returnStatus);
}





