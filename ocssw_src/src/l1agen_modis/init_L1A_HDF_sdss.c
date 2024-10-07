#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "mapi.h"
#include "PGS_SMF.h"


PGSt_SMF_status  init_L1A_HDF_sdss ( MODFILE *L1A_file_ptr,
                                     int16   nscans )

/*
!C*****************************************************************************

!Description:  This function calls other procedures to create all of the SDSs 
               that will be written to by the L1A processing.

!Input Parameters:
               MODFILE   *L1A_file_ptr     **  Pointer to the L1A file          **
               int16     nscans            **  Number of scans in this granule  **

!Output Parameters:
               None

Return Values: 
               MODIS_S_SUCCESS           (PGS_MODIS_35005.h)
               MODIS_E_SDS_CREATE_FAILED (PGS_MODIS_35005.h)

Externally Defined:  
               None

Called By:
               create_L1A_granule

Routines Called:
               init_L1A_scan_meta_HDF_sdss
               init_L1A_pix_qual_HDF_sdss
               init_L1A_scan_data_HDF_sdss
               log_fmt_msg  

!Revision History:
               revision 1.0 1997/09/05  17:30:00
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
  /*****************************************************************************/
  /*                                                                           */
  /*              Define and Initialize Local Variables                        */
  /*                                                                           */
  /*****************************************************************************/

  char                     *routine = "init_L1A_HDF_sdss";
  PGSt_SMF_status          returnStatus;


  /*****************************************************************************/
  /*                                                                           */
  /* Set returnStatus to MODIS_S_SUCCESS                                       */
  /*                                                                           */
  /*****************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /*****************************************************************************/
  /*                                                                           */
  /* CALL init_L1A_scan_meta_HDF_sdss to create all of the SDSs for the Scan   */
  /*    Level  Metadata (section 2 of the MODIS Level 1A Data Product Format)  */
  /*   INPUT:  L1A_file_ptr, nscans                                            */
  /*   OUTPUT: None                                                            */
  /*   RETURN: returnStatus                                                    */
  /*                                                                           */
  /*****************************************************************************/

  returnStatus = init_L1A_scan_meta_HDF_sdss(L1A_file_ptr, nscans);


  /*****************************************************************************/
  /*                                                                           */
  /* IF returnStatus is not equal to MODIS_S_SUCCESS                           */
  /* THEN                                                                      */
  /*    CALL log_fmt_msg to report that not all of the Scan Level Metadata     */
  /*       SDSs could be created                                               */
  /*      INPUT:  returnStatus, routine, msg                                   */
  /*      OUTPUT: None                                                         */
  /*      RETURN: None                                                         */
  /* ELSE                                                                      */
  /*                                                                           */
  /*    CALL init_L1A_pix_qual_HDF_sdss to create all of the SDSs for the      */
  /*       Pixel Quality Data (section 3 of the MODIS Level 1A Data Product    */
  /*       Format)                                                             */
  /*      INPUT:  L1A_file_ptr, nscans                                         */
  /*      OUTPUT: None                                                         */
  /*      RETURN: returnStatus                                                 */
  /*                                                                           */
  /*    IF returnStatus is not equal to MODIS_S_SUCCESS                        */
  /*    THEN                                                                   */
  /*       CALL log_fmt_msg to report that not all of the Pixel Quality Data   */
  /*          SDSs could be created                                            */
  /*         INPUT:  returnStatus, routine, msg                                */
  /*         OUTPUT: None                                                      */
  /*         RETURN: None                                                      */
  /*    ELSE                                                                   */
  /*                                                                           */
  /*       CALL init_L1A_scan_data_HDF_sdss to create all of the SDSs for the  */
  /*          Scan Data (section 4 of the MODIS Level 1A Data Product Format)  */
  /*         INPUT:  L1A_file_ptr, nscans                                      */
  /*         OUTPUT: None                                                      */
  /*         RETURN: returnStatus                                              */
  /*                                                                           */
  /*       IF returnStatus is not equal to MODIS_S_SUCCESS                     */
  /*       THEN                                                                */
  /*          CALL log_fmt_msg to report that not all of the Scan Data SDSs    */
  /*             could be created                                              */
  /*            INPUT:  returnStatus, routine, msg                             */
  /*            OUTPUT: None                                                   */
  /*            RETURN: None                                                   */
  /*       ENDIF                                                               */
  /*    ENDIF                                                                  */
  /* ENDIF                                                                     */
  /*                                                                           */
  /*****************************************************************************/

  if (returnStatus != MODIS_S_SUCCESS)
    log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, 
                "Not all of the Scan Level Metadata SDSs could be created");
  else
    {
      returnStatus = init_L1A_pix_qual_HDF_sdss(L1A_file_ptr, nscans);
      if (returnStatus != MODIS_S_SUCCESS)
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, 
                    "Not all of the Pixel Quality Data SDSs could be created");
      else
        {
          returnStatus = init_L1A_scan_data_HDF_sdss(L1A_file_ptr, nscans);
          if (returnStatus != MODIS_S_SUCCESS)
            log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine,
                        "Not all of the Scan Data SDSs could be created");
        }
    }

  /*****************************************************************************/
  /*                                                                           */
  /* RETURN returnStatus                                                       */
  /*                                                                           */
  /*****************************************************************************/

  return (returnStatus);
}
