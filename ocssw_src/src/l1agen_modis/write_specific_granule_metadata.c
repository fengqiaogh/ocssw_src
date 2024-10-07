#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "mapi.h"
#include "hdf.h"
#include "MD_metadata.h"


PGSt_SMF_status write_specific_granule_metadata (MODFILE                *mfile,
                                                 MD_L1A_SPECIFIC_MET_t  *l1a_specific_met)
/*
!C***************************************************************************

!Description:   Function  write_specific_granule_metadata stores the L1A Specific
                Metadata (MD_L1A_SPECIFIC_MET_t) in the L1A file.  
                 
!Input Parameters:
                MODFILE                 mfile               **  Address of MODFILE 
                                                                structure               **
                MD_L1A_SPECIFIC_MET_t   l1a_specific_met    **  L1A Specific Metadata 
                                                                structure               **

!Output Parameters:
                None

Return Values:                
                MODIS_S_SUCCESS                      (PGS_MODIS_35005.h)
                MODIS_E_GATTRIB_FAILED               (PGS_MODIS_35005.h)

Externally Defined:
                MODFILE                              (mapi.h)
                MAPIOK                               (mapi.h)
                MFAIL                                (mapi.h)
                I32                                  (mapi.h)
                TXT                                  (mapi.h)
                MD_L1A_SPECIFIC_MET_t                (MD_metadata.h)
                MD_NUM_SCANS_TXT                     (MD_metadata.h)
                MD_NUM_DAY_SCANS_TXT                 (MD_metadata.h)
                MD_NUM_NIGHT_SCANS_TXT               (MD_metadata.h)
                MD_MAX_TOTAL_FRAMES_TXT              (MD_metadata.h)
                MD_MAX_EARTH_FRAMES_TXT              (MD_metadata.h)
                MD_MAX_SD_FRAMES_TXT                 (MD_metadata.h)
                MD_MAX_SRCA_FRAMES_TXT               (MD_metadata.h)
                MD_MAX_BB_FRAMES_TXT                 (MD_metadata.h)
                MD_MAX_SV_FRAMES_TXT                 (MD_metadata.h)
                MD_SCAN_TYPES_TXT                    (MD_metadata.h)
                MD_INCOMPLETE_SCANS_TXT              (MD_metadata.h)
                MD_MISSING_PKTS_TXT                  (MD_metadata.h)
                MD_BAD_CRC_PKTS_TXT                  (MD_metadata.h)
                MD_DISCARDED_PKTS_TXT                (MD_metadata.h)

Called By:
                write_global_metadata

Routines Called:
                putMODISfileinfo
                log_fmt_msg

!Revision History:

                revision 1.0 1997/08/21  17:30:00
                Qi Huang/RDC    (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support 
                Team (SDST) for the National Aeronautics and Space Administration 
                (NASA), Goddard Space Flight Center (GSFC), under contract 
                NAS5-32373.

!References and Credits:
                None

!Design Notes:
                None

!END***************************************************************************
*/

{
  /****************************************************************************/
  /*                                                                          */
  /*              Define and Initialize Local Variables                       */
  /*                                                                          */
  /****************************************************************************/

  PGSt_SMF_status   returnStatus;
  char              *routine = "write_specific_granule_metadata";
  char              msg[300];
  int               mapiStatus;


  /****************************************************************************/
  /*                                                                          */
  /*  Set routine to "write_specific_granule_metadata" (done during           */
  /*     (declaration)                                                        */
  /*  Set msg to "Unable to set the attribute "                               */
  /*  Set returnStatus to MODIS_S_SUCCESS                                     */
  /*  Set mapiStatus to MAPIOK                                                */
  /*                                                                          */
  /****************************************************************************/

  returnStatus = MODIS_S_SUCCESS;
  mapiStatus = MAPIOK;

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Number of scans attribute            */
  /*    INPUT:  mfile, MD_NUM_SCANS_TXT, I32, 1,                              */
  /*            &MD_L1A_SPECIFIC_MET_t.num_scans                              */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, NUM_SCANS               */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_NUM_SCANS_TXT, I32, 1,
                       &(l1a_specific_met->num_scans))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_NUM_SCANS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg); 

    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Number of Day mode scans attribute   */
  /*    INPUT:  mfile, MD_NUM_DAY_SCANS_TXT, I32, 1,                          */
  /*            &MD_L1A_SPECIFIC_MET_t.num_day_scans                          */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, NUM_DAY_SCANS           */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_NUM_DAY_SCANS_TXT, I32, 1,
                       &(l1a_specific_met->num_day_scans))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_NUM_DAY_SCANS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Number of Night mode scans attribute */
  /*    INPUT:  mfile, MD_NUM_NIGHT_SCANS_TXT, I32, 1,                        */
  /*            &MD_L1A_SPECIFIC_MET_t.num_night_scans                        */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, NUM_NITE_SCANS          */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_NUM_NIGHT_SCANS_TXT, I32, 1,
                       &(l1a_specific_met->num_night_scans))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_NUM_NIGHT_SCANS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max Total Frames attribute           */
  /*    INPUT:  mfile, MD_MAX_TOTAL_FRAMES_TXT, I32, 1,                       */
  /*            &MD_L1A_SPECIFIC_MET_t.max_total_frames                       */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_TOTAL_FRAMES        */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_TOTAL_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_total_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_TOTAL_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max Earth Frames attribute           */
  /*    INPUT:  mfile, MD_MAX_EARTH_FRAMES_TXT, I32, 1,                       */
  /*            &MD_L1A_SPECIFIC_MET_t.max_earth_frames                       */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_EARTH_FRAMES        */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_EARTH_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_earth_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_EARTH_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max SD Frames attribute              */
  /*    INPUT:  mfile, MD_MAX_SD_FRAMES_TXT, I32, 1,                          */
  /*            &MD_L1A_SPECIFIC_MET_t.max_sd_frames                          */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_SD_FRAMES           */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_SD_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_sd_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_SD_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max SRCA Frames attribute            */
  /*    INPUT:  mfile, MD_MAX_SRCA_FRAMES_TXT, I32, 1,                        */
  /*            &MD_L1A_SPECIFIC_MET_t.max_srca_frames                        */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_SRCA_FRAMES         */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_SRCA_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_srca_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_SRCA_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max BB Frames attribute              */
  /*    INPUT:  mfile, MD_MAX_BB_FRAMES_TXT, I32, 1,                          */
  /*            &MD_L1A_SPECIFIC_MET_t.max_bb_frames                          */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_BB_FRAMES           */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_BB_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_bb_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_BB_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Max SV Frames attribute              */
  /*    INPUT:  mfile, MD_MAX_SV_FRAMES_TXT, I32, 1,                          */
  /*            &MD_L1A_SPECIFIC_MET_t.max_sv_frames                          */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MAX_SV_FRAMES           */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MAX_SV_FRAMES_TXT, I32, 1,
                       &(l1a_specific_met->max_sv_frames))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MAX_SV_FRAMES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Scan Types in product attribute      */
  /*    INPUT:  mfile, MD_SCAN_TYPES_TXT, TXT,                                */
  /*            sizeof(MD_L1A_SPECIFIC_MET_t.scan_types_product),             */
  /*            MD_L1A_SPECIFIC_MET_t.scan_types_product                      */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, SCAN_TYPES_PRO          */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_SCAN_TYPES_TXT, TXT, 
                       sizeof(l1a_specific_met->scan_types_product),
                       &(l1a_specific_met->scan_types_product))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_SCAN_TYPES_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Incomplete Scans attribute           */
  /*    INPUT:  mfile, MD_INCOMPLETE_SCANS_TXT, I32, 1,                       */
  /*            &MD_L1A_SPECIFIC_MET_t.incomplete_scans                       */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, INCOMPLETE_SCANS        */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_INCOMPLETE_SCANS_TXT, I32, 1,
                       &(l1a_specific_met->incomplete_scans))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_INCOMPLETE_SCANS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }
 
 
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Missing Packets attribute            */
  /*    INPUT:  mfile, MD_MISSING_PKTS_TXT, I32, 1,                           */
  /*            &MD_L1A_SPECIFIC_MET_t.missing_packets                        */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, MISSING_PKTS            */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_MISSING_PKTS_TXT, I32, 1,
                       &(l1a_specific_met->missing_packets))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_MISSING_PKTS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }
 
 
  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Packets with bad CRC attribute       */
  /*    INPUT:  mfile, MD_BAD_CRC_PKTS_TXT, I32, 1,                           */
  /*            &MD_L1A_SPECIFIC_MET_t.packets_bad_crc                        */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, BAD_CRC_PKTS            */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_BAD_CRC_PKTS_TXT, I32, 1,
                       &(l1a_specific_met->packets_bad_crc))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_BAD_CRC_PKTS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
     }
 

  /****************************************************************************/
  /*                                                                          */
  /*  CALL putMODISfileinfo to store the Discarded packets attribute          */
  /*    INPUT:  mfile, MD_DISCARDED_PKTS_TXT, I32, 1,                         */
  /*            &MD_L1A_SPECIFIC_MET_t.discarded_packets                      */
  /*    OUTPUT: None                                                          */
  /*    RETURN: mapiStatus                                                    */
  /*                                                                          */
  /*  IF mapiStatus is equal to MFAIL                                         */
  /*    Set returnStatus to MODIS_E_GATTRIB_FAILED                            */
  /*    CALL log_fmt_msg to report that the attribute could not be set        */
  /*      INPUT:  mapiStatus, routine, "%s%s\n", msg, DISCARDED_PKTS          */
  /*      OUTPUT: None                                                        */
  /*      RETURN: None                                                        */
  /*  ENDIF                                                                   */
  /*                                                                          */
  /****************************************************************************/

  if ((mapiStatus = putMODISfileinfo(mfile, MD_DISCARDED_PKTS_TXT, I32, 1,
                       &(l1a_specific_met->discarded_packets))) == MFAIL)
    {
      returnStatus = MODIS_E_GATTRIB_FAILED;
      sprintf(msg, "Unable to set the attribute %s", MD_DISCARDED_PKTS_TXT);
      log_fmt_msg(MODIS_E_GATTRIB_FAILED, routine, msg);
    }

  
  /****************************************************************************/
  /*                                                                          */
  /*  RETURN returnStatus                                                     */
  /*                                                                          */
  /****************************************************************************/

  return (returnStatus);

}
  
