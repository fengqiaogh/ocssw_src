#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "mapi.h"
#include "mapiL1A.h"
#include "PGS_SMF.h"
#include "SC_scan.h"


PGSt_SMF_status   init_L1A_scan_meta_HDF_sdss (MODFILE  *mfile, 
                                               int      nscans)

/*
!C*****************************************************************************

!Description:   Function init_L1A_scan_meta_HDF_sdss creates all the Scan Level
                Metadata (section 2 of the MODIS Level 1A Data Product Format
                SDSs (arrays) that will be written to by L1A processing.
              
!Input Parameters:
                MODFILE mfile          **  Pointer to MODFILE structure  **
                int     nscans         **  Number of scans  **

!Output Parameters:
                None

Return Values:        
                MODIS_S_SUCCESS            (PGS_MODIS_35005.h)
                MODIS_E_SDS_CREATE_FAILED  (PGS_MODIS_35005.h)        

Externally Defined:
                PGSt_SMF_status            (PGS_SMF.h)
                M01SCAN_NUMBER             (mapiL1A.h)
                M01FRAME_COUNT_ARRAY       (mapiL1A.h)
                M01SCAN_TYPE               (mapiL1A.h)
                M01SD_START_TIME           (mapiL1A.h)
                M01SRCA_START_TIME         (mapiL1A.h)
                M01BB_START_TIME           (mapiL1A.h)
                M01SV_START_TIME           (mapiL1A.h)
                M01EV_START_TIME           (mapiL1A.h)
                M01SRCA_CALIBRATION_MODE   (mapiL1A.h)
                M01PACKET_SCAN_COUNT       (mapiL1A.h)
                M01CCSDS_APID              (mapiL1A.h)
                M01PACKET_QL               (mapiL1A.h)
                M01MIRROR_SIDE             (mapiL1A.h)
                M01SCAN_QUALITY_ARRAY      (mapiL1A.h)
                MAPIOK                     (mapi.h)
                MFAIL                      (mapi.h)
                MODFILE                    (mapi.h)
                DATATYPELENMAX             (mapic.h)
                MFILL_VALUE                (mapi.h)

Called By:
                init_L1A_HDF_sdss

Routines Called:
                log_fmt_msg
                createMODISarray
                putMODISarinfo

!Revision History:
$Log: init_L1A_scan_meta_HDF_sdss.c,v $
Revision 6.2  2015/08/21 20:13:17  jkuyper
Finished resolving Bug 1492 by setting two additional _FillValue attribute
length to 1.

Revision 6.1  2010/08/24 15:36:55  kuyper
Corrected length of _FillValue attribute to 1.

Revision 5.1  2004/09/23 18:48:46  seaton
Collection 5 changed. Added _FillValue attribute to SDSs.
seaton@saicmodis.com

                Revision 1.0  1997/09/10  17:30:00
                Qi Huang/RDC (qhuang@ltpmail.gsfc.nasa.gov)
                Original development

!Team-unique Header:
                This software is developed by the MODIS Science Data Support 
                Team (SDST) for the National Aeronautics and Space Administration
                (NASA), Goddard Space Flight Center (GSFC), under contract 
                NAS5-32373.

References and Credits:
                None

Design Notes:
                None

!END***************************************************************************
*/
{
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/

  char              *routine = "init_L1A_scan_meta_HDF_sdss";
  char              msg[300];
  char              data_type[14][DATATYPELENMAX] = {I16,I16,TXT,R64,R64,R64,R64,
                                                     R64,I16,I16,I16,I16,I16,I32};
  long int          rank[14] = {1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2};
  long int          dimsizes[14][2] ={{0,0},{0,6},{0,10},{0,0},{0,0},{0,0},{0,0},
                                      {0,0},{0,0},{0,0},{0,3},{0,0},{0,0},{0,4}};
  char              array_name[14][256] = {M01SCAN_NUMBER, M01FRAME_COUNT_ARRAY,
                                           M01SCAN_TYPE, M01SD_START_TIME,
                                           M01SRCA_START_TIME, M01BB_START_TIME,
                                           M01SV_START_TIME, M01EV_START_TIME,
                                           M01SRCA_CALIBRATION_MODE,
                                           M01PACKET_SCAN_COUNT, M01CCSDS_APID,
                                           M01PACKET_QL, M01MIRROR_SIDE,
                                           M01SCAN_QUALITY_ARRAY};
  int               i;
  PGSt_SMF_status   returnStatus;
  int16             fill_zero = 0;
  float64           fill_time = TIME_FILL_VALUE;
  int16             fill_neg_one = SC_FILL_VALUE;

  /*******************************************************************************/
  /*                                                                             */
  /* Set mapi_status to MAPIOK                                                   */
  /*                                                                             */
  /* Set returnStatus to MODIS_S_SUCCESS                                         */
  /*                                                                             */
  /*******************************************************************************/

  returnStatus = MODIS_S_SUCCESS;


  /*******************************************************************************/
  /*                                                                             */
  /* Store all array names (M01SCAN_NUMBER, M01FRAME_COUNT_ARRAY, M01SCAN_TYPE,  */
  /*    M01SD_START_TIME, M01SRCA_START_TIME, M01BB_START_TIME, M01SV_START_TIME,*/
  /*    M01EV_START_TIME, M01SRCA_CALIBRATION_MODE, M01PACKET_SCAN_COUNT,        */
  /*    M01CCSDS_APID, M01PACKET_QL, M01MIRROR_SIDE, M01SCAN_QUALITY_ARRAY) in   */
  /*    array_name   (done during declaration)                                   */
  /*                                                                             */
  /* Store all data types (int16, int16, char, float64, float64, float64,        */
  /*    float64, float64, int16, int16, int16, int16, int16, int32) in data_type */
  /*    (done during declaration)                                                */
  /*                                                                             */
  /* Store all dimension sizes (nscans; nscans,6; nscans,10; nscans; nscans;     */
  /*    nscans; nscans; nscans; nscans; nscans; nscans,3; nscans; nscans;        */
  /*    nscans,4) in dimsizes                                                    */
  /*                                                                             */
  /* Store all ranks (1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2) in rank    (done */
  /*    during declaration)                                                      */
  /*                                                                             */
  /*******************************************************************************/

  dimsizes[0][0] = dimsizes[1][0] = dimsizes[2][0] = dimsizes[3][0] = dimsizes[4][0] =
  dimsizes[5][0] = dimsizes[6][0] = dimsizes[7][0] = dimsizes[8][0] = dimsizes[9][0] = 
  dimsizes[10][0] = dimsizes[11][0] = dimsizes[12][0] = dimsizes[13][0] = nscans;

 
  /*******************************************************************************/
  /*                                                                             */
  /* FOR i equals all elements of Scan Level Metadata (14)                       */
  /*    CALL createMODISarray to create the appropriate SDS in the L1A file      */
  /*      INPUT:  mfile, array_name[i], NULL, data_type[i], rank[i], dimsizes[i] */
  /*      OUTPUT: None                                                           */
  /*      RETURN: mapi_status                                                    */
  /*                                                                             */
  /*    IF mapi_status is equal to MFAIL                                         */
  /*    THEN                                                                     */
  /*      Set returnStatus to MODIS_E_SDS_CREATE_FAILED                          */
  /*      Set msg to "unable to create the 'array_name' SDS"                     */
  /*      Set routine to "init_L1A_scan_meta_HDF_sdss"                           */
  /*      CALL log_fmt_msg to report that the 'array_name' SDS was not created   */
  /*        INPUT:  returnStatus, routine, msg                                   */
  /*        OUTPUT: None                                                         */
  /*        RETURN: None                                                         */
  /*    ENDIF                                                                    */
  /*                                                                             */
  /* ENDFOR                                                                      */
  /*                                                                             */
  /*******************************************************************************/

  for (i=0; i<14; i++) {
    if (createMODISarray(mfile, array_name[i], NULL, data_type[i], rank[i],
                         dimsizes[i]) == MFAIL)
      {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "SDS Name: %s", array_name[i]);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
      }

    switch (i) {
      case 0: if (putMODISarinfo(mfile, array_name[i], NULL, MFILL_VALUE,
		  data_type[i], 1, &fill_zero) != MAPIOK)
	      {
                  returnStatus = MODIS_E_SDS_CREATE_FAILED;
                  sprintf(msg, "Could not initialize SDS to 0 for SDS Name: %s", array_name[i]);
                  log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
              }
              break;
      case 3:
      case 4:
      case 5:
      case 6:
      case 7: if (putMODISarinfo(mfile, array_name[i], NULL, MFILL_VALUE,
		    data_type[i], 1, &fill_time) != MAPIOK)
	      {
                  returnStatus = MODIS_E_SDS_CREATE_FAILED;
                  sprintf(msg, "Could not initialize SDS to 0 for SDS Name: %s", array_name[i]);
                  log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
              }
              break;
      case 8:
      case 9:
      case 10:
      case 11:
      case 12: if (putMODISarinfo(mfile, array_name[i], NULL, MFILL_VALUE,
		    data_type[i], 1, &fill_neg_one) != MAPIOK)
	       {
                  returnStatus = MODIS_E_SDS_CREATE_FAILED;
                  sprintf(msg, "Could not initialize SDS to 0 for SDS Name: %s", array_name[i]);
                  log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
               }
               break;

    } /* end switch */

  } /* end for */

  /*******************************************************************************/
  /*                                                                             */
  /* RETURN returnStatus                                                         */
  /*                                                                             */
  /*******************************************************************************/

  return returnStatus;
}
