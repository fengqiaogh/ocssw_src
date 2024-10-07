#include "L1A_prototype.h"
#include "PGS_MODIS_35005.h"
#include "mapi.h"
#include "mapiL1A.h"
#include "PH_pkt_hdr.h"
#include "PGS_SMF.h"


PGSt_SMF_status   init_L1A_pix_qual_HDF_sdss (MODFILE  *mfile, 
                                              int      nscans)

/*
!C*****************************************************************************

!Description:   Function init_L1A_pix_qual_HDF_sdss creates all the Pixel Quality
                Data (section 3 of the MODIS Level 1A Data Product Format) SDSs 
                (arrays) that will be written to by L1A processing.
              
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
                M01SD_PIX_QUAL             (mapiL1A.h)
                M01SRCA_PIX_QUAL           (mapiL1A.h)
                M01BB_PIX_QUAL             (mapiL1A.h)
                M01SV_PIX_QUAL             (mapiL1A.h)
                M01EV_PIX_QUAL             (mapiL1A.h)
                MFAIL                      (mapi.h)
                MODFILE                    (mapi.h)
                DATATYPELENMAX             (mapi.h)
                I16                        (mapi.h)
                MFILL_VALUE                (mapi.h)

Called By:
                init_L1A_HDF_sdss

Routines Called:
                log_fmt_msg
                createMODISarray
                putMODISarinfo

!Revision History:
$Log: init_L1A_pix_qual_HDF_sdss.c,v $
Revision 5.1  2004/09/23 18:30:14  seaton
Collection 5 changes. Updated SDSs to have a _FillValue attribute in the product output file.
seaton@saicmodis.com


                Revision 1.0  1997/09/05  17:30:00
                Qi Huang/RDC (qhuang@ltpmail.gsfc.nasa.gov)
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
  /*************************************************************************/
  /*                                                                       */
  /*              Define and Initialize Local Variables                    */
  /*                                                                       */
  /*************************************************************************/

  char                *routine = "init_L1A_pix_qual_HDF_sdss";
  char                msg[300];
  char                data_type[DATATYPELENMAX] = I16;
  long int            rank = 3;
  char                array_name[5][256] = {M01SD_PIX_QUAL, M01SRCA_PIX_QUAL,
                                             M01BB_PIX_QUAL, M01SV_PIX_QUAL,
                                             M01EV_PIX_QUAL};
  long int            dimsizes[5][3] = {{0,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2},
                                        {0,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2},
                                        {0,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2},
                                        {0,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2},
                                        {0,PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT,2}};
  PGSt_SMF_status     returnStatus;
  int                 i;
  int16               fill_val[1]={1};

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
  /* Set data_type to "int16"            (done during declaration)               */
  /*                                                                             */
  /* Set rank to 3                       (done during declaration)               */
  /*                                                                             */
  /* Store all array names (M01SD_PIX_QUAL, M01SRCA_PIX_QUAL, M01BB_PIX_QUAL,    */
  /*    M01SV_PIX_QUAL, M01EV_PIX_QUAL) in array_name   (done during declaration)*/
  /*                                                                             */
  /* Store all dimension sizes (nscans,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2;     */
  /*    nscans,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2;                             */
  /*    nscans,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2;                             */
  /*    nscans,PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX,2;                             */
  /*    nscans,PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT,2) in dimsizes             */
  /*                                                                             */
  /*******************************************************************************/

  dimsizes[0][0] = dimsizes[1][0] = dimsizes[2][0] = dimsizes[3][0] = dimsizes[4][0] 
  = nscans;

 
  /*******************************************************************************/
  /*                                                                             */
  /* FOR i equals all elements of Pixel Quality Data (5)                         */
  /*    CALL createMODISarray to create the appropriate SDS in the L1A file      */
  /*      INPUT:  mfile, array_name[i], NULL, data_type, rank, dimsizes[i]       */
  /*      OUTPUT: None                                                           */
  /*      RETURN: mapi_status                                                    */
  /*                                                                             */
  /*    IF mapi_status is equal to MFAIL                                         */
  /*    THEN                                                                     */
  /*      Set returnStatus to MODIS_E_SDS_CREATE_FAILED                          */
  /*      Set msg to "unable to create the 'array_name' SDS"                     */
  /*      Set routine to "init_L1A_HDF_sdss"                                     */
  /*      CALL log_fmt_msg to report that the 'array_name' SDS was not created   */
  /*        INPUT:  returnStatus, routine, msg                                   */
  /*        OUTPUT: None                                                         */
  /*        RETURN: None                                                         */
  /*    ENDIF                                                                    */
  /*                                                                             */
  /* ENDFOR                                                                      */
  /*                                                                             */
  /*******************************************************************************/

  for (i=0; i<5; i++) {
    if (createMODISarray(mfile,array_name[i],NULL,data_type,rank,dimsizes[i])
         == MFAIL)
      {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "SDS name: %s", array_name[i]);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
      }

    if (putMODISarinfo(mfile, array_name[i],
                                NULL, MFILL_VALUE,
                                data_type, 1, &fill_val) != MAPIOK) {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "Could not initialize SDS to %d for SDS Name: %s", *fill_val, array_name[i]);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
      }

  } /* end for */


  /*******************************************************************************/
  /*                                                                             */
  /* RETURN returnStatus                                                         */
  /*                                                                             */
  /*******************************************************************************/

  return (returnStatus);
}
