#include "L1A_prototype.h"
#include "hdf.h"
#include "mapi.h"
#include "mapiL1A.h"
#include "PGS_MODIS_35005.h"
#include "PGS_SMF.h"
#include "SC_scan.h"
#include "PD_pkt_data.h"
#include "PH_pkt_hdr.h"


PGSt_SMF_status   init_L1A_scan_data_HDF_sdss (MODFILE  *mfile, 
                                               int      nscans)

/*
!C***************************************************************************

!Description:   Function init_L1A_scan_data_HDF_sdss creates all the scan 
                SDSs (arrays) that will be written to by L1A processing.
              
!Input Parameters:
                mfile          **  Pointer to MODFILE structure  **
                nscans         **  Number of scans  **
!Output Parameters:
                None
    

Return Values:        
                MODIS_S_SUCCESS                        (PGS_MODIS_35005.h)
                MODIS_E_SDS_CREATE_FAILED              (PGS_MODIS_35005.h)

Externally Defined:
                PGSt_SMF_status                        (PGS_SMF.h)
                MO1SD_250M                             (mapiL1A.h)
                MO1SD_500M                             (mapiL1A.h)
                MO1SD_1KM_DA                           (mapiL1A.h)
                MO1SD_1KM_NITE                         (mapiL1A.h)
                MO1SRCA_250M                           (mapiL1A.h)
                MO1SRCA_500M                           (mapiL1A.h)
                MO1SRCA_1KM_DAY                        (mapiL1A.h)
                MO1SRCA_1KM_NITE                       (mapiL1A.h)
                M01BB_250M                             (mapiL1A.h)
                M01BB_500M                             (mapiL1A.h)
                M01BB_1KM_DAY                          (mapiL1A.h)
                M01BB_1KM_NITE                         (mapiL1A.h)
                M01SV_250M                             (mapiL1A.h)
                M01SV_500M                             (mapiL1A.h)
                M01SV_1KM_DA                           (mapiL1A.h)
                M01SV_1KM_NITE                         (mapiL1A.h)
                M01EV_250M                             (mapiL1A.h)
                M01EV_500M                             (mapiL1A.h)
                M01EV_1KM_DAY                          (mapiL1A.h)
                M01EV_1KM_NITE                         (mapiL1A.h)
                M01RAW_MIR_ENC                         (mapiL1A.h)
                M01RAW_VS_DEF                          (mapiL1A.h)
                M01RAW_VS_ACT                          (mapiL1A.h)
                M01RAW_SC_ANCIL                        (mapiL1A.h)
                M01FPA_AEM_CONFIG                      (mapiL1A.h)
                M01SCIENCE_STATE                       (mapiL1A.h)
                M01SCIENCE_ABNORMAL                    (mapiL1A.h)
                M01FPA_DCR_OFFSET                      (mapiL1A.h)
                M01RAW_SCI_ENG                         (mapiL1A.h)
                M01RAW_HK_TELEM                        (mapiL1A.h)
                M01RAW_PARAM                           (mapiL1A.h)
                M01RAW_PV_GAINS                        (mapiL1A.h)
                MFAIL                                  (mapi.h)
                MODFILE                                (mapi.h)
                DATATYPELENMAX                         (mapi.h)
                I16                                    (mapi.h)
                I8                                     (mapi.h)
                PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX     (PH_pkt_hdr.h)
                PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT (PH_pkt_hdr.h)
                PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS     (PH_pkt_hdr.h)
                SC_NUM_SCI_ENG_BYTES_IN_SDS            (SC_scan.h)
                PD_E1P1_NUM_FPA_DCR_OFFSETS            (PD_pkt_data.h)
                PD_E1P2_NUM_EARTH_ENCODER_TIMES        (PD_pkt_data.h)
                PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS    (PD_pkt_data.h)
                PD_E1P2_NUM_VIEW_SECTOR_ACTUALS        (PD_pkt_data.h)
                PD_E2P1_NUM_HK_TELEM_BYTES             (PD_pkt_data.h)
                PD_E2P1_NUM_SC_ANCIL_WORDS             (PD_pkt_data.h)
                PD_E2P1_NUM_PARAM_BYTES                (PD_pkt_data.h)
                PD_E2P2_NUM_PV_GAINS                   (PD_pkt_data.h)
                PD_DN_NUM_250M_DETECTORS               (PD_pkt_data.h)
                PD_DN_NUM_250M_BANDS                   (PD_pkt_data.h)
                PD_DN_BAND_RATIO_250M                  (PD_pkt_data.h)
                PD_DN_NUM_500M_DETECTORS               (PD_pkt_data.h)
                PD_DN_NUM_500M_BANDS                   (PD_pkt_data.h)
                PD_DN_BAND_RATIO_500M                  (PD_pkt_data.h)
                PD_DN_NUM_1KMDAY_DETECTORS             (PD_pkt_data.h)
                PD_DN_NUM_1KMDAY_BANDS                 (PD_pkt_data.h)
                PD_DN_BAND_RATIO_1KM                   (PD_pkt_data.h)
                PD_DN_NUM_1KMNIGHT_DETECTORS           (PD_pkt_data.h)
                PD_DN_NUM_1KMNIGHT_BANDS               (PD_pkt_data.h)


Called By:
                init_L1A_HDF_sdss

Routines Called:
                log_fmt_msg
                createMODISarray
                putMODISarinfo

!Revision History:
$Log: init_L1A_scan_data_HDF_sdss.c,v $
Revision 5.1  2004/09/23 18:41:06  seaton
Collection 5 changes.
Added code to set a _FillValue attribute for all SDSs in the output product.
seaton@saicmodis.com


                Revision 3.0  2001/04/13
                John Seaton SAIC/GSC
                Set raw_mir_enc to type UI16

                Revision 1.1  1999/10/18  12:03:00
                John Seaton/GSC (seaton@ltpmail.gsfc.nasa.gov)
                Modified to initialize in a loop and to set the fill value for 
                some SDSs to -1 as specified in the L1A filespec.

                Revision 1.0  1997/09/09  17:30:00
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

#define SDS_SECTOR_DATA_INDEX 20
#define NUM_SDSS              32

{
/******************************************************************************/
/*                Define and Initialize Local Variables                       */
/******************************************************************************/

  char              *routine = "init_L1A_HDF_scan_data_sdss";
  char              msg[300];
  PGSt_SMF_status   returnStatus;
  int               i;

/******************************************************************************/
/* Initialize sds structure storing sds names, datatype, rank and dimensions  */
/* for all Scan data as specified in Section 4 of the L1A filespec. The first */
/* index of the dimension member is initially set to zero. It is initialized  */
/* to the correct value later in thie routine. It cannot be initialized here  */
/* because the number of scans is needed to compute that size, and nscans is  */
/* not a constant, therefore, it cannot be used in the initialization of these*/
/* SDSs here.                                                                 */
/******************************************************************************/

  struct {
    char *sds_name;
    char datatype[DATATYPELENMAX];
    long int rank;
    long int dimensions[3];
  } sds_info[] = {
     {M01SD_250M,     I16, 3, {0, PD_DN_NUM_250M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_250M}},
     {M01SD_500M,     I16, 3, {0, PD_DN_NUM_500M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_500M}},
     {M01SD_1KM_DAY,  I16, 3, {0, PD_DN_NUM_1KMDAY_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01SD_1KM_NITE, I16, 3, {0, PD_DN_NUM_1KMNIGHT_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01SRCA_250M,     I16, 3, {0, PD_DN_NUM_250M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_250M}},
     {M01SRCA_500M,     I16, 3, {0, PD_DN_NUM_500M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_500M}},
     {M01SRCA_1KM_DAY,  I16, 3, {0, PD_DN_NUM_1KMDAY_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01SRCA_1KM_NITE, I16, 3, {0, PD_DN_NUM_1KMNIGHT_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01BB_250M,     I16, 3, {0, PD_DN_NUM_250M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_250M}},
     {M01BB_500M,     I16, 3, {0, PD_DN_NUM_500M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_500M}},
     {M01BB_1KM_DAY,  I16, 3, {0, PD_DN_NUM_1KMDAY_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01BB_1KM_NITE, I16, 3, {0, PD_DN_NUM_1KMNIGHT_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01SV_250M,     I16, 3, {0, PD_DN_NUM_250M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_250M}},
     {M01SV_500M,     I16, 3, {0, PD_DN_NUM_500M_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_500M}},
     {M01SV_1KM_DAY,  I16, 3, {0, PD_DN_NUM_1KMDAY_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01SV_1KM_NITE, I16, 3, {0, PD_DN_NUM_1KMNIGHT_BANDS, PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * PD_DN_BAND_RATIO_1KM}},
     {M01EV_250M,     I16, 3, {0, PD_DN_NUM_250M_BANDS, PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_250M}},
     {M01EV_500M,     I16, 3, {0, PD_DN_NUM_500M_BANDS, PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_500M}},
     {M01EV_1KM_DAY,  I16, 3, {0, PD_DN_NUM_1KMDAY_BANDS, PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_1KM}},
     {M01EV_1KM_NITE, I16, 3, {0, PD_DN_NUM_1KMNIGHT_BANDS, PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * PD_DN_BAND_RATIO_1KM}},
     {M01FPA_AEM_CONFIG, I8, 2, {0, PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS, 0}},
     {M01SCIENCE_STATE, I8, 1, {0, 0, 0}},
     {M01SCIENCE_ABNORM, I8, 1, {0, 0, 0}},
     {M01FPA_DCR_OFFSET, I8, 2, {0, PD_E1P1_NUM_FPA_DCR_OFFSETS, 0}},
     {M01RAW_MIR_ENC, I16, 2, {0, PD_E1P2_NUM_EARTH_ENCODER_TIMES, 0}},
     {M01RAW_VS_DEF, I16, 2, {0, PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS, 0}},
     {M01RAW_VS_ACT, I16, 2, {0, PD_E1P2_NUM_VIEW_SECTOR_ACTUALS, 0}},
     {M01RAW_SCI_ENG, I8, 2, {0, SC_NUM_SCI_ENG_BYTES_IN_SDS, 0}},
     {M01RAW_HK_TELEM, I8, 2, {0, PD_E2P1_NUM_HK_TELEM_BYTES, 0}},
     {M01RAW_SC_ANCIL, I16, 2, {0, PD_E2P1_NUM_SC_ANCIL_WORDS, 0}},
     {M01RAW_PARAM, I8, 2, {0, PD_E2P1_NUM_PARAM_BYTES, 0}},
     {M01RAW_PV_GAINS, I8, 2, {0, PD_E2P2_NUM_PV_GAINS, 0}}
  };

/******************************************************************************/
/* Set the fill value variable to -1 for future calls to API routines         */
/******************************************************************************/

  int16 fill16[1] = {-1};
  int8 fill8[1] = {-1};

/******************************************************************************/
/* Set the first dimension for the SDSs as they could not be initialized due  */
/* to the fact that the nscans variable is not a constant.                    */
/******************************************************************************/

  sds_info[0].dimensions[0] = nscans * PD_DN_NUM_250M_DETECTORS;
  sds_info[1].dimensions[0] = nscans * PD_DN_NUM_500M_DETECTORS;
  sds_info[2].dimensions[0] = nscans * PD_DN_NUM_1KMDAY_DETECTORS;
  sds_info[3].dimensions[0] = nscans * PD_DN_NUM_1KMNIGHT_DETECTORS;
  sds_info[4].dimensions[0] = nscans * PD_DN_NUM_250M_DETECTORS;
  sds_info[5].dimensions[0] = nscans * PD_DN_NUM_500M_DETECTORS;
  sds_info[6].dimensions[0] = nscans * PD_DN_NUM_1KMDAY_DETECTORS;
  sds_info[7].dimensions[0] = nscans * PD_DN_NUM_1KMNIGHT_DETECTORS;
  sds_info[8].dimensions[0] = nscans * PD_DN_NUM_250M_DETECTORS;
  sds_info[9].dimensions[0] = nscans * PD_DN_NUM_500M_DETECTORS;
  sds_info[10].dimensions[0] = nscans * PD_DN_NUM_1KMDAY_DETECTORS;
  sds_info[11].dimensions[0] = nscans * PD_DN_NUM_1KMNIGHT_DETECTORS;
  sds_info[12].dimensions[0] = nscans * PD_DN_NUM_250M_DETECTORS;
  sds_info[13].dimensions[0] = nscans * PD_DN_NUM_500M_DETECTORS;
  sds_info[14].dimensions[0] = nscans * PD_DN_NUM_1KMDAY_DETECTORS;
  sds_info[15].dimensions[0] = nscans * PD_DN_NUM_1KMNIGHT_DETECTORS;
  sds_info[16].dimensions[0] = nscans * PD_DN_NUM_250M_DETECTORS;
  sds_info[17].dimensions[0] = nscans * PD_DN_NUM_500M_DETECTORS;
  sds_info[18].dimensions[0] = nscans * PD_DN_NUM_1KMDAY_DETECTORS;
  sds_info[19].dimensions[0] = nscans * PD_DN_NUM_1KMNIGHT_DETECTORS;
  sds_info[20].dimensions[0] = nscans;
  sds_info[21].dimensions[0] = nscans;
  sds_info[22].dimensions[0] = nscans;
  sds_info[23].dimensions[0] = nscans;
  sds_info[24].dimensions[0] = nscans;
  sds_info[25].dimensions[0] = nscans;
  sds_info[26].dimensions[0] = nscans;
  sds_info[27].dimensions[0] = nscans;
  sds_info[28].dimensions[0] = nscans;
  sds_info[29].dimensions[0] = nscans;
  sds_info[30].dimensions[0] = nscans;
  sds_info[31].dimensions[0] = nscans;


  returnStatus = MODIS_S_SUCCESS;

/******************************************************************************/
/* Loop through the above defined structure, creating the SDSs.               */
/******************************************************************************/

  for (i=0; i < NUM_SDSS; i++) {

/******************************************************************************/
/* call createMODISarray for each SDS defined above. If return from           */
/* createMODISarrray fails, then print LogStatus message stating that the     */
/* named SDS could not be created and set returnStatus to                     */
/* MODIS_E_SDS_CREATE_FAILED.                                                 */
/******************************************************************************/

    if ( createMODISarray(mfile, sds_info[i].sds_name, 
                                 NULL, 
                                 sds_info[i].datatype, 
                                 sds_info[i].rank,
                                 sds_info[i].dimensions) != MAPIOK) {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "SDS Name: %s", sds_info[i].sds_name);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
    }


/******************************************************************************/
/* Then call putMODISarinfo to set its fill value to -1. If this call fails   */
/* then print an appropriate error message in the LogStatus file and set the  */
/* returnStatus to MODIS_E_SDS_CREATE_FAILED.                                 */
/******************************************************************************/

    if (strncmp(sds_info[i].datatype, I16, sizeof(I16)) == 0) {
      if (putMODISarinfo(mfile, sds_info[i].sds_name, 
                                NULL, MFILL_VALUE, 
                                sds_info[i].datatype, 1, &fill16) != MAPIOK) {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "Could not initialize SDS to -1 for SDS Name: %s", sds_info[i].sds_name);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
      }

    } /* end if i  = I16 */
    else {
      if (putMODISarinfo(mfile, sds_info[i].sds_name,
                                NULL, MFILL_VALUE,
                                sds_info[i].datatype, 1, &fill8) != MAPIOK) {
        returnStatus = MODIS_E_SDS_CREATE_FAILED;
        sprintf(msg, "Could not initialize SDS to -1 for SDS Name: %s", sds_info[i].sds_name);
        log_fmt_msg(MODIS_E_SDS_CREATE_FAILED, routine, msg);
      }

    } /* end else */

  }  /* end for */

  /****************************************************************************/
  /* RETURN returnStatus                                                      */
  /****************************************************************************/

  return (returnStatus);
}
