#include "L1A_prototype.h"
#include "SC_scan.h"

void   initialize_scan_data(SC_SCAN_DATA_t    *L1A_scan)

/*****************************************************************************
!C

!Description:  This function fills the Scan Data structure (section 4 of
               the MODIS Level 1A Data Product Format) with the appropriate
               fill data (-1).

!Input Parameters: None

!Output Parameters:
               SC_SCAN_DATA_t    *L1A_scan      The Scan Data structure

Return Values: None

Externally Defined:  
               PD_E1P2_NUM_EARTH_ENCODER_TIMES          (PD_pkt_data.h)
               PD_E1P1_NUM_FPA_DCR_OFFSETS              (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_ACTUALS          (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS      (PD_pkt_data.h)
               PD_E2P1_NUM_HK_TELEM_BYTES               (PD_pkt_data.h)
               PD_E2P1_NUM_PARAM_BYTES                  (PD_pkt_data.h)
               PD_E2P1_NUM_SC_ANCIL_WORDS               (PD_pkt_data.h)
               PD_E2P2_NUM_PV_GAINS                     (PD_pkt_data.h)
               PD_DN_BAND_RATIO_250M                    (PD_pkt_data.h)
               PD_DN_BAND_RATIO_500M                    (PD_pkt_data.h)
               PD_DN_BAND_RATIO_1KM                     (PD_pkt_data.h)    
               PD_DN_NUM_1KMDAY_BANDS                   (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_DETECTORS               (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_BANDS                 (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_DETECTORS             (PD_pkt_data.h)
               PD_DN_NUM_250M_BANDS                     (PD_pkt_data.h)
               PD_DN_NUM_250M_DETECTORS                 (PD_pkt_data.h)
               PD_DN_NUM_500M_BANDS                     (PD_pkt_data.h)
               PD_DN_NUM_500M_DETECTORS                 (PD_pkt_data.h)
               PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS       (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX       (PD_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT   (PH_pkt_hdr.h) 
               SC_FILL_VALUE                            (SC_scan.h)
               SC_NUM_SCI_ENG_BYTES_IN_SDS              (SC_scan.h) 
               SC_SCAN_DATA_t                           (SC_scan.h)

Called By:     initialize_scan

Routines Called: None

!Revision History:
  $Log: initialize_scan_data.c,v $
  Revision 4.1  2003/03/07 15:55:05  vlin
  Updated after code walkthrough

  Revision 4.0  2002/12/02 15:14:05  vlin
  called memcmp before using memset
  vlin@saicmodis.com

               Revision 2.2  1997/09/08  18:00 EDT 
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Originated Code.          

!Team-unique Header:

        This software is developed by the MODIS Science Data Support Team 
        for the National Aeronautics and Space Administration, 
        Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:    None

Design Notes:   The CODE below was developed in C language.

!END
**************************************************************************/

{
  int i, j, k;

  L1A_scan->science_state = SC_FILL_VALUE;
  L1A_scan->science_abnormal = SC_FILL_VALUE;
  L1A_scan->fpa_aem_config[0] = SC_FILL_VALUE;
  memset(L1A_scan->fpa_aem_config+1, SC_FILL_VALUE, 
                                     sizeof(L1A_scan->fpa_aem_config[0]));

  if (memcmp(L1A_scan->fpa_aem_config, L1A_scan->fpa_aem_config+1, 
      sizeof(L1A_scan->fpa_aem_config[0]))) {         /* memset didn't work */
      for (i=1; i<PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS; i++)
           L1A_scan->fpa_aem_config[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E1P1_NUM_FPA_DCR_OFFSETS; i++)
           L1A_scan->fpa_dcr_offset[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E1P2_NUM_EARTH_ENCODER_TIMES; i++) 
           L1A_scan->raw_mir_enc[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS; i++)
           L1A_scan->raw_vs_def[i] = SC_FILL_VALUE; 
      for (i=0; i<PD_E1P2_NUM_VIEW_SECTOR_ACTUALS; i++)
           L1A_scan->raw_vs_act[i] = SC_FILL_VALUE;
      for (i=0; i<SC_NUM_SCI_ENG_BYTES_IN_SDS; i++)
           L1A_scan->raw_sci_eng[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E2P1_NUM_HK_TELEM_BYTES; i++)
           L1A_scan->raw_hk_telem[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E2P1_NUM_SC_ANCIL_WORDS; i++)
           L1A_scan->raw_sc_ancil[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E2P1_NUM_PARAM_BYTES; i++)
           L1A_scan->raw_param[i] = SC_FILL_VALUE;
      for (i=0; i<PD_E2P2_NUM_PV_GAINS; i++)
           L1A_scan->raw_pv_gains[i] = SC_FILL_VALUE;

      for (i=0;i<PD_DN_NUM_250M_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_250M_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX*PD_DN_BAND_RATIO_250M;k++) {
           L1A_scan->SD_250m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SRCA_250m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->BB_250m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SV_250m[i][j][k] = SC_FILL_VALUE;
      }

      for (i=0;i<PD_DN_NUM_250M_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_250M_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT*PD_DN_BAND_RATIO_250M;k++)
           L1A_scan->EV_250m[i][j][k] = SC_FILL_VALUE;
          
      for (i=0;i<PD_DN_NUM_500M_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_500M_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX*PD_DN_BAND_RATIO_500M;k++) { 
           L1A_scan->SD_500m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SRCA_500m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->BB_500m[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SV_500m[i][j][k] = SC_FILL_VALUE;
      }

      for (i=0;i<PD_DN_NUM_500M_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_500M_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT*PD_DN_BAND_RATIO_500M;k++)
           L1A_scan->EV_500m[i][j][k] = SC_FILL_VALUE;

      for (i=0;i<PD_DN_NUM_1KMDAY_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_1KMDAY_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX*PD_DN_BAND_RATIO_1KM;k++) { 
           L1A_scan->SD_1km_day[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SRCA_1km_day[i][j][k] = SC_FILL_VALUE;
           L1A_scan->BB_1km_day[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SV_1km_day[i][j][k] = SC_FILL_VALUE;
      }

      for (i=0;i<PD_DN_NUM_1KMDAY_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_1KMDAY_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT*PD_DN_BAND_RATIO_1KM;k++)
           L1A_scan->EV_1km_day[i][j][k] = SC_FILL_VALUE;

      for (i=0;i<PD_DN_NUM_1KMNIGHT_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_1KMNIGHT_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX*PD_DN_BAND_RATIO_1KM;k++) { 
           L1A_scan->SD_1km_night[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SRCA_1km_night[i][j][k] = SC_FILL_VALUE;
           L1A_scan->BB_1km_night[i][j][k] = SC_FILL_VALUE;
           L1A_scan->SV_1km_night[i][j][k] = SC_FILL_VALUE;
      }

      for (i=0;i<PD_DN_NUM_1KMNIGHT_DETECTORS;i++)
      for (j=0;j<PD_DN_NUM_1KMNIGHT_BANDS;j++)
      for (k=0;k<PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT*PD_DN_BAND_RATIO_1KM;k++)
           L1A_scan->EV_1km_night[i][j][k] = SC_FILL_VALUE;
  }
  else {
     memset(L1A_scan->SD_250m, SC_FILL_VALUE, sizeof(L1A_scan->SD_250m));
     memset(L1A_scan->SD_500m, SC_FILL_VALUE, sizeof(L1A_scan->SD_500m));
     memset(L1A_scan->SD_1km_day, SC_FILL_VALUE, sizeof(L1A_scan->SD_1km_day));
     memset(L1A_scan->SD_1km_night, SC_FILL_VALUE, 
                                    sizeof(L1A_scan->SD_1km_night));
     memset(L1A_scan->SRCA_250m, SC_FILL_VALUE, sizeof(L1A_scan->SRCA_250m));
     memset(L1A_scan->SRCA_500m, SC_FILL_VALUE, sizeof(L1A_scan->SRCA_500m));
     memset(L1A_scan->SRCA_1km_day, SC_FILL_VALUE, 
                                    sizeof(L1A_scan->SRCA_1km_day));
     memset(L1A_scan->SRCA_1km_night, SC_FILL_VALUE, 
                                      sizeof(L1A_scan->SRCA_1km_night));
     memset(L1A_scan->BB_250m, SC_FILL_VALUE, sizeof(L1A_scan->BB_250m));
     memset(L1A_scan->BB_500m, SC_FILL_VALUE, sizeof(L1A_scan->BB_500m));
     memset(L1A_scan->BB_1km_day, SC_FILL_VALUE, sizeof(L1A_scan->BB_1km_day));
     memset(L1A_scan->BB_1km_night, SC_FILL_VALUE, 
                                    sizeof(L1A_scan->BB_1km_night));
     memset(L1A_scan->SV_250m, SC_FILL_VALUE, sizeof(L1A_scan->SV_250m));
     memset(L1A_scan->SV_500m, SC_FILL_VALUE, sizeof(L1A_scan->SV_500m));
     memset(L1A_scan->SV_1km_day, SC_FILL_VALUE, sizeof(L1A_scan->SV_1km_day));
     memset(L1A_scan->SV_1km_night, SC_FILL_VALUE, 
                                    sizeof(L1A_scan->SV_1km_night));
     memset(L1A_scan->EV_250m, SC_FILL_VALUE, sizeof(L1A_scan->EV_250m));
     memset(L1A_scan->EV_500m, SC_FILL_VALUE, sizeof(L1A_scan->EV_500m));
     memset(L1A_scan->EV_1km_day, SC_FILL_VALUE, sizeof(L1A_scan->EV_1km_day));
     memset(L1A_scan->EV_1km_night, SC_FILL_VALUE, 
                                    sizeof(L1A_scan->EV_1km_night));
     memset(L1A_scan->fpa_aem_config+2, SC_FILL_VALUE, 
            sizeof(L1A_scan->fpa_aem_config)-2*sizeof(L1A_scan->fpa_aem_config[0]));
     memset(L1A_scan->fpa_dcr_offset, SC_FILL_VALUE,
                                      sizeof(L1A_scan->fpa_dcr_offset));
     memset(L1A_scan->raw_mir_enc, SC_FILL_VALUE,
                                   sizeof(L1A_scan->raw_mir_enc));
     memset(L1A_scan->raw_vs_def, SC_FILL_VALUE, sizeof(L1A_scan->raw_vs_def));
     memset(L1A_scan->raw_vs_act, SC_FILL_VALUE, sizeof(L1A_scan->raw_vs_act));
     memset(L1A_scan->raw_sci_eng, SC_FILL_VALUE, 
                                   sizeof(L1A_scan->raw_sci_eng));
     memset(L1A_scan->raw_hk_telem, SC_FILL_VALUE,
                                    sizeof(L1A_scan->raw_hk_telem));
     memset(L1A_scan->raw_sc_ancil, SC_FILL_VALUE,
                                    sizeof(L1A_scan->raw_sc_ancil));
     memset(L1A_scan->raw_param, SC_FILL_VALUE,
                                 sizeof(L1A_scan->raw_param));
     memset(L1A_scan->raw_pv_gains, SC_FILL_VALUE,
                                    sizeof(L1A_scan->raw_pv_gains));
  }
} 
