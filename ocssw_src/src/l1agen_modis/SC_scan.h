#ifndef SC_SCAN_H 
#define SC_SCAN_H

/*
!C-INC************************************************************************

!Description:  This include file contains the definition of the modis_scan 
               structures.  This includes the Scan Data structure (section
               4 of the MODIS Level 1A Data Product Format) and the Pixel
               Quality Data structure (section 3 of the MODIS Level 1A Data
               Product Format). 

               All definitions in this file will begin with "SC_".  

!Input Parameters: 
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:  
               PD_DN_NUM_250M_DETECTORS                 (PD_pkt_data.h)
               PD_DN_NUM_250M_BANDS                     (PD_pkt_data.h)
               PD_DN_BAND_RATIO_250M                    (PD_pkt_data.h)
               PD_DN_NUM_500M_DETECTORS                 (PD_pkt_data.h)
               PD_DN_NUM_500M_BANDS                     (PD_pkt_data.h)
               PD_DN_BAND_RATIO_500M                    (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_DETECTORS               (PD_pkt_data.h)
               PD_DN_NUM_1KMDAY_BANDS                   (PD_pkt_data.h)
               PD_DN_BAND_RATIO_1KM                     (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_DETECTORS             (PD_pkt_data.h)
               PD_DN_NUM_1KMNIGHT_BANDS                 (PD_pkt_data.h)
               PD_DN_BAND_RATIO_1KM                     (PD_pkt_data.h)
               PD_E1P1_NUM_FPA_DCR_OFFSETS              (PD_pkt_data.h)
               PD_E1P2_NUM_EARTH_ENCODER_TIMES          (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS      (PD_pkt_data.h)
               PD_E1P2_NUM_VIEW_SECTOR_ACTUALS          (PD_pkt_data.h)
               PD_E2P1_NUM_HK_TELEM_BYTES               (PD_pkt_data.h)
               PD_E2P1_NUM_SC_ANCIL_WORDS               (PD_pkt_data.h)
               PD_E2P1_NUM_PARAM_BYTES                  (PD_pkt_data.h)
               PD_E2P2_NUM_PV_GAINS                     (PD_pkt_data.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT   (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX       (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS       (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP        (PH_pkt_hdr.h)
               int8                                     (hdfi.h)
               int16                                    (hdfi.h)

Called By:
               N/A

Routines Called:
               N/A

!Revision History:
$Log: SC_scan.h,v $
Revision 5.1  2004/09/30 18:54:46  seaton
Updated to run Collection 5 L1A code.
seaton@saicmodis.com


               Revision 3.0  2001/04/13
               John Seaton SAIC/GSC
               Changed raw_mir_enc from int16 to uint16
               DDTS GSFcd02139

               Revision 2.1  2000/07/05
               John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
               Added SC_SCAN_PROC_STATE_t type for split scan fix.
               DDTS MODx101733

               Revision 2.0  1997/07/09  14:15
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Created include file from Version 1 include file
               modis_scan.h for Version 2

               Revision 1.0  1996/06/24  14:00 EDT
               Keith Degnan/SAIC/GSC (keith.degnan@gsfc.nasa.gov)
               Created the PDL for module

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
	       HDF portions developed at the National Center for 
               Supercomputing Applications at the University of Illinois 
               at Urbana-Champaign.

!Design Notes: 
               The ".h" file below was specifically written for development
               in C. Any other language choice may require reworking of the
               ".h" file before coding can begin.

!END**********************************************************************
*/

#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "hdfi.h"


/************************************************************************/
/*  The following are the time offset percentages from the start of     */
/*    the scan (ie the SD start time).  These percentages were          */
/*    calculated based on a scan rate of 1.477, but should apply no     */
/*    matter what the scan rate.  Remember, these are only percentages  */
/*    and not the actual time between sectors.                          */
/************************************************************************/

#define  SC_SD_TIME_OFFSET_PERCENTAGE           0.0
#define  SC_SRCA_TIME_OFFSET_PERCENTAGE         0.072742045
#define  SC_BB_TIME_OFFSET_PERCENTAGE           0.13473934
#define  SC_SV_TIME_OFFSET_PERCENTAGE           0.21849695
#define  SC_EV_TIME_OFFSET_PERCENTAGE           0.34462424


#define  SC_SCAN_RATE_TOLERANCE                 0.05


/************************************************************************/
/*  The following are the values of the instrument field of view (IFOV) */
/*  pixel quality flags.                                                */
/************************************************************************/

#define  SC_GOOD_PIXEL_PACKET                   0 
#define  SC_MISSING_PACKET                      1 
#define  SC_BAD_CHECKSUM_PACKET                 2 
#define  SC_DISCARDED_PACKET                    4 


/************************************************************************/
/*  The following are the default values                                */
/************************************************************************/

#define  SC_FILL_VALUE                         -1
#define  TIME_FILL_VALUE                       -2E9

/**********************************************************************/

#define  SC_NUM_SCI_ENG_BYTES_IN_SDS            224


    /******************************************************************/
    /*  The following are the typedefs for calibration target arrays  */
    /******************************************************************/

typedef int16 SC_CAL_250M[PD_DN_NUM_250M_DETECTORS]
                         [PD_DN_NUM_250M_BANDS]
                         [PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * 
                          PD_DN_BAND_RATIO_250M];

typedef int16 SC_CAL_500M[PD_DN_NUM_500M_DETECTORS]
                         [PD_DN_NUM_500M_BANDS]
                         [PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * 
                          PD_DN_BAND_RATIO_500M];

typedef int16 SC_CAL_1KM_DAY[PD_DN_NUM_1KMDAY_DETECTORS]
                            [PD_DN_NUM_1KMDAY_BANDS]
                            [PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * 
                             PD_DN_BAND_RATIO_1KM];

typedef int16 SC_CAL_1KM_NIGHT[PD_DN_NUM_1KMNIGHT_DETECTORS]
                              [PD_DN_NUM_1KMNIGHT_BANDS]
                              [PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX * 
                               PD_DN_BAND_RATIO_1KM];


    /******************************************************************/
    /*  The following are the typedefs for earth target arrays        */
    /******************************************************************/

typedef int16 SC_EV_250M[PD_DN_NUM_250M_DETECTORS]
                        [PD_DN_NUM_250M_BANDS]
                        [PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * 
                         PD_DN_BAND_RATIO_250M];

typedef int16 SC_EV_500M[PD_DN_NUM_500M_DETECTORS]
                        [PD_DN_NUM_500M_BANDS]
                        [PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * 
                         PD_DN_BAND_RATIO_500M];

typedef int16 SC_EV_1KM_DAY[PD_DN_NUM_1KMDAY_DETECTORS]
                           [PD_DN_NUM_1KMDAY_BANDS]
                           [PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * 
                            PD_DN_BAND_RATIO_1KM];

typedef int16 SC_EV_1KM_NIGHT[PD_DN_NUM_1KMNIGHT_DETECTORS]
                             [PD_DN_NUM_1KMNIGHT_BANDS]         
                             [PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT * 
                              PD_DN_BAND_RATIO_1KM];


    /******************************************************************/
    /*  The following is the data structure for the Scan Data         */
    /*    (section 4 of the MODIS Level 1A Data Product Format        */
    /******************************************************************/

typedef struct 
   {
		/* SOLAR DIFFUSER SOURCE (SD) */
   SC_CAL_250M       SD_250m;
   SC_CAL_500M       SD_500m;
   SC_CAL_1KM_DAY    SD_1km_day;
   SC_CAL_1KM_NIGHT  SD_1km_night;

		/* SRCA CALIBRATION SOURCE (SRCA) */
   SC_CAL_250M       SRCA_250m;
   SC_CAL_500M       SRCA_500m;
   SC_CAL_1KM_DAY    SRCA_1km_day;
   SC_CAL_1KM_NIGHT  SRCA_1km_night;

		/* BLACK BODY SOURCE (BB) */
   SC_CAL_250M       BB_250m;
   SC_CAL_500M       BB_500m;
   SC_CAL_1KM_DAY    BB_1km_day;
   SC_CAL_1KM_NIGHT  BB_1km_night;

		/* SPACE VIEW SOURCE (SV) */
   SC_CAL_250M       SV_250m;
   SC_CAL_500M       SV_500m;
   SC_CAL_1KM_DAY    SV_1km_day;
   SC_CAL_1KM_NIGHT  SV_1km_night;

		/* EARTH VIEW SOURCE (EV) */
   SC_EV_250M        EV_250m;
   SC_EV_500M        EV_500m;
   SC_EV_1KM_DAY     EV_1km_day;
   SC_EV_1KM_NIGHT   EV_1km_night;

        /* Data from the engineering packets */

   int8    fpa_aem_config[PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS];
   int8    science_state;
   int8    science_abnormal;
   int8    fpa_dcr_offset[PD_E1P1_NUM_FPA_DCR_OFFSETS];
   int16   raw_mir_enc[PD_E1P2_NUM_EARTH_ENCODER_TIMES];
   int16   raw_vs_def[PD_E1P2_NUM_VIEW_SECTOR_DEFINITIONS];
   int16   raw_vs_act[PD_E1P2_NUM_VIEW_SECTOR_ACTUALS];
   int8    raw_sci_eng[SC_NUM_SCI_ENG_BYTES_IN_SDS];
   int8    raw_hk_telem[PD_E2P1_NUM_HK_TELEM_BYTES];
   int16   raw_sc_ancil[PD_E2P1_NUM_SC_ANCIL_WORDS];
   int8    raw_param[PD_E2P1_NUM_PARAM_BYTES];
   int8    raw_pv_gains[PD_E2P2_NUM_PV_GAINS];

   } SC_SCAN_DATA_t;


typedef struct 
   {
   int16   SD_pix_qual[PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX]
                      [PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP];
   int16   SRCA_pix_qual[PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX]
                      [PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP];
   int16   BB_pix_qual[PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX]
                      [PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP];
   int16   SV_pix_qual[PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX]
                      [PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP];
   int16   EV_pix_qual[PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT]
                      [PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP];

   } SC_PIXEL_QUALITY_DATA_t;

typedef struct
   {
   PGSt_double pkt_TAI_time;
   int8 sector_packet_count;
   } SC_SCAN_PROC_STATE_t;


#endif   /* SC_SCAN_H */

