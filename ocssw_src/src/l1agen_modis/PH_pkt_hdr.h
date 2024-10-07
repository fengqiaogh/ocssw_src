#ifndef PH_PKT_HDR_H
#define PH_PKT_HDR_H

/*
!C-INC************************************************************************

!Description:  This include file contains the definitions of values and byte
               and bit locations for the primary, secondary and MODIS headers.

               All definitions in this file will begin with "PH_".  
               Definitions for the primary header will begin with "PH_PRI".  
               Definitions for the secondary header will begin with "PH_SEC".  
               Definitions for the MODIS header will begin with "PH_MOD".  

!Input Parameters:
               N/A

!Output Parameters:
               N/A

Return Values: 
               N/A

Externally Defined:  
               PGSt_scTime             (PGS_TD.h)
               PGSt_double             (PGS_TD.h)
               int8                    (hdfi.h)
               int16                   (hdfi.h)

Called By:
               N/A

Routines Called:
               N/A

!Revision History:
               Revision 2.1  1999/08/20  10:25
               John Seaton/SAIC/GSC (seaton@ltpmail.gsfc.nasa.gov)
               Added #define PH_SEC_PKT_TYPE_SPARE 3 so all the packet header,
               packet type spare values will be checked.

               Revision 2.0  1997/07/10  11:00
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Created include file from Version 1 include files 
               modis_pkt_hdr_pos.h and modis_pkt_hdr_vals.h

               Revision 1.0  1996/06/21  13:00 EDT
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

#include "PGS_TD.h"
#include "hdfi.h"

/*********************************************************/
/*  The following general information applies to the     */
/*  primary, secondary, and MODIS header.                */
/*********************************************************/

#define  PH_NUM_12BIT_WORDS_IN_HEADER         12


#define  PH_REALISTIC_NUM_SD_PACKETS         100       
#define  PH_REALISTIC_NUM_SRCA_PACKETS        20       
#define  PH_REALISTIC_NUM_BB_PACKETS         100       
#define  PH_REALISTIC_NUM_SV_PACKETS         100       
#define  PH_REALISTIC_NUM_EV_DAY_PACKETS    2708
#define  PH_REALISTIC_NUM_EV_NIGHT_PACKETS  1354
#define  PH_REALISTIC_NUM_ENG_PACKETS          4       


#define  PH_REALISTIC_NUM_SD_FRAMES           50       
#define  PH_REALISTIC_NUM_SRCA_FRAMES         10       
#define  PH_REALISTIC_NUM_BB_FRAMES           50       
#define  PH_REALISTIC_NUM_SV_FRAMES           50       
#define  PH_REALISTIC_NUM_EV_FRAMES         1354


/*********************************************************/
/*  The following byte/bit offset and value information  */
/*  is for the primary header.                           */
/*********************************************************/

#define  PH_PRI_VERSION_BYTE_OFFSET            0
#define  PH_PRI_VERSION_BIT_OFFSET             0
#define  PH_PRI_VERSION_NUM_BITS               3
#define  PH_PRI_VERSION_VALUE                  0   

#define  PH_PRI_TYPE_BYTE_OFFSET               0
#define  PH_PRI_TYPE_BIT_OFFSET                3
#define  PH_PRI_TYPE_NUM_BITS                  1
#define  PH_PRI_TYPE_VALUE                     0
  
#define  PH_PRI_SEC_HDR_FLAG_BYTE_OFFSET       0
#define  PH_PRI_SEC_HDR_FLAG_BIT_OFFSET        4
#define  PH_PRI_SEC_HDR_FLAG_NUM_BITS          1
#define  PH_PRI_SEC_HDR_PRESENT                1

#define  PH_PRI_APID_BYTE_OFFSET               0
#define  PH_PRI_APID_BIT_OFFSET                5
#define  PH_PRI_APID_NUM_BITS                 11
#define  PH_PRI_MIN_MODIS_APID_AM1            64 
#define  PH_PRI_MAX_MODIS_APID_AM1           127
#define  PH_PRI_APID_TEST_PACKET             127

#define  PH_PRI_SEQUENCE_FLAG_BYTE_OFFSET      2
#define  PH_PRI_SEQUENCE_FLAG_BIT_OFFSET       0
#define  PH_PRI_SEQUENCE_FLAG_NUM_BITS         2
#define  PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP    1
#define  PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP   2
#define  PH_PRI_SEQUENCE_ONLY_PKT_IN_GROUP     3
#define  PH_PRI_SEQUENCE_NOT_USED              0

#define  PH_PRI_SOURCE_SEQ_CNT_BYTE_OFFSET     2
#define  PH_PRI_SOURCE_SEQ_CNT_BIT_OFFSET      2
#define  PH_PRI_SOURCE_SEQ_CNT_NUM_BITS       14
#define  PH_PRI_MAX_CCSDS_PKT_SEQ_COUNT    16383

#define  PH_PRI_PKT_LENGTH_BYTE_OFFSET         4
#define  PH_PRI_PKT_LENGTH_BIT_OFFSET          0
#define  PH_PRI_PKT_LENGTH_NUM_BITS           16
#define  PH_PRI_LONG_PKT_LENGTH              635
#define  PH_PRI_SHORT_PKT_LENGTH             269


/*********************************************************/
/*  The following byte/bit offset and value information  */
/*  is for the secondary header.                         */
/*********************************************************/

#define  PH_SEC_TIME_TAG_BYTE_OFFSET           6
#define  PH_SEC_TIME_TAG_NUM_BYTES             8

#define  PH_SEC_QUICK_LOOK_FLAG_BYTE_OFFSET   14
#define  PH_SEC_QUICK_LOOK_FLAG_BIT_OFFSET     0
#define  PH_SEC_QUICK_LOOK_FLAG_NUM_BITS       1
#define  PH_SEC_QUICK_LOOK_FLAG_SET            1
#define  PH_SEC_QUICK_LOOK_FLAG_NOT_SET        0

#define  PH_SEC_PKT_TYPE_BYTE_OFFSET          14
#define  PH_SEC_PKT_TYPE_BIT_OFFSET            1
#define  PH_SEC_PKT_TYPE_NUM_BITS              3
#define  PH_SEC_PKT_TYPE_DAY_GROUP             0
#define  PH_SEC_PKT_TYPE_NIGHT_GROUP           1
#define  PH_SEC_PKT_TYPE_ENG1_GROUP            2
#define  PH_SEC_PKT_TYPE_ENG2_GROUP            4
#define  PH_SEC_PKT_TYPE_SPARE                 3
#define  PH_SEC_PKT_TYPE_MAX_PKTS_IN_GROUP     2

#define  PH_SEC_SCAN_CNT_BYTE_OFFSET          14
#define  PH_SEC_SCAN_CNT_BIT_OFFSET            4
#define  PH_SEC_SCAN_CNT_NUM_BITS              3
#define  PH_SEC_SCAN_CNT_MAX                   7

#define  PH_SEC_MIRROR_SIDE_BYTE_OFFSET       14
#define  PH_SEC_MIRROR_SIDE_BIT_OFFSET         7
#define  PH_SEC_MIRROR_SIDE_NUM_BITS           1
#define  PH_SEC_MIRROR_SIDE_1                  0
#define  PH_SEC_MIRROR_SIDE_2                  1


/*********************************************************/
/*  The following byte/bit offset and value information  */
/*  is for the MODIS header.                             */
/*********************************************************/

#define  PH_MOD_SOURCE_ID_TYPE_FLAG_BYTE_OFFSET          15
#define  PH_MOD_SOURCE_ID_TYPE_FLAG_BIT_OFFSET            0
#define  PH_MOD_SOURCE_ID_TYPE_FLAG_NUM_BITS              1
#define  PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH                 0
#define  PH_MOD_SOURCE_ID_TYPE_FLAG_CAL                   1

#define  PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BYTE_OFFSET    15
#define  PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BIT_OFFSET      1
#define  PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_NUM_BITS       11
#define  PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_MAX          2048
#define  PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT        1400

#define  PH_MOD_SOURCE_ID_CAL_TYPE_BYTE_OFFSET           15
#define  PH_MOD_SOURCE_ID_CAL_TYPE_BIT_OFFSET             1
#define  PH_MOD_SOURCE_ID_CAL_TYPE_NUM_BITS               2
#define  PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE  0
#define  PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE        1
#define  PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE       2
#define  PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE           3

#define  PH_MOD_SOURCE_ID_CAL_MODE_BYTE_OFFSET           15
#define  PH_MOD_SOURCE_ID_CAL_MODE_BIT_OFFSET             3
#define  PH_MOD_SOURCE_ID_CAL_MODE_NUM_BITS               2
#define  PH_MOD_SOURCE_ID_CAL_MODE_RADIOMETRIC_CAL_MODE   0
#define  PH_MOD_SOURCE_ID_CAL_MODE_SPATIAL_CAL_MODE       1
#define  PH_MOD_SOURCE_ID_CAL_MODE_SPECTRAL_CAL_MODE      2
#define  PH_MOD_SOURCE_ID_CAL_MODE_NON_CAL_MODE           3

#define  PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BYTE_OFFSET      15
#define  PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BIT_OFFSET        6
#define  PH_MOD_SOURCE_ID_CAL_FRAME_CNT_NUM_BITS          6
#define  PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX              64

#define  PH_MOD_FPA_AEM_CONFIG_BYTE_OFFSET               16
#define  PH_MOD_FPA_AEM_CONFIG_BIT_OFFSET                 4
#define  PH_MOD_FPA_AEM_CONFIG_NUM_BITS                   1
#define  PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS              10

#define  PH_MOD_SCI_STATE_BYTE_OFFSET                    17
#define  PH_MOD_SCI_STATE_BIT_OFFSET                      6
#define  PH_MOD_SCI_STATE_NUM_BITS                        1
#define  PH_MOD_SCI_STATE_TEST                            0
#define  PH_MOD_SCI_STATE_NORMAL                          1

#define  PH_MOD_SCI_ABNORM_BYTE_OFFSET                   17
#define  PH_MOD_SCI_ABNORM_BIT_OFFSET                     7
#define  PH_MOD_SCI_ABNORM_NUM_BITS                       1
#define  PH_MOD_SCI_ABNORM_TEST                           0
#define  PH_MOD_SCI_ABNORM_NORMAL                         1


/*********************************************************/
/*  The following structure defines the entire packet    */
/*  header.                                              */
/*********************************************************/

typedef struct 
   {
   int8         version;
   int8         type;
   int8         sec_hdr_flag;
   int16        apid;
   int8         sequence_flag;
   int16        pkt_seq_cnt;
   int16        pkt_length;
   PGSt_scTime  pkt_time_code[8];
   PGSt_double  pkt_TAI_time;
   int8         QL_flag;
   int8         pkt_type;
   int8         scan_cnt;
   int8         mirror_side;
   int8         source_ID_type_flag;
   int16        earth_frame_cnt;
   int8         cal_type;
   int8         cal_mode;
   int8         cal_frame_cnt;
   int8         fpa_aem_config[PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS];
   int8         sci_state;
   int8         sci_abnorm;

   } PH_PACKET_HEADER_t;


#endif /* PH_PKT_HDR_H */
