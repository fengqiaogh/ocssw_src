#include "PGS_SMF.h"
#include "PGS_TD.h"
#include "PGS_IO.h"
#include "PGS_IO_L0.h"
#include "PGS_TYPES.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "PGS_MODIS_35005.h"
#include "L1A_prototype.h"

PGSt_SMF_status   unpack_MODIS_header (PGSt_IO_L0_Packet   *pkt,
                                       PH_PACKET_HEADER_t  *packet_header )

/*
!C************************************************************************

!Description:  This function extracts the MODIS header information 
               contained in the MODIS packet and places each item into the 
               appropriate variable in the PH_PACKET_HEADER_t structure.

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt             ** The MODIS packet **

!Output Parameters:
               PH_PACKET_HEADER_t  *packet_header   ** Pointer to structure 
                                                    that contains unpacked 
                                                    contents of the primary
                                                    header, secondary header, 
                                                    and MODIS header    **

Return Values: 
               MODIS_S_SUCCESS                               (PGS_MODIS_35005.h)
               MODIS_E_EARTH_FR_CNT_INVALID                  (PGS_MODIS_35005.h)
               MODIS_E_CAL_FR_CNT_INVALID                    (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_IO_L0_Packet                             (PGS_IO.h)
               PGSt_SMF_status                               (PGS_SMF.h)
               PH_PACKET_HEADER_t                            (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG1_GROUP                    (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP                    (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_BYTE_OFFSET        (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_BIT_OFFSET         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_NUM_BITS           (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH              (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BYTE_OFFSET  (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BIT_OFFSET   (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_NUM_BITS     (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_LIMIT        (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_BYTE_OFFSET         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_BIT_OFFSET          (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_NUM_BITS            (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_MODE_BYTE_OFFSET         (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_MODE_BIT_OFFSET          (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_MODE_NUM_BITS            (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BYTE_OFFSET    (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BIT_OFFSET     (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_NUM_BITS       (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX            (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS            (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_BYTE_OFFSET             (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_BIT_OFFSET              (PH_pkt_hdr.h)
               PH_MOD_FPA_AEM_CONFIG_NUM_BITS                (PH_pkt_hdr.h)
               PH_MOD_SCI_STATE_BYTE_OFFSET                  (PH_pkt_hdr.h)
               PH_MOD_SCI_STATE_BIT_OFFSET                   (PH_pkt_hdr.h)
               PH_MOD_SCI_STATE_NUM_BITS                     (PH_pkt_hdr.h)
               PH_MOD_SCI_ABNORM_BYTE_OFFSET                 (PH_pkt_hdr.h)
               PH_MOD_SCI_ABNORM_BIT_OFFSET                  (PH_pkt_hdr.h)
               PH_MOD_SCI_ABNORM_NUM_BITS                    (PH_pkt_hdr.h)
               PD_NUM_BITS_IN_BYTE                           (PD_pkt_data.h)

Called By:
               unpack_packet_header

Routines Called:
               extr_bits                       
               log_fmt_msg

!Revision History:
  James Kuyper Jr: James.R.Kuyper@nasa.gov

  $Log: unpack_MODIS_header.c,v $
  Revision 6.4  2015/08/22 20:36:00  jkuyper
  Changed to treat frame counts less than 1 as invalid.

  Revision 6.3  2014/10/22 18:33:26  jkuyper
  Corrected second error message the same way as the first one.

  Revision 6.2  2014/10/20 21:27:46  jkuyper
  Resolve Bug 5570 by moving the line that sets the earth_frame_count to a fill value.

  Revision 6.1  2010/08/25 18:26:28  kuyper
  Changed to always extract the source_ID_type_flag, even for engineering data
    packets.

  Revision 5.1  2005/12/30 19:30:32  vlin
  validate frame counts using 1354


!Team-unique Header:
               This software is developed by the MODIS Science 
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

Design Notes: 
               The details for the packet header data locations were taken
               from Hughes Santa Barbara Remote Sensing (SBRS) Contract Data
               Requirements List (CDRL) 305, MODIS Engineering Telemetry 
               Description, Figures 3-7 and 30-8.


!END************************************************************************
*/


{

  PGSt_SMF_status  returnStatus;   /* SMF-style message returned by function */
  char             *routine = "unpack_MODIS_header";
  char             msg[300];
  int              i;


  returnStatus = MODIS_S_SUCCESS;

  packet_header->source_ID_type_flag = extr_bits (pkt,
      PH_MOD_SOURCE_ID_TYPE_FLAG_BIT_OFFSET,
      PH_MOD_SOURCE_ID_TYPE_FLAG_BYTE_OFFSET,
      PH_MOD_SOURCE_ID_TYPE_FLAG_NUM_BITS);

  if ((packet_header->pkt_type != PH_SEC_PKT_TYPE_ENG1_GROUP) &&
     (packet_header->pkt_type != PH_SEC_PKT_TYPE_ENG2_GROUP))
    {
     if (packet_header->source_ID_type_flag == PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH)
       {
        packet_header->earth_frame_cnt = extr_bits (pkt,
           PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BIT_OFFSET,
           PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_BYTE_OFFSET,
           PH_MOD_SOURCE_ID_EARTH_FRAME_CNT_NUM_BITS);

        if (packet_header->earth_frame_cnt < 1 ||
	    packet_header->earth_frame_cnt > PH_REALISTIC_NUM_EV_FRAMES)
          {
           returnStatus = MODIS_E_EARTH_FR_CNT_INVALID;
           sprintf(msg, "Earth Frame Count: %d",
	       packet_header->earth_frame_cnt);
           log_fmt_msg (MODIS_E_EARTH_FR_CNT_INVALID, routine, msg);

           packet_header->earth_frame_cnt = -1;
          }

       }

     else
       {
        packet_header->cal_type = extr_bits (pkt, PH_MOD_SOURCE_ID_CAL_TYPE_BIT_OFFSET,
           PH_MOD_SOURCE_ID_CAL_TYPE_BYTE_OFFSET, PH_MOD_SOURCE_ID_CAL_TYPE_NUM_BITS);

        packet_header->cal_mode = extr_bits (pkt, PH_MOD_SOURCE_ID_CAL_MODE_BIT_OFFSET,
           PH_MOD_SOURCE_ID_CAL_MODE_BYTE_OFFSET, PH_MOD_SOURCE_ID_CAL_MODE_NUM_BITS);

        packet_header->cal_frame_cnt = extr_bits (pkt,
           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BIT_OFFSET,
           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_BYTE_OFFSET,
           PH_MOD_SOURCE_ID_CAL_FRAME_CNT_NUM_BITS);

        if (packet_header->cal_frame_cnt < 1 ||
	    packet_header->cal_frame_cnt > PH_MOD_SOURCE_ID_CAL_FRAME_CNT_MAX)
          {
           returnStatus = MODIS_E_CAL_FR_CNT_INVALID;
           sprintf(msg, "Calibration Frame Count: %d",
	       packet_header->cal_frame_cnt);
           log_fmt_msg (MODIS_E_CAL_FR_CNT_INVALID, routine, msg);
           packet_header->cal_frame_cnt = -1;
          }
       }
    }

  for ( i = 0; i < PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS; i++ )
     packet_header->fpa_aem_config[i] = extr_bits(pkt,
        ((PH_MOD_FPA_AEM_CONFIG_BIT_OFFSET + i) % PD_NUM_BITS_IN_BYTE),
        (PH_MOD_FPA_AEM_CONFIG_BYTE_OFFSET +
        ((PH_MOD_FPA_AEM_CONFIG_BIT_OFFSET + i) / PD_NUM_BITS_IN_BYTE)),
        PH_MOD_FPA_AEM_CONFIG_NUM_BITS);

  packet_header->sci_state = extr_bits (pkt, PH_MOD_SCI_STATE_BIT_OFFSET,
     PH_MOD_SCI_STATE_BYTE_OFFSET, PH_MOD_SCI_STATE_NUM_BITS);

  packet_header->sci_abnorm = extr_bits (pkt, PH_MOD_SCI_ABNORM_BIT_OFFSET,
     PH_MOD_SCI_ABNORM_BYTE_OFFSET, PH_MOD_SCI_ABNORM_NUM_BITS);

  return (returnStatus);

}
