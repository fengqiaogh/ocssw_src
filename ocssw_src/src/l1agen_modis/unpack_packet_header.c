#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "PH_pkt_hdr.h"
#include "L1A_prototype.h"
#include "packet_stats.h"


PGSt_SMF_status   unpack_packet_header ( PGSt_IO_L0_Packet   *pkt,
                                         PH_PACKET_HEADER_t  *packet_header )

/*
!C********************************************************************************

!Description:  This function calls other routines to extract and validate
               the packet's primary, secondary, and MODIS header. 

!Input Parameters:
               PGSt_IO_L0_Packet   *pkt             ** The MODIS packet         **

!Output Parameters:
               PH_PACKET_HEADER_t  *packet_header   ** Unpacked primary, secondary
                                                       and MODIS header         **

Return Values: 
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_E_INV_VERSION                (PGS_MODIS_35005.h)
               MODIS_E_INV_TYPE                   (PGS_MODIS_35005.h)
               MODIS_E_INV_SEC_HDR_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_APID                   (PGS_MODIS_35005.h)
               MODIS_E_INV_APID_TEST              (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_SEQ_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_PKTLEN                 (PGS_MODIS_35005.h)
               MODIS_E_FAILED_TIMECODE_CONV       (PGS_MODIS_35005.h)
               MODIS_E_INV_QL_FLAG                (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TYPE               (PGS_MODIS_35005.h)
               MODIS_E_EARTH_FR_CNT_EXC_LIM       (PGS_MODIS_35005.h)
               MODIS_E_CAL_FR_CNT_EXC_LIM         (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_IO_L0_Packet                  (PGS_IO.h)
               PH_PACKET_HEADER_t                 (PH_pkt_hdr.h)
	       PH_PRI_LONG_PKT_LENGTH		  (PH_pkt_hdr.h)
	       PH_PRI_PKT_LENGTH_BYTE_OFFSET	  (PH_pkt_hdr.h)
	       PH_PRI_SHORT_PKT_LENGTH		  (PH_pkt_hdr.h)
	       PH_SEC_PKT_TYPE_DAY_GROUP	  (PH_pkt_hdr.h)
	       PH_SEC_PKT_TYPE_ENG1_GROUP	  (PH_pkt_hdr.h)
	       PH_SEC_PKT_TYPE_ENG2_GROUP	  (PH_pkt_hdr.h)
	       PH_SEC_PKT_TYPE_NIGHT_GROUP	  (PH_pkt_hdr.h)

Called By:
               load_eng_data
               process_a_packet

Routines Called:
               unpack_primary_header
               unpack_secondary_header
               unpack_MODIS_header

!Revision History:
  $Log: unpack_packet_header.c,v $
  Revision 6.1  2010/08/25 19:10:41  kuyper
  Changed to return MODIS_E_L1A for any failure detected in subroutines.
  Changed from correcting packet lengths to be consistent with the packet type,
    to marking the packet as invalid when they are inconsistent.
  Added two other intra-packet consistency tests.
  Added collection of packet filtering statistics.

  Revision 5.1  2005/10/21 16:12:15  kuyper
  Changed to correct invalid packet lengths based upon the packet type.

  Revision 4.2  2002/12/03 20:47:54  vlin
  Updated after code walk through.

  Revision 4.1  2002/10/16 16:10:21  vlin
  call memcmp before calling memset.

               Revision 2.0  1997/07/15  16:07
               Tom Johnson/SAIC/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Completely overhauled for Version 2

               Revision 1.0  1997/06/18  16:40 EDT
               Timi Adelekan/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:

       This software is developed by the MODIS Science Data Support Team 
       for the National Aeronautics and Space Administration, 
       Goddard Space Flight Center, under contract NAS5-32373.

Design Notes: 
               The details for the packet header data locations were taken
               from Hughes Santa Barbara Remote Sensing (SBRS) Contract Data
               Requirements List (CDRL) 305, MODIS Engineering Telemetry 
               Description, Figures 3-7 and 30-8.


!END
*************************************************************************/

{
  int              i;
  PGSt_SMF_status  returnStatus = MODIS_S_SUCCESS, tempStatus; 
  char	*routine = "unpack_packet_header";
  char	msg[256];

  packet_header->fpa_aem_config[0] = -1;
  memset(packet_header->fpa_aem_config+1, -1, 
         sizeof(packet_header->fpa_aem_config[0]));

  if (memcmp(packet_header->fpa_aem_config, packet_header->fpa_aem_config+1,
      sizeof(packet_header->fpa_aem_config[0])) ) { /* memset didn't work */
      for (i=1; i<PH_MOD_FPA_AEM_CONFIG_NUM_ELEMENTS; i++)
           packet_header->fpa_aem_config[i] = -1; 
  }
  else
      memset(packet_header->fpa_aem_config+2, -1, 
             sizeof(packet_header->fpa_aem_config));

  packet_header->version = -1;
  packet_header->type = -1;
  packet_header->sec_hdr_flag = -1;
  packet_header->apid = -1;
  packet_header->sequence_flag = -1;
  packet_header->pkt_seq_cnt = -1;
  packet_header->pkt_length = -1;
  packet_header->pkt_TAI_time = -1.0;
  packet_header->QL_flag = -1;
  packet_header->pkt_type = -1;
  packet_header->scan_cnt = -1;
  packet_header->mirror_side = -1;
  packet_header->source_ID_type_flag = -1;
  packet_header->earth_frame_cnt = -1;
  packet_header->cal_type = -1;
  packet_header->cal_mode = -1;
  packet_header->cal_frame_cnt = -1;
  packet_header->sci_state = -1;
  packet_header->sci_abnorm = -1;
  memset(packet_header->pkt_time_code, '\0', 
         sizeof(packet_header->pkt_time_code));

  stats.packets++;
  tempStatus = unpack_primary_header(pkt, packet_header);
  if (tempStatus != MODIS_S_SUCCESS)
  {
      switch(tempStatus)
      {
      case MODIS_E_INV_VERSION:		stats.version++;	break;
      case MODIS_E_INV_TYPE:		stats.type++;		break;
      case MODIS_E_INV_SEC_HDR_FLAG:	stats.sec_hdr_flag++;	break;
      case MODIS_E_INV_APID:		/*FALLTHROUGH*/	
      case MODIS_E_INV_APID_TEST:	stats.apid++;		break;
      case MODIS_E_INV_PKT_SEQ_FLAG:	stats.seq_flag++;	break;
      default:
	  log_fmt_msg(MODIS_E_L1A, routine,
	      "unpack_primary_header() returned %ld", (long)tempStatus);
      }

      returnStatus = MODIS_E_L1A;
  }

  tempStatus = unpack_secondary_header(pkt, packet_header);
  if (tempStatus != MODIS_S_SUCCESS)
  {
     switch(tempStatus)
     {
     case MODIS_E_FAILED_TIMECODE_CONV:	stats.time_tag++;	break;
     case MODIS_E_INV_QL_FLAG:		stats.quick_look++;	break;
     case MODIS_E_INV_PKT_TYPE:		stats.pkt_type++;	break;
     default:
	  log_fmt_msg(MODIS_E_L1A, routine,
	      "unpack_secondary_header() returned %ld", (long)tempStatus);
     }
     returnStatus = MODIS_E_L1A;
  }

  if (returnStatus == MODIS_S_SUCCESS)
  {
      if(packet_header->pkt_length !=
	  (packet_header->pkt_type == PH_SEC_PKT_TYPE_NIGHT_GROUP
	  ?  PH_PRI_SHORT_PKT_LENGTH : PH_PRI_LONG_PKT_LENGTH))
      {

	  sprintf(msg, "packet length value: %d Valid Values: long packet %d "
	      "short packet %d", packet_header->pkt_length,
	      PH_PRI_LONG_PKT_LENGTH, PH_PRI_SHORT_PKT_LENGTH);
	  log_fmt_msg(MODIS_W_INV_PKTLEN, routine, msg);

	  stats.length_type++;
	  returnStatus = MODIS_W_INV_PKTLEN;
      }
      else if(
	  (packet_header->sequence_flag == PH_PRI_SEQUENCE_ONLY_PKT_IN_GROUP)
	  != (packet_header->pkt_type == PH_SEC_PKT_TYPE_NIGHT_GROUP))
      {
	  stats.seq_type++;
	  returnStatus = MODIS_E_INV_PKT_SEQ_FLAG;
      }
  }

  tempStatus = unpack_MODIS_header(pkt, packet_header);
  if (tempStatus != MODIS_S_SUCCESS)
  {
     stats.frame_count++;	/* Only failure mode is excess frame count. */
     returnStatus = MODIS_E_L1A;
  }

  if(returnStatus == MODIS_S_SUCCESS &&
      packet_header->source_ID_type_flag != PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH 
      && packet_header->pkt_type != PH_SEC_PKT_TYPE_DAY_GROUP)
  {
      stats.type_flag++;
      returnStatus = MODIS_E_INV_SECTOR;
  }

  return returnStatus;
} 

