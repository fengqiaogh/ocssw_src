#include "L1A_prototype.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "MD_metadata.h"
#include "PH_pkt_hdr.h"
#include "SC_scan.h"


void   update_scan_metadata (PH_PACKET_HEADER_t  packet_header,
                             PGSt_SMF_status     packetStatus,
                             MD_SCAN_MET_t       *scan_meta,
                             int16               *qual_value)

/*
!C************************************************************************

!Description:  This function will set the scan-level metadata fields by 
               looking at the last packet. 

!Input Parameters:
               PH_PACKET_HEADER_t  packet_header  ** Packet header         **
               PGSt_SMF_status     packetStatus   ** Results of the validation
                                                     of this packet        **

!Output Parameters:
               int16               *qual_value    ** Pixel quality value   **

!Input/Output Parameters:
               MD_SCAN_MET_t       *scan_meta     ** Scan metadata         **

Return Values: 
               None

Externally Defined:  
               MODIS_S_SUCCESS                                 (PGS_SMF.h)
               MODIS_E_CHECKSUM_NOT_VALID                      (PGS_MODIS_35005.h)
               MD_SCAN_MET_t                                   (MD_metadata.h)
               MD_SCAN_DATA_PRESENCE                           (MD_metadata.h)
               MD_BAD_CHECKSUM_PACKET                          (MD_metadata.h)
               MD_DISCARDED_PACKET                             (MD_metadata.h)
               MD_SOME_VALID_DATA_IN_SCAN                      (MD_metadata.h)
               MD_SD_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_SRCA_FRAMES_IN_SCAN                          (MD_metadata.h)
               MD_BB_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_SV_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_EV_FRAMES_IN_SCAN                            (MD_metadata.h)
               MD_DAY_SCAN                                     (MD_metadata.h)
               MD_NIGHT_SCAN                                   (MD_metadata.h)
               PH_PACKET_HEADER_t                              (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG1_GROUP                      (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_ENG2_GROUP                      (PH_pkt_hdr.h)
               PH_SEC_PKT_TYPE_DAY_GROUP                       (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH                (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE       (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE      (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE          (PH_pkt_hdr.h)
               SC_GOOD_PIXEL_PACKET                            (SC_scan.h)
               SC_BAD_CHECKSUM_PACKET                          (SC_scan.h) 
               SC_DISCARDED_PACKET                             (SC_scan.h) 
               int16                                           (hdfi.h)

Called By:
               process_a_scan

Routines Called:
               None

!Revision History:
               Revision 2.1  1997/08/25  11:50 EDT
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Original design.

               Revision 2.0  1997/07/25  10:00
               Tom Johnson/GSC (johnson@ltpmail.gsfc.nasa.gov)
               Updated for Version 2 from Version 1

               Revision 1.0  1997/06/18  16:40 EDT  
               Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
               Baseline from Version 1.

!Team-unique Header:
               This software is developed by the MODIS Science
               Data Support Team (SDST) for the National Aeronautics
               and Space Administration (NASA), Goddard Space Flight
               Center (GSFC), under contract NAS5-32373.

!References and Credits:
               None

!Design Notes: 
               This routine assumes that the packet has been validated 
               before this routine is executed.  No validation of data 
               is performed in this routine.

!END*************************************************************************
*/
{
     if (packetStatus == MODIS_S_SUCCESS) {
        *qual_value = SC_GOOD_PIXEL_PACKET;
        scan_meta->packet_scan_count = packet_header.scan_cnt;

        if ((scan_meta->ccsds_apids[0] == SC_FILL_VALUE) ||
            (scan_meta->ccsds_apids[0] == packet_header.apid))
           scan_meta->ccsds_apids[0] = packet_header.apid;
        else if ((scan_meta->ccsds_apids[1] == SC_FILL_VALUE) ||
                 (scan_meta->ccsds_apids[1] == packet_header.apid))
               scan_meta->ccsds_apids[1] = packet_header.apid;
        else if ((scan_meta->ccsds_apids[2] == SC_FILL_VALUE) ||
                 (scan_meta->ccsds_apids[2] == packet_header.apid))
               scan_meta->ccsds_apids[2] = packet_header.apid;
        
        scan_meta->packet_expedited_data_flag = packet_header.QL_flag;
        scan_meta->mirror_side = packet_header.mirror_side;
        scan_meta->scan_qual_array[MD_SCAN_DATA_PRESENCE] =
                                                    MD_SOME_VALID_DATA_IN_SCAN;
   
        if ((packet_header.pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) ||
            (packet_header.pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP)) {
           if (scan_meta->ev_start_time < 0.0)
              scan_meta->ev_start_time = packet_header.pkt_TAI_time;
        }
        else if (packet_header.source_ID_type_flag == 
                                           PH_MOD_SOURCE_ID_TYPE_FLAG_EARTH) {
           scan_meta->ev_start_time = packet_header.pkt_TAI_time;
	   if(scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN] <
	       packet_header.earth_frame_cnt)
	       scan_meta->frame_count_array[MD_EV_FRAMES_IN_SCAN] =
		   packet_header.earth_frame_cnt; 

           if (packet_header.pkt_type == PH_SEC_PKT_TYPE_DAY_GROUP)
             {
              memset(scan_meta->scan_type, '\0', sizeof(scan_meta->scan_type));
              strcpy(scan_meta->scan_type, MD_DAY_SCAN);
             }
           else
             { 
              memset(scan_meta->scan_type, '\0', sizeof(scan_meta->scan_type));
              strcpy(scan_meta->scan_type, MD_NIGHT_SCAN);
             }
        }
        else { 
           switch (packet_header.cal_type) {
              case PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE:
                 scan_meta->sd_start_time = packet_header.pkt_TAI_time;
                 if(scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN] <
		     packet_header.cal_frame_cnt)
		     scan_meta->frame_count_array[MD_SD_FRAMES_IN_SCAN] =
			 packet_header.cal_frame_cnt;
                 break;
 
              case PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE:
                 scan_meta->srca_start_time = packet_header.pkt_TAI_time; 
                 if(scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN] <
		     packet_header.cal_frame_cnt)
		     scan_meta->frame_count_array[MD_SRCA_FRAMES_IN_SCAN] =
			 packet_header.cal_frame_cnt;
                 scan_meta->srca_cal_mode = packet_header.cal_mode;
                 break;

              case PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE:
                 scan_meta->bb_start_time = packet_header.pkt_TAI_time; 
                 if(scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN] <
		     packet_header.cal_frame_cnt)
		     scan_meta->frame_count_array[MD_BB_FRAMES_IN_SCAN] =
			 packet_header.cal_frame_cnt;
                 break;

              case PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE:
                 scan_meta->sv_start_time = packet_header.pkt_TAI_time; 
                 if(scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN] <
		     packet_header.cal_frame_cnt)
		     scan_meta->frame_count_array[MD_SV_FRAMES_IN_SCAN] =
			 packet_header.cal_frame_cnt;
                 break;
           }
        }
     }
     else {
        if (packetStatus == MODIS_E_CHECKSUM_NOT_VALID) {
           *qual_value = (int16)SC_BAD_CHECKSUM_PACKET;
           scan_meta->scan_qual_array[MD_BAD_CHECKSUM_PACKET] = 
                     scan_meta->scan_qual_array[MD_BAD_CHECKSUM_PACKET] + 1;
           scan_meta->scan_qual_array[MD_DISCARDED_PACKET] = 
                     scan_meta->scan_qual_array[MD_DISCARDED_PACKET] + 1;
        }
        else {
           *qual_value = (int16)SC_DISCARDED_PACKET;
           scan_meta->scan_qual_array[MD_DISCARDED_PACKET] = 
                     scan_meta->scan_qual_array[MD_DISCARDED_PACKET] + 1;
        }
     }

}
