#include "PGS_MODIS_35005.h"
#include "L1A_prototype.h"


PGSt_SMF_status packet_of_scan ( PH_PACKET_HEADER_t *pkt_header,
                                 PGSt_double next_scan_start_time,
				 int8 		*previous_scan_count,
                                 SC_SCAN_PROC_STATE_t *scan_proc_state,
				 PGSt_IO_L0_VirtualDataSet  *L0_file )

/*******************************************************************************
!C

!Description:  This routine determines, from the pkt_header, whether a
		packet belongs in the current scan, the next scan, or is
		a misplaced packet because of a corrupted time stamp field.

!Input Parameters: 
       PH_PACKET_HEADER_t   pkt_header             Unpacked packet header
       PGSt_double          next_scan_start_time   Estimated (maximum) start 
                                                   time of the next scan
       PGSt_IO_L0_VirtualDataSet  L0_file          The L0 file that contains the 
                                                   packets

!Output Parameters:         None

!Input/Output Parameters: 
       int8		    previous_scan_count    scan count field of last
				                   valid packet
       SC_SCAN_PROC_STATE_t scan_proc_state        processing state of scan used 
                                                   to identify misplaced packets

Return Values:       
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_M_PKT_NOT_IN_SCAN		  (PGS_MODIS_35005.h)
               MODIS_E_SCANCNT_NOT_VALID          (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TIME		  (PGS_MODIS_35005.h)

Externally Defined:      
               int8                               (hdfi.h)
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_M_PKT_NOT_IN_SCAN            (PGS_MODIS_35005.h)
               MODIS_E_SCANCNT_NOT_VALID          (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TIME               (PGS_MODIS_35005.h)
               PGSt_double                        (PGS_TYPES.h)
               PGSt_IO_L0_VirtualDataSet          (PGS_SMF.h)
               PGSt_SMF_status                    (PGS_SMF.h)
               PH_MOD_SOURCE_ID_TYPE_FLAG_CAL     (PH_pkt_hdr.h)
               PH_MOD_SOURCE_ID_CAL_TYPE_NUM_BITS (PH_pkt_hdr.h)
               PH_PACKET_HEADER_t                 (PH_pkt_hdr.h)
               SC_FILL_VALUE                      (SC_scan.h)
               SC_SCAN_PROC_STATE_t               (SC_scan.h)
               SC_SCAN_RATE_TOLERANCE             (SC_scan.h)
               uint8                              (hdfi.h)

Called By:
               process_a_scan

Routines Called:
	       process_next_packet


!Revision History: 
     $Log: packet_of_scan.c,v $
     Revision 5.1  2007/01/26 22:10:33  kuyper
     Corrected to partially resolve bug 484 by changing the logic for detecting
       misplaced packets and packets with corrupted scan_cnt values.
     A full resolution will require major re-design.

     Revision 4.1  2003/05/09 18:50:04  kuyper
     Changed to recieve next_pkt_header as a pointer to the look-ahead buffer,
       rather than as a copy of it.

     Revision 4.0  2002/11/22 16:13:46  vlin
     Updated according to packet_of_scan.pdl revision 4.0
     vlin@saicmodis.com

     Revision 1.3  2001/04/13 18:09:20  seaton
     Changed raw_mir_enc to type uint16 and fixed numerous prologs.

     Revision 1.2  2000/10/03  15:27:01  seaton
     Split Scan fix
     
!Team-unique Header:

       This software is developed by the MODIS Science Data Support Team 
       for the National Aeronautics and Space Administration, 
       Goddard Space Flight Center, under contract NAS5-32373.

References and Credits:

Design Notes: 

!END
*****************************************************************************/

{

    PGSt_SMF_status            Status, returnStatus = MODIS_S_SUCCESS;
    int                        i;
    uint8                      hash_key;
    PH_PACKET_HEADER_t         *next_pkt_header=NULL;

    enum {NOPACKETPROCESSED, NOPACKETSHARETIME, MORETHANONEPACKETSHARETIME};

    if (pkt_header->source_ID_type_flag == PH_MOD_SOURCE_ID_TYPE_FLAG_CAL) 
	hash_key = pkt_header->cal_type;
    else /* Earth view sector or engineering packets */
	hash_key = (PH_MOD_SOURCE_ID_CAL_TYPE_NUM_BITS << 1);

    if (*previous_scan_count == SC_FILL_VALUE)
    {
	*previous_scan_count = pkt_header->scan_cnt;
	for (i=0; i<5; i++)
	{
	    scan_proc_state[i].pkt_TAI_time = 0.0;
	    scan_proc_state[i].sector_packet_count = 0;
	}
    }

    if ((pkt_header->pkt_TAI_time + SC_SCAN_RATE_TOLERANCE) >=
	next_scan_start_time)
    {    /* packet may have an invalid time for this scan */
	returnStatus = MODIS_M_PKT_NOT_IN_SCAN;
	*previous_scan_count = pkt_header->scan_cnt;
	for (i=0; i<5; i++)
	{
	    scan_proc_state[i].pkt_TAI_time = 0.0;
	    scan_proc_state[i].sector_packet_count = 0;
	}
	scan_proc_state[hash_key].pkt_TAI_time = pkt_header->pkt_TAI_time;
	scan_proc_state[hash_key].sector_packet_count = NOPACKETSHARETIME;
    }
    else if (pkt_header->scan_cnt != *previous_scan_count)
    {
	/* packet may have an invalid scan count for this scan */
	Status = process_next_packet(L0_file, &next_pkt_header);

	if ((Status != MODIS_S_SUCCESS) ||
	    (((*previous_scan_count+1)%8) == pkt_header->scan_cnt &&
	    (next_pkt_header->pkt_TAI_time >= next_scan_start_time ||
	     pkt_header->scan_cnt == next_pkt_header->scan_cnt)))
	{
	    /* Nothing in next packet refutes the assumption that the current
	     * packet belongs to the next scan
	     */
	    returnStatus = MODIS_M_PKT_NOT_IN_SCAN;
	    *previous_scan_count = pkt_header->scan_cnt;
	    for (i=0; i<5; i++) {
		scan_proc_state[i].pkt_TAI_time = 0.0;
		scan_proc_state[i].sector_packet_count = 0;
	    }
	    scan_proc_state[hash_key].pkt_TAI_time = pkt_header->pkt_TAI_time;
	    scan_proc_state[hash_key].sector_packet_count = NOPACKETSHARETIME;
	}
	else /* found a misplaced packet or a packet with a corrupted scan_cnt*/
	    returnStatus = MODIS_E_SCANCNT_NOT_VALID;
    }
    else switch (scan_proc_state[hash_key].sector_packet_count)
    {	/* packet appears to be part of this scan, but may be a rogue packet */
    case NOPACKETPROCESSED :  
	/* no packets processed for sector */
	scan_proc_state[hash_key].pkt_TAI_time = pkt_header->pkt_TAI_time;
	scan_proc_state[hash_key].sector_packet_count = NOPACKETSHARETIME;
	break;
    case NOPACKETSHARETIME :  
	/* no processed packets for sector share the same time */
	if (scan_proc_state[hash_key].pkt_TAI_time == pkt_header->pkt_TAI_time) 
	   scan_proc_state[hash_key].sector_packet_count = 
		  MORETHANONEPACKETSHARETIME;
	else
	   scan_proc_state[hash_key].pkt_TAI_time=pkt_header->pkt_TAI_time;
	break;
    case MORETHANONEPACKETSHARETIME : 
	/*more than one packet for sector processed with the same time */
	if (scan_proc_state[hash_key].pkt_TAI_time != pkt_header->pkt_TAI_time)
	    returnStatus = MODIS_E_INV_PKT_TIME;
	break;
    }

    return returnStatus;
}

