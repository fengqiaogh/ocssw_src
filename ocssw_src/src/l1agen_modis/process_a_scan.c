#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_TYPES.h"
#include "PGS_MODIS_35005.h"
#include "SC_scan.h"
#include "PH_pkt_hdr.h"
#include "PD_pkt_data.h"
#include "MD_metadata.h"
#include "EN_eng_data.h"
#include "FP_failed_pkt_queue.h"
#include "L1A_prototype.h"

PGSt_SMF_status  process_a_scan (
	int                        *scan_number,
	PGSt_IO_L0_Packet          *pkt,
	PGSt_double                *scan_rate,
	PGSt_double                *scan_time,
	SC_SCAN_DATA_t             *L1A_scan,
	MD_SCAN_MET_t              *scan_meta,
	EN_VDATA_TYPE_t            *eng_data,
	FP_QUEUE_t                 *failed_pkts,
	PH_PACKET_HEADER_t         *pkt_header,
	SC_PIXEL_QUALITY_DATA_t    *scan_pixel,
	PGSt_IO_L0_VirtualDataSet  *L0_file
)
/*
!C************************************************************************

!Description: 
	This routine drives the accumulation of scans of MODIS data. It calls
	process_a_packet in a loop until it gets the first packet of a new scan.
	As each packet is received, the contents are placed into the
	appropriate part of the scan data structure, and the scan metadata and
	pixel quality data are updated.  If a data dropout is detected, then
	the routine updates the scan metadata and may return if the dropout
	extends past the end of the next scan. When the next scan is detected,
	the routine returns the scan data structure, scan metadata, pixel
	quality data, and the start time of the next scan.

!Input Parameters:
	PGSt_double              *scan_rate    ** Scan rate            **


!Output Parameters:
	SC_SCAN_DATA_t           L1A_scan      ** Scan data structure  **

	MD_SCAN_MET_t            *scan_meta    ** Scan Level Metadata  **

	FP_QUEUE_t               *failed_pkts  ** Discarded packets    **

	SC_PIXEL_QUALITY_DATA_t  *scan_pixel   ** Pixel Quality Data   **


!Input/Output Parameters:
	EN_VDATA_TYPE_t          *eng_data     ** Eng Data             **

	PGSt_double              *scan_time    ** Start time for the   **
						** scan                 **

	PGSt_IO_L0_VirtualDataSet * L0_file     ** The L0 file that     **
						** contains the packets **

	PGSt_IO_L0_Packet        *pkt          ** Packet of data       **

	PH_PACKET_HEADER_t       *pkt_header   ** Unpacked packet      **
                                                      ** header               **

	int                      *scan_number  ** Scan number          **

Return Values:
	MODIS_E_CHECKSUM_NOT_VALID
	MODIS_E_INV_PKT_TIME
	MODIS_E_SCANCNT_NOT_VALID
	MODIS_F_PKT_READ_FAILED
	MODIS_M_PKT_NOT_IN_SCAN
	MODIS_S_SUCCESS
	MODIS_W_NO_MORE_PACKETS

Externally Defined:
	EN_VDATA_TYPE_t					(EN_eng_data.h)
	FP_QUEUE_t					(FP_failed_pkt_queue.h)
	MD_SCAN_MET_t					(MD_metadata.h)
	MODIS_E_CHECKSUM_NOT_VALID			(PGS_MODIS_35005.h)
	MODIS_E_INV_PKT_TIME				(PGS_MODIS_35005.h)
	MODIS_E_SCANCNT_NOT_VALID			(PGS_MODIS_35005.h)
	MODIS_F_PKT_READ_FAILED				(PGS_MODIS-35005.h)
	MODIS_M_PKT_NOT_IN_SCAN				(PGS_MODIS_35005.h)
	MODIS_S_SUCCESS					(PGS_MODIS_35005.h)
	MODIS_W_NO_MORE_PACKETS				(PGS_MODIS-35005.h)
	PGS_FALSE					(PGS_SMF.h)
	PGS_TRUE					(PGS_SMF.h)
	PGSt_double					(PGS_TYPES.h)
	PGSt_IO_L0_Packet				(PGS_IO.h)
	PGSt_IO_L0_VirtualDataSet			(PGS_SMF.h)
	PGSt_SMF_boolean				(PGS_SMF.h)
	PGSt_SMF_status					(PGS_SMF.h)
	PH_MOD_SOURCE_ID_TYPE_FLAG_CAL			(PH_pkt_hdr.h)
	PH_PACKET_HEADER_t				(PH_pkt_hdr.h)
	PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP		(PH_pkt_hdr.h)
	PH_PRI_SEQUENCE_ONLY_PKT_IN_GROUP		(PH_pkt_hdr.h)
	PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP		(PH_pkt_hdr.h)
	PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE	(PH_pkt_hdr.h)
	PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE	(PH_pkt_hdr.h)
	PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE		(PH_pkt_hdr.h)
	PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE	(PH_pkt_hdr.h)
	PH_SEC_PKT_TYPE_DAY_GROUP			(PH_pkt_hdr.h)
	PH_SEC_PKT_TYPE_ENG1_GROUP			(PH_pkt_hdr.h)
	PH_SEC_PKT_TYPE_ENG2_GROUP			(PH_pkt_hdr.h)
	PH_SEC_PKT_TYPE_NIGHT_GROUP			(PH_pkt_hdr.h)
	SC_PIXEL_QUALITY_DATA_t				(SC_scan.h)
	SC_SCAN_DATA_t					(SC_scan.h)

Called By:
	process_a_granule

Routines Called:
	initialize_scan
	unpack_packet_contents
	check_checksum
	put_pkt_cont_in_scan
	accumulate_failed_packets
	update_pixel_qual_data
	update_scan_metadata
	process_eng_packet
	process_a_packet
	finalize_scan_metadata
	finalize_pixel_qual_data
	compute_SD_start_time
	log_fmt_msg
	packet_of_scan

!Revision History:
$Log: process_a_scan.c,v $
Revision 6.1  2014/11/13 19:05:06  jkuyper
Simplified error message to reduce log file size.

Revision 4.3  2003/04/07 17:32:11  kuyper
Fixed some comments.

Revision 4.2  2003/03/10 23:19:29  kuyper
Cleanup, standardization for ease of reading.

James Kuyper Jr. (kuyper@saicmodis.com)

Revision 2.4  2001/01/05
John Seaton  (seaton@ltpmail.gsfc.nasa.gov)
Added code to correctly count packets in a scan

Revision 2.3  2000/07/05
John Seaton/SAIC/GSC  (seaton@ltpmail.gsfc.nasa.gov)
Added code to fix split scan problem. DDTS MODx101733

Revision 2.2  1997/09/05  16:45 EDT
Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
Originated Code.

Revision 2.1  1997/08/27
Tom Johnson  (johnson@ltpmail.gsfc.nasa.gov)
Incorporate PDL walkthru comments

Revision 2.0  1997/07/24  11:55 EDT
Timi Adelekan/SAIC/GSC (adelekan@ltpmail.gsfc.nasa.gov)
Created Module for version 2 from guideline of version 1
code.

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

	This routine assumes that the data contained in the input parameters
	have been validated before this routine is executed.  No validation of
	data is performed in this routine.

!END*************************************************************************
*/

{
    /*************************************************************************/
    /*                                                                       */
    /*              Define and Initialize Local Variables                    */
    /*                                                                       */
    /*************************************************************************/
    /*                                                                       */
    /* initialize Status to MODIS_S_SUCCESS                                  */
    /*                                                                       */
    /* Set previous_scan_count to -1                                         */
    /*                                                                       */
    /*************************************************************************/

    #define FIRSTPKT 0
    #define SECONDPKT 1
    #define MAX_CALIB_FRAMES 64
    #define MAX_EARTH_FRAMES 1400


    SC_SCAN_PROC_STATE_t scan_proc_state[5];

    char   *routine = "process_a_scan";

    PGSt_SMF_status    tempStatus = MODIS_S_SUCCESS;
				   /* SMF-style message returned by function */

    static PGSt_SMF_status Status = MODIS_S_SUCCESS;
				   /* SMF-style message returned by function */

    PGSt_SMF_boolean current_scan; /* current scan being processed           */

    PGSt_double  next_scan_start_time = 0.0;
				   /* Start time scan after current scan     */
				   /* being processed                        */

    uint16    pkt_contents[PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT +
			   PH_NUM_12BIT_WORDS_IN_HEADER + 1];
				   /* unpacked packet data array             */

    int16     qual_value = 0;      /* The pixel quality value to be written  */
				   /* to the pixel quality array             */

    int16     num_packets=0;         /* Number of packets read from L0 file */

    int8      previous_scan_count = -1;

    char      msg[300];

    int       i;

    int8      pktnum;

    /*************************************************************************/
    /*                                                                       */
    /* define scan count arrays                                              */
    /*                                                                       */
    /*************************************************************************/

    int      calib_cnt_SD[2][MAX_CALIB_FRAMES];
    int      calib_cnt_SRCA[2][MAX_CALIB_FRAMES];
    int      calib_cnt_BB[2][MAX_CALIB_FRAMES];
    int      calib_cnt_SV[2][MAX_CALIB_FRAMES];
    int      earth_cnt_EV[2][MAX_EARTH_FRAMES];
    int      eng_cnt_ENG[2][2];

    /*************************************************************************/
    /*                                                                       */
    /* initialize all count arrays (eng, cal and Earth view) to zero         */
    /*                                                                       */
    /*************************************************************************/

    memset(calib_cnt_SD, 0, sizeof(calib_cnt_SD));
    memset(calib_cnt_SRCA, 0, sizeof(calib_cnt_SRCA));
    memset(calib_cnt_BB, 0, sizeof(calib_cnt_BB));
    memset(calib_cnt_SV, 0, sizeof(calib_cnt_SV));
    memset(earth_cnt_EV, 0, sizeof(earth_cnt_EV));
    memset(eng_cnt_ENG, 0, sizeof(eng_cnt_ENG));

    initialize_scan (L1A_scan, scan_pixel, scan_meta);

    scan_meta->scan_num = ++*scan_number;

    unpack_packet_contents (pkt, pkt_header, pkt_contents);

    Status = check_checksum (*pkt_header, pkt_contents);

    if (Status != MODIS_S_SUCCESS) {
	sprintf(msg, "Scan Number %d: Packet Number %d",
	    scan_meta->scan_num, num_packets);
	log_fmt_msg (MODIS_E_CHECKSUM_NOT_VALID, routine, msg);
    }

    next_scan_start_time = *scan_time + *scan_rate;
    current_scan = PGS_TRUE;

    do {
	if (Status == MODIS_S_SUCCESS)
	{
	   put_pkt_cont_in_scan (*pkt_header, pkt, pkt_contents, L1A_scan);
	   if ((pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG1_GROUP) ||
	       (pkt_header->pkt_type == PH_SEC_PKT_TYPE_ENG2_GROUP))
		process_eng_packet (eng_data, *scan_number, pkt_header, pkt);

	   /* To figure out if the packet is the first or second or only packet
	    * in the group. The default case or a '0' case is not needed. These
	    * checks happen when the packet header is verified.
	    */
	   switch (pkt_header->sequence_flag) {
	       case PH_PRI_SEQUENCE_FIRST_PKT_IN_GROUP:	/* FALLTHROUGH */
	       case PH_PRI_SEQUENCE_ONLY_PKT_IN_GROUP:
		   pktnum = FIRSTPKT;	break;
	       case PH_PRI_SEQUENCE_SECOND_PKT_IN_GROUP:
		   pktnum = SECONDPKT;	break;
	   }

	   if (pkt_header->source_ID_type_flag==PH_MOD_SOURCE_ID_TYPE_FLAG_CAL)
	   {
	      switch (pkt_header->cal_type) {
		 case PH_MOD_SOURCE_ID_CAL_TYPE_SOLAR_DIFFUSER_SOURCE:
		       calib_cnt_SD[pktnum][pkt_header->cal_frame_cnt -1] = 1;
		       break;
		 case PH_MOD_SOURCE_ID_CAL_TYPE_SRCA_CAL_SOURCE:
		       calib_cnt_SRCA[pktnum][pkt_header->cal_frame_cnt -1] = 1;
		       break;
		 case PH_MOD_SOURCE_ID_CAL_TYPE_BLACKBODY_SOURCE:
		       calib_cnt_BB[pktnum][pkt_header->cal_frame_cnt -1] = 1;
		       break;
		 case PH_MOD_SOURCE_ID_CAL_TYPE_SPACE_SOURCE:
		       calib_cnt_SV[pktnum][pkt_header->cal_frame_cnt -1] = 1;
		       break;
		 /* No default case needed, only a 2 bit field */
	      }
	   }
	   else
	   {
		switch (pkt_header->pkt_type)
		{
		    case PH_SEC_PKT_TYPE_DAY_GROUP:	/* FALLTHROUGH */
		    case PH_SEC_PKT_TYPE_NIGHT_GROUP:
			earth_cnt_EV[pktnum][pkt_header->earth_frame_cnt-1] = 1;
			break;
		    case PH_SEC_PKT_TYPE_ENG1_GROUP:
			eng_cnt_ENG[pktnum][FIRSTPKT] = 1;
			break;
		    case PH_SEC_PKT_TYPE_ENG2_GROUP:
			eng_cnt_ENG[pktnum][SECONDPKT] = 1;
			break;
		    /* No default case needed, only a 2 bit field */
		}
	    }
	}
	else
	{
	    tempStatus = accumulate_failed_packets (pkt, *failed_pkts);
	    if (tempStatus != MODIS_S_SUCCESS)
		log_fmt_msg (tempStatus, routine,
		    "Unable to add discarded packet to queue");
	}

	num_packets++;

	update_scan_metadata (*pkt_header, Status, scan_meta, &qual_value);
	update_pixel_qual_data (*pkt_header, qual_value, scan_pixel);
	Status = process_a_packet (L0_file, pkt, pkt_header, pkt_contents);

	if (Status == MODIS_E_CHECKSUM_NOT_VALID)
	{
	    sprintf(msg, "Scan Number = %d: Packet Number %d",
		scan_meta->scan_num, num_packets);
	    log_fmt_msg (MODIS_E_CHECKSUM_NOT_VALID, routine, msg);
	}

	if (Status==MODIS_W_NO_MORE_PACKETS ||
	    Status==MODIS_F_PKT_READ_FAILED)
	{
	    current_scan = PGS_FALSE;
	    log_fmt_msg (Status, routine,
		"Unable to continue packet processing");
	}
	else
	{
	    tempStatus = packet_of_scan(pkt_header, next_scan_start_time,
		&previous_scan_count, scan_proc_state, L0_file);

	    switch (tempStatus)
	    {
		case MODIS_S_SUCCESS :
		    previous_scan_count = pkt_header->scan_cnt;
		    break;

		case MODIS_M_PKT_NOT_IN_SCAN :
		    previous_scan_count = pkt_header->scan_cnt;
		    current_scan = FALSE;
		    break;

		case MODIS_E_SCANCNT_NOT_VALID :
		    if (Status != MODIS_E_CHECKSUM_NOT_VALID)
			Status = tempStatus;
		    log_fmt_msg(MODIS_E_SCANCNT_NOT_VALID, routine,
			"%d: invalid scan count:%d", num_packets,
			pkt_header->scan_cnt);
		    break;

		case MODIS_E_INV_PKT_TIME :
		    previous_scan_count = pkt_header->scan_cnt;
		    if (Status != MODIS_E_CHECKSUM_NOT_VALID)
			Status = tempStatus;
		    sprintf(msg, "\nScan Number = %d: Packet Number %d",
			scan_meta->scan_num,num_packets);
		    log_fmt_msg(MODIS_E_INV_PKT_TIME, routine, msg);
		    break;
	    }
	}
    } while (current_scan == PGS_TRUE); /*  End_while  */

    num_packets = 0;
    num_packets = eng_cnt_ENG[FIRSTPKT][0] + eng_cnt_ENG[FIRSTPKT][1] +
	eng_cnt_ENG[SECONDPKT][0] + eng_cnt_ENG[SECONDPKT][1];
    for(i=0; i<MAX_CALIB_FRAMES; i++) {
	num_packets += calib_cnt_SD[FIRSTPKT][i] +
		       calib_cnt_SD[SECONDPKT][i] +
		       calib_cnt_SRCA[FIRSTPKT][i] +
		       calib_cnt_SRCA[SECONDPKT][i] +
		       calib_cnt_BB[FIRSTPKT][i] +
		       calib_cnt_BB[SECONDPKT][i] +
		       calib_cnt_SV[FIRSTPKT][i] +
		       calib_cnt_SV[SECONDPKT][i];
    }
    for(i=0;i<MAX_EARTH_FRAMES; i++)
	num_packets += earth_cnt_EV[FIRSTPKT][i] + earth_cnt_EV[SECONDPKT][i];

    finalize_scan_metadata (scan_meta, num_packets);
    finalize_pixel_qual_data (scan_pixel, scan_meta);

    if (pkt_header->pkt_TAI_time > 0.0)
	compute_SD_start_time (pkt_header, scan_time);

    return (Status);

}  /* End process_a_scan.c */
