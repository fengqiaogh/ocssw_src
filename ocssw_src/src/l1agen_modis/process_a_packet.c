#include "L1A_prototype.h"
#include "PGS_IO_L0.h"
#include "PGS_IO.h"
#include "PGS_SMF.h"
#include "PGS_MODIS_35005.h"
#include "PH_pkt_hdr.h"

static PGSt_SMF_status		next_return;
static PH_PACKET_HEADER_t	next_packet_header;   
static PGSt_IO_L0_Packet	next_pkt[PD_PKT_BUF_MAX];
static uint16			next_packet_cont[
				    PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT
				    + PH_NUM_12BIT_WORDS_IN_HEADER + 1];
static PGSt_PC_Logical		next_packet_full = FALSE;


PGSt_SMF_status  process_a_packet (
	PGSt_IO_L0_VirtualDataSet	*L0_file,
	PGSt_IO_L0_Packet		*pkt,
	PH_PACKET_HEADER_t		*packet_header,
	uint16				*packet_cont
)

/*
!C************************************************************************

!Description:  This function is a driver that calls all the functions 
               necessary to validate a packet and unpack its contents.

!Input Parameters:
               None

!Output Parameters:
               PH_PACKET_HEADER_t  *packet_header  **  L0 data packet header **
               uint16              *packet_cont    **  Unpacked data packet  **

!Input/Output Parameters:
               PGSt_IO_L0_VirtualDataSet  *L0_file **  L0 virtual data set   **

Return Values: 
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_W_NO_MORE_PACKETS            (PGS_MODIS_35005.h)
               MODIS_F_PKT_PROCESS_FAILED         (PGS_MODIS-35005.h)
               MODIS_E_INV_VERSION                (PGS_MODIS_35005.h)
               MODIS_E_INV_TYPE                   (PGS_MODIS_35005.h)
               MODIS_E_INV_SEC_HDR_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_APID                   (PGS_MODIS_35005.h)
               MODIS_E_INV_APID_TEST              (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_SEQ_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_FAILED_TIMECODE_CONV       (PGS_MODIS_35005.h)
               MODIS_E_INV_QL_FLAG                (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TYPE               (PGS_MODIS_35005.h)
               MODIS_E_EARTH_FR_CNT_EXC_LIM       (PGS_MODIS_35005.h)
               MODIS_E_CAL_FR_CNT_EXC_LIM         (PGS_MODIS_35005.h)
               MODIS_E_CHECKSUM_NOT_VALID         (PGS_MODIS_35005.h)

Externally Defined:
               PGSt_SMF_status                 (PGS_SMF.h)
               PGSt_IO_L0_Packet               (PGS_IO.h)
               PGSt_IO_L0_VirtualDataSet       (PGS_IO.h)
               PH_PACKET_HEADER_t              (PH_pkt_hdr.h)

Called By:
               process_a_scan
               process_next_packet

Routines Called:
               read_a_packet
               unpack_packet_header   
               unpack_packet_contents
               check_checksum  

!Revision History:
$Log: process_a_packet.c,v $
Revision 5.1  2005/10/21 21:38:52  kuyper
No longer has return value indicating an invalid packet length.

Revision 4.1  2003/05/09 19:44:57  kuyper
Corrected process_a_packet() to copy over the packet and it's contents.
Change process_next_packet() to update an externally provided pointer, rather
  than copying the data it points at.

James Kuyper Jr. (kuyper@saicmodis.com)

Revision 4.0  2002/11/13 15:43:23  vlin
Updated according to process_next_packet.pdl revision 3.2

Revision 1.5  2001/11/20 17:01:41  seaton
Bit Flip Bug Fixes.

Revision 1.4  2001/04/13 18:09:20  seaton
Changed raw_mir_enc to type uint16 and fixed numerous prologs.

Revision 1.3  2000/10/03  15:27:01  seaton
Split Scan fix

Revision 2.2 2000/06/30
John Seaton GSC (seaton@ltpmail.gsfc.nasa.gov)
Code to fix split scan bug. DDTS MODx101733

!Team-unique Header:
	This software is developed by the MODIS Science Data Support Team 
	for the National Aeronautics and Space Administration, 
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits: None

Design Notes:
	The details for the packet data locations were taken from Santa Barbara
	Remote Sensing (SBRS) Contract Data Requirement List (CDRL)305, MODIS
	Engineering Telemetry Description, Tables T30-5A, T30-5B, T30-5C,
	T30-5D, and T30-5E.
!END*************************************************************************
*/
{
    PGSt_SMF_status returnStatus; /* SMF-style message returned by function */

    if (next_packet_full)
    {
	memcpy(pkt, next_pkt, sizeof(next_pkt));
	*packet_header = next_packet_header;
	memcpy(packet_cont, next_packet_cont, sizeof(next_packet_cont));
	next_packet_full = FALSE;  
	return next_return;
    }          

    if( (returnStatus = read_a_packet (L0_file, pkt)) == MODIS_S_SUCCESS &&
	(returnStatus = unpack_packet_header (pkt, packet_header))
	== MODIS_S_SUCCESS )
    {
	unpack_packet_contents (pkt, packet_header, packet_cont);
	returnStatus = check_checksum (*packet_header, packet_cont);
    }

    return returnStatus;
}

/*===========================================================================*/

PGSt_SMF_status  process_next_packet(
	PGSt_IO_L0_VirtualDataSet	*L0_file,
	PH_PACKET_HEADER_t		**packet_header
)

/*****************************************************************************
!C

!Description:  'Peeks' at next packet available in the L0_file while
                leaving it in line to be read by process_a_packet.

!Input Parameters:
       PGSt_IO_L0_VirtualDataSet  L0_file         L0 virtual data set

!Output Parameters:
       PH_PACKET_HEADER_t         packet_header   L0 data packet header

Return Values:
               MODIS_E_CAL_FR_CNT_EXC_LIM
               MODIS_E_CHECKSUM_NOT_VALID
               MODIS_E_EARTH_FR_CNT_EXC_LIM
               MODIS_E_FAILED_TIMECODE_CONV
               MODIS_E_INV_APID
               MODIS_E_INV_APID_TEST
               MODIS_E_INV_PKT_SEQ_FLAG
               MODIS_E_INV_PKT_TYPE
               MODIS_E_INV_QL_FLAG
               MODIS_E_INV_SEC_HDR_FLAG
               MODIS_E_INV_TYPE
               MODIS_F_PKT_PROCESS_FAILED
               MODIS_E_INV_VERSION
               MODIS_S_SUCCESS
               MODIS_W_CANT_PKT_PEEK
               MODIS_W_NO_MORE_PACKETS

Externally Defined:
               MODIS_E_CAL_FR_CNT_EXC_LIM         (PGS_MODIS_35005.h)
               MODIS_E_CHECKSUM_NOT_VALID         (PGS_MODIS_35005.h)
               MODIS_E_EARTH_FR_CNT_EXC_LIM       (PGS_MODIS_35005.h)
               MODIS_E_FAILED_TIMECODE_CONV       (PGS_MODIS_35005.h)
               MODIS_E_INV_APID                   (PGS_MODIS_35005.h)
               MODIS_E_INV_APID_TEST              (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_SEQ_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_PKT_TYPE               (PGS_MODIS_35005.h)
               MODIS_E_INV_QL_FLAG                (PGS_MODIS_35005.h)
               MODIS_E_INV_SEC_HDR_FLAG           (PGS_MODIS_35005.h)
               MODIS_E_INV_TYPE                   (PGS_MODIS_35005.h)
               MODIS_E_INV_VERSION                (PGS_MODIS_35005.h)
               MODIS_S_SUCCESS                    (PGS_MODIS_35005.h)
               MODIS_F_PKT_PROCESS_FAILED         (PGS_MODIS_35005.h)
               MODIS_W_CANT_PKT_PEEK              (PGS_MODIS_35005.h)
               MODIS_W_NO_MORE_PACKETS            (PGS_MODIS_35005.h)
               PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT (PD_pkt_data.h)
               PGSt_IO_L0_Packet                  (PGS_IO.h)
               PGSt_IO_L0_VirtualDataSet          (PGS_IO.h)
               PGSt_PC_Logical                    (PGS_PC.h)
               PGSt_SMF_status                    (PGS_SMF.h)
               PH_PACKET_HEADER_t                 (PH_pkt_hdr.h)
               PH_NUM_12BIT_WORDS_IN_HEADER       (PH_pkt_hder.h)
               TRUE                               (hdf.h)
               uint8                              (hdfi.h)

               File scope variables shared only with process_a_packet:
               PGSt_SMF_status          next_return
               PGSt_IO_L0_Packet        next_pkt[PD_PKT_BUF_MAX]
               PH_PACKET_HEADER_t       next_packet_header
               uint16                   next_packet_cont[PD_NUM_ELMTS_IN_DATA_FIELD_DAY_PKT
                                        + PH_NUM_12BIT_WORDS_IN_HEADER + 1]
               PGSt_PC_Logical          next_packet_full initialized to FALSE 

Called By:
               packet_of_scan
                
Routines Called:
               process_a_packet
               
!Revision History:
   Please see prologue in "process_a_packet()"

!Team-unique Header:

      This software is developed by the MODIS Science Data Support Team 
      for the National Aeronautics and Space Administration, 
      Goddard Space Flight Center, under contract NAS5-32373.
                                                               
References and Credits:   None

Design Notes:
               Issue: Globals modified by subroutines of process_a_packet
                (like incrementing global packet counters) could cause
                unintended side-effects from calling this routine.
               
!END
******************************************************************************/

{
    if (next_packet_full)
	return MODIS_E_CANT_PKT_PEEK;
    else
    {       
	next_return = process_a_packet(L0_file, next_pkt, &next_packet_header,
	    next_packet_cont);
	*packet_header = &next_packet_header;
	next_packet_full = TRUE;

	return next_return;
    }
}
