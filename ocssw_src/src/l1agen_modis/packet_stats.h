#ifndef PACKET_STATS
#define PACKET_STATS
/*!C-INC************************************************************************

!Description:
	Header file for packet filtering statistics.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Values:
	N/A

Externally Defined:
	stats		"print_stats.c"

Called by:
	N/A

Routines Called:
	N/A

!Revision History:
$Log: packet_stats.h,v $
Revision 6.1  2010/08/25 19:19:17  kuyper
Initial revision.

James Kuyper Jr.	James.R.Kuyper@NASA.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
extern struct packet_stats
{
    unsigned long packets;
    /* Single field validity checks */
    unsigned long version;	/* Must be 0.	*/
    unsigned long type;		/* Must be 0 (normal).	*/
    unsigned long seq_flag;	/* Must be non-zero.	*/
    unsigned long sec_hdr_flag;	/* Must be 1 (secondary header). */
    unsigned long apid;		/* 64 or 127(?) (MODIS).	*/
    unsigned long time_tag;
    unsigned long quick_look;	/* Must be 0 (not selected)	*/
    unsigned long pkt_type;	/* 0, 1, 2, 4 */

    /* Single packet consistency checks */
    unsigned long seq_type;	/* == 3 only iff pkt_type == night	*/
    unsigned long length_type;	/* short length iff pkt_type == night	*/
    unsigned long type_flag;	/* type_flag non-zero only if pkt_type == Day */
    unsigned long frame_count;	/* Consistent with type_flag, calib_type. */
    unsigned long cksum;	/* Matches checksum of data field.	*/
} stats;

void print_stats(void);
#endif
