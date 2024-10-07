#include <stdio.h>
#include "PGS_SMF.h"
#include "packet_stats.h"
struct packet_stats stats;

void print_stats(void)
/*!C****************************************************************************

!Description:   
	Writes the packet filtering statistics to the Report Log.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Values:
	None

Externally Defined:
	N/A

Called by:
	process_a_granule()		"L1A_prototypes.h"

Routines Called:
	PGS_SMF_GenerateStatusReport()	"PGS_SMF.h"

!Revision History:
James Kuyper Jr.	James.R.Kuyper@NASA.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
{
    char buffer[1024];

    sprintf(buffer, "Packet filtering statistics.\n"
            "\tpackets:                %9lu\n"
            "\tNon-zero version:       %9lu\n"
            "\tTest:                   %9lu\n"
            "\tNon-zero sec_flag:      %9lu\n"
            "\tNon-zero sec_hdr_flag:  %9lu\n"
            "\tBad APID:               %9lu\n"
            "\tInvalid time_tag:       %9lu\n"
            "\tQuick Look:             %9lu\n"
            "\tBad Packet type:        %9lu\n"
            "\tSeq_flag/type mismatch: %9lu\n"
            "\tLength/type mismatch:   %9lu\n"
            "\tType mismatch:          %9lu\n"
            "\tinvalid frame_count:    %9lu\n",
       stats.packets, stats.version, stats.type, stats.seq_flag,
       stats.sec_hdr_flag, stats.apid, stats.time_tag, stats.quick_look,
       stats.pkt_type, stats.seq_type, stats.length_type, stats.type_flag,
       stats.frame_count);

    PGS_SMF_GenerateStatusReport(buffer);
    /* The only way this function can fail is if it can't open the log file.
     * Therefore, there's no way to report the failure of this function.
     */
}

