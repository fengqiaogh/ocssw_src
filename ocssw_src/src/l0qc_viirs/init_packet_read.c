//  Ported from IDL subroutine to initialize the reading of a packet file by finding
//  the first packet of a scan.

//       Arguments
//
//       Name    Type      I/O     Description
//       ----    ----      ---     -----------
//       infile    FILE        I    Input file pointer 
//       epacket byte(9318) O      Byte array containing first packet of scan
//	     len_packet	int		O		packet size
//       endfile long int        O      File Length (bytes); value set to 0 at the end of file

// Liang Hong, July 22, 2015

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "l0qc_viirs.h"

int init_packet_read(FILE *infile, uint8_t epacket[], int *len_packet, long int *endfile) {
    int apid = 0;
    uint8_t packet[30000];

    // Read to the first packet of a scan (APID = 826)
    while (apid != 826 || *len_packet != 9318) {
        read_packet(infile, packet, len_packet, endfile);
        if (*endfile == 0) {
            printf("End of file reached before start of scan\n");
            return 1;
        } else {
            apid = (packet[0] % 8)*256 + packet[1];
        }
    }

    memcpy(epacket, packet, *len_packet);

    return 0;
}



