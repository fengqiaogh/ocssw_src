// Ported from IDL procedure to read a single CCSDS packet from a file.

//       Arguments
//
//       Name    Type    I/O     Description
//       ----    ----    ---     -----------
//       infile    FILE        I    Input file pointer
//       packet  byte(*)  O      Byte array containing packet
//       len     int      O      Length of packet in bytes
//       endfile long int        O      File Length (bytes); value set to 0 at the end of file
// Liang Hong, July 22, 2015
// Liang Hong, Sept. 21, 2015: changed feof check to compare pointer location with file size

#include <stdio.h>
#include <string.h>
#include <stdint.h>

int read_packet(FILE *infile, uint8_t packet[], int *len, long int *endfile) {
    // Check for end of file
    long int currpos = ftell(infile);
    if (currpos >= *endfile) {
        //Check to see if the pointer reaches end of file
        *endfile = 0;
        *len = 0;
        //printf("End of file\n");
        return 0;
    }

    // Read packet header
    uint8_t phead[6];
    fread(phead, 1, 6, infile);

    // Get length of packet body
    *len = phead[4]*256 + phead[5] + 1;
    uint8_t pbod[*len];

    fread(pbod, 1, *len, infile);

    // Concatenate header and body
    memcpy(packet, phead, 6);
    memcpy(packet + 6, pbod, *len);

    *len = *len + 6;

    return 0;
}
