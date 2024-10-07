// Ported from IDL procedure to read all of the packets for a VIIRS scan 
//  into a packet buffer

//       Arguments
//
//       Name    Type      I/O     Description
//       ----    ----      ---     -----------
//       infile    FILE        I    Input file pointer 
//       epacket byte(9318) I/O    First packet already read from
//                                  file
//	     len_packet	int		O		packet size
//       pbuffer byte(32000,513) O Packet buffer array
//       npkts   int        O      Number of packets stored in buffer
//       endfile long int        O      File Length (bytes); value set to 0 at the end of file
// Liang Hong, July 22, 2015

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <timeutils.h>
#include "l0qc_viirs.h"

int read_viirs_scan_packets(FILE *infile, uint8_t epacket[], int *len_packet, uint8_t pbuffer[], int *npkts, long int *endfile) {
    // Get VIIRS scan number and start time from first packet

    int32_t iyear, iday;
    int i;
    double stime;
    uint8_t cctime[8];
    memcpy(cctime, &epacket[6], 8);
    ccsds_to_yds(cctime, &iyear, &iday, &stime);
    int32_t jd1 = jday(iyear, 1, iday);

    long scnum;
    if (epacket[1] == 58) {
        //scnum = swap_endian(long(epacket(48:51),0));
        scnum = (epacket[48] << 24) + (epacket[49] << 16) + (epacket[50] << 8) + epacket[51];
    } else {
        //scnum = swap_endian(long(epacket(34:37),0));
        scnum = (epacket[34] << 24) + (epacket[35] << 16) + (epacket[36] << 8) + epacket[37];
    }

    // Store the engineering packet in the buffer
    for (i = 0; i<*len_packet; i++) pbuffer[i + 0 * 32000] = epacket[i];
    *npkts = 1;

    // Read all of the packets with this start time and scan number 
    //  and store in the packet buffer

    // Read the next packet
    uint8_t packet[30000];
    read_packet(infile, packet, len_packet, endfile);
    if (*endfile == 0) return 0;

    // Check for missing first packet of a group
    int first = (packet[0] & 8) / 8;
    while (!first) {
        for (i = 0; i<*len_packet; i++) pbuffer[i + *npkts * 32000] = packet[i];
        *npkts = *npkts + 1;
        read_packet(infile, packet, len_packet, endfile);
        if (*endfile == 0) return 0;
        first = (packet[0] & 8) / 8;
    }

    // Check scan number and start time

    //scn = swap_endian(long(packet(34:37),0))
    long scn;
    scn = (packet[34] << 24) + (packet[35] << 16) + (packet[36] << 8) + packet[37];
    memcpy(cctime, &packet[6], 8);
    double stm;
    ccsds_to_yds(cctime, &iyear, &iday, &stm);
    stm = stm + 864 * (jday(iyear, 1, iday) - jd1);

    int apd, apid, npkg;
    double usec, sec;
    int16_t iyr, ih, mn, mm, dd;
    while ((scn <= scnum) && (stm <= stime)) {
        // Check for packet out of order
        if ((scn != scnum) || (stm != stime)) {
            usec = yds2unix(iyear, iday, stm);
            unix2ymdhms(usec, &iyr, &mm, &dd, &ih, &mn, &sec);
            printf("Packet out of order at %02d:%02d:%02d\n", ih, mn, (int) sec);
            // Find next packet group
            first = 0;
            while (!first) {
                read_packet(infile, packet, len_packet, endfile);
                if (*endfile == 0) return 0;
                first = (packet[0] & 8) / 8;
            }
        } else {
            // Get APID and number of packets in group and load packet into buffer
            apid = (packet[0] % 8)*256 + packet[1];
            npkg = packet[14];
            for (i = 0; i<*len_packet; i++) pbuffer[i + *npkts * 32000] = packet[i];
            *npkts = *npkts + 1;

            // Read remaining packets in group
            int ipkg = 0;
            while (ipkg <= npkg) {
                read_packet(infile, packet, len_packet, endfile);
                if (*endfile == 0) return 0;
                apd = (packet[0] % 8)*256 + packet[1];

                //  Check APID// if not a match, missing packets
                if (apd == apid) {
                    for (i = 0; i<*len_packet; i++) pbuffer[i + *npkts * 32000] = packet[i];
                    *npkts = *npkts + 1;
                    ipkg++;
                } else {
                    ipkg = npkg + 1;
                }
            }

            //  Check for another first packet of a group
            first = (packet[0] & 8) / 8;
            while (!first) {
                for (i = 0; i<*len_packet; i++) pbuffer[i + *npkts * 32000] = packet[i];
                *npkts = *npkts + 1;
                read_packet(infile, packet, len_packet, endfile);
                if (*endfile == 0) return 0;
                first = (packet[0] & 8) / 8;
            }
        }
        scn = (packet[34] << 24) + (packet[35] << 16) + (packet[36] << 8) + packet[37];

        memcpy(cctime, &packet[6], 8);
        ccsds_to_yds(cctime, &iyear, &iday, &stm);

        stm = stm + 864 * (jday(iyear, 1, iday) - jd1);
    }

    // Load last packet read into engineering packet for next scan
    memcpy(epacket, &packet, *len_packet);

    return 0;
}







