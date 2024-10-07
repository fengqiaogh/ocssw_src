// Ported from IDL procedure to unpack and convert the CCSDS segmented time code
// from the first packet of a VIIRS packet group


//       Arguments
//
//       Name    Type    I/O     Description
//       ----    ----    ---     -----------
//       p1      byte(*)  I      Input array containing first packet
//       pyear   I*4      O      Year of packet time
//       pday    I*4      O      Day of packet time
//       eyear   I*4      O      Year of packet end time
//       eday    I*4      O      Day of packet end time
//       syear   I*4      O      Year of scan start time
//       sday    I*4      O      Day of scan start time
//       stime    R*8      O     Packet scan start time
//       ptime    R*8      O     Packet time
//       etime    R*8      O     Packet scan end time

// Liang Hong, July 22, 2015
// Liang Hong, added year, day for packet, start and end time

#include <stdint.h>
#include <timeutils.h>

int get_viirs_packet_time(uint8_t p1[], int32_t *pyear, int32_t *eyear, int32_t *syear, int32_t *pday, int32_t *eday, int32_t *sday, double *stime, double *ptime, double *etime) {
    uint8_t cctime[8];
    int toff;
    //  Get day count since Jan. 1, 1958 (Julian day 2436205)

    // Packet time
    toff = 18; // changed from 20, Liang H. 8/5/2015
    memcpy(cctime, &p1[toff], 8);
    ccsds_to_yds(cctime, pyear, pday, ptime);

    // Scan end time
    toff = 40; // changed from 38, Liang H. 8/5/2015
    memcpy(cctime, &p1[toff], 8);
    ccsds_to_yds(cctime, eyear, eday, etime);

    // Scan start time
    toff = 6;
    memcpy(cctime, &p1[toff], 8);
    ccsds_to_yds(cctime, syear, sday, stime);

    return 0;
}

