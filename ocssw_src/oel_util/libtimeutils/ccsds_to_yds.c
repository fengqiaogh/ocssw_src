#include <stdio.h>
#include <stdint.h>
#include <timeutils.h>

int ccsds_to_yds(uint8_t *cctime,
        int32_t *iyear, int32_t *iday, double *sec) {

    // Convert CCSDS time code to year, day-of-year and seconds-of-day
    *iday = cctime[0]*256 + cctime[1];
    int32_t jday = *iday + 2436205; // Jan. 1, 1958 is Julian day 2436205
    jdate(jday, iyear, iday);

    // Get milliseconds of day count
    int32_t msec = ((cctime[2]*256 + cctime[3])*256 + cctime[4])*256 + cctime[5];

    // Get microseconds
    int32_t musec = cctime[6]*256 + cctime[7];

    *sec = msec / (double) 1000.0 + musec / (double) 1.e6;

    return 0;
}
