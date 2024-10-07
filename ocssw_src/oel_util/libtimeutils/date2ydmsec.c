#include <timeutils.h>
#include <genutils.h>

void date2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec) {
    int yy, mm, dd, hh, mn, sc, cc;

    int count = sscanf(date, "%4d%2d%2d%2d%2d%2d%2d",
            &yy, &mm, &dd, &hh, &mn, &sc, &cc);
    if(count != 7) {
        *year = BAD_INT;
        *day = BAD_INT;
        *msec = BAD_INT;
        return;
    }

    ymdhms2ydmsec(yy, mm, dd, hh, mn, sc, year, day, msec);
    *msec = *msec + (int32_t) cc * 10;
}

void isodate2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec) {
    int yy, mm, dd, hh, mn, sc;
    int count = sscanf(date, "%4d-%2d-%2dT%2d:%2d:%2dZ",
            &yy, &mm, &dd, &hh, &mn, &sc);
    if(count != 6) {
        *year = BAD_INT;
        *day = BAD_INT;
        *msec = BAD_INT;
        return;
    }

    ymdhms2ydmsec(yy, mm, dd, hh, mn, sc, year, day, msec);
}

