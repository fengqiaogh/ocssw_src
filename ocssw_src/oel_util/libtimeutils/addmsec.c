#include <timeutils.h>

void addmsec(int32_t *year, int32_t *day, int32_t *msec, int32_t delta) {
    int16_t year1, day1;
    double sec1, unixTime;

    unixTime = yds2unix(*year, *day, *msec / 1000.0) + delta / 1000.0;
    unix2yds(unixTime, &year1, &day1, &sec1);
    *year = year1;
    *day = day1;
    *msec = sec1 * 1000.0;
}



