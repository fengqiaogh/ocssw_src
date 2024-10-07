#include <time.h>
#include <timeutils.h>
#include <ctype.h>

/* ------------------------------------------------------------------- *
 * ymds2unix() - converts yr (ie 2002), mon (1-12), and day (1-31) of  *
 * month to secs since 1/1/70                                          *
 * ------------------------------------------------------------------- */
double ymds2unix(int16_t year, int16_t month, int16_t day, double secs) {
    double usec;
    struct tm trec;
    time_t secSince;

    /*                                                  */
    /* calculate seconds from 1/1/70 to input year, day */
    /* and correct to GMT                               */
    /*                                                  */
    trec.tm_year = year - 1900;
    trec.tm_mon = month - 1;
    trec.tm_mday = day;
    trec.tm_hour = 0;
    trec.tm_min = 0;
    trec.tm_sec = 0;
    trec.tm_isdst = 0;
    secSince = mktime(&trec) - gmt_offset();

    /*                                                  */
    /* Now we total the seconds                         */
    /*                                                  */
    usec = secSince + secs;

    return usec;
}
