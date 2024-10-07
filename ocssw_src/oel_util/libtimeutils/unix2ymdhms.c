#include <timeutils.h>
#include <genutils.h>

/* --------------------------------------------------------------------- *
 * unix2ymds() - converts secs since 1/1/1970 to                         *
 * yr(2013), mon(1-12), day(1-31), hour, min, secs                       *
 * --------------------------------------------------------------------- */
void unix2ymdhms(double usec, int16_t *year, int16_t *mon, int16_t *day, int16_t *hour, int16_t *min, double *sec) {

    if (usec == BAD_FLT) {
        *year = BAD_INT;
        *mon = BAD_INT;
        *day = BAD_INT;
        *hour = BAD_INT;
        *min = BAD_INT;
        *sec = BAD_FLT;
    } else {
        struct tm *trec;
        time_t utime = (time_t) usec;

        trec = gmtime(&utime);

        *year = trec->tm_year + 1900;
        *mon = trec->tm_mon + 1;
        *day = trec->tm_mday;
        *hour = trec->tm_hour;
        *min = trec->tm_min;
        *sec = trec->tm_sec + fmod(usec, 1.);
    }
}

