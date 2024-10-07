#include <timeutils.h>
#include <genutils.h>

/* --------------------------------------------------------------------- */
/* unix2ymds() - converts secs since 1/1/70 to yr, mon, day, secs of day */

/* --------------------------------------------------------------------- */
void unix2ymds(double usec, int16_t *year, int16_t *mon, int16_t *day, double *secs) {

    if (usec == BAD_FLT) {
        *mon = BAD_INT;
        *day = BAD_INT;
        *secs = BAD_FLT;
    } else {
        struct tm *trec;
        time_t utime = (time_t) usec;

        trec = gmtime(&utime);

        *year = trec->tm_year + 1900;
        *day = trec->tm_mday;
        *mon = trec->tm_mon + 1;
        *secs = trec->tm_hour * 3600.
                + trec->tm_min * 60.
                + trec->tm_sec
                + fmod(usec, 1.);
    }
}

