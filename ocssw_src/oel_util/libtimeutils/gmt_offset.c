#include <timeutils.h>


/* -------------------------------------------------------------- */
/* gmt_offset) - determines offset from local to GMT time n secs  */

/* -------------------------------------------------------------- */
time_t gmt_offset(void) {
    struct tm trec;
    trec.tm_year = 70;
    trec.tm_mon = 0;
    trec.tm_mday = 1;
    trec.tm_hour = 0;
    trec.tm_min = 0;
    trec.tm_sec = 0;
    trec.tm_isdst = 0;

    return (mktime(&trec));
}

