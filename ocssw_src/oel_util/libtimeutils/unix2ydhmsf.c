#include <timeutils.h>
#include <genutils.h>

/* -------------------------------------------------------------- */
/* unix2ydhmsf() - converts secs since 1/1/70 to YYYYDDDHHMMSSFFF */

/* -------------------------------------------------------------- */
char * unix2ydhmsf(double usec, char zone) {
    struct tm *ts;
    time_t itime;
    static char string[17];

    if(usec == BAD_FLT) {
        return("Undefined time");
    }

    itime = (time_t) usec;
    switch (zone) {
    case 'G': ts = gmtime(&itime);
        break;
    case 'L': ts = localtime(&itime);
        break;
    default:
        fprintf(stderr, "-E- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Bad timezone argument passed to ydhmsf().\n");
        exit(1);

    }
    sprintf(string, "%d%03d%02d%02d%02d%03.0f",
            ts->tm_year + 1900,
            ts->tm_yday + 1,
            ts->tm_hour,
            ts->tm_min,
            ts->tm_sec,
            floor(1000 * (usec - itime)));
    return (string);
}


