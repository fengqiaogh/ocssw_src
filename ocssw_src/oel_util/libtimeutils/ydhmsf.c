#include <timeutils.h>
#include <genutils.h>

/************************************************************************
This function converts its input arguments to a string "YYYYDDDHHMMSSFFF"
where YYYY is the year, DDD is the day of the year, HH is the hour,
MM is the minute, SS is the second, and FFF is the fraction of a second.
The first argument represents the number of seconds elapsed since
1-Jan-1970 00:00:00.000 GMT.  The second argument determines whether
the output string represents local time (L) or GMT (G).
 ************************************************************************/
char * ydhmsf(double dtime, char zone) {
    struct tm *ts;
    time_t itime;
    static char string[17];

    if(dtime == BAD_FLT) {
        return "Undefined time";
    }

    itime = (time_t) dtime;
    switch (zone) {
    case 'G': ts = gmtime(&itime);
        break;
    case 'L': ts = localtime(&itime);
        break;
    default:
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Bad timezone argument passed to ydhmsf().\n");
        exit(EXIT_FAILURE);
    }
    // Add 0.0005 to correct for round off error JMG (08/03/2009)
    // Check for fracSec == 1000 JMG (02/05/2012)
    int fracSec = floor(1000 * (dtime - itime) + 0.0005);
    if (fracSec == 1000) {
        fracSec = 0;
        ts->tm_sec += 1;
    }

    sprintf(string, "%d%03d%02d%02d%02d%03.0f",
            ts->tm_year + 1900,
            ts->tm_yday + 1,
            ts->tm_hour,
            ts->tm_min,
            ts->tm_sec,
            (float) fracSec);
    return (string);
}
