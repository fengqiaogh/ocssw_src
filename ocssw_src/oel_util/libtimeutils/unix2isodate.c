#include <timeutils.h>
#include <genutils.h>

/**
 * creates a time string "1994-11-05T13:15:30.123Z"
 * @param dtime unix time as a double
 * @param zone time zone character 'G'=gmt, 'L'=local
 * @return ISO 8601 time string (internal memory)
 */
char * unix2isodate(double dtime, char zone) {
    struct tm *ts;
    time_t itime;
    static char string[30];
    char tzone[4];

    if(dtime == BAD_FLT) {
        return "Undefined time";
    }

    itime = (time_t) dtime;
    switch (zone) {
    case 'G':
        ts = gmtime(&itime);
        break;
    case 'L':
        ts = localtime(&itime);
        break;
    default:
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Bad timezone argument passed to unix2isodate().\n");
        exit(EXIT_FAILURE);
    }
    // Add 0.0005 to correct for round off error JMG (08/03/2009)
    // Check for fracSec == 1000 JMG (02/05/2012)
    int fracSec = floor(1000 * (dtime - itime) + 0.0005);
    if (fracSec == 1000) {
        fracSec = 0;
        ts->tm_sec += 1;
    }

    sprintf(string, "%d-%02d-%02dT%02d:%02d:%02d.%03.0f", ts->tm_year + 1900, ts->tm_mon + 1,
            ts->tm_mday, ts->tm_hour, ts->tm_min, ts->tm_sec, (float) fracSec);

    switch (zone) {
    case 'G':
        strcat(string, "Z");
        break;
    case 'L':
        sprintf(tzone, "%+03d", (int) ts->tm_gmtoff / 3600);
        strcat(string, tzone);
        break;
    default:
        fprintf(stderr, "-W- %s line %d: ", __FILE__, __LINE__);
        fprintf(stderr, "Bad timezone argument passed to unix2isodate().\n");
        exit(EXIT_FAILURE);
    }

    return (string);
}

double isodate2unix(const char *isodate) {
    struct tm trec = {0};
    strptime(isodate, "%Y-%m-%dT%H:%M:%S", &trec);
    double secSince = mktime(&trec) - gmt_offset();
    char* ptr = strchr(isodate, '.');
    if (ptr) {
        double fracsecs = atof(ptr);
        secSince += fracsecs;
    }
    return secSince;
}

void isodate2ymds(const char* isodate,  int16_t *year, int16_t *month, int16_t *day, double *sec) {
    if (isodate == NULL) {
        *year = BAD_INT;
        *month = BAD_INT;
        *day = BAD_INT;
        *sec = BAD_FLT; 
        return;
    }

    struct tm trec = {0};
    if (strptime(isodate, "%Y-%m-%dT%H:%M:%S", &trec) == NULL) {
        *year = BAD_INT;
        *month = BAD_INT;
        *day = BAD_INT;
        *sec = BAD_FLT; 
    }
    else {
        *year = trec.tm_year + 1900; 
        *month = trec.tm_mon + 1; 
        *day = trec.tm_mday; 
        *sec = trec.tm_hour * 3600 + trec.tm_min * 60 + trec.tm_sec;

        // check for fraction of a second
        char* ptr = strchr(isodate, '.');
        if (ptr) {
            double fracsecs = atof(ptr);
            *sec += fracsecs;
        }

    }

}
