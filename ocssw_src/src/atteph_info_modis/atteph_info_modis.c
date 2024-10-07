#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> 
#include <PGS_EPH.h>
#include <stddef.h>     /* for `size_t' */
#include <timeutils.h>

/* prototype for function used to swap byte order */
extern int byteswap(char *, int);

#define RTRN_ERROR   255
#define RTRN_SUCCEED   0

int getEphemHeaders(char *ephfile) {
    FILE *ephemFilePtr;
    PGSt_ephemHeader header;
    char *ptr, *ptr1, *ptr2, *ptr3, *ptr4;
    int ll, ll1, n, m, q;
    int16_t year, month, day, hour, minute, second;
    double secs;
    double secTAI93 = 725846400.;

    static char buf[255];
    static char shortname[8];
    size_t len = 8;

    strcpy(buf, ephfile);
    strncpy(shortname, buf, len);
    fprintf(stdout, "shortname=%s\n", shortname);
    fprintf(stdout, "localgranuleid=%s\n", buf);

    if ((ephemFilePtr = fopen(buf, "r")) == NULL) {
        fprintf(stderr, "Error opening Ephemeris file\n");
        return RTRN_ERROR;
    } else {
        fread(&header, sizeof (PGSt_ephemHeader), 1, ephemFilePtr);
        fclose(ephemFilePtr);
    }

    ptr = (char *) (&header.startTime);
    ll = sizeof (PGSt_double);
    byteswap(ptr, ll);

    ptr1 = (char *) (&header.endTime);
    ll1 = sizeof (PGSt_double);
    byteswap(ptr1, ll1);

    ptr2 = (char *) (&header.interval);
    for (n = 0; n < 6; n++, ptr2 += 4) byteswap(ptr2, 4);

    ptr3 = (char *) (&header.keplerElements);
    for (m = 0; m < 7; m++, ptr3 += 8) byteswap(ptr3, 8);

    ptr4 = (char *) (&header.qaParameters);
    for (q = 0; q < 20; q++, ptr4 += 4) byteswap(ptr4, 4);

    unix2ymds(secTAI93 + header.startTime, &year, &month, &day, &secs);

    fprintf(stdout, "rangebeginningdate=%4d-%02d-%02d\n", year, month, day);

    hour = (int) (24. * (secs / 86400.));
    minute = (int) ((secs - 3600. * (double) hour) / 60.);
    second = (int) ((secs - 3600. * (double) hour - (double) minute * 60.));
    double fsec = 1000000. * (secs - 3600. * (double) hour - (double) minute * 60. -
            (double) second);

    fprintf(stdout, "rangebeginningtime=%02d:%02d:%02d.%d\n", hour, minute, second, (int) fsec);

    unix2ymds(secTAI93 + header.endTime, &year, &month, &day, &secs);

    fprintf(stdout, "rangeendingdate=%4d-%02d-%02d\n", year, month, day);

    hour = (int) (24. * (secs / 86400.));
    minute = (int) ((secs - 3600. * (double) hour) / 60.);
    second = (int) ((secs - 3600. * (double) hour - (double) minute * 60.));
    fsec = 1000000. * (secs - 3600. * (double) hour - (double) minute * 60. -
            (double) second);

    fprintf(stdout, "rangeendingtime=%02d:%02d:%02d.%d\n", hour, minute, second, (int) fsec);
    fprintf(stdout, "percentMissingData=%lf\n", header.qaStatistics[1]);

    return RTRN_SUCCEED;
}

int getAttHeaders(char *attfile) {
    FILE *attitFilePtr;
    PGSt_attitHeader header;
    char *ptr, *ptr1, *ptr2;
    int ll, ll1, n;
    static char buf[255];
    int16_t year, month, day, hour, minute, second;
    double secs;
    double secTAI93 = 725846400.;
    static char shortname[8];
    size_t len = 8;

    strcpy(buf, attfile);
    strncpy(shortname, buf, len);
    fprintf(stdout, "shortname=%s\n", shortname);
    fprintf(stdout, "localgranuleid=%s\n", buf);

    if ((attitFilePtr = fopen(buf, "r")) == NULL) {
        fprintf(stderr, "Error opening Attitude file\n");
        return RTRN_ERROR;
    } else {
        fread(&header, sizeof (PGSt_attitHeader), 1, attitFilePtr);
        fclose(attitFilePtr);
    }
    ptr = (char *) (&header.startTime);
    ll = sizeof (PGSt_double);
    byteswap(ptr, ll);

    ptr1 = (char *) (&header.endTime);
    ll1 = sizeof (PGSt_double);
    byteswap(ptr1, ll1);

    ptr2 = (char *) (&header.interval);
    for (n = 0; n < 26; n++, ptr2 += 4) byteswap(ptr2, 4);

    unix2ymds(secTAI93 + header.startTime, &year, &month, &day, &secs);

    fprintf(stdout, "rangebeginningdate=%4d-%02d-%02d\n", year, month, day);

    hour = (int) (24. * (secs / 86400.));
    minute = (int) ((secs - 3600. * (double) hour) / 60.);
    second = (int) ((secs - 3600. * (double) hour - (double) minute * 60.));
    double fsec = 1000000. * (secs - 3600. * (double) hour - (double) minute * 60. -
            (double) second);

    fprintf(stdout, "rangebeginningtime=%02d:%02d:%02d.%d\n", hour, minute, second, (int) fsec);

    unix2ymds(secTAI93 + header.endTime, &year, &month, &day, &secs);

    fprintf(stdout, "rangeendingdate=%4d-%02d-%02d\n", year, month, day);

    hour = (int) (24. * (secs / 86400.));
    minute = (int) ((secs - 3600. * (double) hour) / 60.);
    second = (int) ((secs - 3600. * (double) hour - (double) minute * 60.));
    fsec = 1000000. * (secs - 3600. * (double) hour - (double) minute * 60. -
            (double) second);

    fprintf(stdout, "rangeendingtime=%02d:%02d:%02d.%d\n", hour, minute, second, (int) fsec);
    fprintf(stdout, "percentMissingData=%lf\n", header.qaStatistics[1]);

    return RTRN_SUCCEED;
}

int main(int argc, char *argv[]) {
    char *attephfile;
    int rtrn;


    if (2 == argc) {
        attephfile = argv[1];

        if (strstr(attephfile, "ATT")) {
            rtrn = getAttHeaders(attephfile);
        } else {
            rtrn = getEphemHeaders(attephfile);
        }
    } else {
        printf("Usage: %s atteph-file\n", argv[0]);
        rtrn = RTRN_ERROR;
    }

    return rtrn;
}
