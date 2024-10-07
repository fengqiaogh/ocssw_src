#define LEAPSEC_DAT "leapsec.dat"
#define JD_1993 2448988.5
#define SECONDS_IN_DAY 86400
#define MAX_LEAP_SECONDS 32
#define UNIX_TIME_AT_1993 725846400.0

static double leap_seconds_tai93[MAX_LEAP_SECONDS] = {0};
static double leap_seconds_unix[MAX_LEAP_SECONDS] = {0};

#include <timeutils.h>
#include <stdlib.h>
#include <stdio.h>

void init_leapsecond_arrays() {
    char *varRoot;
    char leapsecdat[FILENAME_MAX];

    if (!leap_seconds_tai93[0]) {
        if ((varRoot = getenv("OCVARROOT")) == NULL) {
            printf("-E- OCVARROOT environment variable is not defined.\n");
            exit(EXIT_FAILURE);
        }
        snprintf(leapsecdat, FILENAME_MAX, "%s/common/tai-utc.dat", varRoot);
        FILE *file_h = fopen(leapsecdat, "rb");
        if (file_h == NULL) {
            snprintf(leapsecdat, FILENAME_MAX, "%s/viirsn/IETTime.dat", varRoot);
            file_h = fopen(leapsecdat, "rb");
            if (file_h == NULL) {
                snprintf(leapsecdat, FILENAME_MAX, "%s/modis/leapsec.dat", varRoot);
                file_h = fopen(leapsecdat, "rb");
                if (file_h == NULL) {
                    printf("-E- %s:%d - leap second file (%s/common/tai-utc.dat) not found.\n", __FILE__, __LINE__, varRoot);
                    exit(EXIT_FAILURE);
                }
            }
        }
        int i = 0, year;
        float jd;
        while (fgetc(file_h) != '\n')
            ;
        while (fscanf(file_h, " %d %*s %*d =JD %f%*[^\n]", &year, &jd) == 2) {
            if (year >= 1993) {
                leap_seconds_tai93[i] = ((jd - JD_1993) * SECONDS_IN_DAY) + i;
                leap_seconds_unix[i] =  ((jd - JD_1993) * SECONDS_IN_DAY) + UNIX_TIME_AT_1993;
                i++;
                if (i >= MAX_LEAP_SECONDS) {
                    printf("-E- %s:%d - max leap seconds reached in %s\n", __FILE__, __LINE__, leapsecdat);
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(file_h);
    }
    
}


int leapseconds_since_1993(double tai93time) {

    init_leapsecond_arrays();

    int i;
    for (i = 0; i < MAX_LEAP_SECONDS - 1 && leap_seconds_tai93[i]; i++) {
        if (tai93time < leap_seconds_tai93[i]) {
            return i;
        }
    }
    return i;
}

int leapseconds_since_1993_unix(double unixtime) {

    init_leapsecond_arrays();

    int i;
    for (i = 0; i < MAX_LEAP_SECONDS - 1 && leap_seconds_unix[i]; i++) {
        if (unixtime < leap_seconds_unix[i]) {
            return i;
        }
    }
    return i;
}
