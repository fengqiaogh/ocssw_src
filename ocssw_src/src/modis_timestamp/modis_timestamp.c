#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>

int main(int argc, char *argv[]) {
    int len;
    char *ifile = NULL;
    FILE *fp = NULL;
    FILE *fptr = NULL;
    FILE *testptr = NULL;
    char ncdump[500000] = "";
    char tmpstr[1000] = "";
    char cmd_string[250] = "";
    char *ptr = NULL;
    char range_date[11] = "";
    char range_time[9] = "";
    char platform[4] = "";
    char first_letter;
    char year[5] = "";
    char month[3] = "";
    char day[3] = "";
    int year_int;
    int day_int;
    int mod;
    int add_leap_day;
    int basedays = 0;
    int julian_days;
    char julian_str[4] = "";
    char jul_str[4] = "";
    char hour[3] = "";
    char min[3] = "";
    char sec[3] = "";
    char date_stamp[15] = "";
    char *seadas = NULL;

    if (argc == 1) {
        printf("\nUSAGE: %s MODIS_L1A_L1B_or_GEO_file [start|stop]\n", argv[0]);
        exit(1);
    }
    if (!(fptr = fopen(argv[1], "r"))) {
        perror(argv[1]);
        exit(1);
    }
    fclose(fptr);

    len = strlen(argv[1]);
    ifile = malloc(1 + (len * sizeof (char)));
    strncpy(ifile, argv[1], len);

    seadas = getenv("OCSSWROOT");
    if (seadas == NULL) {
        strcat(cmd_string, "ncdump");
    } else {
        strcpy(cmd_string, seadas);
        strcat(cmd_string, "/run/bin3/ncdump");
        if (!(testptr = fopen(cmd_string, "r"))) {
            strcpy(cmd_string, seadas);
            strcat(cmd_string, "/bin/ncdump");
            if (!(testptr = fopen(cmd_string, "r"))) {
                printf("modis_timestamp: Error - ncdump must exist in $SEADAS/run/bin3/ or $SEADAS/bin/ to run this program.\n");
                exit(1);
            }
        }
    }

    strcat(cmd_string, " -h ");
    strcat(cmd_string, ifile);

    fp = popen(cmd_string, "r");
    while (fgets(tmpstr, sizeof tmpstr, fp)) {
        strcat(ncdump, tmpstr);
    }
    pclose(fp);

    /* find the date */
    if (argc == 2) argv[2] = argv[1]; /* HACK - do this to avoid bus error */
    if (strcmp(argv[2], "stop") == 0) ptr = strstr(ncdump, "RANGEENDINGDATE");
    else ptr = strstr(ncdump, "RANGEBEGINNINGDATE");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    ptr = strstr(ptr, "VALUE");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    strncpy(range_date, ptr + 25, 10);
    range_date[10] = '\0';

    /* find the time */
    if (strcmp(argv[2], "stop") == 0) ptr = strstr(ncdump, "RANGEENDINGTIME");
    else ptr = strstr(ncdump, "RANGEBEGINNINGTIME");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    ptr = strstr(ptr, "VALUE");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    strncpy(range_time, ptr + 25, 8);
    range_time[8] = '\0';

    /* find out if it's aqua or terra */
    ptr = strstr(ncdump, " SHORTNAME");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    ptr = strstr(ptr, "VALUE");
    if (ptr == NULL) {
        printf("modis_timestamp: Error - File doesn't appear to have a valid MODIS HDF header.\n");
        exit(1);
    }
    strncpy(platform, ptr + 25, 3);
    platform[3] = '\0';

    /* build the return string, like AYYYYDDDHHMMSS */
    if (strcmp(platform, "MYD") == 0) {
        first_letter = 'A';
    } else if (strcmp(platform, "MOD") == 0) {
        first_letter = 'T';
    } else {
        printf("modis_timestamp: Error - Could not determine platform.\n");
        exit(1);
    }

    strncpy(year, &range_date[0], 4);
    year[4] = '\0';
    strncpy(month, &range_date[5], 2);
    month[2] = '\0';
    strncpy(day, &range_date[8], 2);
    day[2] = '\0';

    /* convert year string into integer for leap calc */
    sscanf(year, "%d", &year_int);
    mod = year_int % 4;
    if (mod == 0) {
        add_leap_day = 1;
    } else {
        add_leap_day = 0;
    }

    sscanf(day, "%d", &day_int);

    if (strcmp(month, "01") == 0)
        basedays = 0;
    else if (strcmp(month, "02") == 0)
        basedays = 31;
    else if (strcmp(month, "03") == 0)
        basedays = 59 + add_leap_day;
    else if (strcmp(month, "04") == 0)
        basedays = 90 + add_leap_day;
    else if (strcmp(month, "05") == 0)
        basedays = 120 + add_leap_day;
    else if (strcmp(month, "06") == 0)
        basedays = 151 + add_leap_day;
    else if (strcmp(month, "07") == 0)
        basedays = 181 + add_leap_day;
    else if (strcmp(month, "08") == 0)
        basedays = 212 + add_leap_day;
    else if (strcmp(month, "09") == 0)
        basedays = 243 + add_leap_day;
    else if (strcmp(month, "10") == 0)
        basedays = 273 + add_leap_day;
    else if (strcmp(month, "11") == 0)
        basedays = 304 + add_leap_day;
    else if (strcmp(month, "12") == 0)
        basedays = 334 + add_leap_day;

    julian_days = basedays + day_int;

    if (julian_days < 10) {
        strcpy(julian_str, "00\0");
    } else if (julian_days < 100) {
        strcpy(julian_str, "0\0");
    } else {
        strcpy(julian_str, "\0");
    }

    sprintf(jul_str, "%d", julian_days);
    strcat(julian_str, jul_str);

    strncpy(hour, range_time, 2);
    hour[2] = '\0';
    ptr = &range_time[3];
    strncpy(min, ptr, 2);
    min[2] = '\0';
    ptr = &range_time[6];
    strncpy(sec, ptr, 2);
    sec[2] = '\0';

    /* put it all together */
    date_stamp[0] = first_letter;
    date_stamp[1] = '\0';
    strcat(date_stamp, year);
    strcat(date_stamp, julian_str);
    strcat(date_stamp, hour);
    strcat(date_stamp, min);
    strcat(date_stamp, sec);
    date_stamp[14] = '\0';

    printf("%s", date_stamp);

    exit(0);
}
