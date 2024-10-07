#ifndef  _TIME_UTILS_H
#define  _TIME_UTILS_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define SECONDS_IN_DAY 86400

#ifdef __cplusplus
extern "C" {
#endif

int isleap(int year);

time_t gmt_offset(void);
double yds2unix(int16_t year, int16_t day, double secs);
double yds2unix_(int16_t *year, int16_t *day, double *secs);
void unix2yds(double usec, int16_t *year, int16_t *day, double *secs);
char * unix2ydhmsf(double usec, char zone);
double ymds2unix(int16_t year, int16_t month, int16_t day, double secs);
void unix2ymds(double usec, int16_t *year, int16_t *mon, int16_t *day, double *secs);
void unix2ymdhms(double usec, int16_t *year, int16_t *mon, int16_t *day, int16_t *hour, int16_t *min, double *sec);
void yd2md(int16_t year, int16_t doy, int16_t *month, int16_t *dom);
char * ydhmsf(double dtime, char zone);
double now(void);
void ymdhms2ydmsec(int yy, int mm, int dd, int hh, int mn, int sc,
        int32_t *year, int32_t *day, int32_t *msec);
void date2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec);
void isodate2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec);
char * unix2isodate(double dtime, char zone);
double zulu2unix(char *zulu);
double isodate2unix(const char *isodate);
void addmsec(int32_t *year, int32_t *day, int32_t *msec, int32_t delta);

void get_time(char *pr_time);

int leapseconds_since_1993(double tai93time);
int leapseconds_since_1993_unix(double unixtime);
double unix_to_tai93(double unixtime);
double tai93_to_unix(double tai93);
double unix_to_tai58(double unixtime);
double tai58_to_unix(double tai58);

int32_t jday(int16_t year, int16_t month, int16_t day);
int jdate(int32_t julianDay, int32_t *year, int32_t *dayOfYear);
double yds2tai93(int16_t iyr, int16_t idy, double sec);
int ccsds_to_yds(uint8_t *cctime, int32_t *iyear, int32_t *iday, double *sec);

double udunits2unix(const char* str);
const char* unix2udunits_c(double unixtime, int hasZ);

#ifdef __cplusplus
}  // end extern "C"

#include <string>

double udunits2unix(const std::string& str);
std::string unix2udunits(double unixtime, bool hasZ = true);

#endif

#endif
