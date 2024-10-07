/*
 *----------------------------------------------------------------------
 *  @(#) time-utils.c	1.0	05 Jun 98	<shc>
 *  Copyright (c) 1993, CSIRO Division of Oceanography.
 *  Copyright (c) 1998, Datron Transco Inc.
 *----------------------------------------------------------------------
 *
 * time-utils --
 *
 *	A collection of routines for converting between normal time/date
 *	formats and Julian dates.
 *
 * History:
 *	05 Jun 98  <shc>
 *	   Split out from old julian.c
 *
 *----------------------------------------------------------------------
 */

#include <math.h>
#include <time.h>

#include "orbit.h" /* Prototype checking */

/*
 *  This function converts the date and time in a tm struct to the
 *  modified Julian day format used internally in all the astronomical
 *  and orbit prediction code.
 */

double
tmtojul(time)
struct tm *time;
{
    int year, month, day;
    int32_t hour, minute, second;
    int32_t mjd;
    double fday;
    double result;

    day = time->tm_mday;
    month = time->tm_mon + 1; /* tm_mon is in range 0-11 */
    year = 1900 + time->tm_year; /* tm_year is years since 1900 */

    second = time->tm_sec;
    minute = time->tm_min;
    hour = time->tm_hour;

    mjd = julday(year, month, day) - 2400000;

    fday = ((hour * 60 + minute)*60 + second) / 86400.0;

    result = mjd + fday - 0.5;
    return (result);
}

/*
 *  This function does the inverse of TOJUL.
 */

void
jultotm(tjd, time)
double tjd;
struct tm *time;
{
    int32_t mjd;
    double secs;
    int year, month, day;
    int hour, minute, second;

    mjd = (int32_t) (tjd + 0.5);
    secs = (tjd + 0.5 - mjd) * 86400.0 + 0.0005;
    mjd += 2400000;

    caldat(mjd, &year, &month, &day);

    hour = (int) (secs / 3600.0);
    minute = (int) (secs / 60.0) % 60;
    second = (int) (secs) % 60;

    time->tm_sec = second;
    time->tm_min = minute;
    time->tm_hour = hour;
    time->tm_mday = day;
    time->tm_mon = month;
    time->tm_year = year;
    time->tm_isdst = 0;

    return;
}

/*
 *  This function converts a date and time in compact integer form to the
 *  modified Julian day format used internally in all the orbit prediction
 *  code.  Note that a year number of less than 50 is interpreted as being
 *  in the range 2000 - 2049.
 */

double
itojul(date, time)
int32_t date; /* yymmdd */
int32_t time; /* hhmmssfff */
{
    int year, month, day, hour, minute, second, millisec;
    int32_t mjd;
    double fday, result;

    day = date % 100L;
    month = date / 100L % 100L;
    year = date / 10000L % 100L;
    millisec = time % 1000L;
    second = time / 1000L % 100L;
    minute = time / 100000L % 100L;
    hour = time / 10000000L % 100L;

    if (year >= 50)
        year = year + 1900;
    else
        year = year + 2000;

    mjd = julday(year, month, day) - 2400000L;

    fday = (((hour * 60.0 + minute)*60.0 + second) + millisec / 1000.0) / 86400.0;

    result = mjd + fday - 0.5;
    return (result);
}

/*
 *  This function does the inverse of ITOJUL.
 */

void
jultoi(tjd, date, time)
double tjd;
int32_t *date, *time;
{
    int year, month, day;
    int32_t mjd, hour, minute, second, millisec;
    double secs;

    mjd = (int32_t) (tjd + 0.5);
    secs = (tjd + 0.5 - mjd) * 86400.0 + 0.0005;
    caldat(mjd + 2400000L, &year, &month, &day);

    hour = (int32_t) (secs / 3600.0);
    minute = (int32_t) (secs / 60.0) % 60L;
    second = (int32_t) (secs) % 60L;
    millisec = (int32_t) (secs * 1000.0) % 1000L;

    *date = ((int32_t) year % 100L)*10000L + (int32_t) month * 100L + (int32_t) day;
    *time = hour * 10000000L + minute * 100000L + second * 1000L + millisec;

    return;
}
