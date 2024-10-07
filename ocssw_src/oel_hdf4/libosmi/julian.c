/*
 *----------------------------------------------------------------------
 *  @(#) julian.c	1.0	20 Mar 93	<shc>
 *  Copyright (c) 1993, CSIRO Division of Oceanography.
 *  Copyright (c) 1998, Datron Transco Inc..
 *----------------------------------------------------------------------
 *
 * julian --
 *
 *	A collection of routines for converting between normal time/date
 *	formats and Julian dates.
 *
 * History:
 *	25 Oct 93  <shc>
 *	   Initial version
 *
 *	05 Jun 98  <shc>
 *	   Divided into julday/caldat and other utility routines.
 *
 *----------------------------------------------------------------------
 */

#include <math.h>

#include "orbit.h" /* Prototype checking */


/*
 *  This routine returns the Julian Day Number which begins at noon of
 *  the calendar date specified by the arguments.  Positive year signifies
 *  A.D., negative B.C.  Remember that the year after 1 B.C. was 1 A.D.
 */


#define IGREG1 (15L + 31L * (10L + 12L * 1582L))

int32_t julday(year, month, day)
int year, month, day;
{
    int32_t y, m, adj, date;

    if (year == 0)
        return (-1);

    if (year < 0)
        year += 1;

    if (month > 2) {
        y = year;
        m = month + 1;
    } else {
        y = year - 1;
        m = month + 13;
    }

    date = (int32_t) floor(365.25 * y) + (int32_t) (30.6001 * m) +
            (int32_t) day + 1720995L;

    if ((int32_t) day + 31 * ((int32_t) month + 12 * (int32_t) year) > IGREG1) {
        adj = (int32_t) (0.01 * y);
        date = date + 2 - adj + (int32_t) (0.25 * adj);
    }

    return (date);
}


/*
 *  Inverse of the function "julday" above.
 */

#define IGREG2 2299161L

void
caldat(julian, year, month, day)
int32_t julian;
int *year, *month, *day;
{
    int32_t ja, jb, jc, jd, je, jalpha;

    if (julian < IGREG2) {
        ja = julian;
    } else {
        jalpha = (int32_t) (((julian - 1867216L) - 0.25) / 36524.25);
        ja = julian + 1 + jalpha - (int32_t) (0.25 * jalpha);
    }

    jb = ja + 1524;
    jc = (int32_t) (6680.0 + ((jb - 2439870L) - 122.1) / 365.25);
    jd = 365 * jc + (int32_t) (0.25 * jc);
    je = (int32_t) ((jb - jd) / 30.6001);

    *day = jb - jd - (int32_t) (30.6001 * je);
    *month = je - 1;
    if (*month > 12)
        *month = *month - 12;
    *year = jc - 4715;

    if (*month > 2)
        *year = *year - 1;
    if (*year <= 0)
        *year = *year - 1;

    return;
}
