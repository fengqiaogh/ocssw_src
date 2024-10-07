/*****************************************************************
 *
 *  NSSDC/CDF							EPOCHbreakdown.
 *
 *  Version 2.0, 10-Jul-91, ST Systems (STX)
 *
 *  Modification history:
 *
 *   V1.0 ~11-Feb-91, J Love	Original version (for CDF V2.0).  Used
 *				algorithm from EPOCH_BREAKDOWN.F (author
 *				unknown).
 *   V1.1   9-Jul-91, J Love	Added include of 'math.h'.
 *   V2.0  10-Jul-91, J Love	Renamed (was EPOCH_BREAKDOWN.C in CDF V2.0
 *				distribution).
 *
 *******************************************************************/

#include <math.h>

void EPOCHbreakdown(epoch, year, month, day, hour, minute, second, msec)
double epoch;
long *year, *month, *day, *hour, *minute, *second, *msec;
{
    long jd, i, j, k, l, n;
    double msecAD, secondAD, minuteAD, hourAD, dayAD;

    msecAD = epoch;
    secondAD = msecAD / 1000.0;
    minuteAD = secondAD / 60.0;
    hourAD = minuteAD / 60.0;
    dayAD = hourAD / 24.0;

    jd = 1721060 + dayAD;
    l = jd + 68569;
    n = 4 * l / 146097;
    l = l - (146097 * n + 3) / 4;
    i = 4000 * (l + 1) / 1461001;
    l = l - 1461 * i / 4 + 31;
    j = 80 * l / 2447;
    k = l - 2447 * j / 80;
    l = j / 11;
    j = j + 2 - 12 * l;
    i = 100 * (n - 49) + i + l;

    *year = i;
    *month = j;
    *day = k;

    *hour = fmod(hourAD, 24.0);
    *minute = fmod(minuteAD, 60.0);
    *second = fmod(secondAD, 60.0);
    *msec = fmod(msecAD, 1000.0);

    return;
}
