#include <stdio.h>
#include <stdint.h>
#include <timeutils.h>

/**
 * @brief Takes in a julian day and mutates `year` and `doy` to contain the gregorian year and day of that 
 * gregorian year of that julian day
 * 
 * @param julianDay The Julian day to be converted
 * @param year Out-parameter indicating the gregorian year that `julianDay` is in
 * @param dayOfYear Out-parameter indicating the number of days into the gregorian year that `julianDay` is. 
*/
int jdate(int32_t julianDay, int32_t *year, int32_t *dayOfYear) {

    int32_t month, day;
    int32_t ja, jb, jc, jd, je, jalpha;

    jalpha = (int32_t) (((julianDay - 1867216L) - 0.25) / 36524.25);
    ja = julianDay + 1 + jalpha - (int32_t) (0.25 * jalpha);

    jb = ja + 1524;
    jc = (int32_t) (6680.0 + ((jb - 2439870L) - 122.1) / 365.25);
    jd = 365 * jc + (int32_t) (0.25 * jc);
    je = (int32_t) ((jb - jd) / 30.6001);

    day = jb - jd - (int32_t) (30.6001 * je);
    month = je - 1;
    if (month > 12)
        month = month - 12;
    *year = jc - 4715;

    if (month > 2) 
        *year = *year - 1;
    if (*year <= 0)
        *year = *year - 1;

    int numDaysEachMonth[24] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    *dayOfYear = day;
    for (int idxMonth = 1; idxMonth < month; idxMonth++) 
        *dayOfYear = *dayOfYear + numDaysEachMonth[isleap(*year) * 12 + (idxMonth - 1)];

    return 0;
}
