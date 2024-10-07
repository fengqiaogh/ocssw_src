#include <stdio.h>
#include <stdint.h>

/**
 * @brief Converts a calendar date to the corresponding Julian day starting at noon on the calendar date.
 * Originally implemented by Frederick S. Patt, GSC, November 4 1992.
 * @cite  Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, November 1979, p. 400.
 *
 * @param year 4-digit year
 * @param month 1-12
 * @param day 1-31
 * @return The astronomical Julian day corresponding to the one starting at noon on the calendar day
 */
int32_t jday(int16_t year, int16_t month, int16_t day) {
    return 367 * year - 7 * (year + (month + 9) / 12) / 4 + 275 * month / 9 + day + 1721014;
}
