#include <timeutils.h>

/* -------------------------------------------------------------- */
/* yd2md() - converts year and day to month and day-of-month      */
/* -------------------------------------------------------------- */
void yd2md(int16_t year, int16_t doy, int16_t *month, int16_t *dom) {
    double sec = 0.0;
    int16_t yr;

    unix2ymds(yds2unix(year, doy, sec), &yr, month, dom, &sec);

    return;
}

