#include <timeutils.h>

/****************************************************************************
Return the current unix time as a double.
 *****************************************************************************/
double now(void) {
    /* I could call gettimeofday() here to resolve the time down to the
     *  millisecond level, but I can't imagine that anyone really needs
     *  that kind of precision. */
    return ((double) time(NULL));
}


