#include <stdint.h>
#include <timeutils.h>

double yds2tai93(int16_t iyr, int16_t idy, double sec) {
    // add the seconds after the conversion to TAI93 so the leap second 
    // is not lost
    double unixTime = yds2unix(iyr, idy, 0.0);
    return unix_to_tai93(unixTime) + sec;
}

double unix_to_tai93(double unixtime) {
    double tai93 = unixtime - 725846400; // seconds between 1970 and 1993
    return tai93 + leapseconds_since_1993_unix(unixtime); // TAI includes leap seconds
}

double tai93_to_unix(double tai93) {
    double unixtime = tai93 - leapseconds_since_1993(tai93);
    return unixtime + 725846400;
}

// The TAI58 routines will have the wrong number of leap seconds
// for times before 1993.
double unix_to_tai58(double unixtime) {
    // seconds between UTC 1993 and UTC 1958 = 1104537600
    // plus leap seconds between 1993 and 1958 = 27
    return unix_to_tai93(unixtime) + 1104537600.0 + 27.0;
}

double tai58_to_unix(double tai58) {
    return tai93_to_unix(tai58 - 1104537600.0 - 27.0);
}
