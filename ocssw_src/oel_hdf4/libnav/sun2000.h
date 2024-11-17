#ifndef OCSSW_SUN2000_H
#define OCSSW_SUN2000_H


/**
 * @brief
This subroutine computes the Sun vector in geocentric inertial
(equatorial) coodinates.  It uses the model referenced in The
Astronomical Almanac for 1984, Section S (Supplement) and documented
in "Exact closed-form geolocation algorithm for Earth survey
sensors", by F.S. Patt and W.W. Gregg, Int. Journal of Remote
Sensing, 1993.  The accuracy of the Sun vector is approximately 0.1
       Subprograms referenced:
       JD              Computes Julian day from calendar date
       EPHPARMS        Computes mean solar longitude and anomaly and
                        mean lunar lontitude and ascending node
       NUTATE          Compute nutation corrections to lontitude and
                        obliquityc
       Coded by:  Frederick S. Patt, GSC, November 2, 1992
       Modified to include Earth constants subroutine by W. Gregg,
               May 11, 1993.
arcminute.
 * @param iyr Year, four digits (i.e, 1993)
 * @param iday Day of year (1-366)
 * @param sec Seconds of day
 * @param sun  Unit Sun vector in geocentric inertial  coordinates of date
 * @param rs  Magnitude of the Sun vector (AU)
 */
void sun2000(int iyr, int iday, double sec, float sun[3], float* rs);

/**
 @brief
This subroutine computes the nutation in longitude and the obliquity
of the ecliptic corrected for nutation.  It uses the model referenced
in The Astronomical Almanac for 1984, Section S (Supplement) and
documented in "Exact closed-form geolocation algorithm for Earth
survey sensors", by F.S. Patt and W.W. Gregg, Int. Journal of
Remote Sensing, 1993.  These parameters are used to compute the
apparent time correction to the Greenwich Hour Angle and for the
calculation of the geocentric Sun vector.  The input ephemeris
parameters are computed using subroutine ephparms.  Terms are
included to 0.1 arcsecond.
       Program written by:     Frederick S. Patt
                               General Sciences Corporation
                               October 21, 1992

 * @param t Time in days since January 1, 2000 at  12 hours UT
 * @param xls Mean solar longitude (degrees)
 * @param gs Mean solar anomaly   (degrees)
 * @param xlm Mean lunar longitude (degrees)
 * @param omega Ascending node of mean lunar orbit  (degrees)
 * @param dpsi Nutation in longitude (degrees)
 * @param eps Obliquity of the Ecliptic (degrees) (includes nutation in obliquity)
 *
 */
void nutate(double t, double xls, double gs, double xlm, double omega, double* dpsi, double* eps);

/**
 * @brief
This subroutine computes the Greenwich hour angle in degrees for the
input time.  It uses the model referenced in The Astronomical Almanac
for 1984, Section S (Supplement) and documented in "Exact
closed-form geolocation algorithm for Earth survey sensors", by
F.S. Patt and W.W. Gregg, Int. Journal of Remote Sensing, 1993.
It includes the correction to mean sideral time for nutation
as well as precession.
       Subprograms referenced:

       JD              Computes Julian day from calendar date
       EPHPARMS        Computes mean solar longitude and anomaly and
                        mean lunar lontitude and ascending node
       NUTATE          Compute nutation corrections to lontitude and
                        obliquity


       Program written by:     Frederick S. Patt
                               General Sciences Corporation
                               November 2, 1992
 * @param iyr Year (four digits)
 * @param day Day (time of day as fraction)
 * @param gha reenwich hour angle (degrees)
 */
void gha2000(int iyr, double day, double* gha);

/**
 * @brief
This subroutine computes ephemeris parameters used by other Mission
Operations routines:  the solar mean longitude and mean anomaly, and
the lunar mean longitude and mean ascending node.  It uses the model
referenced in The Astronomical Almanac for 1984, Section S
(Supplement) and documented and documented in "Exact closed-form
geolocation algorithm for Earth survey sensors", by F.S. Patt and
W.W. Gregg, Int. Journal of Remote Sensing, 1993.  These parameters
are used to compute the solar longitude and the nutation in
longitude and obliquity.
 * @param t Time in days since January 1, 2000 at 12 hours UT
 * @param xls Mean solar longitude (degrees)
 * @param gs Mean solar anomaly (degrees)
 * @param xlm Mean lunar longitude (degrees)
 * @param omega Ascending node of mean lunar orbit (degrees)
 */
void ephparms(double t, double* xls, double* gs, double* xlm, double* omega);

#endif  // OCSSW_SUN2000_H
