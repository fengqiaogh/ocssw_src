#include "libnav.h"
#include "sun2000.h"
#include "nav.h"
/**
 * @brief c  Computes unit Sun vector in geocentric rotating coodinates, using
          c  subprograms to compute inertial Sun vector and Greenwich hour angle
c       Subprograms referenced:
c
c       SUN2000         Computes inertial Sun vector
c       GHA2000         Computes Greenwich sidereal angle
c
c       Coded by:  Frederick S. Patt, GSC, September 29, 1992
c
c       Modification History:
c
c       Modifified to use new Sun and hour angle routines
c       Frederick S. Patt, November 3, 1992
c
c       Removed internal jd() function, since it is available as an
c       independent module.  B. A. Franz, GSC, November 14, 1997.
 * @param iyr input      Year, four digits (i.e, 1993)
 * @param idoy input  Day of year (1-366)
 * @param sec input Seconds of day
 * @param sunr output Unit Sun vector in geocentric rotating coordinates
 * @param rs Earth-to-Sun distance (AU)
 */
extern "C" void l_sun_(int *iyr, int *iday, double *sec, float sunr[3], float *rs) {
    getUnitRotatingSunVector(*iyr, *iday, *sec, sunr, *rs);
}

/**
 * @brief Computes unit Sun vector in geocentric rotating coodinates, using subprograms to compute inertial
 * Sun vector and Greenwich hour angle 
 * 
 * Subprograms referenced:
 * SUN2000; Computes inertial Sun vector
 * GHA2000; Computes Greenwich sidereal angle
 *
 * Coded by:  Frederick S. Patt, GSC, September 29, 1992
 *
 * Modification history:
 * Modified to use new Sun and hour angle routines
 * Frederick S. Patt, November 3, 1992
 *
 * Removed internal jd() function, since it is available as an
 * independent module.  B. A. Franz, GSC, November 14, 1997.
 * 
 * Conveted to C++.  Jakob C Lindo, SSAI, September 30, 2025
 * @param[in] year Year, four digits (i.e, 1993)
 * @param[in] dayOfYear Day of year (1-366)
 * @param[in] secondsOfDay Seconds of day
 * @param[out] sunRotatingVector Unit Sun vector in geocentric rotating coordinates
 * @param[out] earthToSunDistance Earth-to-Sun distance (in AU)
 */
template <typename T>
void getUnitRotatingSunVector(int &year, int &dayOfYear, double &secondsOfDay, T rotatingSunVector[3],
                              T &earthToSunDistance) {
    double inertialUnitSunVector[3];  // Unit vector
    double esDistance;
    getUnitInertialSunVector(year, dayOfYear, secondsOfDay, inertialUnitSunVector, esDistance);
    earthToSunDistance = esDistance;

    double dayWithFraction = dayOfYear + secondsOfDay / 864.e2;
    double greenwichHourAngle;
    gha2000(year, dayWithFraction, &greenwichHourAngle);
    double greenwichHourAngleRotating = greenwichHourAngle / radeg;

    double sinGhar = sin(greenwichHourAngleRotating);
    double cosGhar = cos(greenwichHourAngleRotating);

    rotatingSunVector[0] = inertialUnitSunVector[0] * cosGhar + inertialUnitSunVector[1] * sinGhar;
    rotatingSunVector[1] = inertialUnitSunVector[1] * cosGhar - inertialUnitSunVector[0] * sinGhar;
    rotatingSunVector[2] = inertialUnitSunVector[2];
}

// force the instanciation of the float and double version
template void getUnitRotatingSunVector<float>(int &year, int &dayOfYear, double &secondsOfDay, float rotatingSunVector[3],
                              float &earthToSunDistance);
template void getUnitRotatingSunVector<double>(int &year, int &dayOfYear, double &secondsOfDay, double rotatingSunVector[3],
                              double &earthToSunDistance);
