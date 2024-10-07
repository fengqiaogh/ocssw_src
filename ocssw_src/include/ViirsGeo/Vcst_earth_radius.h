/*******************************************************************************
 *
 * NAME: Vcst_earth_radius
 *
 * DESCRIPTION: Computes the radius of the earth, in kilometers, from the
 * geodetic latitude.
 *
 * USAGE:  radius = earth_radius_D(rlat);
 *
 * PARAMETER          TYPE  I/O  DESCRIPTION
 * ----------- -----------  ---  --------------------------------
 * radius           double   O   radius of earth in kilometers
 * rlat             double   I   geodetic latitude, radians, positive north
 * ----------- -----------  ---  --------------------------------
 *
 *
 * DESCRIPTION: Computes the radius of the earth, in kilometers, from the
 * geocentric latitude.
 *
 * USAGE:  radius = Vcst_earth_radius_C(clat);
 *
 * PARAMETER          TYPE  I/O  DESCRIPTION
 * ----------- -----------  ---  --------------------------------
 * radius           double   O   radius of earth in kilometers
 * clat             double   I   geocentric latitude, radians, positive north
 * ----------- -----------  ---  --------------------------------
 *
 * REFERENCES: USAF, Orbit Analyst Manuals, 1978
 * The main reference for the WGS84 geoid is NIMA Physical Geodesy
 * web page: 164.214.2.59/GandG/wgs-84/egm96.htm
 *
 * LIMITATIONS: Because this function is called millions of times, there is no
 * error checking. Therefore, it is the responsibility of the calling program
 * to ensure legal input values:  -pio2 < rlat < pio2
 *
 * NOTES (MISCELLANEOUS) SECTION: 
 *
 *******************************************************************************/

#ifndef _EARTH_RADIUS_H
#define _EARTH_RADIUS_H


#ifdef __cplusplus
extern "C" {
#endif

double earth_radius_D(double rlat);

double earth_radius_C(double clat);

#ifdef __cplusplus
}
#endif

#endif

