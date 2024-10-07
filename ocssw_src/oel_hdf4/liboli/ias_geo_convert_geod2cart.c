/*******************************************************************************
NAME: ias_geo_convert_geod2cart

PURPOSE: Convert geodetic coordinate (lat, lon, height) into Cartesian
         coordinates (x, y, z).

RETURN VALUE:
none

NOTES: Input lat & lon are in radians, height, semi-major axis, and 
       output Cartesian position vector are in meters
       flattening is a unitless number.

ALGORITHM REFERENCES:
             Geodesy: The Concept, by Vanicek and Karkivsky, 1986;
             or any other geodesy textbook.

*******************************************************************************/
#include <math.h>
#include "ias_math.h"
#include "ias_structures.h"
#include "ias_geo.h"

void ias_geo_convert_geod2cart
(
    double latitude,    /* I: Lat of geodetic coordinates in radians */
    double longitude,   /* I: Long of geodetic coordinates in radians */
    double height,      /* I: Height (elevation) of geodetic coord in meters */
    double semimajor,   /* I: Reference ellipsoid semi-major axis in meters */
    double flattening,  /* I: Flattening of the ellipsoid 
                           (semimajor-semiminor)/semimajor */
    IAS_VECTOR *cart        /* O: Cartesian vector for the coord */
)
{
    double prime_vertical_radius;/* radius of prime vertical */
    double ecc2;               /* square of eccentricity */
    double coslat, sinlat;     /* cosine and sine of lat */

    /* Calculate the radius of prime vertical of the ellipsoid */ 
    ecc2 = flattening * (2.0 - flattening);
    coslat = cos(latitude);
    sinlat = sin(latitude);

    /* Since sinlat is less than or equal to 1 and ecc2 is less than 1 no
       test is performed for negative square root or division by zero */
    prime_vertical_radius = semimajor / sqrt(1.0 - ecc2 * sinlat * sinlat);

    /* Calculate the cartesian coordinates */
    cart->x = (prime_vertical_radius + height)*coslat*cos(longitude);
    cart->y = (prime_vertical_radius + height)*coslat*sin(longitude);
    cart->z = (prime_vertical_radius * (1.0 - flattening) * (1.0 - flattening) 
        + height) * sinlat;
}
