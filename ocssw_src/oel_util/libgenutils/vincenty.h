/** @file vincenty.h
 @brief Provides a single function to calculate geographical distances.

 The Vincenty algorithm is described at:
 https://en.wikipedia.org/wiki/Vincenty's_formulae#Inverse_problem
 */

#ifndef INCLUDE_VINCENTY_H_
#define INCLUDE_VINCENTY_H_

/** @brief Calculate geographical distances using Vincenty's algorithm.

 @param[in] lat1 Latitude, in degrees, of first point.
 @param[in] lon1 Longitude, in degrees, of first point.
 @param[in] lat2 Latitude, in degrees, of second point.
 @param[in] lon2 Longitude, in degrees, of second point.
 @param[out] meters Distance between the two points, in meters.

 @return 0 on success.
 */
double vincenty_distance(double lat1, double lon1, double lat2, double lon2);

#endif /* INCLUDE_VINCENTY_H_ */
