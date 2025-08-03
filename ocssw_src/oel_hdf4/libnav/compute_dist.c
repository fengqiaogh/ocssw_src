/*
 * compute_dist.c
 *
 *  Created on: Jun 29, 2016
 *      Author: rhealy
 */
#include <math.h>
#include <genutils.h>

double deg2rad(double deg);

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Compute the angular distance                                  :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

double angular_distance(double lat1, double lon1, double lat2, double lon2) {
  double theta, dist;
  theta = lon1 - lon2;
  dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
  dist = dist> 1? 1:dist;
  dist = dist<-1?-1:dist;
  return (acos(dist));
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
  return (deg * OEL_DEGRAD);
}
/* and the inverse */
double rad2deg(double rad) {
  return (rad * OEL_RADEG);
}

