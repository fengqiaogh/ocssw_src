/**

An implementation of a Lon/Lat to UTM zone code written by Chuck Gantz <chuck.gantz@globalstar.com>
...in 1998...found on http://www.gpsy.com/gpsinfo/geotoutm/, a site published by Karen Nakamura, 
last updated 22 June 2000.

Minor modifications to the code were made to eliminate deprecation warnings 

 */

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "LatLong-UTMconversion.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " Lon Lat" << std::endl;
        return 1;
    }

    int ReferenceEllipsoid = 23; //WGS84

    double Lon, Lat;
    sscanf(argv[1], "%lf", &Lon);
    sscanf(argv[2], "%lf", &Lat);
    double UTMNorthing;
    double UTMEasting;
    char UTMZone[10];

    LLtoUTM(ReferenceEllipsoid, Lat, Lon, UTMNorthing, UTMEasting, UTMZone);

    std::cout << "+zone=" << UTMZone << std::endl;
    std::cout << "easting=" << std::fixed << std::setprecision(0)
            << UTMEasting << std::endl;
    std::cout << "northing=" << std::fixed << std::setprecision(0)
            << UTMNorthing << std::endl;
}
