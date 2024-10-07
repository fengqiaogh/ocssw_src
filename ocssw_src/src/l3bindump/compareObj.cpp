#include "compareObj.h"
#include <cstdio>
#include <L3Shape.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>

using namespace l3;

CompareObjLonOnly::CompareObjLonOnly(double east, double west) {
    this->east = constrainLon(east);
    this->west = constrainLon(west);
    if (this->east < this->west)
        this->east += 360.0;
}

bool CompareObjLonOnly::isInside(double lat, double lon) {
    lon = constrainLon(lon);
    if (lon < west)
        lon += 360.0;
    if (lon <= east)
        return true;
    return false;
}

CompareObjCircle::CompareObjCircle(double lat, double lon, double radius) {
    lat0 = lat;
    lon0 = lon;
    this->radius = radius;
}

bool CompareObjCircle::isInside(double lat, double lon) {
    double const earth_radius = 6371.229;

    using namespace boost::geometry;
    typedef model::point<double, 2, cs::spherical_equatorial<degree> > P;
    model::linestring<P> line;
    line.push_back(P(lon0, lat0));
    line.push_back(P(lon, lat));
    double dist = length(line, strategy::distance::haversine<double>(earth_radius));
    if (dist <= radius)
        return true;
    return false;
}
