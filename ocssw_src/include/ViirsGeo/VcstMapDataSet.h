/**************************************************************************
 * 
 * NAME: VcstMapDataSet
 *
 **************************************************************************/

#ifndef _VcstMapDataSet_H_
#define _VcstMapDataSet_H_

#include <algorithm>
#include <VcstCmnConsts.h>
#include <mds.h>
#include <novas.h>

//------------------------------------------------------------------------------
// Constants for instantiating N. & S. Mds.  The values for the lat and lon are
// just in the opposite hemisphere that is being used, so that the MDS
// encapsulates the entire hemisphere.
//------------------------------------------------------------------------------
const double NORTH_UPPER_LEFT_LAT = (DEG2RAD * -20.747206768662e0);
const double NORTH_UPPER_LEFT_LON = (DEG2RAD * 145.0e0);
const double NORTH_LOWER_RIGHT_LAT = NORTH_UPPER_LEFT_LAT;
const double NORTH_LOWER_RIGHT_LON = (NORTH_UPPER_LEFT_LON - PI);
const double SOUTH_UPPER_LEFT_LAT = -NORTH_UPPER_LEFT_LAT;
const double SOUTH_UPPER_LEFT_LON = (DEG2RAD * (-125.0e0));
//const double SOUTH_LOWER_RIGHT_LAT  =  SOUTH_UPPER_LEFT_LAT;
const double SOUTH_LOWER_RIGHT_LON = (SOUTH_UPPER_LEFT_LON + PI);

//------------------------------------------------------------------------------
// Constant used to allow for small float-to-double conversion errors
// when checking values against PI.
//------------------------------------------------------------------------------
const double PI_EPSILON = 1.0e-7;

/**
 * Enum datatype defining Hemisphere
 */
enum Hemisphere {
    NORTHERN, SOUTHERN, EQUATOR
};

/**
 * Enum datatype defining GridType
 */
enum GridType {
    GRIDTYPE_NONE = 0,
    MERCATOR_COAX = 11,
    POLAR_STEREO_NORTH = 21,
    POLAR_STEREO_SOUTH = 25,
    LAMBERT_TANGENT_NORTH = 31,
    LAMBERT_TANGENT_SOUTH = 35,
    LAMBERT_SECANT_NORTH = 41,
    LAMBERT_SECANT_SOUTH = 45,
    CYLINDRICAL_EQUIDISTANT = 51,
    POLAR_AZMUTHAL_NORTH = 61,
    POLAR_AZMUTHAL_SOUTH = 65
};

/**
 * This class encapsulates a Map Data Set (MDS), which consists of an area of
 * interest defined on a map of the earth.  An MDS does not contain any
 * particular data plotted on a map-it is simply a defined region. The map on
 * which the MDS is defined may have any projection, although the two most
 * commonly used are a mercator projection (flat-earth) and a
 * polar-stereographic projection (either northern or southern hemisphere). 
 */
class VcstMapDataSet {
public:

    /**
     * Constructor
     *
     * @param theHemisphere the hemisphere for the map dataset
     */
    VcstMapDataSet(Hemisphere theHemisphere);

    /**
     * Virtual Destructor
     */
    virtual ~VcstMapDataSet();

    /**
     * fills anMds with the values from this instance
     *
     * This method will not change where anMds points,
     * but it will change the values inside the mds_type struct.
     * If anMds == 0, this method is a no-op.
     *
     * @param anMds Ptr to mds_type (map dataset struct)
     */
    void fillMds(mds_type * const anMds) const;

    /**
     * Method returns upper left latitude value
     */
    /* inline */
    double getUpperLeftLatitude() const;

    /**
     * Method returns upper left longitude value
     */
    /* inline */
    double getUpperLeftLongitude() const;

    /**
     * Method returns lower right latitude value
     */
    /* inline */
    double getLowerRightLatitude() const;

    /**
     * Method returns lower right longitude value
     */
    /* inline */
    double getLowerRightLongitude() const;

    /**
     * Method returns largest latitude value
     */
    /* inline */
    double getBiggestLatitude() const;

    /**
     * Method returns largest longitude value
     */
    /* inline */
    double getBiggestLongitude() const;

    /**
     * Method returns smallest latitude value
     */
    /* inline */
    double getSmallestLatitude() const;

    /**
     * Method returns smallest longitude value
     */
    /* inline */
    double getSmallestLongitude() const;

    /**
     * Method returns Hemisphere enum value
     */
    /* inline */
    Hemisphere getHemisphere() const;

    /**
     * Method returns GridType enum value
     */
    /* inline */
    int getGridType() const;

    /**
     * Converts the input latitude and longitude to a grid position
     *
     * @param rlat Latitude 
     * @param rlon Longitude
     * @param row Grid row of latitude/longitude (rlat/rlon)
     * @param column Grid column of latitude/longitude (rlat/rlon)
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int latLonToGrid(
            double rlat,
            double rlon,
            double* row,
            double* column) = 0;

    /**
     * Converts the input grid position to a latitude and longitude.
     *
     * @param row Grid row
     * @param column Grid column
     * @param rlat Latitude of row/col
     * @param rlon Longitude of row/col
     * @return PRO_SUCCESS or PRO_FAIL
     */
    virtual int gridToLatLon(
            double row,
            double column,
            double* rlat,
            double* rlon) = 0;


protected:

    /**
     * GridType value
     */
    GridType gridType_;

    /**
     * Hemisphere value
     */
    Hemisphere hemisphere_;

    /**
     * wedge rotation
     */
    short wedgeRotation_;

    /**
     * Mds number
     */
    int mdsNum_;

    /**
     * standard latitude1
     */
    double stanLat1_;

    /**
     * standard latitude2
     */
    double stanLat2_;

    /**
     * base longitude
     */
    double baseLongitude_;

    /**
     * grid incremental constant
     */
    double gridIncConstant_;

    /**
     * max number of rows
     */
    double maxRow_;

    /**
     * max number of columns
     */
    double maxCol_;

    /**
     * upper left latitude
     */
    double upperLeftLatitude_;

    /**
     * upper left longitude
     */
    double upperLeftLongitude_;

    /**
     * lower right latitude
     */
    double lowerRightLatitude_;

    /**
     * lower right longitude
     */
    double lowerRightLongitude_;

    /**
     * upper left X coordinate
     */
    double upperLeftX_;

    /**
     * upper left Y coordinate
     */
    double upperLeftY_;

    /**
     * lower right X coordinate
     */
    double lowerRightX_;

    /**
     * lower right Y coordinate
     */
    double lowerRightY_;

    /**
     * smallest latitude value
     */
    double smallestLatitude_;

    /**
     * largest latitude value
     */
    double biggestLatitude_;

    /**
     * smallest longitude value
     */
    double smallestLongitude_[2];

    /**
     * largest longitude value
     */
    double biggestLongitude_[2];

    /**
     * number of range values
     */
    short numberOfRanges_;

    /**
     * Checks the given latitude against PI/2.  If the magnitude is
     * greater than PI/2, but is within a given epsilon of PI/2,
     * then the lat value is truncated at +-PI/2 (to allow for float-
     * to-double conversion errors).  If its value lies outside the 
     * epsilon range, then an error is returned.
     * 
     * @param rlat   Latitude in radians
     * @return PRO_FAIL if rlat is outside the allowable range or
     *         PRO_SUCCESS if the value is within the allowable range.
     */
    int checkLat(double &rlat);

    /**
     * Checks the given longitude PI.  If the magnitude is greater 
     * than PI, but is within a given epsilon PI, then the long value 
     * is truncated at +-PI (to allow for float-to-double conversion 
     * errors).  If its values lies outside the epsilon range, then 
     * an error is returned.
     * 
     * @param rlon   Longitude in radians
     * @return PRO_FAIL if rlon is outside the allowable range or
     *         PRO_SUCCESS if the value is within the allowable range.
     */
    int checkLon(double &rlon);

private:

};

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getUpperLeftLatitude() const {
    return upperLeftLatitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getUpperLeftLongitude() const {
    return upperLeftLongitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getLowerRightLatitude() const {
    return lowerRightLatitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getLowerRightLongitude() const {
    return lowerRightLongitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getBiggestLatitude() const {
    return biggestLatitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getSmallestLatitude() const {
    return smallestLatitude_;
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getBiggestLongitude() const {
    if (numberOfRanges_ == 1) {
        return biggestLongitude_[0];
    } else {
        return std::max(biggestLongitude_[0], biggestLongitude_[1]);
    }
}

//-----------------------------------------------------------------------------

inline double
VcstMapDataSet::getSmallestLongitude() const {
    if (numberOfRanges_ == 1) {
        return smallestLongitude_[0];
    } else {
        return std::min(smallestLongitude_[0], smallestLongitude_[1]);
    }
}

//-----------------------------------------------------------------------------

inline Hemisphere
VcstMapDataSet::getHemisphere() const {
    return hemisphere_;
}

//-----------------------------------------------------------------------------

inline int
VcstMapDataSet::getGridType() const {
    return gridType_;
}


#endif



