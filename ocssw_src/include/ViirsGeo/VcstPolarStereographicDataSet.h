/**************************************************************************
 * 
 * NAME: VcstPolarStereographicDataSet
 *
 * DESCRIPTION: See description preceding class declaration below
 *
 *
 **************************************************************************/

#ifndef _VcstPolarStereographicDataSet_H_
#define _VcstPolarStereographicDataSet_H_

//#include <ProCmnDefs.h>
//#include <ProCmnMessage.h>
#include <string>
#include <VcstMapDataSet.h>

/**
 * This class extends the MapDataSet base type for a polar-
 * stereographic projection.  It encapsulates a rectangular area of
 * interest that can be located anywhere on the X-Y plane of the 
 * map.  The implementation is based on two coordinante systems.
 * The first is a standard Cartesian coordinate plane.  The origin
 * of the X,Y axes represents the north or south pole (depending on
 * the hemisphere used), with the X-axis along the base longitude
 * of the MDS (positive to the right) and the Y-axis completing a 
 * right-handed coordinate system (positive upward).  Any latitude
 * line forms a circle around the pole (origin), with the equator
 * being the outmost such circle.  Any longitude line is a straight
 * line from the pole to a point on the equatorial circle.  
 * Eastward motion is counterclockwise (for the northern hemisphere)
 * or clockwise (for the southern hemisphere) along any cirlce of
 * latitude.  The second system consists of the row-column 
 * coordinate plane of the MDS rectangle.  The center of the upper-
 * left pixel is at the origin (row 0, column 0) for the 
 * convenience of other software already in use.  WIth row #'s 
 * increasing downward and cloumn #'s increasing to the right, this
 * is actually another right-handed Cartesian coordinate system in 
 * the map plane.  The common relationship between the two coordinate
 * systems is that the X-axis is parallel to the top of the MDS 
 * rectangle.  Also note that it does not matter whether the long  
 * side of the rectangle is considered left-right or top-down. 
 */
class VcstPolarStereographicDataSet : public VcstMapDataSet {
public:

    /**
     * Constructor - This constructs a polar stereographic MDS given the 
     * hemisphere for which to construct it, the latitude values of the corners
     * in radians, and the number of pixels in the up-down and left-right
     * directions.
     *
     * @param theHemisphere the hemisphere of the dataset
     * @param upperLeftLatitude the upper left latitude value
     * @param upperLeftLongitude the upper left longitude value
     * @param lowerRightLatitude the lower right latitude value
     * @param lowerRightLongitude the lower right longitude value
     * @param heightInPixels the pixal height of the area
     * @param widthInPixels the pixal width of the area
     */
    VcstPolarStereographicDataSet(
            Hemisphere theHemisphere,
            double upperLeftLatitude,
            double upperLeftLongitude,
            double lowerRightLatitude,
            double lowerRightLongitude,
            int heightInPixels,
            int widthInPixels
            );


    /**
     * Converts the input grid position to a latitude and longitude
     *
     * @param row Grid row
     * @param col Grid column
     * @param rlat Latitude of row/col
     * @param rlon Longitude of row/col
     * @retval PRO_SUCCESS or PRO_FAIL
     */
    /* virtual */
    int gridToLatLon(
            double row,
            double col,
            double* rlat,
            double* rlon
            );

    /**
     * Converts the input latitude and longitude to a grid position
     *
     * @param rlat Latitude 
     * @param rlon Longitude
     * @param row Grid row of latitude/longitude (rlat/rlon)
     * @param col Grid column of latitude/longitude (rlat/rlon)
     * @retval PRO_SUCCESS or PRO_FAIL
     */
    /* virtual */
    int latLonToGrid(
            double rlat,
            double rlon,
            double* row,
            double* column
            );



private:

    /**
     * C style message string to hold error messages and other information
     */
    //    std::string msgStr_;
    char msgStr_[255];

    /**
     * Message container used to give status via msgStr
     */
    //    ProCmnMessage* message_;

    /**
     * State holders for the DEBUG message throttle.  Protect against
     * flooding the logs with DEBUG messages during a maneuver.
     */
    int dbgMsgThrottle1_;
    int dbgMsgThrottle2_;
    int dbgMsgThrottle3_;
    int dbgMsgThrottle4_;

    /**
     * Method returns true if this is a northern hemisphere MDS.
     *
     * @retval TRUE or FALSE
     */
    bool isNorthernHemis();

    /**
     * Method finds the lat/lon ranges in the MDS
     *
     * @retval PRO_SUCCESS or PRO_FAIL
     */
    int findLatLonRanges();


};

#endif
