#include <L3Shape.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace l3 {

/**
 * limit latitude values to +-90
 * @param lat latitude to normalize
 * @return lat value clamped to -90 to 90
 */
double constrainLat(double lat) {
    if (lat > 90.0)
        return 90.0;
    if (lat < -90.0)
        return -90.0;
    return lat;
}

/**
 * limit the longitude to +-180
 * @param lon longitude to normalize
 * @return lat value normalized to -180 to 180
 */
double constrainLon(double lon) {
    if(lon > 180.0 || lon <= -180.0) {
        lon = fmod(lon, 360.0);
        if(lon > 180.0)
            lon -= 360.0;
    }
    return lon;
}

/**
 * constructor for a integerized sine shaped l3 bin file
 * @param numRows number of rows in the bin file
 */
L3Shape::L3Shape(int32_t numRows) {
    seamLon = -180.0; /*  this value should be passed in  */
    totalRows = numRows;
    totalBins = 0;
}

/**
 * destructor that frees the memory for this class
 */
L3Shape::~L3Shape() {
}

/**
 * get the number of rows in this shape
 * @return number of rows
 */
int32_t L3Shape::getNumRows() const {
    return totalRows;
}

/**
 * get the total number of bins in this shape
 * @return number of rows
 */
int64_t L3Shape::getNumBins() const {
    return totalBins;
}

/**
 * keep the row within 0 to number of rows
 * @param row to change if necessary
 */
void L3Shape::constrainRow(int32_t &row) const {
    if (row < 0)
        row = 0;
    else if (row >= totalRows)
        row = totalRows - 1;
}

/**
 * get the center lat/lon for the given row/col
 * @param row of bin
 * @param col of bin
 * @param lat center lat
 * @param lon center lon
 */
void L3Shape::rowcol2latlon(int32_t row, int32_t col,
            float &lat, float &lon) const {
    double lat1, lon1;
    rowcol2latlon(row, col, lat1, lon1);
    lat = lat1;
    lon = lon1;
};

/**
 * get the center lat/lon for the given bin number
 * @param bin to look up
 * @param lat center lat
 * @param lon center lon
 */
void L3Shape::bin2latlon(int64_t bin, double &lat, double &lon) {
    int32_t row, col;

    bin2rowcol(bin, row, col);
    rowcol2latlon(row, col, lat, lon);
}

/**
 * get the center lat/lon for the given bin number
 * @param bin to look up
 * @param lat center lat
 * @param lon center lon
 */
void L3Shape::bin2latlon(int64_t bin, float &lat, float &lon) {
    double lat1, lon1;
    bin2latlon(bin, lat1, lon1);
    lat = lat1;
    lon = lon1;
}

/**
 * get the boundaries of the bin at row/col
 * @param row of bin
 * @param col of bin
 * @param north extent of bin
 * @param south extent of bin
 * @param east extent of bin
 * @param west extent of bin
 */
void L3Shape::rowcol2bounds(int32_t row, int32_t col,
            float &north, float &south,
            float &east, float &west) const {
    double n, s, e, w;
    rowcol2bounds(row, col, n, s, e, w);
    north = n;
    south = s;
    east = e;
    west = w;
}


/**
 * get the boundaries of the bin number
 * @param bin to look up
 * @param north extent of bin
 * @param south extent of bin
 * @param east extent of bin
 * @param west extent of bin
 */
void L3Shape::bin2bounds(int64_t bin,
        double &north, double &south,
        double &east, double &west) {
    int32_t row, col;
    bin2rowcol(bin, row, col);
    rowcol2bounds(row, col, north, south, east, west);
}

/**
 * get the boundaries of the bin number
 * @param bin to look up
 * @param north extent of bin
 * @param south extent of bin
 * @param east extent of bin
 * @param west extent of bin
 */
void L3Shape::bin2bounds(int64_t bin,
        float &north, float &south,
        float &east, float &west) {
    double n, s, e, w;
    bin2bounds(bin, n, s, e, w);
    north = n;
    south = s;
    east = e;
    west = w;
}

}
