#include <L3ShapeSMI.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace l3 {

/**
 * constructor for a SMI shaped l3 bin file
 * @param numRows number of rows in the bin file
 */
L3ShapeSMI::L3ShapeSMI(int32_t numRows, int32_t numCols,
        double north, double south,
        double east, double west)
: L3Shape(numRows) {
    totalCols = numCols;
    totalBins = totalRows * totalCols;

    north_ = north;
    south_ = south;
    east_ = east;
    west_ = west;
}

/**
 * destructor that frees the memory for this class
 */
L3ShapeSMI::~L3ShapeSMI() {
}

/**
 * keep row and col in the valid range
 * @param row to modify if necessary
 * @param col to modify if necessary
 */
void L3ShapeSMI::constrainRowCol(int32_t &row, int32_t &col) const {
    if (row < 0)
        row = 0;
    else if (row >= totalRows)
        row = totalRows - 1;

    if (col < 0) {
        col = totalCols + col;
        if(col < 0)
            col = 0;
    } else if (col >= totalCols) {
        col = col - totalCols;
        if(col >= totalCols)
            col = totalCols - 1;
    }
}

/**
 * get the bin number of the first column in the given row
 * @param row to look up
 * @return the first bin in this row
 */
int64_t L3ShapeSMI::getBaseBin(int32_t row) const {
    constrainRow(row);
    return (int64_t) row * (int64_t) totalCols + 1;
}

/**
 * get the number of columns in the given row
 * @param row to look up
 * @return the number of cols in this row
 */
int32_t L3ShapeSMI::getNumCols(int32_t row) const {
    return totalCols;
}

/**
 * find the row that contains the given bin number
 * @param bin to look up
 * @return row this bin is in
 */
int32_t L3ShapeSMI::bin2row(int64_t bin) {
    if (bin < 1)
        return 0; /* south pole */
    if (bin >= totalBins)
        return totalRows - 1; /* north pole */
    return (bin - 1) / totalCols;
}

/**
 * get the row and col for the given bin
 * @param bin to look up
 * @param row that contains the bin
 * @param col for the bin
 */
void L3ShapeSMI::bin2rowcol(int64_t bin, int32_t &row, int32_t &col) {
    if (bin < 1) {
        row = col = 0; /* south pole */
        return;
    }
    if (bin >= totalBins) {
        row = totalRows - 1;
        col = totalCols - 1; /* north pole */
        return;
    }
    bin--;
    row = bin / totalCols;
    col = bin - row * totalCols;
}

/**
 * get the bin at the given row/col
 * @param row of the bin
 * @param col of the bin
 * @return bin number at row/col
 */
int64_t L3ShapeSMI::rowcol2bin(int32_t row, int32_t col) const {
    constrainRowCol(row, col);
    return row * totalCols + col + 1;
}

/**
 * get the center lat/lon for the given row/col
 * @param row of bin
 * @param col of bin
 * @param lat center lat
 * @param lon center lon
 */
void L3ShapeSMI::rowcol2latlon(int32_t row, int32_t col,
        double &lat, double &lon) const {
    constrainRowCol(row, col);
    lat = (north_ - south_) * (row + 0.5) / totalRows + south_;
    lon = (east_ - west_) * (col + 0.5) / totalCols + west_;
}

/**
 * get the row for the given lat
 * @param lat to look up
 * @return row of containing bin
 */
int32_t L3ShapeSMI::lat2row(double lat) const {
    lat = constrainLat(lat);
    if (lat > north_ || lat < south_)
        return -1;

    int32_t row =
            (int32_t) ((lat - south_) * totalRows / (north_ - south_));
    if (row >= totalRows)
        row = totalRows - 1;
    if (row < 0)
        row = 0;
    return row;
}

/**
 * get the center latitude of the given row
 * @param row of bin to look up
 * @return center latitude of the row
 */
double L3ShapeSMI::row2lat(int32_t row) const {
    return (north_ - south_) * (row + 0.5) / totalRows + south_;
}


/**
 * get the row/col that contains the given lat/lon
 * @param lat to look up
 * @param lon to look up
 * @param row of containing bin
 * @param col of containing bin
 */
void L3ShapeSMI::latlon2rowcol(double lat, double lon,
        int32_t &row, int32_t &col) const {
    row = lat2row(lat);
    if (row == -1) {
        col = -1;
        return;
    }
    lon = constrainLon(lon);
    if (lon < west_ || lon > east_) {
        row = col = -1;
        return;
    }
    col = (int32_t) (totalCols * (lon - west_) / (east_ - west_));
    if (col >= totalCols || col < 0)
        row = col = -1;
}

/**
 * get the bin number containing the given lat/lon
 * @param lat to look up
 * @param lon to look up
 * @return  bin number containing lat/lon
 */
int64_t L3ShapeSMI::latlon2bin(double lat, double lon) const {
    int32_t row, col;

    latlon2rowcol(lat, lon, row, col);
    return (int64_t) row * (int64_t) totalCols + (int64_t) col;
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
void L3ShapeSMI::rowcol2bounds(int32_t row, int32_t col,
        double &north, double &south,
        double &east, double &west) const {
    double lat, lon;
    rowcol2latlon(row, col, lat, lon);
    double tmp = (north_ - south_) / totalRows / 2.0;
    north = lat + tmp;
    south = lat - tmp;

    tmp = (east_ - west_) / totalCols / 2.0;
    east = lon + tmp;
    west = lon - tmp;
}

}
