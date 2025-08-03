#include <L3ShapeIsine.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <genutils.h>

namespace l3 {


/**
 * constructor for a integerized sine shaped l3 bin file
 * @param numRows number of rows in the bin file
 */
L3ShapeIsine::L3ShapeIsine(int32_t numRows) : L3Shape(numRows) {
    int32_t i;

    oldRow = 0;

    numBin = (int32_t *) calloc(totalRows, sizeof (int32_t));
    baseBin = (int64_t *) calloc(totalRows + 1, sizeof (int64_t));
    latBin = (double *) calloc(totalRows, sizeof (double));

    for (i = 0; i < totalRows; i++) {
        latBin[i] = (i + 0.5) * (180.0 / totalRows) - 90.0;
        numBin[i] = (int32_t) (cos(latBin[i] * OEL_DEGRAD) * (2.0 * totalRows) + 0.5);
    }

    baseBin[0] = 1;

    for (i = 1; i <= totalRows; i++)
        baseBin[i] = baseBin[i - 1] + numBin[i - 1];
    totalBins = baseBin[totalRows] - 1;
}

/**
 * destructor that frees the memory for this class
 */
L3ShapeIsine::~L3ShapeIsine() {
    free(numBin);
    free(baseBin);
    free(latBin);
}

/**
 * keep row and col in the valid range
 * @param row to modify if necessary
 * @param col to modify if necessary
 */
void L3ShapeIsine::constrainRowCol(int32_t &row, int32_t &col) const {
    if (row < 0) {
        row = 0 - row;
        // if pass a pole even times, need to change col to other side
        if((row / totalRows) % 2 == 0) {
            row %= totalRows;
            col += numBin[row] / 2;
        } else {
            row %= totalRows;
        }
    } else if (row >= totalRows) {
        // if pass a pole odd times, need to change col to other side
        if((row / totalRows) % 2 == 1) {
            row %= totalRows;
            col += numBin[row] / 2;
        } else {
            row %= totalRows;
        }
    }

    if (col < 0 || col >= numBin[row]) {
        col %= numBin[row];
    }
}

/**
 * get the bin number of the first column in the given row
 * @param row to look up
 * @return the first bin in this row
 */
int64_t L3ShapeIsine::getBaseBin(int32_t row) const {
    constrainRow(row);
    return baseBin[row];
}

/**
 * get the number of columns in the given row
 * @param row to look up
 * @return the number of cols in this row
 */
int32_t L3ShapeIsine::getNumCols(int32_t row) const {
    constrainRow(row);
    return numBin[row];
}

/**
 * find the row that contains the given bin number
 * @param bin to look up
 * @return row this bin is in
 */
int32_t L3ShapeIsine::bin2row(int64_t bin) {
    int32_t rlow, rhi, rmid;

    if (baseBin[oldRow] <= bin && baseBin[oldRow + 1] > bin) {
        return oldRow;
    } else {
        if (bin < 1)
            return 0; /* south pole */
        if (bin >= totalBins)
            return totalRows - 1; /* north pole */

        /* binary search for row */
        rlow = 0;
        rhi = totalRows - 1;
        while (1) {
            rmid = (rlow + rhi + 1) / 2;
            if (baseBin[rmid] > bin)
                rhi = rmid - 1;
            else
                rlow = rmid;

            if (rlow == rhi) {
                oldRow = rlow;
                return rlow;
            }
        }
    }
    return 0;
}

/**
 * get the row and col for the given bin
 * @param bin to look up
 * @param row that contains the bin
 * @param col for the bin
 */
void L3ShapeIsine::bin2rowcol(int64_t bin, int32_t &row, int32_t &col) {
    if (bin < 1) {
        row = col = 0; /* south pole */
        return;
    }
    if (bin >= totalBins) {
        row = totalRows - 1;
        col = numBin[row] - 1; /* north pole */
        return;
    }
    row = bin2row(bin);
    col = bin - baseBin[row];
}

/**
 * get the bin at the given row/col
 * @param row of the bin
 * @param col of the bin
 * @return bin number at row/col
 */
int64_t L3ShapeIsine::rowcol2bin(int32_t row, int32_t col) const {
    constrainRowCol(row, col);
    return baseBin[row] + col;
}

/**
 * get the center lat/lon for the given row/col
 * @param row of bin
 * @param col of bin
 * @param lat center lat
 * @param lon center lon
 */
void L3ShapeIsine::rowcol2latlon(int32_t row, int32_t col,
        double &lat, double &lon) const {
    constrainRowCol(row, col);
    lat = latBin[row];
    lon = 360.0 * (col + 0.5) / numBin[row] + seamLon;
}

/**
 * get the row for the given lat
 * @param lat to look up
 * @return row of containing bin
 */
int32_t L3ShapeIsine::lat2row(double lat) const {
    lat = constrainLat(lat);
    int32_t row = (int32_t) ((90.0 + lat) * totalRows / 180.0);
    if (row >= totalRows)
        row = totalRows - 1;
    return row;
}

/**
 * get the center latitude of the given row
 * @param row of bin to look up
 * @return center latitude of the row
 */
double L3ShapeIsine::row2lat(int32_t row) const {
    if(row < 0)
        return latBin[0];
    if(row >= totalRows)
        return latBin[totalRows-1];
    return latBin[row];
}

/**
 * get the row/col that contains the given lat/lon
 * @param lat to look up
 * @param lon to look up
 * @param row of containing bin
 * @param col of containing bin
 */
void L3ShapeIsine::latlon2rowcol(double lat, double lon,
        int32_t &row, int32_t &col) const {
    row = lat2row(lat);
    lon = constrainLon(lon);
    col = (int32_t) ((double) numBin[row] * (lon - seamLon) / 360.0);
    if (col >= numBin[row])
        col = numBin[row] - 1;
}

/**
 * get the bin number containing the given lat/lon
 * @param lat to look up
 * @param lon to look up
 * @return  bin number containing lat/lon
 */
int64_t L3ShapeIsine::latlon2bin(double lat, double lon) const {
    int32_t row, col;

    latlon2rowcol(lat, lon, row, col);
    return baseBin[row] + col;
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
void L3ShapeIsine::rowcol2bounds(int32_t row, int32_t col,
        double &north, double &south,
        double &east, double &west) const {
    double lat, lon;
    rowcol2latlon(row, col, lat, lon);
    north = lat + 90.0 / totalRows;
    south = lat - 90.0 / totalRows;
    east = lon + 180.0 / numBin[row];
    west = lon - 180.0 / numBin[row];
}

}
