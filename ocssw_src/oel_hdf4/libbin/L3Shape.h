/*
 * File:   L3BinShape.h
 * Author: dshea
 *
 * Created on July 21, 2015, 8:45 AM
 */

#ifndef L3SHAPE_H
#define L3SHAPE_H

#include <stdint.h>

namespace l3 {

double constrainLat(double lat);
double constrainLon(double lon);

class L3Shape {
protected:
    int64_t totalBins; // total number of bins in the L3 bin shape
    int32_t totalRows; // total number of rows in the L3 bin shape
    double seamLon; // longitude of the start of the shape

public:
    L3Shape(int32_t numRows);
    virtual ~L3Shape();

    virtual int32_t getNumRows() const;
    virtual int64_t getNumBins() const;
    virtual void constrainRow(int32_t &row) const;
    virtual void constrainRowCol(int32_t &row, int32_t &col) const = 0;
    virtual int64_t getBaseBin(int32_t row) const = 0;
    virtual int32_t getNumCols(int32_t row) const = 0;

    virtual int32_t bin2row(int64_t bin) = 0;
    virtual void bin2rowcol(int64_t bin, int32_t &row, int32_t &col) = 0;
    virtual int64_t rowcol2bin(int32_t row, int32_t col) const = 0;
    virtual void rowcol2latlon(int32_t row, int32_t col,
            double &lat, double &lon) const = 0;
    virtual void rowcol2latlon(int32_t row, int32_t col,
            float &lat, float &lon) const;
    virtual void bin2latlon(int64_t bin, double &lat, double &lon);
    virtual void bin2latlon(int64_t bin, float &lat, float &lon);
    virtual int32_t lat2row(double lat) const = 0;
    virtual double row2lat(int32_t row) const = 0;
    virtual void latlon2rowcol(double lat, double lon,
            int32_t &row, int32_t &col) const = 0;
    virtual int64_t latlon2bin(double lat, double lon) const = 0;
    virtual void rowcol2bounds(int32_t row, int32_t col,
            double &north, double &south,
            double &east, double &west) const = 0;
    virtual void rowcol2bounds(int32_t row, int32_t col,
            float &north, float &south,
            float &east, float &west) const;
    virtual void bin2bounds(int64_t bin,
            double &north, double &south,
            double &east, double &west);
    virtual void bin2bounds(int64_t bin,
            float &north, float &south,
            float &east, float &west);
};


}

#endif /* L3SHAPE_H */
