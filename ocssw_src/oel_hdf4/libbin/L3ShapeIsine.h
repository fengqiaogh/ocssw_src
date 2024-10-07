/*
 * File:   L3BinShape.h
 * Author: dshea
 *
 * Created on July 21, 2015, 8:45 AM
 */

#ifndef L3SHAPE_ISINE_H
#define L3SHAPE_ISINE_H

#include <stdint.h>
#include "L3Shape.h"

namespace l3 {

class L3ShapeIsine : public L3Shape {
protected:
    int32_t *numBin; // number of bins in each row
    int64_t *baseBin; // 1st bin of each row
    double *latBin; // center latitude of each row

    int32_t oldRow; // row that was found on last search

public:
    L3ShapeIsine(int32_t numRows);
    virtual ~L3ShapeIsine();

    virtual void constrainRowCol(int32_t &row, int32_t &col) const;
    virtual int64_t getBaseBin(int32_t row) const;
    virtual int32_t getNumCols(int32_t row) const;

    virtual int32_t bin2row(int64_t bin);
    virtual void bin2rowcol(int64_t bin, int32_t &row, int32_t &col);
    virtual int64_t rowcol2bin(int32_t row, int32_t col) const;
    virtual void rowcol2latlon(int32_t row, int32_t col,
            double &lat, double &lon) const;
    virtual int32_t lat2row(double lat) const;
    virtual double row2lat(int32_t row) const;
    virtual void latlon2rowcol(double lat, double lon,
            int32_t &row, int32_t &col) const;
    virtual int64_t latlon2bin(double lat, double lon) const;
    virtual void rowcol2bounds(int32_t row, int32_t col,
            double &north, double &south,
            double &east, double &west) const;
};

}

#endif /* L3SHAPE_ISINE_H */
