/*
 * File:   L3BinShape.h
 * Author: dshea
 *
 * Created on July 21, 2015, 8:45 AM
 */

#ifndef L3SHAPE_SMI_H
#define L3SHAPE_SMI_H

#include <L3Shape.h>

namespace l3 {

class L3ShapeSMI : public L3Shape {
protected:
    int32_t totalCols;
    double north_, south_, east_, west_;

public:
    L3ShapeSMI(int32_t numRows, int32_t numCols,
            double north, double south,
            double east, double west);
    virtual ~L3ShapeSMI();
    void constrainRowCol(int32_t &row, int32_t &col) const;

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

#endif /* L3SHAPE_SMI_H */
