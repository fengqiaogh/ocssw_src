/*
 * EnvsatTiepoint.h
 *
 *  Created on: Jan 4, 2013
 *      Author: dshea
 */

#ifndef ENVSATTIEPOINT_H_
#define ENVSATTIEPOINT_H_

#include "EnvsatDSR.h"

class TiepointDSR : public EnvsatDSR {
public:
    TiepointDSR(int size);
    virtual ~TiepointDSR();

    virtual void print();

    virtual double getStartTime();

    virtual int getNumPixels() {
        return numPixels;
    }

    virtual double getLat(int pixel);
    virtual double getLon(int pixel);

private:
    int numPixels;

};

#endif /* ENVSATTIEPOINT_H_ */
