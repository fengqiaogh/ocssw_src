/*
 * FlagsDSR.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#ifndef FLAGSDSR_H_
#define FLAGSDSR_H_

#include "MeasurementDSR.h"

class FlagsDSR : public MeasurementDSR {
public:
    FlagsDSR(int size);
    virtual ~FlagsDSR();

    virtual void setRange(int offset, int count, int val);

    virtual void print();

};

#endif /* FLAGSDSR_H_ */
