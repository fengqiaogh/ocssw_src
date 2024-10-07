/*
 * MeasurementDSR.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#ifndef MEASUREMENTDSR_H_
#define MEASUREMENTDSR_H_

#include "EnvsatDSR.h"

class MeasurementDSR : public EnvsatDSR {
public:
    MeasurementDSR(int size);
    virtual ~MeasurementDSR();

    virtual void print();
    virtual double getStartTime();
    virtual unsigned char getQuality();

};

#endif /* MEASUREMENTDSR_H_ */
