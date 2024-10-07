/*
 * RadianceDSR.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#ifndef RADIANCEDSR_H_
#define RADIANCEDSR_H_

#include "MeasurementDSR.h"

class RadianceDSR : public MeasurementDSR {
public:
    RadianceDSR(int size);
    virtual ~RadianceDSR();

    virtual void print();

};

#endif /* RADIANCEDSR_H_ */
