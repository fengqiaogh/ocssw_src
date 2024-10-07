/*
 * RadianceDSR.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#include "RadianceDSR.h"

#include <stdio.h>
#include <stdlib.h>

RadianceDSR::RadianceDSR(int size) :
MeasurementDSR(size) {
}

RadianceDSR::~RadianceDSR() {
    // TODO Auto-generated destructor stub
}

void RadianceDSR::print() {
    MeasurementDSR::print();
    printf("\nRadianceDSR-----\n");
    printf("no extra data\n");
}

