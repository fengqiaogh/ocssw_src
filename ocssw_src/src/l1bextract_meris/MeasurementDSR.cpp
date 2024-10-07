/*
 * MeasurementDSR.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#include "MeasurementDSR.h"

#include "EnvsatUtil.h"

#include <stdio.h>
#include <stdlib.h>

#include <timeutils.h>

MeasurementDSR::MeasurementDSR(int size) :
EnvsatDSR(size) {
}

MeasurementDSR::~MeasurementDSR() {
    // TODO Auto-generated destructor stub
}

void MeasurementDSR::print() {
    EnvsatDSR::print();
    printf("\nMeasurementDSR-----\n");
    printf("startTime = %s\n", unix2merisTime(getStartTime()).c_str());
    printf("quality = %X\n", getQuality());
}

double MeasurementDSR::getStartTime() {
    return getMJD(getBuffer(), 0);
}

unsigned char MeasurementDSR::getQuality() {
    return ((unsigned char*) getBuffer())[12];
}
