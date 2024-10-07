/*
 * FlagsDSR.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: dshea
 */

#include "FlagsDSR.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

FlagsDSR::FlagsDSR(int size) :
MeasurementDSR(size) {

    // 12 bytes for time
    // 1 byte for quality indicator
    // 1 bytes needed for each pixel quality
    // 2 bytes needed for each pixel Detector Index
    pixelSize = 1;
    numPixels = (size - arrayOffset) / 3;
    if (numPixels * 3 + arrayOffset != size) {
        printf("Error %s:%s:%d   size does not divide evenly into pixels\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
}

FlagsDSR::~FlagsDSR() {
    // TODO Auto-generated destructor stub
}

void FlagsDSR::setRange(int offset, int count, int val) {
    if (offset < 0)
        return;
    if (offset + count > numPixels)
        count = numPixels - offset;

    uint8_t* ptr = (uint8_t*) (getBuffer() + arrayOffset + offset);
    int i;
    for (i = 0; i < count; i++) {
        *ptr = val;
        ptr++;
    }
}

void FlagsDSR::print() {
    MeasurementDSR::print();
    printf("\nFlagsDSR-----\n");
}
