/*
 * EnvsatTiepoint.cpp
 *
 *  Created on: Jan 4, 2013
 *      Author: dshea
 */

#include "TiepointDSR.h"

#include "EnvsatUtil.h"
#include <stdlib.h>
#include <stdio.h>

#include <timeutils.h>

TiepointDSR::TiepointDSR(int size) :
EnvsatDSR(size) {

    // 12 bytes for time
    // 1 byte for attachment flag
    // 50 bytes need for each pixel
    numPixels = (size - 13) / (50);
    if (numPixels * 50 + 13 != size) {
        printf("Error %s:%s:%d   size does not divide evenly into pixels\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
}

TiepointDSR::~TiepointDSR() {

}

void TiepointDSR::print() {
    EnvsatDSR::print();
    printf("\nTiepointDSR-----\n");
    printf("numPixels = %d\n", numPixels);
    printf("startTime = %s\n", unix2merisTime(getStartTime()).c_str());
}

double TiepointDSR::getStartTime() {
    return getMJD(getBuffer(), 0);
}

double TiepointDSR::getLat(int pixel) {
    int offset = 13 + 4 * pixel;
    return getRawInt32(getBuffer(), offset) * 1e-6;
}

double TiepointDSR::getLon(int pixel) {
    int offset = 13 + 4 * numPixels + 4 * pixel;
    return getRawInt32(getBuffer(), offset) * 1e-6;
}
