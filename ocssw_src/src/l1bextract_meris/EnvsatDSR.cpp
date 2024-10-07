/* 
 * File:   EnvsatDSR.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "EnvsatDSR.h"

#include "EnvsatUtil.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <genutils.h>

using namespace std;

EnvsatDSR::EnvsatDSR(const EnvsatDSR& orig) {
    init(orig.size);
}

EnvsatDSR::EnvsatDSR(int size) {
    init(size);
}

EnvsatDSR::~EnvsatDSR() {
    if (buffer)
        free(buffer);
}

void EnvsatDSR::init(int size) {
    this->size = size;
    buffer = (char*) malloc(size);
    arrayOffset = 13;
    pixelSize = 2;
    numPixels = (size - arrayOffset) / pixelSize;
}

void EnvsatDSR::print() {
    printf("EnvsatDSR-----\n");
    printf("size = %d\n", size);
    printf("numPixels = %d\n", numPixels);
}

void EnvsatDSR::printRecursive() {
    print();
}

int EnvsatDSR::readData(int fin) {
    ssize_t result;
    int num = 0;

    while (num < size) {
        result = read(fin, buffer + num, size - num);
        if (result <= 0)
            return -1;
        else
            num += result;
    }
    return 0;
}

int EnvsatDSR::writeData(int fout) {
    ssize_t result;
    int num = 0;

    while (num < size) {
        result = write(fout, buffer + num, size - num);
        if (result <= 0)
            return -1;
        else
            num += result;
    }
    return 0;
}

void EnvsatDSR::setRange(int offset, int count, int val) {
    if (offset < 0)
        return;
    if (offset + count > numPixels)
        count = numPixels - offset;

    // the gcc optimizer crashes if accessing 16 bit int on an odd boundary in
    // some cases, so we need to accessing the array by byte.
    char* ptr = buffer + arrayOffset + offset * pixelSize;
    int16_t val16 = val;
    char* valPtr = (char*) &val16;


    if (!endianess())
        swapc_bytes(valPtr, 2, 1);

    int i;
    for (i = 0; i < count; i++) {
        *ptr = valPtr[0];
        ptr++;
        *ptr = valPtr[1];
        ptr++;
    }

}

