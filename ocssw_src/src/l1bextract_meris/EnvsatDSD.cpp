/* 
 * File:   EnvsatDSD.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "EnvsatDSD.h"

#include "EnvsatUtil.h"
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

const string EnvsatDSD::DSD_NAME_KEY_VALUE = "DS_NAME";

EnvsatDSD::EnvsatDSD(const EnvsatDSD& orig) {
    init(orig.sph, orig.buffer, orig.index);
}

EnvsatDSD::EnvsatDSD(EnvsatSPH* sph, char* buffer, int index) {
    init(sph, buffer, index);
}

EnvsatDSD::~EnvsatDSD() {

}

void EnvsatDSD::init(EnvsatSPH* sph, char* buffer, int index) {
    this->sph = sph;
    this->buffer = buffer;
    this->index = index;
}

bool EnvsatDSD::isValid() {
    string name;
    getString(buffer, name, DSD_NAME_KEY_OFFSET, DSD_NAME_KEY_LENGTH);
    return name.compare(DSD_NAME_KEY_VALUE) == 0;
}

string& EnvsatDSD::getName() {
    static string name;
    getString(buffer, name, DS_NAME_OFFSET, DS_NAME_LENGTH);
    return name;
}

char EnvsatDSD::getType() {
    return buffer[DS_TYPE_OFFSET];
}

int64_t EnvsatDSD::getDSOffset() {
    return getInt64(buffer, DS_OFFSET_OFFSET, DS_OFFSET_LENGTH);
}

void EnvsatDSD::setDSOffset(int64_t offset) {
    setInt64(offset, buffer, DS_OFFSET_OFFSET, DS_OFFSET_LENGTH);
}

int64_t EnvsatDSD::getDSSize() {
    return getInt64(buffer, DS_SIZE_OFFSET, DS_SIZE_LENGTH);
}

void EnvsatDSD::setDSSize(int64_t size) {
    setInt64(size, buffer, DS_SIZE_OFFSET, DS_SIZE_LENGTH);
}

int EnvsatDSD::getNumDSRs() {
    return getInt(buffer, NUM_DSR_OFFSET, NUM_DSR_LENGTH);
}

void EnvsatDSD::setNumDSRs(int size) {
    setInt(size, buffer, NUM_DSR_OFFSET, NUM_DSR_LENGTH);
}

int64_t EnvsatDSD::getDSRSize() {
    return getInt64(buffer, DSR_SIZE_OFFSET, DSR_SIZE_LENGTH);
}

void EnvsatDSD::setDSRSize(int64_t size) {
    setInt64(size, buffer, DSR_SIZE_OFFSET, DSR_SIZE_LENGTH);
}

void EnvsatDSD::seekDataStart(int fd) {
    off_t result;
    result = lseek(fd, (off_t) getDSOffset(), SEEK_SET);
    if (result == -1) {
        printf("Error %s:%s:%d   seek to data start failed\n", __FILE__,
                __func__, __LINE__);
        exit(1);
    }
}

void EnvsatDSD::print() {
    printf("EnvsatDSD-----\n");
    if (isValid()) {
        printf("name     = %s\n", getName().c_str());
        printf("index    = %d\n", getIndex());
        printf("Type     = %c\n", getType());
        printf("DSOffset = %ld\n", (long) getDSOffset());
        printf("DSSize   = %ld\n", (long) getDSSize());
        printf("numDSRs  = %d\n", getNumDSRs());
        printf("DSRSize  = %ld\n", (long) getDSRSize());
    } else {
        printf("DSD %d is invalid\n", index);
    }
}

void EnvsatDSD::printRecursive() {
    print();
}
