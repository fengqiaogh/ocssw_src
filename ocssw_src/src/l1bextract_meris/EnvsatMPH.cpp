/* 
 * File:   EnvsatMPH.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "EnvsatMPH.h"

#include "EnvsatUtil.h"
#include "EnvsatFile.h"

#include <string.h>
#include <unistd.h>
#include <stdio.h>

using namespace std;

const string EnvsatMPH::TYPE_MERIS_RR_L1 = "MER_RR__1P";
const string EnvsatMPH::TYPE_MERIS_FR_L1 = "MER_FR__1P";
const string EnvsatMPH::TYPE_MERIS_FRS_L1 = "MER_FRS_1P";
const string EnvsatMPH::OBPG_EXTRACT_STRING = "OBPG_EXTRACT=1";
const string EnvsatMPH::START_PIXEL_LABEL_STRING = "START_PIXEL=";
const string EnvsatMPH::PIXEL_UNITS_STRING = "<pixels>";
const string EnvsatMPH::END_PIXEL_LABEL_STRING = "END_PIXEL=";

//EnvsatMPH::EnvsatMPH() {
//    init(NULL);
//}

EnvsatMPH::EnvsatMPH(EnvsatFile* file) {
    init(file);
}

EnvsatMPH::EnvsatMPH(const EnvsatMPH& orig) {
    init(orig.envsatFile);
    memcpy(this->buffer, orig.buffer, MPH_LENGTH);
}

EnvsatMPH::~EnvsatMPH() {
    delete buffer;
}

void EnvsatMPH::init(EnvsatFile* file) {
    buffer = new char[MPH_LENGTH];
    envsatFile = file;
}

// note that the EnvsatFile is not copied

EnvsatMPH& EnvsatMPH::operator=(const EnvsatMPH& src) {
    memcpy(buffer, src.buffer, MPH_LENGTH);
    return *this;
}

int EnvsatMPH::readHeader(int fin) {
    int result;
    int num = 0;

    while (num < (int) MPH_LENGTH) {
        result = read(fin, buffer + num, MPH_LENGTH - num);
        if (result == -1)
            return -1;
        else
            num += result;
    }
    return 0;
}

int EnvsatMPH::writeHeader(int fout) {
    int result;
    int num = 0;

    while (num < (int) MPH_LENGTH) {
        result = write(fout, buffer + num, MPH_LENGTH - num);
        if (result == -1)
            return -1;
        else
            num += result;
    }
    return 0;
}

EnvsatSPH* EnvsatMPH::createSPH() {
    if (getProductType() == TYPE_MERIS_RR_L1
            || getProductType() == TYPE_MERIS_FR_L1
            || getProductType() == TYPE_MERIS_FRS_L1) {
        return new MerisSPH(envsatFile, getSPHSize(), getNumDSDs(),
                getDSDSize());
    }
    return NULL;
}

string& EnvsatMPH::getProductName() {
    static string name;
    getString(buffer, name, PRODUCTNAME_OFFSET, PRODUCTNAME_LENGTH);
    return name;
}

void EnvsatMPH::setProductName(const string &name) {
    setString(name, buffer, PRODUCTNAME_OFFSET, PRODUCTNAME_LENGTH);
}

string& EnvsatMPH::getProductType() {
    static string name;
    getString(buffer, name, PRODUCTNAME_OFFSET, PRODUCT_TYPE_LENGTH);
    return name;
}

bool EnvsatMPH::getObpgExtract() {
    string str;
    getString(buffer, str, OBPG_EXTRACT_OFFSET, OBPG_EXTRACT_LENGTH);
    if (str.compare(OBPG_EXTRACT_STRING) == 0)
        return true;
    else
        return false;
}

void EnvsatMPH::setObpgExtract(bool val) {
    if (val)
        setString(OBPG_EXTRACT_STRING, buffer, OBPG_EXTRACT_OFFSET, OBPG_EXTRACT_LENGTH);
    else
        setString(" ", buffer, OBPG_EXTRACT_OFFSET, OBPG_EXTRACT_LENGTH);
}

/**
 * get the start pixel, if this is an OBPG extracted file.
 *
 * @return start pixel (1 based), or -1 if not set
 */
int EnvsatMPH::getStartPixel() {
    string str;
    int labelSize = START_PIXEL_LABEL_STRING.size();
    getString(buffer, str, START_PIXEL_LABEL_OFFSET, labelSize);
    if (str.compare(START_PIXEL_LABEL_STRING) == 0)
        return getInt(buffer, START_PIXEL_LABEL_OFFSET + labelSize, START_PIXEL_LENGTH);
    else
        return -1;
}

/**
 * set the start pixel (1 based).  if pix<1 is given then the label and value are over
 * written with spaces
 */
void EnvsatMPH::setStartPixel(int pix) {
    int labelSize = START_PIXEL_LABEL_STRING.size();
    int unitsSize = PIXEL_UNITS_STRING.size();
    if (pix < 1) {
        setString(" ", buffer, START_PIXEL_LABEL_OFFSET, labelSize + START_PIXEL_LENGTH + unitsSize);
        return;
    }
    setString(START_PIXEL_LABEL_STRING, buffer, START_PIXEL_LABEL_OFFSET, labelSize);
    setInt(pix, buffer, START_PIXEL_LABEL_OFFSET + labelSize, START_PIXEL_LENGTH);
    setString(PIXEL_UNITS_STRING, buffer, START_PIXEL_LABEL_OFFSET + labelSize + START_PIXEL_LENGTH, unitsSize);
}

/**
 * get the end pixel, if this is an OBPG extracted file.
 *
 * @return end pixel (1 based), or -1 if not set
 */
int EnvsatMPH::getEndPixel() {
    int labelSize = END_PIXEL_LABEL_STRING.size();
    string str;
    getString(buffer, str, END_PIXEL_LABEL_OFFSET, labelSize);
    if (str.compare(END_PIXEL_LABEL_STRING) == 0) {
        return getInt(buffer, END_PIXEL_LABEL_OFFSET + labelSize, END_PIXEL_LENGTH);
    } else
        return -1;
}

/**
 * set the end pixel (1 based).  if pix<1 is given then the label and value are
 * over written with spaces
 */
void EnvsatMPH::setEndPixel(int pix) {
    int labelSize = END_PIXEL_LABEL_STRING.size();
    int unitsSize = PIXEL_UNITS_STRING.size();
    if (pix < 1) {
        setString(" ", buffer, END_PIXEL_LABEL_OFFSET, labelSize + END_PIXEL_LENGTH + unitsSize);
        return;
    }
    setString(END_PIXEL_LABEL_STRING, buffer, END_PIXEL_LABEL_OFFSET, labelSize);
    setInt(pix, buffer, END_PIXEL_LABEL_OFFSET + labelSize, END_PIXEL_LENGTH);
    setString(PIXEL_UNITS_STRING, buffer, END_PIXEL_LABEL_OFFSET + labelSize + END_PIXEL_LENGTH, unitsSize);
}

double EnvsatMPH::getSensingStart() {
    string startTime;
    getString(buffer, startTime, SENSING_START_OFFSET, UTC_DATE_LENGTH);
    return merisTime2unix(startTime);
}

void EnvsatMPH::setSensingStart(double unixTime) {
    setString(unix2merisTime(unixTime), buffer, SENSING_START_OFFSET,
            UTC_DATE_LENGTH);
}

double EnvsatMPH::getSensingStop() {
    string stopTime;
    getString(buffer, stopTime, SENSING_STOP_OFFSET, UTC_DATE_LENGTH);
    return merisTime2unix(stopTime);
}

void EnvsatMPH::setSensingStop(double unixTime) {
    setString(unix2merisTime(unixTime), buffer, SENSING_STOP_OFFSET,
            UTC_DATE_LENGTH);
}

int64_t EnvsatMPH::getTotalSize() {
    return getInt64(buffer, TOT_SIZE_OFFSET, TOT_SIZE_LENGTH);
}

void EnvsatMPH::setTotalSize(int64_t size) {
    setInt64(size, buffer, TOT_SIZE_OFFSET, TOT_SIZE_LENGTH);
}

int EnvsatMPH::getSPHSize() {
    return getInt(buffer, SPH_SIZE_OFFSET, SPH_SIZE_LENGTH);
}

void EnvsatMPH::setSPHSize(int size) {
    setInt(size, buffer, SPH_SIZE_OFFSET, SPH_SIZE_LENGTH);
}

int EnvsatMPH::getNumDSDs() {
    return getInt(buffer, NUM_DSD_OFFSET, NUM_DSD_LENGTH);
}

void EnvsatMPH::setNumDSDs(int num) {
    setInt(num, buffer, NUM_DSD_OFFSET, NUM_DSD_LENGTH);
}

int EnvsatMPH::getDSDSize() {
    return getInt(buffer, DSD_SIZE_OFFSET, DSD_SIZE_LENGTH);
}

void EnvsatMPH::setDSDSize(int size) {
    setInt(size, buffer, DSD_SIZE_OFFSET, DSD_SIZE_LENGTH);
}

int EnvsatMPH::getNumDSs() {
    return getInt(buffer, NUM_DS_OFFSET, NUM_DS_LENGTH);
}

void EnvsatMPH::setNumDSs(int num) {
    setInt(num, buffer, NUM_DS_OFFSET, NUM_DS_LENGTH);
}

void EnvsatMPH::print() {
    printf("EnvsatMPH-----\n");
    printf("MPHSize     = %d\n", getMPHSize());
    printf("productName = %s\n", getProductName().c_str());
    printf("productType = %s\n", getProductType().c_str());
    if (getObpgExtract())
        printf("OBPG Extract= 1\n");
    int i;
    i = getStartPixel();
    if (i != -1)
        printf("startPixel  = %d\n", i);
    i = getEndPixel();
    if (i != -1)
        printf("endPixel    = %d\n", i);
    printf("startTime   = %s\n", unix2merisTime(getSensingStart()).c_str());
    printf("stopTime    = %s\n", unix2merisTime(getSensingStop()).c_str());
    printf("totalSize   = %ld\n", (long) getTotalSize());
    printf("SPHSize     = %d\n", getSPHSize());
    printf("numDSDs     = %d\n", getNumDSDs());
    printf("DSDSize     = %d\n", getDSDSize());
    printf("numDSs      = %d\n", getNumDSs());
}

void EnvsatMPH::printRecursive() {
    print();
}

