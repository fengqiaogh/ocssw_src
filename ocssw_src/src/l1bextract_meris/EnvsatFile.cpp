/* 
 * File:   EnvsatFile.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:07 PM
 */

#include "EnvsatFile.h"

#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <timeutils.h>

using namespace std;

//EnvsatFile::EnvsatFile() {
//    init("");
//}

EnvsatFile::EnvsatFile(const EnvsatFile& orig) {
    init(orig.filename);
    mph = new EnvsatMPH(*orig.mph);
    sph = orig.mph->createSPH();
    *sph = *(orig.sph);
}

EnvsatFile::EnvsatFile(const string& filenameStr) {
    init(filenameStr);
    mph = new EnvsatMPH(this);
}

EnvsatFile::~EnvsatFile() {
    if (fd != -1)
        close(fd);
    if (sph)
        delete sph;
    if (mph)
        delete mph;
}

void EnvsatFile::init(const string& name) {
    filename = name;
    fd = -1;
    mph = NULL;
    sph = NULL;
}

EnvsatFile& EnvsatFile::operator=(const EnvsatFile& src) {
    filename = src.filename;
    fd = -1;
    if (mph)
        *mph = *(src.mph);
    else
        mph = new EnvsatMPH(*(src.mph));

    if (sph) {
        free(sph);
        sph = NULL;
    }

    if (src.sph) {
        sph = mph->createSPH();
        *sph = *(src.sph);
    }
    return *this;
}

int EnvsatFile::openFile(bool write) {
    if (fd != -1)
        close(fd);
    if (write)
        fd = creat(filename.c_str(), 0666);
    else
        fd = open(filename.c_str(), O_RDONLY);
    if (fd == -1)
        return -1;
    return 0;
}

int EnvsatFile::readHeader() {
    int result = mph->readHeader(fd);
    if (result != 0)
        return result;
    sph = mph->createSPH();
    if (sph == NULL) {
        printf("Error %s:%s:%d   Could not create SPH\n", __FILE__, __func__,
                __LINE__);
        exit(1);
    }
    return sph->readHeader(fd);
}

int EnvsatFile::writeHeader() {
    int result = mph->writeHeader(fd);
    if (result != 0) {
        return result;
    }
    if (sph == NULL) {
        printf("Error %s:%s:%d   Could not create SPH\n", __FILE__, __func__,
                __LINE__);
        exit(1);
    }
    return sph->writeHeader(fd);
}

const string& EnvsatFile::getFileName() {
    return filename;
}

void EnvsatFile::setFileName(const string& name) {
    filename = name;
}

EnvsatMPH* EnvsatFile::getMPH() {
    return mph;
}

EnvsatSPH* EnvsatFile::getSPH() {
    return sph;
}

void EnvsatFile::print() {
    printf("EnvsatFile-----\n");
    printf("File Name = %s\n", getFileName().c_str());
}

void EnvsatFile::printRecursive() {
    print();
    printf("\n");
    mph->printRecursive();
    if (sph != NULL) {
        printf("\n");
        sph->printRecursive();
    }
}

void EnvsatFile::closeFile() {
    if (fd != -1) {
        close(fd);
        fd = -1;
    }
}

int EnvsatFile::getNumScans() {
    if (sph == NULL) {
        printf("Error %s:%s:%d   SPH is NULL\n", __FILE__, __func__, __LINE__);
        exit(1);
    }

    EnvsatDSD* dsd = sph->findFirstDSD('M');
    if (dsd == NULL) {
        printf("Error %s:%s:%d   Could not find first Measurement DSD\n",
                __FILE__, __func__, __LINE__);
        exit(1);
    }
    return dsd->getNumDSRs();
}

int EnvsatFile::getNumPixels() {
    if (sph == NULL) {
        printf("Error %s:%s:%d   SPH is NULL\n", __FILE__, __func__, __LINE__);
        exit(1);
    }
    return sph->getSamplesPerLine();
}

int EnvsatFile::seekData(EnvsatDSD* dsd, int scanLine) {
    off_t offset = dsd->getDSOffset() + dsd->getDSRSize() * scanLine;
    off_t result = lseek(fd, offset, SEEK_SET);
    if (result == -1)
        return -1;
    return 0;
}

int EnvsatFile::readData(EnvsatDSR* dsr) {
    return dsr->readData(fd);
}

int EnvsatFile::writeData(EnvsatDSR* dsr) {
    return dsr->writeData(fd);
}

void EnvsatFile::modifyProductName() {
    if (mph == NULL) {
        printf("Error %s:%s:%d   mph == NULL\n", __FILE__, __func__, __LINE__);
        exit(1);
    }
    if (sph == NULL) {
        printf("Error %s:%s:%d   sph == NULL\n", __FILE__, __func__, __LINE__);
        exit(1);
    }

    char str[32];
    int16_t year, month, day, hour, min;
    double sec;

    string oldName = mph->getProductName();
    string newName = oldName.substr(0, 14);
    unix2ymdhms((double) sph->getFirstLineTime(), &year, &month, &day, &hour, &min,
            &sec);
    sprintf(str, "%04d%02d%02d_%02d%02d%02d_", year, month, day, hour, min,
            (int) sec);
    newName += str;
    int deltaT = abs((int) (sph->getLastLineTime() - sph->getFirstLineTime()));
    sprintf(str, "%08d", deltaT);
    newName += str;
    newName += oldName.substr(38, 17);
    newName += "0000.N1";
    mph->setProductName(newName);
}
