/* 
 * File:   EnvsatSPH.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "EnvsatSPH.h"

#include "EnvsatUtil.h"
#include "EnvsatDSD.h"

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

//EnvsatSPH::EnvsatSPH() {
//    size=0;
//    buffer = NULL;
//    envsatFile = NULL;
//}

EnvsatSPH::EnvsatSPH(const EnvsatSPH& orig) {
    init(orig.envsatFile, orig.size, orig.numDSDs, orig.DSDSize);
    memcpy(this->buffer, orig.buffer, size);
    createDSDs();
}

EnvsatSPH::EnvsatSPH(EnvsatFile* file, int size, int numDSDs, int DSDSize) {
    init(file, size, numDSDs, DSDSize);
}

EnvsatSPH::~EnvsatSPH() {
    delete buffer;
    deleteDSDs();
}

EnvsatSPH& EnvsatSPH::operator=(const EnvsatSPH& src) {
    if (size != src.size) {
        printf("Error %s:%s:%d   size mismatch, this=%d, src=%d\n", __FILE__,
                __func__, __LINE__, size, src.size);
        exit(1);
    }
    memcpy(buffer, src.buffer, size);
    numDSDs = src.numDSDs;
    DSDSize = src.DSDSize;
    deleteDSDs();
    createDSDs();
    return *this;
}

void EnvsatSPH::init(EnvsatFile* file, int size, int numDSDs, int DSDSize) {
    this->size = size;
    buffer = new char[size];
    envsatFile = file;
    this->numDSDs = numDSDs;
    this->DSDSize = DSDSize;
    envsatDSDs = new EnvsatDSD*[numDSDs];
    for (int i = 0; i < numDSDs; i++)
        envsatDSDs[i] = NULL;
}

int EnvsatSPH::readHeader(int fin) {
    int result;
    int num = 0;

    while (num < size) {
        result = read(fin, buffer + num, size - num);
        if (result == -1)
            return -1;
        else
            num += result;
    }

    deleteDSDs();
    createDSDs();
    return 0;
}

int EnvsatSPH::writeHeader(int fout) {
    int result;
    int num = 0;

    while (num < size) {
        result = write(fout, buffer + num, size - num);
        if (result == -1)
            return -1;
        else
            num += result;
    }
    return 0;
}

int EnvsatSPH::getNumDSDs() {
    return numDSDs;
}

EnvsatDSD* EnvsatSPH::getDSD(int index) {
    if (numDSDs == 0)
        return NULL;
    if (index >= numDSDs)
        return NULL;
    return envsatDSDs[index];
}

EnvsatDSD* EnvsatSPH::findFirstDSD(char DSDType) {
    for (int i = 0; i < getNumDSDs(); i++) {
        EnvsatDSD* dsd = getDSD(i);
        if (dsd != NULL)
            if (dsd->getType() == DSDType)
                return dsd;
    }
    return NULL;
}

EnvsatDSD* EnvsatSPH::findDSD(const string& name) {
    for (int i = 0; i < getNumDSDs(); i++) {
        EnvsatDSD* dsd = getDSD(i);
        if (dsd != NULL) {
            if (dsd->getName() == name)
                return dsd;
        }
    }
    return NULL;
}

void EnvsatSPH::createDSDs() {
    if (envsatDSDs == NULL) {
        printf("Error %s:%s:%d   envsatDSDs is NULL\n", __FILE__, __func__,
                __LINE__);
        exit(1);
    }

    int dsdOffset = size - numDSDs * DSDSize;
    for (int i = 0; i < numDSDs; i++) {
        envsatDSDs[i] = new EnvsatDSD(this, buffer + dsdOffset, i);
        dsdOffset += DSDSize;
    }
}

void EnvsatSPH::deleteDSDs() {
    if (envsatDSDs == NULL) {
        printf("Error %s:%s:%d   DSD buffer alignment \n", __FILE__, __func__,
                __LINE__);
        exit(1);
    }

    for (int i = 0; i < numDSDs; i++) {
        if (envsatDSDs[i] != NULL)
            delete envsatDSDs[i];
    }
}

void EnvsatSPH::print() {
    printf("EnvsatSPH-----\n");
    printf("size    = %d\n", getSize());
    printf("numDSDs = %d\n", getNumDSDs());
    printf("DSDSize = %d\n", getDSDSize());
}

void EnvsatSPH::printRecursive() {
    print();
    for (int i = 0; i < numDSDs; i++) {
        if (envsatDSDs[i] != NULL) {
            printf("\n");
            envsatDSDs[i]->printRecursive();
        }
    }
}
