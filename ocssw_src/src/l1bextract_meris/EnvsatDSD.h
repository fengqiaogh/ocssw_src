/* 
 * File:   EnvsatDSD.h
 * Author: dshea
 *
 * Created on December 10, 2012, 12:33 PM
 */

#ifndef ENVSATDSD_H
#define ENVSATDSD_H

class EnvsatSPH;

#include <stdint.h>
#include <string>

class EnvsatDSD {
public:
    //EnvsatDSD();
    EnvsatDSD(const EnvsatDSD& orig);
    EnvsatDSD(EnvsatSPH* sph, char* buffer, int index);
    virtual ~EnvsatDSD();

    virtual EnvsatSPH* getSPH() {
        return sph;
    }

    virtual int getIndex() {
        return index;
    }

    virtual bool isValid();

    virtual std::string& getName();
    virtual char getType();

    virtual int64_t getDSOffset();
    virtual void setDSOffset(int64_t offset);

    virtual int64_t getDSSize();
    virtual void setDSSize(int64_t size);

    virtual int getNumDSRs();
    virtual void setNumDSRs(int size);

    virtual int64_t getDSRSize();
    virtual void setDSRSize(int64_t size);

    virtual void seekDataStart(int fd);

    virtual void print();
    virtual void printRecursive();

private:
    void init(EnvsatSPH* sph, char* buffer, int index);

    static const unsigned int DSD_NAME_KEY_OFFSET = 0;
    static const unsigned int DSD_NAME_KEY_LENGTH = 7;
    static const std::string DSD_NAME_KEY_VALUE;

    static const unsigned int DS_NAME_OFFSET = 9;
    static const unsigned int DS_NAME_LENGTH = 28;

    static const unsigned int DS_TYPE_OFFSET = 47;

    static const unsigned int DS_OFFSET_OFFSET = 133;
    static const unsigned int DS_OFFSET_LENGTH = 21;

    static const unsigned int DS_SIZE_OFFSET = 170;
    static const unsigned int DS_SIZE_LENGTH = 21;

    static const unsigned int NUM_DSR_OFFSET = 207;
    static const unsigned int NUM_DSR_LENGTH = 11;

    static const unsigned int DSR_SIZE_OFFSET = 228;
    static const unsigned int DSR_SIZE_LENGTH = 11;

    EnvsatSPH* sph;
    char* buffer; // pointer into SPH buffer
    int index; // index of this DSD

};

#endif /* ENVSATDSD_H */

