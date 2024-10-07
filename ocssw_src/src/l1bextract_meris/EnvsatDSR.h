/* 
 * File:   EnvsatDSR.h
 * Author: dshea
 *
 * Created on December 10, 2012, 12:33 PM
 */

#ifndef ENVSATDSR_H
#define ENVSATDSR_H

#include <string>

class EnvsatDSR {
public:
    EnvsatDSR(const EnvsatDSR& orig);
    EnvsatDSR(int size);

    virtual ~EnvsatDSR();

    virtual int getSize() {
        return size;
    }

    virtual char* getBuffer() {
        return buffer;
    }

    virtual int getNumPixels() {
        return numPixels;
    }

    virtual void print();
    virtual void printRecursive();

    virtual int readData(int fin);
    virtual int writeData(int fout);

    virtual void setRange(int offset, int count, int val);

protected:
    // num bytes into buffer that the pixel array starts
    int arrayOffset;
    // size of pixel data in bytes
    int pixelSize;
    // number of pixels in the array
    int numPixels;

private:
    void init(int size);

    // size of the buffer in bytes
    int size;
    // memory for the pixel data
    char* buffer;

};

#endif /* ENVSATDSR_H */

