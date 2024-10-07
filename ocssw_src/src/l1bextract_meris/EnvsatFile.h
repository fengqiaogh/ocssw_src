/* 
 * File:   EnvsatFile.h
 * Author: dshea
 *
 * Created on November 28, 2012, 1:07 PM
 */

#ifndef ENVSATFILE_H
#define ENVSATFILE_H

#include "EnvsatMPH.h"
#include "EnvsatSPH.h"
#include "MerisSPH.h"
#include "EnvsatDSD.h"
#include "EnvsatDSR.h"
#include "TiepointDSR.h"
#include "MeasurementDSR.h"
#include "RadianceDSR.h"
#include "FlagsDSR.h"
#include "EnvsatUtil.h"

class EnvsatFile {
public:

    //EnvsatFile();
    EnvsatFile(const EnvsatFile& orig);
    EnvsatFile(const std::string& filename);
    virtual ~EnvsatFile();

    virtual EnvsatFile& operator=(const EnvsatFile& src);

    virtual int openFile(bool write = false);

    /** open file descriptor and read header */
    virtual int readHeader();

    /** open file descriptor and write header */
    virtual int writeHeader();

    virtual const std::string& getFileName();
    virtual void setFileName(const std::string& name);

    virtual EnvsatMPH* getMPH();
    virtual EnvsatSPH* getSPH();

    virtual void print();
    virtual void printRecursive();

    virtual void closeFile();

    virtual int getNumScans();
    virtual int getNumPixels();

    virtual int seekData(EnvsatDSD* dsd, int scanLine = 0);
    virtual int readData(EnvsatDSR* dsr);
    virtual int writeData(EnvsatDSR* dsr);

    virtual void modifyProductName();

private:
    /** name of the file to read and write */
    std::string filename;

    /** low level file descriptor */
    int fd;

    /** Main Product header */
    EnvsatMPH* mph;
    EnvsatSPH* sph;

    void init(const std::string& name);

};

#endif /* ENVSATFILE_H */

