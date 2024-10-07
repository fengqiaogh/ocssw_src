/* 
 * File:   EnvsatMPH.h
 * Author: dshea
 *
 * Created on November 28, 2012, 1:11 PM
 */

#ifndef ENVSATMPH_H
#define ENVSATMPH_H

class EnvsatFile;
class EnvsatSPH;

#include <stdint.h>
#include <string>

class EnvsatMPH {
public:

    // Product type constants
    static const std::string TYPE_MERIS_RR_L1;
    static const std::string TYPE_MERIS_FR_L1;
    static const std::string TYPE_MERIS_FRS_L1;

    //EnvsatMPH();
    EnvsatMPH(const EnvsatMPH& orig);
    EnvsatMPH(EnvsatFile* file);
    virtual ~EnvsatMPH();

    virtual EnvsatMPH& operator=(const EnvsatMPH& src);

    virtual EnvsatFile* getEnvsatFile() {
        return envsatFile;
    }

    virtual int readHeader(int fin);
    virtual int writeHeader(int fout);

    // function to create the correct SPH object
    virtual EnvsatSPH* createSPH();

    virtual int getMPHSize() {
        return MPH_LENGTH;
    }

    virtual std::string& getProductName();
    virtual void setProductName(const std::string& name);

    /** return the first part of the product name which identifies the type
     * of the product.  For example:
     * MER_FR__1PNUPA20030723_105132_000000982018_00223_07291_0388.N1
     * would return MER_FR__1P */
    virtual std::string& getProductType();

    virtual bool getObpgExtract();
    virtual void setObpgExtract(bool val);

    virtual int getStartPixel();
    virtual void setStartPixel(int pix);
    virtual int getEndPixel();
    virtual void setEndPixel(int pix);

    virtual double getSensingStart();
    virtual void setSensingStart(double unixTime);
    virtual double getSensingStop();
    virtual void setSensingStop(double unixTime);

    /** return the size of the whole produce file in bytes*/
    virtual int64_t getTotalSize();
    virtual void setTotalSize(int64_t size);

    virtual int getSPHSize();
    virtual void setSPHSize(int size);

    virtual int getNumDSDs();
    virtual void setNumDSDs(int num);
    virtual int getDSDSize();
    virtual void setDSDSize(int size);

    virtual int getNumDSs();
    virtual void setNumDSs(int num);

    virtual void print();
    virtual void printRecursive();

protected:

private:
    // constants needed to index into the buffer
    static const unsigned int MPH_LENGTH = 1247;

    static const unsigned int PRODUCTNAME_OFFSET = 9;
    static const unsigned int PRODUCTNAME_LENGTH = 62;
    static const unsigned int PRODUCT_TYPE_LENGTH = 10;

    static const unsigned int OBPG_EXTRACT_OFFSET = 120;
    static const unsigned int OBPG_EXTRACT_LENGTH = 14;
    static const std::string OBPG_EXTRACT_STRING;

    static const unsigned int SOFTWARE_VER_OFFSET = 279;
    static const unsigned int SOFTWARE_VER_LENGTH = 14;

    // begin 296
    // START_PIXEL=+1234567890
    static const unsigned int START_PIXEL_LABEL_OFFSET = 295;
    static const std::string START_PIXEL_LABEL_STRING;
    static const std::string PIXEL_UNITS_STRING;
    static const unsigned int START_PIXEL_LENGTH = 11;

    static const unsigned int SENSING_START_OFFSET = 351;
    static const unsigned int SENSING_STOP_OFFSET = 394;

    // begin 424
    // END_PIXEL=+1234567890
    static const unsigned int END_PIXEL_LABEL_OFFSET = 423;
    static const std::string END_PIXEL_LABEL_STRING;
    static const unsigned int END_PIXEL_LENGTH = 11;

    static const unsigned int TOT_SIZE_OFFSET = 1075;
    static const unsigned int TOT_SIZE_LENGTH = 21;

    static const unsigned int SPH_SIZE_OFFSET = 1113;
    static const unsigned int SPH_SIZE_LENGTH = 11;

    static const unsigned int NUM_DSD_OFFSET = 1140;
    static const unsigned int NUM_DSD_LENGTH = 11;

    static const unsigned int DSD_SIZE_OFFSET = 1161;
    static const unsigned int DSD_SIZE_LENGTH = 11;

    static const unsigned int NUM_DS_OFFSET = 1194;
    static const unsigned int NUM_DS_LENGTH = 11;

    // init all of the standard internal structures
    void init(EnvsatFile* file);

    // member variables
    EnvsatFile* envsatFile;
    char* buffer;

};

#endif /* ENVSATMPH_H */

