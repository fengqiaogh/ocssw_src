/* 
 * File:   EnvsatSPH.h
 * Author: dshea
 *
 * Created on December 10, 2012, 12:33 PM
 */

#ifndef ENVSATSPH_H
#define ENVSATSPH_H

#include <string>

class EnvsatFile;
class EnvsatDSD;

class EnvsatSPH {
public:
    //EnvsatSPH();
    EnvsatSPH(const EnvsatSPH& orig);
    EnvsatSPH(EnvsatFile* file, int size, int numDSDs, int DSDSize);
    virtual ~EnvsatSPH();

    virtual EnvsatSPH& operator=(const EnvsatSPH& src);

    virtual int readHeader(int fin);
    virtual int writeHeader(int fout);

    virtual EnvsatFile* getEnvsatFile() {
        return envsatFile;
    }

    virtual int getSize() {
        return size;
    }

    virtual int getNumDSDs();
    virtual EnvsatDSD* getDSD(int index);
    virtual EnvsatDSD* findFirstDSD(char DSDType);
    virtual EnvsatDSD* findDSD(const std::string& name);

    virtual int getDSDSize() {
        return DSDSize;
    }

    virtual std::string& getDescriptor() = 0;

    virtual double getFirstLineTime() = 0;
    virtual void setFirstLineTime(double unixTime) = 0;

    virtual double getLastLineTime() = 0;
    virtual void setLastLineTime(double unixTime) = 0;

    virtual double getFirstFirstLat() = 0;
    virtual void setFirstFirstLat(double lat) = 0;
    virtual double getFirstFirstLon() = 0;
    virtual void setFirstFirstLon(double lon) = 0;

    virtual double getFirstMidLat() = 0;
    virtual void setFirstMidLat(double lat) = 0;
    virtual double getFirstMidLon() = 0;
    virtual void setFirstMidLon(double lon) = 0;

    virtual double getFirstLastLat() = 0;
    virtual void setFirstLastLat(double lat) = 0;
    virtual double getFirstLastLon() = 0;
    virtual void setFirstLastLon(double lon) = 0;

    virtual double getLastFirstLat() = 0;
    virtual void setLastFirstLat(double lat) = 0;
    virtual double getLastFirstLon() = 0;
    virtual void setLastFirstLon(double lon) = 0;

    virtual double getLastMidLat() = 0;
    virtual void setLastMidLat(double lat) = 0;
    virtual double getLastMidLon() = 0;
    virtual void setLastMidLon(double lon) = 0;

    virtual double getLastLastLat() = 0;
    virtual void setLastLastLat(double lat) = 0;
    virtual double getLastLastLon() = 0;
    virtual void setLastLastLon(double lon) = 0;

    virtual int getSamplesPerLine() = 0;
    virtual int getLinesPerTiepoint() = 0;
    virtual int getSamplesPerTiepoint() = 0;

    virtual const std::string& getQualityName() = 0;
    virtual const std::string& getTiepointName() = 0;

    virtual void print();
    virtual void printRecursive();

protected:
    char* buffer;

    virtual void createDSDs();
    virtual void deleteDSDs();

private:
    int size;
    EnvsatFile* envsatFile;
    int numDSDs;
    int DSDSize;
    EnvsatDSD** envsatDSDs;

    void init(EnvsatFile* file, int size, int numDSDs, int DSDSize);

};

#endif /* ENVSATSPH_H */

