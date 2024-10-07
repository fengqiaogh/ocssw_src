/* 
 * File:   MerisSPH.h
 * Author: dshea
 *
 * Created on December 10, 2012, 12:33 PM
 */

#ifndef MERIS_SPH_H
#define MERIS_SPH_H

#include "EnvsatSPH.h"

#include <string>

class MerisSPH : public EnvsatSPH {
public:
    MerisSPH(const MerisSPH& orig);
    MerisSPH(EnvsatFile* file, int size, int numDSDs, int DSDSize);
    //    virtual ~MerisSPH();

    virtual std::string& getDescriptor();

    virtual double getFirstLineTime();
    virtual void setFirstLineTime(double unixTime);

    virtual double getLastLineTime();
    virtual void setLastLineTime(double unixTime);

    virtual double getFirstFirstLat();
    virtual void setFirstFirstLat(double lat);
    virtual double getFirstFirstLon();
    virtual void setFirstFirstLon(double lon);

    virtual double getFirstMidLat();
    virtual void setFirstMidLat(double lat);
    virtual double getFirstMidLon();
    virtual void setFirstMidLon(double lon);

    virtual double getFirstLastLat();
    virtual void setFirstLastLat(double lat);
    virtual double getFirstLastLon();
    virtual void setFirstLastLon(double lon);

    virtual double getLastFirstLat();
    virtual void setLastFirstLat(double lat);
    virtual double getLastFirstLon();
    virtual void setLastFirstLon(double lon);

    virtual double getLastMidLat();
    virtual void setLastMidLat(double lat);
    virtual double getLastMidLon();
    virtual void setLastMidLon(double lon);

    virtual double getLastLastLat();
    virtual void setLastLastLat(double lat);
    virtual double getLastLastLon();
    virtual void setLastLastLon(double lon);

    virtual int getSamplesPerLine();
    virtual int getLinesPerTiepoint();
    virtual int getSamplesPerTiepoint();

    virtual const std::string& getQualityName();
    virtual const std::string& getTiepointName();

    virtual void print();

private:
    static const unsigned int DESCRIPTOR_OFFSET = 16;
    static const unsigned int DESCRIPTOR_LENGTH = 28;

    static const unsigned int FIRST_LINE_TIME_OFFSET = 135;
    static const unsigned int LAST_LINE_TIME_OFFSET = 180;

    static const unsigned int FIRST_FIRST_LAT_OFFSET = 225;
    static const unsigned int FIRST_FIRST_LON_OFFSET = 264;
    static const unsigned int FIRST_MID_LAT_OFFSET = 300;
    static const unsigned int FIRST_MID_LON_OFFSET = 337;
    static const unsigned int FIRST_LAST_LAT_OFFSET = 374;
    static const unsigned int FIRST_LAST_LON_OFFSET = 412;

    static const unsigned int LAST_FIRST_LAT_OFFSET = 449;
    static const unsigned int LAST_FIRST_LON_OFFSET = 487;
    static const unsigned int LAST_MID_LAT_OFFSET = 522;
    static const unsigned int LAST_MID_LON_OFFSET = 558;
    static const unsigned int LAST_LAST_LAT_OFFSET = 594;
    static const unsigned int LAST_LAST_LON_OFFSET = 631;

    static const unsigned int SAMPLES_PER_LINE_OFFSET = 1404;
    static const unsigned int SAMPLES_PER_LINE_LENGTH = 6;
    static const unsigned int LINES_PER_TIEPOINT_OFFSET = 1437;
    static const unsigned int LINES_PER_TIEPOINT_LENGTH = 4;
    static const unsigned int SAMPLES_PER_TIEPOINT_OFFSET = 1461;
    static const unsigned int SAMPLES_PER_TIEPOINT_LENGTH = 4;

    static const unsigned int FIRST_DSD_OFFSET = 1542;

};

#endif /* MERIS_SPH_H */

