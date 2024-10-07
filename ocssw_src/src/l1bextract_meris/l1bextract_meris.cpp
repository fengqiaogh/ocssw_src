/*
 * File:   extract_envsat.cpp
 * Author: dshea
 *
 * Created on November 28, 2012, 1:01 PM
 */

#include "EnvsatFile.h"

#include <string>
#include <cstdlib>
#include <stdio.h>

#define BOUNDS_ERROR 110
#undef DEBUG
// #define DEBUG

using namespace std;

/*
 *
 */
int main(int argc, char** argv) {
    string inFilename;
    int spixl = -1;
    int epixl = -1;
    int sscan;
    int escan;
    string outFilename;

    if (argc == 5) {
        inFilename = argv[1];
        sscan = atoi(argv[2]) - 1;
        escan = atoi(argv[3]) - 1;
        outFilename = argv[4];
    } else if (argc == 7) {
        inFilename = argv[1];
        spixl = atoi(argv[2]) - 1;
        epixl = atoi(argv[3]) - 1;
        sscan = atoi(argv[4]) - 1;
        escan = atoi(argv[5]) - 1;
        outFilename = argv[6];
    } else {
        printf("\nusage: l1bextract_meris infile [spixl epixl] sscan escan outfile\n");
        printf("  where:\n");
        printf("    infile  - input envsat file\n");
        printf("    spixl   - start pixel (1-based)\n");
        printf("    epixl   - end pixel (1-based)\n");
        printf("    sscan   - start scan line (1-based)\n");
        printf("    escan   - end scan line (1-based)\n");
        printf("    outfile - output file name\n");
        exit(EXIT_FAILURE);
    }

    if ((spixl > epixl) || (sscan > escan)) {
        printf("Invalid range requested:\n");
        int npixl = epixl - spixl + 1;
        int nscan = escan - sscan + 1;
        printf("spixl: %5d  epixl: %5d  npixl: %5d\n", spixl+1, epixl+1, npixl);
        printf("sscan: %5d  escan: %5d  nscan: %5d\n", sscan+1, escan+1, nscan);
        exit(BOUNDS_ERROR);
    }

    EnvsatFile* envsatFile = new EnvsatFile(inFilename);
    if (envsatFile->openFile()) {
        printf("Error: Could not open file %s\n", inFilename.c_str());
        exit(EXIT_FAILURE);
    }
    if (envsatFile->readHeader()) {
        printf("Error: Could not read header for file %s\n", inFilename.c_str());
        exit(EXIT_FAILURE);
    }

    EnvsatMPH* mph = envsatFile->getMPH();
    if (mph == NULL) {
        printf("Error: Could not get MPH from file\n");
        exit(EXIT_FAILURE);
    }
    EnvsatSPH* sph = envsatFile->getSPH();
    if (sph == NULL) {
        printf("Error: Could not get SPH from file\n");
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    printf("Original File-------------------------------\n");
    envsatFile->printRecursive();
#endif

    EnvsatFile* envsatFile2 = new EnvsatFile(*envsatFile);
    envsatFile2->setFileName(outFilename);

    EnvsatMPH* mph2 = envsatFile2->getMPH();
    if (mph2 == NULL) {
        printf("Error: Could not get MPH from file2\n");
        exit(EXIT_FAILURE);
    }
    EnvsatSPH* sph2 = envsatFile2->getSPH();
    if (sph2 == NULL) {
        printf("Error: Could not get SPH from file2\n");
        exit(EXIT_FAILURE);
    }

    // fix start / end line
    int linesPerTiepoint = sph->getLinesPerTiepoint();
    int origNumScans = envsatFile->getNumScans();
    int origNumPixels = envsatFile->getNumPixels();
    if (sscan < 0)
        sscan = 0;
    if (sscan >= origNumScans) {
        printf("sscan %d >= scan count (%d)\n",
               sscan, origNumScans);
        exit(BOUNDS_ERROR);
    }
    if (escan < 0)
        escan = origNumScans - 1;
    if (escan >= origNumScans) {
        printf("escan %d >= scan count (%d); adjusting.\n",
               escan, origNumScans);
        escan = origNumScans - 1;
    }

    // calc start end tiepoint line
    int sscanTiepoint = sscan / linesPerTiepoint;
    int escanTiepoint = escan / linesPerTiepoint;
    if (escan % linesPerTiepoint > 0)
        escanTiepoint++;

    // make start scan be on the tiepoint boundary
    sscan = sscanTiepoint * linesPerTiepoint;

    // calc quality ADS start and stop lines
    int sscanQuality = sscanTiepoint / 8;
    int escanQuality = escanTiepoint / 8;
    if (escanTiepoint % 8 > 0)
        escanQuality++;

    // calc the number of scans
    int numScans = escan - sscan + 1;
    int numScansTiepoint = escanTiepoint - sscanTiepoint + 1;
    int numScansQuality = escanQuality - sscanQuality + 1;

    // fix the spixl and epixl
    if (spixl < 0)
        spixl = -1;
    if (spixl >= origNumPixels) {
        printf("spixl %d >= pixel count (%d)\n",
               spixl, origNumPixels);
        exit(BOUNDS_ERROR);
    }
    if (epixl < 0)
        epixl = -1;
    if (epixl >= origNumPixels) {
        printf("epixl %d >= pixel count (%d); adjusting.\n",
               epixl, origNumPixels);
        epixl = -1;
        // epixl = origNumPixels - 1;
    }

    // since the epr api that l2gen uses to read the N1 file reverses
    // the pixels on a line, we need to do the same for epixl and spixl
    // so the command line for the extractor and l2gen match
    int tmp_spixl = spixl;
    if (epixl == -1)
        spixl = -1;
    else
        spixl = origNumPixels - epixl - 1;

    if (tmp_spixl == -1)
        epixl = -1;
    else
        epixl = origNumPixels - tmp_spixl - 1;

    // loop through all of the DSDs setting the offset and size
    int64_t offset = mph2->getMPHSize() + mph2->getSPHSize();

    EnvsatDSD* dsd;
    EnvsatDSD* dsd2;

    for (int i = 0; i < sph2->getNumDSDs(); i++) {
        dsd2 = sph2->getDSD(i);
        if (dsd2 != NULL && dsd2->isValid()) {

            if (dsd2->getName() == sph2->getQualityName()) {
                // Quality ADS
                dsd2->setDSOffset(offset);
                dsd2->setNumDSRs(numScansQuality);
                int64_t dsSize = numScansQuality * dsd2->getDSRSize();
                dsd2->setDSSize(dsSize);
                offset += dsSize;

            } else if (dsd2->getName() == sph2->getTiepointName()) {
                // Tiepoint ADS
                dsd2->setDSOffset(offset);
                dsd2->setNumDSRs(numScansTiepoint);
                int64_t dsSize = numScansTiepoint * dsd2->getDSRSize();
                dsd2->setDSSize(dsSize);
                offset += dsSize;

            } else if (dsd2->getType() == 'M') {
                // Measurement DS
                dsd2->setDSOffset(offset);
                dsd2->setNumDSRs(numScans);
                int64_t dsSize = numScans * dsd2->getDSRSize();
                dsd2->setDSSize(dsSize);
                offset += dsSize;

            } else if (dsd2->getType() == 'R') {
                // Reference DS
                // do nothing

            } else {
                // everything else just set the offset to copy the whole data set
                if (dsd2->getDSOffset() != 0) {
                    dsd2->setDSOffset(offset);
                    offset += dsd2->getDSSize();
                }

            }

        } // dsd good

    } // loop DSDs

    // set the total file size
    mph2->setTotalSize(offset);

    // get tie point for first line
    dsd = sph->findDSD(sph->getTiepointName());
    if (dsd == NULL) {
        printf("Error: Could not find Tiepoint DSD\n");
        exit(EXIT_FAILURE);
    }

    TiepointDSR* tiepointDSR1 = new TiepointDSR(dsd->getDSRSize());
    TiepointDSR* tiepointDSR2 = new TiepointDSR(dsd->getDSRSize());

    // set first line metadata
    if (envsatFile->seekData(dsd, sscanTiepoint)) {
        printf("Error: Could seek to DSR#%d of Tiepoint DSD\n", sscanTiepoint);
        exit(EXIT_FAILURE);
    }
    if (envsatFile->readData(tiepointDSR1)) {
        printf("Error: Could read data for DSR#%d of Tiepoint DSD\n", sscanTiepoint);
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    printf("sscan tiepoint# = %d\n", sscanTiepoint);
    tiepointDSR1->print();
#endif

    int samplesPerTiepoint = sph->getSamplesPerTiepoint();
    int startTiepoint, midTiepoint, endTiepoint;

    if (spixl == -1) {
        startTiepoint = 0;
    } else {
        startTiepoint = spixl / samplesPerTiepoint;
    }
    if (epixl == -1) {
        endTiepoint = tiepointDSR1->getNumPixels() - 1;
    } else {
        endTiepoint = epixl / samplesPerTiepoint;
    }
    midTiepoint = (endTiepoint - startTiepoint) / 2 + startTiepoint;

    mph2->setSensingStart(tiepointDSR1->getStartTime());
    sph2->setFirstLineTime(tiepointDSR1->getStartTime());
    sph2->setFirstFirstLat(tiepointDSR1->getLat(startTiepoint));
    sph2->setFirstFirstLon(tiepointDSR1->getLon(startTiepoint));
    sph2->setFirstMidLat(tiepointDSR1->getLat(midTiepoint));
    sph2->setFirstMidLon(tiepointDSR1->getLon(midTiepoint));
    sph2->setFirstLastLat(tiepointDSR1->getLat(endTiepoint));
    sph2->setFirstLastLon(tiepointDSR1->getLon(endTiepoint));

    // set last line metadata
    int tiepointIndex1, tiepointIndex2;
    int tiepointScan1, tiepointScan2;
    tiepointIndex1 = escan / linesPerTiepoint;
    if (tiepointIndex1 >= (dsd->getNumDSRs() - 1)) {
        tiepointIndex1--;
    }
    tiepointIndex2 = tiepointIndex1 + 1;
    tiepointScan1 = tiepointIndex1 * linesPerTiepoint;
    tiepointScan2 = tiepointIndex2 * linesPerTiepoint;

    if (envsatFile->seekData(dsd, tiepointIndex1)) {
        printf("Error: Could not seek file position for DSR#%d of DSD \"%s\"\n", tiepointIndex1, dsd->getName().c_str());
        exit(EXIT_FAILURE);
    }
    if (envsatFile->readData(tiepointDSR1)) {
        printf("Error: Could not read data for DSR#%d of DSD \"%s\"\n", tiepointIndex1, dsd->getName().c_str());
        exit(EXIT_FAILURE);
    }
    if (envsatFile->readData(tiepointDSR2)) {
        printf("Error: Could not read data for DSR#%d of DSD \"%s\"\n", tiepointIndex1 + 1, dsd->getName().c_str());
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    printf("tiepointScan1 = %d\n", tiepointScan1);
    printf("tiepointScan2 = %d\n", tiepointScan2);
    tiepointDSR1->print();
#endif

    double x1, x2, y1, y2;
    double y;
    x1 = tiepointScan1;
    x2 = tiepointScan2;
    y1 = tiepointDSR1->getStartTime();
    y2 = tiepointDSR2->getStartTime();
    y = envsatInterp(x1, x2, y1, y2, escan);
    mph2->setSensingStop(y);
    sph2->setLastLineTime(y);

    y1 = tiepointDSR1->getLat(startTiepoint);
    y2 = tiepointDSR2->getLat(startTiepoint);
    sph2->setLastFirstLat(envsatInterp(x1, x2, y1, y2, escan));

    y1 = tiepointDSR1->getLon(startTiepoint);
    y2 = tiepointDSR2->getLon(startTiepoint);
    sph2->setLastFirstLon(envsatInterp(x1, x2, y1, y2, escan));

    y1 = tiepointDSR1->getLat(midTiepoint);
    y2 = tiepointDSR2->getLat(midTiepoint);
    sph2->setLastMidLat(envsatInterp(x1, x2, y1, y2, escan));

    y1 = tiepointDSR1->getLon(midTiepoint);
    y2 = tiepointDSR2->getLon(midTiepoint);
    sph2->setLastMidLon(envsatInterp(x1, x2, y1, y2, escan));

    y1 = tiepointDSR1->getLat(endTiepoint);
    y2 = tiepointDSR2->getLat(endTiepoint);
    sph2->setLastLastLat(envsatInterp(x1, x2, y1, y2, escan));

    y1 = tiepointDSR1->getLon(endTiepoint);
    y2 = tiepointDSR2->getLon(endTiepoint);
    sph2->setLastLastLon(envsatInterp(x1, x2, y1, y2, escan));

    envsatFile2->modifyProductName();

    // set the pixel extract information in the header
    if (spixl == -1 && epixl == -1)
        mph2->setObpgExtract(false);
    else
        mph2->setObpgExtract(true);
    mph2->setStartPixel(spixl + 1); // if spixl=-1 then 0 also clears this field
    mph2->setEndPixel(epixl + 1); // if epixl=-1 then 0 also clears this field

#ifdef DEBUG
    printf("\n\nModify File-------------------------------\n");
    envsatFile2->printRecursive();
#endif

    // done with header info, write out the new file header
    if (envsatFile2->openFile(true)) {
        printf("Error: Could not open file \"%s\"\n", envsatFile2->getFileName().c_str());
        exit(EXIT_FAILURE);
    }
    if (envsatFile2->writeHeader()) {
        printf("Error: Could not write header for file \"%s\"\n", envsatFile2->getFileName().c_str());
        exit(EXIT_FAILURE);
    }

    // loop through data sets copying data and
    for (int i = 0; i < sph->getNumDSDs(); i++) {

        // bool firstLine = true;

        EnvsatDSR* dsr;
        EnvsatDSD* dsd = sph->getDSD(i);

#ifdef DEBUG
        printf("copying DSD %s\n", dsd->getName().c_str());
#endif

        if (dsd != NULL && dsd->isValid()) {
            if (dsd->getNumDSRs() > 0) {
                int start, end;
                bool extractPixel = false;
                int fillVal = 0;

                if (dsd->getName() == sph->getQualityName()) {
                    // Quality ADS
                    dsr = new EnvsatDSR(dsd->getDSRSize());
                    start = sscanQuality;
                    end = escanQuality;
                } else if (dsd->getName() == sph->getTiepointName()) {
                    // Tiepoint ADS
                    dsr = new EnvsatDSR(dsd->getDSRSize());
                    start = sscanTiepoint;
                    end = escanTiepoint;
                } else if (dsd->getType() == 'M') {
                    // Measurement DS
                    if (dsd->getName().substr(0, 12).compare("Radiance MDS") == 0)
                        dsr = new RadianceDSR(dsd->getDSRSize());
                    else if (dsd->getName().substr(0, 9).compare("Flags MDS") == 0) {
                        dsr = new FlagsDSR(dsd->getDSRSize());
                        fillVal = 128;
                    } else
                        dsr = new MeasurementDSR(dsd->getDSRSize());
                    start = sscan;
                    end = escan;
                    extractPixel = true;
                } else {
                    // copy all DSRs
                    dsr = new EnvsatDSR(dsd->getDSRSize());
                    start = 0;
                    end = dsd->getNumDSRs() - 1;
                }

                if (envsatFile->seekData(dsd, start)) {
                    printf("Error: Could not seek file position for DSR#%d of DSD \"%s\"\n", start, dsd->getName().c_str());
                    exit(EXIT_FAILURE);
                }
                for (int row = start; row <= end; row++) {
                    if (envsatFile->readData(dsr)) {
                        printf("Error: Could not read data for DSR#%d of DSD \"%s\"\n", row, dsd->getName().c_str());
                        exit(EXIT_FAILURE);
                    }

                    if (extractPixel) {
                        if (spixl > 0)
                            dsr->setRange(0, spixl, fillVal);
                        if (epixl >= 0)
                            dsr->setRange(epixl + 1, dsr->getNumPixels() - epixl - 1, fillVal);
                    }

                    if (envsatFile2->writeData(dsr)) {
                        printf("Error: Could not write data for DSR#%d of DSD \"%s\"\n", row, dsd->getName().c_str());
                        exit(EXIT_FAILURE);
                    }
                }
                free(dsr);
            }
        }
    }

    // close the files
    envsatFile2->closeFile();
    envsatFile->closeFile();

    return EXIT_SUCCESS;
}
