#include "process_science_data.h"
#include <iostream>

using namespace std;

// copy over good l0 data that have a ccd line and band index number over to the science data arr
void copyGoodL0BandDataToSciArr(uint16_t **sciArr, uint16_t **l0BandData, uint16_t numCcdPixels,
                                uint16_t numBands, int16_t *ccdLineNums, int16_t *ccdBandLineIndices) {
    for (size_t i = 0; i < numCcdPixels; i++) {
        int16_t lineNum = ccdLineNums[i];
        int16_t ccdBandIndex = lineNum != -1 ? ccdBandLineIndices[lineNum] : -1;

        // skip if there is no ccdBandIndex because there is no line number
        if (ccdBandIndex == -1) {
            continue;
        }

        for (size_t j = 0; j < numBands; j++) {
            sciArr[j][ccdBandIndex] = l0BandData[i][j];
        }
    }
}

// go through L0 data and id good pixels to copy over to a new array (int16_t version)
void copyGoodL0BandDataToDarkSciArr(uint16_t **darkSciArr, uint16_t **l0BandData, uint16_t numCcdPixels,
                                    uint16_t numDarkCcdPixels, uint16_t numBands, int16_t *ccdLineNums,
                                    int16_t *ccdBandLineIndices) {
    for (size_t i = 0; i < numCcdPixels; i++) {
        int16_t lineNum = ccdLineNums[i];
        int16_t darkLineIndex = lineNum != -1 ? ccdBandLineIndices[lineNum] : -1;
        bool darkLineIndexInRange = darkLineIndex <= numDarkCcdPixels && darkLineIndex >= 0;

        if (darkLineIndex != -1 && darkLineIndexInRange) {
            for (size_t j = 0; j < numBands; j++) {
                darkSciArr[j][darkLineIndex] = l0BandData[i][j];
            }
        }
    }
}

// overloaded version to take uint32_t for SWIR dark science data
void copyGoodSwirL0BandDataToDarkSciArr(uint32_t **darkSciArr, uint32_t **l0BandData, uint16_t numPixels,
                                        uint16_t numDarkPixels, uint16_t numBands, int16_t *ccdLineNums,
                                        int16_t *ccdBandLineIndices) {
    for (size_t i = 0; i < numPixels; i++) {
        int16_t lineNum = ccdLineNums[i];
        int16_t darkLineIndex = lineNum != -1 ? ccdBandLineIndices[lineNum] : -1;
        bool darkLineIndexInRange = darkLineIndex <= numDarkPixels && darkLineIndex >= 0;

        if (darkLineIndex != -1 && darkLineIndexInRange) {
            for (size_t j = 0; j < numBands; j++) {
                darkSciArr[j][darkLineIndex] = l0BandData[i][j];
            }
        }
    }
}

// get the number of pixels that have a bad ccd band line index
int getNumOfInvalidCcdLineIndicies(uint16_t numPixels, int16_t *ccdLineNums, int16_t *ccdBandLineIndices) {
    int numBadIndicies = -1;
    for (size_t i = 0; i < numPixels; i++) {
        int16_t lineNum = ccdLineNums[i];

        // skip if no line number
        if (lineNum == -1) {
            continue;
        }

        // bad band index for line number
        if (ccdBandLineIndices[lineNum] == -1)
            numBadIndicies++;
    }
    return numBadIndicies;
}

// check for ccd line index sequence errors. return 1 if there is one. Otherwise, 0
int isSequenceOutOfOrder(string errorToReport, uint16_t numPixels, int16_t *ccdLineNums,
                         int16_t *ccdBandLineIndices) {
    // check for sequence error by adding the band line index for the first line + ith pixel
    // and comparing it to the current ith pixel band line index.
    int16_t firstLineNum = ccdLineNums[0];

    for (size_t i = 0; i < numPixels; i++) {
        int16_t lineNum = ccdLineNums[i];
        // skip if no line number
        if (lineNum == -1) {
            continue;
        }
        // sequence error, return 1
        if (ccdBandLineIndices[lineNum] != (ccdBandLineIndices[firstLineNum] + int(i))) {
            cout << errorToReport << endl;
            return 1;
        }
    }
    return 0;
}

// count the number of red and blue band pixels and report if any is missing
int checkMissingCcdPixels(std::string errorToReport, bool checkCcdLines, uint16_t numPixels,
                          uint16_t numDarkPixels, int16_t *ccdLineNums, int16_t *ccdBandLineIndices) {
    uint16_t goodPixelCount = 0;
    for (size_t i = 0; i < numPixels; i++) {
        int16_t lineNum = ccdLineNums[i];

        // skip bad lines
        if (lineNum == -1)
            continue;

        if (ccdBandLineIndices[lineNum] > -1)
            goodPixelCount++;
    }

    if ((goodPixelCount < numPixels - numDarkPixels) && checkCcdLines) {
        cout << errorToReport << endl;
        return 2;  // science data status of 2 if missing pixels
    }

    return 0;  // if nothing wrong, or no checks needed, 0
}

// check red, blue and swir for good line numbers and band indicies for all pixels. If it is good, then copy
// it over to a new science data array from l0 array
int checkAndLoadScienceData(short dataType, uint16_t numCcdPixels, uint16_t numSwirPixels,
                            uint16_t numDarkCcdPixels, uint16_t numDarkSwirPixels, uint16_t numBlueBands,
                            uint16_t numRedBands, uint16_t numSwirBands, int16_t *ccdBandLineIndices,
                            int16_t *swirBandLineIndices, int16_t *ccdBandDarkLineIndices,
                            int16_t *swirBandDarkLineIndices, uint16_t **l0BlueBandData,
                            uint16_t **l0RedBandData, uint32_t **l0SwirBandData, int16_t *blueCcdLineNums,
                            int16_t *redCcdLineNums, int16_t *swirLineNums, uint16_t **blueScienceData,
                            uint16_t **redScienceData, uint32_t **swirScienceData, uint16_t **blueCcdDarkData,
                            uint16_t **redCcdDarkData, uint32_t **swirDarkData, int8_t &ccdLineError,
                            int &sciDataStatus) {
    sciDataStatus = 0;

    // CCD line check except in linearity mode
    bool checkCcdLines = (dataType != 5);

    // Check for non null band data and valid line numbers for ccd and swir
    bool hasBlueBandAndCcdLines = ((l0BlueBandData[0] != NULL) && (blueCcdLineNums[0] != -1));
    bool hasRedBandAndCcdLines = ((l0RedBandData[0] != NULL) && (redCcdLineNums[0] != -1));
    bool hasSwirBandAndSwirLines = ((l0SwirBandData[0] != NULL) && (swirLineNums[0] != -1));

    //////////////////// Blue Bands ////////////////////
    int badBlueIndicesCount = -1;

    // Check for invalid line numbers
    if (hasBlueBandAndCcdLines && checkCcdLines) {
        badBlueIndicesCount =
            getNumOfInvalidCcdLineIndicies(numCcdPixels, blueCcdLineNums, ccdBandLineIndices);
        ccdLineError = isSequenceOutOfOrder("Blue CCD line sequence error", numCcdPixels, blueCcdLineNums,
                                            ccdBandLineIndices);
    }

    // Identify pixels and load into output arrays
    if (hasBlueBandAndCcdLines) {
        sciDataStatus = checkMissingCcdPixels("Missing blue band science pixels", checkCcdLines, numCcdPixels,
                                              numDarkCcdPixels, blueCcdLineNums, ccdBandLineIndices);

        // process dark and non-dark blue band data if l0 blue band data is good
        copyGoodL0BandDataToSciArr(blueScienceData, l0BlueBandData, numCcdPixels, numBlueBands,
                                   blueCcdLineNums, ccdBandLineIndices);
        copyGoodL0BandDataToDarkSciArr(blueCcdDarkData, l0BlueBandData, numCcdPixels, numDarkCcdPixels,
                                       numBlueBands, blueCcdLineNums, ccdBandDarkLineIndices);

        // Identify dark count pixels and load into output arrays
        uint16_t nbd = 0;
        for (size_t i = 0; i < numCcdPixels; i++) {
            if (blueCcdLineNums[i] == -1)
                continue;
            if (ccdBandDarkLineIndices[blueCcdLineNums[i]] > -1)
                nbd++;  // LH, 8/24/2020
        }
        if (nbd < numDarkCcdPixels) {
            cout << "Missing blue band dark pixels" << endl;
            sciDataStatus = 2;
        }
    }

    //////////////////// Red Bands ////////////////////
    int badRedIndicesCount = -1;

    // Check for invalid line numbers
    if (hasRedBandAndCcdLines && checkCcdLines) {
        badRedIndicesCount = getNumOfInvalidCcdLineIndicies(numCcdPixels, redCcdLineNums, ccdBandLineIndices);
        ccdLineError = isSequenceOutOfOrder("Red CCD line sequence error", numCcdPixels, redCcdLineNums,
                                            ccdBandLineIndices);
    }

    // Identify pixels and load into output arrays
    if (hasRedBandAndCcdLines) {
        sciDataStatus = checkMissingCcdPixels("Missing red band science pixels", checkCcdLines, numCcdPixels,
                                              numDarkCcdPixels, redCcdLineNums, ccdBandLineIndices);

        // copy over non-dark red band data to science array
        copyGoodL0BandDataToSciArr(redScienceData, l0RedBandData, numCcdPixels, numRedBands, redCcdLineNums,
                                   ccdBandLineIndices);
        copyGoodL0BandDataToDarkSciArr(redCcdDarkData, l0RedBandData, numCcdPixels, numDarkCcdPixels,
                                       numRedBands, redCcdLineNums, ccdBandDarkLineIndices);

        // Identify dark count pixels and load into output arrays
        uint16_t nrd = 0;
        for (size_t i = 0; i < numCcdPixels; i++) {
            if (redCcdLineNums[i] == -1)
                continue;
            if (ccdBandDarkLineIndices[redCcdLineNums[i]] > -1)
                nrd++;
        }

        if (nrd < numDarkCcdPixels) {
            cout << "Missing red band dark pixels" << endl;
            sciDataStatus = 2;
        }
    }

    if (badBlueIndicesCount != -1 || badRedIndicesCount != -1) {
        // Disable SWIR line check for ETU
        cout << "Invalid line numbers" << endl;
        sciDataStatus = 4;
    }

    //////////////////// SWIR Bands ////////////////////

    // Identify pixels and load into output arrays
    // SWIR science data line index has band-by-band offsets, LH, 11/20/2020
    if (hasSwirBandAndSwirLines) {
        uint16_t nss = 0;
        for (size_t i = 0; i < numSwirPixels; i++) {
            if (swirLineNums[i] > -1) {
                nss++;
                for (size_t j = 0; j < numSwirBands; j++) {
                    if (swirBandLineIndices[swirLineNums[i] * numSwirBands + j] > -1) {
                        swirScienceData[j][swirBandLineIndices[swirLineNums[i] * numSwirBands + j]] =
                            l0SwirBandData[i][j];
                    }
                }
            }
        }

        if (nss < numSwirPixels - numDarkSwirPixels) {
            cout << "Missing SWIR band science pixels" << endl;
            sciDataStatus = 2;
        }

        // Identify dark count pixels and load into output arrays
        uint16_t nsd = 0;
        for (size_t i = 0; i < numSwirPixels; i++) {
            if (swirLineNums[i] == -1)
                continue;
            if (swirBandDarkLineIndices[swirLineNums[i]] > -1)
                nsd++;
        }

        copyGoodSwirL0BandDataToDarkSciArr(swirDarkData, l0SwirBandData, numSwirPixels, numDarkSwirPixels,
                                           numSwirBands, swirLineNums, swirBandDarkLineIndices);
    }
    return 0;
}