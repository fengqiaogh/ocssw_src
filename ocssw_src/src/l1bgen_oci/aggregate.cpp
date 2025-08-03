#include "aggregate.hpp"
#include "genutils.h"
#include "calibrations.hpp"

int aggregateBands(const size_t numTaps, size_t *tapAggFactors, size_t binCounts[NUMBER_OF_TAPS],
                   const int16_t aggregationFlags[NUMBER_OF_TAPS], size_t &numInstrumentBands,
                   size_t &numL1bBands) {
    uint32_t numActiveTaps = 0;

    for (size_t i = 0; i < numTaps; i++) {
        if (aggregationFlags[i] > 0) {
            tapAggFactors[numActiveTaps] = i;
            numActiveTaps++;
        }
    }

    if (numActiveTaps == 0) {
        return 2;
    } else {
        for (size_t i = 0; i < NUMBER_OF_TAPS; i++)
            binCounts[i] = 0;
        for (size_t i = 0; i < numActiveTaps; i++)
            binCounts[tapAggFactors[i]] = BANDS_PER_TAP / aggregationFlags[tapAggFactors[i]];
        numInstrumentBands = 0;
        for (size_t i = 0; i < NUMBER_OF_TAPS; i++)
            numInstrumentBands += binCounts[i];
    }

    // Compute number of bands for 8x aggregation with overlapping bands
    numL1bBands = (binCounts[tapAggFactors[0]] * 3) / 4 + 1;
    for (size_t i = 1; i < numActiveTaps; i++) {
        if (aggregationFlags[tapAggFactors[i]] >= aggregationFlags[tapAggFactors[i] - 1])
            numL1bBands += binCounts[tapAggFactors[i]];
        else
            numL1bBands += (binCounts[tapAggFactors[i]] * 3) / 4 + binCounts[tapAggFactors[i] - 1] / 4;
    }

    return 0;
}

int getAggregationMatrices(size_t *taps, int16_t tapAggFactors[NUMBER_OF_TAPS],
                           size_t binCounts[NUMBER_OF_TAPS], uint32_t numInstrumentBands,
                           uint32_t numL1bBands, float **instrumentAggMatrix, float **gainAggMatrix) {
    // tap start/end locations for instrument bands in terms of band indices
    int16_t tapBounds[NUMBER_OF_TAPS][2];

    for (size_t i = 0; i < NUMBER_OF_BANDS; i++) {
        gainAggMatrix[i] = new float[numInstrumentBands];
        for (size_t j = 0; j < numInstrumentBands; j++)
            gainAggMatrix[i][j] = 0.0;
    }

    for (size_t i = 0; i < numInstrumentBands; i++) {
        instrumentAggMatrix[i] = new float[numL1bBands];
        for (size_t j = 0; j < numL1bBands; j++)
            instrumentAggMatrix[i][j] = 0;
    }

    populateGainAggMatrix(gainAggMatrix, tapAggFactors, binCounts, numInstrumentBands, tapBounds);
    populateInstrumentAggMatrix(instrumentAggMatrix, tapAggFactors, binCounts, tapBounds, taps);

    return 0;
}

void populateGainAggMatrix(float **gainAggMatrix, int16_t tapAggFactors[NUMBER_OF_TAPS],
                           size_t binCounts[NUMBER_OF_TAPS], uint32_t numInstrumentBands,
                           int16_t tapBounds[NUMBER_OF_TAPS][2]) {
    int16_t tapBoundBuf = 0;  // Either the start of the current tap or the end of it
    for (size_t currTap = 0; currTap < NUMBER_OF_TAPS; currTap++) {
        tapBounds[currTap][0] = tapBoundBuf;

        if (tapAggFactors[currTap] <= 0) {
            continue;
        }

        for (size_t bin = 0; bin < binCounts[currTap]; bin++) {
            // The index of the first band in this tap
            size_t firstTapBand = BANDS_PER_TAP * currTap;
            // The index of the current bin from the start of the taps
            size_t binIndex = bin * tapAggFactors[currTap];

            for (size_t wavelengthIndex = 0; wavelengthIndex < (size_t)tapAggFactors[currTap];
                 wavelengthIndex++)
                gainAggMatrix[firstTapBand + binIndex + wavelengthIndex][tapBoundBuf + bin] =
                    (1.0 / tapAggFactors[currTap]);
        }
        tapBoundBuf += binCounts[currTap];
        tapBounds[currTap][1] = tapBoundBuf - 1;
    }
}

void fillColInsAggMat(float **instrumentAggMatrix, size_t currTap, size_t currBand,
                      uint16_t endBinsCurrentTap, uint16_t startBinsCurrentTap, int16_t tapBounds[16][2],
                      int16_t tapAggFactors[16], int16_t firstTap) {
    auto previousTap = tapBounds[currTap - 1];
    auto currentTap = tapBounds[currTap];

    for (size_t j = 0; j <= endBinsCurrentTap; j++)
        instrumentAggMatrix[previousTap[1] - endBinsCurrentTap + j][firstTap + currBand] =
            tapAggFactors[currTap - 1] / 8.0;
    for (size_t j = 0; j <= startBinsCurrentTap; j++)
        instrumentAggMatrix[currentTap[0] + j][firstTap + currBand] = tapAggFactors[currTap] / 8.0;
}

void populateInstrumentAggMatrix(float **instrumentAggMatrix, int16_t *tapAggFactors,
                                 size_t binCounts[NUMBER_OF_TAPS], int16_t tapBounds[NUMBER_OF_TAPS][2],
                                 size_t *taps) {
    size_t firstTap = taps[0];

    // First tap
    int16_t currBandIndex = (binCounts[firstTap] * 3) / 4;
    for (int i = 0; i <= currBandIndex; i++) {
        for (size_t j = 0; j <= binCounts[firstTap] / 4 - 1; j++)
            instrumentAggMatrix[i + j][i] = tapAggFactors[firstTap] / 8.0;
    }
    currBandIndex++;

    for (size_t currTap = taps[1]; currTap < NUMBER_OF_TAPS; currTap++) {
        uint16_t numRemainingBands;
        size_t previousBinCount = binCounts[currTap - 1];
        size_t currentBinCount = binCounts[currTap];

        if (currentBinCount <= 0) {
            continue;
        }

        if (previousBinCount <= currentBinCount)  // Transition resolution determined by preceding tap
            numRemainingBands = previousBinCount / 4 - 1;
        else  // Transition resolution determined by current tap
            numRemainingBands = currentBinCount / 4 - 1;

        if (numRemainingBands <= 0)
            continue;

        uint16_t endBinsCurrentTap, startBinsCurrentTap;
        for (size_t currBand = 0; currBand < numRemainingBands; currBand++) {
            if (previousBinCount <= currentBinCount) {
                endBinsCurrentTap = numRemainingBands - currBand - 1;
                startBinsCurrentTap = ((currBand + 1) * currentBinCount) / previousBinCount - 1;
            } else {
                endBinsCurrentTap = ((numRemainingBands - currBand) * previousBinCount) / currentBinCount - 1;
                startBinsCurrentTap = currBand;
            }
            fillColInsAggMat(instrumentAggMatrix, currTap, currBand, endBinsCurrentTap, startBinsCurrentTap,
                             tapBounds, tapAggFactors, currBandIndex);
        }

        currBandIndex += numRemainingBands;

        // Remaining bands using this tap
        for (size_t j = 0; j <= (currentBinCount * 3) / 4; j++) {
            for (size_t k = 0; k < currentBinCount / 4; k++)
                instrumentAggMatrix[tapBounds[currTap][0] + j + k][currBandIndex + j] =
                    tapAggFactors[currTap] / 8.0;
        }
        currBandIndex += (currentBinCount * 3) / 4 + 1;
    }
}

void aggAndCalcRefls(size_t currScan, size_t numBands, size_t numInsBands, const GeoData &geoData,
                     float *preAgg, float *aggregated, float **insAggMat, bool radianceGen,
                     std::vector<double> &solIrrL1a, const float *cosSolZens) {
    for (size_t insBand = 0; insBand < numInsBands; insBand++) {
        for (size_t l1bBand = 0; l1bBand < numBands; l1bBand++) {
            float instrumentAgg = insAggMat[insBand][l1bBand];

            for (size_t pix = 0; pix < geoData.numCcdPix; pix++) {
                size_t dataIndex = insBand * geoData.numCcdPix + pix;
                size_t outIndex = l1bBand * geoData.numCcdPix + pix;

                aggregated[outIndex] += instrumentAgg * preAgg[dataIndex];
            }
        }
    }

    for (size_t l1bBand = 0; l1bBand < numBands; l1bBand++) {
        for (size_t pix = 0; pix < geoData.numCcdPix; pix++) {
            size_t dataIndex = l1bBand * geoData.numCcdPix + pix;

            if (!radianceGen) {
                if (geoData.solarZeniths[currScan * geoData.numCcdPix + pix] < MAX_SOLZ)
                    aggregated[dataIndex] *=
                        OEL_PI * geoData.auCorrection /
                        (solIrrL1a[l1bBand] * cosSolZens[currScan * geoData.numCcdPix + pix]);
                else
                    aggregated[dataIndex] = BAD_FLT;
            }
        }
    }
}