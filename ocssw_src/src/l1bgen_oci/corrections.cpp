/**
 * @brief A collection of functions related to L1B corrections
 *
 * @authors Joel Gales (SAIC) Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#include "corrections.hpp"
#include <allocate2d.h>

using namespace std;

int getDarkCorrection(const size_t scanIndex, const uint32_t numScans, const vector<uint8_t> &hamSides,
                      const uint16_t numScansAvg, const uint16_t numPixSkip, const int16_t sciSpatialAgg,
                      const int16_t darkSpatialAgg, const uint32_t numTaps, const vector<int16_t> &specAgg,
                      const uint32_t fillValue, const int16_t numPixAverage, uint32_t ***darkPixels,
                      vector<double> &darkCorrections) {
    // Determine number of bands per tap for hyperspectral data
    vector<int16_t> numBandsPerTap(numTaps);

    if (numTaps == NUMBER_OF_TAPS) {
        // hyperspectral bands
        for (size_t i = 0; i < numTaps; i++)
            if (specAgg[i] > 0)
                numBandsPerTap[i] = 32 / specAgg[i];
            else
                numBandsPerTap[i] = 0;
    } else {
        for (size_t i = 0; i < numTaps; i++)
            numBandsPerTap[i] = 9;
    }

    // Select data for HAM side and determine scan indices
    vector<int32_t> matchingHamSideIndices(numScans);
    int32_t matchingHamSideCount = 0;
    for (size_t i = 0; i < numScans; i++) {
        if (hamSides[i] == hamSides[scanIndex]) {
            matchingHamSideIndices[i] = (int32_t)i;
            matchingHamSideCount++;
        } else {
            matchingHamSideIndices[i] = -1; // Would never be this value otherwise
        }
    }

    int32_t matchedScanIndex = 0;            // HAM sides index for this scan
    for (size_t i = 0; i < numScans; i++) {  // Find the HAM side index for this scan
        if (matchingHamSideIndices[i] == (int32_t)scanIndex) {
            matchedScanIndex = (int32_t)i;
            break;
        }
    }

    // Check for valid dark collect data within specified range
    uint16_t scanRangeSize = numScansAvg;
    bool validDarkFound = false;
    int32_t averageLowerBound = matchedScanIndex, averageUpperBound = matchedScanIndex;

    while (!validDarkFound && scanRangeSize <= matchingHamSideCount) {
        if (scanRangeSize > 1) {
            averageLowerBound = matchedScanIndex - (scanRangeSize / 2);
            averageUpperBound = matchedScanIndex + (scanRangeSize / 2);
            // Check for start or end of granule
            if (averageLowerBound < 0) {
                averageLowerBound = 0;
                averageUpperBound = scanRangeSize - 1;
            }
            if (averageUpperBound >= matchingHamSideCount) {
                averageLowerBound = matchingHamSideCount - scanRangeSize;
                averageUpperBound = matchingHamSideCount - 1;
            }
        }

        // If no valid dark data, expand scan range
        for (size_t i = averageLowerBound; i <= (size_t)averageUpperBound; i++) {
            size_t hamScanLine = matchingHamSideIndices[i];
            if (hamScanLine == -1) // HAM side doesn't match
                continue;

            for (size_t j = numPixSkip; j < (size_t)numPixAverage; j++) {
                if (darkPixels[hamScanLine][0][j] != fillValue) {
                    validDarkFound = true;
                    break;
                }
            }
        }
        if (!validDarkFound)
            scanRangeSize += 2;
    }

    if (!validDarkFound) {
        return -1;
    }

    // Loop through taps and compute dark correction
    int16_t bandIndex = 0;
    for (size_t tapIndex = 0; tapIndex < numTaps; tapIndex++) {
        if (specAgg[tapIndex] <= 0) {
            continue;
        }

        float darkDivisor = 1.0;
        float darkOffset = 0.0;

        const bool needsAggregationAdjustment = (sciSpatialAgg * specAgg[tapIndex] > 4);
        if (needsAggregationAdjustment) {
            darkDivisor = darkSpatialAgg * specAgg[tapIndex] / 4.0;
            darkOffset = (darkDivisor - 1) / (2 * darkDivisor);
        }

        for (int bandInTap = 0; bandInTap < numBandsPerTap[tapIndex]; bandInTap++) {
            size_t currentBand = bandIndex + bandInTap;
            float sumDarkValues = 0.0;
            int numValidPixels = 0;

            for (int tempScan = averageLowerBound; tempScan <= averageUpperBound; tempScan++) {
                size_t hamScanLine = matchingHamSideIndices[tempScan];
                if (hamScanLine == -1) // HAM side doesn't match
                    continue;

                for (int pixelIndex = numPixSkip; pixelIndex < numPixAverage; pixelIndex++) {
                    float pixelValue = darkPixels[hamScanLine][currentBand][pixelIndex];
                    if (pixelValue != fillValue) {
                        sumDarkValues += pixelValue;
                        numValidPixels++;
                    }
                }
            }

            float averageDarkValue = (numValidPixels > 0) ? (sumDarkValues / numValidPixels) : 0;
            darkCorrections[currentBand] = (averageDarkValue / darkDivisor) - darkOffset;
        }
        bandIndex += numBandsPerTap[tapIndex];
    }

    if (scanRangeSize > numScansAvg)
        return 1;

    return 0;
}

int getTempCorrection(uint32_t numInsBands, const Gains &gains, const float *referenceTemps,
                      const double *calibrationTemps, uint32_t numScans, float *temperatureCorrections) {
    const uint16_t numTemperatures = gains.dimensions[TEMP];
    const uint16_t numCoefficients = gains.dimensions[TEMP_CORR];

    // Initialize all temperature corrections to 1.0
    fill(temperatureCorrections, temperatureCorrections + numInsBands, 1.0f);

    // Calculate temperature corrections for each band
    for (size_t tempIndex = 0; tempIndex < numTemperatures; ++tempIndex) {
        float temperatureDifference = calibrationTemps[tempIndex] - referenceTemps[tempIndex];

        for (size_t coeffIndex = 0; coeffIndex < numCoefficients; ++coeffIndex) {
            float temperatureFactor = std::pow(temperatureDifference, coeffIndex + 1);

            for (size_t bandIndex = 0; bandIndex < numInsBands; ++bandIndex) {
                temperatureCorrections[bandIndex] -=
                    gains.k3Coefs[bandIndex][tempIndex][coeffIndex] * temperatureFactor;
            }
        }
    }
    return EXIT_SUCCESS;
}

void getTempCorrection(const size_t numInsBands, const float *referenceTemps,
                       const vector<double> &calibrationTemps, const Gains &gains,
                       vector<double> &temperatureCorrections) {
    // Initialize temperature corrections for each band to 1.0
    for (double &correction : temperatureCorrections) {
        correction = 1.0;
    }
    for (size_t tempIndex = 0; tempIndex < gains.dimensions[TEMP]; ++tempIndex) {
        float temperatureDifference = calibrationTemps[tempIndex] - referenceTemps[tempIndex];

        for (size_t coeffIndex = 0; coeffIndex < gains.dimensions[TEMP_CORR]; ++coeffIndex) {
            float temperatureFactor = std::pow(temperatureDifference, coeffIndex + 1);

            for (size_t bandIndex = 0; bandIndex < numInsBands; ++bandIndex) {
                temperatureCorrections[bandIndex] -=
                    (gains.k3Coefs[bandIndex][tempIndex][coeffIndex] * temperatureFactor);
            }
        }
    }
}

int getRvsCorrection(uint32_t numInsBands, uint16_t numPixels, uint8_t hamSide, const Gains &gains,
                     const double *scanAngles, float **rvsCorrections) {
    const uint16_t numRvsCoefficients = gains.dimensions[3];

    // Initialize all RVS corrections to 1.0
    for (size_t band = 0; band < numInsBands; ++band) {
        std::fill_n(rvsCorrections[band], numPixels, 1.0f);
    }

    // Apply RVS corrections
    for (size_t coeffIndex = 0; coeffIndex < numRvsCoefficients; ++coeffIndex) {
        for (size_t band = 0; band < numInsBands; ++band) {
            const float coefficient = gains.k4Coefs[band][hamSide][coeffIndex];

            for (size_t pixel = 0; pixel < numPixels; ++pixel) {
                const float scanAnglePower = std::pow(scanAngles[pixel], coeffIndex + 1);
                rvsCorrections[band][pixel] += coefficient * scanAnglePower;
            }
        }
    }
    return EXIT_SUCCESS;
}

void getRvsCorrection(const uint32_t numInsBands, const uint16_t numPixels, const uint8_t hamSide,
                      const Gains &gains, const vector<double> scanAngles, vec2D<double> &rvsCorrections) {
    for (vector<double> &row : rvsCorrections) {
        for (double &element : row) {
            element = 1.0;
        }
    }

    const size_t numRvsCoefs = gains.dimensions[3];

    // Apply RVS corrections
    for (size_t coeffIndex = 0; coeffIndex < numRvsCoefs; ++coeffIndex) {
        for (size_t band = 0; band < numInsBands; ++band) {
            const float coefficient = gains.k4Coefs[band][hamSide][coeffIndex];

            for (size_t pixel = 0; pixel < numPixels; ++pixel) {
                const float scanAnglePower = std::pow(scanAngles[pixel], coeffIndex + 1);
                rvsCorrections[band][pixel] += coefficient * scanAnglePower;
            }
        }
    }
}

int getNonlinearityCorrection(uint32_t numInsBands, uint16_t numPixels, uint32_t numNonlinearTerms,
                              const Gains &gains, float **digitalNumbers, float **nonlinearityCorrections) {
    for (size_t band = 0; band < numInsBands; ++band) {
        const vector<double>& k5CoefsForBand = gains.k5Coefs[band];
        float* nonlinearityCorrectionForBand = nonlinearityCorrections[band];
        float* digitalNumbersForBand = digitalNumbers[band];

        for (size_t pixel = 0; pixel < numPixels; ++pixel) {
            // Initialize with the zeroth-order correction
            float correction = k5CoefsForBand[0];
            float dn = digitalNumbersForBand[pixel];
            float dnPower = dn;

            // Add higher-order terms
            for (size_t term = 1; term < numNonlinearTerms; ++term) {
                correction += k5CoefsForBand[term] * dnPower;
                dnPower *= dn;
            }

            nonlinearityCorrectionForBand[pixel] = correction;
        }
    }

    return EXIT_SUCCESS;
}

void getNonlinearityCorrection(const uint32_t numInsBands, const size_t numPixels,
                               const uint32_t numNonlinearTerms, const vec2D<double> &k5Coefs,
                               const vec2D<float> &digitalNumbers, vec2D<double> &nonlinearityCorrections) {
    for (size_t band = 0; band < numInsBands; ++band) {
        for (size_t pixel = 0; pixel < numPixels; ++pixel) {  // TODO: Loop should be in the caller
            // Initialize with the zeroth-order correction
            float correction = k5Coefs[band][0];

            // Add higher-order terms
            for (size_t term = 1; term < numNonlinearTerms; ++term) {
                float dnPower = std::pow(digitalNumbers[band][pixel], term);
                correction += k5Coefs[band][term] * dnPower;
            }

            nonlinearityCorrections[band][pixel] = correction;
        }
    }
}

void getQualityFlags(size_t numPix, size_t numBands, size_t numInsBands, vec2D<float> &digitalNumbers,
                     float **insAggMat, std::vector<uint32_t> &saturationThresholds,
                     boost::multi_array<uint8_t, 2> &qualityFlags) {
    using namespace boost;

    vector<uint8_t> saturationFlags(numInsBands);

    for (size_t i = 0; i < numPix; i++) {
                                           // Check saturation for all instrument bands
        for (size_t j = 0; j < numInsBands; j++) {
            saturationFlags[j] = (digitalNumbers[j][i] >= saturationThresholds[j]) ? 1 : 0;
        }

        // Process all output bands
        for (size_t j = 0; j < numBands; j++) {
            float sum = 0.0;
            for (size_t l = 0; l < numInsBands; l++) {
                sum += insAggMat[l][j] * saturationFlags[l];
            }
            qualityFlags[j][i] = (sum > 0) ? 1 : 0;
        }
    }
}