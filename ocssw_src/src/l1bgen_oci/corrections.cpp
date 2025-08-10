/**
 * @brief A collection of functions related to L1B corrections
 *
 * @authors Joel Gales (SAIC) Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#include "calibrations.hpp"
#include "corrections.hpp"
#include "types.hpp"
#include <allocate2d.h>
#include <gsl/gsl_cblas.h>

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
            matchingHamSideIndices[i] = -1;  // Would never be this value otherwise
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
            int32_t hamScanLine = matchingHamSideIndices[i];
            if (hamScanLine == -1)  // HAM side doesn't match
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
                int32_t hamScanLine = matchingHamSideIndices[tempScan];
                if (hamScanLine == -1)  // HAM side doesn't match
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

double getTempCorrection(const size_t band, const float *referenceTemps, const vector<double> &measuredTemps,
                         const Gains &gains, Device fpa) {
    double correction = 1.0;
    /* The reference and measured temperatures are both shared by the CCDs and the SWIR FPA, and SWIR temps
     start at index 8 */
    size_t tempIndex = (fpa == SWIR) ? 8 : 0;
    size_t numTemps = gains.dimensions[TEMP];           // For data localization
    size_t numTempCoefs = gains.dimensions[TEMP_CORR];  // For data localization

    for (; tempIndex < numTemps; tempIndex++) {
        double tempDiff = measuredTemps[tempIndex] - referenceTemps[tempIndex];
        for (size_t coefIndex = 0; coefIndex < numTempCoefs; coefIndex++) {
            double tempFactor = pow(tempDiff, coefIndex + 1);
            correction -= (gains.k3Coefs[band][tempIndex][coefIndex]) * tempFactor;
        }
    }

    return correction;
}

double getRvsCorrection(const size_t band, const size_t numPixels, const uint8_t hamSide, const Gains &gains,
                        double scanAngle) {
    const size_t numRvsCoefs = gains.dimensions[3];
    double correction = 1.0;
    double scanAnglePower = scanAngle;

    for (size_t i = 0; i < numRvsCoefs; i++) {
        const double coefficient = gains.k4Coefs[band][hamSide][i];
        correction += coefficient * scanAnglePower;
        scanAnglePower *= scanAngle;
    }

    return correction;
}

double getNonlinearityCorrection(const size_t band, const size_t pixel, const uint32_t numNonlinearTerms,
                                 const vec2D<double> &k5Coefs, const vec2D<float> &digitalNumbers) {
    // Initialize with the zeroth-order correction
    float correction = k5Coefs[band][0];

    // Add higher-order terms
    float dnPower = 1;  // Remove std::pow call
    for (size_t term = 1; term < numNonlinearTerms; ++term) {
        dnPower *= digitalNumbers[band][pixel];
        correction += k5Coefs[band][term] * dnPower;
    }

    return correction;
}

vector<uint8_t> saturationFlags;
vector<float> qualityMatrix;
size_t largestDim{0};
void getQualityFlags(size_t numPix, size_t numBands, size_t numInsBands, vec2D<float> &digitalNumbers,
                     float **insAggMat, std::vector<uint32_t> &saturationThresholds,
                     uint8_t *qualityFlags) {
    using namespace boost;

    if (largestDim < max(numInsBands * numPix, numBands * numPix)) {
        largestDim = max(numInsBands * numPix, numBands * numPix);
        saturationFlags.resize(largestDim);
        qualityMatrix.resize(largestDim);
    }

    fill(qualityMatrix.begin(), qualityMatrix.end(), 0);

    for (size_t insBand = 0; insBand < numInsBands; insBand++) {
        for (size_t pix = 0; pix < numPix; pix++) {
            saturationFlags[insBand * numPix + pix] =
                (digitalNumbers[insBand][pix] >= saturationThresholds[insBand]) ? 1 : 0;
        }
    }

    for (size_t insBand = 0; insBand < numInsBands; insBand++) {
        for (size_t l1bBand = 0; l1bBand < numBands; l1bBand++) {
            for (size_t pix = 0; pix < numPix; pix++) {
                qualityMatrix[l1bBand * numPix + pix] +=
                    insAggMat[insBand][l1bBand] * saturationFlags[insBand * numPix + pix];
            }
        }
    }

    for (size_t l1bBand = 0; l1bBand < numBands; l1bBand++) {
        for (size_t pix = 0; pix < numPix; pix++) {
            if(qualityMatrix[l1bBand * numPix + pix] > 0) {
                qualityFlags[l1bBand * numPix + pix] |= 1;
            }
        }
    }
}

void makeXtalkMat(int16_t spatialAgg, CalibrationData &calData, float ***cmat_in, size_t ncpix,
                  size_t nxbands, size_t nbands) {
    // Program to generate the OCI cross-correlation matrix for a blue band
    // and aggregate the coefficients for the instrument configuration

    /*
     * I/O:
     * size_t ncpix : passed in (possibly modified and sent out)
     * size_t nxbands: should be passed in
     * size_t nbands: should be passed in
     * float ***cmat_in : passed in
     * NUM_BLUE_WAVELENGHS: preprocessor-defined macro constant
     * uint32_t numInstrumentBands: passed in via blueCalData
     * float **gainAggMAtrix: passed in via blueCalData
     * double ***cmat: should be passed out
     * int16_t spatialAgg: passed in
     *
     * local:
     * double gmat_tmp[nbands][numInstrumentBands]
     * double gmatu[nbands][numInstrumentBands]
     * double cmat_tmp: [ncpix][numInsBands][numInsBands]
     *
     */

    // Generate blue band correction matrix for aggregated bands
    // This requires a unitized version of the gain aggregation matrix
    // gmat is sized for 512 (NUM_BLUE_WAVELENGTHS) bands,
    // the crosstalk coefficients are for 480 (nbands) bands.

    double gmat_tmp[nbands][calData.numInsBands];
    double gmatu[nbands][calData.numInsBands];
    double cmat_tmp[ncpix][calData.numInsBands][calData.numInsBands];
    string color = determineColor(calData.color);

    cout << "generating " << color << " band crosstalk matrix ..." << endl;

    for (size_t i = 0; i < nbands; i++) {
        for (uint32_t j = 0; j < calData.numInsBands; j++) {
            gmat_tmp[i][j] = calData.gainAgg[NUM_BLUE_WAVELENGTHS - nbands + i][j];
            if (gmat_tmp[i][j] != 0.) {
                gmatu[i][j] = 1.;
            } else {
                gmatu[i][j] = 0;
            }
        }
    }

    // flatten gmat_tmp and transpose(gmatu)
    size_t gmat_tmpOneDSize = nbands * calData.numInsBands;
    size_t gmatuOneDSize = nbands * calData.numInsBands;
    std::vector<double> gmat_tmpOneD(gmat_tmpOneDSize);
    std::vector<double> gmatuOneDT(gmatuOneDSize);
    for (uint32_t i = 0; i < calData.numInsBands; i++) {
        for (size_t n = 0; n < nbands; n++) {
            gmat_tmpOneD[calData.numInsBands * n + i] = gmat_tmp[n][i];
            gmatuOneDT[n + nbands * i] = gmatu[n][i];
        }
    }

    // flatten cmat_in @ influence pixel
    size_t cmat_inOneDSize = nxbands * nbands;
    size_t cmat_tmpOneDSize = calData.numInsBands * calData.numInsBands;
    for (size_t p = 0; p < ncpix; p++) {
        std::vector<double> cmat_inOneD(cmat_inOneDSize);
        for (size_t m = 0; m < nxbands; m++) {
            for (size_t n = 0; n < nbands; n++) {
                cmat_inOneD[m * nbands + n] = static_cast<double>(cmat_in[p][m][n]);
            }
        }
        // calculate flattened cmat_tmp = transpose(gmatu) * cmat_in * gmat_tmp
        // C = transpose(gmatu) * cmat_in
        // D = C * gmat_tmp
        std::vector<double> C(nbands * calData.numInsBands, 0.0);  // Initialize matrix C with zeros
        std::vector<double> D(calData.numInsBands * calData.numInsBands, 0.0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, calData.numInsBands, nbands, nbands, 1.0,
                    gmatuOneDT.data(), nbands, cmat_inOneD.data(), nbands, 0.0, C.data(), nbands);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, calData.numInsBands, calData.numInsBands,
                    nbands, 1.0, C.data(), nbands, gmat_tmpOneD.data(), calData.numInsBands, 0.0, D.data(),
                    calData.numInsBands);
        // form cmat_tmp
        for (size_t k = 0; k < cmat_tmpOneDSize; k++) {
            size_t i = k / calData.numInsBands;
            size_t j = k - i * calData.numInsBands;
            cmat_tmp[p][i][j] = D[k];
        }
    }

    // If spatial aggregation is not 8

    if (spatialAgg != 8) {
        size_t iagf = 8 / spatialAgg;
        size_t ncpixc = ncpix * iagf;
        // allocate calData.mat and assign modified cmat_tmp to it.
        calData.cmat = allocate3d_double(ncpixc, calData.numInsBands, calData.numInsBands);
        for (size_t i = 0; i < ncpixc; ++i) {
            for (size_t j = 0; j < calData.numInsBands; ++j) {
                for (size_t k = 0; k < calData.numInsBands; k++) {
                    for (size_t ic = 0; ic < iagf; ic++) {
                        calData.cmat[iagf * i + ic][j][k] = cmat_tmp[i][j][k] / iagf;
                    }
                }
            }
        }
        // assign ncpixc to calData.ncpix
        calData.ncpix = ncpixc;
    } else {
        // allocate calData.cmat and assign cmat_tmp to it.
        calData.cmat = allocate3d_double(ncpix, calData.numInsBands, calData.numInsBands);
        for (size_t i = 0; i < ncpix; ++i) {
            for (size_t j = 0; j < calData.numInsBands; ++j) {
                for (size_t k = 0; k < calData.numInsBands; k++) {
                    calData.cmat[i][j][k] = cmat_tmp[i][j][k];
                }
            }
        }

        // assign ncpix to calData.ncpix
        calData.ncpix = ncpix;
    }

    cout << "done." << endl;
}

void getXtalkCorrection(uint32_t numInsBands, size_t numCcdPix, size_t ncpix, double ***cmat,
                        const vec2D<float> &digitalNumbers, vec2D<double> &xtalkCorrection) {
    size_t nco2 = ncpix / 2;
    vec2D<float> bandPad = digitalNumbers;

    // create band array with padding

    for (auto &vec : bandPad) {
        vec.resize(vec.size() + 2 * nco2);
        std::rotate(vec.rbegin(), vec.rbegin() + nco2, vec.rend());
    }

    // make flattened matrices
    size_t crossOneDSize = numInsBands * numCcdPix;
    size_t xbandPadOneDSize = numInsBands * numCcdPix;
    size_t cmatOneDSize = numInsBands * numInsBands;
    std::vector<double> crossOneD(crossOneDSize, 0.0);

    for (size_t p = 0; p < ncpix; p++) {
        std::vector<double> cmatOneD(cmatOneDSize);
        std::vector<double> xbandPadOneDT(xbandPadOneDSize);
        std::vector<double> C(crossOneDSize);
        for (size_t j = 0; j < numInsBands; j++) {
            for (size_t k = 0; k < numInsBands; k++) {
                cmatOneD[numInsBands * j + k] = cmat[p][j][k];
            }
            for (size_t q = 0; q < numCcdPix; q++) {
                xbandPadOneDT[j + numInsBands * q] = bandPad[j][q + p];
            }
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numCcdPix, numInsBands, numInsBands, 1.0,
                    xbandPadOneDT.data(), numInsBands, cmatOneD.data(), numInsBands, 0.0, C.data(),
                    numInsBands);
        for (size_t i = 0; i < crossOneDSize; i++) {
            crossOneD[i] += C[i];
        }
    }

    // form xtalkCorrection
    for (size_t k = 0; k < crossOneDSize; k++) {
        size_t ip = k / numInsBands;
        size_t i = k - numInsBands * ip;
        xtalkCorrection[i][ip] = crossOneD[k];
    }
}
