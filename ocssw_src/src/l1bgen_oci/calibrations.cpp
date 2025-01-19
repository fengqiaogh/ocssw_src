/**
 * @file
 * @brief OCI calibration and gain calculation functions
 *
 * This file contains functions for reading OCI calibration lookup tables,
 * calculating gains, and retrieving calibration temperatures.
 *
 * @author Jakob C Lindo (SSAI)
 * @date Aug 2024
 */

#include "calibrations.hpp"
#include <timeutils.h>
#include "corrections.hpp"
#include <genutils.h>
#include "geo_data.hpp"
#include "aggregate.hpp"
#include "device.hpp"

#include <stdexcept>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

const double J2000 = 2451545;

using namespace std;
using namespace netCDF;

template <typename T>
using matrix2D = boost::multi_array<T, 2>;

static vector<double> darkCorrections;
static vec2D<float> digitalNumbers;
static float **preAggCalData = nullptr;
static size_t mostInsBands;
vector<double> tempCorrections;
vec2D<double> k4;
vec2D<double> k5;
matrix2D<uint8_t> qualityFlags(boost::extents[0][0]);

template <typename T>
bool matrixIsEmpty(boost::multi_array<T, 2> &matrix) {
    return matrix.shape()[0] == 0 || matrix.shape()[1] == 0;
}

void readSciData(const size_t &currScan, CalibrationData &calData, const NcGroup &l1aSciData,
                 const size_t &numCcdPix) {
    string color = determineColor(calData.color);
    vector<size_t> start = {currScan, 0, 0};
    vector<size_t> count = {1, calData.numInsBands, numCcdPix};
    l1aSciData.getVar("sci_" + color).getVar(start, count, &calData.sciData[0][0]);
}

void writeCalData(const Level1bFile &outfile, const CalibrationData &calData, const size_t &currScan,
                  const size_t numCcdPix, bool radianceGenerationEnabled,
                  const matrix2D<uint8_t> &qualityFlags) {
    string color = determineColor(calData.color);
    vector<size_t> start = {0, currScan, 0};
    vector<size_t> count = {calData.numBands, 1, numCcdPix};

    // Output to L1B file
    if (!radianceGenerationEnabled) {
        try {
            outfile.observationData.getVar("rhot_" + color)
                .putVar(start, count, &calData.calibratedData[0][0]);
        } catch (const exception &e) {
            cout  << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color << " rhot values into L1B file" << endl;
            cout << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        try {
            outfile.observationData.getVar("Lt_" + color).putVar(start, count, &calData.calibratedData[0][0]);
        } catch (const exception &e) {
            cout  << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color << " Lt values into L1B file" << endl;
            cout << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    }

    try {
        outfile.observationData.getVar("qual_" + color).putVar(start, count, calData.qualityFlags);
    } catch (const exception &e) {
        cout  << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color << " quality flags into L1B file" << endl;
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void calibrate(CalibrationData &calData, const GeoData &geoData, DarkData &darkData,
               const NcGroup &l1aSciData, const Level1bFile &outfile, const size_t currScan,
               const vector<int16_t> &spatialAgg, const size_t spatialAggregationIndex, const size_t numTaps,
               const size_t numNonlinCoefs, const vector<double> &sciScanAngles, const float *referenceTemps,
               const vec2D<double> &temps, const vector<float> &cosSolZens,
               const bool radianceGenerationEnabled) {
    bool mostInsBandsChanged = false;  // Will probably only be true once or twice

    // Specific to each CalibrationData, so numInsBands and numCcdPix won't change
    if (calData.sciData == nullptr)
        calData.sciData = (uint16_t **)allocate2d_short(calData.numInsBands, geoData.numCcdPix);
    readSciData(currScan, calData, l1aSciData, geoData.numCcdPix);
    if (calData.calibratedData == nullptr)
        calData.calibratedData = allocate2d_float(calData.numBands, geoData.numCcdPix);

    // Setup for some variables used in this function
    if (mostInsBands < calData.numInsBands) {
        if (preAggCalData != nullptr)
            free2d_float(preAggCalData);
        mostInsBands = calData.numInsBands;
        mostInsBandsChanged = true;
    }

    if (preAggCalData == nullptr)
        preAggCalData = allocate2d_float(mostInsBands, geoData.numCcdPix);
    else if (mostInsBandsChanged) {
        preAggCalData = allocate2d_float(mostInsBands, geoData.numCcdPix);  // TODO: Free first
    }

    // Compute dark offset, correct data, and apply absolute and
    // temporal gain and temperature correction
    if (darkData.corrections.size() == 0 || (mostInsBandsChanged))
        darkData.corrections.resize(mostInsBands);

    int darkCorrectionStatus = getDarkCorrection(
        currScan, geoData.numGoodScans, geoData.hamSides, darkData.numScansAvg, darkData.numPixSkip,
        spatialAgg[spatialAggregationIndex], spatialAgg[darkData.darkZone], numTaps, *calData.specAgg,
        calData.fillValue, darkData.numPix, darkData.data, darkData.corrections);

    // Raw numbers out of the CCD
    if (digitalNumbers.empty() || mostInsBandsChanged) {
        digitalNumbers.resize(mostInsBands, vector<float>());
        for (auto &vec : digitalNumbers)
            vec.resize(geoData.numCcdPix);
    }

    if (darkCorrectionStatus != -1) {  // Found valid dark data
        if (tempCorrections.size() == 0 || mostInsBandsChanged)
            tempCorrections.resize(mostInsBands);
        getTempCorrection(calData.numInsBands, referenceTemps, temps[currScan], *calData.gains,
                          tempCorrections);

        for (size_t j = 0; j < calData.numInsBands; j++) {
            for (size_t k = 0; k < geoData.numCcdPix; k++) {
                // Handle fill value
                if (calData.sciData[j][k] == calData.fillValue) {
                    digitalNumbers[j][k] = BAD_FLT;
                    preAggCalData[j][k] = BAD_FLT;
                    continue;
                }

                // Need to save dn for linearity correction
                digitalNumbers[j][k] = calData.sciData[j][k] - darkData.corrections[j];
                preAggCalData[j][k] = tempCorrections[j] *
                                      calData.gains->k1K2[j][geoData.hamSides[currScan]] *
                                      digitalNumbers[j][k];
            }
        }
    } else {
        cout << "-E- Dark correction not found.\n";
    }

    // Compute and apply RVS and linearity corrections
    if (k4.size() == 0 || mostInsBandsChanged)
        k4 = vec2D<double>(mostInsBands, vector<double>(geoData.numCcdPix));
    getRvsCorrection(calData.numInsBands, geoData.numCcdPix, geoData.hamSides[currScan], *calData.gains,
                     sciScanAngles, k4);

    if (k5.size() == 0 || mostInsBandsChanged)
        k5 = vec2D<double>(mostInsBands, vector<double>(geoData.numCcdPix));
    getNonlinearityCorrection(calData.numInsBands, geoData.numCcdPix, numNonlinCoefs, calData.gains->k5Coefs,
                              digitalNumbers, k5);

    for (size_t j = 0; j < calData.numInsBands; j++) {
        for (size_t k = 0; k < geoData.numCcdPix; k++) {
            if (preAggCalData[j][k] != BAD_FLT)
                preAggCalData[j][k] *= k4[j][k] * k5[j][k];
        }
    }

    // Aggregate to L1B bands
    // bcalb = transpose(bamat#transpose(bcal))
    aggToL1b(currScan, calData.numBands, calData.numInsBands, geoData, preAggCalData, calData.calibratedData,
             calData.insAgg, radianceGenerationEnabled, *calData.solIrrL1a, cosSolZens.data());

    // Check for saturation

    if (matrixIsEmpty(qualityFlags) || mostInsBandsChanged) {
        qualityFlags.resize(boost::extents[mostInsBands][geoData.numCcdPix]);
    }

    getQualityFlags(geoData.numCcdPix, calData.numBands, calData.numInsBands, digitalNumbers, calData.insAgg,
                    calData.gains->saturationThresholds, qualityFlags);

    calData.qualityFlags = qualityFlags.data();

    writeCalData(outfile, calData, currScan, geoData.numCcdPix, radianceGenerationEnabled, qualityFlags);
}

CalibrationLut readOciCalLut(const NcFile *calLUTfile, Device device, const NcGroup &lutGroup,
                             uint32_t &numBands, const uint32_t mcedim) {
    string color = determineColor(device);
    numBands = calLUTfile->getDim(color + "_bands").getSize();

    uint32_t numTemps;
    switch (device) {
        case RED:
        case BLUE:
            numTemps = calLUTfile->getDim("number_of_CCD_temperatures").getSize();
            break;
        case SWIR:
            numTemps = calLUTfile->getDim("number_of_SWIR_temperatures").getSize();
            break;
        default:
            throw std::invalid_argument("Unknown device type");
    }

    uint32_t numTempCoeffs = calLUTfile->getDim("number_of_T_coefficients").getSize();
    uint32_t numRvsCoefs = calLUTfile->getDim("number_of_RVS_coefficients").getSize();
    uint32_t numNonlinearityCoefs = calLUTfile->getDim("number_of_nonlinearity_coefficients").getSize();
    uint32_t numPolarizationCoefs = calLUTfile->getDim("number_of_polarization_coefficients").getSize();
    uint32_t numHamSides = 2;  // TODO: find  a better name
    uint32_t numTimes = calLUTfile->getDim("number_of_times").getSize();

    CalibrationLut calLut(numBands, numHamSides, numTimes, numTemps, numTempCoeffs, mcedim, numRvsCoefs,
                          numNonlinearityCoefs, numPolarizationCoefs);

    lutGroup.getVar("K1").getVar(&calLut.k1[0][0]);
    lutGroup.getVar("K2").getVar(&calLut.k2[0][0][0]);
    lutGroup.getVar("K3_coef").getVar(&calLut.k3Coefs[0][0][0]);
    lutGroup.getVar("K4_coef").getVar(&calLut.k4Coefs[0][0][0][0]);
    lutGroup.getVar("K5_coef").getVar(&calLut.k5Coefs[0][0]);
    lutGroup.getVar("sat_thres").getVar(&calLut.saturationThresholds[0]);

    lutGroup.getVar("m12_coef").getVar(&calLut.m12Coefs[0][0][0]);
    lutGroup.getVar("m13_coef").getVar(&calLut.m13Coefs[0][0][0]);

    return calLut;
}

Gains *makeOciGains(uint32_t numInsBands, uint32_t numBands, uint16_t year, uint32_t julianDay,
                    double scanTime, size_t numTimes, double *relGainFactors, int16_t boardId,
                    int16_t spatialAgg, int16_t *gainAgg, CalibrationLut &calLut, float **gainMat) {
    Gains *gains = (Gains *)calloc(1, sizeof(Gains));
    for (size_t i = 0; i < 6; i++)
        gains->dimensions[i] = calLut.dimensions[i];

    numTimes = gains->dimensions[T];
    uint16_t numTemps = gains->dimensions[TEMP];
    uint16_t numTempCorrections = gains->dimensions[TEMP_CORR];
    uint16_t rvsdim = gains->dimensions[RVS];
    uint16_t numNonlinCoefs = gains->dimensions[NL];
    uint16_t mirrorSideDim = gains->dimensions[MS];

    bool hyperspectral = false;
    uint16_t localBoardId;  // Used only when working with a non-hyperspectral instrument
    int16_t *insAdjFactors = NULL;

    // Hyperspectral bands
    if (boardId == -1) {
        hyperspectral = true;
        insAdjFactors = new int16_t[numInsBands];  // TODO: Find a better name
        int ib = 0;
        for (size_t i = 0; i < 16; i++) {
            if (gainAgg[i] <= 0) {
                continue;
            }

            uint32_t nb = 32 / gainAgg[i];
            for (size_t j = 0; j < nb; j++) {
                if (spatialAgg * gainAgg[i] < 4)
                    insAdjFactors[ib + j] = 4 / (spatialAgg * gainAgg[i]);
                else
                    insAdjFactors[ib + j] = 4 / 4;
            }
            ib += nb;
        }
    } else
        localBoardId = boardId % 2;

    // gains->k1K2 = allocate2d_float(numInsBands, mirrorSideDim);
    gains->k1K2 = vec2D<float>(numInsBands, vector<float>(mirrorSideDim));
    gains->k3Coefs = vec3D<float>(numInsBands, vec2D<float>(numTemps, vector<float>(numTempCorrections)));
    gains->k4Coefs = vec3D<float>(numInsBands, vec2D<float>(mirrorSideDim, vector<float>(rvsdim)));
    gains->k5Coefs = vec2D<double>(numInsBands, vector<double>(numNonlinCoefs));
    gains->saturationThresholds = vector<uint32_t>(numInsBands);

    // Mirror-side dependent gains
    double *K2 = new double[numBands];
    for (size_t mirrorSide = 0; mirrorSide < mirrorSideDim; mirrorSide++) {
        // Get temporal gain and combine with absolute gain
        double daysSinceJ2000 = julianDay - J2000 + scanTime / SECONDS_IN_DAY;

        size_t gainFactorIndex = 0;
        for (size_t j = numTimes - 1; j >= 0; j--) {
            if (daysSinceJ2000 > relGainFactors[j]) {
                gainFactorIndex = j;
                break;
            }
        }
        if (gainFactorIndex < (size_t)(numTimes - 1)) {
            double interpFactor = (daysSinceJ2000 - relGainFactors[gainFactorIndex]) /
                                  (relGainFactors[gainFactorIndex + 1] - relGainFactors[gainFactorIndex]);
            for (size_t band = 0; band < numBands; band++)
                K2[band] = calLut.k2[band][mirrorSide][gainFactorIndex] * (1.0 - interpFactor) +
                           calLut.k2[band][mirrorSide][gainFactorIndex + 1] * interpFactor;
        } else {
            for (size_t band = 0; band < numBands; band++)
                K2[band] = calLut.k2[band][mirrorSide][gainFactorIndex];
        }

        int16_t insAdjFactor = 1;  // Scales hyperspectral data, not multispectral
        for (size_t insBand = 0; insBand < numInsBands; insBand++) {
            gains->k1K2[insBand][mirrorSide] = 0;
            if (hyperspectral)
                insAdjFactor = insAdjFactors[insBand];
            for (size_t band = 0; band < numBands; band++)
                gains->k1K2[insBand][mirrorSide] +=
                    gainMat[band][insBand] * calLut.k1[band][mirrorSide] * K2[band] * insAdjFactor;
        }

        // Generate RVS coefficents
        for (size_t insBand = 0; insBand < numInsBands; insBand++) {
            for (size_t rvs = 0; rvs < rvsdim; rvs++) {
                gains->k4Coefs[insBand][mirrorSide][rvs] = 0;
                for (size_t band = 0; band < numBands; band++) {
                    if (hyperspectral)
                        gains->k4Coefs[insBand][mirrorSide][rvs] +=
                            gainMat[band][insBand] * calLut.k4Coefs[band][mirrorSide][0][rvs];
                    else
                        gains->k4Coefs[insBand][mirrorSide][rvs] +=
                            gainMat[band][insBand] * calLut.k4Coefs[band][mirrorSide][localBoardId][rvs];
                }
            }
        }
    }
    delete[] K2;

    // Generate temperature coefficients
    for (size_t i = 0; i < numInsBands; i++) {
        for (size_t k = 0; k < numTempCorrections; k++) {
            for (size_t l = 0; l < numTemps; l++) {
                gains->k3Coefs[i][l][k] = 0;
                for (size_t j = 0; j < numBands; j++)
                    gains->k3Coefs[i][l][k] += gainMat[j][i] * calLut.k3Coefs[j][l][k];
            }
        }
    }

    // Generate linearity coefficients
    for (size_t i = 0; i < numInsBands; i++) {
        for (size_t k = 0; k < numNonlinCoefs; k++) {
            gains->k5Coefs[i][k] = 0;
            for (size_t j = 0; j < numBands; j++) {
                if (hyperspectral)
                    gains->k5Coefs[i][k] +=
                        gainMat[j][i] * calLut.k5Coefs[j][k] * powf(insAdjFactors[i], (float)i);
                else
                    gains->k5Coefs[i][k] += gainMat[j][i] * calLut.k5Coefs[j][k];
            }
        }
    }

    // Generate saturation thresholds
    for (size_t i = 0; i < numInsBands; i++) {
        gains->saturationThresholds[i] = 0;
        for (size_t j = 0; j < numBands; j++) {
            if (hyperspectral)
                gains->saturationThresholds[i] +=
                    gainMat[j][i] * calLut.saturationThresholds[j] / insAdjFactors[i];
            else
                gains->saturationThresholds[i] += gainMat[j][i] * calLut.saturationThresholds[j];
        }
    }

    if (hyperspectral)
        delete[] insAdjFactors;

    return gains;
}

vec2D<double> interpTemps(NcFile *l1aFile, uint16_t numTemps, uint32_t numScans, double evtime[]) {
    using namespace boost::accumulators;
    uint32_t numTlmPackets = l1aFile->getDim("tlm_packets").getSize();
    uint32_t numDaucTemps = l1aFile->getDim("DAUC_temps").getSize();
    uint32_t numIcduThermisters = l1aFile->getDim("ICDU_therm").getSize();
    NcGroup engineeringGroup = l1aFile->getGroup("engineering_data");
    vector<double> dauctime(numTlmPackets);
    vector<double> icdutime(numTlmPackets);
    float **daucTemps = allocate2d_float(numTlmPackets, numDaucTemps);
    // Measurements from the onboard ICDU thermisters
    float **icduTherms = allocate2d_float(numTlmPackets, numIcduThermisters);

    engineeringGroup.getVar("DAUC_temp_time").getVar(dauctime.data());
    engineeringGroup.getVar("DAUC_temperatures").getVar(&daucTemps[0][0]);
    engineeringGroup.getVar("TC_tlm_time").getVar(icdutime.data());
    engineeringGroup.getVar("ICDU_thermisters").getVar(&icduTherms[0][0]);

    // Indices of required temperatures
    // MLA and lens housings (blue CCD, blue grating, red CCD, red grating)
    vector<size_t> icduIndices = {23, 24, 25, 26, 11};
    // Red and blue CCDs, SDA detectors, AOB 1, 2, 3, 4, 7, 8
    vector<size_t> daucIndices = {5,  6,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 36, 37};

    // Loop through temperatures
    vector<double> scanTimes(numScans);  // As an offset from the first earth view time
    for (size_t i = 0; i < numScans; i++)
        scanTimes[i] = evtime[i] - evtime[0];

    // ICDU thermistors
    vector<double> icduTimes(numTlmPackets);  // ICDU measurement times relative to the first earth view time
    vector<double> dummy(numTlmPackets);
    vec2D<double> calibrationTemps(numScans, vector<double>(numTemps));

    double intercept, slope;

    for (size_t i = 0; i <= 3; i++) {
        uint32_t k = 0;
        for (size_t j = 0; j < numTlmPackets; j++) {
            if (icdutime[j] > 0) {
                icduTimes[k] = icdutime[j] - evtime[0];
                dummy[k++] = icduTherms[j][icduIndices[i]];
            }
        }

        // Create accumulators
        accumulator_set<double, stats<tag::mean, tag::variance>> accX, accY;
        accumulator_set<double, stats<tag::mean>> accXY;

        // Accumulate data
        for (uint32_t j = 0; j < k; ++j) {
            accX(icduTimes[j]);
            accY(dummy[j]);
            accXY(icduTimes[j] * dummy[j]);
        }

        // Calculate regression coefficients
        double meanX = mean(accX);
        double meanY = mean(accY);
        double varx = variance(accX);
        double covXY = mean(accXY) - meanX * meanY;

        slope = covXY / varx;
        intercept = meanY - slope * meanX;

        for (size_t j = 0; j < numScans; j++) {
            calibrationTemps[j][i] = intercept + slope * scanTimes[j];
        }
    }

    size_t k = 0;
    for (size_t j = 0; j < numTlmPackets; j++)
        if (icdutime[j] > 0)
            dummy[k++] = icduTherms[j][icduIndices[4]];

    accumulator_set<double, stats<tag::mean, tag::variance>> accX, accY, accXY;
    for (size_t i = 0; i < k; ++i) {
        accX(icduTimes[i]);
        accY(dummy[i]);
        accXY(icduTimes[i] * dummy[i]);
    }

    double meanX = mean(accX);
    double meanY = mean(accY);
    double varX = variance(accX);
    double covXY = mean(accXY) - meanX * meanY;

    slope = covXY / varX;
    intercept = meanY - slope * meanX;

    for (size_t j = 0; j < numScans; j++)
        calibrationTemps[j][30] = intercept + slope * scanTimes[j];

    // DAUC temperatures
    // DAUC measurement times as an offset from the first earth view time
    vector<double> daucTimes(numTlmPackets);
    for (size_t i = 0; i < daucIndices.size(); i++) {
        size_t k = 0;
        for (size_t j = 0; j < numTlmPackets; j++)
            if (dauctime[j] > 0) {
                daucTimes[k] = dauctime[j] - evtime[0];
                dummy[k++] = daucTemps[j][daucIndices[i]];
            }

        accumulator_set<double, stats<tag::mean, tag::variance>> accX, accY, accXY;
        for (size_t j = 0; j < k; ++j) {
            accX(daucTimes[j]);
            accY(dummy[j]);
            accXY(daucTimes[j] * dummy[j]);
        }

        double meanX = mean(accX);
        double meanY = mean(accY);
        double varX = variance(accX);
        double covXY = mean(accXY) - meanX * meanY;

        slope = covXY / varX;
        intercept = meanY - slope * meanX;

        for (size_t j = 0; j < numScans; j++)
            calibrationTemps[j][i + 4] = intercept + slope * scanTimes[j];
    }

    return calibrationTemps;
}
