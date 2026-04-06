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
#include "l1b_options.hpp"
#include "types.hpp"

#include <stdexcept>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

const double J2000 = 2451545;

using namespace std;
using namespace netCDF;

static vector<double> darkCorrections;
static vec2D<double> digitalNumbers;
static double *preAggCalData = nullptr;  // Instrument bands by pixels
static size_t mostInsBands;

template <typename T>
bool matrixIsEmpty(boost::multi_array<T, 2> &matrix) {
    return matrix.shape()[0] == 0 || matrix.shape()[1] == 0;
}

/**
 * @brief Reads scientific data from the L1A file for a specific scan.
 *
 * This function reads scientific data for a given scan from the L1A file and stores it in the CalibrationData
 * structure. It assumes that calData.sciData is already initialized to the proper size.
 *
 * @param currScan The current scan number.
 * @param calData The CalibrationData structure to store the read data.
 * @param l1aSciData The NetCDF group containing the scientific data.
 * @param numCcdPix The number of CCD pixels.
 */
void readSciData(const size_t &currScan, CalibrationData &calData, const NcGroup &l1aSciData,
                 const size_t &numCcdPix) {
    string color = determineColor(calData.color);

    for (size_t band = 0; band < calData.numInsBands; ++band) {
        vector<size_t> bandStart = {currScan, band, 0};
        vector<size_t> bandCount = {1, 1, numCcdPix};
        try {
            l1aSciData.getVar("sci_" + color).getVar(bandStart, bandCount, calData.sciData[band].data());
        } catch (const exceptions::NcException &e) {
            cout << "-E- Could not read science data for scan " << currScan << ". Skipping... " << endl;
            continue;
        } catch (const std::exception &e) {
            cout << "-E- Caught exception: " << e.what() << " for scan " << currScan << ". Skipping... "
                 << endl;
            continue;
        } catch (...) {
            cout << "-E- Caught unknown exception (possibly segmentation fault) for scan " << currScan
                 << ". Skipping... " << endl;
            continue;
        }
    }
}

void writeCalData(const Level1bFile &outfile, const CalibrationData &calData, const size_t &currScan,
                  const size_t numCcdPix, bool radianceGenerationEnabled) {
    string color = determineColor(calData.color);
    vector<size_t> start = {0, currScan, 0};
    vector<size_t> count = {calData.numBands, 1, numCcdPix};

    // Output to L1B file
    if (!radianceGenerationEnabled) {
        try {
            outfile.observationData.getVar("rhot_" + color).putVar(start, count, &calData.calibratedData[0]);
        } catch (const exception &e) {
            cout << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color
                 << " rhot values into L1B file" << endl;
            cout << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        try {
            outfile.observationData.getVar("Lt_" + color).putVar(start, count, &calData.calibratedData[0]);
        } catch (const exception &e) {
            cout << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color
                 << " Lt values into L1B file" << endl;
            cout << e.what() << endl;
            exit(EXIT_FAILURE);
        }
    }

    try {
        outfile.observationData.getVar("qual_" + color).putVar(start, count, calData.qualityFlags);
    } catch (const exception &e) {
        cout << "-E- " << __FILE__ << ":" << __LINE__ << " - Couldn't put " << color
             << " quality flags into L1B file" << endl;
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}

void calibrate(CalibrationData &calData, const GeoData &geoData, DarkData &darkData,
               const NcGroup &l1aSciData, const Level1bFile &outfile, const size_t currScan,
               const vector<int16_t> &spatialAgg, const size_t spatialAggregationIndex, const size_t numTaps,
               const size_t numNonlinCoefs, const vector<double> &sciScanAngles, const float *referenceTemps,
               const vec2D<double> &temps, const vector<double> &cosSolZens,
               const bool radianceGenerationEnabled) {
    bool mostInsBandsChanged = false;  // Will probably only be true once or twice
    vector<double> tempCorrection(calData.numInsBands);

    // Specific to each CalibrationData, so numInsBands and numCcdPix won't change
    if (calData.sciData.size() == 0)
        calData.sciData = vector<vector<uint16_t>>(calData.numInsBands, vector<uint16_t>(geoData.numCcdPix));
    readSciData(currScan, calData, l1aSciData, geoData.numCcdPix);
    if (calData.calibratedData == nullptr)
        calData.calibratedData = new float[calData.numInsBands * geoData.numCcdPix]();
    else
        memset(calData.calibratedData, 0, sizeof(float) * (calData.numInsBands * geoData.numCcdPix));

    // Setup for some variables used in this function
    if (mostInsBands < calData.numInsBands) {
        if (preAggCalData != nullptr)
            delete[] preAggCalData;  // Should only happen once per run
        mostInsBands = calData.numInsBands;
        mostInsBandsChanged = true;
    }

    size_t insBandsByPix = mostInsBands * geoData.numCcdPix;
    if (preAggCalData == nullptr)
        preAggCalData = new double[insBandsByPix]();  // Obviate later need to 0 on sciData fill value
    else if (mostInsBandsChanged) {
        preAggCalData = new double[insBandsByPix]();  // Obviate later need to 0 on sciData fill value
    }

    // Compute dark offset, correct data, and apply absolute and
    // temporal gain and temperature correction
    if (darkData.corrections.size() == 0 || (mostInsBandsChanged))
        darkData.corrections.resize(mostInsBands);

    int darkCorrectionStatus = getDarkCorrection(
        currScan, geoData.numGoodScans, geoData.hamSides, darkData.numScansAvg, darkData.numPixSkip,
        spatialAgg[spatialAggregationIndex], spatialAgg[darkData.darkZone], numTaps, calData.specAggModes,
        calData.fillValue, darkData.numPix, darkData.data, darkData.corrections);

    // Raw numbers out of the CCD
    if (digitalNumbers.empty() || mostInsBandsChanged) {
        digitalNumbers.resize(mostInsBands, vector<double>());
        for (auto &vec : digitalNumbers)
            vec.resize(geoData.numCcdPix);
    }

    if (darkCorrectionStatus != -1) {  // Found valid dark data
        for (size_t j = 0; j < calData.numInsBands; j++) {
            // This is k3[j]
            tempCorrection[j] = getTempCorrection(j, referenceTemps, temps[currScan], *calData.gains, CCD);
            for (size_t k = 0; k < geoData.numCcdPix; k++) {
                // Handle fill value
                if (calData.sciData[j][k] == calData.fillValue) {
                    digitalNumbers[j][k] = BAD_FLT;
                    continue;
                }

                // Need to save dn for linearity correction
                digitalNumbers[j][k] = calData.sciData[j][k] - darkData.corrections[j];
            }
        }

        vec2D<double> xtalkCorrection;
        xtalkCorrection.resize(calData.numInsBands, std::vector<double>());
        for (auto &vec : xtalkCorrection)
            vec.resize(geoData.numCcdPix, 0.0);

        if (calData.color == BLUE && calData.enableCrosstalk) {
            getXtalkCorrection(calData.numInsBands, geoData.numCcdPix, calData.ncpix, calData.cmat,
                               digitalNumbers, xtalkCorrection);
        }

        for (size_t j = 0; j < calData.numInsBands; j++) {
            for (size_t k = 0; k < geoData.numCcdPix; k++) {
                size_t dataIndex = j * geoData.numCcdPix + k;

                if (digitalNumbers[j][k] == BAD_FLT) { // science data was fill
                    preAggCalData[dataIndex] = BAD_FLT;
                    continue; // Bad data, no corrections
                }

                if(preAggCalData[dataIndex] == BAD_FLT)
                    continue;

                preAggCalData[dataIndex] = tempCorrection[j] *
                                           calData.gains->k1K2[j][geoData.hamSides[currScan]] *
                                           (digitalNumbers[j][k] - xtalkCorrection[j][k]);
                // This is k4
                double rvsCorrection =
                    getRvsCorrection(j, k, geoData.hamSides[currScan], *calData.gains, sciScanAngles[k]);
                // This is k5
                double nonlinearityCorrection =
                    getNonlinearityCorrection(j, k, numNonlinCoefs, calData.gains->k5Coefs, digitalNumbers);
                preAggCalData[dataIndex] *= rvsCorrection * nonlinearityCorrection;
            }
        }

        // Aggregate to L1B bands
        // bcalb = transpose(bamat#transpose(bcal))
        aggAndCalcRefls(currScan, calData.numBands, calData.numInsBands, geoData, preAggCalData,
                        calData.calibratedData, calData.insAgg, radianceGenerationEnabled, calData.solIrrL1a,
                        cosSolZens.data());

        // Check for saturation
        getQualityFlags(calData, geoData.numCcdPix, digitalNumbers);

    } else {
        cout << "-W- Dark correction not found, skipping calibration corrections\n";
        for (size_t i = 0; i < geoData.numCcdPix; i++) {
            for (size_t j = 0; j < calData.numBands; j++) {
                size_t dataIndex = i * geoData.numCcdPix + j;
                calData.calibratedData[dataIndex] = BAD_FLT;
                calData.qualityFlags[j * geoData.numCcdPix + i] = 0;
            }
        }
    }

    writeCalData(outfile, calData, currScan, geoData.numCcdPix, radianceGenerationEnabled);
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
    uint32_t numHamSides = 2;
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
    Gains *gains = new Gains();
    for (size_t i = 0; i < GAIN_DIMS_SIZE; i++)
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
        insAdjFactors = new int16_t[numInsBands];
        int bandIndex = 0;
        for (size_t i = 0; i < 16; i++) {
            if (gainAgg[i] <= 0) {  // Prevent a divide by zero error
                continue;
            }

            uint32_t numBands = 32 / gainAgg[i];
            for (size_t j = 0; j < numBands; j++) {
                if (spatialAgg * gainAgg[i] < 4)
                    insAdjFactors[bandIndex + j] = 4 / (spatialAgg * gainAgg[i]);
                else
                    insAdjFactors[bandIndex + j] = 4 / 4;
            }
            bandIndex += numBands;
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

CalibrationData::CalibrationData(const netCDF::NcFile *l1aFile, const GeoData &geoData, Device device,
                                 uint32_t gainBands, size_t k2tTimes, double *k2t,
                                 vector<int16_t> &spatialAgg, size_t spatialAggIndex, CalibrationLut &calLut,
                                 vector<double> &solarIrradiances, const bool aggregationOff) {
    color = device;
    sciData = std::vector<std::vector<uint16_t>>();
    cmat = nullptr; // Only used for blue
    enableCrosstalk = false; // Only used for blue
    ncpix = 0;
    gains = nullptr;
    darkData = nullptr;
    calibratedData = nullptr;

    std::string colorStr = determineColor(color);
    std::string colorBands = colorStr + "_bands";

    {
        bool _;
        l1aFile->getGroup("science_data").getVar("sci_" + colorStr).getFillModeParameters(_, fillValue);
    }

    numBandsL1a = l1aFile->getDim(colorBands).getSize();
    numInsBands = numBandsL1a;
    qualityFlags = (uint8_t *)malloc(numBandsL1a * geoData.numCcdPix);

    insAgg = new float *[numInsBands];
    if (color == SWIR) {
        numBands = NUM_SWIR_WAVELENGTHS;
        gainAgg = allocate2d_float(numBands, numBands);

        for (size_t i = 0; i < numBands; i++) {
            for (size_t j = 0; j < numBands; j++) {
                gainAgg[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

    } else {
        gainAgg = new float *[NUM_CCD_WAVELENGTHS];
    }

    size_t numTaps = l1aFile->getDim("number_of_taps").getSize();

    if (color != SWIR) {  // This var doesn't exist for SWIR
        specAggModes = std::vector<int16_t>(numTaps);
        std::string colorSpectralMode = colorStr + "_spectral_mode";

        l1aFile->getGroup("spatial_spectral_modes").getVar(colorSpectralMode).getVar(specAggModes.data());
    }

    vector<size_t> tapAggFactors(numTaps);  // A list of scalars, 0-16 inclusive
    size_t binCounts[NUMBER_OF_TAPS]{};

    // *********************************************************** //
    // *** Generate matrices for spectral and gain aggregation *** //
    // *********************************************************** //

    int status = 0;

    if (color != SWIR && !aggregationOff) {
        status =
            aggregateBands(numTaps, tapAggFactors.data(), binCounts, specAggModes.data(), numInsBands, numBands);
    } else { // Either SWIR cal data or aggregation has been turned off
        numBands = numBandsL1a;
    }

    if (status == 2) {
        std::cout << "-W- All " << colorStr << " bands disabled" << std::endl;
    }

    if (numInsBands != 1 && color != SWIR) {
        getAggregationMatrices(tapAggFactors.data(), specAggModes.data(), binCounts, numInsBands, numBands,
                               aggregationOff, insAgg, gainAgg, color);
    }

    if (numInsBands != numBandsL1a) {
        std::string error =
            "Number of " + colorStr + " in input file not consistent with spectral aggregation";
        throw std::runtime_error(error);
    } else if (numInsBands < 4) {
        std::cout << "-W- No " + colorStr + " bands in input file" << std::endl;
    }

    // Note: blue & red gain matrix are not necessarily contiguous

    int16_t year, month, dayOfMonth;
    double _;
    unix2ymds(geoData.unixTimeStart, &year, &month, &dayOfMonth, &_);
    int32_t astronomicalJulianDay = jday(year, month, dayOfMonth);

    if (numInsBands >= 4 && color != SWIR) {
        gains = makeOciGains(numInsBands, gainBands, year, astronomicalJulianDay, geoData.earthViewTimes[0],
                             k2tTimes, k2t, -1, spatialAgg[spatialAggIndex], specAggModes.data(), calLut, gainAgg);
    } else if (color == SWIR) {
        gains = makeOciGains(numBands, numBands, year, astronomicalJulianDay, geoData.earthViewTimes[0],
                             k2tTimes, k2t, geoData.mceTelem.mceBoardId, -1, NULL, calLut, gainAgg);
    }

    // Just in case num blue != num red.
    const size_t LOCAL_NUM_WAVELENGTHS = color == BLUE ? NUM_BLUE_WAVELENGTHS : NUM_RED_WAVELENGTHS;

    if (color != SWIR) {  // SWIR doesn't get aggregated
        vector<double> wavelengthsByInsBands =
            aggregateIrradiances(numInsBands, LOCAL_NUM_WAVELENGTHS, solarIrradiances, gainAgg);
        solIrrL1a = aggregateIrradiances(numBands, numInsBands, wavelengthsByInsBands, insAgg);
    } else { // But we still need the irradiances
        solIrrL1a = solarIrradiances;
    }
}

void CalibrationData::writeSaturationPercent(Level1bFile &outfile, const GeoData &geoData) {

    if (!outfile.l1bFile) {
        std::cout << "-E- Outfile doesn't exist";
        exit(EXIT_FAILURE);
    }

    std::vector<float> percentages;
    size_t numPix = geoData.numCcdPix;
    size_t numScans = geoData.numGoodScans;

    if (geoData.numCcdPix != geoData.numSwirPix && color == SWIR) {
        numPix = geoData.numSwirPix;
    }
    
    for (size_t count : saturatedPixelCounts) {
        // Num pix * num scans because this is the percentage of all pixels in this band
        float proportion = static_cast<float>(count) / static_cast<float>(numPix * numScans);
        percentages.push_back(proportion * 100);
    }

    std::string ncVarName = determineColor(color) + "_percent_saturated";

    try {
        NcVar percentSaturatedVar = outfile.l1bFile->getGroup("observation_data").getVar(ncVarName);
        percentSaturatedVar.putVar(percentages.data());
    } catch (const NcException &e) {
        std::cout << "-E- Couldn't put saturation percentage for " << determineColor(color) << ": "
                  << e.what();
        exit(EXIT_FAILURE);
    }
}

void CalibrationData::writeBandInfo(Level1bFile& outfile, const CalibrationLut& calLut,
                                    const vector<float>& l1aWavelengths) {
    vector<float> l1bWavelengths = vector<float>(numBands);
    vector<float> sum = vector<float>(numInsBands);
    string color = determineColor(this->color);

    // Right now, blue wavelengths equals red wavelengths
    const size_t LOCAL_NUM_WAVELENGTHS = this->color == SWIR ? NUM_SWIR_WAVELENGTHS : NUM_BLUE_WAVELENGTHS;

    if (this->color != SWIR) {
        for (size_t i = 0; i < numInsBands; i++) {
            sum[i] = 0.0;
            for (size_t j = 0; j < LOCAL_NUM_WAVELENGTHS; j++) {
                sum[i] += l1aWavelengths[j] * gainAgg[j][i];
            }
        }

        // bamat#sum
        for (size_t i = 0; i < numBands; i++) {
            l1bWavelengths[i] = 0.0;
            for (size_t j = 0; j < numInsBands; j++) {
                l1bWavelengths[i] += insAgg[j][i] * sum[j];
            }
        }
    }

    const size_t HAM_SIDES = 2;
    const size_t POLARIZATION_COEFS = 3;
    float*** l1bM12 = allocate3d_float(numBands, HAM_SIDES, POLARIZATION_COEFS);
    float*** l1bM13 = allocate3d_float(numBands, HAM_SIDES, POLARIZATION_COEFS);

    auto insAggTimesMeullerTimesGainAgg = [&](float*** l1b, float*** out) {
        for (size_t coef = 0; coef < POLARIZATION_COEFS; coef++) {
            for (size_t hamSide = 0; hamSide < HAM_SIDES; hamSide++) {
                for (size_t band = 0; band < numBands; band++) {
                    float& result = l1b[band][hamSide][coef];
                    result = 0.0;
                    for (size_t insBand = 0; insBand < numInsBands; insBand++) {
                        float sum = 0.0;
                        for (size_t wave = 0; wave < LOCAL_NUM_WAVELENGTHS; wave++) {
                            sum += out[wave][hamSide][coef] * gainAgg[wave][insBand];
                        }
                        result += insAgg[insBand][band] * sum;
                    }
                }
            }
        }
    };

    if (this->color != SWIR) {
        insAggTimesMeullerTimesGainAgg(l1bM12, calLut.m12Coefs);
        insAggTimesMeullerTimesGainAgg(l1bM13, calLut.m13Coefs);
    }

    NcGroup sensorBandParameters = outfile.l1bFile->getGroup("sensor_band_parameters");

    sensorBandParameters.getVar(color + "_solar_irradiance").putVar({0}, {numBands}, solIrrL1a.data());

    vector<size_t> start(3, 0);
    vector<size_t> count(3, 0);
    count[0] = numBands;
    count[1] = 2;
    count[2] = 3;

    if (this->color == SWIR) {
        sensorBandParameters.getVar(color + "_wavelength").putVar({0}, {numBands}, l1aWavelengths.data());
        sensorBandParameters.getVar("SWIR_bandpass").putVar({0}, {NUM_SWIR_WAVELENGTHS}, SWIR_BANDPASS);
        sensorBandParameters.getVar(color + "_m12_coef").putVar(start, count, &calLut.m12Coefs[0][0][0]);
        sensorBandParameters.getVar(color + "_m13_coef").putVar(start, count, &calLut.m13Coefs[0][0][0]);
    } else {
        sensorBandParameters.getVar(color + "_wavelength").putVar({0}, {numBands}, l1bWavelengths.data());
        sensorBandParameters.getVar(color + "_m12_coef").putVar(start, count, &l1bM12[0][0][0]);
        sensorBandParameters.getVar(color + "_m13_coef").putVar(start, count, &l1bM13[0][0][0]);

        // Doesn't exist on some inputs
        NcVar spectralMode = sensorBandParameters.getVar(color + "_spectral_mode");
        if (!spectralMode.isNull()) {
            spectralMode.putVar({0}, {outfile.l1bFile->getDim("number_of_taps").getSize()}, specAggModes.data());
        }
    }

    free3d_float(l1bM12);
    free3d_float(l1bM13);

}
