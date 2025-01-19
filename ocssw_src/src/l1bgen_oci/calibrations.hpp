/**
 * @file l1b_calibrations.h
 * @brief Header file for OCI L1B calibration functions
 *
 * This file contains function declarations for reading calibration lookup tables,
 * calculating OCI gains, and retrieving calibration temperatures for the OCI instrument.
 *
 * @authors Jakob C Lindo (SSAI)
 *
 * @date Aug 2024
 */

#ifndef __L1B_CALIBRATIONS_H__
#define __L1B_CALIBRATIONS_H__

#include <netcdf>
#include "gains.hpp"
#include <boost/multi_array.hpp>
#include "allocate2d.h"
#include "allocate3d.h"
#include "allocate4d.h"
#include "dark_data.hpp"
#include "l1b_file.hpp"
#include "device.hpp"
#include "geo_data.hpp"

// NOTE: we are comparing to the int val of solz (scale=0.01)
#define MAX_SOLZ (88 * 100)

struct CalibrationData {  // A collection of data specific to on CCD needed for calibration.

    // Inputs
    Device color;  // Which device (of RED, BLUE, SWIR) owns this data
    size_t numInsBands = 1;
    size_t numBands = 1;
    size_t numBandsL1a;
    uint32_t fillValue;
    // TODO: std container
    uint16_t **sciData;  // Science data from L1A. Bands by pixels
    float **insAgg;      // Instrument aggregation matrix. numInsBands by numBands
    float **gainAgg;     // Gain aggregation matrix
    std::vector<int16_t> *specAgg;
    std::vector<double> *solIrrL1a;  // Aggregated solar irradiances for each L1A band
    Gains *gains;
    DarkData *darkData;

    // Outputs
    float **calibratedData = nullptr;  // Bands by pixels
    uint8_t *qualityFlags = nullptr;   // Bands by pixels

    ~CalibrationData() {
        if (gainAgg)
            delete[] gainAgg;
        if (insAgg)
            delete[] insAgg;
        if (sciData)
            free2d_short((short**) sciData);
        if (calibratedData)
            free2d_float(calibratedData);
    }
};

struct CalibrationLut {
    uint16_t dimensions[7];  // The shape of this LUT
    float **k1;              // Absolute gain factor
    float ***k2;             // Relative gain factor over time
    float ***k3Coefs;        // Temperature correction factor
    float ****k4Coefs;       // Response vs scan (one dim is scan angle, usually replaced by pixel number)
    double **k5Coefs;        // Nonlinearity factor
    uint32_t *saturationThresholds;
    float ***m12Coefs;  // Rotation of polarization, described by the Mueller Matrix at (1, 2)
    float ***m13Coefs;  // Describes how much OCI prefers to attenuate light along x compared to z

    CalibrationLut(const size_t numBands, const size_t numHamSides, const size_t numTimes,
                   const size_t numTemps, const size_t numTempCoefs, const size_t mceDim,
                   const size_t numRvsCoefs, const size_t numNonlinCoefs, const size_t numPolarizationCoefs) {
        k1 = allocate2d_float(numBands, numHamSides);
        k2 = allocate3d_float(numBands, numHamSides, numTimes);  // A function of time
        k3Coefs = allocate3d_float(numBands, numTemps, numTempCoefs);
        k4Coefs = allocate4d_float(numBands, numHamSides, mceDim, numRvsCoefs);  // An angle
        k5Coefs = allocate2d_double(numBands, numNonlinCoefs);
        saturationThresholds = new uint32_t[numBands];
        m12Coefs =
            allocate3d_float(numBands, numHamSides, numPolarizationCoefs);  // Polarization sensitivity
        m13Coefs =
            allocate3d_float(numBands, numHamSides, numPolarizationCoefs);  // Polarization sensitivity

        // Assign cal lut dimensions array
        dimensions[0] = numTimes;
        dimensions[1] = numTemps;
        dimensions[2] = numTempCoefs;
        dimensions[3] = numRvsCoefs;
        dimensions[4] = numNonlinCoefs;
        dimensions[5] = numHamSides;
        dimensions[6] = numPolarizationCoefs;
    }

    ~CalibrationLut() {
        free2d_float(k1);
        free3d_float(k2);
        free3d_float(k3Coefs);
        free4d_float(k4Coefs);
        free2d_double(k5Coefs);
        delete[] saturationThresholds;
        free3d_float(m12Coefs);
        free3d_float(m13Coefs);
    }
};

/**
 * @brief Read a calibration look up table for OCI
 * @param calLUTfile The OCI calibration look up table
 * @param device Indicates which device (between red and blue) is being read for
 * @param gidLUT The NcGroup that contains the LUT for a specific device
 * @param numGainBands The number of bands per pixel
 * @param mcedim A dimension for response vs scan
 * @return The LUT
 */
CalibrationLut readOciCalLut(const netCDF::NcFile *calLutFile, const Device device,
                             const netCDF::NcGroup &calLutGroup, uint32_t &numGainBands,
                             const uint32_t mcedim);
/**
 * @brief Make a Gains object
 * @param numInsBands The number of instrument bands
 * @param banddim The number of gain bands
 * @param year The year in which this file began
 * @param julianDay The astronomical Julian day in which this file began
 * @param scanTime The time of day of the start of the first scan in this file
 * @param numTimes The number of times recorded in this file
 * @param relGainFactors Relative gain factors
 * @param boardId The MCE board ID. -1 if an OCI CCD
 * @param spatialAgg A spatial aggregation factor
 * @param specAgg The spectral aggregation matrix
 * @param calLut A calibration look up table for this a specific device
 * @param gainMat Gain matrix
 * @return A pointer to a Gains object
 */
Gains *makeOciGains(uint32_t numInsBands, uint32_t banddim, uint16_t year, uint32_t julianDay,
                    double scanTime, size_t numTimes, double *relGainFactors, int16_t boardId,
                    int16_t spatialAgg, int16_t *specAgg, CalibrationLut &calLut, float **gainMat);

/**
 * @brief Read and return temperatures from L1A, interpolated against scan times.
 * @param l1afile The L1A file
 * @param numTemps Number of temperatures
 * @param numScans Number of scans whose temps will be read
 * @param earthViewTimes Earth view times
 * @return The temperatures from L1A, interpolated to scan times
 */
vec2D<double> interpTemps(netCDF::NcFile *l1afile, uint16_t numTemps, uint32_t numScans,
                          double earthViewTimes[]);

/**
 * @brief Read L1A science data, perform the L1B calibration/correction algorithm, and write the results to an
 * L1B file. This function calibrates one line at a time.
 * @param calData Information specific to one CCD that is required for calibration
 * @param geoData The output of geolocation which is needed for some corrections
 * @param darkData Information regarding data captured when the RTA was facing the inside of the housing
 * @param l1aSciData The raw data out of the CCD
 * @param outfile The L1B file to write to
 * @param currScan The scan number being processed
 * @param spatialAggregationIndex The index into calData.spatAgg
 * @param numTaps The number of taps in this CCD. Usually 16
 * @param numNonlinCoefs The number of nonlinearity coefficients, used in nonlinearity correction
 * @param referenceTemps The temperature measurements (citation needed)
 * @param sciScanAngles The angles of each scan
 * @param temps The temperatures to compare to during temperature correction
 * @param cosSolZens The cosines of each solar zenith
 * @param radianceGenerationEnabled Whether we are processing radiance
 */
void calibrate(CalibrationData &calData, const GeoData &geoData, DarkData &darkData,
               const netCDF::NcGroup &l1aSciData, const Level1bFile &outfile, const size_t currScan,
               const std::vector<int16_t> &spatialAgg, const size_t spatialAggregationIndex,
               const size_t numTaps, const size_t numNonlinCoefs, const std::vector<double> &sciScanAngles,
               const float *referenceTemps, const vec2D<double> &temps, const std::vector<float> &cosSolZens,
               const bool radianceGenerationEnabled);

#endif
