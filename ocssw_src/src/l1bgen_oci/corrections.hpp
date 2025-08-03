/**
 * @brief Function declarations for L1B corrections
 *
 * @authors Joel Gales (SAIC), Jakob Lindo (SSAI)
 * @date Aug 2024
 */

#ifndef __L1B_CORR__
#define __L1B_CORR__

#include <boost/multi_array.hpp>
#include "gains.hpp"
#include "dark_data.hpp"
#include "device.hpp"
#include "calibrations.hpp"

#define NUMBER_OF_TAPS 16

/**
 * @brief Generate dark corrections for OCI data
 * Averages dark collect data and corrects for bit shift and truncation as necessary. Dark data are taken
 * during the portion of the scan wherein the telescope is looking at the inside of the instrument
 * @param scanIndex Which scan (usually out of ~1710) whose dark data is being corrected
 * @param numScans Total number of scans (usually ~1710)
 * @param hamSides Array of values indicated which side of the half angle mirror was used for this scan
 * @param numScansAvg Number of dark collect scans to average (should be an odd number)
 * @param numPixSkip Number of dark pixels to skip at the start of collect
 * @param sciSpatialAgg Spatial aggregation for science data
 * @param darkSpatialAgg Spatial aggregation for dark data
 * @param numTaps Number of taps
 * @param specAgg Spectral aggregation factors for each tap. One of 0, 2, 4, 8, 16.
 * @param fillValue Fill value
 * @param numPixAverage Number of dark pixels to average at one time
 * @param darkPixels The dark collect pixel values
 * @param unused Unused parameter
 * @param darkCorrections The dark corrections for each pixel
 * @return 0 if OK, 1 if used adjacent scan(s), -1 if unable to find valid dark data
 */
int getDarkCorrection(const size_t scanIndex, const uint32_t numScans, const std::vector<uint8_t> &hamSides,
                      const uint16_t numScansAvg, const uint16_t numPixSkip, const int16_t sciSpatialAgg,
                      const int16_t darkSpatialAgg, const uint32_t numTaps,
                      const std::vector<int16_t> &specAgg, const uint32_t fillValue,
                      const int16_t numPixAverage, uint32_t ***darkPixels,
                      std::vector<double> &darkCorrections);

/**
 * @brief Calculate temperature corrections for each band
 *
 * @param numInsBands The number of instrument bands
 * @param refTemps Pointer to an array of reference temperatures for each temperature sensor
 * @param calibrationTemps Vector of measured temperatures
 * @param gains Object containing gain information
 * @param fps For which device the temperature corrections are being calculated
 * @param temperatureCorrections Vector to store the calculated temperature corrections
 */
double getTempCorrection(const size_t numInsBands, const float *refTemps,
                         const std::vector<double> &calibrationTemps, const Gains &gains, Device fpa);

/**
 * @brief Calculate response versus scan (RVS) corrections for each band and pixel
 *
 * @param numInsBands The number of instrument bands
 * @param numPixels The number of pixels
 * @param hamSide The half angle mirror side (0 or 1)
 * @param gains Object containing gain information
 * @param scanAngles Vector of scan angles
 * @param rvsCorrections Vector to store the calculated RVS corrections
 *
 * @return The response versus scsan correction for the given band and pixel
 */
double getRvsCorrection(const size_t band, const size_t pixel, const uint8_t hamSide, const Gains &gains,
                        double scanAngle);

/**
 * @brief Calculate nonlinearity corrections for one pixel in a single band
 *
 * @param numInsBands The band index for which to find a nonlinearity correction
 * @param numPixels The pixel for which to find a nonlinearity correction
 * @param numNonlinearTerms The number of nonlinear terms in the correction
 * @param k5Coefs The K5 coefficients for nonlinearity correction
 * @param digitalNumbers The digital numbers (raw sensor readings)
 * @param nonlinearityCorrections Vector to store the calculated nonlinearity corrections
 *
 * @return The nonlinearity correction for the given pixel and band
 */

double getNonlinearityCorrection(const size_t band, const size_t pixel, const uint32_t numNonlinearTerms,
                                 const vec2D<double> &k5Coefs, const vec2D<float> &digitalNumbers);

/**
 * @brief Compute and return quality flags for one line of data.
 * @param numPix The number of pixels in to work over
 * @param numBands The number of bands
 * @param numInsBands The number of bands after instrument aggregation
 * @param digitalNumbers The raw data values from the device
 * @param saturationThresholds The values above which each pixel will be considered saturated
 * @return A 2D array of flag values that indicate saturation for each pixel in the current line
 */
void getQualityFlags(size_t numPix, size_t numBands, size_t numInsBands, vec2D<float> &digitalNumbers,
                     float **insAggMat, std::vector<uint32_t> &saturationThresholds,
                     boost::multi_array<uint8_t, 2> &qualityFlags);

/**
 * @brief Read the OCI cross-correlation matrix for a band and aggregate the coefficients for the instrument
 *configuration
 * @param taps. Spatial aggregation factor. Each tap is responsible for 32 bands.
 * @param numInstrumentBands Number of instrument bands based on aggregation factors
 * @param gainAggMatrix Matrix to aggregate gains from CCD to instrument output
 * @param ncpix Number of influence pixels
 * @param cmat crosstalk coefficient matrix
 **/
void makeXtalkMat(int16_t spatialAgg, CalibrationData &calData, float ***cmat_in, size_t ncpix,
                  size_t nxbands, size_t nbands);

/**
 * @brief Compute along-scan stray light and crosstalk for OCI data from one FPA
 * @param numInstrumentBands Number of instrument bands based on aggregation factors
 * @param numCcdPix The number of pixels to work over
 * @param ncpix Number of influence pixels around each pixel
 * @param cmat crosstalk influence coefficient matrix
 * @param digitalNumbers The digital numbers (raw sensor readings)
 * @param The crosstalk correction for sensor bands in question
 **/
void getXtalkCorrection(uint32_t numInsBands, size_t numCcdPix, size_t ncpix, double ***cmat,
                        const vec2D<float> &digitalNumbers, vec2D<double> &xtalkCorrection);

#endif
