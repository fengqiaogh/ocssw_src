#ifndef __L1B_AGG__
#define __L1B_AGG__

#include <stdint.h>
#include <stddef.h>
#include "geolocate_oci.h"

#define NUMBER_OF_TAPS 16
#define BANDS_PER_TAP 32
#define NUMBER_OF_BANDS 512

/**
 * @brief Aggregate bands
 * @param taps The taps. Each tap is responsible for 32 bands.
 * @param tapAggFactors A set of scalar values indicating the aggregation factor for each tap. Each values is
 * expected to be one of 0, 1, 2, 4, or 8, with each number indicating the size of each aggregation in the tap
 * @param binCounts A list whose indices correspond to the taps and whose values correspond to the number of
 * aggregations for that tap
 * @param aggregationFlags Aggregation factor enable flag for each tap
 * @param numInstrumentBands Number of instrument bands based on aggregation factors
 * @param numL1bBands Number of L1B bands after aggregation
 * @return 2 if all taps for this band are disabled, 0 otherwise
 */
int aggregateBands(const size_t taps, size_t *tapAggFactors, size_t binCounts[NUMBER_OF_TAPS],
                   const int16_t aggregationFlags[NUMBER_OF_TAPS], size_t &numInstrumentBands,
                   size_t &numL1bBands);

/**
 * @brief Generate matrices to to help with aggregation of gain corrections
 * @param taps The taps. Each tap is responsible for 32 bands.
 * @param tapAggFactors A set of scalar values indicating the aggregation factor for each tap. Each values is
 * expected to be one of 0, 1, 2, 4, or 8, with each number indicating the size of each aggregation in the tap
 * @param binCounts A list whose indices correspond to the taps and whose values correspond to the number of
 * aggregations for that tap
 * @param numInstrumentBands Number of instrument bands based on aggregation factors
 * @param numL1bBands Number of L1B bands after aggregation
 * @param instrumentAggMatrix Matrix to aggregate instrument bands to L1B bands
 * @param gainAggMatrix Matrix to aggregate gains from CCD to instrument output
 * @return
 */
int getAggregationMatrices(size_t *taps, int16_t tapAggFactors[NUMBER_OF_TAPS],
                           size_t binCounts[NUMBER_OF_TAPS], uint32_t numInstrumentBands,
                           uint32_t numL1bBands, float **instrumentAggMatrix, float **gainAggMatrix);

/**
 * @brief Populate the instrument aggregation matrix with the correct values
 * @param instrumentAggMatrix  The instrument aggregation matrix. Will be mutated on each call
 * @param tapAggFactors A set of scalar values indicating the aggregation factor for each tap. Each values is
 * expected to be one of 0, 1, 2, 4, or 8, with each number indicating the size of each aggregation in the tap
 * @param binCounts A list whose indices correspond to the taps and whose values correspond to the number of
 * aggregations for that tap
 * @param tapBounds Each tap's beginning and end boundaries in terms of band indices
 * @param taps The taps. Each tap is responsible for 32 bands.
 */
void populateInstrumentAggMatrix(float **instrumentAggMatrix, int16_t *tapAggFactors,
                                 size_t binCounts[NUMBER_OF_TAPS], int16_t tapBounds[NUMBER_OF_TAPS][2],
                                 size_t *taps);

/**
 * @brief Populate the gain aggregation matrix with the correct values
 * @param gainAggMatrix The gain aggregation matrix. Will be mutated on each call
 * @param tapAggFactors A set of scalar values indicating the aggregation factor for each tap. Each values is
 * expected to be one of 0, 1, 2, 4, or 8, with each number indicating the size of each aggregation in the tap
 * @param binCounts A list whose indices correspond to the taps and whose values correspond to the number of
 * aggregations for that tap
 * @param numInstrumentBands Number of instrument bands based on aggregation factors
 * @param tapBounds Each tap's beginning and end boundaries in terms of band indices
 */
void populateGainAggMatrix(float **gainAggMatrix, int16_t tapAggFactors[NUMBER_OF_TAPS],
                           size_t binCounts[NUMBER_OF_TAPS], uint32_t numInstrumentBands,
                           int16_t tapBounds[NUMBER_OF_TAPS][2]);

/**
 * @brief Aggregate calibration data to L1B bands
 * @param currScan Scan index into this current file
 * @param numBands Number of bands
 * @param numInsBands Number of bands to aggregate to
 * @param geoData Output of geolocation
 * @param preAgg Pre-aggregated data. Instrument bands X pixels
 * @param aggregated Agreggated & calibrated data. L1B bands X pixels
 * @param insAggMat Instrument aggregation matrix
 * @param radianceGen Indicates whether radiance generation is enabled
 * @param solIrrL1a Solar irradiances from L1A
 * @param cosSolZens Cosine of solar zeniths
 */
void aggAndCalcRefls(size_t currScan, size_t numBands, size_t numInsBands, const GeoData &geoData,
                     float *preAgg, float *aggregated, float **insAggMat, bool radianceGen,
                     std::vector<double> &solIrrL1a, const float *cosSolZens);

#endif