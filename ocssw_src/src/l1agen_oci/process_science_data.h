#ifndef _PROCESS_SCIENCE_DATA_HPP_
#define _PROCESS_SCIENCE_DATA_HPP_

#include <stdint.h>
#include <string>

/** copy and filter out bad band data extracted from L0 file and write the good ones to the science arr
 * @brief
 * @param sciArr to write into
 * @param l0BandData from L0 file
 * @param numCcdPixels
 * @param numBands
 * @param ccdLineNums color specific, red or blue
 * @param ccdBandLineIndices color specific, red or blue
 */
void copyGoodL0BandDataToSciArr(uint16_t **sciArr, uint16_t **l0BandData, uint16_t numCcdPixels,
                                uint16_t numBands, int16_t *ccdLineNums, int16_t *ccdBandLineIndices);

/**
 * @brief copy and filter out bad band data from l0 and write the good ones to dark science data
 * @param darkSciArr
 * @param l0BandData
 * @param numCcdPixels
 * @param numDarkCcdPixels
 * @param numBands
 * @param ccdLineNums
 * @param ccdBandLineIndices
 */
void copyGoodSwirL0BandDataToDarkSciArr(uint32_t **darkSciArr, uint32_t **l0BandData, uint16_t numPixels,
                                        uint16_t numDarkPixels, uint16_t numBands, int16_t *ccdLineNums,
                                        int16_t *ccdBandLineIndices);

/**
 * @brief copy and filter out bad band data from l0 and write the good ones to dark science data
 * @param darkSciArr
 * @param l0BandData
 * @param numCcdPixels
 * @param numDarkCcdPixels
 * @param numBands
 * @param ccdLineNums
 * @param ccdBandLineIndices
 */
void copyGoodL0BandDataToDarkSciArr(uint16_t **darkSciArr, uint16_t **l0BandData, uint16_t numCcdPixels,
                                    uint16_t numDarkCcdPixels, uint16_t numBands, int16_t *ccdLineNums,
                                    int16_t *ccdBandLineIndices);

/**
 * @brief Return the number of ccd band line indices that are -1
 * @param numPixels ccd pixels
 * @param ccdLineNums red or blue ccd line numbers
 * @param ccdBandLineIndicies
 * @return num of ccd band line indicies that are -1
 */
int getNumOfInvalidCcdLineIndicies(uint16_t numPixels, int16_t *ccdLineNums, int16_t *ccdBandLineIndicies);

/**
 * @brief check to see if ccd band line indicies are out of order
 * @param errorToReport message to print when sequence is out of order
 * @param numPixels ccd pixels
 * @param ccdLineNums red or blue line numbers array
 * @param ccdBandLineIndices
 * @return 1 if out of order. 0 if not.
 */
int isSequenceOutOfOrder(std::string errorToReport, uint16_t numPixels, int16_t *ccdLineNums,
                         int16_t *ccdBandLineIndices);

int checkMissingCcdPixels(std::string errorToReport, bool checkCcdLines, uint16_t numPixels,
                          uint16_t numDarkPixels, int16_t *ccdLineNums, int16_t *ccdBandLineIndices);

/**
 * @brief extract good l0 band data for swir, red and blue and put them in a clean array for writing
 * @param dataType
 * @param numCcdPixels
 * @param numSwirPixels
 * @param numDarkCcdPixels
 * @param numDarkSwirPixels
 * @param numBlueBands
 * @param numRedBands
 * @param numSwirBands
 * @param ccdBandLineIndices
 * @param swirBandLineIndices
 * @param ccdBandDarkLineIndices
 * @param swirBandDarkLineIndices
 * @param l0BlueBandData
 * @param l0RedBandData
 * @param l0SwirBandData
 * @param blueCcdLineNums
 * @param redCcdLineNums
 * @param swirLineNums
 * @param blueScienceData
 * @param redScienceData
 * @param swirScienceData
 * @param blueCcdDarkData
 * @param redCcdDarkData
 * @param swirDarkData
 * @param ccdLineError
 * @param sciDataStatus
 * @return
 */
int checkAndLoadScienceData(short dataType, uint16_t numCcdPixels, uint16_t numSwirPixels,
                            uint16_t numDarkCcdPixels, uint16_t numDarkSwirPixels, uint16_t numBlueBands,
                            uint16_t numRedBands, uint16_t numSwirBands, int16_t *ccdBandLineIndices,
                            int16_t *swirBandLineIndices, int16_t *ccdBandDarkLineIndices,
                            int16_t *swirBandDarkLineIndices, uint16_t **l0BlueBandData,
                            uint16_t **l0RedBandData, uint32_t **l0SwirBandData, int16_t *blueCcdLineNums,
                            int16_t *redCcdLineNums, int16_t *swirLineNums, uint16_t **blueScienceData,
                            uint16_t **redScienceData, uint32_t **swirScienceData, uint16_t **blueCcdDarkData,
                            uint16_t **redCcdDarkData, uint32_t **swirDarkData, int8_t &ccdLineError,
                            int &sciDataStatus);

#endif