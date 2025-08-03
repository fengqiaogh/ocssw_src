#ifndef _l0_PACKET_UTIL_OCI_H_
#define _l0_PACKET_UTIL_OCI_H_

#include <stdint.h>
#include <fstream>
#include <timeutils.h>
#include "l0stream.hpp"

// L1A File Constants
#define ANCSIZE 104
#define TLMSIZE 3200

const int MAX_PACKET_SIZE = 2048;  // change from 3200 to 2048 on 03/10/22
const int MAX_LINES = 32768;
const int MAX_NUM_PACKETS = 30000;
const int NUM_TAPS = 16;

// Base spectral aggrigation factor for focal planes
const long BASE_BLUE_SPECTRAL_AGG_FACTOR = 2863311530;
const long BASE_RED_SPECTRAL_AGG_FACTOR = 2774100650;

const int16_t SWIR_LINE_OFFSET_DEFAULT[9] = {-16, 80, 64, 0, -16, 64, 0, 80, 96};
const int16_t SWIR_LINE_OFFSET_ETU[9] = {80, 88, 72, 64, 80, 72, 0, 0, 0};
const uint16_t SPATIAL_AGG_CODES[4] = {1, 2, 4, 8};

// Ancillary packet APID
const uint32_t ANCILLARY_APID = 636;

// Telemetry APIDs with time fields and spin numbers
const uint32_t APIDS_WITH_TIME_SPIN[8] = {
    711, 712, 713, 715,
    716, 717, 721, 723};  //  telemetry APIDs with time fields and spin numbers, ver 0.30, LH, 5/5/2022

// Science APIDs with spin numbers but no time fields
// 720 == SWIR packet
const uint32_t APIDS_WITH_SPIN[2] = {700, 720};

// APID valid range
const uint32_t MIN_APID = 636;
const uint32_t MAX_APID = 745;

typedef struct {
    int32_t year;
    int32_t day;
    double second;
} AncillaryPktTimeStamp;

typedef struct {
    short dataType;
    short spatialAggCode;
    short lines;
} spatialAggTable;

/**
 * @brief Grabs pixel, band and taps information from the ancillary packet
 * @param ancillaryPacket - packet to get the info from
 * @param numCcdPixels
 * @param numBlueBands
 * @param numRedBands
 * @param numSwirPixels
 * @param numDarkCcdPixels
 * @param numDarkSwirPixels
 * @param isBlueCcdTapsEnabled
 * @param isRedCcdTapsEnabled
 * @param spatialAggList
 * @return 0
 */
int getBandDimensions(uint8_t *ancillaryPacket, uint16_t &numCcdPixels, uint16_t &numBlueBands,
                      uint16_t &numRedBands, uint16_t &numSwirPixels, uint16_t &numDarkCcdPixels,
                      uint16_t &numDarkSwirPixels, uint16_t *isBlueCcdTapsEnabled,
                      uint16_t *isRedCcdTapsEnabled, spatialAggTable *spatialAggList);

/**
 * @brief Grab the packet time and save it to the passed in year, day and starttime ref
 * @param ancillaryPacket - packet to get the time from
 * @param year
 * @param day
 * @param startTtime
 * @return
 */
int getAncillaryPacketTime(uint8_t *ancillaryPacket, int32_t &year, int32_t &day, double &startTtime);

/**
 * @brief load all packets that have the same spin num as the first packet into the buffer for processing
 * @param l0Filestream - stream of L0 data
 * @param firstPacket - buffer containing data of the first packet
 * @param allPacketsBuffer - buffer for each packet related to the first packet based on spin number
 * @param numPackets - total packets for the current spin number
 * @param firstPacketSpinNum - spin number for the first packet of the series of packets
 * @param ancillaryIndex
 * @param telemetryIndices
 * @param sequenceErrorFlag
 * @param isEndFile
 * @param isSPW
 * @return
 */
int readScanPackets(L0Stream *l0Filestream, uint8_t *firstPacket,
                    uint8_t (*allPacketsBuffer)[MAX_PACKET_SIZE], uint32_t &numPackets,
                    int32_t &firstPacketSpinNum, int32_t &ancillaryIndex,
                    std::vector<int32_t> &telemetryIndices, uint8_t &sequenceErrorFlag, int32_t &isEndFile,
                    bool isSPW);

/**
 * @brief Compare 2 ancillary packets and determine if the spectral or spatial table were changed
 * @param currAncillaryPacket
 * @param nextAncillaryPacket
 * @return 0 if changed. 1 if not changed.
 */
int compareAncillaryPackets(uint8_t *currAncillaryPacket, uint8_t *nextAncillaryPacket);

/**
 * @brief Grab the spin number from the referenced packet
 * @param packet
 * @return
 */
int32_t getSpinNumFromPacket(uint32_t apid, uint8_t *packet);

/**
 * @brief Get spin number from the packet, but if the packet has a telemetry apid, save it to the telemetry
 *          spin num variable
 * @param spinNum - reference to spinNum holding spinNum for ancillary and science packet apids
 * @param spinNumTelemetry - reference to spinNum for telemetry apid
 * @param apid - apid of the packet being passed in. determines which variable to save the spinNum to
 * @param packet packet data
 */
void getSpinNumFromPacket(int32_t &spinNum, int32_t &spinNumTelemetry, uint32_t apid, uint8_t *packet);

/**
 * @brief Skip a packet from the l0 data file given a number of bytes
 * @param l0FileStream contains l0 data
 * @param packetLen number of bytes to skip
 * @param isEndFile if there is no more bytes, updates it to 1
 * @param isSPW skips the spacewire heaer if it contains it
 * @return
 */
int skipSinglePacket(L0Stream *l0FileStream, uint32_t &packetLen, int32_t &isEndFile, bool isSPW);

/**
 * @brief Reads a single packet from the file stream or just get the packet length, apid and end file flag
 * @param l0FileStream - L0 data stream where the packet will come from
 * @param packetBuffer - save packet into this pointer. NULL to just read in packetLen, apid and endFile flag
 *                          for the next file and not move the packet pointer forward.
 *
 * @param packetLen - extract packet body length from L0 and save it here
 * @param apid - id of the packet
 * @param isEndFile - end of file marker
 * @param isSPW - contains spacewire header or not
 * @return
 */
int readSinglePacket(L0Stream *l0FileStream, uint8_t *packetBuffer, uint32_t &packetLen, uint32_t &apid,
                     int32_t &isEndFile, bool isSPW);

/**
 * @brief Grab the apid and packet length of the next packet, but dont move the stream's pointer forward
 * @param l0FileStream stream to get the data from
 * @param packetLen length of next packet
 * @param apid apid of next packet
 * @param isEndFile if end file flag
 * @param isSPW skip header if there is one
 * @return
 */
int getNextPacketInfo(L0Stream *l0FileStream, uint32_t &packetLen, uint32_t &apid, int32_t &isEndFile,
                      bool isSPW);

/**
 * @brief Get the swir mode data from the packet
 * @param packets
 * @param numPackets
 * @param swirMode
 * @return
 */
int getSwirMode(uint8_t (*packets)[MAX_PACKET_SIZE], uint32_t numPackets, uint16_t &swirMode);

/**
 * @brief Takes group of 9, 20-bit samples that are packed into 23 bytes and unpacks
 * @param swirBandPacket
 * @param outDataArr
 * @return
 */
int extractSwirBandData(uint8_t *swirBandPacket, uint32_t *outDataArrRef);

/**
 * @brief Takes an ancillary packet and use the red or blue tap flags to extact ccd packet data.
 *        It also sets the spectral aggregation factors.
 * @param ancillaryPacket
 * @param isBlueCcdTapsEnabled
 * @param isRedCcdTapsEnabled
 * @param ccdId
 * @param lineNum
 * @param dataType
 * @param spatialAggregation
 * @param spectralAggregation
 * @param numBands
 * @param ccdData
 * @param overScanSumData
 * @return
 */
int unpackCcdPacket(uint8_t *ancillaryPacket, uint16_t isBlueCcdTapsEnabled[NUM_TAPS],
                    uint16_t isRedCcdTapsEnabled[NUM_TAPS], uint16_t &ccdId, uint32_t &lineNum,
                    uint16_t &dataType, uint16_t &spatialAggregation, uint16_t spectralAggregation[NUM_TAPS],
                    uint16_t &numBands, uint16_t **ccdData, uint16_t overScanSumData[NUM_TAPS]);

/** Given a packet, extract the swir line numbers, frame type and data from it
 * @brief
 * @param packet where to extract the data from
 * @param swirLines array to store line numbers
 * @param swirFrameType array to store frame type
 * @param swirData array to save the swir pixel data
 * @return
 */
int unpackSwirPacket(uint8_t *packet, int16_t *swirLineNums, uint8_t *swirFrameType, uint32_t *swirData);

/** Compare the spectral aggregation codes from when it was extracted from the packet to when it was copied
 * @brief
 * @param spectralAggCodes - from packet
 * @param bandSpectralAggCodes - copied from packet to its own array
 * @return true if matching, false otherwise
 */
bool inconsistentAggregationCodes(uint16_t *spectralAggCodes, uint16_t *bandSpectralAggCodes);

/** Compare datatypes and make sure they are the same
 * @brief
 * @param dataType
 * @param bandDataType
 * @return
 */
bool inconsistentDataTypes(uint16_t &dataType, uint16_t &bandDataType);

/**
 * @brief Unpack a science packet and extract red/blue ccd data and swir data from it
 * @param totalNumPackets
 * @param spin
 * @param numCcdPixels
 * @param numSwirPixels
 * @param maxNumOfSwirPixels
 * @param numBandsRef
 * @param isBlueCcdTapsEnabled
 * @param isRedCcdTapsEnabled
 * @param ociPacketBuffer
 * @param l0BlueBandData
 * @param l0RedBandData
 * @param l0SwirBandData
 * @param blueCcdLineNums
 * @param redCcdLineNums
 * @param swirLineNums
 * @param blueDataType
 * @param blueSpectralAggregations
 * @param redDataType
 * @param redSpectralAggregations
 * @param swirFrameTypeArrRef
 * @param returnStatus
 * @return 0 on finish
 */
int unpackScienceData(uint32_t totalNumPackets, int32_t spin, uint16_t numCcdPixels, uint16_t numSwirPixels,
                      uint16_t maxNumOfSwirPixels, uint16_t &numBandsRef,
                      uint16_t isBlueCcdTapsEnabled[NUM_TAPS], uint16_t isRedCcdTapsEnabled[NUM_TAPS],
                      uint8_t (*ociPacketBuffer)[MAX_PACKET_SIZE], uint16_t **l0BlueBandData,
                      uint16_t **l0RedBandData, uint32_t **l0SwirBandData, int16_t *blueCcdLineNums,
                      int16_t *redCcdLineNums, int16_t *swirLineNums, uint16_t &blueDataType,
                      uint16_t blueSpectralAggregations[NUM_TAPS], uint16_t &redDataType,
                      uint16_t redSpectralAggregations[NUM_TAPS], int8_t *swirFrameTypeArrRef,
                      int &returnStatus);

/** generate oci line index array for ccd and swir bands
 * @brief
 * @param spatialAggList table of 10 tables
 * @param ccdBandLineIndices
 * @param swirBandLineIndices
 * @param ccdBandDarkLineIndices
 * @param swirBandDarkLineIndices
 * @param swirLineOffsetArr
 * @return
 */
int makeOciLineIndex(spatialAggTable *spatialAggList, int16_t *ccdBandLineIndices,
                     int16_t *swirBandLineIndices, int16_t *ccdBandDarkLineIndices,
                     int16_t *swirBandDarkLineIndices, int16_t *swirLineOffsetArr);

#endif