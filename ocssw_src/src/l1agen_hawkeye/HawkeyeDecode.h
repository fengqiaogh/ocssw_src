/*
		HawkeyeDecode.h

		Matt Longmire
		Cloudland Instruments

		Contains the interface to the utility functions for
		decoding a Hawkeye Format Image file (stream).


*/

#ifndef _HAWKEYE_DECODE_
#define _HAWKEYE_DECODE_

#include "Hawkeye.h"

#define MAX_FINDERSCOPE_IMAGES	25
#define MAX_TELEMETRY_RECORDS	12
#define INCLUDE_GPS				0		// set to 1 to decode binaries with GPS data

typedef enum { HSE_NO_ERROR, HSE_NO_IMAGE_PARAMETERS, HSE_NO_HEADER_FOUND } HAWKEYE_STREAM_ERROR;

typedef enum { EBT_NULL, EBT_UNCOMPRESSED, EBT_COMPRESSED, EBT_IMAGE_PARAMS, EBT_TELEMETRY, EBT_EOF, EBT_MISSION_LOG } ENCODED_BLOCK_TYPE;

typedef struct {
  uint32_t recordNo;
  ENCODED_BLOCK_TYPE blockType;
  uint16_t blockLen;
} HeaderInfo;

typedef struct {
  uint32_t timeStamp;
  int height, width;
  int noMissingRows;
  int averageDarkRowPresent;
  double electronicGain;
} HawkeyeBandInfo;

typedef struct {
  uint32_t timeStamp;
  uint16_t noChannels;
  GetTelemetryResponse telemetry;
} HawkeyeTelemetryInfo;

typedef struct {
  uint16_t errorCode;
  uint32_t exposureID;
  uint32_t imageID;
  uint64_t epochT0;
  uint32_t hostDeltaEpoch;
  uint32_t hawkeyeDeltaEpoch;
  uint16_t spectralBinning;
  uint16_t finderscopeBinning;
  uint16_t channelBitfield;
  uint16_t ccd1Exposure;
  uint16_t ccd2Exposure;
  uint16_t ccd3Exposure;
  uint16_t ccd4Exposure;
  uint16_t height;
  uint16_t darkHeight;
  uint16_t interval;
  uint16_t oversampling;
  uint16_t finderscopeExposure;
  uint16_t noFinderscopeImages;
  DATA_COMPRESSION compression;
  uint16_t darkSubtracted;
  uint16_t shutterSolenoid;
  uint16_t readoutOrder;
#if INCLUDE_GPS
  uint8_t gpsBinary[116];
#endif
} HawkeyeImageInfo;

typedef struct {
  uint32_t streamLength;
  uint32_t noRecords;
  uint32_t noUnknownRecords;
  uint32_t noMissingRecords;
  uint32_t noBadRecords;
  uint32_t noMissingBytes;
  HawkeyeImageInfo imageInfo;
  int noSpectralImages;
  HawkeyeBandInfo spectralInfo[8];
  int noFinderscopeImages;
  HawkeyeBandInfo finderscopeInfo[MAX_FINDERSCOPE_IMAGES];
  int noTelemetryRecords;
  HawkeyeTelemetryInfo telemetryInfo[MAX_TELEMETRY_RECORDS];
  char missionLog[MISSION_LOG_LENGTH];
  uint32_t hebBufSize;
 } HawkeyeStreamInfo;

#ifdef __cplusplus

extern "C" HAWKEYE_STREAM_ERROR HawkeyeScanStream(uint8_t *stream, uint32_t streamLength, HawkeyeStreamInfo *streamInfo);
extern "C" HAWKEYE_STREAM_ERROR HawkeyeDecodeSpectralImage(uint16_t bandNo, uint16_t *pixels, uint32_t pixelLength,
                                                           uint16_t *averageDarkPixels, uint16_t averageDarkPixelLength);
extern "C" HAWKEYE_STREAM_ERROR HawkeyeDecodeFinderscopeImage(uint16_t imageNo, uint16_t *pixels, uint32_t pixelLength);

#else

extern HAWKEYE_STREAM_ERROR HawkeyeScanStream(uint8_t *stream, uint32_t streamLength, HawkeyeStreamInfo *streamInfo);
extern HAWKEYE_STREAM_ERROR HawkeyeDecodeSpectralImage(uint16_t bandNo, uint16_t *pixels, uint32_t pixelLength,
                                                       uint16_t *averageDarkPixels, uint16_t averageDarkPixelLength);
extern HAWKEYE_STREAM_ERROR HawkeyeDecodeFinderscopeImage(uint16_t imageNo, uint16_t *pixels, uint32_t pixelLength);

#endif

#endif // _HAWKEYE_DECODE_

