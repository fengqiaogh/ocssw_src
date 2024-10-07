/*
                HawkeyeDecode.c

                Matt Longmire
                Cloudland Instruments

                Contains the utility functions for scanning and
                decoding a Hawkeye Image Format file (stream)

*/
#include "HawkeyeDecode.h"
#include <string.h>
#include <stdio.h>

#ifdef _DEBUG
        #define DEBUG_LOG       1
#else
        #define DEBUG_LOG       0
#endif

#define HEADER_DATA_LEN 7


/*
        RowInfo:

        Struct for holding the decoded parameters associated with
        a row of pixel data.

*/
typedef struct {
  uint16_t ls5Bits;
  uint16_t width;
  uint16_t rowNumber;
  uint32_t timeStamp;
  uint16_t noDarksAtStart;
  uint16_t isSpectral;
  uint16_t isDark;
  uint16_t dataLen;
  uint16_t electronicGain;
  uint16_t slope1, slope2;
  uint16_t knee;
  uint8_t *pPixels;
} RowInfo;

/*
        HeaderInfo:

        Struct for holding the parameters found in
        the header to each record.

*/

static HawkeyeStreamInfo scannedStreamInfo;
static uint32_t finderscopeOffsets[MAX_FINDERSCOPE_IMAGES];
static uint8_t *scannedStream;
static uint8_t missionLog[MISSION_LOG_LENGTH];

/*      
        ChecksumRecord:

        Calculate and retyurn the 16-bit uint16_t
        checksum of a block.

*/
static uint16_t ChecksumRecord(uint8_t *stream, uint16_t len)
{
  uint16_t sum = 0, recSum;
  int i;

  for (i = 0; i < len; i++)
    sum += stream[i];
  recSum = stream[i] << 8;
  recSum += stream[i + 1];
  return sum == recSum;
}

/*
        FindHeader:

        Scan the stream and return the starting position of the next
        found header. Verifies the checksum and fills in the
        HeaderInfo struct.

*/
static uint8_t HAWKEYE_HEADER_START[] = { 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF };
uint32_t FindHeader(uint8_t *stream, uint32_t streamLength, uint32_t startPos, HeaderInfo *header)
{
  uint32_t i, j;
  uint8_t uc;
  uint32_t ul;
  uint16_t us;

  memset(header, 0, sizeof(HeaderInfo));
  for (i = startPos; i <= streamLength - sizeof(HAWKEYE_HEADER_START); i++)
    {
      uc = stream[i];                 //get the next byte
      if (uc == HAWKEYE_HEADER_START[0]) {                      //does it match the fisrt byte of the header start sequence?
        for (j = 1; j < sizeof(HAWKEYE_HEADER_START); j++) {    //check the full header start sequence length
          if ((i+j) >= scannedStreamInfo.hebBufSize) {
            // Can't find header
            return(0);
          }
          if (stream[i + j] != HAWKEYE_HEADER_START[j])         //does it match?
            break;                    //if not start over at the next byte in the stream
        }
        if (j == sizeof(HAWKEYE_HEADER_START)) {                //have we checked all the header start bytes?
          j += i;                     //move past the start sequence
          if (j + HEADER_DATA_LEN > streamLength)               //we need 7 bytes from the header
            return streamLength;      //too short so bail
          header->blockType = stream[j++];                      //decode the block type
          ul = ((uint32_t)stream[j++]) << 24;                   //decode the record number
          ul += ((uint32_t)stream[j++]) << 16;
          ul += ((uint32_t)stream[j++]) << 8;
          ul += stream[j++];
          ul = (ul & 0x0FFFFFF0) >> 4;//24-bits in 32-bit field
          header->recordNo = ul;
          us = ((uint16_t)stream[j++]) << 8;                    //decode the block len
          us += stream[j++];          //this includes the start, header, data and checksum
          header->blockLen = us;
          if (i + us > streamLength)  //is there enough data left for the payload and checksum?
            return streamLength;      //too short, bail
          if (!ChecksumRecord(stream + i, us - 2))              //does the checksum match?
            continue;                 //if not start at next byte in stream
          return i;                   //otherwise return the count to the 1st byte in the header
        }
      }
    }
  return streamLength;
}

/*
        DecodeImageParams:

        Decode an Image Parameters record, filling in the
        passed HawkeyeImageInfo struct.

*/
static int DecodeImageParams(uint8_t *stream, uint32_t blockLen, uint32_t startPos, HawkeyeImageInfo *pImgInfo)
{
  uint32_t ul;
  uint16_t us;
  uint8_t uc;
  uint64_t ll;
#if INCLUDE_GPS
  uint64_t inBits, ull;
  int i, noInBits, gpsCount;
#endif  
  memset(pImgInfo, 0, sizeof(HawkeyeImageInfo));
#if INCLUDE_GPS
  if (blockLen != 44 + 34 * 4)
    return -1;
#else
  if (blockLen != 44)
    return -1;
#endif
  ul = (uint32_t)stream[startPos++] << 24;      //Decode Exposure ID
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += stream[startPos++];
  pImgInfo->exposureID = ul;
  uc = stream[startPos++];                         //Decode Binnings
  pImgInfo->spectralBinning = uc & 0x0F;
  pImgInfo->finderscopeBinning = (uc >> 4);
  ul = (uint32_t)stream[startPos++] << 24;      //Decode Image ID
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += stream[startPos++];
  pImgInfo->imageID = ul;
  us = (uint16_t)stream[startPos++] << 8;      //Decode Channels
  us += stream[startPos++];
  pImgInfo->channelBitfield = us;
  ll =  (uint64_t)stream[startPos++] << 40;                 //Decode Epoch Time
  ll += (uint64_t)stream[startPos++] << 32;
  ll += (uint64_t)stream[startPos++] << 24;
  ll += (uint64_t)stream[startPos++] << 16;
  ll += (uint64_t)stream[startPos++] << 8;
  ll += stream[startPos++];
  ll = (ll & 0x0FFFFFFFFFF0) >> 4;
  pImgInfo->epochT0 = ll;
  ul = (uint32_t)stream[startPos++] << 24;      //Decode Delta-Host time
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += stream[startPos++];
  pImgInfo->hostDeltaEpoch = (ul >> 4) & 0xFFFFFF;
  ul = (uint32_t)stream[startPos++] << 24;      //Decode Delta-Hawkeye time
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += stream[startPos++];
  pImgInfo->hawkeyeDeltaEpoch = (ul >> 4) & 0xFFFFFF;
  us = (uint16_t)stream[startPos++] << 8;      //Decode CCD 1 Exposure
  us += stream[startPos++];
  pImgInfo->ccd1Exposure = us;
  us = (uint16_t)stream[startPos++] << 8;      //Decode CCD 2 Exposure
  us += stream[startPos++];
  pImgInfo->ccd2Exposure = us;
  us = (uint16_t)stream[startPos++] << 8;      //Decode CCD 3 Exposure
  us += stream[startPos++];
  pImgInfo->ccd3Exposure = us;
  us = (uint16_t)stream[startPos++] << 8;      //Decode CCD 4 Exposure
  us += stream[startPos++];
  pImgInfo->ccd4Exposure = us;
  us = (uint16_t)stream[startPos++] << 8;      //Decode Spectral Height
  us += stream[startPos++];
  pImgInfo->height = us;
  us = (uint16_t)stream[startPos++] << 8;      //Decode Spectral Interval
  us += stream[startPos++];
  pImgInfo->interval = us;
  uc = stream[startPos++];                         //Decode Dark Height / Oversampling
  pImgInfo->oversampling = (uc & 0x03) + 1;                  // Oversampling (0-3 => 1-4)
  pImgInfo->darkHeight = (uc >> 2);        // Dark Height
  if (pImgInfo->darkHeight == 63)          //  only room for 6 bits so 64
    pImgInfo->darkHeight = 64;               //  was encoded as 63
  us = (uint16_t)stream[startPos++] << 8;      //Decode Finderscope Exposure
  us += stream[startPos++];
  pImgInfo->finderscopeExposure = us;
  uc = stream[startPos++];                         //Decode No Finderscope Images
  pImgInfo->noFinderscopeImages = uc;
  uc = stream[startPos++];                         //Decode bitfield
  if ((uc & 03) == 1)
    pImgInfo->compression = DC_PACKED;
  else if ((uc & 03) == 2)
    pImgInfo->compression = DC_DELTA;
  else
    pImgInfo->compression = DC_UNCOMPRESSED;
  if (uc & 0x04)  pImgInfo->darkSubtracted = 1;
  pImgInfo->readoutOrder = (uc & 0x08) ? RO_BLUE_FIRST : RO_GREEN_FIRST;
  pImgInfo->shutterSolenoid = (uc >> 4) & 0x03;      //Shutter solenoid used
  us = (uint16_t)stream[startPos++] << 8;      //Decode Error Code
  us += stream[startPos++];
  pImgInfo->errorCode = us;
#if INCLUDE_GPS
  inBits = 0;
  gpsCount = noInBits = 0;
  for (i = 0; i < 34; i++)                         //Decode 34 32-bit GPS Binary fields
    {
      ul  = (uint32_t)stream[startPos++] << 24;  // Get the next 32-bit field
      ul += (uint32_t)stream[startPos++] << 16;
      ul += (uint32_t)stream[startPos++] << 8;
      ul += stream[startPos++];
      ull = (ul & 0x3FFFFFFC) >> 2;    //extract 28 bits, throwing away leading 0 and trailing 10
      inBits = inBits | (ull << (12 - noInBits));        //add 28-bits to 40-bit bit field
      noInBits += 28;                                  //just added 28 bits
      while (noInBits >= 8) {
        pImgInfo->gpsBinary[gpsCount++] = (uint8_t)(inBits >> 32);   //peel off ms 8 bits
        inBits = (inBits & 0x00FFFFFFFF) << 8;                           //keep ls 32 bits and move to ms bits
        noInBits -= 8;                                                                           //just ate 8 bits
        if (gpsCount == 116)                                                             //have we decoded 116 bytes?
          goto gps_done;                                                                   // if so we're done
      }
    }
gps_done:
#endif
  return 0;
}

/*
        DecodeTelemetry:

        Decode a Telemetry Record, filling in the
        HawkeyeTelemetryInfo struct.

*/
static int DecodeTelemetry(uint8_t *stream, uint32_t blockLen, uint32_t startPos, HawkeyeTelemetryInfo *pTelem)
{
  int i;
  uint32_t ul;
  uint16_t us;
  int noChannels;

  memset(pTelem, 0, sizeof(HawkeyeTelemetryInfo));                   //zero the struct
  noChannels = (blockLen - 4) / 3;                         //calculate no channels reported
  if (noChannels > TC_NO_CHANNELS)                         //is it more than we know about?
    noChannels = TC_NO_CHANNELS;                     // if so limit it
  ul =  (uint32_t)stream[startPos++] << 24;                     //decode the time stamp
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += (uint32_t)stream[startPos++];
  ul = (ul & 0x0FFFFFF0) >> 4;                             //24-bit in 32-bit field
  pTelem->timeStamp = ul;
  for (i = 0; i < noChannels; i++) {                       //for each channel
    us = (uint32_t)stream[startPos++] << 8;               //decode the value
    us += (uint32_t)stream[startPos++];
    pTelem->telemetry.channel[i].channelValue = us;
    us = (uint32_t)stream[startPos++];  //and decode the interpretation
    pTelem->telemetry.channel[i].channelInterp = us;
  }
  pTelem->noChannels = noChannels;                         //record the no channels decoded
  return 0;
}


/*
  DecodeMissionLog:

  Decode a Mission Log Record

*/
static int DecodeMissionLog(uint8_t *stream, uint32_t blockLen, uint32_t startPos, uint8_t *dest)
{
  memset(dest, 0, MISSION_LOG_LENGTH);             //zero the results
  if ( blockLen > MISSION_LOG_LENGTH - 1)          //leave at least 1 byte for NULL terminator
    blockLen = MISSION_LOG_LENGTH - 1;
  memcpy(dest, stream + startPos, blockLen);       //extract text string
  return 0;
}

/*
        ScanRow:

        Scan the Umcompressed or Compressed Pixels record,
        filling in the passed RowInfo struct.  No pixel
        decoding is done here.

*/
static int ScanRow(uint8_t *stream, uint32_t blockLen, uint32_t startPos, RowInfo *pRow)
{
  uint8_t uc, ls5Bits, isSpectral, isDark;
  uint16_t us;
  uint32_t ul;
  int pos = startPos;

  memset(pRow, 0, sizeof(RowInfo));
  if (blockLen < 10)
    return -1;
  uc = stream[pos++];                              //get the flags
  isDark     = (uc & 0x80) ? 1 : 0;//decode dark row flag
  isSpectral = (uc & 0x40) ? 1 : 0;//decode spectral flag
  ls5Bits = (uc & 0x01F);                   //decode ls 5 bits
  if (isSpectral && (ls5Bits < 1 || ls5Bits > 8))
    return -1;
  if (!isSpectral && (ls5Bits < 1 || ls5Bits > FINDERSCOPE_MAX_IMAGES))
    return -1;
  pRow->ls5Bits = ls5Bits;
  pRow->isDark = isDark;
  pRow->isSpectral = isSpectral;
  ul =  (uint32_t)stream[pos++] << 24;          //decode timestamp
  ul += (uint32_t)stream[pos++] << 16;
  ul += (uint32_t)stream[pos++] << 8;
  ul += stream[pos++];
  ul = (ul & 0x0FFFFFF0) >> 4;
  pRow->timeStamp = ul;
  us = (uint16_t)stream[pos++] << 8;           //decode Row Number
  us += stream[pos++];
  pRow->rowNumber = us;
  us = (uint16_t)stream[pos++] << 8;           //decode Row Width
  us += stream[pos++];
  pRow->width = us;
  uc = stream[pos++];                              //decode no darks at start
  pRow->noDarksAtStart = uc;
  uc = stream[pos++];                              //decode electronic gain
  pRow->electronicGain = uc;
  if (isSpectral) {
    us = (uint16_t)stream[pos++] << 8;   //decode slope 1
    us += stream[pos++];
    pRow->slope1 = us;
    us = (uint16_t)stream[pos++] << 8;   //decode slope 2
    us += stream[pos++];
    pRow->slope2 = us;
    us = (uint16_t)stream[pos++] << 8;   //decode knee
    us += stream[pos++];
    pRow->knee = us;
  } else {
    pRow->slope1 = us;
    pRow->slope2 = us;
    pRow->knee = 65535;
  }
  pRow->pPixels = stream + pos;
  pRow->dataLen = (uint16_t)(blockLen - (pos - startPos));
  return 0;
}

/*
        DecodeEOF:

        Decode the End of File record, returning
        the expected file length.

*/
static int DecodeEOF(uint8_t *stream, uint32_t blockLen, uint32_t startPos, int32_t *pExpLen)
{
  uint32_t ul;
  *pExpLen = 0;
  if (blockLen != 4)
    return -1;
  ul = (uint32_t)stream[startPos++] << 24;
  ul += (uint32_t)stream[startPos++] << 16;
  ul += (uint32_t)stream[startPos++] << 8;
  ul += (uint32_t)stream[startPos++];
  *pExpLen = ul;
  return 0;
}

/*
        HawkeyeScanStream:

        Scan the Hawkeye Image File Stream (data bytes) and fill in
        the HawkeyeStreamInfo struct based on what we find. No images
        are stored here, just scanned and detected.
*/
HAWKEYE_STREAM_ERROR HawkeyeScanStream(uint8_t *stream, uint32_t streamLength, HawkeyeStreamInfo *streamInfo)
{
  HAWKEYE_STREAM_ERROR res = HSE_NO_ERROR;
  uint32_t streamPos;
  uint32_t headerPos, nextPos;
  uint32_t blockLen;
  ENCODED_BLOCK_TYPE blockType;
  HawkeyeTelemetryInfo telem;
  HawkeyeImageInfo imgInfo;
  RowInfo rowInfo;
  int blockRes;
  int lastFSI = -1, missingRows;
  int noFSI = 0, band;
  int32_t expectedLen;
  int32_t missing;
  int imageParamsFound = 0;
  uint32_t recNo, nextRecNo = 0;
  HeaderInfo header;
  int firstLight[8];
  int darkHeight[8];
  uint32_t maxBlockLen = 0;
  //ENCODED_BLOCK_TYPE maxBlockType = EBT_NULL;

#if DEBUG_LOG
  FILE *fp;
  fp = fopen("HawkeyeDecodeLog.txt", "wt");
#endif

  memset(streamInfo, 0, sizeof(HawkeyeStreamInfo));                //zero the results
  memset(finderscopeOffsets, 0, sizeof(finderscopeOffsets));       //and where we found each finderscope image
  memset(firstLight, 0, sizeof(firstLight));                       //init the dark counters
  memset(darkHeight, 0, sizeof(darkHeight));
  scannedStream = stream;                                          //remember the stream for when we extract components later
  streamInfo->streamLength = streamLength;
  streamInfo->hebBufSize = streamLength;
  scannedStreamInfo.hebBufSize = streamLength;
  streamInfo->noMissingBytes = -1;

  // Find the first valid Image Parameters
  streamPos = 0;
  headerPos = FindHeader(stream, streamLength, streamPos, &header);//find the first valid header
  while (streamPos < streamLength) {                               //process until we are done
    blockType = header.blockType;                                  //remember what record type we found
    blockLen = header.blockLen;                                    //remember the length of the records data block
    streamPos = headerPos + sizeof(HAWKEYE_HEADER_START)+HEADER_DATA_LEN;      //skip passed this recordheader
    nextPos = FindHeader(stream, streamLength, streamPos, &header);            //and find the next
    if (blockType != EBT_EOF && nextPos - headerPos != blockLen) {             //is the next headr at the right position or is this an EOF
      streamPos = headerPos = nextPos;                                         // and try the next record
      continue;
    }
    blockLen -= 9 + sizeof(HAWKEYE_HEADER_START);                            //skip the record headr bytes
    if (blockType == EBT_IMAGE_PARAMS) {
      if ((blockRes = DecodeImageParams(stream, blockLen, streamPos, &imgInfo)) == 0) {  //if so decode it
        streamInfo->imageInfo = imgInfo;                                                                 //and record it
        imageParamsFound = 1;
        break;
      }
    }
    streamPos = headerPos = nextPos;                                                         //skip to next header
  }
  if ( !imageParamsFound )
    res = HSE_NO_IMAGE_PARAMETERS;
  else {
    stream = scannedStream;                                                                          //start over at the beginning
    streamPos = 0;
    headerPos = FindHeader(stream, streamLength, streamPos, &header);                  //find the first valid header
    while (streamPos < streamLength) {                                                       //process until we are done
      blockType = header.blockType;                                                    //remember what record type we found
      recNo = header.recordNo;                                                                 //and the record number
      if (recNo > nextRecNo)                                                                   //have we missed any records?
        streamInfo->noMissingRecords += recNo - nextRecNo;       // if so keep track of that
      nextRecNo = recNo + 1;                                                                   //the next record number we expect to find
      blockLen = header.blockLen;                                                              //remember the length of the records data block
      streamPos = headerPos + sizeof(HAWKEYE_HEADER_START)+HEADER_DATA_LEN;   // sko passed this recordheader
      nextPos = FindHeader(stream, streamLength, streamPos, &header);            //and find the next
      if (blockType != EBT_EOF && nextPos - headerPos != blockLen) {             //is the next headr at the right position or is this an EOF
        streamInfo->noBadRecords++;                                                      // if no recod this a a bad record
        streamPos = headerPos = nextPos;                                         // and try the next record
        continue;
      }
      streamInfo->noRecords = recNo + 1;                                               //we've found a good record
      if (blockLen > maxBlockLen)     {                                                        //keep track of the longest block length
        maxBlockLen = blockLen;
        //maxBlockType = blockType;
      }
      blockLen -= 9 + sizeof(HAWKEYE_HEADER_START);                    //skip the record headr bytes
      blockRes = 0;                                                    //assume it's a good block
      switch (blockType) {                                             //what type of block did we find?
      case EBT_UNCOMPRESSED:                                           //Uncompressed or
      case EBT_COMPRESSED:                                             //Compressed Pixel Data?
        if ((blockRes = ScanRow(stream, blockLen, streamPos, &rowInfo)) == 0) {                    //if so scan the record
          if (!rowInfo.isSpectral) {                                                                               //is this Finderscope data?
            noFSI = rowInfo.ls5Bits - 1;                                                             // if so get the finderscope image no
            if (noFSI < MAX_FINDERSCOPE_IMAGES) {                                            //do we have room for it?
              if (noFSI != lastFSI) {                                                                  //have we jumped to the next image?
#if DEBUG_LOG
                fprintf(fp, "Finderscope Image %d\n", noFSI + 1);
#endif
                finderscopeOffsets[noFSI] = headerPos;
                //                printf("finderscopeOffsets: %d\n", headerPos);
                  
                //remember it's position so it's faster to retrieve
                streamInfo->noFinderscopeImages = noFSI + 1;             //count it
                //                printf("# FS2: %d\n", streamInfo->noFinderscopeImages);
                streamInfo->finderscopeInfo[noFSI].timeStamp = rowInfo.timeStamp;  //and record the timestamp
              }
              if (streamInfo->finderscopeInfo[noFSI].width == 0) {     //have we detected the width yet?
                streamInfo->finderscopeInfo[noFSI].width = rowInfo.width;                  //if not record it and the electronic gain
                streamInfo->finderscopeInfo[noFSI].electronicGain = rowInfo.electronicGain / 10.0;
              }
              missingRows = rowInfo.rowNumber - streamInfo->finderscopeInfo[noFSI].height;
              if (missingRows > 0)                                                                     //any missing rows?
                streamInfo->finderscopeInfo[noFSI].noMissingRows += missingRows;   //if so count them
              if (missingRows >= 0)                                                                    //is this a new row
                streamInfo->finderscopeInfo[noFSI].height = rowInfo.rowNumber + 1; // if so add it to the height
            }
            lastFSI = noFSI;                                                       //remember this image number
          }
          else {                                                                   //otherwise this is Spectral Band pixel data
#if DEBUG_LOG
            fprintf(fp, "Spectral Band %d %s Row %d\n", rowInfo.ls5Bits,
                    rowInfo.isDark ? "Dark " : "Light",
                    rowInfo.rowNumber);
#endif
            band = rowInfo.ls5Bits - 1;                                            //remember which band it came from
            if (streamInfo->spectralInfo[band].width == 0) {                       //have we not seen this band before
              streamInfo->spectralInfo[band].width = rowInfo.width;                //if so record the width
              streamInfo->spectralInfo[band].timeStamp = rowInfo.timeStamp;        //and the timestamp and electronic gain
              streamInfo->spectralInfo[band].electronicGain = rowInfo.electronicGain / 100.0;
            }
            if (rowInfo.isDark)                                                    //is it a dark row?
              streamInfo->spectralInfo[band].averageDarkRowPresent = 1;            //let them know we found it
            else {
              if (streamInfo->spectralInfo[band].height == 0)                      //is this the first light we've seen
                firstLight[band] = rowInfo.rowNumber;                              // if so remember it
              missingRows = rowInfo.rowNumber - streamInfo->spectralInfo[band].height;
              if (missingRows > 0)                                                 //did we miss any rows?
                streamInfo->spectralInfo[band].noMissingRows += missingRows;       //if so keep track
              if (missingRows >= 0)                                                //is this a new row?
                streamInfo->spectralInfo[band].height = rowInfo.rowNumber + 1;     //if so it adds to the height
            }
          }
        };
        break;
      case EBT_IMAGE_PARAMS:                                                       //is the record Image Parameters?
        if ((blockRes = DecodeImageParams(stream, blockLen, streamPos, &imgInfo)) == 0) {  //if so decode it
#if DEBUG_LOG
          fprintf(fp, "Image Parameters\n");
#endif
        }
        break;
      case EBT_TELEMETRY:                                                          //is the record Telemetry Data
        if ((blockRes = DecodeTelemetry(stream, blockLen, streamPos, &telem)) == 0) {      //if so Decode it
#if DEBUG_LOG
          fprintf(fp, "Telemetry Record\n");
#endif
          if (streamInfo->noTelemetryRecords < MAX_TELEMETRY_RECORDS) {    //do we have room for it?
            streamInfo->telemetryInfo[streamInfo->noTelemetryRecords] = telem;                 //if so record it
            streamInfo->noTelemetryRecords++;                                                        //and count it
          }
        }
        break;
      case EBT_MISSION_LOG:                                                        //is the record the Mission Log?
        if ((blockRes = DecodeMissionLog(stream, blockLen, streamPos, missionLog)) == 0) { //if so Decode it
#if DEBUG_LOG
          fprintf(fp, "Mission Log Record\n");
#endif
          if ( strlen(streamInfo->missionLog) == 0)                                //do we have it already
            memcpy(streamInfo->missionLog, missionLog, MISSION_LOG_LENGTH);        //if not record it
        }
        break;
      case EBT_EOF:                                                                //is the record and EOF?
        if ((blockRes = DecodeEOF(stream, blockLen, streamPos, &expectedLen)) == 0) {      //if so decode it
#if DEBUG_LOG
          fprintf(fp, "EOF Record\n");
#endif
          missing = expectedLen - (streamPos + 6);
          if (missing >= 0)                                                        //any bytes missing
            streamInfo->noMissingBytes = missing;                                  //  if so record them
          streamInfo->streamLength = expectedLen;                                  //and record the expected length
        }
        break;
      default:                                                                     //it was an unknown record
#if DEBUG_LOG
        fprintf(fp, "Unknown Record\n");
#endif
        streamInfo->noUnknownRecords++;                                            //count it as such
        break;
      }
      if (blockRes != 0)                                                           //was there an error in decoding the block
        streamInfo->noBadRecords++;                                                // if so count it as such
      streamPos = headerPos = nextPos;
      if (blockType == EBT_EOF)
        break;
    } // while not EOF or more bytes to process
    for (band = 0; band < 8; band++) {                               //for each of the spectral bands
      if (streamInfo->spectralInfo[band].height != 0) {              //did we see any row data?
        streamInfo->noSpectralImages++;                  //if so then count it
      }
    }
    scannedStreamInfo = *streamInfo;                                 //remember the stream info so we return so we can reuse it
  }
#if DEBUG_LOG
  fclose(fp);
#endif
  return res;
}

/*
        UncompressKLIPixels:

        Uncompress a row of Spectral Pixels that were compressed
        using the dual-slope + knee method, then otionally
        delta encoded and then finally packed.
*/
static void UncompressKLIPixels(uint16_t *pixels, uint8_t *stream, RowInfo rowInfo)
{
  uint16_t i, p0=0, p1, p2, vid;
  int d1, d2;
  uint32_t bits=0;
  int noBits;
  uint32_t lvid, s1, s2, v1;
  int width = rowInfo.width;
  uint32_t eGain;

  // Init compression params
  s1 = rowInfo.slope1 * 10;          //1st slope, e-/ADU x 100
  s2 = rowInfo.slope2 * 10;          //2nd slope, e-/ADU x 100
  eGain = rowInfo.electronicGain; // KLI EGain, e-/ADU x 100
  v1 = rowInfo.knee;                         //Knee ADU

  if (scannedStreamInfo.imageInfo.compression == DC_PACKED) {        //is the data jsut packed (no delta compression)
    // Unpack the data 12-bit fields in the 8-bit byte stream
    noBits = 0;
    for (i = 0; i < width;) {                                        //for each pixel in the row
      while (noBits < 12) {                                    //make sure we have at least 12 bits
        bits = bits | (*stream++ << (24 - noBits));                //if not get 8 more bits
        noBits += 8;                                             //and keep count
      }
      vid = bits >> 20;                                                //peel off the 12-bit pixle value
      bits = bits << 12;                                               //keep the others
      noBits -= 12;                                                    //keep count
      pixels[i++] = vid & 0xFFF;                               //and store the pixel
    }
  } else {                                                                                 //other wise the data is packed and delat encoded
    // Unpack the data 13-bit fields in the 8-bit byte stream
    noBits = 0;
    for (i = 0; i < width;) {                                        //for each pixel in the row
      while (noBits < 13) {                                    //make sure we have at least 13 bits
        bits = bits | (*stream++ << (24 - noBits));                //if not get another 8 bits
        noBits += 8;                                             //and keep count
      }
      vid = bits >> 19;                                                //peel off the 13-bit encoded value
      bits = bits << 13;                                               //keep the others
      noBits -= 13;                                                    //and keep count
      if (vid & 0x1000) {                                              //is the flag bit (msb) set?
        p0 = vid & 0xFFF;                                        // if so it's just 1 12-bit pixle value
        pixels[i++] = p0;                                        // so save it
      } else {                                                                 //otherwise the flag bit (msb) is clear
        d1 = (vid >> 6) & 0x3F;                          //so we have two 6-bit delats in the lower 12-bits
        if (d1 & 0x20)  d1 = -64 + d1;           //sign exten the first delat
        p1 = p0 + d1;                                            //calculate the resulting pixel value
        d2 = (vid & 0x3F);                                       //extract the 2nd 6-bit delta
        if (d2 & 0x20)  d2 = -64 + d2;           //sign extend it
        p0 = p2 = p1 + d2;                                       //calculate the 2nd pixel value
        pixels[i++] = p1;                                        //store the first
        pixels[i++] = p2;                                        //and the 2nd
      }
    }
  }

  // decompress the data
  for (i = 0; i < width; i++) {                                    //for each pixel in the row
    lvid = pixels[i];                                                        //get the dual-slope encoded value
    lvid = (lvid * s1) / eGain;                                      //decode the 1st slope encoding
    if (lvid > v1) {                                                         //was it above the knee
      lvid -= v1;                                                              // if so subtract the knee
      lvid = lvid * s2 / s1;                                   // decode the 2nd slope encoding
      lvid += v1;                                                              // and add the 2nd slope offset
    }
    pixels[i] = (uint16_t)(lvid > 65535 ? 65535 : lvid); //finally clip and save the decoded pixel value
  }
}

/*
        HawkeyeDecodeSpectralImage:

        Decode a Spectral Band image, saving the data as an array of
        uint16_ts in to the passed buffer.

        For each row, the first two unsigend shorts represent the time
        code offset from start of image with the ms 16-bits first
        followed by the ls 16-bits.

        After the time code offset the row pixel data is stored.

*/
HAWKEYE_STREAM_ERROR HawkeyeDecodeSpectralImage(uint16_t bandNo, uint16_t *pixels, uint32_t pixelLength,
                                                uint16_t *averageDarkPixels, uint16_t averageDarkPixelLength)
{
  HAWKEYE_STREAM_ERROR res = HSE_NO_ERROR;
  uint32_t streamPos;
  uint32_t headerPos, nextPos;
  uint32_t blockLen;
  ENCODED_BLOCK_TYPE blockType;
  RowInfo rowInfo;
  uint32_t pixelPos = 0;
  int i;
  uint16_t pix;
  uint8_t *pPix;
  HeaderInfo header;
  uint16_t darkStored = 0;

  memset(pixels, 0, 2 * pixelLength);                                                                //zero the results
  streamPos = 0;
  headerPos = FindHeader(scannedStream, scannedStreamInfo.streamLength, streamPos, &header);         //find the first record header
  while (streamPos < scannedStreamInfo.streamLength) {                                               //while there's bytes left
    blockType = header.blockType;                                                                    //remember the block type
    blockLen = header.blockLen;                                                                      //and length
    streamPos = headerPos + sizeof(HAWKEYE_HEADER_START)+HEADER_DATA_LEN;            //move passed this record header
    nextPos = FindHeader(scannedStream, scannedStreamInfo.streamLength, streamPos, &header);   //and find the next
    if (nextPos<=streamPos) { 																// avoid infinite loop, LH, Jan. 19, 2021
      printf("HawkeyeDecodeSpectralImage error.\n");
      break;
    }
    if (nextPos - headerPos != blockLen) {                                                     //is the next header at the right postiton?
      streamPos = headerPos = nextPos;                                                         // if not then skip this one
      continue;
    }
    blockLen -= 9 + sizeof(HAWKEYE_HEADER_START);                                              //account fo the header bytes
    if ( blockType== EBT_UNCOMPRESSED || blockType == EBT_COMPRESSED ) {             //did we find a Pixel Data record?
      if (ScanRow(scannedStream, blockLen, streamPos, &rowInfo) == 0) {                //if so then scan it
        if (rowInfo.isSpectral && rowInfo.ls5Bits - 1 == bandNo) {                       //is it Spectral Data and the correct Band?
          if (rowInfo.isDark) {                                                                                    //is it dark data?
            if (averageDarkPixels != 0 && !darkStored) {                             //is there a buffer and have we not stored it already?
              if (rowInfo.width + 2 <= averageDarkPixelLength ) {            //will it fit?
                pixelPos = 0;                                                //only one row of average darks, at the start of the passed buffer
                pPix = rowInfo.pPixels;                                      //get the pointer to the stream's pixel data
                //                printf( "Dark Data: %d Saving Time Stamp\n", streamPos);
                pixels[pixelPos++] = (uint16_t)(rowInfo.timeStamp >> 16);    //save the time code ms 16-bits
                pixels[pixelPos++] = (uint16_t)(rowInfo.timeStamp);                  // and ls 16-bits
                if (blockType == EBT_UNCOMPRESSED || rowInfo.dataLen ==                    //is the data uncompressed
                    2 * rowInfo.width) {                                                     //or stored as uncompressed?
                  for (i = 0; i < rowInfo.width; i++) {                    // if so then for every pixel
                    pix = *pPix++ << 8;                                                      // get the ms byte
                    pix += *pPix++;                                                          // and ls byte
                    averageDarkPixels[i + pixelPos] = pix;           // and save the pixel value
                  }
                }
                else {
                  printf( "Dark Data Uncompressing\n");
                  UncompressKLIPixels(averageDarkPixels + pixelPos, pPix, rowInfo);          //so uncompress it
                }
                darkStored = 1;
              } // pixels would fit in buffer
              else
                break;
            }
          }
          else {
            // Non-dark data
            //printf("row number: %d\n", rowInfo.rowNumber);
            pixelPos = rowInfo.rowNumber * (scannedStreamInfo.spectralInfo[bandNo].width + 2);
            if (pixelPos + rowInfo.width + 2 <= pixelLength) {                                       //is there room for the whole row?
              pPix = rowInfo.pPixels;                                                      //get the pointer to the stream's pixel data
              //              printf( "%d Saving Time Stamp\n", streamPos);
              pixels[pixelPos++] = (uint16_t)(rowInfo.timeStamp >> 16);                    //save the time code ms 16-bits
              pixels[pixelPos++] = (uint16_t)(rowInfo.timeStamp);        // and ls 16-bits
              if (blockType == EBT_UNCOMPRESSED || rowInfo.dataLen == 2 * rowInfo.width) {// is the data uncompressed
                for (i = 0; i < rowInfo.width; i++) {                                    // if so then for every pixel
                  pix = *pPix++ << 8;                                                                      // get the ms byte
                  pix += *pPix++;                                                                          // and ls byte
                  pixels[i + pixelPos] = pix;                                                      // and save the pixel value
                }
              }
              else {
                //otherwise the data is compressed
                //                printf( "Uncompressing\n");
                UncompressKLIPixels(pixels + pixelPos, pPix, rowInfo);   //so uncompress it
              }
            }  // pixels would fit into buffer
            else
              break;
          }
        } // is the correct band
      } // row scanned OK
    } // block is pixel data
    streamPos = headerPos = nextPos;                                     //go the the next header
    if (blockType == EBT_EOF)                                            //did we find an EOF
      break;                                                             //if so we are done
  }
  return res;
}

/*
        HawkeyeDecodeFinderscopeImage:

        Decode and save an image from the Finderscope a an
        array of unsigend shorts into the passed buffer.
*/
HAWKEYE_STREAM_ERROR HawkeyeDecodeFinderscopeImage(uint16_t imageNo, uint16_t *pixels, uint32_t pixelLength)
{
  HAWKEYE_STREAM_ERROR res = HSE_NO_ERROR;
  uint32_t streamPos;
  uint32_t headerPos, nextPos;
  uint32_t blockLen;
  ENCODED_BLOCK_TYPE blockType;
  RowInfo rowInfo;
  uint32_t pixelPos = 0;
  int i;
  uint16_t pix;
  uint8_t *pPix;
  int fsiNo;
  uint32_t bits, pv;
  int noBits;
  HeaderInfo header;

  memset(pixels, 0, 2 * pixelLength);                                 //xero the results
  streamPos = finderscopeOffsets[imageNo];                            //start at the previously found location
  headerPos = FindHeader(scannedStream, scannedStreamInfo.streamLength, streamPos, &header);    //and fid the first record header
  while (streamPos < scannedStreamInfo.streamLength) {                //while there's bytes left to process
    blockType = header.blockType;                                     //remember the block type
    blockLen = header.blockLen;                                       //and length
    streamPos = headerPos + sizeof(HAWKEYE_HEADER_START)+HEADER_DATA_LEN;                       //skip passed this header
    nextPos = FindHeader(scannedStream, scannedStreamInfo.streamLength, streamPos, &header);    //and find the next
    if (nextPos == 0) {
      if (headerPos == FindHeader(scannedStream, scannedStreamInfo.streamLength, finderscopeOffsets[imageNo], &header)) {
        printf("No FinderScope Header Found.\n");					  // only when the loop has not found any data yet
        return(HSE_NO_HEADER_FOUND);								  // LH, 4/19/2022, l1agen_hawkeye.cpp, v1.0.4
      } else {
        return res;
      }
    }
    if (nextPos - headerPos != blockLen ) {                           //is the next in the right position?
      streamPos = headerPos = nextPos;                                // if not ship this record
      continue;
    }
    blockLen -= 9 + sizeof(HAWKEYE_HEADER_START);                     //account for the record header
    if (blockType == EBT_UNCOMPRESSED || blockType == EBT_COMPRESSED) {                         //did we find Pixel data?
      if (ScanRow(scannedStream, blockLen, streamPos, &rowInfo) == 0) {                         //if so scan the record
        if (!rowInfo.isSpectral) {                                    //is it Finderscope data
          fsiNo = rowInfo.ls5Bits - 1;                                //decode the image number
          if (fsiNo == imageNo) {                                     //is that the one we want?
            pixelPos = rowInfo.rowNumber * scannedStreamInfo.finderscopeInfo[imageNo].width;
            if (pixelPos + rowInfo.width <= pixelLength) {            //is there room for al the pixels?
              pPix = rowInfo.pPixels;                                 //get a poiner to the record's pixel data
              if (blockType == EBT_UNCOMPRESSED) {                    //is it uncompressed data
                for (i = 0; i < rowInfo.width; i++) {                 //  if so then for each pixel
                  pix = *pPix++ << 8;                                 //  get the ms 8-bits
                  pix += *pPix++;                                     //  and the ls 8-bitd
                  pixels[i+pixelPos] = pix;                           //  and save the pixel value
                }
              } else {                                                //otherwise the pixel data is packed 
                bits = 0;
                noBits = 0;
                for (i = 0; i < rowInfo.width; i++) {                 //for each pixel in the row
                  while (noBits < 10) {                               //we need at least 10-bits
                    pv = *pPix++;                                     // if we dont have it
                    bits = bits | (pv << (24 - noBits));              // then get 8 more
                    noBits += 8;                                      // and account for it
                  }
                  pix = (bits >> 22);                                 //extract the 10-bit pixel value
                  bits = bits << 10;                                  //keep the others
                  noBits -= 10;                                       //account for it
                  pixels[i + pixelPos] = pix;                         //and save the pixel value
                }
              } // compressed data
            }  // pixels would fit into buffer
            else
              break;
          } // is the correct image no
        } // is finderscope data
      } // row scanned OK
    } // block is pixel data
    streamPos = headerPos = nextPos;                                  //go to the next record header
    if (blockType == EBT_EOF)                                         //did we find an EOF?
      break;                                                          //if so we're done
  }
  return res;
}
