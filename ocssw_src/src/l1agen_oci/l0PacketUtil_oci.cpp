#include <algorithm>
#include <iterator>

#include <genutils.h>
#include "l0PacketUtil_oci.h"

using namespace std;

/// Given an ancillary packet, extract the number of pixels, bands and tap information from it.
/// The parameters are references
int getBandDimensions(uint8_t *ancillaryPacket, uint16_t &numCcdPixels, uint16_t &numBlueBands,
                      uint16_t &numRedBands, uint16_t &numSwirPixels, uint16_t &numDarkCcdPixels,
                      uint16_t &numDarkSwirPixels, uint16_t *isBlueCcdTapsEnabled,
                      uint16_t *isRedCcdTapsEnabled, spatialAggTable *spatialAggList) {
    // Extract spatial aggregation table and compute numbers of pixels
    short ioff = 36;
    numCcdPixels = 0;
    numSwirPixels = 0;
    numDarkCcdPixels = 0;
    numDarkSwirPixels = 0;

    for (size_t i = 0; i < 10; i++) {
        spatialAggList[i].dataType = ancillaryPacket[ioff + 3] % 16;
        spatialAggList[i].spatialAggCode = ancillaryPacket[ioff + 2] % 4;
        if (spatialAggList[i].dataType == 5)
            spatialAggList[i].spatialAggCode = 0;
        spatialAggList[i].lines = ancillaryPacket[ioff] * 256 + ancillaryPacket[ioff + 1];
        ioff += 4;

        if (spatialAggList[i].dataType > 0 && spatialAggList[i].dataType <= 14 &&
            spatialAggList[i].dataType != 10) {  // Changed dtype<=12 to 14, LH, 3/30/2022
            if (spatialAggList[i].dataType == 2) {
                numDarkCcdPixels +=
                    spatialAggList[i].lines / SPATIAL_AGG_CODES[spatialAggList[i].spatialAggCode];
                numDarkSwirPixels += spatialAggList[i].lines / 8;
            }
            numCcdPixels += spatialAggList[i].lines / SPATIAL_AGG_CODES[spatialAggList[i].spatialAggCode];
            numSwirPixels += spatialAggList[i].lines / 8;
        }
    }

    // to ensure that the science and dark count arrays will be created with at least one pixel
    if (numCcdPixels == 0) {
        numCcdPixels = 1;
        numSwirPixels = 1;
    }

    if (numDarkCcdPixels == 0) {
        numDarkCcdPixels = 1;
        numDarkSwirPixels = 1;
    }

    ioff += 4;

    // Extract spectral aggregation and compute numbers of bands
    // Tap enable flags
    short blueTapFlag = ancillaryPacket[ioff + 2] * 256 + ancillaryPacket[ioff + 3];
    short redTapFlag = ancillaryPacket[ioff] * 256 + ancillaryPacket[ioff + 1];

    // Tap aggregation factors
    uint32_t blueSpectralAggFactor;
    uint32_t redSpectralAggFactor;
    memcpy(&blueSpectralAggFactor, &ancillaryPacket[ioff + 8], sizeof(uint32_t));
    memcpy(&redSpectralAggFactor, &ancillaryPacket[ioff + 4], sizeof(uint32_t));
    swapc_bytes((char *)&blueSpectralAggFactor, sizeof(uint32_t), 1);
    swapc_bytes((char *)&redSpectralAggFactor, sizeof(uint32_t), 1);

    if (spatialAggList[1].dataType == 1 && (blueSpectralAggFactor != BASE_BLUE_SPECTRAL_AGG_FACTOR ||
                                            redSpectralAggFactor != BASE_RED_SPECTRAL_AGG_FACTOR)) {
        spatialAggList[1].dataType = 13;
    }

    // Compute number of bands for enabled taps
    // Bands are in reverse spectral order
    numBlueBands = 0;
    numRedBands = 0;

    // there is a total of 16 taps, currTapBit is used as bitwise& to check each tap from the least sig bit to
    // the most significant (right to left). After each tap iteration, multiply by 2 to get the next bit
    // 1  = 0000 0000 0000 0001 : checks the right most bit for blue or red tap and checks if it is enabled
    // 2  = 0000 0000 0000 0010
    // 4  = 0000 0000 0000 0100
    // ...
    // 65536 = 1000 0000 0000 0001
    uint16_t currTap = 1;

    // works the same as currTapBit, except you are extracting the aggrigation factor which are 2 bits and bc
    // the aggrigation factors are 32 bits. Each iteration gets multiplied by 4.
    // 3  = 0000 0000 0000 0000 0000 0000 0000 0011
    // 12 = 0000 0000 0000 0000 0000 0000 0000 1100
    // 48 = 0000 0000 0000 0000 0000 0000 0011 0000
    // ...
    // 12884901888 = 1100 0000 0000 0000 0000 0000 0000 0000
    uint32_t currAggFactorBit = 3;

    // works with currAggBit. After grabbing the aggigration factor, divide that number by this to extract the
    // bits. This will be factors of 4, so it will move 2 bits to the right
    // so if we have the aggigration factor at 12, the 2nd iteration, extract bit is 4.12/4 = 3
    // 12 = 0000 0000 0000 0000 0000 0000 0000 1100
    // 3  = 0000 0000 0000 0000 0000 0000 0000 0011 : after shifting it 2 bytes to the right
    uint32_t currAggFactorExtractBit = 1;

    for (int i = 15; i >= 0; i--) {
        isBlueCcdTapsEnabled[i] = (blueTapFlag & currTap) / currTap;
        if (isBlueCcdTapsEnabled[i])
            numBlueBands +=
                32 / SPATIAL_AGG_CODES[(blueSpectralAggFactor & currAggFactorBit) / currAggFactorExtractBit];
        isRedCcdTapsEnabled[i] = (redTapFlag & currTap) / currTap;
        if (isRedCcdTapsEnabled[i])
            numRedBands +=
                32 / SPATIAL_AGG_CODES[(redSpectralAggFactor & currAggFactorBit) / currAggFactorExtractBit];
        currTap *= 2;
        currAggFactorBit *= 4;
        currAggFactorExtractBit *= 4;
    }

    return 0;
}

int getAncillaryPacketTime(uint8_t *ancillaryPacket, int32_t &year, int32_t &day, double &startTime) {
    // Unpack and convert the CCSDS segmented time code
    // from the OCI ancillary packet

    // time offset from start of ancillary packet that gives the time
    short int timeOffSet = 28;

    // Get day count since Jan. 1, 1958 (Julian day 2436205)
    // sec58 = seconds since 0 UT Jan. 1, 1958
    uint32_t secSinceJan1958Tai;
    memcpy(&secSinceJan1958Tai, &ancillaryPacket[timeOffSet], sizeof(uint32_t));
    swapc_bytes((char *)&secSinceJan1958Tai, sizeof(uint32_t), 1);

    int leap = leapseconds_since_1993((double)secSinceJan1958Tai) + 27;
    secSinceJan1958Tai -= leap;

    // Convert to year and day
    day = secSinceJan1958Tai / 86400;
    int32_t julianDay = day + 2436205;  // Jan. 1, 1958 is Julian day 2436205
    jdate(julianDay, &year, &day);

    // Get microseconds
    uint32_t startTimeMicroseconds;
    memcpy(&startTimeMicroseconds, &ancillaryPacket[timeOffSet + 4], sizeof(uint32_t));
    swapc_bytes((char *)&startTimeMicroseconds, sizeof(uint32_t), 1);
    startTimeMicroseconds = startTimeMicroseconds / 4096;  // 20 MSBs used, last 12 are spares

    // get seconds and add microseconds
    uint32_t isec = secSinceJan1958Tai % 86400;
    startTime = isec + ((double)startTimeMicroseconds) * 1e-6;

    return 0;
}

/**
 * @brief Get spin number from the packet based on the apid
 * @param apid
 * @param packet
 * @return byte swapped spinNum from the packet
 */
int32_t getSpinNumFromPacket(uint32_t apid, uint8_t *packet) {
    int32_t spinNum = 0;

    // Ancillary Packet APID
    if (apid == ANCILLARY_APID) {
        memcpy(&spinNum, &packet[24], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);

        // APIDs that only have spin number
    } else if (apid == APIDS_WITH_SPIN[0] || apid == APIDS_WITH_SPIN[1]) {
        // Science Packet without time field, ignore SWIR packets for now
        memcpy(&spinNum, &packet[6], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);

        // APIDS that have both time field and spin number field
    } else if (std::find(std::begin(APIDS_WITH_TIME_SPIN), std::end(APIDS_WITH_TIME_SPIN), apid) !=
               std::end(APIDS_WITH_TIME_SPIN)) {
        // Packet with time field
        memcpy(&spinNum, &packet[12], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);
    }

    return spinNum;
}

// Get spin number from packet and save it to to appropriate spin var
void getSpinNumFromPacket(int32_t &spinNum, int32_t &spinNumTelemetry, uint32_t apid, uint8_t *packet) {
    if (apid == ANCILLARY_APID) {
        memcpy(&spinNum, &packet[24], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);

    } else if (apid == APIDS_WITH_SPIN[0] || apid == APIDS_WITH_SPIN[1]) {  // ver 1.00.01
        memcpy(&spinNum, &packet[6], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);

    } else if (std::find(std::begin(APIDS_WITH_TIME_SPIN), std::end(APIDS_WITH_TIME_SPIN), apid) !=
               std::end(APIDS_WITH_TIME_SPIN)) {
        memcpy(&spinNumTelemetry, &packet[12], 4);
        swapc_bytes((char *)&spinNumTelemetry, sizeof(uint32_t), 1);
    }
}

// get the next packets apid and length without moving the streams pointer
int getNextPacketInfo(L0Stream *l0FileStream, uint32_t &packetLen, uint32_t &apid, int32_t &isEndFile,
                      bool isSPW) {
    // no more packets to read
    if (l0FileStream->eof()) {
        isEndFile = 1;
        packetLen = 0;
        cout << "End of packet file" << endl;
        return 0;
    }

    // Read spacewire header if needed
    if (isSPW) {
        uint8_t spwhead[2];
        l0FileStream->read((char *)&spwhead, 2);
    }

    // Read packet header
    uint8_t phead[6];

    l0FileStream->read((char *)&phead, 6);

    // Get length of packet body and APID
    packetLen = phead[4] * 256 + phead[5] + 1 + 6;
    apid = (phead[0] % 8) * 256 + phead[1];

    if (l0FileStream->tellg() == -1)
        isEndFile = 1;
    l0FileStream->seekg(-6, ios_base::cur);
    if (isSPW)
        l0FileStream->seekg(-2, ios_base::cur);
    return 0;
}

// Read all packets with the same spin number from the first packet and loads it into the all packet buffer
// The main l1agen_oci program will unpack science data based on the first packets spin number and the total
// number of packets that share the same spin number
int readScanPackets(L0Stream *l0FileStream, uint8_t *packet, uint8_t (*allPacketsBuffer)[MAX_PACKET_SIZE],
                    uint32_t &numPackets, int32_t &firstPacketSpinNum, int32_t &ancillaryIndex,
                    vector<int32_t> &telemetryIndices, uint8_t &sequenceErrorFlag, int32_t &isEndFile,
                    bool isSPW) {
    if (isEndFile)
        return 0;

    // these are referenced variables for the main l1agen_oci loop. reset them to store new data from
    // reading all packets that have the firstPacket's spin number
    telemetryIndices.clear();
    ancillaryIndex = -1;
    numPackets = 0;

    // initial packet length and apid from the passed in argument, packet.
    uint32_t packetLen = packet[4] * 256 + packet[5] + 7;
    uint32_t currPacketApid = (packet[0] % 8) * 256 + packet[1];

    /// update the firstPacketSpinNum to reflect the first packets spin number that was passed in.
    /// all packets after this packet should have the same spin number.
    firstPacketSpinNum = getSpinNumFromPacket(currPacketApid, packet);

    // track spin number for each packet locally when going through the packets
    // do not modify firstPacketSpinNum since it will be used in main l1agen_oci loop
    int32_t currPacketSpinNum = firstPacketSpinNum;
    int32_t currPacketSpinNumTelemetry = firstPacketSpinNum;  // always 1 less than than the initial spinNum
    // uint8_t *packet = packet;

    // science packet sequence number of the previous packet. used to make make sure sequence is in order
    int prevSequenceNum = 0;

    // Read all of the packets that have the same spin number as the first packet and store it in the all
    // packets buffer in the order that it was seen
    while (currPacketSpinNum <= firstPacketSpinNum &&
           currPacketSpinNumTelemetry <= (firstPacketSpinNum + 1) && (numPackets < MAX_NUM_PACKETS)) {
        // copy packet data into the all packets buffer
        memcpy(&allPacketsBuffer[numPackets][0], packet, packetLen);
        currPacketApid = (packet[0] % 8) * 256 + packet[1];

        // Check science packet sequence numbers and check for out of order, LH 11/20/2020
        if (currPacketApid == APIDS_WITH_SPIN[0]) {
            int currSequenceNum =
                (allPacketsBuffer[numPackets][2] % 64) * 256 + allPacketsBuffer[numPackets][3];

            // if previous sequence is not the first and if the difference between the current and last seq
            // number is not 1, then there is a sequence jump, report the error.
            if ((prevSequenceNum > 0) &&
                ((currSequenceNum - prevSequenceNum) != 1 && (currSequenceNum - prevSequenceNum) != -16383)) {
                sequenceErrorFlag = 1;
            }
            prevSequenceNum = currSequenceNum;
        }

        // note the index where the ancillary packet can be found in the buffer
        if (currPacketApid == ANCILLARY_APID)
            ancillaryIndex = numPackets;

        // note the indicies where telemetry packets are located in the buffer
        if (currPacketApid >= MIN_APID & currPacketApid <= MAX_APID && currPacketApid != APIDS_WITH_SPIN[0] &&
            currPacketApid != APIDS_WITH_SPIN[1]) {
            telemetryIndices.push_back(numPackets);
        }

        numPackets++;

        if (numPackets == MAX_NUM_PACKETS) {
            cout << "Maximum number of packets: " << MAX_NUM_PACKETS << " exceeded." << endl;
        }

        // Read next packet's length, apid and if its end file flag (NOTE: l0 file pointer does not advance)
        getNextPacketInfo(l0FileStream, packetLen, currPacketApid, isEndFile, isSPW);

        // if next packet is the end, then stop reading packets and return to the main l1agen_oci loop
        if (isEndFile)
            return 0;

        // if next packet is larger than the default packet size, move the packet pointer forward and skip it
        if (packetLen > MAX_PACKET_SIZE) {
            cout << "Packet size > " << MAX_PACKET_SIZE << " (" << packetLen << ")" << endl;
            skipSinglePacket(l0FileStream, packetLen, isEndFile, isSPW);

            // if everything is okay with the next packet, read it
        } else {
            readSinglePacket(l0FileStream, packet, packetLen, currPacketApid, isEndFile, isSPW);
        }

        // current apid is not in valid range, read through the packets until one is in range so it can be
        // added to the buffer list when it goes back to the outter while loop
        while (currPacketApid < MIN_APID || currPacketApid > MAX_APID) {
            getNextPacketInfo(l0FileStream, packetLen, currPacketApid, isEndFile, isSPW);
            if (isEndFile)
                return 0;
            if (packetLen > MAX_PACKET_SIZE) {
                skipSinglePacket(l0FileStream, packetLen, isEndFile, isSPW);
            } else {
                readSinglePacket(l0FileStream, packet, packetLen, currPacketApid, isEndFile, isSPW);
            }
        }

        // Get spin number for the current packet, but separate the num to 2 different variables
        getSpinNumFromPacket(currPacketSpinNum, currPacketSpinNumTelemetry, currPacketApid, packet);

    }  // while (currPacketSpinNum <= firstPacketSpinNum && currPacketSpinNumTelemetry <= firstPacketSpinNum)

    return 0;
}

int compareAncillaryPackets(uint8_t *currAncillaryPacket, uint8_t *nextAncillaryPacket) {
    // 0 if spatial or spectral data tables were changed from the current packet to the next
    int tableNotChanged = 1;

    // ancillary packet without valid mode table, v0.20 == no table change
    if (nextAncillaryPacket[15] > 0)
        return tableNotChanged;

    // ---SPATIAL DATA COLLECTION FIELD---
    const int SPATIAL_DATA_FIELD_OFFSET = 36;
    const size_t NUM_SPATIAL_DATA_FIELDS = 40;
    bool hasSpatialDataFieldDiff = false;  // tracks if the data is different from the next and current packet

    // ---SPECTRAL DATA COLLECTION FIELD---
    const int SPECTRAL_DATA_FIELD_OFFSET = 80;
    const size_t NUM_SPECTRAL_DATA_FIELDS = 12;
    bool hasSpectralDataDiff = false;

    // go through all spatial data field offsets to find a difference
    for (size_t i = 0; i < NUM_SPATIAL_DATA_FIELDS; i++) {
        if (nextAncillaryPacket[SPATIAL_DATA_FIELD_OFFSET + i] !=
            currAncillaryPacket[SPATIAL_DATA_FIELD_OFFSET + i]) {
            hasSpatialDataFieldDiff = true;
            break;
        }
    }

    // go through spectral data fields and find if there is at least 1 difference
    for (size_t i = 0; i < NUM_SPECTRAL_DATA_FIELDS; i++) {
        if (nextAncillaryPacket[SPECTRAL_DATA_FIELD_OFFSET + i] !=
            currAncillaryPacket[SPECTRAL_DATA_FIELD_OFFSET + i]) {
            hasSpectralDataDiff = true;
            break;
        }
    }

    // if either detected a difference with its fields, get the spin number and time and print it
    if (hasSpatialDataFieldDiff or hasSpectralDataDiff) {
        // ---TIME---
        double startTime = 0.0;
        int32_t year = 0;
        int32_t day = 0;
        uint32_t spinNum = 0;

        getAncillaryPacketTime(nextAncillaryPacket, year, day, startTime);
        memcpy(&spinNum, &nextAncillaryPacket[24], 4);
        swapc_bytes((char *)&spinNum, sizeof(uint32_t), 1);

        uint16_t hour = (uint16_t)floor(startTime / 3600);
        uint16_t mins = (uint16_t)floor((startTime - hour * 3600) / 60);
        uint16_t sec = (uint16_t)floor(startTime - hour * 3600 - mins * 60);

        // report based on which has has a diff
        if (hasSpatialDataFieldDiff) {
            cout << "Spatial table change at: spin=" << spinNum << ", " << hour << ":" << mins << ":" << sec
                 << endl;
        }

        if (hasSpectralDataDiff) {
            cout << "Spectral table change at: spin=" << spinNum << ", " << hour << ":" << mins << ":" << sec
                 << endl;
        }

        tableNotChanged = 0;
    }

    return tableNotChanged;
}

// Skips a packet from the l0 file stream so you do not need to read into a buffer and then discard it
int skipSinglePacket(L0Stream *l0FileStream, uint32_t &packetLen, int32_t &isEndFile, bool isSPW) {
    // no more packets to read
    if (l0FileStream->eof()) {
        isEndFile = 1;
        packetLen = 0;
        cout << "End of packet file" << endl;
        return 0;
    }

    // Read spacewire header if needed
    if (isSPW) {
        uint8_t spwhead[2];
        l0FileStream->read((char *)&spwhead, 2);
    }

    l0FileStream->skip(packetLen);
    return 0;
}
// Loads a single packet into the buffer and extract the apid, packet length and note if its the end of file
// if packetBuffer is null, don't load anything and just read in the packet length, apid and end file flag
int readSinglePacket(L0Stream *l0FileStream, uint8_t *packetBuffer, uint32_t &packetLen, uint32_t &apid,
                     int32_t &isEndFile, bool isSPW) {
    // no more packets to read
    if (l0FileStream->eof()) {
        isEndFile = 1;
        packetLen = 0;
        cout << "End of packet file" << endl;
        return 0;
    }

    // Read spacewire header if needed
    if (isSPW) {
        uint8_t spwhead[2];
        l0FileStream->read((char *)&spwhead, 2);
    }

    // read into buffer
    l0FileStream->read((char *)packetBuffer, packetLen);

    // extract the apid from the buffer and save it
    apid = (packetBuffer[0] % 8) * 256 + packetBuffer[1];

    return 0;
}

// Given a pointer to multiple packets, find a Swir packet and get the mode from its meta data
int getSwirMode(uint8_t (*packets)[MAX_PACKET_SIZE], uint32_t numPackets, uint16_t &swirMode) {
    // go through each packet and find the SWIR band packet with apid 720
    for (size_t i = 0; i < numPackets; i++) {
        uint32_t apid = (packets[i][0] % 8) * 256 + packets[i][1];

        if (apid == 720) {
            //  Get mode from packet
            uint8_t swirMetaData = packets[i][12];
            swirMode = (swirMetaData % 64) / 16;
            return 0;
        }
    }
    return 0;
}

// Get the swir band data from the packet buffer and save it to the outDataArrRef
int extractSwirBandData(uint8_t *packetBuffer, uint32_t *outDataArrRef) {
    // This routine takes groups of 9 20-bit samples that are packed into
    // 23 bytes and unpacks them into a 4-byte integer array

    for (size_t i = 0; i < 5; i++) {
        outDataArrRef[i * 2] =
            4096 * packetBuffer[i * 5] + 16 * packetBuffer[i * 5 + 1] + packetBuffer[i * 5 + 2] / 16;
    }

    for (size_t i = 0; i < 4; i++) {
        outDataArrRef[i * 2 + 1] =
            65536 * (packetBuffer[i * 5 + 2] % 16) + 256 * packetBuffer[i * 5 + 3] + packetBuffer[i * 5 + 4];
    }

    return 0;
}

// given a packet reference, extract swir line numbers, pixel frame type and data
int unpackSwirPacket(uint8_t *packet, int16_t *swirLineNums, uint8_t *swirFrameTypeArr, uint32_t *swirData) {
    // Program to unpack one OCI SWIR data packet
    // Each packet contains data for eight science pixels
    // Reference: OCI-ELEC-SPEC-0028

    // Allocate output buffers (9 bands for 8 pixels)

    int offset = 10;
    uint8_t swirMetaData;

    for (size_t i = 0; i < 8; i++) {
        // Get line number and metadata
        if ((unsigned(packet[offset]) == 255) && (unsigned(packet[offset + 1]) == 255)) {
            // handle fill values in SWIR, LH 10/28/2020
            swirLineNums[i] = 0;
        } else {
            swirLineNums[i] = packet[offset] * 256 + packet[offset + 1];
        }

        swirMetaData = packet[offset + 2];
        swirFrameTypeArr[i] = swirMetaData % 8;

        // Extract SWIR band data
        extractSwirBandData(&packet[offset + 3], &swirData[9 * i]);

        offset += 26;
    }

    return 0;
}

// Used when unpacking science data. Passes in an ancillary packet and depending on the CCD Id, use
// either the red or blue taps to determine which taps has valid data to extract band data from
int unpackCcdPacket(uint8_t *ancillaryPacket, uint16_t isBlueCcdTapsEnabled[NUM_TAPS],
                    uint16_t isRedCcdTapsEnabled[NUM_TAPS], uint16_t &ccdId, uint32_t &lineNum,
                    uint16_t &dataType, uint16_t &spatialAggregation, uint16_t spectralAggregation[NUM_TAPS],
                    uint16_t &numBands, uint16_t **ccdData, uint16_t overScanSumData[NUM_TAPS]) {
    // Get CCD ID, line number, data type and spatial aggregation and save it to the referenced variable
    ccdId = (ancillaryPacket[12] & 64) / 64;
    lineNum = ancillaryPacket[10] * 256 + ancillaryPacket[11];
    dataType = (ancillaryPacket[12] % 64) / 4;
    uint16_t overScanSum = ancillaryPacket[17] % 16;
    spatialAggregation = ancillaryPacket[12] % 4;
    uint16_t spectralAggregationCodes[4] = {1, 2, 4, 8};
    spatialAggregation = spectralAggregationCodes[spatialAggregation];

    // Get tap enable flags for the current focal plane. If the id is 1, use blue and if it is 0,
    // use red
    uint16_t *currCcdTapEnabledFlags;
    if (ccdId)
        currCcdTapEnabledFlags = isBlueCcdTapsEnabled;
    else
        currCcdTapEnabledFlags = isRedCcdTapsEnabled;

    // Get spectral aggregation factors for all taps. The indicies indicate the range to write the data into
    // So:
    //      spectraAggIndices[0] will access spectralAggregation array from index 0 to 3.
    //      spectraAggIndices[2] ... index 4 to 7
    //      ....
    //      spectraAggIndices[n] ... index n to n+3
    uint16_t spectraAggIndices[4] = {0, 4, 8, 12};
    uint32_t taps[16];
    for (size_t i = 0; i < 4; i++) {
        spectralAggregation[spectraAggIndices[i]] =
            spectralAggregationCodes[(ancillaryPacket[13 + i] & 192) / 64];
        spectralAggregation[spectraAggIndices[i] + 1] =
            spectralAggregationCodes[(ancillaryPacket[13 + i] & 48) / 16];
        spectralAggregation[spectraAggIndices[i] + 2] =
            spectralAggregationCodes[(ancillaryPacket[13 + i] & 12) / 4];
        spectralAggregation[spectraAggIndices[i] + 3] =
            spectralAggregationCodes[(ancillaryPacket[13 + i] & 3)];
    }

    // set red or blue tap flags
    for (size_t i = 0; i < 16; i++)
        taps[i] = 32 * currCcdTapEnabledFlags[i] / spectralAggregation[i];

    // Allocate output buffer
    numBands = 0;
    for (size_t i = 0; i < 16; i++)
        numBands += taps[i];

    *ccdData = new uint16_t[numBands];

    // offset from the start of the ancillary packet to get band data
    int bandDataPacketOffset = 18;
    // the current band that is being extracted from the packet
    uint16_t currExtractingBandNum = numBands - 1;
    // For each tap (1-16):
    for (size_t j = 0; j < NUM_TAPS; j++) {
        // Copy band data from ancillaryPacket to output buffer
        // Packet data are in reverse spectral order, so swap order here
        if (currCcdTapEnabledFlags[j]) {
            for (size_t i = 0; i < taps[j]; i++) {
                uint16_t temp = 0;
                // save data to temp, swap the order then save it to ccdData array
                memcpy(&temp, &ancillaryPacket[bandDataPacketOffset + 2 * i], 2);
                swapc_bytes((char *)&temp, sizeof(int16_t), 1);
                memcpy(&(*ccdData)[currExtractingBandNum - i], &temp, 2);
            }

            // update the band number to get next and the offset
            currExtractingBandNum -= taps[j];
            bandDataPacketOffset += 2 * taps[j];
        }
    }
    // over scanned, extract the over scanned data
    if (overScanSum > 0) {
        for (size_t j = 0; j < NUM_TAPS; j++) {
            uint16_t temp = 0;
            memcpy(&temp, &ancillaryPacket[bandDataPacketOffset + 2 * j], 2);
            swapc_bytes((char *)&temp, sizeof(int16_t), 1);
            memcpy(&overScanSumData[j], &temp, 2);
        }
    }
    return 0;
}

/** Compare the spectral aggregation codes from when it was extracted from the packet to when it was copied
 * @brief
 * @param spectralAggCodes - from packet aggregation code array of size 16
 * @param bandSpectralAggCodes - copied from packet to its own array
 * @return true if matching, false otherwise
 */
bool inconsistentAggregationCodes(uint16_t *spectralAggCodes, uint16_t *bandSpectralAggCodes) {
    uint8_t diffCount = 0;
    for (size_t i = 0; i < 16; i++) {
        if (spectralAggCodes[i] != bandSpectralAggCodes[i]) {
            diffCount++;
        }
    }
    return diffCount > 0 ? true : false;
}

/** Compare datatypes and make sure they are the same
 * @brief
 * @param dataType
 * @param bandDataType
 * @return
 */
bool inconsistentDataTypes(uint16_t &dataType, uint16_t &bandDataType) {
    return (dataType != bandDataType && dataType != 2 && dataType != 5);
}

int unpackScienceData(uint32_t totalNumPackets, int32_t spin, uint16_t numCcdPixels, uint16_t numSwirPixels,
                      uint16_t maxNumOfSwirPixels, uint16_t &numBandsRef,
                      uint16_t isBlueCcdTapsEnabled[NUM_TAPS], uint16_t isRedCcdTapsEnabled[NUM_TAPS],
                      uint8_t (*ociPacketBuffer)[MAX_PACKET_SIZE], uint16_t **l0BlueBandData,
                      uint16_t **l0RedBandData, uint32_t **l0SwirBandData, int16_t *blueCcdLineNums,
                      int16_t *redCcdLineNums, int16_t *swirLineNums, uint16_t &blueDataType,
                      uint16_t blueSpectralAggregations[NUM_TAPS], uint16_t &redDataType,
                      uint16_t redSpectralAggregations[NUM_TAPS], int8_t *swirFrameTypeArrRef,
                      int &returnStatus) {
    // Unpack all of the science data for an OCI scan.
    // The SWIR bands are re-ordered in ascending wavelength order.
    // Order of dual-gain bands is (SG, HG).

    for (size_t i = 0; i < numCcdPixels; i++) {
        blueCcdLineNums[i] = -1;
        redCcdLineNums[i] = -1;
    }
    for (size_t i = 0; i < maxNumOfSwirPixels; i++)
        swirLineNums[i] = -1;

    // lists indices for bands. this is the default band order
    std::vector<int> bandArrIndices = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    // SWIR band sorting indices for science mode
    uint16_t swirMode = 0;
    getSwirMode((uint8_t(*)[MAX_PACKET_SIZE]) & ociPacketBuffer[0][0], totalNumPackets, swirMode);

    if (swirMode == 0) {
        // indicies to access the band in acending wavelength order
        bandArrIndices = {3, 0, 1, 2, 8, 6, 7, 5, 4};
    }

    uint16_t ccdId;        // store the current packets ccd id
    uint32_t currLineNum;  // curr packets line number to report error if detected
    uint16_t spatialAggCode;
    uint16_t spectralAggCodePerTap[NUM_TAPS];

    uint16_t **ccdData = new uint16_t *;
    uint16_t overScanSumData[NUM_TAPS];

    uint16_t numBlueBands;
    uint16_t numRedBands;
    uint16_t numBands;  // num bands for this local function and not a ref from l1agen_oci main

    uint32_t apid;
    int32_t spinNum;

    int bluePixelCount = 0;
    int redPixelCount = 0;
    int swirPixelCount = 0;
    returnStatus = 0;

    for (size_t ipkt = 0; ipkt < totalNumPackets; ipkt++) {
        uint32_t uint32Temp = 0;  // temp to store uint32 for byte swapping
        apid = (ociPacketBuffer[ipkt][0] % 8) * 256 + ociPacketBuffer[ipkt][1];
        uint16_t dataType = (ociPacketBuffer[ipkt][12] % 64) / 4;
        memcpy(&uint32Temp, &ociPacketBuffer[ipkt][6], 4);
        swapc_bytes((char *)&uint32Temp, sizeof(uint32_t), 1);
        spinNum = uint32Temp;

        // If CCD science APID
        if (apid == 700 && dataType > 0 && spinNum == spin) {
            unpackCcdPacket(ociPacketBuffer[ipkt], isBlueCcdTapsEnabled, isRedCcdTapsEnabled, ccdId,
                            currLineNum, dataType, spatialAggCode, spectralAggCodePerTap, numBandsRef,
                            ccdData, overScanSumData);

            // If blue
            if (ccdId) {
                if (bluePixelCount == 0) {
                    numBlueBands = numBandsRef;
                    blueDataType = dataType;
                    memcpy(blueSpectralAggregations, spectralAggCodePerTap, 16 * sizeof(uint16_t));
                    for (size_t i = 0; i < numCcdPixels; i++) {
                        if (l0BlueBandData[i])
                            delete[] l0BlueBandData[i];
                        l0BlueBandData[i] = new uint16_t[numBandsRef];
                        for (size_t j = 0; j < numBandsRef; j++)
                            l0BlueBandData[i][j] = 65535;  // LH, 11/23/2020
                    }

                } else {
                    // Check for inconsistent non-dark data type or spectral aggregation
                    bool hasInconsistentSpectralAgg =
                        inconsistentAggregationCodes(spectralAggCodePerTap, blueSpectralAggregations);
                    bool hasInconsistentDataType = inconsistentDataTypes(dataType, redDataType);

                    if (hasInconsistentDataType || hasInconsistentSpectralAgg) {  // Ignore dark and linearity
                        cout << "Data type or spectral aggregation error, CCDID: " << ccdId
                             << "  Line: " << currLineNum << endl;
                        returnStatus = 4;
                    }
                }
                if (numBandsRef <= numBlueBands)
                    numBands = numBandsRef;
                else
                    numBands = numBlueBands;
                if (bluePixelCount < numCcdPixels) {
                    memcpy(l0BlueBandData[bluePixelCount], &(*ccdData)[0], numBands * sizeof(uint16_t));
                    blueCcdLineNums[bluePixelCount] = (currLineNum / spatialAggCode) * spatialAggCode;
                } else {
                    cout << "Number of blue pixels exceeded in spin: " << spin
                         << " in packet (0-based): " << ipkt << endl;
                }
                bluePixelCount++;

                // Red
            } else {
                if (redPixelCount == 0) {
                    numRedBands = numBandsRef;
                    redDataType = dataType;
                    memcpy(redSpectralAggregations, spectralAggCodePerTap, 16 * sizeof(uint16_t));
                    for (size_t i = 0; i < numCcdPixels; i++) {
                        if (l0RedBandData[i])
                            delete[] l0RedBandData[i];
                        l0RedBandData[i] = new uint16_t[numBandsRef];
                        for (size_t j = 0; j < numBandsRef; j++)
                            l0RedBandData[i][j] = 65535;  // LH, 11/23/2020
                    }
                } else {
                    // Check for inconsistent non-dark data type or spectral aggregation
                    bool hasInconsistentSpectralAgg =
                        inconsistentAggregationCodes(spectralAggCodePerTap, redSpectralAggregations);
                    bool hasInconsistentDataType = inconsistentDataTypes(dataType, redDataType);

                    if (hasInconsistentDataType || hasInconsistentSpectralAgg) {  // Ignore dark and linearity
                        cout << "Data type or spectral aggregation error, CCDID: " << ccdId
                             << "  Line: " << currLineNum << endl;
                        returnStatus = 4;
                    }
                }

                if (numBandsRef <= numRedBands)
                    numBands = numBandsRef;
                else
                    numBands = numRedBands;

                if (redPixelCount < numCcdPixels) {
                    memcpy(l0RedBandData[redPixelCount], &(*ccdData)[0], numBands * sizeof(uint16_t));
                    redCcdLineNums[redPixelCount] = (currLineNum / spatialAggCode) * spatialAggCode;
                } else {
                    cout << "Number of red pixels exceeded in spin: " << spin
                         << " in packet (0-based): " << ipkt << endl;
                }
                redPixelCount++;
            }  // if (ccdId)

            delete[] *ccdData;

        }  // if (apid == 700 && dataType > 0)

        // If SWIR science APID
        if (apid == 720 && (swirPixelCount + 7) < maxNumOfSwirPixels && spinNum == spin) {
            int16_t swirLineNumArr[8];
            uint8_t swirFrameTypeArr[8];
            uint32_t swirData[8 * 9];

            // swirData: 8 rows of 9 columns
            unpackSwirPacket(ociPacketBuffer[ipkt], swirLineNumArr, swirFrameTypeArr, swirData);

            for (size_t i = 0; i < 8; i++) {
                swirLineNums[swirPixelCount + i] = (swirLineNumArr[i] / 8) * 8;
                for (size_t j = 0; j < 9; j++)
                    l0SwirBandData[swirPixelCount + i][j] =
                        swirData[bandArrIndices[j] + 9 * i];  // bands in ascending order
            }

            // copy to the the reference array from l1agen_oci
            memcpy(&swirFrameTypeArrRef[swirPixelCount], swirFrameTypeArr, 8 * sizeof(uint8_t));

            swirPixelCount += 8;
        }  // if (apid == 720 && (swirPixelCount+7) < maxNumOfSwirPixels)

        if ((bluePixelCount <= 0) || (bluePixelCount > numCcdPixels))
            blueCcdLineNums[0] = -1;
        if ((redPixelCount <= 0) || (redPixelCount > numCcdPixels))
            redCcdLineNums[0] = -1;
        if ((swirPixelCount <= 0) || (swirPixelCount > maxNumOfSwirPixels))
            swirLineNums[0] = -1;

    }  // ipkt loop

    delete ccdData;

    return 0;
}

// generate oci line index array for ccd and swir bands
int makeOciLineIndex(spatialAggTable *spatialAggList, int16_t *ccdBandLineIndices,
                     int16_t *swirBandLineIndices, int16_t *ccdBandDarkLineIndices,
                     int16_t *swirBandDarkLineIndices, int16_t *swirLineOffsetArr) {
    const uint16_t NUM_SWIR_BANDS = 9;
    const size_t TABLE_SIZE = 10;
    const size_t AGGRIGATION_SIZE = 8;

    // Loop through data zones
    uint16_t lineOffset = 0;
    uint16_t ccdPixel = 0;
    uint16_t swirPixel = 0;
    uint16_t ccdDarkPixel = 0;
    uint16_t swirDarkPixel = 0;

    for (size_t i = 0; i < TABLE_SIZE; i++) {
        // If not "no data" type, add indices to array
        if ((spatialAggList[i].dataType != 0) && (spatialAggList[i].dataType != TABLE_SIZE)) {
            if (spatialAggList[i].lines > 0) {
                // CCD pixel index
                if (spatialAggList[i].spatialAggCode > 3 || spatialAggList[i].spatialAggCode < 0) {
                    printf(
                        "--Error-- : the value of spatialAggList at i = %d is an out of range index %d. See "
                        "%s at "
                        "%d\n\n",
                        int(i), spatialAggList[i].spatialAggCode, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                uint16_t spatialAggregation = SPATIAL_AGG_CODES[spatialAggList[i].spatialAggCode];
                uint16_t lines = spatialAggList[i].lines;
                // Check for total lines within limit
                if ((lineOffset + lines) > MAX_LINES) {
                    cout << "Mode table entry " << i << " exceeds max lines" << endl;
                    lines = MAX_LINES - lineOffset - spatialAggregation;
                }

                uint16_t numCcdLinesAfterAgg = lines / spatialAggregation;
                for (size_t j = 0; j < numCcdLinesAfterAgg; j++) {
                    uint16_t ccdIndex = lineOffset + j * spatialAggregation;
                    ccdBandLineIndices[ccdIndex] = ccdPixel + j;
                }
                ccdPixel += numCcdLinesAfterAgg;

                // SWIR pixel index (fixed aggregation of AGGRIGATION_SIZE)
                uint16_t numSwirLinesAfterAgg = lines / AGGRIGATION_SIZE;
                for (size_t j = 0; j < numSwirLinesAfterAgg; j++) {
                    uint16_t swirIndex = (lineOffset / AGGRIGATION_SIZE + j) * AGGRIGATION_SIZE;
                    for (size_t k = 0; k < NUM_SWIR_BANDS; k++) {
                        int16_t swirLineOffset = 0;
                        // If not dark view, use line offsets
                        if (spatialAggList[i].dataType != 2) {
                            swirLineOffset = swirLineOffsetArr[k];
                        }
                        if ((swirIndex - swirLineOffset) < 0) {
                            printf("--Error-- : %s:%d : SWIR line offset caused a negative index\n", __FILE__,
                                   __LINE__);
                            exit(EXIT_FAILURE);
                        }
                        swirBandLineIndices[(swirIndex - swirLineOffset) * NUM_SWIR_BANDS + k] =
                            swirPixel + j;
                    }
                }
                swirPixel += numSwirLinesAfterAgg;

                // Dark collect ( if dataType equals 2)
                if (spatialAggList[i].dataType == 2) {
                    for (size_t j = 0; j < numCcdLinesAfterAgg; j++) {
                        uint16_t ccdIndex = lineOffset + j * spatialAggregation;
                        ccdBandDarkLineIndices[ccdIndex] = ccdDarkPixel + j;
                    }
                    ccdDarkPixel += numCcdLinesAfterAgg;

                    for (size_t j = 0; j < numSwirLinesAfterAgg; j++) {
                        uint16_t swirIndex = lineOffset + j * AGGRIGATION_SIZE;
                        swirBandDarkLineIndices[swirIndex] = swirDarkPixel + j;
                    }
                    swirDarkPixel += numSwirLinesAfterAgg;
                }
            } else {
                cout << "Data zone " << i << " type " << spatialAggList[i].dataType << " has zero lines"
                     << endl;
            }
        }
        lineOffset += spatialAggList[i].lines;
        if (lineOffset >= MAX_LINES)
            break;
    }

    return 0;
}
