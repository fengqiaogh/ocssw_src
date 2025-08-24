#include <stdio.h>
#include <math.h>
#include <string.h>
#include <regex>
#include <getopt.h>
#include <clo.h>
#include <boost/algorithm/string.hpp>
#include <libgen.h>

#include "l0stream.hpp"
#include "l1afile_oci.hpp"
#include "process_science_data.h"
#include "l1afile_oci_manager.hpp"

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


/**
 * @brief Reference to dark calibration **uint16_t** array of size maxScans for blue, red or swir
 * @param arr blue, red or swir array that needs to be initalized
 * @param numPixels dark ccd or normal ccd pixels
 * @param numBands  blue, red or swir pixels
 * @param maxScans number of scans for this instance
 */
void initializeDarkCalibrationArray(uint16_t **&arr, uint16_t numPixels, uint16_t numBands,
                                    uint16_t maxScans) {
    // initialize the first index of the array to hold (numBands * numPixels) of elements for all scans
    // if numbands is 0, use 1.
    if (numBands > 0) {
        arr[0] = new uint16_t[numPixels * numBands * maxScans];
    } else {
        arr[0] = new uint16_t[numPixels * maxScans];
    }

    //  define the start index for each scan relative to arr[0]. So:
    //      arr[1]'s start index should be: arr[0] + (numPixels * numBands)
    //      arr[2]'s  start index should be: arr[1] + (numPixels * numBands)
    for (size_t i = 1; i < maxScans; i++) {
        arr[i] = arr[i - 1] + numPixels * numBands;
    }
}

/**
 * @brief Reference to dark calibration **uint32_t** array of size maxScans for blue, red or swir
 * @param arr overloaded uint32_t version instead of uint16_t (short)
 * @param numPixels dark ccd or normal ccd pixels
 * @param numBands  blue, red or swir pixels
 * @param maxScans number of scans for this instance
 */
void initializeDarkCalibrationArray(uint32_t **&arr, uint16_t numPixels, uint16_t numBands,
                                    uint16_t maxScans) {
    // initialize the first index of the array to hold (numBands * numPixels) of elements for all scans
    // if numbands is 0, use 1.
    if (numBands > 0) {
        arr[0] = new uint32_t[numPixels * numBands * maxScans];
    } else {
        arr[0] = new uint32_t[numPixels * maxScans];
    }

    //  define the start index for each scan relative to arr[0]. So:
    //      arr[1]'s start index should be: arr[0] + (numPixels * numBands)
    //      arr[2]'s  start index should be: arr[1] + (numPixels * numBands)
    for (size_t i = 1; i < maxScans; i++) {
        arr[i] = arr[i - 1] + numPixels * numBands;
    }
}

/**
 * @brief Reference to science **uint16_t** array of size numBands
 * @param arr array that needs to be initalized
 * @param numPixels dark ccd or normal ccd pixels
 * @param numBands  blue, red or swir pixels
 */
void initializeScienceDataArray(uint16_t **&arr, uint16_t numPixels, uint16_t numBands) {
    // initialize the first index of the array to hold (numBands * numPixels) of elements for all scans
    // if numbands is 0, use 1.
    if (numBands > 0) {
        arr[0] = new uint16_t[numPixels * numBands];
    } else {
        arr[0] = new uint16_t[numPixels];
    }

    //  define the start index for each scan relative to arr[0]. So:
    //      arr[1]'s start index should be: arr[0] + numPixels
    //      arr[2]'s start index should be: arr[1] + numPixels
    for (size_t i = 1; i < numBands; i++) {
        arr[i] = arr[i - 1] + numPixels;
    }
}

/**
 * @brief Reference to science **uint32_t** array of size numBands
 * @param arr array that needs to be initalized
 * @param numPixels dark ccd or normal ccd pixels
 * @param numBands  blue, red or swir pixels
 */
void initializeScienceDataArray(uint32_t **&arr, uint16_t numPixels, uint16_t numBands) {
    // initialize the first index of the array to hold (numBands * numPixels) of elements for all scans
    // if numbands is 0, use 1.
    if (numBands > 0) {
        arr[0] = new uint32_t[numPixels * numBands];
    } else {
        arr[0] = new uint32_t[numPixels];
    }

    //  define the start index for each scan relative to arr[0]. So:
    //      arr[1]'s start index should be: arr[0] + numPixels
    //      arr[2]'s start index should be: arr[1] + numPixels
    for (size_t i = 1; i < numBands; i++) {
        arr[i] = arr[i - 1] + numPixels;
    }
}

void verify_packet_len(L0Stream *tfileStream, uint32_t &currPacketLength, uint32_t apid, int32_t isEndFile,
                       bool isSPW) {
    getNextPacketInfo(tfileStream, currPacketLength, apid, isEndFile, isSPW);
    if (currPacketLength > MAX_PACKET_SIZE) {
        cout << "Packet too big (" << currPacketLength << ") for buffer (" << MAX_PACKET_SIZE << ")" << endl;
        exit(EXIT_FAILURE);
    }
}

//! Makes clo option and alias at the same time. Wraps around the clo_addOption and clo_addAlias.
void makeCloOptionAndAlias(clo_optionList_t *list, const char *optionName, const char *alias,
                           enum clo_dataType_t dataType, const char *defaultVal, const char *desc) {
    // make option
    clo_addOption(list, optionName, dataType, defaultVal, desc);

    // make alias
    clo_addAlias(list, optionName, alias);
}

// take time make hours, mins and seconds. Return a time formatted string with Z
stringstream makeTimeCoverageString(AncillaryPktTimeStamp &time) {
    // Write start, end, create time attributes
    int16_t month = 0;
    int16_t day = 0;
    int32_t hour = 0;
    int32_t mins = 0;
    stringstream timeCoveragetStr;  // format time coverage start and end using this stream then write it
    
    yd2md((int16_t)time.year, (int16_t)time.day, &month, &day);

    time.second = floor(time.second * 1000) / 1000;  // LH, 11/18/2020
    hour = (int)(time.second / 3600);
    time.second -= hour * 3600;
    mins = (int)(time.second / 60);
    time.second -= mins * 60;
    // yyyy-mn-dyThr:mn:ss.sss
    // for each, set the width of how long it should be and set a fill value if it is less than the width set.
    timeCoveragetStr = stringstream();
    timeCoveragetStr << setw(4) << to_string(time.year) << "-";  // year is 4 digits always, no fill
    timeCoveragetStr << setw(2) << setfill('0') << month << "-";      // month can be single digit, fill w/ 0
    timeCoveragetStr << setw(2) << setfill('0') << day << "T";        // day is the same as month

    timeCoveragetStr << setw(2) << setfill('0') << hour << ":";
    timeCoveragetStr << setw(2) << setfill('0') << mins << ":";
    timeCoveragetStr << fixed << setw(6) << setprecision(3) << setfill('0') << time.second;

    return timeCoveragetStr;
}

int main(int argc, char *argv[]) {
    // version of l1agen_oci
    const string VERSION = "02.13.00_2025-08-20";

    cout << "l1agen_oci " << VERSION << " (" << __DATE__ << " " << __TIME__ << ")" << endl;

    // Set up buffer to read in command line options
    clo_optionList_t *commandLineOptionsList = clo_createList();
    clo_setHelpStr(
        "\nUsage: l1agen_oci [options] \n\nRetuns status codes:\n\
        0 - good\n\
        110 - no ancillary packets found in file\n\
        120 - no L1A file was generated\n\
        130 - duplicate l1a file being generated\n"
    );

    // L0 file list and granule length, no alias.
    clo_addOption(commandLineOptionsList, "ifile", CLO_TYPE_IFILE, NULL,
                  "Single .oci file or a text file that lists many .oci files");
    makeCloOptionAndAlias(commandLineOptionsList, "ofile", "f", CLO_TYPE_STRING, NULL,
                          "Name of L1A output file if not following system naming.");
    clo_addOption(commandLineOptionsList, "granule_len", CLO_TYPE_STRING, NULL, "Length of the L1A file");

    makeCloOptionAndAlias(commandLineOptionsList, "start_time", "t", CLO_TYPE_STRING, "",
                          "Granule Start Time. Format: YYYYmmddTHHMMSS or YYYY-mm-ddTHH:MM:SS");
    makeCloOptionAndAlias(commandLineOptionsList, "maxgap", "g", CLO_TYPE_INT, "65535",
                          "Max Gap between Spins");
    makeCloOptionAndAlias(commandLineOptionsList, "max_file_gap", "m", CLO_TYPE_DOUBLE, "300",
                        "Max Gap between Files in Seconds");
    makeCloOptionAndAlias(commandLineOptionsList, "hktlist", "k", CLO_TYPE_IFILE, NULL,
                          "List of Housekeeping Telemetry (HKT) files");
    makeCloOptionAndAlias(commandLineOptionsList, "swir_loff_set", "s", CLO_TYPE_STRING, NULL,
                          "SWIR Line Offset Config. list of 9 comma separated integers.");
    makeCloOptionAndAlias(commandLineOptionsList, "outlist", "o", CLO_TYPE_STRING, NULL,
                          "Lists output files generated by l1agen_oci with timestamps.");
    makeCloOptionAndAlias(commandLineOptionsList, "nametag", "p", CLO_TYPE_STRING, "PACE_OCI",
                          "Prefix to output file names");
    makeCloOptionAndAlias(commandLineOptionsList, "doi", "d", CLO_TYPE_STRING, NULL,
                          "Digital Object Identifier (DOI)");
    makeCloOptionAndAlias(commandLineOptionsList, "pversion", "v", CLO_TYPE_STRING, "Unspecified",
                          "Program Processing Version");
    makeCloOptionAndAlias(commandLineOptionsList, "noSPW", "n", CLO_TYPE_BOOL, "false",
                          "No Spacewire header Flag");

    // if no arguments in the command line, print the useage to the user and exit.
    if (argc == 1) {
        clo_printHelpString();
        clo_printOptions(commandLineOptionsList);
        exit(EXIT_FAILURE);
    }

    // record down the command line parameters used to run this instance of the program. This will be
    // saved in the global attributes "history"
    string history = call_sequence(argc, argv);

    // flag to determine if you need to read the advanace the packet pointer until it finds a packet where
    // the start time is >= granuleStartTime
    bool isGranuleStartTimeGiven = false;
    int maxGap = 10;            // by default, the max allowed missing scan  = 10
    double maxFileGap = 300;    // default == 300 seconds == 5 mins 
    bool isSPW = true;          // by default, OCI data is from DSB and has space wire header
    int returnStatus = 0;

    // command line input parameters
    string ifileName = "";            // name of text file that will list out all L0 files being used
    string granuleLengthStr = "";  // how long the granule should be in mins
    string hktList = "";           // name of text file that will list out all the hkt files being used
    string swirLineOffsetString = "";
    string granuleStartTimeStr = "";
    string outlist = "";  // name of file that will list out all files made with timestamps
    string ofileName = ""; // Optional
    string outfilePrefixTag = "";     // prefix to output files. if empty, default to "PACE_OCI"
    string doi = "";                  // digital object identifier
    string pversion = "Unspecified";  // processing version

    // parse all the command line arguments and read them
    clo_readArgs(commandLineOptionsList, argc, argv);

    // list of L0 files to open
    if (clo_isSet(commandLineOptionsList, "ifile")) {
        ifileName = clo_getString(commandLineOptionsList, "ifile");
    } else {
        invalid_argument up = invalid_argument("Input L0 file must be specified");
        throw up;
    }

    // granule length
    if (clo_isSet(commandLineOptionsList, "granule_len")) {
        granuleLengthStr = clo_getString(commandLineOptionsList, "granule_len");
    }

    // maxgap
    if (clo_isSet(commandLineOptionsList, "maxgap")) {
        maxGap = clo_getInt(commandLineOptionsList, "maxgap");
        if (maxGap == 0) {
            maxGap = 65535;
        }
    }

    // max_file_gap
    if (clo_isSet(commandLineOptionsList, "max_file_gap")) {
        maxFileGap = clo_getDouble(commandLineOptionsList, "max_file_gap");
    }

    // swir line offset
    if (clo_isSet(commandLineOptionsList, "swir_loff_set")) {
        swirLineOffsetString = clo_getString(commandLineOptionsList, "swir_loff_set");
    }

    // start_time
    if (clo_isSet(commandLineOptionsList, "start_time")) {
        granuleStartTimeStr = clo_getString(commandLineOptionsList, "start_time");
        isGranuleStartTimeGiven = true;
    }

    // get HKT List
    if (clo_isSet(commandLineOptionsList, "hktlist")) {
        hktList = clo_getString(commandLineOptionsList, "hktlist");
    } else {
        invalid_argument up = invalid_argument("List of HKT files must be specified");
        throw up;
    }

    // file to write outputs generated by l1agen_oci
    if (clo_isSet(commandLineOptionsList, "outlist")) {
        outlist = clo_getString(commandLineOptionsList, "outlist");
    }

    // optional output file name
    if (clo_isSet(commandLineOptionsList, "ofile")) {
        ofileName = clo_getString(commandLineOptionsList, "ofile");
    }

    // diff nametag that is not PACE_OCI
    if (clo_isSet(commandLineOptionsList, "nametag")) {
        outfilePrefixTag = clo_getString(commandLineOptionsList, "nametag");
    }

    // doi
    if (clo_isSet(commandLineOptionsList, "doi")) {
        doi = clo_getString(commandLineOptionsList, "doi");
    }

    // pversion
    if (clo_isSet(commandLineOptionsList, "pversion")) {
        pversion = clo_getString(commandLineOptionsList, "pversion");
    }

    // has spacewire header flag
    if (clo_isSet(commandLineOptionsList, "noSPW")) {
        isSPW = true;
    }

    // free list sincei t wont be used anymore from here on out
    free(commandLineOptionsList);

    nc_set_chunk_cache(CHUNK_CACHE_SIZE, CHUNK_CACHE_NELEMS, CHUNK_CACHE_PREEMPTION);

    // Error if output is not just one granule and an output file name is specified
    if (granuleLengthStr != "0" && granuleStartTimeStr.empty() && !ofileName.empty()) {
        cout << "Set granule length parameter to 0 or provide start time when specifying output file name."
             << endl;
        exit(EXIT_FAILURE);
    }

    // Break the given start time and date into individual int componets. Required tm struct to pass into
    // the strptime function
    struct tm granuleStartTimeIntStruct;

    if (isGranuleStartTimeGiven) {
        // regular expression that highlights date time pattern YYYY-mm-ddTHH:MM:SS

        const std::regex timeHasPunctuationRegEx(
            "([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2})");

        // regular expression that highlights date time pattern YYYYmmddTHHMMSS
        const std::regex timeNoPunctuationRegEx(
            "([0-9]{4})([0-9]{2})([0-9]{2})T([0-9]{2})([0-9]{2})([0-9]{2})");

        // take the input start time string and parse it into an int time structure
        if (regex_match(granuleStartTimeStr, timeHasPunctuationRegEx)) {
            strptime(granuleStartTimeStr.c_str(), "%Y-%m-%dT%H:%M:%S", &granuleStartTimeIntStruct);

        } else if (regex_match(granuleStartTimeStr, timeNoPunctuationRegEx)) {
            strptime(granuleStartTimeStr.c_str(), "%Y%m%dT%H%M%S", &granuleStartTimeIntStruct);
        } else {
            cout << "Start time format is wrong. Use  YYYY-mm-ddTHH:MM:SS  or YYYYmmddTHHMMSS." << endl;
            exit(EXIT_FAILURE);
        }
    }

    vector<string> fileNames;
    fileNames = readFileList(ifileName);
    L0Stream tfileStream(fileNames);

    // if the input is empty
    if (fileNames.size() == 0) {
        cout << "Error - No L0 files found" << endl;
        return 0;
    }
    if (tfileStream.fail())
        cout << "Failed to open tfileStream at line " << __LINE__ << endl;

    uint32_t apid = 0;
    uint32_t currPacketLength = 0;  // length of packet
    int32_t isEndFile = 0;
    uint8_t packetBuffer[MAX_PACKET_SIZE];

    // when reading and writing Ancillary packets, you get some from the next packet
    uint8_t currAncillaryPkt[ANCSIZE];
    uint8_t nextAncillaryPkt[ANCSIZE];

    uint8_t **nextOciPacketBuffer = new uint8_t *[MAX_NUM_PACKETS];
    nextOciPacketBuffer[0] = new uint8_t[MAX_PACKET_SIZE * MAX_NUM_PACKETS];
    for (size_t i = 1; i < MAX_NUM_PACKETS; i++)
        nextOciPacketBuffer[i] = nextOciPacketBuffer[i - 1] + MAX_PACKET_SIZE;

    uint8_t **currOciPacketBuffer = new uint8_t *[MAX_NUM_PACKETS];
    currOciPacketBuffer[0] = new uint8_t[MAX_PACKET_SIZE * MAX_NUM_PACKETS];
    for (size_t i = 1; i < MAX_NUM_PACKETS; i++)
        currOciPacketBuffer[i] = currOciPacketBuffer[i - 1] + MAX_PACKET_SIZE;

    // get the number of packets for the next packet
    uint32_t nextNumPackets;
    int32_t nextPktAncillaryIndex;
    vector<int32_t> nextTelemetryIndices;

    //! updates when spin contains science data
    int32_t lastGoodSpinNum;
    int32_t nextSpinNum;
    int32_t currSpinNum;

    // Get first science or ancillary packet

    // open up the first packet and check the length of the packet body.
    // if the packet can't fit into the buffer, it will exit the program
    verify_packet_len(&tfileStream, currPacketLength, apid, isEndFile, isSPW);

    // after verifying that the packet will fit into the buffer, load it into the buffer
    readSinglePacket(&tfileStream, packetBuffer, currPacketLength, apid, isEndFile, isSPW);

    int numPacketsSkipped = 0;

    // Continue to read until the right apid is found
    while (apid != 636 && apid != 700 && apid != 720 && !isEndFile) {
        numPacketsSkipped++;

        // previous read moved the pointer to a new location, verify it again before reading
        verify_packet_len(&tfileStream, currPacketLength, apid, isEndFile, isSPW);
        readSinglePacket(&tfileStream, packetBuffer, currPacketLength, apid, isEndFile, isSPW);
    }
    if (isEndFile) {
        cout << "No science packets found in file" << endl;
        exit(EXIT_FAILURE);
    }
    if (numPacketsSkipped > 0)
        cout << numPacketsSkipped << " packets skipped" << endl;

    // Read first scan and check for ancillary packet
    nextPktAncillaryIndex = -1;
    nextTelemetryIndices.clear();

    int32_t ancillaryPktYear, ancillaryPktDay;
    double ancillaryPktTime;

    // Get granule duration in minutes
    // takes the granule mins from the commandline (a string) and save it
    // as an int to be used
    int32_t granuleMins;
    istringstream(granuleLengthStr) >> granuleMins;

    spatialAggTable spatialAggList[10];

    
    /**
     * when evaluating L0 packets, there can be a data type change. These are all the possible data types.
     * 
     * index 0 no data type. used in the instrument configuration to indicate segments of the scan where 
     * data are not collected. 
     * 
     * index 1 is science data (Earth viewing data). No suffix bc all of our l1a files are this.
     * 
     * index 10 is filler and we don't use it. Included so the other data type ids line up.
     * 
     * NBSB (Non-Baseline Spectral Bands) is NOT a instrument data type but we have it to note if the science data is collected
     * using the non standard spectral config. 
     * 
     */
    string const DATA_TYPES[] = {"",      "",      "_DARK", "_SOL", "_SPCA",   "_LIN",    "_LUN",
                                 "_DIAG", "_STAT", "_SPEC", "",     "_SNAP-X", "_SNAP-I", "_NBSB"};

    string const SWIR_DATA_MODES[] = {"", "_SDIAG", "_SRAW", "_STEST"};

    //! a copy that gets unpacked and thrown away later
    uint8_t **tempOciPacketBuffer = new uint8_t *[MAX_NUM_PACKETS];
    tempOciPacketBuffer[0] = new uint8_t[MAX_PACKET_SIZE * MAX_NUM_PACKETS];
    for (size_t i = 1; i < MAX_NUM_PACKETS; i++)
        tempOciPacketBuffer[i] = tempOciPacketBuffer[i - 1] + MAX_PACKET_SIZE;

    int8_t *ccdLineError = new int8_t[MAX_NUM_PACKETS];
    uint8_t *sciPacketSequenceError = new uint8_t[MAX_NUM_PACKETS];
    for (size_t i = 0; i < MAX_NUM_PACKETS; i++) {
        ccdLineError[i] = 0;
        sciPacketSequenceError[i] = 0;
    }
    uint8_t noSequenceErrorFlag = 255;
    uint16_t swirMode = 0;

    // l1afile manager so if there's a data type change, dont close the files
    L1aFileManager l1aFileManager = L1aFileManager();

    // Before processing, find the next good spin number that contains an ancillary index
    // Using "next" ancilliary packet var because we have not started processessing
    while ((nextPktAncillaryIndex == -1) && !isEndFile) {
        readScanPackets(&tfileStream, packetBuffer, (uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0],
                        nextNumPackets, nextSpinNum, nextPktAncillaryIndex, nextTelemetryIndices,
                        noSequenceErrorFlag, isEndFile, isSPW);
    }

    if (isEndFile) {
        cout << "No ancillary packets found in file" << endl;
        l1aFileManager.dumpOutlistBuffer(outlist);
        exit(110);  // ver 0.99.21
    }

    // time variables
    int32_t granuleStartTime = -1;  // time used to name the file
    int32_t granuleMaxTime = -1;
    uint16_t maxScans;
    double prevAncillaryPktTime = 0.0;

    // Each granule is ~5 mins = 300 seconds and it has ~1710 scans or lines.
    // 1710/300 ~= 5.7 scans per second. 1/5.7 ~= 0.1755 seconds per scan.
    const double SECONDS_PER_SCAN = 0.1755;

    /** Max swir pixels allowed  */
    uint16_t maxNumOfSwirPixels = 4060;

    // Pixel Variables
    uint16_t numCcdPixels = 0;
    uint16_t numDarkCcdPixels = 0;
    uint16_t numSwirPixels = 0;
    uint16_t numDarkSwirPixels = 0;

    uint16_t numBlueBands = 0;
    uint16_t numRedBands = 0;

    // 1 = enabled
    uint16_t isBlueCcdTapsEnabled[16] = {};
    uint16_t isRedCcdTapsEnabled[16] = {};

    L1aFile* currentL1aFile = nullptr; 

    // Packets will continue to be read as long as the ancillary data collection fields did not change.
    // Resets to 1 when processing scans for a new packet. Initializes to 1 in the main loop.
    int ancillaryDataTbleNotModified; 

    // records the time of each new ancillary packet
    double lastAncillaryPktTime = 0.0; 

    ////////////////////// Main Loop ///////////////////
    while (!isEndFile) {
        // Check for zero science pixels
        if (nextPktAncillaryIndex != -1) {
            memcpy(currAncillaryPkt, &nextOciPacketBuffer[nextPktAncillaryIndex][0], ANCSIZE);
            getBandDimensions(currAncillaryPkt, numCcdPixels, numBlueBands, numRedBands, numSwirPixels,
                              numDarkCcdPixels, numDarkSwirPixels, isBlueCcdTapsEnabled, isRedCcdTapsEnabled,
                              spatialAggList);
        }

        while ((numCcdPixels == 1 || nextPktAncillaryIndex == -1) && !isEndFile) {
            if (nextPktAncillaryIndex != -1)
                cout << "Ancillary packet has zero science pixels at spin " << nextSpinNum << endl;
            readScanPackets(&tfileStream, packetBuffer,
                            (uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0], nextNumPackets,
                            nextSpinNum, nextPktAncillaryIndex, nextTelemetryIndices, noSequenceErrorFlag,
                            isEndFile, isSPW);
            if (nextPktAncillaryIndex != -1) {
                memcpy(currAncillaryPkt, &nextOciPacketBuffer[nextPktAncillaryIndex][0], ANCSIZE);
                getBandDimensions(currAncillaryPkt, numCcdPixels, numBlueBands, numRedBands, numSwirPixels,
                                  numDarkCcdPixels, numDarkSwirPixels, isBlueCcdTapsEnabled,
                                  isRedCcdTapsEnabled, spatialAggList);
            }
        }

        if (nextPktAncillaryIndex == -1) {
            cout << "No ancillary packets for last granule" << endl;
            break;
        }

        cout << endl;
        cout << "number of CCD band pixels (ccd_pixels): " << numCcdPixels << endl;
        cout << "number of SWIR band pixels:             " << numSwirPixels << endl;
        cout << "number of dark collect pixels:          " << numDarkCcdPixels << endl;
        cout << "number of dark SWIR pixels:             " << numDarkSwirPixels << endl;
        cout << "number of blue bands:                   " << numBlueBands << endl;
        cout << "number of red bands:                    " << numRedBands << endl;

        // get the first scan time of the current packet and save it to prevAncillaryPktTime
        getAncillaryPacketTime(currAncillaryPkt, ancillaryPktYear, ancillaryPktDay, ancillaryPktTime);

        // Assign the first ancillary packet time
        prevAncillaryPktTime = ancillaryPktTime - SECONDS_PER_SCAN;

        // make object to hold all time info so this object can be passed instead of 3+ variables
        AncillaryPktTimeStamp ancillaryPktStartTime, ancillaryPktEndTime;
        ancillaryPktStartTime.year = ancillaryPktYear;
        ancillaryPktStartTime.day = ancillaryPktDay;
        ancillaryPktStartTime.second = ancillaryPktTime;

        /*
            --start_time was given, set the file time name (ltime) to be this time.
            granuleStartTime string will be cleared after so that if there's a non-science
            data inbetween the science data, the non-science data won't use --start_time
            and the next science data file will also not use it

        */
        if (isGranuleStartTimeGiven) {
            granuleStartTime =
                (int32_t)(granuleStartTimeIntStruct.tm_hour * 3600 + granuleStartTimeIntStruct.tm_min * 60 +
                          granuleStartTimeIntStruct.tm_sec);

            // adjust for different dates

            // passed in day is Day of Year (DOY), so ignore month ie: 229 DOY = Aug 17th
            int32_t startTimeAsJulianDay = jday(granuleStartTimeIntStruct.tm_year + 1900, 1, granuleStartTimeIntStruct.tm_yday + 1);
            int32_t ancillaryPktTimeAsJulianDay = jday(ancillaryPktYear, 1, ancillaryPktDay);
            int32_t julianDayDiff = startTimeAsJulianDay - ancillaryPktTimeAsJulianDay;
            granuleStartTime += julianDayDiff * 86400;

            // set max time for the L0 files and max scans depending on granule length.
            // if minutes is 0, then max time is +600 by default to read the entire file.
            // else, add (granule length * 60) to start_time to get max time for this file.
            granuleMaxTime = granuleMins > 0 ? granuleStartTime + granuleMins * 60 : granuleStartTime + 600;

            // update max scans based on granule's duration in minutes
            maxScans = granuleMins > 0 ? (uint16_t)((granuleMins * 60 / SECONDS_PER_SCAN) + 2) : 3600;

            // **flip so this does not run again **
            isGranuleStartTimeGiven = false;

            // moves packet pointer forward if the packets actual start time is less than what we want our
            // granule start time to be.
            while (((ancillaryPktTime < granuleStartTime) || (numCcdPixels == 1)) && !isEndFile) {

                // note previous packet time
                lastAncillaryPktTime = ancillaryPktTime;

                readScanPackets(&tfileStream, packetBuffer,
                                (uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0], nextNumPackets,
                                nextSpinNum, nextPktAncillaryIndex, nextTelemetryIndices, noSequenceErrorFlag,
                                isEndFile, isSPW);
                if (nextPktAncillaryIndex != -1) {
                    memcpy(currAncillaryPkt, &nextOciPacketBuffer[nextPktAncillaryIndex][0], ANCSIZE);
                    getAncillaryPacketTime(currAncillaryPkt, ancillaryPktYear, ancillaryPktDay,
                                           ancillaryPktTime);

                    getBandDimensions(currAncillaryPkt, numCcdPixels, numBlueBands, numRedBands,
                                      numSwirPixels, numDarkCcdPixels, numDarkSwirPixels,
                                      isBlueCcdTapsEnabled, isRedCcdTapsEnabled, spatialAggList);
                    if (numCcdPixels == 1)
                        cout << "Ancillary packet has zero science pixels at spin " << nextSpinNum << endl;
                }

                // check for crossing date line (ancillary time becomes smaller bc of a new day) and if the 
                // granuleStartTime > 86400 (max seconds in a day). when the adjusted granuleStartTime via jday()
                // is > 86400, then the start of the ancillary packet time is in the previous day by 
                // (granduleStartTime - 86400) seconds. 
                if (ancillaryPktTime < lastAncillaryPktTime && granuleStartTime > 86400) {
                    granuleStartTime = granuleStartTime % 86400;
                }
            }

            if (isEndFile) {
                cout << "No science data in time range" << endl;
                l1aFileManager.dumpOutlistBuffer(outlist);
                return 110;
            }

            ancillaryPktStartTime.second = ancillaryPktTime;
        }

        /*  --start_time is not set, use the first scantime of the packet.
            This runs if there is a **non-science data** packet in between OCI data.
            
            If the science data gets split because a non-science data packet is inbewteen it,
            when it gets back to the science data, it will append to the same file with the same
            scan time. 

            So if the L0 has the following data types and the following times: 
                OCI (03:00) --> SOL (03:01) --> SPCA (03:02) --> OCI (03:04)
                
            The first and second OCI file will use time (03:00) because the first file will be
            kept open. SOL, SPCA will have their own file with their time as well. Note that 
            not all non-science will append to the same file like OCI.
        */
        else {
            // no start time given by user, use packet's first ancillary packet time 
            granuleStartTime = (int32_t)ancillaryPktTime;
            
            // user did not speciify how long the granuel is suppose to be for the set of L0 files,
            // so  process everything. When there is a switch in the data type, update the granule
            // max time for every new ancillary packet
            if (granuleMins == 0) {
                granuleMaxTime = granuleStartTime + (15 * 60); // 15 mins
            }
            // otherwise, the granule length is specified. Only calculate the max time once
            // so the program ends when that time is reached even if a new ancillary packet is read
            else if (granuleMins > 0 && granuleMaxTime == -1) {
                granuleMaxTime = granuleStartTime + (granuleMins * 60);
            }
            // update max scans based on granule_len
            maxScans = granuleMins > 0 ? (uint16_t)((granuleMins * 60 / SECONDS_PER_SCAN) + 2) : 3600;
        }

        // number of SWIR bands
        numSwirPixels = ((numSwirPixels + 7) / 8) *
                        8;  // Round up SWIR number of pixels to a multiple of 8,  LH, 4/15/2022, v1.03.00
        unsigned short numSwirBands = 9;

        // reset to not modified for the next packet
        ancillaryDataTbleNotModified = 1;

        int16_t ccdBandLineIndices[MAX_LINES];
        int16_t swirBandLineIndices[MAX_LINES * numSwirBands];
        int16_t ccdBandDarkLineIndices[MAX_LINES];
        int16_t swirBandDarkLineIndices[MAX_LINES];
        int16_t swirLineOffset[9];

        for (size_t i = 0; i < MAX_LINES; i++) {
            ccdBandLineIndices[i] = -1;
            for (size_t j = 0; j < numSwirBands; j++)
                swirBandLineIndices[i * numSwirBands + j] = -1;
            ccdBandDarkLineIndices[i] = -1;
            swirBandDarkLineIndices[i] = -1;
        }

        if (!swirLineOffsetString.empty()) {
            if (swirLineOffsetString.compare("ETU") == 0) {
                memcpy(swirLineOffset, SWIR_LINE_OFFSET_ETU, 9 * sizeof(int16_t));
            } else {
                vector<string> parts;
                boost::split(parts, swirLineOffsetString, boost::is_any_of(","));
                if (parts.size() != 9) {
                    cout << "Processing with single file option" << endl;
                    exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < 9; i++) {
                    swirLineOffset[i] = atoi(parts[i].c_str());
                }
            }
        } else {
            memcpy(swirLineOffset, SWIR_LINE_OFFSET_DEFAULT, 9 * sizeof(int16_t));
        }

        makeOciLineIndex(spatialAggList, ccdBandLineIndices, swirBandLineIndices, ccdBandDarkLineIndices,
                         swirBandDarkLineIndices, swirLineOffset);

        // Get SWIR band data mode
        // uint16_t smode = 0;      // LH, 09/10/2020
        getSwirMode((uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0], nextNumPackets, swirMode);
        uint16_t initialSwirMode = swirMode;
        int isSameSwirMode = 1;

        // (0 = enabled, 1 = reset, 2 = video)
        uint16_t cdsMode;

        // Determine start and end time of granule
        uint32_t jd0 = jday(ancillaryPktYear, 1, ancillaryPktDay);
        int16_t yr16 = (int16_t)ancillaryPktYear;
        int16_t doy = (int16_t)ancillaryPktDay;
        int16_t month, day;
        yd2md(yr16, doy, &month, &day);

        // ---- SETTING TIME NAMING ----
        int32_t granuelStartTimeHour = (int32_t)(granuleStartTime / 3600);
        int32_t granuelStartTimeMins = (int32_t)((granuleStartTime - granuelStartTimeHour * 3600) / 60);
        int32_t granuelStartTimeSecs =
            (int32_t)(granuleStartTime - granuelStartTimeHour * 3600 - granuelStartTimeMins * 60);

        stringstream timeString, dateString;

        // set time in the format: HHMMSS
        timeString << setfill('0') << setw(2) << granuelStartTimeHour << setfill('0') << setw(2)
                   << granuelStartTimeMins << setfill('0') << setw(2) << granuelStartTimeSecs;

        // set date in the format: YYYYMMDD
        dateString << setfill('0') << setw(4) << ancillaryPktYear << setfill('0') << setw(2) << month
                   << setfill('0') << setw(2) << day;

        string baseName;
        char *tmpFileStr = strdup(ifileName.c_str());  // need this since basename may modify the char* passed in
        baseName.assign(basename(tmpFileStr));
        free(tmpFileStr);

        if (outfilePrefixTag == "")
            outfilePrefixTag.assign("PACE_OCI");
        short dataType = spatialAggList[1].dataType;  // ;  Get data type for file name and metadata
        short maxdtype = 2;
        for (size_t i = 0; i < 10; i++) {
            if (spatialAggList[i].dataType > maxdtype)
                maxdtype = spatialAggList[i].dataType;
        }
        // if (maxdtype > 2) dtype = maxdtype;
        if ((maxdtype != 2) && (maxdtype != 10)) {
            dataType = maxdtype;
            cout << "\nWARNING: Non-Science Data is now being processed. Type: " << DATA_TYPES[dataType]
                 << "\n"
                 << endl;
        }

        // data type mod for ETU before June 2020
        if ((jd0 < 2459000) && dataType == 11)
            dataType = 9;

        // set the l1a output file name. If there is not custom name given by the user, then name it
        // normally as: {tag}.{date}T{time}.L1A.nc (ie. PACE_OCI.20240323T003546.L1A.nc)
        // or with a type: {tag}_{type}.{date}T{time}.L1A.nc (ie. PACE_OCI_SPEC.20240323T003908.L1A.nc)
        string l1aFileName = ofileName;
        if (ofileName.compare("") == 0) {
            l1aFileName = outfilePrefixTag + DATA_TYPES[dataType] + SWIR_DATA_MODES[swirMode];
            // keep input filename substrings for OCI instrument test data
            if (!isSPW)
                l1aFileName += "_" + baseName.substr(0, 4) + baseName.substr(5, 3) + baseName.substr(9, 3);
            l1aFileName += "." + dateString.str() + "T" + timeString.str() + ".L1A.nc";
        }


        // Initialize data arrays

        // blue, red, swir band dark collect data for granule (Dark Calibration (DC) in NetCDF file)
        uint16_t **blueDarkCalibrationData = new uint16_t *[maxScans];
        uint16_t **redDarkCalibrationData = new uint16_t *[maxScans];
        uint32_t **swirDarkCalibrationData = new uint32_t *[maxScans];
        initializeDarkCalibrationArray(blueDarkCalibrationData, numDarkCcdPixels, numBlueBands, maxScans);
        initializeDarkCalibrationArray(redDarkCalibrationData, numDarkCcdPixels, numRedBands, maxScans);
        initializeDarkCalibrationArray(swirDarkCalibrationData, numDarkSwirPixels, numSwirBands, maxScans);

        uint16_t **blueCcdDarkData = new uint16_t *[numBlueBands];
        uint16_t **redCcdDarkData = new uint16_t *[numRedBands];
        uint32_t **swirDarkData = new uint32_t *[numSwirBands];

        // make and initialize science data arrays
        uint16_t **blueScienceData = new uint16_t *[numBlueBands];
        uint16_t **redScienceData = new uint16_t *[numRedBands];
        uint32_t **swirScienceData = new uint32_t *[numSwirBands];
        initializeScienceDataArray(blueScienceData, numCcdPixels, numBlueBands);
        initializeScienceDataArray(redScienceData, numCcdPixels, numRedBands);
        initializeScienceDataArray(swirScienceData, numSwirPixels, numSwirBands);

        uint8_t **ancillaryData = new uint8_t *[maxScans + 1];
        ancillaryData[0] = new uint8_t[ANCSIZE * (maxScans + 1)];
        for (size_t i = 1; i < (size_t)(maxScans + 1); i++)
            ancillaryData[i] = ancillaryData[i - 1] + ANCSIZE;
        
        // initialize to be 0 so when writing, it does not grab garbage for instances
        // where ancillary apid is not found. Cant use BAD_INT because it is negative
        for (int i = 0; i < ANCSIZE * (maxScans + 1); i++) {
            ancillaryData[0][i] = 0;
        }

        uint32_t const maxNumTelemetryPkt = maxScans * 10;
        uint8_t **tlmdata = new uint8_t *[maxNumTelemetryPkt];
        tlmdata[0] = new uint8_t[TLMSIZE * (maxNumTelemetryPkt)];
        for (size_t i = 1; i < (size_t)(maxNumTelemetryPkt); i++)
            tlmdata[i] = tlmdata[i - 1] + TLMSIZE;

        // Round up SWIR number of pixels to a multiple of 8
        uint16_t numSwirPixelsRounded = ((numSwirPixels + 7) / 8) * 8;

        // Note: "bands" arrays here are the reverse order as IDL versions

        uint16_t **l0BlueBandData = new uint16_t *[numCcdPixels];
        uint16_t **l0RedBandData = new uint16_t *[numCcdPixels];
        for (size_t i = 0; i < numCcdPixels; i++) {
            l0BlueBandData[i] = NULL;
            l0RedBandData[i] = NULL;
        }

        // Blue pixel line number from spacecraft
        int16_t *blueCcdLineNums = new int16_t[numCcdPixels];
        // Red pixel line number from spacecraft
        int16_t *redCcdLineNums = new int16_t[numCcdPixels];

        for (size_t i = 0; i < numCcdPixels; i++) 
            blueCcdLineNums[i] = -1;
        for (size_t i = 0; i < numCcdPixels; i++)
            redCcdLineNums[i] = -1;

        // May be full-scan SWIR pixels in ETU data
        maxNumOfSwirPixels = 4096;

        uint32_t **l0SwirBandData = new uint32_t *[maxNumOfSwirPixels];
        for (size_t i = 0; i < maxNumOfSwirPixels; i++) {
            l0SwirBandData[i] = new uint32_t[numSwirBands];
        }

        // Swir Pixel line number. Swir's line number == ccd line number 8x spatial aggregation
        int16_t *swirLineNums = new int16_t[maxNumOfSwirPixels];

        // swir frame type for dark collect
        int8_t **swirDarkCalFrameTypeData = new int8_t *[maxScans];
        swirDarkCalFrameTypeData[0] = new int8_t[numDarkSwirPixels * maxScans];
        for (size_t i = 1; i < maxScans; i++)
            swirDarkCalFrameTypeData[i] = swirDarkCalFrameTypeData[i - 1] + numDarkSwirPixels;
        memset(swirDarkCalFrameTypeData[0], -1, sizeof(int8_t) * numDarkSwirPixels * maxScans);

        // directly from the L0 packet, reordered in ascending wavelength
        int8_t *l0SwirPixelFrameTypes = new int8_t[maxNumOfSwirPixels];
        int8_t *swirFrameTypesSciData = new int8_t[maxNumOfSwirPixels];


        string maxgap_string = to_string(maxGap);
        string noSPW_string = to_string(isSPW);

        // track if the file is being merge into another file of the same data type
        bool isFileBeingAppended = false;

        // if the file manager does not have a l1a for this data type, initialize output file and get object IDs for EV data
        // -- for lunear calibration, the file is removed after it is done. So Lunar Calibration will not be grouped because
        // it contains different ccd and swir pixel dimensions that will cause an indexing error
        if (!l1aFileManager.contains(dataType)) {
            cout << "Creating: " << l1aFileName.c_str() << endl;
            cout << "Starting at spin number " << nextSpinNum << endl;

            // create a new l1a file and initialize it
            currentL1aFile = l1aFileManager.createL1aOutputFile(dataType);
            currentL1aFile->initializeL1aFile((char *)l1aFileName.c_str(), maxScans, numCcdPixels, numBlueBands,
                                     numRedBands, numSwirPixels, numDarkCcdPixels);

            // write in what command and files were used to make this file
            currentL1aFile->writeProcessingControl(hktList, ifileName, granuleStartTimeStr, maxgap_string,
                                          outfilePrefixTag, swirLineOffsetString, outlist, l1aFileName, doi,
                                          pversion, noSPW_string, VERSION);
        }

        // file manager has the file, but the time gap between this packet and last is larger than requested
        // close the current one and create a new 
        else if (l1aFileManager.contains(dataType) && (ancillaryPktTime - lastAncillaryPktTime > maxFileGap)) {
            currentL1aFile = l1aFileManager.getL1aFile(dataType);

            cout << "Closing: " << currentL1aFile->getFileName() << " because max gap between packets exceeds " << maxFileGap/60 << " mins" << endl;
            l1aFileManager.closeAndRemoveFile(dataType);

            // update the outlist buffer. removed to prevent outlist merging
            l1aFileManager.processPrevFile();

            cout << "Creating: " << l1aFileName.c_str() << endl;
            cout << "Starting at spin number " << nextSpinNum << endl;

            // create a new l1a file and initialize it
            currentL1aFile = l1aFileManager.createL1aOutputFile(dataType);
            currentL1aFile->initializeL1aFile((char *)l1aFileName.c_str(), maxScans, numCcdPixels, numBlueBands,
                                     numRedBands, numSwirPixels, numDarkCcdPixels);

            // write in what command and files were used to make this file
            currentL1aFile->writeProcessingControl(hktList, ifileName, granuleStartTimeStr, maxgap_string,
                                          outfilePrefixTag, swirLineOffsetString, outlist, l1aFileName, doi,
                                          pversion, noSPW_string, VERSION);
        }

        // if file manager has the file, then switch the current working l1afile to be the current dataType
        else {
            cout << "Creating: " << l1aFileName.c_str() << endl;
            cout << "          " << l1aFileName.c_str() << " will be appened into " << currentL1aFile->getFileName() 
                << " because it shares the same data-type:  " << DATA_TYPES[dataType] << "." << endl;
            cout << "Starting at spin number " << nextSpinNum << endl;
            isFileBeingAppended = true;
        }

        // Read and process OCI scans
        int completeFlag = 0;
        int scienceDataStatus = 0;
        uint32_t currScan = 0;
        //! current packet + next packet
        uint32_t totalNumPackets;
        uint32_t currNumPackets;
        int32_t currPktAncillaryIndex;
        vector<int32_t> currTelemetryIndices;

        // nextSpinNum is the first good spin number when searching through packets
        // it has a nextPktAncillaryIndex, etc.
        lastGoodSpinNum = nextSpinNum;
        int32_t endData = 0;
        int spinGap = 1;
        uint32_t numTelemetryPackets = 0;  // number of telemtry packets

        uint16_t btype, rtype;
        uint16_t blueSpectralAggregations[16];
        uint16_t redSpectralAggregations[16];
        uint16_t nbands;

        // Dimension sizes for the current file
        DimensionShape *currFileDimShape = l1aFileManager.getCurrL1aFileDimShape(dataType);

        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------

        while (ancillaryPktTime < granuleMaxTime && ancillaryDataTbleNotModified && !endData &&
               isSameSwirMode && (spinGap <= maxGap) && currScan < maxScans) {
            // sometimes, data contain more scans than possible, set currScan<maxsc to avoid overflow, LH
            // 8/26/2020
            if ((currScan % 100) == 0) {
                cout << "Processing scan " << currScan << endl;
            }

            // load next packet into the current packet
            memcpy(&currOciPacketBuffer[0][0], &nextOciPacketBuffer[0][0], MAX_PACKET_SIZE * MAX_NUM_PACKETS);

            // update all current variables to be next so the next packet can be processed and a new
            // packet can be read in from the L0 file stream
            currPktAncillaryIndex = nextPktAncillaryIndex;
            currTelemetryIndices = nextTelemetryIndices;
            currNumPackets = nextNumPackets;
            currSpinNum = nextSpinNum;
            endData = isEndFile;

            // after making the next packets to be the current, read and load in the next, next packets
            readScanPackets(&tfileStream, packetBuffer,
                            (uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0], nextNumPackets,
                            nextSpinNum, nextPktAncillaryIndex, nextTelemetryIndices,
                            sciPacketSequenceError[currScan], isEndFile, isSPW);

            // Disabled maxsc check for I&T
            if (spinGap >= 1 && currNumPackets > 1) {  // ver 1.03.00
                // Save ancillary packet in array
                if (currPktAncillaryIndex != -1) {
                    memcpy(ancillaryData[currScan], currOciPacketBuffer[currPktAncillaryIndex], ANCSIZE);
                } else {
                    // insert spin # if missing anc packet
                    uint32_t ui32 = 0;
                    swapc_bytes2((char *)&currSpinNum, (char *)&ui32, sizeof(int32_t), 1);
                    memcpy(&ancillaryData[currScan][24], &ui32, 4);

                    // set apid to 0 where there is missing anc packet, LH 4/28/2020
                    ui32 = 0;
                    memcpy(&ancillaryData[currScan][0], &ui32, 4);
                }

                // Save telemetry packet in array
                int numTelemetryIndices = currTelemetryIndices.size();
                if (numTelemetryIndices > 0 && numTelemetryPackets < maxNumTelemetryPkt) {
                    int telemetryIndex = 0;
                    while (telemetryIndex < numTelemetryIndices && numTelemetryPackets < maxNumTelemetryPkt) {
                        memcpy(tlmdata[numTelemetryPackets],
                               currOciPacketBuffer[currTelemetryIndices[telemetryIndex]], TLMSIZE);
                        numTelemetryPackets++;
                        telemetryIndex++;
                    }
                    if (numTelemetryPackets >= maxNumTelemetryPkt)
                        cout << "Maximum number of telemetry packets exceeded at spin: " << currSpinNum
                             << endl;
                }

                // Zero out swirLineNums, l0SwirPixelFrameTypes, and l0SwirBandData arrays
                for (size_t i = 0; i < maxNumOfSwirPixels; i++) {
                    swirLineNums[i] = -1;
                    l0SwirPixelFrameTypes[i] = 0;
                    for (size_t j = 0; j < numSwirBands; j++)
                        l0SwirBandData[i][j] = 0;
                }

                // Unpack science data from packets
                totalNumPackets = currNumPackets + nextNumPackets;
                if (totalNumPackets > MAX_NUM_PACKETS) {  // LH, 5/27/2022, ver 1.04.01
                    cout << "Packets exceed max [" << MAX_NUM_PACKETS << "]." << endl;
                    totalNumPackets = MAX_NUM_PACKETS;
                    nextNumPackets = totalNumPackets - currNumPackets;
                }
                memcpy(&tempOciPacketBuffer[0][0], &currOciPacketBuffer[0][0],
                       MAX_PACKET_SIZE * currNumPackets);
                if (nextNumPackets > 0)
                    memcpy(&tempOciPacketBuffer[currNumPackets][0], &nextOciPacketBuffer[0][0],
                           MAX_PACKET_SIZE * nextNumPackets);  // tempOciPacketBuffer[*,currNumPackets:totalNumPackets-1]
                                                               // = pbuffer2[*,0:npkt2-1]

                unpackScienceData(totalNumPackets, currSpinNum, numCcdPixels, numSwirPixelsRounded,
                                  maxNumOfSwirPixels, nbands, isBlueCcdTapsEnabled, isRedCcdTapsEnabled,
                                  (uint8_t(*)[MAX_PACKET_SIZE]) & tempOciPacketBuffer[0][0], l0BlueBandData,
                                  l0RedBandData, l0SwirBandData, blueCcdLineNums, redCcdLineNums,
                                  swirLineNums, btype, blueSpectralAggregations, rtype,
                                  redSpectralAggregations, l0SwirPixelFrameTypes, scienceDataStatus);
                if (scienceDataStatus != 0)
                    cout << "Science data unpacking error in spin: " << currSpinNum << endl;
                completeFlag = completeFlag | scienceDataStatus;

                blueCcdDarkData[0] = blueDarkCalibrationData[currScan];
                for (size_t i = 1; i < numBlueBands; i++)
                    blueCcdDarkData[i] = blueCcdDarkData[i - 1] + numDarkCcdPixels;

                redCcdDarkData[0] = redDarkCalibrationData[currScan];
                for (size_t i = 1; i < numRedBands; i++)
                    redCcdDarkData[i] = redCcdDarkData[i - 1] + numDarkCcdPixels;

                swirDarkData[0] = swirDarkCalibrationData[currScan];
                for (size_t i = 1; i < numSwirBands; i++)
                    swirDarkData[i] = swirDarkData[i - 1] + numDarkSwirPixels;

                // Initialize science/dark array with fill value
                std::fill(blueScienceData[0],
                          blueScienceData[0] + numCcdPixels * ((numBlueBands > 0) ? numBlueBands : 1), 65535);
                std::fill(redScienceData[0],
                          redScienceData[0] + numCcdPixels * ((numRedBands > 0) ? numRedBands : 1), 65535);
                std::fill(swirScienceData[0],
                          swirScienceData[0] + numSwirPixels * ((numSwirBands > 0) ? numSwirBands : 1),
                          1048575);
                memset(swirFrameTypesSciData, -1, sizeof(int8_t) * maxNumOfSwirPixels);

                // Initialize with fill value, LH 8/6/2020
                std::fill(blueCcdDarkData[0],
                          blueCcdDarkData[0] + numDarkCcdPixels * ((numBlueBands > 0) ? numBlueBands : 1),
                          65535);
                std::fill(redCcdDarkData[0],
                          redCcdDarkData[0] + numDarkCcdPixels * ((numRedBands > 0) ? numRedBands : 1),
                          65535);
                std::fill(swirDarkData[0],
                          swirDarkData[0] + numDarkSwirPixels * ((numSwirBands > 0) ? numSwirBands : 1),
                          1048575);

                // If science data in scan, check data for gaps or inconsistencies and determine data types
                if ((blueCcdLineNums[0] != -1) || (redCcdLineNums[0] != -1) || (swirLineNums[0] != -1)) {
                    checkAndLoadScienceData(
                        dataType, numCcdPixels, numSwirPixels, numDarkCcdPixels, numDarkSwirPixels,
                        numBlueBands, numRedBands,
                        numSwirBands,  // maxNumOfSwirPixels -> numSwirPixels, LH, ver 1.03.00
                        ccdBandLineIndices, swirBandLineIndices, ccdBandDarkLineIndices,
                        swirBandDarkLineIndices, l0BlueBandData, l0RedBandData, l0SwirBandData,
                        blueCcdLineNums, redCcdLineNums, swirLineNums, blueScienceData, redScienceData,
                        swirScienceData, blueCcdDarkData, redCcdDarkData, swirDarkData,
                        ccdLineError[currScan], scienceDataStatus);
                    if (scienceDataStatus >= 2)
                        cout << "Science data checking error in spin: " << currSpinNum << endl;
                    completeFlag = completeFlag | scienceDataStatus;

                    // Store SWIR band frame types
                    for (size_t i = 0; i < maxNumOfSwirPixels; i++) {
                        // LH, 8/24/2020; use SWIR band0, LH 11/20/2020
                        if (swirBandLineIndices[swirLineNums[i] * numSwirBands] > -1) {
                            swirFrameTypesSciData[swirBandLineIndices[swirLineNums[i] * numSwirBands]] =
                                l0SwirPixelFrameTypes[i];
                        }
                    }

                    for (size_t i = 0; i < maxNumOfSwirPixels; i++) {
                        if (swirBandDarkLineIndices[swirLineNums[i]] > -1) {  // LH, 8/24/2020
                            swirDarkCalFrameTypeData[currScan][swirBandDarkLineIndices[swirLineNums[i]]] =
                                l0SwirPixelFrameTypes[i];
                        }
                    }

                    // Write science data to file
                    currentL1aFile->writeScienceData(currFileDimShape, currScan, numBlueBands, numRedBands, numSwirBands,
                                                numCcdPixels, numSwirPixels, blueScienceData, redScienceData,
                                                swirScienceData, swirFrameTypesSciData);
                    //  Check for instrument HKT packets in scan (placeholder for now)

                    currScan++;
                    prevAncillaryPktTime = ancillaryPktTime;
                    ancillaryPktEndTime.year = ancillaryPktYear;
                    ancillaryPktEndTime.day = ancillaryPktDay;
                    ancillaryPktEndTime.second = ancillaryPktTime;
                    while (ancillaryPktEndTime.second > 86400)
                        ancillaryPktEndTime.second -= 86400;
                    lastGoodSpinNum = currSpinNum;
                    lastGoodSpinNum = currSpinNum;
                }

            }
            // No spin gap and has more than 1 packet
            else {
                int32_t packetHour = (int32_t)(ancillaryPktTime / 3600);
                int32_t packetMins = (int32_t)((ancillaryPktTime - packetHour * 3600) / 60);
                int32_t packetSecs = (int32_t)(ancillaryPktTime - packetHour * 3600 - packetMins * 60);

                // if the current ancillary packet time is less than the previous,
                // then the scans are out of order
                if (ancillaryPktTime <= prevAncillaryPktTime && currPktAncillaryIndex != -1) {
                    cout << "Scan " << currSpinNum << " out of order at" << packetHour << " " << packetMins
                         << " " << packetSecs << endl;
                }
            }  // if (dspn >= 1 && totalNumPackets >

            // Get scan time and band dimensions for next scan
            if (!endData && nextPktAncillaryIndex != -1) {
                memcpy(nextAncillaryPkt, &nextOciPacketBuffer[nextPktAncillaryIndex][0], ANCSIZE);
                getAncillaryPacketTime(nextAncillaryPkt, ancillaryPktYear, ancillaryPktDay, ancillaryPktTime);
                uint32_t jd = jday(ancillaryPktYear, 1, ancillaryPktDay);
                ancillaryPktTime += (jd - jd0) * 86400;

                getSwirMode((uint8_t(*)[MAX_PACKET_SIZE]) & nextOciPacketBuffer[0][0], nextNumPackets,
                            swirMode);

                // check that SWIR mode has not changed since getting the mode on the first packet and
                // report the change if it did
                isSameSwirMode = (swirMode == initialSwirMode);
                if (!isSameSwirMode) {
                    int32_t packetHour = (int32_t)(ancillaryPktTime / 3600);
                    int32_t packetMins = (int32_t)((ancillaryPktTime - packetHour * 3600) / 60);
                    int32_t packetSecs = (int32_t)(ancillaryPktTime - packetHour * 3600 - packetMins * 60);
                    cout << "SWIR data mode change at: " << packetHour << " " << packetMins << " "
                         << packetSecs << endl;
                }
                ancillaryDataTbleNotModified = compareAncillaryPackets(currAncillaryPkt, nextAncillaryPkt);
            }

            // next spin was read at the start of the loop.
            // lastGoodSpin updates to the currSpinNum **if the curr spin has science data.**
            // otherwise, it doesnt update; curr and nextSpin gets incremeneted and there will be a gap
            spinGap = nextSpinNum - lastGoodSpinNum;
            if (spinGap > maxGap)
                cout << "Spin number gap at spin  " << lastGoodSpinNum << endl;

        }  // while (ancillaryPktTime < mtime && acomp && !endData && isSameSwirMode && (dspn <= maxGap) &&
           // currScan < maxsc)

        if (granuleMins > 0 && currScan < (size_t)(maxScans - 10) && ancillaryDataTbleNotModified) {
            completeFlag = completeFlag | 1;
        }

        cout << "Scans in file: " << currScan << endl;
        cout << "Complete flag: " << completeFlag << endl;

        if (currScan > 0) {
            // Need to include ancillary packet from next scan
            if (nextPktAncillaryIndex != -1)
                memcpy(ancillaryData[currScan], nextAncillaryPkt, ANCSIZE);

            // Write calibration data to file
            currentL1aFile->writeCalibrationData(currFileDimShape, currScan, numBlueBands, numRedBands, numSwirBands,
                                            numDarkCcdPixels, numDarkSwirPixels, blueDarkCalibrationData[0],
                                            redDarkCalibrationData[0], swirDarkCalibrationData[0],
                                            swirDarkCalFrameTypeData[0]);

            // Get scan metadata and write to file
            // Save spid id computed in this function for later use
            int32_t *spinID = new int32_t[currScan + 100];
            for (size_t i = 0; i < (currScan + 100); i++)
                spinID[i] = BAD_INT;


            currentL1aFile->writeScanMetaData(currFileDimShape, currScan, ancillaryData[0], sciPacketSequenceError, ccdLineError,
                                         spinID, ancillaryPktStartTime);

            currentL1aFile->writeAncillaryData(currFileDimShape, currScan, ancillaryData[0]);

            // Unpack OCI telemetry data and write to file
            cdsMode = 0;
            int numTelemetryIndices = nextTelemetryIndices.size();
            if (nextTelemetryIndices.size() != 0 && spinGap <= maxGap && !isEndFile) {  // 0.99.20
                int telemetryIndex = 0;
                while (telemetryIndex < numTelemetryIndices && numTelemetryPackets < maxNumTelemetryPkt) {
                    memcpy(tlmdata[numTelemetryPackets],
                           nextOciPacketBuffer[nextTelemetryIndices[telemetryIndex]], TLMSIZE);
                    numTelemetryPackets++;
                    telemetryIndex++;
                }
            }

            if (numTelemetryPackets > 0) {
                cout << numTelemetryPackets << " HKT packets" << endl;
                currentL1aFile->writeTelemetryData(currFileDimShape, spatialAggList, numTelemetryPackets,
                                              (uint8_t(*)[TLMSIZE]) & tlmdata[0][0], spinID, cdsMode,
                                              currScan, ancillaryPktStartTime);
            }
            delete[] spinID;

            // Locate navigation data and write to file
            if (hktList.compare("") != 0) {
                currentL1aFile->writeNavigationData(currFileDimShape, hktList, ancillaryPktStartTime, ancillaryPktEndTime);
            }

            // Generate granule metadata and write to file
            string sdir = "Ascending";  // Hard-code until we have spacecraft data
            string edir = "Ascending";  // Hard-code until we have spacecraft data
               
            // get ancillary packet start and end times as a string
            stringstream timeCoverageStart = makeTimeCoverageString(ancillaryPktStartTime);
            stringstream timeCoverageEnd = makeTimeCoverageString(ancillaryPktEndTime);

            // also writes the time converage start into the outlist file for the current file
            currentL1aFile->writeGlobalMetaData(timeCoverageStart, timeCoverageEnd, isFileBeingAppended, l1aFileName,
                                           sdir, edir,
                                           dataType,  // spatialAggList[1].dataType changed to dataType
                                           swirMode, cdsMode);
            
            // write the file information to the outlist buffer
            if (swirMode > 0) {
                l1aFileManager.addFileToOutlistBuffer(l1aFileName, timeCoverageStart.str(), 
                timeCoverageEnd.str(), completeFlag, "DIAG", dataType
                );
            } 
            else if (dataType == 0) {
                l1aFileManager.addFileToOutlistBuffer(l1aFileName, timeCoverageStart.str(), 
                timeCoverageEnd.str(), completeFlag, "_", dataType
                ); 
            } 
            // named data types
            else if (DATA_TYPES[dataType].substr(0, 1) == "_") {
                l1aFileManager.addFileToOutlistBuffer(l1aFileName, timeCoverageStart.str(), 
                timeCoverageEnd.str(), completeFlag, DATA_TYPES[dataType].substr(1), dataType
                );
            } 
            // normal OCI l1a file and any other data types without a name in DATA_TYPES
            else {
                l1aFileManager.addFileToOutlistBuffer(l1aFileName, timeCoverageStart.str(), 
                timeCoverageEnd.str(), completeFlag, "", dataType
                );
            }

            // note the last datatype used
            l1aFileManager.updateLastDataTypeSeen(dataType);

            // Write common global metadata
            currentL1aFile->writeGlobalAttributes(history, doi, pversion);
            
            // update the number_of_scans dim. This is used in multiple write methods
            // so it is done outside the write method, unlike the other dims.
            currFileDimShape->incrementNumScansShape(currScan);

        } else {
            // Remove 0-scan file
            int status = remove(l1aFileName.c_str());
            if (status == 0) {
                cout << "Removing 0-scan file: " << l1aFileName.c_str() << endl;
            } else {
                cout << "Error removing " << l1aFileName.c_str() << endl;
                exit(EXIT_FAILURE);
            }
        }

        // after writing all the data, increment each UNLIMITED dimensions so the next
        // write will know the index to put the next batch of data 
        
        // tlm_packets is incremented in writeTelemetryPackets because the amt of data being
        // written is within the function and you can't get it from outside of it without modifying
        // the function itself 


        // -- LUNAR CALIBRATION --
        //  Do not group lunar calibration files. So close and remove it so sebsequent lunar cal files
        //  will be in a file of its own and not be appended into the first one. 
        if (DATA_TYPES[dataType] == "_LUN") {
            l1aFileManager.closeAndRemoveFile(dataType);
        }




        delete[] blueDarkCalibrationData[0];
        delete[] blueDarkCalibrationData;
        delete[] redDarkCalibrationData[0];
        delete[] redDarkCalibrationData;
        delete[] swirDarkCalibrationData[0];
        delete[] swirDarkCalibrationData;

        delete[] blueCcdDarkData;
        delete[] redCcdDarkData;
        delete[] swirDarkData;

        delete[] blueScienceData[0];
        delete[] blueScienceData;
        delete[] redScienceData[0];
        delete[] redScienceData;
        delete[] swirScienceData[0];
        delete[] swirScienceData;

        for (size_t i = 0; i < numCcdPixels; i++)
            if (l0BlueBandData[i] != NULL)
                delete[] l0BlueBandData[i];
        delete[] l0BlueBandData;
        l0BlueBandData = NULL;
        for (size_t i = 0; i < numCcdPixels; i++)
            if (l0RedBandData[i] != NULL)
                delete[] l0RedBandData[i];
        delete[] l0RedBandData;
        l0RedBandData = NULL;
        for (size_t i = 0; i < maxNumOfSwirPixels; i++)
            if (l0SwirBandData[i] != NULL)
                delete[] l0SwirBandData[i];
        delete[] l0SwirBandData;
        l0SwirBandData = NULL;

        delete[] blueCcdLineNums;
        delete[] redCcdLineNums;
        delete[] swirLineNums;

        delete[] swirDarkCalFrameTypeData[0];
        delete[] swirDarkCalFrameTypeData;

        delete[] l0SwirPixelFrameTypes;
        delete[] swirFrameTypesSciData;

        delete[] ancillaryData[0];
        delete[] ancillaryData;

        delete[] tlmdata[0];
        delete[] tlmdata;

        // update last seen time for the next packet to make sure the gap is
        // not too big
        lastAncillaryPktTime = ancillaryPktTime;

        if (ancillaryPktTime > granuleMaxTime)
            break;

    }  // while (!isEndFile) main loop end curly brace

    /////////////////// End Main Loop ///////////////////

    l1aFileManager.closeAllL1aFiles();

    // no l1a made, do not dump outlist buffer so no outlist is generated
    if (!l1aFileManager.filesWereGenerated()) {
        returnStatus = 120;
    }
    else {
        // dump all the file info into the file and close it
        l1aFileManager.dumpOutlistBuffer(outlist);
    }
    

    // Deallocate
    delete[] nextOciPacketBuffer[0];
    delete[] nextOciPacketBuffer;
    delete[] currOciPacketBuffer[0];
    delete[] currOciPacketBuffer;
    delete[] tempOciPacketBuffer[0];
    delete[] tempOciPacketBuffer;

    delete[] ccdLineError;
    delete[] sciPacketSequenceError;

    return returnStatus;
}
