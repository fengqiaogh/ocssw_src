#include "l1afile_oci_manager.hpp"
using namespace std;


// initialize previous packet time of size 14 (number of DATA_TYPES for l1agen_coi)
// set to -1 to indiciate it has not been seen yet
L1aFileManager::L1aFileManager() : prevAncillaryPktTimes(14, -1.0) {}

L1aFileManager::~L1aFileManager() {}


// update time for dataType
void L1aFileManager::updatePrevAncillaryPktTimeFor(short dataType, double time) {
    prevAncillaryPktTimes[dataType] = time;
}

// request for the previous ancillary packet time with dataType
double L1aFileManager::getPrevAncillaryPktTimeFor(short dataType) {
    
    // this is called when the dataType has been seen and has a time set
    // if there's no time set, something is wrong. Issue a warning 
    if (prevAncillaryPktTimes[dataType] == -1) {
        cout << "--WARNING-- using previous ancillary packet time, but it was never updated at least once" << endl;
    }
    return prevAncillaryPktTimes[dataType];
}

// when l1a for said datatype has not been created, make one and return a reference to it
L1aFile* L1aFileManager::createL1aOutputFile(short dataType) {
    if (fileCabinet.count(dataType) == 0) {
        fileCabinet.emplace(dataType, L1aFile());              // make L1A file object for data type
        fileDimShapes.emplace(dataType, DimensionShape());     // make corresponding shape tracker for it
    }
    return &fileCabinet[dataType];
}


// checks if the dataType already exists
bool L1aFileManager::contains(short dataType) {
    if (fileCabinet.find(dataType) != fileCabinet.end()) {
        return true;
    }
    return false;
}

// update the the number_of_scans dimension so the next write knows where to start
bool L1aFileManager::incrementNumScansShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        fileDimShapes[dataType].incrementNumScansShape(incrementBy);
        return true;
    }
    return false;
}

// update the the number_of_mce_scans dimension so the next write knows where to start
bool L1aFileManager::incrementNumMceScanShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].numMceScans;
        fileDimShapes[dataType].numMceScans = currSize + incrementBy;
        return true;
    }
    return false;
}


// update the the number_of_sca_scans dimension so the next write knows where to start
bool L1aFileManager::incrementNumScaScanShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].numScaScans;
        fileDimShapes[dataType].numScaScans = currSize + incrementBy;
        return true;
    }
    return false;
}

// update the the att_records dimension so the next write knows where to start
bool L1aFileManager::incrementAttRecordsShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].attRecords;
        fileDimShapes[dataType].attRecords = currSize + incrementBy;
        return true;
    }
    return false;
}

// update the the orb_records dimension so the next write knows where to start
bool L1aFileManager::incrementOrbRecordsShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].orbRecords;
        fileDimShapes[dataType].orbRecords = currSize + incrementBy;
        return true;
    }
    return false;
}

// update the the tlm_packets dimension so the next write knows where to start
bool L1aFileManager::incrementTlmPacketsShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].tlmPackets;
        fileDimShapes[dataType].tlmPackets = currSize + incrementBy;
        return true;
    }
    return false;
}

// update the the tilt_sample dimension so the next write knows where to start
bool L1aFileManager::incrementTiltSampleShape(short dataType, size_t incrementBy) {
    if (contains(dataType)) {
        size_t currSize = fileDimShapes[dataType].tiltSamples;
        fileDimShapes[dataType].tiltSamples = currSize + incrementBy;
        return true;
    }
    return false;
}


// close all l1afiles
void L1aFileManager::closeAllL1aFiles() {
    for (auto file: fileCabinet) {
        file.second.close();
        file.second.freeFile();
    }
    fileCabinet.clear();
    fileDimShapes.clear();
}

// close and remove file for dataType from the map and also its cooresponding DimShape
void L1aFileManager::closeAndRemoveFile(short dataType) {
    fileCabinet[dataType].close();
    fileCabinet.erase(dataType);
    fileDimShapes.erase(dataType);
}

DimensionShape* L1aFileManager::getCurrL1aFileDimShape(short dataType) {
    if (fileCabinet.find(dataType) != fileCabinet.end()) {
        return &fileDimShapes[dataType];
    }
    return nullptr;
}

// for datatype, find the l1afile for said
L1aFile* L1aFileManager::getL1aFile(short dataType) {
    if (fileCabinet.find(dataType) != fileCabinet.end()) {
        return &fileCabinet[dataType];
    }
    return nullptr;
}

//////////////////////////////////////////////
//////////// outlist functions////////////////
/////////////////////////////////////////////

// add file into the 
void L1aFileManager::addFileToOutlistBuffer(string fileName, string startTime, string endTime,
                                        int completeFlag, string dataType, short dataTypeId) {


    // lunar files, format the string as it and keep it in a vector for writing later
    if (dataType == "LUN") {
        string outString = fileName + " " + startTime + " " + endTime + " " \
                             + to_string(completeFlag) + " " + dataType + "\n";
        noMergeOutlistStrings.push_back(outString);
        return;
    }

    // not lunar file and it is already in the map merge meta data
    if (outlistData.find(dataTypeId) != outlistData.end()) {

        // end time will get extended
        outlistData[dataTypeId]["end_time"] = endTime;

        // merging complete and incomplete file results in incomplete
        if (completeFlag == 1) {
            outlistData[dataTypeId]["complete_flag"] = "1";
        }
        return;
    }

    // otherwise, new data type to be added 
    outlistData[dataTypeId] = {
        {"data_type", dataType},
        {"file_name", fileName},
        {"start_time", startTime},
        {"end_time", endTime},
        {"complete_flag", to_string(completeFlag)},
    };
}


// take all the data from the buffer and write it to the file
void L1aFileManager::dumpOutlistBuffer(string outlist) {

    // open outlist file to write into
    ofstream outFile;
    if (outlist.compare("") != 0)
        outFile.open(outlist.c_str());

    // no outlist provided, dont write anything and just close
    else {
        outFile.close();
        return; 
    }

    // empty case, no data or no l1a files made, keep outlist empty
    if (outlistData.size() == 0 && noMergeOutlistStrings.size() == 0) {
        outFile.close();
        return;
    }

    // check if there's lunar files or files that were closed early due to 
    // file gap being too large
    if (noMergeOutlistStrings.size() > 0) {
        for (string str: noMergeOutlistStrings) {
            // dont need new line bc string should already end in one when added to
            // the list 
            outFile << str; 
        }
    }
    
    // check all the other files that can be merged
    if (outlistData.size() > 0) {
        for (const auto& item : outlistData) {
            short dataTypeId = item.first;
            string dataType = outlistData[dataTypeId]["data_type"];
            string fileName = outlistData[dataTypeId]["file_name"];
            string startTime = outlistData[dataTypeId]["start_time"];
            string endTime = outlistData[dataTypeId]["end_time"];
            string completeFlag = outlistData[dataTypeId]["complete_flag"];

            // check navigation coverage if Science Fiile
            if (dataTypeId == SCIENCE_FILE && !navTimeFrame.isGoodNavigationCoverage()) {
                completeFlag = "-1";
            }

            outFile << fileName << " " << startTime << " " << endTime << " " \
                << completeFlag << " " << dataType << "\n";
        }
    }
    outFile.close();

}

// // update last file type seen
// void L1aFileManager::updateLastDataTypeSeen(short dataType) {
//     lastDataTypeSeen = dataType;
// }

// add the previous file to the noMergeOutlistString
void L1aFileManager::processOutlistFor(short dataType) {
    string type = outlistData[dataType]["data_type"];
    string fileName = outlistData[dataType]["file_name"];
    string startTime = outlistData[dataType]["start_time"];
    string endTime = outlistData[dataType]["end_time"];
    string completeFlag = outlistData[dataType]["complete_flag"];

    // bad navigation coverage for a science file, replace complete flag
    if (dataType == SCIENCE_FILE && !navTimeFrame.isGoodNavigationCoverage()) {
        completeFlag = "-1";
    }

    string outString = fileName + " " + startTime + " " + endTime + " " \
                             + completeFlag + " " + type + "\n";

    noMergeOutlistStrings.push_back(outString);

    // remove it from the map 
    outlistData.erase(dataType);
}

bool L1aFileManager::filesWereGenerated() {
    return outlistData.size() > 0 || noMergeOutlistStrings.size() > 0;
}


//////////////////////////////////////////////
///////// navigation coverage functions//////
/////////////////////////////////////////////

NavigationTimeFrame* L1aFileManager::getNavTimeFrameRef() {
    return &navTimeFrame;
}




