#include "l1afile_oci_manager.hpp"
using namespace std;


L1aFileManager::L1aFileManager() {}

L1aFileManager::~L1aFileManager() {}

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