/* 
 * File:   MerisSPH.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "MerisSPH.h"

#include "EnvsatUtil.h"

#include <stdio.h>

using namespace std;

MerisSPH::MerisSPH(const MerisSPH& orig) :
EnvsatSPH(orig) {

}

MerisSPH::MerisSPH(EnvsatFile* file, int size, int numDSDs, int DSDSize) :
EnvsatSPH(file, size, numDSDs, DSDSize) {

}

//MerisSPH::~MerisSPH() {
//
//}

string& MerisSPH::getDescriptor() {
    static string name;
    getString(buffer, name, DESCRIPTOR_OFFSET, DESCRIPTOR_LENGTH);
    return name;
}

double MerisSPH::getFirstLineTime() {
    string str;
    getString(buffer, str, FIRST_LINE_TIME_OFFSET, UTC_DATE_LENGTH);
    return merisTime2unix(str);
}

void MerisSPH::setFirstLineTime(double unixTime) {
    setString(unix2merisTime(unixTime), buffer, FIRST_LINE_TIME_OFFSET,
            UTC_DATE_LENGTH);
}

double MerisSPH::getLastLineTime() {
    string str;
    getString(buffer, str, LAST_LINE_TIME_OFFSET, UTC_DATE_LENGTH);
    return merisTime2unix(str);
}

void MerisSPH::setLastLineTime(double unixTime) {
    setString(unix2merisTime(unixTime), buffer, LAST_LINE_TIME_OFFSET,
            UTC_DATE_LENGTH);
}

double MerisSPH::getFirstFirstLat() {
    return getLatLon(buffer, FIRST_FIRST_LAT_OFFSET);
}

void MerisSPH::setFirstFirstLat(double lat) {
    setLatLon(lat, buffer, FIRST_FIRST_LAT_OFFSET);
}

double MerisSPH::getFirstFirstLon() {
    return getLatLon(buffer, FIRST_FIRST_LON_OFFSET);
}

void MerisSPH::setFirstFirstLon(double lon) {
    setLatLon(lon, buffer, FIRST_FIRST_LON_OFFSET);
}

double MerisSPH::getFirstMidLat() {
    return getLatLon(buffer, FIRST_MID_LAT_OFFSET);
}

void MerisSPH::setFirstMidLat(double lat) {
    setLatLon(lat, buffer, FIRST_MID_LAT_OFFSET);
}

double MerisSPH::getFirstMidLon() {
    return getLatLon(buffer, FIRST_MID_LON_OFFSET);
}

void MerisSPH::setFirstMidLon(double lon) {
    setLatLon(lon, buffer, FIRST_MID_LON_OFFSET);
}

double MerisSPH::getFirstLastLat() {
    return getLatLon(buffer, FIRST_LAST_LAT_OFFSET);
}

void MerisSPH::setFirstLastLat(double lat) {
    setLatLon(lat, buffer, FIRST_LAST_LAT_OFFSET);
}

double MerisSPH::getFirstLastLon() {
    return getLatLon(buffer, FIRST_LAST_LON_OFFSET);
}

void MerisSPH::setFirstLastLon(double lon) {
    setLatLon(lon, buffer, FIRST_LAST_LON_OFFSET);
}

double MerisSPH::getLastFirstLat() {
    return getLatLon(buffer, LAST_FIRST_LAT_OFFSET);
}

void MerisSPH::setLastFirstLat(double lat) {
    setLatLon(lat, buffer, LAST_FIRST_LAT_OFFSET);
}

double MerisSPH::getLastFirstLon() {
    return getLatLon(buffer, LAST_FIRST_LON_OFFSET);
}

void MerisSPH::setLastFirstLon(double lon) {
    setLatLon(lon, buffer, LAST_FIRST_LON_OFFSET);
}

double MerisSPH::getLastMidLat() {
    return getLatLon(buffer, LAST_MID_LAT_OFFSET);
}

void MerisSPH::setLastMidLat(double lat) {
    setLatLon(lat, buffer, LAST_MID_LAT_OFFSET);
}

double MerisSPH::getLastMidLon() {
    return getLatLon(buffer, LAST_MID_LON_OFFSET);
}

void MerisSPH::setLastMidLon(double lon) {
    setLatLon(lon, buffer, LAST_MID_LON_OFFSET);
}

double MerisSPH::getLastLastLat() {
    return getLatLon(buffer, LAST_LAST_LAT_OFFSET);
}

void MerisSPH::setLastLastLat(double lat) {
    setLatLon(lat, buffer, LAST_LAST_LAT_OFFSET);
}

double MerisSPH::getLastLastLon() {
    return getLatLon(buffer, LAST_LAST_LON_OFFSET);
}

void MerisSPH::setLastLastLon(double lon) {
    setLatLon(lon, buffer, LAST_LAST_LON_OFFSET);
}

int MerisSPH::getSamplesPerLine() {
    return getInt(buffer, SAMPLES_PER_LINE_OFFSET, SAMPLES_PER_LINE_LENGTH);
}

int MerisSPH::getLinesPerTiepoint() {
    return getInt(buffer, LINES_PER_TIEPOINT_OFFSET, LINES_PER_TIEPOINT_LENGTH);
}

int MerisSPH::getSamplesPerTiepoint() {
    return getInt(buffer, SAMPLES_PER_TIEPOINT_OFFSET,
            SAMPLES_PER_TIEPOINT_LENGTH);
}

const string& MerisSPH::getQualityName() {
    static string result("Quality ADS");
    return result;
}

const string& MerisSPH::getTiepointName() {
    static string result("Tie points ADS");
    return result;
}

void MerisSPH::print() {
    EnvsatSPH::print();
    printf("\nMeris SPH\n");
    printf("firstLineTime  = %s\n", unix2merisTime(getFirstLineTime()).c_str());
    printf("lastLineTime   = %s\n", unix2merisTime(getLastLineTime()).c_str());

    printf("firstFirstLat  = %f\n", getFirstFirstLat());
    printf("firstFirstLon  = %f\n", getFirstFirstLon());
    printf("firstMidLat    = %f\n", getFirstMidLat());
    printf("firstMidLon    = %f\n", getFirstMidLon());
    printf("firstLastLat   = %f\n", getFirstLastLat());
    printf("firstLastLon   = %f\n", getFirstLastLon());

    printf("lastFirstLat   = %f\n", getLastFirstLat());
    printf("lastFirstLon   = %f\n", getLastFirstLon());
    printf("lastMidLat     = %f\n", getLastMidLat());
    printf("lastMidLon     = %f\n", getLastMidLon());
    printf("lastLastLat    = %f\n", getLastLastLat());
    printf("lastLastLon    = %f\n", getLastLastLon());

    printf("samplesPerLine = %d\n", getSamplesPerLine());
    printf("linesPerTiepoint = %d\n", getLinesPerTiepoint());
    printf("samplesPerTiepoint = %d\n", getSamplesPerTiepoint());
}
