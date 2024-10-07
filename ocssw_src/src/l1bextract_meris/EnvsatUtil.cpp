/* 
 * File:   EnvsatUtil.cpp
 * Author: dshea
 * 
 * Created on November 28, 2012, 1:11 PM
 */

#include "EnvsatUtil.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <timeutils.h>

using namespace std;

const char* monthArray[12] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
    "AUG", "SEP", "OCT", "NOV", "DEC"};

string& trim_right(string& s, const string& delimiters) {
    return s.erase(s.find_last_not_of(delimiters) + 1);
}

string& trim_left(string& s, const string& delimiters) {
    return s.erase(0, s.find_first_not_of(delimiters));
}

string& trim(string& s, const string& delimiters) {
    return trim_left(trim_right(s, delimiters), delimiters);
}

void setString(const string &str, char* buffer, unsigned int offset, int size) {
    int i;
    for (i = 0; i < (int) str.size(); i++) {
        if (i >= size)
            break;
        buffer[offset + i] = str[i];
    }
    for (; i < size; i++) {
        buffer[offset + i] = ' ';
    }
}

void getString(const char* buffer, string &str, unsigned int offset, int size) {
    str.assign(buffer + offset, size);
    trim(str);
}

int getInt(const char* buffer, unsigned int offset, int size) {
    int val;
    char formatStr[32];

    sprintf(formatStr, "%%%dd", size);
    sscanf(buffer + offset, formatStr, &val);
    return val;
}

void setInt(int val, char* buffer, unsigned int offset, int size) {
    char tmpBuff[size + 32];

    sprintf(tmpBuff, "%+0*d", size, val);
    memcpy(buffer + offset, tmpBuff, size);
}

int64_t getInt64(const char* buffer, unsigned int offset, int size) {
    long long val;
    char formatStr[32];

    sprintf(formatStr, "%%%dlld", size);
    sscanf(buffer + offset, formatStr, &val);
    return val;
}

void setInt64(int64_t val, char* buffer, unsigned int offset, int size) {
    char tmpBuff[size + 32];

    sprintf(tmpBuff, "%+0*lld", size, (long long) val);
    memcpy(buffer + offset, tmpBuff, size);
}

double getLatLon(const char* buffer, unsigned int offset) {
    int64_t val = getInt64(buffer, offset, LAT_LON_LENGTH);
    return ((double) val) * 1.0e-6;
}

void setLatLon(double val, char* buffer, unsigned int offset) {
    int64_t i = round(val * 1.0e6);
    setInt64(i, buffer, offset, LAT_LON_LENGTH);
}

int getRawInt16(const char* buffer, unsigned int offset) {
    return (((buffer[offset + 0] & 0xff) << 8) | (buffer[offset + 1] & 0xff));
}

void setRawInt16(int val, char* buffer, unsigned int offset) {
    ((unsigned char*) buffer)[offset + 0] = (val >> 8) & 0xff;
    ((unsigned char*) buffer)[offset + 1] = val & 0xff;
}

int getRawInt32(const char* buffer, unsigned int offset) {
    return (((buffer[offset + 0] & 0xff) << 24)
            | ((buffer[offset + 1] & 0xff) << 16)
            | ((buffer[offset + 2] & 0xff) << 8) | (buffer[offset + 3] & 0xff));
}

void setRawInt32(int val, char* buffer, unsigned int offset) {
    ((unsigned char*) buffer)[offset + 0] = (val >> 24) & 0xff;
    ((unsigned char*) buffer)[offset + 1] = (val >> 16) & 0xff;
    ((unsigned char*) buffer)[offset + 2] = (val >> 8) & 0xff;
    ((unsigned char*) buffer)[offset + 3] = val & 0xff;
}

double getMJD(const char* buffer, unsigned int offset) {
    int day, sec, microsec;
    struct tm trec;
    time_t secSince;

    day = getRawInt32(buffer, offset);
    sec = getRawInt32(buffer, offset + 4);
    microsec = getRawInt32(buffer, offset + 8);

    trec.tm_year = 2000 - 1900;
    trec.tm_mon = 0;
    trec.tm_mday = day + 1;
    trec.tm_hour = 0;
    trec.tm_min = 0;
    trec.tm_sec = 0;
    trec.tm_isdst = 0;
    secSince = mktime(&trec) - gmt_offset();

    /*                                                  */
    /* Now we total the seconds                         */
    /*                                                  */
    return secSince + sec + microsec * 1.0e-6;
}

void setMJD(int day, int sec, int microsec, char* buffer, unsigned int offset) {
    setRawInt32(day, buffer, offset);
    setRawInt32(sec, buffer, offset + 4);
    setRawInt32(microsec, buffer, offset + 8);
}

double merisTime2unix(const string& timeStr) {
    // 02-AUG-2002 18:45:14.058298
    int16_t year, month, day, hour, min;
    char monthStr[5];
    double sec;

    sscanf(timeStr.c_str(), "%2hd-%3s-%4hd %2hd:%2hd:%lf", &day, monthStr,
            &year, &hour, &min, &sec);
    for (month = 0; month < 12; month++) {
        if (strcmp(monthStr, monthArray[month]) == 0)
            break;
    }
    sec += hour * 3600.0 + min * 60.0;
    return ymds2unix(year, month + 1, day, sec);
}

const string& unix2merisTime(double unixTime) {
    // 02-AUG-2002 18:45:14.058298
    int16_t year, month, day, hour, min;
    double sec;
    char str[32];
    static string timeStr;

    unix2ymdhms(unixTime, &year, &month, &day, &hour, &min, &sec);
    sprintf(str, "%02d-%s-%04d %02d:%02d:%09.6f", day, monthArray[month - 1],
            year, hour, min, sec);
    timeStr.assign(str);
    return timeStr;
}

double envsatInterp(double x1, double x2, double y1, double y2, double xin) {
    double m = (y2 - y1) / (x2 - x1);
    double x = xin - x1;
    return m * x + y1;
}

