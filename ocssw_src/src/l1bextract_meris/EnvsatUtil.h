/* 
 * File:   EnvsatUtil.h
 * Author: dshea
 *
 * Created on December 10, 2012, 12:33 PM
 */

#ifndef ENVSATUTIL_H
#define ENVSATUTIL_H

#include <string>
#include <stdint.h>

/** size of the UTC date fields */
static const unsigned int UTC_DATE_LENGTH = 27;
static const unsigned int LAT_LON_LENGTH = 11;

std::string& trim_right(std::string& s, const std::string& delimiters =
        " \n\r\t");
std::string& trim_left(std::string& s,
        const std::string& delimiters = " \n\r\t");
std::string& trim(std::string& s, const std::string& delimiters = " \n\r\t");

// generic set and get of string out of buffer
void setString(const std::string &str, char* buffer, unsigned int offset,
        int size);
void getString(const char* buffer, std::string &str, unsigned int offset,
        int size);

int getInt(const char* buffer, unsigned int offset, int size);
void setInt(int val, char* buffer, unsigned int offset, int size);
int64_t getInt64(const char* buffer, unsigned int offset, int size);
void setInt64(int64_t val, char* buffer, unsigned int offset, int size);

/** get a lat or lon value out of the character buffer
 * These assume the lat and lon are written as ints * 1e6
 */
double getLatLon(const char* buffer, unsigned int offset);
void setLatLon(double val, char* buffer, unsigned int offset);

int getRawInt16(const char* buffer, unsigned int offset);
void setRawInt16(int val, char* buffer, unsigned int offset);
int getRawInt32(const char* buffer, unsigned int offset);
void setRawInt32(int val, char* buffer, unsigned int offset);

double getMJD(const char* buffer, unsigned int offset);
void setMJD(int day, int sec, int microsec, char* buffer, unsigned int offset);

double merisTime2unix(const std::string& timeStr);
const std::string& unix2merisTime(double unixTime);

double envsatInterp(double x1, double x2, double y1, double y2, double xin);

#endif  /* ENVSATUTIL_H */
