#include <timeutils.h>
#include <string>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;
/* ---------------------------------------------------------------
    
    Given a string in the format:
        1. yyyy-mm-dd HH:MM:SSZ
        2. yyyy-mm-dd HH:MM:SS
        3. yyyy-mm-dd HH:MM:SS.fffZ
        4. yyyy-mm-dd HH:MM:SS.fff
        5. or use a T between date and time

    covert it into a double.
------------------------------------------------------------------*/
double udunits2unix(const string& str) {

    // split the string into 2, index 0 containing the date and index 1 with the time
    vector<string> datetime, date, time;
    boost::split(datetime, str, boost::is_any_of(" T"));
    if (datetime.size() != 2) {
        cout << "-E- Error: Udunits needs a space or 'T' between the date and time - " << str << endl;
        exit(1);
    }

    boost::split(date, datetime[0], boost::is_any_of("-"));
    if (date.size() != 3) {
        cout << "-E- Error: Udunits needs a '-' between the date fields - " << datetime[0] << endl;
        exit(1);
    }

    boost::split(time, datetime[1], boost::is_any_of(":"));
    if (time.size() != 3) {
        cout << "-E- Error: Udunits needs a ':' between the time fields - " << datetime[1] << endl;
        exit(1);
    }

    int year, month, day;
    int hours, mins;
    double secs;

    year = stoi(date[0]);
    month = stoi(date[1]);
    day = stoi(date[2]);

    hours = stoi(time[0]);
    mins = stoi(time[1]);
    
    // seconds may have a Z at the end.
    // remove the last char if there is a Z
    if (time[2][time[2].size()-1] == 'Z') {
        time[2] = time[2].substr(0, time[2].size()-1);
    }
    secs = stod(time[2]);

    double totalSecs =  hours * 3600.0 + mins * 60.0 + secs;
    return ymds2unix((int16_t)year, (int16_t)month, (int16_t)day, totalSecs);
}

// given a number as a string, if it's a single digit, pad 0 to
// the left
string padZero(const string& str) {
    string newStr = str;
    if (str.size() == 1)
        newStr = "0" + str;
    return newStr;
}

/* ---------------------------------------------------------------
    
    Given a unix time as a double, convert it into udunit string:
        1. yyyy-mm-dd HH:MM:SSZ
        2. yyyy-mm-dd HH:MM:SS

    @param time - unix time as a double
    @param hasZ - inclue the Z or not at the end
------------------------------------------------------------------*/
string unix2udunits(double unixtime, bool hasZ) {

    if (unixtime < 0) {
        cout << "-E- Error: time can't be less than 0 to convert to udunits in unix2udunits_cpp().\n" << endl;
        exit(1);
    }

    // convert into individual units
    int16_t year, month, day, hour, min;
    double sec;
    unix2ymdhms(unixtime, &year, &month, &day, &hour, &min, &sec);

    // format the month, day, and time to be 2 digits if it's 1
    string monthStr = padZero(std::to_string(month));
    string dayStr = padZero(std::to_string(day));
    string hourStr = padZero(std::to_string(hour));
    string minStr = padZero(std::to_string(min));
    string secStr = padZero(std::to_string((int)sec));

    // format the time
    ostringstream udunitStr;
    udunitStr << year << "-" 
                << monthStr << "-"
                << dayStr << " "
                << hourStr << ":"
                << minStr << ":"
                << secStr; 
    if (hasZ) {
        udunitStr << "Z";
    }
    return udunitStr.str();
}

// use the C++ function to do the work for the C version
extern "C" double udunits2unix(const char* str) {
    return udunits2unix((string) str);
}

// use the C++ function to do the work for the C version
extern "C" const char* unix2udunits_c(double unixtime, int hasZ) {
    bool wantZ;
    if(hasZ)
        wantZ = true;
    else
        wantZ = false;
    static string tmpStr = unix2udunits(unixtime, wantZ);
    return tmpStr.c_str();
}
