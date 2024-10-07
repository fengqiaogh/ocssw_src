#include <genutils.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

bool is_digits(const std::string &str)
{
    return std::all_of(str.begin(), str.end(), ::isdigit);
}

double string2resolution(std::string resolutionStr) {
    double resolution = BAD_FLT;
    boost::trim(resolutionStr);
    boost::to_lower(resolutionStr);
    resolutionStr.erase(remove_if(resolutionStr.begin(), resolutionStr.end(), ::isspace), resolutionStr.end());

    if (resolutionStr.compare("90km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 432.0;
    else if (resolutionStr.compare("36km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if (resolutionStr.compare("18km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 2160.0;
    else if (resolutionStr.compare("9km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4320.0;
    else if (resolutionStr.compare("4km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if (resolutionStr.compare("2km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 17280.0;
    else if (resolutionStr.compare("1km") == 0)
        resolution = EARTH_CIRCUMFERENCE / 34560.0;
    else if (resolutionStr.compare("hkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 69120.0;
    else if (resolutionStr.compare("qkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 138240.0;
    else if (resolutionStr.compare("hqkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 276480.0;
    else if (resolutionStr.compare("hhkm") == 0)
        resolution = EARTH_CIRCUMFERENCE / 552960.0;
    else if (resolutionStr.compare("smi") == 0)
        resolution = EARTH_CIRCUMFERENCE / 4096.0;
    else if (resolutionStr.compare("smi4") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8192.0;
    else if (resolutionStr.compare("land") == 0)
        resolution = EARTH_CIRCUMFERENCE / 8640.0;
    else if (resolutionStr.compare("thirddeg") == 0)
        resolution = EARTH_CIRCUMFERENCE / 1080.0;
    else if (boost::ends_with(resolutionStr, "km")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 2);
        resolution = atof(val.c_str()) * 1000.0;
    }  else if (boost::ends_with(resolutionStr, "m")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 1);
        resolution = atof(val.c_str());
    } else if (boost::ends_with(resolutionStr, "deg")) {
        std::string val = resolutionStr.substr(0, resolutionStr.length() - 3);
        resolution = atof(val.c_str()) / 360.0 * EARTH_CIRCUMFERENCE;
    } else if (is_digits(resolutionStr)) {
        resolution = atof(resolutionStr.c_str());
    } else {
        resolution = BAD_FLT;
    }
    return resolution;
}

std::string resolution2string(double resolution) {
    std::string resolutionStr = "Unknown";
    if (resolution > 0) {
        std::string geoUnits = " m";
        if (resolution > 1000.0) {
            resolution = resolution / 1000.0;
            geoUnits = " km";
        } 
        resolutionStr = std::to_string(resolution) + geoUnits;
    }
    return resolutionStr;
}

extern "C" double str2resolution(char const * resolutionStr) {
    std::string str(resolutionStr);
    return string2resolution(str);
}

extern "C" const char* resolution2str(double resolution) {
    static std::string stringRes = resolution2string(resolution);
    return stringRes.c_str();
}

extern "C" double resolution2degrees(double resolution) {
    return resolution / EARTH_CIRCUMFERENCE * 360.0;
}

extern "C" double degrees2resolution(double degrees) {
    return degrees / 360.0 * EARTH_CIRCUMFERENCE;
}

// this is called resolve since the string is different that the resolution string above.
// This is the input string (resolve) to l2bin
void resolve2binRows(std::string resolve, int32_t &nrows, double &resolution) {
    if (resolve == "1D") {
        nrows = 180;
        resolution = str2resolution("1deg");
    } else if (resolve == "HD") {
        nrows = 360;
        resolution = str2resolution("0.5deg");
    } else if (resolve == "QD") {
        nrows = 720;
        resolution = str2resolution("0.25deg");
    } else if (resolve == "36") {
        nrows = 2160 / 4;
        resolution = str2resolution("36km");
    } else if (resolve == "18") {
        nrows = 2160 / 2;
        resolution = str2resolution("18km");
    } else if (resolve == "9") {
        nrows = 2160;
        resolution = str2resolution("9km");
    } else if (resolve == "4") {
        nrows = 2160 * 2;
        resolution = str2resolution("4km");
    } else if (resolve == "2") {
        nrows = 2160 * 4;
        resolution = str2resolution("2km");
    } else if (resolve == "1") {
        nrows = 2160 * 8;
        resolution = str2resolution("1km");
    } else if (resolve == "H") {
        nrows = 2160 * 16;
        resolution = str2resolution("hkm");
    } else if (resolve == "Q") {
        nrows = 2160 * 32;
        resolution = str2resolution("qkm");
    } else if (resolve == "HQ") {
        nrows = 2160 * 80;
        resolution = str2resolution("hqkm");
    } else if (resolve == "HH") {
        nrows = 2160 * 160;
        resolution = str2resolution("hhkm");
    } else {
        nrows = -1;
        resolution = -1;
    }
}

extern "C" int32_t resolve2binRows(const char *resolve) {
    int32_t nrows;
    double resolution;
    resolve2binRows(resolve, nrows, resolution);
    return nrows;
}

extern "C" double resolve2resolution(const char *resolve) {
    int32_t nrows;
    double resolution;
    resolve2binRows(resolve, nrows, resolution);
    return resolution;
}


