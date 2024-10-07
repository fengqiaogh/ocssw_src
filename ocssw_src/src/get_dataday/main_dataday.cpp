#include "get_dataday.h"

#include <iostream>
#include <string>

#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(0);
    }

    std::string inputFile = argv[1];
    if (inputFile == "-h" || inputFile == "--help") {
        printUsage(0);
    }
    // Open the file for read access
    int32_t day0, day1;
    set_verbosity(1);
    get_datadays(inputFile.c_str(), &day0, &day1);
    // const auto scan_time_earliest = boost::posix_time::ptime(
    //     boost::gregorian::date(1970, 1, 1) + boost::gregorian::days(day0));
    // const auto scan_time_latest = boost::posix_time::ptime(
    //     boost::gregorian::date(1970, 1, 1) + boost::gregorian::days(day1));
    // std::cout << "Day 0 is " << scan_time_earliest << std::endl;
    // std::cout << "Day 1 is " << scan_time_latest << std::endl;
    return 0;
}