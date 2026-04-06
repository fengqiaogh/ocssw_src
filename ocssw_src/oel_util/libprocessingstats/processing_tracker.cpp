
#include "processing_tracker.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

/**
 * @brief Constructs a ProcessingTracker
 * @param startingScan The index upon which processing should start
 * @param endingScan The index upon which processing should end
 * @param loggingRate The logging rate in percent of the delta between start line and end line. Should be
 * between 0 and 100
 */
oel::ProcessingTracker::ProcessingTracker(int32_t startingScan, int32_t endingScan, double loggingRate) {
    if (endingScan < 0 || endingScan < startingScan) {
        // Add one to the parameters so they match what the user input
        std::invalid_argument up("Invalid starting line: " + std::to_string(startingScan + 1) +
                                 " is not less than ending line of " + std::to_string(endingScan + 1));
        throw up;
    }

    if (loggingRate < 0.0) {
        std::invalid_argument up("Invalid logging rate: Must be greater than 0");
        throw up;
    }

    this->endingScan = endingScan;

    this->startingScan = std::max(0, startingScan);
    currScan = std::max(0, startingScan);

    delta = endingScan - startingScan;

    this->loggingRate = loggingRate * 100.0;
    lastLoggedScan = 0;
}

/**
 * @brief Initializes an invalid object. The user is expected to use the setters after using this constructor
 */
oel::ProcessingTracker::ProcessingTracker() {
    startingScan = 0;
    currScan = startingScan;
    endingScan = -1;
    loggingRate = 10.0;
    lastLoggedScan = 0;
    delta = 0;
    percentComplete = 0.0;
}

bool oel::ProcessingTracker::shouldLogLine() {

    if (delta == 0 || endingScan < 0 || endingScan < startingScan) {
        return false;
    }

    percentComplete = ((currScan - startingScan) / delta) * 100;  // Precision between 0 & 1 isn't good

    if (currScan == startingScan) {
        return true;
    }

    if (currScan == endingScan - 1) {
        currScan++;
        percentComplete = ((currScan - startingScan) / delta) * 100;  // Precision between 0 & 1 isn't good
        return true;
    }

    size_t currentInterval = static_cast<int>(percentComplete / loggingRate);

    if (currentInterval > lastLoggedScan) {
        lastLoggedScan = currentInterval;
        return true;
    } else {
        return false;
    }
}

std::string oel::ProcessingTracker::getLogString() {

    std::string currScanString = std::to_string(currScan);
    std::string percentString = std::to_string((int)std::round(percentComplete));

    return "    Scan " + currScanString + ": " + percentString + "\%";
}

void oel::ProcessingTracker::reset() {
    lastLoggedScan = 0;
    percentComplete = 0.0;
    currScan = 0;
}

void oel::ProcessingTracker::update() {
    if (shouldLogLine()) {
        std::cout << getLogString() << std::endl;
    }

    currScan++;
}