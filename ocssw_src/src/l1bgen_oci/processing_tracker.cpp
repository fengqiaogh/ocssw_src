
#include "processing_tracker.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>

/**
 * @brief Constructs a ProcessingTracker
 * @param startingLine The index upon which processing should start
 * @param endingLine The index upon which processing should end
 * @param interval The logging interval in percent of the delta between start line and end line. Should be
 * between 0 and 1
 */
oel::ProcessingTracker::ProcessingTracker(int32_t startingLine, int32_t endingLine, double interval) {
    if (0 < endingLine && endingLine < startingLine) {
        // Add one to the parameters so they match what the user input
        std::invalid_argument up("Invalid starting line: " + std::to_string(startingLine + 1) +
                                 " is not less than ending line of " + std::to_string(endingLine + 1));
        throw up;
    }

    if (interval < 0.0) {
        std::invalid_argument up("Invalid interval: Must be greater than 0");
        throw up;
    }

    this->endingLine = endingLine;

    if (startingLine < 0) {
        this->startingLine = 0;
        currScan = 0;
    } else {
        this->startingLine = startingLine;
        currScan = startingLine;
    }

    delta = endingLine - startingLine;

    this->interval = interval * 100.0;
    lastLoggedScan = 0;
}

oel::ProcessingTracker::ProcessingTracker() {
    startingLine = 0;
    currScan = startingLine;
    endingLine = -1;
    interval = 10.0;
    lastLoggedScan = 0;
}

bool oel::ProcessingTracker::shouldLogLine(size_t currLine) {
    percentComplete = ((currLine - startingLine) / delta) * 100;  // Precision between 0 & 1 isn't good

    if (currLine == startingLine) {
        return true;
    }

    if (currLine == endingLine - 1) {
        currLine++;
        percentComplete = ((currLine - startingLine) / delta) * 100;  // Precision between 0 & 1 isn't good
        return true;
    }

    size_t currentInterval = static_cast<int>(percentComplete / interval);

    if (currentInterval > lastLoggedScan) {
        lastLoggedScan = currentInterval;
        return true;
    } else {
        return false;
    }
}

std::string oel::ProcessingTracker::getLogString(size_t currScan) {
    return "    Scan " + std::to_string(currScan) + ": " + std::to_string((int)std::round(percentComplete)) +
           "\%";
}

void oel::ProcessingTracker::reset() {
    lastLoggedScan = 0;
    percentComplete = 0.0;
    currScan = 0;
}

void oel::ProcessingTracker::update() {
    if (shouldLogLine(currScan)) {
        std::cout << getLogString(currScan) << std::endl;
    }

    currScan++;
}