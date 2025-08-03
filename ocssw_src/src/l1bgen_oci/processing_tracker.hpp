
#ifndef L1B_SCAN_RANGE_H
#define L1B_SCAN_RANGE_H

#include <string>

namespace oel {
class ProcessingTracker {
   public:
    ProcessingTracker(int32_t startingLine, int32_t endingLine, double interval);
    ProcessingTracker();

    /**
     * @brief Updates state variables and performs operations as necessary
     * @note May print to stdout or stderr
     */
    void update();

    inline size_t getStartingLine() {
        return startingLine;
    }

    inline size_t getEndingLine() {
        return endingLine;
    }

    /**
     * @brief Resets the state to 0
     */
    void reset();

   private:
    size_t currScan;
    size_t startingLine;  // starting line (first line = 0)
    size_t endingLine;    // ending line inclusive (first line = 0)
    double delta;  // Difference between endingLine and startingLine. Double so logging calc is in floating
                   // point space
    double interval;
    size_t lastLoggedScan;
    double percentComplete;

    bool shouldLogLine(size_t currLine);        // Determines whether the current line should be logged
    std::string getLogString(size_t currScan);  // Constructs and returns a string to be used for logging
};
}  // namespace oel

#endif