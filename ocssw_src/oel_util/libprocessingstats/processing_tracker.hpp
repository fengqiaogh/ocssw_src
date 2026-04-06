
#ifndef OEL_PROCESSING_TRACKER
#define OEL_PROCESSING_TRACKER

#include <string>
#include <stdexcept>

namespace oel {
class ProcessingTracker {
   public:
    ProcessingTracker(int32_t startingScan, int32_t endingScan, double loggingRate);
    ProcessingTracker();

    /**
     * @brief Updates state variables and performs operations as necessary
     * @note May print to stdout or stderr
     */
    void update();

    inline size_t getStartingLine() const {
        return startingScan;
    }

    inline size_t getEndingLine() const {
        return endingScan;
    }

    inline double getInterval() const {
        return loggingRate;
    }

    inline size_t getCurrScan() const {
        return currScan;
    }

    inline size_t getPercentComplete() const {
        return percentComplete;
    }
    
    /**
     * @brief Allows the user to set the current scan manually.
     * @param newCurrScan The desired current scan
     * @throws std::invalid argument when passed scan number is outside the range defined by the current
     * starting line and ending line
     */
    inline void setCurrScan(size_t newCurrScan) {

        // If user passed a negative number, newCurrScan will likely be of a higher value than endingLine
        if (newCurrScan < startingScan || endingScan < newCurrScan) {
            throw std::invalid_argument("Current scan must be between starting and ending line");
        }

        currScan = newCurrScan;
    }

    /**
     * @brief Sets the logging loggingRate
     * @param newInterval The desired loggingRate as a percentage 
     * @throws std::invalid_argument if passed loggingRate is negative
     * @throws std::invalid_argument if passed loggingRate is over 1.0
     */
    inline void setInterval(double newInterval) {
        if (newInterval < 0.0) {
            throw std::invalid_argument("Interval must be positive");
        }

        if (newInterval > 1.0) {
            throw std::invalid_argument("Interval must be less than 1.0");
        }

        loggingRate = newInterval * 100.0;
    }

    /**
     * @brief Resets the state to 0
     */
    void reset();

   private:
    size_t currScan;
    size_t startingScan;  // starting line (first line = 0)
    size_t endingScan;    // ending line inclusive (first line = 0)
    double delta;  // Difference between endingLine and startingLine. Double so logging calc is in floating
                   // point space
    double loggingRate;
    size_t lastLoggedScan;
    double percentComplete;

    bool shouldLogLine();        // Determines whether the current line should be logged
    std::string getLogString();  // Constructs and returns a string to be used for logging
};
}  // namespace oel

#endif