/**************************************************************************
 *
 * NAME: VcstLogger.h
 *
 * DESCRIPTION: provides logging capabilities for the runtime environment.
 * Debug messages can be set to off, low, medium, or high and can be sent
 * to the console or to a log file.
 *
 * Adapted directly from the IDPS file ProCmnLogger.h published
 * by Raytheon Company.
 *
 **************************************************************************/

#ifndef _VcstLogger_H_
#define _VcstLogger_H_

#include <fstream>
#include <VcstTime.h>

/**
 * This class provides logging capabilities for the ADL environment.  Debug
 * messages can be set to off, low, medium, or high and can be sent to the
 * console or to a log file.
 */
class VcstLogger {
public:

    /**
     * Logger detail levels.
     */
    enum Level {
        DBG_OFF = 0, // Disable logging
        DBG_HIGH, // Most basic details
        DBG_MED, // Debugging level
        DBG_LOW // Low-level detail
    };

    /**
     * Logger destination.
     */
    enum LoggingDest {
        STDOUT = 0, LOGFILE
    };

    /**
     * Destructor
     */
    virtual ~VcstLogger();

    /**
     * Initialize my logger.  Sets the logging detail level retrieved
     * from the algorithm execution XML file.
     *
     * @param level An enum value representing the log level.
     * @param dest An enum value representing the message destination
     * (STDOUT or a log file).
     * @param processName The name of the executable that is running.
     * @param path The path where the log file will be written.
     */
    void init(VcstLogger::Level level, VcstLogger::LoggingDest dest,
            std::string processName, std::string path = ".");

    /**
     * Initialize my logger.  Sets the logging detail level retrieved
     * from the algorithm execution XML file.
     *
     * This init method calls the version of init which takes enum values to
     * do the actual initialization.
     *
     * @param level A string value representing the log level.
     * @param dest A string value representing the message destination
     * (STDOUT or a log file).
     * @param processName The name of the executable that is running.
     * @param path The path where the log file will be written.
     */
    void init(std::string level, std::string dest, std::string processName,
            std::string path = ".");

    /**
     * Get visibility to the reference to the INF logging
     * service implementation.
     */
    static VcstLogger& getLogger();

    /**
     * Setter for our logger instance.  My implementation of the
     * VcstLoggerGCIF interface method.
     *
     * @param aLogger A new logger instance.  My garbage collector
     * deletes this memory.
     * @retval The old logger instance.  My caller is responsible for
     * deleting this memory.
     */
    static VcstLogger* setLogger(VcstLogger* aLogger);

    /**
     * Translate a return error number into a string.
     *
     * @param status, one of PRO_FAIL, PRO_SUCCESS, ... (see Typedefs.h>
     * @retval string representation of the error code.
     */
    static const char *errorAsStr(int status);

    /**
     * Determine if a given error level is an error case or not.
     *
     * @param status, one of PRO_FAIL, PRO_SUCCESS, ... (see Typedefs.h>
     * @retval false if the error was PRO_SUCCESS, else return false.
     */
    static bool wasError(int status);

    /**
     * Returns whether or not DBG_HIGH is enabled.
     *
     * @retval true if DBG_HIGH is enabled.
     */
    bool isDebugHighEnabled();

    /**
     * Returns whether or not DBG_MED is enabled.
     *
     * @retval true if DBG_MED is enabled.
     */
    bool isDebugMedEnabled();

    /**
     * Returns whether or not DBG_LOW is enabled.
     *
     * @retval true if DBG_LOW is enabled.
     */
    bool isDebugLowEnabled();

    /**
     * Logs a HIGH debug message.
     *
     * @param msg The debug message to send.
     * @param file The name of the file this message was sent from.
     * @param line The line number this message was sent from.
     */
    void debugHigh(std::string msg, const char* file, const int line);

    /**
     * Logs a MED debug message.
     *
     * @param msg The debug message to send.
     * @param file The name of the file this message was sent from.
     * @param line The line number this message was sent from.
     */
    void debugMed(std::string msg, const char* file, const int line);

    /**
     * Logs a LOW debug message.
     *
     * @param msg The debug message to send.
     * @param file The name of the file this message was sent from.
     * @param line The line number this message was sent from.
     */
    void debugLow(std::string msg, const char* file, const int line);

    /**
     * Logs a debug message.
     *
     */
    void log(std::string strLog, VcstLogger::Level level, std::string msg,
            const char* file, const int line);

    /**
     * Returns whether or not level is enabled
     *
     */
    bool isEnabledFor(VcstLogger::Level level);

private:

    /**
     * Default Constructor
     *
     * @throws AdlCmnException if log file cannot be opened.
     */
    VcstLogger();

    /**
     * Disable Copy constructor
     *
     * @param orig object to be copied
     */
    VcstLogger(const VcstLogger &orig);

    /**
     * Disable Assignment operator
     *
     * @param rhs object whose value is to be assumed
     */
    VcstLogger &operator=(const VcstLogger &rhs);

    /**
     * Initialize the log file.
     *
     * @param path The directory the log file should be written to.
     * @param processName The name of the current process be executed.
     */
    void initLogFile(const std::string& path, const std::string& processName);

    /**
     * Generates a human-readable timestamp from the current GDT.
     *
     * @returns a string with the GDT in human-readable format.
     */
    void getCurrentHRTime(std::ostringstream& inStream);

    /**
     * Pointer to the time utility object.
     */

    VcstTime* timApi_;

    /**
     * Pointer to the output stream.  This may be cout or a log file, but
     * the writing to the stream is the same:  with the << operator.
     */
    std::ostream* outStreamPtr_;

    /**
     * Our log file.
     */
    std::ofstream logFile_;

    /**
     * Static reference to our logger instance.
     */
    static VcstLogger* ourLogger_;

    /**
     * Logging detail level.
     */
    VcstLogger::Level level_;

    /**
     * Logging output destination.
     */
    VcstLogger::LoggingDest dest_;

    /**
     * Tells whether or not DBG_HIGH is enabled.
     */
    bool dbgHighEn_;

    /**
     * Tells whether or not DBG_MED is enabled.
     */
    bool dbgMedEn_;

    /**
     * Tells whether or not DBG_LOW is enabled.
     */
    bool dbgLowEn_;

    /**
     * Strings for the different debug levels
     */
    static const std::string DEBUG_OFF_STRING;
    static const std::string DEBUG_HIGH_STRING;
    static const std::string DEBUG_MEDIUM_STRING;
    static const std::string DEBUG_LOW_STRING;

    /**
     * Strings for the different debug destinations
     */
    static const std::string DEST_STDOUT_STRING;
    static const std::string DEST_LOGFILE_STRING;
};

#define DEBUG_HIGH(logger, message) \
                logger.debugHigh(message, __FILE__, __LINE__);

#define DEBUG_MED(logger, message) \
                logger.debugMed(message, __FILE__, __LINE__);

#define DEBUG_LOW(logger, message) \
                logger.debugLow(message, __FILE__, __LINE__);

#endif
