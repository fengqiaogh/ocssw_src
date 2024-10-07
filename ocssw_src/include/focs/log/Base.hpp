
#ifndef FOCS_LOG_BASE
#define FOCS_LOG_BASE

#include <limits>
#include <memory>
#include <string>

namespace focs {

/*! \enum focs::LogSeverity

    \brief Log severities, for use with focs::Log

    https://en.wikipedia.org/wiki/Syslog#Severity_level

    Value   Severity        Keyword   Description/Examples
    0       Emergency       emerg     System is unusable
                                      This level should not be used by applications.
    1       Alert           alert     Should be corrected immediately
                                      Loss of the primary ISP connection.
                                      E.g., Ski Haus Delta has not reported status within status_timeout (120)
    2       Critical        crit      Critical conditions
                                      A failure in the system's primary application.
                                      E.g., Ski Haus Delta reports temperature < low_critical (30)
    3       Error           err       Error conditions
                                      E.g., An application has exceeded its file storage limit and attempts to write are failing.
                                      E.g., Ski Haus Delta reports temperature < low_error (32)
    4       Warning         warn      May indicate that an error will occur if action is not taken.
                                      E.g., A non-root file system has only 2GB remaining.
                                      E.g., Ski Haus Delta reports temperature < low_warning (36)
    5       Notice          notice    Events that are unusual, but not error conditions.
                                      E.g., Ski Haus Delta reports temperature < low_notice(50)
    6       Informational   info      Normal operational messages that require no action.
                                      An application has started, paused or ended successfully.
                                      E.g., Ski Haus Delta reports temperature 60
    7       Debug           debug     Information useful to developers for debugging the application.

    The meaning of severity levels other than Emergency and Debug are relative to the application.
    For example, if the purpose of the system is to process transactions to update customer account
    balance information, an error in the final step should be assigned Alert level. However, an error
    occurring in an attempt to display the ZIP code of the customer may be assigned Error or even
    Warning level.
*/
enum class LogSeverity : int {
    debug = 7,      //!< debug severity
    info = 6,       //!< info severity
    notice = 5,     //!< notice severity
    warning = 4,    //!< warning severity
    error = 3,      //!< error severity
    critical = 2,   //!< critical severity
    alert = 1,      //!< alert severity
    emergency = 0,  //!< emergency severity

    warn = warning,    //!< warning severity (alias for warning)
    err = error,       //!< error severity (alias for error)
    crit = critical,   //!< critical severity (alias for critical)
    emerg = emergency, //!< emergency severity (alias for emergency)

    min = std::numeric_limits<int>::max(), //!< minimum severity possible
    max = std::numeric_limits<int>::min()  //!< maximum severity possible
};
/*!
    \class focs::LogFacility

    \brief Base class to create new loggers
*/
class LogFacility {
    public:
        /*!  \brief Empty, virtual destructor */
        virtual ~LogFacility(){}

        /*!
            \brief Print a message to this log facility, with a given severity

            \param severity numeric severity level
            \param s message to log
        */
        virtual void log(int severity, const std::string& s) = 0;
        /*!
            \brief Wrapper for the above function, converting the LogSeverity enum
                to its underlying numeric severity.  This shouldn't need to be
                overridden manually.
        
            \param severity numeric severity level
            \param s message to log
        */
        virtual void log(LogSeverity severity, const std::string& s){ this->log(static_cast<int>(severity), s); }
};
/*!
    \class focs::LogTarget

    \brief Simple class containing a log facility and the severity ranges it should receive

    This is used internally for focs::Log and isn't meant for end users.
*/
//! \cond
class LogTarget {
    public:
        LogTarget(int min_severity, int max_severity, std::unique_ptr<LogFacility> f) : min_severity_{min_severity}, max_severity_{max_severity}, log_facility_{std::move(f)} {}
        LogTarget(int min_severity, std::unique_ptr<LogFacility> f) : min_severity_{min_severity}, log_facility_{std::move(f)}{}
        LogTarget(std::unique_ptr<LogFacility> f) : log_facility_{std::move(f)}{}

        // move
        LogTarget(LogTarget&& other){
            swap(std::move(other));
        }
        LogTarget& operator=(LogTarget&& other){
            swap(std::move(other));
            return *this;
        }

        // no copying
        LogTarget(const LogTarget&) = delete;
        LogTarget& operator=(const LogTarget&) = delete;

        int min_severity(){ return min_severity_; }
        int max_severity(){ return max_severity_; }
        LogFacility& log_facility(){ return *log_facility_; }

    private:
        int min_severity_{static_cast<int>(LogSeverity::info)};
        int max_severity_{static_cast<int>(LogSeverity::max)};
        std::unique_ptr<LogFacility> log_facility_;

        void swap(LogTarget&& other){
            auto min = min_severity_;
            auto max = max_severity_;
            auto log = std::move(log_facility_);
            min_severity_ = other.min_severity_;
            max_severity_ = other.max_severity_;
            log_facility_ = std::move(other.log_facility_);
            other.min_severity_ = min;
            other.max_severity_ = max;
            other.log_facility_ = std::move(log);
        }
};
//! \endcond
} // namespace focs

#endif // FOCS_LOG_BASE

