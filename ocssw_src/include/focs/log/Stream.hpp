

#ifndef FOCS_LOG_STREAM
#define FOCS_LOG_STREAM

#include "focs/Log.hpp"

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

namespace focs {
    /*!
        \class focs::FileStreamLogger
    
        \brief Log facility for printing to a std::ostream
    */
class FileStreamLogger : public LogFacility {
    public:
        /*!  \brief Default constructor, print to std::cout */
        FileStreamLogger() = default;
        /*!
            \brief Print to a given std::ostream
            \param stream Stream to which to print
        */
        explicit FileStreamLogger(std::ostream& stream) : stream_{stream} {}
        ~FileStreamLogger() override;

        void log(int severity, const std::string& s) override;
    private:
        std::ostream& stream_{std::cout};
};
}

#endif // FOCS_LOG_STREAM

