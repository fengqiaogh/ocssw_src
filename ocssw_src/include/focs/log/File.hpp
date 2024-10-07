
#ifndef FOCS_LOG_FILE
#define FOCS_LOG_FILE

#include "focs/Log.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace focs {
/*!
    \class focs::FileLogger

    \brief Log facility for printing to a file
*/
class FileLogger : public LogFacility {
    public:
        /*!
            \brief Open a file at the given path, optionally with a specific file mode
        
            \param filename Path to log file
            \param mode Mode with which to open the file
        */
        FileLogger(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out) : path_{filename}, handle_{filename, mode} {}
        /*!  \brief Clean-up destructor to flush and close file handle */
        ~FileLogger() override;

        void log(int severity, const std::string& s) override;
    private:
        std::string path_;
        std::ofstream handle_;
};
}

#endif // FOCS_LOG_FILE
