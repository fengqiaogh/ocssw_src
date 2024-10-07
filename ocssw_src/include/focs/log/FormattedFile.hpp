
#ifndef FOCS_LOG_FORMATTEDFILE
#define FOCS_LOG_FORMATTEDFILE

#include "focs/Log.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace focs {
/*!
    \class focs::FormattedFileLogger

    \brief Log facility for printing to a file, allowing text to be generated around messages

    For the detailed description of the format parameter, see the documentation
    for focs::FormattedStreamLogger.
*/
class FormattedFileLogger : public FormattedStreamLogger {
    using named_function = std::function<std::string(const std::string& name, std::ostream& stream, const int severity, const std::string& str)>;
    using map_of_functions = std::unordered_map<std::string, named_function>;

    using indexed_function = std::function<std::string(const size_t index, std::ostream& stream, const int severity, const std::string& str)>;
    using vector_of_functions = std::vector<indexed_function>;

    public:
        /*!
            \brief Create a new log facility without external function calls

            \param filename string with which to format messages
            \param format string with which to format messages
            \param named_functions map of functions referred to by the format string
            \param numbered_functions vector of functions referred to by the format string
        */
        FormattedFileLogger(const std::string& filename, const std::string& format, map_of_functions named_functions={}, vector_of_functions numbered_functions={}) : FormattedStreamLogger(format, named_functions, numbered_functions, handle_), path_{filename}, handle_{filename, std::ios_base::out | std::ios_base::app}  {}

        /*!
            \brief Create a new log facility without external function calls

            \param filename string with which to format messages
            \param mode file mode to use when opening the file
            \param format string with which to format messages
            \param named_functions map of functions referred to by the format string
            \param numbered_functions vector of functions referred to by the format string
        */
        FormattedFileLogger(const std::string& filename, std::ios_base::openmode mode, const std::string& format, map_of_functions named_functions={}, vector_of_functions numbered_functions={}) : FormattedStreamLogger(format, named_functions, numbered_functions, handle_), path_{filename}, handle_{filename, mode}  {}
//std::ios_base::out | std::ios_base::app)

        ~FormattedFileLogger() override {
            handle_.close();
        }
    private:
        std::string path_;
        std::ofstream handle_;
};
}

#endif // FOCS_LOG_FORMATTEDFILE
