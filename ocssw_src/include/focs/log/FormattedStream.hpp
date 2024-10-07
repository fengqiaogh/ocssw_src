
#ifndef FOCS_LOG_FORMATTEDSTREAM
#define FOCS_LOG_FORMATTEDSTREAM

#include "focs/Log.hpp"

#include <functional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace focs {
/*!
    \class focs::FormatCommand

    \brief Class for storing logging functions (internal, not for end users)
*/
//! \cond
class FormatCommand {
    public:
        virtual ~FormatCommand() {}
        virtual void log(std::ostream& stream, const int severity, const std::string& s) const = 0;
        virtual bool is_user_string() { return false; }
    private:
};
//! \endcond

/*!
    \class focs::FormattedStreamLogger

    \brief Log facility allowing text to be generated around messages

    \section formattedstreamlogger-format Format String

    \code{md}
    %s - message from user

    %% - literal %

    %{function_name}f - function from map_of_functions
    %[function_index]f - function from vector_of_functions

    %l - severity lowercase
    %u - severity uppercase
    %m - severity mixed case
    %n - severity number

    %t - time
    %{time_string} - time
    %{local|gmt|utc,time_string} - time
    \endcode

    If the user message (\%s) is omitted, the entire format string is treated
    as a prefix and the user message is automatically added to the end.
*/
class FormattedStreamLogger : public LogFacility {
    using named_function = std::function<std::string(const std::string& name, std::ostream& stream, const int severity, const std::string& str)>;
    using map_of_functions = std::unordered_map<std::string, named_function>;

    using indexed_function = std::function<std::string(const size_t index, std::ostream& stream, const int severity, const std::string& str)>;
    using vector_of_functions = std::vector<indexed_function>;

    public:
        /*!
            \brief Create a new log facility with the given format, with optional function containers and a specific output stream
        
            \param format string with which to format messages
            \param named_functions map of functions referred to by the format string
            \param numbered_functions vector of functions referred to by the format string
            \param stream stream to which to print
        */
        explicit FormattedStreamLogger(const std::string& format, map_of_functions named_functions={}, vector_of_functions numbered_functions={}, std::ostream& stream=std::cout);
        /*!
            \brief Create a new log facility without external function calls

            \param format string with which to format messages
            \param stream stream to which to print
        */
        FormattedStreamLogger(const std::string& format, std::ostream& stream) : FormattedStreamLogger(format, {}, {}, stream) {}

        ~FormattedStreamLogger() override;

        void log(int severity, const std::string& s) override;
    private:
        // const std::string& format_;
        std::vector<std::unique_ptr<FormatCommand>> commands_{};
        std::vector<std::unique_ptr<FormatCommand>>::const_iterator after_user_string_{commands_.cend()};
        std::ostream& stream_{std::cout};

        map_of_functions named_functions_{};
        vector_of_functions numbered_functions_{};

        void parse_format_string(const std::string& format);

        bool user_string_is_last_{true};
        bool last_format_char_is_newline_{false};
        bool last_character_was_newline_{true};
};
}

#endif //  FOCS_LOG_FORMATTEDSTREAM

