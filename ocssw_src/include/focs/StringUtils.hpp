#ifndef FOCS_STRINGUTILS
#define FOCS_STRINGUTILS


#include <functional>

#include <string>
#include <vector>

namespace focs {

/*!
    \class focs::StringUtils

    \brief Static functions for parsing/processing/modifying strings

    \section string_utils-ex1 String Replacement
    \snippet tests/string_utils.cpp focs::StringUtils example1

    \section string_utils-ex2 String-to-Vector Examples
    \snippet tests/string_utils.cpp focs::StringUtils example2

    \section string_utils-ex3 Quote/Bracket Stripping
    \snippet tests/string_utils.cpp focs::StringUtils example3
*/
class StringUtils {
    public:
        /*!
            \brief Replace $OCSSWROOT, $OCDATAROOT, etc, in string (returns a copy)

            \param str String containing environment variables to replace

            \return A new string with the known environment variables replaced
        */
        static std::string replace_oc_roots(const std::string& str);
        /*!
            \brief Replace $OCSSWROOT, $OCDATAROOT, etc, in string, (in-place)

            \param str String containing environment variables to replace

            \return Reference to the input string with environment variables replaced
        */
        static std::string& replace_oc_roots(std::string& str);
        /*!
            \brief Replace paths in the input with $OCSSWROOT, $OCDATAROOT, etc, (returns a copy)

            \param str String containing paths to replace

            \return A new string with the known environment variables inserted in place of absolute paths
        */
        static std::string insert_oc_roots(const std::string& str);
        /*!
            \brief Replace paths in the input with $OCSSWROOT, $OCDATAROOT, etc, (in-place)

            \param str String containing paths to replace

            \return reference to the input string with paths replaced
        */
        static std::string& insert_oc_roots(std::string& str);
        /*!
            \brief Replace all occurrences of a string (returns a copy)

            \param haystack String with text to replace
            \param needle What to find within haystack
            \param replacement What to insert into haystack

            \return Copy of the input string with substrings replaced
        */
        static std::string replace_all(const std::string& haystack, const char *needle, const char *replacement);
        /*!
            \brief Replace all occurrences of a string (in-place)

            \param haystack String with text to replace
            \param needle What to find within haystack
            \param replacement What to insert into haystack

            \return Reference to input string
        */
        static std::string& replace_all(std::string& haystack, const char *needle, const char *replacement);

        /*!
            \brief Remove start/end characters from ends of a string (returns a copy)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify
            \param first Expected first character of haystack
            \param last Expected last character of haystack

            \return Copy of the input string with enclosure, if found, removed
        */
        static std::string strip_enclosure(const std::string& haystack, const char first, const char last);
        /*!
            \brief Remove start/end characters from ends of a string (in-place)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify
            \param first Expected first character of haystack
            \param last Expected last character of haystack

            \return reference to the input string with enclosure, if found, removed
        */
        static std::string& strip_enclosure(std::string& haystack, const char first, const char last);

        /*!
            \brief Remove square brackets from ends of input string (returns a copy)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify

            \return reference to the input string with square brackets, if found, removed
        */
        static std::string strip_brackets(const std::string& haystack);
        /*!
            \brief Remove square brackets from ends of input string (in-place)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify

            \return reference to the input string with square brackets, if found, removed
        */
        static std::string& strip_brackets(std::string& haystack);
        /*!
            \brief Remove single- and double-quotes from ends of input string (returns a copy)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify

            \return copy of the input string with quotes brackets, if found, removed
        */
        static std::string strip_quotes(const std::string& haystack);
        /*!
            \brief Remove single- and double-quotes from ends of input string (in-place)

            Only removes enclosure if first and last characters match as expected.

            \param haystack String to modify

            \return reference to the input string with quotes, if found, removed
        */
        static std::string& strip_quotes(std::string& haystack);

        /*!
            \brief Convert delimited string to vector

            This is declared to provide template specialization.  The provided
            specializations are below.  All default specializations do not
            check for validity, e.g., entire string matched, string wasn't
            empty, value wasn't narrowed.

            + std::string (no conversion)
            + float (std::stof)
            + double (std::stod)
            + int (std::stoi)
            + long (std::stol)
            + long double (std::stold)
            + long long (std::stoll)
            + unsigned (std::stoul)
            + unsigned long (std::stoul)
            + unsigned long long (std::stoull)

            \tparam T type of elements in output vector
            \param str Delimited string
            \param delims The substring separating elements with the string (defaults to comma)
            \param merge_multiple_delims Whether multiple, consecutive
                delimiters (i.e., an empty value) should be considered one
                delimiter.
            \param default_value What to insert if merging multiple delimiters.

            \return Newly allocated vector of the template type
        */
        template<typename T>
        static std::vector<T> stov(std::string& str, const std::string& delims=",", bool merge_multiple_delims=false, T default_value=T{});

        /*!
            \brief Convert delimited string to vector

            This should be implemented for custom types not already specialized.

            \tparam T type of elements in output vector
            \param str Delimited string
            \param delims The substring separating elements with the string (defaults to comma)
            \param merge_multiple_delims Whether multiple, consecutive
                delimiters (i.e., an empty value) should be considered one
                delimiter.
            \param default_value What to insert if merging multiple delimiters.
            \param parser The function to use to convert the elements' strings to the desired type

            \return Newly allocated vector of the template type
        */
        template<typename T>
        static std::vector<T> stov(std::string& str, const std::string& delims, bool merge_multiple_delims, T default_value, std::function<T(const std::string&)> parser);

        // TODO: provide things like this, to make it easier to customize? (All this does is reorder the params and provide some defaults.
        // template<typename T>
        // static std::vector<T> stov(std::string& str, const std::string& delims, std::function<T(const std::string&)> parser, bool merge_multiple_delims=false, T default_value=T{});
    private:

};


} // namespace focs

#endif // FOCS_STRINGUTILS
