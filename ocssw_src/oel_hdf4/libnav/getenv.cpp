#include <boost/algorithm/string.hpp>
#include <vector>
#include <regex>
#include <iostream>

#ifdef WINPATH
#define SEP "\\"
#define REGEX ".*%(\\w+)%.*"
#else
#define SEP "/"
#define REGEX ".*\\$(\\w+).*"
#endif

std::string filenv(const std::string &file_env_path) {
    std::vector<std::string> path_sep;
    boost::split(path_sep,file_env_path,boost::is_any_of(SEP));
    std::string out_path;
    for (const auto &el : path_sep) {
        std::string s = el;
        std::regex rgx(REGEX);
        std::smatch match;
        if (std::regex_search(s, match, rgx)) {
            std::string env_ = match[1];
            out_path += std::getenv(env_.c_str());
        } else {
            out_path += SEP + s;
        }
    }
    return out_path;
}

extern "C" void filenv_(char *inp, char *out) {
    strcpy(out, filenv(inp).c_str());
}