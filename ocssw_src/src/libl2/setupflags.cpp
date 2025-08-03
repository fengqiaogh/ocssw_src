#include <iostream>
#include <setupflags.h>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
void setupflags(char *flagdef, char *flaguse, uint32_t *flagusemask, uint32_t *required, int *status) {
    int bitNum;
    char *tmpFlags;
    char *ptr, *ptr2;

    *status = 0;
    *flagusemask = 0;
    *required = 0;

    bitNum = 0;
    tmpFlags = strdup(flagdef);
    ptr = strtok(tmpFlags, ",");
    while (ptr) {
        if (ptr) {
            if ((ptr2 = strstr(flaguse, ptr))) {
                ptr2--;
                if (*ptr2 == '~')
                    *required = *required | 1<<bitNum;
                else
                    *flagusemask = *flagusemask | 1<<bitNum;
            }
        }
        ptr = strtok(NULL, ",");
        bitNum++;
        if (bitNum > 33) {
            *status = -1;
            break;
        }
    }

    free(tmpFlags);
}

int setupflags(std::string &flaguse, const std::unordered_map<std::string, int> &flag_l2_meaning_bit_dict,
                uint32_t &flagusemask, uint32_t &required) {
    std::vector<std::string> flaguse_separated;
    boost::split(flaguse_separated, flaguse, boost::is_any_of(","), boost::algorithm::token_compress_on);
    std::vector<std::string> flagusemask_list;
    std::vector<std::string> required_list;
    for (auto &word : flaguse_separated) {
        std::string flag = word;
        size_t pos = word.find("~");
        bool found = pos != std::string::npos;
        if (found) {
            flag = word.substr(pos+1);
        }
        if (flag_l2_meaning_bit_dict.count(flag) == 0) {
            std::cerr << "ERROR: Unknown flag " << flag << std::endl;
            return 1;
        }
        if (found)
            required_list.push_back(flag);
        else
            flagusemask_list.push_back(flag);
    }
    flagusemask = 0;
    required = 0;
    for (auto &flag : flagusemask_list) {
        flagusemask = flagusemask | flag_l2_meaning_bit_dict.at(flag);
    }
    for (auto &flag : required_list) {
        required = required | flag_l2_meaning_bit_dict.at(flag);
    }
    return 0;
}