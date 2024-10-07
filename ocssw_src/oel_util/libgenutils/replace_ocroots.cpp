#include "genutils.h"

#include <stdlib.h>
#include <string.h>
#include <string>

/**
 * modify the input string with environment variable names substituted for the 
 * values for:  OCVARROOT, OCDATAROOT and OCSSWROOT 
 * @param str string to do substitutions on
 */
void replaceOCroots(std::string& str) {
    constexpr size_t numSubstitutes = 3;
    char *envSubstitutes[numSubstitutes][2] = {
        {"$OCVARROOT", getenv("OCVARROOT")},
        {"$OCDATAROOT", getenv("OCDATAROOT")},
        {"$OCSSWROOT", getenv("OCSSWROOT")}
    };
    
    for(size_t sub = 0; sub < numSubstitutes; sub++) {
        if(envSubstitutes[sub][1] == NULL)
            continue;
        size_t pos;
        while((pos = str.find(envSubstitutes[sub][1])) != std::string::npos) {
            str.replace(pos, strlen(envSubstitutes[sub][1]), envSubstitutes[sub][0]);
        }
    }
}


/**
 * allocate a new string with environment variable names substituted for the 
 * values for:  OCVARROOT, OCDATAROOT and OCSSWROOT 
 * @param inStr string to do substitutions on
 * @return sanitized string that must be freed by the caller.
 */
extern "C" char* replace_ocroots(const char* inStr) {
    std::string str = inStr;
    replaceOCroots(str);
    return strdup(str.c_str());
}

