#include "l2_utils.hpp"
static char arr[256];
/**
 * @brief Convert a std::string to a char*
 *  strdup is not used because it uses malloc and we are using a static array.
  * @param data 
 * @return char* 
 */
char* convert_string(const std::string& data) {
    strncpy(arr, data.c_str(), sizeof(arr) - 1);
    return arr;
}

