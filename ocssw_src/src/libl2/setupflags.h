/*
 * File:   setupflags.h
 * Author: dshea
 *
 * Created on January 31, 2012, 8:47 AM
 */

#ifndef SETUPFLAGS_H
#define SETUPFLAGS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/**
 * @brief Sets up flag masks based on flag definitions and usage
 * @param flagdef Comma-separated string of flag definitions
 * @param flaguse String specifying which flags to use and how
 * @param flagusemask Output bitmask for flags to use
 * @param required Output bitmask for required flags
 * @param status Output status code (0 for success, -1 for error)
 * @param BITS Array of bit values for each flag position
 * 
 * Parses flag definitions and usage to create bitmasks indicating which flags
 * should be used and which are required. A '~' prefix in flaguse indicates
 * a required flag.
 */
void setupflags(char *flagsdef, char *flaguse, uint32_t *flagusemask, uint32_t *required, int *status);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
#include <unordered_map>
#include <string>
/**
 * @brief Sets up flag masks based on flag usage and meaning dictionary
 * @param flaguse String specifying which flags to use and how
 * @param flag_l2_meaning_bit_dict Dictionary mapping flag names to bit positions
 * @param flagusemask Output bitmask for flags to use
 * @param required Output bitmask for required flags
 * @return 0 on success, 1 if unknown flag encountered
 * 
 * Modern C++ version that uses a dictionary to look up flag bit positions.
 * Parses comma-separated flag usage string where '~' prefix indicates
 * a required flag. Validates all flags exist in dictionary.
 */
int setupflags(std::string &flaguse, const std::unordered_map<std::string, int> & flag_l2_meaning_bit_dict,
                uint32_t &flagusemask, uint32_t &required);
#endif
#endif /* SETUPFLAGS_H */
