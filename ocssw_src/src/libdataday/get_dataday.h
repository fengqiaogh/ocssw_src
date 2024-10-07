#include <stdint.h>

#ifndef GET_DATADAY_H
#define GET_DATADAY_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Get the plus day
 * @return int 
 */
int32_t get_plusday();

/**
 * @brief Get the datadays object
 * 
 * @param path - path to the input granule
 * @param day0 - the earliset day the granule covers
 * @param day1 - the latest day the granule covers
 * @return int 
 */
int get_datadays(const char* path, int32_t* day0, int32_t* day1);
/**
 * @brief 
 * prints error code
 * @param exitStatus 
 */
void printUsage(int32_t exitStatus);
/**
 * @brief Get the equatorial crossing time 
 * 
 * @return float 
 */
float get_equatorial_crossingTime();
/**
 * @brief Set the verbosity 
 * 
 * @param val if val > 0 verbosity is set to 1 (true)
 */
void set_verbosity(int val);
#ifdef __cplusplus
}
#endif
#endif
