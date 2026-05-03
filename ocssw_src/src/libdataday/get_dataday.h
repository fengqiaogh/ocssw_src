#include <stdint.h>
#include <time.h>
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
 * @param sensorID - the sensor ID of the input granule
 * @param dayNight - whether the granule is a day or night granule
 * @param starttime - the start time of the granule in seconds since 1 Jan. 1970 00:00:00 GMT
 * @param equatorialCrossingTime - the equatorial crossing time in hours (output)
 * @param plusDay - whether to add a day to the dataday output of get_datadays (output)
 * @return float 
 */
void getEquatorCrossingTime(int32_t sensorID, bool dayNight, time_t starttime, float* equatorialCrossingTime,
                            int32_t* plusDay);
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
