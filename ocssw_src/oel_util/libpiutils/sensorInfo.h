/*
 * sensorInfo.h
 *
 *  Created on: Oct 29, 2013
 *      Author: dshea
 */

#ifndef SENSOR_INFO_H_
#define SENSOR_INFO_H_

#include <sensorDefs.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" { 
#endif


int32_t rdsensorinfo(int32_t sensorID, int32_t evalmask, const char *pname, void **pval);

const char* sensorId2SensorName(int sensorId);
const char* sensorId2InstrumentName(int sensorId);
const char* sensorId2PlatformName(int sensorId);
const char* sensorId2SensorDir(int sensorId);
const char* subsensorId2SubsensorDir(int subsensorId);

int sensorId2InstrumentId(int sensorId);
const char* instrumentId2InstrumentName(int instrumentId);

int sensorName2SensorId(const char* sensorName);
int instrumentPlatform2SensorId(const char* instrument, const char* platform);
int sensorId2SubsensorId(int sensorId);

const char* instrumentPlatform2SensorName(const char* instrument, const char* platform);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <string>

int sensorName2SensorId(const std::string &sensorName);
int instrumentPlatform2SensorId(const std::string &instrument, const std::string &platform);

std::string instrumentPlatform2SensorName(const std::string &instrument, const std::string &platform);

#endif

#endif /* SENSOR_INFO_H_ */
