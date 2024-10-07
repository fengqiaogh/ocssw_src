/* 
 * File:   GetStationInfo.h
 * Author: dshea
 *
 * Created on October 13, 2011, 4:01 PM
 */

#ifndef GETSTATIONINFO_H
#define GETSTATIONINFO_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char code[5];
    char *data_center;
    char *station_name;
    float station_latitude;
    float station_longitude;
} StationInfo;

int GetStationInfo(char *stationInfoFile, StationInfo *stationInfo);


#ifdef __cplusplus
}
#endif

#endif /* GETSTATIONINFO_H */

