/* 
 * File:   sst_dsdi.h
 * Author: swbaile1
 *
 * Created on March 29, 2019, 1:09 PM
 */

#ifndef SST_DSDI_H
#define SST_DSDI_H

#ifdef __cplusplus
extern "C" {
#endif
    float dust_correction(float dustExtinction, float csenz, float BT39, float BT85, float BT11, 
            float BT12, int32_t sensorID);


#ifdef __cplusplus
}
#endif

#endif /* SST_DSDI_H */

