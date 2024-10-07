#ifndef CALIBRATE_VIIRS_H
#define CALIBRATE_VIIRS_H

#ifdef __cplusplus
extern "C" {
#endif

int load_fcal_lut(char* calfile, int64_t UTC58usec, double ****ftable);

#ifdef __cplusplus
}
#endif
#endif