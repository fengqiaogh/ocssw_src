#ifndef CALL1APROTO_H
#define CALL1APROTO_H

#include "l1a.h"
#include <stdint.h>

int32_t calibrate_l1a(char *cal_path, int16_t syear, int16_t sday, int32_t smsec,
        int16_t eday, int32_t msec, char *dtype, int32_t nsta, int32_t ninc,
        int32_t nsamp, float *dark_mean, int16_t *gain, int16_t *tdi,
        int16_t *scan_temp, float *inst_temp, int16_t side,
        int16_t *l1a_data, float *l1b_data, cal_mod_struc *cal_mod);

int jul2jday(int year, int yday);

int cal2jday(int year, int month, int mday);

#endif /* CALL1APROTO_H */
