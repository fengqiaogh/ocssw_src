#ifndef CALL1APROTO_H
#define CALL1APROTO_H

int32 calibrate_l1a(char *cal_path, int16 syear, int16 sday, int32 smsec,
        int16 eday, int32 msec, char *dtype, int32 st_samp,
        int32 nsamp, int16 *dark_rest, float32 *dark_mean, int16 *gain,
        int16 *tdi, int16 *scan_temp, int16 side, int16 *l1a_data,
        float32 *l1b_data, cal_mod_struc *cal_mod);

int jul2jday(int year, int yday);

int cal2jday(int year, int month, int mday);

#endif /* CALL1APROTO_H */
