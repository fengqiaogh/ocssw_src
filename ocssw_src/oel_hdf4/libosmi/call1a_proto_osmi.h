#ifndef CALL1APROTO_OSMI_H
#define CALL1APROTO_OSMI_H

#include <stdint.h>

int32_t calibrate_l1a_osmi(char *cal_path,
        int16_t syear,
        int16_t sday,
        int32_t smsec,
        int16_t eday,
        int32_t msec,
        char *dtype, /* "","SOL","TDI","IGC","LUN" */
        int32_t st_samp, /* 1 */
        int32_t nsamp, /* 1044 */
        int32_t fpixel, /* 0-95 */
        int16_t gain[4], /* quadrant gain settings from modified ISD */
        int16_t offset, /* quadrant offset settings from mod ISD */
        int16_t scan_temp, /* sensor temps from mod ISD */
        int16_t *l1a_data,
        float *l1b_data,
        cal_mod_struc *cal_mod);

int jul2jday(int year, int yday);

int cal2jday(int year, int month, int mday);

#endif /* CALL1APROTO_H */
