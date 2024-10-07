#ifndef CALPROTOI_H /* avoid re-inclusion */
#define CALPROTOI_H

#include <stdint.h>

int32_t get_cal_osmi(char *cal_path, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day, int32_t *cal_msec,
        float *eoffset, float *egain, float *temp,
        float *mirror, float *t_const,
        float *t_linear, float *t_quadratic,
        float *slope, float *dc,
        float *sm);


int32_t get_index_osmi(int32_t fid, int16_t syear, int16_t sday, int16_t eday, int32_t msec,
        int16_t *cal_year, int16_t *cal_day, int32_t *cal_msec);

int32_t read_parm_data_osmi(int32_t fid, int32_t sdfid, int32_t index, float *eoffset,
        float *egain, float *temp, float *mirror,
        float *t_const, float *t_linear, float *t_quadratic,
        float *slope, float *dc, float *sm);

#endif  /* CALPROTOI_H */
