#ifndef CALPROTOI_H /* avoid re-inclusion */
#define CALPROTOI_H

int32 get_cal(char *cal_path, int16 syear, int16 sday, int16 eday,
        int32 msec, char *dtype, int16 *tdi, int16 *cal_year,
        int16 *cal_day, int16 *ref_year, int16 *ref_day,
        int16 *ref_min, float32 temps[256][BANDS],
        float32 scan_mod[2][1285], float32 mirror[2][BANDS],
        float64 *tfactor_const, float64 *tfactor_linear,
        float64 *tfactor_quadratic, float32 *cal_offset,
        float32 counts[BANDS][4][5], float32 rads[BANDS][4][5]);

int32 read_time_data(int32 fid, int16 *year, int16 *day, int32 *msec,
        int32 *elts);

int32 get_index(int32 fid, int16 syear, int16 sday, int16 eday, int32 msec,
        int16 *cal_year, int16 *cal_day);

int32 read_parm_data(int32 fid, int32 sdfid, int32 index, int32 idoffs[8][16],
        float32 gains[8][16], float32 temps[256][8],
        float32 scan_mod[2][1285], float64 *tfactor_const,
        float64 *tfactor_linear, float64 *tfactor_quadratic,
        float32 *cal_offset, float32 mirror[2][8],
        int16 tdi_list[8][4]);

void calc_knees(int16 *tdi, int16 tdi_list[8][4], int32 idoffs[8][16],
        float32 gains[8][16], float32 counts[8][4][5],
        float32 rads[8][4][5]);

void setup_scanmod(char *dtype, float32 scan_mod[2][1285]);

void sort_srads(float32 *srads, int32 *oindex);

int32 get_ref_time(int32 sdfid, int16 *ref_year, int16 *ref_day, int16 *ref_min);

#endif  /* CALPROTOI_H */
